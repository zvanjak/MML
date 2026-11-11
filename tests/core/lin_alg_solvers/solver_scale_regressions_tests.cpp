#include <catch2/catch_all.hpp>

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#elif defined(MML_USE_SINGLE_HEADER)
#include <MML.h>
#else
#include <mml/core/LinAlgEqSolvers/LinAlgDirect.h>
#include <mml/core/LinAlgEqSolvers/LinAlgQR.h>
#include <mml/core/LinAlgEqSolvers/LinAlgSVD.h>
#include <mml/core/LinAlgEqSolvers/LinAlgComplexSVD.h>
#include <mml/core/LinAlgEqSolvers/LinAlgEqSolvers_iterative.h>
#include <mml/core/SparseSolvers/GMRES.h>
#include <mml/core/SparseSolvers/ConjugateGradient.h>
#include <mml/core/SparseSolvers/Preconditioners.h>
#include <mml/algorithms/MatrixAlg/Decompositions.h>
#endif

#include <cmath>
#include <limits>
#include <memory>

namespace MML::Tests::Core::SolverScaleRegressions
{
	namespace
	{
		Real RelativeResidual(const Matrix<Real>& a, const Vector<Real>& x, const Vector<Real>& b)
		{
			Real scale = 0;
			for (int row = 0; row < a.rows(); ++row) {
				scale = std::max(scale, std::abs(b[row]));
				for (int col = 0; col < a.cols(); ++col)
					scale = std::max(scale, std::abs(a(row, col)));
			}
			Real residual = 0, rhs = 0;
			for (int row = 0; row < a.rows(); ++row) {
				Real value = -b[row] / scale;
				for (int col = 0; col < a.cols(); ++col)
					value += (a(row, col) / scale) * x[col];
				residual = std::hypot(residual, value);
				rhs = std::hypot(rhs, b[row] / scale);
			}
			return rhs == 0 ? residual : residual / rhs;
		}

		SparseMatrix::SparseMatrixCSR<Real> Sparse(const Matrix<Real>& a)
		{
			SparseMatrix::SparseMatrixCOO<Real> coo(a.rows(), a.cols());
			for (int row = 0; row < a.rows(); ++row)
				for (int col = 0; col < a.cols(); ++col)
					if (a(row, col) != 0)
						coo.addEntry(row, col, a(row, col));
			return SparseMatrix::SparseMatrixCSR<Real>(coo);
		}

		class ScaledIdentity : public SparseSolvers::Preconditioner<Real>
		{
		public:
			void apply(const Vector<Real>& r, Vector<Real>& z) const override
			{
				z = r;
				for (auto& value : z)
					value *= Real(1e-20);
			}
		};
	}

	TEST_CASE("L01 GMRES breakdown is independent of RHS and operator scale", "[ke4v][GMRES][scale]")
	{
		for (const Real operatorScale : {Real(1), Real(1e-20)}) {
			for (const Real rhsScale : {Real(1), Real(1e16)}) {
				CAPTURE(operatorScale, rhsScale);
				const Matrix<Real> a(2, 2, {operatorScale, 0, 0, 2 * operatorScale});
				const Vector<Real> b{rhsScale, rhsScale};
				Vector<Real> x(2, Real(0));
				SparseSolvers::SolverConfig<Real> config;
				config.tolerance = Precision::EqualityTolerance;
				config.absTolerance = 0;
				config.maxIterations = 2;
				const auto result = SparseSolvers::solveGMRES(Sparse(a), b, x, 2, config);
				CHECK(result.converged);
				CHECK(result.iterations_used == 2);
				CHECK(RelativeResidual(a, x, b) < Precision::EqualityTolerance);
				CHECK(std::abs(result.relativeResidual - RelativeResidual(a, x, b)) < Precision::EqualityTolerance);
			}
		}
	}

	TEST_CASE("L01 GMRES spends actual restart work and preserves zero-operator iterate", "[ke4v][GMRES]")
	{
		SECTION("early inner exits do not consume fictitious restart work") {
			const Matrix<Real> a(3, 3, {1, 0, 0, 0, 2, 0, 0, 0, 3});
			const Vector<Real> b{1, 1, 1};
			Vector<Real> x(3, Real(0));
			SparseSolvers::SolverConfig<Real> config;
			config.tolerance = Precision::Eps;
			config.absTolerance = Precision::RoundoffZeroThreshold;
			config.maxIterations = 4;
			SparseSolvers::GMRES<Real> solver(2);
			solver.setConfig(config);
			solver.setPreconditioner(std::make_shared<ScaledIdentity>());
			const auto result = solver.solve(Sparse(a), b, x);
			CHECK(result.status == AlgorithmStatus::MaxIterationsExceeded);
			CHECK(result.iterations_used == config.maxIterations);
			CHECK(std::abs(result.relativeResidual - RelativeResidual(a, x, b)) < Precision::EqualityTolerance);
		}
		SECTION("a singular projected solve cannot invent a correction") {
			const Matrix<Real> a(2, 2);
			const Vector<Real> b{1, 2};
			Vector<Real> x(2, Real(0));
			const auto result = SparseSolvers::solveGMRES(Sparse(a), b, x, 2);
			CHECK_FALSE(result.converged);
			CHECK(result.status == AlgorithmStatus::NumericalInstability);
			CHECK(x[0] == 0);
			CHECK(x[1] == 0);
			CHECK(result.relativeResidual == 1);
		}
		SECTION("a zero work budget still checks the initial guess") {
			const Matrix<Real> a = Matrix<Real>::Identity(2);
			const Vector<Real> b{1, 2};
			Vector<Real> x = b;
			SparseSolvers::SolverConfig<Real> config;
			config.maxIterations = 0;
			const auto result = SparseSolvers::solveGMRES(Sparse(a), b, x, 2, config);
			CHECK(result.converged);
			CHECK(result.iterations_used == 0);
			CHECK(RelativeResidual(a, x, b) == 0);
		}
	}

	TEST_CASE("L02 finite condition-one extreme LU has finite accurate solution", "[ke4v][LUSolver][scale]")
	{
		const Real scale = std::numeric_limits<Real>::max() * Real(0.75);
		const Matrix<Real> a(2, 2, {scale, scale, -scale, scale});
		const Vector<Real> b{scale, 0};
		SECTION("owning LU") {
			LUSolver<Real> solver(a);
			const auto x = solver.Solve(b);
			CHECK(std::isfinite(x[0]));
			CHECK(std::isfinite(x[1]));
			CHECK(std::abs(x[0] - Real(0.5)) < Precision::EqualityTolerance);
			CHECK(std::abs(x[1] - Real(0.5)) < Precision::EqualityTolerance);
			CHECK(RelativeResidual(a, x, b) < Precision::EqualityTolerance);
			REQUIRE_THROWS_AS(solver.GetDecomposition(), MatrixNumericalError);
		}
		SECTION("in-place LU") {
			Matrix<Real> factors(a);
			LUSolverInPlace<Real> solver(factors);
			const auto x = solver.Solve(b);
			CHECK(RelativeResidual(a, x, b) < Precision::EqualityTolerance);
		}
	}

	TEST_CASE("L02 LU exports original-scale factors when they are representable", "[ke4v][LUSolver][scale]")
	{
		if (std::numeric_limits<Real>::max_exponent10 < 200)
			SKIP("Double-range witness");
		for (const Real scale : {Real(1), Real(1e200), Real(1e-200)}) {
			const Matrix<Real> a(2, 2, {scale, 2 * scale, 3 * scale, 4 * scale});
			LUSolver<Real> solver(a);
			const auto& factors = solver.GetDecomposition();
			auto permuted = a;
			const auto& pivots = solver.GetPivotIndices();
			for (int row = 0; row < 2; ++row)
				for (int col = 0; col < 2; ++col)
					if (row != pivots[row])
						std::swap(permuted(row, col), permuted(pivots[row], col));
			for (int row = 0; row < 2; ++row)
				for (int col = 0; col < 2; ++col) {
					Real product = 0;
					for (int k = 0; k < 2; ++k) {
						const Real l = row == k ? 1 : row > k ? factors(row, k) : 0;
						const Real u = k <= col ? factors(k, col) / scale : 0;
						product += l * u;
					}
					CHECK(std::abs(product - permuted(row, col) / scale) < Precision::EqualityTolerance);
				}
		}
	}

	TEST_CASE("L02 LU equilibration does not erase representable small rows", "[ke4v][LUSolver][scale]")
	{
		if (std::numeric_limits<Real>::max_exponent10 < 200)
			SKIP("Double-range witness");
		const Matrix<Real> a(2, 2, {Real(1e200), 0, 0, Real(1e-200)});
		LUSolver<Real> solver(a);
		const auto x = solver.Solve(Vector<Real>{Real(1e200), Real(1e-200)});
		CHECK(std::abs(x[0] - 1) < Precision::EqualityTolerance);
		CHECK(std::abs(x[1] - 1) < Precision::EqualityTolerance);
		CHECK(std::abs(solver.det() - 1) < Precision::EqualityTolerance);
	}

	TEST_CASE("L03 GaussJordan rejects all-NaN input before selecting a pivot", "[ke4v][GaussJordanSolver][unsafe-red]")
	{
		const Real nan = std::numeric_limits<Real>::quiet_NaN();
		Matrix<Real> a(2, 2, {nan, nan, nan, nan});
		Vector<Real> b{1, 2};
		REQUIRE_THROWS_AS(GaussJordanSolver<Real>::SolveInPlace(a, b), MatrixNumericalError);
	}

	TEST_CASE("L04 QR rejects empty input and recognizes a zero terminal pivot", "[ke4v][QRSolver][unsafe-red]")
	{
		SECTION("empty") {
			const Matrix<Real> empty;
			REQUIRE_THROWS_AS(QRSolver<Real>(empty), MatrixDimensionError);
		}
		SECTION("zero 1x1") {
			QRSolver<Real> solver(Matrix<Real>(1, 1, {0}));
			CHECK(solver.isSingular());
			REQUIRE_THROWS_AS(solver.Solve(Vector<Real>{1}), SingularMatrixError);
		}
	}

	TEST_CASE("L06 Detailed solvers reject nonfinite inputs and never invent residuals", "[ke4v][Detailed]")
	{
		const Real nan = std::numeric_limits<Real>::quiet_NaN();
		LinearSolverConfig config;
		config.check_finite = true;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		const Matrix<Real> identity = Matrix<Real>::Identity(2);
		const Vector<Real> badRhs{nan, 1};
		CHECK(LUSolveDetailed(identity, badRhs, config).status == AlgorithmStatus::NumericalInstability);
		CHECK(GaussJordanSolveDetailed(identity, badRhs, config).status == AlgorithmStatus::NumericalInstability);
		CHECK(CholeskySolveDetailed(identity, badRhs, config).status == AlgorithmStatus::NumericalInstability);
		const BandDiagonalMatrix bandIdentity(2, 0, 0, Matrix<Real>(2, 1, {1, 1}));
		CHECK(BandDiagonalSolveDetailed(bandIdentity, badRhs, config).status == AlgorithmStatus::NumericalInstability);
		const BandDiagonalMatrix zero(2, 0, 0, Matrix<Real>(2, 1));
		CHECK(BandDiagonalSolveDetailed(zero, Vector<Real>{1, 1}, config).status == AlgorithmStatus::SingularMatrix);
		const BandDiagonalMatrix nonfinite(2, 0, 0, Matrix<Real>(2, 1, {nan, 1}));
		CHECK(BandDiagonalSolveDetailed(nonfinite, Vector<Real>{1, 1}, config).status == AlgorithmStatus::NumericalInstability);
		CHECK(CholeskySolveDetailed(Matrix<Real>(1, 1, {nan}), Vector<Real>{1}, config).status == AlgorithmStatus::NumericalInstability);
		REQUIRE_THROWS_AS(CholeskySolver<Real>(Matrix<Real>(1, 1, {nan})), MatrixNumericalError);
		config.estimate_error = false;
		const auto completed = LUSolveDetailed(identity, Vector<Real>{1, 2}, config);
		CHECK(completed.status == AlgorithmStatus::Completed);
		CHECK(std::isnan(completed.residual_norm));
		CHECK(std::isnan(GaussJordanSolveDetailed(identity, Vector<Real>{1, 2}, config).residual_norm));
		CHECK(std::isnan(CholeskySolveDetailed(identity, Vector<Real>{1, 2}, config).residual_norm));
		CHECK(std::isnan(BandDiagonalSolveDetailed(bandIdentity, Vector<Real>{1, 2}, config).residual_norm));
		config.estimate_error = true;
		const auto measured = LUSolveDetailed(identity, Vector<Real>{1, 2}, config);
		CHECK(measured.IsSuccess());
		CHECK(measured.residual_norm == 0);
		config.estimate_error = false;
		CHECK(LUSolveDetailed(Matrix<Real>(1, 1, {std::numeric_limits<Real>::min()}),
			Vector<Real>{std::numeric_limits<Real>::max()}, config).status == AlgorithmStatus::NumericalInstability);
	}

	TEST_CASE("L07 condition-one tiny systems retain scale-independent solvability", "[ke4v][scale]")
	{
		const Real scale = Real(1e-20);
		const Matrix<Real> a(2, 2, {scale, 0, 0, scale});
		const Vector<Real> b{scale, -2 * scale};
		SECTION("Cholesky") {
			CholeskySolver<Real> solver(a);
			CHECK(RelativeResidual(a, solver.Solve(b), b) < Precision::EqualityTolerance);
			const auto& l = solver.L();
			CHECK(std::abs(l(0, 0) * l(0, 0) / scale - 1) < Precision::EqualityTolerance);
			const auto reconstruction = l * l.Transposed();
			for (int row = 0; row < 2; ++row)
				for (int col = 0; col < 2; ++col)
					CHECK(std::abs((reconstruction(row, col) - a(row, col)) / scale) < Precision::EqualityTolerance);
		}
		SECTION("QR") {
			QRSolver<Real> solver(a);
			CHECK_FALSE(solver.isSingular());
			CHECK(RelativeResidual(a, solver.Solve(b), b) < Precision::EqualityTolerance);
			const auto reconstruction = solver.GetQ() * solver.GetR();
			for (int row = 0; row < 2; ++row)
				for (int col = 0; col < 2; ++col)
					CHECK(std::abs((reconstruction(row, col) - a(row, col)) / scale) < Precision::EqualityTolerance);
		}
		SECTION("Jacobi") {
			const auto result = JacobiSolver::Solve(a, b);
			CHECK(result.converged);
			CHECK(RelativeResidual(a, result.solution, b) < Precision::EqualityTolerance);
		}
		SECTION("Gauss-Seidel") {
			const auto result = GaussSeidelSolver::Solve(a, b);
			CHECK(result.converged);
			CHECK(RelativeResidual(a, result.solution, b) < Precision::EqualityTolerance);
		}
		SECTION("SOR") {
			const auto result = SORSolver::Solve(a, b, 1);
			CHECK(result.converged);
			CHECK(RelativeResidual(a, result.solution, b) < Precision::EqualityTolerance);
		}
		SECTION("CG optional SPD validation") {
			SparseSolvers::SolverConfig<Real> config;
			config.absTolerance = 0;
			config.validateSPD = true;
			Vector<Real> x(2, Real(0));
			const auto result = SparseSolvers::solveCG(Sparse(a), b, x, config);
			CHECK(result.converged);
			CHECK(RelativeResidual(a, x, b) < Precision::EqualityTolerance);
		}
	}

	TEST_CASE("Sparse norms preserve scalar precision and extreme magnitudes", "[SparseSolvers][norm][scale]")
	{
		static_assert(std::is_same_v<decltype(SparseMatrix::norm2(Vector<double>{})), double>);
		static_assert(std::is_same_v<decltype(SparseMatrix::norm2(std::vector<double>{})), double>);
		for (const double magnitude : {1e-200, 1.0, 1e200}) {
			CAPTURE(magnitude);
			const double expected = std::sqrt(2.0);
			const double dynamic = SparseMatrix::norm2(Vector<double>{magnitude, -magnitude});
			const double standard = SparseMatrix::norm2(std::vector<double>{magnitude, -magnitude});
			REQUIRE(std::isfinite(dynamic));
			REQUIRE(std::isfinite(standard));
			REQUIRE(std::abs(dynamic / magnitude - expected) < PrecisionValues<double>::RoundoffZeroThreshold);
			REQUIRE(std::abs(standard / magnitude - expected) < PrecisionValues<double>::RoundoffZeroThreshold);
		}
		CHECK(SparseMatrix::norm2(Vector<double>{}) == 0);
		CHECK(SparseMatrix::norm2(std::vector<double>{0, 0}) == 0);
		CHECK(SparseMatrix::norm2(Vector<int>{3, 4}) == 5);
		CHECK(SparseMatrix::norm2(std::vector<int>{3, 4}) == 5);
		CHECK(SparseMatrix::norm2(Vector<unsigned>{3, 4}) == 5);
		const int lowest = std::numeric_limits<int>::lowest();
		CHECK(SparseMatrix::norm2(Vector<int>{lowest}) == -double(lowest));
		CHECK(SparseMatrix::norm2(std::vector<int>{lowest}) == -double(lowest));
		const Vector<std::complex<double>> complexValues{{1e200, 1e200}, {1e200, 0}};
		CHECK(std::abs(SparseMatrix::norm2(complexValues) / 1e200 - std::sqrt(3.0)) <
			PrecisionValues<double>::RoundoffZeroThreshold);
		CHECK(std::abs(SparseMatrix::norm2(std::vector<std::complex<double>>(complexValues.begin(), complexValues.end())) /
			1e200 - std::sqrt(3.0)) < PrecisionValues<double>::RoundoffZeroThreshold);
	}

	TEST_CASE("CG retains scale-independent convergence with and without preconditioning", "[SparseSolvers][CG][scale]")
	{
		const Real large = std::sqrt(std::numeric_limits<Real>::max()) * Real(4);
		for (const Real scale : {Real(1) / large, Real(1), large})
			for (bool preconditioned : {false, true}) {
				CAPTURE(scale, preconditioned);
				const Matrix<Real> a(2, 2, {scale, 0, 0, scale});
				const auto sparse = Sparse(a);
				const Vector<Real> b{scale, -2 * scale};
				Vector<Real> x(2, Real(0));
				SparseSolvers::SolverConfig<Real> config;
				config.absTolerance = 0;
				config.tolerance = Precision::EqualityTolerance;
				config.validateSPD = true;
				std::shared_ptr<SparseSolvers::Preconditioner<Real>> preconditioner;
				if (preconditioned) {
					preconditioner = std::make_shared<SparseSolvers::JacobiPreconditioner<Real>>();
					preconditioner->setup(sparse);
				}
				const auto result = preconditioned ? SparseSolvers::solveCG(sparse, b, x, preconditioner, config) :
					SparseSolvers::solveCG(sparse, b, x, config);
				INFO(result.error_message);
				REQUIRE(result.status == AlgorithmStatus::Success);
				REQUIRE(result.converged);
				REQUIRE(result.iterations_used <= 2);
				CHECK(std::abs(x[0] - 1) < Precision::EqualityTolerance);
				CHECK(std::abs(x[1] + 2) < Precision::EqualityTolerance);
				CHECK(RelativeResidual(a, x, b) < Precision::EqualityTolerance);
				CHECK(result.relativeResidual < Precision::EqualityTolerance);
			}
	}

	TEST_CASE("L08 real SVD reconstructs and solves representable extreme scales", "[ke4v][SVD][scale]")
	{
		if (std::numeric_limits<Real>::max_exponent10 < 200)
			SKIP("Double-range witness");
		for (const Real scale : {Real(1e200), Real(1e-200)}) {
			CAPTURE(scale);
			const Matrix<Real> a(2, 2, {scale, 2 * scale, 3 * scale, 4 * scale});
			SVDecompositionSolver<Real> solver(a);
			const auto u = solver.getU(), v = solver.getV();
			const auto w = solver.getW();
			for (int row = 0; row < 2; ++row)
				for (int col = 0; col < 2; ++col) {
					Real reconstruction = 0, orthogonalU = 0, orthogonalV = 0;
					for (int k = 0; k < 2; ++k) {
						reconstruction += u(row, k) * (w[k] / scale) * v(col, k);
						orthogonalU += u(k, row) * u(k, col);
						orthogonalV += v(k, row) * v(k, col);
					}
					CHECK(std::abs(reconstruction - a(row, col) / scale) < Precision::EqualityTolerance);
					CHECK(std::abs(orthogonalU - (row == col ? 1 : 0)) < Precision::EqualityTolerance);
					CHECK(std::abs(orthogonalV - (row == col ? 1 : 0)) < Precision::EqualityTolerance);
				}
			const Vector<Real> b{3 * scale, 7 * scale};
			CHECK(RelativeResidual(a, solver.Solve(b), b) < Precision::EqualityTolerance);
			CHECK(solver.Rank() == 2);
		}
	}

	TEST_CASE("L14 complex SVD retains raw tiny singular values independently of numerical rank", "[ke4v][ComplexSVD][scale]")
	{
		if (std::numeric_limits<Real>::min_exponent10 > -200)
			SKIP("Double-range witness");
		const Real tiny = Real(1e-200);
		// Orthogonal columns with tiny components in both rows exercise scaled column norms.
		const Real c = 1 / std::sqrt(Real(2));
		const Matrix<Complex> a(2, 2, {Complex(c), Complex(0, tiny * c),
			Complex(c), Complex(0, -tiny * c)});
		const auto result = ComplexSVDecompositionSolver::Decompose(a);
		REQUIRE(result.singularValues[1] > 0);
		CHECK(std::abs(result.singularValues[1] / tiny - 1) < Precision::EqualityTolerance);
		CHECK(result.singularValues[1] < Precision::Eps * result.singularValues[0]);
		const auto defaultRank = MatrixAlg::SVDDecompose(a);
		const auto explicitRank = MatrixAlg::SVDDecompose(a, Real(0));
		CHECK(defaultRank.rank == 1);
		CHECK(explicitRank.rank == 2);
		CHECK(defaultRank.singularValues[1] == result.singularValues[1]);
		CHECK(explicitRank.singularValues[1] == result.singularValues[1]);
		for (int row = 0; row < 2; ++row)
			for (int col = 0; col < 2; ++col) {
				Complex reconstructed{}, unitaryU{}, unitaryV{};
				for (int k = 0; k < 2; ++k) {
					reconstructed += result.U(row, k) * result.singularValues[k] * std::conj(result.V(col, k));
					unitaryU += std::conj(result.U(k, row)) * result.U(k, col);
					unitaryV += std::conj(result.V(k, row)) * result.V(k, col);
				}
				CHECK(std::abs((reconstructed - a(row, col)) / std::abs(a(row, col))) < Precision::EqualityTolerance);
				CHECK(std::abs(unitaryU - Complex(row == col ? Real(1) : Real(0))) < Precision::EqualityTolerance);
				CHECK(std::abs(unitaryV - Complex(row == col ? Real(1) : Real(0))) < Precision::EqualityTolerance);
			}
	}

	TEST_CASE("L14 complex SVD diagonal and correlated tiny columns retain representable singular values", "[ke4v][ComplexSVD][scale]")
	{
		if (std::numeric_limits<Real>::min_exponent10 > -200)
			SKIP("Double-range witness");
		const Real tiny = Real(1e-200);
		SECTION("the literal diagonal witness") {
			const Matrix<Complex> a(2, 2, {Complex(1), Complex(0), Complex(0), Complex(tiny)});
			const auto result = ComplexSVDecompositionSolver::Decompose(a);
			CHECK(std::abs(result.singularValues[1] / tiny - 1) < Precision::EqualityTolerance);
		}
		SECTION("nonzero cross-correlation") {
			const Matrix<Complex> a(2, 2, {Complex(1), Complex(tiny), Complex(0), Complex(tiny)});
			const auto result = ComplexSVDecompositionSolver::Decompose(a);
			REQUIRE(result.singularValues[1] > 0);
			CHECK(std::abs(result.singularValues[1] / tiny - 1) < Precision::EqualityTolerance);
			for (int row = 0; row < 2; ++row) {
				Complex reconstructed{};
				for (int k = 0; k < 2; ++k)
					reconstructed += result.U(row, k) * result.singularValues[k] * std::conj(result.V(1, k));
				CHECK(std::abs((reconstructed - a(row, 1)) / tiny) < Precision::EqualityTolerance);
			}
		}
	}
}
