#include <MML.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEST_CASE("Single header preserves storage boundary checks and sparse aliases", "[single_header][storage-regression]") {
	using namespace MML;
	CHECK_THROWS_AS((SparseMatrix::SparseMatrixCSR<Real>(2, 2, {1}, {0}, {0, 2, 1})), ArgumentError);
	CHECK_THROWS_AS((Matrix3D<Real>(65536, 65536, 1)), MatrixAllocationError);
	SparseMatrix::SparseMatrixCSR<Real> matrix(2, 2, {2, 3}, {1, 0}, {0, 1, 2});
	Vector<Real> input{1, 2};
	matrix.Multiply(input, input);
	CHECK(input[0] == 4);
	CHECK(input[1] == 3);
	CHECK(SparseMatrix::SparseMatrixCSR<Real>().Norm1() == 0);
}

TEST_CASE("Single header verifies scale-safe residuals", "[single_header][verification-regression]") {
	using namespace MML;
	const Real magnitude = std::numeric_limits<Real>::max() * Real(.75);
	const Systems::LinearSystem<Real> system(Matrix<Real>::Identity(4), Vector<Real>(4, magnitude));
	const auto result = system.Verify(Vector<Real>{0, magnitude, magnitude, magnitude});
	CHECK_THAT(result.relativeResidual, Catch::Matchers::WithinRel(Real(.5), Precision::EqualityTolerance));
	CHECK_FALSE(result.isAccurate);
}

TEST_CASE("Single header exposes range-safe factorization and Krylov solving", "[single_header][solver-scale-regression]") {
	using namespace MML;
	const Real large = std::numeric_limits<Real>::max() * Real(.75);
	LUSolver<Real> lu(Matrix<Real>(2, 2, {large, large, -large, large}));
	const auto x = lu.Solve(Vector<Real>{large, 0});
	CHECK_THAT(x[0], Catch::Matchers::WithinRel(Real(.5), Precision::EqualityTolerance));
	CHECK_THAT(x[1], Catch::Matchers::WithinRel(Real(.5), Precision::EqualityTolerance));
	SparseMatrix::SparseMatrixCOO<Real> coo(2, 2);
	coo.addEntry(0, 0, 1);
	coo.addEntry(1, 1, 2);
	SparseMatrix::SparseMatrixCSR<Real> matrix(coo);
	const Real rhsScale = std::sqrt(std::numeric_limits<Real>::max()) / 2;
	Vector<Real> rhs{rhsScale, rhsScale}, solution(2, Real(0));
	SparseSolvers::SolverConfig<Real> config;
	config.absTolerance = 0;
	config.maxIterations = 2;
	const auto result = SparseSolvers::solveGMRES(matrix, rhs, solution, 2, config);
	CHECK(result.converged);
	CHECK_THAT(solution[0] / rhsScale, Catch::Matchers::WithinRel(Real(1), Precision::EqualityTolerance));
	CHECK_THAT(solution[1] / rhsScale, Catch::Matchers::WithinRel(Real(.5), Precision::EqualityTolerance));
}

TEST_CASE("Single header preserves full-range exact arithmetic", "[single_header][exact-integer-regression]") {
	using namespace MML;
	CHECK(NumberTheory::IsPrime(73));
	CHECK(NumberTheory::IsPrime(193));
	CHECK(NumberTheory::Factorize(0).empty());
	CHECK(Combinatorics::Fibonacci(92) == 7540113804746346429LL);
	CHECK(Combinatorics::Catalan(35) == 3116285494907301262LL);
	const auto minimum = std::numeric_limits<long long>::min();
	const auto maximum = std::numeric_limits<long long>::max();
	CHECK(Rational(minimum) / Rational(minimum) == Rational(1));
	CHECK(Rational(maximum, 2) > Rational(maximum - 2, 3));
}
