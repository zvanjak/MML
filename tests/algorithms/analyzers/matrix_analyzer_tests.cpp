#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/algorithms/Analyzers/MatrixAnalyzer.h>
#endif

#include <type_traits>
#include <utility>

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Algorithms::MatrixAnalyzerTests
{
	template<typename Scalar>
	concept SupportsSymmetricAnalyzerEigensystem = requires(const MatrixAnalyzer<Scalar>& analyzer) {
		analyzer.SymmetricEigensystem();
	};

	template<typename Scalar>
	concept SupportsHermitianAnalyzerEigensystem = requires(const MatrixAnalyzer<Scalar>& analyzer) {
		analyzer.HermitianEigensystem();
	};

	template<typename Scalar>
	concept SupportsAnalyzerHessenberg = requires(const MatrixAnalyzer<Scalar>& analyzer) {
		analyzer.ReduceToHessenberg();
	};

	template<typename Scalar>
	concept SupportsLegacyAnalyzerGetEigen = requires(const MatrixAnalyzer<Scalar>& analyzer) {
		analyzer.GetEigen();
	};

	template<typename Scalar>
	concept SupportsLegacyAnalyzerEigenvaluesSymmetric = requires(const MatrixAnalyzer<Scalar>& analyzer) {
		analyzer.EigenvaluesSymmetric();
	};

	TEST_CASE("MatrixAnalyzer exposes scalar-correct Real and Complex APIs", "[MatrixAnalyzer][Complex][Contract][MML2]")
	{
		STATIC_REQUIRE(std::is_copy_constructible_v<MatrixAnalyzer<Complex>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<const MatrixAnalyzer<Complex>&>().Trace()), Complex>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<const MatrixAnalyzer<Complex>&>().ConditionNumber()), Real>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<const MatrixAnalyzer<Complex>&>().Eigenvalues()), Vector<Complex>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<const MatrixAnalyzer<Complex>&>().Eigensystem()),
			const MatrixAlg::EigensystemResult<Complex>&>);
		STATIC_REQUIRE(SupportsSymmetricAnalyzerEigensystem<Real>);
		STATIC_REQUIRE_FALSE(SupportsSymmetricAnalyzerEigensystem<Complex>);
		STATIC_REQUIRE_FALSE(SupportsHermitianAnalyzerEigensystem<Real>);
		STATIC_REQUIRE(SupportsHermitianAnalyzerEigensystem<Complex>);
		STATIC_REQUIRE(SupportsAnalyzerHessenberg<Real>);
		STATIC_REQUIRE_FALSE(SupportsAnalyzerHessenberg<Complex>);
		STATIC_REQUIRE_FALSE(SupportsLegacyAnalyzerGetEigen<Real>);
		STATIC_REQUIRE_FALSE(SupportsLegacyAnalyzerEigenvaluesSymmetric<Real>);
	}

	TEST_CASE("MatrixAnalyzer Complex caches canonical decompositions and eigensystems", "[MatrixAnalyzer][Complex][Cache][MML2]")
	{
		const Matrix<Complex> matrix{2, 2, {
			Complex{4.0, 0.0}, Complex{1.0, 2.0},
			Complex{1.0, -2.0}, Complex{3.0, 0.0}
		}};
		MatrixAnalyzer<Complex> analyzer(matrix);

		REQUIRE(&analyzer.LUDecompose() == &analyzer.LUDecompose());
		REQUIRE(&analyzer.QRDecompose() == &analyzer.QRDecompose());
		REQUIRE(&analyzer.SVDDecompose() == &analyzer.SVDDecompose());
		REQUIRE(&analyzer.CholeskyDecompose() == &analyzer.CholeskyDecompose());
		REQUIRE(&analyzer.Eigensystem() == &analyzer.Eigensystem());
		REQUIRE(&analyzer.HermitianEigensystem() == &analyzer.HermitianEigensystem());
		REQUIRE(&analyzer.Eigensystem(REAL(1e-8), 200) != &analyzer.Eigensystem(REAL(1e-8), 400));
		REQUIRE(&analyzer.CholeskyDecompose(REAL(1e-8)) != &analyzer.CholeskyDecompose(REAL(1e-6)));

		const auto cholesky = analyzer.CholeskyDecompose();
		Matrix<Complex> adjoint(cholesky.L.cols(), cholesky.L.rows());
		for (int row = 0; row < cholesky.L.rows(); ++row)
			for (int col = 0; col < cholesky.L.cols(); ++col)
				adjoint(col, row) = std::conj(cholesky.L(row, col));
		REQUIRE((cholesky.L * adjoint).IsEqualTo(matrix, TOL(1e-9, 1e-4)));
		REQUIRE(analyzer.Eigensystem().algorithmName == "HermitianJacobi");
		REQUIRE(analyzer.HermitianEigenvalues().IsEqualTo(
			Vector<Real>{REAL(1.2087121525220803), REAL(5.79128784747792)}, TOL(1e-9, 1e-4)));
	}

	TEST_CASE("MatrixAnalyzer Complex forwards SVD subspaces and general eigenanalysis", "[MatrixAnalyzer][Complex][Forwarding][MML2]")
	{
		const Matrix<Complex> matrix{2, 3, {
			Complex{1.0, 1.0}, Complex{2.0, 0.0}, Complex{0.0, -1.0},
			Complex{2.0, 2.0}, Complex{4.0, 0.0}, Complex{0.0, -2.0}
		}};
		MatrixAnalyzer<Complex> analyzer(matrix);
		REQUIRE(analyzer.Rank() == 1);
		REQUIRE(analyzer.Nullity() == 2);
		REQUIRE((matrix * analyzer.PseudoInverse() * matrix).IsEqualTo(matrix, TOL(1e-8, 1e-4)));
		REQUIRE(MatrixAlg::FrobeniusNorm(matrix * analyzer.NullSpace()) < TOL(1e-8, 1e-4));

		const Matrix<Complex> general{2, 2, {
			Complex{1.0, 2.0}, Complex{3.0, -1.0},
			Complex{}, Complex{-2.0, 0.5}
		}};
		MatrixAnalyzer<Complex> generalAnalyzer(general);
		REQUIRE(generalAnalyzer.Eigensystem().algorithmName == "ComplexQR");
		REQUIRE(generalAnalyzer.HasComplexEigenvalues());
		REQUIRE(std::abs(generalAnalyzer.SpectralRadius() - std::abs(Complex{1.0, 2.0})) < TOL(1e-9, 1e-4));
	}

	TEST_CASE("MatrixAnalyzer Complex analysis uses Hermitian definiteness", "[MatrixAnalyzer][Complex][Analysis][MML2]")
	{
		MatrixAnalyzer<Complex> hermitian(Matrix<Complex>{2, 2, {
			Complex{4.0, 0.0}, Complex{1.0, 2.0},
			Complex{1.0, -2.0}, Complex{3.0, 0.0}
		}});
		const auto analysis = hermitian.Analyze();
		REQUIRE(analysis.isHermitian);
		REQUIRE_FALSE(analysis.isSymmetric);
		REQUIRE(analysis.definiteness == MatrixAlg::Definiteness::PositiveDefinite);

		MatrixAnalyzer<Complex> nonHermitian(Matrix<Complex>{2, 2, {
			Complex{1.0, 1.0}, Complex{2.0, 0.0},
			Complex{0.0, 1.0}, Complex{3.0, -1.0}
		}});
		REQUIRE_FALSE(nonHermitian.Analyze().definiteness.has_value());
	}

	TEST_CASE("MatrixAnalyzer owns an immutable matrix snapshot", "[MatrixAnalyzer][Ownership][MML2]")
	{
		STATIC_REQUIRE(std::is_copy_constructible_v<MatrixAnalyzer<Real>>);
		STATIC_REQUIRE(std::is_move_constructible_v<MatrixAnalyzer<Real>>);

		Matrix<Real> source{2, 2, {
			REAL(2.0), REAL(0.0),
			REAL(0.0), REAL(3.0)
		}};
		MatrixAnalyzer<Real> analyzer(source);
		source(0, 0) = REAL(100.0);

		REQUIRE(analyzer.GetMatrix()(0, 0) == REAL(2.0));
		REQUIRE(analyzer.Trace() == REAL(5.0));
		REQUIRE(analyzer.Rows() == 2);
		REQUIRE(analyzer.Cols() == 2);
		REQUIRE(analyzer.IsSquare());

		MatrixAnalyzer<Real> moved(std::move(analyzer));
		REQUIRE(moved.GetMatrix()(0, 0) == REAL(2.0));
		REQUIRE(moved.Trace() == REAL(5.0));
	}

	TEST_CASE("MatrixAnalyzer caches decompositions by complete request key", "[MatrixAnalyzer][Cache][MML2]")
	{
		MatrixAnalyzer<Real> analyzer(Matrix<Real>{2, 2, {
			REAL(10.0), REAL(0.0),
			REAL(0.0), REAL(0.1)
		}});

		const auto& automatic = analyzer.SVDDecompose();
		const auto& automaticAgain = analyzer.SVDDecompose();
		const auto& truncated = analyzer.SVDDecompose(REAL(0.5));
		REQUIRE(&automatic == &automaticAgain);
		REQUIRE(&automatic != &truncated);
		REQUIRE(automatic.rank == 2);
		REQUIRE(truncated.rank == 1);
		REQUIRE(&automatic == &analyzer.SVDDecompose());
		for (int index = 1; index <= 32; ++index)
			(void)analyzer.SVDDecompose(REAL(0.01) * index);
		REQUIRE(&automatic == &analyzer.SVDDecompose());

		const auto& spaces = analyzer.FundamentalSubspacesOf();
		const auto& spacesAgain = analyzer.FundamentalSubspacesOf();
		const auto& truncatedSpaces = analyzer.FundamentalSubspacesOf(REAL(0.5));
		REQUIRE(&spaces == &spacesAgain);
		REQUIRE(&spaces != &truncatedSpaces);
		REQUIRE(spaces.rank == 2);
		REQUIRE(truncatedSpaces.rank == 1);

		REQUIRE(&analyzer.LUDecompose() == &analyzer.LUDecompose());
		REQUIRE(&analyzer.QRDecompose() == &analyzer.QRDecompose());
		REQUIRE(&analyzer.CholeskyDecompose() == &analyzer.CholeskyDecompose());
	}

	TEST_CASE("MatrixAnalyzer property caches respect tolerance keys", "[MatrixAnalyzer][Cache][Tolerance][MML2]")
	{
		const Matrix<Real> almostSymmetric{2, 2, {
			REAL(2.0), REAL(1.0),
			REAL(1.0001), REAL(3.0)
		}};
		MatrixAnalyzer<Real> analyzer(almostSymmetric);
		const MatrixAnalyzer<Real>::ComparisonTolerance strict{REAL(1e-8), REAL(0.0)};
		const MatrixAnalyzer<Real>::ComparisonTolerance loose{REAL(1e-3), REAL(0.0)};

		REQUIRE_FALSE(analyzer.IsSymmetric(strict));
		REQUIRE(analyzer.IsSymmetric(loose));
		REQUIRE_FALSE(analyzer.IsSymmetric(strict));

		MatrixAnalyzer<Real> nearZero(Matrix<Real>{2, 2, {
			REAL(5e-5), REAL(0.0),
			REAL(0.0), REAL(2.0)
		}});
		REQUIRE(nearZero.ClassifyDefiniteness(REAL(1e-4)) == MatrixAlg::Definiteness::PositiveSemidefinite);
		REQUIRE(nearZero.IsPositiveSemiDefinite(REAL(1e-4)));
		REQUIRE_FALSE(nearZero.IsNegativeSemiDefinite(REAL(1e-4)));
	}

	TEST_CASE("MatrixAnalyzer caches failed expensive requests", "[MatrixAnalyzer][Cache][Errors][MML2]")
	{
		MatrixAnalyzer<Real> singular(Matrix<Real>{2, 2, {
			REAL(1.0), REAL(2.0),
			REAL(2.0), REAL(4.0)
		}});
		REQUIRE_THROWS_AS(singular.LUDecompose(), SingularMatrixError);
		REQUIRE_THROWS_AS(singular.LUDecompose(), SingularMatrixError);

		MatrixAnalyzer<Real> wide(Matrix<Real>{2, 3});
		REQUIRE_THROWS_AS(wide.QRDecompose(), MatrixDimensionError);
		REQUIRE_THROWS_AS(wide.QRDecompose(), MatrixDimensionError);
	}

	TEST_CASE("MatrixAnalyzer forwards canonical matrix operations", "[MatrixAnalyzer][Forwarding][MML2]")
	{
		const Matrix<Real> matrix{3, 2, {
			REAL(1.0), REAL(2.0),
			REAL(2.0), REAL(4.0),
			REAL(3.0), REAL(6.0)
		}};
		MatrixAnalyzer<Real> analyzer(matrix);

		REQUIRE(analyzer.IsTall());
		REQUIRE(analyzer.Rank() == MatrixAlg::Rank(matrix));
		REQUIRE(analyzer.Nullity() == MatrixAlg::Nullity(matrix));
		REQUIRE(analyzer.ConditionNumber() == MatrixAlg::ConditionNumber(matrix));
		REQUIRE(analyzer.PseudoInverse().IsEqualTo(MatrixAlg::PseudoInverse(matrix), TOL(1e-8, 1e-4)));
		REQUIRE(analyzer.NullSpace().IsEqualTo(MatrixAlg::NullSpace(matrix), TOL(1e-8, 1e-4)));
		REQUIRE(analyzer.ColumnSpace().IsEqualTo(MatrixAlg::ColumnSpace(matrix), TOL(1e-8, 1e-4)));
		REQUIRE(analyzer.RowSpace().IsEqualTo(MatrixAlg::RowSpace(matrix), TOL(1e-8, 1e-4)));
		REQUIRE(analyzer.LeftNullSpace().IsEqualTo(MatrixAlg::LeftNullSpace(matrix), TOL(1e-8, 1e-4)));
	}

	TEST_CASE("MatrixAnalyzer Analyze returns the approved matrix-only report", "[MatrixAnalyzer][Analysis][MML2]")
	{
		MatrixAnalyzer<Real> analyzer(Matrix<Real>{3, 3, {
			REAL(4.0), REAL(1.0), REAL(0.0),
			REAL(1.0), REAL(4.0), REAL(1.0),
			REAL(0.0), REAL(1.0), REAL(4.0)
		}});
		const auto analysis = analyzer.Analyze();
		const auto* cachedSVD = &analyzer.SVDDecompose();
		const auto repeatedAnalysis = analyzer.Analyze();
		REQUIRE(cachedSVD == &analyzer.SVDDecompose());
		REQUIRE(repeatedAnalysis.rank == analysis.rank);

		REQUIRE(analysis.rows == 3);
		REQUIRE(analysis.cols == 3);
		REQUIRE(analysis.isSquare);
		REQUIRE(analysis.isSymmetric);
		REQUIRE(analysis.isHermitian);
		REQUIRE(analysis.isDiagonallyDominant);
		REQUIRE(analysis.rank == 3);
		REQUIRE(analysis.nullity == 0);
		REQUIRE(analysis.determinant.has_value());
		REQUIRE(analysis.stability == MatrixAlg::MatrixStability::WellConditioned);
		REQUIRE(analysis.expectedDigitsLost == std::optional<int>{0});
		REQUIRE(analysis.definiteness == MatrixAlg::Definiteness::PositiveDefinite);
		REQUIRE(analysis.report.find("Matrix: 3x3") != std::string::npos);
		REQUIRE(analysis.report.find("Rank: 3") != std::string::npos);

		const MatrixAnalyzer<Real>::ComparisonTolerance looseSymmetry{REAL(1e-3), REAL(0.0)};
		const auto explicitToleranceAnalysis = analyzer.Analyze(std::nullopt, {}, looseSymmetry, TOL(1e-12, 1e-5));
		REQUIRE(explicitToleranceAnalysis.definiteness == MatrixAlg::Definiteness::PositiveDefinite);

		MatrixAnalyzer<Real> almostSymmetric(Matrix<Real>{2, 2, {
			REAL(5e-5), REAL(1e-5),
			REAL(1.01e-5), REAL(2.0)
		}});
		const MatrixAnalyzer<Real>::ComparisonTolerance strictSymmetry{REAL(1e-9), REAL(0.0)};
		const auto acceptedStrictEigen = almostSymmetric.Analyze(std::nullopt, {}, looseSymmetry, REAL(1e-6));
		const auto acceptedLooseEigen = almostSymmetric.Analyze(std::nullopt, {}, looseSymmetry, REAL(1e-4));
		const auto rejectedSymmetry = almostSymmetric.Analyze(std::nullopt, {}, strictSymmetry, REAL(1e-6));
		const auto acceptedStrictEigenAgain = almostSymmetric.Analyze(std::nullopt, {}, looseSymmetry, REAL(1e-6));

		REQUIRE(acceptedStrictEigen.isSymmetric);
		REQUIRE(acceptedStrictEigen.definiteness == MatrixAlg::Definiteness::PositiveDefinite);
		REQUIRE(acceptedLooseEigen.definiteness == MatrixAlg::Definiteness::PositiveSemidefinite);
		REQUIRE_FALSE(rejectedSymmetry.isSymmetric);
		REQUIRE_FALSE(rejectedSymmetry.definiteness.has_value());
		REQUIRE(acceptedStrictEigenAgain.definiteness == acceptedStrictEigen.definiteness);
	}

	TEST_CASE("MatrixAnalyzer rejects comprehensive analysis of empty matrices", "[MatrixAnalyzer][Errors][MML2]")
	{
		MatrixAnalyzer<Real> analyzer(Matrix<Real>{});
		REQUIRE(analyzer.Rows() == 0);
		REQUIRE(analyzer.Cols() == 0);
		REQUIRE(analyzer.IsSquare());
		REQUIRE_THROWS_AS(analyzer.Analyze(), MatrixDimensionError);
	}
}
