#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Matrix/MatrixSym.h>
#include <mml/algorithms/MatrixAnalysisTypes.h>
#include <mml/algorithms/MatrixAlg.h>
#include <mml/algorithms/Eigen/HermitianMatEigenSolverJacobi.h>
#endif

#include <concepts>
#include <utility>

#include "../test_beds/linear_alg_eq_systems_test_bed.h"

using namespace MML;
using namespace MML::Testing;
using namespace MML::MatrixAlg;

namespace MML::Tests::Algorithms::MatrixAlgTests
{
	template<typename Scalar>
	concept SupportsSVD = requires(const Matrix<Scalar>& matrix) {
		SVDecompositionSolver<Scalar>{matrix};
	};

	template<typename Magnitude>
	concept SupportsMatrixComparisonTolerance = requires {
		typename MatrixComparisonTolerance<Magnitude>;
	};

	template<typename Scalar>
	concept SupportsCanonicalSVD = requires(const Matrix<Scalar>& matrix) {
		SVDDecompose(matrix);
	};

	template<typename Scalar>
	concept SupportsCanonicalCholesky = requires(const Matrix<Scalar>& matrix) {
		CholeskyDecompose(matrix);
	};

	template<MMLScalar Scalar>
	void VerifyMatrixAnalysisResultTypes()
	{
		using Magnitude = MatrixMagnitude<Scalar>;
		using ComplexScalar = MatrixComplexScalar<Scalar>;

		STATIC_REQUIRE(std::same_as<decltype(std::declval<LUDecomposition<Scalar>>().L), Matrix<Scalar>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<LUDecomposition<Scalar>>().U), Matrix<Scalar>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<LUDecomposition<Scalar>>().permutation), std::vector<int>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<LUDecomposition<Scalar>>().determinant), Scalar>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<QRDecomposition<Scalar>>().Q), Matrix<Scalar>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<QRDecomposition<Scalar>>().R), Matrix<Scalar>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<SVDDecomposition<Scalar>>().U), Matrix<Scalar>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<SVDDecomposition<Scalar>>().singularValues), Vector<Magnitude>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<SVDDecomposition<Scalar>>().V), Matrix<Scalar>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<SVDDecomposition<Scalar>>().rank), int>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<SVDDecomposition<Scalar>>().threshold), Magnitude>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<CholeskyDecomposition<Scalar>>().L), Matrix<Scalar>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAlg::FundamentalSubspaces<Scalar>>().columnSpace), Matrix<Scalar>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAlg::FundamentalSubspaces<Scalar>>().rowSpace), Matrix<Scalar>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAlg::FundamentalSubspaces<Scalar>>().nullSpace), Matrix<Scalar>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAlg::FundamentalSubspaces<Scalar>>().leftNullSpace), Matrix<Scalar>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAlg::FundamentalSubspaces<Scalar>>().rank), int>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAlg::FundamentalSubspaces<Scalar>>().threshold), Magnitude>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<FaddeevLeVerrierResult<Scalar>>().characteristicCoefficients), Vector<Scalar>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<FaddeevLeVerrierResult<Scalar>>().determinant), Scalar>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<FaddeevLeVerrierResult<Scalar>>().inverse), std::optional<Matrix<Scalar>>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<EigensystemResult<Scalar>>().eigenvalues), Vector<ComplexScalar>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<EigensystemResult<Scalar>>().eigenvectors), Matrix<ComplexScalar>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<EigensystemResult<Scalar>>().converged), bool>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<EigensystemResult<Scalar>>().iterations), int>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<EigensystemResult<Scalar>>().maxResidual), Magnitude>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<EigensystemResult<Scalar>>().status), AlgorithmStatus>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<EigensystemResult<Scalar>>().algorithmName), std::string>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<EigensystemResult<Scalar>>().errorMessage), std::string>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<SelfAdjointEigensystemResult<Scalar>>().eigenvalues), Vector<Magnitude>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<SelfAdjointEigensystemResult<Scalar>>().eigenvectors), Matrix<Scalar>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<SelfAdjointEigensystemResult<Scalar>>().converged), bool>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<SelfAdjointEigensystemResult<Scalar>>().iterations), int>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<SelfAdjointEigensystemResult<Scalar>>().maxResidual), Magnitude>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<SelfAdjointEigensystemResult<Scalar>>().status), AlgorithmStatus>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<SelfAdjointEigensystemResult<Scalar>>().algorithmName), std::string>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<SelfAdjointEigensystemResult<Scalar>>().errorMessage), std::string>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().rows), int>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().cols), int>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().isSquare), bool>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().isTall), bool>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().isWide), bool>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().isUpperTriangular), bool>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().isLowerTriangular), bool>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().isDiagonal), bool>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().isUpperHessenberg), bool>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().isSymmetric), bool>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().isSkewSymmetric), bool>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().isHermitian), bool>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().isSkewHermitian), bool>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().isDiagonallyDominant), bool>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().sparsity), Magnitude>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().rank), int>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().nullity), int>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().determinant), std::optional<Scalar>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().conditionNumber), Magnitude>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().stability), MatrixStability>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().expectedDigitsLost), std::optional<int>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().definiteness), std::optional<Definiteness>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<MatrixAnalysis<Scalar>>().report), std::string>);
	}

	TEST_CASE("MML 2.0 canonical matrix analysis result types", "[MatrixAlg][Contract][MML2][Types]")
	{
		STATIC_REQUIRE(std::same_as<MatrixMagnitude<float>, float>);
		STATIC_REQUIRE(std::same_as<MatrixMagnitude<double>, double>);
		STATIC_REQUIRE(std::same_as<MatrixMagnitude<long double>, long double>);
		STATIC_REQUIRE(std::same_as<MatrixMagnitude<std::complex<float>>, float>);
		STATIC_REQUIRE(std::same_as<MatrixMagnitude<std::complex<double>>, double>);
		STATIC_REQUIRE(std::same_as<MatrixMagnitude<std::complex<long double>>, long double>);

		STATIC_REQUIRE(std::same_as<MatrixComplexScalar<float>, std::complex<float>>);
		STATIC_REQUIRE(std::same_as<MatrixComplexScalar<double>, std::complex<double>>);
		STATIC_REQUIRE(std::same_as<MatrixComplexScalar<long double>, std::complex<long double>>);
		STATIC_REQUIRE(std::same_as<MatrixComplexScalar<std::complex<double>>, std::complex<double>>);

		STATIC_REQUIRE(SupportsMatrixComparisonTolerance<Real>);
		STATIC_REQUIRE_FALSE(SupportsMatrixComparisonTolerance<Complex>);
		VerifyMatrixAnalysisResultTypes<Real>();
		VerifyMatrixAnalysisResultTypes<Complex>();

		const MatrixAnalysis<Real> analysis;
		REQUIRE(analysis.isSquare);
		REQUIRE(analysis.stability == MatrixStability::Singular);
		REQUIRE_FALSE(analysis.determinant.has_value());
		REQUIRE_FALSE(analysis.expectedDigitsLost.has_value());
		REQUIRE_FALSE(analysis.definiteness.has_value());
		REQUIRE(Definiteness::ZeroSemidefinite != Definiteness::Indefinite);

		const EigensystemResult<Real> eigensystem;
		REQUIRE_FALSE(eigensystem.converged);
		REQUIRE(eigensystem.status == AlgorithmStatus::AlgorithmSpecificFailure);
	}

	TEST_CASE("MML 2.0 current matrix foundation scalar contracts", "[MatrixAlg][Characterization][Contract][MML2]")
	{
		STATIC_REQUIRE(MMLScalar<Real>);
		STATIC_REQUIRE(MMLScalar<Complex>);

		STATIC_REQUIRE(std::same_as<decltype(Trace(std::declval<const Matrix<Real>&>())), Real>);
		STATIC_REQUIRE(std::same_as<decltype(Trace(std::declval<const Matrix<Complex>&>())), Complex>);
		STATIC_REQUIRE(std::same_as<decltype(Determinant(std::declval<const Matrix<Real>&>())), Real>);
		STATIC_REQUIRE(std::same_as<decltype(Determinant(std::declval<const Matrix<Complex>&>())), Complex>);

		STATIC_REQUIRE(std::same_as<decltype(FrobeniusNorm(std::declval<const Matrix<Real>&>())), Real>);
		STATIC_REQUIRE(std::same_as<decltype(FrobeniusNorm(std::declval<const Matrix<Complex>&>())), Real>);
		STATIC_REQUIRE(std::same_as<decltype(OneNorm(std::declval<const Matrix<Complex>&>())), Real>);
		STATIC_REQUIRE(std::same_as<decltype(InfinityNorm(std::declval<const Matrix<Complex>&>())), Real>);

		STATIC_REQUIRE(requires(const Matrix<Complex>& matrix) { LUSolver<Complex>{matrix}; });
		STATIC_REQUIRE(requires(const Matrix<Complex>& matrix) { QRSolver<Complex>{matrix}; });
		STATIC_REQUIRE(SupportsSVD<Real>);
		STATIC_REQUIRE_FALSE(SupportsSVD<Complex>);

		using HermitianResult = HermitianMatEigenSolverJacobi::Result;
		STATIC_REQUIRE(std::same_as<decltype(std::declval<HermitianResult>().eigenvalues), Vector<Real>>);
		STATIC_REQUIRE(std::same_as<decltype(std::declval<HermitianResult>().eigenvectors), Matrix<Complex>>);
	}

	TEST_CASE("MML 2.0 current matrix foundation complex solver contracts", "[MatrixAlg][Characterization][Contract][MML2]")
	{
		const Matrix<Complex> matrix{2, 2, {
			Complex{2.0, 1.0}, Complex{1.0, -1.0},
			Complex{0.0, 2.0}, Complex{3.0, 0.0}
		}};
		const Vector<Complex> expected{Complex{1.0, 1.0}, Complex{-2.0, 0.5}};
		const Vector<Complex> rightHandSide = matrix * expected;

		LUSolver<Complex> lu(matrix);
		QRSolver<Complex> qr(matrix);
		REQUIRE(lu.Solve(rightHandSide).IsEqualTo(expected, TOL(1e-10, 1e-5)));
		REQUIRE(qr.Solve(rightHandSide).IsEqualTo(expected, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("MML 2.0 current matrix foundation singular value contracts", "[MatrixAlg][Characterization][Contract][MML2]")
	{
		const Matrix<Real> singularReal{2, 2, {REAL(1.0), REAL(0.0), REAL(0.0), REAL(0.0)}};
		const Matrix<Complex> singularComplex{2, 2, {Complex{1.0, 1.0}, Complex{}, Complex{}, Complex{}}};

		REQUIRE(Determinant(singularReal) == REAL(0.0));
		REQUIRE(Determinant(singularComplex) == Complex{});
		REQUIRE(std::isinf(ConditionNumber(singularReal)));
	}

	TEST_CASE("MML 2.0 MatrixAlg shape and scale-aware structure", "[MatrixAlg][Structure][MML2]")
	{
		const Matrix<Real> empty;
		REQUIRE(Rows(empty) == 0);
		REQUIRE(Cols(empty) == 0);
		REQUIRE(IsSquare(empty));
		REQUIRE_FALSE(IsTall(empty));
		REQUIRE_FALSE(IsWide(empty));

		const Matrix<Real> tall{3, 2, {
			REAL(1.0), REAL(2.0),
			REAL(0.0), REAL(3.0),
			REAL(0.0), REAL(0.0)
		}};
		REQUIRE(IsTall(tall));
		REQUIRE(IsUpperTriangular(tall));
		REQUIRE_FALSE(IsLowerTriangular(tall));
		REQUIRE_FALSE(IsSymmetric(tall));

		const Matrix<Real> wideLower{2, 3, {
			REAL(1.0), REAL(0.0), REAL(0.0),
			REAL(2.0), REAL(3.0), REAL(0.0)
		}};
		REQUIRE(IsWide(wideLower));
		REQUIRE(IsLowerTriangular(wideLower));
		REQUIRE(IsDiagonal(Matrix<Real>{2, 3, {
			REAL(1.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(2.0), REAL(0.0)
		}}));

		const Matrix<Real> scaled{2, 2, {
			REAL(1e8), REAL(2.0),
			REAL(0.005), REAL(1e8)
		}};
		const MatrixComparisonTolerance<Real> tolerance{REAL(1e-10), REAL(1e-10)};
		REQUIRE(IsUpperTriangular(scaled, tolerance));
		REQUIRE_FALSE(IsUpperTriangular(scaled, MatrixComparisonTolerance<Real>{REAL(1e-10), REAL(0.0)}));
		REQUIRE_THROWS_AS(IsDiagonal(scaled, MatrixComparisonTolerance<Real>{REAL(-1.0), REAL(0.0)}), DomainError);
	}

	TEST_CASE("MML 2.0 MatrixAlg distinguishes complex symmetry families", "[MatrixAlg][Structure][Complex][MML2]")
	{
		const Matrix<Complex> hermitian{2, 2, {
			Complex{2.0, 0.0}, Complex{0.0, 1.0},
			Complex{0.0, -1.0}, Complex{3.0, 0.0}
		}};
		REQUIRE(IsHermitian(hermitian));
		REQUIRE_FALSE(IsSymmetric(hermitian));

		const Matrix<Complex> symmetric{2, 2, {
			Complex{1.0, 0.0}, Complex{0.0, 1.0},
			Complex{0.0, 1.0}, Complex{2.0, 0.0}
		}};
		REQUIRE(IsSymmetric(symmetric));
		REQUIRE_FALSE(IsHermitian(symmetric));

		const Matrix<Complex> skewHermitian{2, 2, {
			Complex{0.0, 1.0}, Complex{1.0, 2.0},
			Complex{-1.0, 2.0}, Complex{0.0, -3.0}
		}};
		REQUIRE(IsSkewHermitian(skewHermitian));
		REQUIRE_FALSE(IsSkewSymmetric(skewHermitian));
	}

	TEST_CASE("MML 2.0 MatrixAlg scalar quantities preserve scalar domains", "[MatrixAlg][Scalar][Complex][MML2]")
	{
		const Matrix<Complex> matrix{2, 2, {
			Complex{3.0, 4.0}, Complex{},
			Complex{}, Complex{1.0, -2.0}
		}};

		STATIC_REQUIRE(std::same_as<decltype(Trace(matrix)), Complex>);
		STATIC_REQUIRE(std::same_as<decltype(FrobeniusNorm(matrix)), Real>);
		STATIC_REQUIRE(std::same_as<decltype(OneNorm(matrix)), Real>);
		STATIC_REQUIRE(std::same_as<decltype(InfinityNorm(matrix)), Real>);
		STATIC_REQUIRE(std::same_as<decltype(Sparsity(matrix)), Real>);

		REQUIRE(Trace(matrix) == Complex{4.0, 2.0});
		REQUIRE(std::abs(FrobeniusNorm(matrix) - std::sqrt(REAL(30.0))) < TOL(1e-12, 1e-5));
		REQUIRE(std::abs(OneNorm(matrix) - REAL(5.0)) < TOL(1e-12, 1e-5));
		REQUIRE(std::abs(InfinityNorm(matrix) - REAL(5.0)) < TOL(1e-12, 1e-5));
		REQUIRE(std::abs(Sparsity(matrix) - REAL(0.5)) < TOL(1e-12, 1e-5));

		const Matrix<Complex> dominant{2, 2, {
			Complex{5.0, 0.0}, Complex{1.0, 1.0},
			Complex{0.5, 0.0}, Complex{4.0, 0.0}
		}};
		REQUIRE(IsDiagonallyDominant(dominant));
		REQUIRE_FALSE(IsDiagonallyDominant(Matrix<Real>{2, 2, {
			REAL(1.0), REAL(1.0),
			REAL(0.0), REAL(1.0)
		}}));
		REQUIRE_THROWS_AS(Sparsity(matrix, REAL(-1.0)), DomainError);
		REQUIRE_THROWS_AS(Trace(Matrix<Real>{2, 3}), MatrixDimensionError);
		REQUIRE_THROWS_AS(FrobeniusNorm(Matrix<Real>{}), MatrixDimensionError);

		const Real large = std::numeric_limits<Real>::max() / REAL(4.0);
		const Real largeNorm = FrobeniusNorm(Matrix<Real>{1, 2, {large, large}});
		REQUIRE(std::isfinite(largeNorm));
		REQUIRE(std::abs(largeNorm / large - std::sqrt(REAL(2.0))) < TOL(1e-12, 1e-5));
	}

	TEST_CASE("MML 2.0 MatrixAlg orthogonal and unitary predicates", "[MatrixAlg][Structure][Complex][MML2]")
	{
		const Matrix<Real> rotation{2, 2, {
			REAL(0.0), REAL(-1.0),
			REAL(1.0), REAL(0.0)
		}};
		REQUIRE(IsOrthogonal(rotation));
		REQUIRE(IsUnitary(rotation));

		const Real inverseSqrt2 = REAL(1.0) / std::sqrt(REAL(2.0));
		const Matrix<Complex> unitary{2, 2, {
			Complex{inverseSqrt2, 0.0}, Complex{0.0, inverseSqrt2},
			Complex{0.0, inverseSqrt2}, Complex{inverseSqrt2, 0.0}
		}};
		REQUIRE(IsUnitary(unitary));
		REQUIRE_FALSE(IsOrthogonal(unitary));
		REQUIRE_FALSE(IsUnitary(Matrix<Complex>{2, 3}));
	}

	TEST_CASE("MML 2.0 MatrixAlg algebraic properties and direct operations", "[MatrixAlg][Algebraic][MML2]")
	{
		const Matrix<Real> nilpotent{2, 2, {
			REAL(0.0), REAL(1.0),
			REAL(0.0), REAL(0.0)
		}};
		REQUIRE(IsNilpotent(nilpotent));
		REQUIRE(IsNilpotent(Matrix<Real>{1, 1, {REAL(0.0)}}));
		REQUIRE_FALSE(IsNilpotent(Matrix<Real>{1, 1, {REAL(1.0)}}));

		const Matrix<Complex> unipotent{2, 2, {
			Complex{1.0, 0.0}, Complex{0.0, 1.0},
			Complex{}, Complex{1.0, 0.0}
		}};
		REQUIRE(IsUnipotent(unipotent));

		const Matrix<Complex> diagonal{2, 2, {
			Complex{2.0, 1.0}, Complex{},
			Complex{}, Complex{3.0, -1.0}
		}};
		STATIC_REQUIRE(std::same_as<decltype(Determinant(diagonal)), Complex>);
		REQUIRE(Determinant(diagonal) == Complex{7.0, 1.0});
		REQUIRE(Determinant(Matrix<Real>{2, 2, {REAL(1.0), REAL(2.0), REAL(2.0), REAL(4.0)}}) == REAL(0.0));
		REQUIRE((diagonal * Inverse(diagonal)).IsEqualTo(Matrix<Complex>::Identity(2), TOL(1e-12, 1e-5)));
		REQUIRE_THROWS_AS(Inverse(Matrix<Real>{2, 2, {REAL(1.0), REAL(2.0), REAL(2.0), REAL(4.0)}}), SingularMatrixError);

		const Matrix<Real> badlyScaled{2, 2, {
			REAL(1e20), REAL(0.0),
			REAL(0.0), REAL(1.0)
		}};
		REQUIRE(std::abs(Determinant(badlyScaled) / REAL(1e20) - REAL(1.0)) < TOL(1e-12, 1e-5));
		REQUIRE((badlyScaled * Inverse(badlyScaled)).IsEqualTo(Matrix<Real>::Identity(2), TOL(1e-12, 1e-5)));

		const Matrix<Real> rankOne{2, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(2.0), REAL(4.0), REAL(6.0)
		}};
		REQUIRE(RankGaussian(rankOne) == 1);
		REQUIRE(RankGaussian(Matrix<Complex>{2, 2, {
			Complex{1.0, 1.0}, Complex{2.0, 0.0},
			Complex{2.0, 2.0}, Complex{4.0, 0.0}
		}}) == 1);

		const auto result = FaddeevLeVerrier(Matrix<Real>{2, 2, {
			REAL(2.0), REAL(0.0),
			REAL(0.0), REAL(3.0)
		}});
		REQUIRE(result.characteristicCoefficients.IsEqualTo(Vector<Real>{REAL(1.0), REAL(-5.0), REAL(6.0)}, TOL(1e-12, 1e-5)));
		REQUIRE(result.determinant == REAL(6.0));
		REQUIRE(result.inverse.has_value());
		REQUIRE(std::abs((*result.inverse)(0, 0) - REAL(0.5)) < TOL(1e-12, 1e-5));
		REQUIRE(std::abs((*result.inverse)(1, 1) - REAL(1.0 / 3.0)) < TOL(1e-12, 1e-5));

		const auto singularResult = FaddeevLeVerrier(Matrix<Real>{2, 2, {
			REAL(1.0), REAL(2.0),
			REAL(2.0), REAL(4.0)
		}});
		REQUIRE(singularResult.determinant == REAL(0.0));
		REQUIRE_FALSE(singularResult.inverse.has_value());

		const auto complexSingularResult = FaddeevLeVerrier(Matrix<Complex>{2, 2, {
			Complex{1.0, 1.0}, Complex{2.0, 0.0},
			Complex{2.0, 2.0}, Complex{4.0, 0.0}
		}});
		REQUIRE_FALSE(complexSingularResult.inverse.has_value());

		const Complex scalar{2.0, -1.0};
		const auto scalarResult = FaddeevLeVerrier(Matrix<Complex>{1, 1, {scalar}});
		REQUIRE(scalarResult.characteristicCoefficients.IsEqualTo(Vector<Complex>{Complex{1.0, 0.0}, -scalar}, TOL(1e-12, 1e-5)));
		REQUIRE(scalarResult.determinant == scalar);
		REQUIRE(scalarResult.inverse.has_value());
	}

	TEST_CASE("MML 2.0 MatrixAlg LU QR and Cholesky decompositions", "[MatrixAlg][Decomposition][MML2]")
	{
		const Matrix<Real> square{3, 3, {
			REAL(0.0), REAL(2.0), REAL(1.0),
			REAL(1.0), REAL(1.0), REAL(0.0),
			REAL(2.0), REAL(0.0), REAL(1.0)
		}};
		const auto lu = LUDecompose(square);
		Matrix<Real> permuted(square.rows(), square.cols());
		for (int row = 0; row < square.rows(); ++row)
			for (int col = 0; col < square.cols(); ++col)
				permuted(row, col) = square(lu.permutation[row], col);
		REQUIRE((lu.L * lu.U).IsEqualTo(permuted, TOL(1e-10, 1e-5)));
		REQUIRE(std::abs(lu.determinant - Determinant(square)) < TOL(1e-10, 1e-5));

		const Matrix<Real> tall{4, 2, {
			REAL(1.0), REAL(2.0),
			REAL(2.0), REAL(0.0),
			REAL(0.0), REAL(3.0),
			REAL(1.0), REAL(1.0)
		}};
		const auto qr = QRDecompose(tall);
		REQUIRE(qr.Q.rows() == 4);
		REQUIRE(qr.Q.cols() == 2);
		REQUIRE(qr.R.rows() == 2);
		REQUIRE(qr.R.cols() == 2);
		REQUIRE((qr.Q * qr.R).IsEqualTo(tall, TOL(1e-9, 1e-4)));
		REQUIRE(IsOrthogonal(qr.Q.transpose() * qr.Q));
		REQUIRE_THROWS_AS(QRDecompose(Matrix<Real>{2, 3}), MatrixDimensionError);

		const Matrix<Real> spd{2, 2, {
			REAL(4.0), REAL(2.0),
			REAL(2.0), REAL(3.0)
		}};
		const auto cholesky = CholeskyDecompose(spd);
		REQUIRE((cholesky.L * cholesky.L.transpose()).IsEqualTo(spd, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("MML 2.0 MatrixAlg full SVD handles tall and wide matrices", "[MatrixAlg][Decomposition][SVD][MML2]")
	{
		STATIC_REQUIRE(SupportsCanonicalSVD<Real>);
		STATIC_REQUIRE(SupportsCanonicalSVD<Complex>);
		STATIC_REQUIRE(SupportsCanonicalCholesky<Real>);
		STATIC_REQUIRE(SupportsCanonicalCholesky<Complex>);

		const auto verify = [](const Matrix<Real>& matrix) {
			const auto svd = SVDDecompose(matrix);
			const int minDimension = std::min(matrix.rows(), matrix.cols());

			REQUIRE(svd.U.rows() == matrix.rows());
			REQUIRE(svd.U.cols() == matrix.rows());
			REQUIRE(svd.V.rows() == matrix.cols());
			REQUIRE(svd.V.cols() == matrix.cols());
			REQUIRE(svd.singularValues.size() == minDimension);
			REQUIRE(IsOrthogonal(svd.U));
			REQUIRE(IsOrthogonal(svd.V));

			Matrix<Real> sigma(matrix.rows(), matrix.cols());
			for (int index = 0; index < minDimension; ++index)
				sigma(index, index) = svd.singularValues[index];
			REQUIRE((svd.U * sigma * svd.V.transpose()).IsEqualTo(matrix, TOL(1e-8, 1e-4)));
		};

		verify(Matrix<Real>{4, 2, {
			REAL(1.0), REAL(2.0),
			REAL(2.0), REAL(0.0),
			REAL(0.0), REAL(3.0),
			REAL(1.0), REAL(1.0)
		}});
		verify(Matrix<Real>{2, 4, {
			REAL(1.0), REAL(2.0), REAL(0.0), REAL(1.0),
			REAL(0.0), REAL(1.0), REAL(3.0), REAL(1.0)
		}});
	}

	TEST_CASE("MML 2.0 MatrixAlg shares one SVD threshold", "[MatrixAlg][SVD][Conditioning][MML2]")
	{
		const Matrix<Real> diagonal{2, 2, {
			REAL(10.0), REAL(0.0),
			REAL(0.0), REAL(0.1)
		}};
		const auto automatic = SVDDecompose(diagonal);
		REQUIRE(automatic.rank == 2);
		REQUIRE(Rank(diagonal) == 2);
		REQUIRE(Nullity(diagonal) == 0);
		REQUIRE(std::abs(ConditionNumber(diagonal) - REAL(100.0)) < TOL(1e-9, 1e-4));
		REQUIRE(std::abs(ConditionNumber1(diagonal) - REAL(100.0)) < TOL(1e-9, 1e-4));
		REQUIRE(std::abs(ConditionNumberInfinity(diagonal) - REAL(100.0)) < TOL(1e-9, 1e-4));

		const std::optional<Real> threshold{REAL(0.5)};
		const auto truncated = SVDDecompose(diagonal, threshold);
		REQUIRE(truncated.threshold == *threshold);
		REQUIRE(truncated.rank == 1);
		REQUIRE(Rank(diagonal, threshold) == 1);
		REQUIRE(Nullity(diagonal, threshold) == 1);
		REQUIRE(std::isinf(ConditionNumber(diagonal, threshold)));
		REQUIRE_THROWS_AS(SVDDecompose(diagonal, REAL(-1.0)), DomainError);
		REQUIRE_THROWS_AS(SVDDecompose(diagonal, std::numeric_limits<Real>::quiet_NaN()), DomainError);
		REQUIRE_THROWS_AS(SVDDecompose(diagonal, std::numeric_limits<Real>::infinity()), DomainError);

		REQUIRE(AssessStability(Matrix<Real>::Identity(2)) == MatrixStability::WellConditioned);
		REQUIRE(AssessStability(Matrix<Real>{2, 2, {REAL(1.0), REAL(0.0), REAL(0.0), REAL(0.0)}}) == MatrixStability::Singular);
		REQUIRE(ExpectedDigitsLost(Matrix<Real>::Identity(2)) == std::optional<int>{0});
		REQUIRE_FALSE(ExpectedDigitsLost(Matrix<Real>{2, 2, {REAL(1.0), REAL(0.0), REAL(0.0), REAL(0.0)}}).has_value());

		const Matrix<Real> rectangularIdentity{3, 2, {
			REAL(1.0), REAL(0.0),
			REAL(0.0), REAL(1.0),
			REAL(0.0), REAL(0.0)
		}};
		REQUIRE(std::abs(ConditionNumber1(rectangularIdentity) - REAL(1.0)) < TOL(1e-9, 1e-4));
		REQUIRE(std::abs(ConditionNumberInfinity(rectangularIdentity) - REAL(1.0)) < TOL(1e-9, 1e-4));
	}

	TEST_CASE("MML 2.0 MatrixAlg pseudoinverse and fundamental subspaces", "[MatrixAlg][SVD][Subspaces][MML2]")
	{
		const Matrix<Real> matrix{2, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(2.0), REAL(4.0), REAL(6.0)
		}};
		const Matrix<Real> pseudoinverse = PseudoInverse(matrix);
		REQUIRE((matrix * pseudoinverse * matrix).IsEqualTo(matrix, TOL(1e-8, 1e-4)));
		REQUIRE((pseudoinverse * matrix * pseudoinverse).IsEqualTo(pseudoinverse, TOL(1e-8, 1e-4)));

		const auto spaces = MatrixAlg::FundamentalSubspacesOf(matrix);
		REQUIRE(spaces.rank == 1);
		REQUIRE(spaces.columnSpace.rows() == 2);
		REQUIRE(spaces.columnSpace.cols() == 1);
		REQUIRE(spaces.rowSpace.rows() == 3);
		REQUIRE(spaces.rowSpace.cols() == 1);
		REQUIRE(spaces.nullSpace.rows() == 3);
		REQUIRE(spaces.nullSpace.cols() == 2);
		REQUIRE(spaces.leftNullSpace.rows() == 2);
		REQUIRE(spaces.leftNullSpace.cols() == 1);
		REQUIRE(FrobeniusNorm(matrix * spaces.nullSpace) < TOL(1e-8, 1e-4));
		REQUIRE(FrobeniusNorm(matrix.transpose() * spaces.leftNullSpace) < TOL(1e-8, 1e-4));
		REQUIRE((spaces.columnSpace.transpose() * spaces.columnSpace).IsEqualTo(Matrix<Real>::Identity(spaces.rank), TOL(1e-8, 1e-4)));
		REQUIRE((spaces.rowSpace.transpose() * spaces.rowSpace).IsEqualTo(Matrix<Real>::Identity(spaces.rank), TOL(1e-8, 1e-4)));
		REQUIRE(NullSpace(matrix).IsEqualTo(spaces.nullSpace, TOL(1e-8, 1e-4)));
		REQUIRE(ColumnSpace(matrix).IsEqualTo(spaces.columnSpace, TOL(1e-8, 1e-4)));
		REQUIRE(RowSpace(matrix).IsEqualTo(spaces.rowSpace, TOL(1e-8, 1e-4)));
		REQUIRE(LeftNullSpace(matrix).IsEqualTo(spaces.leftNullSpace, TOL(1e-8, 1e-4)));

		const Matrix<Real> identity = Matrix<Real>::Identity(3);
		const auto fullRankSpaces = MatrixAlg::FundamentalSubspacesOf(identity);
		REQUIRE(fullRankSpaces.nullSpace.rows() == 3);
		REQUIRE(fullRankSpaces.nullSpace.cols() == 0);
		REQUIRE(fullRankSpaces.leftNullSpace.rows() == 3);
		REQUIRE(fullRankSpaces.leftNullSpace.cols() == 0);
		const Matrix<Real> zeroProduct = identity * fullRankSpaces.nullSpace;
		REQUIRE(zeroProduct.rows() == 3);
		REQUIRE(zeroProduct.cols() == 0);
	}

	TEST_CASE("LinearAlgEqTestBed exposes condition metadata filters", "[MatrixAlg][TestBeds]")
	{
		const auto containsSystem = [](const auto& systems, const std::string& systemName) {
			for (const auto& system : systems)
				if (system.first == systemName)
					return true;
			return false;
		};

		const auto wellConditioned = TestBeds::LinearAlgEqTestBed::getWellConditioned();
		const auto moderatelyConditioned = TestBeds::LinearAlgEqTestBed::getModeratelyConditioned();
		const auto illConditioned = TestBeds::LinearAlgEqTestBed::getIllConditioned();

		REQUIRE(containsSystem(wellConditioned, "mat_3x3"));
		REQUIRE(containsSystem(wellConditioned, "diag_dominant_4x4"));
		REQUIRE(containsSystem(moderatelyConditioned, "hilbert_3x3"));
		REQUIRE(containsSystem(moderatelyConditioned, "frank_4x4"));
		REQUIRE(containsSystem(illConditioned, "hilbert_8x8"));
		REQUIRE(containsSystem(illConditioned, "vandermonde_5x5"));
	}

	/*********************************************************************/
	/*****               MATRIX SYMMETRY TESTS                       *****/
	/*********************************************************************/
	
	TEST_CASE("IsSymmetric_SymmetricMatrix", "[MatrixAlg][Symmetry]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A{3, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(2.0), REAL(4.0), REAL(5.0),
			REAL(3.0), REAL(5.0), REAL(6.0)
		}};
		
		REQUIRE(IsSymmetric(A));
	}
	
	TEST_CASE("IsSymmetric_NonSymmetricMatrix", "[MatrixAlg][Symmetry]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A{3, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(4.0), REAL(5.0), REAL(6.0),  // A[1][0] = 4 != 2 = A[0][1]
			REAL(7.0), REAL(8.0), REAL(9.0)
		}};
		
		REQUIRE_FALSE(IsSymmetric(A));
	}
	
	TEST_CASE("IsSkewSymmetric_True", "[MatrixAlg][Symmetry]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A{3, 3, {
			 REAL(0.0),  REAL(2.0), -REAL(3.0),
			-REAL(2.0),  REAL(0.0),  REAL(5.0),
			 REAL(3.0), -REAL(5.0),  REAL(0.0)
		}};
		
		REQUIRE(IsSkewSymmetric(A));
	}
	
	TEST_CASE("IsSkewSymmetric_NonZeroDiagonal", "[MatrixAlg][Symmetry]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A{3, 3, {
			 REAL(1.0),  REAL(2.0), -REAL(3.0),  // Non-zero diagonal
			-REAL(2.0),  REAL(0.0),  REAL(5.0),
			 REAL(3.0), -REAL(5.0),  REAL(0.0)
		}};
		
		REQUIRE_FALSE(IsSkewSymmetric(A));
	}
	
	/*********************************************************************/
	/*****               MATRIX DEFINITENESS TESTS                   *****/
	/*********************************************************************/
	
	TEST_CASE("IsPositiveDefinite_SPD_3x3", "[MatrixAlg][Definiteness]")
	{
			TEST_PRECISION_INFO();
		// SPD matrix from test bed
		auto sys = TestBeds::spd_3x3();
		
		REQUIRE(IsPositiveDefinite(sys._mat));
		REQUIRE(IsPositiveSemiDefinite(sys._mat));
		REQUIRE_FALSE(IsNegativeDefinite(sys._mat));
		REQUIRE_FALSE(IsIndefinite(sys._mat));
	}
	
	TEST_CASE("IsPositiveDefinite_SPD_CorrelationMatrix", "[MatrixAlg][Definiteness]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::spd_4x4_correlation();
		
		REQUIRE(IsPositiveDefinite(sys._mat));
	}
	
	TEST_CASE("IsPositiveDefinite_SPD_MassMatrix", "[MatrixAlg][Definiteness]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::spd_5x5_mass_matrix();
		
		REQUIRE(IsPositiveDefinite(sys._mat));
	}

	TEST_CASE("IsPositiveDefinite_respects_strict_tolerance", "[MatrixAlg][Definiteness]")
	{
		TEST_PRECISION_INFO();
		const Real tol = TOL(1e-8, 1e-4);
		MatrixSym<Real> A{2, {tol * REAL(0.5), REAL(0.0), REAL(2.0)}};

		REQUIRE_FALSE(IsPositiveDefinite(A, tol));
		REQUIRE(IsPositiveDefinite(A, tol * REAL(0.1)));
	}
	
	TEST_CASE("IsPositiveSemiDefinite_GraphLaplacian", "[MatrixAlg][Definiteness]")
	{
			TEST_PRECISION_INFO();
		// Graph Laplacian is PSD (has zero eigenvalue for constant vector)
		auto sys = TestBeds::spd_6x6_graph_laplacian();
		
		REQUIRE(IsPositiveSemiDefinite(sys._mat));
		// Graph Laplacians have a zero eigenvalue, so not strictly positive definite
		// (though with floating point, it might appear as very small positive)
	}
	
	TEST_CASE("IsNegativeDefinite_NegativeDiagonal", "[MatrixAlg][Definiteness]")
	{
			TEST_PRECISION_INFO();
		MatrixSym<Real> A{3, {
			-REAL(4.0),
			 REAL(0.0), -REAL(5.0),
			 REAL(0.0),  REAL(0.0), -REAL(6.0)
		}};
		
		REQUIRE(IsNegativeDefinite(A));
		REQUIRE(IsNegativeSemiDefinite(A));
		REQUIRE_FALSE(IsPositiveDefinite(A));
	}
	
	TEST_CASE("IsIndefinite_MixedEigenvalues", "[MatrixAlg][Definiteness]")
	{
			TEST_PRECISION_INFO();
		// Indefinite matrix: has both positive and negative eigenvalues
		MatrixSym<Real> A{3, {
			 REAL(2.0),
			 REAL(0.0), -REAL(1.0),
			 REAL(0.0),  REAL(0.0),  REAL(3.0)
		}};
		
		REQUIRE(IsIndefinite(A));
		REQUIRE_FALSE(IsPositiveDefinite(A));
		REQUIRE_FALSE(IsNegativeDefinite(A));
	}
	
	TEST_CASE("ClassifyDefiniteness_All", "[MatrixAlg][Definiteness]")
	{
			TEST_PRECISION_INFO();
		// Positive definite
		MatrixSym<Real> pd{2, {REAL(2.0), REAL(0.5), REAL(2.0)}};
		REQUIRE(ClassifyDefiniteness(pd) == Definiteness::PositiveDefinite);
		
		// Negative definite
		MatrixSym<Real> nd{2, {-REAL(2.0), REAL(0.5), -REAL(2.0)}};
		REQUIRE(ClassifyDefiniteness(nd) == Definiteness::NegativeDefinite);
		
		// Positive semi-definite (singular)
		MatrixSym<Real> psd{2, {REAL(1.0), REAL(1.0), REAL(1.0)}};  // rank 1 matrix
		REQUIRE(ClassifyDefiniteness(psd) == Definiteness::PositiveSemidefinite);
		
		// Indefinite
		MatrixSym<Real> indef{2, {REAL(1.0), REAL(0.0), -REAL(1.0)}};
		REQUIRE(ClassifyDefiniteness(indef) == Definiteness::Indefinite);

		MatrixSym<Real> zero{2, {REAL(0.0), REAL(0.0), REAL(0.0)}};
		REQUIRE(ClassifyDefiniteness(zero) == Definiteness::ZeroSemidefinite);
		REQUIRE(IsPositiveSemiDefinite(zero));
		REQUIRE(IsNegativeSemiDefinite(zero));
		REQUIRE_FALSE(IsPositiveDefinite(zero));
		REQUIRE_FALSE(IsNegativeDefinite(zero));
	}
	
	TEST_CASE("Definiteness_GeneralMatrix", "[MatrixAlg][Definiteness]")
	{
			TEST_PRECISION_INFO();
		// Symmetric matrices may use the general storage type.
		Matrix<Real> A{3, 3, {
			REAL(4.0), REAL(1.0), REAL(0.0),
			REAL(1.0), REAL(4.0), REAL(1.0),
			REAL(0.0), REAL(1.0), REAL(4.0)
		}};
		
		REQUIRE(IsPositiveDefinite(A));

		Matrix<Real> nonsymmetric{2, 2, {
			REAL(2.0), REAL(1.0),
			REAL(0.0), REAL(2.0)
		}};
		REQUIRE_THROWS_AS(ClassifyDefiniteness(nonsymmetric), MatrixDimensionError);
		REQUIRE_THROWS_AS(IsPositiveDefinite(nonsymmetric), MatrixDimensionError);
	}
	
	/*********************************************************************/
	/*****                SVD-BASED UTILITIES                        *****/
	/*********************************************************************/
	
	TEST_CASE("SingularValues_DiagonalMatrix", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A{3, 3, {
			REAL(5.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(3.0), REAL(0.0),
			REAL(0.0), REAL(0.0), REAL(1.0)
		}};
		
		Vector<Real> sv = SingularValues(A);
		
		// Singular values should be sorted in descending order
		REQUIRE(sv.size() == 3);
		REQUIRE(std::abs(sv[0] - REAL(5.0)) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(sv[1] - REAL(3.0)) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(sv[2] - REAL(1.0)) < TOL(1e-10, 1e-5));
	}
	
	TEST_CASE("Rank_FullRank", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_4x4();
		
		int rank = RankGaussian(sys._mat);
		REQUIRE(rank == 4);
	}
	
	TEST_CASE("Rank_RankDeficient", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		// Rank 2 matrix (3x3)
		Matrix<Real> A{3, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(2.0), REAL(4.0), REAL(6.0),   // Row 2 = 2 * Row 1
			REAL(3.0), REAL(6.0), REAL(9.0)    // Row 3 = 3 * Row 1
		}};
		
		int rank = RankGaussian(A);
		REQUIRE(rank == 1);
	}
	
	TEST_CASE("Nullity_FullRank", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_4x4();
		
		int nullity = Nullity(sys._mat);
		REQUIRE(nullity == 0);
	}
	
	TEST_CASE("Nullity_RankDeficient", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		// Rank 1 matrix → nullity = 2
		Matrix<Real> A{3, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(2.0), REAL(4.0), REAL(6.0),
			REAL(3.0), REAL(6.0), REAL(9.0)
		}};
		
		int nullity = Nullity(A);
		REQUIRE(nullity == 2);
	}
	
	TEST_CASE("ConditionNumber_WellConditioned", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		// Diagonal matrix with condition number = 5/1 = 5
		Matrix<Real> A{3, 3, {
			REAL(5.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(3.0), REAL(0.0),
			REAL(0.0), REAL(0.0), REAL(1.0)
		}};
		
		Real cond = ConditionNumber(A);
		REQUIRE(std::abs(cond - REAL(5.0)) < TOL(1e-10, 1e-5));
	}
	
	TEST_CASE("ConditionNumber_Identity", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		auto I = Matrix<Real>::Identity(5);
		
		Real cond = ConditionNumber(I);
		REQUIRE(std::abs(cond - REAL(1.0)) < TOL(1e-10, 1e-5));
	}
	
	TEST_CASE("NullSpace_FullRank", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_4x4();
		
		Matrix<Real> ns = NullSpace(sys._mat);
		
		// Full rank matrix has trivial null space
		REQUIRE(ns.cols() == 0);
	}
	
	TEST_CASE("NullSpace_RankDeficient", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		// Rank 1 matrix
		Matrix<Real> A{3, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(2.0), REAL(4.0), REAL(6.0),
			REAL(3.0), REAL(6.0), REAL(9.0)
		}};
		
		Matrix<Real> ns = NullSpace(A);
		
		// Should have 2-dimensional null space
		REQUIRE(ns.cols() == 2);
		REQUIRE(ns.rows() == 3);
		
		// Verify: A * x should be zero for each null space vector
		for (int col = 0; col < ns.cols(); col++)
		{
			Vector<Real> x(3);
			for (int i = 0; i < 3; i++)
				x[i] = ns(i, col);
			
			Vector<Real> Ax = A * x;
			REQUIRE(Ax.NormL2() < TOL(1e-10, 1e-5));
		}
	}
	
	TEST_CASE("ColumnSpace_FullRank", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		auto sys = TestBeds::diag_dominant_4x4();
		
		Matrix<Real> cs = ColumnSpace(sys._mat);
		
		// Full rank 4x4 → column space is R^4
		REQUIRE(cs.cols() == 4);
		REQUIRE(cs.rows() == 4);
	}
	
	TEST_CASE("ColumnSpace_RankDeficient", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		// Rank 1 matrix
		Matrix<Real> A{3, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(2.0), REAL(4.0), REAL(6.0),
			REAL(3.0), REAL(6.0), REAL(9.0)
		}};
		
		Matrix<Real> cs = ColumnSpace(A);
		
		// Column space is 1-dimensional
		REQUIRE(cs.cols() == 1);
		REQUIRE(cs.rows() == 3);
	}
	
	TEST_CASE("FundamentalSubspaces_RankNullityTheorem", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		// For m×n matrix: rank + nullity = n
		Matrix<Real> A{4, 3, {
			REAL(1.0), REAL(2.0), REAL(3.0),
			REAL(2.0), REAL(4.0), REAL(6.0),
			REAL(3.0), REAL(6.0), REAL(9.0),
			REAL(4.0), REAL(8.0), REAL(12.0)
		}};
		
		auto fs = FundamentalSubspacesOf(A);
		
		int m = A.rows();  // 4
		int n = A.cols();  // 3
		
		// Rank-nullity theorem: rank + nullity = n
		REQUIRE(fs.rank + fs.nullSpace.cols() == n);
		
		// dim(col space) = rank
		REQUIRE(fs.columnSpace.cols() == fs.rank);
		
		// dim(row space) = rank
		REQUIRE(fs.rowSpace.cols() == fs.rank);
		
		// dim(left null space) = m - rank
		REQUIRE(fs.leftNullSpace.cols() == m - fs.rank);
	}
	
	TEST_CASE("PseudoInverse_FullRankSquare", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		// For invertible matrix, pseudoinverse = inverse
		Matrix<Real> A{3, 3, {
			REAL(2.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(3.0), REAL(0.0),
			REAL(0.0), REAL(0.0), REAL(4.0)
		}};
		
		Matrix<Real> Ainv = PseudoInverse(A);
		
		// A * A⁺ should be identity
		Matrix<Real> I = A * Ainv;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
			{
				Real expected = (i == j) ? REAL(1.0) : REAL(0.0);
				REQUIRE(std::abs(I(i, j) - expected) < TOL(1e-10, 1e-5));
			}
	}
	
	TEST_CASE("PseudoInverse_LeastSquares", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		// Overdetermined system: pseudoinverse gives least squares solution
		Matrix<Real> A{4, 2, {
			REAL(1.0), REAL(0.0),
			REAL(0.0), REAL(1.0),
			REAL(1.0), REAL(0.0),
			REAL(0.0), REAL(1.0)
		}};
		
		Vector<Real> b{REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)};
		
		Matrix<Real> Aplus = PseudoInverse(A);
		
		// x = A⁺ * b is least squares solution
		Vector<Real> x(2);
		for (int i = 0; i < 2; i++)
		{
			x[i] = REAL(0.0);
			for (int j = 0; j < 4; j++)
				x[i] += Aplus(i, j) * b[j];
		}
		
		// Expected: x = [2, 3] (average of (1,3) and (2,4))
		REQUIRE(std::abs(x[0] - REAL(2.0)) < TOL(1e-10, 1e-5));
		REQUIRE(std::abs(x[1] - REAL(3.0)) < TOL(1e-10, 1e-5));
	}
	
	TEST_CASE("ComputeSVD_Basic", "[MatrixAlg][SVD]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A{3, 3, {
			REAL(1.0), REAL(0.0), REAL(0.0),
			REAL(0.0), REAL(2.0), REAL(0.0),
			REAL(0.0), REAL(0.0), REAL(3.0)
		}};
		
		auto svd = SVDDecompose(A);
		
		REQUIRE(svd.rank == 3);
		REQUIRE(std::abs(ConditionNumber(A) - REAL(3.0)) < TOL(1e-10, 1e-5));
		REQUIRE(svd.singularValues.size() == 3);
	}

} // namespace MML::Tests::Algorithms::MatrixAlgTests
