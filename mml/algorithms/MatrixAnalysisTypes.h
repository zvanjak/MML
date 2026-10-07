///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MatrixAnalysisTypes.h                                               ///
///  Description: Scalar traits and result types for matrix analysis                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_MATRIX_ANALYSIS_TYPES_H
#define MML_MATRIX_ANALYSIS_TYPES_H

#include <mml/MMLConcepts.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Vector/Vector.h>

#include <cmath>
#include <complex>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace MML::MatrixAlg
{
	template<MMLScalar Scalar>
	using MatrixMagnitude = decltype(std::abs(std::declval<Scalar>()));

	template<MMLScalar Scalar>
	using MatrixComplexScalar = std::conditional_t<
		MMLComplex<Scalar>,
		std::remove_cvref_t<Scalar>,
		std::complex<std::remove_cvref_t<Scalar>>>;

	template<MMLReal Magnitude>
	struct MatrixComparisonTolerance
	{
		Magnitude absolute{};
		Magnitude relative{};

		bool operator==(const MatrixComparisonTolerance&) const = default;
	};

	enum class MatrixStability
	{
		WellConditioned,
		ModeratelyConditioned,
		IllConditioned,
		Singular
	};

	enum class Definiteness
	{
		PositiveDefinite,
		PositiveSemidefinite,
		NegativeDefinite,
		NegativeSemidefinite,
		Indefinite,
		ZeroSemidefinite
	};

	template<MMLScalar Scalar>
	struct LUDecomposition
	{
		Matrix<Scalar> L;
		Matrix<Scalar> U;
		std::vector<int> permutation;
		Scalar determinant{};
	};

	template<MMLScalar Scalar>
	struct QRDecomposition
	{
		Matrix<Scalar> Q;
		Matrix<Scalar> R;
	};

	template<MMLScalar Scalar>
	struct SVDDecomposition
	{
		Matrix<Scalar> U;
		Vector<MatrixMagnitude<Scalar>> singularValues;
		Matrix<Scalar> V;
		int rank = 0;
		MatrixMagnitude<Scalar> threshold{};
	};

	template<MMLScalar Scalar>
	struct CholeskyDecomposition
	{
		Matrix<Scalar> L;
	};

	template<MMLScalar Scalar>
	struct FundamentalSubspaces
	{
		Matrix<Scalar> columnSpace;
		Matrix<Scalar> rowSpace;
		Matrix<Scalar> nullSpace;
		Matrix<Scalar> leftNullSpace;
		int rank = 0;
		MatrixMagnitude<Scalar> threshold{};
	};

	template<MMLScalar Scalar>
	struct FaddeevLeVerrierResult
	{
		Vector<Scalar> characteristicCoefficients;
		Scalar determinant{};
		std::optional<Matrix<Scalar>> inverse;
	};

	template<MMLScalar Scalar>
	struct EigensystemResult
	{
		Vector<MatrixComplexScalar<Scalar>> eigenvalues;
		Matrix<MatrixComplexScalar<Scalar>> eigenvectors;
		bool converged = false;
		int iterations = 0;
		MatrixMagnitude<Scalar> maxResidual{};
		AlgorithmStatus status = AlgorithmStatus::AlgorithmSpecificFailure;
		std::string algorithmName;
		std::string errorMessage;
	};

	template<MMLScalar Scalar>
	struct SelfAdjointEigensystemResult
	{
		Vector<MatrixMagnitude<Scalar>> eigenvalues;
		Matrix<Scalar> eigenvectors;
		bool converged = false;
		int iterations = 0;
		MatrixMagnitude<Scalar> maxResidual{};
		AlgorithmStatus status = AlgorithmStatus::AlgorithmSpecificFailure;
		std::string algorithmName;
		std::string errorMessage;
	};

	template<MMLScalar Scalar>
	struct MatrixAnalysis
	{
		int rows = 0;
		int cols = 0;
		bool isSquare = true;
		bool isTall = false;
		bool isWide = false;
		bool isUpperTriangular = false;
		bool isLowerTriangular = false;
		bool isDiagonal = false;
		bool isUpperHessenberg = false;
		bool isSymmetric = false;
		bool isSkewSymmetric = false;
		bool isHermitian = false;
		bool isSkewHermitian = false;
		bool isDiagonallyDominant = false;
		MatrixMagnitude<Scalar> sparsity{};
		int rank = 0;
		int nullity = 0;
		std::optional<Scalar> determinant;
		MatrixMagnitude<Scalar> conditionNumber{};
		MatrixStability stability = MatrixStability::Singular;
		std::optional<int> expectedDigitsLost;
		std::optional<Definiteness> definiteness;
		std::string report;
	};
}

#endif // MML_MATRIX_ANALYSIS_TYPES_H
