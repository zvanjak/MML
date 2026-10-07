///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MatrixAnalyzer.h                                                    ///
///  Description: Owning cached facade for canonical matrix analysis                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_MATRIX_ANALYZER_H
#define MML_MATRIX_ANALYZER_H

#include <mml/algorithms/MatrixAlg.h>

#include <algorithm>
#include <deque>
#include <exception>
#include <optional>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace MML
{
	template<MMLScalar Scalar = Real>
	class MatrixAnalyzer
	{
	public:
		using Magnitude = MatrixAlg::MatrixMagnitude<Scalar>;
		using ComparisonTolerance = MatrixAlg::MatrixComparisonTolerance<Magnitude>;
		using Threshold = std::optional<Magnitude>;

		explicit MatrixAnalyzer(const Matrix<Scalar>& matrix) : _matrix(matrix) {}
		explicit MatrixAnalyzer(Matrix<Scalar>&& matrix) : _matrix(std::move(matrix)) {}

		const Matrix<Scalar>& GetMatrix() const noexcept { return _matrix; }

		int Rows() const noexcept { return MatrixAlg::Rows(_matrix); }
		int Cols() const noexcept { return MatrixAlg::Cols(_matrix); }
		bool IsSquare() const noexcept { return MatrixAlg::IsSquare(_matrix); }
		bool IsTall() const noexcept { return MatrixAlg::IsTall(_matrix); }
		bool IsWide() const noexcept { return MatrixAlg::IsWide(_matrix); }

		bool IsUpperHessenberg(ComparisonTolerance tolerance = MatrixAlg::Detail::DefaultDiagonalTolerance<Scalar>()) const {
			return MatrixAlg::IsUpperHessenberg(_matrix, tolerance);
		}
		bool IsUpperTriangular(ComparisonTolerance tolerance = MatrixAlg::Detail::DefaultDiagonalTolerance<Scalar>()) const {
			return MatrixAlg::IsUpperTriangular(_matrix, tolerance);
		}
		bool IsLowerTriangular(ComparisonTolerance tolerance = MatrixAlg::Detail::DefaultDiagonalTolerance<Scalar>()) const {
			return MatrixAlg::IsLowerTriangular(_matrix, tolerance);
		}
		bool IsDiagonal(ComparisonTolerance tolerance = MatrixAlg::Detail::DefaultDiagonalTolerance<Scalar>()) const {
			return MatrixAlg::IsDiagonal(_matrix, tolerance);
		}
		bool IsSymmetric(ComparisonTolerance tolerance = MatrixAlg::Detail::DefaultSymmetryTolerance<Scalar>()) const {
			return Cached(_symmetryCache, tolerance, [&] { return MatrixAlg::IsSymmetric(_matrix, tolerance); });
		}
		bool IsSkewSymmetric(ComparisonTolerance tolerance = MatrixAlg::Detail::DefaultSymmetryTolerance<Scalar>()) const {
			return MatrixAlg::IsSkewSymmetric(_matrix, tolerance);
		}
		bool IsHermitian(ComparisonTolerance tolerance = MatrixAlg::Detail::DefaultSymmetryTolerance<Scalar>()) const {
			return MatrixAlg::IsHermitian(_matrix, tolerance);
		}
		bool IsSkewHermitian(ComparisonTolerance tolerance = MatrixAlg::Detail::DefaultSymmetryTolerance<Scalar>()) const {
			return MatrixAlg::IsSkewHermitian(_matrix, tolerance);
		}
		bool IsDiagonallyDominant(ComparisonTolerance tolerance = MatrixAlg::Detail::DefaultDiagonalTolerance<Scalar>()) const {
			return MatrixAlg::IsDiagonallyDominant(_matrix, tolerance);
		}

		Magnitude Sparsity(Magnitude threshold = PrecisionValues<Magnitude>::EigenSolverZeroThreshold) const {
			return MatrixAlg::Sparsity(_matrix, threshold);
		}
		Scalar Trace() const { return MatrixAlg::Trace(_matrix); }
		Magnitude FrobeniusNorm() const { return MatrixAlg::FrobeniusNorm(_matrix); }
		Magnitude OneNorm() const { return MatrixAlg::OneNorm(_matrix); }
		Magnitude InfinityNorm() const { return MatrixAlg::InfinityNorm(_matrix); }
		bool IsOrthogonal(ComparisonTolerance tolerance = MatrixAlg::Detail::DefaultSymmetryTolerance<Scalar>()) const {
			return MatrixAlg::IsOrthogonal(_matrix, tolerance);
		}
		bool IsUnitary(ComparisonTolerance tolerance = MatrixAlg::Detail::DefaultSymmetryTolerance<Scalar>()) const {
			return MatrixAlg::IsUnitary(_matrix, tolerance);
		}
		bool IsNilpotent(ComparisonTolerance tolerance = MatrixAlg::Detail::DefaultDiagonalTolerance<Scalar>()) const {
			return MatrixAlg::IsNilpotent(_matrix, tolerance);
		}
		bool IsUnipotent(ComparisonTolerance tolerance = MatrixAlg::Detail::DefaultDiagonalTolerance<Scalar>()) const {
			return MatrixAlg::IsUnipotent(_matrix, tolerance);
		}

		Scalar Determinant() const { return MatrixAlg::Determinant(_matrix); }
		Matrix<Scalar> Inverse() const { return MatrixAlg::Inverse(_matrix); }
		int RankGaussian(Magnitude tolerance = static_cast<Magnitude>(Defaults::RankAlgEPS)) const {
			return MatrixAlg::RankGaussian(_matrix, tolerance);
		}
		MatrixAlg::FaddeevLeVerrierResult<Scalar> FaddeevLeVerrier() const {
			return MatrixAlg::FaddeevLeVerrier(_matrix);
		}

		const MatrixAlg::LUDecomposition<Scalar>& LUDecompose() const {
			return Cached(_luCache, true, [&] { return MatrixAlg::LUDecompose(_matrix); });
		}
		const MatrixAlg::QRDecomposition<Scalar>& QRDecompose() const {
			return Cached(_qrCache, true, [&] { return MatrixAlg::QRDecompose(_matrix); });
		}
		const MatrixAlg::CholeskyDecomposition<Scalar>& CholeskyDecompose(
				Magnitude tolerance = PrecisionValues<Magnitude>::IsMatrixSymmetricTolerance) const {
			return Cached(_choleskyCache, tolerance, [&] { return MatrixAlg::CholeskyDecompose(_matrix, tolerance); });
		}
		const MatrixAlg::SVDDecomposition<Scalar>& SVDDecompose(Threshold threshold = std::nullopt) const {
			return Cached(_svdCache, threshold, [&] { return MatrixAlg::SVDDecompose(_matrix, threshold); });
		}

		Vector<Magnitude> SingularValues() const { return SVDDecompose().singularValues; }
		int Rank(Threshold threshold = std::nullopt) const { return SVDDecompose(threshold).rank; }
		int Nullity(Threshold threshold = std::nullopt) const { return Cols() - SVDDecompose(threshold).rank; }
		Magnitude ConditionNumber(Threshold threshold = std::nullopt) const {
			return MatrixAlg::Detail::ConditionNumberFromSVD(SVDDecompose(threshold));
		}
		Magnitude ConditionNumber1(Threshold threshold = std::nullopt) const {
			const auto& svd = SVDDecompose(threshold);
			if (svd.rank < svd.singularValues.size())
				return std::numeric_limits<Magnitude>::infinity();
			return OneNorm() * MatrixAlg::OneNorm(MatrixAlg::Detail::PseudoInverseFromSVD(svd, Rows(), Cols()));
		}
		Magnitude ConditionNumberInfinity(Threshold threshold = std::nullopt) const {
			const auto& svd = SVDDecompose(threshold);
			if (svd.rank < svd.singularValues.size())
				return std::numeric_limits<Magnitude>::infinity();
			return InfinityNorm() * MatrixAlg::InfinityNorm(MatrixAlg::Detail::PseudoInverseFromSVD(svd, Rows(), Cols()));
		}
		MatrixAlg::MatrixStability AssessStability(Threshold threshold = std::nullopt) const {
			return MatrixAlg::Detail::StabilityFromConditionNumber(ConditionNumber(threshold));
		}
		std::optional<int> ExpectedDigitsLost(Threshold threshold = std::nullopt) const {
			return MatrixAlg::Detail::DigitsLostFromConditionNumber(ConditionNumber(threshold));
		}
		Matrix<Scalar> PseudoInverse(Threshold threshold = std::nullopt) const {
			return MatrixAlg::Detail::PseudoInverseFromSVD(SVDDecompose(threshold), Rows(), Cols());
		}

		const MatrixAlg::FundamentalSubspaces<Scalar>& FundamentalSubspacesOf(Threshold threshold = std::nullopt) const {
			return Cached(_subspacesCache, threshold, [&] {
				return MatrixAlg::Detail::FundamentalSubspacesFromSVD(SVDDecompose(threshold), Rows(), Cols());
			});
		}
		Matrix<Scalar> NullSpace(Threshold threshold = std::nullopt) const { return FundamentalSubspacesOf(threshold).nullSpace; }
		Matrix<Scalar> ColumnSpace(Threshold threshold = std::nullopt) const { return FundamentalSubspacesOf(threshold).columnSpace; }
		Matrix<Scalar> RowSpace(Threshold threshold = std::nullopt) const { return FundamentalSubspacesOf(threshold).rowSpace; }
		Matrix<Scalar> LeftNullSpace(Threshold threshold = std::nullopt) const { return FundamentalSubspacesOf(threshold).leftNullSpace; }

		MatrixAlg::Definiteness ClassifyDefiniteness(Magnitude tolerance = PrecisionValues<Magnitude>::EigenSolverConvergenceTolerance) const {
			return Cached(_definitenessCache, tolerance, [&] { return MatrixAlg::ClassifyDefiniteness(_matrix, tolerance); });
		}
		bool IsPositiveDefinite(Magnitude tolerance = PrecisionValues<Magnitude>::EigenSolverConvergenceTolerance) const {
			return ClassifyDefiniteness(tolerance) == MatrixAlg::Definiteness::PositiveDefinite;
		}
		bool IsPositiveDefinite(ComparisonTolerance symmetryTolerance, Magnitude definitenessTolerance) const {
			return IsSelfAdjoint(symmetryTolerance) &&
				   ClassifyAcceptedSelfAdjointDefiniteness(symmetryTolerance, definitenessTolerance) ==
					   MatrixAlg::Definiteness::PositiveDefinite;
		}
		bool IsPositiveSemiDefinite(Magnitude tolerance = PrecisionValues<Magnitude>::EigenSolverConvergenceTolerance) const {
			const auto value = ClassifyDefiniteness(tolerance);
			return value == MatrixAlg::Definiteness::PositiveDefinite ||
				   value == MatrixAlg::Definiteness::PositiveSemidefinite ||
				   value == MatrixAlg::Definiteness::ZeroSemidefinite;
		}
		bool IsNegativeDefinite(Magnitude tolerance = PrecisionValues<Magnitude>::EigenSolverConvergenceTolerance) const {
			return ClassifyDefiniteness(tolerance) == MatrixAlg::Definiteness::NegativeDefinite;
		}
		bool IsNegativeSemiDefinite(Magnitude tolerance = PrecisionValues<Magnitude>::EigenSolverConvergenceTolerance) const {
			const auto value = ClassifyDefiniteness(tolerance);
			return value == MatrixAlg::Definiteness::NegativeDefinite ||
				   value == MatrixAlg::Definiteness::NegativeSemidefinite ||
				   value == MatrixAlg::Definiteness::ZeroSemidefinite;
		}
		bool IsIndefinite(Magnitude tolerance = PrecisionValues<Magnitude>::EigenSolverConvergenceTolerance) const {
			return ClassifyDefiniteness(tolerance) == MatrixAlg::Definiteness::Indefinite;
		}

		MatrixAlg::HessenbergResult ReduceToHessenberg() const requires MMLReal<Scalar> {
			return MatrixAlg::ReduceToHessenberg(_matrix);
		}

		const MatrixAlg::EigensystemResult<Scalar>& Eigensystem(
				Magnitude tolerance = PrecisionValues<Magnitude>::EigenSolverConvergenceTolerance,
				int maxIterations = 1000) const {
			const EigenRequest request{tolerance, maxIterations};
			return Cached(_eigenCache, request, [&] { return MatrixAlg::Eigensystem(_matrix, tolerance, maxIterations); });
		}

		Vector<MatrixAlg::MatrixComplexScalar<Scalar>> Eigenvalues(
				Magnitude tolerance = PrecisionValues<Magnitude>::EigenSolverConvergenceTolerance,
				int maxIterations = 1000) const {
			return Eigensystem(tolerance, maxIterations).eigenvalues;
		}

		const MatrixAlg::SelfAdjointEigensystemResult<Scalar>& SymmetricEigensystem(
				Magnitude tolerance = PrecisionValues<Magnitude>::EigenSolverConvergenceTolerance,
				int maxIterations = 100) const requires MMLReal<Scalar> {
			const EigenRequest request{tolerance, maxIterations};
			return Cached(_selfAdjointEigenCache, request, [&] {
				return MatrixAlg::SymmetricEigensystem(_matrix, tolerance, maxIterations);
			});
		}

		Vector<Magnitude> SymmetricEigenvalues(
				Magnitude tolerance = PrecisionValues<Magnitude>::EigenSolverConvergenceTolerance,
				int maxIterations = 100) const requires MMLReal<Scalar> {
			return SymmetricEigensystem(tolerance, maxIterations).eigenvalues;
		}

		const MatrixAlg::SelfAdjointEigensystemResult<Scalar>& HermitianEigensystem(
				Magnitude tolerance = PrecisionValues<Magnitude>::EigenSolverConvergenceTolerance,
				int maxIterations = 100) const requires MMLComplex<Scalar> {
			const EigenRequest request{tolerance, maxIterations};
			return Cached(_selfAdjointEigenCache, request, [&] {
				return MatrixAlg::HermitianEigensystem(_matrix, tolerance, maxIterations);
			});
		}

		Vector<Magnitude> HermitianEigenvalues(
				Magnitude tolerance = PrecisionValues<Magnitude>::EigenSolverConvergenceTolerance,
				int maxIterations = 100) const requires MMLComplex<Scalar> {
			return HermitianEigensystem(tolerance, maxIterations).eigenvalues;
		}

		Magnitude SpectralRadius(
				Magnitude tolerance = PrecisionValues<Magnitude>::EigenSolverConvergenceTolerance,
				int maxIterations = 1000) const {
			Magnitude maximum{};
			for (const auto& eigenvalue : Eigensystem(tolerance, maxIterations).eigenvalues)
				maximum = std::max(maximum, static_cast<Magnitude>(std::abs(eigenvalue)));
			return maximum;
		}

		bool HasComplexEigenvalues(
				Magnitude tolerance = PrecisionValues<Magnitude>::EigenSolverConvergenceTolerance,
				int maxIterations = 1000) const {
			for (const auto& eigenvalue : Eigensystem(tolerance, maxIterations).eigenvalues)
				if (std::abs(eigenvalue.imag()) > tolerance)
					return true;
			return false;
		}

		MatrixAlg::MatrixAnalysis<Scalar> Analyze(
				Threshold threshold = std::nullopt,
				ComparisonTolerance structuralTolerance = MatrixAlg::Detail::DefaultDiagonalTolerance<Scalar>(),
				ComparisonTolerance symmetryTolerance = MatrixAlg::Detail::DefaultSymmetryTolerance<Scalar>(),
				Magnitude definitenessTolerance = PrecisionValues<Magnitude>::EigenSolverConvergenceTolerance) const {
			const auto& svd = SVDDecompose(threshold);
			MatrixAlg::MatrixAnalysis<Scalar> result;
			result.rows = Rows();
			result.cols = Cols();
			result.isSquare = IsSquare();
			result.isTall = IsTall();
			result.isWide = IsWide();
			result.isUpperTriangular = IsUpperTriangular(structuralTolerance);
			result.isLowerTriangular = IsLowerTriangular(structuralTolerance);
			result.isDiagonal = IsDiagonal(structuralTolerance);
			result.isUpperHessenberg = IsUpperHessenberg(structuralTolerance);
			result.isSymmetric = IsSymmetric(symmetryTolerance);
			result.isSkewSymmetric = IsSkewSymmetric(symmetryTolerance);
			result.isHermitian = IsHermitian(symmetryTolerance);
			result.isSkewHermitian = IsSkewHermitian(symmetryTolerance);
			result.isDiagonallyDominant = IsDiagonallyDominant(structuralTolerance);
			result.sparsity = Sparsity();
			result.rank = svd.rank;
			result.nullity = Cols() - svd.rank;
			if (result.isSquare)
				result.determinant = Determinant();
			result.conditionNumber = MatrixAlg::Detail::ConditionNumberFromSVD(svd);
			result.stability = MatrixAlg::Detail::StabilityFromConditionNumber(result.conditionNumber);
			result.expectedDigitsLost = MatrixAlg::Detail::DigitsLostFromConditionNumber(result.conditionNumber);
			if (IsSelfAdjointResult(result))
				result.definiteness = ClassifyAcceptedSelfAdjointDefiniteness(symmetryTolerance, definitenessTolerance);
			result.report = GenerateReport(result);
			return result;
		}

	private:
		struct DefinitenessRequest
		{
			ComparisonTolerance symmetryTolerance;
			Magnitude eigenvalueTolerance;

			bool operator==(const DefinitenessRequest&) const = default;
		};

		struct EigenRequest
		{
			Magnitude tolerance;
			int maxIterations;

			bool operator==(const EigenRequest&) const = default;
		};

		template<typename Key, typename Value>
		struct CacheEntry
		{
			Key key;
			std::optional<Value> value;
			std::exception_ptr error;
		};

		template<typename Key, typename Value, typename Factory>
		static const Value& Cached(std::deque<CacheEntry<Key, Value>>& cache, const Key& key, Factory&& factory) {
			const auto found = std::find_if(cache.begin(), cache.end(), [&](const auto& entry) { return entry.key == key; });
			if (found != cache.end()) {
				if (found->error)
					std::rethrow_exception(found->error);
				return *found->value;
			}

			try {
				cache.push_back({key, std::forward<Factory>(factory)(), {}});
				return *cache.back().value;
			}
			catch (...) {
				cache.push_back({key, std::nullopt, std::current_exception()});
				std::rethrow_exception(cache.back().error);
			}
		}

		static std::string GenerateReport(const MatrixAlg::MatrixAnalysis<Scalar>& analysis) {
			std::ostringstream report;
			report << "Matrix: " << analysis.rows << 'x' << analysis.cols << '\n';
			report << "Rank: " << analysis.rank << ", nullity: " << analysis.nullity << '\n';
			report << "Condition number: " << analysis.conditionNumber << '\n';
			if (analysis.definiteness.has_value())
				report << "Definiteness: " << static_cast<int>(*analysis.definiteness) << '\n';
			return report.str();
		}

		bool IsSelfAdjoint(ComparisonTolerance tolerance) const {
			if constexpr (MMLComplex<Scalar>)
				return IsHermitian(tolerance);
			else
				return IsSymmetric(tolerance);
		}

		static bool IsSelfAdjointResult(const MatrixAlg::MatrixAnalysis<Scalar>& analysis) {
			if constexpr (MMLComplex<Scalar>)
				return analysis.isHermitian;
			else
				return analysis.isSymmetric;
		}

		MatrixAlg::Definiteness ClassifyAcceptedSelfAdjointDefiniteness(
				ComparisonTolerance symmetryTolerance, Magnitude eigenvalueTolerance) const {
			const DefinitenessRequest request{symmetryTolerance, eigenvalueTolerance};
			return Cached(_acceptedDefinitenessCache, request, [&] {
				const int size = Rows();
				Matrix<Scalar> selfAdjoint(size, size);
				for (int row = 0; row < size; ++row)
					for (int col = 0; col < size; ++col)
						selfAdjoint(row, col) = Scalar{0.5} *
							(_matrix(row, col) + MatrixAlg::Detail::Conjugate(_matrix(col, row)));
				return MatrixAlg::ClassifyDefiniteness(selfAdjoint, eigenvalueTolerance);
			});
		}

		Matrix<Scalar> _matrix;
		mutable std::deque<CacheEntry<bool, MatrixAlg::LUDecomposition<Scalar>>> _luCache;
		mutable std::deque<CacheEntry<bool, MatrixAlg::QRDecomposition<Scalar>>> _qrCache;
		mutable std::deque<CacheEntry<Magnitude, MatrixAlg::CholeskyDecomposition<Scalar>>> _choleskyCache;
		mutable std::deque<CacheEntry<Threshold, MatrixAlg::SVDDecomposition<Scalar>>> _svdCache;
		mutable std::deque<CacheEntry<Threshold, MatrixAlg::FundamentalSubspaces<Scalar>>> _subspacesCache;
		mutable std::deque<CacheEntry<ComparisonTolerance, bool>> _symmetryCache;
		mutable std::deque<CacheEntry<Magnitude, MatrixAlg::Definiteness>> _definitenessCache;
		mutable std::deque<CacheEntry<DefinitenessRequest, MatrixAlg::Definiteness>> _acceptedDefinitenessCache;
		mutable std::deque<CacheEntry<EigenRequest, MatrixAlg::EigensystemResult<Scalar>>> _eigenCache;
		mutable std::deque<CacheEntry<EigenRequest, MatrixAlg::SelfAdjointEigensystemResult<Scalar>>> _selfAdjointEigenCache;
	};
}

#endif // MML_MATRIX_ANALYZER_H
