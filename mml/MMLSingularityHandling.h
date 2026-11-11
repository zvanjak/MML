///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MMLSingularityHandling.h                                            ///
///  Description: Shared singularity policy and safe arithmetic helpers               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_SINGULARITY_HANDLING_H
#define MML_SINGULARITY_HANDLING_H

#include <mml/MMLBase.h>
#include <mml/MMLExceptions.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

namespace MML
{
	/// Policy for handling mathematical singularities.
	enum class SingularityPolicy
	{
		Throw,
		ReturnNaN,
		ReturnInf,
		Clamp,
		ReturnZero
	};

	inline const char* SingularityPolicyToString(SingularityPolicy policy) {
		switch (policy) {
			case SingularityPolicy::Throw:      return "Throw";
			case SingularityPolicy::ReturnNaN:  return "ReturnNaN";
			case SingularityPolicy::ReturnInf:  return "ReturnInf";
			case SingularityPolicy::Clamp:      return "Clamp";
			case SingularityPolicy::ReturnZero: return "ReturnZero";
			default:                            return "Unknown";
		}
	}

	namespace Singularity
	{
		/// Default tolerance for singularity detection, adapted to the Real type.
		static constexpr Real DEFAULT_SINGULARITY_TOL = Precision::NumericalZeroThreshold;

		/// Default behavior when a singular operation is encountered.
		static constexpr SingularityPolicy DEFAULT_POLICY = SingularityPolicy::Throw;

		inline void RequireMatrixPolicy(SingularityPolicy policy) {
			switch (policy) {
			case SingularityPolicy::Throw:
			case SingularityPolicy::ReturnNaN:
			case SingularityPolicy::ReturnInf:
			case SingularityPolicy::ReturnZero:
			case SingularityPolicy::Clamp:
				return;
			default:
				throw ArgumentError("Unknown matrix singularity policy");
			}
		}

		template<int N, class TVector>
		void RequireFiniteCoordinates(const TVector& pos, const char* context) {
			for (int i = 0; i < N; ++i)
				if (!std::isfinite(pos[i]))
					throw DomainError(std::string(context) + ": nonfinite coordinates");
		}

		template<class TMatrix>
		void RequireFiniteMatrix(const TMatrix& matrix, const char* context) {
			for (int i = 0; i < matrix.rows(); ++i)
				for (int j = 0; j < matrix.cols(); ++j)
					if (!std::isfinite(matrix(i, j)))
						throw DomainError(std::string(context) + ": nonfinite matrix");
		}

		/// Explicit sentinel for a singular fixed-size real matrix, never an inverse.
		/// Clamp cannot select a unique inverse and must not silently regularize geometry.
		template<class TMatrix>
		TMatrix MatrixAtSingularity(SingularityPolicy policy, const char* context) {
			Real value;
			switch (policy) {
			case SingularityPolicy::Throw:
				throw SingularMatrixError(context);
			case SingularityPolicy::Clamp:
				throw DomainError(std::string(context) + ": Clamp cannot define a singular matrix inverse");
			case SingularityPolicy::ReturnNaN:
				value = std::numeric_limits<Real>::quiet_NaN();
				break;
			case SingularityPolicy::ReturnInf:
				value = std::numeric_limits<Real>::infinity();
				break;
			case SingularityPolicy::ReturnZero:
				value = 0;
				break;
			default:
				throw ArgumentError("Unknown matrix singularity policy");
			}
			TMatrix result;
			for (int i = 0; i < result.rows(); ++i)
				for (int j = 0; j < result.cols(); ++j)
					result(i, j) = value;
			return result;
		}



		/// Log absolute determinant via scaled partial-pivot elimination.
		/// Uses the LU determinant identity (Numerical Recipes, chapter 2).
		/// Avoids determinant underflow/overflow; pivot cutoff matches MatrixNM's norm-scaled rule.
		template<class TMatrix>
		Real LogAbsDeterminant(const TMatrix& matrix) {
			RequireFiniteMatrix(matrix, "Metric determinant");
			const int n = matrix.rows();
			if (n != matrix.cols() || n == 0)
				throw ArgumentError("LogAbsDeterminant requires a nonempty square matrix");
			Real scale = 0;
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					scale = std::max(scale, std::abs(matrix(i, j)));
			if (scale == 0)
				throw SingularMatrixError("Metric determinant: zero matrix");
			auto work = matrix;
			for (int i = 0; i < n; ++i)
				for (int j = 0; j < n; ++j)
					work(i, j) /= scale;
			Real result = n * std::log(scale);
			const Real cutoff = Precision::Eps * n;
			for (int col = 0; col < n; ++col) {
				int pivot = col;
				for (int row = col + 1; row < n; ++row)
					if (std::abs(work(row, col)) > std::abs(work(pivot, col)))
						pivot = row;
				if (!(std::abs(work(pivot, col)) > cutoff))
					throw SingularMatrixError("Metric determinant: rank-deficient or numerically singular matrix");
				if (pivot != col)
					for (int j = col; j < n; ++j)
						std::swap(work(col, j), work(pivot, j));
				result += std::log(std::abs(work(col, col)));
				for (int row = col + 1; row < n; ++row) {
					const Real factor = work(row, col) / work(col, col);
					for (int j = col + 1; j < n; ++j)
						work(row, j) -= factor * work(col, j);
				}
			}
			if (!std::isfinite(result))
				throw NumericalMethodError("Metric log determinant is nonfinite");
			return result;
		}

		inline bool IsNearZero(Real value, Real tol = DEFAULT_SINGULARITY_TOL) {
			return std::abs(value) < tol;
		}

		inline Real SafeDivide(Real numerator, Real denominator,
		                       SingularityPolicy policy = DEFAULT_POLICY,
		                       const char* context = nullptr,
		                       Real tol = DEFAULT_SINGULARITY_TOL)
		{
			if (!IsNearZero(denominator, tol)) {
				return numerator / denominator;
			}

			switch (policy)
			{
				case SingularityPolicy::Throw: {
					std::string msg = "Division by near-zero denominator";
					if (context) {
						msg += " in ";
						msg += context;
					}
					msg += " (denominator = " + std::to_string(denominator) + ")";
					throw DomainError(msg);
				}

				case SingularityPolicy::ReturnNaN:
					return std::numeric_limits<Real>::quiet_NaN();

				case SingularityPolicy::ReturnInf:
					if (numerator >= 0)
						return std::numeric_limits<Real>::infinity();
					else
						return -std::numeric_limits<Real>::infinity();

				case SingularityPolicy::Clamp: {
					Real clampedDenom = (denominator >= 0) ? tol : -tol;
					return numerator / clampedDenom;
				}

				case SingularityPolicy::ReturnZero:
					return Real(0);

				default:
					return std::numeric_limits<Real>::quiet_NaN();
			}
		}
	} // namespace Singularity
} // namespace MML

#endif // MML_SINGULARITY_HANDLING_H