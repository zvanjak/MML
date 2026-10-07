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

#include <cmath>
#include <limits>
#include <string>

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