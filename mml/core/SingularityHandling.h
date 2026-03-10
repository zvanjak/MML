///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        SingularityHandling.h                                               ///
///  Description: Singularity detection and handling for curvilinear coordinates     ///
///               Provides configurable policies for operations at singular points   ///
///                                                                                   ///
///  SINGULARITY POLICY                                                               ///
///  ==================                                                               ///
///  Many curvilinear coordinate systems have unavoidable coordinate singularities:   ///
///    - Spherical: r=0 (origin), θ=0 or π (poles)                                    ///
///    - Cylindrical: r=0 (z-axis)                                                    ///
///                                                                                   ///
///  At these points, metric components become zero or infinite, and operations       ///
///  like gradient, divergence, Laplacian involve division by zero.                   ///
///                                                                                   ///
///  This header defines:                                                             ///
///    1. Detection functions to check if a point is near a singularity              ///
///    2. Policy enum to configure behavior at singularities                          ///
///    3. Safe arithmetic helpers that apply the configured policy                    ///
///                                                                                   ///
///  USAGE                                                                            ///
///  =====                                                                            ///
///    // Check before computing                                                      ///
///    if (Singularity::IsNearSphericalPole(pos)) { ... }                             ///
///                                                                                   ///
///    // Or use safe division with automatic policy application                      ///
///    Real result = Singularity::SafeDivide(a, r, policy);                           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_SINGULARITY_HANDLING_H
#define MML_SINGULARITY_HANDLING_H

#include "MMLBase.h"
#include "MMLExceptions.h"

#include <cmath>
#include <limits>
#include <string>

namespace MML
{
	///////////////////////////////////////////////////////////////////////////////////////////
	///                         SINGULARITY HANDLING POLICY                                ///
	///////////////////////////////////////////////////////////////////////////////////////////
	
	/// @brief Policy for handling mathematical singularities in coordinate operations
	///
	/// Defines behavior when an operation encounters a coordinate singularity
	/// (division by zero, undefined limits, etc.)
	enum class SingularityPolicy
	{
		/// Throw a DomainError exception with detailed message
		/// Use when: Caller must handle singular cases explicitly
		Throw,
		
		/// Return NaN (std::numeric_limits<Real>::quiet_NaN())
		/// Use when: Propagating undefined values for later detection
		ReturnNaN,
		
		/// Return signed infinity (+/-inf based on sign of numerator)
		/// Use when: Mathematically correct limit behavior is needed
		ReturnInf,
		
		/// Clamp denominator to a small epsilon, preserving sign
		/// Use when: Approximate but finite values are acceptable near singularities
		Clamp,
		
		/// Return zero (useful for some physical contexts)
		/// Use when: Zero is a sensible default (e.g., field strength at origin)
		ReturnZero
	};

	/// @brief Convert singularity policy to string for diagnostics
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

	///////////////////////////////////////////////////////////////////////////////////////////
	///                         SINGULARITY DETECTION AND SAFE OPERATIONS                  ///
	///////////////////////////////////////////////////////////////////////////////////////////
	
	namespace Singularity
	{
		/// @brief Default tolerance for singularity detection
		/// 
		/// Points within this distance of a singular coordinate are considered singular.
		/// Chosen to be above machine epsilon but small enough for most practical uses.
		static constexpr Real DEFAULT_SINGULARITY_TOL = 1e-12;

		/// @brief Default policy when none is specified
		static constexpr SingularityPolicy DEFAULT_POLICY = SingularityPolicy::Throw;

		///////////////////////////////////////////////////////////////////////////////////
		//                           SINGULARITY DETECTION                               //
		///////////////////////////////////////////////////////////////////////////////////

		/// @brief Check if a value is effectively zero (potential singularity)
		/// @param value Value to check
		/// @param tol Tolerance for near-zero detection
		/// @return true if |value| < tol
		inline bool IsNearZero(Real value, Real tol = DEFAULT_SINGULARITY_TOL) {
			return std::abs(value) < tol;
		}

		/// @brief Check if spherical radius r is at origin singularity
		/// @param r Radial coordinate
		/// @param tol Tolerance
		/// @return true if r < tol
		inline bool IsAtSphericalOrigin(Real r, Real tol = DEFAULT_SINGULARITY_TOL) {
			return r < tol;
		}

		/// @brief Check if spherical theta is at a pole (θ=0 or θ=π)
		/// 
		/// The poles are singular because sin(θ)=0, causing division by zero
		/// in φ-components of operators.
		/// 
		/// @param theta Polar angle in radians [0, π]
		/// @param tol Tolerance
		/// @return true if sin(theta) < tol
		inline bool IsAtSphericalPole(Real theta, Real tol = DEFAULT_SINGULARITY_TOL) {
			return std::abs(std::sin(theta)) < tol;
		}

		/// @brief Check if cylindrical radius is at z-axis singularity
		/// @param r Radial distance from z-axis
		/// @param tol Tolerance
		/// @return true if r < tol
		inline bool IsAtCylindricalAxis(Real r, Real tol = DEFAULT_SINGULARITY_TOL) {
			return r < tol;
		}

		/// @brief Combined check: is spherical position (r,θ,φ) at ANY singularity?
		/// @param r Radial coordinate
		/// @param theta Polar angle
		/// @param tol Tolerance
		/// @return true if at origin OR at pole
		inline bool IsAtSphericalSingularity(Real r, Real theta, 
		                                     Real tol = DEFAULT_SINGULARITY_TOL) {
			return IsAtSphericalOrigin(r, tol) || IsAtSphericalPole(theta, tol);
		}

		/// @brief Get descriptive string for spherical singularity type
		/// @return "origin" if r≈0, "pole" if θ≈0 or π, "none" otherwise
		inline const char* DescribeSphericalSingularity(Real r, Real theta,
		                                                Real tol = DEFAULT_SINGULARITY_TOL) {
			if (IsAtSphericalOrigin(r, tol)) return "origin (r=0)";
			if (IsAtSphericalPole(theta, tol)) return "pole (sin(θ)=0)";
			return "none";
		}

		///////////////////////////////////////////////////////////////////////////////////
		//                           SAFE ARITHMETIC OPERATIONS                          //
		///////////////////////////////////////////////////////////////////////////////////

		/// @brief Safe division with singularity policy
		///
		/// Performs a/b with configurable behavior when b is near zero.
		///
		/// @param numerator Numerator value
		/// @param denominator Denominator value (checked for near-zero)
		/// @param policy What to do if denominator is singular
		/// @param context Optional context string for error messages
		/// @param tol Tolerance for near-zero detection
		/// @return Result of division or policy-determined value
		/// @throws DomainError if policy is Throw and denominator is singular
		inline Real SafeDivide(Real numerator, Real denominator, 
		                       SingularityPolicy policy = DEFAULT_POLICY,
		                       const char* context = nullptr,
		                       Real tol = DEFAULT_SINGULARITY_TOL)
		{
			if (!IsNearZero(denominator, tol)) {
				return numerator / denominator;
			}

			// Denominator is near zero - apply policy
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
					// Clamp denominator to tol with same sign
					Real clampedDenom = (denominator >= 0) ? tol : -tol;
					return numerator / clampedDenom;
				}

				case SingularityPolicy::ReturnZero:
					return Real(0);

				default:
					return std::numeric_limits<Real>::quiet_NaN();
			}
		}

		/// @brief Safe 1/r operation common in spherical/cylindrical coordinates
		/// @param r Radial coordinate
		/// @param policy Singularity handling policy
		/// @param context Optional context for error messages
		/// @param tol Near-zero tolerance
		/// @return 1/r or policy-determined value
		inline Real SafeInverseR(Real r, 
		                         SingularityPolicy policy = DEFAULT_POLICY,
		                         const char* context = "1/r",
		                         Real tol = DEFAULT_SINGULARITY_TOL) {
			return SafeDivide(Real(1), r, policy, context, tol);
		}

		/// @brief Safe 1/r² operation
		/// @param r Radial coordinate  
		/// @param policy Singularity handling policy
		/// @param context Optional context for error messages
		/// @param tol Near-zero tolerance
		/// @return 1/r² or policy-determined value
		inline Real SafeInverseR2(Real r, 
		                          SingularityPolicy policy = DEFAULT_POLICY,
		                          const char* context = "1/r²",
		                          Real tol = DEFAULT_SINGULARITY_TOL) {
			return SafeDivide(Real(1), r * r, policy, context, tol);
		}

		/// @brief Safe 1/(r*sin(θ)) for spherical azimuthal component
		/// @param r Radial coordinate
		/// @param theta Polar angle
		/// @param policy Singularity handling policy
		/// @param context Optional context for error messages
		/// @param tol Near-zero tolerance
		/// @return 1/(r*sin(θ)) or policy-determined value
		inline Real SafeInverseRSinTheta(Real r, Real theta,
		                                 SingularityPolicy policy = DEFAULT_POLICY,
		                                 const char* context = "1/(r·sinθ)",
		                                 Real tol = DEFAULT_SINGULARITY_TOL) {
			return SafeDivide(Real(1), r * std::sin(theta), policy, context, tol);
		}

		/// @brief Safe 1/(r²*sin²(θ)) for spherical Laplacian azimuthal term
		/// @param r Radial coordinate
		/// @param theta Polar angle
		/// @param policy Singularity handling policy
		/// @param context Optional context for error messages
		/// @param tol Near-zero tolerance
		/// @return 1/(r²·sin²θ) or policy-determined value
		inline Real SafeInverseR2Sin2Theta(Real r, Real theta,
		                                   SingularityPolicy policy = DEFAULT_POLICY,
		                                   const char* context = "1/(r²·sin²θ)",
		                                   Real tol = DEFAULT_SINGULARITY_TOL) {
			Real sinTheta = std::sin(theta);
			return SafeDivide(Real(1), r * r * sinTheta * sinTheta, policy, context, tol);
		}

		/// @brief Safe cot(θ)/r² for spherical Laplacian polar term
		/// @param r Radial coordinate
		/// @param theta Polar angle
		/// @param policy Singularity handling policy
		/// @param context Optional context for error messages
		/// @param tol Near-zero tolerance
		/// @return cot(θ)/r² or policy-determined value
		inline Real SafeCotThetaOverR2(Real r, Real theta,
		                               SingularityPolicy policy = DEFAULT_POLICY,
		                               const char* context = "cotθ/r²",
		                               Real tol = DEFAULT_SINGULARITY_TOL) {
			Real sinTheta = std::sin(theta);
			Real cosTheta = std::cos(theta);
			return SafeDivide(cosTheta, r * r * sinTheta, policy, context, tol);
		}

	} // namespace Singularity

} // namespace MML

#endif // MML_SINGULARITY_HANDLING_H
