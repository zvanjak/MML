///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        CoordinateSingularities.h                                           ///
///  Description: Detection and safe operations for coordinate singularities          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COORDINATE_SINGULARITIES_H
#define MML_COORDINATE_SINGULARITIES_H

#include <mml/MMLSingularityHandling.h>

#include <cmath>

namespace MML::Singularity
{
	inline bool IsAtSphericalOrigin(Real r, Real tol = DEFAULT_SINGULARITY_TOL) {
		return r < tol;
	}

	inline bool IsAtSphericalPole(Real theta, Real tol = DEFAULT_SINGULARITY_TOL) {
		return std::abs(std::sin(theta)) < tol;
	}

	inline bool IsAtCylindricalAxis(Real r, Real tol = DEFAULT_SINGULARITY_TOL) {
		return r < tol;
	}

	inline bool IsAtSphericalSingularity(Real r, Real theta,
	                                     Real tol = DEFAULT_SINGULARITY_TOL) {
		return IsAtSphericalOrigin(r, tol) || IsAtSphericalPole(theta, tol);
	}

	inline const char* DescribeSphericalSingularity(Real r, Real theta,
	                                                Real tol = DEFAULT_SINGULARITY_TOL) {
		if (IsAtSphericalOrigin(r, tol)) return "origin (r=0)";
		if (IsAtSphericalPole(theta, tol)) return "pole (sin(θ)=0)";
		return "none";
	}

	inline Real SafeInverseR(Real r,
	                         SingularityPolicy policy = DEFAULT_POLICY,
	                         const char* context = "1/r",
	                         Real tol = DEFAULT_SINGULARITY_TOL) {
		return SafeDivide(Real(1), r, policy, context, tol);
	}

	inline Real SafeInverseR2(Real r,
	                          SingularityPolicy policy = DEFAULT_POLICY,
	                          const char* context = "1/r²",
	                          Real tol = DEFAULT_SINGULARITY_TOL) {
		return SafeDivide(Real(1), r * r, policy, context, tol);
	}

	inline Real SafeInverseRSinTheta(Real r, Real theta,
	                                 SingularityPolicy policy = DEFAULT_POLICY,
	                                 const char* context = "1/(r·sinθ)",
	                                 Real tol = DEFAULT_SINGULARITY_TOL) {
		return SafeDivide(Real(1), r * std::sin(theta), policy, context, tol);
	}

	inline Real SafeInverseR2Sin2Theta(Real r, Real theta,
	                                   SingularityPolicy policy = DEFAULT_POLICY,
	                                   const char* context = "1/(r²·sin²θ)",
	                                   Real tol = DEFAULT_SINGULARITY_TOL) {
		Real sinTheta = std::sin(theta);
		return SafeDivide(Real(1), r * r * sinTheta * sinTheta, policy, context, tol);
	}

	inline Real SafeCotThetaOverR2(Real r, Real theta,
	                               SingularityPolicy policy = DEFAULT_POLICY,
	                               const char* context = "cotθ/r²",
	                               Real tol = DEFAULT_SINGULARITY_TOL) {
		Real sinTheta = std::sin(theta);
		Real cosTheta = std::cos(theta);
		return SafeDivide(cosTheta, r * r * sinTheta, policy, context, tol);
	}
} // namespace MML::Singularity

#endif // MML_COORDINATE_SINGULARITIES_H