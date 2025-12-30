///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FieldAnalyzers.h                                                    ///
///  Description: Vector field analysis tools (streamlines, potential, flux)          ///
///               Divergence, curl, gradient visualizations                           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FIELDS_ANALYZERS_H
#define MML_FIELDS_ANALYZERS_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/VectorN.h"
#include "base/Geometry3D.h"
#include "core/FieldOperations.h"
#include "core/Integration.h"
#include "core/Integration/PathIntegration.h"
#include "core/FunctionHelpers.h"

namespace MML
{
	/////////////////////////////////////////////////////////////////////////////////////
	///                           SCALAR FIELD ANALYZER                               ///
	/////////////////////////////////////////////////////////////////////////////////////
	
	// Analyzes properties of scalar fields in 3D space
	class ScalarFieldAnalyzer
	{
		const IScalarFunction<3>& _field;
		
	public:
		ScalarFieldAnalyzer(const IScalarFunction<3>& f) : _field(f) {}

		// Get the gradient at a point
		Vec3Cart GradientAt(const Vec3Cart& pos) const
		{
			return ScalarFieldOperations::GradientCart(_field, pos);
		}

		// Check if gradient is zero at a point (critical point)
		bool IsCriticalPoint(const Vec3Cart& pos, Real tol = 1e-6) const
		{
			Vec3Cart grad = GradientAt(pos);
			return grad.NormL2() < tol;
		}

		// Find approximate critical points by sampling a region
		std::vector<Vec3Cart> FindCriticalPoints(const Vec3Cart& minCorner, const Vec3Cart& maxCorner, 
		                                          int samplesPerDim = 10, Real tol = 1e-4) const
		{
			std::vector<Vec3Cart> criticalPoints;
			
			Real dx = (maxCorner.X() - minCorner.X()) / samplesPerDim;
			Real dy = (maxCorner.Y() - minCorner.Y()) / samplesPerDim;
			Real dz = (maxCorner.Z() - minCorner.Z()) / samplesPerDim;
			
			for (int i = 0; i <= samplesPerDim; ++i) {
				for (int j = 0; j <= samplesPerDim; ++j) {
					for (int k = 0; k <= samplesPerDim; ++k) {
						Vec3Cart pos(minCorner.X() + i * dx,
						             minCorner.Y() + j * dy,
						             minCorner.Z() + k * dz);
						if (IsCriticalPoint(pos, tol)) {
							criticalPoints.push_back(pos);
						}
					}
				}
			}
			return criticalPoints;
		}

		// Compute Laplacian at a point (∇²f = ∂²f/∂x² + ∂²f/∂y² + ∂²f/∂z²)
		Real LaplacianAt(const Vec3Cart& pos) const
		{
			return ScalarFieldOperations::LaplacianCart(_field, pos);
		}

		// Check if field satisfies Laplace equation (harmonic) in a region
		bool IsHarmonic(const Vec3Cart& minCorner, const Vec3Cart& maxCorner, 
		                int samplesPerDim = 5, Real tol = 1e-4) const
		{
			Real dx = (maxCorner.X() - minCorner.X()) / samplesPerDim;
			Real dy = (maxCorner.Y() - minCorner.Y()) / samplesPerDim;
			Real dz = (maxCorner.Z() - minCorner.Z()) / samplesPerDim;
			
			for (int i = 0; i <= samplesPerDim; ++i) {
				for (int j = 0; j <= samplesPerDim; ++j) {
					for (int k = 0; k <= samplesPerDim; ++k) {
						Vec3Cart pos(minCorner.X() + i * dx,
						             minCorner.Y() + j * dy,
						             minCorner.Z() + k * dz);
						if (std::abs(LaplacianAt(pos)) > tol) {
							return false;
						}
					}
				}
			}
			return true;
		}
	};

	/////////////////////////////////////////////////////////////////////////////////////
	///                           VECTOR FIELD ANALYZER                               ///
	/////////////////////////////////////////////////////////////////////////////////////
	
	// Analyzes properties of vector fields in 3D space
	class VectorFieldAnalyzer
	{
		const IVectorFunction<3>& _field;
		
	public:
		VectorFieldAnalyzer(const IVectorFunction<3>& f) : _field(f) {}

		// Get divergence at a point
		Real DivergenceAt(const Vec3Cart& pos) const
		{
			return VectorFieldOperations::DivCart(_field, pos);
		}

		// Get curl at a point
		Vec3Cart CurlAt(const Vec3Cart& pos) const
		{
			return VectorFieldOperations::CurlCart(_field, pos);
		}

		// Check if field is solenoidal (divergence-free) at a point
		bool IsSolenoidalAt(const Vec3Cart& pos, Real tol = 1e-6) const
		{
			return std::abs(DivergenceAt(pos)) < tol;
		}

		// Check if field is irrotational (curl-free) at a point
		bool IsIrrotationalAt(const Vec3Cart& pos, Real tol = 1e-6) const
		{
			return CurlAt(pos).NormL2() < tol;
		}

		// Check if field is solenoidal throughout a region (sampling-based)
		// Solenoidal fields have ∇·F = 0 (no sources or sinks)
		bool IsSolenoidal(const Vec3Cart& minCorner, const Vec3Cart& maxCorner, 
		                  int samplesPerDim = 5, Real tol = 1e-4) const
		{
			Real dx = (maxCorner.X() - minCorner.X()) / samplesPerDim;
			Real dy = (maxCorner.Y() - minCorner.Y()) / samplesPerDim;
			Real dz = (maxCorner.Z() - minCorner.Z()) / samplesPerDim;
			
			for (int i = 0; i <= samplesPerDim; ++i) {
				for (int j = 0; j <= samplesPerDim; ++j) {
					for (int k = 0; k <= samplesPerDim; ++k) {
						Vec3Cart pos(minCorner.X() + i * dx,
						             minCorner.Y() + j * dy,
						             minCorner.Z() + k * dz);
						if (!IsSolenoidalAt(pos, tol)) {
							return false;
						}
					}
				}
			}
			return true;
		}

		// Check if field is irrotational throughout a region (sampling-based)
		// Irrotational fields have ∇×F = 0
		bool IsIrrotational(const Vec3Cart& minCorner, const Vec3Cart& maxCorner, 
		                    int samplesPerDim = 5, Real tol = 1e-4) const
		{
			Real dx = (maxCorner.X() - minCorner.X()) / samplesPerDim;
			Real dy = (maxCorner.Y() - minCorner.Y()) / samplesPerDim;
			Real dz = (maxCorner.Z() - minCorner.Z()) / samplesPerDim;
			
			for (int i = 0; i <= samplesPerDim; ++i) {
				for (int j = 0; j <= samplesPerDim; ++j) {
					for (int k = 0; k <= samplesPerDim; ++k) {
						Vec3Cart pos(minCorner.X() + i * dx,
						             minCorner.Y() + j * dy,
						             minCorner.Z() + k * dz);
						if (!IsIrrotationalAt(pos, tol)) {
							return false;
						}
					}
				}
			}
			return true;
		}

		// Check if field is conservative (path-independent)
		// For simply-connected regions, conservative ⟺ irrotational
		bool IsConservative(const Vec3Cart& minCorner, const Vec3Cart& maxCorner, 
		                    int samplesPerDim = 5, Real tol = 1e-4) const
		{
			// In simply-connected domains, conservative = irrotational
			return IsIrrotational(minCorner, maxCorner, samplesPerDim, tol);
		}

		// Get maximum divergence magnitude in a region
		Real MaxDivergence(const Vec3Cart& minCorner, const Vec3Cart& maxCorner, 
		                   int samplesPerDim = 5) const
		{
			Real maxDiv = 0.0;
			Real dx = (maxCorner.X() - minCorner.X()) / samplesPerDim;
			Real dy = (maxCorner.Y() - minCorner.Y()) / samplesPerDim;
			Real dz = (maxCorner.Z() - minCorner.Z()) / samplesPerDim;
			
			for (int i = 0; i <= samplesPerDim; ++i) {
				for (int j = 0; j <= samplesPerDim; ++j) {
					for (int k = 0; k <= samplesPerDim; ++k) {
						Vec3Cart pos(minCorner.X() + i * dx,
						             minCorner.Y() + j * dy,
						             minCorner.Z() + k * dz);
						maxDiv = std::max(maxDiv, std::abs(DivergenceAt(pos)));
					}
				}
			}
			return maxDiv;
		}

		// Get maximum curl magnitude in a region
		Real MaxCurlMagnitude(const Vec3Cart& minCorner, const Vec3Cart& maxCorner, 
		                      int samplesPerDim = 5) const
		{
			Real maxCurl = 0.0;
			Real dx = (maxCorner.X() - minCorner.X()) / samplesPerDim;
			Real dy = (maxCorner.Y() - minCorner.Y()) / samplesPerDim;
			Real dz = (maxCorner.Z() - minCorner.Z()) / samplesPerDim;
			
			for (int i = 0; i <= samplesPerDim; ++i) {
				for (int j = 0; j <= samplesPerDim; ++j) {
					for (int k = 0; k <= samplesPerDim; ++k) {
						Vec3Cart pos(minCorner.X() + i * dx,
						             minCorner.Y() + j * dy,
						             minCorner.Z() + k * dz);
						maxCurl = std::max(maxCurl, CurlAt(pos).NormL2());
					}
				}
			}
			return maxCurl;
		}

		// Compute circulation around a closed circular path
		Real CirculationAroundCircle(const Vec3Cart& center, const Vec3Cart& normal, 
		                             Real radius, int numPoints = 100) const
		{
			// Create orthonormal basis for the plane
			Vec3Cart n = normal.GetAsUnitVector();
			Vec3Cart u, v;
			
			// Find a vector not parallel to n
			if (std::abs(n.X()) < 0.9) {
				u = VectorProduct(n, Vec3Cart(1, 0, 0)).GetAsUnitVector();
			} else {
				u = VectorProduct(n, Vec3Cart(0, 1, 0)).GetAsUnitVector();
			}
			v = VectorProduct(n, u);
			
			// Integrate F·dr around the circle
			Real circulation = 0.0;
			Real dtheta = 2.0 * Constants::PI / numPoints;
			
			for (int i = 0; i < numPoints; ++i) {
				Real theta = i * dtheta;
				Real theta_next = (i + 1) * dtheta;
				
				// Current and next points on circle
				Vec3Cart p1 = center + u * (radius * std::cos(theta)) + v * (radius * std::sin(theta));
				Vec3Cart p2 = center + u * (radius * std::cos(theta_next)) + v * (radius * std::sin(theta_next));
				
				// Tangent direction (dr)
				Vec3Cart dr = p2 - p1;
				
				// Midpoint for field evaluation
				Vec3Cart midpoint = (p1 + p2) * 0.5;
				VectorN<Real, 3> F = _field(midpoint);
				
				// Add contribution F·dr
				circulation += F[0] * dr.X() + F[1] * dr.Y() + F[2] * dr.Z();
			}
			
			return circulation;
		}
	};

} // end namespace MML

#endif // MML_FIELDS_ANALYZERS_H
