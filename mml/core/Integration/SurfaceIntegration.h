///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        SurfaceIntegration.h                                                ///
///  Description: Surface integrals over parametric surfaces                          ///
///               Flux integrals and scalar surface integration                       ///
///                                                                                   ///
///  Features:    - Flux integrals ∬_S F·n dS (vector field through surface)         ///
///               - Parametric surface integration with numerical derivatives         ///
///               - Triangle mesh surface integration with adaptive refinement        ///
///               - Rectangular mesh surface integration with adaptive refinement     ///
///                                                                                   ///
///  Usage:                                                                           ///
///    // Flux through parametric surface                                             ///
///    Real flux = SurfaceIntegration::SurfaceIntegral(field, surface,                ///
///                                                    u1, u2, v1, v2);               ///
///                                                                                   ///
///    // Flux through closed triangular mesh (e.g., tetrahedron)                     ///
///    Real flux = SurfaceIntegration::SurfaceIntegral(field, solidBody);             ///
///                                                                                   ///
///  Mathematical Background:                                                         ///
///    Surface integral:  ∬_S F·dS = ∬_D F(r(u,v))·(∂r/∂u × ∂r/∂v) du dv            ///
///    Flux:              Φ = ∬_S F·n dS  (measures flow through surface)            ///
///                                                                                   ///
///  Gauss's Divergence Theorem:                                                      ///
///    ∬_S F·n dS = ∭_V ∇·F dV  (flux through closed surface = divergence in volume) ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_SURFACE_INTEGRATION_H
#define MML_SURFACE_INTEGRATION_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/Vector/Vector.h"
#include "base/Vector/VectorN.h"
#include "mml/base/Geometry/Geometry3D.h"

#include "core/Derivation.h"
#include "core/Integration.h"
#include "core/FieldOperations.h"

namespace MML
{
	///////////////////////////////////////////////////////////////////////////
	///                       SurfaceIntegration                            ///
	///////////////////////////////////////////////////////////////////////////
	/// @brief Static class for computing surface integrals (flux integrals)
	/// 
	/// Provides methods for computing flux of vector fields through surfaces:
	/// - Parametric surfaces with rectangular parameter domains
	/// - Triangular mesh surfaces (BodyWithTriangleSurfaces)
	/// - Rectangular mesh surfaces (BodyWithRectSurfaces)
	///
	/// @par Supported surface types:
	/// - IParametricSurfaceRect: Surfaces defined by r(u,v) over [u1,u2]×[v1,v2]
	/// - Triangle3D: Individual triangular patches
	/// - RectSurface3D: Individual rectangular (quadrilateral) patches
	/// - BodyWithTriangleSurfaces: Closed surfaces made of triangles
	/// - BodyWithRectSurfaces: Closed surfaces made of rectangles
	///
	/// @par Integration methods:
	/// - Midpoint rule for parametric surfaces
	/// - Adaptive recursive refinement for mesh surfaces
	///
	/// @see PathIntegration for line integrals
	/// @see VolumeIntegration for volume integrals
	///////////////////////////////////////////////////////////////////////////
	class SurfaceIntegration
	{
	public:
		///////////////////////////////////////////////////////////////////////
		///              Parametric Surface Integration                     ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Flux integral over a parametric surface
		/// 
		/// Computes ∬_S F·dS = ∬_D F(r(u,v))·(∂r/∂u × ∂r/∂v) du dv
		/// using midpoint rule with numerical derivatives.
		///
		/// @param vectorField Vector field F: ℝ³ → ℝ³ to integrate
		/// @param surface Parametric surface r(u,v): D ⊂ ℝ² → ℝ³
		/// @param x1 Starting u-parameter value
		/// @param x2 Ending u-parameter value
		/// @param y1 Starting v-parameter value
		/// @param y2 Ending v-parameter value
		/// @param numU Number of subdivisions in u-direction (default: 20)
		/// @param numV Number of subdivisions in v-direction (default: 20)
		/// @return Flux of vector field through the surface
		///
		/// @par Algorithm:
		/// 1. Divide parameter domain into numU × numV patches
		/// 2. For each patch, compute surface normal via cross product:
		///    n = (∂r/∂u × ∂r/∂v) / |∂r/∂u × ∂r/∂v|
		/// 3. Evaluate F·n at patch center, multiply by area element
		/// 4. Sum all contributions
		///
		/// @par Example:
		/// @code
		/// // Flux of radial field through sphere
		/// VectorFunction<3> radialField([](const VectorN<3>& p) {
		///     return p;  // F = r (pointing outward)
		/// });
		/// // Sphere parametrization: r(θ,φ) = (sinθ cosφ, sinθ sinφ, cosθ)
		/// ParametricSurfaceRect<3> sphere([](Real theta, Real phi) {
		///     return VectorN<3>({sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)});
		/// });
		/// Real flux = SurfaceIntegration::SurfaceIntegral(radialField, sphere,
		///                                                  0, PI, 0, 2*PI, 40, 40);
		/// // By Gauss's theorem: flux = ∭_V ∇·F dV = ∭_V 3 dV = 4π (for unit sphere)
		/// @endcode
		///
		/// @note Uses central differences for derivative computation
		/// @warning Accuracy depends on numU, numV; increase for complex surfaces
		static IntegrationResult SurfaceIntegral(const IVectorFunction<3>& vectorField, const IParametricSurfaceRect<3>& surface, 
		                            const Real x1, const Real x2, const Real y1, const Real y2, 
		                            int numU = 20, int numV = 20)
		{
			if (x2 <= x1 || y2 <= y1)
				return IntegrationResult(0.0, 0.0, 0, true);  // Invalid parameter range

			Real du = (x2 - x1) / numU;
			Real dv = (y2 - y1) / numV;
			Real total = 0.0;

			// Integrate using midpoint rule with adaptive normal computation
			for (int i = 0; i < numU; i++)
			{
				for (int j = 0; j < numV; j++)
				{
					// Midpoint of current patch
					Real u = x1 + (i + 0.5) * du;
					Real v = y1 + (j + 0.5) * dv;

					// Compute surface point at midpoint
					VectorN<Real, 3> point = surface(u, v);

					// Compute partial derivatives using central differences
					Real h = PrecisionValues<Real>::DerivativeStepSize;

					VectorN<Real, 3> point_u_plus = surface(u + h, v);
					VectorN<Real, 3> point_u_minus = surface(u - h, v);
					VectorN<Real, 3> dS_du = (point_u_plus - point_u_minus) / (2.0 * h);

					VectorN<Real, 3> point_v_plus = surface(u, v + h);
					VectorN<Real, 3> point_v_minus = surface(u, v - h);
					VectorN<Real, 3> dS_dv = (point_v_plus - point_v_minus) / (2.0 * h);

					// Compute cross product: ∂r/∂u × ∂r/∂v (normal vector, not normalized)
					Vector3Cartesian partial_u(dS_du[0], dS_du[1], dS_du[2]);
					Vector3Cartesian partial_v(dS_dv[0], dS_dv[1], dS_dv[2]);
					Vector3Cartesian crossProd = VectorProduct(partial_u, partial_v);
					
					// Area element magnitude
					Real dS_magnitude = crossProd.NormL2();
					if (dS_magnitude < PrecisionValues<Real>::SurfaceNormalThreshold)
						continue;  // Skip degenerate patch

					// Unit normal vector
					Vector3Cartesian normal = crossProd / dS_magnitude;

					// Evaluate vector field at surface point
					Vector3Cartesian fieldValue = vectorField(Vector3Cartesian(point[0], point[1], point[2]));

					// Flux contribution: (F·n) dS where dS = |∂r/∂u × ∂r/∂v| du dv
					Real dotProduct = ScalarProduct(normal, fieldValue);
					Real patchArea = dS_magnitude * du * dv;

					total += dotProduct * patchArea;
				}
			}

			return IntegrationResult(total, 0.0, numU * numV, true);
		}

		///////////////////////////////////////////////////////////////////////
		///              Triangle Mesh Surface Integration                  ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Flux integral over a body with triangular surfaces
		/// 
		/// Computes ∬_S F·n dS for a closed surface composed of triangles.
		/// Sums the flux contribution from each triangular face.
		///
		/// @param vectorField Vector field F: ℝ³ → ℝ³
		/// @param solid Body defined by triangular surfaces
		/// @param eps Convergence tolerance for adaptive refinement (default: 0.001)
		/// @return Total flux through all surfaces
		///
		/// @par Physical interpretation:
		/// For a closed surface enclosing volume V, by Gauss's theorem:
		/// ∬_S F·n dS = ∭_V ∇·F dV
		///
		/// @par Example:
		/// @code
		/// // Verify Gauss's theorem for a tetrahedron
		/// Tetrahedron tetra(p1, p2, p3, p4);
		/// VectorFunction<3> F([](const VectorN<3>& p) {
		///     return VectorN<3>({p[0], p[1], p[2]});  // F = r, ∇·F = 3
		/// });
		/// Real flux = SurfaceIntegration::SurfaceIntegral(F, tetra);
		/// // Should equal 3 * volume of tetrahedron
		/// @endcode
		static IntegrationResult SurfaceIntegral(const IVectorFunction<3>& vectorField, const BodyWithTriangleSurfaces& solid, Real eps = 0.001)
		{
			Real total = 0.0;
			for (int i = 0; i < solid.GetSurfaceCount(); i++)
				total += SurfaceIntegral(vectorField, solid.GetSurface(i), eps).value;
			return IntegrationResult(total, 0.0, solid.GetSurfaceCount(), true);
		}

		/// @brief Flux integral over a single triangle with adaptive refinement
		/// 
		/// Uses recursive subdivision to achieve desired accuracy.
		/// Each triangle is subdivided into 4 smaller triangles by connecting
		/// edge midpoints until convergence or max level is reached.
		///
		/// @param vectorField Vector field F: ℝ³ → ℝ³
		/// @param triangle The triangular surface patch
		/// @param eps Convergence tolerance
		/// @param maxLevel Maximum recursion depth (default: 7)
		/// @return Flux through the triangle
		///
		/// @par Subdivision pattern:
		/// @code
		///       A                    A
		///      / \                  /|\
		///     /   \     --->      m1-+-m3
		///    /     \              |\ | /|
		///   B-------C             B--m2--C
		/// @endcode
		/// where m1, m2, m3 are edge midpoints.
		static IntegrationResult SurfaceIntegral(const IVectorFunction<3>& vectorField, const Triangle3D& triangle, Real eps, int maxLevel = 7)
		{
			Real result = CalcSurfaceContrib(vectorField, triangle);
			Real refined = SurfaceIntegralImproveRecursively(vectorField, triangle, result, eps, maxLevel);
			return IntegrationResult(refined, 0.0, 1, true);
		}

	private:
		/// @brief Recursive refinement for triangle surface integration
		/// @details Subdivides triangle into 4 smaller triangles and checks convergence
		static Real SurfaceIntegralImproveRecursively(const IVectorFunction<3>& vectorField, const Triangle3D& triangle, Real prev_value, Real eps, int level)
		{
			// Subdivide triangle into 4 by connecting edge midpoints
			Pnt3Cart m12 = (triangle.Pnt1() + triangle.Pnt2()) / 2.0;
			Pnt3Cart m23 = (triangle.Pnt2() + triangle.Pnt3()) / 2.0;
			Pnt3Cart m31 = (triangle.Pnt3() + triangle.Pnt1()) / 2.0;

			Triangle3D t1(triangle.Pnt1(), m12, m31);  // Corner triangle at vertex 1
			Triangle3D t2(m12, triangle.Pnt2(), m23);  // Corner triangle at vertex 2
			Triangle3D t3(m31, m23, triangle.Pnt3());  // Corner triangle at vertex 3
			Triangle3D t4(m12, m23, m31);              // Central triangle

			Real c1 = CalcSurfaceContrib(vectorField, t1);
			Real c2 = CalcSurfaceContrib(vectorField, t2);
			Real c3 = CalcSurfaceContrib(vectorField, t3);
			Real c4 = CalcSurfaceContrib(vectorField, t4);
			Real result = c1 + c2 + c3 + c4;

			// Check convergence
			if (fabs(result - prev_value) < eps || level == 0)
			{
				return result;
			}
			else
			{
				// Recurse on each sub-triangle with its own contribution as prev_value
				Real new_result = 0.0;
				new_result += SurfaceIntegralImproveRecursively(vectorField, t1, c1, eps, level - 1);
				new_result += SurfaceIntegralImproveRecursively(vectorField, t2, c2, eps, level - 1);
				new_result += SurfaceIntegralImproveRecursively(vectorField, t3, c3, eps, level - 1);
				new_result += SurfaceIntegralImproveRecursively(vectorField, t4, c4, eps, level - 1);
				return new_result;
			}
		}

		/// @brief Calculate flux contribution of a single triangle using centroid rule
		/// @details Evaluates F·n at triangle center, multiplies by triangle area
		static Real CalcSurfaceContrib(const IVectorFunction<3>& vectorField, const Triangle3D& triangle)
		{
			// Compute outward normal using right-hand rule: edge1 × edge2
			Vec3Cart edge1(triangle.Pnt1(), triangle.Pnt2());
			Vec3Cart edge2(triangle.Pnt1(), triangle.Pnt3());
			Vec3Cart normal = VectorProduct(edge1, edge2).GetAsUnitVector();

			// Triangle area
			Real area = triangle.Area();

			// Centroid (center of mass)
			Pnt3Cart center = (triangle.Pnt1() + triangle.Pnt2() + triangle.Pnt3()) / 3.0;

			// Evaluate field and compute dot product
			Vector3Cartesian value = vectorField(Vector3Cartesian(center.X(), center.Y(), center.Z()));
			Real dotProduct = ScalarProduct(normal, value);

			return dotProduct * area;
		}

	public:
		///////////////////////////////////////////////////////////////////////
		///            Rectangular Mesh Surface Integration                 ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Flux integral over a body with rectangular surfaces
		/// 
		/// Computes ∬_S F·n dS for a closed surface composed of quadrilaterals.
		///
		/// @param vectorField Vector field F: ℝ³ → ℝ³
		/// @param solid Body defined by rectangular surfaces (e.g., Box3D)
		/// @param eps Convergence tolerance for adaptive refinement (default: 0.001)
		/// @return Total flux through all surfaces
		///
		/// @par Example:
		/// @code
		/// // Flux through a box
		/// Box3D box(Point3Cartesian(0,0,0), Point3Cartesian(1,1,1));
		/// VectorFunction<3> F([](const VectorN<3>& p) {
		///     return VectorN<3>({1, 0, 0});  // Uniform field in x-direction
		/// });
		/// Real flux = SurfaceIntegration::SurfaceIntegral(F, box);
		/// // Net flux = 0 (same flux in through left face, out through right)
		/// @endcode
		static IntegrationResult SurfaceIntegral(const IVectorFunction<3>& vectorField, const BodyWithRectSurfaces& solid, Real eps = 0.001)
		{
			Real total = 0.0;
			for (int i = 0; i < solid.GetSurfaceCount(); i++)
				total += SurfaceIntegral(vectorField, solid.GetSurface(i), eps).value;
			return IntegrationResult(total, 0.0, solid.GetSurfaceCount(), true);
		}

		/// @brief Flux integral over a single rectangular surface with adaptive refinement
		/// 
		/// Uses recursive subdivision to achieve desired accuracy.
		/// Each rectangle is subdivided into 4 smaller rectangles.
		///
		/// @param vectorField Vector field F: ℝ³ → ℝ³
		/// @param surface The rectangular surface patch
		/// @param eps Convergence tolerance
		/// @param maxLevel Maximum recursion depth (default: 7)
		/// @return IntegrationResult with flux value (implicitly converts to Real)
		static IntegrationResult SurfaceIntegral(const IVectorFunction<3>& vectorField, const RectSurface3D& surface, Real eps, int maxLevel = 7)
		{
			Real result = CalcSurfaceContrib(vectorField, surface);
			Real refined = SurfaceIntegralImproveRecursively(vectorField, surface, result, eps, maxLevel);
			return IntegrationResult(refined, 0.0, 1, true);
		}

	private:
		/// @brief Recursive refinement for rectangular surface integration
		/// @details Subdivides rectangle into 4 by connecting edge midpoints to center
		static Real SurfaceIntegralImproveRecursively(const IVectorFunction<3>& vectorField, const RectSurface3D& surface, Real prev_value, Real eps, int level)
		{
			// Compute midpoints of edges and center
			Point3Cartesian pnt_mid12 = (surface._pnt1 + surface._pnt2) / 2.0;
			Point3Cartesian pnt_mid23 = (surface._pnt2 + surface._pnt3) / 2.0;
			Point3Cartesian pnt_mid34 = (surface._pnt3 + surface._pnt4) / 2.0;
			Point3Cartesian pnt_mid41 = (surface._pnt4 + surface._pnt1) / 2.0;
			Point3Cartesian center = surface.getCenter();

			// Create 4 sub-rectangles
			RectSurface3D s1(surface._pnt1, pnt_mid12, center, pnt_mid41);
			RectSurface3D s2(pnt_mid12, surface._pnt2, pnt_mid23, center);
			RectSurface3D s3(center, pnt_mid23, surface._pnt3, pnt_mid34);
			RectSurface3D s4(pnt_mid41, center, pnt_mid34, surface._pnt4);

			// Calculate contributions for each subsurface
			Real contrib1 = CalcSurfaceContrib(vectorField, s1);
			Real contrib2 = CalcSurfaceContrib(vectorField, s2);
			Real contrib3 = CalcSurfaceContrib(vectorField, s3);
			Real contrib4 = CalcSurfaceContrib(vectorField, s4);
			Real result = contrib1 + contrib2 + contrib3 + contrib4;

			// Check convergence
			if (fabs(result - prev_value) < eps || level == 0)
			{
				return result;
			}
			else
			{
				// Recurse on each sub-rectangle
				Real new_result = 0.0;
				new_result += SurfaceIntegralImproveRecursively(vectorField, s1, contrib1, eps, level - 1);
				new_result += SurfaceIntegralImproveRecursively(vectorField, s2, contrib2, eps, level - 1);
				new_result += SurfaceIntegralImproveRecursively(vectorField, s3, contrib3, eps, level - 1);
				new_result += SurfaceIntegralImproveRecursively(vectorField, s4, contrib4, eps, level - 1);
				return new_result;
			}
		}

		/// @brief Calculate flux contribution of a single rectangle using center rule
		/// @details Evaluates F·n at rectangle center, multiplies by rectangle area
		static Real CalcSurfaceContrib(const IVectorFunction<3>& vectorField, const RectSurface3D& surface)
		{
			Vector3Cartesian normal = surface.getNormal();
			Real area = surface.getArea();
			Point3Cartesian center = surface.getCenter();

			Vector3Cartesian value = vectorField(Vector3Cartesian(center.X(), center.Y(), center.Z()));
			Real dotProduct = ScalarProduct(normal, value);

			return dotProduct * area;
		}
	};

} // end namespace MML

#endif // MML_SURFACE_INTEGRATION_H
