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
///    Real flux = SurfaceIntegration::SurfaceIntegralDetailed(field, surface,                ///
///                                                    u1, u2, v1, v2);               ///
///                                                                                   ///
///    // Flux through closed triangular mesh (e.g., tetrahedron)                     ///
///    Real flux = SurfaceIntegration::SurfaceIntegralDetailed(field, solidBody);             ///
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

#include <mml/MMLBase.h>

#include <mml/interfaces/IFunction.h>

#include <mml/base/Vector/Vector.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Geometry/Geometry3D.h>
#include <mml/base/Geometry/Geometry3DBodies.h>

#include <mml/core/Derivation.h>
#include <mml/core/Integration.h>
#include <mml/core/Fields/FieldOperations.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

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
		/// Real flux = SurfaceIntegration::SurfaceIntegralDetailed(radialField, sphere,
		///                                                  0, PI, 0, 2*PI, 40, 40);
		/// // By Gauss's theorem: flux = ∭_V ∇·F dV = ∭_V 3 dV = 4π (for unit sphere)
		/// @endcode
		///
		/// @note Uses central differences for derivative computation
		/// @warning Accuracy depends on numU, numV; increase for complex surfaces
		static IntegrationResult SurfaceIntegralDetailed(const IVectorFunction<3>& vectorField, const IParametricSurfaceRect<3>& surface, 
		                            const Real x1, const Real x2, const Real y1, const Real y2, 
		                            int numU = 20, int numV = 20)
		{
			ValidateBounds(x1, x2, "SurfaceIntegral u");
			ValidateBounds(y1, y2, "SurfaceIntegral v");
			if (x2 <= x1 || y2 <= y1)
				throw ArgumentError("SurfaceIntegral: parameter range must satisfy x1 < x2 and y1 < y2");
			if (numU <= 0 || numV <= 0)
				throw ArgumentError("SurfaceIntegral: numU and numV must be positive");
			if (numV > std::numeric_limits<int>::max() / numU)
				throw ResourceLimitError("SurfaceIntegral: midpoint grid exceeds evaluation counter");
			const Real minU = surface.getMinU(), maxU = surface.getMaxU();
			const Real minV = surface.getMinW(), maxV = surface.getMaxW();
			if (!(minU < maxU && minV < maxV) || x1 < minU || x2 > maxU || y1 < minV || y2 > maxV)
				throw ArgumentError("SurfaceIntegral: requested range must lie inside the declared surface domain");

			AlgorithmTimer timer;
			const int patches = numU * numV;
			const Real halfU = IntegrationDetail::HalfWidth(x1, x2);
			const Real halfV = IntegrationDetail::HalfWidth(y1, y2);
			IntegrationDetail::Sum total;
			auto sample = [&](Real u, Real v) {
				auto point = surface(u, v);
				RequireFiniteVector(point, "Surface point");
				return point;
			};
			for (int i = 0; i < numU; i++)
			{
				for (int j = 0; j < numV; j++)
				{
					const Real u = std::lerp(x1, x2, (Real(i) + Real(0.5)) / numU);
					const Real v = std::lerp(y1, y2, (Real(j) + Real(0.5)) / numV);
					const auto point = sample(u, v);
					const auto partialU = BoundedPartial([&](Real t) { return sample(t, v); }, point, u, x1, x2);
					const auto partialV = BoundedPartial([&](Real t) { return sample(u, t); }, point, v, y1, y2);
					const auto field = vectorField(point);
					total.Add(ScaledFlux(field, partialU, partialV, halfU, halfV, Real(4) / patches));
				}
			}

			IntegrationDetail::RequireFinite(total.Value(), "SurfaceIntegral");
			IntegrationResult result(total.Value(), NoEstimate, patches, false, false);
			result.function_evaluations = patches;
			result.algorithm_name = "SurfaceIntegralMidpoint";
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}

		static Real SurfaceIntegral(const IVectorFunction<3>& vectorField, const IParametricSurfaceRect<3>& surface, 
		                            const Real x1, const Real x2, const Real y1, const Real y2, 
		                            int numU = 20, int numV = 20)
		{
			return SurfaceIntegralDetailed(vectorField, surface, x1, x2, y1, y2, numU, numV).ValueOrThrow("SurfaceIntegral");
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
		/// Real flux = SurfaceIntegration::SurfaceIntegralDetailed(F, tetra);
		/// // Should equal 3 * volume of tetrahedron
		/// @endcode
		static IntegrationResult SurfaceIntegralDetailed(const IVectorFunction<3>& vectorField, const BodyWithTriangleSurfaces& solid, Real eps = Precision::SurfaceIntegrationTolerance)
		{
			std::vector<Triangle3D> faces;
			for (int i = 0; i < solid.getSurfaceCount(); ++i) faces.push_back(solid.getSurface(i));
			return IntegrateSurfaces(vectorField, faces, eps, Precision::SurfaceIntegrationMaxDepth);
		}

		static Real SurfaceIntegral(const IVectorFunction<3>& vectorField, const BodyWithTriangleSurfaces& solid, Real eps = Precision::SurfaceIntegrationTolerance)
		{
			return SurfaceIntegralDetailed(vectorField, solid, eps).ValueOrThrow("SurfaceIntegral");
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
		/** @par Subdivision pattern:
		 * @code
		 *       A                    A
		 *      / \                  /|\
		 *     /   \     --->      m1-+-m3
		 *    /     \              |\ | /|
		 *   B-------C             B--m2--C
		 * @endcode
		 */
		/// where m1, m2, m3 are edge midpoints.
		static IntegrationResult SurfaceIntegralDetailed(const IVectorFunction<3>& vectorField, const Triangle3D& triangle, Real eps, int maxLevel = Precision::SurfaceIntegrationMaxDepth)
		{
			return IntegrateSurfaces(vectorField, std::vector<Triangle3D>{triangle}, eps, maxLevel);
		}

		static Real SurfaceIntegral(const IVectorFunction<3>& vectorField, const Triangle3D& triangle, Real eps, int maxLevel = Precision::SurfaceIntegrationMaxDepth)
		{
			return SurfaceIntegralDetailed(vectorField, triangle, eps, maxLevel).ValueOrThrow("SurfaceIntegral");
		}

	private:
		static void RequireFiniteVector(const VectorN<Real, 3>& value, const char* context) {
			for (int i = 0; i < 3; ++i)
				IntegrationDetail::RequireFinite(value[i], context);
		}

		template<class Sample>
		static VectorN<Real, 3> BoundedPartial(Sample sample, const VectorN<Real, 3>& center,
			Real t, Real lower, Real upper) {
			const Real step = Derivation::ScaleStep(PrecisionValues<Real>::DerivativeStepSize, t);
			const auto derivative = Derivation::Detail::BoundedFirstDerivative(sample, center, t, lower, upper, step);
			RequireFiniteVector(derivative, "Surface derivative");
			return derivative;
		}

		static Real ScaledFlux(const VectorN<Real, 3>& field, const VectorN<Real, 3>& first,
			const VectorN<Real, 3>& second, Real scaleU = 1, Real scaleV = 1, Real fraction = 1) {
			RequireFiniteVector(field, "Surface integrand");
			RequireFiniteVector(first, "Surface tangent");
			RequireFiniteVector(second, "Surface tangent");
			IntegrationDetail::Sum flux;
			for (int i = 0; i < 3; ++i) {
				const int j = (i + 1) % 3, k = (i + 2) % 3;
				// Retain the field and quadrature weights before reconstructing each product's exponent.
				flux.Add(IntegrationDetail::Product({field[i], first[j], second[k], scaleU, scaleV, fraction}));
				flux.Add(-IntegrationDetail::Product({field[i], first[k], second[j], scaleU, scaleV, fraction}));
			}
			IntegrationDetail::RequireFinite(flux.Value(), "Surface flux");
			return flux.Value();
		}

		/// @brief Subdivide a triangle while retaining its orientation.
		static std::vector<Triangle3D> SplitSurface(const Triangle3D& triangle)
		{
			// Subdivide triangle into 4 by connecting edge midpoints
			Pnt3Cart m12 = (triangle.Pnt1() + triangle.Pnt2()) / 2.0;
			Pnt3Cart m23 = (triangle.Pnt2() + triangle.Pnt3()) / 2.0;
			Pnt3Cart m31 = (triangle.Pnt3() + triangle.Pnt1()) / 2.0;

			Triangle3D t1(triangle.Pnt1(), m12, m31);  // Corner triangle at vertex 1
			Triangle3D t2(m12, triangle.Pnt2(), m23);  // Corner triangle at vertex 2
			Triangle3D t3(m31, m23, triangle.Pnt3());  // Corner triangle at vertex 3
			Triangle3D t4(m12, m23, m31);              // Central triangle

			return {t1, t2, t3, t4};
		}

		/// @brief Calculate flux contribution of a single triangle using centroid rule
		/// @details Evaluates F·n at triangle center, multiplies by triangle area
		static Real CalcSurfaceContrib(const IVectorFunction<3>& vectorField, const Triangle3D& triangle, int* evaluations = nullptr)
		{
			Vec3Cart edge1(triangle.Pnt1(), triangle.Pnt2());
			Vec3Cart edge2(triangle.Pnt1(), triangle.Pnt3());
			// Centroid (center of mass)
			Pnt3Cart center = (triangle.Pnt1() + triangle.Pnt2() + triangle.Pnt3()) / 3.0;

			// Evaluate field and compute dot product
			if (evaluations) {
				if (*evaluations == std::numeric_limits<int>::max()) throw ResourceLimitError("Surface integration counter exhausted");
				++*evaluations;
			}
			Vector3Cartesian value = vectorField(Vector3Cartesian(center.X(), center.Y(), center.Z()));
			return ScaledFlux(value, edge1, edge2, 1, 1, Real(0.5));
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
		/// Real flux = SurfaceIntegration::SurfaceIntegralDetailed(F, box);
		/// // Net flux = 0 (same flux in through left face, out through right)
		/// @endcode
		static IntegrationResult SurfaceIntegralDetailed(const IVectorFunction<3>& vectorField, const BodyWithRectSurfaces& solid, Real eps = Precision::SurfaceIntegrationTolerance)
		{
			std::vector<RectSurface3D> faces;
			for (int i = 0; i < solid.getSurfaceCount(); ++i) faces.push_back(solid.getSurface(i));
			return IntegrateSurfaces(vectorField, faces, eps, Precision::SurfaceIntegrationMaxDepth);
		}

		static Real SurfaceIntegral(const IVectorFunction<3>& vectorField, const BodyWithRectSurfaces& solid, Real eps = Precision::SurfaceIntegrationTolerance)
		{
			return SurfaceIntegralDetailed(vectorField, solid, eps).ValueOrThrow("SurfaceIntegral");
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
		/// @return IntegrationResult with flux value (use .value; check converged)
		static IntegrationResult SurfaceIntegralDetailed(const IVectorFunction<3>& vectorField, const RectSurface3D& surface, Real eps, int maxLevel = Precision::SurfaceIntegrationMaxDepth)
		{
			return IntegrateSurfaces(vectorField, std::vector<RectSurface3D>{surface}, eps, maxLevel);
		}

		static Real SurfaceIntegral(const IVectorFunction<3>& vectorField, const RectSurface3D& surface, Real eps, int maxLevel = Precision::SurfaceIntegrationMaxDepth)
		{
			return SurfaceIntegralDetailed(vectorField, surface, eps, maxLevel).ValueOrThrow("SurfaceIntegral");
		}

	private:
		/// @brief Subdivide a rectangular patch while retaining its orientation.
		static std::vector<RectSurface3D> SplitSurface(const RectSurface3D& surface)
		{
			// Compute midpoints of edges and center
			Point3Cartesian pnt_mid12 = detail::GeometryMidpoint(surface._pnt1, surface._pnt2);
			Point3Cartesian pnt_mid23 = detail::GeometryMidpoint(surface._pnt2, surface._pnt3);
			Point3Cartesian pnt_mid34 = detail::GeometryMidpoint(surface._pnt3, surface._pnt4);
			Point3Cartesian pnt_mid41 = detail::GeometryMidpoint(surface._pnt4, surface._pnt1);
			Point3Cartesian center = surface.getCenter();

			// Create 4 sub-rectangles
			RectSurface3D s1(surface._pnt1, pnt_mid12, center, pnt_mid41);
			RectSurface3D s2(pnt_mid12, surface._pnt2, pnt_mid23, center);
			RectSurface3D s3(center, pnt_mid23, surface._pnt3, pnt_mid34);
			RectSurface3D s4(pnt_mid41, center, pnt_mid34, surface._pnt4);

			// Calculate contributions for each subsurface
			return {s1, s2, s3, s4};
		}

		/// @brief Calculate flux contribution of a single rectangle using center rule
		/// @details Evaluates F·n at rectangle center, multiplies by rectangle area
		static Real CalcSurfaceContrib(const IVectorFunction<3>& vectorField, const RectSurface3D& surface, int* evaluations = nullptr)
		{
			Vector3Cartesian normal = surface.getNormal();
			Real area = surface.getArea();
			Point3Cartesian center = surface.getCenter();

			if (evaluations) {
				if (*evaluations == std::numeric_limits<int>::max()) throw ResourceLimitError("Surface integration counter exhausted");
				++*evaluations;
			}
			Vector3Cartesian value = vectorField(Vector3Cartesian(center.X(), center.Y(), center.Z()));
			for (int i = 0; i < 3; ++i) IntegrationDetail::RequireFinite(value[i], "Surface integrand");
			Real dotProduct = ScalarProduct(normal, value);

			return dotProduct * area;
		}

		template<class Surface>
		static IntegrationResult IntegrateSurfaces(const IVectorFunction<3>& field, const std::vector<Surface>& surfaces,
			Real tolerance, int maxDepth) {
			IntegrationResult result;
			result.algorithm_name = "SurfaceIntegral";
			AlgorithmTimer timer;
			auto evaluate = [&](const Surface& surface) {
				const Real coarse = CalcSurfaceContrib(field, surface, &result.function_evaluations);
				IntegrationDetail::Sum fine;
				for (const auto& child : SplitSurface(surface))
					fine.Add(CalcSurfaceContrib(field, child, &result.function_evaluations));
				IntegrationResult estimate(fine.Value(), std::abs(fine.Value() - coarse), 1, false);
				estimate.status = AlgorithmStatus::Completed;
				return estimate;
			};
			auto split = [](const Surface& surface) { return SplitSurface(surface); };
			IntegrationDetail::RefineCells(result, surfaces, evaluate, split, tolerance, 0, maxDepth,
				std::numeric_limits<int>::max(), 5);
			if (result.function_evaluations > 0) ++result.iterations;
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}
	};

} // end namespace MML

#endif // MML_SURFACE_INTEGRATION_H
