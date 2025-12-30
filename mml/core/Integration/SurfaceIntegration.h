///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        SurfaceIntegration.h                                                ///
///  Description: Surface integrals over parametric surfaces                          ///
///               Flux integrals and scalar surface integration                       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_SURFACE_INTEGRATION_H
#define MML_SURFACE_INTEGRATION_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/Vector.h"
#include "base/VectorN.h"
#include "base/Geometry3D.h"

#include "core/Derivation.h"
#include "core/Integration.h"
#include "core/FieldOperations.h"

namespace MML
{
	class SurfaceIntegration
	{
	public:
		// povrsinski integral preko
		// - ParametricSurface
		// - ParametricSurfaceRect
		// - SolidSurfaces3D -> 
		//		- SurfacesWithTriangles
		// 		- SurfacesWithRects
    // IParametricSurface integration on rectangular patch
	  // Uses finite differences to compute surface normals and integrates over parameter domain [x1,x2] x [y1,y2]
	  static Real SurfaceIntegral(const IVectorFunction<3>& vectorField, const IParametricSurfaceRect<3>& surface, 
                                const Real x1, const Real x2, const Real y1, const Real y2, 
                                int numU = 20, int numV = 20)
	  {
      if (x2 <= x1 || y2 <= y1)
        return 0.0;  // Invalid parameter range

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

          // Compute partial derivatives using finite differences
          Real h = PrecisionValues<Real>::DerivativeStepSize;  // Small step for numerical derivative

          VectorN<Real, 3> point_u_plus = surface(u + h, v);
          VectorN<Real, 3> point_u_minus = surface(u - h, v);
          VectorN<Real, 3> dS_du = (point_u_plus - point_u_minus) / (2.0 * h);

          VectorN<Real, 3> point_v_plus = surface(u, v + h);
          VectorN<Real, 3> point_v_minus = surface(u, v - h);
          VectorN<Real, 3> dS_dv = (point_v_plus - point_v_minus) / (2.0 * h);

          // Compute cross product: dS/du x dS/dv
          Vector3Cartesian partial_u(dS_du[0], dS_du[1], dS_du[2]);
          Vector3Cartesian partial_v(dS_dv[0], dS_dv[1], dS_dv[2]);
          Vector3Cartesian crossProd = VectorProduct(partial_u, partial_v);
          
          // Area element magnitude (not normalized yet)
          Real dS_magnitude = crossProd.NormL2();
          if (dS_magnitude < PrecisionValues<Real>::SurfaceNormalThreshold)
            continue;  // Degenerate patch

          // Unit normal
          Vector3Cartesian normal = crossProd / dS_magnitude;

          // Evaluate vector field at surface point
          Vector3Cartesian fieldValue = vectorField(Vector3Cartesian(point[0], point[1], point[2]));

          // Flux contribution: (F · n) * dS where dS = |∂S/∂u × ∂S/∂v| * du * dv
          Real dotProduct = ScalarProduct(normal, fieldValue);
          Real patchArea = dS_magnitude * du * dv;

          total += dotProduct * patchArea;
        }
      }

      return total;
		}

    ////////////////////////////////////////////////////////////////////////////////////////////
		// Surface integral over a solid with triangular surfaces (using Triangle3D)
		static Real SurfaceIntegral(const IVectorFunction<3>& vectorField, const BodyWithTriangleSurfaces& solid, Real eps = 0.001)
		{
			Real total = 0.0;
			for (int i = 0; i < solid.GetSurfaceCount(); i++)
				total += SurfaceIntegral(vectorField, solid.GetSurface(i), eps);
			return total;
		}

		// Surface integral over a single triangle (with recursive refinement)
		static Real SurfaceIntegral(const IVectorFunction<3>& vectorField, const Triangle3D& triangle, Real eps, int maxLevel = 7)
		{
			Real result = CalcSurfaceContrib(vectorField, triangle);

			return SurfaceIntegralImproveRecursively(vectorField, triangle, result, eps, maxLevel);
		}

		// Recursive refinement for triangle surface integration
		static Real SurfaceIntegralImproveRecursively(const IVectorFunction<3>& vectorField, const Triangle3D& triangle, Real prev_value, Real eps, int level)
		{
			// Subdivide triangle into 4 smaller triangles by connecting edge midpoints
			Pnt3Cart m12 = (triangle.Pnt1() + triangle.Pnt2()) / 2.0;
			Pnt3Cart m23 = (triangle.Pnt2() + triangle.Pnt3()) / 2.0;
			Pnt3Cart m31 = (triangle.Pnt3() + triangle.Pnt1()) / 2.0;

			Triangle3D t1(triangle.Pnt1(), m12, m31);
			Triangle3D t2(m12, triangle.Pnt2(), m23);
			Triangle3D t3(m31, m23, triangle.Pnt3());
			Triangle3D t4(m12, m23, m31);

			Real result = 0.0;
			result += CalcSurfaceContrib(vectorField, t1);
			result += CalcSurfaceContrib(vectorField, t2);
			result += CalcSurfaceContrib(vectorField, t3);
			result += CalcSurfaceContrib(vectorField, t4);

			if (fabs(result - prev_value) < eps || level == 0)
			{
				return result;
			}
			else
			{
				Real new_result = 0.0;
				new_result += SurfaceIntegralImproveRecursively(vectorField, t1, result, eps, level - 1);
				new_result += SurfaceIntegralImproveRecursively(vectorField, t2, result, eps, level - 1);
				new_result += SurfaceIntegralImproveRecursively(vectorField, t3, result, eps, level - 1);
				new_result += SurfaceIntegralImproveRecursively(vectorField, t4, result, eps, level - 1);
				return new_result;
			}
		}

		// Calculate the contribution of a single triangle (center-point rule)
		static Real CalcSurfaceContrib(const IVectorFunction<3>& vectorField, const Triangle3D& triangle)
		{
			// Compute normal using cross product of two edges
			Vec3Cart edge1(triangle.Pnt1(), triangle.Pnt2());
			Vec3Cart edge2(triangle.Pnt1(), triangle.Pnt3());
			Vec3Cart normal = VectorProduct(edge1, edge2).GetAsUnitVector();

			// Area using Heron's formula (already available as triangle.Area())
			Real area = triangle.Area();

			// Center (centroid) of the triangle
			Pnt3Cart center = (triangle.Pnt1() + triangle.Pnt2() + triangle.Pnt3()) / 3.0;

			Vector3Cartesian value = vectorField(Vector3Cartesian(center.X(), center.Y(), center.Z()));
			Real dotProduct = ScalarProduct(normal, value);

			return dotProduct * area;
		}

    ////////////////////////////////////////////////////////////////////////////////////////////
		// surface integral over a solid with rectangular surfaces
		static Real SurfaceIntegral(const IVectorFunction<3>& vectorField, const BodyWithRectSurfaces& solid, Real eps = 0.001)
		{
			Real total = 0.0;
			for (int i = 0; i < solid.GetSurfaceCount(); i++)
				total += SurfaceIntegral(vectorField, solid.GetSurface(i), eps);

			return total;
		}

		// Adaptive surface integration using eps for convergence, maxLevel for safety limit
		static Real SurfaceIntegral(const IVectorFunction<3>& vectorField, const RectSurface3D& surface, Real eps, int maxLevel = 7)
		{
			Real result = CalcSurfaceContrib(vectorField, surface);

			return SurfaceIntegralImproveRecursively(vectorField, surface, result, eps, maxLevel);
		}

		static Real SurfaceIntegralImproveRecursively(const IVectorFunction<3>& vectorField, const RectSurface3D& surface, Real prev_value, Real eps, int level)
		{
			// now we will calculate integral for surface divided to 4 equal parts
			Point3Cartesian pnt_mid12 = (surface._pnt1 + surface._pnt2) / 2.0;
			Point3Cartesian pnt_mid23 = (surface._pnt2 + surface._pnt3) / 2.0;
			Point3Cartesian pnt_mid34 = (surface._pnt3 + surface._pnt4) / 2.0;
			Point3Cartesian pnt_mid41 = (surface._pnt4 + surface._pnt1) / 2.0;
			Point3Cartesian center = surface.getCenter();

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

			// compare to prev_value - if converged or max level reached, return result
			if (fabs(result - prev_value) < eps || level == 0)
			{
				return result;
			}
			else
			{
				// Need more refinement - recursively improve each subsurface
				Real new_result = 0.0;

				new_result += SurfaceIntegralImproveRecursively(vectorField, s1, contrib1, eps, level - 1);
				new_result += SurfaceIntegralImproveRecursively(vectorField, s2, contrib2, eps, level - 1);
				new_result += SurfaceIntegralImproveRecursively(vectorField, s3, contrib3, eps, level - 1);
				new_result += SurfaceIntegralImproveRecursively(vectorField, s4, contrib4, eps, level - 1);

				return new_result;
			}
		}

		static Real CalcSurfaceContrib(const IVectorFunction<3>& vectorField, const RectSurface3D& surface)
		{
			Real result = 0.0;
			Vector3Cartesian    normal = surface.getNormal();
			Real area = surface.getArea();
			Point3Cartesian    center = surface.getCenter();

			Vector3Cartesian    value = vectorField(Vector3Cartesian(center.X(), center.Y(), center.Z()));

			Real dotProduct = ScalarProduct(normal, value);

			result = dotProduct * area;

			return result;
		}
	};
} // end namespace

#endif // MML_SURFACE_INTEGRATION_H
