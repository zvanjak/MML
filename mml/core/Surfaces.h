///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Surfaces.h                                                          ///
///  Description: Parametric surface classes (2D parameter domain to 3D)              ///
///               Normal vectors, surface area, Gaussian/mean curvature               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_CURVES_SURFACES_H
#define MML_CURVES_SURFACES_H

#include "MMLBase.h"

#include "base/VectorN.h"
#include "base/Vector.h"
#include "base/VectorTypes.h"
#include "base/BaseUtils.h"

#include "base/Function.h"
#include "base/Geometry3D.h"

#include "core/Derivation.h"

using namespace MML::Utils;
using namespace MML::Derivation;

namespace MML
{
	namespace Surfaces
	{
		class ISurfaceCartesian : public IParametricSurfaceRect<3>	
		{
		public:
			virtual VectorN<Real, 3> Normal(Real u, Real w) const
			{
				// Compute partial derivatives
				Vec3Cart r_u = NDer1_u(*this, u, w);
				Vec3Cart r_w = NDer1_w(*this, u, w);

			// Compute cross product
			VectorN<Real, 3> n = VectorProduct(r_u, r_w);

			// Normalize
			Real norm = n.NormL2();
			constexpr Real eps = PrecisionValues<Real>::NumericalZeroThreshold;
			if (norm < eps)
				return VectorN<Real, 3>{0.0, 0.0, 0.0}; // Degenerate case

			return n / norm;
			}

			virtual void Tangents(Real u, Real w, VectorN<Real, 3>& tU, VectorN<Real, 3>& tW) const
			{
				// Compute tangent vector in the u direction
				tU = NDer1_u(*this, u, w);
				// Compute tangent vector in the w direction
				tW = NDer1_w(*this, u, w);
			}

			virtual void PrincipalCurvatures(Real u, Real w, Real& k1, Real& k2) const
			{
				// Get first and second fundamental form coefficients
				Real E, F, G;
				GetFirstNormalFormCoefficients(u, w, E, F, G);
				Real L, M, N;
				GetSecondNormalFormCoefficients(u, w, L, M, N);

				// Compute the invariants
				Real denom = E * G - F * F;
				constexpr Real eps = PrecisionValues<Real>::NumericalZeroThreshold;
				if (std::abs(denom) < eps) {
					k1 = k2 = 0.0;
					return;
				}

				// Trace and determinant of the shape operator
				Real H = (E * N + G * L - 2 * F * M) / (2 * denom); // mean curvature
				Real K = (L * N - M * M) / denom;                   // Gaussian curvature

				// Solve quadratic equation: k^2 - 2H k + K = 0
				Real discr = H * H - K;
				if (discr < 0) discr = 0; // Clamp for numerical stability

				Real sqrt_discr = std::sqrt(discr);
				k1 = H + sqrt_discr;
				k2 = H - sqrt_discr;
			}

			virtual void PrincipalDirections(Real u, Real w, VectorN<Real, 3>& dir1, VectorN<Real, 3>& dir2) const
			{
				// Get first and second fundamental form coefficients
				Real E, F, G;
				GetFirstNormalFormCoefficients(u, w, E, F, G);
				Real L, M, N;
				GetSecondNormalFormCoefficients(u, w, L, M, N);

				// Compute the invariants
				Real denom = E * G - F * F;
				constexpr Real eps = PrecisionValues<Real>::NumericalZeroThreshold;
				if (std::abs(denom) < eps) {
					dir1 = dir2 = VectorN<Real, 3>{ 0.0, 0.0, 0.0 };
					return;
				}

				// Get principal curvatures
				Real k1, k2;
				PrincipalCurvatures(u, w, k1, k2);

				// Tangent vectors
				VectorN<Real, 3> tU, tW;
				Tangents(u, w, tU, tW);

				// For each principal curvature, solve (II - k*I) v = 0 for v = (a, b)
				// That is:
				// [L - k*E, M - k*F] [a] = 0
				// [M - k*F, N - k*G] [b]   0

				auto getDirection = [&](Real k) -> VectorN<Real, 3> {
					Real a = M - k * F;
					Real b = L - k * E;
					if (std::abs(a) > std::abs(b)) {
						// Set v = (1, -b/a)
						Real v0 = 1.0;
						Real v1 = (a != 0.0) ? -b / a : 0.0;
						return (v0 * tU + v1 * tW).Normalized();
					}
					else if (std::abs(b) > eps) {
						// Set v = (-a/b, 1)
						Real v0 = (b != 0.0) ? -a / b : 0.0;
						Real v1 = 1.0;
						return (v0 * tU + v1 * tW).Normalized();
					}
					else {
						// Fallback: use tU
						return tU.Normalized();
					}
					};

				dir1 = getDirection(k1);
				dir2 = getDirection(k2);
			}

			virtual void		 GetFirstNormalFormCoefficients(Real u, Real w, Real& outE, Real& outF, Real& outG) const
			{
				// first, we need to calculate necessary derivatives
				VectorN<Real, 3> r_by_u = NDer1_u(*this, u, w);
				VectorN<Real, 3> r_by_w = NDer1_w(*this, u, w);

				outE = ScalarProduct(r_by_u, r_by_u);
				outF = ScalarProduct(r_by_u, r_by_w);
				outG = ScalarProduct(r_by_w, r_by_w);
			}
			VectorN<Real, 3> GetFirstNormalFormCoefficients(Real u, Real w) const
			{
				VectorN<Real, 3> ret;

				GetFirstNormalFormCoefficients(u, w, ret[0], ret[1], ret[2]);
				
				return ret;
			}
			virtual void		 GetSecondNormalFormCoefficients(Real u, Real w, Real &outL, Real &outM, Real &outN) const
			{
				// first, we need to calculate necessary derivatives
				Vec3Cart r_by_u = NDer1_u(*this, u, w);
				Vec3Cart r_by_w = NDer1_w(*this, u, w);
				
			// normal
			Vec3Cart n = VectorProduct(r_by_u, r_by_w).Normalized();				// now we can calculate second derivatives
				Vec3Cart r_by_uu = NDer2_uu(*this, u, w);
				Vec3Cart r_by_uw = NDer2_uw(*this, u, w);
				Vec3Cart r_by_ww = NDer2_ww(*this, u, w);

				outL = ScalarProduct(r_by_uu, n); // L
				outM = ScalarProduct(r_by_uw, n); // M
				outN = ScalarProduct(r_by_ww, n); // N
			}
			VectorN<Real, 3> GetSecondNormalFormCoefficients(Real u, Real w) const
			{
				VectorN<Real, 3> ret;

				GetSecondNormalFormCoefficients(u, w, ret[0], ret[1], ret[2]);

				return ret;
			}

			Real GaussianCurvature(Real u, Real w)
			{
				Real E, F, G;
				GetFirstNormalFormCoefficients(u, w, E, F, G);

				Real L, M, N;
				GetSecondNormalFormCoefficients(u, w, L, M, N);

				// Compute denominator and check for degeneracy
				Real denom = E * G - F * F;
				constexpr Real eps = PrecisionValues<Real>::NumericalZeroThreshold;
				if (std::abs(denom) < eps)
				{
					// Degenerate metric: curvature undefined
					return 0.0;
					// Alternatively: return std::numeric_limits<Real>::quiet_NaN();
				}

				// Gaussian curvature formula: K = (L*N - M^2) / (E*G - F^2)
				return (L * N - M * M) / denom;
			}

			Real MeanCurvature(Real u, Real w)
			{
				Real E, F, G;
				GetFirstNormalFormCoefficients(u, w, E, F, G);

				Real L, M, N;
				GetSecondNormalFormCoefficients(u, w, L, M, N);

				// Compute denominator and check for degeneracy
				Real denom = 2 * (E * G - F * F);
				constexpr Real eps = PrecisionValues<Real>::NumericalZeroThreshold;
				if (std::abs(denom) < eps)
				{
					// Degenerate metric: mean curvature undefined
					return 0.0;
					// Alternatively: return std::numeric_limits<Real>::quiet_NaN();
				}

				// Mean curvature formula: H = (E*N + G*L - 2*F*M) / (2*(E*G - F^2))
				return (E * N + G * L - 2 * F * M) / denom;
			}

			bool isRegular(Real u, Real w) const
			{
				// Compute tangent vectors
				Vec3Cart tU = NDer1_u(*this, u, w);
				Vec3Cart tW = NDer1_w(*this, u, w);

				// Compute normal vector (cross product)
				VectorN<Real, 3> n = VectorProduct(tU, tW);

				// If the norm is above a small threshold, the surface is regular at (u, w)
				constexpr Real eps = PrecisionValues<Real>::NumericalZeroThreshold;
				return n.NormL2() > eps;
			}
			bool isFlat(Real u, Real w, Real eps = PrecisionValues<Real>::DefaultToleranceRelaxed)
			{
				return std::abs(GaussianCurvature(u, w)) < eps;
			}
			bool isParabolic(Real u, Real w, Real eps = PrecisionValues<Real>::DefaultToleranceRelaxed)
			{
				return std::abs(GaussianCurvature(u, w)) < eps && !isFlat(u, w, eps);
			}
			bool isHyperbolic(Real u, Real w, Real eps = PrecisionValues<Real>::DefaultToleranceRelaxed)
			{
				return GaussianCurvature(u, w) < -eps;
			}
		};

		// plane surface, given by point and normal, with local (u, v) coordinate system
		class PlaneSurface : public ISurfaceCartesian
		{
			VectorN<Real, 3> _point, _normal;
			VectorN<Real, 3> _uAxis, _vAxis;

			Real _minU = Constants::NegInf, _maxU = Constants::PosInf;
			Real _minV = Constants::NegInf, _maxV = Constants::PosInf;

		public:
			PlaneSurface(const Vec3Cart& point, const Vec3Cart& normal, const Vec3Cart& uAxis, const Vec3Cart& vAxis,
									 Real minU, Real maxU,	Real minV, Real maxV) 
				: PlaneSurface(point, normal, uAxis, vAxis)
			{
				SetPlaneLimits(minU, maxU, minV, maxV);
			}
			// complete ctor, with all info given
			PlaneSurface(const Vec3Cart& point, const Vec3Cart& normal, const Vec3Cart& uAxis, const Vec3Cart& vAxis)
				: _point(point), _normal(normal), _uAxis(uAxis), _vAxis(vAxis)
			{
				// check that the uAxis and vAxis are orthogonal to the normal
				if (!normal.IsPerpendicularTo(uAxis) || !normal.IsPerpendicularTo(vAxis))
				{
					throw GeometryError("uAxis and vAxis must be orthogonal to the normal.");
				}
			}
			// this ctor gets only one local axis direction, the other is computed as cross product with normal
			PlaneSurface(const Vec3Cart& point, const Vec3Cart& normal, const Vec3Cart& uAxisDirection)
			{
				// first we have to make sure that the uAxisDirection is not parallel to the normal
				if (normal.IsParallelTo(uAxisDirection))
				{
					throw GeometryError("uAxisDirection must not be parallel to the normal.");
				}
				// then we have to project uAxisDirection onto the plane defined by point and normal
				_point = point;
				_normal = normal.Normalized();
				Vec3Cart uAxis = uAxisDirection - (uAxisDirection.ScalarProduct(_normal) * _normal);

				// now we can compute the vAxis as cross product of normal and uAxis
				_vAxis = VectorProduct(_normal, uAxis).Normalized();
			}
			// this ctor only gets point and normal, uAxis and vAxis are computed as orthogonal to normal
			// with uAxis trying to align as much as possible with x axis
			PlaneSurface(const Vec3Cart& point, const Vec3Cart& normal)
			{
				_point = point;
				_normal = normal.Normalized();

				// try to align uAxis with x axis, if normal is not parallel to x axis
				if (!normal.IsParallelTo(Vec3Cart(1,0,0)))
				{
					_uAxis = VectorProduct(normal, Vec3Cart(1,0,0)).Normalized();
				}
				else // if normal is parallel to x axis, align uAxis with y axis
				{
					_uAxis = VectorProduct(normal, Vec3Cart(0,1,0)).Normalized();
				}
				_vAxis = VectorProduct(normal, _uAxis).Normalized();
			}

			Real getMinU() const { return _minU; }
			Real getMaxU() const { return _maxU; }
			Real getMinW() const { return _minV; }
			Real getMaxW() const { return _maxV; }

			void SetPlaneLimits(Real minU, Real maxU, Real minV, Real maxV)
			{
				_minU = minU;
				_maxU = maxU;
				_minV = minV;
				_maxV = maxV;
			}

			VectorN<Real, 3> operator()(Real u, Real w) const
			{
				return _point + u * _uAxis + w * _vAxis;
			}
		};

		class Cylinder : public ISurfaceCartesian
		{
			Real _R, _H;
		public:
			Cylinder() : _R(1), _H(1) {}
			Cylinder(Real R, Real H) : _R(R), _H(H) {}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return 2 * Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return _H; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return VectorN<Real, 3>{_R* cos(u), _R* sin(u), w}; }
		};
		
		class Torus : public ISurfaceCartesian
		{
			Real _R, _r;
		public:
			Torus() : _R(1), _r(0.5) {}
			Torus(Real R, Real r) : _R(R), _r(r) {}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return 2 * Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{(_R + _r * cos(w)) * cos(u), (_R + _r * cos(w)) * sin(u), _r * sin(w)}; }
		};

		class Sphere : public ISurfaceCartesian
		{
			Real _R;
		public:
			Sphere() : _R(1) {}
			Sphere(Real R) : _R(R) {}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{_R * sin(u)* cos(w), _R * sin(u)* sin(w), _R * cos(u)}; }
		};

		class MonkeySaddle : public ISurfaceCartesian
		{
		public:
			Real getMinU() const { return -10; }
			Real getMaxU() const { return 10; }
			Real getMinW() const { return -10; }
			Real getMaxW() const { return 10; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{u, w, u* (u * u - 3 * w * w)}; }
		};

		class MobiusStrip : public ISurfaceCartesian
		{
		public:
			Real getMinU() const { return 0; }
			Real getMaxU() const { return 2 * Constants::PI; }
			Real getMinW() const { return -1; }
			Real getMaxW() const { return 1; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{(1 + w * cos(u / 2)) * cos(u), (1 + w * cos(u / 2)) * sin(u), w * sin(u / 2)}; }
		};

		class Ellipsoid : public ISurfaceCartesian
		{
			Real _a, _b, _c;
		public:
			Ellipsoid() : _a(1), _b(1), _c(1) {}
			Ellipsoid(Real a, Real b, Real c) : _a(a), _b(b), _c(c) {}

		Real getMinU() const { return 0; }
		Real getMaxU() const { return Constants::PI; }
		Real getMinW() const { return 0; }
		Real getMaxW() const { return 2 * Constants::PI; }

		VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{_a * sin(u) * cos(w), _b * sin(u) * sin(w), _c * cos(u)}; }
		};

		class Hyperboloid : public ISurfaceCartesian
		{
			Real _a, _b, _c;
		public:
			Hyperboloid() : _a(1), _b(1), _c(1) {}
			Hyperboloid(Real a, Real b, Real c) : _a(a), _b(b), _c(c) {}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{_a * cosh(u) * cos(w), _b * cosh(u) * sin(w), _c * sinh(u)}; }
		};

		class Paraboloid : public ISurfaceCartesian
		{
			Real _a, _h;
		public:
			Paraboloid() : _a(1), _h(1) {}
			Paraboloid(Real a, Real h) : _a(a), _h(h) {}

			Real getMinU() const { return -10; }
			Real getMaxU() const { return 10; }
			Real getMinW() const { return -10; }
			Real getMaxW() const { return 10; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{_a * sqrt(u/_h) * cos(w), _a * sqrt(u/_h) * sin(w), u}; }
		};
	} // end namespace Surfaces
}

#endif // MML_CURVES_SURFACES_H