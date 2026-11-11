///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Surfaces.h                                                          ///
///  Description: Parametric surface classes (2D parameter domain to 3D)              ///
///               Normal vectors, surface area, Gaussian/mean curvature               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_CURVES_SURFACES_H
#define MML_CURVES_SURFACES_H

#include <mml/MMLBase.h>

#include <mml/base/Vector/VectorN.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Vector/VectorTypes3D.h>
#include <mml/base/BaseUtils.h>

#include <mml/base/Function.h>
#include <mml/base/Geometry/Geometry3D.h>

#include <mml/core/Derivation.h>

namespace MML
{
	namespace Surfaces
	{
		/// @brief Base class for parametric surfaces r(u,w) : R² → R³ in Cartesian coordinates
		/// @note Provides normal vectors, tangents, curvatures (Gaussian/mean), principal directions
		/// @note Uses first/second fundamental forms (E,F,G) and (L,M,N) for curvature calculations
		class ISurfaceCartesian : public IParametricSurfaceRect<3>	
		{
		public:
			/// @brief Compute unit normal vector n = (∂r/∂u × ∂r/∂w) / ||∂r/∂u × ∂r/∂w||
			/// @param u,w Surface parameters
			/// @return Normalized normal vector (zero if degenerate)
			virtual VectorN<Real, 3> Normal(Real u, Real w) const
			{
				// Compute partial derivatives
				Vec3Cart r_u = Derivation::NDer1_u(*this, u, w);
				Vec3Cart r_w = Derivation::NDer1_w(*this, u, w);

			// Compute cross product
			VectorN<Real, 3> n = VectorProduct(r_u, r_w);

			// Normalize
			Real norm = n.NormL2();
			constexpr Real eps = PrecisionValues<Real>::NumericalZeroThreshold;
			if (norm < eps)
				return VectorN<Real, 3>{0.0, 0.0, 0.0}; // Degenerate case

			return n / norm;
			}

			/// @brief Compute tangent vectors along parameter directions
			/// @param u,w Surface parameters
			/// @param tU Output: tangent in u direction (∂r/∂u)
			/// @param tW Output: tangent in w direction (∂r/∂w)
			virtual void Tangents(Real u, Real w, VectorN<Real, 3>& tU, VectorN<Real, 3>& tW) const
			{
				// Compute tangent vector in the u direction
				tU = Derivation::NDer1_u(*this, u, w);
				// Compute tangent vector in the w direction
				tW = Derivation::NDer1_w(*this, u, w);
			}

			/// @brief Compute principal curvatures κ₁ and κ₂ (eigenvalues of shape operator)
			/// @param u,w Surface parameters
			/// @param k1 Output: maximum principal curvature
			/// @param k2 Output: minimum principal curvature
			/// @note Solves κ² - 2H·κ + K = 0 where H=mean curvature, K=Gaussian curvature
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

			/// @brief Compute principal direction vectors (eigenvectors of shape operator)
			/// @param u,w Surface parameters
			/// @param dir1 Output: direction of maximum curvature κ₁
			/// @param dir2 Output: direction of minimum curvature κ₂
			/// @note Solves (II - κ·I) v = 0 for v = (a,b) in tangent space
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

			/// @brief Get first fundamental form coefficients (E,F,G) - measures surface metric
			/// @param u,w Surface parameters
			/// @param outE Output: E = ∂r/∂u · ∂r/∂u
			/// @param outF Output: F = ∂r/∂u · ∂r/∂w
			/// @param outG Output: G = ∂r/∂w · ∂r/∂w
			/// @note ds² = E du² + 2F du dw + G dw² (arc length formula)
			virtual void		 GetFirstNormalFormCoefficients(Real u, Real w, Real& outE, Real& outF, Real& outG) const
			{
				// first, we need to calculate necessary derivatives
				VectorN<Real, 3> r_by_u = Derivation::NDer1_u(*this, u, w);
				VectorN<Real, 3> r_by_w = Derivation::NDer1_w(*this, u, w);

				outE = Utils::ScalarProduct(r_by_u, r_by_u);
				outF = Utils::ScalarProduct(r_by_u, r_by_w);
				outG = Utils::ScalarProduct(r_by_w, r_by_w);
			}
			/// @brief Get first fundamental form as vector (E,F,G)
			/// @return Vector with [E, F, G] coefficients
			VectorN<Real, 3> GetFirstNormalFormCoefficients(Real u, Real w) const
			{
				VectorN<Real, 3> ret;

				GetFirstNormalFormCoefficients(u, w, ret[0], ret[1], ret[2]);
				
				return ret;
			}
			/// @brief Get second fundamental form coefficients (L,M,N) - measures surface curvature
			/// @param u,w Surface parameters
			/// @param outL Output: L = ∂²r/∂u² · n
			/// @param outM Output: M = ∂²r/∂u∂w · n
			/// @param outN Output: N = ∂²r/∂w² · n
			/// @note Shape operator S = I⁻¹·II where I=(E,F,G), II=(L,M,N)
			virtual void		 GetSecondNormalFormCoefficients(Real u, Real w, Real &outL, Real &outM, Real &outN) const
			{
				// first, we need to calculate necessary derivatives
				Vec3Cart r_by_u = Derivation::NDer1_u(*this, u, w);
				Vec3Cart r_by_w = Derivation::NDer1_w(*this, u, w);
				
				// normal — check for degenerate point (parallel tangents / cusp / fold)
				Vec3Cart cross = VectorProduct(r_by_u, r_by_w);
				Real mag = cross.NormL2();
				constexpr Real eps = PrecisionValues<Real>::NumericalZeroThreshold;
				if (mag < eps) {
					// Degenerate point: second fundamental form is undefined
					outL = outM = outN = 0.0;
					return;
				}
				Vec3Cart n = cross / mag;

				// now we can calculate second derivatives
				Vec3Cart r_by_uu = Derivation::NDer2_uu(*this, u, w);
				Vec3Cart r_by_uw = Derivation::NDer2_uw(*this, u, w);
				Vec3Cart r_by_ww = Derivation::NDer2_ww(*this, u, w);

				outL = Utils::ScalarProduct(r_by_uu, n); // L
				outM = Utils::ScalarProduct(r_by_uw, n); // M
				outN = Utils::ScalarProduct(r_by_ww, n); // N
			}
			/// @brief Get second fundamental form as vector (L,M,N)
			/// @return Vector with [L, M, N] coefficients
			VectorN<Real, 3> GetSecondNormalFormCoefficients(Real u, Real w) const
			{
				VectorN<Real, 3> ret;

				GetSecondNormalFormCoefficients(u, w, ret[0], ret[1], ret[2]);

				return ret;
			}

			/// @brief Compute Gaussian curvature K = (LN - M²)/(EG - F²) = κ₁·κ₂
			/// @param u,w Surface parameters
			/// @return K (positive: elliptic, zero: parabolic/flat, negative: hyperbolic)
			/// @note Intrinsic invariant (preserved under isometries)
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

			/// @brief Compute mean curvature H = (EN + GL - 2FM)/(2(EG - F²)) = (κ₁+κ₂)/2
			/// @param u,w Surface parameters
			/// @return H (zero for minimal surfaces)
			/// @note Extrinsic invariant (depends on embedding in R³)
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
				Vec3Cart tU = Derivation::NDer1_u(*this, u, w);
				Vec3Cart tW = Derivation::NDer1_w(*this, u, w);

				// Compute normal vector (cross product)
				VectorN<Real, 3> n = VectorProduct(tU, tW);

				// If the norm is above a small threshold, the surface is regular at (u, w)
				constexpr Real eps = PrecisionValues<Real>::NumericalZeroThreshold;
				return n.NormL2() > eps;
			}
			bool isFlat(Real u, Real w, Real eps = PrecisionValues<Real>::DefaultToleranceStrict)
			{
				return std::abs(GaussianCurvature(u, w)) < eps;
			}
			bool isParabolic(Real u, Real w, Real eps = PrecisionValues<Real>::DefaultToleranceStrict)
			{
				return std::abs(GaussianCurvature(u, w)) < eps && std::abs(MeanCurvature(u, w)) > eps;
			}
			bool isHyperbolic(Real u, Real w, Real eps = PrecisionValues<Real>::DefaultToleranceStrict)
			{
				return GaussianCurvature(u, w) < -eps;
			}
		};

		/// @brief Plane surface r(u,w) = point + u·uAxis + w·vAxis with local (u,v) coordinate system
		/// @note Flat surface: K=0, H=0 everywhere
		/// @note Supports optional limits (minU,maxU,minV,maxV) for bounded plane regions
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

		/// @brief Cylindrical surface r(u,w) = (R·cos(u), R·sin(u), w), radius R, height H
		/// @note u∈[0,2π], w∈[0,H]. Parabolic surface: K=0 (one principal curvature zero)
		class CylinderSurface : public ISurfaceCartesian
		{
			Real _radius, _height;
		public:
			CylinderSurface() : _radius(1), _height(1) {}
			CylinderSurface(Real R, Real H) : _radius(R), _height(H) {}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return 2 * Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return _height; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return VectorN<Real, 3>{static_cast<Real>(_radius* cos(u)), static_cast<Real>(_radius* sin(u)), w}; }
		};
		
		/// @brief Torus r(u,w) = ((R+r·cos(w))·cos(u), (R+r·cos(w))·sin(u), r·sin(w))
		/// @note R=major radius, r=minor radius. u,w∈[0,2π]. K varies (elliptic inner, hyperbolic outer)
		class Torus : public ISurfaceCartesian
		{
			Real _majorRadius, _minorRadius;
		public:
			Torus() : _majorRadius(1), _minorRadius(0.5) {}
			Torus(Real R, Real r) : _majorRadius(R), _minorRadius(r) {}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return 2 * Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{static_cast<Real>((_majorRadius + _minorRadius * cos(w)) * cos(u)), static_cast<Real>((_majorRadius + _minorRadius * cos(w)) * sin(u)), static_cast<Real>(_minorRadius * sin(w))}; }
		};

		/// @brief Spherical surface r(u,w) = R·(sin(u)·cos(w), sin(u)·sin(w), cos(u)), radius R
		/// @note u∈[0,π] (latitude), w∈[0,2π] (longitude). Constant positive curvature K=1/R²
		class Sphere : public ISurfaceCartesian
		{
			Real _radius;
		public:
			Sphere() : _radius(1) {}
			Sphere(Real R) : _radius(R) {}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{static_cast<Real>(_radius * sin(u)* cos(w)), static_cast<Real>(_radius * sin(u)* sin(w)), static_cast<Real>(_radius * cos(u))}; }
		};





		/// @brief Ellipsoid r(u,w) = (a·sin(u)·cos(w), b·sin(u)·sin(w), c·cos(u)), semi-axes (a,b,c)
		/// @note u∈[0,π], w∈[0,2π]. Generalization of sphere (a=b=c gives sphere)
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

		VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{static_cast<Real>(_a * sin(u) * cos(w)), static_cast<Real>(_b * sin(u) * sin(w)), static_cast<Real>(_c * cos(u))}; }
		};

		/// @brief One-sheet hyperboloid r(u,w) = (a·cosh(u)·cos(w), b·cosh(u)·sin(w), c·sinh(u))
		/// @note Ruled surface with negative Gaussian curvature. u∈[0,π], w∈[0,2π]
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

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{static_cast<Real>(_a * cosh(u) * cos(w)), static_cast<Real>(_b * cosh(u) * sin(w)), static_cast<Real>(_c * sinh(u))}; }
		};

		/// @brief Paraboloid of revolution z = u, r(u,w) = (a·√(u/h)·cos(w), a·√(u/h)·sin(w), u)
		/// @note a=radius scaling, h=height scaling. Parabolic shape with positive curvature
		class Paraboloid : public ISurfaceCartesian
		{
			Real _a, _h;
		public:
			Paraboloid() : _a(1), _h(1) {}
			Paraboloid(Real a, Real h) : _a(a), _h(h) {}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return 10; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{static_cast<Real>(_a * sqrt(u/_h) * cos(w)), static_cast<Real>(_a * sqrt(u/_h) * sin(w)), u}; }
		};








	} // end namespace Surfaces
}

#endif // MML_CURVES_SURFACES_H