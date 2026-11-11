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
#include <mml/core/detail/GeometryChecks.h>

namespace MML
{
	namespace Surfaces
	{
		/// @brief Base class for parametric surfaces r(u,w) : R² → R³ in Cartesian coordinates
		/// @note Provides normal vectors, tangents, curvatures (Gaussian/mean), principal directions
		/// @note Uses first/second fundamental forms (E,F,G) and (L,M,N) for curvature calculations
		class ISurfaceCartesian : public IParametricSurfaceRect<3>	
		{
			bool UnitTangents(Real u, Real w, Vec3Cart& first, Vec3Cart& second, Vec3Cart& normal) const
			{
				Tangents(u, w, first, second);
				const Real firstNorm = GeometryDetail::Norm(first), secondNorm = GeometryDetail::Norm(second);
				if (firstNorm == 0 || secondNorm == 0)
					return false;
				first = first / firstNorm;
				second = second / secondNorm;
				const Vec3Cart cross = VectorProduct(first, second);
				const Real sine = GeometryDetail::Norm(cross);
				if (sine * sine <= Precision::RoundoffZeroThreshold)
					return false;
				normal = cross / sine;
				return true;
			}

			bool ShapeOperator(Real u, Real w, Real& a, Real& b, Real& c) const
			{
				Real E, F, G;
				GetFirstNormalFormCoefficients(u, w, E, F, G);
				GeometryDetail::Checked(E);
				GeometryDetail::Checked(F);
				GeometryDetail::Checked(G);
				if (E < 0 || G < 0)
					throw GeometryError("Surface metric - negative squared tangent length");
				Real L, M, N;
				GetSecondNormalFormCoefficients(u, w, L, M, N);
				GeometryDetail::Checked(L);
				GeometryDetail::Checked(M);
				GeometryDetail::Checked(N);
				if (E == 0 || G == 0)
					return false;
				const Real rootE = std::sqrt(E), rootG = std::sqrt(G);
				const Real cosine = (F / rootE) / rootG;
				if (!std::isfinite(cosine) || std::abs(cosine) > 1 + Precision::RoundoffZeroThreshold)
					throw GeometryError("Surface metric - invalid normalized Gram matrix");
				const Real determinant = (1 - cosine) * (1 + cosine);
				if (determinant <= Precision::RoundoffZeroThreshold)
					return false;
				const Real l = GeometryDetail::Checked(L / E);
				const Real m = GeometryDetail::Checked((M / rootE) / rootG);
				const Real n = GeometryDetail::Checked(N / G);
				// Express I^-1 II in the orthonormal basis obtained from the two parameter tangents.
				a = l;
				b = GeometryDetail::Checked((m - cosine * l) / std::sqrt(determinant));
				c = GeometryDetail::Checked((n - 2 * cosine * m + cosine * cosine * l) / determinant);
				return true;
			}

		public:
			/// @brief Compute unit normal vector n = (∂r/∂u × ∂r/∂w) / ||∂r/∂u × ∂r/∂w||
			/// @param u,w Surface parameters
			/// @return Normalized normal vector (zero if degenerate)
			virtual VectorN<Real, 3> Normal(Real u, Real w) const
			{
				Vec3Cart first, second, normal;
				return UnitTangents(u, w, first, second, normal) ? normal : Vec3Cart{};
			}

			/// @brief Compute tangent vectors along parameter directions
			/// @param u,w Surface parameters
			/// @param tU Output: tangent in u direction (∂r/∂u)
			/// @param tW Output: tangent in w direction (∂r/∂w)
			virtual void Tangents(Real u, Real w, VectorN<Real, 3>& tU, VectorN<Real, 3>& tW) const
			{
				GeometryDetail::RequireFinite(u, "Surface tangent - u must be finite");
				GeometryDetail::RequireFinite(w, "Surface tangent - w must be finite");
				tU = GeometryDetail::Checked(Derivation::NDer1_u(*this, u, w));
				tW = GeometryDetail::Checked(Derivation::NDer1_w(*this, u, w));
			}

			/// @brief Compute principal curvatures κ₁ and κ₂ (eigenvalues of shape operator)
			/// @param u,w Surface parameters
			/// @param k1 Output: maximum principal curvature
			/// @param k2 Output: minimum principal curvature
			/// @note Solves κ² - 2H·κ + K = 0 where H=mean curvature, K=Gaussian curvature
			virtual void PrincipalCurvatures(Real u, Real w, Real& k1, Real& k2) const
			{
				Real a, b, c;
				if (!ShapeOperator(u, w, a, b, c)) {
					k1 = k2 = 0;
					return;
				}
				const Real mean = a / 2 + c / 2;
				const Real radius = std::hypot(a / 2 - c / 2, b);
				k1 = GeometryDetail::Checked(mean + radius);
				k2 = GeometryDetail::Checked(mean - radius);
			}

			/// @brief Compute principal direction vectors (eigenvectors of shape operator)
			/// @param u,w Surface parameters
			/// @param dir1 Output: direction of maximum curvature κ₁
			/// @param dir2 Output: direction of minimum curvature κ₂
			/// @note Solves (II - κ·I) v = 0 for v = (a,b) in tangent space
			virtual void PrincipalDirections(Real u, Real w, VectorN<Real, 3>& dir1, VectorN<Real, 3>& dir2) const
			{
				Real a, b, c;
				Vec3Cart first, second, normal;
				if (!ShapeOperator(u, w, a, b, c) || !UnitTangents(u, w, first, second, normal)) {
					dir1 = dir2 = Vec3Cart{};
					return;
				}
				const Vec3Cart perpendicular = GeometryDetail::Unit(VectorProduct(normal, first));
				const Real angle = std::atan2(b, a / 2 - c / 2) / 2;
				dir1 = std::cos(angle) * first + std::sin(angle) * perpendicular;
				dir2 = -std::sin(angle) * first + std::cos(angle) * perpendicular;
			}

			/// @brief Get first fundamental form coefficients (E,F,G) - measures surface metric
			/// @param u,w Surface parameters
			/// @param outE Output: E = ∂r/∂u · ∂r/∂u
			/// @param outF Output: F = ∂r/∂u · ∂r/∂w
			/// @param outG Output: G = ∂r/∂w · ∂r/∂w
			/// @note ds² = E du² + 2F du dw + G dw² (arc length formula)
			virtual void		 GetFirstNormalFormCoefficients(Real u, Real w, Real& outE, Real& outF, Real& outG) const
			{
				VectorN<Real, 3> r_by_u, r_by_w;
				Tangents(u, w, r_by_u, r_by_w);
				outE = GeometryDetail::Checked(Utils::ScalarProduct(r_by_u, r_by_u));
				outF = GeometryDetail::Checked(Utils::ScalarProduct(r_by_u, r_by_w));
				outG = GeometryDetail::Checked(Utils::ScalarProduct(r_by_w, r_by_w));
				if ((outE == 0 && GeometryDetail::Norm(r_by_u) != 0) ||
					(outG == 0 && GeometryDetail::Norm(r_by_w) != 0))
					throw GeometryError("Surface metric - squared tangent length underflow");
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
				const Vec3Cart n = GeometryDetail::Checked(Normal(u, w));
				if (GeometryDetail::Norm(n) == 0) {
					outL = outM = outN = 0.0;
					return;
				}

				Vec3Cart r_by_uu = GeometryDetail::Checked(Derivation::NDer2_uu(*this, u, w));
				Vec3Cart r_by_uw = GeometryDetail::Checked(Derivation::NDer2_uw(*this, u, w));
				Vec3Cart r_by_ww = GeometryDetail::Checked(Derivation::NDer2_ww(*this, u, w));

				outL = GeometryDetail::Checked(Utils::ScalarProduct(r_by_uu, n));
				outM = GeometryDetail::Checked(Utils::ScalarProduct(r_by_uw, n));
				outN = GeometryDetail::Checked(Utils::ScalarProduct(r_by_ww, n));
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
				Real a, b, c;
				return ShapeOperator(u, w, a, b, c) ? GeometryDetail::Checked(a * c - b * b) : 0;
			}

			/// @brief Compute mean curvature H = (EN + GL - 2FM)/(2(EG - F²)) = (κ₁+κ₂)/2
			/// @param u,w Surface parameters
			/// @return H (zero for minimal surfaces)
			/// @note Extrinsic invariant (depends on embedding in R³)
			Real MeanCurvature(Real u, Real w)
			{
				Real a, b, c;
				return ShapeOperator(u, w, a, b, c) ? GeometryDetail::Checked(a / 2 + c / 2) : 0;
			}

			bool isRegular(Real u, Real w) const
			{
				Vec3Cart first, second, normal;
				return UnitTangents(u, w, first, second, normal);
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
				GeometryDetail::RequireFinite(point, "PlaneSurface - point must be finite");
				GeometryDetail::RequireFinite(normal, "PlaneSurface - normal must be finite");
				GeometryDetail::RequireFinite(uAxis, "PlaneSurface - u axis must be finite");
				GeometryDetail::RequireFinite(vAxis, "PlaneSurface - v axis must be finite");
				_normal = GeometryDetail::Unit(normal);
				const Real uNorm = GeometryDetail::Norm(uAxis), vNorm = GeometryDetail::Norm(vAxis);
				if ((uNorm != 0 && std::abs(Utils::ScalarProduct(Vec3Cart(_normal), Vec3Cart(uAxis / uNorm))) > Precision::OrthogonalityTolerance) ||
					(vNorm != 0 && std::abs(Utils::ScalarProduct(Vec3Cart(_normal), Vec3Cart(vAxis / vNorm))) > Precision::OrthogonalityTolerance))
					throw GeometryError("uAxis and vAxis must be orthogonal to the normal.");
			}
			// this ctor gets only one local axis direction, the other is computed as cross product with normal
			PlaneSurface(const Vec3Cart& point, const Vec3Cart& normal, const Vec3Cart& uAxisDirection)
			{
				GeometryDetail::RequireFinite(point, "PlaneSurface - point must be finite");
				GeometryDetail::RequireFinite(normal, "PlaneSurface - normal must be finite");
				GeometryDetail::RequireFinite(uAxisDirection, "PlaneSurface - axis direction must be finite");
				_point = point;
				_normal = GeometryDetail::Unit(normal);
				const Vec3Cart direction = GeometryDetail::Unit(uAxisDirection);
				const Vec3Cart projected = direction - Utils::ScalarProduct(direction, _normal) * _normal;
				if (GeometryDetail::Norm(projected) <= Precision::RoundoffZeroThreshold)
					throw GeometryError("PlaneSurface - axis direction must not be parallel to the normal");
				_uAxis = GeometryDetail::Unit(projected);
				_vAxis = GeometryDetail::Unit(VectorProduct(Vec3Cart(_normal), Vec3Cart(_uAxis)));
			}
			// this ctor only gets point and normal, uAxis and vAxis are computed as orthogonal to normal
			// with uAxis trying to align as much as possible with x axis
			PlaneSurface(const Vec3Cart& point, const Vec3Cart& normal)
			{
				GeometryDetail::RequireFinite(point, "PlaneSurface - point must be finite");
				GeometryDetail::RequireFinite(normal, "PlaneSurface - normal must be finite");
				_point = point;
				_normal = GeometryDetail::Unit(normal);
				Vec3Cart cross = VectorProduct(_normal, Vec3Cart{1, 0, 0});
				if (GeometryDetail::Norm(cross) <= Precision::Vec3CartIsParallelTolerance)
					cross = VectorProduct(_normal, Vec3Cart{0, 1, 0});
				_uAxis = GeometryDetail::Unit(cross);
				_vAxis = GeometryDetail::Unit(VectorProduct(Vec3Cart(_normal), Vec3Cart(_uAxis)));
			}

			Real getMinU() const { return _minU; }
			Real getMaxU() const { return _maxU; }
			Real getMinW() const { return _minV; }
			Real getMaxW() const { return _maxV; }

			void SetPlaneLimits(Real minU, Real maxU, Real minV, Real maxV)
			{
				GeometryDetail::RequireRange(minU, maxU, "PlaneSurface - invalid u interval");
				GeometryDetail::RequireRange(minV, maxV, "PlaneSurface - invalid w interval");
				_minU = minU;
				_maxU = maxU;
				_minV = minV;
				_maxV = maxV;
			}

			VectorN<Real, 3> operator()(Real u, Real w) const
			{
				GeometryDetail::RequireFinite(u, "PlaneSurface - u must be finite");
				GeometryDetail::RequireFinite(w, "PlaneSurface - w must be finite");
				return GeometryDetail::Checked(_point + u * _uAxis + w * _vAxis);
			}
		};

		/// @brief Cylindrical surface r(u,w) = (R·cos(u), R·sin(u), w), radius R, height H
		/// @note u∈[0,2π], w∈[0,H]. Parabolic surface: K=0 (one principal curvature zero)
		class CylinderSurface : public ISurfaceCartesian
		{
			Real _radius, _height;
		public:
			CylinderSurface() : _radius(1), _height(1) {}
			CylinderSurface(Real R, Real H) : _radius(R), _height(H)
			{
				GeometryDetail::RequireNonnegative(R, "CylinderSurface - radius must be finite and nonnegative");
				GeometryDetail::RequireNonnegative(H, "CylinderSurface - height must be finite and nonnegative");
			}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return 2 * Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return _height; }

			VectorN<Real, 3> operator()(Real u, Real w) const
			{
				GeometryDetail::RequireFinite(u, "CylinderSurface - u must be finite");
				GeometryDetail::RequireFinite(w, "CylinderSurface - w must be finite");
				return GeometryDetail::Checked(Vec3Cart{_radius * cos(u), _radius * sin(u), w});
			}
		};
		
		/// @brief Torus r(u,w) = ((R+r·cos(w))·cos(u), (R+r·cos(w))·sin(u), r·sin(w))
		/// @note R=major radius, r=minor radius. u,w∈[0,2π]. K varies (elliptic inner, hyperbolic outer)
		class Torus : public ISurfaceCartesian
		{
			Real _majorRadius, _minorRadius;
		public:
			Torus() : _majorRadius(1), _minorRadius(0.5) {}
			Torus(Real R, Real r) : _majorRadius(R), _minorRadius(r)
			{
				GeometryDetail::RequireNonnegative(R, "Torus - major radius must be finite and nonnegative");
				GeometryDetail::RequireNonnegative(r, "Torus - minor radius must be finite and nonnegative");
			}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return 2 * Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real u, Real w) const
			{
				GeometryDetail::RequireFinite(u, "Torus - u must be finite");
				GeometryDetail::RequireFinite(w, "Torus - w must be finite");
				return GeometryDetail::Checked(Vec3Cart{(_majorRadius + _minorRadius * cos(w)) * cos(u),
					(_majorRadius + _minorRadius * cos(w)) * sin(u), _minorRadius * sin(w)});
			}
		};

		/// @brief Spherical surface r(u,w) = R·(sin(u)·cos(w), sin(u)·sin(w), cos(u)), radius R
		/// @note u∈[0,π] (latitude), w∈[0,2π] (longitude). Constant positive curvature K=1/R²
		class Sphere : public ISurfaceCartesian
		{
			Real _radius;
		public:
			Sphere() : _radius(1) {}
			Sphere(Real R) : _radius(R)
			{
				GeometryDetail::RequireNonnegative(R, "Sphere - radius must be finite and nonnegative");
			}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real u, Real w) const
			{
				GeometryDetail::RequireFinite(u, "Sphere - u must be finite");
				GeometryDetail::RequireFinite(w, "Sphere - w must be finite");
				return GeometryDetail::Checked(Vec3Cart{_radius * sin(u) * cos(w), _radius * sin(u) * sin(w), _radius * cos(u)});
			}
		};





		/// @brief Ellipsoid r(u,w) = (a·sin(u)·cos(w), b·sin(u)·sin(w), c·cos(u)), semi-axes (a,b,c)
		/// @note u∈[0,π], w∈[0,2π]. Generalization of sphere (a=b=c gives sphere)
		class Ellipsoid : public ISurfaceCartesian
		{
			Real _a, _b, _c;
		public:
			Ellipsoid() : _a(1), _b(1), _c(1) {}
			Ellipsoid(Real a, Real b, Real c) : _a(a), _b(b), _c(c)
			{
				GeometryDetail::RequireNonnegative(a, "Ellipsoid - a must be finite and nonnegative");
				GeometryDetail::RequireNonnegative(b, "Ellipsoid - b must be finite and nonnegative");
				GeometryDetail::RequireNonnegative(c, "Ellipsoid - c must be finite and nonnegative");
			}

		Real getMinU() const { return 0; }
		Real getMaxU() const { return Constants::PI; }
		Real getMinW() const { return 0; }
		Real getMaxW() const { return 2 * Constants::PI; }

		VectorN<Real, 3> operator()(Real u, Real w) const
		{
			GeometryDetail::RequireFinite(u, "Ellipsoid - u must be finite");
			GeometryDetail::RequireFinite(w, "Ellipsoid - w must be finite");
			return GeometryDetail::Checked(Vec3Cart{_a * sin(u) * cos(w), _b * sin(u) * sin(w), _c * cos(u)});
		}
		};

		/// @brief One-sheet hyperboloid r(u,w) = (a·cosh(u)·cos(w), b·cosh(u)·sin(w), c·sinh(u))
		/// @note Ruled surface with negative Gaussian curvature. u∈[0,π], w∈[0,2π]
		class Hyperboloid : public ISurfaceCartesian
		{
			Real _a, _b, _c;
		public:
			Hyperboloid() : _a(1), _b(1), _c(1) {}
			Hyperboloid(Real a, Real b, Real c) : _a(a), _b(b), _c(c)
			{
				GeometryDetail::RequireNonnegative(a, "Hyperboloid - a must be finite and nonnegative");
				GeometryDetail::RequireNonnegative(b, "Hyperboloid - b must be finite and nonnegative");
				GeometryDetail::RequireNonnegative(c, "Hyperboloid - c must be finite and nonnegative");
			}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real u, Real w) const
			{
				GeometryDetail::RequireFinite(u, "Hyperboloid - u must be finite");
				GeometryDetail::RequireFinite(w, "Hyperboloid - w must be finite");
				return GeometryDetail::Checked(Vec3Cart{_a * cosh(u) * cos(w), _b * cosh(u) * sin(w), _c * sinh(u)});
			}
		};

		/// @brief Paraboloid of revolution z = u, r(u,w) = (a·√(u/h)·cos(w), a·√(u/h)·sin(w), u)
		/// @note a=radius scaling, h=height scaling. Parabolic shape with positive curvature
		class Paraboloid : public ISurfaceCartesian
		{
			Real _a, _h;
		public:
			Paraboloid() : _a(1), _h(1) {}
			Paraboloid(Real a, Real h) : _a(a), _h(h)
			{
				GeometryDetail::RequireNonnegative(a, "Paraboloid - radius scale must be finite and nonnegative");
				if (!std::isfinite(h) || h <= 0)
					throw DomainError("Paraboloid - height scale must be finite and positive");
			}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return 10; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real u, Real w) const
			{
				GeometryDetail::RequireNonnegative(u, "Paraboloid - u must be finite and nonnegative");
				GeometryDetail::RequireFinite(w, "Paraboloid - w must be finite");
				return GeometryDetail::Checked(Vec3Cart{_a * sqrt(u / _h) * cos(w), _a * sqrt(u / _h) * sin(w), u});
			}
		};








	} // end namespace Surfaces
}

#endif // MML_CURVES_SURFACES_H