///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Curves.h                                                            ///
///  Description: Parametric curve classes (2D, 3D, space curves)                     ///
///               Arc length, curvature, torsion, Frenet frame calculations           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_CURVES_H
#define MML_CURVES_H

#include "MMLBase.h"

#include "base/BaseUtils.h"
#include "base/Vector/VectorN.h"
#include "base/Vector/Vector.h"

#include "base/Function.h"
#include "mml/base/Geometry/Geometry3D.h"

#include "core/Derivation.h"
#include "core/Integration/PathIntegration.h"


namespace MML
{
	namespace Curves
	{
		/////////////////////////////              CARTESIAN PLANAR CURVES                  ///////////////////////////////
		/// @brief Abstract base class for Cartesian planar curves (2D)
		/// @details Provides tangent and normal vector calculations using numerical derivatives
		/// @brief Calculate curvature for any 2D parametric curve with Cartesian output
		/// @details Uses formula κ = |x'y'' - y'x''| / (x'² + y'²)^(3/2)
		/// @note This function assumes the curve returns Cartesian coordinates (x, y).
		///       For curves returning polar coordinates (r, θ), this formula is incorrect;
		///       use the polar curvature formula: κ = (r² + 2r'² - r·r'') / (r² + r'²)^(3/2)
		/// @param curve The 2D parametric curve (must return Cartesian x,y values)
		/// @param t Parameter value
		/// @return Scalar curvature at parameter t
		inline Real getCurvature2D(const IParametricCurve<2>& curve, Real t)
		{
			Vec2Cart r_prime = Derivation::DeriveCurve<2>(curve, t, nullptr);
			Vec2Cart r_double_prime = Derivation::DeriveCurveSec<2>(curve, t, nullptr);

			Real x_prime = r_prime[0];
			Real y_prime = r_prime[1];
			Real x_double_prime = r_double_prime[0];
			Real y_double_prime = r_double_prime[1];

			// κ = |x'y'' - y'x''| / (x'² + y'²)^(3/2)
			Real numerator = std::abs(x_prime * y_double_prime - y_prime * x_double_prime);
			Real speed_squared = x_prime * x_prime + y_prime * y_prime;
			Real denominator = std::pow(speed_squared, 1.5);

			if (denominator < 1e-15)
				return 0.0;  // Avoid division by zero at singular points

			return numerator / denominator;
		}

		/// @brief Calculate signed curvature for any 2D parametric curve with Cartesian output
		/// @details Uses formula κ = (x'y'' - y'x'') / (x'² + y'²)^(3/2)
		///          Positive = counter-clockwise, Negative = clockwise
		/// @note This function assumes the curve returns Cartesian coordinates (x, y).
		///       Not valid for polar-valued curves.
		/// @param curve The 2D parametric curve (must return Cartesian x,y values)
		/// @param t Parameter value
		/// @return Signed curvature at parameter t
		inline Real getSignedCurvature2D(const IParametricCurve<2>& curve, Real t)
		{
			Vec2Cart r_prime = Derivation::DeriveCurve<2>(curve, t, nullptr);
			Vec2Cart r_double_prime = Derivation::DeriveCurveSec<2>(curve, t, nullptr);

			Real x_prime = r_prime[0];
			Real y_prime = r_prime[1];
			Real x_double_prime = r_double_prime[0];
			Real y_double_prime = r_double_prime[1];

			// κ = (x'y'' - y'x'') / (x'² + y'²)^(3/2)  (signed version)
			Real numerator = x_prime * y_double_prime - y_prime * x_double_prime;
			Real speed_squared = x_prime * x_prime + y_prime * y_prime;
			Real denominator = std::pow(speed_squared, 1.5);

			if (denominator < 1e-15)
				return 0.0;  // Avoid division by zero at singular points

			return numerator / denominator;
		}

		/// @brief Calculate curvature for a polar curve r = r(θ) where θ is the parameter
		/// @details Uses formula κ = |r² + 2r'² - r·r''| / (r² + r'²)^(3/2)
		///          where r' = dr/dθ and r'' = d²r/dθ²
		/// @note For curves returning (r, θ) where θ = t (parameter equals angle)
		/// @param curve The polar curve returning (r, θ)
		/// @param t Parameter value (angle θ)
		/// @return Scalar curvature at parameter t
		inline Real getCurvaturePolar2D(const IParametricCurve<2>& curve, Real t)
		{
			// For polar curve, output is (r, θ)
			// We need derivatives of r with respect to θ
			// If θ = t (common case), then dr/dθ = dr/dt
			Vec2Cart vals = curve(t);
			Vec2Cart derivs = Derivation::DeriveCurve<2>(curve, t, nullptr);
			Vec2Cart derivs2 = Derivation::DeriveCurveSec<2>(curve, t, nullptr);

			Real r = vals[0];           // r value
			Real theta_prime = derivs[1]; // dθ/dt
			
			// dr/dθ = (dr/dt) / (dθ/dt)  using chain rule
			Real r_prime = (std::abs(theta_prime) > 1e-15) ? derivs[0] / theta_prime : 0.0;
			
			// d²r/dθ² requires more careful calculation
			// d²r/dθ² = d/dθ(dr/dθ) = (d/dt(dr/dθ)) / (dθ/dt)
			Real theta_double_prime = derivs2[1];
			Real r_double_prime_dt = derivs2[0];  // d²r/dt²
			
			// Using quotient rule: d/dt(r'/θ') = (r''·θ' - r'·θ'')/θ'²
			// Then divide by θ' to get d²r/dθ²
			Real r_double_prime = 0.0;
			if (std::abs(theta_prime) > 1e-15) {
				Real theta_prime_sq = theta_prime * theta_prime;
				r_double_prime = (r_double_prime_dt * theta_prime - derivs[0] * theta_double_prime) / (theta_prime_sq * theta_prime);
			}

			// κ = |r² + 2r'² - r·r''| / (r² + r'²)^(3/2)
			Real r_sq = r * r;
			Real r_prime_sq = r_prime * r_prime;
			Real numerator = std::abs(r_sq + 2.0 * r_prime_sq - r * r_double_prime);
			Real denominator = std::pow(r_sq + r_prime_sq, 1.5);

			if (denominator < 1e-15)
				return 0.0;

			return numerator / denominator;
		}

		/// @brief Abstract base class for Cartesian planar curves (2D)
		/// @details Provides tangent, normal, and curvature calculations using numerical derivatives
		class ICurveCartesian2D : public IParametricCurve<2>
		{
		public:
			Vec2Cart getTangent(Real t)
			{
				return Derivation::DeriveCurve<2>(*this, t, nullptr);
			}
			Vec2Cart getTangentUnit(Real t)
			{
				return getTangent(t).GetAsUnitVector();
			}
			Vec2Cart getNormal(Real t)
			{
				return Derivation::DeriveCurveSec<2>(*this, t, nullptr);
			}
			Vec2Cart getNormalUnit(Real t)
			{
				return getNormal(t).GetAsUnitVector();
			}

			/// @brief Get scalar curvature κ(t) = |x'y'' - y'x''| / (x'² + y'²)^(3/2)
			virtual Real getCurvature(Real t)
			{
				return getCurvature2D(*this, t);
			}

			/// @brief Get signed curvature (positive = CCW, negative = CW)
			virtual Real getSignedCurvature(Real t)
			{
				return getSignedCurvature2D(*this, t);
			}

			/// @brief Get radius of curvature R(t) = 1/κ(t)
			virtual Real getRadiusOfCurvature(Real t)
			{
				Real kappa = getCurvature(t);
				if (kappa < 1e-15)
					return Constants::PosInf;  // Straight line has infinite radius
				return 1.0 / kappa;
			}
		};

		/// @brief Circle curve in 2D Cartesian coordinates
		/// @details Parametrized as (x,y) = center + r(cos(t), sin(t)), t∈[0,2π]
		class Circle2DCurve : public ICurveCartesian2D
		{
			Real _radius;
			Pnt2Cart _center;
		public:
			Circle2DCurve() : _radius(1), _center(0, 0) {}
			Circle2DCurve(Real radius) : _radius(radius), _center(0,0) {}
			Circle2DCurve(Real radius, const Pnt2Cart& center) : _radius(radius), _center(center) {}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			VectorN<Real, 2> operator()(Real t) const { 
				return MML::VectorN<Real, 2>{_center.X()+_radius* cos(t), _center.Y()+_radius* sin(t)}; 
			}
		};

		/// @brief Logarithmic (exponential) spiral curve
		/// @details r(t) = c·e^(λt) with λ < 0, parametrized as (e^(λt)cos(t), e^(λt)sin(t))
		class LogSpiralCurve : public ICurveCartesian2D
		{
			Real _lambda, _c;
		public:
			LogSpiralCurve() : _lambda(-1), _c(1) {}
			LogSpiralCurve(Real lambda) : _lambda(lambda), _c(1) {
				if (lambda >= 0) throw std::invalid_argument("LogSpiralCurve: lambda must be negative.");
			}
			LogSpiralCurve(Real lambda, Real c) : _lambda(lambda), _c(c) {
				if (lambda >= 0) throw std::invalid_argument("LogSpiralCurve: lambda must be negative.");
				if (c == 0) throw std::invalid_argument("LogSpiralCurve: c must not be zero.");
			}

			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{exp(_lambda * t)* cos(t), exp(_lambda * t)* sin(t)}; }
		};

		/// @brief Lemniscate of Gerono curve (figure-eight curve)
		/// @details Parametrized as (cos(t)/(1+sin²(t)), sin(t)cos(t)/(1+sin²(t)))
		class LemniscateCurve : public ICurveCartesian2D
		{
		public:
			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{cos(t) / (1 + sin(t) * sin(t)), sin(t)* cos(t) / (1 + sin(t) * sin(t))}; }
		};

		/// @brief Deltoid (tricuspid hypocycloid) curve
		/// @details Three-cusped curve traced by point on circle rolling inside larger circle
		class DeltoidCurve : public ICurveCartesian2D
		{
			int _n;
		public:
			DeltoidCurve() : _n(1) {}
			DeltoidCurve(int n) : _n(n) {}

			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{2 * _n * cos(t) * (1 + cos(t)), 2 * _n * sin(t) * (1 - cos(t))}; }
		};

		/// @brief Astroid (four-cusped hypocycloid) curve
		/// @details Star-shaped curve with parametrization (c·cos³(t), c·sin³(t)), c > 0
		class AstroidCurve : public ICurveCartesian2D
		{
			Real _c;
		public:
			AstroidCurve() : _c(1) {}
			AstroidCurve(Real c) : _c(c) {
				if (c <= 0) throw std::invalid_argument("AstroidCurve: c must be positive.");
			}

			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{_c* cos(t)* cos(t)* cos(t), _c* sin(t)* sin(t)* sin(t)}; }
		};

		/// @brief Epitrochoid curve (generalized epicycloid)
		/// @details Curve traced by point attached to circle rolling outside another circle
		class EpitrochoidCurve : public ICurveCartesian2D
		{
			Real _radius, _c;
			int _n;
		public:
			EpitrochoidCurve() : _radius(1), _c(1), _n(1) {}
			EpitrochoidCurve(Real radius, Real c, int n) : _radius(radius), _c(c), _n(n) {}

			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{cos(t) - _c * cos(_n * t), sin(t) - _c * sin(_n * t) }; }
		};

		/// @brief Archimedean spiral curve
		/// @details Linear spiral with constant separation between turns, r = a·t
		class ArchimedeanSpiralCurve : public ICurveCartesian2D
		{
			Real _a;
		public:
			ArchimedeanSpiralCurve() : _a(1) {}
			ArchimedeanSpiralCurve(Real a) : _a(a) {}

			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{_a* t* cos(t), _a* t* sin(t)}; }
		};

		/// @brief Butterfly curve (Temple H. Fay)
		/// @details A beautiful transcendental plane curve discovered by Temple H. Fay.
		///          Parametric equations:
		///            x(t) = sin(t) · (e^cos(t) - 2·cos(4t) - sin⁵(t/12))
		///            y(t) = cos(t) · (e^cos(t) - 2·cos(4t) - sin⁵(t/12))
		///          Best visualized for t ∈ [0, 12π]
		class ButterflyCurve : public ICurveCartesian2D
		{
			Real _scale;
		public:
			ButterflyCurve() : _scale(1.0) {}
			ButterflyCurve(Real scale) : _scale(scale) {}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 12.0 * Constants::PI; }

			VectorN<Real, 2> operator()(Real t) const {
				Real sin_t12_5 = std::pow(std::sin(t / 12.0), 5);
				Real r = std::exp(std::cos(t)) - 2.0 * std::cos(4.0 * t) - sin_t12_5;
				return MML::VectorN<Real, 2>{ _scale * std::sin(t) * r, _scale * std::cos(t) * r };
			}
		};


		/////////////////////////////               POLAR PLANAR CURVES                  ///////////////////////////////
		/// @brief Abstract base class for polar planar curves (2D)
		/// @details Curves return (r, θ) values. Provides tangent in Cartesian coordinates
		///          and curvature using the polar formula.
		/// @note The operator() returns VectorN<Real,2>{r, θ} where θ is typically = t
		class ICurvePolar2D : public IParametricCurve<2>
		{
		public:
			/// @brief Convert polar point (r, θ) to Cartesian (x, y)
			Vec2Cart toCartesian(Real t) const
			{
				VectorN<Real, 2> polar = (*this)(t);
				Real r = polar[0];
				Real theta = polar[1];
				return Vec2Cart{r * std::cos(theta), r * std::sin(theta)};
			}

			/// @brief Get tangent vector in Cartesian coordinates
			/// @details Converts polar derivatives to Cartesian tangent:
			///          dx/dt = r'·cos(θ) - r·θ'·sin(θ)
			///          dy/dt = r'·sin(θ) + r·θ'·cos(θ)
			Vec2Cart getTangent(Real t)
			{
				VectorN<Real, 2> polar = (*this)(t);
				Vec2Cart derivs = Derivation::DeriveCurve<2>(*this, t, nullptr);
				
				Real r = polar[0];
				Real theta = polar[1];
				Real r_prime = derivs[0];      // dr/dt
				Real theta_prime = derivs[1];  // dθ/dt

				// Tangent in Cartesian coords
				Real dx_dt = r_prime * std::cos(theta) - r * theta_prime * std::sin(theta);
				Real dy_dt = r_prime * std::sin(theta) + r * theta_prime * std::cos(theta);
				return Vec2Cart{dx_dt, dy_dt};
			}

			Vec2Cart getTangentUnit(Real t)
			{
				return getTangent(t).GetAsUnitVector();
			}

			/// @brief Get normal vector (perpendicular to tangent, rotated 90° CCW)
			Vec2Cart getNormal(Real t)
			{
				Vec2Cart tang = getTangent(t);
				return Vec2Cart{-tang[1], tang[0]};  // Rotate 90° CCW
			}

			Vec2Cart getNormalUnit(Real t)
			{
				return getNormal(t).GetAsUnitVector();
			}

			/// @brief Get curvature using polar formula κ = |r² + 2r'² - r·r''| / (r² + r'²)^(3/2)
			virtual Real getCurvature(Real t)
			{
				return getCurvaturePolar2D(*this, t);
			}

			/// @brief Get radius of curvature R(t) = 1/κ(t)
			virtual Real getRadiusOfCurvature(Real t)
			{
				Real kappa = getCurvature(t);
				if (kappa < 1e-15)
					return Constants::PosInf;
				return 1.0 / kappa;
			}
		};

		/// @brief Circle curve in polar coordinates
		/// @details Returns (r, θ) with constant radius r, varying angle θ∈[0,2π]
		class Circle2DCurvePolar : public ICurvePolar2D
		{
			Real _radius;
			Pnt2Cart _center;

		public:
			Circle2DCurvePolar() : _radius(1), _center(0, 0) {}
			Circle2DCurvePolar(Real radius) : _radius(radius), _center(0, 0) {}
			Circle2DCurvePolar(Real radius, const Pnt2Cart& center) : _radius(radius), _center(center) {}
			
			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }
			
			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{_radius, t}; }
		};

		///////////////////////////////             CARTESIAN SPACE CURVES                  ///////////////////////////////
		/// @brief Abstract base class for 3D Cartesian space curves
		/// @details Provides Frenet-Serret frame (tangent, normal, binormal), curvature, torsion
		class ICurveCartesian3D : public IParametricCurve<3>
		{
		public:
			/// @brief Get tangent vector at parameter t
			/// @param t Parameter value
			/// @return Tangent vector r'(t)
			Vec3Cart getTangent(Real t) const
			{
				return Derivation::DeriveCurve<3>(*this, t, nullptr);
			}
			/// @brief Get unit tangent vector T(t) = r'(t)/|r'(t)|
			Vec3Cart getTangentUnit(Real t) const
			{
				return getTangent(t).GetAsUnitVector();
			}
			/// @brief Get second derivative (acceleration) vector r''(t)
			/// @details This is the acceleration, NOT the Frenet principal normal
			Vec3Cart getSecondDerivative(Real t) const
			{
				return Derivation::DeriveCurveSec<3>(*this, t, nullptr);
			}
			/// @brief Get Frenet principal normal unit vector N(t)
			/// @details N = (T'/|T'|) - perpendicular to tangent, points toward center of curvature
			/// @note Returns unit vector by definition of Frenet frame
			Vec3Cart getNormal(Real t) const
			{
				Vec3Cart r_prime = getTangent(t);
				Vec3Cart r_double_prime = getSecondDerivative(t);
				Vector3Cartesian cross1 = VectorProduct(Vector3Cartesian(r_double_prime), Vector3Cartesian(r_prime));
				Vector3Cartesian result = VectorProduct(Vector3Cartesian(r_prime), cross1);
				return Vec3Cart(result / (r_prime.NormL2() * cross1.NormL2()));
			}
			/// @brief Get principal normal unit vector N(t) (same as getNormal)
			/// @deprecated Use getNormal() instead
			Vec3Cart getNormalUnit(Real t) const
			{
				return getNormal(t);  // Already unit by definition
			}
			/// @brief Get binormal unit vector B(t) = T(t) × N(t)
			Vec3Cart getBinormal(Real t) const
			{
				Vector3Cartesian T = Vector3Cartesian(getTangentUnit(t));
				Vector3Cartesian N = Vector3Cartesian(getNormal(t));
				return Vec3Cart(VectorProduct(T, N));
			}

			/// @brief Get curvature vector κ(t)·N(t)
			Vec3Cart getCurvatureVector(Real t) const
			{
				return getCurvature(t) * getNormal(t);
			}
			
			/// @brief Get scalar curvature κ(t) = |r'×r''| / |r'|³
			virtual Real getCurvature(Real t) const
			{
				Vec3Cart r_prime = getTangent(t);
				Vec3Cart r_double_prime = getSecondDerivative(t);
				Vector3Cartesian cross = VectorProduct(Vector3Cartesian(r_prime), Vector3Cartesian(r_double_prime));
				Real speed_cubed = std::pow(r_prime.NormL2(), 3);
				return cross.NormL2() / speed_cubed;
			}
			
			/// @brief Get radius of curvature ρ(t) = 1/κ(t)
			virtual Real getRadiusOfCurvature(Real t) const
			{
				Real kappa = getCurvature(t);
				return (kappa > 0) ? 1.0 / kappa : std::numeric_limits<Real>::infinity();
			}
			/// @brief Get torsion τ(t) = (r'×r'')·r''' / |r'×r''|² (twist of curve)
			virtual Real getTorsion(Real t) const
			{
				Vec3Cart r_prime = getTangent(t);
				Vec3Cart r_double_prime = getSecondDerivative(t);
				Vec3Cart r_triple_prime = Derivation::DeriveCurveThird<3>(*this, t, nullptr);

				Vec3Cart cross_product = VectorProduct(r_prime, r_double_prime);
				Real numerator = Utils::ScalarProduct(cross_product, r_triple_prime);
				Real denominator = std::pow(cross_product.NormL2(), 2);

				return numerator / denominator;
			}

			/// @brief Get osculating plane at t (contains tangent and principal normal)
			/// @details The osculating plane is perpendicular to the binormal
			virtual Plane3D getOsculationPlane(Real t) const
			{
				return Plane3D(Vector3Cartesian((*this)(t)).getAsPoint(), Vector3Cartesian(getBinormal(t)));
			}
			/// @brief Get normal plane at t (perpendicular to tangent)
			virtual Plane3D getNormalPlane(Real t) const
			{
				return Plane3D(Vector3Cartesian((*this)(t)).getAsPoint(), Vector3Cartesian(getTangentUnit(t)));
			}
			/// @brief Get rectifying plane at t (contains tangent and binormal)
			virtual Plane3D getRectifyingPlane(Real t) const
			{
				return Plane3D(Vector3Cartesian((*this)(t)).getAsPoint(), Vector3Cartesian(getBinormal(t)));
			}

			/// @brief Get Frenet-Serret moving frame (T, N, B) at parameter t
			/// @param t Parameter value
			/// @param[out] tangent Unit tangent vector
			/// @param[out] normal Principal normal vector
			/// @param[out] binormal Binormal vector
			void getMovingTrihedron(Real t, Vector3Cartesian& tangent, Vector3Cartesian& normal, Vector3Cartesian& binormal)
			{
				tangent = Vector3Cartesian(getTangentUnit(t));
				normal = Vector3Cartesian(getNormalUnit(t));
				binormal = Vector3Cartesian(getBinormal(t));
			}

			/// @brief Check if curve is parametrized by arc length
			/// @param t1 Start parameter
			/// @param t2 End parameter
			/// @param numPoints Number of test points
			/// @return true if |r'(t)| ≈ 1 everywhere (natural parametrization)
			bool isArcLengthParametrized(Real t1, Real t2, int numPoints = 100) const
			{
				Real delta = (t2 - t1) / numPoints;
				for (Real t = t1 + delta; t < t2; t += delta)
				{
					Real len = PathIntegration::ParametricCurveLength(*this, t1, t);
					if (fabs(len - (t - t1)) > PrecisionValues<Real>::DefaultToleranceRelaxed)
						return false;
				}
				return true;
			}
		};

		/// @brief Generic 3D Cartesian curve from function pointer or lambda
		/// @details Allows creation of custom curves without subclassing
		class CurveCartesian3D : public ICurveCartesian3D
		{
			Real _minT;
			Real _maxT;
			VectorN<Real, 3>(*_func)(Real);
		public:
			CurveCartesian3D(VectorN<Real, 3>(*inFunc)(Real)) 
				: _func(inFunc), _minT(Constants::NegInf), _maxT(Constants::PosInf) {}
			CurveCartesian3D(Real minT, Real maxT, VectorN<Real, 3>(*inFunc)(Real)) 
				: _func(inFunc), _minT(minT), _maxT(maxT) {}

			Real getMinT() const { return _minT; }
			Real getMaxT() const { return _maxT; }

			virtual VectorN<Real, 3> operator()(Real x) const { return _func(x); }
		};

		/// @brief Straight line segment in 3D
		/// @details Parametrized as r(t) = p + t·v, with starting point and direction
		class LineCurve : public ICurveCartesian3D
		{
			Line3D  _line;
			Real _minT;
			Real _maxT;
		public:
			LineCurve(Real minT, Real maxT, const Point3Cartesian& pnt, const Vector3Cartesian& dir) : _line(pnt, dir), _minT(minT), _maxT(maxT) {}
			LineCurve(Real t1, const Point3Cartesian& pnt1, Real t2, const Point3Cartesian& pnt2)
			{
				// tocno samo ako je t1 = 0.0!!!
				_line.StartPoint() = pnt1;
				Vec3Cart dir = Vec3Cart(pnt1, pnt2);
				_line.Direction() = dir.NormL2() / (t2 - t1) * dir.GetAsUnitVector();
				_minT = t1;
				_maxT = t2;
			}

			Real getMinT() const { return _minT; }
			Real getMaxT() const { return _maxT; }

			VectorN<Real, 3> operator()(Real t) const
			{
				//if (t < _minT || t > _maxT)
				//	throw std::invalid_argument("LineCurve: t is out of range.");

				auto pnt = _line(t);
				return VectorN<Real, 3>{pnt.X(), pnt.Y(), pnt.Z()};
			}
		};

		/// @brief Circle in 3D XY-plane (z=const)
		/// @details Parametrized as (center.x+r·cos(t), center.y+r·sin(t), center.z)
		class Circle3DXYCurve : public ICurveCartesian3D {
			Real _radius;
			Pnt3Cart _center;
		public:
			Circle3DXYCurve() : _radius(1), _center(0, 0, 0) {}
			Circle3DXYCurve(Real radius) : _radius(radius), _center(0, 0, 0) {}
			Circle3DXYCurve(Real radius, const Pnt3Cart& center) : _radius(radius), _center(center) {}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real t) const { 
				return MML::VectorN<Real, 3>{_center.X() + _radius * cos(t), _center.Y() + _radius * sin(t), _center.Z()}; 
			}
		};

		/// @brief Circle in 3D XZ-plane (y=const)
		/// @details Parametrized as (center.x+r·cos(t), center.y, center.z+r·sin(t))
		class Circle3DXZCurve : public ICurveCartesian3D {
			Real _radius;
			Pnt3Cart _center;
		public:
			Circle3DXZCurve() : _radius(1), _center(0, 0, 0) {}
			Circle3DXZCurve(Real radius) : _radius(radius), _center(0, 0, 0) {}
			Circle3DXZCurve(Real radius, const Pnt3Cart& center) : _radius(radius), _center(center) {}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real t) const { 
				return MML::VectorN<Real, 3>{_center.X() + _radius * cos(t), _center.Y(), _center.Z() + _radius * sin(t)}; 
			}
		};

		/// @brief Circle in 3D YZ-plane (x=const)
		/// @details Parametrized as (center.x, center.y+r·cos(t), center.z+r·sin(t))
		class Circle3DYZCurve : public ICurveCartesian3D {
			Real _radius;
			Pnt3Cart _center;
		public:
			Circle3DYZCurve() : _radius(1), _center(0, 0, 0) {}
			Circle3DYZCurve(Real radius) : _radius(radius), _center(0, 0, 0) {}
			Circle3DYZCurve(Real radius, const Pnt3Cart& center) : _radius(radius), _center(center) {}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real t) const { 
				return MML::VectorN<Real, 3>{_center.X(), _center.Y() + _radius * cos(t), _center.Z() + _radius * sin(t)}; 
			}
		};

		/// @brief General circle in 3D lying in arbitrary plane
		/// @details Circle with given radius, center, and plane normal vector
		class Circle3DCurve : public ICurveCartesian3D {
			Real _radius;
			Vec3Cart _normal;
			Pnt3Cart _center;
		public:
			Circle3DCurve(Real radius, const Vec3Cart& normal, const Pnt3Cart& center) 
				: _radius(radius), _normal(normal), _center(center) 
			{
				if (normal.NormL2() == 0.0) throw std::invalid_argument("Circle3D: normal must not be zero vector.");
				_normal = normal.GetAsUnitVector();
			}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real t) const
			{
				Vec3Cart u, v;

				// Find two orthonormal (unit) vectors (u, v) perpendicular to _normal
				_normal.GetPerpendicularVectors(u, v);

				// Parametric equation for the circle in 3D
				Vec3Cart point = _center + _radius * (u * std::cos(t) + v * std::sin(t));

				return VectorN<Real, 3>{point.X(), point.Y(), point.Z()};
			}
		};

		/// @brief Circular helix curve
		/// @details Parametrized as (r·cos(t), r·sin(t), b·t) with constant curvature and torsion
		class HelixCurve : public ICurveCartesian3D
		{
			Real _radius, _b;
		public:
			HelixCurve() : _radius(1.0), _b(1.0) {}
			HelixCurve(Real radius, Real b) : _radius(radius), _b(b) {}

			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 3> operator()(Real t) const { 
				return MML::VectorN<Real, 3>{_radius * cos(t), _radius * sin(t), _b * t}; 
			}

			/// @brief Get constant curvature κ = r/(r²+b²)
			Real getCurvature(Real t) const { return _radius / (POW2(_radius) + POW2(_b)); }
			/// @brief Get constant torsion τ = b/(r²+b²)
			Real getTorsion(Real t) const { return _b / (POW2(_radius) + POW2(_b)); }
		};

		/// @brief Twisted cubic space curve
		/// @details Parametrized as (t, t², t³), simplest non-planar algebraic curve
		class TwistedCubicCurve : public ICurveCartesian3D
		{
		public:
			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 3> operator()(Real t) const { return MML::VectorN<Real, 3>{t, t* t, t* t* t}; }
		};

		/// @brief Toroidal spiral curve winding around a torus
		/// @details Parametrized curve on torus surface with n windings
		class ToroidalSpiralCurve : public ICurveCartesian3D
		{
			int _n;
			Real _scale = 1.0;

		public:
			ToroidalSpiralCurve() : _n(1) {}
			ToroidalSpiralCurve(int n) : _n(n) {}
			ToroidalSpiralCurve(Real scale) : _n(1), _scale(scale) {}
			ToroidalSpiralCurve(int n, Real scale) : _n(n), _scale(scale) {}

			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 3> operator()(Real t) const { return MML::VectorN<Real, 3>{(_scale* (4 + sin(_n * t))* cos(t)), _scale* (4 + sin(_n * t))* sin(t), _scale* cos(_n* t)}; }
		};
	}
}

#endif // MML_CURVES_H