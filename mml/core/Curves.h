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

#include <mml/MMLBase.h>

#include <mml/base/BaseUtils.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Vector/VectorTypes2D.h>
#include <mml/base/Vector/VectorTypes3D.h>

#include <mml/base/Function.h>
#include <mml/base/Geometry/Geometry3D.h>

#include <mml/core/Derivation.h>
#include <mml/core/Derivation/DerivativeProviders.h>
#include <mml/core/Integration/PathIntegration.h>
#include <mml/core/detail/GeometryChecks.h>


namespace MML
{
	namespace Curves
	{
		namespace Detail
		{
			inline Real SignedCurvature(const Vec2Cart& velocity, const Vec2Cart& acceleration)
			{
				const Real speed = GeometryDetail::Norm(velocity);
				if (speed == 0)
					throw GeometryError("Curve curvature - zero speed");
				const Real accelerationNorm = GeometryDetail::Norm(acceleration);
				if (accelerationNorm == 0)
					return 0;
				const Vec2Cart tangent = velocity / speed, direction = acceleration / accelerationNorm;
				return GeometryDetail::Checked(
					((tangent[0] * direction[1] - tangent[1] * direction[0]) * (accelerationNorm / speed)) / speed);
			}

			class PolarCartesianCurve : public IParametricCurve<2>
			{
				const IParametricCurve<2>& _curve;
			public:
				explicit PolarCartesianCurve(const IParametricCurve<2>& curve) : _curve(curve) {}
				Real getMinT() const override { return _curve.getMinT(); }
				Real getMaxT() const override { return _curve.getMaxT(); }
				VectorN<Real, 2> operator()(Real t) const override
				{
					const auto point = GeometryDetail::Checked(_curve(t));
					return GeometryDetail::Checked(Vec2Cart{point[0] * std::cos(point[1]), point[0] * std::sin(point[1])});
				}
			};

			inline Vec2Cart PolarCartesianDerivative(const IParametricCurve<2>& curve, Real t, int order)
			{
				GeometryDetail::RequireFinite(t, "Polar curve - parameter must be finite");
				if (const auto* provider = dynamic_cast<const ICurveDerivativeProvider<2>*>(&curve))
				{
					const auto point = GeometryDetail::Checked(curve(t));
					const auto first = GeometryDetail::Checked(provider->EvaluateDerivative(t, 1));
					const Real c = std::cos(point[1]), s = std::sin(point[1]);
					if (order == 1)
						return GeometryDetail::Checked(Vec2Cart{
							first[0] * c - point[0] * first[1] * s, first[0] * s + point[0] * first[1] * c});
					const auto second = GeometryDetail::Checked(provider->EvaluateDerivative(t, 2));
					const Real radial = second[0] - point[0] * first[1] * first[1];
					const Real angular = 2 * first[0] * first[1] + point[0] * second[1];
					return GeometryDetail::Checked(Vec2Cart{radial * c - angular * s, radial * s + angular * c});
				}
				const PolarCartesianCurve cartesian(curve);
				return GeometryDetail::Checked(Differential::CurveDerivative<2>(cartesian, t, order));
			}
		}
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
			GeometryDetail::RequireFinite(t, "Curve curvature - parameter must be finite");
			return std::abs(Detail::SignedCurvature(Differential::CurveDerivative<2>(curve, t, 1),
				Differential::CurveDerivative<2>(curve, t, 2)));
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
			GeometryDetail::RequireFinite(t, "Curve curvature - parameter must be finite");
			return Detail::SignedCurvature(Differential::CurveDerivative<2>(curve, t, 1),
				Differential::CurveDerivative<2>(curve, t, 2));
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
			return std::abs(Detail::SignedCurvature(Detail::PolarCartesianDerivative(curve, t, 1),
				Detail::PolarCartesianDerivative(curve, t, 2)));
		}

		/// @brief Abstract base class for Cartesian planar curves (2D)
		/// @details Provides tangent, normal, and curvature calculations using numerical derivatives
		class ICurveCartesian2D : public IParametricCurve<2>
		{
		public:
			Vec2Cart getTangent(Real t)
			{
				GeometryDetail::RequireFinite(t, "Curve tangent - parameter must be finite");
				return GeometryDetail::Checked(Differential::CurveDerivative<2>(*this, t, 1));
			}
			Vec2Cart getTangentUnit(Real t)
			{
				return GeometryDetail::Unit(getTangent(t));
			}
			Vec2Cart getSecondDerivative(Real t)
			{
				GeometryDetail::RequireFinite(t, "Curve acceleration - parameter must be finite");
				return GeometryDetail::Checked(Differential::CurveDerivative<2>(*this, t, 2));
			}
			Vec2Cart getNormal(Real t)
			{
				const auto tangent = getTangentUnit(t);
				const auto acceleration = getSecondDerivative(t);
				const Real perpendicular = tangent[0] * acceleration[1] - tangent[1] * acceleration[0];
				return GeometryDetail::Checked(Vec2Cart{-tangent[1] * perpendicular, tangent[0] * perpendicular});
			}
			Vec2Cart getNormalUnit(Real t)
			{
				return GeometryDetail::Unit(getNormal(t));
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
				if (kappa == 0)
					return Constants::PosInf;  // Straight line has infinite radius
				return MML::Real(1.0) / kappa;
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
			Circle2DCurve(Real radius) : Circle2DCurve(radius, Pnt2Cart{}) {}
			Circle2DCurve(Real radius, const Pnt2Cart& center) : _radius(radius), _center(center)
			{
				GeometryDetail::RequireNonnegative(radius, "Circle2DCurve - radius must be finite and nonnegative");
				GeometryDetail::RequireFinite(Vec2Cart{center.X(), center.Y()}, "Circle2DCurve - center must be finite");
			}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			VectorN<Real, 2> operator()(Real t) const { 
				GeometryDetail::RequireFinite(t, "Circle2DCurve - parameter must be finite");
				return GeometryDetail::Checked(Vec2Cart{_center.X() + _radius * std::cos(t), _center.Y() + _radius * std::sin(t)});
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
				if (!std::isfinite(lambda) || lambda >= 0) throw ArgumentError("LogSpiralCurve: lambda must be finite and negative.");
			}
			LogSpiralCurve(Real lambda, Real c) : _lambda(lambda), _c(c) {
				if (!std::isfinite(lambda) || lambda >= 0) throw ArgumentError("LogSpiralCurve: lambda must be finite and negative.");
				if (!std::isfinite(c) || c == 0) throw ArgumentError("LogSpiralCurve: c must be finite and nonzero.");
			}

			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 2> operator()(Real t) const
			{
				GeometryDetail::RequireFinite(t, "LogSpiralCurve - parameter must be finite");
				const Real radius = GeometryDetail::Checked(_c * std::exp(_lambda * t));
				return GeometryDetail::Checked(Vec2Cart{radius * std::cos(t), radius * std::sin(t)});
			}
		};









		/// @brief Archimedean spiral curve
		/// @details Linear spiral with constant separation between turns, r = a·t
		class ArchimedeanSpiralCurve : public ICurveCartesian2D
		{
			Real _a;
		public:
			ArchimedeanSpiralCurve() : _a(1) {}
			ArchimedeanSpiralCurve(Real a) : _a(a)
			{
				GeometryDetail::RequireFinite(a, "ArchimedeanSpiralCurve - scale must be finite");
			}

			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 2> operator()(Real t) const
			{
				GeometryDetail::RequireFinite(t, "ArchimedeanSpiralCurve - parameter must be finite");
				return GeometryDetail::Checked(Vec2Cart{_a * t * std::cos(t), _a * t * std::sin(t)});
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
				GeometryDetail::RequireFinite(t, "Polar curve - parameter must be finite");
				VectorN<Real, 2> polar = GeometryDetail::Checked((*this)(t));
				Real r = polar[0];
				Real theta = polar[1];
				return GeometryDetail::Checked(Vec2Cart{r * std::cos(theta), r * std::sin(theta)});
			}

			/// @brief Get tangent vector in Cartesian coordinates
			/// @details Converts polar derivatives to Cartesian tangent:
			///          dx/dt = r'·cos(θ) - r·θ'·sin(θ)
			///          dy/dt = r'·sin(θ) + r·θ'·cos(θ)
			Vec2Cart getTangent(Real t)
			{
				return Detail::PolarCartesianDerivative(*this, t, 1);
			}

			Vec2Cart getTangentUnit(Real t)
			{
				return GeometryDetail::Unit(getTangent(t));
			}

			/// @brief Get normal vector (perpendicular to tangent, rotated 90° CCW)
			Vec2Cart getNormal(Real t)
			{
				Vec2Cart tang = getTangent(t);
				return Vec2Cart{-tang[1], tang[0]};  // Rotate 90° CCW
			}

			Vec2Cart getNormalUnit(Real t)
			{
				return GeometryDetail::Unit(getNormal(t));
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
				if (kappa == 0)
					return Constants::PosInf;
				return MML::Real(1.0) / kappa;
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
			Circle2DCurvePolar(Real radius) : Circle2DCurvePolar(radius, Pnt2Cart{}) {}
			Circle2DCurvePolar(Real radius, const Pnt2Cart& center) : _radius(radius), _center(center)
			{
				GeometryDetail::RequireNonnegative(radius, "Circle2DCurvePolar - radius must be finite and nonnegative");
				GeometryDetail::RequireFinite(Vec2Cart{center.X(), center.Y()}, "Circle2DCurvePolar - center must be finite");
			}
			
			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }
			
			VectorN<Real, 2> operator()(Real t) const {
				GeometryDetail::RequireFinite(t, "Circle2DCurvePolar - parameter must be finite");
				Real x = _center.X() + _radius * cos(t);
				Real y = _center.Y() + _radius * sin(t);
				return GeometryDetail::Checked(Vec2Cart{std::hypot(x, y), std::atan2(y, x)});
			}
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
				GeometryDetail::RequireFinite(t, "Curve tangent - parameter must be finite");
				return GeometryDetail::Checked(Differential::CurveDerivative<3>(*this, t, 1));
			}
			/// @brief Get unit tangent vector T(t) = r'(t)/|r'(t)|
			Vec3Cart getTangentUnit(Real t) const
			{
				return GeometryDetail::Unit(getTangent(t));
			}
			/// @brief Get second derivative (acceleration) vector r''(t)
			/// @details This is the acceleration, NOT the Frenet principal normal
			Vec3Cart getSecondDerivative(Real t) const
			{
				GeometryDetail::RequireFinite(t, "Curve acceleration - parameter must be finite");
				return GeometryDetail::Checked(Differential::CurveDerivative<3>(*this, t, 2));
			}
			/// @brief Get Frenet principal normal unit vector N(t)
			/// @details N = (T'/|T'|) - perpendicular to tangent, points toward center of curvature
			/// @note Returns unit vector by definition of Frenet frame
			Vec3Cart getNormal(Real t) const
			{
				const Vec3Cart tangent = getTangentUnit(t);
				const Vec3Cart acceleration = GeometryDetail::Unit(getSecondDerivative(t));
				const Vec3Cart binormal = GeometryDetail::Unit(VectorProduct(tangent, acceleration));
				return GeometryDetail::Unit(VectorProduct(binormal, tangent));
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
				const Real curvature = getCurvature(t);
				return curvature == 0 ? Vec3Cart{} : Vec3Cart(GeometryDetail::Checked(curvature * getNormal(t)));
			}
			
			/// @brief Get scalar curvature κ(t) = |r'×r''| / |r'|³
			virtual Real getCurvature(Real t) const
			{
				const Vec3Cart velocity = getTangent(t), acceleration = getSecondDerivative(t);
				const Real speed = GeometryDetail::Norm(velocity);
				if (speed == 0)
					throw GeometryError("Curve curvature - zero speed");
				const Real accelerationNorm = GeometryDetail::Norm(acceleration);
				if (accelerationNorm == 0)
					return 0;
				const Vec3Cart cross = VectorProduct(Vec3Cart(velocity / speed), Vec3Cart(acceleration / accelerationNorm));
				return GeometryDetail::Checked((GeometryDetail::Norm(cross) * (accelerationNorm / speed)) / speed);
			}
			
			/// @brief Get radius of curvature ρ(t) = 1/κ(t)
			virtual Real getRadiusOfCurvature(Real t) const
			{
				Real kappa = getCurvature(t);
				return (kappa > 0) ? MML::Real(1.0) / kappa : std::numeric_limits<Real>::infinity();
			}
			/// @brief Get torsion τ(t) = (r'×r'')·r''' / |r'×r''|² (twist of curve)
			virtual Real getTorsion(Real t) const
			{
				const Vec3Cart velocity = getTangent(t), acceleration = getSecondDerivative(t);
				const Real speed = GeometryDetail::Norm(velocity), accelerationNorm = GeometryDetail::Norm(acceleration);
				const Vec3Cart cross = VectorProduct(Vec3Cart(GeometryDetail::Unit(velocity)),
					Vec3Cart(GeometryDetail::Unit(acceleration)));
				const Real crossNorm = GeometryDetail::Norm(cross);
				if (crossNorm == 0)
					throw GeometryError("Curve torsion - undefined Frenet frame");
				const Vec3Cart third = GeometryDetail::Checked(Differential::CurveDerivative<3>(*this, t, 3));
				const Real thirdNorm = GeometryDetail::Norm(third);
				if (thirdNorm == 0)
					return 0;
				const Real projection = Utils::ScalarProduct(Vec3Cart(cross / crossNorm), Vec3Cart(third / thirdNorm));
				return GeometryDetail::Checked(((projection / crossNorm) * (thirdNorm / accelerationNorm)) / speed);
			}

			/// @brief Get osculating plane at t (contains tangent and principal normal)
			/// @details The osculating plane is perpendicular to the binormal
			virtual Plane3D getOsculationPlane(Real t) const
			{
				return Plane3D(Vector3Cartesian((*this)(t)).ToPoint(), Vector3Cartesian(getBinormal(t)));
			}
			/// @brief Get normal plane at t (perpendicular to tangent)
			virtual Plane3D getNormalPlane(Real t) const
			{
				return Plane3D(Vector3Cartesian((*this)(t)).ToPoint(), Vector3Cartesian(getTangentUnit(t)));
			}
			/// @brief Get rectifying plane at t (contains tangent and binormal)
			/// @details The rectifying plane is perpendicular to the principal normal
			virtual Plane3D getRectifyingPlane(Real t) const
			{
				return Plane3D(Vector3Cartesian((*this)(t)).ToPoint(), Vector3Cartesian(getNormal(t)));
			}

			/// @brief Get Frenet-Serret moving frame (T, N, B) at parameter t
			/// @param t Parameter value
			/// @param[out] tangent Unit tangent vector
			/// @param[out] normal Principal normal vector
			/// @param[out] binormal Binormal vector
			void getMovingTrihedron(Real t, Vector3Cartesian& tangent, Vector3Cartesian& normal, Vector3Cartesian& binormal)
			{
				tangent = Vector3Cartesian(getTangentUnit(t));
				normal = Vector3Cartesian(getNormal(t));
				binormal = Vector3Cartesian(getBinormal(t));
			}

			/// @brief Check if curve is parametrized by arc length
			/// @param t1 Start parameter
			/// @param t2 End parameter
			/// @param numPoints Number of test points
			/// @return true if |r'(t)| ≈ 1 everywhere (natural parametrization)
			bool isArcLengthParametrized(Real t1, Real t2, int numPoints = 100) const
			{
				if (numPoints <= 0 || !std::isfinite(t1) || !std::isfinite(t2) || t1 > t2)
					throw DomainError("Arc-length parameter check - positive sample count and finite ordered endpoints required");
				for (int index = 0; index < numPoints; ++index)
				{
					const Real fraction = (static_cast<Real>(index) + 1) / (static_cast<Real>(numPoints) + 1);
					const Real t = std::signbit(t1) != std::signbit(t2) ?
						(1 - fraction) * t1 + fraction * t2 : t1 + fraction * (t2 - t1);
					const Real speed = GeometryDetail::Norm(getTangent(t));
					if (std::abs(speed - 1) > Precision::ArcLengthSpeedTolerance)
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
				: CurveCartesian3D(Constants::NegInf, Constants::PosInf, inFunc) {}
			CurveCartesian3D(Real minT, Real maxT, VectorN<Real, 3>(*inFunc)(Real)) 
				: _minT(minT), _maxT(maxT), _func(inFunc)
			{
				if (!inFunc)
					throw DomainError("CurveCartesian3D - callback must not be null");
				GeometryDetail::RequireRange(minT, maxT, "CurveCartesian3D - invalid parameter interval");
			}

			Real getMinT() const { return _minT; }
			Real getMaxT() const { return _maxT; }

			virtual VectorN<Real, 3> operator()(Real x) const
			{
				GeometryDetail::RequireFinite(x, "CurveCartesian3D - parameter must be finite");
				return GeometryDetail::Checked(_func(x));
			}
		};

		/// @brief Straight line segment in 3D
		/// @details Parametrized as r(t) = p + t·v, with starting point and direction
		class LineCurve : public ICurveCartesian3D
		{
			Line3D  _line;
			Real _minT;
			Real _maxT;
		public:
			LineCurve(Real minT, Real maxT, const Point3Cartesian& pnt, const Vector3Cartesian& dir) : _minT(minT), _maxT(maxT)
			{
				GeometryDetail::RequireRange(minT, maxT, "LineCurve - invalid parameter interval");
				GeometryDetail::RequireFinite(Vec3Cart{pnt.X(), pnt.Y(), pnt.Z()}, "LineCurve - point must be finite");
				GeometryDetail::RequireFinite(dir, "LineCurve - direction must be finite");
				_line.StartPoint() = pnt;
				_line.Direction() = GeometryDetail::Unit(dir);
			}
			LineCurve(Real t1, const Point3Cartesian& pnt1, Real t2, const Point3Cartesian& pnt2)
			{
				if (!std::isfinite(t1) || !std::isfinite(t2) || t1 >= t2)
					throw DomainError("LineCurve - finite strictly ordered endpoint parameters required");
				GeometryDetail::RequireFinite(Vec3Cart{pnt1.X(), pnt1.Y(), pnt1.Z()}, "LineCurve - first point must be finite");
				GeometryDetail::RequireFinite(Vec3Cart{pnt2.X(), pnt2.Y(), pnt2.Z()}, "LineCurve - second point must be finite");
				Vec3Cart dir = Vec3Cart(pnt1, pnt2);
				if (GeometryDetail::Norm(dir) == 0)
					throw GeometryError("LineCurve - distinct endpoints required");
				_line.Direction() = GeometryDetail::Checked(dir / GeometryDetail::Checked(t2 - t1));
				_line.StartPoint() = pnt1 + (-t1) * _line.Direction();
				GeometryDetail::Checked(Vec3Cart{_line.StartPoint().X(), _line.StartPoint().Y(), _line.StartPoint().Z()});
				_minT = t1;
				_maxT = t2;
			}

			Real getMinT() const { return _minT; }
			Real getMaxT() const { return _maxT; }

			VectorN<Real, 3> operator()(Real t) const
			{
				GeometryDetail::RequireFinite(t, "LineCurve - parameter must be finite");
				auto pnt = _line(t);
				return GeometryDetail::Checked(Vec3Cart{pnt.X(), pnt.Y(), pnt.Z()});
			}
		};

		/// @brief Circle in 3D XY-plane (z=const)
		/// @details Parametrized as (center.x+r·cos(t), center.y+r·sin(t), center.z)
		class Circle3DXYCurve : public ICurveCartesian3D {
			Real _radius;
			Pnt3Cart _center;
		public:
			Circle3DXYCurve() : _radius(1), _center(0, 0, 0) {}
			Circle3DXYCurve(Real radius) : Circle3DXYCurve(radius, Pnt3Cart{}) {}
			Circle3DXYCurve(Real radius, const Pnt3Cart& center) : _radius(radius), _center(center)
			{
				GeometryDetail::RequireNonnegative(radius, "Circle3DXYCurve - radius must be finite and nonnegative");
				GeometryDetail::RequireFinite(Vec3Cart{center.X(), center.Y(), center.Z()}, "Circle3DXYCurve - center must be finite");
			}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real t) const { 
				GeometryDetail::RequireFinite(t, "Circle3DXYCurve - parameter must be finite");
				return GeometryDetail::Checked(Vec3Cart{_center.X() + _radius * cos(t), _center.Y() + _radius * sin(t), _center.Z()});
			}
		};

		/// @brief Circle in 3D XZ-plane (y=const)
		/// @details Parametrized as (center.x+r·cos(t), center.y, center.z+r·sin(t))
		class Circle3DXZCurve : public ICurveCartesian3D {
			Real _radius;
			Pnt3Cart _center;
		public:
			Circle3DXZCurve() : _radius(1), _center(0, 0, 0) {}
			Circle3DXZCurve(Real radius) : Circle3DXZCurve(radius, Pnt3Cart{}) {}
			Circle3DXZCurve(Real radius, const Pnt3Cart& center) : _radius(radius), _center(center)
			{
				GeometryDetail::RequireNonnegative(radius, "Circle3DXZCurve - radius must be finite and nonnegative");
				GeometryDetail::RequireFinite(Vec3Cart{center.X(), center.Y(), center.Z()}, "Circle3DXZCurve - center must be finite");
			}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real t) const { 
				GeometryDetail::RequireFinite(t, "Circle3DXZCurve - parameter must be finite");
				return GeometryDetail::Checked(Vec3Cart{_center.X() + _radius * cos(t), _center.Y(), _center.Z() + _radius * sin(t)});
			}
		};

		/// @brief Circle in 3D YZ-plane (x=const)
		/// @details Parametrized as (center.x, center.y+r·cos(t), center.z+r·sin(t))
		class Circle3DYZCurve : public ICurveCartesian3D {
			Real _radius;
			Pnt3Cart _center;
		public:
			Circle3DYZCurve() : _radius(1), _center(0, 0, 0) {}
			Circle3DYZCurve(Real radius) : Circle3DYZCurve(radius, Pnt3Cart{}) {}
			Circle3DYZCurve(Real radius, const Pnt3Cart& center) : _radius(radius), _center(center)
			{
				GeometryDetail::RequireNonnegative(radius, "Circle3DYZCurve - radius must be finite and nonnegative");
				GeometryDetail::RequireFinite(Vec3Cart{center.X(), center.Y(), center.Z()}, "Circle3DYZCurve - center must be finite");
			}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real t) const { 
				GeometryDetail::RequireFinite(t, "Circle3DYZCurve - parameter must be finite");
				return GeometryDetail::Checked(Vec3Cart{_center.X(), _center.Y() + _radius * cos(t), _center.Z() + _radius * sin(t)});
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
				GeometryDetail::RequireNonnegative(radius, "Circle3DCurve - radius must be finite and nonnegative");
				GeometryDetail::RequireFinite(Vec3Cart{center.X(), center.Y(), center.Z()}, "Circle3DCurve - center must be finite");
				GeometryDetail::RequireFinite(normal, "Circle3DCurve - normal must be finite");
				_normal = GeometryDetail::Unit(normal);
			}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real t) const
			{
				GeometryDetail::RequireFinite(t, "Circle3DCurve - parameter must be finite");
				Vec3Cart u, v;

				// Find two orthonormal (unit) vectors (u, v) perpendicular to _normal
				_normal.GetPerpendicularVectors(u, v);

				// Parametric equation for the circle in 3D
				Vec3Cart point = _center + _radius * (u * std::cos(t) + v * std::sin(t));

				return GeometryDetail::Checked(point);
			}
		};

		/// @brief Circular helix curve
		/// @details Parametrized as (r·cos(t), r·sin(t), b·t) with constant curvature and torsion
		class HelixCurve : public ICurveCartesian3D
		{
			Real _radius, _b;
		public:
			HelixCurve() : _radius(1.0), _b(1.0) {}
			HelixCurve(Real radius, Real b) : _radius(radius), _b(b)
			{
				GeometryDetail::RequireNonnegative(radius, "HelixCurve - radius must be finite and nonnegative");
				GeometryDetail::RequireFinite(b, "HelixCurve - pitch must be finite");
			}

			Real getMinT() const { return Constants::NegInf; }
			Real getMaxT() const { return Constants::PosInf; }

			VectorN<Real, 3> operator()(Real t) const { 
				GeometryDetail::RequireFinite(t, "HelixCurve - parameter must be finite");
				return GeometryDetail::Checked(Vec3Cart{_radius * cos(t), _radius * sin(t), _b * t});
			}

			/// @brief Get constant curvature κ = r/(r²+b²)
			Real getCurvature(Real t) const
			{
				GeometryDetail::RequireFinite(t, "HelixCurve - parameter must be finite");
				const Real speed = GeometryDetail::Checked(std::hypot(_radius, _b));
				if (speed == 0) throw GeometryError("HelixCurve curvature - zero speed");
				return GeometryDetail::Checked((_radius / speed) / speed);
			}
			/// @brief Get constant torsion τ = b/(r²+b²)
			Real getTorsion(Real t) const
			{
				GeometryDetail::RequireFinite(t, "HelixCurve - parameter must be finite");
				if (_radius == 0) throw GeometryError("HelixCurve torsion - undefined Frenet frame");
				const Real speed = GeometryDetail::Checked(std::hypot(_radius, _b));
				return GeometryDetail::Checked((_b / speed) / speed);
			}
		};




	}
}

#endif // MML_CURVES_H