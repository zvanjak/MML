#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include <type_traits>
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/core/Curves.h>
#include <mml/core/Surfaces.h>
#endif

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Core::CurveSurfaceContracts
{
	class Parabola : public Curves::ICurveCartesian2D
	{
	public:
		Real getMinT() const override { return Constants::NegInf; }
		Real getMaxT() const override { return Constants::PosInf; }
		VectorN<Real, 2> operator()(Real t) const override { return {t, t * t}; }
	};

	class ExactHelix : public Curves::ICurveCartesian3D, public ICurveDerivativeProvider<3>
	{
		Real _scale, _rate;
	public:
		ExactHelix(Real scale, Real rate) : _scale(scale), _rate(rate) {}
		Real getMinT() const override { return Constants::NegInf; }
		Real getMaxT() const override { return Constants::PosInf; }
		VectorN<Real, 3> operator()(Real t) const override
		{
			return _scale * Vec3Cart{std::cos(_rate * t), std::sin(_rate * t), 2 * _rate * t};
		}
		VectorN<Real, 3> EvaluateDerivative(Real t, int order) const override
		{
			const Real angle = _rate * t;
			if (order == 1) return (_scale * _rate) * Vec3Cart{-std::sin(angle), std::cos(angle), 2};
			if (order == 2) return (_scale * _rate * _rate) * Vec3Cart{-std::cos(angle), -std::sin(angle), 0};
			return (_scale * _rate * _rate * _rate) * Vec3Cart{std::sin(angle), -std::cos(angle), 0};
		}
	};

	class ExactLine : public Curves::ICurveCartesian3D, public ICurveDerivativeProvider<3>
	{
		Real _speed;
	public:
		mutable int calls = 0;
		mutable Real lastT = 0;
		explicit ExactLine(Real speed) : _speed(speed) {}
		Real getMinT() const override { return Constants::NegInf; }
		Real getMaxT() const override { return Constants::PosInf; }
		VectorN<Real, 3> operator()(Real t) const override { return {_speed * t, 0, 0}; }
		VectorN<Real, 3> EvaluateDerivative(Real t, int order) const override
		{
			++calls;
			lastT = t;
			return order == 1 ? Vec3Cart{_speed, 0, 0} : Vec3Cart{};
		}
	};

	class CollapsedSurface : public Surfaces::ISurfaceCartesian
	{
	public:
		Real getMinU() const override { return -1; }
		Real getMaxU() const override { return 1; }
		Real getMinW() const override { return -1; }
		Real getMaxW() const override { return 1; }
		VectorN<Real, 3> operator()(Real u, Real) const override { return {u, 0, 0}; }
	};

	class QuadraticPatch : public Surfaces::ISurfaceCartesian
	{
		Real _scale, _a, _b, _shear;
	public:
		QuadraticPatch(Real scale, Real a, Real b, Real shear) : _scale(scale), _a(a), _b(b), _shear(shear) {}
		Real getMinU() const override { return -1; }
		Real getMaxU() const override { return 1; }
		Real getMinW() const override { return -1; }
		Real getMaxW() const override { return 1; }
		VectorN<Real, 3> operator()(Real u, Real w) const override
		{
			const Real x = _a * u + _shear * _b * w, y = _b * w;
			return _scale * Vec3Cart{x, y, (2 * x * x + 3 * y * y) / 2};
		}
		void Tangents(Real u, Real w, VectorN<Real, 3>& first, VectorN<Real, 3>& second) const override
		{
			const Real x = _a * u + _shear * _b * w, y = _b * w;
			first = _scale * Vec3Cart{_a, 0, 2 * _a * x};
			second = _scale * Vec3Cart{_shear * _b, _b, 2 * _shear * _b * x + 3 * _b * y};
		}
	};

	class PolarProvider : public Curves::ICurvePolar2D, public ICurveDerivativeProvider<2>
	{
	public:
		mutable int calls = 0;
		Real getMinT() const override { return Constants::NegInf; }
		Real getMaxT() const override { return Constants::PosInf; }
		VectorN<Real, 2> operator()(Real t) const override { return {2, std::atan2(std::sin(t), std::cos(t))}; }
		VectorN<Real, 2> EvaluateDerivative(Real, int order) const override
		{
			++calls;
			return order == 1 ? Vec2Cart{0, 1} : Vec2Cart{};
		}
	};

	TEST_CASE("Plane constructors retain rank orientation and explicit axis scales", "[geometry-contract][surfaces]")
	{
		const Vec3Cart origin{1, 2, 3}, normal{0, 0, 7};
		Surfaces::PlaneSurface projected(origin, normal, Vec3Cart{2, 0, 3});
		const Vec3Cart u = projected(1, 0) - projected(0, 0);
		const Vec3Cart v = projected(0, 1) - projected(0, 0);
		CHECK(std::abs(u.NormL2() - 1) < Precision::RoundoffZeroThreshold);
		CHECK(std::abs(v.NormL2() - 1) < Precision::RoundoffZeroThreshold);
		CHECK(std::abs(Utils::ScalarProduct(u, v)) < Precision::RoundoffZeroThreshold);
		CHECK(u[0] == Real{1});
		CHECK(v[1] == Real{1});
		CHECK(projected.isRegular(0, 0));
		Surfaces::PlaneSurface automatic(origin, normal);
		CHECK((automatic(1, 0) - origin).NormL2() == Real{1});
		Surfaces::PlaneSurface explicitAxes(origin, normal, Vec3Cart{2, 0, 0}, Vec3Cart{0, 3, 0}, -2, 4, -3, 5);
		CHECK(explicitAxes(2, 3)[0] == Real{5});
		CHECK(explicitAxes(2, 3)[1] == Real{11});
		CHECK(explicitAxes.getMinU() == Real{-2});
		CHECK(explicitAxes.getMaxW() == Real{5});
	}

	TEST_CASE("Log spiral uses its finite signed scale", "[geometry-contract][curves]")
	{
		const Real c = Real(GENERATE(-3, 4));
		Curves::LogSpiralCurve spiral(-1, c);
		CHECK(spiral(0)[0] == c);
		const Real t = Real(0.7);
		CHECK(std::abs(spiral(t)[0] - c * std::exp(-t) * std::cos(t)) <= Precision::RoundoffZeroThreshold * std::abs(c));
		CHECK(std::abs(spiral(t)[1] - c * std::exp(-t) * std::sin(t)) <= Precision::RoundoffZeroThreshold * std::abs(c));
	}

	TEST_CASE("Cartesian planar normal is perpendicular acceleration", "[geometry-contract][curves]")
	{
		Parabola curve;
		const auto tangent = curve.getTangent(1);
		const auto normal = curve.getNormal(1);
		CHECK(std::abs(Utils::ScalarProduct(tangent, normal)) < TOL(1e-7, 1e-3));
		CHECK(std::abs(normal[0] + Real(0.8)) < TOL(1e-7, 1e-3));
		CHECK(std::abs(normal[1] - Real(0.4)) < TOL(1e-7, 1e-3));
		CHECK(std::abs(curve.getNormalUnit(1).NormL2() - 1) < Precision::RoundoffZeroThreshold);
	}

	TEST_CASE("Numerical circle curvature and frames survive unit changes", "[geometry-contract][curves]")
	{
		const Real radius = GENERATE(Real(1e-8), Real(1e-5), Real{1}, Real(1e5), Real(1e16));
		CAPTURE(radius);
		Curves::Circle2DCurve planar(radius);
		Curves::Circle3DXYCurve spatial(radius);
		CHECK(std::abs(planar.getCurvature(Real(0.5)) * radius - 1) < TOL(1e-6, 1e-3));
		CHECK(std::abs(planar.getSignedCurvature(Real(0.5)) * radius - 1) < TOL(1e-6, 1e-3));
		CHECK(std::abs(planar.getRadiusOfCurvature(Real(0.5)) / radius - 1) < TOL(1e-6, 1e-3));
		CHECK(std::abs(spatial.getCurvature(Real(0.5)) * radius - 1) < TOL(1e-6, 1e-3));
		CHECK(std::abs(spatial.getNormal(Real(0.5)).NormL2() - 1) < Precision::RoundoffZeroThreshold);
		CHECK(std::abs(spatial.getBinormal(Real(0.5)).NormL2() - 1) < TOL(1e-6, 1e-3));
	}

	TEST_CASE("Analytic Frenet frames respect scaling and reparameterization", "[geometry-contract][curves]")
	{
		const Real scale = GENERATE(TOL(1e-100, 1e-12), Real(1e-8), Real{1}, Real(1e8), TOL(1e100, 1e12));
		const Real rate = GENERATE(Real(1e-8), Real{1}, Real(1e8));
		CAPTURE(scale, rate);
		ExactHelix curve(scale, rate);
		const auto t = curve.getTangentUnit(Real(0.5));
		const auto n = curve.getNormal(Real(0.5));
		const auto b = curve.getBinormal(Real(0.5));
		CHECK(std::abs(t.NormL2() - 1) < Precision::RoundoffZeroThreshold);
		CHECK(std::abs(n.NormL2() - 1) < Precision::RoundoffZeroThreshold);
		CHECK(std::abs(b.NormL2() - 1) < Precision::RoundoffZeroThreshold);
		CHECK(std::abs(Utils::ScalarProduct(t, n)) < Precision::RoundoffZeroThreshold);
		CHECK(std::abs(Utils::ScalarProduct(n, b)) < Precision::RoundoffZeroThreshold);
		CHECK(std::abs(curve.getCurvature(Real(0.5)) * scale - Real(1) / 5) < Precision::RoundoffZeroThreshold);
		CHECK(std::abs(curve.getTorsion(Real(0.5)) * scale - Real(2) / 5) < Precision::RoundoffZeroThreshold);
	}

	TEST_CASE("Polar Cartesian derivatives cross the angular seam and the polar origin", "[geometry-contract][curves]")
	{
		const Real center = Real(GENERATE(-1, 0, 2));
		const Real offset = GENERATE(Real(-0.001), Real{0}, Real(0.001));
		Curves::Circle2DCurvePolar circle(2, Pnt2Cart{center, 0});
		const Real t = Constants::PI + offset;
		const auto tangent = circle.getTangent(t);
		CHECK(std::abs(tangent[0] + 2 * std::sin(t)) < TOL(1e-6, 1e-3));
		CHECK(std::abs(tangent[1] - 2 * std::cos(t)) < TOL(1e-6, 1e-3));
		CHECK(std::abs(circle.getCurvature(t) * 2 - 1) < TOL(1e-5, 1e-2));
	}

	TEST_CASE("Surface normals and shape operator are invariant to uniform size", "[geometry-contract][surfaces]")
	{
		const Real radius = GENERATE(Real(1e-8), Real(1e-4), Real{1}, Real(1e4), Real(1e8));
		CAPTURE(radius);
		Surfaces::Sphere sphere(radius);
		const Real u = Constants::PI / 3, w = Real(0.7);
		CHECK(sphere.isRegular(u, w));
		CHECK(std::abs(sphere.Normal(u, w).NormL2() - 1) < Precision::RoundoffZeroThreshold);
		CHECK(std::abs(sphere.GaussianCurvature(u, w) * radius * radius - 1) < TOL(1e-5, 1e-2));
		CHECK(std::abs(sphere.MeanCurvature(u, w) * radius + 1) < TOL(1e-5, 1e-2));
		Real k1, k2;
		sphere.PrincipalCurvatures(u, w, k1, k2);
		CHECK(std::abs(k1 * radius + 1) < TOL(1e-5, 1e-2));
		CHECK(std::abs(k2 * radius + 1) < TOL(1e-5, 1e-2));
		Vec3Cart d1, d2;
		sphere.PrincipalDirections(u, w, d1, d2);
		CHECK(std::abs(d1.NormL2() - 1) < Precision::RoundoffZeroThreshold);
		CHECK(std::abs(d2.NormL2() - 1) < Precision::RoundoffZeroThreshold);
		CHECK(std::abs(Utils::ScalarProduct(d1, d2)) < TOL(1e-6, 1e-3));
	}

	TEST_CASE("Small and anisotropic plane charts retain regularity", "[geometry-contract][surfaces]")
	{
		const Real scaleU = GENERATE(Real(1e-8), Real{1}, Real(1e8));
		const Real scaleW = GENERATE(Real(1e-8), Real{1}, Real(1e8));
		Surfaces::PlaneSurface plane(Vec3Cart{}, Vec3Cart{0, 0, 1}, Vec3Cart{scaleU, 0, 0}, Vec3Cart{0, scaleW, 0});
		CHECK(plane.isRegular(0, 0));
		CHECK(std::abs(plane.Normal(0, 0).NormL2() - 1) < Precision::RoundoffZeroThreshold);
		CHECK(plane.GaussianCurvature(0, 0) == Real{0});
		CHECK(plane.MeanCurvature(0, 0) == Real{0});
	}

	TEST_CASE("Shape operator gives independent principal directions on anisotropic charts", "[geometry-contract][surfaces]")
	{
		const Real scale = GENERATE(Real(1e-6), Real{1}, Real(1e6));
		const Real a = GENERATE(Real(0.01), Real{1}, Real{100});
		const Real b = GENERATE(Real(0.01), Real{1}, Real{100});
		const Real shear = GENERATE(Real(-0.75), Real{0}, Real(0.75));
		QuadraticPatch patch(scale, a, b, shear);
		// Mixed-partial cancellation grows with the chart axis ratio.
		const Real derivativeBudget = Precision::RoundoffZeroThreshold * (a / b + b / a);
		Real k1, k2;
		patch.PrincipalCurvatures(0, 0, k1, k2);
		CHECK(std::abs(k1 * scale - 3) < derivativeBudget);
		CHECK(std::abs(k2 * scale - 2) < derivativeBudget);
		CHECK(std::abs(patch.GaussianCurvature(0, 0) * scale * scale - 6) < derivativeBudget);
		CHECK(std::abs(patch.MeanCurvature(0, 0) * scale - Real(2.5)) < derivativeBudget);
		Vec3Cart first, second;
		patch.PrincipalDirections(0, 0, first, second);
		CHECK(std::abs(std::abs(first[1]) - 1) < Precision::RoundoffZeroThreshold);
		CHECK(std::abs(std::abs(second[0]) - 1) < Precision::RoundoffZeroThreshold);
	}

	TEST_CASE("Polar curves retain analytic derivative dispatch at seams", "[geometry-contract][curves]")
	{
		PolarProvider curve;
		CHECK(std::abs(curve.getTangent(Constants::PI)[1] + 2) < Precision::RoundoffZeroThreshold);
		CHECK(std::abs(curve.getCurvature(Constants::PI) - Real(0.5)) < Precision::RoundoffZeroThreshold);
		CHECK(curve.calls == 4);
	}

	TEST_CASE("Degenerate surface compatibility remains explicit", "[geometry-contract][surfaces]")
	{
		CollapsedSurface surface;
		CHECK_FALSE(surface.isRegular(0, 0));
		CHECK(surface.Normal(0, 0).NormL2() == Real{0});
		const auto second = surface.GetSecondNormalFormCoefficients(0, 0);
		CHECK(second.NormL2() == Real{0});
		CHECK(surface.GaussianCurvature(0, 0) == Real{0});
		CHECK(surface.MeanCurvature(0, 0) == Real{0});
		Real k1, k2;
		surface.PrincipalCurvatures(0, 0, k1, k2);
		CHECK(k1 == Real{0});
		CHECK(k2 == Real{0});
		Vec3Cart d1, d2;
		surface.PrincipalDirections(0, 0, d1, d2);
		CHECK(d1.NormL2() == Real{0});
		CHECK(d2.NormL2() == Real{0});
	}

	TEST_CASE("Undefined curve frames and stationary curvature throw explicitly", "[geometry-contract][curves]")
	{
		ExactLine stationary(0), line(1);
		CHECK_THROWS_AS(stationary.getCurvature(0), GeometryError);
		CHECK_THROWS_AS(stationary.getTangentUnit(0), GeometryError);
		CHECK_THROWS_AS(line.getNormal(0), GeometryError);
		CHECK_THROWS_AS(line.getBinormal(0), GeometryError);
		CHECK_THROWS_AS(line.getTorsion(0), GeometryError);
		CHECK(line.getCurvature(0) == Real{0});
		CHECK(std::isinf(line.getRadiusOfCurvature(0)));
		CHECK(line.getCurvatureVector(0).NormL2() == Real{0});
		CHECK_THROWS_AS(Curves::Circle2DCurve(0).getCurvature(0), GeometryError);
		CHECK_THROWS_AS(Curves::HelixCurve(0, 1).getTorsion(0), GeometryError);
	}

	TEST_CASE("Arc length check samples speed even with one sample", "[geometry-contract][curves][arc-sampling]")
	{
		ExactLine unit(1), fast(2);
		CHECK_FALSE(fast.isArcLengthParametrized(0, 1, 1));
		CHECK(unit.isArcLengthParametrized(0, 1, 17));
		CHECK(unit.calls == 17);
		CHECK_THROWS_AS(unit.isArcLengthParametrized(0, 1, 0), DomainError);
		CHECK_THROWS_AS(unit.isArcLengthParametrized(0, 1, -1), DomainError);
		CHECK_THROWS_AS(unit.isArcLengthParametrized(1, 0, 2), DomainError);
		CHECK_THROWS_AS(unit.isArcLengthParametrized(0, Constants::PosInf, 2), DomainError);
	}

	TEST_CASE("Arc length check remains bounded at adjacent huge and opposite endpoints", "[geometry-contract][curves][arc-sampling]")
	{
		ExactLine line(1);
		CHECK(line.isArcLengthParametrized(-3, 9, 1));
		CHECK(line.lastT == Real{3});
		line.calls = 0;
		const Real huge = Real(1e20), adjacent = std::nextafter(huge, Constants::PosInf);
		CHECK(line.isArcLengthParametrized(huge, adjacent, 100));
		CHECK(line.calls == 100);
		line.calls = 0;
		CHECK(line.isArcLengthParametrized(-std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max(), 17));
		CHECK(line.calls == 17);
		CHECK(line.isArcLengthParametrized(2, 2, 1));
		CHECK(line.lastT == Real{2});
		CHECK_THROWS_AS(line.isArcLengthParametrized(std::numeric_limits<Real>::quiet_NaN(), 2, 1), DomainError);
	}

	TEST_CASE("Hyperboloid evaluation preserves the configured scalar precision", "[geometry-contract][hyperboloid-precision]")
	{
		static_assert(std::is_same_v<decltype(MML::cosh(Real{})), Real>);
		static_assert(std::is_same_v<decltype(MML::sinh(Real{})), Real>);
		static_assert(Precision::GeometryEpsilon == Precision::EqualityTolerance);
		const Surfaces::Hyperboloid surface(2, 3, 4);
		const Real u = GENERATE(Real(-1) / 2, Real(0), Real(1) / 2, Real(1) + Precision::Eps);
		const Real w = Real(1) / 3;
		const auto point = surface(u, w);
		CHECK(point[0] == 2 * std::cosh(u) * std::cos(w));
		CHECK(point[1] == 3 * std::cosh(u) * std::sin(w));
		CHECK(point[2] == 4 * std::sinh(u));
		const Real x = point[0] / 2, y = point[1] / 3, z = point[2] / 4;
		CHECK(IsNearlyEqual(x * x + y * y - z * z, Real(1), Precision::RoundoffZeroThreshold));
	}

	TEST_CASE("Invalid callbacks geometry data and plane frames surface errors", "[geometry-contract][domains]")
	{
		Curves::CurveCartesian3D bad(+[](Real) -> VectorN<Real, 3> {
			return {std::numeric_limits<Real>::quiet_NaN(), 0, 0};
		});
		CHECK_THROWS_AS(bad(0), GeometryError);
		CHECK_THROWS_AS(bad.getTangent(0), GeometryError);
		CHECK_THROWS_AS(bad(Constants::PosInf), DomainError);
		CHECK_THROWS_AS(Surfaces::PlaneSurface(Vec3Cart{}, Vec3Cart{}), GeometryError);
		CHECK_THROWS_AS(Surfaces::PlaneSurface(Vec3Cart{}, Vec3Cart{0, 0, 1}, Vec3Cart{0, 0, 2}), GeometryError);
		const Real small = TOL(1e-100, std::sqrt(std::numeric_limits<Real>::min()) / Real(4));
		REQUIRE(small > 0);
		Surfaces::PlaneSurface tinyDirections(Vec3Cart{}, Vec3Cart{0, 0, small}, Vec3Cart{small, 0, 0});
		CHECK(tinyDirections(1, 0)[0] == Real{1});
		CHECK(tinyDirections(0, 1)[1] == Real{1});
		CHECK(Surfaces::Sphere(0)(1, 1).NormL2() == Real{0});
		CHECK(Curves::Circle2DCurve(0)(1).NormL2() == Real{0});
		CHECK_THROWS_AS(Surfaces::Paraboloid{}(-1, 0), DomainError);
	}

	TEST_CASE("Primitive dimensions reject negative or nonfinite sizes", "[geometry-contract][domains]")
	{
		const Real bad = GENERATE(Real{-1}, Constants::PosInf, Constants::NegInf, std::numeric_limits<Real>::quiet_NaN());
		CHECK_THROWS_AS(Curves::Circle2DCurve(bad), DomainError);
		CHECK_THROWS_AS(Curves::Circle2DCurvePolar(bad), DomainError);
		CHECK_THROWS_AS(Curves::Circle3DXYCurve(bad), DomainError);
		CHECK_THROWS_AS(Curves::Circle3DXZCurve(bad), DomainError);
		CHECK_THROWS_AS(Curves::Circle3DYZCurve(bad), DomainError);
		CHECK_THROWS_AS(Curves::Circle3DCurve(bad, Vec3Cart{0, 0, 1}, Pnt3Cart{}), DomainError);
		CHECK_THROWS_AS(Curves::HelixCurve(bad, 1), DomainError);
		CHECK_THROWS_AS(Surfaces::Sphere(bad), DomainError);
		CHECK_THROWS_AS(Surfaces::CylinderSurface(bad, 1), DomainError);
		CHECK_THROWS_AS(Surfaces::CylinderSurface(1, bad), DomainError);
		CHECK_THROWS_AS(Surfaces::Torus(bad, 1), DomainError);
		CHECK_THROWS_AS(Surfaces::Torus(1, bad), DomainError);
		CHECK_THROWS_AS(Surfaces::Ellipsoid(bad, 1, 1), DomainError);
		CHECK_THROWS_AS(Surfaces::Ellipsoid(1, bad, 1), DomainError);
		CHECK_THROWS_AS(Surfaces::Ellipsoid(1, 1, bad), DomainError);
		CHECK_THROWS_AS(Surfaces::Hyperboloid(bad, 1, 1), DomainError);
		CHECK_THROWS_AS(Surfaces::Paraboloid(bad, 1), DomainError);
		CHECK_THROWS_AS(Surfaces::Paraboloid(1, bad), DomainError);
	}

	TEST_CASE("Callbacks bounds and divisors reject invalid construction", "[geometry-contract][domains]")
	{
		CHECK_THROWS_AS(Curves::CurveCartesian3D(nullptr), DomainError);
		CHECK_THROWS_AS(Surfaces::Paraboloid(1, 0), DomainError);
		CHECK_THROWS_AS(Curves::LogSpiralCurve(-1, 0), ArgumentError);
		CHECK_THROWS_AS(Curves::LogSpiralCurve(std::numeric_limits<Real>::quiet_NaN()), ArgumentError);
		CHECK_THROWS_AS(Curves::LogSpiralCurve(-1, Constants::PosInf), ArgumentError);
		CHECK_THROWS_AS(Curves::LineCurve(0, Pnt3Cart{}, 0, Pnt3Cart{1, 0, 0}), DomainError);
		auto callback = +[](Real t) -> VectorN<Real, 3> { return {t, 0, 0}; };
		CHECK_THROWS_AS(Curves::CurveCartesian3D(2, 1, callback), DomainError);
		Surfaces::PlaneSurface plane(Vec3Cart{}, Vec3Cart{0, 0, 1});
		CHECK_THROWS_AS(plane.SetPlaneLimits(1, -1, 0, 1), DomainError);
		CHECK(plane.getMinU() == Constants::NegInf);
	}
}
