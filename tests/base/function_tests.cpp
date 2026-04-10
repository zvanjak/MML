#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Function.h"
#endif

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Core::FunctionTests
{
	/*********************************************************************/
	/*****                   RealFunction tests                      *****/
	/*********************************************************************/

	Real FunctionTests_TestFunc(Real x)
	{
		return sin(x) * (REAL(1.0) + REAL(0.5) * x * x);
	}

	TEST_CASE("Test_Function_func_pointer", "[simple]") {
			TEST_PRECISION_INFO();
		// creating a function object from a already existing (standalone) function
		RealFunction f1(FunctionTests_TestFunc);

		// or creating a function object directly
		RealFunction f2{ [](Real x) -> Real { return sin(x) * (REAL(1.0) + x * x / REAL(2.0)); } };

		REQUIRE(f1(REAL(1.0)) == FunctionTests_TestFunc(REAL(1.0)));
		REQUIRE(f2(REAL(1.0)) == FunctionTests_TestFunc(REAL(1.0)));
	}

	class ClassProvidingFuncToDerive
	{
	private:
		double _param;
	public:
		ClassProvidingFuncToDerive(double param) : _param(param) { }
		void setParam(double inParam) { _param = inParam; }

		Real operator()(Real x) const { return _param * sin(x); }
	};

	TEST_CASE("Test_Function_class_obj_overload_op()", "[simple]")
	{
		ClassProvidingFuncToDerive   funcObj(REAL(3.0));

		RealFunctionFromStdFunc f1(std::function<Real(Real)>{funcObj});

		REQUIRE(f1(REAL(1.0)) == Catch::Approx(REAL(3.0) * sin(REAL(1.0))));
		REQUIRE(funcObj(REAL(1.0)) == Catch::Approx(REAL(3.0) * sin(REAL(1.0))));

		REQUIRE(f1(-REAL(2.0)) == funcObj(-REAL(2.0)));
		REQUIRE(f1(REAL(1.0)) == funcObj(REAL(1.0)));
		REQUIRE(f1(REAL(5.0)) == funcObj(REAL(5.0)));

		// VERY IMPORTANT - our f1 has a COPY of funcObj, so if we change funcObj, f1 will not change
		funcObj.setParam(REAL(6.0));

		REQUIRE(REAL(REAL(3.0) * sin(REAL(1.0))) == f1(REAL(1.0)));   // WE HAVE AN OLD VALUE OF PARAMETAR!
		REQUIRE(REAL(REAL(6.0) * sin(REAL(1.0))) != f1(REAL(1.0)));
		REQUIRE(REAL(REAL(6.0) * sin(REAL(1.0))) == funcObj(REAL(1.0)));

		REQUIRE(f1(REAL(1.0)) != funcObj(REAL(1.0)));
		REQUIRE(f1(REAL(5.0)) != funcObj(REAL(5.0)));

		// if we create a new copy f2, it gets the new value of parametar
		RealFunctionFromStdFunc f2(std::function<Real(Real)>{funcObj});

		REQUIRE(REAL(REAL(6.0) * sin(REAL(1.0))) == f2(REAL(1.0)));
		REQUIRE(f2(REAL(1.0)) == funcObj(REAL(1.0)));
		REQUIRE(f2(REAL(5.0)) == funcObj(REAL(5.0)));
	}

	class ClassProvidingFuncToDerive2 : public IRealFunction
	{
	private:
		Real _param;
	public:
		ClassProvidingFuncToDerive2(Real param) : _param(param) { }
		void setParam(Real inParam) { _param = inParam; }

		Real operator()(Real x) const { return _param * sin(x); }
	};

	TEST_CASE("Test_Function_class_inheriting_IRealFunction", "[simple]")
	{
			TEST_PRECISION_INFO();
		ClassProvidingFuncToDerive2   funcObj(REAL(3.0));

		IRealFunction& f1 = funcObj;

		REQUIRE(f1(REAL(1.0)) == funcObj(REAL(1.0)));
		REQUIRE(f1(REAL(5.0)) == funcObj(REAL(5.0)));

		funcObj.setParam(REAL(6.0));

		REQUIRE(f1(REAL(1.0)) == funcObj(REAL(1.0)));
		REQUIRE(f1(REAL(5.0)) == funcObj(REAL(5.0)));
	}

	class BigComplexClassYouCantChange
	{
		// Has data and functionality for calculating function value we want
		// BUT, we have to manually implement exact calculation we need
	public:
		double _param;
		Real complexCalc(Real x) const { return cos(x); }
	};

	class BigComplexDerivFunc
	{
		const BigComplexClassYouCantChange& _ref;
	public:
		BigComplexDerivFunc(const BigComplexClassYouCantChange& bigClass) : _ref(bigClass) { }

		// here we implement function with behavior we want
		Real operator()(Real x)
		{
			return _ref._param * _ref.complexCalc(x);
		}
	};

	TEST_CASE("Test_Function_class_unchangeable_overload_op()", "[simple]")
	{
		BigComplexClassYouCantChange    ref;
		ref._param = REAL(3.0);

		BigComplexDerivFunc             funcObj(ref);
		RealFunctionFromStdFunc         f1(std::function<Real(Real)>{funcObj});

		REQUIRE(f1(REAL(1.0)) == funcObj(REAL(1.0)));  // f1 and funcObj must match exactly

		REQUIRE(f1(REAL(1.0)) == funcObj(REAL(1.0)));
		REQUIRE(f1(REAL(5.0)) == funcObj(REAL(5.0)));

		// f1 object is valid even if you changed the referenced class object
		ref._param = REAL(6.0);
		REQUIRE(f1(REAL(1.0)) == funcObj(REAL(1.0)));
		REQUIRE(f1(REAL(5.0)) == funcObj(REAL(5.0)));
	}

	/*********************************************************************/
	/*****                 ScalarFunction tests                    *****/
	/*********************************************************************/

	Real ScalarFunc3D(const VectorN<Real, 3>& v)
	{
		return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];  // x² + y² + z²
	}

	TEST_CASE("Function_ScalarFunction_func_pointer", "[Function]")
	{
		TEST_PRECISION_INFO();
		
		ScalarFunction<3> f(ScalarFunc3D);
		
		VectorN<Real, 3> v1({ REAL(1.0), REAL(2.0), REAL(3.0) });
		REQUIRE_THAT(f(v1), RealWithinAbs(REAL(14.0), TOL(1e-10, 1e-5)));  // 1 + 4 + 9 = 14
		
		VectorN<Real, 3> v2({ REAL(0.0), REAL(0.0), REAL(0.0) });
		REQUIRE_THAT(f(v2), RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		
		VectorN<Real, 3> v3({ REAL(1.0), REAL(1.0), REAL(1.0) });
		REQUIRE_THAT(f(v3), RealWithinAbs(REAL(3.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Function_ScalarFunctionFromStdFunc", "[Function]")
	{
		TEST_PRECISION_INFO();
		
		// Using lambda
		ScalarFunctionFromStdFunc<3> f([](const VectorN<Real, 3>& v) {
			return v[0] + v[1] + v[2];  // x + y + z
		});
		
		VectorN<Real, 3> v1({ REAL(1.0), REAL(2.0), REAL(3.0) });
		REQUIRE_THAT(f(v1), RealWithinAbs(REAL(6.0), TOL(1e-10, 1e-5)));
		
		VectorN<Real, 3> v2({ -REAL(1.0), REAL(1.0), REAL(0.0) });
		REQUIRE_THAT(f(v2), RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Function_ScalarFunction_2D", "[Function]")
	{
		TEST_PRECISION_INFO();
		
		// 2D scalar function: f(x,y) = x*y
		ScalarFunctionFromStdFunc<2> f([](const VectorN<Real, 2>& v) {
			return v[0] * v[1];
		});
		
		VectorN<Real, 2> v1({ REAL(3.0), REAL(4.0) });
		REQUIRE_THAT(f(v1), RealWithinAbs(REAL(12.0), TOL(1e-10, 1e-5)));
		
		VectorN<Real, 2> v2({ -REAL(2.0), REAL(5.0) });
		REQUIRE_THAT(f(v2), RealWithinAbs(-REAL(10.0), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****                 VectorFunction tests                    *****/
	/*********************************************************************/

	VectorN<Real, 3> VectorFunc3D(const VectorN<Real, 3>& v)
	{
		return VectorN<Real, 3>({ v[0] * REAL(2.0), v[1] * REAL(2.0), v[2] * REAL(2.0) });
	}

	TEST_CASE("Function_VectorFunction_func_pointer", "[Function]")
	{
		TEST_PRECISION_INFO();
		
		VectorFunction<3> f(VectorFunc3D);
		
		VectorN<Real, 3> v1({ REAL(1.0), REAL(2.0), REAL(3.0) });
		auto result = f(v1);
		
		REQUIRE_THAT(result[0], RealWithinAbs(REAL(2.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result[1], RealWithinAbs(REAL(4.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result[2], RealWithinAbs(REAL(6.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Function_VectorFunctionFromStdFunc", "[Function]")
	{
		TEST_PRECISION_INFO();
		
		// Rotation-like transformation
		VectorFunctionFromStdFunc<2> f([](const VectorN<Real, 2>& v) {
			return VectorN<Real, 2>({ -v[1], v[0] });  // 90-degree rotation
		});
		
		VectorN<Real, 2> v1({ REAL(1.0), REAL(0.0) });
		auto result = f(v1);
		
		REQUIRE_THAT(result[0], RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result[1], RealWithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****                 VectorNMFunction tests                    *****/
	/*********************************************************************/

	VectorN<Real, 2> VectorFuncNM(const VectorN<Real, 3>& v)
	{
		return VectorN<Real, 2>({ v[0] + v[1], v[1] + v[2] });
	}

	TEST_CASE("Function_VectorFunctionNM_func_pointer", "[Function]")
	{
		TEST_PRECISION_INFO();
		
		VectorFunctionNM<3, 2> f(VectorFuncNM);
		
		VectorN<Real, 3> v1({ REAL(1.0), REAL(2.0), REAL(3.0) });
		auto result = f(v1);
		
		REQUIRE(result.size() == 2);
		REQUIRE_THAT(result[0], RealWithinAbs(REAL(3.0), TOL(1e-10, 1e-5)));  // 1 + 2
		REQUIRE_THAT(result[1], RealWithinAbs(REAL(5.0), TOL(1e-10, 1e-5)));  // 2 + 3
	}

	TEST_CASE("Function_VectorFunctionNMFromStdFunc", "[Function]")
	{
		TEST_PRECISION_INFO();
		
		// R^2 -> R^3: (x,y) -> (x, y, x+y)
		VectorFunctionNMFromStdFunc<2, 3> f([](const VectorN<Real, 2>& v) {
			return VectorN<Real, 3>({ v[0], v[1], v[0] + v[1] });
		});
		
		VectorN<Real, 2> v1({ REAL(2.0), REAL(5.0) });
		auto result = f(v1);
		
		REQUIRE(result.size() == 3);
		REQUIRE_THAT(result[0], RealWithinAbs(REAL(2.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result[1], RealWithinAbs(REAL(5.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(result[2], RealWithinAbs(REAL(7.0), TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****                 ParametricCurve tests                    *****/
	/*********************************************************************/

	VectorN<Real, 3> HelixCurve(Real t)
	{
		return VectorN<Real, 3>({ std::cos(t), std::sin(t), t });
	}

	TEST_CASE("Function_ParametricCurve_func_pointer", "[Function]")
	{
		TEST_PRECISION_INFO();
		
		ParametricCurve<3> helix(HelixCurve);
		
		auto p0 = helix(REAL(0.0));
		REQUIRE_THAT(p0[0], RealWithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));  // cos(0) = 1
		REQUIRE_THAT(p0[1], RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));  // sin(0) = 0
		REQUIRE_THAT(p0[2], RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));  // z = 0
		
		auto pPi = helix(Constants::PI);
		REQUIRE_THAT(pPi[0], RealWithinAbs(-REAL(1.0), TOL(1e-10, 1e-5)));  // cos(π) = -1
		REQUIRE_THAT(pPi[1], RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));   // sin(π) ≈ 0
		REQUIRE_THAT(pPi[2], RealWithinAbs(Constants::PI, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Function_ParametricCurve_bounded", "[Function]")
	{
		TEST_PRECISION_INFO();
		
		ParametricCurve<3> helix(REAL(0.0), REAL(2.0) * Constants::PI, HelixCurve);
		
		REQUIRE_THAT(helix.getMinT(), RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(helix.getMaxT(), RealWithinAbs(REAL(2.0) * Constants::PI, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Function_ParametricCurveFromStdFunc", "[Function]")
	{
		TEST_PRECISION_INFO();
		
		// Circle in 2D
		ParametricCurveFromStdFunc<2> circle([](Real t) {
			return VectorN<Real, 2>({ std::cos(t), std::sin(t) });
		});
		
		auto p0 = circle(REAL(0.0));
		REQUIRE_THAT(p0[0], RealWithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p0[1], RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		
		auto pHalfPi = circle(Constants::PI / REAL(2.0));
		REQUIRE_THAT(pHalfPi[0], RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(pHalfPi[1], RealWithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Function_ParametricCurveFromStdFunc_bounded", "[Function]")
	{
		TEST_PRECISION_INFO();
		
		ParametricCurveFromStdFunc<2> arc(
			REAL(0.0), Constants::PI,
			[](Real t) { return VectorN<Real, 2>({ std::cos(t), std::sin(t) }); }
		);
		
		REQUIRE_THAT(arc.getMinT(), RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(arc.getMaxT(), RealWithinAbs(Constants::PI, TOL(1e-10, 1e-5)));
	}

	/*********************************************************************/
	/*****                 ParametricSurface tests                    *****/
	/*********************************************************************/

	VectorN<Real, 3> SphereSurface(Real theta, Real phi)
	{
		return VectorN<Real, 3>({ 
			std::sin(theta) * std::cos(phi),
			std::sin(theta) * std::sin(phi),
			std::cos(theta)
		});
	}

	Real ConstantZero(Real) { return REAL(0.0); }
	Real ConstantTwoPi(Real) { return REAL(2.0) * Constants::PI; }

	TEST_CASE("Function_ParametricSurface_func_pointer", "[Function]")
	{
		TEST_PRECISION_INFO();
		
		// Unit sphere: θ ∈ [0, π], φ ∈ [0, 2π]
		ParametricSurface<3> sphere(SphereSurface, 
			REAL(0.0), Constants::PI, 
			ConstantZero, ConstantTwoPi);
		
		// North pole: θ = 0
		auto northPole = sphere(REAL(0.0), REAL(0.0));
		REQUIRE_THAT(northPole[0], RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(northPole[1], RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(northPole[2], RealWithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
		
		// Equator point: θ = π/2, φ = 0
		auto equator = sphere(Constants::PI / REAL(2.0), REAL(0.0));
		REQUIRE_THAT(equator[0], RealWithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(equator[1], RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(equator[2], RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Function_ParametricSurface_domain_bounds", "[Function]")
	{
		TEST_PRECISION_INFO();
		
		ParametricSurface<3> sphere(SphereSurface, 
			REAL(0.0), Constants::PI, 
			ConstantZero, ConstantTwoPi);
		
		REQUIRE_THAT(sphere.getMinU(), RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(sphere.getMaxU(), RealWithinAbs(Constants::PI, TOL(1e-10, 1e-5)));
		REQUIRE_THAT(sphere.getMinW(REAL(0.5)), RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(sphere.getMaxW(REAL(0.5)), RealWithinAbs(REAL(2.0) * Constants::PI, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Function_ParametricSurface_domain_error", "[Function]")
	{
		TEST_PRECISION_INFO();
		
		ParametricSurface<3> sphere(SphereSurface, 
			REAL(0.0), Constants::PI, 
			ConstantZero, ConstantTwoPi);
		
		// x out of domain
		REQUIRE_THROWS_AS(sphere(-REAL(0.1), REAL(0.0)), DomainError);
		REQUIRE_THROWS_AS(sphere(Constants::PI + REAL(0.1), REAL(0.0)), DomainError);
		
		// y out of domain  
		REQUIRE_THROWS_AS(sphere(REAL(0.5), -REAL(0.1)), DomainError);
		REQUIRE_THROWS_AS(sphere(REAL(0.5), REAL(2.0) * Constants::PI + REAL(0.1)), DomainError);
	}

	VectorN<Real, 3> PlaneSurface(Real u, Real w)
	{
		return VectorN<Real, 3>({ u, w, REAL(0.0) });
	}

	TEST_CASE("Function_ParametricSurfaceRect_func_pointer", "[Function]")
	{
		TEST_PRECISION_INFO();
		
		ParametricSurfaceRect<3> plane(PlaneSurface, REAL(0.0), REAL(1.0), REAL(0.0), REAL(1.0));
		
		auto p00 = plane(REAL(0.0), REAL(0.0));
		REQUIRE_THAT(p00[0], RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p00[1], RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p00[2], RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		
		auto p11 = plane(REAL(1.0), REAL(1.0));
		REQUIRE_THAT(p11[0], RealWithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p11[1], RealWithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(p11[2], RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Function_ParametricSurfaceRect_bounds", "[Function]")
	{
		TEST_PRECISION_INFO();
		
		ParametricSurfaceRect<3> plane(PlaneSurface, -REAL(1.0), REAL(2.0), -REAL(3.0), REAL(4.0));
		
		REQUIRE_THAT(plane.getMinU(), RealWithinAbs(-REAL(1.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(plane.getMaxU(), RealWithinAbs(REAL(2.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(plane.getMinW(), RealWithinAbs(-REAL(3.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(plane.getMaxW(), RealWithinAbs(REAL(4.0), TOL(1e-10, 1e-5)));
	}

	TEST_CASE("Function_ParametricSurfaceRect_unbounded", "[Function]")
	{
		TEST_PRECISION_INFO();
		
		ParametricSurfaceRect<3> plane(PlaneSurface);
		
		REQUIRE(plane.getMinU() == Constants::NegInf);
		REQUIRE(plane.getMaxU() == Constants::PosInf);
		REQUIRE(plane.getMinW() == Constants::NegInf);
		REQUIRE(plane.getMaxW() == Constants::PosInf);
	}

	TEST_CASE("Function_ParametricSurfaceFromStdFunc", "[Function]")
	{
		TEST_PRECISION_INFO();
		
		// Saddle surface: z = x² - y²
		std::function<VectorN<Real, 3>(Real, Real)> saddleFunc = 
			[](Real u, Real w) { return VectorN<Real, 3>({ u, w, u*u - w*w }); };
		
		ParametricSurfaceFromStdFunc<3> saddle(saddleFunc, -REAL(1.0), REAL(1.0), -REAL(1.0), REAL(1.0));
		
		auto origin = saddle(REAL(0.0), REAL(0.0));
		REQUIRE_THAT(origin[0], RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(origin[1], RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		REQUIRE_THAT(origin[2], RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));
		
		auto p11 = saddle(REAL(1.0), REAL(1.0));
		REQUIRE_THAT(p11[2], RealWithinAbs(REAL(0.0), TOL(1e-10, 1e-5)));  // 1² - 1² = 0
		
		auto p10 = saddle(REAL(1.0), REAL(0.0));
		REQUIRE_THAT(p10[2], RealWithinAbs(REAL(1.0), TOL(1e-10, 1e-5)));  // 1² - 0² = 1
	}

} // namespace Tests::Core::FunctionTests
