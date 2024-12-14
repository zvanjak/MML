#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Integration.h"

#endif

#include "../test_data/real_functions_test_bed.h"

using namespace MML;

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Core::IntegrationTests
{
	TEST_CASE("Test_Integration_Trap_const_func", "[simple]")
	{
		RealFunction constFunc0{ [](double x) { return 0.0; } };
		RealFunction constFunc1{ [](double x) { return 2.0; } };
		RealFunction constFunc2{ [](double x) { return -5.0; } };
		RealFunction constFunc3{ [](double x) { return 1e-10; } };

		REQUIRE_THAT(0.0, WithinAbs(IntegrateTrap(constFunc0, 0.0, 1.0), 1e-15));
		REQUIRE_THAT(0.0, WithinAbs(IntegrateTrap(constFunc0, -3.0, 1.0), 1e-15));

		REQUIRE_THAT(2.0, WithinAbs(IntegrateTrap(constFunc1, 0.0, 1.0), 1e-15));
		REQUIRE_THAT(8.0, WithinAbs(IntegrateTrap(constFunc1, -3.0, 1.0), 1e-15));

		REQUIRE_THAT(-5.0, WithinAbs(IntegrateTrap(constFunc2, 0.0, 1.0), 1e-15));
		REQUIRE_THAT(-20.0, WithinAbs(IntegrateTrap(constFunc2, -3.0, 1.0), 1e-15));

		REQUIRE_THAT(1e-10, WithinAbs(IntegrateTrap(constFunc3, 0.0, 1.0), 1e-15));
		REQUIRE_THAT(4e-10, WithinAbs(IntegrateTrap(constFunc3, -3.0, 1.0), 1e-15));
	}

	// TODO
	//TEST_CASE("Test_Integration_almost_const_func", "[simple]") 
	//{
	//	RealFunction linearFun{ [](double x) { return 2.0; } };

	//	double int_trap = IntegrateTrap(linearFun, 0.0, 1.0);

	//	REQUIRE_THAT(1.0, Catch::Matchers::WithinAbs(IntegrateTrap(linearFun, 0.0, 1.0), 1e-16));
	//}

	TEST_CASE("Test_Integration_Trap_linear_func", "[simple]")
	{
		RealFunction linearFunc1{ [](double x) { return 2 * x; } };
		RealFunction linearFunc2{ [](double x) { return -3 * x; } };
		RealFunction linearFunc3{ [](double x) { return 1e-10 * x; } };
		RealFunction linearFunc4{ [](double x) { return 2 * (x - 0.5); } };

		REQUIRE_THAT(1.0, WithinAbs(IntegrateTrap(linearFunc1, 0.0, 1.0), 1e-15));
		REQUIRE_THAT(-1.5, WithinAbs(IntegrateTrap(linearFunc2, 0.0, 1.0), 1e-15));
		REQUIRE_THAT(0.5e-10, WithinAbs(IntegrateTrap(linearFunc3, 0.0, 1.0), 1e-15));
		REQUIRE_THAT(0.0, WithinAbs(IntegrateTrap(linearFunc4, 0.0, 1.0), 1e-15));

		REQUIRE_THAT(-5.0, WithinAbs(IntegrateTrap(linearFunc1, -3.0, 2.0), 1e-15));
		REQUIRE_THAT(7.5, WithinAbs(IntegrateTrap(linearFunc2, -3.0, 2.0), 1e-15));
		REQUIRE_THAT(-2.5e-10, WithinAbs(IntegrateTrap(linearFunc3, -3.0, 2.0), 1e-15));
		REQUIRE_THAT(-10.0, WithinAbs(IntegrateTrap(linearFunc4, -3.0, 2.0), 1e-15));
	}

	TEST_CASE("Test_Integration_Trap_standard_func_precision", "[simple]")
	{
		RealFunction sinFunc = TestBeds::RealFunctionsTestBed::getTestFunctionReal("Sin")._func;
		RealFunction sinFunc_int = TestBeds::RealFunctionsTestBed::getTestFunctionReal("Sin")._funcIntegrated;

		double a = 0.0, b = 1.0;
		REQUIRE_THAT(sinFunc_int(b) - sinFunc_int(a), WithinAbs(IntegrateTrap(sinFunc, a, b), 1e-6));
		REQUIRE_THAT(sinFunc_int(b) - sinFunc_int(a), !WithinAbs(IntegrateTrap(sinFunc, a, b), 1e-7));

		// row - integration method
		// column - req prec
		// testirati mogucu preciznost za dvije funkcije
		// sa svim metodama
		// provjeriti exception to many steps za preveliku preciznost

	}

	TEST_CASE("Test_Integration_Trap_integral_test_bed1", "[simple]")
	{
		RealFunction testFunc = TestBeds::RealFunctionsTestBed::getTestFunctionRealWithIntegral("TestInt1")._func;
		RealFunction testFunc_int = TestBeds::RealFunctionsTestBed::getTestFunctionRealWithIntegral("TestInt1")._funcIntegrated;

		double a = 0.0, b = 1.0;
		double integralExactVal = testFunc_int(b) - testFunc_int(a);
		REQUIRE_THAT(integralExactVal,  WithinAbs(IntegrateTrap(testFunc, a, b), 1e-6));
		REQUIRE_THAT(integralExactVal, !WithinAbs(IntegrateTrap(testFunc, a, b), 1e-7));

		REQUIRE_THAT(integralExactVal,  WithinRel(IntegrateTrap(testFunc, a, b), 1e-5));
		REQUIRE_THAT(integralExactVal, !WithinRel(IntegrateTrap(testFunc, a, b), 1e-6));

		a = 0.0, b = 5.0;
		integralExactVal = testFunc_int(b) - testFunc_int(a);
		REQUIRE_THAT(integralExactVal,  WithinAbs(IntegrateTrap(testFunc, a, b), 1e-3));
		REQUIRE_THAT(integralExactVal, !WithinAbs(IntegrateTrap(testFunc, a, b), 1e-4));

		REQUIRE_THAT(integralExactVal,  WithinRel(IntegrateTrap(testFunc, a, b), 1e-5));
		REQUIRE_THAT(integralExactVal, !WithinRel(IntegrateTrap(testFunc, a, b), 1e-6));

		a = -5.0, b = 2.0;
		integralExactVal = testFunc_int(b) - testFunc_int(a);
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinAbs(IntegrateTrap(testFunc, a, b), 1e-2));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinAbs(IntegrateTrap(testFunc, a, b), 1e-3));

		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinRel(IntegrateTrap(testFunc, a, b), 1e-5));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinRel(IntegrateTrap(testFunc, a, b), 1e-6));
	}

	TEST_CASE("Test_Integration_Trap_integral_test_bed2", "[simple]")
	{
		RealFunction testFunc = TestBeds::RealFunctionsTestBed::getTestFunctionRealWithIntegral("TestInt2")._func;
		RealFunction testFunc_int = TestBeds::RealFunctionsTestBed::getTestFunctionRealWithIntegral("TestInt2")._funcIntegrated;

		double a = 0.0, b = 1.0;
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinAbs(IntegrateTrap(testFunc, a, b), 1e-6));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinAbs(IntegrateTrap(testFunc, a, b), 1e-7));

		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinRel(IntegrateTrap(testFunc, a, b), 1e-5));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinRel(IntegrateTrap(testFunc, a, b), 1e-6));

		a = 0.0, b = 5.0;
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinAbs(IntegrateTrap(testFunc, a, b), 1e-5));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinAbs(IntegrateTrap(testFunc, a, b), 1e-6));

		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinRel(IntegrateTrap(testFunc, a, b), 1e-5));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinRel(IntegrateTrap(testFunc, a, b), 1e-6));

		a = -5.0, b = 2.0;
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinAbs(IntegrateTrap(testFunc, a, b), 1e-5));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinAbs(IntegrateTrap(testFunc, a, b), 1e-6));

		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), WithinRel(IntegrateTrap(testFunc, a, b), 1e-5));
		REQUIRE_THAT(testFunc_int(b) - testFunc_int(a), !WithinRel(IntegrateTrap(testFunc, a, b), 1e-6));
	}
}