#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Integration.h"

#endif

#include "../test_data/real_functions_test_bed.h"

using namespace MML;


TEST_CASE("Test_Integration", "[simple]") {
    RealFunction testFunc = TestBeds::RealFunctionsTestBed::getTestFunctionRealWithIntegral("TestInt1")._func;
    RealFunction integratedTestFunc = TestBeds::RealFunctionsTestBed::getTestFunctionRealWithIntegral("TestInt1")._funcIntegrated;

    double a = 0.0;
    double b = 1.0;
    double int_trap = Integration::IntegrateTrap(testFunc,a,b);
    double int_simp = Integration::IntegrateSimpson(testFunc,a,b);
    double int_romb = Integration::IntegrateRomberg(testFunc,a,b);

	REQUIRE(int_trap == Approx(integratedTestFunc(b) - integratedTestFunc(a)).epsilon(1e-5));
	REQUIRE(int_simp == Approx(integratedTestFunc(b) - integratedTestFunc(a)).epsilon(1e-7));
	REQUIRE(int_romb == Approx(integratedTestFunc(b) - integratedTestFunc(a)).epsilon(1e-9));
    // row - integration method
    // column - req prec
    // testirati mogucu preciznost za dvije funkcije
    // sa svim metodama
    // provjeriti exception to many steps za preveliku preciznost
    
}