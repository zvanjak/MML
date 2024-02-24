#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "core/Function.h"
#include "core/Derivation.h"
#endif

#include "../test_data/real_functions_test_bed.h"
#include "../test_data/scalar_functions_test_bed.h"
#include "../test_data/vector_functions_test_bed.h"

using namespace MML;


////////////////////////////////////////////////////////////////////////////////////////////////////
///////         Derivation of REAL functions
////////////////////////////////////////////////////////////////////////////////////////////////////

// NDer1 achieves precision of at least 5 decimal places on a set of functions
// definirati parove <RealFunction, Interval>, init iz TestBeda
//  onda za svaki ekvidistan covering intervala, izracunati NDer1 i usporediti sa analitickim derivacijama
TEST_CASE("Test_NDer_sin_func", "[simple]") 
{
    RealFunction sinFunc = TestBeds::RealFunctionsTestBed::getTestFunctionReal("Sin")._func;
    RealFunction sinFuncDer = TestBeds::RealFunctionsTestBed::getTestFunctionReal("Sin")._funcDerived;

    double der1 = Derivation::NDer1(sinFunc, 0.5);
    REQUIRE(sinFuncDer(0.5) == Approx(der1).epsilon(1e-8));

    double der2 = Derivation::NDer2(sinFunc, 0.5);
    REQUIRE(sinFuncDer(0.5) == Approx(der2).epsilon(1e-11));

    double der4 = Derivation::NDer4(sinFunc, 0.5);
    REQUIRE(sinFuncDer(0.5) == Approx(der4).epsilon(1e-13));

    double der6 = Derivation::NDer6(sinFunc, 0.5);
    REQUIRE(sinFuncDer(0.5) == Approx(der6).epsilon(1e-13));

    double der8 = Derivation::NDer8(sinFunc, 0.5);
    REQUIRE(sinFuncDer(0.5) == Approx(der8).epsilon(1e-14));
}

// TODO 0.9 - see results for complex and hard test function
TEST_CASE("Test_NDer_TestDer1", "[simple]") 
{
    RealFunction testFunc = TestBeds::RealFunctionsTestBed::getTestFunctionRealWithDerivation("TestDer1")._func;
    RealFunction testFuncDer = TestBeds::RealFunctionsTestBed::getTestFunctionRealWithDerivation("TestDer1")._funcDerived;
}

TEST_CASE("Test_NSecDer", "[simple]") 
{
    RealFunction sinFunc = TestBeds::RealFunctionsTestBed::getTestFunctionReal("Sin")._func;
    RealFunction sinFuncSecDer = TestBeds::RealFunctionsTestBed::getTestFunctionReal("Sin")._funcSecDer;

    double secder2 = Derivation::NSecDer2(sinFunc, 0.5);
    REQUIRE(sinFuncSecDer(0.5) == Approx(secder2).epsilon(1e-9));

    double secder4 = Derivation::NSecDer4(sinFunc, 0.5);
    REQUIRE(sinFuncSecDer(0.5) == Approx(secder4).epsilon(1e-11));

    double secder6 = Derivation::NSecDer6(sinFunc, 0.5);
    REQUIRE(sinFuncSecDer(0.5) == Approx(secder6).epsilon(1e-11));

    double secder8 = Derivation::NSecDer8(sinFunc, 0.5);
    REQUIRE(sinFuncSecDer(0.5) == Approx(secder8).epsilon(1e-13));    
}

TEST_CASE("Test_NThirdDer", "[simple]") 
{
    RealFunction sinFunc = TestBeds::RealFunctionsTestBed::getTestFunctionReal("Sin")._func;
    RealFunction sinFuncThirdDer = TestBeds::RealFunctionsTestBed::getTestFunctionReal("Sin")._funcThirdDer;

    double thirdder1 = Derivation::NThirdDer2(sinFunc, 0.5);
    REQUIRE(sinFuncThirdDer(0.5) == Approx(thirdder1).epsilon(1e-6));

    double thirdder4 = Derivation::NThirdDer4(sinFunc, 0.5);
    REQUIRE(sinFuncThirdDer(0.5) == Approx(thirdder4).epsilon(1e-9));

    double thirdder6 = Derivation::NThirdDer6(sinFunc, 0.5);
    REQUIRE(sinFuncThirdDer(0.5) == Approx(thirdder6).epsilon(1e-10));

    double thirdder8 = Derivation::NThirdDer8(sinFunc, 0.5);
    REQUIRE(sinFuncThirdDer(0.5) == Approx(thirdder8).epsilon(1e-12)); 
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///////         Derivation of SCALAR functions
////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Test_NDerPartial", "[simple]") 
{
    VectorN<Real, 3> point{1.0, 1.0, 3.0 };

    ScalarFunction<3> f = TestBeds::ScalarFunctionsTestBed::getTestFunctionScalar3("Scalar func 2")._func;
    Real (*fDer)(const VectorN<Real, 3> &, int ind) = TestBeds::ScalarFunctionsTestBed::getTestFunctionScalar3("Scalar func 2")._funcDerived;

    Real der_x, der_y, der_z;

    der_x = Derivation::NDer1Partial(f, 0, point);
    der_y = Derivation::NDer1Partial(f, 1, point);
    der_z = Derivation::NDer1Partial(f, 2, point);

    REQUIRE(fDer(point, 0) == Approx(der_x).epsilon(1e-7));
    REQUIRE(fDer(point, 1) == Approx(der_y).epsilon(1e-7));
    REQUIRE(fDer(point, 2) == Approx(der_z).epsilon(1e-7));

    der_x = Derivation::NDer2Partial(f, 0, point);
    der_y = Derivation::NDer2Partial(f, 1, point);
    der_z = Derivation::NDer2Partial(f, 2, point);

    REQUIRE(fDer(point, 0) == Approx(der_x).epsilon(1e-9));
    REQUIRE(fDer(point, 1) == Approx(der_y).epsilon(1e-11));
    REQUIRE(fDer(point, 2) == Approx(der_z).epsilon(1e-11));

    der_x = Derivation::NDer4Partial(f, 0, point);
    der_y = Derivation::NDer4Partial(f, 1, point);
    der_z = Derivation::NDer4Partial(f, 2, point);

    REQUIRE(fDer(point, 0) == Approx(der_x).epsilon(1e-12));
    REQUIRE(fDer(point, 1) == Approx(der_y).epsilon(1e-10));
    REQUIRE(fDer(point, 2) == Approx(der_z).epsilon(1e-14));

    der_x = Derivation::NDer6Partial(f, 0, point);
    der_y = Derivation::NDer6Partial(f, 1, point);
    der_z = Derivation::NDer6Partial(f, 2, point);

    REQUIRE(fDer(point, 0) == Approx(der_x).epsilon(1e-12));
    REQUIRE(fDer(point, 1) == Approx(der_y).epsilon(1e-13));
    REQUIRE(fDer(point, 2) == Approx(der_z).epsilon(1e-14));

    der_x = Derivation::NDer8Partial(f, 0, point);
    der_y = Derivation::NDer8Partial(f, 1, point);
    der_z = Derivation::NDer8Partial(f, 2, point);

    REQUIRE(fDer(point, 0) == Approx(der_x).epsilon(1e-11));
    REQUIRE(fDer(point, 1) == Approx(der_y).epsilon(1e-8));
    REQUIRE(fDer(point, 2) == Approx(der_z).epsilon(1e-13));
}

// TODO - test for 2nd partial derivation of scalar function
TEST_CASE("Test_NSecDerPartial", "[simple]") 
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///////         Derivation of VECTOR functions
////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Test_NDer1VecPartial", "[simple]") 
{
    VectorN<Real, 3> point{1.0, 3.0, -1.0 };

    VectorFunction<3> f = TestBeds::VectorFunctionsTestBed::getTestFunctionVector(0)._func;
    VectorN<Real, 3> (*fDer)(const VectorN<Real, 3> &, int ind) = TestBeds::VectorFunctionsTestBed::getTestFunctionVector(0)._funcDerived;

    double der_x_x, der_x_y, der_x_z;
    double der_y_x, der_y_y, der_y_z;
    double der_z_x, der_z_y, der_z_z;

    der_x_x = Derivation::NDer1Partial(f, 0, 0, point);
    der_x_y = Derivation::NDer1Partial(f, 0, 1, point);
    der_x_z = Derivation::NDer1Partial(f, 0, 2, point);

    REQUIRE(fDer(point, 0)[0] == Approx(der_x_x).epsilon(1e-8));
    REQUIRE(fDer(point, 0)[1] == Approx(der_x_y).epsilon(1e-7));
    REQUIRE(fDer(point, 0)[2] == Approx(der_x_z).epsilon(1e-7));

    der_y_x = Derivation::NDer1Partial(f, 1, 0, point);
    der_y_y = Derivation::NDer1Partial(f, 1, 1, point);
    der_y_z = Derivation::NDer1Partial(f, 1, 2, point);

    REQUIRE(fDer(point, 1)[0] == Approx(der_y_x).epsilon(1e-7));
    REQUIRE(fDer(point, 1)[1] == Approx(der_y_y).epsilon(1e-8));
    REQUIRE(fDer(point, 1)[2] == Approx(der_y_z).epsilon(1e-8));

    der_z_x = Derivation::NDer1Partial(f, 2, 0, point);
    der_z_y = Derivation::NDer1Partial(f, 2, 1, point);
    der_z_z = Derivation::NDer1Partial(f, 2, 2, point);

    REQUIRE(fDer(point, 2)[0] == Approx(der_z_x).epsilon(1e-7));
    REQUIRE(fDer(point, 2)[1] == Approx(der_z_y).epsilon(1e-8));
    REQUIRE(fDer(point, 2)[2] == Approx(der_z_z).epsilon(1e-7));
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///////         Derivation of Tensor field
////////////////////////////////////////////////////////////////////////////////////////////////////
// TODO - tests for tensor field derivation

////////////////////////////////////////////////////////////////////////////////////////////////////
///////         Derivation of Parametric curves
////////////////////////////////////////////////////////////////////////////////////////////////////

// TODO - test for derivation of Parametric curves
TEST_CASE("Test_NDer_curve", "[simple]") 
{
}

TEST_CASE("Test_NSecDer_curve", "[simple]") 
{
}

TEST_CASE("Test_NThirdDer_curve", "[simple]") 
{
}

TEST_CASE("Verify_vector_equation1")
{

}