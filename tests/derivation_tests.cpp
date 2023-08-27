#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Vector.h"
#include "basic_types/Function.h"
#include "algorithms/Derivation.h"
#endif

#include "../test_data/real_functions_test_bed.h"
#include "../test_data/scalar_functions_test_bed.h"
#include "../test_data/vector_functions_test_bed.h"

using namespace MML;

TEST_CASE("Test_sin_derivation_precision", "[simple]") 
{
    RealFunction  f_sin = TestData::RealFunctionsTestBed::getTestFunctionReal(0)._func;
    RealFunction  f_der = TestData::RealFunctionsTestBed::getTestFunctionReal(0)._funcDerived;

    RealFunction  f_cos = TestData::RealFunctionsTestBed::getTestFunctionReal(1)._func;

    Vector<Real> x_val{0.0, 0.5, 1.0, 2.0, 3.0, 5.0 };
    
    Vector<Real> nder1_prec{1e-16,  1e-8,  1e-7,  1e-7,  1e-8,  1e-7 };
    Vector<Real> nder2_prec{1e-10, 1e-11, 1e-10, 1e-10, 1e-10, 1e-10 };
    Vector<Real> nder4_prec{1e-13, 1e-13, 1e-12, 1e-13, 1e-13, 1e-12 };
    Vector<Real> nder6_prec{1e-15, 1e-13, 1e-14, 1e-13, 1e-13, 1e-13 };
    Vector<Real> nder8_prec{1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14 };
    
    for( int i=0; i<x_val.size(); i++ )
    {
        double x = x_val[i];

        double num_der = Derivation::NDer1(f_sin, x);
        REQUIRE(f_der(x) == Approx(num_der).epsilon(nder1_prec[i]));
        REQUIRE(f_der(x) != Approx(num_der).epsilon(nder1_prec[i]*0.1));

        num_der = Derivation::NDer2(f_sin, x);
        REQUIRE(f_der(x) == Approx(num_der).epsilon(nder2_prec[i]));
        REQUIRE(f_der(x) != Approx(num_der).epsilon(nder2_prec[i]*0.1));

        num_der = Derivation::NDer4(f_sin, x);
        REQUIRE(f_der(x) == Approx(num_der).epsilon(nder4_prec[i]));
        REQUIRE(f_der(x) != Approx(num_der).epsilon(nder4_prec[i]*0.1));  

        num_der = Derivation::NDer6(f_sin, x);
        REQUIRE(f_der(x) == Approx(num_der).epsilon(nder6_prec[i]));
        REQUIRE(f_der(x) != Approx(num_der).epsilon(nder6_prec[i]*0.1));              

        num_der = Derivation::NDer8(f_sin, x);
        REQUIRE(f_der(x) == Approx(num_der).epsilon(nder8_prec[i]));
        REQUIRE(f_der(x) != Approx(num_der).epsilon(nder8_prec[i]*0.1));              
    }   
}

// TODO - implement derivation precision test as above
////////////////////////////////////////////////////////////////////////////////////////////////////
///////         Derivation of REAL functions
////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Test_NDer1_with_different_h", "[simple]") 
{
    RealFunction  f = TestData::RealFunctionsTestBed::getTestFunctionReal(0)._func;
    RealFunction  fDer = TestData::RealFunctionsTestBed::getTestFunctionReal(0)._funcDerived;

    double x = 1.0;

    double der = Derivation::NDer1(f, x, 1e-4);
    REQUIRE(fDer(x) == Approx(der).epsilon(1e-4));
    REQUIRE(fDer(x) != Approx(der).epsilon(1e-5));

    der = Derivation::NDer1(f, x, 1e-5);
    REQUIRE(fDer(x) == Approx(der).epsilon(1e-5));
    REQUIRE(fDer(x) != Approx(der).epsilon(1e-6));

    der = Derivation::NDer1(f, x, 1e-6);
    REQUIRE(fDer(x) == Approx(der).epsilon(1e-6));
    REQUIRE(fDer(x) != Approx(der).epsilon(1e-7));
    
    der = Derivation::NDer1(f, x, 1e-7);
    REQUIRE(fDer(x) == Approx(der).epsilon(1e-7));
    REQUIRE(fDer(x) != Approx(der).epsilon(1e-8));
    
    der = Derivation::NDer1(f, x, 1e-8);
    REQUIRE(fDer(x) == Approx(der).epsilon(1e-8));
    REQUIRE(fDer(x) != Approx(der).epsilon(1e-9));    
}

TEST_CASE("Test_NDer2_with_different_h", "[simple]") 
{
    RealFunction  f = TestData::RealFunctionsTestBed::getTestFunctionReal(0)._func;
    RealFunction  fDer = TestData::RealFunctionsTestBed::getTestFunctionReal(0)._funcDerived;

    // different step sizes for h!
    double der = Derivation::NDer2(f, 1.0, 1e-4);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-8));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-9));

    der = Derivation::NDer2(f, 1.0, 1e-5);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-10));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-11));

    der = Derivation::NDer2(f, 1.0, 1e-6);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-10));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-11));
    
    der = Derivation::NDer2(f, 1.0, 1e-7);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-9));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-10));

    der = Derivation::NDer2(f, 1.0, 1e-8);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-8));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-9));
}

TEST_CASE("Test_NDer4_with_different_h", "[simple]") 
{
    RealFunction  f = TestData::RealFunctionsTestBed::getTestFunctionReal(0)._func;
    RealFunction  fDer = TestData::RealFunctionsTestBed::getTestFunctionReal(0)._funcDerived;

    // different step sizes for h!
    double der = Derivation::NDer4(f, 1.0, 1e-4);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-8));
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-11));

    der = Derivation::NDer4(f, 1.0, 1e-5);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-11));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-12));

    der = Derivation::NDer4(f, 1.0, 1e-6);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-10));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-11));
    
    der = Derivation::NDer4(f, 1.0, 1e-7);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-9));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-10));

    der = Derivation::NDer4(f, 1.0, 1e-8);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-8));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-9));
}

TEST_CASE("Test_NDer6_with_different_h", "[simple]") 
{
    RealFunction  f = TestData::RealFunctionsTestBed::getTestFunctionReal(0)._func;
    RealFunction  fDer = TestData::RealFunctionsTestBed::getTestFunctionReal(0)._funcDerived;

    // different step sizes for h!
    double der = Derivation::NDer6(f, 1.0, 1e-4);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-8));
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-11));

    der = Derivation::NDer6(f, 1.0, 1e-5);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-11));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-12));

    der = Derivation::NDer6(f, 1.0, 1e-6);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-10));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-11));
    
    der = Derivation::NDer6(f, 1.0, 1e-7);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-9));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-10));

    der = Derivation::NDer6(f, 1.0, 1e-8);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-8));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-9));
}

TEST_CASE("Test_NDer8_with_different_h", "[simple]") 
{
    RealFunction  f = TestData::RealFunctionsTestBed::getTestFunctionReal(0)._func;
    RealFunction  fDer = TestData::RealFunctionsTestBed::getTestFunctionReal(0)._funcDerived;

    // different step sizes for h!
    double der = Derivation::NDer8(f, 1.0, 1e-4);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-8));
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-11));

    der = Derivation::NDer8(f, 1.0, 1e-5);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-11));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-12));

    der = Derivation::NDer8(f, 1.0, 1e-6);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-10));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-11));
    
    der = Derivation::NDer8(f, 1.0, 1e-7);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-9));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-10));

    der = Derivation::NDer8(f, 1.0, 1e-8);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-8));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-9));
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///////         Derivation of SCALAR functions
////////////////////////////////////////////////////////////////////////////////////////////////////
TEST_CASE("Test_NDer1Partial", "[simple]") 
{
    VectorN<Real, 3> point{1.0, 1.0, 3.0 };

    ScalarFunction<3> f = TestData::ScalarFunctionsTestBed::_listFuncScalar3[0]._func;
    double (*fDer)(const VectorN<Real, 3> &, int ind) = TestData::ScalarFunctionsTestBed::_listFuncScalar3[0]._funcDerived;

    // prva komponenta je cos() pa očekujemo derivaciju -sin()
    double der1 = Derivation::NDer1Partial(f, 0, point, 1e-6);
    REQUIRE(fDer(point, 0) == Approx(der1).epsilon(1e-6));
    REQUIRE(fDer(point, 0) != Approx(der1).epsilon(1e-7));

    // druga  komponenta je sin() pa očekujemo derivaciju cos()
    double der2 = Derivation::NDer1Partial(f, 1, point, 1e-6);
    REQUIRE(fDer(point, 1) == Approx(der2).epsilon(1e-6));
    REQUIRE(fDer(point, 1) != Approx(der2).epsilon(1e-7));

    // treća komponenta je exp() pa očekujemo derivaciju exp()
    double der3 = Derivation::NDer1Partial(f, 2, point, 1e-6);
    REQUIRE(fDer(point, 2) == Approx(der3).epsilon(1e-6));
    REQUIRE(fDer(point, 2) != Approx(der3).epsilon(1e-7));
}

TEST_CASE("Test_NDer2Partial", "[simple]") 
{
    VectorN<Real, 3> point{1.0, 1.0, 3.0 };

    ScalarFunction<3> f = TestData::ScalarFunctionsTestBed::_listFuncScalar3[0]._func;
    double (*fDer)(const VectorN<Real, 3> &, int ind) = TestData::ScalarFunctionsTestBed::_listFuncScalar3[0]._funcDerived;

    // prva komponenta je cos() pa očekujemo derivaciju -sin()
    double der1 = Derivation::NDer2Partial(f, 0, point, 1e-6);
    REQUIRE(fDer(point, 0) == Approx(der1).epsilon(1e-9));
    REQUIRE(fDer(point, 0) != Approx(der1).epsilon(1e-10));

    // druga  komponenta je sin() pa očekujemo derivaciju cos()
    double der2 = Derivation::NDer2Partial(f, 1, point, 1e-6);
    REQUIRE(fDer(point, 1) == Approx(der2).epsilon(1e-9));
    REQUIRE(fDer(point, 1) != Approx(der2).epsilon(1e-10));

    // treća komponenta je exp() pa očekujemo derivaciju exp()
    double der3 = Derivation::NDer2Partial(f, 2, point, 1e-6);
    REQUIRE(fDer(point, 2) == Approx(der3).epsilon(1e-9));
    REQUIRE(fDer(point, 2) != Approx(der3).epsilon(1e-10));
}

TEST_CASE("Test_NDer4Partial", "[simple]") 
{
    VectorN<Real, 3> point{1.0, 1.0, 3.0 };

    ScalarFunction<3> f = TestData::ScalarFunctionsTestBed::_listFuncScalar3[0]._func;
    double (*fDer)(const VectorN<Real, 3> &, int ind) = TestData::ScalarFunctionsTestBed::_listFuncScalar3[0]._funcDerived;

    // prva komponenta je cos() pa očekujemo derivaciju -sin()
    double der1 = Derivation::NDer4Partial(f, 0, point, 1e-6);
    REQUIRE(fDer(point, 0) == Approx(der1).epsilon(1e-9));
    REQUIRE(fDer(point, 0) != Approx(der1).epsilon(1e-10));

    // druga  komponenta je sin() pa očekujemo derivaciju cos()
    double der2 = Derivation::NDer4Partial(f, 1, point, 1e-6);
    REQUIRE(fDer(point, 1) == Approx(der2).epsilon(1e-9));
    REQUIRE(fDer(point, 1) != Approx(der2).epsilon(1e-10));

    // treća komponenta je exp() pa očekujemo derivaciju exp()
    double der3 = Derivation::NDer4Partial(f, 2, point, 1e-6);
    REQUIRE(fDer(point, 2) == Approx(der3).epsilon(1e-9));
    REQUIRE(fDer(point, 2) != Approx(der3).epsilon(1e-10));
}

TEST_CASE("Test_NDer6Partial", "[simple]") 
{
    VectorN<Real, 3> point{1.0, 1.0, 3.0 };

    ScalarFunction<3> f = TestData::ScalarFunctionsTestBed::_listFuncScalar3[0]._func;
    double (*fDer)(const VectorN<Real, 3> &, int ind) = TestData::ScalarFunctionsTestBed::_listFuncScalar3[0]._funcDerived;

    // prva komponenta je cos() pa očekujemo derivaciju -sin()
    double der1 = Derivation::NDer6Partial(f, 0, point, 1e-6);
    REQUIRE(fDer(point, 0) == Approx(der1).epsilon(1e-8));
    REQUIRE(fDer(point, 0) != Approx(der1).epsilon(1e-9));

    // druga  komponenta je sin() pa očekujemo derivaciju cos()
    double der2 = Derivation::NDer6Partial(f, 1, point, 1e-6);
    REQUIRE(fDer(point, 1) == Approx(der2).epsilon(1e-9));
    REQUIRE(fDer(point, 1) != Approx(der2).epsilon(1e-10));

    // treća komponenta je exp() pa očekujemo derivaciju exp()
    double der3 = Derivation::NDer6Partial(f, 2, point, 1e-6);
    REQUIRE(fDer(point, 2) == Approx(der3).epsilon(1e-9));
    REQUIRE(fDer(point, 2) != Approx(der3).epsilon(1e-10));
}

TEST_CASE("Test_NDer8Partial", "[simple]") 
{
    VectorN<Real, 3> point{1.0, 1.0, 3.0 };

    ScalarFunction<3> f = TestData::ScalarFunctionsTestBed::_listFuncScalar3[0]._func;
    double (*fDer)(const VectorN<Real, 3> &, int ind) = TestData::ScalarFunctionsTestBed::_listFuncScalar3[0]._funcDerived;

    // prva komponenta je cos() pa očekujemo derivaciju -sin()
    double der1 = Derivation::NDer8Partial(f, 0, point, 1e-6);
    REQUIRE(fDer(point, 0) == Approx(der1).epsilon(1e-8));
    REQUIRE(fDer(point, 0) != Approx(der1).epsilon(1e-9));

    // druga  komponenta je sin() pa očekujemo derivaciju cos()
    double der2 = Derivation::NDer8Partial(f, 1, point, 1e-6);
    REQUIRE(fDer(point, 1) == Approx(der2).epsilon(1e-9));
    REQUIRE(fDer(point, 1) != Approx(der2).epsilon(1e-10));

    // treća komponenta je exp() pa očekujemo derivaciju exp()
    double der3 = Derivation::NDer8Partial(f, 2, point, 1e-6);
    REQUIRE(fDer(point, 2) == Approx(der3).epsilon(1e-9));
    REQUIRE(fDer(point, 2) != Approx(der3).epsilon(1e-10));
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///////         Derivation of VECTOR functions
////////////////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Test_NDer1VecPartial", "[simple]") 
{
    VectorN<Real, 3> point{1.0, 1.0, 3.0 };

    VectorFunction<3> f = TestData::VectorFunctionsTestBed::_listFuncVector3[0]._func;
    VectorN<Real, 3> (*fDer)(const VectorN<Real, 3> &, int ind) = TestData::VectorFunctionsTestBed::_listFuncVector3[0]._funcDerived;

    VectorN<Real, 3> der1 = Derivation::NDer1PartialByAll(f, 0, point, 1e-6);
    REQUIRE(fDer(point, 0).IsEqual(der1, 1e-6));
    REQUIRE(!fDer(point, 0).IsEqual(der1, 1e-7));

    VectorN<Real, 3> der2 = Derivation::NDer1PartialByAll(f, 1, point, 1e-6);
    REQUIRE(fDer(point, 1).IsEqual(der2, 1e-6));
    REQUIRE(!fDer(point, 1).IsEqual(der2, 1e-7));

    // treća komponenta je exp() pa očekujemo derivaciju exp()
    VectorN<Real, 3> der3 = Derivation::NDer1PartialByAll(f, 2, point, 1e-6);
    REQUIRE(fDer(point, 2).IsEqual(der3, 1e-4));
    REQUIRE(!fDer(point, 2).IsEqual(der3, 1e-5));
}

TEST_CASE("Test_NDer2VecPartial", "[simple]") 
{
    VectorN<Real, 3> point{1.0, 1.0, 3.0 };

    VectorFunction<3> f = TestData::VectorFunctionsTestBed::_listFuncVector3[0]._func;
    VectorN<Real, 3> (*fDer)(const VectorN<Real, 3> &, int ind) = TestData::VectorFunctionsTestBed::_listFuncVector3[0]._funcDerived;

    VectorN<Real, 3> der1 = Derivation::NDer2PartialByAll(f, 0, point, 1e-6);
    REQUIRE(fDer(point, 0).IsEqual(der1, 1e-10));
    REQUIRE(!fDer(point, 0).IsEqual(der1, 1e-11));

    VectorN<Real, 3> der2 = Derivation::NDer2PartialByAll(f, 1, point, 1e-6);
    REQUIRE(fDer(point, 1).IsEqual(der2, 1e-10));
    REQUIRE(!fDer(point, 1).IsEqual(der2, 1e-11));

    // treća komponenta je exp() pa očekujemo derivaciju exp()
    VectorN<Real, 3> der3 = Derivation::NDer2PartialByAll(f, 2, point, 1e-6);
    REQUIRE(fDer(point, 2).IsEqual(der3, 1e-8));
    REQUIRE(!fDer(point, 2).IsEqual(der3, 1e-9));
}

TEST_CASE("Test_NDer4VecPartial", "[simple]") 
{
    VectorN<Real, 3> point{1.0, 1.0, 3.0 };

    VectorFunction<3> f = TestData::VectorFunctionsTestBed::_listFuncVector3[0]._func;
    VectorN<Real, 3> (*fDer)(const VectorN<Real, 3> &, int ind) = TestData::VectorFunctionsTestBed::_listFuncVector3[0]._funcDerived;

    VectorN<Real, 3> der1 = Derivation::NDer4PartialByAll(f, 0, point, 1e-6);
    REQUIRE(fDer(point, 0).IsEqual(der1, 1e-10));
    REQUIRE(!fDer(point, 0).IsEqual(der1, 1e-11));

    VectorN<Real, 3> der2 = Derivation::NDer4PartialByAll(f, 1, point, 1e-6);
    REQUIRE(fDer(point, 1).IsEqual(der2, 1e-10));
    REQUIRE(!fDer(point, 1).IsEqual(der2, 1e-11));

    // treća komponenta je exp() pa očekujemo derivaciju exp()
    VectorN<Real, 3> der3 = Derivation::NDer4PartialByAll(f, 2, point, 1e-6);
    REQUIRE(fDer(point, 2).IsEqual(der3, 1e-8));
    REQUIRE(!fDer(point, 2).IsEqual(der3, 1e-9));
}

TEST_CASE("Test_NDer6VecPartial", "[simple]") 
{
    VectorN<Real, 3> point{1.0, 1.0, 3.0 };

    VectorFunction<3> f = TestData::VectorFunctionsTestBed::_listFuncVector3[0]._func;
    VectorN<Real, 3> (*fDer)(const VectorN<Real, 3> &, int ind) = TestData::VectorFunctionsTestBed::_listFuncVector3[0]._funcDerived;

    VectorN<Real, 3> der1 = Derivation::NDer6PartialByAll(f, 0, point, 1e-6);
    REQUIRE(fDer(point, 0).IsEqual(der1, 1e-10));
    REQUIRE(!fDer(point, 0).IsEqual(der1, 1e-11));

    VectorN<Real, 3> der2 = Derivation::NDer6PartialByAll(f, 1, point, 1e-6);
    REQUIRE(fDer(point, 1).IsEqual(der2, 1e-10));
    REQUIRE(!fDer(point, 1).IsEqual(der2, 1e-11));

    // treća komponenta je exp() pa očekujemo derivaciju exp()
    VectorN<Real, 3> der3 = Derivation::NDer6PartialByAll(f, 2, point, 1e-6);
    REQUIRE(fDer(point, 2).IsEqual(der3, 1e-8));
    REQUIRE(!fDer(point, 2).IsEqual(der3, 1e-9));
}

TEST_CASE("Test_NDer8VecPartial", "[simple]") 
{
    VectorN<Real, 3> point{1.0, 1.0, 3.0 };

    VectorFunction<3> f = TestData::VectorFunctionsTestBed::_listFuncVector3[0]._func;
    VectorN<Real, 3> (*fDer)(const VectorN<Real, 3> &, int ind) = TestData::VectorFunctionsTestBed::_listFuncVector3[0]._funcDerived;

    VectorN<Real, 3> der1 = Derivation::NDer8PartialByAll(f, 0, point, 1e-6);
    REQUIRE(fDer(point, 0).IsEqual(der1, 1e-10));
    REQUIRE(!fDer(point, 0).IsEqual(der1, 1e-11));

    VectorN<Real, 3> der2 = Derivation::NDer8PartialByAll(f, 1, point, 1e-6);
    REQUIRE(fDer(point, 1).IsEqual(der2, 1e-10));
    REQUIRE(!fDer(point, 1).IsEqual(der2, 1e-11));

    // treća komponenta je exp() pa očekujemo derivaciju exp()
    VectorN<Real, 3> der3 = Derivation::NDer8PartialByAll(f, 2, point, 1e-6);
    REQUIRE(fDer(point, 2).IsEqual(der3, 1e-8));
    REQUIRE(!fDer(point, 2).IsEqual(der3, 1e-9));
}


TEST_CASE("Verify_vector_equation1")
{

}