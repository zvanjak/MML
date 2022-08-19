#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/Function.h"
#include "algorithms/Derivation.h"
#endif

#include "test_data/functions_test_bed.h"

TEST_CASE("Test_sin_derivation_precision", "[simple]") 
{
    MML::RealFunction  f = MML::Tests::FunctionsTestBed::_listFuncReal[0]._func;
    MML::RealFunction  fDer = MML::Tests::FunctionsTestBed::_listFuncReal[0]._funcDerived;

    double x = 1.0;

    double der = MML::Derivation::NDer1(f, x);
    REQUIRE(fDer(x) == Approx(der).epsilon(1e-7));
    REQUIRE(fDer(x) != Approx(der).epsilon(1e-8));

    der = MML::Derivation::NDer2(f, x);
    REQUIRE(fDer(x) == Approx(der).epsilon(1e-10));
    REQUIRE(fDer(x) != Approx(der).epsilon(1e-11));

    der = MML::Derivation::NDer4(f, x);
    REQUIRE(fDer(x) == Approx(der).epsilon(1e-12));
    REQUIRE(fDer(x) != Approx(der).epsilon(1e-13));
    
    der = MML::Derivation::NDer6(f, x);
    REQUIRE(fDer(x) == Approx(der).epsilon(1e-14));
    REQUIRE(fDer(x) != Approx(der).epsilon(1e-15));
    
    der = MML::Derivation::NDer8(f, x);
    REQUIRE(fDer(x) == Approx(der).epsilon(1e-14));
    REQUIRE(fDer(x) != Approx(der).epsilon(1e-15));    
}

TEST_CASE("Test_NDer1", "[simple]") 
{
    MML::RealFunction  f = MML::Tests::FunctionsTestBed::_listFuncReal[0]._func;
    MML::RealFunction  fDer = MML::Tests::FunctionsTestBed::_listFuncReal[0]._funcDerived;

    double x = 1.0;

    double der = MML::Derivation::NDer1(f, x, 1e-4);
    REQUIRE(fDer(x) == Approx(der).epsilon(1e-4));
    REQUIRE(fDer(x) != Approx(der).epsilon(1e-5));

    der = MML::Derivation::NDer1(f, x, 1e-5);
    REQUIRE(fDer(x) == Approx(der).epsilon(1e-5));
    REQUIRE(fDer(x) != Approx(der).epsilon(1e-6));

    der = MML::Derivation::NDer1(f, x, 1e-6);
    REQUIRE(fDer(x) == Approx(der).epsilon(1e-6));
    REQUIRE(fDer(x) != Approx(der).epsilon(1e-7));
    
    der = MML::Derivation::NDer1(f, x, 1e-7);
    REQUIRE(fDer(x) == Approx(der).epsilon(1e-7));
    REQUIRE(fDer(x) != Approx(der).epsilon(1e-8));
    
    der = MML::Derivation::NDer1(f, x, 1e-8);
    REQUIRE(fDer(x) == Approx(der).epsilon(1e-8));
    REQUIRE(fDer(x) != Approx(der).epsilon(1e-9));    
}

TEST_CASE("Test_NDer2", "[simple]") 
{
    MML::RealFunction  f = MML::Tests::FunctionsTestBed::_listFuncReal[0]._func;
    MML::RealFunction  fDer = MML::Tests::FunctionsTestBed::_listFuncReal[0]._funcDerived;

    // different step sizes for h!
    double der = MML::Derivation::NDer2(f, 1.0, 1e-4);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-8));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-9));

    der = MML::Derivation::NDer2(f, 1.0, 1e-5);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-10));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-11));

    der = MML::Derivation::NDer2(f, 1.0, 1e-6);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-10));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-11));
    
    der = MML::Derivation::NDer2(f, 1.0, 1e-7);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-9));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-10));

    der = MML::Derivation::NDer2(f, 1.0, 1e-8);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-8));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-9));
}

TEST_CASE("Test_NDer4", "[simple]") 
{
    MML::RealFunction  f = MML::Tests::FunctionsTestBed::_listFuncReal[0]._func;
    MML::RealFunction  fDer = MML::Tests::FunctionsTestBed::_listFuncReal[0]._funcDerived;

    // different step sizes for h!
    double der = MML::Derivation::NDer4(f, 1.0, 1e-4);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-8));
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-11));

    der = MML::Derivation::NDer4(f, 1.0, 1e-5);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-11));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-12));

    der = MML::Derivation::NDer4(f, 1.0, 1e-6);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-10));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-11));
    
    der = MML::Derivation::NDer4(f, 1.0, 1e-7);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-9));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-10));

    der = MML::Derivation::NDer4(f, 1.0, 1e-8);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-8));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-9));
}

TEST_CASE("Test_NDer6", "[simple]") 
{
    MML::RealFunction  f = MML::Tests::FunctionsTestBed::_listFuncReal[0]._func;
    MML::RealFunction  fDer = MML::Tests::FunctionsTestBed::_listFuncReal[0]._funcDerived;

    // different step sizes for h!
    double der = MML::Derivation::NDer6(f, 1.0, 1e-4);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-8));
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-11));

    der = MML::Derivation::NDer6(f, 1.0, 1e-5);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-11));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-12));

    der = MML::Derivation::NDer6(f, 1.0, 1e-6);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-10));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-11));
    
    der = MML::Derivation::NDer6(f, 1.0, 1e-7);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-9));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-10));

    der = MML::Derivation::NDer6(f, 1.0, 1e-8);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-8));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-9));
}

TEST_CASE("Test_NDer8", "[simple]") 
{
    MML::RealFunction  f = MML::Tests::FunctionsTestBed::_listFuncReal[0]._func;
    MML::RealFunction  fDer = MML::Tests::FunctionsTestBed::_listFuncReal[0]._funcDerived;

    // different step sizes for h!
    double der = MML::Derivation::NDer8(f, 1.0, 1e-4);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-8));
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-11));

    der = MML::Derivation::NDer8(f, 1.0, 1e-5);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-11));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-12));

    der = MML::Derivation::NDer8(f, 1.0, 1e-6);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-10));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-11));
    
    der = MML::Derivation::NDer8(f, 1.0, 1e-7);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-9));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-10));

    der = MML::Derivation::NDer8(f, 1.0, 1e-8);
    REQUIRE(fDer(1.0) == Approx(der).epsilon(1e-8));
    REQUIRE(fDer(1.0) != Approx(der).epsilon(1e-9));
}

TEST_CASE("Test_NDer1Partial", "[simple]") 
{
    MML::VectorN<3> point{1.0, 1.0, 3.0 };

    MML::ScalarFunctionFromFuncPtr<3> f = MML::Tests::FunctionsTestBed::_listFuncScalar[0]._func;
    double (*fDer)(const MML::VectorN<3> &, int ind) = MML::Tests::FunctionsTestBed::_listFuncScalar[0]._funcDerived;

    // prva komponenta je cos() pa očekujemo derivaciju -sin()
    double der1 = MML::Derivation::NDer1Partial(f, 0, point, 1e-6);
    REQUIRE(fDer(point, 0) == Approx(der1).epsilon(1e-6));
    REQUIRE(fDer(point, 0) != Approx(der1).epsilon(1e-7));

    // druga  komponenta je sin() pa očekujemo derivaciju cos()
    double der2 = MML::Derivation::NDer1Partial(f, 1, point, 1e-6);
    REQUIRE(fDer(point, 1) == Approx(der2).epsilon(1e-6));
    REQUIRE(fDer(point, 1) != Approx(der2).epsilon(1e-7));

    // treća komponenta je exp() pa očekujemo derivaciju exp()
    double der3 = MML::Derivation::NDer1Partial(f, 2, point, 1e-6);
    REQUIRE(fDer(point, 2) == Approx(der3).epsilon(1e-6));
    REQUIRE(fDer(point, 2) != Approx(der3).epsilon(1e-7));
}

TEST_CASE("Test_NDer2Partial", "[simple]") 
{
    MML::VectorN<3> point{1.0, 1.0, 3.0 };

    MML::ScalarFunctionFromFuncPtr<3> f = MML::Tests::FunctionsTestBed::_listFuncScalar[0]._func;
    double (*fDer)(const MML::VectorN<3> &, int ind) = MML::Tests::FunctionsTestBed::_listFuncScalar[0]._funcDerived;

    // prva komponenta je cos() pa očekujemo derivaciju -sin()
    double der1 = MML::Derivation::NDer2Partial(f, 0, point, 1e-6);
    REQUIRE(fDer(point, 0) == Approx(der1).epsilon(1e-9));
    REQUIRE(fDer(point, 0) != Approx(der1).epsilon(1e-10));

    // druga  komponenta je sin() pa očekujemo derivaciju cos()
    double der2 = MML::Derivation::NDer2Partial(f, 1, point, 1e-6);
    REQUIRE(fDer(point, 1) == Approx(der2).epsilon(1e-9));
    REQUIRE(fDer(point, 1) != Approx(der2).epsilon(1e-10));

    // treća komponenta je exp() pa očekujemo derivaciju exp()
    double der3 = MML::Derivation::NDer2Partial(f, 2, point, 1e-6);
    REQUIRE(fDer(point, 2) == Approx(der3).epsilon(1e-9));
    REQUIRE(fDer(point, 2) != Approx(der3).epsilon(1e-10));
}

TEST_CASE("Test_NDer4Partial", "[simple]") 
{
    MML::VectorN<3> point{1.0, 1.0, 3.0 };

    MML::ScalarFunctionFromFuncPtr<3> f = MML::Tests::FunctionsTestBed::_listFuncScalar[0]._func;
    double (*fDer)(const MML::VectorN<3> &, int ind) = MML::Tests::FunctionsTestBed::_listFuncScalar[0]._funcDerived;

    // prva komponenta je cos() pa očekujemo derivaciju -sin()
    double der1 = MML::Derivation::NDer4Partial(f, 0, point, 1e-6);
    REQUIRE(fDer(point, 0) == Approx(der1).epsilon(1e-9));
    REQUIRE(fDer(point, 0) != Approx(der1).epsilon(1e-10));

    // druga  komponenta je sin() pa očekujemo derivaciju cos()
    double der2 = MML::Derivation::NDer4Partial(f, 1, point, 1e-6);
    REQUIRE(fDer(point, 1) == Approx(der2).epsilon(1e-9));
    REQUIRE(fDer(point, 1) != Approx(der2).epsilon(1e-10));

    // treća komponenta je exp() pa očekujemo derivaciju exp()
    double der3 = MML::Derivation::NDer4Partial(f, 2, point, 1e-6);
    REQUIRE(fDer(point, 2) == Approx(der3).epsilon(1e-9));
    REQUIRE(fDer(point, 2) != Approx(der3).epsilon(1e-10));
}

TEST_CASE("Test_NDer6Partial", "[simple]") 
{
    MML::VectorN<3> point{1.0, 1.0, 3.0 };

    MML::ScalarFunctionFromFuncPtr<3> f = MML::Tests::FunctionsTestBed::_listFuncScalar[0]._func;
    double (*fDer)(const MML::VectorN<3> &, int ind) = MML::Tests::FunctionsTestBed::_listFuncScalar[0]._funcDerived;

    // prva komponenta je cos() pa očekujemo derivaciju -sin()
    double der1 = MML::Derivation::NDer6Partial(f, 0, point, 1e-6);
    REQUIRE(fDer(point, 0) == Approx(der1).epsilon(1e-8));
    REQUIRE(fDer(point, 0) != Approx(der1).epsilon(1e-9));

    // druga  komponenta je sin() pa očekujemo derivaciju cos()
    double der2 = MML::Derivation::NDer6Partial(f, 1, point, 1e-6);
    REQUIRE(fDer(point, 1) == Approx(der2).epsilon(1e-9));
    REQUIRE(fDer(point, 1) != Approx(der2).epsilon(1e-10));

    // treća komponenta je exp() pa očekujemo derivaciju exp()
    double der3 = MML::Derivation::NDer6Partial(f, 2, point, 1e-6);
    REQUIRE(fDer(point, 2) == Approx(der3).epsilon(1e-9));
    REQUIRE(fDer(point, 2) != Approx(der3).epsilon(1e-10));
}

TEST_CASE("Test_NDer8Partial", "[simple]") 
{
    MML::VectorN<3> point{1.0, 1.0, 3.0 };

    MML::ScalarFunctionFromFuncPtr<3> f = MML::Tests::FunctionsTestBed::_listFuncScalar[0]._func;
    double (*fDer)(const MML::VectorN<3> &, int ind) = MML::Tests::FunctionsTestBed::_listFuncScalar[0]._funcDerived;

    // prva komponenta je cos() pa očekujemo derivaciju -sin()
    double der1 = MML::Derivation::NDer8Partial(f, 0, point, 1e-6);
    REQUIRE(fDer(point, 0) == Approx(der1).epsilon(1e-8));
    REQUIRE(fDer(point, 0) != Approx(der1).epsilon(1e-9));

    // druga  komponenta je sin() pa očekujemo derivaciju cos()
    double der2 = MML::Derivation::NDer8Partial(f, 1, point, 1e-6);
    REQUIRE(fDer(point, 1) == Approx(der2).epsilon(1e-9));
    REQUIRE(fDer(point, 1) != Approx(der2).epsilon(1e-10));

    // treća komponenta je exp() pa očekujemo derivaciju exp()
    double der3 = MML::Derivation::NDer8Partial(f, 2, point, 1e-6);
    REQUIRE(fDer(point, 2) == Approx(der3).epsilon(1e-9));
    REQUIRE(fDer(point, 2) != Approx(der3).epsilon(1e-10));
}