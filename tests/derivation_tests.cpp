#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MMLBasicTypes.h"
#else
#include "basic_types/Function.h"
#include "algorithms/Derivation.h"
#endif

#include "functions_test_bed.h"

double TestSin(double x)
{
    return sin(x);
}

// TODO - definirati još par vsti funkcija - member, + razno

TEST_CASE("Test_Derivative", "[simple]") {
    std::function<double(double)> func = TestSin;

    MML::RealFunctionFromStdFunc g{func};
    //MML::RealFunction  f{TestSin};
    MML::RealFunction  f = MML::Tests::FunctionsTestBed::_listFunc[0]._func;


    // different step sizes for h!
    double der = MML::Derivation::NDer1(f, 1.0, 1e-4);
    REQUIRE(cos(1.0) == Approx(der).epsilon(1e-3));
    REQUIRE(cos(1.0) == Approx(der).epsilon(1e-4));
    REQUIRE(cos(1.0) != Approx(der).epsilon(1e-5));
    REQUIRE(cos(1.0) != Approx(der).epsilon(1e-6));
    REQUIRE(cos(1.0) != Approx(der).epsilon(1e-7));

    der = MML::Derivation::NDer1(f, 1.0, 1e-5);
    REQUIRE(cos(1.0) == Approx(der).epsilon(1e-4));
    REQUIRE(cos(1.0) == Approx(der).epsilon(1e-5));
    REQUIRE(cos(1.0) != Approx(der).epsilon(1e-6));
    REQUIRE(cos(1.0) != Approx(der).epsilon(1e-7));

    der = MML::Derivation::NDer1(f, 1.0, 1e-6);
    REQUIRE(cos(1.0) == Approx(der).epsilon(1e-5));
    REQUIRE(cos(1.0) == Approx(der).epsilon(1e-6));
    REQUIRE(cos(1.0) != Approx(der).epsilon(1e-7));
    REQUIRE(cos(1.0) != Approx(der).epsilon(1e-8));
    
    der = MML::Derivation::NDer1(f, 1.0, 1e-7);
    REQUIRE(cos(1.0) == Approx(der).epsilon(1e-6));
    REQUIRE(cos(1.0) == Approx(der).epsilon(1e-7));
    REQUIRE(cos(1.0) != Approx(der).epsilon(1e-8));
    REQUIRE(cos(1.0) != Approx(der).epsilon(1e-9));
}

double TestScalarFunction(MML::VectorN<3> &x)
{
    return cos(x[0]) + sin(x[1]) + exp(x[2]);
}

TEST_CASE("Test_Partial_Derivative", "[simple]") {
    std::function<double(MML::VectorN<3> &)> func = TestScalarFunction;

    MML::ScalarFunctionFromStdFunc<3> f{func};

    MML::VectorN<3> point{1.0, 1.0, 1.0 };
    double c = func(point);

    double d = f(point);

    // prva komponenta je cos() pa očekujemo derivaciju -sin()
    double der1 = MML::Derivation::NDer1Partial(f, 0, point, 1e-10);
    REQUIRE(-sin(1.0) == Approx(der1).epsilon(1e-8));

    // druga  komponenta je sin() pa očekujemo derivaciju cos()
    double der2 = MML::Derivation::NDer1Partial(f, 1, point, 1e-10);
    REQUIRE(cos(1.0) == Approx(der2).epsilon(1e-10));

    // treća komponenta je exp() pa očekujemo derivaciju exp()
    double der3 = MML::Derivation::NDer1Partial(f, 2, point, 1e-10);
    REQUIRE(exp(1.0) == Approx(der3).epsilon(1e-10));

}