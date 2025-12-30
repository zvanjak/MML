#include <catch2/catch_all.hpp>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/Function.h"
#include "core/Derivation.h"
#endif

#include "../test_data/real_functions_test_bed.h"
#include "../test_data/scalar_functions_test_bed.h"
#include "../test_data/vector_functions_test_bed.h"

using namespace MML;

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;


////////////////////////////////////////////////////////////////////////////////////////////////////
///////         Derivation of REAL functions
////////////////////////////////////////////////////////////////////////////////////////////////////

// TODO - Test_..._derivation_precision_default_step_size - za 5 funkcija, i prosiriti interval

TEST_CASE("Test_1st_derivation_precision_sin_func_default_step_size", "[simple]") 
{
    RealFunction  f_sin = TestBeds::RealFunctionsTestBed::getFunc(0)._func;
    RealFunction  f_der = TestBeds::RealFunctionsTestBed::getFunc(0)._funcDerived;

    RealFunction  f_cos = TestBeds::RealFunctionsTestBed::getFunc(1)._func;

    Vector<Real> x_val{0.0, 0.5, 1.0, 2.0, 3.0, 5.0 };

#if defined(__clang__)
    Vector<Real> nder1_prec{1e-16,  1e-8,  1e-7,  1e-7,  1e-8,  1e-7 };
    Vector<Real> nder2_prec{1e-10, 1e-11, 1e-10, 1e-10, 1e-10, 1e-10 };
    Vector<Real> nder4_prec{1e-13, 1e-13, 1e-12, 1e-13, 1e-13, 1e-12 };
    Vector<Real> nder6_prec{1e-15, 1e-13, 1e-14, 1e-13, 1e-13, 1e-13 };
    Vector<Real> nder8_prec{1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14 };
#elif defined(__GNUC__) || defined(__GNUG__)
    Vector<Real> nder1_prec{1e-16,  1e-8,  1e-7,  1e-7,  1e-8,  1e-7 };
    Vector<Real> nder2_prec{1e-10, 1e-11, 1e-10, 1e-10, 1e-10, 1e-10 };
    Vector<Real> nder4_prec{1e-13, 1e-13, 1e-12, 1e-13, 1e-13, 1e-12 };
    Vector<Real> nder6_prec{1e-15, 1e-13, 1e-14, 1e-13, 1e-13, 1e-13 };
    Vector<Real> nder8_prec{1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14 };
#elif defined(_MSC_VER)
    Vector<Real> nder1_prec{1e-16,  1e-8,  1e-7,  1e-7,  1e-8,  1e-7 };
    Vector<Real> nder2_prec{1e-10, 1e-11, 1e-10, 1e-10, 1e-10, 1e-10 };
    Vector<Real> nder4_prec{1e-13, 1e-13, 1e-12, 1e-13, 1e-13, 1e-12 };
    Vector<Real> nder6_prec{1e-15, 1e-13, 1e-14, 1e-13, 1e-13, 1e-13 };
    Vector<Real> nder8_prec{1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14 };
#endif

    for( int i=0; i<x_val.size(); i++ )
    {
        double x = x_val[i];

        double num_der = Derivation::NDer1(f_sin, x);
        REQUIRE_THAT(f_der(x) , WithinRel(num_der, nder1_prec[i]));
        REQUIRE_THAT(f_der(x) , !WithinRel(num_der, nder1_prec[i]*0.1));

        num_der = Derivation::NDer2(f_sin, x);
        REQUIRE_THAT(f_der(x) , WithinRel(num_der, nder2_prec[i]));
        REQUIRE_THAT(f_der(x) , !WithinRel(num_der, nder2_prec[i]*0.1));

        num_der = Derivation::NDer4(f_sin, x);
        REQUIRE_THAT(f_der(x) , WithinRel(num_der, nder4_prec[i]));
        REQUIRE_THAT(f_der(x) , !WithinRel(num_der, nder4_prec[i]*0.1));  

        num_der = Derivation::NDer6(f_sin, x);
        REQUIRE_THAT(f_der(x) , WithinRel(num_der, nder6_prec[i]));
        REQUIRE_THAT(f_der(x) , !WithinRel(num_der, nder6_prec[i]*0.1));              

        num_der = Derivation::NDer8(f_sin, x);
        REQUIRE_THAT(f_der(x) , WithinRel(num_der, nder8_prec[i]));
        REQUIRE_THAT(f_der(x) , !WithinRel(num_der, nder8_prec[i]*0.1));              
    }   
}

TEST_CASE("Test_2nd_derivation_precision_sin_func_default_step_size", "[simple]") 
{
}

TEST_CASE("Test_3rd_derivation_precision_sin_func_default_step_size", "[simple]") 
{
}
