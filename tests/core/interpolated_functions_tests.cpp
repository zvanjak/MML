#if !defined __MML_INTERPOLATED_FUNCTIONS_H
#define __MML_INTERPOLATED_FUNCTIONS_H

#include "../catch/catch.hpp"

#include <vector>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/InterpolatedFunction.h"
#endif

using namespace MML;

namespace Tests
{
    static double Func1(double x) { return sin(x); }
    static double Func1_derived(double x) { return cos(x); }
    static double Func1_integrated(double x) { return -cos(x); }
}

// TODO - HIGH, EMPTY!!! - implement tests for interpolation methods and verify precision

// TEST_CASE("Test_TabulatedValues1D", "[simple]") {
//     std::vector<double> values(101);

//     for(int i=0; i<=100; i++ )
//         values[i] = sin(i / 10.0);

//     TabulatedValues1DEqualSpacing tabValues{0.0, 10.0, values};

//     REQUIRE(sin(0.0) == tabValues.Value(0));
//     REQUIRE(sin(5.0) == tabValues.Value(50));
//     REQUIRE(sin(10.0) == tabValues.Value(100));
// }

#endif