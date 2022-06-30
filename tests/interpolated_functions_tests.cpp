#if !defined __MML_INTERPOLATED_FUNCTIONS_H
#define __MML_INTERPOLATED_FUNCTIONS_H

#include "../catch/catch.hpp"

#include <vector>

#ifdef MML_USE_SINGLE_HEADER
#include "MMLBasicTypes.h"
#else
#include "basic_types/InterpolatedFunction.h"
#endif

namespace MML::Tests
{
    static double Func1(double x) { return sin(x); }
    static double Func1_derived(double x) { return cos(x); }
    static double Func1_integrated(double x) { return -cos(x); }
}

TEST_CASE("Test_TabulatedValues1D", "[simple]") {
    std::vector<double> values(101);

    for(int i=0; i<=100; i++ )
        values[i] = sin(i / 10.0);

    MML::TabulatedValues1DEqualSpacing tabValues{0.0, 10.0, values};

    REQUIRE(sin(0.0) == tabValues.Value(0));
    REQUIRE(sin(5.0) == tabValues.Value(50));
    REQUIRE(sin(10.0) == tabValues.Value(100));
}

#endif