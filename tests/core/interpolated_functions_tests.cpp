#if !defined __MML_INTERPOLATED_FUNCTIONS_H
#define __MML_INTERPOLATED_FUNCTIONS_H

#include "../catch/catch.hpp"

#include <vector>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/InterpolatedFunction.h"

#include "algorithms/FunctionAnalyzers.h"
#endif

using namespace MML;

namespace Tests
{
double test_func(const double x)
{
    const double eps = 1.0;
    return x*exp(-x)/(POW2(x-1.0)+eps*eps);
}

void CreateInterpolatedValues(RealFunction f, Real x1, Real x2, int numPnt, Vector<Real> &outX, Vector<Real> &outY)
{
    outX.Resize(numPnt);
    outY.Resize(numPnt);

    for (int i=0;i<numPnt;i++) {
        outX[i] = x1 + i * (x2 - x1) / (numPnt - 1);
        outY[i] = f(outX[i]);
    } 
}

// TODO - HIGH, EMPTY!!! - implement tests for interpolation methods and verify precision
TEST_CASE("Test_LinearInterp_sin_func", "[simple]") {
    RealFunction f{[](Real x) { return sin(x); } };

    Vector<Real> vec_x, vec_y;
    CreateInterpolatedValues(f, 0.0, 2.0*Constants::PI, 10, vec_x, vec_y);

    LinearInterpRealFunc    linear_interp(vec_x, vec_y);

    RealFunctionComparer      comparer(f, linear_interp);

    double totalAbsDiff = comparer.getAbsDiffSum(0.0, 2.0*Constants::PI, 100);
    double avgAbsDiff = comparer.getAbsDiffAvg(0.0, 2.0*Constants::PI, 100);
    double maxAbsDiff = comparer.getAbsDiffMax(0.0, 2.0*Constants::PI, 100);

    double totalRelDiff = comparer.getRelDiffSum(0.0, 2.0*Constants::PI, 100);
    double avgRelDiff = comparer.getRelDiffAvg(0.0, 2.0*Constants::PI, 100);
    double maxRelDiff = comparer.getRelDiffMax(0.0, 2.0*Constants::PI, 100);

    // REQUIRE(totalAbsDiff == Approx(0.0).margin(0.0001));
    // REQUIRE(avgAbsDiff == Approx(0.0).margin(0.0001));
    // REQUIRE(maxAbsDiff == Approx(0.0).margin(0.0001));

    // proci kroz 100 tocaka, i usporediti apsolutne vrijednosti
    // avg razlike, max razlike, min razlike
    // relativna greška u odnosu na točnu vrijednost
}


} // end namespace
#endif