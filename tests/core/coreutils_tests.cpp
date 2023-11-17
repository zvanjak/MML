#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/CoreUtils.h"
#endif

using namespace MML;

// TODO - LOW, EMPTY!!! finish these tests
TEST_CASE("Test Abs")
{
    Complex a(1, -3);

    double b = Abs(a);

    REQUIRE(b == sqrt(10));
}

TEST_CASE("Test_Utils_Vector_operations", "[simple]") {
    Vector<Complex> a({Complex(1,0), Complex(0,1)});
    Vector<Real> b({1,2});
    
    auto c = Utils::AddVec(a,b);

    REQUIRE(2 == c.size());
    REQUIRE(Complex(2,0) == c[0]);
    REQUIRE(Complex(2,1) == c[1]);
}

TEST_CASE("Test_Utils_Matrix_operations", "[simple]") {
}