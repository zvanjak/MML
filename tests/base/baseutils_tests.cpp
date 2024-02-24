#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/BaseUtils.h"
#endif

using namespace MML;

// TODO 0.9 - LOW - finish these tests
TEST_CASE("Test Abs")
{
    Complex a(1, -3);

    double b = Abs(a);

    REQUIRE(b == sqrt(10));
}

TEST_CASE("Test_Complex_AreEqual", "[simple]") {
}

TEST_CASE("Test_Vector<Complex>_AreEqual", "[simple]") {
}

TEST_CASE("Test_VectorProjectionParallelTo", "[simple]") {
}

TEST_CASE("Test_VectorProjectionPerpendicularTo", "[simple]") {
}

TEST_CASE("Test_OuterProduct", "[simple]") {
}

TEST_CASE("Test_Utils_AddVec", "[simple]") {
    Vector<Complex> a({Complex(1,0), Complex(0,1)});
    Vector<Real>    b({1,2});
    
    auto c = Utils::AddVec(a,b);

    REQUIRE(2 == c.size());
    REQUIRE(Complex(2,0) == c[0]);
    REQUIRE(Complex(2,1) == c[1]);

    auto d = Utils::AddVec(b, a);
    REQUIRE(2 == d.size());
    REQUIRE(Complex(2,0) == d[0]);
    REQUIRE(Complex(2,1) == d[1]);
}

TEST_CASE("Test_Utils_SubVec", "[simple]") {
    Vector<Complex> a({Complex(1,0), Complex(0,1)});
    Vector<Real>    b({1,2});
    
    auto c = Utils::SubVec(a,b);

    REQUIRE(2 == c.size());
    REQUIRE(Complex(0,0) == c[0]);
    REQUIRE(Complex(-2,1) == c[1]);

    auto d = Utils::SubVec(b, a);
    REQUIRE(2 == d.size());
    REQUIRE(Complex(0,0) == d[0]);
    REQUIRE(Complex(2,-1) == d[1]);

}
