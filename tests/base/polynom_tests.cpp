#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Polynom.h"
#endif

using namespace MML;

namespace MML::Tests::Base::PolynomTests
{
  TEST_CASE("Polynom coef.order")
  {
    RealPolynom poly4({1.0, 0.0, 0.0, 0.0, 5.0});   // polynom 1 + 5 * x^4
    
    REQUIRE(poly4.GetDegree() == 4);
    REQUIRE(poly4[0] == 1.0);
    REQUIRE(poly4[4] == 5.0);

    REQUIRE(poly4(2.0) == 5 * 2*2*2*2 + 1);
  }

  TEST_CASE("Polynom Init, GetDegree, operator[]", "[simple]")
  {
    RealPolynom poly_empty;
    REQUIRE(poly_empty.GetDegree() == -1);

    RealPolynom poly3(3);
    REQUIRE(poly3.GetDegree() == 3);

    RealPolynom poly4({1.0, 2.0, 3.0, 4.0, 5.0});
    REQUIRE(poly4.GetDegree() == 4);
    REQUIRE(poly4[0] == 1.0);
    REQUIRE(poly4[1] == 2.0);
    REQUIRE(poly4[2] == 3.0);
    REQUIRE(poly4[3] == 4.0);
    REQUIRE(poly4[4] == 5.0);

    std::vector<Real> coefs{-1.0, -2.0, -3.0, -4.0, -5.0};
    RealPolynom poly4b(coefs);
    REQUIRE(poly4b.GetDegree() == 4);
    REQUIRE(poly4b[0] == -1.0);
    REQUIRE(poly4b[1] == -2.0);
    REQUIRE(poly4b[2] == -3.0);
    REQUIRE(poly4b[3] == -4.0);
    REQUIRE(poly4b[4] == -5.0);

    // negativan argument
    // prazan vektor
    // TODO - polynom test - implement for Complex and Matrix
  }

  TEST_CASE("Polynom SetDegree", "[simple]")
  {
    RealPolynom poly;
    REQUIRE(poly.GetDegree() == -1);

    poly.SetDegree(3);
    REQUIRE(poly.GetDegree() == 3);
  }

  TEST_CASE("Polynom operator()", "[simple]")
  {
    RealPolynom pol_constant({3.5});

    REQUIRE(pol_constant(0.0) == 3.5);
    REQUIRE(pol_constant(1.0) == Approx(3.5));
    REQUIRE(pol_constant(-200.0) == Approx(3.5));
    REQUIRE(pol_constant(1.0e-10) == Approx(3.5));
    REQUIRE(pol_constant(1.0e+10) == Approx(3.5));

    RealPolynom pol_lin({3.5, 1.0});

    REQUIRE(pol_lin(0.0) == Approx(3.5));
    REQUIRE(pol_lin(1.0) == 4.5);
    REQUIRE(pol_lin(-100.0) == Approx(-96.5));
  }

  TEST_CASE("Polynom operator==", "[simple]")
  {
    RealPolynom poly3(3);

    RealPolynom poly4({1.0, 2.0, 3.0, 4.0, 5.0});
    RealPolynom poly4a({1.0, 2.0, 3.0, 4.0, 5.0});

    std::vector<Real> coefs{-1.0, -2.0, -3.0, -4.0, -5.0};
    RealPolynom poly4b(coefs);

    REQUIRE(poly4 == poly4a);
    REQUIRE(poly4 != poly4b);
    REQUIRE(poly4 == poly4b * (-1));
  }

  TEST_CASE("Polynom operator+", "[simple]")
  {
    RealPolynom poly1({3.5});
    RealPolynom poly2({1.0, 0.5, 2.0});

    RealPolynom poly3 = poly1 + poly2;

    REQUIRE(poly3.GetDegree() == 2);
    REQUIRE(poly3[0] == 4.5);
    REQUIRE(poly3[1] == 0.5);
    REQUIRE(poly3[2] == 2.0);
  }
  TEST_CASE("Polynom operator-", "[simple]")
  {
    RealPolynom poly1({3.5});
    RealPolynom poly2({1.0, 0.5, 2.0});

    RealPolynom poly3 = poly1 - poly2;

    REQUIRE(poly3.GetDegree() == 2);
    REQUIRE(poly3[0] == 2.5);
    REQUIRE(poly3[1] == -0.5);
    REQUIRE(poly3[2] == -2.0);
  }
  TEST_CASE("Polynom operator*", "[simple]")
  {
    RealPolynom poly1({2, 1, 3});
    RealPolynom poly2({1, 2});

    RealPolynom poly3 = poly1 * poly2;
    REQUIRE(poly3.GetDegree() == 3);
    REQUIRE(poly3[0] == 2);
    REQUIRE(poly3[1] == 5);
    REQUIRE(poly3[2] == 5);
    REQUIRE(poly3[3] == 6);
  }
  TEST_CASE("Polynom operator*(Poly, CoefType)", "[simple]")
  {
    RealPolynom poly({1, 2, 3});

    RealPolynom poly2 = poly * 3.0;
    REQUIRE(poly2.GetDegree() == 2);
    REQUIRE(poly2[0] == 3.0);
    REQUIRE(poly2[1] == 6.0);
    REQUIRE(poly2[2] == 9.0);
  }
  TEST_CASE("Test_Polynom_operator*(CoefType, Poly)", "[simple]")
  {
    RealPolynom poly({1, 2, 3});

    RealPolynom poly2 = 3.0 * poly;
    REQUIRE(poly2.GetDegree() == 2);
    REQUIRE(poly2[0] == 3.0);
    REQUIRE(poly2[1] == 6.0);
    REQUIRE(poly2[2] == 9.0);
  }
  TEST_CASE("Polynom operator/(Poly, CoefType)", "[simple]")
  {
    RealPolynom poly({2, 4, 6});

    RealPolynom poly2 = poly / 2.0;
    REQUIRE(poly2.GetDegree() == 2);
    REQUIRE(poly2[0] == 1.0);
    REQUIRE(poly2[1] == 2.0);
    REQUIRE(poly2[2] == 3.0);
  }
  
  TEST_CASE("Polynom to_string", "[simple]")
  {
    RealPolynom poly({1, 2, 3});

    REQUIRE(poly.to_string(10, 5) ==   "   3.00000*x^2 + 2.00000*x + 1.00000");    
  }

  TEST_CASE("Polynom poldiv1", "[simple]")
  {
    RealPolynom p1({-1.0, 1.0}); // x - 1
    RealPolynom p2({-2.0, 1.0}); // x - 2
    RealPolynom p3({-3.0, 1.0}); // x - 3

    RealPolynom p = p1 * p2 * p3;

    REQUIRE(p.GetDegree() == 3);
    REQUIRE(p[0] == -6.0);
    REQUIRE(p[1] == 11.0);
    REQUIRE(p[2] == -6.0);
    REQUIRE(p[3] == 1.0);

    RealPolynom res, rem;

    RealPolynom::poldiv(p, p1, res, rem);
    REQUIRE(p.IsEqual(res * p1 + rem));
    REQUIRE(res == (p2 * p3));

    RealPolynom::poldiv(p, p2, res, rem);
    REQUIRE(p.IsEqual(res * p2 + rem));
    REQUIRE(res == (p1 * p3));

    RealPolynom::poldiv(p, p3, res, rem);
    REQUIRE(p.IsEqual(res * p3 + rem));
    REQUIRE(res == (p1 * p2));
  }
  
  TEST_CASE("Polynom poldiv2", "[simple]")
  {
    RealPolynom p_1({2.0, 0.0, 1.0});
    RealPolynom p_2({-4.0, 1.0});
    RealPolynom p_quot({4.0, 1.0});
    RealPolynom p_rem({18.0});

    RealPolynom res, rem;

    RealPolynom::poldiv(p_1, p_2, res, rem);

    REQUIRE(p_1 == res * p_2 + rem);
    REQUIRE(res == p_quot);
    REQUIRE(rem == p_rem);

    RealPolynom r_1({-1.0, -5.0, 2.0});
    RealPolynom r_2({-3.0, 1.0});
    RealPolynom r_quot({1.0, 2.0});
    RealPolynom r_rem({2.0});

    RealPolynom::poldiv(r_1, r_2, res, rem);

    REQUIRE(r_1 == res * r_2 + rem);
    REQUIRE(res == r_quot);
    REQUIRE(rem == r_rem);
  }
  
  TEST_CASE("Polynom poldiv3", "[simple]")
  {
    RealPolynom p3_u({-1.0, -5.0 });
    RealPolynom p3_v({-3.0, 1.0, 2.0});

    RealPolynom res, rem;

    RealPolynom::poldiv(p3_u, p3_v, res, rem);

    REQUIRE(res.GetDegree() == -1);
    REQUIRE(res.IsNullPolynom() == true);
    REQUIRE(rem == p3_u);
  }
}