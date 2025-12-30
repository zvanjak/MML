#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Polynom.h"
#endif

using namespace MML;
using namespace MML::Testing;
using namespace MML::Testing;

using Catch::Matchers::WithinAbs;
// using statement removed - RealWithinRel is a function, not a type

namespace MML::Tests::Base::PolynomTests
{
  TEST_CASE("Polynom::coef.order")
  {
    PolynomRealFunc poly4({REAL(REAL(1.0)), REAL(REAL(0.0)), REAL(REAL(0.0)), REAL(REAL(0.0)), REAL(REAL(5.0))});   // polynom 1 + 5 * x^4
    
    REQUIRE(poly4.GetDegree() == 4);
    REQUIRE(poly4[0] == REAL(REAL(1.0)));
    REQUIRE(poly4[4] == REAL(REAL(5.0)));

    REQUIRE(poly4(REAL(REAL(2.0))) == 5 * 2*2*2*2 + 1);
  }

  TEST_CASE("Polynom::Init, GetDegree, operator[]", "[simple]")
  {
    PolynomRealFunc poly_empty;
    REQUIRE(poly_empty.GetDegree() == -1);

    PolynomRealFunc poly3(3);
    REQUIRE(poly3.GetDegree() == 3);

    PolynomRealFunc poly4({REAL(REAL(1.0)), REAL(REAL(2.0)), REAL(REAL(3.0)), REAL(REAL(4.0)), REAL(REAL(5.0))});
    REQUIRE(poly4.GetDegree() == 4);
    REQUIRE(poly4[0] == REAL(REAL(1.0)));
    REQUIRE(poly4[1] == REAL(REAL(2.0)));
    REQUIRE(poly4[2] == REAL(REAL(3.0)));
    REQUIRE(poly4[3] == REAL(REAL(4.0)));
    REQUIRE(poly4[4] == REAL(REAL(5.0)));

    std::vector<Real> coefs{-REAL(REAL(1.0)), -REAL(REAL(2.0)), -REAL(REAL(3.0)), -REAL(REAL(4.0)), -REAL(REAL(5.0))};
    PolynomRealFunc poly4b(coefs);
    REQUIRE(poly4b.GetDegree() == 4);
    REQUIRE(poly4b[0] == -REAL(REAL(1.0)));
    REQUIRE(poly4b[1] == -REAL(REAL(2.0)));
    REQUIRE(poly4b[2] == -REAL(REAL(3.0)));
    REQUIRE(poly4b[3] == -REAL(REAL(4.0)));
    REQUIRE(poly4b[4] == -REAL(REAL(5.0)));

    // negativan argument
    // prazan vektor
    // TODO - polynom test - implement for Complex and Matrix
  }

  TEST_CASE("Polynom::SetDegree", "[simple]")
  {
    PolynomRealFunc poly;
    REQUIRE(poly.GetDegree() == -1);

    poly.SetDegree(3);
    REQUIRE(poly.GetDegree() == 3);
  }

  TEST_CASE("Polynom::operator()", "[simple]")
  {
    PolynomRealFunc pol_constant({REAL(REAL(3.5))});

    REQUIRE(pol_constant(REAL(REAL(0.0))) == REAL(REAL(3.5)));
    REQUIRE_THAT(pol_constant(REAL(REAL(1.0))) , RealWithinRel(REAL(REAL(3.5)), 1e-5));
    REQUIRE_THAT(pol_constant(-REAL(REAL(200.0))) , RealWithinRel(REAL(REAL(3.5)), 1e-5));
    REQUIRE_THAT(pol_constant(1.0e-10) , RealWithinRel(REAL(REAL(3.5)), 1e-5));
    REQUIRE_THAT(pol_constant(1.0e+10) , RealWithinRel(REAL(REAL(3.5)), 1e-5));

    PolynomRealFunc pol_lin({REAL(REAL(3.5)), REAL(REAL(1.0))});

    REQUIRE_THAT(pol_lin(REAL(REAL(0.0))) , RealWithinRel(REAL(REAL(3.5)), 1e-5));
    REQUIRE(pol_lin(REAL(REAL(1.0))) == REAL(REAL(4.5)));
    REQUIRE_THAT(pol_lin(-REAL(REAL(100.0))) , RealWithinRel(-REAL(REAL(96.5)), REAL(1e-5)));
  }

  TEST_CASE("Polynom::operator==", "[simple]")
  {
    PolynomRealFunc poly3(3);

    PolynomRealFunc poly4({REAL(REAL(1.0)), REAL(REAL(2.0)), REAL(REAL(3.0)), REAL(REAL(4.0)), REAL(REAL(5.0))});
    PolynomRealFunc poly4a({REAL(REAL(1.0)), REAL(REAL(2.0)), REAL(REAL(3.0)), REAL(REAL(4.0)), REAL(REAL(5.0))});

    std::vector<Real> coefs{-REAL(REAL(1.0)), -REAL(REAL(2.0)), -REAL(REAL(3.0)), -REAL(REAL(4.0)), -REAL(REAL(5.0))};
    PolynomRealFunc poly4b(coefs);

    REQUIRE(poly4 == poly4a);
    REQUIRE(poly4 != poly4b);
    REQUIRE(poly4 == poly4b * (-1));
  }

  TEST_CASE("Polynom::operator== with different degrees (trailing zeros)", "[equality]")
  {
    // Test that polynomials with trailing zeros are considered equal
    PolynomRealFunc p1({REAL(REAL(1.0)), REAL(REAL(2.0)), REAL(REAL(3.0))});           // 1 + 2x + 3x²
    PolynomRealFunc p2({REAL(REAL(1.0)), REAL(REAL(2.0)), REAL(REAL(3.0)), REAL(REAL(0.0))});      // 1 + 2x + 3x² + 0x³
    PolynomRealFunc p3({REAL(REAL(1.0)), REAL(REAL(2.0)), REAL(REAL(3.0)), REAL(REAL(0.0)), REAL(REAL(0.0))}); // 1 + 2x + 3x² + 0x³ + 0x⁴

    REQUIRE(p1 == p2);
    REQUIRE(p2 == p3);
    REQUIRE(p1 == p3);
    REQUIRE(p1.IsEqual(p2));
    REQUIRE(p2.IsEqual(p3));

    // Test with single coefficient
    PolynomRealFunc c1({REAL(REAL(5.0))});
    PolynomRealFunc c2({REAL(REAL(5.0)), REAL(REAL(0.0)), REAL(REAL(0.0))});
    REQUIRE(c1 == c2);

    // Test that different polynomials are still not equal
    PolynomRealFunc p4({REAL(REAL(1.0)), REAL(REAL(2.0)), REAL(REAL(3.0)), REAL(REAL(1.0))}); // 1 + 2x + 3x² + x³
    REQUIRE(p1 != p4);
    REQUIRE(p2 != p4);

    // Test with zero polynomial
    PolynomRealFunc zero1({REAL(REAL(0.0))});
    PolynomRealFunc zero2({REAL(REAL(0.0)), REAL(REAL(0.0)), REAL(REAL(0.0))});
    REQUIRE(zero1 == zero2);
  }

  TEST_CASE("Polynom::operator+", "[simple]")
  {
    PolynomRealFunc poly1({REAL(REAL(3.5))});
    PolynomRealFunc poly2({REAL(REAL(1.0)), REAL(REAL(0.5)), REAL(REAL(2.0))});

    PolynomRealFunc poly3 = poly1 + poly2;

    REQUIRE(poly3.GetDegree() == 2);
    REQUIRE(poly3[0] == REAL(REAL(4.5)));
    REQUIRE(poly3[1] == REAL(REAL(0.5)));
    REQUIRE(poly3[2] == REAL(REAL(2.0)));
  }
  TEST_CASE("Polynom::operator-", "[simple]")
  {
    PolynomRealFunc poly1({REAL(REAL(3.5))});
    PolynomRealFunc poly2({REAL(REAL(1.0)), REAL(REAL(0.5)), REAL(REAL(2.0))});

    PolynomRealFunc poly3 = poly1 - poly2;

    REQUIRE(poly3.GetDegree() == 2);
    REQUIRE(poly3[0] == REAL(REAL(2.5)));
    REQUIRE(poly3[1] == -REAL(REAL(0.5)));
    REQUIRE(poly3[2] == -REAL(REAL(2.0)));
  }
  TEST_CASE("Polynom::operator*", "[simple]")
  {
    PolynomRealFunc poly1({2, 1, 3});
    PolynomRealFunc poly2({1, 2});

    PolynomRealFunc poly3 = poly1 * poly2;
    REQUIRE(poly3.GetDegree() == 3);
    REQUIRE(poly3[0] == 2);
    REQUIRE(poly3[1] == 5);
    REQUIRE(poly3[2] == 5);
    REQUIRE(poly3[3] == 6);
  }
  TEST_CASE("Polynom::operator*(Poly, CoefType)", "[simple]")
  {
    PolynomRealFunc poly({1, 2, 3});

    PolynomRealFunc poly2 = poly * REAL(REAL(3.0));
    REQUIRE(poly2.GetDegree() == 2);
    REQUIRE(poly2[0] == REAL(REAL(3.0)));
    REQUIRE(poly2[1] == REAL(REAL(6.0)));
    REQUIRE(poly2[2] == REAL(REAL(9.0)));
  }
  TEST_CASE("Polynom::operator*(CoefType, Poly)", "[simple]")
  {
    PolynomRealFunc poly({1, 2, 3});

    PolynomRealFunc poly2 = REAL(REAL(3.0)) * poly;
    REQUIRE(poly2.GetDegree() == 2);
    REQUIRE(poly2[0] == REAL(REAL(3.0)));
    REQUIRE(poly2[1] == REAL(REAL(6.0)));
    REQUIRE(poly2[2] == REAL(REAL(9.0)));
  }
  TEST_CASE("Polynom::operator/(Poly, CoefType)", "[simple]")
  {
    PolynomRealFunc poly({2, 4, 6});

    PolynomRealFunc poly2 = poly / REAL(REAL(2.0));
    REQUIRE(poly2.GetDegree() == 2);
    REQUIRE(poly2[0] == REAL(REAL(1.0)));
    REQUIRE(poly2[1] == REAL(REAL(2.0)));
    REQUIRE(poly2[2] == REAL(REAL(3.0)));
  }
  
  TEST_CASE("Polynom::to_string", "[simple]")
  {
    PolynomRealFunc poly({1, 2, 3});

    REQUIRE(poly.to_string(10, 5) ==   "   3.00000x^2 + 2.00000x + 1.00000");    
  }

  TEST_CASE("Polynom::poldiv_1", "[simple]")
  {
    PolynomRealFunc p1({-REAL(REAL(1.0)), REAL(REAL(1.0))}); // x - 1
    PolynomRealFunc p2({-REAL(REAL(2.0)), REAL(REAL(1.0))}); // x - 2
    PolynomRealFunc p3({-REAL(REAL(3.0)), REAL(REAL(1.0))}); // x - 3

    PolynomRealFunc p = p1 * p2 * p3;

    REQUIRE(p.GetDegree() == 3);
    REQUIRE(p[0] == -REAL(REAL(6.0)));
    REQUIRE(p[1] == REAL(REAL(11.0)));
    REQUIRE(p[2] == -REAL(REAL(6.0)));
    REQUIRE(p[3] == REAL(REAL(1.0)));

    PolynomRealFunc res, rem;

    PolynomRealFunc::poldiv(p, p1, res, rem);
    REQUIRE(p.IsEqual(res * p1 + rem));
    REQUIRE(res == (p2 * p3));

    PolynomRealFunc::poldiv(p, p2, res, rem);
    REQUIRE(p.IsEqual(res * p2 + rem));
    REQUIRE(res == (p1 * p3));

    PolynomRealFunc::poldiv(p, p3, res, rem);
    REQUIRE(p.IsEqual(res * p3 + rem));
    REQUIRE(res == (p1 * p2));
  }
  
  TEST_CASE("Polynom::poldiv_2", "[simple]")
  {
    PolynomRealFunc p_1({REAL(REAL(2.0)), REAL(REAL(0.0)), REAL(REAL(1.0))});
    PolynomRealFunc p_2({-REAL(REAL(4.0)), REAL(REAL(1.0))});
    PolynomRealFunc p_quot({REAL(REAL(4.0)), REAL(REAL(1.0))});
    PolynomRealFunc p_rem({REAL(REAL(18.0))});

    PolynomRealFunc res, rem;

    PolynomRealFunc::poldiv(p_1, p_2, res, rem);

    REQUIRE(p_1 == res * p_2 + rem);
    REQUIRE(res == p_quot);
    REQUIRE(rem == p_rem);

    PolynomRealFunc r_1({-REAL(REAL(1.0)), -REAL(REAL(5.0)), REAL(REAL(2.0))});
    PolynomRealFunc r_2({-REAL(REAL(3.0)), REAL(REAL(1.0))});
    PolynomRealFunc r_quot({REAL(REAL(1.0)), REAL(REAL(2.0))});
    PolynomRealFunc r_rem({REAL(REAL(2.0))});

    PolynomRealFunc::poldiv(r_1, r_2, res, rem);

    REQUIRE(r_1 == res * r_2 + rem);
    REQUIRE(res == r_quot);
    REQUIRE(rem == r_rem);
  }
  
  TEST_CASE("Polynom::poldiv_3", "[simple]")
  {
    PolynomRealFunc p3_u({-REAL(REAL(1.0)), -REAL(REAL(5.0)) });
    PolynomRealFunc p3_v({-REAL(REAL(3.0)), REAL(REAL(1.0)), REAL(REAL(2.0))});

    PolynomRealFunc res, rem;

    PolynomRealFunc::poldiv(p3_u, p3_v, res, rem);

    REQUIRE(res.GetDegree() == -1);
    REQUIRE(res.IsNullPolynom() == true);
    REQUIRE(rem == p3_u);
  }

  /*********************************************************************/
  /*****                 Polynomial Derivatives                    *****/
  /*********************************************************************/
  TEST_CASE("Polynom::Derivatives", "[derivatives]")
  {
    // Test simple polynomial: 3x^2 + 2x + 1
    PolynomRealFunc poly({REAL(REAL(1.0)), REAL(REAL(2.0)), REAL(REAL(3.0))});
    Vector<Real> pd(4);  // Value + 3 derivatives
    
    // Evaluate at x = 2
    poly.Derive(REAL(REAL(2.0)), pd);
    
    // f(2) = 3*4 + 2*2 + 1 = 17
    REQUIRE_THAT(pd[0] , RealWithinRel(REAL(REAL(17.0)), 1e-5));
    
    // f'(2) = 6*2 + 2 = 14
    REQUIRE_THAT(pd[1] , RealWithinRel(REAL(REAL(14.0)), 1e-5));
    
    // f''(2) = 6
    REQUIRE_THAT(pd[2] , RealWithinRel(REAL(REAL(6.0)), 1e-5));
    
    // f'''(2) = 0
    REQUIRE_THAT(pd[3] , RealWithinRel(REAL(REAL(0.0)), 1e-5));
  }

  TEST_CASE("Polynom::Derivatives_quartic", "[derivatives]")
  {
    // Test quartic: x^4 + 2x^3 - 3x^2 + 4x - 5
    PolynomRealFunc poly({-REAL(REAL(5.0)), REAL(REAL(4.0)), -REAL(REAL(3.0)), REAL(REAL(2.0)), REAL(REAL(1.0))});
    Vector<Real> pd(5);  // Value + 4 derivatives
    
    // Evaluate at x = 1
    poly.Derive(REAL(REAL(1.0)), pd);
    
    // f(1) = 1 + 2 - 3 + 4 - 5 = -1
    REQUIRE_THAT(pd[0] , RealWithinRel(-REAL(REAL(1.0)), REAL(1e-5)));
    
    // f'(1) = 4 + 6 - 6 + 4 = 8
    REQUIRE_THAT(pd[1] , RealWithinRel(REAL(REAL(8.0)), 1e-5));
    
    // f''(1) = 12 + 12 - 6 = 18
    REQUIRE_THAT(pd[2] , RealWithinRel(REAL(REAL(18.0)), 1e-5));
    
    // f'''(1) = 24 + 12 = 36
    REQUIRE_THAT(pd[3] , RealWithinRel(REAL(REAL(36.0)), 1e-5));
    
    // f''''(1) = 24
    REQUIRE_THAT(pd[4] , RealWithinRel(REAL(REAL(24.0)), 1e-5));
  }

  TEST_CASE("Polynom::Derivatives_constant", "[derivatives]")
  {
    // Constant polynomial
    PolynomRealFunc poly({REAL(REAL(5.0))});
    Vector<Real> pd(3);
    
    poly.Derive(REAL(REAL(10.0)), pd);
    
    REQUIRE_THAT(pd[0] , RealWithinRel(REAL(REAL(5.0)), 1e-5));  // f(x) = 5
    REQUIRE_THAT(pd[1] , RealWithinRel(REAL(REAL(0.0)), 1e-5));  // f'(x) = 0
    REQUIRE_THAT(pd[2] , RealWithinRel(REAL(REAL(0.0)), 1e-5));  // f''(x) = 0
  }

  /*********************************************************************/
  /*****           Polynomial Special Constructors                 *****/
  /*********************************************************************/
  TEST_CASE("Polynom::Static_constructors", "[constructors]")
  {
    // Test Zero polynomial
    auto zero = PolynomRealFunc::Zero();
    REQUIRE(zero.GetDegree() == 0);
    REQUIRE_THAT(zero(REAL(REAL(5.0))) , RealWithinRel(REAL(REAL(0.0)), 1e-5));
    
    // Test Monomial x^3
    auto monomial = PolynomRealFunc::Monomial(3);
    REQUIRE(monomial.GetDegree() == 3);
    REQUIRE_THAT(monomial(REAL(REAL(2.0))) , RealWithinRel(REAL(REAL(8.0)), 1e-5));   // 2^3 = 8
    REQUIRE_THAT(monomial(REAL(REAL(3.0))) , RealWithinRel(REAL(REAL(27.0)), 1e-5));  // 3^3 = 27
    
    // Test Constant polynomial
    auto constant = PolynomRealFunc::Constant(REAL(REAL(7.5)));
    REQUIRE(constant.GetDegree() == 0);
    REQUIRE_THAT(constant(REAL(REAL(100.0))) , RealWithinRel(REAL(REAL(7.5)), 1e-5));
    REQUIRE_THAT(constant(-REAL(REAL(50.0))) , RealWithinRel(REAL(REAL(7.5)), 1e-5));
    
    // Test Linear polynomial: 3x + 5
    auto linear = PolynomRealFunc::Linear(REAL(REAL(3.0)), REAL(REAL(5.0)));
    REQUIRE(linear.GetDegree() == 1);
    REQUIRE_THAT(linear(REAL(REAL(0.0))) , RealWithinRel(REAL(REAL(5.0)), 1e-5));
    REQUIRE_THAT(linear(REAL(REAL(1.0))) , RealWithinRel(REAL(REAL(8.0)), 1e-5));
    REQUIRE_THAT(linear(REAL(REAL(2.0))) , RealWithinRel(REAL(REAL(11.0)), 1e-5));
  }

  /*********************************************************************/
  /*****              Polynomial Root Finding                      *****/
  /*********************************************************************/
  TEST_CASE("Polynom::Roots_linear", "[roots]")
  {
    // Linear: 2x + 4 = 0, root at x = -2
    PolynomRealFunc linear({REAL(REAL(4.0)), REAL(REAL(2.0))});
    
    REQUIRE_THAT(linear(-REAL(REAL(2.0))) , RealWithinRel(REAL(REAL(0.0)), 1e-5));
  }

  TEST_CASE("Polynom::Roots_quadratic", "[roots]")
  {
    // Quadratic: x^2 - 5x + 6 = (x-2)(x-3), roots at x=2 and x=3
    PolynomRealFunc quad({REAL(REAL(6.0)), -REAL(REAL(5.0)), REAL(REAL(1.0))});
    
    REQUIRE_THAT(quad(REAL(REAL(2.0))) , RealWithinRel(REAL(REAL(0.0)), 1e-5));
    REQUIRE_THAT(quad(REAL(REAL(3.0))) , RealWithinRel(REAL(REAL(0.0)), 1e-5));
  }

  TEST_CASE("Polynom::Roots_cubic", "[roots]")
  {
    // Cubic: (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6
    PolynomRealFunc cubic({-REAL(REAL(6.0)), REAL(REAL(11.0)), -REAL(REAL(6.0)), REAL(REAL(1.0))});
    
    REQUIRE_THAT(cubic(REAL(REAL(1.0))) , RealWithinRel(REAL(REAL(0.0)), 1e-5));
    REQUIRE_THAT(cubic(REAL(REAL(2.0))) , RealWithinRel(REAL(REAL(0.0)), 1e-5));
    REQUIRE_THAT(cubic(REAL(REAL(3.0))) , RealWithinRel(REAL(REAL(0.0)), 1e-5));
  }

  /*********************************************************************/
  /*****           Polynomial Arithmetic Properties                *****/
  /*********************************************************************/
  TEST_CASE("Polynom::Arithmetic_properties", "[properties]")
  {
    PolynomRealFunc p1({REAL(REAL(1.0)), REAL(REAL(2.0)), REAL(REAL(3.0))});  // 3x^2 + 2x + 1
    PolynomRealFunc p2({REAL(REAL(4.0)), REAL(REAL(5.0))});       // 5x + 4
    
    // Test commutativity of addition
    auto sum1 = p1 + p2;
    auto sum2 = p2 + p1;
    REQUIRE(sum1 == sum2);
    
    // Test commutativity of multiplication
    auto prod1 = p1 * p2;
    auto prod2 = p2 * p1;
    REQUIRE(prod1 == prod2);
    
    // Test associativity of addition
    PolynomRealFunc p3({REAL(REAL(1.0)), REAL(REAL(1.0))});  // x + 1
    auto sum_assoc1 = (p1 + p2) + p3;
    auto sum_assoc2 = p1 + (p2 + p3);
    REQUIRE(sum_assoc1 == sum_assoc2);
  }

  TEST_CASE("Polynom::Scalar_operations", "[scalar]")
  {
    PolynomRealFunc poly({REAL(REAL(1.0)), REAL(REAL(2.0)), REAL(REAL(3.0))});  // 3x^2 + 2x + 1
    
    // Test addition with scalar
    auto p_plus_5 = poly + REAL(REAL(5.0));
    REQUIRE_THAT(p_plus_5[0] , RealWithinRel(REAL(REAL(6.0)), 1e-5));  // constant term becomes 6
    REQUIRE_THAT(p_plus_5[1] , RealWithinRel(REAL(REAL(2.0)), 1e-5));  // other terms unchanged
    REQUIRE_THAT(p_plus_5[2] , RealWithinRel(REAL(REAL(3.0)), 1e-5));
    
    // Test scalar + polynomial
    auto five_plus_p = REAL(REAL(5.0)) + poly;
    REQUIRE(five_plus_p == p_plus_5);
    
    // Test subtraction with scalar
    auto p_minus_3 = poly - REAL(REAL(3.0));
    REQUIRE_THAT(p_minus_3[0] , RealWithinRel(-REAL(REAL(2.0)), REAL(1e-5)));  // 1 - 3 = -2
    REQUIRE_THAT(p_minus_3[1] , RealWithinRel(REAL(REAL(2.0)), 1e-5));
    
    // Test scalar - polynomial (using multiplication by -1 and addition)
    auto neg_p = poly * (-REAL(REAL(1.0)));
    auto three_minus_p = REAL(REAL(3.0)) + neg_p;
    REQUIRE_THAT(three_minus_p[0] , RealWithinRel(REAL(REAL(2.0)), 1e-5));  // 3 + (-1) = 2
    REQUIRE_THAT(three_minus_p[1] , RealWithinRel(-REAL(REAL(2.0)), REAL(1e-5))); // 0 + (-2) = -2
    REQUIRE_THAT(three_minus_p[2] , RealWithinRel(-REAL(REAL(3.0)), REAL(1e-5))); // 0 + (-3) = -3
  }

  TEST_CASE("Polynom::Leading_and_constant_terms", "[properties]")
  {
    PolynomRealFunc poly({REAL(REAL(5.0)), REAL(REAL(3.0)), REAL(REAL(2.0)), REAL(REAL(7.0))});  // 7x^3 + 2x^2 + 3x + 5
    
    REQUIRE_THAT(poly.leadingTerm() , RealWithinRel(REAL(REAL(7.0)), 1e-5));
    REQUIRE_THAT(poly.constantTerm() , RealWithinRel(REAL(REAL(5.0)), 1e-5));
    
    // Test empty polynomial
    PolynomRealFunc empty;
    REQUIRE_THAT(empty.leadingTerm() , RealWithinRel(REAL(REAL(0.0)), 1e-5));
    REQUIRE_THAT(empty.constantTerm() , RealWithinRel(REAL(REAL(0.0)), 1e-5));
  }

  TEST_CASE("Polynom::Reduce_degree", "[properties]")
  {
    // Polynomial with trailing zeros: 2x^2 + 3x + 1 + 0x^3 + 0x^4
    PolynomRealFunc poly({REAL(REAL(1.0)), REAL(REAL(3.0)), REAL(REAL(2.0)), REAL(REAL(0.0)), REAL(REAL(0.0))});
    
    REQUIRE(poly.GetDegree() == 4);  // Before reduction
    
    poly.Reduce();
    
    REQUIRE(poly.GetDegree() == 2);  // After reduction
    REQUIRE_THAT(poly[0] , RealWithinRel(REAL(REAL(1.0)), 1e-5));
    REQUIRE_THAT(poly[1] , RealWithinRel(REAL(REAL(3.0)), 1e-5));
    REQUIRE_THAT(poly[2] , RealWithinRel(REAL(REAL(2.0)), 1e-5));
  }

  TEST_CASE("Polynom::FromValues - Quadratic (degree 2)", "[interpolation]")
  {
    // Test polynomial: P(x) = 2 + 3x - x^2
    // P(0) = 2, P(1) = 4, P(2) = 4
    std::vector<Real> x_vals = {REAL(REAL(0.0)), REAL(REAL(1.0)), REAL(REAL(2.0))};
    std::vector<Real> y_vals = {REAL(REAL(2.0)), REAL(REAL(4.0)), REAL(REAL(4.0))};
    
    PolynomReal poly = PolynomReal::FromValues(x_vals, y_vals);
    
    REQUIRE(poly.GetDegree() == 2);
    
    // Verify coefficients: 2, 3, -1
    REQUIRE_THAT(poly[0], WithinAbs(REAL(REAL(2.0)), 1e-10));
    REQUIRE_THAT(poly[1], WithinAbs(REAL(REAL(3.0)), 1e-10));
    REQUIRE_THAT(poly[2], WithinAbs(-REAL(REAL(1.0)), REAL(1e-10)));
    
    // Verify polynomial passes through all points
    REQUIRE_THAT(poly(REAL(REAL(0.0))), WithinAbs(REAL(REAL(2.0)), 1e-10));
    REQUIRE_THAT(poly(REAL(REAL(1.0))), WithinAbs(REAL(REAL(4.0)), 1e-10));
    REQUIRE_THAT(poly(REAL(REAL(2.0))), WithinAbs(REAL(REAL(4.0)), 1e-10));
    
    // Test intermediate point
    REQUIRE_THAT(poly(REAL(REAL(0.5))), WithinAbs(REAL(REAL(3.25)), 1e-10));  // 2 + REAL(REAL(1.5)) - REAL(REAL(0.25)) = REAL(REAL(3.25))
  }

  TEST_CASE("Polynom::FromValues - Cubic (degree 3)", "[interpolation]")
  {
    // Test: Interpolate through 4 points to get a cubic polynomial
    // P(0) = 1, P(1) = 1, P(-1) = 1, P(2) = 9
    std::vector<Real> x_vals = {REAL(REAL(0.0)), REAL(REAL(1.0)), -REAL(REAL(1.0)), REAL(REAL(2.0))};
    std::vector<Real> y_vals = {REAL(REAL(1.0)), REAL(REAL(1.0)), REAL(REAL(1.0)), REAL(REAL(9.0))};
    
    PolynomReal poly = PolynomReal::FromValues(x_vals, y_vals);
    
    REQUIRE(poly.GetDegree() == 3);
    
    // The key requirement is that it passes through all points
    // (The specific coefficients depend on the interpolation algorithm)
    REQUIRE_THAT(poly(REAL(REAL(0.0))), WithinAbs(REAL(REAL(1.0)), 1e-10));
    REQUIRE_THAT(poly(REAL(REAL(1.0))), WithinAbs(REAL(REAL(1.0)), 1e-10));
    REQUIRE_THAT(poly(-REAL(REAL(1.0))), WithinAbs(REAL(REAL(1.0)), 1e-10));
    REQUIRE_THAT(poly(REAL(REAL(2.0))), WithinAbs(REAL(REAL(9.0)), 1e-10));
    
    // Test another point to verify interpolation
    // The polynomial passing through these points should give consistent intermediate values
    Real y_half = poly(REAL(REAL(0.5)));
    REQUIRE(y_half > REAL(REAL(0.0)));  // Should be positive based on the data
    REQUIRE(y_half < REAL(REAL(9.0)));  // Should be less than maximum value
  }

  TEST_CASE("Polynom::FromValues - Quintic (degree 5)", "[interpolation]")
  {
    // Test polynomial: P(x) = 1 + x + x^2 + x^3 + x^4 + x^5
    // Sample at 6 equally spaced points
    std::vector<Real> x_vals = {REAL(REAL(0.0)), REAL(REAL(0.5)), REAL(REAL(1.0)), REAL(REAL(1.5)), REAL(REAL(2.0)), REAL(REAL(2.5))};
    std::vector<Real> y_vals;
    
    // Calculate y values for the polynomial
    for (Real x : x_vals) {
      Real y = REAL(REAL(1.0)) + x + x*x + x*x*x + x*x*x*x + x*x*x*x*x;
      y_vals.push_back(y);
    }
    
    PolynomReal poly = PolynomReal::FromValues(x_vals, y_vals);
    
    REQUIRE(poly.GetDegree() == 5);
    
    // Verify coefficients (all should be REAL(REAL(1.0)))
    REQUIRE_THAT(poly[0], WithinAbs(REAL(REAL(1.0)), 1e-8));
    REQUIRE_THAT(poly[1], WithinAbs(REAL(REAL(1.0)), 1e-8));
    REQUIRE_THAT(poly[2], WithinAbs(REAL(REAL(1.0)), 1e-8));
    REQUIRE_THAT(poly[3], WithinAbs(REAL(REAL(1.0)), 1e-8));
    REQUIRE_THAT(poly[4], WithinAbs(REAL(REAL(1.0)), 1e-8));
    REQUIRE_THAT(poly[5], WithinAbs(REAL(REAL(1.0)), 1e-8));
    
    // Verify polynomial passes through all sample points
    for (size_t i = 0; i < x_vals.size(); i++) {
      REQUIRE_THAT(poly(x_vals[i]), WithinAbs(y_vals[i], REAL(1e-8)));
    }
    
    // Test intermediate point
    Real x_test = REAL(REAL(1.25));
    Real y_expected = REAL(REAL(1.0)) + x_test + x_test*x_test + x_test*x_test*x_test 
                      + x_test*x_test*x_test*x_test + x_test*x_test*x_test*x_test*x_test;
    REQUIRE_THAT(poly(x_test), WithinAbs(y_expected, REAL(1e-7)));
  }

  TEST_CASE("Polynom::FromValues - Degree 8", "[interpolation]")
  {
    // Test polynomial: P(x) = 3 - 2x + 4x^2 - x^3 + 2x^4 + x^5 - 3x^6 + x^7 + x^8
    std::vector<Real> coeffs = {REAL(REAL(3.0)), -REAL(REAL(2.0)), REAL(REAL(4.0)), -REAL(REAL(1.0)), REAL(REAL(2.0)), REAL(REAL(1.0)), -REAL(REAL(3.0)), REAL(REAL(1.0)), REAL(REAL(1.0))};
    
    // Helper to evaluate the original polynomial
    auto eval_poly = [&coeffs](Real x) -> Real {
      Real y = coeffs[0];
      Real x_pow = x;
      for (size_t i = 1; i < coeffs.size(); i++) {
        y += coeffs[i] * x_pow;
        x_pow *= x;
      }
      return y;
    };
    
    // Sample at 9 points (need n+1 points for degree n)
    std::vector<Real> x_vals = {-REAL(REAL(2.0)), -REAL(REAL(1.5)), -REAL(REAL(1.0)), -REAL(REAL(0.5)), REAL(REAL(0.0)), REAL(REAL(0.5)), REAL(REAL(1.0)), REAL(REAL(1.5)), REAL(REAL(2.0))};
    std::vector<Real> y_vals;
    
    // Calculate y values using the polynomial
    for (Real x : x_vals) {
      y_vals.push_back(eval_poly(x));
    }
    
    PolynomReal poly = PolynomReal::FromValues(x_vals, y_vals);
    
    REQUIRE(poly.GetDegree() == 8);
    
    // Verify coefficients (with reasonable tolerance for high-degree polynomials)
    for (size_t i = 0; i < coeffs.size(); i++) {
      REQUIRE_THAT(poly[i], WithinAbs(coeffs[i], REAL(1e-6)));
    }
    
    // Verify polynomial passes through all sample points
    for (size_t i = 0; i < x_vals.size(); i++) {
      REQUIRE_THAT(poly(x_vals[i]), WithinAbs(y_vals[i], REAL(1e-6)));
    }
    
    // Test intermediate points using the original polynomial for comparison
    Real x_test1 = -REAL(REAL(0.75));
    Real y_expected1 = eval_poly(x_test1);
    REQUIRE_THAT(poly(x_test1), WithinAbs(y_expected1, REAL(1e-5)));
    
    Real x_test2 = REAL(REAL(1.25));
    Real y_expected2 = eval_poly(x_test2);
    REQUIRE_THAT(poly(x_test2), WithinAbs(y_expected2, REAL(1e-5)));
  }

  TEST_CASE("Polynom::FromValues - Edge cases", "[interpolation]")
  {
    SECTION("Single point (constant polynomial)")
    {
      std::vector<Real> x_vals = {REAL(REAL(1.0))};
      std::vector<Real> y_vals = {REAL(REAL(5.0))};
      
      PolynomReal poly = PolynomReal::FromValues(x_vals, y_vals);
      
      REQUIRE(poly.GetDegree() == 0);
      REQUIRE_THAT(poly[0], WithinAbs(REAL(REAL(5.0)), 1e-10));
      REQUIRE_THAT(poly(REAL(REAL(1.0))), WithinAbs(REAL(REAL(5.0)), 1e-10));
      REQUIRE_THAT(poly(REAL(REAL(10.0))), WithinAbs(REAL(REAL(5.0)), 1e-10));  // Constant everywhere
    }
    
    SECTION("Linear polynomial")
    {
      std::vector<Real> x_vals = {REAL(REAL(0.0)), REAL(REAL(1.0))};
      std::vector<Real> y_vals = {REAL(REAL(2.0)), REAL(REAL(5.0))};  // y = 2 + 3x
      
      PolynomReal poly = PolynomReal::FromValues(x_vals, y_vals);
      
      REQUIRE(poly.GetDegree() == 1);
      REQUIRE_THAT(poly[0], WithinAbs(REAL(REAL(2.0)), 1e-10));
      REQUIRE_THAT(poly[1], WithinAbs(REAL(REAL(3.0)), 1e-10));
      REQUIRE_THAT(poly(REAL(REAL(0.5))), WithinAbs(REAL(REAL(3.5)), 1e-10));
    }
    
    SECTION("Invalid inputs - mismatched sizes")
    {
      std::vector<Real> x_vals = {REAL(REAL(0.0)), REAL(REAL(1.0)), REAL(REAL(2.0))};
      std::vector<Real> y_vals = {REAL(REAL(1.0)), REAL(REAL(2.0))};
      
      REQUIRE_THROWS_AS(PolynomReal::FromValues(x_vals, y_vals), std::invalid_argument);
    }
    
    SECTION("Invalid inputs - empty arrays")
    {
      std::vector<Real> x_vals;
      std::vector<Real> y_vals;
      
      REQUIRE_THROWS_AS(PolynomReal::FromValues(x_vals, y_vals), std::invalid_argument);
    }
  }

  // ========== EDGE CASE TESTS ==========

  TEST_CASE("Polynom::EdgeCase_zero_polynomial", "[Polynom][edge_cases]")
  {
    // Zero polynomial: p(x) = 0 for all x
    PolynomReal zero({REAL(0.0)});
    
    REQUIRE(zero.GetDegree() == 0);
    REQUIRE(zero(REAL(0.0)) == REAL(0.0));
    REQUIRE(zero(REAL(100.0)) == REAL(0.0));
    REQUIRE(zero(REAL(-50.0)) == REAL(0.0));
    
    // Zero + anything = anything
    PolynomReal other({REAL(1.0), REAL(2.0), REAL(3.0)});
    PolynomReal sum = zero + other;
    REQUIRE(sum(REAL(1.0)) == other(REAL(1.0)));
    
    // Zero * anything = zero
    PolynomReal prod = zero * other;
    REQUIRE(prod(REAL(5.0)) == REAL(0.0));
  }

  TEST_CASE("Polynom::EdgeCase_constant_polynomial", "[Polynom][edge_cases]")
  {
    // Constant polynomial: p(x) = 7 for all x
    PolynomReal constant({REAL(7.0)});
    
    REQUIRE(constant.GetDegree() == 0);
    REQUIRE(constant(REAL(0.0)) == REAL(7.0));
    REQUIRE(constant(REAL(1000.0)) == REAL(7.0));
    REQUIRE(constant(REAL(-999.0)) == REAL(7.0));
    
    // Derivative of constant is zero (using Derive method)
    Vector<Real> derivs(2);  // value and 1st derivative
    constant.Derive(REAL(5.0), derivs);
    REQUIRE(derivs[1] == REAL(0.0));  // f'(5) = 0
  }

  TEST_CASE("Polynom::EdgeCase_high_degree", "[Polynom][edge_cases]")
  {
    // High degree polynomial: x^10
    std::vector<Real> coeffs(11, REAL(0.0));
    coeffs[10] = REAL(1.0);  // x^10
    PolynomReal highDeg(coeffs);
    
    REQUIRE(highDeg.GetDegree() == 10);
    REQUIRE(highDeg(REAL(1.0)) == REAL(1.0));     // 1^10 = 1
    REQUIRE(highDeg(REAL(2.0)) == REAL(1024.0));  // 2^10 = 1024
    REQUIRE(highDeg(REAL(0.0)) == REAL(0.0));     // 0^10 = 0
  }

  TEST_CASE("Polynom::EdgeCase_negative_coefficients", "[Polynom][edge_cases]")
  {
    // p(x) = -x^2 + 2x - 1 = -(x-1)^2, always <= 0
    PolynomReal negQuad({REAL(-1.0), REAL(2.0), REAL(-1.0)});
    
    REQUIRE(negQuad(REAL(1.0)) == REAL(0.0));   // Root at x=1
    REQUIRE(negQuad(REAL(0.0)) == REAL(-1.0));
    REQUIRE(negQuad(REAL(2.0)) == REAL(-1.0));
    REQUIRE(negQuad(REAL(3.0)) < REAL(0.0));    // Always negative except at x=1
  }
}