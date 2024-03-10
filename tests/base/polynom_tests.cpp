#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Polynom.h"
#endif

using namespace MML;

TEST_CASE("Test_Polynom_init_GetDegree_operator[]", "[simple]") {
	RealPolynom poly_empty;
	REQUIRE(poly_empty.GetDegree() == -1);

	RealPolynom poly3(3);
	REQUIRE(poly3.GetDegree() == 3);

	RealPolynom poly4({ 1.0, 2.0, 3.0, 4.0, 5.0 });
	REQUIRE(poly4.GetDegree() == 4);
	REQUIRE(poly4[0] == 1.0);
	REQUIRE(poly4[1] == 2.0);
	REQUIRE(poly4[2] == 3.0);
	REQUIRE(poly4[3] == 4.0);
	REQUIRE(poly4[4] == 5.0);

	std::vector<Real> coefs{ -1.0, -2.0, -3.0, -4.0, -5.0 };
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

TEST_CASE("Test_Polynom_SetDegree", "[simple]") {
	RealPolynom poly;
	REQUIRE(poly.GetDegree() == -1);

	poly.SetDegree(3);
	REQUIRE(poly.GetDegree() == 3);
}

TEST_CASE("Test_Polynom_operator()", "[simple]") {
	RealPolynom pol_constant({ 3.5 });

	REQUIRE(pol_constant(0.0) == Approx(3.5));
	REQUIRE(pol_constant(1.0) == Approx(3.5));
	REQUIRE(pol_constant(-200.0) == Approx(3.5));
	REQUIRE(pol_constant(1.0e-10) == Approx(3.5));
	REQUIRE(pol_constant(1.0e+10) == Approx(3.5));

	RealPolynom pol_lin({ 3.5, 1.0 });

	REQUIRE(pol_lin(0.0) == Approx(3.5));
	REQUIRE(pol_lin(1.0) == Approx(4.5));
	REQUIRE(pol_lin(-100.0) == Approx(-96.5));
}

TEST_CASE("Test_Polynom_operator==", "[simple]") {
	RealPolynom poly3(3);

	RealPolynom poly4({ 1.0, 2.0, 3.0, 4.0, 5.0 });
	RealPolynom poly4a({ 1.0, 2.0, 3.0, 4.0, 5.0 });

	std::vector<Real> coefs{ -1.0, -2.0, -3.0, -4.0, -5.0 };
	RealPolynom poly4b(coefs);

	REQUIRE(poly4 == poly4a);
	REQUIRE(poly4 != poly4b);
	REQUIRE(poly4 == poly4b * (-1));
}

// TODO 0.9  - finish these
TEST_CASE("Test_Polynom_operator+", "[simple]") {
	RealPolynom poly({ 3.5 });
}
TEST_CASE("Test_Polynom_operator-", "[simple]") {
	RealPolynom poly({ 3.5 });
}
TEST_CASE("Test_Polynom_operator*", "[simple]") {
	RealPolynom poly({ 3.5 });
}
TEST_CASE("Test_Polynom_operator*(Poly, CoefType)", "[simple]") {
	RealPolynom poly({ 3.5 });
}
TEST_CASE("Test_Polynom_operator*(CoefType, Poly)", "[simple]") {
	RealPolynom poly({ 3.5 });
}
TEST_CASE("Test_Polynom_operator/(Poly, CoefType)", "[simple]") {
	RealPolynom poly({ 3.5 });
}
TEST_CASE("Test_Polynom_to_string", "[simple]") {
	RealPolynom poly({ 3.5 });
}

TEST_CASE("Test_Polynom_poldiv", "[simple]") {
	RealPolynom p1({ -1.0, 1.0 });  // x - 1
	RealPolynom p2({ -2.0, 1.0 });  // x - 2
	RealPolynom p3({ -3.0, 1.0 });  // x - 3

	RealPolynom p = p1 * p2 * p3;

	REQUIRE(p.GetDegree() == 3);
	REQUIRE(p[0] == -6.0);
	REQUIRE(p[1] == 11.0);
	REQUIRE(p[2] == -6.0);
	REQUIRE(p[3] == 1.0);

	RealPolynom res, rem;

	RealPolynom::poldiv(p, p1, res, rem);

	REQUIRE(p == res * p1 + rem);
	REQUIRE(res == p2 * p3);

	RealPolynom p_1({ 2.0, 0.0, 1.0 });
	RealPolynom p_2({ -4.0, 1.0 });
	RealPolynom p_quot({ 4.0, 1.0 });
	RealPolynom p_rem({ 18.0 });

	RealPolynom::poldiv(p_1, p_2, res, rem);

	REQUIRE(p_1 == res * p_2 + rem);
	REQUIRE(res == p_quot);
	REQUIRE(rem == p_rem);

	RealPolynom r_1({ -1.0, -5.0, 2.0 });
	RealPolynom r_2({ -3.0, 1.0 });
	RealPolynom r_quot({ 1.0, 2.0 });
	RealPolynom r_rem({ 2.0 });

	RealPolynom::poldiv(r_1, r_2, res, rem);

	REQUIRE(r_1 == res * r_2 + rem);
	REQUIRE(res == r_quot);
	REQUIRE(rem == r_rem);
}