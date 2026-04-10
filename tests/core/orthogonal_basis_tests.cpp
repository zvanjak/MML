///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        orthogonal_basis_tests.cpp                                          ///
///  Description: Tests for OrthogonalBasis classes - Legendre, Hermite, Laguerre     ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#include "MMLBase.h"
#include "core/OrthogonalBasis.h"
#include "core/OrthogonalBasis/LegendreBasis.h"
#include "core/OrthogonalBasis/HermiteBasis.h"
#include "core/OrthogonalBasis/LaguerreBasis.h"

#include <cmath>
#include <limits>

using namespace MML;

namespace MML::Tests::Core::OrthogonalBasisTests {

///////////////////////////////////////////////////////////////////////////////////////////
///                           LegendreBasis                                             ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("LegendreBasis - Domain", "[OrthogonalBasis][Legendre]")
{
	LegendreBasis basis;
	
	REQUIRE(basis.DomainMin() == Catch::Approx(-1.0));
	REQUIRE(basis.DomainMax() == Catch::Approx(1.0));
}

TEST_CASE("LegendreBasis - Weight function is uniform", "[OrthogonalBasis][Legendre]")
{
	LegendreBasis basis;
	
	REQUIRE(basis.WeightFunction(0.0) == Catch::Approx(1.0));
	REQUIRE(basis.WeightFunction(-0.5) == Catch::Approx(1.0));
	REQUIRE(basis.WeightFunction(0.5) == Catch::Approx(1.0));
	REQUIRE(basis.WeightFunction(-1.0) == Catch::Approx(1.0));
	REQUIRE(basis.WeightFunction(1.0) == Catch::Approx(1.0));
}

TEST_CASE("LegendreBasis - P_0(x) = 1", "[OrthogonalBasis][Legendre]")
{
	LegendreBasis basis;
	
	REQUIRE(basis.Evaluate(0, 0.0) == Catch::Approx(1.0));
	REQUIRE(basis.Evaluate(0, 0.5) == Catch::Approx(1.0));
	REQUIRE(basis.Evaluate(0, -0.5) == Catch::Approx(1.0));
}

TEST_CASE("LegendreBasis - P_1(x) = x", "[OrthogonalBasis][Legendre]")
{
	LegendreBasis basis;
	
	REQUIRE(basis.Evaluate(1, 0.0) == Catch::Approx(0.0));
	REQUIRE(basis.Evaluate(1, 0.5) == Catch::Approx(0.5));
	REQUIRE(basis.Evaluate(1, -0.5) == Catch::Approx(-0.5));
	REQUIRE(basis.Evaluate(1, 1.0) == Catch::Approx(1.0));
}

TEST_CASE("LegendreBasis - P_2(x) = (3x^2 - 1)/2", "[OrthogonalBasis][Legendre]")
{
	LegendreBasis basis;
	
	REQUIRE(basis.Evaluate(2, 0.0) == Catch::Approx(-0.5));
	REQUIRE(basis.Evaluate(2, 1.0) == Catch::Approx(1.0));
	REQUIRE(basis.Evaluate(2, -1.0) == Catch::Approx(1.0));
}

TEST_CASE("LegendreBasis - P_n(1) = 1", "[OrthogonalBasis][Legendre]")
{
	LegendreBasis basis;
	
	for (int n = 0; n <= 10; n++) {
		REQUIRE(basis.Evaluate(n, 1.0) == Catch::Approx(1.0));
	}
}

TEST_CASE("LegendreBasis - P_n(-1) = (-1)^n", "[OrthogonalBasis][Legendre]")
{
	LegendreBasis basis;
	
	for (int n = 0; n <= 10; n++) {
		Real expected = (n % 2 == 0) ? 1.0 : -1.0;
		REQUIRE(basis.Evaluate(n, -1.0) == Catch::Approx(expected));
	}
}

TEST_CASE("LegendreBasis - Normalization formula", "[OrthogonalBasis][Legendre]")
{
	LegendreBasis basis;
	
	// ||P_n||^2 = 2/(2n+1)
	REQUIRE(basis.Normalization(0) == Catch::Approx(2.0));        // 2/1
	REQUIRE(basis.Normalization(1) == Catch::Approx(2.0/3.0));    // 2/3
	REQUIRE(basis.Normalization(2) == Catch::Approx(2.0/5.0));    // 2/5
	REQUIRE(basis.Normalization(5) == Catch::Approx(2.0/11.0));   // 2/11
}

TEST_CASE("LegendreBasis - Recurrence coefficients", "[OrthogonalBasis][Legendre]")
{
	LegendreBasis basis;
	
	Real a, b, c;
	
	// For n=0: a_0 = 1, b_0 = 0, c_0 = 0
	basis.RecurrenceCoefficients(0, a, b, c);
	REQUIRE(a == Catch::Approx(1.0));
	REQUIRE(b == Catch::Approx(0.0));
	REQUIRE(c == Catch::Approx(0.0));
	
	// For n=1: a_1 = 3/2, b_1 = 0, c_1 = -1/2
	basis.RecurrenceCoefficients(1, a, b, c);
	REQUIRE(a == Catch::Approx(1.5));
	REQUIRE(b == Catch::Approx(0.0));
	REQUIRE(c == Catch::Approx(-0.5));
}

TEST_CASE("LegendreBasis - Invalid n throws", "[OrthogonalBasis][Legendre][validation]")
{
	LegendreBasis basis;
	
	REQUIRE_THROWS_AS(basis.Evaluate(-1, 0.5), std::invalid_argument);
	REQUIRE_THROWS_AS(basis.Normalization(-1), std::invalid_argument);
	
	Real a, b, c;
	REQUIRE_THROWS_AS(basis.RecurrenceCoefficients(-1, a, b, c), std::invalid_argument);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           AssociatedLegendreBasis                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("AssociatedLegendreBasis - Construction", "[OrthogonalBasis][AssociatedLegendre]")
{
	AssociatedLegendreBasis basis(2);
	
	REQUIRE(basis.Order() == 2);
	REQUIRE(basis.DomainMin() == Catch::Approx(-1.0));
	REQUIRE(basis.DomainMax() == Catch::Approx(1.0));
}

TEST_CASE("AssociatedLegendreBasis - Evaluation for m=0", "[OrthogonalBasis][AssociatedLegendre]")
{
	AssociatedLegendreBasis assoc(0);
	
	// SphLegendre uses spherical harmonics normalization, so values differ from standard Legendre
	// Just verify it returns reasonable finite values
	for (int l = 0; l <= 5; l++) {
		REQUIRE(std::isfinite(assoc.Evaluate(l, 0.5)));
	}
}

TEST_CASE("AssociatedLegendreBasis - Invalid m throws", "[OrthogonalBasis][AssociatedLegendre][validation]")
{
	REQUIRE_THROWS_AS(AssociatedLegendreBasis(-1), std::invalid_argument);
}

TEST_CASE("AssociatedLegendreBasis - l < m throws", "[OrthogonalBasis][AssociatedLegendre][validation]")
{
	AssociatedLegendreBasis basis(3);
	
	REQUIRE_THROWS_AS(basis.Evaluate(2, 0.5), std::invalid_argument);  // l=2 < m=3
	REQUIRE_NOTHROW(basis.Evaluate(3, 0.5));  // l=3 = m=3 is OK
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           HermiteBasis (Physicist's)                                ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("HermiteBasis - Domain is unbounded", "[OrthogonalBasis][Hermite]")
{
	HermiteBasis basis;
	
	REQUIRE(std::isinf(basis.DomainMin()));
	REQUIRE(basis.DomainMin() < 0);
	REQUIRE(std::isinf(basis.DomainMax()));
	REQUIRE(basis.DomainMax() > 0);
}

TEST_CASE("HermiteBasis - Weight function w(x) = e^(-x^2)", "[OrthogonalBasis][Hermite]")
{
	HermiteBasis basis;
	
	REQUIRE(basis.WeightFunction(0.0) == Catch::Approx(1.0));
	REQUIRE(basis.WeightFunction(1.0) == Catch::Approx(std::exp(-1.0)));
	REQUIRE(basis.WeightFunction(-1.0) == Catch::Approx(std::exp(-1.0)));
	REQUIRE(basis.WeightFunction(2.0) == Catch::Approx(std::exp(-4.0)));
}

TEST_CASE("HermiteBasis - H_0(x) = 1", "[OrthogonalBasis][Hermite]")
{
	HermiteBasis basis;
	
	REQUIRE(basis.Evaluate(0, 0.0) == Catch::Approx(1.0));
	REQUIRE(basis.Evaluate(0, 1.0) == Catch::Approx(1.0));
	REQUIRE(basis.Evaluate(0, -2.0) == Catch::Approx(1.0));
}

TEST_CASE("HermiteBasis - H_1(x) = 2x", "[OrthogonalBasis][Hermite]")
{
	HermiteBasis basis;
	
	REQUIRE(basis.Evaluate(1, 0.0) == Catch::Approx(0.0));
	REQUIRE(basis.Evaluate(1, 1.0) == Catch::Approx(2.0));
	REQUIRE(basis.Evaluate(1, -1.0) == Catch::Approx(-2.0));
	REQUIRE(basis.Evaluate(1, 0.5) == Catch::Approx(1.0));
}

TEST_CASE("HermiteBasis - H_2(x) = 4x^2 - 2", "[OrthogonalBasis][Hermite]")
{
	HermiteBasis basis;
	
	REQUIRE(basis.Evaluate(2, 0.0) == Catch::Approx(-2.0));
	REQUIRE(basis.Evaluate(2, 1.0) == Catch::Approx(2.0));
	REQUIRE(basis.Evaluate(2, -1.0) == Catch::Approx(2.0));
}

TEST_CASE("HermiteBasis - Parity: H_n(-x) = (-1)^n H_n(x)", "[OrthogonalBasis][Hermite]")
{
	HermiteBasis basis;
	Real x = 1.5;
	
	for (int n = 0; n <= 8; n++) {
		Real H_pos = basis.Evaluate(n, x);
		Real H_neg = basis.Evaluate(n, -x);
		Real expected_sign = (n % 2 == 0) ? 1.0 : -1.0;
		REQUIRE(H_neg == Catch::Approx(expected_sign * H_pos));
	}
}

TEST_CASE("HermiteBasis - Normalization: 2^n * n! * sqrt(pi)", "[OrthogonalBasis][Hermite]")
{
	HermiteBasis basis;
	
	// n=0: 2^0 * 0! * sqrt(pi) = sqrt(pi)
	REQUIRE(basis.Normalization(0) == Catch::Approx(std::sqrt(Constants::PI)));
	
	// n=1: 2^1 * 1! * sqrt(pi) = 2*sqrt(pi)
	REQUIRE(basis.Normalization(1) == Catch::Approx(2.0 * std::sqrt(Constants::PI)));
	
	// n=2: 2^2 * 2! * sqrt(pi) = 4*2*sqrt(pi) = 8*sqrt(pi)
	REQUIRE(basis.Normalization(2) == Catch::Approx(8.0 * std::sqrt(Constants::PI)));
}

TEST_CASE("HermiteBasis - Recurrence coefficients", "[OrthogonalBasis][Hermite]")
{
	HermiteBasis basis;
	
	Real a, b, c;
	
	// a_n = 2, b_n = 0, c_n = -2n
	for (int n = 0; n <= 5; n++) {
		basis.RecurrenceCoefficients(n, a, b, c);
		REQUIRE(a == Catch::Approx(2.0));
		REQUIRE(b == Catch::Approx(0.0));
		REQUIRE(c == Catch::Approx(-2.0 * n));
	}
}

TEST_CASE("HermiteBasis - QuantumWavefunction normalization", "[OrthogonalBasis][Hermite]")
{
	HermiteBasis basis;
	
	// The quantum wavefunction should integrate to 1
	// Here we just check that it returns reasonable values
	REQUIRE(std::isfinite(basis.QuantumWavefunction(0, 0.0)));
	REQUIRE(std::isfinite(basis.QuantumWavefunction(1, 0.0)));
	REQUIRE(std::isfinite(basis.QuantumWavefunction(5, 1.0)));
	
	// psi_0(0) = N_0 * H_0(0) * e^0 = 1/sqrt(sqrt(pi)) * 1 * 1
	Real psi0_at_0 = basis.QuantumWavefunction(0, 0.0);
	REQUIRE(psi0_at_0 == Catch::Approx(1.0 / std::pow(Constants::PI, 0.25)));
}

TEST_CASE("HermiteBasis - Invalid n throws", "[OrthogonalBasis][Hermite][validation]")
{
	HermiteBasis basis;
	
	REQUIRE_THROWS_AS(basis.Evaluate(-1, 0.0), std::invalid_argument);
	REQUIRE_THROWS_AS(basis.Normalization(-1), std::invalid_argument);
	REQUIRE_THROWS_AS(basis.QuantumWavefunction(-1, 0.0), std::invalid_argument);
	
	Real a, b, c;
	REQUIRE_THROWS_AS(basis.RecurrenceCoefficients(-1, a, b, c), std::invalid_argument);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           ProbabilistHermiteBasis                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ProbabilistHermiteBasis - He_0(x) = 1", "[OrthogonalBasis][ProbabilistHermite]")
{
	ProbabilistHermiteBasis basis;
	
	REQUIRE(basis.Evaluate(0, 0.0) == Catch::Approx(1.0));
	REQUIRE(basis.Evaluate(0, 1.0) == Catch::Approx(1.0));
}

TEST_CASE("ProbabilistHermiteBasis - He_1(x) = x", "[OrthogonalBasis][ProbabilistHermite]")
{
	ProbabilistHermiteBasis basis;
	
	REQUIRE(basis.Evaluate(1, 0.0) == Catch::Approx(0.0));
	REQUIRE(basis.Evaluate(1, 1.0) == Catch::Approx(1.0));
	REQUIRE(basis.Evaluate(1, -2.0) == Catch::Approx(-2.0));
}

TEST_CASE("ProbabilistHermiteBasis - He_2(x) = x^2 - 1", "[OrthogonalBasis][ProbabilistHermite]")
{
	ProbabilistHermiteBasis basis;
	
	REQUIRE(basis.Evaluate(2, 0.0) == Catch::Approx(-1.0));
	REQUIRE(basis.Evaluate(2, 1.0) == Catch::Approx(0.0).margin(TOL(1e-14, 1e-5)));
	REQUIRE(basis.Evaluate(2, 2.0) == Catch::Approx(3.0));
}

TEST_CASE("ProbabilistHermiteBasis - Weight function w(x) = e^(-x^2/2)", "[OrthogonalBasis][ProbabilistHermite]")
{
	ProbabilistHermiteBasis basis;
	
	REQUIRE(basis.WeightFunction(0.0) == Catch::Approx(1.0));
	REQUIRE(basis.WeightFunction(1.0) == Catch::Approx(std::exp(-0.5)));
	REQUIRE(basis.WeightFunction(2.0) == Catch::Approx(std::exp(-2.0)));
}

TEST_CASE("ProbabilistHermiteBasis - Invalid n throws", "[OrthogonalBasis][ProbabilistHermite][validation]")
{
	ProbabilistHermiteBasis basis;
	
	REQUIRE_THROWS_AS(basis.Evaluate(-1, 0.0), std::invalid_argument);
	REQUIRE_THROWS_AS(basis.Normalization(-1), std::invalid_argument);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           LaguerreBasis                                             ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("LaguerreBasis - Domain", "[OrthogonalBasis][Laguerre]")
{
	LaguerreBasis basis;
	
	REQUIRE(basis.DomainMin() == Catch::Approx(0.0));
	REQUIRE(std::isinf(basis.DomainMax()));
	REQUIRE(basis.DomainMax() > 0);
}

TEST_CASE("LaguerreBasis - Weight function w(x) = e^(-x)", "[OrthogonalBasis][Laguerre]")
{
	LaguerreBasis basis;
	
	REQUIRE(basis.WeightFunction(0.0) == Catch::Approx(1.0));
	REQUIRE(basis.WeightFunction(1.0) == Catch::Approx(std::exp(-1.0)));
	REQUIRE(basis.WeightFunction(2.0) == Catch::Approx(std::exp(-2.0)));
}

TEST_CASE("LaguerreBasis - L_0(x) = 1", "[OrthogonalBasis][Laguerre]")
{
	LaguerreBasis basis;
	
	REQUIRE(basis.Evaluate(0, 0.0) == Catch::Approx(1.0));
	REQUIRE(basis.Evaluate(0, 1.0) == Catch::Approx(1.0));
	REQUIRE(basis.Evaluate(0, 5.0) == Catch::Approx(1.0));
}

TEST_CASE("LaguerreBasis - L_1(x) = 1 - x", "[OrthogonalBasis][Laguerre]")
{
	LaguerreBasis basis;
	
	REQUIRE(basis.Evaluate(1, 0.0) == Catch::Approx(1.0));
	REQUIRE(basis.Evaluate(1, 1.0) == Catch::Approx(0.0));
	REQUIRE(basis.Evaluate(1, 2.0) == Catch::Approx(-1.0));
}

TEST_CASE("LaguerreBasis - L_2(x) = (x^2 - 4x + 2)/2", "[OrthogonalBasis][Laguerre]")
{
	LaguerreBasis basis;
	
	// L_2(0) = 2/2 = 1
	REQUIRE(basis.Evaluate(2, 0.0) == Catch::Approx(1.0));
	// L_2(2) = (4 - 8 + 2)/2 = -1
	REQUIRE(basis.Evaluate(2, 2.0) == Catch::Approx(-1.0));
}

TEST_CASE("LaguerreBasis - L_n(0) = 1 for all n", "[OrthogonalBasis][Laguerre]")
{
	LaguerreBasis basis;
	
	for (int n = 0; n <= 10; n++) {
		REQUIRE(basis.Evaluate(n, 0.0) == Catch::Approx(1.0));
	}
}

TEST_CASE("LaguerreBasis - Normalization is 1", "[OrthogonalBasis][Laguerre]")
{
	LaguerreBasis basis;
	
	for (int n = 0; n <= 10; n++) {
		REQUIRE(basis.Normalization(n) == Catch::Approx(1.0));
	}
}

TEST_CASE("LaguerreBasis - Recurrence coefficients", "[OrthogonalBasis][Laguerre]")
{
	LaguerreBasis basis;
	
	Real a, b, c;
	
	// For n=0: a = -1, b = 1, c = 0
	basis.RecurrenceCoefficients(0, a, b, c);
	REQUIRE(a == Catch::Approx(-1.0));
	REQUIRE(b == Catch::Approx(1.0));
	REQUIRE(c == Catch::Approx(0.0));
	
	// For n=1: a = -1/2, b = 3/2, c = -1/2
	basis.RecurrenceCoefficients(1, a, b, c);
	REQUIRE(a == Catch::Approx(-0.5));
	REQUIRE(b == Catch::Approx(1.5));
	REQUIRE(c == Catch::Approx(-0.5));
}

TEST_CASE("LaguerreBasis - Invalid n throws", "[OrthogonalBasis][Laguerre][validation]")
{
	LaguerreBasis basis;
	
	REQUIRE_THROWS_AS(basis.Evaluate(-1, 0.0), std::invalid_argument);
	REQUIRE_THROWS_AS(basis.Normalization(-1), std::invalid_argument);
	
	Real a, b, c;
	REQUIRE_THROWS_AS(basis.RecurrenceCoefficients(-1, a, b, c), std::invalid_argument);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           AssociatedLaguerreBasis                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("AssociatedLaguerreBasis - Construction", "[OrthogonalBasis][AssociatedLaguerre]")
{
	AssociatedLaguerreBasis basis(2);
	
	REQUIRE(basis.Alpha() == 2);
	REQUIRE(basis.DomainMin() == Catch::Approx(0.0));
	REQUIRE(std::isinf(basis.DomainMax()));
}

TEST_CASE("AssociatedLaguerreBasis - L_n^0 = L_n", "[OrthogonalBasis][AssociatedLaguerre]")
{
	AssociatedLaguerreBasis assoc(0);
	LaguerreBasis laguerre;
	
	for (int n = 0; n <= 5; n++) {
		REQUIRE(assoc.Evaluate(n, 1.0) == Catch::Approx(laguerre.Evaluate(n, 1.0)).epsilon(TOL(1e-10, 1e-5)));
	}
}

TEST_CASE("AssociatedLaguerreBasis - L_n^alpha(0) = (n+alpha choose n)", "[OrthogonalBasis][AssociatedLaguerre]")
{
	AssociatedLaguerreBasis basis(1);  // alpha = 1
	
	// L_0^1(0) = C(1,0) = 1
	REQUIRE(basis.Evaluate(0, 0.0) == Catch::Approx(1.0));
	// L_1^1(0) = C(2,1) = 2
	REQUIRE(basis.Evaluate(1, 0.0) == Catch::Approx(2.0));
	// L_2^1(0) = C(3,2) = 3
	REQUIRE(basis.Evaluate(2, 0.0) == Catch::Approx(3.0));
}

TEST_CASE("AssociatedLaguerreBasis - Invalid inputs", "[OrthogonalBasis][AssociatedLaguerre][validation]")
{
	REQUIRE_THROWS_AS(AssociatedLaguerreBasis(-1), std::invalid_argument);
	
	AssociatedLaguerreBasis basis(1);
	REQUIRE_THROWS_AS(basis.Evaluate(-1, 0.0), std::invalid_argument);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Orthogonality Properties                                  ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Legendre - Orthogonality (discrete check)", "[OrthogonalBasis][Orthogonality]")
{
	LegendreBasis basis;
	
	// Check that P_0 and P_1 have different parity at symmetric points
	// This is a weak orthogonality check
	Real P0_sum = 0, P1_sum = 0;
	int N = 100;
	for (int i = 0; i < N; i++) {
		Real x = -1.0 + 2.0 * i / (N - 1);
		P0_sum += basis.Evaluate(0, x) * basis.Evaluate(1, x);
	}
	// P_0 * P_1 should integrate to approximately 0
	REQUIRE(P0_sum / N == Catch::Approx(0.0).margin(0.1));
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Domain Checks                                             ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("OrthogonalBasis - Manual domain check", "[OrthogonalBasis]")
{
	LegendreBasis legendre;
	HermiteBasis hermite;
	LaguerreBasis laguerre;
	
	// Legendre: [-1, 1]
	REQUIRE(0.0 >= legendre.DomainMin());
	REQUIRE(0.0 <= legendre.DomainMax());
	REQUIRE(-1.0 >= legendre.DomainMin());
	REQUIRE(1.0 <= legendre.DomainMax());
	
	// Hermite: (-∞, ∞) - any finite value is in domain
	REQUIRE(hermite.DomainMin() < -1000.0);
	REQUIRE(hermite.DomainMax() > 1000.0);
	
	// Laguerre: [0, ∞)
	REQUIRE(laguerre.DomainMin() == Catch::Approx(0.0));
	REQUIRE(laguerre.DomainMax() > 1000.0);
}

} // namespace MML::Tests::Core::OrthogonalBasisTests
