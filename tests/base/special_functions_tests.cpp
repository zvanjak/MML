///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML) Tests                            ///
///                                                                                   ///
///  File:        special_functions_tests.cpp                                         ///
///  Description: Tests for SpecialFunctions.h                                         ///
///               Incomplete gamma/beta (+inverses), digamma/trigamma/polygamma,       ///
///               erf/erfc inverses                                                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#include "../TestPrecision.h"
#include "../../mml/base/SpecialFunctions.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <cmath>

using namespace MML;
using namespace MML::SpecialFunctions;
using Catch::Approx;

namespace MML::Tests::Base::SpecialFunctionsTests {

TEST_CASE("SpecialFunctions::Beta_and_Gamma", "[special][gamma][beta]") {
	TEST_PRECISION_INFO();
	REQUIRE(Beta(2.0, 3.0) == Approx(1.0 / 12.0).epsilon(TOL(1e-12, 1e-6)));
	REQUIRE(LnBeta(2.0, 3.0) == Approx(std::log(1.0 / 12.0)).epsilon(TOL(1e-12, 1e-6)));
	REQUIRE(Gamma(5.0) == Approx(24.0).epsilon(TOL(1e-12, 1e-6)));
}

TEST_CASE("SpecialFunctions::RegularizedGammaP", "[special][gamma]") {
	TEST_PRECISION_INFO();
	// P(1,x) = 1 - e^{-x}
	REQUIRE(RegularizedGammaP(1.0, 1.0) == Approx(1.0 - std::exp(-1.0)).epsilon(TOL(1e-10, 2e-6)));
	// P(2,2) = 1 - 3 e^{-2}
	REQUIRE(RegularizedGammaP(2.0, 2.0) == Approx(1.0 - 3.0 * std::exp(-2.0)).epsilon(TOL(1e-10, 2e-6)));
	// P(1/2, x) = erf(sqrt(x))
	REQUIRE(RegularizedGammaP(0.5, 0.5) == Approx(std::erf(std::sqrt(0.5))).epsilon(TOL(1e-10, 2e-6)));
	// Q = 1 - P
	REQUIRE(RegularizedGammaQ(2.0, 2.0) == Approx(3.0 * std::exp(-2.0)).epsilon(TOL(1e-10, 2e-6)));
	// Unregularized lower: gamma(1,x) = 1 - e^{-x}
	REQUIRE(IncompleteGammaLower(1.0, 1.0) == Approx(1.0 - std::exp(-1.0)).epsilon(TOL(1e-10, 2e-6)));
	REQUIRE(RegularizedGammaP(3.0, 0.0) == Approx(0.0).margin(1e-14));
}

TEST_CASE("SpecialFunctions::RegularizedBetaI", "[special][beta]") {
	TEST_PRECISION_INFO();
	REQUIRE(RegularizedBetaI(0.5, 2.0, 2.0) == Approx(0.5).epsilon(TOL(1e-10, 2e-6)));          // symmetric
	REQUIRE(RegularizedBetaI(0.3, 2.0, 3.0) == Approx(0.3483).epsilon(TOL(1e-6, 2e-6)));        // closed-form CDF
	REQUIRE(RegularizedBetaI(0.0, 2.0, 3.0) == Approx(0.0).margin(1e-14));
	REQUIRE(RegularizedBetaI(1.0, 2.0, 3.0) == Approx(1.0).epsilon(TOL(1e-12, 1e-6)));
	// I_x(1,1) = x (uniform CDF)
	REQUIRE(RegularizedBetaI(0.37, 1.0, 1.0) == Approx(0.37).epsilon(TOL(1e-10, 2e-6)));
}

TEST_CASE("SpecialFunctions::Digamma", "[special][digamma]") {
	TEST_PRECISION_INFO();
	const double euler = 0.5772156649015329;
	REQUIRE(Digamma(1.0) == Approx(-euler).epsilon(TOL(1e-10, 2e-6)));
	REQUIRE(Digamma(2.0) == Approx(1.0 - euler).epsilon(TOL(1e-10, 2e-6)));
	REQUIRE(Digamma(0.5) == Approx(-euler - 2.0 * std::log(2.0)).epsilon(TOL(1e-10, 2e-6)));
	// Recurrence check: psi(x+1) - psi(x) = 1/x
	REQUIRE(Digamma(3.7) - Digamma(2.7) == Approx(1.0 / 2.7).epsilon(TOL(1e-9, 2e-6)));
	REQUIRE_THROWS_AS(Digamma(0.0), DomainError);
	REQUIRE_THROWS_AS(Digamma(-2.0), DomainError);
}

TEST_CASE("SpecialFunctions::Trigamma_and_Polygamma", "[special][polygamma]") {
	TEST_PRECISION_INFO();
	const double pi = Constants::PI;
	REQUIRE(Trigamma(1.0) == Approx(pi * pi / 6.0).epsilon(TOL(1e-9, 2e-6)));
	REQUIRE(Trigamma(0.5) == Approx(pi * pi / 2.0).epsilon(TOL(1e-9, 2e-6)));
	// Polygamma(1,x) == Trigamma(x)
	REQUIRE(Polygamma(1, 2.3) == Approx(Trigamma(2.3)).epsilon(TOL(1e-9, 2e-6)));
	// psi''(1) = -2 zeta(3)
	REQUIRE(Polygamma(2, 1.0) == Approx(-2.404113806319188).epsilon(TOL(1e-8, 5e-6)));
	// psi'''(1) = pi^4/15
	REQUIRE(Polygamma(3, 1.0) == Approx(std::pow(pi, 4) / 15.0).epsilon(TOL(1e-8, 5e-6)));
}

TEST_CASE("SpecialFunctions::ErfInv_ErfcInv", "[special][erf]") {
	TEST_PRECISION_INFO();
	REQUIRE(ErfInv(0.0) == Approx(0.0).margin(1e-14));
	REQUIRE(ErfInv(0.5) == Approx(0.4769362762044699).epsilon(TOL(1e-10, 2e-6)));
	// Round trip: erf(erfinv(x)) = x
	for (double x : { -0.9, -0.3, 0.1, 0.75, 0.99 })
		REQUIRE(std::erf(ErfInv(x)) == Approx(x).epsilon(TOL(1e-10, 3e-6)).margin(TOL(1e-12, 1e-7)));
	// erfcinv(erfc(z)) = z
	REQUIRE(ErfcInv(std::erfc(1.0)) == Approx(1.0).epsilon(TOL(1e-9, 2e-6)));
	REQUIRE(ErfcInv(1.0) == Approx(0.0).margin(1e-10));
}

TEST_CASE("SpecialFunctions::RegularizedGammaPInv_roundtrip", "[special][gamma][inverse]") {
	TEST_PRECISION_INFO();
	for (double a : { 0.5, 1.0, 2.5, 7.0 }) {
		for (double x : { 0.2, 1.0, 3.0, 9.0 }) {
			double p = RegularizedGammaP(a, x);
			CAPTURE(a, x, p);
			REQUIRE(RegularizedGammaPInv(a, p) == Approx(x).epsilon(TOL(1e-7, 2e-4)));
		}
	}
}

TEST_CASE("SpecialFunctions::RegularizedBetaIInv_roundtrip", "[special][beta][inverse]") {
	TEST_PRECISION_INFO();
	for (double a : { 0.7, 1.0, 2.0, 5.0 }) {
		for (double b : { 0.9, 1.0, 3.0, 8.0 }) {
			for (double x : { 0.15, 0.4, 0.8 }) {
				double p = RegularizedBetaI(x, a, b);
				CAPTURE(a, b, x, p);
				REQUIRE(RegularizedBetaIInv(p, a, b) == Approx(x).epsilon(TOL(1e-7, 2e-3)).margin(TOL(1e-10, 1e-6)));
			}
		}
	}
}

TEST_CASE("SpecialFunctions::ExponentialIntegrals", "[special][expint]") {
	TEST_PRECISION_INFO();
	REQUIRE(ExpIntegralEi(1.0) == Approx(1.8951178163559368).epsilon(TOL(1e-9, 3e-6)));
	REQUIRE(ExpIntegralEi(2.0) == Approx(4.954234356001891).epsilon(TOL(1e-9, 3e-6)));
	REQUIRE(ExpIntegralEi(-1.0) == Approx(-0.21938393439552029).epsilon(TOL(1e-9, 3e-6))); // = -E1(1)
	REQUIRE(ExpIntegralE1(1.0) == Approx(0.21938393439552027).epsilon(TOL(1e-9, 3e-6)));
	REQUIRE(ExpIntegralEn(2, 1.0) == Approx(0.14849550677592205).epsilon(TOL(1e-9, 3e-6)));
	REQUIRE(ExpIntegralEn(0, 2.0) == Approx(std::exp(-2.0) / 2.0).epsilon(TOL(1e-12, 1e-6))); // E0(x)=e^-x/x
	REQUIRE_THROWS_AS(ExpIntegralEi(0.0), DomainError);
	REQUIRE_THROWS_AS(ExpIntegralEn(1, 0.0), DomainError);
}

TEST_CASE("SpecialFunctions::SineCosineIntegrals", "[special][sici]") {
	TEST_PRECISION_INFO();
	REQUIRE(SineIntegral(1.0) == Approx(0.9460830703671830).epsilon(TOL(1e-9, 3e-6)));
	REQUIRE(SineIntegral(10.0) == Approx(1.6583475942188740).epsilon(TOL(1e-9, 3e-6)));
	REQUIRE(CosineIntegral(1.0) == Approx(0.3374039229009681).epsilon(TOL(1e-9, 3e-6)));
	REQUIRE(CosineIntegral(10.0) == Approx(-0.04545643300445537).epsilon(TOL(1e-7, 5e-6)).margin(TOL(1e-9, 1e-6)));
	REQUIRE(SineIntegral(0.0) == Approx(0.0).margin(1e-14));
	REQUIRE(SineIntegral(-1.0) == Approx(-0.9460830703671830).epsilon(TOL(1e-9, 3e-6))); // Si is odd
	// Si(+inf) = pi/2
	REQUIRE(SineIntegral(1000.0) == Approx(Constants::PI / 2.0).epsilon(1e-3));
}

TEST_CASE("SpecialFunctions::Dawson", "[special][dawson]") {
	TEST_PRECISION_INFO();
	REQUIRE(Dawson(0.1) == Approx(0.09933599239785287).epsilon(1e-6)); // small-x branch
	REQUIRE(Dawson(1.0) == Approx(0.5380795069127684).epsilon(1e-6));
	REQUIRE(Dawson(2.0) == Approx(0.30134038892379196).epsilon(1e-6));
	REQUIRE(Dawson(-1.0) == Approx(-0.5380795069127684).epsilon(1e-6)); // odd
	REQUIRE(Dawson(0.0) == Approx(0.0).margin(1e-14));
}

TEST_CASE("SpecialFunctions::Fresnel", "[special][fresnel]") {
	TEST_PRECISION_INFO();
	REQUIRE(FresnelC(1.0) == Approx(0.7798934003768228).epsilon(TOL(1e-9, 3e-6)));
	REQUIRE(FresnelS(1.0) == Approx(0.4382591473903548).epsilon(TOL(1e-9, 3e-6)));
	REQUIRE(FresnelC(0.5) == Approx(0.4923442258714463).epsilon(TOL(1e-9, 3e-6))); // series branch
	REQUIRE(FresnelC(2.0) == Approx(0.4882534060753407).epsilon(TOL(1e-9, 3e-6))); // CF branch
	REQUIRE(FresnelS(2.0) == Approx(0.3434156783636982).epsilon(TOL(1e-9, 3e-6)));
	REQUIRE(FresnelC(-1.0) == Approx(-0.7798934003768228).epsilon(TOL(1e-9, 3e-6))); // odd
	REQUIRE(FresnelS(0.0) == Approx(0.0).margin(1e-14));
	// C(+inf) = S(+inf) = 1/2
	REQUIRE(FresnelC(30.0) == Approx(0.5).margin(0.02));
	REQUIRE(FresnelS(30.0) == Approx(0.5).margin(0.02));
}

TEST_CASE("SpecialFunctions::LambertW0", "[special][lambertw]") {
	TEST_PRECISION_INFO();
	REQUIRE(LambertW0(0.0) == Approx(0.0).margin(1e-14));
	REQUIRE(LambertW0(1.0) == Approx(0.5671432904097838).epsilon(TOL(1e-10, 2e-6))); // Omega constant
	REQUIRE(LambertW0(std::exp(1.0)) == Approx(1.0).epsilon(TOL(1e-10, 2e-6)));       // W0(e) = 1
	REQUIRE(LambertW0(10.0) == Approx(1.7455280027406994).epsilon(TOL(1e-10, 2e-6)));
	REQUIRE(LambertW0(-0.2) == Approx(-0.2591711018190738).epsilon(TOL(1e-9, 2e-6)));
	REQUIRE(LambertW0(-0.36787944117144232) == Approx(-1.0).margin(1e-6)); // branch point -1/e
	// Defining identity W(x) e^{W(x)} = x
	for (double x : { -0.3, -0.1, 0.5, 2.0, 5.0, 50.0 }) {
		double w = LambertW0(x);
		REQUIRE(w * std::exp(w) == Approx(x).epsilon(TOL(1e-9, 3e-6)).margin(TOL(1e-12, 1e-6)));
	}
	REQUIRE_THROWS_AS(LambertW0(-0.5), DomainError);
}

TEST_CASE("SpecialFunctions::LambertWm1", "[special][lambertw]") {
	TEST_PRECISION_INFO();
	REQUIRE(LambertWm1(-0.1) == Approx(-3.577152063957297).epsilon(TOL(1e-9, 2e-6)));
	REQUIRE(LambertWm1(-0.2) == Approx(-2.542641357773526).epsilon(TOL(1e-9, 2e-6)));
	REQUIRE(LambertWm1(-0.3) == Approx(-1.781337023421627).epsilon(TOL(1e-9, 2e-6)));
	for (double x : { -0.35, -0.2, -0.05, -0.001 }) {
		double w = LambertWm1(x);
		REQUIRE(w * std::exp(w) == Approx(x).epsilon(TOL(1e-8, 3e-6)).margin(TOL(1e-12, 1e-6)));
	}
	REQUIRE_THROWS_AS(LambertWm1(0.1), DomainError);
	REQUIRE_THROWS_AS(LambertWm1(-0.5), DomainError);
}

TEST_CASE("SpecialFunctions::Airy_values", "[special][airy]") {
	TEST_PRECISION_INFO();
	REQUIRE(AiryAi(0.0) == Approx(0.3550280538878172).epsilon(TOL(1e-10, 5e-6)));
	REQUIRE(AiryAi(1.0) == Approx(0.13529241631288147).epsilon(TOL(1e-9, 5e-6)));
	REQUIRE(AiryAi(-1.0) == Approx(0.5355608832923521).epsilon(TOL(1e-9, 5e-6)));
	REQUIRE(AiryAi(2.0) == Approx(0.03492413042327235).epsilon(TOL(1e-9, 5e-6)));
	REQUIRE(AiryBi(0.0) == Approx(0.6149266274460007).epsilon(TOL(1e-10, 5e-6)));
	REQUIRE(AiryBi(1.0) == Approx(1.2074235949528713).epsilon(TOL(1e-9, 5e-6)));
	REQUIRE(AiryBi(-1.0) == Approx(0.10399738949694459).epsilon(TOL(1e-9, 5e-6)));
	REQUIRE(AiryAiPrime(0.0) == Approx(-0.2588194037928068).epsilon(TOL(1e-10, 5e-6)));
	REQUIRE(AiryBiPrime(0.0) == Approx(0.4482883573538264).epsilon(TOL(1e-10, 5e-6)));
}

TEST_CASE("SpecialFunctions::Airy_Wronskian_identity", "[special][airy]") {
	TEST_PRECISION_INFO();
	// Ai(x) Bi'(x) - Ai'(x) Bi(x) = 1/pi everywhere (validates series + both asymptotics)
	const double invPi = 1.0 / Constants::PI;
	for (double x : { -12.0, -8.0, -3.0, 0.0, 3.0, 8.0, 12.0 }) {
		CAPTURE(x);
		Real ai, bi, aip, bip;
		Airy(x, ai, bi, aip, bip);
		REQUIRE(ai * bip - aip * bi == Approx(invPi).epsilon(TOL(1e-7, 1e-4)).margin(TOL(1e-9, 1e-6)));
	}
}

TEST_CASE("SpecialFunctions::Airy_zeros", "[special][airy]") {
	TEST_PRECISION_INFO();
	REQUIRE(AiryAiZero(1) == Approx(-2.338107410459767).epsilon(TOL(1e-8, 3e-6)));
	REQUIRE(AiryAiZero(2) == Approx(-4.087949444130970).epsilon(TOL(1e-8, 3e-6)));
	REQUIRE(AiryAiZero(3) == Approx(-5.520559828095551).epsilon(TOL(1e-8, 3e-6)));
	REQUIRE(AiryAiZero(5) == Approx(-7.944133587120853).epsilon(TOL(1e-8, 3e-6)));
	REQUIRE(AiryAiZero(10) == Approx(-12.828776742769930).epsilon(TOL(1e-7, 3e-6))); // asymptotic region
	REQUIRE(AiryBiZero(1) == Approx(-1.173713222709128).epsilon(TOL(1e-8, 3e-6)));
	REQUIRE(AiryBiZero(2) == Approx(-3.271093302836357).epsilon(TOL(1e-8, 3e-6)));
	// Verify Ai actually vanishes at its computed zero
	REQUIRE(AiryAi(AiryAiZero(4)) == Approx(0.0).margin(TOL(1e-9, 1e-5)));
	REQUIRE_THROWS_AS(AiryAiZero(0), DomainError);
}

} // namespace MML::Tests::Base::SpecialFunctionsTests
