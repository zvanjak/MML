#include "../TestPrecision.h"
#include "../ProbabilityTailTestUtils.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/SpecialFunctions.h>
#endif
#include <catch2/catch_all.hpp>
#include <type_traits>

using namespace MML;
using namespace MML::SpecialFunctions;
using Catch::Matchers::WithinRel;
using MML::Testing::GammaLogRoundoff;
using MML::Testing::GammaTailTolerance;

namespace {
	void CheckTail(Real actual, Real expected, Real tolerance = Precision::EqualityTolerance) {
		REQUIRE(std::isfinite(expected));
		REQUIRE(expected > 0);
		INFO("actual/expected=" << actual / expected << ", relative error=" << std::abs(actual / expected - 1)
			<< ", log(expected)=" << std::log(expected));
		CHECK(std::isfinite(actual));
		CHECK(actual > 0);
		CHECK_THAT(actual, WithinRel(expected, tolerance));
	}
	Real BetaLogRoundoff(Real a, Real b, Real x) {
		return 4 * Precision::Eps * (1 + std::abs(a * std::log(x)) + std::abs(b * std::log1p(-x))
			+ std::abs(std::lgamma(a)) + std::abs(std::lgamma(b)) + std::abs(std::lgamma(a + b)));
	}
	Real InverseLogBudget(Real p, Real arithmetic) {
		const Real logTarget = std::log(std::min(p, 1 - p));
		return Precision::RoundoffZeroThreshold * (1 + std::abs(logTarget)) + arithmetic;
	}
}

TEST_CASE("Gamma lower series retains its working precision beyond the float range",
	"[special][forward-shape-contract]") {
	using S = SpecialFunctionsDetail::ShapeReal;
	const S x = 100;
	const S expected = x + std::log1p(-std::exp(-x)) - std::log(x);
	const S actual = SpecialFunctionsDetail::GammaLowerSeriesLog(S(1), x);
	REQUIRE(std::isfinite(actual));
	CHECK(std::abs(actual - expected) <= PrecisionValues<S>::RoundoffZeroThreshold * std::abs(expected));
}

TEST_CASE("Forward gamma resolves large central shapes and both tails", "[special][forward-shape-contract]") {
	struct Reference { Real a, x, lower, upper; };
	for (const auto& ref : {
		Reference{1000, 1000, Real(0.504205244180215508503777843602L), Real(0.495794755819784491496222156398L)},
		Reference{2000, 2000, Real(0.502973548444202534658120005029L), Real(0.497026451555797465341879994971L)},
		Reference{10000, 10000, Real(0.50132980833995520038274225129985137L), Real(0.49867019166004479961725774870014863L)},
		Reference{1000000, 1000000, Real(0.50013298076087259124432281750320671L), Real(0.49986701923912740875567718249679329L)},
		Reference{1000000, 998000, Real(0.022696114006736802806015042359720038L), Real(0.97730388599326319719398495764027996L)},
		Reference{1000000, 1002000, Real(0.97719590410123013724183986891560052L), Real(0.022804095898769862758160131084399481L)}
	}) {
		CAPTURE(ref.a, ref.x);
		CheckTail(RegularizedGammaP(ref.a, ref.x), ref.lower);
		CheckTail(RegularizedGammaQ(ref.a, ref.x), ref.upper);
	}
	if constexpr (!std::is_same_v<Real, float>) {
		CheckTail(RegularizedGammaP(Real(1e12), Real(1e12)), Real(0.50000013298076013381163142864942839L));
		CheckTail(RegularizedGammaQ(Real(1e12), Real(1000002000000.L)), Real(0.022750185939118724884581480465989594L));
	}
}

TEST_CASE("Forward shape scaling survives extreme finite normalization ranges", "[special][forward-shape-contract]") {
	const Real largest = std::numeric_limits<Real>::max();
	CheckTail(RegularizedGammaP(largest, largest), Real(0.5));
	CheckTail(RegularizedGammaQ(largest, largest), Real(0.5));
	CheckTail(RegularizedBetaI(16 / largest, 16, largest), RegularizedGammaP(16, 16));
	for (Real a : {std::nextafter(Real(Precision::ForwardShapeQuadratureSwitch), Real(0)),
		Real(Precision::ForwardShapeQuadratureSwitch),
		std::nextafter(Real(Precision::ForwardShapeQuadratureSwitch), largest)}) {
		CAPTURE(a);
		CheckTail(RegularizedGammaP(a, a), Real(0.504205244180215508503777843602L));
	}
}

TEST_CASE("Forward gamma preserves tiny shape complements", "[special][forward-shape-contract]") {
	const Real a = Real(1e-20);
	CheckTail(RegularizedGammaQ(a, 1), Real(2.193839343955202736794085258678966e-21L));
	CheckTail(IncompleteGammaUpper(a, 1), Real(0.21938393439552027367916377546012165L));
	CHECK(RegularizedGammaP(a, 1) == 1);
	CheckTail(RegularizedGammaQ(a, Real(0.5)), a * Real(0.5597735947761608117467959393150852352L));
	CHECK(RegularizedGammaQ(std::numeric_limits<Real>::denorm_min(), Real(0.5)) == std::numeric_limits<Real>::denorm_min());
}

TEST_CASE("Forward beta resolves central large shapes and tiny normalizations", "[special][forward-shape-contract]") {
	CHECK_THAT(RegularizedBetaI(Real(0.5), Real(1e12), Real(1e12)), WithinRel(Real(0.5), Precision::EqualityTolerance));
	if constexpr (!std::is_same_v<Real, float>) {
		const Real lower = RegularizedBetaI(Real(0.499999L), Real(1e12), Real(1e12));
		CHECK_THAT(lower, Catch::Matchers::WithinAbs(Real(0.0023388674905055493067805613442109366L), Precision::EqualityTolerance));
		CHECK_THAT(RegularizedBetaI(Real(0.500001L), Real(1e12), Real(1e12)),
			Catch::Matchers::WithinAbs(Real(0.99766113250949445069321943865578906L), Precision::EqualityTolerance));
	}
	CheckTail(RegularizedBetaI(Real(0.001), Real(0.5), 1000), Real(0.84275272573327598758675922467219181L));
	for (Real a : {Real(1e-20), TOL(Real(1e-300), Real(1e-30))})
		CheckTail(RegularizedBetaI(Real(0.25), a, 2 * a), Real(2) / 3);
}

TEST_CASE("Log beta avoids cancellation in disparate positive shapes", "[special][forward-shape-contract]") {
	CHECK_THAT(LnBeta(Real(0.5), 1000),
		Catch::Matchers::WithinAbs(Real(-2.8813876965715767707261080255075613L), Precision::RoundoffZeroThreshold));
	for (Real b : {Real(1000), Real(1e12), std::numeric_limits<Real>::max()}) {
		CAPTURE(b);
		CHECK_THAT(LnBeta(1, b), Catch::Matchers::WithinAbs(-std::log(b), Precision::RoundoffZeroThreshold * (1 + std::log(b))));
		CheckTail(Beta(1, b), 1 / b);
	}
}

TEST_CASE("Incomplete beta scales tiny probabilities without zero times infinity", "[special][forward-shape-contract]") {
	for (Real b : {std::numeric_limits<Real>::min(), std::numeric_limits<Real>::denorm_min()}) {
		CHECK(IncompleteBeta(0, 2, b) == 0);
		CHECK(IncompleteBeta(1, 2, b) == Beta(2, b));
		for (Real x : {Real(0.125), Real(0.5), Real(0.875)}) {
			CAPTURE(b, x);
			// Independent limits subtract up to three nearly equal terms.
			// Evaluate the reference before rounding to Real, especially for float.
			using ReferenceReal = std::common_type_t<Real, double>;
			const ReferenceReal refX = x, logarithm = -std::log1p(-refX);
			CheckTail(IncompleteBeta(x, 1, b), Real(logarithm));
			CheckTail(IncompleteBeta(x, 2, b), Real(logarithm - refX));
			CheckTail(IncompleteBeta(x, 3, b), Real(logarithm - refX - refX * refX / 2));
		}
	}
}

TEST_CASE("Large shape inverses consume the accurate forward kernels", "[special][forward-shape-contract]") {
	const Real a = Real(1000000);
	const Real median = RegularizedGammaPInv(a, Real(0.5));
	CHECK(median < a);
	CHECK(median > a - 1);
	// A large-shape median lies on the Real grid. Bound the probability change
	// by one x ulp times the central normal density (~1/sqrt(2*pi*a)).
	const Real medianSpacing = std::nextafter(median, a) - median;
	const Real medianBudget = Precision::EqualityTolerance + 2 * medianSpacing /
		std::sqrt(2 * Constants::PI * a);
	CHECK_THAT(RegularizedGammaP(a, median), WithinRel(Real(0.5), medianBudget));
	const Real betaShape = TOL(Real(1e12), Real(10000));
	CHECK_THAT(RegularizedBetaIInv(Real(0.5), betaShape, betaShape), WithinRel(Real(0.5), Precision::EqualityTolerance));
	for (Real p : {Real(0.1), Real(0.9)}) {
		const Real x = RegularizedBetaIInv(p, betaShape, betaShape);
		const Real densityScale = std::sqrt(betaShape);
		CHECK_THAT(RegularizedBetaI(x, betaShape, betaShape),
			Catch::Matchers::WithinAbs(p, Precision::EqualityTolerance + 4 * Precision::Eps * densityScale));
	}
}

TEST_CASE("Tiny beta shape inverse retains a resolvable complementary tail", "[special][forward-shape-contract]") {
	for (Real a : {std::numeric_limits<Real>::denorm_min(), Real(1e-20), Real(1e12)})
		CHECK(RegularizedBetaIInv(Real(0.5), a, a) == Real(0.5));
	if constexpr (!std::is_same_v<Real, float>) {
		const Real a = Real(1e-12), p = 1 - Real(1e-10);
		const Real expected = std::exp((std::log(p) - std::log1p(a)) / a);
		CheckTail(RegularizedBetaIInv(p, a, 2), expected);
	}
}

TEST_CASE("Closed-form shape tails remain probabilities at rounded endpoints", "[special][forward-shape-contract]") {
	for (Real x : {std::numeric_limits<Real>::denorm_min(), std::numeric_limits<Real>::min(),
		Real(0.125), Real(0.5), std::nextafter(Real(1), Real(0))}) {
		for (Real b : {Real(0.5), Real(1), Real(2.5), Real(1000)}) {
			CAPTURE(x, b);
			for (Real p : {RegularizedBetaI(x, 1, b), RegularizedBetaI(x, b, 1),
				RegularizedGammaP(1, x), RegularizedGammaQ(1, x)}) {
				CHECK(std::isfinite(p));
				CHECK(p >= 0);
				CHECK(p <= 1);
			}
		}
	}
}

TEST_CASE("Beta inverse retains representable extreme lower quantiles", "[special][beta-inverse-contract]") {
	const Real tolerance = std::max(Precision::EqualityTolerance, Precision::RoundoffZeroThreshold);
	for (Real p : {TOL(1e-100, 1e-20), TOL(1e-300, 1e-36)}) {
		CAPTURE(p);
		REQUIRE(p > 0);
		const Real root = std::sqrt(p);
		CheckTail(RegularizedBetaIInv(p, 2, 1), root, tolerance);
		CheckTail(RegularizedBetaIInv(p, 2, 2), root / std::sqrt(Real(3)), tolerance);
	}
	const Real p = TOL(1e-100, 1e-10);
	CheckTail(RegularizedBetaIInv(p, Real(0.5), Real(0.5)),
		Constants::PI * Constants::PI / 4 * TOL(1e-200, 1e-20), tolerance);
}

TEST_CASE("Beta inverse preserves exact uniform and stable unit shape identities", "[special][beta-inverse-contract]") {
	for (Real p : {std::numeric_limits<Real>::denorm_min(), std::numeric_limits<Real>::min(),
		Real(1e-300), Real(1e-100), Real(0.01), Real(0.5), Real(0.99), std::nextafter(Real(1), Real(0))}) {
		CAPTURE(p);
		CHECK(RegularizedBetaIInv(p, 1, 1) == p);
		for (Real b : {Real(0.1), Real(0.5), Real(2), Real(1000)}) {
			const Real expected = p == std::numeric_limits<Real>::denorm_min() && b == 2
				? p : -std::expm1(std::log1p(-p) / b);
			CHECK(RegularizedBetaIInv(p, 1, b) == expected);
		}
	}
	const Real p = std::nextafter(Real(1), Real(0));
	CHECK(RegularizedBetaIInv(p, 1, 2) == 1 - std::sqrt(1 - p));
	for (Real shape : {Real(2), Real(4), Real(16)}) {
		CHECK(RegularizedBetaIInv(1 - shape * Precision::Eps / 4, shape, 1) == std::nextafter(Real(1), Real(0)));
		CHECK(RegularizedBetaIInv(shape / 2 * std::numeric_limits<Real>::denorm_min(), 1, shape) ==
			std::numeric_limits<Real>::denorm_min());
	}
}

TEST_CASE("Beta inverse rounds genuine boundary range loss without an interior floor", "[special][beta-inverse-contract]") {
	const Real unit = std::numeric_limits<Real>::denorm_min();
	const Real rootUnit = std::sqrt(unit);
	for (Real p : {rootUnit / 8, rootUnit / 4}) {
		REQUIRE(p > 0);
		const Real x = RegularizedBetaIInv(p, Real(0.5), Real(0.5));
		CHECK(x == 0);
		CHECK_FALSE(std::signbit(x));
	}
	CHECK(RegularizedBetaIInv(rootUnit / 2, Real(0.5), Real(0.5)) == unit);
	CHECK(RegularizedBetaIInv(rootUnit, Real(0.5), Real(0.5)) == 2 * unit);
	CHECK(RegularizedBetaIInv(2 * rootUnit, Real(0.5), Real(0.5)) == 10 * unit);
	// sin^2(pi*p/2) is asymptotically (pi^2/4)*p^2, giving 39.478... subnormal units.
	CHECK(RegularizedBetaIInv(4 * rootUnit, Real(0.5), Real(0.5)) == 39 * unit);
	CHECK(RegularizedBetaIInv(Real(0.99), 2, Real(0.1)) == 1);
	CHECK(RegularizedBetaIInv(Real(0.99), Real(0.1), 2) < 1);
}

TEST_CASE("Beta inverse validates strict domains before exact endpoints", "[special][beta-inverse-contract]") {
	const Real infinity = std::numeric_limits<Real>::infinity();
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	for (Real invalid : {Real(0), Real(-1), infinity, -infinity, nan})
		for (Real p : {Real(0), Real(0.5), Real(1)}) {
			CHECK_THROWS_AS(RegularizedBetaIInv(p, invalid, 2), DomainError);
			CHECK_THROWS_AS(RegularizedBetaIInv(p, 2, invalid), DomainError);
		}
	for (Real p : {-std::numeric_limits<Real>::epsilon(), 1 + std::numeric_limits<Real>::epsilon(), infinity, -infinity, nan})
		CHECK_THROWS_AS(RegularizedBetaIInv(p, 2, 3), DomainError);
	for (Real a : {Real(0.1), Real(1), Real(2), Real(1000)}) {
		for (Real p : {Real(0), -Real(0)}) {
			const Real x = RegularizedBetaIInv(p, a, 2);
			CHECK(x == 0);
			CHECK_FALSE(std::signbit(x));
		}
		CHECK(RegularizedBetaIInv(1, a, 2) == 1);
	}
}

TEST_CASE("Beta inverse is monotone and preserves ordinary roundtrips and symmetry", "[special][beta-inverse-contract]") {
	for (Real a : {Real(0.1), Real(0.5), Real(2), Real(10), Real(1000)})
		for (Real b : {Real(0.1), Real(0.5), Real(2), Real(10), Real(1000)}) {
			Real previous = 0;
			for (Real p : {Real(1e-100), Real(1e-20), Real(0.01), Real(0.1), Real(0.5), Real(0.9), Real(0.99), std::nextafter(Real(1), Real(0))}) {
				CAPTURE(a, b, p);
				Real x = 0;
				REQUIRE_NOTHROW(x = RegularizedBetaIInv(p, a, b));
				CHECK(std::isfinite(x));
				CHECK(x >= previous);
				CHECK(x <= 1);
				if (x > 0 && x < 1 && p >= Real(0.01) && p <= Real(0.99)) {
					// Near one, quantization of x limits the accuracy of its forward probability.
					const Real propagatedUlp = std::numeric_limits<Real>::epsilon() /
						std::min(x, 1 - x);
					const Real logBudget = InverseLogBudget(p, 2 * BetaLogRoundoff(a, b, x));
					const Real inverseBudget = std::min(p, 1 - p) * std::expm1(logBudget);
					const Real budget = TOL(Precision::EqualityTolerance,
						std::max(Precision::EqualityTolerance, inverseBudget));
					CHECK(IsWithinAbs(RegularizedBetaI(x, a, b), p, budget + propagatedUlp));
				}
				previous = x;
			}
			for (Real p : {Real(0.125), Real(0.5), Real(0.875)}) {
				CAPTURE(a, b, p);
				const Real x = RegularizedBetaIInv(p, a, b);
				const Real reflected = 1 - RegularizedBetaIInv(1 - p, b, a);
				Real budget = Precision::RoundoffZeroThreshold;
				if (x > 0 && x < 1) {
					const Real logDensity = (a - 1) * std::log(x) + (b - 1) * std::log1p(-x)
						- (std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b));
					const Real logBudget = InverseLogBudget(p, 2 * BetaLogRoundoff(a, b, x));
					budget += 2 * std::min(p, 1 - p) * std::expm1(logBudget) * std::exp(-logDensity);
				}
				CAPTURE(x, reflected, budget);
				CHECK(IsWithinAbs(x, reflected, budget));
			}
		}
}

TEST_CASE("Beta inverse preserves explicit forward arithmetic and convergence failures", "[special][beta-inverse-contract]") {
	CHECK_THROWS_AS(RegularizedBetaIInv(Real(0.5), std::numeric_limits<Real>::max(),
		std::numeric_limits<Real>::max()), NumericalMethodError);
	try {
		SpecialFunctionsDetail::BetaFraction(Real(0.5), Real(1e12), Real(1e12));
		FAIL("The legacy fraction must still report exhausted work");
	}
	catch (const ConvergenceError& error) {
		CHECK(error.iterations() == Precision::IncompleteBetaMaxIterations);
		CHECK(error.residual() > Precision::Eps);
	}
}

TEST_CASE("Beta inverse reports its own exhausted refinement budget", "[special][beta-inverse-contract]") {
	try {
		SpecialFunctionsDetail::BetaIInverse(Real(0.25), 2, 3, 1);
		FAIL("The inverse must not return an unchecked result after its own budget");
	}
	catch (const ConvergenceError& error) {
		CHECK(error.iterations() == 1);
		CHECK(error.residual() > Precision::RoundoffZeroThreshold);
	}
}

TEST_CASE("Gamma inverse retains representable extreme lower quantiles", "[special][gamma-inverse-contract]") {
	const Real tolerance = std::max(Precision::EqualityTolerance, Precision::RoundoffZeroThreshold);
	CheckTail(RegularizedGammaPInv(2, TOL(1e-100, 1e-20)),
		TOL(Real(1.414213562373095048801688724209698079e-50L), Real(1.4142135623730950488e-10L)), tolerance);
	CheckTail(RegularizedGammaPInv(2, TOL(1e-300, 1e-36)),
		TOL(Real(1.414213562373095048801688724209698079e-150L), Real(1.4142135623730950488e-18L)), tolerance);
	CheckTail(RegularizedGammaPInv(3, TOL(1e-300, 1e-36)),
		TOL(Real(1.817120592832139658891211756327260502e-100L), Real(1.8171205928321396589e-12L)), tolerance);
}

TEST_CASE("Inverse tails retain finite quantiles at the minimum positive probabilities", "[special][beta-inverse-contract][gamma-inverse-contract]") {
	const Real tolerance = std::max(Precision::EqualityTolerance, Precision::RoundoffZeroThreshold);
	for (Real p : {std::numeric_limits<Real>::denorm_min(), std::numeric_limits<Real>::min()}) {
		CAPTURE(p);
		REQUIRE(p > 0);
		const Real root = std::sqrt(p);
		CheckTail(RegularizedBetaIInv(p, 2, 1), root, tolerance);
		CheckTail(RegularizedBetaIInv(p, 2, 2), root / std::sqrt(Real(3)), tolerance);
		CheckTail(RegularizedGammaPInv(2, p), std::sqrt(Real(2)) * root, tolerance);
	}
}

TEST_CASE("Asymptotic inverse paths respect correction and iteration boundaries", "[special][beta-inverse-contract][gamma-inverse-contract]") {
	for (Real factor : {Real(0.5), Real(2)}) {
		const Real x = factor * 3 * Precision::Eps;
		// Independent polynomial/series identities for I_x(2,2) and P(2,x).
		const Real betaP = x * x * (3 - 2 * x);
		const Real gammaP = x * x / 2 * (1 - 2 * x / 3 + x * x / 4);
		CAPTURE(factor, x);
		if (factor < 1) {
			CheckTail(SpecialFunctionsDetail::BetaIInverse(betaP, 2, 2, 2), x);
			CheckTail(SpecialFunctionsDetail::GammaPInverse(2, gammaP, 1), x);
		}
		else {
			CHECK_THROWS_AS(SpecialFunctionsDetail::BetaIInverse(betaP, 2, 2, 2), ConvergenceError);
			CHECK_THROWS_AS(SpecialFunctionsDetail::GammaPInverse(2, gammaP, 1), ConvergenceError);
		}
	}
}

TEST_CASE("Unit shape gamma inverse follows the exact exponential quantile", "[special][gamma-inverse-contract]") {
	for (Real p : {std::numeric_limits<Real>::denorm_min(), std::numeric_limits<Real>::min(),
		Real(1e-100), Real(1e-20), Real(0.01), Real(0.5), Real(0.99), std::nextafter(Real(1), Real(0))}) {
		CAPTURE(p);
		CHECK(RegularizedGammaPInv(1, p) == -std::log1p(-p));
	}
}

TEST_CASE("Gamma inverse rounds genuine lower range loss without an artificial floor", "[special][gamma-inverse-contract]") {
	for (Real a : {Real(0.1), Real(0.25), Real(0.5)}) {
		const Real result = RegularizedGammaPInv(a, std::numeric_limits<Real>::denorm_min());
		CHECK(result == 0);
		CHECK_FALSE(std::signbit(result));
	}
	const Real unit = std::numeric_limits<Real>::denorm_min();
	const Real rootUnit = std::sqrt(unit);
	CHECK(RegularizedGammaPInv(Real(0.5), rootUnit / 2) == 0);
	CHECK(RegularizedGammaPInv(Real(0.5), rootUnit) == unit);
	CHECK(RegularizedGammaPInv(Real(0.5), 2 * rootUnit) == 3 * unit);
	CHECK(RegularizedGammaPInv(Real(0.5), 4 * rootUnit) == 13 * unit);
}

TEST_CASE("Gamma inverse resolves upper quantiles against independent erfc", "[special][gamma-inverse-contract]") {
	for (Real p : {Real(0.1), Real(0.5), Real(0.99), std::nextafter(Real(1), Real(0))}) {
		const Real x = RegularizedGammaPInv(Real(0.5), p);
		CAPTURE(p, x);
		CheckTail(std::erfc(std::sqrt(x)), 1 - p);
	}
}

TEST_CASE("Gamma inverse validates strict domains before exact endpoints", "[special][gamma-inverse-contract]") {
	const Real infinity = std::numeric_limits<Real>::infinity();
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	for (Real a : {Real(0), Real(-1), infinity, -infinity, nan})
		for (Real p : {Real(0), Real(0.5), Real(1)})
			CHECK_THROWS_AS(RegularizedGammaPInv(a, p), DomainError);
	for (Real p : {-std::numeric_limits<Real>::epsilon(), 1 + std::numeric_limits<Real>::epsilon(), infinity, -infinity, nan})
		for (Real a : {Real(0.5), Real(1), Real(2)})
			CHECK_THROWS_AS(RegularizedGammaPInv(a, p), DomainError);
	for (Real a : {Real(0.1), Real(0.5), Real(1), Real(2), Real(1000)}) {
		for (Real p : {Real(0), -Real(0)}) {
			const Real x = RegularizedGammaPInv(a, p);
			CHECK(x == 0);
			CHECK_FALSE(std::signbit(x));
		}
		CHECK(RegularizedGammaPInv(a, 1) == infinity);
	}
}

TEST_CASE("Gamma inverse remains monotone and preserves ordinary forward values", "[special][gamma-inverse-contract]") {
	for (Real a : {Real(0.1), Real(0.25), Real(0.5), Real(1), Real(2), Real(10), Real(100), Real(1000)}) {
		Real previous = 0;
		for (Real p : {Real(1e-100), Real(1e-20), Real(0.01), Real(0.1), Real(0.5), Real(0.9), Real(0.99), std::nextafter(Real(1), Real(0))}) {
			CAPTURE(a, p);
			Real x = 0;
			REQUIRE_NOTHROW(x = RegularizedGammaPInv(a, p));
			CAPTURE(x);
			CHECK(std::isfinite(x));
			CHECK(x >= previous);
			if (p >= Real(0.01) && p <= Real(0.99)) {
				REQUIRE(x > 0);
				const Real logBudget = InverseLogBudget(p, 2 * GammaLogRoundoff(a, x));
				const Real density = std::exp((a - 1) * std::log(x) - x - std::lgamma(a));
				const Real ulp = std::nextafter(x, std::numeric_limits<Real>::infinity()) - x;
				const Real budget = std::min(p, 1 - p) * std::expm1(logBudget) + ulp * density;
				CHECK_THAT(RegularizedGammaP(a, x), Catch::Matchers::WithinAbs(p,
					TOL(p * Precision::EqualityTolerance, std::max(p * Precision::EqualityTolerance, budget))));
			}
			previous = x;
		}
	}
}

TEST_CASE("Gamma inverse preserves explicit normalization failures", "[special][gamma-inverse-contract]") {
	CHECK_THROWS_AS(RegularizedGammaPInv(std::numeric_limits<Real>::max(), Real(0.5)), NumericalMethodError);
}

TEST_CASE("Gamma inverse reports its own exhausted refinement budget", "[special][gamma-inverse-contract]") {
	try {
		SpecialFunctionsDetail::GammaPInverse(2, Real(0.5), 1);
		FAIL("The inverse must not return an unchecked result after its own budget");
	}
	catch (const ConvergenceError& error) {
		CHECK(error.iterations() == 1);
		CHECK(error.residual() > Precision::RoundoffZeroThreshold);
	}
}

TEST_CASE("Shared sine cosine integrals expose exact signed zero limits", "[special][sici-endpoints]") {
	const Real negativeInfinity = -std::numeric_limits<Real>::infinity();
	for (Real zero : {Real(0), -Real(0)}) {
		Real si = Real(0.25), ci = Real(-0.75);
		SinCosIntegral(zero, si, ci);
		CHECK(si == 0);
		CHECK(std::signbit(si) == std::signbit(zero));
		CHECK(ci == negativeInfinity);
		const Real sine = SineIntegral(zero);
		CHECK(sine == 0);
		CHECK(std::signbit(sine) == std::signbit(zero));
		CHECK(CosineIntegral(zero) == negativeInfinity);
	}
}

TEST_CASE("Shared sine cosine integrals expose exact infinite argument limits", "[special][sici-endpoints]") {
	const Real infinity = std::numeric_limits<Real>::infinity();
	for (Real x : {infinity, -infinity}) {
		const Real expectedSine = std::copysign(Constants::PI / 2, x);
		Real si, ci;
		SinCosIntegral(x, si, ci);
		CHECK(si == expectedSine);
		CHECK(ci == 0);
		CHECK_FALSE(std::signbit(ci));
		CHECK(SineIntegral(x) == expectedSine);
		const Real cosine = CosineIntegral(x);
		CHECK(cosine == 0);
		CHECK_FALSE(std::signbit(cosine));
	}
}

TEST_CASE("Shared sine cosine integrals reject NaN before modifying outputs", "[special][sici-endpoints]") {
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	for (Real x : {nan, -nan}) {
		Real si = Real(0.25), ci = Real(-0.75);
		CHECK_THROWS_AS(SinCosIntegral(x, si, ci), DomainError);
		CHECK(si == Real(0.25));
		CHECK(ci == Real(-0.75));
		CHECK_THROWS_AS(SineIntegral(x), DomainError);
		CHECK_THROWS_AS(CosineIntegral(x), DomainError);
	}
}

TEST_CASE("Nonzero sine cosine inputs retain their finite logarithmic limit", "[special][sici-compatibility]") {
	const Real tinyMethodBoundary = std::sqrt(4 * std::numeric_limits<Real>::min());
	const Real eulerGamma = Real(0.577215664901532860606512090082402431L);
	for (Real magnitude : {std::numeric_limits<Real>::denorm_min(), std::numeric_limits<Real>::min(),
		std::sqrt(std::numeric_limits<Real>::min()), std::nextafter(tinyMethodBoundary, Real(0)),
		tinyMethodBoundary, std::nextafter(tinyMethodBoundary, std::numeric_limits<Real>::infinity()),
		std::numeric_limits<Real>::epsilon() / 2}) {
		REQUIRE(magnitude > 0);
		for (Real x : {magnitude, -magnitude}) {
			CAPTURE(x);
			Real si, ci;
			SinCosIntegral(x, si, ci);
			CHECK(si == x);
			CHECK(std::signbit(si) == std::signbit(x));
			CHECK(std::isfinite(ci));
			CHECK(IsNearlyEqual(ci, std::log(magnitude) + eulerGamma, Precision::RoundoffZeroThreshold));
			CHECK(SineIntegral(x) == si);
			CHECK(CosineIntegral(x) == ci);
		}
	}
}

TEST_CASE("Finite sine cosine parity and wrapper consistency are preserved", "[special][sici-compatibility]") {
	for (Real x : {Real(0.25), Real(1), std::nextafter(Real(2), Real(0)), Real(2),
		std::nextafter(Real(2), std::numeric_limits<Real>::infinity()), Real(10), Real(1000)}) {
		CAPTURE(x);
		Real si, ci, negativeSi, negativeCi;
		SinCosIntegral(x, si, ci);
		SinCosIntegral(-x, negativeSi, negativeCi);
		CHECK(std::isfinite(si));
		CHECK(std::isfinite(ci));
		CHECK(negativeSi == -si);
		CHECK(negativeCi == ci);
		CHECK(SineIntegral(x) == si);
		CHECK(CosineIntegral(x) == ci);
	}
	CHECK(IsNearlyEqual(SineIntegral(1), Real(0.94608307036718301494135331382318L), Precision::RoundoffZeroThreshold));
	CHECK(IsNearlyEqual(CosineIntegral(1), Real(0.33740392290096813466264620388915L), Precision::RoundoffZeroThreshold));
	CHECK(IsNearlyEqual(SineIntegral(10), Real(1.65834759421887404933097187938967L), Precision::RoundoffZeroThreshold));
	CHECK(IsNearlyEqual(CosineIntegral(10), Real(-0.04545643300445537263453282995263L), Precision::RoundoffZeroThreshold));
}

TEST_CASE("Shared gamma upper tail matches independent exponential references", "[special][probability-tails]") {
	for (Real x : {Real(1), Real(2), Real(40), Real(100), Real(500), Real(700)}) {
		Real term = std::exp(-x);
		if (term < std::numeric_limits<Real>::min())
			continue;
		Real sum = term;
		for (int a = 1; a <= 6; ++a) {
			if (a > 1) {
				term *= x / (a - 1);
				sum += term;
			}
			CAPTURE(a, x);
			CheckTail(RegularizedGammaQ(Real(a), x), sum, GammaTailTolerance(Real(a), x));
		}
	}
}

TEST_CASE("Shared gamma half-shape upper tail matches the standard library erfc", "[special][probability-tails]") {
	for (Real x : {Real(0.1), Real(1), Real(2), Real(10), Real(40), Real(100), Real(500)}) {
		Real expected = std::erfc(std::sqrt(x));
		if (expected >= std::numeric_limits<Real>::min())
			CheckTail(RegularizedGammaQ(Real(0.5), x), expected, GammaTailTolerance(Real(0.5), x));
	}
}

TEST_CASE("Shared gamma lower tail and complementary branch retain ordinary values", "[special][probability-tails]") {
	for (Real x : {Real(1e-20), Real(0.1), Real(1), Real(2), Real(40)}) {
		CheckTail(RegularizedGammaP(1, x), -std::expm1(-x), GammaTailTolerance(1, x));
		CheckTail(RegularizedGammaQ(1, x), std::exp(-x), GammaTailTolerance(1, x));
		CHECK_THAT(RegularizedGammaP(1, x) + RegularizedGammaQ(1, x), WithinRel(Real(1), Precision::EqualityTolerance));
	}
}

TEST_CASE("Shared gamma kernels preserve endpoints and reject invalid domains", "[special][probability-tails]") {
	CHECK(RegularizedGammaP(2, 0) == 0);
	CHECK(RegularizedGammaQ(2, 0) == 1);
	CHECK(RegularizedGammaP(2, std::numeric_limits<Real>::infinity()) == 1);
	CHECK(RegularizedGammaQ(2, std::numeric_limits<Real>::infinity()) == 0);
	for (Real a : {Real(0), Real(-1), std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		CHECK_THROWS_AS(RegularizedGammaP(a, 1), DomainError);
		CHECK_THROWS_AS(RegularizedGammaQ(a, 1), DomainError);
	}
	for (Real x : {Real(-1), std::numeric_limits<Real>::quiet_NaN()}) {
		CHECK_THROWS_AS(RegularizedGammaP(1, x), DomainError);
		CHECK_THROWS_AS(RegularizedGammaQ(1, x), DomainError);
	}
}

TEST_CASE("Legacy gamma series explicitly reports its exhausted work budget", "[special][probability-tails]") {
	try {
		SpecialFunctionsDetail::GammaLowerSeriesLog(Real(1e6), Real(1e6));
		FAIL("The exhausted gamma series must not return an unconverged probability");
	}
	catch (const ConvergenceError& error) {
		CHECK(error.iterations() == Precision::IncompleteGammaMaxIterations);
		CHECK(error.residual() > Precision::Eps);
	}
}

TEST_CASE("Unregularized lower gamma retains finite values beyond complete gamma range", "[special][gamma-scaling]") {
	// Independent mpmath gammainc references at 90 decimal digits.
	const Real shapes[] = {200, 500, 200, 200, 1000};
	const Real limits[] = {1, 2, Real(0.1), 20, 1};
	const Real expected[] = {
		Real(0.001848593963122709083131623858),
		Real(8.895616059109884099857857384e146),
		Real(4.526439044354275457088826613e-203),
		Real(1.838953687713328083545091318e249),
		Real(0.0003682473202451029346381449979)
	};
	for (int i = 0; i < 5; ++i) {
		if (std::isfinite(expected[i]) && expected[i] > 0) {
			CAPTURE(shapes[i], limits[i]);
			CheckTail(IncompleteGammaLower(shapes[i], limits[i]), expected[i]);
		}
	}
}

TEST_CASE("Unregularized upper gamma retains finite tails after regularization underflow", "[special][gamma-scaling]") {
	const Real shapes[] = {200, 200, 170, 200};
	const Real limits[] = {2000, 1000, 1000, 1600};
	const Real expected[] = {
		Real(2.298767256466362778315741555e-212),
		Real(6.335068535476873445332340430e162),
		Real(6.106763595778072306920042577e72),
		Real(6.402408688546389022318648866e-58)
	};
	for (int i = 0; i < 4; ++i) {
		if (std::isfinite(expected[i]) && expected[i] > 0) {
			CAPTURE(shapes[i], limits[i]);
			CheckTail(IncompleteGammaUpper(shapes[i], limits[i]), expected[i]);
		}
	}
}

TEST_CASE("Unregularized gamma preserves complements and recurrences on both branches", "[special][gamma-scaling]") {
	for (Real a : {Real(0.5), Real(1), Real(2), Real(6)}) {
		for (Real x : {Real(0.1), Real(1), a + 1, Real(20)}) {
			CAPTURE(a, x);
			const Real lower = IncompleteGammaLower(a, x);
			const Real upper = IncompleteGammaUpper(a, x);
			const Real tolerance = GammaTailTolerance(a + 1, x);
			CheckTail(lower + upper, std::tgamma(a), tolerance);
			CheckTail(lower, RegularizedGammaP(a, x) * std::tgamma(a), tolerance);
			CheckTail(upper, RegularizedGammaQ(a, x) * std::tgamma(a), tolerance);
			const Real boundary = std::exp(a * std::log(x) - x);
			CheckTail(IncompleteGammaLower(a + 1, x) + boundary, a * lower, tolerance);
			CheckTail(IncompleteGammaUpper(a + 1, x), a * upper + boundary, tolerance);
		}
	}
}

TEST_CASE("Unregularized gamma endpoints avoid zero times overflowing complete gamma", "[special][gamma-scaling]") {
	const Real infinity = std::numeric_limits<Real>::infinity();
	for (Real a : {Real(0.5), Real(2), Real(200)}) {
		CAPTURE(a);
		CHECK(IncompleteGammaLower(a, 0) == 0);
		CHECK(IncompleteGammaUpper(a, 0) == std::tgamma(a));
		CHECK(IncompleteGammaLower(a, infinity) == std::tgamma(a));
		CHECK(IncompleteGammaUpper(a, infinity) == 0);
	}
}

TEST_CASE("Unregularized gamma reports genuine IEEE range loss without NaNs", "[special][gamma-scaling]") {
	const Real infinity = std::numeric_limits<Real>::infinity();
	if (std::isinf(std::tgamma(Real(200)))) {
		CHECK(IncompleteGammaLower(200, 200) == infinity);
		CHECK(IncompleteGammaUpper(200, 200) == infinity);
	}
	const Real farTail = Real(16) * std::numeric_limits<Real>::max_exponent;
	if (std::isinf(std::tgamma(Real(200))))
		CHECK(IncompleteGammaLower(200, farTail) == infinity);
	else
		CheckTail(IncompleteGammaLower(200, farTail), std::tgamma(Real(200)));
	CHECK(IncompleteGammaUpper(200, farTail) == 0);
	if constexpr (std::is_same_v<Real, double>) {
		const Real subnormal = std::exp(Real(-740));
		REQUIRE(subnormal > 0);
		CHECK(IncompleteGammaUpper(1, 740) == subnormal);
	}
	const Real maximumShape = std::numeric_limits<Real>::max();
	CheckTail(maximumShape * IncompleteGammaLower(maximumShape, 1), std::exp(Real(-1)));
	CHECK(IncompleteGammaLower(std::numeric_limits<Real>::denorm_min(), 1) == infinity);
}

TEST_CASE("Unregularized gamma rejects invalid domains before endpoint handling", "[special][gamma-scaling]") {
	for (auto evaluate : {&IncompleteGammaLower, &IncompleteGammaUpper}) {
		for (Real a : {Real(0), Real(-1), std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
			CHECK_THROWS_AS(evaluate(a, 0), DomainError);
			CHECK_THROWS_AS(evaluate(a, std::numeric_limits<Real>::infinity()), DomainError);
		}
		for (Real x : {Real(-1), std::numeric_limits<Real>::quiet_NaN()})
			CHECK_THROWS_AS(evaluate(2, x), DomainError);
	}
}

TEST_CASE("Shape quadrature reports exhausted work instead of returning an iterate", "[special][gamma-scaling][forward-shape-contract]") {
	try {
		SpecialFunctionsDetail::PositiveShapeIntegral([](SpecialFunctionsDetail::ShapeReal x) { return std::exp(-x); },
			SpecialFunctionsDetail::ShapeReal(1), 8);
		FAIL("Quadrature must not return an unchecked value");
	}
	catch (const ConvergenceError& error) {
		CHECK(error.iterations() == 8);
		CHECK(std::isinf(error.residual()));
	}
	for (auto evaluate : {&IncompleteGammaLower, &IncompleteGammaUpper}) {
		CHECK(evaluate(Real(1e6), Real(1e6)) == std::numeric_limits<Real>::infinity());
	}
}

TEST_CASE("Unregularized gamma retains complementary tails at vanishing shapes", "[special][gamma-scaling]") {
	for (Real a : {std::numeric_limits<Real>::min(), std::numeric_limits<Real>::denorm_min()})
		CheckTail(IncompleteGammaUpper(a, Real(0.5)), Real(0.5597735947761608117467959393150852352L));
}

TEST_CASE("Lambert real branches distinguish the exact endpoint from its neighbors", "[special][lambertw-contracts]") {
	const Real endpoint = Real(-0.367879441171442321595523770161460867L);
	const Real below = std::nextafter(endpoint, -std::numeric_limits<Real>::infinity());
	const Real above = std::nextafter(endpoint, Real(0));
	CHECK(LambertW0(endpoint) == -1);
	CHECK(LambertWm1(endpoint) == -1);
	CHECK_THROWS_AS(LambertW0(below), DomainError);
	CHECK_THROWS_AS(LambertWm1(below), DomainError);
	CHECK(LambertW0(above) > -1);
	CHECK(LambertWm1(above) < -1);
}

TEST_CASE("Lambert branch-point values match independent high-precision references", "[special][lambertw-contracts]") {
	if constexpr (std::numeric_limits<Real>::digits == 53) {
		// mpmath at 100 decimal digits, using the exact binary inputs rather than rounded decimal offsets.
		const Real inputs[] = {
			Real(-0x1.78b56362cef37p-2), Real(-0x1.78b56361f30c8p-2), Real(-0x1.78b56362ca8dap-2),
			Real(-0x1.78b562b7027c7p-2), Real(-0x1.78b2c44c1dd71p-2), Real(-0x1.789b2c7fe3d74p-2)
		};
		const Real principal[] = {
			Real(-0.9999999846957458715), Real(-0.9999835128792687693), Real(-0.9999976683981105763),
			Real(-0.9997668537219889770), Real(-0.9926447551971230044), Real(-0.9768628655744246052)
		};
		const Real secondary[] = {
			Real(-1.0000000153042542846), Real(-1.0000164873019499892), Real(-1.0000023316055136743),
			Real(-1.0002331825217691424), Real(-1.0073914890313083653), Real(-1.0234996190820795719)
		};
		for (int i = 0; i < 6; ++i) {
			CAPTURE(inputs[i]);
			CHECK_THAT(LambertW0(inputs[i]), Catch::Matchers::WithinAbs(principal[i], Precision::RoundoffZeroThreshold));
			CHECK_THAT(LambertWm1(inputs[i]), Catch::Matchers::WithinAbs(secondary[i], Precision::RoundoffZeroThreshold));
		}
	}
}

TEST_CASE("Lambert branches are strictly monotone immediately above the endpoint", "[special][lambertw-contracts]") {
	Real x = Real(-0.367879441171442321595523770161460867L);
	Real principal = -1;
	Real secondary = -1;
	for (int i = 0; i < 32; ++i) {
		x = std::nextafter(x, Real(0));
		const Real nextPrincipal = LambertW0(x);
		const Real nextSecondary = LambertWm1(x);
		CAPTURE(x);
		CHECK(std::isfinite(nextPrincipal));
		CHECK(std::isfinite(nextSecondary));
		CHECK(nextPrincipal > principal);
		CHECK(nextSecondary < secondary);
		principal = nextPrincipal;
		secondary = nextSecondary;
	}
}

TEST_CASE("Lambert branches preserve ordinary values and logarithmic defining identities", "[special][lambertw-contracts]") {
	for (Real x : {Real(-0.35), Real(-0.3), Real(-0.2), Real(-1e-30), Real(1e-30),
	               Real(0.1), Real(1.49), Real(1.5), Real(1.51), Real(2), Real(100),
	               std::numeric_limits<Real>::max()}) {
		const Real w = LambertW0(x);
		CAPTURE(x);
		REQUIRE(std::isfinite(w));
		CHECK(w > -1);
		CHECK((w > 0) == (x > 0));
		CHECK(IsNearlyEqual(w + std::log(std::abs(w)), std::log(std::abs(x)), Precision::EqualityTolerance));
	}
	CheckTail(LambertW0(Real(1.5)), Real(0.7258613577662262570));
	for (Real x : {Real(-0.35), Real(-0.3), Real(-0.1), Real(-0.001),
	               -std::numeric_limits<Real>::min(), -std::numeric_limits<Real>::denorm_min()}) {
		const Real w = LambertWm1(x);
		CAPTURE(x);
		REQUIRE(std::isfinite(w));
		CHECK(w < -1);
		CHECK(IsNearlyEqual(w + std::log(-w), std::log(-x), Precision::EqualityTolerance));
	}
}

TEST_CASE("Lambert real branches enforce their non-finite and zero domain contracts", "[special][lambertw-contracts]") {
	const Real infinity = std::numeric_limits<Real>::infinity();
	const Real nan = std::numeric_limits<Real>::quiet_NaN();
	CHECK(LambertW0(0) == 0);
	CHECK(LambertW0(-Real(0)) == 0);
	CHECK(LambertW0(infinity) == infinity);
	CHECK_THROWS_AS(LambertW0(nan), DomainError);
	CHECK_THROWS_AS(LambertWm1(nan), DomainError);
	CHECK_THROWS_AS(LambertW0(-infinity), DomainError);
	for (Real x : {Real(0), -Real(0), Real(1), infinity, -infinity})
		CHECK_THROWS_AS(LambertWm1(x), DomainError);
}

TEST_CASE("Ei retains finite positive values past the exponential overflow boundary", "[special][ei-range]") {
	// Independent mpmath Ei references at 100 decimal digits.
	const Real inputs[] = {36, 40, 100, 700, 709, 710, 715};
	const Real expected[] = {
		Real(1.232852079912097685430891908e14), Real(6.039718263611241578359231419e15),
		Real(2.715552744853879821914014642e41), Real(1.450978736052560852620882522e301),
		Real(1.160794336657263602480355632e305), Real(3.150915688206201214938641192e305),
		Real(4.643625670370571100383392157e307)
	};
	for (int i = 0; i < 7; ++i) {
		if (std::isfinite(expected[i])) {
			CAPTURE(inputs[i]);
			CheckTail(ExpIntegralEi(inputs[i]), expected[i]);
		}
	}
}

TEST_CASE("Ei distinguishes the last finite double input from genuine overflow", "[special][ei-range]") {
	if constexpr (std::numeric_limits<Real>::digits == 53 && std::numeric_limits<Real>::max_exponent == 1024) {
		const Real lastFinite = Real(0x1.662d80b6ceef0p+9);
		CheckTail(ExpIntegralEi(lastFinite), Real(1.797693134862150202646889473e308));
		CHECK(ExpIntegralEi(std::nextafter(lastFinite, std::numeric_limits<Real>::infinity())) ==
			std::numeric_limits<Real>::infinity());
	}
	CHECK(ExpIntegralEi(std::numeric_limits<Real>::max()) == std::numeric_limits<Real>::infinity());
}

TEST_CASE("Ei preserves its method boundary and negative finite tails", "[special][ei-range]") {
	const Real boundary = -std::log(Precision::Eps);
	const Real below = ExpIntegralEi(std::nextafter(boundary, Real(0)));
	const Real above = ExpIntegralEi(std::nextafter(boundary, std::numeric_limits<Real>::infinity()));
	CheckTail(below, above);
	for (Real x : {Real(1), Real(10), Real(100)}) {
		CAPTURE(x);
		CHECK(ExpIntegralEi(-x) == -ExpIntegralE1(x));
	}
	const Real tail = ExpIntegralEi(-std::log(std::numeric_limits<Real>::max()));
	CHECK(std::isfinite(tail));
	CHECK(tail < 0);
}

TEST_CASE("Ei enforces explicit NaN zero and infinite endpoint contracts", "[special][ei-range]") {
	const Real infinity = std::numeric_limits<Real>::infinity();
	CHECK_THROWS_AS(ExpIntegralEi(0), DomainError);
	CHECK_THROWS_AS(ExpIntegralEi(-Real(0)), DomainError);
	CHECK_THROWS_AS(ExpIntegralEi(std::numeric_limits<Real>::quiet_NaN()), DomainError);
	CHECK(ExpIntegralEi(infinity) == infinity);
	CHECK(ExpIntegralEi(-infinity) == 0);
	CHECK(std::signbit(ExpIntegralEi(-infinity)));
}

TEST_CASE("Dawson ordinary sampling remains compatible with independent references", "[special][dawson-safe]") {
	const Real inputs[] = {Real(0.1), Real(0.2), 1, 2, Real(9.999)};
	const Real expected[] = {
		Real(0.09933599239785286659), Real(0.19475103336802805989), Real(0.53807950691276841914),
		Real(0.30134038892379196603), Real(0.05025892464699347300)
	};
	for (int i = 0; i < 5; ++i) {
		CAPTURE(inputs[i]);
		CHECK_THAT(Dawson(inputs[i]), WithinRel(expected[i], Precision::LooseTolerance));
		CHECK(Dawson(-inputs[i]) == -Dawson(inputs[i]));
	}
}

TEST_CASE("Dawson asymptotic tail matches independent representable references", "[special][dawson-safe]") {
	const Real inputs[] = {10, 20, 100, 10000};
	const Real expected[] = {
		Real(0.0502538471875985280327), Real(0.0250313679264036719470),
		Real(0.0050002500375093782827), Real(0.00005000000025000000375)
	};
	for (int i = 0; i < 4; ++i) {
		CAPTURE(inputs[i]);
		CheckTail(Dawson(inputs[i]), expected[i]);
		CHECK(Dawson(-inputs[i]) == -Dawson(inputs[i]));
	}
}

TEST_CASE("Dawson method switches preserve the ordinary approximation accuracy", "[special][dawson-safe]") {
	for (Real boundary : {Precision::DawsonSmallLimit, Precision::DawsonAsymptoticLimit}) {
		CAPTURE(boundary);
		const Real below = Dawson(std::nextafter(boundary, Real(0)));
		const Real above = Dawson(std::nextafter(boundary, std::numeric_limits<Real>::infinity()));
		CHECK(std::isfinite(below));
		CHECK(std::isfinite(above));
		CHECK_THAT(below, WithinRel(above, Precision::LooseTolerance));
	}
}

TEST_CASE("Dawson large finite arguments never require an integer grid index", "[special][dawson-extreme]") {
	for (Real x : {Real(std::numeric_limits<int>::max()) / 2, Real(std::numeric_limits<int>::max()),
	               Real(std::numeric_limits<int>::max()) * 4,
	               std::sqrt(std::numeric_limits<Real>::max()), std::numeric_limits<Real>::max()}) {
		CAPTURE(x);
		CheckTail(Dawson(x), Real(0.5) / x);
		CHECK(Dawson(-x) == -Dawson(x));
	}
	const Real maximum = std::numeric_limits<Real>::max();
	CHECK(Dawson(maximum) > 0);
	CHECK(Dawson(maximum) < std::numeric_limits<Real>::min());
}

TEST_CASE("Dawson enforces explicit NaN signed-zero and infinite endpoint contracts", "[special][dawson-extreme]") {
	const Real infinity = std::numeric_limits<Real>::infinity();
	CHECK_THROWS_AS(Dawson(std::numeric_limits<Real>::quiet_NaN()), DomainError);
	CHECK(Dawson(infinity) == 0);
	CHECK_FALSE(std::signbit(Dawson(infinity)));
	CHECK(Dawson(-infinity) == 0);
	CHECK(std::signbit(Dawson(-infinity)));
	CHECK(Dawson(0) == 0);
	CHECK(Dawson(-Real(0)) == 0);
	CHECK(std::signbit(Dawson(-Real(0))));
	const Real tiny = std::numeric_limits<Real>::denorm_min();
	CHECK(Dawson(tiny) == tiny);
	CHECK(Dawson(-tiny) == -tiny);
}

TEST_CASE("Complementary inverse error function retains tiny representable inputs", "[special][probability-tails]") {
	// Independently generated with mpmath at 400 decimal digits.
	const Real inputs[] = {Real(1e-20), Real(1e-100), Real(1e-300)};
	const Real expected[] = {Real(6.6015806223551425615), Real(15.065574702592645704), Real(26.209469960516123886)};
	for (int i = 0; i < 3; ++i) {
		if (inputs[i] > 0) {
			CAPTURE(inputs[i]);
			CheckTail(ErfcInv(inputs[i]), expected[i]);
			CheckTail(std::erfc(ErfcInv(inputs[i])), inputs[i]);
		}
	}
}

TEST_CASE("Complementary inverse error function retains ordinary symmetry and endpoints", "[special][probability-tails]") {
	for (Real x : {Real(0.01), Real(0.3), Real(1), Real(2), Real(5), Real(10)}) {
		Real y = std::erfc(x);
		if (y >= std::numeric_limits<Real>::min()) {
			// Rounding erfc(x) to Real loses absolute input accuracy near y=1.
			const Real inputUlp = std::nextafter(y, std::numeric_limits<Real>::infinity()) - y;
			const Real inverseSensitivity = std::sqrt(Constants::PI) / 2 * std::exp(x * x);
			CHECK_THAT(ErfcInv(y), Catch::Matchers::WithinAbs(x,
				Precision::EqualityTolerance * x + inputUlp * inverseSensitivity));
		}
		if (2 - y < 2)
			CHECK_THAT(ErfcInv(2 - y), WithinRel(-ErfcInv(2 - (2 - y)), Precision::EqualityTolerance));
	}
	CHECK(ErfcInv(1) == 0);
	CHECK(ErfcInv(0) == std::numeric_limits<Real>::infinity());
	CHECK(ErfcInv(2) == -std::numeric_limits<Real>::infinity());
	CHECK_THROWS_AS(ErfcInv(-1), DomainError);
	CHECK_THROWS_AS(ErfcInv(3), DomainError);
	CHECK_THROWS_AS(ErfcInv(std::numeric_limits<Real>::quiet_NaN()), DomainError);
}

TEST_CASE("Complementary inverse error function resolves the smallest double input", "[special][probability-tails]") {
	if constexpr (std::numeric_limits<Real>::digits == 53 && std::numeric_limits<Real>::max_exponent == 1024)
		CheckTail(ErfcInv(std::numeric_limits<Real>::denorm_min()), Real(27.213293210812948815));
}
