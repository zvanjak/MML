///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML) Tests                            ///
///                                                                                   ///
///  File:        quasi_random_tests.cpp                                              ///
///  Description: Tests for QuasiRandom.h (Halton, Sobol) and quasi-Monte Carlo        ///
///               integration                                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

#include "../TestPrecision.h"
#include "../../mml/base/QuasiRandom.h"
#include "../../mml/core/Integration/MonteCarloIntegration.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <functional>
#include <vector>
#include <cmath>

using namespace MML;
using Catch::Approx;

namespace MML::Tests::Base::QuasiRandomTests {

template<int N>
struct ScalarFn : public IScalarFunction<N> {
	std::function<Real(const VectorN<Real, N>&)> f;
	explicit ScalarFn(std::function<Real(const VectorN<Real, N>&)> fn) : f(std::move(fn)) {}
	Real operator()(const VectorN<Real, N>& x) const override { return f(x); }
};

TEST_CASE("QuasiRandom::Halton_radical_inverse", "[quasirandom][halton]") {
	TEST_PRECISION_INFO();
	REQUIRE(HaltonSequence::RadicalInverse(1, 2) == Approx(0.5).epsilon(TOL(1e-12, 1e-6)));
	REQUIRE(HaltonSequence::RadicalInverse(2, 2) == Approx(0.25).epsilon(TOL(1e-12, 1e-6)));
	REQUIRE(HaltonSequence::RadicalInverse(3, 2) == Approx(0.75).epsilon(TOL(1e-12, 1e-6)));
	REQUIRE(HaltonSequence::RadicalInverse(1, 3) == Approx(1.0 / 3.0).epsilon(TOL(1e-12, 1e-6)));
	REQUIRE(HaltonSequence::RadicalInverse(4, 2) == Approx(0.125).epsilon(TOL(1e-12, 1e-6)));
}

TEST_CASE("QuasiRandom::Halton_range_and_determinism", "[quasirandom][halton]") {
	TEST_PRECISION_INFO();
	HaltonSequence a(3), b(3);
	for (int i = 0; i < 500; ++i) {
		VectorN<Real, 3> pa = a.Next<3>();
		VectorN<Real, 3> pb = b.Next<3>();
		for (int d = 0; d < 3; ++d) {
			REQUIRE(pa[d] >= 0.0);
			REQUIRE(pa[d] < 1.0);
			REQUIRE(pa[d] == pb[d]);  // deterministic
		}
	}
	REQUIRE_THROWS_AS(HaltonSequence(0), ArgumentError);
	REQUIRE_THROWS_AS(HaltonSequence(41), ArgumentError);
}

TEST_CASE("QuasiRandom::Sobol_range_and_determinism", "[quasirandom][sobol]") {
	TEST_PRECISION_INFO();
	SobolSequence a(4), b(4);
	for (int i = 0; i < 1000; ++i) {
		VectorN<Real, 4> pa = a.Next<4>();
		VectorN<Real, 4> pb = b.Next<4>();
		for (int d = 0; d < 4; ++d) {
			REQUIRE(pa[d] >= 0.0);
			REQUIRE(pa[d] < 1.0);
			REQUIRE(pa[d] == pb[d]);  // deterministic
		}
	}
	SobolSequence c(2);
	c.Next<2>();
	c.Reset();
	VectorN<Real, 2> first = c.Next<2>();
	SobolSequence d(2);
	VectorN<Real, 2> firstAgain = d.Next<2>();
	REQUIRE(first[0] == firstAgain[0]);
	REQUIRE(first[1] == firstAgain[1]);
	REQUIRE_THROWS_AS(SobolSequence(0), ArgumentError);
	REQUIRE_THROWS_AS(SobolSequence(7), ArgumentError);
}

TEST_CASE("QuasiRandom::Sobol_uniformity_1D_mean", "[quasirandom][sobol]") {
	TEST_PRECISION_INFO();
	SobolSequence s(1);
	double sum = 0.0;
	const int n = 1 << 16;
	for (int i = 0; i < n; ++i) { Real u; s.Next(&u); sum += static_cast<double>(u); }
	REQUIRE(sum / n == Approx(0.5).margin(1e-3));  // mean of low-discrepancy points ~ 0.5
}

TEST_CASE("QMC::Sobol_integration_2D", "[quasirandom][qmc]") {
	TEST_PRECISION_INFO();
	// ∫∫_[0,1]^2 x*y dx dy = 1/4
	ScalarFn<2> f([](const VectorN<Real, 2>& x) { return x[0] * x[1]; });
	MonteCarloIntegrator<2> integ;
	VectorN<Real, 2> lo{ 0.0, 0.0 }, hi{ 1.0, 1.0 };
	auto r = integ.integrateQuasiRandom(f, lo, hi, 1 << 16, /*useSobol=*/true);
	REQUIRE(r.value == Approx(0.25).margin(5e-4));
	REQUIRE(r.algorithm_name == std::string("QuasiMonteCarloIntegrator"));
}

TEST_CASE("QMC::Sobol_integration_3D_sum", "[quasirandom][qmc]") {
	TEST_PRECISION_INFO();
	// ∫_[0,1]^3 (x+y+z) = 3/2
	ScalarFn<3> f([](const VectorN<Real, 3>& x) { return x[0] + x[1] + x[2]; });
	MonteCarloIntegrator<3> integ;
	VectorN<Real, 3> lo{ 0.0, 0.0, 0.0 }, hi{ 1.0, 1.0, 1.0 };
	auto r = integ.integrateQuasiRandom(f, lo, hi, 1 << 16, true);
	REQUIRE(r.value == Approx(1.5).margin(1e-3));
}

TEST_CASE("QMC::Halton_fallback_high_dimension", "[quasirandom][qmc][halton]") {
	TEST_PRECISION_INFO();
	// N = 8 > Sobol max dim (6): integrateQuasiRandom falls back to Halton.
	// ∫_[0,1]^8 (sum x_i) = 8 * 1/2 = 4.
	ScalarFn<8> f([](const VectorN<Real, 8>& x) {
		Real s = 0.0; for (int d = 0; d < 8; ++d) s += x[d]; return s; });
	MonteCarloIntegrator<8> integ;
	VectorN<Real, 8> lo, hi;
	for (int d = 0; d < 8; ++d) { lo[d] = 0.0; hi[d] = 1.0; }
	auto r = integ.integrateQuasiRandom(f, lo, hi, 1 << 16, true);  // useSobol ignored for N>6
	REQUIRE(r.value == Approx(4.0).margin(2e-2));
}

TEST_CASE("QMC::beats_plain_MC_on_smooth_integrand", "[quasirandom][qmc]") {
	TEST_PRECISION_INFO();
	// Smooth integrand over [0,1]^2 with known integral; QMC (Sobol) should be at least
	// as accurate as plain MC at the same sample count.
	// ∫∫ exp(x+y) = (e-1)^2 ≈ 2.952492442...
	const double truth = (std::exp(1.0) - 1.0) * (std::exp(1.0) - 1.0);
	ScalarFn<2> f([](const VectorN<Real, 2>& x) { return std::exp(x[0] + x[1]); });
	MonteCarloIntegrator<2> integ;
	VectorN<Real, 2> lo{ 0.0, 0.0 }, hi{ 1.0, 1.0 };

	const size_t n = 1 << 15;
	auto q = integ.integrateQuasiRandom(f, lo, hi, n, true);
	MonteCarloConfig cfg; cfg.num_samples = n; cfg.seed = 12345;
	auto m = integ.integrate(f, lo, hi, cfg);

	double qErr = std::abs(static_cast<double>(q.value) - truth);
	double mErr = std::abs(static_cast<double>(m.value) - truth);
	REQUIRE(qErr <= mErr + 1e-12);   // QMC no worse than plain MC here
	REQUIRE(qErr < 1e-3);            // and QMC is accurate in absolute terms
}

} // namespace MML::Tests::Base::QuasiRandomTests
