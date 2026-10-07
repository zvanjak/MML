#include <catch2/catch_all.hpp>

#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#include <mml/MMLBase.h>

#include <mml/algorithms/Fourier/FourierConvolution.h>
#include <mml/algorithms/Fourier/Fourier.h>
#include <mml/algorithms/Fourier/FourierRealFFT.h>
#include <mml/algorithms/Fourier/FourierSpectrum.h>
#include <mml/algorithms/Fourier/FourierWindowing.h>

#include "../../test_beds/fourier_test_bed.h"

using namespace MML;
using namespace MML::Fourier;
using namespace MML::TestBeds;

namespace MML::Tests::Algorithms::FourierTests {

	///////////////////////////////////////////////////////////////////////////////////////////
	//                          WINDOW FUNCTIONS AND METRICS                                   //
	///////////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Fourier windows - additional cosine families", "[Fourier][Windowing]") {
		SECTION("Blackman-Harris has expected endpoints and symmetry") {
			auto window = Windows::BlackmanHarris(17);

			REQUIRE(window.size() == 17);
			REQUIRE(std::abs(window[0] - REAL(0.00006)) < TOL(1e-12, 1e-6));
			REQUIRE(std::abs(window[8] - REAL(1.0)) < TOL(1e-12, 1e-6));
			for (int i = 0; i < 8; i++)
				REQUIRE(std::abs(window[i] - window[16 - i]) < TOL(1e-12, 1e-6));
		}

		SECTION("FlatTop has expected endpoints and symmetry") {
			auto window = Windows::FlatTop(17);

			REQUIRE(window.size() == 17);
			REQUIRE(std::abs(window[0]) < TOL(1e-12, 1e-6));
			REQUIRE(std::abs(window[8] - REAL(4.636)) < TOL(1e-12, 1e-5));
			for (int i = 0; i < 8; i++)
				REQUIRE(std::abs(window[i] - window[16 - i]) < TOL(1e-12, 1e-5));
		}
	}

	TEST_CASE("Fourier windows - Tukey behavior", "[Fourier][Windowing]") {
		auto rectangular = Windows::Rectangular(8);
		auto hann = Windows::Hann(8);
		auto tukey_rect = Windows::Tukey(8, 0.0);
		auto tukey_hann = Windows::Tukey(8, 1.0);

		for (int i = 0; i < 8; i++) {
			REQUIRE(std::abs(tukey_rect[i] - rectangular[i]) < TOL(1e-12, 1e-6));
			REQUIRE(std::abs(tukey_hann[i] - hann[i]) < TOL(1e-12, 1e-6));
		}

		auto tapered = Windows::Tukey(9, 0.5);
		REQUIRE(std::abs(tapered[0]) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(tapered[4] - REAL(1.0)) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(tapered[8]) < TOL(1e-12, 1e-6));
		REQUIRE_THROWS_AS(Windows::Tukey(8, -0.1), std::invalid_argument);
		REQUIRE_THROWS_AS(Windows::Tukey(8, 1.1), std::invalid_argument);
	}

	TEST_CASE("Fourier windows - apply helper and metrics", "[Fourier][Windowing]") {
		Vector<Real> signal(4);
		signal[0] = 1.0;
		signal[1] = 2.0;
		signal[2] = 3.0;
		signal[3] = 4.0;

		auto hann = Windows::Hann(4);
		auto windowed = Windows::ApplyWindow(signal, hann);
		REQUIRE(std::abs(windowed[0]) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(windowed[1] - REAL(1.5)) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(windowed[2] - REAL(2.25)) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(windowed[3]) < TOL(1e-12, 1e-6));

		Vector<Complex> complex_signal(4);
		for (int i = 0; i < 4; i++)
			complex_signal[i] = Complex(signal[i], -signal[i]);
		auto complex_windowed = Windows::ApplyWindow(complex_signal, hann);
		REQUIRE(std::abs(complex_windowed[1] - Complex(REAL(1.5), REAL(-1.5))) < TOL(1e-12, 1e-6));

		auto metrics = Windows::Metrics(hann);
		REQUIRE(std::abs(metrics.sum - REAL(1.5)) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(metrics.sum_squares - REAL(1.125)) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(metrics.coherent_gain - REAL(0.375)) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(metrics.rms_gain - std::sqrt(REAL(0.28125))) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(metrics.enbw_bins - REAL(2.0)) < TOL(1e-12, 1e-6));

		Vector<Real> bad_window(3, 1.0);
		REQUIRE_THROWS_AS(Windows::ApplyWindow(signal, bad_window), std::invalid_argument);
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	//                          SPECTRUM HELPERS                                               //
	///////////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Fourier spectrum helpers - magnitude power and phase", "[Fourier][Spectrum]") {
		Vector<Complex> spectrum(4);
		spectrum[0] = Complex(3.0, 4.0);
		spectrum[1] = Complex(1.0, -1.0);
		spectrum[2] = Complex(-2.0, 0.0);
		spectrum[3] = Complex(0.0, 2.0);

		auto magnitude = Magnitude(spectrum);
		auto magnitude_squared = MagnitudeSquared(spectrum);
		auto power = Power(spectrum);
		auto phase = Phase(spectrum);

		REQUIRE(std::abs(magnitude[0] - REAL(5.0)) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(magnitude[1] - std::sqrt(REAL(2.0))) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(magnitude_squared[0] - REAL(25.0)) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(power[1] - magnitude_squared[1]) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(phase[0] - std::atan2(REAL(4.0), REAL(3.0))) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(phase[2] - Constants::PI) < TOL(1e-12, 1e-6));

		Vector<Complex> empty;
		REQUIRE_THROWS_AS(Magnitude(empty), std::invalid_argument);
		REQUIRE_THROWS_AS(Phase(empty), std::invalid_argument);
	}

	TEST_CASE("Fourier spectrum helpers - phase unwrap", "[Fourier][Spectrum]") {
		Vector<Real> wrapped(5);
		wrapped[0] = 0.0;
		wrapped[1] = 0.75 * Constants::PI;
		wrapped[2] = -0.75 * Constants::PI;
		wrapped[3] = -0.5 * Constants::PI;
		wrapped[4] = 0.5 * Constants::PI;

		auto unwrapped = UnwrapPhase(wrapped);

		REQUIRE(std::abs(unwrapped[0] - REAL(0.0)) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(unwrapped[1] - 0.75 * Constants::PI) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(unwrapped[2] - 1.25 * Constants::PI) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(unwrapped[3] - 1.5 * Constants::PI) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(unwrapped[4] - 2.5 * Constants::PI) < TOL(1e-12, 1e-6));
		REQUIRE_THROWS_AS(UnwrapPhase(wrapped, 0.0), std::invalid_argument);
	}

	TEST_CASE("Fourier spectrum helpers - frequency axis layouts", "[Fourier][Spectrum]") {
		auto one_sided = FrequencyAxis(8, 80.0, FrequencyAxisLayout::OneSided);
		REQUIRE(one_sided.size() == 5);
		for (int i = 0; i < one_sided.size(); i++)
			REQUIRE(std::abs(one_sided[i] - i * REAL(10.0)) < TOL(1e-12, 1e-6));

		auto two_sided = FrequencyAxis(8, 80.0);
		std::array<Real, 8> expected_two_sided = {0.0, 10.0, 20.0, 30.0, -40.0, -30.0, -20.0, -10.0};
		for (int i = 0; i < two_sided.size(); i++)
			REQUIRE(std::abs(two_sided[i] - expected_two_sided[i]) < TOL(1e-12, 1e-6));

		auto shifted = FrequencyAxis(8, 80.0, FrequencyAxisLayout::TwoSided, true);
		std::array<Real, 8> expected_shifted = {-40.0, -30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0};
		for (int i = 0; i < shifted.size(); i++)
			REQUIRE(std::abs(shifted[i] - expected_shifted[i]) < TOL(1e-12, 1e-6));

		REQUIRE_THROWS_AS(FrequencyAxis(0, 80.0), std::invalid_argument);
		REQUIRE_THROWS_AS(FrequencyAxis(8, 0.0), std::invalid_argument);
	}

	TEST_CASE("Fourier spectrum helpers - FFTShift and IFFTShift", "[Fourier][Spectrum]") {
		Vector<Real> even(6);
		for (int i = 0; i < even.size(); i++)
			even[i] = i;

		auto even_shifted = FFTShift(even);
		std::array<Real, 6> expected_even = {3.0, 4.0, 5.0, 0.0, 1.0, 2.0};
		for (int i = 0; i < even_shifted.size(); i++)
			REQUIRE(std::abs(even_shifted[i] - expected_even[i]) < TOL(1e-12, 1e-6));

		auto even_round_trip = IFFTShift(even_shifted);
		for (int i = 0; i < even.size(); i++)
			REQUIRE(std::abs(even_round_trip[i] - even[i]) < TOL(1e-12, 1e-6));

		Vector<Real> odd(5);
		for (int i = 0; i < odd.size(); i++)
			odd[i] = i;

		auto odd_shifted = FFTShift(odd);
		std::array<Real, 5> expected_odd = {3.0, 4.0, 0.0, 1.0, 2.0};
		for (int i = 0; i < odd_shifted.size(); i++)
			REQUIRE(std::abs(odd_shifted[i] - expected_odd[i]) < TOL(1e-12, 1e-6));

		auto odd_round_trip = IFFTShift(odd_shifted);
		for (int i = 0; i < odd.size(); i++)
			REQUIRE(std::abs(odd_round_trip[i] - odd[i]) < TOL(1e-12, 1e-6));
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	//                          CONVOLUTION HELPERS                                            //
	///////////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Fourier convolution - real linear modes", "[Fourier][Convolution]") {
		Vector<Real> x(4);
		x[0] = 1.0;
		x[1] = 2.0;
		x[2] = 3.0;
		x[3] = 4.0;
		Vector<Real> h(3);
		h[0] = 0.5;
		h[1] = 1.0;
		h[2] = 0.5;

		auto full = Convolve(x, h, ConvolutionMode::Full, ConvolutionMethod::Direct);
		std::array<Real, 6> expected_full = {0.5, 2.0, 4.0, 6.0, 5.5, 2.0};
		REQUIRE(full.size() == 6);
		for (int i = 0; i < full.size(); i++)
			REQUIRE(std::abs(full[i] - expected_full[i]) < TOL(1e-12, 1e-6));

		auto same = Convolve(x, h, ConvolutionMode::Same, ConvolutionMethod::Direct);
		std::array<Real, 4> expected_same = {2.0, 4.0, 6.0, 5.5};
		REQUIRE(same.size() == 4);
		for (int i = 0; i < same.size(); i++)
			REQUIRE(std::abs(same[i] - expected_same[i]) < TOL(1e-12, 1e-6));

		auto valid = Convolve(x, h, ConvolutionMode::Valid, ConvolutionMethod::Direct);
		std::array<Real, 2> expected_valid = {4.0, 6.0};
		REQUIRE(valid.size() == 2);
		for (int i = 0; i < valid.size(); i++)
			REQUIRE(std::abs(valid[i] - expected_valid[i]) < TOL(1e-12, 1e-6));
	}

	TEST_CASE("Fourier convolution - circular mode", "[Fourier][Convolution]") {
		Vector<Real> x(4);
		x[0] = 1.0;
		x[1] = 2.0;
		x[2] = 3.0;
		x[3] = 4.0;
		Vector<Real> h(3);
		h[0] = 1.0;
		h[1] = -1.0;
		h[2] = 0.5;

		auto circular = Convolve(x, h, ConvolutionMode::Circular, ConvolutionMethod::Direct);
		std::array<Real, 4> expected = {-1.5, 3.0, 1.5, 2.0};
		REQUIRE(circular.size() == 4);
		for (int i = 0; i < circular.size(); i++)
			REQUIRE(std::abs(circular[i] - expected[i]) < TOL(1e-12, 1e-6));
	}

	TEST_CASE("Fourier convolution - direct and FFT paths agree", "[Fourier][Convolution]") {
		Vector<Real> x(5);
		Vector<Real> h(4);
		for (int i = 0; i < x.size(); i++)
			x[i] = std::sin(REAL(0.3) * i) + i;
		for (int i = 0; i < h.size(); i++)
			h[i] = std::cos(REAL(0.4) * i) - REAL(0.25) * i;

		for (auto mode : {ConvolutionMode::Full, ConvolutionMode::Same, ConvolutionMode::Valid, ConvolutionMode::Circular}) {
			auto direct = Convolve(x, h, mode, ConvolutionMethod::Direct);
			auto fft = Convolve(x, h, mode, ConvolutionMethod::FFT);
			REQUIRE(direct.size() == fft.size());
			for (int i = 0; i < direct.size(); i++)
				REQUIRE(std::abs(direct[i] - fft[i]) < TOL(1e-10, 1e-5));
		}
	}

	TEST_CASE("Fourier convolution - complex values and validation", "[Fourier][Convolution]") {
		Vector<Complex> x(2);
		x[0] = Complex(1.0, 1.0);
		x[1] = Complex(2.0, -1.0);
		Vector<Complex> h(2);
		h[0] = Complex(0.5, 0.0);
		h[1] = Complex(0.0, 1.0);

		auto full = Convolve(x, h, ConvolutionMode::Full, ConvolutionMethod::Direct);
		REQUIRE(full.size() == 3);
		REQUIRE(std::abs(full[0] - Complex(0.5, 0.5)) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(full[1] - Complex(0.0, 0.5)) < TOL(1e-12, 1e-6));
		REQUIRE(std::abs(full[2] - Complex(1.0, 2.0)) < TOL(1e-12, 1e-6));

		Vector<Real> empty;
		Vector<Real> non_empty(1, 1.0);
		REQUIRE_THROWS_AS(Convolve(empty, non_empty), std::invalid_argument);
		REQUIRE_THROWS_AS(Convolve(non_empty, empty), std::invalid_argument);
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	//                          TEST 11: DCT (DISCRETE COSINE TRANSFORM)                      //
	///////////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("DCT-II Forward/Inverse Round-Trip", "[Fourier][DCT]") {
		// Test with various signal sizes
		std::vector<int> sizes = {4, 8, 16, 32, 64};

		for (int N : sizes) {
			Vector<Real> signal(N);
			for (int i = 0; i < N; i++) {
				signal[i] = std::sin(2.0 * Constants::PI * i / N) + 0.5 * std::cos(4.0 * Constants::PI * i / N);
			}

			auto coeffs = DCT::ForwardII(signal);
			auto recovered = DCT::InverseII(coeffs);

			REQUIRE(recovered.size() == signal.size());

			for (int i = 0; i < N; i++) {
				REQUIRE(std::abs(recovered[i] - signal[i]) < TOL(1e-12, 1e-5));
			}
		}
	}

	TEST_CASE("DCT-II Known Signal Test", "[Fourier][DCT]") {
		// Test with DC signal (constant)
		Vector<Real> dc(8, 1.0);
		auto dctDC = DCT::ForwardII(dc);

		// DC signal should have all energy in first coefficient
		REQUIRE(std::abs(dctDC[0] - 2.0) < TOL(1e-10, 1e-5)); // DCT-II of constant 1 is 2
		for (int i = 1; i < 8; i++) {
			REQUIRE(std::abs(dctDC[i]) < TOL(1e-10, 1e-5));
		}
	}

	TEST_CASE("DCT-II Energy Conservation", "[Fourier][DCT]") {
		int N = 32;
		Vector<Real> signal(N);
		for (int i = 0; i < N; i++) {
			signal[i] = std::sin(2.0 * Constants::PI * 5.0 * i / N);
		}

		Real signalEnergy = 0.0;
		for (int i = 0; i < N; i++) {
			signalEnergy += signal[i] * signal[i];
		}

		auto coeffs = DCT::ForwardII(signal);

		// Parseval's theorem for DCT-II with 2/N normalization:
		// ||x||² = (N/2) * [ (1/2)|X[0]|² + Σ_{k=1}^{N-1} |X[k]|² ]
		Real dctEnergy = (N / 2.0) * (0.5 * coeffs[0] * coeffs[0]);
		for (int k = 1; k < N; k++) {
			dctEnergy += (N / 2.0) * coeffs[k] * coeffs[k];
		}

		Real relativeError = std::abs(dctEnergy - signalEnergy) / signalEnergy;
		REQUIRE(relativeError < TOL(1e-10, 1e-5));
	}

	TEST_CASE("DCT Helper - VerifyRoundTrip", "[Fourier][DCT]") {
		Vector<Real> signal(16);
		for (int i = 0; i < 16; i++) {
			signal[i] = std::exp(-0.1 * i) * std::cos(2.0 * Constants::PI * i / 16);
		}

		REQUIRE(DCT::VerifyRoundTrip(signal, TOL(1e-12, 1e-5)));
	}

	TEST_CASE("DCT Fast vs Reference - Non-zero DC", "[Fourier][DCT]") {
		// Verify fast implementation matches reference for non-zero mean signal
		// This exposes the DC coefficient double-halving bug (0.5 * 0.5 = 0.25 instead of 0.5)
		int N = 32;
		Vector<Real> signal(N);
		for (int i = 0; i < N; i++) {
			signal[i] = 3.0 + std::sin(2.0 * Constants::PI * i / N); // DC offset of 3.0
		}

		// Inverse comparison: fast should match reference
		auto ref_fwd = DCT::ForwardII_Reference(signal);
		auto inv_fast = DCT::InverseII_Fast(ref_fwd);
		auto inv_ref = DCT::InverseII_Reference(ref_fwd);
		REQUIRE(inv_fast.size() == inv_ref.size());
		for (int i = 0; i < N; i++) {
			REQUIRE(std::abs(inv_fast[i] - inv_ref[i]) < TOL(1e-10, 1e-5));
		}
	}

	TEST_CASE("DCT Fast vs Reference", "[Fourier][DCT]") {
		// Verify fast implementation matches reference for both forward and inverse
		int N = 32; // Above threshold, uses fast
		Vector<Real> signal(N);
		for (int i = 0; i < N; i++) {
			signal[i] = std::sin(2.0 * Constants::PI * i / N);
		}

		// Forward comparison
		auto fast_fwd = DCT::ForwardII_Fast(signal);
		auto ref_fwd = DCT::ForwardII_Reference(signal);
		REQUIRE(fast_fwd.size() == ref_fwd.size());
		for (int i = 0; i < N; i++) {
			REQUIRE(std::abs(fast_fwd[i] - ref_fwd[i]) < TOL(1e-10, 1e-5));
		}

		// Inverse comparison
		auto inv_fast = DCT::InverseII_Fast(ref_fwd);
		auto inv_ref = DCT::InverseII_Reference(ref_fwd);
		REQUIRE(inv_fast.size() == inv_ref.size());
		for (int i = 0; i < N; i++) {
			REQUIRE(std::abs(inv_fast[i] - inv_ref[i]) < TOL(1e-10, 1e-5));
		}
	}

	TEST_CASE("DCT Fast supports non-power-of-two sizes", "[Fourier][DCT][Bluestein]") {
		const int N = 45;
		Vector<Real> signal(N);
		for (int i = 0; i < N; i++)
			signal[i] = REAL(1.5) + std::sin(REAL(0.21) * i) - REAL(0.3) * std::cos(REAL(0.47) * i);

		const auto fastCoefficients = DCT::ForwardII_Fast(signal);
		const auto referenceCoefficients = DCT::ForwardII_Reference(signal);
		for (int i = 0; i < N; i++)
			REQUIRE(std::abs(fastCoefficients[i] - referenceCoefficients[i]) < TOL(1e-10, 2e-4));

		const auto recovered = DCT::InverseII_Fast(fastCoefficients);
		for (int i = 0; i < N; i++)
			REQUIRE(std::abs(recovered[i] - signal[i]) < TOL(1e-10, 2e-4));
	}

	TEST_CASE("DST-I Forward/Inverse Round-Trip", "[Fourier][DST]") {
		// Test DST with boundary conditions (x[0] = x[N+1] = 0)
		int N = 16;
		Vector<Real> signal(N);
		for (int i = 0; i < N; i++) {
			// Sine function naturally satisfies boundary conditions
			signal[i] = std::sin(Constants::PI * (i + 1) / (N + 1));
		}

		auto coeffs = DCT::ForwardDST(signal);
		auto recovered = DCT::InverseDST(coeffs);

		REQUIRE(recovered.size() == signal.size());

		for (int i = 0; i < N; i++) {
			REQUIRE(std::abs(recovered[i] - signal[i]) < TOL(1e-12, 1e-5));
		}
	}

	TEST_CASE("DST-I Known Signal Test", "[Fourier][DST]") {
		// Pure sine mode: sin(π*k*(n+1)/(N+1)) should give single peak at coefficient k-1
		int N = 8;
		int mode = 3; // Third mode

		Vector<Real> signal(N);
		for (int n = 0; n < N; n++) {
			signal[n] = std::sin(Constants::PI * mode * (n + 1) / (N + 1));
		}

		auto coeffs = DCT::ForwardDST(signal);

		// Should have peak at coefficient (mode-1)
		Real maxCoeff = 0.0;
		int maxIdx = 0;
		for (int k = 0; k < N; k++) {
			if (std::abs(coeffs[k]) > maxCoeff) {
				maxCoeff = std::abs(coeffs[k]);
				maxIdx = k;
			}
		}

		REQUIRE(maxIdx == mode - 1);
		REQUIRE(maxCoeff > 0.9); // Should be close to 1
	}

	TEST_CASE("DCT Performance - Fast vs Reference", "[Fourier][DCT][performance]") {
		// Verify O(N log N) vs O(N^2) behavior
		// This test ensures the fast implementation is significantly faster for large N

		std::vector<int> sizes = {64, 256, 1024, 4096};

		for (int N : sizes) {
			Vector<Real> signal(N);
			for (int i = 0; i < N; i++) {
				signal[i] = std::sin(2.0 * Constants::PI * i / N) + 0.5 * std::cos(4.0 * Constants::PI * i / N);
			}

			// Time fast implementation
			auto start_fast = std::chrono::high_resolution_clock::now();
			auto dct_fast = DCT::ForwardII_Fast(signal);
			auto end_fast = std::chrono::high_resolution_clock::now();
			auto fast_us = std::chrono::duration_cast<std::chrono::microseconds>(end_fast - start_fast).count();

			// Time reference implementation
			auto start_ref = std::chrono::high_resolution_clock::now();
			auto dct_ref = DCT::ForwardII_Reference(signal);
			auto end_ref = std::chrono::high_resolution_clock::now();
			auto ref_us = std::chrono::duration_cast<std::chrono::microseconds>(end_ref - start_ref).count();

			// Results should match
			for (int i = 0; i < N; i++) {
				REQUIRE(std::abs(dct_fast[i] - dct_ref[i]) < TOL(1e-10, 1e-5));
			}

			// For N >= 256, fast should be noticeably faster
			// (For smaller N, overhead may dominate)
			if (N >= 256) {
				// Fast should be at least 2x faster for large N
				// Use INFO to report times but don't fail on timing (platform-dependent)
				INFO("N=" << N << " Fast: " << fast_us << "us, Reference: " << ref_us
									<< "us, Speedup: " << (ref_us > 0 ? (double)ref_us / fast_us : 0));
			}
		}
	}

	TEST_CASE("DCT Edge Cases", "[Fourier][DCT][edge-case]") {
		// Single element
		Vector<Real> single(1);
		single[0] = 5.0;
		auto dctSingle = DCT::ForwardII(single);
		REQUIRE(dctSingle.size() == 1);
		REQUIRE(std::abs(dctSingle[0] - 10.0) < TOL(1e-10, 1e-5)); // 2/1 * 5 = 10

		// Two elements
		Vector<Real> two(2);
		two[0] = 1.0;
		two[1] = 2.0;
		auto dctTwo = DCT::ForwardII(two);
		REQUIRE(dctTwo.size() == 2);
		if constexpr (std::is_same_v<Real, float>) {
			// DCT round-trip accumulates too much error with float precision
			auto recovered = DCT::InverseII(dctTwo);
			for (int i = 0; i < two.size(); ++i)
				REQUIRE_THAT(recovered[i], Catch::Matchers::WithinAbs(two[i], REAL(1e-4)));
		} else {
			REQUIRE(DCT::VerifyRoundTrip(two));
		}
	}


	///////////////////////////////////////////////////////////////////////////////////////////
	//                          TEST 1: DFT vs FFT CORRECTNESS                                //
	///////////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("DFT vs FFT - Pure Tone", "[Fourier][DFT][FFT]") {
		auto signal = getPureToneSignal();
		Vector<Complex> data(signal.data.size());
		for (size_t i = 0; i < signal.data.size(); i++) {
			data[i] = Complex(signal.data[i], 0.0);
		}

		auto dftResult = DFT::Forward(data);
		auto fftResult = FFT::Forward(data);

		REQUIRE(dftResult.size() == fftResult.size());

		for (size_t i = 0; i < dftResult.size(); i++) {
			REQUIRE(std::abs(dftResult[i] - fftResult[i]) < TOL(1e-10, 5e-3));
		}
	}

	TEST_CASE("DFT vs FFT - Dual Tone", "[Fourier][DFT][FFT]") {
		auto signal = getDualToneSignal();
		Vector<Complex> data(signal.data.size());
		for (size_t i = 0; i < signal.data.size(); i++) {
			data[i] = Complex(signal.data[i], 0.0);
		}

		auto dftResult = DFT::Forward(data);
		auto fftResult = FFT::Forward(data);

		REQUIRE(dftResult.size() == fftResult.size());

		for (size_t i = 0; i < dftResult.size(); i++) {
			REQUIRE(std::abs(dftResult[i] - fftResult[i]) < TOL(1e-9, 0.05));
		}
	}

	TEST_CASE("DFT vs FFT - All Test Signals", "[Fourier][DFT][FFT][comprehensive]") {
		auto signals = getAllSignals();

		for (const auto& signal : signals) {
			DYNAMIC_SECTION("Signal: " << signal.name) {
				Vector<Complex> data(signal.data.size());
				for (size_t i = 0; i < signal.data.size(); i++) {
					data[i] = Complex(signal.data[i], 0.0);
				}

				auto dftResult = DFT::Forward(data);
				auto fftResult = FFT::Forward(data);

				REQUIRE(dftResult.size() == fftResult.size());

				Real maxError = 0.0;
				for (size_t i = 0; i < dftResult.size(); i++) {
					Real error = std::abs(dftResult[i] - fftResult[i]);
					maxError = std::max(maxError, error);
				}

				REQUIRE(maxError < TOL(1e-8, 0.05));
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	//                          TEST 2: FFT INVERSE RECOVERY                                  //
	///////////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("FFT Forward/Inverse Round-Trip", "[Fourier][FFT][inverse]") {
		auto signal = getHarmonicSignal();
		Vector<Complex> data(signal.data.size());
		for (size_t i = 0; i < signal.data.size(); i++) {
			data[i] = Complex(signal.data[i], 0.0);
		}

		auto original = data;
		auto spectrum = FFT::Forward(data);
		auto recovered = FFT::Inverse(spectrum);

		REQUIRE(original.size() == recovered.size());

		for (size_t i = 0; i < original.size(); i++) {
			REQUIRE(std::abs(original[i] - recovered[i]) < TOL(1e-10, 1e-5));
		}
	}

	TEST_CASE("FFT Inverse Recovery - All Signals", "[Fourier][FFT][inverse][comprehensive]") {
		auto signals = getAllSignals();

		for (const auto& signal : signals) {
			DYNAMIC_SECTION("Signal: " << signal.name) {
				Vector<Complex> data(signal.data.size());
				for (size_t i = 0; i < signal.data.size(); i++) {
					data[i] = Complex(signal.data[i], 0.0);
				}

				auto original = data;
				auto spectrum = FFT::Forward(data);
				auto recovered = FFT::Inverse(spectrum);

				Real maxError = 0.0;
				for (size_t i = 0; i < original.size(); i++) {
					Real error = std::abs(original[i] - recovered[i]);
					maxError = std::max(maxError, error);
				}

				REQUIRE(maxError < TOL(1e-9, 1e-4));
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	//                          TEST 3: PARSEVAL'S THEOREM                                    //
	///////////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Parseval's Theorem - Energy Conservation", "[Fourier][FFT][Parseval]") {
		auto signal = getPureToneSignal();
		Vector<Complex> data(signal.data.size());
		for (size_t i = 0; i < signal.data.size(); i++) {
			data[i] = Complex(signal.data[i], 0.0);
		}

		// Time domain energy
		Real timeEnergy = 0.0;
		for (size_t i = 0; i < data.size(); i++) {
			timeEnergy += std::norm(data[i]);
		}

		// Frequency domain energy
		auto spectrum = FFT::Forward(data);
		Real freqEnergy = 0.0;
		for (size_t i = 0; i < spectrum.size(); i++) {
			freqEnergy += std::norm(spectrum[i]);
		}
		freqEnergy /= spectrum.size(); // Normalization

		REQUIRE(std::abs(timeEnergy - freqEnergy) / timeEnergy < TOL(1e-10, 1e-5));
	}

	TEST_CASE("Parseval's Theorem - All Signals", "[Fourier][FFT][Parseval][comprehensive]") {
		auto signals = getAllSignals();

		for (const auto& signal : signals) {
			DYNAMIC_SECTION("Signal: " << signal.name) {
				Vector<Complex> data(signal.data.size());
				for (size_t i = 0; i < signal.data.size(); i++) {
					data[i] = Complex(signal.data[i], 0.0);
				}

				Real timeEnergy = 0.0;
				for (size_t i = 0; i < data.size(); i++) {
					timeEnergy += std::norm(data[i]);
				}

				auto spectrum = FFT::Forward(data);
				Real freqEnergy = 0.0;
				for (size_t i = 0; i < spectrum.size(); i++) {
					freqEnergy += std::norm(spectrum[i]);
				}
				freqEnergy /= spectrum.size();

				Real relativeError = std::abs(timeEnergy - freqEnergy) / timeEnergy;
				REQUIRE(relativeError < TOL(1e-9, 1e-4));
			}
		}
	}

	TEST_CASE("FFT - Power of Two Requirement", "[Fourier][FFT][edge-case]") {
		REQUIRE(FFT::IsPowerOfTwo(128));
		REQUIRE(FFT::IsPowerOfTwo(256));
		REQUIRE(!FFT::IsPowerOfTwo(100));
		REQUIRE(!FFT::IsPowerOfTwo(127));

		REQUIRE(FFT::NextPowerOfTwo(100) == 128);
		REQUIRE(FFT::NextPowerOfTwo(128) == 128);
		REQUIRE(FFT::NextPowerOfTwo(129) == 256);
	}

	TEST_CASE("FFT - Empty Input", "[Fourier][FFT][edge-case]") {
		Vector<Complex> empty;
		REQUIRE_THROWS(FFT::Forward(empty));
	}

	TEST_CASE("FFT - Single Element", "[Fourier][FFT][edge-case]") {
		Vector<Complex> single(1);
		single[0] = Complex(1.0, 0.0);
		auto result = FFT::Forward(single);
		REQUIRE(result.size() == 1);
		REQUIRE(std::abs(result[0] - Complex(1.0, 0.0)) < TOL(1e-10, 1e-5));
	}

	TEST_CASE("FFT Bluestein matches DFT for odd and prime sizes", "[Fourier][FFT][Bluestein]") {
		for (int size : {3, 5, 7, 11, 25}) {
			DYNAMIC_SECTION("N=" << size) {
				Vector<Complex> data(size);
				for (int i = 0; i < size; i++)
					data[i] = Complex(std::sin(REAL(0.37) * i) + REAL(0.1) * i,
					                  std::cos(REAL(0.23) * i) - REAL(0.05) * i);

				const auto expectedForward = DFT::Forward(data, TransformNormalization::None);
				const auto actualForward = FFT::Forward(data, TransformNormalization::None);
				for (int i = 0; i < size; i++)
					REQUIRE(std::abs(actualForward[i] - expectedForward[i]) < TOL(1e-10, 2e-4));

				const auto expectedInverse = DFT::Inverse(data, TransformNormalization::None);
				const auto actualInverse = FFT::Inverse(data, TransformNormalization::None);
				for (int i = 0; i < size; i++)
					REQUIRE(std::abs(actualInverse[i] - expectedInverse[i]) < TOL(1e-10, 2e-4));
			}
		}
	}

	TEST_CASE("FFT arbitrary sizes honor every normalization mode", "[Fourier][FFT][Bluestein][Normalization]") {
		for (int size : {6, 13, 31}) {
			Vector<Complex> data(size);
			for (int i = 0; i < size; i++)
				data[i] = Complex(std::sin(REAL(0.19) * i), std::cos(REAL(0.31) * i));

			for (auto normalization : {TransformNormalization::Legacy, TransformNormalization::None,
			                           TransformNormalization::Forward, TransformNormalization::Inverse,
			                           TransformNormalization::Orthonormal}) {
				DYNAMIC_SECTION("N=" << size << " normalization=" << static_cast<int>(normalization)) {
					const auto spectrum = FFT::Forward(data, normalization);
					const auto recovered = FFT::Inverse(spectrum, normalization);
					const Real expectedScale = normalization == TransformNormalization::None ? size : 1.0;
					for (int i = 0; i < size; i++)
						REQUIRE(std::abs(recovered[i] - expectedScale * data[i]) < TOL(1e-9, 5e-4));
				}
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	//                       NORMALIZATION API TESTS                                          //
	///////////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("DFT normalization modes preserve explicit round trips", "[Fourier][Normalization]") {
		Vector<Complex> data(4);
		data[0] = Complex(1.0, 0.5);
		data[1] = Complex(2.0, -1.0);
		data[2] = Complex(-0.5, 0.25);
		data[3] = Complex(0.75, 1.5);

		for (auto normalization : {TransformNormalization::Forward, TransformNormalization::Inverse, TransformNormalization::Orthonormal}) {
			auto spectrum = DFT::Forward(data, normalization);
			auto recovered = DFT::Inverse(spectrum, normalization);
			for (int i = 0; i < data.size(); i++)
				REQUIRE(std::abs(recovered[i] - data[i]) < TOL(1e-12, 1e-6));
		}

		auto raw_spectrum = DFT::Forward(data, TransformNormalization::None);
		auto raw_inverse = DFT::Inverse(raw_spectrum, TransformNormalization::None);
		for (int i = 0; i < data.size(); i++)
			REQUIRE(std::abs(raw_inverse[i] - data[i] * REAL(data.size())) < TOL(1e-12, 1e-6));
	}

	TEST_CASE("FFT normalization modes preserve legacy behavior and Parseval scaling", "[Fourier][Normalization]") {
		int N = 8;
		Vector<Complex> data(N);
		for (int i = 0; i < N; i++)
			data[i] = Complex(std::sin(2.0 * Constants::PI * i / N), std::cos(4.0 * Constants::PI * i / N));

		auto legacy = FFT::Forward(data);
		auto explicit_inverse = FFT::Forward(data, TransformNormalization::Inverse);
		for (int i = 0; i < N; i++)
			REQUIRE(std::abs(legacy[i] - explicit_inverse[i]) < TOL(1e-12, 1e-6));

		for (auto normalization : {TransformNormalization::Forward, TransformNormalization::Inverse, TransformNormalization::Orthonormal}) {
			auto spectrum = FFT::Forward(data, normalization);
			auto recovered = FFT::Inverse(spectrum, normalization);
			for (int i = 0; i < N; i++)
				REQUIRE(std::abs(recovered[i] - data[i]) < TOL(1e-10, 1e-6));
		}

		Real time_energy = 0.0;
		for (int i = 0; i < N; i++)
			time_energy += std::norm(data[i]);

		auto orthonormal_spectrum = FFT::Forward(data, TransformNormalization::Orthonormal);
		Real frequency_energy = 0.0;
		for (int i = 0; i < N; i++)
			frequency_energy += std::norm(orthonormal_spectrum[i]);

		REQUIRE(std::abs(time_energy - frequency_energy) < TOL(1e-10, 1e-5));
	}

	TEST_CASE("RealFFT normalization modes round trip real data", "[Fourier][Normalization][RealFFT]") {
		int N = 16;
		Vector<Real> data(N);
		for (int i = 0; i < N; i++)
			data[i] = std::sin(2.0 * Constants::PI * i / N) + REAL(0.25) * std::cos(6.0 * Constants::PI * i / N);

		for (auto normalization : {TransformNormalization::Forward, TransformNormalization::Inverse, TransformNormalization::Orthonormal}) {
			auto spectrum = RealFFT::Forward(data, normalization);
			auto recovered = RealFFT::Inverse(spectrum, normalization);
			for (int i = 0; i < N; i++)
				REQUIRE(std::abs(recovered[i] - data[i]) < TOL(1e-10, 1e-5));
		}
	}

	TEST_CASE("Fourier detailed APIs report normalization metadata", "[Fourier][Normalization][Detailed]") {
		Vector<Complex> data(8);
		for (int i = 0; i < data.size(); i++)
			data[i] = Complex(std::sin(2.0 * Constants::PI * i / data.size()), 0.0);

		FourierConfig config;
		config.normalization = TransformNormalization::Orthonormal;
		auto result = FFT::ForwardDetailed(data, config);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.input_size == data.size());
		REQUIRE(result.output_size == data.size());
		REQUIRE(result.normalization == TransformNormalization::Orthonormal);
		REQUIRE_FALSE(result.padding_applied);
	}

	///////////////////////////////////////////////////////////////////////////////////////////
	//                       DETAILED API TESTS                                               //
	///////////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("DFT ForwardDetailed - Success", "[Fourier][DFT][Detailed]") {
		Vector<Complex> data(4);
		for (int i = 0; i < 4; i++)
			data[i] = Complex(std::cos(2.0 * Constants::PI * i / 4), 0.0);

		auto result = DFT::ForwardDetailed(data);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.value.size() == 4);
		REQUIRE(result.input_size == 4);
		REQUIRE(result.algorithm_name == "DFT::Forward");
		REQUIRE(result.elapsed_time_ms >= 0.0);
	}

	TEST_CASE("DFT ForwardDetailed - Real Input", "[Fourier][DFT][Detailed]") {
		Vector<Real> data(8);
		for (int i = 0; i < 8; i++)
			data[i] = std::sin(2.0 * Constants::PI * i / 8);

		auto result = DFT::ForwardDetailed(data);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.value.size() == 8);
		REQUIRE(result.input_size == 8);
	}

	TEST_CASE("DFT ForwardDetailed - Empty Input Returns Failure", "[Fourier][DFT][Detailed]") {
		FourierConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		Vector<Complex> empty;

		auto result = DFT::ForwardDetailed(empty, config);

		REQUIRE_FALSE(result.IsSuccess());
		REQUIRE(result.status == AlgorithmStatus::InvalidInput);
		REQUIRE_FALSE(result.error_message.empty());
	}

	TEST_CASE("DFT InverseDetailed - Round Trip", "[Fourier][DFT][Detailed]") {
		Vector<Complex> data(4);
		data[0] = Complex(1.0, 0.0);
		data[1] = Complex(2.0, 0.0);
		data[2] = Complex(3.0, 0.0);
		data[3] = Complex(4.0, 0.0);

		auto fwd = DFT::ForwardDetailed(data);
		REQUIRE(fwd.IsSuccess());

		auto inv = DFT::InverseDetailed(fwd.value);
		REQUIRE(inv.IsSuccess());

		for (int i = 0; i < 4; i++)
			REQUIRE(std::abs(inv.value[i] - data[i]) < TOL(1e-12, 1e-5));
	}

	TEST_CASE("DFT InverseRealDetailed - Round Trip", "[Fourier][DFT][Detailed]") {
		Vector<Real> data(8);
		for (int i = 0; i < 8; i++)
			data[i] = std::sin(2.0 * Constants::PI * i / 8);

		auto fwd = DFT::ForwardDetailed(data);
		REQUIRE(fwd.IsSuccess());

		auto inv = DFT::InverseRealDetailed(fwd.value);
		REQUIRE(inv.IsSuccess());
		REQUIRE(inv.value.size() == 8);

		for (int i = 0; i < 8; i++)
			REQUIRE(std::abs(inv.value[i] - data[i]) < TOL(1e-12, 1e-5));
	}

	TEST_CASE("FFT ForwardDetailed - Success", "[Fourier][FFT][Detailed]") {
		int N = 16;
		Vector<Complex> data(N);
		for (int i = 0; i < N; i++)
			data[i] = Complex(std::sin(2.0 * Constants::PI * i / N), 0.0);

		auto result = FFT::ForwardDetailed(data);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.value.size() == N);
		REQUIRE(result.input_size == N);
		REQUIRE(result.algorithm_name == "FFT::Forward");
	}

	TEST_CASE("FFT Detailed - Round Trip", "[Fourier][FFT][Detailed]") {
		int N = 32;
		Vector<Complex> data(N);
		for (int i = 0; i < N; i++)
			data[i] = Complex(std::cos(4.0 * Constants::PI * i / N), std::sin(2.0 * Constants::PI * i / N));

		auto fwd = FFT::ForwardDetailed(data);
		REQUIRE(fwd.IsSuccess());

		auto inv = FFT::InverseDetailed(fwd.value);
		REQUIRE(inv.IsSuccess());

		for (int i = 0; i < N; i++)
			REQUIRE(std::abs(inv.value[i] - data[i]) < TOL(1e-10, 1e-5));
	}

	TEST_CASE("FFT ForwardDetailed - Empty Input ConvertToStatus", "[Fourier][FFT][Detailed]") {
		FourierConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		Vector<Complex> empty;

		auto result = FFT::ForwardDetailed(empty, config);

		REQUIRE_FALSE(result.IsSuccess());
		REQUIRE(result.status == AlgorithmStatus::InvalidInput);
	}

	TEST_CASE("DCT ForwardIIDetailed - Success", "[Fourier][DCT][Detailed]") {
		Vector<Real> signal(16);
		for (int i = 0; i < 16; i++)
			signal[i] = std::sin(2.0 * Constants::PI * i / 16) + 0.5;

		auto result = DCT::ForwardIIDetailed(signal);

		REQUIRE(result.IsSuccess());
		REQUIRE(result.value.size() == 16);
		REQUIRE(result.input_size == 16);
		REQUIRE(result.algorithm_name == "DCT::ForwardII");
	}

	TEST_CASE("DCT Detailed - Round Trip", "[Fourier][DCT][Detailed]") {
		Vector<Real> signal(8);
		for (int i = 0; i < 8; i++)
			signal[i] = std::cos(2.0 * Constants::PI * i / 8) + 0.3;

		auto fwd = DCT::ForwardIIDetailed(signal);
		REQUIRE(fwd.IsSuccess());

		auto inv = DCT::InverseIIDetailed(fwd.value);
		REQUIRE(inv.IsSuccess());

		for (int i = 0; i < 8; i++)
			REQUIRE(std::abs(inv.value[i] - signal[i]) < TOL(1e-12, 1e-5));
	}

	TEST_CASE("DST Detailed - Round Trip", "[Fourier][DCT][Detailed]") {
		Vector<Real> signal(8);
		for (int i = 0; i < 8; i++)
			signal[i] = std::sin(Constants::PI * (i + 1) / 9);

		auto fwd = DCT::ForwardDSTDetailed(signal);
		REQUIRE(fwd.IsSuccess());

		auto inv = DCT::InverseDSTDetailed(fwd.value);
		REQUIRE(inv.IsSuccess());

		for (int i = 0; i < 8; i++)
			REQUIRE(std::abs(inv.value[i] - signal[i]) < TOL(1e-12, 1e-5));
	}

	TEST_CASE("Detailed API - Propagate Exception Policy", "[Fourier][Detailed]") {
		FourierConfig config;
		config.exception_policy = EvaluationExceptionPolicy::Propagate;
		Vector<Complex> empty;

		REQUIRE_THROWS_AS(DFT::ForwardDetailed(empty, config), std::invalid_argument);
	}

} // namespace MML::Tests::Algorithms::FourierTests
