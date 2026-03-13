#include <catch2/catch_all.hpp>

#include "../TestPrecision.h"
#include "../TestMatchers.h"

#include "MMLBase.h"

#include "algorithms/Fourier.h"

#include "../../test_beds/fourier_test_bed.h"

using namespace MML;
using namespace MML::Fourier;
using namespace MML::TestBeds;

namespace MML::Tests::Algorithms::FourierTests {

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
				REQUIRE(std::abs(recovered[i] - signal[i]) < 1e-12);
			}
		}
	}

	TEST_CASE("DCT-II Known Signal Test", "[Fourier][DCT]") {
		// Test with DC signal (constant)
		Vector<Real> dc(8, 1.0);
		auto dctDC = DCT::ForwardII(dc);

		// DC signal should have all energy in first coefficient
		REQUIRE(std::abs(dctDC[0] - 2.0) < 1e-10); // DCT-II of constant 1 is 2
		for (int i = 1; i < 8; i++) {
			REQUIRE(std::abs(dctDC[i]) < 1e-10);
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
		REQUIRE(relativeError < 1e-10);
	}

	TEST_CASE("DCT Helper - VerifyRoundTrip", "[Fourier][DCT]") {
		Vector<Real> signal(16);
		for (int i = 0; i < 16; i++) {
			signal[i] = std::exp(-0.1 * i) * std::cos(2.0 * Constants::PI * i / 16);
		}

		REQUIRE(DCT::VerifyRoundTrip(signal, 1e-12));
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
			REQUIRE(std::abs(fast_fwd[i] - ref_fwd[i]) < 1e-10);
		}

		// Inverse comparison
		auto inv_fast = DCT::InverseII_Fast(ref_fwd);
		auto inv_ref = DCT::InverseII_Reference(ref_fwd);
		REQUIRE(inv_fast.size() == inv_ref.size());
		for (int i = 0; i < N; i++) {
			REQUIRE(std::abs(inv_fast[i] - inv_ref[i]) < 1e-10);
		}
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
			REQUIRE(std::abs(recovered[i] - signal[i]) < 1e-12);
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
				REQUIRE(std::abs(dct_fast[i] - dct_ref[i]) < 1e-10);
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
		REQUIRE(std::abs(dctSingle[0] - 10.0) < 1e-10); // 2/1 * 5 = 10

		// Two elements
		Vector<Real> two(2);
		two[0] = 1.0;
		two[1] = 2.0;
		auto dctTwo = DCT::ForwardII(two);
		REQUIRE(dctTwo.size() == 2);
		REQUIRE(DCT::VerifyRoundTrip(two));
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
			REQUIRE(std::abs(dftResult[i] - fftResult[i]) < 1e-10);
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
			REQUIRE(std::abs(dftResult[i] - fftResult[i]) < 1e-9);
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

				REQUIRE(maxError < 1e-8);
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
			REQUIRE(std::abs(original[i] - recovered[i]) < 1e-10);
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

				REQUIRE(maxError < 1e-9);
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

		REQUIRE(std::abs(timeEnergy - freqEnergy) / timeEnergy < 1e-10);
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
				REQUIRE(relativeError < 1e-9);
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
		REQUIRE(std::abs(result[0] - Complex(1.0, 0.0)) < 1e-10);
	}

} // namespace MML::Tests::Algorithms::FourierTests
