///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FFT.h                                                               ///
///  Description: Fast Fourier Transform (Cooley-Tukey radix-2 algorithm)             ///
///               O(n log n) forward and inverse transforms                           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FFT_H
#define MML_FFT_H

#include "MMLBase.h"

#include "base/Vector.h"

namespace MML {
	///////////////////////////////////////////////////////////////////////////////////////////
	// FFT - Fast Fourier Transform
	//
	// Cooley-Tukey radix-2 FFT algorithm - O(n log n) complexity
	//
	// Forward transform:  X[k] = Σ_{n=0}^{N-1} x[n] * exp(-2πi*k*n/N)
	// Inverse transform:  x[n] = (1/N) * Σ_{k=0}^{N-1} X[k] * exp(2πi*k*n/N)
	//
	// Requirements:
	// - Input size must be a power of 2 (use NextPowerOfTwo to pad if needed)
	// - For arbitrary sizes, use DFT or zero-pad to next power of 2
	//
	// Algorithm: Decimation-in-time (DIT) with bit-reversal permutation
	//
	// Performance: ~100x faster than DFT for n=1024, ~1000x for n=4096
	///////////////////////////////////////////////////////////////////////////////////////////
	class FFT {
	public:
		///////////////////////////////////////////////////////////////////////////////////////////
		// Core FFT algorithms
		///////////////////////////////////////////////////////////////////////////////////////////

		// In-place FFT (modifies input array)
		// isign: +1 for forward, -1 for inverse (caller must normalize by 1/N)
		static void Transform(Vector<Complex>& data, int isign = 1) {
			int n = data.size();

			if (n == 0) {
				throw std::invalid_argument("FFT::Transform - empty input vector");
			}

			if (!IsPowerOfTwo(n)) {
				throw std::invalid_argument("FFT::Transform - size must be power of 2");
			}

			// Bit-reversal permutation
			BitReversalPermute(data);

			// Cooley-Tukey FFT with iterative decimation-in-time
			int mmax = 1;
			while (n > mmax) {
				int istep = mmax << 1;						// istep = 2 * mmax
				Real theta = -isign * Constants::PI / mmax; // Negative sign for DIT

				// Complex exponential for twiddle factor
				Complex wp(std::cos(theta), std::sin(theta));
				Complex w(1.0, 0.0);

				for (int m = 0; m < mmax; m++) {
					for (int i = m; i < n; i += istep) {
						int j = i + mmax;

						// Butterfly operation
						Complex temp = w * data[j];
						Complex u = data[i];
						data[j] = u - temp;
						data[i] = u + temp;
					}
					w *= wp; // Update twiddle factor
				}
				mmax = istep;
			}
		}

		// Forward FFT: time domain -> frequency domain
		static Vector<Complex> Forward(const Vector<Complex>& data) {
			Vector<Complex> result = data; // Copy input
			Transform(result, 1);
			return result;
		}

		// Inverse FFT: frequency domain -> time domain
		static Vector<Complex> Inverse(const Vector<Complex>& spectrum) {
			Vector<Complex> result = spectrum; // Copy input
			Transform(result, -1);			   // Inverse transform

			// Normalize by 1/N
			Real norm = 1.0 / result.size();
			for (int i = 0; i < result.size(); i++) {
				result[i] *= norm;
			}

			return result;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Utility functions
		///////////////////////////////////////////////////////////////////////////////////////////

		// Check if n is a power of 2
		static bool IsPowerOfTwo(int n) { return (n > 0) && ((n & (n - 1)) == 0); }

		// Find next power of 2 >= n
		static int NextPowerOfTwo(int n) {
			if (n <= 0)
				return 1;
			if (IsPowerOfTwo(n))
				return n;

			int power = 1;
			while (power < n) {
				power <<= 1;
			}
			return power;
		}

		// Zero-pad vector to next power of 2
		static Vector<Complex> ZeroPad(const Vector<Complex>& data) {
			int n = data.size();
			int padded_size = NextPowerOfTwo(n);

			if (padded_size == n) {
				return data; // Already power of 2
			}

			Vector<Complex> padded(padded_size, Complex(0.0, 0.0));
			for (int i = 0; i < n; i++) {
				padded[i] = data[i];
			}

			return padded;
		}

		// Compute log2 of n (assumes n is power of 2)
		static int Log2(int n) {
			int log = 0;
			while (n > 1) {
				n >>= 1;
				log++;
			}
			return log;
		}

	private:
		///////////////////////////////////////////////////////////////////////////////////////////
		// Internal helper functions
		///////////////////////////////////////////////////////////////////////////////////////////

		// Bit-reversal permutation (in-place)
		static void BitReversalPermute(Vector<Complex>& data) {
			int n = data.size();
			int j = 0;

			for (int i = 0; i < n - 1; i++) {
				if (i < j) {
					std::swap(data[i], data[j]);
				}

				// Bit-reversal increment
				int k = n >> 1;
				while (k <= j) {
					j -= k;
					k >>= 1;
				}
				j += k;
			}
		}
	};

} // namespace MML

#endif // MML_FFT_H
