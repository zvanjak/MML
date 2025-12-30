///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        RealFFT.h                                                           ///
///  Description: Optimized FFT for real-valued signals                               ///
///               Exploits Hermitian symmetry for 2x performance                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_REAL_FFT_H
#define MML_REAL_FFT_H

#include "MMLBase.h"

#include "FFT.h"

namespace MML 
{
	///////////////////////////////////////////////////////////////////////////////////////////
	// RealFFT - Optimized FFT for real-valued data using conjugate symmetry
	//
	// PERFORMANCE:
	//   - ~2x faster than complex FFT for real data
	//   - Uses packed storage format (n real -> n/2+1 complex)
	//
	// THEORY:
	//   Real signals have conjugate-symmetric spectra: X[k] = X*[N-k]
	//   Only need to compute/store half the spectrum + DC and Nyquist terms
	//
	// REFERENCE: Numerical Recipes fourier.h realft()
	///////////////////////////////////////////////////////////////////////////////////////////

	class RealFFT {
	public:
		///////////////////////////////////////////////////////////////////////////////////////////
		// Forward - Real to complex transform
		//
		// INPUT:  n real samples
		// OUTPUT: n/2+1 complex frequency components
		//         spectrum[0] = DC component (k=0)
		//         spectrum[k] = positive frequencies (k=1 to n/2-1)
		//         spectrum[n/2] = Nyquist frequency (k=n/2)
		//
		// CONJUGATE SYMMETRY:
		//   For k > n/2: X[k] = conj(X[n-k])
		//   This is implicit - we don't store negative frequencies
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Complex> Forward(const Vector<Real>& data) {
			int n = data.size();
			if (n < 2 || !FFT::IsPowerOfTwo(n))
				throw std::invalid_argument("RealFFT::Forward requires power-of-2 size >= 2");

			// Pack real data into complex vector (treating pairs as complex numbers)
			int n2 = n / 2;
			Vector<Complex> temp(n2);
			for (int i = 0; i < n2; i++)
				temp[i] = Complex(data[2 * i], data[2 * i + 1]);

			// Complex FFT of packed data
			FFT::Transform(temp, 1);

			// Unpack using conjugate symmetry to get real FFT result
			Vector<Complex> result(n2 + 1);

			// DC component (sum of even + i*sum of odd)
			result[0] = Complex(temp[0].real() + temp[0].imag(), 0);

			// Nyquist component (difference of even - i*odd)
			result[n2] = Complex(temp[0].real() - temp[0].imag(), 0);

			// Middle components using symmetry relations
			Real theta = -Constants::PI / n2; // -pi/n for positive frequencies
			for (int k = 1; k < n2; k++) {
				Complex w = std::exp(Complex(0, k * theta)); // exp(-i*2*pi*k/n)

				// Get symmetric pairs
				Complex z1 = temp[k];
				Complex z2 = std::conj(temp[n2 - k]);

				// Even part (real DFT of even samples)
				Complex even = Real(0.5) * (z1 + z2);

				// Odd part (real DFT of odd samples)
				Complex odd = Complex(0, -0.5) * (z1 - z2);

				// Combine with twiddle factor
				result[k] = even + w * odd;
			}

			return result;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Inverse - Complex to real transform
		//
		// INPUT:  n/2+1 complex frequency components (conjugate-symmetric spectrum)
		// OUTPUT: n real samples
		//
		// NOTE: Input must satisfy conjugate symmetry for real output
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> Inverse(const Vector<Complex>& spectrum) {
			int n2p1 = spectrum.size();
			int n2 = n2p1 - 1;
			int n = 2 * n2;

			if (n2 < 1)
				throw std::invalid_argument("RealFFT::Inverse requires spectrum size >= 2");

			// Pack spectrum back into complex vector using reverse process
			Vector<Complex> temp(n2);

			// Reconstruct DC and Nyquist into temp[0]
			temp[0] =
				Complex(Real(0.5) * (spectrum[0].real() + spectrum[n2].real()), Real(0.5) * (spectrum[0].real() - spectrum[n2].real()));

			// Reconstruct middle components
			Real theta = Constants::PI / n2; // pi/n for inverse
			for (int k = 1; k < n2; k++) {
				Complex w = std::exp(Complex(0, k * theta)); // exp(i*2*pi*k/n)

				// Get spectrum values and conjugate pair
				Complex z1 = spectrum[k];
				Complex z2 = std::conj(spectrum[n2 - k]);

				// Even and odd parts
				Complex even = Real(0.5) * (z1 + z2);
				Complex odd = Complex(0, 0.5) * (z1 - z2);

				// Recombine
				temp[k] = even + w * odd;
			}

			// Inverse complex FFT
			FFT::Transform(temp, -1);

			// Unpack complex to real (interleaved real/imag become consecutive reals)
			Vector<Real> result(n);
			for (int i = 0; i < n2; i++) {
				result[2 * i] = temp[i].real();
				result[2 * i + 1] = temp[i].imag();
			}

			// Apply normalization
			for (int i = 0; i < n; i++)
				result[i] /= n2;

			return result;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// In-place transform (memory-efficient packed format)
		//
		// INPUT/OUTPUT: data vector
		//   isign =  1: Real->Complex (input: n real, output: packed format)
		//   isign = -1: Complex->Real (input: packed format, output: n real)
		//
		// PACKED FORMAT (NR convention):
		//   data[0] = DC component
		//   data[1] = Nyquist component
		//   data[2*k] = real part of X[k] for k=1..n/2-1
		//   data[2*k+1] = imag part of X[k] for k=1..n/2-1
		///////////////////////////////////////////////////////////////////////////////////////////
		static void Transform(Vector<Real>& data, int isign = 1) {
			int n = data.size();
			if (n < 2 || !FFT::IsPowerOfTwo(n))
				throw std::invalid_argument("RealFFT::Transform requires power-of-2 size >= 2");

			int n2 = n / 2;

			// Pack into complex and do complex FFT
			Vector<Complex> temp(n2);
			for (int i = 0; i < n2; i++)
				temp[i] = Complex(data[2 * i], data[2 * i + 1]);

			FFT::Transform(temp, isign);

			if (isign == 1) // Real to complex
			{
				// Unpack into NR format
				data[0] = temp[0].real() + temp[0].imag(); // DC
				data[1] = temp[0].real() - temp[0].imag(); // Nyquist

				Real theta = -Constants::PI / n2;
				for (int k = 1; k < n2; k++) {
					Complex w = std::exp(Complex(0, k * theta));
					Complex z1 = temp[k];
					Complex z2 = std::conj(temp[n2 - k]);

					Complex even = Real(0.5) * (z1 + z2);
					Complex odd = Complex(0, -0.5) * (z1 - z2);
					Complex result = even + w * odd;

					data[2 * k] = result.real();
					data[2 * k + 1] = result.imag();
				}
			} else // Complex to real (inverse)
			{
				// Reconstruct temp from packed format
				temp[0] = Complex(Real(0.5) * (data[0] + data[1]), Real(0.5) * (data[0] - data[1]));

				Real theta = Constants::PI / n2;
				for (int k = 1; k < n2; k++) {
					Complex spectrum_k(data[2 * k], data[2 * k + 1]);
					Complex spectrum_nk(data[n - 2 * k], -data[n - 2 * k + 1]); // Conjugate

					Complex w = std::exp(Complex(0, k * theta));
					Complex even = Real(0.5) * (spectrum_k + spectrum_nk);
					Complex odd = Complex(0, 0.5) * (spectrum_k - spectrum_nk);

					temp[k] = even + w * odd;
				}

				FFT::Transform(temp, -1);

				// Unpack and normalize
				for (int i = 0; i < n2; i++) {
					data[2 * i] = temp[i].real() / n2;
					data[2 * i + 1] = temp[i].imag() / n2;
				}
			}
		}
	};

} // namespace MML

#endif // MML_REAL_FFT_H
