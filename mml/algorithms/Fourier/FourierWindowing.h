///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Fourier/FourierWindowing.h                                          ///
///  Description: Window functions for spectral analysis and FFT                      ///
///               Hamming, Hann, Blackman, Kaiser windows for leakage reduction       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                                       ///
///                                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_WINDOWING_H
#define MML_WINDOWING_H

#include <mml/MMLBase.h>

#include <mml/base/Vector/Vector.h>

#include <mml/algorithms/Fourier/Fourier.h>

#include <cmath>
#include <stdexcept>

namespace MML::Fourier 
{
	///////////////////////////////////////////////////////////////////////////////////////////
	// Windows - Collection of window functions for spectral analysis
	//
	// PURPOSE:
	//   Window functions reduce spectral leakage in FFT by tapering signal edges.
	//   They smooth the transition to zero at signal boundaries, minimizing artifacts.
	//
	// APPLICATIONS:
	//   - FFT analysis (reduce spectral leakage)
	//   - Filter design (finite impulse response filters)
	//   - Audio processing (smooth transitions)
	//   - Signal segmentation (overlapping windows)
	//
	// TRADE-OFFS:
	//   - Main lobe width: Narrower = better frequency resolution
	//   - Side lobe level: Lower = less spectral leakage
	//   - Each window optimizes different characteristics
	//
	// REFERENCE: Numerical Recipes spectrum.h, Harris (1978) "On the use of windows"
	///////////////////////////////////////////////////////////////////////////////////////////
	namespace Windows {
		///////////////////////////////////////////////////////////////////////////////////////////
		// Rectangular - No windowing (box window)
		//
		// PROPERTIES:
		//   - Narrowest main lobe (best frequency resolution)
		//   - Highest side lobes (worst spectral leakage: -13 dB)
		//   - Use when signal exactly matches FFT period
		//
		// FORMULA: w[n] = 1 for all n
		///////////////////////////////////////////////////////////////////////////////////////////
		inline Vector<Real> Rectangular(int n) {
			FourierValidation::ValidateWindowSize(n);

			return Vector<Real>(n, 1.0);
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Hann (Hanning) - Raised cosine window
		//
		// PROPERTIES:
		//   - Moderate main lobe width
		//   - Good side lobe suppression (-31 dB)
		//   - General-purpose window, widely used
		//
		// FORMULA: w[n] = 0.5 * (1 - cos(2πn/(N-1)))
		///////////////////////////////////////////////////////////////////////////////////////////
		inline Vector<Real> Hann(int n) {
			FourierValidation::ValidateWindowSize(n);

			if (n == 1)
				return Vector<Real>(1, 1.0);

			Vector<Real> window(n);
			for (int i = 0; i < n; i++)
				window[i] = 0.5 * (1.0 - std::cos(2.0 * Constants::PI * i / (n - 1)));

			return window;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Hamming - Optimized raised cosine window
		//
		// PROPERTIES:
		//   - Slightly wider main lobe than Hann
		//   - Better side lobe suppression (-43 dB)
		//   - Optimized for minimizing first side lobe
		//
		// FORMULA: w[n] = 0.54 - 0.46 * cos(2πn/(N-1))
		///////////////////////////////////////////////////////////////////////////////////////////
		inline Vector<Real> Hamming(int n) {
			FourierValidation::ValidateWindowSize(n);

			if (n == 1)
				return Vector<Real>(1, 1.0);

			Vector<Real> window(n);
			for (int i = 0; i < n; i++)
				window[i] = 0.54 - 0.46 * std::cos(2.0 * Constants::PI * i / (n - 1));

			return window;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Blackman - Three-term cosine window
		//
		// PROPERTIES:
		//   - Wide main lobe (worst frequency resolution)
		//   - Excellent side lobe suppression (-58 dB)
		//   - Use when minimizing spectral leakage is critical
		//
		// FORMULA: w[n] = 0.42 - 0.5*cos(2πn/(N-1)) + 0.08*cos(4πn/(N-1))
		///////////////////////////////////////////////////////////////////////////////////////////
		inline Vector<Real> Blackman(int n) {
			FourierValidation::ValidateWindowSize(n);

			if (n == 1)
				return Vector<Real>(1, 1.0);

			Vector<Real> window(n);
			for (int i = 0; i < n; i++) {
				Real angle = 2.0 * Constants::PI * i / (n - 1);
				window[i] = 0.42 - 0.5 * std::cos(angle) + 0.08 * std::cos(2.0 * angle);
			}

			return window;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Blackman-Harris - Four-term cosine window
		//
		// FORMULA: w[n] = a0 - a1*cos(θ) + a2*cos(2θ) - a3*cos(3θ)
		//   where θ = 2πn/(N-1), a0=0.35875, a1=0.48829, a2=0.14128, a3=0.01168
		///////////////////////////////////////////////////////////////////////////////////////////
		inline Vector<Real> BlackmanHarris(int n) {
			FourierValidation::ValidateWindowSize(n);

			if (n == 1)
				return Vector<Real>(1, 1.0);

			Vector<Real> window(n);
			for (int i = 0; i < n; i++) {
				Real angle = 2.0 * Constants::PI * i / (n - 1);
				window[i] = 0.35875 - 0.48829 * std::cos(angle) + 0.14128 * std::cos(2.0 * angle) - 0.01168 * std::cos(3.0 * angle);
			}

			return window;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// FlatTop - Five-term cosine window for amplitude-accurate spectral peaks
		//
		// FORMULA: w[n] = a0 - a1*cos(θ) + a2*cos(2θ) - a3*cos(3θ) + a4*cos(4θ)
		//   where θ = 2πn/(N-1), a0=1, a1=1.93, a2=1.29, a3=0.388, a4=0.028
		///////////////////////////////////////////////////////////////////////////////////////////
		inline Vector<Real> FlatTop(int n) {
			FourierValidation::ValidateWindowSize(n);

			if (n == 1)
				return Vector<Real>(1, 1.0);

			Vector<Real> window(n);
			for (int i = 0; i < n; i++) {
				Real angle = 2.0 * Constants::PI * i / (n - 1);
				window[i] = 1.0 - 1.93 * std::cos(angle) + 1.29 * std::cos(2.0 * angle) - 0.388 * std::cos(3.0 * angle) + 0.028 * std::cos(4.0 * angle);
			}

			return window;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Bartlett (Triangular) - Linear taper window
		//
		// PROPERTIES:
		//   - Moderate main lobe width
		//   - Moderate side lobe suppression (-25 dB)
		//   - Simple triangular shape
		//
		// FORMULA: w[n] = 1 - |2n/(N-1) - 1|
		///////////////////////////////////////////////////////////////////////////////////////////
		inline Vector<Real> Bartlett(int n) {
			FourierValidation::ValidateWindowSize(n);

			if (n == 1)
				return Vector<Real>(1, 1.0);

			Vector<Real> window(n);
			for (int i = 0; i < n; i++)
				window[i] = 1.0 - std::abs(2.0 * i / (n - 1) - 1.0);

			return window;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Welch - Parabolic window
		//
		// PROPERTIES:
		//   - Similar to Bartlett but smoother
		//   - Continuous first derivative (Bartlett has kink at peak)
		//   - Side lobe suppression: -21 dB
		//
		// FORMULA: w[n] = 1 - (2n/(N-1) - 1)²
		///////////////////////////////////////////////////////////////////////////////////////////
		inline Vector<Real> Welch(int n) {
			FourierValidation::ValidateWindowSize(n);

			if (n == 1)
				return Vector<Real>(1, 1.0);

			Vector<Real> window(n);
			for (int i = 0; i < n; i++) {
				Real x = 2.0 * i / (n - 1) - 1.0;
				window[i] = 1.0 - x * x;
			}

			return window;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Kaiser - Parametric window using modified Bessel function
		//
		// PARAMETERS:
		//   beta: Shape parameter controlling window properties
		//     - beta = 0: Rectangular window
		//     - beta = 5: Hamming-like
		//     - beta = 8.6: Blackman-like (-60 dB side lobes)
		//     - Higher beta: Wider main lobe, better side lobe suppression
		//
		// PROPERTIES:
		//   - Adjustable trade-off between main lobe width and side lobe level
		//   - Near-optimal in terms of energy concentration
		//
		// FORMULA: w[n] = I₀(β√(1-(2n/(N-1)-1)²)) / I₀(β)
		//   where I₀ is the zeroth-order modified Bessel function of the first kind
		///////////////////////////////////////////////////////////////////////////////////////////
		inline Vector<Real> Kaiser(int n, Real beta) {
			FourierValidation::ValidateWindowSize(n);
			FourierValidation::ValidateNonNegative(beta, "beta", "Windows::Kaiser");

			// Modified Bessel function I0(x) using series expansion
			auto bessel_i0 = [](Real x) -> Real {
				Real sum = 1.0;
				Real term = 1.0;
				Real x_half = x / 2.0;

				for (int k = 1; k <= 50; k++) // Sufficient for typical accuracy
				{
					term *= (x_half / k);
					term *= (x_half / k);
					sum += term;

					if (term < 1e-12 * sum) // Convergence check
						break;
				}

				return sum;
			};

			Real i0_beta = bessel_i0(beta);
			Vector<Real> window(n);

			for (int i = 0; i < n; i++) {
				Real x = 2.0 * i / (n - 1) - 1.0; // Normalize to [-1, 1]
				Real arg = beta * std::sqrt(1.0 - x * x);
				window[i] = bessel_i0(arg) / i0_beta;
			}

			return window;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Gaussian - Gaussian-shaped window
		//
		// PARAMETERS:
		//   sigma: Standard deviation (larger = narrower window)
		//     - Typical: sigma = 0.4 to 0.5
		//     - sigma = 0.4: Side lobes at -60 dB
		//
		// PROPERTIES:
		//   - Smooth, bell-shaped taper
		//   - Minimal overshoot in time and frequency domains
		//   - Good for signals with Gaussian statistics
		//
		// FORMULA: w[n] = exp(-0.5 * ((n - (N-1)/2) / (σ * (N-1)/2))²)
		///////////////////////////////////////////////////////////////////////////////////////////
		inline Vector<Real> Gaussian(int n, Real sigma) {
			FourierValidation::ValidateWindowSize(n);
			FourierValidation::ValidatePositive(sigma, "sigma", "Windows::Gaussian");

			Vector<Real> window(n);
			Real center = (n - 1) / 2.0;
			Real denominator = sigma * (n - 1) / 2.0;

			for (int i = 0; i < n; i++) {
				Real x = (i - center) / denominator;
				window[i] = std::exp(-0.5 * x * x);
			}

			return window;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Tukey - Tapered cosine window
		//
		// PARAMETERS:
		//   alpha: Taper fraction in [0, 1]. alpha=0 is rectangular, alpha=1 is Hann.
		///////////////////////////////////////////////////////////////////////////////////////////
		inline Vector<Real> Tukey(int n, Real alpha = 0.5) {
			FourierValidation::ValidateWindowSize(n);
			if (alpha < 0.0 || alpha > 1.0)
				throw FourierError("Windows::Tukey: alpha must be in [0, 1]");

			if (alpha == 0.0)
				return Rectangular(n);
			if (alpha == 1.0)
				return Hann(n);
			if (n == 1)
				return Vector<Real>(1, 1.0);

			Vector<Real> window(n, 1.0);
			Real denominator = alpha * (n - 1);
			Real left_edge = alpha * (n - 1) / 2.0;
			Real right_edge = (n - 1) * (1.0 - alpha / 2.0);

			for (int i = 0; i < n; i++) {
				if (i < left_edge) {
					window[i] = 0.5 * (1.0 + std::cos(Constants::PI * (2.0 * i / denominator - 1.0)));
				} else if (i > right_edge) {
					window[i] = 0.5 * (1.0 + std::cos(Constants::PI * (2.0 * i / denominator - 2.0 / alpha + 1.0)));
				}
			}

			return window;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// ApplyWindow - Element-wise window multiplication
		///////////////////////////////////////////////////////////////////////////////////////////
		inline Vector<Real> ApplyWindow(const Vector<Real>& signal, const Vector<Real>& window) {
			FourierValidation::ValidateNonEmpty(signal.size(), "Windows::ApplyWindow (signal)");
			if (signal.size() != window.size())
				throw FourierError("Windows::ApplyWindow: signal and window sizes must match");

			Vector<Real> result(signal.size());
			for (int i = 0; i < signal.size(); i++)
				result[i] = signal[i] * window[i];

			return result;
		}

		inline Vector<Complex> ApplyWindow(const Vector<Complex>& signal, const Vector<Real>& window) {
			FourierValidation::ValidateNonEmpty(signal.size(), "Windows::ApplyWindow (signal)");
			if (signal.size() != window.size())
				throw FourierError("Windows::ApplyWindow: signal and window sizes must match");

			Vector<Complex> result(signal.size());
			for (int i = 0; i < signal.size(); i++)
				result[i] = signal[i] * window[i];

			return result;
		}

		struct WindowMetrics {
			Real sum = 0.0;
			Real sum_squares = 0.0;
			Real coherent_gain = 0.0;
			Real rms_gain = 0.0;
			Real enbw_bins = 0.0;
		};

		inline WindowMetrics Metrics(const Vector<Real>& window) {
			FourierValidation::ValidateNonEmpty(window.size(), "Windows::Metrics");

			WindowMetrics metrics;
			for (int i = 0; i < window.size(); i++) {
				metrics.sum += window[i];
				metrics.sum_squares += window[i] * window[i];
			}

			if (metrics.sum == 0.0)
				throw FourierError("Windows::Metrics: window sum must be non-zero");

			Real n = static_cast<Real>(window.size());
			metrics.coherent_gain = metrics.sum / n;
			metrics.rms_gain = std::sqrt(metrics.sum_squares / n);
			metrics.enbw_bins = n * metrics.sum_squares / (metrics.sum * metrics.sum);

			return metrics;
		}
	} // namespace Windows
} // namespace MML::Fourier

#endif // MML_WINDOWING_H
