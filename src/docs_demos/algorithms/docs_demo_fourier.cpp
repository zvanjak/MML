#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector.h"
#include "algorithms/Fourier/DFT.h"
#include "algorithms/Fourier/FFT.h"
#include "algorithms/Fourier/RealFFT.h"
#include "algorithms/Fourier/Spectrum.h"
#include "algorithms/Fourier/Convolution.h"
#include "algorithms/Fourier/Windowing.h"
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                    DISCRETE FOURIER TRANSFORM (DFT)                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Fourier_DFT_Basic()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Discrete Fourier Transform (DFT)\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nDFT: X[k] = Σ x[n] * exp(-2πi*k*n/N)\n";
	std::cout << "Reference O(n²) implementation - use FFT for production!\n";
	
	// Simple signal: DC + single frequency
	std::cout << "\n--- Pure Sine Wave at 2 Hz ---\n";
	int N = 8;
	Real fs = 8.0;  // Sampling frequency
	Real f0 = 2.0;  // Signal frequency
	
	Vector<Complex> signal(N);
	for (int n = 0; n < N; n++) {
		Real t = n / fs;
		signal[n] = Complex(std::sin(2.0 * Constants::PI * f0 * t), 0.0);
	}
	
	std::cout << "Signal (time domain): ";
	for (int n = 0; n < N; n++)
		std::cout << signal[n].real() << " ";
	std::cout << std::endl;
	
	// Forward DFT
	Vector<Complex> spectrum = DFT::Forward(signal);
	
	std::cout << "\nSpectrum magnitudes: ";
	for (int k = 0; k < N; k++)
		std::cout << std::abs(spectrum[k]) << " ";
	std::cout << std::endl;
	
	std::cout << "\nNote: Peak at k=2 (and N-2=6) corresponds to 2 Hz\n";
	
	// Inverse DFT
	Vector<Complex> recovered = DFT::Inverse(spectrum);
	
	std::cout << "\n--- Round-Trip Test ---\n";
	Real maxError = 0;
	for (int n = 0; n < N; n++)
		maxError = std::max(maxError, std::abs(signal[n] - recovered[n]));
	std::cout << "Max reconstruction error: " << maxError << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    FAST FOURIER TRANSFORM (FFT)                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Fourier_FFT_Basic()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Fast Fourier Transform (FFT)\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nFFT: Cooley-Tukey radix-2 algorithm, O(n log n)\n";
	std::cout << "Requires input size to be power of 2.\n";
	
	// Create a test signal: sum of two sinusoids
	int N = 64;
	Real fs = 64.0;
	
	Vector<Complex> signal(N);
	for (int n = 0; n < N; n++) {
		Real t = n / fs;
		// 5 Hz + 12 Hz components
		signal[n] = Complex(
			std::sin(2.0 * Constants::PI * 5.0 * t) + 
			0.5 * std::sin(2.0 * Constants::PI * 12.0 * t), 
			0.0
		);
	}
	
	std::cout << "\n--- Signal: sin(2π*5t) + 0.5*sin(2π*12t) ---\n";
	
	// FFT
	Vector<Complex> spectrum = FFT::Forward(signal);
	
	std::cout << "Significant frequency components:\n";
	for (int k = 0; k < N/2; k++) {
		Real mag = std::abs(spectrum[k]) / (N/2);
		if (mag > 0.1) {
			Real freq = k * fs / N;
			std::cout << "  f = " << freq << " Hz, amplitude = " << mag << std::endl;
		}
	}
	
	// Inverse FFT
	Vector<Complex> recovered = FFT::Inverse(spectrum);
	
	Real maxError = 0;
	for (int n = 0; n < N; n++)
		maxError = std::max(maxError, std::abs(signal[n] - recovered[n]));
	std::cout << "\nReconstruction error: " << maxError << std::endl;
}

void Docs_Demo_Fourier_FFT_vs_DFT()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: FFT vs DFT Performance\n";
	std::cout << "==========================================================================\n";
	
	// Create test signal
	int N = 256;
	Vector<Complex> signal(N);
	for (int n = 0; n < N; n++)
		signal[n] = Complex(std::sin(2.0 * Constants::PI * n / N), 0.0);
	
	std::cout << "\nComparing FFT and DFT results for N = " << N << "\n";
	
	// Both transforms
	Vector<Complex> fft_result = FFT::Forward(signal);
	Vector<Complex> dft_result = DFT::Forward(signal);
	
	// Compare
	Real maxDiff = 0;
	for (int k = 0; k < N; k++)
		maxDiff = std::max(maxDiff, std::abs(fft_result[k] - dft_result[k]));
	
	std::cout << "Max difference between FFT and DFT: " << maxDiff << std::endl;
	std::cout << "FFT is O(n log n) = O(" << N * std::log2(N) << ")\n";
	std::cout << "DFT is O(n²) = O(" << N*N << ")\n";
	std::cout << "Speedup factor: ~" << (N*N)/(N * std::log2(N)) << "x\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    REAL-VALUED FFT                                                  ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Fourier_RealFFT()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Real-Valued FFT\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nRealFFT exploits conjugate symmetry for real signals.\n";
	std::cout << "Twice as efficient as complex FFT for real data.\n";
	
	// Real signal
	int N = 32;
	Vector<Real> signal(N);
	for (int n = 0; n < N; n++)
		signal[n] = std::sin(2.0 * Constants::PI * 3.0 * n / N) + 0.5;  // Sine + DC offset
	
	std::cout << "\n--- Real signal: sin(2π*3n/N) + 0.5 ---\n";
	
	// Forward transform
	Vector<Complex> spectrum = RealFFT::Forward(signal);
	
	std::cout << "Spectrum (first half + Nyquist):\n";
	for (int k = 0; k <= N/2; k++) {
		Real mag = std::abs(spectrum[k]) / N;
		if (mag > 0.01)
			std::cout << "  k=" << k << ": magnitude = " << mag << std::endl;
	}
	
	// Inverse
	Vector<Real> recovered = RealFFT::Inverse(spectrum);
	
	Real maxError = 0;
	for (int n = 0; n < N; n++)
		maxError = std::max(maxError, std::abs(signal[n] - recovered[n]));
	std::cout << "\nReconstruction error: " << maxError << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    SPECTRUM ANALYSIS                                                ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Fourier_Spectrum()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Spectrum Analysis\n";
	std::cout << "==========================================================================\n";
	
	// Create signal with known components
	int N = 128;
	Real fs = 1000.0;  // 1 kHz sampling rate
	
	Vector<Real> signal(N);
	for (int n = 0; n < N; n++) {
		Real t = n / fs;
		// 50 Hz, 120 Hz, 200 Hz components
		signal[n] = std::sin(2.0 * Constants::PI * 50.0 * t) +
		            0.5 * std::sin(2.0 * Constants::PI * 120.0 * t) +
		            0.3 * std::sin(2.0 * Constants::PI * 200.0 * t);
	}
	
	std::cout << "\n--- Signal: sin(50t) + 0.5*sin(120t) + 0.3*sin(200t) ---\n";
	std::cout << "Sampling rate: " << fs << " Hz, Duration: " << N/fs*1000 << " ms\n";
	
	// Compute spectrum
	Vector<Complex> fft = RealFFT::Forward(signal);
	
	// Find peaks
	std::cout << "\nDetected frequency components:\n";
	for (int k = 1; k < N/2; k++) {
		Real mag = 2.0 * std::abs(fft[k]) / N;
		Real freq = k * fs / N;
		if (mag > 0.2)  // Threshold
			std::cout << "  f = " << freq << " Hz, amplitude = " << mag << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    WINDOWING                                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Fourier_Windowing()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Window Functions\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nWindowing reduces spectral leakage in FFT analysis.\n";
	std::cout << "Different windows trade main lobe width vs sidelobe level.\n";
	
	int N = 16;
	
	// Generate windows
	auto hann = Windows::Hann(N);
	auto hamming = Windows::Hamming(N);
	auto blackman = Windows::Blackman(N);
	
	std::cout << "\n--- Window Functions (N=" << N << ") ---\n";
	std::cout << "n\tHann\t\tHamming\t\tBlackman\n";
	for (int n = 0; n < N; n += 2) {
		std::cout << n << "\t" 
		          << hann[n] << "\t" 
		          << hamming[n] << "\t" 
		          << blackman[n] << std::endl;
	}
	
	std::cout << "\nWindow properties:\n";
	std::cout << "  Hann:     Smooth, good sidelobe suppression (-31 dB)\n";
	std::cout << "  Hamming:  Optimized for minimum sidelobe level (-43 dB)\n";
	std::cout << "  Blackman: Excellent sidelobe suppression (-57 dB)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    CONVOLUTION                                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Fourier_Convolution()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: FFT-Based Convolution\n";
	std::cout << "==========================================================================\n";
	
	std::cout << "\nConvolution theorem: conv(f,g) = IFFT(FFT(f) * FFT(g))\n";
	std::cout << "O(n log n) vs O(n²) for direct convolution.\n";
	
	// Signal
	Vector<Real> signal({1, 2, 3, 4, 5, 4, 3, 2});
	
	// Simple moving average filter (kernel)
	Vector<Real> kernel({0.25, 0.5, 0.25});
	
	std::cout << "\n--- Signal: [1,2,3,4,5,4,3,2] ---\n";
	std::cout << "--- Kernel: [0.25, 0.5, 0.25] (smoothing) ---\n";
	
	// FFT-based linear convolution
	auto result = Convolution::Linear(signal, kernel);
	
	std::cout << "\nConvolved result:\n  [";
	for (int i = 0; i < result.size(); i++)
		std::cout << result[i] << (i < result.size()-1 ? ", " : "");
	std::cout << "]\n";
	
	std::cout << "\nNote: Output length = len(signal) + len(kernel) - 1\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         MAIN DEMO FUNCTION                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Fourier()
{
	std::cout << "\n##########################################################################\n";
	std::cout << "#                   FOURIER TRANSFORM DEMOS                              #\n";
	std::cout << "##########################################################################\n";
	
	Docs_Demo_Fourier_DFT_Basic();
	Docs_Demo_Fourier_FFT_Basic();
	Docs_Demo_Fourier_FFT_vs_DFT();
	Docs_Demo_Fourier_RealFFT();
	Docs_Demo_Fourier_Spectrum();
	Docs_Demo_Fourier_Windowing();
	Docs_Demo_Fourier_Convolution();
	
	std::cout << "\n=== All Fourier Transform Demos Complete ===\n";
}
