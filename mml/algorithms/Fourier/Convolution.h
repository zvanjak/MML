///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Convolution.h                                                       ///
///  Description: Fast convolution and deconvolution via FFT                          ///
///               Circular and linear convolution, Wiener deconvolution               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CONVOLUTION_H
#define MML_CONVOLUTION_H

#include "../../MMLBase.h"
#include "FFT.h"

namespace MML
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    // Convolution - Fast convolution and deconvolution via FFT
    //
    // CONVOLUTION THEOREM:
    //   x * y = IFFT(FFT(x) ⊙ FFT(y))
    //
    // PERFORMANCE:
    //   - Direct convolution: O(n²)
    //   - FFT convolution: O(n log n)
    //   - Speedup: ~100x for n=1024
    //
    // REFERENCE: Numerical Recipes convlv.h
    ///////////////////////////////////////////////////////////////////////////////////////////

    class Convolution
    {
    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // Linear - Linear convolution (standard discrete convolution)
        //
        // OUTPUT LENGTH: n + m - 1
        //   - For signal of length n and kernel of length m
        //   - Produces complete convolution with no wraparound
        //
        // ALGORITHM:
        //   1. Zero-pad both to next power of 2 >= n+m-1
        //   2. FFT both
        //   3. Element-wise multiply
        //   4. IFFT
        //   5. Extract result (first n+m-1 samples)
        //
        // APPLICATIONS:
        //   - Filtering (FIR filters)
        //   - Smoothing
        //   - Edge detection
        //   - Pattern matching
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Vector<Real> Linear(const Vector<Real>& signal, const Vector<Real>& kernel)
        {
            int n = signal.size();
            int m = kernel.size();
            if (n == 0 || m == 0)
                throw std::invalid_argument("Convolution::Linear requires non-empty inputs");
            
            int result_length = n + m - 1;
            int fft_size = FFT::NextPowerOfTwo(result_length);
            
            // Convert to complex and zero-pad
            Vector<Complex> signal_complex(fft_size, Complex(0, 0));
            Vector<Complex> kernel_complex(fft_size, Complex(0, 0));
            for (int i = 0; i < n; i++)
                signal_complex[i] = Complex(signal[i], 0);
            for (int i = 0; i < m; i++)
                kernel_complex[i] = Complex(kernel[i], 0);
            
            // FFT both
            FFT::Transform(signal_complex, 1);
            FFT::Transform(kernel_complex, 1);
            
            // Element-wise multiply
            for (int i = 0; i < fft_size; i++)
                signal_complex[i] *= kernel_complex[i];
            
            // IFFT
            FFT::Transform(signal_complex, -1);
            
            // Normalize by 1/N
            Real norm = 1.0 / fft_size;
            for (int i = 0; i < fft_size; i++)
                signal_complex[i] *= norm;
            
            // Extract real part (first n+m-1 samples)
            Vector<Real> result(result_length);
            for (int i = 0; i < result_length; i++)
                result[i] = signal_complex[i].real();
            
            return result;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // Linear (complex version)
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Vector<Complex> Linear(const Vector<Complex>& x, const Vector<Complex>& y)
        {
            int n = x.size();
            int m = y.size();
            if (n == 0 || m == 0)
                throw std::invalid_argument("Convolution::Linear requires non-empty inputs");
            
            int result_length = n + m - 1;
            int fft_size = FFT::NextPowerOfTwo(result_length);
            
            // Zero-pad
            Vector<Complex> x_padded(fft_size, Complex(0, 0));
            Vector<Complex> y_padded(fft_size, Complex(0, 0));
            for (int i = 0; i < n; i++)
                x_padded[i] = x[i];
            for (int i = 0; i < m; i++)
                y_padded[i] = y[i];
            
            // FFT
            FFT::Transform(x_padded, 1);
            FFT::Transform(y_padded, 1);
            
            // Multiply
            for (int i = 0; i < fft_size; i++)
                x_padded[i] *= y_padded[i];
            
            // IFFT
            FFT::Transform(x_padded, -1);
            
            // Normalize by 1/N
            Real norm = 1.0 / fft_size;
            for (int i = 0; i < fft_size; i++)
                x_padded[i] *= norm;
            
            // Extract result
            Vector<Complex> result(result_length);
            for (int i = 0; i < result_length; i++)
                result[i] = x_padded[i];
            
            return result;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // Circular - Circular convolution (periodic convolution)
        //
        // OUTPUT LENGTH: max(n, m)
        //   - Same length as longer input
        //   - Wraparound at boundaries (periodic)
        //
        // ALGORITHM:
        //   1. Pad both to same length (next power of 2)
        //   2. FFT both
        //   3. Element-wise multiply
        //   4. IFFT
        //
        // APPLICATIONS:
        //   - Periodic signals
        //   - Spectral analysis
        //   - Fast polynomial multiplication (mod x^n - 1)
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Vector<Real> Circular(const Vector<Real>& x, const Vector<Real>& y)
        {
            int n = x.size();
            int m = y.size();
            if (n == 0 || m == 0)
                throw std::invalid_argument("Convolution::Circular requires non-empty inputs");
            
            int max_len = std::max(n, m);
            int fft_size = FFT::NextPowerOfTwo(max_len);
            
            // Convert to complex and zero-pad
            Vector<Complex> x_complex(fft_size, Complex(0, 0));
            Vector<Complex> y_complex(fft_size, Complex(0, 0));
            for (int i = 0; i < n; i++)
                x_complex[i] = Complex(x[i], 0);
            for (int i = 0; i < m; i++)
                y_complex[i] = Complex(y[i], 0);
            
            // FFT both
            FFT::Transform(x_complex, 1);
            FFT::Transform(y_complex, 1);
            
            // Element-wise multiply
            for (int i = 0; i < fft_size; i++)
                x_complex[i] *= y_complex[i];
            
            // IFFT
            FFT::Transform(x_complex, -1);
            
            // Normalize by 1/N
            Real norm = 1.0 / fft_size;
            for (int i = 0; i < fft_size; i++)
                x_complex[i] *= norm;
            
            // Extract real part (length = max_len)
            Vector<Real> result(max_len);
            for (int i = 0; i < max_len; i++)
                result[i] = x_complex[i].real();
            
            return result;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // Deconvolve - Inverse filtering (deconvolution)
        //
        // PROBLEM: Given y = x * h, recover x given y and h
        // SOLUTION: x = IFFT(FFT(y) / FFT(h))
        //
        // WARNING:
        //   - Unstable if kernel has small frequency components
        //   - May amplify noise
        //   - Consider regularization for real applications
        //
        // PARAMETERS:
        //   signal - Convolved signal y
        //   kernel - Known kernel h
        //   epsilon - Regularization parameter (prevents division by zero)
        //
        // APPLICATIONS:
        //   - Image deblurring
        //   - Channel equalization
        //   - System identification
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Vector<Real> Deconvolve(const Vector<Real>& signal, 
                                      const Vector<Real>& kernel,
                                      Real epsilon = 1e-10)
        {
            int n = signal.size();
            int m = kernel.size();
            if (n == 0 || m == 0)
                throw std::invalid_argument("Convolution::Deconvolve requires non-empty inputs");
            
            // Use FFT size matching signal length
            int fft_size = FFT::NextPowerOfTwo(n);
            
            // Convert to complex and zero-pad
            Vector<Complex> signal_complex(fft_size, Complex(0, 0));
            Vector<Complex> kernel_complex(fft_size, Complex(0, 0));
            for (int i = 0; i < n; i++)
                signal_complex[i] = Complex(signal[i], 0);
            for (int i = 0; i < m; i++)
                kernel_complex[i] = Complex(kernel[i], 0);
            
            // FFT both
            FFT::Transform(signal_complex, 1);
            FFT::Transform(kernel_complex, 1);
            
            // Divide with regularization
            for (int i = 0; i < fft_size; i++)
            {
                Real kernel_mag_sq = std::norm(kernel_complex[i]);
                if (kernel_mag_sq > epsilon)
                    signal_complex[i] /= kernel_complex[i];
                else
                    signal_complex[i] = Complex(0, 0);  // Zero out unstable frequencies
            }
            
            // IFFT
            FFT::Transform(signal_complex, -1);
            
            // Normalize by 1/N
            Real norm = 1.0 / fft_size;
            for (int i = 0; i < fft_size; i++)
                signal_complex[i] *= norm;
            
            // Extract real part
            Vector<Real> result(n);
            for (int i = 0; i < n; i++)
                result[i] = signal_complex[i].real();
            
            return result;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // Deconvolve (complex version)
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Vector<Complex> Deconvolve(const Vector<Complex>& signal,
                                         const Vector<Complex>& kernel,
                                         Real epsilon = 1e-10)
        {
            int n = signal.size();
            int m = kernel.size();
            if (n == 0 || m == 0)
                throw std::invalid_argument("Convolution::Deconvolve requires non-empty inputs");
            
            int fft_size = FFT::NextPowerOfTwo(n);
            
            // Zero-pad
            Vector<Complex> signal_padded(fft_size, Complex(0, 0));
            Vector<Complex> kernel_padded(fft_size, Complex(0, 0));
            for (int i = 0; i < n; i++)
                signal_padded[i] = signal[i];
            for (int i = 0; i < m; i++)
                kernel_padded[i] = kernel[i];
            
            // FFT
            FFT::Transform(signal_padded, 1);
            FFT::Transform(kernel_padded, 1);
            
            // Divide with regularization
            for (int i = 0; i < fft_size; i++)
            {
                Real kernel_mag_sq = std::norm(kernel_padded[i]);
                if (kernel_mag_sq > epsilon)
                    signal_padded[i] /= kernel_padded[i];
                else
                    signal_padded[i] = Complex(0, 0);
            }
            
            // IFFT
            FFT::Transform(signal_padded, -1);
            
            // Normalize by 1/N
            Real norm = 1.0 / fft_size;
            for (int i = 0; i < fft_size; i++)
                signal_padded[i] *= norm;
            
            // Extract result
            Vector<Complex> result(n);
            for (int i = 0; i < n; i++)
                result[i] = signal_padded[i];
            
            return result;
        }
    };

} // namespace MML

#endif // MML_CONVOLUTION_H
