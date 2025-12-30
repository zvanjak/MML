///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Correlation.h                                                       ///
///  Description: Cross-correlation and auto-correlation via FFT                      ///
///               Signal similarity analysis and lag detection                        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_CORRELATION_H
#define MML_CORRELATION_H

#include "FFT.h"

#include "base/Vector.h"

#include <complex>
#include <algorithm>
#include <stdexcept>

namespace MML
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    // Correlation - Fast correlation using FFT
    //
    // CORRELATION THEOREM:
    //   corr(x,y)[k] = Σ x[i] * conj(y[i+k])  (discrete cross-correlation)
    //   corr(x,y) = IFFT(FFT(x) ⊙ conj(FFT(y)))
    //
    // COMPLEXITY: O(n log n) using FFT vs O(n²) naive
    //
    // APPLICATIONS:
    //   - Signal processing (template matching, echo detection)
    //   - Time series analysis (lag detection, similarity)
    //   - Image processing (pattern recognition)
    //   - Neuroscience (spike train analysis)
    ///////////////////////////////////////////////////////////////////////////////////////////
    class Correlation
    {
    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // Cross - Cross-correlation of two real signals
        //
        // OUTPUT LENGTH: n + m - 1 (same as linear convolution)
        //   - Negative lags: [0, m-2]
        //   - Zero lag: [m-1]
        //   - Positive lags: [m, n+m-2]
        //
        // ALGORITHM:
        //   Correlation = Convolution with time-reversed signal
        //   1. Reverse y → y_rev
        //   2. Convolve(x, y_rev) using FFT
        //
        // NOTE: Result is NOT symmetric unless x == y
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Vector<Real> Cross(const Vector<Real>& x, const Vector<Real>& y)
        {
            int n = x.size();
            int m = y.size();
            if (n == 0 || m == 0)
                throw std::invalid_argument("Correlation::Cross requires non-empty inputs");
            
            // Reverse y for correlation (correlation = convolution with reversed y)
            Vector<Real> y_reversed(m);
            for (int i = 0; i < m; i++)
                y_reversed[i] = y[m - 1 - i];
            
            // Output length same as linear convolution
            int result_length = n + m - 1;
            int fft_size = FFT::NextPowerOfTwo(result_length);
            
            // Convert to complex and zero-pad
            Vector<Complex> x_complex(fft_size, Complex(0, 0));
            Vector<Complex> y_complex(fft_size, Complex(0, 0));
            for (int i = 0; i < n; i++)
                x_complex[i] = Complex(x[i], 0);
            for (int i = 0; i < m; i++)
                y_complex[i] = Complex(y_reversed[i], 0);
            
            // FFT both signals
            FFT::Transform(x_complex, 1);
            FFT::Transform(y_complex, 1);
            
            // Element-wise multiply (no conjugate needed - already reversed y)
            for (int i = 0; i < fft_size; i++)
                x_complex[i] *= y_complex[i];
            
            // IFFT
            FFT::Transform(x_complex, -1);
            
            // Normalize by 1/N
            Real norm = 1.0 / fft_size;
            for (int i = 0; i < fft_size; i++)
                x_complex[i] *= norm;
            
            // Extract real part (result_length elements)
            Vector<Real> result(result_length);
            for (int i = 0; i < result_length; i++)
                result[i] = x_complex[i].real();
            
            return result;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // Cross (complex version)
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Vector<Complex> Cross(const Vector<Complex>& x, const Vector<Complex>& y)
        {
            int n = x.size();
            int m = y.size();
            if (n == 0 || m == 0)
                throw std::invalid_argument("Correlation::Cross requires non-empty inputs");
            
            // Reverse y for correlation
            Vector<Complex> y_reversed(m);
            for (int i = 0; i < m; i++)
                y_reversed[i] = y[m - 1 - i];
            
            int result_length = n + m - 1;
            int fft_size = FFT::NextPowerOfTwo(result_length);
            
            // Zero-pad
            Vector<Complex> x_padded(fft_size, Complex(0, 0));
            Vector<Complex> y_padded(fft_size, Complex(0, 0));
            for (int i = 0; i < n; i++)
                x_padded[i] = x[i];
            for (int i = 0; i < m; i++)
                y_padded[i] = y_reversed[i];
            
            // FFT
            FFT::Transform(x_padded, 1);
            FFT::Transform(y_padded, 1);
            
            // Multiply (no conjugate needed - already reversed)
            for (int i = 0; i < fft_size; i++)
                x_padded[i] *= y_padded[i];
            
            // IFFT
            FFT::Transform(x_padded, -1);
            
            // Normalize
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
        // Auto - Auto-correlation (correlation of signal with itself)
        //
        // PROPERTIES:
        //   - Symmetric around zero lag
        //   - Maximum at zero lag (lag = n-1 in output)
        //   - Decays as lag increases (for non-periodic signals)
        //
        // APPLICATIONS:
        //   - Periodicity detection
        //   - Signal energy estimation
        //   - Power spectral density (Wiener-Khinchin theorem)
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Vector<Real> Auto(const Vector<Real>& x)
        {
            return Cross(x, x);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // CrossNormalized - Normalized cross-correlation coefficient
        //
        // RANGE: [-1, 1]
        //   -1: Perfect anti-correlation
        //    0: No correlation
        //   +1: Perfect correlation
        //
        // FORMULA:
        //   ρ[k] = corr(x,y)[k] / sqrt(E_x * E_y)
        //   where E_x = Σ x[i]², E_y = Σ y[i]²
        //
        // APPLICATIONS:
        //   - Template matching (find best lag for similarity)
        //   - Signal comparison (measure similarity independent of scale)
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Vector<Real> CrossNormalized(const Vector<Real>& x, const Vector<Real>& y)
        {
            // Compute raw cross-correlation
            Vector<Real> corr = Cross(x, y);
            
            // Compute signal energies
            Real energy_x = 0.0;
            Real energy_y = 0.0;
            for (int i = 0; i < x.size(); i++)
                energy_x += x[i] * x[i];
            for (int i = 0; i < y.size(); i++)
                energy_y += y[i] * y[i];
            
            // Normalize
            Real norm = std::sqrt(energy_x * energy_y);
            if (norm > 1e-10)  // Avoid division by zero
            {
                for (int i = 0; i < corr.size(); i++)
                    corr[i] /= norm;
            }
            
            return corr;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // FindPeakLag - Find lag with maximum correlation
        //
        // RETURNS: Lag index where correlation is maximum
        //
        // NOTE: For output of length n+m-1:
        //   - Lag = 0 is at index n-1
        //   - Negative lags: indices [0, n-2]
        //   - Positive lags: indices [n, n+m-2]
        ///////////////////////////////////////////////////////////////////////////////////////////
        static int FindPeakLag(const Vector<Real>& correlation)
        {
            if (correlation.size() == 0)
                throw std::invalid_argument("Correlation::FindPeakLag requires non-empty input");
            
            int max_idx = 0;
            Real max_val = correlation[0];
            for (int i = 1; i < correlation.size(); i++)
            {
                if (correlation[i] > max_val)
                {
                    max_val = correlation[i];
                    max_idx = i;
                }
            }
            
            return max_idx;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // ConvertLagToIndex - Convert lag value to index in correlation output
        //
        // For signals x[n] and y[m]:
        //   - Output length: n + m - 1
        //   - Zero lag at: index n - 1
        //   - Negative lag k: index (n-1) - k
        //   - Positive lag k: index (n-1) + k
        ///////////////////////////////////////////////////////////////////////////////////////////
        static int ConvertLagToIndex(int lag, int x_size)
        {
            return (x_size - 1) + lag;
        }
    };

} // namespace MML

#endif // MML_CORRELATION_H
