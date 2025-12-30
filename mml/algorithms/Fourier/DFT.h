///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DFT.h                                                               ///
///  Description: Discrete Fourier Transform (reference O(n²) implementation)         ///
///               Forward and inverse DFT for complex signals                         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DFT_H
#define MML_DFT_H

#include "MMLBase.h"

#include "base/Vector.h"

namespace MML
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    // DFT - Discrete Fourier Transform
    //
    // Naive O(n²) implementation of DFT for reference and testing purposes.
    // 
    // Forward transform:  X[k] = Σ_{n=0}^{N-1} x[n] * exp(-2πi*k*n/N)
    // Inverse transform:  x[n] = (1/N) * Σ_{k=0}^{N-1} X[k] * exp(2πi*k*n/N)
    //
    // Features:
    // - Works for any size (not just power of 2)
    // - Reference for FFT correctness testing
    // - Educational value - clear implementation of DFT formula
    //
    // WARNING: This is O(n²) - use FFT for production code!
    ///////////////////////////////////////////////////////////////////////////////////////////
    class DFT
    {
    public:
        // Forward DFT: time domain -> frequency domain
        static Vector<Complex> Forward(const Vector<Complex>& data)
        {
            int N = data.size();
            if (N == 0) {
                throw std::invalid_argument("DFT::Forward - empty input vector");
            }

            Vector<Complex> result(N);
            const Real two_pi = 2.0 * Constants::PI;

            for (int k = 0; k < N; k++)
            {
                Complex sum(0.0, 0.0);
                for (int n = 0; n < N; n++)
                {
                    Real angle = -two_pi * k * n / N;
                    Complex twiddle(std::cos(angle), std::sin(angle));
                    sum += data[n] * twiddle;
                }
                result[k] = sum;
            }

            return result;
        }

        // Forward DFT for real-valued input (converts to complex internally)
        static Vector<Complex> Forward(const Vector<Real>& data)
        {
            int N = data.size();
            if (N == 0) {
                throw std::invalid_argument("DFT::Forward - empty input vector");
            }

            Vector<Complex> complex_data(N);
            for (int i = 0; i < N; i++) {
                complex_data[i] = Complex(data[i], 0.0);
            }

            return Forward(complex_data);
        }

        // Inverse DFT: frequency domain -> time domain
        static Vector<Complex> Inverse(const Vector<Complex>& spectrum)
        {
            int N = spectrum.size();
            if (N == 0) {
                throw std::invalid_argument("DFT::Inverse - empty input vector");
            }

            Vector<Complex> result(N);
            const Real two_pi = 2.0 * Constants::PI;

            for (int n = 0; n < N; n++)
            {
                Complex sum(0.0, 0.0);
                for (int k = 0; k < N; k++)
                {
                    Real angle = two_pi * k * n / N;  // Positive sign for inverse
                    Complex twiddle(std::cos(angle), std::sin(angle));
                    sum += spectrum[k] * twiddle;
                }
                result[n] = sum / static_cast<Real>(N);  // Normalization by 1/N
            }

            return result;
        }

        // Extract real part from inverse DFT (for real-valued signals)
        static Vector<Real> InverseReal(const Vector<Complex>& spectrum)
        {
            Vector<Complex> complex_result = Inverse(spectrum);
            int N = complex_result.size();
            
            Vector<Real> result(N);
            for (int i = 0; i < N; i++) {
                result[i] = complex_result[i].real();
            }

            return result;
        }

        // Utility: Check if DFT size is reasonable (warn if too large)
        static bool IsReasonableSize(int n)
        {
            // DFT is O(n²), so even moderate sizes get slow
            // n=1024 -> ~1M operations
            // n=4096 -> ~16M operations
            return n <= 4096;
        }

        // Utility: Compute single frequency bin (useful for testing)
        static Complex ComputeBin(const Vector<Complex>& data, int k)
        {
            int N = data.size();
            if (k < 0 || k >= N) {
                throw std::invalid_argument("DFT::ComputeBin - invalid frequency bin index");
            }

            Complex sum(0.0, 0.0);
            const Real two_pi = 2.0 * Constants::PI;

            for (int n = 0; n < N; n++)
            {
                Real angle = -two_pi * k * n / N;
                Complex twiddle(std::cos(angle), std::sin(angle));
                sum += data[n] * twiddle;
            }

            return sum;
        }
    };

} // namespace MML

#endif // MML_DFT_H
