///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DCT.h                                                               ///
///  Description: Discrete Cosine Transform (DCT-II, DCT-III, DST-I)                  ///
///               Used in signal processing, JPEG, and Chebyshev approximation        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DCT_H
#define MML_DCT_H

#include "../../MMLBase.h"

namespace MML
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    // DCT - Discrete Cosine Transform
    //
    // The DCT is widely used in signal processing, image compression (JPEG), and spectral
    // methods. It's closely related to the Discrete Fourier Transform but operates on real
    // values and has better energy compaction properties for many real-world signals.
    //
    // KEY MATHEMATICAL CONNECTION:
    // DCT-II is the transform used to compute Chebyshev polynomial coefficients when
    // sampling at Chebyshev nodes: x_k = cos(π(k+0.5)/N). This is because Chebyshev
    // polynomials T_n(cos θ) = cos(nθ) form an orthogonal basis, and the DCT implements
    // the discrete orthogonality relation.
    //
    // TYPES IMPLEMENTED:
    // - DCT-II: Standard "DCT" used in JPEG, spectral methods, Chebyshev approximation
    // - DCT-III: Inverse of DCT-II (scaled)
    // - DST-I: Discrete Sine Transform (boundary conditions for PDEs)
    //
    // APPLICATIONS:
    // - Image/audio compression (JPEG, MP3)
    // - Spectral methods (Chebyshev approximation, PDE solving)
    // - Feature extraction and dimensionality reduction
    // - Lossy data compression
    //
    // PERFORMANCE:
    // Current implementation is O(N²) naive evaluation.
    // TODO: Implement fast DCT via FFT (O(N log N)) for large N.
    //
    // REFERENCE: 
    // - Numerical Recipes 3rd Ed., Chapter 5 (Chebyshev approximation)
    // - Discrete Cosine Transform (Wikipedia)
    // - ISO/IEC JPEG standard
    ///////////////////////////////////////////////////////////////////////////////////////////
    class DCT
    {
    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // ForwardII - DCT Type-II (standard "DCT")
        //
        // FORMULA:
        //   X[k] = (2/N) * Σ_{n=0}^{N-1} x[n] * cos(π * k * (n + 0.5) / N)
        //
        // OUTPUT LENGTH: N (same as input)
        //
        // PROPERTIES:
        // - Most commonly used DCT variant
        // - Used in JPEG compression (8×8 blocks)
        // - Computes Chebyshev polynomial coefficients at Chebyshev nodes
        // - Real-to-real transform
        // - Good energy compaction (most energy in low frequencies)
        //
        // NOTE: Some definitions omit the 2/N factor or apply it differently.
        //       This implementation matches Chebyshev approximation convention.
        //
        // EXAMPLE:
        //   Vector<Real> signal = {1.0, 2.0, 3.0, 4.0};
        //   Vector<Real> dctCoeffs = DCT::ForwardII(signal);
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Vector<Real> ForwardII(const Vector<Real>& data)
        {
            int N = data.size();
            if (N == 0) {
                throw std::invalid_argument("DCT::ForwardII - empty input");
            }

            const Real pi = Constants::PI;
            Real fac = 2.0 / N;
            
            Vector<Real> result(N);
            for (int k = 0; k < N; k++)
            {
                Real sum = 0.0;
                for (int n = 0; n < N; n++)
                {
                    sum += data[n] * std::cos(pi * k * (n + 0.5) / N);
                }
                result[k] = fac * sum;
            }
            
            return result;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // InverseII - DCT Type-III (inverse of DCT-II)
        //
        // FORMULA:
        //   x[n] = X[0]/2 + Σ_{k=1}^{N-1} X[k] * cos(π * k * (n + 0.5) / N)
        //
        // OUTPUT LENGTH: N (same as input)
        //
        // PROPERTIES:
        // - Inverse transform of DCT-II
        // - Used in JPEG decompression
        // - Reconstructs signal from DCT coefficients
        // - Real-to-real transform
        //
        // NOTE: This is DCT-III scaled to be the exact inverse of our ForwardII.
        //       The first coefficient is halved to maintain orthogonality.
        //
        // EXAMPLE:
        //   Vector<Real> coeffs = DCT::ForwardII(signal);
        //   Vector<Real> recovered = DCT::InverseII(coeffs);
        //   // recovered ≈ signal (within numerical precision)
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Vector<Real> InverseII(const Vector<Real>& coefficients)
        {
            int N = coefficients.size();
            if (N == 0) {
                throw std::invalid_argument("DCT::InverseII - empty input");
            }

            const Real pi = Constants::PI;
            
            Vector<Real> result(N);
            for (int n = 0; n < N; n++)
            {
                // First term: DC component (halved for orthogonality)
                Real sum = coefficients[0] / 2.0;
                
                // Remaining terms
                for (int k = 1; k < N; k++)
                {
                    sum += coefficients[k] * std::cos(pi * k * (n + 0.5) / N);
                }
                result[n] = sum;
            }
            
            return result;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // ForwardDST - DST Type-I (Discrete Sine Transform)
        //
        // FORMULA:
        //   X[k] = sqrt(2/(N+1)) * Σ_{n=0}^{N-1} x[n] * sin(π * (k+1) * (n+1) / (N+1))
        //
        // OUTPUT LENGTH: N (same as input)
        //
        // PROPERTIES:
        // - Uses orthonormal basis (self-adjoint transform)
        // - Used for problems with Dirichlet boundary conditions (x[0] = x[N+1] = 0)
        // - Common in PDE solving (finite difference methods)
        // - Sine basis instead of cosine
        // - Real-to-real transform
        //
        // APPLICATIONS:
        // - Solving Poisson equation with boundary conditions
        // - Heat equation on finite domains
        // - Wave equation discretization
        //
        // EXAMPLE:
        //   Vector<Real> u = {0.5, 1.0, 0.5};  // Interior points
        //   Vector<Real> dstCoeffs = DCT::ForwardDST(u);
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Vector<Real> ForwardDST(const Vector<Real>& data)
        {
            int N = data.size();
            if (N == 0) {
                throw std::invalid_argument("DCT::ForwardDST - empty input");
            }

            const Real pi = Constants::PI;
            Real norm = std::sqrt(2.0 / (N + 1));  // Orthonormal basis
            
            Vector<Real> result(N);
            for (int k = 0; k < N; k++)
            {
                Real sum = 0.0;
                for (int n = 0; n < N; n++)
                {
                    sum += data[n] * std::sin(pi * (k + 1) * (n + 1) / (N + 1));
                }
                result[k] = norm * sum;
            }
            
            return result;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // InverseDST - DST Type-I inverse
        //
        // FORMULA:
        //   x[n] = sqrt(2/(N+1)) * Σ_{k=0}^{N-1} X[k] * sin(π * (k+1) * (n+1) / (N+1))
        //
        // PROPERTIES:
        // - DST-I with orthonormal basis is self-adjoint (its own inverse)
        // - Used to reconstruct signal from DST coefficients
        //
        // EXAMPLE:
        //   Vector<Real> coeffs = DCT::ForwardDST(data);
        //   Vector<Real> recovered = DCT::InverseDST(coeffs);
        //   // recovered ≈ data
        ///////////////////////////////////////////////////////////////////////////////////////////
        static Vector<Real> InverseDST(const Vector<Real>& coefficients)
        {
            // DST-I with orthonormal basis is self-adjoint - inverse is the same as forward
            return ForwardDST(coefficients);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // HELPER: Verify DCT-II/III round-trip
        //
        // Tests that InverseII(ForwardII(x)) ≈ x within tolerance
        ///////////////////////////////////////////////////////////////////////////////////////////
        static bool VerifyRoundTrip(const Vector<Real>& data, Real tolerance = 1e-12)
        {
            Vector<Real> coeffs = ForwardII(data);
            Vector<Real> recovered = InverseII(coeffs);
            
            if (recovered.size() != data.size()) return false;
            
            for (int i = 0; i < data.size(); i++)
            {
                if (std::abs(recovered[i] - data[i]) > tolerance)
                    return false;
            }
            return true;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////
        // FUTURE ENHANCEMENTS:
        //
        // 1. Fast DCT via FFT (O(N log N))
        //    DCT can be computed using FFT with pre/post-processing
        //    See: Makhoul, "A fast cosine transform in one and two dimensions" (1980)
        //
        // 2. 2D DCT for image compression
        //    static Matrix<Real> Forward2D(const Matrix<Real>& image);
        //    static Matrix<Real> Inverse2D(const Matrix<Real>& dctImage);
        //    Apply 1D DCT to rows, then columns (separable transform)
        //
        // 3. Other DCT types (DCT-IV used in MDCT for audio)
        //    static Vector<Real> ForwardIV(const Vector<Real>& data);
        //
        // 4. Block DCT (8×8 blocks for JPEG)
        //    static Matrix<Real> BlockDCT(const Matrix<Real>& image, int blockSize = 8);
        ///////////////////////////////////////////////////////////////////////////////////////////
    };

} // namespace MML

#endif // MML_DCT_H