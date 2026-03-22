///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Fourier.h                                                           ///
///  Description: Fourier transforms for MML core algorithms layer                    ///
///               DFT  - Naive O(n²) reference implementation                        ///
///               FFT  - Cooley-Tukey radix-2 O(n log n)                             ///
///               DCT  - Discrete Cosine Transform (DCT-II/III, DST-I)               ///
///               Fast O(N log N) DCT via FFT for power-of-2 sizes                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                             ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FOURIER_H
#define MML_FOURIER_H

#include "../MMLBase.h"
#include "../base/Vector/Vector.h"
#include "../core/AlgorithmTypes.h"

#include <stdexcept>
#include <string>

namespace MML::Fourier 
{
	///////////////////////////////////////////////////////////////////////////////////////////
	// FourierConfig - Configuration for Fourier transform detailed APIs
	///////////////////////////////////////////////////////////////////////////////////////////
	struct FourierConfig : public EvaluationConfigBase {
		// Inherits: estimate_error, check_finite, exception_policy
		// No additional parameters needed for Fourier transforms
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// FourierResult - Result type for Fourier transform detailed APIs
	///////////////////////////////////////////////////////////////////////////////////////////
	template<typename T>
	struct FourierResult : public EvaluationResultBase {
		/// The computed transform output
		T value{};

		/// Input size that was processed
		int input_size = 0;
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// FourierDetail - Internal helpers for Detailed API execution
	///////////////////////////////////////////////////////////////////////////////////////////
	namespace FourierDetail
	{
		/// Execute a Fourier Detailed operation with timing and exception handling.
		template<typename ResultType, typename ComputeFn>
		ResultType ExecuteFourierDetailed(const char* algorithm_name,
		                                 const FourierConfig& config,
		                                 ComputeFn&& compute)
		{
			auto execute = [&]() {
				AlgorithmTimer timer;

				ResultType result = MakeEvaluationSuccessResult<ResultType>(algorithm_name);

				compute(result);

				result.elapsed_time_ms = timer.elapsed_ms();
				return result;
			};

			if (config.exception_policy == EvaluationExceptionPolicy::Propagate)
				return execute();

			try {
				return execute();
			}
			catch (const std::invalid_argument& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::InvalidInput, ex.what(), algorithm_name);
			}
			catch (const std::exception& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::AlgorithmSpecificFailure, ex.what(), algorithm_name);
			}
		}
	} // namespace FourierDetail

	///////////////////////////////////////////////////////////////////////////////////////////
	// FourierValidation - Reusable validation helpers for Fourier transform functions
	///////////////////////////////////////////////////////////////////////////////////////////
	namespace FourierValidation 
	{
		///////////////////////////////////////////////////////////////////////////////////////////
		// Basic size validations
		///////////////////////////////////////////////////////////////////////////////////////////

		/// @brief Validate that a container is not empty
		inline void ValidateNonEmpty(int size, const std::string& functionName) {
			if (size == 0) {
				throw std::invalid_argument(functionName + " - empty input");
			}
		}

		/// @brief Validate that a size is positive (> 0)
		inline void ValidatePositiveSize(int size, const std::string& paramName, const std::string& functionName) {
			if (size <= 0) {
				throw std::invalid_argument(functionName + " - " + paramName + " must be positive");
			}
		}

		/// @brief Validate that a size meets a minimum requirement
		inline void ValidateMinSize(int size, int minSize, const std::string& functionName) {
			if (size < minSize) {
				throw std::invalid_argument(functionName + " - size must be at least " + std::to_string(minSize));
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Power-of-2 validations (critical for FFT)
		///////////////////////////////////////////////////////////////////////////////////////////

		inline bool IsPowerOfTwo(int n) { return (n > 0) && ((n & (n - 1)) == 0); }

		inline void ValidatePowerOfTwo(int size, const std::string& functionName) {
			if (!IsPowerOfTwo(size)) {
				throw std::invalid_argument(functionName + " - size must be power of 2");
			}
		}

		inline void ValidateFFTSize(int size, const std::string& functionName) {
			ValidateNonEmpty(size, functionName);
			ValidatePowerOfTwo(size, functionName);
		}

		inline void ValidateRealFFTSize(int size, int minSize, const std::string& functionName) {
			if (size < minSize || !IsPowerOfTwo(size)) {
				throw std::invalid_argument(functionName + " - requires power-of-2 size >= " + std::to_string(minSize));
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Window validations
		///////////////////////////////////////////////////////////////////////////////////////////

		inline void ValidateWindowSize(int n) {
			if (n <= 0) {
				throw std::invalid_argument("Window size must be positive");
			}
		}

		inline void ValidateWindowMatch(int windowSize, int dataSize, const std::string& functionName) {
			if (windowSize != dataSize) {
				throw std::invalid_argument(functionName + " - window size mismatch");
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Parameter validations
		///////////////////////////////////////////////////////////////////////////////////////////

		template<typename T>
		inline void ValidatePositive(T value, const std::string& paramName, const std::string& functionName) {
			if (value <= T(0)) {
				throw std::invalid_argument(functionName + " - " + paramName + " must be positive");
			}
		}

		template<typename T>
		inline void ValidateNonNegative(T value, const std::string& paramName, const std::string& functionName) {
			if (value < T(0)) {
				throw std::invalid_argument(functionName + " - " + paramName + " must be non-negative");
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Index validations
		///////////////////////////////////////////////////////////////////////////////////////////

		inline void ValidateIndex(int index, int maxIndex, const std::string& functionName) {
			if (index < 0 || index >= maxIndex) {
				throw std::invalid_argument(functionName + " - invalid index");
			}
		}

	} // namespace FourierValidation

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
	class DFT {
	public:
		// Forward DFT: time domain -> frequency domain
		static Vector<Complex> Forward(const Vector<Complex>& data) {
			int N = data.size();
			FourierValidation::ValidateNonEmpty(N, "DFT::Forward");

			Vector<Complex> result(N);
			const Real two_pi = 2.0 * Constants::PI;

			for (int k = 0; k < N; k++) {
				Complex sum(0.0, 0.0);
				for (int n = 0; n < N; n++) {
					Real angle = -two_pi * k * n / N;
					Complex twiddle(std::cos(angle), std::sin(angle));
					sum += data[n] * twiddle;
				}
				result[k] = sum;
			}

			return result;
		}

		// Forward DFT for real-valued input (converts to complex internally)
		static Vector<Complex> Forward(const Vector<Real>& data) {
			int N = data.size();
			FourierValidation::ValidateNonEmpty(N, "DFT::Forward");

			Vector<Complex> complex_data(N);
			for (int i = 0; i < N; i++) {
				complex_data[i] = Complex(data[i], 0.0);
			}

			return Forward(complex_data);
		}

		// Inverse DFT: frequency domain -> time domain
		static Vector<Complex> Inverse(const Vector<Complex>& spectrum) {
			int N = spectrum.size();
			FourierValidation::ValidateNonEmpty(N, "DFT::Inverse");

			Vector<Complex> result(N);
			const Real two_pi = 2.0 * Constants::PI;

			for (int n = 0; n < N; n++) {
				Complex sum(0.0, 0.0);
				for (int k = 0; k < N; k++) {
					Real angle = two_pi * k * n / N; // Positive sign for inverse
					Complex twiddle(std::cos(angle), std::sin(angle));
					sum += spectrum[k] * twiddle;
				}
				result[n] = sum / static_cast<Real>(N); // Normalization by 1/N
			}

			return result;
		}

		// Extract real part from inverse DFT (for real-valued signals)
		static Vector<Real> InverseReal(const Vector<Complex>& spectrum) {
			Vector<Complex> complex_result = Inverse(spectrum);
			int N = complex_result.size();

			Vector<Real> result(N);
			for (int i = 0; i < N; i++) {
				result[i] = complex_result[i].real();
			}

			return result;
		}

		// Utility: Check if DFT size is reasonable (warn if too large)
		static bool IsReasonableSize(int n) {
			// DFT is O(n²), so even moderate sizes get slow
			// n=1024 -> ~1M operations
			// n=4096 -> ~16M operations
			return n <= 4096;
		}

		// Utility: Compute single frequency bin (useful for testing)
		static Complex ComputeBin(const Vector<Complex>& data, int k) {
			int N = data.size();
			FourierValidation::ValidateIndex(k, N, "DFT::ComputeBin");

			Complex sum(0.0, 0.0);
			const Real two_pi = 2.0 * Constants::PI;

			for (int n = 0; n < N; n++) {
				Real angle = -two_pi * k * n / N;
				Complex twiddle(std::cos(angle), std::sin(angle));
				sum += data[n] * twiddle;
			}

			return sum;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Detailed API variants - return FourierResult with AlgorithmStatus + timing
		///////////////////////////////////////////////////////////////////////////////////////////

		/// Forward DFT with full diagnostics (complex input)
		static FourierResult<Vector<Complex>> ForwardDetailed(
			const Vector<Complex>& data, const FourierConfig& config = {})
		{
			using ResultType = FourierResult<Vector<Complex>>;
			return FourierDetail::ExecuteFourierDetailed<ResultType>(
				"DFT::Forward", config,
				[&](ResultType& result) {
					result.input_size = data.size();
					result.value = Forward(data);
					result.function_evaluations = data.size() * data.size(); // O(n²)
				});
		}

		/// Forward DFT with full diagnostics (real input)
		static FourierResult<Vector<Complex>> ForwardDetailed(
			const Vector<Real>& data, const FourierConfig& config = {})
		{
			using ResultType = FourierResult<Vector<Complex>>;
			return FourierDetail::ExecuteFourierDetailed<ResultType>(
				"DFT::Forward", config,
				[&](ResultType& result) {
					result.input_size = data.size();
					result.value = Forward(data);
					result.function_evaluations = data.size() * data.size();
				});
		}

		/// Inverse DFT with full diagnostics
		static FourierResult<Vector<Complex>> InverseDetailed(
			const Vector<Complex>& spectrum, const FourierConfig& config = {})
		{
			using ResultType = FourierResult<Vector<Complex>>;
			return FourierDetail::ExecuteFourierDetailed<ResultType>(
				"DFT::Inverse", config,
				[&](ResultType& result) {
					result.input_size = spectrum.size();
					result.value = Inverse(spectrum);
					result.function_evaluations = spectrum.size() * spectrum.size();
				});
		}

		/// Inverse DFT returning real part, with full diagnostics
		static FourierResult<Vector<Real>> InverseRealDetailed(
			const Vector<Complex>& spectrum, const FourierConfig& config = {})
		{
			using ResultType = FourierResult<Vector<Real>>;
			return FourierDetail::ExecuteFourierDetailed<ResultType>(
				"DFT::InverseReal", config,
				[&](ResultType& result) {
					result.input_size = spectrum.size();
					result.value = InverseReal(spectrum);
					result.function_evaluations = spectrum.size() * spectrum.size();
				});
		}
	};

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

			FourierValidation::ValidateFFTSize(n, "FFT::Transform");

			// Bit-reversal permutation
			BitReversalPermute(data);

			// Cooley-Tukey FFT with iterative decimation-in-time
			int mmax = 1;
			while (n > mmax) {
				int istep = mmax << 1;											// istep = 2 * mmax
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
			Transform(result, -1);						 // Inverse transform

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

		///////////////////////////////////////////////////////////////////////////////////////////
		// Detailed API variants - return FourierResult with AlgorithmStatus + timing
		///////////////////////////////////////////////////////////////////////////////////////////

		/// Forward FFT with full diagnostics
		static FourierResult<Vector<Complex>> ForwardDetailed(
			const Vector<Complex>& data, const FourierConfig& config = {})
		{
			using ResultType = FourierResult<Vector<Complex>>;
			return FourierDetail::ExecuteFourierDetailed<ResultType>(
				"FFT::Forward", config,
				[&](ResultType& result) {
					result.input_size = data.size();
					result.value = Forward(data);
					result.function_evaluations = data.size() * Log2(data.size()); // O(n log n)
				});
		}

		/// Inverse FFT with full diagnostics
		static FourierResult<Vector<Complex>> InverseDetailed(
			const Vector<Complex>& spectrum, const FourierConfig& config = {})
		{
			using ResultType = FourierResult<Vector<Complex>>;
			return FourierDetail::ExecuteFourierDetailed<ResultType>(
				"FFT::Inverse", config,
				[&](ResultType& result) {
					result.input_size = spectrum.size();
					result.value = Inverse(spectrum);
					result.function_evaluations = spectrum.size() * Log2(spectrum.size());
				});
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
	// PERFORMANCE:
	// - Fast DCT via FFT: O(N log N) for power-of-2 sizes
	// - Reference O(N²) fallback for non-power-of-2 sizes
	// - Auto-selects fastest algorithm based on input size
	//
	// APPLICATIONS:
	// - Image/audio compression (JPEG, MP3)
	// - Spectral methods (Chebyshev approximation, PDE solving)
	// - Feature extraction and dimensionality reduction
	// - Lossy data compression
	//
	// REFERENCE:
	// - Numerical Recipes 3rd Ed., Chapter 5 (Chebyshev approximation)
	// - Makhoul, "A fast cosine transform in one and two dimensions" (1980)
	// - FFTW documentation (www.fftw.org)
	///////////////////////////////////////////////////////////////////////////////////////////
	class DCT {
	private:
		/// @brief Threshold for using fast FFT-based DCT vs naive O(N²)
		/// Below this threshold, naive may be faster due to FFT overhead
		static constexpr int FAST_DCT_THRESHOLD = 32;

	public:
		///////////////////////////////////////////////////////////////////////////////////////////
		// ForwardII - DCT Type-II (standard "DCT") - AUTO-SELECTS FASTEST ALGORITHM
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> ForwardII(const Vector<Real>& data) {
			int N = data.size();
			if (N == 0) {
				throw std::invalid_argument("DCT::ForwardII - empty input");
			}

			// Small sizes: use reference implementation (lower overhead)
			// Non-power-of-2: use reference implementation
			if (N < FAST_DCT_THRESHOLD || !FFT::IsPowerOfTwo(N)) {
				return ForwardII_Reference(data);
			}

			// Fast DCT-II via FFT for power-of-2 sizes
			return ForwardII_Fast(data);
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// InverseII - DCT Type-III (inverse of DCT-II) - AUTO-SELECTS FASTEST ALGORITHM
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> InverseII(const Vector<Real>& coefficients) {
			int N = coefficients.size();
			if (N == 0) {
				throw std::invalid_argument("DCT::InverseII - empty input");
			}

			// Small sizes or non-power-of-2: use reference implementation
			if (N < FAST_DCT_THRESHOLD || !FFT::IsPowerOfTwo(N)) {
				return InverseII_Reference(coefficients);
			}

			// Fast DCT-III via FFT
			return InverseII_Fast(coefficients);
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// ForwardII_Fast - Fast DCT-II via FFT (O(N log N))
		//
		// ALGORITHM:
		// 1. Create mirrored extension: y[n] = x[n] for n<N, y[n] = x[2N-1-n] for n>=N
		// 2. Compute FFT of length 2N
		// 3. Extract DCT coefficients: X[k] = Re(exp(-iπk/2N) * Y[k]) * (2/N)
		//
		// REQUIREMENTS: N must be a power of 2
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> ForwardII_Fast(const Vector<Real>& data) {
			int N = data.size();
			if (N == 0) {
				throw std::invalid_argument("DCT::ForwardII_Fast - empty input");
			}
			if (!FFT::IsPowerOfTwo(N)) {
				throw std::invalid_argument("DCT::ForwardII_Fast - size must be power of 2");
			}

			// Create extended sequence of length 2N with even symmetry
			int N2 = 2 * N;
			Vector<Complex> y(N2);
			for (int n = 0; n < N; ++n) {
				y[n] = Complex(data[n], 0.0);
			}
			for (int n = N; n < N2; ++n) {
				y[n] = Complex(data[N2 - 1 - n], 0.0);
			}

			// Compute FFT of extended sequence
			FFT::Transform(y, 1);

			// Extract DCT coefficients
			const Real pi = Constants::PI;
			Real scale = 1.0 / N;

			Vector<Real> result(N);
			for (int k = 0; k < N; ++k) {
				Real theta = -pi * k / N2;
				Complex twiddle(std::cos(theta), std::sin(theta));
				Complex product = twiddle * y[k];
				result[k] = scale * product.real();
			}

			return result;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// InverseII_Fast - Fast DCT-III (inverse DCT-II) via FFT (O(N log N))
		//
		// REQUIREMENTS: N must be a power of 2
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> InverseII_Fast(const Vector<Real>& coefficients) {
			int N = coefficients.size();
			if (N == 0) {
				throw std::invalid_argument("DCT::InverseII_Fast - empty input");
			}
			if (!FFT::IsPowerOfTwo(N)) {
				throw std::invalid_argument("DCT::InverseII_Fast - size must be power of 2");
			}

			const Real pi = Constants::PI;
			int N2 = 2 * N;

			Vector<Complex> Y(N2, Complex(0.0, 0.0));

			for (int k = 0; k < N; ++k) {
				Real coeff = coefficients[k];

				Real theta = pi * k / N2;
				Complex twiddle(std::cos(theta), std::sin(theta));
				Y[k] = coeff * twiddle;
			}

			// Set up conjugate symmetry for real output
			for (int k = 1; k < N; ++k) {
				Y[N2 - k] = std::conj(Y[k]);
			}

			// Inverse FFT (unnormalized)
			FFT::Transform(Y, -1);

			// Extract first N samples, divide by 2
			Vector<Real> result(N);
			for (int n = 0; n < N; ++n) {
				result[n] = Y[n].real() * 0.5;
			}

			return result;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// ForwardII_Reference - Reference O(N²) implementation
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> ForwardII_Reference(const Vector<Real>& data) {
			int N = data.size();
			if (N == 0) {
				throw std::invalid_argument("DCT::ForwardII_Reference - empty input");
			}

			const Real pi = Constants::PI;
			Real fac = 2.0 / N;

			Vector<Real> result(N);
			for (int k = 0; k < N; k++) {
				Real sum = 0.0;
				for (int n = 0; n < N; n++) {
					sum += data[n] * std::cos(pi * k * (n + 0.5) / N);
				}
				result[k] = fac * sum;
			}

			return result;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// InverseII_Reference - Reference O(N²) implementation
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> InverseII_Reference(const Vector<Real>& coefficients) {
			int N = coefficients.size();
			if (N == 0) {
				throw std::invalid_argument("DCT::InverseII_Reference - empty input");
			}

			const Real pi = Constants::PI;

			Vector<Real> result(N);
			for (int n = 0; n < N; n++) {
				// First term: DC component (halved for orthogonality)
				Real sum = coefficients[0] / 2.0;

				// Remaining terms
				for (int k = 1; k < N; k++) {
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
		// PROPERTIES:
		// - Uses orthonormal basis (self-adjoint transform)
		// - Used for problems with Dirichlet boundary conditions
		// - Common in PDE solving (finite difference methods)
		//
		// PERFORMANCE: O(N²) reference implementation
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> ForwardDST(const Vector<Real>& data) {
			int N = data.size();
			if (N == 0) {
				throw std::invalid_argument("DCT::ForwardDST - empty input");
			}

			const Real pi = Constants::PI;
			Real norm = std::sqrt(2.0 / (N + 1)); // Orthonormal basis

			Vector<Real> result(N);
			for (int k = 0; k < N; k++) {
				Real sum = 0.0;
				for (int n = 0; n < N; n++) {
					sum += data[n] * std::sin(pi * (k + 1) * (n + 1) / (N + 1));
				}
				result[k] = norm * sum;
			}

			return result;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// InverseDST - DST Type-I inverse
		// DST-I with orthonormal basis is self-adjoint (its own inverse)
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> InverseDST(const Vector<Real>& coefficients) {
			return ForwardDST(coefficients);
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// HELPER: Verify DCT-II/III round-trip
		///////////////////////////////////////////////////////////////////////////////////////////
		static bool VerifyRoundTrip(const Vector<Real>& data, Real tolerance = 1e-12) {
			Vector<Real> coeffs = ForwardII(data);
			Vector<Real> recovered = InverseII(coeffs);

			if (recovered.size() != data.size())
				return false;

			for (int i = 0; i < data.size(); i++) {
				if (std::abs(recovered[i] - data[i]) > tolerance)
					return false;
			}
			return true;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// HELPER: Verify fast implementation matches reference
		///////////////////////////////////////////////////////////////////////////////////////////
		static bool VerifyFastMatchesReference(const Vector<Real>& data, Real tolerance = 1e-10) {
			int N = data.size();
			if (!FFT::IsPowerOfTwo(N)) {
				return true; // Not applicable for non-power-of-2
			}

			Vector<Real> fast = ForwardII_Fast(data);
			Vector<Real> ref = ForwardII_Reference(data);

			if (fast.size() != ref.size())
				return false;

			for (int i = 0; i < N; ++i) {
				if (std::abs(fast[i] - ref[i]) > tolerance)
					return false;
			}
			return true;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Detailed API variants - return FourierResult with AlgorithmStatus + timing
		///////////////////////////////////////////////////////////////////////////////////////////

		/// DCT-II forward with full diagnostics (auto-selects fast/reference)
		static FourierResult<Vector<Real>> ForwardIIDetailed(
			const Vector<Real>& data, const FourierConfig& config = {})
		{
			using ResultType = FourierResult<Vector<Real>>;
			return FourierDetail::ExecuteFourierDetailed<ResultType>(
				"DCT::ForwardII", config,
				[&](ResultType& result) {
					result.input_size = data.size();
					result.value = ForwardII(data);
					result.function_evaluations = FFT::IsPowerOfTwo(data.size()) ?
						data.size() * FFT::Log2(data.size()) : data.size() * data.size();
				});
		}

		/// DCT-II inverse (DCT-III) with full diagnostics
		static FourierResult<Vector<Real>> InverseIIDetailed(
			const Vector<Real>& coefficients, const FourierConfig& config = {})
		{
			using ResultType = FourierResult<Vector<Real>>;
			return FourierDetail::ExecuteFourierDetailed<ResultType>(
				"DCT::InverseII", config,
				[&](ResultType& result) {
					result.input_size = coefficients.size();
					result.value = InverseII(coefficients);
					result.function_evaluations = FFT::IsPowerOfTwo(coefficients.size()) ?
						coefficients.size() * FFT::Log2(coefficients.size()) : coefficients.size() * coefficients.size();
				});
		}

		/// DST-I forward with full diagnostics
		static FourierResult<Vector<Real>> ForwardDSTDetailed(
			const Vector<Real>& data, const FourierConfig& config = {})
		{
			using ResultType = FourierResult<Vector<Real>>;
			return FourierDetail::ExecuteFourierDetailed<ResultType>(
				"DCT::ForwardDST", config,
				[&](ResultType& result) {
					result.input_size = data.size();
					result.value = ForwardDST(data);
					result.function_evaluations = data.size() * data.size();
				});
		}

		/// DST-I inverse with full diagnostics
		static FourierResult<Vector<Real>> InverseDSTDetailed(
			const Vector<Real>& coefficients, const FourierConfig& config = {})
		{
			using ResultType = FourierResult<Vector<Real>>;
			return FourierDetail::ExecuteFourierDetailed<ResultType>(
				"DCT::InverseDST", config,
				[&](ResultType& result) {
					result.input_size = coefficients.size();
					result.value = InverseDST(coefficients);
					result.function_evaluations = coefficients.size() * coefficients.size();
				});
		}
	};

} // namespace MML::Fourier
#endif // MML_FOURIER_H
