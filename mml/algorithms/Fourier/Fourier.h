///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Fourier/Fourier.h                                                   ///
///  Description: Fourier transforms for MML core algorithms layer                    ///
///               DFT  - Naive O(n²) reference implementation                        ///
///               FFT  - Radix-2 and Bluestein arbitrary-size O(n log n)             ///
///               DCT  - Discrete Cosine Transform (DCT-II/III, DST-I)               ///
///               Fast O(N log N) DCT via FFT                                        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                             ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FOURIER_H
#define MML_FOURIER_H

#include <mml/MMLBase.h>
#include <mml/MMLExceptions.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/AlgorithmTypes.h>

#include <algorithm>
#include <limits>
#include <cstdint>
#include <new>
#include <type_traits>
#include <stdexcept>
#include <string>

namespace MML::Fourier 
{
	enum class TransformNormalization {
		Legacy,
		None,
		Forward,
		Inverse,
		Orthonormal
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// FourierConfig - Configuration for Fourier transform detailed APIs
	///////////////////////////////////////////////////////////////////////////////////////////
	struct FourierConfig : public EvaluationConfigBase {
		// Inherits: estimate_error, check_finite, exception_policy
		TransformNormalization normalization = TransformNormalization::Legacy;
		bool allow_padding = false;
		int target_size = 0;
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

		/// Output size produced by the transform
		int output_size = 0;

		/// Normalization convention used for this result
		TransformNormalization normalization = TransformNormalization::Legacy;

		/// Whether padding was applied before computing the result
		bool padding_applied = false;

		/// Direct, Radix2, Bluestein, FFT/Radix2, or FFT/Bluestein.
		std::string kernel_name;

		/// Largest direct or FFT kernel length (Bluestein's convolution length).
		int work_size = 0;

		/// Summation terms, FFT butterflies and complex products; excludes scaling/copies.
		std::uint64_t work_units = 0;
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// FourierDetail - Internal helpers for Detailed API execution
	///////////////////////////////////////////////////////////////////////////////////////////
	namespace FourierDetail
	{
		inline void ValidateNormalization(TransformNormalization normalization)
		{
			switch (normalization) {
			case TransformNormalization::Legacy:
			case TransformNormalization::None:
			case TransformNormalization::Forward:
			case TransformNormalization::Inverse:
			case TransformNormalization::Orthonormal:
				return;
			}
			throw FourierError("Fourier transform - unsupported normalization");
		}

		template<typename T>
		void ValidateFinite(const Vector<T>& values, bool input)
		{
			for (int i = 0; i < values.size(); ++i) {
				bool finite;
				if constexpr (std::is_same_v<T, Complex>)
					finite = std::isfinite(values[i].real()) && std::isfinite(values[i].imag());
				else
					finite = std::isfinite(values[i]);
				if (!finite) {
					if (input)
						throw NumericInputError("Fourier transform - nonfinite input at index " + std::to_string(i));
					throw NumericalMethodError("Fourier transform - nonfinite output at index " + std::to_string(i));
				}
			}
		}

		inline int NextPowerOfTwo(int size)
		{
			int power = 1;
			while (power < size) {
				if (power > std::numeric_limits<int>::max() / 2)
					throw FourierError("FFT::NextPowerOfTwo - size exceeds supported range");
				power *= 2;
			}
			return power;
		}

		inline int DoubledSize(int size)
		{
			if (size <= 0 || size > std::numeric_limits<int>::max() / 2)
				throw ResourceLimitError("Fourier transform - doubled work size exceeds supported range");
			return 2 * size;
		}

		inline int BluesteinSize(int size)
		{
			const int required = DoubledSize(size) - 1;
			if (required > std::numeric_limits<int>::max() / 2 + 1)
				throw ResourceLimitError("FFT - Bluestein convolution size exceeds supported range");
			return NextPowerOfTwo(required);
		}

		struct Work {
			const char* kernel;
			int size;
			std::uint64_t units;
		};

		inline Work DirectWork(int size)
		{
			if (size <= 0)
				throw FourierError("Fourier transform - direct work size must be positive");
			return {"Direct", size, std::uint64_t(size) * std::uint64_t(size)};
		}

		inline Work FFTWork(int size)
		{
			if (size <= 0)
				throw FourierError("FFT - empty input");
			const bool radix2 = (size & (size - 1)) == 0;
			const int workSize = radix2 ? size : BluesteinSize(size);
			int stages = 0;
			for (int n = workSize; n > 1; n /= 2) ++stages;
			const auto butterflies = std::uint64_t(workSize / 2) * std::uint64_t(stages);
			return {radix2 ? "Radix2" : "Bluestein", workSize,
				radix2 ? butterflies : 3 * butterflies + std::uint64_t(workSize) + 2 * std::uint64_t(size)};
		}

		template<class Result>
		void SetWork(Result& result, const Work& work)
		{
			result.kernel_name = work.kernel;
			result.work_size = work.size;
			result.work_units = work.units;
		}

		/// Execute a Fourier Detailed operation with timing and exception handling.
		template<typename ResultType, typename T, typename ComputeFn>
		ResultType ExecuteFourierDetailed(const char* algorithm_name,
		                                 const Vector<T>& input,
		                                 const FourierConfig& config,
		                                 ComputeFn&& compute,
		                                 bool supportsPadding = false)
		{
			AlgorithmTimer timer;
			ResultType result = MakeEvaluationSuccessResult<ResultType>(algorithm_name);
			result.input_size = input.size();
			result.normalization = config.normalization;
			auto execute = [&]() {
				ValidateNormalization(config.normalization);
				if (input.size() == 0)
					throw FourierError(std::string(algorithm_name) + " - empty input");
				if (config.target_size < 0 || (!config.allow_padding && config.target_size != 0) ||
					(!supportsPadding && (config.allow_padding || config.target_size != 0)))
					throw FourierError(std::string(algorithm_name) + " - unsupported padding options");
				if (config.check_finite) ValidateFinite(input, true);
				compute(result);
				result.output_size = result.value.size();
				if (config.check_finite) ValidateFinite(result.value, false);
			};
			auto fail = [&](AlgorithmStatus status, const std::exception& error) {
				result.status = status;
				result.error_message = error.what();
			};
			if (config.exception_policy == EvaluationExceptionPolicy::Propagate) {
				execute();
			} else {
				try { execute(); }
				catch (const std::invalid_argument& ex) { fail(AlgorithmStatus::InvalidInput, ex); }
				catch (const NumericInputError& ex) { fail(AlgorithmStatus::InvalidInput, ex); }
				catch (const NumericalMethodError& ex) { fail(AlgorithmStatus::NumericalInstability, ex); }
				catch (const ResourceLimitError& ex) { fail(AlgorithmStatus::ResourceLimitExceeded, ex); }
				catch (const std::bad_alloc& ex) { fail(AlgorithmStatus::ResourceLimitExceeded, ex); }
				catch (const std::length_error& ex) { fail(AlgorithmStatus::ResourceLimitExceeded, ex); }
			}
			result.elapsed_time_ms = timer.elapsed_ms();
			return result;
		}

		inline Real NormalizationScale(int size, TransformNormalization normalization, bool forward_direction)
		{
			if (normalization == TransformNormalization::Legacy)
				normalization = forward_direction ? TransformNormalization::None : TransformNormalization::Inverse;

			switch (normalization) {
			case TransformNormalization::None:
					return 1.0;
			case TransformNormalization::Forward:
				return forward_direction ? MML::Real(1.0) / size : MML::Real(1.0);
			case TransformNormalization::Inverse:
				return forward_direction ? MML::Real(1.0) : MML::Real(1.0) / size;
			case TransformNormalization::Orthonormal:
				return MML::Real(1.0) / std::sqrt(static_cast<Real>(size));
			case TransformNormalization::Legacy:
				return 1.0;
			}

			throw FourierError("Fourier transform - unsupported normalization");
		}

		template<typename T>
		void ApplyScale(Vector<T>& values, Real scale)
		{
			if (scale == 1.0)
				return;

			for (int i = 0; i < values.size(); i++)
				values[i] *= scale;
		}

		inline int CeilLog2(int size)
		{
			int exponent = 0;
			for (int remaining = size - 1; remaining > 0; remaining /= 2)
				++exponent;
			return exponent;
		}

		inline Real ScalePowerOfTwo(Real value, int exponent)
		{
			if (exponent == 0) return value;
			return std::scalbn(value, exponent);
		}

		inline Complex ScalePowerOfTwo(Complex value, int exponent)
		{
			if (exponent == 0) return value;
			return {std::scalbn(value.real(), exponent), std::scalbn(value.imag(), exponent)};
		}

		template<typename T>
		int TransformInputShift(const Vector<T>& values, int growthExponent)
		{
			Real largest = 0;
			for (const auto& value : values) {
				if constexpr (std::is_same_v<T, Complex>) {
					// Preserve the existing nonfinite policy, including unchecked Detailed calls.
					if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
						return 0;
					largest = std::max({largest, std::abs(value.real()), std::abs(value.imag())});
				} else {
					if (!std::isfinite(value))
						return 0;
					largest = std::max(largest, std::abs(value));
				}
			}
			int exponent = 0;
			std::frexp(largest, &exponent);
			return std::max(0, exponent + growthExponent - std::numeric_limits<Real>::max_exponent);
		}

		template<typename T>
		void RestoreTransformScale(Vector<T>& values, Real scale, int inputShift)
		{
			if (inputShift == 0) {
				ApplyScale(values, scale);
				return;
			}
			int scaleExponent = 0;
			const Real scaleFraction = std::frexp(scale, &scaleExponent);
			auto restore = [&](Real value) {
				int exponent = 0;
				const Real fraction = std::frexp(value, &exponent);
				// Neither an overflowing reciprocal input scale nor a premature tiny
				// product is needed to combine normalization with exponent restoration.
				return std::scalbn(fraction * scaleFraction, exponent + scaleExponent + inputShift);
			};
			for (auto& value : values) {
				if constexpr (std::is_same_v<T, Complex>)
					value = Complex(restore(value.real()), restore(value.imag()));
				else
					value = restore(value);
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
				throw FourierError(functionName + " - empty input");
			}
		}

		/// @brief Validate that a size is positive (> 0)
		inline void ValidatePositiveSize(int size, const std::string& paramName, const std::string& functionName) {
			if (size <= 0) {
				throw FourierError(functionName + " - " + paramName + " must be positive");
			}
		}

		/// @brief Validate that a size meets a minimum requirement
		inline void ValidateMinSize(int size, int minSize, const std::string& functionName) {
			if (size < minSize) {
				throw FourierError(functionName + " - size must be at least " + std::to_string(minSize));
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Power-of-2 validations (critical for FFT)
		///////////////////////////////////////////////////////////////////////////////////////////

		inline bool IsPowerOfTwo(int n) { return (n > 0) && ((n & (n - 1)) == 0); }

		inline void ValidatePowerOfTwo(int size, const std::string& functionName) {
			if (!IsPowerOfTwo(size)) {
				throw FourierError(functionName + " - size must be power of 2");
			}
		}

		inline void ValidateFFTSize(int size, const std::string& functionName) {
			ValidateNonEmpty(size, functionName);
			ValidatePowerOfTwo(size, functionName);
		}

		inline void ValidateRealFFTSize(int size, int minSize, const std::string& functionName) {
			if (size < minSize || !IsPowerOfTwo(size)) {
				throw FourierError(functionName + " - requires power-of-2 size >= " + std::to_string(minSize));
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Window validations
		///////////////////////////////////////////////////////////////////////////////////////////

		inline void ValidateWindowSize(int n) {
			if (n <= 0) {
				throw FourierError("Window size must be positive");
			}
		}

		inline void ValidateWindowMatch(int windowSize, int dataSize, const std::string& functionName) {
			if (windowSize != dataSize) {
				throw FourierError(functionName + " - window size mismatch");
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Parameter validations
		///////////////////////////////////////////////////////////////////////////////////////////

		template<typename T>
		inline void ValidatePositive(T value, const std::string& paramName, const std::string& functionName) {
			if (!std::isfinite(value) || value <= T(0)) {
				throw FourierError(functionName + " - " + paramName + " must be finite and positive");
			}
		}

		template<typename T>
		inline void ValidateNonNegative(T value, const std::string& paramName, const std::string& functionName) {
			if (!std::isfinite(value) || value < T(0)) {
				throw FourierError(functionName + " - " + paramName + " must be finite and non-negative");
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Index validations
		///////////////////////////////////////////////////////////////////////////////////////////

		inline void ValidateIndex(int index, int maxIndex, const std::string& functionName) {
			if (index < 0 || index >= maxIndex) {
				throw FourierError(functionName + " - invalid index");
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
			return Forward(data, TransformNormalization::Legacy);
		}

		// Forward DFT with explicit normalization
		static Vector<Complex> Forward(const Vector<Complex>& data, TransformNormalization normalization) {
			int N = data.size();
			FourierValidation::ValidateNonEmpty(N, "DFT::Forward");
			const Real scale = FourierDetail::NormalizationScale(N, normalization, true);
			const int inputShift = FourierDetail::TransformInputShift(data, FourierDetail::CeilLog2(N) + 2);

			Vector<Complex> result(N);
			const Real two_pi = 2.0 * Constants::PI;

			for (int k = 0; k < N; k++) {
				Complex sum(0.0, 0.0);
				for (int n = 0; n < N; n++) {
					Real angle = -two_pi * k * n / N;
					Complex twiddle(std::cos(angle), std::sin(angle));
					sum += FourierDetail::ScalePowerOfTwo(data[n], -inputShift) * twiddle;
				}
				result[k] = sum;
			}

			FourierDetail::RestoreTransformScale(result, scale, inputShift);

			return result;
		}

		// Forward DFT for real-valued input (converts to complex internally)
		static Vector<Complex> Forward(const Vector<Real>& data) {
			return Forward(data, TransformNormalization::Legacy);
		}

		// Forward DFT for real-valued input with explicit normalization
		static Vector<Complex> Forward(const Vector<Real>& data, TransformNormalization normalization) {
			int N = data.size();
			FourierValidation::ValidateNonEmpty(N, "DFT::Forward");

			Vector<Complex> complex_data(N);
			for (int i = 0; i < N; i++) {
				complex_data[i] = Complex(data[i], 0.0);
			}

			return Forward(complex_data, normalization);
		}

		// Inverse DFT: frequency domain -> time domain
		static Vector<Complex> Inverse(const Vector<Complex>& spectrum) {
			return Inverse(spectrum, TransformNormalization::Legacy);
		}

		// Inverse DFT with explicit normalization
		static Vector<Complex> Inverse(const Vector<Complex>& spectrum, TransformNormalization normalization) {
			int N = spectrum.size();
			FourierValidation::ValidateNonEmpty(N, "DFT::Inverse");
			const Real scale = FourierDetail::NormalizationScale(N, normalization, false);
			const int inputShift = FourierDetail::TransformInputShift(spectrum, FourierDetail::CeilLog2(N) + 2);

			Vector<Complex> result(N);
			const Real two_pi = 2.0 * Constants::PI;

			for (int n = 0; n < N; n++) {
				Complex sum(0.0, 0.0);
				for (int k = 0; k < N; k++) {
					Real angle = two_pi * k * n / N; // Positive sign for inverse
					Complex twiddle(std::cos(angle), std::sin(angle));
					sum += FourierDetail::ScalePowerOfTwo(spectrum[k], -inputShift) * twiddle;
				}
				result[n] = sum;
			}

			FourierDetail::RestoreTransformScale(result, scale, inputShift);

			return result;
		}

		// Extract real part from inverse DFT (for real-valued signals)
		static Vector<Real> InverseReal(const Vector<Complex>& spectrum) {
			return InverseReal(spectrum, TransformNormalization::Legacy);
		}

		// Extract real part from inverse DFT with explicit normalization
		static Vector<Real> InverseReal(const Vector<Complex>& spectrum, TransformNormalization normalization) {
			Vector<Complex> complex_result = Inverse(spectrum, normalization);
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
				"DFT::Forward", data, config,
				[&](ResultType& result) {
					result.input_size = data.size();
					result.value = Forward(data, config.normalization);
					result.output_size = result.value.size();
					result.normalization = config.normalization;
					FourierDetail::SetWork(result, FourierDetail::DirectWork(data.size()));
				});
		}

		/// Forward DFT with full diagnostics (real input)
		static FourierResult<Vector<Complex>> ForwardDetailed(
			const Vector<Real>& data, const FourierConfig& config = {})
		{
			using ResultType = FourierResult<Vector<Complex>>;
			return FourierDetail::ExecuteFourierDetailed<ResultType>(
				"DFT::Forward", data, config,
				[&](ResultType& result) {
					result.input_size = data.size();
					result.value = Forward(data, config.normalization);
					result.output_size = result.value.size();
					result.normalization = config.normalization;
					FourierDetail::SetWork(result, FourierDetail::DirectWork(data.size()));
				});
		}

		/// Inverse DFT with full diagnostics
		static FourierResult<Vector<Complex>> InverseDetailed(
			const Vector<Complex>& spectrum, const FourierConfig& config = {})
		{
			using ResultType = FourierResult<Vector<Complex>>;
			return FourierDetail::ExecuteFourierDetailed<ResultType>(
				"DFT::Inverse", spectrum, config,
				[&](ResultType& result) {
					result.input_size = spectrum.size();
					result.value = Inverse(spectrum, config.normalization);
					result.output_size = result.value.size();
					result.normalization = config.normalization;
					FourierDetail::SetWork(result, FourierDetail::DirectWork(spectrum.size()));
				});
		}

		/// Inverse DFT returning real part, with full diagnostics
		static FourierResult<Vector<Real>> InverseRealDetailed(
			const Vector<Complex>& spectrum, const FourierConfig& config = {})
		{
			using ResultType = FourierResult<Vector<Real>>;
			return FourierDetail::ExecuteFourierDetailed<ResultType>(
				"DFT::InverseReal", spectrum, config,
				[&](ResultType& result) {
					result.input_size = spectrum.size();
					result.value = InverseReal(spectrum, config.normalization);
					result.output_size = result.value.size();
					result.normalization = config.normalization;
					FourierDetail::SetWork(result, FourierDetail::DirectWork(spectrum.size()));
				});
		}
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// FFT - Fast Fourier Transform
	//
	// Cooley-Tukey radix-2 FFT with Bluestein dispatch - O(n log n) complexity
	//
	// Forward transform:  X[k] = Σ_{n=0}^{N-1} x[n] * exp(-2πi*k*n/N)
	// Inverse transform:  x[n] = (1/N) * Σ_{k=0}^{N-1} X[k] * exp(2πi*k*n/N)
	//
	// Algorithm: radix-2 decimation-in-time for powers of two, Bluestein chirp-z otherwise
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
			TransformScaled(data, isign, Real(1));
		}

	private:
		static void TransformScaled(Vector<Complex>& data, int isign, Real scale) {
			const int n = data.size();
			FourierValidation::ValidateNonEmpty(n, "FFT::Transform");
			if (isign != 1 && isign != -1)
				throw FourierError("FFT::Transform - sign must be +1 or -1");
			const bool radix2 = IsPowerOfTwo(n);
			const int workSize = radix2 ? n : FourierDetail::BluesteinSize(n);
			// Radix-2 grows by at most O(N). Bluestein additionally needs headroom
			// for two forward FFTs, their complex products and the inverse FFT.
			const int stages = FourierDetail::CeilLog2(workSize);
			const int growthExponent = radix2 ? stages + 2 : 3 * stages + 6;
			const int inputShift = FourierDetail::TransformInputShift(data, growthExponent);
			if (inputShift != 0)
				for (auto& value : data)
					value = FourierDetail::ScalePowerOfTwo(value, -inputShift);
			if (radix2)
				TransformRadix2(data, isign);
			else
				TransformBluestein(data, isign);
			FourierDetail::RestoreTransformScale(data, scale, inputShift);
		}

		// Existing iterative Cooley-Tukey kernel. Input size must be a power of two.
		static void TransformRadix2(Vector<Complex>& data, int isign) {
			int n = data.size();

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

		// Bluestein rewrites an N-point DFT as a linear convolution of chirp sequences.
		static void TransformBluestein(Vector<Complex>& data, int isign) {
			const int n = data.size();
			const int convolutionSize = FourierDetail::BluesteinSize(n);
			Vector<Complex> chirpedData(convolutionSize, Complex(0.0, 0.0));
			Vector<Complex> chirpKernel(convolutionSize, Complex(0.0, 0.0));

			for (int index = 0; index < n; index++) {
				const Real angle = isign * Constants::PI *
					(static_cast<Real>(index) * static_cast<Real>(index) / n);
				const Complex inputChirp(std::cos(angle), -std::sin(angle));
				const Complex convolutionChirp(std::cos(angle), std::sin(angle));
				chirpedData[index] = data[index] * inputChirp;
				chirpKernel[index] = convolutionChirp;
				if (index != 0)
					chirpKernel[convolutionSize - index] = convolutionChirp;
			}

			TransformRadix2(chirpedData, 1);
			TransformRadix2(chirpKernel, 1);
			for (int index = 0; index < convolutionSize; index++)
				chirpedData[index] *= chirpKernel[index];
			TransformRadix2(chirpedData, -1);

			const Real inverseConvolutionSize = MML::Real(1.0) / convolutionSize;
			for (int index = 0; index < n; index++) {
				const Real angle = isign * Constants::PI *
					(static_cast<Real>(index) * static_cast<Real>(index) / n);
				const Complex outputChirp(std::cos(angle), -std::sin(angle));
				data[index] = chirpedData[index] * outputChirp * inverseConvolutionSize;
			}
		}

	public:

		// Forward FFT: time domain -> frequency domain
		static Vector<Complex> Forward(const Vector<Complex>& data) {
			return Forward(data, TransformNormalization::Legacy);
		}

		// Forward FFT with explicit normalization
		static Vector<Complex> Forward(const Vector<Complex>& data, TransformNormalization normalization) {
			FourierValidation::ValidateNonEmpty(data.size(), "FFT::Forward");
			const Real scale = FourierDetail::NormalizationScale(data.size(), normalization, true);
			Vector<Complex> result = data; // Copy input
			TransformScaled(result, 1, scale);
			return result;
		}

		// Inverse FFT: frequency domain -> time domain
		static Vector<Complex> Inverse(const Vector<Complex>& spectrum) {
			return Inverse(spectrum, TransformNormalization::Legacy);
		}

		// Inverse FFT with explicit normalization
		static Vector<Complex> Inverse(const Vector<Complex>& spectrum, TransformNormalization normalization) {
			FourierValidation::ValidateNonEmpty(spectrum.size(), "FFT::Inverse");
			const Real scale = FourierDetail::NormalizationScale(spectrum.size(), normalization, false);
			Vector<Complex> result = spectrum; // Copy input
			TransformScaled(result, -1, scale);

			return result;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// Utility functions
		///////////////////////////////////////////////////////////////////////////////////////////

		// Check if n is a power of 2
		static bool IsPowerOfTwo(int n) { return (n > 0) && ((n & (n - 1)) == 0); }

		// Find next power of 2 >= n
		static int NextPowerOfTwo(int n) {
			return FourierDetail::NextPowerOfTwo(n);
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
			return TransformDetailed(data, config, true);
		}

		/// Inverse FFT with full diagnostics
		static FourierResult<Vector<Complex>> InverseDetailed(
			const Vector<Complex>& spectrum, const FourierConfig& config = {})
		{
			return TransformDetailed(spectrum, config, false);
		}

	private:
		static FourierResult<Vector<Complex>> TransformDetailed(
			const Vector<Complex>& input, const FourierConfig& config, bool forward)
		{
			using ResultType = FourierResult<Vector<Complex>>;
			return FourierDetail::ExecuteFourierDetailed<ResultType>(
				forward ? "FFT::Forward" : "FFT::Inverse", input, config,
				[&](ResultType& result) {
					const int size = !config.allow_padding ? input.size() :
						(config.target_size == 0 ? NextPowerOfTwo(input.size()) : config.target_size);
					if (size < input.size())
						throw FourierError("FFT - padding target cannot truncate input");
					const auto work = FourierDetail::FFTWork(size);
					Vector<Complex> padded;
					const Vector<Complex>* effective = &input;
					if (size != input.size()) {
						padded = Vector<Complex>(size, Complex(0, 0));
						for (int i = 0; i < input.size(); ++i) padded[i] = input[i];
						effective = &padded;
						result.padding_applied = true;
					}
					result.value = forward ? Forward(*effective, config.normalization) : Inverse(*effective, config.normalization);
					FourierDetail::SetWork(result, work);
				}, true);
		}

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
	// - Fast DCT via FFT: O(N log N)
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

		static FourierDetail::Work DescribeWork(int size)
		{
			if (size < FAST_DCT_THRESHOLD)
				return FourierDetail::DirectWork(size);
			auto work = FourierDetail::FFTWork(FourierDetail::DoubledSize(size));
			work.kernel = FFT::IsPowerOfTwo(size) ? "FFT/Radix2" : "FFT/Bluestein";
			work.units += std::uint64_t(size);
			return work;
		}

		static Real SineScale(int size, TransformNormalization normalization, bool forward)
		{
			const Real rawScale = std::sqrt(2 * (Real(size) + 1));
			switch (normalization) {
			case TransformNormalization::Legacy:
			case TransformNormalization::Orthonormal:
				return 1;
			case TransformNormalization::None:
				return rawScale;
			case TransformNormalization::Forward:
				return forward ? 1 / rawScale : rawScale;
			case TransformNormalization::Inverse:
				return forward ? rawScale : 1 / rawScale;
			}
			throw FourierError("DST - unsupported normalization");
		}

	public:
		///////////////////////////////////////////////////////////////////////////////////////////
		// ForwardII - DCT Type-II (standard "DCT") - AUTO-SELECTS FASTEST ALGORITHM
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> ForwardII(const Vector<Real>& data) {
			int N = data.size();
			if (N == 0) {
				throw FourierError("DCT::ForwardII - empty input");
			}

			// Small sizes use the reference implementation to avoid FFT setup overhead.
			if (N < FAST_DCT_THRESHOLD) {
				return ForwardII_Reference(data);
			}

			// Fast DCT-II via arbitrary-size FFT
			return ForwardII_Fast(data);
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// InverseII - DCT Type-III (inverse of DCT-II) - AUTO-SELECTS FASTEST ALGORITHM
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> InverseII(const Vector<Real>& coefficients) {
			int N = coefficients.size();
			if (N == 0) {
				throw FourierError("DCT::InverseII - empty input");
			}

			// Small sizes use the reference implementation to avoid FFT setup overhead.
			if (N < FAST_DCT_THRESHOLD) {
				return InverseII_Reference(coefficients);
			}

			// Fast DCT-III via FFT
			return InverseII_Fast(coefficients);
		}

		static Vector<Real> ForwardII(const Vector<Real>& data, TransformNormalization normalization)
		{
			FourierDetail::ValidateNormalization(normalization);
			if ((normalization == TransformNormalization::Forward || normalization == TransformNormalization::Orthonormal) &&
				data.size() == 1)
				return data;
			auto result = ForwardII(data);
			switch (normalization) {
			case TransformNormalization::Legacy:
				break;
			case TransformNormalization::None:
			case TransformNormalization::Inverse:
				FourierDetail::ApplyScale(result, Real(data.size()));
				break;
			case TransformNormalization::Orthonormal:
				FourierDetail::ApplyScale(result, std::sqrt(Real(data.size()) / 2));
				result[0] /= std::sqrt(Real(2));
				break;
			case TransformNormalization::Forward:
				FourierDetail::ApplyScale(result, Real(0.5));
				break;
			}
			return result;
		}

		static Vector<Real> InverseII(const Vector<Real>& coefficients, TransformNormalization normalization)
		{
			FourierDetail::ValidateNormalization(normalization);
			if (normalization == TransformNormalization::Orthonormal) {
				FourierValidation::ValidateNonEmpty(coefficients.size(), "DCT::InverseII");
				auto scaled = coefficients;
				scaled[0] = 0;
				FourierDetail::ApplyScale(scaled, std::sqrt(Real(2) / coefficients.size()));
				auto result = InverseII(scaled);
				const Real dc = coefficients[0] / std::sqrt(Real(coefficients.size()));
				for (int i = 0; i < result.size(); ++i) result[i] += dc;
				return result;
			}
			if (normalization == TransformNormalization::Inverse) {
				FourierValidation::ValidateNonEmpty(coefficients.size(), "DCT::InverseII");
				auto scaled = coefficients;
				FourierDetail::ApplyScale(scaled, Real(1) / coefficients.size());
				return InverseII(scaled);
			}
			auto result = InverseII(coefficients);
			if (normalization != TransformNormalization::Legacy)
				FourierDetail::ApplyScale(result, Real(2));
			return result;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// ForwardII_Fast - Fast DCT-II via FFT (O(N log N))
		//
		// ALGORITHM:
		// 1. Create mirrored extension: y[n] = x[n] for n<N, y[n] = x[2N-1-n] for n>=N
		// 2. Compute FFT of length 2N
		// 3. Extract legacy coefficients: X[k] = Re(exp(-iπk/2N) * Y[k]) / N
		//
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> ForwardII_Fast(const Vector<Real>& data) {
			int N = data.size();
			if (N == 0) {
				throw FourierError("DCT::ForwardII_Fast - empty input");
			}

			// Create extended sequence of length 2N with even symmetry
			const int N2 = FourierDetail::DoubledSize(N);
			FourierDetail::FFTWork(N2);
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
			Real scale = MML::Real(1.0) / N;

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
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> InverseII_Fast(const Vector<Real>& coefficients) {
			int N = coefficients.size();
			if (N == 0) {
				throw FourierError("DCT::InverseII_Fast - empty input");
			}

			const Real pi = Constants::PI;
			const int N2 = FourierDetail::DoubledSize(N);
			FourierDetail::FFTWork(N2);

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
				result[n] = Y[n].real() * MML::Real(0.5);
			}

			return result;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// ForwardII_Reference - Reference O(N²) implementation
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> ForwardII_Reference(const Vector<Real>& data) {
			int N = data.size();
			if (N == 0) {
				throw FourierError("DCT::ForwardII_Reference - empty input");
			}

			using AccumulationReal = std::common_type_t<Real, double>;
			const AccumulationReal pi = std::acos(AccumulationReal(-1));
			const AccumulationReal fac = AccumulationReal(2) / N;

			Vector<Real> result(N);
			for (int k = 0; k < N; k++) {
				AccumulationReal sum = 0;
				for (int n = 0; n < N; n++) {
					sum += data[n] * std::cos(pi * k * (n + AccumulationReal(0.5)) / N);
				}
				result[k] = Real(fac * sum);
			}

			return result;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// InverseII_Reference - Reference O(N²) implementation
		///////////////////////////////////////////////////////////////////////////////////////////
		static Vector<Real> InverseII_Reference(const Vector<Real>& coefficients) {
			int N = coefficients.size();
			if (N == 0) {
				throw FourierError("DCT::InverseII_Reference - empty input");
			}

			using AccumulationReal = std::common_type_t<Real, double>;
			const AccumulationReal pi = std::acos(AccumulationReal(-1));

			Vector<Real> result(N);
			for (int n = 0; n < N; n++) {
				// First term: DC component (halved for orthogonality)
				AccumulationReal sum = coefficients[0] / AccumulationReal(2);

				// Remaining terms
				for (int k = 1; k < N; k++) {
					sum += coefficients[k] * std::cos(pi * k * (n + AccumulationReal(0.5)) / N);
				}
				result[n] = Real(sum);
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
				throw FourierError("DCT::ForwardDST - empty input");
			}

			const Real pi = Constants::PI;
			const Real denominator = Real(N) + 1;
			Real norm = std::sqrt(MML::Real(2.0) / denominator); // Orthonormal basis

			Vector<Real> result(N);
			for (int k = 0; k < N; k++) {
				Real sum = 0.0;
				for (int n = 0; n < N; n++) {
					sum += data[n] * std::sin(pi * (k + 1) * (n + 1) / denominator);
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

		static Vector<Real> ForwardDST(const Vector<Real>& data, TransformNormalization normalization)
		{
			FourierDetail::ValidateNormalization(normalization);
			auto result = ForwardDST(data);
			FourierDetail::ApplyScale(result, SineScale(data.size(), normalization, true));
			return result;
		}

		static Vector<Real> InverseDST(const Vector<Real>& coefficients, TransformNormalization normalization)
		{
			FourierDetail::ValidateNormalization(normalization);
			auto result = InverseDST(coefficients);
			FourierDetail::ApplyScale(result, SineScale(coefficients.size(), normalization, false));
			return result;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// HELPER: Verify DCT-II/III round-trip
		///////////////////////////////////////////////////////////////////////////////////////////
		static bool VerifyRoundTrip(const Vector<Real>& data, Real tolerance = Precision::ZeroTolerance) {
			FourierValidation::ValidateNonEmpty(data.size(), "DCT::VerifyRoundTrip");
			FourierValidation::ValidateNonNegative(tolerance, "tolerance", "DCT::VerifyRoundTrip");
			FourierDetail::ValidateFinite(data, true);
			Vector<Real> coeffs = ForwardII(data);
			for (int i = 0; i < coeffs.size(); ++i)
				if (!std::isfinite(coeffs[i]))
					return false;
			Vector<Real> recovered = InverseII(coeffs);

			if (recovered.size() != data.size())
				return false;

			for (int i = 0; i < data.size(); i++) {
				if (!std::isfinite(recovered[i]) || !(std::abs(recovered[i] - data[i]) <= tolerance))
					return false;
			}
			return true;
		}

		///////////////////////////////////////////////////////////////////////////////////////////
		// HELPER: Verify fast implementation matches reference
		///////////////////////////////////////////////////////////////////////////////////////////
		static bool VerifyFastMatchesReference(const Vector<Real>& data, Real tolerance = Precision::EqualityTolerance) {
			FourierValidation::ValidateNonEmpty(data.size(), "DCT::VerifyFastMatchesReference");
			FourierValidation::ValidateNonNegative(tolerance, "tolerance", "DCT::VerifyFastMatchesReference");
			FourierDetail::ValidateFinite(data, true);
			int N = data.size();

			Vector<Real> fast = ForwardII_Fast(data);
			Vector<Real> ref = ForwardII_Reference(data);

			if (fast.size() != ref.size())
				return false;

			for (int i = 0; i < N; ++i) {
				if (!std::isfinite(fast[i]) || !std::isfinite(ref[i]) || !(std::abs(fast[i] - ref[i]) <= tolerance))
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
				"DCT::ForwardII", data, config,
				[&](ResultType& result) {
					result.input_size = data.size();
					const auto work = DescribeWork(data.size());
					result.value = ForwardII(data, config.normalization);
					FourierDetail::SetWork(result, work);
				});
		}

		/// DCT-II inverse (DCT-III) with full diagnostics
		static FourierResult<Vector<Real>> InverseIIDetailed(
			const Vector<Real>& coefficients, const FourierConfig& config = {})
		{
			using ResultType = FourierResult<Vector<Real>>;
			return FourierDetail::ExecuteFourierDetailed<ResultType>(
				"DCT::InverseII", coefficients, config,
				[&](ResultType& result) {
					result.input_size = coefficients.size();
					const auto work = DescribeWork(coefficients.size());
					result.value = InverseII(coefficients, config.normalization);
					FourierDetail::SetWork(result, work);
				});
		}

		/// DST-I forward with full diagnostics
		static FourierResult<Vector<Real>> ForwardDSTDetailed(
			const Vector<Real>& data, const FourierConfig& config = {})
		{
			using ResultType = FourierResult<Vector<Real>>;
			return FourierDetail::ExecuteFourierDetailed<ResultType>(
				"DCT::ForwardDST", data, config,
				[&](ResultType& result) {
					result.input_size = data.size();
					result.value = ForwardDST(data, config.normalization);
					FourierDetail::SetWork(result, FourierDetail::DirectWork(data.size()));
				});
		}

		/// DST-I inverse with full diagnostics
		static FourierResult<Vector<Real>> InverseDSTDetailed(
			const Vector<Real>& coefficients, const FourierConfig& config = {})
		{
			using ResultType = FourierResult<Vector<Real>>;
			return FourierDetail::ExecuteFourierDetailed<ResultType>(
				"DCT::InverseDST", coefficients, config,
				[&](ResultType& result) {
					result.input_size = coefficients.size();
					result.value = InverseDST(coefficients, config.normalization);
					FourierDetail::SetWork(result, FourierDetail::DirectWork(coefficients.size()));
				});
		}
	};

} // namespace MML::Fourier
#endif // MML_FOURIER_H
