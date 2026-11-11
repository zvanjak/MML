#include "../../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/Fourier/Fourier.h>
#include <mml/algorithms/Fourier/FourierRealFFT.h>
#endif
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <limits>

using namespace MML;
using namespace MML::Fourier;

namespace {
	using ComplexOperation = FourierResult<Vector<Complex>> (*)(const Vector<Complex>&, const FourierConfig&);
	using RealOperation = FourierResult<Vector<Real>> (*)(const Vector<Real>&, const FourierConfig&);
	const std::array<ComplexOperation, 4> complexOperations{
		DFT::ForwardDetailed, DFT::InverseDetailed, FFT::ForwardDetailed, FFT::InverseDetailed};
	const std::array<RealOperation, 4> realOperations{
		DCT::ForwardIIDetailed, DCT::InverseIIDetailed, DCT::ForwardDSTDetailed, DCT::InverseDSTDetailed};
	const std::array modes{TransformNormalization::Legacy, TransformNormalization::None,
		TransformNormalization::Forward, TransformNormalization::Inverse, TransformNormalization::Orthonormal};

	FourierConfig ConvertFailures() {
		FourierConfig config;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		return config;
	}

	Real SumRoundoffBound(int n, Real inputL1) {
		// O(n) phase argument error plus summation/butterfly rounding. The
		// independent trigonometric oracle below is evaluated in long double.
		return std::max(Precision::EqualityTolerance,
			Real(8) * n * std::numeric_limits<Real>::epsilon() * inputL1);
	}
}

TEST_CASE("Fourier Detailed finite policy covers every input component and transform", "[Fourier][FourierDetailedContract]") {
	auto config = ConvertFailures();
	for (const Real invalid : {std::numeric_limits<Real>::quiet_NaN(), Constants::PosInf, Constants::NegInf}) {
		for (auto operation : complexOperations) {
			for (const auto bad : {Complex(invalid, 0), Complex(0, invalid)}) {
				const auto result = operation(Vector<Complex>{Complex(1, 0), bad}, config);
				CHECK(result.status == AlgorithmStatus::InvalidInput);
				CHECK_FALSE(result.IsSuccess());
				CHECK_FALSE(result.error_message.empty());
				CHECK(result.input_size == 2);
				CHECK(result.output_size == 0);
			}
		}
		for (auto operation : realOperations)
			CHECK(operation(Vector<Real>{1, invalid}, config).status == AlgorithmStatus::InvalidInput);
		CHECK(DFT::ForwardDetailed(Vector<Real>{1, invalid}, config).status == AlgorithmStatus::InvalidInput);
		CHECK(DFT::InverseRealDetailed(Vector<Complex>{Complex(1, 0), Complex(0, invalid)}, config).status == AlgorithmStatus::InvalidInput);
	}
	config.check_finite = false;
	CHECK(FFT::ForwardDetailed(Vector<Complex>{Complex(std::numeric_limits<Real>::quiet_NaN(), 0)}, config).IsSuccess());
	CHECK(DCT::ForwardIIDetailed(Vector<Real>{Constants::PosInf}, config).IsSuccess());
	config.check_finite = true;
	config.exception_policy = EvaluationExceptionPolicy::Propagate;
	CHECK_THROWS_AS(FFT::ForwardDetailed(Vector<Complex>{Complex(Constants::PosInf, 0)}, config), NumericInputError);
}

TEST_CASE("Fourier Detailed detects nonfinite outputs of finite inputs", "[Fourier][FourierDetailedContract]") {
	auto config = ConvertFailures();
	config.normalization = TransformNormalization::None;
	const Real huge = std::numeric_limits<Real>::max();
	for (auto operation : complexOperations) {
		const auto result = operation(Vector<Complex>{Complex(huge, 0), Complex(huge, 0)}, config);
		CHECK(result.status == AlgorithmStatus::NumericalInstability);
		CHECK_FALSE(result.error_message.empty());
	}
	for (auto operation : realOperations)
		CHECK(operation(Vector<Real>{huge, huge}, config).status == AlgorithmStatus::NumericalInstability);
	CHECK(DFT::InverseRealDetailed(Vector<Complex>{Complex(huge, 0), Complex(huge, 0)}, config).status == AlgorithmStatus::NumericalInstability);
}

TEST_CASE("Complex Fourier normalization avoids overflow before representable coefficients", "[Fourier][FourierRange]") {
	const Real dc = std::numeric_limits<Real>::max() * Real(0.75);
	for (int n : {1, 2, 3, 5, 8, 31, 32, 33, 64}) {
		for (auto mode : modes) {
			for (bool forward : {true, false}) {
				CAPTURE(n, mode, forward);
				const Real scale = FourierDetail::NormalizationScale(n, mode, forward);
				const Real amplitude = dc / (n * scale);
				const Vector<Complex> input(n, Complex(amplitude, -amplitude));
				const Real bound = SumRoundoffBound(n, Real(1));
				for (bool fast : {false, true}) {
					CAPTURE(fast);
					const auto output = fast ?
						(forward ? FFT::Forward(input, mode) : FFT::Inverse(input, mode)) :
						(forward ? DFT::Forward(input, mode) : DFT::Inverse(input, mode));
					REQUIRE(output.size() == n);
					for (int k = 0; k < n; ++k) {
						CAPTURE(k);
						REQUIRE(std::isfinite(output[k].real()));
						REQUIRE(std::isfinite(output[k].imag()));
						CHECK(IsWithinAbs(output[k].real() / dc, k == 0 ? Real(1) : Real(0), bound));
						CHECK(IsWithinAbs(output[k].imag() / dc, k == 0 ? Real(-1) : Real(0), bound));
					}
				}
				FourierConfig config;
				config.normalization = mode;
				config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
				for (auto operation : forward ?
					std::array<ComplexOperation, 2>{DFT::ForwardDetailed, FFT::ForwardDetailed} :
					std::array<ComplexOperation, 2>{DFT::InverseDetailed, FFT::InverseDetailed}) {
					const auto result = operation(input, config);
					CHECK(result.IsSuccess());
					CHECK(result.output_size == n);
					CHECK(result.normalization == mode);
				}
			}
		}
	}
}

TEST_CASE("Bluestein range scaling retains phases and convolution normalization", "[Fourier][FourierRange]") {
	const Real amplitude = std::numeric_limits<Real>::max() / 16;
	for (int n : {3, 5, 31, 33}) {
		Vector<Complex> unit(n);
		Vector<Complex> huge(n);
		for (int i = 0; i < n; ++i) {
			unit[i] = Complex(Real(i % 3 - 1) / n, Real(i % 5 - 2) / n);
			huge[i] = amplitude * unit[i];
		}
		for (auto mode : modes) {
			for (bool forward : {true, false}) {
				CAPTURE(n, mode, forward);
				const auto expected = forward ? DFT::Forward(unit, mode) : DFT::Inverse(unit, mode);
				const auto actual = forward ? FFT::Forward(huge, mode) : FFT::Inverse(huge, mode);
				for (int k = 0; k < n; ++k) {
					REQUIRE(std::isfinite(actual[k].real()));
					REQUIRE(std::isfinite(actual[k].imag()));
					CHECK(IsWithinAbs(actual[k].real() / amplitude, expected[k].real(), SumRoundoffBound(n, Real(2))));
					CHECK(IsWithinAbs(actual[k].imag() / amplitude, expected[k].imag(), SumRoundoffBound(n, Real(2))));
				}
			}
		}
	}
}

TEST_CASE("Real-input Fourier surfaces retain normalized extreme coefficients", "[Fourier][FourierRange]") {
	const Real dc = std::numeric_limits<Real>::max() * Real(0.75);
	for (int n : {2, 8, 32}) {
		for (auto mode : modes) {
			CAPTURE(n, mode);
			const Real forwardScale = FourierDetail::NormalizationScale(n, mode, true);
			const Real inverseScale = FourierDetail::NormalizationScale(n, mode, false);
			const Vector<Real> input(n, dc / (n * forwardScale));
			for (const auto& spectrum : {DFT::Forward(input, mode), RealFFT::Forward(input, mode)}) {
				for (int k = 0; k < spectrum.size(); ++k) {
					REQUIRE(std::isfinite(spectrum[k].real()));
					REQUIRE(std::isfinite(spectrum[k].imag()));
					CHECK(IsWithinAbs(spectrum[k].real() / dc, k == 0 ? Real(1) : Real(0), SumRoundoffBound(n, Real(1))));
				}
			}
			Vector<Complex> full(n, Complex(dc / (n * inverseScale), 0));
			Vector<Complex> half(n / 2 + 1, full[0]);
			for (const auto& restored : {DFT::InverseReal(full, mode), RealFFT::Inverse(half, mode)}) {
				for (int i = 0; i < n; ++i) {
					REQUIRE(std::isfinite(restored[i]));
					CHECK(IsWithinAbs(restored[i] / dc, i == 0 ? Real(1) : Real(0), SumRoundoffBound(n, Real(1))));
				}
			}
		}
	}
}

TEST_CASE("Fourier range protection does not predivide tiny inputs or change raw Transform", "[Fourier][FourierRange]") {
	const Real tiny = std::numeric_limits<Real>::denorm_min();
	for (int n : {2, 8}) {
		const Vector<Complex> input(n, Complex(tiny, 0));
		CHECK(FFT::Forward(input, TransformNormalization::Forward)[0] == Complex(tiny, 0));
		CHECK(FFT::Inverse(input)[0] == Complex(tiny, 0));
		CHECK(DFT::Forward(input, TransformNormalization::Forward)[0] == Complex(tiny, 0));
		CHECK(DFT::Inverse(input)[0] == Complex(tiny, 0));
		CHECK(RealFFT::Forward(Vector<Real>(n, tiny), TransformNormalization::Forward)[0] == Complex(tiny, 0));
		const auto inverse = RealFFT::Inverse(Vector<Complex>(n / 2 + 1, Complex(tiny, 0)), TransformNormalization::None);
		CHECK(inverse[0] == n * tiny);
		Vector<Complex> impulse(n / 2 + 1, Complex(0, 0));
		impulse[0] = Complex(2 * tiny, 0);
		const auto rawInverse = RealFFT::Inverse(impulse, TransformNormalization::None);
		for (Real sample : rawInverse) CHECK(sample == 2 * tiny);

		const Real largest = std::numeric_limits<Real>::max();
		const Vector<Complex> maximum(n, Complex(largest, 0));
		CHECK(FFT::Forward(maximum, TransformNormalization::Forward)[0] == Complex(largest, 0));
		CHECK(FFT::Inverse(maximum)[0] == Complex(largest, 0));
		CHECK(DFT::Forward(maximum, TransformNormalization::Forward)[0] == Complex(largest, 0));
		CHECK(DFT::Inverse(maximum)[0] == Complex(largest, 0));
		CHECK(RealFFT::Forward(Vector<Real>(n, largest), TransformNormalization::Forward)[0] == Complex(largest, 0));
	}
	for (int n : {2, 3, 8}) {
		const Real dc = std::numeric_limits<Real>::max() / 4;
		Vector<Complex> input(n, Complex(dc / n, 0));
		for (int sign : {1, -1}) {
			auto output = input;
			FFT::Transform(output, sign);
			CHECK(std::isfinite(output[0].real()));
			CHECK(IsWithinAbs(output[0].real() / dc, Real(1), SumRoundoffBound(n, Real(1))));
		}
	}

}

TEST_CASE("Range protection retains the existing Detailed nonfinite policy", "[Fourier][FourierRange]") {
	for (auto operation : complexOperations) {
		auto config = ConvertFailures();
		config.normalization = TransformNormalization::None;
		const Vector<Complex> unrepresentable(2, Complex(std::numeric_limits<Real>::max(), 0));
		CHECK(operation(unrepresentable, config).status == AlgorithmStatus::NumericalInstability);
		config.exception_policy = EvaluationExceptionPolicy::Propagate;
		CHECK_THROWS_AS(operation(unrepresentable, config), NumericalMethodError);
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		config.check_finite = false;
		const auto unchecked = operation(unrepresentable, config);
		CHECK(unchecked.IsSuccess());
		CHECK_FALSE(std::isfinite(unchecked.value[0].real()));
	}
}

TEST_CASE("Range-protected FFT Detailed padding uses its effective transform length", "[Fourier][FourierRange]") {
	const Real huge = std::numeric_limits<Real>::max() * Real(0.75);
	const Vector<Complex> input{Complex(huge, 0), Complex(huge, 0)};
	for (int target : {4, 5}) {
		FourierConfig config;
		config.allow_padding = true;
		config.target_size = target;
		config.normalization = TransformNormalization::Forward;
		config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
		const auto result = FFT::ForwardDetailed(input, config);
		REQUIRE(result.IsSuccess());
		CHECK(result.padding_applied);
		CHECK(result.input_size == 2);
		CHECK(result.output_size == target);
		CHECK(IsWithinAbs(result.value[0].real() / huge, Real(2) / target, SumRoundoffBound(target, Real(1))));
	}
}

TEST_CASE("FFT Detailed padding is opt-in explicit and normalized at the actual output length", "[Fourier][FourierDetailedContract]") {
	const Vector<Complex> input{Complex(1, 2), Complex(3, -1), Complex(-2, 1)};
	CHECK(FFT::ForwardDetailed(input).output_size == 3);
	for (auto mode : modes) {
		FourierConfig config;
		config.allow_padding = true;
		config.target_size = 8;
		config.normalization = mode;
		Vector<Complex> padded(8, Complex(0, 0));
		for (int i = 0; i < input.size(); ++i) padded[i] = input[i];
		const auto forward = FFT::ForwardDetailed(input, config);
		const auto inverse = FFT::InverseDetailed(input, config);
		REQUIRE(forward.value.size() == 8);
		REQUIRE(inverse.value.size() == 8);
		CHECK(forward.input_size == 3);
		CHECK(forward.output_size == 8);
		CHECK(forward.padding_applied);
		CHECK(inverse.padding_applied);
		CHECK(forward.normalization == mode);
		const auto expectedForward = DFT::Forward(padded, mode);
		const auto expectedInverse = DFT::Inverse(padded, mode);
		Real inputL1 = 0;
		for (const auto& value : input) inputL1 += std::abs(value);
		const Real comparisonBound = SumRoundoffBound(8, inputL1);
		for (int i = 0; i < 8; ++i) {
			CHECK(std::abs(forward.value[i] - expectedForward[i]) < comparisonBound);
			CHECK(std::abs(inverse.value[i] - expectedInverse[i]) < comparisonBound);
		}
	}
	FourierConfig config;
	config.allow_padding = true;
	CHECK(FFT::ForwardDetailed(input, config).output_size == 4);
	config.target_size = 5;
	CHECK(FFT::ForwardDetailed(input, config).output_size == 5);
	config.target_size = 3;
	CHECK_FALSE(FFT::ForwardDetailed(input, config).padding_applied);
	CHECK(input.size() == 3);
	CHECK(input[0] == Complex(1, 2));
}

TEST_CASE("Fourier Detailed rejects unsupported options instead of ignoring them", "[Fourier][FourierDetailedContract]") {
	auto config = ConvertFailures();
	config.allow_padding = true;
	config.target_size = 8;
	CHECK(DFT::ForwardDetailed(Vector<Real>{1, 2, 3}, config).status == AlgorithmStatus::InvalidInput);
	CHECK(DFT::ForwardDetailed(Vector<Complex>{Complex(1, 0)}, config).status == AlgorithmStatus::InvalidInput);
	CHECK(DFT::InverseDetailed(Vector<Complex>{Complex(1, 0)}, config).status == AlgorithmStatus::InvalidInput);
	CHECK(DFT::InverseRealDetailed(Vector<Complex>{Complex(1, 0)}, config).status == AlgorithmStatus::InvalidInput);
	for (auto operation : realOperations)
		CHECK(operation(Vector<Real>{1, 2, 3}, config).status == AlgorithmStatus::InvalidInput);
	config.allow_padding = false;
	CHECK(FFT::ForwardDetailed(Vector<Complex>{Complex(1, 0)}, config).status == AlgorithmStatus::InvalidInput);
	config.allow_padding = true;
	config.target_size = -1;
	CHECK(FFT::InverseDetailed(Vector<Complex>{Complex(1, 0)}, config).status == AlgorithmStatus::InvalidInput);
	config.target_size = 1;
	CHECK(FFT::ForwardDetailed(Vector<Complex>{Complex(1, 0), Complex(2, 0)}, config).status == AlgorithmStatus::InvalidInput);
	config.target_size = 0;
	config.allow_padding = false;
	config.normalization = static_cast<TransformNormalization>(-1);
	for (auto operation : complexOperations)
		CHECK(operation(Vector<Complex>{Complex(1, 0)}, config).status == AlgorithmStatus::InvalidInput);
	for (auto operation : realOperations)
		CHECK(operation(Vector<Real>{1}, config).status == AlgorithmStatus::InvalidInput);
}

TEST_CASE("DCT orthonormal Detailed coefficients use unit-norm cosine rows", "[Fourier][FourierDetailedContract]") {
	FourierConfig config;
	config.normalization = TransformNormalization::Orthonormal;
	const auto result = DCT::ForwardIIDetailed(Vector<Real>(8, Real(1)), config);
	CHECK(result.output_size == 8);
	CHECK(result.normalization == TransformNormalization::Orthonormal);
	CHECK(IsWithinAbs(result.value[0], std::sqrt(Real(8)), Precision::EqualityTolerance));
	for (int i = 1; i < 8; ++i)
		CHECK(std::abs(result.value[i]) < Precision::EqualityTolerance);
	CHECK(IsWithinAbs(DCT::ForwardII(Vector<Real>(8, Real(1)))[0], Real(2), Precision::EqualityTolerance));
}

TEST_CASE("Normalized singleton DCT preserves finite extremes and subnormal coefficients", "[Fourier][FourierDetailedContract]") {
	for (Real value : {std::numeric_limits<Real>::max(), std::numeric_limits<Real>::denorm_min()}) {
		const Vector<Real> data{value};
		CHECK(DCT::ForwardII(data, TransformNormalization::Forward)[0] == value);
		CHECK(DCT::ForwardII(data, TransformNormalization::Orthonormal)[0] == value);
		CHECK(DCT::InverseII(data, TransformNormalization::Orthonormal)[0] == value);
	}
}

TEST_CASE("DCT and DST Detailed modes match independent sums across direct and FFT dispatch", "[Fourier][FourierDetailedContract]") {
	for (int n : {1, 2, 7, 8, 31, 32, 33, 64}) {
		Vector<Real> data(n);
		for (int i = 0; i < n; ++i) data[i] = std::sin(Real(i + 1)) + Real(i) / n;
		for (const auto mode : modes) {
			CAPTURE(n, mode);
			FourierConfig config;
			config.normalization = mode;
			const auto cosine = DCT::ForwardIIDetailed(data, config);
			const auto sine = DCT::ForwardDSTDetailed(data, config);
			REQUIRE(cosine.output_size == n);
			REQUIRE(sine.output_size == n);
			CHECK(cosine.normalization == mode);
			CHECK(sine.normalization == mode);
			Real energy = 0, cosineEnergy = 0, sineEnergy = 0;
			Real inputL1 = 0;
			for (const Real value : data) inputL1 += std::abs(value);
			for (int k = 0; k < n; ++k) {
				long double rawCosine = 0, rawSine = 0;
				for (int j = 0; j < n; ++j) {
					const long double pi = std::acos(-1.L);
					rawCosine += 2.L * data[j] * std::cos(pi * k * (j + 0.5L) / n);
					rawSine += 2.L * data[j] * std::sin(pi * (k + 1) * (j + 1.L) / (n + 1.L));
				}

				Real cosineScale = 1, sineScale = 1;
				if (mode == TransformNormalization::Legacy) cosineScale = Real(1) / n;
				if (mode == TransformNormalization::Forward) {
					cosineScale = Real(1) / (2 * Real(n));
					sineScale = Real(1) / (2 * (Real(n) + 1));
				}
				if (mode == TransformNormalization::Orthonormal)
					cosineScale = Real(1) / std::sqrt((k == 0 ? Real(4) : Real(2)) * n);
				if (mode == TransformNormalization::Orthonormal || mode == TransformNormalization::Legacy)
					sineScale = Real(1) / std::sqrt(2 * (Real(n) + 1));
				CHECK(IsWithinAbs(cosine.value[k], Real(rawCosine * cosineScale),
					SumRoundoffBound(n, 2 * inputL1) * cosineScale));
				CHECK(IsWithinAbs(sine.value[k], Real(rawSine * sineScale),
					SumRoundoffBound(n, 2 * inputL1) * sineScale));
				energy += data[k] * data[k];
				cosineEnergy += cosine.value[k] * cosine.value[k];
				sineEnergy += sine.value[k] * sine.value[k];
			}
			if (mode == TransformNormalization::Orthonormal) {
				const Real energyBound = std::max(n * Precision::EqualityTolerance,
					Real(16) * n * std::numeric_limits<Real>::epsilon() * energy);
				CHECK(IsWithinAbs(cosineEnergy, energy, energyBound));
				CHECK(IsWithinAbs(sineEnergy, energy, energyBound));
			}
			const auto cosineInverse = DCT::InverseIIDetailed(cosine.value, config);
			const auto sineInverse = DCT::InverseDSTDetailed(sine.value, config);
			REQUIRE(cosineInverse.output_size == n);
			REQUIRE(sineInverse.output_size == n);
			for (int i = 0; i < n; ++i) {
				const Real cosineFactor = mode == TransformNormalization::None ? 2 * Real(n) : Real(1);
				const Real sineFactor = mode == TransformNormalization::None ? 2 * (Real(n) + 1) : Real(1);
				CHECK(IsWithinAbs(cosineInverse.value[i], cosineFactor * data[i],
					cosineFactor * SumRoundoffBound(n, inputL1)));
				CHECK(IsWithinAbs(sineInverse.value[i], sineFactor * data[i],
					sineFactor * SumRoundoffBound(n, inputL1)));
			}
		}
	}
}

TEST_CASE("Fourier Detailed work metadata follows direct radix2 and Bluestein kernels", "[Fourier][FourierDetailedContract]") {
	for (int size : {1, 8, 31, 32, 33}) {
		Vector<Complex> complex(size, Complex(1, 0));
		const auto dft = DFT::ForwardDetailed(complex);
		const auto fft = FFT::ForwardDetailed(complex);
		CHECK(dft.kernel_name == "Direct");
		CHECK(dft.work_size == size);
		CHECK(dft.work_units == std::uint64_t(size) * size);
		CHECK(dft.function_evaluations == 0);
		CHECK(fft.function_evaluations == 0);
		const bool radix2 = size == 1 || size == 8 || size == 32;
		CHECK(fft.kernel_name == (radix2 ? "Radix2" : "Bluestein"));
		const int workSize = size == 31 ? 64 : size == 33 ? 128 : size;
		int stages = 0;
		for (int n = workSize; n > 1; n /= 2) ++stages;
		const auto butterflies = std::uint64_t(workSize / 2) * stages;
		CHECK(fft.work_size == workSize);
		CHECK(fft.work_units == (radix2 ? butterflies : 3 * butterflies + workSize + 2 * std::uint64_t(size)));
		const auto inverseFFT = FFT::InverseDetailed(complex);
		CHECK(inverseFFT.work_size == fft.work_size);
		CHECK(inverseFFT.work_units == fft.work_units);
		Vector<Real> real(size, Real(1));
		const auto cosine = DCT::ForwardIIDetailed(real);
		const auto inverseCosine = DCT::InverseIIDetailed(real);
		const auto sine = DCT::ForwardDSTDetailed(real);
		CHECK(cosine.function_evaluations == 0);
		CHECK(sine.function_evaluations == 0);
		CHECK(sine.kernel_name == "Direct");
		CHECK(sine.work_units == std::uint64_t(size) * size);
		if (size < 32) {
			CHECK(cosine.kernel_name == "Direct");
			CHECK(cosine.work_size == size);
			CHECK(cosine.work_units == std::uint64_t(size) * size);
		} else if (size == 32) {
			CHECK(cosine.kernel_name == "FFT/Radix2");
			CHECK(cosine.work_size == 64);
			CHECK(cosine.work_units == 224);
		} else {
			CHECK(cosine.kernel_name == "FFT/Bluestein");
			CHECK(cosine.work_size == 256);
			CHECK(cosine.work_units == 3493);
		}
		CHECK(inverseCosine.kernel_name == cosine.kernel_name);
		CHECK(inverseCosine.work_units == cosine.work_units);
	}
}

TEST_CASE("Fourier work-size and counter boundaries are checked without allocating huge inputs", "[Fourier][FourierDetailedContract]") {
	const int max = std::numeric_limits<int>::max();
	CHECK(FourierDetail::DoubledSize(max / 2) == max - 1);
	CHECK_THROWS_AS(FourierDetail::DoubledSize(max / 2 + 1), ResourceLimitError);
	const int largestPower = max / 2 + 1;
	CHECK(FFT::NextPowerOfTwo(largestPower) == largestPower);
	CHECK_THROWS_AS(FFT::NextPowerOfTwo(largestPower + 1), FourierError);
	CHECK(FourierDetail::BluesteinSize(largestPower / 2) == largestPower);
	CHECK_THROWS_AS(FourierDetail::BluesteinSize(largestPower / 2 + 1), ResourceLimitError);
	CHECK(FourierDetail::DirectWork(max).units == std::uint64_t(max) * std::uint64_t(max));
	CHECK(FourierDetail::FFTWork(largestPower).units > std::uint64_t(max));
	auto config = ConvertFailures();
	config.allow_padding = true;
	config.target_size = max;
	const auto rejected = FFT::ForwardDetailed(Vector<Complex>{Complex(1, 0)}, config);
	CHECK(rejected.status == AlgorithmStatus::ResourceLimitExceeded);
	CHECK(rejected.input_size == 1);
	CHECK(rejected.output_size == 0);
	CHECK(rejected.work_units == 0);
	CHECK_FALSE(rejected.padding_applied);
	CHECK_FALSE(rejected.error_message.empty());
	Vector<Complex> data{Complex(1, 0), Complex(2, 0)};
	CHECK_THROWS_AS(FFT::Transform(data, std::numeric_limits<int>::min()), FourierError);
	CHECK(data[0] == Complex(1, 0));
}

TEST_CASE("Fourier empty failures retain requested metadata and no invented work", "[Fourier][FourierDetailedContract]") {
	auto config = ConvertFailures();
	config.normalization = TransformNormalization::Orthonormal;
	auto check = [](const auto& result) {
		CHECK(result.status == AlgorithmStatus::InvalidInput);
		CHECK(result.normalization == TransformNormalization::Orthonormal);
		CHECK(result.input_size == 0);
		CHECK(result.output_size == 0);
		CHECK(result.work_size == 0);
		CHECK(result.work_units == 0);
		CHECK(result.function_evaluations == 0);
		CHECK_FALSE(result.algorithm_name.empty());
		CHECK_FALSE(result.error_message.empty());
	};
	for (auto operation : complexOperations) check(operation(Vector<Complex>{}, config));
	for (auto operation : realOperations) check(operation(Vector<Real>{}, config));
	check(DFT::ForwardDetailed(Vector<Real>{}, config));
	check(DFT::InverseRealDetailed(Vector<Complex>{}, config));
}
