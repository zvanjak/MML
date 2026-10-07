///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Fourier/FourierConvolution.h                                        ///
///  Description: Core linear and circular convolution helpers                        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                             ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FOURIER_CONVOLUTION_H
#define MML_FOURIER_CONVOLUTION_H

#include <mml/MMLBase.h>

#include <mml/algorithms/Fourier/Fourier.h>

#include <algorithm>
#include <stdexcept>

namespace MML::Fourier
{
	enum class ConvolutionMode {
		Full,
		Same,
		Valid,
		Circular
	};

	enum class ConvolutionMethod {
		Auto,
		Direct,
		FFT
	};

	namespace ConvolutionDetail
	{
		template<typename T>
		Vector<T> DirectFull(const Vector<T>& x, const Vector<T>& h, const char* function_name) {
			FourierValidation::ValidateNonEmpty(x.size(), std::string(function_name) + " (x)");
			FourierValidation::ValidateNonEmpty(h.size(), std::string(function_name) + " (h)");

			Vector<T> result(x.size() + h.size() - 1, T{});
			for (int i = 0; i < x.size(); i++)
				for (int j = 0; j < h.size(); j++)
					result[i + j] += x[i] * h[j];

			return result;
		}

		inline Vector<Complex> FFTFullComplex(const Vector<Complex>& x, const Vector<Complex>& h, const char* function_name) {
			FourierValidation::ValidateNonEmpty(x.size(), std::string(function_name) + " (x)");
			FourierValidation::ValidateNonEmpty(h.size(), std::string(function_name) + " (h)");

			int result_size = x.size() + h.size() - 1;
			int fft_size = FFT::NextPowerOfTwo(result_size);
			Vector<Complex> x_padded(fft_size, Complex(0.0, 0.0));
			Vector<Complex> h_padded(fft_size, Complex(0.0, 0.0));

			for (int i = 0; i < x.size(); i++)
				x_padded[i] = x[i];
			for (int i = 0; i < h.size(); i++)
				h_padded[i] = h[i];

			FFT::Transform(x_padded, 1);
			FFT::Transform(h_padded, 1);
			for (int i = 0; i < fft_size; i++)
				x_padded[i] *= h_padded[i];

			FFT::Transform(x_padded, -1);
			Real norm = 1.0 / fft_size;
			Vector<Complex> result(result_size);
			for (int i = 0; i < result_size; i++)
				result[i] = x_padded[i] * norm;

			return result;
		}

		inline Vector<Complex> ToComplex(const Vector<Real>& data) {
			Vector<Complex> result(data.size());
			for (int i = 0; i < data.size(); i++)
				result[i] = Complex(data[i], 0.0);

			return result;
		}

		inline Vector<Real> RealPart(const Vector<Complex>& data) {
			Vector<Real> result(data.size());
			for (int i = 0; i < data.size(); i++)
				result[i] = data[i].real();

			return result;
		}

		template<typename T>
		Vector<T> SliceLinearResult(const Vector<T>& full, int x_size, int h_size, ConvolutionMode mode) {
			if (mode == ConvolutionMode::Full)
				return full;

			int start = 0;
			int result_size = 0;
			if (mode == ConvolutionMode::Same) {
				start = (h_size - 1) / 2;
				result_size = x_size;
			} else if (mode == ConvolutionMode::Valid) {
				int min_size = std::min(x_size, h_size);
				int max_size = std::max(x_size, h_size);
				start = min_size - 1;
				result_size = max_size - min_size + 1;
			} else {
				throw FourierError("SliceLinearResult - circular mode is not a linear slice");
			}

			Vector<T> result(result_size);
			for (int i = 0; i < result_size; i++)
				result[i] = full[start + i];

			return result;
		}

		template<typename T>
		Vector<T> FoldCircular(const Vector<T>& full, int size) {
			Vector<T> result(size, T{});
			for (int i = 0; i < full.size(); i++)
				result[i % size] += full[i];

			return result;
		}

		template<typename T>
		Vector<T> DirectCircular(const Vector<T>& x, const Vector<T>& h, const char* function_name) {
			FourierValidation::ValidateNonEmpty(x.size(), std::string(function_name) + " (x)");
			FourierValidation::ValidateNonEmpty(h.size(), std::string(function_name) + " (h)");

			int result_size = std::max(x.size(), h.size());
			Vector<T> result(result_size, T{});
			for (int i = 0; i < x.size(); i++)
				for (int j = 0; j < h.size(); j++)
					result[(i + j) % result_size] += x[i] * h[j];

			return result;
		}

		inline bool ShouldUseFFT(int x_size, int h_size, ConvolutionMethod method) {
			if (method == ConvolutionMethod::FFT)
				return true;
			if (method == ConvolutionMethod::Direct)
				return false;

			return static_cast<long long>(x_size) * static_cast<long long>(h_size) > 4096;
		}
	} // namespace ConvolutionDetail

	inline Vector<Complex> Convolve(const Vector<Complex>& x,
	                               const Vector<Complex>& h,
	                               ConvolutionMode mode = ConvolutionMode::Full,
	                               ConvolutionMethod method = ConvolutionMethod::Auto) {
		if (mode == ConvolutionMode::Circular) {
			if (ConvolutionDetail::ShouldUseFFT(x.size(), h.size(), method)) {
				auto full = ConvolutionDetail::FFTFullComplex(x, h, "Convolve");
				return ConvolutionDetail::FoldCircular(full, std::max(x.size(), h.size()));
			}

			return ConvolutionDetail::DirectCircular(x, h, "Convolve");
		}

		Vector<Complex> full = ConvolutionDetail::ShouldUseFFT(x.size(), h.size(), method)
			? ConvolutionDetail::FFTFullComplex(x, h, "Convolve")
			: ConvolutionDetail::DirectFull(x, h, "Convolve");

		return ConvolutionDetail::SliceLinearResult(full, x.size(), h.size(), mode);
	}

	inline Vector<Real> Convolve(const Vector<Real>& x,
	                            const Vector<Real>& h,
	                            ConvolutionMode mode = ConvolutionMode::Full,
	                            ConvolutionMethod method = ConvolutionMethod::Auto) {
		if (mode == ConvolutionMode::Circular && method != ConvolutionMethod::FFT)
			return ConvolutionDetail::DirectCircular(x, h, "Convolve");

		Vector<Complex> complex_result = Convolve(ConvolutionDetail::ToComplex(x), ConvolutionDetail::ToComplex(h), mode, method);
		return ConvolutionDetail::RealPart(complex_result);
	}
} // namespace MML::Fourier

#endif // MML_FOURIER_CONVOLUTION_H