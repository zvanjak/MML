///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Fourier/FourierSpectrum.h                                           ///
///  Description: Small spectrum helper utilities for Fourier transform outputs       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///               Copyright (c) 2024-2026 Zvonimir Vanjak                             ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FOURIER_SPECTRUM_H
#define MML_FOURIER_SPECTRUM_H

#include <mml/MMLBase.h>

#include <mml/algorithms/Fourier/Fourier.h>

#include <cmath>
#include <stdexcept>

namespace MML::Fourier
{
	enum class FrequencyAxisLayout {
		TwoSided,
		OneSided
	};

	inline Vector<Real> Magnitude(const Vector<Complex>& spectrum) {
		FourierValidation::ValidateNonEmpty(spectrum.size(), "Magnitude");

		Vector<Real> result(spectrum.size());
		for (int i = 0; i < spectrum.size(); i++)
			result[i] = std::abs(spectrum[i]);

		return result;
	}

	inline Vector<Real> MagnitudeSquared(const Vector<Complex>& spectrum) {
		FourierValidation::ValidateNonEmpty(spectrum.size(), "MagnitudeSquared");

		Vector<Real> result(spectrum.size());
		for (int i = 0; i < spectrum.size(); i++)
			result[i] = std::norm(spectrum[i]);

		return result;
	}

	inline Vector<Real> Power(const Vector<Complex>& spectrum) {
		return MagnitudeSquared(spectrum);
	}

	inline Vector<Real> Phase(const Vector<Complex>& spectrum) {
		FourierValidation::ValidateNonEmpty(spectrum.size(), "Phase");

		Vector<Real> result(spectrum.size());
		for (int i = 0; i < spectrum.size(); i++)
			result[i] = std::atan2(spectrum[i].imag(), spectrum[i].real());

		return result;
	}

	inline Vector<Real> UnwrapPhase(const Vector<Real>& phase, Real discontinuity = Constants::PI) {
		FourierValidation::ValidateNonEmpty(phase.size(), "UnwrapPhase");
		FourierValidation::ValidatePositive(discontinuity, "discontinuity", "UnwrapPhase");

		Vector<Real> result = phase;
		Real offset = 0.0;
		const Real period = 2.0 * Constants::PI;

		for (int i = 1; i < result.size(); i++) {
			Real delta = phase[i] - phase[i - 1];
			if (delta > discontinuity)
				offset -= period;
			else if (delta < -discontinuity)
				offset += period;

			result[i] = phase[i] + offset;
		}

		return result;
	}

	inline Vector<Real> FrequencyAxis(int n,
	                                  Real sample_rate,
	                                  FrequencyAxisLayout layout = FrequencyAxisLayout::TwoSided,
	                                  bool shifted = false) {
		FourierValidation::ValidatePositiveSize(n, "n", "FrequencyAxis");
		FourierValidation::ValidatePositive(sample_rate, "sample_rate", "FrequencyAxis");

		Real df = sample_rate / n;
		if (layout == FrequencyAxisLayout::OneSided) {
			int half_size = n / 2 + 1;
			Vector<Real> frequencies(half_size);
			for (int i = 0; i < half_size; i++)
				frequencies[i] = i * df;

			return frequencies;
		}

		Vector<Real> frequencies(n);
		if (shifted) {
			int first_bin = -n / 2;
			for (int i = 0; i < n; i++)
				frequencies[i] = (first_bin + i) * df;
		} else {
			int positive_count = (n + 1) / 2;
			for (int i = 0; i < n; i++) {
				int bin = (i < positive_count) ? i : i - n;
				frequencies[i] = bin * df;
			}
		}

		return frequencies;
	}

	template<typename T>
	Vector<T> FFTShift(const Vector<T>& data) {
		FourierValidation::ValidateNonEmpty(data.size(), "FFTShift");

		int n = data.size();
		int shift = (n + 1) / 2;
		Vector<T> result(n);
		for (int i = 0; i < n; i++)
			result[i] = data[(i + shift) % n];

		return result;
	}

	template<typename T>
	Vector<T> IFFTShift(const Vector<T>& data) {
		FourierValidation::ValidateNonEmpty(data.size(), "IFFTShift");

		int n = data.size();
		int shift = n / 2;
		Vector<T> result(n);
		for (int i = 0; i < n; i++)
			result[i] = data[(i + shift) % n];

		return result;
	}
} // namespace MML::Fourier

#endif // MML_FOURIER_SPECTRUM_H