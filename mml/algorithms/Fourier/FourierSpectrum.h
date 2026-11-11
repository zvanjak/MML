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

	/// Raw squared coefficient magnitude; no window, sample-rate or one-sided power normalization.
	inline Vector<Real> MagnitudeSquared(const Vector<Complex>& spectrum) {
		FourierValidation::ValidateNonEmpty(spectrum.size(), "MagnitudeSquared");

		Vector<Real> result(spectrum.size());
		for (int i = 0; i < spectrum.size(); i++)
			result[i] = std::norm(spectrum[i]);

		return result;
	}

	inline Vector<Real> Phase(const Vector<Complex>& spectrum) {
		FourierValidation::ValidateNonEmpty(spectrum.size(), "Phase");

		Vector<Real> result(spectrum.size());
		for (int i = 0; i < spectrum.size(); i++)
			result[i] = std::atan2(spectrum[i].imag(), spectrum[i].real());

		return result;
	}

	/// Removes the minimum whole turns needed to cross the selected jump threshold.
	/// Accepts arbitrary finite radians; thresholds below pi can produce opposite-sign jumps exceeding the threshold.
	inline Vector<Real> UnwrapPhase(const Vector<Real>& phase, Real discontinuity = Constants::PI) {
		FourierValidation::ValidateNonEmpty(phase.size(), "UnwrapPhase");
		FourierValidation::ValidatePositive(discontinuity, "discontinuity", "UnwrapPhase");
		FourierDetail::ValidateFinite(phase, true);

		Vector<Real> result = phase;
		const Real period = 2.0 * Constants::PI;

		for (int i = 1; i < result.size(); i++) {
			Real delta = phase[i] - phase[i - 1];
			if (delta > discontinuity || delta < -discontinuity) {
				// Reduce before subtracting: huge phase differences can overflow or lose whole turns.
				Real reduced = std::fmod(phase[i], period) - std::fmod(phase[i - 1], period);
				if (delta < 0)
					reduced = -reduced;
				Real excess = std::fmod(reduced - std::fmod(discontinuity, period), period);
				if (excess < 0)
					excess += period;
				const Real adjusted = discontinuity - (excess == 0 ? Real{0} : period - excess);
				delta = delta > 0 ? adjusted : -adjusted;
			}
			result[i] = result[i - 1] + delta;
			if (!std::isfinite(result[i]))
				throw NumericalMethodError("UnwrapPhase - unwrapped output exceeds the finite range");
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