///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DerivationBase.h                                                    ///
///  Description: Base definitions for numerical differentiation                      ///
///               Step size control and derivative types                              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_BASE_H
#define MML_DERIVATION_BASE_H

#include <mml/MMLBase.h>
#include <mml/MMLNumericValidation.h>
#include <algorithm>

#include <mml/interfaces/IFunction.h>

namespace MML
{
	namespace Derivation
	{
		/// @brief Base step for a first derivative with first-order accuracy (forward/backward difference)
		/// @note h = 2√ε minimizes roundoff + truncation error for O(h) methods
		static inline const Real NDer1_h = 2 * std::sqrt(Constants::Eps);
		/// @brief Base step for a first derivative with second-order accuracy (central difference)
		/// @note h = (3ε)^(1/3) for O(h²) methods, balances O(ε/h) roundoff with O(h²) truncation
		static inline const Real NDer2_h = std::pow(3 * Constants::Eps, MML::Real(1.0) / MML::Real(3.0));   // tol-ok: optimal-step constants from the truncation/roundoff balance (NR 5.7)
		/// @brief Base step for a first derivative with fourth-order accuracy (5-point stencil)
		/// @note h = (11.25ε)^(1/5) ≈ 0.0012 for double precision
		static inline const Real NDer4_h = std::pow(MML::Real(11.25) * Constants::Eps, MML::Real(1.0) / MML::Real(5.0));   // tol-ok: optimal-step constant     // 0.0012009323661373839 for double!
		/// @brief Base step for a first derivative with sixth-order accuracy (7-point stencil)
		/// @note h = (ε/168)^(1/7) for O(h⁶) methods
		static inline const Real NDer6_h = std::pow(Constants::Eps / MML::Real(168.0), MML::Real(1.0) / MML::Real(7.0));
		/// @brief Base step for a first derivative with eighth-order accuracy (9-point stencil)
		/// @note h = (551.25ε)^(1/9) for O(h⁸) methods
		static inline const Real NDer8_h = std::pow(MML::Real(551.25) * Constants::Eps, MML::Real(1.0) / MML::Real(9.0));   // tol-ok: optimal-step constant

		/// @brief Scale step size relative to |x| to avoid cancellation at large x
		/// @details For |x| > 1 scales h proportionally; for |x| <= 1 returns h unchanged
		static inline Real ScaleStep(Real h, Real x) {
			return h * std::max(Real(1), std::abs(x));
		}

		namespace Detail
		{
			template<class Sample, class Value>
			Value BoundedFirstDerivative(Sample sample, const Value& center,
				Real t, Real lower, Real upper, Real step)
			{
				if (!std::isfinite(step) || step <= 0)
					throw NumericalMethodError("Bounded derivative: step must be positive and finite");
				const Real centralStep = std::min({step, t - lower, upper - t});
				Real first = std::max(lower, t - centralStep);
				Real second = std::min(upper, t + centralStep);
				if (!(first < t && t < second)) {
					const bool forward = upper - t >= t - lower;
					const Real end = forward ? upper : lower;
					const Real h = std::min(step, std::abs(end - t) / 2);
					const Real direction = forward ? Real(1) : Real(-1);
					first = std::clamp(t + direction * h, lower, upper);
					second = std::clamp(t + direction * (2 * h), lower, upper);
					if (first == t) first = std::nextafter(t, end);
					if (second == first || second == t) second = std::nextafter(first, end);
				}
				if (first == t || second == t || first == second ||
					first < lower || first > upper || second < lower || second > upper ||
					!std::isfinite(first - t) || !std::isfinite(second - t) || !std::isfinite(second - first))
					throw NumericalMethodError("Bounded derivative: insufficient distinct in-domain nodes");
				const auto firstSlope = (sample(first) - center) / (first - t);
				const auto secondSlope = (sample(second) - center) / (second - t);
				// Quadratic interpolation at actual representable nodes, including one-sided stencils.
				return firstSlope * ((second - t) / (second - first))
					- secondSlope * ((first - t) / (second - first));
			}

			// Taylor moments: these residuals approximate f^(m+p) h^(m+p), leaving h^p after division by h^m.
			template<typename Value, typename Magnitude>
			Real EstimateSecondTruncation2(const Value& ym2, const Value& ym, const Value& y0,
				const Value& yp, const Value& yp2, Real h, Magnitude magnitude) {
				return magnitude(ym2 - 4 * ym + 6 * y0 - 4 * yp + yp2) / h / h / 12;
			}

			template<typename Value, typename Magnitude>
			Real EstimateSecondTruncation4(const Value& ym3, const Value& ym2, const Value& ym,
				const Value& y0, const Value& yp, const Value& yp2, const Value& yp3, Real h, Magnitude magnitude) {
				return magnitude(ym3 - 6 * ym2 + 15 * ym - 20 * y0 + 15 * yp - 6 * yp2 + yp3) / h / h / 90;
			}

			template<typename Value, typename Magnitude>
			Real EstimateThirdTruncation2(const Value& ym3, const Value& ym2, const Value& ym,
				const Value& yp, const Value& yp2, const Value& yp3, Real h, Magnitude magnitude) {
				return magnitude(-ym3 + 4 * ym2 - 5 * ym + 5 * yp - 4 * yp2 + yp3) / h / h / h / 8;
			}

			template<typename Value, typename Magnitude>
			Real EstimateThirdTruncation4(const Value& ym4, const Value& ym3, const Value& ym2,
				const Value& ym, const Value& yp, const Value& yp2, const Value& yp3, const Value& yp4,
				Real h, Magnitude magnitude) {
				return magnitude(-ym4 + 6 * ym3 - 14 * ym2 + 14 * ym - 14 * yp + 14 * yp2 - 6 * yp3 + yp4)
					/ h / h / h * (Real(7) / 240);
			}

			inline Real RefineMixedDerivativeStep(Real coordinate, Real step) {
				ValidateFinite(coordinate, "Mixed derivative coordinate");
				const Real half = step / 2;
				if (!std::isfinite(coordinate - step) || !std::isfinite(coordinate + step) ||
					!(coordinate - step < coordinate - half && coordinate - half < coordinate &&
					  coordinate < coordinate + half && coordinate + half < coordinate + step))
					throw NumericalMethodError("Mixed derivative error estimate: half-step coordinates collapse or overflow");
				return half;
			}

			template<int AccuracyOrder, typename Value, typename Magnitude>
			Real EstimateMixedDerivativeError(const Value& coarse, const Value& fine,
				Real coarseRoundoff, Real fineRoundoff, Magnitude magnitude) {
				const Real factor = AccuracyOrder == 2 ? Real(4) / 3 : Real(16) / 15;
				const Real estimate = factor * (magnitude(coarse - fine) + coarseRoundoff + fineRoundoff);
				if (!std::isfinite(estimate))
					throw NumericalMethodError("Mixed derivative: non-finite error estimate");
				return estimate;
			}
		}
	}
}

#endif // MML_DERIVATION_BASE_H
