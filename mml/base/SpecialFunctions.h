///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        SpecialFunctions.h                                                  ///
///  Description: Public special-function library (gamma/beta family, digamma, erf-1) ///
///               Portable implementations independent of std:: special math          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_SPECIAL_FUNCTIONS_H
#define MML_SPECIAL_FUNCTIONS_H

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <type_traits>

#include <mml/MMLBase.h>

namespace MML
{
	/// @brief Portable special functions (gamma/beta family, polygamma, error-function inverses).
	/// @details These are real-valued implementations using the configured Real type, independent of the
	///          C++17 `std::` special-math functions (which are unavailable on some platforms).
	///          The regularized incomplete gamma/beta routines and their inverses are the shared,
	///          public home for machinery previously duplicated privately inside the statistical
	///          distributions.
	namespace SpecialFunctions
	{
		namespace SpecialFunctionsDetail
		{
			using ShapeReal = std::common_type_t<Real, double>;
			// DLMF 5.11.1, through the x^-15 Bernoulli term.
			inline constexpr std::array<ShapeReal, 8> StirlingCoefficients = {ShapeReal(1) / 12, -ShapeReal(1) / 360,
				ShapeReal(1) / 1260, -ShapeReal(1) / 1680, ShapeReal(1) / 1188, -ShapeReal(691) / 360360,
				ShapeReal(1) / 156, -ShapeReal(3617) / 122400};

			inline ShapeReal LogOnePlusMinusX(ShapeReal x)
			{
				if (std::abs(x) > Precision::ForwardShapeLogSeriesLimit)
					return std::log1p(x) - x;
				ShapeReal term = -x * x / 2, sum = term;
				for (int n = 3; n <= Precision::ForwardShapeLogSeriesIterations; ++n) {
					term *= -x * ShapeReal(n - 1) / n;
					const ShapeReal next = sum + term;
					if (next == sum) return sum;
					sum = next;
				}
				throw ConvergenceError("SpecialFunctions logarithmic deviance series exhausted its budget",
					Precision::ForwardShapeLogSeriesIterations, static_cast<Real>(std::abs(term)));
			}

			inline ShapeReal StirlingCorrection(ShapeReal x)
			{
				using S = ShapeReal;
				const S inv = 1 / x, squared = inv * inv;
				S result = StirlingCoefficients[7];
				for (int i = 6; i >= 0; --i) result = StirlingCoefficients[i] + squared * result;
				return inv * result;
			}

			inline ShapeReal LogOnePlusOverX(ShapeReal x) { return x == 0 ? ShapeReal(1) : std::log1p(x) / x; }
			inline ShapeReal Expm1OverX(ShapeReal x) { return x == 0 ? ShapeReal(1) : std::expm1(x) / x; }

			inline ShapeReal GammaLogOnePlusRatio(ShapeReal a)
			{
				using S = ShapeReal;
				if (a > Precision::GammaInverseSmallShapeLimit)
					return std::lgamma(1 + a) / a;
				const S coefficients[] = {
					S(-0.577215664901532860606512090082402431L), S(0.822467033424113218236207583323012595L),
					S(-0.400685634386531428466579387170483330L), S(0.270580808427784547879000924135291976L),
					S(-0.207385551028673985266273097291406834L), S(0.169557176997408189952419654965153421L),
					S(-0.144049896768846118119971078549970966L), S(0.125509669524743042422335654813581558L),
					S(-0.111334265869564690490872529914712451L), S(0.100099457512781808533714595890031902L)
				};
				S result = coefficients[9];
				for (int i = 8; i >= 0; --i) result = coefficients[i] + a * result;
				return result;
			}

			inline ShapeReal PositiveLogBeta(ShapeReal a, ShapeReal b)
			{
				using S = ShapeReal;
				if (a > b) std::swap(a, b);
				if (b < Precision::ForwardShapeStirlingMin) {
					if (a + b <= Precision::GammaInverseSmallShapeLimit)
						return a * GammaLogOnePlusRatio(a) + b * GammaLogOnePlusRatio(b)
							- (a + b) * GammaLogOnePlusRatio(a + b)
							+ std::log1p(a / b) - std::log(a);
					return std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
				}
				const S ratio = a / b, logRatio = std::log1p(ratio);
				if (a < Precision::ForwardShapeStirlingMin) {
					return std::lgamma(a) - a * std::log(b)
						- (b + a - S(0.5)) * LogOnePlusMinusX(ratio) - (a - S(0.5)) * ratio
						+ StirlingCorrection(b) - StirlingCorrection(a + b);
				}
				const S logP = std::log(a) - std::log(b) - logRatio;
				return a * logP - b * logRatio
					+ (std::log(S(2) * Constants::PI) - std::log(a) + logRatio) / 2
					+ StirlingCorrection(a) + StirlingCorrection(b) - StirlingCorrection(a + b);
			}

			// Positive, mode-scaled densities keep quadrature conditioning independent of shape.
			// An eight-point Gauss rule is compared with its two half-interval rules.
			template<class Function>
			inline ShapeReal PositiveShapeIntegral(Function&& function, ShapeReal endpoint,
				int maxEvaluations = Precision::ForwardShapeQuadratureMaxEvaluations)
			{
				using S = ShapeReal;
				const S nodes[] = {S(0.183434642495649804939476142360183981L), S(0.525532409916328985817739049189246349L),
					S(0.796666477413626739591553936475830437L), S(0.960289856497536231683560868569472991L)};
				const S weights[] = {S(0.362683783378361982965150449277195612L), S(0.313706645877887287337962201986601313L),
					S(0.222381034453374470544355994426240884L), S(0.101228536290376259152531354309962190L)};
				int evaluations = 0;
				auto rule = [&](S low, S high) {
					const S middle = (low + high) / 2, radius = (high - low) / 2;
					S sum = 0;
					for (int i = 0; i < 4; ++i) {
						if (evaluations + 2 > maxEvaluations)
							throw ConvergenceError("SpecialFunctions shape quadrature exhausted its evaluation budget",
								evaluations, std::numeric_limits<Real>::infinity());
						evaluations += 2;
						const S pair = function(middle - radius * nodes[i]) + function(middle + radius * nodes[i]);
						if (!(std::isfinite(pair) && pair >= 0))
							throw NumericalMethodError("SpecialFunctions shape quadrature density is invalid");
						sum += weights[i] * pair;
					}
					return radius * sum;
				};
				auto refine = [&](auto&& self, S low, S high, S whole, S budget, int depth) -> S {
					const S middle = (low + high) / 2;
					const S left = rule(low, middle), right = rule(middle, high), sum = left + right;
					// Subdivision halves the truncation budget but not the relative
					// rounding error of the positive Gauss sums. Account for both
					// rules' arithmetic locally instead of refining below that floor.
					if (std::abs(sum - whole) <= budget + PrecisionValues<S>::RoundoffZeroThreshold * sum)
						return sum;
					if (depth == Precision::ForwardShapeQuadratureMaxDepth)
						throw ConvergenceError("SpecialFunctions shape quadrature exhausted its subdivision budget",
							evaluations, static_cast<Real>(std::abs(sum - whole)));
					return self(self, low, middle, left, budget / 2, depth + 1)
						+ self(self, middle, high, right, budget / 2, depth + 1);
				};
				const S initial = rule(0, endpoint);
				const S result = refine(refine, 0, endpoint, initial,
					PrecisionValues<S>::RoundoffZeroThreshold * initial, 0);
				if (!(std::isfinite(result) && result > 0))
					throw NumericalMethodError("SpecialFunctions shape quadrature integral is invalid");
				return result;
			}

			inline ShapeReal ShapeIntegralWindow()
			{
				return -Precision::ForwardShapeQuadratureTailLogs * std::log(PrecisionValues<ShapeReal>::Eps);
			}
		}

		/// @brief Gamma function Γ(x).
		inline Real Gamma(Real x) { return std::tgamma(x); }

		/// @brief Natural log of |Γ(x)|.
		inline Real LnGamma(Real x) { return std::lgamma(x); }

		/// @brief Beta function B(a,b) = Γ(a)Γ(b)/Γ(a+b).
		inline Real Beta(Real a, Real b) {
			if (std::isfinite(a) && a > 0 && std::isfinite(b) && b > 0)
				return static_cast<Real>(std::exp(SpecialFunctionsDetail::PositiveLogBeta(a, b)));
			return std::exp(std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b));
		}

		/// @brief Natural log of the Beta function.
		inline Real LnBeta(Real a, Real b) {
			if (std::isfinite(a) && a > 0 && std::isfinite(b) && b > 0)
				return static_cast<Real>(SpecialFunctionsDetail::PositiveLogBeta(a, b));
			return std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
		}

		///////////////////////////////////////////////////////////////////////
		///                  Incomplete gamma family                        ///
		///////////////////////////////////////////////////////////////////////

		namespace SpecialFunctionsDetail
		{
			// Numerical Recipes modified Lentz fraction for the upper gamma, without its exponential prefactor.
			inline ShapeReal GammaUpperFraction(ShapeReal a, ShapeReal x)
			{
				const ShapeReal fpmin = std::numeric_limits<ShapeReal>::min() / PrecisionValues<ShapeReal>::Eps;
				ShapeReal b = x + 1 - a;
				ShapeReal c = 1 / fpmin;
				ShapeReal d = 1 / b;
				ShapeReal h = d;
				ShapeReal residual = 0;
				for (int n = 1; n <= Precision::IncompleteGammaMaxIterations; ++n)
				{
					ShapeReal an = -ShapeReal(n) * (ShapeReal(n) - a);
					b += 2;
					d = an * d + b;
					if (std::abs(d) < fpmin) d = std::copysign(fpmin, d);
					c = b + an / c;
					if (std::abs(c) < fpmin) c = std::copysign(fpmin, c);
					d = 1 / d;
					ShapeReal delta = d * c;
					h *= delta;
					if (!std::isfinite(h) || h <= 0)
						throw NumericalMethodError("SpecialFunctions gamma continued fraction became non-finite or non-positive");
					residual = std::abs(delta - 1);
					if (residual <= PrecisionValues<ShapeReal>::Eps)
						return h;
				}
				throw ConvergenceError("SpecialFunctions gamma continued fraction exhausted its iteration budget",
					Precision::IncompleteGammaMaxIterations, static_cast<Real>(residual));
			}

			inline void ValidateGammaArguments(Real a, Real x)
			{
				if (!std::isfinite(a) || a <= 0)
					throw DomainError("SpecialFunctions incomplete gamma - a must be finite and > 0");
				if (std::isnan(x) || x < 0)
					throw DomainError("SpecialFunctions incomplete gamma - x must be >= 0");
			}

			inline ShapeReal GammaLowerSeriesLog(ShapeReal a, ShapeReal x)
			{
				// Factor 1/a out of the Numerical Recipes lower series before summation.
				ShapeReal sum = 1;
				ShapeReal delta = sum;
				for (int n = 1; n <= Precision::IncompleteGammaMaxIterations; ++n)
				{
					delta *= x / (a + n);
					sum += delta;
					if (!std::isfinite(sum))
						throw NumericalMethodError("SpecialFunctions gamma series became non-finite");
					if (std::abs(delta) <= PrecisionValues<ShapeReal>::Eps * std::abs(sum))
						return std::log(sum);
				}
				throw ConvergenceError("SpecialFunctions gamma series exhausted its iteration budget",
					Precision::IncompleteGammaMaxIterations, static_cast<Real>(std::abs(delta / sum)));
			}

			inline ShapeReal GammaLogScale(ShapeReal a, ShapeReal x)
			{
				using S = ShapeReal;
				if (a < Precision::ForwardShapeStirlingMin)
					return a * std::log(x) - x - a * GammaLogOnePlusRatio(a) + std::log(a);
				const S relative = (x - a) / a;
				const S deviance = relative > -1 && std::isfinite(relative)
					? a * LogOnePlusMinusX(relative) : a * (std::log(x) - std::log(a)) + a - x;
				return deviance + (std::log(a) - std::log(S(2) * Constants::PI)) / 2 - StirlingCorrection(a);
			}

			inline ShapeReal GammaSmallLowerCoefficient(ShapeReal a, ShapeReal x)
			{
				// gamma(a,x)=x^a*(1/a + sum((-x)^n/(n!*(a+n)))); retain log(P)/a before multiplying by a.
				using S = ShapeReal;
				S power = -x, sum = power / (a + 1);
				for (int n = 2; n <= Precision::IncompleteGammaMaxIterations; ++n) {
					power *= -x / n;
					const S next = sum + power / (a + n);
					if (next == sum) {
						const S z = a * sum;
						return std::log(x) - GammaLogOnePlusRatio(a)
							+ sum * LogOnePlusOverX(z);
					}
					sum = next;
				}
				throw ConvergenceError("SpecialFunctions small-shape gamma series exhausted its budget",
					Precision::IncompleteGammaMaxIterations, static_cast<Real>(std::abs(power)));
			}

			inline ShapeReal GammaQuadratureLog(ShapeReal a, ShapeReal x, bool lower, bool normalized)
			{
				using S = ShapeReal;
				const S difference = (a - x) - 1;
				const S peakWidth = std::sqrt(a) / x;
				const S relativeWidth = difference == 0 ? peakWidth : std::min(peakWidth, 1 / std::abs(difference));
				const S direction = lower ? -S(1) : S(1);
				const S endpoint = lower ? std::min(ShapeIntegralWindow(), 1 / relativeWidth) : ShapeIntegralWindow();
				const S integral = PositiveShapeIntegral([&](S y) {
					const S r = direction * relativeWidth * y;
					return std::exp((a - 1) * LogOnePlusMinusX(r) + difference * r);
				}, endpoint);
				const S scale = normalized ? GammaLogScale(a, x) : a * std::log(x) - x;
				return scale + std::log(relativeWidth) + std::log(integral);
			}

			inline Real GammaDirectTailLog(Real a, Real x, Real logScale)
			{
				if (a >= Precision::ForwardShapeQuadratureSwitch)
					return static_cast<Real>(GammaQuadratureLog(a, x, x < a + 1, false));
				if (x < a + 1)
					return static_cast<Real>((GammaLowerSeriesLog(a, x) - std::log(a)) + logScale);
				return static_cast<Real>(std::log(GammaUpperFraction(a, x)) + logScale);
			}

			inline ShapeReal CheckedShapeLogProbability(ShapeReal logTail)
			{
				if (std::isnan(logTail) || logTail == std::numeric_limits<Real>::infinity())
					throw NumericalMethodError("SpecialFunctions log probability is invalid");
				if (logTail > 0)
				{
					if (logTail > Precision::RoundoffZeroThreshold)
						throw NumericalMethodError("SpecialFunctions probability exceeds one");
					logTail = 0;
				}
				return logTail;
			}

			inline Real CheckedGammaLogProbability(Real logTail)
			{
				return static_cast<Real>(CheckedShapeLogProbability(logTail));
			}

			inline ShapeReal LogScaledComplement(ShapeReal factor, ShapeReal coefficient)
			{
				const ShapeReal value = factor * coefficient;
				if (value < -Precision::ForwardShapeLogSeriesLimit)
					return std::log(-std::expm1(value));
				return CheckedShapeLogProbability(std::log(factor) + std::log(-coefficient) + std::log(Expm1OverX(value)));
			}

			inline ShapeReal GammaLogTail(Real a, Real x, bool upper)
			{
				if (a <= 1 && x <= 1) {
					const ShapeReal coefficient = GammaSmallLowerCoefficient(a, x), z = ShapeReal(a) * coefficient;
					return upper ? LogScaledComplement(a, coefficient) : z;
				}
				const bool lower = x < a + 1;
				ShapeReal direct;
				if (a >= Precision::ForwardShapeQuadratureSwitch)
					direct = GammaQuadratureLog(a, x, lower, true);
				else
					direct = GammaLogScale(a, x) + (lower
						? ShapeReal(GammaLowerSeriesLog(a, x)) - std::log(ShapeReal(a))
						: std::log(ShapeReal(GammaUpperFraction(a, x))));
				const ShapeReal logTail = CheckedShapeLogProbability(direct);
				return upper == lower ? std::log(-std::expm1(logTail)) : logTail;
			}

			inline Real RegularizedGamma(Real a, Real x, bool upper)
			{
				ValidateGammaArguments(a, x);
				if (x == 0) return upper ? Real(1) : Real(0);
				if (std::isinf(x)) return upper ? Real(0) : Real(1);

				return static_cast<Real>(std::exp(GammaLogTail(a, x, upper)));
			}

			inline Real IncompleteGamma(Real a, Real x, bool upper)
			{
				ValidateGammaArguments(a, x);
				if (x == 0) return upper ? std::tgamma(a) : Real(0);
				if (std::isinf(x)) return upper ? Real(0) : std::tgamma(a);
				if (a <= 1 && x <= 1 && upper) {
					const ShapeReal coefficient = GammaSmallLowerCoefficient(a, x);
					const ShapeReal z = ShapeReal(a) * coefficient;
					const ShapeReal ratio = Expm1OverX(z);
					return static_cast<Real>(std::exp(ShapeReal(a) * GammaLogOnePlusRatio(a)) * -coefficient * ratio);
				}

				// Keep normalization and complementary-tail subtraction in log-
				// probability units. Subtracting two rounded O(a log a) logs
				// destroys a central complement in float.
				if (a >= Precision::ForwardShapeQuadratureSwitch) {
					const bool lower = x < a + 1;
					if (upper != lower)
						return static_cast<Real>(std::exp(GammaQuadratureLog(a, x, lower, false)));
					const ShapeReal logProbability = GammaLogTail(a, x, upper);
					const ShapeReal logGamma = std::lgamma(ShapeReal(a));
					return static_cast<Real>(std::exp(logGamma + logProbability));
				}
				const Real logTail = GammaDirectTailLog(a, x, a * std::log(x) - x);
				if (std::isnan(logTail))
					throw NumericalMethodError("SpecialFunctions incomplete gamma log tail is NaN");
				if (upper != (x < a + 1))
					return std::exp(logTail);

				const Real logGamma = std::lgamma(a);
				if (logGamma == std::numeric_limits<Real>::infinity())
					return logGamma;
				if (logTail == -std::numeric_limits<Real>::infinity())
					return std::exp(logGamma);
				const Real logProbability = CheckedGammaLogProbability(logTail - logGamma);
				if (logProbability == 0)
					throw NumericalMethodError("SpecialFunctions incomplete gamma complement lost to roundoff");
				return std::exp(logGamma + std::log(-std::expm1(logProbability)));
			}
		} // namespace SpecialFunctionsDetail

		/// @brief Regularized lower incomplete gamma P(a,x) = γ(a,x)/Γ(a).
		/// @details Checked series/fraction kernels and mode-scaled density quadrature for large shapes.
		/// @throws DomainError for invalid parameters; ConvergenceError on exhausted iteration budget.
		inline Real RegularizedGammaP(Real a, Real x) { return SpecialFunctionsDetail::RegularizedGamma(a, x, false); }

		/// @brief Regularized upper incomplete gamma Q(a,x) = Γ(a,x)/Γ(a) = 1 - P(a,x).
		/// @details Direct upper tails and factored small-shape complements preserve representable probabilities.
		inline Real RegularizedGammaQ(Real a, Real x) { return SpecialFunctionsDetail::RegularizedGamma(a, x, true); }

		/// @brief Lower incomplete gamma γ(a,x) = P(a,x)·Γ(a).
		/// @details Log-domain scaling avoids forming P or Γ separately. Requires finite a > 0 and x >= 0.
		///          Accepts x = +infinity; genuine final overflow/underflow returns +infinity/rounded zero.
		/// @throws DomainError for invalid parameters; ConvergenceError on exhausted iteration budget;
		///         NumericalMethodError on arithmetic failure or a complement lost to roundoff.
		inline Real IncompleteGammaLower(Real a, Real x) { return SpecialFunctionsDetail::IncompleteGamma(a, x, false); }

		/// @brief Upper incomplete gamma Γ(a,x) = Q(a,x)·Γ(a).
		/// @details Uses the same domain, range and error contract as IncompleteGammaLower.
		///          For x >= a+1, scales the upper continued fraction without forming Q or Γ separately.
		inline Real IncompleteGammaUpper(Real a, Real x) { return SpecialFunctionsDetail::IncompleteGamma(a, x, true); }

		///////////////////////////////////////////////////////////////////////
		///                  Incomplete beta family                         ///
		///////////////////////////////////////////////////////////////////////

		namespace SpecialFunctionsDetail
		{
			// Numerical Recipes modified Lentz fraction, without the exponential prefactor.
			inline Real BetaFraction(Real x, Real a, Real b)
			{
				const Real qab = a + b;
				const Real fpmin = std::numeric_limits<Real>::min() / std::numeric_limits<Real>::epsilon();
				const Real qap = a + MML::Real(1.0);
				const Real qam = a - MML::Real(1.0);
				Real c = 1.0;
				Real d = MML::Real(1.0) - qab * x / qap;
				if (std::abs(d) < fpmin) d = fpmin;
				d = MML::Real(1.0) / d;
				Real h = d;
				Real residual = 0;
				for (int m = 1; m <= Precision::IncompleteBetaMaxIterations; m++)
				{
					int m2 = 2 * m;
					Real aa = (Real(m) / (a + m2)) * ((b - m) / (qam + m2)) * x;
					d = MML::Real(1.0) + aa * d;
					if (std::abs(d) < fpmin) d = fpmin;
					c = MML::Real(1.0) + aa / c;
					if (std::abs(c) < fpmin) c = fpmin;
					d = MML::Real(1.0) / d;
					h *= d * c;

					aa = -((a + m) / (a + m2)) * (((qab + m) / (qap + m2)) * x);
					d = MML::Real(1.0) + aa * d;
					if (std::abs(d) < fpmin) d = fpmin;
					c = MML::Real(1.0) + aa / c;
					if (std::abs(c) < fpmin) c = fpmin;
					d = MML::Real(1.0) / d;
					Real del = d * c;
					h *= del;
					if (!std::isfinite(h) || h <= 0)
						throw NumericalMethodError("SpecialFunctions incomplete beta fraction became non-finite or non-positive");
					residual = std::abs(del - 1);
					if (residual <= Precision::Eps)
						return h;
				}
				throw ConvergenceError("SpecialFunctions incomplete beta exhausted its iteration budget",
					Precision::IncompleteBetaMaxIterations, residual);
			}

			inline ShapeReal SmallBetaNormalizationRatio(ShapeReal a, ShapeReal b)
			{
				// log(a*B(a,b))/a: shift b with beta's recurrence, then difference Stirling terms without subtracting them.
				using S = ShapeReal;
				S recurrence = 0;
				while (b < Precision::ForwardShapeStirlingMin) {
					recurrence += LogOnePlusOverX(a / b) / b;
					b += 1;
				}
				const S ratio = a / b, logRatio = std::log1p(ratio);
				const S inv = 1 / b;
				S power = inv, difference = 0;
				for (int i = 0; i < 8; ++i) {
					const S exponent = S(2 * i + 1);
					difference += StirlingCoefficients[i] * power * exponent * inv * LogOnePlusOverX(ratio)
						* Expm1OverX(-exponent * logRatio);
					power *= inv * inv;
				}
				const S devianceRatio = ratio == 0 ? S(0) : LogOnePlusMinusX(ratio) / ratio;
				return GammaLogOnePlusRatio(a) - std::log(b) - (1 + (a - S(0.5)) / b) * devianceRatio
					- (a - S(0.5)) / b + difference + recurrence;
			}

			inline ShapeReal SmallBetaLogTail(ShapeReal x, ShapeReal logX, ShapeReal a, ShapeReal b, bool upper)
			{
				using S = ShapeReal;
				S term = (1 - b) * x, sum = term / (a + 1);
				for (int n = 2; n <= Precision::IncompleteBetaMaxIterations; ++n) {
					term *= (S(n) - b) * (x / n);
					const S next = sum + term / (a + n);
					if (next == sum) {
						if (b >= Precision::GammaInverseSmallShapeLimit) {
							const S coefficient = logX - SmallBetaNormalizationRatio(a, b) + sum * LogOnePlusOverX(a * sum);
							if (!(coefficient < 0))
								throw NumericalMethodError("SpecialFunctions small-shape beta tail is invalid");
							const S value = a * coefficient;
							return upper ? LogScaledComplement(a, coefficient) : value;
						}
						const S ratio = a <= b ? std::log1p(a / b) : std::log(a) - std::log(b) + std::log1p(b / a);
						const S normalization = a * GammaLogOnePlusRatio(a) + b * GammaLogOnePlusRatio(b)
							- (a + b) * GammaLogOnePlusRatio(a + b) + ratio;
						const S value = CheckedShapeLogProbability(a * logX - normalization + std::log1p(a * sum));
						return upper ? std::log(-std::expm1(value)) : value;
					}
					sum = next;
				}
				throw ConvergenceError("SpecialFunctions small-shape beta series exhausted its budget",
					Precision::IncompleteBetaMaxIterations, static_cast<Real>(std::abs(term)));
			}

			inline ShapeReal BetaLogFrontOverA(ShapeReal x, ShapeReal logX, ShapeReal logY, ShapeReal a, ShapeReal b)
			{
				using S = ShapeReal;
				const S sum = a + b;
				if (a <= Precision::GammaInverseSmallShapeLimit && b >= Precision::GammaInverseSmallShapeLimit)
					return a * (logX - SmallBetaNormalizationRatio(a, b)) + b * logY;
				if (a >= Precision::ForwardShapeStirlingMin && b >= Precision::ForwardShapeStirlingMin) {
					const S smaller = std::min(a, b), larger = std::max(a, b);
					const S sumError = (larger - sum) + smaller;
					const S numerator = std::fma(sum, x, -a) + sumError * x;
					const S deviance = a * LogOnePlusMinusX(numerator / a) + b * LogOnePlusMinusX(-numerator / b);
					return deviance + (std::log(smaller) - std::log1p(smaller / larger)
						- std::log(S(2) * Constants::PI)) / 2
						- StirlingCorrection(a) - StirlingCorrection(b) + StirlingCorrection(sum) - std::log(a);
				}
				if (a >= Precision::ForwardShapeStirlingMin)
					return BetaLogFrontOverA(1 - x, logY, logX, b, a) + std::log(b) - std::log(a);
				if (b >= Precision::ForwardShapeStirlingMin) {
					const S ratio = a / b;
					const S logProduct = x > 0 ? std::log(b * x) : logX + std::log(b);
					return a * logProduct + b * logY - a * GammaLogOnePlusRatio(a)
						+ (b + a - S(0.5)) * LogOnePlusMinusX(ratio) + (a - S(0.5)) * ratio
						- StirlingCorrection(b) + StirlingCorrection(sum);
				}
				if (sum <= Precision::GammaInverseSmallShapeLimit) {
					const S logRatio = a <= b ? std::log1p(a / b) : std::log(a) - std::log(b) + std::log1p(b / a);
					const S normalization = a * GammaLogOnePlusRatio(a) + b * GammaLogOnePlusRatio(b)
						- sum * GammaLogOnePlusRatio(sum) + logRatio;
					return a * logX + b * logY - normalization;
				}
				return a * logX + b * logY - PositiveLogBeta(a, b) - std::log(a);
			}

			inline ShapeReal BetaQuadratureLog(ShapeReal x, ShapeReal a, ShapeReal b, bool lower)
			{
				using S = ShapeReal;
				const S sum = a + b, y = 1 - x;
				const S sumError = (std::max(a, b) - sum) + std::min(a, b);
				const S numerator = std::fma(-sum, x, a) - sumError * x + (2 * x - 1);
				const S peakWidth = std::sqrt(y) / std::sqrt(x) / std::sqrt(sum - 2);
				const S relativeX = numerator == 0 ? peakWidth : std::min(peakWidth, y / std::abs(numerator));
				const S relativeY = relativeX * x / y;
				const S direction = lower ? -S(1) : S(1);
				const S relativeWidth = lower ? relativeX : relativeY;
				const S endpoint = relativeWidth == 0 ? ShapeIntegralWindow() : std::min(ShapeIntegralWindow(), 1 / relativeWidth);
				const S integral = PositiveShapeIntegral([&](S t) {
					const S dx = direction * relativeX * t, dy = -direction * relativeY * t;
					return std::exp((a - 1) * LogOnePlusMinusX(dx) + (b - 1) * LogOnePlusMinusX(dy)
						+ direction * (numerator * relativeX / y) * t);
				}, endpoint);
				const S logX = std::log(x), logY = std::log1p(-x);
				return BetaLogFrontOverA(x, logX, logY, a, b) + std::log(a)
					+ std::log(relativeX) - logY + std::log(integral);
			}

			inline ShapeReal BetaLogTail(Real x, ShapeReal logX, Real a, Real b, bool upper)
			{
				using S = ShapeReal;
				if (!std::isfinite(a + b))
					throw NumericalMethodError("SpecialFunctions incomplete beta shape sum is non-finite");
				const S sum = S(a) + b, logY = std::log1p(-S(x));
				if (a == b && x == S(0.5)) return -std::log(S(2));
				if (b == 1) {
					const S logP = S(a) * logX;
					if (!upper || std::isinf(logP)) return upper ? S(0) : logP;
					return LogScaledComplement(a, logX);
				}
				if (a == 1) {
					const S logQ = S(b) * logY;
					if (upper || std::isinf(logQ)) return upper ? logQ : S(0);
					return LogScaledComplement(b, logY);
				}
				if (a <= Precision::GammaInverseSmallShapeLimit && x <= S(0.5) && S(b) * x <= 1)
					return SmallBetaLogTail(x, logX, a, b, upper);
				if (b <= Precision::GammaInverseSmallShapeLimit && x >= S(0.5) && S(a) * (1 - S(x)) <= 1)
					return SmallBetaLogTail(1 - S(x), logY, b, a, !upper);
				const bool lower = x <= (S(a) + 1) / (sum + 2);
				S direct;
				if (std::min(a, b) >= Precision::ForwardShapeStirlingMin && sum >= Precision::ForwardShapeQuadratureSwitch && x > 0)
					direct = BetaQuadratureLog(x, a, b, lower);
				else if (lower)
					direct = BetaLogFrontOverA(x, logX, logY, a, b) + std::log(S(BetaFraction(x, a, b)));
				else
					direct = BetaLogFrontOverA(1 - S(x), logY, logX, b, a) + std::log(S(BetaFraction(1 - x, b, a)));
				if (upper == lower && direct >= -Precision::SqrtEps) {
					const S opposite = lower
						? BetaLogFrontOverA(1 - S(x), logY, logX, b, a) + std::log(S(BetaFraction(1 - x, b, a)))
						: BetaLogFrontOverA(x, logX, logY, a, b) + std::log(S(BetaFraction(x, a, b)));
					return CheckedShapeLogProbability(opposite);
				}
				const S logTail = CheckedShapeLogProbability(direct);
				return upper == lower ? std::log(-std::expm1(logTail)) : logTail;
			}

			inline void ValidateBetaArguments(Real x, Real a, Real b)
			{
				if (!(std::isfinite(a) && a > 0 && std::isfinite(b) && b > 0))
					throw DomainError("SpecialFunctions incomplete beta - a,b must be finite and > 0");
				if (!(x >= 0 && x <= 1))
					throw DomainError("SpecialFunctions incomplete beta - x must be in [0,1]");
			}
		} // namespace SpecialFunctionsDetail

		/// @brief Regularized incomplete beta I_x(a,b) = B(x;a,b)/B(a,b).
		/// @details Stable normalization, small-shape series, ratio-ordered Lentz fraction, and large-shape density quadrature.
		/// @throws DomainError for non-finite/non-positive shapes or x outside [0,1].
		/// @throws ConvergenceError if the continued-fraction budget is exhausted.
		inline Real RegularizedBetaI(Real x, Real a, Real b)
		{
			SpecialFunctionsDetail::ValidateBetaArguments(x, a, b);
			if (x == 0.0) return 0.0;
			if (x == 1.0) return 1.0;

			return static_cast<Real>(std::exp(SpecialFunctionsDetail::BetaLogTail(x, std::log(SpecialFunctionsDetail::ShapeReal(x)), a, b, false)));
		}

		/// @brief Incomplete beta B(x;a,b) = I_x(a,b)·B(a,b).
		/// @details Log scaling avoids 0*infinity; domains and convergence errors match RegularizedBetaI.
		inline Real IncompleteBeta(Real x, Real a, Real b) {
			SpecialFunctionsDetail::ValidateBetaArguments(x, a, b);
			if (x == 0) return 0;
			if (x == 1) return Beta(a, b);
			const auto logTail = SpecialFunctionsDetail::BetaLogTail(x, std::log(SpecialFunctionsDetail::ShapeReal(x)), a, b, false);
			return static_cast<Real>(std::exp(SpecialFunctionsDetail::PositiveLogBeta(a, b) + logTail));
		}

		///////////////////////////////////////////////////////////////////////
		///                     Polygamma family                            ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Digamma ψ(x) = d/dx ln Γ(x).
		/// @details Reflection for x < 0.5, recurrence up to x >= 6, then asymptotic Bernoulli series.
		/// @throws DomainError at non-positive integers (poles).
		inline Real Digamma(Real x)
		{
			if (x <= 0.0 && x == std::floor(x)) throw DomainError("SpecialFunctions::Digamma - pole at non-positive integer");

			Real result = 0.0;
			if (x < 0.5)
			{
				result -= Constants::PI / std::tan(Constants::PI * x);
				x = MML::Real(1.0) - x;
			}
			while (x < MML::Real(6.0)) { result -= MML::Real(1.0) / x; x += MML::Real(1.0); }

			const Real inv = MML::Real(1.0) / x;
			const Real inv2 = inv * inv;
			result += std::log(x) - MML::Real(0.5) * inv;
			result -= inv2 * (static_cast<MML::Real>(1.0) / 12 - inv2 * (static_cast<MML::Real>(1.0) / 120 - inv2 * (static_cast<MML::Real>(1.0) / 252
				- inv2 * (static_cast<MML::Real>(1.0) / 240 - inv2 * (static_cast<MML::Real>(1.0) / 132)))));
			return result;
		}

		/// @brief Trigamma ψ'(x).
		/// @details Reflection for x < 0.5, recurrence up to x >= 6, then asymptotic series.
		inline Real Trigamma(Real x)
		{
			if (x < 0.5)
			{
				const Real s = std::sin(Constants::PI * x);
				return (Constants::PI * Constants::PI) / (s * s) - Trigamma(MML::Real(1.0) - x);
			}
			Real result = 0.0;
			while (x < MML::Real(6.0)) { result += MML::Real(1.0) / (x * x); x += MML::Real(1.0); }

			const Real inv = MML::Real(1.0) / x;
			const Real inv2 = inv * inv;
			result += inv + MML::Real(0.5) * inv2 + inv * inv2 * (static_cast<MML::Real>(1.0) / 6 - inv2 * (static_cast<MML::Real>(1.0) / 30
				- inv2 * (static_cast<MML::Real>(1.0) / 42 - inv2 * (static_cast<MML::Real>(1.0) / 30))));
			return result;
		}

		/// @brief Polygamma ψ^(n)(x), the n-th derivative of the digamma function.
		/// @details n==0 → Digamma, n==1 → Trigamma; for n≥1 uses
		///          ψ^(n)(x) = (-1)^(n+1) n! ζ(n+1, x) with Euler–Maclaurin on the Hurwitz zeta.
		/// @throws DomainError if n < 0 or x <= 0.
		inline Real Polygamma(int n, Real x)
		{
			if (n < 0) throw DomainError("SpecialFunctions::Polygamma - n must be >= 0");
			if (n == 0) return Digamma(x);
			if (n == 1) return Trigamma(x);
			if (x <= 0.0) throw DomainError("SpecialFunctions::Polygamma - x must be > 0 for n >= 2");

			const int s = n + 1;
			Real nfact = 1.0;
			for (int k = 2; k <= n; k++) nfact *= k;
			const Real sign = (n % 2 == 0) ? -MML::Real(1.0) : MML::Real(1.0);   // (-1)^(n+1)

			// Sum (x+k)^{-s} while shifting x up to a safe base for the asymptotic tail.
			Real sum = 0.0;
			const Real target = MML::Real(10.0) + n;
			while (x < target) { sum += std::pow(x, -(Real)s); x += 1.0; }

			// Euler–Maclaurin tail of the Hurwitz zeta ζ(s, x).
			Real tail = std::pow(x, MML::Real(1.0) - s) / (s - 1) + MML::Real(0.5) * std::pow(x, -(Real)s);
			static const Real emC[] = { static_cast<MML::Real>(1.0) / 12, static_cast<MML::Real>(-1.0) / 720, static_cast<MML::Real>(1.0) / 30240,
										static_cast<MML::Real>(-1.0) / 1209600, static_cast<MML::Real>(1.0) / 47900160 };
			Real poch = Real(s);                    // (s)_{2j-1}, starts at (s)_1 = s
			Real xpow = std::pow(x, -(Real)(s + 1)); // x^{-(s+2j-1)}, starts at x^{-(s+1)}
			const Real xinv2 = MML::Real(1.0) / (x * x);
			for (int j = 1; j <= 5; j++)
			{
				tail += emC[j - 1] * poch * xpow;
				poch *= (Real)(s + 2 * j - 1) * (Real)(s + 2 * j);
				xpow *= xinv2;
			}

			return sign * nfact * (sum + tail);
		}

		///////////////////////////////////////////////////////////////////////
		///                  Error-function inverses                        ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Inverse error function erf⁻¹(x), x in (-1,1).
		/// @details Giles (2010) rational initial guess refined by two Halley iterations against erf.
		inline Real ErfInv(Real x)
		{
			if (x <= -1.0) { if (x == -1.0) return -std::numeric_limits<Real>::infinity(); throw DomainError("SpecialFunctions::ErfInv - |x| must be <= 1"); }   // tol-ok: exact special value at the domain edge
			if (x >= 1.0) { if (x == 1.0) return std::numeric_limits<Real>::infinity(); throw DomainError("SpecialFunctions::ErfInv - |x| must be <= 1"); }
			if (x == 0.0) return 0.0;

			Real w = -std::log((MML::Real(1.0) - x) * (MML::Real(1.0) + x));
			Real p;
			if (w < 5.0)
			{
				w -= 2.5;
				p = static_cast<MML::Real>(2.81022636e-08);   // tol-ok: Giles erfinv coefficients
				p = static_cast<MML::Real>(3.43273939e-07) + p * w;
				p = static_cast<MML::Real>(-3.5233877e-06) + p * w;
				p = static_cast<MML::Real>(-4.39150654e-06) + p * w;
				p = static_cast<MML::Real>(0.00021858087) + p * w;
				p = static_cast<MML::Real>(-0.00125372503) + p * w;
				p = static_cast<MML::Real>(-0.00417768164) + p * w;
				p = static_cast<MML::Real>(0.246640727) + p * w;
				p = static_cast<MML::Real>(1.50140941) + p * w;
			}
			else
			{
				w = std::sqrt(w) - MML::Real(3.0);
				p = static_cast<MML::Real>(-0.000200214257);
				p = static_cast<MML::Real>(0.000100950558) + p * w;
				p = static_cast<MML::Real>(0.00134934322) + p * w;
				p = static_cast<MML::Real>(-0.00367342844) + p * w;
				p = static_cast<MML::Real>(0.00573950773) + p * w;
				p = static_cast<MML::Real>(-0.0076224613) + p * w;
				p = static_cast<MML::Real>(0.00943887047) + p * w;
				p = static_cast<MML::Real>(1.00167406) + p * w;
				p = static_cast<MML::Real>(2.83297682) + p * w;
			}
			Real res = p * x;

			// Halley refinement against erf (2 iterations reach full double precision).
			const Real twoOverSqrtPi = static_cast<MML::Real>(2.0) / std::sqrt(Constants::PI);
			for (int i = 0; i < 2; i++)
			{
				Real err = std::erf(res) - x;
				Real deriv = twoOverSqrtPi * std::exp(-res * res);
				res -= err / (deriv - res * err);   // Halley step
			}
			return res;
		}

		/// @brief Inverse complementary error function erfc⁻¹(y), y in (0,2).
		inline Real ErfcInv(Real y)
		{
			if (std::isnan(y))
				throw DomainError("SpecialFunctions::ErfcInv - y must be in [0,2]");
			if (y <= 0.0) { if (y == 0.0) return std::numeric_limits<Real>::infinity(); throw DomainError("SpecialFunctions::ErfcInv - y must be in (0,2)"); }
			if (y >= 2.0) { if (y == 2.0) return -std::numeric_limits<Real>::infinity(); throw DomainError("SpecialFunctions::ErfcInv - y must be in (0,2)"); }
			if (y > 1) return -ErfcInv(2 - y);
			if (y >= Real(0.5)) return ErfInv(1 - y);

			const Real logTarget = std::log(y);
			const Real sqrtPi = std::sqrt(Constants::PI);
			Real result = std::sqrt(-logTarget);
			Real step = 0;
			for (int i = 0; i < Precision::ErfcInverseMaxIterations; ++i)
			{
				Real squared = result * result;
				Real scaledErfc, logTail;
				Real tail = std::erfc(result);
				if (std::isnormal(tail))
				{
					scaledErfc = tail * std::exp(squared);
					logTail = std::log(tail);
				}
				else
				{
					// erfc(x) = Q(1/2,x^2): reuse its fraction when the library tail loses precision to underflow.
					scaledErfc = static_cast<Real>(
						SpecialFunctionsDetail::GammaUpperFraction(Real(0.5), squared) * result / sqrtPi);
					logTail = std::log(scaledErfc) - squared;
				}
				step = (logTail - logTarget) * (sqrtPi / 2) * scaledErfc;
				Real next = result + step;
				if (!std::isfinite(next) || next < 0)
					throw NumericalMethodError("SpecialFunctions::ErfcInv - non-finite or negative tail iterate");
				if (std::abs(step) <= Precision::Eps * std::max(Real(1), next))
					return next;
				result = next;
			}
			throw ConvergenceError("SpecialFunctions::ErfcInv - exhausted log-tail refinement budget",
				Precision::ErfcInverseMaxIterations, std::abs(step));
		}

		///////////////////////////////////////////////////////////////////////
		///             Inverses of the regularized functions               ///
		///////////////////////////////////////////////////////////////////////

		namespace SpecialFunctionsDetail
		{
			inline Real GammaLogOnePlus(Real a)
			{
				return static_cast<Real>(ShapeReal(a) * GammaLogOnePlusRatio(a));
			}

			inline Real GammaPInverse(Real a, Real p, int maxIterations)
			{
				if (!(std::isfinite(a) && a > 0))
					throw DomainError("SpecialFunctions::RegularizedGammaPInv - a must be finite and > 0");
				if (!(std::isfinite(p) && p >= 0 && p <= 1))
					throw DomainError("SpecialFunctions::RegularizedGammaPInv - p must be finite and in [0,1]");
				if (p == 0) return 0;
				if (p == 1) return std::numeric_limits<Real>::infinity();
				if (a == 1) return -std::log1p(-p);
				if (a == Real(0.5)) {
					const Real root = p <= Real(0.5) ? ErfInv(p) : ErfcInv(1 - p);
					return root * root;
				}

				const Real logGamma = std::lgamma(a);
				const Real logGammaOnePlus = GammaLogOnePlus(a);
				if (!(std::isfinite(logGamma) && std::isfinite(logGammaOnePlus)))
					throw NumericalMethodError("SpecialFunctions gamma inverse normalization is non-finite");
				const Real initial = std::exp((std::log(p) + logGammaOnePlus) / a);
				// P(a,x) ~ x^a/Gamma(a+1); the inverse's first relative correction is x/(a+1).
				// Refining below this scale only adds quantization error from the log-tail arithmetic.
				if (initial <= std::numeric_limits<Real>::min() || initial / (a + 1) <= Precision::Eps) // tol-ok: IEEE range boundary
					return initial;

				const bool upper = p > Real(0.5);
				const Real logTarget = upper ? std::log1p(-p) : std::log(p);
				int iterations = 0;
				Real residual = std::numeric_limits<Real>::infinity();
				auto evaluate = [&](Real x) {
					if (iterations >= maxIterations)
						throw ConvergenceError("SpecialFunctions gamma inverse exhausted its iteration budget",
							iterations, residual);
					++iterations;
					return static_cast<Real>(GammaLogTail(a, x, upper));
				};
				auto belowRoot = [&](Real logTail) { return upper ? logTail > logTarget : logTail < logTarget; };

				Real low = 0, high = std::max(a, Real(1));
				Real lowLog = upper ? Real(0) : -std::numeric_limits<Real>::infinity();
				Real highLog = evaluate(high);
				residual = std::abs(highLog - logTarget);
				while (belowRoot(highLog)) {
					low = high;
					lowLog = highLog;
					high = high > std::numeric_limits<Real>::max() / 2
						? std::numeric_limits<Real>::max() : 2 * high;
					if (high == low)
						throw NumericalMethodError("SpecialFunctions gamma inverse cannot bracket a finite result");
					highLog = evaluate(high);
					residual = std::abs(highLog - logTarget);
				}

				Real x = initial > low && initial < high ? initial : low + (high - low) / 2;
				const Real targetLogTolerance = Precision::RoundoffZeroThreshold * (1 + std::abs(logTarget));
				for (;;) {
					const Real logTail = evaluate(x);
					const Real error = logTail - logTarget;
					const Real logScale = static_cast<Real>(GammaLogScale(a, x));
					const Real rate = std::exp(logScale - logTail);
					const Real logTolerance = targetLogTolerance + Precision::Eps *
						(1 + std::abs(logScale) + std::abs(logTail) + 2 * rate);
					residual = std::abs(error);
					if (belowRoot(logTail)) { low = x; lowLog = logTail; }
					else { high = x; highLog = logTail; }
					if (residual <= logTolerance && high - low <= Precision::RoundoffZeroThreshold * high)
						return x;
					if (std::nextafter(low, high) == high) {
						if (std::min(std::abs(lowLog - logTarget), std::abs(highLog - logTarget)) <= logTolerance)
							return std::abs(lowLog - logTarget) < std::abs(highLog - logTarget) ? low : high;
						throw ConvergenceError("SpecialFunctions gamma inverse stagnated before reaching tail accuracy",
							iterations, residual);
					}

					// Newton on log P/log Q: x*d(log tail)/dx avoids forming an underflowed density.
					const Real correction = error / (upper ? -rate : rate);
					Real next = x - x * correction;
					if (!(std::isfinite(next) && next > low && next < high))
						next = low > 0 ? std::sqrt(low) * std::sqrt(high) : high / 2;
					if (!(next > low && next < high))
						next = low + (high - low) / 2;
					x = next;
				}
			}
		} // namespace SpecialFunctionsDetail

		/// @brief Inverse of the regularized lower incomplete gamma: solves P(a,x)=p for x.
		/// @details Bracketed log-tail Newton refinement; exact shape-1 and error-function shape-1/2 identities.
		///          Requires finite a > 0 and finite p in [0,1]. Endpoints are 0 and +infinity; tiny results round to zero.
		/// @throws DomainError for invalid inputs; ConvergenceError for exhausted inverse or forward budgets.
		inline Real RegularizedGammaPInv(Real a, Real p)
		{
			return SpecialFunctionsDetail::GammaPInverse(a, p, Precision::IncompleteGammaInverseMaxIterations);
		}

		namespace SpecialFunctionsDetail
		{
			inline Real BetaIInverse(Real p, Real a, Real b, int maxIterations)
			{
				if (!(std::isfinite(a) && a > 0 && std::isfinite(b) && b > 0))
					throw DomainError("SpecialFunctions::RegularizedBetaIInv - a,b must be finite and > 0");
				if (!(std::isfinite(p) && p >= 0 && p <= 1))
					throw DomainError("SpecialFunctions::RegularizedBetaIInv - p must be finite and in [0,1]");
				if (p == 0) return 0;
				if (p == 1) return 1;
				if (a == 1 && b == 1) return p;
				if (b == 1) {
					const Real result = std::exp(std::log(p) / a);
					// Bernoulli's inequality: this linear distance already reaches half an endpoint ULP.
					if (result == 1 && a > 1 && (1 - p) / (Precision::Eps / 4) >= a)
						return std::nextafter(Real(1), Real(0));
					return result;
				}
				if (a == 1) {
					const Real result = -std::expm1(std::log1p(-p) / b);
					// The positive quadratic correction breaks an exact half-subnormal tie upward.
					if (result == 0 && b > 1 && p / std::numeric_limits<Real>::denorm_min() == b / 2)
						return std::numeric_limits<Real>::denorm_min();
					return result;
				}

				const Real sum = a + b;
				const Real logBeta = LnBeta(a, b);
				if (!(std::isfinite(sum) && std::isfinite(logBeta)))
					throw NumericalMethodError("SpecialFunctions beta inverse normalization is non-finite");
				if (a == b && p == Real(0.5)) return Real(0.5);
				bool upper = p > Real(0.5);
				const Real q = upper ? 1 - p : p;
				const Real logTarget = std::log(q);
				int iterations = 0;
				Real residual = std::numeric_limits<Real>::infinity();
				auto evaluate = [&](Real logX) {
					if (iterations >= maxIterations)
						throw ConvergenceError("SpecialFunctions beta inverse exhausted its iteration budget",
							iterations, residual);
					++iterations;
					return static_cast<Real>(BetaLogTail(std::exp(logX), logX, a, b, upper));
				};
				auto belowRoot = [&](Real logTail) { return upper ? logTail > logTarget : logTail < logTarget; };

				// Work in log distance from the nearer support boundary; never subtract two near-one iterates.
				Real high = std::log(Real(0.5));
				const Real middleLog = evaluate(high);
				const bool reflect = belowRoot(middleLog);
				if (reflect) {
					std::swap(a, b);
					upper = !upper;
				}
				// Half an endpoint ULP: these are rounding boundaries, not artificial interior floors.
				Real low = reflect ? std::log(Precision::Eps / 4)
					: std::log(std::numeric_limits<Real>::denorm_min()) - std::log(Real(2));
				Real lowLog = evaluate(low), highLog = middleLog;
				if (!belowRoot(lowLog))
					return reflect ? Real(1) : Real(0);

				const Real logCdfTarget = upper ? std::log1p(-q) : logTarget;
				const Real initial = (logCdfTarget + std::log(a) + logBeta) / a;
				if (!upper && !reflect && initial > low && initial < high) {
					const Real distance = std::exp(initial);
					// I_x(a,b) ~ x^a/(a*B(a,b)); the inverse correction is (b-1)*x/(a+1).
					// Leave subnormal and reflected endpoint rounding to the checked bracket below.
					if (std::isnormal(distance) && std::abs(b - 1) * distance / (a + 1) <= Precision::Eps)
						return distance;
				}
				Real logX = initial > low && initial < high ? initial : high - std::log(Real(2));
				const Real targetTolerance = Precision::RoundoffZeroThreshold * (1 + std::abs(logTarget));
				auto resultAt = [&](Real position) {
					const Real distance = std::exp(position);
					return reflect ? 1 - distance : distance;
				};
				for (;;) {
					const Real logTail = evaluate(logX);
					const Real error = logTail - logTarget;
					const Real x = std::exp(logX), logY = std::log1p(-x);
					const Real logFront = static_cast<Real>(BetaLogFrontOverA(x, logX, logY, a, b)) + std::log(a);
					const Real rate = std::exp(logFront - logY - logTail);
					const Real tolerance = targetTolerance + Precision::Eps *
						(1 + std::abs(logFront) + std::abs(logTail) + 2 * rate);
					residual = std::abs(error);
					if (belowRoot(logTail)) { low = logX; lowLog = logTail; }
					else { high = logX; highLog = logTail; }
					if (residual <= tolerance && high - low <= Precision::RoundoffZeroThreshold)
						return resultAt(logX);
					if (std::nextafter(low, high) == high) {
						if (std::min(std::abs(lowLog - logTarget), std::abs(highLog - logTarget)) <= tolerance)
							return resultAt(std::abs(lowLog - logTarget) < std::abs(highLog - logTarget) ? low : high);
						throw ConvergenceError("SpecialFunctions beta inverse stagnated before reaching tail accuracy",
							iterations, residual);
					}
					Real next = logX - error / (upper ? -rate : rate);
					if (!(std::isfinite(next) && next > low && next < high))
						next = low + (high - low) / 2;
					logX = next;
				}
			}
		} // namespace SpecialFunctionsDetail

		/// @brief Inverse of the regularized incomplete beta: solves I_x(a,b)=p for x.
		/// @details Checked log-tail Newton refinement in log distance from the nearer endpoint.
		///          Requires finite positive a,b and finite p in [0,1]. Genuine range loss rounds to 0 or 1.
		/// @throws DomainError for invalid inputs; ConvergenceError for exhausted inverse or forward budgets;
		///         NumericalMethodError for invalid normalization or tail arithmetic.
		inline Real RegularizedBetaIInv(Real p, Real a, Real b)
		{
			return SpecialFunctionsDetail::BetaIInverse(p, a, b, Precision::IncompleteBetaInverseMaxIterations);
		}

		///////////////////////////////////////////////////////////////////////
		///                    Exponential integrals                        ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Exponential integral E_n(x) = ∫_1^∞ e^{-x t} / t^n dt.
		/// @details Continued fraction for x > 1, power series otherwise (Numerical Recipes).
		/// @throws DomainError for n < 0, x < 0, or x == 0 with n in {0,1}.
		inline Real ExpIntegralEn(int n, Real x)
		{
			static const Real EulerGamma = static_cast<MML::Real>(0.5772156649015328606);
			const Real eps = std::numeric_limits<Real>::epsilon();
			const Real big = std::numeric_limits<Real>::max() * eps;
			const int MAXIT = 100;
			const int nm1 = n - 1;

			if (n < 0 || x < 0.0 || (x == 0.0 && (n == 0 || n == 1)))
				throw DomainError("SpecialFunctions::ExpIntegralEn - bad arguments");

			if (n == 0) return std::exp(-x) / x;
			if (x == MML::Real(0.0)) return MML::Real(1.0) / nm1;

			if (x > 1.0)
			{
				Real b = x + n;
				Real c = big;
				Real d = MML::Real(1.0) / b;
				Real h = d;
				for (int i = 1; i <= MAXIT; i++)
				{
					Real a = -Real(i) * (nm1 + i);
					b += 2.0;
					d = MML::Real(1.0) / (a * d + b);
					c = b + a / c;
					Real del = c * d;
					h *= del;
					if (std::abs(del - 1.0) <= eps)
						return h * std::exp(-x);
				}
				return h * std::exp(-x);
			}

			Real ans = (nm1 != 0) ? MML::Real(1.0) / nm1 : -std::log(x) - EulerGamma;
			Real fact = 1.0;
			for (int i = 1; i <= MAXIT; i++)
			{
				fact *= -x / i;
				Real del;
				if (i != nm1)
					del = -fact / (i - nm1);
				else
				{
					Real psi = -EulerGamma;
					for (int ii = 1; ii <= nm1; ii++) psi += MML::Real(1.0) / ii;
					del = fact * (-std::log(x) + psi);
				}
				ans += del;
				if (std::abs(del) < std::abs(ans) * eps) return ans;
			}
			return ans;
		}

		/// @brief Exponential integral E_1(x) = ∫_x^∞ e^{-t} / t dt, x > 0.
		inline Real ExpIntegralE1(Real x) { return ExpIntegralEn(1, x); }

		/// @brief Exponential integral Ei(x) = -PV ∫_{-x}^∞ e^{-t}/t dt, x ≠ 0.
		/// @details Series for moderate x > 0, least-term asymptotic for large x; Ei(x) = -E_1(-x) for x < 0.
		///          Split exponential scaling retains finite large positive values. Genuine final overflow is +infinity.
		///          Ei(+infinity) = +infinity and Ei(-infinity) = -0.
		/// @throws DomainError for NaN or either signed zero; ConvergenceError if a positive-input budget is exhausted.
		inline Real ExpIntegralEi(Real x)
		{
			static const Real EulerGamma = static_cast<MML::Real>(0.5772156649015328606);
			const Real eps = std::numeric_limits<Real>::epsilon();
			const Real fpmin = std::numeric_limits<Real>::min() / eps;
			if (std::isnan(x) || x == 0.0)
				throw DomainError("SpecialFunctions::ExpIntegralEi - x must be nonzero and not NaN");
			if (std::isinf(x)) return x > 0 ? x : -Real(0);
			if (x < 0.0) return -ExpIntegralEn(1, -x);
			if (x < fpmin) return std::log(x) + EulerGamma;

			if (x <= -std::log(eps))
			{
				Real sum = 0.0, fact = 1.0;
				for (int k = 1; k <= Precision::ExpIntegralEiMaxIterations; k++)
				{
					fact *= x / k;
					Real term = fact / k;
					sum += term;
					if (IsNegligibleRelativeTo(term, sum, eps))
						return sum + std::log(x) + EulerGamma;
				}
				throw ConvergenceError("SpecialFunctions::ExpIntegralEi - positive series exhausted its iteration budget",
					Precision::ExpIntegralEiMaxIterations, std::abs((fact / Precision::ExpIntegralEiMaxIterations) / sum));
			}

			Real sum = 0.0, term = 1.0;
			bool completed = false;
			for (int k = 1; k <= Precision::ExpIntegralEiMaxIterations; k++)
			{
				Real prev = term;
				term *= Real(k) / x;
				if (term < eps) { completed = true; break; }
				if (term < prev) sum += term;
				else { sum -= prev; completed = true; break; }
			}
			if (!completed)
				throw ConvergenceError("SpecialFunctions::ExpIntegralEi - asymptotic expansion exhausted its iteration budget",
					Precision::ExpIntegralEiMaxIterations, term / (1 + sum));
			// Each half exponential stays representable throughout the finite Ei range.
			const Real halfExp = std::exp(x / 2);
			return halfExp * (((1 + sum) / x) * halfExp);
		}

		///////////////////////////////////////////////////////////////////////
		///                   Sine and cosine integrals                     ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Computes both Si(x) = ∫_0^x sin t/t dt and Ci(x) = γ + ln|x| + ∫_0^x (cos t − 1)/t dt.
		/// @details Power series for |x| ≤ 2, complex continued fraction otherwise (Numerical Recipes).
		///          Real extension: Si is odd and Ci is even. Si preserves signed zero; Ci(±0) = -∞.
		///          At ±∞, Si = ±π/2 and Ci = positive zero.
		/// @throws DomainError for NaN, before modifying either output.
		inline void SinCosIntegral(Real x, Real& si, Real& ci)
		{
			if (std::isnan(x)) throw DomainError("SpecialFunctions::SinCosIntegral - x must not be NaN");
			static const Real EulerGamma = static_cast<MML::Real>(0.5772156649015328606);
			const Real PIBY2 = Constants::PI / 2.0;
			const Real eps = std::numeric_limits<Real>::epsilon();
			const Real fpmin = std::numeric_limits<Real>::min() * static_cast<MML::Real>(4.0);
			const Real big = std::numeric_limits<Real>::max() * eps;
			const Real TMIN = static_cast<MML::Real>(2.0);
			const int MAXIT = 100;

			Real t = std::abs(x);
			if (t == 0.0) { si = x; ci = -std::numeric_limits<Real>::infinity(); return; }
			if (std::isinf(t)) { si = std::copysign(PIBY2, x); ci = Real(0); return; }

			if (t > TMIN)
			{
				std::complex<Real> b(1.0, t), c(big, 0.0), d, h, del;
				d = h = static_cast<MML::Real>(1.0) / b;
				for (int i = 1; i < MAXIT; i++)
				{
					Real a = -Real(i) * Real(i);
					b += std::complex<Real>(2.0, 0.0);
					d = static_cast<MML::Real>(1.0) / (a * d + b);
					c = b + a / c;
					del = c * d;
					h *= del;
					if (std::abs(std::real(del) - 1.0) + std::abs(std::imag(del)) <= eps) break;
				}
				h = std::complex<Real>(std::cos(t), -std::sin(t)) * h;
				ci = -std::real(h);
				si = PIBY2 + std::imag(h);
			}
			else
			{
				Real sum = 0.0, sums = 0.0, sumc = 0.0, sign = 1.0, fact = 1.0;
				bool odd = true;
				if (t < std::sqrt(fpmin)) { sumc = 0.0; sums = t; }
				else
				{
					for (int k = 1; k <= MAXIT; k++)
					{
						fact *= t / k;
						Real term = fact / k;
						sum += sign * term;
						Real err = term / std::abs(sum);
						if (odd) { sign = -sign; sums = sum; sum = sumc; }
						else { sumc = sum; sum = sums; }
						if (err < eps) break;
						odd = !odd;
					}
				}
				si = sums;
				ci = sumc + std::log(t) + EulerGamma;
			}
			if (x < 0.0) si = -si;
		}

		/// @brief Sine integral Si(x); preserves signed zero and has limits ±π/2 at ±∞.
		/// @throws DomainError for NaN.
		inline Real SineIntegral(Real x) { Real si, ci; SinCosIntegral(x, si, ci); return si; }

		/// @brief Real, even cosine integral Ci(x); Ci(±0) = -∞ and Ci(±∞) = positive zero.
		/// @throws DomainError for NaN.
		inline Real CosineIntegral(Real x) { Real si, ci; SinCosIntegral(x, si, ci); return ci; }

		///////////////////////////////////////////////////////////////////////
		///                   Dawson and Fresnel integrals                  ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Dawson's integral D(x) = e^{-x²} ∫_0^x e^{t²} dt.
		/// @details Maclaurin series near zero, sampling-theorem approximation at ordinary magnitudes
		///          (Numerical Recipes), and a reciprocal asymptotic expansion for |x| >= DawsonAsymptoticLimit.
		///          Preserves signed zero, including the limits at +/-infinity.
		/// @throws DomainError for NaN; ConvergenceError if the asymptotic tail cannot meet its stopping criterion.
		inline Real Dawson(Real x)
		{
			if (std::isnan(x)) throw DomainError("SpecialFunctions::Dawson - x must not be NaN");
			if (x == 0) return x;
			if (std::isinf(x)) return std::copysign(Real(0), x);
			const Real xx = std::abs(x);
			if (xx >= Precision::DawsonAsymptoticLimit)
			{
				// DLMF 7.12: D(x) ~ (1/(2*x)) * sum (2k-1)!!/(2*x*x)^k.
				// Reciprocals avoid forming x*x or 2*x, either of which can overflow.
				const Real inverse = 1 / xx;
				Real sum = 1, term = 1;
				for (int k = 1; k <= Precision::DawsonMaxTerms; ++k)
				{
					const Real previous = term;
					term *= (Real(2 * k - 1) / 2) * inverse * inverse;
					if (term > previous)
						throw ConvergenceError("SpecialFunctions::Dawson - asymptotic terms grew before convergence",
							k, term / sum);
					sum += term;
					if (IsNegligibleRelativeTo(term, sum, Precision::Eps))
						return std::copysign((Real(0.5) / xx) * sum, x);
				}
				throw ConvergenceError("SpecialFunctions::Dawson - asymptotic expansion exhausted its term budget",
					Precision::DawsonMaxTerms, term / sum);
			}

			const int NMAX = 6;
			const Real H = Precision::DawsonSamplingStep;
			const Real A1 = static_cast<MML::Real>(2.0) / static_cast<MML::Real>(3.0), A2 = static_cast<MML::Real>(0.4), A3 = static_cast<MML::Real>(2.0) / static_cast<MML::Real>(7.0);
			static const std::array<Real, NMAX> ctab = [] {
				std::array<Real, NMAX> c{};
				for (int i = 0; i < NMAX; i++) { Real t = (MML::Real(2.0) * i + MML::Real(1.0)) * Precision::DawsonSamplingStep; c[i] = std::exp(-t * t); }
				return c;
			}();

			if (xx < Precision::DawsonSmallLimit)
			{
				Real x2 = x * x;
				return x * (MML::Real(1.0) - A1 * x2 * (MML::Real(1.0) - A2 * x2 * (MML::Real(1.0) - A3 * x2)));
			}

			int n0 = 2 * int(0.5 * xx / H + 0.5);
			Real xp = xx - n0 * H;
			Real e1 = std::exp(MML::Real(2.0) * xp * H);
			Real e2 = e1 * e1;
			Real d1 = Real(n0 + 1);
			Real d2 = d1 - MML::Real(2.0);
			Real sum = 0.0;
			for (int i = 0; i < NMAX; i++, d1 += 2.0, d2 -= 2.0, e1 *= e2)
				sum += ctab[i] * (e1 / d1 + MML::Real(1.0) / (d2 * e1));

			Real val = (static_cast<MML::Real>(1.0) / std::sqrt(Constants::PI)) * std::exp(-xp * xp) * sum;
			return (x >= 0.0) ? val : -val;
		}

		/// @brief Computes both Fresnel integrals S(x) = ∫_0^x sin(π t²/2) dt and C(x) = ∫_0^x cos(π t²/2) dt.
		/// @details Power series for |x| ≤ 1.5, complex continued fraction otherwise (Numerical Recipes).
		inline void FresnelSC(Real x, Real& s, Real& c)
		{
			const Real eps = std::numeric_limits<Real>::epsilon();
			const Real fpmin = std::numeric_limits<Real>::min();
			const Real big = std::numeric_limits<Real>::max() * eps;
			const Real XMIN = static_cast<MML::Real>(1.5);
			const Real PIBY2 = Constants::PI / 2.0;
			const int MAXIT = 100;

			Real ax = std::abs(x);
			if (ax < std::sqrt(fpmin)) { s = 0.0; c = ax; }
			else if (ax <= XMIN)
			{
				Real sum = 0.0, sums = 0.0, sumc = ax, sign = 1.0;
				Real fact = PIBY2 * ax * ax;
				bool odd = true;
				Real term = ax;
				int n = 3;
				for (int k = 1; k <= MAXIT; k++)
				{
					term *= fact / k;
					sum += sign * term / n;
					Real test = std::abs(sum) * eps;
					if (odd) { sign = -sign; sums = sum; sum = sumc; }
					else { sumc = sum; sum = sums; }
					if (term < test) break;
					odd = !odd;
					n += 2;
				}
				s = sums; c = sumc;
			}
			else
			{
				Real pix2 = Constants::PI * ax * ax;
				std::complex<Real> b(1.0, -pix2), cc(big, 0.0), d, h, del, cs;
				d = h = static_cast<MML::Real>(1.0) / b;
				int n = -1;
				for (int k = 2; k <= MAXIT; k++)
				{
					n += 2;
					Real a = -Real(n) * Real(n + 1);
					b += std::complex<Real>(4.0, 0.0);
					d = static_cast<MML::Real>(1.0) / (a * d + b);
					cc = b + a / cc;
					del = cc * d;
					h *= del;
					if (std::abs(std::real(del) - 1.0) + std::abs(std::imag(del)) <= eps) break;
				}
				h *= std::complex<Real>(ax, -ax);
				cs = std::complex<Real>(MML::Real(0.5), MML::Real(0.5)) * (static_cast<MML::Real>(1.0) - std::complex<Real>(std::cos(MML::Real(0.5) * pix2), std::sin(MML::Real(0.5) * pix2)) * h);
				c = std::real(cs); s = std::imag(cs);
			}
			if (x < 0.0) { c = -c; s = -s; }
		}

		/// @brief Fresnel integral C(x) = ∫_0^x cos(π t²/2) dt.
		inline Real FresnelC(Real x) { Real s, c; FresnelSC(x, s, c); return c; }

		/// @brief Fresnel integral S(x) = ∫_0^x sin(π t²/2) dt.
		inline Real FresnelS(Real x) { Real s, c; FresnelSC(x, s, c); return s; }

		///////////////////////////////////////////////////////////////////////
		///                        Lambert W function                       ///
		///////////////////////////////////////////////////////////////////////

		namespace SpecialFunctionsDetail
		{
			inline Real LambertWBranchSeries(Real p)
			{
				// Corless et al. (1996): Puiseux expansion in p = +/-sqrt(2(1+e*x)).
				return -1 + p * (1 + p * (-Real(1) / 3 + p * (Real(11) / 72 +
					p * (-Real(43) / 540 + p * (Real(769) / 17280 + p * (-Real(221) / 8505 +
					p * (Real(680863) / 43545600 + p * (-Real(1963) / 204120 +
					p * (Real(226287557) / 37623398400LL)))))))));
			}

			inline Real LambertW(Real x, bool secondary)
			{
				constexpr Real minX = Real(-0.367879441171442321595523770161460867L);
				if (std::isnan(x) || x < minX || (secondary && x >= 0))
					throw DomainError(secondary ? "SpecialFunctions::LambertWm1 - x must be in [-1/e, 0)"
						: "SpecialFunctions::LambertW0 - x must be >= -1/e");
				if (x == minX) return -1;
				if (x == 0 || std::isinf(x)) return x;

				Real w;
				if (x < Precision::LambertWBranchInitialLimit)
				{
					// Split -1/e: the first two terms are exact in float too. Subtractions
					// retain the true distance even at the first valid neighbor of the rounded endpoint.
					const Real distance = ((x - Real(-0.367879450321197509765625L)) -
						Real(0.00000000914975561983055740711279213428497314453125L)) -
						Real(-4.316604561772742530017307842755630178345078368e-16L);
					const Real p = std::sqrt(2 * Real(2.718281828459045235360287471352662498L) * distance);
					w = LambertWBranchSeries(secondary ? -p : p);
					if (p <= Precision::LambertWSeriesLimit)
						return w;
				}
				else if (!secondary && x < Precision::LambertWLogInitialLimit)
					w = x / (1 + x);
				else
				{
					const Real l1 = std::log(std::abs(x));
					const Real l2 = std::log(std::abs(l1));
					w = l1 - l2 + l2 / l1;
				}

				const Real logX = std::log(std::abs(x));
				Real residual = std::numeric_limits<Real>::infinity();
				for (int i = 0; i < Precision::LambertWMaxIterations; ++i)
				{
					// Halley on w + log(|w|) = log(|x|) avoids exponential range failures.
					const Real g = w + std::log(std::abs(w)) - logX;
					const Real wp1 = w + 1;
					const Real dw = (g * (w / wp1)) / (1 + g / (2 * wp1 * wp1));
					w -= dw;
					if (!std::isfinite(w) || !std::isfinite(dw) ||
						(secondary ? w >= -1 : w <= -1) || ((w > 0) != (x > 0)))
						throw NumericalMethodError("SpecialFunctions Lambert W refinement left its finite real branch");
					residual = std::abs(dw);
					if (IsNegligibleRelativeTo(dw, 1 + std::abs(w), Precision::Eps))
						return w;
				}
				throw ConvergenceError("SpecialFunctions Lambert W exhausted its iteration budget",
					Precision::LambertWMaxIterations, residual);
			}
		} // namespace SpecialFunctionsDetail

		/// @brief Principal real branch W_0(x), solving w*exp(w) = x for x >= -1/e.
		/// @details Only the rounded -1/e endpoint returns exactly -1. Nearby inputs use a
		///          cancellation-safe Puiseux expansion; other inputs use logarithmic Halley refinement.
		///          Zero returns zero; positive infinity returns positive infinity.
		/// @throws DomainError for NaN or x below the endpoint; NumericalMethodError on arithmetic failure;
		///         ConvergenceError on exhausted refinement budget.
		inline Real LambertW0(Real x) { return SpecialFunctionsDetail::LambertW(x, false); }

		/// @brief Secondary real branch W_{-1}(x), for -1/e <= x < 0.
		/// @details Shares the endpoint, branch expansion and explicit failure contracts of LambertW0.
		/// @throws DomainError for NaN or x outside [-1/e, 0), including either signed zero.
		inline Real LambertWm1(Real x) { return SpecialFunctionsDetail::LambertW(x, true); }

		///////////////////////////////////////////////////////////////////////
		///                        Airy functions                           ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Computes all four Airy functions Ai(x), Bi(x), Ai'(x), Bi'(x) at once.
		/// @details Maclaurin series for |x| < 6, DLMF asymptotic expansions beyond.
		inline void Airy(Real x, Real& ai, Real& bi, Real& aip, Real& bip)
		{
			const Real sqrtPi = std::sqrt(Constants::PI);
			const Real c1 = static_cast<MML::Real>(0.3550280538878172);   // Ai(0)
			const Real c2 = static_cast<MML::Real>(0.2588194037928068);   // -Ai'(0)
			const Real eps = std::numeric_limits<Real>::epsilon();

			if (std::abs(x) < static_cast<MML::Real>(6.0))
			{
				const Real sqrt3 = std::sqrt(static_cast<MML::Real>(3.0));
				if (x == 0.0) { ai = c1; bi = sqrt3 * c1; aip = -c2; bip = sqrt3 * c2; return; }
				Real x3 = x * x * x;
				Real F = 1.0, G = x, Fp = 0.0, Gp = 1.0, tF = 1.0, tG = x;
				for (int k = 1; k <= 300; k++)
				{
					tF *= x3 / (Real(3 * k - 1) * Real(3 * k));
					tG *= x3 / (Real(3 * k) * Real(3 * k + 1));
					F += tF; G += tG;
					Fp += Real(3 * k) * tF / x;
					Gp += Real(3 * k + 1) * tG / x;
					if (std::abs(tF) <= eps * std::abs(F) && std::abs(tG) <= eps * std::abs(G)) break;
				}
				ai = c1 * F - c2 * G;
				bi = sqrt3 * (c1 * F + c2 * G);
				aip = c1 * Fp - c2 * Gp;
				bip = sqrt3 * (c1 * Fp + c2 * Gp);
				return;
			}

			const Real t = std::abs(x);
			const Real zeta = (static_cast<MML::Real>(2.0) / static_cast<MML::Real>(3.0)) * std::pow(t, static_cast<MML::Real>(1.5));
			const Real t14 = std::pow(t, static_cast<MML::Real>(0.25));
			const Real zinv = static_cast<MML::Real>(1.0) / zeta;

			if (x > 0.0)
			{
				Real Sc = 0, ScAlt = 0, Sv = 0, SvAlt = 0;
				Real cm = 1.0, zpow = 1.0, prev = Constants::PosInf;
				for (int m = 0; m <= 40; m++)
				{
					Real vm = (m == 0) ? static_cast<MML::Real>(1.0) : -cm * Real(6 * m + 1) / Real(6 * m - 1);
					Real termC = cm * zpow, termV = vm * zpow;
					if (m > 0 && std::abs(termC) > prev) break;   // asymptotic divergence
					Real sgn = (m % 2 == 0) ? static_cast<MML::Real>(1.0) : static_cast<MML::Real>(-1.0);
					Sc += termC; ScAlt += sgn * termC;
					Sv += termV; SvAlt += sgn * termV;
					prev = std::abs(termC);
					if (std::abs(termC) < eps) break;
					int k = m + 1;
					cm *= Real(6 * k - 5) * Real(6 * k - 3) * Real(6 * k - 1) / (Real(2 * k - 1) * static_cast<MML::Real>(216.0) * Real(k));
					zpow *= zinv;
				}
				Real emz = std::exp(-zeta), epz = std::exp(zeta);
				ai = emz / (static_cast<MML::Real>(2.0) * sqrtPi * t14) * ScAlt;
				bi = epz / (sqrtPi * t14) * Sc;
				aip = -t14 * emz / (static_cast<MML::Real>(2.0) * sqrtPi) * SvAlt;
				bip = t14 * epz / sqrtPi * Sv;
				return;
			}

			// x < 0
			const Real theta = zeta - Constants::PI / static_cast<MML::Real>(4.0);
			const Real cth = std::cos(theta), sth = std::sin(theta);
			Real P = 0, Q = 0, R = 0, S = 0;
			Real cm = 1.0, zpow = 1.0, prev = Constants::PosInf;
			for (int m = 0; m <= 40; m++)
			{
				Real vm = (m == 0) ? static_cast<MML::Real>(1.0) : -cm * Real(6 * m + 1) / Real(6 * m - 1);
				Real termC = cm * zpow, termV = vm * zpow;
				if (m > 0 && std::abs(termC) > prev) break;
				Real sgn = ((m / 2) % 2 == 0) ? static_cast<MML::Real>(1.0) : static_cast<MML::Real>(-1.0);
				if (m % 2 == 0) { P += sgn * termC; R += sgn * termV; }
				else { Q += sgn * termC; S += sgn * termV; }
				prev = std::abs(termC);
				if (std::abs(termC) < eps) break;
				int k = m + 1;
				cm *= Real(6 * k - 5) * Real(6 * k - 3) * Real(6 * k - 1) / (Real(2 * k - 1) * static_cast<MML::Real>(216.0) * Real(k));
				zpow *= zinv;
			}
			Real inv = static_cast<MML::Real>(1.0) / (sqrtPi * t14);
			ai = inv * (cth * P + sth * Q);
			bi = inv * (-sth * P + cth * Q);
			Real invd = t14 / sqrtPi;
			aip = invd * (sth * R - cth * S);
			bip = invd * (cth * R + sth * S);
		}

		/// @brief Airy function Ai(x).
		inline Real AiryAi(Real x) { Real ai, bi, aip, bip; Airy(x, ai, bi, aip, bip); return ai; }
		/// @brief Airy function Bi(x).
		inline Real AiryBi(Real x) { Real ai, bi, aip, bip; Airy(x, ai, bi, aip, bip); return bi; }
		/// @brief Derivative Ai'(x).
		inline Real AiryAiPrime(Real x) { Real ai, bi, aip, bip; Airy(x, ai, bi, aip, bip); return aip; }
		/// @brief Derivative Bi'(x).
		inline Real AiryBiPrime(Real x) { Real ai, bi, aip, bip; Airy(x, ai, bi, aip, bip); return bip; }

		namespace detail
		{
			/// @brief Asymptotic locator T(t) for Airy zeros (DLMF 9.9.6).
			inline Real AiryZeroT(Real t)
			{
				Real u = static_cast<MML::Real>(1.0) / (t * t);
				return std::pow(t, static_cast<MML::Real>(2.0) / static_cast<MML::Real>(3.0)) *
					(static_cast<MML::Real>(1.0) + u * (static_cast<MML::Real>(5.0) / 48 - u * (static_cast<MML::Real>(5.0) / 36 - u * (static_cast<MML::Real>(77125.0) / 82944))));
			}
		}

		/// @brief The n-th real (negative) zero a_n of Ai (n ≥ 1).
		inline Real AiryAiZero(int n)
		{
			if (n < 1) throw DomainError("SpecialFunctions::AiryAiZero - n must be >= 1");
			Real t = static_cast<MML::Real>(3.0) * Constants::PI * (static_cast<MML::Real>(4.0) * n - static_cast<MML::Real>(1.0)) / static_cast<MML::Real>(8.0);
			Real x = -detail::AiryZeroT(t);
			for (int i = 0; i < 60; i++)
			{
				Real ai, bi, aip, bip; Airy(x, ai, bi, aip, bip);
				Real dx = ai / aip;
				x -= dx;
				if (IsNegligibleRelativeTo(dx, static_cast<MML::Real>(1.0) + std::abs(x), Precision::RoundoffZeroThreshold)) break;
			}
			return x;
		}

		/// @brief The n-th real (negative) zero b_n of Bi (n ≥ 1).
		inline Real AiryBiZero(int n)
		{
			if (n < 1) throw DomainError("SpecialFunctions::AiryBiZero - n must be >= 1");
			Real t = static_cast<MML::Real>(3.0) * Constants::PI * (static_cast<MML::Real>(4.0) * n - static_cast<MML::Real>(3.0)) / static_cast<MML::Real>(8.0);
			Real x = -detail::AiryZeroT(t);
			for (int i = 0; i < 60; i++)
			{
				Real ai, bi, aip, bip; Airy(x, ai, bi, aip, bip);
				Real dx = bi / bip;
				x -= dx;
				if (IsNegligibleRelativeTo(dx, static_cast<MML::Real>(1.0) + std::abs(x), Precision::RoundoffZeroThreshold)) break;
			}
			return x;
		}
	} // namespace SpecialFunctions
}
#endif
