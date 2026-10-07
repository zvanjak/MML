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

#include <mml/MMLBase.h>

namespace MML
{
	/// @brief Portable special functions (gamma/beta family, polygamma, error-function inverses).
	/// @details These are real-valued, double-precision implementations that do not depend on the
	///          C++17 `std::` special-math functions (which are unavailable on some platforms).
	///          The regularized incomplete gamma/beta routines and their inverses are the shared,
	///          public home for machinery previously duplicated privately inside the statistical
	///          distributions.
	namespace SpecialFunctions
	{
		/// @brief Gamma function Γ(x).
		inline Real Gamma(Real x) { return std::tgamma(x); }

		/// @brief Natural log of |Γ(x)|.
		inline Real LnGamma(Real x) { return std::lgamma(x); }

		/// @brief Beta function B(a,b) = Γ(a)Γ(b)/Γ(a+b).
		inline Real Beta(Real a, Real b) { return std::exp(std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b)); }

		/// @brief Natural log of the Beta function.
		inline Real LnBeta(Real a, Real b) { return std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b); }

		///////////////////////////////////////////////////////////////////////
		///                  Incomplete gamma family                        ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Regularized lower incomplete gamma P(a,x) = γ(a,x)/Γ(a).
		/// @details Series representation for x < a+1, continued fraction otherwise (Numerical Recipes).
		/// @throws DomainError if a <= 0.
		inline Real RegularizedGammaP(Real a, Real x)
		{
			if (a <= 0.0) throw DomainError("SpecialFunctions::RegularizedGammaP - a must be > 0");
			if (x < 0.0) throw DomainError("SpecialFunctions::RegularizedGammaP - x must be >= 0");
			if (x == 0.0) return 0.0;

			const Real gln = std::lgamma(a);
			if (x < a + 1.0)
			{
				Real ap = a;
				Real sum = 1.0 / a;
				Real del = sum;
				for (int n = 1; n <= 300; n++)
				{
					ap += 1.0;
					del *= x / ap;
					sum += del;
					if (std::abs(del) < std::abs(sum) * REAL(1e-16)) break;
				}
				return sum * std::exp(-x + a * std::log(x) - gln);
			}
			else
			{
				const Real fpmin = std::numeric_limits<Real>::min() / std::numeric_limits<Real>::epsilon();
				Real b = x + 1.0 - a;
				Real c = 1.0 / fpmin;
				Real d = 1.0 / b;
				Real h = d;
				for (int n = 1; n <= 300; n++)
				{
					Real an = -n * (n - a);
					b += 2.0;
					d = an * d + b;
					if (std::abs(d) < fpmin) d = fpmin;
					c = b + an / c;
					if (std::abs(c) < fpmin) c = fpmin;
					d = 1.0 / d;
					Real del = d * c;
					h *= del;
					if (std::abs(del - 1.0) < REAL(1e-16)) break;
				}
				return 1.0 - h * std::exp(-x + a * std::log(x) - gln);
			}
		}

		/// @brief Regularized upper incomplete gamma Q(a,x) = Γ(a,x)/Γ(a) = 1 - P(a,x).
		inline Real RegularizedGammaQ(Real a, Real x) { return 1.0 - RegularizedGammaP(a, x); }

		/// @brief Lower incomplete gamma γ(a,x) = P(a,x)·Γ(a).
		inline Real IncompleteGammaLower(Real a, Real x) { return RegularizedGammaP(a, x) * std::tgamma(a); }

		/// @brief Upper incomplete gamma Γ(a,x) = Q(a,x)·Γ(a).
		inline Real IncompleteGammaUpper(Real a, Real x) { return RegularizedGammaQ(a, x) * std::tgamma(a); }

		///////////////////////////////////////////////////////////////////////
		///                  Incomplete beta family                         ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Regularized incomplete beta I_x(a,b) = B(x;a,b)/B(a,b).
		/// @details Continued fraction via modified Lentz's method (Numerical Recipes),
		///          using the symmetry relation I_x(a,b) = 1 - I_{1-x}(b,a) for fast convergence.
		/// @throws DomainError if a <= 0 or b <= 0, or x outside [0,1].
		inline Real RegularizedBetaI(Real x, Real a, Real b)
		{
			if (a <= 0.0 || b <= 0.0) throw DomainError("SpecialFunctions::RegularizedBetaI - a,b must be > 0");
			if (x < 0.0 || x > 1.0) throw DomainError("SpecialFunctions::RegularizedBetaI - x must be in [0,1]");
			if (x == 0.0) return 0.0;
			if (x == 1.0) return 1.0;

			if (x > (a + 1.0) / (a + b + 2.0))
				return 1.0 - RegularizedBetaI(1.0 - x, b, a);

			const Real logBeta = std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
			const Real front = std::exp(a * std::log(x) + b * std::log(1.0 - x) - logBeta) / a;

			const Real fpmin = std::numeric_limits<Real>::min() / std::numeric_limits<Real>::epsilon();
			const Real qab = a + b;
			const Real qap = a + 1.0;
			const Real qam = a - 1.0;
			Real c = 1.0;
			Real d = 1.0 - qab * x / qap;
			if (std::abs(d) < fpmin) d = fpmin;
			d = 1.0 / d;
			Real h = d;
			for (int m = 1; m <= 300; m++)
			{
				int m2 = 2 * m;
				Real aa = m * (b - m) * x / ((qam + m2) * (a + m2));
				d = 1.0 + aa * d;
				if (std::abs(d) < fpmin) d = fpmin;
				c = 1.0 + aa / c;
				if (std::abs(c) < fpmin) c = fpmin;
				d = 1.0 / d;
				h *= d * c;

				aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
				d = 1.0 + aa * d;
				if (std::abs(d) < fpmin) d = fpmin;
				c = 1.0 + aa / c;
				if (std::abs(c) < fpmin) c = fpmin;
				d = 1.0 / d;
				Real del = d * c;
				h *= del;
				if (std::abs(del - 1.0) < REAL(1e-16)) break;
			}
			return front * h;
		}

		/// @brief Incomplete beta B(x;a,b) = I_x(a,b)·B(a,b).
		inline Real IncompleteBeta(Real x, Real a, Real b) { return RegularizedBetaI(x, a, b) * Beta(a, b); }

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
				x = 1.0 - x;
			}
			while (x < 6.0) { result -= 1.0 / x; x += 1.0; }

			const Real inv = 1.0 / x;
			const Real inv2 = inv * inv;
			result += std::log(x) - 0.5 * inv;
			result -= inv2 * (REAL(1.0) / 12 - inv2 * (REAL(1.0) / 120 - inv2 * (REAL(1.0) / 252
				- inv2 * (REAL(1.0) / 240 - inv2 * (REAL(1.0) / 132)))));
			return result;
		}

		/// @brief Trigamma ψ'(x).
		/// @details Reflection for x < 0.5, recurrence up to x >= 6, then asymptotic series.
		inline Real Trigamma(Real x)
		{
			if (x < 0.5)
			{
				const Real s = std::sin(Constants::PI * x);
				return (Constants::PI * Constants::PI) / (s * s) - Trigamma(1.0 - x);
			}
			Real result = 0.0;
			while (x < 6.0) { result += 1.0 / (x * x); x += 1.0; }

			const Real inv = 1.0 / x;
			const Real inv2 = inv * inv;
			result += inv + 0.5 * inv2 + inv * inv2 * (REAL(1.0) / 6 - inv2 * (REAL(1.0) / 30
				- inv2 * (REAL(1.0) / 42 - inv2 * (REAL(1.0) / 30))));
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
			const Real sign = (n % 2 == 0) ? -1.0 : 1.0;   // (-1)^(n+1)

			// Sum (x+k)^{-s} while shifting x up to a safe base for the asymptotic tail.
			Real sum = 0.0;
			const Real target = 10.0 + n;
			while (x < target) { sum += std::pow(x, -(Real)s); x += 1.0; }

			// Euler–Maclaurin tail of the Hurwitz zeta ζ(s, x).
			Real tail = std::pow(x, 1.0 - s) / (s - 1) + 0.5 * std::pow(x, -(Real)s);
			static const Real emC[] = { REAL(1.0) / 12, REAL(-1.0) / 720, REAL(1.0) / 30240,
										REAL(-1.0) / 1209600, REAL(1.0) / 47900160 };
			Real poch = s;                          // (s)_{2j-1}, starts at (s)_1 = s
			Real xpow = std::pow(x, -(Real)(s + 1)); // x^{-(s+2j-1)}, starts at x^{-(s+1)}
			const Real xinv2 = 1.0 / (x * x);
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
			if (x <= -1.0) { if (x == -1.0) return -std::numeric_limits<Real>::infinity(); throw DomainError("SpecialFunctions::ErfInv - |x| must be <= 1"); }
			if (x >= 1.0) { if (x == 1.0) return std::numeric_limits<Real>::infinity(); throw DomainError("SpecialFunctions::ErfInv - |x| must be <= 1"); }
			if (x == 0.0) return 0.0;

			Real w = -std::log((1.0 - x) * (1.0 + x));
			Real p;
			if (w < 5.0)
			{
				w -= 2.5;
				p = REAL(2.81022636e-08);
				p = REAL(3.43273939e-07) + p * w;
				p = REAL(-3.5233877e-06) + p * w;
				p = REAL(-4.39150654e-06) + p * w;
				p = REAL(0.00021858087) + p * w;
				p = REAL(-0.00125372503) + p * w;
				p = REAL(-0.00417768164) + p * w;
				p = REAL(0.246640727) + p * w;
				p = REAL(1.50140941) + p * w;
			}
			else
			{
				w = std::sqrt(w) - 3.0;
				p = REAL(-0.000200214257);
				p = REAL(0.000100950558) + p * w;
				p = REAL(0.00134934322) + p * w;
				p = REAL(-0.00367342844) + p * w;
				p = REAL(0.00573950773) + p * w;
				p = REAL(-0.0076224613) + p * w;
				p = REAL(0.00943887047) + p * w;
				p = REAL(1.00167406) + p * w;
				p = REAL(2.83297682) + p * w;
			}
			Real res = p * x;

			// Halley refinement against erf (2 iterations reach full double precision).
			const Real twoOverSqrtPi = REAL(2.0) / std::sqrt(Constants::PI);
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
			if (y <= 0.0) { if (y == 0.0) return std::numeric_limits<Real>::infinity(); throw DomainError("SpecialFunctions::ErfcInv - y must be in (0,2)"); }
			if (y >= 2.0) { if (y == 2.0) return -std::numeric_limits<Real>::infinity(); throw DomainError("SpecialFunctions::ErfcInv - y must be in (0,2)"); }
			return ErfInv(1.0 - y);
		}

		///////////////////////////////////////////////////////////////////////
		///             Inverses of the regularized functions               ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Inverse of the regularized lower incomplete gamma: solves P(a,x)=p for x.
		/// @details Initial guess (Wilson–Hilferty / Numerical Recipes) refined by Halley iteration.
		inline Real RegularizedGammaPInv(Real a, Real p)
		{
			if (a <= 0.0) throw DomainError("SpecialFunctions::RegularizedGammaPInv - a must be > 0");
			if (p <= 0.0) return 0.0;
			if (p >= 1.0) return std::numeric_limits<Real>::infinity();

			const Real gln = std::lgamma(a);
			const Real a1 = a - 1.0;
			Real x;
			if (a > 1.0)
			{
				Real pp = (p < 0.5) ? p : 1.0 - p;
				Real t = std::sqrt(-2.0 * std::log(pp));
				Real x0 = (2.30753 + t * 0.27061) / (1.0 + t * (0.99229 + t * 0.04481)) - t;
				if (p < 0.5) x0 = -x0;
				x = a * std::pow(1.0 - 1.0 / (9.0 * a) - x0 / (3.0 * std::sqrt(a)), 3.0);
				x = (x <= 0.0) ? REAL(1e-8) : x;
			}
			else
			{
				Real t = 1.0 - a * (0.253 + a * 0.12);
				if (p < t) x = std::pow(p / t, 1.0 / a);
				else x = 1.0 - std::log(1.0 - (p - t) / (1.0 - t));
			}

			for (int j = 0; j < 15; j++)
			{
				if (x <= 0.0) { x = REAL(1e-8); break; }
				Real err = RegularizedGammaP(a, x) - p;
				Real deriv = std::exp(-x + a1 * std::log(x) - gln);   // dP/dx
				Real u = err / deriv;
				Real correction = u / (1.0 - 0.5 * std::min<Real>(REAL(1.0), u * (a1 / x - REAL(1.0))));
				x -= correction;
				if (x <= 0.0) x = 0.5 * (x + correction);
				if (std::abs(correction) < REAL(1e-14) * std::abs(x)) break;
			}
			return x;
		}

		/// @brief Inverse of the regularized incomplete beta: solves I_x(a,b)=p for x.
		/// @details Initial guess (Numerical Recipes) refined by Halley iteration.
		inline Real RegularizedBetaIInv(Real p, Real a, Real b)
		{
			if (a <= 0.0 || b <= 0.0) throw DomainError("SpecialFunctions::RegularizedBetaIInv - a,b must be > 0");
			if (p <= 0.0) return 0.0;
			if (p >= 1.0) return 1.0;

			const Real a1 = a - 1.0;
			const Real b1 = b - 1.0;
			Real x, t, w;
			if (a >= 1.0 && b >= 1.0)
			{
				Real pp = (p < 0.5) ? p : 1.0 - p;
				t = std::sqrt(-2.0 * std::log(pp));
				Real xg = (2.30753 + t * 0.27061) / (1.0 + t * (0.99229 + t * 0.04481)) - t;
				if (p < 0.5) xg = -xg;
				Real al = (xg * xg - 3.0) / 6.0;
				Real h = 2.0 / (1.0 / (2.0 * a - 1.0) + 1.0 / (2.0 * b - 1.0));
				w = (xg * std::sqrt(al + h) / h) - (1.0 / (2.0 * b - 1.0) - 1.0 / (2.0 * a - 1.0)) * (al + 5.0 / 6.0 - 2.0 / (3.0 * h));
				x = a / (a + b * std::exp(2.0 * w));
			}
			else
			{
				Real lna = std::log(a / (a + b));
				Real lnb = std::log(b / (a + b));
				t = std::exp(a * lna) / a;
				Real u = std::exp(b * lnb) / b;
				w = t + u;
				if (p < t / w) x = std::pow(a * w * p, 1.0 / a);
				else x = 1.0 - std::pow(b * w * (1.0 - p), 1.0 / b);
			}

			const Real afac = -std::lgamma(a) - std::lgamma(b) + std::lgamma(a + b);
			for (int j = 0; j < 12; j++)
			{
				if (x <= 0.0 || x >= 1.0) return (x <= 0.0) ? 0.0 : 1.0;
				Real err = RegularizedBetaI(x, a, b) - p;
				Real deriv = std::exp(a1 * std::log(x) + b1 * std::log(1.0 - x) + afac);   // dI/dx
				Real u = err / deriv;
				Real correction = u / (1.0 - 0.5 * std::min<Real>(REAL(1.0), u * (a1 / x - b1 / (REAL(1.0) - x))));
				x -= correction;
				if (x <= 0.0) x = 0.5 * (x + correction);
				if (x >= 1.0) x = 0.5 * (x + correction + 1.0);
				if (std::abs(correction) < REAL(1e-14) * x) break;
			}
			return x;
		}

		///////////////////////////////////////////////////////////////////////
		///                    Exponential integrals                        ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Exponential integral E_n(x) = ∫_1^∞ e^{-x t} / t^n dt.
		/// @details Continued fraction for x > 1, power series otherwise (Numerical Recipes).
		/// @throws DomainError for n < 0, x < 0, or x == 0 with n in {0,1}.
		inline Real ExpIntegralEn(int n, Real x)
		{
			static const Real EulerGamma = REAL(0.5772156649015328606);
			const Real eps = std::numeric_limits<Real>::epsilon();
			const Real big = std::numeric_limits<Real>::max() * eps;
			const int MAXIT = 100;
			const int nm1 = n - 1;

			if (n < 0 || x < 0.0 || (x == 0.0 && (n == 0 || n == 1)))
				throw DomainError("SpecialFunctions::ExpIntegralEn - bad arguments");

			if (n == 0) return std::exp(-x) / x;
			if (x == 0.0) return 1.0 / nm1;

			if (x > 1.0)
			{
				Real b = x + n;
				Real c = big;
				Real d = 1.0 / b;
				Real h = d;
				for (int i = 1; i <= MAXIT; i++)
				{
					Real a = -Real(i) * (nm1 + i);
					b += 2.0;
					d = 1.0 / (a * d + b);
					c = b + a / c;
					Real del = c * d;
					h *= del;
					if (std::abs(del - 1.0) <= eps)
						return h * std::exp(-x);
				}
				return h * std::exp(-x);
			}

			Real ans = (nm1 != 0) ? 1.0 / nm1 : -std::log(x) - EulerGamma;
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
					for (int ii = 1; ii <= nm1; ii++) psi += 1.0 / ii;
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
		/// @details Series for moderate x > 0, asymptotic for large x; Ei(x) = -E_1(-x) for x < 0.
		/// @throws DomainError if x == 0.
		inline Real ExpIntegralEi(Real x)
		{
			static const Real EulerGamma = REAL(0.5772156649015328606);
			const Real eps = std::numeric_limits<Real>::epsilon();
			const Real fpmin = std::numeric_limits<Real>::min() / eps;
			const int MAXIT = 100;

			if (x == 0.0) throw DomainError("SpecialFunctions::ExpIntegralEi - x must be nonzero");
			if (x < 0.0) return -ExpIntegralEn(1, -x);
			if (x < fpmin) return std::log(x) + EulerGamma;

			if (x <= -std::log(eps))
			{
				Real sum = 0.0, fact = 1.0;
				for (int k = 1; k <= MAXIT; k++)
				{
					fact *= x / k;
					Real term = fact / k;
					sum += term;
					if (term < eps * sum) break;
				}
				return sum + std::log(x) + EulerGamma;
			}

			Real sum = 0.0, term = 1.0;
			for (int k = 1; k <= MAXIT; k++)
			{
				Real prev = term;
				term *= Real(k) / x;
				if (term < eps) break;
				if (term < prev) sum += term;
				else { sum -= prev; break; }
			}
			return std::exp(x) * (1.0 + sum) / x;
		}

		///////////////////////////////////////////////////////////////////////
		///                   Sine and cosine integrals                     ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Computes both Si(x) = ∫_0^x sin t/t dt and Ci(x) = γ + ln|x| + ∫_0^x (cos t − 1)/t dt.
		/// @details Power series for |x| ≤ 2, complex continued fraction otherwise (Numerical Recipes).
		inline void SinCosIntegral(Real x, Real& si, Real& ci)
		{
			static const Real EulerGamma = REAL(0.5772156649015328606);
			const Real PIBY2 = Constants::PI / 2.0;
			const Real eps = std::numeric_limits<Real>::epsilon();
			const Real fpmin = std::numeric_limits<Real>::min() * REAL(4.0);
			const Real big = std::numeric_limits<Real>::max() * eps;
			const Real TMIN = REAL(2.0);
			const int MAXIT = 100;

			Real t = std::abs(x);
			if (t == 0.0) { si = 0.0; ci = -big; return; }

			if (t > TMIN)
			{
				std::complex<Real> b(1.0, t), c(big, 0.0), d, h, del;
				d = h = REAL(1.0) / b;
				for (int i = 1; i < MAXIT; i++)
				{
					Real a = -Real(i) * Real(i);
					b += std::complex<Real>(2.0, 0.0);
					d = REAL(1.0) / (a * d + b);
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

		/// @brief Sine integral Si(x).
		inline Real SineIntegral(Real x) { Real si, ci; SinCosIntegral(x, si, ci); return si; }

		/// @brief Cosine integral Ci(x). Ci(0) = -∞.
		inline Real CosineIntegral(Real x) { Real si, ci; SinCosIntegral(x, si, ci); return ci; }

		///////////////////////////////////////////////////////////////////////
		///                   Dawson and Fresnel integrals                  ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Dawson's integral D(x) = e^{-x²} ∫_0^x e^{t²} dt.
		/// @details Maclaurin series for |x| < 0.2, sampling-theorem evaluation otherwise (Numerical Recipes).
		inline Real Dawson(Real x)
		{
			const int NMAX = 6;
			const Real H = REAL(0.4);
			const Real A1 = REAL(2.0) / REAL(3.0), A2 = REAL(0.4), A3 = REAL(2.0) / REAL(7.0);
			static const std::array<Real, NMAX> ctab = [] {
				std::array<Real, NMAX> c{};
				for (int i = 0; i < NMAX; i++) { Real t = (2.0 * i + 1.0) * REAL(0.4); c[i] = std::exp(-t * t); }
				return c;
			}();

			if (std::abs(x) < REAL(0.2))
			{
				Real x2 = x * x;
				return x * (1.0 - A1 * x2 * (1.0 - A2 * x2 * (1.0 - A3 * x2)));
			}

			Real xx = std::abs(x);
			int n0 = 2 * int(0.5 * xx / H + 0.5);
			Real xp = xx - n0 * H;
			Real e1 = std::exp(2.0 * xp * H);
			Real e2 = e1 * e1;
			Real d1 = n0 + 1;
			Real d2 = d1 - 2.0;
			Real sum = 0.0;
			for (int i = 0; i < NMAX; i++, d1 += 2.0, d2 -= 2.0, e1 *= e2)
				sum += ctab[i] * (e1 / d1 + 1.0 / (d2 * e1));

			Real val = (REAL(1.0) / std::sqrt(Constants::PI)) * std::exp(-xp * xp) * sum;
			return (x >= 0.0) ? val : -val;
		}

		/// @brief Computes both Fresnel integrals S(x) = ∫_0^x sin(π t²/2) dt and C(x) = ∫_0^x cos(π t²/2) dt.
		/// @details Power series for |x| ≤ 1.5, complex continued fraction otherwise (Numerical Recipes).
		inline void FresnelSC(Real x, Real& s, Real& c)
		{
			const Real eps = std::numeric_limits<Real>::epsilon();
			const Real fpmin = std::numeric_limits<Real>::min();
			const Real big = std::numeric_limits<Real>::max() * eps;
			const Real XMIN = REAL(1.5);
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
				d = h = REAL(1.0) / b;
				int n = -1;
				for (int k = 2; k <= MAXIT; k++)
				{
					n += 2;
					Real a = -Real(n) * Real(n + 1);
					b += std::complex<Real>(4.0, 0.0);
					d = REAL(1.0) / (a * d + b);
					cc = b + a / cc;
					del = cc * d;
					h *= del;
					if (std::abs(std::real(del) - 1.0) + std::abs(std::imag(del)) <= eps) break;
				}
				h *= std::complex<Real>(ax, -ax);
				cs = std::complex<Real>(0.5, 0.5) * (REAL(1.0) - std::complex<Real>(std::cos(0.5 * pix2), std::sin(0.5 * pix2)) * h);
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

		/// @brief Principal branch W_0(x) of the Lambert W function, x ≥ -1/e.
		/// @details Solves w·e^w = x. Series near the branch point, log approximation for
		///          large x, refined by Halley iteration.
		/// @throws DomainError if x < -1/e.
		inline Real LambertW0(Real x)
		{
			const Real EE = REAL(2.718281828459045235360287);
			const Real minX = REAL(-0.36787944117144232159553);  // -1/e
			const Real eps = std::numeric_limits<Real>::epsilon();

			if (x < minX - REAL(1e-12)) throw DomainError("SpecialFunctions::LambertW0 - x must be >= -1/e");
			if (x <= minX + REAL(1e-10)) return REAL(-1.0);   // at/near the branch point -1/e
			if (x == 0.0) return 0.0;

			Real w;
			if (x < REAL(-0.3))
			{
				Real arg = REAL(2.0) * (EE * x + REAL(1.0));
				Real p = std::sqrt(arg < 0.0 ? REAL(0.0) : arg);
				w = REAL(-1.0) + p - p * p / REAL(3.0) + REAL(11.0) * p * p * p / REAL(72.0);
			}
			else if (x < REAL(1.5))
				w = x / (REAL(1.0) + x);   // simple, Halley refines
			else
			{
				Real l1 = std::log(x);
				Real l2 = std::log(l1);
				w = l1 - l2 + l2 / l1;
			}

			for (int i = 0; i < 100; i++)
			{
				Real ew = std::exp(w);
				Real f = w * ew - x;
				Real wp1 = w + REAL(1.0);
				Real dw = f / (ew * wp1 - (w + REAL(2.0)) * f / (REAL(2.0) * wp1));
				w -= dw;
				if (std::abs(dw) <= eps * (REAL(1.0) + std::abs(w))) break;
			}
			return w;
		}

		/// @brief Secondary branch W_{-1}(x) of the Lambert W function, -1/e ≤ x < 0.
		/// @throws DomainError if x is outside [-1/e, 0).
		inline Real LambertWm1(Real x)
		{
			const Real EE = REAL(2.718281828459045235360287);
			const Real minX = REAL(-0.36787944117144232159553);  // -1/e
			const Real eps = std::numeric_limits<Real>::epsilon();

			if (x >= 0.0 || x < minX - REAL(1e-12)) throw DomainError("SpecialFunctions::LambertWm1 - x must be in [-1/e, 0)");
			if (x <= minX + REAL(1e-10)) return REAL(-1.0);   // at/near the branch point -1/e

			Real w;
			if (x < REAL(-0.3))
			{
				Real arg = REAL(2.0) * (EE * x + REAL(1.0));
				Real p = -std::sqrt(arg < 0.0 ? REAL(0.0) : arg);  // negative branch
				w = REAL(-1.0) + p - p * p / REAL(3.0) + REAL(11.0) * p * p * p / REAL(72.0);
			}
			else
			{
				Real l1 = std::log(-x);
				Real l2 = std::log(-l1);
				w = l1 - l2 + l2 / l1;
			}

			for (int i = 0; i < 100; i++)
			{
				Real ew = std::exp(w);
				Real f = w * ew - x;
				Real wp1 = w + REAL(1.0);
				Real dw = f / (ew * wp1 - (w + REAL(2.0)) * f / (REAL(2.0) * wp1));
				w -= dw;
				if (std::abs(dw) <= eps * (REAL(1.0) + std::abs(w))) break;
			}
			return w;
		}

		///////////////////////////////////////////////////////////////////////
		///                        Airy functions                           ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Computes all four Airy functions Ai(x), Bi(x), Ai'(x), Bi'(x) at once.
		/// @details Maclaurin series for |x| < 6, DLMF asymptotic expansions beyond.
		inline void Airy(Real x, Real& ai, Real& bi, Real& aip, Real& bip)
		{
			const Real sqrtPi = std::sqrt(Constants::PI);
			const Real c1 = REAL(0.3550280538878172);   // Ai(0)
			const Real c2 = REAL(0.2588194037928068);   // -Ai'(0)
			const Real eps = std::numeric_limits<Real>::epsilon();

			if (std::abs(x) < REAL(6.0))
			{
				const Real sqrt3 = std::sqrt(REAL(3.0));
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
			const Real zeta = (REAL(2.0) / REAL(3.0)) * std::pow(t, REAL(1.5));
			const Real t14 = std::pow(t, REAL(0.25));
			const Real zinv = REAL(1.0) / zeta;

			if (x > 0.0)
			{
				Real Sc = 0, ScAlt = 0, Sv = 0, SvAlt = 0;
				Real cm = 1.0, zpow = 1.0, prev = REAL(1e300);
				for (int m = 0; m <= 40; m++)
				{
					Real vm = (m == 0) ? REAL(1.0) : -cm * Real(6 * m + 1) / Real(6 * m - 1);
					Real termC = cm * zpow, termV = vm * zpow;
					if (m > 0 && std::abs(termC) > prev) break;   // asymptotic divergence
					Real sgn = (m % 2 == 0) ? REAL(1.0) : REAL(-1.0);
					Sc += termC; ScAlt += sgn * termC;
					Sv += termV; SvAlt += sgn * termV;
					prev = std::abs(termC);
					if (std::abs(termC) < eps) break;
					int k = m + 1;
					cm *= Real(6 * k - 5) * Real(6 * k - 3) * Real(6 * k - 1) / (Real(2 * k - 1) * REAL(216.0) * Real(k));
					zpow *= zinv;
				}
				Real emz = std::exp(-zeta), epz = std::exp(zeta);
				ai = emz / (REAL(2.0) * sqrtPi * t14) * ScAlt;
				bi = epz / (sqrtPi * t14) * Sc;
				aip = -t14 * emz / (REAL(2.0) * sqrtPi) * SvAlt;
				bip = t14 * epz / sqrtPi * Sv;
				return;
			}

			// x < 0
			const Real theta = zeta - Constants::PI / REAL(4.0);
			const Real cth = std::cos(theta), sth = std::sin(theta);
			Real P = 0, Q = 0, R = 0, S = 0;
			Real cm = 1.0, zpow = 1.0, prev = REAL(1e300);
			for (int m = 0; m <= 40; m++)
			{
				Real vm = (m == 0) ? REAL(1.0) : -cm * Real(6 * m + 1) / Real(6 * m - 1);
				Real termC = cm * zpow, termV = vm * zpow;
				if (m > 0 && std::abs(termC) > prev) break;
				Real sgn = ((m / 2) % 2 == 0) ? REAL(1.0) : REAL(-1.0);
				if (m % 2 == 0) { P += sgn * termC; R += sgn * termV; }
				else { Q += sgn * termC; S += sgn * termV; }
				prev = std::abs(termC);
				if (std::abs(termC) < eps) break;
				int k = m + 1;
				cm *= Real(6 * k - 5) * Real(6 * k - 3) * Real(6 * k - 1) / (Real(2 * k - 1) * REAL(216.0) * Real(k));
				zpow *= zinv;
			}
			Real inv = REAL(1.0) / (sqrtPi * t14);
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
				Real u = REAL(1.0) / (t * t);
				return std::pow(t, REAL(2.0) / REAL(3.0)) *
					(REAL(1.0) + u * (REAL(5.0) / 48 - u * (REAL(5.0) / 36 - u * (REAL(77125.0) / 82944))));
			}
		}

		/// @brief The n-th real (negative) zero a_n of Ai (n ≥ 1).
		inline Real AiryAiZero(int n)
		{
			if (n < 1) throw DomainError("SpecialFunctions::AiryAiZero - n must be >= 1");
			Real t = REAL(3.0) * Constants::PI * (REAL(4.0) * n - REAL(1.0)) / REAL(8.0);
			Real x = -detail::AiryZeroT(t);
			for (int i = 0; i < 60; i++)
			{
				Real ai, bi, aip, bip; Airy(x, ai, bi, aip, bip);
				Real dx = ai / aip;
				x -= dx;
				if (std::abs(dx) <= REAL(1e-14) * (REAL(1.0) + std::abs(x))) break;
			}
			return x;
		}

		/// @brief The n-th real (negative) zero b_n of Bi (n ≥ 1).
		inline Real AiryBiZero(int n)
		{
			if (n < 1) throw DomainError("SpecialFunctions::AiryBiZero - n must be >= 1");
			Real t = REAL(3.0) * Constants::PI * (REAL(4.0) * n - REAL(3.0)) / REAL(8.0);
			Real x = -detail::AiryZeroT(t);
			for (int i = 0; i < 60; i++)
			{
				Real ai, bi, aip, bip; Airy(x, ai, bi, aip, bip);
				Real dx = bi / bip;
				x -= dx;
				if (std::abs(dx) <= REAL(1e-14) * (REAL(1.0) + std::abs(x))) break;
			}
			return x;
		}
	} // namespace SpecialFunctions
}
#endif
