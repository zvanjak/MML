///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        StandardFunctions.h                                                 ///
///  Description: Standard mathematical functions (Bessel, Gamma, Legendre, etc.)     ///
///               Special functions from C++ math spec                                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#if !defined  MML_FUNCTIONS_H
#define MML_FUNCTIONS_H

#include "MMLBase.h"

// Feature detection for C++17 special math functions
// macOS libc++ doesn't support these yet, so we need fallbacks
#if defined(__cpp_lib_math_special_functions) || \
    (defined(_MSC_VER) && _MSC_VER >= 1910) || \
    (defined(__GNUC__) && !defined(__clang__) && __GNUC__ >= 7) || \
    (defined(__clang__) && !defined(__APPLE__) && __clang_major__ >= 10)
    #define MML_HAS_STD_SPECIAL_FUNCTIONS 1
#else
    #define MML_HAS_STD_SPECIAL_FUNCTIONS 0
#endif

namespace MML
{
	namespace Functions
	{
		// Basic math functions - templated to support all numeric types (Real, Complex, float, long double, etc.)
		// These delegate to std:: functions which handle type-specific implementations
		template<typename T>
		static inline T Sin(T x) { return std::sin(x); }
		template<typename T>
		static inline T Cos(T x) { return std::cos(x); }
		template<typename T>
		static inline T Sec(T x) { T c = std::cos(x); if (std::abs(c) < std::numeric_limits<T>::epsilon()) throw std::domain_error("Sec: cos(x) is zero"); return T{1} / c; }
		template<typename T>
		static inline T Csc(T x) { T s = std::sin(x); if (std::abs(s) < std::numeric_limits<T>::epsilon()) throw std::domain_error("Csc: sin(x) is zero"); return T{1} / s; }
		template<typename T>
		static inline T Tan(T x) { return std::tan(x); }
		template<typename T>
		static inline T Ctg(T x) { T t = std::tan(x); if (std::abs(t) < std::numeric_limits<T>::epsilon()) throw std::domain_error("Ctg: tan(x) is zero"); return T{1} / t; }

		template<typename T>
		static inline T Exp(T x) { return std::exp(x); }
		template<typename T>
		static inline T Log(T x) { return std::log(x); }
		template<typename T>
		static inline T Log10(T x) { return std::log10(x); }
		template<typename T>
		static inline T Sqrt(T x) { return std::sqrt(x); }
		template<typename T>
		static inline T Pow(T x, T y) { return std::pow(x, y); }

		template<typename T>
		static inline T Sinh(T x) { return std::sinh(x); }
		template<typename T>
		static inline T Cosh(T x) { return std::cosh(x); }
		template<typename T>
		static inline T Sech(T x) { return T{1} / std::cosh(x); }
		template<typename T>
		static inline T Csch(T x) { T s = std::sinh(x); if (std::abs(s) < std::numeric_limits<T>::epsilon()) throw std::domain_error("Csch: sinh(x) is zero"); return T{1} / s; }
		template<typename T>
		static inline T Tanh(T x) { return std::tanh(x); }
		template<typename T>
		static inline T Ctgh(T x) { T t = std::tanh(x); if (std::abs(t) < std::numeric_limits<T>::epsilon()) throw std::domain_error("Ctgh: tanh(x) is zero"); return T{1} / t; }

		template<typename T>
		static inline T Asin(T x) { return std::asin(x); }
		template<typename T>
		static inline T Acos(T x) { return std::acos(x); }
		template<typename T>
		static inline T Atan(T x) { return std::atan(x); }

		template<typename T>
		static inline T Asinh(T x) { return std::asinh(x); }
		template<typename T>
		static inline T Acosh(T x) { return std::acosh(x); }
		template<typename T>
		static inline T Atanh(T x) { return std::atanh(x); }

		// Special functions - Real only (std:: library limitation)
		static inline Real Erf(Real x) { return std::erf(x); }
		static inline Real Erfc(Real x) { return std::erfc(x); }

		static inline Real TGamma(Real x) { return std::tgamma(x); }
		static inline Real LGamma(Real x) { return std::lgamma(x); }

#if MML_HAS_STD_SPECIAL_FUNCTIONS
		// C++17 special math functions (not available on macOS libc++)
		// NOTE: These are ONLY wrappers around std:: functions. They are NOT used by MML internally.
		//       MML uses its own implementations in mml/algorithms/ which work on all platforms.
		static inline Real RiemannZeta(Real x) { return std::riemann_zeta(x); }

		static inline Real Hermite(unsigned int n, Real x) { return std::hermite(n, x); }
		static inline Real Legendre(unsigned int n, Real x) { return std::legendre(n, x); }
		static inline Real Laguerre(unsigned int n, Real x) { return std::laguerre(n, x); }
		static inline Real SphBessel(unsigned int n, Real x) { return std::sph_bessel(n, x); }
		static inline Real SphLegendre(int n1, int n2, Real x) { return std::sph_legendre(n1, n2, x); }

		// Elliptic integrals (incomplete and complete)
		static inline Real Ellint_1(Real k, Real phi) { return std::ellint_1(k, phi); }
		static inline Real Ellint_2(Real k, Real phi) { return std::ellint_2(k, phi); }
		static inline Real Ellint_3(Real k, Real n, Real phi) { return std::ellint_3(k, n, phi); }
		static inline Real Comp_ellint_1(Real x) { return std::comp_ellint_1(x); }
		static inline Real Comp_ellint_2(Real x) { return std::comp_ellint_2(x); }
		static inline Real Comp_ellint_3(Real k, Real n) { return std::comp_ellint_3(k, n); }

		// Bessel functions
		static inline Real CylBesselJ(Real n, Real x) { return std::cyl_bessel_j(n, x); }
		static inline Real CylBesselI(Real n, Real x) { return std::cyl_bessel_i(n, x); }
		static inline Real CylBesselK(Real n, Real x) { return std::cyl_bessel_k(n, x); }
		static inline Real CylNeumann(Real n, Real x) { return std::cyl_neumann(n, x); }

		// Spherical Neumann function
		static inline Real SphNeumann(unsigned int n, Real x) { return std::sph_neumann(n, x); }

		// Associated Legendre and Laguerre polynomials
		static inline Real AssocLegendre(unsigned int l, unsigned int m, Real x) { return std::assoc_legendre(l, m, x); }
		static inline Real AssocLaguerre(unsigned int n, unsigned int m, Real x) { return std::assoc_laguerre(n, m, x); }

		// Beta and exponential integral functions
		static inline Real Beta(Real x, Real y) { return std::beta(x, y); }
		static inline Real Expint(Real x) { return std::expint(x); }

		// Alias for Bessel function of the second kind (same as CylNeumann)
		static inline Real CylBesselY(Real n, Real x) { return std::cyl_neumann(n, x); }
#else
		// =====================================================================================
		// FALLBACK IMPLEMENTATIONS for platforms without C++17 special math (e.g., macOS libc++)
		// These provide the same interface as the std:: functions above.
		// =====================================================================================

		// Legendre polynomial Pₙ(x) using Bonnet's recurrence
		static inline Real Legendre(unsigned int n, Real x) {
			if (n == 0) return 1.0;
			if (n == 1) return x;
			Real P_prev2 = 1.0, P_prev1 = x, P_n = 0.0;
			for (unsigned int k = 2; k <= n; ++k) {
				P_n = ((2.0 * k - 1.0) * x * P_prev1 - (k - 1.0) * P_prev2) / k;
				P_prev2 = P_prev1;
				P_prev1 = P_n;
			}
			return P_n;
		}

		// Hermite polynomial Hₙ(x) (physicist's convention) using recurrence
		static inline Real Hermite(unsigned int n, Real x) {
			if (n == 0) return 1.0;
			if (n == 1) return 2.0 * x;
			Real H_prev2 = 1.0, H_prev1 = 2.0 * x, H_n = 0.0;
			for (unsigned int k = 2; k <= n; ++k) {
				H_n = 2.0 * x * H_prev1 - 2.0 * (k - 1.0) * H_prev2;
				H_prev2 = H_prev1;
				H_prev1 = H_n;
			}
			return H_n;
		}

		// Laguerre polynomial Lₙ(x) using recurrence
		static inline Real Laguerre(unsigned int n, Real x) {
			if (n == 0) return 1.0;
			if (n == 1) return 1.0 - x;
			Real L_prev2 = 1.0, L_prev1 = 1.0 - x, L_n = 0.0;
			for (unsigned int k = 2; k <= n; ++k) {
				L_n = ((2.0 * k - 1.0 - x) * L_prev1 - (k - 1.0) * L_prev2) / k;
				L_prev2 = L_prev1;
				L_prev1 = L_n;
			}
			return L_n;
		}

		// Associated Legendre polynomial Pₗᵐ(x) using recurrence
		static inline Real AssocLegendre(unsigned int l, unsigned int m, Real x) {
			if (m > l) return 0.0;
			// Start with P_m^m
			Real P_mm = 1.0;
			if (m > 0) {
				Real sqrt1mx2 = std::sqrt(1.0 - x * x);
				Real fact = 1.0;
				for (unsigned int i = 1; i <= m; ++i) {
					P_mm *= -fact * sqrt1mx2;
					fact += 2.0;
				}
			}
			if (l == m) return P_mm;
			// P_{m+1}^m
			Real P_mp1_m = x * (2.0 * m + 1.0) * P_mm;
			if (l == m + 1) return P_mp1_m;
			// Recurrence for higher l
			Real P_l_m = 0.0;
			for (unsigned int ll = m + 2; ll <= l; ++ll) {
				P_l_m = (x * (2.0 * ll - 1.0) * P_mp1_m - (ll + m - 1.0) * P_mm) / (ll - m);
				P_mm = P_mp1_m;
				P_mp1_m = P_l_m;
			}
			return P_l_m;
		}

		// Associated Laguerre polynomial Lₙᵐ(x) using recurrence
		static inline Real AssocLaguerre(unsigned int n, unsigned int m, Real x) {
			if (n == 0) return 1.0;
			if (n == 1) return 1.0 + m - x;
			Real L_prev2 = 1.0;
			Real L_prev1 = 1.0 + m - x;
			Real L_n = 0.0;
			for (unsigned int k = 2; k <= n; ++k) {
				L_n = ((2.0 * k - 1.0 + m - x) * L_prev1 - (k - 1.0 + m) * L_prev2) / k;
				L_prev2 = L_prev1;
				L_prev1 = L_n;
			}
			return L_n;
		}

		// Spherical Legendre function Y_l^m(theta) (real, normalized)
		static inline Real SphLegendre(unsigned int l, unsigned int m, Real theta) {
			// Y_l^m(theta) = K_l^m * P_l^m(cos(theta))
			// Normalization: K_l^m = sqrt((2l+1)/(4π) * (l-m)!/(l+m)!)
			Real x = std::cos(theta);
			Real P_lm = AssocLegendre(l, m, x);
			// Compute (l-m)! / (l+m)!
			Real ratio = 1.0;
			for (unsigned int i = l - m + 1; i <= l + m; ++i)
				ratio /= i;
			Real K = std::sqrt((2.0 * l + 1.0) / (4.0 * Constants::PI) * ratio);
			return K * P_lm;
		}

		// Complete elliptic integral of the first kind K(k) using AGM
		static inline Real Comp_ellint_1(Real k) {
			Real a = 1.0, b = std::sqrt(1.0 - k * k);
			while (std::abs(a - b) > 1e-15) {
				Real temp = (a + b) / 2.0;
				b = std::sqrt(a * b);
				a = temp;
			}
			return Constants::PI / (2.0 * a);
		}

		// Complete elliptic integral of the second kind E(k) using AGM
		static inline Real Comp_ellint_2(Real k) {
			Real a = 1.0, b = std::sqrt(1.0 - k * k);
			Real c = k, sum = k * k / 2.0;
			Real power_of_2 = 1.0;
			while (std::abs(c) > 1e-15) {
				Real a_new = (a + b) / 2.0;
				Real b_new = std::sqrt(a * b);
				c = (a - b) / 2.0;
				power_of_2 *= 2.0;
				sum += power_of_2 * c * c;
				a = a_new;
				b = b_new;
			}
			return Constants::PI / (2.0 * a) * (1.0 - sum);
		}

		// Note: More special functions (Bessel, Beta, etc.) can be added here if needed.
		// For now, use MML's own implementations from mml/algorithms/ for those.
#endif

		// Factorial functions
		static inline Real Factorial(int n) {
			if (n < 0) throw std::domain_error("Factorial: negative argument");
			Real fact = 1.0;
			for (int i = 2; i <= n; i++)
				fact *= i;
			return fact;
		}
		static inline long long FactorialInt(int n) {
			if (n < 0) throw std::domain_error("FactorialInt: negative argument");
			long long fact = 1;
			for (int i = 2; i <= n; i++)
				fact *= i;
			return fact;
		}
		static inline Real FactorialStirling(int n)
		{
			if (n < 0) throw std::domain_error("FactorialStirling: negative argument");
			if (n == 0 || n == 1) return 1.0;
			return sqrt(2 * Constants::PI * n) * std::pow(n / Constants::E, n);
		}
	}
}

#endif