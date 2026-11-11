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
#if !defined  MML_FUNCTIONS_H
#define MML_FUNCTIONS_H

#include <mml/MMLBase.h>
#include <mml/base/ExactInteger.h>

#include <cmath>

// Detect the standard library's C++17 special math support, not the compiler's version.
#ifndef MML_HAS_STD_SPECIAL_FUNCTIONS
#if defined(__cpp_lib_math_special_functions) && __cpp_lib_math_special_functions >= 201603L
    #define MML_HAS_STD_SPECIAL_FUNCTIONS 1
#else
    #define MML_HAS_STD_SPECIAL_FUNCTIONS 0
#endif
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
		static inline T Sec(T x) { T c = std::cos(x); if (std::abs(c) < std::numeric_limits<decltype(std::abs(c))>::epsilon()) throw DomainError("Sec: cos(x) is zero"); return T{1} / c; }
		template<typename T>
		static inline T Csc(T x) { T s = std::sin(x); if (std::abs(s) < std::numeric_limits<decltype(std::abs(s))>::epsilon()) throw DomainError("Csc: sin(x) is zero"); return T{1} / s; }
		template<typename T>
		static inline T Tan(T x) { return std::tan(x); }
		template<typename T>
		static inline T Ctg(T x) { T t = std::tan(x); if (std::abs(t) < std::numeric_limits<decltype(std::abs(t))>::epsilon()) throw DomainError("Ctg: tan(x) is zero"); return T{1} / t; }

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
		static inline T Csch(T x) { T s = std::sinh(x); if (std::abs(s) < std::numeric_limits<decltype(std::abs(s))>::epsilon()) throw DomainError("Csch: sinh(x) is zero"); return T{1} / s; }
		template<typename T>
		static inline T Tanh(T x) { return std::tanh(x); }
		template<typename T>
		static inline T Ctgh(T x) { T t = std::tanh(x); if (std::abs(t) < std::numeric_limits<decltype(std::abs(t))>::epsilon()) throw DomainError("Ctgh: tanh(x) is zero"); return T{1} / t; }

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
		// C++17 special math functions supplied by the selected standard library.
		// Complete elliptic E uses the shared portable implementation below.
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
		static inline Real Comp_ellint_3(Real k, Real n) { return std::comp_ellint_3(k, n); }

		// Bessel functions
		static inline Real CylBesselJ(Real n, Real x) { return std::cyl_bessel_j(n, x); }
		static inline Real CylBesselI(Real n, Real x) { return std::cyl_bessel_i(n, x); }
		static inline Real CylBesselK(Real n, Real x) { return std::cyl_bessel_k(n, x); }
		static inline Real CylNeumann(Real n, Real x) { return std::cyl_neumann(n, x); }

		// Spherical Neumann function
		static inline Real SphNeumann(unsigned int n, Real x) { return std::sph_neumann(n, x); }

		// Associated Legendre and Laguerre polynomials
		static inline Real AssocLegendre(unsigned int l, unsigned int m, Real x) {
			if (std::isnan(x)) return x;
			if (std::abs(x) > 1) throw DomainError("AssocLegendre: |x| must not exceed 1");
			return std::assoc_legendre(l, m, x);
		}
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
			if (std::isnan(x)) return x;
			if (std::abs(x) > 1) throw DomainError("AssocLegendre: |x| must not exceed 1");
			if (m > l) return 0.0;
			// Start with P_m^m
			Real P_mm = 1.0;
			if (m > 0) {
				Real sqrt1mx2 = std::sqrt((1 - x) * (1 + x));
				Real fact = 1.0;
				for (unsigned int i = 1; i <= m; ++i) {
					P_mm *= fact * sqrt1mx2;
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
			if (std::isnan(theta)) return theta;
			if (m > l) return 0;
			// Unlike AssocLegendre, the spherical function includes the Condon-Shortley phase.
			// Normalization: K_l^m = sqrt((2l+1)/(4π) * (l-m)!/(l+m)!)
			Real x = std::cos(theta);
			Real P_lm = AssocLegendre(l, m, x);
			if (m % 2 != 0) P_lm = -P_lm;
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
			while (std::abs(a - b) > PrecisionValues<Real>::AGMConvergenceTolerance) {
				Real temp = (a + b) / 2.0;
				b = std::sqrt(a * b);
				a = temp;
			}
			return Constants::PI / (2.0 * a);
		}

		// Note: More special functions (Bessel, Beta, etc.) can be added here if needed.
		// For now, use MML's own implementations from mml/algorithms/ for those.
#endif

		// Shared AGM implementation: native libstdc++ E loses accuracy near |k|=1.
		static inline Real Comp_ellint_2(Real k) {
			if (std::isnan(k)) return k;
			if (std::abs(k) > 1) throw DomainError("Comp_ellint_2: |k| must not exceed 1");
			if (std::abs(k) == 1) return 1;
			Real a = 1, b = std::sqrt((1 - k) * (1 + k));
			Real sum = k * k / 2, power_of_2 = 1;
			// DLMF 19.8(i): E/K = 1 - k^2/2 - c_1^2 - 2*c_2^2 - 4*c_3^2 - ...
			for (int iteration = 0; iteration < Precision::EllipticAGMMaxIterations; ++iteration) {
				const Real c = (a - b) / 2;
				const Real a_new = (a + b) / 2;
				const Real b_new = std::sqrt(a * b);
				sum += power_of_2 * c * c;
				a = a_new;
				b = b_new;
				if (std::abs(c) <= Precision::AGMConvergenceTolerance * a)
					return Constants::PI / (2 * a) * (1 - sum);
				power_of_2 *= 2;
			}
			throw ConvergenceError("Comp_ellint_2: AGM iteration budget exhausted");
		}

		// Factorial functions
		/// @brief n! as Real; negative n throws DomainError, overflow returns positive infinity.
		static inline Real Factorial(int n) {
			if (n < 0) throw DomainError("Factorial: negative argument");
			Real fact = 1.0;
			for (int i = 2; i <= n; i++) {
				if (fact > std::numeric_limits<Real>::max() / i)
					return Constants::PosInf;
				fact *= i;
			}
			return fact;
		}
		/// @brief Exact factorial; throws DomainError for negative n or a result outside signed 64-bit range.
		static inline long long FactorialInt(int n) {
			if (n < 0) throw DomainError("FactorialInt: negative argument");
			long long fact = 1;
			for (int i = 2; i <= n; i++)
				fact = ExactInteger::Multiply(fact, i);
			return fact;
		}
		static inline Real FactorialStirling(int n)
		{
			if (n < 0) throw DomainError("FactorialStirling: negative argument");
			if (n == 0 || n == 1) return 1.0;
			return sqrt(2 * Constants::PI * n) * std::pow(n / Constants::E, Real(n));
		}
	}
}

#endif