///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        StandardFunctions.h                                                 ///
///  Description: Standard mathematical functions (Bessel, Gamma, Legendre, etc.)     ///
///               Special functions from C++ math spec                                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1

#if !defined  MML_FUNCTIONS_H
#define MML_FUNCTIONS_H

#include "MMLBase.h"

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
		static inline T Sec(T x) { return T{1} / std::cos(x); }
		template<typename T>
		static inline T Csc(T x) { return T{1} / std::sin(x); }
		template<typename T>
		static inline T Tan(T x) { return std::tan(x); }
		template<typename T>
		static inline T Ctg(T x) { return T{1} / std::tan(x); }

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
		static inline T Csch(T x) { return T{1} / std::sinh(x); }
		template<typename T>
		static inline T Tanh(T x) { return std::tanh(x); }
		template<typename T>
		static inline T Ctgh(T x) { return T{1} / std::tanh(x); }

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

    // Note: Basic math functions (Sin, Cos, Exp, Log, etc.) now templated above.
    // They work with Real, Complex, float, long double, etc. - no separate Complex overloads needed.

		// Factorial functions
		static inline Real Factorial(int n) {
			Real fact = 1.0;
			for (int i = 2; i <= n; i++)
				fact *= i;
			return fact;
		}
		static inline long long FactorialInt(int n) {
			long long fact = 1;
			for (int i = 2; i <= n; i++)
				fact *= i;
			return fact;
		}
		static inline Real FactorialStirling(int n)
		{
			if (n < 0) return 0.0; // Factorial is not defined for negative integers
			if (n == 0 || n == 1) return 1.0;
			return sqrt(2 * Constants::PI * n) * std::pow(n / Constants::E, n);
		}
	}
}

#endif