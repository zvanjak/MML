///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ComparisonUtils.h                                                   ///
///  Description: Equality comparison utilities for Complex, Vector types             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_COMPARISON_UTILS_H
#define MML_COMPARISON_UTILS_H

#include <cmath>
#include <complex>

#include "mml/MMLBase.h"
#include "mml/base/Vector/Vector.h"

namespace MML
{
	namespace Utils
	{
		// ============================================================================
		// Complex Number Comparisons
		// ============================================================================

		/// @brief Tests if two complex numbers are equal within tolerance (per component).
		/// @param a First complex number
		/// @param b Second complex number
		/// @param eps Tolerance for real and imaginary parts separately
		/// @return True if both real and imaginary differences are within tolerance
		static bool AreEqual(const Complex& a, const Complex& b, double eps = Defaults::ComplexAreEqualTolerance)
		{
			if (std::abs(a.real() - b.real()) > eps || std::abs(a.imag() - b.imag()) > eps)
				return false;
			return true;
		}

		/// @brief Tests if two complex numbers are equal within tolerance (absolute value).
		/// @param a First complex number
		/// @param b Second complex number
		/// @param eps Tolerance for |a - b|
		/// @return True if |a - b| <= eps
		static bool AreEqualAbs(const Complex& a, const Complex& b, double eps = Defaults::ComplexAreEqualAbsTolerance)
		{
			if (Abs(a - b) > eps)
				return false;
			return true;
		}

		// ============================================================================
		// Real Vector Comparisons
		// ============================================================================

		/// @brief Tests if two real vectors are equal within tolerance.
		/// @param a First vector
		/// @param b Second vector
		/// @param eps Tolerance for element-wise comparison
		/// @return True if vectors are equal within tolerance
		static bool AreEqual(const Vector<Real> &a, const Vector<Real> &b, 
		                     Real eps = Defaults::VectorIsEqualTolerance)
		{
			return a.IsEqualTo(b, eps);
		}

		// ============================================================================
		// Complex Vector Comparisons
		// ============================================================================

		/// @brief Tests if two complex vectors are equal within tolerance (per component).
		/// @param a First complex vector
		/// @param b Second complex vector
		/// @param eps Tolerance for each complex element
		/// @return True if all elements are equal within tolerance
		static bool AreEqual(const Vector<Complex>& a, const Vector<Complex>& b, 
		                     Real eps = Defaults::ComplexAreEqualTolerance)
		{
			if (a.size() != b.size())
				return false;

			for (int i = 0; i < a.size(); i++)
				if( !AreEqual(a[i], b[i], eps) )
					return false;

			return true;
		}

		/// @brief Tests if two complex vectors are equal within tolerance (absolute value).
		/// @param a First complex vector
		/// @param b Second complex vector
		/// @param eps Tolerance for |a[i] - b[i]|
		/// @return True if all elements differ by at most eps in absolute value
		static bool AreEqualAbs(const Vector<Complex>& a, const Vector<Complex>& b, 
		                        Real eps = Defaults::ComplexAreEqualAbsTolerance)
		{
			if (a.size() != b.size())
				return false;

			for (int i = 0; i < a.size(); i++)
				if ( !AreEqualAbs(a[i], b[i], eps) )
					return false;

			return true;
		}

	} // namespace Utils
} // namespace MML

#endif // MML_COMPARISON_UTILS_H
