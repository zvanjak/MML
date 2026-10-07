///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        InterpolatedRealFunctionPolynomial.h                                ///
///  Description: Polynomial interpolation for real-valued functions                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file InterpolatedRealFunctionPolynomial.h
/// @brief Polynomial interpolation using Neville's algorithm.
/// @ingroup Interpolation

#if !defined MML_INTERPOLATED_REAL_FUNCTION_POLYNOMIAL_H
#define MML_INTERPOLATED_REAL_FUNCTION_POLYNOMIAL_H

#include <mml/base/InterpolatedFunctions/InterpolatedRealFunctionLinear.h>

namespace MML {
	/////////////////////////////////////////////////////////////////////////////////////
	///                      POLYNOMIAL INTERPOLATION                                 ///
	/////////////////////////////////////////////////////////////////////////////////////

	/// @brief Polynomial interpolation using Neville's algorithm.
	/// PolynomInterpRealFunc constructs an interpolating polynomial of degree mm-1
	/// through mm consecutive data points. The polynomial passes exactly through
	/// all used points and provides smooth, infinitely differentiable interpolation.
	/// **Algorithm:** Neville's algorithm builds the interpolating polynomial
	/// iteratively, computing successively higher-degree approximations. It also
	/// provides an error estimate as a byproduct.
	/// **Advantages:**
	/// - Exact interpolation through all data points
	/// - Built-in error estimation
	/// - C^∞ smooth within each local region
	/// **Limitations:**
	/// - Runge phenomenon: High-degree polynomials can oscillate wildly
	/// - Not recommended for more than ~10 points
	/// - Error can be large between data points for high degrees
	/// @warning For many data points, use SplineInterpRealFunc or BarycentricRationalInterp instead.
	/// @see BarycentricRationalInterp for oscillation-free alternative
	/// @ingroup Interpolation

	class PolynomInterpRealFunc : public RealFunctionInterpolated {
	private:
		mutable Real _errorEst;

	public:
		/// @brief Construct a polynomial interpolation function.
		/// @param xv Vector of x-values (abscissas), must be sorted
		/// @param yv Vector of y-values (ordinates)
		/// @param m Number of points to use in interpolation (polynomial degree = m-1)

		PolynomInterpRealFunc(const Vector<Real>& xv, const Vector<Real>& yv, int m)
			: RealFunctionInterpolated(xv, yv, m)
			, _errorEst(0.) {}

		/// /** @brief Get the error estimate from the last interpolation. */

		Real getLastErrorEst() const { return _errorEst; }
		Real DetailedErrorEstimate() const override { return std::abs(_errorEst); }
		const char* InterpolationMethodName() const override { return "NevillePolynomial"; }

		/// @brief Calculate polynomial interpolated value using Neville's algorithm.
		/// @param startInd Starting index in the data array
		/// @param x Point at which to interpolate
		/// @return Interpolated value; error estimate stored in _errorEst
		/// Uses mm-point polynomial interpolation on the subrange [startInd, startInd+mm-1].

		Real calcInterpValue(int startInd, Real x) const override {
			int i, m, ns = 0;
			Real y, den, dif, dift, ho, hp, w, dy;
			Vector<Real> c(getInterpOrder()), d(getInterpOrder());

			dif = std::abs(x - X(startInd));

			// First we find the index ns of the closest table entry
			for (i = 0; i < getInterpOrder(); i++) {
				if ((dift = std::abs(x - X(startInd + i))) < dif) {
					ns = i;
					dif = dift;
				}
				c[i] = Y(startInd + i);
				d[i] = Y(startInd + i);
			}
			y = Y(startInd + ns--); // This is the initial approximation to y

			for (m = 1; m < getInterpOrder(); m++) // For each column of the tableau
			{
				for (i = 0; i < getInterpOrder() - m; i++) {
					// we loop over the current c's and d's and update
					ho = X(startInd + i) - x;
					hp = X(startInd + i + m) - x;
					w = c[i + 1] - d[i];

					if ((den = ho - hp) == 0.0)
						throw RealFuncInterpRuntimeError("Polynomial interpolation: identical x-values encountered");

					den = w / den;
					d[i] = hp * den;
					c[i] = ho * den;
				}
				dy = (2 * (ns + 1) < (getInterpOrder() - m) ? c[ns + 1] : d[ns--]);

				y += dy;
			}
			_errorEst = dy;
			return y;
		}

	};

} // namespace MML

#endif // MML_INTERPOLATED_REAL_FUNCTION_POLYNOMIAL_H