///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        InterpolatedRealFunctionRational.h                                  ///
///  Description: Rational interpolation for real-valued functions                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file InterpolatedRealFunctionRational.h
/// @brief Rational function interpolation using the Bulirsch-Stoer algorithm.
/// @ingroup Interpolation

#if !defined MML_INTERPOLATED_REAL_FUNCTION_RATIONAL_H
#define MML_INTERPOLATED_REAL_FUNCTION_RATIONAL_H

#include <limits>

#include <mml/base/InterpolatedFunctions/InterpolatedRealFunctionLinear.h>

namespace MML {
	/////////////////////////////////////////////////////////////////////////////////////
	///                      RATIONAL INTERPOLATION                                   ///
	/////////////////////////////////////////////////////////////////////////////////////

	/// @brief Rational function interpolation using Bulirsch-Stoer algorithm.
	/// RationalInterpRealFunc constructs a diagonal rational function (equal degree
	/// numerator and denominator) through mm data points. This is superior to
	/// polynomial interpolation when the underlying function has poles or
	/// asymptotic behavior.
	/// **Algorithm:** Uses the Bulirsch-Stoer algorithm which builds up a continued
	/// fraction representation of the rational interpolant. Provides error estimation.
	/// **Advantages:**
	/// - Handles functions with poles or asymptotes
	/// - Better extrapolation behavior than polynomials
	/// - Built-in error estimation
	/// **Limitations:**
	/// - May detect poles during interpolation (throws exception)
	/// - More complex than polynomial methods
	/// @warning Will throw if a pole is detected during interpolation.
	/// @see PolynomInterpRealFunc for simpler polynomial alternative
	/// @see BarycentricRationalInterp for pole-free rational interpolation
	/// @ingroup Interpolation

	class RationalInterpRealFunc : public RealFunctionInterpolated {
	private:
		mutable Real _errorEst;
		static constexpr Real TINY = std::numeric_limits<Real>::min(); ///< Small positive value to prevent 0/0

	public:
		/// @brief Construct a rational interpolation function.
		/// @param xv Vector of x-values (abscissas), must be sorted
		/// @param yv Vector of y-values (ordinates)
		/// @param m Number of points to use (rational function has degree (m-1)/2 in each)

		RationalInterpRealFunc(const Vector<Real>& xv, const Vector<Real>& yv, int m)
			: RealFunctionInterpolated(xv, yv, m)
			, _errorEst(0.) {}

		/// /** @brief Get the error estimate from the last interpolation. */

		Real getLastErrorEst() const { return _errorEst; }
		Real DetailedErrorEstimate() const override { return std::abs(_errorEst); }
		const char* InterpolationMethodName() const override { return "BulirschStoerRational"; }

		/// @brief Calculate rational interpolated value using Bulirsch-Stoer algorithm.
		/// @param startInd Starting index in the data array
		/// @param x Point at which to interpolate
		/// @return Interpolated value; error estimate stored in _errorEst
		/// @throws RealFuncInterpRuntimeError if a pole is detected

		Real calcInterpValue(int startInd, Real x) const override {
			int m, i, ns = 0;
			Real y, w, t, hh, h, dd;
			int mm = getInterpOrder();
			Vector<Real> c(mm), d(mm);

			hh = std::abs(x - X(startInd));
			for (i = 0; i < mm; i++) {
				h = std::abs(x - X(startInd + i));
				const Real nodeScale = std::max({Real{1}, std::abs(x), std::abs(X(startInd + i))});
				const Real nodeTolerance = Real{4} * std::numeric_limits<Real>::epsilon() * nodeScale;
				if (h <= nodeTolerance) {
					_errorEst = 0.0;
					return Y(startInd + i);
				} else if (h < hh) {
					ns = i;
					hh = h;
				}
				c[i] = Y(startInd + i);
				d[i] = Y(startInd + i) + TINY; // TINY prevents rare zero-over-zero condition
			}
			y = Y(startInd + ns--);

			for (m = 1; m < mm; m++) {
				for (i = 0; i < mm - m; i++) {
					w = c[i + 1] - d[i];
					h = X(startInd + i + m) - x;
					if (std::abs(h) < TINY)
						throw RealFuncInterpRuntimeError("Error in RationalInterpRealFunc: x coincides with a node");
					t = (X(startInd + i) - x) * d[i] / h;
					dd = t - c[i + 1];
					if (dd == 0.0)
						throw RealFuncInterpRuntimeError("Error in RationalInterpRealFunc: pole detected");
					dd = w / dd;
					d[i] = c[i + 1] * dd;
					c[i] = t * dd;
				}
				_errorEst = (2 * (ns + 1) < (mm - m) ? c[ns + 1] : d[ns--]);
				y += _errorEst;
			}
			return y;
		}
	};

} // namespace MML

#endif // MML_INTERPOLATED_REAL_FUNCTION_RATIONAL_H