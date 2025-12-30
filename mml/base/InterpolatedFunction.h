///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        InterpolatedFunction.h                                              ///
///  Description: Interpolation classes (Linear, Polynomial, Spline, Rational)        ///
///               Table-driven function approximation and smoothing                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_INTERPOLATEDFUNCTION_H
#define MML_INTERPOLATEDFUNCTION_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/Vector.h"
#include "base/Matrix.h"

#include "base/Function.h"

namespace MML
{
	class RealFunctionInterpolated : public IRealFunction
	{
	private:
		int		_numPoints, _usedPoints;
		Vector<Real> _x, _y;						// we are storing copies of given values!

	public:
		RealFunctionInterpolated(const Vector<Real> &x, const Vector<Real> &y, 
														 int usedPointsInInterpolation)
							:  _x(x), _y(y), _numPoints(x.size()), _usedPoints(usedPointsInInterpolation)
		{
			// throw if not enough points
			if (_numPoints < 2 || _usedPoints < 2 || _usedPoints > _numPoints)
				throw RealFuncInterpInitError("RealFunctionInterpolated size error");
		}

		virtual ~RealFunctionInterpolated() {}

		Real virtual calcInterpValue(int startInd, Real x) const = 0;

		inline Real MinX() const { return X(0); }
		inline Real MaxX() const { return X(_numPoints-1); }

		inline Real X(int i) const { return _x[i]; }
		inline Real Y(int i) const { return _y[i]; }

		inline int	getNumPoints() const { return _numPoints; }
		inline int  getInterpOrder() const { return _usedPoints; }

		Real operator()(Real x) const
		{
			int startInd = locate(x);
			return calcInterpValue(startInd, x);
		}

		// Given a value x, return a value j such that x is (insofar as possible) centered in the subrange
		// xx[j..j+mm-1], where xx is the stored pointer. The values in xx must be monotonic, either
		// increasing or decreasing. The returned value is not less than 0, nor greater than _numPoints-1.
		int locate(const Real x) const
		{
			int  indUpper, indMid, indLower;
			bool isAscending = (_x[_numPoints - 1] >= _x[0]);

			indLower = 0;
			indUpper = _numPoints - 1;

			while (indUpper - indLower > 1) {
				indMid = (indUpper + indLower) / 2;
				if ((x >= _x[indMid]) == isAscending)
					indLower = indMid;
				else
					indUpper = indMid;
			}
			return std::max(0, std::min(_numPoints - _usedPoints, indLower - ((_usedPoints - 2) / 2)));
		}
	};

	////////////////////////             LINEAR INTERPOLATION						  ///////////////////////
	class LinearInterpRealFunc : public RealFunctionInterpolated
	{
		bool _extrapolateOutsideOfRange;
	public:
		LinearInterpRealFunc(const Vector<Real>& xv, const Vector<Real>& yv, bool extrapolateOutsideOfRange = false) 
						: RealFunctionInterpolated(xv, yv, 2), _extrapolateOutsideOfRange(extrapolateOutsideOfRange) {}

		Real calcInterpValue(int j, Real x) const override {
			if(_extrapolateOutsideOfRange == false && (x < MinX() || x > MaxX()) )
				return 0.0;

			if (X(j) == X(j + 1)) 
				return Y(j);
			else 
				return Y(j) + ((x - X(j)) / (X(j + 1) - X(j)) * (Y(j + 1) - Y(j)));
		}
	};

	////////////////////////           POLYNOMIAL INTERPOLATION						///////////////////////
	// Polynomial interpolation object.Construct with x and y vectors, and the number M of points
	// to be used locally(polynomial order plus one), then call interp for interpolated values.
	class PolynomInterpRealFunc : public RealFunctionInterpolated
	{
	private:
		mutable Real _errorEst;
	public:
		PolynomInterpRealFunc(const Vector<Real>& xv, const Vector<Real>& yv, int m)
						: RealFunctionInterpolated(xv, yv, m), _errorEst(0.) { }

		Real getLastErrorEst() const { return _errorEst; }	

		// Given a value x, and using pointers to data xx and yy, this routine returns an interpolated
		// value y, and stores an error estimate _errorEst.The returned value is obtained by mm-point polynomial
		// interpolation on the subrange xx[startInd..startInd + mm - 1].
		Real calcInterpValue(int startInd, Real x) const override
		{
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
			y = Y(startInd + ns--);				// This is the initial approximation to y
			
			for (m = 1; m < getInterpOrder(); m++)  	// For each column of the tableau
			{
				for (i = 0; i < getInterpOrder() - m; i++) 
				{
					// we loop over the current c's and d's and update
					ho = X(startInd + i) - x;
					hp = X(startInd + i + m) - x;
					w = c[i + 1] - d[i];
					
					if ((den = ho - hp) == 0.0) 
						throw("Poly_interp error");
					
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

	///////////////////////           CUBIC SPLINE INTERPOLATION						//////////////////////
	// Construct with x and y vectors, and (optionally) values of
	// the first derivative at the endpoints, then call interp for interpolated values.
	//
	// Natural spline: yp1 = ypn = 1e99 (default) - zero second derivative at endpoints
	// Clamped spline: specify yp1 and/or ypn as first derivatives at endpoints
	class SplineInterpRealFunc : public RealFunctionInterpolated
	{
	private:
		Vector<Real> _secDerY;

	public:

		SplineInterpRealFunc(const Vector<Real>& xv, const Vector<Real>& yv, Real yp1 = 1.e99, Real ypn = 1.e99)
			: RealFunctionInterpolated(xv, yv, 2), _secDerY(xv.size())
		{
			initSecDerivs(&xv[0], &yv[0], yp1, ypn);
		}

		// This routine stores an array _secDerY[0..numPoints-1] with second derivatives of the interpolating function
		// at the tabulated points pointed to by xv, using function values pointed to by yv. If yp1 and/or
		// ypn are equal to 10e99 or larger, the routine is signaled to set the corresponding boundary
		// condition for a natural spline, with zero second derivative on that boundary; otherwise, they are
		// the values of the first derivatives at the endpoints.
		void initSecDerivs(const Real* xv, const Real* yv, Real yp1, Real ypn)
		{
			int		i, k;
			Real	p, qn, sig, un;

			int numPoints = (int)_secDerY.size();
			Vector<Real> u(numPoints - 1);

			if (yp1 > 0.99e99)
				_secDerY[0] = u[0] = 0.0;
			else 
			{
				_secDerY[0] = -0.5;
				u[0] = (3.0 / (xv[1] - xv[0])) * ((yv[1] - yv[0]) / (xv[1] - xv[0]) - yp1);
			}
			for (i = 1; i < numPoints - 1; i++) 
			{
				sig = (xv[i] - xv[i - 1]) / (xv[i + 1] - xv[i - 1]);
				p = sig * _secDerY[i - 1] + 2.0;

				_secDerY[i] = (sig - 1.0) / p;

				u[i] = (yv[i + 1] - yv[i]) / (xv[i + 1] - xv[i]) - (yv[i] - yv[i - 1]) / (xv[i] - xv[i - 1]);
				u[i] = (6.0 * u[i] / (xv[i + 1] - xv[i - 1]) - sig * u[i - 1]) / p;
			}
			if (ypn > 0.99e99)
				qn = un = 0.0;
			else 
			{
				qn = 0.5;
				un = (3.0 / (xv[numPoints - 1] - xv[numPoints - 2])) * (ypn - (yv[numPoints - 1] - yv[numPoints - 2]) / (xv[numPoints - 1] - xv[numPoints - 2]));
			}

			_secDerY[numPoints - 1] = (un - qn * u[numPoints - 2]) / (qn * _secDerY[numPoints - 2] + 1.0);

			for (k = numPoints - 2; k >= 0; k--)
				_secDerY[k] = _secDerY[k] * _secDerY[k + 1] + u[k];
		}

		// Given a value x, and using pointers to data xx and yy, and the stored vector of second derivatives
		// _secDerY, this routine returns the cubic spline interpolated value y.        
		Real calcInterpValue(int startInd, Real x) const override
		{
			int indLow = startInd, indUpp = startInd + 1;
			Real y, h, b, a;
			h = X(indUpp) - X(indLow);
			
			if (h == 0.0) 
				throw("Bad input to routine splint");
			
			a = (X(indUpp) - x) / h;
			b = (x - X(indLow)) / h;

			y = a * Y(indLow) + b * Y(indUpp) + 
				  ((POW3(a) - a) * _secDerY[indLow] + (POW3(b) - b) * _secDerY[indUpp]) * (h * h) / 6.0;
			
			return y;
		}

		// Evaluate the first derivative of the cubic spline at point x
		// Returns dy/dx at the given point
		Real Derivative(Real x) const
		{
			int startInd = locate(x);
			int indLow = startInd, indUpp = startInd + 1;
			Real h, b, a;
			h = X(indUpp) - X(indLow);
			
			if (h == 0.0) 
				throw("Bad input to routine splint derivative");
			
			a = (X(indUpp) - x) / h;
			b = (x - X(indLow)) / h;

			// Derivative of cubic spline:
			// dy/dx = (y[hi] - y[lo])/h - (3a² - 1)/6 * h * y2[lo] + (3b² - 1)/6 * h * y2[hi]
			Real dydx = (Y(indUpp) - Y(indLow)) / h
				- (3.0 * a * a - 1.0) / 6.0 * h * _secDerY[indLow]
				+ (3.0 * b * b - 1.0) / 6.0 * h * _secDerY[indUpp];
			
			return dydx;
		}

		// Evaluate the second derivative of the cubic spline at point x
		// Returns d²y/dx² at the given point (linearly interpolated between knots)
		Real SecondDerivative(Real x) const
		{
			int startInd = locate(x);
			int indLow = startInd, indUpp = startInd + 1;
			Real h, b, a;
			h = X(indUpp) - X(indLow);
			
			if (h == 0.0) 
				throw("Bad input to routine splint second derivative");
			
			a = (X(indUpp) - x) / h;
			b = (x - X(indLow)) / h;

			// Second derivative is linear interpolation of stored second derivatives
			return a * _secDerY[indLow] + b * _secDerY[indUpp];
		}

		// Compute definite integral of spline from a to b
		// Uses exact integration of cubic polynomials in each segment
		Real Integrate(Real a, Real b) const
		{
			if (a > b) return -Integrate(b, a);
			if (a < MinX() || b > MaxX())
				throw("Integration bounds outside spline domain");

			Real integral = 0.0;
			int iLow = locate(a);
			int iHigh = locate(b);

			// Helper lambda to integrate one segment from x0 to x1 within [X(i), X(i+1)]
			auto integrateSegment = [this](int i, Real x0, Real x1) -> Real {
				Real h = X(i + 1) - X(i);
				if (h == 0.0) return 0.0;

				// Compute a and b at both endpoints
				Real a0 = (X(i + 1) - x0) / h;
				Real b0 = (x0 - X(i)) / h;
				Real a1 = (X(i + 1) - x1) / h;
				Real b1 = (x1 - X(i)) / h;

				// Integral of spline segment: ∫[a*y_lo + b*y_hi + ((a³-a)*y2_lo + (b³-b)*y2_hi)*h²/6] dx
				// Using substitution and exact integration
				Real y_lo = Y(i), y_hi = Y(i + 1);
				Real y2_lo = _secDerY[i], y2_hi = _secDerY[i + 1];

				// Antiderivative evaluated at x1 minus at x0
				// For term a*y_lo: integral is -h/2 * a² * y_lo
				// For term b*y_hi: integral is h/2 * b² * y_hi
				// For (a³-a)*h²/6*y2_lo: integral is -h/6 * (a⁴/4 - a²/2) * h * y2_lo = -h³/6 * (a⁴/4 - a²/2) * y2_lo
				// For (b³-b)*h²/6*y2_hi: integral is h/6 * (b⁴/4 - b²/2) * h * y2_hi = h³/6 * (b⁴/4 - b²/2) * y2_hi

				Real term1 = h * y_lo * 0.5 * (a0*a0 - a1*a1);
				Real term2 = h * y_hi * 0.5 * (b1*b1 - b0*b0);
				Real term3 = h*h*h / 6.0 * y2_lo * ((a0*a0*a0*a0/4.0 - a0*a0/2.0) - (a1*a1*a1*a1/4.0 - a1*a1/2.0));
				Real term4 = h*h*h / 6.0 * y2_hi * ((b1*b1*b1*b1/4.0 - b1*b1/2.0) - (b0*b0*b0*b0/4.0 - b0*b0/2.0));

				return term1 + term2 + term3 + term4;
			};

			if (iLow == iHigh) {
				// Both endpoints in same segment
				integral = integrateSegment(iLow, a, b);
			}
			else {
				// Multiple segments
				integral += integrateSegment(iLow, a, X(iLow + 1));
				for (int i = iLow + 1; i < iHigh; i++) {
					integral += integrateSegment(i, X(i), X(i + 1));
				}
				integral += integrateSegment(iHigh, X(iHigh), b);
			}

			return integral;
		}

		// Access second derivative at node i
		Real GetSecondDerivative(int i) const {
			if (i < 0 || i >= getNumPoints())
				throw std::out_of_range("Index out of range in GetSecondDerivative");
			return _secDerY[i];
		}
	};


	///////////////////////           RATIONAL INTERPOLATION                      //////////////////////
	// Rational function interpolation using diagonal rational function through mm points.
	// Better than polynomial interpolation when data has poles or asymptotic behavior.
	// Uses Bulirsch-Stoer algorithm (continued fraction representation).
	class RationalInterpRealFunc : public RealFunctionInterpolated
	{
	private:
		mutable Real _errorEst;
		static constexpr Real TINY = 1.0e-99;  // Small number to prevent division by zero
	
	public:
		RationalInterpRealFunc(const Vector<Real>& xv, const Vector<Real>& yv, int m)
			: RealFunctionInterpolated(xv, yv, m), _errorEst(0.) { }

		Real getLastErrorEst() const { return _errorEst; }

		// Rational function interpolation using Bulirsch-Stoer algorithm.
		// Returns interpolated value and stores error estimate in _errorEst.
		Real calcInterpValue(int startInd, Real x) const override
		{
			int m, i, ns = 0;
			Real y, w, t, hh, h, dd;
			int mm = getInterpOrder();
			Vector<Real> c(mm), d(mm);

			hh = std::abs(x - X(startInd));
			for (i = 0; i < mm; i++) {
				h = std::abs(x - X(startInd + i));
				if (h == 0.0) {
					_errorEst = 0.0;
					return Y(startInd + i);  // Exact hit on a data point
				}
				else if (h < hh) {
					ns = i;
					hh = h;
				}
				c[i] = Y(startInd + i);
				d[i] = Y(startInd + i) + TINY;  // TINY prevents rare zero-over-zero condition
			}
			y = Y(startInd + ns--);

			for (m = 1; m < mm; m++) {
				for (i = 0; i < mm - m; i++) {
					w = c[i + 1] - d[i];
					h = X(startInd + i + m) - x;
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


	///////////////////////        BARYCENTRIC RATIONAL INTERPOLATION             //////////////////////
	// Barycentric rational interpolation with no poles on the real line.
	// Floater-Hormann algorithm - provides smooth interpolation without the
	// oscillations of high-degree polynomial interpolation (Runge phenomenon).
	// Parameter d controls the "degree" - higher d gives smoother interpolation
	// but may be less accurate for rapidly varying functions.
	class BarycentricRationalInterp : public IRealFunction
	{
	private:
		int _n;             // Number of data points
		int _d;             // Blending parameter (0 to n-1)
		Vector<Real> _x;    // Abscissas
		Vector<Real> _y;    // Ordinates
		Vector<Real> _w;    // Barycentric weights

		void computeWeights() {
			for (int k = 0; k < _n; k++) {
				int imin = std::max(k - _d, 0);
				int imax = k >= _n - _d ? _n - _d - 1 : k;
				Real temp = (imin & 1) ? -1.0 : 1.0;  // Sign alternates
				Real sum = 0.0;
				
				for (int i = imin; i <= imax; i++) {
					int jmax = std::min(i + _d, _n - 1);
					Real term = 1.0;
					for (int j = i; j <= jmax; j++) {
						if (j == k) continue;
						term *= (_x[k] - _x[j]);
					}
					term = temp / term;
					temp = -temp;
					sum += term;
				}
				_w[k] = sum;
			}
		}

	public:
		// Construct with data points and blending parameter d.
		// d = 0 gives piecewise constant, d = 1 gives piecewise linear,
		// d = n-1 gives polynomial interpolation (equivalent to Lagrange).
		// Recommended: d = 3 to 5 for most applications.
		BarycentricRationalInterp(const Vector<Real>& xv, const Vector<Real>& yv, int d = 3)
			: _n(xv.size()), _d(d), _x(xv), _y(yv), _w(_n)
		{
			if (_n <= _d)
				throw RealFuncInterpInitError("BarycentricRationalInterp: d too large for number of points");
			if (_d < 0)
				throw RealFuncInterpInitError("BarycentricRationalInterp: d must be non-negative");
			computeWeights();
		}

		Real MinX() const { return _x[0]; }
		Real MaxX() const { return _x[_n - 1]; }
		int  getNumPoints() const { return _n; }
		int  getBlendingParameter() const { return _d; }

		Real operator()(Real x) const override
		{
			Real num = 0.0, den = 0.0;
			for (int i = 0; i < _n; i++) {
				Real h = x - _x[i];
				if (h == 0.0) {
					return _y[i];  // Exact hit on a data point
				}
				Real temp = _w[i] / h;
				num += temp * _y[i];
				den += temp;
			}
			return num / den;
		}

		// Get barycentric weight at index i
		Real getWeight(int i) const {
			if (i < 0 || i >= _n)
				throw std::out_of_range("Index out of range in getWeight");
			return _w[i];
		}
	};


	///////////////////////            2D INTERPOLATION CLASSES                   //////////////////////

	///////////////////////             BILINEAR INTERPOLATION                    //////////////////////
	// Bilinear interpolation on a regular 2D grid.
	// Given data z[i][j] at grid points (x1[i], x2[j]), interpolates to any point (x1, x2).
	class BilinearInterp2D
	{
	private:
		int _m, _n;                    // Grid dimensions
		Vector<Real> _x1, _x2;         // Grid coordinates
		const Matrix<Real>& _z;        // Reference to data matrix z[m][n]

		// Binary search to find interval containing x
		int locate(const Vector<Real>& arr, Real x) const {
			int n = arr.size();
			int lo = 0, hi = n - 1;
			bool ascending = (arr[n - 1] >= arr[0]);

			while (hi - lo > 1) {
				int mid = (hi + lo) / 2;
				if ((x >= arr[mid]) == ascending)
					lo = mid;
				else
					hi = mid;
			}
			return std::max(0, std::min(n - 2, lo));
		}

	public:
		BilinearInterp2D(const Vector<Real>& x1v, const Vector<Real>& x2v, const Matrix<Real>& zm)
			: _m(x1v.size()), _n(x2v.size()), _x1(x1v), _x2(x2v), _z(zm)
		{
			if (zm.RowNum() != _m || zm.ColNum() != _n)
				throw RealFuncInterpInitError("BilinearInterp2D: matrix dimensions don't match grid");
		}

		Real operator()(Real x1, Real x2) const
		{
			int i = locate(_x1, x1);
			int j = locate(_x2, x2);

			// Compute normalized coordinates t, u in [0, 1]
			Real t = (_x1[i + 1] - _x1[i] != 0.0) ? 
					 (x1 - _x1[i]) / (_x1[i + 1] - _x1[i]) : 0.0;
			Real u = (_x2[j + 1] - _x2[j] != 0.0) ? 
					 (x2 - _x2[j]) / (_x2[j + 1] - _x2[j]) : 0.0;

			// Bilinear interpolation formula
			return (1.0 - t) * (1.0 - u) * _z(i, j) 
				 + t * (1.0 - u) * _z(i + 1, j)
				 + (1.0 - t) * u * _z(i, j + 1) 
				 + t * u * _z(i + 1, j + 1);
		}

		// Evaluate with partial derivatives
		void interpWithDerivatives(Real x1, Real x2, Real& z, Real& dz_dx1, Real& dz_dx2) const
		{
			int i = locate(_x1, x1);
			int j = locate(_x2, x2);

			Real dx1 = _x1[i + 1] - _x1[i];
			Real dx2 = _x2[j + 1] - _x2[j];

			Real t = (dx1 != 0.0) ? (x1 - _x1[i]) / dx1 : 0.0;
			Real u = (dx2 != 0.0) ? (x2 - _x2[j]) / dx2 : 0.0;

			// Function value
			z = (1.0 - t) * (1.0 - u) * _z(i, j) 
			  + t * (1.0 - u) * _z(i + 1, j)
			  + (1.0 - t) * u * _z(i, j + 1) 
			  + t * u * _z(i + 1, j + 1);

			// Partial derivatives
			if (dx1 != 0.0) {
				dz_dx1 = ((1.0 - u) * (_z(i + 1, j) - _z(i, j)) 
						+ u * (_z(i + 1, j + 1) - _z(i, j + 1))) / dx1;
			} else {
				dz_dx1 = 0.0;
			}

			if (dx2 != 0.0) {
				dz_dx2 = ((1.0 - t) * (_z(i, j + 1) - _z(i, j)) 
						+ t * (_z(i + 1, j + 1) - _z(i + 1, j))) / dx2;
			} else {
				dz_dx2 = 0.0;
			}
		}

		int getGridDim1() const { return _m; }
		int getGridDim2() const { return _n; }
	};


	///////////////////////           BICUBIC SPLINE INTERPOLATION                //////////////////////
	// 2D spline interpolation on a regular grid using cubic splines in both directions.
	// More accurate than bilinear for smooth data, with continuous first derivatives.
	class BicubicSplineInterp2D
	{
	private:
		int _m, _n;
		Vector<Real> _x1, _x2;
		Matrix<Real> _z;               // Local copy of data
		std::vector<SplineInterpRealFunc*> _rowSplines;  // Splines along x2 for each x1[i]

	public:
		BicubicSplineInterp2D(const Vector<Real>& x1v, const Vector<Real>& x2v, const Matrix<Real>& zm)
			: _m(x1v.size()), _n(x2v.size()), _x1(x1v), _x2(x2v), _z(zm), _rowSplines(_m)
		{
			if (zm.RowNum() != _m || zm.ColNum() != _n)
				throw RealFuncInterpInitError("BicubicSplineInterp2D: matrix dimensions don't match grid");

			// Create a spline along x2 for each row (fixed x1[i])
			for (int i = 0; i < _m; i++) {
				Vector<Real> row(_n);
				for (int j = 0; j < _n; j++)
					row[j] = _z(i, j);
				_rowSplines[i] = new SplineInterpRealFunc(_x2, row);
			}
		}

		~BicubicSplineInterp2D() {
			for (int i = 0; i < _m; i++)
				delete _rowSplines[i];
		}

		// Disable copy to avoid double-free
		BicubicSplineInterp2D(const BicubicSplineInterp2D&) = delete;
		BicubicSplineInterp2D& operator=(const BicubicSplineInterp2D&) = delete;

		Real operator()(Real x1, Real x2) const
		{
			// First interpolate along x2 for each x1[i] to get values at fixed x2
			Vector<Real> yv(_m);
			for (int i = 0; i < _m; i++)
				yv[i] = (*_rowSplines[i])(x2);

			// Then interpolate along x1
			SplineInterpRealFunc colSpline(_x1, yv);
			return colSpline(x1);
		}

		// Evaluate with partial derivatives
		void interpWithDerivatives(Real x1, Real x2, Real& z, Real& dz_dx1, Real& dz_dx2) const
		{
			// Get values and x2-derivatives along rows
			Vector<Real> yv(_m), dyv_dx2(_m);
			for (int i = 0; i < _m; i++) {
				yv[i] = (*_rowSplines[i])(x2);
				dyv_dx2[i] = _rowSplines[i]->Derivative(x2);
			}

			// Interpolate values along x1
			SplineInterpRealFunc colSpline(_x1, yv);
			z = colSpline(x1);
			dz_dx1 = colSpline.Derivative(x1);

			// Interpolate x2-derivatives along x1
			SplineInterpRealFunc dSpline(_x1, dyv_dx2);
			dz_dx2 = dSpline(x1);
		}

		int getGridDim1() const { return _m; }
		int getGridDim2() const { return _n; }
	};


	//////////////////          PARAMETRIC CURVE LINEAR INTERPOLATION						/////////////////
	// Object for interpolating a curve specified by _numPoints points in N dimensions.
	template<int N>
	class LinInterpParametricCurve : public IParametricCurve<N>
	{
		int _numPoints;
		Matrix<Real> _curvePoints;
		Vector<Real> _arcLen; // normalized arc length parameter [0,1] for each point
		bool _isCurveClosed;

	public:
		LinInterpParametricCurve(const Matrix<Real>& ptsin, bool close = false)
			: _numPoints(ptsin.RowNum()), _curvePoints(ptsin), _arcLen(_numPoints), _isCurveClosed(close)
		{
			// Compute normalized arc length parameterization
			_arcLen[0] = 0.0;
			Real totalLen = 0.0;
			for (int i = 1; i < _numPoints; ++i) {
				Real segLen = 0.0;
				for (int j = 0; j < N; ++j)
					segLen += std::pow(_curvePoints[i][j] - _curvePoints[i - 1][j], 2);
				totalLen += std::sqrt(segLen);
				_arcLen[i] = totalLen;
			}
			// Normalize to [0,1]
			for (int i = 0; i < _numPoints; ++i)
				_arcLen[i] /= totalLen > 0 ? totalLen : 1.0;
		}

		Real getMinT() const override { return 0.0; }
		Real getMaxT() const override { return 1.0; }

		VectorN<Real, N> operator()(Real t) const override
		{
			// Clamp or wrap t
			if (_isCurveClosed) {
				t = t - std::floor(t);
			}
			else {
				if (t <= 0.0) t = 0.0;
				if (t >= 1.0) t = 1.0;
			}

			// Find segment
			int i = 0;
			while (i < _numPoints - 2 && t > _arcLen[i + 1])
				++i;

			// Linear interpolation
			Real t0 = _arcLen[i], t1 = _arcLen[i + 1];
			Real alpha = (t1 - t0 > 0) ? (t - t0) / (t1 - t0) : 0.0;

			VectorN<Real, N> result;
			for (int j = 0; j < N; ++j)
				result[j] = (1 - alpha) * _curvePoints[i][j] + alpha * _curvePoints[i + 1][j];
			return result;
		}
	};

	//////////////////          PARAMETRIC CURVE SPLINE INTERPOLATION						/////////////////
	// Object for interpolating a curve specified by _numPoints points in N dimensions.
	template<int N>
	class SplineInterpParametricCurve : public IParametricCurve<N>
	{
		Real _minT, _maxT;
		int  _dim, _numPoints, _bemba;
		bool _isCurveClosed;

		Matrix<Real> _curvePoints;
		Vector<Real> s;
		Vector<Real> ans;

		std::vector<SplineInterpRealFunc*> srp;
	
	public:
		// Constructor. The _numPoints _dim matrix ptsin inputs the data points. Input close as 0 for
		// an open curve, 1 for a closed curve. (For a closed curve, the last data point should not
		// duplicate the first - the algorithm will connect them.)
		SplineInterpParametricCurve(Real minT, Real maxT, const Matrix<Real>& ptsin, bool close = 0)
			: _numPoints(ptsin.RowNum()), _dim(ptsin.ColNum()), _bemba(close ? 2 * _numPoints : _numPoints),
			_isCurveClosed(close), _curvePoints(_dim, _bemba), s(_bemba), ans(_dim), srp(_dim),
			_minT(minT), _maxT(maxT)
		{
			// check N == dim
			if (N != _dim)
				throw("SplineInterpParametricCurve: N != dim");

			int i, ii, im, j, ofs;
			Real ss, soff, db, de;

			ofs = close ? _numPoints / 2 : 0;
			s[0] = 0.;
			for (i = 0; i < _bemba; i++) 
			{
				ii = (i - ofs + _numPoints) % _numPoints;
				im = (ii - 1 + _numPoints) % _numPoints;

				for (j = 0; j < _dim; j++) 
					_curvePoints[j][i] = ptsin[ii][j];
				
				if (i > 0) 
				{
					s[i] = s[i - 1] + rad(&ptsin[ii][0], &ptsin[im][0]);
					
					if (s[i] == s[i - 1]) 
						throw("error in Curve_interp");
					// Consecutive points may not be identical. For a closed curve, the last data
					// point should not duplicate the first.                    
				}
			}
			ss = close ? s[ofs + _numPoints] - s[ofs] : s[_numPoints - 1] - s[0];
			soff = s[ofs];
			
			for (i = 0; i < _bemba; i++) 
				s[i] = (s[i] - soff) / ss;
			
			for (j = 0; j < _dim; j++) 
			{
				db = _bemba < 4 ? 1.e99 : fprime(&s[0], &_curvePoints[j][0], 1);
				de = _bemba < 4 ? 1.e99 : fprime(&s[_bemba - 1], &_curvePoints[j][_bemba - 1], -1);

				Vector<Real> vec = _curvePoints.VectorFromRow(j);

				srp[j] = new SplineInterpRealFunc(s, vec, db, de);
			}
		}

		SplineInterpParametricCurve(const Matrix<Real>& ptsin, bool close = 0)
			: SplineInterpParametricCurve(0.0, 1.0, ptsin, close) {	}

		~SplineInterpParametricCurve() {
			for (int j = 0; j < _dim; j++) 
				delete srp[j];
		}
		
		Real getMinT() const { return _minT; }
		Real getMaxT() const { return _maxT; }

		// Interpolate a point on the stored curve. The point is parameterized by t, in the range [0,1].
		// For open curves, values of t outside this range will return extrapolations (dangerous!). For
		// closed curves, t is periodic with period 1
		VectorN<Real, N> operator()(Real t) const
		{
			if (_isCurveClosed == true && (t < _minT || t > _maxT))
				throw("SplineInterpParametricCurve: t outside interval");

			VectorN<Real, N> ans;

			if (_isCurveClosed)
				t = t - floor(t);

			// we have to map t from [minT, maxT] to [0, 1]
			t = (t - _minT) / (_maxT - _minT);
			for (int j = 0; j < _dim; j++)
				ans[j] = (*srp[j])(t);

			return ans;
		}

		// Utility for estimating the derivatives at the endpoints. x and y point to the abscissa and
		// ordinate of the endpoint. If pm is C1, points to the right will be used (left endpoint); if it
		// is -1, points to the left will be used (right endpoint). 
		Real fprime(Real* x, Real* y, int pm) {
			Real s1 = x[0] - x[pm * 1], s2 = x[0] - x[pm * 2], s3 = x[0] - x[pm * 3],
				s12 = s1 - s2, s13 = s1 - s3, s23 = s2 - s3;
			return -(s1 * s2 / (s13 * s23 * s3)) * y[pm * 3] + (s1 * s3 / (s12 * s2 * s23)) * y[pm * 2]
				- (s2 * s3 / (s1 * s12 * s13)) * y[pm * 1] + (1. / s1 + 1. / s2 + 1. / s3) * y[0];
		}

		Real rad(const Real* p1, const Real* p2) {
			Real sum = 0.;
			for (int i = 0; i < _dim; i++)
				sum += POW2(p1[i] - p2[i]);
			return sqrt(sum);
		}
	};


	template<int N>
	class InterpolatedSurface : public IParametricSurfaceRect<N>
	{
	public:
		InterpolatedSurface() {}

		VectorN<Real, N> operator()(const VectorN<Real, 2>& x) const { return VectorN<Real, N>{}; }
	};
}

#endif