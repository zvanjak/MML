///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Polynom.h                                                           ///
///  Description: Polynomial class with arithmetic, evaluation, differentiation       ///
///               Root finding and algebraic operations                               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_POLYNOM_H
#define MML_POLYNOM_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/Vector.h"
#include "base/VectorN.h"
#include "base/MatrixNM.h"

namespace MML
{
	template <typename CoefT, typename FieldT = CoefT>
	class Polynom
	{
	protected:
		// polynom coefficients - with _vecCoef[0] being the constant term
		std::vector<CoefT> _vecCoef;
	public:
		Polynom() {}
		Polynom(int n) { _vecCoef.resize(n + 1); }
		Polynom(const std::vector<CoefT>& vecCoef) : _vecCoef(vecCoef) {}
		Polynom(std::initializer_list<CoefT> list) : _vecCoef(list) {}
		Polynom(const Polynom& Copy) = default;

		Polynom(Polynom&&) noexcept = default;
		Polynom& operator=(Polynom&&) noexcept = default;

		// Static constructors for common polynomials
		// Returns the zero polynomial (all coefficients zero)
		static Polynom Zero() {
			return Polynom(std::initializer_list<CoefT>{CoefT(0)});
		}

		// Returns the monomial x^degree (coefficient 1 at degree, 0 elsewhere)
		static Polynom Monomial(int degree) {
			std::vector<CoefT> coef(degree + 1, CoefT(0));
			if (degree >= 0) coef[degree] = CoefT(1);
			return Polynom(coef);
		}

		// Returns the constant polynomial (all coefficients zero except constant term)
		static Polynom Constant(const CoefT& value) {
			return Polynom(std::initializer_list<CoefT>{value});
		}

    // Given arrays x[0..n-1] and y[0..n-1] containing a tabulated function y[i] = f(x[i]), 
    // this routine returns a polynomial such that y[i] = P(x[i]) for all i.
    // Uses Lagrange interpolation via Newton's divided differences method.
    static Polynom FromValues(const std::vector<FieldT>& x, const std::vector<FieldT>& y)
    {
      int n = (int)x.size();
      if (n != (int)y.size())
        throw std::invalid_argument("FromValues: x and y arrays must have the same size");
      if (n == 0)
        throw std::invalid_argument("FromValues: arrays cannot be empty");

      std::vector<FieldT> cof(n, FieldT(0));
      std::vector<FieldT> s(n);
      
      // Initialize
      for (int i = 0; i < n; i++)
        s[i] = FieldT(0);
      
      s[n-1] = -x[0];
      
      // Build the polynomial basis
      for (int i = 1; i < n; i++) {
        for (int j = n-1-i; j < n-1; j++)
          s[j] -= x[i] * s[j+1];
        s[n-1] -= x[i];
      }
      
      // Compute coefficients using the basis
      for (int j = 0; j < n; j++) {
        FieldT phi = FieldT(n);
        for (int k = n-1; k > 0; k--)
          phi = FieldT(k) * s[k] + x[j] * phi;
      
        FieldT ff = y[j] / phi;
        FieldT b = FieldT(1);
        
        for (int k = n-1; k >= 0; k--) {
          cof[k] += b * ff;
          b = s[k] + x[j] * b;
        }
      }
      
      return Polynom(cof);
    }    

		// Returns the linear polynomial a*x + b
		static Polynom Linear(const CoefT& a, const CoefT& b) {
			return Polynom(std::initializer_list<CoefT>{b, a});
		}

		int  GetDegree() const { return (int)_vecCoef.size() - 1; }
		void SetDegree(int newDeg) { _vecCoef.resize(newDeg + 1); }
		bool IsNullPolynom() const { return _vecCoef.size() == 0; }
		void Reduce() { while (!_vecCoef.empty() && _vecCoef.back() == CoefT(0)) _vecCoef.pop_back(); }

		CoefT leadingTerm() const { return _vecCoef.empty() ? CoefT(0) : _vecCoef.back(); }
		CoefT constantTerm() const { return _vecCoef.empty() ? CoefT(0) : _vecCoef[0]; }

		const CoefT& operator[] (int i) const { return _vecCoef[i]; }
		CoefT& operator[] (int i) { return _vecCoef[i]; }

		///////////////////////////            Iterators              ///////////////////////////

		using iterator = typename std::vector<CoefT>::iterator;
		using const_iterator = typename std::vector<CoefT>::const_iterator;

		iterator begin() { return _vecCoef.begin(); }
		iterator end() { return _vecCoef.end(); }
		const_iterator begin() const { return _vecCoef.begin(); }
		const_iterator end() const { return _vecCoef.end(); }
		const_iterator cbegin() const { return _vecCoef.cbegin(); }
		const_iterator cend() const { return _vecCoef.cend(); }

		///////////////////////////            Operators              ///////////////////////////
		Polynom& operator=(const Polynom& Copy) { _vecCoef = Copy._vecCoef; return *this; }

		FieldT operator() (const FieldT& x) const {
			int j = GetDegree();
			FieldT p = _vecCoef[j] * FieldT(1);

			while (j > 0)
				p = p * x + _vecCoef[--j];
			return p;
		}

		bool IsEqual(const Polynom& b) const
		{
			// Compare polynomials by effective degree (ignoring trailing zeros)
			int maxSize = (int)std::max(_vecCoef.size(), b._vecCoef.size());
			
			for (int i = 0; i < maxSize; i++)
			{
				CoefT a_coef = (i < (int)_vecCoef.size()) ? _vecCoef[i] : CoefT(0);
				CoefT b_coef = (i < (int)b._vecCoef.size()) ? b._vecCoef[i] : CoefT(0);
				
				if (a_coef != b_coef)
					return false;
			}
			return true;
		}

		bool operator==(const Polynom& b) const
		{
			return IsEqual(b);
		}

		bool operator!=(const Polynom& b) const
		{
			return !IsEqual(b);
		}

		Polynom operator+(const Polynom& b) const
		{
			Polynom result;
			int n = (int)std::max(_vecCoef.size(), b._vecCoef.size());
			result._vecCoef.resize(n);
			for (int i = 0; i < n; i++)
			{
				if (i < _vecCoef.size())
					result._vecCoef[i] += _vecCoef[i];
				if (i < b._vecCoef.size())
					result._vecCoef[i] += b._vecCoef[i];
			}
			return result;
		}

		Polynom operator-(const Polynom& b) const
		{
			Polynom result;
			int n = (int)std::max(_vecCoef.size(), b._vecCoef.size());
			result._vecCoef.resize(n);
			for (int i = 0; i < n; i++)
			{
				if (i < _vecCoef.size())
					result._vecCoef[i] += _vecCoef[i];
				if (i < b._vecCoef.size())
					result._vecCoef[i] -= b._vecCoef[i];
			}
			return result;
		}

		Polynom operator*(const Polynom& b) const
		{
			Polynom result;

			int n = (int)(_vecCoef.size() + b._vecCoef.size() - 1);
			result._vecCoef.resize(n);
			for (int i = 0; i < _vecCoef.size(); i++)
				for (int j = 0; j < b._vecCoef.size(); j++)
					result._vecCoef[i + j] += _vecCoef[i] * b._vecCoef[j];
			return result;
		}

		// Addition by scalar (adds scalar to the constant term)
		Polynom operator+(const CoefT& scalar) const {
			Polynom result = *this;
			if (result._vecCoef.empty())
				result._vecCoef.resize(1, CoefT(0));
			result._vecCoef[0] += scalar;
			return result;
		}
		friend Polynom operator+(const CoefT& scalar, const Polynom& poly) {
			return poly + scalar;
		}

		// Subtraction by scalar (subtracts scalar from the constant term)
		Polynom operator-(const CoefT& scalar) const {
			Polynom result = *this;
			if (result._vecCoef.empty())
				result._vecCoef.resize(1, CoefT(0));
			result._vecCoef[0] -= scalar;
			return result;
		}
		friend Polynom operator-(const CoefT& scalar, const Polynom& poly) {
			Polynom result = -poly;
			if (result._vecCoef.empty())
				result._vecCoef.resize(1, CoefT(0));
			result._vecCoef[0] += scalar;
			return result;
		}

		// Multiplication by scalar
		Polynom operator*(const CoefT& scalar) const {
			Polynom result = *this;
			for (auto& coef : result._vecCoef)
				coef *= scalar;
			return result;
		}
		friend Polynom operator*(const CoefT& scalar, const Polynom& poly) {
			return poly * scalar;
		}

		// Division by scalar
		Polynom operator/(const CoefT& scalar) const {
			Polynom result = *this;
			for (auto& coef : result._vecCoef)
				coef /= scalar;
			return result;
		}

		static void poldiv(const Polynom& u, const Polynom& v, Polynom& qout, Polynom& rout)
		{
			int k, j, n = u.GetDegree(), nv = v.GetDegree();

			// find real degree of v
			while (nv >= 0 && v._vecCoef[nv] == CoefT(0))
				nv--;

			if (nv < 0)
				throw std::domain_error("poldiv divide by zero polynomial");

			Polynom r = u;
			Polynom q(u.GetDegree());
			for (k = n - nv; k >= 0; k--)
			{
				q[k] = r[nv + k] / v[nv];

				for (j = nv + k - 1; j >= k; j--)
					r[j] -= q[k] * v[j - k];
			}
			for (j = nv; j <= n; j++)
				r[j] = CoefT(0);


			int nq = q.GetDegree();
			while (nq >= 0 && q[nq] == CoefT(0))
				nq--;

			// setting exact size for quotient
			qout.SetDegree(nq);
			for (j = 0; j <= nq; j++)
				qout[j] = q[j];

			// setting exact size for remainder
			rout.SetDegree(nv - 1);
			for (j = 0; j < nv; j++)
				rout[j] = r[j];
		}

		///////////////////////////           Operations              ///////////////////////////
		// Given the coefficients of a polynomial of degree nc as an array c[0..nc] of size nc+1 (with
		// c[0] being the constant term), and given a value x, this routine fills an output array pd of size
		// nd+1 with the value of the polynomial evaluated at x in pd[0], and the first nd derivatives at
		// x in pd[1..nd].
		void Derive(const Real x, Vector<Real>& pd)
		{
			int  nnd, j, i;
			int  nc = GetDegree();
			int  nd = pd.size() - 1;
			Real cnst = 1.0;

			pd[0] = (*this)[nc];
			for (j = 1; j < nd + 1; j++)
				pd[j] = 0.0;

			for (i = nc - 1; i >= 0; i--)
			{
				nnd = (nd < (nc - i) ? nd : nc - i);
				for (j = nnd; j > 0; j--)
					pd[j] = pd[j] * x + pd[j - 1];

				pd[0] = pd[0] * x + (*this)[i];
			}
			for (i = 2; i < nd + 1; i++) {
				cnst *= i;
				pd[i] *= cnst;
			}
		}

		// Returns the first derivative of the polynomial as a new Polynom
		Polynom Derive() const {
			Polynom result;
			int deg = GetDegree();
			if (deg <= 0) {
				result._vecCoef.clear();
				return result;
			}
			result._vecCoef.resize(deg);
			for (int i = 1; i <= deg; ++i) {
				result._vecCoef[i - 1] = _vecCoef[i] * CoefT(i);
			}
			return result;
		}

		// Returns the indefinite integral of the polynomial as a new Polynom.
		// The constant of integration is set to zero.
		Polynom Integrate() const {
			Polynom result;
			int deg = GetDegree();
			result._vecCoef.resize(deg + 2, CoefT(0));
			for (int i = 0; i <= deg; ++i) {
				result._vecCoef[i + 1] = _vecCoef[i] / CoefT(i + 1);
			}
			return result;
		}

		///////////////////////////               I/O                 ///////////////////////////
		std::string to_string(int width, int precision) const
		{
			std::stringstream str;

			Print(str, width, precision);

			return str.str();
		}

		std::ostream& Print(std::ostream& stream) const
		{
			return Print(stream, 10, 6); // Default width and precision
		}

		std::ostream& Print(std::ostream& stream, int width, int precision) const
		{
			using std::abs;
			using std::to_string;

			std::ios_base::fmtflags f(stream.flags());
			stream << std::fixed << std::setw(width) << std::setprecision(precision);

			bool first = true;
			for (int i = (int)_vecCoef.size() - 1; i >= 0; --i)
			{
				// Skip zero coefficients
				bool is_zero = false;
				if constexpr (std::is_floating_point_v<CoefT>) {
					is_zero = (std::abs(_vecCoef[i]) < PrecisionValues<Real>::PolynomialCoeffZeroThreshold);
				}
				else if constexpr (std::is_same_v<CoefT, std::complex<double>> || std::is_same_v<CoefT, std::complex<float>>) {
					is_zero = (std::abs(_vecCoef[i]) < PrecisionValues<Real>::PolynomialCoeffZeroThreshold);
				}
				else {
					is_zero = (_vecCoef[i] == CoefT(0));
				}
				if (is_zero) continue;

				// Print sign and separator
				if (!first) {
					if constexpr (std::is_arithmetic_v<CoefT>) {
						stream << (_vecCoef[i] >= CoefT(0) ? " + " : " - ");
					}
					else {
						stream << " + "; // Always use '+' for non-ordered types
					}
				}
				else {
					if constexpr (std::is_arithmetic_v<CoefT>) {
						if (_vecCoef[i] < CoefT(0)) stream << "-";
					}
					first = false;
				}

				// Print coefficient (absolute value for arithmetic types, as-is for others)
				if constexpr (std::is_arithmetic_v<CoefT>) {
					CoefT abs_coef = _vecCoef[i] < CoefT(0) ? -_vecCoef[i] : _vecCoef[i];
					bool print_coef = (i == 0) || (abs_coef != CoefT(1));
					if (print_coef) stream << abs_coef;
				}
				else {
					stream << _vecCoef[i];
				}

				// Variable and exponent
				if (i > 0) {
					stream << "x";
					if (i > 1) stream << "^" << i;
				}
			}

			if (first) stream << CoefT(0);

			stream.flags(f);
			return stream;
		}

		friend std::ostream& operator<<(std::ostream& stream, const Polynom& poly)
		{
			return poly.Print(stream);
		}
	};

	class PolynomRealFunc : public Polynom<Real>, public IRealFunction
	{
	public:
		PolynomRealFunc() {}
		PolynomRealFunc(int n) { _vecCoef.resize(n + 1); }
		PolynomRealFunc(const std::vector<Real>& vecCoef) : Polynom(vecCoef) {}
		PolynomRealFunc(std::initializer_list<Real> list) : Polynom(list) {}
		PolynomRealFunc(const Polynom& Copy) : Polynom(Copy) {}
		~PolynomRealFunc() {}

		Real operator()(Real x) const { return Polynom::operator()(x); }
	};

	typedef Polynom<Real>      PolynomReal;
	typedef Polynom<Complex>   PolynomComplex;
	
	typedef Polynom<Real, Complex>   PolynomComplexRealCoef;

	typedef Polynom<Real, MatrixNM<Real, 2, 2>>       Matrix2Polynom;
	typedef Polynom<Real, MatrixNM<Real, 3, 3>>       Matrix3Polynom;
	typedef Polynom<Real, MatrixNM<Real, 4, 4>>       Matrix4Polynom;
}

#endif