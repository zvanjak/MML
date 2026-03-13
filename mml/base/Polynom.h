///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Polynom.h                                                           ///
///  Description: Polynomial class with arithmetic, evaluation, differentiation       ///
///               Root finding and algebraic operations                               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_POLYNOM_H
#define MML_POLYNOM_H

#include "MMLBase.h"
#include "MMLExceptions.h"

#include "interfaces/IFunction.h"

#include "base/Vector/Vector.h"
#include "base/Vector/VectorN.h"
#include "base/Matrix/MatrixNM.h"

// Standard headers - include what we use
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace MML
{
	/// @brief Polynomial class with arithmetic, evaluation, differentiation, and root finding.
	/// @tparam CoefT Coefficient type (Real, Complex, etc.)
	/// @tparam FieldT Field type for evaluation (defaults to CoefT)
	template <typename CoefT, typename FieldT = CoefT>
	class Polynom
	{
	protected:
		/// @brief Polynomial coefficients with _vecCoef[0] being the constant term
		std::vector<CoefT> _vecCoef;
	public:
		/// @brief Default constructor.
		Polynom() {}
		/// @brief Constructs polynomial of degree n with zero coefficients.
		/// @param n Degree of polynomial
		Polynom(int n) { _vecCoef.resize(n + 1); }
		/// @brief Constructs polynomial from coefficient vector.
		/// @param vecCoef Coefficients (vecCoef[0] = constant term)
		Polynom(const std::vector<CoefT>& vecCoef) : _vecCoef(vecCoef) {}
		/// @brief Constructs polynomial from initializer list.
		/// @param list Coefficients in ascending order
		Polynom(std::initializer_list<CoefT> list) : _vecCoef(list) {}
		/// @brief Copy constructor.
		Polynom(const Polynom& Copy) = default;

		/// @brief Constructs a constant polynomial from a single coefficient value.
		/// @param value Constant term value
		/// @note Uses explicit to avoid ambiguity with Polynom(int)
		explicit Polynom(const CoefT& value) : _vecCoef{value} {}

		/// @brief Move constructor.
		Polynom(Polynom&&) noexcept = default;
		/// @brief Move assignment.
		Polynom& operator=(Polynom&&) noexcept = default;

		/// @brief Returns the zero polynomial (all coefficients zero).
		static Polynom Zero() {
			return Polynom(std::initializer_list<CoefT>{CoefT(0)});
		}

		/// @brief Returns the monomial x^degree (coefficient 1 at degree, 0 elsewhere).
		/// @param degree Degree of monomial
		static Polynom Monomial(int degree) {
			std::vector<CoefT> coef(degree + 1, CoefT(0));
			if (degree >= 0) coef[degree] = CoefT(1);
			return Polynom(coef);
		}

		/// @brief Returns the constant polynomial.
		/// @param value Constant value
		static Polynom Constant(const CoefT& value) {
			return Polynom(std::initializer_list<CoefT>{value});
		}

		/// @brief Constructs interpolating polynomial from data points using Lagrange/Newton method.
		/// @param x X-coordinates of data points
		/// @param y Y-coordinates of data points
		/// @return Polynomial such that P(x[i]) = y[i]
    static Polynom FromValues(const std::vector<FieldT>& x, const std::vector<FieldT>& y)
    {
      int n = (int)x.size();
      if (n != (int)y.size())
        throw ArgumentError("FromValues: x and y arrays must have the same size");
      if (n == 0)
        throw ArgumentError("FromValues: arrays cannot be empty");

      // Check for duplicate x-values (causes division by zero in interpolation)
      for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
          if (x[i] == x[j])
            throw ArgumentError("FromValues: duplicate x-value at indices " + std::to_string(i) + " and " + std::to_string(j));

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

		/// @brief Returns the linear polynomial a*x + b.
		/// @param a Coefficient of x
		/// @param b Constant term
		static Polynom Linear(const CoefT& a, const CoefT& b) {
			return Polynom(std::initializer_list<CoefT>{b, a});
		}

		/// @brief Returns degree of polynomial (preferred API).
		int  degree() const noexcept { return (int)_vecCoef.size() - 1; }
		/// @brief Sets degree of polynomial (resizes coefficient vector).
		/// @param newDeg New degree
		void SetDegree(int newDeg) { _vecCoef.resize(newDeg + 1); }
		/// @brief Checks if polynomial is null (no coefficients).
		bool isNull() const noexcept { return _vecCoef.size() == 0; }
		/// @brief Removes trailing zero coefficients.
		void Reduce() { while (!_vecCoef.empty() && _vecCoef.back() == CoefT(0)) _vecCoef.pop_back(); }
		/// @brief Removes trailing coefficients smaller than eps in absolute value.
		/// @param eps Tolerance threshold for near-zero coefficients
		void Reduce(CoefT eps) { while (!_vecCoef.empty() && std::abs(_vecCoef.back()) < eps) _vecCoef.pop_back(); }
		
		/// @brief Returns leading coefficient.
		CoefT leadingTerm() const noexcept { return _vecCoef.empty() ? CoefT(0) : _vecCoef.back(); }
		/// @brief Returns constant term.
		CoefT constantTerm() const noexcept { return _vecCoef.empty() ? CoefT(0) : _vecCoef[0]; }

		/// @brief Coefficient access (const).
		/// @param i Coefficient index (0 = constant term)
		const CoefT& operator[] (int i) const noexcept { return _vecCoef[i]; }
		/// @brief Coefficient access (non-const).
		/// @param i Coefficient index (0 = constant term)
		CoefT& operator[] (int i) noexcept { return _vecCoef[i]; }

		///////////////////////////            Iterators              ///////////////////////////

		using iterator = typename std::vector<CoefT>::iterator;
		using const_iterator = typename std::vector<CoefT>::const_iterator;

		iterator begin() { return _vecCoef.begin(); }
		iterator end() { return _vecCoef.end(); }
		const_iterator begin() const { return _vecCoef.begin(); }
		const_iterator end() const { return _vecCoef.end(); }
		const_iterator cbegin() const { return _vecCoef.cbegin(); }
		const_iterator cend() const { return _vecCoef.cend(); }

		//////////////////   Additional convenience methods   /////////////////
		
		/// @brief Alias for operator() - evaluates polynomial at x.
		/// @param x Evaluation point
		FieldT eval(const FieldT& x) const { return operator()(x); }
		
		/// @brief Returns coefficient vector (const reference).
		const std::vector<CoefT>& coefficients() const { return _vecCoef; }
		
		/// @brief Returns coefficient at index i.
		/// @param i Coefficient index (0 = constant term)
		CoefT coeff(int i) const { return (i >= 0 && i < (int)_vecCoef.size()) ? _vecCoef[i] : CoefT(0); }
		
		/// @brief Returns string representation with default formatting.
		std::string to_string() const { return to_string(10, 6); }

		///////////////////////////            Operators              ///////////////////////////
		/// @brief Copy assignment.
		Polynom& operator=(const Polynom& Copy) { _vecCoef = Copy._vecCoef; return *this; }

		/// @brief Evaluates polynomial at given value using Horner's method.
		/// @param x Evaluation point
		/// @return P(x)
		FieldT operator() (const FieldT& x) const {
			int j = degree();
			FieldT p = _vecCoef[j] * FieldT(1);

			while (j > 0)
				p = p * x + _vecCoef[--j];
			return p;
		}

		/// @brief Checks exact equality with another polynomial (ignoring trailing zeros).
		/// @param b Polynomial to compare
		/// @note For tolerance-based comparison, use IsEqualTo()
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

		/// @brief Checks tolerance-based equality with another polynomial.
		/// @param b Polynomial to compare
		/// @param eps Tolerance for coefficient comparison
		/// @return True if all coefficients are within eps of each other
		template<typename T = CoefT>
		typename std::enable_if<std::is_floating_point<T>::value, bool>::type
		IsEqualTo(const Polynom& b, T eps = Defaults::DefaultTolerance) const
		{
			int maxSize = (int)std::max(_vecCoef.size(), b._vecCoef.size());
			
			for (int i = 0; i < maxSize; i++)
			{
				T a_coef = (i < (int)_vecCoef.size()) ? _vecCoef[i] : T(0);
				T b_coef = (i < (int)b._vecCoef.size()) ? b._vecCoef[i] : T(0);
				
				if (std::abs(a_coef - b_coef) > eps)
					return false;
			}
			return true;
		}

		/// @brief Equality operator.
		bool operator==(const Polynom& b) const
		{
			return IsEqual(b);
		}

		/// @brief Inequality operator.
		bool operator!=(const Polynom& b) const
		{
			return !IsEqual(b);
		}

		/// @brief Unary negation operator.
		Polynom operator-() const
		{
			Polynom result;
			result._vecCoef.resize(_vecCoef.size());
			for (size_t i = 0; i < _vecCoef.size(); i++)
				result._vecCoef[i] = -_vecCoef[i];
			return result;
		}

		/// @brief Polynomial addition.
		/// @param b Polynomial to add
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

		/// @brief Polynomial subtraction.
		/// @param b Polynomial to subtract
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

		/// @brief Polynomial multiplication.
		/// @param b Polynomial to multiply
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

		/// @brief Addition by scalar (adds to constant term).
		/// @param scalar Value to add
		Polynom operator+(const CoefT& scalar) const {
			Polynom result = *this;
			if (result._vecCoef.empty())
				result._vecCoef.resize(1, CoefT(0));
			result._vecCoef[0] += scalar;
			return result;
		}
		/// @brief Scalar + polynomial addition.
		friend Polynom operator+(const CoefT& scalar, const Polynom& poly) {
			return poly + scalar;
		}

		/// @brief Subtraction by scalar (subtracts from constant term).
		/// @param scalar Value to subtract
		Polynom operator-(const CoefT& scalar) const {
			Polynom result = *this;
			if (result._vecCoef.empty())
				result._vecCoef.resize(1, CoefT(0));
			result._vecCoef[0] -= scalar;
			return result;
		}
		/// @brief Scalar - polynomial subtraction.
		friend Polynom operator-(const CoefT& scalar, const Polynom& poly) {
			Polynom result = -poly;
			if (result._vecCoef.empty())
				result._vecCoef.resize(1, CoefT(0));
			result._vecCoef[0] += scalar;
			return result;
		}

		/// @brief Multiplication by scalar.
		/// @param scalar Multiplier
		Polynom operator*(const CoefT& scalar) const {
			Polynom result = *this;
			for (auto& coef : result._vecCoef)
				coef *= scalar;
			return result;
		}
		/// @brief Scalar * polynomial multiplication.
		friend Polynom operator*(const CoefT& scalar, const Polynom& poly) {
			return poly * scalar;
		}

		/// @brief Division by scalar.
		/// @param scalar Divisor
		Polynom operator/(const CoefT& scalar) const {
			Polynom result = *this;
			for (auto& coef : result._vecCoef)
				coef /= scalar;
			return result;
		}

		/// @brief Polynomial division with quotient and remainder.
		/// @param u Dividend polynomial
		/// @param v Divisor polynomial
		/// @param qout Quotient (output)
		/// @param rout Remainder (output)
		static void poldiv(const Polynom& u, const Polynom& v, Polynom& qout, Polynom& rout)
		{
			int k, j, n = u.degree(), nv = v.degree();

			// find real degree of v
			while (nv >= 0 && v._vecCoef[nv] == CoefT(0))
				nv--;

			if (nv < 0)
				throw DivisionByZeroError("poldiv: divide by zero polynomial");

			Polynom r = u;
			Polynom q(u.degree());
			for (k = n - nv; k >= 0; k--)
			{
				q[k] = r[nv + k] / v[nv];

				for (j = nv + k - 1; j >= k; j--)
					r[j] -= q[k] * v[j - k];
			}
			for (j = nv; j <= n; j++)
				r[j] = CoefT(0);


			int nq = q.degree();
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
		/// @brief Evaluates polynomial and its first n derivatives at point x.
		/// @param x Evaluation point
		/// @param pd Output vector: pd[0]=P(x), pd[1]=P'(x), ..., pd[n]=P^(n)(x)
		void Derive(const Real x, Vector<Real>& pd)
		{
			int  nnd, j, i;
			int  nc = degree();
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

		/// @brief Returns the first derivative as a new polynomial (preferred API).
		Polynom derivative() const {
			Polynom result;
			int deg = degree();
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
		
		/// @brief Returns the indefinite integral (antiderivative) with constant = 0 (preferred API).
		Polynom integral() const {
			Polynom result;
			int deg = degree();
			result._vecCoef.resize(deg + 2, CoefT(0));
			for (int i = 0; i <= deg; ++i) {
				result._vecCoef[i + 1] = _vecCoef[i] / CoefT(i + 1);
			}
			return result;
		}
		
		///////////////////////////               I/O                 ///////////////////////////
		/// @brief Converts polynomial to string representation.
		/// @param width Field width
		/// @param precision Decimal precision
		std::string to_string(int width, int precision) const
		{
			std::stringstream str;

			Print(str, width, precision);

			return str.str();
		}

		/// @brief Prints polynomial to stream with default formatting.
		std::ostream& Print(std::ostream& stream) const
		{
			return Print(stream, 10, 6); // Default width and precision
		}

		/// @brief Prints polynomial to stream with custom formatting.
		/// @param stream Output stream
		/// @param width Field width
		/// @param precision Decimal precision
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

		/// @brief Stream output operator.
		friend std::ostream& operator<<(std::ostream& stream, const Polynom& poly)
		{
			return poly.Print(stream);
		}
	};

	/// @brief Real polynomial that implements IRealFunction interface.
	class PolynomRealFunc : public Polynom<Real>, public IRealFunction
	{
	public:
		/// @brief Default constructor.
		PolynomRealFunc() {}
		/// @brief Constructs polynomial of degree n.
		PolynomRealFunc(int n) { _vecCoef.resize(n + 1); }
		/// @brief Constructs from coefficient vector.
		PolynomRealFunc(const std::vector<Real>& vecCoef) : Polynom(vecCoef) {}
		/// @brief Constructs from initializer list.
		PolynomRealFunc(std::initializer_list<Real> list) : Polynom(list) {}
		/// @brief Constructs from existing polynomial.
		PolynomRealFunc(const Polynom& Copy) : Polynom(Copy) {}
		~PolynomRealFunc() {}

		/// @brief Evaluates polynomial at x (IRealFunction interface).
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