///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Polynomial.h                                                        ///
///  Description: Canonical exact polynomials over coefficient fields                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ALGEBRA_POLYNOMIAL_H
#define MML_ALGEBRA_POLYNOMIAL_H

#include <mml/MMLExceptions.h>

#include <algorithm>
#include <initializer_list>
#include <utility>
#include <vector>

namespace MML
{
	namespace Algebra
	{
		/// @brief Exact polynomial with coefficients stored in ascending power order.
		template<class Coeff>
		class Polynomial
		{
			std::vector<Coeff> _coefficients;

			void normalize()
			{
				while (!_coefficients.empty() && _coefficients.back() == Coeff(0))
					_coefficients.pop_back();
			}

		public:
			using coefficient_type = Coeff;

			Polynomial() = default;
			Polynomial(std::initializer_list<Coeff> coefficients) : _coefficients(coefficients) { normalize(); }
			explicit Polynomial(std::vector<Coeff> coefficients)
				: _coefficients(std::move(coefficients)) { normalize(); }

			static Polynomial zero() { return Polynomial(); }
			static Polynomial one() { return Polynomial({Coeff(1)}); }
			static Polynomial monomial(int degree, Coeff coefficient = Coeff(1))
			{
				if (degree < 0)
					throw ArgumentError("Polynomial monomial degree cannot be negative");
				std::vector<Coeff> coefficients(static_cast<std::size_t>(degree + 1), Coeff(0));
				coefficients[static_cast<std::size_t>(degree)] = coefficient;
				return Polynomial(std::move(coefficients));
			}

			bool is_zero() const noexcept { return _coefficients.empty(); }
			int degree() const noexcept { return static_cast<int>(_coefficients.size()) - 1; }
			Coeff leading_coefficient() const { return is_zero() ? Coeff(0) : _coefficients.back(); }
			Coeff coefficient(int degree) const
			{
				return degree >= 0 && static_cast<std::size_t>(degree) < _coefficients.size()
					? _coefficients[static_cast<std::size_t>(degree)] : Coeff(0);
			}
			const std::vector<Coeff>& coefficients() const noexcept { return _coefficients; }

			Coeff evaluate(Coeff value) const
			{
				Coeff result(0);
				for (auto coefficient = _coefficients.rbegin(); coefficient != _coefficients.rend(); ++coefficient)
					result = result * value + *coefficient;
				return result;
			}

			Polynomial& operator+=(const Polynomial& other)
			{
				if (_coefficients.size() < other._coefficients.size())
					_coefficients.resize(other._coefficients.size(), Coeff(0));
				for (std::size_t index = 0; index < other._coefficients.size(); index++)
					_coefficients[index] += other._coefficients[index];
				normalize();
				return *this;
			}

			Polynomial& operator-=(const Polynomial& other)
			{
				if (_coefficients.size() < other._coefficients.size())
					_coefficients.resize(other._coefficients.size(), Coeff(0));
				for (std::size_t index = 0; index < other._coefficients.size(); index++)
					_coefficients[index] -= other._coefficients[index];
				normalize();
				return *this;
			}

			Polynomial& operator*=(const Polynomial& other)
			{
				if (is_zero() || other.is_zero()) {
					_coefficients.clear();
					return *this;
				}
				std::vector<Coeff> product(_coefficients.size() + other._coefficients.size() - 1, Coeff(0));
				for (std::size_t left = 0; left < _coefficients.size(); left++)
					for (std::size_t right = 0; right < other._coefficients.size(); right++)
						product[left + right] += _coefficients[left] * other._coefficients[right];
				_coefficients = std::move(product);
				normalize();
				return *this;
			}

			Polynomial& operator*=(Coeff scalar)
			{
				for (auto& coefficient : _coefficients)
					coefficient *= scalar;
				normalize();
				return *this;
			}

			Polynomial monic() const
			{
				if (is_zero())
					return *this;
				return *this * leading_coefficient().inverse();
			}

			std::pair<Polynomial, Polynomial> divmod(const Polynomial& divisor) const
			{
				if (divisor.is_zero())
					throw DivisionByZeroError("Polynomial division by zero");
				Polynomial quotient;
				Polynomial remainder(*this);
				std::vector<Coeff> quotientCoefficients(
					remainder.degree() >= divisor.degree()
						? static_cast<std::size_t>(remainder.degree() - divisor.degree() + 1) : 0,
					Coeff(0));

				while (!remainder.is_zero() && remainder.degree() >= divisor.degree()) {
					const int degreeDifference = remainder.degree() - divisor.degree();
					const Coeff factor = remainder.leading_coefficient() / divisor.leading_coefficient();
					quotientCoefficients[static_cast<std::size_t>(degreeDifference)] += factor;
					remainder -= divisor * monomial(degreeDifference, factor);
				}
				quotient = Polynomial(std::move(quotientCoefficients));
				return {quotient, remainder};
			}

			friend Polynomial operator+(Polynomial left, const Polynomial& right) { return left += right; }
			friend Polynomial operator-(Polynomial left, const Polynomial& right) { return left -= right; }
			friend Polynomial operator-(Polynomial polynomial) { return polynomial *= Coeff(-1); }
			friend Polynomial operator*(Polynomial left, const Polynomial& right) { return left *= right; }
			friend Polynomial operator*(Polynomial polynomial, Coeff scalar) { return polynomial *= scalar; }
			friend Polynomial operator*(Coeff scalar, Polynomial polynomial) { return polynomial *= scalar; }
			friend Polynomial operator/(const Polynomial& dividend, const Polynomial& divisor)
			{
				return dividend.divmod(divisor).first;
			}
			friend Polynomial operator%(const Polynomial& dividend, const Polynomial& divisor)
			{
				return dividend.divmod(divisor).second;
			}

			friend bool operator==(const Polynomial& left, const Polynomial& right) noexcept
			{
				return left._coefficients == right._coefficients;
			}
			friend bool operator!=(const Polynomial& left, const Polynomial& right) noexcept { return !(left == right); }
		};
	}
}

#endif // MML_ALGEBRA_POLYNOMIAL_H