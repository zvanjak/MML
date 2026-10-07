///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FiniteField.h                                                       ///
///  Description: Prime-power extension field value types and finite adapters         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FINITE_FIELD_H
#define MML_FINITE_FIELD_H

#include <mml/MMLExceptions.h>
#include <mml/base/Algebra/ModInt.h>
#include <mml/base/Algebra/Polynomial.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <vector>

namespace MML
{
	namespace Algebra
	{
		namespace Detail
		{
			constexpr std::uint64_t IntegerPower(std::uint64_t base, int exponent) noexcept
			{
				std::uint64_t result = 1;
				for (int index = 0; index < exponent; index++)
					result *= base;
				return result;
			}
		}

		/// @brief Element of GF(P^N) represented modulo a supplied monic polynomial.
		/// @tparam ModulusPolynomial Provider with static modulus() returning
		///         Polynomial<PrimeFieldElement<P>> of degree N.
		template<int P, int N, class ModulusPolynomial>
		class FiniteFieldElement
		{
			static_assert(IsPrime<P>::value, "FiniteFieldElement characteristic must be prime");
			static_assert(N > 1, "Use PrimeFieldElement<P> for degree-one fields");

		public:
			using coefficient_type = PrimeFieldElement<P>;
			using polynomial_type = Polynomial<coefficient_type>;
			static constexpr int characteristic = P;
			static constexpr int extension_degree = N;
			static constexpr std::uint64_t field_order = Detail::IntegerPower(P, N);

		private:
			std::array<coefficient_type, N> _coefficients{};

			static const polynomial_type& modulus_polynomial()
			{
				static const polynomial_type modulus = [] {
					polynomial_type value = ModulusPolynomial::modulus();
					if (value.degree() != N || value.leading_coefficient() != coefficient_type(1))
						throw ArgumentError("Finite-field modulus must be monic with degree N");
					return value;
				}();
				return modulus;
			}

			static FiniteFieldElement from_reduced_polynomial(const polynomial_type& polynomial)
			{
				FiniteFieldElement result;
				for (int degree = 0; degree < N; degree++)
					result._coefficients[static_cast<std::size_t>(degree)] = polynomial.coefficient(degree);
				return result;
			}

		public:
			FiniteFieldElement() = default;
			FiniteFieldElement(std::int64_t constant) { _coefficients[0] = coefficient_type(constant); }
			explicit FiniteFieldElement(const std::array<coefficient_type, N>& coefficients)
				: _coefficients(coefficients) { (void)modulus_polynomial(); }
			FiniteFieldElement(std::initializer_list<coefficient_type> coefficients)
			{
				if (coefficients.size() > static_cast<std::size_t>(N))
					throw ArgumentError("Too many finite-field coefficients");
				std::copy(coefficients.begin(), coefficients.end(), _coefficients.begin());
				(void)modulus_polynomial();
			}

			static const polynomial_type& modulus() { return modulus_polynomial(); }
			const std::array<coefficient_type, N>& coefficients() const noexcept { return _coefficients; }
			coefficient_type coefficient(int degree) const
			{
				return degree >= 0 && degree < N ? _coefficients[static_cast<std::size_t>(degree)] : coefficient_type(0);
			}
			bool is_zero() const noexcept
			{
				for (const auto& coefficient : _coefficients)
					if (coefficient != coefficient_type(0))
						return false;
				return true;
			}

			polynomial_type polynomial() const
			{
				return polynomial_type(std::vector<coefficient_type>(_coefficients.begin(), _coefficients.end()));
			}

			FiniteFieldElement& operator+=(const FiniteFieldElement& other) noexcept
			{
				for (int index = 0; index < N; index++)
					_coefficients[static_cast<std::size_t>(index)] += other._coefficients[static_cast<std::size_t>(index)];
				return *this;
			}

			FiniteFieldElement& operator-=(const FiniteFieldElement& other) noexcept
			{
				for (int index = 0; index < N; index++)
					_coefficients[static_cast<std::size_t>(index)] -= other._coefficients[static_cast<std::size_t>(index)];
				return *this;
			}

			FiniteFieldElement& operator*=(const FiniteFieldElement& other)
			{
				*this = from_reduced_polynomial((polynomial() * other.polynomial()) % modulus_polynomial());
				return *this;
			}

			FiniteFieldElement pow(std::uint64_t exponent) const
			{
				FiniteFieldElement result(1);
				FiniteFieldElement factor(*this);
				while (exponent > 0) {
					if ((exponent & 1U) != 0)
						result *= factor;
					factor *= factor;
					exponent >>= 1U;
				}
				return result;
			}

			FiniteFieldElement inverse() const
			{
				if (is_zero())
					throw DivisionByZeroError("Zero has no finite-field inverse");
				return pow(field_order - 2);
			}

			FiniteFieldElement& operator/=(const FiniteFieldElement& other) { return *this *= other.inverse(); }

			friend FiniteFieldElement operator+(FiniteFieldElement left, const FiniteFieldElement& right) { return left += right; }
			friend FiniteFieldElement operator-(FiniteFieldElement left, const FiniteFieldElement& right) { return left -= right; }
			friend FiniteFieldElement operator-(FiniteFieldElement value) { return FiniteFieldElement() - value; }
			friend FiniteFieldElement operator*(FiniteFieldElement left, const FiniteFieldElement& right) { return left *= right; }
			friend FiniteFieldElement operator/(FiniteFieldElement left, const FiniteFieldElement& right) { return left /= right; }
			friend bool operator==(const FiniteFieldElement& left, const FiniteFieldElement& right) noexcept
			{
				return left._coefficients == right._coefficients;
			}
			friend bool operator!=(const FiniteFieldElement& left, const FiniteFieldElement& right) noexcept { return !(left == right); }
		};

		template<int P, int N, class ModulusPolynomial>
		class ExtensionField
		{
		public:
			using element_type = FiniteFieldElement<P, N, ModulusPolynomial>;
			using scalar_type = element_type;

		private:
			std::vector<element_type> _elements;

		public:
			ExtensionField()
			{
				static_assert(element_type::field_order <= static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max()),
					"Extension field is too large to enumerate");
				_elements.reserve(static_cast<std::size_t>(element_type::field_order));
				for (std::uint64_t encoded = 0; encoded < element_type::field_order; encoded++) {
					std::uint64_t remaining = encoded;
					std::array<typename element_type::coefficient_type, N> coefficients{};
					for (int degree = 0; degree < N; degree++) {
						coefficients[static_cast<std::size_t>(degree)] = typename element_type::coefficient_type(remaining % P);
						remaining /= P;
					}
					_elements.emplace_back(coefficients);
				}
			}

			const std::vector<element_type>& elements() const noexcept { return _elements; }
			element_type zero() const { return element_type(0); }
			element_type one() const { return element_type(1); }
			element_type add(const element_type& left, const element_type& right) const { return left + right; }
			element_type negate(const element_type& element) const { return -element; }
			element_type multiply(const element_type& left, const element_type& right) const { return left * right; }
			element_type inverse(const element_type& element) const { return element.inverse(); }
			std::size_t order() const noexcept { return _elements.size(); }
		};
	}
}

#endif // MML_FINITE_FIELD_H