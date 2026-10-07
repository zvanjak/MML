///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ModInt.h                                                            ///
///  Description: Exact modular integers and prime-field scalar elements              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_MOD_INT_H
#define MML_MOD_INT_H

#include <mml/MMLExceptions.h>

#include <cstdint>
#include <type_traits>

namespace MML
{
	namespace Algebra
	{
		constexpr bool IsPrimeValue(int value) noexcept
		{
			if (value < 2)
				return false;
			if (value % 2 == 0)
				return value == 2;
			for (int divisor = 3; static_cast<std::int64_t>(divisor) * divisor <= value; divisor += 2)
				if (value % divisor == 0)
					return false;
			return true;
		}

		template<int Value>
		struct IsPrime : std::bool_constant<IsPrimeValue(Value)>
		{
		};

		namespace Detail
		{
			inline std::int64_t ExtendedGcd(std::int64_t left, std::int64_t right,
				std::int64_t& leftCoefficient, std::int64_t& rightCoefficient)
			{
				if (right == 0) {
					leftCoefficient = 1;
					rightCoefficient = 0;
					return left;
				}

				std::int64_t nextLeftCoefficient = 0;
				std::int64_t nextRightCoefficient = 0;
				const std::int64_t divisor = ExtendedGcd(right, left % right,
					nextLeftCoefficient, nextRightCoefficient);
				leftCoefficient = nextRightCoefficient;
				rightCoefficient = nextLeftCoefficient - (left / right) * nextRightCoefficient;
				return divisor;
			}
		}

		/// @brief Exact integer arithmetic modulo a positive compile-time modulus.
		template<int Modulus>
		class ModInt
		{
			static_assert(Modulus > 1, "ModInt modulus must be greater than one");
			int _value = 0;

			static constexpr int normalize(std::int64_t value) noexcept
			{
				const std::int64_t remainder = value % Modulus;
				return static_cast<int>(remainder < 0 ? remainder + Modulus : remainder);
			}

		public:
			using value_type = int;
			static constexpr int modulus = Modulus;
			static constexpr bool is_field = IsPrime<Modulus>::value;

			constexpr ModInt() noexcept = default;
			constexpr ModInt(std::int64_t value) noexcept : _value(normalize(value)) { }

			constexpr int value() const noexcept { return _value; }
			explicit constexpr operator int() const noexcept { return _value; }

			constexpr ModInt operator+() const noexcept { return *this; }
			constexpr ModInt operator-() const noexcept { return ModInt(-static_cast<std::int64_t>(_value)); }

			constexpr ModInt& operator+=(ModInt other) noexcept
			{
				_value = normalize(static_cast<std::int64_t>(_value) + other._value);
				return *this;
			}

			constexpr ModInt& operator-=(ModInt other) noexcept
			{
				_value = normalize(static_cast<std::int64_t>(_value) - other._value);
				return *this;
			}

			constexpr ModInt& operator*=(ModInt other) noexcept
			{
				_value = normalize(static_cast<std::int64_t>(_value) * other._value);
				return *this;
			}

			ModInt inverse() const
			{
				if (_value == 0)
					throw DivisionByZeroError("Zero has no modular inverse");
				std::int64_t coefficient = 0;
				std::int64_t unused = 0;
				const std::int64_t divisor = Detail::ExtendedGcd(_value, Modulus, coefficient, unused);
				if (divisor != 1)
					throw DomainError("Modular integer is not a unit");
				return ModInt(coefficient);
			}

			ModInt& operator/=(ModInt other)
			{
				return *this *= other.inverse();
			}

			constexpr ModInt pow(std::uint64_t exponent) const noexcept
			{
				ModInt result(1);
				ModInt factor(*this);
				while (exponent > 0) {
					if ((exponent & 1U) != 0)
						result *= factor;
					factor *= factor;
					exponent >>= 1U;
				}
				return result;
			}

			friend constexpr ModInt operator+(ModInt left, ModInt right) noexcept { return left += right; }
			friend constexpr ModInt operator-(ModInt left, ModInt right) noexcept { return left -= right; }
			friend constexpr ModInt operator*(ModInt left, ModInt right) noexcept { return left *= right; }
			friend ModInt operator/(ModInt left, ModInt right) { return left /= right; }

			friend constexpr bool operator==(ModInt left, ModInt right) noexcept { return left._value == right._value; }
			friend constexpr bool operator!=(ModInt left, ModInt right) noexcept { return !(left == right); }
		};

		template<int Prime>
		struct PrimeFieldType
		{
			static_assert(IsPrime<Prime>::value, "PrimeFieldElement modulus must be prime");
			using type = ModInt<Prime>;
		};

		template<int Prime>
		using PrimeFieldElement = typename PrimeFieldType<Prime>::type;
	}
}

#endif // MML_MOD_INT_H