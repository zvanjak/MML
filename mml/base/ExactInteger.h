///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///  File:        ExactInteger.h                                                      ///
///  Description: Prechecked signed arithmetic and fixed two-limb rational scratch     ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                         ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_EXACT_INTEGER_H
#define MML_EXACT_INTEGER_H

#include <cstdint>
#include <limits>
#include <mml/MMLExceptions.h>

namespace MML::ExactInteger
{
	// Unsigned conversion is modulo 2^64, including for LLONG_MIN.
	inline std::uint64_t Magnitude(long long value)
	{
		const auto u = static_cast<std::uint64_t>(value);
		return value < 0 ? std::uint64_t{0} - u : u;
	}

	inline std::uint64_t Gcd(std::uint64_t a, std::uint64_t b)
	{
		while (b != 0) { const auto r = a % b; a = b; b = r; }
		return a;
	}

	inline long long FromMagnitude(std::uint64_t value, bool negative)
	{
		constexpr auto max = static_cast<std::uint64_t>((std::numeric_limits<long long>::max)());
		if (value <= max) {
			const auto signedValue = static_cast<long long>(value);
			return negative ? -signedValue : signedValue;
		}
		if (negative && value == max + 1) return (std::numeric_limits<long long>::min)();
		throw DomainError("Exact integer: result exceeds signed 64-bit range");
	}

	inline long long Add(long long a, long long b)
	{
		constexpr auto max = (std::numeric_limits<long long>::max)();
		constexpr auto min = (std::numeric_limits<long long>::min)();
		if ((b > 0 && a > max - b) || (b < 0 && a < min - b))
			throw DomainError("Exact integer: addition exceeds signed 64-bit range");
		return a + b;
	}

	inline long long Subtract(long long a, long long b)
	{
		constexpr auto max = (std::numeric_limits<long long>::max)();
		constexpr auto min = (std::numeric_limits<long long>::min)();
		if ((b > 0 && a < min + b) || (b < 0 && a > max + b))
			throw DomainError("Exact integer: subtraction exceeds signed 64-bit range");
		return a - b;
	}

	inline long long Multiply(long long a, long long b)
	{
		const auto ua = Magnitude(a), ub = Magnitude(b);
		const bool negative = (a < 0) != (b < 0);
		const auto limit = static_cast<std::uint64_t>((std::numeric_limits<long long>::max)())
			+ (negative ? 1ULL : 0ULL);
		if (ub != 0 && ua > limit / ub)
			throw DomainError("Exact integer: product exceeds signed 64-bit range");
		return FromMagnitude(ua * ub, negative);
	}

	// Fixed two-limb scratch space for rational cross-products, not arbitrary precision.
	// Callers keep sums within 128 bits; all limb arithmetic is unsigned.
	struct WideUnsigned
	{
		std::uint64_t high = 0, low = 0;

		static WideUnsigned Product(std::uint64_t a, std::uint64_t b)
		{
			constexpr std::uint64_t mask = 0xffffffffULL;
			const auto a0 = a & mask, a1 = a >> 32;
			const auto b0 = b & mask, b1 = b >> 32;
			const auto w0 = a0 * b0;
			const auto t = a1 * b0 + (w0 >> 32);
			const auto w1 = (t & mask) + a0 * b1;
			return {a1 * b1 + (t >> 32) + (w1 >> 32), (w1 << 32) | (w0 & mask)};
		}

		bool operator<(const WideUnsigned& r) const
		{
			return high < r.high || (high == r.high && low < r.low);
		}

		WideUnsigned operator+(const WideUnsigned& r) const
		{
			const auto sum = low + r.low;
			return {high + r.high + (sum < low ? 1ULL : 0ULL), sum};
		}

		// Precondition: *this >= r.
		WideUnsigned operator-(const WideUnsigned& r) const
		{
			return {high - r.high - (low < r.low ? 1ULL : 0ULL), low - r.low};
		}

		WideUnsigned Divide(std::uint64_t divisor, std::uint64_t& remainder) const
		{
			if (divisor == 0) throw DomainError("Exact integer: zero divisor");
			if (high == 0) {
				remainder = low % divisor;
				return {0, low / divisor};
			}
			WideUnsigned quotient;
			remainder = 0;
			for (int bit = 127; bit >= 0; --bit) {
				const bool carry = (remainder >> 63) != 0;
				const auto digit = bit >= 64 ? (high >> (bit - 64)) & 1 : (low >> bit) & 1;
				remainder = (remainder << 1) | digit;
				if (carry || remainder >= divisor) {
					remainder -= divisor;
					if (bit >= 64) quotient.high |= std::uint64_t{1} << (bit - 64);
					else quotient.low |= std::uint64_t{1} << bit;
				}
			}
			return quotient;
		}
	};
}

#endif // MML_EXACT_INTEGER_H
