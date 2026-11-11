///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Rational.h                                                          ///
///  Description: Exact rational number backed by 64-bit integers                     ///
///               Always stored in lowest terms with a positive denominator          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_RATIONAL_H
#define MML_RATIONAL_H

#include <ostream>

#include <mml/MMLBase.h>
#include <limits>
#include <mml/base/NumberTheory.h>

namespace MML
{
	/// @brief Exact rational number p/q backed by 64-bit integers.
	/// @details The value is always kept in lowest terms with a strictly positive
	///          denominator, so equality is a plain member comparison. Arithmetic is
	///          exact, including LLONG_MIN numerators; DomainError means the reduced
	///          numerator or positive denominator cannot fit long long.
	class Rational
	{
		long long num_ = 0;   // numerator, carries the sign
		long long den_ = 1;   // denominator, always > 0

		void normalize()
		{
			if (den_ == 0) throw DivisionByZeroError("Rational: zero denominator");
			auto n = ExactInteger::Magnitude(num_), d = ExactInteger::Magnitude(den_);
			const auto g = ExactInteger::Gcd(n, d);
			const bool negative = (num_ < 0) != (den_ < 0);
			num_ = ExactInteger::FromMagnitude(n / g, negative);
			den_ = ExactInteger::FromMagnitude(d / g, false);
		}

		static Rational fromMagnitudes(std::uint64_t n, std::uint64_t d, bool negative)
		{
			const auto g = ExactInteger::Gcd(n, d);
			Rational result;
			result.num_ = ExactInteger::FromMagnitude(n / g, negative);
			result.den_ = ExactInteger::FromMagnitude(d / g, false);
			return result;
		}

		Rational add(const Rational& r, bool subtract) const
		{
			using ExactInteger::WideUnsigned;
			const auto b = static_cast<std::uint64_t>(den_);
			const auto d = static_cast<std::uint64_t>(r.den_);
			const auto g = ExactInteger::Gcd(b, d);
			const auto left = WideUnsigned::Product(ExactInteger::Magnitude(num_), d / g);
			const auto right = WideUnsigned::Product(ExactInteger::Magnitude(r.num_), b / g);
			const bool leftNegative = num_ < 0;
			const bool rightNegative = (r.num_ < 0) != subtract;
			const bool negative = leftNegative == rightNegative || !(left < right)
				? leftNegative : rightNegative;
			const auto sum = leftNegative == rightNegative ? left + right
				: (left < right ? right - left : left - right);
			// Only a factor of g can still cancel. Keep the at-most-127-bit
			// numerator exact until that cancellation, even without __int128.
			std::uint64_t remainder = 0;
			sum.Divide(g, remainder);
			const auto h = ExactInteger::Gcd(remainder, g);
			const auto numerator = sum.Divide(h, remainder);
			if (numerator.high != 0)
				throw DomainError("Rational: reduced numerator exceeds signed 64-bit range");
			Rational result;
			result.num_ = ExactInteger::FromMagnitude(numerator.low, negative);
			result.den_ = ExactInteger::Multiply(static_cast<long long>(b / g),
				static_cast<long long>(d / h));
			return result;
		}

		Rational multiply(const Rational& r, bool divide) const
		{
			auto a = ExactInteger::Magnitude(num_), b = static_cast<std::uint64_t>(den_);
			auto c = divide ? static_cast<std::uint64_t>(r.den_) : ExactInteger::Magnitude(r.num_);
			auto d = divide ? ExactInteger::Magnitude(r.num_) : static_cast<std::uint64_t>(r.den_);
			if (d == 0) throw DivisionByZeroError("Rational: division by zero");
			const auto g1 = ExactInteger::Gcd(a, d), g2 = ExactInteger::Gcd(c, b);
			a /= g1; d /= g1; c /= g2; b /= g2;
			const auto n = ExactInteger::WideUnsigned::Product(a, c);
			const auto den = ExactInteger::WideUnsigned::Product(b, d);
			if (n.high != 0 || den.high != 0)
				throw DomainError("Rational: reduced result exceeds signed 64-bit range");
			return fromMagnitudes(n.low, den.low, (num_ < 0) != (r.num_ < 0));
		}

	public:
		Rational() = default;
		Rational(long long n) : num_(n), den_(1) {}
		Rational(long long n, long long d) : num_(n), den_(d) { normalize(); }

		long long Num() const { return num_; }
		long long Den() const { return den_; }
		Real ToReal() const { return static_cast<Real>(num_) / static_cast<Real>(den_); }

		Rational operator-() const { return Rational(ExactInteger::Subtract(0, num_), den_); }

		Rational operator+(const Rational& r) const
		{
			return add(r, false);
		}

		Rational operator-(const Rational& r) const { return add(r, true); }

		Rational operator*(const Rational& r) const
		{
			return multiply(r, false);
		}

		Rational operator/(const Rational& r) const
		{
			return multiply(r, true);
		}

		Rational& operator+=(const Rational& r) { *this = *this + r; return *this; }
		Rational& operator-=(const Rational& r) { *this = *this - r; return *this; }
		Rational& operator*=(const Rational& r) { *this = *this * r; return *this; }
		Rational& operator/=(const Rational& r) { *this = *this / r; return *this; }

		bool operator==(const Rational& r) const { return num_ == r.num_ && den_ == r.den_; }
		bool operator!=(const Rational& r) const { return !(*this == r); }

		bool operator<(const Rational& r) const
		{
			if ((num_ < 0) != (r.num_ < 0)) return num_ < 0;
			// Compare magnitudes exactly in fixed-width unsigned scratch space.
			const auto left = ExactInteger::WideUnsigned::Product(
				ExactInteger::Magnitude(num_), static_cast<std::uint64_t>(r.den_));
			const auto right = ExactInteger::WideUnsigned::Product(
				ExactInteger::Magnitude(r.num_), static_cast<std::uint64_t>(den_));
			return num_ < 0 ? right < left : left < right;
		}
		bool operator>(const Rational& r) const { return r < *this; }
		bool operator<=(const Rational& r) const { return !(r < *this); }
		bool operator>=(const Rational& r) const { return !(*this < r); }

		friend std::ostream& operator<<(std::ostream& os, const Rational& r)
		{
			if (r.den_ == 1) os << r.num_;
			else os << r.num_ << "/" << r.den_;
			return os;
		}
	};
} // namespace MML

#endif // MML_RATIONAL_H
