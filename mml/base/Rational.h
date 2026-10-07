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
#include <mml/base/NumberTheory.h>

namespace MML
{
	/// @brief Exact rational number p/q backed by 64-bit integers.
	/// @details The value is always kept in lowest terms with a strictly positive
	///          denominator, so equality is a plain member comparison. Arithmetic is
	///          overflow-checked and throws DomainError if a 64-bit result would overflow.
	class Rational
	{
		long long num_ = 0;   // numerator, carries the sign
		long long den_ = 1;   // denominator, always > 0

		static long long mulChecked(long long a, long long b)
		{
			if (a == 0 || b == 0) return 0;
			long long r = a * b;
			if (r / a != b) throw DomainError("Rational: 64-bit integer overflow");
			return r;
		}

		static long long addChecked(long long a, long long b)
		{
			long long r = a + b;
			if (((a ^ r) & (b ^ r)) < 0) throw DomainError("Rational: 64-bit integer overflow");
			return r;
		}

		void normalize()
		{
			if (den_ == 0) throw DivisionByZeroError("Rational: zero denominator");
			if (den_ < 0) { num_ = -num_; den_ = -den_; }
			long long g = NumberTheory::Gcd(num_, den_);
			if (g > 1) { num_ /= g; den_ /= g; }
		}

	public:
		Rational() = default;
		Rational(long long n) : num_(n), den_(1) {}
		Rational(long long n, long long d) : num_(n), den_(d) { normalize(); }

		long long Num() const { return num_; }
		long long Den() const { return den_; }
		Real ToReal() const { return static_cast<Real>(num_) / static_cast<Real>(den_); }

		Rational operator-() const { return Rational(-num_, den_); }

		Rational operator+(const Rational& r) const
		{
			long long g = NumberTheory::Gcd(den_, r.den_);
			long long den = mulChecked(den_ / g, r.den_);
			long long num = addChecked(mulChecked(num_, r.den_ / g), mulChecked(r.num_, den_ / g));
			return Rational(num, den);
		}

		Rational operator-(const Rational& r) const { return *this + (-r); }

		Rational operator*(const Rational& r) const
		{
			// Cross-reduce before multiplying to keep the operands small.
			long long g1 = NumberTheory::Gcd(num_, r.den_);
			long long g2 = NumberTheory::Gcd(r.num_, den_);
			long long num = mulChecked(num_ / g1, r.num_ / g2);
			long long den = mulChecked(den_ / g2, r.den_ / g1);
			return Rational(num, den);
		}

		Rational operator/(const Rational& r) const
		{
			if (r.num_ == 0) throw DivisionByZeroError("Rational: division by zero");
			return *this * Rational(r.den_, r.num_);
		}

		Rational& operator+=(const Rational& r) { *this = *this + r; return *this; }
		Rational& operator-=(const Rational& r) { *this = *this - r; return *this; }
		Rational& operator*=(const Rational& r) { *this = *this * r; return *this; }
		Rational& operator/=(const Rational& r) { *this = *this / r; return *this; }

		bool operator==(const Rational& r) const { return num_ == r.num_ && den_ == r.den_; }
		bool operator!=(const Rational& r) const { return !(*this == r); }

		bool operator<(const Rational& r) const
		{
#if defined(__SIZEOF_INT128__)
			return static_cast<__int128>(num_) * r.den_ < static_cast<__int128>(r.num_) * den_;
#else
			return mulChecked(num_, r.den_) < mulChecked(r.num_, den_);
#endif
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
