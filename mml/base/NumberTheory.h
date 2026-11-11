///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        NumberTheory.h                                                      ///
///  Description: Runtime integer number theory (gcd, modular arithmetic, primes)     ///
///               64-bit routines: Miller-Rabin, Pollard rho, CRT, totient, sieve     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_NUMBER_THEORY_H
#define MML_NUMBER_THEORY_H

#include <algorithm>
#include <cstdint>
#include <limits>
#include <random>
#include <utility>
#include <vector>

#include <mml/MMLBase.h>
#include <mml/base/ExactInteger.h>

namespace MML
{
	/// @brief Runtime integer number theory on 64-bit integers.
	/// @details General-purpose routines operating on `long long` / `std::uint64_t`, distinct
	///          from the compile-time algebraic-structure module in `mml/base/Algebra`. All
	///          modular routines are correct across the full unsigned 64-bit range.
	namespace NumberTheory
	{
		// ---------------------------------------------------------------------------------
		//  gcd / lcm / extended Euclid
		// ---------------------------------------------------------------------------------

		/// @brief Greatest common divisor of |a| and |b| (Gcd(0,0)=0).
		/// @throws DomainError when the nonnegative result is 2^63.
		inline long long Gcd(long long a, long long b)
		{
			return ExactInteger::FromMagnitude(
				ExactInteger::Gcd(ExactInteger::Magnitude(a), ExactInteger::Magnitude(b)), false);
		}

		/// @brief Least common multiple of |a| and |b| (0 if either is 0).
		inline long long Lcm(long long a, long long b)
		{
			if (a == 0 || b == 0) return 0;
			auto ua = ExactInteger::Magnitude(a), ub = ExactInteger::Magnitude(b);
			ua /= ExactInteger::Gcd(ua, ub);
			const auto max = static_cast<std::uint64_t>((std::numeric_limits<long long>::max)());
			if (ua > max / ub) throw DomainError("Lcm: result exceeds signed 64-bit range");
			return static_cast<long long>(ua * ub);
		}

		/// @brief Extended Euclid: returns g = gcd(a,b) and sets x,y with a*x + b*y = g.
		inline long long ExtendedGcd(long long a, long long b, long long& x, long long& y)
		{
			auto old_r = ExactInteger::Magnitude(a), r = ExactInteger::Magnitude(b);
			std::uint64_t old_s = 1, s = 0, old_t = 0, t = 1;
			bool negative = false;
			while (r != 0) {
				const auto q = old_r / r, next_r = old_r % r;
				// The final unused coefficients can reach 2^63, so retain unsigned
				// magnitudes. Euclid bounds them by |b|/g and |a|/g respectively.
				const auto next_s = old_s + q * s, next_t = old_t + q * t;
				old_r = r; r = next_r;
				old_s = s; s = next_s;
				old_t = t; t = next_t;
				negative = !negative;
			}
			const auto g = ExactInteger::FromMagnitude(old_r, false);
			const auto new_x = ExactInteger::FromMagnitude(old_s, negative != (a < 0));
			const auto new_y = ExactInteger::FromMagnitude(old_t, (!negative) != (b < 0));
			x = new_x; y = new_y;
			return g;
		}

		namespace detail
		{
			inline std::uint64_t GcdU(std::uint64_t a, std::uint64_t b)
			{
				return ExactInteger::Gcd(a, b);
			}

			// a and b are residues in [0,m); the subtraction form cannot wrap.
			inline std::uint64_t AddMod(std::uint64_t a, std::uint64_t b, std::uint64_t m)
			{
				return a >= m - b ? a - (m - b) : a + b;
			}

			inline long long Residue(long long a, long long m)
			{
				if (m <= 0) throw DomainError("Modular arithmetic: modulus must be positive");
				const auto r = a % m;
				return r < 0 ? r + m : r;
			}
		}

		// ---------------------------------------------------------------------------------
		//  modular multiply / power / inverse (full 64-bit range)
		// ---------------------------------------------------------------------------------

		/// @brief (a * b) mod m, correct for the full unsigned 64-bit range.
		inline std::uint64_t MulMod(std::uint64_t a, std::uint64_t b, std::uint64_t m)
		{
			if (m == 0) throw DomainError("MulMod: modulus must be positive");
#if defined(__SIZEOF_INT128__)
			return static_cast<std::uint64_t>((static_cast<unsigned __int128>(a) * b) % m);
#else
			// Overflow-safe binary (double-and-add) multiply, valid for any m up to 2^64-1.
			std::uint64_t result = 0;
			a %= m;
			b %= m;
			while (b != 0) {
				if (b & 1u) result = detail::AddMod(result, a, m);
				a = detail::AddMod(a, a, m);
				b >>= 1;
			}
			return result;
#endif
		}

		/// @brief (base^exp) mod m by binary exponentiation.
		inline std::uint64_t ModPow(std::uint64_t base, std::uint64_t exp, std::uint64_t m)
		{
			if (m == 0) throw DomainError("ModPow: modulus must be positive");
			if (m == 1) return 0;
			std::uint64_t result = 1 % m;
			base %= m;
			while (exp != 0) {
				if (exp & 1u) result = MulMod(result, base, m);
				base = MulMod(base, base, m);
				exp >>= 1;
			}
			return result;
		}

		/// @brief Modular inverse of a modulo m; throws DomainError if gcd(a,m) != 1.
		inline long long ModInverse(long long a, long long m)
		{
			if (m <= 0) throw DomainError("ModInverse: modulus must be positive");
			long long aa = detail::Residue(a, m);
			long long x, y;
			if (ExtendedGcd(aa, m, x, y) != 1)
				throw DomainError("ModInverse: arguments are not coprime");
			return detail::Residue(x, m);
		}

		// ---------------------------------------------------------------------------------
		//  primality (deterministic Miller-Rabin) and factorization (Pollard rho)
		// ---------------------------------------------------------------------------------

		/// @brief Deterministic Miller-Rabin primality test, exact for all 64-bit integers.
		inline bool IsPrime(std::uint64_t n)
		{
			if (n < 2) return false;
			for (std::uint64_t p : { 2ull, 3ull, 5ull, 7ull, 11ull, 13ull, 17ull, 19ull,
									 23ull, 29ull, 31ull, 37ull }) {
				if (n % p == 0) return n == p;
			}
			std::uint64_t d = n - 1;
			int s = 0;
			while ((d & 1u) == 0) { d >>= 1; ++s; }
			// These seven bases give a deterministic test over the whole 64-bit range.
			for (std::uint64_t a : { 2ull, 325ull, 9375ull, 28178ull, 450775ull,
									 9780504ull, 1795265022ull }) {
				a %= n;
				if (a == 0) continue;
				std::uint64_t x = ModPow(a, d, n);
				if (x == 1 || x == n - 1) continue;
				bool composite = true;
				for (int r = 1; r < s; ++r) {
					x = MulMod(x, x, n);
					if (x == n - 1) { composite = false; break; }
				}
				if (composite) return false;
			}
			return true;
		}

		namespace detail
		{
			/// @brief One nontrivial factor of composite n via Pollard's rho (Floyd variant).
			inline std::uint64_t PollardRho(std::uint64_t n)
			{
				if (n % 2 == 0) return 2;
				if (n % 3 == 0) return 3;
				std::mt19937_64 rng(0x9E3779B97F4A7C15ull ^ n);
				while (true) {
					std::uint64_t c = rng() % (n - 1) + 1;
					std::uint64_t x = rng() % n, y = x, d = 1;
					auto f = [&](std::uint64_t v) { return AddMod(MulMod(v, v, n), c, n); };
					while (d == 1) {
						x = f(x);
						y = f(f(y));
						std::uint64_t diff = (x > y) ? (x - y) : (y - x);
						if (diff == 0) { d = n; break; }
						d = GcdU(diff, n);
					}
					if (d != n && d != 0) return d;
				}
			}

			inline void FactorizeInto(std::uint64_t n, std::vector<std::uint64_t>& out)
			{
				if (n < 2) return;
				if (IsPrime(n)) { out.push_back(n); return; }
				std::uint64_t d = PollardRho(n);
				FactorizeInto(d, out);
				FactorizeInto(n / d, out);
			}
		}

		/// @brief Sorted list of prime factors of n with multiplicity (empty for n<2).
		inline std::vector<std::uint64_t> Factorize(std::uint64_t n)
		{
			std::vector<std::uint64_t> out;
			detail::FactorizeInto(n, out);
			std::sort(out.begin(), out.end());
			return out;
		}

		/// @brief Prime factorization of n as (prime, exponent) pairs, ascending by prime.
		inline std::vector<std::pair<std::uint64_t, int>> FactorizePowers(std::uint64_t n)
		{
			std::vector<std::pair<std::uint64_t, int>> out;
			for (std::uint64_t p : Factorize(n)) {
				if (!out.empty() && out.back().first == p) out.back().second++;
				else out.push_back({ p, 1 });
			}
			return out;
		}

		// ---------------------------------------------------------------------------------
		//  multiplicative functions
		// ---------------------------------------------------------------------------------

		/// @brief Euler's totient phi(n): count of integers in [1,n] coprime to n.
		inline long long EulerTotient(long long n)
		{
			if (n <= 0) throw DomainError("EulerTotient: argument must be positive");
			if (n == 1) return 1;
			long long result = n;
			for (auto [p, e] : FactorizePowers(static_cast<std::uint64_t>(n))) {
				(void)e;
				result -= result / static_cast<long long>(p);
			}
			return result;
		}

		/// @brief Mobius function mu(n): 0 if n has a squared prime factor, else (-1)^k.
		inline int MobiusMu(long long n)
		{
			if (n <= 0) throw DomainError("MobiusMu: argument must be positive");
			if (n == 1) return 1;
			int k = 0;
			for (auto [p, e] : FactorizePowers(static_cast<std::uint64_t>(n))) {
				(void)p;
				if (e > 1) return 0;
				++k;
			}
			return (k & 1) ? -1 : 1;
		}

		/// @brief Number of positive divisors of n (tau).
		inline long long DivisorCount(long long n)
		{
			if (n <= 0) throw DomainError("DivisorCount: argument must be positive");
			long long result = 1;
			for (auto [p, e] : FactorizePowers(static_cast<std::uint64_t>(n))) {
				(void)p;
				result *= (e + 1);
			}
			return result;
		}

		/// @brief Sum of the positive divisors of n (sigma).
		inline long long DivisorSum(long long n)
		{
			if (n <= 0) throw DomainError("DivisorSum: argument must be positive");
			long long result = 1;
			for (auto [p, e] : FactorizePowers(static_cast<std::uint64_t>(n))) {
				long long term = 1, power = 1;
				for (int i = 0; i < e; ++i) {
					power = ExactInteger::Multiply(power, static_cast<long long>(p));
					term = ExactInteger::Add(term, power);
				}
				result = ExactInteger::Multiply(result, term);
			}
			return result;
		}

		// ---------------------------------------------------------------------------------
		//  prime enumeration
		// ---------------------------------------------------------------------------------

		/// @brief All primes <= limit via the sieve of Eratosthenes.
		inline std::vector<std::uint64_t> PrimesUpTo(std::uint64_t limit)
		{
			std::vector<std::uint64_t> primes;
			if (limit < 2) return primes;
			if (limit >= std::vector<bool>().max_size() ||
				limit >= (std::numeric_limits<std::size_t>::max)())
				throw DomainError("PrimesUpTo: sieve size is not representable");
			std::vector<bool> composite(static_cast<std::size_t>(limit) + 1, false);
			for (std::uint64_t i = 2; i <= limit / i; ++i) {
				if (composite[static_cast<std::size_t>(i)]) continue;
				for (std::uint64_t j = i * i; j <= limit; ) {
					composite[static_cast<std::size_t>(j)] = true;
					if (limit - j < i) break;
					j += i;
				}
			}
			for (std::uint64_t i = 2; i <= limit; ++i)
				if (!composite[static_cast<std::size_t>(i)]) primes.push_back(i);
			return primes;
		}

		/// @brief Smallest prime strictly greater than n; DomainError if none fits uint64_t.
		inline std::uint64_t NextPrime(std::uint64_t n)
		{
			constexpr auto max = (std::numeric_limits<std::uint64_t>::max)();
			if (n == max) throw DomainError("NextPrime: unsigned 64-bit range exhausted");
			std::uint64_t candidate = n + 1;
			while (!IsPrime(candidate)) {
				if (candidate == max) throw DomainError("NextPrime: unsigned 64-bit range exhausted");
				++candidate;
			}
			return candidate;
		}

		// ---------------------------------------------------------------------------------
		//  Chinese remainder theorem (non-coprime moduli supported)
		// ---------------------------------------------------------------------------------

		namespace detail
		{
			// (a * b) mod m for signed operands, with 0 <= result < m and m > 0.
			inline long long MulModS(long long a, long long b, long long m)
			{
				long long ra = Residue(a, m);
				long long rb = Residue(b, m);
				return static_cast<long long>(MulMod(static_cast<std::uint64_t>(ra),
					static_cast<std::uint64_t>(rb), static_cast<std::uint64_t>(m)));
			}
		}

		/// @brief Solve x ≡ rem[i] (mod mod[i]) for all i.
		/// @return Pair (x, L) with x the least non-negative solution modulo L = lcm(mod[i]).
		/// @throws DomainError on invalid input, inconsistent congruences, or L exceeding LLONG_MAX.
		inline std::pair<long long, long long> ChineseRemainder(
			const std::vector<long long>& rem, const std::vector<long long>& mod)
		{
			if (rem.size() != mod.size() || rem.empty())
				throw DomainError("ChineseRemainder: rem and mod must be non-empty and equal length");

			long long r0 = 0, m0 = 1;
			for (std::size_t i = 0; i < mod.size(); ++i) {
				if (mod[i] <= 0) throw DomainError("ChineseRemainder: moduli must be positive");
				long long ri = detail::Residue(rem[i], mod[i]);
				long long mi = mod[i];

				long long p, q;
				long long g = ExtendedGcd(m0, mi, p, q);
				long long diff = ri - r0;
				if (diff % g != 0)
					throw DomainError("ChineseRemainder: congruences are inconsistent");

				long long lcm = ExactInteger::Multiply(m0 / g, mi);
				long long modReduced = mi / g;
				long long t = detail::Residue(diff / g, modReduced);
				long long mult = detail::MulModS(t, p, modReduced);
				long long add = detail::MulModS(m0 % lcm, mult, lcm);
				long long newr = static_cast<long long>(detail::AddMod(
					static_cast<std::uint64_t>(r0), static_cast<std::uint64_t>(add),
					static_cast<std::uint64_t>(lcm)));

				r0 = newr;
				m0 = lcm;
			}
			return { r0, m0 };
		}
	} // namespace NumberTheory
} // namespace MML

#endif // MML_NUMBER_THEORY_H
