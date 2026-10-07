///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Combinatorics.h                                                     ///
///  Description: Combinatorial counting, integer sequences, partitions, enumeration  ///
///               Binomials, Catalan/Fibonacci/Bell/Stirling, Bernoulli (exact/real)  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_COMBINATORICS_H
#define MML_COMBINATORICS_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

#include <mml/MMLBase.h>
#include <mml/base/NumberTheory.h>
#include <mml/base/Rational.h>

namespace MML
{
	/// @brief Combinatorial counting, integer sequences, partitions, and basic enumeration.
	/// @details Exact integer results use `long long` and throw DomainError once a value would
	///          exceed the 64-bit range; floating alternatives (e.g. BinomialCoefficientReal)
	///          extend the usable range. Enumeration helpers yield index vectors over {0..n-1}.
	namespace Combinatorics
	{
		namespace detail
		{
			inline long long MulCheckedLL(long long a, long long b)
			{
				if (a == 0 || b == 0) return 0;
				long long r = a * b;
				if (r / a != b) throw DomainError("Combinatorics: 64-bit integer overflow");
				return r;
			}
		}

		// ---------------------------------------------------------------------------------
		//  counting
		// ---------------------------------------------------------------------------------

		/// @brief Exact binomial coefficient C(n,k); throws DomainError if it exceeds 64 bits.
		inline long long BinomialCoefficient(int n, int k)
		{
			if (n < 0) throw DomainError("BinomialCoefficient: negative n");
			if (k < 0 || k > n) return 0;
			k = std::min(k, n - k);
			long long result = 1;
			for (int i = 1; i <= k; ++i) {
				long long num = n - k + i;
				long long ii = i;
				// Reduce num/ii against ii and the running result so the product stays exact.
				long long g = NumberTheory::Gcd(num, ii);
				num /= g; ii /= g;
				long long g2 = NumberTheory::Gcd(result, ii);
				result /= g2; ii /= g2;   // ii is now 1 (C(n,k) is integral)
				if (num != 0 && result > (std::numeric_limits<long long>::max)() / num)
					throw DomainError("BinomialCoefficient: exceeds 64-bit range; use BinomialCoefficientReal");
				result *= num;
			}
			return result;
		}

		/// @brief Binomial coefficient as a double via log-gamma (usable for large n).
		inline double BinomialCoefficientReal(int n, int k)
		{
			if (n < 0) throw DomainError("BinomialCoefficientReal: negative n");
			if (k < 0 || k > n) return 0.0;
			double v = std::lgamma((double)n + 1.0) - std::lgamma((double)k + 1.0)
				- std::lgamma((double)(n - k) + 1.0);
			return std::round(std::exp(v));
		}

		/// @brief Multinomial coefficient (sum counts)! / prod(counts!); throws on overflow.
		inline long long Multinomial(const std::vector<int>& counts)
		{
			long long total = 0, result = 1;
			for (int c : counts) {
				if (c < 0) throw DomainError("Multinomial: negative count");
				total += c;
				result = detail::MulCheckedLL(result, BinomialCoefficient((int)total, c));
			}
			return result;
		}

		/// @brief Falling factorial x(x-1)...(x-k+1); throws on overflow.
		inline long long FallingFactorial(long long x, int k)
		{
			if (k < 0) throw DomainError("FallingFactorial: negative k");
			long long r = 1;
			for (int i = 0; i < k; ++i) r = detail::MulCheckedLL(r, x - i);
			return r;
		}

		/// @brief Rising factorial x(x+1)...(x+k-1); throws on overflow.
		inline long long RisingFactorial(long long x, int k)
		{
			if (k < 0) throw DomainError("RisingFactorial: negative k");
			long long r = 1;
			for (int i = 0; i < k; ++i) r = detail::MulCheckedLL(r, x + i);
			return r;
		}

		// ---------------------------------------------------------------------------------
		//  named integer sequences
		// ---------------------------------------------------------------------------------

		/// @brief nth Catalan number C(2n,n)/(n+1).
		inline long long Catalan(int n)
		{
			if (n < 0) throw DomainError("Catalan: negative n");
			return BinomialCoefficient(2 * n, n) / (n + 1);
		}

		namespace detail
		{
			// (F(n), F(n+1)) via fast doubling.
			inline std::pair<long long, long long> FibPair(int n)
			{
				if (n == 0) return { 0, 1 };
				auto hp = FibPair(n / 2);
				long long a = hp.first, b = hp.second;
				long long c = a * (2 * b - a);
				long long d = a * a + b * b;
				if (n & 1) return { d, c + d };
				return { c, d };
			}
		}

		/// @brief nth Fibonacci number (F(0)=0, F(1)=1); exact up to F(92).
		inline long long Fibonacci(int n)
		{
			if (n < 0) throw DomainError("Fibonacci: negative n");
			return detail::FibPair(n).first;
		}

		/// @brief nth Lucas number (L(0)=2, L(1)=1).
		inline long long Lucas(int n)
		{
			if (n < 0) throw DomainError("Lucas: negative n");
			auto fp = detail::FibPair(n);
			return 2 * fp.second - fp.first;
		}

		/// @brief nth Bell number (number of set partitions of an n-element set).
		inline long long BellNumber(int n)
		{
			if (n < 0) throw DomainError("BellNumber: negative n");
			std::vector<long long> row{ 1 };
			for (int i = 1; i <= n; ++i) {
				std::vector<long long> next(i + 1);
				next[0] = row.back();
				for (int j = 1; j <= i; ++j) next[j] = next[j - 1] + row[j - 1];
				row = std::move(next);
			}
			return row[0];
		}

		/// @brief Unsigned Stirling number of the first kind [n k] (permutations of n with k cycles).
		inline long long StirlingFirst(int n, int k)
		{
			if (n < 0 || k < 0) throw DomainError("StirlingFirst: negative argument");
			if (n == 0) return k == 0 ? 1 : 0;
			if (k > n) return 0;
			std::vector<long long> s(k + 1, 0);
			s[0] = 1;   // c(0,0) = 1
			for (int i = 1; i <= n; ++i) {
				std::vector<long long> t(k + 1, 0);
				for (int j = 1; j <= std::min(i, k); ++j)
					t[j] = s[j - 1] + (long long)(i - 1) * s[j];
				s = std::move(t);
			}
			return s[k];
		}

		/// @brief Stirling number of the second kind S(n,k) (partitions of n into k blocks).
		inline long long StirlingSecond(int n, int k)
		{
			if (n < 0 || k < 0) throw DomainError("StirlingSecond: negative argument");
			if (n == 0) return k == 0 ? 1 : 0;
			if (k > n) return 0;
			std::vector<long long> s(k + 1, 0);
			s[0] = 1;   // S(0,0) = 1
			for (int i = 1; i <= n; ++i) {
				std::vector<long long> t(k + 1, 0);
				for (int j = 1; j <= std::min(i, k); ++j)
					t[j] = (long long)j * s[j] + s[j - 1];
				s = std::move(t);
			}
			return s[k];
		}

		/// @brief nth Euler zigzag (up/down) number: alternating permutations of {1..n}.
		inline long long ZigzagNumber(int n)
		{
			if (n < 0) throw DomainError("ZigzagNumber: negative n");
			std::vector<long long> e(n + 1, 0);
			e[0] = 1;                       // E(0,0)
			long long zig = 1;
			for (int i = 1; i <= n; ++i) {
				std::vector<long long> next(i + 1, 0);
				// next[k] = next[k-1] + e[i-k]  (boustrophedon / Entringer triangle)
				for (int k = 1; k <= i; ++k)
					next[k] = next[k - 1] + e[i - k];
				e = std::move(next);
				zig = e[i];
			}
			return zig;
		}

		/// @brief nth tangent (zag) number: coefficient in tan(x) = sum T_n x^(2n-1)/(2n-1)!.
		inline long long Tangent(int n)
		{
			if (n < 1) throw DomainError("Tangent: n must be >= 1");
			return ZigzagNumber(2 * n - 1);
		}

		/// @brief nth Bernoulli number as an exact Rational (convention B_1 = +1/2).
		/// @details Uses the Akiyama-Tanigawa algorithm; throws DomainError once a numerator
		///          or denominator would exceed the 64-bit range (roughly n > 34).
		inline Rational Bernoulli(int n)
		{
			if (n < 0) throw DomainError("Bernoulli: negative n");
			std::vector<Rational> a(n + 1);
			for (int m = 0; m <= n; ++m) {
				a[m] = Rational(1, m + 1);
				for (int j = m; j >= 1; --j)
					a[j - 1] = Rational(j) * (a[j - 1] - a[j]);
			}
			return a[0];
		}

		/// @brief nth Bernoulli number as a double (convention B_1 = +1/2); wide range via zeta.
		inline double BernoulliReal(int n)
		{
			if (n < 0) throw DomainError("BernoulliReal: negative n");
			if (n == 1) return 0.5;
			if (n % 2 == 1) return 0.0;
			// Prefer the exact rational value while it fits 64 bits (fast and precise).
			try {
				const Rational value = Bernoulli(n);
				return static_cast<double>(value.Num()) / static_cast<double>(value.Den());
			}
			catch (const DomainError&) { /* fall through to the asymptotic zeta formula */ }
			// B_2m = (-1)^(m+1) * 2 * (2m)! * zeta(2m) / (2*pi)^(2m)
			double zeta = 0.0, term = 0.0;
			for (int j = 1; j < 1000000; ++j) {
				term = std::pow((double)j, -(double)n);
				zeta += term;
				if (term < 1e-18 * zeta) break;
			}
			const double pi = std::acos(-1.0);
			double logMag = std::log(2.0) + std::lgamma((double)n + 1.0) + std::log(zeta)
				- (double)n * std::log(2.0 * pi);
			int m = n / 2;
			double sign = (m % 2 == 1) ? 1.0 : -1.0;
			return sign * std::exp(logMag);
		}

		// ---------------------------------------------------------------------------------
		//  integer partitions
		// ---------------------------------------------------------------------------------

		/// @brief Number of integer partitions p(n) via Euler's pentagonal number theorem.
		inline long long PartitionCount(int n)
		{
			if (n < 0) return 0;
			std::vector<long long> p(n + 1, 0);
			p[0] = 1;
			for (int i = 1; i <= n; ++i) {
				long long sum = 0;
				for (int k = 1; ; ++k) {
					long long g1 = (long long)k * (3 * k - 1) / 2;   // generalized pentagonal numbers
					long long g2 = (long long)k * (3 * k + 1) / 2;
					if (g1 > i && g2 > i) break;
					long long sign = (k & 1) ? 1 : -1;
					if (g1 <= i) sum += sign * p[i - g1];
					if (g2 <= i) sum += sign * p[i - g2];
				}
				p[i] = sum;
			}
			return p[n];
		}

		namespace detail
		{
			template<class Fn>
			void PartitionRec(int remaining, int maxPart, std::vector<int>& current, Fn& fn)
			{
				if (remaining == 0) { fn(current); return; }
				for (int part = std::min(remaining, maxPart); part >= 1; --part) {
					current.push_back(part);
					PartitionRec(remaining - part, part, current, fn);
					current.pop_back();
				}
			}
		}

		/// @brief Enumerate every partition of n (as non-increasing part vectors), calling fn on each.
		template<class Fn>
		void ForEachPartition(int n, Fn&& fn)
		{
			if (n < 0) throw DomainError("ForEachPartition: negative n");
			std::vector<int> current;
			if (n == 0) { fn(current); return; }
			detail::PartitionRec(n, n, current, fn);
		}

		// ---------------------------------------------------------------------------------
		//  basic enumeration over {0, 1, ..., n-1}
		// ---------------------------------------------------------------------------------

		/// @brief Advance c to the next k-combination of {0..n-1} in lexicographic order.
		/// @return false when c already holds the last combination.
		inline bool NextCombination(std::vector<int>& c, int n)
		{
			int k = (int)c.size();
			if (k == 0) return false;
			int i = k - 1;
			while (i >= 0 && c[i] == n - k + i) --i;
			if (i < 0) return false;
			++c[i];
			for (int j = i + 1; j < k; ++j) c[j] = c[j - 1] + 1;
			return true;
		}

		/// @brief Enumerate every k-subset of {0..n-1} in lexicographic order, calling fn on each.
		template<class Fn>
		void ForEachCombination(int n, int k, Fn&& fn)
		{
			if (k < 0 || k > n) return;
			std::vector<int> c(k);
			for (int i = 0; i < k; ++i) c[i] = i;
			do { fn(c); } while (NextCombination(c, n));
		}

		/// @brief Advance c to the next length-k multiset (combination with repetition) of {0..n-1}.
		/// @return false when c already holds the last multiset.
		inline bool NextMultiset(std::vector<int>& c, int n)
		{
			int k = (int)c.size();
			if (k == 0) return false;
			int i = k - 1;
			while (i >= 0 && c[i] == n - 1) --i;
			if (i < 0) return false;
			int v = c[i] + 1;
			for (int j = i; j < k; ++j) c[j] = v;
			return true;
		}

		/// @brief Enumerate every length-k multiset over {0..n-1}, calling fn on each.
		template<class Fn>
		void ForEachMultiset(int n, int k, Fn&& fn)
		{
			if (k < 0 || (n <= 0 && k > 0)) return;
			std::vector<int> c(k, 0);
			do { fn(c); } while (NextMultiset(c, n));
		}

		/// @brief Enumerate every permutation of {0..n-1} in lexicographic order, calling fn on each.
		template<class Fn>
		void ForEachPermutation(int n, Fn&& fn)
		{
			if (n < 0) throw DomainError("ForEachPermutation: negative n");
			std::vector<int> p(n);
			for (int i = 0; i < n; ++i) p[i] = i;
			do { fn(p); } while (std::next_permutation(p.begin(), p.end()));
		}
	} // namespace Combinatorics
} // namespace MML

#endif // MML_COMBINATORICS_H
