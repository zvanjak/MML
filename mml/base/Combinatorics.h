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
				return ExactInteger::Multiply(a, b);
			}

			inline long long BinomialLL(long long n, long long k)
			{
				k = std::min(k, n - k);
				long long result = 1;
				for (long long i = 1; i <= k; ++i) {
					long long num = n - k + i, divisor = i;
					const auto g = NumberTheory::Gcd(num, divisor);
					num /= g; divisor /= g;
					const auto g2 = NumberTheory::Gcd(result, divisor);
					result /= g2;
					result = MulCheckedLL(result, num);
				}
				return result;
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
			return detail::BinomialLL(n, k);
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
				total = ExactInteger::Add(total, c);
				result = detail::MulCheckedLL(result, detail::BinomialLL(total, c));
			}
			return result;
		}

		/// @brief Falling factorial x(x-1)...(x-k+1); throws on overflow.
		inline long long FallingFactorial(long long x, int k)
		{
			if (k < 0) throw DomainError("FallingFactorial: negative k");
			if (x >= 0 && x < k) return 0;
			long long r = 1;
			for (int i = 0; i < k; ++i)
				r = detail::MulCheckedLL(r, ExactInteger::Subtract(x, i));
			return r;
		}

		/// @brief Rising factorial x(x+1)...(x+k-1); throws on overflow.
		inline long long RisingFactorial(long long x, int k)
		{
			if (k < 0) throw DomainError("RisingFactorial: negative k");
			if (x <= 0 && x > -static_cast<long long>(k)) return 0;
			long long r = 1;
			for (int i = 0; i < k; ++i)
				r = detail::MulCheckedLL(r, ExactInteger::Add(x, i));
			return r;
		}

		// ---------------------------------------------------------------------------------
		//  named integer sequences
		// ---------------------------------------------------------------------------------

		/// @brief nth Catalan number C(2n,n)/(n+1).
		inline long long Catalan(int n)
		{
			if (n < 0) throw DomainError("Catalan: negative n");
			if (n > 35) throw DomainError("Catalan: result exceeds signed 64-bit range");
			long long result = 1;
			for (int i = 1; i <= n; ++i) {
				// C_i = C_(i-1) * 2(2i-1)/(i+1); cancel before multiplying.
				long long numerator = 2LL * (2LL * i - 1), denominator = i + 1;
				const auto g = NumberTheory::Gcd(numerator, denominator);
				numerator /= g; denominator /= g;
				result /= denominator;
				result = detail::MulCheckedLL(result, numerator);
			}
			return result;
		}

		/// @brief nth Fibonacci number (F(0)=0, F(1)=1); exact up to F(92).
		inline long long Fibonacci(int n)
		{
			if (n < 0) throw DomainError("Fibonacci: negative n");
			if (n > 92) throw DomainError("Fibonacci: result exceeds signed 64-bit range");
			if (n == 0) return 0;
			long long previous = 0, current = 1;
			for (int i = 1; i < n; ++i) {
				const auto next = ExactInteger::Add(previous, current);
				previous = current; current = next;
			}
			return current;
		}

		/// @brief nth Lucas number (L(0)=2, L(1)=1).
		inline long long Lucas(int n)
		{
			if (n < 0) throw DomainError("Lucas: negative n");
			if (n > 90) throw DomainError("Lucas: result exceeds signed 64-bit range");
			if (n == 0) return 2;
			long long previous = 2, current = 1;
			for (int i = 1; i < n; ++i) {
				const auto next = ExactInteger::Add(previous, current);
				previous = current; current = next;
			}
			return current;
		}

		/// @brief nth Bell number (number of set partitions of an n-element set).
		inline long long BellNumber(int n)
		{
			if (n < 0) throw DomainError("BellNumber: negative n");
			if (n > 25) throw DomainError("BellNumber: result exceeds signed 64-bit range");
			if (n == 0) return 1;
			std::vector<long long> row{ 1 };
			// The last entry of row i is B_(i+1); do not construct B_(n+1).
			for (int i = 1; i < n; ++i) {
				std::vector<long long> next(i + 1);
				next[0] = row.back();
				for (int j = 1; j <= i; ++j)
					next[j] = ExactInteger::Add(next[j - 1], row[j - 1]);
				row = std::move(next);
			}
			return row.back();
		}

		namespace detail
		{
			inline long long Stirling(int n, int k, bool firstKind)
			{
				if (n == 0) return k == 0 ? 1 : 0;
				if (k > n || k == 0) return 0;
				if (k == n || (!firstKind && k == 1)) return 1;
				if (k == n - 1) return BinomialLL(n, 2);
				// Track the smaller coordinate, cycles/blocks j or deficit i-j.
				// Only states which can still reach (n,k) are evaluated: overflowing
				// irrelevant triangle entries must not reject a representable answer.
				const bool byDeficit = n - k < k;
				const int width = byDeficit ? n - k : k;
				std::vector<long long> s{1};
				for (long long i = 1; i <= n; ++i) {
					const int lower = byDeficit ? 0 : static_cast<int>((std::max)(1LL, k - (n - i)));
					const int upper = static_cast<int>((std::min)(i, static_cast<long long>(width)));
					std::vector<long long> t(static_cast<std::size_t>(upper) + 1, 0);
					for (int j = lower; j <= upper; ++j) {
						const auto same = static_cast<std::size_t>(j) < s.size() ? s[j] : 0;
						if (byDeficit) {
							const auto coefficient = firstKind ? i - 1 : i - j;
							t[j] = ExactInteger::Add(same, j == 0 ? 0 :
								MulCheckedLL(coefficient, s[j - 1]));
						}
						else {
							const auto coefficient = firstKind ? i - 1 : j;
							t[j] = ExactInteger::Add(s[j - 1], MulCheckedLL(coefficient, same));
						}
					}
					s = std::move(t);
				}
				return s[width];
			}
		}

		/// @brief Unsigned Stirling number of the first kind [n k] (permutations of n with k cycles).
		inline long long StirlingFirst(int n, int k)
		{
			if (n < 0 || k < 0) throw DomainError("StirlingFirst: negative argument");
			return detail::Stirling(n, k, true);
		}

		/// @brief Stirling number of the second kind S(n,k) (partitions of n into k blocks).
		inline long long StirlingSecond(int n, int k)
		{
			if (n < 0 || k < 0) throw DomainError("StirlingSecond: negative argument");
			return detail::Stirling(n, k, false);
		}

		/// @brief nth Euler zigzag (up/down) number: alternating permutations of {1..n}.
		inline long long ZigzagNumber(int n)
		{
			if (n < 0) throw DomainError("ZigzagNumber: negative n");
			if (n > 23) throw DomainError("ZigzagNumber: result exceeds signed 64-bit range");
			std::vector<long long> e(n + 1, 0);
			e[0] = 1;                       // E(0,0)
			long long zig = 1;
			for (int i = 1; i <= n; ++i) {
				std::vector<long long> next(i + 1, 0);
				// next[k] = next[k-1] + e[i-k]  (boustrophedon / Entringer triangle)
				for (int k = 1; k <= i; ++k)
					next[k] = ExactInteger::Add(next[k - 1], e[i - k]);
				e = std::move(next);
				zig = e[i];
			}
			return zig;
		}

		/// @brief nth tangent (zag) number: coefficient in tan(x) = sum T_n x^(2n-1)/(2n-1)!.
		inline long long Tangent(int n)
		{
			if (n < 1) throw DomainError("Tangent: n must be >= 1");
			if (n > 12) throw DomainError("Tangent: result exceeds signed 64-bit range");
			return ZigzagNumber(2 * n - 1);
		}

		/// @brief nth Bernoulli number as an exact Rational (convention B_1 = +1/2).
		/// @details Uses the Akiyama-Tanigawa algorithm; throws DomainError once a numerator
		///          or denominator would exceed the 64-bit range (roughly n > 34).
		inline Rational Bernoulli(int n)
		{
			if (n < 0) throw DomainError("Bernoulli: negative n");
			if (n > 1 && n % 2 == 1) return Rational(0);
			// Grow only as far as the checked exact recurrence can reach.
			std::vector<Rational> a;
			for (long long m = 0; m <= n; ++m) {
				a.emplace_back(1, m + 1);
				for (long long j = m; j >= 1; --j)
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
				if (term < PrecisionValues<double>::Eps * zeta) break;   // series converged to double roundoff
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

		/// @brief Number of integer partitions p(n) via the positive generating-function recurrence.
		inline long long PartitionCount(int n)
		{
			if (n < 0) return 0;
			if (n > 405) throw DomainError("PartitionCount: result exceeds signed 64-bit range");
			std::vector<long long> p(n + 1, 0);
			p[0] = 1;
			// Positive generating-function recurrence avoids oversized alternating
			// partial sums in the pentagonal formula near p(405).
			for (int part = 1; part <= n; ++part) {
				for (int i = part; i <= n; ++i) {
					p[i] = ExactInteger::Add(p[i], p[i - part]);
				}
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
