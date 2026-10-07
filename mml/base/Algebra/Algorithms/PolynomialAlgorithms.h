///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        PolynomialAlgorithms.h                                              ///
///  Description: Exact polynomial GCD, extended GCD, and irreducibility checks       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_POLYNOMIAL_ALGORITHMS_H
#define MML_POLYNOMIAL_ALGORITHMS_H

#include <mml/MMLExceptions.h>
#include <mml/base/Algebra/Polynomial.h>

#include <cstdint>
#include <vector>

namespace MML
{
	namespace Algebra
	{
		template<class Coeff>
		struct PolynomialExtendedGcdResult
		{
			Polynomial<Coeff> gcd;
			Polynomial<Coeff> left_coefficient;
			Polynomial<Coeff> right_coefficient;
		};

		template<class Coeff>
		Polynomial<Coeff> PolynomialGcd(Polynomial<Coeff> left, Polynomial<Coeff> right)
		{
			while (!right.is_zero()) {
				auto remainder = left % right;
				left = std::move(right);
				right = std::move(remainder);
			}
			return left.monic();
		}

		template<class Coeff>
		PolynomialExtendedGcdResult<Coeff> PolynomialExtendedGcd(
			Polynomial<Coeff> left, Polynomial<Coeff> right)
		{
			Polynomial<Coeff> oldLeft = Polynomial<Coeff>::one();
			Polynomial<Coeff> currentLeft = Polynomial<Coeff>::zero();
			Polynomial<Coeff> oldRight = Polynomial<Coeff>::zero();
			Polynomial<Coeff> currentRight = Polynomial<Coeff>::one();

			while (!right.is_zero()) {
				auto division = left.divmod(right);
				left = std::move(right);
				right = std::move(division.second);

				auto nextLeft = oldLeft - division.first * currentLeft;
				oldLeft = std::move(currentLeft);
				currentLeft = std::move(nextLeft);

				auto nextRight = oldRight - division.first * currentRight;
				oldRight = std::move(currentRight);
				currentRight = std::move(nextRight);
			}

			if (left.is_zero())
				return {left, oldLeft, oldRight};
			const Coeff scale = left.leading_coefficient().inverse();
			return {left * scale, oldLeft * scale, oldRight * scale};
		}

		template<class Coeff>
		Polynomial<Coeff> PolynomialPowMod(Polynomial<Coeff> base, std::uint64_t exponent,
			const Polynomial<Coeff>& modulus)
		{
			Polynomial<Coeff> result = Polynomial<Coeff>::one();
			base = base % modulus;
			while (exponent > 0) {
				if ((exponent & 1U) != 0)
					result = (result * base) % modulus;
				base = (base * base) % modulus;
				exponent >>= 1U;
			}
			return result;
		}

		template<class Coeff>
		bool IsIrreducible(const Polynomial<Coeff>& polynomial)
		{
			static_assert(Coeff::is_field && IsPrime<Coeff::modulus>::value,
				"Irreducibility test requires a prime field coefficient type");
			constexpr int Prime = Coeff::modulus;
			if (polynomial.degree() <= 0)
				return false;
			const auto monic = polynomial.monic();
			const auto x = Polynomial<Coeff>::monomial(1);
			Polynomial<Coeff> power = x;

			for (int step = 1; step <= monic.degree() / 2; step++) {
				power = PolynomialPowMod(power, Prime, monic);
				if (PolynomialGcd(monic, power - x).degree() > 0)
					return false;
			}

			power = x;
			for (int step = 0; step < monic.degree(); step++)
				power = PolynomialPowMod(power, Prime, monic);
			return power == x;
		}
	}
}

#endif // MML_POLYNOMIAL_ALGORITHMS_H