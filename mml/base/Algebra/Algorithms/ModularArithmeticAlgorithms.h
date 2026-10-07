///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ModularArithmeticAlgorithms.h                                       ///
///  Description: Number-theory helpers for prime modular fields                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_MODULAR_ARITHMETIC_ALGORITHMS_H
#define MML_MODULAR_ARITHMETIC_ALGORITHMS_H

#include <mml/MMLExceptions.h>
#include <mml/base/Algebra/ModInt.h>

#include <cstdint>
#include <vector>

namespace MML
{
	namespace Algebra
	{
		namespace Detail
		{
			inline std::vector<int> DistinctPrimeFactors(int value)
			{
				std::vector<int> factors;
				for (int divisor = 2; static_cast<std::int64_t>(divisor) * divisor <= value; divisor++) {
					if (value % divisor != 0)
						continue;
					factors.push_back(divisor);
					while (value % divisor == 0)
						value /= divisor;
				}
				if (value > 1)
					factors.push_back(value);
				return factors;
			}
		}

		template<int Prime>
		PrimeFieldElement<Prime> PrimitiveRoot()
		{
			static_assert(IsPrime<Prime>::value, "PrimitiveRoot requires a prime modulus");
			if constexpr (Prime == 2)
				return PrimeFieldElement<Prime>(1);

			const auto factors = Detail::DistinctPrimeFactors(Prime - 1);
			for (int candidate = 2; candidate < Prime; candidate++) {
				const PrimeFieldElement<Prime> element(candidate);
				bool primitive = true;
				for (int factor : factors)
					if (element.pow(static_cast<std::uint64_t>((Prime - 1) / factor)) ==
						PrimeFieldElement<Prime>(1)) {
						primitive = false;
						break;
					}
				if (primitive)
					return element;
			}
			throw DomainError("Prime field has no primitive root");
		}

		template<int Prime>
		int LegendreSymbol(PrimeFieldElement<Prime> value)
		{
			static_assert(IsPrime<Prime>::value, "LegendreSymbol requires a prime modulus");
			if constexpr (Prime == 2)
				return value == PrimeFieldElement<Prime>(0) ? 0 : 1;
			if (value == PrimeFieldElement<Prime>(0))
				return 0;
			const auto result = value.pow(static_cast<std::uint64_t>((Prime - 1) / 2));
			return result == PrimeFieldElement<Prime>(1) ? 1 : -1;
		}

		template<int Prime>
		bool IsQuadraticResidue(PrimeFieldElement<Prime> value)
		{
			return LegendreSymbol<Prime>(value) >= 0;
		}
	}
}

#endif // MML_MODULAR_ARITHMETIC_ALGORITHMS_H