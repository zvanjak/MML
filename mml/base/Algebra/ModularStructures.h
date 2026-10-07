///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ModularStructures.h                                                 ///
///  Description: Finite ring and prime-field structure adapters                     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_MODULAR_STRUCTURES_H
#define MML_MODULAR_STRUCTURES_H

#include <mml/base/Algebra/ModInt.h>

#include <vector>

namespace MML
{
	namespace Algebra
	{
		template<int Modulus>
		class ModularRing
		{
			std::vector<ModInt<Modulus>> _elements;

		public:
			using element_type = ModInt<Modulus>;
			using scalar_type = element_type;

			ModularRing()
			{
				_elements.reserve(Modulus);
				for (int value = 0; value < Modulus; value++)
					_elements.emplace_back(value);
			}

			const std::vector<element_type>& elements() const noexcept { return _elements; }
			constexpr element_type zero() const noexcept { return element_type(0); }
			constexpr element_type one() const noexcept { return element_type(1); }
			constexpr element_type add(element_type left, element_type right) const noexcept { return left + right; }
			constexpr element_type negate(element_type element) const noexcept { return -element; }
			constexpr element_type multiply(element_type left, element_type right) const noexcept { return left * right; }
			std::size_t order() const noexcept { return _elements.size(); }
		};

		template<int Prime>
		class PrimeField : public ModularRing<Prime>
		{
			static_assert(IsPrime<Prime>::value, "PrimeField modulus must be prime");

		public:
			using element_type = PrimeFieldElement<Prime>;
			using scalar_type = element_type;

			element_type inverse(element_type element) const { return element.inverse(); }
		};
	}
}

#endif // MML_MODULAR_STRUCTURES_H