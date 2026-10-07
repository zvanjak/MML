///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        CyclicGroup.h                                                       ///
///  Description: Finite cyclic groups and polygon rotation actions                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CYCLIC_GROUP_H
#define MML_CYCLIC_GROUP_H

#include <mml/MMLExceptions.h>
#include <mml/base/Algebra/Permutation.h>

#include <cstddef>
#include <utility>
#include <vector>

namespace MML
{
	namespace Algebra
	{
		struct CyclicElement
		{
			int exponent = 0;

			bool operator==(const CyclicElement& other) const noexcept { return exponent == other.exponent; }
			bool operator!=(const CyclicElement& other) const noexcept { return !(*this == other); }
		};

		/// @brief The cyclic group C_n represented by powers of one generator.
		class CyclicGroup
		{
			int _degree;
			std::vector<CyclicElement> _elements;

			int normalize(int exponent) const noexcept
			{
				int result = exponent % _degree;
				return result < 0 ? result + _degree : result;
			}

		public:
			using element_type = CyclicElement;

			explicit CyclicGroup(int degree)
				: _degree(degree)
			{
				if (degree <= 0)
					throw ArgumentError("CyclicGroup degree must be positive");
				_elements.reserve(static_cast<std::size_t>(degree));
				for (int exponent = 0; exponent < degree; exponent++)
					_elements.push_back({exponent});
			}

			const std::vector<CyclicElement>& elements() const noexcept { return _elements; }
			CyclicElement identity() const noexcept { return {0}; }
			CyclicElement generator() const noexcept { return rotation(1); }
			CyclicElement rotation(int exponent) const noexcept { return {normalize(exponent)}; }
			std::size_t order() const noexcept { return _elements.size(); }
			int degree() const noexcept { return _degree; }

			CyclicElement compose(CyclicElement left, CyclicElement right) const noexcept
			{
				return rotation(left.exponent + right.exponent);
			}

			CyclicElement inverse(CyclicElement element) const noexcept
			{
				return rotation(-element.exponent);
			}

			bool contains(CyclicElement element) const noexcept
			{
				return element.exponent >= 0 && element.exponent < _degree;
			}

			DynamicPermutation permutation(CyclicElement element) const
			{
				if (!contains(element))
					throw ArgumentError("CyclicGroup element is not canonical for this group");
				std::vector<int> images(static_cast<std::size_t>(_degree));
				for (int vertex = 0; vertex < _degree; vertex++)
					images[static_cast<std::size_t>(vertex)] = normalize(vertex + element.exponent);
				return DynamicPermutation(std::move(images));
			}
		};
	}
}

#endif // MML_CYCLIC_GROUP_H