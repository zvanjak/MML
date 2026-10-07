///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DihedralGroup.h                                                     ///
///  Description: Finite dihedral groups and regular-polygon symmetry actions         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DIHEDRAL_GROUP_H
#define MML_DIHEDRAL_GROUP_H

#include <mml/MMLExceptions.h>
#include <mml/base/Algebra/Permutation.h>

#include <cstddef>
#include <utility>
#include <vector>

namespace MML
{
	namespace Algebra
	{
		/// @brief Canonical dihedral element r^rotation s^reflected.
		struct DihedralElement
		{
			int rotation = 0;
			bool reflected = false;

			bool operator==(const DihedralElement& other) const noexcept
			{
				return rotation == other.rotation && reflected == other.reflected;
			}
			bool operator!=(const DihedralElement& other) const noexcept { return !(*this == other); }
		};

		/// @brief The order-2n symmetry group D_n of a regular n-gon.
		class DihedralGroup
		{
			int _degree;
			std::vector<DihedralElement> _elements;

			int normalize(int exponent) const noexcept
			{
				int result = exponent % _degree;
				return result < 0 ? result + _degree : result;
			}

		public:
			using element_type = DihedralElement;

			explicit DihedralGroup(int degree)
				: _degree(degree)
			{
				if (degree < 2)
					throw ArgumentError("DihedralGroup degree must be at least two");
				_elements.reserve(static_cast<std::size_t>(2 * degree));
				for (int exponent = 0; exponent < degree; exponent++)
					_elements.push_back({exponent, false});
				for (int exponent = 0; exponent < degree; exponent++)
					_elements.push_back({exponent, true});
			}

			const std::vector<DihedralElement>& elements() const noexcept { return _elements; }
			DihedralElement identity() const noexcept { return {0, false}; }
			DihedralElement generator() const noexcept { return rotation(1); }
			DihedralElement reflection() const noexcept { return {0, true}; }
			DihedralElement rotation(int exponent) const noexcept { return {normalize(exponent), false}; }
			DihedralElement reflection(int exponent) const noexcept { return {normalize(exponent), true}; }
			std::size_t order() const noexcept { return _elements.size(); }
			int degree() const noexcept { return _degree; }

			DihedralElement compose(DihedralElement left, DihedralElement right) const noexcept
			{
				const int signedRightRotation = left.reflected ? -right.rotation : right.rotation;
				return {normalize(left.rotation + signedRightRotation), left.reflected != right.reflected};
			}

			DihedralElement inverse(DihedralElement element) const noexcept
			{
				return element.reflected ? element : rotation(-element.rotation);
			}

			bool contains(DihedralElement element) const noexcept
			{
				return element.rotation >= 0 && element.rotation < _degree;
			}

			DynamicPermutation permutation(DihedralElement element) const
			{
				if (!contains(element))
					throw ArgumentError("DihedralGroup element is not canonical for this group");
				std::vector<int> images(static_cast<std::size_t>(_degree));
				for (int vertex = 0; vertex < _degree; vertex++) {
					const int transformed = element.reflected ? element.rotation - vertex : element.rotation + vertex;
					images[static_cast<std::size_t>(vertex)] = normalize(transformed);
				}
				return DynamicPermutation(std::move(images));
			}
		};
	}
}

#endif // MML_DIHEDRAL_GROUP_H