///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FiniteGroup.h                                                       ///
///  Description: Runtime finite-group value wrapper                                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FINITE_GROUP_H
#define MML_FINITE_GROUP_H

#include <mml/MMLExceptions.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

namespace MML
{
	namespace Algebra
	{
		/// @brief A small finite set together with exact group operations.
		/// @details This Base type stores structure only. Use Core algebra law checks to
		///          validate closure and the group axioms.
		template<class Element, class Equal = std::equal_to<Element>>
		class FiniteGroup
		{
		public:
			using element_type = Element;
			using compose_type = std::function<Element(const Element&, const Element&)>;
			using inverse_type = std::function<Element(const Element&)>;

		private:
			std::vector<Element> _elements;
			Element _identity;
			compose_type _compose;
			inverse_type _inverse;
			Equal _equal;

		public:
			FiniteGroup(std::vector<Element> elements, Element identity,
				compose_type composeOperation, inverse_type inverseOperation, Equal equal = Equal())
				: _elements(std::move(elements)), _identity(std::move(identity)),
				  _compose(std::move(composeOperation)), _inverse(std::move(inverseOperation)),
				  _equal(std::move(equal))
			{
				if (_elements.empty())
					throw ArgumentError("FiniteGroup requires at least one element");
				if (!_compose || !_inverse)
					throw ArgumentError("FiniteGroup operations must be callable");
				if (!contains(_identity))
					throw ArgumentError("FiniteGroup elements must contain the identity");

				for (std::size_t left = 0; left < _elements.size(); left++)
					for (std::size_t right = left + 1; right < _elements.size(); right++)
						if (_equal(_elements[left], _elements[right]))
							throw ArgumentError("FiniteGroup elements must be unique");
			}

			const std::vector<Element>& elements() const noexcept { return _elements; }
			const Element& identity() const noexcept { return _identity; }
			std::size_t order() const noexcept { return _elements.size(); }

			Element compose(const Element& left, const Element& right) const
			{
				return _compose(left, right);
			}

			Element inverse(const Element& element) const
			{
				return _inverse(element);
			}

			bool contains(const Element& element) const
			{
				return std::any_of(_elements.begin(), _elements.end(),
					[this, &element](const Element& candidate) { return _equal(candidate, element); });
			}
		};
	}
}

#endif // MML_FINITE_GROUP_H