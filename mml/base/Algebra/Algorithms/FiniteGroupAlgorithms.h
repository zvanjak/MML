///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FiniteGroupAlgorithms.h                                             ///
///  Description: Generic algorithms for small finite groups                          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FINITE_GROUP_ALGORITHMS_H
#define MML_FINITE_GROUP_ALGORITHMS_H

#include <mml/MMLExceptions.h>
#include <mml/base/Algebra/AlgebraTraits.h>

#include <cstddef>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace MML
{
	namespace Algebra
	{
		template<class Element>
		struct CayleyGraphEdge
		{
			std::size_t from = 0;
			std::size_t to = 0;
			std::size_t generator = 0;
		};

		template<class Element>
		struct CayleyGraphData
		{
			std::vector<Element> elements;
			std::vector<Element> generators;
			std::vector<CayleyGraphEdge<Element>> edges;
		};

		namespace Detail
		{
			template<class Elements, class Element, class Equal>
			std::size_t FindElementIndex(const Elements& elements, const Element& element, Equal equal)
			{
				for (std::size_t index = 0; index < elements.size(); index++)
					if (equal(elements[index], element))
						return index;
				throw DomainError("Finite-group operation produced an element outside the group");
			}

			template<class Element, class Equal>
			bool AppendUnique(std::vector<Element>& elements, const Element& element, Equal equal)
			{
				for (const auto& existing : elements)
					if (equal(existing, element))
						return false;
				elements.push_back(element);
				return true;
			}
		}

		template<class Group,
			class Equal = std::equal_to<typename AlgebraTraits<Group>::element_type>>
		std::size_t element_order(const Group& group,
			const typename AlgebraTraits<Group>::element_type& element, Equal equal = Equal())
		{
			const auto& elements = group.elements();
			Detail::FindElementIndex(elements, element, equal);
			const auto identity = group.identity();
			auto power = identity;
			for (std::size_t exponent = 1; exponent <= elements.size(); exponent++) {
				power = group.compose(power, element);
				if (equal(power, identity))
					return exponent;
			}
			throw DomainError("Element order exceeds the finite group order");
		}

		template<class Group,
			class Equal = std::equal_to<typename AlgebraTraits<Group>::element_type>>
		std::vector<std::vector<std::size_t>> CayleyTable(const Group& group, Equal equal = Equal())
		{
			const auto& elements = group.elements();
			std::vector<std::vector<std::size_t>> table(elements.size(),
				std::vector<std::size_t>(elements.size()));
			for (std::size_t row = 0; row < elements.size(); row++)
				for (std::size_t column = 0; column < elements.size(); column++)
					table[row][column] = Detail::FindElementIndex(elements,
						group.compose(elements[row], elements[column]), equal);
			return table;
		}

		template<class Group,
			class Equal = std::equal_to<typename AlgebraTraits<Group>::element_type>>
		std::vector<typename AlgebraTraits<Group>::element_type> GeneratedSubgroup(const Group& group,
			const std::vector<typename AlgebraTraits<Group>::element_type>& generators, Equal equal = Equal())
		{
			using Element = typename AlgebraTraits<Group>::element_type;
			const auto& groupElements = group.elements();
			std::vector<Element> subgroup{group.identity()};

			for (const auto& generator : generators) {
				Detail::FindElementIndex(groupElements, generator, equal);
				Detail::AppendUnique(subgroup, generator, equal);
				Detail::AppendUnique(subgroup, group.inverse(generator), equal);
			}

			bool changed = true;
			while (changed) {
				changed = false;
				const std::vector<Element> snapshot = subgroup;
				for (const auto& left : snapshot)
					for (const auto& right : snapshot) {
						const auto product = group.compose(left, right);
						Detail::FindElementIndex(groupElements, product, equal);
						changed = Detail::AppendUnique(subgroup, product, equal) || changed;
					}
			}
			return subgroup;
		}

		template<class Group,
			class Equal = std::equal_to<typename AlgebraTraits<Group>::element_type>>
		std::vector<std::vector<typename AlgebraTraits<Group>::element_type>> ConjugacyClasses(
			const Group& group, Equal equal = Equal())
		{
			using Element = typename AlgebraTraits<Group>::element_type;
			const auto& elements = group.elements();
			std::vector<std::vector<Element>> classes;
			std::vector<Element> assigned;

			for (const auto& element : elements) {
				bool alreadyAssigned = false;
				for (const auto& existing : assigned)
					if (equal(existing, element)) {
						alreadyAssigned = true;
						break;
					}
				if (alreadyAssigned)
					continue;

				std::vector<Element> conjugacyClass;
				for (const auto& conjugator : elements) {
					const auto conjugate = group.compose(
						group.compose(conjugator, element), group.inverse(conjugator));
					Detail::FindElementIndex(elements, conjugate, equal);
					Detail::AppendUnique(conjugacyClass, conjugate, equal);
					Detail::AppendUnique(assigned, conjugate, equal);
				}
				classes.push_back(std::move(conjugacyClass));
			}
			return classes;
		}

		template<class Group,
			class Equal = std::equal_to<typename AlgebraTraits<Group>::element_type>>
		CayleyGraphData<typename AlgebraTraits<Group>::element_type> MakeCayleyGraphData(
			const Group& group,
			const std::vector<typename AlgebraTraits<Group>::element_type>& generators,
			Equal equal = Equal())
		{
			using Element = typename AlgebraTraits<Group>::element_type;
			CayleyGraphData<Element> data;
			data.elements.assign(group.elements().begin(), group.elements().end());
			data.generators = generators;

			for (const auto& generator : generators)
				Detail::FindElementIndex(data.elements, generator, equal);

			for (std::size_t from = 0; from < data.elements.size(); from++)
				for (std::size_t generator = 0; generator < generators.size(); generator++) {
					const auto destination = group.compose(generators[generator], data.elements[from]);
					data.edges.push_back({from, Detail::FindElementIndex(data.elements, destination, equal), generator});
				}
			return data;
		}
	}
}

#endif // MML_FINITE_GROUP_ALGORITHMS_H