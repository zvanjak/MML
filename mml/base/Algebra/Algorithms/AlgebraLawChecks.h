///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        AlgebraLawChecks.h                                                  ///
///  Description: Exhaustive law checks for small finite algebraic structures         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ALGEBRA_LAW_CHECKS_H
#define MML_ALGEBRA_LAW_CHECKS_H

#include <mml/base/Algebra/AlgebraTraits.h>

#include <functional>

namespace MML
{
	namespace Algebra
	{
		///////////////////////////////////////////////////////////////////////////////////
		// Finite structure protocol
		//
		// Group: elements(), identity(), compose(left, right), inverse(element)
		// Ring:  elements(), zero(), one(), add(left, right), negate(element),
		//        multiply(left, right)
		// Field: ring protocol plus inverse(element) for nonzero elements
		//
		// Composition convention: compose(left, right) returns left * right. For
		// transformations this means that right is applied first, then left.
		//
		// These checks are exhaustive and therefore O(n^3). They are intended for
		// small finite structures and tests, not large-group verification.
		///////////////////////////////////////////////////////////////////////////////////

		namespace Detail
		{
			template<class Elements, class Element, class Equal>
			bool Contains(const Elements& elements, const Element& value, Equal equal)
			{
				for (const auto& element : elements)
					if (equal(element, value))
						return true;
				return false;
			}

			template<class Elements, class Operation, class Equal>
			bool CheckClosureForOperation(const Elements& elements, Operation operation, Equal equal)
			{
				for (const auto& left : elements)
					for (const auto& right : elements)
						if (!Contains(elements, operation(left, right), equal))
							return false;
				return true;
			}

			template<class Elements, class Operation, class Equal>
			bool CheckAssociativityForOperation(const Elements& elements, Operation operation, Equal equal)
			{
				for (const auto& first : elements)
					for (const auto& second : elements)
						for (const auto& third : elements)
							if (!equal(operation(operation(first, second), third),
								operation(first, operation(second, third))))
								return false;
				return true;
			}
		}

		template<class Structure,
			class Equal = std::equal_to<typename AlgebraTraits<Structure>::element_type>>
		bool CheckClosure(const Structure& structure, Equal equal = Equal())
		{
			const auto& elements = structure.elements();
			return Detail::CheckClosureForOperation(elements,
				[&structure](const auto& left, const auto& right) { return structure.compose(left, right); }, equal);
		}

		template<class Structure,
			class Equal = std::equal_to<typename AlgebraTraits<Structure>::element_type>>
		bool CheckAssociativity(const Structure& structure, Equal equal = Equal())
		{
			const auto& elements = structure.elements();
			return Detail::CheckAssociativityForOperation(elements,
				[&structure](const auto& left, const auto& right) { return structure.compose(left, right); }, equal);
		}

		template<class Structure,
			class Equal = std::equal_to<typename AlgebraTraits<Structure>::element_type>>
		bool CheckIdentity(const Structure& structure, Equal equal = Equal())
		{
			const auto& elements = structure.elements();
			const auto identity = structure.identity();
			if (!Detail::Contains(elements, identity, equal))
				return false;

			for (const auto& element : elements)
				if (!equal(structure.compose(identity, element), element) ||
					!equal(structure.compose(element, identity), element))
					return false;
			return true;
		}

		template<class Structure,
			class Equal = std::equal_to<typename AlgebraTraits<Structure>::element_type>>
		bool CheckInverses(const Structure& structure, Equal equal = Equal())
		{
			const auto& elements = structure.elements();
			const auto identity = structure.identity();
			for (const auto& element : elements) {
				const auto inverse = structure.inverse(element);
				if (!Detail::Contains(elements, inverse, equal) ||
					!equal(structure.compose(element, inverse), identity) ||
					!equal(structure.compose(inverse, element), identity))
					return false;
			}
			return true;
		}

		template<class Structure,
			class Equal = std::equal_to<typename AlgebraTraits<Structure>::element_type>>
		bool CheckGroupLaws(const Structure& structure, Equal equal = Equal())
		{
			return CheckClosure(structure, equal) &&
				CheckAssociativity(structure, equal) &&
				CheckIdentity(structure, equal) &&
				CheckInverses(structure, equal);
		}

		template<class Structure,
			class Equal = std::equal_to<typename AlgebraTraits<Structure>::element_type>>
		bool CheckRingLaws(const Structure& structure, Equal equal = Equal())
		{
			const auto& elements = structure.elements();
			const auto zero = structure.zero();
			const auto one = structure.one();
			auto add = [&structure](const auto& left, const auto& right) { return structure.add(left, right); };
			auto multiply = [&structure](const auto& left, const auto& right) { return structure.multiply(left, right); };

			if (!Detail::Contains(elements, zero, equal) || !Detail::Contains(elements, one, equal) ||
				!Detail::CheckClosureForOperation(elements, add, equal) ||
				!Detail::CheckClosureForOperation(elements, multiply, equal) ||
				!Detail::CheckAssociativityForOperation(elements, add, equal) ||
				!Detail::CheckAssociativityForOperation(elements, multiply, equal))
				return false;

			for (const auto& left : elements) {
				if (!equal(add(left, zero), left) || !equal(add(zero, left), left) ||
					!equal(multiply(left, one), left) || !equal(multiply(one, left), left))
					return false;

				const auto negative = structure.negate(left);
				if (!Detail::Contains(elements, negative, equal) ||
					!equal(add(left, negative), zero) || !equal(add(negative, left), zero))
					return false;

				for (const auto& right : elements) {
					if (!equal(add(left, right), add(right, left)))
						return false;
					for (const auto& third : elements) {
						if (!equal(multiply(left, add(right, third)),
							add(multiply(left, right), multiply(left, third))) ||
							!equal(multiply(add(left, right), third),
								add(multiply(left, third), multiply(right, third))))
							return false;
					}
				}
			}
			return true;
		}

		template<class Structure,
			class Equal = std::equal_to<typename AlgebraTraits<Structure>::element_type>>
		bool CheckFieldLaws(const Structure& structure, Equal equal = Equal())
		{
			if (!CheckRingLaws(structure, equal))
				return false;

			const auto& elements = structure.elements();
			const auto zero = structure.zero();
			const auto one = structure.one();
			for (const auto& left : elements) {
				for (const auto& right : elements)
					if (!equal(structure.multiply(left, right), structure.multiply(right, left)))
						return false;

				if (!equal(left, zero)) {
					const auto inverse = structure.inverse(left);
					if (!Detail::Contains(elements, inverse, equal) ||
						!equal(structure.multiply(left, inverse), one) ||
						!equal(structure.multiply(inverse, left), one))
						return false;
				}
			}
			return true;
		}
	}
}

#endif // MML_ALGEBRA_LAW_CHECKS_H