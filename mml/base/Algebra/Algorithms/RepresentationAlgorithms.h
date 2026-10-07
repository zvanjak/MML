///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        RepresentationAlgorithms.h                                          ///
///  Description: Characters, homomorphism checks, and invariant averaging            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_REPRESENTATION_ALGORITHMS_H
#define MML_REPRESENTATION_ALGORITHMS_H

#include <mml/MMLExceptions.h>
#include <mml/base/Algebra/AlgebraTraits.h>
#include <mml/base/Algebra/Algorithms/FiniteGroupAlgorithms.h>

#include <cstddef>
#include <functional>
#include <utility>

namespace MML
{
	namespace Algebra
	{
		template<class Group, class RepresentationType, class MatrixEqual>
		bool VerifyRepresentation(const Group& group, const RepresentationType& representation,
			MatrixEqual matrixEqual)
		{
			using Matrix = typename RepresentationType::matrix_type;
			if (!matrixEqual(representation.matrix(group.identity()), Matrix::Identity()))
				return false;

			for (const auto& left : group.elements())
				for (const auto& right : group.elements())
					if (!matrixEqual(representation.matrix(group.compose(left, right)),
						representation.matrix(left) * representation.matrix(right)))
						return false;
			return true;
		}

		template<class Group, class RepresentationType>
		bool VerifyRepresentation(const Group& group, const RepresentationType& representation)
		{
			return VerifyRepresentation(group, representation,
				[](const auto& left, const auto& right) { return left == right; });
		}

		template<class RepresentationType>
		typename RepresentationType::scalar_type Character(const RepresentationType& representation,
			const typename RepresentationType::group_element_type& element)
		{
			return representation.matrix(element).Trace();
		}

		template<class Group, class RepresentationType, class ScalarEqual,
			class ElementEqual = std::equal_to<typename AlgebraTraits<Group>::element_type>>
		bool IsCharacterConstantOnConjugacyClasses(const Group& group,
			const RepresentationType& representation, ScalarEqual scalarEqual,
			ElementEqual elementEqual = ElementEqual())
		{
			for (const auto& conjugacyClass : ConjugacyClasses(group, elementEqual)) {
				if (conjugacyClass.empty())
					continue;
				const auto expected = Character(representation, conjugacyClass.front());
				for (const auto& element : conjugacyClass)
					if (!scalarEqual(Character(representation, element), expected))
						return false;
			}
			return true;
		}

		template<class Group, class RepresentationType>
		bool IsCharacterConstantOnConjugacyClasses(const Group& group,
			const RepresentationType& representation)
		{
			return IsCharacterConstantOnConjugacyClasses(group, representation,
				[](const auto& left, const auto& right) { return left == right; });
		}

		template<class Group, class RepresentationType>
		typename RepresentationType::matrix_type InvariantProjection(
			const Group& group, const RepresentationType& representation)
		{
			using Matrix = typename RepresentationType::matrix_type;
			using Scalar = typename RepresentationType::scalar_type;
			if (group.elements().empty())
				throw ArgumentError("InvariantProjection requires a non-empty finite group");
			Matrix projection;
			for (const auto& element : group.elements())
				projection = projection + representation.matrix(element);
			return projection / Scalar(group.elements().size());
		}

		template<class Group, class Action, class Object, class Add, class Scale>
		Object GroupAverage(const Group& group, const Action& action, const Object& object,
			Object zero, Add add, Scale scale)
		{
			if (group.elements().empty())
				throw ArgumentError("GroupAverage requires a non-empty finite group");
			Object sum = std::move(zero);
			for (const auto& element : group.elements())
				sum = add(std::move(sum), action(element, object));
			return scale(std::move(sum), group.elements().size());
		}

		template<class Group, class RepresentationType>
		typename RepresentationType::vector_type SymmetrizeVector(const Group& group,
			const RepresentationType& representation,
			const typename RepresentationType::vector_type& vector)
		{
			return InvariantProjection(group, representation) * vector;
		}
	}
}

#endif // MML_REPRESENTATION_ALGORITHMS_H