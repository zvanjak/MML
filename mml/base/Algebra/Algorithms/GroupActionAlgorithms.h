///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        GroupActionAlgorithms.h                                             ///
///  Description: Orbits, stabilizers, invariants, and Burnside counting              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GROUP_ACTION_ALGORITHMS_H
#define MML_GROUP_ACTION_ALGORITHMS_H

#include <mml/MMLExceptions.h>
#include <mml/base/Algebra/AlgebraTraits.h>

#include <cstddef>
#include <functional>
#include <vector>

namespace MML
{
	namespace Algebra
	{
		namespace Detail
		{
			template<class Object, class Equal>
			bool ActionContains(const std::vector<Object>& objects, const Object& object, Equal equal)
			{
				for (const auto& candidate : objects)
					if (equal(candidate, object))
						return true;
				return false;
			}

			template<class Object, class Equal>
			void AppendActionUnique(std::vector<Object>& objects, const Object& object, Equal equal)
			{
				if (!ActionContains(objects, object, equal))
					objects.push_back(object);
			}
		}

		template<class Group, class Action, class Object,
			class Equal = std::equal_to<Object>>
		std::vector<Object> Orbit(const Group& group, const Action& action,
			const Object& object, Equal equal = Equal())
		{
			std::vector<Object> orbit;
			for (const auto& element : group.elements())
				Detail::AppendActionUnique(orbit, action.apply(element, object), equal);
			return orbit;
		}

		template<class Group, class Action, class Object,
			class Equal = std::equal_to<Object>>
		std::vector<typename AlgebraTraits<Group>::element_type> Stabilizer(
			const Group& group, const Action& action, const Object& object, Equal equal = Equal())
		{
			std::vector<typename AlgebraTraits<Group>::element_type> stabilizer;
			for (const auto& element : group.elements())
				if (equal(action.apply(element, object), object))
					stabilizer.push_back(element);
			return stabilizer;
		}

		template<class Group, class Action, class Object,
			class Relation = std::equal_to<Object>>
		bool IsInvariant(const Group& group, const Action& action,
			const Object& object, Relation relation = Relation())
		{
			for (const auto& element : group.elements())
				if (!relation(action.apply(element, object), object))
					return false;
			return true;
		}

		template<class Group, class Action, class Object, class Predicate>
		bool IsPredicateInvariant(const Group& group, const Action& action,
			const std::vector<Object>& objects, Predicate predicate)
		{
			for (const auto& object : objects)
				for (const auto& element : group.elements())
					if (static_cast<bool>(predicate(action.apply(element, object))) !=
						static_cast<bool>(predicate(object)))
						return false;
			return true;
		}

		template<class Group, class Action, class Object,
			class Equal = std::equal_to<Object>>
		std::vector<Object> FixedPoints(const Group& group, const Action& action,
			const typename AlgebraTraits<Group>::element_type& element,
			const std::vector<Object>& objects, Equal equal = Equal())
		{
			std::vector<Object> fixedPoints;
			for (const auto& object : objects)
				if (equal(action.apply(element, object), object))
					fixedPoints.push_back(object);
			return fixedPoints;
		}

		template<class Group, class Action, class Object,
			class Equal = std::equal_to<Object>>
		std::size_t FixedPointCount(const Group& group, const Action& action,
			const typename AlgebraTraits<Group>::element_type& element,
			const std::vector<Object>& objects, Equal equal = Equal())
		{
			return FixedPoints(group, action, element, objects, equal).size();
		}

		template<class Group, class Action, class Object,
			class Equal = std::equal_to<Object>>
		std::vector<std::vector<Object>> OrbitPartition(const Group& group, const Action& action,
			const std::vector<Object>& objects, Equal equal = Equal())
		{
			std::vector<std::vector<Object>> partition;
			std::vector<Object> assigned;
			for (const auto& object : objects) {
				if (Detail::ActionContains(assigned, object, equal))
					continue;
				auto orbit = Orbit(group, action, object, equal);
				for (const auto& member : orbit) {
					if (!Detail::ActionContains(objects, member, equal))
						throw DomainError("Group action is not closed on the supplied object set");
					Detail::AppendActionUnique(assigned, member, equal);
				}
				partition.push_back(std::move(orbit));
			}
			return partition;
		}

		template<class Group, class Action, class Object,
			class Equal = std::equal_to<Object>>
		std::size_t BurnsideCount(const Group& group, const Action& action,
			const std::vector<Object>& objects, Equal equal = Equal())
		{
			const auto& elements = group.elements();
			if (elements.empty())
				throw ArgumentError("BurnsideCount requires a non-empty finite group");

			std::size_t fixedPointTotal = 0;
			for (const auto& element : elements) {
				for (const auto& object : objects) {
					const auto transformed = action.apply(element, object);
					if (!Detail::ActionContains(objects, transformed, equal))
						throw DomainError("Group action is not closed on the supplied object set");
					if (equal(transformed, object))
						fixedPointTotal++;
				}
			}

			if (fixedPointTotal % elements.size() != 0)
				throw DomainError("Burnside fixed-point average is not integral");
			return fixedPointTotal / elements.size();
		}
	}
}

#endif // MML_GROUP_ACTION_ALGORITHMS_H