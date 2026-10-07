///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        GroupAction.h                                                       ///
///  Description: Callable value adapter for group actions                            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GROUP_ACTION_H
#define MML_GROUP_ACTION_H

#include <mml/MMLExceptions.h>

#include <functional>
#include <utility>

namespace MML
{
	namespace Algebra
	{
		/// @brief Type-erased left action of group elements on objects.
		/// @details apply(g, object) computes g.object. A valid action satisfies
		///          apply(identity, x) = x and apply(compose(g, h), x) =
		///          apply(g, apply(h, x)).
		template<class GroupElement, class Object>
		class GroupAction
		{
		public:
			using group_element_type = GroupElement;
			using object_type = Object;
			using apply_type = std::function<Object(const GroupElement&, const Object&)>;

		private:
			apply_type _apply;

		public:
			explicit GroupAction(apply_type applyOperation)
				: _apply(std::move(applyOperation))
			{
				if (!_apply)
					throw ArgumentError("GroupAction operation must be callable");
			}

			Object apply(const GroupElement& element, const Object& object) const
			{
				return _apply(element, object);
			}

			Object operator()(const GroupElement& element, const Object& object) const
			{
				return apply(element, object);
			}
		};

		template<class GroupElement, class Object, class Apply>
		GroupAction<GroupElement, Object> MakeGroupAction(Apply applyOperation)
		{
			return GroupAction<GroupElement, Object>(std::move(applyOperation));
		}
	}
}

#endif // MML_GROUP_ACTION_H