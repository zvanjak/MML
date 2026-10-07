///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        AlgebraTraits.h                                                     ///
///  Description: Metadata and equality policy for algebraic structures               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ALGEBRA_TRAITS_H
#define MML_ALGEBRA_TRAITS_H

#include <cstddef>
#include <limits>
#include <type_traits>

namespace MML
{
	namespace Algebra
	{
		inline constexpr std::size_t DynamicAlgebraExtent = std::numeric_limits<std::size_t>::max();

		enum class AlgebraEquality
		{
			Exact,
			Approximate
		};

		namespace Detail
		{
			template<class Structure, class = void>
			struct AlgebraScalarType
			{
				using type = typename Structure::element_type;
			};

			template<class Structure>
			struct AlgebraScalarType<Structure, std::void_t<typename Structure::scalar_type>>
			{
				using type = typename Structure::scalar_type;
			};
		}

		/// @brief Customization point for algebraic structure metadata.
		/// @details Concrete structures must expose element_type. They may expose scalar_type,
		///          or specialize this trait to provide fixed order, dimension, and equality policy.
		template<class Structure>
		struct AlgebraTraits
		{
			using structure_type = Structure;
			using element_type = typename Structure::element_type;
			using scalar_type = typename Detail::AlgebraScalarType<Structure>::type;

			static constexpr std::size_t static_order = DynamicAlgebraExtent;
			static constexpr std::size_t dimension = DynamicAlgebraExtent;
			static constexpr AlgebraEquality equality = AlgebraEquality::Exact;
		};
	}
}

#endif // MML_ALGEBRA_TRAITS_H