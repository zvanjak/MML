///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ICoordTransf.h                                                      ///
///  Description: Coordinate transformation interface                                 ///
///               Abstract base for all coordinate system conversions                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @file ICoordTransf.h
 * @brief Coordinate transformation interfaces for N-dimensional spaces.
 * 
 * This file defines abstract interfaces for coordinate system transformations.
 * These interfaces provide a consistent protocol for converting between
 * different coordinate representations (e.g., Cartesian to spherical).
 * 
 * @note All coordinate transformations are implemented as vector functions
 *       where each output component can be accessed individually.
 */

#if !defined MML_ICoordSystemTransf_H
#define MML_ICoordSystemTransf_H

#include <mml/MMLBase.h>

#include "IFunction.h"

#include <mml/base/Vector/VectorN.h>

#include <type_traits>
#include <utility>

namespace MML
{
	namespace Detail
	{
		/// Deduces N for any type derived from VectorN<T, N> (Vector3Cartesian, Vector2Polar, ...)
		template<typename T, int N>
		std::integral_constant<int, N> DeduceVectorNDimension(const VectorN<T, N>&);

		template<typename V>
		inline constexpr int VectorDimension = decltype(DeduceVectorNDimension(std::declval<const V&>()))::value;
	}

	/**
	 * @brief Interface for coordinate transformations between vector spaces.
	 *
	 * A transformation is a map `VectorFrom -> VectorTo`; its dimensions are **deduced from the
	 * vector types** (`FromDim`, `ToDim`) rather than declared separately, so the base signature
	 * can no longer contradict the template parameters. The base is `IFunction<VectorTo, const
	 * VectorFrom&>`: a transformation is callable as `transf(in)` or `operator()(in)` and a plain
	 * `VectorN<Real, FromDim>` is accepted through the typed vectors' converting constructors.
	 *
	 * The component functions `coordTransfFunc(i)`, `i in [0, ToDim)`, are scalar functions of the
	 * `FromDim` source coordinates; the Jacobian of such a map is `ToDim x FromDim`.
	 *
	 * @tparam VectorFrom Type of the source coordinate vector (derived from VectorN<Real, FromDim>)
	 * @tparam VectorTo   Type of the destination coordinate vector (derived from VectorN<Real, ToDim>)
	 * @tparam N          Optional, defaults to FromDim; kept so existing `ICoordTransf<From, To, 3>`
	 *                    spellings compile. It must equal FromDim (checked at compile time) - it is
	 *                    no longer an independent degree of freedom.
	 *
	 * Example implementations:
	 * - Cartesian to spherical coordinates
	 * - Cartesian to cylindrical coordinates
	 * - General curvilinear coordinate systems
	 */
	template<typename VectorFrom, typename VectorTo, int N = Detail::VectorDimension<VectorFrom>>
	class ICoordTransf : public IFunction<VectorTo, const VectorFrom&>
	{
	public:
		static constexpr int FromDim = Detail::VectorDimension<VectorFrom>;   ///< Dimension of the source space
		static constexpr int ToDim   = Detail::VectorDimension<VectorTo>;     ///< Dimension of the destination space

		static_assert(N == FromDim,
			"ICoordTransf: the explicit dimension N must equal the dimension of VectorFrom "
			"(the dimension is deduced from the vector types; pass only <VectorFrom, VectorTo>)");

		/**
		 * @brief Transform coordinates from source to destination representation.
		 * @param in Input coordinates in source representation
		 * @return Transformed coordinates in destination representation
		 */
		virtual       VectorTo            transf(const VectorFrom& in) const = 0;

		/// Function-call form of transf() (IFunction contract).
		VectorTo operator()(const VectorFrom& in) const override { return transf(in); }

		/**
		 * @brief Get the scalar function for a specific destination coordinate component.
		 * @param i Component index (0 to ToDim-1)
		 * @return Reference to the scalar function of the FromDim source coordinates computing
		 *         the i-th destination component
		 *
		 * @note **Lifetime contract**: The returned reference remains valid for the lifetime of
		 *       this ICoordTransf object. Callers must not cache the reference beyond the
		 *       parent object's lifetime. The transform object owns the underlying functions.
		 */
		virtual const IScalarFunction<FromDim>& coordTransfFunc(int i) const = 0;

		virtual ~ICoordTransf() {}
	};

	/**
	 * @brief Interface for invertible coordinate transformations.
	 *
	 * Extends ICoordTransf with inverse transformation capabilities.
	 * Use this interface when the coordinate transformation is bijective
	 * and both forward and inverse mappings are needed.
	 *
	 * @tparam VectorFrom Type of the source coordinate vector
	 * @tparam VectorTo   Type of the destination coordinate vector
	 * @tparam N          Optional, defaults to FromDim (see ICoordTransf)
	 */
	template<typename VectorFrom, typename VectorTo, int N = Detail::VectorDimension<VectorFrom>>
	class ICoordTransfWithInverse : public virtual ICoordTransf<VectorFrom, VectorTo, N>
	{
	public:
		using ICoordTransf<VectorFrom, VectorTo, N>::FromDim;
		using ICoordTransf<VectorFrom, VectorTo, N>::ToDim;

		/**
		 * @brief Inverse transformation from destination back to source coordinates.
		 * @param in Input coordinates in destination representation
		 * @return Transformed coordinates in source representation
		 */
		virtual       VectorFrom          transfInverse(const VectorTo& in) const = 0;

		/**
		 * @brief Get the scalar function for a specific inverse coordinate component.
		 * @param i Component index (0 to FromDim-1)
		 * @return Reference to the scalar function of the ToDim destination coordinates computing
		 *         the i-th source component
		 *
		 * @note **Lifetime contract**: The returned reference remains valid for the lifetime of
		 *       this ICoordTransfWithInverse object. Callers must not cache the reference beyond
		 *       the parent object's lifetime. The transform object owns the underlying functions.
		 */
		virtual const IScalarFunction<ToDim>& inverseCoordTransfFunc(int i) const = 0;

		virtual ~ICoordTransfWithInverse() {}
	};
}
#endif