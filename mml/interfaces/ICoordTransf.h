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

#include "MMLBase.h"

#include "IFunction.h"

#include "base/Vector/VectorN.h"

namespace MML
{
	/**
	 * @brief Interface for coordinate transformations between vector spaces.
	 * 
	 * This template class defines the protocol for transforming coordinates
	 * from one representation to another. Implementations must provide both
	 * the transformation and access to individual component functions.
	 * 
	 * @tparam VectorFrom Type of the source coordinate vector
	 * @tparam VectorTo   Type of the destination coordinate vector
	 * @tparam N          Dimension of the coordinate space
	 * 
	 * Example implementations:
	 * - Cartesian to spherical coordinates
	 * - Cartesian to cylindrical coordinates
	 * - General curvilinear coordinate systems
	 */
	template<typename VectorFrom, typename VectorTo, int N>
	class ICoordTransf : public IVectorFunction<N>
	{
	public:
		/**
		 * @brief Transform coordinates from source to destination representation.
		 * @param in Input coordinates in source representation
		 * @return Transformed coordinates in destination representation
		 */
		virtual       VectorTo            transf(const VectorFrom& in) const = 0;
		
		/**
		 * @brief Get the scalar function for a specific coordinate component.
		 * @param i Component index (0 to N-1)
		 * @return Reference to the scalar function computing the i-th output component
		 * 
		 * @note **Lifetime contract**: The returned reference remains valid for the lifetime of
		 *       this ICoordTransf object. Callers must not cache the reference beyond the
		 *       parent object's lifetime. The transform object owns the underlying functions.
		 */
		virtual const IScalarFunction<N>& coordTransfFunc(int i) const = 0;

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
	 * @tparam N          Dimension of the coordinate space
	 */
	template<typename VectorFrom, typename VectorTo, int N>
	class ICoordTransfWithInverse : public virtual ICoordTransf<VectorFrom, VectorTo, N>
	{
	public:
		/**
		 * @brief Inverse transformation from destination back to source coordinates.
		 * @param in Input coordinates in destination representation
		 * @return Transformed coordinates in source representation
		 */
		virtual       VectorFrom          transfInverse(const VectorTo& in) const = 0;
		
		/**
		 * @brief Get the scalar function for a specific inverse coordinate component.
		 * @param i Component index (0 to N-1)
		 * @return Reference to the scalar function computing the i-th inverse component
		 * 
		 * @note **Lifetime contract**: The returned reference remains valid for the lifetime of
		 *       this ICoordTransfWithInverse object. Callers must not cache the reference beyond
		 *       the parent object's lifetime. The transform object owns the underlying functions.
		 */
		virtual const IScalarFunction<N>& inverseCoordTransfFunc(int i) const = 0;

		virtual ~ICoordTransfWithInverse() {}
	};
}
#endif