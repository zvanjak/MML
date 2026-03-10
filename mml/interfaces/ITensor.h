///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ITensor.h                                                           ///
///  Description: Interface for tensor objects (contravariant/covariant indices)      ///
///               Base class for tensor algebra operations                            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @file ITensor.h
 * @brief Interfaces for tensor objects with index variance tracking.
 * 
 * Defines abstract interfaces for tensors of rank 2 through 5 in N-dimensional
 * space. Each tensor interface tracks the number of contravariant (upper) and
 * covariant (lower) indices, enabling proper transformation behavior under
 * coordinate changes.
 * 
 * **Tensor Notation:**
 * - Contravariant indices (superscripts): Transform like coordinate differentials dx^i
 * - Covariant indices (subscripts): Transform like gradient components ∂f/∂x_i
 * - Mixed tensors: Have both types, e.g., T^i_j (1 contravariant, 1 covariant)
 * 
 * **Rank Hierarchy:**
 * - ITensor2<N>: Rank-2 tensors (matrices), e.g., metric tensor g_ij
 * - ITensor3<N>: Rank-3 tensors, e.g., Christoffel symbols Γ^i_jk
 * - ITensor4<N>: Rank-4 tensors, e.g., Riemann curvature R^i_jkl
 * - ITensor5<N>: Rank-5 tensors for specialized applications
 * 
 * @see ITensorField.h for position-dependent tensor fields
 */

#if !defined MML_ITENSOR_H
#define MML_ITENSOR_H

#include "MMLBase.h"
#include "base/Vector/VectorN.h"

namespace MML
{
	// Forward declarations for contraction return types
	template<int N> class Tensor2;
	template<int N> class Tensor3;

	/**
	 * @brief Enumeration for tensor index variance type.
	 */
	enum TensorIndexType { 
		CONTRAVARIANT,  ///< Upper index, transforms with inverse Jacobian
		COVARIANT       ///< Lower index, transforms with Jacobian
	};
    
	/**
	 * @brief Interface for rank-2 tensors in N-dimensional space.
	 * 
	 * Represents tensors with 2 indices, such as:
	 * - Metric tensor g_ij (0 contravariant, 2 covariant)
	 * - Inverse metric g^ij (2 contravariant, 0 covariant)
	 * - Mixed tensor T^i_j (1 contravariant, 1 covariant)
	 * 
	 * @tparam N Dimension of the underlying space
	 */
	template<int N>
	class ITensor2
	{
	public:
		virtual ~ITensor2() = default;

		/** @brief Get the number of contravariant (upper) indices. */
		virtual int   NumContravar() const = 0;
		
		/** @brief Get the number of covariant (lower) indices. */
		virtual int   NumCovar() const = 0;

		/**
		 * @brief Access tensor component (const).
		 * @param i First index (0 to N-1)
		 * @param j Second index (0 to N-1)
		 * @return The tensor component T_ij or T^ij
		 */
		virtual Real  operator()(int i, int j) const = 0;
		
		/**
		 * @brief Access tensor component (mutable).
		 * @param i First index (0 to N-1)
		 * @param j Second index (0 to N-1)
		 * @return Reference to the tensor component
		 */
		virtual Real& operator()(int i, int j) = 0;

		/**
		 * @brief Contract the tensor (trace operation for mixed tensors).
		 * 
		 * For a mixed tensor T^i_j, returns the trace: sum of T^i_i over i.
		 * This is the fundamental contraction operation that produces a scalar.
		 * 
		 * @return Scalar result of contraction (trace)
		 * @throws TensorCovarContravarNumError if tensor is not mixed type (1 up, 1 down)
		 */
		virtual Real Contract() const = 0;

		/**
		 * @brief Apply tensor as a multilinear form to two vectors.
		 * 
		 * Computes the contraction T(v1, v2) = sum_{i,j} T_ij * v1_i * v2_j.
		 * This evaluates the tensor as a bilinear form on the given vectors.
		 * 
		 * @param v1 First vector argument
		 * @param v2 Second vector argument
		 * @return Scalar result of the bilinear form evaluation
		 */
		virtual Real operator()(const VectorN<Real, N>& v1, const VectorN<Real, N>& v2) const = 0;
	};

	/**
	 * @brief Interface for rank-3 tensors in N-dimensional space.
	 * 
	 * Represents tensors with 3 indices, commonly used for:
	 * - Christoffel symbols Γ^i_jk (connection coefficients)
	 * - Torsion tensors T^i_jk
	 * - Structure constants of Lie algebras
	 * 
	 * @tparam N Dimension of the underlying space
	 */
	template<int N>
	class ITensor3
	{
	public:
		virtual ~ITensor3() = default;

		virtual int   NumContravar() const = 0;
		virtual int   NumCovar() const = 0;

		/**
		 * @brief Access tensor component.
		 * @param i First index
		 * @param j Second index
		 * @param k Third index
		 */
		virtual Real  operator()(int i, int j, int k) const = 0;
		virtual Real& operator()(int i, int j, int k) = 0;

		/**
		 * @brief Contract two indices, producing a rank-1 tensor (vector).
		 * 
		 * Sums over the contracted index: result_a = sum_s T_{...s...s...}
		 * where the two 's' positions are determined by ind1 and ind2.
		 * 
		 * @param ind1 First index to contract (0, 1, or 2)
		 * @param ind2 Second index to contract (0, 1, or 2), must differ from ind1
		 * @return VectorN<Real, N> - the contracted rank-1 result
		 * @throws TensorIndexError if indices out of range or equal
		 */
		virtual VectorN<Real, N> Contract(int ind1, int ind2) const = 0;

		/**
		 * @brief Apply tensor as a multilinear form to three vectors.
		 * 
		 * Computes T(v1, v2, v3) = sum_{i,j,k} T_ijk * v1_i * v2_j * v3_k.
		 * This evaluates the tensor as a trilinear form on the given vectors.
		 * 
		 * @param v1 First vector argument
		 * @param v2 Second vector argument
		 * @param v3 Third vector argument
		 * @return Scalar result of the trilinear form evaluation
		 */
		virtual Real operator()(const VectorN<Real, N>& v1, const VectorN<Real, N>& v2,
		                         const VectorN<Real, N>& v3) const = 0;
	};

	/**
	 * @brief Interface for rank-4 tensors in N-dimensional space.
	 * 
	 * Represents tensors with 4 indices, most notably:
	 * - Riemann curvature tensor R^i_jkl
	 * - Elasticity tensor C_ijkl (in continuum mechanics)
	 * - Weyl conformal tensor
	 * 
	 * @tparam N Dimension of the underlying space
	 */
	template<int N>
	class ITensor4
	{
	public:
		virtual ~ITensor4() = default;

		virtual int   NumContravar() const = 0;
		virtual int   NumCovar() const = 0;

		/**
		 * @brief Access tensor component.
		 * @param i First index
		 * @param j Second index
		 * @param k Third index
		 * @param l Fourth index
		 */
		virtual Real  operator()(int i, int j, int k, int l) const = 0;
		virtual Real& operator()(int i, int j, int k, int l) = 0;

		/**
		 * @brief Contract two indices, producing a rank-2 tensor.
		 * 
		 * Sums over the contracted index: result_{ab} = sum_s T_{...s...s...}
		 * For example, contracting indices 0 and 2 of Riemann R^i_jkl gives Ricci R_jl.
		 * 
		 * @param ind1 First index to contract (0, 1, 2, or 3)
		 * @param ind2 Second index to contract (0, 1, 2, or 3), must differ from ind1
		 * @return Tensor2<N> - the contracted rank-2 result
		 * @throws TensorIndexError if indices out of range or equal
		 */
		virtual Tensor2<N> Contract(int ind1, int ind2) const = 0;

		/**
		 * @brief Apply tensor as a multilinear form to four vectors.
		 * 
		 * Computes T(v1, v2, v3, v4) = sum_{i,j,k,l} T_ijkl * v1_i * v2_j * v3_k * v4_l.
		 * This evaluates the tensor as a quadrilinear form on the given vectors.
		 * 
		 * @param v1 First vector argument
		 * @param v2 Second vector argument
		 * @param v3 Third vector argument
		 * @param v4 Fourth vector argument
		 * @return Scalar result of the quadrilinear form evaluation
		 */
		virtual Real operator()(const VectorN<Real, N>& v1, const VectorN<Real, N>& v2,
		                         const VectorN<Real, N>& v3, const VectorN<Real, N>& v4) const = 0;
	};

	/**
	 * @brief Interface for rank-5 tensors in N-dimensional space.
	 * 
	 * Represents tensors with 5 indices for specialized applications
	 * such as higher-order derivatives or multi-field theories.
	 * 
	 * @tparam N Dimension of the underlying space
	 */
	template<int N>
	class ITensor5
	{
	public:
		virtual ~ITensor5() = default;

		virtual int   NumContravar() const = 0;
		virtual int   NumCovar() const = 0;

		/**
		 * @brief Access tensor component.
		 * @param i First index
		 * @param j Second index
		 * @param k Third index
		 * @param l Fourth index
		 * @param m Fifth index
		 */
		virtual Real  operator()(int i, int j, int k, int l, int m) const = 0;
		virtual Real& operator()(int i, int j, int k, int l, int m) = 0;

		/**
		 * @brief Contract two indices, producing a rank-3 tensor.
		 * 
		 * Sums over the contracted index: result_{abc} = sum_s T_{...s...s...}
		 * Contraction reduces rank by 2, so rank-5 becomes rank-3.
		 * 
		 * @param ind1 First index to contract (0-4)
		 * @param ind2 Second index to contract (0-4), must differ from ind1
		 * @return Tensor3<N> - the contracted rank-3 result
		 * @throws TensorIndexError if indices out of range or equal
		 */
		virtual Tensor3<N> Contract(int ind1, int ind2) const = 0;

		/**
		 * @brief Apply tensor as a multilinear form to five vectors.
		 * 
		 * Computes T(v1, v2, v3, v4, v5) = sum_{i,j,k,l,m} T_ijklm * v1_i * v2_j * v3_k * v4_l * v5_m.
		 * This evaluates the tensor as a 5-linear form on the given vectors.
		 * 
		 * @param v1 First vector argument
		 * @param v2 Second vector argument
		 * @param v3 Third vector argument
		 * @param v4 Fourth vector argument
		 * @param v5 Fifth vector argument
		 * @return Scalar result of the 5-linear form evaluation
		 */
		virtual Real operator()(const VectorN<Real, N>& v1, const VectorN<Real, N>& v2,
		                         const VectorN<Real, N>& v3, const VectorN<Real, N>& v4,
		                         const VectorN<Real, N>& v5) const = 0;
	};
}
#endif


