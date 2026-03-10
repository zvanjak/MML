///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ITensorField.h                                                      ///
///  Description: Interface for tensor-valued fields (scalar, vector, tensor fields)  ///
///               Base classes for field operations on manifolds                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @file ITensorField.h
 * @brief Interfaces for position-dependent tensor fields.
 * 
 * Defines abstract interfaces for tensor fields - functions that assign a tensor
 * to each point in N-dimensional space. These are essential for differential
 * geometry and physics applications.
 * 
 * **Field Types:**
 * - ITensorField2<N>: Rank-2 tensor field, e.g., metric field g_ij(x)
 * - ITensorField3<N>: Rank-3 tensor field, e.g., Christoffel symbols Γ^i_jk(x)
 * - ITensorField4<N>: Rank-4 tensor field, e.g., Riemann curvature R^i_jkl(x)
 * - ITensorField5<N>: Rank-5 tensor field for specialized applications
 * 
 * Each field tracks its index variance (number of contravariant/covariant indices)
 * and provides both full tensor evaluation via operator() and component-wise
 * access via Component().
 * 
 * **Usage Examples:**
 * - Metric tensor fields for curved spacetime
 * - Electromagnetic field tensor
 * - Stress-energy tensor in general relativity
 * - Connection coefficients for parallel transport
 * 
 * @see ITensor.h for constant (position-independent) tensors
 * @see RiemannianManifold for applications
 */

#if !defined  MML_ITENSOR_FIELD_H
#define MML_ITENSOR_FIELD_H

#include "MMLBase.h"

#include "IFunction.h"

#include "base/Vector/Vector.h"
#include "base/Vector/VectorN.h"
#include "base/Tensor.h"

namespace MML
{
	/**
	 * @brief Interface for rank-2 tensor fields in N-dimensional space.
	 * 
	 * Maps each position x ∈ R^N to a rank-2 tensor T_ij(x) or T^ij(x).
	 * The most common example is the metric tensor field g_ij(x).
	 * 
	 * Implementations must provide the Component() method for individual
	 * component evaluation, which is often more efficient than computing
	 * the full tensor.
	 * 
	 * @tparam N Dimension of the base space
	 */
	template<int N>
	class ITensorField2 : public IFunction<Tensor2<N>, const VectorN<Real, N>& >
	{
		int _numContravar;  ///< Number of contravariant (upper) indices
		int _numCovar;      ///< Number of covariant (lower) indices
	public:
		/**
		 * @brief Construct a tensor field with specified index structure.
		 * @param numContra Number of contravariant indices (0, 1, or 2)
		 * @param numCo Number of covariant indices (2-numContra)
		 */
		ITensorField2(int numContra, int numCo) : _numContravar(numContra), _numCovar(numCo) {}

		/** @brief Get the number of contravariant indices. */
		int getNumContravar() const { return _numContravar; }
		
		/** @brief Get the number of covariant indices. */
		int getNumCovar()			const { return _numCovar; }

		/**
		 * @brief Evaluate a single tensor component at a position.
		 * 
		 * More efficient than operator() when only specific components are needed.
		 * 
		 * @param i First tensor index
		 * @param j Second tensor index
		 * @param pos Position in N-dimensional space
		 * @return The component T_ij(pos) or T^ij(pos)
		 */
		virtual Real Component(int i, int j, const VectorN<Real, N>& pos) const = 0;

		virtual ~ITensorField2() {}
	};

	/**
	 * @brief Interface for rank-3 tensor fields in N-dimensional space.
	 * 
	 * Maps each position to a rank-3 tensor. Primary application is the
	 * Christoffel symbol field Γ^i_jk(x) representing connection coefficients
	 * on a manifold with metric.
	 * 
	 * @tparam N Dimension of the base space
	 */
	template<int N>
	class ITensorField3 : public IFunction<Tensor3<N>, const VectorN<Real, N>& >
	{
	protected:
		int _numContravar;
		int _numCovar;
	public:
		ITensorField3(int numContra, int numCo) : _numContravar(numContra), _numCovar(numCo) {}
		ITensorField3() : _numContravar(0), _numCovar(3) {}

		int getNumContravar() const { return _numContravar; }
		int getNumCovar() const { return _numCovar; }

		/**
		 * @brief Evaluate a single tensor component at a position.
		 * @param i First tensor index
		 * @param j Second tensor index
		 * @param k Third tensor index
		 * @param pos Position in N-dimensional space
		 * @return The component value at pos
		 */
		virtual Real    Component(int i, int j, int k, const VectorN<Real, N>& pos) const = 0;

		virtual ~ITensorField3() {}
	};

	/**
	 * @brief Interface for rank-4 tensor fields in N-dimensional space.
	 * 
	 * Maps each position to a rank-4 tensor. Primary application is the
	 * Riemann curvature tensor field R^i_jkl(x) describing intrinsic
	 * curvature of a manifold.
	 * 
	 * @tparam N Dimension of the base space
	 */
	template<int N>
	class ITensorField4 : public IFunction<Tensor4<N>, const VectorN<Real, N>& >
	{
	protected:
		int _numContravar;
		int _numCovar;
	public:
		ITensorField4(int numContra, int numCo) : _numContravar(numContra), _numCovar(numCo) {}
		ITensorField4() : _numContravar(0), _numCovar(4) {}

		int getNumContravar() const { return _numContravar; }
		int getNumCovar() const { return _numCovar; }

		/**
		 * @brief Evaluate a single tensor component at a position.
		 * @param i First tensor index
		 * @param j Second tensor index
		 * @param k Third tensor index
		 * @param l Fourth tensor index
		 * @param pos Position in N-dimensional space
		 * @return The component value at pos
		 */
		virtual Real    Component(int i, int j, int k, int l, const VectorN<Real, N>& pos) const = 0;

		virtual ~ITensorField4() {}
	};

	/**
	 * @brief Interface for rank-5 tensor fields in N-dimensional space.
	 * 
	 * Maps each position to a rank-5 tensor for specialized applications
	 * such as covariant derivatives of the Riemann tensor or higher-order
	 * curvature invariants.
	 * 
	 * @tparam N Dimension of the base space
	 */
	template<int N>
	class ITensorField5 : public IFunction<Tensor5<N>, const VectorN<Real, N>& >
	{
	protected:
		int _numContravar;
		int _numCovar;
	public:
		ITensorField5(int numContra, int numCo) : _numContravar(numContra), _numCovar(numCo) {}
		ITensorField5() : _numContravar(0), _numCovar(5) {}

		int getNumContravar() const { return _numContravar; }
		int getNumCovar() const { return _numCovar; }

		/**
		 * @brief Evaluate a single tensor component at a position.
		 * @param i First tensor index
		 * @param j Second tensor index
		 * @param k Third tensor index
		 * @param l Fourth tensor index
		 * @param m Fifth tensor index
		 * @param pos Position in N-dimensional space
		 * @return The component value at pos
		 */
		virtual Real    Component(int i, int j, int k, int l, int m, const VectorN<Real, N>& pos) const = 0;
		virtual ~ITensorField5() {}
	};
}
#endif