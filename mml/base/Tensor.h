///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Tensor.h                                                            ///
///  Description: Tensor classes for multi-linear algebra                             ///
///               Rank-2, Rank-3, Rank-4, Rank-5 tensors with index notation          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @file Tensor.h
 * @brief Concrete tensor classes implementing ITensor interfaces.
 * 
 * Provides tensor classes of rank 2 through 5 in N-dimensional space, with full
 * support for covariant/contravariant index tracking and tensor operations.
 * 
 * **Index Variance Convention:**
 * - Constructor takes (nCovar, nContravar) - number of each index type
 * - By default, indices are ordered: covariant first, then contravariant
 * - The _isContravar[i] array tracks the variance of each index position
 * - Example: Tensor2<3>(1, 1) creates a mixed tensor T_i^j (1 down, 1 up)
 * 
 * **Tensor Notation Reference:**
 * - Contravariant (upper) indices: T^i transform with inverse Jacobian
 * - Covariant (lower) indices: T_i transform with Jacobian
 * - Contraction: Summing over one upper and one lower index (Einstein notation)
 * 
 * **Class Hierarchy:**
 * - Tensor2<N>: Rank-2 tensors (metric, stress, etc.)
 * - Tensor3<N>: Rank-3 tensors (Christoffel symbols, etc.)
 * - Tensor4<N>: Rank-4 tensors (Riemann curvature, elasticity, etc.)
 * - Tensor5<N>: Rank-5 tensors (specialized applications)
 * 
 * @see ITensor.h for abstract interfaces
 * @see ITensorField.h for position-dependent tensor fields
 */

#if !defined MML_TENSORS_H
#define MML_TENSORS_H

#include <cassert>

#include "MMLBase.h"

#include "interfaces/ITensor.h"

#include "base/Vector/VectorN.h"
#include "base/Matrix/MatrixNM.h"

// Standard headers - include what we use
#include <iomanip>
#include <initializer_list>
#include <iostream>

namespace MML
{
	///////////////////////////////////////////////////////////////////////////
	///                            Tensor2                                  ///
	///////////////////////////////////////////////////////////////////////////

	/// @brief Rank-2 tensor in N-dimensional space with covariant/contravariant tracking.
	/// @details Represents tensors with 2 indices, such as:
	///          - Metric tensor g_ij (2 covariant, 0 contravariant)
	///          - Inverse metric g^ij (0 covariant, 2 contravariant)
	///          - Mixed tensor T^i_j (1 covariant, 1 contravariant)
	///          - Stress tensor, electromagnetic field tensor, etc.
	///
	///          Components are stored in an N×N matrix. Index variance is tracked
	///          per-position via the _isContravar array, enabling proper tensor
	///          algebra including contraction (trace) operations.
	///
	/// @tparam N Dimension of the underlying vector space (e.g., 3 for 3D, 4 for spacetime)
	///
	/// @par Example Usage:
	/// @code
	///     // Create a metric tensor g_ij (fully covariant)
	///     Tensor2<3> metric(2, 0);  // 2 covariant, 0 contravariant
	///     metric(0,0) = 1.0; metric(1,1) = 1.0; metric(2,2) = 1.0;
	///
	///     // Create a mixed tensor T^i_j for contraction
	///     Tensor2<3> mixed(1, 1);  // 1 covariant, 1 contravariant
	///     Real trace = mixed.Contract();  // Sum of diagonal elements
	///
	///     // Evaluate as bilinear form
	///     VectorN<Real, 3> v1{1, 0, 0}, v2{0, 1, 0};
	///     Real result = metric(v1, v2);  // = g_ij * v1^i * v2^j
	/// @endcode
	template <int N>
	class Tensor2 : public ITensor2<N>
	{
		MatrixNM<Real, N, N> _coeff;   ///< Storage for tensor components
		int _numContravar = 0;         ///< Number of contravariant (upper) indices
		int _numCovar = 0;             ///< Number of covariant (lower) indices
		bool _isContravar[2];          ///< Per-index variance: true = contravariant, false = covariant
	public:

		///////////////////////////////////////////////////////////////////////
		///                        Constructors                             ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Constructs a zero-initialized rank-2 tensor with specified index variance.
		/// @param nCovar Number of covariant (lower) indices, must be non-negative
		/// @param nContraVar Number of contravariant (upper) indices, must be non-negative
		/// @throws TensorCovarContravarNumError if nCovar + nContraVar != 2 or either is negative
		/// @note The first nCovar indices are covariant, the remaining are contravariant.
		///       For example, (1, 1) creates T_i^j where index 0 is covariant and index 1 is contravariant.
		Tensor2(int nCovar, int nContraVar) : _numContravar(nContraVar), _numCovar(nCovar)
		{
			if ( _numContravar < 0 || _numCovar  < 0 || _numContravar + _numCovar != 2 )
				throw TensorCovarContravarNumError("Tensor2 ctor, wrong number of contravariant and covariant indices", nCovar, nContraVar);

			for (int i = 0; i < _numCovar; i++)
				_isContravar[i] = false;

			for (int i = _numCovar; i < _numCovar + _numContravar; i++)
				_isContravar[i] = true;
		}

		/// @brief Constructs a rank-2 tensor with specified values and index variance.
		/// @param nCovar Number of covariant (lower) indices
		/// @param nContraVar Number of contravariant (upper) indices
		/// @param values Initializer list of N×N values in row-major order
		/// @throws TensorCovarContravarNumError if nCovar + nContraVar != 2
		/// @note Values are filled row-by-row. Missing values are zero-initialized.
		Tensor2(int nCovar, int nContraVar, std::initializer_list<Real> values) : _numContravar(nContraVar), _numCovar(nCovar)
		{
			if ( _numContravar < 0 || _numCovar  < 0 || _numContravar + _numCovar != 2 )
				throw TensorCovarContravarNumError("Tensor2 ctor, wrong number of covariant and contravariant indices", nCovar, nContraVar);

			for (int i = 0; i < _numCovar; i++)
				_isContravar[i] = false;

			for (int i = _numCovar; i < _numCovar + _numContravar; i++)
				_isContravar[i] = true;

			if (values.size() > N * N)
				throw TensorCovarContravarNumError("Tensor2 ctor, initializer_list has more than N*N elements", static_cast<int>(values.size()), N * N);

			auto val = values.begin();
			for (size_t i = 0; i < N; ++i)
				for (size_t j = 0; j < N; ++j)
					if (val != values.end())
						_coeff[i][j] = *val++;
					else
						_coeff[i][j] = 0.0;
		}

		///////////////////////////////////////////////////////////////////////
		///                     Index Variance Query                        ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Returns the number of contravariant (upper) indices.
		int  NumContravar() const override { return _numContravar; }

		/// @brief Returns the number of covariant (lower) indices.
		int  NumCovar()     const override { return _numCovar; }

		/// @brief Returns the underlying coefficient matrix.
		MatrixNM<Real, N, N> GetMatrix() const {	return _coeff; }

		/// @brief Checks if index i is contravariant (upper index).
		/// @param i Index position (0 or 1)
		bool IsContravar(int i) const override { return _isContravar[i]; }

		/// @brief Checks if index i is covariant (lower index).
		/// @param i Index position (0 or 1)
		bool IsCovar(int i) const			{ return !_isContravar[i]; }

		/// @brief Sets the variance of index i.
		void setContravar(int i, bool val) { _isContravar[i] = val; }

		///////////////////////////////////////////////////////////////////////
		///                        Element Access                           ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Access tensor element (const) with assertion bounds checking.
		/// @param i First index (0 to N-1)
		/// @param j Second index (0 to N-1)
		/// @return The tensor component T(i,j)
		Real  operator()(int i, int j) const override { 
			assert(i >= 0 && i < N && "Tensor2: index i out of bounds");
			assert(j >= 0 && j < N && "Tensor2: index j out of bounds");
			return _coeff[i][j]; 
		}

		/// @brief Access tensor element (mutable) with assertion bounds checking.
		/// @param i First index (0 to N-1)
		/// @param j Second index (0 to N-1)
		/// @return Reference to the tensor component T(i,j)
		Real& operator()(int i, int j) override { 
			assert(i >= 0 && i < N && "Tensor2: index i out of bounds");
			assert(j >= 0 && j < N && "Tensor2: index j out of bounds");
			return _coeff[i][j]; 
		}
		
		/// @brief Checked element access (const) with exception on out-of-bounds.
		/// @param i First index (0 to N-1)
		/// @param j Second index (0 to N-1)
		/// @return The tensor component T(i,j)
		/// @throws TensorIndexError if any index is out of bounds
		Real at(int i, int j) const {
			if (i < 0 || i >= N || j < 0 || j >= N)
				throw TensorIndexError("Tensor2::at - index out of bounds");
			return _coeff[i][j];
		}

		/// @brief Checked element access (mutable) with exception on out-of-bounds.
		/// @param i First index (0 to N-1)
		/// @param j Second index (0 to N-1)
		/// @return Reference to the tensor component T(i,j)
		/// @throws TensorIndexError if any index is out of bounds
		Real& at(int i, int j) {
			if (i < 0 || i >= N || j < 0 || j >= N)
				throw TensorIndexError("Tensor2::at - index out of bounds");
			return _coeff[i][j];
		}

		///////////////////////////////////////////////////////////////////////
		///                      Arithmetic Operations                      ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Add two tensors of the same variance type.
		/// @throws TensorCovarContravarArithmeticError if index variance doesn't match
		Tensor2 operator+(const Tensor2& other) const
		{
			if (_numContravar != other._numContravar || _numCovar != other._numCovar)
				throw TensorCovarContravarArithmeticError("Tensor2 operator+, wrong number of contravariant and covariant indices", _numContravar, _numCovar, other._numContravar, other._numCovar);

			Tensor2 result(_numCovar, _numContravar);

			result._coeff = _coeff + other._coeff;

			return result;
		}

		/// @brief Subtract two tensors of the same variance type.
		/// @throws TensorCovarContravarArithmeticError if index variance doesn't match
		Tensor2 operator-(const Tensor2& other) const
		{
			if (_numContravar != other._numContravar || _numCovar != other._numCovar)
				throw TensorCovarContravarArithmeticError("Tensor2 operator-, wrong number of contravariant and covariant indices", _numContravar, _numCovar, other._numContravar, other._numCovar);

			Tensor2 result(_numCovar, _numContravar);

			result._coeff = _coeff - other._coeff;

			return result;
		}

		/// @brief Multiply tensor by a scalar.
		Tensor2 operator*(Real scalar) const
		{
			Tensor2 result(_numCovar, _numContravar);

			result._coeff = _coeff * scalar;

			return result;
		}

		/// @brief Divide tensor by a scalar.
		Tensor2 operator/(Real scalar) const
		{
			Tensor2 result(_numCovar, _numContravar);

			result._coeff = _coeff / scalar;

			return result;
		}

		/// @brief Scalar multiplication from the left (scalar * tensor).
		friend Tensor2 operator*(Real scalar, const Tensor2& b)
		{
			Tensor2 result(b.NumCovar(), b.NumContravar());

			result._coeff = b._coeff * scalar;

			return result;
		}

		///////////////////////////////////////////////////////////////////////
		///                       Tensor Operations                         ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Contract (trace) the tensor over its two indices.
		/// @details Computes the trace: sum_i T^i_i for a mixed tensor.
		///          This operation requires exactly one covariant and one contravariant index.
		/// @return The scalar trace value
		/// @throws TensorCovarContravarNumError if tensor is not mixed (1,1)
		/// @note Contraction is only meaningful for mixed tensors. For a (2,0) or (0,2) tensor,
		///       you would first need to raise/lower an index using a metric tensor.
		Real Contract() const override
		{
			if (_numContravar != 1 || _numCovar != 1)
				throw TensorCovarContravarNumError("Tensor2 Contract, wrong number of contravariant and covariant indices", _numContravar, _numCovar);

			Real result = 0.0;
			for (int i = 0; i < N; i++)
				result += _coeff[i][i];

			return result;
		}

		/// @brief Evaluate tensor as a bilinear form on two vectors.
		/// @details Computes T_ij * v1^i * v2^j (or appropriate variant based on variance).
		/// @param v1 First vector argument
		/// @param v2 Second vector argument
		/// @return The scalar result of the bilinear form evaluation
		Real operator()(const VectorN<Real, N>& v1, const VectorN<Real, N>& v2) const override
		{
			Real sum = 0.0;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					sum += _coeff[i][j] * v1[i] * v2[j];

			return sum;
		}

		///////////////////////////////////////////////////////////////////////
		///                            Output                               ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Print tensor components to a stream.
		void   Print(std::ostream& stream, int width, int precision) const
		{
			stream << std::fixed << "(N = " << N << ")" << std::endl;

			for (size_t i = 0; i < N; i++)
			{
				stream << "[ ";
				for (size_t j = 0; j < N; j++)
					stream << std::setw(width) << std::setprecision(precision) << _coeff[i][j] << ", ";
				stream << " ]" << std::endl;
			}
		}
		friend std::ostream& operator<<(std::ostream& stream, const Tensor2& a)
		{
			a.Print(stream, 15, 10);

			return stream;
		}
	};

	///////////////////////////////////////////////////////////////////////////
	///                            Tensor3                                  ///
	///////////////////////////////////////////////////////////////////////////

	/// @brief Rank-3 tensor in N-dimensional space with covariant/contravariant tracking.
	/// @details Represents tensors with 3 indices, commonly used for:
	///          - Christoffel symbols Γ^k_ij (connection coefficients)
	///          - Structure constants of Lie algebras
	///          - Completely antisymmetric tensors (Levi-Civita symbol)
	///          - Third-order material tensors
	///
	///          Index variance is tracked per-position. Contraction reduces rank by 2,
	///          returning a VectorN (rank-1 tensor).
	///
	/// @tparam N Dimension of the underlying vector space
	///
	/// @par Example Usage:
	/// @code
	///     // Create Christoffel-like tensor Γ^k_ij (1 up, 2 down)
	///     Tensor3<3> christoffel(2, 1);  // 2 covariant, 1 contravariant
	///     
	///     // Contract indices 1 and 2 (trace over the covariant indices)
	///     VectorN<Real, 3> contracted = christoffel.Contract(1, 2);
	/// @endcode
	template <int N>
	class Tensor3 : public ITensor3<N>
	{
		Real _coeff[N][N][N] = { 0 };  ///< Storage for N³ tensor components
		int _numContravar;             ///< Number of contravariant (upper) indices
		int _numCovar;                 ///< Number of covariant (lower) indices
		bool _isContravar[3];          ///< Per-index variance: true = contravariant, false = covariant
	public:

		///////////////////////////////////////////////////////////////////////
		///                        Constructors                             ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Constructs a zero-initialized rank-3 tensor with specified index variance.
		/// @param nCovar Number of covariant (lower) indices
		/// @param nContraVar Number of contravariant (upper) indices
		/// @throws TensorCovarContravarNumError if nCovar + nContraVar != 3
		Tensor3(int nCovar, int nContraVar) : _numContravar(nContraVar), _numCovar(nCovar)
		{
			if (_numContravar + _numCovar != 3)
				throw TensorCovarContravarNumError("Tensor3 ctor, wrong number of contravariant and covariant indices", nCovar, nContraVar);

			for (int i = 0; i < _numCovar; i++)
				_isContravar[i] = false;

			for (int i = _numCovar; i < _numCovar + _numContravar; i++)
				_isContravar[i] = true;
		}

		/// @brief Constructs a rank-3 tensor with specified values and index variance.
		/// @param nCovar Number of covariant (lower) indices
		/// @param nContraVar Number of contravariant (upper) indices  
		/// @param values Initializer list of N³ values in lexicographic order
		/// @throws TensorCovarContravarNumError if nCovar + nContraVar != 3
		Tensor3(int nCovar, int nContraVar, std::initializer_list<Real> values) : _numContravar(nContraVar), _numCovar(nCovar)
		{
			if (_numContravar + _numCovar != 3)
				throw TensorCovarContravarNumError("Tensor3 ctor, wrong number of contravariant and covariant indices", nCovar, nContraVar);
			
			for (int i = 0; i < _numCovar; i++)
				_isContravar[i] = false;
			
			for (int i = _numCovar; i < _numCovar + _numContravar; i++)
				_isContravar[i] = true;
			
			if (values.size() > N * N * N)
				throw TensorCovarContravarNumError("Tensor3 ctor, initializer_list has more than N*N*N elements", static_cast<int>(values.size()), N * N * N);

			auto val = values.begin();
			for (size_t i = 0; i < N; ++i)
				for (size_t j = 0; j < N; ++j)
					for (size_t k = 0; k < N; ++k)
						if (val != values.end())
							_coeff[i][j][k] = *val++;
						else
							_coeff[i][j][k] = 0.0;
		}

		///////////////////////////////////////////////////////////////////////
		///                     Index Variance Query                        ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Returns the number of contravariant (upper) indices.
		int   NumContravar() const override { return _numContravar; }

		/// @brief Returns the number of covariant (lower) indices.
		int   NumCovar()     const override { return _numCovar; }

		/// @brief Checks if index i is contravariant (upper index).
		bool IsContravar(int i) const override { return _isContravar[i]; }

		/// @brief Checks if index i is covariant (lower index).
		bool IsCovar(int i) const { return !_isContravar[i]; }

		/// @brief Sets the variance of index i.
		void setContravar(int i, bool val) { _isContravar[i] = val; }

		///////////////////////////////////////////////////////////////////////
		///                        Element Access                           ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Access tensor element (const) with assertion bounds checking.
		Real  operator()(int i, int j, int k) const override { 
			assert(i >= 0 && i < N && "Tensor3: index i out of bounds");
			assert(j >= 0 && j < N && "Tensor3: index j out of bounds");
			assert(k >= 0 && k < N && "Tensor3: index k out of bounds");
			return _coeff[i][j][k]; 
		}

		/// @brief Access tensor element (mutable) with assertion bounds checking.
		Real& operator()(int i, int j, int k) override { 
			assert(i >= 0 && i < N && "Tensor3: index i out of bounds");
			assert(j >= 0 && j < N && "Tensor3: index j out of bounds");
			assert(k >= 0 && k < N && "Tensor3: index k out of bounds");
			return _coeff[i][j][k]; 
		}
		
		/// @brief Checked element access (const) with exception on out-of-bounds.
		/// @throws TensorIndexError if any index is out of bounds
		Real at(int i, int j, int k) const {
			if (i < 0 || i >= N || j < 0 || j >= N || k < 0 || k >= N)
				throw TensorIndexError("Tensor3::at - index out of bounds");
			return _coeff[i][j][k];
		}

		/// @brief Checked element access (mutable) with exception on out-of-bounds.
		/// @throws TensorIndexError if any index is out of bounds
		Real& at(int i, int j, int k) {
			if (i < 0 || i >= N || j < 0 || j >= N || k < 0 || k >= N)
				throw TensorIndexError("Tensor3::at - index out of bounds");
			return _coeff[i][j][k];
		}

		///////////////////////////////////////////////////////////////////////
		///                       Tensor Operations                         ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Evaluate tensor as a trilinear form on three vectors.
		/// @details Computes T_ijk * v1^i * v2^j * v3^k (or appropriate variant).
		/// @return The scalar result of the trilinear form evaluation
		Real operator()(const VectorN<Real, N>& v1, const VectorN<Real, N>& v2, const VectorN<Real, N>& v3) const override
		{
			Real sum = 0.0;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
						sum += _coeff[i][j][k] * v1[i] * v2[j] * v3[k];

			return sum;
		}

		/// @brief Contract two indices, reducing rank from 3 to 1.
		/// @details Returns a VectorN by summing over the diagonal of the specified indices.
		///          For example, contracting indices 0 and 1 computes sum_s T(s,s,k) for each k.
		/// @param ind1 First index to contract (0-2)
		/// @param ind2 Second index to contract (0-2), must differ from ind1
		/// @return VectorN containing the contracted result (rank-1 tensor)
		/// @throws TensorIndexError if indices are out of range or equal
		/// @note The result inherits the variance of the remaining (non-contracted) index.
		VectorN<Real, N> Contract(int ind1, int ind2) const override
		{
			if (ind1 < 0 || ind1 > 2 || ind2 < 0 || ind2 > 2)
				throw TensorIndexError("Tensor3 Contract, index out of range [0,2]");
			if (ind1 == ind2)
				throw TensorIndexError("Tensor3 Contract, indices must be different");
			if (ind1 > ind2) std::swap(ind1, ind2);  // Normalize: ind1 < ind2

			VectorN<Real, N> result;

			// Find the remaining index (the one not being contracted)
			int remaining = -1;
			for (int i = 0; i < 3; i++)
				if (i != ind1 && i != ind2)
					remaining = i;

			// Contract: sum over the contracted index
			for (int a = 0; a < N; a++)  // result index
			{
				Real sum = 0.0;
				for (int s = 0; s < N; s++)  // contracted index
				{
					int idx[3];
					idx[ind1] = s;
					idx[ind2] = s;
					idx[remaining] = a;
					sum += _coeff[idx[0]][idx[1]][idx[2]];
				}
				result[a] = sum;
			}

			return result;
		}
	};

	///////////////////////////////////////////////////////////////////////////
	///                            Tensor4                                  ///
	///////////////////////////////////////////////////////////////////////////

	/// @brief Rank-4 tensor in N-dimensional space with covariant/contravariant tracking.
	/// @details Represents tensors with 4 indices, commonly used for:
	///          - Riemann curvature tensor R^a_bcd
	///          - Elasticity tensor C_ijkl
	///          - Electromagnetic constitutive tensors
	///          - Fourth-order material property tensors
	///
	///          Index variance is tracked per-position. Contraction reduces rank by 2,
	///          returning a Tensor2 with properly propagated index variance.
	///
	/// @tparam N Dimension of the underlying vector space
	///
	/// @par Example Usage:
	/// @code
	///     // Create Riemann tensor R^a_bcd (1 up, 3 down)
	///     Tensor4<4> riemann(3, 1);  // 3 covariant, 1 contravariant
	///     
	///     // Contract indices to get Ricci tensor R_bd = R^a_bad
	///     Tensor2<4> ricci = riemann.Contract(0, 2);  // Contract indices 0 and 2
	/// @endcode
	template <int N>
	class Tensor4 : public ITensor4<N>
	{
		Real _coeff[N][N][N][N] = { 0 };  ///< Storage for N⁴ tensor components
		int _numContravar;                 ///< Number of contravariant (upper) indices
		int _numCovar;                     ///< Number of covariant (lower) indices
		bool _isContravar[4];              ///< Per-index variance: true = contravariant, false = covariant
	public:

		///////////////////////////////////////////////////////////////////////
		///                        Constructors                             ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Constructs a zero-initialized rank-4 tensor with specified index variance.
		/// @param nCovar Number of covariant (lower) indices
		/// @param nContraVar Number of contravariant (upper) indices
		/// @throws TensorCovarContravarNumError if nCovar + nContraVar != 4
		Tensor4(int nCovar, int nContraVar) : _numContravar(nContraVar), _numCovar(nCovar) 
		{
			if (_numContravar + _numCovar != 4)
				throw TensorCovarContravarNumError("Tensor4 ctor, wrong number of contravariant and covariant indices", nCovar, nContraVar);

			for (int i = 0; i < _numCovar; i++)
				_isContravar[i] = false;

			for (int i = _numCovar; i < _numCovar + _numContravar; i++)
				_isContravar[i] = true;
		}

		///////////////////////////////////////////////////////////////////////
		///                     Index Variance Query                        ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Returns the number of contravariant (upper) indices.
		int   NumContravar() const override { return _numContravar; }

		/// @brief Returns the number of covariant (lower) indices.
		int   NumCovar()     const override { return _numCovar; }

		/// @brief Checks if index i is contravariant (upper index).
		bool IsContravar(int i) const override { return _isContravar[i]; }

		/// @brief Checks if index i is covariant (lower index).
		bool IsCovar(int i) const { return !_isContravar[i]; }

		/// @brief Sets the variance of index i.
		void setContravar(int i, bool val) { _isContravar[i] = val; }

		///////////////////////////////////////////////////////////////////////
		///                        Element Access                           ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Access tensor element (const) with assertion bounds checking.
		Real  operator()(int i, int j, int k, int l) const override { 
			assert(i >= 0 && i < N && "Tensor4: index i out of bounds");
			assert(j >= 0 && j < N && "Tensor4: index j out of bounds");
			assert(k >= 0 && k < N && "Tensor4: index k out of bounds");
			assert(l >= 0 && l < N && "Tensor4: index l out of bounds");
			return _coeff[i][j][k][l]; 
		}

		/// @brief Access tensor element (mutable) with assertion bounds checking.
		Real& operator()(int i, int j, int k, int l) override { 
			assert(i >= 0 && i < N && "Tensor4: index i out of bounds");
			assert(j >= 0 && j < N && "Tensor4: index j out of bounds");
			assert(k >= 0 && k < N && "Tensor4: index k out of bounds");
			assert(l >= 0 && l < N && "Tensor4: index l out of bounds");
			return _coeff[i][j][k][l]; 
		}
		
		/// @brief Checked element access (const) with exception on out-of-bounds.
		/// @throws TensorIndexError if any index is out of bounds
		Real at(int i, int j, int k, int l) const {
			if (i < 0 || i >= N || j < 0 || j >= N || k < 0 || k >= N || l < 0 || l >= N)
				throw TensorIndexError("Tensor4::at - index out of bounds");
			return _coeff[i][j][k][l];
		}

		/// @brief Checked element access (mutable) with exception on out-of-bounds.
		/// @throws TensorIndexError if any index is out of bounds
		Real& at(int i, int j, int k, int l) {
			if (i < 0 || i >= N || j < 0 || j >= N || k < 0 || k >= N || l < 0 || l >= N)
				throw TensorIndexError("Tensor4::at - index out of bounds");
			return _coeff[i][j][k][l];
		}

		///////////////////////////////////////////////////////////////////////
		///                       Tensor Operations                         ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Evaluate tensor as a quadrilinear form on four vectors.
		/// @details Computes T_ijkl * v1^i * v2^j * v3^k * v4^l (or appropriate variant).
		/// @return The scalar result of the quadrilinear form evaluation
		Real operator()(const VectorN<Real, N>& v1, const VectorN<Real, N>& v2, const VectorN<Real, N>& v3, const VectorN<Real, N>& v4) const override
		{
			Real sum = 0.0;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
						for (int l = 0; l < N; l++)
							sum += _coeff[i][j][k][l] * v1[i] * v2[j] * v3[k] * v4[l];

			return sum;
		}

		/// @brief Contract two indices, reducing rank from 4 to 2.
		/// @details Returns a Tensor2 by summing over the diagonal of the specified indices.
		///          The result tensor correctly inherits the variance of non-contracted indices.
		///
		///          For example, given tensor T^i_j^k_l and contracting indices 0 and 1:
		///          - Result has indices from positions 2 and 3 (originally k and l)
		///          - Result's _isContravar is copied from original positions [2] and [3]
		///          - So result would have variance pattern matching ^k_l
		///
		/// @param ind1 First index to contract (0-3)
		/// @param ind2 Second index to contract (0-3), must differ from ind1
		/// @return Tensor2 containing the contracted result with proper index variance
		/// @throws TensorIndexError if indices are out of range or equal
		/// @note Index variance is explicitly propagated from the original tensor,
		///       not reconstructed from covariant/contravariant counts.
		Tensor2<N> Contract(int ind1, int ind2) const override
		{
			if (ind1 < 0 || ind1 > 3 || ind2 < 0 || ind2 > 3)
				throw TensorIndexError("Tensor4 Contract, index out of range [0,3]");
			if (ind1 == ind2)
				throw TensorIndexError("Tensor4 Contract, indices must be different");
			if (ind1 > ind2) std::swap(ind1, ind2);  // Normalize: ind1 < ind2

			// Determine covariant/contravariant counts for result
			int newCovar = _numCovar;
			int newContravar = _numContravar;
			if (_isContravar[ind1]) newContravar--; else newCovar--;
			if (_isContravar[ind2]) newContravar--; else newCovar--;

			Tensor2<N> result(newCovar, newContravar);

			// Build mapping: which original indices become result indices
			int map[2];  // map[result_idx] = original_idx
			int m_idx = 0;
			for (int i = 0; i < 4; i++)
				if (i != ind1 && i != ind2)
					map[m_idx++] = i;

			// Properly propagate index variance from original tensor
			for (int i = 0; i < 2; i++)
				result.setContravar(i, _isContravar[map[i]]);

			// Contract: sum over the contracted index
			for (int a = 0; a < N; a++)        // result index 0
				for (int b = 0; b < N; b++)    // result index 1
				{
					Real sum = 0.0;
					for (int s = 0; s < N; s++)  // contracted index
					{
						int idx[4];
						idx[ind1] = s;
						idx[ind2] = s;
						idx[map[0]] = a;
						idx[map[1]] = b;
						sum += _coeff[idx[0]][idx[1]][idx[2]][idx[3]];
					}
					result(a, b) = sum;
				}

			return result;
		}
	};

	///////////////////////////////////////////////////////////////////////////
	///                            Tensor5                                  ///
	///////////////////////////////////////////////////////////////////////////

	/// @brief Rank-5 tensor in N-dimensional space with covariant/contravariant tracking.
	/// @details Represents tensors with 5 indices, used in specialized applications:
	///          - Higher-order material property tensors
	///          - Derivatives of fourth-order tensors
	///          - Theoretical physics applications
	///
	///          Index variance is tracked per-position. Contraction reduces rank by 2,
	///          returning a Tensor3 with properly propagated index variance.
	///
	/// @tparam N Dimension of the underlying vector space
	///
	/// @par Example Usage:
	/// @code
	///     // Create a rank-5 tensor with 3 covariant, 2 contravariant indices
	///     Tensor5<3> t(3, 2);  // T_ijk^lm
	///     
	///     // Contract to get a Tensor3
	///     Tensor3<3> contracted = t.Contract(0, 3);  // Contract first covariant with first contravariant
	/// @endcode
	template <int N>
	class Tensor5 : public ITensor5<N>
	{
		Real _coeff[N][N][N][N][N] = { 0 };  ///< Storage for N⁵ tensor components
		int _numContravar;                    ///< Number of contravariant (upper) indices
		int _numCovar;                        ///< Number of covariant (lower) indices
		bool _isContravar[5];                 ///< Per-index variance: true = contravariant, false = covariant
	public:

		///////////////////////////////////////////////////////////////////////
		///                        Constructors                             ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Constructs a zero-initialized rank-5 tensor with specified index variance.
		/// @param nCovar Number of covariant (lower) indices
		/// @param nContraVar Number of contravariant (upper) indices
		/// @throws TensorCovarContravarNumError if nCovar + nContraVar != 5
		Tensor5(int nCovar, int nContraVar) : _numContravar(nContraVar), _numCovar(nCovar) 
		{
			if (_numContravar + _numCovar != 5)
				throw TensorCovarContravarNumError("Tensor5 ctor, wrong number of contravariant and covariant indices", nCovar, nContraVar);

			for (int i = 0; i < _numCovar; i++)
				_isContravar[i] = false;

			for (int i = _numCovar; i < _numCovar + _numContravar; i++)
				_isContravar[i] = true;

		}

		///////////////////////////////////////////////////////////////////////
		///                     Index Variance Query                        ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Returns the number of contravariant (upper) indices.
		int   NumContravar() const override { return _numContravar; }

		/// @brief Returns the number of covariant (lower) indices.
		int   NumCovar()     const override { return _numCovar; }

		/// @brief Checks if index i is contravariant (upper index).
		bool IsContravar(int i) const { return _isContravar[i]; }

		/// @brief Checks if index i is covariant (lower index).
		bool IsCovar(int i) const { return !_isContravar[i]; }

		/// @brief Sets the variance of index i.
		void setContravar(int i, bool val) { _isContravar[i] = val; }

		///////////////////////////////////////////////////////////////////////
		///                        Element Access                           ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Access tensor element (const) with assertion bounds checking.
		Real  operator()(int i, int j, int k, int l, int m) const override { 
			assert(i >= 0 && i < N && "Tensor5: index i out of bounds");
			assert(j >= 0 && j < N && "Tensor5: index j out of bounds");
			assert(k >= 0 && k < N && "Tensor5: index k out of bounds");
			assert(l >= 0 && l < N && "Tensor5: index l out of bounds");
			assert(m >= 0 && m < N && "Tensor5: index m out of bounds");
			return _coeff[i][j][k][l][m]; 
		}

		/// @brief Access tensor element (mutable) with assertion bounds checking.
		Real& operator()(int i, int j, int k, int l, int m) override { 
			assert(i >= 0 && i < N && "Tensor5: index i out of bounds");
			assert(j >= 0 && j < N && "Tensor5: index j out of bounds");
			assert(k >= 0 && k < N && "Tensor5: index k out of bounds");
			assert(l >= 0 && l < N && "Tensor5: index l out of bounds");
			assert(m >= 0 && m < N && "Tensor5: index m out of bounds");
			return _coeff[i][j][k][l][m]; 
		}
		
		/// @brief Checked element access (const) with exception on out-of-bounds.
		/// @throws TensorIndexError if any index is out of bounds
		Real at(int i, int j, int k, int l, int m) const {
			if (i < 0 || i >= N || j < 0 || j >= N || k < 0 || k >= N || l < 0 || l >= N || m < 0 || m >= N)
				throw TensorIndexError("Tensor5::at - index out of bounds");
			return _coeff[i][j][k][l][m];
		}

		/// @brief Checked element access (mutable) with exception on out-of-bounds.
		/// @throws TensorIndexError if any index is out of bounds
		Real& at(int i, int j, int k, int l, int m) {
			if (i < 0 || i >= N || j < 0 || j >= N || k < 0 || k >= N || l < 0 || l >= N || m < 0 || m >= N)
				throw TensorIndexError("Tensor5::at - index out of bounds");
			return _coeff[i][j][k][l][m];
		}

		///////////////////////////////////////////////////////////////////////
		///                       Tensor Operations                         ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Evaluate tensor as a quintilinear form on five vectors.
		/// @details Computes T_ijklm * v1^i * v2^j * v3^k * v4^l * v5^m (or appropriate variant).
		/// @return The scalar result of the quintilinear form evaluation
		Real operator()(const VectorN<Real, N>& v1, const VectorN<Real, N>& v2, 
		                const VectorN<Real, N>& v3, const VectorN<Real, N>& v4, 
		                const VectorN<Real, N>& v5) const override
		{
			Real sum = 0.0;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
						for (int l = 0; l < N; l++)
							for (int m = 0; m < N; m++)
								sum += _coeff[i][j][k][l][m] * v1[i] * v2[j] * v3[k] * v4[l] * v5[m];
			return sum;
		}

		/// @brief Contract two indices, reducing rank from 5 to 3.
		/// @details Returns a Tensor3 by summing over the diagonal of the specified indices.
		///          The result tensor correctly inherits the variance of non-contracted indices.
		///
		///          For example, given tensor T_ij^k_l^m and contracting indices 1 and 2:
		///          - Result has indices from positions 0, 3, 4 (originally i, l, m)
		///          - Result's _isContravar is copied from original positions [0], [3], [4]
		///          - So result would have variance pattern matching _i_l^m
		///
		/// @param ind1 First index to contract (0-4)
		/// @param ind2 Second index to contract (0-4), must differ from ind1
		/// @return Tensor3 containing the contracted result with proper index variance
		/// @throws TensorIndexError if indices are out of range or equal
		/// @note Index variance is explicitly propagated from the original tensor,
		///       not reconstructed from covariant/contravariant counts.
		Tensor3<N> Contract(int ind1, int ind2) const override
		{
			if (ind1 < 0 || ind1 > 4 || ind2 < 0 || ind2 > 4)
				throw TensorIndexError("Tensor5 Contract, index out of range [0,4]");
			if (ind1 == ind2)
				throw TensorIndexError("Tensor5 Contract, indices must be different");
			if (ind1 > ind2) std::swap(ind1, ind2);  // Normalize: ind1 < ind2

			// Determine covariant/contravariant counts for result
			int newCovar = _numCovar;
			int newContravar = _numContravar;
			if (_isContravar[ind1]) newContravar--; else newCovar--;
			if (_isContravar[ind2]) newContravar--; else newCovar--;

			Tensor3<N> result(newCovar, newContravar);

			// Build mapping: which original indices become result indices
			int map[3];  // map[result_idx] = original_idx
			int m_idx = 0;
			for (int i = 0; i < 5; i++)
				if (i != ind1 && i != ind2)
					map[m_idx++] = i;

			// Properly propagate index variance from original tensor
			for (int i = 0; i < 3; i++)
				result.setContravar(i, _isContravar[map[i]]);

			// Contract: sum over the contracted index
			for (int a = 0; a < N; a++)        // result index 0
				for (int b = 0; b < N; b++)    // result index 1
					for (int c = 0; c < N; c++)// result index 2
					{
						Real sum = 0.0;
						for (int s = 0; s < N; s++)  // contracted index
						{
							// Build the 5 indices for _coeff access
							int idx[5];
							idx[ind1] = s;
							idx[ind2] = s;
							idx[map[0]] = a;
							idx[map[1]] = b;
							idx[map[2]] = c;
							sum += _coeff[idx[0]][idx[1]][idx[2]][idx[3]][idx[4]];
						}
						result(a, b, c) = sum;
					}

			return result;
		}
	};
}
#endif