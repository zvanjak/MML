///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Tensor5.h                                                          ///
///  Description: Rank-5 concrete tensor and rank-specific operations                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_TENSOR5_H
#define MML_TENSOR5_H

#include <cassert>
#include <utility>

#include <mml/MMLBase.h>
#include <mml/interfaces/ITensor.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Tensor/MatrixNDim.h>
#include <mml/base/Tensor/Tensor3.h>

namespace MML
{
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
		MatrixNDim<Real, N, 5> _coeff;        ///< Storage for N⁵ tensor components
		int _numContravar;                    ///< Number of contravariant (upper) indices
		int _numCovar;                        ///< Number of covariant (lower) indices
		bool _isContravar[5];                 ///< Per-index variance: true = contravariant, false = covariant

		/// @brief Zero tensor with this tensor's exact per-index variance arrangement.
		Tensor5 sameShapeZero() const
		{
			Tensor5 r(_numCovar, _numContravar);
			for (int i = 0; i < 5; i++) r.setContravar(i, _isContravar[i]);
			return r;
		}

		/// @brief Throw if `other` does not share this tensor's per-index variance.
		void checkSameVariance(const Tensor5& other, const char* where) const
		{
			for (int i = 0; i < 5; i++)
				if (_isContravar[i] != other._isContravar[i])
					throw TensorCovarContravarArithmeticError(where, _numContravar, _numCovar, other._numContravar, other._numCovar);
		}
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
		bool IsContravar(int i) const override { return _isContravar[i]; }

		/// @brief Checks if index i is covariant (lower index).
		bool IsCovar(int i) const { return !_isContravar[i]; }

		/// @brief Sets the variance of index i, keeping covariant/contravariant counts in sync.
		void setContravar(int i, bool val)
		{
			_isContravar[i] = val;
			_numContravar = 0;
			for (bool c : _isContravar)
				if (c) _numContravar++;
			_numCovar = static_cast<int>(sizeof(_isContravar) / sizeof(_isContravar[0])) - _numContravar;
		}

		///////////////////////////////////////////////////////////////////////
		///                        Element Access                           ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Unchecked element access (const) — fast, assertion-only bounds checking.
		/// @warning Assertions are disabled in Release builds. For runtime safety, use at().
		Real  operator()(int i, int j, int k, int l, int m) const override { 
			assert(i >= 0 && i < N && "Tensor5: index i out of bounds");
			assert(j >= 0 && j < N && "Tensor5: index j out of bounds");
			assert(k >= 0 && k < N && "Tensor5: index k out of bounds");
			assert(l >= 0 && l < N && "Tensor5: index l out of bounds");
			assert(m >= 0 && m < N && "Tensor5: index m out of bounds");
			return _coeff(i, j, k, l, m); 
		}

		/// @brief Unchecked element access (mutable) — fast, assertion-only bounds checking.
		/// @warning Assertions are disabled in Release builds. For runtime safety, use at().
		Real& operator()(int i, int j, int k, int l, int m) override { 
			assert(i >= 0 && i < N && "Tensor5: index i out of bounds");
			assert(j >= 0 && j < N && "Tensor5: index j out of bounds");
			assert(k >= 0 && k < N && "Tensor5: index k out of bounds");
			assert(l >= 0 && l < N && "Tensor5: index l out of bounds");
			assert(m >= 0 && m < N && "Tensor5: index m out of bounds");
			return _coeff(i, j, k, l, m); 
		}
		
		/// @brief Checked element access (const) — safe, throws on out-of-bounds.
		/// @details Use this instead of operator() when runtime bounds safety is needed.
		/// @throws TensorIndexError if any index is out of bounds
		Real at(int i, int j, int k, int l, int m) const {
			if (i < 0 || i >= N || j < 0 || j >= N || k < 0 || k >= N || l < 0 || l >= N || m < 0 || m >= N)
				throw TensorIndexError("Tensor5::at - index out of bounds");
			return _coeff(i, j, k, l, m);
		}

		/// @brief Checked element access (mutable) — safe, throws on out-of-bounds.
		/// @details Use this instead of operator() when runtime bounds safety is needed.
		/// @throws TensorIndexError if any index is out of bounds
		Real& at(int i, int j, int k, int l, int m) {
			if (i < 0 || i >= N || j < 0 || j >= N || k < 0 || k >= N || l < 0 || l >= N || m < 0 || m >= N)
				throw TensorIndexError("Tensor5::at - index out of bounds");
			return _coeff(i, j, k, l, m);
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
								sum += _coeff(i, j, k, l, m) * v1[i] * v2[j] * v3[k] * v4[l] * v5[m];
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
							sum += _coeff(idx[0], idx[1], idx[2], idx[3], idx[4]);
						}
						result(a, b, c) = sum;
					}

			return result;
		}

		///////////////////////////////////////////////////////////////////////
		///                      Arithmetic Operations                      ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Add two tensors sharing the same per-index variance.
		/// @throws TensorCovarContravarArithmeticError if index variance doesn't match
		Tensor5 operator+(const Tensor5& other) const
		{
			checkSameVariance(other, "Tensor5 operator+, index variance mismatch");
			Tensor5 result = sameShapeZero();
			result._coeff = _coeff + other._coeff;
			return result;
		}

		/// @brief Subtract two tensors sharing the same per-index variance.
		/// @throws TensorCovarContravarArithmeticError if index variance doesn't match
		Tensor5 operator-(const Tensor5& other) const
		{
			checkSameVariance(other, "Tensor5 operator-, index variance mismatch");
			Tensor5 result = sameShapeZero();
			result._coeff = _coeff - other._coeff;
			return result;
		}

		/// @brief Unary negation.
		Tensor5 operator-() const
		{
			Tensor5 result = sameShapeZero();
			result._coeff = -_coeff;
			return result;
		}

		/// @brief Multiply tensor by a scalar.
		Tensor5 operator*(Real scalar) const
		{
			Tensor5 result = sameShapeZero();
			result._coeff = _coeff * scalar;
			return result;
		}

		/// @brief Divide tensor by a scalar.
		Tensor5 operator/(Real scalar) const
		{
			Tensor5 result = sameShapeZero();
			result._coeff = _coeff / scalar;
			return result;
		}

		/// @brief Scalar multiplication from the left (scalar * tensor).
		friend Tensor5 operator*(Real scalar, const Tensor5& t)
		{
			Tensor5 result = t.sameShapeZero();
			result._coeff = t._coeff * scalar;
			return result;
		}

		Tensor5& operator+=(const Tensor5& other) { checkSameVariance(other, "Tensor5 operator+=, index variance mismatch"); _coeff += other._coeff; return *this; }
		Tensor5& operator-=(const Tensor5& other) { checkSameVariance(other, "Tensor5 operator-=, index variance mismatch"); _coeff -= other._coeff; return *this; }
		Tensor5& operator*=(Real scalar) { _coeff *= scalar; return *this; }
		Tensor5& operator/=(Real scalar) { _coeff /= scalar; return *this; }

		/// @brief Exact equality: same per-index variance and identical components.
		bool operator==(const Tensor5& other) const
		{
			for (int i = 0; i < 5; i++)
				if (_isContravar[i] != other._isContravar[i])
					return false;
			return _coeff == other._coeff;
		}
		bool operator!=(const Tensor5& other) const { return !(*this == other); }
	};
}
#endif
