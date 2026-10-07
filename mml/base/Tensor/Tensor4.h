///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Tensor4.h                                                          ///
///  Description: Rank-4 concrete tensor and rank-specific operations                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_TENSOR4_H
#define MML_TENSOR4_H

#include <array>
#include <cassert>
#include <cmath>
#include <utility>

#include <mml/MMLBase.h>
#include <mml/interfaces/ITensor.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/BaseUtils/SymbolUtils.h>
#include <mml/base/Tensor/TensorCommon.h>
#include <mml/base/Tensor/MatrixNDim.h>
#include <mml/base/Tensor/Tensor2.h>

namespace MML
{
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
		MatrixNDim<Real, N, 4> _coeff;     ///< Storage for N⁴ tensor components
		int _numContravar;                 ///< Number of contravariant (upper) indices
		int _numCovar;                     ///< Number of covariant (lower) indices
		bool _isContravar[4];              ///< Per-index variance: true = contravariant, false = covariant

		/// @brief Zero tensor with this tensor's exact per-index variance arrangement.
		Tensor4 sameShapeZero() const
		{
			Tensor4 r(_numCovar, _numContravar);
			for (int i = 0; i < 4; i++) r.setContravar(i, _isContravar[i]);
			return r;
		}

		/// @brief Throw if `other` does not share this tensor's per-index variance.
		void checkSameVariance(const Tensor4& other, const char* where) const
		{
			for (int i = 0; i < 4; i++)
				if (_isContravar[i] != other._isContravar[i])
					throw TensorCovarContravarArithmeticError(where, _numContravar, _numCovar, other._numContravar, other._numCovar);
		}
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
		Real  operator()(int i, int j, int k, int l) const override { 
			assert(i >= 0 && i < N && "Tensor4: index i out of bounds");
			assert(j >= 0 && j < N && "Tensor4: index j out of bounds");
			assert(k >= 0 && k < N && "Tensor4: index k out of bounds");
			assert(l >= 0 && l < N && "Tensor4: index l out of bounds");
			return _coeff(i, j, k, l); 
		}

		/// @brief Unchecked element access (mutable) — fast, assertion-only bounds checking.
		/// @warning Assertions are disabled in Release builds. For runtime safety, use at().
		Real& operator()(int i, int j, int k, int l) override { 
			assert(i >= 0 && i < N && "Tensor4: index i out of bounds");
			assert(j >= 0 && j < N && "Tensor4: index j out of bounds");
			assert(k >= 0 && k < N && "Tensor4: index k out of bounds");
			assert(l >= 0 && l < N && "Tensor4: index l out of bounds");
			return _coeff(i, j, k, l); 
		}
		
		/// @brief Checked element access (const) — safe, throws on out-of-bounds.
		/// @details Use this instead of operator() when runtime bounds safety is needed.
		/// @throws TensorIndexError if any index is out of bounds
		Real at(int i, int j, int k, int l) const {
			if (i < 0 || i >= N || j < 0 || j >= N || k < 0 || k >= N || l < 0 || l >= N)
				throw TensorIndexError("Tensor4::at - index out of bounds");
			return _coeff(i, j, k, l);
		}

		/// @brief Checked element access (mutable) — safe, throws on out-of-bounds.
		/// @details Use this instead of operator() when runtime bounds safety is needed.
		/// @throws TensorIndexError if any index is out of bounds
		Real& at(int i, int j, int k, int l) {
			if (i < 0 || i >= N || j < 0 || j >= N || k < 0 || k >= N || l < 0 || l >= N)
				throw TensorIndexError("Tensor4::at - index out of bounds");
			return _coeff(i, j, k, l);
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
							sum += _coeff(i, j, k, l) * v1[i] * v2[j] * v3[k] * v4[l];

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
						sum += _coeff(idx[0], idx[1], idx[2], idx[3]);
					}
					result(a, b) = sum;
				}

			return result;
		}

		///////////////////////////////////////////////////////////////////////
		///                      Arithmetic Operations                      ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Add two tensors sharing the same per-index variance.
		/// @throws TensorCovarContravarArithmeticError if index variance doesn't match
		Tensor4 operator+(const Tensor4& other) const
		{
			checkSameVariance(other, "Tensor4 operator+, index variance mismatch");
			Tensor4 result = sameShapeZero();
			result._coeff = _coeff + other._coeff;
			return result;
		}

		/// @brief Subtract two tensors sharing the same per-index variance.
		/// @throws TensorCovarContravarArithmeticError if index variance doesn't match
		Tensor4 operator-(const Tensor4& other) const
		{
			checkSameVariance(other, "Tensor4 operator-, index variance mismatch");
			Tensor4 result = sameShapeZero();
			result._coeff = _coeff - other._coeff;
			return result;
		}

		/// @brief Unary negation.
		Tensor4 operator-() const
		{
			Tensor4 result = sameShapeZero();
			result._coeff = -_coeff;
			return result;
		}

		/// @brief Multiply tensor by a scalar.
		Tensor4 operator*(Real scalar) const
		{
			Tensor4 result = sameShapeZero();
			result._coeff = _coeff * scalar;
			return result;
		}

		/// @brief Divide tensor by a scalar.
		Tensor4 operator/(Real scalar) const
		{
			Tensor4 result = sameShapeZero();
			result._coeff = _coeff / scalar;
			return result;
		}

		/// @brief Scalar multiplication from the left (scalar * tensor).
		friend Tensor4 operator*(Real scalar, const Tensor4& t)
		{
			Tensor4 result = t.sameShapeZero();
			result._coeff = t._coeff * scalar;
			return result;
		}

		Tensor4& operator+=(const Tensor4& other) { checkSameVariance(other, "Tensor4 operator+=, index variance mismatch"); _coeff += other._coeff; return *this; }
		Tensor4& operator-=(const Tensor4& other) { checkSameVariance(other, "Tensor4 operator-=, index variance mismatch"); _coeff -= other._coeff; return *this; }
		Tensor4& operator*=(Real scalar) { _coeff *= scalar; return *this; }
		Tensor4& operator/=(Real scalar) { _coeff /= scalar; return *this; }

		/// @brief Exact equality: same per-index variance and identical components.
		bool operator==(const Tensor4& other) const
		{
			for (int i = 0; i < 4; i++)
				if (_isContravar[i] != other._isContravar[i])
					return false;
			return _coeff == other._coeff;
		}
		bool operator!=(const Tensor4& other) const { return !(*this == other); }
	};

	/// @brief Permute rank-4 tensor indices.
	/// @details permutation[resultSlot] gives the source slot that becomes that result slot.
	template<int N>
	Tensor4<N> PermuteIndices(const Tensor4<N>& tensor, const std::array<int, 4>& permutation)
	{
		detail::ValidateTensorPermutation<4>(permutation);

		bool resultContravar[4] = {
			tensor.IsContravar(permutation[0]),
			tensor.IsContravar(permutation[1]),
			tensor.IsContravar(permutation[2]),
			tensor.IsContravar(permutation[3])
		};
		int numContravar = 0;
		for (int i = 0; i < 4; i++)
			if (resultContravar[i])
				numContravar++;

		Tensor4<N> result(4 - numContravar, numContravar);
		for (int i = 0; i < 4; i++)
			result.setContravar(i, resultContravar[i]);

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				for (int k = 0; k < N; k++)
					for (int l = 0; l < N; l++)
					{
						int sourceIndex[4];
						sourceIndex[permutation[0]] = i;
						sourceIndex[permutation[1]] = j;
						sourceIndex[permutation[2]] = k;
						sourceIndex[permutation[3]] = l;
						result(i, j, k, l) = tensor(sourceIndex[0], sourceIndex[1], sourceIndex[2], sourceIndex[3]);
					}

		return result;
	}

	/// @brief Covariant Levi-Civita tensor in 4D: E_abcd = sqrt(|det(g)|) epsilon_abcd.
	template<int N>
	Tensor4<N> LeviCivitaTensor4(const MatrixNM<Real, N, N>& covariantMetric)
	{
		static_assert(N == 4, "LeviCivitaTensor4 is available only for N=4");

		Real weight = std::sqrt(std::abs(detail::DeterminantForLeviCivitaTensor<N>(covariantMetric)));
		Tensor4<N> result(4, 0);

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				for (int k = 0; k < N; k++)
					for (int l = 0; l < N; l++)
						result(i, j, k, l) = weight * Utils::LeviCivita(i + 1, j + 1, k + 1, l + 1);

		return result;
	}

	/// @brief Outer product of two rank-2 tensors, producing a rank-4 tensor.
	/// @details Result indices are ordered as the left tensor's two indices followed
	///          by the right tensor's two indices, with variance preserved per slot.
	template<int N>
	Tensor4<N> OuterProduct(const Tensor2<N>& left, const Tensor2<N>& right)
	{
		bool resultContravar[4] = {
			left.IsContravar(0),
			left.IsContravar(1),
			right.IsContravar(0),
			right.IsContravar(1)
		};
		int numContravar = 0;
		for (int i = 0; i < 4; i++)
			if (resultContravar[i])
				numContravar++;

		Tensor4<N> result(4 - numContravar, numContravar);
		for (int i = 0; i < 4; i++)
			result.setContravar(i, resultContravar[i]);

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				for (int k = 0; k < N; k++)
					for (int l = 0; l < N; l++)
						result(i, j, k, l) = left(i, j) * right(k, l);

		return result;
	}
}
#endif
