///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Tensor3.h                                                          ///
///  Description: Rank-3 concrete tensor and rank-specific operations                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_TENSOR3_H
#define MML_TENSOR3_H

#include <array>
#include <cassert>
#include <cmath>
#include <initializer_list>
#include <utility>

#include <mml/MMLBase.h>
#include <mml/interfaces/ITensor.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/BaseUtils/SymbolUtils.h>
#include <mml/base/Tensor/TensorCommon.h>
#include <mml/base/Tensor/MatrixNDim.h>

namespace MML
{
	///////////////////////////////////////////////////////////////////////////
	///                            Tensor3                                  ///
	///////////////////////////////////////////////////////////////////////////

	/// @brief Rank-3 tensor in N-dimensional space with covariant/contravariant tracking.
	/// @details Represents tensors with 3 indices, commonly used for:
	///          - Structure constants of Lie algebras
	///          - Completely antisymmetric tensors (Levi-Civita symbol)
	///          - Third-order material tensors
	///
	/// @note Christoffel symbols are often stored in rank-3 arrays but are not
	///       tensors; do not transform them with tensor transformation APIs.
	///
	///          Index variance is tracked per-position. Contraction reduces rank by 2,
	///          returning a VectorN (rank-1 tensor).
	///
	/// @tparam N Dimension of the underlying vector space
	///
	/// @par Example Usage:
	/// @code
	///     // Create a true rank-3 tensor with 1 up, 2 down indices
	///     Tensor3<3> torsion(2, 1);  // 2 covariant, 1 contravariant
	///     
	///     // Contract indices 1 and 2 (trace over the covariant indices)
	///     VectorN<Real, 3> contracted = torsion.Contract(1, 2);
	/// @endcode
	template <int N>
	class Tensor3 : public ITensor3<N>
	{
		MatrixNDim<Real, N, 3> _coeff;  ///< Storage for N³ tensor components
		int _numContravar;             ///< Number of contravariant (upper) indices
		int _numCovar;                 ///< Number of covariant (lower) indices
		bool _isContravar[3];          ///< Per-index variance: true = contravariant, false = covariant

		/// @brief Zero tensor with this tensor's exact per-index variance arrangement.
		Tensor3 sameShapeZero() const
		{
			Tensor3 r(_numCovar, _numContravar);
			for (int i = 0; i < 3; i++) r.setContravar(i, _isContravar[i]);
			return r;
		}

		/// @brief Throw if `other` does not share this tensor's per-index variance.
		void checkSameVariance(const Tensor3& other, const char* where) const
		{
			for (int i = 0; i < 3; i++)
				if (_isContravar[i] != other._isContravar[i])
					throw TensorCovarContravarArithmeticError(where, _numContravar, _numCovar, other._numContravar, other._numCovar);
		}
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
			for (int i = 0; i < N; ++i)
				for (int j = 0; j < N; ++j)
					for (int k = 0; k < N; ++k)
						if (val != values.end())
							_coeff(i, j, k) = *val++;
						else
							_coeff(i, j, k) = 0.0;
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
		Real  operator()(int i, int j, int k) const override { 
			assert(i >= 0 && i < N && "Tensor3: index i out of bounds");
			assert(j >= 0 && j < N && "Tensor3: index j out of bounds");
			assert(k >= 0 && k < N && "Tensor3: index k out of bounds");
			return _coeff(i, j, k); 
		}

		/// @brief Unchecked element access (mutable) — fast, assertion-only bounds checking.
		/// @warning Assertions are disabled in Release builds. For runtime safety, use at().
		Real& operator()(int i, int j, int k) override { 
			assert(i >= 0 && i < N && "Tensor3: index i out of bounds");
			assert(j >= 0 && j < N && "Tensor3: index j out of bounds");
			assert(k >= 0 && k < N && "Tensor3: index k out of bounds");
			return _coeff(i, j, k); 
		}
		
		/// @brief Checked element access (const) — safe, throws on out-of-bounds.
		/// @details Use this instead of operator() when runtime bounds safety is needed.
		/// @throws TensorIndexError if any index is out of bounds
		Real at(int i, int j, int k) const {
			if (i < 0 || i >= N || j < 0 || j >= N || k < 0 || k >= N)
				throw TensorIndexError("Tensor3::at - index out of bounds");
			return _coeff(i, j, k);
		}

		/// @brief Checked element access (mutable) — safe, throws on out-of-bounds.
		/// @details Use this instead of operator() when runtime bounds safety is needed.
		/// @throws TensorIndexError if any index is out of bounds
		Real& at(int i, int j, int k) {
			if (i < 0 || i >= N || j < 0 || j >= N || k < 0 || k >= N)
				throw TensorIndexError("Tensor3::at - index out of bounds");
			return _coeff(i, j, k);
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
						sum += _coeff(i, j, k) * v1[i] * v2[j] * v3[k];

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
					sum += _coeff(idx[0], idx[1], idx[2]);
				}
				result[a] = sum;
			}

			return result;
		}

		///////////////////////////////////////////////////////////////////////
		///                      Arithmetic Operations                      ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Add two tensors sharing the same per-index variance.
		/// @throws TensorCovarContravarArithmeticError if index variance doesn't match
		Tensor3 operator+(const Tensor3& other) const
		{
			checkSameVariance(other, "Tensor3 operator+, index variance mismatch");
			Tensor3 result = sameShapeZero();
			result._coeff = _coeff + other._coeff;
			return result;
		}

		/// @brief Subtract two tensors sharing the same per-index variance.
		/// @throws TensorCovarContravarArithmeticError if index variance doesn't match
		Tensor3 operator-(const Tensor3& other) const
		{
			checkSameVariance(other, "Tensor3 operator-, index variance mismatch");
			Tensor3 result = sameShapeZero();
			result._coeff = _coeff - other._coeff;
			return result;
		}

		/// @brief Unary negation.
		Tensor3 operator-() const
		{
			Tensor3 result = sameShapeZero();
			result._coeff = -_coeff;
			return result;
		}

		/// @brief Multiply tensor by a scalar.
		Tensor3 operator*(Real scalar) const
		{
			Tensor3 result = sameShapeZero();
			result._coeff = _coeff * scalar;
			return result;
		}

		/// @brief Divide tensor by a scalar.
		Tensor3 operator/(Real scalar) const
		{
			Tensor3 result = sameShapeZero();
			result._coeff = _coeff / scalar;
			return result;
		}

		/// @brief Scalar multiplication from the left (scalar * tensor).
		friend Tensor3 operator*(Real scalar, const Tensor3& t)
		{
			Tensor3 result = t.sameShapeZero();
			result._coeff = t._coeff * scalar;
			return result;
		}

		Tensor3& operator+=(const Tensor3& other) { checkSameVariance(other, "Tensor3 operator+=, index variance mismatch"); _coeff += other._coeff; return *this; }
		Tensor3& operator-=(const Tensor3& other) { checkSameVariance(other, "Tensor3 operator-=, index variance mismatch"); _coeff -= other._coeff; return *this; }
		Tensor3& operator*=(Real scalar) { _coeff *= scalar; return *this; }
		Tensor3& operator/=(Real scalar) { _coeff /= scalar; return *this; }

		/// @brief Exact equality: same per-index variance and identical components.
		bool operator==(const Tensor3& other) const
		{
			for (int i = 0; i < 3; i++)
				if (_isContravar[i] != other._isContravar[i])
					return false;
			return _coeff == other._coeff;
		}
		bool operator!=(const Tensor3& other) const { return !(*this == other); }
	};

	/// @brief Permute rank-3 tensor indices.
	/// @details permutation[resultSlot] gives the source slot that becomes that result slot.
	template<int N>
	Tensor3<N> PermuteIndices(const Tensor3<N>& tensor, const std::array<int, 3>& permutation)
	{
		detail::ValidateTensorPermutation<3>(permutation);

		bool resultContravar[3] = {
			tensor.IsContravar(permutation[0]),
			tensor.IsContravar(permutation[1]),
			tensor.IsContravar(permutation[2])
		};
		int numContravar = 0;
		for (int i = 0; i < 3; i++)
			if (resultContravar[i])
				numContravar++;

		Tensor3<N> result(3 - numContravar, numContravar);
		for (int i = 0; i < 3; i++)
			result.setContravar(i, resultContravar[i]);

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				for (int k = 0; k < N; k++)
				{
					int sourceIndex[3];
					sourceIndex[permutation[0]] = i;
					sourceIndex[permutation[1]] = j;
					sourceIndex[permutation[2]] = k;
					result(i, j, k) = tensor(sourceIndex[0], sourceIndex[1], sourceIndex[2]);
				}

		return result;
	}

	/// @brief Covariant Levi-Civita tensor in 3D: E_ijk = sqrt(|det(g)|) epsilon_ijk.
	template<int N>
	Tensor3<N> LeviCivitaTensor(const MatrixNM<Real, N, N>& covariantMetric)
	{
		static_assert(N == 3, "Tensor3 LeviCivitaTensor is available only for N=3");

		Real weight = std::sqrt(std::abs(detail::DeterminantForLeviCivitaTensor<N>(covariantMetric)));
		Tensor3<N> result(3, 0);

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				for (int k = 0; k < N; k++)
					result(i, j, k) = weight * Utils::LeviCivita(i + 1, j + 1, k + 1);

		return result;
	}
}
#endif
