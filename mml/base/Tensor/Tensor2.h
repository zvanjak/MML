///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Tensor2.h                                                          ///
///  Description: Rank-2 concrete tensor and rank-specific operations                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_TENSOR2_H
#define MML_TENSOR2_H

#include <array>
#include <cassert>
#include <iomanip>
#include <initializer_list>
#include <iostream>

#include <mml/MMLBase.h>
#include <mml/interfaces/ITensor.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/Tensor/TensorCommon.h>
#include <mml/base/Tensor/MatrixNDim.h>

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
		MatrixNDim<Real, N, 2> _coeff;   ///< Storage for tensor components
		int _numContravar = 0;         ///< Number of contravariant (upper) indices
		int _numCovar = 0;             ///< Number of covariant (lower) indices
		bool _isContravar[2];          ///< Per-index variance: true = contravariant, false = covariant

		/// @brief Zero tensor with this tensor's exact per-index variance arrangement.
		Tensor2 sameShapeZero() const
		{
			Tensor2 r(_numCovar, _numContravar);
			for (int i = 0; i < 2; i++) r.setContravar(i, _isContravar[i]);
			return r;
		}

		/// @brief Throw if `other` does not share this tensor's per-index variance.
		void checkSameVariance(const Tensor2& other, const char* where) const
		{
			for (int i = 0; i < 2; i++)
				if (_isContravar[i] != other._isContravar[i])
					throw TensorCovarContravarArithmeticError(where, _numContravar, _numCovar, other._numContravar, other._numCovar);
		}
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
			for (int i = 0; i < N; ++i)
				for (int j = 0; j < N; ++j)
					if (val != values.end())
						_coeff(i, j) = *val++;
					else
						_coeff(i, j) = 0.0;
		}

		///////////////////////////////////////////////////////////////////////
		///                     Index Variance Query                        ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Returns the number of contravariant (upper) indices.
		int  NumContravar() const override { return _numContravar; }

		/// @brief Returns the number of covariant (lower) indices.
		int  NumCovar()     const override { return _numCovar; }

		/// @brief Returns the underlying coefficient matrix.
		MatrixNM<Real, N, N> GetMatrix() const {
			MatrixNM<Real, N, N> m;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					m[i][j] = _coeff(i, j);
			return m;
		}

		/// @brief Checks if index i is contravariant (upper index).
		/// @param i Index position (0 or 1)
		bool IsContravar(int i) const override { return _isContravar[i]; }

		/// @brief Checks if index i is covariant (lower index).
		/// @param i Index position (0 or 1)
		bool IsCovar(int i) const			{ return !_isContravar[i]; }

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
		/// @param i First index (0 to N-1)
		/// @param j Second index (0 to N-1)
		/// @return The tensor component T(i,j)
		Real  operator()(int i, int j) const override { 
			assert(i >= 0 && i < N && "Tensor2: index i out of bounds");
			assert(j >= 0 && j < N && "Tensor2: index j out of bounds");
			return _coeff(i, j); 
		}

		/// @brief Unchecked element access (mutable) — fast, assertion-only bounds checking.
		/// @warning Assertions are disabled in Release builds. For runtime safety, use at().
		/// @param i First index (0 to N-1)
		/// @param j Second index (0 to N-1)
		/// @return Reference to the tensor component T(i,j)
		Real& operator()(int i, int j) override { 
			assert(i >= 0 && i < N && "Tensor2: index i out of bounds");
			assert(j >= 0 && j < N && "Tensor2: index j out of bounds");
			return _coeff(i, j); 
		}
		
		/// @brief Checked element access (const) — safe, throws on out-of-bounds.
		/// @details Use this instead of operator() when runtime bounds safety is needed.
		/// @param i First index (0 to N-1)
		/// @param j Second index (0 to N-1)
		/// @return The tensor component T(i,j)
		/// @throws TensorIndexError if any index is out of bounds
		Real at(int i, int j) const {
			if (i < 0 || i >= N || j < 0 || j >= N)
				throw TensorIndexError("Tensor2::at - index out of bounds");
			return _coeff(i, j);
		}

		/// @brief Checked element access (mutable) — safe, throws on out-of-bounds.
		/// @details Use this instead of operator() when runtime bounds safety is needed.
		/// @param i First index (0 to N-1)
		/// @param j Second index (0 to N-1)
		/// @return Reference to the tensor component T(i,j)
		/// @throws TensorIndexError if any index is out of bounds
		Real& at(int i, int j) {
			if (i < 0 || i >= N || j < 0 || j >= N)
				throw TensorIndexError("Tensor2::at - index out of bounds");
			return _coeff(i, j);
		}

		///////////////////////////////////////////////////////////////////////
		///                      Arithmetic Operations                      ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Add two tensors sharing the same per-index variance.
		/// @throws TensorCovarContravarArithmeticError if index variance doesn't match
		Tensor2 operator+(const Tensor2& other) const
		{
			checkSameVariance(other, "Tensor2 operator+, index variance mismatch");
			Tensor2 result = sameShapeZero();
			result._coeff = _coeff + other._coeff;
			return result;
		}

		/// @brief Subtract two tensors sharing the same per-index variance.
		/// @throws TensorCovarContravarArithmeticError if index variance doesn't match
		Tensor2 operator-(const Tensor2& other) const
		{
			checkSameVariance(other, "Tensor2 operator-, index variance mismatch");
			Tensor2 result = sameShapeZero();
			result._coeff = _coeff - other._coeff;
			return result;
		}

		/// @brief Unary negation.
		Tensor2 operator-() const
		{
			Tensor2 result = sameShapeZero();
			result._coeff = -_coeff;
			return result;
		}

		/// @brief Multiply tensor by a scalar.
		Tensor2 operator*(Real scalar) const
		{
			Tensor2 result = sameShapeZero();
			result._coeff = _coeff * scalar;
			return result;
		}

		/// @brief Divide tensor by a scalar.
		Tensor2 operator/(Real scalar) const
		{
			Tensor2 result = sameShapeZero();
			result._coeff = _coeff / scalar;
			return result;
		}

		/// @brief Scalar multiplication from the left (scalar * tensor).
		friend Tensor2 operator*(Real scalar, const Tensor2& b)
		{
			Tensor2 result = b.sameShapeZero();
			result._coeff = b._coeff * scalar;
			return result;
		}

		Tensor2& operator+=(const Tensor2& other) { checkSameVariance(other, "Tensor2 operator+=, index variance mismatch"); _coeff += other._coeff; return *this; }
		Tensor2& operator-=(const Tensor2& other) { checkSameVariance(other, "Tensor2 operator-=, index variance mismatch"); _coeff -= other._coeff; return *this; }
		Tensor2& operator*=(Real scalar) { _coeff *= scalar; return *this; }
		Tensor2& operator/=(Real scalar) { _coeff /= scalar; return *this; }

		/// @brief Exact equality: same per-index variance and identical components.
		bool operator==(const Tensor2& other) const
		{
			for (int i = 0; i < 2; i++)
				if (_isContravar[i] != other._isContravar[i])
					return false;
			return _coeff == other._coeff;
		}
		bool operator!=(const Tensor2& other) const { return !(*this == other); }

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
				result += _coeff(i, i);

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
					sum += _coeff(i, j) * v1[i] * v2[j];

			return sum;
		}

		///////////////////////////////////////////////////////////////////////
		///                    Index Raising / Lowering                      ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Raise the first index using a metric tensor.
		/// @details Converts T_ij -> T^i_j using metric inverse: T^i_j = g^ik * T_kj
		///          Or converts T_i^j -> T^ij using: T^ij = g^ik * T_k^j
		/// @param metricInverse The inverse metric tensor g^{ij} (must be (2,0) type)
		/// @return New tensor with first index raised
		/// @throws TensorCovarContravarNumError if first index is already contravariant
		///         or metric is not (2,0) type
		Tensor2 RaiseFirstIndex(const Tensor2& metricInverse) const
		{
			if (_isContravar[0])
				throw TensorCovarContravarNumError("Tensor2::RaiseFirstIndex - first index is already contravariant", _numCovar, _numContravar);
			if (metricInverse.NumContravar() != 2)
				throw TensorCovarContravarNumError("Tensor2::RaiseFirstIndex - metric must be fully contravariant (2,0)", metricInverse.NumCovar(), metricInverse.NumContravar());

			Tensor2 result(_numCovar - 1, _numContravar + 1);
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++) {
					Real sum = 0;
					for (int k = 0; k < N; k++)
						sum += metricInverse(i, k) * _coeff(k, j);
					result(i, j) = sum;
				}
			return result;
		}

		/// @brief Raise the second index using a metric tensor.
		/// @details Converts T_ij -> T_i^j using: T_i^j = g^jk * T_ik
		///          Or converts T^i_j -> T^ij using: T^ij = g^jk * T^i_k
		/// @param metricInverse The inverse metric tensor g^{ij} (must be (2,0) type)
		/// @return New tensor with second index raised
		/// @throws TensorCovarContravarNumError if second index is already contravariant
		Tensor2 RaiseSecondIndex(const Tensor2& metricInverse) const
		{
			if (_isContravar[1])
				throw TensorCovarContravarNumError("Tensor2::RaiseSecondIndex - second index is already contravariant", _numCovar, _numContravar);
			if (metricInverse.NumContravar() != 2)
				throw TensorCovarContravarNumError("Tensor2::RaiseSecondIndex - metric must be fully contravariant (2,0)", metricInverse.NumCovar(), metricInverse.NumContravar());

			Tensor2 result(_numCovar - 1, _numContravar + 1);
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++) {
					Real sum = 0;
					for (int k = 0; k < N; k++)
						sum += _coeff(i, k) * metricInverse(k, j);
					result(i, j) = sum;
				}
			return result;
		}

		/// @brief Lower the first index using a metric tensor.
		/// @details Converts T^ij -> T_i^j using: T_i^j = g_ik * T^kj
		///          Or converts T^i_j -> T_ij using: T_ij = g_ik * T^k_j
		/// @param metric The metric tensor g_{ij} (must be (0,2) type)
		/// @return New tensor with first index lowered
		/// @throws TensorCovarContravarNumError if first index is already covariant
		Tensor2 LowerFirstIndex(const Tensor2& metric) const
		{
			if (!_isContravar[0])
				throw TensorCovarContravarNumError("Tensor2::LowerFirstIndex - first index is already covariant", _numCovar, _numContravar);
			if (metric.NumCovar() != 2)
				throw TensorCovarContravarNumError("Tensor2::LowerFirstIndex - metric must be fully covariant (0,2)", metric.NumCovar(), metric.NumContravar());

			Tensor2 result(_numCovar + 1, _numContravar - 1);
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++) {
					Real sum = 0;
					for (int k = 0; k < N; k++)
						sum += metric(i, k) * _coeff(k, j);
					result(i, j) = sum;
				}
			return result;
		}

		/// @brief Lower the second index using a metric tensor.
		/// @details Converts T^ij -> T^i_j using: T^i_j = T^ik * g_kj
		///          Or converts T_i^j -> T_ij using: T_ij = T_ik * g_kj  (with appropriate index match)
		/// @param metric The metric tensor g_{ij} (must be (0,2) type)
		/// @return New tensor with second index lowered
		/// @throws TensorCovarContravarNumError if second index is already covariant
		Tensor2 LowerSecondIndex(const Tensor2& metric) const
		{
			if (!_isContravar[1])
				throw TensorCovarContravarNumError("Tensor2::LowerSecondIndex - second index is already covariant", _numCovar, _numContravar);
			if (metric.NumCovar() != 2)
				throw TensorCovarContravarNumError("Tensor2::LowerSecondIndex - metric must be fully covariant (0,2)", metric.NumCovar(), metric.NumContravar());

			Tensor2 result(_numCovar + 1, _numContravar - 1);
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++) {
					Real sum = 0;
					for (int k = 0; k < N; k++)
						sum += _coeff(i, k) * metric(k, j);
					result(i, j) = sum;
				}
			return result;
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
					stream << std::setw(width) << std::setprecision(precision) << _coeff(i, j) << ", ";
				stream << " ]" << std::endl;
			}
		}
		friend std::ostream& operator<<(std::ostream& stream, const Tensor2& a)
		{
			a.Print(stream, 15, 10);

			return stream;
		}
	};

	/// @brief Contract one index from each of two rank-2 tensors, producing a rank-2 tensor.
	/// @details Computes the single-index contraction of two distinct tensors. The remaining
	///          index from the left tensor becomes result index 0; the remaining index from
	///          the right tensor becomes result index 1.
	/// @param left Left tensor
	/// @param leftIndex Index of left tensor to contract (0 or 1)
	/// @param right Right tensor
	/// @param rightIndex Index of right tensor to contract (0 or 1)
	/// @return Rank-2 tensor containing sum_s left(...s...) * right(...s...)
	/// @throws TensorIndexError if either index is outside [0,1]
	/// @throws TensorCovarContravarArithmeticError if contracted indices have the same variance
	template<int N>
	Tensor2<N> Contract(const Tensor2<N>& left, int leftIndex, const Tensor2<N>& right, int rightIndex)
	{
		if (leftIndex < 0 || leftIndex > 1 || rightIndex < 0 || rightIndex > 1)
			throw TensorIndexError("Contract(Tensor2,Tensor2), index out of range [0,1]");

		if (left.IsContravar(leftIndex) == right.IsContravar(rightIndex))
			throw TensorCovarContravarArithmeticError("Contract(Tensor2,Tensor2), contracted indices must have opposite variance",
				left.NumContravar(), left.NumCovar(), right.NumContravar(), right.NumCovar());

		int leftRemaining = 1 - leftIndex;
		int rightRemaining = 1 - rightIndex;
		bool resultContravar[2] = { left.IsContravar(leftRemaining), right.IsContravar(rightRemaining) };
		int numContravar = (resultContravar[0] ? 1 : 0) + (resultContravar[1] ? 1 : 0);
		Tensor2<N> result(2 - numContravar, numContravar);
		result.setContravar(0, resultContravar[0]);
		result.setContravar(1, resultContravar[1]);

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				Real sum = REAL(0.0);
				for (int s = 0; s < N; s++)
				{
					int leftIdx[2];
					int rightIdx[2];
					leftIdx[leftIndex] = s;
					leftIdx[leftRemaining] = i;
					rightIdx[rightIndex] = s;
					rightIdx[rightRemaining] = j;

					sum += left(leftIdx[0], leftIdx[1]) * right(rightIdx[0], rightIdx[1]);
				}
				result(i, j) = sum;
			}

		return result;
	}

	/// @brief Outer product of two rank-1 vectors, producing a rank-2 tensor.
	/// @details VectorN does not encode index variance, so the caller may specify
	///          whether each vector should become a contravariant or covariant index.
	///          By default, vectors are treated as contravariant components.
	template<int N>
	Tensor2<N> OuterProduct(const VectorN<Real, N>& left, const VectorN<Real, N>& right,
		TensorIndexType leftIndexType = CONTRAVARIANT, TensorIndexType rightIndexType = CONTRAVARIANT)
	{
		bool resultContravar[2] = { leftIndexType == CONTRAVARIANT, rightIndexType == CONTRAVARIANT };
		int numContravar = (resultContravar[0] ? 1 : 0) + (resultContravar[1] ? 1 : 0);
		Tensor2<N> result(2 - numContravar, numContravar);
		result.setContravar(0, resultContravar[0]);
		result.setContravar(1, resultContravar[1]);

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				result(i, j) = left[i] * right[j];

		return result;
	}

	/// @brief Symmetric part of a rank-2 tensor: (T_ij + T_ji) / 2.
	/// @details Symmetrization is only defined when both index slots have the same variance.
	template<int N>
	Tensor2<N> SymmetricPart(const Tensor2<N>& tensor)
	{
		if (tensor.IsContravar(0) != tensor.IsContravar(1))
			throw TensorCovarContravarArithmeticError("SymmetricPart requires both tensor indices to have the same variance",
				tensor.NumContravar(), tensor.NumCovar(), tensor.NumContravar(), tensor.NumCovar());

		Tensor2<N> result(tensor.NumCovar(), tensor.NumContravar());
		result.setContravar(0, tensor.IsContravar(0));
		result.setContravar(1, tensor.IsContravar(1));

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				result(i, j) = REAL(0.5) * (tensor(i, j) + tensor(j, i));

		return result;
	}

	/// @brief Antisymmetric part of a rank-2 tensor: (T_ij - T_ji) / 2.
	/// @details Antisymmetrization is only defined when both index slots have the same variance.
	template<int N>
	Tensor2<N> AntisymmetricPart(const Tensor2<N>& tensor)
	{
		if (tensor.IsContravar(0) != tensor.IsContravar(1))
			throw TensorCovarContravarArithmeticError("AntisymmetricPart requires both tensor indices to have the same variance",
				tensor.NumContravar(), tensor.NumCovar(), tensor.NumContravar(), tensor.NumCovar());

		Tensor2<N> result(tensor.NumCovar(), tensor.NumContravar());
		result.setContravar(0, tensor.IsContravar(0));
		result.setContravar(1, tensor.IsContravar(1));

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				result(i, j) = REAL(0.5) * (tensor(i, j) - tensor(j, i));

		return result;
	}

	/// @brief Permute rank-2 tensor indices.
	/// @details permutation[resultSlot] gives the source slot that becomes that result slot.
	template<int N>
	Tensor2<N> PermuteIndices(const Tensor2<N>& tensor, const std::array<int, 2>& permutation)
	{
		detail::ValidateTensorPermutation<2>(permutation);

		bool resultContravar[2] = { tensor.IsContravar(permutation[0]), tensor.IsContravar(permutation[1]) };
		int numContravar = (resultContravar[0] ? 1 : 0) + (resultContravar[1] ? 1 : 0);
		Tensor2<N> result(2 - numContravar, numContravar);
		result.setContravar(0, resultContravar[0]);
		result.setContravar(1, resultContravar[1]);

		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				int sourceIndex[2];
				sourceIndex[permutation[0]] = i;
				sourceIndex[permutation[1]] = j;
				result(i, j) = tensor(sourceIndex[0], sourceIndex[1]);
			}

		return result;
	}

	/// @brief Transpose a rank-2 tensor by swapping its two index slots.
	template<int N>
	Tensor2<N> Transpose(const Tensor2<N>& tensor)
	{
		return PermuteIndices(tensor, std::array<int, 2>{ 1, 0 });
	}
}
#endif
