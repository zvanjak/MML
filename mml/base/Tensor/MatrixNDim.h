///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MatrixNDim.h                                                        ///
///  Description: Fixed-size rank-R, N-per-axis dense storage for concrete tensors   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_MATRIX_NDIM_H
#define MML_MATRIX_NDIM_H

#include <array>
#include <cassert>
#include <cmath>

#include <mml/MMLBase.h>

namespace MML
{
	namespace detail
	{
		/// @brief Compile-time integer power base^exp (exp >= 0).
		constexpr int IntPow(int base, int exp)
		{
			int result = 1;
			for (int i = 0; i < exp; ++i)
				result *= base;
			return result;
		}
	}

	///////////////////////////////////////////////////////////////////////////
	///                            MatrixNDim                               ///
	///////////////////////////////////////////////////////////////////////////

	/// @brief Fixed-size rank-`Rank` dense array with `N` elements along every axis.
	/// @details Provides the shared component storage for the concrete tensor classes
	///          (Tensor2..Tensor5). Components live in a single row-major flat buffer of
	///          `N^Rank` elements, so element-wise arithmetic is written once here instead
	///          of being duplicated across each tensor rank.
	///
	///          Multi-index access uses a variadic `operator()` taking exactly `Rank`
	///          integer indices; the row-major offset matches the layout of a C array
	///          `T[N][N]...[N]`, so initializer order is preserved when migrating.
	///
	/// @tparam T    Component type (Real, float, Complex, ...)
	/// @tparam N    Number of components along each axis (dimension of the space)
	/// @tparam Rank Number of indices / axes (tensor rank)
	template<typename T, int N, int Rank>
	class MatrixNDim
	{
		static_assert(N > 0, "MatrixNDim requires N > 0");
		static_assert(Rank > 0, "MatrixNDim requires Rank > 0");

	public:
		static constexpr int Dim = N;                                ///< Elements per axis
		static constexpr int TensorRank = Rank;                      ///< Number of axes
		static constexpr int TotalSize = detail::IntPow(N, Rank);    ///< Flat element count (N^Rank)

		typedef T value_type;                                        ///< STL-style element alias

	private:
		std::array<T, TotalSize> _data{};   ///< Row-major flat storage (zero-initialized)

		/// @brief Row-major flat offset from `Rank` indices.
		template<typename... Idx>
		static int offset(Idx... idx)
		{
			static_assert(sizeof...(Idx) == Rank, "MatrixNDim: wrong number of indices");
			const int indices[Rank] = { static_cast<int>(idx)... };
			int off = 0;
			for (int d = 0; d < Rank; ++d)
			{
				assert(indices[d] >= 0 && indices[d] < N && "MatrixNDim: index out of bounds");
				off = off * N + indices[d];
			}
			return off;
		}

	public:
		///////////////////////////////////////////////////////////////////////
		///                        Element Access                           ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Multi-index access (mutable). Requires exactly `Rank` indices.
		template<typename... Idx>
		T& operator()(Idx... idx) { return _data[offset(idx...)]; }

		/// @brief Multi-index access (const). Requires exactly `Rank` indices.
		template<typename... Idx>
		const T& operator()(Idx... idx) const { return _data[offset(idx...)]; }

		/// @brief Flat access into the underlying buffer (mutable).
		T& flat(int i) { assert(i >= 0 && i < TotalSize && "MatrixNDim::flat out of bounds"); return _data[i]; }

		/// @brief Flat access into the underlying buffer (const).
		const T& flat(int i) const { assert(i >= 0 && i < TotalSize && "MatrixNDim::flat out of bounds"); return _data[i]; }

		/// @brief Number of stored components (N^Rank).
		constexpr int size() const { return TotalSize; }

		///////////////////////////////////////////////////////////////////////
		///                      Arithmetic Operations                      ///
		///////////////////////////////////////////////////////////////////////

		MatrixNDim operator+(const MatrixNDim& other) const
		{
			MatrixNDim result;
			for (int i = 0; i < TotalSize; ++i)
				result._data[i] = _data[i] + other._data[i];
			return result;
		}

		MatrixNDim operator-(const MatrixNDim& other) const
		{
			MatrixNDim result;
			for (int i = 0; i < TotalSize; ++i)
				result._data[i] = _data[i] - other._data[i];
			return result;
		}

		MatrixNDim operator-() const
		{
			MatrixNDim result;
			for (int i = 0; i < TotalSize; ++i)
				result._data[i] = -_data[i];
			return result;
		}

		MatrixNDim operator*(T scalar) const
		{
			MatrixNDim result;
			for (int i = 0; i < TotalSize; ++i)
				result._data[i] = _data[i] * scalar;
			return result;
		}

		MatrixNDim operator/(T scalar) const
		{
			MatrixNDim result;
			for (int i = 0; i < TotalSize; ++i)
				result._data[i] = _data[i] / scalar;
			return result;
		}

		friend MatrixNDim operator*(T scalar, const MatrixNDim& m)
		{
			MatrixNDim result;
			for (int i = 0; i < TotalSize; ++i)
				result._data[i] = scalar * m._data[i];
			return result;
		}

		MatrixNDim& operator+=(const MatrixNDim& other)
		{
			for (int i = 0; i < TotalSize; ++i)
				_data[i] += other._data[i];
			return *this;
		}

		MatrixNDim& operator-=(const MatrixNDim& other)
		{
			for (int i = 0; i < TotalSize; ++i)
				_data[i] -= other._data[i];
			return *this;
		}

		MatrixNDim& operator*=(T scalar)
		{
			for (int i = 0; i < TotalSize; ++i)
				_data[i] *= scalar;
			return *this;
		}

		MatrixNDim& operator/=(T scalar)
		{
			for (int i = 0; i < TotalSize; ++i)
				_data[i] /= scalar;
			return *this;
		}

		///////////////////////////////////////////////////////////////////////
		///                          Comparison                             ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Exact component-wise equality.
		bool operator==(const MatrixNDim& other) const
		{
			for (int i = 0; i < TotalSize; ++i)
				if (_data[i] != other._data[i])
					return false;
			return true;
		}

		bool operator!=(const MatrixNDim& other) const { return !(*this == other); }

		/// @brief Approximate component-wise equality within an absolute tolerance.
		bool IsEqualTo(const MatrixNDim& other, T tolerance) const
		{
			for (int i = 0; i < TotalSize; ++i)
				if (std::abs(_data[i] - other._data[i]) > tolerance)
					return false;
			return true;
		}
	};
}
#endif
