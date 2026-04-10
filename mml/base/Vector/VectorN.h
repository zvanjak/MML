///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        VectorN.h                                                           ///
///  Description: Fixed-size vector template for compile-time dimension vectors       ///
///               Static allocation, optimized for small dimensions                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_VECTORN_H
#define MML_VECTORN_H

#include "MMLBase.h"
#include "mml/base/Geometry/Geometry.h"

// Standard headers - include what we use
#include <cassert>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace MML
{
	/// @brief Fixed-size vector with compile-time dimension (stack-allocated).
	/// @tparam Type Element type
	/// @tparam N Vector dimension (fixed at compile time)
	template<class Type, int N>
	class VectorN
	{
	private:
		Type  _val[N] = { 0 };

	protected:
		/// @brief Protected element accessor (non-const) for derived classes.
		Type& val(int i) noexcept { return _val[i]; }
		/// @brief Protected element accessor (const) for derived classes.
		const Type& val(int i) const noexcept { return _val[i]; }

	public:
		typedef Type value_type;      // make T available externally

		///////////////////////          Constructors and destructor       //////////////////////
		/// @brief Default constructor (zero vector).
		VectorN() {}
		/// @brief Constructs vector with all elements set to init_val.
		/// @param init_val Initialization value for all elements
		explicit VectorN(const Type& init_val) {
			for (int i = 0; i < N; ++i)
				_val[i] = init_val;
		}
		/// @brief Constructs vector from initializer list.
		/// @param list Initializer list (must have <= N elements)
		/// @throws VectorDimensionError if list.size() > N
		VectorN(std::initializer_list<Type> list)
		{
			if (static_cast<int>(list.size()) > N)
				throw VectorDimensionError("VectorN initializer list has more elements than vector dimension", N, static_cast<int>(list.size()));
			int count = 0;
			for (auto element : list)
			{
				_val[count] = element;
				++count;
			}
		}
		/// @brief Constructs vector from std::vector.
		/// @param list Vector of values (must have <= N elements)
		/// @throws VectorDimensionError if list.size() > N
		VectorN(std::vector<Type> list)
		{
			if (static_cast<int>(list.size()) > N)
				throw VectorDimensionError("VectorN std::vector has more elements than vector dimension", N, static_cast<int>(list.size()));
			int count{ 0 };
			for (auto element : list)
			{
				_val[count] = element;
				++count;
			}
		}
		/// @brief Constructs vector from raw array.
		/// @param vals Array of N elements
		explicit VectorN(Type* vals)
		{
			for (int i = 0; i < N; ++i)
				_val[i] = vals[i];
		}
		
		/// @brief Returns unit vector (1 at indUnit, 0 elsewhere).
		/// @param indUnit Index for the unit component
		static VectorN UnitVector(int indUnit)
		{
			VectorN ret;
			ret[indUnit] = 1.0;
			return ret;
		}
		/// @brief Returns normalized vector (unit length).
		VectorN Normalized() const
		{
			VectorN ret;
			Real norm = NormL2();
			if (norm < std::numeric_limits<Real>::epsilon() * 100)
				throw VectorDimensionError("VectorN::Normalized - cannot normalize near-zero vector", N, 0);
			for (int i = 0; i < N; ++i)
				ret._val[i] = _val[i] / norm;
			return ret;
		}
		////////////////////////            Standard stuff             ////////////////////////
		/// @brief Returns vector dimension.
		int  size() const noexcept { return N; }
		/// @brief Sets all elements to zero.
		void clear() noexcept {
			for (int i = 0; i < N; ++i)
				_val[i] = Type{ 0 };
		}

		///////////////////////            Accessing elements             ///////////////////////
		/// @brief Element access (non-const).
		/// @param n Index
		inline Type& operator[](int n) noexcept {
			assert(n >= 0 && n < N && "VectorN::operator[] - index out of bounds");
			return _val[n];
		}
		/// @brief Element access (const).
		/// @param n Index
		inline const Type& operator[](int n) const noexcept {
			assert(n >= 0 && n < N && "VectorN::operator[] - index out of bounds");
			return _val[n];
		}

		/// @brief Checked element access (non-const).
		/// @param n Index (throws if out of bounds)
		Type& at(int n) {
			if (n < 0 || n >= N)
				throw VectorDimensionError("VectorN::at - index out of bounds", N, n);
			else
				return _val[n];
		}
		/// @brief Checked element access (const).
		/// @param n Index (throws if out of bounds)
		const Type& at(int n) const {
			if (n < 0 || n >= N)
				throw VectorDimensionError("VectorN::at - index out of bounds", N, n);
			else
				return _val[n];
		}

		///////////////////////            Iterator support                ///////////////////////
		Type*       begin()        noexcept { return _val; }
		Type*       end()          noexcept { return _val + N; }
		const Type* begin()  const noexcept { return _val; }
		const Type* end()    const noexcept { return _val + N; }
		const Type* cbegin() const noexcept { return _val; }
		const Type* cend()   const noexcept { return _val + N; }
		Type*       data()         noexcept { return _val; }
		const Type* data()   const noexcept { return _val; }

		///////////////////////            Testing equality                ///////////////////////
		/// @brief Exact equality comparison.
		/// @param b Vector to compare
		bool operator==(const VectorN& b) const
		{
			for (int i = 0; i < size(); i++)
				if ((*this)[i] != b[i])
					return false;

			return true;
		}
		/// @brief Checks equality within tolerance (static version).
		/// @param a First vector
		/// @param b Second vector
		/// @param eps Tolerance
		static bool AreEqual(const VectorN& a, const VectorN& b, Type eps = Defaults::VectorIsEqualTolerance)
		{
			return a.IsEqualTo(b, eps);
		}
		/// @brief Checks equality within tolerance.
		/// @param b Vector to compare
		/// @param eps Tolerance
		bool IsEqualTo(const VectorN& b, Type eps = Defaults::VectorIsEqualTolerance) const
		{
			for (int i = 0; i < N; i++)
			{
				if (Abs((*this)[i] - b[i]) > eps)
					return false;
			}
			return true;
		}
		/// @brief Inequality operator.
		/// @param b Vector to compare
		bool operator!=(const VectorN& b) const
		{
			return !(*this == b);
		}

		/// @brief Checks if vector is zero (all elements = 0).
		bool isZero() const noexcept
		{
			for (int i = 0; i < N; i++)
				if (std::abs(_val[i]) > std::numeric_limits<Real>::epsilon() * 100)
					return false;

			return true;
		}

		/////////////////////             Arithmetic operators             ///////////////////////
		/// @brief Unary minus (negation).
		[[nodiscard]] VectorN  operator-() const
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = Type{ -1 } *_val[i];
			return ret;
		}
		
		/// @brief Vector addition.
		/// @param b Vector to add
		[[nodiscard]] VectorN  operator+(const VectorN& b) const
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = _val[i] + b._val[i];
			return ret;
		}
		/// @brief In-place addition.
		/// @param b Vector to add
		VectorN& operator+=(const VectorN& b)
		{
			for (int i = 0; i < N; i++)
				_val[i] += b._val[i];
			return *this;
		}
		/// @brief Vector subtraction.
		/// @param b Vector to subtract
		[[nodiscard]] VectorN  operator-(const VectorN& b) const
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = _val[i] - b._val[i];
			return ret;
		}
		/// @brief In-place subtraction.
		/// @param b Vector to subtract
		VectorN& operator-=(const VectorN& b)
		{
			for (int i = 0; i < N; i++)
				_val[i] -= b._val[i];
			return *this;
		}

		/// @brief Scalar multiplication (vector * scalar).
		/// @param b Scalar multiplier
		[[nodiscard]] VectorN  operator*(const Type &b) const
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = _val[i] * b;
			return ret;
		}
		/// @brief In-place scalar multiplication.
		/// @param b Scalar multiplier
		VectorN& operator*=(const Type& b)
		{
			for (int i = 0; i < N; i++)
				_val[i] *= b;
			return *this;
		}
		/// @brief Scalar division (vector / scalar).
		/// @param b Scalar divisor
		[[nodiscard]] VectorN  operator/(const Type &b) const
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = _val[i] / b;
			return ret;
		}		
		/// @brief In-place scalar division.
		/// @param b Scalar divisor
		VectorN& operator/=(const Type& b)
		{
			for (int i = 0; i < N; i++)
				_val[i] /= b;
			return *this;
		}
		
		/// @brief Scalar multiplication (scalar * vector).
		/// @param a Scalar multiplier
		/// @param b Vector
		friend VectorN operator*(Type a, const VectorN<Type, N>& b)
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = a * b[i];
			return ret;
		}

		//////////////////////                 Operations                 ///////////////////////
		/// @brief Computes L1 norm (Manhattan norm).
		Real NormL1() const
		{
			Real norm{ 0.0 };
			for (int i = 0; i < N; i++)
				norm += Abs((*this)[i]);
			return norm;
		}
		/// @brief Computes L2 norm (Euclidean norm): ||v||_2 = √(Σ|v_i|²)
		/// @details For complex vectors, uses std::norm(z) = |z|² = real² + imag²
		Real NormL2() const
		{
			Real norm{ 0.0 };
			for (int i = 0; i < N; i++) {
				if constexpr (std::is_same_v<Type, Complex> || 
				              std::is_same_v<Type, std::complex<float>> ||
				              std::is_same_v<Type, std::complex<long double>>) {
					norm += std::norm((*this)[i]);  // |z|² for complex
				} else {
					norm += (*this)[i] * (*this)[i];  // x² for real
				}
			}
			return std::sqrt(norm);
		}
		/// @brief Computes L-infinity norm (maximum absolute value).
		Real NormLInf() const
		{
			Real norm{ 0.0 };
			for (int i = 0; i < N; i++)
				norm = std::max(norm, Abs((*this)[i]));
			return norm;
		}

		///////////////////////////               I/O                 ///////////////////////////
		/// @brief Prints vector to stream.
		/// @param stream Output stream
		/// @param width Field width
		/// @param precision Decimal precision
		std::ostream& Print(std::ostream& stream, int width, int precision) const
		{
			stream << "[";
			bool first = true;
			for (const Type& x : _val)
			{
				if (!first)
					stream << ", ";
				else
					first = false;

				stream << std::setw(width) << std::setprecision(precision) << x;
			}
			stream << "]";

			return stream;
		}
		/// @brief Prints vector with values below threshold shown as zero.
		/// @param stream Output stream
		/// @param width Field width
		/// @param precision Decimal precision
		/// @param zeroThreshold Threshold for displaying as zero
		std::ostream& Print(std::ostream& stream, int width, int precision, Real zeroThreshold) const
		{
			stream << "[";
			bool first = true;
			for (const Type& x : _val)
			{
				if (!first)
					stream << ", ";
				else
					first = false;

				if (Abs(x) > zeroThreshold)
					stream << std::setw(width) << std::setprecision(precision) << x;
				else
					stream << std::setw(width) << std::setprecision(precision) << 0.0;
			}
			stream << "]";

			return stream;
		}
		/// @brief Prints vector with message prefix on single line.
		/// @param stream Output stream
		/// @param msg Message prefix
		/// @param width Field width
		/// @param precision Decimal precision
		std::ostream& PrintLine(std::ostream& stream, std::string msg, int width, int precision) const
		{
			stream << msg;
			Print(stream, width, precision);
			stream << std::endl;

			return stream;
		}

		/// @brief Converts vector to string.
		/// @param width Field width
		/// @param precision Decimal precision
		std::string to_string(int width, int precision) const
		{
			std::stringstream str;

			Print(str, width, precision);

			return str.str();
		}
		/// @brief Stream output operator.
		friend std::ostream& operator<<(std::ostream& stream, const VectorN<Type, N>& a)
		{
			a.Print(stream, Defaults::VectorNPrintWidth, Defaults::VectorNPrintPrecision);

			return stream;
		}
	};

	typedef VectorN<Real, 2> Vec2;
	typedef VectorN<Real, 3> Vec3;
	typedef VectorN<Real, 4> Vec4;

	typedef VectorN<float, 2> Vec2Flt;
	typedef VectorN<float, 3> Vec3Flt;
	typedef VectorN<float, 4> Vec4Flt;

	typedef VectorN<double, 2> Vec2Dbl;
	typedef VectorN<double, 3> Vec3Dbl;
	typedef VectorN<double, 4> Vec4Dbl;

	typedef VectorN<Complex, 2> Vec2Complex;
	typedef VectorN<Complex, 3> Vec3Complex;
	typedef VectorN<Complex, 4> Vec4Complex;

	typedef VectorN<float, 2> Vec2F;
	typedef VectorN<float, 3> Vec3F;
	typedef VectorN<float, 4> Vec4F;

	typedef VectorN<double, 2> Vec2D;
	typedef VectorN<double, 3> Vec3D;
	typedef VectorN<double, 4> Vec4D;

	typedef VectorN<Complex, 2> Vec2C;
	typedef VectorN<Complex, 3> Vec3C;
	typedef VectorN<Complex, 4> Vec4C;
}

#endif
