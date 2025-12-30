///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        VectorN.h                                                           ///
///  Description: Fixed-size vector template for compile-time dimension vectors       ///
///               Static allocation, optimized for small dimensions                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_VECTORN_H
#define MML_VECTORN_H

#include "MMLBase.h"
#include "base/Geometry.h"

#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <vector>

namespace MML
{
	template<class Type, int N>
	class VectorN
	{
	protected:
		Type  _val[N] = { 0 };

	public:
		typedef Type value_type;      // make T available externally

		///////////////////////          Constructors and destructor       //////////////////////
		VectorN() {}
		explicit VectorN(const Type& init_val) {
			for (int i = 0; i < N; ++i)
				_val[i] = init_val;
		}
		VectorN(std::initializer_list<Type> list)
		{
			int count = 0;
			for (auto element : list)
			{
				_val[count] = element;
				++count;

				if (count >= N)
					break;
			}
		}
		VectorN(std::vector<Type> list)
		{
			int count{ 0 };
			for (auto element : list)
			{
				_val[count] = element;
				++count;

				if (count >= N)
					break;
			}
		}
		explicit VectorN(Type* vals)
		{
			for (int i = 0; i < N; ++i)
				_val[i] = vals[i];
		}
		
		static VectorN GetUnitVector(int indUnit)
		{
			VectorN ret;
			ret[indUnit] = 1.0;
			return ret;
		}
		VectorN Normalized() const
		{
			VectorN ret;
			Real norm = NormL2();
			if (norm == 0.0)
				throw VectorDimensionError("VectorN::Normalized - cannot normalize zero vector", N, 0);
			for (int i = 0; i < N; ++i)
				ret._val[i] = _val[i] / norm;
			return ret;
		}
		////////////////////////            Standard stuff             ////////////////////////
		int  size() const { return N; }
		void clear() {
			for (int i = 0; i < N; ++i)
				_val[i] = Type{ 0 };
		}

		///////////////////////            Accessing elements             ///////////////////////
		inline Type& operator[](int n) { return _val[n]; }
		inline Type  operator[](int n) const { return _val[n]; }

		// checked access
		Type& at(int n) {
			if (n < 0 || n >= N)
				throw VectorDimensionError("VectorN::at - index out of bounds", N, n);
			else
				return _val[n];
		}
		Type  at(int n) const {
			if (n < 0 || n >= N)
				throw VectorDimensionError("VectorN::at - index out of bounds", N, n);
			else
				return _val[n];
		}

		///////////////////////            Testing equality                ///////////////////////
		bool operator==(const VectorN& b) const
		{
			for (int i = 0; i < size(); i++)
				if ((*this)[i] != b[i])
					return false;

			return true;
		}
		static bool AreEqual(const VectorN& a, const VectorN& b, Type eps = Defaults::VectorIsEqualTolerance)
		{
			return a.IsEqualTo(b, eps);
		}
		bool IsEqualTo(const VectorN& b, Type eps = Defaults::VectorIsEqualTolerance) const
		{
			for (int i = 0; i < N; i++)
			{
				if (Abs((*this)[i] - b[i]) > eps)
					return false;
			}
			return true;
		}
		bool operator!=(const VectorN& b) const
		{
			return !(*this == b);
		}

		bool IsNullVec() const
		{
			for (int i = 0; i < N; i++)
				if (_val[i] != 0.0)
					return false;

			return true;
		}
		
		/////////////////////             Arithmetic operators             ///////////////////////
		VectorN  operator-() const        // unary minus
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = Type{ -1 } *_val[i];
			return ret;
		}
		
		VectorN  operator+(const VectorN& b) const
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = _val[i] + b._val[i];
			return ret;
		}
		VectorN& operator+=(const VectorN& b)
		{
			for (int i = 0; i < N; i++)
				_val[i] += b._val[i];
			return *this;
		}
		VectorN  operator-(const VectorN& b) const
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = _val[i] - b._val[i];
			return ret;
		}
		VectorN& operator-=(const VectorN& b)
		{
			for (int i = 0; i < N; i++)
				_val[i] -= b._val[i];
			return *this;
		}

		VectorN  operator*(const Type &b) const
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = _val[i] * b;
			return ret;
		}
		VectorN& operator*=(const Type& b)
		{
			for (int i = 0; i < N; i++)
				_val[i] *= b;
			return *this;
		}
		VectorN  operator/(const Type &b) const
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = _val[i] / b;
			return ret;
		}		
		VectorN& operator/=(const Type& b)
		{
			for (int i = 0; i < N; i++)
				_val[i] /= b;
			return *this;
		}
		
		friend VectorN operator*(Type a, const VectorN<Type, N>& b)
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = a * b[i];
			return ret;
		}

		//////////////////////                 Operations                 ///////////////////////
		Real NormL1() const
		{
			Real norm{ 0.0 };
			for (int i = 0; i < N; i++)
				norm += Abs((*this)[i]);
			return norm;
		}
		Real NormL2() const
		{
			Real norm{ 0.0 };
			for (int i = 0; i < N; i++)
				norm += (*this)[i] * (*this)[i];
			return std::sqrt(norm);
		}
		Real NormLInf() const
		{
			Real norm{ 0.0 };
			for (int i = 0; i < N; i++)
				norm = std::max(norm, Abs((*this)[i]));
			return norm;
		}

		///////////////////////////               I/O                 ///////////////////////////
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
		std::ostream& PrintLine(std::ostream& stream, std::string msg, int width, int precision) const
		{
			stream << msg;
			Print(stream, width, precision);
			stream << std::endl;

			return stream;
		}

		std::string to_string(int width, int precision) const
		{
			std::stringstream str;

			Print(str, width, precision);

			return str.str();
		}
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