///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        VectorTypes4D.h                                                     ///
///  Description: 4D Minkowski spacetime vector type and alias                        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_VECTOR_TYPES_4D_H
#define MML_VECTOR_TYPES_4D_H

#include <mml/MMLBase.h>

#include <mml/base/Vector/VectorN.h>

#include <cmath>

namespace MML
{
	/// @brief 4D Minkowski spacetime vector with T (time), X, Y, Z components.
	class Vector4Minkowski : public VectorN<Real, 4>
	{
	public:
		/// @brief Returns T component (const).
		Real  T() const { return val(0); }
		/// @brief Returns T component (non-const).
		Real& T()				{ return val(0); }
		/// @brief Returns X component (const).
		Real  X() const { return val(1); }
		/// @brief Returns X component (non-const).
		Real& X()				{ return val(1); }
		/// @brief Returns Y component (const).
		Real  Y() const { return val(2); }
		/// @brief Returns Y component (non-const).
		Real& Y()				{ return val(2); }
		/// @brief Returns Z component (const).
		Real  Z() const { return val(3); }
		/// @brief Returns Z component (non-const).
		Real& Z()				{ return val(3); }

		/// @brief Default constructor (zero 4-vector).
		Vector4Minkowski() : VectorN<Real, 4>{ 0.0, 0.0, 0.0, 0.0 } {}
		/// @brief Constructs from initializer list.
		/// @param list Initializer list of values
		Vector4Minkowski(std::initializer_list<Real> list) : VectorN<Real, 4>(list) { }

		/// @brief Minkowski scalar product using the library signature (-,+,+,+).
		/// @param a First 4-vector
		/// @param b Second 4-vector
		friend Real ScalarProduct(const Vector4Minkowski& a, const Vector4Minkowski& b)
		{
			return -a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
		}

		/// @brief Computes Minkowski norm (spacetime interval).
		Real Norm() const
		{
			if(T() == 0.0 && X() == 0.0 && Y() == 0.0 && Z() == 0.0)
				return 0.0;

			Real interval_squared = ScalarProduct(*this, *this);
			if (isTimelike())
				return sqrt(-interval_squared);
			else if (isSpacelike())
				return -sqrt(interval_squared);
			else 
				return 0.0; // Lightlike vectors have zero Minkowski norm
		}

		/// @brief Computes Minkowski distance to another 4-vector.
		/// @param b Other 4-vector
		Real Distance(const Vector4Minkowski& b) const
		{
			Vector4Minkowski diff{ T() - b.T(), X() - b.X(), Y() - b.Y(), Z() - b.Z() };
			return diff.Norm();
		}
		
		/// @brief Checks if 4-vector is timelike (inside light cone).
		bool isTimelike() const
		{
			return ScalarProduct(*this, *this) < 0;
		}
		/// @brief Checks if 4-vector is spacelike (outside light cone).
		bool isSpacelike() const
		{
			return ScalarProduct(*this, *this) > 0;
		}
		/// @brief Checks if 4-vector is lightlike (on light cone).
		bool isLightlike() const
		{
			return ScalarProduct(*this, *this) == 0;
		}
	};

	typedef Vector4Minkowski Vec4Mink;
}

#endif