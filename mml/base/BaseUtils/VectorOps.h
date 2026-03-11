///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        VectorOps.h                                                         ///
///  Description: Vector operation utilities - scalar products, angles, projections   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_VECTOR_OPS_H
#define MML_VECTOR_OPS_H

#include <cmath>
#include <complex>

#include "mml/MMLBase.h"
#include "mml/MMLExceptions.h"
#include "mml/base/Vector/Vector.h"
#include "mml/base/Vector/VectorN.h"
#include "mml/base/Matrix/Matrix.h"

namespace MML
{
	namespace Utils
	{
		// ============================================================================
		// Vector<Real> Operations
		// ============================================================================

		/// @brief Computes scalar (dot) product of two real vectors.
		/// @param a First vector
		/// @param b Second vector
		/// @return Sum of element-wise products
		/// @throws VectorDimensionError if vectors have different sizes
		static Real ScalarProduct(const Vector<Real>& a, const Vector<Real>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("ScalarProduct - must be same dim", a.size(), b.size());

			Real ret = 0;
			for (int i = 0; i < a.size(); i++)
				ret += a[i] * b[i];

			return ret;
		}

		/// @brief Computes scalar (dot) product of two complex vectors.
		/// @details Computes sum of a[i] * conj(b[i]) for proper inner product.
		/// @param a First complex vector
		/// @param b Second complex vector (conjugated in product)
		/// @return Complex inner product
		/// @throws VectorDimensionError if vectors have different sizes
		static Complex ScalarProduct(const Vector<Complex>& a, const Vector<Complex>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("ScalarProduct - must be same dim", a.size(), b.size());

			Complex ret = 0;
			for (int i = 0; i < a.size(); i++)
				ret += a[i] * std::conj(b[i]);

			return ret;
		}

		/// @brief Computes angle between two real vectors.
		/// @param a First vector
		/// @param b Second vector
		/// @return Angle in radians [0, π]
		/// @throws VectorDimensionError if vectors have different sizes
		static Real VectorsAngle(const Vector<Real>& a, const Vector<Real>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("Vector::AngleToVector - vectors must be equal size", a.size(), b.size());

			Real cosAngle = ScalarProduct(a, b) / (a.NormL2() * b.NormL2());
			cosAngle = std::clamp(cosAngle, Real(-1), Real(1));
			return std::acos(cosAngle);
		}

		/// @brief Computes projection of vector onto direction of another vector.
		/// @param orig Vector to project
		/// @param b Direction vector to project onto
		/// @return Component of orig parallel to b
		static Vector<Real> VectorProjectionParallelTo(const Vector<Real>& orig, const Vector<Real>& b)
		{
			return ScalarProduct(orig, b) / ScalarProduct(b, b) * b;
		}

		/// @brief Computes rejection of vector from direction of another vector.
		/// @param orig Vector to project
		/// @param b Direction vector to reject from
		/// @return Component of orig perpendicular to b
		static Vector<Real> VectorProjectionPerpendicularTo(const Vector<Real>& orig, const Vector<Real>& b)
		{
			return orig - VectorProjectionParallelTo(orig, b);
		}

		// ============================================================================
		// VectorN<Real/Complex, N> Operations
		// ============================================================================

		/// @brief Computes scalar (dot) product of two fixed-size real vectors.
		/// @tparam N Vector dimension
		/// @param a First vector
		/// @param b Second vector
		/// @return Sum of element-wise products
		template<int N>
		static Real ScalarProduct(const VectorN<Real, N> &a, const VectorN<Real, N> &b)
		{
			Real ret = 0;
			for (int i = 0; i < N; i++)
				ret += a[i] * b[i];
			return ret;
		}

		/// @brief Computes scalar (dot) product of two fixed-size complex vectors.
		/// @details Computes sum of a[i] * conj(b[i]) for proper inner product.
		/// @tparam N Vector dimension
		/// @param a First complex vector
		/// @param b Second complex vector (conjugated in product)
		/// @return Complex inner product
		template<int N> 
		static Complex ScalarProduct(const VectorN<Complex, N>& a, const VectorN<Complex, N>& b)
		{
			Complex ret = 0;
			for (int i = 0; i < N; i++)
				ret += a[i] * std::conj(b[i]);
			return ret;
		}

		/// @brief Computes angle between two fixed-size real vectors.
		/// @tparam N Vector dimension
		/// @param a First vector
		/// @param b Second vector
		/// @return Angle in radians [0, π]
		template<int N> 
		static Real VectorsAngle(const VectorN<Real, N> &a, const VectorN<Real, N> &b)
		{
			Real cosAngle = ScalarProduct(a, b) / (a.NormL2() * b.NormL2());
			cosAngle = std::clamp(cosAngle, Real(-1), Real(1));
			return std::acos(cosAngle);
		}

		// ============================================================================
		// Outer Product
		// ============================================================================

		/// @brief Computes outer product of two vectors.
		/// @details Result is matrix C where C[i][j] = a[i] * b[j].
		/// @tparam Type Element type
		/// @param a First vector (determines rows)
		/// @param b Second vector (determines columns)
		/// @return Matrix of size a.size() x b.size()
		template<class Type> 
		static Matrix<Type> OuterProduct(const Vector<Type>& a, const Vector<Type>& b)
		{
			Matrix<Type> ret(a.size(), b.size());

			for (int i = 0; i < a.size(); i++)
				for (int j = 0; j < b.size(); j++)
					ret[i][j] = a[i] * b[j];

			return ret;
		}

	} // namespace Utils
} // namespace MML

#endif // MML_VECTOR_OPS_H
