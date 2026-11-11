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
#include <algorithm>
#include <array>
#include <utility>

namespace MML
{
	/// @brief 4D Minkowski spacetime vector with T (time), X, Y, Z components.
	class Vector4Minkowski : public VectorN<Real, 4>
	{
		struct ScaledProduct {
			Real mantissa = 0;
			int exponent = 0;
		};

		bool HasFiniteComponents() const
		{
			return std::isfinite(T()) && std::isfinite(X()) && std::isfinite(Y()) && std::isfinite(Z());
		}

		static ScaledProduct ComputeProduct(const Vector4Minkowski& a, const Vector4Minkowski& b)
		{
			if (!a.HasFiniteComponents() || !b.HasFiniteComponents())
				throw NumericInputError("Minkowski products require finite components");
			struct Term { Real high, low; int exponent; };
			std::array<Term, 4> terms{};
			size_t count = 0;
			for (int i = 0; i < 4; ++i) {
				if (a[i] == 0 || b[i] == 0)
					continue;
				int firstExponent, secondExponent;
				const Real first = std::frexp(a[i], &firstExponent);
				const Real second = std::frexp(b[i], &secondExponent);
				Real high = first * second;
				Real low = std::fma(first, second, -high);
				if (i == 0) { high = -high; low = -low; }
				terms[count++] = {high, low, firstExponent + secondExponent};
			}
			// Bound the four-term sort explicitly to avoid GCC's inlined std::sort bounds diagnostic.
			for (size_t i = 1; i < count; ++i) {
				const Term term = terms[i];
				size_t position = i;
				while (position > 0 && term.exponent > terms[position - 1].exponent) {
					terms[position] = terms[position - 1];
					--position;
				}
				terms[position] = term;
			}
			const auto twoSum = [](Real first, Real second) {
				const Real sum = first + second;
				const Real part = sum - first;
				return std::pair<Real, Real>{sum, (first - (sum - part)) + (second - part)};
			};
			Real high = 0, low = 0;
			int exponent = 0;
			// Sum largest products first, retaining FMA product errors. Rebase after cancellation
			// so a small surviving spatial term is not lost beneath an enormous lightlike pair.
			for (size_t i = 0; i < count; ++i) {
				const auto& term = terms[i];
				if (high == 0 && low == 0) {
					high = term.high;
					low = term.low;
					exponent = term.exponent;
				} else {
					const int commonExponent = std::max(exponent, term.exponent);
					high = std::scalbn(high, exponent - commonExponent);
					low = std::scalbn(low, exponent - commonExponent);
					const auto sum = twoSum(high, std::scalbn(term.high, term.exponent - commonExponent));
					const auto normalized = twoSum(sum.first, low + sum.second +
						std::scalbn(term.low, term.exponent - commonExponent));
					high = normalized.first;
					low = normalized.second;
					exponent = commonExponent;
				}
				if (high != 0) {
					int shift;
					high = std::frexp(high, &shift);
					low = std::scalbn(low, -shift);
					exponent += shift;
				}
			}
			return {high + low, exponent};
		}

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
			const auto product = ComputeProduct(a, b);
			return std::scalbn(product.mantissa, product.exponent);
		}

		/// @brief Computes Minkowski norm (spacetime interval).
		Real Norm() const
		{
			auto interval = ComputeProduct(*this, *this);
			if (interval.mantissa == 0)
				return 0;
			const bool timelike = interval.mantissa < 0;
			interval.mantissa = std::abs(interval.mantissa);
			if (interval.exponent % 2 != 0) {
				interval.mantissa *= Real{2};
				--interval.exponent;
			}
			const Real norm = std::scalbn(std::sqrt(interval.mantissa), interval.exponent / 2);
			return timelike ? norm : -norm;
		}

		/// @brief Computes Minkowski distance to another 4-vector.
		/// @param b Other 4-vector
		Real Distance(const Vector4Minkowski& b) const
		{
			Vector4Minkowski diff{ T() - b.T(), X() - b.X(), Y() - b.Y(), Z() - b.Z() };
			return diff.Norm();
		}
		
		/// @brief Checks if 4-vector is timelike (inside light cone).
		bool IsTimelike() const
		{
			return HasFiniteComponents() && ComputeProduct(*this, *this).mantissa < 0;
		}
		/// @brief Checks if 4-vector is spacelike (outside light cone).
		bool IsSpacelike() const
		{
			return HasFiniteComponents() && ComputeProduct(*this, *this).mantissa > 0;
		}
		/// @brief Checks if 4-vector is lightlike (on light cone).
		bool IsLightlike() const
		{
			return HasFiniteComponents() && ComputeProduct(*this, *this).mantissa == 0;
		}

		/// @deprecated Use IsTimelike() (CONTRACT_POLICY §11: mathematical predicates are PascalCase).
		[[deprecated("use IsTimelike()")]]
		bool isTimelike() const { return IsTimelike(); }
		/// @deprecated Use IsSpacelike().
		[[deprecated("use IsSpacelike()")]]
		bool isSpacelike() const { return IsSpacelike(); }
		/// @deprecated Use IsLightlike().
		[[deprecated("use IsLightlike()")]]
		bool isLightlike() const { return IsLightlike(); }
	};

	typedef Vector4Minkowski Vec4Mink;
}

#endif