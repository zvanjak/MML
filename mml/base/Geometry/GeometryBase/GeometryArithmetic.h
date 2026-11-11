///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///  File:        GeometryArithmetic.h                                                ///
///  Description: Scaled direction, cross-product and midpoint arithmetic              ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                         ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_GEOMETRY_ARITHMETIC_H
#define MML_GEOMETRY_ARITHMETIC_H

#include <mml/base/Geometry/GeometryBase/GeometryPoints.h>
#include <mml/base/Vector/VectorN.h>

#include <array>
#include <numeric>

namespace MML::detail {
	template<int N>
	VectorN<Real, N> GeometryUnitDirection(const VectorN<Real, N>& value,
		const char* nonfinite = "Geometry - non-finite direction",
		const char* zero = "Geometry - undefined direction at a zero vector") {
		Real scale = 0;
		for (int index = 0; index < N; ++index) {
			if (!std::isfinite(value[index])) throw GeometryError(nonfinite);
			scale = std::max(scale, std::abs(value[index]));
		}
		if (scale == 0) throw GeometryError(zero);
		const auto scaled = value / scale;
		return scaled / scaled.NormL2();
	}

	struct GeometryBinary {
		Real fraction = 0;
		int exponent = 0;

		GeometryBinary() = default;
		explicit GeometryBinary(Real value) {
			ValidateGeometryCoordinates({value}, "Geometry arithmetic");
			fraction = std::frexp(value, &exponent);
		}

		static GeometryBinary DifferenceOfProducts(const GeometryBinary& a, const GeometryBinary& b,
			const GeometryBinary& c, const GeometryBinary& d) {
			const bool firstZero = a.fraction == 0 || b.fraction == 0;
			const bool secondZero = c.fraction == 0 || d.fraction == 0;
			if (firstZero && secondZero) return {};
			const int firstExponent = a.exponent + b.exponent;
			const int secondExponent = c.exponent + d.exponent;
			const int exponent = firstZero ? secondExponent : secondZero ? firstExponent
				: std::max(firstExponent, secondExponent);
			const Real first = firstZero ? 0 : std::ldexp(b.fraction, firstExponent - exponent);
			const Real second = secondZero ? 0 : std::ldexp(d.fraction, secondExponent - exponent);
			const Real product = c.fraction * second;
			// Compensated difference of products (Kahan): retain cancellation below one product ulp.
			const Real difference = std::fma(a.fraction, first, -product)
				+ std::fma(-c.fraction, second, product);
			GeometryBinary result;
			int shift = 0;
			result.fraction = std::frexp(difference, &shift);
			if (result.fraction != 0) result.exponent = exponent + shift;
			return result;
		}

		static GeometryBinary Difference(Real a, Real b) {
			return DifferenceOfProducts(GeometryBinary(a), GeometryBinary(1), GeometryBinary(b), GeometryBinary(1));
		}

		static GeometryBinary Add(const GeometryBinary& a, GeometryBinary b) {
			b.fraction = -b.fraction;
			return DifferenceOfProducts(a, GeometryBinary(1), b, GeometryBinary(1));
		}

		Real Value() const { return std::ldexp(fraction, exponent); }

		static GeometryBinary Scaled(Real fraction, int exponent) {
			GeometryBinary result;
			int shift = 0;
			result.fraction = std::frexp(fraction, &shift);
			if (result.fraction != 0) result.exponent = exponent + shift;
			return result;
		}

		static GeometryBinary Product(const GeometryBinary& a, const GeometryBinary& b) {
			return Scaled(a.fraction * b.fraction, a.exponent + b.exponent);
		}

		static std::array<GeometryBinary, 2> ProductParts(const GeometryBinary& a, const GeometryBinary& b) {
			const Real product = a.fraction * b.fraction;
			return {Scaled(product, a.exponent + b.exponent),
				Scaled(std::fma(a.fraction, b.fraction, -product), a.exponent + b.exponent)};
		}

		static GeometryBinary Quotient(const GeometryBinary& a, const GeometryBinary& b) {
			if (b.fraction == 0) throw GeometryError("Geometry arithmetic - zero divisor");
			return Scaled(a.fraction / b.fraction, a.exponent - b.exponent);
		}

		GeometryBinary Sqrt() const {
			if (fraction < 0) throw GeometryError("Geometry arithmetic - negative squared length");
			const int half = exponent / 2;
			return Scaled(std::sqrt(std::ldexp(fraction, exponent - 2 * half)), half);
		}

		Real CheckedValue() const {
			const Real result = Value();
			if (!std::isfinite(result) || (result == 0 && fraction != 0))
				throw GeometryError("Geometry query - result is not representable");
			return result;
		}

		bool AbsLessThan(Real tolerance) const {
			if (tolerance == 0) return false;
			if (fraction == 0) return true;
			const GeometryBinary bound(tolerance);
			return exponent < bound.exponent ||
				(exponent == bound.exponent && std::abs(fraction) < bound.fraction);
		}

		static std::array<GeometryBinary, 2> TwoSum(GeometryBinary a, GeometryBinary b) {
			if (a.fraction == 0) return {b, {}};
			if (b.fraction == 0) return {a, {}};
			if (a.exponent < b.exponent) std::swap(a, b);
			const int gap = a.exponent - b.exponent;
			if (gap > std::numeric_limits<Real>::digits + 1) return {a, b};
			const Real small = std::ldexp(b.fraction, -gap);
			const Real sum = a.fraction + small;
			const Real virtualSmall = sum - a.fraction;
			const Real error = (a.fraction - (sum - virtualSmall)) + (small - virtualSmall);
			return {Scaled(sum, a.exponent), Scaled(error, a.exponent)};
		}

		template<std::size_t N>
		static GeometryBinary SumProducts(const std::array<GeometryBinary, N>& a,
			const std::array<GeometryBinary, N>& b, GeometryBinary* remainder = nullptr) {
			// Error-free products and growing expansions (Shewchuk, 1997), with
			// separate binary exponents so cancellation cannot overflow or erase tiny terms.
			std::array<GeometryBinary, 2 * N> expansion{};
			std::size_t count = 0;
			for (std::size_t index = 0; index < N; ++index) {
				const auto parts = ProductParts(a[index], b[index]);
				for (auto term : {parts[1], parts[0]}) {
					std::array<GeometryBinary, 2 * N> next{};
					std::size_t length = 0;
					for (std::size_t j = 0; j < count; ++j) {
						const auto sum = TwoSum(term, expansion[j]);
						if (sum[1].fraction != 0) next[length++] = sum[1];
						term = sum[0];
					}
					if (term.fraction != 0) next[length++] = term;
					expansion = next;
					count = length;
				}
			}
			GeometryBinary low;
			for (std::size_t index = 0; index + 1 < count; ++index) low = Add(low, expansion[index]);
			const auto result = count == 0 ? std::array<GeometryBinary, 2>{} : TwoSum(expansion[count - 1], low);
			if (remainder) *remainder = result[1];
			return result[0];
		}
	};

	inline GeometryBinary GeometryPlaneEvaluation(const Pnt3Cart& point,
		const VectorN<Real, 3>& normal, Real offset, GeometryBinary* remainder = nullptr) {
		using B = GeometryBinary;
		return B::SumProducts<4>({B(normal[0]), B(normal[1]), B(normal[2]), B(offset)},
			{B(point.X()), B(point.Y()), B(point.Z()), B(1)}, remainder);
	}

	inline Real GeometryPlaneOffset(const Pnt3Cart& point, const VectorN<Real, 3>& normal) {
		const auto sum = GeometryPlaneEvaluation(point, normal, 0);
		const Real offset = -sum.Value();
		if (!std::isfinite(offset) || (offset == 0 && sum.fraction != 0))
			throw GeometryError("Plane3D ctor - offset is not representable");
		return offset;
	}

	class GeometryCrossProduct {
		VectorN<Real, 3> _scaled;
		int _exponent = 0;
	public:
		GeometryCrossProduct(const Pnt3Cart& a, const Pnt3Cart& b, const Pnt3Cart& c, const Pnt3Cart& d) {
			const std::array<GeometryBinary, 3> first{GeometryBinary::Difference(b.X(), a.X()),
				GeometryBinary::Difference(b.Y(), a.Y()), GeometryBinary::Difference(b.Z(), a.Z())};
			const std::array<GeometryBinary, 3> second{GeometryBinary::Difference(d.X(), c.X()),
				GeometryBinary::Difference(d.Y(), c.Y()), GeometryBinary::Difference(d.Z(), c.Z())};
			std::array<GeometryBinary, 3> cross;
			bool nonzero = false;
			for (int index = 0; index < 3; ++index) {
				const int next = (index + 1) % 3, last = (index + 2) % 3;
				cross[index] = GeometryBinary::DifferenceOfProducts(first[next], second[last], first[last], second[next]);
				if (cross[index].fraction != 0) {
					_exponent = nonzero ? std::max(_exponent, cross[index].exponent) : cross[index].exponent;
					nonzero = true;
				}
			}
			for (int index = 0; index < 3; ++index)
				_scaled[index] = cross[index].fraction == 0 ? 0
					: std::ldexp(cross[index].fraction, cross[index].exponent - _exponent);
		}

		VectorN<Real, 3> Unit() const { return GeometryUnitDirection(_scaled); }

		Real HalfMagnitude() const {
			const Real area = std::ldexp(_scaled.NormL2() / 2, _exponent);
			if (!std::isfinite(area) || area == 0)
				throw GeometryError("RectSurface3D::getArea - area is not representable or is zero");
			return area;
		}
	};

	inline Pnt3Cart GeometryMidpoint(const Pnt3Cart& a, const Pnt3Cart& b) {
		return {std::midpoint(a.X(), b.X()), std::midpoint(a.Y(), b.Y()), std::midpoint(a.Z(), b.Z())};
	}
}

#endif // MML_GEOMETRY_ARITHMETIC_H
