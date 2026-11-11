///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///  File:        GeometryPredicates.h                                                 ///
///  Description: Exact geometry determinants for finite binary coordinates           ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                         ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GEOMETRY_PREDICATES_H
#define MML_GEOMETRY_PREDICATES_H

#include <mml/base/Geometry/GeometryBase/GeometryPoints.h>
#include <mml/base/Geometry/GeometryBase/GeometryArithmetic.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>

namespace MML::detail {
	class GeometryExactArithmetic {
	public:
		static_assert(std::numeric_limits<Real>::radix == 2, "Exact predicates require binary floating-point Real");
		static constexpr int WordBits = 32;
		static constexpr int CoordinateWords = (std::numeric_limits<Real>::digits + WordBits - 1) / WordBits;
		static constexpr int MinimumCoordinateExponent = std::numeric_limits<Real>::min_exponent
			- std::numeric_limits<Real>::digits - CoordinateWords * WordBits;

		struct BinaryCoordinate {
			std::array<std::uint32_t, CoordinateWords> words{};
			int exponent = 0;
			bool negative = false;

			BinaryCoordinate() = default;
			explicit BinaryCoordinate(Real value) : negative(std::signbit(value)) {
				Real fraction = std::frexp(std::abs(value), &exponent);
				for (int i = CoordinateWords; i > 0; --i) {
					fraction = std::ldexp(fraction, WordBits);
					words[i - 1] = static_cast<std::uint32_t>(fraction);
					fraction -= words[i - 1];
				}
				exponent -= CoordinateWords * WordBits;
			}
		};

		template<std::size_t N>
		static void AddWord(std::array<std::uint32_t, N>& accumulator, std::size_t index, std::uint64_t carry) {
			while (carry != 0) {
				carry += accumulator[index];
				accumulator[index++] = static_cast<std::uint32_t>(carry);
				carry >>= WordBits;
			}
		}

		template<std::size_t N>
		static int Compare(const std::array<std::uint32_t, N>& positive, const std::array<std::uint32_t, N>& negative) {
			for (std::size_t i = N; i > 0; --i)
				if (positive[i - 1] != negative[i - 1])
					return positive[i - 1] > negative[i - 1] ? 1 : -1;
			return 0;
		}

		template<std::size_t Degree>
		class PolynomialAccumulator {
			static_assert(Degree >= 1 && Degree <= 5, "Geometry determinant degree must be between one and five");
			static constexpr int MinimumExponent = static_cast<int>(Degree) * MinimumCoordinateExponent;
			// One spare word covers carries from up to 360 signed determinant terms.
			static constexpr int Words = (static_cast<int>(Degree) * std::numeric_limits<Real>::max_exponent
				- MinimumExponent + WordBits - 1) / WordBits + 1;
			std::array<std::uint32_t, Words> positive{}, negative{};

		public:
			void AddProduct(const std::array<const BinaryCoordinate*, Degree>& factors, bool subtract) {
				for (const auto* factor : factors)
					if (std::all_of(factor->words.begin(), factor->words.end(), [](std::uint32_t word) { return word == 0; })) return;
				std::array<std::uint32_t, Degree * CoordinateWords> product{};
				std::copy(factors[0]->words.begin(), factors[0]->words.end(), product.begin());
				int count = CoordinateWords, exponent = factors[0]->exponent;
				bool negativeTerm = subtract != factors[0]->negative;
				for (std::size_t f = 1; f < Degree; ++f) {
					const auto& factor = *factors[f];
					std::array<std::uint32_t, Degree * CoordinateWords> next{};
					for (int i = 0; i < count; ++i) {
						std::uint64_t carry = 0;
						for (int j = 0; j < CoordinateWords; ++j) {
							const std::uint64_t term = std::uint64_t(product[i]) * factor.words[j] + next[i + j] + carry;
							next[i + j] = static_cast<std::uint32_t>(term);
							carry = term >> WordBits;
						}
						next[i + CoordinateWords] = static_cast<std::uint32_t>(carry);
					}
					count += CoordinateWords;
					product = next;
					exponent += factor.exponent;
					negativeTerm = negativeTerm != factor.negative;
				}
				auto& accumulator = negativeTerm ? negative : positive;
				const int shift = exponent - MinimumExponent;
				for (std::size_t i = 0; i < product.size(); ++i)
					AddWord(accumulator, static_cast<std::size_t>(shift / WordBits) + i,
						std::uint64_t(product[i]) << (shift % WordBits));
			}

			int Sign() const { return Compare(positive, negative); }
		};
	};

	class GeometryExact2D : private GeometryExactArithmetic {
		static constexpr int MinimumProductExponent = 2 * MinimumCoordinateExponent;
		// Three carry bits cover six products across Real's full binary exponent range.
		static constexpr int AccumulatorWords = (2 * std::numeric_limits<Real>::max_exponent
			- MinimumProductExponent + 3 + WordBits - 1) / WordBits;
		using Accumulator = std::array<std::uint32_t, AccumulatorWords>;

		static void AddProduct(Accumulator& positive, Accumulator& negative,
			const BinaryCoordinate& left, const BinaryCoordinate& right, bool subtract) {
			auto& accumulator = (left.negative != right.negative) != subtract ? negative : positive;
			const int shift = left.exponent + right.exponent - MinimumProductExponent;
			const int wordShift = shift / WordBits;
			const int bitShift = shift % WordBits;
			for (int i = 0; i < CoordinateWords; ++i) {
				for (int j = 0; j < CoordinateWords; ++j) {
					const auto product = std::uint64_t(left.words[i]) * right.words[j];
					const auto index = static_cast<std::size_t>(wordShift + i + j);
					AddWord(accumulator, index, std::uint64_t(static_cast<std::uint32_t>(product)) << bitShift);
					AddWord(accumulator, index + 1, (product >> WordBits) << bitShift);
				}
			}
		}

		static GeometryBinary Value(const Accumulator& positive, const Accumulator& negative) {
			const int sign = Compare(positive, negative);
			if (sign == 0) return {};
			const auto& larger = sign > 0 ? positive : negative;
			const auto& smaller = sign > 0 ? negative : positive;
			Accumulator magnitude{};
			std::uint64_t borrow = 0;
			for (std::size_t i = 0; i < magnitude.size(); ++i) {
				const std::uint64_t subtrahend = std::uint64_t(smaller[i]) + borrow;
				magnitude[i] = static_cast<std::uint32_t>(std::uint64_t(larger[i]) - subtrahend);
				borrow = std::uint64_t(larger[i]) < subtrahend;
			}
			int top = AccumulatorWords - 1;
			while (magnitude[top] == 0) --top;
			Real fraction = 0;
			for (int i = top; i >= 0 && top - i <= CoordinateWords; --i)
				fraction += std::ldexp(static_cast<Real>(magnitude[i]), WordBits * (i - top));
			return GeometryBinary::Scaled(sign * fraction, MinimumProductExponent + WordBits * top);
		}

		static GeometryBinary ExactDeterminant(const Point2Cartesian& a, const Point2Cartesian& b, const Point2Cartesian& c) {
			const std::array<BinaryCoordinate, 3> x{
				BinaryCoordinate(a.X()), BinaryCoordinate(b.X()), BinaryCoordinate(c.X())};
			const std::array<BinaryCoordinate, 3> y{
				BinaryCoordinate(a.Y()), BinaryCoordinate(b.Y()), BinaryCoordinate(c.Y())};
			Accumulator positive{}, negative{};
			// Six signed coordinate products avoid rounded/overflowing coordinate differences.
			// Fixed binary accumulators also retain products below Real's subnormal range.
			for (std::size_t i = 0; i < x.size(); ++i) {
				AddProduct(positive, negative, x[i], y[(i + 1) % 3], false);
				AddProduct(positive, negative, x[i], y[(i + 2) % 3], true);
			}
			return Value(positive, negative);
		}

	public:
		// Magnitudes are rounded to Real precision, but retain a separate exponent.
		static GeometryBinary Determinant(const Point2Cartesian& a, const Point2Cartesian& b, const Point2Cartesian& c) {
			ValidateGeometryCoordinates({a.X(), a.Y(), b.X(), b.Y(), c.X(), c.Y()}, "Planar determinant");
			return ExactDeterminant(a, b, c);
		}

		static int Sign(const Point2Cartesian& a, const Point2Cartesian& b, const Point2Cartesian& c) {
			ValidateGeometryCoordinates({a.X(), a.Y(), b.X(), b.Y(), c.X(), c.Y()}, "Orientation2D");
			const Real acx = a.X() - c.X(), acy = a.Y() - c.Y();
			const Real bcx = b.X() - c.X(), bcy = b.Y() - c.Y();
			const Real left = acx * bcy, right = acy * bcx;
			const Real determinant = left - right;
			const Real determinantSum = std::abs(left) + std::abs(right);
			const bool leftSafe = std::isnormal(left) || acx == 0 || bcy == 0;
			const bool rightSafe = std::isnormal(right) || acy == 0 || bcx == 0;
			// Shewchuk (1997), Adaptive Precision Floating-Point Arithmetic: orient2d filter.
			if (leftSafe && rightSafe && std::isnormal(determinant) && std::isfinite(determinantSum)
				&& std::abs(determinant) > Precision::GeometryOrientationErrorBound * determinantSum)
				return determinant > 0 ? 1 : -1;
			const Real fraction = ExactDeterminant(a, b, c).fraction;
			return fraction > 0 ? 1 : fraction < 0 ? -1 : 0;
		}

		static int InCircleSign(const Point2Cartesian& a, const Point2Cartesian& b,
			const Point2Cartesian& c, const Point2Cartesian& d) {
			PolynomialAccumulator<4> accumulator;
			const std::array<BinaryCoordinate, 4> x{
				BinaryCoordinate(a.X()), BinaryCoordinate(b.X()), BinaryCoordinate(c.X()), BinaryCoordinate(d.X())};
			const std::array<BinaryCoordinate, 4> y{
				BinaryCoordinate(a.Y()), BinaryCoordinate(b.Y()), BinaryCoordinate(c.Y()), BinaryCoordinate(d.Y())};
			const auto add = [&](const BinaryCoordinate& lift, const BinaryCoordinate& first,
				const BinaryCoordinate& second, bool subtract) {
				accumulator.AddProduct({&lift, &lift, &first, &second}, subtract);
			};
			// Expand det[x,y,x*x+y*y,1] along its lifted column, without
			// coordinate subtraction or floating products (Shewchuk's incircle determinant).
			for (int lifted = 0; lifted < 4; ++lifted) {
				std::array<int, 3> others{};
				int count = 0;
				for (int i = 0; i < 4; ++i)
					if (i != lifted) others[count++] = i;
				for (const auto& lift : {x[lifted], y[lifted]}) {
					for (int j = 0; j < 3; ++j) {
						add(lift, x[others[j]], y[others[(j + 1) % 3]], lifted % 2 != 0);
						add(lift, x[others[j]], y[others[(j + 2) % 3]], lifted % 2 == 0);
					}
				}
			}
			return accumulator.Sign();
		}
	};

	class GeometryExact3D : private GeometryExactArithmetic {
		template<std::size_t N>
		static int DeterminantSign(const std::array<Point3Cartesian, N>& points) {
			std::array<std::array<BinaryCoordinate, N>, 3> coordinates{};
			std::array<int, N> permutation{};
			for (std::size_t i = 0; i < N; ++i) {
				coordinates[0][i] = BinaryCoordinate(points[i].X());
				coordinates[1][i] = BinaryCoordinate(points[i].Y());
				coordinates[2][i] = BinaryCoordinate(points[i].Z());
				permutation[i] = static_cast<int>(i);
			}
			PolynomialAccumulator<N == 4 ? 3 : 5> accumulator;
			// Shewchuk (1997): determinants [x,y,z,1] and [x,y,z,x*x+y*y+z*z,1].
			// Original-coordinate integer products retain even mixed-scale tails.
			do {
				bool subtract = false;
				for (std::size_t i = 0; i < N; ++i)
					for (std::size_t j = i + 1; j < N; ++j)
						if (permutation[i] > permutation[j]) subtract = !subtract;
				const auto* x = &coordinates[0][permutation[0]];
				const auto* y = &coordinates[1][permutation[1]];
				const auto* z = &coordinates[2][permutation[2]];
				if constexpr (N == 4) accumulator.AddProduct({x, y, z}, subtract);
				else
					for (const auto& axis : coordinates) {
						const auto* lift = &axis[permutation[3]];
						accumulator.AddProduct({x, y, z, lift, lift}, subtract);
					}
			} while (std::next_permutation(permutation.begin(), permutation.end()));
			return accumulator.Sign();
		}

	public:
		static int OrientationSign(const std::array<Point3Cartesian, 4>& points) { return DeterminantSign(points); }
		static int InSphereSign(const std::array<Point3Cartesian, 5>& points) { return DeterminantSign(points); }
	};

	class GeometryOrientation2D {
	public:
		static int Sign(const Point2Cartesian& a, const Point2Cartesian& b, const Point2Cartesian& c) {
			return GeometryExact2D::Sign(a, b, c);
		}
	};
}

#endif // MML_GEOMETRY_PREDICATES_H
