///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        RobustPredicates.h                                                  ///
///  Description: Adaptive exact-sign orientation and in-circle/in-sphere predicates ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COMP_GEOMETRY_ROBUST_PREDICATES_H
#define MML_COMP_GEOMETRY_ROBUST_PREDICATES_H

#include <mml/MMLBase.h>
#include <mml/base/Geometry/Geometry2D.h>
#include <mml/base/Geometry/Geometry3D.h>

#include <array>
#include <cmath>
#include <limits>
#include <utility>

namespace MML::CompGeometry {

	/// Adaptive predicates use a fast error-bounded evaluation and fall back to
	/// floating-point expansions whose sum is the exact value for the input values.
	class RobustPredicates {
		template<std::size_t Capacity>
		struct Expansion {
			std::array<Real, Capacity> components{};
			std::size_t size = 0;
		};

		static int Sign(Real value) {
			return (value > REAL(0.0)) - (value < REAL(0.0));
		}

		static void TwoSum(Real a, Real b, Real& sum, Real& error) {
			sum = a + b;
			const Real bVirtual = sum - a;
			const Real aVirtual = sum - bVirtual;
			error = (a - aVirtual) + (b - bVirtual);
		}

		static Expansion<2> Difference(Real a, Real b) {
			const Real head = a - b;
			const Real bVirtual = a - head;
			const Real tail = (a - (head + bVirtual)) + (bVirtual - b);
			Expansion<2> result;
			if (tail != REAL(0.0)) result.components[result.size++] = tail;
			if (head != REAL(0.0) || result.size == 0) result.components[result.size++] = head;
			return result;
		}

		template<std::size_t Capacity>
		static void Grow(const Expansion<Capacity>& expansion, Real value,
			Expansion<Capacity>& result) {
			result.size = 0;
			Real accumulator = value;
			for (std::size_t i = 0; i < expansion.size; ++i) {
				Real sum, error;
				TwoSum(accumulator, expansion.components[i], sum, error);
				if (error != REAL(0.0)) result.components[result.size++] = error;
				accumulator = sum;
			}
			if (accumulator != REAL(0.0) || result.size == 0)
				result.components[result.size++] = accumulator;
		}

		template<std::size_t LeftCapacity, std::size_t RightCapacity>
		static Expansion<LeftCapacity + RightCapacity> Add(
			const Expansion<LeftCapacity>& left, const Expansion<RightCapacity>& right) {
			Expansion<LeftCapacity + RightCapacity> result;
			for (std::size_t i = 0; i < left.size; ++i)
				result.components[result.size++] = left.components[i];
			Expansion<LeftCapacity + RightCapacity> scratch;
			for (std::size_t i = 0; i < right.size; ++i) {
				Grow(result, right.components[i], scratch);
				std::swap(result, scratch);
			}
			return result;
		}

		template<std::size_t Capacity>
		static Expansion<Capacity> Negate(Expansion<Capacity> value) {
			for (std::size_t i = 0; i < value.size; ++i)
				value.components[i] = -value.components[i];
			return value;
		}

		template<std::size_t Capacity>
		static Expansion<2 * Capacity> Scale(const Expansion<Capacity>& expansion, Real factor) {
			Expansion<2 * Capacity> result;
			result.components[result.size++] = REAL(0.0);
			Expansion<2 * Capacity> scratch;
			for (std::size_t i = 0; i < expansion.size; ++i) {
				const Real product = expansion.components[i] * factor;
				const Real error = std::fma(expansion.components[i], factor, -product);
				if (error != REAL(0.0)) {
					Grow(result, error, scratch);
					std::swap(result, scratch);
				}
				Grow(result, product, scratch);
				std::swap(result, scratch);
			}
			return result;
		}

		template<std::size_t LeftCapacity, std::size_t RightCapacity>
		static Expansion<2 * LeftCapacity * RightCapacity> Multiply(
			const Expansion<LeftCapacity>& left, const Expansion<RightCapacity>& right) {
			Expansion<2 * LeftCapacity * RightCapacity> result;
			result.components[result.size++] = REAL(0.0);
			Expansion<2 * LeftCapacity * RightCapacity> scratch;
			for (std::size_t i = 0; i < left.size; ++i) {
				for (std::size_t j = 0; j < right.size; ++j) {
					const Real product = left.components[i] * right.components[j];
					const Real error = std::fma(left.components[i], right.components[j], -product);
					if (error != REAL(0.0)) {
						Grow(result, error, scratch);
						std::swap(result, scratch);
					}
					Grow(result, product, scratch);
					std::swap(result, scratch);
				}
			}
			return result;
		}

		template<std::size_t Capacity>
		static int ExpansionSign(const Expansion<Capacity>& expansion) {
			for (std::size_t i = expansion.size; i > 0; --i)
				if (expansion.components[i - 1] != REAL(0.0)) return Sign(expansion.components[i - 1]);
			return 0;
		}

		template<std::size_t Capacity>
		static Expansion<4 * Capacity * Capacity> ExactCross(
			const Expansion<Capacity>& ax, const Expansion<Capacity>& ay,
			const Expansion<Capacity>& bx, const Expansion<Capacity>& by) {
			return Add(Multiply(ax, by), Negate(Multiply(ay, bx)));
		}

		static Expansion<1> Scalar(Real value) {
			Expansion<1> result;
			result.components[result.size++] = value;
			return result;
		}

		template<std::size_t Capacity>
		static void Append(Expansion<Capacity>& result, const Real value,
			Expansion<Capacity>& scratch) {
			Grow(result, value, scratch);
			std::swap(result, scratch);
		}

	public:
		/// Sign of the oriented area: +1 counter-clockwise, -1 clockwise, 0 collinear.
		static int Orientation2D(const Point2Cartesian& a, const Point2Cartesian& b,
			const Point2Cartesian& c) {
			const Real acx = a.X() - c.X();
			const Real bcx = b.X() - c.X();
			const Real acy = a.Y() - c.Y();
			const Real bcy = b.Y() - c.Y();
			const Real determinant = acx * bcy - acy * bcx;
			const Real determinantSum = std::abs(acx * bcy) + std::abs(acy * bcx);
			const Real epsilon = std::numeric_limits<Real>::epsilon();
			if (std::abs(determinant) > (REAL(3.0) + REAL(16.0) * epsilon) * epsilon * determinantSum)
				return Sign(determinant);

			return ExpansionSign(ExactCross(
				Difference(a.X(), c.X()), Difference(a.Y(), c.Y()),
				Difference(b.X(), c.X()), Difference(b.Y(), c.Y())));
		}

		/// Sign of the in-circle determinant for oriented triangle (a,b,c).
		/// For counter-clockwise (a,b,c), +1 means d is inside and 0 cocircular.
		static int InCircle2D(const Point2Cartesian& a, const Point2Cartesian& b,
			const Point2Cartesian& c, const Point2Cartesian& d) {
			const Real adx = a.X() - d.X(), ady = a.Y() - d.Y();
			const Real bdx = b.X() - d.X(), bdy = b.Y() - d.Y();
			const Real cdx = c.X() - d.X(), cdy = c.Y() - d.Y();
			const Real abdet = adx * bdy - bdx * ady;
			const Real bcdet = bdx * cdy - cdx * bdy;
			const Real cadet = cdx * ady - adx * cdy;
			const Real alift = adx * adx + ady * ady;
			const Real blift = bdx * bdx + bdy * bdy;
			const Real clift = cdx * cdx + cdy * cdy;
			const Real determinant = alift * bcdet + blift * cadet + clift * abdet;
			const Real permanent = alift * (std::abs(bdx * cdy) + std::abs(cdx * bdy))
				+ blift * (std::abs(cdx * ady) + std::abs(adx * cdy))
				+ clift * (std::abs(adx * bdy) + std::abs(bdx * ady));
			const Real epsilon = std::numeric_limits<Real>::epsilon();
			if (std::abs(determinant) > (REAL(10.0) + REAL(96.0) * epsilon) * epsilon * permanent)
				return Sign(determinant);

			const auto adxExact = Difference(a.X(), d.X());
			const auto adyExact = Difference(a.Y(), d.Y());
			const auto bdxExact = Difference(b.X(), d.X());
			const auto bdyExact = Difference(b.Y(), d.Y());
			const auto cdxExact = Difference(c.X(), d.X());
			const auto cdyExact = Difference(c.Y(), d.Y());
			const auto abdetExact = ExactCross(adxExact, adyExact, bdxExact, bdyExact);
			const auto bcdetExact = ExactCross(bdxExact, bdyExact, cdxExact, cdyExact);
			const auto cadetExact = ExactCross(cdxExact, cdyExact, adxExact, adyExact);
			const auto aliftExact = Add(Multiply(adxExact, adxExact), Multiply(adyExact, adyExact));
			const auto bliftExact = Add(Multiply(bdxExact, bdxExact), Multiply(bdyExact, bdyExact));
			const auto cliftExact = Add(Multiply(cdxExact, cdxExact), Multiply(cdyExact, cdyExact));

			return ExpansionSign(Add(Add(Multiply(aliftExact, bcdetExact),
				Multiply(bliftExact, cadetExact)), Multiply(cliftExact, abdetExact)));
		}

		/// Sign of the oriented volume: +1 when d is below oriented plane (a,b,c),
		/// -1 when above, and 0 when the four points are coplanar.
		static int Orientation3D(const Point3Cartesian& a, const Point3Cartesian& b,
			const Point3Cartesian& c, const Point3Cartesian& d) {
			const std::array<Point3Cartesian, 4> points = {a, b, c, d};
			Real maxCoordinate = REAL(0.0);
			for (const auto& point : points) {
				maxCoordinate = std::max(maxCoordinate, std::abs(point.X()));
				maxCoordinate = std::max(maxCoordinate, std::abs(point.Y()));
				maxCoordinate = std::max(maxCoordinate, std::abs(point.Z()));
			}
			int coordinateExponent = 0;
			if (std::isfinite(maxCoordinate) && maxCoordinate != REAL(0.0))
				std::frexp(maxCoordinate, &coordinateExponent);
			const int safeCoordinateExponent = (std::numeric_limits<Real>::max_exponent - 4) / 3;
			const int scaleExponent = std::min(0, safeCoordinateExponent - coordinateExponent);
			const auto scaled = [scaleExponent](Real value) {
				return std::scalbn(value, scaleExponent);
			};

			const Real ax = scaled(a.X()), ay = scaled(a.Y()), az = scaled(a.Z());
			const Real bx = scaled(b.X()), by = scaled(b.Y()), bz = scaled(b.Z());
			const Real cx = scaled(c.X()), cy = scaled(c.Y()), cz = scaled(c.Z());
			const Real dx = scaled(d.X()), dy = scaled(d.Y()), dz = scaled(d.Z());
			const Real adx = ax - dx, ady = ay - dy, adz = az - dz;
			const Real bdx = bx - dx, bdy = by - dy, bdz = bz - dz;
			const Real cdx = cx - dx, cdy = cy - dy, cdz = cz - dz;
			const Real bdxcdy = bdx * cdy, cdxbdy = cdx * bdy;
			const Real cdxady = cdx * ady, adxcdy = adx * cdy;
			const Real adxbdy = adx * bdy, bdxady = bdx * ady;
			const Real determinant = adz * (bdxcdy - cdxbdy)
				+ bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);
			const Real permanent = (std::abs(bdxcdy) + std::abs(cdxbdy)) * std::abs(adz)
				+ (std::abs(cdxady) + std::abs(adxcdy)) * std::abs(bdz)
				+ (std::abs(adxbdy) + std::abs(bdxady)) * std::abs(cdz);
			const Real epsilon = std::numeric_limits<Real>::epsilon();
			if (std::abs(determinant) > (REAL(7.0) + REAL(56.0) * epsilon) * epsilon * permanent)
				return Sign(determinant);

			const auto adxExact = Difference(ax, dx);
			const auto adyExact = Difference(ay, dy);
			const auto adzExact = Difference(az, dz);
			const auto bdxExact = Difference(bx, dx);
			const auto bdyExact = Difference(by, dy);
			const auto bdzExact = Difference(bz, dz);
			const auto cdxExact = Difference(cx, dx);
			const auto cdyExact = Difference(cy, dy);
			const auto cdzExact = Difference(cz, dz);
			const auto bc = ExactCross(bdxExact, bdyExact, cdxExact, cdyExact);
			const auto ca = ExactCross(cdxExact, cdyExact, adxExact, adyExact);
			const auto ab = ExactCross(adxExact, adyExact, bdxExact, bdyExact);
			return ExpansionSign(Add(Add(Multiply(adzExact, bc), Multiply(bdzExact, ca)),
				Multiply(cdzExact, ab)));
		}

		/// Sign of the in-sphere determinant for positively oriented tetrahedron
		/// (a,b,c,d): +1 inside, -1 outside, and 0 when e is cospherical.
		static int InSphere3D(const Point3Cartesian& a, const Point3Cartesian& b,
			const Point3Cartesian& c, const Point3Cartesian& d, const Point3Cartesian& e) {
			const std::array<Point3Cartesian, 5> inputPoints = {a, b, c, d, e};
			Real maxCoordinate = REAL(0.0);
			for (const auto& point : inputPoints) {
				maxCoordinate = std::max(maxCoordinate, std::abs(point.X()));
				maxCoordinate = std::max(maxCoordinate, std::abs(point.Y()));
				maxCoordinate = std::max(maxCoordinate, std::abs(point.Z()));
			}
			int coordinateExponent = 0;
			if (std::isfinite(maxCoordinate) && maxCoordinate != REAL(0.0))
				std::frexp(maxCoordinate, &coordinateExponent);
			const int safeCoordinateExponent = (std::numeric_limits<Real>::max_exponent - 6) / 5;
			const int scaleExponent = std::min(0, safeCoordinateExponent - coordinateExponent);
			const auto scaled = [scaleExponent](Real value) {
				return std::scalbn(value, scaleExponent);
			};
			const std::array<Point3Cartesian, 5> points = {
				Point3Cartesian(scaled(a.X()), scaled(a.Y()), scaled(a.Z())),
				Point3Cartesian(scaled(b.X()), scaled(b.Y()), scaled(b.Z())),
				Point3Cartesian(scaled(c.X()), scaled(c.Y()), scaled(c.Z())),
				Point3Cartesian(scaled(d.X()), scaled(d.Y()), scaled(d.Z())),
				Point3Cartesian(scaled(e.X()), scaled(e.Y()), scaled(e.Z()))
			};
			const Real aex = points[0].X() - points[4].X();
			const Real aey = points[0].Y() - points[4].Y();
			const Real aez = points[0].Z() - points[4].Z();
			const Real bex = points[1].X() - points[4].X();
			const Real bey = points[1].Y() - points[4].Y();
			const Real bez = points[1].Z() - points[4].Z();
			const Real cex = points[2].X() - points[4].X();
			const Real cey = points[2].Y() - points[4].Y();
			const Real cez = points[2].Z() - points[4].Z();
			const Real dex = points[3].X() - points[4].X();
			const Real dey = points[3].Y() - points[4].Y();
			const Real dez = points[3].Z() - points[4].Z();
			const Real ab = aex * bey - bex * aey;
			const Real bc = bex * cey - cex * bey;
			const Real cd = cex * dey - dex * cey;
			const Real da = dex * aey - aex * dey;
			const Real ac = aex * cey - cex * aey;
			const Real bd = bex * dey - dex * bey;
			const Real abc = aez * bc - bez * ac + cez * ab;
			const Real bcd = bez * cd - cez * bd + dez * bc;
			const Real cda = cez * da + dez * ac + aez * cd;
			const Real dab = dez * ab + aez * bd + bez * da;
			const Real alift = aex * aex + aey * aey + aez * aez;
			const Real blift = bex * bex + bey * bey + bez * bez;
			const Real clift = cex * cex + cey * cey + cez * cez;
			const Real dlift = dex * dex + dey * dey + dez * dez;
			const Real determinant = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);
			const Real permanent = ((std::abs(cex * dey) + std::abs(dex * cey)) * std::abs(bez)
				+ (std::abs(dex * bey) + std::abs(bex * dey)) * std::abs(cez)
				+ (std::abs(bex * cey) + std::abs(cex * bey)) * std::abs(dez)) * alift
				+ ((std::abs(dex * aey) + std::abs(aex * dey)) * std::abs(cez)
				+ (std::abs(aex * cey) + std::abs(cex * aey)) * std::abs(dez)
				+ (std::abs(cex * dey) + std::abs(dex * cey)) * std::abs(aez)) * blift
				+ ((std::abs(aex * bey) + std::abs(bex * aey)) * std::abs(dez)
				+ (std::abs(bex * dey) + std::abs(dex * bey)) * std::abs(aez)
				+ (std::abs(dex * aey) + std::abs(aex * dey)) * std::abs(bez)) * clift
				+ ((std::abs(bex * cey) + std::abs(cex * bey)) * std::abs(aez)
				+ (std::abs(cex * aey) + std::abs(aex * cey)) * std::abs(bez)
				+ (std::abs(aex * bey) + std::abs(bex * aey)) * std::abs(cez)) * dlift;
			const Real epsilon = std::numeric_limits<Real>::epsilon();
			if (std::abs(determinant) > (REAL(16.0) + REAL(224.0) * epsilon) * epsilon * permanent)
				return Sign(determinant);

			std::array<Expansion<6>, 5> lifts;
			for (std::size_t i = 0; i < points.size(); ++i) {
				const auto x = Scalar(points[i].X());
				const auto y = Scalar(points[i].Y());
				const auto z = Scalar(points[i].Z());
				lifts[i] = Add(Add(Multiply(x, x), Multiply(y, y)), Multiply(z, z));
			}

			Expansion<5760> exactDeterminant;
			exactDeterminant.components[exactDeterminant.size++] = REAL(0.0);
			Expansion<5760> scratch;
			for (int px = 0; px < 5; ++px)
				for (int py = 0; py < 5; ++py) if (py != px)
					for (int pz = 0; pz < 5; ++pz) if (pz != px && pz != py)
						for (int pl = 0; pl < 5; ++pl) if (pl != px && pl != py && pl != pz)
							for (int po = 0; po < 5; ++po) if (po != px && po != py && po != pz && po != pl) {
								const std::array<int, 5> permutation = {px, py, pz, pl, po};
								int inversions = 0;
								for (int i = 0; i < 5; ++i)
									for (int j = i + 1; j < 5; ++j)
										if (permutation[i] > permutation[j]) ++inversions;
								const Real termSign = (inversions % 2 == 0) ? REAL(1.0) : -REAL(1.0);
								const auto termX = Scale(lifts[pl], points[px].X());
								const auto termXY = Scale(termX, points[py].Y());
								const auto termXYZ = Scale(termXY, points[pz].Z());
								for (std::size_t i = 0; i < termXYZ.size; ++i)
									Append(exactDeterminant, termSign * termXYZ.components[i], scratch);
							}
			return ExpansionSign(exactDeterminant);
		}
	};
}

#endif