///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry3DPlane.h                                                   ///
///  Description: Plane3D class for 3D plane representation and operations            ///
///               Point/line/plane tests, distances, and intersections               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GEOMETRY_3D_PLANE_H
#define MML_GEOMETRY_3D_PLANE_H

#include <mml/MMLBase.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Vector/VectorTypes3D.h>
#include <mml/base/Geometry/Geometry.h>
#include <mml/base/Geometry/GeometryBase/GeometryArithmetic.h>
#include <mml/base/BaseUtils.h>
#include <mml/base/Geometry/Geometry3D/Geometry3DLines.h>

namespace MML
{

	/// @brief Plane in 3D space
	/// @details Represented in general form Ax + By + Cz + D = 0 with a finite nonzero normal.
	///          Provides multiple construction methods, point/line/plane operations,
	///          distance calculations, and intersection computations.
	class Plane3D
	{
	private:
		Real _coefA, _coefB, _coefC, _coefD;
		using Binary = detail::GeometryBinary;
		using BinaryVector = std::array<Binary, 3>;

		BinaryVector CheckedNormal() const {
			const BinaryVector normal{Binary(_coefA), Binary(_coefB), Binary(_coefC)};
			detail::ValidateGeometryCoordinates({_coefD}, "Plane3D query");
			if (normal[0].fraction == 0 && normal[1].fraction == 0 && normal[2].fraction == 0)
				throw GeometryError("Plane3D query - zero normal");
			return normal;
		}
		static BinaryVector Coordinates(const Pnt3Cart& p) { return {Binary(p.X()), Binary(p.Y()), Binary(p.Z())}; }
		static BinaryVector Direction(const Vec3Cart& v) {
			(void)detail::GeometryUnitDirection<3>(v);
			return {Binary(v.X()), Binary(v.Y()), Binary(v.Z())};
		}
		static Binary Dot(const BinaryVector& a, const BinaryVector& b, Binary* remainder = nullptr) {
			return Binary::SumProducts<3>(a, b, remainder);
		}
		static Binary Negated(Binary value) { value.fraction = -value.fraction; return value; }
		static std::array<Binary, 2> CorrectedQuotient(const Binary& numerator, const Binary& denominator,
			const Binary& numeratorLow, const Binary& denominatorLow) {
			const auto quotient = Binary::Quotient(numerator, denominator);
			const auto error = Binary::SumProducts<4>({numerator, numeratorLow, Negated(denominator), Negated(denominatorLow)},
				{Binary(1), Binary(1), quotient, quotient});
			return {quotient, Binary::Quotient(error, denominator)};
		}
		static BinaryVector Cross(const BinaryVector& a, const BinaryVector& b) {
			BinaryVector result;
			for (int i = 0; i < 3; ++i) {
				const int j = (i + 1) % 3, k = (i + 2) % 3;
				result[i] = Binary::DifferenceOfProducts(a[j], b[k], a[k], b[j]);
			}
			return result;
		}
		static void CheckTolerance(Real tolerance) {
			if (!std::isfinite(tolerance) || tolerance < 0)
				throw GeometryError("Plane3D query - tolerance must be finite and nonnegative");
		}
		Vec3Cart UnitNormal() const {
			(void)CheckedNormal();
			return detail::GeometryUnitDirection<3>(Normal());
		}
		Binary Evaluation(const Pnt3Cart& point, Binary* remainder = nullptr) const {
			return detail::GeometryPlaneEvaluation(point, Normal(), _coefD, remainder);
		}

	public:
		/// @brief Constructs a plane from a point and normal vector
		/// @param a A point on the plane
		/// @param normal The plane normal (will be normalized)
		/// @throws GeometryError for a zero/nonfinite normal, nonfinite point or unrepresentable offset
		Plane3D(const Pnt3Cart& a, const Vec3Cart& normal)
		{
			const Vec3Cart unitNormal = detail::GeometryUnitDirection<3>(normal);

			_coefA = unitNormal.X();
			_coefB = unitNormal.Y();
			_coefC = unitNormal.Z();
			_coefD = detail::GeometryPlaneOffset(a, unitNormal);
		}
		
		/// @brief Constructs a plane from three points
		/// @param a First point
		/// @param b Second point
		/// @param c Third point
		Plane3D(const Pnt3Cart& a, const Pnt3Cart& b, const Pnt3Cart& c)
			: Plane3D(a, detail::GeometryCrossProduct(a, b, a, c).Unit())
		{
		}
		
		/// @brief Constructs a plane in Hesse normal form
		/// @param alpha Direction cosine angle for x
		/// @param beta Direction cosine angle for y
		/// @param gamma Direction cosine angle for z
		/// @param d Distance from origin
		Plane3D(Real alpha, Real beta, Real gamma, Real d)
		{
			_coefA = cos(alpha);
			_coefB = cos(beta);
			_coefC = cos(gamma);
			_coefD = -d;
		}
		
		/// @brief Constructs a plane from intercepts on coordinate axes
		/// @param seg_x X-axis intercept
		/// @param seg_y Y-axis intercept
		/// @param seg_z Z-axis intercept
		/// @throws GeometryError if any intercept is zero
		Plane3D(Real seg_x, Real seg_y, Real seg_z)
		{
			if (seg_x == 0 || seg_y == 0 || seg_z == 0)
				throw GeometryError("Plane3D ctor - zero segment");

			_coefA = 1 / seg_x;
			_coefB = 1 / seg_y;
			_coefC = 1 / seg_z;
			_coefD = -1;
		}

		/// @brief Gets the XY plane (z = 0)
		static Plane3D XYPlane() { return Plane3D(Pnt3Cart(0, 0, 0), Vec3Cart(0, 0, 1)); }
		/// @deprecated Use XYPlane() (CONTRACT_POLICY §11: factories are PascalCase nouns).
		[[deprecated("use XYPlane()")]]
		static Plane3D GetXYPlane()  { return XYPlane(); }
		/// @brief Gets the XZ plane (y = 0)
		static Plane3D XZPlane() { return Plane3D(Pnt3Cart(0, 0, 0), Vec3Cart(0, 1, 0)); }
		/// @deprecated Use XZPlane() (CONTRACT_POLICY §11: factories are PascalCase nouns).
		[[deprecated("use XZPlane()")]]
		static Plane3D GetXZPlane()  { return XZPlane(); }
		/// @brief Gets the YZ plane (x = 0)
		static Plane3D YZPlane() { return Plane3D(Pnt3Cart(0, 0, 0), Vec3Cart(1, 0, 0)); }
		/// @deprecated Use YZPlane() (CONTRACT_POLICY §11: factories are PascalCase nouns).
		[[deprecated("use YZPlane()")]]
		static Plane3D GetYZPlane()  { return YZPlane(); }

		/// @brief Gets coefficient A (const)
		Real  A() const { return _coefA; }
		/// @brief Gets coefficient A (mutable)
		Real& A() { return _coefA; }
		/// @brief Gets coefficient B (const)
		Real  B() const { return _coefB; }
		/// @brief Gets coefficient B (mutable)
		Real& B() { return _coefB; }
		/// @brief Gets coefficient C (const)
		Real  C() const { return _coefC; }
		/// @brief Gets coefficient C (mutable)
		Real& C() { return _coefC; }
		/// @brief Gets coefficient D (const)
		Real  D() const { return _coefD; }
		/// @brief Gets coefficient D (mutable)
		Real& D() { return _coefD; }

		/// @brief Gets the stored normal coefficients (not necessarily unit length)
		Vec3Cart	Normal() const { return Vec3Cart(_coefA, _coefB, _coefC); }
		
		/// @brief Gets the point on the plane closest to the origin
		Pnt3Cart	getPointOnPlane() const {
			return ProjectionToPlane(Pnt3Cart(0, 0, 0));
		}
		/// @deprecated Use getPointOnPlane() (CONTRACT_POLICY §11: accessors are getX()).
		[[deprecated("use getPointOnPlane()")]]
		Pnt3Cart	GetPointOnPlane() const { return getPointOnPlane(); }

		/// @brief Gets the intercepts on coordinate axes
		/// @param outseg_x Output: X-axis intercept (PosInf if parallel)
		/// @param outseg_y Output: Y-axis intercept (PosInf if parallel)
		/// @param outseg_z Output: Z-axis intercept (PosInf if parallel)
		void getCoordAxisSegments(Real& outseg_x, Real& outseg_y, Real& outseg_z) const
		{
			const auto normal = CheckedNormal();
			std::array<Real, 3> segments;
			for (int i = 0; i < 3; ++i)
				segments[i] = normal[i].fraction == 0 ? Constants::PosInf
					: Binary::Quotient(Binary(-_coefD), normal[i]).CheckedValue();
			outseg_x = segments[0]; outseg_y = segments[1]; outseg_z = segments[2];
		}
		/// @deprecated Use getCoordAxisSegments().
		[[deprecated("use getCoordAxisSegments()")]]
		void GetCoordAxisSegments(Real& outseg_x, Real& outseg_y, Real& outseg_z) const { getCoordAxisSegments(outseg_x, outseg_y, outseg_z); }

		/// @brief Tests if a point lies on this plane
		/// @param pnt The point to test
		/// @param defEps Tolerance for the test
		bool IsPointOnPlane(const Pnt3Cart& pnt, Real defEps = Precision::Plane3DIsPointOnPlaneTolerance) const
		{
			(void)CheckedNormal();
			CheckTolerance(defEps);
			return Evaluation(pnt).AbsLessThan(defEps);
		}
		
		/// @brief Calculates the distance from a point to this plane
		Real DistToPoint(const Pnt3Cart& pnt) const
		{
			const auto normal = CheckedNormal();
			return std::abs(Binary::Quotient(Evaluation(pnt), Dot(normal, normal).Sqrt()).CheckedValue());
		}

		/// @brief Projects a point onto this plane
		/// @param pnt The point to project
		/// @return The projected point on the plane
		/// @throws GeometryError for invalid geometry or an unrepresentable output coordinate
		Pnt3Cart ProjectionToPlane(const Pnt3Cart& pnt) const
		{
			const auto normal = CheckedNormal();
			const auto point = Coordinates(pnt);
			Binary squareLow;
			const auto square = Dot(normal, normal, &squareLow);
			std::array<Real, 3> result;
			for (int i = 0; i < 3; ++i) {
				// Expand p_i*|n|^2 - n_i*(n.p+D) before dividing. Exact product
				// cancellation preserves true zero coordinates even at subnormal scales.
				std::array<Binary, 13> left{}, right{};
				for (int j = 0; j < 3; ++j) {
					const auto first = Binary::ProductParts(normal[j], normal[j]);
					const auto second = Binary::ProductParts(Negated(normal[i]), normal[j]);
					for (int part = 0; part < 2; ++part) {
						left[4 * j + part] = first[part]; right[4 * j + part] = point[i];
						left[4 * j + 2 + part] = second[part]; right[4 * j + 2 + part] = point[j];
					}
				}
				left[12] = Negated(normal[i]); right[12] = Binary(_coefD);
				Binary low;
				const auto numerator = Binary::SumProducts(left, right, &low);
				const auto quotient = CorrectedQuotient(numerator, square, low, squareLow);
				result[i] = Binary::Add(quotient[0], quotient[1]).CheckedValue();
			}
			return {result[0], result[1], result[2]};
		}

		/// @brief Tests if a line lies entirely on this plane
		bool IsLineOnPlane(const Line3D& line) const
		{
			const auto normal = CheckedNormal();
			const auto direction = Direction(line.Direction());
			const auto start = Evaluation(line.StartPoint());
			return start.AbsLessThan(Precision::Plane3DIsPointOnPlaneTolerance) &&
				Binary::Add(start, Dot(normal, direction)).AbsLessThan(Precision::Plane3DIsPointOnPlaneTolerance);
		}
		
		/// @brief Calculates the angle between a line and this plane
		/// @return Angle in radians
		Real AngleToLine(const Line3D& line) const
		{
			// angle between line and normal to plane
			auto direction = Vec3Cart(detail::GeometryUnitDirection<3>(line.Direction()));
			return Constants::PI / MML::Real(2.0) - direction.AngleToVector(UnitNormal());
		}
		
		/// @brief Finds the intersection point of a line with this plane
		/// @param line The line
		/// @param out_inter_pnt Output: intersection point if exists
		/// @return false if line is parallel to plane
		bool IntersectionWithLine(const Line3D& line, Pnt3Cart& out_inter_pnt, Real eps = Precision::Line3DIsParallelTolerance) const
		{
			CheckTolerance(eps);
			const auto normal = CheckedNormal(), direction = Direction(line.Direction());
			const auto point = Coordinates(line.StartPoint());
			Binary dotLow;
			const auto dot = Dot(normal, direction, &dotLow);
			const auto lengths = Binary::Product(Dot(normal, normal).Sqrt(), Dot(direction, direction).Sqrt());
			if (dot.fraction == 0 || Binary::Quotient(dot, lengths).AbsLessThan(eps)) return false;
			Binary residualLow;
			const auto residual = Evaluation(line.StartPoint(), &residualLow);
			const auto parameter = CorrectedQuotient(residual, dot, residualLow, dotLow);
			std::array<Real, 3> result;
			for (int i = 0; i < 3; ++i)
				result[i] = Binary::SumProducts<3>({point[i], Negated(direction[i]), Negated(direction[i])},
					{Binary(1), parameter[0], parameter[1]}).CheckedValue();
			out_inter_pnt = Pnt3Cart(result[0], result[1], result[2]);
			return true;
		}

		/// @brief Tests if this plane is parallel to another plane
		bool IsParallelToPlane(const Plane3D& plane) const
		{
			return UnitNormal().IsParallelTo(plane.UnitNormal());
		}
		
		/// @brief Tests if this plane is perpendicular to another plane
		bool IsPerpendicularToPlane(const Plane3D& plane) const
		{
			return UnitNormal().IsPerpendicularTo(plane.UnitNormal());
		}
		
		/// @brief Calculates the angle between two planes
		/// @return Angle in radians
		Real AngleToPlane(const Plane3D& plane) const
		{
			// angle between normals of two planes
			auto normal = UnitNormal();
			return normal.AngleToVector(plane.UnitNormal());
		}
		
		/// @brief Calculates the distance between two planes
		/// @return Distance (0 if planes intersect)
		Real DistToPlane(const Plane3D& plane) const
		{
			if (IsParallelToPlane(plane))
			{
				const auto n1 = CheckedNormal(), n2 = plane.CheckedNormal();
				const auto d1 = Binary::Quotient(Binary(_coefD), Dot(n1, n1).Sqrt());
				auto d2 = Binary::Quotient(Binary(plane._coefD), Dot(n2, n2).Sqrt());
				if (Dot(n1, n2).fraction > 0) d2 = Negated(d2);
				return std::abs(Binary::Add(d1, d2).CheckedValue());
			}
			else
				return 0.0;
		}

		/// @brief Finds the intersection line of two planes
		/// @param plane The other plane
		/// @param out_inter_line Output: intersection line if exists
		/// @return false if planes are parallel
		bool IntersectionWithPlane(const Plane3D& plane, Line3D& out_inter_line) const
		{
			const auto n1 = CheckedNormal(), n2 = plane.CheckedNormal();
			const auto cross = Cross(n1, n2);
			const auto square = Dot(cross, cross);
			const auto lengths = Binary::Product(Dot(n1, n1).Sqrt(), Dot(n2, n2).Sqrt());
			if (square.fraction == 0 ||
				Binary::Quotient(square.Sqrt(), lengths).AbsLessThan(Precision::Vec3CartIsParallelTolerance))
				return false;
			BinaryVector weighted;
			for (int i = 0; i < 3; ++i)
				weighted[i] = Binary::DifferenceOfProducts(Binary(plane._coefD), n1[i], Binary(_coefD), n2[i]);
			const auto numerator = Cross(weighted, cross);
			std::array<Real, 3> point;
			int exponent = cross[0].fraction != 0 ? cross[0].exponent
				: cross[1].fraction != 0 ? cross[1].exponent : cross[2].exponent;
			for (int i = 0; i < 3; ++i) {
				point[i] = Binary::Quotient(numerator[i], square).CheckedValue();
				if (cross[i].fraction != 0) exponent = std::max(exponent, cross[i].exponent);
			}
			Vec3Cart direction;
			for (int i = 0; i < 3; ++i)
				direction[i] = cross[i].fraction == 0 ? 0 : std::ldexp(cross[i].fraction, cross[i].exponent - exponent);
			out_inter_line = Line3D(Pnt3Cart(point[0], point[1], point[2]), direction);
			return true;
		}
	};

} // namespace MML

#endif // MML_GEOMETRY_3D_PLANE_H
