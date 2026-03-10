///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry2DShapes.h                                                  ///
///  Description: Pure 2D geometric shapes - coordinate-free intrinsic geometry       ///
///               (Triangle, Circle, Ellipse, Sector, Segment, Annulus, etc.)         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GEOMETRY_2D_SHAPES_H
#define MML_GEOMETRY_2D_SHAPES_H

#include <algorithm>
#include <cmath>

#include "mml/MMLBase.h"


namespace MML {

	//=============================================================================
	// PURE 2D GEOMETRIC SHAPES - Coordinate-free intrinsic geometry
	//=============================================================================

	/// @brief Triangle defined by three side lengths (coordinate-free)
	/// Represents a triangle purely by its side lengths a, b, c.
	/// All properties are intrinsic - independent of position or orientation.

	class Triangle {
	private:
		Real _a, _b, _c;

	public:
		// Accessors
		Real A() const { return _a; }
		Real& A() { return _a; }
		Real B() const { return _b; }
		Real& B() { return _b; }
		Real C() const { return _c; }
		Real& C() { return _c; }

		// Constructors
		Triangle()
			: _a(0.0)
			, _b(0.0)
			, _c(0.0) {}
		Triangle(Real a, Real b, Real c)
			: _a(a)
			, _b(b)
			, _c(c) {}

		// Validity
		bool IsValid() const { return _a > 0 && _b > 0 && _c > 0 && _a + _b > _c && _a + _c > _b && _b + _c > _a; }

		// Basic properties
		Real Perimeter() const { return _a + _b + _c; }
		Real Semiperimeter() const { return (_a + _b + _c) / 2.0; }

		Real Area() const {
			Real s = Semiperimeter();
			return sqrt(s * (s - _a) * (s - _b) * (s - _c)); // Heron's formula
		}

		// Angles (in radians) - using law of cosines
		Real AngleA() const { return acos((_b * _b + _c * _c - _a * _a) / (2 * _b * _c)); } // opposite to side a
		Real AngleB() const { return acos((_a * _a + _c * _c - _b * _b) / (2 * _a * _c)); } // opposite to side b
		Real AngleC() const { return acos((_a * _a + _b * _b - _c * _c) / (2 * _a * _b)); } // opposite to side c

		// Altitudes (heights)
		Real AltitudeToA() const { return 2.0 * Area() / _a; }
		Real AltitudeToB() const { return 2.0 * Area() / _b; }
		Real AltitudeToC() const { return 2.0 * Area() / _c; }

		// Medians
		Real MedianToA() const { return 0.5 * sqrt(2 * _b * _b + 2 * _c * _c - _a * _a); }
		Real MedianToB() const { return 0.5 * sqrt(2 * _a * _a + 2 * _c * _c - _b * _b); }
		Real MedianToC() const { return 0.5 * sqrt(2 * _a * _a + 2 * _b * _b - _c * _c); }

		// Inscribed and circumscribed circles
		Real Inradius() const { return Area() / Semiperimeter(); }
		Real Circumradius() const { return (_a * _b * _c) / (4.0 * Area()); }

		// Classification
		bool IsRight(Real eps = Defaults::ShapePropertyTolerance) const {
			Real a2 = _a * _a, b2 = _b * _b, c2 = _c * _c;
			return std::abs(a2 + b2 - c2) < eps || std::abs(a2 + c2 - b2) < eps || std::abs(b2 + c2 - a2) < eps;
		}
		bool IsIsosceles(Real eps = Defaults::ShapePropertyTolerance) const { return std::abs(_a - _b) < eps || std::abs(_a - _c) < eps || std::abs(_b - _c) < eps; }
		bool IsEquilateral(Real eps = Defaults::ShapePropertyTolerance) const { return std::abs(_a - _b) < eps && std::abs(_a - _c) < eps; }
		bool IsAcute() const {
			Real a2 = _a * _a, b2 = _b * _b, c2 = _c * _c;
			return (a2 + b2 > c2) && (a2 + c2 > b2) && (b2 + c2 > a2);
		}
		bool IsObtuse() const {
			Real a2 = _a * _a, b2 = _b * _b, c2 = _c * _c;
			return (a2 + b2 < c2) || (a2 + c2 < b2) || (b2 + c2 < a2);
		}
		bool IsScalene(Real eps = Defaults::ShapePropertyTolerance) const { return std::abs(_a - _b) > eps && std::abs(_a - _c) > eps && std::abs(_b - _c) > eps; }
	};

	/// @brief Circle defined by radius (coordinate-free)

	class Circle {
	private:
		Real _radius;

	public:
		Real Radius() const { return _radius; }
		Real& Radius() { return _radius; }
		Real Diameter() const { return 2.0 * _radius; }

		Circle()
			: _radius(0.0) {}
		explicit Circle(Real radius)
			: _radius(radius) {}

		bool IsValid() const { return _radius > 0; }

		Real Area() const { return Constants::PI * _radius * _radius; }
		Real Circumference() const { return 2.0 * Constants::PI * _radius; }
		Real Perimeter() const { return Circumference(); }

		// Arc and sector for given central angle (radians)
		Real ArcLength(Real angle) const { return _radius * angle; }
		Real SectorArea(Real angle) const { return 0.5 * _radius * _radius * angle; }

		// Chord length for given central angle
		Real ChordLength(Real angle) const { return 2.0 * _radius * sin(angle / 2.0); }

		// Segment area (region between chord and arc)
		Real SegmentArea(Real angle) const { return 0.5 * _radius * _radius * (angle - sin(angle)); }

		// Relationships
		static Circle FromArea(Real area) { return Circle(sqrt(area / Constants::PI)); }
		static Circle FromCircumference(Real circ) { return Circle(circ / (2.0 * Constants::PI)); }
		static Circle Inscribed(const Triangle& t) { return Circle(t.Inradius()); }
		static Circle Circumscribed(const Triangle& t) { return Circle(t.Circumradius()); }
	};

	/// @brief Ellipse defined by semi-major and semi-minor axes (coordinate-free)

	class Ellipse {
	private:
		Real _a; // semi-major axis
		Real _b; // semi-minor axis

	public:
		Real SemiMajor() const { return _a; }
		Real& SemiMajor() { return _a; }
		Real SemiMinor() const { return _b; }
		Real& SemiMinor() { return _b; }
		Real MajorAxis() const { return 2.0 * _a; }
		Real MinorAxis() const { return 2.0 * _b; }

		Ellipse()
			: _a(0.0)
			, _b(0.0) {}
		Ellipse(Real semiMajor, Real semiMinor)
			: _a(semiMajor)
			, _b(semiMinor) {}

		bool IsValid() const { return _a > 0 && _b > 0; }
		bool IsCircle(Real eps = Defaults::ShapePropertyTolerance) const { return std::abs(_a - _b) < eps; }

		Real Area() const { return Constants::PI * _a * _b; }

		// Circumference approximation (Ramanujan's formula)
		Real Circumference() const {
			Real h = POW2(_a - _b) / POW2(_a + _b);
			return Constants::PI * (_a + _b) * (1.0 + 3.0 * h / (10.0 + sqrt(4.0 - 3.0 * h)));
		}
		Real Perimeter() const { return Circumference(); }

		// Eccentricity e = sqrt(1 - b²/a²), ranges from 0 (circle) to 1 (parabola limit)
		Real Eccentricity() const {
			if (_a == 0)
				return 0;
			return sqrt(1.0 - (_b * _b) / (_a * _a));
		}

		// Linear eccentricity (distance from center to focus)
		Real LinearEccentricity() const { return sqrt(_a * _a - _b * _b); }

		// Focal distance (distance between foci)
		Real FocalDistance() const { return 2.0 * LinearEccentricity(); }

		// Semi-latus rectum (half the chord through focus perpendicular to major axis)
		Real SemiLatusRectum() const { return (_b * _b) / _a; }

		// Directrix distance from center
		Real DirectrixDistance() const { return _a / Eccentricity(); }

		static Ellipse FromFoci(Real focalDist, Real majorAxis) {
			Real a = majorAxis / 2.0;
			Real c = focalDist / 2.0;
			Real b = sqrt(a * a - c * c);
			return Ellipse(a, b);
		}
	};

	/// @brief Circular sector (pie slice) defined by radius and central angle
	/// @details A circular sector is the portion of a disk enclosed by two radii
	/// and an arc. Sometimes called a "pie slice" or "wedge".

	class CircularSector {
	private:
		Real _radius;
		Real _angle; // central angle in radians

	public:
		/// @brief Gets the radius (const version)
		Real Radius() const { return _radius; }
		/// @brief Gets the radius (mutable reference)
		Real& Radius() { return _radius; }
		/// @brief Gets the central angle in radians (const version)
		Real Angle() const { return _angle; }
		/// @brief Gets the central angle in radians (mutable reference)
		Real& Angle() { return _angle; }
		/// @brief Gets the central angle in degrees
		Real AngleDegrees() const { return _angle * 180.0 / Constants::PI; }

		/// @brief Default constructor - creates invalid sector
		CircularSector()
			: _radius(0.0)
			, _angle(0.0) {}
		/// @brief Constructs a sector from radius and central angle
		/// @param radius The radius of the sector
		/// @param angle The central angle in radians
		CircularSector(Real radius, Real angle)
			: _radius(radius)
			, _angle(angle) {}

		/// @brief Checks if the sector parameters are valid
		/// @return True if radius > 0 and 0 < angle ≤ 2π
		bool IsValid() const { return _radius > 0 && _angle > 0 && _angle <= 2.0 * Constants::PI; }

		/// @brief Calculates the area of the sector
		/// @return Area = ½r²θ
		Real Area() const { return 0.5 * _radius * _radius * _angle; }
		/// @brief Calculates the arc length
		/// @return Arc length = rθ
		Real ArcLength() const { return _radius * _angle; }
		/// @brief Calculates the perimeter (two radii plus arc)
		/// @return Perimeter = 2r + arc length
		Real Perimeter() const { return 2.0 * _radius + ArcLength(); }
		/// @brief Calculates the chord length (straight line between arc endpoints)
		/// @return Chord length = 2r·sin(θ/2)
		Real ChordLength() const { return 2.0 * _radius * sin(_angle / 2.0); }
	};

	/// @brief Circular segment (region between chord and arc)

	class CircularSegment {
	private:
		Real _radius;
		Real _angle; // central angle in radians

	public:
		Real Radius() const { return _radius; }
		Real& Radius() { return _radius; }
		Real Angle() const { return _angle; }
		Real& Angle() { return _angle; }

		CircularSegment()
			: _radius(0.0)
			, _angle(0.0) {}
		CircularSegment(Real radius, Real angle)
			: _radius(radius)
			, _angle(angle) {}

		// Construct from radius and chord height (sagitta)
		static CircularSegment FromSagitta(Real radius, Real sagitta) {
			Real angle = 2.0 * acos(1.0 - sagitta / radius);
			return CircularSegment(radius, angle);
		}

		bool IsValid() const { return _radius > 0 && _angle > 0 && _angle <= 2.0 * Constants::PI; }

		Real Area() const { return 0.5 * _radius * _radius * (_angle - sin(_angle)); }
		Real ArcLength() const { return _radius * _angle; }
		Real ChordLength() const { return 2.0 * _radius * sin(_angle / 2.0); }
		Real Sagitta() const { return _radius * (1.0 - cos(_angle / 2.0)); } // height of segment
		Real Perimeter() const { return ArcLength() + ChordLength(); }
	};

	/// @brief Annulus (ring) defined by inner and outer radii

	class Annulus {
	private:
		Real _innerRadius;
		Real _outerRadius;

	public:
		Real InnerRadius() const { return _innerRadius; }
		Real& InnerRadius() { return _innerRadius; }
		Real OuterRadius() const { return _outerRadius; }
		Real& OuterRadius() { return _outerRadius; }
		Real Width() const { return _outerRadius - _innerRadius; }

		Annulus()
			: _innerRadius(0.0)
			, _outerRadius(0.0) {}
		Annulus(Real inner, Real outer)
			: _innerRadius(inner)
			, _outerRadius(outer) {}

		bool IsValid() const { return _innerRadius >= 0 && _outerRadius > _innerRadius; }

		Real Area() const { return Constants::PI * (_outerRadius * _outerRadius - _innerRadius * _innerRadius); }

		Real InnerCircumference() const { return 2.0 * Constants::PI * _innerRadius; }
		Real OuterCircumference() const { return 2.0 * Constants::PI * _outerRadius; }
		Real MeanCircumference() const { return Constants::PI * (_innerRadius + _outerRadius); }

		Circle InnerCircle() const { return Circle(_innerRadius); }
		Circle OuterCircle() const { return Circle(_outerRadius); }
	};

	/// @brief Regular polygon defined by number of sides and side length

	class RegularPolygon {
	private:
		int _n; // number of sides
		Real _sideLength;

	public:
		int NumSides() const { return _n; }
		int& NumSides() { return _n; }
		Real SideLength() const { return _sideLength; }
		Real& SideLength() { return _sideLength; }

		RegularPolygon()
			: _n(3)
			, _sideLength(0.0) {}
		RegularPolygon(int n, Real sideLength)
			: _n(n)
			, _sideLength(sideLength) {}

		bool IsValid() const { return _n >= 3 && _sideLength > 0; }

		Real InteriorAngle() const { return (_n - 2) * Constants::PI / _n; }
		Real ExteriorAngle() const { return 2.0 * Constants::PI / _n; }
		Real CentralAngle() const { return 2.0 * Constants::PI / _n; }

		Real Perimeter() const { return _n * _sideLength; }

		Real Area() const { return 0.25 * _n * _sideLength * _sideLength / tan(Constants::PI / _n); }

		// Apothem (distance from center to midpoint of side)
		Real Apothem() const { return _sideLength / (2.0 * tan(Constants::PI / _n)); }

		// Circumradius (distance from center to vertex)
		Real Circumradius() const { return _sideLength / (2.0 * sin(Constants::PI / _n)); }

		// Inradius = apothem
		Real Inradius() const { return Apothem(); }

		Circle InscribedCircle() const { return Circle(Inradius()); }
		Circle CircumscribedCircle() const { return Circle(Circumradius()); }

		// Factory methods
		static RegularPolygon FromCircumradius(int n, Real R) {
			Real side = 2.0 * R * sin(Constants::PI / n);
			return RegularPolygon(n, side);
		}

		static RegularPolygon FromInradius(int n, Real r) {
			Real side = 2.0 * r * tan(Constants::PI / n);
			return RegularPolygon(n, side);
		}

		static RegularPolygon FromArea(int n, Real area) {
			Real side = sqrt(4.0 * area * tan(Constants::PI / n) / n);
			return RegularPolygon(n, side);
		}

		// Named polygons
		static RegularPolygon EquilateralTriangle(Real side) { return RegularPolygon(3, side); }
		static RegularPolygon Square(Real side) { return RegularPolygon(4, side); }
		static RegularPolygon Pentagon(Real side) { return RegularPolygon(5, side); }
		static RegularPolygon Hexagon(Real side) { return RegularPolygon(6, side); }
		static RegularPolygon Heptagon(Real side) { return RegularPolygon(7, side); }
		static RegularPolygon Octagon(Real side) { return RegularPolygon(8, side); }
		static RegularPolygon Nonagon(Real side) { return RegularPolygon(9, side); }
		static RegularPolygon Decagon(Real side) { return RegularPolygon(10, side); }
		static RegularPolygon Dodecagon(Real side) { return RegularPolygon(12, side); }
	};

	/// @brief Rectangle defined by width and height (coordinate-free)

	class Rectangle {
	private:
		Real _width;
		Real _height;

	public:
		Real Width() const { return _width; }
		Real& Width() { return _width; }
		Real Height() const { return _height; }
		Real& Height() { return _height; }

		Rectangle()
			: _width(0.0)
			, _height(0.0) {}
		Rectangle(Real width, Real height)
			: _width(width)
			, _height(height) {}

		bool IsValid() const { return _width > 0 && _height > 0; }
		bool IsSquare(Real eps = Defaults::ShapePropertyTolerance) const { return std::abs(_width - _height) < eps; }

		Real Area() const { return _width * _height; }
		Real Perimeter() const { return 2.0 * (_width + _height); }
		Real Diagonal() const { return sqrt(_width * _width + _height * _height); }

		// Aspect ratio (width/height)
		Real AspectRatio() const { return _width / _height; }

		// Inscribed/circumscribed circles
		Real InscribedCircleRadius() const { return std::min(_width, _height) / 2.0; }
		Real CircumscribedCircleRadius() const { return Diagonal() / 2.0; }

		static Rectangle Square(Real side) { return Rectangle(side, side); }
		static Rectangle GoldenRectangle(Real width) { return Rectangle(width, width / Constants::GoldenRatio); }
		static Rectangle FromDiagonalAndRatio(Real diagonal, Real ratio) {
			Real h = diagonal / sqrt(1.0 + ratio * ratio);
			Real w = h * ratio;
			return Rectangle(w, h);
		}
	};

	/// @brief Parallelogram defined by two sides and included angle

	class Parallelogram {
	private:
		Real _a;	 // one pair of parallel sides
		Real _b;	 // other pair of parallel sides
		Real _angle; // included angle in radians

	public:
		Real SideA() const { return _a; }
		Real& SideA() { return _a; }
		Real SideB() const { return _b; }
		Real& SideB() { return _b; }
		Real Angle() const { return _angle; }
		Real& Angle() { return _angle; }
		Real AngleDegrees() const { return _angle * 180.0 / Constants::PI; }

		Parallelogram()
			: _a(0.0)
			, _b(0.0)
			, _angle(0.0) {}
		Parallelogram(Real a, Real b, Real angle)
			: _a(a)
			, _b(b)
			, _angle(angle) {}

		bool IsValid() const { return _a > 0 && _b > 0 && _angle > 0 && _angle < Constants::PI; }
		bool IsRectangle(Real eps = Defaults::ShapePropertyTolerance) const { return std::abs(_angle - Constants::PI / 2) < eps; }
		bool IsRhombus(Real eps = Defaults::ShapePropertyTolerance) const { return std::abs(_a - _b) < eps; }
		bool IsSquare(Real eps = Defaults::ShapePropertyTolerance) const { return IsRectangle(eps) && IsRhombus(eps); }

		Real Area() const { return _a * _b * sin(_angle); }
		Real Perimeter() const { return 2.0 * (_a + _b); }

		// Heights
		Real HeightToA() const { return _b * sin(_angle); }
		Real HeightToB() const { return _a * sin(_angle); }

		// Diagonals (using law of cosines)
		Real DiagonalP() const { return sqrt(_a * _a + _b * _b - 2 * _a * _b * cos(_angle)); }
		Real DiagonalQ() const { return sqrt(_a * _a + _b * _b + 2 * _a * _b * cos(_angle)); }

		// The two angles
		Real AcuteAngle() const { return std::min(_angle, Constants::PI - _angle); }
		Real ObtuseAngle() const { return std::max(_angle, Constants::PI - _angle); }

		static Parallelogram FromDiagonals(Real p, Real q, Real angleBetween) {
			// Diagonals bisect each other; use parallelogram law
			Real a = 0.5 * sqrt(p * p + q * q + 2 * p * q * cos(angleBetween));
			Real b = 0.5 * sqrt(p * p + q * q - 2 * p * q * cos(angleBetween));
			Real angle = acos((a * a + b * b - 0.25 * p * p) / (2 * a * b));
			return Parallelogram(a, b, angle);
		}
	};

	/// @brief Rhombus defined by side length and one diagonal (or both diagonals)

	class Rhombus {
	private:
		Real _side;
		Real _d1; // first diagonal
		Real _d2; // second diagonal

	public:
		Real Side() const { return _side; }
		Real& Side() { return _side; }
		Real Diagonal1() const { return _d1; }
		Real& Diagonal1() { return _d1; }
		Real Diagonal2() const { return _d2; }
		Real& Diagonal2() { return _d2; }

		Rhombus()
			: _side(0.0)
			, _d1(0.0)
			, _d2(0.0) {}

		// From diagonals (diagonals of rhombus are perpendicular bisectors)
		Rhombus(Real d1, Real d2)
			: _d1(d1)
			, _d2(d2) {
			_side = sqrt(d1 * d1 + d2 * d2) / 2.0;
		}

		static Rhombus FromSideAndDiagonal(Real side, Real diagonal) {
			Real halfD1 = diagonal / 2.0;
			Real halfD2 = sqrt(side * side - halfD1 * halfD1);
			return Rhombus(diagonal, 2.0 * halfD2);
		}

		static Rhombus FromSideAndAngle(Real side, Real angle) {
			Real d1 = 2.0 * side * sin(angle / 2.0);
			Real d2 = 2.0 * side * cos(angle / 2.0);
			return Rhombus(d1, d2);
		}

		bool IsValid() const { return _side > 0 && _d1 > 0 && _d2 > 0; }
		bool IsSquare(Real eps = Defaults::ShapePropertyTolerance) const { return std::abs(_d1 - _d2) < eps; }

		Real Area() const { return 0.5 * _d1 * _d2; }
		Real Perimeter() const { return 4.0 * _side; }

		// Angles (acute and obtuse)
		Real AcuteAngle() const { return 2.0 * atan(std::min(_d1, _d2) / std::max(_d1, _d2)); }
		Real ObtuseAngle() const { return Constants::PI - AcuteAngle(); }

		// Height (perpendicular distance between parallel sides)
		Real Height() const { return Area() / _side; }

		// Inscribed circle (rhombus has incircle tangent to all sides)
		Real Inradius() const { return (_d1 * _d2) / (2.0 * sqrt(_d1 * _d1 + _d2 * _d2)); }
		Circle InscribedCircle() const { return Circle(Inradius()); }
	};

	/// @brief Trapezoid defined by parallel sides and height

	class Trapezoid {
	private:
		Real _a;	  // one parallel side
		Real _b;	  // other parallel side
		Real _height; // perpendicular distance between parallel sides

	public:
		Real ParallelSideA() const { return _a; }
		Real& ParallelSideA() { return _a; }
		Real ParallelSideB() const { return _b; }
		Real& ParallelSideB() { return _b; }
		Real Height() const { return _height; }
		Real& Height() { return _height; }

		Trapezoid()
			: _a(0.0)
			, _b(0.0)
			, _height(0.0) {}
		Trapezoid(Real a, Real b, Real height)
			: _a(a)
			, _b(b)
			, _height(height) {}

		bool IsValid() const { return _a > 0 && _b > 0 && _height > 0; }
		bool IsParallelogram(Real eps = Defaults::ShapePropertyTolerance) const { return std::abs(_a - _b) < eps; }

		Real Area() const { return 0.5 * (_a + _b) * _height; }

		// Median (midsegment) - line connecting midpoints of non-parallel sides
		Real Median() const { return 0.5 * (_a + _b); }

		// For isosceles trapezoid with known leg length
		static Trapezoid Isosceles(Real base, Real top, Real leg) {
			Real height = sqrt(leg * leg - 0.25 * POW2(base - top));
			return Trapezoid(base, top, height);
		}
	};

} // namespace MML
#endif
