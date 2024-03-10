# Geometry 2D & 3D classes

Defined geometrical objects:
- Line2D
- SegmentLine2D
- Triangle2D
- Polygon2D
- Line3D
- SegmentLine3D
- Plane3D
- Triangle3D

## Geometry 2D

**Line2D**

~~~C++
class Line2D
{
private:
	Point2Cartesian _point;
	Vector2Cartesian _direction; // unit vector in line direction

public:
	Line2D(const Point2Cartesian& pnt, const Vector2Cartesian dir);
	Line2D(const Point2Cartesian& a, const Point2Cartesian& b);

	Point2Cartesian     StartPoint() const { return _point; }
	Point2Cartesian& StartPoint() { return _point; }

	Vector2Cartesian    Direction() const { return _direction; }
	Vector2Cartesian& Direction() { return _direction; }

	Point2Cartesian PointOnLine(Real t);
};
~~~

**SegmentLine2D**

~~~C++
class SegmentLine2D
{
private:
	Point2Cartesian _point1;
	Point2Cartesian _point2;

public:
	SegmentLine2D(Point2Cartesian pnt1, Point2Cartesian pnt2) : _point1(pnt1), _point2(pnt2);
	SegmentLine2D(const Point2Cartesian& pnt1, const Vector2Cartesian& direction, Real t) : _point1(pnt1);

	Point2Cartesian     StartPoint() const { return _point1; }
	Point2Cartesian& StartPoint() { return _point1; }

	Point2Cartesian     EndPoint()  const { return _point2; }
	Point2Cartesian& EndPoint() { return _point2; }

	Point2Cartesian PointOnSegment(Real t);

	Real                Length()    const { return _point1.Dist(_point2); }
	Vector2Cartesian    Direction() const { return Vector2Cartesian(_point1, _point2); }
};
~~~

**Triangle2D**

~~~C++
class Triangle2D
{
private:
	Point2Cartesian _pnt1, _pnt2, _pnt3;

public:
	Triangle2D(Point2Cartesian pnt1, Point2Cartesian pnt2, Point2Cartesian pnt3) : _pnt1(pnt1), _pnt2(pnt2), _pnt3(pnt3)

	Point2Cartesian  Pnt1() const { return _pnt1; }
	Point2Cartesian& Pnt1() { return _pnt1; }
	Point2Cartesian  Pnt2() const { return _pnt2; }
	Point2Cartesian& Pnt2() { return _pnt2; }
	Point2Cartesian  Pnt3() const { return _pnt3; }
	Point2Cartesian& Pnt3() { return _pnt3; }

	Real Area() const;
};
~~~

**Polygon2D**

~~~C++
~~~

## Geometry 3D

**Line3D**

~~~C++
class Line3D
{
private:
	Point3Cartesian  _point;
	Vector3Cartesian _direction;

public:
	Line3D() {}
	Line3D(const Point3Cartesian& pnt, const Vector3Cartesian dir)
	Line3D(const Point3Cartesian& a, const Point3Cartesian& b)

	Point3Cartesian   StartPoint() const { return _point; }
	Point3Cartesian&  StartPoint() { return _point; }

	Vector3Cartesian  Direction() const { return _direction; }
	Vector3Cartesian& Direction() { return _direction; }

	Point3Cartesian PointOnLine(Real t) const { return _point + t * _direction; }

	bool IsPerpendicular(const Line3D& b) const;
	bool IsParallel(const Line3D& b) const;

	Real Dist(const Point3Cartesian& pnt) const;

	Point3Cartesian NearestPoint(const Point3Cartesian& pnt) const;
};
~~~

**SegmentLine3D**

~~~C++
class SegmentLine3D
{
private:
	Point3Cartesian _point1;
	Point3Cartesian _point2;

public:
	SegmentLine3D(Point3Cartesian pnt1, Point3Cartesian pnt2) : _point1(pnt1), _point2(pnt2)
	SegmentLine3D(Point3Cartesian pnt1, Vector3Cartesian direction, Real t)

	Point3Cartesian   StartPoint() const { return _point1; }
	Point3Cartesian& StartPoint() { return _point1; }

	Point3Cartesian   EndPoint() const { return _point2; }
	Point3Cartesian& EndPoint() { return _point2; }

	Point3Cartesian PointOnSegment(Real t);

	Real                Length()    const { return _point1.Dist(_point2); }
	Vector3Cartesian    Direction() const { return Vector3Cartesian(_point1, _point2); }
};
~~~

**Plane3D**

~~~C++
class Plane3D
{
private:
	Real _A, _B, _C, _D;

public:
	Plane3D(const Point3Cartesian& a, const Vector3Cartesian& normal)
	Plane3D(const Point3Cartesian& a, const Point3Cartesian& b, const Point3Cartesian& c)
		: Plane3D(a, VectorProd(Vector3Cartesian(a, b), Vector3Cartesian(a, c))) { }

	Plane3D(Real alpha, Real beta, Real gamma, Real d)      // Hesseov (normalni) oblik
	Plane3D(Real seg_x, Real seg_y, Real seg_z)

	static Plane3D GetXYPlane() { return Plane3D(Point3Cartesian(0, 0, 0), Vector3Cartesian(0, 0, 1)); }
	static Plane3D GetXZPlane() { return Plane3D(Point3Cartesian(0, 0, 0), Vector3Cartesian(0, 1, 0)); }
	static Plane3D GetYZPlane() { return Plane3D(Point3Cartesian(0, 0, 0), Vector3Cartesian(1, 0, 0)); }

	Real  A() const { return _A; }
	Real& A() { return _A; }
	Real  B() const { return _B; }
	Real& B() { return _B; }
	Real  C() const { return _C; }
	Real& C() { return _C; }
	Real  D() const { return _D; }
	Real& D() { return _D; }

	Point3Cartesian GetPointOnPlane() const;

	Vector3Cartesian Normal() const { return Vector3Cartesian(_A, _B, _C); }

	void GetCoordAxisSegments(Real& outseg_x, Real& outseg_y, Real& outseg_z);

	bool IsPointOnPlane(const Point3Cartesian& pnt, Real defEps = 1e-15) const;
	Real DistToPoint(const Point3Cartesian& pnt) const;

	Point3Cartesian ProjectionToPlane(const Point3Cartesian& pnt) const;

	bool IsLineOnPlane(const Line3D& line) const;

	Real AngleToLine(const Line3D& line) const;

	bool IntersectionWithLine(const Line3D& line, Point3Cartesian& out_inter_pnt) const;

	bool IsParallelToPlane(const Plane3D& plane) const;
	bool IsPerpendicularToPlane(const Plane3D& plane) const;
	Real AngleToPlane(const Plane3D& plane) const;
	Real DistToPlane(const Plane3D& plane) const;
	bool IntersectionWithPlane(const Plane3D& plane, Line3D& out_inter_line) const;
};
~~~

**Triangle3D**

~~~C++
class Triangle3D
{
private:
	Point3Cartesian _pnt1, _pnt2, _pnt3;
public:
	Triangle3D(Point3Cartesian pnt1, Point3Cartesian pnt2, Point3Cartesian pnt3)
		: _pnt1(pnt1), _pnt2(pnt2), _pnt3(pnt3)
	{}

	Point3Cartesian& Pnt1() { return _pnt1; }
	Point3Cartesian& Pnt2() { return _pnt2; }
	Point3Cartesian& Pnt3() { return _pnt3; }
};
~~~

## Set of classes used for modeling surfaces and 3D bodies

**TriangleSurface3D**

**RectSurface3D**

**IntegrableVolume3D**

**SolidSurfaces3D**

**ComposedSolidSurfaces3D**

**TriangleSurface3D**


