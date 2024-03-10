# General geometry classes

Set of simple geometry classes

## Point2Cartesian

~~~C++
class Point2Cartesian
{
private:
	Real _x, _y;

public:
	Real  X() const { return _x; }
	Real& X() { return _x; }
	Real  Y() const { return _y; }
	Real& Y() { return _y; }

	Point2Cartesian() {}
	Point2Cartesian(Real x, Real y) : _x(x), _y(y) {}

	Real Dist(const Point2Cartesian& b) const { return sqrt(POW2(b._x - _x) + POW2(b._y - _y)); }
};	
~~~

## Point2Polar

~~~C++
class Point2Polar
{
private:
	Real _r, _phi;

public:
	Real  R() const { return _r; }
	Real& R() { return _r; }
	Real  Phi() const { return _phi; }
	Real& Phi() { return _phi; }

	Point2Polar() {}
	Point2Polar(Real r, Real phi) : _r(r), _phi(phi) {}

	Real Dist(const Point2Polar& b) const { return sqrt(R() * R() + b.R() * b.R() - 2 * R() * b.R() * cos(b.Phi() - Phi())); }
};
~~~

## Point3Cartesian

~~~C++
class Point3Cartesian
{
private:
	Real _x, _y, _z;

public:
	Real  X() const { return _x; }
	Real& X() { return _x; }
	Real  Y() const { return _y; }
	Real& Y() { return _y; }
	Real  Z() const { return _z; }
	Real& Z() { return _z; }

	Point3Cartesian() {}
	Point3Cartesian(Real x, Real y, Real z) : _x(x), _y(y), _z(z) {}

	Real Dist(const Point3Cartesian& b) const { return sqrt(POW2(b._x - _x) + POW2(b._y - _y) + POW2(b._z - _z)); }

	Point3Cartesian operator+(const Point3Cartesian& b) const { return Point3Cartesian(_x + b._x, _y + b._y, _z + b._z); }
	Point3Cartesian operator-(const Point3Cartesian& b) const { return Point3Cartesian(_x - b._x, _y - b._y, _z - b._z); }

	friend Point3Cartesian operator*(const Point3Cartesian& a, Real b) { return Point3Cartesian(a._x * b, a._y * b, a._z * b); }
	friend Point3Cartesian operator*(Real a, const Point3Cartesian& b) { return Point3Cartesian(a * b._x, a * b._y, a * b._z); }
	friend Point3Cartesian operator/(const Point3Cartesian& a, Real b) { return Point3Cartesian(a._x / b, a._y / b, a._z / b); }
};
~~~

## Triangle

~~~C++
class Triangle
{
private:
	Real _a, _b, _c;

public:
	Real  A() const { return _a; }
	Real& A() { return _a; }
	Real  B() const { return _b; }
	Real& B() { return _b; }
	Real  C() const { return _c; }
	Real& C() { return _c; }

	Triangle() {}
	Triangle(Real a, Real b, Real c) : _a(a), _b(b), _c(c) {}

	Real Area() const
	{
		Real s = (_a + _b + _c) / 2.0;
		return sqrt(s * (s - _a) * (s - _b) * (s - _c));
	}
	bool IsRight() const
	{
		return ( POW2(_a) + POW2(_b) == POW2(_c)) || 
						(POW2(_a) + POW2(_c) == POW2(_b)) || 
						(POW2(_b) + POW2(_c) == POW2(_a) );
	}
	bool IsIsosceles() const
	{
		return (_a == _b) || (_a == _c) || (_b == _c);
	}
	bool IsEquilateral() const
	{
		return (_a == _b) && (_a == _c);
	}
};
~~~

