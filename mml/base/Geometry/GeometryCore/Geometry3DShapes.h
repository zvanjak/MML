///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry3DShapes.h                                                  ///
///  Description: Pure 3D geometric shapes - coordinate-free intrinsic geometry       ///
///               (Sphere, Cylinder, Cone, Frustum, Tetrahedron, Spheroid, Torus)     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GEOMETRY_3D_SHAPES_H
#define MML_GEOMETRY_3D_SHAPES_H

#include <algorithm>
#include <cmath>

#include "mml/MMLBase.h"
#include "mml/base/Geometry/GeometryCore/Geometry2DShapes.h"


namespace MML {

	//=============================================================================
	// 3D PURE GEOMETRIC SHAPES - Coordinate-free intrinsic geometry
	//=============================================================================

	/// @brief Sphere defined by radius (coordinate-free pure geometry)
	/// Note: Named SphereGeom to avoid conflict with Surfaces::Sphere parametric surface.

	class SphereGeom {
	private:
		Real _radius;

	public:
		Real Radius() const { return _radius; }
		Real& Radius() { return _radius; }
		Real Diameter() const { return 2.0 * _radius; }

		SphereGeom()
			: _radius(0.0) {}
		explicit SphereGeom(Real radius)
			: _radius(radius) {}

		bool IsValid() const { return _radius > 0; }

		Real Volume() const { return (4.0 / 3.0) * Constants::PI * _radius * _radius * _radius; }
		Real SurfaceArea() const { return 4.0 * Constants::PI * _radius * _radius; }

		// Great circle (equator)
		Circle GreatCircle() const { return Circle(_radius); }

		// Spherical cap for given height h
		Real CapArea(Real h) const { return 2.0 * Constants::PI * _radius * h; }
		Real CapVolume(Real h) const { return Constants::PI * h * h * (3.0 * _radius - h) / 3.0; }

		// Spherical zone (between two parallel planes)
		Real ZoneArea(Real h) const { return 2.0 * Constants::PI * _radius * h; }

		static SphereGeom FromVolume(Real vol) { return SphereGeom(cbrt(3.0 * vol / (4.0 * Constants::PI))); }
		static SphereGeom FromSurfaceArea(Real area) { return SphereGeom(sqrt(area / (4.0 * Constants::PI))); }
		static SphereGeom Circumscribed(const Rectangle& rect, Real depth) {
			Real d = sqrt(rect.Width() * rect.Width() + rect.Height() * rect.Height() + depth * depth);
			return SphereGeom(d / 2.0);
		}
	};

	/// @brief CylinderSurface defined by radius and height (coordinate-free pure geometry)
	/// Note: Named CylinderGeom to avoid conflict with Surfaces::CylinderSurface parametric surface.

	class CylinderGeom {
	private:
		Real _radius;
		Real _height;

	public:
		Real Radius() const { return _radius; }
		Real& Radius() { return _radius; }
		Real Height() const { return _height; }
		Real& Height() { return _height; }
		Real Diameter() const { return 2.0 * _radius; }

		CylinderGeom()
			: _radius(0.0)
			, _height(0.0) {}
		CylinderGeom(Real radius, Real height)
			: _radius(radius)
			, _height(height) {}

		bool IsValid() const { return _radius > 0 && _height > 0; }

		Real Volume() const { return Constants::PI * _radius * _radius * _height; }

		Real LateralArea() const { return 2.0 * Constants::PI * _radius * _height; }
		Real BaseArea() const { return Constants::PI * _radius * _radius; }
		Real SurfaceArea() const { return LateralArea() + 2.0 * BaseArea(); }

		Circle Base() const { return Circle(_radius); }

		static CylinderGeom FromVolume(Real volume, Real radius) {
			Real h = volume / (Constants::PI * radius * radius);
			return CylinderGeom(radius, h);
		}
	};

	/// @brief Cone defined by base radius and height (coordinate-free pure geometry)
	/// Note: Named ConeGeom to avoid conflict with Surfaces::Cone parametric surface.

	class ConeGeom {
	private:
		Real _radius;
		Real _height;

	public:
		Real Radius() const { return _radius; }
		Real& Radius() { return _radius; }
		Real Height() const { return _height; }
		Real& Height() { return _height; }

		ConeGeom()
			: _radius(0.0)
			, _height(0.0) {}
		ConeGeom(Real radius, Real height)
			: _radius(radius)
			, _height(height) {}

		bool IsValid() const { return _radius > 0 && _height > 0; }

		Real SlantHeight() const { return sqrt(_radius * _radius + _height * _height); }
		Real Volume() const { return Constants::PI * _radius * _radius * _height / 3.0; }

		Real LateralArea() const { return Constants::PI * _radius * SlantHeight(); }
		Real BaseArea() const { return Constants::PI * _radius * _radius; }
		Real SurfaceArea() const { return LateralArea() + BaseArea(); }

		// Half-angle at apex
		Real ApexAngle() const { return atan(_radius / _height); }

		Circle Base() const { return Circle(_radius); }

		static ConeGeom FromSlantHeight(Real radius, Real slant) {
			Real h = sqrt(slant * slant - radius * radius);
			return ConeGeom(radius, h);
		}
	};

	/// @brief Frustum (truncated cone) defined by two radii and height

	class Frustum {
	private:
		Real _r1; // bottom radius
		Real _r2; // top radius
		Real _height;

	public:
		Real BottomRadius() const { return _r1; }
		Real& BottomRadius() { return _r1; }
		Real TopRadius() const { return _r2; }
		Real& TopRadius() { return _r2; }
		Real Height() const { return _height; }
		Real& Height() { return _height; }

		Frustum()
			: _r1(0.0)
			, _r2(0.0)
			, _height(0.0) {}
		Frustum(Real r1, Real r2, Real height)
			: _r1(r1)
			, _r2(r2)
			, _height(height) {}

		bool IsValid() const { return _r1 > 0 && _r2 > 0 && _height > 0; }

		Real SlantHeight() const {
			Real dr = _r1 - _r2;
			return sqrt(dr * dr + _height * _height);
		}

		Real Volume() const { return Constants::PI * _height * (_r1 * _r1 + _r1 * _r2 + _r2 * _r2) / 3.0; }

		Real LateralArea() const { return Constants::PI * (_r1 + _r2) * SlantHeight(); }

		Real SurfaceArea() const { return LateralArea() + Constants::PI * (_r1 * _r1 + _r2 * _r2); }

		Circle BottomBase() const { return Circle(_r1); }
		Circle TopBase() const { return Circle(_r2); }
	};

	/// @brief Tetrahedron defined by six edge lengths (coordinate-free)
	/// Vertices: A, B, C, D
	/// Edges: AB=a, AC=b, AD=c, BC=d, BD=e, CD=f

	class Tetrahedron {
	private:
		Real _ab, _ac, _ad; // edges from vertex A
		Real _bc, _bd, _cd; // edges opposite to edges from A

	public:
		Real AB() const { return _ab; }
		Real& AB() { return _ab; }
		Real AC() const { return _ac; }
		Real& AC() { return _ac; }
		Real AD() const { return _ad; }
		Real& AD() { return _ad; }
		Real BC() const { return _bc; }
		Real& BC() { return _bc; }
		Real BD() const { return _bd; }
		Real& BD() { return _bd; }
		Real CD() const { return _cd; }
		Real& CD() { return _cd; }

		Tetrahedron()
			: _ab(0)
			, _ac(0)
			, _ad(0)
			, _bc(0)
			, _bd(0)
			, _cd(0) {}
		Tetrahedron(Real ab, Real ac, Real ad, Real bc, Real bd, Real cd)
			: _ab(ab)
			, _ac(ac)
			, _ad(ad)
			, _bc(bc)
			, _bd(bd)
			, _cd(cd) {}

		// Regular tetrahedron (all edges equal)
		static Tetrahedron Regular(Real edge) { return Tetrahedron(edge, edge, edge, edge, edge, edge); }

		bool IsRegular(Real eps = Defaults::ShapePropertyTolerance) const {
			return std::abs(_ab - _ac) < eps && std::abs(_ab - _ad) < eps && std::abs(_ab - _bc) < eps && std::abs(_ab - _bd) < eps &&
				   std::abs(_ab - _cd) < eps;
		}

		// Volume using Cayley-Menger determinant
		Real Volume() const {
			Real a2 = _ab * _ab, b2 = _ac * _ac, c2 = _ad * _ad;
			Real d2 = _bc * _bc, e2 = _bd * _bd, f2 = _cd * _cd;

			// Cayley-Menger determinant / 288
			Real term1 = a2 * f2 * (b2 + c2 + d2 + e2 - a2 - f2);
			Real term2 = b2 * e2 * (a2 + c2 + d2 + f2 - b2 - e2);
			Real term3 = c2 * d2 * (a2 + b2 + e2 + f2 - c2 - d2);
			Real term4 = -a2 * b2 * f2 - a2 * c2 * e2 - b2 * c2 * d2;
			Real term5 = -d2 * e2 * f2;

			Real det = term1 + term2 + term3 + term4 + term5;
			return sqrt(std::abs(det)) / 12.0;
		}

		// Surface area (sum of four triangular faces)
		Real SurfaceArea() const {
			Triangle abc(_ab, _ac, _bc); // face ABC
			Triangle abd(_ab, _ad, _bd); // face ABD
			Triangle acd(_ac, _ad, _cd); // face ACD
			Triangle bcd(_bc, _bd, _cd); // face BCD
			return abc.Area() + abd.Area() + acd.Area() + bcd.Area();
		}

		// The four faces as triangles
		Triangle FaceABC() const { return Triangle(_ab, _ac, _bc); }
		Triangle FaceABD() const { return Triangle(_ab, _ad, _bd); }
		Triangle FaceACD() const { return Triangle(_ac, _ad, _cd); }
		Triangle FaceBCD() const { return Triangle(_bc, _bd, _cd); }

		// For regular tetrahedron
		Real Inradius() const { return 3.0 * Volume() / SurfaceArea(); }

		Real Circumradius() const {
			if (!IsRegular())
				return 0; // Only valid for regular
			return _ab * sqrt(6.0) / 4.0;
		}
	};

	/// @brief Spheroid (ellipsoid of revolution) - oblate or prolate

	class Spheroid {
	private:
		Real _equatorial; // equatorial radius (a)
		Real _polar;	  // polar radius (c)

	public:
		Real EquatorialRadius() const { return _equatorial; }
		Real& EquatorialRadius() { return _equatorial; }
		Real PolarRadius() const { return _polar; }
		Real& PolarRadius() { return _polar; }

		Spheroid()
			: _equatorial(0.0)
			, _polar(0.0) {}
		Spheroid(Real equatorial, Real polar)
			: _equatorial(equatorial)
			, _polar(polar) {}

		bool IsValid() const { return _equatorial > 0 && _polar > 0; }
		bool IsSphere(Real eps = Defaults::ShapePropertyTolerance) const { return std::abs(_equatorial - _polar) < eps; }
		bool IsOblate() const { return _equatorial > _polar; }	// flattened at poles (like Earth)
		bool IsProlate() const { return _polar > _equatorial; } // elongated at poles

		Real Volume() const { return (4.0 / 3.0) * Constants::PI * _equatorial * _equatorial * _polar; }

		// Flattening f = (a - c) / a
		Real Flattening() const { return (_equatorial - _polar) / _equatorial; }

		// Eccentricity
		Real Eccentricity() const {
			if (IsOblate())
				return sqrt(1.0 - (_polar * _polar) / (_equatorial * _equatorial));
			else
				return sqrt(1.0 - (_equatorial * _equatorial) / (_polar * _polar));
		}

		// Surface area (exact formula involves elliptic integrals, using approximation)
		Real SurfaceArea() const {
			if (IsSphere())
				return 4.0 * Constants::PI * _equatorial * _equatorial;

			Real a = _equatorial, c = _polar;
			if (IsOblate()) {
				Real e = Eccentricity();
				return 2.0 * Constants::PI * a * a * (1.0 + (1.0 - e * e) / e * atanh(e));
			} else // Prolate
			{
				Real e = Eccentricity();
				return 2.0 * Constants::PI * a * a * (1.0 + c / (a * e) * asin(e));
			}
		}

		static Spheroid Earth() {
			// WGS84 ellipsoid
			return Spheroid(6378137.0, 6356752.3142);
		}
	};

	/// @brief Torus defined by major and minor radii (coordinate-free pure geometry)
	/// Note: Named TorusGeom to avoid conflict with Surfaces::Torus parametric surface.

	class TorusGeom {
	private:
		Real _R; // major radius (center of tube to center of torus)
		Real _r; // minor radius (radius of tube)

	public:
		Real MajorRadius() const { return _R; }
		Real& MajorRadius() { return _R; }
		Real MinorRadius() const { return _r; }
		Real& MinorRadius() { return _r; }

		TorusGeom()
			: _R(0.0)
			, _r(0.0) {}
		TorusGeom(Real major, Real minor)
			: _R(major)
			, _r(minor) {}

		bool IsValid() const { return _R > 0 && _r > 0; }
		bool IsRingTorus() const { return _R > _r; }	// donut shape, hole in middle
		bool IsHornTorus() const { return _R == _r; }	// tangent to axis
		bool IsSpindleTorus() const { return _R < _r; } // self-intersecting

		Real Volume() const { return 2.0 * Constants::PI * Constants::PI * _R * _r * _r; }
		Real SurfaceArea() const { return 4.0 * Constants::PI * Constants::PI * _R * _r; }

		// Inner and outer radii (for ring torus)
		Real InnerRadius() const { return _R - _r; }
		Real OuterRadius() const { return _R + _r; }

		// Cross-section
		Circle TubeCrossSection() const { return Circle(_r); }
	};

} // namespace MML
#endif
