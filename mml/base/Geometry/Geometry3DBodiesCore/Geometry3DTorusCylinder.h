///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Geometry3DTorusCylinder.h                                           ///
///  Description: Torus and Cylinder 3D solid classes                                 ///
///               Parametric surface mesh representations                             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file Geometry3DTorusCylinder.h
/// @brief Torus and Cylinder solid body representations.
/// - Torus3D: Donut shape with rectangular surface mesh
/// - Cylinder3D: Cylinder with triangular surface mesh
/// @see Geometry3DBodies.h for the aggregate header

#if !defined MML_GEOMETRY_3D_TORUS_CYLINDER_H
#define MML_GEOMETRY_3D_TORUS_CYLINDER_H

#include <cmath>
#include <sstream>

#include "mml/MMLBase.h"
#include "mml/base/Vector/VectorTypes.h"
#include "mml/base/Geometry/Geometry3D.h"
#include "mml/base/Geometry/Geometry3DBodiesCore/Geometry3DBounding.h"
#include "mml/base/Geometry/Geometry3DBodiesCore/Geometry3DBodyBase.h"

namespace MML {

	/// @brief Torus (donut shape) with rectangular surface mesh.
	/// Parametric surface: (R + r·cos(v))·cos(u), (R + r·cos(v))·sin(u), r·sin(v)
	/// where R is major radius and r is minor (tube) radius.
	/// @note Formulas:
	/// - Volume = 2π²Rr²
	/// - Surface Area = 4π²Rr

	class Torus3D : public BodyWithRectSurfaces {
		Real _R;		  ///< Major radius (center to tube center)
		Real _r;		  ///< Minor radius (tube radius)
		Pnt3Cart _center; ///< Center of the torus
		int _numU;		  ///< Divisions around major circle
		int _numV;		  ///< Divisions around tube
	public:
		/// @name Constructors
		/// @{

		/// @brief Construct torus centered at origin.
		/// @param R Major radius
		/// @param r Minor radius (tube)
		/// @param numU Divisions around major circle (default 20)
		/// @param numV Divisions around tube (default 12)

		Torus3D(Real R, Real r, int numU = 20, int numV = 12)
			: _R(R)
			, _r(r)
			, _center(0, 0, 0)
			, _numU(numU)
			, _numV(numV) {
			constructSurfaces();
		}

		/// @brief Construct torus centered at specified point.
		/// @param R Major radius
		/// @param r Minor radius (tube)
		/// @param center Center point
		/// @param numU Divisions around major circle (default 20)
		/// @param numV Divisions around tube (default 12)

		Torus3D(Real R, Real r, const Pnt3Cart& center, int numU = 20, int numV = 12)
			: _R(R)
			, _r(r)
			, _center(center)
			, _numU(numU)
			, _numV(numV) {
			constructSurfaces();
		}
		/// @}

		/// @name Surface Generation
		/// @{

		/// /** @brief Generate rectangular patches for the torus surface. */

		void constructSurfaces() {
			_surfaces.clear();
			Real du = 2.0 * Constants::PI / _numU;
			Real dv = 2.0 * Constants::PI / _numV;

			for (int i = 0; i < _numU; ++i) {
				for (int j = 0; j < _numV; ++j) {
					Real u1 = i * du;
					Real u2 = (i + 1) * du;
					Real v1 = j * dv;
					Real v2 = (j + 1) * dv;

					// Four corners of this patch
					Pnt3Cart p1 = torusPoint(u1, v1);
					Pnt3Cart p2 = torusPoint(u2, v1);
					Pnt3Cart p3 = torusPoint(u2, v2);
					Pnt3Cart p4 = torusPoint(u1, v2);

					_surfaces.push_back(RectSurface3D(p1, p2, p3, p4));
				}
			}
		}

		/// @brief Compute point on torus surface.
		/// @param u Angle around major circle [0, 2π]
		/// @param v Angle around tube [0, 2π]

		Pnt3Cart torusPoint(Real u, Real v) const {
			Real x = (_R + _r * std::cos(v)) * std::cos(u);
			Real y = (_R + _r * std::cos(v)) * std::sin(u);
			Real z = _r * std::sin(v);
			return Pnt3Cart(_center.X() + x, _center.Y() + y, _center.Z() + z);
		}
		/// @}

		/// @name Containment Test
		/// @{

		bool IsInside(const Pnt3Cart& pnt) const override {
			// Translate to torus-centered coordinates
			Real x = pnt.X() - _center.X();
			Real y = pnt.Y() - _center.Y();
			Real z = pnt.Z() - _center.Z();

			// Distance from point to Z-axis
			Real d = std::sqrt(x * x + y * y);

			// Distance from point to the tube center circle
			Real dist = std::sqrt((d - _R) * (d - _R) + z * z);

			return dist <= _r;
		}
		/// @}

		/// @name IBody Interface
		/// @{

		/// /** @brief Volume = 2π²Rr² */

		Real Volume() const override { return 2.0 * Constants::PI * Constants::PI * _R * _r * _r; }

		/// /** @brief Surface area = 4π²Rr */

		Real SurfaceArea() const override { return 4.0 * Constants::PI * Constants::PI * _R * _r; }

		Pnt3Cart GetCenter() const override { return _center; }

		Box3D GetBoundingBox() const override {
			Real outerRadius = _R + _r;
			Pnt3Cart min(_center.X() - outerRadius, _center.Y() - outerRadius, _center.Z() - _r);
			Pnt3Cart max(_center.X() + outerRadius, _center.Y() + outerRadius, _center.Z() + _r);
			return Box3D(min, max);
		}

		/// /** @brief Bounding sphere radius = R + r (outer radius) */

		BoundingSphere3D GetBoundingSphere() const override {
			Real radius = _R + _r;
			return BoundingSphere3D(_center, radius);
		}

		std::string ToString() const override {
			std::ostringstream oss;
			oss << "Torus3D{Center=(" << _center.X() << ", " << _center.Y() << ", " << _center.Z() << "), MajorRadius=" << _R
				<< ", MinorRadius=" << _r << ", Volume=" << Volume() << ", SurfaceArea=" << SurfaceArea() << "}";
			return oss.str();
		}
		/// @}

		/// @name Accessors
		/// @{

		Real GetMajorRadius() const { return _R; }
		Real GetMinorRadius() const { return _r; }
		int GetNumU() const { return _numU; }
		int GetNumV() const { return _numV; }
		/// @}
	};

	/// @brief Cylinder with triangular surface mesh.
	/// Cylinder aligned with Z-axis, with circular caps at top and bottom.
	/// Surface includes:
	/// - Bottom cap (triangular fan)
	/// - Top cap (triangular fan)
	/// - Lateral surface (rectangular strips as triangle pairs)
	/// @note Formulas:
	/// - Volume = πR²H
	/// - Surface Area = 2πR² + 2πRH

	class Cylinder3D : public BodyWithTriangleSurfaces {
		Real _R;		  ///< Radius
		Real _H;		  ///< Height
		Pnt3Cart _center; ///< Bottom center
		int _numSegments; ///< Divisions around circumference
	public:
		/// @name Constructors
		/// @{

		/// @brief Construct cylinder centered at origin (bottom at z = -H/2).
		/// @param R Radius
		/// @param H Height
		/// @param numSegments Circumferential divisions (default 20)

		Cylinder3D(Real R, Real H, int numSegments = 20)
			: _R(R)
			, _H(H)
			, _center(0, 0, -H / 2)
			, _numSegments(numSegments) {
			constructSurfaces();
		}

		/// @brief Construct cylinder with specified bottom center.
		/// @param R Radius
		/// @param H Height
		/// @param center Bottom center point
		/// @param numSegments Circumferential divisions (default 20)

		Cylinder3D(Real R, Real H, const Pnt3Cart& center, int numSegments = 20)
			: _R(R)
			, _H(H)
			, _center(center)
			, _numSegments(numSegments) {
			constructSurfaces();
		}
		/// @}

		/// @name Surface Generation
		/// @{

		/// /** @brief Generate triangles for caps and lateral surface. */

		void constructSurfaces() {
			_surfaces.clear();
			Real angleStep = 2.0 * Constants::PI / _numSegments;

			Pnt3Cart bottomCenter(_center.X(), _center.Y(), _center.Z());
			Pnt3Cart topCenter(_center.X(), _center.Y(), _center.Z() + _H);

			for (int i = 0; i < _numSegments; ++i) {
				Real angle1 = i * angleStep;
				Real angle2 = (i + 1) * angleStep;

				// Points on bottom circle
				Pnt3Cart b1(_center.X() + _R * std::cos(angle1), _center.Y() + _R * std::sin(angle1), _center.Z());
				Pnt3Cart b2(_center.X() + _R * std::cos(angle2), _center.Y() + _R * std::sin(angle2), _center.Z());

				// Points on top circle
				Pnt3Cart t1(_center.X() + _R * std::cos(angle1), _center.Y() + _R * std::sin(angle1), _center.Z() + _H);
				Pnt3Cart t2(_center.X() + _R * std::cos(angle2), _center.Y() + _R * std::sin(angle2), _center.Z() + _H);

				// Bottom cap (CCW when viewed from above)
				_surfaces.emplace_back(bottomCenter, b2, b1);

				// Top cap (CCW when viewed from below)
				_surfaces.emplace_back(topCenter, t1, t2);

				// Lateral surface (two triangles per segment)
				_surfaces.emplace_back(b1, b2, t2);
				_surfaces.emplace_back(b1, t2, t1);
			}
		}
		/// @}

		/// @name Containment Test
		/// @{

		bool IsInside(const Pnt3Cart& pnt) const override {
			Real x = pnt.X() - _center.X();
			Real y = pnt.Y() - _center.Y();
			Real z = pnt.Z() - _center.Z();

			// Height check
			if (z < 0.0 || z > _H)
				return false;

			// Circular cross-section check
			Real distSq = x * x + y * y;
			return distSq <= _R * _R;
		}
		/// @}

		/// @name IBody Interface
		/// @{

		/// /** @brief Volume = πR²H */

		Real Volume() const override { return Constants::PI * _R * _R * _H; }

		/// /** @brief Surface area = 2πR² + 2πRH */

		Real SurfaceArea() const override { return 2.0 * Constants::PI * _R * _R + 2.0 * Constants::PI * _R * _H; }

		/// /** @brief Return geometric center (middle of cylinder). */

		Pnt3Cart GetCenter() const override { return Pnt3Cart(_center.X(), _center.Y(), _center.Z() + _H / 2.0); }

		Box3D GetBoundingBox() const override {
			Pnt3Cart min(_center.X() - _R, _center.Y() - _R, _center.Z());
			Pnt3Cart max(_center.X() + _R, _center.Y() + _R, _center.Z() + _H);
			return Box3D(min, max);
		}

		/// /** @brief Bounding sphere radius = √(R² + (H/2)²) */

		BoundingSphere3D GetBoundingSphere() const override {
			Pnt3Cart geomCenter = GetCenter();
			Real halfH = _H / 2.0;
			Real radius = std::sqrt(_R * _R + halfH * halfH);
			return BoundingSphere3D(geomCenter, radius);
		}

		std::string ToString() const override {
			std::ostringstream oss;
			Pnt3Cart geomCenter = GetCenter();
			oss << "Cylinder3D{Center=(" << geomCenter.X() << ", " << geomCenter.Y() << ", " << geomCenter.Z() << "), Radius=" << _R
				<< ", Height=" << _H << ", Volume=" << Volume() << ", SurfaceArea=" << SurfaceArea() << "}";
			return oss.str();
		}
		/// @}

		/// @name Accessors
		/// @{

		Real GetRadius() const { return _R; }
		Real GetHeight() const { return _H; }
		Pnt3Cart GetBaseCenter() const { return _center; }
		int GetNumSegments() const { return _numSegments; }
		/// @}
	};

} // namespace MML

#endif
