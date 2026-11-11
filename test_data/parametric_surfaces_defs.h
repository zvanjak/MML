#if !defined MML_PARAMETRIC_SURFACES_DEFS_H
#define MML_PARAMETRIC_SURFACES_DEFS_H

#include <cmath>

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Vector/VectorN.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Function.h>
#include <mml/core/Surfaces.h>
#endif

namespace MML::TestBeds
{
    /*******************************************************************************************************************
     * PARAMETRIC SURFACE DEFINITIONS FOR TESTING
     * 
     * This file contains surface functions and shared ShowcaseSurfaces classes.
     * Showcase classes retain their original parameterizations outside the core MML distribution.
     * Each surface includes:
     *   - Position function r(u,v)
     *   - Geometric properties (curvature, area, etc.)
     * 
     * Categories:
     *   - Quadric surfaces: Sphere, Ellipsoid, Paraboloid, Hyperboloid
     *   - Ruled surfaces: Cylinder, Cone, Hyperboloid of one sheet
     *   - Minimal surfaces: Helicoid, Catenoid, Enneper
     *   - Topologically interesting: Torus, Möbius strip, Klein bottle
     *   - Explicit surfaces: z = f(x,y) forms
     *******************************************************************************************************************/

    //==================================================================================================================
    // PARAMETRIC SURFACE FUNCTIONS
    //==================================================================================================================

    // ========== 1. UNIT SPHERE ==========
    // Classic test surface with constant positive curvature
    // x = cos(u)sin(v), y = sin(u)sin(v), z = cos(v)
    // K = 1, H = 1, Area = 4π
    static VectorN<Real,3> UnitSphere(Real u, Real v) {
        return VectorN<Real,3>{ std::cos(u) * std::sin(v), 
                                std::sin(u) * std::sin(v), 
                                std::cos(v) };
    }

    // ========== 2. TORUS ==========
    // (R + r*cos(v))cos(u), (R + r*cos(v))sin(u), r*sin(v)
    // R=2 (major), r=1 (minor), Area = 4π²Rr = 8π²
    // K varies: positive on outside, negative on inside
    static VectorN<Real,3> Torus(Real u, Real v) {
        const Real R = static_cast<MML::Real>(2.0), r = static_cast<MML::Real>(1.0);
        return VectorN<Real,3>{ (R + r * std::cos(v)) * std::cos(u),
                                (R + r * std::cos(v)) * std::sin(u),
                                r * std::sin(v) };
    }

    // ========== 3. CYLINDER ==========
    // Developable surface (K=0), H≠0
    // x = cos(u), y = sin(u), z = v
    // Area = 2π*h (for height h)
    static VectorN<Real,3> Cylinder(Real u, Real v) {
        return VectorN<Real,3>{ std::cos(u), std::sin(u), v };
    }

    // ========== 4. CONE ==========
    // Developable surface (K=0)
    // x = v*cos(u), y = v*sin(u), z = v
    static VectorN<Real,3> Cone(Real u, Real v) {
        return VectorN<Real,3>{ v * std::cos(u), v * std::sin(u), v };
    }

    // ========== 5. HELICOID ==========
    // Minimal surface (H=0), K<0
    // x = v*cos(u), y = v*sin(u), z = u
    // One of only two ruled minimal surfaces (other is plane)
    static VectorN<Real,3> Helicoid(Real u, Real v) {
        return VectorN<Real,3>{ v * std::cos(u), v * std::sin(u), u };
    }

    // ========== 6. CATENOID ==========
    // Minimal surface (H=0), K<0
    // x = cosh(v)*cos(u), y = cosh(v)*sin(u), z = v
    // Surface of revolution of catenary
    static VectorN<Real,3> Catenoid(Real u, Real v) {
        return VectorN<Real,3>{ std::cosh(v) * std::cos(u), 
                                std::cosh(v) * std::sin(u), 
                                v };
    }

    // ========== 7. ENNEPER SURFACE ==========
    // Self-intersecting minimal surface
    // Classic example for testing minimal surface algorithms
    static VectorN<Real,3> EnneperSurface(Real u, Real v) {
        return VectorN<Real,3>{ u - u*u*u/3 + u*v*v,
                                v - v*v*v/3 + v*u*u,
                                u*u - v*v };
    }

    // ========== 8. MÖBIUS STRIP ==========
    // Non-orientable surface
    // x = (1 + v/2*cos(u/2))cos(u), y = (1 + v/2*cos(u/2))sin(u), z = v/2*sin(u/2)
    static VectorN<Real,3> MobiusStrip(Real u, Real v) {
        Real halfU = u / static_cast<MML::Real>(2.0);
        return VectorN<Real,3>{ (static_cast<MML::Real>(1.0) + v/static_cast<MML::Real>(2.0) * std::cos(halfU)) * std::cos(u),
                                (static_cast<MML::Real>(1.0) + v/static_cast<MML::Real>(2.0) * std::cos(halfU)) * std::sin(u),
                                v/static_cast<MML::Real>(2.0) * std::sin(halfU) };
    }

    // ========== 9. HYPERBOLOID OF ONE SHEET ==========
    // Doubly ruled surface, K<0
    // x = cosh(v)*cos(u), y = cosh(v)*sin(u), z = sinh(v)
    static VectorN<Real,3> HyperboloidOneSheet(Real u, Real v) {
        return VectorN<Real,3>{ std::cosh(v) * std::cos(u),
                                std::cosh(v) * std::sin(u),
                                std::sinh(v) };
    }

    // ========== 10. PARABOLOID ==========
    // Bowl shape, K>0 everywhere
    // x = u, y = v, z = u² + v²
    static VectorN<Real,3> Paraboloid(Real u, Real v) {
        return VectorN<Real,3>{ u, v, u*u + v*v };
    }

    // ========== 11. HYPERBOLIC PARABOLOID (SADDLE) ==========
    // Doubly ruled, K<0 everywhere
    // x = u, y = v, z = u² - v²  (saddle point at origin)
    static VectorN<Real,3> HyperbolicParaboloid(Real u, Real v) {
        return VectorN<Real,3>{ u, v, u*u - v*v };
    }

    // ========== 12. KLEIN BOTTLE ==========
    // Non-orientable, self-intersecting immersion (Figure-8 immersion in R³)
    static VectorN<Real,3> KleinBottle(Real u, Real v) {
        Real r = static_cast<MML::Real>(4.0) * (static_cast<MML::Real>(1.0) - std::cos(u) / static_cast<MML::Real>(2.0));
        if (u < Constants::PI) {
            return VectorN<Real,3>{ 
                static_cast<MML::Real>(6.0) * std::cos(u) * (static_cast<MML::Real>(1.0) + std::sin(u)) + r * std::cos(u) * std::cos(v),
                static_cast<MML::Real>(16.0) * std::sin(u) + r * std::sin(u) * std::cos(v),
                r * std::sin(v) };
        } else {
            return VectorN<Real,3>{ 
                static_cast<MML::Real>(6.0) * std::cos(u) * (static_cast<MML::Real>(1.0) + std::sin(u)) + r * std::cos(v + static_cast<Real>(Constants::PI)),
                static_cast<MML::Real>(16.0) * std::sin(u),
                r * std::sin(v) };
        }
    }

    //==================================================================================================================
    // NEW SURFACES - Additional mathematically interesting surfaces
    //==================================================================================================================

    // ========== 13. SCHERK'S SURFACE ==========
    // Doubly periodic minimal surface
    // z = ln(cos(y)/cos(x)) for |x|, |y| < π/2
    static VectorN<Real,3> ScherkSurface(Real u, Real v) {
        // Domain: |u|, |v| < π/2 - ε to avoid singularities
        Real cosU = std::cos(u);
        Real cosV = std::cos(v);
        Real z = (std::abs(cosU) > static_cast<MML::Real>(0.01) && std::abs(cosV) > static_cast<MML::Real>(0.01))
                 ? std::log(std::abs(cosV / cosU)) : static_cast<MML::Real>(0.0);
        return VectorN<Real,3>{ u, v, z };
    }

    // ========== 14. DINI'S SURFACE ==========
    // Twisted pseudosphere - surface of constant negative curvature
    // x = a*cos(u)*sin(v), y = a*sin(u)*sin(v), z = a*(cos(v) + ln(tan(v/2))) + b*u
    static VectorN<Real,3> DiniSurface(Real u, Real v) {
        const Real a = static_cast<MML::Real>(1.0), b = static_cast<MML::Real>(0.2);
        Real sinV = std::sin(v);
        Real cosV = std::cos(v);
        Real tanHalfV = std::tan(v / static_cast<MML::Real>(2.0));
        Real logTerm = (tanHalfV > static_cast<MML::Real>(0.001)) ? std::log(tanHalfV) : static_cast<MML::Real>(-6.9);  // clamp for small v
        return VectorN<Real,3>{ a * std::cos(u) * sinV,
                                a * std::sin(u) * sinV,
                                a * (cosV + logTerm) + b * u };
    }

    // ========== 15. BOY'S SURFACE ==========
    // Immersion of real projective plane in R³ (non-orientable, self-intersecting)
    static VectorN<Real,3> BoysSurface(Real u, Real v) {
        Real x = std::cos(u) * std::sin(v);
        Real y = std::sin(u) * std::sin(v);
        Real z = std::cos(v);
        
        Real x2 = x * x, y2 = y * y, z2 = z * z;
        Real xyz = x * y * z;
        
        // Bryant-Kusner parameterization (simplified)
        return VectorN<Real,3>{ 
            (x2 - y2) / static_cast<MML::Real>(2.0) + xyz,
            xyz + (y2 - z2) / static_cast<MML::Real>(2.0),
            (z2 - x2) / static_cast<MML::Real>(2.0) + xyz
        };
    }

    // ========== 16. ROMAN SURFACE (STEINER SURFACE) ==========
    // Another immersion of RP² with tetrahedral symmetry
    static VectorN<Real,3> RomanSurface(Real u, Real v) {
        Real sinU = std::sin(u), cosU = std::cos(u);
        Real sinV = std::sin(v), cosV = std::cos(v);
        Real sin2V = std::sin(static_cast<MML::Real>(2.0) * v);
        
        return VectorN<Real,3>{ 
            sin2V * cosU * cosU,
            sin2V * sinU * cosU,
            sinV * sinV * cosV
        };
    }

    // ========== 17. CROSS-CAP ==========
    // Yet another non-orientable surface (simplest immersion of RP²)
    static VectorN<Real,3> CrossCap(Real u, Real v) {
        Real sinU = std::sin(u), cosU = std::cos(u);
        Real sinV = std::sin(v), cosV = std::cos(v);
        Real sin2U = std::sin(static_cast<MML::Real>(2.0) * u);
        
        return VectorN<Real,3>{ 
            sinU * sin2U / static_cast<MML::Real>(2.0),
            sinU * cosU,
            cosU * cosV
        };
    }

    // ========== 18. HENNEBERG'S SURFACE ==========
    // Minimal surface discovered by Henneberg (1875)
    static VectorN<Real,3> HennebergSurface(Real u, Real v) {
        Real sinhU = std::sinh(u), coshU = std::cosh(u);
        Real cosV = std::cos(v), sinV = std::sin(v);
        Real cos2V = std::cos(static_cast<MML::Real>(2.0) * v), sin2V = std::sin(static_cast<MML::Real>(2.0) * v);
        Real cos3V = std::cos(static_cast<MML::Real>(3.0) * v), sin3V = std::sin(static_cast<MML::Real>(3.0) * v);
        
        return VectorN<Real,3>{ 
            static_cast<MML::Real>(2.0) * sinhU * cosV - (static_cast<MML::Real>(2.0)/static_cast<MML::Real>(3.0)) * std::sinh(static_cast<MML::Real>(3.0)*u) * cos3V,
            static_cast<MML::Real>(2.0) * sinhU * sinV - (static_cast<MML::Real>(2.0)/static_cast<MML::Real>(3.0)) * std::sinh(static_cast<MML::Real>(3.0)*u) * sin3V,
            static_cast<MML::Real>(2.0) * coshU * cos2V
        };
    }


    // Shared showcase types are fixtures, not part of the MML core API.
    namespace ShowcaseSurfaces {
		/// @brief Monkey saddle z = u·(u²-3w²), classic example with 3 critical directions
		/// @note Hyperbolic surface with vanishing Gaussian curvature at origin
		class MonkeySaddle : public MML::Surfaces::ISurfaceCartesian
		{
		public:
			Real getMinU() const { return -10; }
			Real getMaxU() const { return 10; }
			Real getMinW() const { return -10; }
			Real getMaxW() const { return 10; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{u, w, u* (u * u - 3 * w * w)}; }
		};

		/// @brief Möbius strip - non-orientable surface with one side and one edge
		/// @note u∈[0,2π], w∈[-1,1]. Famous example of non-orientability
		class MobiusStrip : public MML::Surfaces::ISurfaceCartesian
		{
			Real _scale;  // Scale factor for visibility
		public:
			MobiusStrip() : _scale(1.0) {}
			MobiusStrip(Real scale) : _scale(scale) {}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return 2 * Constants::PI; }
			Real getMinW() const { return -1; }
			Real getMaxW() const { return 1; }

			VectorN<Real, 3> operator()(Real u, Real w) const {
				return VectorN<Real, 3>{
					static_cast<Real>(_scale * (1 + w * cos(u / 2)) * cos(u)),
					static_cast<Real>(_scale * (1 + w * cos(u / 2)) * sin(u)),
					static_cast<Real>(_scale * w * sin(u / 2))
				};
			}
		};

		/// @brief Helicoid - minimal ruled surface like a spiral staircase
		/// @note u is radial distance, w is angle. Gaussian curvature K = -c²/(u² + c²)²
		class Helicoid : public MML::Surfaces::ISurfaceCartesian
		{
			Real _pitch;  // Controls vertical rise per revolution
		public:
			Helicoid() : _pitch(1.0) {}
			Helicoid(Real pitch) : _pitch(pitch) {}

			Real getMinU() const { return -30; }
			Real getMaxU() const { return 30; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return 4 * Constants::PI; }

			VectorN<Real, 3> operator()(Real u, Real w) const {
				return VectorN<Real, 3>{static_cast<Real>(u * cos(w)), static_cast<Real>(u * sin(w)), _pitch * w};
			}
		};

		/// @brief Enneper surface - self-intersecting minimal surface with elegant shape
		/// @note Classic example from differential geometry with K = -16/(1+u²+v²)⁴
		class EnneperSurface : public MML::Surfaces::ISurfaceCartesian
		{
			Real _scale;
		public:
			EnneperSurface() : _scale(1.0) {}
			EnneperSurface(Real scale) : _scale(scale) {}

			Real getMinU() const { return -2; }
			Real getMaxU() const { return 2; }
			Real getMinW() const { return -2; }
			Real getMaxW() const { return 2; }

			VectorN<Real, 3> operator()(Real u, Real v) const {
				return VectorN<Real, 3>{
					static_cast<Real>(_scale * (u - u*u*u/3.0 + u*v*v)),
					static_cast<Real>(_scale * (v - v*v*v/3.0 + v*u*u)),
					static_cast<Real>(_scale * (u*u - v*v))
				};
			}
		};

		/// @brief Klein Bottle - non-orientable closed surface (Pinched Torus immersion)
		/// @note Famous topological surface with no distinct "inside" or "outside"
		/// @note u∈[0,2π], v∈[0,2π]. This immersion creates the classic bottle shape
		class KleinBottle : public MML::Surfaces::ISurfaceCartesian
		{
			Real _scale;
		public:
			KleinBottle() : _scale(1.0) {}
			KleinBottle(Real scale) : _scale(scale) {}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return 2 * Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real u, Real v) const {
				// Classic "bottle" immersion parametrization
				// u goes around the main loop, v goes around the tube
				Real cosu = cos(u), sinu = sin(u);
				Real cosv = cos(v), sinv = sin(v);
				Real cos_half_u = cos(u / 2), sin_half_u = sin(u / 2);

				// Radius varies to create the bottle shape
				Real r = 2 - cosv;

				Real x, y, z;
				if (u < Constants::PI) {
					// First half: normal torus-like shape
					x = 6 * cosu * (1 + sinu) + r * cos_half_u * cosv;
					y = 16 * sinu + r * sin_half_u * cosv;
				} else {
					// Second half: tube goes through itself
					x = 6 * cosu * (1 + sinu) - r * cos_half_u * cosv;
					y = 16 * sinu;
				}
				z = r * sinv;

				return VectorN<Real, 3>{ _scale * x, _scale * y, _scale * z };
			}
		};

		/// @brief Dini's Surface - twisted pseudosphere with constant negative curvature
		/// @note Elegant spiral shape, u is length along surface, v is angle
		class DiniSurface : public MML::Surfaces::ISurfaceCartesian
		{
			Real _a;  // Scaling factor
			Real _b;  // Controls twist rate
		public:
			DiniSurface() : _a(1.0), _b(0.15) {}
			DiniSurface(Real a, Real b) : _a(a), _b(b) {}

			Real getMinU() const { return 0.01; }  // Avoid singularity at 0
			Real getMaxU() const { return 4 * Constants::PI; }
			Real getMinW() const { return 0.01; }  // Avoid singularity at 0
			Real getMaxW() const { return 2; }

			VectorN<Real, 3> operator()(Real u, Real v) const {
				return VectorN<Real, 3>{
					static_cast<Real>(_a * cos(u) * sin(v)),
					static_cast<Real>(_a * sin(u) * sin(v)),
					static_cast<Real>(_a * (cos(v) + log(tan(v / 2))) + _b * u)
				};
			}
		};
    } // namespace ShowcaseSurfaces

}  // namespace MML::TestBeds

#endif  // MML_PARAMETRIC_SURFACES_DEFS_H
