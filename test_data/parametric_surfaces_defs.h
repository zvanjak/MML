#if !defined MML_PARAMETRIC_SURFACES_DEFS_H
#define MML_PARAMETRIC_SURFACES_DEFS_H

#include <cmath>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/VectorN.h"
#include "base/Matrix.h"
#include "base/Function.h"
#endif

namespace MML::TestBeds
{
    /*******************************************************************************************************************
     * PARAMETRIC SURFACE DEFINITIONS FOR TESTING
     * 
     * This file contains definitions of 3D parametric surfaces for testing surface analysis algorithms.
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
        const Real R = REAL(2.0), r = REAL(1.0);
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
        Real halfU = u / REAL(2.0);
        return VectorN<Real,3>{ (REAL(1.0) + v/REAL(2.0) * std::cos(halfU)) * std::cos(u),
                                (REAL(1.0) + v/REAL(2.0) * std::cos(halfU)) * std::sin(u),
                                v/REAL(2.0) * std::sin(halfU) };
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
        Real r = REAL(4.0) * (REAL(1.0) - std::cos(u) / REAL(2.0));
        if (u < Constants::PI) {
            return VectorN<Real,3>{ 
                REAL(6.0) * std::cos(u) * (REAL(1.0) + std::sin(u)) + r * std::cos(u) * std::cos(v),
                REAL(16.0) * std::sin(u) + r * std::sin(u) * std::cos(v),
                r * std::sin(v) };
        } else {
            return VectorN<Real,3>{ 
                REAL(6.0) * std::cos(u) * (REAL(1.0) + std::sin(u)) + r * std::cos(v + static_cast<Real>(Constants::PI)),
                REAL(16.0) * std::sin(u),
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
        Real z = (std::abs(cosU) > REAL(0.01) && std::abs(cosV) > REAL(0.01)) 
                 ? std::log(std::abs(cosV / cosU)) : REAL(0.0);
        return VectorN<Real,3>{ u, v, z };
    }

    // ========== 14. DINI'S SURFACE ==========
    // Twisted pseudosphere - surface of constant negative curvature
    // x = a*cos(u)*sin(v), y = a*sin(u)*sin(v), z = a*(cos(v) + ln(tan(v/2))) + b*u
    static VectorN<Real,3> DiniSurface(Real u, Real v) {
        const Real a = REAL(1.0), b = REAL(0.2);
        Real sinV = std::sin(v);
        Real cosV = std::cos(v);
        Real tanHalfV = std::tan(v / REAL(2.0));
        Real logTerm = (tanHalfV > REAL(0.001)) ? std::log(tanHalfV) : REAL(-6.9);  // clamp for small v
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
            (x2 - y2) / REAL(2.0) + xyz,
            xyz + (y2 - z2) / REAL(2.0),
            (z2 - x2) / REAL(2.0) + xyz
        };
    }

    // ========== 16. ROMAN SURFACE (STEINER SURFACE) ==========
    // Another immersion of RP² with tetrahedral symmetry
    static VectorN<Real,3> RomanSurface(Real u, Real v) {
        Real sinU = std::sin(u), cosU = std::cos(u);
        Real sinV = std::sin(v), cosV = std::cos(v);
        Real sin2V = std::sin(REAL(2.0) * v);
        
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
        Real sin2U = std::sin(REAL(2.0) * u);
        
        return VectorN<Real,3>{ 
            sinU * sin2U / REAL(2.0),
            sinU * cosU,
            cosU * cosV
        };
    }

    // ========== 18. HENNEBERG'S SURFACE ==========
    // Minimal surface discovered by Henneberg (1875)
    static VectorN<Real,3> HennebergSurface(Real u, Real v) {
        Real sinhU = std::sinh(u), coshU = std::cosh(u);
        Real cosV = std::cos(v), sinV = std::sin(v);
        Real cos2V = std::cos(REAL(2.0) * v), sin2V = std::sin(REAL(2.0) * v);
        Real cos3V = std::cos(REAL(3.0) * v), sin3V = std::sin(REAL(3.0) * v);
        
        return VectorN<Real,3>{ 
            REAL(2.0) * sinhU * cosV - (REAL(2.0)/REAL(3.0)) * std::sinh(REAL(3.0)*u) * cos3V,
            REAL(2.0) * sinhU * sinV - (REAL(2.0)/REAL(3.0)) * std::sinh(REAL(3.0)*u) * sin3V,
            REAL(2.0) * coshU * cos2V
        };
    }

}  // namespace MML::TestBeds

#endif  // MML_PARAMETRIC_SURFACES_DEFS_H
