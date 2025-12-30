#if !defined MML_PARAMETRIC_CURVES_DEFS_H
#define MML_PARAMETRIC_CURVES_DEFS_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/VectorN.h"
#include "base/Function.h"
#include "core/Curves.h"
#endif

namespace MML::TestBeds
{
    /*******************************************************************************************************************
     * PARAMETRIC CURVE DEFINITIONS FOR TESTING
     * 
     * This file contains definitions of 3D parametric curves for testing curve analysis algorithms.
     * Each curve includes:
     *   - Position function r(t)
     *   - First derivative r'(t)
     *   - Second derivative r''(t)
     *   - Curvature κ(t)
     *   - Torsion τ(t)
     * 
     * Categories:
     *   - Basic curves: Helix, Circle, Lines
     *   - Polynomial curves: Twisted cubic, Schaums examples
     *   - Trigonometric curves: Viviani, Lissajous
     *   - Knots: Torus knots, Trefoil knot
     *   - Special curves: Catenary, Cycloid variants
     *******************************************************************************************************************/

    using namespace MML::Curves;

    //==================================================================================================================
    // HELPER: Static functions for curves (avoids lambda-to-function-pointer issues with Real=float)
    //==================================================================================================================

    // ========== TWISTED CUBIC: r(t) = (t, t², t³) ==========
    // Classic space curve used in algebraic geometry
    static VectorN<Real,3> TwistedCubic_pos(Real t) { return VectorN<Real,3>{ t, t*t, t*t*t }; }
    static VectorN<Real,3> TwistedCubic_deriv(Real t) { return VectorN<Real,3>{ REAL(1.0), REAL(2.0)*t, REAL(3.0)*t*t }; }
    static VectorN<Real,3> TwistedCubic_deriv2(Real t) { return VectorN<Real,3>{ REAL(0.0), REAL(2.0), REAL(6.0)*t }; }
    static Real TwistedCubic_curvature(Real t) {
        Real denom = std::pow(REAL(1.0) + REAL(4.0)*t*t + REAL(9.0)*t*t*t*t, REAL(1.5));
        return REAL(2.0) * std::abs(REAL(3.0)*t*t - REAL(1.0)) / denom;
    }
    static Real TwistedCubic_torsion(Real t) {
        Real denom = std::pow(REAL(1.0) + REAL(4.0)*t*t + REAL(9.0)*t*t*t*t, REAL(2.0));
        return REAL(12.0) * t / denom;
    }

    // ========== VIVIANI CURVE ==========
    // Intersection of sphere and cylinder - beautiful closed curve
    static VectorN<Real,3> Viviani_pos(Real t) { return VectorN<Real,3>{ REAL(1.0) + std::cos(t), std::sin(t), REAL(2.0)*std::sin(t/REAL(2.0)) }; }
    static VectorN<Real,3> Viviani_deriv(Real t) { return VectorN<Real,3>{ -std::sin(t), std::cos(t), std::cos(t/REAL(2.0)) }; }
    static VectorN<Real,3> Viviani_deriv2(Real t) { return VectorN<Real,3>{ -std::cos(t), -std::sin(t), REAL(-0.5)*std::sin(t/REAL(2.0)) }; }
    static Real Viviani_curvature(Real t) {
        Real num = std::sqrt(REAL(2.0) * (REAL(1.0) + std::cos(t)));
        Real denom = std::pow(REAL(2.0) + std::cos(t), REAL(1.5));
        return num / denom;
    }
    static Real Viviani_torsion(Real t) {
        Real num = -std::cos(t/REAL(2.0)) * std::sin(t/REAL(2.0));
        Real denom = std::pow(REAL(2.0) + std::cos(t), REAL(2.0));
        return num / denom;
    }

    // ========== TORUS KNOT (3,2) ==========
    // Winds 3 times around the torus tube while going 2 times around the axis
    static VectorN<Real,3> TorusKnot32_pos(Real t) {
        return VectorN<Real,3>{ 
            (REAL(2.0) + std::cos(REAL(2.0)*t)) * std::cos(REAL(3.0)*t), 
            (REAL(2.0) + std::cos(REAL(2.0)*t)) * std::sin(REAL(3.0)*t), 
            std::sin(REAL(2.0)*t) 
        };
    }
    static VectorN<Real,3> TorusKnot32_deriv(Real t) {
        Real dx = -REAL(2.0)*std::sin(REAL(2.0)*t)*std::cos(REAL(3.0)*t) - REAL(3.0)*(REAL(2.0) + std::cos(REAL(2.0)*t))*std::sin(REAL(3.0)*t);
        Real dy = -REAL(2.0)*std::sin(REAL(2.0)*t)*std::sin(REAL(3.0)*t) + REAL(3.0)*(REAL(2.0) + std::cos(REAL(2.0)*t))*std::cos(REAL(3.0)*t);
        Real dz = REAL(2.0)*std::cos(REAL(2.0)*t);
        return VectorN<Real,3>{ dx, dy, dz };
    }
    static VectorN<Real,3> TorusKnot32_deriv2(Real t) {
        Real ddx = -REAL(4.0)*std::cos(REAL(2.0)*t)*std::cos(REAL(3.0)*t) + REAL(6.0)*std::sin(REAL(2.0)*t)*std::sin(REAL(3.0)*t) - REAL(9.0)*(REAL(2.0)+std::cos(REAL(2.0)*t))*std::cos(REAL(3.0)*t);
        Real ddy = -REAL(4.0)*std::cos(REAL(2.0)*t)*std::sin(REAL(3.0)*t) - REAL(6.0)*std::sin(REAL(2.0)*t)*std::cos(REAL(3.0)*t) - REAL(9.0)*(REAL(2.0)+std::cos(REAL(2.0)*t))*std::sin(REAL(3.0)*t);
        Real ddz = -REAL(4.0)*std::sin(REAL(2.0)*t);
        return VectorN<Real,3>{ ddx, ddy, ddz };
    }
    static Real TorusKnot32_curvature(Real t) {
        Real dx = -REAL(2.0)*std::sin(REAL(2.0)*t)*std::cos(REAL(3.0)*t) - REAL(3.0)*(REAL(2.0)+std::cos(REAL(2.0)*t))*std::sin(REAL(3.0)*t);
        Real dy = -REAL(2.0)*std::sin(REAL(2.0)*t)*std::sin(REAL(3.0)*t) + REAL(3.0)*(REAL(2.0)+std::cos(REAL(2.0)*t))*std::cos(REAL(3.0)*t);
        Real dz = REAL(2.0)*std::cos(REAL(2.0)*t);
        Real ddx = -REAL(4.0)*std::cos(REAL(2.0)*t)*std::cos(REAL(3.0)*t) + REAL(6.0)*std::sin(REAL(2.0)*t)*std::sin(REAL(3.0)*t) - REAL(9.0)*(REAL(2.0)+std::cos(REAL(2.0)*t))*std::cos(REAL(3.0)*t);
        Real ddy = -REAL(4.0)*std::cos(REAL(2.0)*t)*std::sin(REAL(3.0)*t) - REAL(6.0)*std::sin(REAL(2.0)*t)*std::cos(REAL(3.0)*t) - REAL(9.0)*(REAL(2.0)+std::cos(REAL(2.0)*t))*std::sin(REAL(3.0)*t);
        Real ddz = -REAL(4.0)*std::sin(REAL(2.0)*t);
        Real v1 = std::sqrt(dx*dx + dy*dy + dz*dz);
        Real num = std::sqrt(POW2(dx*ddy - dy*ddx) + POW2(dx*ddz - dz*ddx) + POW2(dy*ddz - dz*ddy));
        return num / POW3(v1);
    }
    static Real TorusKnot32_torsion(Real t) {
        Real dx = -REAL(2.0)*std::sin(REAL(2.0)*t)*std::cos(REAL(3.0)*t) - REAL(3.0)*(REAL(2.0)+std::cos(REAL(2.0)*t))*std::sin(REAL(3.0)*t);
        Real dy = -REAL(2.0)*std::sin(REAL(2.0)*t)*std::sin(REAL(3.0)*t) + REAL(3.0)*(REAL(2.0)+std::cos(REAL(2.0)*t))*std::cos(REAL(3.0)*t);
        Real dz = REAL(2.0)*std::cos(REAL(2.0)*t);
        Real ddx = -REAL(4.0)*std::cos(REAL(2.0)*t)*std::cos(REAL(3.0)*t) + REAL(6.0)*std::sin(REAL(2.0)*t)*std::sin(REAL(3.0)*t) - REAL(9.0)*(REAL(2.0)+std::cos(REAL(2.0)*t))*std::cos(REAL(3.0)*t);
        Real ddy = -REAL(4.0)*std::cos(REAL(2.0)*t)*std::sin(REAL(3.0)*t) - REAL(6.0)*std::sin(REAL(2.0)*t)*std::cos(REAL(3.0)*t) - REAL(9.0)*(REAL(2.0)+std::cos(REAL(2.0)*t))*std::sin(REAL(3.0)*t);
        Real ddz = -REAL(4.0)*std::sin(REAL(2.0)*t);
        Real dddx = REAL(8.0)*std::sin(REAL(2.0)*t)*std::cos(REAL(3.0)*t) + REAL(12.0)*std::cos(REAL(2.0)*t)*std::sin(REAL(3.0)*t) + REAL(27.0)*(REAL(2.0)+std::cos(REAL(2.0)*t))*std::sin(REAL(3.0)*t);
        Real dddy = REAL(8.0)*std::sin(REAL(2.0)*t)*std::sin(REAL(3.0)*t) - REAL(12.0)*std::cos(REAL(2.0)*t)*std::cos(REAL(3.0)*t) + REAL(27.0)*(REAL(2.0)+std::cos(REAL(2.0)*t))*std::cos(REAL(3.0)*t);
        Real dddz = -REAL(8.0)*std::cos(REAL(2.0)*t);
        Real tx = dy*ddz - dz*ddy;
        Real ty = dz*ddx - dx*ddz;
        Real tz = dx*ddy - dy*ddx;
        Real cross_norm2 = tx*tx + ty*ty + tz*tz;
        return (tx*dddx + ty*dddy + tz*dddz) / cross_norm2;
    }

    // ========== LISSAJOUS (2,3,5) ==========
    // 3D Lissajous figure with frequencies 2:3:5
    static VectorN<Real,3> Lissajous235_pos(Real t) { return VectorN<Real,3>{ std::sin(REAL(2.0)*t), std::sin(REAL(3.0)*t), std::sin(REAL(5.0)*t) }; }
    static VectorN<Real,3> Lissajous235_deriv(Real t) { return VectorN<Real,3>{ REAL(2.0)*std::cos(REAL(2.0)*t), REAL(3.0)*std::cos(REAL(3.0)*t), REAL(5.0)*std::cos(REAL(5.0)*t) }; }
    static VectorN<Real,3> Lissajous235_deriv2(Real t) { return VectorN<Real,3>{ REAL(-4.0)*std::sin(REAL(2.0)*t), REAL(-9.0)*std::sin(REAL(3.0)*t), REAL(-25.0)*std::sin(REAL(5.0)*t) }; }
    static Real Lissajous235_curvature(Real t) {
        Real dx = REAL(2.0)*std::cos(REAL(2.0)*t), dy = REAL(3.0)*std::cos(REAL(3.0)*t), dz = REAL(5.0)*std::cos(REAL(5.0)*t);
        Real ddx = REAL(-4.0)*std::sin(REAL(2.0)*t), ddy = REAL(-9.0)*std::sin(REAL(3.0)*t), ddz = REAL(-25.0)*std::sin(REAL(5.0)*t);
        Real v1 = std::sqrt(dx*dx + dy*dy + dz*dz);
        Real num = std::sqrt(POW2(dx*ddy - dy*ddx) + POW2(dx*ddz - dz*ddx) + POW2(dy*ddz - dz*ddy));
        return num / POW3(v1);
    }
    static Real Lissajous235_torsion(Real t) {
        Real dx = REAL(2.0)*std::cos(REAL(2.0)*t), dy = REAL(3.0)*std::cos(REAL(3.0)*t), dz = REAL(5.0)*std::cos(REAL(5.0)*t);
        Real ddx = REAL(-4.0)*std::sin(REAL(2.0)*t), ddy = REAL(-9.0)*std::sin(REAL(3.0)*t), ddz = REAL(-25.0)*std::sin(REAL(5.0)*t);
        Real dddx = REAL(-8.0)*std::cos(REAL(2.0)*t), dddy = REAL(-27.0)*std::cos(REAL(3.0)*t), dddz = REAL(-125.0)*std::cos(REAL(5.0)*t);
        Real tx = dy*ddz - dz*ddy;
        Real ty = dz*ddx - dx*ddz;
        Real tz = dx*ddy - dy*ddx;
        Real cross_norm2 = tx*tx + ty*ty + tz*tz;
        return (tx*dddx + ty*dddy + tz*dddz) / cross_norm2;
    }

    // ========== TREFOIL KNOT ==========
    // The simplest non-trivial knot - iconic in knot theory
    // r(t) = ((sin(t) + 2*sin(2t)), (cos(t) - 2*cos(2t)), -sin(3t))
    static VectorN<Real,3> TrefoilKnot_pos(Real t) {
        return VectorN<Real,3>{ 
            std::sin(t) + REAL(2.0)*std::sin(REAL(2.0)*t),
            std::cos(t) - REAL(2.0)*std::cos(REAL(2.0)*t),
            -std::sin(REAL(3.0)*t)
        };
    }
    static VectorN<Real,3> TrefoilKnot_deriv(Real t) {
        return VectorN<Real,3>{
            std::cos(t) + REAL(4.0)*std::cos(REAL(2.0)*t),
            -std::sin(t) + REAL(4.0)*std::sin(REAL(2.0)*t),
            REAL(-3.0)*std::cos(REAL(3.0)*t)
        };
    }
    static VectorN<Real,3> TrefoilKnot_deriv2(Real t) {
        return VectorN<Real,3>{
            -std::sin(t) - REAL(8.0)*std::sin(REAL(2.0)*t),
            -std::cos(t) + REAL(8.0)*std::cos(REAL(2.0)*t),
            REAL(9.0)*std::sin(REAL(3.0)*t)
        };
    }
    static Real TrefoilKnot_curvature(Real t) {
        VectorN<Real,3> r1 = TrefoilKnot_deriv(t);
        VectorN<Real,3> r2 = TrefoilKnot_deriv2(t);
        Real dx = r1[0], dy = r1[1], dz = r1[2];
        Real ddx = r2[0], ddy = r2[1], ddz = r2[2];
        Real v1 = std::sqrt(dx*dx + dy*dy + dz*dz);
        Real num = std::sqrt(POW2(dx*ddy - dy*ddx) + POW2(dx*ddz - dz*ddx) + POW2(dy*ddz - dz*ddy));
        return num / POW3(v1);
    }
    static Real TrefoilKnot_torsion(Real t) {
        VectorN<Real,3> r1 = TrefoilKnot_deriv(t);
        VectorN<Real,3> r2 = TrefoilKnot_deriv2(t);
        Real dx = r1[0], dy = r1[1], dz = r1[2];
        Real ddx = r2[0], ddy = r2[1], ddz = r2[2];
        Real dddx = -std::cos(t) - REAL(16.0)*std::cos(REAL(2.0)*t);
        Real dddy = std::sin(t) + REAL(16.0)*std::sin(REAL(2.0)*t);
        Real dddz = REAL(27.0)*std::cos(REAL(3.0)*t);
        Real tx = dy*ddz - dz*ddy;
        Real ty = dz*ddx - dx*ddz;
        Real tz = dx*ddy - dy*ddx;
        Real cross_norm2 = tx*tx + ty*ty + tz*tz;
        if (cross_norm2 < REAL(1e-20)) return REAL(0.0);
        return (tx*dddx + ty*dddy + tz*dddz) / cross_norm2;
    }

    // ========== CONICAL HELIX ==========
    // Helix that spirals outward like a cone: r(t) = (t*cos(t), t*sin(t), t)
    static VectorN<Real,3> ConicalHelix_pos(Real t) {
        return VectorN<Real,3>{ t*std::cos(t), t*std::sin(t), t };
    }
    static VectorN<Real,3> ConicalHelix_deriv(Real t) {
        return VectorN<Real,3>{
            std::cos(t) - t*std::sin(t),
            std::sin(t) + t*std::cos(t),
            REAL(1.0)
        };
    }
    static VectorN<Real,3> ConicalHelix_deriv2(Real t) {
        return VectorN<Real,3>{
            -REAL(2.0)*std::sin(t) - t*std::cos(t),
            REAL(2.0)*std::cos(t) - t*std::sin(t),
            REAL(0.0)
        };
    }
    static Real ConicalHelix_curvature(Real t) {
        VectorN<Real,3> r1 = ConicalHelix_deriv(t);
        VectorN<Real,3> r2 = ConicalHelix_deriv2(t);
        Real dx = r1[0], dy = r1[1], dz = r1[2];
        Real ddx = r2[0], ddy = r2[1], ddz = r2[2];
        Real v1 = std::sqrt(dx*dx + dy*dy + dz*dz);
        Real num = std::sqrt(POW2(dx*ddy - dy*ddx) + POW2(dx*ddz - dz*ddx) + POW2(dy*ddz - dz*ddy));
        return num / POW3(v1);
    }
    static Real ConicalHelix_torsion(Real t) {
        VectorN<Real,3> r1 = ConicalHelix_deriv(t);
        VectorN<Real,3> r2 = ConicalHelix_deriv2(t);
        Real dx = r1[0], dy = r1[1], dz = r1[2];
        Real ddx = r2[0], ddy = r2[1], ddz = r2[2];
        Real dddx = -REAL(3.0)*std::cos(t) + t*std::sin(t);
        Real dddy = -REAL(3.0)*std::sin(t) - t*std::cos(t);
        Real dddz = REAL(0.0);
        Real tx = dy*ddz - dz*ddy;
        Real ty = dz*ddx - dx*ddz;
        Real tz = dx*ddy - dy*ddx;
        Real cross_norm2 = tx*tx + ty*ty + tz*tz;
        if (cross_norm2 < REAL(1e-20)) return REAL(0.0);
        return (tx*dddx + ty*dddy + tz*dddz) / cross_norm2;
    }

    // ========== SPHERICAL SPIRAL (Loxodrome on sphere) ==========
    // Spiral that winds around a sphere from pole to pole
    // r(t) = (cos(t)/cosh(at), sin(t)/cosh(at), tanh(at)) with a = 0.3
    static VectorN<Real,3> SphericalSpiral_pos(Real t) {
        Real a = REAL(0.3);
        Real ch = std::cosh(a*t);
        Real th = std::tanh(a*t);
        return VectorN<Real,3>{ std::cos(t)/ch, std::sin(t)/ch, th };
    }
    static VectorN<Real,3> SphericalSpiral_deriv(Real t) {
        Real a = REAL(0.3);
        Real ch = std::cosh(a*t);
        Real sh = std::sinh(a*t);
        Real ch2 = ch*ch;
        // d/dt[cos(t)/cosh(at)] = -sin(t)/cosh(at) - a*cos(t)*sinh(at)/cosh²(at)
        Real dx = -std::sin(t)/ch - a*std::cos(t)*sh/ch2;
        Real dy = std::cos(t)/ch - a*std::sin(t)*sh/ch2;
        Real dz = a/ch2;  // d/dt[tanh(at)] = a*sech²(at)
        return VectorN<Real,3>{ dx, dy, dz };
    }
    static VectorN<Real,3> SphericalSpiral_deriv2(Real t) {
        Real a = REAL(0.3);
        Real ch = std::cosh(a*t);
        Real sh = std::sinh(a*t);
        Real ch2 = ch*ch;
        Real ch3 = ch2*ch;
        // Second derivatives (computed analytically)
        Real ddx = -std::cos(t)/ch + REAL(2.0)*a*std::sin(t)*sh/ch2 - a*a*std::cos(t)*(ch2 - REAL(2.0)*sh*sh)/ch3;
        Real ddy = -std::sin(t)/ch - REAL(2.0)*a*std::cos(t)*sh/ch2 - a*a*std::sin(t)*(ch2 - REAL(2.0)*sh*sh)/ch3;
        Real ddz = -REAL(2.0)*a*a*sh/ch3;
        return VectorN<Real,3>{ ddx, ddy, ddz };
    }
    static Real SphericalSpiral_curvature(Real t) {
        VectorN<Real,3> r1 = SphericalSpiral_deriv(t);
        VectorN<Real,3> r2 = SphericalSpiral_deriv2(t);
        Real dx = r1[0], dy = r1[1], dz = r1[2];
        Real ddx = r2[0], ddy = r2[1], ddz = r2[2];
        Real v1 = std::sqrt(dx*dx + dy*dy + dz*dz);
        Real num = std::sqrt(POW2(dx*ddy - dy*ddx) + POW2(dx*ddz - dz*ddx) + POW2(dy*ddz - dz*ddy));
        return num / POW3(v1);
    }
    static Real SphericalSpiral_torsion(Real t) {
        // Numerical approximation for torsion (third derivative is complex)
        Real h = REAL(1e-6);
        VectorN<Real,3> r2_p = SphericalSpiral_deriv2(t + h);
        VectorN<Real,3> r2_m = SphericalSpiral_deriv2(t - h);
        VectorN<Real,3> r1 = SphericalSpiral_deriv(t);
        VectorN<Real,3> r2 = SphericalSpiral_deriv2(t);
        Real dx = r1[0], dy = r1[1], dz = r1[2];
        Real ddx = r2[0], ddy = r2[1], ddz = r2[2];
        Real dddx = (r2_p[0] - r2_m[0])/(REAL(2.0)*h);
        Real dddy = (r2_p[1] - r2_m[1])/(REAL(2.0)*h);
        Real dddz = (r2_p[2] - r2_m[2])/(REAL(2.0)*h);
        Real tx = dy*ddz - dz*ddy;
        Real ty = dz*ddx - dx*ddz;
        Real tz = dx*ddy - dy*ddx;
        Real cross_norm2 = tx*tx + ty*ty + tz*tz;
        if (cross_norm2 < REAL(1e-20)) return REAL(0.0);
        return (tx*dddx + ty*dddy + tz*dddz) / cross_norm2;
    }

    // ========== FIGURE-EIGHT KNOT (4₁ knot) ==========
    // Second simplest prime knot after trefoil
    static VectorN<Real,3> FigureEightKnot_pos(Real t) {
        Real c1 = std::cos(t);
        Real c2 = std::cos(REAL(2.0)*t);
        Real s1 = std::sin(t);
        Real s2 = std::sin(REAL(2.0)*t);
        Real s3 = std::sin(REAL(3.0)*t);
        return VectorN<Real,3>{
            (REAL(2.0) + c2) * c1,
            (REAL(2.0) + c2) * s1,
            s3
        };
    }
    static VectorN<Real,3> FigureEightKnot_deriv(Real t) {
        Real c1 = std::cos(t);
        Real c2 = std::cos(REAL(2.0)*t);
        Real s1 = std::sin(t);
        Real s2 = std::sin(REAL(2.0)*t);
        Real c3 = std::cos(REAL(3.0)*t);
        return VectorN<Real,3>{
            -REAL(2.0)*s2*c1 - (REAL(2.0) + c2)*s1,
            -REAL(2.0)*s2*s1 + (REAL(2.0) + c2)*c1,
            REAL(3.0)*c3
        };
    }
    static VectorN<Real,3> FigureEightKnot_deriv2(Real t) {
        Real c1 = std::cos(t);
        Real c2 = std::cos(REAL(2.0)*t);
        Real s1 = std::sin(t);
        Real s2 = std::sin(REAL(2.0)*t);
        Real s3 = std::sin(REAL(3.0)*t);
        return VectorN<Real,3>{
            -REAL(4.0)*c2*c1 + REAL(2.0)*s2*s1 + REAL(2.0)*s2*s1 - (REAL(2.0) + c2)*c1,
            -REAL(4.0)*c2*s1 - REAL(2.0)*s2*c1 - REAL(2.0)*s2*c1 - (REAL(2.0) + c2)*s1,
            REAL(-9.0)*s3
        };
    }
    static Real FigureEightKnot_curvature(Real t) {
        VectorN<Real,3> r1 = FigureEightKnot_deriv(t);
        VectorN<Real,3> r2 = FigureEightKnot_deriv2(t);
        Real dx = r1[0], dy = r1[1], dz = r1[2];
        Real ddx = r2[0], ddy = r2[1], ddz = r2[2];
        Real v1 = std::sqrt(dx*dx + dy*dy + dz*dz);
        Real num = std::sqrt(POW2(dx*ddy - dy*ddx) + POW2(dx*ddz - dz*ddx) + POW2(dy*ddz - dz*ddy));
        return num / POW3(v1);
    }
    static Real FigureEightKnot_torsion(Real t) {
        // Numerical approximation
        Real h = REAL(1e-6);
        VectorN<Real,3> r2_p = FigureEightKnot_deriv2(t + h);
        VectorN<Real,3> r2_m = FigureEightKnot_deriv2(t - h);
        VectorN<Real,3> r1 = FigureEightKnot_deriv(t);
        VectorN<Real,3> r2 = FigureEightKnot_deriv2(t);
        Real dx = r1[0], dy = r1[1], dz = r1[2];
        Real ddx = r2[0], ddy = r2[1], ddz = r2[2];
        Real dddx = (r2_p[0] - r2_m[0])/(REAL(2.0)*h);
        Real dddy = (r2_p[1] - r2_m[1])/(REAL(2.0)*h);
        Real dddz = (r2_p[2] - r2_m[2])/(REAL(2.0)*h);
        Real tx = dy*ddz - dz*ddy;
        Real ty = dz*ddx - dx*ddz;
        Real tz = dx*ddy - dy*ddx;
        Real cross_norm2 = tx*tx + ty*ty + tz*tz;
        if (cross_norm2 < REAL(1e-20)) return REAL(0.0);
        return (tx*dddx + ty*dddy + tz*dddz) / cross_norm2;
    }

}  // namespace MML::TestBeds

#endif  // MML_PARAMETRIC_CURVES_DEFS_H
