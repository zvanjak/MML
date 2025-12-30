#if !defined MML_PARAMETRIC_SURFACE_TEST_BED_H
#define MML_PARAMETRIC_SURFACE_TEST_BED_H

#include <cmath>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/VectorN.h"
#include "base/Matrix.h"
#include "base/Function.h"
#endif

#include "parametric_surfaces_defs.h"

using namespace MML;

namespace MML::TestBeds
{
    /*******************************************************************************************************************
     * PARAMETRIC SURFACES TEST BED
     * 
     * Provides test infrastructure for parametric surface analysis:
     *   - Test data structures (TestParametricSurfaceRect3)
     *   - Collection of test surfaces with known properties
     *   - Explicit surface test bed (z = f(x,y))
     *   - Access by index or name
     * 
     * Uses surface definitions from parametric_surfaces_defs.h
     *******************************************************************************************************************/

    //==================================================================================================================
    // TEST DATA STRUCTURES
    //==================================================================================================================

    /**
     * @brief Test data for a 3D parametric surface over rectangular domain
     * 
     * Contains the surface, its parameter ranges, and geometric properties.
     */
    struct TestParametricSurfaceRect3
    {
        std::string _surfaceName;
        std::string _surfaceExpr;
        Real _u1, _u2;   // u parameter range
        Real _v1, _v2;   // v parameter range

        MML::ParametricSurfaceRect<3> _surface;
        
        // Surface properties for validation
        Real _expectedArea;           // Analytical area (NaN if unknown/complex)
        Real _gaussianCurvature;      // Constant Gaussian curvature (NaN if varies)
        bool _isMinimal;              // Mean curvature H = 0 everywhere
        bool _isDevelopable;          // K = 0 everywhere (can flatten without distortion)

        TestParametricSurfaceRect3( std::string surfaceName, std::string surfaceExpr, 
                                    Real u1, Real u2, Real v1, Real v2,
                                    VectorN<Real,3> (*f1)(Real, Real),
                                    Real expectedArea = std::numeric_limits<Real>::quiet_NaN(),
                                    Real gaussianCurvature = std::numeric_limits<Real>::quiet_NaN(),
                                    bool isMinimal = false, bool isDevelopable = false
                              ) : _surfaceName(surfaceName), _surfaceExpr(surfaceExpr), 
                                  _u1(u1), _u2(u2), _v1(v1), _v2(v2),
                                  _surface(f1),
                                  _expectedArea(expectedArea), _gaussianCurvature(gaussianCurvature),
                                  _isMinimal(isMinimal), _isDevelopable(isDevelopable)
        {}
    };      

    //==================================================================================================================
    // PARAMETRIC SURFACE TEST BED
    //==================================================================================================================

    class ParametricSurfaceTestBed
    {
    public:
        static int getNumTestSurfaces3() { return 18; }  // 12 original + 6 new

        const static TestParametricSurfaceRect3& getTestSurface3(int i) { return _listSurfaces[i]; }

        const static TestParametricSurfaceRect3& getTestSurface3(const std::string &surfaceName)
        {
            for (int i = 0; i < getNumTestSurfaces3(); i++)
            {
                if (_listSurfaces[i]._surfaceName == surfaceName)
                    return _listSurfaces[i];
            }
            throw std::runtime_error("ParametricSurface " + surfaceName + " not found!");
        }

    private:
        const static inline TestParametricSurfaceRect3 _listSurfaces[] = { 
            //==========================================================================================================
            // CLASSIC SURFACES (using static functions from defs.h)
            //==========================================================================================================

            // 1. Unit Sphere - constant positive curvature baseline
            { "Unit Sphere", 
              "r(u,v) = (cos(u)sin(v), sin(u)sin(v), cos(v))", 
              0.0, 2.0 * Constants::PI,  // u: azimuthal angle
              0.0, Constants::PI,        // v: polar angle
              UnitSphere,
              4.0 * Constants::PI,       // Area = 4π
              1.0,                        // K = 1/R² = 1
              false, false },

            // 2. Torus - mixed curvature (positive outside, negative inside)
            { "Torus", 
              "r(u,v) = ((R+r*cos(v))cos(u), (R+r*cos(v))sin(u), r*sin(v)), R=2, r=1", 
              0.0, 2.0 * Constants::PI, 
              0.0, 2.0 * Constants::PI, 
              Torus,
              4.0 * Constants::PI * Constants::PI * 2.0,  // Area = 4π²Rr ≈ 78.96
              std::numeric_limits<Real>::quiet_NaN(),     // K varies
              false, false },

            // 3. Cylinder - developable (K=0), non-minimal
            { "Cylinder", 
              "r(u,v) = (cos(u), sin(u), v)", 
              0.0, 2.0 * Constants::PI, 
              0.0, 2.0,                  // height = 2
              Cylinder,
              4.0 * Constants::PI,       // Area = 2πr*h = 4π
              0.0,                        // K = 0 (developable)
              false, true },

            // 4. Cone - developable (K=0)
            { "Cone", 
              "r(u,v) = (v*cos(u), v*sin(u), v)", 
              0.0, 2.0 * Constants::PI, 
              0.0, 1.0,                  // v from 0 to 1
              Cone,
              std::numeric_limits<Real>::quiet_NaN(),  // Area depends on cone angle
              0.0,                        // K = 0 (developable)
              false, true },

            // 5. Helicoid - minimal surface (H=0), ruled
            { "Helicoid", 
              "r(u,v) = (v*cos(u), v*sin(u), u)", 
              0.0, 2.0 * Constants::PI, 
              -1.0, 1.0, 
              Helicoid,
              std::numeric_limits<Real>::quiet_NaN(),
              std::numeric_limits<Real>::quiet_NaN(),  // K = -1/(1+v²)² varies
              true, false },             // minimal!

            // 6. Catenoid - minimal surface (H=0), surface of revolution
            { "Catenoid", 
              "r(u,v) = (cosh(v)cos(u), cosh(v)sin(u), v)", 
              0.0, 2.0 * Constants::PI, 
              -1.0, 1.0, 
              Catenoid,
              std::numeric_limits<Real>::quiet_NaN(),
              std::numeric_limits<Real>::quiet_NaN(),  // K = -1/cosh⁴(v) varies
              true, false },             // minimal!

            // 7. Enneper Surface - minimal, self-intersecting
            { "Enneper Surface", 
              "r(u,v) = (u - u³/3 + uv², v - v³/3 + vu², u² - v²)", 
              -1.0, 1.0, 
              -1.0, 1.0, 
              EnneperSurface,
              std::numeric_limits<Real>::quiet_NaN(),
              std::numeric_limits<Real>::quiet_NaN(),  // K varies
              true, false },             // minimal!

            // 8. Mobius Strip - non-orientable
            { "Mobius Strip", 
              "r(u,v) = ((1 + v/2*cos(u/2))cos(u), (1 + v/2*cos(u/2))sin(u), v/2*sin(u/2))", 
              0.0, 2.0 * Constants::PI, 
              -0.5, 0.5,                 // width of strip
              MobiusStrip,
              std::numeric_limits<Real>::quiet_NaN(),
              std::numeric_limits<Real>::quiet_NaN(),
              false, false },

            // 9. Hyperboloid of One Sheet - doubly ruled, K<0
            { "Hyperboloid One Sheet", 
              "r(u,v) = (cosh(v)cos(u), cosh(v)sin(u), sinh(v))", 
              0.0, 2.0 * Constants::PI, 
              -1.0, 1.0, 
              HyperboloidOneSheet,
              std::numeric_limits<Real>::quiet_NaN(),
              std::numeric_limits<Real>::quiet_NaN(),  // K < 0 everywhere
              false, false },

            // 10. Paraboloid - bowl shape, K>0 everywhere
            { "Paraboloid", 
              "r(u,v) = (u, v, u² + v²)", 
              -1.0, 1.0, 
              -1.0, 1.0, 
              Paraboloid,
              std::numeric_limits<Real>::quiet_NaN(),
              std::numeric_limits<Real>::quiet_NaN(),  // K > 0 but varies
              false, false },

            // 11. Hyperbolic Paraboloid (Saddle) - doubly ruled, K<0
            { "Hyperbolic Paraboloid", 
              "r(u,v) = (u, v, u² - v²)", 
              -1.0, 1.0, 
              -1.0, 1.0, 
              HyperbolicParaboloid,
              std::numeric_limits<Real>::quiet_NaN(),
              std::numeric_limits<Real>::quiet_NaN(),  // K < 0 everywhere
              false, false },

            // 12. Klein Bottle - non-orientable, self-intersecting immersion
            { "Klein Bottle", 
              "Figure-8 immersion", 
              0.0, 2.0 * Constants::PI, 
              0.0, 2.0 * Constants::PI, 
              KleinBottle,
              std::numeric_limits<Real>::quiet_NaN(),
              std::numeric_limits<Real>::quiet_NaN(),
              false, false },

            //==========================================================================================================
            // NEW SURFACES (added for expanded test coverage)
            //==========================================================================================================

            // 13. Scherk's Surface - doubly periodic minimal surface
            { "Scherk Surface", 
              "z = ln(cos(v)/cos(u))", 
              -1.2, 1.2,                 // avoid singularities at ±π/2
              -1.2, 1.2, 
              ScherkSurface,
              std::numeric_limits<Real>::quiet_NaN(),
              std::numeric_limits<Real>::quiet_NaN(),
              true, false },             // minimal!

            // 14. Dini's Surface - twisted pseudosphere, constant K<0
            { "Dini Surface", 
              "Twisted pseudosphere with constant negative curvature", 
              0.0, 4.0 * Constants::PI, 
              0.1, Constants::PI - 0.1,  // avoid singularities at 0 and π
              DiniSurface,
              std::numeric_limits<Real>::quiet_NaN(),
              -1.0,                       // K = -1 (constant negative curvature)
              false, false },

            // 15. Boy's Surface - non-orientable immersion of RP²
            { "Boy Surface", 
              "Bryant-Kusner parameterization of RP²", 
              0.0, Constants::PI, 
              0.0, Constants::PI, 
              BoysSurface,
              std::numeric_limits<Real>::quiet_NaN(),
              std::numeric_limits<Real>::quiet_NaN(),
              false, false },

            // 16. Roman Surface (Steiner) - tetrahedral symmetry RP² immersion
            { "Roman Surface", 
              "Steiner's Roman surface - RP² with tetrahedral symmetry", 
              0.0, 2.0 * Constants::PI, 
              0.0, Constants::PI, 
              RomanSurface,
              std::numeric_limits<Real>::quiet_NaN(),
              std::numeric_limits<Real>::quiet_NaN(),
              false, false },

            // 17. Cross-Cap - simplest RP² immersion
            { "Cross Cap", 
              "Cross-cap parameterization of RP²", 
              0.0, 2.0 * Constants::PI, 
              0.0, Constants::PI, 
              CrossCap,
              std::numeric_limits<Real>::quiet_NaN(),
              std::numeric_limits<Real>::quiet_NaN(),
              false, false },

            // 18. Henneberg's Surface - classical minimal surface
            { "Henneberg Surface", 
              "Henneberg's minimal surface (1875)", 
              -1.0, 1.0, 
              0.0, 2.0 * Constants::PI, 
              HennebergSurface,
              std::numeric_limits<Real>::quiet_NaN(),
              std::numeric_limits<Real>::quiet_NaN(),
              true, false }              // minimal!
        };     
    };

    //==================================================================================================================
    // EXPLICIT SURFACE TEST BED (z = f(x,y))
    //==================================================================================================================

    class ExplicitSurfaceTestBed
    {
    public:
        static int getNumExplicitSurfaces3() { return 4; }

        const static TestParametricSurfaceRect3& getExplicitSurface3(int i) { return _listSurfaces[i]; }

        const static TestParametricSurfaceRect3& getExplicitSurface3(const std::string &surfaceName)
        {
            for (int i = 0; i < getNumExplicitSurfaces3(); i++)
            {
                if (_listSurfaces[i]._surfaceName == surfaceName)
                    return _listSurfaces[i];
            }
            throw std::runtime_error("ExplicitSurface " + surfaceName + " not found!");
        }

    private:
        const static inline TestParametricSurfaceRect3 _listSurfaces[] = { 
            // 1. Upper Unit Hemisphere: z = sqrt(1 - x² - y²)
            { "Upper Unit Hemisphere", 
              "z = sqrt(1 - x² - y²)", 
              -0.99, 0.99,   // avoid singularity at edge
              -0.99, 0.99, 
              [](Real x, Real y) { 
                  Real r2 = x*x + y*y;
                  return VectorN<Real,3>{ x, y, r2 < REAL(1.0) ? std::sqrt(REAL(1.0) - r2) : REAL(0.0) }; 
              },
              Constants::PI,              // Area of hemisphere = 2πR² = 2π for R=1 (but domain is square)
              1.0,                         // K = 1
              false, false },

            // 2. Monkey Saddle: z = x(x² - 3y²), K<0
            { "Monkey Saddle", 
              "z = x(x² - 3y²)", 
              -1.0, 1.0, 
              -1.0, 1.0, 
              [](Real x, Real y) { return VectorN<Real,3>{ x, y, x * (x*x - 3*y*y) }; },
              std::numeric_limits<Real>::quiet_NaN(),
              std::numeric_limits<Real>::quiet_NaN(),  // K < 0 near origin
              false, false },

            // 3. Egg Carton: z = sin(x)sin(y), periodic saddles
            { "Egg Carton", 
              "z = sin(x)*sin(y)", 
              -Constants::PI, Constants::PI, 
              -Constants::PI, Constants::PI, 
              [](Real x, Real y) { return VectorN<Real,3>{ x, y, std::sin(x) * std::sin(y) }; },
              std::numeric_limits<Real>::quiet_NaN(),
              std::numeric_limits<Real>::quiet_NaN(),
              false, false },

            // 4. Gaussian Bump: z = exp(-(x² + y²)), smooth with K>0 at peak
            { "Gaussian Bump", 
              "z = exp(-(x² + y²))", 
              -2.0, 2.0, 
              -2.0, 2.0, 
              [](Real x, Real y) { return VectorN<Real,3>{ x, y, std::exp(-(x*x + y*y)) }; },
              std::numeric_limits<Real>::quiet_NaN(),
              std::numeric_limits<Real>::quiet_NaN(),
              false, false }
        };     
    };    

}  // namespace MML::TestBeds

#endif  // MML_PARAMETRIC_SURFACE_TEST_BED_H
