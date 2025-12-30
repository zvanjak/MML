#if !defined MML_PARAMETRIC_CURVES_TEST_BED_H
#define MML_PARAMETRIC_CURVES_TEST_BED_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/VectorN.h"
#include "base/Function.h"
#include "core/Curves.h"
#endif

#include "parametric_curves_defs.h"

using namespace MML::Curves;

namespace MML::TestBeds
{
    /*******************************************************************************************************************
     * PARAMETRIC CURVES TEST BED
     * 
     * Provides test infrastructure for parametric curve analysis:
     *   - Test data structures (TestCartesianCurve3D)
     *   - Collection of test curves with known properties
     *   - Arc-length parameterized variants
     *   - Access by index or name
     * 
     * Uses curve definitions from parametric_curves_defs.h
     *******************************************************************************************************************/

    //==================================================================================================================
    // TEST DATA STRUCTURES
    //==================================================================================================================

    /**
     * @brief Test data for a 3D Cartesian parametric curve
     * 
     * Contains the curve, its derivatives, and analytical curvature/torsion functions.
     */
    struct TestCartesianCurve3D
    {
        std::string _curveName;
        std::string _curveExpr;
        Real _start, _end;

        CurveCartesian3D _curve;
        CurveCartesian3D _curveDerived;
        CurveCartesian3D _curveDerSecond;
        RealFunction _curvatureFunc;
        RealFunction _torsionFunc;

        TestCartesianCurve3D( std::string curveName, std::string curveExpr, Real x1, Real x2,
                              VectorN<Real, 3>(*f1)(Real), VectorN<Real, 3>(*f2)(Real), VectorN<Real, 3>(*f3)(Real),
                              Real(*f_curv)(Real), Real(*f_tors)(Real) ) 
            : _start(x1), _end(x2), _curveName(curveName), _curveExpr(curveExpr),
              _curve(f1), _curveDerived(f2), _curveDerSecond(f3),
              _curvatureFunc(f_curv), _torsionFunc(f_tors)
        {
        }
    };

    //==================================================================================================================
    // PARAMETRIC CURVES TEST BED
    //==================================================================================================================

    class ParametricCurvesTestBed
    {
    public:
        static int getNumTestCurves() { return 13; }  // Total number of test curves
        static int getNumTestCurvesArcLenParam() { return 1; }

        const static TestCartesianCurve3D& getTestCurve(int i) { return _listCurves[i]; }
        const static TestCartesianCurve3D& getTestCurveArcLenParam(int i) { return _listCurvesArcLenParam[i]; }

        const static TestCartesianCurve3D& getTestCurve(const std::string& curveName)
        {
            for (int i = 0; i < getNumTestCurves(); i++)
            {
                if (_listCurves[i]._curveName == curveName)
                    return _listCurves[i];
            }
            throw std::runtime_error("TestCartesianCurve3D " + curveName + " not found!");
        }

        const static TestCartesianCurve3D& getTestCurveArcLenParam(const std::string& curveName)
        {
            for (int i = 0; i < getNumTestCurvesArcLenParam(); i++)
            {
                if (_listCurvesArcLenParam[i]._curveName == curveName)
                    return _listCurvesArcLenParam[i];
            }
            throw std::runtime_error("TestSpaceCurveArcLenParam " + curveName + " not found!");
        }

    private:
        //==============================================================================================================
        // MAIN TEST CURVE COLLECTION
        //==============================================================================================================
        const static inline TestCartesianCurve3D _listCurves[] = {
            // --- Basic curves (using simple lambdas) ---
            {"Helix", "{cos(t), sin(t), t}",
                Real(0.0), Real(2.0) * Constants::PI,
                [](Real t) { return VectorN<Real,3>{ cos(t), sin(t), t}; },
                [](Real t) { return VectorN<Real,3>{-sin(t), cos(t), 1.0}; },
                [](Real t) { return VectorN<Real,3>{-cos(t),-sin(t), 0}; },
                [](Real t) { return (Real)0.5; },
                [](Real t) { return (Real)0.5; }
            },
            {"Helix2", "{2*cos(t), 2*sin(t), 0.5*t}",
                Real(0.0), Real(2.0) * Constants::PI,
                [](Real t) -> VectorN<Real,3> { return VectorN<Real,3>{ Real(2)*cos(t), Real(2)*sin(t), Real(0.5)*t}; },
                [](Real t) -> VectorN<Real,3> { return VectorN<Real,3>{Real(-2)*sin(t), Real(2)*cos(t), Real(0.5)}; },
                [](Real t) -> VectorN<Real,3> { return VectorN<Real,3>{Real(-2)*cos(t),Real(-2)*sin(t), Real(0)}; },
                [](Real t) -> Real { return Real(2.0)/(Real(4.0) + Real(0.25)); },
                [](Real t) -> Real { return Real(0.5)/(Real(4.0) + Real(0.25)); }
            },
            {"Schaums1", "{ 3*t - t*t*t, 3*t*t, 3*t + t*t*t}",
                Real(0.0), Real(2.0) * Constants::PI,
                [](Real t) { return VectorN<Real,3>{ 3 * t - t * t * t, 3 * t * t, 3 * t + t * t * t}; },
                [](Real t) { return VectorN<Real,3>{-sin(t), cos(t), 1.0}; },
                [](Real t) { return VectorN<Real,3>{-cos(t), -sin(t), 0}; },
                [](Real t) { return (Real)1.0 / (3 * POW2(1 + t * t)); },
                [](Real t) { return (Real)1.0 / (3 * POW2(1 + t * t)); }
            },
            {"Schaums2", "{ t - sin(t), 1 - cos(t), t}",
                Real(0.0), Real(2.0) * Constants::PI,
                [](Real t) { return VectorN<Real,3>{ t - sin(t), 1 - cos(t), t}; },
                [](Real t) { return VectorN<Real,3>{0, 0, 0}; },
                [](Real t) { return VectorN<Real,3>{0, 0, 0}; },
                [](Real t) { return sqrt(1 + 4 * (Real)pow(sin(t / 2), 4)) / (Real)pow((1 + 4 * pow(sin(t / 2), 2)), 1.5); },
                [](Real t) { return -1 / (1 + 4 * (Real)pow(sin(t / 2), 4)); }
            },
            {"Circle3DXY", "{cos(t), sin(t), 0}",
                Real(0.0), Real(2.0) * Constants::PI,
                [](Real t) { return VectorN<Real,3>{ cos(t), sin(t), 0.0}; },
                [](Real t) { return VectorN<Real,3>{-sin(t), cos(t), 0.0}; },
                [](Real t) { return VectorN<Real,3>{-cos(t), -sin(t), 0.0}; },
                [](Real t) { return Real{1.0}; },  // curvature = 1 for unit circle
                [](Real t) { return Real{0.0}; }   // torsion = 0 for planar curve
            },

            // --- Polynomial curves (static functions) ---
            {"TwistedCubic", "{t, t², t³}",
                Real(-2.0), Real(2.0),
                TwistedCubic_pos, TwistedCubic_deriv, TwistedCubic_deriv2,
                TwistedCubic_curvature, TwistedCubic_torsion
            },

            // --- Trigonometric curves (static functions) ---
            {"Viviani", "{1 + cos(t), sin(t), 2*sin(t/2)}",
                Real(0.0), Real(2.0) * Constants::PI,
                Viviani_pos, Viviani_deriv, Viviani_deriv2,
                Viviani_curvature, Viviani_torsion
            },
            {"Lissajous_2_3_5", "{sin(2t), sin(3t), sin(5t)}",
                Real(0.0), Real(2.0) * Constants::PI,
                Lissajous235_pos, Lissajous235_deriv, Lissajous235_deriv2,
                Lissajous235_curvature, Lissajous235_torsion
            },

            // --- Knots (static functions) ---
            {"TorusKnot_3_2", "{ (2 + cos(2t)) cos(3t), (2 + cos(2t)) sin(3t), sin(2t) }",
                Real(0.0), Real(2.0) * Constants::PI,
                TorusKnot32_pos, TorusKnot32_deriv, TorusKnot32_deriv2,
                TorusKnot32_curvature, TorusKnot32_torsion
            },
            {"TrefoilKnot", "{ sin(t) + 2sin(2t), cos(t) - 2cos(2t), -sin(3t) }",
                Real(0.0), Real(2.0) * Constants::PI,
                TrefoilKnot_pos, TrefoilKnot_deriv, TrefoilKnot_deriv2,
                TrefoilKnot_curvature, TrefoilKnot_torsion
            },
            {"FigureEightKnot", "{ (2 + cos(2t)) cos(t), (2 + cos(2t)) sin(t), sin(3t) }",
                Real(0.0), Real(2.0) * Constants::PI,
                FigureEightKnot_pos, FigureEightKnot_deriv, FigureEightKnot_deriv2,
                FigureEightKnot_curvature, FigureEightKnot_torsion
            },

            // --- Special curves (static functions) ---
            {"ConicalHelix", "{ t*cos(t), t*sin(t), t }",
                Real(0.1), Real(4.0) * Constants::PI,  // Start from 0.1 to avoid singularity
                ConicalHelix_pos, ConicalHelix_deriv, ConicalHelix_deriv2,
                ConicalHelix_curvature, ConicalHelix_torsion
            },
            {"SphericalSpiral", "{ cos(t)/cosh(0.3t), sin(t)/cosh(0.3t), tanh(0.3t) }",
                Real(-10.0), Real(10.0),
                SphericalSpiral_pos, SphericalSpiral_deriv, SphericalSpiral_deriv2,
                SphericalSpiral_curvature, SphericalSpiral_torsion
            }
        };

        //==============================================================================================================
        // ARC-LENGTH PARAMETERIZED CURVES
        //==============================================================================================================
        const static inline TestCartesianCurve3D _listCurvesArcLenParam[] = {
            {"Helix", "1/sqrt(2) * {cos(t), sin(t), t}",
                Real(0.0), Real(2.0) * Constants::PI,
                [](Real t) { return VectorN<Real,3>{ cos(t) / sqrt((Real)2.0), sin(t) / sqrt((Real)2.0), t / sqrt((Real)2.0)}; },
                [](Real t) { return VectorN<Real,3>{-sin(t) / sqrt((Real)2.0), cos(t) / sqrt((Real)2.0), 1 / sqrt((Real)2.0)}; },
                [](Real t) { return VectorN<Real,3>{-cos(t) / sqrt((Real)2.0), -sin(t) / sqrt((Real)2.0), 0}; },
                [](Real t) { return Real{0.5}; },
                [](Real t) { return Real{1.0}; }
            }
        };
    };

}  // namespace MML::TestBeds

#endif  // MML_PARAMETRIC_CURVES_TEST_BED_H
