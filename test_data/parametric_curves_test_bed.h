#if !defined MML_PARAMETRIC_CURVES_TEST_BED_H
#define MML_PARAMETRIC_CURVES_TEST_BED_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/VectorN.h"
#include "core/Function.h"
#endif

namespace MML:: TestBeds
{
    // TODO 1.0 - add couple more curves
     struct TestSpaceCurve
    {
        std::string _curveName;
        std::string _curveExpr;
        Real _start, _end;

        ParametricCurve<3> _curve;
        ParametricCurve<3> _curveDerived;
        ParametricCurve<3> _curveDerSecond;
        RealFunction _curvatureFunc;
        RealFunction _torsionFunc;

        TestSpaceCurve(std::string curveName, std::string curveExpr, Real x1, Real x2,
                        VectorN<Real,3> (*f1)(Real), VectorN<Real,3> (*f2)(Real), VectorN<Real,3> (*f3)(Real), 
                        Real (*f_curv)(Real), Real (*f_tors)(Real)  
                        ) : _start(x1), _end(x2), _curveName(curveName), _curveExpr(curveExpr), 
                            _curve(f1), _curveDerived(f2), _curveDerSecond(f3), 
                            _curvatureFunc(f_curv), _torsionFunc(f_tors)
        {}
    };   

    class ParametricCurvesTestBed
    {
    public:
        static int getNumTestCurves()            { return 3; }
        static int getNumTestCurvesArcLenParam() { return 1; }

        const static TestSpaceCurve& getTestCurve(int i)            { return _listCurves[i]; }
        const static TestSpaceCurve& getTestCurveArcLenParam(int i) { return _listCurvesArcLenParam[i]; }

        const static TestSpaceCurve& getTestCurve(const std::string &curveName)
        {
            for (int i = 0; i < getNumTestCurves(); i++)
            {
                if (_listCurves[i]._curveName == curveName)
                    return _listCurves[i];
            }
            throw std::runtime_error("TestSpaceCurve " + curveName + " not found!");
        }
        const static TestSpaceCurve& getTestCurveArcLenParam(const std::string &curveName)
        {
            for (int i = 0; i < getNumTestCurvesArcLenParam(); i++)
            {
                if (_listCurvesArcLenParam[i]._curveName == curveName)
                    return _listCurvesArcLenParam[i];
            }
            throw std::runtime_error("TestSpaceCurveArcLenParam " + curveName + " not found!");
        } 
           
    private:
        const static inline TestSpaceCurve _listCurves[] = { 
                {"Helix", "{cos(t), sin(t), t}", 
                        0.0, 2.0 * Constants::PI,  
                        [](Real t) { return VectorN<Real,3>{ cos(t), sin(t), t}; }, 
                        [](Real t) { return VectorN<Real,3>{-sin(t), cos(t), 1.0};}, 
                        [](Real t) { return VectorN<Real,3>{-cos(t), -sin(t), 0};},
                        [](Real t) { return (Real) 0.5; }, 
                        [](Real t) { return (Real) 1.0; } 
                },
                {"Schaums", "{ 3*t - t*t*t, 3*t*t, 3*t + t*t*t}", 
                        0.0, 2.0 * Constants::PI,  
                        [](Real t) { return VectorN<Real,3>{ 3*t - t*t*t, 3*t*t, 3*t + t*t*t}; }, 
                        [](Real t) { return VectorN<Real,3>{-sin(t), cos(t), 1.0};}, 
                        [](Real t) { return VectorN<Real,3>{-cos(t), -sin(t), 0};},
                        [](Real t) { return (Real) 1.0 / (3 * POW2(1 + t*t)); }, 
                        [](Real t) { return (Real) 2.0 / (3 * POW2(1 + t*t)); } 
                },
                {"Schaums2", "{ t - sin(t), 1 - cos(t), t}", 
                        0.0, 2.0 * Constants::PI,  
                        [](Real t) { return VectorN<Real,3>{ t - sin(t), 1 - cos(t), t}; }, 
                        [](Real t) { return VectorN<Real,3>{0, 0, 0};}, 
                        [](Real t) { return VectorN<Real,3>{0, 0, 0};},
                        [](Real t) { return sqrt(1 + 4 * (Real) pow(sin(t/2), 4)) / (Real) pow((1 + 4 * pow(sin(t/2), 2)), 1.5); }, 
                        [](Real t) { return -1 / (1 + 4 * (Real) pow(sin(t/2), 4)); } 
                }                        
            };        

        const static inline TestSpaceCurve _listCurvesArcLenParam[] = { 
                {"Helix", "1/sqrt(2) * {cos(t), sin(t), t}", 
                        0.0, 2.0 * Constants::PI,  
                        [](Real t) { return VectorN<Real,3>{ cos(t) / sqrt((Real) 2.0), sin(t) / sqrt((Real) 2.0), t / sqrt((Real) 2.0)}; }, 
                        [](Real t) { return VectorN<Real,3>{-sin(t) / sqrt((Real) 2.0), cos(t) / sqrt((Real) 2.0), 1 / sqrt((Real) 2.0)};}, 
                        [](Real t) { return VectorN<Real,3>{-cos(t) / sqrt((Real) 2.0), -sin(t) / sqrt((Real) 2.0), 0};},
                        [](Real t) { return Real{0.5}; }, 
                        [](Real t) { return Real{1.0}; } 
                }
            };               
    };
}
#endif