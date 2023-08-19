#if !defined MML_PARAMETRIC_CURVES_TEST_BED_H
#define MML_PARAMETRIC_CURVES_TEST_BED_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/VectorN.h"
#include "basic_types/Function.h"
#endif

namespace MML::TestData
{
    // TODO - add planar curves
    // TODO - add planar polar curves
    // TODO - add 5 curves with arc length parametrization
    // arbitrary parametrization
    struct TestSpaceCurve
    {
        std::string _curveName;
        std::string _curveExpr;
        double _start, _end;

        MML::ParametricCurve<3> _curve;
        MML::ParametricCurve<3> _curveDerived;
        MML::ParametricCurve<3> _curveDerSecond;
        MML::RealFunction _curvatureFunc;
        MML::RealFunction _torsionFunc;

        TestSpaceCurve(std::string curveName, std::string curveExpr, double x1, double x2,
                        VectorN<Real,3> (*f1)(double), VectorN<Real,3> (*f2)(double), VectorN<Real,3> (*f3)(double), 
                        double (*f_curv)(double), double (*f_tors)(double)  
                        ) : _start(x1), _end(x2), _curveName(curveName), _curveExpr(curveExpr), 
                            _curve(f1), _curveDerived(f2), _curveDerSecond(f3), 
                            _curvatureFunc(f_curv), _torsionFunc(f_tors)
        {}
    };   

    struct TestCurveArcLengthParametrization
    {

    };

    class ParametricCurvesTestBed
    {
    public:
        const static inline TestSpaceCurve _listCurves[] = { 
                {"Helix", "{cos(t), sin(t), t}", 0.0, 2.0 * Constants::PI,  
                        [](double t) { return VectorN<Real,3>{ cos(t), sin(t), t}; }, 
                        [](double t) { return VectorN<Real,3>{-sin(t), cos(t), 1.0};}, 
                        [](double t) { return VectorN<Real,3>{-cos(t), -sin(t), 0};},
                        [](double t) { return 0.5; }, 
                        [](double t) { return 1.0; } 
                },
                {"Schaums", "{ 3*t - t*t*t, 3*t*t, 3*t + t*t*t}", 0.0, 2.0 * Constants::PI,  
                        [](double t) { return VectorN<Real,3>{ 3*t - t*t*t, 3*t*t, 3*t + t*t*t}; }, 
                        [](double t) { return VectorN<Real,3>{-sin(t), cos(t), 1.0};}, 
                        [](double t) { return VectorN<Real,3>{-cos(t), -sin(t), 0};},
                        [](double t) { return 1.0 / (3 * SQR(1 + t*t)); }, 
                        [](double t) { return 2.0 / (3 * SQR(1 + t*t)); } 
                },
                {"Schaums2 ", "{ t - sin(t), 1 - cos(t), t}", 0.0, 2.0 * Constants::PI,  
                        [](double t) { return VectorN<Real,3>{ t - sin(t), 1 - cos(t), t}; }, 
                        [](double t) { return VectorN<Real,3>{0, 0, 0};}, 
                        [](double t) { return VectorN<Real,3>{0, 0, 0};},
                        [](double t) { return sqrt(1 + 4 * pow(sin(t/2), 4)) / pow((1 + 4 * pow(sin(t/2), 2)), 1.5); }, 
                        [](double t) { return -1.0 / (1 + 4 * pow(sin(t/2), 4)); } 
                }                        
            };        

        const static inline TestSpaceCurve _listCurvesArcLenParam[] = { 
                {"Helix", "1/sqrt(2) * {cos(t), sin(t), t}", 0.0, 2.0 * Constants::PI,  
                        [](double t) { return VectorN<Real,3>{ cos(t) / sqrt(2.0), sin(t) / sqrt(2.0), t / sqrt(2.0)}; }, 
                        [](double t) { return VectorN<Real,3>{-sin(t), cos(t), 1.0};}, 
                        [](double t) { return VectorN<Real,3>{-cos(t), -sin(t), 0};},
                        [](double t) { return 0.5; }, 
                        [](double t) { return 1.0; } 
                }
            };               
    };
}
#endif