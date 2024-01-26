#if !defined MML_PARAMETRIC_SURFACE_TEST_BED_H
#define MML_PARAMETRIC_SURFACE_TEST_BED_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/VectorN.h"
#include "base/Matrix.h"

#include "core/Function.h"
#endif

using namespace MML;

namespace MML::TestBeds
{
    struct TestParametricSurface3
    {
        std::string _surfaceName;
        std::string _surfaceExpr;
        double _x1, _x2;
        double _y1, _y2;

        MML::ParametricSurface<3> _surface;

        TestParametricSurface3( std::string surfaceName, std::string surfaceExpr, double x1, double x2, double y1, double y2,
                                VectorN<Real,3> (*f1)(double, double)
                              ) : _x1(x1), _x2(x2), _surfaceName(surfaceName), _surfaceExpr(surfaceExpr), 
                                  _surface(f1)
        {}
    };      

    class ParametricSurfaceTestBed
    {
    public:
        static int getNumTestSurfaces3() { return 1; }

        const static TestParametricSurface3& getTestSurface3(int i)            { return _listSurfaces[i]; }

        const static TestParametricSurface3& getTestSurface3(const std::string &surfaceName)
        {
            for (int i = 0; i < getNumTestSurfaces3(); i++)
            {
                if (_listSurfaces[i]._surfaceName == surfaceName)
                    return _listSurfaces[i];
            }
            throw std::runtime_error("ParametricSurface " + surfaceName + " not found!");
        }

    private:
        const static inline TestParametricSurface3 _listSurfaces[] = { 
                {"Test", "{cos(t), sin(t), t}", 0.0, 2.0 * Constants::PI, 0.0, 2.0 * Constants::PI, 
                        [](double u, double v) { return VectorN<Real,3>{ cos(u), sin(u), u}; } 
                }
        };     
    };

    class ExplicitSurfaceTestBed
    {
    public:
        static int getNumExplicitSurfaces3() { return 2; }

        const static TestParametricSurface3& getExplicitSurface3(int i)            { return _listSurfaces[i]; }

        const static TestParametricSurface3& getExplicitSurface3(const std::string &surfaceName)
        {
            for (int i = 0; i < getNumExplicitSurfaces3(); i++)
            {
                if (_listSurfaces[i]._surfaceName == surfaceName)
                    return _listSurfaces[i];
            }
            throw std::runtime_error("ExplicitSurface " + surfaceName + " not found!");
        }

    private:
        const static inline TestParametricSurface3 _listSurfaces[] = { 
                // TODO 0.8 - treba modelirati preciznije domenu u x-y ravnini, za sferu
                {"UpperUnitSphere", "todo", -1.0, 1.0, -1.0, 1.0, 
                                    [](double x, double y) { return VectorN<Real,3>{ x, y, x*x+y*y>1.0 ? 0.0 : sqrt(1 - x*x - y*y)}; } 
                },
                {"Monkey saddle", "todo", -1.0, 1.0, -1.0, 1.0, 
                                    [](double x, double y) { return VectorN<Real,3>{ x, y, x * (x*x - 3 * y*y)}; } 
                }
        };     
    };    
}
#endif