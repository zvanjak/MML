#if !defined __MML_VECTOR_FUNCTIONS_TEST_BED_H
#define __MML_VECTOR_FUNCTIONS_TEST_BED_H

#include <string>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Function.h"
#endif

namespace MML::TestBeds
{
    template<int N>
    struct TestFunctionVector
    {
        std::string _funcName;

        VectorFunction<N> _func;
        VectorN<Real, N> (*_funcDerived)(const MML::VectorN<Real, N> &, int ind);

        std::string _funcExpr;
        std::string _funcDerivedExpr;
        // curl

        TestFunctionVector( std::string funcName,
                            VectorN<Real, N> (*f1)(const VectorN<Real, N> &), std::string funcExpr, 
                            VectorN<Real, N> (*f2)(const VectorN<Real, N> &, int ind), std::string funcDerivedExpr
                            ) : _funcName(funcName),
                                _func(f1), _funcDerived(f2), 
                                _funcExpr(funcExpr), _funcDerivedExpr(funcDerivedExpr)
        {}
    };    

    static MML::VectorN<Real, 3> TestVectorFunc1(const VectorN<Real, 3> &xVal) 
    {
        Real x = xVal[0];
        Real y = xVal[1];
        Real z = xVal[2];

        Real valx = x*cos(y)*z*z;
        Real valy = sin(x)*(y*y + z*z);
        Real valz = exp(x*y/(z*z+1));

        return VectorN<Real, 3>{valx, valy, valz};
    }
    static MML::VectorN<Real, 3> TestVectorFunc1_derived(const VectorN<Real, 3> &xVal, int ind) 
    { 
        Real x = xVal[0];
        Real y = xVal[1];
        Real z = xVal[2];
// d/dx(x cos(y) z z) = z^2 cos(y)
// d/dy(x cos(y) z z) = -x z^2 sin(y)
// d/dz(x cos(y) z z) = 2 x z cos(y)

// d/dx(sin(x) (y y + z z)) = cos(x) (y^2 + z^2)
// d/dy(sin(x) (y y + z z)) = 2 y sin(x)
// d/dz(sin(x) (y y + z z)) = 2 z sin(x)

// d/dx(exp((x y)/(z z + 1))) = (y e^((x y)/(z^2 + 1)))/(z^2 + 1)
// d/dy(exp((x y)/(z z + 1))) = (x e^((x y)/(z^2 + 1)))/(z^2 + 1)
// d/dz(exp((x y)/(z z + 1))) = -(2 x y z e^((x y)/(z^2 + 1)))/(z^2 + 1)^2

        if( ind == 0 ) 
            return VectorN<Real, 3>{z*z*cos(y)                            , -x * z*z * sin(y)                   , 2 * x * z * cos(y)};
        else if( ind == 1 ) 
            return VectorN<Real, 3>{cos(x) * (y*y + z*z)                  , 2 * y * sin(x)                      , 2 * z * sin(x)};
        else 
            return VectorN<Real, 3>{(y * exp((x * y)/(z*z + 1)))/(z*z + 1), x * exp((x * y)/(z*z + 1))/(z*z + 1), -(2 * x * y * z * exp((x * y)/(z*z + 1))) / (Real) pow((z*z + 1),2)};
    }
    
    class VectorFunctionsTestBed
    {
    public:
        static int getNumTestFunctionVector() { return 1; }

        const static TestFunctionVector<3>& getTestFunctionVector(int i)  { return _listFuncVector3[i]; }

        const static TestFunctionVector<3>& getTestFunctionVector(const std::string &funcName)
        {
            for (int i = 0; i < getNumTestFunctionVector(); i++)
            {
                if (_listFuncVector3[i]._funcName == funcName)
                    return _listFuncVector3[i];
            }
            throw std::runtime_error("TestFunctionVector " + funcName + " not found!");
        }
    private:
        const static inline TestFunctionVector<3> _listFuncVector3[] = { 
            { "Simple vector func", 
              TestVectorFunc1, "( x*cos(y)*z*z , sin(x)*(y*y + z*z) , exp(x*y/(z*z+1)) )", 
              TestVectorFunc1_derived, "( -sin(x[0]), 0, 0 ); ( 0, cos(x[1]), 0 ); ( 0, 0, exp(x[2] )" } 
        };
    };

}

#endif