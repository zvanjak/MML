#if !defined __MML_SCALAR_FUNCTIONS_TEST_BED_H
#define __MML_SCALAR_FUNCTIONS_TEST_BED_H

#include <string>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Function.h"
#endif

namespace MML::TestBeds
{
    template<int N>
    struct TestFunctionScalar
    {
        std::string _funcName;

        ScalarFunction<N> _func;
        Real (*_funcDerived)(const VectorN<Real, N> &, int ind);

        std::string _funcExpr;
        std::string _funcDerivedExpr;
        // gradijent

        TestFunctionScalar(std::string funcName,
                            Real (*f1)(const VectorN<Real, N> &), std::string funcExpr, 
                            Real (*f2)(const VectorN<Real, N> &, int ind), std::string funcDerivedExpr
                            ) : _funcName(funcName),
                                _func(f1), _funcDerived(f2), 
                                _funcExpr(funcExpr), _funcDerivedExpr(funcDerivedExpr)
        {}
    };    

    // TODO - nekoliko primjera različitih realnih skalarnih polja - potencijali
    
    static Real TestScalarFunc1(const VectorN<Real, 3> &x) { return cos(x[0]) + sin(x[1]) + exp(x[2]); }
    static Real TestScalarFunc1_derived(const VectorN<Real, 3> &x, int ind) 
    { 
        if( ind == 0 ) return -sin(x[0]);
        else if( ind == 1 ) return cos(x[1]);
        else return exp(x[2]);
    }
    static Real TestScalarFunc2(const VectorN<Real, 3> &x) { return sin(x[0] * x[1]) * exp(x[2] / (x[1] * x[1] +1 )) / (1 + x[0] * x[0]); }
    static Real TestScalarFunc2_derived(const VectorN<Real, 3> &xVal, int ind) 
    { 
        Real x = xVal[0];
        Real y = xVal[1];
        Real z = xVal[2];
        if( ind == 0 ) return (exp(z / (y*y + 1)) * ((x*x + 1) * y * cos(x * y) - 2 * x * sin(x * y)))/pow((x*x + 1),2);
        else if( ind == 1 ) return (exp(z/(y*y + 1)) * (x * pow((y*y + 1),2) * cos(x * y) - 2 * y * z * sin(x * y)))/((x*x + 1) * pow((y*y + 1),2));
        else return (exp(z/(y*y + 1)) * sin(x * y))/((x*x + 1) * (y*y + 1));
    }    
    class ScalarFunctionsTestBed
    {
    public:
        static int getNumTestFunctionScalar3() { return 2; }

        const static TestFunctionScalar<3>& getTestFunctionScalar3(int i)  { return _listFuncScalar3[i]; }

        const static TestFunctionScalar<3>& getTestFunctionScalar3(const std::string &funcName)
        {
            for (int i = 0; i < getNumTestFunctionScalar3(); i++)
            {
                if (_listFuncScalar3[i]._funcName == funcName)
                    return _listFuncScalar3[i];
            }
            throw std::runtime_error("TestFunctionScalar " + funcName + " not found!");
        }
    
    private:
        const static inline TestFunctionScalar<3> _listFuncScalar3[] = { 
            { "Scalar func 1", 
              TestScalarFunc1, "cos(x[0]) + sin(x[1]) + exp(x[2])", 
              TestScalarFunc1_derived, "-sin(x[0]); cos(x[1]); exp(x[2]" },
            { "Scalar func 2", 
              TestScalarFunc2, "sin(x*y)*exp(z/(y*y+1))/(1+x*x)", 
            // d/dx((sin(x y) exp(z/(y y + 1)))/(1 + x x)) = (e^(z/(y^2 + 1)) ((x^2 + 1) y cos(x y) - 2 x sin(x y)))/(x^2 + 1)^2
            // d/dy((sin(x y) exp(z/(y y + 1)))/(1 + x x)) = (e^(z/(y^2 + 1)) (x (y^2 + 1)^2 cos(x y) - 2 y z sin(x y)))/((x^2 + 1) (y^2 + 1)^2)
            // d/dz((sin(x y) exp(z/(y y + 1)))/(1 + x x)) = (e^(z/(y^2 + 1)) sin(x y))/((x^2 + 1) (y^2 + 1))
              TestScalarFunc2_derived, "(e^(z/(y^2 + 1)) ((x^2 + 1) y cos(x y) - 2 x sin(x y)))/(x^2 + 1)^2; (e^(z/(y^2 + 1)) (x (y^2 + 1)^2 cos(x y) - 2 y z sin(x y)))/((x^2 + 1) (y^2 + 1)^2); (e^(z/(y^2 + 1)) sin(x y))/((x^2 + 1) (y^2 + 1))" }               
        };
        
    };
}

#endif