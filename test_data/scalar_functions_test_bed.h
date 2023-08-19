#if !defined __MML_SCALAR_FUNCTIONS_TEST_BED_H
#define __MML_SCALAR_FUNCTIONS_TEST_BED_H

#include <string>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/Function.h"
#endif

namespace MML::TestData
{
    template<int N>
    struct TestFunctionScalar
    {
        std::string _funcName;

        MML::ScalarFunction<N> _func;
        double (*_funcDerived)(const MML::VectorN<Real, N> &, int ind);

        std::string _funcExpr;
        std::string _funcDerivedExpr;
        // gradijent

        TestFunctionScalar(std::string funcName,
                            double (*f1)(const MML::VectorN<Real, N> &), std::string funcExpr, 
                            double (*f2)(const MML::VectorN<Real, N> &, int ind), std::string funcDerivedExpr
                            ) : _funcName(funcName),
                                _func(f1), _funcDerived(f2), 
                                _funcExpr(funcExpr), _funcDerivedExpr(funcDerivedExpr)
        {}
    };    

    // primjeri različitih realnih skalarnih polja - potencijali
    
    static double TestScalarFunc1(const MML::VectorN<Real, 3> &x) { return cos(x[0]) + sin(x[1]) + exp(x[2]); }
    static double TestScalarFunc1_derived(const MML::VectorN<Real, 3> &x, int ind) 
    { 
        if( ind == 0 ) return -sin(x[0]);
        else if( ind == 1 ) return cos(x[1]);
        else return exp(x[2]);
    }
    class ScalarFunctionsTestBed
    {
    public:
        const static inline TestFunctionScalar<3> _listFuncScalar3[] = { 
            { "Scalar func 1", 
              TestScalarFunc1, "cos(x[0]) + sin(x[1]) + exp(x[2])", 
              TestScalarFunc1_derived, "-sin(x[0]); cos(x[1]); exp(x[2]" } 
        };
    };

}

#endif