#if !defined __MML_VECTOR_FUNCTIONS_TEST_BED_H
#define __MML_VECTOR_FUNCTIONS_TEST_BED_H

#include <string>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/Function.h"
#endif

namespace MML::TestData
{
    template<int N>
    struct TestFunctionVector
    {
        std::string _funcName;

        MML::VectorFunction<N> _func;
        MML::VectorN<Real, N> (*_funcDerived)(const MML::VectorN<Real, N> &, int ind);

        std::string _funcExpr;
        std::string _funcDerivedExpr;
        // curl

        TestFunctionVector(std::string funcName,
                            MML::VectorN<Real, N> (*f1)(const MML::VectorN<Real, N> &), std::string funcExpr, 
                            MML::VectorN<Real, N> (*f2)(const MML::VectorN<Real, N> &, int ind), std::string funcDerivedExpr
                            ) : _funcName(funcName),
                                _func(f1), _funcDerived(f2), 
                                _funcExpr(funcExpr), _funcDerivedExpr(funcDerivedExpr)
        {}
    };    

    static MML::VectorN<Real, 3> TestVectorFunc1(const MML::VectorN<Real, 3> &x) { return VectorN<Real, 3>{cos(x[0]), sin(x[1]), exp(x[2])}; }
    static MML::VectorN<Real, 3> TestVectorFunc1_derived(const MML::VectorN<Real, 3> &x, int ind) 
    { 
        if( ind == 0 ) return VectorN<Real, 3>{-sin(x[0]), 0.0, 0.0};
        else if( ind == 1 ) return VectorN<Real, 3>{0, cos(x[1]), 0};
        else return VectorN<Real, 3>{0, 0, exp(x[2])};
    }
    
    class VectorFunctionsTestBed
    {
    public:
        const static inline TestFunctionVector<3> _listFuncVector3[] = { 
            { "Simple vector func", 
              TestVectorFunc1, "( cos(x[0]) , sin(x[1]) , exp(x[2]) )", 
              TestVectorFunc1_derived, "( -sin(x[0]), 0, 0 ); ( 0, cos(x[1]), 0 ); ( 0, 0, exp(x[2] )" } 
        };
    };

}

#endif