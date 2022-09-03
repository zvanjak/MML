#if !defined __MML_FUNCTIONS_TEST_BED_H
#define __MML_FUNCTIONS_TEST_BED_H

#include <string>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/Function.h"
#endif

namespace MML::Tests
{
    static double TestFunc1_Sin(double x) { return sin(x); }
    static double TestFunc1_Sin_derived(double x) { return cos(x); }
    static double TestFunc1_Sin_integrated(double x) { return -cos(x); }

    static double TestFunc2_Cos(double x) { return cos(x); }
    static double TestFunc2_Cos_derived(double x) { return -sin(x); }
    static double TestFunc2_Cos_integrated(double x) { return sin(x); }

    static double TestScalarFunc1(const MML::VectorN<3> &x) { return cos(x[0]) + sin(x[1]) + exp(x[2]); }
    static double TestScalarFunc1_derived(const MML::VectorN<3> &x, int ind) 
    { 
        if( ind == 0 ) return -sin(x[0]);
        else if( ind == 1 ) return cos(x[1]);
        else return exp(x[2]);
    }

    static MML::VectorN<3> TestVectorFunc1(const MML::VectorN<3> &x) { return VectorN<3>{cos(x[0]), sin(x[1]), exp(x[2])}; }
    static MML::VectorN<3> TestVectorFunc1_derived(const MML::VectorN<3> &x, int ind) 
    { 
        if( ind == 0 ) return VectorN<3>{-sin(x[0]), 0.0, 0.0};
        else if( ind == 1 ) return VectorN<3>{0, cos(x[1]), 0};
        else return VectorN<3>{0, 0, exp(x[2])};
    }

    struct TestFunctionReal
    {
        double _start, _end;
        std::string _funcName;

        MML::RealFunction _func;
        MML::RealFunction _funcDerived;
        MML::RealFunction _funcIntegrated;

        std::string _funcExpr;
        std::string _funcDerivedExpr;
        std::string _funcIntegratedExpr;

        TestFunctionReal(double x1, double x2, std::string funcName,
                        double (*f1)(double), std::string funcExpr, double (*f2)(double), std::string funcDerivedExpr, double (*f3)(double), std::string funcIntegratedExpr
                        ) : _start(x1), _end(x2), _funcName(funcName),
                            _func(f1), _funcDerived(f2), _funcIntegrated(f3),
                            _funcExpr(funcExpr), _funcDerivedExpr(funcDerivedExpr), _funcIntegratedExpr(funcIntegratedExpr)
        {}
    };

    struct TestFunctionScalar3
    {
        std::string _funcName;

        MML::ScalarFunctionFromFuncPtr<3> _func;
        double (*_funcDerived)(const MML::VectorN<3> &, int ind);

        std::string _funcExpr;
        std::string _funcDerivedExpr;

        TestFunctionScalar3(std::string funcName,
                            double (*f1)(const MML::VectorN<3> &), double (*f2)(const MML::VectorN<3> &, int ind), 
                            std::string funcExpr, std::string funcDerivedExpr
                            ) : _funcName(funcName),
                                _func(f1), _funcDerived(f2), 
                                _funcExpr(funcExpr), _funcDerivedExpr(funcDerivedExpr)
        {}
    };    

    struct TestFunctionVector3
    {
        std::string _funcName;

        MML::VectorFunctionFromFuncPtr<3> _func;
        MML::VectorN<3> (*_funcDerived)(const MML::VectorN<3> &, int ind);

        std::string _funcExpr;
        std::string _funcDerivedExpr;

        TestFunctionVector3(std::string funcName,
                            MML::VectorN<3> (*f1)(const MML::VectorN<3> &), MML::VectorN<3> (*f2)(const MML::VectorN<3> &, int ind), 
                            std::string funcExpr, std::string funcDerivedExpr
                            ) : _funcName(funcName),
                                _func(f1), _funcDerived(f2), 
                                _funcExpr(funcExpr), _funcDerivedExpr(funcDerivedExpr)
        {}
    };    

    class FunctionsTestBed
    {
    public:
        const static inline TestFunctionReal _listFuncReal[] = { 
                {0.0, 1, "Sin(x)", [](double x) { return sin(x);},  "sin(x)", 
                                   [](double x) { return cos(x);},  "cos(x)",  
                                   [](double x) { return -cos(x);}, "-cos(x)"},
                {0.0, 1, "Cos(x)", [](double x) { return cos(x);}, "cos(x)", 
                                   [](double x) { return -sin(x);}, "-sin(x)", 
                                   [](double x) { return sin(x);}, "sin(x)"},
                {0.0, 1, "x^3",    [](double x) { return x*x*x;}, "x^3",    
                                   [](double x) { return 3*x*x;}, "3*x^2",
                                   [](double x) { return 1.0/4.0*x*x*x*x;}, "1/4 * x^4"},
                {0.0, 1, "Exp(x)", [](double x) { return exp(x);}, "exp(x)", 
                                   [](double x) { return exp(x);}, "exp(x)",
                                   [](double x) { return exp(x);}, "exp(x)"}
            };

        const static inline TestFunctionScalar3 _listFuncScalar[] = { {"Simple func", TestScalarFunc1, TestScalarFunc1_derived, "cos(x[0]) + sin(x[1]) + exp(x[2])", "-sin(x[0]); cos(x[1]); exp(x[2]"} };
        const static inline TestFunctionVector3 _listFuncVector[] = { {"Simple vector func", TestVectorFunc1, TestVectorFunc1_derived, "( cos(x[0]) , sin(x[1]) , exp(x[2]) )", "( -sin(x[0]), 0, 0 ); ( 0, cos(x[1]), 0 ); ( 0, 0, exp(x[2] )"} };

    };

}

#endif