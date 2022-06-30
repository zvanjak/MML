#if !defined __MML_FUNCTIONS_TEST_BED_H
#define __MML_FUNCTIONS_TEST_BED_H

#include <string>

#ifdef MML_USE_SINGLE_HEADER
#include "MMLBasicTypes.h"
#else
#include "basic_types/Function.h"
#endif

namespace MML::Tests
{
    static double TestFunc1_Sin(double x) { return sin(x); }
    static double TestFunc1_Sin_derived(double x) { return cos(x); }
    static double TestFunc1_Sin_integrated(double x) { return -cos(x); }

    struct TestFunction
    {
        double _start, _end;
        std::string _funcName;

        MML::RealFunction _func;
        MML::RealFunction _funcDerived;
        MML::RealFunction _funcIntegrated;

        std::string _funcExpr;
        std::string _funcDerivedExpr;
        std::string _funcIntegratedExpr;

        TestFunction(double x1, double x2, std::string funcName,
                     double (*f1)(double), double (*f2)(double), double (*f3)(double) ,
                     std::string funcExpr, std::string funcDerivedExpr, std::string funcIntegratedExpr
                     ) : _start(x1), _end(x2), _funcName(funcName),
                         _func(f1), _funcDerived(f2), _funcIntegrated(f3),
                         _funcExpr(funcExpr), _funcDerivedExpr(funcDerivedExpr), _funcIntegratedExpr(funcIntegratedExpr)
        {}
    };

    class FunctionsTestBed
    {
    public:
        //const static inline TestFunction _f1{ 0.0, 1, Func1, Func1_derived, Func1_integrated };

        const static inline TestFunction _listFunc[] = { {0.0, 1, "Sin(x)", TestFunc1_Sin, TestFunc1_Sin_derived, TestFunc1_Sin_integrated, "sin(x)", "cos(x)", "-cos(x)"} };
    };

}

#endif