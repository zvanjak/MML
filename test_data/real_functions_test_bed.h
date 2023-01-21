#if !defined __MML_REAL_FUNCTIONS_TEST_BED_H
#define __MML_REAL_FUNCTIONS_TEST_BED_H

#include <string>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/Function.h"
#endif

namespace MML::TestData
{
    // dodati test bed za standard functions
    
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
                        double (*f1)(double), std::string funcExpr, 
                        double (*f2)(double), std::string funcDerivedExpr, 
                        double (*f3)(double), std::string funcIntegratedExpr
                        ) : _start(x1), _end(x2), _funcName(funcName),
                            _func(f1), _funcDerived(f2), _funcIntegrated(f3),
                            _funcExpr(funcExpr), _funcDerivedExpr(funcDerivedExpr), _funcIntegratedExpr(funcIntegratedExpr)
        {}
    };

    struct TestFunctionRealWithDerivation
    {
        double _start, _end;
        std::string _funcName;

        MML::RealFunction _func;
        MML::RealFunction _funcDerived;

        std::string _funcExpr;
        std::string _funcDerivedExpr;

        TestFunctionRealWithDerivation( double x1, double x2, std::string funcName,
                                        double (*f1)(double), std::string funcExpr, 
                                        double (*f2)(double), std::string funcDerivedExpr
                                        ) : _start(x1), _end(x2), _funcName(funcName),
                                            _func(f1), _funcDerived(f2),
                                            _funcExpr(funcExpr), _funcDerivedExpr(funcDerivedExpr)
        {}
    };

    struct TestFunctionRealWithIntegral
    {
        double _start, _end;
        std::string _funcName;

        MML::RealFunction _func;
        MML::RealFunction _funcIntegrated;

        std::string _funcExpr;
        std::string _funcIntegratedExpr;

        TestFunctionRealWithIntegral(double x1, double x2, std::string funcName,
                                    double (*f1)(double), std::string funcExpr, 
                                    double (*f2)(double), std::string funcIntegratedExpr
                                    ) : _start(x1), _end(x2), _funcName(funcName),
                                        _func(f1), _funcIntegrated(f2),
                                        _funcExpr(funcExpr), _funcIntegratedExpr(funcIntegratedExpr)
        {}
    };

    class RealFunctionsTestBed
    {
    public:
        const static inline TestFunctionReal _listFuncReal[] = { 
                {0.0, 1, "Sin(x)", [](double x) { return sin(x);},  "sin(x)", 
                                   [](double x) { return cos(x);},  "cos(x)",  
                                   [](double x) { return -cos(x);}, "-cos(x)"},
                {0.0, 1, "Cos(x)", [](double x) { return cos(x);}, "cos(x)", 
                                   [](double x) { return -sin(x);}, "-sin(x)", 
                                   [](double x) { return sin(x);}, "sin(x)"},
                {0.0, 1, "Sqrt(x)",[](double x) { return sqrt(x);}, "sqrt(x)",    
                                   [](double x) { return 0.5/sqrt(x);}, "1/(2 * sqrt(x))",
                                   [](double x) { return 2.0/3.0*x*sqrt(x);}, "2/3 * x^(3/2)"},
                {0.0, 1, "x^2",    [](double x) { return x*x;}, "x^2",    
                                   [](double x) { return 2*x;}, "2*x",
                                   [](double x) { return 1.0/3.0*x*x*x;}, "1/3 * x^3"},
                {0.0, 1, "x^3",    [](double x) { return x*x*x;}, "x^3",    
                                   [](double x) { return 3*x*x;}, "3*x^2",
                                   [](double x) { return 1.0/4.0*x*x*x*x;}, "1/4 * x^4"},
                {0.0, 1, "Exp(x)", [](double x) { return exp(x);}, "exp(x)", 
                                   [](double x) { return exp(x);}, "exp(x)",
                                   [](double x) { return exp(x);}, "exp(x)"}
            };

        const static inline TestFunctionRealWithDerivation _listFuncRealWithDerivation[] = { 
                {0.0, 1, "Sin(x)", [](double x) { return sin(x);},  "sin(x)", 
                                   [](double x) { return cos(x);},  "cos(x)"}
            };

        const static inline TestFunctionRealWithIntegral _listFuncRealWithIntegral[] = { 
                {0.0, 1, "Sin(x)", [](double x) { return sin(x);},  "sin(x)", 
                                   [](double x) { return -cos(x);},  "-cos(x)"}
            };
    };

}

#endif