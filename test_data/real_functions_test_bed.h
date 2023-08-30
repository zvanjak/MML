#if !defined __MML_REAL_FUNCTIONS_TEST_BED_H
#define __MML_REAL_FUNCTIONS_TEST_BED_H

#include <string>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Constants.h"
#include "basic_types/Function.h"
#endif

namespace MML::TestData
{
    // dodati test beds - u svakom po 15ak funkcija
    // TODO - base test bed - do trece derivacije i integral poznati
    // TODO - derivation test bed - 
    // TODO - integration test bed
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

        TestFunctionReal(std::string funcName, double x1, double x2,
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

        TestFunctionRealWithDerivation( std::string funcName, double x1, double x2, 
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

        TestFunctionRealWithIntegral(std::string funcName, double x1, double x2, 
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
        static int getNumTestFunctionReal()                 { return sizeof(_listFuncReal)/sizeof(TestFunctionReal); }
        static int getNumTestFunctionRealWithDerivation()   { return 1; }
        static int getNumTestFunctionRealWithIntegral()     { return 2; }
        
        const static TestFunctionReal& getTestFunctionReal(int i)                               { return _listFuncReal[i]; }
        const static TestFunctionRealWithDerivation& getTestFunctionRealWithDerivation(int i)   { return _listFuncRealWithDerivation[i]; }
        const static TestFunctionRealWithIntegral& getTestFunctionRealWithIntegral(int i)       { return _listFuncRealWithIntegral[i]; }

        const static TestFunctionReal& getTestFunctionReal(const std::string &funcName)
        {
            for (int i = 0; i < getNumTestFunctionReal(); i++)
            {
                if (_listFuncReal[i]._funcName == funcName)
                    return _listFuncReal[i];
            }
            throw std::runtime_error("TestFunctionReal " + funcName + " not found!");
        }
        const static TestFunctionRealWithDerivation& getTestFunctionRealWithDerivation(const std::string &funcName)
        {
            for (int i = 0; i < getNumTestFunctionRealWithDerivation(); i++)
            {
                if (_listFuncRealWithDerivation[i]._funcName == funcName)
                    return _listFuncRealWithDerivation[i];
            }
            throw std::runtime_error("TestFunctionRealWithDerivation " + funcName + " not found!");
        }
        const static TestFunctionRealWithIntegral& getTestFunctionRealWithIntegral(const std::string &funcName)
        {
            for (int i = 0; i < getNumTestFunctionRealWithIntegral(); i++)
            {
                if (_listFuncRealWithIntegral[i]._funcName == funcName)
                    return _listFuncRealWithIntegral[i];
            }
            throw std::runtime_error("TestFunctionRealWithIntegral " + funcName + " not found!");
        }

    private:
        const static inline TestFunctionReal _listFuncReal[] = { 
                {"Sin", 0.0, 2.0 * MML::Constants::PI,  
                        [](double x) { return sin(x);},  "sin(x)", 
                        [](double x) { return cos(x);},  "cos(x)",  
                        [](double x) { return -cos(x);}, "-cos(x)"},
                {"Cos", 0.0, 2.0 * MML::Constants::PI,  
                        [](double x) { return cos(x);}, "cos(x)", 
                        [](double x) { return -sin(x);}, "-sin(x)", 
                        [](double x) { return sin(x);}, "sin(x)"},
                {"Tan", 0.0, 2.0 * MML::Constants::PI,  
                        [](double x) { return tan(x);}, "tan(x)", 
                        [](double x) { return -sin(x);}, "???", 
                        [](double x) { return -log(std::abs(cos(x))); }, "-log(std::abs(cos(x)))"},
                {"Sinh", 0.0, 2.0 * MML::Constants::PI,  
                        [](double x) { return sinh(x);}, "sinh(x)", 
                        [](double x) { return cosh(x);}, "cosh(x)",  
                        [](double x) { return cosh(x);}, "cosh(x)"},
                {"Cosh", 0.0, 2.0 * MML::Constants::PI,  
                        [](double x) { return cosh(x);}, "cosh(x)", 
                        [](double x) { return sinh(x);}, "sinh(x)", 
                        [](double x) { return sinh(x);}, "sinh(x)"},
                {"Tanh", 0.0, 2.0 * MML::Constants::PI,  
                        [](double x) { return tanh(x);}, "tanh(x)", 
                        [](double x) { return -sin(x);}, "???", 
                        [](double x) { return log(cosh(x)); }, "log(cosh(x))"},
                {"Sqrt", 0.0, 1e6, 
                        [](double x) { return sqrt(x);}, "sqrt(x)",    
                        [](double x) { return 0.5/sqrt(x);}, "1/(2 * sqrt(x))",
                        [](double x) { return 2.0/3.0*x*sqrt(x);}, "2/3 * x^(3/2)"},
                {"x^2", -1e3, 1e3,
                        [](double x) { return x*x;}, "x^2",    
                        [](double x) { return 2*x;}, "2*x",
                        [](double x) { return 1.0/3.0*x*x*x;}, "1/3 * x^3"},
                {"x^3", -1e3, 1e3,
                        [](double x) { return x*x*x;}, "x^3",    
                        [](double x) { return 3*x*x;}, "3*x^2",
                        [](double x) { return 1.0/4.0*x*x*x*x;}, "1/4 * x^4"},
                {"Exp", -20, 20,
                        [](double x) { return exp(x);}, "exp(x)", 
                        [](double x) { return exp(x);}, "exp(x)",
                        [](double x) { return exp(x);}, "exp(x)"}
            };

        // TODO - naći još bar 5 složenih funkcija za derivaciju
        const static inline TestFunctionRealWithDerivation _listFuncRealWithDerivation[] = { 
                {"Sin", 0.0, 1,  [](double x) { return sin(x);},  "sin(x)", 
                                 [](double x) { return cos(x);},  "cos(x)"}
            };

        // TODO - naći još bar 5 složenih funkcija za integraciju
        const static inline TestFunctionRealWithIntegral _listFuncRealWithIntegral[] = { 
                {"Sin", 0.0, 1, [](double x) { return sin(x);},  "sin(x)", 
                                [](double x) { return -cos(x);},  "-cos(x)"},
                {"Test1", 0.0, 10.0, [](double x) { return x*x*(x*x-2.0)*sin(x);},  "x*x*(x*x-2.0)*sin(x)", 
                                     [](double x) { return 4.0*x*(x*x-7.0)*sin(x)-(pow(x,4.0)-14.0*x*x+28.0)*cos(x);},  "4.0*x*(x*x-7.0)*sin(x)-(pow(x,4.0)-14.0*x*x+28.0)*cos(x)"}                                   
            };
    };

}

#endif