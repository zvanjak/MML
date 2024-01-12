#if !defined __MML_REAL_FUNCTIONS_TEST_BED_H
#define __MML_REAL_FUNCTIONS_TEST_BED_H

#include <string>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "utilities/Intervals.h"
#include "utilities/StdFunctions.h"
#include "utilities/Constants.h"

#include "core/Function.h"
#endif

using namespace MML::Functions;

namespace MML::TestBeds
{
    struct TestFunctionReal
    {
        std::string _funcName;

        const IInterval &_intervalDef;
        const IInterval &_intervalTest;

        MML::RealFunction _func;
        MML::RealFunction _funcDerived;
        MML::RealFunction _funcSecDer;
        MML::RealFunction _funcThirdDer;
        MML::RealFunction _funcIntegrated;

        std::string _funcExpr;
        std::string _funcDerivedExpr;
        std::string _funcSecDerExpr;
        std::string _funcThirdDerExpr;
        std::string _funcIntegratedExpr;

        TestFunctionReal(std::string funcName, 
                        const IInterval &inIntervalDef,
                        double (*f1)(double), std::string funcExpr, 
                        double (*f2)(double), std::string funcDerivedExpr, 
                        double (*f3)(double), std::string funcSecDerExpr, 
                        double (*f4)(double), std::string funcThirdDerExpr, 
                        double (*f5)(double), std::string funcIntegratedExpr
                        ) : _funcName(funcName),
                            _intervalDef(inIntervalDef), _intervalTest(inIntervalDef),
                            _func(f1), _funcExpr(funcExpr),
                            _funcDerived(f2), _funcDerivedExpr(funcDerivedExpr),
                            _funcSecDer(f3), _funcSecDerExpr(funcSecDerExpr),
                            _funcThirdDer(f4), _funcThirdDerExpr(funcThirdDerExpr),
                            _funcIntegrated(f5), _funcIntegratedExpr(funcIntegratedExpr)
        {}

        TestFunctionReal(std::string funcName, 
                        const IInterval &inIntervalDef, const IInterval &inIntervalTest,
                        double (*f1)(double), std::string funcExpr, 
                        double (*f2)(double), std::string funcDerivedExpr, 
                        double (*f3)(double), std::string funcSecDerExpr, 
                        double (*f4)(double), std::string funcThirdDerExpr, 
                        double (*f5)(double), std::string funcIntegratedExpr
                        ) : _funcName(funcName),
                            _intervalDef(inIntervalDef), _intervalTest(inIntervalDef),
                            _func(f1), _funcExpr(funcExpr),
                            _funcDerived(f2), _funcDerivedExpr(funcDerivedExpr),
                            _funcSecDer(f3), _funcSecDerExpr(funcSecDerExpr),
                            _funcThirdDer(f4), _funcThirdDerExpr(funcThirdDerExpr),
                            _funcIntegrated(f5), _funcIntegratedExpr(funcIntegratedExpr)
        {}
    };

    struct TestFunctionRealWithDerivation
    {
        double _start, _end;
        std::string _funcName;

        const IInterval &_intervalDef;
        const IInterval &_intervalTest;        

        MML::RealFunction _func;
        MML::RealFunction _funcDerived;

        std::string _funcExpr;
        std::string _funcDerivedExpr;

        TestFunctionRealWithDerivation( std::string funcName, const IInterval &inIntervalDef,
                                        double (*f1)(double), std::string funcExpr, 
                                        double (*f2)(double), std::string funcDerivedExpr
                                        ) : _funcName(funcName),
                                            _intervalDef(inIntervalDef), _intervalTest(inIntervalDef),
                                            _func(f1), _funcDerived(f2),
                                            _funcExpr(funcExpr), _funcDerivedExpr(funcDerivedExpr)
        {}

        TestFunctionRealWithDerivation( std::string funcName, const IInterval &inIntervalDef, const IInterval &inIntervalTest,
                                        double (*f1)(double), std::string funcExpr, 
                                        double (*f2)(double), std::string funcDerivedExpr
                                        ) : _funcName(funcName),
                                            _intervalDef(inIntervalDef), _intervalTest(inIntervalTest),
                                            _func(f1), _funcDerived(f2),
                                            _funcExpr(funcExpr), _funcDerivedExpr(funcDerivedExpr)
        {}        
    };

    struct TestFunctionRealWithIntegral
    {
        double _start, _end;
        std::string _funcName;
        
        const IInterval &_intervalDef;
        const IInterval &_intervalTest;

        MML::RealFunction _func;
        MML::RealFunction _funcIntegrated;

        std::string _funcExpr;
        std::string _funcIntegratedExpr;

        TestFunctionRealWithIntegral(std::string funcName, const IInterval &inIntervalDef, 
                                    double (*f1)(double), std::string funcExpr, 
                                    double (*f2)(double), std::string funcIntegratedExpr
                                    ) : _funcName(funcName),
                                        _intervalDef(inIntervalDef), _intervalTest(inIntervalDef),
                                        _func(f1), _funcIntegrated(f2),
                                        _funcExpr(funcExpr), _funcIntegratedExpr(funcIntegratedExpr)
        {}

        TestFunctionRealWithIntegral(std::string funcName, const IInterval &inIntervalDef, const IInterval &inIntervalTest, 
                                    double (*f1)(double), std::string funcExpr, 
                                    double (*f2)(double), std::string funcIntegratedExpr
                                    ) : _funcName(funcName),
                                        _intervalDef(inIntervalDef), _intervalTest(inIntervalTest),
                                        _func(f1), _funcIntegrated(f2),
                                        _funcExpr(funcExpr), _funcIntegratedExpr(funcIntegratedExpr)
        {}
    };

    class RealFunctionsTestBed
    {
    public:
        static int getNumTestFunctionReal()                 { return 19; }
        static int getNumTestFunctionRealWithDerivation()   { return 5; }
        static int getNumTestFunctionRealWithIntegral()     { return 5; }
        
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
                {"Sin", CompleteRInterval(), 
                        ClosedInterval(0.0, 2.0 * MML::Constants::PI),  
                        [](double x) { return sin(x);},  "sin(x)", 
                        [](double x) { return cos(x);},  "cos(x)",  
                        [](double x) { return -sin(x);}, "-sin(x)",
                        [](double x) { return -cos(x);}, "-cos(x)",
                        [](double x) { return -cos(x);}, "-cos(x)"},
                {"Cos", CompleteRInterval(), 
                        ClosedInterval(0.0, 2.0 * MML::Constants::PI),  
                        [](double x) { return cos(x);}, "cos(x)", 
                        [](double x) { return -sin(x);}, "-sin(x)", 
                        [](double x) { return -cos(x);}, "-cos(x)", 
                        [](double x) { return sin(x);}, "sin(x)", 
                        [](double x) { return sin(x);}, "sin(x)"},
                {"Tan", Interval::Union(ClosedOpenInterval(-2.0*Constants::PI, -1.5*Constants::PI), 
                                        OpenInterval(-1.5*Constants::PI, -0.5*Constants::PI))
                                    .PerformUnion(OpenInterval(-0.5*Constants::PI, 0.5*Constants::PI)) 
                                    .PerformUnion(OpenInterval(0.5*Constants::PI, 1.5*Constants::PI))
                                    .PerformUnion(OpenClosedInterval(1.5*Constants::PI, 2.0*Constants::PI)), 
                        [](double x) { return tan(x);}, "tan(x)", 
                        [](double x) { return 1.0 / (cos(x) * cos(x));}, "1 / cos(x)^2", 
                        [](double x) { return 2.0 * tan(x) / (cos(x) * cos(x));}, "2*tan(x) / cos(x)^2", 
                        [](double x) { return -2 * (cos(2*x)-2) / pow(cos(x), 4);}, "-2*(cos(2*x)-2) / cos(x)^4)", 
                        [](double x) { return -log(std::abs(cos(x))); }, "-log(std::abs(cos(x)))"},
                {"Sinh", CompleteRInterval(), 
                        ClosedInterval(-10, 10),   
                        [](double x) { return sinh(x);}, "sinh(x)", 
                        [](double x) { return cosh(x);}, "cosh(x)",  
                        [](double x) { return sinh(x);}, "sinh(x)",  
                        [](double x) { return cosh(x);}, "cosh(x)",  
                        [](double x) { return cosh(x);}, "cosh(x)"},
                {"Cosh", CompleteRInterval(), 
                        ClosedInterval(-10, 10),   
                        [](double x) { return cosh(x);}, "cosh(x)", 
                        [](double x) { return sinh(x);}, "sinh(x)", 
                        [](double x) { return cosh(x);}, "cosh(x)", 
                        [](double x) { return sinh(x);}, "sinh(x)", 
                        [](double x) { return sinh(x);}, "sinh(x)"},
                {"Tanh", CompleteRInterval(), 
                        ClosedInterval(-10, 10),   
                        [](double x) { return tanh(x);}, "tanh(x)", 
                        [](double x) { return 1.0 / (cosh(x) * cosh(x));}, "1.0 / cosh(x)^2", 
                        [](double x) { return -2.0 * sinh(x) / pow(cosh(x), 3.0);}, "-2.0 * sinh(x) / cosh(x)^3", 
                        [](double x) { return -2.0 / pow(cosh(x), 4.0) + 4.0 * pow(sinh(x),2.0) / pow(cosh(x), 4.0) ;}, "-2.0 / cosh(x)^4 + 4 * sinh(x)^2 / cosh(x)^4", 
                        [](double x) { return log(cosh(x)); }, "log(cosh(x))"},
                {"Sqrt",OpenToInfInterval(0.0), 
                        OpenClosedInterval(0.0, 1e6 ),
                        [](double x) { return sqrt(x);}, "sqrt(x)",    
                        [](double x) { return 0.5/sqrt(x);}, "1 / (2 * sqrt(x))",
                        [](double x) { return -0.25 / pow(x, 1.5);}, "-1 / (4 * x^(3/2))",
                        [](double x) { return 3.0/(8.0 * pow(x, 2.5));}, "3 / (8 * x^(5/2)))",
                        [](double x) { return 2.0/3.0*x*sqrt(x);}, "2/3 * x^(3/2)"},
                {"x^2", CompleteRInterval(), 
                        ClosedInterval(-100.0, 100.0),  
                        [](double x) { return x*x;},  "x^2", 
                        [](double x) { return 2*x;},  "2*x",  
                        [](double x) { return 2.0;}, "2.0",
                        [](double x) { return 0.0;}, "0.0",
                        [](double x) { return 1.0/3.0 * x*x*x;}, "1.0/3.0 * x^3"},
                {"x^3", CompleteRInterval(), 
                        ClosedInterval(-100.0, 100.0),  
                        [](double x) { return x*x*x;},  "x^3", 
                        [](double x) { return 3*x*x;},  "3*x^2",  
                        [](double x) { return 6*x;}, "6*x",
                        [](double x) { return 6.0;}, "0.0",
                        [](double x) { return 0.25 * x*x*x*x;}, "0.25 * x^4"},
                {"x^4", CompleteRInterval(), 
                        ClosedInterval(-100.0, 100.0),  
                        [](double x) { return x*x*x*x;},  "x^4", 
                        [](double x) { return 4*x*x*x;},  "4*x^3",  
                        [](double x) { return 12*x*x;}, "12*x^2",
                        [](double x) { return 24*x;}, "24*x",
                        [](double x) { return 0.2 * x*x*x*x*x;}, "0.2 * x^5"},
                {"x^5", CompleteRInterval(), 
                        ClosedInterval(-100.0, 100.0),  
                        [](double x) { return x*x*x*x*x;},  "x^5", 
                        [](double x) { return 5*x*x*x*x;},  "5*x^4",  
                        [](double x) { return 20*x*x*x;}, "20*x^3",
                        [](double x) { return 60*x*x;}, "60*x^2",
                        [](double x) { return 1.0/6.0 * x*x*x*x*x*x;}, "1.0/6.0 * x^6"},
                {"Ln",  OpenToInfInterval(0.0), 
                        OpenClosedInterval(0.0, 1e6 ),
                        [](double x) { return log(x);}, "ln(x)", 
                        [](double x) { return 1.0 / x;}, "1 / x",
                        [](double x) { return 1.0 / (x * x);}, "1 / x^2",
                        [](double x) { return 2.0 / (x * x * x);}, "2 / x^3",
                        [](double x) { return x * (log(x) - 1);}, "x * (log(x) - 1)"},
                {"Exp", CompleteRInterval(), 
                        OpenInterval(-20.0, 20.0 ),
                        [](double x) { return exp(x);}, "exp(x)", 
                        [](double x) { return exp(x);}, "exp(x)",
                        [](double x) { return exp(x);}, "exp(x)",
                        [](double x) { return exp(x);}, "exp(x)",
                        [](double x) { return exp(x);}, "exp(x)"},
                {"Asin", OpenInterval(-1.0, 1.0), 
                        OpenInterval(-1.0, 1.0),
                        [](double x) { return asin(x);}, "asin(x)", 
                        [](double x) { return 1.0 / sqrt(1 - x*x);}, "1.0 / sqrt(1 - xˇ2)",
                        [](double x) { return x / pow(1 - x*x, 1.5);}, "x / (1 - x*x)^(3/2)",
                        [](double x) { return (2*x*x + 1) / pow(1 - x*x, 2.5);}, "(2*x^2 + 1) / (1 - x*x)^(5/2)",
                        [](double x) { return sqrt(1.0 - x*x) + x * asin(x);}, "sqrt(1.0 - x^2) + x * asin(x)"},
                {"Acos", OpenInterval(-1.0, 1.0), 
                        OpenInterval(-1.0, 1.0),
                        [](double x) { return acos(x);}, "", 
                        [](double x) { return -1.0 / sqrt(1 - x*x);}, "-1.0 / sqrt(1 - x^2)",
                        [](double x) { return -x / pow(1-x*x, 1.5);}, "-x / (1-x*x)^3/2",
                        [](double x) { return (-2*x*x - 1) / pow(1 - x*x, 2.5);}, "(-2*x^2 - 1) / (1-x*x)^5/2",
                        [](double x) { return x*acos(x) - sqrt(1 - x*x);}, "x*acos(x) - sqrt(1 - x^2)"}, 
                {"Atan", CompleteRInterval(), 
                        ClosedInterval(-100.0, 100.0 ),
                        [](double x) { return atan(x);}, "", 
                        [](double x) { return 1.0 / (1.0 + x*x);}, "1.0 / (1.0 + x^2)",
                        [](double x) { return -2.0*x / pow((1 + x*x),2);}, "-2.0*x / (1 + x*x)^2",
                        [](double x) { return (6.0*x*x - 2) / pow((1 + x*x),3);}, "(6.0*x*x - 2) / (1 + x*x)^3",
                        [](double x) { return x*atan(x) - 0.5 * log(x*x -1);}, "x*atan(x) - 0.5 * log(x*x -1)"}, 
                {"Asinh", CompleteRInterval(), 
                        ClosedInterval(-10.0, 10.0 ),
                        [](double x) { return asinh(x);}, "", 
                        [](double x) { return 1.0 / sqrt(x*x + 1);}, "1.0 / sqrt(x*x + 1)",
                        [](double x) { return -x / pow(x*x + 1, 1.5);}, "-x / (x*x + 1)^3/2",
                        [](double x) { return (2*x*x + 1) / pow(x*x + 1, 3.2);}, "(2*x*x + 1) / (x*x + 1)^5/2",
                        [](double x) { return x*asinh(x) - sqrt(x*x + 1);}, "x*asinh(x) - sqrt(x*x + 1)"},
                {"Acosh", Interval::Union(NegInfToOpenInterval(-1.0), 
                                          OpenToInfInterval(1.0)),
                        [](double x) { return acosh(x);}, "", 
                        [](double x) { return 1.0 / ( sqrt(x - 1) * sqrt(x + 1) );}, "1.0 / ( sqrt(x - 1) * sqrt(x + 1) )",
                        [](double x) { return -x / ( pow(x - 1, 1.5) * pow(x + 1, 1.5) );}, "-x / ( (x - 1)^3/2 * (x + 1)^3/2 )",
                        [](double x) { return (2*x*x + 1) / ( pow(x - 1, 2.5) * pow(x + 1, 2.5) );}, "(2*x^2 + 1) / ( (x - 1)^5/2 * (x + 1)^5/2 )",
                        [](double x) { return x*acosh(x) - sqrt(x-1) * sqrt(x+1);}, "x*acosh(x) - sqrt(x-1) * sqrt(x+1)"},
                {"Atanh", Interval::Union(NegInfToOpenInterval(-1.0), 
                                          OpenInterval(-1.0, 1.0))
                                          .PerformUnion(OpenToInfInterval(1.0)),
                        [](double x) { return atanh(x);}, "", 
                        [](double x) { return 1.0 / (1 - x*x);}, "1.0 / (1 - x^2)",
                        [](double x) { return 2*x / pow(1 - x*x, 2);}, "2*x / (1 - x^2)^2)",
                        [](double x) { return -2.0 * (3*x*x + 1) / pow(x*x - 1, 3);}, "-2.0 * (3*x^2 + 1) / (x*x - 1)^3",
                        [](double x) { return 0.5 * log(1 - x*x) + x * atanh(x);}, "0.5 * log(1 - x^2) + x * atanh(x)"}
            };

        // TODO - HIGH naći još bar 5 složenih funkcija za derivaciju
        const static inline TestFunctionRealWithDerivation _listFuncRealWithDerivation[] = { 
                {"TestDer1", CompleteRInterval(), 
                            OpenInterval(-20.0, 20.0 ),
                            [](double x) { return sin(x);},  "sin(x)", 
                            [](double x) { return cos(x);},  "cos(x)"},
                {"TestDer2", CompleteRInterval(), 
                            OpenInterval(-20.0, 20.0 ),
                            [](double x) { return sin(x);},  "sin(x)", 
                            [](double x) { return cos(x);},  "cos(x)"},
                {"TestDer3", CompleteRInterval(), 
                            OpenInterval(-20.0, 20.0 ),
                            [](double x) { return sin(x);},  "sin(x)", 
                            [](double x) { return cos(x);},  "cos(x)"},
                {"TestDer4", CompleteRInterval(), 
                            OpenInterval(-20.0, 20.0 ),
                            [](double x) { return sin(x);},  "sin(x)", 
                            [](double x) { return cos(x);},  "cos(x)"},
                {"TestDer5", CompleteRInterval(), 
                            OpenInterval(-20.0, 20.0 ),
                            [](double x) { return sin(x);},  "sin(x)", 
                            [](double x) { return cos(x);},  "cos(x)"}
            };

        // TODO - HIGH naći još bar 5 složenih funkcija za integraciju
        const static inline TestFunctionRealWithIntegral _listFuncRealWithIntegral[] = { 
                {"TestInt1", CompleteRInterval(), 
                            OpenInterval(-20.0, 20.0 ),
                            [](double x) { return x*x*(x*x-2.0)*sin(x);},  "x*x*(x*x-2.0)*sin(x)", 
                            [](double x) { return 4.0*x*(x*x-7.0)*sin(x)-(pow(x,4.0)-14.0*x*x+28.0)*cos(x);},  "4.0*x*(x*x-7.0)*sin(x)-(pow(x,4.0)-14.0*x*x+28.0)*cos(x)"},                                   
                {"TestInt2", CompleteRInterval(), 
                            OpenInterval(-20.0, 20.0 ),
                            [](double x) { return sin(x);},  "sin(x)", 
                            [](double x) { return -cos(x);},  "-cos(x)"},
                {"TestInt3", CompleteRInterval(), 
                            OpenInterval(-20.0, 20.0 ),
                            [](double x) { return sin(x);},  "sin(x)", 
                            [](double x) { return -cos(x);},  "-cos(x)"},
                {"TestInt14", CompleteRInterval(), 
                            OpenInterval(-20.0, 20.0 ),
                            [](double x) { return sin(x);},  "sin(x)", 
                            [](double x) { return -cos(x);},  "-cos(x)"},
                {"TestInt5", CompleteRInterval(), 
                            OpenInterval(-20.0, 20.0 ),
                            [](double x) { return sin(x);},  "sin(x)", 
                            [](double x) { return -cos(x);},  "-cos(x)"}                                                                        
            };
    };
}
#endif
