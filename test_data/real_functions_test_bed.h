#if !defined __MML_REAL_FUNCTIONS_TEST_BED_H
#define __MML_REAL_FUNCTIONS_TEST_BED_H

#include <string>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/StdFunctions.h"

#include "base/Intervals.h"

#include "core/Function.h"
#endif

using namespace MML::Functions;

namespace MML::TestBeds
{
    struct TestFunctionReal
    {
        std::string _funcName;

        std::shared_ptr<IInterval> _intervalDef;
        std::shared_ptr<IInterval> _intervalTest;

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

        // TestFunctionReal(std::string funcName, 
        //                 IInterval *inIntervalDef,
        //                 Real (*f1)(Real), std::string funcExpr, 
        //                 Real (*f2)(Real), std::string funcDerivedExpr, 
        //                 Real (*f3)(Real), std::string funcSecDerExpr, 
        //                 Real (*f4)(Real), std::string funcThirdDerExpr, 
        //                 Real (*f5)(Real), std::string funcIntegratedExpr
        //                 ) : _funcName(funcName),
        //                     _intervalDef(inIntervalDef), _intervalTest(inIntervalDef),
        //                     _func(f1), _funcExpr(funcExpr),
        //                     _funcDerived(f2), _funcDerivedExpr(funcDerivedExpr),
        //                     _funcSecDer(f3), _funcSecDerExpr(funcSecDerExpr),
        //                     _funcThirdDer(f4), _funcThirdDerExpr(funcThirdDerExpr),
        //                     _funcIntegrated(f5), _funcIntegratedExpr(funcIntegratedExpr)
        // {
        //  moze ako se za intervalTest napravi KOPIJA ulaznog intervala!
        // }

        TestFunctionReal(std::string funcName, 
                        IInterval *inIntervalDef, IInterval *inIntervalTest,
                        Real (*f1)(Real), std::string funcExpr, 
                        Real (*f2)(Real), std::string funcDerivedExpr, 
                        Real (*f3)(Real), std::string funcSecDerExpr, 
                        Real (*f4)(Real), std::string funcThirdDerExpr, 
                        Real (*f5)(Real), std::string funcIntegratedExpr
                        ) : _funcName(funcName),
                            _func(f1), _funcExpr(funcExpr),
                            _funcDerived(f2), _funcDerivedExpr(funcDerivedExpr),
                            _funcSecDer(f3), _funcSecDerExpr(funcSecDerExpr),
                            _funcThirdDer(f4), _funcThirdDerExpr(funcThirdDerExpr),
                            _funcIntegrated(f5), _funcIntegratedExpr(funcIntegratedExpr)
        {
            _intervalDef = std::shared_ptr<IInterval>(inIntervalDef);
            _intervalTest = std::shared_ptr<IInterval>(inIntervalTest);
        }
    };

    struct TestFunctionRealWithDerivation
    {
        double _start, _end;
        std::string _funcName;

        std::shared_ptr<IInterval> _intervalDef;
        std::shared_ptr<IInterval> _intervalTest;        

        MML::RealFunction _func;
        MML::RealFunction _funcDerived;

        std::string _funcExpr;
        std::string _funcDerivedExpr;

        TestFunctionRealWithDerivation( std::string funcName, IInterval *inIntervalDef,
                                        Real (*f1)(Real), std::string funcExpr, 
                                        Real (*f2)(Real), std::string funcDerivedExpr
                                        ) : _funcName(funcName),
                                            _intervalDef(inIntervalDef), _intervalTest(inIntervalDef),
                                            _func(f1), _funcDerived(f2),
                                            _funcExpr(funcExpr), _funcDerivedExpr(funcDerivedExpr)
        {}

        TestFunctionRealWithDerivation( std::string funcName, IInterval *inIntervalDef, IInterval *inIntervalTest,
                                        Real (*f1)(Real), std::string funcExpr, 
                                        Real (*f2)(Real), std::string funcDerivedExpr
                                        ) : _funcName(funcName),
                                            _intervalDef(inIntervalDef), _intervalTest(inIntervalTest),
                                            _func(f1), _funcDerived(f2),
                                            _funcExpr(funcExpr), _funcDerivedExpr(funcDerivedExpr)
        {}        
    };

    struct TestFunctionRealWithIntegral
    {
        Real _start, _end;
        std::string _funcName;
        
        std::shared_ptr<IInterval> _intervalDef;
        std::shared_ptr<IInterval> _intervalTest;

        MML::RealFunction _func;
        MML::RealFunction _funcIntegrated;

        std::string _funcExpr;
        std::string _funcIntegratedExpr;

        TestFunctionRealWithIntegral(std::string funcName, IInterval *inIntervalDef, 
                                    Real (*f1)(Real), std::string funcExpr, 
                                    Real (*f2)(Real), std::string funcIntegratedExpr
                                    ) : _funcName(funcName),
                                        _intervalDef(inIntervalDef), _intervalTest(inIntervalDef),
                                        _func(f1), _funcIntegrated(f2),
                                        _funcExpr(funcExpr), _funcIntegratedExpr(funcIntegratedExpr)
        {}
       
        TestFunctionRealWithIntegral(std::string funcName, IInterval *inIntervalDef, IInterval *inIntervalTest, 
                                    Real (*f1)(Real), std::string funcExpr, 
                                    Real (*f2)(Real), std::string funcIntegratedExpr
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
                {"Sin", new CompleteRInterval(), 
                        new ClosedInterval(-2.0*Constants::PI, 2.0 * MML::Constants::PI),  
                        [](Real x) { return sin(x);},  "sin(x)", 
                        [](Real x) { return cos(x);},  "cos(x)",  
                        [](Real x) { return -sin(x);}, "-sin(x)",
                        [](Real x) { return -cos(x);}, "-cos(x)",
                        [](Real x) { return -cos(x);}, "-cos(x)"},
                {"Cos", new CompleteRInterval(), 
                        new ClosedInterval(-2.0*Constants::PI, 2.0 * MML::Constants::PI),  
                        [](Real x) { return cos(x);}, "cos(x)", 
                        [](Real x) { return -sin(x);}, "-sin(x)", 
                        [](Real x) { return -cos(x);}, "-cos(x)", 
                        [](Real x) { return sin(x);}, "sin(x)", 
                        [](Real x) { return sin(x);}, "sin(x)"},
                {"Tan", new CompleteRWithReccuringPointHoles(0.5*Constants::PI, Constants::PI),
                        new ClosedIntervalWithReccuringPointHoles(-5*Constants::PI, 5*Constants::PI, 0.5*Constants::PI, Constants::PI),
                        [](Real x) { return tan(x);}, "tan(x)", 
                        [](Real x) { return 1 / (cos(x) * cos(x));}, "1 / cos(x)^2", 
                        [](Real x) { return 2 * tan(x) / (cos(x) * cos(x));}, "2*tan(x) / cos(x)^2", 
                        [](Real x) { return -2 * (cos(2*x)-2) / (Real) pow(cos(x), 4);}, "-2*(cos(2*x)-2) / cos(x)^4)", 
                        [](Real x) { return -log(std::abs(cos(x))); }, "-log(std::abs(cos(x)))"},
                {"Sinh", new CompleteRInterval(), 
                        new ClosedInterval(-10, 10),   
                        [](Real x) { return sinh(x);}, "sinh(x)", 
                        [](Real x) { return cosh(x);}, "cosh(x)",  
                        [](Real x) { return sinh(x);}, "sinh(x)",  
                        [](Real x) { return cosh(x);}, "cosh(x)",  
                        [](Real x) { return cosh(x);}, "cosh(x)"},
                {"Cosh", new CompleteRInterval(), 
                        new ClosedInterval(-10, 10),   
                        [](Real x) { return cosh(x);}, "cosh(x)", 
                        [](Real x) { return sinh(x);}, "sinh(x)", 
                        [](Real x) { return cosh(x);}, "cosh(x)", 
                        [](Real x) { return sinh(x);}, "sinh(x)", 
                        [](Real x) { return sinh(x);}, "sinh(x)"},
                {"Tanh", new CompleteRInterval(), 
                        new ClosedInterval(-10, 10),   
                        [](Real x) { return tanh(x);}, "tanh(x)", 
                        [](Real x) { return 1 / (cosh(x) * cosh(x));}, "1.0 / cosh(x)^2", 
                        [](Real x) { return -2 * sinh(x) / (Real) pow(cosh(x), 3.0);}, "-2.0 * sinh(x) / cosh(x)^3", 
                        [](Real x) { return -2 / (Real) pow(cosh(x), 4) + 4 * (Real) pow(sinh(x),2.0) / (Real) pow(cosh(x), 4.0) ;}, "-2.0 / cosh(x)^4 + 4 * sinh(x)^2 / cosh(x)^4", 
                        [](Real x) { return log(cosh(x)); }, "log(cosh(x))"},
                {"Sqrt", new OpenToInfInterval(0.0), 
                        new OpenClosedInterval(0.0, 1e6 ),
                        [](Real x) { return sqrt(x);}, "sqrt(x)",    
                        [](Real x) { return 1 / (2 * sqrt(x));}, "1 / (2 * sqrt(x))",
                        [](Real x) { return Real{-0.25} / (Real) pow(x, 1.5);}, "-1 / (4 * x^(3/2))",
                        [](Real x) { return 3/(8 * (Real) pow(x, 2.5));}, "3 / (8 * x^(5/2)))",
                        [](Real x) { return Real{2/3.0} * x*sqrt(x);}, "2/3 * x^(3/2)"},
                {"x^2", new CompleteRInterval(), 
                        new ClosedInterval(-100.0, 100.0),  
                        [](Real x) { return x*x;},  "x^2", 
                        [](Real x) { return 2*x;},  "2*x",  
                        [](Real x) { return Real{2};}, "2.0",
                        [](Real x) { return Real{0};}, "0.0",
                        [](Real x) { return Real{1}/3 * x*x*x;}, "1.0/3.0 * x^3"},
                {"x^3", new CompleteRInterval(), 
                        new ClosedInterval(-100.0, 100.0),  
                        [](Real x) { return x*x*x;},  "x^3", 
                        [](Real x) { return 3*x*x;},  "3*x^2",  
                        [](Real x) { return 6*x;}, "6*x",
                        [](Real x) { return Real{6};}, "0.0",
                        [](Real x) { return Real{0.25} * x*x*x*x;}, "0.25 * x^4"},
                {"x^4", new CompleteRInterval(), 
                        new ClosedInterval(-100.0, 100.0),  
                        [](Real x) { return x*x*x*x;},  "x^4", 
                        [](Real x) { return 4*x*x*x;},  "4*x^3",  
                        [](Real x) { return 12*x*x;}, "12*x^2",
                        [](Real x) { return 24*x;}, "24*x",
                        [](Real x) { return Real{0.2} * x*x*x*x*x;}, "0.2 * x^5"},
                {"x^5", new CompleteRInterval(), 
                        new ClosedInterval(-100.0, 100.0),  
                        [](Real x) { return x*x*x*x*x;},  "x^5", 
                        [](Real x) { return 5*x*x*x*x;},  "5*x^4",  
                        [](Real x) { return 20*x*x*x;}, "20*x^3",
                        [](Real x) { return 60*x*x;}, "60*x^2",
                        [](Real x) { return Real{1.0/6.0} * x*x*x*x*x*x;}, "1.0/6.0 * x^6"},
                {"Ln",  new OpenToInfInterval(0.0), 
                        new OpenClosedInterval(0.0, 1000.0 ),
                        [](Real x) { return log(x);}, "ln(x)", 
                        [](Real x) { return 1 / x;}, "1 / x",
                        [](Real x) { return 1 / (x * x);}, "1 / x^2",
                        [](Real x) { return 2 / (x * x * x);}, "2 / x^3",
                        [](Real x) { return x * (log(x) - 1);}, "x * (log(x) - 1)"},
                {"Exp", new CompleteRInterval(), 
                        new OpenInterval(-20.0, 20.0 ),
                        [](Real x) { return exp(x);}, "exp(x)", 
                        [](Real x) { return exp(x);}, "exp(x)",
                        [](Real x) { return exp(x);}, "exp(x)",
                        [](Real x) { return exp(x);}, "exp(x)",
                        [](Real x) { return exp(x);}, "exp(x)"},
                {"Asin", new OpenInterval(-1.0, 1.0), 
                        new OpenInterval(-1.0, 1.0),
                        [](Real x) { return asin(x);}, "asin(x)", 
                        [](Real x) { return 1 / sqrt(1 - x*x);}, "1.0 / sqrt(1 - xˇ2)",
                        [](Real x) { return x / (Real) pow(1 - x*x, 1.5);}, "x / (1 - x*x)^(3/2)",
                        [](Real x) { return (2*x*x + 1) / (Real) pow(1 - x*x, 2.5);}, "(2*x^2 + 1) / (1 - x*x)^(5/2)",
                        [](Real x) { return sqrt(1 - x*x) + x * asin(x);}, "sqrt(1.0 - x^2) + x * asin(x)"},
                {"Acos", new OpenInterval(-1.0, 1.0), 
                        new OpenInterval(-1.0, 1.0),
                        [](Real x) { return acos(x);}, "", 
                        [](Real x) { return -1 / sqrt(1 - x*x);}, "-1.0 / sqrt(1 - x^2)",
                        [](Real x) { return -x / (Real) pow(1-x*x, 1.5);}, "-x / (1-x*x)^3/2",
                        [](Real x) { return (-2*x*x - 1) / (Real) pow(1 - x*x, 2.5);}, "(-2*x^2 - 1) / (1-x*x)^5/2",
                        [](Real x) { return x*acos(x) - sqrt(1 - x*x);}, "x*acos(x) - sqrt(1 - x^2)"}, 
                {"Atan", new CompleteRInterval(), 
                        new ClosedInterval(-100.0, 100.0 ),
                        [](Real x) { return atan(x);}, "", 
                        [](Real x) { return 1 / (1 + x*x);}, "1.0 / (1.0 + x^2)",
                        [](Real x) { return -2*x / (Real) pow((1 + x*x),2);}, "-2.0*x / (1 + x*x)^2",
                        [](Real x) { return (6*x*x - 2) / (Real) pow((1 + x*x),3);}, "(6.0*x*x - 2) / (1 + x*x)^3",
                        [](Real x) { return x*atan(x) - log(x*x -1) / 2;}, "x*atan(x) - 0.5 * log(x*x -1)"}, 
                {"Asinh", new CompleteRInterval(), 
                        new ClosedInterval(-10.0, 10.0 ),
                        [](Real x) { return asinh(x);}, "", 
                        [](Real x) { return 1 / sqrt(x*x + 1);}, "1.0 / sqrt(x*x + 1)",
                        [](Real x) { return -x / (Real) pow(x*x + 1, 1.5);}, "-x / (x*x + 1)^3/2",
                        [](Real x) { return (2*x*x + 1) / (Real) pow(x*x + 1, 3.2);}, "(2*x*x + 1) / (x*x + 1)^5/2",
                        [](Real x) { return x*asinh(x) - sqrt(x*x + 1);}, "x*asinh(x) - sqrt(x*x + 1)"},
                {"Acosh", new Interval({ new NegInfToOpenInterval(-1.0), 
                                         new OpenToInfInterval(1.0) } ),
                          new Interval({ new ClosedOpenInterval(-100.0, -1.0), 
                                         new OpenClosedInterval(1.0, 100.0) } ),                        
                        [](Real x) { return acosh(x);}, "", 
                        [](Real x) { return 1 / ( sqrt(x - 1) * sqrt(x + 1) );}, "1.0 / ( sqrt(x - 1) * sqrt(x + 1) )",
                        [](Real x) { return -x / (Real) ( pow(x - 1, 1.5) * pow(x + 1, 1.5) );}, "-x / ( (x - 1)^3/2 * (x + 1)^3/2 )",
                        [](Real x) { return (2*x*x + 1) / (Real) ( pow(x - 1, 2.5) * pow(x + 1, 2.5) );}, "(2*x^2 + 1) / ( (x - 1)^5/2 * (x + 1)^5/2 )",
                        [](Real x) { return x*acosh(x) - sqrt(x-1) * sqrt(x+1);}, "x*acosh(x) - sqrt(x-1) * sqrt(x+1)"},
                {"Atanh", new Interval({ new NegInfToOpenInterval(-1.0), 
                                          new OpenInterval(-1.0, 1.0),
                                          new OpenToInfInterval(1.0) } ),
                          new Interval({ new ClosedOpenInterval(-100.0, -1.0), 
                                          new OpenInterval(-1.0, 1.0),
                                          new OpenClosedInterval(1.0, 100.0) } ),
                        [](Real x) { return atanh(x);}, "", 
                        [](Real x) { return 1 / (1 - x*x);}, "1.0 / (1 - x^2)",
                        [](Real x) { return 2*x / (Real) pow(1 - x*x, 2);}, "2*x / (1 - x^2)^2)",
                        [](Real x) { return -2 * (3*x*x + 1) / (Real) pow(x*x - 1, 3);}, "-2.0 * (3*x^2 + 1) / (x*x - 1)^3",
                        [](Real x) { return log(1 - x*x) / 2 + x * atanh(x);}, "0.5 * log(1 - x^2) + x * atanh(x)"}
            };

        // TODO 0.9 - HIGH naći još bar 5 složenih funkcija za derivaciju
        const static inline TestFunctionRealWithDerivation _listFuncRealWithDerivation[] = { 
                {"TestDer1", new CompleteRInterval(), 
                            new OpenInterval(-20.0, 20.0 ),
                            [](Real x) { return sin(x);},  "sin(x)", 
                            [](Real x) { return cos(x);},  "cos(x)"},
                {"TestDer2", new CompleteRInterval(), 
                            new OpenInterval(-20.0, 20.0 ),
                            [](Real x) { return sin(x);},  "sin(x)", 
                            [](Real x) { return cos(x);},  "cos(x)"},
                {"TestDer3", new CompleteRInterval(), 
                            new OpenInterval(-20.0, 20.0 ),
                            [](Real x) { return sin(x);},  "sin(x)", 
                            [](Real x) { return cos(x);},  "cos(x)"},
                {"TestDer4", new CompleteRInterval(), 
                            new OpenInterval(-20.0, 20.0 ),
                            [](Real x) { return sin(x);},  "sin(x)", 
                            [](Real x) { return cos(x);},  "cos(x)"},
                {"TestDer5", new CompleteRInterval(), 
                            new OpenInterval(-20.0, 20.0 ),
                            [](Real x) { return sin(x);},  "sin(x)", 
                            [](Real x) { return cos(x);},  "cos(x)"}
            };

        // TODO 0.9 - HIGH naći još bar 5 složenih funkcija za integraciju
        const static inline TestFunctionRealWithIntegral _listFuncRealWithIntegral[] = { 
                {"TestInt1", new CompleteRInterval(), 
                            new OpenInterval(0.0, 5.0 ),
                            [](Real x) { return x*x*(x*x-2)*sin(x);},  "x*x*(x*x-2.0)*sin(x)", 
                            [](Real x) { return 4*x*(x*x-7)*sin(x)-((Real) pow(x,4.0)-14*x*x+28)*cos(x);},  "4.0*x*(x*x-7.0)*sin(x)-(pow(x,4.0)-14.0*x*x+28.0)*cos(x)"},                                   
                {"TestInt2", new CompleteRInterval(), 
                            new OpenInterval(-20.0, 20.0 ),
                            [](Real x) { return sin(x);},  "sin(x)", 
                            [](Real x) { return -cos(x);},  "-cos(x)"},
                {"TestInt3", new CompleteRInterval(), 
                            new OpenInterval(-20.0, 20.0 ),
                            [](Real x) { return sin(x);},  "sin(x)", 
                            [](Real x) { return -cos(x);},  "-cos(x)"},
                {"TestInt14", new CompleteRInterval(), 
                            new OpenInterval(-20.0, 20.0 ),
                            [](Real x) { return sin(x);},  "sin(x)", 
                            [](Real x) { return -cos(x);},  "-cos(x)"},
                {"TestInt5", new CompleteRInterval(), 
                            new OpenInterval(-20.0, 20.0 ),
                            [](Real x) { return sin(x);},  "sin(x)", 
                            [](Real x) { return -cos(x);},  "-cos(x)"}                                                                        
            };
    };
}
#endif
