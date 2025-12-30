#if !defined __MML_REAL_FUNCTIONS_TEST_BED_H
#define __MML_REAL_FUNCTIONS_TEST_BED_H

#include <string>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/StandardFunctions.h"

#include "base/Intervals.h"

#include "base/Function.h"
#endif

using namespace MML::Functions;

namespace MML::TestBeds
{
    /*******************************************************************************************************************
     * FUNCTION CATEGORY ENUMERATION
     *******************************************************************************************************************/
    
    enum class RealFunctionCategory
    {
        Standard,           // Standard mathematical functions (sin, cos, exp, log, polynomials)
        NumericalAnalysis,  // Functions interesting for numerical analysis testing
        Derivation,         // Functions specifically for derivative testing
        Integration         // Functions specifically for integral testing
    };

    /*******************************************************************************************************************
     * TEST FUNCTION STRUCTURES
     *******************************************************************************************************************/

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
        std::string _funcName;

        std::shared_ptr<IInterval> _intervalDef;
        std::shared_ptr<IInterval> _intervalTest;        

        MML::RealFunctionFromStdFunc _func;
        MML::RealFunctionFromStdFunc _funcDerived;

        std::string _funcExpr;
        std::string _funcDerivedExpr;

        TestFunctionRealWithDerivation( std::string funcName, IInterval *inIntervalDef,
                                        const std::function<Real(Real)>& f1, std::string funcExpr, 
                                        const std::function<Real(Real)>& f2, std::string funcDerivedExpr
                                        ) : _funcName(funcName),
                                            _intervalDef(inIntervalDef), _intervalTest(inIntervalDef),
                                            _func(f1), _funcDerived(f2),
                                            _funcExpr(funcExpr), _funcDerivedExpr(funcDerivedExpr)
        {}

        TestFunctionRealWithDerivation( std::string funcName, IInterval *inIntervalDef, IInterval *inIntervalTest,
                                        const std::function<Real(Real)>& f1, std::string funcExpr, 
                                        const std::function<Real(Real)>& f2, std::string funcDerivedExpr
                                        ) : _funcName(funcName),
                                            _intervalDef(inIntervalDef), _intervalTest(inIntervalTest),
                                            _func(f1), _funcDerived(f2),
                                            _funcExpr(funcExpr), _funcDerivedExpr(funcDerivedExpr)
        {}        
    };

    struct TestFunctionRealWithIntegral
    {
        std::string _funcName;
        
        std::shared_ptr<IInterval> _intervalDef;
        std::shared_ptr<IInterval> _intervalTest;

        MML::RealFunctionFromStdFunc _func;
        MML::RealFunctionFromStdFunc _funcIntegrated;

        std::string _funcExpr;
        std::string _funcIntegratedExpr;

        TestFunctionRealWithIntegral(std::string funcName, IInterval *inIntervalDef, 
                                    const std::function<Real(Real)>& f1, std::string funcExpr, 
                                    const std::function<Real(Real)>& f2, std::string funcIntegratedExpr
                                    ) : _funcName(funcName),
                                        _intervalDef(inIntervalDef), _intervalTest(inIntervalDef),
                                        _func(f1), _funcIntegrated(f2),
                                        _funcExpr(funcExpr), _funcIntegratedExpr(funcIntegratedExpr)
        {}
       
        TestFunctionRealWithIntegral(std::string funcName, IInterval *inIntervalDef, IInterval *inIntervalTest, 
                                    const std::function<Real(Real)>& f1, std::string funcExpr, 
                                    const std::function<Real(Real)>& f2, std::string funcIntegratedExpr
                                    ) : _funcName(funcName),
                                        _intervalDef(inIntervalDef), _intervalTest(inIntervalTest),
                                        _func(f1), _funcIntegrated(f2),
                                        _funcExpr(funcExpr), _funcIntegratedExpr(funcIntegratedExpr)
        {}
    };

    /*******************************************************************************************************************
     * REAL FUNCTIONS TEST BED
     * 
     * Provides test functions organized by category:
     *   - Standard: Basic mathematical functions (sin, cos, exp, log, polynomials, etc.)
     *   - Numerical Analysis: Functions specifically chosen for their challenging numerical properties
     *   - Derivation: Functions with known derivatives for testing differentiation algorithms
     *   - Integration: Functions with known antiderivatives for testing integration algorithms
     *******************************************************************************************************************/

    class RealFunctionsTestBed
    {
    public:
        //--------------------------------------------------------------------------------------------------------------
        // COUNTS
        //--------------------------------------------------------------------------------------------------------------
        
        static int getNumFunc()                     { return 19; }
        static int getNumFuncExtended()             { return 8; }
        static int getNumFuncWithDerivation()       { return 11; }
        static int getNumFuncWithIntegral()         { return 16; }
        static int getNumFuncExtendedWithIntegral() { return 6; }
        
        //--------------------------------------------------------------------------------------------------------------
        // STANDARD FUNCTIONS ACCESSORS
        //--------------------------------------------------------------------------------------------------------------
        
        const static TestFunctionReal& getFunc(int i)                               { return _listFuncReal[i]; }
        const static TestFunctionRealWithDerivation& getFuncWithDerivation(int i)   { return _listFuncRealWithDerivation[i]; }
        const static TestFunctionRealWithIntegral& getFuncWithIntegral(int i)       { return _listFuncRealWithIntegral[i]; }

        const static TestFunctionReal& getFunc(const std::string &funcName)
        {
          for (int i = 0; i < getNumFunc(); i++)
          {
              if (_listFuncReal[i]._funcName == funcName)
                  return _listFuncReal[i];
          }
          throw std::runtime_error("TestFunctionReal " + funcName + " not found!");
        }
        const static TestFunctionRealWithDerivation& getFuncWithDerivation(const std::string &funcName)
        {
            for (int i = 0; i < getNumFuncWithDerivation(); i++)
            {
                if (_listFuncRealWithDerivation[i]._funcName == funcName)
                    return _listFuncRealWithDerivation[i];
            }
            throw std::runtime_error("TestFunctionRealWithDerivation " + funcName + " not found!");
        }
        const static TestFunctionRealWithIntegral& getFuncWithIntegral(const std::string &funcName)
        {
            for (int i = 0; i < getNumFuncWithIntegral(); i++)
            {
                if (_listFuncRealWithIntegral[i]._funcName == funcName)
                    return _listFuncRealWithIntegral[i];
            }
            throw std::runtime_error("TestFunctionRealWithIntegral " + funcName + " not found!");
        }

        //--------------------------------------------------------------------------------------------------------------
        // EXTENDED (NUMERICAL ANALYSIS) FUNCTIONS ACCESSORS
        //--------------------------------------------------------------------------------------------------------------
        
        const static TestFunctionRealWithDerivation& getFuncExtended(int i) { return _listFuncExtended[i]; }
        const static TestFunctionRealWithIntegral& getFuncExtendedWithIntegral(int i) { return _listFuncExtendedWithIntegral[i]; }
        
        const static TestFunctionRealWithDerivation& getFuncExtended(const std::string &funcName)
        {
            for (int i = 0; i < getNumFuncExtended(); i++)
            {
                if (_listFuncExtended[i]._funcName == funcName)
                    return _listFuncExtended[i];
            }
            throw std::runtime_error("Extended TestFunctionReal " + funcName + " not found!");
        }
        
        const static TestFunctionRealWithIntegral& getFuncExtendedWithIntegral(const std::string &funcName)
        {
            for (int i = 0; i < getNumFuncExtendedWithIntegral(); i++)
            {
                if (_listFuncExtendedWithIntegral[i]._funcName == funcName)
                    return _listFuncExtendedWithIntegral[i];
            }
            throw std::runtime_error("Extended TestFunctionRealWithIntegral " + funcName + " not found!");
        }
        
        //--------------------------------------------------------------------------------------------------------------
        // CATEGORY-BASED ACCESSORS
        //--------------------------------------------------------------------------------------------------------------
        
        /** Get all standard functions */
        static std::vector<const TestFunctionReal*> getStandardFunctions()
        {
            std::vector<const TestFunctionReal*> result;
            for (int i = 0; i < getNumFunc(); i++)
                result.push_back(&_listFuncReal[i]);
            return result;
        }
        
        /** Get all numerically interesting functions (extended set) */
        static std::vector<const TestFunctionRealWithDerivation*> getNumericalAnalysisFunctions()
        {
            std::vector<const TestFunctionRealWithDerivation*> result;
            for (int i = 0; i < getNumFuncExtended(); i++)
                result.push_back(&_listFuncExtended[i]);
            return result;
        }
        
        /** Get challenging integration test functions */
        static std::vector<const TestFunctionRealWithIntegral*> getChallengingIntegrationFunctions()
        {
            std::vector<const TestFunctionRealWithIntegral*> result;
            for (int i = 0; i < getNumFuncExtendedWithIntegral(); i++)
                result.push_back(&_listFuncExtendedWithIntegral[i]);
            return result;
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

        const static inline TestFunctionRealWithDerivation _listFuncRealWithDerivation[] = { 
                {"TestDer1", new CompleteRInterval(), 
                            new OpenInterval(-20.0, 20.0 ),
                            [](Real x) { return sin(x);},  "sin(x)", 
                            [](Real x) { return cos(x);},  "cos(x)"},
                // Example 2: x^3, derivative 3*x^2
                {"TestDer2", new CompleteRInterval(),
                            new OpenInterval(-20.0, 20.0),
                            [](Real x) { return x * x * x; },  "x^3",
                            [](Real x) { return 3 * x * x; },  "3*x^2"
                },
                // Example 3: tanh(x), derivative 1/cosh(x)^2
                {"TestDer3", new CompleteRInterval(),
                            new OpenInterval(-10.0, 10.0),
                            [](Real x) { return tanh(x); },  "tanh(x)",
                            [](Real x) { return 1.0 / (cosh(x) * cosh(x)); },  "1/cosh(x)^2"
                },
                // Gaussian: f(x) = exp(-x^2), f'(x) = -2x * exp(-x^2)
                {"Gaussian", new CompleteRInterval(),
                            new ClosedInterval(-5.0, 5.0),
                            [](Real x) { return exp(-x * x); }, "exp(-x^2)",
                            [](Real x) { return -2.0 * x * exp(-x * x); }, "-2x * exp(-x^2)"
                },
                // Logistic: f(x) = 1 / (1 + exp(-x)), f'(x) = exp(-x) / (1 + exp(-x))^2
                {"Logistic", new CompleteRInterval(),
                            new ClosedInterval(-10.0, 10.0),
                            [](Real x) { return 1.0 / (1.0 + exp(-x)); }, "1 / (1 + exp(-x))",
                            [](Real x) { Real e = exp(-x); return e / ((1.0 + e) * (1.0 + e)); }, "exp(-x) / (1 + exp(-x))^2"
                },
                // Error function: f(x) = erf(x), f'(x) = 2/sqrt(pi) * exp(-x^2)
                {"ErrorFunction", new CompleteRInterval(),
                                  new ClosedInterval(-3.0, 3.0),
                                  [](Real x) { return std::erf(x); }, "erf(x)",
                                  [](Real x) { return 2.0 / sqrt(MML::Constants::PI) * exp(-x * x); }, "2/sqrt(pi) * exp(-x^2)"
                },
                // Reciprocal quadratic: f(x) = 1 / (1 + x^2), f'(x) = -2x / (1 + x^2)^2
                {"ReciprocalQuadratic", new CompleteRInterval(),
                                        new ClosedInterval(-10.0, 10.0),
                                        [](Real x) { return 1.0 / (1.0 + x * x); }, "1 / (1 + x^2)",
                                        [](Real x) { return -2.0 * x / pow(1.0 + x * x, 2); }, "-2x / (1 + x^2)^2"
                },
                // Polynomial with Mixed Powers: f(x) = x^5 - 3x^3 + 2x, f'(x) = 5x^4 - 9x^2 + 2
                {"PolyMixedPowers", new CompleteRInterval(),
                                    new ClosedInterval(-10.0, 10.0),
                                    [](Real x) { return pow(x, 5) - 3.0 * pow(x, 3) + 2.0 * x; }, "x^5 - 3x^3 + 2x",
                                    [](Real x) { return 5.0 * pow(x, 4) - 9.0 * pow(x, 2) + 2.0; }, "5x^4 - 9x^2 + 2"
                },
                // Composite Trigonometric: f(x) = sin(x) * cos(x), f'(x) = cos^2(x) - sin^2(x)
                {"SinCosComposite", new CompleteRInterval(),
                                    new ClosedInterval(-2.0 * MML::Constants::PI, 2.0 * MML::Constants::PI),
                                    [](Real x) { return sin(x) * cos(x); }, "sin(x) * cos(x)",
                                    [](Real x) { return cos(x) * cos(x) - sin(x) * sin(x); }, "cos^2(x) - sin^2(x)"
                },
                // Exponential Times Polynomial: f(x) = x^2 * exp(x), f'(x) = exp(x) * (x^2 + 2x)
                {"ExpPoly", new CompleteRInterval(),
                            new ClosedInterval(-10.0, 10.0),
                            [](Real x) { return x * x * exp(x); }, "x^2 * exp(x)",
                            [](Real x) { return exp(x) * (x * x + 2.0 * x); }, "exp(x) * (x^2 + 2x)"
                },
                // Rational Function: f(x) = (x^2 + 1) / (x^3 + 2), f'(x) = [2x(x^3 + 2) - 3x^2(x^2 + 1)] / (x^3 + 2)^2
                {"RationalFunc", new CompleteRInterval(),
                                new ClosedInterval(-10.0, 10.0),
                                [](Real x) { return (x * x + 1.0) / (pow(x, 3) + 2.0); }, "(x^2 + 1) / (x^3 + 2)",
                                [](Real x) {
                                    Real num = 2.0 * x * (pow(x, 3) + 2.0) - 3.0 * pow(x, 2) * (x * x + 1.0);
                                    Real denom = pow(pow(x, 3) + 2.0, 2);
                                    return num / denom;
                                }, "[2x(x^3 + 2) - 3x^2(x^2 + 1)] / (x^3 + 2)^2"
                }
            };

        const static inline TestFunctionRealWithIntegral _listFuncRealWithIntegral[] = { 
                {"TestInt1", new CompleteRInterval(), 
                            new OpenInterval(0.0, 5.0 ),
                            [](Real x) { return x*x*(x*x-2)*sin(x);},  "x*x*(x*x-2.0)*sin(x)", 
                            [](Real x) { return 4*x*(x*x-7)*sin(x)-((Real) pow(x,4.0)-14*x*x+28)*cos(x);},  "4.0*x*(x*x-7.0)*sin(x)-(pow(x,4.0)-14.0*x*x+28.0)*cos(x)"},                                   
                {"TestInt2", new CompleteRInterval(), 
                            new OpenInterval(-20.0, 20.0 ),
                            [](Real x) { return 1.0 / (4.0 + x*x);},  "1 / (4 + x*x)", 
                            [](Real x) { return 0.5 * atan(x/2);},  "1 / 2 * atan(x/2)"},
                // x^2, integral is x^3 / 3
                {"XSquared", new CompleteRInterval(),
                            new OpenInterval(-10.0, 10.0),
                            [](Real x) { return x * x; }, "x^2",
                            [](Real x) { return x * x * x / 3.0; }, "x^3 / 3"
                },
                // 1/x, integral is ln|x|
                {"Reciprocal", new OpenInterval(-10.0, 0.0), // avoid x=0
                              new OpenInterval(0.1, 10.0),
                              [](Real x) { return 1.0 / x; }, "1/x",
                              [](Real x) { return log(std::abs(x)); }, "ln|x|"
                },
                // Exponential Decay: f(x) = exp(-a*x), integral = -1/a * exp(-a*x)
                {"ExpDecay", new CompleteRInterval(),
                            new ClosedInterval(-10.0, 10.0),
                            [](Real x) { constexpr Real a = 2.0; return exp(-a * x); }, "exp(-2x)",
                            [](Real x) { constexpr Real a = 2.0; return -1.0 / a * exp(-a * x); }, "-1/2 * exp(-2x)"
                },
                // Gaussian: f(x) = exp(-x^2), integral = sqrt(pi)/2 * erf(x)
                {"Gaussian", new CompleteRInterval(),
                            new ClosedInterval(-3.0, 3.0),
                            [](Real x) { return exp(-x * x); }, "exp(-x^2)",
                            [](Real x) { return 0.5 * sqrt(MML::Constants::PI) * std::erf(x); }, "0.5*sqrt(pi)*erf(x)"
                },
                // Logistic: f(x) = 1/(1+exp(-x)), integral = x - ln(1+exp(-x))
                {"Logistic", new CompleteRInterval(),
                            new ClosedInterval(-10.0, 10.0),
                            [](Real x) { return 1.0 / (1.0 + exp(-x)); }, "1/(1+exp(-x))",
                            [](Real x) { return x - log(1.0 + exp(-x)); }, "x - ln(1+exp(-x))"
                },
                // Rational: f(x) = x/(x^2+1), integral = 0.5*ln(x^2+1)
                {"RationalXOverXSquaredPlus1", new CompleteRInterval(),
                                              new ClosedInterval(-10.0, 10.0),
                                              [](Real x) { return x / (x * x + 1.0); }, "x/(x^2+1)",
                                              [](Real x) { return 0.5 * log(x * x + 1.0); }, "0.5*ln(x^2+1)"
                },
                // Arctangent: f(x) = atan(x), integral = x*atan(x) - 0.5*ln(1+x^2)
                {"Arctangent", new CompleteRInterval(),
                              new ClosedInterval(-10.0, 10.0),
                              [](Real x) { return atan(x); }, "atan(x)",
                              [](Real x) { return x * atan(x) - 0.5 * log(1.0 + x * x); }, "x*atan(x) - 0.5*ln(1+x^2)"
                },
                // Sine Squared: f(x) = sin^2(x), integral = x/2 - (1/4) * sin(2x)
                {"SineSquared", new CompleteRInterval(),
                                new ClosedInterval(-10.0, 10.0),
                                [](Real x) { return sin(x) * sin(x); }, "sin^2(x)",
                                [](Real x) { return x / 2.0 - 0.25 * sin(2.0 * x); }, "x/2 - 0.25*sin(2x)"
                },
                // Cosine Cubed: f(x) = cos^3(x), integral = (sin(x)/3) + (sin(3x)/9)
                {"CosineCubed", new CompleteRInterval(),
                                new ClosedInterval(-10.0, 10.0),
                                [](Real x) { return pow(cos(x), 3); }, "cos^3(x)",
                                [](Real x) { return (sin(x) / 3.0) + (sin(3.0 * x) / 9.0); }, "sin(x)/3 + sin(3x)/9"
                },
                // Reciprocal Quadratic: f(x) = 1/(x^2 + a^2), integral = (1/a) * atan(x/a)
                {"ReciprocalQuadratic", new CompleteRInterval(),
                                        new ClosedInterval(-10.0, 10.0),
                                        [](Real x) { constexpr Real a = 2.0; return 1.0 / (x * x + a * a); }, "1/(x^2+4)",
                                        [](Real x) { constexpr Real a = 2.0; return (1.0 / a) * atan(x / a); }, "0.5*atan(x/2)"
                },
                // Exponential Times Polynomial: f(x) = x*exp(x), integral = (x-1)*exp(x)
                {"ExpTimesX", new CompleteRInterval(),
                              new ClosedInterval(-10.0, 10.0),
                              [](Real x) { return x * exp(x); }, "x*exp(x)",
                              [](Real x) { return (x - 1.0) * exp(x); }, "(x-1)*exp(x)"
                },
                // Logarithmic Composite: f(x) = ln(1 + x^2), integral = x*ln(1 + x^2) - 2x + 2*atan(x)
                {"LogComposite", new CompleteRInterval(),
                                new ClosedInterval(-10.0, 10.0),
                                [](Real x) { return log(1.0 + x * x); }, "ln(1 + x^2)",
                                [](Real x) { return x * log(1.0 + x * x) - 2.0 * x + 2.0 * atan(x); }, "x*ln(1 + x^2) - 2x + 2*atan(x)"
                },
                // Hyperbolic Tangent: f(x) = tanh(x), integral = ln(cosh(x))
                {"Tanh", new CompleteRInterval(),
                        new ClosedInterval(-10.0, 10.0),
                        [](Real x) { return tanh(x); }, "tanh(x)",
                        [](Real x) { return log(cosh(x)); }, "ln(cosh(x))"
                },
                // Exponential of Quadratic: f(x) = exp(a*x^2), integral = sqrt(pi)/(2*sqrt(a)) * erf(sqrt(a)*x), a > 0
                {"ExpQuadratic", new CompleteRInterval(),
                                new ClosedInterval(-3.0, 3.0),
                                [](Real x) { constexpr Real a = 1.0; return exp(a * x * x); }, "exp(x^2)",
                                [](Real x) { constexpr Real a = 1.0; return 0.5 * sqrt(MML::Constants::PI / a) * std::erf(sqrt(a) * x); }, "0.5*sqrt(pi)*erf(x)"
                }
            };

        //==============================================================================================================
        // EXTENDED FUNCTIONS - Numerically Interesting Functions
        // 
        // These functions are specifically chosen for their challenging numerical properties:
        // - Runge function: Classic example of polynomial interpolation failure
        // - Steep gradients: Test adaptive algorithms
        // - Highly oscillatory: Challenge integration routines
        // - Peaked functions: Test resolution capabilities
        // - Near-singular behavior: Test boundary handling
        //==============================================================================================================
        
        const static inline TestFunctionRealWithDerivation _listFuncExtended[] = {
            // 1. Runge Function: f(x) = 1/(1 + 25x²)
            // Classic example showing Runge's phenomenon - polynomial interpolation fails at edges
            // f'(x) = -50x / (1 + 25x²)²
            {"Runge", new CompleteRInterval(),
                      new ClosedInterval(-1.0, 1.0),
                      [](Real x) { return 1.0 / (1.0 + 25.0 * x * x); }, "1/(1 + 25x^2)",
                      [](Real x) { Real d = 1.0 + 25.0 * x * x; return -50.0 * x / (d * d); }, "-50x/(1 + 25x^2)^2"
            },
            
            // 2. Steep Gradient: f(x) = tanh(kx) with large k
            // Tests adaptive methods with rapid transitions
            // f'(x) = k / cosh²(kx)
            {"SteepGradient_k10", new CompleteRInterval(),
                                  new ClosedInterval(-2.0, 2.0),
                                  [](Real x) { constexpr Real k = 10.0; return tanh(k * x); }, "tanh(10x)",
                                  [](Real x) { constexpr Real k = 10.0; Real c = cosh(k * x); return k / (c * c); }, "10/cosh^2(10x)"
            },
            
            // 3. Peaked Gaussian: f(x) = exp(-x²/ε) for small ε
            // Narrow peak challenges quadrature rules
            // f'(x) = -2x/ε * exp(-x²/ε)
            {"PeakedGaussian_eps0.01", new CompleteRInterval(),
                                       new ClosedInterval(-1.0, 1.0),
                                       [](Real x) { constexpr Real eps = 0.01; return exp(-x * x / eps); }, "exp(-x^2/0.01)",
                                       [](Real x) { constexpr Real eps = 0.01; return -2.0 * x / eps * exp(-x * x / eps); }, "-200x*exp(-x^2/0.01)"
            },
            
            // 4. Witch of Agnesi: f(x) = 1/(1 + x²)
            // Classical benchmark function, well-behaved but interesting
            // f'(x) = -2x / (1 + x²)²
            {"WitchOfAgnesi", new CompleteRInterval(),
                             new ClosedInterval(-10.0, 10.0),
                             [](Real x) { return 1.0 / (1.0 + x * x); }, "1/(1 + x^2)",
                             [](Real x) { Real d = 1.0 + x * x; return -2.0 * x / (d * d); }, "-2x/(1 + x^2)^2"
            },
            
            // 5. Oscillatory Decay: f(x) = exp(-x) * sin(10x)
            // Combines decay with oscillation - challenges many algorithms
            // f'(x) = exp(-x) * (10*cos(10x) - sin(10x))
            {"OscillatoryDecay", new CompleteRInterval(),
                                 new ClosedInterval(0.0, 10.0),
                                 [](Real x) { return exp(-x) * sin(10.0 * x); }, "exp(-x)*sin(10x)",
                                 [](Real x) { return exp(-x) * (10.0 * cos(10.0 * x) - sin(10.0 * x)); }, "exp(-x)*(10*cos(10x) - sin(10x))"
            },
            
            // 6. Near-Singular Power: f(x) = x^0.1 for x > 0
            // Tests handling of functions with infinite derivative at boundary
            // f'(x) = 0.1 * x^(-0.9)
            {"NearSingularPower", new OpenToInfInterval(0.0),
                                  new OpenClosedInterval(0.001, 10.0),  // Avoid x=0
                                  [](Real x) { return pow(x, 0.1); }, "x^0.1",
                                  [](Real x) { return 0.1 * pow(x, -0.9); }, "0.1*x^(-0.9)"
            },
            
            // 7. Bessel-like oscillation: f(x) = sin(x)/x (sinc function)
            // Removable singularity at x=0, oscillatory decay
            // f'(x) = (x*cos(x) - sin(x))/x²
            {"Sinc", new CompleteRInterval(),
                     new ClosedInterval(0.1, 20.0),  // Avoid singularity at 0 for testing
                     [](Real x) { return x != 0.0 ? sin(x) / x : 1.0; }, "sin(x)/x",
                     [](Real x) { return x != 0.0 ? (x * cos(x) - sin(x)) / (x * x) : 0.0; }, "(x*cos(x) - sin(x))/x^2"
            },
            
            // 8. Bump-like function: f(x) = exp(-1/(1-x²)) for |x| < 1
            // Smooth function with very flat behavior near boundaries
            // Used in smooth approximations
            {"SmoothBump", new OpenInterval(-1.0, 1.0),
                           new OpenInterval(-0.99, 0.99),
                           [](Real x) { 
                               Real t = 1.0 - x * x;
                               return t > 0 ? exp(-1.0 / t) : 0.0; 
                           }, "exp(-1/(1-x^2))",
                           [](Real x) { 
                               Real t = 1.0 - x * x;
                               if (t <= 0) return 0.0;
                               return -2.0 * x * exp(-1.0 / t) / (t * t); 
                           }, "-2x*exp(-1/(1-x^2))/(1-x^2)^2"
            }
        };

        //==============================================================================================================
        // EXTENDED INTEGRATION FUNCTIONS - Challenging Integrals
        // 
        // These are particularly challenging for numerical integration:
        // - Fresnel-type integrands
        // - Peaked integrands
        // - Oscillatory integrands
        //==============================================================================================================
        
        const static inline TestFunctionRealWithIntegral _listFuncExtendedWithIntegral[] = {
            // 1. Fresnel-like: f(x) = sin(x²)
            // Integral is Fresnel S function: sqrt(π/2) * S(sqrt(2/π) * x)
            // Highly oscillatory, challenging for standard quadrature
            {"FresnelSine", new CompleteRInterval(),
                            new ClosedInterval(0.0, 5.0),
                            [](Real x) { return sin(x * x); }, "sin(x^2)",
                            [](Real x) { 
                                // Fresnel S integral approximation using series or erf
                                // For simplicity, we use a reference computation
                                Real u = x * sqrt(2.0 / MML::Constants::PI);
                                // S(u) approximation via series for small u
                                Real sum = 0.0;
                                Real term = u;
                                for (int n = 0; n < 50 && std::abs(term) > 1e-15; n++) {
                                    sum += term;
                                    term *= -MML::Constants::PI * MML::Constants::PI * u * u * u * u / 
                                            ((4.0 * n + 3.0) * (4.0 * n + 4.0) * (2.0 * n + 1.0) * (2.0 * n + 2.0));
                                }
                                return sqrt(MML::Constants::PI / 2.0) * sum;
                            }, "sqrt(pi/2)*FresnelS(sqrt(2/pi)*x)"
            },
            
            // 2. Fresnel-like: f(x) = cos(x²)
            // Integral is Fresnel C function
            {"FresnelCosine", new CompleteRInterval(),
                              new ClosedInterval(0.0, 5.0),
                              [](Real x) { return cos(x * x); }, "cos(x^2)",
                              [](Real x) { 
                                  Real u = x * sqrt(2.0 / MML::Constants::PI);
                                  Real sum = 0.0;
                                  Real term = 1.0;
                                  for (int n = 0; n < 50 && std::abs(term) > 1e-15; n++) {
                                      sum += term;
                                      term *= -MML::Constants::PI * MML::Constants::PI * u * u * u * u / 
                                              ((4.0 * n + 1.0) * (4.0 * n + 2.0) * (2.0 * n + 1.0) * (2.0 * n + 2.0));
                                  }
                                  return sqrt(MML::Constants::PI / 2.0) * u * sum;
                              }, "sqrt(pi/2)*FresnelC(sqrt(2/pi)*x)"
            },
            
            // 3. Runge function integral
            // ∫ 1/(1+25x²) dx = (1/5)*atan(5x)
            {"RungeIntegral", new CompleteRInterval(),
                              new ClosedInterval(-1.0, 1.0),
                              [](Real x) { return 1.0 / (1.0 + 25.0 * x * x); }, "1/(1+25x^2)",
                              [](Real x) { return 0.2 * atan(5.0 * x); }, "(1/5)*atan(5x)"
            },
            
            // 4. Peaked Gaussian integral
            // ∫ exp(-x²/ε) dx = sqrt(π*ε) * erf(x/sqrt(ε)) / 2
            {"PeakedGaussianIntegral", new CompleteRInterval(),
                                       new ClosedInterval(-1.0, 1.0),
                                       [](Real x) { constexpr Real eps = 0.01; return exp(-x * x / eps); }, "exp(-x^2/0.01)",
                                       [](Real x) { constexpr Real eps = 0.01; return 0.5 * sqrt(MML::Constants::PI * eps) * std::erf(x / sqrt(eps)); }, "0.5*sqrt(0.01*pi)*erf(10x)"
            },
            
            // 5. Oscillatory decay integral
            // ∫ exp(-x)*sin(ωx) dx = exp(-x)*(-sin(ωx) - ω*cos(ωx))/(1+ω²)
            {"OscillatoryDecayIntegral", new CompleteRInterval(),
                                         new ClosedInterval(0.0, 10.0),
                                         [](Real x) { constexpr Real w = 10.0; return exp(-x) * sin(w * x); }, "exp(-x)*sin(10x)",
                                         [](Real x) { 
                                             constexpr Real w = 10.0; 
                                             return exp(-x) * (-sin(w * x) - w * cos(w * x)) / (1.0 + w * w); 
                                         }, "exp(-x)*(-sin(10x) - 10*cos(10x))/101"
            },
            
            // 6. Near-singular power integral
            // ∫ x^α dx = x^(α+1)/(α+1) for α > -1
            {"NearSingularPowerIntegral", new OpenToInfInterval(0.0),
                                          new OpenClosedInterval(0.001, 10.0),
                                          [](Real x) { return pow(x, 0.1); }, "x^0.1",
                                          [](Real x) { return pow(x, 1.1) / 1.1; }, "x^1.1/1.1"
            }
        };
    };
}
#endif
