#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Function.h"
#include "core/InterpolatedFunction.h"

#include "core/Integration.h"
#endif

using namespace std;
using namespace MML;


//////////////////////          REAL FUNCTION DERIVATION           ////////////////////////
double DemoIntRealFunc_TestFunc(double x) 
{ 
    return sin(x)*(1.0 + 0.5*x*x); 
}

void Demo_Integration_func_ptr()
{
    // creating a function object from a already existing (standalone) function
    RealFunction f1(DemoIntRealFunc_TestFunc);

    // or creating a function object directly
    RealFunction f2{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };

    // creating a function object from a std::function object
    std::function<double(double)> h{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };
    RealFunctionFromStdFunc f3(h);

    double a = 0.0;
    double b = 1.0;
    double int_trap = Integration::IntegrateTrap(f1,a,b);
    double int_simp = Integration::IntegrateSimpson(f1,a,b);
    double int_romb = Integration::IntegrateRomberg(f1,a,b);

    // we can use default Integrate routine (default set to IntegrateSimpson)
    double int_def = Integration::Integrate(f1, a, b, 1e-04);
}

// If you CAN change the class where your data for calculation is
class ClassProvidingFuncToIntegrate
{
    private:
        double _param;
    public:
        ClassProvidingFuncToIntegrate(double param) : _param(param) { }
    
        // just provide operator() for your class
        double operator()(double x ) { return _param * sin(x); }
};

void Demo_Integration_member_fun()
{
    ClassProvidingFuncToIntegrate   funcObj(3.0);
    
    RealFunctionFromStdFunc g(std::function<double(double)>{funcObj});

    double a = 0.0;
    double b = 1.0;
    double int_trap = Integration::IntegrateTrap(g,a,b);
    double int_simp = Integration::IntegrateSimpson(g,a,b);
    double int_romb = Integration::IntegrateRomberg(g,a,b);
}

// If you CAN'T change the class where your data for calculation is
class BigComplexClassYouCantChangeInt
{
    // has data for calculating function you want to derive
};

// Create a helper wrapper class that will provide operator() for your class
// and implement calculation of function value
class BigComplexIntegrateFunc
{
    const BigComplexClassYouCantChangeInt &_ref;
public:
    BigComplexIntegrateFunc(const BigComplexClassYouCantChangeInt &bigClass) : _ref(bigClass) { }

    double operator()(double x ) 
    {
        double val = 0.0; 
        // complex calculations using _ref
        return val; 
    }
};

void Demo_Integration_member_fun2(const BigComplexClassYouCantChangeInt &ref)
{
    BigComplexIntegrateFunc funcObj(ref);  
    RealFunctionFromStdFunc func_to_integrate(std::function<double(double)>{funcObj});
    
    double a = 0.0;
    double b = 1.0;
    double int_trap = Integration::IntegrateTrap(func_to_integrate,a,b);
    double int_simp = Integration::IntegrateSimpson(func_to_integrate,a,b);
    double int_romb = Integration::IntegrateRomberg(func_to_integrate,a,b);
}

void Demo_Integration_Interpolated_RealFunc()
{
    // creating interpolated function
    Vector<double> x_val(100);
    Vector<double> y_val(100);
    for( int i=0; i<100; i++ )
    {
        x_val[i] = i / 100.0;
        y_val[i] = sin(x_val[i])*(1.0 + 0.5*x_val[i]*x_val[i]);
    }

    LinearInterpRealFunc f_linear(x_val, y_val);

    double a = 0.0;
    double b = 1.0;
    double int_trap = Integration::IntegrateTrap(f_linear,a,b);
    double int_simp = Integration::IntegrateSimpson(f_linear,a,b);
    double int_romb = Integration::IntegrateRomberg(f_linear,a,b);
}

void Demo_Integration()
{
    std::cout << endl;
    std::cout << "***********************************************************************" << endl;
    std::cout << "****                         INTEGRATION                           ****" << endl;
    std::cout << "***********************************************************************" << endl;

    Demo_Integration_func_ptr();
    Demo_Integration_member_fun();
    Demo_Integration_member_fun2(BigComplexClassYouCantChangeInt{});
    Demo_Integration_Interpolated_RealFunc();
}