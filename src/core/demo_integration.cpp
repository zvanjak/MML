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


//////////////////////          REAL FUNCTION INTEGRATION           ////////////////////////
Real DemoIntRealFunc_TestFunc(Real x) 
{ 
    return sin(x)*(1 + x*x / 2); 
}

void Demo_Integration_func_ptr()
{
    // creating a function object from a already existing (standalone) function
    RealFunction f1(DemoIntRealFunc_TestFunc);

    // or creating a function object directly
    RealFunction f2{[](Real x) { return sin(x)*(1 + x*x / 2); } };

    // creating a function object from a std::function object
    std::function<Real(Real)> h{[](Real x) { return sin(x)*(1 + x*x / 2); } };
    RealFunctionFromStdFunc f3(h);

    Real a = 0.0;
    Real b = 1.0;
    Real int_trap = IntegrateTrap(f1,a,b);
    Real int_simp = IntegrateSimpson(f1,a,b);
    Real int_romb = IntegrateRomberg(f1,a,b);

    // we can use default Integrate routine (default set to IntegrateSimpson)
    Real int_def = Integrate(f1, a, b, 1e-04);
}

// If you CAN change the class where your data for calculation is
class ClassProvidingFuncToIntegrate
{
    private:
        Real _param;
    public:
        ClassProvidingFuncToIntegrate(Real param) : _param(param) { }
    
        // just provide operator() for your class
        Real operator()(Real x ) { return _param * sin(x); }
};

void Demo_Integration_member_fun()
{
    ClassProvidingFuncToIntegrate   funcObj(3.0);
    
    RealFunctionFromStdFunc g(std::function<double(double)>{funcObj});

    Real a = 0.0;
    Real b = 1.0;
    Real int_trap = IntegrateTrap(g,a,b);
    Real int_simp = IntegrateSimpson(g,a,b);
    Real int_romb = IntegrateRomberg(g,a,b);
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

    Real operator()(Real x ) 
    {
        Real val = 0.0; 
        // complex calculations using _ref
        return val; 
    }
};

void Demo_Integration_member_fun2(const BigComplexClassYouCantChangeInt &ref)
{
    BigComplexIntegrateFunc funcObj(ref);  
    RealFunctionFromStdFunc func_to_integrate(std::function<double(double)>{funcObj});
    
    Real a = 0.0;
    Real b = 1.0;
    Real int_trap = IntegrateTrap(func_to_integrate,a,b);
    Real int_simp = IntegrateSimpson(func_to_integrate,a,b);
    Real int_romb = IntegrateRomberg(func_to_integrate,a,b);
}

void Demo_Integration_Interpolated_RealFunc()
{
    // creating interpolated function
    Vector<Real> x_val(100);
    Vector<Real> y_val(100);
    for( int i=0; i<100; i++ )
    {
        x_val[i] = i / 100.0;
        y_val[i] = sin(x_val[i])*(1.0 + 0.5*x_val[i]*x_val[i]);
    }

    LinearInterpRealFunc f_linear(x_val, y_val);

    Real a = 0.0;
    Real b = 1.0;
    Real int_trap = IntegrateTrap(f_linear,a,b);
    Real int_simp = IntegrateSimpson(f_linear,a,b);
    Real int_romb = IntegrateRomberg(f_linear,a,b);
}

void Demo_Integration2D()
{
    ScalarFunction<2> f([](const VectorN<Real, 2> &x) { return Real{1.0}; });
    
    Real integral     = IntegrateSurface(f, IntegrationMethod::GAUSS10, 0, 2, [](Real x) { return Real{1.0};}, [](Real x) { return Real{4.0};}); 
    Real integralCircle = IntegrateSurface(f, IntegrationMethod::TRAP, -2, 2, [](Real x) { return -sqrt(4 - x*x);}, [](Real x) { return sqrt(4 - x*x);}); 

    std::cout << "Integral = " << integral << std::endl;
    std::cout << "Integral circle = " << integralCircle << "4 * PI = " << 4 * Constants::PI << std::endl;
}

void Demo_Integration3D()
{
    ScalarFunction<3> f([](const VectorN<Real, 3> &x) { return Real{1.0}; });
    
    Real integral = IntegrateVolume(f, 0, 1, [](Real x) { return Real{1.0};}, [](Real x) { return Real{4.0};}, [](Real x, Real y) { return Real{1.0};}, [](Real x, Real y) { return Real{5.0};});

    std::cout << "Integral = " << integral << std::endl;
}

void Demo_Integration()
{
    std::cout << endl;
    std::cout << "***********************************************************************" << endl;
    std::cout << "****                         INTEGRATION                           ****" << endl;
    std::cout << "***********************************************************************" << endl;

    // Demo_Integration_func_ptr();
    // Demo_Integration_member_fun();
    // Demo_Integration_member_fun2(BigComplexClassYouCantChangeInt{});
    // Demo_Integration_Interpolated_RealFunc();

    Demo_Integration2D();
    Demo_Integration3D();
}