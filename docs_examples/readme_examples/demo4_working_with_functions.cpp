#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Function.h"
#include "core/FunctionHelpers.h"
#include "core/Derivation.h"
#include "core/Integration.h"
#endif


using namespace MML;

void Readme_deriving_functions()
{
    RealFunction       f1{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };

    // numerical derivation of real function (available orders - 1, 2, 4, 6, 8)
    double der_f1 = Derivation::NDer1(f1, 0.5);
    double der_f4 = Derivation::NDer2(f1, 0.5, 1e-6);   // setting explicit step size
    double err;
    double der_f6 = Derivation::NDer6(f1, 0.5, &err);   // if we need error estimate    
    // we can use default Derive routine (set to NDer4), but it requires error estimate
    double num_der4 = Derivation::Derive(f1, 0.5, nullptr);

    // second and third derivatives
    double sec_der_f1   = Derivation::NSecDer2(f1, 0.5);
    double third_der_f1 = Derivation::NThirdDer2(f1, 0.5);

    // creating new function that is derivation of existing function
    RealFuncDerived4    f1_der4(f1);        // 4th order derivation

    // scalar and vector functions
    ScalarFunction<3>   f2Scal([](const VectorN<Real, 3> &x) { return 1.0 / pow(x.NormL2(), 2); });
    VectorFunction<3>   f3Vec([](const VectorN<Real, 3> &x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });
    VectorN<Real, 3>    der_point{1.0, 1.0, 1.0};

    double der_f2               = Derivation::NDer1Partial(f2Scal, 1, der_point);
    VectorN<Real, 3> der_f2_all = Derivation::NDer1PartialByAll(f2Scal, der_point);

    double der_f3 = Derivation::NDer1Partial(f3Vec, 1, 1, der_point);
    VectorN<Real, 3>     der_f3_by1    = Derivation::NDer2PartialByAll(f3Vec, 1, der_point);
    MatrixNM<Real, 3, 3> der_f3_by_all = Derivation::NDer4PartialAllByAll(f3Vec, der_point);

    // TODO 0.8 - s Function Analyzerom vidjeti koliko je dobra derivacija
}

void Readme_integrating_functions()
{
    RealFunction f1{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };

    double a = 0.0;
    double b = 1.0;
    double int_trap = IntegrateTrap(f1,a,b);
    double int_simp = IntegrateSimpson(f1,a,b);
    double int_romb = IntegrateRomberg(f1,a,b);
    // we can use default Integrate routine (set to IntegrateSimpson), requires precision
    double int_def = Integrate(f1, a, b, 1e-04);

    // 2D integration of constant scalar 2D function (ie. we'll get the area of the surface)
    ScalarFunction<2> f2([](const VectorN<Real, 2> &x) { return 1.0; });
    
    // we integrate over circle with radius 2
    Real val = IntegrateSurface(f2, IntegrationMethod::GAUSS10, 
                                    -2, 2,              // x range
                                    [](Real x) { return -sqrt(4 - x*x);},   // y range lower limit
                                    [](Real x) { return sqrt(4 - x*x);});   // y range upper limit

    std::cout << "Calc. area = " << val << ", exact value: 4 * PI = " << 4 * Constants::PI << std::endl;

    // 3D integration of constant scalar 3D function (ie. we'll get the volume of the solid)
    ScalarFunction<3> f3([](const VectorN<Real, 3> &x) { return 1.0; });
    
    // integration over sphere of radius 1
    Real vol = IntegrateVolume( f3, 
                                -1, 1,              
                                [](Real x) { return -sqrt(1 - x*x);}, 
                                [](Real x) { return sqrt(1 - x*x);}, 
                                [](Real x, Real y) { return -sqrt(1 - x*x - y*y);}, 
                                [](Real x, Real y) { return sqrt(1 - x*x - y*y);});

    std::cout << "Calc. vol. = " << vol << ", exact value: 4/3 * PI = " << 4.0/3.0 * Constants::PI << std::endl;

/* OUTPUT 
    Calc. area = 12.57211164, exact value: 4 * PI = 12.56637061
    Calc. vol. = 4.190703882, exact value: 4/3 * PI = 4.188790205
*/
}

void Readme_working_with_functions()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                README - working with functions                ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Readme_deriving_functions();
    Readme_integrating_functions();
}