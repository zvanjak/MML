#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Function.h"
#include "core/Derivation.h"
#include "core/Integration.h"
#endif


using namespace MML;

void Readme_deriving_functions()
{
    RealFunction       f1{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };
    ScalarFunction<3>  f2Scal([](const VectorN<Real, 3> &x) { return 1.0 / pow(x.NormL2(), 2); });
    VectorFunction<3>  f3Vec([](const VectorN<Real, 3> &x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });

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

    VectorN<Real, 3> der_point{1.0, 1.0, 1.0};

    double der_f2               = Derivation::NDer1Partial(f2Scal, 1, der_point);
    VectorN<Real, 3> der_f2_all = Derivation::NDer1PartialByAll(f2Scal, der_point);

    double der_f3 = Derivation::NDer1Partial(f3Vec, 1, 1, der_point);
    VectorN<Real, 3>     der_f3_by1    = Derivation::NDer2PartialByAll(f3Vec, 1, der_point);
    MatrixNM<Real, 3, 3> der_f3_by_all = Derivation::NDer4PartialAllByAll(f3Vec, der_point);
}

void Readme_integrating_functions()
{
    RealFunction f1{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };

    double a = 0.0;
    double b = 1.0;
    double int_trap = Integration::IntegrateTrap(f1,a,b);
    double int_simp = Integration::IntegrateSimpson(f1,a,b);
    double int_romb = Integration::IntegrateRomberg(f1,a,b);
    // we can use default Integrate routine (set to IntegrateSimpson), requires precision
    double int_def = Integration::Integrate(f1, a, b, 1e-04);  
}

void Readme_interpolating_functions()
{

}

void Readme_working_with_functions()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                README - working with functions                ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    Readme_deriving_functions();
    Readme_integrating_functions();
    Readme_interpolating_functions();
}