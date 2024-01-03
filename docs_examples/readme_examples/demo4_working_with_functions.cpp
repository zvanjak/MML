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
    RealFunction f1{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };

    // numerical derivation of real function (different orders)
    double der_f1 = Derivation::NDer1(f1, 0.5);
    double der_f2 = Derivation::NDer2(f1, 0.5);
    double der_f4 = Derivation::NDer4(f1, 0.5, 1e-6);   // setting explicit step size
    double der_f6 = Derivation::NDer6(f1, 0.5);    
    double der_f8 = Derivation::NDer8(f1, 0.5);
    // we can use default Derive routine (default set to NDer4)
    double num_der4 = Derivation::Derive(f1, 0.5, nullptr);
}

void Readme_integrating_functions()
{
    RealFunction f1{[](double x) { return sin(x)*(1.0 + 0.5*x*x); } };
        
    double a = 0.0;
    double b = 1.0;
    double int_trap = Integration::IntegrateTrap(f1,a,b);
    double int_simp = Integration::IntegrateSimpson(f1,a,b);
    double int_romb = Integration::IntegrateRomberg(f1,a,b);
    // we can use default Integrate routine (default set to IntegrateSimpson)
    double int_def = Integration::Integrate(f1, a, b, 1e-04);  
}

void Readme_interpolating_functions()
{

}

void Readme_working_with_functions()
{

}