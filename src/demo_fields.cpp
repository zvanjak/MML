#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/VectorN.h"
#include "basic_types/Function.h"
#include "basic_types/Fields.h"
#endif

using namespace MML;

double GravityPotential(double r)
{
    return 1 / r;
}

void Demo_Fields()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                          FIELDS                               ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

/*
    - definirati različite potencijalne funkcije (ScalarFunction)
        - gravitacija
    - definirati vektorske funcije polja
        - polje EM naboja u pokretu
        - polje vodiča kroz koji teče struja
*/
    InverseRadialFieldCart gravityPotentialField(10);
    
    RealFunction func2([](double x) { return x; });
    ScalarFunction<3> funcScalar([](const VectorN<Real, 3> &x) { return x[0]; });
    VectorFunction<3> funcVector([](const VectorN<Real, 3> &x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });

    ParametricCurve<3> paramCurve([](double x) { return VectorN<Real, 3>{x, 2 * x, 3 * x}; });
}