#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/Vector.h"
#include "algorithms/DiffEqSolvers.h"
#endif

void VanDerPol(double eps, const double x, MML::Vector<Real> &y, MML::Vector<Real> &dydx) {
    dydx[0]= y[1];
    dydx[1]=((1.0-y[0]*y[0])*y[1]-y[0])/eps;
}

void Demo_DiffEqSolvers()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                    DIFF.EQUATIONS SOLVERS                     ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
}