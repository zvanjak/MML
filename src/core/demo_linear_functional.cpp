// TODO - BIG, EMPTY!!! - Demo functionals
#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"    

#include "base/LinearFunctional.h"
#endif

using namespace MML;

// TODO - BIG!
void Demo_LinearFunctionals()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                    LINEAR FUNCTIONALS                         ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    LinearFunctionalN<3> f({1, 2, 3});

    VectorN<Real, 3> v1({1, 2, 3});
    Real r = f(v1);

    LinearFunctionalN<2, Complex> f_complex_2({Complex(1,1), Complex(1,1)});
    VectorN<Complex, 2> v2({Complex(1,1), Complex(1,1)});
    Complex c = f_complex_2(v2);
}
