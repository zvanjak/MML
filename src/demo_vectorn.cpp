#include <iostream>
#include <iomanip>
#include <cmath>

#ifdef MML_USE_SINGLE_HEADER
#include "MMLBasicTypes.h"
#else
#include "basic_types/VectorN.h"
#endif

void Demo_VectorN()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                           VECTOR N                            ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    MML::VectorN<3>   a({1.0, 1.0, 1.0});
    MML::VectorN<3>   b({2.0, 2.0, 2.0});
    MML::VectorN<3>   c({3.0, 3.0, 3.0});

    MML::VectorN<3> d = a + b + c;

    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "c = " << c << std::endl;

    std::cout << "a + b = " << a + b << std::endl;
    std::cout << "a - b = " << a - b << std::endl;
    std::cout << "2.0 * c = " << 2.0 * c << std::endl;
    std::cout << "c * 2.0 = " << c * 2.0 << std::endl;
    std::cout << "c / 2.0 = " << c / 2.0 << std::endl;

    std::cout << "NormCartesian(a) = " << a.NormL2() << std::endl;

    MML::Vector3Cartesian c1(b), c2(c);

    std::cout << "b.ScalarProductCartesian(c) = " << ScalarProd(c1, c2) << std::endl;   
}