#include <iostream>
#include <iomanip>
#include <cmath>

#ifdef MML_USE_SINGLE_HEADER
#include "MMLBasicTypes.h"
#else
#endif

void Demo_Vector()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                            VECTOR                             ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    
    MML::Vector   a({1.0, 1.0, 1.0});
    MML::Vector   b({2.0, 2.0, 2.0});
    MML::Vector   c({3.0, 3.0, 3.0});

    MML::Vector d = a + b + c;

    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "c = " << c << std::endl;

    std::cout << "a + b = " << a + b << std::endl;
    std::cout << "a - b = " << a - b << std::endl;
    std::cout << "2.0 * c = " << 2.0 * c << std::endl;
    std::cout << "c * 2.0 = " << c * 2.0 << std::endl;
    std::cout << "c / 2.0 = " << c / 2.0 << std::endl;

    std::cout << "NormCartesian(a) = " << a.NormCartesian() << std::endl;
    std::cout << "b.ScalarProductCartesian(c) = " << b.ScalarProductCartesian(c) << std::endl;

}