#include <iostream>
#include <iomanip>
#include <cmath>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/Vector.h"
#endif

void Demo_Vector()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                            VECTOR                             ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    
    std::vector<double> vec{1.0, 2.0, 3.0, 4.0, 5.0};
    double  arr[5] = {0.0, 0.0, 0.0, 10.0, 0.0};

    MML::Vector   a(5);
    MML::Vector   b(5, 3.14159);
    MML::Vector   c(vec);
    MML::Vector   d({3.0, 3.0, 3.0});
    MML::Vector   e(5, arr);

    MML::Vector f = a + b + c;

    std::cout << "a = " << a << std::endl;
    std::cout << "b = " << b << std::endl;
    std::cout << "c = " << c << std::endl;
    std::cout << "d = " << d << std::endl;
    std::cout << "e = " << e << std::endl;
    std::cout << "f = " << f << std::endl;

    std::cout << "a + b = " << a + b << std::endl;
    std::cout << "a - b = " << a - b << std::endl;
    std::cout << "2.0 * c = " << 2.0 * c << std::endl;
    std::cout << "c * 2.0 = " << c * 2.0 << std::endl;
    std::cout << "c / 2.0 = " << c / 2.0 << std::endl;

    std::cout << "NormCartesian(a) = " << a.NormL2() << std::endl;
    std::cout << "b.ScalarProductCartesian(c) = " << b.ScalarProductCartesian(c) << std::endl;

}