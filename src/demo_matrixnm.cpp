#include <iostream>
#include <iomanip>
#include <cmath>

#ifdef MML_USE_SINGLE_HEADER
#include "MMLBasicTypes.h"
#else
#include "basic_types/MatrixNM.h"
#endif

void Basic_MatrixNM_operations() 
{
    MML::MatrixNM<2,2> m1({1.0, -1.0, 1.5, 3.0}), m2;
    m2.MakeUnitMatrix();

    std::cout << "m1 = " << m1 << std::endl;
    std::cout << "m2 = " << m2 << std::endl;

    std::cout << "m1 + m2 = " << m1 + m2  << std::endl;
    std::cout << "m1 - m2 = " << m1 - m2 << std::endl;
    std::cout << "m1 * m2 = " << m1 * m2 << std::endl;

    MML::MatrixNM<1,3> m3{1.0, 1.0, 1.0};
    MML::MatrixNM<3,4> m4{1.0, 0.0, 0.0, 0.0,
                            0.0, 1.0, 0.0, 0.0, 
                            0.0, 0.0, 1.0, 1.0};
    MML::MatrixNM<1,4> m5  = m3 * m4;

    std::cout << "m3 = " << m3 << std::endl;
    std::cout << "m4 = " << m4 << std::endl;
    std::cout << "m3 * m4 = " << m5  << std::endl;

}

void Demo_MatrixNM()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                         MATRIX N_M                            ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    
    MML::MatrixNM<2,2> b;
    MML::MatrixNM<2,2> c({1.0, 0.0, 0.0, 1.0});
    MML::MatrixNM<2,2> d(c);
    MML::MatrixNM<2,2> e = c;

    std::cout << "b = " << b << std::endl;
    std::cout << "c = " << c << std::endl;
    std::cout << "d = " << d << std::endl;
    std::cout << "e = " << e << std::endl;

    Basic_MatrixNM_operations();
}