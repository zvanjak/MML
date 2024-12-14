#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/VectorN.h"
#endif

using namespace MML;

void Demo_VectorN()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                           VECTOR N                            ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    std::vector<double> std_vec{-1.0, 5.0, -2.0};
    float  arr[3] = {-1.0, 5.0, -2.0};

    VectorN<double, 5> vecN_dbl_1;                       // init vector with 5 elements
    Vec3Dbl            vecN_dbl_2(3.14159);              // init with constant value
    Vec3D              vecN_dbl_3({ 1.0, 2.0, 3.0 });    // init with list of values
    VectorN<double, 5> vecN_dbl_4(vecN_dbl_1);           // init with copy ctor
    VectorN<double, 3> vecN_dbl_5 = vecN_dbl_2;          // init with assignment
    VectorN<double, 3> vecN_dbl_6(std_vec);              // init with std::vector<>
    Vec3Flt            vecN_flt_1(arr);                  // init with C/C++ array
    VectorN<double, 5> vecN_unit(3);                     // unit vector

    VectorN<Complex, 3>    vecN_cmplx_1({ 1.0, 2.0, 3.0 });     // init with list of double values
    Vec3C                  vecN_cmplx_2({ Complex(1,1), 
                                          Complex(-1,2), 
                                          Complex(2, -0.5) });  // init with list of complex values
    std::cout << "\nVectorN - declarations and initializations:" << std::endl;

    std::cout << "std::vector<double> std_vec{-1.0, 5.0, -2.0};" << std::endl;
    std::cout << "float  arr[3] = {-1.0, 5.0, -2.0};" << std::endl;

    std::cout << "VectorN<double, 5> vecN_dbl_1;                    vecN_dbl_1   = " << vecN_dbl_1 << std::endl;
    std::cout << "Vec3Dbl            vecN_dbl_2(3.14159);           vecN_dbl_2   = " << vecN_dbl_2 << std::endl;
    std::cout << "Vec3D              vecN_dbl_3({ 1.0, 2.0, 3.0 }); vecN_dbl_3   = " << vecN_dbl_3 << std::endl;
    std::cout << "VectorN<double, 5> vecN_dbl_4(vecN_dbl_1);        vecN_dbl_4   = " << vecN_dbl_4 << std::endl;
    std::cout << "VectorN<double, 3> vecN_dbl_5 = vecN_dbl_2;       vecN_dbl_5   = " << vecN_dbl_5 << std::endl;
    std::cout << "VectorN<double, 3> vecN_dbl_6(std_vec);           vecN_dbl_6   = " << vecN_dbl_6 << std::endl;
    std::cout << "Vec3Flt            vecN_flt_1(arr);               vecN_flt_1   = " << vecN_flt_1 << std::endl;
    std::cout << "VectorN<double, 5> vecN_unit(3);                  vecN_unit    = " << vecN_unit << std::endl;
    std::cout << "vecN_cmplx_1 = " << vecN_cmplx_1 << std::endl;
    std::cout << "vecN_cmplx_2 = " << vecN_cmplx_2 << std::endl;

    std::cout << "\nVectorN - arithmetic operations:" << std::endl;
    std::cout << "vecN_dbl_2 + vecN_dbl_6 = " << vecN_dbl_2 + vecN_dbl_6 << std::endl;
    std::cout << "vecN_dbl_2 - vecN_dbl_6 = " << vecN_dbl_2 - vecN_dbl_6 << std::endl;
    std::cout << "       2.0 * vecN_dbl_6 = " << 2.0 * vecN_dbl_6 << std::endl;
    std::cout << "vecN_flt_1 * 2.0        = " << vecN_flt_1 * 2.0 << std::endl;
    std::cout << "vecN_flt_1 / 2.0        = " << vecN_flt_1 / 2.0 << std::endl;
    std::cout << "vecN_cmplx_1 + vecN_cmplx_2 = " << vecN_cmplx_1 + vecN_cmplx_2 << std::endl;
    std::cout << "vecN_cmplx_2 * 3.0          = " << vecN_cmplx_2 * 3.0 << std::endl;
    std::cout << "         3.0 * vecN_cmplx_2 = " << 3.0 * vecN_cmplx_2 << std::endl;

    std::cout << "\nVectorN operations:" << std::endl;
    VectorN<double, 5> vecN_dbl_4_almost_equal(vecN_dbl_4);
    vecN_dbl_4_almost_equal[0] = vecN_dbl_4_almost_equal[0] + 1e-6;
    std::cout << "vecN_dbl_4              = " << vecN_dbl_4 << std::endl;
    std::cout << "vecN_dbl_4_almost_equal = " << vecN_dbl_4_almost_equal << std::endl;
    std::cout << "IsEqual(vecN_dbl_4, vecN_dbl_4_almost_equal, 1e-07) = " << vecN_dbl_4.IsEqualTo(vecN_dbl_4_almost_equal, 1e-07) << std::endl;
    std::cout << "IsEqual(vecN_dbl_4, vecN_dbl_4_almost_equal, 1e-05) = " << vecN_dbl_4.IsEqualTo(vecN_dbl_4_almost_equal, 1e-05) << std::endl;

    std::cout << "NormL2(vecN_dbl_3) = " << vecN_dbl_3.NormL2() << std::endl;
    std::cout << "vecN_dbl_2.ScalarProductCartesian(vecN_dbl_6) = " << vecN_dbl_2.ScalarProductCartesian(vecN_dbl_6) << std::endl;
    std::cout << "vecN_dbl_3.AngleToVector(vecN_dbl_6)          = " << vecN_dbl_3.AngleToVector(vecN_dbl_6) << std::endl;

    std::cout << "\ndouble vector output:\n";
    std::cout << "std::cout << vecN_dbl_6                  => " << vecN_dbl_6 << std::endl;
    std::cout << "std::cout << vecN_dbl_6.to_string(10, 5) => " << vecN_dbl_6.to_string(10, 5) << std::endl;
    std::cout << "vecN_dbl_6.Print(std::cout, 7, 3)        => "; vecN_dbl_6.Print(std::cout, 7, 3);

    std::cout << "\nComplex vector output:\n";
    std::cout << vecN_cmplx_2 << std::endl;
    std::cout << vecN_cmplx_2.to_string(10, 5) << std::endl;
    vecN_cmplx_2.Print(std::cout, 7, 3);  
}