#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector.h"
#include "base/BaseUtils.h"
#endif

using namespace MML;

void Demo_Vector()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                            VECTOR                             ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    
    std::vector<Real> std_vec{-1.0, 5.0, -2.0, 10.0, 4.0};
    float  arr[5] = {-1.0, 5.0, -2.0, 10.0, 4.0};

    Vector<Real>    vec_dbl_1(5);                       // init vector with 5 elements
    Vector<Real>    vec_dbl_2(5, 3.14159);              // init with constant value
    Vector<Real>    vec_dbl_3({ 1.0, 2.0, 3.0 });       // init with list of values
    Vector<Real>    vec_dbl_4(vec_dbl_3);               // init with copy ctor
    Vector<Real>    vec_dbl_5 = vec_dbl_2;              // init with assignment
    Vector<Real>    vec_dbl_6(std_vec);                 // init with std::vector<>
    VectorFlt       vec_flt_1(5, arr);                  // init with C/C++ array
    Vector<Real>    vec_unit(5, 3);                     // unit vector with 1.0 at index 3

    Vector<Complex>  vec_cmplx_1({ 1.0, 2.0, 3.0 });     // init with list of real values
    VecC             vec_cmplx_2({ Complex(1,1), 
                                        Complex(-1,2), 
                                        Complex(2, -0.5) });  // init with list of complex values

    std::cout << "\nVector - declarations and initializations:" << std::endl;

    std::cout << "std::vector<Real> std_vec{-1.0, 5.0, -2.0, 10.0, 4.0};" << std::endl;
    std::cout << "float  arr[5] = {-1.0, 5.0, -2.0, 10.0, 4.0};" << std::endl;

    std::cout << "Vector<Real>    vec_dbl_1(5);                 vec_dbl_1 = " << vec_dbl_1 << std::endl;
    std::cout << "Vector<Real>    vec_dbl_2(5, 3.14159);        vec_dbl_2 = " << vec_dbl_2 << std::endl;
    std::cout << "Vector<Real>    vec_dbl_3({ 1.0, 2.0, 3.0 }); vec_dbl_3 = " << vec_dbl_3 << std::endl;
    std::cout << "Vector<Real>    vec_dbl_4(vec_dbl_3);         vec_dbl_4 = " << vec_dbl_4 << std::endl;
    std::cout << "Vector<Real>    vec_dbl_5 = vec_dbl_2;        vec_dbl_5 = " << vec_dbl_5 << std::endl;
    std::cout << "Vector<Real>    vec_dbl_6(std_vec);           vec_dbl_6 = " << vec_dbl_6 << std::endl;
    std::cout << "VectorFlt       vec_flt_1(5, arr);            vec_flt_1 = " << vec_flt_1 << std::endl;
    std::cout << "Vector<Real>    vec_unit(5, 3);               vec_unit  = " << vec_unit << std::endl;
    std::cout << "vec_cmplx_1 = " << vec_cmplx_1 << std::endl;
    std::cout << "vec_cmplx_2 = " << vec_cmplx_2 << std::endl;

    std::cout << "\nVector - arithmetic operations:" << std::endl;
    std::cout << "vec_dbl_2 + vec_dbl_6 = " << vec_dbl_2 + vec_dbl_6 << std::endl;
    std::cout << "vec_dbl_2 - vec_dbl_6 = " << vec_dbl_2 - vec_dbl_6 << std::endl;
    std::cout << "      2.0 * vec_dbl_6 = " << 2.0 * vec_dbl_6 << std::endl;
    std::cout << "vec_flt_1 * 2.0       = " << vec_flt_1 * 2.0 << std::endl;
    std::cout << "vec_flt_1 / 2.0       = " << vec_flt_1 / 2.0 << std::endl;
    std::cout << "vec_cmplx_1 + vec_cmplx_2 = " << vec_cmplx_1 + vec_cmplx_2 << std::endl;
    std::cout << "vec_cmplx_2 * 3.0         = " << vec_cmplx_2 * 3.0 << std::endl;
    std::cout << "        3.0 * vec_cmplx_2 = " << 3.0 * vec_cmplx_2 << std::endl;

    std::cout << "\nHandling exceptions:" << std::endl;
    try {
      std::cout << "vec_dbl_1 + vec_dbl_3 = " << vec_dbl_1 + vec_dbl_3 << std::endl;
    }
    catch (VectorDimensionError& ) {
      std::cout << "\nCan't add vectors of different dimension!\n";
    }

    std::cout << "\nVector operations:" << std::endl;
    Vector<Real> vec_dbl_4_almost_equal(vec_dbl_4);
    vec_dbl_4_almost_equal[0] = vec_dbl_4_almost_equal[0] + 1e-6;
    std::cout << "vec_dbl_4              = " << vec_dbl_4 << std::endl;
    std::cout << "vec_dbl_4_almost_equal = " << vec_dbl_4_almost_equal << std::endl;
    std::cout << "IsEqual(vec_dbl_4, vec_dbl_4_almost_equal, 1e-07) = " << vec_dbl_4.IsEqualTo(vec_dbl_4_almost_equal, 1e-07) << std::endl;
    std::cout << "IsEqual(vec_dbl_4, vec_dbl_4_almost_equal, 1e-05) = " << vec_dbl_4.IsEqualTo(vec_dbl_4_almost_equal, 1e-05) << std::endl;

    std::cout << "NormL2(vec_dbl_3) = " << vec_dbl_3.NormL2() << std::endl;
    std::cout << "Utils::ScalarProduct(vec_dbl_2, vec_dbl_6) = " << Utils::ScalarProduct(vec_dbl_2, vec_dbl_6) << std::endl;
    std::cout << "Utils::VectorsAngle(vec_dbl_2, vecN_dbl_6) = " << Utils::VectorsAngle(vec_dbl_2, vec_dbl_6) << std::endl;

    std::cout << "\nReal vector output:\n";
    std::cout << "std::cout << vec_dbl_6                  => " << vec_dbl_6 << std::endl;
    std::cout << "std::cout << vec_dbl_6.to_string(10, 5) => " << vec_dbl_6.to_string(10, 5) << std::endl;
    std::cout << "vec_dbl_6.Print(std::cout, 7, 3)        => "; vec_dbl_6.Print(std::cout, 7, 3);

    std::cout << "\n\nComplex vector output:\n";
    std::cout << vec_cmplx_2 << std::endl;
    std::cout << vec_cmplx_2.to_string(10, 5) << std::endl;
    vec_cmplx_2.Print(std::cout, 7, 3);

}