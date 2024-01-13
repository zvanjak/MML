# VectorN objects

## Initializing/creating VectorN objects
Vector class
~~~ c++
    std::vector<double> std_vec{-1.0, 5.0, -2.0};
    float  arr[3] = {-1.0, 5.0, -2.0};

    MML::VectorN<double, 5> vec_dbl_1;                       // init vector with 5 elements
    MML::Vector3Dbl         vec_dbl_2(3.14159);              // init with constant value
    MML::Vec3D              vec_dbl_3({ 1.0, 2.0, 3.0 });    // init with list of values
    MML::VectorN<double, 5> vec_dbl_4(vec_dbl_1);            // init with copy ctor
    MML::VectorN<double, 3> vec_dbl_5 = vec_dbl_2;           // init with assignment
    MML::VectorN<double, 3> vec_dbl_6(std_vec);              // init with std::vector<>
    MML::Vector3Flt         vec_flt_1(arr);                  // init with C/C++ array

    MML::VectorN<Complex, 3>    vec_cmplx_1({ 1.0, 2.0, 3.0 });     // init with list of double values
    MML::Vec3C                  vec_cmplx_2({ Complex(1,1), 
                                              Complex(-1,2), 
                                              Complex(2, -0.5) });  // init with list of complex values
~~~

## Math operations with VectorN objects

~~~ c++
    std::cout << "vec_dbl_2 + vec_dbl_6 = " << vec_dbl_2 + vec_dbl_6 << std::endl;
    std::cout << "vec_dbl_2 - vec_dbl_6 = " << vec_dbl_2 - vec_dbl_6 << std::endl;
    std::cout << "      2.0 * vec_dbl_6 = " << 2.0 * vec_dbl_6 << std::endl;
    std::cout << "vec_flt_1 * 2.0       = " << vec_flt_1 * 2.0 << std::endl;
    std::cout << "vec_flt_1 / 2.0       = " << vec_flt_1 / 2.0 << std::endl;
    std::cout << "vec_cmplx_1 + vec_cmplx_2 = " << vec_cmplx_1 + vec_cmplx_2 << std::endl;
    std::cout << "vec_cmplx_2 * 3.0         = " << vec_cmplx_2 * 3.0 << std::endl;
    std::cout << "        3.0 * vec_cmplx_2 = " << 3.0 * vec_cmplx_2 << std::endl;
 ~~~

 ## Operations on VectorN objects
~~~ c++
    std::cout << "\nVector operations:" << std::endl;
    MML::VectorN<double, 5> vec_dbl_4_almost_equal(vec_dbl_4);
    vec_dbl_4_almost_equal[0] = vec_dbl_4_almost_equal[0] + 1e-6;
    std::cout << "vec_dbl_4              = " << vec_dbl_4 << std::endl;
    std::cout << "vec_dbl_4_almost_equal = " << vec_dbl_4_almost_equal << std::endl;
    std::cout << "IsEqual(vec_dbl_4, vec_dbl_4_almost_equal, 1e-07) = " << vec_dbl_4.IsEqual(vec_dbl_4_almost_equal, 1e-07) << std::endl;
    std::cout << "IsEqual(vec_dbl_4, vec_dbl_4_almost_equal, 1e-05) = " << vec_dbl_4.IsEqual(vec_dbl_4_almost_equal, 1e-05) << std::endl;

    std::cout << "NormL2(vec_dbl_3) = " << vec_dbl_3.NormL2() << std::endl;
    std::cout << "vec_dbl_2.ScalarProductCartesian(vec_dbl_6) = " << vec_dbl_2.ScalarProductCartesian(vec_dbl_6) << std::endl;
~~~

## VectorN IO

~~~ c++
    std::cout << "\nDouble vector output:\n";
    std::cout << vec_dbl_6 << std::endl;
    std::cout << vec_dbl_6.to_string(10, 5) << std::endl;
    vec_dbl_6.Print(std::cout, 7, 3);

    std::cout << "\nComplex vector output:\n";
    std::cout << vec_cmplx_2 << std::endl;
    std::cout << vec_cmplx_2.to_string(10, 5) << std::endl;
    vec_cmplx_2.Print(std::cout, 7, 3);  
~~~