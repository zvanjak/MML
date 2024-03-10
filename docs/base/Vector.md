# Vector class

Class representing vector of given type, with dynamic size.
Based on std::vector, supports all expected arithmetic operations.

**Interesting operations**
~~~c++
    void Resize(int newLen)
    
    VectorN GetAsUnitVector() const
    bool IsEqual(const VectorN &b, Type eps) const

    Type ScalarProductCartesian(const VectorN &b) const
    Type AngleToVector(const VectorN &b) const
    Type NormL2() const

    void Print(std::ostream &out, int width, int precision) const
    std::string to_string(int width, int precision) const
~~~

## Initializing/creating vectors
Vector class
~~~ c++
    std::vector<double> std_vec{-1.0, 5.0, -2.0, 10.0, 4.0};
    float  arr[5] = {-1.0, 5.0, -2.0, 10.0, 4.0};

    MML::Vector<double>   vec_dbl_1(5);                       // init vector with 5 elements
    MML::VectorDbl        vec_dbl_2(5, 3.14159);              // init with constant value
    MML::VecD             vec_dbl_3({ 1.0, 2.0, 3.0 });       // init with list of values
    MML::Vector<double>   vec_dbl_4(vec_dbl_3);               // init with copy ctor
    MML::Vector<double>   vec_dbl_5 = vec_dbl_2;              // init with assignment
    MML::Vector<double>   vec_dbl_6(std_vec);                 // init with std::vector<>
    MML::VectorFlt        vec_flt_1(5, arr);                  // init with C/C++ array

    MML::Vector<Complex>  vec_cmplx_1({ 1.0, 2.0, 3.0 });     // init with list of real values
    MML::VecC             vec_cmplx_2({ Complex(1,1), 
                                        Complex(-1,2), 
                                        Complex(2, -0.5) });  // init with list of complex values
~~~

## Math operations with vector objects

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

 ## Operations on vector objects
~~~ c++
    MML::Vector<double> vec_dbl_4_almost_equal(vec_dbl_4);
    vec_dbl_4_almost_equal[0] = vec_dbl_4_almost_equal[0] + 1e-6;
    std::cout << "vec_dbl_4              = " << vec_dbl_4 << std::endl;
    std::cout << "vec_dbl_4_almost_equal = " << vec_dbl_4_almost_equal << std::endl;

    std::cout << "IsEqual(vec_dbl_4, vec_dbl_4_almost_equal, 1e-07) = " << vec_dbl_4.IsEqual(vec_dbl_4_almost_equal, 1e-07) << std::endl;
    std::cout << "IsEqual(vec_dbl_4, vec_dbl_4_almost_equal, 1e-05) = " << vec_dbl_4.IsEqual(vec_dbl_4_almost_equal, 1e-05) << std::endl;

    std::cout << "NormL2(vec_dbl_3) = " << vec_dbl_3.NormL2() << std::endl;
    std::cout << "vec_dbl_2.ScalarProductCartesian(vec_dbl_6) = " << vec_dbl_2.ScalarProductCartesian(vec_dbl_6) << std::endl;
~~~
## Handling exceptions
~~~ c++
    try {
      std::cout << "a + b = " << vec_dbl_1 + vec_dbl_3 << std::endl;
    }
    catch (MML::VectorDimensionError& ) {
      std::cout << "\nCan't add vectors of different dimension!\n";
    }
~~~

## Vector IO

~~~ c++
    std::cout << "\nReal vector output:\n";
    std::cout << vec_dbl_6 << std::endl;
    std::cout << vec_dbl_6.to_string(10, 5) << std::endl;
    vec_dbl_6.Print(std::cout, 7, 3);

    std::cout << "\nComplex vector output:\n";
    std::cout << vec_cmplx_2 << std::endl;
    std::cout << vec_cmplx_2.to_string(10, 5) << std::endl;
    vec_cmplx_2.Print(std::cout, 7, 3);
~~~

## Defined typedefs for easier use

Vector class
~~~ c++
    typedef Vector<int>     VectorInt;
    typedef Vector<float>   VectorFlt;
    typedef Vector<double>  VectorDbl;
    typedef Vector<Complex> VectorComplex;

    typedef Vector<int>     VecI;
    typedef Vector<float>   VecF;
    typedef Vector<double>  VecD;
    typedef Vector<Complex> VecC; 
~~~