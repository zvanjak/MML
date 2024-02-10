# MatrixNM class

Class representing matrix of statically define type and size (as template parameters).

## Defined typedefs for easier use

MatrixNM class
~~~ c++
    typedef MatrixNM<float, 2, 2> Matrix22Flt;
    typedef MatrixNM<float, 3, 3> Matrix33Flt;
    typedef MatrixNM<float, 4, 4> Matrix44Flt;

    typedef MatrixNM<double, 2, 2> Matrix22Dbl;
    typedef MatrixNM<double, 3, 3> Matrix33Dbl;
    typedef MatrixNM<double, 4, 4> Matrix44Dbl;

    typedef MatrixNM<Complex, 2, 2> Matrix22Complex;
    typedef MatrixNM<Complex, 3, 3> Matrix33Complex;
    typedef MatrixNM<Complex, 4, 4> Matrix44Complex;

    typedef MatrixNM<float, 2, 2> Mat22F;
    typedef MatrixNM<float, 3, 3> Mat33F;
    typedef MatrixNM<float, 4, 4> Mat44F;

    typedef MatrixNM<double, 2, 2> Mat22D;
    typedef MatrixNM<double, 3, 3> Mat33D;
    typedef MatrixNM<double, 4, 4> Mat44D;

    typedef MatrixNM<Complex, 2, 2> Mat22C;
    typedef MatrixNM<Complex, 3, 3> Mat33C;
    typedef MatrixNM<Complex, 4, 4> Mat44C;
~~~