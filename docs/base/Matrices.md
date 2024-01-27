# Matrix types

## Defined matrix types

Basic types representing matrices:
- [Matrix\<class T>](/docs/base/Matrix.md) - matrix with run-time defined size
- [MatrixNM\<class T, int N, int M>>](/docs/base/MatrixNM.md) - matrix with compile-time defined size

Additional types
- MatrixSymm
- MatrixTridiag
- MatrixBandDiag


## Defined typedefs for easier use

Matrix class
~~~ c++
    typedef Matrix<int>     MatrixInt;
    typedef Matrix<float>   MatrixFlt;
    typedef Matrix<double>  MatrixDbl;
    typedef Matrix<Complex> MatrixComplex;

    typedef Matrix<int>     MatI;
    typedef Matrix<float>   MatF;
    typedef Matrix<double>  MatD;
    typedef Matrix<Complex> MatC; 
~~~
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
