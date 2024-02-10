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

