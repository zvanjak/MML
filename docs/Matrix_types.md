# Matrix types

## Defined matrix types

Basic types representing matrices:
- Matrix\<class T> - matrix with run-time defined size
- MatrixNM\<class T, int N, int M> - matrix with compile-time defined size

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
VectorN class
~~~ c++
    typedef VectorN<float, 2> Vector2Flt;
    typedef VectorN<float, 3> Vector3Flt;
    typedef VectorN<float, 4> Vector4Flt;

    typedef VectorN<double, 2> Vector2Dbl;
    typedef VectorN<double, 3> Vector3Dbl;
    typedef VectorN<double, 4> Vector4Dbl;

    typedef VectorN<Complex, 2> Vector2Complex;
    typedef VectorN<Complex, 3> Vector3Complex;
    typedef VectorN<Complex, 4> Vector4Complex;

    typedef VectorN<float, 2> Vec2F;
    typedef VectorN<float, 3> Vec3F;
    typedef VectorN<float, 4> Vec4F;

    typedef VectorN<double, 2> Vec2D;
    typedef VectorN<double, 3> Vec3D;
    typedef VectorN<double, 4> Vec4D;

    typedef VectorN<Complex, 2> Vec2C;
    typedef VectorN<Complex, 3> Vec3C;
    typedef VectorN<Complex, 4> Vec4C;
~~~
