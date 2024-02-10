# Numerical derivation

Functions for performing numerical derivation of order 1, 2, 4, 6, 8.

Types of functions that can be derived are:
- IRealFunction
- IScalarFunction\<int N>
- IVectorFunction\<int N>
- IParametricCurve\<int N>
- ITnesorField\<int N> - of order 2, 3, 4 and 5

Example declarations for NDer1 variant:
~~~c++
//////////////////////////              RealFunction            //////////////////////////
static Real NDer1(const IRealFunction &f, Real x, Real* error = nullptr);
static Real NDer1(const IRealFunction &f, Real x, Real h, Real* error = nullptr);
static Real NSecDer1(const IRealFunction &f, Real x, Real* error = nullptr)
static Real NSecDer1(const IRealFunction &f, Real x, Real h, Real* error = nullptr)
static Real NThirdDer1(const IRealFunction &f, Real x, Real* error = nullptr)
static Real NThirdDer1(const IRealFunction &f, Real x, Real h, Real* error = nullptr)

//////////////////////////             ScalarFunction           //////////////////////////
template <int N> static Real NDer1Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, Real *error = nullptr)
template <int N> static Real NDer1Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, Real h, Real *error = nullptr)
template <int N> static Real NSecDer1Partial(const IScalarFunction<N> &f, int der_ind1, int der_ind2, const VectorN<Real, N> &point, Real *error = nullptr)
template <int N> static Real NSecDer1Partial(const IScalarFunction<N> &f, int der_ind1, int der_ind2, const VectorN<Real, N> &point, Real h, Real *error = nullptr)

template <int N> static VectorN<Real, N> NDer1PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
template <int N> static VectorN<Real, N> NDer1PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, Real h, VectorN<Real, N> *error = nullptr)

//////////////////////////             VectorFunction           //////////////////////////
template <int N> static Real NDer1Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, Real *error = nullptr)
template <int N> static Real NDer1Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, Real h, Real *error = nullptr)
template <int N> static VectorN<Real, N> NDer1PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
template <int N> static VectorN<Real, N> NDer1PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, Real h, VectorN<Real, N> *error = nullptr)
template <int N> static MatrixNM<Real, N,N> NDer1PartialAllByAll(const IVectorFunction<N> &f, const VectorN<Real, N> &point, MatrixNM<Real, N,N> *error = nullptr)
template <int N> static MatrixNM<Real, N,N> NDer1PartialAllByAll(const IVectorFunction<N> &f, const VectorN<Real, N> &point, Real h, MatrixNM<Real, N,N> *error = nullptr)

//////////////////////////             TensorField           //////////////////////////
template <int N>
static Real NDer1Partial(const ITensorField2<N> &f, int i, int j, int deriv_index, const VectorN<Real, N> &point, Real *error = nullptr)
template <int N>
static Real NDer1Partial(const ITensorField2<N> &f, int i, int j, int deriv_index, const VectorN<Real, N> &point, Real h, Real *error = nullptr)
// also for TensorField3, TensorField4, TensorField5

/////////////////////////             ParametricCurve           /////////////////////////
template <int N> static VectorN<Real, N> NDer1(const IParametricCurve<N> &f, Real t, Real* error = nullptr)
template <int N> static VectorN<Real, N> NDer1(const IParametricCurve<N> &f, Real t, Real h, Real* error = nullptr)
template <int N> static VectorN<Real, N> NSecDer1(const IParametricCurve<N> &f, Real x, Real* error = nullptr)
template <int N> static VectorN<Real, N> NSecDer1(const IParametricCurve<N> &f, Real x, Real h, Real* error = nullptr)
template <int N> static VectorN<Real, N> NThirdDer1(const IParametricCurve<N> &f, Real x, Real* error = nullptr)
template <int N> static VectorN<Real, N> NThirdDer1(const IParametricCurve<N> &f, Real x, Real h, Real* error = nullptr)
~~~
