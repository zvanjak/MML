# Scalar and vector field operations

Numerical calculation of standard scalar and vector field operations - gradient, divergence, curl and Laplacian.
The operations are available for Cartesian, spherical and cylindrical coordinate systems. 
Gradient and divergence can be calculated for fields defined on general coordinate systems, if provided with metric tensor.

## Gradient

~~~C++
// General gradient (requires metric tensor)
template<int N>
static VectorN<Real, N> Gradient(IScalarFunction<N>& scalarField, const VectorN<Real, N>& pos, const MetricTensorField<N>& metricTensorField)

// Cartesian gradient
template<int N>
static VectorN<Real, N> GradientCart(const IScalarFunction<N>& scalarField, const VectorN<Real, N>& pos)
template<int N>
static VectorN<Real, N> GradientCart(const IScalarFunction<N>& scalarField, const VectorN<Real, N>& pos, int der_order)

// Spherical gradient
static Vector3Spherical GradientSpher(const IScalarFunction<3>& scalarField, const Vector3Spherical& pos)
static Vector3Spherical GradientSpher(const IScalarFunction<3>& scalarField, const Vector3Spherical& pos, int der_order)

// Cylindrical gradient
static Vector3Cylindrical GradientCyl(const IScalarFunction<3>& scalarField, const Vector3Cylindrical& pos)
static Vector3Cylindrical GradientCyl(const IScalarFunction<3>& scalarField, const Vector3Cylindrical& pos, int der_order)
~~~

## Divergence

~~~C++
// General divergence (requires metric tensor)
template<int N>
static Real Divergence(const IVectorFunction<N>& vectorField, const VectorN<Real, N>& pos, const MetricTensorField<N>& metricTensorField)

// Cartesian divergence
template<int N>
static Real DivCart(const IVectorFunction<N>& vectorField, const VectorN<Real, N>& pos)

// Spherical divergence
static Real DivSpher(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& x)

// Cylindrical divergence
static Real DivCyl(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& x)
~~~

## Curl

~~~C++
static Vector3Cartesian CurlCart(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& pos)

static Vector3Spherical CurlSpher(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& pos)

static Vector3Cylindrical CurlCyl(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& pos)
~~~

## Laplacian

~~~C++
template<int N>
static Real LaplacianCart(const IScalarFunction<N>& scalarField, const VectorN<Real, N>& pos)

static Real LaplacianSpher(const IScalarFunction<3>& scalarField, const Vector3Spherical& pos)

static Real LaplacianCyl(const IScalarFunction<3>& scalarField, const Vector3Cylindrical& pos)
~~~

## Example usage

