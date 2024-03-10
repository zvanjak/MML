# Function types

MML introduces functions as full fledged objects, which can be passed as arguments to other functions (and are expected as arguments for many MML algorithms).
This is achieved by defining a set of interfaces for different types of functions:

- IRealFunction
- IScalarFunction\<int N>
- IVectorFunction\<int N>
- IVectorFunctionNM\<int N, int M>
- IRealToVectorFunction\<int N>
- IParametricCurve\<int N>
- IParametricSurface\<int N>
- ITensorField\<int N> - of order 2, 3, 4 and 5

In addition to these interfaces, MML provides a set of basic function objects, which implement these interfaces, 
and enable user to create function objects in a simple way, whether by using function pointers, lambda expressions, class members, or by interpolating from a set of values.

## Interfaces for function types
Interfaces for function types are defined in the following way:
~~~C++
template<typename _RetType, typename _ArgType>
class IFunction
{
public:
	virtual _RetType operator()(_ArgType) const = 0;
};

class IRealFunction : public IFunction<Real, Real>
{
public:
	virtual Real operator()(Real) const = 0;
};

template<int N>
class IScalarFunction : public IFunction<Real, const VectorN<Real, N>&>
{
public:
	virtual Real operator()(const VectorN<Real, N>& x) const = 0;
};

template<int N>
class IRealToVectorFunction : public IFunction<VectorN<Real, N>, Real>
{
public:
	virtual VectorN<Real, N> operator()(Real x) const = 0;
};

template<int N>
class IVectorFunction : public IFunction<VectorN<Real, N>, const VectorN<Real, N>&>
{
public:
	virtual VectorN<Real, N> operator()(const VectorN<Real, N>& x) const = 0;

	virtual Real operator()(const VectorN<Real, N>& x, int component) const
	{
		VectorN<Real, N> val = (*this)(x);
		return val[component];
	}
};

template<int N, int M>
class IVectorFunctionNM : public IFunction<VectorN<Real, M>, const VectorN<Real, N>&>
{
public:
	virtual VectorN<Real, M> operator()(const VectorN<Real, N>& x) const = 0;
	virtual Real operator()(const VectorN<Real, N>& x, int component) const
	{
		VectorN<Real, M> val = (*this)(x);
		return val[component];
	}
};
~~~

Parametric curve and surface interfaces are defined as follows:

~~~C++
template<int N>
class IParametricCurve : public IRealToVectorFunction<N>
{
public:
	virtual VectorN<Real, N> operator()(Real x) const = 0;

	virtual Real getMinT() const = 0;
	virtual Real getMaxT() const = 0;

	std::vector<VectorN<Real, N>> GetTrace(double t1, double t2, int numPoints) const
	{
		std::vector<VectorN<Real, N>> ret;
		double deltaT = (t2 - t1) / (numPoints - 1);
		for (Real t = t1; t <= t2; t += deltaT)
			ret.push_back((*this)(t));
		return ret;
	}
};

// simple regular surface, defined on rectangular coordinate patch
template<int N>
class IParametricSurface : public IFunction<VectorN<Real, N>, const VectorN<Real, 2>&>
{
public:
	virtual VectorN<Real, N> operator()(Real u, Real w) const = 0;

	virtual Real getMinX() const = 0;
	virtual Real getMaxX() const = 0;
	virtual Real getMinY() const = 0;
	virtual Real getMaxY() const = 0;

	virtual VectorN<Real, N> operator()(const VectorN<Real, 2>& coord) const
	{
		return operator()(coord[0], coord[1]);
	}
	bool Serialize2DCartesian(Real x1, Real x2, int numPointsX, Real y1, Real y2, int numPointsY, std::string fileName) const
	{
		return false;
	}
};
~~~

Tensor fields are defined as follows:

~~~C++
template<int N>
class ITensorField2 : public IFunction<Tensor2<N>, const VectorN<Real, N>& >
{
	int _numContravar;
	int _numCovar;
public:
	ITensorField2(int numContra, int numCo) : _numContravar(numContra), _numCovar(numCo) { }

	int getNumContravar() const { return _numContravar; }
	int getNumCovar() const { return _numCovar; }

	virtual Real    Component(int i, int j, const VectorN<Real, N>& pos) const = 0;
};

template<int N>
class ITensorField3 : public IFunction<Tensor3<N>, const VectorN<Real, N>& >
{
	int _numContravar;
	int _numCovar;
public:
	int getNumContravar() const { return _numContravar; }
	int getNumCovar() const { return _numCovar; }

	virtual Real    Component(int i, int j, int k, const VectorN<Real, N>& pos) const = 0;
};

template<int N>
class ITensorField4 : public IFunction<Tensor4<N>, const VectorN<Real, N>& >
{
	int _numContravar;
	int _numCovar;
public:
	int getNumContravar() const { return _numContravar; }
	int getNumCovar() const { return _numCovar; }

	virtual Real    Component(int i, int j, int k, int l, const VectorN<Real, N>& pos) const = 0;
};

template<int N>
class ITensorField5 : public IFunction<Tensor5<N>, const VectorN<Real, N>& >
{
	int _numContravar;
	int _numCovar;
public:
	int getNumContravar() const { return _numContravar; }
	int getNumCovar() const { return _numCovar; }

	virtual Real    Component(int i, int j, int k, int l, int m, const VectorN<Real, N>& pos) const = 0;
};
~~~

## Basic function objects, implementing above interfaces

For each interface, there are two variants of basic function objects - one that accepts function pointers, and one that accepts std::function objects.
Both variants are concrete classes, meaning that they can be instantiated (providing correct parameters), and used as function objects.
First variant enables simple usage scenarios, with function pointers or lambda expression providing the function.
Second variant enables more complex usage scenarios, with class members providing the function.

~~~C++
///////////////////////////     REAL FUNCTION      ////////////////////////////////////
class RealFunction : public IRealFunction
{
	Real(*_func)(const Real);
public:
	RealFunction(Real(*inFunc)(const Real)) : _func(inFunc) {}

	Real operator()(const Real x) const { return _func(x); }
};
class RealFunctionFromStdFunc : public IRealFunction
{
	std::function<Real(const Real)> _func;
public:
	RealFunctionFromStdFunc(std::function<Real(const Real)> inFunc) : _func(inFunc) {}

	Real operator()(const Real x) const { return _func(x); }
};

///////////////////////////     SCALAR FUNCTION       //////////////////////////////////
template<int N>
class ScalarFunction : public IScalarFunction<N>
{
	Real(*_func)(const VectorN<Real, N>&);
public:
	ScalarFunction(Real(*inFunc)(const VectorN<Real, N>&)) : _func(inFunc) {}

	Real operator()(const VectorN<Real, N>& x) const { return _func(x); }
};

template<int N>
class ScalarFunctionFromStdFunc : public IScalarFunction<N>
{
	std::function<Real(const VectorN<Real, N>&)> _func;
public:
	ScalarFunctionFromStdFunc(std::function<Real(const VectorN<Real, N>&)> inFunc) : _func(inFunc) {}

	Real operator()(const VectorN<Real, N>& x) const { return _func(x); }
};

/////////////////////////    VECTOR FUNCTION N -> N      ///////////////////////////////////
template<int N>
class VectorFunction : public IVectorFunction<N>
{
	VectorN<Real, N>(*_func)(const VectorN<Real, N>&);
public:
	VectorFunction(VectorN<Real, N>(*inFunc)(const VectorN<Real, N>&)) : _func(inFunc) {}

	VectorN<Real, N>     operator()(const VectorN<Real, N>& x) const { return _func(x); }

	MatrixNM<Real, N, N> jacobian(const VectorN<Real, N>& x) const { return Jacobian<VectorN<Real, N>, VectorN<Real, N>, N>::calc(*this, x); }
};

template<int N>
class VectorFunctionFromStdFunc : public IVectorFunction<N>
{
	std::function<VectorN<Real, N>(const VectorN<Real, N>&)> _func;
public:
	VectorFunctionFromStdFunc(std::function<VectorN<Real, N>(const VectorN<Real, N>&)>& inFunc) : _func(inFunc) {}

	VectorN<Real, N>     operator()(const VectorN<Real, N>& x) const { return _func(x); }

	MatrixNM<Real, N, N> jacobian(const VectorN<Real, N>& x) const { return Jacobian<VectorN<Real, N>, VectorN<Real, N>, N>::calc(*this, x); }
};

/////////////////////////    VECTOR FUNCTION N -> M      ///////////////////////////////////
template<int N, int M>
class VectorFunctionNM : public IVectorFunctionNM<N, M>
{
	VectorN<Real, M>(*_func)(const VectorN<Real, N>&);
public:
	VectorFunctionNM(VectorN<Real, M>(*inFunc)(const VectorN<Real, N>&)) : _func(inFunc) {}

	VectorN<Real, M>     operator()(const VectorN<Real, N>& x) const { return _func(x); }

	MatrixNM<Real, M, N> jacobian(const VectorN<Real, N>& x) const { return Jacobian<VectorN<Real, N>, VectorN<Real, M>, N>::calc(*this, x); }
};

template<int N, int M>
class VectorFunctionNMFromStdFunc : public IVectorFunctionNM<N, M>
{
	std::function<VectorN<Real, M>(const VectorN<Real, N>&)> _func;
public:
	VectorFunctionNMFromStdFunc(std::function<VectorN<Real, M>(const VectorN<Real, N>&)>& inFunc) : _func(inFunc) {}

	VectorN<Real, M>     operator()(const VectorN<Real, N>& x) const { return _func(x); }

	MatrixNM<Real, M, N> jacobian(const VectorN<Real, N>& x) const { return Jacobian<VectorN<Real, N>, VectorN<Real, M>, N>::calc(*this, x); }
};

//////////////////////     PARAMETRIC CURVE             ///////////////////////////////////
template<int N>
class   ParametricCurve : public IParametricCurve<N>
{
	Real _minT;
	Real _maxT;
	VectorN<Real, N>(*_func)(Real);
public:
	ParametricCurve(VectorN<Real, N>(*inFunc)(Real)) : _func(inFunc), _minT(Constants::NegativeInf), _maxT(Constants::PositiveInf) {}
	ParametricCurve(Real minT, Real maxT, VectorN<Real, N>(*inFunc)(Real)) : _func(inFunc), _minT(minT), _maxT(maxT) {}

	Real getMinT() const { return _minT; }
	Real getMaxT() const { return _maxT; }

	virtual VectorN<Real, N> operator()(Real x) const { return _func(x); }
};

template<int N>
class ParametricCurveFromStdFunc : public IParametricCurve<N>
{
	Real _minT;
	Real _maxT;
	std::function<VectorN<Real, N>(Real)> _func;
public:
	ParametricCurveFromStdFunc(std::function<VectorN<Real, N>(Real)>& inFunc) : _func(inFunc), _minT(Constants::NegativeInf), _maxT(Constants::PositiveInf) {}
	ParametricCurveFromStdFunc(Real minT, Real maxT, std::function<VectorN<Real, N>(Real)>& inFunc) : _func(inFunc), _minT(minT), _maxT(maxT) {}

	Real getMinT() const { return _minT; }
	Real getMaxT() const { return _maxT; }

	VectorN<Real, N> operator()(Real x) const { return _func(x); }
};

/////////////////////       PARAMETRIC SURFACE         //////////////////////////////////
template<int N>
class ParametricSurface : public IParametricSurface<N>
{
	// TODO - ensure that N is at least 3!!!
	Real _minX;
	Real _maxX;
	Real _minY;
	Real _maxY;
	VectorN<Real, N>(*_func)(Real u, Real w);

public:
	ParametricSurface(VectorN<Real, N>(*inFunc)(Real u, Real w)) : _func(inFunc), _minX(Constants::NegativeInf), _maxX(Constants::PositiveInf), _minY(Constants::NegativeInf), _maxY(Constants::PositiveInf) {}
	ParametricSurface(VectorN<Real, N>(*inFunc)(Real u, Real w), Real minX, Real maxX, Real minY, Real maxY) : _func(inFunc), _minX(minX), _maxX(maxX), _minY(minY), _maxY(maxY) {}

	VectorN<Real, N> operator()(Real u, Real w) const { return _func(u, w); }

	virtual Real getMinX() const { return _minX; }
	virtual Real getMaxX() const { return _maxX; }
	virtual Real getMinY() const { return _minY; }
	virtual Real getMaxY() const { return _maxY; }

	// TODO - double getStartY(double x) const;     // ako surface patch  nije kvadratni
};

// imati cemo i surface Discrete, kreiran od triangles?

template<int N>
class ParametricSurfaceFromStdFunc : public IParametricSurface<N>
{
	Real _minX;
	Real _maxX;
	Real _minY;
	Real _maxY;
	std::function<VectorN<Real, N>(Real u, Real w)> _func;
public:
	ParametricSurfaceFromStdFunc(std::function<VectorN<Real, N>(Real u, Real w)>& inFunc) : _func(inFunc), _minX(Constants::NegativeInf), _maxX(Constants::PositiveInf), _minY(Constants::NegativeInf), _maxY(Constants::PositiveInf) {}
	ParametricSurfaceFromStdFunc(std::function<VectorN<Real, N>(Real u, Real w)>& inFunc, Real minX, Real maxX, Real minY, Real maxY) : _func(inFunc), _minX(minX), _maxX(maxX), _minY(minY), _maxY(maxY) {}

	VectorN<Real, N> operator()(Real u, Real w) const { return _func(u, w); }

	virtual Real getMinX() const { return _minX; }
	virtual Real getMaxX() const { return _maxX; }
	virtual Real getMinY() const { return _minY; }
	virtual Real getMaxY() const { return _maxY; }
};
~~~

## How to create Function objects?

### From already existing (standalone/static) functions

~~~C++
// examples of already existing function
Real Docs_Demo_functions_RealFunc(Real x) {
  return sin(x) * (1 + x * x / 2);
}
Real Docs_Demo_functions_ScalarFunc(const VectorN<Real, 3>& x) {
	return 1 / x.NormL2();
}
// little bit more real scalar function :)
Real Docs_Demo_functions_two_masses_gravity_potential(const VectorN<Real, 3>& x)
{
  const VectorN<Real, 3> x1{ 10.0, 0.0, 0.0 };
  const VectorN<Real, 3> x2{ -10.0, 0.0, 0.0 };
  const Real m1 = 1000.0;
  const Real m2 = 1000.0;
  const Real G = 1.0;
  return -G * m1 / (x - x1).NormL2() - G * m2 / (x - x2).NormL2();
}
VectorN<Real, 3> Docs_Demo_functions_VectorFunc(const VectorN<Real, 3>& x) {
	return VectorN<Real, 3>{0, x[0] * x[1], 0};
}
VectorN<Real, 3> Docs_Demo_functions_VectorFuncNM(const VectorN<Real, 2>& x) {
	return VectorN<Real, 3>{0, x[0] * x[1], 0};
}
VectorN<Real, 3> Docs_Demo_functions_ParamCurve(Real t) {
	return VectorN<Real, 3>{t, 2 * t, 3 * t};
}
VectorN<Real, 3> Docs_Demo_functions_ParamSurface(Real x, Real y) {
	return VectorN<Real, 3>{x * y, 2 * x * y, 3 * x};
}

void Docs_Demo_functions_case_1_usage()
{
  // creating a function object from an already existing (standalone) function
  RealFunction            fReal(Docs_Demo_functions_RealFunc);
  ScalarFunction<3>       fScalar(Docs_Demo_functions_ScalarFunc);
  ScalarFunction<3>       fScalar2(Docs_Demo_functions_two_masses_gravity_potential);
  VectorFunction<3>       fVector(Docs_Demo_functions_VectorFunc);
  VectorFunctionNM<2, 3>  fVectorNM(Docs_Demo_functions_VectorFuncNM);
  ParametricCurve<3>      paramCurve(Docs_Demo_functions_ParamCurve);
  ParametricSurface<3>    paramSurface(Docs_Demo_functions_ParamSurface);
}
~~~

### Creating function directly from lambda

~~~C++
  std::cout << "**********************************************************************" << std::endl;
  std::cout << "***  CASE 2 - creating function directly from lambda" << std::endl;

  RealFunction f2{ [](Real x) { return (Real)sin(x) * (1 + x * x / 2); } };

  // creating directly different types of functions
  ScalarFunction<3> fScalar([](const VectorN<Real, 3>& x) { return x[0]; });
  ScalarFunction<3> two_masses_gravity_field_potential{ [](const VectorN<Real, 3>& x)
  {
      const VectorN<Real, 3> x1{ 10.0, 0.0, 0.0 };
      const VectorN<Real, 3> x2{ -10.0, 0.0, 0.0 };
      const Real m1 = 1000.0;
      const Real m2 = 1000.0;
      const Real G = 1.0;
      return -G * m1 / (x - x1).NormL2() - G * m2 / (x - x2).NormL2();
  } };
  VectorFunction<3>       fVector([](const VectorN<Real, 3>& x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });
  VectorFunctionNM<2, 3>  fVectorNM([](const VectorN<Real, 2>& x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });
  ParametricCurve<3>      paramCurve([](Real x) { return VectorN<Real, 3>{x, 2 * x, 3 * x}; });
  ParametricSurface<3>    paramSurface([](Real x, Real y) { return VectorN<Real, 3>{x* y, 2 * x * y, 3 * x}; });

  // evaluating the functions
  std::cout << "f2(0.0) = " << f2(0.0) << std::endl;
  std::cout << "fScalar({1, 1, 1})   = " << fScalar(VectorN<Real, 3>{1, 1, 1}) << std::endl;
  std::cout << "fScalar({-1, -1, 1}) = " << fScalar(VectorN<Real, 3>{-1, -1, 1}) << std::endl;
  std::cout << "fScalar({2, 2, 2})   = " << fScalar(VectorN<Real, 3>{2, 2, 2}) << std::endl;
  std::cout << "fVector({1, 1, 1})   = " << fVector(VectorN<Real, 3>{1, 1, 1}) << std::endl;
  std::cout << "fVector({-1, -1, 1}) = " << fVector(VectorN<Real, 3>{-1, -1, 1}) << std::endl;
  std::cout << "fVector({2, 2, 2})   = " << fVector(VectorN<Real, 3>{2, 2, 2}) << std::endl;
  std::cout << "fVectorNM({1, 1})    = " << fVectorNM(VectorN<Real, 2>{1, 1}) << std::endl;
  std::cout << "fVectorNM({-1, -1})  = " << fVectorNM(VectorN<Real, 2>{-1, -1}) << std::endl;
  std::cout << "fVectorNM({2, 2})    = " << fVectorNM(VectorN<Real, 2>{2, 2}) << std::endl;
  std::cout << "paramCurve(1.0)        = " << paramCurve(1.0) << std::endl;
  std::cout << "paramSurface(1.0, 1.0) = " << paramSurface(1.0, 1.0) << std::endl;

/* OUTPUT
f2(0.0) = 0
fScalar({1, 1, 1})   = 1
fScalar({-1, -1, 1}) = -1
fScalar({2, 2, 2})   = 2
fVector({1, 1, 1})   = [              0,               1,               0]
fVector({-1, -1, 1}) = [              0,               1,               0]
fVector({2, 2, 2})   = [              0,               4,               0]
fVectorNM({1, 1})    = [              0,               1,               0]
fVectorNM({-1, -1})  = [              0,               1,               0]
fVectorNM({2, 2})    = [              0,               4,               0]
paramCurve(1.0)        = [              1,               2,               3]
paramSurface(1.0, 1.0) = [              1,               2,               3]
*/
~~~

### Simple function, but with some parameters

~~~C++
// CASE 3 - we have simple function, but that has some additional parameters that it depend on (and that we might wanna vary easily)

class TwoMassesGravityPotential : public IScalarFunction<3>
{
  Real _m1, _m2, _G;
  VectorN<Real, 3> _x1, _x2;

public:
    TwoMassesGravityPotential(Real m1, Real m2, Real G, const VectorN<Real, 3>& x1, const VectorN<Real, 3>& x2) : _m1(m1), _m2(m2), _G(G), _x1(x1), _x2(x2) {}  

    void SetM1(Real m1) { _m1 = m1; }
    void SetM2(Real m2) { _m2 = m2; }
    void SetX1(const VectorN<Real, 3>& x1) { _x1 = x1; }
    void SetX2(const VectorN<Real, 3>& x2) { _x2 = x2; }

    Real operator()(const VectorN<Real, 3> & x) const { return -_G * (_m1 / (x - _x1).NormL2() + _m2 / (x - _x2).NormL2()); }
};

void Docs_Demo_functions_case_3_usage()
{
  std::cout << "**********************************************************************" << std::endl;
  std::cout << "***  CASE 3 - we have simple function, but that has some additional parameters that it depend on (and that we might wanna vary easily)" << std::endl;

  TwoMassesGravityPotential f3(1000.0, 1000.0, 1.0, VectorN<Real, 3>{10.0, 0.0, 0.0}, VectorN<Real, 3>{-10.0, 0.0, 0.0});

  std::cout << "f3({1, 1, 1})   = " << f3(VectorN<Real, 3>{1, 1, 1}) << std::endl;
  std::cout << "f3({-1, -1, 1}) = " << f3(VectorN<Real, 3>{-1, -1, 1}) << std::endl;
  std::cout << "f3({2, 2, 2})   = " << f3(VectorN<Real, 3>{2, 2, 2}) << std::endl;

/* OUTPUT
f3({1, 1, 1})   = -199.9312235
f3({-1, -1, 1}) = -199.9312235
f3({2, 2, 2})   = -198.9618408
*/
}
~~~

### There is an external class that has EXACT FUNCTION you want to use

~~~C++
class BigComplexClassYouCantChange1 {
public:
  double _param1, _param2;
  double Weight(double x) const { return 58 * _param1; }

  double UsefulRealFunc(double x) const { return (_param1 * cos(x) + _param2 * sin(x)) / Weight(x); }
};

void Docs_Demo_functions_case_4_usage()
{
  std::cout << "**********************************************************************" << std::endl;
  std::cout << "***  CASE 4 - there is an external class that has EXACT FUNCTION you want to use (but you can't change the class)" << std::endl;
  
  BigComplexClassYouCantChange1 bigObj;

  bigObj._param1 = 1.0;
  bigObj._param2 = 2.0;

  RealFunctionFromStdFunc f4(std::function<Real(Real)>{[&bigObj](Real x) { return bigObj.UsefulRealFunc(x); }});

  std::cout << "f4(0.0) = " << f4(0.0) << std::endl;
  std::cout << "f4(1.0) = " << f4(1.0) << std::endl;
  std::cout << "f4(2.0) = " << f4(2.0) << std::endl;

  bigObj._param1 = -2.0;
  bigObj._param2 = 5.0;

  std::cout << "After change of object parameters:" << std::endl;

  std::cout << "f4(0.0) = " << f4(0.0) << std::endl;
  std::cout << "f4(1.0) = " << f4(1.0) << std::endl;
  std::cout << "f4(2.0) = " << f4(2.0) << std::endl;

/* OUTPUT
f4(0.0) = 0.01724137931
f4(1.0) = 0.03833179785
f4(2.0) = 0.02418013823
After change of object parameters:
f4(0.0) = 0.01724137931
f4(1.0) = -0.02695474407
f4(2.0) = -0.04636880006
*/
~~~

### There is an external class that has data relevant to calculation

~~~C++
class BigComplexClassYouCantChange2 {
public:
  double _param1, _param2;
  double Weight(double x) const { return 58*_param1; }
};

// Create a helper wrapper class, inherit it from IRealFunction and use it as RealFunction
class BigComplexRealFunc2 : public IRealFunction {
  const BigComplexClassYouCantChange2& _ref;
public:
  BigComplexRealFunc2(const BigComplexClassYouCantChange2& bigClass) : _ref(bigClass) { }

  Real operator()(Real x) const {
    return (_ref._param1 * cos(x) + _ref._param2 * sin(x)) / _ref.Weight(x) ;  
  }
};

void Docs_Demo_functions_case_5_usage()
{
  std::cout << "**********************************************************************" << std::endl;
  std::cout << "***  CASE 5 - there is an external class that has data relevant to calculation (but you can't change the class)" << std::endl;
  BigComplexClassYouCantChange2 bigObj;

  bigObj._param1 = 1.0;
  bigObj._param2 = 2.0;

  BigComplexRealFunc2    f5(bigObj);     // usable RealFunction object

  std::cout << "f5(0.0) = " << f5(0.0) << std::endl;
  std::cout << "f5(1.0) = " << f5(1.0) << std::endl;
  std::cout << "f5(2.0) = " << f5(2.0) << std::endl;

  bigObj._param1 = -2.0;
  bigObj._param2 = 5.0;

  // our f5 function is now CHNAGED!!!
  std::cout << "After change of object parameters:" << std::endl;

  std::cout << "f5(0.0) = " << f5(0.0) << std::endl;
  std::cout << "f5(1.0) = " << f5(1.0) << std::endl;
  std::cout << "f5(2.0) = " << f5(2.0) << std::endl;

/* OUTPUT
f5(0.0) = 0.01724137931
f5(1.0) = 0.03833179785
f5(2.0) = 0.02418013823
After change of object parameters:
f5(0.0) = 0.01724137931
f5(1.0) = -0.02695474407
f5(2.0) = -0.04636880006
*/
~~~


### Creating interpolated from given values

Currently, you can create inteprolated RealFunction, ParametricCurve and ScalaFunction<2> as function objects from interpolated values.
Detailed explanation is given in - [Interpolated functions](/docs/core/Interpolated_functions.md).



