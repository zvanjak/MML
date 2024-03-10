# ODE system

Class representing system of ordinary differential equations (ODEs).

## Interface definition

~~~C++
class IODESystem
{
public:
	virtual int     getDim() const = 0;
	virtual void    derivs(const Real, const Vector<Real>&, Vector<Real>&) const = 0;
	void operator()(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const { return derivs(t, x, dxdt); }
};
~~~

## Base class

~~~C++
class ODESystem : public IODESystem
{
protected:
	int _dim;
	void (*_func)(Real, const Vector<Real>&, Vector<Real>&);

public:
	ODESystem() : _dim(0), _func(nullptr) { }
	ODESystem(int n, void (*inFunc)(Real, const Vector<Real>&, Vector<Real>&)) : _dim(n), _func(inFunc) { }

	int getDim() const { return _dim; }
	void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const
	{
		_func(t, x, dxdt);
	}
};
~~~

## Systems with Jacobian

~~~C++
class ODESystemWithJacobian : public ODESystem
{
private:
	void (*_funcJac)(const Real, const Vector<Real>&, Vector<Real>&, Matrix<Real>&);

public:
	ODESystemWithJacobian() { }
	ODESystemWithJacobian(int n,
		void (*inFunc)(Real, const Vector<Real>&, Vector<Real>&),
		void (*inFuncJac)(const Real t, const Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx)
	) : ODESystem(n, inFunc), _funcJac(inFuncJac) { }

	void jacobian(const Real t, Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx)
	{
		_funcJac(t, x, dxdt, dydx);
	}
};

// calculates needed Jacobian numerically
class ODESystemWithNumJacobian : public ODESystem
{
public:
	ODESystemWithNumJacobian() { }
	ODESystemWithNumJacobian(int n, void (*inFunc)(Real, const Vector<Real>&, Vector<Real>&) ) 
		: ODESystem(n, inFunc) { }

	void jacobian(const Real t, Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx)
	{
		// formirati lokalnu vektorsku funkciju na bazi derivs()
//			_funcJac(t, x, dxdt, dydx);
	}
};
~~~

## ODE solution class

~~~C++
class ODESystemSolution
{
	// first values are initial conditions
public:
	int  _sys_dim;
	int  _count;
	Real _x1, _x2;
	Vector<Real> _xval;
	Matrix<Real> _yval;

	ODESystemSolution() {}
	ODESystemSolution(Real x1, Real x2, int dim, int maxSteps) : _sys_dim(dim), _count(maxSteps + 1), _x1(x1), _x2(x2)

	template<int N>
	ParametricCurveInterpolated<N> getSolutionAsParametricCurve() const;

	template<int N>
	SplineInterpParametricCurve<N> getSolutionAsSplineParametricCurve() const;

	LinearInterpRealFunc getSolutionAsLinearInterp(int component) const;
	PolynomInterpRealFunc getSolutionAsPolynomInterp(int component, int polynomDegree) const;
	SplineInterpRealFunc getSolutionAsSplineInterp(int component) const;

	bool Serialize(std::string fileName, std::string title) const;
	bool SerializeAsParametricCurve3D(std::string fileName, std::string title) const;
};
~~~
