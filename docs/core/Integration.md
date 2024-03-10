# Numerical integration

Set of classes and functions to perform numerical integration of functions.

## Implemented algorithms for real function integration

enum IntegrationMethod { TRAP, SIMPSON, ROMBERG, GAUSS10 };

### Trapezoidal integration

~~~c++
static Real IntegrateTrap(const IRealFunction& func, const Real a, const Real b, Real req_eps)

// Returns the integral of the function func from a to b. The parameters EPS can be set to the
// desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum allowed
// number of steps. Integration is performed by the trapezoidal rule.

// Unsophisticated as it is, routine qtrap is in fact a fairly robust way of doing
// integrals of functions that are not very smooth. Increased sophistication will usually
// translate into a higher-order method whose efficiency will be greater only for
// sufficiently smooth integrands. qtrap is the method of choice, e.g., for an integrand
// which is a function of a variable that is linearly interpolated between measured data
// points. Be sure that you do not require too stringent an EPS, however: If qtrap takes
// too many steps in trying to achieve your required accuracy, accumulated roundoff
// errors may start increasing, and the routine may never converge. 
// Value 1e-6 is just on the edge of trouble for most 32-bit machines; it is achievable when the
// convergence is moderately rapid, but not otherwise.
~~~

### Simpson's rule

~~~C++
static Real IntegrateSimpson(const IRealFunction& func, const Real a, const Real b, Real req_eps)

// Returns the integral of the function func from a to b. The parameters EPS can be set to the
// desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum allowed
// number of steps. Integration is performed by Simpson’s rule.

// The routine qsimp will in general be more efficient than qtrap (i.e., require
// fewer function evaluations) when the function to be integrated has a finite 4th
// derivative (i.e., a continuous 3rd derivative). The combination of qsimp and its
// necessary workhorse TrapRefine is a good one for light-duty work.
~~~

### Romberg integration

~~~C++
static Real IntegrateRomberg(const IRealFunction& func, const Real a, const Real b, Real req_eps)

// Returns the integral of the function func from a to b. Integration is performed by Romberg’s
// method of order 2K, where, e.g., K=2 is Simpson’s rule.

// The routine IntegrateRomberg, along with its required TrapRefine and polint, is quite
// powerful for sufficiently smooth (e.g., analytic) integrands, integrated over intervals
// which contain no singularities, and where the enRealoints are also nonsingular. qromb,
// in such circumstances, takes many, many fewer function evaluations than either of
// the routines in x4.2
~~~

### Gaussian quadrature

~~~C++
static Real IntegrateGauss10(const IRealFunction& func, const Real a, const Real b)

// Returns the integral of the function func between a and b, by ten-point GaussLegendre integration: 
// the function is evaluated exactly ten times at interior points in the range of integration.  
~~~

## Example usage

~~~C++
RealFunction f1{ [](Real x) { return (Real)(sin(x) * (1.0 + 0.5 * x * x)); } };
RealFunction f1_integral{ [](Real x) { return (Real)(x * (-0.5 * x * cos(x) + sin(x))); } };

double a = 0.0;
double b = 10.0;
double int_trap  = IntegrateTrap(f1, a, b);
double int_simp  = IntegrateSimpson(f1, a, b);
double int_romb  = IntegrateRomberg(f1, a, b);
double int_gauss = IntegrateGauss10(f1, a, b);
// we can use default Integrate routine (set to IntegrateSimpson), requires precision
double int_def = Integrate(f1, a, b, 1e-04);

std::cout << "Integrating function f1 from " << a << " to " << b << std::endl;
std::cout << "Exact integral   = " << f1_integral(b) - f1_integral(a) << std::endl;
std::cout << "IntegrateTrap    = " << int_trap << std::endl;
std::cout << "IntegrateSimpson = " << int_simp << std::endl;
std::cout << "IntegrateRomberg = " << int_romb << std::endl;
std::cout << "IntegrateGauss10 = " << int_gauss << std::endl;

/* OUTPUT
Integrating function f1 from 0 to 10
Exact integral   = 36.5134
IntegrateTrap    = 36.5133
IntegrateSimpson = 36.5134
IntegrateRomberg = 36.5134
IntegrateGauss10 = 36.5134
*/
~~~
