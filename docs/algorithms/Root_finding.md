# Root finding

List of functions for root finding:
- zbrac
- zbrak
- FindRootBisection
- FindRootFalsePosition
- FindRootSecant
- FindRootRidders
- FindRootBrent
- FindRootNewton
- FindRootSafe

## Functions

~~~C++
// Given a function func and an initial guessed range x1 to x2, the routine expands
// the range geometrically until a root is bracketed by the returned values x1 and x2(in which
// case zbrac returns true) or until the range becomes unacceptably large(in which case zbrac
// returns false).
static bool zbrac(const IRealFunction& func, double& x1, double& x2)
~~~

~~~C++
// Given a function or functor fx defined on the interval[x1, x2], subdivide the interval into
// n equally spaced segments, and search for zero crossings of the function.nroot will be set
// to the number of bracketing pairs found.If it is positive, the arrays xb1[0..nroot - 1] and
// xb2[0..nroot - 1] will be filled sequentially with any bracketing pairs that are found.On input,
// these vectors may have any size, including zero; they will be resized to   nroot.
static void zbrak(const IRealFunction& fx, const Real x1, const Real x2, const int n, Vector<Real>& xb1,
	                Vector<Real>& xb2, int& nroot)
~~~

~~~C++
// Using bisection, return the root of a function or functor func known to lie between x1 and x2.
// The root will be refined until its accuracy is xacc.
static Real FindRootBisection(const IRealFunction& func, const Real x1, const Real x2, const Real xacc)
~~~

~~~C++
// Using the false - position method, return the root of a function or functor func known to lie
// between x1 and x2.The root is refined until its accuracy is xacc
static Real FindRootFalsePosition(const IRealFunction& func, const Real x1, const Real x2, const Real xacc)
~~~

~~~C++
// Using the secant method, return the root of a function or functor func thought to lie between
// x1 and x2.The root is refined until its accuracy is ˙xacc.
static Real FindRootSecant(const IRealFunction& func, const Real x1, const Real x2, const Real xacc)
~~~

~~~C++
// Using Ridders’ method, return the root of a function or functor func known to lie between x1
// and x2.The root will be refined to an approximate accuracy xacc.
static Real FindRootRidders(const IRealFunction& func, const Real x1, const Real x2, const Real xacc)
~~~

~~~C++
// Using Brent’s method, return the root of a function or functor func known to lie between x1
// and x2.The root will be refined until its accuracy is tol.
static Real FindRootBrent(IRealFunction& func, const Real x1, const Real x2, const Real tol)
~~~

~~~C++
// Using the Newton-Raphson method, return the root of a function known to lie in the interval
// x1; x2. The root will be refined until its accuracy is known within ˙xacc.
static Real FindRootNewton(const IRealFunction& funcd, const Real x1, const Real x2, const Real xacc)
~~~

~~~C++
//Using a combination of Newton - Raphson and bisection, return the root of a function bracketed
//between x1 and x2.The root will be refined until its accuracy is known within ˙xacc.funcd
//is a user - supplied struct that returns the function value as a functor and the first derivative of
//the function at the point x as the function df(see text).
static Real FindRootSafe(const IRealFunction& funcd, const Real x1, const Real x2, const Real xacc)
~~~
