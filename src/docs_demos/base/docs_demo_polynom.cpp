#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Polynom.h"
#include "algorithms/RootFindingPolynoms.h"
#endif

using namespace MML;

// Demo functions for Polynoms.md runnable examples
void Docs_Demo_Polynom_Calculus();
void Docs_Demo_Polynom_RootFinding();

void Docs_Demo_Polynom()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                           POLYNOM                             ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// Initialization of polynomials
	PolynomRealFunc pol_constant({ 1 });
	PolynomRealFunc pol_linear({ 1, 2 });					
	PolynomRealFunc pol_quadratic({ 1, 2, 3 });			
	PolynomRealFunc pol_cubic({ 1, 2, 3, 4 });			
	PolynomRealFunc pol_quartic({ 1, 2, 3, 4, 5 });

	PolynomRealFunc    poly_real({ -1, 0.25, 1.3, -2, 0.4 });
	PolynomComplex poly_cmplx({ 1, 2, 3, 4 });
	Matrix2Polynom poly_mat({ 1, -1.0/2, 1./6 });            // matrix polynomial of 3rd order

	// Basic output of polynomials
	std::cout << "pol_constant  : " << pol_constant << std::endl;
	std::cout << "pol_linear    : " << pol_linear << std::endl;
	std::cout << "pol_quadratic : " << pol_quadratic << std::endl;
	std::cout << "pol_cubic     : " << pol_cubic << std::endl;
	std::cout << "pol_quartic   : " << pol_quartic << std::endl;
	std::cout << "poly_real     : " << poly_real << std::endl;
	std::cout << "poly_cmplx    : " << poly_cmplx << std::endl;
	std::cout << "poly_mat      : " << poly_mat << std::endl << std::endl;

	// Evaluation of polynomials
	std::cout << "poly_real(5.0)   = " << poly_real(1.0) << std::endl;
	std::cout << "poly_real(-2.0)  = " << poly_real(-2.0) << std::endl;
	std::cout << "poly_real(1.276) = " << poly_real(1.276) << std::endl << std::endl;
	
	std::cout << "poly_cmplx( 2 + 3i) = " << poly_cmplx(Complex(2, 3)) << std::endl;
	std::cout << "poly_cmplx(-1 - i)  = " << poly_cmplx(Complex(-1, -1)) << std::endl;
	std::cout << "poly_cmplx( 1 - 2i) = " << poly_cmplx(Complex(1, -2)) << std::endl << std::endl;

	// evaluate the polynomial at the given matrix value
	std::cout << "poly_mat({ 1, 0.5, -1.4, 2.8 }) = " << poly_mat(MatrixNM<Real, 2, 2>({ 1, 0.5, -1.4, 2.8 })) << std::endl << std::endl;

	// evaluation of matrix polynom in detail
	MatrixNM<Real, 2, 2> eval_mat({ 1, 2, 3, 4 });     // matrix to evaluate the polynomial at
	MatrixNM<Real, 2, 2> m2_2 = poly_mat(eval_mat);    // evaluate the polynomial at the given matrix value

	// Operations on polynomials
	PolynomRealFunc pol_sum = pol_quadratic + pol_cubic;
	std::cout << "( " << pol_quadratic << " )  +  ( " << pol_cubic << " )  =  " << pol_sum << std::endl;
	PolynomRealFunc pol_diff = pol_quadratic - pol_cubic;
	std::cout << "( " << pol_quadratic << " )  -  ( " << pol_cubic << " )  =  " << pol_diff << std::endl;
	PolynomRealFunc pol_prod = pol_quadratic * pol_cubic;
	std::cout << "( " << pol_quadratic << " )  *  ( " << pol_cubic << " )  =  " << pol_prod << std::endl;
	PolynomRealFunc pol_div, pol_rem;
	PolynomRealFunc::poldiv(pol_prod, pol_cubic, pol_div, pol_rem);
	std::cout << "( " << pol_prod << " )  /  ( " << pol_cubic << " )  =  " << pol_div << "  remainder: " << pol_rem << std::endl;

	PolynomRealFunc pol_sum2 = pol_quadratic * 2.0;
	PolynomRealFunc pol_diff2 = pol_quadratic / 2.0;
	PolynomRealFunc pol_prod2 = 2.0 * pol_quadratic;

	std::cout << "\nReal coef. polynom output:\n";
	std::cout << poly_real << std::endl;
	std::cout << poly_real.to_string(10, 5) << std::endl;
	poly_real.Print(std::cout, 7, 3);

	std::cout << "\nComplex coef. polynom output:\n";
	std::cout << poly_cmplx << std::endl;
	std::cout << poly_cmplx.to_string(10, 5) << std::endl;
	poly_cmplx.Print(std::cout, 7, 3);

/* OUTPUT
pol_constant  : 1
pol_linear    : 2 * x1 + 1
pol_quadratic : 3 * x^2 + 2 * x1 + 1
pol_cubic     : 4 * x^3 + 3 * x^2 + 2 * x1 + 1
pol_quartic   : 5 * x^4 + 4 * x^3 + 3 * x^2 + 2 * x1 + 1
poly_real     : 0.4 * x^4 + -2 * x^3 + 1.3 * x^2 + 0.25 * x1 + -1
poly_cmplx_r  : (4,0) * x^3 + (3,0) * x^2 + (2,0) * x1 + (1,0)
poly_cmplx    : (-2,3) * x1 + (0.5,-1)
poly_mat      : 0.166667 * x^2 + -0.5 * x1 + 1

poly_real(5.0)   = -1.05
poly_real(-2.0)  = 26.1
poly_real(1.276) = -1.65909

poly_cmplx( 2 + 3i) = (-12.5,-1)
poly_cmplx(-1 - i)  = (5.5,-2)
poly_cmplx( 1 - 2i) = (4.5,6)

poly_mat({ 1, 0.5, -1.4, 2.8 }) = Rows: 2  Cols: 2
[       0.55,     0.0667,  ]
[     -0.187,       0.79,  ]


( 3 * x^2 + 2 * x1 + 1 )  +  ( 4 * x^3 + 3 * x^2 + 2 * x1 + 1 ) = 4 * x^3 + 6 * x^2 + 4 * x1 + 2
( 3 * x^2 + 2 * x1 + 1 )  -  ( 4 * x^3 + 3 * x^2 + 2 * x1 + 1 ) = -4 * x^3
( 3 * x^2 + 2 * x1 + 1 )  *  ( 4 * x^3 + 3 * x^2 + 2 * x1 + 1 ) = 12 * x^5 + 17 * x^4 + 16 * x^3 + 10 * x^2 + 4 * x1 + 1
( 12 * x^5 + 17 * x^4 + 16 * x^3 + 10 * x^2 + 4 * x1 + 1 )  /  ( 4 * x^3 + 3 * x^2 + 2 * x1 + 1 ) = 2 * x1 + 1  remainder:

Real coef. polynom output:
0.4 * x^4 + -2 * x^3 + 1.3 * x^2 + 0.25 * x1 + -1
			 0.4 * x^4 +         -2 * x^3 +        1.3 * x^2 +       0.25 * x1 +         -1
		0.4 * x^4 +      -2 * x^3 +     1.3 * x^2 +    0.25 * x1 +      -1
Complex coef. polynom output:
(-2,3) * x1 + (0.5,-1)
		(-2,3) * x1 +   (0.5,-1)
 (-2,3) * x1 + (0.5,-1)
*/

	// Run additional demos
	Docs_Demo_Polynom_Calculus();
	Docs_Demo_Polynom_RootFinding();
}

// Polynom - Calculus operations demo
void Docs_Demo_Polynom_Calculus()
{
	std::cout << "\n*** Polynom Calculus Operations ***" << std::endl;

	// Create a polynomial: p(x) = x^3 - 3x^2 + 2x - 1
	PolynomReal p({-1, 2, -3, 1});  // coefficients from lowest to highest degree
	std::cout << "p(x)  = " << p << std::endl;

	// Derivative - returns a new polynomial
	PolynomReal dp = p.Derive();
	std::cout << "p'(x) = " << dp << std::endl;

	// Second derivative
	PolynomReal ddp = dp.Derive();
	std::cout << "p''(x)= " << ddp << std::endl;

	// Integration - returns a polynomial with constant term = 0
	PolynomReal int_p = p.Integrate();
	std::cout << "∫p(x) = " << int_p << " (+ C)" << std::endl;

	// Derive() with position: evaluates polynomial and derivatives at a point
	Vector<Real> derivs_at_2(4);  // [p(2), p'(2), p''(2), p'''(2)]
	p.Derive(2.0, derivs_at_2);
	std::cout << "At x = 2:" << std::endl;
	std::cout << "  p(2)    = " << derivs_at_2[0] << std::endl;
	std::cout << "  p'(2)   = " << derivs_at_2[1] << std::endl;
	std::cout << "  p''(2)  = " << derivs_at_2[2] << std::endl;
	std::cout << "  p'''(2) = " << derivs_at_2[3] << std::endl;
}

// Polynom - Root Finding demo
void Docs_Demo_Polynom_RootFinding()
{
	std::cout << "\n*** Polynom Root Finding ***" << std::endl;

	// Create a polynomial with known roots: (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6
	PolynomReal p({-6, 11, -6, 1});
	std::cout << "p(x) = " << p << std::endl;

	// Find roots using Laguerre's method
	Vector<Complex> roots = RootFinding::LaguerreRoots(p, 1e-10, 100);
	
	std::cout << "Roots found:" << std::endl;
	for (int i = 0; i < roots.size(); ++i) {
		std::cout << "  root " << i+1 << " = " << roots[i];
		if (std::abs(roots[i].imag()) < 1e-10) {
			std::cout << " (real: " << roots[i].real() << ")";
		}
		std::cout << std::endl;
	}

	// Verify roots by using complex polynomial evaluation
	PolynomComplex pc({-6, 11, -6, 1});  // same polynomial, but accepts Complex input
	std::cout << "Verification (should be ~0):" << std::endl;
	for (int i = 0; i < roots.size(); ++i) {
		Complex val = pc(roots[i]);
		std::cout << "  p(root " << i+1 << ") = " << val << std::endl;
	}
}