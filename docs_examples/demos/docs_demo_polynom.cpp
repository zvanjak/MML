#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Polynom.h"
#endif

using namespace MML;

void Docs_Demo_Polynom()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                           POLYNOM                             ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// Initialization of polynomials
	RealPolynom pol_constant({ 1 });
	RealPolynom pol_linear({ 1, 2 });					
	RealPolynom pol_quadratic({ 1, 2, 3 });			
	RealPolynom pol_cubic({ 1, 2, 3, 4 });			
	RealPolynom pol_quartic({ 1, 2, 3, 4, 5 });

	RealPolynom    poly_real({ -1, 0.25, 1.3, -2, 0.4 });
	ComplexPolynom poly_cmplx_r({ 1, 2, 3, 4 });
	ComplexPolynom poly_cmplx({ Complex(0.5,-1), Complex(-2,3) });
	Matrix2Polynom poly_mat({ 1, -1.0/2, 1./6 });            // matrix polynomial of 3rd order

	// Basic output of polynomials
	std::cout << "pol_constant  : " << pol_constant << std::endl;
	std::cout << "pol_linear    : " << pol_linear << std::endl;
	std::cout << "pol_quadratic : " << pol_quadratic << std::endl;
	std::cout << "pol_cubic     : " << pol_cubic << std::endl;
	std::cout << "pol_quartic   : " << pol_quartic << std::endl;
	std::cout << "poly_real     : " << poly_real << std::endl;
	std::cout << "poly_cmplx_r  : " << poly_cmplx_r << std::endl;
	std::cout << "poly_cmplx    : " << poly_cmplx << std::endl;
	std::cout << "m2            : " << poly_mat << std::endl << std::endl;

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
	RealPolynom pol_sum = pol_quadratic + pol_cubic;
	std::cout << "( " << pol_quadratic << " )  +  ( " << pol_cubic << " ) = " << pol_sum << std::endl;
	RealPolynom pol_diff = pol_quadratic - pol_cubic;
	std::cout << "( " << pol_quadratic << " )  -  ( " << pol_cubic << " ) = " << pol_diff << std::endl;
	RealPolynom pol_prod = pol_quadratic * pol_cubic;
	std::cout << "( " << pol_quadratic << " )  *  ( " << pol_cubic << " ) = " << pol_prod << std::endl;
	RealPolynom pol_div, pol_rem;
	RealPolynom::poldiv(pol_prod, pol_cubic, pol_div, pol_rem);
	std::cout << "( " << pol_prod << " )  /  ( " << pol_cubic << " ) = " << pol_div << "  remainder: " << pol_rem << std::endl;

	RealPolynom pol_sum2 = pol_quadratic * 2.0;
	RealPolynom pol_diff2 = pol_quadratic / 2.0;
	RealPolynom pol_prod2 = 2.0 * pol_quadratic;

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
m2            : 0.166667 * x^2 + -0.5 * x1 + 1

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
( 3 * x^2 + 2 * x1 + 1 )  *  ( 4 * x^3 + 3 * x^2 + 2 * x1 + 1 ) = 12 * x^5 + 17 * x^4 + 16 * x^3 + 10 * x^2 + 4 * x1 + 1( 12 * x^5 + 17 * x^4 + 16 * x^3 + 10 * x^2 + 4 * x1 + 1 )  /  ( 4 * x^3 + 3 * x^2 + 2 * x1 + 1 ) = 2 * x1 + 1  remainder:

Real coef. polynom output:
0.4 * x^4 + -2 * x^3 + 1.3 * x^2 + 0.25 * x1 + -1
			 0.4 * x^4 +         -2 * x^3 +        1.3 * x^2 +       0.25 * x1 +         -1
		0.4 * x^4 +      -2 * x^3 +     1.3 * x^2 +    0.25 * x1 +      -1
Complex coef. polynom output:
(-2,3) * x1 + (0.5,-1)
		(-2,3) * x1 +   (0.5,-1)
 (-2,3) * x1 + (0.5,-1)
*/
}