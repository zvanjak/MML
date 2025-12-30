#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Polynom.h"
#endif

using namespace MML;

// TODO - evaluate exp(x) series for matrix and compare with direct calculation
void Demo_Polynom()
{
	std::cout << std::endl;
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                           POLYNOM                             ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// Initialization of polynomials
	PolynomRealFunc pol_constant({ -2.5 });
	PolynomRealFunc pol_linear({ -1, 2 });
	PolynomRealFunc pol_quadratic({ -2, -1, -1 });
	PolynomRealFunc pol_cubic({ -1, 2, 3, -4 });
	PolynomRealFunc pol_quartic({ 1, -2, 3, 4, 1 });

	// I/O of polynomials
	std::cout << "Polynomial initialization:\n";
	std::cout << "pol_constant  : " << pol_constant << std::endl;
	std::cout << "pol_linear    : " << pol_linear << std::endl;
	std::cout << "pol_quadratic : " << pol_quadratic << std::endl;
	std::cout << "pol_cubic     : " << pol_cubic << std::endl;
	std::cout << "pol_quartic   : " << pol_quartic << std::endl;

	PolynomRealFunc		 poly_real({ 1, 2, 3, 4, 5 });
	PolynomComplex poly_cmplx({ 2, -1, 3 });
	Matrix2Polynom poly_mat({ 1, -1, 2 });            // matrix polynomial of 3rd order

	// Evaluation of polynomials
	double  v = poly_real(5.0);
	Complex q = poly_cmplx(Complex(2, 3));
	MatrixNM<Real, 2, 2> m2_3 = poly_mat(MatrixNM<Real, 2, 2>({ 1, -1, 0.5, 2 }));     // evaluate the polynomial at the given matrix value

	// in detail
	MatrixNM<Real, 2, 2> eval_mat({ 1, -1, 0.5, 2 });          // matrix to evaluate the polynomial at
	MatrixNM<Real, 2, 2> m2_2 = poly_mat(eval_mat);       // evaluate the polynomial at the given matrix value

	// Operations on polynomials
	PolynomRealFunc pol_sum = pol_quadratic + pol_cubic;
	std::cout << "\nAdd : (" << pol_quadratic << ")  +  (" << pol_cubic << ") = " << pol_quadratic + pol_cubic << std::endl;

	PolynomRealFunc pol_diff = pol_quadratic - pol_cubic;
	std::cout << "Sub : (" << pol_quadratic << ")  -  (" << pol_cubic << ") = " << pol_quadratic - pol_cubic << std::endl;

	PolynomRealFunc pol_prod = pol_quadratic * pol_cubic;
	std::cout << "Mul : (" << pol_quadratic << ")  *  (" << pol_cubic << ") = " << pol_quadratic * pol_cubic << std::endl;

	PolynomRealFunc pol_sum2 = pol_quadratic * 2.0;
	PolynomRealFunc pol_diff2 = pol_quadratic / 2.0;
	PolynomRealFunc pol_prod2 = 2.0 * pol_quadratic;

	PolynomRealFunc pol_quot, pol_rem;
	PolynomRealFunc::poldiv(pol_quartic, pol_quadratic, pol_quot, pol_rem);
	std::cout << "Div : (" << pol_quartic << ")  /  (" << pol_quadratic << ") = " << pol_quot << "  remainder: " << pol_rem << std::endl;

	PolynomRealFunc pol_calc = pol_quadratic * pol_quot + pol_rem;
	std::cout << "Div check : (" << pol_quadratic << ")  *  (" << pol_quot << ") + (" << pol_rem << ") = " << pol_quadratic * pol_quot + pol_rem << std::endl;

	// create polynom respresenting sinus and cosinus functions up to x^5
	PolynomRealFunc pol_sin5({1./120,     0, -1./6,     0, 1, 0});
	PolynomRealFunc pol_cos5({     0, 1./24,     0, -1./2, 0, 1});

	PolynomRealFunc pol_sin7({-1./5040,    0   , 1./120,   0  , -1./6,   0  , 1, 0});
	PolynomRealFunc pol_cos7({    0   , -1./720,    0  , 1./24,   0  , -1./2, 0, 1});

	std::cout << "\nPolynom output:\n";
	std::cout << "std::cout << poly_real                  => " << poly_real << std::endl;
	std::cout << "std::cout << poly_real.to_string(10, 5) => " << poly_real.to_string(10, 5) << std::endl;
	std::cout << "poly_real.Print(std::cout, 7, 3);       => ";	poly_real.Print(std::cout, 7, 3); std::cout << std::endl;
	std::cout << "poly_real.Print(std::cout);             => "; poly_real.Print(std::cout); std::cout << std::endl;
}
