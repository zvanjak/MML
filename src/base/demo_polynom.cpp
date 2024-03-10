#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Polynom.h"
#include "base/LinearFunctional.h"
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
	RealPolynom pol_constant({ 1 });
	RealPolynom pol_linear({ 1, 2 });
	RealPolynom pol_quadratic({ 1, 2, 3 });
	RealPolynom pol_cubic({ 1, 2, 3, 4 });
	RealPolynom pol_quartic({ 1, 2, 3, 4, 5 });

	ComplexPolynom p2({ 1, 2, 3, 4, 5 });

	RealPolynom poly_real({ 1, 2, 3, 4, 5 });
	ComplexPolynom poly_cmplx({ Complex(1,1), Complex(1,1) });
	Matrix2Polynom poly_mat({ 1, 2, 3 });            // matrix polynomial of 3rd order

	// I/O of polynomials
	std::cout << "pol_constant  : " << pol_constant << std::endl;
	std::cout << "pol_linear    : " << pol_linear << std::endl;
	std::cout << "pol_quadratic : " << pol_quadratic << std::endl;
	std::cout << "pol_cubic     : " << pol_cubic << std::endl;
	std::cout << "pol_quartic   : " << pol_quartic << std::endl;
	std::cout << "p2            : " << p2 << std::endl;
	std::cout << "poly_real     : " << poly_real << std::endl;
	std::cout << "poly_cmplx    : " << poly_cmplx << std::endl;
	std::cout << "m2            : " << poly_mat << std::endl;

	// Evaluation of polynomials
	double  v = poly_real(5.0);
	Complex q = poly_cmplx(Complex(2, 3));
	MatrixNM<Real, 2, 2> m2_3 = poly_mat(MatrixNM<Real, 2, 2>({ 1, 2, 3, 4 }));     // evaluate the polynomial at the given matrix value

	// in detail
	MatrixNM<Real, 2, 2> eval_mat({ 1, 2, 3, 4 });          // matrix to evaluate the polynomial at
	MatrixNM<Real, 2, 2> m2_2 = poly_mat(eval_mat);       // evaluate the polynomial at the given matrix value

	// Operations on polynomials
	RealPolynom pol_sum = pol_quadratic + pol_cubic;
	RealPolynom pol_diff = pol_quadratic - pol_cubic;
	RealPolynom pol_prod = pol_quadratic * pol_cubic;

	RealPolynom pol_sum2 = pol_quadratic * 2.0;
	RealPolynom pol_diff2 = pol_quadratic / 2.0;
	RealPolynom pol_prod2 = 2.0 * pol_quadratic;

	std::cout << "\nReal polynom output:\n";
	std::cout << poly_real << std::endl;
	std::cout << poly_real.to_string(10, 5) << std::endl;
	poly_real.Print(std::cout, 7, 3);

	std::cout << "\nComplex polynom output:\n";
	std::cout << poly_cmplx << std::endl;
	std::cout << poly_cmplx.to_string(10, 5) << std::endl;
	poly_cmplx.Print(std::cout, 7, 3);

	std::cout << "\nReal matrix output:\n";
	std::cout << poly_mat << std::endl;
	std::cout << poly_mat.to_string(10, 5) << std::endl;
	poly_mat.Print(std::cout, 7, 3);
}
