#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Function.h"
#include "core/FunctionHelpers.h"
#include "core/Derivation.h"
#include "core/Integration.h"

#include "tools/ConsolePrinter.h"
#endif


using namespace MML;

void Readme_deriving_functions()
{
	RealFunction       f1{ [](Real x) { return (Real)(sin(x) * (1.0 + 0.5 * x * x)); } };

	// numerical derivation of real function (available orders - 1, 2, 4, 6, 8)
	double der_f1 = Derivation::NDer1(f1, 0.5);
	double der_f4 = Derivation::NDer2(f1, 0.5, 1e-6);   // setting explicit step size
	Real err;
	double der_f6 = Derivation::NDer6(f1, 0.5, &err);   // if we need error estimate    
	// we can use default Derive routine (set to NDer4)
	double num_der4 = Derivation::Derive(f1, 0.5);

	// second and third derivatives
	double sec_der_f1 = Derivation::NSecDer2(f1, 0.5);
	double third_der_f1 = Derivation::NThirdDer2(f1, 0.5);

	// creating new function that is derivation of existing function
	RealFuncDerived4    f1_der4(f1);        // 4th order derivation

	// scalar and vector functions
	ScalarFunction<3>   f2Scal([](const VectorN<Real, 3>& x) { return (Real)(1.0 / pow(x.NormL2(), 2)); });
	VectorFunction<3>   f3Vec([](const VectorN<Real, 3>& x) { return VectorN<Real, 3>{0, x[0] * x[1], 0}; });
	VectorN<Real, 3>    der_point{ 1.0, 1.0, 1.0 };

	double der_f2 = Derivation::NDer1Partial(f2Scal, 1, der_point);
	VectorN<Real, 3> der_f2_all = Derivation::NDer1PartialByAll(f2Scal, der_point);

	double der_f3 = Derivation::NDer1Partial(f3Vec, 1, 1, der_point);
	VectorN<Real, 3>     der_f3_by1 = Derivation::NDer2PartialByAll(f3Vec, 1, der_point);
	MatrixNM<Real, 3, 3> der_f3_by_all = Derivation::NDer4PartialAllByAll(f3Vec, der_point);
}

void Readme_Derivation_precision()
{
	// our test function
	RealFunction       f{ [](Real x) { return (Real)(sin(x) * (1.0 + 0.5 * x * x)); } };
	// and its exact derivation
	RealFunction       f_der{ [](Real x) { return (Real)(cos(x) * (1.0 + 0.5 * x * x) + sin(x) * x); } };

	double x1 = -7.0, x2 = 7.0;
	double err_sum1 = 0.0, err_sum2 = 0.0, err_sum4 = 0.0, err_sum6 = 0.0, err_sum8 = 0.0;
	int    numPntForEval = 20;

	std::cout << "\nAVERAGE DERIVATION ERROR FOR DIFFERENT ORDERS" << std::endl;

	TablePrinter<double, double> print_data("x", 8, 3,
		{ "Exact der.",
			"Nder1", "Nder1 err.", "Nder2", "Nder2 err.", "Nder4", "Nder4 err.", "Nder6", "Nder6 err.", "Nder8", "Nder8 err."
		},
		{ {12,7,'F'},
			{13,7,'F'}, {15,6,'S'}, {13,7,'F'}, {15,6,'S'}, {13,7,'F'}, {15,6,'S'}, {13,7,'F'}, {15,6,'S'}, {13,7,'F'}, {15,6,'S'}
		}
	);

	for (int i = 0; i < numPntForEval; i++) {
		double x = x1 + (x2 - x1) * i / (numPntForEval - 1);

		double exact_der = f_der(x);

		double num_der1 = MML::Derivation::NDer1(f, x);
		double num_der2 = MML::Derivation::NDer2(f, x);
		double num_der4 = MML::Derivation::NDer4(f, x);
		double num_der6 = MML::Derivation::NDer6(f, x);
		double num_der8 = MML::Derivation::NDer8(f, x);

		double err1 = num_der1 - exact_der;
		double err2 = num_der2 - exact_der;
		double err4 = num_der4 - exact_der;
		double err6 = num_der6 - exact_der;
		double err8 = num_der8 - exact_der;

		print_data.addRow(x, { exact_der, num_der1, err1, num_der2, err2, num_der4, err4, num_der6, err6, num_der8, err8 });

		err_sum1 += std::abs(err1);
		err_sum2 += std::abs(err2);
		err_sum4 += std::abs(err4);
		err_sum6 += std::abs(err6);
		err_sum8 += std::abs(err8);
	}
	print_data.Print();

	std::cout << "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << "Total abs error =                      " << std::scientific << err_sum1 << "                  "
		<< err_sum2 << "                  "
		<< err_sum4 << "                  "
		<< err_sum6 << "                  "
		<< err_sum8 << std::endl;

	/* OUTPUT
 AVERAGE DERIVATION ERROR FOR DIFFERENT ORDERS
			 x   Exact der.         Nder1      Nder1 err.         Nder2      Nder2 err.         Nder4      Nder4 err.         Nder6      Nder6 err.         Nder8      Nder8 err.
	-7.000   23.8234137    23.8234137    5.288187e-08    23.8234137   -8.782806e-10    23.8234137   -1.159250e-11    23.8234137    3.481659e-13    23.8234137   -8.881784e-14
	-6.263   20.4840131    20.4840129   -1.939339e-07    20.4840131   -3.526814e-10    20.4840131   -8.846257e-12    20.4840131    4.369838e-13    20.4840131    1.634248e-13
	-5.526    8.0335341     8.0335339   -2.569286e-07     8.0335341    2.451923e-10     8.0335341   -2.472689e-12     8.0335341    3.339551e-13     8.0335341    2.540190e-13
	-4.789   -3.8149928    -3.8149929   -1.382276e-07    -3.8149928    2.549019e-10    -3.8149928    3.900436e-12    -3.8149928   -4.605205e-13    -3.8149928    2.198242e-13
	-4.053   -8.8483626    -8.8483626   -6.771076e-09    -8.8483626    1.632490e-10    -8.8483626    4.400036e-12    -8.8483626    1.030287e-13    -8.8483626    1.598721e-14
	-3.316   -6.9735845    -6.9735844    8.421828e-08    -6.9735845    1.058416e-10    -6.9735845   -1.277201e-12    -6.9735845   -2.193801e-13    -6.9735845   -1.927347e-13
	-2.579   -2.2830219    -2.2830218    9.390760e-08    -2.2830219   -4.583889e-11    -2.2830219   -1.086242e-12    -2.2830219   -9.459100e-14    -2.2830219   -2.216005e-13
	-1.842    1.0520333     1.0520333    4.414189e-08     1.0520333   -8.739320e-11     1.0520333   -5.551115e-13     1.0520333   -1.398881e-13     1.0520333   -1.421085e-13
	-1.105    1.7107321     1.7107321   -6.558459e-09     1.7107321   -5.454659e-11     1.7107321   -3.994582e-13     1.7107321    5.573320e-14     1.7107321    3.264056e-14
	-0.368    1.1288943     1.1288943   -8.092415e-09     1.1288943    1.287215e-11     1.1288943    5.462297e-13     1.1288943    6.439294e-15     1.1288943    1.654232e-13
	 0.368    1.1288943     1.1288944    1.053404e-08     1.1288943    1.287215e-11     1.1288943    5.462297e-13     1.1288943    6.439294e-15     1.1288943    1.654232e-13
	 1.105    1.7107321     1.7107321    1.579328e-08     1.7107321   -4.183387e-11     1.7107321   -3.990142e-13     1.7107321    1.276756e-13     1.7107321    3.730349e-14
	 1.842    1.0520333     1.0520332   -3.036391e-08     1.0520333   -8.739098e-11     1.0520333   -1.045830e-12     1.0520333   -2.109424e-14     1.0520333   -1.376677e-13
	 2.579   -2.2830219    -2.2830220   -7.000517e-08    -2.2830219    5.014655e-12    -2.2830219   -8.353318e-13    -2.2830219   -6.616929e-14    -2.2830219   -2.358114e-13
	 3.316   -6.9735845    -6.9735846   -8.714507e-08    -6.9735845    1.058416e-10    -6.9735845   -1.277201e-12    -6.9735845   -2.193801e-13    -6.9735845   -1.927347e-13
	 4.053   -8.8483626    -8.8483625    8.263589e-08    -8.8483626    1.632490e-10    -8.8483626    4.400036e-12    -8.8483626    1.030287e-13    -8.8483626    1.598721e-14
	 4.789   -3.8149928    -3.8149926    1.597957e-07    -3.8149928    2.549019e-10    -3.8149928    3.900436e-12    -3.8149928   -4.605205e-13    -3.8149928    2.198242e-13
	 5.526    8.0335341     8.0335344    2.795132e-07     8.0335341    4.177991e-11     8.0335341   -3.844036e-12     8.0335341    2.238210e-13     8.0335341    3.073097e-13
	 6.263   20.4840131    20.4840133    1.934963e-07    20.4840131   -3.526814e-10    20.4840131   -8.846257e-12    20.4840131    4.369838e-13    20.4840131    1.634248e-13
	 7.000   23.8234137    23.8234136   -6.632742e-08    23.8234137   -8.782806e-10    23.8234137   -1.159250e-11    23.8234137    3.481659e-13    23.8234137   -8.881784e-14
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Total abs error =                      1.881272e-06                  4.144644e-09                  7.176304e-11                  4.211964e-12                  3.060885e-12
*/
}


void Readme_derivation()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                       README - derivation                     ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Readme_deriving_functions();
	Readme_Derivation_precision();
}