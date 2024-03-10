#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Function.h"
#include "core/Integration.h"
#endif

#include "../test_data/real_functions_test_bed.h"
#include "../test_data/scalar_functions_test_bed.h"
#include "../test_data/vector_functions_test_bed.h"

using namespace MML;

void Test_Precision_Integration_Single_Func()
{
	const TestBeds::TestFunctionRealWithIntegral& f_wrap = TestBeds::RealFunctionsTestBed::getTestFunctionRealWithIntegral(0);

	const RealFunction& f = f_wrap._func;
	const RealFunction& f_int = f_wrap._funcIntegrated;
	double x1 = f_wrap._intervalTest->getLowerBound();
	double x2 = f_wrap._intervalTest->getUpperBound();

	const int numIntervals = 20;

	double err_sum1 = 0.0;
	double err_sum2 = 0.0;
	double err_sum3 = 0.0;

	std::cout << "\nAVERAGE INTEGRATION ERROR FOR DIFFERENT INTEGRATORS" << std::endl;
	std::cout << "   Interval        Exact int.         Trap         Trap err.           Simpson      Simpson  err.       Romberg       Romberg err.          " << std::endl;
	std::cout << "--------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;

	for (int i = 1; i < numIntervals; i++) {
		double x = x1 + (x2 - x1) * i / (numIntervals - 1);

		double integral = f_int(x) - f_int(x1);

		double int_trap = MML::IntegrateTrap(f, x1, x, 1e-3);
		double int_simp = MML::IntegrateSimpson(f, x1, x, 1e-3);
		double int_romb = MML::IntegrateRomberg(f, x1, x);

		double err1 = int_trap - integral;
		double err2 = int_simp - integral;
		double err3 = int_romb - integral;

		std::cout << "[" << std::fixed
			<< std::setw(6) << std::setprecision(3) << x1 << ", "
			<< std::setw(6) << std::setprecision(3) << x << "]  "
			<< std::setw(13) << std::setprecision(8) << integral << " "
			<< std::setw(13) << std::setprecision(8) << int_trap << "   "
			<< std::scientific << std::setw(15) << err1 << "   " << std::fixed
			<< std::setw(13) << std::setprecision(8) << int_simp << "   "
			<< std::scientific << std::setw(15) << err2 << "   " << std::fixed
			<< std::setw(13) << std::setprecision(8) << int_romb << "   "
			<< std::scientific << std::setw(15) << err3 << "   " << std::fixed
			<< std::endl;

		err_sum1 += std::abs(err1);
		err_sum2 += std::abs(err2);
		err_sum3 += std::abs(err3);
	}

	std::cout << "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << "Total abs error =                    " << std::scientific << err_sum1 << "                    "
		<< err_sum2 << "                    "
		<< err_sum3 << std::endl;
}

void Test_Precision_Integration()
{
	Test_Precision_Integration_Single_Func();
}