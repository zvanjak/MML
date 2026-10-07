///////////////////////////////////////////////////////////////////////////////////////////
// docs_demo_fundamentals.cpp - runnable counterpart of docs/Fundamentals.md
// Every snippet in that document must compile and execute here (docs<->demos contract).
///////////////////////////////////////////////////////////////////////////////////////////
#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>
#include <mml/MMLConcepts.h>
#include <mml/MMLExceptions.h>
#include <mml/base/Function.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/algorithms/RootFinding.h>
#endif

#include <functional>
#include <iostream>

using namespace MML;

namespace
{
	// §2 - concepts are compile-time facts; assert them right here
	static_assert(MMLScalar<Real> && MMLScalar<Complex>);
	static_assert(!MMLScalar<int>);
	static_assert(Field<Complex>);

	void Fundamentals_Section2_ConceptsAndContainers(std::ostream& out)
	{
		Matrix<Real>    A(3, 3);
		Matrix<Complex> C(2, 2);
		C(0, 0) = Complex(1, 1);
		C(1, 1) = Complex(0, -1);

		out << "  [2] Matrix<Real> 3x3 and Matrix<Complex> 2x2 constructed; C(0,0) = "
		    << C(0, 0) << "\n";
	}

	void Fundamentals_Section3_FunctionForms(std::ostream& out)
	{
		// 1) Concrete adapter over a plain function (or capture-less lambda)
		RealFunction f1([](Real x) { return x * x - 2; });

		// 2) Adapter over std::function - lambdas WITH captures
		Real shift = 2.0;
		RealFunctionFromStdFunc f2(std::function<Real(Real)>(
			[shift](Real x) { return x * x - shift; }));

		// 3) Raw lambda straight into an algorithm (RealFunctionCallable overload);
		//    simple overloads return the answer directly
		Real root = RootFinding::FindRootBisection([](Real x) { return x * x - 2; },
		                                           0.0, 2.0, 1e-12);

		out << "  [3] f1(2) = " << f1(2.0) << ", f2(2) = " << f2(2.0)
		    << ", bisection root of x^2-2 via raw lambda = " << root << "\n";
	}

	void Fundamentals_Section4_ConfigResult(std::ostream& out)
	{
		RealFunction f1([](Real x) { return x * x - 2; });

		// Simple call - plain answer, sensible defaults
		Real simpleRoot = RootFinding::FindRootBrent(f1, 0.0, 2.0, 1e-10);

		// Full control - a Config struct
		RootFinding::RootFindingConfig config;
		config.tolerance      = 1e-14;
		config.max_iterations = 200;
		auto result = RootFinding::FindRootBrent(f1, 0.0, 2.0, config);

		if (result.converged) {
			out << "  [4] sqrt(2) = " << result.root
			    << "  (" << result.algorithm_name
			    << ", iterations: " << result.iterations_used
			    << ", f(root): " << result.function_value
			    << ", f-evals: " << result.function_evaluations << ")\n";
		} else {
			out << "  [4] FAILED: status = " << ToString(result.status)
			    << ", " << result.error_message << "\n";
		}
		out << "  [4] simple overload agrees: " << simpleRoot << "\n";
	}

	void Fundamentals_Section5_Contracts(std::ostream& out)
	{
		Vector<Real> v({ 1e-14, -1e-14, 0.0 });
		bool exact = v.isZero();          // false - exact query
		bool near  = v.isNearZero();      // true  - tolerance query (Defaults-scaled)

		out << "  [5] isZero: " << std::boolalpha << exact
		    << ", isNearZero: " << near << "\n";

		try {
			Real x = v.at(17);            // checked access - throws
			(void)x;
			out << "  [5] ERROR: at(17) did not throw!\n";
		} catch (const MMLException& e) { // one base class catches any MML error
			out << "  [5] at(17) threw as contracted: " << e.message() << "\n";
		}
	}
} // namespace

void Docs_Demo_Fundamentals()
{
	std::ostream& out = std::cout;
	out << "=== Docs_Demo_Fundamentals (docs/Fundamentals.md) ===\n";
	Fundamentals_Section2_ConceptsAndContainers(out);
	Fundamentals_Section3_FunctionForms(out);
	Fundamentals_Section4_ConfigResult(out);
	Fundamentals_Section5_Contracts(out);
	out << "=== Fundamentals demo done ===\n";
}
