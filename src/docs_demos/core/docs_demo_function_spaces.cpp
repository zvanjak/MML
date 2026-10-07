///////////////////////////////////////////////////////////////////////////////////////////
// docs_demo_function_spaces.cpp - Demonstrations for Function_spaces.md
//
// This file contains deterministic examples for the finite FunctionSpaces layer:
//   - L2 interval descriptors and weighted inner products
//   - modal projection into an orthogonal basis
//   - Chebyshev-Lobatto collocation nodes and derivative matrices
//   - dense BVP solve for a Poisson benchmark
//   - matrix-free finite-difference application
///////////////////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <iomanip>
#include <iostream>

#include <mml/MMLBase.h>
#include <mml/core/FunctionSpaces.h>
#include <mml/core/OrthogonalBasis/LegendreBasis.h>

using namespace MML;
using namespace MML::FunctionSpaces;

namespace
{
	class DemoConstantFunction : public IRealFunction
	{
		Real _value;

	public:
		explicit DemoConstantFunction(Real value) : _value(value) { }
		Real operator()(Real) const override { return _value; }
	};

	class DemoLegendrePolynomial : public IRealFunction
	{
	public:
		Real operator()(Real x) const override
		{
			Real p2 = (REAL(3.0) * x * x - REAL(1.0)) / REAL(2.0);
			return REAL(1.0) + REAL(2.0) * x + REAL(3.0) * p2;
		}
	};

	class DemoPoissonRhs : public IRealFunction
	{
	public:
		Real operator()(Real x) const override
		{
			return Constants::PI * Constants::PI * std::sin(Constants::PI * x);
		}
	};
}

void Docs_Demo_FunctionSpaces()
{
	std::cout << "========================================\n";
	std::cout << "   FUNCTION SPACES DEMO\n";
	std::cout << "========================================\n\n";

	std::cout << "--- L2 interval and inner product ---\n";
	L2IntervalSpace l2(-REAL(1.0), REAL(1.0));
	DemoConstantFunction one(REAL(1.0));
	std::cout << "Domain: [" << l2.domainMin() << ", " << l2.domainMax() << "]\n";
	std::cout << "<1,1> = " << l2.innerProduct(one, one) << "\n\n";

	std::cout << "--- L2 projection into Legendre basis ---\n";
	LegendreBasis legendre;
	OrthogonalBasisTrialSpace1D legendreSpace(legendre, 3);
	DemoLegendrePolynomial polynomial;
	FunctionExpansion1D projection = ProjectL2(polynomial, legendreSpace);
	for (int i = 0; i < projection.dimension(); ++i)
		std::cout << "c[" << i << "] = " << std::setprecision(10) << projection.coefficients()[i] << "\n";
	std::cout << "u_N(0.25) = " << projection.evaluate(REAL(0.25)) << "\n\n";

	std::cout << "--- Chebyshev-Lobatto collocation ---\n";
	ChebyshevCollocationSpace1D cheb(-REAL(1.0), REAL(1.0), 8);
	std::cout << "Nodes:";
	for (int i = 0; i < cheb.nodeCount(); ++i)
		std::cout << " " << std::setprecision(5) << cheb.node(i);
	std::cout << "\n";
	std::cout << "D1(0,0) = " << cheb.firstDerivativeMatrix()(0, 0) << "\n\n";

	std::cout << "--- Dense BVP solve: -u'' = pi^2 sin(pi x) ---\n";
	DemoPoissonRhs rhs;
	auto L = LinearDifferentialOperator1D::SecondOrder(
		[](Real) { return -REAL(1.0); },
		[](Real) { return REAL(0.0); },
		[](Real) { return REAL(0.0); });
	BoundaryConditions1D bc{
		BoundaryCondition1D::Dirichlet(-REAL(1.0), REAL(0.0)),
		BoundaryCondition1D::Dirichlet(REAL(1.0), REAL(0.0))
	};
	auto solve = SolveDenseCollocationBVP(L, rhs, cheb, bc);
	std::cout << "Converged: " << (solve.success() ? "yes" : "no") << "\n";
	if (solve.success())
		std::cout << "u_N(0.5) = " << solve.solution->evaluate(REAL(0.5)) << "\n\n";

	std::cout << "--- Matrix-free finite-difference operator ---\n";
	DirichletSecondDerivativeOperator1D d2(3, REAL(0.25));
	Vector<Real> values{ REAL(1.0), REAL(2.0), REAL(1.0) };
	Vector<Real> applied = d2.apply(values);
	std::cout << "D2*[1,2,1] = [" << applied[0] << ", " << applied[1] << ", " << applied[2] << "]\n\n";
}
