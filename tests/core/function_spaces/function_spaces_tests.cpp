#include <catch2/catch_all.hpp>

#include "../../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/core/FunctionSpaces.h>
#endif

#include <limits>
#include <type_traits>

#ifndef MML_USE_SINGLE_HEADER
#include <mml/core/OrthogonalBasis/LegendreBasis.h>
#endif

using namespace MML;
using namespace MML::FunctionSpaces;

namespace MML::Tests::Core::FunctionSpacesTests
{
	constexpr Real QuadratureTolerance = Testing::Tol(REAL(1e-8), REAL(3e-5));
	constexpr Real ProjectionTolerance = Testing::Tol(REAL(1e-8), REAL(3e-5));
	constexpr Real BvpTolerance = Testing::Tol(REAL(1e-8), REAL(5e-5));
	constexpr Real BvpResidualTolerance = Testing::Tol(REAL(1e-8), REAL(3e-2));
	constexpr Real NodeTolerance = Testing::Tol(REAL(1e-14), REAL(1e-6));
	constexpr Real FirstDerivativeTolerance = Testing::Tol(REAL(1e-9), REAL(2e-4));
	constexpr Real SecondDerivativeTolerance = Testing::Tol(REAL(1e-8), REAL(3e-4));
	constexpr Real SmoothDerivativeTolerance = Testing::Tol(REAL(1e-10), REAL(4e-5));

	class ConstantFunction : public IRealFunction
	{
		Real _value;

	public:
		explicit ConstantFunction(Real value) : _value(value) { }

		Real operator()(Real) const override { return _value; }
	};

	class LinearFunction : public IRealFunction
	{
	public:
		Real operator()(Real x) const override { return x; }
	};

	class LegendreCombinationFunction : public IRealFunction
	{
	public:
		Real operator()(Real x) const override
		{
			Real p2 = (REAL(3.0) * x * x - REAL(1.0)) / REAL(2.0);
			return REAL(3.0) + REAL(2.0) * x + REAL(5.0) * p2;
		}
	};

	class PoissonSineRhs : public IRealFunction
	{
	public:
		Real operator()(Real x) const override
		{
			return Constants::PI * Constants::PI * std::sin(Constants::PI * x);
		}
	};

	class MonomialTrialSpace3 : public TrialSpace1D
	{
		const FunctionSpace1D& _space;

	public:
		explicit MonomialTrialSpace3(const FunctionSpace1D& space) : _space(space) { }

		const FunctionSpace1D& functionSpace() const noexcept override { return _space; }
		int dimension() const noexcept override { return 3; }

		Real basisValue(int basisIndex, Real x) const override
		{
			validateBasisIndex(basisIndex, "MonomialTrialSpace3::basisValue");
			if (basisIndex == 0)
				return REAL(1.0);
			if (basisIndex == 1)
				return x;
			return x * x;
		}

		bool hasBasisDerivative(int order = 1) const noexcept override
		{
			return order >= 0 && order <= 2;
		}

		Real basisDerivative(int basisIndex, int order, Real x) const override
		{
			validateBasisIndex(basisIndex, "MonomialTrialSpace3::basisDerivative");
			validateDerivativeOrder(order, "MonomialTrialSpace3::basisDerivative");
			if (order == 0)
				return basisValue(basisIndex, x);
			if (order == 1) {
				if (basisIndex == 0)
					return REAL(0.0);
				if (basisIndex == 1)
					return REAL(1.0);
				return REAL(2.0) * x;
			}
			if (basisIndex == 2)
				return REAL(2.0);
			return REAL(0.0);
		}

		bool hasNodes() const noexcept override { return true; }
		int nodeCount() const noexcept override { return 3; }

		Real node(int nodeIndex) const override
		{
			validateNodeIndex(nodeIndex, "MonomialTrialSpace3::node");
			return functionSpace().domainMin() + (functionSpace().domainLength() * nodeIndex) / REAL(2.0);
		}
	};

	class IdentityLinearOperator3 : public LinearOperator
	{
	public:
		using LinearOperator::apply;

		int rows() const noexcept override { return 3; }
		int cols() const noexcept override { return 3; }

		void apply(const Vector<Real>& input, Vector<Real>& output) const override
		{
			validateInputOutput(input, output, "IdentityLinearOperator3::apply");
			for (int i = 0; i < 3; ++i)
				output[i] = input[i];
		}
	};

	TEST_CASE("FunctionSpaces - Aggregate header exposes base declarations", "[FunctionSpaces][Core]")
	{
		static_assert(std::is_enum_v<DiscretizationMethod>);

		DiscretizationMethod method = DiscretizationMethod::Collocation;
		REQUIRE(method == DiscretizationMethod::Collocation);
	}

	TEST_CASE("FunctionSpace1D - L2 interval reports domain and unit weight", "[FunctionSpaces][FunctionSpace1D]")
	{
		L2IntervalSpace space(-REAL(2.0), REAL(3.0));

		REQUIRE(space.domainMin() == -REAL(2.0));
		REQUIRE(space.domainMax() == REAL(3.0));
		REQUIRE(space.domainLength() == REAL(5.0));
		REQUIRE(space.contains(REAL(0.0)));
		REQUIRE(space.contains(-REAL(2.0)));
		REQUIRE(space.contains(REAL(3.0)));
		REQUIRE_FALSE(space.contains(REAL(4.0)));
		REQUIRE(space.weight(REAL(0.25)) == REAL(1.0));
	}

	TEST_CASE("FunctionSpace1D - Weighted interval evaluates custom weight", "[FunctionSpaces][FunctionSpace1D]")
	{
		WeightedL2IntervalSpace space(REAL(0.0), REAL(1.0), [](Real x) { return REAL(1.0) + x; });

		REQUIRE(space.weight(REAL(0.0)) == REAL(1.0));
		REQUIRE(space.weight(REAL(0.5)) == REAL(1.5));
		REQUIRE(space.weight(REAL(1.0)) == REAL(2.0));
	}

	TEST_CASE("FunctionSpace1D - Inner product uses interval weight", "[FunctionSpaces][FunctionSpace1D]")
	{
		WeightedL2IntervalSpace space(REAL(0.0), REAL(1.0), [](Real x) { return REAL(1.0) + x; });
		ConstantFunction one(REAL(1.0));
		LinearFunction linear;

		REQUIRE(space.innerProduct(one, one) == Catch::Approx(REAL(1.5)).epsilon(QuadratureTolerance));
		REQUIRE(space.innerProduct(linear, one) == Catch::Approx(REAL(5.0) / REAL(6.0)).epsilon(QuadratureTolerance));
	}

	TEST_CASE("FunctionSpace1D - Invalid descriptors throw", "[FunctionSpaces][FunctionSpace1D]")
	{
		REQUIRE_THROWS_AS(L2IntervalSpace(REAL(1.0), REAL(1.0)), FunctionSpaceInputError);
		REQUIRE_THROWS_AS(L2IntervalSpace(REAL(2.0), REAL(1.0)), FunctionSpaceInputError);
		REQUIRE_THROWS_AS(WeightedL2IntervalSpace(REAL(0.0), REAL(1.0), nullptr), FunctionSpaceInputError);

		WeightedL2IntervalSpace badWeight(REAL(0.0), REAL(1.0), [](Real) { return std::numeric_limits<Real>::quiet_NaN(); });
		REQUIRE_THROWS_AS(badWeight.weight(REAL(0.5)), FunctionSpaceInputError);
	}

	TEST_CASE("TrialSpace1D - Exposes basis derivative and node metadata", "[FunctionSpaces][TrialSpace1D]")
	{
		L2IntervalSpace space(-REAL(1.0), REAL(1.0));
		MonomialTrialSpace3 trial(space);

		REQUIRE(&trial.functionSpace() == &space);
		REQUIRE(trial.dimension() == 3);
		REQUIRE(trial.basisValue(0, REAL(0.25)) == REAL(1.0));
		REQUIRE(trial.basisValue(1, REAL(0.25)) == REAL(0.25));
		REQUIRE(trial.basisValue(2, REAL(0.25)) == REAL(0.0625));

		REQUIRE(trial.hasBasisDerivative(1));
		REQUIRE(trial.hasBasisDerivative(2));
		REQUIRE_FALSE(trial.hasBasisDerivative(3));
		REQUIRE(trial.basisDerivative(2, 1, REAL(0.25)) == REAL(0.5));
		REQUIRE(trial.basisDerivative(2, 2, REAL(0.25)) == REAL(2.0));

		REQUIRE(trial.hasNodes());
		REQUIRE(trial.nodeCount() == 3);
		REQUIRE(trial.node(0) == -REAL(1.0));
		REQUIRE(trial.node(1) == REAL(0.0));
		REQUIRE(trial.node(2) == REAL(1.0));

		REQUIRE_THROWS_AS(trial.basisValue(3, REAL(0.0)), FunctionSpaceInputError);
		REQUIRE_THROWS_AS(trial.basisDerivative(0, -1, REAL(0.0)), FunctionSpaceInputError);
		REQUIRE_THROWS_AS(trial.node(3), FunctionSpaceInputError);
	}

	TEST_CASE("FunctionExpansion1D - Evaluates coefficient expansion and derivatives", "[FunctionSpaces][FunctionExpansion1D]")
	{
		L2IntervalSpace space(-REAL(1.0), REAL(1.0));
		MonomialTrialSpace3 trial(space);
		FunctionExpansion1D expansion(trial, Vector<Real>{ REAL(1.0), REAL(2.0), REAL(3.0) });

		REQUIRE(&expansion.trialSpace() == &trial);
		REQUIRE(&expansion.functionSpace() == &space);
		REQUIRE(expansion.dimension() == 3);
		REQUIRE(expansion.coefficients().size() == 3);
		REQUIRE(expansion.evaluate(REAL(0.5)) == REAL(2.75));
		REQUIRE(expansion.derivative(REAL(0.5), 0) == REAL(2.75));
		REQUIRE(expansion.derivative(REAL(0.5), 1) == REAL(5.0));
		REQUIRE(expansion.derivative(REAL(0.5), 2) == REAL(6.0));
	}

	TEST_CASE("FunctionExpansion1D - Rejects invalid coefficient and derivative requests", "[FunctionSpaces][FunctionExpansion1D]")
	{
		L2IntervalSpace space(-REAL(1.0), REAL(1.0));
		MonomialTrialSpace3 trial(space);

		REQUIRE_THROWS_AS(FunctionExpansion1D(trial, Vector<Real>{ REAL(1.0), REAL(2.0) }), FunctionSpaceInputError);

		FunctionExpansion1D expansion(trial, Vector<Real>{ REAL(1.0), REAL(2.0), REAL(3.0) });
		REQUIRE_THROWS_AS(expansion.derivative(REAL(0.0), -1), FunctionSpaceInputError);
		REQUIRE_THROWS_AS(expansion.derivative(REAL(0.0), 3), FunctionSpaceError);
	}

	TEST_CASE("Projection - L2 projection recovers Legendre coefficients", "[FunctionSpaces][Projection][OrthogonalBasis]")
	{
		LegendreBasis basis;
		OrthogonalBasisTrialSpace1D trial(basis, 3);
		LegendreCombinationFunction function;

		FunctionExpansion1D expansion = ProjectL2(function, trial, ProjectionTolerance);

		REQUIRE(expansion.coefficients()[0] == Catch::Approx(REAL(3.0)).epsilon(ProjectionTolerance));
		REQUIRE(expansion.coefficients()[1] == Catch::Approx(REAL(2.0)).epsilon(ProjectionTolerance));
		REQUIRE(expansion.coefficients()[2] == Catch::Approx(REAL(5.0)).epsilon(ProjectionTolerance));
		REQUIRE(expansion.evaluate(REAL(0.25)) == Catch::Approx(function(REAL(0.25))).epsilon(ProjectionTolerance));
	}

	TEST_CASE("Projection - OrthogonalBasis trial space validates dimension", "[FunctionSpaces][Projection][OrthogonalBasis]")
	{
		LegendreBasis basis;
		REQUIRE_THROWS_AS(OrthogonalBasisTrialSpace1D(basis, 0), FunctionSpaceInputError);
	}

	TEST_CASE("Interpolation - Nodal interpolation recovers monomial coefficients", "[FunctionSpaces][Interpolation]")
	{
		L2IntervalSpace space(-REAL(1.0), REAL(1.0));
		MonomialTrialSpace3 trial(space);
		LegendreCombinationFunction function;

		FunctionExpansion1D expansion = Interpolate(function, trial);

		REQUIRE(expansion.coefficients()[0] == Catch::Approx(REAL(0.5)).epsilon(REAL(1e-10)));
		REQUIRE(expansion.coefficients()[1] == Catch::Approx(REAL(2.0)).epsilon(REAL(1e-10)));
		REQUIRE(expansion.coefficients()[2] == Catch::Approx(REAL(7.5)).epsilon(REAL(1e-10)));
		REQUIRE(expansion.evaluate(trial.node(0)) == Catch::Approx(function(trial.node(0))).epsilon(REAL(1e-10)));
		REQUIRE(expansion.evaluate(trial.node(1)) == Catch::Approx(function(trial.node(1))).epsilon(REAL(1e-10)));
		REQUIRE(expansion.evaluate(trial.node(2)) == Catch::Approx(function(trial.node(2))).epsilon(REAL(1e-10)));
	}

	TEST_CASE("Interpolation - Rejects trial spaces without nodes", "[FunctionSpaces][Interpolation]")
	{
		LegendreBasis basis;
		OrthogonalBasisTrialSpace1D trial(basis, 3);
		LegendreCombinationFunction function;

		REQUIRE_THROWS_AS(Interpolate(function, trial), FunctionSpaceInputError);
	}

	TEST_CASE("FunctionSpaces diagnostics - TryProjectL2 reports success", "[FunctionSpaces][Diagnostics]")
	{
		LegendreBasis basis;
		OrthogonalBasisTrialSpace1D trial(basis, 3);
		LegendreCombinationFunction function;

		auto result = TryProjectL2(function, trial, ProjectionTolerance);

		REQUIRE(result.success());
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.algorithm_name == "ProjectL2");
		REQUIRE(result.expansion.has_value());
		REQUIRE(result.expansion->coefficients()[2] == Catch::Approx(REAL(5.0)).epsilon(ProjectionTolerance));
	}

	TEST_CASE("FunctionSpaces diagnostics - TryInterpolate reports unsupported spaces", "[FunctionSpaces][Diagnostics]")
	{
		LegendreBasis basis;
		OrthogonalBasisTrialSpace1D trial(basis, 3);
		LegendreCombinationFunction function;

		auto result = TryInterpolate(function, trial);

		REQUIRE_FALSE(result.success());
		REQUIRE_FALSE(result.expansion.has_value());
		REQUIRE(result.status == AlgorithmStatus::InvalidInput);
		REQUIRE(result.algorithm_name == "Interpolate");
		REQUIRE(result.error_message.find("nodes") != std::string::npos);
	}

	TEST_CASE("DenseBVPSolveResult1D - Reports success and failure states", "[FunctionSpaces][BVP]")
	{
		L2IntervalSpace space(-REAL(1.0), REAL(1.0));
		MonomialTrialSpace3 trial(space);
		FunctionExpansion1D expansion(trial, Vector<Real>{ REAL(1.0), REAL(0.0), REAL(0.0) });

		DenseBVPSolveResult1D success = DenseBVPSolveResult1D::Success(expansion, REAL(1e-12), 2);
		REQUIRE(success.success());
		REQUIRE(success.status == AlgorithmStatus::Success);
		REQUIRE(success.solution.has_value());
		REQUIRE(success.dimension == 3);
		REQUIRE(success.boundary_conditions_applied == 2);
		REQUIRE(success.residual_norm == REAL(1e-12));

		DenseBVPSolveResult1D failure = DenseBVPSolveResult1D::Failure(AlgorithmStatus::SingularMatrix, "singular collocation system");
		REQUIRE_FALSE(failure.success());
		REQUIRE_FALSE(failure.solution.has_value());
		REQUIRE(failure.status == AlgorithmStatus::SingularMatrix);
		REQUIRE(failure.error_message == "singular collocation system");
	}

	TEST_CASE("DenseBVPSolver1D - Solves Poisson Dirichlet benchmark", "[FunctionSpaces][BVP]")
	{
		const int dimension = std::is_same_v<Real, float> ? 12 : 24;
		ChebyshevCollocationSpace1D trial(-REAL(1.0), REAL(1.0), dimension);
		auto minusSecondDerivative = LinearDifferentialOperator1D::SecondOrder(
			[](Real) { return -REAL(1.0); },
			[](Real) { return REAL(0.0); },
			[](Real) { return REAL(0.0); });
		PoissonSineRhs rhs;
		BoundaryConditions1D boundaryConditions{
			BoundaryCondition1D::Dirichlet(-REAL(1.0), REAL(0.0)),
			BoundaryCondition1D::Dirichlet(REAL(1.0), REAL(0.0))
		};

		DenseBVPSolveResult1D result = SolveDenseCollocationBVP(minusSecondDerivative, rhs, trial, boundaryConditions);

		INFO(result.error_message);
		REQUIRE(result.success());
		REQUIRE(result.status == AlgorithmStatus::Success);
		REQUIRE(result.solution.has_value());
		REQUIRE(result.dimension == trial.dimension());
		REQUIRE(result.boundary_conditions_applied == 2);

		for (Real x : { -REAL(0.75), -REAL(0.25), REAL(0.0), REAL(0.25), REAL(0.75) })
			REQUIRE(result.solution->evaluate(x) == Catch::Approx(std::sin(Constants::PI * x)).margin(BvpTolerance));

		Real residual = CollocationResidualMaxNorm(*result.solution, minusSecondDerivative, rhs, trial);
		REQUIRE(residual < BvpResidualTolerance);

		Vector<Real> samplePoints = UniformSamplePoints(trial.functionSpace(), 5);
		Vector<Real> samples = SampleExpansion(*result.solution, samplePoints);
		REQUIRE(samplePoints[0] == -REAL(1.0));
		REQUIRE(samplePoints[2] == REAL(0.0));
		REQUIRE(samplePoints[4] == REAL(1.0));
		for (int index = 0; index < samplePoints.size(); ++index)
			REQUIRE(samples[index] == Catch::Approx(std::sin(Constants::PI * samplePoints[index])).margin(BvpTolerance));
	}

	TEST_CASE("BVPDiagnostics1D - Uniform sample points validates count", "[FunctionSpaces][BVP]")
	{
		L2IntervalSpace space(REAL(2.0), REAL(6.0));
		Vector<Real> single = UniformSamplePoints(space, 1);

		REQUIRE(single.size() == 1);
		REQUIRE(single[0] == REAL(4.0));
		REQUIRE_THROWS_AS(UniformSamplePoints(space, 0), FunctionSpaceInputError);
	}

	TEST_CASE("ChebyshevCollocationSpace1D - Generates Lobatto nodes left to right", "[FunctionSpaces][ChebyshevCollocation]")
	{
		ChebyshevCollocationSpace1D space(-REAL(1.0), REAL(1.0), 5);

		REQUIRE(space.dimension() == 5);
		REQUIRE(space.hasNodes());
		REQUIRE(space.nodeCount() == 5);
		REQUIRE(space.node(0) == Catch::Approx(-REAL(1.0)).margin(NodeTolerance));
		REQUIRE(space.node(1) == Catch::Approx(-std::sqrt(REAL(0.5))).margin(NodeTolerance));
		REQUIRE(space.node(2) == Catch::Approx(REAL(0.0)).margin(NodeTolerance));
		REQUIRE(space.node(3) == Catch::Approx(std::sqrt(REAL(0.5))).margin(NodeTolerance));
		REQUIRE(space.node(4) == Catch::Approx(REAL(1.0)).margin(NodeTolerance));
		REQUIRE(space.nodes().size() == 5);
	}

	TEST_CASE("ChebyshevCollocationSpace1D - Scales Lobatto nodes to custom interval", "[FunctionSpaces][ChebyshevCollocation]")
	{
		ChebyshevCollocationSpace1D space(REAL(2.0), REAL(6.0), 3);

		REQUIRE(space.node(0) == Catch::Approx(REAL(2.0)).margin(REAL(1e-14)));
		REQUIRE(space.node(1) == Catch::Approx(REAL(4.0)).margin(REAL(1e-14)));
		REQUIRE(space.node(2) == Catch::Approx(REAL(6.0)).margin(REAL(1e-14)));
		REQUIRE(space.basisValue(0, REAL(4.0)) == REAL(1.0));
		REQUIRE(space.basisValue(1, REAL(4.0)) == Catch::Approx(REAL(0.0)).margin(REAL(1e-14)));
		REQUIRE(space.basisValue(2, REAL(4.0)) == Catch::Approx(-REAL(1.0)).margin(REAL(1e-14)));
	}

	TEST_CASE("ChebyshevCollocationSpace1D - Rejects invalid dimensions", "[FunctionSpaces][ChebyshevCollocation]")
	{
		REQUIRE_THROWS_AS(ChebyshevCollocationSpace1D(-REAL(1.0), REAL(1.0), 1), FunctionSpaceInputError);
	}

	TEST_CASE("ChebyshevCollocationSpace1D - First derivative matrix differentiates polynomials", "[FunctionSpaces][ChebyshevCollocation]")
	{
		ChebyshevCollocationSpace1D space(-REAL(1.0), REAL(1.0), 6);
		Matrix<Real> derivative = space.firstDerivativeMatrix();

		Vector<Real> constant(space.dimension());
		Vector<Real> quadratic(space.dimension());
		for (int i = 0; i < space.dimension(); ++i) {
			Real x = space.node(i);
			constant[i] = REAL(1.0);
			quadratic[i] = x * x;
		}

		for (int row = 0; row < space.dimension(); ++row) {
			Real constantDerivative = REAL(0.0);
			Real quadraticDerivative = REAL(0.0);
			for (int column = 0; column < space.dimension(); ++column) {
				constantDerivative += derivative(row, column) * constant[column];
				quadraticDerivative += derivative(row, column) * quadratic[column];
			}

			REQUIRE(constantDerivative == Catch::Approx(REAL(0.0)).margin(FirstDerivativeTolerance));
			REQUIRE(quadraticDerivative == Catch::Approx(REAL(2.0) * space.node(row)).margin(FirstDerivativeTolerance));
		}
	}

	TEST_CASE("ChebyshevCollocationSpace1D - First derivative matrix scales to interval", "[FunctionSpaces][ChebyshevCollocation]")
	{
		ChebyshevCollocationSpace1D space(REAL(2.0), REAL(6.0), 6);
		Matrix<Real> derivative = space.firstDerivativeMatrix();
		Vector<Real> linear(space.dimension());

		for (int i = 0; i < space.dimension(); ++i)
			linear[i] = space.node(i);

		for (int row = 0; row < space.dimension(); ++row) {
			Real value = REAL(0.0);
			for (int column = 0; column < space.dimension(); ++column)
				value += derivative(row, column) * linear[column];

			REQUIRE(value == Catch::Approx(REAL(1.0)).margin(FirstDerivativeTolerance));
		}
	}

	TEST_CASE("ChebyshevCollocationSpace1D - Second derivative matrix differentiates polynomials", "[FunctionSpaces][ChebyshevCollocation]")
	{
		ChebyshevCollocationSpace1D space(-REAL(1.0), REAL(1.0), 7);
		Matrix<Real> second = space.secondDerivativeMatrix();
		Vector<Real> cubic(space.dimension());

		for (int i = 0; i < space.dimension(); ++i) {
			Real x = space.node(i);
			cubic[i] = x * x * x;
		}

		for (int row = 0; row < space.dimension(); ++row) {
			Real value = REAL(0.0);
			for (int column = 0; column < space.dimension(); ++column)
				value += second(row, column) * cubic[column];

			REQUIRE(value == Catch::Approx(REAL(6.0) * space.node(row)).margin(FirstDerivativeTolerance));
		}
	}

	TEST_CASE("ChebyshevCollocationSpace1D - Second derivative matrix scales to interval", "[FunctionSpaces][ChebyshevCollocation]")
	{
		ChebyshevCollocationSpace1D space(REAL(2.0), REAL(6.0), 7);
		Matrix<Real> second = space.secondDerivativeMatrix();
		Vector<Real> quadratic(space.dimension());

		for (int i = 0; i < space.dimension(); ++i) {
			Real x = space.node(i);
			quadratic[i] = x * x;
		}

		for (int row = 0; row < space.dimension(); ++row) {
			Real value = REAL(0.0);
			for (int column = 0; column < space.dimension(); ++column)
				value += second(row, column) * quadratic[column];

			REQUIRE(value == Catch::Approx(REAL(2.0)).margin(SecondDerivativeTolerance));
		}
	}

	TEST_CASE("ChebyshevCollocationSpace1D - Smooth derivative error decreases with N", "[FunctionSpaces][ChebyshevCollocation]")
	{
		auto maxDerivativeError = [](int dimension) {
			ChebyshevCollocationSpace1D space(-REAL(1.0), REAL(1.0), dimension);
			Matrix<Real> derivative = space.firstDerivativeMatrix();
			Vector<Real> values(space.dimension());
			for (int i = 0; i < space.dimension(); ++i)
				values[i] = std::sin(Constants::PI * space.node(i));

			Real maxError = REAL(0.0);
			for (int row = 0; row < space.dimension(); ++row) {
				Real numerical = REAL(0.0);
				for (int column = 0; column < space.dimension(); ++column)
					numerical += derivative(row, column) * values[column];
				Real exact = Constants::PI * std::cos(Constants::PI * space.node(row));
				maxError = std::max(maxError, std::abs(numerical - exact));
			}
			return maxError;
		};

		Real error10 = maxDerivativeError(10);
		Real error20 = maxDerivativeError(20);
		Real error40 = maxDerivativeError(40);

		REQUIRE(error20 < error10);
		REQUIRE(error40 < error10);
		REQUIRE(error20 < SmoothDerivativeTolerance);
		REQUIRE(error40 < SmoothDerivativeTolerance);
	}

	TEST_CASE("LinearDifferentialOperator1D - Applies coefficients to trial basis", "[FunctionSpaces][DifferentialOperator1D]")
	{
		L2IntervalSpace space(-REAL(1.0), REAL(1.0));
		MonomialTrialSpace3 trial(space);
		auto op = LinearDifferentialOperator1D::SecondOrder(
			[](Real x) { return -REAL(1.0) + REAL(0.0) * x; },
			[](Real x) { return REAL(2.0) * x; },
			[](Real) { return REAL(3.0); });

		REQUIRE(op.order() == 2);
		REQUIRE(op.coefficient(0, REAL(0.25)) == REAL(3.0));
		REQUIRE(op.coefficient(1, REAL(0.25)) == REAL(0.5));
		REQUIRE(op.coefficient(2, REAL(0.25)) == -REAL(1.0));

		// phi_2=x^2: 3*x^2 + 2*x*(2*x) - 1*2 = 7*x^2 - 2
		REQUIRE(op.applyToBasis(trial, 2, REAL(0.5)) == Catch::Approx(-REAL(0.25)).margin(REAL(1e-14)));
	}

	TEST_CASE("LinearDifferentialOperator1D - Rejects invalid coefficients and missing derivatives", "[FunctionSpaces][DifferentialOperator1D]")
	{
		L2IntervalSpace space(-REAL(1.0), REAL(1.0));
		LegendreBasis basis;
		OrthogonalBasisTrialSpace1D trial(basis, 3);

		REQUIRE_THROWS_AS(LinearDifferentialOperator1D({}), FunctionSpaceInputError);
		REQUIRE_THROWS_AS(LinearDifferentialOperator1D({ nullptr }), FunctionSpaceInputError);

		auto secondOrder = LinearDifferentialOperator1D::SecondOrder(
			[](Real) { return REAL(1.0); },
			[](Real) { return REAL(0.0); },
			[](Real) { return REAL(0.0); });
		REQUIRE_THROWS_AS(secondOrder.coefficient(3, REAL(0.0)), FunctionSpaceInputError);
		REQUIRE_THROWS_AS(secondOrder.applyToBasis(trial, 0, REAL(0.0)), FunctionSpaceError);
	}

	TEST_CASE("BoundaryCondition1D - Factory methods define standard conditions", "[FunctionSpaces][BoundaryCondition1D]")
	{
		auto dirichlet = BoundaryCondition1D::Dirichlet(-REAL(1.0), REAL(2.0));
		auto neumann = BoundaryCondition1D::Neumann(REAL(1.0), -REAL(3.0));
		auto robin = BoundaryCondition1D::Robin(REAL(0.0), REAL(4.0), REAL(5.0), REAL(6.0));
		auto periodic = BoundaryCondition1D::Periodic();

		REQUIRE(dirichlet.kind == BoundaryConditionKind::Dirichlet);
		REQUIRE(dirichlet.alpha == REAL(1.0));
		REQUIRE(dirichlet.beta == REAL(0.0));
		REQUIRE(dirichlet.value == REAL(2.0));
		REQUIRE(neumann.kind == BoundaryConditionKind::Neumann);
		REQUIRE(neumann.alpha == REAL(0.0));
		REQUIRE(neumann.beta == REAL(1.0));
		REQUIRE(robin.kind == BoundaryConditionKind::Robin);
		REQUIRE(robin.alpha == REAL(4.0));
		REQUIRE(robin.beta == REAL(5.0));
		REQUIRE(periodic.kind == BoundaryConditionKind::Periodic);
	}

	TEST_CASE("BoundaryConditions1D - Validates condition collections", "[FunctionSpaces][BoundaryCondition1D]")
	{
		L2IntervalSpace space(-REAL(1.0), REAL(1.0));
		BoundaryConditions1D conditions{
			BoundaryCondition1D::Dirichlet(-REAL(1.0), REAL(0.0)),
			BoundaryCondition1D::Neumann(REAL(1.0), REAL(0.0))
		};

		REQUIRE(conditions.size() == 2);
		REQUIRE_FALSE(conditions.empty());
		REQUIRE(conditions[0].kind == BoundaryConditionKind::Dirichlet);
		REQUIRE_NOTHROW(conditions.validate(space));

		conditions.Add(BoundaryCondition1D::Robin(REAL(0.0), REAL(1.0), REAL(1.0), REAL(0.0)));
		REQUIRE(conditions.size() == 3);
		REQUIRE_NOTHROW(conditions.validate(space));
	}

	TEST_CASE("BoundaryCondition1D - Rejects invalid descriptors", "[FunctionSpaces][BoundaryCondition1D]")
	{
		L2IntervalSpace space(-REAL(1.0), REAL(1.0));
		REQUIRE_THROWS_AS(BoundaryCondition1D::Dirichlet(REAL(2.0), REAL(0.0)).validate(space), FunctionSpaceInputError);
		REQUIRE_THROWS_AS(BoundaryCondition1D::Robin(REAL(0.0), REAL(0.0), REAL(0.0), REAL(0.0)).validate(space), FunctionSpaceInputError);

		BoundaryCondition1D invalid = BoundaryCondition1D::Dirichlet(REAL(0.0), std::numeric_limits<Real>::quiet_NaN());
		REQUIRE_THROWS_AS(invalid.validate(space), FunctionSpaceInputError);
	}

	TEST_CASE("BoundaryCondition1D - Builds Dirichlet Neumann and Robin rows", "[FunctionSpaces][BoundaryCondition1D]")
	{
		L2IntervalSpace space(-REAL(1.0), REAL(1.0));
		MonomialTrialSpace3 trial(space);

		Vector<Real> dirichlet = BoundaryCondition1D::Dirichlet(REAL(0.5), REAL(0.0)).row(trial);
		REQUIRE(dirichlet[0] == REAL(1.0));
		REQUIRE(dirichlet[1] == REAL(0.5));
		REQUIRE(dirichlet[2] == REAL(0.25));

		Vector<Real> neumann = BoundaryCondition1D::Neumann(REAL(0.5), REAL(0.0)).row(trial);
		REQUIRE(neumann[0] == REAL(0.0));
		REQUIRE(neumann[1] == REAL(1.0));
		REQUIRE(neumann[2] == REAL(1.0));

		Vector<Real> robin = BoundaryCondition1D::Robin(REAL(0.5), REAL(2.0), REAL(3.0), REAL(0.0)).row(trial);
		REQUIRE(robin[0] == REAL(2.0));
		REQUIRE(robin[1] == REAL(4.0));
		REQUIRE(robin[2] == Catch::Approx(REAL(3.5)).margin(REAL(1e-14)));
	}

	TEST_CASE("BoundaryCondition1D - Row helper rejects periodic paired assembly", "[FunctionSpaces][BoundaryCondition1D]")
	{
		L2IntervalSpace space(-REAL(1.0), REAL(1.0));
		MonomialTrialSpace3 trial(space);

		REQUIRE_THROWS_AS(BoundaryCondition1D::Periodic().row(trial), FunctionSpaceError);
	}

	TEST_CASE("OperatorAssembly - Builds dense collocation matrix", "[FunctionSpaces][OperatorAssembly]")
	{
		L2IntervalSpace space(-REAL(1.0), REAL(1.0));
		MonomialTrialSpace3 trial(space);
		auto op = LinearDifferentialOperator1D::SecondOrder(
			[](Real) { return -REAL(1.0); },
			[](Real x) { return REAL(2.0) * x; },
			[](Real) { return REAL(3.0); });

		Matrix<Real> matrix = BuildCollocationMatrix(op, trial);

		REQUIRE(matrix.rows() == trial.dimension());
		REQUIRE(matrix.cols() == trial.dimension());
		for (int row = 0; row < trial.dimension(); ++row) {
			Real x = trial.node(row);
			for (int column = 0; column < trial.dimension(); ++column)
				REQUIRE(matrix(row, column) == Catch::Approx(op.applyToBasis(trial, column, x)).margin(REAL(1e-14)));
		}
	}

	TEST_CASE("OperatorAssembly - Rejects trial spaces without nodes", "[FunctionSpaces][OperatorAssembly]")
	{
		LegendreBasis basis;
		OrthogonalBasisTrialSpace1D trial(basis, 3);
		auto op = LinearDifferentialOperator1D::Multiplication([](Real) { return REAL(1.0); });

		REQUIRE_THROWS_AS(BuildCollocationMatrix(op, trial), FunctionSpaceInputError);
	}

	TEST_CASE("OperatorAssembly - Applies boundary row replacement", "[FunctionSpaces][OperatorAssembly][BoundaryCondition1D]")
	{
		L2IntervalSpace space(-REAL(1.0), REAL(1.0));
		MonomialTrialSpace3 trial(space);
		Matrix<Real> matrix(3, 3, { REAL(9.0), REAL(9.0), REAL(9.0), REAL(8.0), REAL(8.0), REAL(8.0), REAL(7.0), REAL(7.0), REAL(7.0) });
		Vector<Real> rhs{ REAL(1.0), REAL(2.0), REAL(3.0) };

		ApplyBoundaryRow(matrix, rhs, 0, trial, BoundaryCondition1D::Dirichlet(-REAL(1.0), REAL(5.0)));
		ApplyBoundaryRow(matrix, rhs, 2, trial, BoundaryCondition1D::Neumann(REAL(1.0), -REAL(2.0)));

		REQUIRE(matrix(0, 0) == REAL(1.0));
		REQUIRE(matrix(0, 1) == -REAL(1.0));
		REQUIRE(matrix(0, 2) == REAL(1.0));
		REQUIRE(rhs[0] == REAL(5.0));

		REQUIRE(matrix(2, 0) == REAL(0.0));
		REQUIRE(matrix(2, 1) == REAL(1.0));
		REQUIRE(matrix(2, 2) == REAL(2.0));
		REQUIRE(rhs[2] == -REAL(2.0));

		REQUIRE(matrix(1, 0) == REAL(8.0));
		REQUIRE(rhs[1] == REAL(2.0));
	}

	TEST_CASE("OperatorAssembly - Boundary row replacement validates dimensions", "[FunctionSpaces][OperatorAssembly][BoundaryCondition1D]")
	{
		L2IntervalSpace space(-REAL(1.0), REAL(1.0));
		MonomialTrialSpace3 trial(space);
		Matrix<Real> badMatrix(2, 2);
		Matrix<Real> matrix(3, 3);
		Vector<Real> badRhs(2);
		Vector<Real> rhs(3);
		auto condition = BoundaryCondition1D::Dirichlet(REAL(0.0), REAL(0.0));

		REQUIRE_THROWS_AS(ApplyBoundaryRow(badMatrix, rhs, 0, trial, condition), MatrixDimensionError);
		REQUIRE_THROWS_AS(ApplyBoundaryRow(matrix, badRhs, 0, trial, condition), VectorDimensionError);
		REQUIRE_THROWS_AS(ApplyBoundaryRow(matrix, rhs, 3, trial, condition), FunctionSpaceInputError);
	}

	TEST_CASE("OperatorAssembly - Builds simple Galerkin stiffness matrix", "[FunctionSpaces][OperatorAssembly][Galerkin]")
	{
		L2IntervalSpace space(-REAL(1.0), REAL(1.0));
		MonomialTrialSpace3 trial(space);
		ConstantFunction one(REAL(1.0));
		ConstantFunction zero(REAL(0.0));

		Matrix<Real> stiffness = BuildGalerkinStiffnessMatrix(trial, one, zero, REAL(1e-9));

		REQUIRE(stiffness(0, 0) == Catch::Approx(REAL(0.0)).margin(REAL(1e-10)));
		REQUIRE(stiffness(0, 1) == Catch::Approx(REAL(0.0)).margin(REAL(1e-10)));
		REQUIRE(stiffness(0, 2) == Catch::Approx(REAL(0.0)).margin(REAL(1e-10)));
		REQUIRE(stiffness(1, 1) == Catch::Approx(REAL(2.0)).epsilon(QuadratureTolerance));
		REQUIRE(stiffness(1, 2) == Catch::Approx(REAL(0.0)).margin(REAL(1e-10)));
		REQUIRE(stiffness(2, 2) == Catch::Approx(REAL(8.0) / REAL(3.0)).epsilon(QuadratureTolerance));
		REQUIRE(stiffness(2, 1) == Catch::Approx(stiffness(1, 2)).margin(REAL(1e-12)));
	}

	TEST_CASE("LinearOperator - Matrix-free apply validates dimensions", "[FunctionSpaces][LinearOperator]")
	{
		IdentityLinearOperator3 op;
		Vector<Real> input{ REAL(1.0), REAL(2.0), REAL(3.0) };

		Vector<Real> output = op.apply(input);
		REQUIRE(output.size() == 3);
		REQUIRE(output[0] == REAL(1.0));
		REQUIRE(output[1] == REAL(2.0));
		REQUIRE(output[2] == REAL(3.0));

		Vector<Real> badInput(2);
		Vector<Real> badOutput(2);
		REQUIRE_THROWS_AS(op.apply(badInput), VectorDimensionError);
		REQUIRE_THROWS_AS(op.apply(input, badOutput), VectorDimensionError);
	}

	TEST_CASE("MatrixFreeOperators1D - Dirichlet second derivative matches dense stencil", "[FunctionSpaces][LinearOperator][MatrixFree]")
	{
		DirichletSecondDerivativeOperator1D op(4, REAL(0.25));
		Vector<Real> input{ REAL(1.0), REAL(2.0), REAL(4.0), REAL(8.0) };
		Vector<Real> output = op.apply(input);

		Matrix<Real> dense(4, 4, {
			-REAL(2.0), REAL(1.0), REAL(0.0), REAL(0.0),
			REAL(1.0), -REAL(2.0), REAL(1.0), REAL(0.0),
			REAL(0.0), REAL(1.0), -REAL(2.0), REAL(1.0),
			REAL(0.0), REAL(0.0), REAL(1.0), -REAL(2.0)
		});
		Real h2inv = REAL(16.0);

		for (int row = 0; row < 4; ++row) {
			Real expected = REAL(0.0);
			for (int column = 0; column < 4; ++column)
				expected += h2inv * dense(row, column) * input[column];
			REQUIRE(output[row] == Catch::Approx(expected).margin(REAL(1e-14)));
		}
	}

	TEST_CASE("MatrixFreeOperators1D - Dirichlet second derivative includes boundary values", "[FunctionSpaces][LinearOperator][MatrixFree]")
	{
		DirichletSecondDerivativeOperator1D op(2, REAL(0.5), REAL(1.0), REAL(3.0));
		Vector<Real> input{ REAL(2.0), REAL(4.0) };
		Vector<Real> output = op.apply(input);

		REQUIRE(output[0] == Catch::Approx(REAL(4.0) * (REAL(1.0) - REAL(4.0) + REAL(4.0))).margin(REAL(1e-14)));
		REQUIRE(output[1] == Catch::Approx(REAL(4.0) * (REAL(2.0) - REAL(8.0) + REAL(3.0))).margin(REAL(1e-14)));
		REQUIRE(op.gridSpacing() == REAL(0.5));
		REQUIRE(op.leftBoundaryValue() == REAL(1.0));
		REQUIRE(op.rightBoundaryValue() == REAL(3.0));
	}

	TEST_CASE("MatrixFreeOperators1D - Rejects invalid construction", "[FunctionSpaces][LinearOperator][MatrixFree]")
	{
		REQUIRE_THROWS_AS(DirichletSecondDerivativeOperator1D(0, REAL(1.0)), FunctionSpaceInputError);
		REQUIRE_THROWS_AS(DirichletSecondDerivativeOperator1D(4, REAL(0.0)), FunctionSpaceInputError);
	}
}
