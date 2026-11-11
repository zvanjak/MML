#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/core/OrthogonalBasis/ChebyshevBasis.h>
#include <mml/core/OrthogonalBasis/HermiteBasis.h>
#include <mml/core/OrthogonalBasis/LaguerreBasis.h>
#include <mml/core/OrthogonalBasis/LegendreBasis.h>
#endif
#include <array>
#include <functional>

using namespace MML;

namespace MML::Tests::Core::OrthogonalBasisContracts
{
	Real CoefficientBudget(Real inputNormSquared, Real basisNormSquared)
	{
		// Cauchy-Schwarz bounds integral(|f * basis|) / ||basis||^2.
		const Real condition = std::max(Real(1), std::sqrt(inputNormSquared / basisNormSquared));
		return std::max(Precision::DefaultToleranceStrict, 2 * Precision::IntegrationRoundoffFloorScale * condition);
	}

	class Function : public IRealFunction
	{
		std::function<Real(Real)> _function;
	public:
		mutable int calls = 0;
		explicit Function(std::function<Real(Real)> function) : _function(std::move(function)) {}
		Real operator()(Real x) const override
		{
			++calls;
			return _function(x);
		}
	};

	struct TrapBudget
	{
		int saved = Defaults::TrapezoidIntegrationMaxSteps;
		TrapBudget() { Defaults::TrapezoidIntegrationMaxSteps = 8; }
		~TrapBudget() { Defaults::TrapezoidIntegrationMaxSteps = saved; }
	};

	struct Bases
	{
		LegendreBasis legendre;
		ChebyshevBasis chebyshev;
		HermiteBasis hermite;
		LaguerreBasis laguerre;
		std::array<const OrthogonalBasis*, 4> all{&legendre, &chebyshev, &hermite, &laguerre};
	};

	class TestBasis : public OrthogonalBasis
	{
		Real _norm, _a, _b;
	public:
		TestBasis(Real norm, Real a, Real b) : _norm(norm), _a(a), _b(b) {}
		Real Evaluate(int, Real) const override { return 1; }
		Real WeightFunction(Real x) const override { return std::exp(-std::abs(x)); }
		Real Normalization(int) const override { return _norm; }
		Real DomainMin() const override { return _a; }
		Real DomainMax() const override { return _b; }
	};

	TEST_CASE("Weighted basis coefficients integrate the full domain", "[OrthogonalBasis][basis-contract][coefficients]")
	{
		TrapBudget budget;
		Bases bases;
		for (const auto* basis : bases.all)
		{
			Function constant([](Real x) {
				if (!std::isfinite(x)) throw DomainError("coefficient sampled a nonfinite coordinate");
				return Real{1};
			});
			CHECK(basis->ComputeCoefficient(constant, 0) == Catch::Approx(1).margin(Precision::DefaultToleranceStrict));
			CHECK(constant.calls < 10000);
		}
	}

	TEST_CASE("Weighted basis projections recover their Gram matrix", "[OrthogonalBasis][basis-contract][coefficients]")
	{
		TrapBudget budget;
		Bases bases;
		for (const auto* basis : bases.all)
			for (int m = 0; m <= 4; ++m)
				for (int n = 0; n <= 4; ++n)
				{
					CAPTURE(m, n, basis->DomainMin(), basis->DomainMax());
					Function function([&](Real x) { return basis->Evaluate(m, x); });
					const Real eps = CoefficientBudget(basis->Normalization(m), basis->Normalization(n));
					CHECK(basis->ComputeCoefficient(function, n, eps) ==
						Catch::Approx(m == n ? 1 : 0).margin(eps));
					CHECK(function.calls < 10000);
				}
	}

	TEST_CASE("Weighted basis coefficients reconstruct a quadratic", "[OrthogonalBasis][basis-contract][coefficients]")
	{
		TrapBudget budget;
		Bases bases;
		const std::array<std::array<Real, 5>, 4> expected{{
			{Real(1) / 3, 0, Real(2) / 3, 0, 0},
			{Real(0.5), 0, Real(0.5), 0, 0},
			{Real(0.5), 0, Real(0.25), 0, 0},
			{2, -4, 2, 0, 0}
		}};
		Function quadratic([](Real x) { return x * x; });
		for (std::size_t i = 0; i < bases.all.size(); ++i)
		{
			Real inputNormSquared = 0;
			for (int n = 0; n <= 4; ++n)
				inputNormSquared += expected[i][n] * expected[i][n] * bases.all[i]->Normalization(n);
			for (int n = 0; n <= 4; ++n)
			{
				const Real eps = CoefficientBudget(inputNormSquared, bases.all[i]->Normalization(n));
				CHECK(bases.all[i]->ComputeCoefficient(quadratic, n, eps) ==
					Catch::Approx(expected[i][n]).margin(eps));
			}
		}
	}

	TEST_CASE("Laguerre coefficients retain the exponential tail", "[OrthogonalBasis][basis-contract][coefficients]")
	{
		TrapBudget budget;
		const Real beta = GENERATE(Real(0.1), Real{1}, Real{10});
		LaguerreBasis basis;
		Function function([&](Real x) { return std::exp(-beta * x); });
		for (int n = 0; n <= 6; ++n)
		{
			const Real expected = std::pow(beta / (1 + beta), MML::Real(n)) / (1 + beta);
			CHECK(basis.ComputeCoefficient(function, n) ==
				Catch::Approx(expected).epsilon(Precision::DefaultToleranceStrict).margin(Precision::DefaultToleranceStrict));
		}
	}

	TEST_CASE("Hermite coefficients agree with the Gaussian Fourier transform", "[OrthogonalBasis][basis-contract][coefficients]")
	{
		TrapBudget budget;
		const Real frequency = GENERATE(Real(0.2), Real{1}, Real{3});
		HermiteBasis basis;
		Function function([&](Real x) { return std::cos(frequency * x); });
		Real expected = std::exp(-frequency * frequency / 4);
		for (int n = 0; n <= 6; ++n)
		{
			if (n != 0 && n % 2 == 0)
				expected *= -frequency * frequency / (4 * n * (n - 1));
			CHECK(basis.ComputeCoefficient(function, n) ==
				Catch::Approx(n % 2 == 0 ? expected : 0).margin(Precision::DefaultToleranceStrict));
		}
	}

	TEST_CASE("Associated Legendre uses raw x polynomials without phase", "[OrthogonalBasis][basis-contract][conventions]")
	{
		AssociatedLegendreBasis m0(0), m1(1), m2(2), m3(3);
		CHECK(m0.Evaluate(0, 0) == Catch::Approx(1));
		CHECK(m0.Evaluate(1, 0) == Catch::Approx(0).margin(Precision::RoundoffZeroThreshold));
		CHECK(m0.Evaluate(2, 0) == Catch::Approx(-0.5));
		CHECK(m1.Evaluate(1, 0) == Catch::Approx(1));
		CHECK(m1.Evaluate(2, Real(0.5)) == Catch::Approx(Real(1.5) * std::sqrt(Real(0.75))));
		CHECK(m2.Evaluate(2, 0) == Catch::Approx(3));
		CHECK(m3.Evaluate(3, 0) == Catch::Approx(15));
		CHECK(m1.Evaluate(3, -1) == Real{0});
		CHECK(m1.Evaluate(3, 1) == Real{0});
		CHECK_THROWS_AS(m0.Evaluate(0, Real(1.1)), DomainError);
		CHECK_THROWS_AS(m0.Evaluate(0, std::numeric_limits<Real>::quiet_NaN()), DomainError);
	}

	TEST_CASE("Constant Chebyshev extrema are rejected explicitly", "[OrthogonalBasis][basis-contract][domains]")
	{
		ChebyshevBasis basis;
		CHECK_THROWS_AS(basis.GetExtrema(0), ArgumentError);
		const auto points = basis.GetExtrema(4);
		REQUIRE(points.size() == 5);
		CHECK(points.front() == Real{1});
		CHECK(points.back() == Real{-1});
		CHECK(std::abs(points[2]) < Precision::RoundoffZeroThreshold);
		CHECK(std::isinf(basis.WeightFunction(-1)));
		CHECK(std::isinf(basis.WeightFunction(1)));
		CHECK_THROWS_AS(basis.WeightFunction(Real(1.1)), DomainError);
		CHECK_THROWS_AS(basis.WeightFunction(std::numeric_limits<Real>::quiet_NaN()), DomainError);
	}

	TEST_CASE("Spherical associated Legendre values retain the explicit angular convention", "[OrthogonalBasis][basis-contract][conventions]")
	{
		AssociatedLegendreBasis m0(0), m1(1);
		CHECK(m0.EvaluateSpherical(0, Real(0.7)) == Catch::Approx(1 / std::sqrt(4 * Constants::PI)));
		CHECK(m0.EvaluateSpherical(1, 0) == Catch::Approx(std::sqrt(3 / (4 * Constants::PI))));
		CHECK(std::abs(m0.EvaluateSpherical(1, Constants::PI / 2)) < Precision::RoundoffZeroThreshold);
		CHECK(m1.EvaluateSpherical(1, Constants::PI / 2) == Catch::Approx(-std::sqrt(3 / (8 * Constants::PI))));
		CHECK_THROWS_AS(m1.EvaluateSpherical(0, 0), ArgumentError);
		CHECK_THROWS_AS(m0.EvaluateSpherical(0, -1), DomainError);
		CHECK_THROWS_AS(m0.EvaluateSpherical(0, Constants::PI + 1), DomainError);
		CHECK_THROWS_AS(m0.EvaluateSpherical(0, Constants::PosInf), DomainError);
	}

	TEST_CASE("Generic weighted coefficients handle each infinite domain and bad normalizations", "[OrthogonalBasis][basis-contract][coefficients]")
	{
		Function constant([](Real) { return Real{1}; });
		CHECK(TestBasis(1 - std::exp(Real(-1)), 0, 1).ComputeCoefficient(constant, 0) ==
			Catch::Approx(1).margin(Precision::DefaultToleranceStrict));
		CHECK(TestBasis(std::exp(Real(-2)), 2, Constants::PosInf).ComputeCoefficient(constant, 0) ==
			Catch::Approx(1).margin(Precision::DefaultToleranceStrict));
		CHECK(TestBasis(std::exp(Real(-2)), Constants::NegInf, -2).ComputeCoefficient(constant, 0) ==
			Catch::Approx(1).margin(Precision::DefaultToleranceStrict));
		CHECK(TestBasis(2, Constants::NegInf, Constants::PosInf).ComputeCoefficient(constant, 0) ==
			Catch::Approx(1).margin(Precision::DefaultToleranceStrict));
		Function huge([](Real) { return std::numeric_limits<Real>::max(); });
		CHECK(TestBasis(std::numeric_limits<Real>::max(), 0, Constants::PosInf).ComputeCoefficient(huge, 0) ==
			Catch::Approx(1).margin(Precision::DefaultToleranceStrict));
		Function tiny([](Real) { return std::numeric_limits<Real>::denorm_min(); });
		CHECK(TestBasis(std::numeric_limits<Real>::denorm_min(), 0, Constants::PosInf).ComputeCoefficient(tiny, 0) ==
			Catch::Approx(1).margin(Precision::DefaultToleranceStrict));
		for (Real norm : {Real{0}, Real{-1}, Constants::PosInf, std::numeric_limits<Real>::quiet_NaN()})
			CHECK_THROWS_AS(TestBasis(norm, 0, 1).ComputeCoefficient(constant, 0), NumericalMethodError);
		CHECK_THROWS_AS(TestBasis(1, 2, 1).ComputeCoefficient(constant, 0), ArgumentError);
		CHECK_THROWS_AS(TestBasis(1, Constants::PosInf, Constants::PosInf).ComputeCoefficient(constant, 0), ArgumentError);
	}

	TEST_CASE("Unachievable coefficient accuracy fails rather than returning a success", "[OrthogonalBasis][basis-contract][coefficients]")
	{
		LegendreBasis basis;
		Function constant([](Real) { return Real{1}; });
		CHECK_THROWS_AS(basis.ComputeCoefficient(constant, 0, Precision::Eps * Precision::Eps), ConvergenceError);
		CHECK(constant.calls < 10000);
	}

	TEST_CASE("Associated Laguerre normalization avoids factorial overflow", "[OrthogonalBasis][basis-contract][normalization]")
	{
		const int n = GENERATE(0, 20, 171, 256);
		CHECK(AssociatedLaguerreBasis(0).Normalization(n) == Real{1});
		CHECK(AssociatedLaguerreBasis(1).Normalization(n) == Catch::Approx(n + 1));
		CHECK(AssociatedLaguerreBasis(2).Normalization(n) == Catch::Approx(Real(n + 1) * (n + 2)));
		CHECK(AssociatedLaguerreBasis(3).Normalization(n) == Catch::Approx(Real(n + 1) * (n + 2) * (n + 3)));
		CHECK(AssociatedLaguerreBasis(Real(0.5)).Normalization(171) ==
			Catch::Approx(Real(13.1053493638911462348)).epsilon(Precision::DefaultToleranceStrict));
		CHECK(AssociatedLaguerreBasis(Real(-0.5)).Normalization(171) ==
			Catch::Approx(Real(0.076416031276333214197)).epsilon(Precision::DefaultToleranceStrict));
		if constexpr (std::numeric_limits<Real>::max_exponent10 > 234)
			CHECK(AssociatedLaguerreBasis(100).Normalization(171) ==
				Catch::Approx(1.45479986608789967123e234).epsilon(Precision::DefaultToleranceStrict));
		else
			CHECK_THROWS_AS(AssociatedLaguerreBasis(100).Normalization(171), NumericalMethodError);
	}

	TEST_CASE("Normalized Hermite waves survive intermediate overflow and underflow", "[OrthogonalBasis][basis-contract][normalization]")
	{
		HermiteBasis basis;
		CHECK(basis.QuantumWavefunction(152, 0) ==
			Catch::Approx(Real(0.190925755629646892434)).epsilon(Precision::DefaultToleranceStrict));
		CHECK(basis.QuantumWavefunction(1000, 40) ==
			Catch::Approx(Real(0.172250520732792269826)).epsilon(Testing::Tol(Precision::DefaultToleranceStrict, Real(0.003))));
		for (int n : {0, 1, 2, 50, 151, 152, 500, 1000})
		{
			const Real positive = basis.QuantumWavefunction(n, 3);
			CHECK(std::isfinite(positive));
			CHECK(basis.QuantumWavefunction(n, -3) == (n % 2 == 0 ? positive : -positive));
		}
		CHECK(basis.QuantumWavefunction(0, std::numeric_limits<Real>::max()) == Real{0});
		CHECK_THROWS_AS(basis.QuantumWavefunction(0, Constants::PosInf), DomainError);
		CHECK_THROWS_AS(basis.QuantumWavefunction(0, std::numeric_limits<Real>::quiet_NaN()), DomainError);
	}

	TEST_CASE("Weighted coefficient domains and callback failures remain explicit", "[OrthogonalBasis][basis-contract][domains]")
	{
		TrapBudget budget;
		Bases bases;
		Function constant([](Real) { return Real{1}; });
		Function bad([](Real) { return std::numeric_limits<Real>::quiet_NaN(); });
		struct CallbackFailure : std::runtime_error { CallbackFailure() : std::runtime_error("intentional callback failure") {} };
		Function throwing([](Real) -> Real { throw CallbackFailure(); });
		for (const auto* basis : bases.all)
		{
			for (Real eps : {Real{0}, Real{-1}, Constants::PosInf, std::numeric_limits<Real>::quiet_NaN()})
				CHECK_THROWS_AS(basis->ComputeCoefficient(constant, 0, eps), ArgumentError);
			CHECK_THROWS_AS(basis->ComputeCoefficient(constant, -1), ArgumentError);
			CHECK_THROWS_AS(basis->ComputeCoefficient(bad, 0), NumericalMethodError);
			CHECK_THROWS_AS(basis->ComputeCoefficient(throwing, 0), CallbackFailure);
		}
	}

	TEST_CASE("Associated Laguerre parameter and weight domains are checked", "[OrthogonalBasis][basis-contract][domains]")
	{
		for (Real alpha : {Real{-1}, Constants::NegInf, Constants::PosInf, std::numeric_limits<Real>::quiet_NaN()})
			CHECK_THROWS_AS(AssociatedLaguerreBasis(alpha), ArgumentError);
		AssociatedLaguerreBasis basis(2);
		CHECK_THROWS_AS(basis.WeightFunction(-1), DomainError);
		CHECK_THROWS_AS(basis.WeightFunction(std::numeric_limits<Real>::quiet_NaN()), DomainError);
		CHECK(basis.WeightFunction(0) == Real{0});
		CHECK(AssociatedLaguerreBasis(0).WeightFunction(0) == Real{1});
		CHECK(std::isinf(AssociatedLaguerreBasis(Real(-0.5)).WeightFunction(0)));
		CHECK(basis.WeightFunction(Real(2) * std::numeric_limits<Real>::max_exponent) == Real{0});
		CHECK_THROWS_AS(basis.Evaluate(0, std::numeric_limits<Real>::quiet_NaN()), DomainError);
		CHECK_THROWS_AS(AssociatedLaguerreBasis(std::numeric_limits<Real>::max_exponent).Normalization(0), NumericalMethodError);
		if constexpr (std::numeric_limits<Real>::max_exponent10 > 135)
			CHECK(AssociatedLaguerreBasis(100).WeightFunction(1000) ==
				Catch::Approx(5.07595889754945676529e-135).epsilon(Precision::DefaultToleranceStrict));
	}
}
