#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/algorithms/RootFinding/RootFindingPolynoms.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace {
	using namespace MML;

	void CheckRoots(std::vector<Complex> roots, std::vector<Complex> expected)
	{
		REQUIRE(roots.size() == expected.size());
		for (const auto& wanted : expected) {
			const auto closest = std::min_element(roots.begin(), roots.end(), [wanted](Complex a, Complex b) {
				return std::abs(a - wanted) < std::abs(b - wanted);
			});
			REQUIRE(closest != roots.end());
			CHECK(std::isfinite(std::abs(*closest)));
			// Root reconstruction includes several rounded cubic/quartic operations.
			CHECK(std::abs(*closest - wanted) <= std::max(Precision::LooseTolerance, 16 * Precision::Eps) *
				std::max(Real(1), std::abs(wanted)));
			roots.erase(closest);
		}
	}

	std::array<Real, 6> Scales()
	{
		const Real large = std::numeric_limits<Real>::max() / 32;
		return {Real(1), Real(1e-20), Real(-1e-20), large, Real(1) / large, -large};
	}

	// Exercise runtime inputs without MSVC C4723 from unreachable constant-folded branches.
	auto* volatile quadratic = static_cast<int (*)(Real, Real, Real, Complex&, Complex&)>(&SolveQuadratic);
	auto* volatile complexQuadratic = static_cast<void (*)(const Complex&, const Complex&, const Complex&, Complex&, Complex&)>(&SolveQuadratic);
	auto* volatile cubic = &SolveCubic;
	auto* volatile quartic = &SolveQuartic;
}

TEST_CASE("Quadratic degree and stable roots survive coefficient scaling", "[RootFinding][Polynomial][polynomial-scale]")
{
	for (Real scale : Scales()) {
		CAPTURE(scale);
		Complex first, second;
		CHECK(quadratic(scale, 0, -scale, first, second) == 2);
		CheckRoots({first, second}, {-1, 1});
		CHECK(quadratic(scale, -3 * scale, 2 * scale, first, second) == 2);
		CheckRoots({first, second}, {1, 2});
		CHECK(quadratic(scale, 0, scale, first, second) == 0);
		CheckRoots({first, second}, {Complex{0, -1}, Complex{0, 1}});
		CHECK(quadratic(0, scale, -2 * scale, first, second) == 1);
		CHECK(std::abs(first - Complex{2}) <= Precision::EqualityTolerance);
	}
}

TEST_CASE("Complex quadratic degree survives nonzero coefficient scaling", "[RootFinding][Polynomial][polynomial-scale]")
{
	for (Real scale : Scales()) {
		const Complex factor{scale, scale / 2};
		Complex first, second;
		complexQuadratic(factor, Complex{}, -factor, first, second);
		CheckRoots({first, second}, {-1, 1});
		complexQuadratic(Complex{}, factor, -Real(2) * factor, first, second);
		CHECK(std::abs(first - Complex{2}) <= Precision::EqualityTolerance);
	}
}

TEST_CASE("Cubic and quartic reductions preserve scaled polynomial degree", "[RootFinding][Polynomial][polynomial-scale]")
{
	for (Real scale : Scales()) {
		CAPTURE(scale);
		Complex a, b, c, d;
		CHECK(cubic(scale, -6 * scale, 11 * scale, -6 * scale, a, b, c) == 3);
		CheckRoots({a, b, c}, {1, 2, 3});
		CHECK(cubic(0, scale, 0, -scale, a, b, c) == 2);
		CheckRoots({a, b}, {-1, 1});
		quartic(scale, 0, -5 * scale, 0, 4 * scale, a, b, c, d);
		CheckRoots({a, b, c, d}, {-2, -1, 1, 2});
		quartic(0, scale, -6 * scale, 11 * scale, -6 * scale, a, b, c, d);
		CheckRoots({a, b, c}, {1, 2, 3});
	}
}

TEST_CASE("General polynomial methods are independent of coefficient units", "[RootFinding][Polynomial][polynomial-scale]")
{
	for (auto method : {RootFinding::PolynomialRootMethod::Laguerre,
		RootFinding::PolynomialRootMethod::CompanionMatrix, RootFinding::PolynomialRootMethod::Bairstow}) {
		for (Real scale : Scales()) {
			CAPTURE(RootFinding::ToString(method), scale);
			RootFinding::PolynomialRootConfig config;
			config.method = method;
			config.max_iterations = 500;
			config.polish = false;
			const PolynomReal polynomial{4 * scale, 0, -5 * scale, 0, scale};
			const auto result = RootFinding::FindPolynomialRootsDetailed(polynomial, config);
			INFO(result.error_message);
			REQUIRE(result.converged);
			CHECK(result.polynomial_degree == 4);
			CHECK(result.max_relative_residual <= config.tolerance);
			std::vector<Complex> roots;
			for (const auto& root : result.roots) roots.push_back(root.root);
			CheckRoots(roots, {-2, -1, 1, 2});
		}
	}
}

TEST_CASE("Original polynomial polishing and verification retain coefficient scale", "[RootFinding][Polynomial][polynomial-scale]")
{
	for (Real scale : Scales()) {
		CAPTURE(scale);
		const PolynomReal polynomial{2 * scale, -3 * scale, scale};
		const Vector<Complex> estimates{Complex{Real(0.9)}, Complex{Real(2.1)}};
		const auto polished = RootFinding::PolishPolynomialRoots(polynomial, estimates);
		INFO(polished.error_message);
		REQUIRE(polished.converged);
		std::vector<Complex> roots;
		for (const auto& root : polished.roots) {
			roots.push_back(root.root);
			CHECK(root.relative_residual <= Precision::PolynomialRootTolerance);
		}
		CheckRoots(roots, {1, 2});
		RootFinding::PolynomialRootConfig config;
		config.polish = false;
		const auto rejected = RootFinding::PolishPolynomialRoots(polynomial, Vector<Complex>{Complex{}, Complex{}}, config);
		REQUIRE_FALSE(rejected.converged);
		REQUIRE(rejected.roots.size() == 2);
		CHECK(rejected.max_relative_residual == Real(1));
		CHECK_FALSE(rejected.error_message.empty());
	}
}

TEST_CASE("Tiny constants and subnormal polynomials are not zero polynomials", "[RootFinding][Polynomial][polynomial-scale]")
{
	const Real tiny = std::numeric_limits<Real>::denorm_min() * 16;
	const auto constant = RootFinding::FindPolynomialRootsDetailed(PolynomReal::Constant(tiny));
	CHECK(constant.status == AlgorithmStatus::Completed);
	const auto result = RootFinding::FindPolynomialRootsDetailed(PolynomReal{-tiny, 0, tiny});
	INFO(result.error_message);
	REQUIRE(result.converged);
	std::vector<Complex> roots;
	for (const auto& root : result.roots) roots.push_back(root.root);
	CheckRoots(roots, {-1, 1});
	Complex first, second;
	CHECK(quadratic(tiny, 0, -tiny, first, second) == 2);
	CheckRoots({first, second}, {-1, 1});
}

TEST_CASE("Companion polynomial solver honors configured iteration budget", "[RootFinding][Polynomial][polynomial-scale]")
{
	const PolynomReal polynomial{-1, 1, 1, 1, 1, 1};
	Matrix<Real> companion(5, 5);
	for (int index = 0; index < 4; ++index) companion(index, index + 1) = 1;
	for (int index = 0; index < 5; ++index) companion(4, index) = -polynomial[index];
	REQUIRE_FALSE(EigenSolver::Solve(companion, Precision::EigenSolverConvergenceTolerance, 1).converged);
	RootFinding::PolynomialRootConfig config;
	config.method = RootFinding::PolynomialRootMethod::CompanionMatrix;
	config.max_iterations = 1;
	const auto limited = RootFinding::FindPolynomialRootsDetailed(polynomial, config);
	CHECK_FALSE(limited.converged);
	CHECK_FALSE(limited.IsSuccess());
	CHECK_FALSE(limited.error_message.empty());
	config.max_iterations = 500;
	CHECK(RootFinding::FindPolynomialRootsDetailed(polynomial, config).converged);
}

TEST_CASE("Raw polynomial methods retain scaled odd and linear terms", "[RootFinding][Polynomial][polynomial-scale]")
{
	for (Real scale : Scales()) {
		CAPTURE(scale);
		const PolynomReal linear{-2 * scale, scale};
		for (const auto& roots : {RootFinding::LaguerreRoots(linear), RootFinding::EigenvalueRoots(linear),
			RootFinding::BairstowRoots(linear)})
			CheckRoots({roots.begin(), roots.end()}, {2});
		const PolynomReal odd{0, 4 * scale, 0, -5 * scale, 0, scale};
		for (const auto& roots : {RootFinding::LaguerreRoots(odd), RootFinding::EigenvalueRoots(odd),
			RootFinding::BairstowRoots(odd)})
			CheckRoots({roots.begin(), roots.end()}, {-2, -1, 0, 1, 2});
	}
}

TEST_CASE("Cardano and Ferrari branches preserve roots under common scaling", "[RootFinding][Polynomial][polynomial-scale]")
{
	for (Real scale : Scales()) {
		CAPTURE(scale);
		Complex a, b, c, d;
		CHECK(cubic(scale, -3 * scale, 3 * scale, -scale, a, b, c) == 3);
		CheckRoots({a, b, c}, {1, 1, 1});
		const Real factor = scale / 1001;
		CHECK(cubic(factor, 3 * factor, 3 * factor, scale, a, b, c) == 1);
		CheckRoots({a, b, c}, {-11, Complex{4, Real(5) * std::sqrt(Real(3))},
			Complex{4, -Real(5) * std::sqrt(Real(3))}});
		quartic(scale, -6 * scale, 3 * scale, 26 * scale, -24 * scale, a, b, c, d);
		CheckRoots({a, b, c, d}, {-2, 1, 3, 4});
	}
}

TEST_CASE("Nonzero leading coefficients never trigger lower degree reduction", "[RootFinding][Polynomial][polynomial-scale]")
{
	const Real small = Real(1e-20);
	Complex a, b, c, d;
	CHECK(quadratic(small, 0, -1, a, b) == 2);
	CHECK(std::abs(a * std::sqrt(small)) == Catch::Approx(1));
	CHECK(std::abs(b * std::sqrt(small)) == Catch::Approx(1));
	CHECK(cubic(small, 0, 0, -1, a, b, c) == 1);
	for (Complex root : {a, b, c})
		CHECK(std::abs(small * root * root * root - Real(1)) <= Precision::LooseTolerance);
	quartic(small, 0, 0, 0, -1, a, b, c, d);
	for (Complex root : {a, b, c, d})
		CHECK(std::abs(small * root * root * root * root - Real(1)) <= Precision::LooseTolerance);
}

TEST_CASE("Stable real and complex quadratic branches retain the small root", "[RootFinding][Polynomial][polynomial-scale]")
{
	const Real separation = Real(1) / std::sqrt(Precision::Eps);
	for (Real scale : Scales()) {
		// Keep every input coefficient finite even at the largest common scale.
		// A subnormal factor quantizes the *input polynomial*, changing the expected
		// root by much more than epsilon. Retain wide scaling with normal coefficients.
		const Real resolvedScale = std::max(std::abs(scale), std::numeric_limits<Real>::min() * separation);
		const Real signedScale = std::copysign(resolvedScale, scale);
		const Real factor = signedScale / separation;
		Complex first, second;
		quadratic(factor, -signedScale, factor, first, second);
		if (std::abs(first) > std::abs(second)) std::swap(first, second);
		CHECK(std::abs(first * separation - Real(1)) <= Precision::LooseTolerance);
		CHECK(std::abs(second / separation - Real(1)) <= Precision::LooseTolerance);
		complexQuadratic(Complex{factor}, Complex{-signedScale}, Complex{factor}, first, second);
		if (std::abs(first) > std::abs(second)) std::swap(first, second);
		CHECK(std::abs(first * separation - Real(1)) <= Precision::LooseTolerance);
		CHECK(std::abs(second / separation - Real(1)) <= Precision::LooseTolerance);
	}
}

TEST_CASE("Polynomial coefficient validation cannot publish nonfinite roots", "[RootFinding][Polynomial][polynomial-scale]")
{
	for (Real invalid : {std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		Complex a, b, c, d;
		REQUIRE_THROWS_AS(quadratic(invalid, 0, -1, a, b), ArgumentError);
		REQUIRE_THROWS_AS(complexQuadratic(Complex{1}, Complex{0, invalid}, Complex{1}, a, b), ArgumentError);
		REQUIRE_THROWS_AS(cubic(1, invalid, 0, -1, a, b, c), ArgumentError);
		REQUIRE_THROWS_AS(quartic(1, 0, invalid, 0, 1, a, b, c, d), ArgumentError);
		const auto result = RootFinding::FindPolynomialRootsDetailed(PolynomReal{invalid, 0, 1});
		CHECK(result.status == AlgorithmStatus::InvalidInput);
		CHECK_FALSE(result.error_message.empty());
	}
	const PolynomReal storedZeroLeading{1, 1, 0};
	CHECK(RootFinding::FindPolynomialRootsDetailed(storedZeroLeading).status == AlgorithmStatus::InvalidInput);
	const auto nonfiniteEstimate = RootFinding::PolishPolynomialRoots(
		PolynomReal{-1, 0, 1}, Vector<Complex>{Complex{std::numeric_limits<Real>::infinity()}, Complex{1}});
	CHECK_FALSE(nonfiniteEstimate.converged);
	CHECK(nonfiniteEstimate.status == AlgorithmStatus::NumericalInstability);
	CHECK(std::isinf(nonfiniteEstimate.max_relative_residual));
	CHECK_FALSE(nonfiniteEstimate.error_message.empty());
	const PolynomReal mixed{std::numeric_limits<Real>::denorm_min(), std::numeric_limits<Real>::max()};
	const auto unsupported = RootFinding::FindPolynomialRootsDetailed(mixed);
	CHECK_FALSE(unsupported.IsSuccess());
	CHECK_FALSE(unsupported.error_message.empty());
	REQUIRE_THROWS_AS(RootFinding::EigenvalueRoots(PolynomReal{-1, 0, 1}, 0), ArgumentError);
}
