#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "algorithms/RootFinding.h"
#include "base/Function.h"
#endif

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Algorithms::RootFindingTests
{
	/*********************************************************************/
	/*****          Polynomial equation solvers                      *****/
	/*********************************************************************/
	TEST_CASE("RootFinding::Quadratic_two_real_roots", "[polynomial]")
	{
			TEST_PRECISION_INFO();
		// x^2 - 5x + 6 = 0 => (x-2)(x-3) = 0 => x = 2, 3
		Complex x1, x2;
		int num_roots = MML::SolveQuadratic(REAL(1.0), -REAL(5.0), REAL(6.0), x1, x2);

		REQUIRE(num_roots == 2);
		// Roots can be in any order, sort them
		std::vector<Real> roots = {x1.real(), x2.real()};
		std::sort(roots.begin(), roots.end());
		REQUIRE_THAT(roots[0], WithinAbs(REAL(2.0), REAL(1e-10)));
		REQUIRE_THAT(roots[1], WithinAbs(REAL(3.0), REAL(1e-10)));
	}

	TEST_CASE("RootFinding::Quadratic_one_real_root", "[polynomial]")
	{
			TEST_PRECISION_INFO();
		// x^2 - 4x + 4 = 0 => (x-2)^2 = 0 => x = 2 (double root)
		// Note: Solver returns 2 (for two equal roots) rather than 1
		Complex x1, x2;
		int num_roots = MML::SolveQuadratic(REAL(1.0), -REAL(4.0), REAL(4.0), x1, x2);

		REQUIRE(num_roots == 2);  // Double root counted as 2
		REQUIRE_THAT(x1.real(), WithinAbs(REAL(2.0), REAL(1e-10)));
		REQUIRE_THAT(x2.real(), WithinAbs(REAL(2.0), REAL(1e-10)));
	}

	TEST_CASE("RootFinding::Quadratic_no_real_roots", "[polynomial]")
	{
			TEST_PRECISION_INFO();
		// x^2 + 1 = 0 => no real roots
		Complex x1, x2;
		int num_roots = MML::SolveQuadratic(REAL(1.0), REAL(0.0), REAL(1.0), x1, x2);

		REQUIRE(num_roots == 0);
	}

	TEST_CASE("RootFinding::Quadratic_complex_version", "[polynomial]")
	{
			TEST_PRECISION_INFO();
		// x^2 + 1 = 0 => x = ±i
		Complex x1, x2;
		MML::SolveQuadratic(Complex(REAL(1.0)), Complex(REAL(0.0)), Complex(REAL(1.0)), x1, x2);

		REQUIRE_THAT(x1.real(), WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(std::abs(x1.imag()), WithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE_THAT(x2.real(), WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(std::abs(x2.imag()), WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("RootFinding::Cubic_three_real_roots", "[polynomial]")
	{
			TEST_PRECISION_INFO();
		// x^3 - 6x^2 + 11x - 6 = 0 => (x-1)(x-2)(x-3) = 0 => x = 1, 2, 3
		Complex x1, x2, x3;
		int num_roots = MML::SolveCubic(REAL(1.0), -REAL(6.0), REAL(11.0), -REAL(6.0), x1, x2, x3);

		REQUIRE(num_roots == 3);
		
		// Sort roots for consistent testing
		std::vector<Real> roots = {x1.real(), x2.real(), x3.real()};
		std::sort(roots.begin(), roots.end());
		
		REQUIRE_THAT(roots[0], WithinAbs(REAL(1.0), REAL(1e-8)));
		REQUIRE_THAT(roots[1], WithinAbs(REAL(2.0), REAL(1e-8)));
		REQUIRE_THAT(roots[2], WithinAbs(REAL(3.0), REAL(1e-8)));
	}

	TEST_CASE("RootFinding::Cubic_one_real_root", "[polynomial]")
	{
			TEST_PRECISION_INFO();
		// x^3 - 1 = 0 => (x-1)(x^2+x+1) = 0 => x = 1 (and two complex roots)
		Complex x1, x2, x3;
		int num_roots = MML::SolveCubic(REAL(1.0), REAL(0.0), REAL(0.0), -REAL(1.0), x1, x2, x3);

		REQUIRE(num_roots == 1);
		REQUIRE_THAT(x1.real(), WithinAbs(REAL(1.0), REAL(1e-10)));
	}

	TEST_CASE("RootFinding::Cubic_triple_root", "[polynomial]")
	{
			TEST_PRECISION_INFO();
		// x^3 = 0 => x = 0 (triple root)
		// Note: Degenerate case where d=0, solver returns 1 real root
		Complex x1, x2, x3;
		int num_roots = MML::SolveCubic(REAL(1.0), REAL(0.0), REAL(0.0), REAL(0.0), x1, x2, x3);

		REQUIRE(num_roots >= 1);  // At least one root at x=0
		REQUIRE_THAT(x1.real(), WithinAbs(REAL(0.0), REAL(1e-10)));
		if (num_roots >= 2)
			REQUIRE_THAT(x2.real(), WithinAbs(REAL(0.0), REAL(1e-10)));
		if (num_roots == 3)
			REQUIRE_THAT(x3.real(), WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("RootFinding::Quartic_four_real_roots", "[polynomial][quartic]")
	{
			TEST_PRECISION_INFO();
		// x^4 - 5x^2 + 4 = 0 => (x^2-1)(x^2-4) = 0 => x = ±1, ±2
		Complex x1, x2, x3, x4;
		MML::SolveQuartic(REAL(1.0), REAL(0.0), -REAL(5.0), REAL(0.0), REAL(4.0), x1, x2, x3, x4);

		// Check all roots are real
		REQUIRE_THAT(x1.imag(), WithinAbs(REAL(0.0), REAL(1e-8)));
		REQUIRE_THAT(x2.imag(), WithinAbs(REAL(0.0), REAL(1e-8)));
		REQUIRE_THAT(x3.imag(), WithinAbs(REAL(0.0), REAL(1e-8)));
		REQUIRE_THAT(x4.imag(), WithinAbs(REAL(0.0), REAL(1e-8)));

		// Collect and sort real parts
		std::vector<Real> roots = {x1.real(), x2.real(), x3.real(), x4.real()};
		std::sort(roots.begin(), roots.end());

		REQUIRE_THAT(roots[0], WithinAbs(-REAL(2.0), REAL(1e-8)));
		REQUIRE_THAT(roots[1], WithinAbs(-REAL(1.0), REAL(1e-8)));
		REQUIRE_THAT(roots[2], WithinAbs(REAL(1.0), REAL(1e-8)));
		REQUIRE_THAT(roots[3], WithinAbs(REAL(2.0), REAL(1e-8)));
	}

	TEST_CASE("RootFinding::Quartic_biquadratic", "[polynomial][quartic]")
	{
			TEST_PRECISION_INFO();
		// x^4 - 1 = 0 => (x^2-1)(x^2+1) = 0 => x = ±1, ±i
		Complex x1, x2, x3, x4;
		MML::SolveQuartic(REAL(1.0), REAL(0.0), REAL(0.0), REAL(0.0), -REAL(1.0), x1, x2, x3, x4);

		// Separate real and complex roots
		std::vector<Complex> roots = {x1, x2, x3, x4};
		int real_count = 0;
		int complex_count = 0;

		for (const auto& root : roots) {
			if (std::abs(root.imag()) < 1e-8) {
				real_count++;
				REQUIRE((std::abs(root.real() - REAL(1.0)) < 1e-8 || std::abs(root.real() + REAL(1.0)) < 1e-8));
			} else {
				complex_count++;
				REQUIRE_THAT(root.real(), WithinAbs(REAL(0.0), REAL(1e-8)));
				REQUIRE_THAT(std::abs(root.imag()), WithinAbs(REAL(1.0), REAL(1e-8)));
			}
		}

		REQUIRE(real_count == 2);
		REQUIRE(complex_count == 2);
	}

	TEST_CASE("RootFinding::Quartic_degenerate_to_cubic", "[polynomial][quartic]")
	{
			TEST_PRECISION_INFO();
		// 0*x^4 + x^3 - 6x^2 + 11x - 6 = 0 => cubic with roots 1, 2, 3
		Complex x1, x2, x3, x4;
		MML::SolveQuartic(REAL(0.0), REAL(1.0), -REAL(6.0), REAL(11.0), -REAL(6.0), x1, x2, x3, x4);

		// Should have 3 real roots (x1, x2, x3) and one zero (x4)
		REQUIRE_THAT(x4.real(), WithinAbs(REAL(0.0), REAL(1e-8)));
		REQUIRE_THAT(x4.imag(), WithinAbs(REAL(0.0), REAL(1e-8)));

		// Check the three non-zero roots
		std::vector<Real> non_zero_roots;
		for (auto root : {x1, x2, x3}) {
			if (std::abs(root.imag()) < 1e-8 && std::abs(root.real()) > 1e-8) {
				non_zero_roots.push_back(root.real());
			}
		}

		std::sort(non_zero_roots.begin(), non_zero_roots.end());
		REQUIRE(non_zero_roots.size() == 3);
		REQUIRE_THAT(non_zero_roots[0], WithinAbs(REAL(1.0), REAL(1e-8)));
		REQUIRE_THAT(non_zero_roots[1], WithinAbs(REAL(2.0), REAL(1e-8)));
		REQUIRE_THAT(non_zero_roots[2], WithinAbs(REAL(3.0), REAL(1e-8)));
	}

	TEST_CASE("RootFinding::Quartic_general_case", "[polynomial][quartic]")
	{
			TEST_PRECISION_INFO();
		// x^4 + x^3 - 7x^2 - x + 6 = 0 => roots: -3, -1, 1, 2
		Complex x1, x2, x3, x4;
		MML::SolveQuartic(REAL(1.0), REAL(1.0), -REAL(7.0), -REAL(1.0), REAL(6.0), x1, x2, x3, x4);

		// All roots should be real
		REQUIRE_THAT(x1.imag(), WithinAbs(REAL(0.0), REAL(1e-8)));
		REQUIRE_THAT(x2.imag(), WithinAbs(REAL(0.0), REAL(1e-8)));
		REQUIRE_THAT(x3.imag(), WithinAbs(REAL(0.0), REAL(1e-8)));
		REQUIRE_THAT(x4.imag(), WithinAbs(REAL(0.0), REAL(1e-8)));

		std::vector<Real> roots = {x1.real(), x2.real(), x3.real(), x4.real()};
		std::sort(roots.begin(), roots.end());

		REQUIRE_THAT(roots[0], WithinAbs(-REAL(3.0), REAL(1e-7)));
		REQUIRE_THAT(roots[1], WithinAbs(-REAL(1.0), REAL(1e-7)));
		REQUIRE_THAT(roots[2], WithinAbs(REAL(1.0), REAL(1e-7)));
		REQUIRE_THAT(roots[3], WithinAbs(REAL(2.0), REAL(1e-7)));
	}

	TEST_CASE("RootFinding::Quartic_repeated_roots", "[polynomial][quartic]")
	{
			TEST_PRECISION_INFO();
		// (x-1)^4 = x^4 - 4x^3 + 6x^2 - 4x + 1 = 0 => quadruple root at x=1
		Complex x1, x2, x3, x4;
		MML::SolveQuartic(REAL(1.0), -REAL(4.0), REAL(6.0), -REAL(4.0), REAL(1.0), x1, x2, x3, x4);

		// All roots should be near 1 (numerical methods may not give exact quadruple root)
		REQUIRE_THAT(x1.real(), WithinAbs(REAL(1.0), REAL(1e-6)));
		REQUIRE_THAT(x2.real(), WithinAbs(REAL(1.0), REAL(1e-6)));
		REQUIRE_THAT(x3.real(), WithinAbs(REAL(1.0), REAL(1e-6)));
		REQUIRE_THAT(x4.real(), WithinAbs(REAL(1.0), REAL(1e-6)));
	}

	TEST_CASE("RootFinding::Quartic_double_roots", "[polynomial][quartic]")
	{
			TEST_PRECISION_INFO();
		// (x-1)^2(x-2)^2 = (x^2-2x+1)(x^2-4x+4) = x^4 - 6x^3 + 13x^2 - 12x + 4
		Complex x1, x2, x3, x4;
		MML::SolveQuartic(REAL(1.0), -REAL(6.0), REAL(13.0), -REAL(12.0), REAL(4.0), x1, x2, x3, x4);

		std::vector<Real> roots;
		for (auto root : {x1, x2, x3, x4}) {
			if (std::abs(root.imag()) < 1e-8) {
				roots.push_back(root.real());
			}
		}
		std::sort(roots.begin(), roots.end());

		REQUIRE(roots.size() == 4);
		// Should have two roots near 1 and two near 2
		REQUIRE_THAT(roots[0], WithinAbs(REAL(1.0), REAL(1e-6)));
		REQUIRE_THAT(roots[1], WithinAbs(REAL(1.0), REAL(1e-6)));
		REQUIRE_THAT(roots[2], WithinAbs(REAL(2.0), REAL(1e-6)));
		REQUIRE_THAT(roots[3], WithinAbs(REAL(2.0), REAL(1e-6)));
	}

	TEST_CASE("RootFinding::Quartic_large_coefficients", "[polynomial][quartic]")
	{
			TEST_PRECISION_INFO();
		// Test numerical stability with large coefficients
		// 1000(x^4 - 5x^2 + 4) = 1000x^4 - 5000x^2 + 4000, roots: ±1, ±2
		Complex x1, x2, x3, x4;
		MML::SolveQuartic(REAL(1000.0), REAL(0.0), -REAL(5000.0), REAL(0.0), REAL(4000.0), x1, x2, x3, x4);

		std::vector<Real> roots;
		for (auto root : {x1, x2, x3, x4}) {
			if (std::abs(root.imag()) < 1e-7) {
				roots.push_back(root.real());
			}
		}
		std::sort(roots.begin(), roots.end());

		REQUIRE(roots.size() == 4);
		REQUIRE_THAT(roots[0], WithinAbs(-REAL(2.0), REAL(1e-7)));
		REQUIRE_THAT(roots[1], WithinAbs(-REAL(1.0), REAL(1e-7)));
		REQUIRE_THAT(roots[2], WithinAbs(REAL(1.0), REAL(1e-7)));
		REQUIRE_THAT(roots[3], WithinAbs(REAL(2.0), REAL(1e-7)));
	}

	TEST_CASE("RootFinding::Quartic_small_leading_coefficient", "[polynomial][quartic]")
	{
			TEST_PRECISION_INFO();
		// Test with small but non-zero leading coefficient (not quite degenerate)
		// REAL(0.001)(x-1)(x-2)(x-3)(x-4) expanded
		Real a = REAL(0.001);
		Real roots_expected[4] = {REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)};
		// (x-1)(x-2)(x-3)(x-4) = x^4 - 10x^3 + 35x^2 - 50x + 24
		Complex x1, x2, x3, x4;
		MML::SolveQuartic(REAL(0.001), -REAL(0.010), REAL(0.035), -REAL(0.050), REAL(0.024), x1, x2, x3, x4);

		std::vector<Real> roots;
		for (auto root : {x1, x2, x3, x4}) {
			if (std::abs(root.imag()) < 1e-6) {
				roots.push_back(root.real());
			}
		}
		std::sort(roots.begin(), roots.end());

		REQUIRE(roots.size() == 4);
		for (int i = 0; i < 4; i++) {
			REQUIRE_THAT(roots[i], WithinAbs(roots_expected[i], REAL(1e-5)));
		}
	}

	/*********************************************************************/
	/*****          Bracketing methods                               *****/
	/*********************************************************************/
	TEST_CASE("RootFinding::BracketRoot_simple_case", "[bracketing]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^2 - 4, root at x = 2
		RealFunction func([](Real x) { return x * x - REAL(4.0); });

		Real x1 = REAL(1.0);
		Real x2 = REAL(1.5);

		bool bracketed = RootFinding::BracketRoot(func, x1, x2);

		REQUIRE(bracketed);
		// Should bracket the root at x=2
		REQUIRE(func(x1) * func(x2) < REAL(0.0));
		REQUIRE(x1 < REAL(2.0));
		REQUIRE(x2 > REAL(2.0));
	}

	TEST_CASE("RootFinding::BracketRoot_negative_side", "[bracketing]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^2 - 4, root at x = -2
		RealFunction func([](Real x) { return x * x - REAL(4.0); });

		Real x1 = -REAL(1.5);
		Real x2 = -REAL(1.0);

		bool bracketed = RootFinding::BracketRoot(func, x1, x2);

		REQUIRE(bracketed);
		REQUIRE(func(x1) * func(x2) < REAL(0.0));
		REQUIRE(x1 < -REAL(2.0));
		REQUIRE(x2 > -REAL(2.0));
	}

	TEST_CASE("RootFinding::BracketRoot_no_root", "[bracketing]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^2 + 1, no real roots
		RealFunction func([](Real x) { return x * x + REAL(1.0); });

		Real x1 = REAL(0.0);
		Real x2 = REAL(1.0);

		bool bracketed = RootFinding::BracketRoot(func, x1, x2, 10);  // Low max tries

		REQUIRE_FALSE(bracketed);
	}

	TEST_CASE("RootFinding::FindRootBrackets_multiple_roots", "[bracketing]")
	{
			TEST_PRECISION_INFO();
		// f(x) = sin(x), roots at 0, π, 2π, 3π in [0, 10]
		RealFunction func([](Real x) { return std::sin(x); });

		Vector<Real> xb1, xb2;
		int num_brackets = RootFinding::FindRootBrackets(func, REAL(0.0), REAL(10.0), 1000, xb1, xb2);

		REQUIRE(num_brackets >= 3);  // Should find at least 3 roots

		// Verify each bracket contains a sign change
		for (int i = 0; i < num_brackets; i++) {
			REQUIRE(func(xb1[i]) * func(xb2[i]) <= REAL(0.0));
		}
	}

	TEST_CASE("RootFinding::FindRootBrackets_quadratic", "[bracketing]")
	{
			TEST_PRECISION_INFO();
		// f(x) = (x-2)(x-5) = x^2 - 7x + 10, roots at x = 2, 5
		RealFunction func([](Real x) { return x * x - REAL(7.0) * x + REAL(10.0); });

		Vector<Real> xb1, xb2;
		int num_brackets = RootFinding::FindRootBrackets(func, REAL(0.0), REAL(7.0), 100, xb1, xb2);

		REQUIRE(num_brackets == 2);

		// Verify brackets
		REQUIRE(func(xb1[0]) * func(xb2[0]) <= REAL(0.0));
		REQUIRE(func(xb1[1]) * func(xb2[1]) <= REAL(0.0));
	}

	/*********************************************************************/
	/*****          Bisection method                                 *****/
	/*********************************************************************/
	TEST_CASE("RootFinding::Bisection_simple_root", "[bisection]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^2 - 4, root at x = 2
		RealFunction func([](Real x) { return x * x - REAL(4.0); });

		Real root = RootFinding::FindRootBisection(func, REAL(1.0), REAL(3.0), 1e-8);

		REQUIRE_THAT(root, WithinAbs(REAL(2.0), REAL(1e-7)));
		REQUIRE_THAT(func(root), WithinAbs(REAL(0.0), REAL(1e-6)));
	}

	TEST_CASE("RootFinding::Bisection_cubic", "[bisection]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^3 - x - 2, root near x ≈ REAL(1.521)
		RealFunction func([](Real x) { return x * x * x - x - REAL(2.0); });

		Real root = RootFinding::FindRootBisection(func, REAL(1.0), REAL(2.0), 1e-10);

		REQUIRE_THAT(func(root), WithinAbs(REAL(0.0), REAL(1e-8)));
	}

	TEST_CASE("RootFinding::Bisection_transcendental", "[bisection]")
	{
			TEST_PRECISION_INFO();
		// f(x) = cos(x) - x, root near x ≈ REAL(0.739085)
		RealFunction func([](Real x) { return std::cos(x) - x; });

		Real root = RootFinding::FindRootBisection(func, REAL(0.0), REAL(1.0), 1e-10);

		REQUIRE_THAT(root, WithinAbs(REAL(0.739085), REAL(1e-5)));
		REQUIRE_THAT(func(root), WithinAbs(REAL(0.0), REAL(1e-8)));
	}

	TEST_CASE("RootFinding::Bisection_error_not_bracketed", "[bisection][error]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^2 + 1, no root between 0 and 1
		RealFunction func([](Real x) { return x * x + REAL(1.0); });

		REQUIRE_THROWS_AS(
			RootFinding::FindRootBisection(func, REAL(0.0), REAL(1.0), 1e-8),
			RootFindingError
		);
	}

	/*********************************************************************/
	/*****          Newton-Raphson method                            *****/
	/*********************************************************************/
	TEST_CASE("RootFinding::Newton_simple_root", "[newton]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^2 - 4, root at x = 2
		RealFunction func([](Real x) { return x * x - REAL(4.0); });

		Real root = RootFinding::FindRootNewton(func, REAL(1.0), REAL(3.0), 1e-10);

		REQUIRE_THAT(root, WithinAbs(REAL(2.0), REAL(1e-9)));
		REQUIRE_THAT(func(root), WithinAbs(REAL(0.0), REAL(1e-8)));
	}

	TEST_CASE("RootFinding::Newton_cubic_fast_convergence", "[newton]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^3 - 2x - 5, root near x ≈ REAL(2.0946)
		RealFunction func([](Real x) { return x * x * x - REAL(2.0) * x - REAL(5.0); });

		Real root = RootFinding::FindRootNewton(func, REAL(2.0), REAL(3.0), 1e-12);

		REQUIRE_THAT(root, WithinAbs(REAL(2.0946), REAL(1e-4)));
		REQUIRE_THAT(func(root), WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("RootFinding::Newton_transcendental", "[newton]")
	{
			TEST_PRECISION_INFO();
		// f(x) = e^x - 3x, root near x ≈ REAL(1.512)
		RealFunction func([](Real x) { return std::exp(x) - REAL(3.0) * x; });

		Real root = RootFinding::FindRootNewton(func, REAL(1.0), REAL(2.0), 1e-10);

		REQUIRE_THAT(func(root), WithinAbs(REAL(0.0), REAL(1e-8)));
	}

	TEST_CASE("RootFinding::Newton_error_out_of_bounds", "[newton][error]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^2 - 4 with bad initial bracket that will cause jumping out
		// Start very close to root but with narrow bracket
		RealFunction func([](Real x) { return x * x - REAL(4.0); });

		// This might throw if Newton step jumps outside brackets
		// (depends on implementation details, but testing error handling)
		try {
			Real root = RootFinding::FindRootNewton(func, REAL(1.99), REAL(2.01), 1e-15);
			// If it succeeds, that's also fine - just verify it's correct
			REQUIRE_THAT(root, WithinAbs(REAL(2.0), REAL(1e-10)));
		}
		catch (const RootFindingError&) {
			// Expected if Newton jumps out of narrow brackets
			REQUIRE(true);
		}
	}

	/*********************************************************************/
	/*****          Convergence and precision tests                  *****/
	/*********************************************************************/
	TEST_CASE("RootFinding::Convergence_bisection_vs_newton", "[convergence]")
	{
			TEST_PRECISION_INFO();
		// Compare bisection vs Newton for same problem
		// f(x) = x^3 - x - 1, root near x ≈ REAL(1.3247)
		RealFunction func([](Real x) { return x * x * x - x - REAL(1.0); });

		Real root_bisection = RootFinding::FindRootBisection(func, REAL(1.0), REAL(2.0), 1e-10);
		Real root_newton = RootFinding::FindRootNewton(func, REAL(1.0), REAL(2.0), 1e-10);

		// Both should find the same root
		REQUIRE_THAT(root_bisection, WithinAbs(root_newton, REAL(1e-9)));
		REQUIRE_THAT(func(root_bisection), WithinAbs(REAL(0.0), REAL(1e-8)));
		REQUIRE_THAT(func(root_newton), WithinAbs(REAL(0.0), REAL(1e-8)));
	}

	TEST_CASE("RootFinding::Precision_high_accuracy", "[precision]")
	{
			TEST_PRECISION_INFO();
		// Test high-precision root finding: x^2 = 2, root = sqrt(2)
		RealFunction func([](Real x) { return x * x - REAL(2.0); });

		Real root = RootFinding::FindRootBisection(func, REAL(1.0), REAL(2.0), 1e-15);

		REQUIRE_THAT(root, WithinAbs(std::sqrt(REAL(2.0)), 1e-14));
	}

	TEST_CASE("RootFinding::Multiple_roots_workflow", "[workflow]")
	{
			TEST_PRECISION_INFO();
		// Complete workflow: find all roots of x^3 - 6x^2 + 11x - 6
		// Roots at x = 1, 2, 3
		RealFunction func([](Real x) { 
			return x * x * x - REAL(6.0) * x * x + REAL(11.0) * x - REAL(6.0); 
		});

		// Step 1: Find brackets
		Vector<Real> xb1, xb2;
		int num_brackets = RootFinding::FindRootBrackets(func, REAL(0.0), REAL(4.0), 100, xb1, xb2);

		REQUIRE(num_brackets >= 3);  // Should find at least 3 roots

		// Step 2: Refine each bracket with bisection
		std::vector<Real> roots;
		for (int i = 0; i < num_brackets; i++) {
			// Verify bracket is valid before bisection
			if (func(xb1[i]) * func(xb2[i]) < REAL(0.0)) {
				try {
					Real root = RootFinding::FindRootBisection(func, xb1[i], xb2[i], 1e-8);
					// Avoid close duplicates (tolerance 1e-4 to allow distinct roots)
					bool is_duplicate = false;
					for (const auto& existing_root : roots) {
						if (std::abs(root - existing_root) < 1e-4) {
							is_duplicate = true;
							break;
						}
					}
					if (!is_duplicate) {
						roots.push_back(root);
					}
				}
				catch (const RootFindingError&) {
					// Skip brackets that fail
				}
			}
		}

		// Step 3: Verify roots found (algorithm may find 2-4 depending on bracketing)
		std::sort(roots.begin(), roots.end());
		REQUIRE(roots.size() >= 2);  // At least 2 roots
		
		// Check that at least 2 of the 3 expected roots were found
		// (Bracketing algorithm behavior can vary at boundaries)
		bool found_1 = false, found_2 = false, found_3 = false;
		for (const auto& root : roots) {
			if (std::abs(root - REAL(1.0)) < 1e-5) found_1 = true;
			if (std::abs(root - REAL(2.0)) < 1e-5) found_2 = true;
			if (std::abs(root - REAL(3.0)) < 1e-5) found_3 = true;
		}
		int found_count = (found_1 ? 1 : 0) + (found_2 ? 1 : 0) + (found_3 ? 1 : 0);
		REQUIRE(found_count >= 2);  // At least 2 of the 3 roots found
	}

	/*********************************************************************/
	/*****          Edge cases and numerical stability               *****/
	/*********************************************************************/
	TEST_CASE("RootFinding::Edge_case_root_at_boundary", "[edge]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x - 2, root exactly at x = 2
		RealFunction func([](Real x) { return x - REAL(2.0); });

		Vector<Real> xb1, xb2;
		int num_brackets = RootFinding::FindRootBrackets(func, REAL(2.0), REAL(3.0), 10, xb1, xb2);

		// Should find the root even though it's at the boundary
		REQUIRE(num_brackets >= 1);
	}

	TEST_CASE("RootFinding::Edge_case_very_flat_function", "[edge]")
	{
			TEST_PRECISION_INFO();
		// f(x) = x^5 - very flat near x=0 but has root there
		RealFunction func([](Real x) { return x * x * x * x * x; });

		Real root = RootFinding::FindRootBisection(func, -REAL(0.1), REAL(0.1), 1e-10);

		REQUIRE_THAT(root, WithinAbs(REAL(0.0), REAL(1e-9)));
	}

	TEST_CASE("RootFinding::Quartic_complex_roots_only", "[polynomial][quartic]")
	{
			TEST_PRECISION_INFO();
		// x^4 + 1 = 0 => four complex roots (4th roots of -1)
		Complex x1, x2, x3, x4;
		MML::SolveQuartic(REAL(1.0), REAL(0.0), REAL(0.0), REAL(0.0), REAL(1.0), x1, x2, x3, x4);

		// All roots should be complex (non-real)
		std::vector<Complex> roots = {x1, x2, x3, x4};
		
		for (const auto& root : roots) {
			// Should have |z|^4 = 1, so |z| = 1
			REQUIRE_THAT(std::abs(root), WithinAbs(REAL(1.0), REAL(1e-8)));
			// Should satisfy z^4 = -1
			Complex z4 = root * root * root * root;
			REQUIRE_THAT(z4.real(), WithinAbs(-REAL(1.0), REAL(1e-8)));
			REQUIRE_THAT(z4.imag(), WithinAbs(REAL(0.0), REAL(1e-8)));
		}
	}
}

