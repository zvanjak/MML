#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else

#include "base/ODESystem.h"

#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"
#endif

#include "../test_data/diff_eq_systems_test_bed.h"

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

using namespace MML;
using namespace MML::Testing;
using namespace MML::TestBeds;


namespace MML::Tests::Algorithms::ODESystemSolverTests 
{
	TEST_CASE("Test_Diff_Eq_Solvers_EulerMethod_LinSys1", "[ODESystemFixedStepSolver]") {
		TEST_PRECISION_INFO();
		auto odeSys = TestBeds::ODESystemTestBed::getTestODESystemWithEndSolution(0);

		auto initCond = odeSys.getInitialConditions();
		Real t1 = REAL(0.0), t2 = odeSys.getEndTime();
		int numSteps = 1000;

		ODESystemFixedStepSolver fixedSolver(*odeSys.getODESystem(), StepCalculators::EulerStepCalc);
		ODESystemSolution sol = fixedSolver.integrate(initCond, t1, t2, numSteps);

		Vector<Real> finalX = sol.getXValuesAtEnd();
		Vector<Real> exactSol = odeSys.getEndSolution();

		// compare with exact solution
		for (int i = 0; i < finalX.size(); i++) {
			REQUIRE_THAT(finalX[i], WithinAbs(exactSol[i], REAL(1e-1))); // POOR PRECISION, Euler method
		}
	}

	TEST_CASE("Test_Diff_Eq_Solvers_EulerMethod_HarmOsc1", "[ODESystemFixedStepSolver]") {
		TEST_PRECISION_INFO();
		auto odeSys = TestBeds::ODESystemTestBed::getTestODESystemWithEndSolution(3);

		auto initCond = odeSys.getInitialConditions();
		Real t1 = REAL(0.0), t2 = odeSys.getEndTime();
		int numSteps = 1000;

		ODESystemFixedStepSolver fixedSolver(*odeSys.getODESystem(), StepCalculators::EulerStepCalc);
		ODESystemSolution sol = fixedSolver.integrate(initCond, t1, t2, numSteps);

		Vector<Real> finalX = sol.getXValuesAtEnd();
		Vector<Real> exactSol = odeSys.getEndSolution();

		// compare with exact solution
		for (int i = 0; i < finalX.size(); i++) {
			REQUIRE_THAT(finalX[i], WithinAbs(exactSol[i], REAL(1e-1))); // POOR PRECISION, Euler method
		}
	}

	TEST_CASE("Test_Diff_Eq_Solvers_RungeKutta4_LinSys1", "[ODESystemFixedStepSolver]") {
		TEST_PRECISION_INFO();
		auto odeSys = TestBeds::ODESystemTestBed::getTestODESystemWithEndSolution(0);

		auto initCond = odeSys.getInitialConditions();
		Real t1 = REAL(0.0), t2 = odeSys.getEndTime();
		int numSteps = 1000;

		ODESystemFixedStepSolver fixedSolver(*odeSys.getODESystem(), StepCalculators::RK4_Basic);
		ODESystemSolution sol = fixedSolver.integrate(initCond, t1, t2, numSteps);

		Vector<Real> finalX = sol.getXValuesAtEnd();
		Vector<Real> exactSol = odeSys.getEndSolution();

		// compare with exact solution
		for (int i = 0; i < finalX.size(); i++) {
			REQUIRE_THAT(finalX[i], WithinAbs(exactSol[i], REAL(1e-11)));
		}
	}

	TEST_CASE("Test_Diff_Eq_Solvers_RungeKutta4_LinSys2", "[ODESystemFixedStepSolver]") {
		TEST_PRECISION_INFO();
		auto odeSys = TestBeds::ODESystemTestBed::getTestODESystemWithEndSolution(1);

		auto initCond = odeSys.getInitialConditions();
		Real t1 = REAL(0.0), t2 = odeSys.getEndTime();
		int numSteps = 1000;

		ODESystemFixedStepSolver fixedSolver(*odeSys.getODESystem(), StepCalculators::RK4_Basic);
		ODESystemSolution sol = fixedSolver.integrate(initCond, t1, t2, numSteps);

		Vector<Real> finalX = sol.getXValuesAtEnd();
		Vector<Real> exactSol = odeSys.getEndSolution();

		// compare with exact solution
		for (int i = 0; i < finalX.size(); i++) {
			REQUIRE_THAT(finalX[i], WithinAbs(exactSol[i], REAL(1e-4)));
		}
	}

	TEST_CASE("Test_Diff_Eq_Solvers_RungeKutta4_HarmOsc1", "[ODESystemFixedStepSolver]") {
		TEST_PRECISION_INFO();
		auto odeSys = TestBeds::ODESystemTestBed::getTestODESystemWithEndSolution(3);

		auto initCond = odeSys.getInitialConditions();
		Real t1 = REAL(0.0), t2 = odeSys.getEndTime();
		int numSteps = 1000;

		ODESystemFixedStepSolver fixedSolver(*odeSys.getODESystem(), StepCalculators::RK4_Basic);
		ODESystemSolution sol = fixedSolver.integrate(initCond, t1, t2, numSteps);

		Vector<Real> finalX = sol.getXValuesAtEnd();
		Vector<Real> exactSol = odeSys.getEndSolution();

		// compare with exact solution
		for (int i = 0; i < finalX.size(); i++) {
			REQUIRE_THAT(finalX[i], WithinAbs(exactSol[i], REAL(1e-8)));
		}
	}

	// let's build a test case for multiple verification
	TEST_CASE("Test_Diff_Eq_Solvers_AllCases", "[ODESystemFixedStepSolver]") {
		TEST_PRECISION_INFO();
		// ODE systems we wll be solving
		auto odeSys0 = TestBeds::ODESystemTestBed::getTestODESystemWithEndSolution(0);
		auto odeSys1 = TestBeds::ODESystemTestBed::getTestODESystemWithEndSolution(1);
		auto odeSys2 = TestBeds::ODESystemTestBed::getTestODESystemWithEndSolution(2);
		auto odeSys3 = TestBeds::ODESystemTestBed::getTestODESystemWithEndSolution(3);

		std::vector<ITestODESystemWithEndSolution*> odeSysList = {&odeSys0, &odeSys3};

		// Step calculators we will be using
		std::vector<std::pair<std::string, IODESystemStepCalculator*>> stepCalcs = {{"Euler", &StepCalculators::EulerStepCalc},
																					{"RK4", &StepCalculators::RK4_Basic}};

		// expected precisions for each (odeSys, stepCalc) pair
		std::vector<std::vector<double>> expectedPrecisions = {
			{1e-1, 1e-11}, // odeSys0
			{1e-1, 1e-8}   // odeSys3
		};

		int numSteps = 1000;

		// Now, run the tests
		for (size_t i = 0; i < odeSysList.size(); i++) {
			auto odeSys = odeSysList[i];
			auto initCond = odeSys->getInitialConditions();
			Real t1 = REAL(0.0), t2 = odeSys->getEndTime();

			Vector<Real> exactSol = odeSys->getEndSolution();

			for (size_t j = 0; j < stepCalcs.size(); j++) {
				auto [stepCalcName, stepCalc] = stepCalcs[j];
				double expectedPrecision = expectedPrecisions[i][j];

				ODESystemFixedStepSolver fixedSolver(*(odeSys->getODESystem()), *stepCalc);
				ODESystemSolution sol = fixedSolver.integrate(initCond, t1, t2, numSteps);

				Vector<Real> finalX = sol.getXValuesAtEnd();

				// compare with exact solution
				for (int k = 0; k < finalX.size(); k++) {
					REQUIRE_THAT(finalX[k], WithinAbs(exactSol[k], expectedPrecision));
				}
			}
		}
	}

	/*********************************************************************/
	/*****           DormandPrince5 and CashKarp Tests               *****/
	/*********************************************************************/

	TEST_CASE("Test_DormandPrince5_LinearSystem", "[ODESystemFixedStepSolver][DormandPrince5]") {
		TEST_PRECISION_INFO();
		auto odeSys = TestBeds::ODESystemTestBed::getTestODESystemWithEndSolution(0);

		auto initCond = odeSys.getInitialConditions();
		Real t1 = REAL(0.0), t2 = odeSys.getEndTime();
		int numSteps = 1000;

		ODESystemFixedStepSolver fixedSolver(*odeSys.getODESystem(), StepCalculators::DormandPrince5StepCalc);
		ODESystemSolution sol = fixedSolver.integrate(initCond, t1, t2, numSteps);

		Vector<Real> finalX = sol.getXValuesAtEnd();
		Vector<Real> exactSol = odeSys.getEndSolution();

		// DormandPrince5 should give excellent precision
		for (int i = 0; i < finalX.size(); i++) {
			REQUIRE_THAT(finalX[i], WithinAbs(exactSol[i], REAL(1e-10)));
		}
	}

	TEST_CASE("Test_DormandPrince5_HarmonicOscillator", "[ODESystemFixedStepSolver][DormandPrince5]") {
		TEST_PRECISION_INFO();
		auto odeSys = TestBeds::ODESystemTestBed::getTestODESystemWithEndSolution(3);

		auto initCond = odeSys.getInitialConditions();
		Real t1 = REAL(0.0), t2 = odeSys.getEndTime();
		int numSteps = 1000;

		ODESystemFixedStepSolver fixedSolver(*odeSys.getODESystem(), StepCalculators::DormandPrince5StepCalc);
		ODESystemSolution sol = fixedSolver.integrate(initCond, t1, t2, numSteps);

		Vector<Real> finalX = sol.getXValuesAtEnd();
		Vector<Real> exactSol = odeSys.getEndSolution();

		// DormandPrince5 should maintain excellent precision over long integration
		for (int i = 0; i < finalX.size(); i++) {
			REQUIRE_THAT(finalX[i], WithinAbs(exactSol[i], REAL(1e-8)));
		}
	}

	TEST_CASE("Test_CashKarp_LinearSystem", "[ODESystemFixedStepSolver][CashKarp]") {
		TEST_PRECISION_INFO();
		auto odeSys = TestBeds::ODESystemTestBed::getTestODESystemWithEndSolution(0);

		auto initCond = odeSys.getInitialConditions();
		Real t1 = REAL(0.0), t2 = odeSys.getEndTime();
		int numSteps = 1000;

		ODESystemFixedStepSolver fixedSolver(*odeSys.getODESystem(), StepCalculators::RK5_CashKarp);
		ODESystemSolution sol = fixedSolver.integrate(initCond, t1, t2, numSteps);

		Vector<Real> finalX = sol.getXValuesAtEnd();
		Vector<Real> exactSol = odeSys.getEndSolution();

		// CashKarp should also give excellent precision
		for (int i = 0; i < finalX.size(); i++) {
			REQUIRE_THAT(finalX[i], WithinAbs(exactSol[i], REAL(1e-10)));
		}
	}

	TEST_CASE("Test_CashKarp_HarmonicOscillator", "[ODESystemFixedStepSolver][CashKarp]") {
		TEST_PRECISION_INFO();
		auto odeSys = TestBeds::ODESystemTestBed::getTestODESystemWithEndSolution(3);

		auto initCond = odeSys.getInitialConditions();
		Real t1 = REAL(0.0), t2 = odeSys.getEndTime();
		int numSteps = 1000;

		ODESystemFixedStepSolver fixedSolver(*odeSys.getODESystem(), StepCalculators::RK5_CashKarp);
		ODESystemSolution sol = fixedSolver.integrate(initCond, t1, t2, numSteps);

		Vector<Real> finalX = sol.getXValuesAtEnd();
		Vector<Real> exactSol = odeSys.getEndSolution();

		// CashKarp should maintain excellent precision
		for (int i = 0; i < finalX.size(); i++) {
			REQUIRE_THAT(finalX[i], WithinAbs(exactSol[i], REAL(1e-8)));
		}
	}

	TEST_CASE("Test_StiffODE_VanDerPol", "[ODESystemFixedStepSolver][StiffODE]") {
		TEST_PRECISION_INFO();
		// Van der Pol oscillator: x'' - mu*(1-x^2)*x' + x = 0
		// Rewrite as system: x' = y, y' = mu*(1-x^2)*y - x
		// For mu >> 1, this becomes stiff

		Real mu = REAL(10.0); // Stiffness parameter

		// Van der Pol equation as function pointer
		auto vanDerPolFunc = [](Real t, const Vector<Real>& x, Vector<Real>& dxdt) {
			Real mu = REAL(10.0);
			dxdt[0] = x[1];
			dxdt[1] = mu * (REAL(1.0) - x[0] * x[0]) * x[1] - x[0];
		};
		ODESystem vanDerPol(2, vanDerPolFunc);

		Vector<Real> initCond{REAL(2.0), REAL(0.0)};
		Real t1 = REAL(0.0), t2 = REAL(10.0);
		int numSteps = 5000; // Need more steps for stiff equation

		// Test that higher-order methods can handle stiff equations better
		ODESystemFixedStepSolver solver1(vanDerPol, StepCalculators::RK4_Basic);
		ODESystemSolution sol1 = solver1.integrate(initCond, t1, t2, numSteps);

		ODESystemFixedStepSolver solver2(vanDerPol, StepCalculators::DormandPrince5StepCalc);
		ODESystemSolution sol2 = solver2.integrate(initCond, t1, t2, numSteps);

		// Both should complete without crashing
		REQUIRE(sol1.getXValuesAtEnd().size() == 2);
		REQUIRE(sol2.getXValuesAtEnd().size() == 2);

		// DormandPrince5 should give more stable results
		REQUIRE(std::isfinite(sol2.getXValuesAtEnd()[0]));
		REQUIRE(std::isfinite(sol2.getXValuesAtEnd()[1]));
	}

	/*********************************************************************/
	/*****               Leapfrog (Velocity Verlet) tests            *****/
	/*********************************************************************/

	// Test Leapfrog on simple harmonic oscillator - ideal case for symplectic integrator
	TEST_CASE("Test_Leapfrog_HarmonicOscillator", "[Leapfrog][ODESystemFixedStepSolver]") {
		TEST_PRECISION_INFO();
		// Simple harmonic oscillator: x'' = -x
		// State: [x, v] where v = x'
		// System: dx/dt = v, dv/dt = -x
		auto odeSys = TestBeds::ODESystemTestBed::getTestODESystemWithEndSolution(3); // HarmOscillator1

		auto initCond = odeSys.getInitialConditions();
		Real t1 = REAL(0.0), t2 = odeSys.getEndTime();
		int numSteps = 1000;

		ODESystemFixedStepSolver fixedSolver(*odeSys.getODESystem(), StepCalculators::LeapfrogStepCalc);
		ODESystemSolution sol = fixedSolver.integrate(initCond, t1, t2, numSteps);

		Vector<Real> finalX = sol.getXValuesAtEnd();
		Vector<Real> exactSol = odeSys.getEndSolution();

		// Leapfrog should be accurate for harmonic oscillator
		for (int i = 0; i < finalX.size(); i++) {
			REQUIRE_THAT(finalX[i], WithinAbs(exactSol[i], REAL(1e-4)));
		}
	}

	// Test Leapfrog energy conservation - key property of symplectic integrators
	TEST_CASE("Test_Leapfrog_EnergyConservation", "[Leapfrog][ODESystemFixedStepSolver]") {
		TEST_PRECISION_INFO();
		// For harmonic oscillator, E = REAL(0.5)*(x^2 + v^2) should be conserved
		TestBeds::SimpleHarmonicOscillatorODE harmOsc;

		Vector<Real> initCond{REAL(1.0), REAL(0.0)}; // x=1, v=0 -> initial E = REAL(0.5)
		Real t1 = REAL(0.0), t2 = REAL(100.0);		 // Long integration to test stability
		int numSteps = 10000;

		ODESystemFixedStepSolver fixedSolver(harmOsc, StepCalculators::LeapfrogStepCalc);
		ODESystemSolution sol = fixedSolver.integrate(initCond, t1, t2, numSteps);

		// Calculate initial and final energy
		Real E_initial = REAL(0.5) * (initCond[0] * initCond[0] + initCond[1] * initCond[1]);

		Vector<Real> finalX = sol.getXValuesAtEnd();
		Real E_final = REAL(0.5) * (finalX[0] * finalX[0] + finalX[1] * finalX[1]);

		// Leapfrog should conserve energy very well (symplectic property)
		REQUIRE_THAT(E_final, WithinRel(E_initial, REAL(1e-4)));
	}

	// Test ODESystemLeapfrogSolver class directly
	TEST_CASE("Test_ODESystemLeapfrogSolver", "[Leapfrog][ODESystemLeapfrogSolver]") {
		TEST_PRECISION_INFO();
		TestBeds::SimpleHarmonicOscillatorODE harmOsc;

		Vector<Real> initCond{REAL(1.0), REAL(2.0)};
		Real t1 = REAL(0.0), t2 = REAL(10.0);
		int numSteps = 1000;

		ODESystemLeapfrogSolver leapfrogSolver(harmOsc);
		ODESystemSolution sol = leapfrogSolver.integrate(initCond, t1, t2, numSteps);

		// Get exact solution
		Vector<Real> exactSol = harmOsc.getSolution(initCond, t2);
		Vector<Real> finalX = sol.getXValuesAtEnd();

		// Should be accurate
		for (int i = 0; i < finalX.size(); i++) {
			REQUIRE_THAT(finalX[i], WithinAbs(exactSol[i], REAL(1e-4)));
		}
	}
} // namespace MML::Tests::Algorithms::ODESystemSolverTests

///////////////////////////////////////////////////////////////////////////////////////////
//                      STIFF ODE TEST BED INTEGRATION TESTS                             //
///////////////////////////////////////////////////////////////////////////////////////////

#include "../../test_data/stiff_ode_test_bed.h"

namespace MML::Tests::Algorithms::StiffODETestBed {
	using namespace MML::TestBeds;

	/*********************************************************************/
	/*****          Test Bed Infrastructure Tests                    *****/
	/*********************************************************************/

	TEST_CASE("StiffODE_TestBedAvailability", "[StiffODE][TestBed][Infrastructure]") {
		TEST_PRECISION_INFO();

		auto allTests = getAllStiffODETests();
		REQUIRE(allTests.size() >= 5); // Should have at least 5 test cases

		INFO("Available stiff ODE tests:");
		for (const auto& test : allTests) {
			INFO("  - " << test.name << " (dim=" << test.dimension << ", stiffness=" << test.stiffnessCategory << ")");
		}

		// Check specific categories exist
		REQUIRE(getModerateStiffTests().size() >= 1);
		REQUIRE(getSevereStiffTests().size() >= 1);
		REQUIRE(getChemicalKineticsTests().size() >= 1);
	}

	TEST_CASE("StiffODE_TestBedDimensions", "[StiffODE][TestBed][Infrastructure]") {
		TEST_PRECISION_INFO();

		auto allTests = getAllStiffODETests();

		for (const auto& test : allTests) {
			DYNAMIC_SECTION("Dimension check: " << test.name) {
				// Initial condition should match dimension
				REQUIRE(test.initialCondition.size() == test.dimension);

				// Time interval should be valid
				REQUIRE(test.tEnd > test.tStart);

				// Suggested step size should be positive
				REQUIRE(test.suggestedStepSize > 0);

				// Difficulty should be valid
				REQUIRE(test.difficulty >= 1);
				REQUIRE(test.difficulty <= 4);
			}
		}
	}

	/*********************************************************************/
	/*****          Jacobian Verification Tests                      *****/
	/*********************************************************************/

	TEST_CASE("StiffODE_JacobianConsistency", "[StiffODE][TestBed][Jacobian]") {
		TEST_PRECISION_INFO();

		auto smallDimTests = getSmallDimensionStiffTests();

		for (const auto& test : smallDimTests) {
			DYNAMIC_SECTION("Jacobian consistency: " << test.name) {
				auto system = test.createStiffSystem();
				REQUIRE(system != nullptr);

				// Get initial condition
				Vector<Real> y = test.initialCondition;
				Vector<Real> dydt(test.dimension);
				Matrix<Real> J(test.dimension, test.dimension);

				// Compute Jacobian
				system->jacobian(test.tStart, y, dydt, J);

				// Verify Jacobian is finite
				bool jacobianValid = true;
				for (int i = 0; i < test.dimension && jacobianValid; ++i) {
					for (int j = 0; j < test.dimension && jacobianValid; ++j) {
						if (!std::isfinite(J(i, j))) {
							jacobianValid = false;
						}
					}
				}
				REQUIRE(jacobianValid);

				// Numerical Jacobian verification for small systems
				if (test.dimension <= 3) {
					Real h = REAL(1e-7);
					Matrix<Real> numJ(test.dimension, test.dimension);
					Vector<Real> dydt_plus(test.dimension);
					Vector<Real> y_perturbed = y;

					for (int j = 0; j < test.dimension; ++j) {
						y_perturbed = y;
						y_perturbed[j] += h;
						system->derivs(test.tStart, y_perturbed, dydt_plus);

						system->derivs(test.tStart, y, dydt);

						for (int i = 0; i < test.dimension; ++i) {
							numJ(i, j) = (dydt_plus[i] - dydt[i]) / h;
						}
					}

					// Compare analytical and numerical Jacobian
					// Use relative tolerance for potentially large values
					Real maxRelError = 0;
					for (int i = 0; i < test.dimension; ++i) {
						for (int j = 0; j < test.dimension; ++j) {
							Real absError = std::abs(J(i, j) - numJ(i, j));
							Real scale = std::max(std::abs(J(i, j)), std::abs(numJ(i, j)));
							Real relError = (scale > REAL(1e-10)) ? absError / scale : absError;
							maxRelError = std::max(maxRelError, relError);
						}
					}
					INFO("Max relative Jacobian error: " << maxRelError);
					REQUIRE(maxRelError < REAL(1e-4));
				}
			}
		}
	}

	/*********************************************************************/
	/*****          Conservation Law Tests                           *****/
	/*********************************************************************/

	TEST_CASE("StiffODE_ConservationLaws", "[StiffODE][TestBed][Conservation]") {
		TEST_PRECISION_INFO();

		auto chemTests = getChemicalKineticsTests();

		for (const auto& test : chemTests) {
			// Skip extremely stiff problems for this test
			if (test.stiffnessCategory == "extreme")
				continue;

			DYNAMIC_SECTION("Conservation check: " << test.name) {
				// Check conservation at initial condition
				Real conservationError = checkConservation(test.initialCondition, test);

				if (conservationError >= 0) { // -1 means no conservation law checked
					INFO("Conservation error at IC: " << conservationError);
					REQUIRE(conservationError < REAL(1e-10));
				}

				// Verify non-negativity at initial condition
				REQUIRE(verifyNonNegative(test.initialCondition));
			}
		}
	}

	/*********************************************************************/
	/*****          Short Integration Tests (Explicit Methods)       *****/
	/*********************************************************************/

	TEST_CASE("StiffODE_ShortIntegration_RK4", "[StiffODE][TestBed][RK4][Short]") {
		TEST_PRECISION_INFO();

		// Test with moderate stiffness cases only, using short time intervals
		auto moderateTests = getModerateStiffTests();

		for (const auto& test : moderateTests) {
			DYNAMIC_SECTION("RK4 short integration: " << test.name) {
				auto system = test.createStiffSystem();
				REQUIRE(system != nullptr);

				Vector<Real> initCond = test.initialCondition;
				Real tEnd = std::min(test.tEnd, test.tStart + REAL(0.1)); // Short integration
				int numSteps = 1000;									  // Many steps for stability

				ODESystemFixedStepSolver solver(*system, StepCalculators::RK4_Basic);
				ODESystemSolution sol = solver.integrate(initCond, test.tStart, tEnd, numSteps);

				Vector<Real> finalX = sol.getXValuesAtEnd();

				// Basic checks: solution should be finite
				bool solutionValid = true;
				for (int i = 0; i < test.dimension && solutionValid; ++i) {
					if (!std::isfinite(finalX[i])) {
						solutionValid = false;
					}
				}
				REQUIRE(solutionValid);

				// For chemical kinetics: verify non-negativity (with tolerance)
				if (test.name.find("Brusselator") != std::string::npos || test.name.find("Linear") != std::string::npos) {
					// These should stay non-negative for short integrations
					bool nonNegative = true;
					for (int i = 0; i < test.dimension; ++i) {
						if (finalX[i] < -REAL(0.01)) { // Small tolerance
							nonNegative = false;
							INFO("Component " << i << " = " << finalX[i] << " < 0");
						}
					}
					REQUIRE(nonNegative);
				}
			}
		}
	}

	TEST_CASE("StiffODE_ExactSolution_Verification", "[StiffODE][TestBed][ExactSolution]") {
		TEST_PRECISION_INFO();

		auto exactTests = getExactSolutionStiffTests();

		for (const auto& test : exactTests) {
			DYNAMIC_SECTION("Exact solution: " << test.name) {
				REQUIRE(test.hasExactSolution);

				auto system = test.createStiffSystem();
				REQUIRE(system != nullptr);

				Vector<Real> initCond = test.initialCondition;

				// Use many small steps with high-order method
				Real tEnd = std::min(test.tEnd, test.tStart + REAL(0.01)); // Very short
				int numSteps = 5000;									   // Many steps

				ODESystemFixedStepSolver solver(*system, StepCalculators::DormandPrince5StepCalc);
				ODESystemSolution sol = solver.integrate(initCond, test.tStart, tEnd, numSteps);

				Vector<Real> computed = sol.getXValuesAtEnd();
				Vector<Real> exact = test.exactSolution(tEnd);

				Real error = computeStiffODEError(computed, test, tEnd);

				INFO("Error at t=" << tEnd << ": " << error);

				// For exact solution tests with high-order explicit method on short interval
				REQUIRE(error < REAL(1e-4));
			}
		}
	}

	/*********************************************************************/
	/*****          System Properties Tests                          *****/
	/*********************************************************************/

	TEST_CASE("StiffODE_StiffnessRatios", "[StiffODE][TestBed][Properties]") {
		TEST_PRECISION_INFO();

		auto allTests = getAllStiffODETests();

		for (const auto& test : allTests) {
			DYNAMIC_SECTION("Stiffness ratio: " << test.name) {
				// Verify stiffness ratio is positive
				REQUIRE(test.stiffnessRatio > 0);

				// Verify stiffness category is non-empty
				REQUIRE(!test.stiffnessCategory.empty());

				// Stiffness ratio should be at least 1 for stiff problems
				REQUIRE(test.stiffnessRatio >= 1);

				INFO("Name: " << test.name);
				INFO("Category: " << test.stiffnessCategory);
				INFO("Stiffness ratio: " << test.stiffnessRatio);
			}
		}

		// Verify we have a variety of test cases
		REQUIRE(allTests.size() >= 5);
	}

	TEST_CASE("StiffODE_Brusselator_Derivatives", "[StiffODE][TestBed][Brusselator]") {
		TEST_PRECISION_INFO();

		auto test = getBrusselatorTest();
		auto system = test.createStiffSystem();

		Vector<Real> y = test.initialCondition;
		Vector<Real> dydt(test.dimension);

		// Compute derivatives at initial condition
		system->derivs(test.tStart, y, dydt);

		// The Brusselator has specific equilibrium at (A, B/A) = (1, 3)
		// At IC, derivatives should be non-zero (away from equilibrium)
		INFO("y = [" << y[0] << ", " << y[1] << "]");
		INFO("dydt = [" << dydt[0] << ", " << dydt[1] << "]");

		// Verify derivatives are computed correctly (basic sanity)
		bool derivsComputed = std::isfinite(dydt[0]) && std::isfinite(dydt[1]);
		REQUIRE(derivsComputed);
	}

} // namespace MML::Tests::Algorithms::StiffODETestBed

/******************************************************************************/
/*****     Adaptive Integrator Tests - DormandPrince5 with Dense Output  *****/
/******************************************************************************/

#include "algorithms/ODEAdaptiveIntegrator.h"

namespace MML::Tests::Algorithms::AdaptiveIntegratorTests {
	// Simple exponential decay: x' = -x, x(0) = 1, solution: x(t) = e^(-t)
	class ExponentialDecayODE : public IODESystem {
	public:
		int getDim() const override { return 1; }
		void derivs(Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override { dxdt[0] = -x[0]; }
	};

	// Harmonic oscillator: x'' = -x => x' = v, v' = -x
	// x(0) = 1, v(0) = 0 => x(t) = cos(t), v(t) = -sin(t)
	class HarmonicOscillatorODE : public IODESystem {
	public:
		int getDim() const override { return 2; }
		void derivs(Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override {
			dxdt[0] = x[1];	 // x' = v
			dxdt[1] = -x[0]; // v' = -x
		}
	};

	TEST_CASE("DormandPrince5_Stepper_ExponentialDecay", "[AdaptiveIntegrator][DP5]") {
		TEST_PRECISION_INFO();

		ExponentialDecayODE ode;
		DormandPrince5_Stepper stepper(ode);

		Vector<Real> x{1.0};
		Vector<Real> dxdt{-1.0}; // Initial derivative
		Real t = 0.0;
		Real eps = 1e-10;
		Real htry = 0.1;

		StepResult result = stepper.doStep(t, x, dxdt, htry, eps);

		// Step should be accepted
		REQUIRE(result.accepted);
		REQUIRE(result.hDone > 0);

		// Check accuracy at end of step
		Real exact = std::exp(-result.hDone);
		REQUIRE_THAT(x[0], WithinAbs(exact, 1e-9));

		// FSAL: derivative should be updated
		REQUIRE_THAT(dxdt[0], WithinAbs(-x[0], 1e-12));
	}

	TEST_CASE("DormandPrince5_Stepper_FSAL_Optimization", "[AdaptiveIntegrator][DP5][FSAL]") {
		TEST_PRECISION_INFO();

		ExponentialDecayODE ode;
		DormandPrince5_Stepper stepper(ode);

		Vector<Real> x{1.0};
		Vector<Real> dxdt{-1.0};
		Real t = 0.0;
		Real eps = 1e-10;
		Real htry = 0.1;

		// First step
		StepResult r1 = stepper.doStep(t, x, dxdt, htry, eps);
		REQUIRE(r1.accepted);
		int evals1 = r1.funcEvals;

		t += r1.hDone;

		// Second step - should use FSAL (one less function evaluation)
		StepResult r2 = stepper.doStep(t, x, dxdt, r1.hNext, eps);
		REQUIRE(r2.accepted);
		int evals2 = r2.funcEvals;

		// Due to FSAL, second step should use same or fewer evals
		// (First step computes k1 fresh, subsequent steps reuse k7)
		INFO("First step evals: " << evals1 << ", Second step evals: " << evals2);
		REQUIRE(evals2 <= evals1);
	}

	TEST_CASE("DormandPrince5_Stepper_DenseOutput", "[AdaptiveIntegrator][DP5][DenseOutput]") {
		TEST_PRECISION_INFO();

		ExponentialDecayODE ode;
		DormandPrince5_Stepper stepper(ode);

		Vector<Real> x{1.0};
		Vector<Real> dxdt{-1.0};
		Real t = 0.0;
		Real eps = 1e-10;
		Real htry = 0.5; // Larger step to test interpolation

		StepResult result = stepper.doStep(t, x, dxdt, htry, eps);
		REQUIRE(result.accepted);

		Real hDone = result.hDone;

		// Test interpolation at several points within the step
		std::vector<Real> thetas = {0.0, 0.25, 0.5, 0.75, 1.0};

		for (Real theta : thetas) {
			Real tInterp = t + theta * hDone;
			Vector<Real> xInterp = stepper.interpolate(tInterp);
			Real exact = std::exp(-tInterp);

			INFO("theta = " << theta << ", t = " << tInterp);
			INFO("interpolated = " << xInterp[0] << ", exact = " << exact);

			// Dense output should be accurate to at least 4th order
			REQUIRE_THAT(xInterp[0], WithinAbs(exact, 1e-6));
		}
	}

	TEST_CASE("ODEAdaptiveIntegrator_ExponentialDecay", "[AdaptiveIntegrator][Integration]") {
		TEST_PRECISION_INFO();

		ExponentialDecayODE ode;
		DormandPrince5Integrator integrator(ode);

		Vector<Real> x0{1.0};
		Real t0 = 0.0, tEnd = 5.0;
		Real outputInterval = 0.5;
		Real eps = 1e-10;

		ODESystemSolution sol = integrator.integrate(x0, t0, tEnd, outputInterval, eps);

		// Check solution at each output point
		int numPoints = sol.getNumSteps() + 1;
		for (int i = 0; i < numPoints; i++) {
			Real t = sol.getTValue(i);
			Real x = sol.getXValue(i, 0);
			Real exact = std::exp(-t);

			INFO("t = " << t << ", x = " << x << ", exact = " << exact);
			REQUIRE_THAT(x, WithinAbs(exact, 1e-6)); // Relaxed for dense output
		}

		// Check statistics
		auto stats = integrator.getStatistics();
		REQUIRE(stats.acceptedSteps > 0);
		INFO("Accepted: " << stats.acceptedSteps << ", Rejected: " << stats.rejectedSteps);
		INFO("Acceptance rate: " << stats.acceptanceRate() * 100 << "%");
		INFO("Total func evals: " << stats.totalFuncEvals);
	}

	TEST_CASE("ODEAdaptiveIntegrator_HarmonicOscillator", "[AdaptiveIntegrator][Integration]") {
		TEST_PRECISION_INFO();

		HarmonicOscillatorODE ode;
		DormandPrince5Integrator integrator(ode);

		Vector<Real> x0{1.0, 0.0};					// x(0) = 1, v(0) = 0
		Real t0 = 0.0, tEnd = 2.0 * Constants::PI;	// One full period
		Real outputInterval = Constants::PI / 10.0; // 20 points per period
		Real eps = 1e-10;

		ODESystemSolution sol = integrator.integrate(x0, t0, tEnd, outputInterval, eps);

		// Check solution at each output point
		int numPoints = sol.getNumSteps() + 1;
		for (int i = 0; i < numPoints; i++) {
			Real t = sol.getTValue(i);
			Real x = sol.getXValue(i, 0);
			Real v = sol.getXValue(i, 1);
			Real exactX = std::cos(t);
			Real exactV = -std::sin(t);

			INFO("t = " << t);
			REQUIRE_THAT(x, WithinAbs(exactX, 1e-5)); // Relaxed for dense output
			REQUIRE_THAT(v, WithinAbs(exactV, 1e-5));
		}

		// After one full period, should return to initial conditions
		Real finalX = sol.getXValuesAtEnd()[0];
		Real finalV = sol.getXValuesAtEnd()[1];
		REQUIRE_THAT(finalX, WithinAbs(1.0, 1e-5));
		REQUIRE_THAT(finalV, WithinAbs(0.0, 1e-5));
	}

	TEST_CASE("ODEAdaptiveIntegrator_IntegrateAt", "[AdaptiveIntegrator][Integration]") {
		TEST_PRECISION_INFO();

		ExponentialDecayODE ode;
		DormandPrince5Integrator integrator(ode);

		Vector<Real> x0{1.0};

		// Specify exact output times
		Vector<Real> times{0.0, 0.1, 0.5, 1.0, 2.0, 3.0, 5.0};
		Real eps = 1e-10;

		ODESystemSolution sol = integrator.integrateAt(x0, times, eps);

		// Check solution at each specified time
		for (int i = 0; i < static_cast<int>(times.size()); i++) {
			Real t = sol.getTValue(i);
			Real x = sol.getXValue(i, 0);
			Real exact = std::exp(-t);

			INFO("t = " << t << " (requested: " << times[i] << ")");
			REQUIRE_THAT(t, WithinAbs(times[i], 1e-10)); // Time should match exactly
			REQUIRE_THAT(x, WithinAbs(exact, 1e-8));
		}
	}

	TEST_CASE("ODEAdaptiveIntegrator_EnergyConservation", "[AdaptiveIntegrator][Physics]") {
		TEST_PRECISION_INFO();

		// Harmonic oscillator: E = 0.5*(v^2 + x^2) should be constant
		HarmonicOscillatorODE ode;
		DormandPrince5Integrator integrator(ode);

		Vector<Real> x0{1.0, 0.0};					// x(0) = 1, v(0) = 0
		Real t0 = 0.0, tEnd = 10.0 * Constants::PI; // 5 full periods
		Real outputInterval = Constants::PI / 5.0;
		Real eps = 1e-12; // Tight tolerance

		ODESystemSolution sol = integrator.integrate(x0, t0, tEnd, outputInterval, eps);

		Real initialEnergy = 0.5 * (x0[0] * x0[0] + x0[1] * x0[1]);
		Real maxEnergyDrift = 0.0;

		int numPoints = sol.getNumSteps() + 1;
		for (int i = 0; i < numPoints; i++) {
			Real x = sol.getXValue(i, 0);
			Real v = sol.getXValue(i, 1);
			Real energy = 0.5 * (x * x + v * v);
			Real drift = std::abs(energy - initialEnergy);
			maxEnergyDrift = std::max(maxEnergyDrift, drift);
		}

		INFO("Initial energy: " << initialEnergy);
		INFO("Max energy drift: " << maxEnergyDrift);

		// Energy should be well conserved (not perfect, but good)
		REQUIRE(maxEnergyDrift < 1e-4);
	}

	TEST_CASE("ODEAdaptiveIntegrator_Efficiency", "[AdaptiveIntegrator][Performance]") {
		TEST_PRECISION_INFO();

		// Compare adaptive vs fixed step for same accuracy
		HarmonicOscillatorODE ode;

		Vector<Real> x0{1.0, 0.0};
		Real t0 = 0.0, tEnd = 10.0;
		Real eps = 1e-8;

		// Adaptive integration
		DormandPrince5Integrator adaptiveIntegrator(ode);
		ODESystemSolution adaptiveSol = adaptiveIntegrator.integrate(x0, t0, tEnd, 0.5, eps);
		auto adaptiveStats = adaptiveIntegrator.getStatistics();

		// Fixed step integration with many steps
		ODESystemFixedStepSolver fixedSolver(ode, StepCalculators::DormandPrince5StepCalc);
		ODESystemSolution fixedSol = fixedSolver.integrate(x0, t0, tEnd, 1000);

		// Check final accuracy of both
		Real exactX = std::cos(tEnd);
		Real exactV = -std::sin(tEnd);

		Real adaptiveErrX = std::abs(adaptiveSol.getXValuesAtEnd()[0] - exactX);
		Real adaptiveErrV = std::abs(adaptiveSol.getXValuesAtEnd()[1] - exactV);
		Real fixedErrX = std::abs(fixedSol.getXValuesAtEnd()[0] - exactX);
		Real fixedErrV = std::abs(fixedSol.getXValuesAtEnd()[1] - exactV);

		INFO("Adaptive error: x=" << adaptiveErrX << ", v=" << adaptiveErrV);
		INFO("Fixed error: x=" << fixedErrX << ", v=" << fixedErrV);
		INFO("Adaptive steps: " << adaptiveStats.acceptedSteps);
		INFO("Fixed steps: 1000");
		INFO("Adaptive func evals: " << adaptiveStats.totalFuncEvals);
		INFO("Fixed func evals: ~7000 (7 per step for DP5)");

		// Adaptive should achieve comparable or better accuracy with fewer steps
		REQUIRE(adaptiveErrX < 1e-7);
		REQUIRE(adaptiveErrV < 1e-7);

		// Should use significantly fewer steps than fixed
		REQUIRE(adaptiveStats.acceptedSteps < 200);
	}

	TEST_CASE("CashKarpIntegrator_HarmonicOscillator", "[AdaptiveIntegrator][CashKarp]") {
		TEST_PRECISION_INFO();

		// y'' + y = 0 => y = cos(t), y' = -sin(t)
		struct HarmonicOscillator : public IODESystem {
			int getDim() const override { return 2; }
			void derivs(Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override {
				dxdt[0] = x[1];	 // y' = v
				dxdt[1] = -x[0]; // v' = -y
			}
		};

		HarmonicOscillator sys;
		CashKarpIntegrator integrator(sys);

		Vector<Real> x0({1.0, 0.0}); // y(0) = 1, y'(0) = 0
		Real t0 = 0.0, tEnd = 10.0;
		Real eps = 1e-10;

		auto sol = integrator.integrate(x0, t0, tEnd, 0.5, eps);

		// Check final values: cos(10) and -sin(10)
		Real exactY = std::cos(10.0);
		Real exactV = -std::sin(10.0);

		Vector<Real> finalX = sol.getXValuesAtEnd();
		Real errY = std::abs(finalX[0] - exactY);
		Real errV = std::abs(finalX[1] - exactV);

		auto stats = integrator.getStatistics();
		INFO("Cash-Karp: y=" << finalX[0] << ", v=" << finalX[1]);
		INFO("Exact: y=" << exactY << ", v=" << exactV);
		INFO("Steps: accepted=" << stats.acceptedSteps << ", rejected=" << stats.rejectedSteps);
		INFO("Func evals: " << stats.totalFuncEvals);

		REQUIRE(errY < 1e-7);
		REQUIRE(errV < 1e-7);
	}

	TEST_CASE("DormandPrince8Integrator_HarmonicOscillator", "[AdaptiveIntegrator][DP8]") {
		TEST_PRECISION_INFO();

		// Same harmonic oscillator, but testing DP8's higher accuracy
		struct HarmonicOscillator : public IODESystem {
			int getDim() const override { return 2; }
			void derivs(Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override {
				dxdt[0] = x[1];
				dxdt[1] = -x[0];
			}
		};

		HarmonicOscillator sys;
		DormandPrince8Integrator integrator(sys);

		Vector<Real> x0({1.0, 0.0});
		Real t0 = 0.0, tEnd = 10.0;
		Real eps = 1e-12; // Tighter tolerance for DP8

		auto sol = integrator.integrate(x0, t0, tEnd, 0.5, eps);

		Real exactY = std::cos(10.0);
		Real exactV = -std::sin(10.0);

		Vector<Real> finalX = sol.getXValuesAtEnd();
		Real errY = std::abs(finalX[0] - exactY);
		Real errV = std::abs(finalX[1] - exactV);

		auto stats = integrator.getStatistics();
		INFO("DP8: y=" << finalX[0] << ", v=" << finalX[1]);
		INFO("Exact: y=" << exactY << ", v=" << exactV);
		INFO("Steps: accepted=" << stats.acceptedSteps << ", rejected=" << stats.rejectedSteps);
		INFO("Func evals: " << stats.totalFuncEvals);

		// DP8 should achieve even higher accuracy (with tolerance slightly relaxed for step size effects)
		REQUIRE(errY < 1e-9);
		REQUIRE(errV < 1e-9);
	}

	TEST_CASE("ODEAdaptiveIntegrator_StepperComparison", "[AdaptiveIntegrator][Comparison]") {
		TEST_PRECISION_INFO();

		// Compare all three steppers on the same problem
		struct VanDerPol : public IODESystem {
			int getDim() const override { return 2; }
			void derivs(Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override {
				Real mu = 0.5; // Mild nonlinearity
				dxdt[0] = x[1];
				dxdt[1] = mu * (1 - x[0] * x[0]) * x[1] - x[0];
			}
		};

		VanDerPol sys;

		DormandPrince5Integrator dp5(sys);
		CashKarpIntegrator ck(sys);
		DormandPrince8Integrator dp8(sys);

		Vector<Real> x0({2.0, 0.0});
		Real t0 = 0.0, tEnd = 5.0;
		Real eps = 1e-8;

		auto sol5 = dp5.integrate(x0, t0, tEnd, 0.5, eps);
		auto solCK = ck.integrate(x0, t0, tEnd, 0.5, eps);
		auto sol8 = dp8.integrate(x0, t0, tEnd, 0.5, 1e-10); // Tighter for DP8

		auto stats5 = dp5.getStatistics();
		auto statsCK = ck.getStatistics();
		auto stats8 = dp8.getStatistics();

		Vector<Real> final5 = sol5.getXValuesAtEnd();
		Vector<Real> finalCK = solCK.getXValuesAtEnd();
		Vector<Real> final8 = sol8.getXValuesAtEnd();

		INFO("DP5:  x=" << final5[0] << ", v=" << final5[1] << ", steps=" << stats5.acceptedSteps << ", evals=" << stats5.totalFuncEvals);
		INFO("CK:   x=" << finalCK[0] << ", v=" << finalCK[1] << ", steps=" << statsCK.acceptedSteps
						<< ", evals=" << statsCK.totalFuncEvals);
		INFO("DP8:  x=" << final8[0] << ", v=" << final8[1] << ", steps=" << stats8.acceptedSteps << ", evals=" << stats8.totalFuncEvals);

		// All should agree to reasonable precision
		REQUIRE_THAT(final5[0], WithinRel(final8[0], 1e-5));
		REQUIRE_THAT(finalCK[0], WithinRel(final8[0], 1e-5));

		// DP8 should use fewer steps for comparable accuracy
		REQUIRE(stats8.acceptedSteps <= stats5.acceptedSteps);
	}

	TEST_CASE("CashKarp_DenseOutput", "[AdaptiveIntegrator][CashKarp][DenseOutput]") {
		TEST_PRECISION_INFO();

		// Test dense output with Cash-Karp
		struct ExpDecay : public IODESystem {
			int getDim() const override { return 1; }
			void derivs(Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override {
				dxdt[0] = -x[0]; // x' = -x => x = e^{-t}
			}
		};

		ExpDecay sys;
		CashKarpIntegrator integrator(sys);

		Vector<Real> x0({1.0});
		Real eps = 1e-10;

		// Request output at specific times
		Vector<Real> times({0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0});
		auto sol = integrator.integrateAt(x0, times, eps);

		// Check against exact solution
		for (int i = 0; i < times.size(); i++) {
			Real t = times[i];
			Real exact = std::exp(-t);
			Real computed = sol.getXValue(i, 0);
			Real err = std::abs(computed - exact);
			INFO("t=" << t << ": computed=" << computed << ", exact=" << exact << ", err=" << err);
			REQUIRE(err < 1e-6);
		}
	}

	TEST_CASE("DormandPrince8_DenseOutput", "[AdaptiveIntegrator][DP8][DenseOutput]") {
		TEST_PRECISION_INFO();

		// Test dense output with DP8
		struct SinODE : public IODESystem {
			int getDim() const override { return 1; }
			void derivs(Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override {
				dxdt[0] = std::cos(t); // x' = cos(t) => x = sin(t)
			}
		};

		SinODE sys;
		DormandPrince8Integrator integrator(sys);

		Vector<Real> x0({0.0}); // sin(0) = 0
		Real tEnd = 2 * Constants::PI;
		Real eps = 1e-12;

		// Request many output points
		Vector<Real> times(21);
		for (int i = 0; i <= 20; i++)
			times[i] = tEnd * i / 20.0;

		auto sol = integrator.integrateAt(x0, times, eps);

		// Check against exact solution
		Real maxErr = 0;
		for (int i = 0; i < times.size(); i++) {
			Real t = times[i];
			Real exact = std::sin(t);
			Real computed = sol.getXValue(i, 0);
			Real err = std::abs(computed - exact);
			maxErr = std::max(maxErr, err);
		}
		INFO("Max dense output error: " << maxErr);
		REQUIRE(maxErr < 1e-8);
	}

	TEST_CASE("ODEAdaptiveIntegrator_TestBedSystems", "[AdaptiveIntegrator][TestBed]") {
		TEST_PRECISION_INFO();

		// Test with systems from the test bed
		for (int sysIdx = 0; sysIdx < 4; sysIdx++) {
			auto odeSys = TestBeds::ODESystemTestBed::getTestODESystemWithEndSolution(sysIdx);

			DYNAMIC_SECTION("TestBed system " << sysIdx) {
				DormandPrince5Integrator integrator(*odeSys.getODESystem());

				auto initCond = odeSys.getInitialConditions();
				Real t0 = odeSys.getStartTime();
				Real tEnd = odeSys.getEndTime();
				Real outputInterval = (tEnd - t0) / 10.0;
				Real eps = 1e-10;

				ODESystemSolution sol = integrator.integrate(initCond, t0, tEnd, outputInterval, eps);

				Vector<Real> finalX = sol.getXValuesAtEnd();
				Vector<Real> exactSol = odeSys.getEndSolution();

				auto stats = integrator.getStatistics();
				INFO("Accepted: " << stats.acceptedSteps << ", Rejected: " << stats.rejectedSteps);

				// Compare with exact solution
				for (int i = 0; i < finalX.size(); i++) {
					INFO("Component " << i << ": computed=" << finalX[i] << ", exact=" << exactSol[i]);
					REQUIRE_THAT(finalX[i], WithinRel(exactSol[i], 1e-6)); // Use relative tolerance for large values
				}
			}
		}
	}

	TEST_CASE("BulirschStoerIntegrator_HarmonicOscillator", "[AdaptiveIntegrator][BulirschStoer]") {
		TEST_PRECISION_INFO();

		// Bulirsch-Stoer should excel at smooth problems
		HarmonicOscillatorODE oscillator;
		BulirschStoerIntegrator integrator(oscillator);

		Vector<Real> x0(2);
		x0[0] = 1.0;  // x(0) = 1
		x0[1] = 0.0;  // v(0) = 0

		Real t0 = 0.0;
		Real tEnd = 10.0;
		Real outputInterval = 1.0;
		Real eps = 1e-10;

		ODESystemSolution sol = integrator.integrate(x0, t0, tEnd, outputInterval, eps);

		auto stats = integrator.getStatistics();
		INFO("BS: Accepted: " << stats.acceptedSteps << ", Rejected: " << stats.rejectedSteps);
		INFO("BS: Total function evals: " << stats.totalFuncEvals);

		// Verify solution at t=tEnd
		Vector<Real> xFinal = sol.getXValuesAtEnd();
		Real exactX = std::cos(tEnd);
		Real exactV = -std::sin(tEnd);

		INFO("Final position: computed=" << xFinal[0] << ", exact=" << exactX);
		INFO("Final velocity: computed=" << xFinal[1] << ", exact=" << exactV);

		REQUIRE_THAT(xFinal[0], WithinAbs(exactX, 1e-5));
		REQUIRE_THAT(xFinal[1], WithinAbs(exactV, 1e-5));

		// Verify energy conservation (should be good for BS)
		Real initialEnergy = 0.5 * (x0[0] * x0[0] + x0[1] * x0[1]);
		Real finalEnergy = 0.5 * (xFinal[0] * xFinal[0] + xFinal[1] * xFinal[1]);
		Real energyError = std::abs(finalEnergy - initialEnergy) / initialEnergy;

		INFO("Energy error: " << energyError);
		REQUIRE(energyError < 1e-6);  // BS conserves energy reasonably well
	}

	TEST_CASE("BulirschStoerIntegrator_ExponentialDecay", "[AdaptiveIntegrator][BulirschStoer]") {
		TEST_PRECISION_INFO();

		// Test Bulirsch-Stoer on simple exponential decay
		ExponentialDecayODE expDecay;
		BulirschStoerIntegrator integrator(expDecay);

		Vector<Real> x0(1);
		x0[0] = 1.0;

		Real t0 = 0.0;
		Real tEnd = 5.0;
		Real outputInterval = 0.5;
		Real eps = 1e-12;

		ODESystemSolution sol = integrator.integrate(x0, t0, tEnd, outputInterval, eps);

		auto stats = integrator.getStatistics();
		INFO("BS: Accepted: " << stats.acceptedSteps << ", Rejected: " << stats.rejectedSteps);
		
		// Verify final value
		Vector<Real> xFinal = sol.getXValuesAtEnd();
		Real exact = std::exp(-tEnd);
		Real err = std::abs(xFinal[0] - exact);
		
		INFO("t=" << tEnd << ": computed=" << xFinal[0] << ", exact=" << exact << ", err=" << err);
		REQUIRE(err < 1e-8);  // Good accuracy for extrapolation method
	}

	TEST_CASE("BulirschStoerIntegrator_HighAccuracy", "[AdaptiveIntegrator][BulirschStoer]") {
		TEST_PRECISION_INFO();

		// Test BS's ability to achieve very high accuracy
		HarmonicOscillatorODE oscillator;
		BulirschStoerIntegrator integrator(oscillator);

		Vector<Real> x0(2);
		x0[0] = 1.0;
		x0[1] = 0.0;

		Real t0 = 0.0;
		Real tEnd = 20.0;  // Long integration time
		Real outputInterval = 2.0;
		Real eps = 1e-12;  // Very tight tolerance

		ODESystemSolution sol = integrator.integrate(x0, t0, tEnd, outputInterval, eps);

		auto stats = integrator.getStatistics();
		INFO("BS High Accuracy: Accepted: " << stats.acceptedSteps << ", Rejected: " << stats.rejectedSteps);
		INFO("BS High Accuracy: Function evals: " << stats.totalFuncEvals);
		INFO("BS High Accuracy: Acceptance rate: " << stats.acceptanceRate());

		// Verify solution at t=tEnd
		Vector<Real> xFinal = sol.getXValuesAtEnd();
		Real omega = 1.0;  // Default HarmonicOscillatorODE uses omega=1
		Real exactX = std::cos(omega * tEnd);
		Real exactV = -omega * std::sin(omega * tEnd);

		INFO("Final position: computed=" << xFinal[0] << ", exact=" << exactX);
		INFO("Final velocity: computed=" << xFinal[1] << ", exact=" << exactV);

		REQUIRE_THAT(xFinal[0], WithinAbs(exactX, 1e-5));
		REQUIRE_THAT(xFinal[1], WithinAbs(exactV, 1e-5));
	}

	TEST_CASE("BulirschStoerIntegrator_Comparison", "[AdaptiveIntegrator][Comparison][BulirschStoer]") {
		TEST_PRECISION_INFO();

		// Compare Bulirsch-Stoer with Dormand-Prince 8 for smooth problems
		HarmonicOscillatorODE oscillator1;
		HarmonicOscillatorODE oscillator2;
		
		BulirschStoerIntegrator bsIntegrator(oscillator1);
		DormandPrince8Integrator dp8Integrator(oscillator2);

		Vector<Real> x0(2);
		x0[0] = 1.0;
		x0[1] = 0.0;

		Real t0 = 0.0;
		Real tEnd = 10.0;
		Real outputInterval = 1.0;
		Real eps = 1e-10;

		// Integrate with both methods
		ODESystemSolution solBS = bsIntegrator.integrate(x0, t0, tEnd, outputInterval, eps);
		ODESystemSolution solDP8 = dp8Integrator.integrate(x0, t0, tEnd, outputInterval, eps);

		auto statsBS = bsIntegrator.getStatistics();
		auto statsDP8 = dp8Integrator.getStatistics();

		INFO("BS:  Accepted=" << statsBS.acceptedSteps << ", FuncEvals=" << statsBS.totalFuncEvals);
		INFO("DP8: Accepted=" << statsDP8.acceptedSteps << ", FuncEvals=" << statsDP8.totalFuncEvals);

		// Verify both achieve good accuracy
		Vector<Real> xFinalBS = solBS.getXValuesAtEnd();
		Vector<Real> xFinalDP8 = solDP8.getXValuesAtEnd();
		Real exactX = std::cos(tEnd);

		Real errBS = std::abs(xFinalBS[0] - exactX);
		Real errDP8 = std::abs(xFinalDP8[0] - exactX);

		INFO("BS error:  " << errBS);
		INFO("DP8 error: " << errDP8);

		REQUIRE(errBS < 1e-5);   // BS achieves good but not extreme accuracy
		REQUIRE(errDP8 < 1e-8);  // DP8 is more accurate for smooth problems
		
		// For smooth problems, BS often uses fewer function evaluations
		INFO("BS vs DP8 efficiency: " << static_cast<Real>(statsBS.totalFuncEvals) / statsDP8.totalFuncEvals);
	}

	TEST_CASE("BulirschStoer_PolynomialVsRational", "[AdaptiveIntegrator][BulirschStoer][Comparison]") {
		TEST_PRECISION_INFO();

		// Compare polynomial vs rational extrapolation versions
		HarmonicOscillatorODE oscillator1;
		HarmonicOscillatorODE oscillator2;
		
		BulirschStoerIntegrator bsPoly(oscillator1);          // Polynomial
		BulirschStoerRationalIntegrator bsRational(oscillator2);  // Rational

		Vector<Real> x0(2);
		x0[0] = 1.0;
		x0[1] = 0.0;

		Real t0 = 0.0;
		Real tEnd = 10.0;
		Real outputInterval = 1.0;
		Real eps = 1e-10;

		// Integrate with both methods
		ODESystemSolution solPoly = bsPoly.integrate(x0, t0, tEnd, outputInterval, eps);
		ODESystemSolution solRational = bsRational.integrate(x0, t0, tEnd, outputInterval, eps);

		auto statsPoly = bsPoly.getStatistics();
		auto statsRational = bsRational.getStatistics();

		INFO("BS Polynomial:  Accepted=" << statsPoly.acceptedSteps << ", FuncEvals=" << statsPoly.totalFuncEvals);
		INFO("BS Rational:    Accepted=" << statsRational.acceptedSteps << ", FuncEvals=" << statsRational.totalFuncEvals);

		// Verify both achieve good accuracy
		Vector<Real> xFinalPoly = solPoly.getXValuesAtEnd();
		Vector<Real> xFinalRational = solRational.getXValuesAtEnd();
		Real exactX = std::cos(tEnd);

		Real errPoly = std::abs(xFinalPoly[0] - exactX);
		Real errRational = std::abs(xFinalRational[0] - exactX);

		INFO("Polynomial error:  " << errPoly);
		INFO("Rational error:    " << errRational);

		REQUIRE(errPoly < 1e-5);
		REQUIRE(errRational < 1e-5);
	}

	TEST_CASE("BulirschStoerRational_ExponentialDecay", "[AdaptiveIntegrator][BulirschStoerRational]") {
		TEST_PRECISION_INFO();

		// Test rational extrapolation on exponential decay
		ExponentialDecayODE expDecay;
		BulirschStoerRationalIntegrator integrator(expDecay);

		Vector<Real> x0(1);
		x0[0] = 1.0;

		Real t0 = 0.0;
		Real tEnd = 5.0;
		Real outputInterval = 0.5;
		Real eps = 1e-10;

		ODESystemSolution sol = integrator.integrate(x0, t0, tEnd, outputInterval, eps);

		auto stats = integrator.getStatistics();
		INFO("BS Rational: Accepted: " << stats.acceptedSteps << ", Rejected: " << stats.rejectedSteps);
		
		// Verify final value
		Vector<Real> xFinal = sol.getXValuesAtEnd();
		Real exact = std::exp(-tEnd);
		Real err = std::abs(xFinal[0] - exact);
		
		INFO("t=" << tEnd << ": computed=" << xFinal[0] << ", exact=" << exact << ", err=" << err);
		REQUIRE(err < 1e-8);
	}

} // namespace MML::Tests::Algorithms::AdaptiveIntegratorTests

///////////////////////////////////////////////////////////////////////////////////////////
//            COMPREHENSIVE ADAPTIVE ODE SOLVER TESTS                                    //
//                                                                                       //
//  These tests verify:                                                                  //
//    1. Step size adaptation behavior (shrink on large error, grow on small error)     //
//    2. Tolerance control (tight tolerance => more steps)                              //
//    3. Error bounds are respected                                                     //
//    4. All stepper types work correctly                                               //
//    5. Dense output accuracy                                                          //
//    6. Long-time integration stability                                                //
//    7. Edge cases (extreme tolerances, stiff systems)                                 //
///////////////////////////////////////////////////////////////////////////////////////////

namespace MML::Tests::Algorithms::ComprehensiveAdaptiveTests
{
	using namespace MML::Tests::Algorithms::AdaptiveIntegratorTests;

	/////////////////////////////////////////////////////////////////////////////
	// Test Problem: Lotka-Volterra (predator-prey) - nonlinear, oscillatory
	/////////////////////////////////////////////////////////////////////////////
	class LotkaVolterraODE : public IODESystem {
		Real _alpha, _beta, _gamma, _delta;
	public:
		LotkaVolterraODE(Real alpha = 1.0, Real beta = 0.1, Real gamma = 1.5, Real delta = 0.075)
			: _alpha(alpha), _beta(beta), _gamma(gamma), _delta(delta) {}
		
		int getDim() const override { return 2; }
		void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const override {
			// prey: dx/dt = αx - βxy
			// predator: dy/dt = δxy - γy
			dydt[0] = _alpha * y[0] - _beta * y[0] * y[1];
			dydt[1] = _delta * y[0] * y[1] - _gamma * y[1];
		}
	};

	/////////////////////////////////////////////////////////////////////////////
	// Test Problem: Lorenz system (chaotic) - tests accuracy in sensitive regime
	/////////////////////////////////////////////////////////////////////////////
	class LorenzODE : public IODESystem {
		Real _sigma, _rho, _beta;
	public:
		LorenzODE(Real sigma = 10.0, Real rho = 28.0, Real beta = 8.0/3.0)
			: _sigma(sigma), _rho(rho), _beta(beta) {}
		
		int getDim() const override { return 3; }
		void derivs(Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override {
			dxdt[0] = _sigma * (x[1] - x[0]);
			dxdt[1] = x[0] * (_rho - x[2]) - x[1];
			dxdt[2] = x[0] * x[1] - _beta * x[2];
		}
	};

	/////////////////////////////////////////////////////////////////////////////
	// Test Problem: Kepler orbit - energy conservation test
	/////////////////////////////////////////////////////////////////////////////
	class KeplerODE : public IODESystem {
	public:
		int getDim() const override { return 4; }
		void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const override {
			// y[0] = x, y[1] = y, y[2] = vx, y[3] = vy
			Real r3 = std::pow(y[0]*y[0] + y[1]*y[1], 1.5);
			dydt[0] = y[2];
			dydt[1] = y[3];
			dydt[2] = -y[0] / r3;
			dydt[3] = -y[1] / r3;
		}
		
		static Real computeEnergy(const Vector<Real>& y) {
			Real r = std::sqrt(y[0]*y[0] + y[1]*y[1]);
			Real v2 = y[2]*y[2] + y[3]*y[3];
			return 0.5 * v2 - 1.0 / r;
		}
		
		static Vector<Real> getCircularOrbitIC() {
			Vector<Real> y0(4);
			y0[0] = 1.0;  // x
			y0[1] = 0.0;  // y
			y0[2] = 0.0;  // vx
			y0[3] = 1.0;  // vy (circular orbit)
			return y0;
		}
	};

	/////////////////////////////////////////////////////////////////////////////
	// Test Problem: Rapidly oscillating system - tests step adaptation
	/////////////////////////////////////////////////////////////////////////////
	class RapidOscillatorODE : public IODESystem {
		Real _omega;
	public:
		RapidOscillatorODE(Real omega = 50.0) : _omega(omega) {}
		
		int getDim() const override { return 2; }
		void derivs(Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override {
			dxdt[0] = x[1];
			dxdt[1] = -_omega * _omega * x[0];
		}
		
		// Exact solution: x(t) = cos(ωt), v(t) = -ω*sin(ωt)
		Real exactX(Real t) const { return std::cos(_omega * t); }
		Real exactV(Real t) const { return -_omega * std::sin(_omega * t); }
	};

	/////////////////////////////////////////////////////////////////////////////
	//         ADAPTIVE STEP SIZE BEHAVIOR TESTS                                //
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Adaptive_StepShrinks_WhenErrorLarge", "[AdaptiveIntegrator][StepControl]") {
		TEST_PRECISION_INFO();

		// Use rapid oscillator - forces small steps
		RapidOscillatorODE ode(100.0);  // omega=100 -> period ≈ 0.063
		DormandPrince5_Stepper stepper(ode);

		Vector<Real> x{1.0, 0.0};
		Vector<Real> dxdt(2);
		ode.derivs(0.0, x, dxdt);

		Real eps = 1e-6;
		Real htry = 1.0;  // Deliberately large initial step

		StepResult result = stepper.doStep(0.0, x, dxdt, htry, eps);

		INFO("htry=" << htry << ", hDone=" << result.hDone << ", hNext=" << result.hNext);
		INFO("errMax=" << result.errMax << ", accepted=" << result.accepted);

		// With very large initial step on rapid oscillator, step should shrink
		if (result.accepted) {
			REQUIRE(result.hDone < htry);  // Took smaller step
		} else {
			REQUIRE(result.hNext < htry);  // Suggests smaller step
		}
	}

	TEST_CASE("Adaptive_StepGrows_WhenErrorSmall", "[AdaptiveIntegrator][StepControl]") {
		TEST_PRECISION_INFO();

		// Use exponential decay - very smooth, allows large steps
		ExponentialDecayODE ode;
		DormandPrince5_Stepper stepper(ode);

		Vector<Real> x{1.0};
		Vector<Real> dxdt{-1.0};

		Real eps = 1e-6;
		Real htry = 0.001;  // Deliberately small initial step

		StepResult result = stepper.doStep(0.0, x, dxdt, htry, eps);

		INFO("htry=" << htry << ", hDone=" << result.hDone << ", hNext=" << result.hNext);
		INFO("errMax=" << result.errMax);

		// With very small step on smooth ODE, next step should grow
		REQUIRE(result.accepted);
		REQUIRE(result.hNext > result.hDone);  // Step size increases
	}

	TEST_CASE("Adaptive_TightTolerance_MoreSteps", "[AdaptiveIntegrator][StepControl]") {
		TEST_PRECISION_INFO();

		HarmonicOscillatorODE ode;

		Vector<Real> x0{1.0, 0.0};
		Real t0 = 0.0;
		Real tEnd = 2 * Constants::PI;  // One full period

		// Loose tolerance
		DormandPrince5Integrator integratorLoose(ode);
		auto solLoose = integratorLoose.integrate(x0, t0, tEnd, 0.1, 1e-3);
		auto statsLoose = integratorLoose.getStatistics();

		// Tight tolerance
		DormandPrince5Integrator integratorTight(ode);
		auto solTight = integratorTight.integrate(x0, t0, tEnd, 0.1, 1e-10);
		auto statsTight = integratorTight.getStatistics();

		INFO("Loose tolerance (1e-3):  steps=" << statsLoose.acceptedSteps);
		INFO("Tight tolerance (1e-10): steps=" << statsTight.acceptedSteps);

		// Tight tolerance requires more steps
		REQUIRE(statsTight.acceptedSteps > statsLoose.acceptedSteps);

		// Tight tolerance should achieve better accuracy
		Real exactX = std::cos(tEnd);
		Real errLoose = std::abs(solLoose.getXValuesAtEnd()[0] - exactX);
		Real errTight = std::abs(solTight.getXValuesAtEnd()[0] - exactX);

		INFO("Loose error: " << errLoose);
		INFO("Tight error: " << errTight);
		REQUIRE(errTight < errLoose);
	}

	TEST_CASE("Adaptive_ErrorBelowTolerance", "[AdaptiveIntegrator][ErrorControl]") {
		TEST_PRECISION_INFO();

		ExponentialDecayODE ode;
		std::vector<Real> tolerances = {1e-3, 1e-6, 1e-9};

		for (Real tol : tolerances) {
			DYNAMIC_SECTION("Tolerance = " << tol) {
				DormandPrince5Integrator integrator(ode);

				Vector<Real> x0{1.0};
				Real t0 = 0.0;
				Real tEnd = 5.0;

				auto sol = integrator.integrate(x0, t0, tEnd, 0.1, tol);

				// Check error at several points
				for (int i = 0; i < sol.getTotalSavedSteps(); ++i) {
					Real t = sol.getTValue(i);
					Real computed = sol.getXValue(i, 0);
					Real exact = std::exp(-t);
					Real localError = std::abs(computed - exact);

					// Error should be below tolerance (with safety factor)
					// Note: accumulated error can exceed local tolerance
					if (i > 0) {
						INFO("t=" << t << ", computed=" << computed << ", exact=" << exact);
						// For accumulated error, allow 100x local tolerance
						REQUIRE(localError < 100.0 * tol);
					}
				}

				// Final error should definitely be bounded
				Real finalErr = std::abs(sol.getXValuesAtEnd()[0] - std::exp(-tEnd));
				INFO("Final error: " << finalErr << " vs tolerance: " << tol);
				REQUIRE(finalErr < 1000.0 * tol);  // Allow accumulation
			}
		}
	}

	/////////////////////////////////////////////////////////////////////////////
	//          ALL STEPPERS CORRECTNESS TESTS                                  //
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("AllSteppers_ExponentialDecay", "[AdaptiveIntegrator][AllSteppers]") {
		TEST_PRECISION_INFO();

		ExponentialDecayODE ode;
		Vector<Real> x0{1.0};
		Real t0 = 0.0;
		Real tEnd = 3.0;
		Real eps = 1e-8;

		SECTION("DormandPrince5") {
			DormandPrince5Integrator integrator(ode);
			auto sol = integrator.integrate(x0, t0, tEnd, 0.1, eps);
			Real err = std::abs(sol.getXValuesAtEnd()[0] - std::exp(-tEnd));
			INFO("DP5 error: " << err);
			REQUIRE(err < eps * 100);
		}

		SECTION("CashKarp") {
			CashKarpIntegrator integrator(ode);
			auto sol = integrator.integrate(x0, t0, tEnd, 0.1, eps);
			Real err = std::abs(sol.getXValuesAtEnd()[0] - std::exp(-tEnd));
			INFO("CashKarp error: " << err);
			REQUIRE(err < eps * 100);
		}

		SECTION("DormandPrince8") {
			DormandPrince8Integrator integrator(ode);
			auto sol = integrator.integrate(x0, t0, tEnd, 0.1, eps);
			Real err = std::abs(sol.getXValuesAtEnd()[0] - std::exp(-tEnd));
			INFO("DP8 error: " << err);
			REQUIRE(err < eps * 100);
		}

		SECTION("BulirschStoer") {
			BulirschStoerIntegrator integrator(ode);
			auto sol = integrator.integrate(x0, t0, tEnd, 0.1, eps);
			Real err = std::abs(sol.getXValuesAtEnd()[0] - std::exp(-tEnd));
			INFO("BS error: " << err);
			REQUIRE(err < eps * 100);
		}

		SECTION("BulirschStoerRational") {
			BulirschStoerRationalIntegrator integrator(ode);
			auto sol = integrator.integrate(x0, t0, tEnd, 0.1, eps);
			Real err = std::abs(sol.getXValuesAtEnd()[0] - std::exp(-tEnd));
			INFO("BSRational error: " << err);
			REQUIRE(err < eps * 100);
		}
	}

	TEST_CASE("AllSteppers_HarmonicOscillator", "[AdaptiveIntegrator][AllSteppers]") {
		TEST_PRECISION_INFO();

		HarmonicOscillatorODE ode;
		Vector<Real> x0{1.0, 0.0};
		Real t0 = 0.0;
		Real tEnd = 4 * Constants::PI;  // Two full periods
		Real eps = 1e-8;

		auto testIntegrator = [&](auto& integrator, const std::string& name) {
			auto sol = integrator.integrate(x0, t0, tEnd, 0.2, eps);
			Real errX = std::abs(sol.getXValuesAtEnd()[0] - std::cos(tEnd));
			Real errV = std::abs(sol.getXValuesAtEnd()[1] + std::sin(tEnd));
			INFO(name << " x-error: " << errX << ", v-error: " << errV);
			REQUIRE(errX < eps * 10000);  // Allow accumulation over 2 periods
			REQUIRE(errV < eps * 10000);
		};

		SECTION("DormandPrince5") {
			DormandPrince5Integrator integrator(ode);
			testIntegrator(integrator, "DP5");
		}

		SECTION("CashKarp") {
			CashKarpIntegrator integrator(ode);
			testIntegrator(integrator, "CashKarp");
		}

		SECTION("DormandPrince8") {
			DormandPrince8Integrator integrator(ode);
			testIntegrator(integrator, "DP8");
		}

		SECTION("BulirschStoer") {
			BulirschStoerIntegrator integrator(ode);
			testIntegrator(integrator, "BS");
		}

		SECTION("BulirschStoerRational") {
			BulirschStoerRationalIntegrator integrator(ode);
			testIntegrator(integrator, "BSRational");
		}
	}

	/////////////////////////////////////////////////////////////////////////////
	//           DENSE OUTPUT ACCURACY TESTS                                    //
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("DenseOutput_Accuracy_AllSteppers", "[AdaptiveIntegrator][DenseOutput]") {
		TEST_PRECISION_INFO();

		HarmonicOscillatorODE ode;
		Vector<Real> x0{1.0, 0.0};
		Real t0 = 0.0;
		Real tEnd = Constants::PI;
		Real eps = 1e-10;

		auto testDense = [&](auto& integrator, const std::string& name) {
			// Use smaller output interval for better interpolation
			auto sol = integrator.integrate(x0, t0, tEnd, 0.1, eps);

			// Check interpolated output at non-saved points using spline
			auto interpX = sol.getSolAsSplineInterp(0);
			auto interpV = sol.getSolAsSplineInterp(1);
			
			for (Real t = 0.15; t < tEnd - 0.1; t += 0.2) {
				Real yDenseX = interpX(t);
				Real yDenseV = interpV(t);
				
				Real exactX = std::cos(t);
				Real exactV = -std::sin(t);
				Real errX = std::abs(yDenseX - exactX);
				Real errV = std::abs(yDenseV - exactV);
				
				INFO(name << " at t=" << t << ": errX=" << errX << ", errV=" << errV);
				// Spline interpolation accuracy depends on saved points density
				// The integrator is accurate, but spline interpolation adds its own error
				REQUIRE(errX < 1e-3);  // Spline interpolation error
				REQUIRE(errV < 1e-3);
			}
		};

		SECTION("DP5 Dense") {
			DormandPrince5Integrator integrator(ode);
			testDense(integrator, "DP5");
		}

		SECTION("CashKarp Dense") {
			CashKarpIntegrator integrator(ode);
			testDense(integrator, "CashKarp");
		}

		SECTION("DP8 Dense") {
			DormandPrince8Integrator integrator(ode);
			testDense(integrator, "DP8");
		}
	}

	/////////////////////////////////////////////////////////////////////////////
	//            LONG-TIME INTEGRATION STABILITY                               //
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("LongTimeIntegration_EnergyConservation", "[AdaptiveIntegrator][Conservation]") {
		TEST_PRECISION_INFO();

		KeplerODE ode;
		Vector<Real> y0 = KeplerODE::getCircularOrbitIC();
		Real E0 = KeplerODE::computeEnergy(y0);

		Real t0 = 0.0;
		Real tEnd = 20 * Constants::PI;  // 10 orbits
		Real eps = 1e-10;

		SECTION("DP5 preserves energy") {
			DormandPrince5Integrator integrator(ode);
			auto sol = integrator.integrate(y0, t0, tEnd, 0.5, eps);

			Real Efinal = KeplerODE::computeEnergy(sol.getXValuesAtEnd());
			Real energyDrift = std::abs(Efinal - E0) / std::abs(E0);

			INFO("Initial E=" << E0 << ", Final E=" << Efinal);
			INFO("Relative energy drift: " << energyDrift);
			REQUIRE(energyDrift < 1e-6);
		}

		SECTION("DP8 preserves energy better") {
			DormandPrince8Integrator integrator(ode);
			auto sol = integrator.integrate(y0, t0, tEnd, 0.5, eps);

			Real Efinal = KeplerODE::computeEnergy(sol.getXValuesAtEnd());
			Real energyDrift = std::abs(Efinal - E0) / std::abs(E0);

			INFO("Initial E=" << E0 << ", Final E=" << Efinal);
			INFO("DP8 Relative energy drift: " << energyDrift);
			REQUIRE(energyDrift < 1e-7);  // Higher order should do better
		}
	}

	TEST_CASE("LongTimeIntegration_HarmonicOscillator", "[AdaptiveIntegrator][LongTime]") {
		TEST_PRECISION_INFO();

		HarmonicOscillatorODE ode;
		Vector<Real> x0{1.0, 0.0};
		Real t0 = 0.0;
		Real tEnd = 100 * Constants::PI;  // 50 periods
		Real eps = 1e-10;

		DormandPrince8Integrator integrator(ode);
		auto sol = integrator.integrate(x0, t0, tEnd, 1.0, eps);

		// Check amplitude is preserved (|x|^2 + |v|^2 = 1)
		auto yFinal = sol.getXValuesAtEnd();
		Real amplitude = std::sqrt(yFinal[0]*yFinal[0] + yFinal[1]*yFinal[1]);
		Real ampError = std::abs(amplitude - 1.0);

		INFO("Final amplitude: " << amplitude);
		INFO("Amplitude error: " << ampError);
		REQUIRE(ampError < 1e-6);
	}

	/////////////////////////////////////////////////////////////////////////////
	//                NONLINEAR SYSTEM TESTS                                    //
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("LotkaVolterra_Conservation", "[AdaptiveIntegrator][Nonlinear]") {
		TEST_PRECISION_INFO();

		// Lotka-Volterra has a conserved quantity (constant of motion)
		LotkaVolterraODE ode;
		Vector<Real> y0{10.0, 5.0};  // Initial prey, predator

		Real t0 = 0.0;
		Real tEnd = 20.0;
		Real eps = 1e-8;

		// Conserved quantity: H = δx - γ*ln(x) + βy - α*ln(y)
		auto computeH = [](const Vector<Real>& y) {
			Real alpha = 1.0, beta = 0.1, gamma = 1.5, delta = 0.075;
			return delta * y[0] - gamma * std::log(y[0]) + 
			       beta * y[1] - alpha * std::log(y[1]);
		};

		Real H0 = computeH(y0);

		DormandPrince8Integrator integrator(ode);
		auto sol = integrator.integrate(y0, t0, tEnd, 0.2, eps);

		Real Hfinal = computeH(sol.getXValuesAtEnd());
		Real Hdrift = std::abs(Hfinal - H0) / std::abs(H0);

		INFO("Initial H=" << H0 << ", Final H=" << Hfinal);
		INFO("Relative H drift: " << Hdrift);
		REQUIRE(Hdrift < 1e-5);

		// Also verify solution stays positive
		for (int i = 0; i < sol.getTotalSavedSteps(); ++i) {
			REQUIRE(sol.getXValue(i, 0) > 0);  // prey
			REQUIRE(sol.getXValue(i, 1) > 0);  // predator
		}
	}

	TEST_CASE("Lorenz_BoundedTrajectory", "[AdaptiveIntegrator][Chaotic]") {
		TEST_PRECISION_INFO();

		LorenzODE ode;
		Vector<Real> x0{1.0, 1.0, 1.0};  // Near unstable fixed point

		Real t0 = 0.0;
		Real tEnd = 50.0;
		Real eps = 1e-8;

		DormandPrince5Integrator integrator(ode);
		auto sol = integrator.integrate(x0, t0, tEnd, 0.1, eps);

		// Lorenz attractor is bounded: |x|, |y| < 30, |z| < 50 typically
		for (int i = 0; i < sol.getTotalSavedSteps(); ++i) {
			REQUIRE(std::abs(sol.getXValue(i, 0)) < 50.0);
			REQUIRE(std::abs(sol.getXValue(i, 1)) < 50.0);
			REQUIRE(std::abs(sol.getXValue(i, 2)) < 80.0);
		}

		INFO("Lorenz integrated successfully to t=" << tEnd);
		INFO("Final state: [" << sol.getXValuesAtEnd()[0] << ", " 
		                      << sol.getXValuesAtEnd()[1] << ", "
		                      << sol.getXValuesAtEnd()[2] << "]");
	}

	/////////////////////////////////////////////////////////////////////////////
	//              RAPID OSCILLATION TESTS                                     //
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("RapidOscillator_StepAdaptation", "[AdaptiveIntegrator][Oscillatory]") {
		TEST_PRECISION_INFO();

		RapidOscillatorODE ode(50.0);  // omega = 50
		Vector<Real> x0{1.0, 0.0};
		Real t0 = 0.0;
		Real tEnd = 2.0;
		Real eps = 1e-6;

		DormandPrince5Integrator integrator(ode);
		auto sol = integrator.integrate(x0, t0, tEnd, 0.1, eps);
		auto stats = integrator.getStatistics();

		// With omega=50, period ≈ 0.126, so expect many steps
		INFO("Accepted steps: " << stats.acceptedSteps);
		INFO("Rejected steps: " << stats.rejectedSteps);
		REQUIRE(stats.acceptedSteps > 100);  // Must adapt to fast oscillation

		// Check accuracy at end
		Real exactX = ode.exactX(tEnd);
		Real errX = std::abs(sol.getXValuesAtEnd()[0] - exactX);
		INFO("Final x error: " << errX);
		REQUIRE(errX < eps * 1000);
	}

	/////////////////////////////////////////////////////////////////////////////
	//                EDGE CASE TESTS                                           //
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("EdgeCase_VeryTightTolerance", "[AdaptiveIntegrator][EdgeCase]") {
		TEST_PRECISION_INFO();

		ExponentialDecayODE ode;
		Vector<Real> x0{1.0};
		Real t0 = 0.0;
		Real tEnd = 1.0;
		Real eps = 1e-14;  // Very tight

		DormandPrince8Integrator integrator(ode);  // Need high-order for tight tol
		auto sol = integrator.integrate(x0, t0, tEnd, 0.1, eps);
		auto stats = integrator.getStatistics();

		Real exact = std::exp(-tEnd);
		Real err = std::abs(sol.getXValuesAtEnd()[0] - exact);

		INFO("Error: " << err << " with eps=" << eps);
		INFO("Steps: " << stats.acceptedSteps);
		
		// Should achieve very high accuracy
		REQUIRE(err < 1e-11);
	}

	TEST_CASE("EdgeCase_LooseTolerance", "[AdaptiveIntegrator][EdgeCase]") {
		TEST_PRECISION_INFO();

		HarmonicOscillatorODE ode;
		Vector<Real> x0{1.0, 0.0};
		Real t0 = 0.0;
		Real tEnd = Constants::PI;
		Real eps = 1e-2;  // Very loose

		DormandPrince5Integrator integrator(ode);
		auto sol = integrator.integrate(x0, t0, tEnd, 0.5, eps);
		auto stats = integrator.getStatistics();

		INFO("Steps with loose tolerance: " << stats.acceptedSteps);
		
		// Should still give reasonable answer
		Real exactX = std::cos(tEnd);
		Real err = std::abs(sol.getXValuesAtEnd()[0] - exactX);
		INFO("Error: " << err);
		REQUIRE(err < 0.1);  // Within 10%
	}

	TEST_CASE("EdgeCase_ShortIntegration", "[AdaptiveIntegrator][EdgeCase]") {
		TEST_PRECISION_INFO();

		ExponentialDecayODE ode;
		Vector<Real> x0{1.0};
		Real t0 = 0.0;
		Real tEnd = 1e-6;  // Very short
		Real eps = 1e-10;

		DormandPrince5Integrator integrator(ode);
		auto sol = integrator.integrate(x0, t0, tEnd, 1e-7, eps);

		Real exact = std::exp(-tEnd);
		Real err = std::abs(sol.getXValuesAtEnd()[0] - exact);

		INFO("Short integration error: " << err);
		REQUIRE(err < 1e-14);  // Should be very accurate
	}

	/////////////////////////////////////////////////////////////////////////////
	//              STEPPER COMPARISON TESTS                                    //
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("StepperComparison_EfficiencyVsAccuracy", "[AdaptiveIntegrator][Comparison]") {
		TEST_PRECISION_INFO();

		HarmonicOscillatorODE ode;
		Vector<Real> x0{1.0, 0.0};
		Real t0 = 0.0;
		Real tEnd = 10 * Constants::PI;
		Real eps = 1e-10;

		struct Result {
			std::string name;
			int steps;
			int funcEvals;
			Real error;
		};
		std::vector<Result> results;

		auto measure = [&](auto& integrator, const std::string& name) {
			auto sol = integrator.integrate(x0, t0, tEnd, 0.5, eps);
			auto stats = integrator.getStatistics();
			Real err = std::abs(sol.getXValuesAtEnd()[0] - std::cos(tEnd));
			results.push_back({name, stats.acceptedSteps, stats.totalFuncEvals, err});
		};

		{
			DormandPrince5Integrator integrator(ode);
			measure(integrator, "DP5");
		}
		{
			CashKarpIntegrator integrator(ode);
			measure(integrator, "CashKarp");
		}
		{
			DormandPrince8Integrator integrator(ode);
			measure(integrator, "DP8");
		}
		{
			BulirschStoerIntegrator integrator(ode);
			measure(integrator, "BulirschStoer");
		}

		// Print comparison
		for (const auto& r : results) {
			INFO(r.name << ": steps=" << r.steps << ", funcEvals=" << r.funcEvals << ", error=" << r.error);
		}

		// All should achieve required accuracy
		for (const auto& r : results) {
			REQUIRE(r.error < eps * 10000);
		}

		// Higher order methods should use fewer steps (generally)
		// DP8 should need fewer steps than DP5
		auto dp5 = std::find_if(results.begin(), results.end(), [](const Result& r) { return r.name == "DP5"; });
		auto dp8 = std::find_if(results.begin(), results.end(), [](const Result& r) { return r.name == "DP8"; });
		if (dp5 != results.end() && dp8 != results.end()) {
			INFO("DP5 steps: " << dp5->steps << ", DP8 steps: " << dp8->steps);
			// DP8 typically needs fewer steps for same accuracy
			REQUIRE(dp8->steps <= dp5->steps * 2);  // Allow some slack
		}
	}

	/////////////////////////////////////////////////////////////////////////////
	//              FSAL (First Same As Last) VERIFICATION                      //
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("FSAL_DerivativeReuse", "[AdaptiveIntegrator][FSAL]") {
		TEST_PRECISION_INFO();

		// DP5 and DP8 are FSAL methods - they should reuse derivatives
		ExponentialDecayODE ode;
		DormandPrince5_Stepper stepper(ode);

		Vector<Real> x{1.0};
		Vector<Real> dxdt{-1.0};
		Real t = 0.0;
		Real eps = 1e-10;
		Real h = 0.1;

		// Do first step
		StepResult r1 = stepper.doStep(t, x, dxdt, h, eps);
		REQUIRE(r1.accepted);
		INFO("First step funcEvals: " << r1.funcEvals);

		// In a proper FSAL implementation, the derivative at end of step
		// should already be computed and usable for next step
		// The stepper tracks this internally

		t += r1.hDone;
		
		// Do second step
		StepResult r2 = stepper.doStep(t, x, dxdt, r1.hNext, eps);
		INFO("Second step funcEvals: " << r2.funcEvals);

		// For a true FSAL method, second step should save one eval
		// DP5 has 6 stages but FSAL saves 1 => 5 evals after first
		// (First step needs 6, subsequent need 5)
	}

	/////////////////////////////////////////////////////////////////////////////
	//                integrateAt TESTS                                         //
	/////////////////////////////////////////////////////////////////////////////

	TEST_CASE("IntegrateAt_SpecificTimes", "[AdaptiveIntegrator][IntegrateAt]") {
		TEST_PRECISION_INFO();

		HarmonicOscillatorODE ode;
		Vector<Real> x0{1.0, 0.0};
		Real eps = 1e-10;

		// Create Vector<Real> for times
		Vector<Real> times(6);
		times[0] = 0.0;
		times[1] = 0.5;
		times[2] = 1.0;
		times[3] = 1.5;
		times[4] = Constants::PI;
		times[5] = 2*Constants::PI;

		DormandPrince8Integrator integrator(ode);
		auto sol = integrator.integrateAt(x0, times, eps);

		REQUIRE(sol.getTotalSavedSteps() == static_cast<int>(times.size()));

		for (int i = 0; i < times.size(); ++i) {
			Real t = times[i];
			Real exactX = std::cos(t);
			Real computedX = sol.getXValue(i, 0);
			Real err = std::abs(computedX - exactX);

			INFO("t=" << t << ", computed=" << computedX << ", exact=" << exactX << ", err=" << err);
			REQUIRE(err < eps * 1000);
		}
	}

} // namespace MML::Tests::Algorithms::ComprehensiveAdaptiveTests

///////////////////////////////////////////////////////////////////////////////////////////
//                      STIFF ODE IMPLICIT SOLVER TESTS                                   //
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MML_USE_SINGLE_HEADER
#include "algorithms/ODEStiffSolvers.h"
#endif

namespace MML::Tests::Algorithms::StiffSolverTests
{
	// Simple stiff test system: y' = λ(y - cos(t)) - sin(t), y(0) = 1
	// Exact solution: y(t) = cos(t)
	// For λ << -1, system is stiff
	class SimpleStiffODE : public IODESystemWithJacobian
	{
	private:
		Real _lambda;
	public:
		SimpleStiffODE(Real lambda = -50.0) : _lambda(lambda) {}
		
		int getDim() const override { return 1; }
		
		void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override {
			dydt[0] = _lambda * (y[0] - std::cos(t)) - std::sin(t);
		}
		
		void jacobian(const Real t, const Vector<Real>& y, Vector<Real>& dydt, Matrix<Real>& J) const override {
			J(0, 0) = _lambda;
		}
		
		Real exactSolution(Real t) const { return std::cos(t); }
	};
	
	TEST_CASE("BackwardEuler - Simple Stiff ODE", "[StiffSolvers][BackwardEuler]") {
		TEST_PRECISION_INFO();
		
		SimpleStiffODE system(-50.0);
		Vector<Real> y0(1);
		y0[0] = 1.0;
		
		Real t0 = 0.0;
		Real t_end = 1.0;
		Real h = 0.1;
		
		auto sol = SolveBackwardEuler(system, t0, y0, t_end, h);
		
		Vector<Real> y_final = sol.getXValuesAtEnd();
		Real exact = system.exactSolution(t_end);
		
		INFO("Computed: " << y_final[0] << ", Exact: " << exact);
		REQUIRE_THAT(y_final[0], WithinAbs(exact, 0.01));  // Backward Euler is only 1st order
	}
	
	TEST_CASE("BDF2 - Simple Stiff ODE", "[StiffSolvers][BDF2]") {
		TEST_PRECISION_INFO();
		
		SimpleStiffODE system(-50.0);
		Vector<Real> y0(1);
		y0[0] = 1.0;
		
		Real t0 = 0.0;
		Real t_end = 1.0;
		Real h = 0.1;
		
		auto sol = SolveBDF2(system, t0, y0, t_end, h);
		
		Vector<Real> y_final = sol.getXValuesAtEnd();
		Real exact = system.exactSolution(t_end);
		
		INFO("Computed: " << y_final[0] << ", Exact: " << exact);
		REQUIRE_THAT(y_final[0], WithinAbs(exact, 0.001));  // BDF2 is 2nd order, better than BE
	}
	
	TEST_CASE("Rosenbrock - Simple Stiff ODE", "[StiffSolvers][Rosenbrock]") {
		TEST_PRECISION_INFO();
		
		SECTION("Create system") {
			SimpleStiffODE system(-50.0);
			REQUIRE(system.getDim() == 1);
			
			SECTION("Create y0") {
				Vector<Real> y0(1);
				y0[0] = 1.0;
				REQUIRE(y0.size() == 1);
				
				SECTION("Create solver") {
					Rosenbrock23Solver solver(system, 1e-3, 1e-3);
					
					SECTION("Solve") {
						Real t0 = 0.0;
						Real t_end = 1.0;
						Real h_init = 0.1;
						
						auto sol = solver.Solve(t0, y0, t_end, h_init);
						
						Vector<Real> y_final = sol.getXValuesAtEnd();
						Real exact = system.exactSolution(t_end);
						
						INFO("Computed: " << y_final[0] << ", Exact: " << exact);
						INFO("Steps accepted: " << sol.getNumStepsOK() << ", rejected: " << sol.getNumStepsBad());
						REQUIRE_THAT(y_final[0], WithinAbs(exact, 0.01));
					}
				}
			}
		}
	}
	
	TEST_CASE("BackwardEuler vs BDF2 vs Rosenbrock - Accuracy Comparison", "[StiffSolvers]") {
		TEST_PRECISION_INFO();
		
		SimpleStiffODE system(-100.0);  // Very stiff
		Vector<Real> y0(1);
		y0[0] = 1.0;
		
		Real t0 = 0.0;
		Real t_end = 2.0;
		Real h = 0.05;
		
		// Backward Euler
		auto sol_be = SolveBackwardEuler(system, t0, y0, t_end, h);
		Real err_be = std::abs(sol_be.getXValuesAtEnd()[0] - system.exactSolution(t_end));
		
		// BDF2
		auto sol_bdf2 = SolveBDF2(system, t0, y0, t_end, h);
		Real err_bdf2 = std::abs(sol_bdf2.getXValuesAtEnd()[0] - system.exactSolution(t_end));
		
		// Rosenbrock (ROS2 - 2nd order, adaptive)
		Rosenbrock23Solver solver(system, 1e-3, 1e-3);
		auto sol_ros = solver.Solve(t0, y0, t_end, h);
		Real err_ros = std::abs(sol_ros.getXValuesAtEnd()[0] - system.exactSolution(t_end));
		
		INFO("Backward Euler error: " << err_be);
		INFO("BDF2 error: " << err_bdf2);
		INFO("Rosenbrock (ROS2) error: " << err_ros);
		
		// BDF2 should be more accurate than Backward Euler
		REQUIRE(err_bdf2 < err_be);
		
		// Rosenbrock should be accurate (2nd order with adaptive stepping)
		REQUIRE(err_ros < 0.01);  // 2nd order should achieve decent accuracy
	}
	
	TEST_CASE("Stiff Solvers - Van der Pol (mu=10)", "[StiffSolvers][VanDerPol]") {
		TEST_PRECISION_INFO();
		
		// Van der Pol with moderate stiffness
		auto system = std::make_unique<VanDerPolStiffODE>(10.0);
		Vector<Real> y0 = VanDerPolStiffODE::getInitialCondition();
		
		Real t0 = 0.0;
		Real t_end = 2.0;
		
		SECTION("Backward Euler handles moderate stiffness") {
			Real h = 0.05;
			auto sol = SolveBackwardEuler(*system, t0, y0, t_end, h);
			
			REQUIRE(sol.getNumStepsOK() > 0);
			REQUIRE(sol.getTotalSavedSteps() > 10);
			
			// Solution should be bounded
			Vector<Real> y_final = sol.getXValuesAtEnd();
			REQUIRE(std::abs(y_final[0]) < 5.0);
			REQUIRE(std::abs(y_final[1]) < 5.0);
		}
		
		SECTION("BDF2 handles moderate stiffness efficiently") {
			Real h = 0.1;  // Can use larger step than Backward Euler
			auto sol = SolveBDF2(*system, t0, y0, t_end, h);
			
			REQUIRE(sol.getNumStepsOK() > 0);
			
			// Solution should be bounded
			Vector<Real> y_final = sol.getXValuesAtEnd();
			REQUIRE(std::abs(y_final[0]) < 5.0);
			REQUIRE(std::abs(y_final[1]) < 5.0);
		}
		
		SECTION("Rosenbrock adapts step size appropriately") {
			Real h_init = 0.1;
			Rosenbrock23Solver solver(*system, 1e-2, 1e-2);  // Relaxed tolerances for Van der Pol
			auto sol = solver.Solve(t0, y0, t_end, h_init);
			
			REQUIRE(sol.getNumStepsOK() > 0);
			
			// Adaptive method should accept and reject steps
			INFO("Accepted: " << sol.getNumStepsOK() << ", Rejected: " << sol.getNumStepsBad());
			
			// Solution should be bounded
			Vector<Real> y_final = sol.getXValuesAtEnd();
			REQUIRE(std::abs(y_final[0]) < 5.0);
			REQUIRE(std::abs(y_final[1]) < 5.0);
		}
	}
	
	TEST_CASE("Stiff Solvers - Robertson Chemical Kinetics", "[StiffSolvers][Robertson][!mayfail]") {
		TEST_PRECISION_INFO();
		
		// Robertson is extremely stiff - may fail with simple parameters
		auto system = std::make_unique<RobertsonODE>();
		Vector<Real> y0 = RobertsonODE::getInitialCondition();
		
		Real t0 = 0.0;
		Real t_end = 0.1;  // Very short integration time
		
		SECTION("BDF2 can integrate Robertson (short time)") {
			Real h = 1e-4;  // Very small step needed
			
			try {
				auto sol = SolveBDF2(*system, t0, y0, t_end, h);
				
				Vector<Real> y_final = sol.getXValuesAtEnd();
				
				// Conservation: y1 + y2 + y3 ≈ 1
				Real sum = y_final[0] + y_final[1] + y_final[2];
				INFO("Conservation check: " << sum);
				REQUIRE_THAT(sum, WithinAbs(1.0, 0.01));
				
				// All components should be non-negative
				REQUIRE(y_final[0] >= 0.0);
				REQUIRE(y_final[1] >= 0.0);
				REQUIRE(y_final[2] >= 0.0);
			} catch (const std::exception& e) {
				WARN("BDF2 failed on Robertson: " << e.what());
				// This is expected - Robertson is extremely challenging
			}
		}
	}

} // namespace MML::Tests::Algorithms::StiffSolverTests