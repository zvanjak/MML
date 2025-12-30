#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/DiracDeltaFunction.h"
#endif

using namespace MML;
using namespace MML::Testing;
using namespace Catch::Matchers;

namespace MML::Tests::Base::DiracDeltaFunctionTests
{
	/*********************************************************************/
	/*****                DiracStep Tests                           *****/
	/*********************************************************************/
	TEST_CASE("DiracStep::Constructor", "[DiracDeltaFunction][DiracStep]")
	{
			TEST_PRECISION_INFO();
		DiracStep delta10(10);
		DiracStep delta100(100);
		
		// Constructor should accept N parameter
		REQUIRE_NOTHROW(DiracStep(1));
		REQUIRE_NOTHROW(DiracStep(1000));
	}

	TEST_CASE("DiracStep::ZeroOutsideInterval", "[DiracDeltaFunction][DiracStep]")
	{
			TEST_PRECISION_INFO();
		DiracStep delta(10);
		Real halfInterval = REAL(1.0) / (2 * 10);  // REAL(0.05)
		
		// Should be zero outside [-1/(2N), 1/(2N)]
		REQUIRE_THAT(delta(-REAL(0.1)), WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(delta(-halfInterval - REAL(0.01)), WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(delta(halfInterval + REAL(0.01)), WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(delta(REAL(0.1)), WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(delta(REAL(1.0)), WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(delta(-REAL(1.0)), WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("DiracStep::NonzeroInsideInterval", "[DiracDeltaFunction][DiracStep]")
	{
			TEST_PRECISION_INFO();
		DiracStep delta(10);
		
		// Should be N inside the interval
		REQUIRE_THAT(delta(REAL(0.0)), WithinAbs(REAL(10.0), REAL(1e-10)));
		REQUIRE_THAT(delta(REAL(0.01)), WithinAbs(REAL(10.0), REAL(1e-10)));
		REQUIRE_THAT(delta(-REAL(0.01)), WithinAbs(REAL(10.0), REAL(1e-10)));
		REQUIRE_THAT(delta(REAL(0.049)), WithinAbs(REAL(10.0), REAL(1e-10)));  // Just inside
		REQUIRE_THAT(delta(-REAL(0.049)), WithinAbs(REAL(10.0), REAL(1e-10))); // Just inside
	}

	TEST_CASE("DiracStep::IncreasingN", "[DiracDeltaFunction][DiracStep]")
	{
			TEST_PRECISION_INFO();
		// As N increases, the function becomes taller and narrower
		DiracStep delta10(10);
		DiracStep delta100(100);
		
		// At x=0, height increases with N
		REQUIRE_THAT(delta10(REAL(0.0)), WithinAbs(REAL(10.0), REAL(1e-10)));
		REQUIRE_THAT(delta100(REAL(0.0)), WithinAbs(REAL(100.0), REAL(1e-10)));
		
		// Width decreases with N
		REQUIRE_THAT(delta10(REAL(0.06)), WithinAbs(REAL(0.0), REAL(1e-10)));   // Outside for N=10
		REQUIRE_THAT(delta100(REAL(0.006)), WithinAbs(REAL(0.0), REAL(1e-10)));  // Outside for N=100
	}

	/*********************************************************************/
	/*****                DiracExp Tests                            *****/
	/*********************************************************************/
	TEST_CASE("DiracExp::Constructor", "[DiracDeltaFunction][DiracExp]")
	{
			TEST_PRECISION_INFO();
		REQUIRE_NOTHROW(DiracExp(1));
		REQUIRE_NOTHROW(DiracExp(10));
		REQUIRE_NOTHROW(DiracExp(100));
	}

	TEST_CASE("DiracExp::MaximumAtOrigin", "[DiracDeltaFunction][DiracExp]")
	{
			TEST_PRECISION_INFO();
		DiracExp delta(10);
		
		// Maximum should be at x=0
		Real valCenter = delta(REAL(0.0));
		Real valLeft = delta(-REAL(0.01));
		Real valRight = delta(REAL(0.01));
		
		REQUIRE(valCenter > valLeft);
		REQUIRE(valCenter > valRight);
		
		// Value at origin: N / sqrt(2*PI) * exp(0) = N / sqrt(2*PI)
		Real expected = REAL(10.0) / std::sqrt(2 * Constants::PI);
		REQUIRE_THAT(valCenter, WithinAbs(expected, REAL(1e-10)));
	}

	TEST_CASE("DiracExp::Symmetry", "[DiracDeltaFunction][DiracExp]")
	{
			TEST_PRECISION_INFO();
		DiracExp delta(10);
		
		// Should be symmetric around x=0
		REQUIRE_THAT(delta(REAL(0.1)), WithinAbs(delta(-REAL(0.1)), REAL(1e-10)));
		REQUIRE_THAT(delta(REAL(0.5)), WithinAbs(delta(-REAL(0.5)), REAL(1e-10)));
		REQUIRE_THAT(delta(REAL(1.0)), WithinAbs(delta(-REAL(1.0)), REAL(1e-10)));
	}

	TEST_CASE("DiracExp::DecaysAwayFromOrigin", "[DiracDeltaFunction][DiracExp]")
	{
			TEST_PRECISION_INFO();
		DiracExp delta(10);
		
		// Should decay exponentially
		Real val0 = delta(REAL(0.0));
		Real val1 = delta(REAL(0.1));
		Real val2 = delta(REAL(0.2));
		Real val3 = delta(REAL(0.5));
		
		REQUIRE(val0 > val1);
		REQUIRE(val1 > val2);
		REQUIRE(val2 > val3);
		REQUIRE(val3 > REAL(0.0));
	}

	TEST_CASE("DiracExp::IncreasingN", "[DiracDeltaFunction][DiracExp]")
	{
			TEST_PRECISION_INFO();
		DiracExp delta10(10);
		DiracExp delta100(100);
		
		// Peak height increases with N
		REQUIRE(delta100(REAL(0.0)) > delta10(REAL(0.0)));
		
		// Width decreases with N (narrower peak)
		// At same distance from origin, larger N should have smaller value
		REQUIRE(delta100(REAL(0.1)) < delta10(REAL(0.1)));
	}

	/*********************************************************************/
	/*****                DiracSqr Tests                            *****/
	/*********************************************************************/
	TEST_CASE("DiracSqr::Constructor", "[DiracDeltaFunction][DiracSqr]")
	{
			TEST_PRECISION_INFO();
		REQUIRE_NOTHROW(DiracSqr(1));
		REQUIRE_NOTHROW(DiracSqr(10));
		REQUIRE_NOTHROW(DiracSqr(100));
	}

	TEST_CASE("DiracSqr::MaximumAtOrigin", "[DiracDeltaFunction][DiracSqr]")
	{
			TEST_PRECISION_INFO();
		DiracSqr delta(10);
		
		// Maximum at x=0: N / PI / (1 + 0) = N / PI
		Real valCenter = delta(REAL(0.0));
		Real expected = REAL(10.0) / Constants::PI;
		REQUIRE_THAT(valCenter, WithinAbs(expected, REAL(1e-10)));
		
		// Should be maximum
		REQUIRE(valCenter > delta(REAL(0.01)));
		REQUIRE(valCenter > delta(REAL(0.1)));
	}

	TEST_CASE("DiracSqr::Symmetry", "[DiracDeltaFunction][DiracSqr]")
	{
			TEST_PRECISION_INFO();
		DiracSqr delta(10);
		
		// Lorentzian is symmetric
		REQUIRE_THAT(delta(REAL(0.1)), WithinAbs(delta(-REAL(0.1)), REAL(1e-10)));
		REQUIRE_THAT(delta(REAL(0.5)), WithinAbs(delta(-REAL(0.5)), REAL(1e-10)));
		REQUIRE_THAT(delta(REAL(1.0)), WithinAbs(delta(-REAL(1.0)), REAL(1e-10)));
	}

	TEST_CASE("DiracSqr::DecaysAwayFromOrigin", "[DiracDeltaFunction][DiracSqr]")
	{
			TEST_PRECISION_INFO();
		DiracSqr delta(10);
		
		// Lorentzian decay: 1/(1+x^2)
		Real val0 = delta(REAL(0.0));
		Real val1 = delta(REAL(0.1));
		Real val2 = delta(REAL(0.5));
		Real val3 = delta(REAL(1.0));
		
		REQUIRE(val0 > val1);
		REQUIRE(val1 > val2);
		REQUIRE(val2 > val3);
		REQUIRE(val3 > REAL(0.0));
	}

	TEST_CASE("DiracSqr::IncreasingN", "[DiracDeltaFunction][DiracSqr]")
	{
			TEST_PRECISION_INFO();
		DiracSqr delta10(10);
		DiracSqr delta100(100);
		
		// Peak height increases with N
		REQUIRE(delta100(REAL(0.0)) > delta10(REAL(0.0)));
		
		// Narrower peak with larger N
		REQUIRE(delta100(REAL(0.1)) < delta10(REAL(0.1)));
	}

	/*********************************************************************/
	/*****                DiracSin Tests                            *****/
	/*********************************************************************/
	TEST_CASE("DiracSin::Constructor", "[DiracDeltaFunction][DiracSin]")
	{
			TEST_PRECISION_INFO();
		REQUIRE_NOTHROW(DiracSin(1));
		REQUIRE_NOTHROW(DiracSin(10));
		REQUIRE_NOTHROW(DiracSin(100));
	}

	TEST_CASE("DiracSin::CentralPeak", "[DiracDeltaFunction][DiracSin]")
	{
			TEST_PRECISION_INFO();
		DiracSin delta(10);
		
		// At x=0, sinc function has limit: sin(Nx)/(PI*x) -> N/PI as x->0
		// Need to test near zero due to singularity
		Real valNearZero = delta(REAL(0.001));
		Real expected = REAL(10.0) / Constants::PI;
		
		// Should be close to N/PI for small x
		REQUIRE_THAT(valNearZero, WithinAbs(expected, REAL(0.01)));
	}

	TEST_CASE("DiracSin::Oscillatory", "[DiracDeltaFunction][DiracSin]")
	{
			TEST_PRECISION_INFO();
		DiracSin delta(10);
		
		// Sinc function oscillates and decays
		// Zero crossings at x = k*PI/N for integer k
		Real zeroCrossing = Constants::PI / REAL(10.0);
		
		REQUIRE_THAT(delta(zeroCrossing), WithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(delta(2 * zeroCrossing), WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("DiracSin::DecaysAwayFromOrigin", "[DiracDeltaFunction][DiracSin]")
	{
			TEST_PRECISION_INFO();
		DiracSin delta(10);
		
		// Overall envelope should decay like 1/x
		Real val1 = std::abs(delta(REAL(0.1)));
		Real val2 = std::abs(delta(REAL(0.5)));
		Real val3 = std::abs(delta(REAL(1.0)));
		
		// Magnitude should generally decrease (ignoring oscillations)
		REQUIRE(val1 > val3);
	}

	TEST_CASE("DiracSin::IncreasingN", "[DiracDeltaFunction][DiracSin]")
	{
			TEST_PRECISION_INFO();
		DiracSin delta10(10);
		DiracSin delta100(100);
		
		// Peak height increases with N (near origin)
		Real val10 = delta10(REAL(0.001));
		Real val100 = delta100(REAL(0.001));
		
		REQUIRE(val100 > val10);
	}

	/*********************************************************************/
	/*****              Integration Tests                           *****/
	/*********************************************************************/
	TEST_CASE("Integration::AllFunctionsSmoothAtOrigin", "[DiracDeltaFunction][Integration]")
	{
			TEST_PRECISION_INFO();
		// All functions except DiracStep should be continuous
		DiracExp deltaExp(10);
		DiracSqr deltaSqr(10);
		
		Real smallStep = REAL(0.001);
		Real valLeftExp = deltaExp(-smallStep);
		Real valRightExp = deltaExp(smallStep);
		Real valCenterExp = deltaExp(REAL(0.0));
		
		// Exponential should be smooth
		REQUIRE(std::abs(valCenterExp - valLeftExp) < REAL(0.5));
		REQUIRE(std::abs(valCenterExp - valRightExp) < REAL(0.5));
		
		Real valLeftSqr = deltaSqr(-smallStep);
		Real valRightSqr = deltaSqr(smallStep);
		Real valCenterSqr = deltaSqr(REAL(0.0));
		
		// Lorentzian should be smooth
		REQUIRE(std::abs(valCenterSqr - valLeftSqr) < REAL(0.5));
		REQUIRE(std::abs(valCenterSqr - valRightSqr) < REAL(0.5));
	}

	TEST_CASE("Integration::AllFunctionsPositive", "[DiracDeltaFunction][Integration]")
	{
			TEST_PRECISION_INFO();
		// All Dirac approximations should be positive (except DiracSin can be negative)
		DiracStep deltaStep(10);
		DiracExp deltaExp(10);
		DiracSqr deltaSqr(10);
		
		for (Real x = -REAL(1.0); x <= REAL(1.0); x += REAL(0.1)) {
			REQUIRE(deltaStep(x) >= REAL(0.0));
			REQUIRE(deltaExp(x) >= REAL(0.0));
			REQUIRE(deltaSqr(x) >= REAL(0.0));
		}
	}

	TEST_CASE("Integration::NarrowingWithIncreasingN", "[DiracDeltaFunction][Integration]")
	{
			TEST_PRECISION_INFO();
		// All functions should concentrate near origin as N increases
		DiracExp deltaExp1(10);
		DiracExp deltaExp2(100);
		
		DiracSqr deltaSqr1(10);
		DiracSqr deltaSqr2(100);
		
		Real testPoint = REAL(0.2);
		
		// Larger N should give smaller values away from origin
		REQUIRE(deltaExp2(testPoint) < deltaExp1(testPoint));
		REQUIRE(deltaSqr2(testPoint) < deltaSqr1(testPoint));
	}
}
