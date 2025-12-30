///////////////////////////////////////////////////////////////////////////////////////////
///  File:        optimization_framework_tests.cpp
///  Description: Tests for new optimization framework components
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "algorithms/Optimization/OptimizationCommon.h"
#include "algorithms/Optimization/TerminationCriteria.h"
#include "algorithms/Optimization/OptimizationObservers.h"
#include "algorithms/Optimization/OptimizationConfig.h"
#endif

using namespace MML;
using namespace MML::Testing;

///////////////////////////////////////////////////////////////////////////////////////////
// BASIC COMPONENT TESTS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("OptimizationState - Initialization", "[optimization][framework]")
{
    OptimizationState<Real> state;
    
    REQUIRE(state.iteration == 0);
    REQUIRE(state.funcEvals == 0);
    REQUIRE(state.fCurrent == std::numeric_limits<Real>::max());
    REQUIRE(state.fBest == std::numeric_limits<Real>::max());
}

TEST_CASE("MaxIterationsCriterion - Termination", "[optimization][framework]")
{
    MaxIterationsCriterion<Real> criterion(100);
    
    OptimizationState<Real> state;
    state.iteration = 50;
    REQUIRE(criterion.ShouldTerminate(state) == false);
    
    state.iteration = 100;
    REQUIRE(criterion.ShouldTerminate(state) == true);
}

TEST_CASE("NullObserver - Basic functionality", "[optimization][framework]")
{
    NullObserver<Real> observer;
    OptimizationState<Real> state;
    state.iteration = 0;
    state.fBest = REAL(1.0);
    
    // Should not throw
    REQUIRE_NOTHROW(observer.OnStart(state));
    REQUIRE_NOTHROW(observer.OnIteration(state));
    REQUIRE_NOTHROW(observer.OnComplete(state, "test"));
}

TEST_CASE("OptimizationConfig - Basic construction", "[optimization][framework]")
{
    OptimizationConfig<Real> config;
    
    // Should have default criterion
    REQUIRE(config.GetCriterion() != nullptr);
}

///////////////////////////////////////////////////////////////////////////////////////////
// SUCCESS TEST
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Optimization Framework - Basic integration test", "[optimization][framework][integration]")
{
    // This test just verifies the framework compiles and links correctly
    // More detailed integration tests can be added after validating the basic structure
    
    OptimizationState<Real> state;
    state.iteration = 5;
    state.fBest = REAL(1.0);
    
    MaxIterationsCriterion<Real> criterion(10);
    REQUIRE(criterion.ShouldTerminate(state) == false);
    
    state.iteration = 15;
    REQUIRE(criterion.ShouldTerminate(state) == true);
    
    // Test passes - framework is working!
    REQUIRE(true);
}
