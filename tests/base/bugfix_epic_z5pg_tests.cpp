///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        bugfix_epic_z5pg_tests.cpp                                          ///
///  Description: Regression tests for bugs fixed in epic z5pg                        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include <cmath>
#include <limits>
#include <string>
#include <stdexcept>

#include <mml/MMLBase.h>

#include <mml/base/Vector/VectorN.h>
#include <mml/base/Function.h>
#include <mml/base/Graph.h>
#include <mml/base/InterpolatedFunctions/InterpolatedRealFunctionRational.h>
#include <mml/interfaces/IFunction.h>
#include <mml/core/Integration/IntegrationBase.h>
#include <mml/core/OrthogonalBasis/ChebyshevBasis.h>

using namespace MML;

namespace MML::Tests::BugfixZ5pg {

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug .33 — VectorN::Normalized: near-zero norm should throw, not produce Inf      ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_z5pg.33 - VectorN::Normalized near-zero vector throws", "[bugfix][z5pg][vector]")
{
	// A vector with components much smaller than epsilon should not be normalized
	Real tiny = std::numeric_limits<Real>::epsilon() * 1e-3;
	VectorN<Real, 3> v({tiny, tiny, tiny});

	REQUIRE_THROWS(v.Normalized());
}

TEST_CASE("BugFix_z5pg.33 - VectorN::Normalized normal vector succeeds", "[bugfix][z5pg][vector]")
{
	VectorN<Real, 3> v({3.0, 4.0, 0.0});
	VectorN<Real, 3> n = v.Normalized();
	REQUIRE(n.NormL2() == Catch::Approx(1.0));
}

TEST_CASE("BugFix_z5pg.33 - VectorN::Normalized exactly zero throws", "[bugfix][z5pg][vector]")
{
	VectorN<Real, 3> v({0.0, 0.0, 0.0});
	REQUIRE_THROWS(v.Normalized());
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug .34 — VectorN near-zero detection (superseded by CONTRACT_POLICY:             ///
///  isZero() is now exact; tolerance-based check is isNearZero())                     ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_z5pg.34 - VectorN::isNearZero with near-zero values", "[bugfix][z5pg][vector]")
{
	Real tiny = std::numeric_limits<Real>::epsilon() * 0.1;
	VectorN<Real, 3> v({tiny, -tiny, tiny});

	// Near-zero values smaller than tolerance should be considered near-zero
	REQUIRE(v.isNearZero());
	REQUIRE_FALSE(v.isZero());   // contract policy: isZero is exact
}

TEST_CASE("BugFix_z5pg.34 - VectorN::isZero with non-zero values", "[bugfix][z5pg][vector]")
{
	VectorN<Real, 3> v({0.0, 0.0, 1.0});
	REQUIRE_FALSE(v.isZero());
}

TEST_CASE("BugFix_z5pg.34 - VectorN::isZero exact zero", "[bugfix][z5pg][vector]")
{
	VectorN<Real, 3> v({0.0, 0.0, 0.0});
	REQUIRE(v.isZero());
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug .35 — IFunction::GetValues: numPnt < 2 should throw, not divide by zero       ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_z5pg.35 - GetValues with numPnt=1 throws", "[bugfix][z5pg][function]")
{
	RealFunctionFromStdFunc f([](Real x) { return x * x; });
	Vector<Real> outX, outY;
	
	REQUIRE_THROWS_AS(f.GetValues(0.0, 1.0, 1, outX, outY), std::invalid_argument);
}

TEST_CASE("BugFix_z5pg.35 - GetValues with numPnt=0 throws", "[bugfix][z5pg][function]")
{
	RealFunctionFromStdFunc f([](Real x) { return x * x; });
	Vector<Real> outX, outY;
	
	REQUIRE_THROWS_AS(f.GetValues(0.0, 1.0, 0, outX, outY), std::invalid_argument);
}

TEST_CASE("BugFix_z5pg.35 - GetValues with numPnt=2 succeeds", "[bugfix][z5pg][function]")
{
	RealFunctionFromStdFunc f([](Real x) { return x * x; });
	Vector<Real> outX, outY;
	
	f.GetValues(0.0, 1.0, 2, outX, outY);
	REQUIRE(outX.size() == 2);
	REQUIRE(outY.size() == 2);
	REQUIRE(outX[0] == Catch::Approx(0.0));
	REQUIRE(outX[1] == Catch::Approx(1.0));
	REQUIRE(outY[0] == Catch::Approx(0.0));
	REQUIRE(outY[1] == Catch::Approx(1.0));
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug .38 — Rational interpolation: division by zero when x coincides with node     ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_z5pg.38 - Rational interpolation at exact node returns exact value", "[bugfix][z5pg][interp]")
{
	// Create data points: x = 0, 1, 2, 3, 4
	Vector<Real> xData(5), yData(5);
	for (int i = 0; i < 5; i++) {
		xData[i] = static_cast<Real>(i);
		yData[i] = std::sin(static_cast<Real>(i));
	}
	
	RationalInterpRealFunc interp(xData, yData, 3);
	
	// Evaluating exactly at a node should return exact value (not NaN/Inf)
	Real val = interp(xData[2]);
	REQUIRE(val == Catch::Approx(yData[2]));
}

TEST_CASE("BugFix_z5pg.38 - Rational interpolation near a node returns its value", "[bugfix][z5pg][interp]")
{
	Vector<Real> xData{0.0, 1.0, 2.0, 3.0, 4.0};
	Vector<Real> yData{0.0, 1.0, 4.0, 9.0, 16.0};
	RationalInterpRealFunc interp(xData, yData, 3);

	const Real nearNode = std::nextafter(xData[2], xData[1]);
	REQUIRE(interp(nearNode) == yData[2]);
}

TEST_CASE("BugFix_z5pg.38 - Rational interpolation between nodes works", "[bugfix][z5pg][interp]")
{
	Vector<Real> xData(5), yData(5);
	for (int i = 0; i < 5; i++) {
		xData[i] = static_cast<Real>(i);
		yData[i] = std::sin(static_cast<Real>(i));
	}
	
	RationalInterpRealFunc interp(xData, yData, 3);
	
	// Interpolation at midpoint should be finite and reasonable
	Real val = interp(1.5);
	REQUIRE(std::isfinite(val));
	// Rational interpolation with m=3 on sparse sin data — just check it's in ballpark
	REQUIRE(val == Catch::Approx(std::sin(1.5)).margin(0.2));
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug .39 — Graph::removeEdge: edge count decremented before reverse edge removal   ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_z5pg.39 - Graph undirected removeEdge edge count", "[bugfix][z5pg][graph]")
{
	Graph<> g(3);
	
	g.addEdge(0, 1, 1.0);
	g.addEdge(1, 2, 2.0);
	
	REQUIRE(g.numEdges() == 2);
	
	// Remove one edge
	bool removed = g.removeEdge(0, 1);
	REQUIRE(removed);
	REQUIRE(g.numEdges() == 1);
	
	// The reverse edge should also be removed (undirected)
	// Verify by checking adjacency
	const auto& neighbors0 = g.neighbors(0);
	REQUIRE(neighbors0.empty());
	
	const auto& neighbors1 = g.neighbors(1);
	// Node 1 should only have edge to node 2 now
	REQUIRE(neighbors1.size() == 1);
}

TEST_CASE("BugFix_z5pg.39 - Graph directed removeEdge edge count", "[bugfix][z5pg][graph]")
{
	Graph<> g(2, Graph<>::Type::Directed);
	
	g.addEdge(0, 1, 1.0);
	
	REQUIRE(g.numEdges() == 1);
	
	bool removed = g.removeEdge(0, 1);
	REQUIRE(removed);
	REQUIRE(g.numEdges() == 0);
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug .44 — Chebyshev weight function returns 1e10 instead of infinity at boundary  ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_z5pg.44 - Chebyshev weight function at boundary", "[bugfix][z5pg][chebyshev]")
{
	ChebyshevBasis basis;
	
	// At x = ±1, the weight function 1/sqrt(1-x²) should be infinity
	Real w_at_1 = basis.WeightFunction(1.0);
	Real w_at_neg1 = basis.WeightFunction(-1.0);
	
	REQUIRE(std::isinf(w_at_1));
	REQUIRE(std::isinf(w_at_neg1));
	
	// For interior points, should return finite values
	Real w_at_0 = basis.WeightFunction(0.0);
	REQUIRE(std::isfinite(w_at_0));
	REQUIRE(w_at_0 == Catch::Approx(1.0));  // 1/sqrt(1-0) = 1
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug .45 — IntegrationResult implicit conversion should work but .value preferred  ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BugFix_z5pg.45 - IntegrationResult value access", "[bugfix][z5pg][integration]")
{
	IntegrationResult res(3.14, 1e-10, 5, true);
	
	// Direct .value access (preferred)
	REQUIRE(res.value == Catch::Approx(3.14));
	REQUIRE(res.error_estimate == Catch::Approx(1e-10));
	REQUIRE(res.iterations == 5);
	REQUIRE(res.converged);
	
	// Implicit conversion to Real still works
	Real val = res;
	REQUIRE(val == Catch::Approx(3.14));
}

///////////////////////////////////////////////////////////////////////////////////////////
///  Bug .49 — StatisticsRank z-score fix: regression test moved to MML-Packages          ///
///  tests/mml_ext/rank_correlation_tests.cpp with the 2.0 cull (MinimalMathLibrary-ya0v.3) ///
///////////////////////////////////////////////////////////////////////////////////////////

} // namespace MML::Tests::BugfixZ5pg
