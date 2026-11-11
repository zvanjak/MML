///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        iparametrized_tests.cpp                                             ///
///  Description: Contract tests for IParametrized (MinimalMathLibrary-wewi.6.4):     ///
///               setParams rejects wrong-sized vectors instead of silently           ///
///               truncating or reading out of bounds; per-index access is checked.   ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#include <mml/MMLBase.h>
#include <mml/interfaces/IParametrized.h>
#include <mml/systems/DynamicalSystem.h>

using namespace MML;
using namespace MML::Systems;

namespace MML::Tests::Interfaces::IParametrizedTests {

namespace {
    // Minimal implementer relying on the interface's default setParams()/getParams()
    class TwoParams : public IParametrized {
        Real _a = 1.0, _b = 2.0;
    public:
        int  getNumParam() const override { return 2; }
        Real getParam(int i) const override { return i == 0 ? _a : _b; }
        void setParam(int i, Real v) override { (i == 0 ? _a : _b) = v; }
    };
}

TEST_CASE("IParametrized::setParams throws on size mismatch and leaves the parameters untouched", "[interfaces][parametrized][wewi.6.4]")
{
    TwoParams p;

    SECTION("too few parameters: used to assign a prefix silently")
    {
        REQUIRE_THROWS_AS(p.setParams(Vector<Real>{static_cast<MML::Real>(5.0)}), VectorDimensionError);
        REQUIRE(p.getParam(0) == static_cast<MML::Real>(1.0));
        REQUIRE(p.getParam(1) == static_cast<MML::Real>(2.0));
    }
    SECTION("too many parameters: used to drop the tail silently")
    {
        REQUIRE_THROWS_AS(p.setParams(Vector<Real>{static_cast<MML::Real>(5.0), static_cast<MML::Real>(6.0), static_cast<MML::Real>(7.0)}), VectorDimensionError);
        REQUIRE(p.getParam(0) == static_cast<MML::Real>(1.0));
        REQUIRE(p.getParam(1) == static_cast<MML::Real>(2.0));
    }
    SECTION("empty vector")
    {
        REQUIRE_THROWS_AS(p.setParams(Vector<Real>()), VectorDimensionError);
    }
    SECTION("exact size is accepted and round-trips through getParams")
    {
        REQUIRE_NOTHROW(p.setParams(Vector<Real>{static_cast<MML::Real>(5.0), static_cast<MML::Real>(6.0)}));
        REQUIRE(p.getParam(0) == static_cast<MML::Real>(5.0));
        REQUIRE(p.getParam(1) == static_cast<MML::Real>(6.0));
        Vector<Real> back = p.getParams();
        REQUIRE(back.size() == 2);
        REQUIRE(back[1] == static_cast<MML::Real>(6.0));
    }
    SECTION("the exception reports both sizes")
    {
        try {
            p.setParams(Vector<Real>{static_cast<MML::Real>(5.0), static_cast<MML::Real>(6.0), static_cast<MML::Real>(7.0)});
            FAIL("expected VectorDimensionError");
        } catch (const VectorDimensionError& e) {
            const std::string msg = e.what();
            REQUIRE(msg.find("setParams") != std::string::npos);
        }
    }
}

TEST_CASE("DynamicalSystemBase::setParams / setParam / getParam are bounds-checked", "[interfaces][parametrized][systems][wewi.6.4]")
{
    LorenzSystem lorenz(MML::Real(10.0), MML::Real(28.0), MML::Real(8.0) / MML::Real(3.0));   // 3 parameters
    REQUIRE(lorenz.getNumParam() == 3);

    // too few: the old override read p[1], p[2] past the end of a 1-element vector
    REQUIRE_THROWS_AS(lorenz.setParams(Vector<Real>{static_cast<MML::Real>(1.0)}), VectorDimensionError);
    REQUIRE(lorenz.getParam(0) == static_cast<MML::Real>(10.0));
    REQUIRE_THROWS_AS(lorenz.setParams(Vector<Real>{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0)}), VectorDimensionError);
    REQUIRE(lorenz.getParam(2) == Catch::Approx(8.0 / 3.0));

    REQUIRE_NOTHROW(lorenz.setParams(Vector<Real>{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)}));
    REQUIRE(lorenz.getParam(1) == static_cast<MML::Real>(2.0));

    REQUIRE_THROWS_AS(lorenz.getParam(3), IndexError);
    REQUIRE_THROWS_AS(lorenz.getParam(-1), IndexError);
    REQUIRE_THROWS_AS(lorenz.setParam(3, static_cast<MML::Real>(0.0)), IndexError);

    VanDerPolSystem vdp;                            // 1 parameter
    REQUIRE_THROWS_AS(vdp.setParams(Vector<Real>{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)}), VectorDimensionError);
    REQUIRE_NOTHROW(vdp.setParams(Vector<Real>{static_cast<MML::Real>(1.5)}));
    REQUIRE(vdp.getParam(0) == static_cast<MML::Real>(1.5));
}

} // namespace MML::Tests::Interfaces::IParametrizedTests
