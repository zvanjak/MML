///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        idiscrete_map_tests.cpp                                             ///
///  Description: Contract tests for IDiscreteMap after it joined the interface       ///
///               family (MinimalMathLibrary-wewi.6.2): namespace MML, IFunction +    ///
///               IParametrized bases, checked iterate(), parameter metadata defaults.///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#include <mml/MMLBase.h>
#include <mml/interfaces/IDiscreteMap.h>
#include <mml/systems/DiscreteMaps.h>

#include <type_traits>

using namespace MML;

namespace MML::Tests::Interfaces::IDiscreteMapTests {

namespace {
    // Minimal implementer that relies on every default the interface provides
    class Doubling : public IDiscreteMap<1> {
        Real _k = 2.0;
    public:
        void map(const Vector<Real>& x, Vector<Real>& xNext) const override { xNext[0] = _k * x[0]; }
        void jacobian(const Vector<Real>&, Matrix<Real>& J) const override { J.resize(1, 1); J(0, 0) = _k; }
        int  getNumParam() const override { return 1; }
        Real getParam(int) const override { return _k; }
        void setParam(int, Real v) override { _k = v; }
    };
}

TEST_CASE("IDiscreteMap lives in namespace MML and is an IFunction + IParametrized", "[interfaces][discrete-map][wewi.6.2]")
{
    STATIC_REQUIRE(std::is_base_of_v<IParametrized, IDiscreteMap<2>>);
    STATIC_REQUIRE((std::is_base_of_v<IFunction<Vector<Real>, const Vector<Real>&>, IDiscreteMap<2>>));
    STATIC_REQUIRE(IDiscreteMap<2>::Dimension == 2);
    // the old spelling still names the same type
    STATIC_REQUIRE(std::is_same_v<Systems::IDiscreteMap<2>, IDiscreteMap<2>>);
    STATIC_REQUIRE(std::is_base_of_v<IDiscreteMap<2>, Systems::HenonMap>);
    STATIC_REQUIRE(std::is_base_of_v<IDiscreteMap<1>, Systems::LogisticMap>);
}

TEST_CASE("IDiscreteMap: operator() == iterate(), dimension is checked, defaults for names", "[interfaces][discrete-map][wewi.6.2]")
{
    Doubling d;
    Vector<Real> x{static_cast<MML::Real>(1.5)};
    REQUIRE(d.iterate(x)[0] == static_cast<MML::Real>(3.0));
    REQUIRE(d(x)[0] == static_cast<MML::Real>(3.0));
    const IFunction<Vector<Real>, const Vector<Real>&>& f = d;
    REQUIRE(f(x)[0] == static_cast<MML::Real>(3.0));
    REQUIRE(d.getDim() == 1);

    REQUIRE_THROWS_AS(d.iterate(Vector<Real>{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)}), VectorDimensionError);
    REQUIRE_THROWS_AS(d.iterate(Vector<Real>()), VectorDimensionError);

    REQUIRE(d.getParamName(0) == "p0");
    REQUIRE(d.getStateName(0) == "x0");
}

TEST_CASE("Concrete maps use the IParametrized mixin: getNumParam, bulk get/set, checked indices", "[interfaces][discrete-map][systems][wewi.6.2]")
{
    Systems::HenonMap henon(1.4, 0.3);
    REQUIRE(henon.getNumParam() == 2);
    Vector<Real> p = henon.getParams();
    REQUIRE(p.size() == 2);
    REQUIRE(p[0] == Catch::Approx(1.4));
    REQUIRE(p[1] == Catch::Approx(0.3));
    REQUIRE_NOTHROW(henon.setParams(Vector<Real>{static_cast<MML::Real>(1.2), static_cast<MML::Real>(0.2)}));
    REQUIRE(henon.getParam(0) == Catch::Approx(1.2));
    REQUIRE_THROWS_AS(henon.setParams(Vector<Real>{static_cast<MML::Real>(1.0)}), VectorDimensionError);     // IParametrized contract (wewi.6.4)
    REQUIRE_THROWS_AS(henon.getParam(2), IndexError);                                        // used to return _b silently
    REQUIRE_THROWS_AS(henon.setParam(-1, static_cast<MML::Real>(0.0)), IndexError);
    REQUIRE(henon.getParamName(0) == "a");
    REQUIRE(henon.getStateName(1) == "y");

    Systems::LogisticMap logistic(3.9);
    REQUIRE(logistic.getNumParam() == 1);
    REQUIRE_THROWS_AS(logistic.getParam(1), IndexError);                                     // used to return 0 silently
    REQUIRE(logistic.getParams()[0] == Catch::Approx(3.9));
    // IFunction form on a library map
    Vector<Real> x0{static_cast<MML::Real>(0.5)};
    REQUIRE(logistic(x0)[0] == Catch::Approx(3.9 * 0.5 * 0.5));

    Systems::StandardMap standard(1.0);
    Systems::TentMap tent(2.0);
    REQUIRE(standard.getNumParam() == 1);
    REQUIRE(tent.getNumParam() == 1);
    REQUIRE_THROWS_AS(standard.setParam(1, static_cast<MML::Real>(2.0)), IndexError);                        // used to overwrite K for any index
    REQUIRE_THROWS_AS(tent.setParam(1, static_cast<MML::Real>(2.0)), IndexError);
}

} // namespace MML::Tests::Interfaces::IDiscreteMapTests
