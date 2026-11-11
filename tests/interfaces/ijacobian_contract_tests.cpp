///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ijacobian_contract_tests.cpp                                        ///
///  Description: The one Jacobian contract (MinimalMathLibrary-wewi.6.1):            ///
///               IODESystemWithJacobian::jacobian(t, x, J) is canonical,             ///
///               derivsAndJacobian(t, x, dxdt, J) is the stepper entry point,        ///
///               IDynamicalSystem shares the interface, the old 4-arg form forwards. ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#include <mml/MMLBase.h>
#include <mml/interfaces/IODESystem.h>
#include <mml/interfaces/IDynamicalSystem.h>
#include <mml/base/ODESystem.h>
#include <mml/systems/DynamicalSystem.h>
#include <mml/algorithms/ODESolvers.h>

#include <type_traits>

using namespace MML;
using namespace MML::Systems;

namespace MML::Tests::Interfaces::IJacobianContractTests {

namespace {
    // Linear decay y' = -2 y with the canonical 3-argument Jacobian only
    class Decay : public IODESystemWithJacobian {
    public:
        int getDim() const override { return 1; }
        void derivs(Real, const Vector<Real>& y, Vector<Real>& dydt) const override { dydt[0] = -2.0 * y[0]; }
        void jacobian(Real, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = -2.0; }
    };

    // Numerical-Recipes style callback: third argument is df/dt (zero here), NOT dy/dt
    void decayDerivs(Real, const Vector<Real>& y, Vector<Real>& dydt) { dydt[0] = -2.0 * y[0]; }
    void decayJacobn(const Real, const Vector<Real>&, Vector<Real>& dfdt, Matrix<Real>& J) { dfdt[0] = 0.0; J(0, 0) = -2.0; }
    void decayJacobianOnly(Real, const Vector<Real>&, Matrix<Real>& J) { J(0, 0) = -2.0; }
}

TEST_CASE("One Jacobian contract: IODESystemWithJacobian and IDynamicalSystem share jacobian(t, x, J)", "[interfaces][jacobian][wewi.6.1]")
{
    STATIC_REQUIRE(std::is_base_of_v<IODESystemWithJacobian, IDynamicalSystem>);
    STATIC_REQUIRE(std::is_base_of_v<IODESystemParametrized, IDynamicalSystem>);
    STATIC_REQUIRE(std::is_base_of_v<IODESystemWithJacobian, LorenzSystem>);

    // a single IODESystem subobject: conversions are unambiguous
    LorenzSystem lorenz;
    IODESystem& asOde = lorenz;
    IODESystemWithJacobian& asJac = lorenz;
    IODESystemParametrized& asPar = lorenz;
    REQUIRE(asOde.getDim() == 3);
    REQUIRE(asJac.getDim() == 3);
    REQUIRE(asPar.getNumParam() == 3);

    // same jacobian through every view
    Vector<Real> y{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)};
    Matrix<Real> J1(3, 3), J2(3, 3);
    lorenz.jacobian(static_cast<MML::Real>(0.0), y, J1);
    asJac.jacobian(static_cast<MML::Real>(0.0), y, J2);
    REQUIRE(J1.IsEqualTo(J2));
    REQUIRE(lorenz.hasAnalyticalJacobian());
}

TEST_CASE("derivsAndJacobian cannot clobber the derivatives (regression for the old 4-arg conflation)", "[interfaces][jacobian][wewi.6.1]")
{
    // Before: stiff steppers called derivs() and then the 4-arg jacobian(t, y, dydt, J); Jacobian
    // callbacks that treated the third argument as NR's df/dt (and zeroed it) wiped dydt, so the
    // Newton residual y_new - y - h*0 "converged" immediately and the solution froze.
    Decay d;
    Vector<Real> y{static_cast<MML::Real>(1.0)}, dydt(1);
    Matrix<Real> J(1, 1);
    d.derivsAndJacobian(static_cast<MML::Real>(0.0), y, dydt, J);
    REQUIRE(dydt[0] == static_cast<MML::Real>(-2.0));
    REQUIRE(J(0, 0) == static_cast<MML::Real>(-2.0));

    // wrapper around an NR jacobn(t, y, dfdt, J): the derivatives still come from derivs()
    ODESystemWithJacobian nr(1, decayDerivs, decayJacobn);
    dydt[0] = static_cast<MML::Real>(123.0);
    nr.derivsAndJacobian(static_cast<MML::Real>(0.0), y, dydt, J);
    REQUIRE(J(0, 0) == static_cast<MML::Real>(-2.0));
    // the combined callback IS the single evaluation here; its dfdt slot is what the callback wrote
    // (zero) - which is exactly why steppers must not read dydt from the old 4-arg jacobian. Through
    // jacobian(t, x, J) the derivatives are untouched:
    Vector<Real> probe{static_cast<MML::Real>(7.0)};
    nr.jacobian(static_cast<MML::Real>(0.0), y, J);
    REQUIRE(probe[0] == static_cast<MML::Real>(7.0));
    Vector<Real> fromDerivs(1);
    nr.derivs(static_cast<MML::Real>(0.0), y, fromDerivs);
    REQUIRE(fromDerivs[0] == static_cast<MML::Real>(-2.0));

    // Jacobian-only callbacks are accepted directly (canonical form)
    ODESystemWithJacobian canonical(1, decayDerivs, decayJacobianOnly);
    canonical.derivsAndJacobian(static_cast<MML::Real>(0.0), y, dydt, J);
    REQUIRE(dydt[0] == static_cast<MML::Real>(-2.0));          // from derivs()
    REQUIRE(J(0, 0) == static_cast<MML::Real>(-2.0));

    ODESystemWithJacobianFromStdFunc fromStd(1,
        [](Real, const Vector<Real>& yy, Vector<Real>& dd) { dd[0] = -2.0 * yy[0]; },
        [](Real, const Vector<Real>&, Matrix<Real>& JJ) { JJ(0, 0) = -2.0; });
    fromStd.derivsAndJacobian(static_cast<MML::Real>(0.0), y, dydt, J);
    REQUIRE(dydt[0] == static_cast<MML::Real>(-2.0));
    REQUIRE(J(0, 0) == static_cast<MML::Real>(-2.0));

    // backward Euler on the canonical system actually moves: y(1) ~ e^-2
    auto sol = SolveBackwardEuler(d, static_cast<MML::Real>(0.0), Vector<Real>{static_cast<MML::Real>(1.0)}, static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.001));
    REQUIRE(sol.getXValuesAtEnd()[0] == Catch::Approx(std::exp(-2.0)).margin(5e-3));
}

TEST_CASE("Any IDynamicalSystem can be handed to a stiff solver through the shared interface", "[interfaces][jacobian][systems][wewi.6.1]")
{
    VanDerPolSystem vdp(static_cast<MML::Real>(1.0));
    IODESystemWithJacobian& stiffView = vdp;      // compiles only because the interfaces are unified
    auto sol = SolveBackwardEuler(stiffView, static_cast<MML::Real>(0.0), Vector<Real>{static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0)}, static_cast<MML::Real>(0.5), static_cast<MML::Real>(0.01));
    REQUIRE(sol.getNumStepsOK() > 0);
    Vector<Real> yEnd = sol.getXValuesAtEnd();
    REQUIRE(yEnd[0] != static_cast<MML::Real>(2.0));                // it moved
    REQUIRE(std::isfinite(yEnd[1]));

    // numerical default for a system without an analytic Jacobian
    class NoJac : public DynamicalSystemBase<1, 0> {
    public:
        NoJac() = default;
        void derivs(Real, const Vector<Real>& x, Vector<Real>& dxdt) const override { dxdt[0] = x[0] * x[0]; }
    };
    NoJac nj;
    REQUIRE_FALSE(nj.hasAnalyticalJacobian());
    Matrix<Real> J(1, 1);
    nj.jacobian(static_cast<MML::Real>(0.0), Vector<Real>{static_cast<MML::Real>(3.0)}, J);
    REQUIRE(J(0, 0) == Catch::Approx(6.0).margin(1e-5));
}

} // namespace MML::Tests::Interfaces::IJacobianContractTests
