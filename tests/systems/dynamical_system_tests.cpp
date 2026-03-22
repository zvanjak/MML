///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML) Tests                            ///
///                                                                                   ///
///  File:        dynamical_system_tests.cpp                                          ///
///  Description: Tests for DynamicalSystem framework                                 ///
///               Covers Lorenz, Rössler, Van der Pol, fixed points, Lyapunov         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "systems/DynamicalSystem.h"

using namespace MML;
using namespace MML::Systems;
using Catch::Approx;

namespace MML::Tests::Systems::DynamicalSystemTests {

///////////////////////////////////////////////////////////////////////////////////////////
//                              LORENZ SYSTEM TESTS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("LorenzSystem - Basic construction and parameters", "[DynamicalSystem][Lorenz]")
{
    LorenzSystem lorenz;
    
    REQUIRE(lorenz.getDim() == 3);
    REQUIRE(lorenz.getNumParam() == 3);
    
    // Default parameters: σ=10, ρ=28, β=8/3
    REQUIRE(lorenz.getParam(0) == Approx(10.0));
    REQUIRE(lorenz.getParam(1) == Approx(28.0));
    REQUIRE(lorenz.getParam(2) == Approx(8.0/3.0));
    
    // Metadata
    REQUIRE(lorenz.getStateName(0) == "x");
    REQUIRE(lorenz.getStateName(1) == "y");
    REQUIRE(lorenz.getStateName(2) == "z");
    
    REQUIRE(lorenz.getParamName(0) == "sigma");
    REQUIRE(lorenz.getParamName(1) == "rho");
    REQUIRE(lorenz.getParamName(2) == "beta");
}

TEST_CASE("LorenzSystem - Derivatives computation", "[DynamicalSystem][Lorenz]")
{
    LorenzSystem lorenz(10.0, 28.0, 8.0/3.0);
    
    Vector<Real> y({1.0, 2.0, 3.0});
    Vector<Real> dydt(3);
    lorenz.derivs(0.0, y, dydt);
    
    // dx/dt = σ(y - x) = 10(2 - 1) = 10
    REQUIRE(dydt[0] == Approx(10.0));
    
    // dy/dt = x(ρ - z) - y = 1(28 - 3) - 2 = 23
    REQUIRE(dydt[1] == Approx(23.0));
    
    // dz/dt = xy - βz = 1*2 - (8/3)*3 = 2 - 8 = -6
    REQUIRE(dydt[2] == Approx(-6.0));
}

TEST_CASE("LorenzSystem - Analytical Jacobian", "[DynamicalSystem][Lorenz]")
{
    LorenzSystem lorenz;
    REQUIRE(lorenz.hasAnalyticalJacobian() == true);
    
    Vector<Real> y({1.0, 2.0, 3.0});
    Matrix<Real> J;
    lorenz.jacobian(0.0, y, J);
    
    // J = [[-σ, σ, 0], [ρ-z, -1, -x], [y, x, -β]]
    REQUIRE(J(0,0) == Approx(-10.0));   // -σ
    REQUIRE(J(0,1) == Approx(10.0));    // σ
    REQUIRE(J(0,2) == Approx(0.0));
    
    REQUIRE(J(1,0) == Approx(25.0));    // ρ - z = 28 - 3
    REQUIRE(J(1,1) == Approx(-1.0));
    REQUIRE(J(1,2) == Approx(-1.0));    // -x
    
    REQUIRE(J(2,0) == Approx(2.0));     // y
    REQUIRE(J(2,1) == Approx(1.0));     // x
    REQUIRE(J(2,2) == Approx(-8.0/3.0)); // -β
}

TEST_CASE("LorenzSystem - Properties", "[DynamicalSystem][Lorenz]")
{
    LorenzSystem lorenz;
    
    REQUIRE(lorenz.isAutonomous() == true);
    REQUIRE(lorenz.isDissipative() == true);
    
    // Divergence = -(σ + 1 + β) = -(10 + 1 + 8/3) ≈ -13.67
    REQUIRE(lorenz.getDivergence() == Approx(-(10.0 + 1.0 + 8.0/3.0)));
}

///////////////////////////////////////////////////////////////////////////////////////////
//                              ROSSLER SYSTEM TESTS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("RosslerSystem - Basic construction", "[DynamicalSystem][Rossler]")
{
    RosslerSystem rossler;
    
    REQUIRE(rossler.getDim() == 3);
    REQUIRE(rossler.getNumParam() == 3);
    
    // Default: a=0.2, b=0.2, c=5.7
    REQUIRE(rossler.getParam(0) == Approx(0.2));
    REQUIRE(rossler.getParam(1) == Approx(0.2));
    REQUIRE(rossler.getParam(2) == Approx(5.7));
}

TEST_CASE("RosslerSystem - Derivatives", "[DynamicalSystem][Rossler]")
{
    RosslerSystem rossler(0.2, 0.2, 5.7);
    
    Vector<Real> y({1.0, 2.0, 3.0});
    Vector<Real> dydt(3);
    rossler.derivs(0.0, y, dydt);
    
    // dx/dt = -y - z = -2 - 3 = -5
    REQUIRE(dydt[0] == Approx(-5.0));
    
    // dy/dt = x + a*y = 1 + 0.2*2 = 1.4
    REQUIRE(dydt[1] == Approx(1.4));
    
    // dz/dt = b + z(x - c) = 0.2 + 3(1 - 5.7) = 0.2 - 14.1 = -13.9
    REQUIRE(dydt[2] == Approx(-13.9));
}

///////////////////////////////////////////////////////////////////////////////////////////
//                              VAN DER POL TESTS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("VanDerPolSystem - Basic construction", "[DynamicalSystem][VanDerPol]")
{
    VanDerPolSystem vdp;
    
    REQUIRE(vdp.getDim() == 2);
    REQUIRE(vdp.getNumParam() == 1);
    REQUIRE(vdp.getParam(0) == Approx(1.0));  // Default μ
}

TEST_CASE("VanDerPolSystem - Fixed point at origin", "[DynamicalSystem][VanDerPol]")
{
    VanDerPolSystem vdp(1.0);
    
    // Origin is a fixed point
    Vector<Real> origin({0.0, 0.0});
    Vector<Real> dydt(2);
    vdp.derivs(0.0, origin, dydt);
    
    REQUIRE(std::abs(dydt[0]) < 1e-14);
    REQUIRE(std::abs(dydt[1]) < 1e-14);
}

TEST_CASE("VanDerPolSystem - Jacobian at origin", "[DynamicalSystem][VanDerPol]")
{
    VanDerPolSystem vdp(2.0);  // μ = 2
    
    Vector<Real> origin({0.0, 0.0});
    Matrix<Real> J;
    vdp.jacobian(0.0, origin, J);
    
    // At origin: J = [[0, 1], [-1, μ]]
    REQUIRE(J(0,0) == Approx(0.0));
    REQUIRE(J(0,1) == Approx(1.0));
    REQUIRE(J(1,0) == Approx(-1.0));
    REQUIRE(J(1,1) == Approx(2.0));  // μ(1 - 0) = μ
}

///////////////////////////////////////////////////////////////////////////////////////////
//                              DUFFING SYSTEM TESTS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("DuffingSystem - Construction", "[DynamicalSystem][Duffing]")
{
    DuffingSystem duffing;
    
    REQUIRE(duffing.getDim() == 3);
    REQUIRE(duffing.getNumParam() == 5);
    REQUIRE(duffing.isAutonomous() == false);  // Has driving term
}

///////////////////////////////////////////////////////////////////////////////////////////
//                              FIXED POINT FINDER TESTS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FixedPointFinder - Lorenz origin", "[DynamicalSystem][FixedPoint]")
{
    LorenzSystem lorenz(10.0, 0.5, 8.0/3.0);  // ρ < 1 → only origin is stable
    
    Vector<Real> guess({0.1, 0.1, 0.1});
    auto fp = FixedPointFinder::Find(lorenz, guess);
    
    // Should find origin
    REQUIRE(fp.convergenceResidual < 1e-8);
    REQUIRE(fp.location[0] == Approx(0.0).margin(1e-8));
    REQUIRE(fp.location[1] == Approx(0.0).margin(1e-8));
    REQUIRE(fp.location[2] == Approx(0.0).margin(1e-8));
}

TEST_CASE("FixedPointFinder - Van der Pol origin", "[DynamicalSystem][FixedPoint]")
{
    VanDerPolSystem vdp(1.0);
    
    Vector<Real> guess({0.1, 0.1});
    auto fp = FixedPointFinder::Find(vdp, guess);
    
    REQUIRE(fp.convergenceResidual < 1e-8);
    REQUIRE(fp.location[0] == Approx(0.0).margin(1e-8));
    REQUIRE(fp.location[1] == Approx(0.0).margin(1e-8));
    
    // Origin is unstable for μ > 0
    REQUIRE(fp.isStable == false);
    REQUIRE(fp.type == FixedPointType::UnstableFocus);
}

TEST_CASE("FixedPointFinder - Classification", "[DynamicalSystem][FixedPoint]")
{
    // Stable node: all eigenvalues negative real
    // Create a simple 2D system with stable node at origin
    class StableNodeSystem : public DynamicalSystemBase<2, 0>
    {
    public:
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            // dx/dt = -x, dy/dt = -2y → eigenvalues -1, -2
            dydt[0] = -y[0];
            dydt[1] = -2.0 * y[1];
        }
    };
    
    StableNodeSystem sys;
    auto fp = FixedPointFinder::Find(sys, Vector<Real>({0.1, 0.1}));
    
    REQUIRE(fp.type == FixedPointType::StableNode);
    REQUIRE(fp.isStable == true);
}

TEST_CASE("FixedPointFinder - Saddle classification", "[DynamicalSystem][FixedPoint]")
{
    class SaddleSystem : public DynamicalSystemBase<2, 0>
    {
    public:
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            // dx/dt = x, dy/dt = -y → eigenvalues +1, -1
            dydt[0] = y[0];
            dydt[1] = -y[1];
        }
    };
    
    SaddleSystem sys;
    auto fp = FixedPointFinder::Find(sys, Vector<Real>({0.0, 0.0}));
    
    REQUIRE(fp.type == FixedPointType::Saddle);
    REQUIRE(fp.isStable == false);
}

///////////////////////////////////////////////////////////////////////////////////////////
//                              LYAPUNOV ANALYZER TESTS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("LyapunovAnalyzer - Simple stable system", "[DynamicalSystem][Lyapunov]")
{
    class DecaySystem : public DynamicalSystemBase<2, 0>
    {
    public:
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            dydt[0] = -y[0];
            dydt[1] = -2.0 * y[1];
        }
        
        bool hasAnalyticalJacobian() const override { return true; }
        
        void jacobian(Real t, const Vector<Real>& y, Matrix<Real>& J) const override
        {
            J.Resize(2, 2);
            J(0,0) = -1; J(0,1) = 0;
            J(1,0) = 0;  J(1,1) = -2;
        }
    };
    
    DecaySystem sys;
    Vector<Real> x0({1.0, 1.0});
    
    auto result = LyapunovAnalyzer::Compute(sys, x0, 50.0, 1.0, 0.01);
    
    // Lyapunov exponents should be approximately -1 and -2
    REQUIRE(result.exponents[0] == Approx(-1.0).margin(0.1));
    REQUIRE(result.exponents[1] == Approx(-2.0).margin(0.1));
    
    REQUIRE(result.isChaotic == false);
    REQUIRE(result.maxExponent < 0);
}

TEST_CASE("LyapunovAnalyzer - Lorenz chaotic regime", "[DynamicalSystem][Lyapunov][Slow]")
{
    LorenzSystem lorenz;  // Default: σ=10, ρ=28, β=8/3
    Vector<Real> x0 = lorenz.getDefaultInitialCondition();
    
    // Short integration for test speed (longer = more accurate)
    auto result = LyapunovAnalyzer::Compute(lorenz, x0, 100.0, 1.0, 0.01);
    
    // Lorenz should have one positive exponent (≈0.9)
    REQUIRE(result.isChaotic == true);
    REQUIRE(result.maxExponent > 0);
    
    // Sum of exponents should be negative (dissipative)
    REQUIRE(result.sum < 0);
    
    // Kaplan-Yorke dimension for Lorenz ≈ 2.06
    REQUIRE(result.kaplanYorkeDimension > 2.0);
    REQUIRE(result.kaplanYorkeDimension < 2.5);
}

///////////////////////////////////////////////////////////////////////////////////////////
//                              PHASE SPACE ANALYZER TESTS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("PhaseSpaceAnalyzer - Trajectory integration", "[DynamicalSystem][PhaseSpace]")
{
    LorenzSystem lorenz;
    Vector<Real> x0({1.0, 1.0, 1.0});
    
    auto traj = PhaseSpaceAnalyzer::IntegrateTrajectory(lorenz, x0, 10.0, 0.1, 0.01);
    
    // Should have approximately 100 points (10.0 / 0.1)
    REQUIRE(traj.size() >= 100);
    
    // Trajectory should stay bounded (characteristic of Lorenz)
    for (const auto& pt : traj)
    {
        REQUIRE(std::abs(pt[0]) < 50);
        REQUIRE(std::abs(pt[1]) < 50);
        REQUIRE(pt[2] > 0);  // z always positive after initial transient
        REQUIRE(pt[2] < 60);
    }
}

TEST_CASE("PhaseSpaceAnalyzer - Poincare section", "[DynamicalSystem][PhaseSpace]")
{
    LorenzSystem lorenz;
    Vector<Real> x0({1.0, 1.0, 1.0});
    
    // Poincaré section at z = 27 (near attractor center)
    PoincareSection<Real> section(2, 27.0, 1);  // z = 27, positive crossing
    
    auto intersections = PhaseSpaceAnalyzer::ComputePoincareSection(
        lorenz, x0, section, 50, 0.001);
    
    REQUIRE(intersections.size() == 50);
    
    // All z values should be approximately 27
    for (const auto& pt : intersections)
    {
        REQUIRE(pt[2] == Approx(27.0).margin(0.1));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
//                              BIFURCATION ANALYZER TESTS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("BifurcationAnalyzer - Lorenz rho sweep", "[DynamicalSystem][Bifurcation][Slow]")
{
    LorenzSystem lorenz;
    Vector<Real> x0({1.0, 1.0, 1.0});
    
    // Sweep ρ from 20 to 30
    auto diagram = BifurcationAnalyzer::Sweep(
        lorenz,
        1,            // param index (rho)
        20.0, 30.0,   // range
        5,            // steps (small for test)
        x0,
        2,            // component to record (z)
        50.0,         // transient
        20.0,         // record time
        0.01);
    
    REQUIRE(diagram.parameterValues.size() == 5);
    REQUIRE(diagram.attractorValues.size() == 5);
    REQUIRE(diagram.parameterName == "rho");
    
    // At ρ=28, should see chaotic attractor with multiple maxima
    int chaotic_index = 3;  // ρ ≈ 27.5
    REQUIRE(diagram.attractorValues[chaotic_index].size() > 1);  // Multiple maxima = chaos
}

///////////////////////////////////////////////////////////////////////////////////////////
//                              NUMERICAL JACOBIAN TESTS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Numerical Jacobian vs Analytical", "[DynamicalSystem][Jacobian]")
{
    LorenzSystem lorenz;
    Vector<Real> y({5.0, 10.0, 20.0});
    
    Matrix<Real> J_analytical;
    lorenz.jacobian(0.0, y, J_analytical);
    
    // Create wrapper that uses numerical Jacobian
    class LorenzNumerical : public LorenzSystem
    {
    public:
        bool hasAnalyticalJacobian() const override { return false; }
        // Uses base class numerical jacobian
    };
    
    LorenzNumerical lorenzNum;
    lorenzNum.setParams(lorenz.getParams());
    
    Matrix<Real> J_numerical;
    // Call parent's numerical jacobian
    static_cast<IDynamicalSystem&>(lorenzNum).jacobian(0.0, y, J_numerical);
    
    // Should match to good precision
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            REQUIRE(J_numerical(i,j) == Approx(J_analytical(i,j)).margin(1e-6));
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
//                              FIXED POINT TYPE TO STRING
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FixedPointType ToString", "[DynamicalSystem][Utilities]")
{
    REQUIRE(ToString(FixedPointType::StableNode) == "Stable Node");
    REQUIRE(ToString(FixedPointType::UnstableNode) == "Unstable Node");
    REQUIRE(ToString(FixedPointType::Saddle) == "Saddle");
    REQUIRE(ToString(FixedPointType::StableFocus) == "Stable Focus");
    REQUIRE(ToString(FixedPointType::UnstableFocus) == "Unstable Focus");
    REQUIRE(ToString(FixedPointType::Center) == "Center");
    REQUIRE(ToString(FixedPointType::Unknown) == "Unknown");
}

///////////////////////////////////////////////////////////////////////////////////////////
//                           ADDITIONAL CONTINUOUS SYSTEMS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ChuaCircuit - Basic construction", "[DynamicalSystem][ChuaCircuit]")
{
    ChuaCircuit chua;
    
    REQUIRE(chua.getDim() == 3);
    REQUIRE(chua.getStateName(0) == "x");
    REQUIRE(chua.getParamName(0) == "alpha");
    
    // Test default IC
    auto ic = chua.getDefaultInitialCondition();
    REQUIRE(ic.size() == 3);
}

TEST_CASE("ChuaCircuit - Derivatives computation", "[DynamicalSystem][ChuaCircuit]")
{
    ChuaCircuit chua;
    Vector<Real> y({0.7, 0.0, 0.0});
    Vector<Real> dydt(3);
    
    chua.derivs(0, y, dydt);
    
    // Derivatives should be well-defined
    REQUIRE(std::isfinite(dydt[0]));
    REQUIRE(std::isfinite(dydt[1]));
    REQUIRE(std::isfinite(dydt[2]));
}

TEST_CASE("HenonHeilesSystem - Hamiltonian properties", "[DynamicalSystem][HenonHeiles]")
{
    HenonHeilesSystem hh;
    
    REQUIRE(hh.getDim() == 4);
    REQUIRE(hh.isHamiltonian() == true);
    REQUIRE(hh.isDissipative() == false);
    REQUIRE(hh.getNumInvariants() == 1);
    REQUIRE(hh.getInvariantName(0) == "Energy");
}

TEST_CASE("HenonHeilesSystem - Energy conservation", "[DynamicalSystem][HenonHeiles]")
{
    HenonHeilesSystem hh;
    auto ic = hh.getDefaultInitialCondition();
    
    Real E0 = hh.computeInvariant(0, ic);
    REQUIRE(std::isfinite(E0));
    REQUIRE(E0 > 0);  // Positive energy for default IC
    
    // Test derivatives
    Vector<Real> dydt(4);
    hh.derivs(0, ic, dydt);
    
    REQUIRE(dydt[0] == Approx(ic[2]));  // dx/dt = px
    REQUIRE(dydt[1] == Approx(ic[3]));  // dy/dt = py
}

TEST_CASE("DoublePendulumSystem - Basic construction", "[DynamicalSystem][DoublePendulum]")
{
    DoublePendulumSystem dp;
    
    REQUIRE(dp.getDim() == 4);
    REQUIRE(dp.isDissipative() == false);  // Conservative
    REQUIRE(dp.getNumInvariants() == 1);
}

TEST_CASE("DoublePendulumSystem - Energy computation", "[DynamicalSystem][DoublePendulum]")
{
    DoublePendulumSystem dp;
    auto ic = dp.getDefaultInitialCondition();
    
    Real E = dp.computeInvariant(0, ic);
    REQUIRE(std::isfinite(E));
    
    // Small angles at rest
    Vector<Real> rest({0.1, 0.1, 0.0, 0.0});
    Real E_rest = dp.computeInvariant(0, rest);
    REQUIRE(std::isfinite(E_rest));
}

TEST_CASE("DoublePendulumSystem - Acceleration at non-trivial state", "[DynamicalSystem][DoublePendulum]")
{
    // m1=m2=1, L1=L2=1, g=10 — matches analysis verification point
    DoublePendulumSystem dp(1.0, 1.0, 1.0, 1.0, 10.0);

    // State with non-zero angular velocities (where the bug manifests)
    Vector<Real> y({0.5, 0.3, 1.0, 2.0});  // theta1, theta2, omega1, omega2
    Vector<Real> dydt(4);
    dp.derivs(0.0, y, dydt);

    // dtheta/dt = omega (trivial check)
    REQUIRE(dydt[0] == Approx(1.0).epsilon(1e-12));
    REQUIRE(dydt[1] == Approx(2.0).epsilon(1e-12));

    // alpha1 ≈ -7.39 from Euler-Lagrange derivation (Cramer's rule solution)
    REQUIRE(dydt[2] == Approx(-7.39).epsilon(0.01));
    // alpha2 ≈ +4.49
    REQUIRE(dydt[3] == Approx(4.49).epsilon(0.02));
}

TEST_CASE("DoublePendulumSystem - Energy conservation (dE/dt = 0)", "[DynamicalSystem][DoublePendulum]")
{
    DoublePendulumSystem dp(1.0, 1.0, 1.0, 1.0, 10.0);

    // Non-trivial state with angular velocities
    Vector<Real> y({0.5, 0.3, 1.0, 2.0});
    Vector<Real> dydt(4);
    dp.derivs(0.0, y, dydt);

    // Compute dE/dt = grad(E) . dydt via central finite differences
    Real eps = 1e-7;
    Real dEdt = 0.0;
    for (int i = 0; i < 4; ++i) {
        Vector<Real> yp = y, ym = y;
        yp[i] += eps;
        ym[i] -= eps;
        Real dEdyi = (dp.getEnergy(yp) - dp.getEnergy(ym)) / (2 * eps);
        dEdt += dEdyi * dydt[i];
    }

    // dE/dt must be zero for correct Euler-Lagrange equations
    REQUIRE(std::abs(dEdt) < 1e-4);
}

///////////////////////////////////////////////////////////////////////////////////////////
//                                   DISCRETE MAPS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("LogisticMap - Basic iteration", "[DynamicalSystem][DiscreteMap][LogisticMap]")
{
    LogisticMap lmap(3.9);
    
    Vector<Real> x({0.5});
    auto x1 = lmap.iterate(x);
    
    // x1 = r * x * (1-x) = 3.9 * 0.5 * 0.5 = 0.975
    REQUIRE(x1[0] == Approx(0.975));
    
    // Test multiple iterations
    Vector<Real> x5 = x;
    for (int i = 0; i < 5; ++i)
        x5 = lmap.iterate(x5);
    REQUIRE(x5.size() == 1);
    REQUIRE(x5[0] >= 0);
    REQUIRE(x5[0] <= 1);
}

TEST_CASE("LogisticMap - Chaos at r=4", "[DynamicalSystem][DiscreteMap][LogisticMap]")
{
    LogisticMap lmap(4.0);
    
    // Known exact Lyapunov exponent
    REQUIRE(lmap.analyticalLyapunov() == Approx(std::log(2.0)));
    
    // Generate orbit - should stay in [0,1]
    Vector<Real> x0({0.3});
    auto orbit = lmap.orbit(x0, 100);
    
    REQUIRE(orbit.size() == 100);
    for (const auto& x : orbit)
    {
        REQUIRE(x[0] >= 0);
        REQUIRE(x[0] <= 1);
    }
}

TEST_CASE("HenonMap - Basic iteration", "[DynamicalSystem][DiscreteMap][HenonMap]")
{
    HenonMap hmap(1.4, 0.3);
    
    REQUIRE(hmap.getA() == 1.4);
    REQUIRE(hmap.getB() == 0.3);
    
    Vector<Real> x({0.0, 0.0});
    auto x1 = hmap.iterate(x);
    
    // x1 = (1 - 1.4*0 + 0, 0.3*0) = (1, 0)
    REQUIRE(x1[0] == Approx(1.0));
    REQUIRE(x1[1] == Approx(0.0));
    
    // Constant area contraction
    REQUIRE(hmap.jacobianDeterminant() == Approx(-0.3));
}

TEST_CASE("HenonMap - Strange attractor", "[DynamicalSystem][DiscreteMap][HenonMap]")
{
    HenonMap hmap(1.4, 0.3);
    
    Vector<Real> x0({0.1, 0.1});
    auto orbit = hmap.orbit(x0, 1000);
    
    REQUIRE(orbit.size() == 1000);
    
    // Attractor is bounded (approximately |x| < 2, |y| < 1)
    for (int i = 500; i < 1000; ++i)  // After transient
    {
        REQUIRE(std::abs(orbit[i][0]) < 3);
        REQUIRE(std::abs(orbit[i][1]) < 2);
    }
}

TEST_CASE("StandardMap - Area preservation", "[DynamicalSystem][DiscreteMap][StandardMap]")
{
    StandardMap smap(0.5);
    
    REQUIRE(smap.isAreaPreserving() == true);
    
    Vector<Real> x({1.0, 1.0});
    auto x1 = smap.iterate(x);
    
    // Jacobian determinant should be 1
    Matrix<Real> J;
    smap.jacobian(x, J);
    Real det = J(0,0)*J(1,1) - J(0,1)*J(1,0);
    REQUIRE(det == Approx(1.0));
}

TEST_CASE("TentMap - Basic iteration", "[DynamicalSystem][DiscreteMap][TentMap]")
{
    TentMap tmap(2.0);
    
    // Below 0.5
    Vector<Real> x1({0.3});
    auto y1 = tmap.iterate(x1);
    REQUIRE(y1[0] == Approx(0.6));  // 2 * 0.3
    
    // Above 0.5
    Vector<Real> x2({0.7});
    auto y2 = tmap.iterate(x2);
    REQUIRE(y2[0] == Approx(0.6));  // 2 * (1 - 0.7)
    
    // Exact Lyapunov exponent
    REQUIRE(tmap.analyticalLyapunov() == Approx(std::log(2.0)));
}

TEST_CASE("TentMap - Lyapunov exponent", "[DynamicalSystem][DiscreteMap][TentMap]")
{
    TentMap tmap(2.0);
    
    // Numerical Lyapunov should match analytical
    Vector<Real> x0({0.1234});
    auto result = DiscreteMapLyapunov<1>::Compute(tmap, x0, 5000, 500);
    
    REQUIRE(result.isChaotic == true);
    REQUIRE(result.maxExponent == Approx(std::log(2.0)).margin(0.1));
}

TEST_CASE("DiscreteMapLyapunov - Henon map chaos detection", "[DynamicalSystem][DiscreteMap][Lyapunov]")
{
    HenonMap hmap(1.4, 0.3);
    
    Vector<Real> x0({0.1, 0.1});
    auto result = DiscreteMapLyapunov<2>::Compute(hmap, x0, 10000, 1000);
    
    REQUIRE(result.exponents.size() == 2);
    REQUIRE(result.isChaotic == true);
    REQUIRE(result.maxExponent > 0);  // Positive Lyapunov = chaos
    
    // For Henon with a=1.4, b=0.3: λ₁ ≈ 0.42, λ₂ ≈ -1.62
    REQUIRE(result.exponents[0] > 0.3);
    REQUIRE(result.exponents[1] < 0);
}

} // namespace MML::Tests::Systems::DynamicalSystemTests

///////////////////////////////////////////////////////////////////////////////////////////
//                      DYNAMICAL SYSTEM DETAILED API TESTS
///////////////////////////////////////////////////////////////////////////////////////////
namespace MML::Tests::Systems::DynSysDetailedTests {

using namespace MML;
using namespace MML::Systems;
using Catch::Approx;

TEST_CASE("FindFixedPointDetailed - converged fixed point", "[DynamicalSystem][Detailed]")
{
    VanDerPolSystem vdp(1.0);
    Vector<Real> guess({0.1, 0.1});

    auto result = FindFixedPointDetailed(vdp, guess);

    REQUIRE(result.IsSuccess());
    REQUIRE(result.algorithm_name == "FixedPointFinder");
    REQUIRE(result.elapsed_time_ms >= 0.0);
    REQUIRE(result.function_evaluations > 0);
    REQUIRE(result.fixed_point.convergenceResidual < 1e-8);
    REQUIRE(result.fixed_point.location[0] == Approx(0.0).margin(1e-8));
    REQUIRE(result.fixed_point.location[1] == Approx(0.0).margin(1e-8));
    REQUIRE(result.fixed_point.isStable == false);
}

TEST_CASE("FindFixedPointDetailed - non-convergence sets status", "[DynamicalSystem][Detailed]")
{
    LorenzSystem lorenz;  // Default params: chaotic, far-off guess won't converge well with 3 iters
    Vector<Real> guess({100.0, 100.0, 100.0});

    auto result = FindFixedPointDetailed(lorenz, guess, 1e-10, 3);

    // With only 3 iterations and a far-off guess, Newton likely won't converge
    if (!result.IsSuccess()) {
        REQUIRE(result.status == AlgorithmStatus::MaxIterationsExceeded);
        REQUIRE(!result.error_message.empty());
    }
    REQUIRE(result.algorithm_name == "FixedPointFinder");
    REQUIRE(result.elapsed_time_ms >= 0.0);
}

TEST_CASE("FindFixedPointDetailed - error suppressed with ConvertToStatus", "[DynamicalSystem][Detailed]")
{
    VanDerPolSystem vdp(1.0);
    Vector<Real> guess({0.1, 0.1});

    DynSysConfig config;
    config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;

    auto result = FindFixedPointDetailed(vdp, guess, 1e-10, 50, config);

    REQUIRE(result.IsSuccess());
    REQUIRE(result.fixed_point.convergenceResidual < 1e-8);
}

TEST_CASE("ComputeLyapunovDetailed - Lorenz chaotic regime", "[DynamicalSystem][Detailed]")
{
    LorenzSystem lorenz;
    Vector<Real> x0({1.0, 1.0, 1.0});

    // Short integration for test speed
    auto result = ComputeLyapunovDetailed(lorenz, x0, 10.0, 1.0, 0.01);

    REQUIRE(result.IsSuccess());
    REQUIRE(result.algorithm_name == "LyapunovAnalyzer");
    REQUIRE(result.elapsed_time_ms >= 0.0);
    REQUIRE(result.function_evaluations > 0);
    REQUIRE(result.lyapunov.exponents.size() == 3);
    // Even with short integration, largest exponent should be positive (chaos)
    REQUIRE(result.lyapunov.maxExponent > 0.0);
}

TEST_CASE("ComputeLyapunovDetailed - error suppressed", "[DynamicalSystem][Detailed]")
{
    LorenzSystem lorenz;
    Vector<Real> x0({1.0, 1.0, 1.0});

    DynSysConfig config;
    config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;

    auto result = ComputeLyapunovDetailed(lorenz, x0, 10.0, 1.0, 0.01, config);

    REQUIRE(result.IsSuccess());
    REQUIRE(result.lyapunov.exponents.size() == 3);
}

TEST_CASE("SweepBifurcationDetailed - basic parameter sweep", "[DynamicalSystem][Detailed]")
{
    LorenzSystem lorenz;
    Vector<Real> x0({1.0, 1.0, 1.0});

    auto result = SweepBifurcationDetailed(
        lorenz, 1, 20.0, 30.0, 5, x0, 2, 50.0, 20.0, 0.01);

    REQUIRE(result.IsSuccess());
    REQUIRE(result.algorithm_name == "BifurcationAnalyzer");
    REQUIRE(result.elapsed_time_ms >= 0.0);
    REQUIRE(result.function_evaluations == 5);
    REQUIRE(result.diagram.parameterValues.size() == 5);
    REQUIRE(result.diagram.attractorValues.size() == 5);
    REQUIRE(result.diagram.parameterName == "rho");
}

TEST_CASE("ComputePoincareSectionDetailed - Lorenz section", "[DynamicalSystem][Detailed]")
{
    LorenzSystem lorenz;
    Vector<Real> x0({1.0, 1.0, 20.0});

    PoincareSection<Real> section;
    section.variable = 2;   // z-axis
    section.value = 25.0;
    section.direction = 1;  // crossing upward

    auto result = ComputePoincareSectionDetailed(lorenz, x0, section, 5, 0.01);

    REQUIRE(result.IsSuccess());
    REQUIRE(result.algorithm_name == "PoincareSection");
    REQUIRE(result.elapsed_time_ms >= 0.0);
    REQUIRE(result.intersections.size() > 0);
}

TEST_CASE("IntegrateTrajectoryDetailed - simple trajectory", "[DynamicalSystem][Detailed]")
{
    LorenzSystem lorenz;
    Vector<Real> x0({1.0, 1.0, 1.0});

    auto result = IntegrateTrajectoryDetailed(lorenz, x0, 5.0, 0.1, 0.01);

    REQUIRE(result.IsSuccess());
    REQUIRE(result.algorithm_name == "TrajectoryIntegration");
    REQUIRE(result.elapsed_time_ms >= 0.0);
    REQUIRE(result.points.size() > 0);
    REQUIRE(result.function_evaluations == static_cast<int>(result.points.size()));
}

} // namespace MML::Tests::Systems::DynSysDetailedTests

