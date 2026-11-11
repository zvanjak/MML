///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        result_honesty_tests.cpp                                            ///
///  Description: Populator audit for docs/RESULT_POLICY.md (MinimalMathLibrary-      ///
///               wewi.3.2): every converged / status / estimate field must be earned.///
///               Each case here guards a site that previously claimed Success,       ///
///               converged = true, or a literal 0.0 error without checking anything. ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#include <mml/MMLBase.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/base/Geometry/Geometry3D.h>
#include <mml/base/Geometry/Geometry3DBodies.h>
#include <mml/core/Integration/SurfaceIntegration.h>
#include <mml/core/Integration/MonteCarloIntegration.h>
#include <mml/core/Surfaces.h>
#include <mml/base/ComplexFunction.h>
#include <mml/core/ComplexAnalysis.h>
#include <mml/algorithms/RootFinding.h>
#include <mml/algorithms/Eigen/EigenSolver.h>
#include <mml/algorithms/Eigen/SymmMatEigenSolverQR.h>
#include <mml/algorithms/Eigen/ComplexEigenSolver.h>

using namespace MML;

namespace MML::Tests::Core::ResultHonestyTests {

namespace {
    // F = (x^2, y^2, z^2): flux through a tilted triangle is not resolved by one centroid sample
    class QuadraticField : public IVectorFunction<3> {
    public:
        VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const override {
            return VectorN<Real, 3>{x[0] * x[0], x[1] * x[1], x[2] * x[2]};
        }
    };
    class ConstantZField : public IVectorFunction<3> {
    public:
        VectorN<Real, 3> operator()(const VectorN<Real, 3>&) const override { return VectorN<Real, 3>{0.0, 0.0, 1.0}; }
    };
    // Highly oscillatory field so that the depth budget runs out before the tolerance is met
    class OscillatoryField : public IVectorFunction<3> {
    public:
        VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const override {
            return VectorN<Real, 3>{MML::Real(0.0), MML::Real(0.0), std::sin(MML::Real(200.0) * x[0]) * std::cos(MML::Real(170.0) * x[1])};
        }
    };
    template<int N>
    struct ScalarFn : public IScalarFunction<N> {
        std::function<Real(const VectorN<Real, N>&)> f;
        explicit ScalarFn(std::function<Real(const VectorN<Real, N>&)> fn) : f(std::move(fn)) {}
        Real operator()(const VectorN<Real, N>& x) const override { return f(x); }
    };
}

///////////////////////////////////////////////////////////////////////////////////////////
//                               Surface integration
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Adaptive surface integrals report convergence and error estimate honestly", "[result-honesty][surface][wewi.3.2]")
{
    Triangle3D triangle(Point3Cartesian(0, 0, 0), Point3Cartesian(1, 0, 0), Point3Cartesian(0, 1, 0));
    RectSurface3D square(Point3Cartesian(0, 0, 0), Point3Cartesian(1, 0, 0), Point3Cartesian(1, 1, 0), Point3Cartesian(0, 1, 0));

    SECTION("constant field: exact after one refinement, converged with a zero change")
    {
        ConstantZField f;
        auto r = SurfaceIntegration::SurfaceIntegralDetailed(f, square, static_cast<MML::Real>(1e-8));
        REQUIRE(r.converged);
        REQUIRE(r.error_estimated);
        REQUIRE(HasEstimate(r.error_estimate));
        REQUIRE(r.error_estimate <= static_cast<MML::Real>(1e-8));
        REQUIRE_THAT(r.value, Catch::Matchers::WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-12)));
    }

    SECTION("depth budget exhausted (maxLevel 0) with a tolerance the single step cannot meet -> not converged")
    {
        QuadraticField f;
        // tilted triangle in the plane x + y + z = 1 so F.n varies quadratically over it
        Triangle3D tilted(Point3Cartesian(1, 0, 0), Point3Cartesian(0, 1, 0), Point3Cartesian(0, 0, 1));
        auto r = SurfaceIntegration::SurfaceIntegralDetailed(f, tilted, static_cast<MML::Real>(1e-15), /*maxLevel*/ 0);
        INFO("value " << r.value << " error " << r.error_estimate);
        REQUIRE_FALSE(r.converged);                     // used to be hard-coded true
        REQUIRE(HasEstimate(r.error_estimate));
        REQUIRE(r.error_estimate > static_cast<MML::Real>(0.0));          // used to be a literal 0.0
    }

    SECTION("oscillatory field with a shallow depth budget on a square -> not converged, body sum inherits it")
    {
        OscillatoryField f;
        auto r = SurfaceIntegration::SurfaceIntegralDetailed(f, square, static_cast<MML::Real>(1e-12), /*maxLevel*/ 2);
        REQUIRE_FALSE(r.converged);
        REQUIRE(r.error_estimate > static_cast<MML::Real>(0.0));

        Cube3D cube(static_cast<MML::Real>(1.0), Point3Cartesian(0.5, 0.5, 0.5));
        auto body = SurfaceIntegration::SurfaceIntegralDetailed(f, cube, static_cast<MML::Real>(1e-12));   // default maxLevel 7
        auto bodyLoose = SurfaceIntegration::SurfaceIntegralDetailed(ConstantZField(), cube, static_cast<MML::Real>(1e-6));
        REQUIRE(bodyLoose.converged);                   // every face converged
        REQUIRE(HasEstimate(body.error_estimate));      // aggregated, never a placeholder
    }

    SECTION("parametric midpoint rule is a fixed rule: not converged, no estimate; invalid ranges throw")
    {
        ConstantZField f;
        Surfaces::PlaneSurface plane(Vec3Cart(0, 0, 0), Vec3Cart(0, 0, 1), Vec3Cart(1, 0, 0), Vec3Cart(0, 1, 0), 0.0, 1.0, 0.0, 1.0);
        auto r = SurfaceIntegration::SurfaceIntegralDetailed(f, plane, 0.0, 1.0, 0.0, 1.0, 8, 8);
        REQUIRE_FALSE(r.converged);
        REQUIRE_FALSE(r.error_estimated);
        REQUIRE_FALSE(HasEstimate(r.error_estimate));
        REQUIRE(r.iterations == 64);
        REQUIRE(r.function_evaluations == 64);
        REQUIRE_THAT(r.value, Catch::Matchers::WithinAbs(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1e-10)));
        REQUIRE_THROWS_AS(SurfaceIntegration::SurfaceIntegral(f, plane, 1.0, 0.0, 0.0, 1.0), ArgumentError);   // was: value 0, converged
        REQUIRE_THROWS_AS(SurfaceIntegration::SurfaceIntegral(f, plane, 0.0, 1.0, 0.0, 1.0, 0, 4), ArgumentError);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
//                               Monte Carlo
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Fixed-sample Monte Carlo variants are Completed with a statistical estimate", "[result-honesty][montecarlo][wewi.3.2]")
{
    ScalarFn<2> f([](const VectorN<Real, 2>& x) { return x[0] * x[1]; });
    VectorN<Real, 2> lo{0.0, 0.0}, hi{1.0, 1.0};
    MonteCarloIntegrator<2> mc(7);

    auto qmc = mc.integrateQuasiRandom(f, lo, hi, 1 << 12, true);
    REQUIRE(qmc.status == AlgorithmStatus::Completed);     // no tolerance was checked
    REQUIRE(qmc.IsSuccess());
    REQUIRE_FALSE(qmc.converged);
    REQUIRE(HasEstimate(qmc.error_estimate));
    REQUIRE(qmc.error_estimate > static_cast<MML::Real>(0.0));

    // measure-zero domain (upper == lower): exactly 0, nothing sampled -> Completed, exact zero error
    auto zero = mc.integrateQuasiRandom(f, lo, lo, 100, true);
    REQUIRE(zero.status == AlgorithmStatus::Completed);
    REQUIRE_FALSE(zero.converged);
    REQUIRE(zero.value == static_cast<MML::Real>(0.0));
    REQUIRE(zero.error_estimate == static_cast<MML::Real>(0.0));
}

///////////////////////////////////////////////////////////////////////////////////////////
//                               Complex analysis
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ContourIntegrationResult does not claim convergence by default", "[result-honesty][complex][wewi.3.2]")
{
    ComplexAnalysis::ContourIntegrationResult r;
    REQUIRE_FALSE(r.converged);                            // default used to be true
    REQUIRE_FALSE(HasEstimate(r.error_estimate));          // default used to be 0.0

    ComplexFunctionFromStdFunc f([](Complex z) -> Complex { return static_cast<MML::Real>(1.0) / z; });
    ComplexAnalysis::CircleContour unit(Complex(0.0, 0.0), 1.0);
    auto I = ComplexAnalysis::ContourIntegral(f, unit);
    REQUIRE(I.converged);                                  // adaptive Simpson checked and met its tolerance
    REQUIRE(HasEstimate(I.error_estimate));
    REQUIRE(std::abs(I.value - Complex(0.0, 2.0 * Constants::PI)) <= std::max(Real(1e-8), I.error_estimate));
}

///////////////////////////////////////////////////////////////////////////////////////////
//                               Root finding
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Endpoint roots report the residual as achieved tolerance and no x-error estimate", "[result-honesty][rootfinding][wewi.3.2]")
{
    RealFunctionFromStdFunc f([](Real x) { return x * x - static_cast<MML::Real>(4.0); });
    RootFinding::RootFindingConfig cfg;
    cfg.tolerance = MML::Real(1e-12);
    auto r = RootFinding::FindRootBrentDetailed(f, static_cast<MML::Real>(2.0), static_cast<MML::Real>(5.0), cfg);   // f(2) == 0 exactly
    REQUIRE(r.converged);
    REQUIRE(r.root == static_cast<MML::Real>(2.0));
    REQUIRE(r.iterations_used == 0);
    REQUIRE(HasEstimate(r.achieved_tolerance));
    REQUIRE(r.achieved_tolerance == std::abs(r.function_value));         // was literal 0.0
    REQUIRE_FALSE(HasEstimate(r.x_error));                               // no bracket was shrunk
}

TEST_CASE("Polishing a non-zero constant polynomial is Completed, not Success", "[result-honesty][rootfinding][wewi.3.2]")
{
    auto r = RootFinding::PolishPolynomialRoots(PolynomReal::Constant(static_cast<MML::Real>(5.0)), Vector<Complex>());
    REQUIRE(r.status == AlgorithmStatus::Completed);
    REQUIRE(r.IsSuccess());
    REQUIRE_FALSE(r.converged);
    REQUIRE(r.roots.empty());
    REQUIRE_FALSE(HasEstimate(r.max_residual));
}

///////////////////////////////////////////////////////////////////////////////////////////
//                               Eigen solvers
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Eigen solver trivial sizes: n = 0 is Completed, n = 1 is a checked exact result", "[result-honesty][eigen][wewi.3.2]")
{
    auto e0 = EigenSolver::Solve(Matrix<Real>());
    REQUIRE(e0.status == AlgorithmStatus::Completed);
    REQUIRE_FALSE(e0.converged);
    REQUIRE_FALSE(HasEstimate(e0.maxResidual));

    auto e1 = EigenSolver::Solve(Matrix<Real>(1, 1, {static_cast<MML::Real>(3.5)}));
    REQUIRE(e1.converged);
    REQUIRE(e1.status == AlgorithmStatus::Success);
    REQUIRE(e1.maxResidual == static_cast<MML::Real>(0.0));
    REQUIRE(e1.eigenvalues[0].real == static_cast<MML::Real>(3.5));

    auto s0 = SymmMatEigenSolverQR::Solve(MatrixSym<Real>(0));
    REQUIRE(s0.status == AlgorithmStatus::Completed);
    REQUIRE_FALSE(s0.converged);
    MatrixSym<Real> one(1);
    one(0, 0) = static_cast<MML::Real>(-2.0);
    auto s1 = SymmMatEigenSolverQR::Solve(one);
    REQUIRE(s1.converged);
    REQUIRE(s1.residual == static_cast<MML::Real>(0.0));
    REQUIRE(s1.eigenvalues[0] == static_cast<MML::Real>(-2.0));

    auto z = ComplexEigenSolver::Solve(Matrix<Complex>(2, 2));
    REQUIRE(z.converged);
    REQUIRE(z.status == AlgorithmStatus::Success);
    REQUIRE(z.maxResidual == static_cast<MML::Real>(0.0));
}

///////////////////////////////////////////////////////////////////////////////////////////
//                               Function spaces
///////////////////////////////////////////////////////////////////////////////////////////



TEST_CASE("IterativeResultBase defaults do not claim an estimate", "[result-honesty][wewi.3.2]")
{
    IterativeResultBase r;
    REQUIRE_FALSE(HasEstimate(r.achieved_tolerance));       // default used to be 0.0
    IntegrationDetailedResult d;
    REQUIRE_FALSE(HasEstimate(d.error_estimate));
}

} // namespace MML::Tests::Core::ResultHonestyTests
