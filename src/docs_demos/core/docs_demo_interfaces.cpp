///////////////////////////////////////////////////////////////////////////////////////////
/// MML Documentation Demo: Interfaces - Abstract Base Classes
///////////////////////////////////////////////////////////////////////////////////////////
///
/// This demo shows how to implement and use the abstract interfaces defined in mml/interfaces/
/// Each interface is demonstrated with a simple concrete implementation.
///
/// Interfaces covered:
///   - IFunction, IRealFunction, IScalarFunction, IVectorFunction
///   - ICoordTransf, ICoordTransfInv
///   - IDynamicalSystem
///   - IInterval
///   - IODESystem, IODESystemWithJacobian
///   - IODESystemDAE
///   - IODESystemStepCalculator, IODESystemStepper
///   - ITensor, ITensorField
///////////////////////////////////////////////////////////////////////////////////////////

#include <mml/MMLBase.h>
#include <mml/interfaces/IFunction.h>
#include <mml/interfaces/ICoordTransf.h>
#include <mml/base/Vector/VectorTypes3D.h>
#include <mml/interfaces/IDynamicalSystem.h>
#include <mml/interfaces/IInterval.h>
#include <mml/interfaces/IODESystem.h>
#include <mml/interfaces/IODESystemDAE.h>
#include <mml/algorithms/DAESolvers.h>

#include <iostream>
#include <cmath>
#include <memory>

using namespace MML;

namespace MML::docs_demos::interfaces
{
    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              IFunction Interface                                 ///
    ///////////////////////////////////////////////////////////////////////////////////////
    // Base template: IFunction<RetType, ArgType>
    // Specialized interfaces:
    //   - IRealFunction (Real -> Real)
    //   - IScalarFunction<N> (R^N -> Real)
    //   - IVectorFunction<N> (R^N -> R^N)
    //   - IParametricCurve<N> (Real -> R^N)

    /// @brief Example: Custom real function implementing IRealFunction
    class QuadraticFunction : public IRealFunction
    {
    private:
        Real a, b, c;  // ax^2 + bx + c
    public:
        QuadraticFunction(Real a_, Real b_, Real c_) : a(a_), b(b_), c(c_) {}

        Real operator()(Real x) const override {
            return a * x * x + b * x + c;
        }
    };

    void demo_ifunction()
    {
        std::cout << "\n=== IFunction Interface ===\n";

        // Using custom implementation
        QuadraticFunction quad(1.0, -3.0, 2.0);  // x^2 - 3x + 2

        std::cout << "f(x) = x² - 3x + 2\n";
        for (double x = 0; x <= 3; x += 0.5) {
            std::cout << "  f(" << x << ") = " << quad(x) << "\n";
        }

        // Roots are at x=1 and x=2
        std::cout << "Expected roots: x=1 (f=" << quad(1.0) << "), x=2 (f=" << quad(2.0) << ")\n";
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              IScalarFunction Interface                           ///
    ///////////////////////////////////////////////////////////////////////////////////////

    /// @brief Example: 2D scalar function implementing IScalarFunction<2>
    class DistanceFunction : public IScalarFunction<2>
    {
    public:
        Real operator()(const VectorN<Real, 2>& x) const override {
            return std::sqrt(x[0]*x[0] + x[1]*x[1]);  // Distance from origin
        }
    };

    void demo_iscalar_function()
    {
        std::cout << "\n=== IScalarFunction<N> Interface ===\n";

        DistanceFunction dist;

        std::vector<VectorN<Real, 2>> points = {
            VectorN<Real, 2>({0.0, 0.0}),
            VectorN<Real, 2>({1.0, 0.0}),
            VectorN<Real, 2>({3.0, 4.0}),
            VectorN<Real, 2>({1.0, 1.0})
        };

        std::cout << "Distance from origin:\n";
        for (const auto& p : points) {
            std::cout << "  |(" << p[0] << ", " << p[1] << ")| = " << dist(p) << "\n";
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              IVectorFunction Interface                           ///
    ///////////////////////////////////////////////////////////////////////////////////////

    /// @brief Example: 2D rotation implementing IVectorFunction<2>
    class Rotation2D : public IVectorFunction<2>
    {
    private:
        Real cosTheta, sinTheta;
    public:
        Rotation2D(Real theta) : cosTheta(std::cos(theta)), sinTheta(std::sin(theta)) {}

        VectorN<Real, 2> operator()(const VectorN<Real, 2>& x) const override {
            return VectorN<Real, 2>({
                cosTheta * x[0] - sinTheta * x[1],
                sinTheta * x[0] + cosTheta * x[1]
            });
        }
    };

    void demo_ivector_function()
    {
        std::cout << "\n=== IVectorFunction<N> Interface ===\n";

        Rotation2D rotate90(Constants::PI / 2);  // 90 degrees

        VectorN<Real, 2> p({1.0, 0.0});
        VectorN<Real, 2> rotated = rotate90(p);

        std::cout << "Rotating (1, 0) by 90°:\n";
        std::cout << "  Original: (" << p[0] << ", " << p[1] << ")\n";
        std::cout << "  Rotated:  (" << rotated[0] << ", " << rotated[1] << ")\n";
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              IParametricCurve Interface                          ///
    ///////////////////////////////////////////////////////////////////////////////////////

    /// @brief Example: Helix curve implementing IParametricCurve<3>
    class Helix : public IParametricCurve<3>
    {
    private:
        Real radius, pitch;
    public:
        Helix(Real r, Real p) : radius(r), pitch(p) {}

        VectorN<Real, 3> operator()(Real t) const override {
            return VectorN<Real, 3>({
                radius * std::cos(t),
                radius * std::sin(t),
                pitch * t / (2 * Constants::PI)
            });
        }

        Real getMinT() const override { return 0.0; }
        Real getMaxT() const override { return 4 * Constants::PI; }  // Two turns
    };

    void demo_iparametric_curve()
    {
        std::cout << "\n=== IParametricCurve<N> Interface ===\n";

        Helix helix(1.0, 1.0);  // radius=1, pitch=1

        std::cout << "Helix with radius=1, pitch=1:\n";
        for (Real t = 0; t <= 2*Constants::PI; t += Constants::PI/2) {
            auto p = helix(t);
            std::cout << "  t=" << t/Constants::PI << "π: ("
                      << std::fixed << std::setprecision(3)
                      << p[0] << ", " << p[1] << ", " << p[2] << ")\n";
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              ICoordTransf Interface                              ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_icoord_transf()
    {
        std::cout << "\n=== ICoordTransf Interface ===\n";

        std::cout << "ICoordTransf<VectorFrom, VectorTo> (dimensions deduced from the vector types) provides:\n";
        std::cout << "  - FromDim / ToDim: static dimensions, e.g. " << ICoordTransf<Vector3Spherical, Vector3Cartesian>::FromDim
                  << " -> " << ICoordTransf<Vector3Spherical, Vector3Cartesian>::ToDim << "\n";
        std::cout << "  - transf(in) / operator()(in): Transform coordinates (IFunction<VectorTo, const VectorFrom&>)\n";
        std::cout << "  - coordTransfFunc(i): i-th target component as IScalarFunction<FromDim>\n\n";

        std::cout << "Implementations include:\n";
        std::cout << "  - CoordTransfCart2DSpher (Cartesian to spherical)\n";
        std::cout << "  - CoordTransfCart2DCyl (Cartesian to cylindrical)\n";
        std::cout << "  - See docs_demo_coord_transf.cpp for full examples\n";
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              IDynamicalSystem Interface                          ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_idynamical_system()
    {
        std::cout << "\n=== IDynamicalSystem Interface ===\n";

        std::cout << "IDynamicalSystem defines the protocol for dynamical systems:\n";
        std::cout << "  - getDim(): System dimension\n";
        std::cout << "  - getNumParams(): Number of parameters\n";
        std::cout << "  - getParam(i), setParam(i, val): Parameter access\n";
        std::cout << "  - derivs(t, x, dxdt): Compute derivatives\n\n";

        std::cout << "See docs_demo_dynamical_system.cpp for full examples.\n";
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              IInterval Interface                                 ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_iinterval()
    {
        std::cout << "\n=== IInterval Interface ===\n";

        std::cout << "IInterval defines 1D interval operations:\n";
        std::cout << "  - getLower(), getUpper(): Bounds\n";
        std::cout << "  - contains(x): Point containment test\n";
        std::cout << "  - length(): Interval length\n";
        std::cout << "  - midpoint(): Center of interval\n\n";

        std::cout << "See docs_demo_intervals.cpp for full examples.\n";
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              IODESystem Interface                                ///
    ///////////////////////////////////////////////////////////////////////////////////////

    /// @brief Example: Simple harmonic oscillator implementing IODESystem
    class HarmonicOscillator : public IODESystem
    {
    private:
        Real omega;  // Angular frequency
    public:
        HarmonicOscillator(Real w) : omega(w) {}

        int getDim() const override { return 2; }

        void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override {
            // x'' = -omega^2 * x
            // Let x[0] = position, x[1] = velocity
            dxdt[0] = x[1];              // dx/dt = v
            dxdt[1] = -omega*omega * x[0]; // dv/dt = -omega^2 * x
        }

        std::string getVarName(int i) const override {
            return i == 0 ? "x" : "v";
        }
    };

    void demo_iode_system()
    {
        std::cout << "\n=== IODESystem Interface ===\n";

        HarmonicOscillator sho(1.0);  // omega = 1

        std::cout << "Simple Harmonic Oscillator (omega=1):\n";
        std::cout << "  Dimension: " << sho.getDim() << "\n";
        std::cout << "  Variables: " << sho.getVarName(0) << ", " << sho.getVarName(1) << "\n";

        // Evaluate at initial condition
        Vector<Real> x({1.0, 0.0});  // x=1, v=0
        Vector<Real> dxdt(2);
        sho.derivs(0.0, x, dxdt);

        std::cout << "  At (x=1, v=0): dx/dt=" << dxdt[0] << ", dv/dt=" << dxdt[1] << "\n";
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              IODESystemDAE Interface                             ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_iode_system_dae()
    {
        std::cout << "\n=== IODESystemDAE Interface ===\n";

        std::cout << "IODESystemDAE describes semi-explicit index-1 systems:\n";
        std::cout << "  - x' = f(t, x, y), 0 = g(t, x, y)\n";
        std::cout << "  - getDiffDim(), getAlgDim(), diffEqs(), algConstraints()\n";
        std::cout << "  - Requires a nonsingular g_y, not a mass-matrix accessor\n\n";

        std::cout << "See docs_demo_dae_system.cpp for full examples.\n";
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              IODESystemStepper Interface                         ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_iode_stepper()
    {
        std::cout << "\n=== IODESystemStepper Interface ===\n";

        std::cout << "IODESystemStepper defines the integration step protocol:\n";
        std::cout << "  - doStep(sys, x, t, h): Take one integration step\n";
        std::cout << "  - getOrder(): Method order\n";
        std::cout << "  - getName(): Stepper name\n\n";

        std::cout << "IODESystemStepCalculator provides error estimation:\n";
        std::cout << "  - doStepCalc(sys, x, t, h, err): Step with error estimate\n\n";

        std::cout << "Implementations: EulerStepper, RK4Stepper, RKCashKarpStepper, etc.\n";
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              ITensor Interface                                   ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_itensor()
    {
        std::cout << "\n=== ITensor Interface ===\n";

        std::cout << "ITensor<N> defines the interface for N-dimensional tensors:\n";
        std::cout << "  - getRank(): Tensor rank (order)\n";
        std::cout << "  - getDim(): Dimension of each index\n";
        std::cout << "  - operator()(indices...): Component access\n\n";

        std::cout << "See docs_demo_tensors.cpp for full examples.\n";
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              ITensorField Interface                              ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_itensor_field()
    {
        std::cout << "\n=== ITensorField Interface ===\n";

        std::cout << "ITensorField<N> defines tensor-valued functions over space:\n";
        std::cout << "  - operator()(position): Evaluate tensor at position\n";
        std::cout << "  - getRank(): Tensor rank\n";
        std::cout << "  - getDim(): Spatial dimension\n\n";

        std::cout << "Used for metric tensors, stress tensors, etc.\n";
        std::cout << "See docs_demo_tensors.cpp and docs_demo_metric_tensor.cpp.\n";
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              INTERFACE SUMMARY                                   ///
    ///////////////////////////////////////////////////////////////////////////////////////

    void demo_interface_summary()
    {
        std::cout << "\n=== Interface Summary ===\n";

        std::cout << "\n┌─────────────────────────────────────────────────────────────┐\n";
        std::cout << "│ Interface              │ Purpose                            │\n";
        std::cout << "├─────────────────────────────────────────────────────────────┤\n";
        std::cout << "│ IFunction<R,A>         │ Generic function template          │\n";
        std::cout << "│ IRealFunction          │ Real → Real function               │\n";
        std::cout << "│ IScalarFunction<N>     │ R^N → Real function                │\n";
        std::cout << "│ IVectorFunction<N>     │ R^N → R^N function                 │\n";
        std::cout << "│ IParametricCurve<N>    │ Real → R^N curve                   │\n";
        std::cout << "│ ICoordTransf           │ Coordinate transformation          │\n";
        std::cout << "│ IDynamicalSystem       │ Parametrized dynamical system      │\n";
        std::cout << "│ IInterval              │ 1D interval operations             │\n";
        std::cout << "│ IODESystem             │ ODE system for solvers             │\n";
        std::cout << "│ IODESystemDAE          │ Differential-algebraic equations   │\n";
        std::cout << "│ IODESystemStepper      │ Integration step method            │\n";
        std::cout << "│ ITensor<N>             │ N-dimensional tensor               │\n";
        std::cout << "│ ITensorField<N>        │ Tensor field over space            │\n";
        std::cout << "└─────────────────────────────────────────────────────────────┘\n";
    }

    void Run()
    {
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "       MML Documentation Demo: Interfaces\n";
        std::cout << std::string(70, '=') << "\n";

        demo_ifunction();
        demo_iscalar_function();
        demo_ivector_function();
        demo_iparametric_curve();
        demo_icoord_transf();
        demo_idynamical_system();
        demo_iinterval();
        demo_iode_system();
        demo_iode_system_dae();
        demo_iode_stepper();
        demo_itensor();
        demo_itensor_field();
        demo_interface_summary();

        std::cout << "\n" << std::string(70, '=') << "\n";
    }
} // namespace

namespace MML::docs_demos::interfaces {
    class IdentitySampler : public IRealFunction {
    public:
        Real operator()(Real x) const override { return x; }
    };

    class ParametrizedDecay : public IODESystemDAEWithJacobian, public IODESystemDAEParametrized {
        Real rate = 1;
    public:
        int getDiffDim() const override { return 1; }
        int getAlgDim() const override { return 1; }
        int getNumParam() const override { return 1; }
        Real getParam(int i) const override {
            if (i != 0) throw IndexError("Invalid parameter index");
            return rate;
        }
        void setParam(int i, Real value) override {
            if (i != 0) throw IndexError("Invalid parameter index");
            rate = value;
        }
        void diffEqs(Real, const Vector<Real>&, const Vector<Real>& y, Vector<Real>& f) const override { f[0] = -y[0]; }
        void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override { g[0] = y[0] - rate * x[0]; }
        void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 0; }
        void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = -1; }
        void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = -rate; }
        void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& j) const override { j(0, 0) = 1; }
    };
}

void Docs_Demo_InterfaceContracts() {
    using namespace MML::docs_demos::interfaces;
    IdentitySampler identity;
    Vector<Real> x, y;
    const Real maximum = std::numeric_limits<Real>::max();
    identity.getValues(-maximum, maximum, 5, x, y);
    if (x[0] != -maximum || x[2] != 0 || x[4] != maximum || y != x)
        throw NumericalMethodError("Finite extreme sampling contract failed");
    Helix helix(1, 1);
    const auto trace = helix.getTrace(0, Constants::PI, 5);
    const auto singleton = helix.getTrace(0, std::numeric_limits<Real>::quiet_NaN(), 1);
    if (trace.size() != 5 || singleton.size() != 1 || !helix.getTrace(0, 1, -1).empty())
        throw NumericalMethodError("Curve count compatibility failed");

    ParametrizedDecay system;
    system.setParams(Vector<Real>{2});
    IODESystemDAE& base = system;
    IODESystemDAEParametrized& parameters = system;
    if (&base != static_cast<IODESystemDAE*>(&parameters))
        throw InvalidStateError("DAE capabilities do not share one base");
    Vector<Real> x0{1}, y0{2};
    Matrix<Real> gy(1, 1);
    system.jacobian_gy(0, x0, y0, gy);
    if (gy(0, 0) != 1) throw NumericalMethodError("Example is not index-1");
    DAESolverConfig config;
    config.step_size = Real(1) / 8;
    const auto result = SolveDAEBackwardEuler(system, 0, x0, y0, 1, config);
    if (result.status != AlgorithmStatus::Success) throw NumericalMethodError(result.error_message);
    const Real expected = std::pow(1 + 2 * config.step_size, -8);
    if (result.total_steps != 8 ||
        !IsNearlyEqual(result.solution.getXValuesAtEnd()[0], expected, Precision::DefaultToleranceStrict) ||
        !IsNearlyEqual(result.solution.getYValuesAtEnd()[0], 2 * expected, Precision::DefaultToleranceStrict))
        throw NumericalMethodError("Index-1 example violates the independent BE recurrence");
    std::cout << "Finite extreme sampling, trace counts, shared DAE capabilities and g_y=1 verified.\n";
    std::cout << "Decay endpoint: " << result.solution.getXValuesAtEnd()[0] << "\n";
}

void Docs_Demo_Interfaces() {
    MML::docs_demos::interfaces::Run();
}
