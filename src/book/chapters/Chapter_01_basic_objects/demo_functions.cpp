///////////////////////////////////////////////////////////////////////////////////////////
///  File:        demo_functions.cpp
///  Description: Comprehensive Function abstraction demonstrations for Chapter 01
///               Shows MML's powerful function type system for numerical computing
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Function.h"
#include "mml/base/Vector/VectorTypes.h"
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                        Function Type Overview                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Function_Types()
{
    std::cout << "\n=== MML Function Types ===\n\n";
    
    std::cout << "MML provides a unified function abstraction system:\n\n";
    std::cout << "  RealFunction              f: R -> R        (function pointer)\n";
    std::cout << "  RealFunctionFromStdFunc   f: R -> R        (std::function/lambda)\n";
    std::cout << "  ScalarFunction<N>         f: R^N -> R      (N-dim input, scalar out)\n";
    std::cout << "  VectorFunction<N>         f: R^N -> R^N    (N-dim to N-dim)\n";
    std::cout << "  VectorFunctionNM<N,M>     f: R^N -> R^M    (N-dim to M-dim)\n";
    std::cout << "  ParametricCurve<N>        g: R -> R^N      (curves in N-space)\n";
    std::cout << "  ParametricSurfaceRect<N>  s: R^2 -> R^N    (surfaces in N-space)\n";
    std::cout << "\nTwo flavors for each:\n";
    std::cout << "  - Function pointer version (stateless lambdas ok)\n";
    std::cout << "  - *FromStdFunc version (captures, closures, bound functions)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Real Functions (R -> R)                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Real_Functions()
{
    std::cout << "\n=== Real Functions (R -> R) ===\n\n";

    // Creating from non-capturing lambda (converts to function pointer)
    std::cout << "--- Creating from Stateless Lambda ---\n";
    RealFunction f1([](Real x) { return x * x - 2.0; });           // Quadratic
    RealFunction f2([](Real x) { return std::sin(x) * std::exp(-x * x / 10.0); }); // Damped sine
    RealFunction f3([](Real x) { return 1.0 / (1.0 + std::exp(-x)); }); // Logistic/sigmoid

    std::cout << "f1(x) = x^2 - 2\n";
    std::cout << "  f1(0) = " << f1(0.0) << "\n";
    std::cout << "  f1(1) = " << f1(1.0) << "\n";
    std::cout << "  f1(sqrt(2)) = " << f1(std::sqrt(2.0)) << " (root!)\n\n";

    std::cout << "f2(x) = sin(x) * exp(-x^2/10)  (damped oscillation)\n";
    for (Real x = 0; x <= 6; x += 2) {
        std::cout << "  f2(" << x << ") = " << f2(x) << "\n";
    }

    std::cout << "\nf3(x) = 1/(1+e^(-x))  (sigmoid/logistic)\n";
    std::cout << "  f3(-5) = " << f3(-5.0) << " (approx 0)\n";
    std::cout << "  f3(0)  = " << f3(0.0) << " (= 0.5)\n";
    std::cout << "  f3(5)  = " << f3(5.0) << " (approx 1)\n";

    // Using RealFunctionFromStdFunc for capturing lambdas
    std::cout << "\n--- Capturing Lambdas (RealFunctionFromStdFunc) ---\n";
    Real g = 9.81;   // gravity (m/s^2)
    Real v0 = 20.0;  // initial velocity (m/s)
    Real theta = 45.0 * Constants::PI / 180.0;  // launch angle
    
    RealFunctionFromStdFunc projectile_height([=](Real t) {
        return v0 * std::sin(theta) * t - 0.5 * g * t * t;
    });
    
    std::cout << "Projectile motion: h(t) = v0*sin(theta)*t - 0.5*g*t^2\n";
    std::cout << "  v0 = " << v0 << " m/s, theta = 45 deg\n";
    for (Real t = 0; t <= 3; t += 0.5) {
        Real h = projectile_height(t);
        std::cout << "  h(" << t << "s) = " << h << " m\n";
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Scalar Functions (R^N -> R)                                  ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Scalar_Functions()
{
    std::cout << "\n=== Scalar Functions (R^N -> R) ===\n\n";

    // 3D scalar field - inverse distance (stateless lambda)
    std::cout << "--- 3D Scalar Field: f(x) = 1/|x| ---\n";
    ScalarFunction<3> inverse_r([](const VectorN<Real, 3>& x) {
        return 1.0 / x.NormL2();
    });

    std::cout << "f(1,0,0) = " << inverse_r(VectorN<Real, 3>{1, 0, 0}) << "\n";
    std::cout << "f(1,1,1) = " << inverse_r(VectorN<Real, 3>{1, 1, 1}) << " (1/sqrt(3) = 0.577)\n";
    std::cout << "f(2,0,0) = " << inverse_r(VectorN<Real, 3>{2, 0, 0}) << "\n";

    // Gravitational potential (capturing lambda - use FromStdFunc)
    std::cout << "\n--- Physics: Two-Body Gravitational Potential ---\n";
    VectorN<Real, 3> mass1_pos{-5.0, 0.0, 0.0};
    VectorN<Real, 3> mass2_pos{5.0, 0.0, 0.0};
    Real M1 = 1000.0, M2 = 1000.0, G = 1.0;
    
    ScalarFunctionFromStdFunc<3> gravity_potential([=](const VectorN<Real, 3>& x) {
        Real r1 = (x - mass1_pos).NormL2();
        Real r2 = (x - mass2_pos).NormL2();
        return -G * (M1 / r1 + M2 / r2);
    });

    std::cout << "Two masses at x = +/-5, M = 1000 each\n";
    std::cout << "phi(0,0,0) = " << gravity_potential(VectorN<Real, 3>{0, 0, 0}) << " (midpoint)\n";
    std::cout << "phi(0,5,0) = " << gravity_potential(VectorN<Real, 3>{0, 5, 0}) << " (off-axis)\n";
    std::cout << "phi(10,0,0) = " << gravity_potential(VectorN<Real, 3>{10, 0, 0}) << " (beyond mass2)\n";

    // Temperature field (uses captures - use FromStdFunc)
    std::cout << "\n--- Engineering: Gaussian Temperature Distribution ---\n";
    Real T_max = 100.0;  // Peak temperature
    Real sigma = 2.0;    // Spread
    
    ScalarFunctionFromStdFunc<3> temperature([=](const VectorN<Real, 3>& x) {
        Real r2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
        return T_max * std::exp(-r2 / (2.0 * sigma * sigma));
    });

    std::cout << "T(x) = 100*exp(-r^2/8)  (Gaussian heat source)\n";
    std::cout << "T(0,0,0) = " << temperature(VectorN<Real, 3>{0, 0, 0}) << " C (peak)\n";
    std::cout << "T(2,0,0) = " << temperature(VectorN<Real, 3>{2, 0, 0}) << " C (at sigma)\n";
    std::cout << "T(4,0,0) = " << temperature(VectorN<Real, 3>{4, 0, 0}) << " C (at 2*sigma)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Vector Functions (R^N -> R^N / R^M)                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Vector_Functions()
{
    std::cout << "\n=== Vector Functions (R^N -> R^M) ===\n\n";

    // 2D rotation field (stateless)
    std::cout << "--- 2D Vector Field: Rotation ---\n";
    VectorFunction<2> rotation_field([](const VectorN<Real, 2>& x) {
        return VectorN<Real, 2>{-x[1], x[0]};  // 90 deg rotation
    });

    std::cout << "F(x,y) = (-y, x)  (counterclockwise rotation)\n";
    std::cout << "F(1,0) = " << rotation_field(VectorN<Real, 2>{1, 0}) << "\n";
    std::cout << "F(0,1) = " << rotation_field(VectorN<Real, 2>{0, 1}) << "\n";
    std::cout << "F(1,1) = " << rotation_field(VectorN<Real, 2>{1, 1}) << "\n";

    // 3D magnetic field (uses captures - use FromStdFunc)
    std::cout << "\n--- Physics: Magnetic Field Around Wire (z-axis) ---\n";
    Real I = 10.0;   // Current (A)
    Real mu0 = 4.0 * Constants::PI * 1e-7;  // Permeability
    
    VectorFunctionFromStdFunc<3> B_field([=](const VectorN<Real, 3>& x) {
        Real rho = std::sqrt(x[0]*x[0] + x[1]*x[1]);  // Distance from z-axis
        if (rho < 1e-10) return VectorN<Real, 3>{0, 0, 0};  // Avoid singularity
        Real B_mag = mu0 * I / (2.0 * Constants::PI * rho);
        // B points in phi direction: (-y/rho, x/rho, 0)
        return VectorN<Real, 3>{-x[1]/rho * B_mag, x[0]/rho * B_mag, 0.0};
    });

    std::cout << "B = (mu0*I/2*pi*rho) phi_hat  for wire along z-axis, I = 10A\n";
    auto B1 = B_field(VectorN<Real, 3>{1, 0, 0});
    std::cout << "B(1,0,0) = " << B1 << " T  (points in +y)\n";
    auto B2 = B_field(VectorN<Real, 3>{0, 2, 0});
    std::cout << "B(0,2,0) = " << B2 << " T  (points in -x)\n";

    // VectorFunctionNM: different input/output dimensions (stateless)
    std::cout << "\n--- VectorFunctionNM: R^2 -> R^3 (surface embedding) ---\n";
    VectorFunctionNM<2, 3> sphere_patch([](const VectorN<Real, 2>& uv) {
        Real u = uv[0], v = uv[1];
        return VectorN<Real, 3>{
            std::cos(u) * std::sin(v),
            std::sin(u) * std::sin(v),
            std::cos(v)
        };
    });

    std::cout << "Sphere patch: (u,v) -> (cos(u)sin(v), sin(u)sin(v), cos(v))\n";
    auto p1 = sphere_patch(VectorN<Real, 2>{0, Constants::PI/2});
    std::cout << "At (0, pi/2): " << p1 << " (equator, x=1)\n";
    auto p2 = sphere_patch(VectorN<Real, 2>{Constants::PI/2, Constants::PI/2});
    std::cout << "At (pi/2, pi/2): " << p2 << " (equator, y=1)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Parametric Curves (R -> R^N)                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Parametric_Curves()
{
    std::cout << "\n=== Parametric Curves (R -> R^N) ===\n\n";

    // Helix (capturing lambda - use FromStdFunc)
    std::cout << "--- 3D Curve: Helix ---\n";
    Real radius = 1.0;
    Real pitch = 0.5;  // Rise per revolution
    
    ParametricCurveFromStdFunc<3> helix([=](Real t) {
        return VectorN<Real, 3>{
            radius * std::cos(t),
            radius * std::sin(t),
            pitch * t / (2.0 * Constants::PI)
        };
    });

    std::cout << "gamma(t) = (cos(t), sin(t), t/(4*pi))  (helix, radius=1, pitch=0.5)\n";
    for (Real t = 0; t <= 2 * Constants::PI; t += Constants::PI / 2) {
        auto p = helix(t);
        std::cout << "  gamma(" << t << ") = (" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";
    }

    // Lissajous curve (stateless lambda)
    std::cout << "\n--- 2D Curve: Lissajous Figure ---\n";
    ParametricCurve<2> lissajous([](Real t) {
        return VectorN<Real, 2>{std::sin(3*t), std::sin(2*t)};
    });

    std::cout << "gamma(t) = (sin(3t), sin(2t))  (3:2 Lissajous)\n";
    std::cout << "  gamma(0) = " << lissajous(0) << "\n";
    std::cout << "  gamma(pi/4) = " << lissajous(Constants::PI/4) << "\n";
    std::cout << "  gamma(pi/2) = " << lissajous(Constants::PI/2) << "\n";

    // Cycloid (capturing - use FromStdFunc)
    std::cout << "\n--- Physics: Cycloid (Brachistochrone) ---\n";
    Real R = 1.0;  // Wheel radius
    ParametricCurveFromStdFunc<2> cycloid([=](Real t) {
        return VectorN<Real, 2>{R * (t - std::sin(t)), R * (1 - std::cos(t))};
    });

    std::cout << "Cycloid: gamma(t) = (R(t-sin(t)), R(1-cos(t)))\n";
    std::cout << "(Path traced by point on rolling wheel - also the brachistochrone!)\n";
    for (Real t = 0; t <= 2 * Constants::PI; t += Constants::PI / 2) {
        auto p = cycloid(t);
        std::cout << "  gamma(" << t << ") = (" << p[0] << ", " << p[1] << ")\n";
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Parametric Surfaces (R^2 -> R^N)                             ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Parametric_Surfaces()
{
    std::cout << "\n=== Parametric Surfaces (R^2 -> R^N) ===\n\n";

    // Torus (capturing lambda - use ParametricSurfaceFromStdFunc)
    std::cout << "--- 3D Surface: Torus ---\n";
    Real R = 3.0;  // Major radius (center of tube to center of torus)
    Real r = 1.0;  // Minor radius (tube radius)
    
    std::function<VectorN<Real, 3>(Real, Real)> torus_func = [=](Real u, Real v) {
        return VectorN<Real, 3>{
            (R + r * std::cos(v)) * std::cos(u),
            (R + r * std::cos(v)) * std::sin(u),
            r * std::sin(v)
        };
    };
    ParametricSurfaceFromStdFunc<3> torus(torus_func);

    std::cout << "Torus: sigma(u,v), R=" << R << " (major), r=" << r << " (minor)\n";
    std::cout << "  sigma(0, 0) = " << torus(0, 0) << " (outer edge)\n";
    std::cout << "  sigma(0, pi) = " << torus(0, Constants::PI) << " (inner edge)\n";
    std::cout << "  sigma(pi/2, 0) = " << torus(Constants::PI/2, 0) << " (90 deg around)\n";

    // Saddle surface (stateless)
    std::cout << "\n--- 3D Surface: Hyperbolic Paraboloid (Saddle) ---\n";
    ParametricSurfaceRect<3> saddle([](Real x, Real y) {
        return VectorN<Real, 3>{x, y, x*x - y*y};
    });

    std::cout << "Saddle: z = x^2 - y^2\n";
    std::cout << "  sigma(0,0) = " << saddle(0, 0) << " (saddle point)\n";
    std::cout << "  sigma(1,0) = " << saddle(1, 0) << " (rises along x)\n";
    std::cout << "  sigma(0,1) = " << saddle(0, 1) << " (falls along y)\n";
    std::cout << "  sigma(1,1) = " << saddle(1, 1) << " (1-1=0)\n";

    // Mobius strip (stateless lambda)
    std::cout << "\n--- 3D Surface: Mobius Strip ---\n";
    ParametricSurfaceRect<3> mobius([](Real u, Real v) {
        // u in [0, 2*pi), v in [-0.5, 0.5]
        Real half_twist = u / 2.0;
        return VectorN<Real, 3>{
            (1 + v * std::cos(half_twist)) * std::cos(u),
            (1 + v * std::cos(half_twist)) * std::sin(u),
            v * std::sin(half_twist)
        };
    });

    std::cout << "Mobius strip (non-orientable, one-sided surface)\n";
    std::cout << "  sigma(0, 0) = " << mobius(0, 0) << " (center of strip at u=0)\n";
    std::cout << "  sigma(pi, 0.5) = " << mobius(Constants::PI, 0.5) << " (edge at u=pi)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Parameterized Functions (with adjustable parameters)         ///
///////////////////////////////////////////////////////////////////////////////////////////

// Class-based scalar function with parameters
class GaussianPotential : public IScalarFunction<3>
{
    Real _amplitude;
    Real _sigma;
    VectorN<Real, 3> _center;

public:
    GaussianPotential(Real amplitude, Real sigma, const VectorN<Real, 3>& center)
        : _amplitude(amplitude), _sigma(sigma), _center(center) {}

    void SetAmplitude(Real a) { _amplitude = a; }
    void SetSigma(Real s) { _sigma = s; }
    void SetCenter(const VectorN<Real, 3>& c) { _center = c; }

    Real operator()(const VectorN<Real, 3>& x) const override
    {
        VectorN<Real, 3> diff = x - _center;
        Real r2 = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2];
        return _amplitude * std::exp(-r2 / (2.0 * _sigma * _sigma));
    }
};

void Demo_Parameterized_Functions()
{
    std::cout << "\n=== Parameterized Functions ===\n\n";

    std::cout << "For complex functions with adjustable parameters,\n";
    std::cout << "inherit from I*Function interfaces.\n\n";

    GaussianPotential gaussian(100.0, 2.0, VectorN<Real, 3>{0, 0, 0});
    
    std::cout << "--- GaussianPotential: A*exp(-r^2/2*sigma^2) ---\n";
    std::cout << "Initial: A=100, sigma=2, center=(0,0,0)\n";
    std::cout << "  phi(0,0,0) = " << gaussian(VectorN<Real, 3>{0, 0, 0}) << "\n";
    std::cout << "  phi(2,0,0) = " << gaussian(VectorN<Real, 3>{2, 0, 0}) << "\n";

    gaussian.SetAmplitude(50.0);
    gaussian.SetCenter(VectorN<Real, 3>{1, 0, 0});
    
    std::cout << "\nAfter: A=50, center=(1,0,0)\n";
    std::cout << "  phi(0,0,0) = " << gaussian(VectorN<Real, 3>{0, 0, 0}) << "\n";
    std::cout << "  phi(1,0,0) = " << gaussian(VectorN<Real, 3>{1, 0, 0}) << " (now at peak)\n";
    std::cout << "  phi(2,0,0) = " << gaussian(VectorN<Real, 3>{2, 0, 0}) << "\n";

    // Lambda capture for simpler cases (using RealFunctionFromStdFunc for captures)
    std::cout << "\n--- Lambda Capture Alternative ---\n";
    Real k = 1.0;
    RealFunctionFromStdFunc spring_force([&k](Real x) { return -k * x; });  // Hooke's law
    
    std::cout << "Spring force: F = -kx, k=" << k << "\n";
    std::cout << "  F(1) = " << spring_force(1.0) << "\n";
    
    k = 5.0;  // Change spring constant
    std::cout << "After k=" << k << ":\n";
    std::cout << "  F(1) = " << spring_force(1.0) << " (updated automatically!)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Practical Applications                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Functions_Applications()
{
    std::cout << "\n=== Practical Applications ===\n\n";

    // Electric field from point charges (capturing - use FromStdFunc)
    std::cout << "--- Electrostatics: Electric Field from Point Charge ---\n";
    Real q = 1.0e-6;  // 1 uC
    Real k_e = 8.99e9;  // Coulomb constant
    VectorN<Real, 3> charge_pos{0, 0, 0};
    
    VectorFunctionFromStdFunc<3> E_field([=](const VectorN<Real, 3>& x) {
        VectorN<Real, 3> r = x - charge_pos;
        Real r_mag = r.NormL2();
        if (r_mag < 1e-10) return VectorN<Real, 3>{0, 0, 0};
        Real E_mag = k_e * q / (r_mag * r_mag);
        return r * (E_mag / r_mag);  // E points radially outward
    });

    std::cout << "Point charge q = 1 uC at origin\n";
    auto E1 = E_field(VectorN<Real, 3>{1, 0, 0});
    std::cout << "  E(1m, 0, 0) = " << E1 << " N/C\n";
    std::cout << "  |E| = " << E1.NormL2() << " N/C (should be ~9e3)\n";

    // Orbital mechanics (capturing - use FromStdFunc)
    std::cout << "\n--- Orbital Mechanics: Elliptical Orbit ---\n";
    Real a = 1.5e11;  // Semi-major axis (m) - roughly Earth-Sun
    Real e = 0.5;     // Eccentricity (exaggerated for demo)
    Real b = a * std::sqrt(1 - e*e);
    Real T = 365.25 * 24 * 3600;  // Period (s)
    
    ParametricCurveFromStdFunc<2> ellipse([=](Real t) {
        Real theta = 2 * Constants::PI * t / T;
        return VectorN<Real, 2>{a * std::cos(theta), b * std::sin(theta)};
    });

    std::cout << "Elliptical orbit: a = 1.5e11 m, e = 0.5\n";
    for (int month = 0; month <= 12; month += 3) {
        Real t = month * 30.44 * 24 * 3600;  // Approximate days to seconds
        auto pos = ellipse(t);
        std::cout << "  Month " << month << ": pos = (" << pos[0]/1e11 << ", " << pos[1]/1e11 << ") x 10^11 m\n";
    }

    // Machine learning: Activation functions (stateless)
    std::cout << "\n--- ML: Neural Network Activation Functions ---\n";
    RealFunction sigmoid([](Real x) { return 1.0 / (1.0 + std::exp(-x)); });
    RealFunction tanh_fn([](Real x) { return std::tanh(x); });
    RealFunction relu([](Real x) { return x > 0 ? x : 0.0; });
    RealFunction leaky_relu([](Real x) { return x > 0 ? x : 0.01 * x; });

    std::cout << "Common activation functions at x = -2, 0, 2:\n";
    std::cout << "  sigmoid: " << sigmoid(-2) << ", " << sigmoid(0) << ", " << sigmoid(2) << "\n";
    std::cout << "  tanh:    " << tanh_fn(-2) << ", " << tanh_fn(0) << ", " << tanh_fn(2) << "\n";
    std::cout << "  ReLU:    " << relu(-2) << ", " << relu(0) << ", " << relu(2) << "\n";
    std::cout << "  Leaky:   " << leaky_relu(-2) << ", " << leaky_relu(0) << ", " << leaky_relu(2) << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Main Demo Entry Point                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Functions()
{
    std::cout << "\n";
    std::cout << "***********************************************************************\n";
    std::cout << "****                    FUNCTIONS IN MML                           ****\n";
    std::cout << "****         Unified Abstraction for Numerical Computing           ****\n";
    std::cout << "***********************************************************************\n";

    Demo_Function_Types();
    Demo_Real_Functions();
    Demo_Scalar_Functions();
    Demo_Vector_Functions();
    Demo_Parametric_Curves();
    Demo_Parametric_Surfaces();
    Demo_Parameterized_Functions();
    Demo_Functions_Applications();
}
