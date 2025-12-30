#if !defined MML_ODE_SYSTEM_DEFS_H
#define MML_ODE_SYSTEM_DEFS_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/ODESystem.h"
#endif

namespace MML::TestBeds
{
    /*******************************************************************************************************************
     * ODE SYSTEM DEFINITIONS FOR TESTING
     * 
     * This file contains definitions of ODE systems for testing ODE solvers.
     * Systems are categorized by:
     *   - Type: Linear, Nonlinear, Oscillatory, Chaotic
     *   - Dimensionality: 1D, 2D, 3D, higher
     *   - Solution availability: Analytical, Numerical reference
     * 
     * NOTE: Stiff systems are in a separate test bed (stiff_systems_test_bed.h)
     *******************************************************************************************************************/

    //==================================================================================================================
    // 1. SIMPLE LINEAR SYSTEMS (with analytical solutions)
    //==================================================================================================================
    
    /**
     * @brief Exponential decay: y' = -ky, y(0) = y0
     * Solution: y(t) = y0 * exp(-k*t)
     * 
     * Simplest possible ODE - good for basic verification
     */
    class ExponentialDecayODE : public IODESystem
    {
        Real _k;  // decay constant
    public:
        ExponentialDecayODE(Real k = 1.0) : _k(k) {}
        
        int getDim() const override { return 1; }
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            dydt[0] = -_k * y[0];
        }
        
        Vector<Real> getSolution(const Vector<Real>& y0, Real t) const
        {
            return Vector<Real>{ y0[0] * std::exp(-_k * t) };
        }
        
        Real getDecayConstant() const { return _k; }
    };

    /**
     * @brief Linear growth: y' = ky, y(0) = y0
     * Solution: y(t) = y0 * exp(k*t)
     */
    class ExponentialGrowthODE : public IODESystem
    {
        Real _k;  // growth constant
    public:
        ExponentialGrowthODE(Real k = 1.0) : _k(k) {}
        
        int getDim() const override { return 1; }
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            dydt[0] = _k * y[0];
        }
        
        Vector<Real> getSolution(const Vector<Real>& y0, Real t) const
        {
            return Vector<Real>{ y0[0] * std::exp(_k * t) };
        }
    };

    /**
     * @brief 2D Linear system with distinct real eigenvalues
     * y1' = -y1 + y2
     * y2' = -4*y2
     * 
     * Eigenvalues: λ1 = -1, λ2 = -4
     * Solution: y1(t) = (y10 - y20/3)*exp(-t) + (y20/3)*exp(-4t)
     *           y2(t) = y20*exp(-4t)
     */
    class Linear2DDistinctEigenODE : public IODESystem
    {
    public:
        int getDim() const override { return 2; }
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            dydt[0] = -y[0] + y[1];
            dydt[1] = -4.0 * y[1];
        }
        
        Vector<Real> getSolution(const Vector<Real>& y0, Real t) const
        {
            Real c1 = y0[0] - y0[1] / 3.0;
            Real c2 = y0[1] / 3.0;
            Vector<Real> sol(2);
            sol[0] = c1 * std::exp(-t) + c2 * std::exp(-4.0 * t);
            sol[1] = y0[1] * std::exp(-4.0 * t);
            return sol;
        }
    };

    /**
     * @brief 3D Linear system: y' = Ay where A has eigenvalues 1, 2, 2
     * y1' = y1 + y2 - y3
     * y2' = -y1 + 3*y2 - y3
     * y3' = -y1 + y2 + y3
     * 
     * For initial condition [1, 1, 1]: solution is [exp(t), exp(t)+exp(2t), exp(t)+exp(2t)]
     */
    class Linear3DSystemODE : public IODESystem
    {
    public:
        int getDim() const override { return 3; }
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            dydt[0] = y[0] + y[1] - y[2];
            dydt[1] = -y[0] + 3.0 * y[1] - y[2];
            dydt[2] = -y[0] + y[1] + y[2];
        }
        
        // Specific solution for y0 = [1, 1, 1]
        Vector<Real> getSolutionForUnitIC(Real t) const
        {
            Real e1 = std::exp(t);
            Real e2 = std::exp(2.0 * t);
            return Vector<Real>{ e1, e1 + e2, e1 + e2 };
        }
    };

    //==================================================================================================================
    // 2. HARMONIC OSCILLATORS (with analytical solutions)
    //==================================================================================================================

    /**
     * @brief Simple harmonic oscillator: x'' + ω²x = 0
     * Written as system: y1' = y2, y2' = -ω²y1
     * 
     * Solution: y1(t) = y10*cos(ωt) + (y20/ω)*sin(ωt)
     *           y2(t) = -ω*y10*sin(ωt) + y20*cos(ωt)
     */
    class SimpleHarmonicOscillatorODE : public IODESystem
    {
        Real _omega;  // angular frequency
    public:
        SimpleHarmonicOscillatorODE(Real omega = 1.0) : _omega(omega) {}
        
        int getDim() const override { return 2; }
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            dydt[0] = y[1];
            dydt[1] = -_omega * _omega * y[0];
        }
        
        Vector<Real> getSolution(const Vector<Real>& y0, Real t) const
        {
            Real c = std::cos(_omega * t);
            Real s = std::sin(_omega * t);
            Vector<Real> sol(2);
            sol[0] = y0[0] * c + (y0[1] / _omega) * s;
            sol[1] = -_omega * y0[0] * s + y0[1] * c;
            return sol;
        }
        
        Real getOmega() const { return _omega; }
    };

    /**
     * @brief Damped harmonic oscillator: x'' + 2ζωx' + ω²x = 0
     * Written as system: y1' = y2, y2' = -2ζωy2 - ω²y1
     * 
     * For underdamped case (ζ < 1):
     * ωd = ω*sqrt(1-ζ²)
     * Solution involves exp(-ζωt)*(A*cos(ωd*t) + B*sin(ωd*t))
     */
    class DampedHarmonicOscillatorODE : public IODESystem
    {
        Real _omega;  // natural frequency
        Real _zeta;   // damping ratio
    public:
        DampedHarmonicOscillatorODE(Real omega = 1.0, Real zeta = 0.1) : _omega(omega), _zeta(zeta) {}
        
        int getDim() const override { return 2; }
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            dydt[0] = y[1];
            dydt[1] = -2.0 * _zeta * _omega * y[1] - _omega * _omega * y[0];
        }
        
        // Analytical solution for underdamped case (zeta < 1)
        Vector<Real> getSolution(const Vector<Real>& y0, Real t) const
        {
            if (_zeta >= 1.0) {
                throw std::runtime_error("Analytical solution only implemented for underdamped case (zeta < 1)");
            }
            Real omega_d = _omega * std::sqrt(1.0 - _zeta * _zeta);
            Real decay = std::exp(-_zeta * _omega * t);
            Real c = std::cos(omega_d * t);
            Real s = std::sin(omega_d * t);
            
            // A = y0[0], B = (y0[1] + zeta*omega*y0[0]) / omega_d
            Real A = y0[0];
            Real B = (y0[1] + _zeta * _omega * y0[0]) / omega_d;
            
            Vector<Real> sol(2);
            sol[0] = decay * (A * c + B * s);
            sol[1] = decay * ((-_zeta * _omega * A + omega_d * B) * c + (-_zeta * _omega * B - omega_d * A) * s);
            return sol;
        }
        
        Real getOmega() const { return _omega; }
        Real getZeta() const { return _zeta; }
        bool isUnderdamped() const { return _zeta < 1.0; }
    };

    //==================================================================================================================
    // 3. NONLINEAR OSCILLATORS (no general analytical solutions)
    //==================================================================================================================

    /**
     * @brief Van der Pol oscillator: x'' - μ(1-x²)x' + x = 0
     * Written as system: y1' = y2, y2' = μ(1-y1²)y2 - y1
     * 
     * For small μ (≤ 0.5), the system is non-stiff and oscillates.
     * For large μ (≥ 10), becomes stiff - use stiff test bed.
     */
    class VanDerPolODE : public IODESystem
    {
        Real _mu;
    public:
        VanDerPolODE(Real mu = 0.5) : _mu(mu) {}
        
        int getDim() const override { return 2; }
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            dydt[0] = y[1];
            dydt[1] = _mu * (1.0 - y[0] * y[0]) * y[1] - y[0];
        }
        
        Real getMu() const { return _mu; }
    };

    /**
     * @brief Duffing oscillator: x'' + δx' + αx + βx³ = γcos(ωt)
     * 
     * Classic nonlinear oscillator showing chaos for certain parameters.
     * Default parameters give non-chaotic behavior for testing.
     */
    class DuffingODE : public IODESystem
    {
        Real _delta;  // damping
        Real _alpha;  // linear stiffness
        Real _beta;   // nonlinear stiffness
        Real _gamma;  // forcing amplitude
        Real _omega;  // forcing frequency
    public:
        DuffingODE(Real delta = 0.3, Real alpha = 1.0, Real beta = 0.2, Real gamma = 0.5, Real omega = 1.0)
            : _delta(delta), _alpha(alpha), _beta(beta), _gamma(gamma), _omega(omega) {}
        
        int getDim() const override { return 2; }
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            dydt[0] = y[1];
            dydt[1] = -_delta * y[1] - _alpha * y[0] - _beta * y[0] * y[0] * y[0] + _gamma * std::cos(_omega * t);
        }
    };

    /**
     * @brief Simple pendulum: θ'' + (g/L)sin(θ) = 0
     * Written as system: y1' = y2, y2' = -(g/L)*sin(y1)
     * 
     * For small angles, approximates harmonic oscillator.
     */
    class SimplePendulumODE : public IODESystem
    {
        Real _g_over_L;  // g/L ratio
    public:
        SimplePendulumODE(Real g_over_L = 9.81) : _g_over_L(g_over_L) {}
        
        int getDim() const override { return 2; }
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            dydt[0] = y[1];
            dydt[1] = -_g_over_L * std::sin(y[0]);
        }
        
        // Small angle approximation solution
        Vector<Real> getSmallAngleSolution(const Vector<Real>& y0, Real t) const
        {
            Real omega = std::sqrt(_g_over_L);
            Real c = std::cos(omega * t);
            Real s = std::sin(omega * t);
            Vector<Real> sol(2);
            sol[0] = y0[0] * c + (y0[1] / omega) * s;
            sol[1] = -omega * y0[0] * s + y0[1] * c;
            return sol;
        }
    };

    //==================================================================================================================
    // 4. POPULATION DYNAMICS AND ECOLOGY
    //==================================================================================================================

    /**
     * @brief Lotka-Volterra predator-prey model
     * y1' = αy1 - βy1y2   (prey)
     * y2' = δy1y2 - γy2   (predator)
     * 
     * Classic ecological model with periodic orbits.
     */
    class LotkaVolterraODE : public IODESystem
    {
        Real _alpha;  // prey growth rate
        Real _beta;   // predation rate
        Real _delta;  // predator growth from predation
        Real _gamma;  // predator death rate
    public:
        LotkaVolterraODE(Real alpha = 1.1, Real beta = 0.4, Real delta = 0.1, Real gamma = 0.4)
            : _alpha(alpha), _beta(beta), _delta(delta), _gamma(gamma) {}
        
        int getDim() const override { return 2; }
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            dydt[0] = _alpha * y[0] - _beta * y[0] * y[1];
            dydt[1] = _delta * y[0] * y[1] - _gamma * y[1];
        }
        
        // First integral (constant of motion): V = δx - γln(x) + βy - αln(y)
        Real getFirstIntegral(const Vector<Real>& y) const
        {
            return _delta * y[0] - _gamma * std::log(y[0]) + _beta * y[1] - _alpha * std::log(y[1]);
        }
    };

    /**
     * @brief Logistic growth: y' = ry(1 - y/K)
     * Solution: y(t) = K / (1 + ((K-y0)/y0)*exp(-rt))
     */
    class LogisticGrowthODE : public IODESystem
    {
        Real _r;  // growth rate
        Real _K;  // carrying capacity
    public:
        LogisticGrowthODE(Real r = 1.0, Real K = 10.0) : _r(r), _K(K) {}
        
        int getDim() const override { return 1; }
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            dydt[0] = _r * y[0] * (1.0 - y[0] / _K);
        }
        
        Vector<Real> getSolution(const Vector<Real>& y0, Real t) const
        {
            Real ratio = (_K - y0[0]) / y0[0];
            return Vector<Real>{ _K / (1.0 + ratio * std::exp(-_r * t)) };
        }
    };

    //==================================================================================================================
    // 5. CHAOTIC SYSTEMS (for long-time integration testing)
    //==================================================================================================================

    /**
     * @brief Lorenz system - classic chaotic attractor
     * x' = σ(y - x)
     * y' = x(ρ - z) - y
     * z' = xy - βz
     * 
     * Standard parameters: σ=10, ρ=28, β=8/3 give chaotic behavior
     */
    class LorenzSystemODE : public IODESystemParametrized
    {
        Real _sigma, _rho, _beta;
    public:
        LorenzSystemODE() : _sigma(10.0), _rho(28.0), _beta(8.0 / 3.0) {}
        LorenzSystemODE(Real sigma, Real rho, Real beta) : _sigma(sigma), _rho(rho), _beta(beta) {}
        
        int getNumParam() const override { return 3; }
        Vector<Real> getParams() const override { return Vector<Real>{_sigma, _rho, _beta}; }
        void setParams(const Vector<Real>& params) override { _sigma = params[0]; _rho = params[1]; _beta = params[2]; }
        Real getParam(int i) const override { return i == 0 ? _sigma : (i == 1 ? _rho : _beta); }
        void setParam(int i, Real val) override { if (i == 0) _sigma = val; else if (i == 1) _rho = val; else _beta = val; }
        
        int getDim() const override { return 3; }
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            dydt[0] = _sigma * (y[1] - y[0]);
            dydt[1] = y[0] * (_rho - y[2]) - y[1];
            dydt[2] = y[0] * y[1] - _beta * y[2];
        }
    };

    /**
     * @brief Rössler system - simpler chaotic system
     * x' = -y - z
     * y' = x + ay
     * z' = b + z(x - c)
     * 
     * Standard: a=0.2, b=0.2, c=5.7 gives chaos
     */
    class RosslerSystemODE : public IODESystem
    {
        Real _a, _b, _c;
    public:
        RosslerSystemODE(Real a = 0.2, Real b = 0.2, Real c = 5.7) : _a(a), _b(b), _c(c) {}
        
        int getDim() const override { return 3; }
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            dydt[0] = -y[1] - y[2];
            dydt[1] = y[0] + _a * y[1];
            dydt[2] = _b + y[2] * (y[0] - _c);
        }
    };

    //==================================================================================================================
    // 6. PHYSICAL SYSTEMS
    //==================================================================================================================

    /**
     * @brief Projectile motion with air drag: F = -cv²
     * x'' = -(c/m)*vx*|v|
     * y'' = -g - (c/m)*vy*|v|
     * 
     * Written as 4D system: [x, vx, y, vy]
     */
    class ProjectileWithDragODE : public IODESystem
    {
        Real _c_over_m;  // drag coefficient / mass
        Real _g;         // gravity
    public:
        ProjectileWithDragODE(Real c_over_m = 0.1, Real g = 9.81) : _c_over_m(c_over_m), _g(g) {}
        
        int getDim() const override { return 4; }
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            Real vx = y[1];
            Real vy = y[3];
            Real v = std::sqrt(vx * vx + vy * vy);
            
            dydt[0] = vx;                           // dx/dt = vx
            dydt[1] = -_c_over_m * vx * v;          // dvx/dt
            dydt[2] = vy;                           // dy/dt = vy
            dydt[3] = -_g - _c_over_m * vy * v;     // dvy/dt
        }
    };

    /**
     * @brief Kepler problem (2-body gravitational)
     * Written in polar coordinates: r'' - r*θ'² = -μ/r²
     * 
     * State: [r, r', θ, θ']
     * Conserves angular momentum L = r²θ'
     */
    class KeplerProblemODE : public IODESystem
    {
        Real _mu;  // gravitational parameter (GM)
    public:
        KeplerProblemODE(Real mu = 1.0) : _mu(mu) {}
        
        int getDim() const override { return 4; }
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            Real r = y[0];
            Real rdot = y[1];
            Real thetadot = y[3];
            
            dydt[0] = rdot;                                    // dr/dt
            dydt[1] = r * thetadot * thetadot - _mu / (r * r); // d²r/dt²
            dydt[2] = thetadot;                                // dθ/dt
            dydt[3] = -2.0 * rdot * thetadot / r;              // dθ'/dt (from angular momentum conservation)
        }
        
        // Angular momentum (should be conserved)
        Real getAngularMomentum(const Vector<Real>& y) const
        {
            return y[0] * y[0] * y[3];  // L = r² * θ'
        }
        
        // Total energy (should be conserved)
        Real getTotalEnergy(const Vector<Real>& y) const
        {
            Real r = y[0];
            Real rdot = y[1];
            Real thetadot = y[3];
            Real v2 = rdot * rdot + r * r * thetadot * thetadot;
            return 0.5 * v2 - _mu / r;
        }
    };

    //==================================================================================================================
    // 7. SPECIAL FUNCTION ODEs (Legendre, Laguerre, Hermite, Bessel)
    //==================================================================================================================

    /**
     * @brief Legendre ODE: (1-x²)y'' - 2xy' + n(n+1)y = 0
     * Rewritten: y'' = 2x/(1-x²) y' - n(n+1)/(1-x²) y
     */
    class LegendreODE : public IODESystem
    {
        int _n;  // order
    public:
        LegendreODE(int n) : _n(n) {}
        
        int getDim() const override { return 2; }
        void derivs(const Real x, const Vector<Real>& y, Vector<Real>& dydx) const override
        {
            Real denom = 1.0 - x * x;
            dydx[0] = y[1];
            dydx[1] = 2.0 * x / denom * y[1] - _n * (_n + 1) / denom * y[0];
        }
        
        int getOrder() const { return _n; }
    };

    /**
     * @brief Laguerre ODE: xy'' + (1-x)y' + ny = 0
     * Rewritten: y'' = (x-1)/x y' - n/x y
     */
    class LaguerreODE : public IODESystem
    {
        int _n;  // order
    public:
        LaguerreODE(int n) : _n(n) {}
        
        int getDim() const override { return 2; }
        void derivs(const Real x, const Vector<Real>& y, Vector<Real>& dydx) const override
        {
            dydx[0] = y[1];
            dydx[1] = (x - 1.0) / x * y[1] - static_cast<Real>(_n) / x * y[0];
        }
        
        int getOrder() const { return _n; }
    };

    /**
     * @brief Hermite ODE: y'' - 2xy' + 2ny = 0
     * Rewritten: y'' = 2xy' - 2ny
     */
    class HermiteODE : public IODESystem
    {
        int _n;  // order
    public:
        HermiteODE(int n) : _n(n) {}
        
        int getDim() const override { return 2; }
        void derivs(const Real x, const Vector<Real>& y, Vector<Real>& dydx) const override
        {
            dydx[0] = y[1];
            dydx[1] = 2.0 * x * y[1] - 2.0 * _n * y[0];
        }
        
        int getOrder() const { return _n; }
    };

    /**
     * @brief Bessel ODE: x²y'' + xy' + (x² - n²)y = 0
     * Rewritten: y'' = -y'/x - (1 - n²/x²)y
     */
    class BesselODE : public IODESystem
    {
        int _n;  // order
    public:
        BesselODE(int n) : _n(n) {}
        
        int getDim() const override { return 2; }
        void derivs(const Real x, const Vector<Real>& y, Vector<Real>& dydx) const override
        {
            dydx[0] = y[1];
            dydx[1] = -y[1] / x - (1.0 - static_cast<Real>(_n * _n) / (x * x)) * y[0];
        }
        
        int getOrder() const { return _n; }
    };

    //==================================================================================================================
    // LEGACY FUNCTION-BASED DEFINITIONS (for backward compatibility)
    //==================================================================================================================

    static void VanDerPol(Real mu, const Real t, const Vector<Real>& y, Vector<Real>& dydt) 
    {
        dydt[0] = y[1];
        dydt[1] = mu * (1.0 - y[0] * y[0]) * y[1] - y[0];
    }
    
    static void VanDerPolMju0_1(const Real t, const Vector<Real>& y, Vector<Real>& dydt) 
    { 
        VanDerPol(0.1, t, y, dydt); 
    }

    static void TestLinODESys(Real t, const Vector<Real>& x, Vector<Real>& ret)
    {
        ret[0] = x[0] + x[1] - x[2];
        ret[1] = -x[0] + 3 * x[1] - x[2];
        ret[2] = -x[0] + x[1] + x[2];
    }
    
    static Vector<Real> TestLinODESys_sol(Real t)
    {
        Real e1 = std::exp(t);
        Real e2 = std::exp(2 * t);
        return Vector<Real>{ e1, e1 + e2, e1 + e2 };
    }
    
    static Vector<Real> TestLinODESys_sol1(Vector<Real> initCond, Real t)
    {
        // For IC = [1, 1, 1], solution is [exp(t), exp(t)+exp(2t), exp(t)+exp(2t)]
        Real e1 = std::exp(t);
        Real e2 = std::exp(2 * t);
        return Vector<Real>{ e1, e1 + e2, e1 + e2 };
    }

    //==================================================================================================================
    // STIFF SYSTEM (single example - more in stiff_systems_test_bed.h)
    //==================================================================================================================

    static void stiff_sys1_derivs(Real t, const Vector<Real>& x, Vector<Real>& dydx)
    {
        dydx[0] = -0.013 * x[0] - 1000.0 * x[0] * x[2];
        dydx[1] = -2500.0 * x[1] * x[2];
        dydx[2] = -0.013 * x[0] - 1000.0 * x[0] * x[2] - 2500.0 * x[1] * x[2];
    }
    
    static void stiff_sys1_jac(const Real t, const Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx)
    {
        int n = (int)x.size();
        for (int i = 0; i < n; i++)
            dxdt[i] = 0.0;

        dydx[0][0] = -0.013 - 1000.0 * x[2];
        dydx[0][1] = 0.0;
        dydx[0][2] = -1000.0 * x[0];
        dydx[1][0] = 0.0;
        dydx[1][1] = -2500.0 * x[2];
        dydx[1][2] = -2500.0 * x[1];
        dydx[2][0] = -0.013 - 1000.0 * x[2];
        dydx[2][1] = -2500.0 * x[2];
        dydx[2][2] = -1000.0 * x[0] - 2500.0 * x[1];
    }

    class TestStiffSysWithJacobian1 : public IODESystemWithJacobian
    {
    public:
        int getDim() const override { return 3; }
        
        void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dydx) const override
        {
            stiff_sys1_derivs(t, x, dydx);
        }

        void jacobian(const Real t, const Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx) const override
        {
            stiff_sys1_jac(t, x, dxdt, dydx);
        }
    };
}

#endif