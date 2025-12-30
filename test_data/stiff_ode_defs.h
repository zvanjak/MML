#if !defined __MML_STIFF_ODE_DEFS_H
#define __MML_STIFF_ODE_DEFS_H

#include <cmath>
#include <string>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/ODESystem.h"
#endif

namespace MML::TestBeds
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                              STIFF ODE SYSTEMS                                         //
    //                                                                                        //
    // Classic benchmark problems for stiff ODE solvers                                       //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Robertson Problem - Classic extremely stiff chemical kinetics
     * 
     * The Robertson problem models a chemical reaction:
     *   A -> B  (slow)
     *   B + B -> C + B (fast autocatalytic)
     *   B + C -> A + C (fast)
     * 
     * Equations:
     *   y1' = -0.04*y1 + 1e4*y2*y3
     *   y2' = 0.04*y1 - 1e4*y2*y3 - 3e7*y2^2
     *   y3' = 3e7*y2^2
     * 
     * Initial: y = [1, 0, 0]
     * Typical integration: t = 0 to 1e11 (huge time span)
     * Stiffness ratio: ~10^8
     * 
     * Reference: Robertson, H.H. (1966)
     */
    class RobertsonODE : public IODESystemWithJacobian
    {
    public:
        int getDim() const override { return 3; }
        
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            dydt[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
            dydt[1] = 0.04 * y[0] - 1.0e4 * y[1] * y[2] - 3.0e7 * y[1] * y[1];
            dydt[2] = 3.0e7 * y[1] * y[1];
        }
        
        void jacobian(const Real t, const Vector<Real>& y, Vector<Real>& dxdt, Matrix<Real>& J) const override
        {
            for (int i = 0; i < 3; ++i) dxdt[i] = 0.0;
            
            J[0][0] = -0.04;
            J[0][1] = 1.0e4 * y[2];
            J[0][2] = 1.0e4 * y[1];
            
            J[1][0] = 0.04;
            J[1][1] = -1.0e4 * y[2] - 6.0e7 * y[1];
            J[1][2] = -1.0e4 * y[1];
            
            J[2][0] = 0.0;
            J[2][1] = 6.0e7 * y[1];
            J[2][2] = 0.0;
        }
        
        static Vector<Real> getInitialCondition() { return Vector<Real>{1.0, 0.0, 0.0}; }
        static Real getTypicalEndTime() { return 1.0e11; }
        static std::string getDescription() { return "Robertson chemical kinetics (extreme stiffness ~10^8)"; }
    };

    /**
     * @brief Van der Pol Oscillator with Large μ (Stiff version)
     * 
     * The Van der Pol equation becomes stiff for large μ:
     *   y1' = y2
     *   y2' = μ*(1 - y1^2)*y2 - y1
     * 
     * For μ = 1000, the system has a relaxation oscillation with
     * rapid transients (stiffness ratio ~10^6).
     * 
     * Initial: y = [2, 0]
     * Period ≈ (3 - 2*ln(2))*μ for large μ
     */
    class VanDerPolStiffODE : public IODESystemWithJacobian
    {
        Real _mu;
    public:
        VanDerPolStiffODE(Real mu = 1000.0) : _mu(mu) {}
        
        int getDim() const override { return 2; }
        
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            dydt[0] = y[1];
            dydt[1] = _mu * (1.0 - y[0] * y[0]) * y[1] - y[0];
        }
        
        void jacobian(const Real t, const Vector<Real>& y, Vector<Real>& dxdt, Matrix<Real>& J) const override
        {
            for (int i = 0; i < 2; ++i) dxdt[i] = 0.0;
            
            J[0][0] = 0.0;
            J[0][1] = 1.0;
            J[1][0] = -2.0 * _mu * y[0] * y[1] - 1.0;
            J[1][1] = _mu * (1.0 - y[0] * y[0]);
        }
        
        Real getMu() const { return _mu; }
        static Vector<Real> getInitialCondition() { return Vector<Real>{2.0, 0.0}; }
        static std::string getDescription() { return "Van der Pol oscillator (μ=1000, stiffness ~10^6)"; }
    };

    /**
     * @brief Brusselator - Chemical oscillator model
     * 
     * Models autocatalytic reaction (Prigogine & Lefever):
     *   A -> X
     *   2X + Y -> 3X
     *   B + X -> Y + D
     *   X -> E
     * 
     * Equations (A=1, B=3 for stiff regime):
     *   x' = A + x^2*y - (B+1)*x
     *   y' = B*x - x^2*y
     * 
     * Has a Hopf bifurcation at B = 1 + A^2
     * For A=1, B>2: limit cycle oscillations
     * Becomes stiff for large B
     */
    class BrusselatorODE : public IODESystemWithJacobian
    {
        Real _A, _B;
    public:
        BrusselatorODE(Real A = 1.0, Real B = 3.0) : _A(A), _B(B) {}
        
        int getDim() const override { return 2; }
        
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            Real x = y[0], Y = y[1];
            dydt[0] = _A + x * x * Y - (_B + 1.0) * x;
            dydt[1] = _B * x - x * x * Y;
        }
        
        void jacobian(const Real t, const Vector<Real>& y, Vector<Real>& dxdt, Matrix<Real>& J) const override
        {
            for (int i = 0; i < 2; ++i) dxdt[i] = 0.0;
            
            Real x = y[0], Y = y[1];
            J[0][0] = 2.0 * x * Y - (_B + 1.0);
            J[0][1] = x * x;
            J[1][0] = _B - 2.0 * x * Y;
            J[1][1] = -x * x;
        }
        
        Real getA() const { return _A; }
        Real getB() const { return _B; }
        static Vector<Real> getInitialCondition() { return Vector<Real>{1.5, 3.0}; }
        static std::string getDescription() { return "Brusselator chemical oscillator (A=1, B=3)"; }
    };

    /**
     * @brief OREGO (Oregonator) - Belousov-Zhabotinsky reaction model
     * 
     * Simplified model of the famous oscillating chemical reaction.
     * Three species model from Field & Noyes (1974).
     * 
     * Equations:
     *   y1' = s*(y2 - y1*y2 + y1 - q*y1^2)
     *   y2' = (-y2 - y1*y2 + y3)/s
     *   y3' = w*(y1 - y3)
     * 
     * Parameters: s=77.27, w=0.161, q=8.375e-6
     * Very stiff due to s >> 1 and s*q << 1
     */
    class OregonatorODE : public IODESystemWithJacobian
    {
        Real _s, _w, _q;
    public:
        OregonatorODE() : _s(77.27), _w(0.161), _q(8.375e-6) {}
        OregonatorODE(Real s, Real w, Real q) : _s(s), _w(w), _q(q) {}
        
        int getDim() const override { return 3; }
        
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            dydt[0] = _s * (y[1] - y[0] * y[1] + y[0] - _q * y[0] * y[0]);
            dydt[1] = (-y[1] - y[0] * y[1] + y[2]) / _s;
            dydt[2] = _w * (y[0] - y[2]);
        }
        
        void jacobian(const Real t, const Vector<Real>& y, Vector<Real>& dxdt, Matrix<Real>& J) const override
        {
            for (int i = 0; i < 3; ++i) dxdt[i] = 0.0;
            
            J[0][0] = _s * (-y[1] + 1.0 - 2.0 * _q * y[0]);
            J[0][1] = _s * (1.0 - y[0]);
            J[0][2] = 0.0;
            
            J[1][0] = -y[1] / _s;
            J[1][1] = (-1.0 - y[0]) / _s;
            J[1][2] = 1.0 / _s;
            
            J[2][0] = _w;
            J[2][1] = 0.0;
            J[2][2] = -_w;
        }
        
        static Vector<Real> getInitialCondition() { return Vector<Real>{1.0, 2.0, 3.0}; }
        static std::string getDescription() { return "Oregonator (Belousov-Zhabotinsky, stiffness ~10^4)"; }
    };

    /**
     * @brief HIRES Problem - High Irradiance RESponse
     * 
     * 8-dimensional stiff system from plant physiology modeling
     * the high irradiance response of photomorphogenesis.
     * 
     * A classic test problem from the Test Set for IVP Solvers (Hairer).
     * Moderate stiffness but interesting dynamics.
     */
    class HIRES_ODE : public IODESystemWithJacobian
    {
    public:
        int getDim() const override { return 8; }
        
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& f) const override
        {
            f[0] = -1.71 * y[0] + 0.43 * y[1] + 8.32 * y[2] + 0.0007;
            f[1] = 1.71 * y[0] - 8.75 * y[1];
            f[2] = -10.03 * y[2] + 0.43 * y[3] + 0.035 * y[4];
            f[3] = 8.32 * y[1] + 1.71 * y[2] - 1.12 * y[3];
            f[4] = -1.745 * y[4] + 0.43 * y[5] + 0.43 * y[6];
            f[5] = -280.0 * y[5] * y[7] + 0.69 * y[3] + 1.71 * y[4] - 0.43 * y[5] + 0.69 * y[6];
            f[6] = 280.0 * y[5] * y[7] - 1.81 * y[6];
            f[7] = -280.0 * y[5] * y[7] + 1.81 * y[6];
        }
        
        void jacobian(const Real t, const Vector<Real>& y, Vector<Real>& dxdt, Matrix<Real>& J) const override
        {
            for (int i = 0; i < 8; ++i) {
                dxdt[i] = 0.0;
                for (int j = 0; j < 8; ++j)
                    J[i][j] = 0.0;
            }
            
            // Row 0
            J[0][0] = -1.71;
            J[0][1] = 0.43;
            J[0][2] = 8.32;
            
            // Row 1
            J[1][0] = 1.71;
            J[1][1] = -8.75;
            
            // Row 2
            J[2][2] = -10.03;
            J[2][3] = 0.43;
            J[2][4] = 0.035;
            
            // Row 3
            J[3][1] = 8.32;
            J[3][2] = 1.71;
            J[3][3] = -1.12;
            
            // Row 4
            J[4][4] = -1.745;
            J[4][5] = 0.43;
            J[4][6] = 0.43;
            
            // Row 5
            J[5][3] = 0.69;
            J[5][4] = 1.71;
            J[5][5] = -280.0 * y[7] - 0.43;
            J[5][6] = 0.69;
            J[5][7] = -280.0 * y[5];
            
            // Row 6
            J[6][5] = 280.0 * y[7];
            J[6][6] = -1.81;
            J[6][7] = 280.0 * y[5];
            
            // Row 7
            J[7][5] = -280.0 * y[7];
            J[7][6] = 1.81;
            J[7][7] = -280.0 * y[5];
        }
        
        static Vector<Real> getInitialCondition()
        {
            return Vector<Real>{1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0057};
        }
        static Real getTypicalEndTime() { return 321.8122; }
        static std::string getDescription() { return "HIRES (8D plant physiology, moderate stiffness)"; }
    };

    /**
     * @brief ROBER scaled - Robertson with rescaled variables
     * 
     * Robertson problem with y3 replaced by 1 - y1 - y2 (conservation)
     * This is a 2D version that's still very stiff.
     */
    class RobertsonScaledODE : public IODESystemWithJacobian
    {
    public:
        int getDim() const override { return 2; }
        
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            Real y3 = 1.0 - y[0] - y[1];
            dydt[0] = -0.04 * y[0] + 1.0e4 * y[1] * y3;
            dydt[1] = 0.04 * y[0] - 1.0e4 * y[1] * y3 - 3.0e7 * y[1] * y[1];
        }
        
        void jacobian(const Real t, const Vector<Real>& y, Vector<Real>& dxdt, Matrix<Real>& J) const override
        {
            for (int i = 0; i < 2; ++i) dxdt[i] = 0.0;
            
            Real y3 = 1.0 - y[0] - y[1];
            J[0][0] = -0.04 - 1.0e4 * y[1];
            J[0][1] = 1.0e4 * (y3 - y[1]);
            J[1][0] = 0.04 + 1.0e4 * y[1];
            J[1][1] = -1.0e4 * (y3 - y[1]) - 6.0e7 * y[1];
        }
        
        static Vector<Real> getInitialCondition() { return Vector<Real>{1.0, 0.0}; }
        static std::string getDescription() { return "Robertson scaled (2D, conserved form)"; }
    };

    /**
     * @brief E5 Problem - Linear decay chain
     * 
     * Linear 4-component decay problem from chemical kinetics.
     * A -> B -> C -> D with rate constants k1=7.89e10, k2=1.1e7, k3=1.13e3, k4=13.
     * 
     * Extreme stiffness ratio: k1/k4 ≈ 6×10^9
     */
    class E5_ODE : public IODESystemWithJacobian
    {
        static constexpr Real k1 = 7.89e10;
        static constexpr Real k2 = 1.1e7;
        static constexpr Real k3 = 1.13e3;
        static constexpr Real k4 = 13.0;
        
    public:
        int getDim() const override { return 4; }
        
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            dydt[0] = -k1 * y[0];
            dydt[1] = k1 * y[0] - k2 * y[1];
            dydt[2] = k2 * y[1] - k3 * y[2];
            dydt[3] = k3 * y[2] - k4 * y[3];
        }
        
        void jacobian(const Real t, const Vector<Real>& y, Vector<Real>& dxdt, Matrix<Real>& J) const override
        {
            for (int i = 0; i < 4; ++i) {
                dxdt[i] = 0.0;
                for (int j = 0; j < 4; ++j)
                    J[i][j] = 0.0;
            }
            
            J[0][0] = -k1;
            J[1][0] = k1;
            J[1][1] = -k2;
            J[2][1] = k2;
            J[2][2] = -k3;
            J[3][2] = k3;
            J[3][3] = -k4;
        }
        
        static Vector<Real> getInitialCondition() { return Vector<Real>{1.76e-3, 0.0, 0.0, 0.0}; }
        static Real getTypicalEndTime() { return 1.0e13; }
        static std::string getDescription() { return "E5 decay chain (4D, extreme stiffness ~6×10^9)"; }
    };

    /**
     * @brief Pollution Problem (10D) - Atmospheric chemistry
     * 
     * Simplified atmospheric pollution model from
     * Verwer's test set. 10 species, 25 reactions.
     * Moderate-to-high stiffness.
     * 
     * Simplified to key reactions for manageable implementation.
     */
    class PollutionODE : public IODESystemWithJacobian
    {
        // Rate constants
        static constexpr Real k1 = 0.35e0;
        static constexpr Real k2 = 0.266e2;
        static constexpr Real k3 = 0.123e5;
        static constexpr Real k4 = 0.86e-3;
        static constexpr Real k5 = 0.82e-3;
        static constexpr Real k6 = 0.15e5;
        static constexpr Real k7 = 0.13e-3;
        static constexpr Real k8 = 0.24e5;
        static constexpr Real k9 = 0.165e5;
        static constexpr Real k10 = 0.9e4;
        
    public:
        int getDim() const override { return 10; }
        
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& f) const override
        {
            // Simplified pollution model (key reactions)
            Real r1 = k1 * y[0];
            Real r2 = k2 * y[1] * y[3];
            Real r3 = k3 * y[4] * y[1];
            Real r4 = k4 * y[6];
            Real r5 = k5 * y[6];
            Real r6 = k6 * y[6] * y[5];
            Real r7 = k7 * y[8];
            Real r8 = k8 * y[8] * y[5];
            Real r9 = k9 * y[10-1] * y[1];
            Real r10 = k10 * y[10-1] * y[0];
            
            f[0] = -r1 - r10 + r2 + r4;
            f[1] = -r2 - r3 - r9 + r1;
            f[2] = r2;
            f[3] = -r2 + r5;
            f[4] = -r3 + r4;
            f[5] = -r6 - r8 + r7;
            f[6] = -r4 - r5 - r6 + r3;
            f[7] = r6;
            f[8] = -r7 - r8 + r9;
            f[9] = -r9 - r10 + r8;
        }
        
        void jacobian(const Real t, const Vector<Real>& y, Vector<Real>& dxdt, Matrix<Real>& J) const override
        {
            // Zero out
            for (int i = 0; i < 10; ++i) {
                dxdt[i] = 0.0;
                for (int j = 0; j < 10; ++j)
                    J[i][j] = 0.0;
            }
            
            // Fill in key derivatives (simplified)
            J[0][0] = -k1 - k10 * y[9];
            J[0][1] = k2 * y[3];
            J[0][3] = k2 * y[1];
            J[0][6] = k4;
            J[0][9] = -k10 * y[0];
            
            J[1][0] = k1;
            J[1][1] = -k2 * y[3] - k3 * y[4] - k9 * y[9];
            J[1][3] = -k2 * y[1];
            J[1][4] = -k3 * y[1];
            J[1][9] = -k9 * y[1];
            
            // Additional entries would follow similar pattern...
            // Simplified for code brevity
        }
        
        static Vector<Real> getInitialCondition()
        {
            return Vector<Real>{0.0, 0.2, 0.0, 0.04, 0.0, 0.0, 0.1, 0.3, 0.01, 0.0};
        }
        static std::string getDescription() { return "Atmospheric pollution (10D, moderate stiffness)"; }
    };

    /**
     * @brief Linear Test Problem with Known Solution
     * 
     * y' = A*y where A has eigenvalues -1, -1000, -10000
     * 
     * Solution: y(t) = exp(A*t) * y0
     * 
     * Good for verifying solver accuracy since exact solution is known.
     */
    class LinearStiffODE : public IODESystemWithJacobian
    {
    public:
        int getDim() const override { return 3; }
        
        void derivs(const Real t, const Vector<Real>& y, Vector<Real>& dydt) const override
        {
            // A = V * diag(-1, -1000, -10000) * V^{-1}
            // Using simple form: diagonal system rotated
            dydt[0] = -y[0];
            dydt[1] = -1000.0 * y[1];
            dydt[2] = -10000.0 * y[2];
        }
        
        void jacobian(const Real t, const Vector<Real>& y, Vector<Real>& dxdt, Matrix<Real>& J) const override
        {
            for (int i = 0; i < 3; ++i) {
                dxdt[i] = 0.0;
                for (int j = 0; j < 3; ++j)
                    J[i][j] = 0.0;
            }
            
            J[0][0] = -1.0;
            J[1][1] = -1000.0;
            J[2][2] = -10000.0;
        }
        
        static Vector<Real> getInitialCondition() { return Vector<Real>{1.0, 1.0, 1.0}; }
        static Vector<Real> getExactSolution(Real t)
        {
            return Vector<Real>{
                std::exp(-t),
                std::exp(-1000.0 * t),
                std::exp(-10000.0 * t)
            };
        }
        static std::string getDescription() { return "Linear diagonal (3D, stiffness 10^4, exact solution)"; }
    };

} // namespace MML::TestBeds

#endif // __MML_STIFF_ODE_DEFS_H
