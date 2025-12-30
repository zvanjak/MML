/******************************************************************************
 * DoublePendulum.h - Self-contained Double Pendulum Physics
 * ============================================================================
 * 
 * This header provides all necessary physics for double pendulum simulation,
 * making the example self-contained without MPL dependencies.
 * 
 * Physics: Lagrangian mechanics for coupled pendulum system
 *   - Four state variables: theta1, omega1, theta2, omega2
 *   - Highly nonlinear coupling → deterministic chaos
 *   - Total energy (T + V) should be conserved
 * 
 * The double pendulum is a classic example of chaos theory:
 *   - Deterministic (same initial conditions → same trajectory)
 *   - Yet unpredictable (tiny changes → exponentially different outcomes)
 * 
 *****************************************************************************/

#ifndef DOUBLE_PENDULUM_SELF_CONTAINED_H
#define DOUBLE_PENDULUM_SELF_CONTAINED_H

#include "MMLBase.h"
#include "interfaces/IODESystem.h"

namespace Pendulum
{
    using namespace MML;

    /**************************************************************************
     * DOUBLE PENDULUM ODE SYSTEM
     * 
     * State vector: [theta1, omega1, theta2, omega2]
     *   theta1 = angle of first arm (from vertical)
     *   omega1 = angular velocity of first arm
     *   theta2 = angle of second arm (from vertical)
     *   omega2 = angular velocity of second arm
     * 
     * Equations derived from Lagrangian mechanics:
     *   L = T - V (kinetic minus potential energy)
     *   d/dt(∂L/∂ω) - ∂L/∂θ = 0 (Euler-Lagrange equations)
     **************************************************************************/

    class DoublePendulumODE : public IODESystem
    {
        Real _m1, _m2;  // masses (kg)
        Real _l1, _l2;  // lengths (m)

    public:
        static constexpr Real g = 9.81;  // gravitational acceleration (m/s²)

        DoublePendulumODE(Real m1, Real m2, Real l1, Real l2) 
            : _m1(m1), _m2(m2), _l1(l1), _l2(l2) {}

        // Accessors for physical parameters
        Real getMass1() const { return _m1; }
        Real getMass2() const { return _m2; }
        Real getLength1() const { return _l1; }
        Real getLength2() const { return _l2; }

        int getDim() const override { return 4; }

        void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
        {
            Real th1 = x[0];  // theta1
            Real w1 = x[1];   // omega1
            Real th2 = x[2];  // theta2
            Real w2 = x[3];   // omega2

            // Common denominator in equations of motion
            Real divisor = (2 * _m1 + _m2 - _m2 * cos(2 * th1 - 2 * th2));

            // d(theta1)/dt = omega1
            dxdt[0] = w1;

            // d(omega1)/dt = [complex expression from Lagrangian]
            dxdt[1] = (-g * (2 * _m1 + _m2) * sin(th1) 
                       - _m2 * g * sin(th1 - 2 * th2)
                       - 2 * sin(th1 - th2) * _m2 * (w2 * w2 * _l2 + w1 * w1 * _l1 * cos(th1 - th2))
                      ) / (_l1 * divisor);

            // d(theta2)/dt = omega2
            dxdt[2] = w2;

            // d(omega2)/dt = [complex expression from Lagrangian]
            dxdt[3] = (2 * sin(th1 - th2) * (w1 * w1 * _l1 * (_m1 + _m2)
                       + g * (_m1 + _m2) * cos(th1)
                       + w2 * w2 * _l2 * _m2 * cos(th1 - th2))
                      ) / (_l2 * divisor);
        }

        std::string getVarName(int ind) const override
        {
            switch (ind)
            {
            case 0: return "theta1";   // angle of first pendulum
            case 1: return "omega1";   // angular velocity of first pendulum
            case 2: return "theta2";   // angle of second pendulum
            case 3: return "omega2";   // angular velocity of second pendulum
            default: return IODESystem::getVarName(ind);
            }
        }

        /***********************************************************************
         * ENERGY CALCULATIONS
         * 
         * Kinetic energy: T = (1/2) m1 v1² + (1/2) m2 v2²
         * Potential energy: V = m1 g y1 + m2 g y2 (relative to pivot)
         * Total energy: E = T + V (should be conserved)
         ***********************************************************************/

        Real calcKineticEnergy(Real theta1, Real omega1, Real theta2, Real omega2) const
        {
            // v1² = l1² ω1²
            // v2² = l1² ω1² + l2² ω2² + 2 l1 l2 ω1 ω2 cos(θ1 - θ2)
            Real T = 0.5 * _m1 * _l1 * _l1 * omega1 * omega1
                   + 0.5 * _m2 * (_l1 * _l1 * omega1 * omega1 
                                + _l2 * _l2 * omega2 * omega2 
                                + 2 * _l1 * _l2 * omega1 * omega2 * cos(theta1 - theta2));
            return T;
        }

        Real calcPotentialEnergy(Real theta1, Real theta2) const
        {
            // y1 = -l1 cos(θ1)
            // y2 = -l1 cos(θ1) - l2 cos(θ2)
            Real V = -(_m1 + _m2) * g * _l1 * cos(theta1) - _m2 * g * _l2 * cos(theta2);
            return V;
        }

        Real calcTotalEnergy(Real theta1, Real omega1, Real theta2, Real omega2) const
        {
            return calcKineticEnergy(theta1, omega1, theta2, omega2) 
                 + calcPotentialEnergy(theta1, theta2);
        }

        /***********************************************************************
         * COORDINATE TRANSFORMATIONS
         * 
         * Convert from angles to Cartesian coordinates for visualization
         ***********************************************************************/

        void anglesToCartesian(Real theta1, Real theta2,
                               Real& x1, Real& y1, Real& x2, Real& y2) const
        {
            // First bob position
            x1 = _l1 * sin(theta1);
            y1 = -_l1 * cos(theta1);

            // Second bob position (relative to first)
            x2 = x1 + _l2 * sin(theta2);
            y2 = y1 - _l2 * cos(theta2);
        }
    };

} // namespace Pendulum

#endif // DOUBLE_PENDULUM_SELF_CONTAINED_H
