/******************************************************************************
 * Projectiles2D.h - Self-contained Projectile Motion Physics
 * ============================================================================
 * 
 * This header provides all necessary physics for 2D projectile simulation,
 * making the example self-contained without MPL dependencies.
 * 
 * Contents:
 *   - IProjectile2DODE: Base interface for projectile ODE systems
 *   - Projectile2DInVacuumODE: Ideal parabolic motion (analytical solutions)
 *   - Projectile2DWithAirResistanceODE: Quadratic drag model
 *   - Projectile2DChangingAirDensityODE: Altitude-dependent drag
 *   - BaseballWithDragCoeffDependentOnSpeedODE: Speed-dependent drag
 *   - ProjectileTrajectory2D: Trajectory storage and interpolation
 *   - ProjectileMotionSolver2D: ODE solver wrapper with ground detection
 * 
 * Physics: F = mg - kv² (quadratic drag model)
 * 
 *****************************************************************************/

#ifndef PROJECTILES_2D_SELF_CONTAINED_H
#define PROJECTILES_2D_SELF_CONTAINED_H

#include "MMLBase.h"
#include "interfaces/IODESystem.h"
#include "base/InterpolatedFunction.h"
#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"

namespace Projectile
{
    using namespace MML;

    /**************************************************************************
     * BASE INTERFACE
     **************************************************************************/

    class IProjectile2DODE : public IODESystem
    {
    public:
        virtual std::string getVarName(int ind) const override
        {
            switch (ind)
            {
            case 0: return "x";     // x position
            case 1: return "y";     // y position
            case 2: return "vx";    // x velocity
            case 3: return "vy";    // y velocity
            default: return IODESystem::getVarName(ind);
            }
        }

        Vector<Real> getInitCond(Real angle, Real initHeight, Real velocity) const
        {
            Vector<Real> initCond(4);
            initCond[0] = 0;                            // initial x position
            initCond[1] = initHeight;                   // initial y position
            initCond[2] = velocity * cos(angle);        // initial x velocity (m/s)
            initCond[3] = velocity * sin(angle);        // initial y velocity (m/s)
            return initCond;
        }
    };

    /**************************************************************************
     * VACUUM MODEL (Ideal Parabolic Motion)
     **************************************************************************/

    class Projectile2DInVacuumODE : public IProjectile2DODE
    {
    public:
        static constexpr Real g = 9.81;  // gravitational acceleration (m/s²)

        int getDim() const override { return 4; }

        void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
        {
            dxdt[0] = x[2];
            dxdt[1] = x[3];
            dxdt[2] = 0;
            dxdt[3] = -g;

            // Stop motion when projectile hits the ground
            if (x[1] < 0)
            {
                dxdt[0] = 0;
                dxdt[1] = 0;
                dxdt[2] = 0;
                dxdt[3] = 0;
            }
        }

        // Analytical solutions for vacuum motion
        Real CalcRange(Real angle, Real initHeight, Real velocity) const
        {
            return (velocity * cos(angle) / g) * 
                   (velocity * sin(angle) + sqrt(velocity * velocity * sin(angle) * sin(angle) + 2 * g * initHeight));
        }

        Real TimeOfFlight(Real angle, Real initHeight, Real velocity) const
        {
            Real v0y = velocity * sin(angle);
            return (v0y + sqrt(v0y * v0y + 2 * g * initHeight)) / g;
        }

        Real MaxHeight(Real angle, Real initHeight, Real velocity) const
        {
            Real v0y = velocity * sin(angle);
            return initHeight + (v0y * v0y) / (2 * g);
        }
    };

    /**************************************************************************
     * CONSTANT AIR RESISTANCE MODEL (Quadratic Drag)
     **************************************************************************/

    class Projectile2DWithAirResistanceODE : public IProjectile2DODE
    {
        Real _dragCoefficient;

    public:
        Projectile2DWithAirResistanceODE(Real dragCoefficient = 0.1) 
            : _dragCoefficient(dragCoefficient) {}

        int getDim() const override { return 4; }

        void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
        {
            Real v = sqrt(x[2] * x[2] + x[3] * x[3]);  // speed

            dxdt[0] = x[2];
            dxdt[1] = x[3];
            dxdt[2] = -_dragCoefficient * v * x[2];
            dxdt[3] = -9.81 - _dragCoefficient * v * x[3];

            if (x[1] < 0)
            {
                dxdt[0] = 0;
                dxdt[1] = 0;
                dxdt[2] = 0;
                dxdt[3] = 0;
            }
        }
    };

    /**************************************************************************
     * ALTITUDE-DEPENDENT AIR DENSITY MODEL
     **************************************************************************/

    enum class AirDensityModel
    {
        Isothermal,    // Exponential decrease with altitude
        Adiabatic      // Standard atmosphere model
    };

    class Projectile2DChangingAirDensityODE : public IProjectile2DODE
    {
        Real _dragCoefficient0;              // drag coefficient at sea level
        Real _airDensity0 = 1.225;           // air density at sea level (kg/m³)
        AirDensityModel _airDensityModel = AirDensityModel::Isothermal;

    public:
        Projectile2DChangingAirDensityODE(Real dragCoefficient0 = 0.1) 
            : _dragCoefficient0(dragCoefficient0) {}

        Projectile2DChangingAirDensityODE(Real dragCoefficient0, AirDensityModel airDensityModel) 
            : _dragCoefficient0(dragCoefficient0), _airDensityModel(airDensityModel) {}

        int getDim() const override { return 4; }

        void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
        {
            Real v = sqrt(x[2] * x[2] + x[3] * x[3]);

            dxdt[0] = x[2];
            dxdt[1] = x[3];
            dxdt[2] = -dragCoefficient(x[3], v) * v * x[2];
            dxdt[3] = -9.81 - dragCoefficient(x[3], v) * v * x[3];

            if (x[1] < 0)
            {
                dxdt[0] = 0;
                dxdt[1] = 0;
                dxdt[2] = 0;
                dxdt[3] = 0;
            }
        }

        // Isothermal air density model (exponential decrease)
        Real airDensity(Real altitude) const
        {
            return 1.225 * exp(-altitude / 8000);  // scale height ~8000 m
        }

        // Adiabatic (standard atmosphere) model
        Real airDensityAdiabatic(Real altitude) const
        {
            return 1.225 * pow((1 - (0.0065 * altitude / 288.15)), 5.2561);
        }

        Real dragCoefficient(Real altitude, Real speed) const
        {
            Real density = (_airDensityModel == AirDensityModel::Adiabatic) 
                         ? airDensityAdiabatic(altitude) 
                         : airDensity(altitude);
            return _dragCoefficient0 * density / _airDensity0;
        }
    };

    /**************************************************************************
     * BASEBALL PHYSICS (Speed-Dependent Drag Coefficient)
     **************************************************************************/

    enum class BaseballType
    {
        Smooth,
        Normal,
        Rough,
        NormalApprox
    };

    class BaseballWithDragCoeffDependentOnSpeedODE : public IProjectile2DODE
    {
    private:
        BaseballType _baseballType;
        PolynomInterpRealFunc* _dragCoefficientSmooth;
        PolynomInterpRealFunc* _dragCoefficientNormal;
        PolynomInterpRealFunc* _dragCoefficientRough;

    public:
        BaseballWithDragCoeffDependentOnSpeedODE(BaseballType baseballType) 
            : _baseballType(baseballType) 
        {
            // Drag coefficient curves for different ball surfaces
            Vector<Real> smoothXVal{ 0.0, 50.0, 100.0, 150.0 };
            Vector<Real> smoothYVal{ 0.005, 0.0045, 0.004, 0.0035 }; 
            _dragCoefficientSmooth = new PolynomInterpRealFunc(smoothXVal, smoothYVal, 3);

            Vector<Real> normalXVal{ 0.0, 50.0, 100.0, 150.0 };
            Vector<Real> normalYVal{ 0.0039, 0.0035, 0.0032, 0.003 };
            _dragCoefficientNormal = new PolynomInterpRealFunc(normalXVal, normalYVal, 3);

            Vector<Real> roughXVal{ 0.0, 50.0, 100.0, 150.0 };
            Vector<Real> roughYVal{ 0.0045, 0.004, 0.0035, 0.003 };
            _dragCoefficientRough = new PolynomInterpRealFunc(roughXVal, roughYVal, 3);
        }

        BaseballWithDragCoeffDependentOnSpeedODE() 
            : BaseballWithDragCoeffDependentOnSpeedODE(BaseballType::NormalApprox) {}

        ~BaseballWithDragCoeffDependentOnSpeedODE()
        {
            delete _dragCoefficientSmooth;
            delete _dragCoefficientNormal;
            delete _dragCoefficientRough;
        }

        int getDim() const override { return 4; }

        void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
        {
            Real v = sqrt(x[2] * x[2] + x[3] * x[3]);

            dxdt[0] = x[2];
            dxdt[1] = x[3];
            dxdt[2] = -dragCoefficient(v) * v * x[2];
            dxdt[3] = -9.81 - dragCoefficient(v) * v * x[3];

            if (x[1] < 0)
            {
                dxdt[0] = 0;
                dxdt[1] = 0;
                dxdt[2] = 0;
                dxdt[3] = 0;
            }
        }

        Real dragCoefficient(Real v) const
        {
            switch (_baseballType)
            {
            case BaseballType::Smooth:       return (*_dragCoefficientSmooth)(v);
            case BaseballType::Normal:       return (*_dragCoefficientNormal)(v);
            case BaseballType::Rough:        return (*_dragCoefficientRough)(v);
            case BaseballType::NormalApprox: return dragCoefficientNormalApprox(v);
            default:                         return (*_dragCoefficientNormal)(v);
            }
        }

        // Analytical approximation for normal baseball
        Real dragCoefficientNormalApprox(Real v) const
        {
            Real vd = 35;
            Real delta = 5;
            return 0.0039 + 0.0058 / (1 + exp((v - vd) / delta));
        }
    };

    /**************************************************************************
     * TRAJECTORY STORAGE AND INTERPOLATION
     **************************************************************************/

    struct ProjectileTrajectory2D
    {
        Real _angle;            // launch angle (radians)
        Real _initHeight;       // initial height (m)
        Real _velocity;         // initial velocity (m/s)
        Real _timeOfFlight;     // total flight time (s)
        Real _range;            // horizontal range (m)

        Vector<Real> _tValues;  // time points
        Vector<Real> _xValues;  // x positions
        Vector<Real> _yValues;  // y positions
        Vector<Real> _vxValues; // x velocities
        Vector<Real> _vyValues; // y velocities

        ProjectileTrajectory2D() 
            : _angle(0.0), _initHeight(0.0), _velocity(0.0), _timeOfFlight(0.0), _range(0.0)
        {}

        ProjectileTrajectory2D(Real angle, Real initHeight, Real velocity)
            : _angle(angle), _initHeight(initHeight), _velocity(velocity), 
              _timeOfFlight(0.0), _range(0.0)
        {}
        
        // Interpolation functions
        LinearInterpRealFunc getXOfT() const
        {
            if (_tValues.size() != _xValues.size())
                throw std::runtime_error("t and x values must have the same size");
            return LinearInterpRealFunc(_tValues, _xValues);
        }

        LinearInterpRealFunc getYOfT() const
        {
            if (_tValues.size() != _yValues.size())
                throw std::runtime_error("t and y values must have the same size");
            return LinearInterpRealFunc(_tValues, _yValues);
        }

        LinearInterpRealFunc getVxOfT() const
        {
            if (_tValues.size() != _vxValues.size())
                throw std::runtime_error("t and vx values must have the same size");
            return LinearInterpRealFunc(_tValues, _vxValues);
        }

        LinearInterpRealFunc getVyOfT() const
        {
            if (_tValues.size() != _vyValues.size())
                throw std::runtime_error("t and vy values must have the same size");
            return LinearInterpRealFunc(_tValues, _vyValues);
        }

        LinearInterpRealFunc getYOfX() const
        {
            if (_xValues.size() != _yValues.size())
                throw std::runtime_error("x and y values must have the same size");
            return LinearInterpRealFunc(_xValues, _yValues);
        }

        PolynomInterpRealFunc getYOfXPoly() const
        {
            if (_xValues.size() != _yValues.size())
                throw std::runtime_error("x and y values must have the same size");
            return PolynomInterpRealFunc(_xValues, _yValues, 3);
        }
    };

    /**************************************************************************
     * PROJECTILE MOTION SOLVER
     **************************************************************************/

    class ProjectileMotionSolver2D
    {
    private:
        IProjectile2DODE& _projectileODESystem;
        ODESystemSolution _lastSol;

    public:
        ProjectileMotionSolver2D(IProjectile2DODE& odeSystem) 
            : _projectileODESystem(odeSystem), _lastSol(0.0, 1.0, 4) {}

        Real getRange(Real angle, Real initHeight, Real velocity) const
        {
            if (auto* vacuumODE = dynamic_cast<const Projectile2DInVacuumODE*>(&_projectileODESystem))
                return vacuumODE->CalcRange(angle, initHeight, velocity);
            return 0.0;  // numerical method needed for non-vacuum
        }

        Real getTimeOfFlight(Real angle, Real initHeight, Real velocity) const
        {
            if (auto* vacuumODE = dynamic_cast<const Projectile2DInVacuumODE*>(&_projectileODESystem))
                return vacuumODE->TimeOfFlight(angle, initHeight, velocity);
            return 0.0;
        }

        Real getTimeOfFlightVacuum(Real angle, Real initHeight, Real velocity) const
        {
            Real g = 9.81;
            Real v0y = velocity * sin(angle);
            return (v0y + sqrt(v0y * v0y + 2 * g * initHeight)) / g;
        }

        ODESystemSolution& getLastODESolution() { return _lastSol; }

        // Helper: Prepare results with ground-hit interpolation
        static ProjectileTrajectory2D prepareResults(
            Real angle, Real initHeight, Real velocity,
            const Vector<Real>& sol_t, const Vector<Real>& sol_x, const Vector<Real>& sol_y,
            const Vector<Real>& sol_vx, const Vector<Real>& sol_vy)
        {
            ProjectileTrajectory2D result(angle, initHeight, velocity);
            int numPoints = sol_t.size();

            // Check if projectile hit the ground
            if (sol_y.back() < 0)
            {
                // Find last point above ground
                int lastIndex = sol_y.size() - 1;
                while (lastIndex > 0 && sol_y[lastIndex] < 0)
                    lastIndex--;

                if (lastIndex < 0 || lastIndex + 1 >= (int)sol_y.size())
                {
                    // Fallback: use raw data
                    result._tValues = sol_t;
                    result._xValues = sol_x;
                    result._yValues = sol_y;
                    result._vxValues = sol_vx;
                    result._vyValues = sol_vy;
                    result._timeOfFlight = sol_t.back();
                    result._range = sol_x.back();
                    return result;
                }

                // Interpolate exact ground hit
                Real t1 = sol_t[lastIndex], y1 = sol_y[lastIndex];
                Real t2 = sol_t[lastIndex + 1], y2 = sol_y[lastIndex + 1];
                Real x1 = sol_x[lastIndex], x2 = sol_x[lastIndex + 1];

                Real alpha = (0 - y1) / (y2 - y1);
                Real tHitGround = t1 + alpha * (t2 - t1);
                Real xHitGround = x1 + alpha * (x2 - x1);
                Real vxHitGround = sol_vx[lastIndex] + alpha * (sol_vx[lastIndex + 1] - sol_vx[lastIndex]);
                Real vyHitGround = sol_vy[lastIndex] + alpha * (sol_vy[lastIndex + 1] - sol_vy[lastIndex]);

                result._timeOfFlight = tHitGround;
                result._range = xHitGround;

                result._tValues.Resize(numPoints);
                result._xValues.Resize(numPoints);
                result._yValues.Resize(numPoints);
                result._vxValues.Resize(numPoints);
                result._vyValues.Resize(numPoints);

                // Copy up to ground hit
                for (int i = 0; i <= lastIndex; i++)
                {
                    result._tValues[i] = sol_t[i];
                    result._xValues[i] = sol_x[i];
                    result._yValues[i] = sol_y[i];
                    result._vxValues[i] = sol_vx[i];
                    result._vyValues[i] = sol_vy[i];
                }

                // Add exact hit point
                result._tValues[lastIndex + 1] = tHitGround;
                result._xValues[lastIndex + 1] = xHitGround;
                result._yValues[lastIndex + 1] = 0;
                result._vxValues[lastIndex + 1] = vxHitGround;
                result._vyValues[lastIndex + 1] = vyHitGround;

                // Fill rest with zeros
                for (int i = lastIndex + 2; i < numPoints; i++)
                {
                    result._tValues[i] = sol_t[i];
                    result._xValues[i] = xHitGround;
                    result._yValues[i] = 0;
                    result._vxValues[i] = 0;
                    result._vyValues[i] = 0;
                }
            }
            else
            {
                result._tValues = sol_t;
                result._xValues = sol_x;
                result._yValues = sol_y;
                result._vxValues = sol_vx;
                result._vyValues = sol_vy;
                result._timeOfFlight = sol_t.back();
                result._range = sol_x.back();
            }

            return result;
        }

        // Multi-angle solving
        Vector<ProjectileTrajectory2D> solveForAnglesEuler(
            const Vector<Real>& angles, Real initHeight, Real velocity, Real tMax, Real dt)
        {
            Vector<ProjectileTrajectory2D> results;
            for (Real angle : angles)
                results.push_back(solveEuler(angle, initHeight, velocity, tMax, dt));
            return results;
        }

        Vector<ProjectileTrajectory2D> solveForAnglesEuler(
            const Vector<Real>& angles, Real initHeight, Real velocity, Real dt)
        {
            // Find max time of flight
            Real tMax = 0.0;
            for (Real angle : angles)
            {
                Real t = getTimeOfFlightVacuum(angle, initHeight, velocity);
                if (t > tMax) tMax = t;
            }
            tMax *= 1.1;  // margin

            return solveForAnglesEuler(angles, initHeight, velocity, tMax, dt);
        }

        // Euler method solver
        ProjectileTrajectory2D solveEuler(
            Real angle, Real initHeight, Real velocity, Real tMax, Real dt)
        {
            Vector<Real> initCond = _projectileODESystem.getInitCond(angle, initHeight, velocity);
            int numSteps = (int)(tMax / dt);

            ODESystemFixedStepSolver odeSolver(_projectileODESystem, StepCalculators::EulerStepCalc);
            ODESystemSolution sol = odeSolver.integrate(initCond, 0.0, tMax, numSteps);
            _lastSol = sol;

            Vector<Real> sol_t{ sol.getTValues() }, sol_x{ sol.getXValues(0) };
            Vector<Real> sol_y{ sol.getXValues(1) }, sol_vx{ sol.getXValues(2) };
            Vector<Real> sol_vy{ sol.getXValues(3) };

            return prepareResults(angle, initHeight, velocity, sol_t, sol_x, sol_y, sol_vx, sol_vy);
        }

        // RK4 method (placeholder)
        ProjectileTrajectory2D solveRK4(
            Real angle, Real initHeight, Real velocity, Real tMax, Real dt)
        {
            ProjectileTrajectory2D result;
            result._angle = angle;
            result._initHeight = initHeight;
            result._velocity = velocity;
            return result;
        }

        // RK5 Cash-Karp adaptive method
        ProjectileTrajectory2D solveRK5(
            Real angle, Real initHeight, Real velocity, Real tMax, Real eps = 1e-7)
        {
            Vector<Real> initCond = _projectileODESystem.getInitCond(angle, initHeight, velocity);

            ODESystemSolver<RK5_CashKarp_Stepper> odeSolver(_projectileODESystem);
            ODESystemSolution sol = odeSolver.integrate(initCond, 0.0, tMax, 0.01, eps, 0.01);
            _lastSol = sol;

            Vector<Real> sol_t{ sol.getTValues() }, sol_x{ sol.getXValues(0) };
            Vector<Real> sol_y{ sol.getXValues(1) }, sol_vx{ sol.getXValues(2) };
            Vector<Real> sol_vy{ sol.getXValues(3) };

            return prepareResults(angle, initHeight, velocity, sol_t, sol_x, sol_y, sol_vx, sol_vy);
        }
    };

} // namespace Projectile

#endif // PROJECTILES_2D_SELF_CONTAINED_H
