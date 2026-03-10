/******************************************************************************
 * Projectiles2D.h - Self-contained Projectile Motion Physics
 * ============================================================================
 * 
 * This header provides all necessary physics for 2D projectile simulation,
 * making the example self-contained without MPL dependencies.
 * 
 * Contents:
 *   - IProjectile2DODE: Base interface (extends IODESystemWithEvents for
 *                        automatic ground-hit detection via zero-crossing)
 *   - Projectile2DInVacuumODE: Ideal parabolic motion (analytical solutions)
 *   - Projectile2DWithAirResistanceODE: Quadratic drag model
 *   - Projectile2DChangingAirDensityODE: Altitude-dependent drag
 *   - BaseballWithDragCoeffDependentOnSpeedODE: Speed-dependent drag
 *   - ProjectileTrajectory2D: Trajectory storage and interpolation
 *   - ProjectileMotionSolver2D: ODE solver with event-detection ground hit
 * 
 * Physics: F = mg - kv² (quadratic drag model)
 * 
 * Ground Detection:
 *   Uses MML's event detection (IODESystemWithEvents + integrateWithEvents)
 *   to precisely locate ground impact via zero-crossing of y(t) to 1e-12
 *   tolerance. No more crude "if(y<0) zero derivatives" hack!
 * 
 *****************************************************************************/

#ifndef PROJECTILES_2D_SELF_CONTAINED_H
#define PROJECTILES_2D_SELF_CONTAINED_H

#include "MMLBase.h"
#include "interfaces/IODESystem.h"
#include "base/InterpolatedFunction.h"
#include "mml/algorithms/ODESolvers/ODEAdaptiveIntegrator.h"
#include "mml/algorithms/ODESolvers/ODEFixedStepIntegrators.h"
#include "mml/algorithms/ODESolvers/ODESystemStepCalculators.h"

namespace Projectile
{
    using namespace MML;

    /**************************************************************************
     * BASE INTERFACE
     * 
     * Extends IODESystemWithEvents so all projectile models get automatic
     * ground-hit detection. Event 0 monitors y position crossing zero
     * in the decreasing direction and stops integration on contact.
     **************************************************************************/

    class IProjectile2DODE : public IODESystemWithEvents
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

        // ---- Event detection: ground contact (y = 0) ----
        int getNumEvents() const override { return 1; }

        Real eventFunction(int /*eventIndex*/, Real /*t*/, const Vector<Real>& x) const override
        {
            return x[1];  // y position — zero means ground contact
        }

        EventDirection getEventDirection(int /*eventIndex*/) const override
        {
            return EventDirection::Decreasing;  // falling toward ground
        }

        EventAction getEventAction(int /*eventIndex*/) const override
        {
            return EventAction::Stop;  // terminate integration on ground hit
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
     * 
     * Primary solver: solveAdaptive() — uses DormandPrince5 with event
     * detection for precise ground-hit location (1e-12 tolerance).
     * Legacy: solveEuler() kept for accuracy comparison.
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

        //---------------------------------------------------------------------
        // PRIMARY: Adaptive solver with event-detected ground impact
        //---------------------------------------------------------------------

        /// @brief Solve using DormandPrince5 adaptive integrator with ground-hit event detection.
        ///
        /// Integration stops automatically when y(t) crosses zero (decreasing).
        /// The ground-hit time and position are located to ~1e-12 precision
        /// using the Illinois root-finding method on the stepper's dense output.
        ///
        /// @param angle        Launch angle (radians)
        /// @param initHeight   Initial height (m)
        /// @param velocity     Launch speed (m/s)
        /// @param outputDt     Output interval — controls trajectory resolution (s)
        /// @param eps          ODE solver tolerance (default 1e-10)
        /// @return             Trajectory truncated precisely at ground impact
        ProjectileTrajectory2D solveAdaptive(
            Real angle, Real initHeight, Real velocity, 
            Real outputDt = 0.01, Real eps = 1e-10)
        {
            Vector<Real> initCond = _projectileODESystem.getInitCond(angle, initHeight, velocity);

            // Generous upper bound — event detection will stop us at ground hit
            Real tMax = 2.0 * getTimeOfFlightVacuum(angle, initHeight, velocity) + 10.0;

            DormandPrince5Integrator integrator(_projectileODESystem);
            auto eventResult = integrator.integrateWithEvents(
                _projectileODESystem, initCond, 0.0, tMax, outputDt, eps);

            return extractTrajectory(eventResult, angle, initHeight, velocity);
        }

        /// @brief Solve multiple launch angles using adaptive event detection.
        Vector<ProjectileTrajectory2D> solveForAngles(
            const Vector<Real>& angles, Real initHeight, Real velocity,
            Real outputDt = 0.01, Real eps = 1e-10)
        {
            Vector<ProjectileTrajectory2D> results;
            for (Real angle : angles)
                results.push_back(solveAdaptive(angle, initHeight, velocity, outputDt, eps));
            return results;
        }

        //---------------------------------------------------------------------
        // LEGACY: Euler (1st order) — kept for accuracy comparison
        //---------------------------------------------------------------------

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

            return prepareResultsLegacy(angle, initHeight, velocity, sol_t, sol_x, sol_y, sol_vx, sol_vy);
        }

        Vector<ProjectileTrajectory2D> solveForAnglesEuler(
            const Vector<Real>& angles, Real initHeight, Real velocity, Real dt)
        {
            Real tMax = 0.0;
            for (Real angle : angles)
            {
                Real t = getTimeOfFlightVacuum(angle, initHeight, velocity);
                if (t > tMax) tMax = t;
            }
            tMax *= 1.1;

            Vector<ProjectileTrajectory2D> results;
            for (Real angle : angles)
                results.push_back(solveEuler(angle, initHeight, velocity, tMax, dt));
            return results;
        }

    private:
        /// @brief Extract trajectory from event result (clean, no post-ground zeros)
        static ProjectileTrajectory2D extractTrajectory(
            const DormandPrince5Integrator::EventResult& eventResult,
            Real angle, Real initHeight, Real velocity)
        {
            ProjectileTrajectory2D result(angle, initHeight, velocity);
            const auto& sol = eventResult.solution;

            // Count valid output points (up to event termination).
            // The solution is pre-allocated for full tMax, but event detection
            // stops early — unfilled entries remain as 0. Detect them by
            // checking that t values remain strictly increasing.
            int numPoints = 0;
            auto tVals = sol.getTValues();
            for (int i = 0; i < (int)tVals.size(); i++)
            {
                if (tVals[i] > eventResult.finalTime + 1e-10)
                    break;
                // After the first point, detect unfilled zeros (t not increasing)
                if (i > 0 && tVals[i] <= tVals[i - 1])
                    break;
                numPoints++;
            }

            // Include the precise ground-hit point if event detected
            bool addEventPoint = eventResult.terminatedByEvent;
            int totalPoints = numPoints + (addEventPoint ? 1 : 0);

            result._tValues.Resize(totalPoints);
            result._xValues.Resize(totalPoints);
            result._yValues.Resize(totalPoints);
            result._vxValues.Resize(totalPoints);
            result._vyValues.Resize(totalPoints);

            for (int i = 0; i < numPoints; i++)
            {
                result._tValues[i]  = sol.getTValues()[i];
                result._xValues[i]  = sol.getXValues(0)[i];
                result._yValues[i]  = sol.getXValues(1)[i];
                result._vxValues[i] = sol.getXValues(2)[i];
                result._vyValues[i] = sol.getXValues(3)[i];
            }

            if (addEventPoint)
            {
                int idx = numPoints;
                result._tValues[idx]  = eventResult.finalTime;
                result._xValues[idx]  = eventResult.finalState[0];
                result._yValues[idx]  = 0.0;  // exact ground
                result._vxValues[idx] = eventResult.finalState[2];
                result._vyValues[idx] = eventResult.finalState[3];

                result._timeOfFlight = eventResult.finalTime;
                result._range = eventResult.finalState[0];
            }
            else
            {
                result._timeOfFlight = sol.getTValues()[numPoints - 1];
                result._range = sol.getXValues(0)[numPoints - 1];
            }

            return result;
        }

        /// @brief Legacy ground-hit interpolation for Euler solver
        static ProjectileTrajectory2D prepareResultsLegacy(
            Real angle, Real initHeight, Real velocity,
            const Vector<Real>& sol_t, const Vector<Real>& sol_x, const Vector<Real>& sol_y,
            const Vector<Real>& sol_vx, const Vector<Real>& sol_vy)
        {
            ProjectileTrajectory2D result(angle, initHeight, velocity);
            int numPoints = sol_t.size();

            if (sol_y.back() < 0)
            {
                int lastIndex = sol_y.size() - 1;
                while (lastIndex > 0 && sol_y[lastIndex] < 0) lastIndex--;

                if (lastIndex < 0 || lastIndex + 1 >= (int)sol_y.size())
                {
                    result._tValues = sol_t;  result._xValues = sol_x;
                    result._yValues = sol_y;  result._vxValues = sol_vx;
                    result._vyValues = sol_vy;
                    result._timeOfFlight = sol_t.back();
                    result._range = sol_x.back();
                    return result;
                }

                Real alpha = (0 - sol_y[lastIndex]) / (sol_y[lastIndex + 1] - sol_y[lastIndex]);
                Real tHit = sol_t[lastIndex] + alpha * (sol_t[lastIndex + 1] - sol_t[lastIndex]);
                Real xHit = sol_x[lastIndex] + alpha * (sol_x[lastIndex + 1] - sol_x[lastIndex]);

                result._timeOfFlight = tHit;
                result._range = xHit;

                // Store only up to ground hit + exact hit point
                int trimmedSize = lastIndex + 2;
                result._tValues.Resize(trimmedSize);
                result._xValues.Resize(trimmedSize);
                result._yValues.Resize(trimmedSize);
                result._vxValues.Resize(trimmedSize);
                result._vyValues.Resize(trimmedSize);

                for (int i = 0; i <= lastIndex; i++)
                {
                    result._tValues[i] = sol_t[i];
                    result._xValues[i] = sol_x[i];
                    result._yValues[i] = sol_y[i];
                    result._vxValues[i] = sol_vx[i];
                    result._vyValues[i] = sol_vy[i];
                }
                result._tValues[lastIndex + 1] = tHit;
                result._xValues[lastIndex + 1] = xHit;
                result._yValues[lastIndex + 1] = 0.0;
                result._vxValues[lastIndex + 1] = sol_vx[lastIndex];
                result._vyValues[lastIndex + 1] = sol_vy[lastIndex];
            }
            else
            {
                result._tValues = sol_t;  result._xValues = sol_x;
                result._yValues = sol_y;  result._vxValues = sol_vx;
                result._vyValues = sol_vy;
                result._timeOfFlight = sol_t.back();
                result._range = sol_x.back();
            }
            return result;
        }
    };

} // namespace Projectile

#endif // PROJECTILES_2D_SELF_CONTAINED_H
