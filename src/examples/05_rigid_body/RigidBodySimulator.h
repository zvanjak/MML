///////////////////////////////////////////////////////////////////////////////////////////
/// @file RigidBodySimulator.h
/// @brief Main simulator engine - Example 05: Rigid Body Collision Simulator
/// @details Self-contained header with simulation engine:
///          - SimulationConfig struct
///          - SimulationSnapshot struct
///          - RigidBodySimulator class
///
/// SIMULATION LOOP:
/// 1. Integrate rigid body dynamics (Euler equations) using RK4
/// 2. Detect collisions (box-box, box-wall, sphere-sphere, sphere-wall)
/// 3. Resolve collisions (apply impulses)
/// 4. Re-normalize quaternions
/// 5. Optionally log state for visualization
///
/// This file is extracted from the MPL library for use in Example 05.
/// Namespace changed from MPL to RigidBodySim to avoid conflicts.
///
/// @author MinimalMathLibrary - Example 05
/// @date February 2026
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef RIGID_BODY_SIMULATOR_H
#define RIGID_BODY_SIMULATOR_H

#include "RigidBodyCore.h"
#include "RigidBodyDynamics.h"
#include "RigidBodyCollision.h"
#include "mml/algorithms/ODESolvers.h"

#include <vector>
#include <functional>
#include <memory>
#include <algorithm>

namespace RigidBodySim
{
    using namespace MML;

    // ========================= Simulation Configuration =========================

    /// @brief Simulation configuration
    struct SimulationConfig
    {
        Real timeStep = 0.001;                  ///< Integration time step (seconds)
        Real totalTime = 10.0;                  ///< Total simulation time (seconds)
        Real containerHalfSize = 5.0;           ///< Half-size of cubic container
        Real coeffRestitution = 1.0;            ///< 1.0 = perfectly elastic
        int maxCollisionIterations = 10;        ///< Max iterations per time step
        bool normalizeQuaternions = true;       ///< Re-normalize quaternions each step
    };

    // ========================= Simulation Snapshot =========================

    /// @brief Snapshot of simulation state at a time point
    struct SimulationSnapshot
    {
        Real time;
        std::vector<RigidBodyState> bodyStates;
        Real totalEnergy;
        Vec3Cart totalMomentum;
        Vec3Cart totalAngularMomentum;
    };

    // ========================= Rigid Body Simulator =========================

    /// @brief Main rigid body simulator engine
    /// @details Simulates multiple rigid bodies with collisions in a cubic container
    class RigidBodySimulator
    {
    public:
        /// @brief Callback for simulation progress (time, bodies)
        using BodyContainer = std::vector<std::unique_ptr<RigidBody>>;
        using ProgressCallback = std::function<void(Real, const BodyContainer&)>;

    private:
        BodyContainer _bodies;
        SimulationConfig _config;
        std::vector<SimulationSnapshot> _history;
        ProgressCallback _progressCallback;
        
        // Step calculators
        RungeKutta4_StepCalculator _rk4;

    public:
        // ========================= Setup =========================

        /// @brief Constructor with default configuration
        RigidBodySimulator() = default;

        /// @brief Constructor with custom configuration
        explicit RigidBodySimulator(const SimulationConfig& config)
            : _config(config) {}

        /// @brief Add a rigid body to the simulation (takes ownership)
        void AddBody(std::unique_ptr<RigidBody> body)
        {
            _bodies.push_back(std::move(body));
        }

        /// @brief Add a rigid body by cloning (convenience for concrete types)
        template<typename T, typename = std::enable_if_t<std::is_base_of_v<RigidBody, T>>>
        void AddBody(const T& body)
        {
            _bodies.push_back(body.Clone());
        }

        /// @brief Clear all bodies
        void ClearBodies()
        {
            _bodies.clear();
        }

        /// @brief Set simulation configuration
        void SetConfig(const SimulationConfig& config)
        {
            _config = config;
        }

        /// @brief Set progress callback (called each timestep)
        void SetProgressCallback(ProgressCallback callback)
        {
            _progressCallback = std::move(callback);
        }

        // ========================= Access =========================

        /// @brief Get number of bodies
        size_t NumBodies() const { return _bodies.size(); }

        /// @brief Get body by index (pointer access)
        RigidBody& GetBody(size_t i) { return *_bodies[i]; }
        const RigidBody& GetBody(size_t i) const { return *_bodies[i]; }

        /// @brief Get body pointer by index
        RigidBody* GetBodyPtr(size_t i) { return _bodies[i].get(); }
        const RigidBody* GetBodyPtr(size_t i) const { return _bodies[i].get(); }

        /// @brief Get all bodies container (for advanced access)
        const BodyContainer& GetBodies() const { return _bodies; }

        /// @brief Get simulation history
        const std::vector<SimulationSnapshot>& GetHistory() const { return _history; }

        /// @brief Get current configuration
        const SimulationConfig& GetConfig() const { return _config; }

        // ========================= Simulation =========================

        /// @brief Run the full simulation
        void Run()
        {
            _history.clear();
            
            Real t = 0.0;
            int numSteps = static_cast<int>(_config.totalTime / _config.timeStep);
            int recordInterval = std::max(1, static_cast<int>(0.01 / _config.timeStep)); // Record every 10ms
            
            // Record initial state
            RecordSnapshot(t);
            
            for (int step = 0; step < numSteps; step++)
            {
                // 1. Integrate dynamics
                IntegrateStep(_config.timeStep);
                
                // 2. Detect and resolve collisions
                ResolveCollisions();
                
                // 3. Normalize quaternions
                if (_config.normalizeQuaternions)
                {
                    for (auto& body : _bodies)
                        body->NormalizeOrientation();
                }
                
                t += _config.timeStep;
                
                // Record state periodically
                if ((step + 1) % recordInterval == 0)
                    RecordSnapshot(t);
                
                // Progress callback
                if (_progressCallback)
                    _progressCallback(t, _bodies);
            }
            
            // Record final state if not just recorded
            if (numSteps % recordInterval != 0)
                RecordSnapshot(t);
        }

        /// @brief Take a single time step (for manual stepping)
        void Step()
        {
            IntegrateStep(_config.timeStep);
            ResolveCollisions();
            
            if (_config.normalizeQuaternions)
            {
                for (auto& body : _bodies)
                    body->NormalizeOrientation();
            }
        }

        // ========================= Conservation =========================

        /// @brief Compute total kinetic energy
        Real TotalKineticEnergy() const
        {
            Real total = 0.0;
            for (const auto& body : _bodies)
                total += body->KineticEnergy();
            return total;
        }

        /// @brief Compute total linear momentum
        Vec3Cart TotalLinearMomentum() const
        {
            Vec3Cart total(0, 0, 0);
            for (const auto& body : _bodies)
                total = total + body->LinearMomentum();
            return total;
        }

        /// @brief Compute total angular momentum (about origin)
        Vec3Cart TotalAngularMomentum() const
        {
            Vec3Cart total(0, 0, 0);
            for (const auto& body : _bodies)
            {
                // L = L_spin + r × p
                Vec3Cart L_spin = body->AngularMomentum();
                Vec3Cart r_cross_p = VectorProduct(body->Position(), body->LinearMomentum());
                total = total + L_spin + r_cross_p;
            }
            return total;
        }

    private:
        /// @brief Integrate all bodies forward by dt
        void IntegrateStep(Real dt)
        {
            if (_bodies.size() == 1)
            {
                IntegrateSingleBody(dt);
            }
            else if (_bodies.size() == 2)
            {
                IntegrateTwoBodies(dt);
            }
            else
            {
                // For more bodies, integrate each independently
                for (auto& body : _bodies)
                {
                    SingleBodyODESystem ode(*body);
                    Vector<Real> state(13), dxdt(13), x_out(13), x_err(13);
                    PackBodyState(state, body->_state, 0);
                    
                    ode.derivs(0, state, dxdt);
                    _rk4.calcStep(ode, 0, state, dxdt, dt, x_out, x_err);
                    
                    body->_state = ExtractBodyState(x_out, 0);
                }
            }
        }

        void IntegrateSingleBody(Real dt)
        {
            SingleBodyODESystem ode(*_bodies[0]);
            Vector<Real> state(13), dxdt(13), x_out(13), x_err(13);
            PackBodyState(state, _bodies[0]->_state, 0);
            
            ode.derivs(0, state, dxdt);
            _rk4.calcStep(ode, 0, state, dxdt, dt, x_out, x_err);
            
            _bodies[0]->_state = ExtractBodyState(x_out, 0);
        }

        void IntegrateTwoBodies(Real dt)
        {
            TwoBodyODESystem ode(*_bodies[0], *_bodies[1]);
            Vector<Real> state(26), dxdt(26), x_out(26), x_err(26);
            
            PackBodyState(state, _bodies[0]->_state, 0);
            PackBodyState(state, _bodies[1]->_state, 1);  // Body index 1, not offset 13
            
            ode.derivs(0, state, dxdt);
            _rk4.calcStep(ode, 0, state, dxdt, dt, x_out, x_err);
            
            _bodies[0]->_state = ExtractBodyState(x_out, 0);
            _bodies[1]->_state = ExtractBodyState(x_out, 1);  // Body index 1
        }

        /// @brief Detect and resolve all collisions
        void ResolveCollisions()
        {
            CollisionParams params;
            params.coeffRestitution = _config.coeffRestitution;
            
            for (int iter = 0; iter < _config.maxCollisionIterations; iter++)
            {
                bool anyCollision = false;
                
                // Body-wall collisions (unified dispatcher handles box vs sphere)
                for (auto& body : _bodies)
                {
                    auto wallCollisions = CollisionDetector::DetectBodyWallCollisions(
                        *body, _config.containerHalfSize);
                    
                    for (const auto& collision : wallCollisions)
                    {
                        if (CollisionResponse::ApplyWallImpulse(*body, collision, params))
                            anyCollision = true;
                    }
                }
                
                // Body-body collisions (unified dispatcher handles all shape combinations)
                for (size_t i = 0; i < _bodies.size(); i++)
                {
                    for (size_t j = i + 1; j < _bodies.size(); j++)
                    {
                        auto collision = CollisionDetector::DetectBodyBodyCollision(
                            *_bodies[i], *_bodies[j]);
                        
                        if (CollisionResponse::ApplyBodyBodyImpulse(
                            *_bodies[i], *_bodies[j], collision, params))
                        {
                            anyCollision = true;
                        }
                    }
                }
                
                if (!anyCollision)
                    break;
            }
        }

        /// @brief Record current state to history
        void RecordSnapshot(Real t)
        {
            SimulationSnapshot snapshot;
            snapshot.time = t;
            snapshot.totalEnergy = TotalKineticEnergy();
            snapshot.totalMomentum = TotalLinearMomentum();
            snapshot.totalAngularMomentum = TotalAngularMomentum();
            
            for (const auto& body : _bodies)
                snapshot.bodyStates.push_back(body->_state);
            
            _history.push_back(snapshot);
        }
    };

} // namespace RigidBodySim

#endif // RIGID_BODY_SIMULATOR_H
