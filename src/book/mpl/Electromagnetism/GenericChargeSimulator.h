///////////////////////////////////////////////////////////////////////////////////////////
/// @file GenericChargeSimulator.h
/// @brief Generic charge distribution simulator using MML's IBody geometry interface
/// @details Simulates repulsive Coulomb forces between charged particles confined to
///          ANY geometry that implements IBody (Sphere3D, Cube3D, Cylinder3D, etc.)
///
/// PHYSICS:
/// - Coulomb force: F = k * q1 * q2 / r^2
/// - Simple Euler integration with velocity damping
/// - Rejection sampling for initial particle placement
/// - Binary search for boundary constraint enforcement
///
/// @author MinimalMathLibrary
/// @date February 2026
///////////////////////////////////////////////////////////////////////////////////////////

#ifndef MPL_GENERIC_CHARGE_SIMULATOR_H
#define MPL_GENERIC_CHARGE_SIMULATOR_H

#include "MMLBase.h"
#include "base/Vector/VectorTypes.h"
#include "base/Geometry/Geometry3DBodies.h"

#include <random>
#include <memory>
#include <vector>
#include <iostream>
#include <iomanip>

namespace MPL
{
    using namespace MML;

    /// @brief Charged particle with position, velocity, charge and mass
    struct ChargedParticle
    {
        Point3Cartesian pos;
        Vector3Cartesian velocity;
        Real charge;
        Real mass;
        
        ChargedParticle() : pos(0,0,0), velocity(0,0,0), charge(1.0), mass(1.0) {}
        ChargedParticle(const Point3Cartesian& p, Real q, Real m)
            : pos(p), velocity(0,0,0), charge(q), mass(m) {}
    };

    /// @brief Generic charge distribution simulator for any IBody geometry
    /// 
    /// Simulates charged particles repelling each other via Coulomb forces
    /// while constrained to remain inside a specified geometry.
    ///
    /// @example
    /// @code
    /// auto sphere = std::make_unique<Sphere3D>(100.0);
    /// GenericChargeSimulator sim(std::move(sphere));
    /// sim.initializeRandomParticles(100, 1.0, 1.0);
    /// auto history = sim.runToEquilibrium(2000);
    /// @endcode
    class GenericChargeSimulator
    {
    public:
        using Particle = ChargedParticle;

    private:
        std::unique_ptr<IBody> _geometry;  ///< The confining geometry
        std::vector<Particle> _particles;
        
        // Physical constants (using normalized units)
        Real _coulombConstant = 1.0;  ///< k in F = k*q1*q2/r^2
        Real _damping = 0.9;          ///< Velocity damping factor per step
        Real _dt = 0.01;              ///< Time step
        
    public:
        /// @brief Constructor takes ownership of any IBody geometry
        template<typename GeometryType>
        GenericChargeSimulator(std::unique_ptr<GeometryType> geometry)
            : _geometry(std::move(geometry))
        {
        }

        /// @brief Set simulation time step
        void setTimeStep(Real dt) { _dt = dt; }
        
        /// @brief Set velocity damping factor (0-1, applied each step)
        void setDamping(Real damping) { _damping = damping; }
        
        /// @brief Set Coulomb constant (k in F = k*q1*q2/r^2)
        void setCoulombConstant(Real k) { _coulombConstant = k; }

        /// @brief Initialize particles randomly within the geometry using rejection sampling
        /// @param numParticles Number of particles to create
        /// @param charge Charge of each particle
        /// @param mass Mass of each particle
        /// @param seed Random seed for reproducibility
        void initializeRandomParticles(int numParticles, Real charge, Real mass, unsigned int seed = 42)
        {
            std::mt19937 rng(seed);
            
            // Get bounding box for sampling region
            Box3D bbox = _geometry->GetBoundingBox();
            std::uniform_real_distribution<Real> xDist(bbox.MinX(), bbox.MaxX());
            std::uniform_real_distribution<Real> yDist(bbox.MinY(), bbox.MaxY());
            std::uniform_real_distribution<Real> zDist(bbox.MinZ(), bbox.MaxZ());

            _particles.clear();
            _particles.reserve(numParticles);

            // Rejection sampling: generate points in bounding box, keep only those inside geometry
            int count = 0;
            int maxAttempts = numParticles * 100;  // Safety limit
            int attempts = 0;
            
            while (count < numParticles && attempts < maxAttempts)
            {
                Pnt3Cart candidate(xDist(rng), yDist(rng), zDist(rng));
                attempts++;
                
                if (_geometry->IsInside(candidate))
                {
                    Particle p;
                    p.pos = Point3Cartesian(candidate.X(), candidate.Y(), candidate.Z());
                    p.velocity = Vector3Cartesian(0, 0, 0);
                    p.charge = charge;
                    p.mass = mass;
                    _particles.push_back(p);
                    count++;
                }
            }
            
            if (count < numParticles)
            {
                std::cerr << "Warning: Only placed " << count << " of " << numParticles 
                          << " particles (rejection sampling limit reached)" << std::endl;
            }
        }

        /// @brief Compute Coulomb force on particle i from all other particles
        Vector3Cartesian computeCoulombForce(int i) const
        {
            Vector3Cartesian force(0, 0, 0);
            const Particle& pi = _particles[i];

            for (size_t j = 0; j < _particles.size(); ++j)
            {
                if (i == static_cast<int>(j)) continue;

                const Particle& pj = _particles[j];
                
                // Vector from j to i
                Vector3Cartesian r_ij(
                    pi.pos.X() - pj.pos.X(),
                    pi.pos.Y() - pj.pos.Y(),
                    pi.pos.Z() - pj.pos.Z()
                );

                Real dist = r_ij.NormL2();
                if (dist < 1e-6) dist = 1e-6;  // Avoid singularity

                // Coulomb force: F = k * q1 * q2 / r^2, direction = r_hat
                Real forceMag = _coulombConstant * pi.charge * pj.charge / (dist * dist);
                
                // Add force (positive = repulsive for same-sign charges)
                force = force + r_ij * (forceMag / dist);
            }

            return force;
        }

        /// @brief Constrain particle to stay within geometry
        /// Uses binary search to find boundary when particle escapes
        void constrainToGeometry(Particle& p) const
        {
            Pnt3Cart pos(p.pos.X(), p.pos.Y(), p.pos.Z());
            
            if (!_geometry->IsInside(pos))
            {
                // Get geometry center
                Pnt3Cart center = _geometry->GetCenter();
                
                // Direction from particle toward center
                Real dx = center.X() - p.pos.X();
                Real dy = center.Y() - p.pos.Y();
                Real dz = center.Z() - p.pos.Z();
                Real dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                
                if (dist < 1e-10) {
                    // Particle at center but outside? Shouldn't happen, but handle gracefully
                    p.pos = Point3Cartesian(center.X(), center.Y(), center.Z());
                    p.velocity = Vector3Cartesian(0, 0, 0);
                    return;
                }
                
                // Normalize direction
                dx /= dist;
                dy /= dist;
                dz /= dist;
                
                // Binary search to find boundary point
                Real lo = 0.0, hi = dist;
                for (int iter = 0; iter < 20; ++iter)  // ~1e-6 precision
                {
                    Real mid = (lo + hi) / 2.0;
                    Pnt3Cart test(p.pos.X() + dx * mid, p.pos.Y() + dy * mid, p.pos.Z() + dz * mid);
                    if (_geometry->IsInside(test))
                        hi = mid;
                    else
                        lo = mid;
                }
                
                // Place particle just inside boundary
                Real pushDist = hi + 1e-6;
                p.pos = Point3Cartesian(
                    p.pos.X() + dx * pushDist,
                    p.pos.Y() + dy * pushDist,
                    p.pos.Z() + dz * pushDist
                );
                
                // Damp velocity when hitting boundary
                p.velocity = p.velocity * 0.5;
            }
        }

        /// @brief Perform one simulation step (Euler integration)
        void step()
        {
            // Compute forces
            std::vector<Vector3Cartesian> forces(_particles.size());
            for (size_t i = 0; i < _particles.size(); ++i)
            {
                forces[i] = computeCoulombForce(static_cast<int>(i));
            }

            // Update velocities and positions
            for (size_t i = 0; i < _particles.size(); ++i)
            {
                Particle& p = _particles[i];
                
                // a = F/m
                Vector3Cartesian accel = forces[i] * (1.0 / p.mass);
                
                // Update velocity with damping
                p.velocity = (p.velocity + accel * _dt) * _damping;
                
                // Update position
                p.pos = Point3Cartesian(
                    p.pos.X() + p.velocity.X() * _dt,
                    p.pos.Y() + p.velocity.Y() * _dt,
                    p.pos.Z() + p.velocity.Z() * _dt
                );
                
                // Apply geometry constraint
                constrainToGeometry(p);
            }
        }

        /// @brief Compute total kinetic energy
        Real totalKineticEnergy() const
        {
            Real ke = 0;
            for (const auto& p : _particles)
            {
                Real v2 = p.velocity.X()*p.velocity.X() + 
                          p.velocity.Y()*p.velocity.Y() + 
                          p.velocity.Z()*p.velocity.Z();
                ke += 0.5 * p.mass * v2;
            }
            return ke;
        }

        /// @brief Compute total potential energy
        Real totalPotentialEnergy() const
        {
            Real pe = 0;
            for (size_t i = 0; i < _particles.size(); ++i)
            {
                for (size_t j = i + 1; j < _particles.size(); ++j)
                {
                    Real dx = _particles[i].pos.X() - _particles[j].pos.X();
                    Real dy = _particles[i].pos.Y() - _particles[j].pos.Y();
                    Real dz = _particles[i].pos.Z() - _particles[j].pos.Z();
                    Real dist = std::sqrt(dx*dx + dy*dy + dz*dz);
                    if (dist < 1e-10) dist = 1e-10;
                    
                    pe += _coulombConstant * _particles[i].charge * _particles[j].charge / dist;
                }
            }
            return pe;
        }

        /// @brief Run simulation until equilibrium (kinetic energy drops below threshold)
        /// @param maxSteps Maximum number of simulation steps
        /// @param keThreshold Kinetic energy threshold for equilibrium detection
        /// @param saveEveryN Save position every N steps (for animation)
        /// @param verbose Print progress information
        /// @return History of positions for each particle
        std::vector<std::vector<Point3Cartesian>> runToEquilibrium(
            int maxSteps, 
            Real keThreshold = 1e-8,
            int saveEveryN = 1,
            bool verbose = true)
        {
            std::vector<std::vector<Point3Cartesian>> history(_particles.size());
            
            // Initialize history with initial positions
            for (size_t i = 0; i < _particles.size(); ++i)
            {
                history[i].push_back(_particles[i].pos);
            }

            Real initialPE = totalPotentialEnergy();
            if (verbose)
            {
                std::cout << "Geometry: " << _geometry->ToString() << std::endl;
                std::cout << "Initial potential energy: " << initialPE << std::endl;
                std::cout << "Running simulation..." << std::endl;
            }

            int equilibriumCount = 0;
            const int equilibriumSteps = 100;  // Need KE below threshold for this many steps

            for (int s = 0; s < maxSteps; ++s)
            {
                step();

                if ((s + 1) % saveEveryN == 0)
                {
                    for (size_t i = 0; i < _particles.size(); ++i)
                    {
                        history[i].push_back(_particles[i].pos);
                    }
                }

                Real ke = totalKineticEnergy();
                
                if (ke < keThreshold)
                    equilibriumCount++;
                else
                    equilibriumCount = 0;

                if (equilibriumCount >= equilibriumSteps)
                {
                    if (verbose)
                    {
                        std::cout << "Equilibrium reached at step " << s + 1 << std::endl;
                        std::cout << "Final kinetic energy: " << ke << std::endl;
                        std::cout << "Final potential energy: " << totalPotentialEnergy() << std::endl;
                    }
                    break;
                }

                if (verbose && (s + 1) % 100 == 0)
                {
                    std::cout << "Step " << std::setw(4) << s + 1 << ": KE = " << std::setw(12) << ke 
                              << ", PE = " << std::setw(12) << totalPotentialEnergy() << std::endl;
                }
            }

            return history;
        }

        /// @brief Get number of particles
        int numParticles() const { return static_cast<int>(_particles.size()); }
        
        /// @brief Get reference to the confining geometry
        const IBody& geometry() const { return *_geometry; }
        
        /// @brief Get bounding box of the geometry
        Box3D getBoundingBox() const { return _geometry->GetBoundingBox(); }
        
        /// @brief Get particle positions
        const std::vector<Particle>& particles() const { return _particles; }
    };

} // namespace MPL

#endif // MPL_GENERIC_CHARGE_SIMULATOR_H
