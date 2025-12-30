/******************************************************************************
 * NBodyGravity.h - Self-contained N-Body Gravity Simulation for MML Example
 * ============================================================================
 * 
 * This header provides everything needed for N-body gravitational simulation
 * without requiring the separate MPL (MinimalPhysicsLibrary) dependency.
 * 
 * Extracted and simplified from MPL for use as a standalone MML example.
 * 
 * Contents:
 *   - Physical constants (gravity constant, solar system data)
 *   - Body state classes (position, velocity, mass)
 *   - N-body system configuration and state
 *   - Multiple ODE integrators for different accuracy/performance tradeoffs
 *   - Pre-configured scenarios (Solar System, Many Bodies, Star Cluster Collision)
 * 
 * Available Solvers (NBodyGravitySimulator):
 *   - SolveEuler(dT, numSteps)      - 1st order, fast, terrible energy conservation
 *   - SolveRK4(dT, numSteps)        - 4th order, classic, energy drifts over time
 *   - SolveVerlet(dT, numSteps)     - Symplectic, excellent energy conservation!
 *   - SolveLeapfrog(dT, numSteps)   - Symplectic, equivalent to Verlet
 *   - SolveRK5(duration, eps, ...)  - Adaptive 5th order Cash-Karp
 *   - SolveDP5(duration, eps, ...)  - Adaptive Dormand-Prince 5th order
 *   - SolveDP8(duration, eps, ...)  - Adaptive Dormand-Prince 8th order (highest accuracy)
 * 
 * Energy Conservation Comparison:
 *   For long-term orbital simulations (10+ years), use Verlet/Leapfrog.
 *   For high-accuracy short runs, use DP5 or DP8.
 *   Euler and RK4 will show significant energy drift in orbital simulations.
 * 
 *****************************************************************************/

#ifndef NBODY_GRAVITY_H
#define NBODY_GRAVITY_H

#include "MMLBase.h"

#include "base/BaseUtils.h"
#include "base/Random.h"
#include "base/Function.h"
#include "base/VectorTypes.h"
#include "base/Geometry3D.h"

#include "interfaces/IODESystem.h"

#include "algorithms/ODEAdaptiveIntegrator.h"
#include "algorithms/ODESystemSolver.h"
#include "algorithms/ODESystemStepCalculators.h"
#include "algorithms/ODESystemSteppers.h"

#include "tools/Serializer.h"
#include "tools/Visualizer.h"

namespace NBody
{
    using namespace MML;

    /**************************************************************************
     * PHYSICAL CONSTANTS
     **************************************************************************/
    
    // Gravitational constant in SI units: m³ kg⁻¹ s⁻²
    constexpr double G_SI = 6.67430e-11;
    
    /**************************************************************************
     * SOLAR SYSTEM DATA
     * 
     * Real astronomical data for simulation scenarios
     **************************************************************************/
    
    struct PlanetData {
        const char* name;
        double mass_kg;           // Mass in kg
        double radius_km;         // Radius in km
        double orbital_dist_km;   // Semi-major axis in km (distance from Sun)
        double mass_jupiter;      // Mass relative to Jupiter
    };
    
    // Solar system bodies with real data
    namespace SolarSystem {
        constexpr double SunMass_kg = 1.989e30;
        constexpr double JupiterMass_kg = 1.898e27;
        
        // Planet data: mass (kg), radius (km), orbital distance (million km)
        const PlanetData Planets[] = {
            {"Mercury", 3.285e23,  2439.7,   57.91e6,  1.6601e-7},
            {"Venus",   4.867e24,  6051.8,  108.21e6,  2.447e-6},
            {"Earth",   5.972e24,  6371.0,  149.60e6,  3.003e-6},
            {"Mars",    6.417e23,  3389.5,  227.92e6,  3.227e-7},
            {"Jupiter", 1.898e27, 69911.0,  778.57e6,  1.0},
            {"Saturn",  5.683e26, 58232.0, 1433.53e6,  0.299},
            {"Uranus",  8.681e25, 25362.0, 2872.46e6,  0.0457},
            {"Neptune", 1.024e26, 24622.0, 4495.06e6,  0.0539}
        };
        constexpr int NumPlanets = 8;
    }
    
    /**************************************************************************
     * GRAVITY BODY STATE
     * 
     * Represents the state of a single gravitating body
     **************************************************************************/
    
    class GravityBodyState
    {
        Real _mass;
        Real _radius;
        std::string _color;
        Vec3Cart _position;
        Vec3Cart _velocity;
        
    public:
        GravityBodyState() 
            : _mass(1), _position{0,0,0}, _velocity{0,0,0}, _color("Black"), _radius(10) {}
        
        GravityBodyState(Real mass, Vec3Cart position, Vec3Cart velocity)
            : _mass(mass), _position(position), _velocity(velocity), _color("Black"), _radius(10) {}
        
        GravityBodyState(Real mass, Vec3Cart position, Vec3Cart velocity, 
                        const std::string& color, Real radius)
            : _mass(mass), _position(position), _velocity(velocity), _color(color), _radius(radius) {}
        
        // Accessors
        Real Mass() const { return _mass; }
        Real& Mass() { return _mass; }
        Real Radius() const { return _radius; }
        Real& Radius() { return _radius; }
        std::string Color() const { return _color; }
        std::string& Color() { return _color; }
        
        Vec3Cart R() const { return _position; }
        Vec3Cart& R() { return _position; }
        Vec3Cart V() const { return _velocity; }
        Vec3Cart& V() { return _velocity; }
    };
    
    /**************************************************************************
     * N-BODY STATE
     * 
     * Complete state of an N-body system at a single instant
     **************************************************************************/
    
    class NBodyState
    {
    public:
        std::vector<GravityBodyState> _bodies;
        
        void addBody(const GravityBodyState& body) {
            _bodies.push_back(body);
        }
        void addBody(Real mass, Vec3Cart position, Vec3Cart velocity) {
            _bodies.push_back(GravityBodyState(mass, position, velocity));
        }
        
        int NumBodies() const { return static_cast<int>(_bodies.size()); }
        
        Real Mass(int i) const { return _bodies[i].Mass(); }
        Vec3Cart Pos(int i) const { return _bodies[i].R(); }
        Vec3Cart Vel(int i) const { return _bodies[i].V(); }
        
        void SetPosition(int i, Vec3Cart pos) { _bodies[i].R() = pos; }
        void SetVelocity(int i, Vec3Cart vel) { _bodies[i].V() = vel; }
        
        // Center of mass calculations
        Vec3Cart CentreOfMassPos() const {
            Vec3Cart com(0,0,0);
            Real totalMass = 0.0;
            for (const auto& body : _bodies) {
                com += body.Mass() * body.R();
                totalMass += body.Mass();
            }
            return com / totalMass;
        }
        
        Vec3Cart CentreOfMassVel() const {
            Vec3Cart com(0,0,0);
            Real totalMass = 0.0;
            for (const auto& body : _bodies) {
                com += body.Mass() * body.V();
                totalMass += body.Mass();
            }
            return com / totalMass;
        }
        
        // Conservation quantities
        Vec3Cart LinearMomentum() const {
            Vec3Cart lm(0,0,0);
            for (const auto& body : _bodies)
                lm += body.Mass() * body.V();
            return lm;
        }
        
        Vec3Cart AngularMomentumCM() const {
            Vec3Cart am(0,0,0);
            Vec3Cart com = CentreOfMassPos();
            for (const auto& body : _bodies) {
                Vec3Cart r = body.R() - com;
                am += VectorProduct(r, body.Mass() * body.V());
            }
            return am;
        }
        
        Real TotalKineticEnergy() const {
            Real ke = 0.0;
            for (const auto& body : _bodies)
                ke += 0.5 * body.Mass() * POW2(body.V().NormL2());
            return ke;
        }
        
        // NOTE: Returns potential energy WITHOUT G factor - multiply by G for actual energy!
        Real TotalPotentialEnergyNoG() const {
            Real pe = 0.0;
            for (size_t i = 0; i < _bodies.size(); ++i) {
                for (size_t j = i + 1; j < _bodies.size(); ++j) {
                    Vec3Cart r_ij = _bodies[j].R() - _bodies[i].R();
                    pe -= (_bodies[i].Mass() * _bodies[j].Mass()) / r_ij.NormL2();
                }
            }
            return pe;
        }
        
        Real TotalPotentialEnergy(Real G) const {
            return G * TotalPotentialEnergyNoG();
        }
        
        Real TotalEnergy(Real G) const {
            return TotalKineticEnergy() + TotalPotentialEnergy(G);
        }
    };
    
    /**************************************************************************
     * N-BODY SIMULATION CONFIGURATION
     * 
     * Initial setup for an N-body simulation
     **************************************************************************/
    
    class NBodyGravitySimConfig
    {
        Real _G;
        NBodyState _initState;
        
    public:
        NBodyGravitySimConfig() : _G(G_SI) {}
        NBodyGravitySimConfig(Real G) : _G(G) {}
        
        Real G() const { return _G; }
        int NumBodies() const { return static_cast<int>(_initState._bodies.size()); }
        
        const NBodyState& InitState() const { return _initState; }
        
        void AddBody(Real mass, Vec3Cart position, Vec3Cart velocity) {
            _initState.addBody(mass, position, velocity);
        }
        void AddBody(Real mass, Vec3Cart position, Vec3Cart velocity, 
                    std::string color, Real radius) {
            _initState._bodies.push_back(GravityBodyState(mass, position, velocity, color, radius));
        }
        
        Real Mass(int i) const { return _initState.Mass(i); }
        Real Radius(int i) const { return _initState._bodies[i].Radius(); }
        std::string Color(int i) const { return _initState._bodies[i].Color(); }
        
        Vec3Cart Position(int i) const { return _initState.Pos(i); }
        Vec3Cart Velocity(int i) const { return _initState.Vel(i); }
        
        // Convert to ODE initial conditions vector
        Vector<Real> getInitCond() {
            Vector<Real> initCond(6 * NumBodies());
            for (int i = 0; i < NumBodies(); i++) {
                initCond[3 * i]     = Position(i).X();
                initCond[3 * i + 1] = Position(i).Y();
                initCond[3 * i + 2] = Position(i).Z();
                initCond[3 * NumBodies() + 3 * i]     = Velocity(i).X();
                initCond[3 * NumBodies() + 3 * i + 1] = Velocity(i).Y();
                initCond[3 * NumBodies() + 3 * i + 2] = Velocity(i).Z();
            }
            return initCond;
        }
    };
    
    /**************************************************************************
     * PRE-CONFIGURED SIMULATION SCENARIOS
     **************************************************************************/
    
    class NBodyGravityConfigGenerator
    {
    public:
        // Scenario 1: Real Solar System with planets
        static NBodyGravitySimConfig Config1_Solar_system()
        {
            // Unit conversions for convenient simulation units
            constexpr double million_km = 1.0e9;  // meters
            constexpr double year = 3.15576e7;    // seconds
            
            double M_jup = SolarSystem::JupiterMass_kg;
            double M_sun = SolarSystem::SunMass_kg;
            
            // G in units: [(million km)³ / (Jupiter mass × year²)]
            double G = G_SI * M_jup * (year * year) / (million_km * million_km * million_km);
            
            NBodyGravitySimConfig config(G);
            
            // Add Sun at origin
            config.AddBody(M_sun / M_jup, Vec3Cart{0,0,0}, Vec3Cart{0,0,0}, "Yellow", 30);
            
            // Planet visualization colors and sizes
            const char* colors[] = {"Gray", "Orange", "Blue", "Red", "Brown", "Gold", "LightBlue", "Blue"};
            double radii[] = {2.5, 6, 6, 4, 50, 40, 20, 20};
            
            // Add planets with circular orbit velocities
            for (int i = 0; i < SolarSystem::NumPlanets; i++) {
                const auto& p = SolarSystem::Planets[i];
                double dist_million_km = p.orbital_dist_km / 1e6;
                
                // Circular orbit velocity: v = sqrt(G * M_sun / r)
                double v = std::sqrt(G * (M_sun / M_jup) / dist_million_km);
                
                // Place planet at (r, 0, 0), velocity (0, v, small_z)
                config.AddBody(
                    p.mass_jupiter,
                    Vec3Cart{Real(dist_million_km), Real(0), Real(0)},
                    Vec3Cart{Real(0), Real(v), Real(10)},
                    colors[i],
                    radii[i]
                );
            }
            
            return config;
        }
        
        // Scenario 2: 100 random bodies orbiting a massive central object
        static NBodyGravitySimConfig Config2_N_bodies_around_massive_object()
        {
            Real G = 1.0;
            NBodyGravitySimConfig config(G);
            
            // Central massive body
            Real centralMass = 10000;
            config.AddBody(centralMass, Vec3Cart{0,0,0}, Vec3Cart{0,0,0}, "Yellow", 15);
            
            // Create 100 bodies in BOUND orbits
            for (int i = 0; i < 100; i++) {
                Real mass = 1 + Random::UniformReal(0.0, 5.0);  // Small masses
                Real radius = 2;  // Visual size
                
                // Random spherical position - varied distances
                Real rad = Random::UniformReal(50.0, 200.0);
                Real theta = Random::UniformReal(0.0, 2 * Constants::PI);
                Real phi = Random::UniformReal(0.3, Constants::PI - 0.3);  // Avoid poles
                Vec3Cart pos = Vec3Cart{
                    rad * sin(phi) * cos(theta),
                    rad * sin(phi) * sin(theta),
                    rad * cos(phi)
                };
                
                // Circular orbital velocity: v = sqrt(G*M/r)
                Real v_circular = std::sqrt(G * centralMass / rad);
                
                // Add some randomness (0.7 to 1.0 of circular velocity for elliptical bound orbits)
                Real v_factor = Random::UniformReal(0.7, 1.0);
                Real v_mag = v_circular * v_factor;
                
                // Velocity perpendicular to radius (tangential)
                // Pick a random tangent direction in the plane perpendicular to pos
                Vec3Cart radial = pos.GetAsUnitVector();
                
                // Create orthogonal basis
                Vec3Cart arbitrary = (std::abs(radial.X()) < 0.9) ? Vec3Cart{1,0,0} : Vec3Cart{0,1,0};
                Vec3Cart tangent1 = VectorProduct(radial, arbitrary).GetAsUnitVector();
                Vec3Cart tangent2 = VectorProduct(radial, tangent1).GetAsUnitVector();
                
                // Random angle in tangent plane
                Real angle = Random::UniformReal(0.0, 2 * Constants::PI);
                Vec3Cart vel = (tangent1 * cos(angle) + tangent2 * sin(angle)) * v_mag;
                
                config.AddBody(mass, pos, vel, "Cyan", radius);
            }
            
            return config;
        }
        
        // Scenario 3: COLLISION OF STAR CLUSTERS
        // Two 100-body clusters approach and swing past each other - spectacular gravitational dynamics!
        static NBodyGravitySimConfig Config3_StarClusterCollision(
            int bodiesPerCluster = 100,
            Real separation = 800.0,
            Real impactParameter = 100.0,  // Vertical offset - clusters pass by, not head-on!
            Real approachSpeed = 2.5,
            Real G = 1.0)
        {
            NBodyGravitySimConfig config(G);
            
            // Helper lambda: Generate a cluster at given offset with bulk velocity
            auto generateCluster = [&](Vec3Cart offset, Vec3Cart bulkVel, 
                                       std::string starColor, Real centralMass, Real centralRadius) 
            {
                // Add central massive body (like a black hole or dense core)
                config.AddBody(centralMass, offset, bulkVel, "White", centralRadius);
                
                // Generate orbiting bodies
                for (int i = 0; i < bodiesPerCluster; i++) {
                    Real mass = 1 + Random::UniformReal(0.0, 3.0);
                    Real radius = 1.5;
                    
                    // Random spherical position around cluster center
                    Real r = Random::UniformReal(25.0, 120.0);
                    Real theta = Random::UniformReal(0.0, 2 * Constants::PI);
                    Real phi = Random::UniformReal(0.4, Constants::PI - 0.4);  // Avoid poles
                    
                    Vec3Cart localPos{
                        r * sin(phi) * cos(theta),
                        r * sin(phi) * sin(theta),
                        r * cos(phi)
                    };
                    
                    // Circular orbit velocity magnitude
                    Real v_circ = std::sqrt(G * centralMass / r);
                    Real v_mag = v_circ * Random::UniformReal(0.7, 0.95);  // Bound orbits
                    
                    // Velocity perpendicular to radial (tangential orbit)
                    Vec3Cart radial = localPos.GetAsUnitVector();
                    Vec3Cart arbitrary = (std::abs(radial.X()) < 0.9) ? Vec3Cart{1,0,0} : Vec3Cart{0,1,0};
                    Vec3Cart tangent1 = VectorProduct(radial, arbitrary).GetAsUnitVector();
                    Vec3Cart tangent2 = VectorProduct(radial, tangent1).GetAsUnitVector();
                    
                    Real angle = Random::UniformReal(0.0, 2 * Constants::PI);
                    Vec3Cart localVel = (tangent1 * cos(angle) + tangent2 * sin(angle)) * v_mag;
                    
                    // Apply cluster offset and bulk motion
                    Vec3Cart pos = offset + localPos;
                    Vec3Cart vel = bulkVel + localVel;
                    
                    config.AddBody(mass, pos, vel, starColor, radius);
                }
            };
            
            // Cluster A (Cyan) - starts on left (offset DOWN in Y), moves right
            Vec3Cart offsetA{-separation / 2, -impactParameter / 2, 0};
            Vec3Cart velA{approachSpeed / 2, 0, 0};
            generateCluster(offsetA, velA, "Cyan", 5000.0, 8.0);
            
            // Cluster B (Magenta) - starts on right (offset UP in Y), moves left  
            Vec3Cart offsetB{+separation / 2, +impactParameter / 2, 0};
            Vec3Cart velB{-approachSpeed / 2, 0, 0};
            generateCluster(offsetB, velB, "Magenta", 5000.0, 8.0);
            
            return config;
        }
    };
    
    /**************************************************************************
     * SIMULATION RESULTS
     * 
     * Stores the complete trajectory of an N-body simulation
     **************************************************************************/
    
    class NBodyGravitySimulationResults
    {
    public:
        const NBodyGravitySimConfig& _config;
        Real _duration;
        Vector<Real> _vecTimes;
        Vector<NBodyState> _vecStates;
        
        NBodyGravitySimulationResults(const NBodyGravitySimConfig& config)
            : _config(config), _duration(0) {}
        
        NBodyGravitySimulationResults(const NBodyGravitySimConfig& config, Real duration)
            : _config(config), _duration(duration) {}
        
        int NumBodies() const { return _config.NumBodies(); }
        int NumSteps() const { return _vecStates.size(); }
        
        Real Time(int i) const { return _vecTimes[i]; }
        NBodyState State(int i) const { return _vecStates[i]; }
        
        Vector<Real> getTimes() const { return _vecTimes; }
        
        // Position getters for visualization
        Vector<Real> getPosX(int i) const {
            Vector<Real> res(_vecStates.size());
            for (int j = 0; j < _vecStates.size(); j++)
                res[j] = _vecStates[j].Pos(i).X();
            return res;
        }
        Vector<Real> getPosY(int i) const {
            Vector<Real> res(_vecStates.size());
            for (int j = 0; j < _vecStates.size(); j++)
                res[j] = _vecStates[j].Pos(i).Y();
            return res;
        }
        Vector<Real> getPosZ(int i) const {
            Vector<Real> res(_vecStates.size());
            for (int j = 0; j < _vecStates.size(); j++)
                res[j] = _vecStates[j].Pos(i).Z();
            return res;
        }
        
        // Energy getters for conservation checks
        Vector<Real> getTotalKineticEnergy() const {
            Vector<Real> res(_vecStates.size());
            for (int j = 0; j < _vecStates.size(); j++)
                res[j] = _vecStates[j].TotalKineticEnergy();
            return res;
        }
        Vector<Real> getTotalPotentialEnergy() const {
            Vector<Real> res(_vecStates.size());
            Real G = _config.G();
            for (int j = 0; j < _vecStates.size(); j++)
                res[j] = _vecStates[j].TotalPotentialEnergy(G);
            return res;
        }
        Vector<Real> getTotalEnergy() const {
            Vector<Real> res(_vecStates.size());
            Real G = _config.G();
            for (int j = 0; j < _vecStates.size(); j++)
                res[j] = _vecStates[j].TotalEnergy(G);
            return res;
        }
        
        // Energy conservation metrics
        Real getInitialEnergy() const {
            return _vecStates[0].TotalEnergy(_config.G());
        }
        Real getFinalEnergy() const {
            return _vecStates[_vecStates.size() - 1].TotalEnergy(_config.G());
        }
        Real getEnergyDrift() const {
            return getFinalEnergy() - getInitialEnergy();
        }
        Real getRelativeEnergyError() const {
            Real E0 = getInitialEnergy();
            if (std::abs(E0) < 1e-15) return 0.0;
            return std::abs(getEnergyDrift() / E0);
        }
        Real getMaxEnergyDeviation() const {
            Real E0 = getInitialEnergy();
            Real maxDev = 0.0;
            Real G = _config.G();
            for (int j = 0; j < _vecStates.size(); j++) {
                Real dev = std::abs(_vecStates[j].TotalEnergy(G) - E0);
                if (dev > maxDev) maxDev = dev;
            }
            return maxDev;
        }
        Real getMaxRelativeEnergyError() const {
            Real E0 = getInitialEnergy();
            if (std::abs(E0) < 1e-15) return 0.0;
            return getMaxEnergyDeviation() / std::abs(E0);
        }
        
        // Visualization: 3D parametric curves for trajectories
        void VisualizeAsParamCurve(std::string baseFileName, Vector<int> vecBodiesIndexToVisualize)
        {
            Vector<Real> t_vals = getTimes();
            std::vector<std::string> fileNames;
            
            for (int i : vecBodiesIndexToVisualize) {
                Vector<Real> body_x = getPosX(i);
                Vector<Real> body_y = getPosY(i);
                Vector<Real> body_z = getPosZ(i);
                
                // Build parametric curve data
                std::vector<VectorN<Real, 3>> res;
                for (int j = 0; j < t_vals.size(); j++)
                    res.push_back(VectorN<Real, 3>{body_x[j], body_y[j], body_z[j]});
                
                std::string fileName = baseFileName + std::to_string(i) + ".txt";
                fileNames.push_back(fileName);
                
                Real t1 = t_vals[0];
                Real t2 = t_vals[t_vals.size() - 1];
                
                Serializer::SaveAsParamCurve<3>(res, "PARAMETRIC_CURVE_CARTESIAN_3D",
                    baseFileName + std::to_string(i), t1, t2, t_vals.size(),
                    GetResultFilesPath() + fileName);
            }
            
            Visualizer::VisualizeMultiParamCurve3D(fileNames);
        }
        
        // Visualization: Particle animation
        void VisualizeAsParticleSimulation(std::string baseFileName, 
                                          Vector<int> vecBodiesIndexToVisualize, double dT)
        {
            std::vector<std::string> vecColors(NumBodies());
            std::vector<Real> vecRad(NumBodies());
            for (int i = 0; i < NumBodies(); i++) {
                vecColors[i] = _config.Color(i);
                vecRad[i] = _config.Radius(i);
            }
            
            std::vector<std::vector<Pnt3Cart>> res2;
            res2.resize(NumBodies());
            
            // Compute bounding box from trajectory data
            Real minX = 1e30, maxX = -1e30;
            Real minY = 1e30, maxY = -1e30;
            Real minZ = 1e30, maxZ = -1e30;
            
            for (int i = 0; i < NumBodies(); i++) {
                res2[i].resize(NumSteps());
                for (int j = 0; j < NumSteps(); j++) {
                    Real x = _vecStates[j].Pos(i).X();
                    Real y = _vecStates[j].Pos(i).Y();
                    Real z = _vecStates[j].Pos(i).Z();
                    res2[i][j] = Pnt3Cart(x, y, z);
                    
                    minX = std::min(minX, x); maxX = std::max(maxX, x);
                    minY = std::min(minY, y); maxY = std::max(maxY, y);
                    minZ = std::min(minZ, z); maxZ = std::max(maxZ, z);
                }
            }
            
            // Add 10% margin and compute dimensions
            Real margin = 0.1;
            Real width  = (maxX - minX) * (1 + margin);
            Real height = (maxY - minY) * (1 + margin);
            Real depth  = (maxZ - minZ) * (1 + margin);
            
            // Ensure minimum dimensions (avoid zero for 2D orbits)
            width  = std::max(width, 1.0);
            height = std::max(height, 1.0);
            depth  = std::max(depth, 1.0);
            
            std::string fullFileName = baseFileName + ".txt";
            Serializer::SaveParticleSimulation3D(GetResultFilesPath() + fullFileName,
                NumBodies(), width, height, depth, res2, vecColors, vecRad, dT);
            
            Visualizer::VisualizeParticleSimulation3D(fullFileName);
        }
    };
    
    /**************************************************************************
     * N-BODY ODE SYSTEM
     * 
     * Formulates N-body gravity as an ODE system for use with RK5 solver
     **************************************************************************/
    
    class NBodyGravitySystemODE : public IODESystem
    {
        NBodyGravitySimConfig _config;
        
    public:
        NBodyGravitySystemODE(NBodyGravitySimConfig config) : _config(config) {}
        
        int getDim() const override { return 6 * _config.NumBodies(); }
        
        void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const override
        {
            int N = _config.NumBodies();
            
            for (int i = 0; i < N; i++) {
                Vec3Cart force(0, 0, 0);
                
                // Calculate gravitational force on body i from all other bodies
                for (int j = 0; j < N; j++) {
                    if (i != j) {
                        Vec3Cart vec_dist(x[3*j] - x[3*i],
                                         x[3*j + 1] - x[3*i + 1],
                                         x[3*j + 2] - x[3*i + 2]);
                        
                        force = force + _config.G() * _config.Mass(i) * _config.Mass(j) 
                                      / POW3(vec_dist.NormL2()) * vec_dist;
                    }
                }
                
                // Position derivatives = velocities
                int half = N * 3;
                dxdt[3*i]     = x[half + 3*i];
                dxdt[3*i + 1] = x[half + 3*i + 1];
                dxdt[3*i + 2] = x[half + 3*i + 2];
                
                // Velocity derivatives = accelerations (F/m)
                dxdt[half + 3*i]     = force.X() / _config.Mass(i);
                dxdt[half + 3*i + 1] = force.Y() / _config.Mass(i);
                dxdt[half + 3*i + 2] = force.Z() / _config.Mass(i);
            }
        }
    };
    
    /**************************************************************************
     * N-BODY GRAVITY SIMULATOR
     * 
     * Main simulator class with multiple integration methods.
     * 
     * Fixed-step methods (specify timestep and number of steps):
     *   - SolveEuler()    - Simple Euler (illustrative, NOT recommended)
     *   - SolveRK4()      - Classic Runge-Kutta 4th order
     *   - SolveVerlet()   - Velocity Verlet (SYMPLECTIC - best for orbits!)
     *   - SolveLeapfrog() - Leapfrog (SYMPLECTIC - equivalent to Verlet)
     * 
     * Adaptive methods (specify duration, tolerance, save interval):
     *   - SolveRK5()  - Cash-Karp 5th order adaptive
     *   - SolveDP5()  - Dormand-Prince 5th order adaptive
     *   - SolveDP8()  - Dormand-Prince 8th order adaptive (highest accuracy)
     * 
     * For long-term orbital simulations, symplectic methods (Verlet/Leapfrog)
     * are strongly recommended as they conserve energy much better than
     * non-symplectic methods like Euler or RK4.
     **************************************************************************/
    
    class NBodyGravitySimulator
    {
    protected:
        NBodyGravitySimConfig _config;
        
    public:
        NBodyGravitySimulator(const NBodyGravitySimConfig& config) : _config(config) {}
        
        // Simple Euler integration (first-order, illustrative)
        NBodyGravitySimulationResults SolveEuler(Real dT, int numSteps)
        {
            NBodyGravitySimulationResults results(_config, dT * numSteps);
            
            results._vecTimes.Resize(numSteps);
            results._vecStates.Resize(numSteps);
            
            // Get initial state from configuration
            NBodyState sysState = _config.InitState();
            
            // Save initial state
            results._vecTimes[0] = 0.0;
            for (int i = 0; i < sysState.NumBodies(); i++)
                results._vecStates[0].addBody(sysState.Mass(i), sysState.Pos(i), sysState.Vel(i));
            
            // Perform Euler integration
            for (int step = 1; step < numSteps; step++) {
                results._vecTimes[step] = results._vecTimes[step - 1] + dT;
                
                // Calculate forces
                Vector<Vec3Cart> force(sysState.NumBodies(), Vec3Cart(0, 0, 0));
                for (int j = 0; j < sysState.NumBodies() - 1; j++) {
                    for (int k = j + 1; k < sysState.NumBodies(); k++) {
                        Vec3Cart radialVec = sysState.Pos(j) - sysState.Pos(k);
                        Vec3Cart F = _config.G() * sysState.Mass(j) * sysState.Mass(k) 
                                   / POW3(radialVec.NormL2()) * radialVec;
                        force[j] = force[j] - F;
                        force[k] = force[k] + F;  // Newton's 3rd law
                    }
                }
                
                // Update velocities and positions
                for (int j = 0; j < _config.NumBodies(); j++) {
                    sysState.SetVelocity(j, sysState.Vel(j) + force[j] * dT / sysState.Mass(j));
                    sysState.SetPosition(j, sysState.Pos(j) + sysState.Vel(j) * dT);
                }
                
                // Save state
                for (int j = 0; j < _config.NumBodies(); j++)
                    results._vecStates[step].addBody(sysState.Mass(j), sysState.Pos(j), sysState.Vel(j));
            }
            
            return results;
        }
        
        // Adaptive RK5 Cash-Karp integration (production quality)
        NBodyGravitySimulationResults SolveRK5(Real duration, Real eps, Real minSaveInterval, Real hStart)
        {
            NBodyGravitySystemODE ode(_config);
            
            ODESystemSolver<RK5_CashKarp_Stepper> rk5solver(ode);
            ODESystemSolution sol = rk5solver.integrate(_config.getInitCond(), 0, duration, 
                                                       minSaveInterval, eps, hStart);
            
            return convertODESolutionToResults(sol, duration);
        }

        // Dormand-Prince 5th order adaptive integrator (very popular, good accuracy)
        NBodyGravitySimulationResults SolveDP5(Real duration, Real eps, Real minSaveInterval, Real hStart)
        {
            NBodyGravitySystemODE ode(_config);
            
            DormandPrince5Integrator dp5solver(ode);
            ODESystemSolution sol = dp5solver.integrate(_config.getInitCond(), 0, duration, 
                                                        minSaveInterval, eps, hStart);
            
            return convertODESolutionToResults(sol, duration);
        }

        // Dormand-Prince 8th order adaptive integrator (highest accuracy, expensive)
        NBodyGravitySimulationResults SolveDP8(Real duration, Real eps, Real minSaveInterval, Real hStart)
        {
            NBodyGravitySystemODE ode(_config);
            
            DormandPrince8Integrator dp8solver(ode);
            ODESystemSolution sol = dp8solver.integrate(_config.getInitCond(), 0, duration, 
                                                        minSaveInterval, eps, hStart);
            
            return convertODESolutionToResults(sol, duration);
        }

        // Velocity Verlet integration (symplectic - excellent energy conservation!)
        // Perfect for long-term orbital simulations where energy conservation matters
        NBodyGravitySimulationResults SolveVerlet(Real dT, int numSteps)
        {
            NBodyGravitySystemODE ode(_config);
            
            VelocityVerlet_StepCalculator verletCalc;
            ODESystemFixedStepSolver solver(ode, verletCalc);
            ODESystemSolution sol = solver.integrate(_config.getInitCond(), 0, dT * numSteps, numSteps);
            
            return convertODESolutionToResults(sol, dT * numSteps);
        }

        // Leapfrog integration (symplectic - equivalent to Velocity Verlet)
        // Another symplectic integrator with bounded energy error
        NBodyGravitySimulationResults SolveLeapfrog(Real dT, int numSteps)
        {
            NBodyGravitySystemODE ode(_config);
            
            Leapfrog_StepCalculator leapfrogCalc;
            ODESystemFixedStepSolver solver(ode, leapfrogCalc);
            ODESystemSolution sol = solver.integrate(_config.getInitCond(), 0, dT * numSteps, numSteps);
            
            return convertODESolutionToResults(sol, dT * numSteps);
        }

        // RK4 fixed-step integration (classic, but energy drifts over time)
        NBodyGravitySimulationResults SolveRK4(Real dT, int numSteps)
        {
            NBodyGravitySystemODE ode(_config);
            
            RungeKutta4_StepCalculator rk4Calc;
            ODESystemFixedStepSolver solver(ode, rk4Calc);
            ODESystemSolution sol = solver.integrate(_config.getInitCond(), 0, dT * numSteps, numSteps);
            
            return convertODESolutionToResults(sol, dT * numSteps);
        }

    private:
        // Helper to convert ODESystemSolution to NBodyGravitySimulationResults
        NBodyGravitySimulationResults convertODESolutionToResults(const ODESystemSolution& sol, Real duration)
        {
            NBodyGravitySimulationResults results(_config, duration);
            
            results._duration = duration;
            results._vecTimes = sol.getTValues();
            
            int numBodies = _config.NumBodies();
            int numSteps = sol.getTValues().size();
            results._vecStates.Resize(numSteps);
            
            for (int i = 0; i < numSteps; i++) {
                NBodyState state;
                for (int j = 0; j < numBodies; j++) {
                    Vec3Cart pos(sol.getXValues()[3*j][i], 
                                sol.getXValues()[3*j + 1][i], 
                                sol.getXValues()[3*j + 2][i]);
                    Vec3Cart vel(sol.getXValues()[numBodies*3 + 3*j][i],
                                sol.getXValues()[numBodies*3 + 3*j + 1][i],
                                sol.getXValues()[numBodies*3 + 3*j + 2][i]);
                    state.addBody(_config.Mass(j), pos, vel);
                }
                results._vecStates[i] = state;
            }
            
            return results;
        }
    };

} // namespace NBody

#endif // NBODY_GRAVITY_H
