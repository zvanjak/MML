#if !defined MML_ODE_SYSTEM_TEST_BED_H
#define MML_ODE_SYSTEM_TEST_BED_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/Matrix.h"
#include "base/ODESystem.h"
#endif

#include "diff_eq_systems_defs.h"

#include <vector>
#include <functional>

namespace MML::TestBeds
{
    /*******************************************************************************************************************
     * ODE TEST BED - Comprehensive Testing Framework for ODE Solvers
     * 
     * Provides:
     *   1. Systems with full analytical solutions (verify at any time)
     *   2. Systems with end-point solutions (perfect for unit testing)
     *   3. Category-based organization (linear, oscillatory, nonlinear, chaotic)
     *   4. Filtering by dimension, stiffness, etc.
     * 
     * NOTE: Stiff systems are primarily in stiff_systems_test_bed.h
     *******************************************************************************************************************/

    //==================================================================================================================
    // INTERFACES FOR TEST SYSTEMS
    //==================================================================================================================

    /**
     * @brief Interface for ODE systems with known analytical solutions
     * 
     * These allow verification at any time point within the integration domain.
     */
    class ITestODESystemWithSolution
    {
    public:
        virtual ~ITestODESystemWithSolution() = default;
        
        virtual const IODESystem* getODESystem() const = 0;
        
        virtual Vector<Real> getInitCondLow() const = 0;
        virtual Vector<Real> getInitCondHigh() const = 0;
        
        virtual Vector<Real> getSolution(const Vector<Real>& initCond, Real t) const = 0;
    };

    /**
     * @brief Interface for ODE systems with known end-point solutions
     * 
     * Perfect for unit testing: integrate from t0 to t_end and compare.
     */
    class ITestODESystemWithEndSolution
    {
    public:
        virtual ~ITestODESystemWithEndSolution() = default;
        
        virtual const IODESystem* getODESystem() const = 0;
        virtual Vector<Real> getInitialConditions() const = 0;
        virtual Real         getStartTime() const = 0;
        virtual Real         getEndTime() const = 0;
        virtual Vector<Real> getEndSolution() const = 0;
    };

    //==================================================================================================================
    // CONCRETE WRAPPER CLASSES
    //==================================================================================================================

    /**
     * @brief Wrapper for systems with analytical solutions
     */
    template<typename ODESystemType>
    class TestODESystemWithSolutionWrapper : public ITestODESystemWithSolution
    {
        ODESystemType _odeSys;
        Vector<Real> _initCondLow;
        Vector<Real> _initCondHigh;
        
    public:
        template<typename... Args>
        TestODESystemWithSolutionWrapper(Vector<Real> icLow, Vector<Real> icHigh, Args&&... args)
            : _odeSys(std::forward<Args>(args)...), _initCondLow(icLow), _initCondHigh(icHigh) {}
        
        const IODESystem* getODESystem() const override { return &_odeSys; }
        Vector<Real> getInitCondLow() const override { return _initCondLow; }
        Vector<Real> getInitCondHigh() const override { return _initCondHigh; }
        
        Vector<Real> getSolution(const Vector<Real>& initCond, Real t) const override
        {
            return _odeSys.getSolution(initCond, t);
        }
        
        const ODESystemType& getTypedODESystem() const { return _odeSys; }
    };

    /**
     * @brief Wrapper for systems with end-point solutions only
     */
    class TestODESystemWithEndSolution : public ITestODESystemWithEndSolution
    {
        std::shared_ptr<IODESystem> _ownedOdeSystem;
        const IODESystem* _odeSystem;
        
        Vector<Real> _initial_conditions;
        Real         _start_time;
        Real         _end_time;
        Vector<Real> _end_solution;
        
    public:
        // Constructor with function pointer (creates ODESystem internally)
        TestODESystemWithEndSolution(int n, void (*inFunc)(Real, const Vector<Real>&, Vector<Real>&),
                                     Vector<Real> init_cond, Real start_time, Real end_time, Vector<Real> end_solution)
            : _ownedOdeSystem(std::make_shared<ODESystem>(n, inFunc)),
              _odeSystem(_ownedOdeSystem.get()),
              _initial_conditions(init_cond), _start_time(start_time), 
              _end_time(end_time), _end_solution(end_solution) {}
        
        // Constructor with shared_ptr to IODESystem
        TestODESystemWithEndSolution(std::shared_ptr<IODESystem> odeSystem,
                                     Vector<Real> init_cond, Real start_time, Real end_time, Vector<Real> end_solution)
            : _ownedOdeSystem(odeSystem),
              _odeSystem(odeSystem.get()),
              _initial_conditions(init_cond), _start_time(start_time),
              _end_time(end_time), _end_solution(end_solution) {}

        // Constructor with pointer to existing IODESystem (no ownership)
        TestODESystemWithEndSolution(const IODESystem* odeSystem,
                                     Vector<Real> init_cond, Real start_time, Real end_time, Vector<Real> end_solution)
            : _ownedOdeSystem(nullptr),
              _odeSystem(odeSystem),
              _initial_conditions(init_cond), _start_time(start_time),
              _end_time(end_time), _end_solution(end_solution) {}
        
        // Legacy constructor (t0 = 0)
        TestODESystemWithEndSolution(int n, void (*inFunc)(Real, const Vector<Real>&, Vector<Real>&),
                                     Vector<Real> init_cond, Real end_time, Vector<Real> end_solution)
            : TestODESystemWithEndSolution(n, inFunc, init_cond, 0.0, end_time, end_solution) {}
        
        TestODESystemWithEndSolution(const IODESystem* odeSystem,
                                     Vector<Real> init_cond, Real end_time, Vector<Real> end_solution)
            : TestODESystemWithEndSolution(odeSystem, init_cond, 0.0, end_time, end_solution) {}
        
        const IODESystem* getODESystem() const override { return _odeSystem; }
        Vector<Real> getInitialConditions() const override { return _initial_conditions; }
        Real         getStartTime() const override { return _start_time; }
        Real         getEndTime() const override { return _end_time; }
        Vector<Real> getEndSolution() const override { return _end_solution; }
    };

    //==================================================================================================================
    // SYSTEM CATEGORY ENUMERATION
    //==================================================================================================================

    enum class ODECategory
    {
        Linear,           // Linear systems (exact analytical solutions)
        Oscillatory,      // Harmonic oscillators (simple, damped)
        Nonlinear,        // Nonlinear oscillators, pendulum
        Population,       // Population dynamics (Lotka-Volterra, logistic)
        Chaotic,          // Lorenz, Rossler
        Physical,         // Projectile, Kepler
        SpecialFunction,  // Legendre, Bessel, Hermite, Laguerre
        Stiff             // Stiff systems (single example here)
    };

    //==================================================================================================================
    // MAIN TEST BED CLASS
    //==================================================================================================================

    class ODESystemTestBed
    {
    private:
        //--------------------------------------------------------------------------------------------------------------
        // STATIC INSTANCES OF ODE SYSTEMS
        //--------------------------------------------------------------------------------------------------------------
        
        // Linear systems (with analytical solutions)
        static inline ExponentialDecayODE _expDecay{1.0};
        static inline ExponentialGrowthODE _expGrowth{0.5};
        static inline Linear2DDistinctEigenODE _linear2D;
        static inline Linear3DSystemODE _linear3D;
        
        // Oscillators
        static inline SimpleHarmonicOscillatorODE _sho{1.0};
        static inline SimpleHarmonicOscillatorODE _shoOmega2{2.0};
        static inline DampedHarmonicOscillatorODE _dampedSHO{1.0, 0.1};
        
        // Nonlinear oscillators
        static inline VanDerPolODE _vanDerPol{0.3};
        static inline DuffingODE _duffing;
        static inline SimplePendulumODE _pendulum{9.81};
        
        // Population dynamics
        static inline LotkaVolterraODE _lotkaVolterra;
        static inline LogisticGrowthODE _logistic{1.0, 10.0};
        
        // Chaotic systems
        static inline LorenzSystemODE _lorenz;
        static inline RosslerSystemODE _rossler;
        
        // Physical systems
        static inline ProjectileWithDragODE _projectile{0.05, 9.81};
        static inline KeplerProblemODE _kepler{1.0};
        
        // Special functions
        static inline LegendreODE _legendre2{2};
        static inline LaguerreODE _laguerre3{3};
        static inline HermiteODE _hermite2{2};
        static inline BesselODE _bessel0{0};
        
        // Stiff (example only)
        static inline TestStiffSysWithJacobian1 _stiff1;

        //--------------------------------------------------------------------------------------------------------------
        // SYSTEMS WITH ANALYTICAL SOLUTIONS
        //--------------------------------------------------------------------------------------------------------------
        
        struct ODESystemWithSolutionEntry
        {
            std::string name;
            ODECategory category;
            const IODESystem* system;
            std::function<Vector<Real>(const Vector<Real>&, Real)> solution;
            Vector<Real> icLow;
            Vector<Real> icHigh;
        };
        
        static inline std::vector<ODESystemWithSolutionEntry> _systemsWithSolution = {
            // Linear systems
            { "ExpDecay_k1", ODECategory::Linear, &_expDecay,
              [](const Vector<Real>& y0, Real t) { return _expDecay.getSolution(y0, t); },
              Vector<Real>{0.1}, Vector<Real>{10.0} },
            
            { "ExpGrowth_k0.5", ODECategory::Linear, &_expGrowth,
              [](const Vector<Real>& y0, Real t) { return _expGrowth.getSolution(y0, t); },
              Vector<Real>{0.1}, Vector<Real>{5.0} },
            
            { "Linear2D_DistinctEigen", ODECategory::Linear, &_linear2D,
              [](const Vector<Real>& y0, Real t) { return _linear2D.getSolution(y0, t); },
              Vector<Real>{-5.0, -5.0}, Vector<Real>{5.0, 5.0} },
            
            // Oscillators
            { "SHO_omega1", ODECategory::Oscillatory, &_sho,
              [](const Vector<Real>& y0, Real t) { return _sho.getSolution(y0, t); },
              Vector<Real>{0.0, -10.0}, Vector<Real>{Constants::PI, 10.0} },
            
            { "SHO_omega2", ODECategory::Oscillatory, &_shoOmega2,
              [](const Vector<Real>& y0, Real t) { return _shoOmega2.getSolution(y0, t); },
              Vector<Real>{0.0, -10.0}, Vector<Real>{Constants::PI, 10.0} },
            
            { "DampedSHO_zeta0.1", ODECategory::Oscillatory, &_dampedSHO,
              [](const Vector<Real>& y0, Real t) { return _dampedSHO.getSolution(y0, t); },
              Vector<Real>{0.0, -10.0}, Vector<Real>{Constants::PI, 10.0} },
            
            // Population dynamics
            { "Logistic_r1_K10", ODECategory::Population, &_logistic,
              [](const Vector<Real>& y0, Real t) { return _logistic.getSolution(y0, t); },
              Vector<Real>{0.1}, Vector<Real>{9.0} }
        };

        //--------------------------------------------------------------------------------------------------------------
        // SYSTEMS WITH END-POINT SOLUTIONS (computed from analytical solutions)
        //--------------------------------------------------------------------------------------------------------------
        
        struct ODESystemWithEndSolutionEntry
        {
            std::string name;
            ODECategory category;
            const IODESystem* system;
            Vector<Real> ic;
            Real t0;
            Real tEnd;
            Vector<Real> endSolution;
        };
        
        static inline std::vector<ODESystemWithEndSolutionEntry> _systemsWithEndSolution = {
            // Linear systems (short, medium, long integration)
            { "ExpDecay_t5", ODECategory::Linear, &_expDecay,
              Vector<Real>{1.0}, 0.0, 5.0, _expDecay.getSolution(Vector<Real>{1.0}, 5.0) },
            
            { "ExpGrowth_t2", ODECategory::Linear, &_expGrowth,
              Vector<Real>{1.0}, 0.0, 2.0, _expGrowth.getSolution(Vector<Real>{1.0}, 2.0) },
            
            { "Linear2D_t3", ODECategory::Linear, &_linear2D,
              Vector<Real>{1.0, 3.0}, 0.0, 3.0, _linear2D.getSolution(Vector<Real>{1.0, 3.0}, 3.0) },
            
            { "Linear3D_t1", ODECategory::Linear, &_linear3D,
              Vector<Real>{1.0, 1.0, 1.0}, 0.0, 1.0, _linear3D.getSolutionForUnitIC(1.0) },
            
            // Oscillators (multiple periods)
            { "SHO_1period", ODECategory::Oscillatory, &_sho,
              Vector<Real>{1.0, 0.0}, 0.0, 2*Constants::PI, _sho.getSolution(Vector<Real>{1.0, 0.0}, 2*Constants::PI) },
            
            { "SHO_5periods", ODECategory::Oscillatory, &_sho,
              Vector<Real>{1.0, 0.0}, 0.0, 10*Constants::PI, _sho.getSolution(Vector<Real>{1.0, 0.0}, 10*Constants::PI) },
            
            { "SHO_omega2_2periods", ODECategory::Oscillatory, &_shoOmega2,
              Vector<Real>{1.0, 0.0}, 0.0, 2*Constants::PI, _shoOmega2.getSolution(Vector<Real>{1.0, 0.0}, 2*Constants::PI) },
            
            { "DampedSHO_t10", ODECategory::Oscillatory, &_dampedSHO,
              Vector<Real>{1.0, 0.0}, 0.0, 10.0, _dampedSHO.getSolution(Vector<Real>{1.0, 0.0}, 10.0) },
            
            { "DampedSHO_t20", ODECategory::Oscillatory, &_dampedSHO,
              Vector<Real>{1.0, 0.0}, 0.0, 20.0, _dampedSHO.getSolution(Vector<Real>{1.0, 0.0}, 20.0) },
            
            // Population dynamics
            { "Logistic_t5", ODECategory::Population, &_logistic,
              Vector<Real>{1.0}, 0.0, 5.0, _logistic.getSolution(Vector<Real>{1.0}, 5.0) },
            
            { "Logistic_t10", ODECategory::Population, &_logistic,
              Vector<Real>{0.5}, 0.0, 10.0, _logistic.getSolution(Vector<Real>{0.5}, 10.0) },
            
            // Legacy compatibility
            { "TestLinODESys_t1", ODECategory::Linear, &_linear3D,
              Vector<Real>{1.0, 1.0, 1.0}, 0.0, 1.0, TestLinODESys_sol(1.0) },
            
            { "TestLinODESys_t5", ODECategory::Linear, &_linear3D,
              Vector<Real>{1.0, 1.0, 1.0}, 0.0, 5.0, TestLinODESys_sol(5.0) }
        };

        //--------------------------------------------------------------------------------------------------------------
        // SYSTEMS WITHOUT ANALYTICAL SOLUTIONS (for qualitative testing)
        //--------------------------------------------------------------------------------------------------------------
        
        struct ODESystemEntry
        {
            std::string name;
            ODECategory category;
            const IODESystem* system;
            Vector<Real> typicalIC;
        };
        
        static inline std::vector<ODESystemEntry> _generalSystems = {
            // Nonlinear oscillators
            { "VanDerPol_mu0.3", ODECategory::Nonlinear, &_vanDerPol, Vector<Real>{2.0, 0.0} },
            { "Duffing", ODECategory::Nonlinear, &_duffing, Vector<Real>{0.0, 1.0} },
            { "SimplePendulum", ODECategory::Nonlinear, &_pendulum, Vector<Real>{0.5, 0.0} },
            
            // Population dynamics
            { "LotkaVolterra", ODECategory::Population, &_lotkaVolterra, Vector<Real>{10.0, 5.0} },
            
            // Chaotic systems
            { "Lorenz", ODECategory::Chaotic, &_lorenz, Vector<Real>{1.0, 1.0, 1.0} },
            { "Rossler", ODECategory::Chaotic, &_rossler, Vector<Real>{0.1, 0.1, 0.1} },
            
            // Physical systems
            { "Projectile", ODECategory::Physical, &_projectile, Vector<Real>{0.0, 10.0, 0.0, 10.0} },
            { "Kepler", ODECategory::Physical, &_kepler, Vector<Real>{1.0, 0.0, 0.0, 1.0} },
            
            // Special functions
            { "Legendre_n2", ODECategory::SpecialFunction, &_legendre2, Vector<Real>{1.0, 0.0} },
            { "Laguerre_n3", ODECategory::SpecialFunction, &_laguerre3, Vector<Real>{1.0, 0.0} },
            { "Hermite_n2", ODECategory::SpecialFunction, &_hermite2, Vector<Real>{1.0, 0.0} },
            { "Bessel_n0", ODECategory::SpecialFunction, &_bessel0, Vector<Real>{1.0, 0.0} }
        };

        //--------------------------------------------------------------------------------------------------------------
        // SYSTEMS WITH JACOBIAN (for implicit methods)
        //--------------------------------------------------------------------------------------------------------------
        
        static inline std::vector<std::pair<std::string, const IODESystemWithJacobian*>> _systemsWithJacobian = {
            { "StiffExample1", &_stiff1 }
        };

    public:
        //--------------------------------------------------------------------------------------------------------------
        // COUNTS AND ITERATION
        //--------------------------------------------------------------------------------------------------------------
        
        static int numSystemsWithSolution() { return static_cast<int>(_systemsWithSolution.size()); }
        static int numSystemsWithEndSolution() { return static_cast<int>(_systemsWithEndSolution.size()); }
        static int numGeneralSystems() { return static_cast<int>(_generalSystems.size()); }
        static int numSystemsWithJacobian() { return static_cast<int>(_systemsWithJacobian.size()); }
        
        //--------------------------------------------------------------------------------------------------------------
        // ACCESSORS BY INDEX
        //--------------------------------------------------------------------------------------------------------------
        
        static const IODESystem* getSystemWithSolution(int index) 
        { 
            return _systemsWithSolution[index].system; 
        }
        
        static Vector<Real> getSolution(int index, const Vector<Real>& ic, Real t)
        {
            return _systemsWithSolution[index].solution(ic, t);
        }
        
        static const ODESystemWithEndSolutionEntry& getSystemWithEndSolution(int index)
        {
            return _systemsWithEndSolution[index];
        }
        
        static const ODESystemEntry& getGeneralSystem(int index)
        {
            return _generalSystems[index];
        }
        
        static const IODESystemWithJacobian* getSystemWithJacobian(int index)
        {
            return _systemsWithJacobian[index].second;
        }

        //--------------------------------------------------------------------------------------------------------------
        // ACCESSORS BY NAME
        //--------------------------------------------------------------------------------------------------------------
        
        static const IODESystem* getSystemWithSolutionByName(const std::string& name)
        {
            for (const auto& entry : _systemsWithSolution)
                if (entry.name == name)
                    return entry.system;
            throw std::runtime_error("System with solution not found: " + name);
        }
        
        static const ODESystemWithEndSolutionEntry* getSystemWithEndSolutionByName(const std::string& name)
        {
            for (const auto& entry : _systemsWithEndSolution)
                if (entry.name == name)
                    return &entry;
            throw std::runtime_error("System with end solution not found: " + name);
        }
        
        static const IODESystem* getGeneralSystemByName(const std::string& name)
        {
            for (const auto& entry : _generalSystems)
                if (entry.name == name)
                    return entry.system;
            throw std::runtime_error("General system not found: " + name);
        }

        //--------------------------------------------------------------------------------------------------------------
        // CATEGORY-BASED ACCESSORS
        //--------------------------------------------------------------------------------------------------------------
        
        static std::vector<const ODESystemWithSolutionEntry*> getSystemsWithSolutionByCategory(ODECategory cat)
        {
            std::vector<const ODESystemWithSolutionEntry*> result;
            for (const auto& entry : _systemsWithSolution)
                if (entry.category == cat)
                    result.push_back(&entry);
            return result;
        }
        
        static std::vector<const ODESystemWithEndSolutionEntry*> getSystemsWithEndSolutionByCategory(ODECategory cat)
        {
            std::vector<const ODESystemWithEndSolutionEntry*> result;
            for (const auto& entry : _systemsWithEndSolution)
                if (entry.category == cat)
                    result.push_back(&entry);
            return result;
        }
        
        static std::vector<const ODESystemEntry*> getGeneralSystemsByCategory(ODECategory cat)
        {
            std::vector<const ODESystemEntry*> result;
            for (const auto& entry : _generalSystems)
                if (entry.category == cat)
                    result.push_back(&entry);
            return result;
        }

        //--------------------------------------------------------------------------------------------------------------
        // FILTER HELPERS
        //--------------------------------------------------------------------------------------------------------------
        
        /** Get systems with solution filtered by predicate */
        template<typename Predicate>
        static std::vector<const ODESystemWithSolutionEntry*> filterSystemsWithSolution(Predicate pred)
        {
            std::vector<const ODESystemWithSolutionEntry*> result;
            for (const auto& entry : _systemsWithSolution)
                if (pred(entry))
                    result.push_back(&entry);
            return result;
        }
        
        /** Get systems with end solution filtered by predicate */
        template<typename Predicate>
        static std::vector<const ODESystemWithEndSolutionEntry*> filterSystemsWithEndSolution(Predicate pred)
        {
            std::vector<const ODESystemWithEndSolutionEntry*> result;
            for (const auto& entry : _systemsWithEndSolution)
                if (pred(entry))
                    result.push_back(&entry);
            return result;
        }
        
        //--------------------------------------------------------------------------------------------------------------
        // DIMENSION-BASED ACCESSORS
        //--------------------------------------------------------------------------------------------------------------
        
        static std::vector<const ODESystemWithEndSolutionEntry*> getSystemsWithEndSolutionByDimension(int dim)
        {
            return filterSystemsWithEndSolution([dim](const ODESystemWithEndSolutionEntry& e) {
                return e.system->getDim() == dim;
            });
        }
        
        static std::vector<const ODESystemEntry*> getGeneralSystemsByDimension(int dim)
        {
            std::vector<const ODESystemEntry*> result;
            for (const auto& entry : _generalSystems)
                if (entry.system->getDim() == dim)
                    result.push_back(&entry);
            return result;
        }

        //--------------------------------------------------------------------------------------------------------------
        // CONVENIENCE: GET ALL LINEAR SYSTEMS (best for initial solver testing)
        //--------------------------------------------------------------------------------------------------------------
        
        static std::vector<const ODESystemWithEndSolutionEntry*> getLinearSystems()
        {
            return getSystemsWithEndSolutionByCategory(ODECategory::Linear);
        }
        
        //--------------------------------------------------------------------------------------------------------------
        // CONVENIENCE: GET ALL OSCILLATORY SYSTEMS
        //--------------------------------------------------------------------------------------------------------------
        
        static std::vector<const ODESystemWithEndSolutionEntry*> getOscillatorySystems()
        {
            return getSystemsWithEndSolutionByCategory(ODECategory::Oscillatory);
        }
        
        //--------------------------------------------------------------------------------------------------------------
        // CONVENIENCE: GET ALL CHAOTIC SYSTEMS
        //--------------------------------------------------------------------------------------------------------------
        
        static std::vector<const ODESystemEntry*> getChaoticSystems()
        {
            return getGeneralSystemsByCategory(ODECategory::Chaotic);
        }
        
        //--------------------------------------------------------------------------------------------------------------
        // CONVENIENCE: GET ALL 1D SYSTEMS (simplest tests)
        //--------------------------------------------------------------------------------------------------------------
        
        static std::vector<const ODESystemWithEndSolutionEntry*> get1DSystems()
        {
            return getSystemsWithEndSolutionByDimension(1);
        }
        
        //--------------------------------------------------------------------------------------------------------------
        // CONVENIENCE: GET ALL 2D SYSTEMS 
        //--------------------------------------------------------------------------------------------------------------
        
        static std::vector<const ODESystemWithEndSolutionEntry*> get2DSystems()
        {
            return getSystemsWithEndSolutionByDimension(2);
        }

        //--------------------------------------------------------------------------------------------------------------
        // LEGACY COMPATIBILITY (old API)
        //--------------------------------------------------------------------------------------------------------------
        
        // Legacy static arrays for old API
        static inline ODESystem _legacySystems[3] = {
            ODESystem{3, TestLinODESys},
            ODESystem{2, VanDerPolMju0_1},
            ODESystem{2, [](Real t, const Vector<Real>& x, Vector<Real>& dxdt) { VanDerPol(0.1, t, x, dxdt); }}
        };
        
        // Legacy end solution systems (for backward compatibility with tests)
        static inline TestODESystemWithEndSolution _legacyEndSolSystems[4] = {
            TestODESystemWithEndSolution{3, TestLinODESys, Vector<Real>{1.0, 2.0, 2.0}, 1.0, TestLinODESys_sol(1.0)},
            TestODESystemWithEndSolution{3, TestLinODESys, Vector<Real>{1.0, 2.0, 2.0}, 5.0, TestLinODESys_sol(5.0)},
            TestODESystemWithEndSolution{3, TestLinODESys, Vector<Real>{1.0, 2.0, 2.0}, 10.0, TestLinODESys_sol(10.0)},
            TestODESystemWithEndSolution{new SimpleHarmonicOscillatorODE(), Vector<Real>{1.0, 2.0}, 10.0, 
                                         SimpleHarmonicOscillatorODE().getSolution(Vector<Real>{1.0, 2.0}, 10.0)}
        };
        
        [[deprecated("Use numGeneralSystems() instead")]]
        static int numODESystem() { return 3; }
        
        [[deprecated("Use getGeneralSystem(index).system instead")]]
        static const ODESystem& getODESystem(int index)
        {
            return _legacySystems[index];
        }
        
        [[deprecated("Use numSystemsWithJacobian() instead")]]
        static int numODESystemWithJacobian() { return 1; }
        
        [[deprecated("Use numSystemsWithSolution() instead")]]
        static int numODESystemWithSolution() { return numSystemsWithSolution(); }
        
        [[deprecated("Use numSystemsWithEndSolution() instead")]]
        static int numODESystemWithEndSolution() { return 4; }
        
        [[deprecated("Use getSystemWithEndSolution(index) instead")]]
        static const TestODESystemWithEndSolution& getTestODESystemWithEndSolution(int index)
        {
            return _legacyEndSolSystems[index];
        }
    };
}

#endif