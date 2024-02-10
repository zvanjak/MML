#if !defined MML_ODE_SYSTEM_TEST_BED_H
#define MML_ODE_SYSTEM_TEST_BED_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/Matrix.h"
#include "core/ODESystem.h"
#endif

#include "diff_eq_systems_defs.h"

namespace MML::TestBeds
{
    // TODO - extend with at least 10 ODE systems with known solutionas
    class TestODESystemWithSolution : public ODESystem
    {
        Vector<Real> _initial_conditions;
        Vector<Real> (*_solution)(Real);

    public:
        TestODESystemWithSolution(int n, 
                                  void (*inFunc)(Real, const Vector<Real>&, Vector<Real> &),  
                                  Vector<Real> init_cond,
                                  Vector<Real> (*inSol)(Real) ) : ODESystem(n, inFunc), _initial_conditions(init_cond), _solution(inSol) { }

        Vector<Real> getSolution(Real t)
        {
            return _solution(t);
        }
        Vector<Real> getInitialConditions() const
        {
            return _initial_conditions;
        }
    };

    class TestODESystemWithEndSolution : public ODESystem
    {
        Vector<Real> _initial_conditions;
        Real       _start_time;
        Real       _end_time;
        Vector<Real> _end_solution;

    public:
        TestODESystemWithEndSolution(int n, 
                                     void (*inFunc)(Real, const Vector<Real>&, Vector<Real> &),  
                                     Vector<Real> init_cond, Real start_time, Real end_time,
                                     Vector<Real> end_solution ) 
                                        : ODESystem(n, inFunc), 
                                          _initial_conditions(init_cond), _start_time(start_time), _end_time(end_time), _end_solution(end_solution) 
        { }

        Vector<Real> getInitialConditions() const
        {
            return _initial_conditions;
        }
    };

    class ODESystemTestBed
    {
    private:
       const static inline std::pair<std::string, ODESystem> _listODESystems[] = 
            { 
                { "simple 1", ODESystem{3, fnc12} },
                { "simple 2", ODESystem{2, VanDerPolMju0_1 } },
                { "VanDerPol 0.1", ODESystem{2, [](Real t, const Vector<Real> &x, Vector<Real> &dxdt) { return VanDerPol(0.1, t, x, dxdt); } } }
            };
        const static inline std::pair<std::string, ODESystemWithJacobian> _listODESystemsWithJac[] = 
            { 
                { "stiff 1", ODESystemWithJacobian{3, stiff_sys1_derivs, stiff_sys1_jac} }
            };
        const static inline std::pair<std::string, TestODESystemWithSolution> _listODESystemsWithSol[] = 
            { 
                { "simple 1", TestODESystemWithSolution{3, fnc12, {1, 2, 2}, fnc12_sol} },
                { "simple 2", TestODESystemWithSolution{3, 
                                          [](Real t, const Vector<Real> &x, Vector<Real> &dxdt)  { dxdt = Vector<Real>{ x[0] + x[1] - x[2], -x[0] + 3*x[1] - x[2], -x[0] + x[1] + x[2] }; },
                                          Vector<Real>{1, 2, 2}, 
                                          [](Real t)  { return Vector<Real>{ exp(t), exp(t) + exp(2*t), exp(t) + exp(2*t) }; } } }
            };        
        const static inline std::pair<std::string, TestODESystemWithEndSolution> _listODESystemsWithEndSol[] = 
            { 
                { "simple", TestODESystemWithEndSolution{3, fnc12, {1, 2, 2}, 0.0, 10.0, fnc12_sol(10.0)} },
            };          
    public:
        static int   numODESystem() { return 3; }     
        static const ODESystem& getODESystem(int index) { return _listODESystems[index].second; }
        static const ODESystem& getODESystem(std::string sysName)
        {
            for (auto &sys : _listODESystems)
            {
                if (sys.first == sysName)
                    return sys.second;
            }
            throw std::runtime_error("ODESystemTestBed::getODESystem: system with name " + sysName + " not found");
        }

        static int   numODESystemWithJacobian() { return 1; }
        static const ODESystemWithJacobian& getODESystemWithJacobian(int index) { return _listODESystemsWithJac[index].second; }
        static const ODESystemWithJacobian& getODESystemWithJacobian(std::string sysName)
        {
            for (auto &sys : _listODESystemsWithJac)
            {
                if (sys.first == sysName)
                    return sys.second;
            }
            throw std::runtime_error("ODESystemTestBed::getODESystemWithJacobian: system with name " + sysName + " not found");
        }

        static int   numODESystemWithSolution() { return 2; }
        static const TestODESystemWithSolution& getTestODESystemWithSolution(int index) { return _listODESystemsWithSol[index].second; }
        static const TestODESystemWithSolution& getTestODESystemWithSolution(std::string sysName)
        {
            for (auto &sys : _listODESystemsWithSol)
            {
                if (sys.first == sysName)
                    return sys.second;
            }
            throw std::runtime_error("ODESystemTestBed::getTestODESystemWithSolution: system with name " + sysName + " not found");
        }

        static int   numODESystemWithEndSolution() { return 1; }
        static const TestODESystemWithEndSolution& getTestODESystemWithEndSolution(int index) { return _listODESystemsWithEndSol[index].second; }
        static const TestODESystemWithEndSolution& getTestODESystemWithEndSolution(std::string sysName)
        {
            for (auto &sys : _listODESystemsWithEndSol)
            {
                if (sys.first == sysName)
                    return sys.second;
            }
            throw std::runtime_error("ODESystemTestBed::getTestODESystemWithEndSolution: system with name " + sysName + " not found");
        } 
    };
}
#endif