#if !defined __MML_STIFF_ODE_TEST_BED_H
#define __MML_STIFF_ODE_TEST_BED_H

#include <string>
#include <vector>
#include <functional>
#include <memory>
#include <cmath>

#include "stiff_ode_defs.h"

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
    //                        STIFF ODE TEST DATA STRUCTURE                                   //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Represents a test case for stiff ODE solvers
     * 
     * Contains the system, initial conditions, and metadata for testing
     * stiff integrators like BDF, implicit Runge-Kutta, etc.
     */
    struct TestStiffODE
    {
        std::string name;                               ///< Descriptive name
        int dimension;                                  ///< System dimension
        
        // Factory functions
        std::function<std::unique_ptr<IODESystem>()> createSystem;           ///< For non-stiff interface
        std::function<std::unique_ptr<IODESystemWithJacobian>()> createStiffSystem;  ///< For stiff solvers
        
        // Problem setup
        Vector<Real> initialCondition;                 ///< y(0)
        Real tStart = 0.0;                             ///< Integration start time
        Real tEnd = 1.0;                               ///< Integration end time
        
        // Reference solution (if known)
        bool hasExactSolution = false;                 ///< Whether exact solution exists
        std::function<Vector<Real>(Real)> exactSolution;  ///< y(t) if known
        
        // Stiffness characteristics
        Real stiffnessRatio = 1.0;                     ///< Ratio of fastest/slowest eigenvalue
        std::string stiffnessCategory;                  ///< mild, moderate, severe, extreme
        
        // Physical interpretation
        std::string description;                       ///< What the system models
        std::string reference;                         ///< Literature reference
        
        // Testing hints
        Real suggestedStepSize = 0.01;                 ///< Reasonable starting step
        Real expectedAccuracy = 1e-6;                  ///< Target accuracy for adaptive
        int difficulty = 2;                            ///< 1=easy, 2=medium, 3=hard, 4=extreme
    };

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                            TEST CASE GENERATORS                                         //
    ///////////////////////////////////////////////////////////////////////////////////////////

    inline TestStiffODE getRobertsonTest()
    {
        TestStiffODE test;
        test.name = "Robertson";
        test.dimension = 3;
        test.createStiffSystem = []() { return std::make_unique<RobertsonODE>(); };
        test.initialCondition = RobertsonODE::getInitialCondition();
        test.tStart = 0.0;
        test.tEnd = 1e11;  // Extremely long integration
        test.stiffnessRatio = 1e8;
        test.stiffnessCategory = "extreme";
        test.description = RobertsonODE::getDescription();
        test.reference = "Robertson, H.H. (1966)";
        test.suggestedStepSize = 1e-6;
        test.expectedAccuracy = 1e-6;
        test.difficulty = 4;
        return test;
    }

    inline TestStiffODE getVanDerPolStiffTest(Real mu = 1000.0)
    {
        TestStiffODE test;
        test.name = "VanDerPol_mu" + std::to_string(static_cast<int>(mu));
        test.dimension = 2;
        test.createStiffSystem = [mu]() { return std::make_unique<VanDerPolStiffODE>(mu); };
        test.initialCondition = VanDerPolStiffODE::getInitialCondition();
        test.tStart = 0.0;
        test.tEnd = 2.0 * mu;  // About one period
        test.stiffnessRatio = mu * mu;
        test.stiffnessCategory = (mu > 100) ? "severe" : "moderate";
        test.description = VanDerPolStiffODE::getDescription();
        test.reference = "Hairer & Wanner, Solving ODEs II";
        test.suggestedStepSize = 0.01 / mu;
        test.expectedAccuracy = 1e-6;
        test.difficulty = (mu > 100) ? 3 : 2;
        return test;
    }

    inline TestStiffODE getBrusselatorTest(Real A = 1.0, Real B = 3.0)
    {
        TestStiffODE test;
        test.name = "Brusselator";
        test.dimension = 2;
        test.createStiffSystem = [A, B]() { return std::make_unique<BrusselatorODE>(A, B); };
        test.initialCondition = BrusselatorODE::getInitialCondition();
        test.tStart = 0.0;
        test.tEnd = 20.0;
        test.stiffnessRatio = B + 1 + A*A;  // Approximate
        test.stiffnessCategory = "moderate";
        test.description = BrusselatorODE::getDescription();
        test.reference = "Prigogine & Lefever";
        test.suggestedStepSize = 0.01;
        test.expectedAccuracy = 1e-6;
        test.difficulty = 2;
        return test;
    }

    inline TestStiffODE getOregonatorTest()
    {
        TestStiffODE test;
        test.name = "Oregonator";
        test.dimension = 3;
        test.createStiffSystem = []() { return std::make_unique<OregonatorODE>(); };
        test.initialCondition = OregonatorODE::getInitialCondition();
        test.tStart = 0.0;
        test.tEnd = 360.0;
        test.stiffnessRatio = 77.27 / 8.375e-6;  // s/q ≈ 10^7
        test.stiffnessCategory = "severe";
        test.description = OregonatorODE::getDescription();
        test.reference = "Field & Noyes (1974)";
        test.suggestedStepSize = 0.01;
        test.expectedAccuracy = 1e-5;
        test.difficulty = 3;
        return test;
    }

    inline TestStiffODE getHIRESTest()
    {
        TestStiffODE test;
        test.name = "HIRES";
        test.dimension = 8;
        test.createStiffSystem = []() { return std::make_unique<HIRES_ODE>(); };
        test.initialCondition = HIRES_ODE::getInitialCondition();
        test.tStart = 0.0;
        test.tEnd = HIRES_ODE::getTypicalEndTime();
        test.stiffnessRatio = 280.0 / 0.43;  // Approximate
        test.stiffnessCategory = "moderate";
        test.description = HIRES_ODE::getDescription();
        test.reference = "Hairer Test Set for IVP Solvers";
        test.suggestedStepSize = 0.001;
        test.expectedAccuracy = 1e-6;
        test.difficulty = 2;
        return test;
    }

    inline TestStiffODE getRobertsonScaledTest()
    {
        TestStiffODE test;
        test.name = "Robertson_Scaled";
        test.dimension = 2;
        test.createStiffSystem = []() { return std::make_unique<RobertsonScaledODE>(); };
        test.initialCondition = RobertsonScaledODE::getInitialCondition();
        test.tStart = 0.0;
        test.tEnd = 1e4;
        test.stiffnessRatio = 1e8;
        test.stiffnessCategory = "extreme";
        test.description = RobertsonScaledODE::getDescription();
        test.reference = "Robertson (conserved form)";
        test.suggestedStepSize = 1e-6;
        test.expectedAccuracy = 1e-6;
        test.difficulty = 4;
        return test;
    }

    inline TestStiffODE getE5Test()
    {
        TestStiffODE test;
        test.name = "E5_Decay";
        test.dimension = 4;
        test.createStiffSystem = []() { return std::make_unique<E5_ODE>(); };
        test.initialCondition = E5_ODE::getInitialCondition();
        test.tStart = 0.0;
        test.tEnd = E5_ODE::getTypicalEndTime();
        test.stiffnessRatio = 7.89e10 / 13.0;  // k1/k4 ≈ 6×10^9
        test.stiffnessCategory = "extreme";
        test.description = E5_ODE::getDescription();
        test.reference = "Enright et al. Test Set";
        test.suggestedStepSize = 1e-12;
        test.expectedAccuracy = 1e-6;
        test.difficulty = 4;
        return test;
    }

    inline TestStiffODE getPollutionTest()
    {
        TestStiffODE test;
        test.name = "Pollution";
        test.dimension = 10;
        test.createStiffSystem = []() { return std::make_unique<PollutionODE>(); };
        test.initialCondition = PollutionODE::getInitialCondition();
        test.tStart = 0.0;
        test.tEnd = 60.0;
        test.stiffnessRatio = 0.123e5 / 0.82e-3;  // Approximate
        test.stiffnessCategory = "moderate";
        test.description = PollutionODE::getDescription();
        test.reference = "Verwer Test Set";
        test.suggestedStepSize = 0.001;
        test.expectedAccuracy = 1e-5;
        test.difficulty = 3;
        return test;
    }

    inline TestStiffODE getLinearStiffTest()
    {
        TestStiffODE test;
        test.name = "Linear_Stiff";
        test.dimension = 3;
        test.createStiffSystem = []() { return std::make_unique<LinearStiffODE>(); };
        test.initialCondition = LinearStiffODE::getInitialCondition();
        test.tStart = 0.0;
        test.tEnd = 1.0;
        test.hasExactSolution = true;
        test.exactSolution = [](Real t) { return LinearStiffODE::getExactSolution(t); };
        test.stiffnessRatio = 10000.0;
        test.stiffnessCategory = "severe";
        test.description = LinearStiffODE::getDescription();
        test.reference = "Standard test problem";
        test.suggestedStepSize = 1e-4;
        test.expectedAccuracy = 1e-8;
        test.difficulty = 2;
        return test;
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                        CATEGORY RETRIEVAL FUNCTIONS                                    //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Get all stiff ODE test cases
     */
    inline std::vector<TestStiffODE> getAllStiffODETests()
    {
        return {
            getRobertsonTest(),
            getVanDerPolStiffTest(1000.0),
            getBrusselatorTest(),
            getOregonatorTest(),
            getHIRESTest(),
            getRobertsonScaledTest(),
            getE5Test(),
            getPollutionTest(),
            getLinearStiffTest()
        };
    }

    /**
     * @brief Get moderate stiffness tests (good starting point)
     */
    inline std::vector<TestStiffODE> getModerateStiffTests()
    {
        return {
            getBrusselatorTest(),
            getHIRESTest(),
            getLinearStiffTest()
        };
    }

    /**
     * @brief Get severe/extreme stiffness tests (challenging)
     */
    inline std::vector<TestStiffODE> getSevereStiffTests()
    {
        return {
            getRobertsonTest(),
            getVanDerPolStiffTest(1000.0),
            getOregonatorTest(),
            getE5Test()
        };
    }

    /**
     * @brief Get tests with known exact solutions
     */
    inline std::vector<TestStiffODE> getExactSolutionStiffTests()
    {
        return {
            getLinearStiffTest()
        };
    }

    /**
     * @brief Get chemical kinetics tests
     */
    inline std::vector<TestStiffODE> getChemicalKineticsTests()
    {
        return {
            getRobertsonTest(),
            getRobertsonScaledTest(),
            getBrusselatorTest(),
            getOregonatorTest(),
            getE5Test()
        };
    }

    /**
     * @brief Get smaller dimension tests (2-3D)
     */
    inline std::vector<TestStiffODE> getSmallDimensionStiffTests()
    {
        return {
            getVanDerPolStiffTest(1000.0),
            getBrusselatorTest(),
            getOregonatorTest(),
            getLinearStiffTest()
        };
    }

    /**
     * @brief Get larger dimension tests (8-10D)
     */
    inline std::vector<TestStiffODE> getLargeDimensionStiffTests()
    {
        return {
            getHIRESTest(),
            getPollutionTest()
        };
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                            VERIFICATION UTILITIES                                       //
    ///////////////////////////////////////////////////////////////////////////////////////////

    /**
     * @brief Compute error against exact solution
     * @param computed Computed solution y(t)
     * @param test Test case with exact solution
     * @param t Time point
     * @return L2 norm of error, or -1 if no exact solution
     */
    inline Real computeStiffODEError(const Vector<Real>& computed, const TestStiffODE& test, Real t)
    {
        if (!test.hasExactSolution) return -1.0;
        
        Vector<Real> exact = test.exactSolution(t);
        Real error = 0.0;
        for (int i = 0; i < test.dimension; ++i)
            error += (computed[i] - exact[i]) * (computed[i] - exact[i]);
        return std::sqrt(error);
    }

    /**
     * @brief Check conservation laws (e.g., mass conservation for Robertson)
     */
    inline Real checkConservation(const Vector<Real>& y, const TestStiffODE& test)
    {
        if (test.name.find("Robertson") != std::string::npos && y.size() == 3) {
            // y1 + y2 + y3 should equal 1 (mass conservation)
            return std::abs(y[0] + y[1] + y[2] - 1.0);
        }
        return -1.0;  // No conservation law checked
    }

    /**
     * @brief Verify non-negativity (for chemical kinetics)
     */
    inline bool verifyNonNegative(const Vector<Real>& y)
    {
        for (int i = 0; i < y.size(); ++i)
            if (y[i] < -1e-10) return false;  // Allow small numerical errors
        return true;
    }

} // namespace MML::TestBeds

#endif // __MML_STIFF_ODE_TEST_BED_H
