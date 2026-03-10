///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        readme15_dynamical_systems.cpp                                      ///
///  Description: README example - Dynamical Systems Analysis (Crown Jewel!)         ///
///               Demonstrates chaos theory, Lyapunov exponents, fixed points,        ///
///               bifurcation analysis, and classic chaotic attractors                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#include "systems/DynamicalSystem.h"
#endif

#include <iostream>
#include <iomanip>

using namespace MML;
using namespace MML::Systems;

void Readme_DynamicalSystems()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****        DYNAMICAL SYSTEMS ANALYSIS - THE CROWN JEWEL          ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    //=========================================================================
    // 1. THE LORENZ SYSTEM - CHAOS THEORY ICON
    //=========================================================================
    std::cout << "\n=== The Lorenz System ===" << std::endl;
    std::cout << "The butterfly effect: dx/dt = sigma(y-x), dy/dt = x(rho-z)-y, dz/dt = xy-beta*z\n" << std::endl;

    // Create Lorenz system with classic parameters: sigma=10, rho=28, beta=8/3
    LorenzSystem lorenz(10.0, 28.0, 8.0/3.0);
    
    std::cout << "Parameters: sigma=" << lorenz.getParam(0) 
              << ", rho=" << lorenz.getParam(1) 
              << ", beta=" << std::setprecision(4) << lorenz.getParam(2) << std::endl;
    std::cout << "Dissipative: " << (lorenz.isDissipative() ? "yes" : "no") << std::endl;
    std::cout << "Flow divergence: " << lorenz.getDivergence() << " (always negative -> attractor exists)" << std::endl;

    //=========================================================================
    // 2. FIXED POINTS - WHERE CHAOS DOESN'T HAPPEN
    //=========================================================================
    std::cout << "\n=== Fixed Point Analysis ===" << std::endl;
    
    // Lorenz has 3 fixed points: origin and two symmetric ones
    std::vector<Vector<Real>> guesses = {
        Vector<Real>{0.0, 0.0, 0.0},                    // Origin
        Vector<Real>{8.0, 8.0, 27.0},                   // Near C+
        Vector<Real>{-8.0, -8.0, 27.0}                  // Near C-
    };
    
    auto fixedPoints = FixedPointFinder::FindMultiple(lorenz, guesses);
    
    std::cout << "Found " << fixedPoints.size() << " fixed points:\n" << std::endl;
    for (size_t i = 0; i < fixedPoints.size(); ++i) {
        const auto& fp = fixedPoints[i];
        std::cout << "  Fixed Point " << (i+1) << ": (" 
                  << std::setprecision(4) << fp.location[0] << ", " 
                  << fp.location[1] << ", " 
                  << fp.location[2] << ")" << std::endl;
        std::cout << "    Type: " << ToString(fp.type) << std::endl;
        std::cout << "    Stable: " << (fp.isStable ? "yes" : "NO (unstable)") << std::endl;
        std::cout << "    Eigenvalues: ";
        for (const auto& ev : fp.eigenvalues)
            std::cout << "(" << std::setprecision(3) << ev.real() << "+" << ev.imag() << "i) ";
        std::cout << std::endl;
    }

    //=========================================================================
    // 3. LYAPUNOV EXPONENTS - QUANTIFYING CHAOS
    //=========================================================================
    std::cout << "\n=== Lyapunov Exponent Analysis ===" << std::endl;
    std::cout << "Measuring sensitivity to initial conditions...\n" << std::endl;
    
    Vector<Real> x0 = lorenz.getDefaultInitialCondition();
    
    // Compute Lyapunov exponents (takes a few seconds for accuracy)
    auto lyapResult = LyapunovAnalyzer::Compute(lorenz, x0, 
                                                 500.0,   // Total integration time
                                                 1.0,     // Orthonormalization interval
                                                 0.01);   // Step size

    std::cout << std::setprecision(6);
    std::cout << "Lyapunov spectrum: [" 
              << lyapResult.exponents[0] << ", "
              << lyapResult.exponents[1] << ", "
              << lyapResult.exponents[2] << "]" << std::endl;
    std::cout << "Maximum exponent: " << lyapResult.maxExponent;
    if (lyapResult.maxExponent > 0)
        std::cout << " (POSITIVE -> CHAOS!)";
    std::cout << std::endl;
    std::cout << "Sum of exponents: " << lyapResult.sum << " (negative for dissipative)" << std::endl;
    std::cout << "Kaplan-Yorke dimension: " << lyapResult.kaplanYorkeDimension 
              << " (fractal attractor!)" << std::endl;
    std::cout << "System is " << (lyapResult.isChaotic ? "CHAOTIC" : "regular") << std::endl;

    //=========================================================================
    // 4. COMPARING CLASSIC CHAOTIC SYSTEMS
    //=========================================================================
    std::cout << "\n=== Classic Chaotic Systems ===" << std::endl;
    
    // Rossler - simpler spiral chaos
    RosslerSystem rossler(0.2, 0.2, 5.7);
    auto rosslerLyap = LyapunovAnalyzer::Compute(rossler, rossler.getDefaultInitialCondition(), 
                                                  300.0, 1.0, 0.01);
    
    // Van der Pol - limit cycle oscillator (not chaotic!)
    VanDerPolSystem vanderpol(1.0);
    auto vdpLyap = LyapunovAnalyzer::Compute(vanderpol, vanderpol.getDefaultInitialCondition(),
                                              200.0, 1.0, 0.01);

    std::cout << std::setprecision(4);
    std::cout << "System          | Max Lyapunov | K-Y Dimension | Chaotic?" << std::endl;
    std::cout << "----------------|--------------|---------------|----------" << std::endl;
    std::cout << "Lorenz          | " << std::setw(12) << lyapResult.maxExponent 
              << " | " << std::setw(13) << lyapResult.kaplanYorkeDimension 
              << " | " << (lyapResult.maxExponent > 0.01 ? "YES" : "no") << std::endl;
    std::cout << "Rossler         | " << std::setw(12) << rosslerLyap.maxExponent 
              << " | " << std::setw(13) << rosslerLyap.kaplanYorkeDimension 
              << " | " << (rosslerLyap.maxExponent > 0.01 ? "YES" : "no") << std::endl;
    std::cout << "Van der Pol     | " << std::setw(12) << vdpLyap.maxExponent 
              << " | " << std::setw(13) << vdpLyap.kaplanYorkeDimension 
              << " | " << (vdpLyap.maxExponent > 0.01 ? "YES" : "no") << " (limit cycle)" << std::endl;

    //=========================================================================
    // 5. DOUBLE PENDULUM - MECHANICAL CHAOS
    //=========================================================================
    std::cout << "\n=== Double Pendulum - Mechanical Chaos ===" << std::endl;
    
    DoublePendulumSystem pendulum(1.0, 1.0, 9.81);  // m=1, L=1, g=9.81
    
    Vector<Real> ic = pendulum.getDefaultInitialCondition();
    std::cout << "Initial: theta1=" << std::setprecision(3) << ic[0] 
              << ", theta2=" << ic[1] << " rad" << std::endl;
    
    // Check energy conservation (Hamiltonian system)
    Real E0 = pendulum.computeInvariant(0, ic);
    std::cout << "Initial energy: " << std::setprecision(6) << E0 << " J" << std::endl;
    
    auto pendLyap = LyapunovAnalyzer::Compute(pendulum, ic, 200.0, 1.0, 0.005);
    std::cout << "Max Lyapunov exponent: " << pendLyap.maxExponent << std::endl;
    std::cout << "Chaos confirmed: " << (pendLyap.isChaotic ? "YES!" : "no") << std::endl;

    //=========================================================================
    // 6. BIFURCATION ANALYSIS - ROUTE TO CHAOS
    //=========================================================================
    std::cout << "\n=== Bifurcation Analysis ===" << std::endl;
    std::cout << "Sweeping Lorenz rho parameter from 20 to 30...\n" << std::endl;
    
    LorenzSystem lorenzSweep;
    Vector<Real> sweepIC{1.0, 1.0, 1.0};
    
    // Quick sweep to show structure
    auto bifurcation = BifurcationAnalyzer::Sweep(
        lorenzSweep,
        1,              // Parameter index (rho)
        20.0, 30.0,     // Parameter range
        6,              // Number of parameter values
        sweepIC,
        2,              // Record z-component maxima
        50.0,           // Transient time
        20.0,           // Recording time
        0.01            // Step size
    );
    
    std::cout << "Attractor structure at different rho values:" << std::endl;
    for (size_t i = 0; i < bifurcation.parameterValues.size(); ++i) {
        std::cout << "  rho = " << std::setprecision(1) << std::fixed 
                  << bifurcation.parameterValues[i] << ": ";
        const auto& maxima = bifurcation.attractorValues[i];
        if (maxima.size() <= 2)
            std::cout << "periodic (limit cycle)";
        else if (maxima.size() <= 4)
            std::cout << "period-doubled";
        else
            std::cout << "CHAOTIC (many local maxima)";
        std::cout << " [" << maxima.size() << " maxima]" << std::endl;
    }

    std::cout << "\n*** Dynamical Systems Analysis Complete ***" << std::endl;
    std::cout << std::endl;
}
