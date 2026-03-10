///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        docs_demo_random.cpp                                                ///
///  Description: Documentation examples for Random class                             ///
///               These examples are referenced from docs/base/Random.md              ///
///                                                                                   ///
///  Usage:       Run MML_DocsDemo application to execute these examples              ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Random.h"
#include "base/Vector/Vector.h"
#endif

#include <iostream>
#include <iomanip>
#include <map>
#include <cmath>

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                           Uniform Real Distribution                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Random_UniformReal()
{
    std::cout << "\n=== Uniform Real Distribution ===\n\n";

    // Generate random numbers in [0, 1)
    std::cout << "Random reals in [0, 1):" << std::endl;
    std::cout << "  ";
    for (int i = 0; i < 10; ++i) {
        std::cout << std::fixed << std::setprecision(4) 
                  << Random::UniformReal(0.0, 1.0) << " ";
    }
    std::cout << std::endl;

    // Generate random numbers in custom range
    std::cout << "\nRandom reals in [-5, 5):" << std::endl;
    std::cout << "  ";
    for (int i = 0; i < 10; ++i) {
        std::cout << std::fixed << std::setprecision(2) 
                  << Random::UniformReal(-5.0, 5.0) << " ";
    }
    std::cout << std::endl;

    // Statistical verification - should average near 0.5 for [0,1)
    std::cout << "\nStatistical verification (1000 samples in [0, 1)):" << std::endl;
    double sum = 0.0;
    double minVal = 1.0, maxVal = 0.0;
    for (int i = 0; i < 1000; ++i) {
        double val = Random::UniformReal(0.0, 1.0);
        sum += val;
        minVal = std::min(minVal, val);
        maxVal = std::max(maxVal, val);
    }
    std::cout << "  Mean:     " << std::fixed << std::setprecision(4) 
              << (sum / 1000.0) << " (expected: 0.5)" << std::endl;
    std::cout << "  Min:      " << minVal << std::endl;
    std::cout << "  Max:      " << maxVal << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Uniform Integer Distribution                              ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Random_UniformInt()
{
    std::cout << "\n=== Uniform Integer Distribution ===\n\n";

    // Generate random integers in [1, 6] (dice roll simulation)
    std::cout << "Dice rolls (10 samples):" << std::endl;
    std::cout << "  ";
    for (int i = 0; i < 10; ++i) {
        std::cout << Random::UniformInt(1, 6) << " ";
    }
    std::cout << std::endl;

    // Histogram verification
    std::cout << "\nDice roll histogram (1000 rolls):" << std::endl;
    std::map<int, int> histogram;
    for (int i = 1; i <= 6; ++i) histogram[i] = 0;
    
    for (int i = 0; i < 1000; ++i) {
        histogram[Random::UniformInt(1, 6)]++;
    }
    
    for (int i = 1; i <= 6; ++i) {
        std::cout << "  " << i << ": ";
        int bars = histogram[i] / 10;  // Scale for display
        for (int j = 0; j < bars; ++j) std::cout << "*";
        std::cout << " (" << histogram[i] << ")" << std::endl;
    }

    // Random indices (common use case)
    std::cout << "\nRandom array indices (array size 5):" << std::endl;
    std::cout << "  ";
    for (int i = 0; i < 10; ++i) {
        std::cout << Random::UniformInt(0, 4) << " ";
    }
    std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           2D Random Direction                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Random_Direction2D()
{
    std::cout << "\n=== 2D Random Direction (Uniform on Circle) ===\n\n";

    // Generate unit vectors on the unit circle
    std::cout << "Random unit directions (|v| = 1.0):" << std::endl;
    for (int i = 0; i < 5; ++i) {
        Real vx, vy;
        Random::UniformVecDirection2(vx, vy, 1.0);
        Real length = std::sqrt(vx*vx + vy*vy);
        std::cout << "  (" << std::setw(8) << std::fixed << std::setprecision(4) << vx 
                  << ", " << std::setw(8) << vy << ")  |v| = " << length << std::endl;
    }

    // Generate vectors with custom magnitude
    std::cout << "\nRandom directions with |v| = 5.0:" << std::endl;
    for (int i = 0; i < 3; ++i) {
        Real vx, vy;
        Random::UniformVecDirection2(vx, vy, 5.0);
        Real length = std::sqrt(vx*vx + vy*vy);
        std::cout << "  (" << std::setw(8) << std::fixed << std::setprecision(4) << vx 
                  << ", " << std::setw(8) << vy << ")  |v| = " << length << std::endl;
    }

    // Verify uniform distribution on circle (check angle distribution)
    std::cout << "\nAngle distribution verification (1000 samples):" << std::endl;
    int quadrants[4] = {0, 0, 0, 0};
    for (int i = 0; i < 1000; ++i) {
        Real vx, vy;
        Random::UniformVecDirection2(vx, vy, 1.0);
        
        // Determine quadrant
        if (vx >= 0 && vy >= 0) quadrants[0]++;
        else if (vx < 0 && vy >= 0) quadrants[1]++;
        else if (vx < 0 && vy < 0) quadrants[2]++;
        else quadrants[3]++;
    }
    std::cout << "  Q1 (+,+): " << quadrants[0] << " (expected ~250)" << std::endl;
    std::cout << "  Q2 (-,+): " << quadrants[1] << " (expected ~250)" << std::endl;
    std::cout << "  Q3 (-,-): " << quadrants[2] << " (expected ~250)" << std::endl;
    std::cout << "  Q4 (+,-): " << quadrants[3] << " (expected ~250)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           3D Random Direction                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Random_Direction3D()
{
    std::cout << "\n=== 3D Random Direction (Uniform on Sphere) ===\n\n";

    // Generate unit vectors on the unit sphere
    std::cout << "Random unit directions (|v| = 1.0):" << std::endl;
    for (int i = 0; i < 5; ++i) {
        Real vx, vy, vz;
        Random::UniformVecDirection3(vx, vy, vz, 1.0);
        Real length = std::sqrt(vx*vx + vy*vy + vz*vz);
        std::cout << "  (" << std::setw(8) << std::fixed << std::setprecision(4) << vx 
                  << ", " << std::setw(8) << vy 
                  << ", " << std::setw(8) << vz << ")  |v| = " << length << std::endl;
    }

    // Generate vectors with custom magnitude
    std::cout << "\nRandom directions with |v| = 10.0:" << std::endl;
    for (int i = 0; i < 3; ++i) {
        Real vx, vy, vz;
        Random::UniformVecDirection3(vx, vy, vz, 10.0);
        Real length = std::sqrt(vx*vx + vy*vy + vz*vz);
        std::cout << "  (" << std::setw(8) << std::fixed << std::setprecision(4) << vx 
                  << ", " << std::setw(8) << vy 
                  << ", " << std::setw(8) << vz << ")  |v| = " << length << std::endl;
    }

    // Verify uniform distribution on sphere (z-coordinate should be uniform in [-1,1])
    std::cout << "\nZ-coordinate distribution (uniform sphere test, 1000 samples):" << std::endl;
    std::cout << "  (For uniform sphere, z should be uniformly distributed)" << std::endl;
    
    int zBins[5] = {0, 0, 0, 0, 0};  // z in [-1,-0.6), [-0.6,-0.2), [-0.2,0.2), [0.2,0.6), [0.6,1]
    for (int i = 0; i < 1000; ++i) {
        Real vx, vy, vz;
        Random::UniformVecDirection3(vx, vy, vz, 1.0);
        
        if (vz < -0.6) zBins[0]++;
        else if (vz < -0.2) zBins[1]++;
        else if (vz < 0.2) zBins[2]++;
        else if (vz < 0.6) zBins[3]++;
        else zBins[4]++;
    }
    std::cout << "  z in [-1.0, -0.6): " << zBins[0] << " (expected ~200)" << std::endl;
    std::cout << "  z in [-0.6, -0.2): " << zBins[1] << " (expected ~200)" << std::endl;
    std::cout << "  z in [-0.2,  0.2): " << zBins[2] << " (expected ~200)" << std::endl;
    std::cout << "  z in [ 0.2,  0.6): " << zBins[3] << " (expected ~200)" << std::endl;
    std::cout << "  z in [ 0.6,  1.0]: " << zBins[4] << " (expected ~200)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Practical Use Cases                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Random_UseCases()
{
    std::cout << "\n=== Practical Use Cases ===\n\n";

    // Monte Carlo integration example: estimate π
    std::cout << "Monte Carlo π estimation (10000 samples):" << std::endl;
    int insideCircle = 0;
    const int samples = 10000;
    for (int i = 0; i < samples; ++i) {
        Real x = Random::UniformReal(-1.0, 1.0);
        Real y = Random::UniformReal(-1.0, 1.0);
        if (x*x + y*y <= 1.0) {
            insideCircle++;
        }
    }
    Real piEstimate = 4.0 * static_cast<Real>(insideCircle) / samples;
    std::cout << "  Estimated π: " << std::fixed << std::setprecision(6) << piEstimate << std::endl;
    std::cout << "  Actual π:    " << Constants::PI << std::endl;
    std::cout << "  Error:       " << std::abs(piEstimate - Constants::PI) << std::endl;

    // Random walk in 2D
    std::cout << "\nRandom walk in 2D (100 steps, step size = 1):" << std::endl;
    Real x = 0.0, y = 0.0;
    for (int step = 0; step < 100; ++step) {
        Real dx, dy;
        Random::UniformVecDirection2(dx, dy, 1.0);
        x += dx;
        y += dy;
    }
    Real distFromOrigin = std::sqrt(x*x + y*y);
    std::cout << "  Final position: (" << std::fixed << std::setprecision(2) 
              << x << ", " << y << ")" << std::endl;
    std::cout << "  Distance from origin: " << distFromOrigin << std::endl;
    std::cout << "  (Expected RMS distance: ~10 = sqrt(100))" << std::endl;

    // Random sampling from vector
    std::cout << "\nRandom element selection from vector:" << std::endl;
    Vector<double> data({10.0, 20.0, 30.0, 40.0, 50.0});
    std::cout << "  Data: " << data << std::endl;
    std::cout << "  Random samples: ";
    for (int i = 0; i < 5; ++i) {
        int idx = Random::UniformInt(0, data.size() - 1);
        std::cout << data[idx] << " ";
    }
    std::cout << std::endl;

    // Shuffle simulation (Fisher-Yates would use this)
    std::cout << "\nRandom permutation indices (for array of size 5):" << std::endl;
    std::cout << "  ";
    for (int i = 4; i >= 0; --i) {
        int swapIdx = Random::UniformInt(0, i);
        std::cout << "swap(" << i << "," << swapIdx << ") ";
    }
    std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Thread Safety Demo                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Random_ThreadSafety()
{
    std::cout << "\n=== Thread Safety Notes ===\n\n";

    std::cout << "The Random class uses thread_local Mersenne Twister generator." << std::endl;
    std::cout << "This provides:" << std::endl;
    std::cout << "  ✓ Thread safety (each thread has its own RNG)" << std::endl;
    std::cout << "  ✓ High performance (no locking required)" << std::endl;
    std::cout << "  ✓ Good randomness quality (MT19937)" << std::endl;
    std::cout << std::endl;
    std::cout << "Note: Results are non-reproducible by design (uses random_device for seeding)." << std::endl;
    std::cout << "For reproducible results, you would need to implement seeding control." << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Main Entry Point                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Random()
{
    std::cout << "╔════════════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║               Random Class Documentation Examples                  ║" << std::endl;
    std::cout << "║                  (from docs/base/Random.md)                        ║" << std::endl;
    std::cout << "╚════════════════════════════════════════════════════════════════════╝" << std::endl;

    Demo_Random_UniformReal();
    Demo_Random_UniformInt();
    Demo_Random_Direction2D();
    Demo_Random_Direction3D();
    Demo_Random_UseCases();
    Demo_Random_ThreadSafety();

    std::cout << "\n=== All Random Examples Complete ===\n" << std::endl;
}
