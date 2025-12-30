#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Function.h"
#include "base/Vector.h"

#include "algorithms/RootFinding.h"
#endif

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
// Demo: BracketRoot - Find Initial Bracket
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_BracketRoot()
{
    std::cout << "--- BracketRoot Demo ---\n";
    
    // f(x) = x² - 2, root at √2 ≈ 1.414
    RealFunction f([](Real x) { return x * x - 2; });
    
    Real x1 = 0.5, x2 = 1.0;
    std::cout << "Function: f(x) = x² - 2\n";
    std::cout << "Initial guess: [" << x1 << ", " << x2 << "]\n";
    
    bool found = RootFinding::BracketRoot(f, x1, x2);
    
    if (found) {
        std::cout << "Bracket found: [" << x1 << ", " << x2 << "]\n";
        std::cout << "f(" << x1 << ") = " << f(x1) << "\n";
        std::cout << "f(" << x2 << ") = " << f(x2) << "\n";
    } else {
        std::cout << "No bracket found!\n";
    }
    std::cout << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
// Demo: FindRootBrackets - Find All Root Brackets
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_FindRootBrackets()
{
    std::cout << "--- FindRootBrackets Demo ---\n";
    
    // f(x) = sin(x), roots at 0, π, 2π, 3π, ...
    RealFunction f([](Real x) { return std::sin(x); });
    
    Vector<Real> xb1, xb2;
    int numRoots = RootFinding::FindRootBrackets(f, 0.0, 10.0, 100, xb1, xb2);
    
    std::cout << "Function: f(x) = sin(x)\n";
    std::cout << "Search interval: [0, 10]\n";
    std::cout << "Number of brackets found: " << numRoots << "\n\n";
    
    for (int i = 0; i < numRoots; ++i) {
        std::cout << "Bracket " << i + 1 << ": [" << std::setprecision(4) << xb1[i] 
                  << ", " << xb2[i] << "]\n";
    }
    std::cout << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
// Demo: FindRootBisection - Guaranteed Convergence
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_FindRootBisection()
{
    std::cout << "--- FindRootBisection Demo ---\n";
    
    // f(x) = x³ - x - 2, root near x = 1.52
    RealFunction f([](Real x) { return x*x*x - x - 2; });
    
    Real root = RootFinding::FindRootBisection(f, 1.0, 2.0, 1e-10);
    
    std::cout << "Function: f(x) = x³ - x - 2\n";
    std::cout << "Bracket: [1.0, 2.0]\n";
    std::cout << "Root found: " << std::setprecision(12) << root << "\n";
    std::cout << "Verification: f(root) = " << std::scientific << f(root) << "\n";
    std::cout << std::fixed << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
// Demo: FindRootFalsePosition - Faster Than Bisection
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_FindRootFalsePosition()
{
    std::cout << "--- FindRootFalsePosition Demo ---\n";
    
    // f(x) = e^x - 3x, roots near 0.62 and 1.51
    RealFunction f([](Real x) { return std::exp(x) - 3*x; });
    
    Real root = RootFinding::FindRootFalsePosition(f, 0.0, 1.0, 1e-10);
    
    std::cout << "Function: f(x) = eˣ - 3x\n";
    std::cout << "Bracket: [0.0, 1.0]\n";
    std::cout << "Root found: " << std::setprecision(12) << root << "\n";
    std::cout << "Verification: f(root) = " << std::scientific << f(root) << "\n";
    std::cout << std::fixed << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
// Demo: FindRootSecant - No Derivative Needed
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_FindRootSecant()
{
    std::cout << "--- FindRootSecant Demo ---\n";
    
    // f(x) = cos(x) - x, root near 0.739
    RealFunction f([](Real x) { return std::cos(x) - x; });
    
    Real root = RootFinding::FindRootSecant(f, 0.0, 1.0, 1e-12);
    
    std::cout << "Function: f(x) = cos(x) - x\n";
    std::cout << "Initial guesses: 0.0, 1.0\n";
    std::cout << "Root found: " << std::setprecision(14) << root << "\n";
    std::cout << "Verification: f(root) = " << std::scientific << f(root) << "\n";
    std::cout << std::fixed << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
// Demo: FindRootNewton - Fast Quadratic Convergence
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_FindRootNewton()
{
    std::cout << "--- FindRootNewton Demo ---\n";
    
    // f(x) = x² - 2, root at √2
    RealFunction f([](Real x) { return x * x - 2; });
    
    Real root = RootFinding::FindRootNewton(f, 1.0, 2.0, 1e-14);
    
    std::cout << "Function: f(x) = x² - 2\n";
    std::cout << "Bracket: [1.0, 2.0]\n";
    std::cout << "Root found: " << std::setprecision(15) << root << "\n";
    std::cout << "Expected:   " << std::sqrt(2.0) << "\n";
    std::cout << "Error:      " << std::scientific << std::abs(root - std::sqrt(2.0)) << "\n";
    std::cout << std::fixed << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
// Demo: FindRootRidders - Quadratic Without Derivatives
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_FindRootRidders()
{
    std::cout << "--- FindRootRidders Demo ---\n";
    
    // f(x) = tan(x) - x, root near 4.49
    RealFunction f([](Real x) { return std::tan(x) - x; });
    
    Real root = RootFinding::FindRootRidders(f, 4.0, 4.6, 1e-12);
    
    std::cout << "Function: f(x) = tan(x) - x\n";
    std::cout << "Bracket: [4.0, 4.6]\n";
    std::cout << "Root found: " << std::setprecision(14) << root << "\n";
    std::cout << "Verification: f(root) = " << std::scientific << f(root) << "\n";
    std::cout << std::fixed << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
// Demo: FindRootBrent - Industry Standard
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_FindRootBrent()
{
    std::cout << "--- FindRootBrent Demo (Recommended Default) ---\n";
    
    // f(x) = x⁵ - 3x⁴ + 2x³ + x² - 5x + 2, root near 0.537
    RealFunction f([](Real x) { 
        return x*x*x*x*x - 3*x*x*x*x + 2*x*x*x + x*x - 5*x + 2; 
    });
    
    Real root = RootFinding::FindRootBrent(f, 0.0, 1.0, 1e-14);
    
    std::cout << "Function: f(x) = x⁵ - 3x⁴ + 2x³ + x² - 5x + 2\n";
    std::cout << "Bracket: [0.0, 1.0]\n";
    std::cout << "Root found: " << std::setprecision(15) << root << "\n";
    std::cout << "Verification: f(root) = " << std::scientific << f(root) << "\n";
    std::cout << std::fixed << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
// Demo: Method Comparison
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_RootFinding_Comparison()
{
    std::cout << "--- Method Comparison ---\n";
    
    // Compare all methods on f(x) = x² - 2
    RealFunction f([](Real x) { return x * x - 2; });
    Real exact = std::sqrt(2.0);
    
    std::cout << "Function: f(x) = x² - 2\n";
    std::cout << "Exact root: √2 = " << std::setprecision(15) << exact << "\n\n";
    
    std::cout << std::setw(20) << "Method" << std::setw(20) << "Root" 
              << std::setw(15) << "Error" << "\n";
    std::cout << std::string(55, '-') << "\n";
    
    Real root;
    
    root = RootFinding::FindRootBisection(f, 1.0, 2.0, 1e-14);
    std::cout << std::setw(20) << "Bisection" << std::setw(20) << std::setprecision(12) << root
              << std::setw(15) << std::scientific << std::abs(root - exact) << "\n";
    
    root = RootFinding::FindRootFalsePosition(f, 1.0, 2.0, 1e-14);
    std::cout << std::setw(20) << "False Position" << std::setw(20) << std::setprecision(12) << std::fixed << root
              << std::setw(15) << std::scientific << std::abs(root - exact) << "\n";
    
    root = RootFinding::FindRootSecant(f, 1.0, 2.0, 1e-14);
    std::cout << std::setw(20) << "Secant" << std::setw(20) << std::setprecision(12) << std::fixed << root
              << std::setw(15) << std::scientific << std::abs(root - exact) << "\n";
    
    root = RootFinding::FindRootNewton(f, 1.0, 2.0, 1e-14);
    std::cout << std::setw(20) << "Newton" << std::setw(20) << std::setprecision(12) << std::fixed << root
              << std::setw(15) << std::scientific << std::abs(root - exact) << "\n";
    
    root = RootFinding::FindRootRidders(f, 1.0, 2.0, 1e-14);
    std::cout << std::setw(20) << "Ridders" << std::setw(20) << std::setprecision(12) << std::fixed << root
              << std::setw(15) << std::scientific << std::abs(root - exact) << "\n";
    
    root = RootFinding::FindRootBrent(f, 1.0, 2.0, 1e-14);
    std::cout << std::setw(20) << "Brent" << std::setw(20) << std::setprecision(12) << std::fixed << root
              << std::setw(15) << std::scientific << std::abs(root - exact) << "\n";
    
    std::cout << std::fixed << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
// Main Demo Entry Point
///////////////////////////////////////////////////////////////////////////////////////////
void Docs_Demo_Root_finding()
{
    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║              ROOT FINDING - DOCUMENTATION DEMOS                ║\n";
    std::cout << "║                                                                ║\n";
    std::cout << "║  Demonstrates all 8 root finding methods:                      ║\n";
    std::cout << "║    • BracketRoot - Find initial bracket                        ║\n";
    std::cout << "║    • FindRootBrackets - Find all root brackets                 ║\n";
    std::cout << "║    • FindRootBisection - Guaranteed convergence                ║\n";
    std::cout << "║    • FindRootFalsePosition - Faster than bisection             ║\n";
    std::cout << "║    • FindRootSecant - No derivative Newton                     ║\n";
    std::cout << "║    • FindRootNewton - Fast quadratic convergence               ║\n";
    std::cout << "║    • FindRootRidders - Quadratic without derivatives           ║\n";
    std::cout << "║    • FindRootBrent - Industry standard (recommended)           ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";

    Docs_Demo_BracketRoot();
    Docs_Demo_FindRootBrackets();
    Docs_Demo_FindRootBisection();
    Docs_Demo_FindRootFalsePosition();
    Docs_Demo_FindRootSecant();
    Docs_Demo_FindRootNewton();
    Docs_Demo_FindRootRidders();
    Docs_Demo_FindRootBrent();
    Docs_Demo_RootFinding_Comparison();

    std::cout << "========================================\n";
    std::cout << "   ROOT FINDING DEMOS COMPLETE\n";
    std::cout << "========================================\n\n";
}