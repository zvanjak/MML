///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  MML Precision Testing Application                                                ///
///  ==================================                                               ///
///                                                                                   ///
///  Comprehensive precision testing for all MML numerical algorithms:               ///
///    - Derivation (numerical derivatives)                                           ///
///    - Integration (1D, 2D, 3D quadrature)                                          ///
///    - Interpolation (polynomial, spline)                                           ///
///    - ODE Solvers (Euler, RK4, adaptive)                                           ///
///    - Root Finding (bisection, Newton, Brent)                                      ///
///    - Linear Algebra (matrix solvers)                                              ///
///                                                                                   ///
///  Usage: MML_PrecisionTestingApp [category] [--format fmt] [--output file]         ///
///         Categories: all, derivation, integration, interpolation, ode, roots, linalg
///         Formats: ascii (default), markdown, csv                                   ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#endif

#include "PrecisionTestFramework.h"

#include <iostream>
#include <string>
#include <vector>
#include <cstring>

using namespace MML::PrecisionTesting;

// Forward declarations for test functions
void Test_Precision_Derivation();
void Test_Precision_Integration();
void Test_Precision_Interpolation();
void Test_Precision_ODE();
void Test_Precision_RootFinding();
void Test_Precision_LinearAlgebra();

//--------------------------------------------------------------------------------------------------
// COMMAND LINE PARSING
//--------------------------------------------------------------------------------------------------

struct CommandLineArgs {
    std::string category = "all";
    std::string format = "ascii";
    std::string output_file;
    bool brief = false;
    bool help = false;
};

void printUsage() {
    std::cout << R"(
MML Precision Testing Application
=================================

Usage: MML_PrecisionTestingApp [options] [category]

Categories:
  all           Run all precision tests (default)
  derivation    Test numerical differentiation
  integration   Test numerical integration
  interpolation Test interpolation methods
  ode           Test ODE solvers
  roots         Test root finding algorithms
  linalg        Test linear algebra solvers

Options:
  --format <fmt>    Output format: ascii (default), markdown, csv
  --output <file>   Write results to file
  --brief           Show summary only, skip detailed tables
  --help, -h        Show this help message

Examples:
  MML_PrecisionTestingApp                     # Run all tests
  MML_PrecisionTestingApp derivation          # Test derivation only
  MML_PrecisionTestingApp --format markdown   # Output as markdown
  MML_PrecisionTestingApp --output report.md --format markdown

)";
}

CommandLineArgs parseArgs(int argc, char* argv[]) {
    CommandLineArgs args;
    
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        
        if (arg == "--help" || arg == "-h") {
            args.help = true;
        }
        else if (arg == "--brief") {
            args.brief = true;
        }
        else if (arg == "--format" && i + 1 < argc) {
            args.format = argv[++i];
        }
        else if (arg == "--output" && i + 1 < argc) {
            args.output_file = argv[++i];
        }
        else if (arg[0] != '-') {
            args.category = arg;
        }
    }
    
    return args;
}

//--------------------------------------------------------------------------------------------------
// MAIN
//--------------------------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    CommandLineArgs args = parseArgs(argc, argv);
    
    if (args.help) {
        printUsage();
        return 0;
    }
    
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                   MML PRECISION TESTING APPLICATION                          ║\n";
    std::cout << "║                   ================================                           ║\n";
    std::cout << "║                                                                              ║\n";
    std::cout << "║   Comprehensive precision verification for numerical algorithms             ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
    
    std::cout << "Configuration:\n";
    std::cout << "  Category: " << args.category << "\n";
    std::cout << "  Format:   " << args.format << "\n";
    if (!args.output_file.empty()) {
        std::cout << "  Output:   " << args.output_file << "\n";
    }
    std::cout << "\n";
    
    // Run requested tests
    bool run_all = (args.category == "all");
    
    if (run_all || args.category == "derivation") {
        std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
        std::cout << "  DERIVATION PRECISION TESTS\n";
        std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
        Test_Precision_Derivation();
    }
    
    if (run_all || args.category == "integration") {
        std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
        std::cout << "  INTEGRATION PRECISION TESTS\n";
        std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
        Test_Precision_Integration();
    }
    
    if (run_all || args.category == "interpolation") {
        std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
        std::cout << "  INTERPOLATION PRECISION TESTS\n";
        std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
        Test_Precision_Interpolation();
    }
    
    if (run_all || args.category == "ode") {
        std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
        std::cout << "  ODE SOLVER PRECISION TESTS\n";
        std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
        Test_Precision_ODE();
    }
    
    if (run_all || args.category == "roots") {
        std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
        std::cout << "  ROOT FINDING PRECISION TESTS\n";
        std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
        Test_Precision_RootFinding();
    }
    
    if (run_all || args.category == "linalg") {
        std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
        std::cout << "  LINEAR ALGEBRA PRECISION TESTS\n";
        std::cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n";
        Test_Precision_LinearAlgebra();
    }
    
    std::cout << "╔══════════════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                     PRECISION TESTING COMPLETE                               ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════════════════╝\n\n";
    
    return 0;
}
