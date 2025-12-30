///////////////////////////////////////////////////////////////////////////////////////////
///  MML Precision Testing - Framework Verification Test
///  
///  Quick test to verify PrecisionTestFramework compiles and works correctly
///////////////////////////////////////////////////////////////////////////////////////////

#include "PrecisionTestFramework.h"

#include <iostream>
#include <cmath>

using namespace MML::PrecisionTesting;

// Simple test: numerical derivative of sin(x) at x=1
void testDerivativeFramework()
{
    PrecisionTestSuite suite("Derivative Framework Test", "Verify framework with simple sin(x) derivative");
    
    double x = 1.0;
    double exact = std::cos(x);  // d/dx sin(x) = cos(x)
    
    // Test different step sizes
    std::vector<double> h_values = {0.1, 0.01, 0.001, 0.0001, 0.00001};
    
    for (double h : h_values) {
        // Central difference: (f(x+h) - f(x-h)) / (2h)
        double computed = (std::sin(x + h) - std::sin(x - h)) / (2 * h);
        
        PrecisionTestResult result("CentralDiff", "sin(x)", exact, computed);
        result.parameters = "h=" + std::to_string(h);
        suite.addResult(result);
    }
    
    // Print results
    suite.printHeader();
    suite.printDetailedTable();
    suite.printSummary();
}

// Test error order matrix with multiple algorithms
void testMultiAlgorithmMatrix()
{
    PrecisionTestSuite suite("Multi-Algorithm Test", "Compare forward, backward, and central differences");
    
    // Test functions: sin at different points
    std::vector<std::pair<std::string, double>> test_points = {
        {"x=0.5", 0.5},
        {"x=1.0", 1.0},
        {"x=2.0", 2.0}
    };
    
    double h = 0.0001;
    
    for (const auto& [name, x] : test_points) {
        double exact = std::cos(x);
        
        // Forward difference: (f(x+h) - f(x)) / h
        double forward = (std::sin(x + h) - std::sin(x)) / h;
        suite.addResult("Forward", name, exact, forward);
        
        // Backward difference: (f(x) - f(x-h)) / h
        double backward = (std::sin(x) - std::sin(x - h)) / h;
        suite.addResult("Backward", name, exact, backward);
        
        // Central difference: (f(x+h) - f(x-h)) / (2h)
        double central = (std::sin(x + h) - std::sin(x - h)) / (2 * h);
        suite.addResult("Central", name, exact, central);
    }
    
    suite.printHeader();
    suite.printErrorOrderMatrix();
    suite.printSummary();
    
    // Test markdown generation
    std::cout << "\n--- MARKDOWN OUTPUT ---\n";
    std::cout << suite.toMarkdownMatrix();
}

// Test CSV export
void testExport()
{
    PrecisionTestSuite suite("Export Test");
    
    suite.addResult("Algo1", "Func1", 1.0, 1.0001);
    suite.addResult("Algo1", "Func2", 2.0, 2.00001);
    suite.addResult("Algo2", "Func1", 1.0, 0.999999);
    
    if (suite.exportCSV("results/precision_test_export.csv")) {
        std::cout << "\n✓ CSV exported to results/precision_test_export.csv\n";
    }
}

// Test precision rating
void testRating()
{
    std::cout << "\n=== PRECISION RATING TEST ===\n";
    
    std::vector<double> errors = {1e-14, 1e-10, 1e-6, 1e-3, 0.1};
    
    for (double err : errors) {
        PrecisionRating rating = ratePrecision(err);
        std::cout << "Error " << std::scientific << err 
                  << " -> " << ratingToEmoji(rating) << " " << ratingToString(rating) << "\n";
    }
}

int main()
{
    std::cout << "======================================================================\n";
    std::cout << "  MML PRECISION TEST FRAMEWORK - VERIFICATION\n";
    std::cout << "======================================================================\n";
    
    testDerivativeFramework();
    testMultiAlgorithmMatrix();
    testRating();
    testExport();
    
    std::cout << "\n======================================================================\n";
    std::cout << "  FRAMEWORK VERIFICATION COMPLETE!\n";
    std::cout << "======================================================================\n";
    
    return 0;
}
