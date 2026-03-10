///////////////////////////////////////////////////////////////////////////////////////////
/// MML Documentation Demo: Serializer - Data Serialization Utilities
///////////////////////////////////////////////////////////////////////////////////////////

#include "MMLBase.h"
// Note: Full serializer includes depend on available modules
// #include "mml/tools/Serializer.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace MML;

namespace MML::docs_demos::serializer
{
    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              SERIALIZATION OVERVIEW                              ///
    ///////////////////////////////////////////////////////////////////////////////////////
    //
    // The MML Serializer module provides utilities for saving mathematical objects:
    //
    // SerializerBase.h       - Base types (SerializeError, SerializeResult) and headers
    // SerializerFunctions.h  - Real function serialization (SaveRealFunc, etc.)
    // SerializerCurves.h     - Parametric curve serialization
    // SerializerSurfaces.h   - Surface and scalar function serialization
    // SerializerVectors.h    - Vector field serialization
    // SerializerODE.h        - ODE solution serialization
    // SerializerSimulation.h - Particle simulation serialization

    void demo_serialization_concepts()
    {
        std::cout << "\n=== Serialization Concepts ===\n";

        std::cout << "\nMML Serializer supports the following formats:\n";
        std::cout << "  - Text files for human-readable output\n";
        std::cout << "  - CSV for spreadsheet compatibility\n";
        std::cout << "  - Custom binary formats for efficiency\n";

        std::cout << "\nKey classes:\n";
        std::cout << "  - SerializeResult: Success/failure with error message\n";
        std::cout << "  - SerializeError: Exception type for serialization failures\n";
    }

    void demo_function_serialization()
    {
        std::cout << "\n=== Function Serialization (Example Pattern) ===\n";

        // Example: Save function values to string stream
        std::ostringstream oss;

        auto f = [](double x) { return std::sin(x); };
        double xMin = 0.0, xMax = Constants::PI;
        int numPoints = 11;

        // Write header
        oss << "# Function: sin(x)\n";
        oss << "# Range: [" << xMin << ", " << xMax << "]\n";
        oss << "# Points: " << numPoints << "\n";
        oss << "x,f(x)\n";

        // Write data points
        for (int i = 0; i < numPoints; ++i) {
            double x = xMin + (xMax - xMin) * i / (numPoints - 1);
            double fx = f(x);
            oss << std::fixed << std::setprecision(6) << x << "," << fx << "\n";
        }

        std::cout << "Generated CSV output:\n";
        std::cout << oss.str().substr(0, 300);  // Print first 300 chars
        if (oss.str().size() > 300) std::cout << "...\n";
    }

    void demo_curve_serialization()
    {
        std::cout << "\n=== Parametric Curve Serialization (Example Pattern) ===\n";

        // Example: Serialize a parametric circle
        std::ostringstream oss;

        double tMin = 0.0, tMax = 2 * Constants::PI;
        int numPoints = 17;

        oss << "# Parametric Curve: Circle\n";
        oss << "t,x,y\n";

        for (int i = 0; i < numPoints; ++i) {
            double t = tMin + (tMax - tMin) * i / (numPoints - 1);
            double x = std::cos(t);
            double y = std::sin(t);
            oss << std::fixed << std::setprecision(6) 
                << t << "," << x << "," << y << "\n";
        }

        std::cout << "Parametric curve output (first 300 chars):\n";
        std::cout << oss.str().substr(0, 300);
        if (oss.str().size() > 300) std::cout << "...\n";
    }

    void demo_ode_solution_serialization()
    {
        std::cout << "\n=== ODE Solution Serialization (Example Pattern) ===\n";

        // Example: Simple exponential decay y' = -y
        std::ostringstream oss;

        double t0 = 0.0, tMax = 5.0;
        double y0 = 1.0;
        int numPoints = 11;

        oss << "# ODE Solution: y' = -y\n";
        oss << "# Initial condition: y(0) = " << y0 << "\n";
        oss << "t,y_numerical,y_exact\n";

        for (int i = 0; i < numPoints; ++i) {
            double t = t0 + (tMax - t0) * i / (numPoints - 1);
            double yExact = y0 * std::exp(-t);
            double yNum = yExact * (1.0 + 0.001 * (rand() % 100 - 50) / 50.0);  // Add noise
            oss << std::fixed << std::setprecision(6) 
                << t << "," << yNum << "," << yExact << "\n";
        }

        std::cout << "ODE solution output:\n";
        std::cout << oss.str();
    }

    void demo_vector_field_serialization()
    {
        std::cout << "\n=== Vector Field Serialization (Example Pattern) ===\n";

        // Example: Serialize a simple 2D vector field
        std::ostringstream oss;

        int gridSize = 5;
        double range = 2.0;

        oss << "# Vector Field: F(x,y) = (-y, x)\n";
        oss << "x,y,Fx,Fy,magnitude\n";

        for (int i = 0; i < gridSize; ++i) {
            for (int j = 0; j < gridSize; ++j) {
                double x = -range + 2 * range * i / (gridSize - 1);
                double y = -range + 2 * range * j / (gridSize - 1);
                double Fx = -y;
                double Fy = x;
                double mag = std::sqrt(Fx*Fx + Fy*Fy);
                oss << std::fixed << std::setprecision(3)
                    << x << "," << y << "," << Fx << "," << Fy << "," << mag << "\n";
            }
        }

        std::cout << "Vector field output (first 400 chars):\n";
        std::cout << oss.str().substr(0, 400);
        if (oss.str().size() > 400) std::cout << "...\n";
    }

    void demo_file_io_pattern()
    {
        std::cout << "\n=== File I/O Pattern ===\n";

        std::cout << R"(
// Typical file output pattern:
std::ofstream file("output.csv");
if (!file) {
    throw SerializeError("Cannot open file for writing");
}

// Write header
file << "# Data description\n";
file << "column1,column2,column3\n";

// Write data
for (const auto& row : data) {
    file << row.col1 << "," << row.col2 << "," << row.col3 << "\n";
}

file.close();

// Check for write errors
if (!file) {
    throw SerializeError("Error writing to file");
}
)";
    }

    void demo_result_handling()
    {
        std::cout << "\n=== SerializeResult Pattern ===\n";

        // Demonstrate result pattern (conceptual)
        struct SerializeResult {
            bool success;
            std::string errorMessage;
        };

        auto serialize = [](bool shouldFail) -> SerializeResult {
            if (shouldFail) {
                return {false, "Simulated failure: disk full"};
            }
            return {true, ""};
        };

        // Success case
        auto result1 = serialize(false);
        if (result1.success) {
            std::cout << "Serialization succeeded\n";
        }

        // Failure case
        auto result2 = serialize(true);
        if (!result2.success) {
            std::cout << "Serialization failed: " << result2.errorMessage << "\n";
        }
    }

    void Run()
    {
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "       MML Documentation Demo: Serializer\n";
        std::cout << std::string(70, '=') << "\n";

        demo_serialization_concepts();
        demo_function_serialization();
        demo_curve_serialization();
        demo_ode_solution_serialization();
        demo_vector_field_serialization();
        demo_file_io_pattern();
        demo_result_handling();

        std::cout << "\n" << std::string(70, '=') << "\n";
    }
} // namespace

void Docs_Demo_Serializer() {
    MML::docs_demos::serializer::Run();
}
