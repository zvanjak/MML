///////////////////////////////////////////////////////////////////////////////////////////
/// MML Documentation Demo: Serializer - Data Serialization Utilities
///////////////////////////////////////////////////////////////////////////////////////////

#include <mml/MMLBase.h>
#include <mml/tools/IOResult.h>
#include <mml/tools/Persistence.h>
// Note: Full serializer includes depend on available modules
// #include <mml/tools/Serializer.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <filesystem>
#include <chrono>

using namespace MML;

void Docs_Demo_VectorPersistence()
{
	auto require = [](const IOResult& result) {
		if (!result) throw std::runtime_error(result.message);
	};
	const auto filename = std::filesystem::temp_directory_path() /
		("mml_vector_demo_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".txt");
	struct Cleanup {
		std::filesystem::path path;
		~Cleanup() {
			std::error_code error;
			std::filesystem::remove(path, error);
			if (error) std::cerr << "Cannot remove vector demo file: " << error.message() << '\n';
		}
	} cleanup{filename};
	Vector<Real> values{1.25, -2.5, std::nextafter(Real(1), Real(2))};
	require(Persistence::Save(values, filename.string()));
	Vector<Real> restored;
	require(Persistence::Load(filename.string(), restored));
	if (restored != values)
		throw std::runtime_error("Vector file persistence lost precision");
	std::stringstream text;
	text << std::fixed << std::setprecision(2);
	require(Persistence::Save(text, values, Persistence::Format::Text));
	require(Persistence::Load(text, restored, Persistence::Format::Text));
	if (restored != values || text.precision() != 2)
		throw std::runtime_error("Vector stream persistence lost precision or changed formatting");

	Vector<int> integers{1, -2, 3};
	std::stringstream integerText;
	require(Persistence::Save(integerText, integers, Persistence::Format::Text));
	std::cout << "Count-prefixed integer vector:\n" << integerText.str();
	Vector<Complex> phasors{Complex(1.25, -0.0),
	                       Complex(std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN())};
	std::stringstream complexText;
	require(Persistence::Save(complexText, phasors, Persistence::Format::Text));
	Vector<Complex> decoded;
	require(Persistence::Load(complexText, decoded, Persistence::Format::Text));
	if (!std::signbit(decoded[0].imag()) || !std::isinf(decoded[1].real()) || !std::isnan(decoded[1].imag()))
		throw std::runtime_error("Vector text failed to preserve special-value classification");
	std::cout << "Complex values and canonical special tokens:\n" << complexText.str();

	Vector<Real> unchanged{9, 8};
	std::istringstream malformed("2\n1\nnot-a-number\n");
	const auto failed = Persistence::Load(malformed, unchanged, Persistence::Format::Text);
	if (failed.success || unchanged != Vector<Real>{9, 8} || failed.offset != 3)
		throw std::runtime_error("Failed vector load did not preserve its destination and diagnostics");
	std::cout << "Transactional load failure at line " << failed.offset << ": " << failed.message << '\n';
	Persistence::LoadOptions limited;
	limited.max_allocation_bytes = sizeof(Real);
	std::istringstream oversized("2\n1\n2\n");
	const auto rejected = Persistence::Load(oversized, unchanged, Persistence::Format::Text, limited);
	if (rejected.error != IOError::ALLOCATION_LIMIT_EXCEEDED)
		throw std::runtime_error("Vector payload allocation limit was not enforced");
	std::cout << "Payload budget: " << rejected.message << '\n';
}

void Docs_Demo_PersistencePortability()
{
	auto require = [](const IOResult& result) {
		if (!result.success)
			throw std::runtime_error(result.message);
	};
	Persistence::SaveOptions options;
	options.precision = 0;
	options.pretty_json = true;
	options.include_metadata = true;
	options.allow_lossy = false;
	options.metadata.title = "test fixture";
	options.metadata.tags["source"] = "unit-test";

	Matrix<double> original(1, 2);
	original(0, 0) = std::nextafter(1.0, 2.0);
	original(0, 1) = 1.0 / 7.0;
	std::stringstream text;
	text << std::fixed << std::setprecision(2);
	require(Persistence::Save(text, original, Persistence::Format::Text, options));
	Matrix<double> restored;
	require(Persistence::Load(text, restored, Persistence::Format::Text));
	if (restored(0, 0) != original(0, 0) || restored(0, 1) != original(0, 1))
		throw std::runtime_error("Matrix text round trip lost precision");
	if (text.precision() != 2 || (text.flags() & std::ios::floatfield) != std::ios::fixed)
		throw std::runtime_error("Matrix persistence changed caller formatting");

	Vector<long double> values{1.0L, -0.0L, std::numeric_limits<long double>::max()};
	std::stringstream binary(std::ios::in | std::ios::out | std::ios::binary);
	require(Persistence::SaveBinary(binary, values));
	const auto envelope = Persistence::ReadBinaryEnvelopeHeader(binary);
	require(envelope.result);
	if (envelope.header.scalar_byte_size != Persistence::BinaryScalarEncoding::LONG_DOUBLE_BYTES)
		throw std::runtime_error("Incorrect canonical long-double width");
	binary.seekg(0);
	Vector<long double> decoded;
	require(Persistence::LoadBinary(binary, decoded));
	if (decoded.size() != values.size() || decoded[0] != values[0] ||
	    !std::signbit(decoded[1]) || decoded[2] != values[2])
		throw std::runtime_error("Canonical long-double round trip failed");

	std::cout << "Exact matrix text round trip; caller precision remains " << text.precision() << '\n'
	          << "Canonical long-double bytes: " << envelope.header.scalar_byte_size
	          << "; host storage bytes: " << sizeof(long double) << '\n';
}

namespace MML::docs_demos::serializer
{
	///////////////////////////////////////////////////////////////////////////////////////
	///                              SERIALIZATION OVERVIEW                              ///
	///////////////////////////////////////////////////////////////////////////////////////
	//
	// The MML Serializer module provides utilities for saving mathematical objects:
	//
	// SerializerBase.h       - Base types (IOError, IOResult) and headers
	// SerializerFunctions.h  - Real function serialization (SaveRealFunc, etc.)
	// SerializerCurves.h     - Parametric curve serialization
	// SerializerSurfaces.h   - Surface and scalar function serialization
	// SerializerVectors.h    - Vector field serialization
	// SerializerODE.h        - ODE solution serialization
	// SerializerSimulation.h - Particle simulation serialization

	void demo_serialization_concepts() {
		std::cout << "\n=== Serialization Concepts ===\n";

		std::cout << "\nMML Serializer supports the following formats:\n";
		std::cout << "  - Text files for human-readable output\n";
		std::cout << "  - CSV for spreadsheet compatibility\n";
		std::cout << "  - Custom binary formats for efficiency\n";

		std::cout << "\nKey classes:\n";
		std::cout << "  - IOResult: Success/failure with error message\n";
		std::cout << "  - IOError: Exception type for serialization failures\n";
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
    throw IOError("Cannot open file for writing");
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
    throw IOError("Error writing to file");
}
)";
		}

    void demo_result_handling()
    {
			std::cout << "\n=== IOResult Pattern ===\n";

			auto serialize = [](bool shouldFail) -> IOResult {
				if (shouldFail) {
					return {false, IOError::WRITE_FAILED, "Simulated failure: disk full"};
				}
				return IOResult::Success();
			};

			// Success case
			auto result1 = serialize(false);
			if (result1.success) {
				std::cout << "Serialization succeeded\n";
        }

        // Failure case
        auto result2 = serialize(true);
        if (!result2.success) {
					std::cout << "Serialization failed: " << result2.message << "\n";
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
