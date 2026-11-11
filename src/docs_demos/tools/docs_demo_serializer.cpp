///////////////////////////////////////////////////////////////////////////////////////////
/// MML Documentation Demo: Serializer - Data Serialization Utilities
///////////////////////////////////////////////////////////////////////////////////////////

#include <mml/MMLBase.h>
#include <mml/base/Function.h>
#include <mml/tools/IOResult.h>
#include <mml/tools/Persistence.h>
#include <mml/tools/serializer/SerializerCurves.h>
#include <mml/tools/serializer/SerializerFunctions.h>
#include <mml/tools/serializer/SerializerSimulation.h>
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

void Docs_Demo_BinarySchema()
{
	const Vector<double> original{1, -2};
	std::ostringstream encoded(std::ios::out | std::ios::binary);
	const auto saved = Persistence::SaveBinary(encoded, original);
	if (!saved.success) throw std::runtime_error(saved.message);
	std::istringstream valid(encoded.str(), std::ios::in | std::ios::binary);
	Vector<double> restored;
	const auto loaded = Persistence::LoadBinary(valid, restored);
	if (!loaded.success || restored != original)
		throw std::runtime_error("Schema v1 binary round trip failed");
	for (bool version : {false, true}) {
		auto bytes = encoded.str();
		if (version) bytes[16] = 2;
		else bytes[10] = 49;
		std::istringstream input(bytes, std::ios::in | std::ios::binary);
		Vector<double> sentinel{9};
		const auto result = Persistence::LoadBinary(input, sentinel);
		const auto expected = version ? IOError::UNSUPPORTED_VERSION : IOError::SCHEMA_MISMATCH;
		if (result.success || result.error != expected || sentinel != Vector<double>{9} || input.tellg() != 48)
			throw std::runtime_error("Invalid binary schema/header was not rejected before payload decoding");
		if (version) {
			input.seekg(0);
			const auto envelope = Persistence::ReadBinaryEnvelopeHeader(input);
			if (!envelope.result.success || envelope.header.object_schema_version != 2)
				throw std::runtime_error("Envelope inspection lost future object schema support");
		}
		std::cout << IOErrorName(result.error) << ": " << result.message << '\n';
	}
	std::cout << "Schema v1 round trip verified; rejected headers leave the destination unchanged.\n";
}

void Docs_Demo_SerializerSampling()
{
	const Real largest = std::numeric_limits<Real>::max();
	std::vector<Real> coordinates;
	RealFunctionFromStdFunc constant([&](Real x) { coordinates.push_back(x); return Real(1); });
	for (bool multi : {false, true}) {
		coordinates.clear();
		std::ostringstream out;
		const auto result = multi
			? Serializer::SaveRealMultiFunc(out, std::vector<IRealFunction*>{&constant}, "constant", {"one"},
				0, largest, 3, std::numeric_limits<Real>::max_digits10)
			: Serializer::SaveRealFunc(out, constant, "constant", -largest, largest, 3,
				std::numeric_limits<Real>::max_digits10);
		const std::vector<Real> expected = multi
			? std::vector<Real>{0, largest / 2, largest}
			: std::vector<Real>{-largest, 0, largest};
		if (!result.success || coordinates != expected)
			throw std::runtime_error("Serializer generated incorrect finite sample coordinates");
		std::istringstream input(out.str());
		std::string line;
		for (int i = 0; i < (multi ? 8 : 6); ++i)
			if (!std::getline(input, line)) throw std::runtime_error("Missing sampling demo header");
		for (Real expectedX : expected) {
			Real x = 0, value = 0;
			if (!(input >> x >> value) || x != expectedX || value != 1)
				throw std::runtime_error("Serialized coordinates did not round-trip");
		}
		std::cout << (multi ? "Multi-function [0, max]" : "Function [-max, max]")
			<< ": all 3 coordinates finite, exact endpoints, round-trip output verified.\n";
	}
}

void Docs_Demo_SerializerMetadata()
{
	const auto prefix = std::filesystem::temp_directory_path() /
		("mml_metadata_demo_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
	const std::filesystem::path curvePath = prefix.string() + "_curve.mml";
	const std::filesystem::path particlesPath = prefix.string() + "_particles.mml";
	struct Cleanup {
		std::vector<std::filesystem::path> paths;
		~Cleanup() {
			for (const auto& path : paths) {
				std::error_code error;
				std::filesystem::remove(path, error);
				if (error) std::cerr << "Cannot remove metadata demo file: " << error.message() << '\n';
			}
		}
	} cleanup{{curvePath, particlesPath}};
	auto read = [](const std::filesystem::path& path) {
		std::ifstream in(path);
		if (!in.is_open()) throw std::runtime_error("Cannot open metadata demo file");
		std::ostringstream out;
		out << in.rdbuf();
		if (in.bad()) throw std::runtime_error("Cannot read metadata demo file");
		return out.str();
	};
	const auto curve = Serializer::SaveAsParamCurve2D(Vector<Real>{0, 1}, Vector<Real>{0, 1}, "", curvePath.string());
	if (!curve.success || read(curvePath).find(
		"MML_PARAMETRIC_CURVE_CARTESIAN_2D\nVERSION: 1\n\nt1: 0\nt2: 1\nNumPoints: 2\n") != 0)
		throw std::runtime_error("Empty curve title shifted the mandatory header positions");
	std::cout << "Empty curve title occupies line 3; t1 remains on line 4.\n";
	for (int dimension : {2, 3}) {
		{
			std::ofstream out(particlesPath, std::ios::binary);
			out << "SENTINEL";
			out.close();
			if (out.fail()) throw std::runtime_error("Cannot seed metadata demo file");
		}
		auto save = [&](const std::string& color) {
			if (dimension == 2)
				return Serializer::SaveParticleSimulation2D(particlesPath.string(), 1, 10, 20, {{{0, 0}}}, {color}, {1}, 1);
			return Serializer::SaveParticleSimulation3D(particlesPath.string(), 1, 10, 20, 30, {{{0, 0, 0}}}, {color}, {1}, 1);
		};
		const auto invalid = save("red\nNumSteps: 999");
		if (invalid.error != IOError::INVALID_PARAMETERS || invalid.filename != particlesPath.string() ||
			read(particlesPath) != "SENTINEL")
			throw std::runtime_error("Malformed particle catalog was not rejected before opening the destination");
		const auto valid = save("red");
		if (!valid.success || read(particlesPath).find("Ball_1 red 1\nNumSteps: 1\nStep 0 0\n") == std::string::npos)
			throw std::runtime_error("Valid particle catalog no longer matches the reader token layout");
		std::cout << dimension << "D particle catalog rejected a multiline color without modifying SENTINEL; valid red token preserved.\n";
	}
}

void Docs_Demo_SerializerPreflight()
{
	const auto path = std::filesystem::temp_directory_path() /
		("mml_preflight_demo_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".mml");
	struct Cleanup {
		std::filesystem::path path;
		~Cleanup() {
			std::error_code error;
			std::filesystem::remove(path, error);
			if (error) std::cerr << "Cannot remove preflight demo file: " << error.message() << '\n';
		}
	} cleanup{path};
	{
		std::ofstream out(path, std::ios::binary);
		out << "SENTINEL";
		out.close();
		if (out.fail()) throw std::runtime_error("Cannot seed preflight demo file");
	}
	auto read = [&]() {
		std::ifstream in(path, std::ios::binary);
		if (!in.is_open()) throw std::runtime_error("Cannot open preflight demo file");
		std::ostringstream out;
		out << in.rdbuf();
		if (in.bad()) throw std::runtime_error("Cannot read preflight demo file");
		return out.str();
	};
	int evaluations = 0;
	RealFunctionFromStdFunc square([&](Real x) { ++evaluations; return x * x; });
	const std::string filename = path.string();
	const auto invalidCount = Serializer::SaveRealFunc(square, "square", 0, 1, 1, filename);
	if (invalidCount.error != IOError::INVALID_PARAMETERS || invalidCount.filename != filename || read() != "SENTINEL")
		throw std::runtime_error("Invalid sampling request modified the destination");
	const std::string invalidName = filename + std::string(1, '\0') + "_ignored_suffix.mml";
	const auto rejectedName = Serializer::SaveRealFuncEquallySpaced(square, "square", 0, 1, 3, invalidName);
	if (rejectedName.error != IOError::INVALID_PARAMETERS || rejectedName.filename != invalidName ||
		read() != "SENTINEL" || evaluations != 0)
		throw std::runtime_error("Embedded NUL filename was not rejected before opening the destination");
	std::cout << "Invalid sample count: " << invalidCount.message << '\n';
	std::cout << "Embedded NUL filename: " << rejectedName.message << '\n';
	std::cout << "Both requests preserved all 8 SENTINEL bytes without evaluating the function.\n";
	const auto valid = Serializer::SaveRealFunc(square, "square", 0, 1, 3, filename);
	if (!valid.success || evaluations != 3 || read().find("MML_REAL_FUNCTION") != 0)
		throw std::runtime_error("Valid serializer request failed");
}

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

void Docs_Demo_JsonConversions()
{
	Vector<float> original{1, 2};
	std::ostringstream saved;
	if (!Persistence::SaveJson(saved, original).success)
		throw std::runtime_error("Could not encode JSON conversion example");
	auto parsed = Persistence::ParseJson(saved.str());
	if (!parsed.result.success)
		throw std::runtime_error(parsed.result.message);
	parsed.value.object_value.at("shape").array_value[0] = Persistence::JsonValue::Number(2.9);
	std::stringstream malformed;
	if (!Persistence::WriteJson(malformed, parsed.value).success)
		throw std::runtime_error("Could not encode malformed-dimension example");
	Vector<float> destination{9};
	const auto rejected = Persistence::LoadJson(malformed, destination);
	if (rejected.error != IOError::SCHEMA_MISMATCH || destination != Vector<float>{9})
		throw std::runtime_error("Fractional JSON dimension was not rejected transactionally");
	std::cout << "Fractional JSON dimension: " << rejected.message << '\n';

	parsed.value.object_value.at("shape").array_value[0] = Persistence::JsonValue::Number(2);
	parsed.value.object_value.at("data").array_value[1] = Persistence::JsonValue::Number(1e100);
	std::stringstream overflow;
	if (!Persistence::WriteJson(overflow, parsed.value).success)
		throw std::runtime_error("Could not encode scalar-range example");
	const auto range = Persistence::LoadJson(overflow, destination);
	if (range.error != IOError::UNSUPPORTED_SCALAR || destination != Vector<float>{9})
		throw std::runtime_error("JSON scalar overflow was not rejected transactionally");
	std::cout << "JSON scalar overflow: " << range.message << '\n';
}

void Docs_Demo_ParticleHistory()
{
	const auto prefix = std::filesystem::temp_directory_path() /
		("mml_particle_demo_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
	const std::string outputPath = prefix.string() + "_2d.mml";
	const std::string outputPath3D = prefix.string() + "_3d.mml";
	struct Cleanup {
		std::string first, second;
		~Cleanup() {
			for (const auto& path : {first, second}) {
				std::error_code error;
				std::filesystem::remove(path, error);
				if (error) std::cerr << "Cannot remove particle demo file: " << error.message() << '\n';
			}
		}
	} cleanup{outputPath, outputPath3D};
	std::vector<std::vector<Pnt2Cart>> trajectory(2);
	for (int particle = 0; particle < 2; ++particle)
		for (int step = 0; step < 5; ++step)
			trajectory[particle].emplace_back(Real(step), Real(particle));
	const auto result = Serializer::SaveParticleSimulation2D(
		outputPath, 2, 10, 20, trajectory, {"red", "blue"}, {Real(0.5), Real(0.75)}, Real(0.25), 2);
	if (!result.success)
		throw std::runtime_error(result.message);
	auto read = [](const std::string& path) {
		std::ifstream in(path);
		if (!in.is_open()) throw std::runtime_error("Cannot read particle demo file");
		std::ostringstream text;
		text << in.rdbuf();
		if (in.bad()) throw std::runtime_error("Cannot read particle demo data");
		return text.str();
	};
	const auto saved = read(outputPath);
	if (saved.find("NumSteps: 3\n") == std::string::npos || saved.find("Step 2 1\n") == std::string::npos)
		throw std::runtime_error("Particle decimation lost its ceiling count or original timestamp");
	trajectory[1].pop_back();
	const auto rejected = Serializer::SaveParticleSimulation2D(
		outputPath, 2, 10, 20, trajectory, {"red", "blue"}, {Real(0.5), Real(0.75)}, Real(0.25), 2);
	if (rejected.error != IOError::INVALID_PARAMETERS || read(outputPath) != saved)
		throw std::runtime_error("Ragged particle history was not rejected before opening the file");

	std::vector<std::vector<Pnt3Cart>> trajectory3D(2);
	for (int particle = 0; particle < 2; ++particle)
		for (int step = 0; step < 5; ++step)
			trajectory3D[particle].emplace_back(Real(step), Real(particle), Real(particle + step));
	const auto sparse = Serializer::SaveParticleSimulation3D(
		outputPath3D, 2, 10, 20, 30, trajectory3D, {"red", "blue"}, {Real(0.5), Real(0.75)}, Real(0.25),
		std::numeric_limits<int>::max());
	if (!sparse.success || read(outputPath3D).find("NumSteps: 1\n") == std::string::npos)
		throw std::runtime_error("Extreme particle stride was not handled correctly");
	std::cout << "Particle export: rectangular histories, original timestamps, pre-open rejection and INT_MAX stride verified.\n";
}

void Docs_Demo_SavePrevalidation()
{
	const auto path = std::filesystem::temp_directory_path() /
		("mml_save_demo_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".mmlj");
	struct Cleanup {
		std::filesystem::path path;
		~Cleanup() {
			std::error_code error;
			std::filesystem::remove(path, error);
			if (error) std::cerr << "Cannot remove save demo file: " << error.message() << '\n';
		}
	} cleanup{path};
	std::ofstream seed(path, std::ios::binary);
	seed << "SENTINEL";
	if (!FinishOutputFile(seed, path.string(), SerializeSuccess()).success)
		throw std::runtime_error("Could not seed save-prevalidation example");
	Vector<Real> invalid{1, std::numeric_limits<Real>::infinity()};
	const auto rejected = Persistence::Save(invalid, path.string());
	if (rejected.error != IOError::UNSUPPORTED_SCALAR || rejected.filename != path.string())
		throw std::runtime_error("Invalid JSON save was not rejected with file context");
	std::ifstream in(path, std::ios::binary);
	std::ostringstream contents;
	contents << in.rdbuf();
	if (contents.str() != "SENTINEL")
		throw std::runtime_error("Invalid JSON save changed the existing file");
	std::cout << "Prevalidated JSON save preserved existing file: " << rejected.message << '\n';
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
	Docs_Demo_JsonConversions();
	Docs_Demo_SavePrevalidation();
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
