///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        serializer_tests.cpp                                                ///
///  Description: Comprehensive unit tests for Serializer                             ///
///                                                                                   ///
///  This test suite validates all serialization formats to catch format changes     ///
///  and regressions. Tests cover:                                                    ///
///    - Real function serialization (equally spaced, specified points)              ///
///    - Multi-function serialization (multiple interpolation types)                 ///
///    - Parametric curve serialization (2D, 3D)                                      ///
///    - Parametric surface serialization                                             ///
///    - Scalar function serialization (2D, 3D grids)                                ///
///    - Vector function serialization (2D, 3D Cartesian and spherical)              ///
///    - ODE solution serialization (component, multi-func, parametric)              ///
///    - Particle simulation serialization (2D, 3D)                                   ///
///    - Error handling for invalid parameters                                        ///
///    - Golden file comparison for content validation                                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <complex>
#include <filesystem>
#include <functional>
#include <limits>
#include <locale>
#include <regex>
#include <string>
#include <type_traits>
#include <vector>

#include <mml/MMLBase.h>
#include <mml/tools/Serializer.h>
#include <mml/tools/serializer/SerializerFieldLines.h>
#include <mml/tools/Persistence.h>
#include <mml/interfaces/IFunction.h>
#include <mml/base/InterpolatedFunction.h>

using namespace MML;

// Alias for cleaner error type access (now at MML namespace level, not nested in class)
using MML::IOError;

namespace MML::Tests::Tools::SerializerTests {
	namespace {
		const char* RealTypeName() {
			if constexpr (std::is_same_v<Real, float>) return "float";
			if constexpr (std::is_same_v<Real, long double>) return "long_double";
			return "double";
		}

		void ReplaceAll(std::string& text, const std::string& from, const std::string& to) {
			std::size_t position = 0;
			while ((position = text.find(from, position)) != std::string::npos) {
				text.replace(position, from.size(), to);
				position += to.size();
			}
		}

		std::string WithRealMetadata(std::string json) {
			ReplaceAll(json, "\"scalar_bytes\": 8", "\"scalar_bytes\": " + std::to_string(sizeof(Real)));
			ReplaceAll(json, "\"scalar_bytes\":8", "\"scalar_bytes\":" + std::to_string(sizeof(Real)));
			ReplaceAll(json, "\"real_type\": \"double\"", "\"real_type\": \"" + std::string(RealTypeName()) + "\"");
			ReplaceAll(json, "\"real_type\":\"double\"", "\"real_type\":\"" + std::string(RealTypeName()) + "\"");
			return json;
		}
	}

	TEST_CASE("Persistence options expose conservative defaults and metadata", "[Persistence][Options]") {
		Persistence::SaveOptions saveOptions;
		REQUIRE(saveOptions.precision == 0);
		REQUIRE(saveOptions.include_metadata);
		REQUIRE(saveOptions.pretty_json);
		REQUIRE_FALSE(saveOptions.allow_lossy);
		REQUIRE(saveOptions.metadata.title.empty());
		REQUIRE(saveOptions.metadata.description.empty());
		REQUIRE(saveOptions.metadata.units.empty());
		REQUIRE(saveOptions.metadata.coordinate_system.empty());
		REQUIRE(saveOptions.metadata.tags.empty());

		saveOptions.metadata.title = "matrix fixture";
		saveOptions.metadata.units = "m";
		saveOptions.metadata.coordinate_system = "cartesian";
		saveOptions.metadata.tags["source"] = "unit-test";

		REQUIRE(saveOptions.metadata.title == "matrix fixture");
		REQUIRE(saveOptions.metadata.units == "m");
		REQUIRE(saveOptions.metadata.coordinate_system == "cartesian");
		REQUIRE(saveOptions.metadata.tags.at("source") == "unit-test");

		Persistence::LoadOptions loadOptions;
		REQUIRE(loadOptions.strict_schema);
		REQUIRE_FALSE(loadOptions.allow_scalar_conversion);
		REQUIRE(loadOptions.max_allocation_bytes == (1ull << 30));

		Persistence::SplineInterpolationOptions splineOptions;
		REQUIRE(splineOptions.first_derivative == SplineInterpRealFunc::NaturalBoundaryDerivative);
		REQUIRE(splineOptions.last_derivative == SplineInterpRealFunc::NaturalBoundaryDerivative);
	}

	TEST_CASE("Persistence error taxonomy exposes stable names and result helpers", "[Persistence][Errors]") {
		REQUIRE(std::string(IOErrorName(IOError::OK)) == "OK");
		REQUIRE(std::string(IOErrorName(IOError::UNSUPPORTED_VERSION)) == "UNSUPPORTED_VERSION");
		REQUIRE(std::string(IOErrorName(IOError::TYPE_MISMATCH)) == "TYPE_MISMATCH");
		REQUIRE(std::string(IOErrorName(IOError::ENDIAN_MISMATCH)) == "ENDIAN_MISMATCH");
		REQUIRE(std::string(IOErrorName(IOError::ALLOCATION_LIMIT_EXCEEDED)) == "ALLOCATION_LIMIT_EXCEEDED");
		REQUIRE(std::string(IOErrorName(IOError::SCHEMA_MISMATCH)) == "SCHEMA_MISMATCH");
		REQUIRE(std::string(IOErrorName(IOError::MALFORMED_INPUT)) == "MALFORMED_INPUT");
		REQUIRE(std::string(IOErrorName(IOError::UNSUPPORTED_SCALAR)) == "UNSUPPORTED_SCALAR");
		REQUIRE(std::string(IOErrorName(IOError::TRUNCATED_INPUT)) == "TRUNCATED_INPUT");

		IOResult ok = SerializeSuccess();
		REQUIRE(ok.success);
		REQUIRE(ok.error == IOError::OK);
		REQUIRE(ok.message == "Success");

		IOResult failure = SerializeFailure(IOError::TRUNCATED_INPUT, "missing binary payload");
		REQUIRE_FALSE(failure.success);
		REQUIRE(failure.error == IOError::TRUNCATED_INPUT);
		REQUIRE(failure.message.find("TRUNCATED_INPUT") != std::string::npos);
		REQUIRE(failure.message.find("missing binary payload") != std::string::npos);
	}

	TEST_CASE("Persistence JSON parser reads MML schema-shaped objects", "[Persistence][JSON]") {
		const std::string json = R"json({
			"mml": {
				"format": "MML_JSON",
				"version": 1,
				"object": "Matrix",
				"object_version": 1,
				"scalar": "Real"
			},
			"shape": [2, 3],
			"data": [1.0, -2.5, 3.25],
			"metadata": {"title": "small matrix"}
		})json";

		Persistence::JsonParseResult parsed = Persistence::ParseJson(json);
		REQUIRE(parsed.result.success);
		REQUIRE(parsed.value.type == Persistence::JsonValueType::Object);

		const auto& root = parsed.value.object_value;
		REQUIRE(root.at("mml").type == Persistence::JsonValueType::Object);
		REQUIRE(root.at("shape").array_value.size() == 2);
		REQUIRE(root.at("data").array_value.size() == 3);
		REQUIRE(root.at("mml").object_value.at("format").string_value == "MML_JSON");
		REQUIRE(root.at("mml").object_value.at("version").number_value == 1.0);
		REQUIRE(root.at("metadata").object_value.at("title").string_value == "small matrix");
	}

	TEST_CASE("Persistence JSON writer emits deterministic escaped JSON", "[Persistence][JSON]") {
		Persistence::JsonValue::Object object;
		object["zeta"] = Persistence::JsonValue::Number(2.5);
		object["alpha"] = Persistence::JsonValue::String("line\nquote\"");
		object["items"] = Persistence::JsonValue::ArrayValue({
			Persistence::JsonValue::Bool(true),
			Persistence::JsonValue::Null()
		});

		Persistence::SaveOptions options;
		options.pretty_json = false;
		options.precision = 6;

		const std::string written = Persistence::ToJsonString(Persistence::JsonValue::ObjectValue(object), options);
		REQUIRE(written == "{\"alpha\":\"line\\nquote\\\"\",\"items\":[true,null],\"zeta\":2.5}");

		Persistence::JsonParseResult parsed = Persistence::ParseJson(written);
		REQUIRE(parsed.result.success);
		REQUIRE(parsed.value.object_value.at("alpha").string_value == "line\nquote\"");
		REQUIRE(parsed.value.object_value.at("items").array_value[0].bool_value);
	}

	TEST_CASE("Persistence JSON parser rejects malformed input", "[Persistence][JSON]") {
		Persistence::JsonParseResult trailingComma = Persistence::ParseJson("{\"a\": 1,}");
		REQUIRE_FALSE(trailingComma.result.success);
		REQUIRE(trailingComma.result.error == IOError::MALFORMED_INPUT);

		Persistence::JsonParseResult badNumber = Persistence::ParseJson("[01]");
		REQUIRE_FALSE(badNumber.result.success);
		REQUIRE(badNumber.result.error == IOError::MALFORMED_INPUT);

		Persistence::JsonParseResult badEscape = Persistence::ParseJson("\"\\x\"");
		REQUIRE_FALSE(badEscape.result.success);
		REQUIRE(badEscape.result.error == IOError::MALFORMED_INPUT);

		Persistence::JsonParseResult loneSurrogate = Persistence::ParseJson("\"\\uD800\"");
		REQUIRE_FALSE(loneSurrogate.result.success);
		REQUIRE(loneSurrogate.result.error == IOError::MALFORMED_INPUT);
		REQUIRE(loneSurrogate.result.message.find("at byte") != std::string::npos);

		std::string deeplyNested(Persistence::MaxJsonNestingDepth + 1, '[');
		deeplyNested += '0';
		deeplyNested.append(Persistence::MaxJsonNestingDepth + 1, ']');
		Persistence::JsonParseResult excessiveDepth = Persistence::ParseJson(deeplyNested);
		REQUIRE_FALSE(excessiveDepth.result.success);
		REQUIRE(excessiveDepth.result.error == IOError::MALFORMED_INPUT);
		REQUIRE(excessiveDepth.result.message.find("Maximum JSON nesting depth exceeded") != std::string::npos);
		REQUIRE(excessiveDepth.result.message.find("at byte") != std::string::npos);
	}

	TEST_CASE("Persistence JSON writer rejects non-finite numbers", "[Persistence][JSON]") {
		std::ostringstream out;
		IOResult result = Persistence::WriteJson(out, Persistence::JsonValue::Number(std::numeric_limits<double>::infinity()));

		REQUIRE_FALSE(result.success);
		REQUIRE(result.error == IOError::UNSUPPORTED_SCALAR);
		REQUIRE(out.str().empty());
	}

	TEST_CASE("Persistence binary envelope writes and reads v1 headers", "[Persistence][Binary]") {
		Persistence::BinaryEnvelopeHeader header;
		header.object_kind = Persistence::BinaryObjectKind::Matrix;
		header.object_schema_version = 2;
		header.scalar_type = Persistence::BinaryScalarType::Float64;
		header.scalar_byte_size = 8;
		header.payload_byte_count = 128;
		header.sidecar_policy = Persistence::BinarySidecarPolicy::Optional;

		std::stringstream buffer(std::ios::in | std::ios::out | std::ios::binary);
		IOResult writeResult = Persistence::WriteBinaryEnvelopeHeader(buffer, header);
		REQUIRE(writeResult.success);
		REQUIRE(buffer.str().size() == Persistence::BinaryEnvelope::HEADER_SIZE);
		REQUIRE(buffer.str().substr(0, 8) == "MML_MATX");

		buffer.seekg(0);
		Persistence::BinaryEnvelopeReadResult read = Persistence::ReadBinaryEnvelopeHeader(buffer);
		REQUIRE(read.result.success);
		REQUIRE(read.header.container_version == Persistence::BinaryEnvelope::VERSION);
		REQUIRE(read.header.header_size == Persistence::BinaryEnvelope::HEADER_SIZE);
		REQUIRE(read.header.object_kind == Persistence::BinaryObjectKind::Matrix);
		REQUIRE(read.header.object_schema_version == 2);
		REQUIRE(read.header.scalar_type == Persistence::BinaryScalarType::Float64);
		REQUIRE(read.header.scalar_byte_size == 8);
		REQUIRE(read.header.endian_marker == Persistence::BinaryEnvelope::ENDIAN_MARKER);
		REQUIRE(read.header.payload_byte_count == 128);
		REQUIRE(read.header.flags == 0);
		REQUIRE(read.header.sidecar_policy == Persistence::BinarySidecarPolicy::Optional);
	}

	TEST_CASE("Persistence binary envelope rejects invalid headers", "[Persistence][Binary]") {
		auto makeHeaderBytes = []() {
			Persistence::BinaryEnvelopeHeader header;
			header.object_kind = Persistence::BinaryObjectKind::Matrix;
			header.scalar_type = Persistence::BinaryScalarType::Float64;
			header.scalar_byte_size = 8;

			std::stringstream buffer(std::ios::in | std::ios::out | std::ios::binary);
			REQUIRE(Persistence::WriteBinaryEnvelopeHeader(buffer, header).success);
			return buffer.str();
		};

		auto writeU16 = [](std::string& bytes, std::size_t offset, uint16_t value) {
			bytes[offset] = static_cast<char>(value & 0xFF);
			bytes[offset + 1] = static_cast<char>((value >> 8) & 0xFF);
		};

		auto writeU32 = [](std::string& bytes, std::size_t offset, uint32_t value) {
			bytes[offset] = static_cast<char>(value & 0xFF);
			bytes[offset + 1] = static_cast<char>((value >> 8) & 0xFF);
			bytes[offset + 2] = static_cast<char>((value >> 16) & 0xFF);
			bytes[offset + 3] = static_cast<char>((value >> 24) & 0xFF);
		};

		{
			std::string bytes = makeHeaderBytes();
			bytes[0] = 'X';
			std::stringstream buffer(bytes, std::ios::in | std::ios::out | std::ios::binary);
			auto read = Persistence::ReadBinaryEnvelopeHeader(buffer);
			REQUIRE_FALSE(read.result.success);
			REQUIRE(read.result.error == IOError::INVALID_FORMAT);
		}

		{
			std::string bytes = makeHeaderBytes();
			writeU16(bytes, 8, Persistence::BinaryEnvelope::VERSION + 1);
			std::stringstream buffer(bytes, std::ios::in | std::ios::out | std::ios::binary);
			auto read = Persistence::ReadBinaryEnvelopeHeader(buffer);
			REQUIRE_FALSE(read.result.success);
			REQUIRE(read.result.error == IOError::UNSUPPORTED_VERSION);
		}

		{
			std::string bytes = makeHeaderBytes();
			writeU32(bytes, 24, 4);
			std::stringstream buffer(bytes, std::ios::in | std::ios::out | std::ios::binary);
			auto read = Persistence::ReadBinaryEnvelopeHeader(buffer);
			REQUIRE_FALSE(read.result.success);
			REQUIRE(read.result.error == IOError::UNSUPPORTED_SCALAR);
		}

		{
			std::string bytes = makeHeaderBytes();
			writeU32(bytes, 28, 0x04030201);
			std::stringstream buffer(bytes, std::ios::in | std::ios::out | std::ios::binary);
			auto read = Persistence::ReadBinaryEnvelopeHeader(buffer);
			REQUIRE_FALSE(read.result.success);
			REQUIRE(read.result.error == IOError::ENDIAN_MISMATCH);
		}

		{
			std::string bytes = makeHeaderBytes().substr(0, Persistence::BinaryEnvelope::HEADER_SIZE - 1);
			std::stringstream buffer(bytes, std::ios::in | std::ios::out | std::ios::binary);
			auto read = Persistence::ReadBinaryEnvelopeHeader(buffer);
			REQUIRE_FALSE(read.result.success);
			REQUIRE(read.result.error == IOError::TRUNCATED_INPUT);
		}
	}

	TEST_CASE("Persistence binary envelope enforces object and allocation validation", "[Persistence][Binary]") {
		Persistence::BinaryEnvelopeHeader header;
		header.object_kind = Persistence::BinaryObjectKind::Matrix;
		header.scalar_type = Persistence::BinaryScalarType::Float64;
		header.scalar_byte_size = 8;
		header.payload_byte_count = 2048;

		Persistence::LoadOptions options;
		options.max_allocation_bytes = 1024;
		IOResult allocation = Persistence::ValidateBinaryEnvelopeHeader(header, options);
		REQUIRE_FALSE(allocation.success);
		REQUIRE(allocation.error == IOError::ALLOCATION_LIMIT_EXCEEDED);

		header.payload_byte_count = 16;
		header.magic = Persistence::BinaryMagic(Persistence::BinaryObjectKind::Vector);
		IOResult mismatch = Persistence::ValidateBinaryEnvelopeHeader(header);
		REQUIRE_FALSE(mismatch.success);
		REQUIRE(mismatch.error == IOError::TYPE_MISMATCH);
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         GOLDEN FILE TEST INFRASTRUCTURE                       ///
	/////////////////////////////////////////////////////////////////////////////////////

	// Path to reference files - works from both:
	// - CTest (runs from build/tests/): ../../tests/tools/serialization_test_data/
	// - VS Code debug (runs from workspace root): tests/tools/serialization_test_data/
	std::string GetReferenceDir() {
		// Try workspace-relative path first (for VS Code debug)
		std::string workspaceRelative = "tests/tools/serialization_test_data/";
		std::ifstream testFile(workspaceRelative + "multi_realfunc_simple.mml");
		if (testFile.good()) {
			return workspaceRelative;
		}
		// Fall back to CTest-relative path (from build/tests/)
		return "../../tests/tools/serialization_test_data/";
	}

	// Helper to get reference file path
	std::string GetReferenceFilePath(const std::string& filename) {
		static std::string refDir = GetReferenceDir();
		return refDir + filename;
	}

	// Parse a line into numeric tokens
	std::vector<Real> ParseNumericLine(const std::string& line) {
		std::vector<Real> values;
		std::istringstream iss(line);
		Real val;
		while (iss >> val) {
			values.push_back(val);
		}
		return values;
	}

	// Check if a line looks like data (starts with number or minus sign)
	bool IsDataLine(const std::string& line) {
		if (line.empty()) return false;
		char c = line[0];
		return std::isdigit(c) || c == '-' || c == '+' || c == '.';
	}

	// Compare two files with hybrid approach:
	// - Exact match for header lines (type, title, metadata)
	// - Numerical tolerance for data lines
	struct CompareResult {
		bool success = true;
		std::string error_message;
		int line_number = 0;
	};

	CompareResult CompareFilesWithTolerance(const std::string& actualPath, 
	                                         const std::string& referencePath,
	                                         Real tolerance = TOL(1e-10, 1e-5)) {
		CompareResult result;
		
		std::ifstream actualFile(actualPath);
		std::ifstream referenceFile(referencePath);
		
		if (!actualFile.is_open()) {
			result.success = false;
			result.error_message = "Cannot open actual file: " + actualPath;
			return result;
		}
		
		if (!referenceFile.is_open()) {
			result.success = false;
			result.error_message = "Cannot open reference file: " + referencePath;
			return result;
		}
		
		std::string actualLine, refLine;
		int lineNum = 0;
		
		while (std::getline(referenceFile, refLine)) {
			lineNum++;
			
			if (!std::getline(actualFile, actualLine)) {
				result.success = false;
				result.line_number = lineNum;
				result.error_message = "Actual file has fewer lines than reference (expected line: " + refLine + ")";
				return result;
			}
			
			if (IsDataLine(refLine) && IsDataLine(actualLine)) {
				// Numerical comparison with tolerance
				auto refValues = ParseNumericLine(refLine);
				auto actualValues = ParseNumericLine(actualLine);
				
				if (refValues.size() != actualValues.size()) {
					result.success = false;
					result.line_number = lineNum;
					result.error_message = "Column count mismatch: expected " + 
					    std::to_string(refValues.size()) + ", got " + 
					    std::to_string(actualValues.size());
					return result;
				}
				
				for (size_t i = 0; i < refValues.size(); ++i) {
					Real diff = std::abs(refValues[i] - actualValues[i]);
					// Use relative tolerance for large values, absolute for small
					Real relativeTol = std::max(tolerance, tolerance * std::abs(refValues[i]));
					if (diff > relativeTol) {
						result.success = false;
						result.line_number = lineNum;
						std::ostringstream oss;
						oss << "Numeric mismatch in column " << i << ": expected " 
						    << refValues[i] << ", got " << actualValues[i] 
						    << " (diff=" << diff << ", tol=" << relativeTol << ")";
						result.error_message = oss.str();
						return result;
					}
				}
			} else {
				// Exact string comparison for headers
				// Trim trailing whitespace for comparison
				auto rtrim = [](std::string& s) {
					s.erase(std::find_if(s.rbegin(), s.rend(), 
					        [](unsigned char ch) { return !std::isspace(ch); }).base(), s.end());
				};
				
				std::string trimmedActual = actualLine;
				std::string trimmedRef = refLine;
				rtrim(trimmedActual);
				rtrim(trimmedRef);
				
				if (trimmedActual != trimmedRef) {
					result.success = false;
					result.line_number = lineNum;
					result.error_message = "Header mismatch: expected '" + trimmedRef + 
					    "', got '" + trimmedActual + "'";
					return result;
				}
			}
		}
		
		// Check if actual file has extra lines
		if (std::getline(actualFile, actualLine)) {
			result.success = false;
			result.line_number = lineNum + 1;
			result.error_message = "Actual file has more lines than reference (extra: " + actualLine + ")";
			return result;
		}
		
		return result;
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         TEST HELPER FUNCTIONS                                 ///
	/////////////////////////////////////////////////////////////////////////////////////

	// Helper to create a temp file path
	std::string GetTempFilePath(const std::string& suffix) {
		std::filesystem::path tempDir = std::filesystem::temp_directory_path();
		return (tempDir / ("mml_serializer_test_" + suffix + ".mml")).string();
	}

	// Helper to read file contents
	std::string ReadFileContents(const std::string& path) {
		std::ifstream file(path);
		std::stringstream buffer;
		buffer << file.rdbuf();
		return buffer.str();
	}

	// Helper to clean up temp file
	void CleanupTempFile(const std::string& path) {
		std::filesystem::remove(path);
	}

	TEST_CASE("Serializer text headers sanitize multiline metadata", "[serializer][headers]") {
		const std::string title = "first line\r\nx1: injected";
		const std::string legend = "series\nNumPoints: injected";
		const std::string sanitizedTitle = "first line  x1: injected";
		const std::string sanitizedLegend = "series NumPoints: injected";

		REQUIRE(Serializer::SanitizeHeaderLine(title) == sanitizedTitle);
		REQUIRE(Serializer::SanitizeHeaderLine("ordinary title") == "ordinary title");

		auto readLines = [](const std::string& content) {
			std::istringstream input(content);
			std::vector<std::string> lines;
			for (std::string line; std::getline(input, line);)
				lines.push_back(line);
			return lines;
		};

		SECTION("real function") {
			std::ostringstream output;
			REQUIRE(Serializer::WriteRealFuncHeader(output, SerializeFormatType::REAL_FUNCTION, title, 0.0, 1.0, 2).success);
			const auto lines = readLines(output.str());
			REQUIRE(lines.size() == 6);
			REQUIRE(lines[2] == sanitizedTitle);
			REQUIRE(lines[3] == "x1: 0");
		}

		SECTION("multiple real functions and legends") {
			std::ostringstream output;
			REQUIRE(Serializer::WriteRealMultiFuncHeader(output, title, 1, {legend}, 0.0, 1.0, 2).success);
			const auto lines = readLines(output.str());
			REQUIRE(lines.size() == 8);
			REQUIRE(lines[2] == sanitizedTitle);
			REQUIRE(lines[4] == sanitizedLegend);
			REQUIRE(lines[5] == "x1: 0");
		}

		SECTION("parametric curve") {
			std::ostringstream output;
			REQUIRE(Serializer::WriteParamCurveHeader(output, SerializeFormatType::PARAMETRIC_CURVE_CARTESIAN_2D, title, 0.0, 1.0, 2).success);
			const auto lines = readLines(output.str());
			REQUIRE(lines.size() == 6);
			REQUIRE(lines[2] == sanitizedTitle);
			REQUIRE(lines[3] == "t1: 0");
		}

		SECTION("vector field") {
			std::ostringstream output;
			REQUIRE(Serializer::WriteVectorFieldHeader(output, SerializeFormatType::VECTOR_FIELD_2D_CARTESIAN, title).success);
			const auto lines = readLines(output.str());
			REQUIRE(lines.size() == 3);
			REQUIRE(lines[2] == sanitizedTitle);
		}

		SECTION("field lines") {
			std::ostringstream output;
			REQUIRE(Serializer::WriteFieldLinesHeader(output, SerializeFormatType::FIELD_LINES_2D, title, 3).success);
			const auto lines = readLines(output.str());
			REQUIRE(lines.size() == 4);
			REQUIRE(lines[2] == "Title: " + sanitizedTitle);
			REQUIRE(lines[3] == "NUM_LINES: 3");
		}

		SECTION("parametric surface") {
			std::ostringstream output;
			REQUIRE(Serializer::WriteParametricSurfaceHeader(output, title, 0.0, 1.0, 2, -1.0, 1.0, 3).success);
			const auto lines = readLines(output.str());
			REQUIRE(lines.size() == 9);
			REQUIRE(lines[2] == sanitizedTitle);
			REQUIRE(lines[8] == "NumPointsW: 3");
		}

		SECTION("2D scalar function") {
			std::ostringstream output;
			REQUIRE(Serializer::WriteScalarFunc2DHeader(output, title, 0.0, 1.0, 2, -1.0, 1.0, 3).success);
			const auto lines = readLines(output.str());
			REQUIRE(lines.size() == 9);
			REQUIRE(lines[2] == sanitizedTitle);
			REQUIRE(lines[8] == "NumPointsY: 3");
		}

		SECTION("3D scalar function") {
			std::ostringstream output;
			REQUIRE(Serializer::WriteScalarFunc3DHeader(output, title, 0.0, 1.0, 2, -1.0, 1.0, 3, 2.0, 4.0, 5).success);
			const auto lines = readLines(output.str());
			REQUIRE(lines.size() == 12);
			REQUIRE(lines[2] == sanitizedTitle);
			REQUIRE(lines[11] == "NumPointsZ: 5");
		}

		SECTION("2D particle simulation") {
			std::ostringstream output;
			REQUIRE(Serializer::WriteParticleSimulation2DHeader(output, 1, 10.0, 20.0, {"red"}, {0.5}, 4).success);
			const auto lines = readLines(output.str());
			REQUIRE(lines.size() == 7);
			REQUIRE(lines[5] == "Ball_1 red 0.5");
			REQUIRE(lines[6] == "NumSteps: 4");
		}

		SECTION("3D particle simulation") {
			std::ostringstream output;
			REQUIRE(Serializer::WriteParticleSimulation3DHeader(output, 1, 10.0, 20.0, 30.0, {"blue"}, {0.25}, 5).success);
			const auto lines = readLines(output.str());
			REQUIRE(lines.size() == 8);
			REQUIRE(lines[6] == "Ball_1 blue 0.25");
			REQUIRE(lines[7] == "NumSteps: 5");
		}
	}

	TEMPLATE_TEST_CASE("Persistence binary round-trips every canonical Vector and Matrix scalar", "[Persistence][Binary][ScalarTypes]",
	                   float, double, long double, (std::complex<float>), (std::complex<double>)) {
		auto value = [](long double real, long double imaginary = 0) {
			if constexpr (std::is_same_v<TestType, std::complex<float>> || std::is_same_v<TestType, std::complex<double>>)
				return TestType(static_cast<typename TestType::value_type>(real), static_cast<typename TestType::value_type>(imaginary));
			else
				return static_cast<TestType>(real);
		};

		Vector<TestType> vec{ value(1.25L, -0.5L), value(-2.5L, 4.0L), value(3.75L, -8.125L) };
		std::stringstream vecStream(std::ios::in | std::ios::out | std::ios::binary);
		REQUIRE(Persistence::SaveBinary(vecStream, vec).success);
		REQUIRE(vecStream.str().substr(0, 8) == "MML_VECT");
		vecStream.seekg(0);
		auto vecEnvelope = Persistence::ReadBinaryEnvelopeHeader(vecStream);
		REQUIRE(vecEnvelope.result.success);
		REQUIRE(vecEnvelope.header.object_kind == Persistence::BinaryObjectKind::Vector);
		REQUIRE(vecEnvelope.header.scalar_type == Persistence::Detail::BinaryScalarTypeFor<TestType>());
		REQUIRE(vecEnvelope.header.scalar_byte_size == Persistence::Detail::BinaryScalarByteSize<TestType>());
		vecStream.seekg(0);

		Vector<TestType> loadedVec;
		REQUIRE(Persistence::LoadBinary(vecStream, loadedVec).success);
		REQUIRE(loadedVec.size() == vec.size());
		for (int i = 0; i < vec.size(); ++i)
			REQUIRE(loadedVec[i] == vec[i]);

		Matrix<TestType> mat(2, 3);
		mat(0, 0) = value(1.0L, -1.5L);
		mat(0, 1) = value(2.0L, 2.25L);
		mat(0, 2) = value(3.0L, -3.5L);
		mat(1, 0) = value(-4.0L, 4.75L);
		mat(1, 1) = value(5.5L, -5.0L);
		mat(1, 2) = value(-6.25L, 6.125L);
		std::stringstream matStream(std::ios::in | std::ios::out | std::ios::binary);
		REQUIRE(Persistence::SaveBinary(matStream, mat).success);
		REQUIRE(matStream.str().substr(0, 8) == "MML_MATX");
		matStream.seekg(0);
		auto matEnvelope = Persistence::ReadBinaryEnvelopeHeader(matStream);
		REQUIRE(matEnvelope.result.success);
		REQUIRE(matEnvelope.header.object_kind == Persistence::BinaryObjectKind::Matrix);
		REQUIRE(matEnvelope.header.scalar_type == Persistence::Detail::BinaryScalarTypeFor<TestType>());
		REQUIRE(matEnvelope.header.scalar_byte_size == Persistence::Detail::BinaryScalarByteSize<TestType>());
		matStream.seekg(0);

		Matrix<TestType> loadedMat;
		REQUIRE(Persistence::LoadBinary(matStream, loadedMat).success);
		REQUIRE(loadedMat.rows() == mat.rows());
		REQUIRE(loadedMat.cols() == mat.cols());
		for (int i = 0; i < mat.rows(); ++i)
			for (int j = 0; j < mat.cols(); ++j)
				REQUIRE(loadedMat(i, j) == mat(i, j));

		matStream.clear();
		matStream.seekg(0);
		if constexpr (std::is_same_v<TestType, float>) {
			Matrix<double> mismatched;
			IOResult mismatch = Persistence::LoadBinary(matStream, mismatched);
			REQUIRE_FALSE(mismatch.success);
			REQUIRE(mismatch.error == IOError::UNSUPPORTED_SCALAR);
		} else {
			Matrix<float> mismatched;
			IOResult mismatch = Persistence::LoadBinary(matStream, mismatched);
			REQUIRE_FALSE(mismatch.success);
			REQUIRE(mismatch.error == IOError::UNSUPPORTED_SCALAR);
		}
	}

	TEST_CASE("Persistence canonical long double preserves edge values", "[Persistence][Binary][LongDouble]") {
		Vector<long double> values{
			0.0L,
			-0.0L,
			std::numeric_limits<long double>::denorm_min(),
			std::numeric_limits<long double>::max(),
			std::numeric_limits<long double>::infinity(),
			-std::numeric_limits<long double>::infinity(),
			std::numeric_limits<long double>::quiet_NaN()
		};
		std::stringstream stream(std::ios::in | std::ios::out | std::ios::binary);
		REQUIRE(Persistence::SaveBinary(stream, values).success);
		stream.seekg(0);

		Vector<long double> loaded;
		REQUIRE(Persistence::LoadBinary(stream, loaded).success);
		REQUIRE(loaded.size() == values.size());
		REQUIRE(loaded[0] == 0.0L);
		REQUIRE_FALSE(std::signbit(loaded[0]));
		REQUIRE(loaded[1] == 0.0L);
		REQUIRE(std::signbit(loaded[1]));
		REQUIRE(loaded[2] == values[2]);
		REQUIRE(loaded[3] == values[3]);
		REQUIRE(loaded[4] == values[4]);
		REQUIRE(loaded[5] == values[5]);
		REQUIRE(std::isnan(loaded[6]));
	}

	TEST_CASE("Persistence converts long double real vectors through complex double", "[Persistence][Binary][LongDouble]") {
		std::vector<long double> standardValues{ 1.25L, -2.5L, 3.75L };
		std::string standardPath = (std::filesystem::temp_directory_path() / "mml_serializer_std_long_double_as_complex.mmlb").string();
		REQUIRE(Persistence::SaveRealAsComplex(standardValues, standardPath).success);

		std::vector<long double> loadedStandardValues;
		REQUIRE(Persistence::LoadComplexAsReal(standardPath, loadedStandardValues).success);
		REQUIRE(loadedStandardValues == standardValues);
		CleanupTempFile(standardPath);

		Vector<long double> values{ 1.25L, -2.5L, 3.75L };
		std::string path = (std::filesystem::temp_directory_path() / "mml_serializer_long_double_as_complex.mmlb").string();
		REQUIRE(Persistence::SaveVectorAsComplex(values, path).success);

		std::vector<std::complex<double>> loaded;
		REQUIRE(Persistence::LoadBinary(path, loaded).success);
		REQUIRE(loaded == std::vector<std::complex<double>>{ { 1.25, 0.0 }, { -2.5, 0.0 }, { 3.75, 0.0 } });

		Vector<long double> loadedValues;
		REQUIRE(Persistence::LoadComplexAsVector(path, loadedValues).success);
		REQUIRE(loadedValues.size() == values.size());
		for (int index = 0; index < values.size(); ++index)
			REQUIRE(loadedValues[index] == values[index]);
		CleanupTempFile(path);
	}

	TEST_CASE("Persistence binary supports generic mmlb path dispatch", "[Persistence][Binary]") {
		Vector<float> vec{ 1.5f, -2.25f, 3.125f };
		std::string vecPath = (std::filesystem::temp_directory_path() / "mml_serializer_vector_binary.mmlb").string();
		REQUIRE(Persistence::Save(vec, vecPath).success);

		Vector<float> loadedVec;
		REQUIRE(Persistence::Load(vecPath, loadedVec).success);
		REQUIRE(loadedVec.size() == vec.size());
		for (int i = 0; i < vec.size(); ++i)
			REQUIRE(loadedVec[i] == vec[i]);
		CleanupTempFile(vecPath);

		Matrix<Real> mat(2, 2);
		mat(0, 0) = 1.0;
		mat(0, 1) = 2.0;
		mat(1, 0) = 3.0;
		mat(1, 1) = 4.0;
		std::string matPath = (std::filesystem::temp_directory_path() / "mml_serializer_matrix_binary.mmlb").string();
		REQUIRE(Persistence::Save(mat, matPath).success);

		Matrix<Real> loadedMat;
		REQUIRE(Persistence::Load(matPath, loadedMat).success);
		REQUIRE(loadedMat.rows() == mat.rows());
		REQUIRE(loadedMat.cols() == mat.cols());
		for (int i = 0; i < mat.rows(); ++i)
			for (int j = 0; j < mat.cols(); ++j)
				REQUIRE(loadedMat(i, j) == mat(i, j));
		CleanupTempFile(matPath);
	}

	TEST_CASE("Persistence binary rejects corrupted payloads", "[Persistence][Binary]") {
		Vector<Real> vec{ 1.0, 2.0 };
		std::stringstream valid(std::ios::in | std::ios::out | std::ios::binary);
		REQUIRE(Persistence::SaveBinary(valid, vec).success);

		std::string badPayloadSize = valid.str();
		for (int i = 0; i < 8; ++i)
			badPayloadSize[32 + i] = 0;
		badPayloadSize[32] = static_cast<char>(9);
		std::stringstream badPayloadStream(badPayloadSize, std::ios::in | std::ios::out | std::ios::binary);
		Vector<Real> payloadRejected;
		IOResult payloadResult = Persistence::LoadBinary(badPayloadStream, payloadRejected);
		REQUIRE_FALSE(payloadResult.success);
		REQUIRE(payloadResult.error == IOError::SCHEMA_MISMATCH);

		std::string truncated = valid.str().substr(0, valid.str().size() - 1);
		std::stringstream truncatedStream(truncated, std::ios::in | std::ios::out | std::ios::binary);
		Vector<Real> truncatedRejected;
		IOResult truncatedResult = Persistence::LoadBinary(truncatedStream, truncatedRejected);
		REQUIRE_FALSE(truncatedResult.success);
		REQUIRE(truncatedResult.error == IOError::TRUNCATED_INPUT);
	}

	TEST_CASE("Persistence dense load guards reject JSON allocation limits", "[Persistence][DenseGuards]") {
		const std::string vectorJson = WithRealMetadata(R"json({
			"mml": {
				"format": "MML_JSON",
				"version": 1,
				"object": "Vector",
				"object_version": 1,
				"scalar": "Real",
				"scalar_bytes": 8,
				"scalar_encoding": "ieee754",
				"real_type": "double"
			},
			"shape": [2],
			"data": [1, 2]
		})json");

		Persistence::LoadOptions tiny;
		tiny.max_allocation_bytes = 2 * sizeof(Real) - 1;
		std::stringstream vectorStream(vectorJson);
		Vector<Real> vectorRejected;
		IOResult vectorResult = Persistence::LoadJson(vectorStream, vectorRejected, tiny);
		REQUIRE_FALSE(vectorResult.success);
		REQUIRE(vectorResult.error == IOError::ALLOCATION_LIMIT_EXCEEDED);
		REQUIRE(vectorRejected.size() == 0);

		const std::string matrixJson = WithRealMetadata(R"json({
			"mml": {
				"format": "MML_JSON",
				"version": 1,
				"object": "Matrix",
				"object_version": 1,
				"scalar": "Real",
				"scalar_bytes": 8,
				"scalar_encoding": "ieee754",
				"layout": "row_major",
				"real_type": "double"
			},
			"shape": [2, 2],
			"data": [1, 2, 3, 4]
		})json");

		tiny.max_allocation_bytes = 4 * sizeof(Real) - 1;
		std::stringstream matrixStream(matrixJson);
		Matrix<Real> matrixRejected;
		IOResult matrixResult = Persistence::LoadJson(matrixStream, matrixRejected, tiny);
		REQUIRE_FALSE(matrixResult.success);
		REQUIRE(matrixResult.error == IOError::ALLOCATION_LIMIT_EXCEEDED);
		REQUIRE(matrixRejected.rows() == 0);
		REQUIRE(matrixRejected.cols() == 0);
	}

	TEST_CASE("Persistence dense load guards reject binary allocation limits and overflow", "[Persistence][DenseGuards]") {
		Vector<Real> vector{ 1.0, 2.0 };
		std::stringstream vectorBuffer(std::ios::in | std::ios::out | std::ios::binary);
		REQUIRE(Persistence::SaveBinary(vectorBuffer, vector).success);

		Persistence::LoadOptions tiny;
		tiny.max_allocation_bytes = 2 * sizeof(Real) - 1;
		vectorBuffer.seekg(0);
		Vector<Real> vectorRejected;
		IOResult vectorResult = Persistence::LoadBinary(vectorBuffer, vectorRejected, tiny);
		REQUIRE_FALSE(vectorResult.success);
		REQUIRE(vectorResult.error == IOError::ALLOCATION_LIMIT_EXCEEDED);
		REQUIRE(vectorRejected.size() == 0);

		Matrix<Real> matrix(2, 2);
		matrix(0, 0) = 1.0;
		matrix(0, 1) = 2.0;
		matrix(1, 0) = 3.0;
		matrix(1, 1) = 4.0;
		std::stringstream matrixBuffer(std::ios::in | std::ios::out | std::ios::binary);
		REQUIRE(Persistence::SaveBinary(matrixBuffer, matrix).success);
		matrixBuffer.seekg(0);

		Matrix<Real> matrixRejected;
		IOResult matrixResult = Persistence::LoadBinary(matrixBuffer, matrixRejected, tiny);
		REQUIRE_FALSE(matrixResult.success);
		REQUIRE(matrixResult.error == IOError::ALLOCATION_LIMIT_EXCEEDED);
		REQUIRE(matrixRejected.rows() == 0);
		REQUIRE(matrixRejected.cols() == 0);

		Persistence::BinaryEnvelopeHeader overflowHeader;
		overflowHeader.object_kind = Persistence::BinaryObjectKind::Matrix;
		overflowHeader.scalar_type = Persistence::Detail::BinaryScalarTypeFor<Real>();
		overflowHeader.scalar_byte_size = Persistence::Detail::BinaryScalarByteSize<Real>();
		overflowHeader.payload_byte_count = 16;
		std::stringstream overflowBuffer(std::ios::in | std::ios::out | std::ios::binary);
		REQUIRE(Persistence::WriteBinaryEnvelopeHeader(overflowBuffer, overflowHeader).success);
		Persistence::WriteUInt64LE(overflowBuffer, static_cast<uint64_t>(std::numeric_limits<int>::max()) + 1);
		Persistence::WriteUInt64LE(overflowBuffer, 1);
		overflowBuffer.seekg(0);

		Persistence::LoadOptions huge;
		huge.max_allocation_bytes = std::numeric_limits<std::size_t>::max();
		Matrix<Real> overflowRejected;
		IOResult overflowResult = Persistence::LoadBinary(overflowBuffer, overflowRejected, huge);
		REQUIRE_FALSE(overflowResult.success);
		REQUIRE(overflowResult.error == IOError::ALLOCATION_LIMIT_EXCEEDED);
	}

	TEST_CASE("Persistence dense load guards reject truncated matrix payload without allocation surprises", "[Persistence][DenseGuards]") {
		Matrix<Real> matrix(2, 2);
		matrix(0, 0) = 1.0;
		matrix(0, 1) = 2.0;
		matrix(1, 0) = 3.0;
		matrix(1, 1) = 4.0;
		std::stringstream valid(std::ios::in | std::ios::out | std::ios::binary);
		REQUIRE(Persistence::SaveBinary(valid, matrix).success);

		std::string truncated = valid.str().substr(0, valid.str().size() - 1);
		std::stringstream truncatedStream(truncated, std::ios::in | std::ios::out | std::ios::binary);
		Matrix<Real> rejected;
		IOResult result = Persistence::LoadBinary(truncatedStream, rejected);
		REQUIRE_FALSE(result.success);
		REQUIRE(result.error == IOError::TRUNCATED_INPUT);
	}

	TEST_CASE("Persistence Vector JSON round-trips through streams and paths", "[Persistence][VectorJSON]") {
		Vector<Real> original{ 1.25, -2.5, 3.75 };

		std::stringstream stream;
		REQUIRE(Persistence::SaveJson(stream, original).success);

		Vector<Real> loaded;
		REQUIRE(Persistence::LoadJson(stream, loaded).success);
		REQUIRE(loaded.size() == original.size());
		for (int i = 0; i < original.size(); ++i)
			REQUIRE(loaded[i] == original[i]);

		std::string path = (std::filesystem::temp_directory_path() / "mml_serializer_vector_roundtrip.mmlj").string();
		REQUIRE(Persistence::SaveJson(original, path).success);

		Vector<Real> fileLoaded;
		REQUIRE(Persistence::LoadJson(path, fileLoaded).success);
		REQUIRE(fileLoaded.size() == original.size());
		for (int i = 0; i < original.size(); ++i)
			REQUIRE(fileLoaded[i] == original[i]);
		CleanupTempFile(path);
	}

	TEST_CASE("Persistence Vector JSON supports float vectors and generic path dispatch", "[Persistence][VectorJSON]") {
		Vector<float> original{ 1.5f, -2.25f, 3.125f };
		std::string path = (std::filesystem::temp_directory_path() / "mml_serializer_vector_float.mmlj").string();

		REQUIRE(Persistence::Save(original, path).success);

		Vector<float> loaded;
		REQUIRE(Persistence::Load(path, loaded).success);
		REQUIRE(loaded.size() == original.size());
		for (int i = 0; i < original.size(); ++i)
			REQUIRE(loaded[i] == original[i]);

		CleanupTempFile(path);
	}

	TEST_CASE("Persistence VectorN JSON round-trips and validates fixed shape", "[Persistence][VectorJSON]") {
		VectorN<Real, 3> original{ 4.0, 5.0, 6.0 };

		std::stringstream stream;
		REQUIRE(Persistence::SaveJson(stream, original).success);

		VectorN<Real, 3> loaded;
		REQUIRE(Persistence::LoadJson(stream, loaded).success);
		for (int i = 0; i < 3; ++i)
			REQUIRE(loaded[i] == original[i]);

		std::string wrongShape = stream.str();
		const std::string from = "\"shape\": [\n    3\n  ]";
		const std::string to = "\"shape\": [\n    2\n  ]";
		const std::size_t pos = wrongShape.find(from);
		REQUIRE(pos != std::string::npos);
		wrongShape.replace(pos, from.size(), to);

		std::stringstream wrongShapeStream(wrongShape);
		VectorN<Real, 3> rejected;
		IOResult result = Persistence::LoadJson(wrongShapeStream, rejected);
		REQUIRE_FALSE(result.success);
		REQUIRE(result.error == IOError::SCHEMA_MISMATCH);
	}

	TEST_CASE("Persistence Vector JSON produces deterministic compact output", "[Persistence][VectorJSON]") {
		Vector<Real> vec{ Real(1.0) / Real(3.0), Real(2.0) };
		Persistence::SaveOptions options;
		options.pretty_json = false;
		options.precision = 4;

		std::stringstream stream;
		REQUIRE(Persistence::SaveJson(stream, vec, options).success);
		REQUIRE(stream.str() == WithRealMetadata("{\"data\":[0.3333,2],\"mml\":{\"format\":\"MML_JSON\",\"object\":\"Vector\",\"object_version\":1,\"real_type\":\"double\",\"scalar\":\"Real\",\"scalar_bytes\":8,\"scalar_encoding\":\"ieee754\",\"version\":1},\"shape\":[2]}"));
	}

	TEST_CASE("Persistence Vector JSON rejects malformed wrong-scalar and wrong-shape input", "[Persistence][VectorJSON]") {
		Vector<Real> original{ 1.0, 2.0 };
		std::stringstream stream;
		REQUIRE(Persistence::SaveJson(stream, original).success);
		const std::string good = stream.str();

		std::string wrongScalar = good;
		const std::size_t scalarPos = wrongScalar.find("\"scalar\": \"Real\"");
		REQUIRE(scalarPos != std::string::npos);
		wrongScalar.replace(scalarPos, std::string("\"scalar\": \"Real\"").size(), "\"scalar\": \"float\"");
		std::stringstream wrongScalarStream(wrongScalar);
		Vector<Real> scalarRejected;
		IOResult scalarResult = Persistence::LoadJson(wrongScalarStream, scalarRejected);
		REQUIRE_FALSE(scalarResult.success);
		REQUIRE(scalarResult.error == IOError::UNSUPPORTED_SCALAR);

		const std::string wrongLength = WithRealMetadata(R"json({
			"mml": {
				"format": "MML_JSON",
				"version": 1,
				"object": "Vector",
				"object_version": 1,
				"scalar": "Real",
				"scalar_bytes": 8,
				"scalar_encoding": "ieee754",
				"real_type": "double"
			},
			"shape": [2],
			"data": [1]
		})json");
		std::stringstream wrongLengthStream(wrongLength);
		Vector<Real> lengthRejected;
		IOResult lengthResult = Persistence::LoadJson(wrongLengthStream, lengthRejected);
		REQUIRE_FALSE(lengthResult.success);
		REQUIRE(lengthResult.error == IOError::SCHEMA_MISMATCH);

		std::stringstream malformed("{\"mml\":");
		Vector<Real> malformedRejected;
		IOResult malformedResult = Persistence::LoadJson(malformed, malformedRejected);
		REQUIRE_FALSE(malformedResult.success);
		REQUIRE(malformedResult.error == IOError::MALFORMED_INPUT);
	}

	TEST_CASE("Persistence Matrix JSON round-trips through streams and paths", "[Persistence][MatrixJSON]") {
		Matrix<Real> original(2, 3);
		original(0, 0) = 1.25;
		original(0, 1) = -2.5;
		original(0, 2) = 3.75;
		original(1, 0) = 4.5;
		original(1, 1) = -5.25;
		original(1, 2) = 6.125;

		std::stringstream stream;
		REQUIRE(Persistence::SaveJson(stream, original).success);

		Matrix<Real> loaded;
		REQUIRE(Persistence::LoadJson(stream, loaded).success);
		REQUIRE(loaded.rows() == original.rows());
		REQUIRE(loaded.cols() == original.cols());
		for (int i = 0; i < original.rows(); ++i)
			for (int j = 0; j < original.cols(); ++j)
				REQUIRE(loaded(i, j) == original(i, j));

		std::string path = (std::filesystem::temp_directory_path() / "mml_serializer_matrix_roundtrip.mmlj").string();
		REQUIRE(Persistence::SaveJson(original, path).success);

		Matrix<Real> fileLoaded;
		REQUIRE(Persistence::LoadJson(path, fileLoaded).success);
		REQUIRE(fileLoaded.rows() == original.rows());
		REQUIRE(fileLoaded.cols() == original.cols());
		for (int i = 0; i < original.rows(); ++i)
			for (int j = 0; j < original.cols(); ++j)
				REQUIRE(fileLoaded(i, j) == original(i, j));
		CleanupTempFile(path);
	}

	TEST_CASE("Persistence MatrixNM JSON round-trips and validates fixed shape", "[Persistence][MatrixJSON]") {
		MatrixNM<Real, 2, 3> original{ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };

		std::stringstream stream;
		REQUIRE(Persistence::SaveJson(stream, original).success);

		MatrixNM<Real, 2, 3> loaded;
		REQUIRE(Persistence::LoadJson(stream, loaded).success);
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 3; ++j)
				REQUIRE(loaded(i, j) == original(i, j));

		const std::string wrongShape = WithRealMetadata(R"json({
			"mml": {
				"format": "MML_JSON",
				"version": 1,
				"object": "MatrixNM",
				"object_version": 1,
				"scalar": "Real",
				"scalar_bytes": 8,
				"scalar_encoding": "ieee754",
				"layout": "row_major",
				"real_type": "double"
			},
			"shape": [2, 2],
			"data": [1, 2, 3, 4]
		})json");
		std::stringstream wrongShapeStream(wrongShape);
		MatrixNM<Real, 2, 3> rejected;
		IOResult result = Persistence::LoadJson(wrongShapeStream, rejected);
		REQUIRE_FALSE(result.success);
		REQUIRE(result.error == IOError::SCHEMA_MISMATCH);
	}

	TEST_CASE("Persistence Matrix JSON produces deterministic compact output", "[Persistence][MatrixJSON]") {
		Matrix<Real> mat(2, 2);
		mat(0, 0) = Real(1.0) / Real(3.0);
		mat(0, 1) = 2.0;
		mat(1, 0) = 3.0;
		mat(1, 1) = 4.0;
		Persistence::SaveOptions options;
		options.pretty_json = false;
		options.precision = 4;

		std::stringstream stream;
		REQUIRE(Persistence::SaveJson(stream, mat, options).success);
		REQUIRE(stream.str() == WithRealMetadata("{\"data\":[0.3333,2,3,4],\"mml\":{\"format\":\"MML_JSON\",\"layout\":\"row_major\",\"object\":\"Matrix\",\"object_version\":1,\"real_type\":\"double\",\"scalar\":\"Real\",\"scalar_bytes\":8,\"scalar_encoding\":\"ieee754\",\"version\":1},\"shape\":[2,2]}"));
	}

	TEST_CASE("Persistence Matrix JSON rejects malformed wrong-scalar wrong-layout and wrong-shape input", "[Persistence][MatrixJSON]") {
		const std::string base = WithRealMetadata(R"json({
			"mml": {
				"format": "MML_JSON",
				"version": 1,
				"object": "Matrix",
				"object_version": 1,
				"scalar": "Real",
				"scalar_bytes": 8,
				"scalar_encoding": "ieee754",
				"layout": "row_major",
				"real_type": "double"
			},
			"shape": [2, 2],
			"data": [1, 2, 3, 4]
		})json");

		std::string wrongScalar = base;
		const std::size_t scalarPos = wrongScalar.find("\"scalar\": \"Real\"");
		REQUIRE(scalarPos != std::string::npos);
		wrongScalar.replace(scalarPos, std::string("\"scalar\": \"Real\"").size(), "\"scalar\": \"float\"");
		std::stringstream wrongScalarStream(wrongScalar);
		Matrix<Real> scalarRejected;
		IOResult scalarResult = Persistence::LoadJson(wrongScalarStream, scalarRejected);
		REQUIRE_FALSE(scalarResult.success);
		REQUIRE(scalarResult.error == IOError::UNSUPPORTED_SCALAR);

		std::string wrongLayout = base;
		const std::size_t layoutPos = wrongLayout.find("\"layout\": \"row_major\"");
		REQUIRE(layoutPos != std::string::npos);
		wrongLayout.replace(layoutPos, std::string("\"layout\": \"row_major\"").size(), "\"layout\": \"column_major\"");
		std::stringstream wrongLayoutStream(wrongLayout);
		Matrix<Real> layoutRejected;
		IOResult layoutResult = Persistence::LoadJson(wrongLayoutStream, layoutRejected);
		REQUIRE_FALSE(layoutResult.success);
		REQUIRE(layoutResult.error == IOError::SCHEMA_MISMATCH);

		const std::string wrongLength = WithRealMetadata(R"json({
			"mml": {
				"format": "MML_JSON",
				"version": 1,
				"object": "Matrix",
				"object_version": 1,
				"scalar": "Real",
				"scalar_bytes": 8,
				"scalar_encoding": "ieee754",
				"layout": "row_major",
				"real_type": "double"
			},
			"shape": [2, 2],
			"data": [1, 2, 3]
		})json");
		std::stringstream wrongLengthStream(wrongLength);
		Matrix<Real> lengthRejected;
		IOResult lengthResult = Persistence::LoadJson(wrongLengthStream, lengthRejected);
		REQUIRE_FALSE(lengthResult.success);
		REQUIRE(lengthResult.error == IOError::SCHEMA_MISMATCH);

		std::stringstream malformed("{\"mml\":");
		Matrix<Real> malformedRejected;
		IOResult malformedResult = Persistence::LoadJson(malformed, malformedRejected);
		REQUIRE_FALSE(malformedResult.success);
		REQUIRE(malformedResult.error == IOError::MALFORMED_INPUT);
	}

	TEST_CASE("Persistence sampled function JSON saves IRealFunction samples and loads data", "[Persistence][FunctionJSON]") {
		RealFunction sinFunction([](const Real x) { return std::sin(x); });
		std::stringstream stream;
		Persistence::SampleGrid grid{0.0, Constants::PI / 2.0, 4};

		REQUIRE(Persistence::SaveSampledFunction(stream, sinFunction, "sin(x)", grid).success);

		Persistence::SampledRealFunctionData loaded;
		REQUIRE(Persistence::LoadSampledFunctionData(stream, loaded).success);
		REQUIRE(loaded.nodes.size() == 4);
		REQUIRE(loaded.values.size() == 4);
		REQUIRE(loaded.label == "sin(x)");
		REQUIRE(loaded.source == "sampled from IRealFunction");
		REQUIRE(std::abs(loaded.nodes[0] - 0.0) < 1e-12);
		REQUIRE(std::abs(loaded.nodes[3] - Constants::PI / 2.0) < 1e-12);
		REQUIRE(std::abs(loaded.values[0] - 0.0) < 1e-12);
		REQUIRE(std::abs(loaded.values[3] - 1.0) < 1e-12);
	}

	TEST_CASE("Persistence sampled function JSON saves interpolated function nodes", "[Persistence][FunctionJSON]") {
		Vector<Real> x(std::vector<Real>{0.0, 1.0, 2.0});
		Vector<Real> y(std::vector<Real>{0.0, 1.0, 4.0});
		LinearInterpRealFunc interp(x, y);

		std::string path = (std::filesystem::temp_directory_path() / "mml_serializer_sampled_function.mmlj").string();
		REQUIRE(Persistence::SaveInterpolatedFunction(interp, path, "linear fixture").success);

		Persistence::SampledRealFunctionData loaded;
		REQUIRE(Persistence::LoadSampledFunctionData(path, loaded).success);
		REQUIRE(loaded.nodes.size() == 3);
		REQUIRE(loaded.values.size() == 3);
		REQUIRE(loaded.label == "linear fixture");
		REQUIRE(loaded.source == "Linear");
		for (int i = 0; i < 3; ++i)
		{
			REQUIRE(loaded.nodes[i] == x[i]);
			REQUIRE(loaded.values[i] == y[i]);
		}

		CleanupTempFile(path);
	}

	TEST_CASE("Persistence sampled function JSON rejects invalid save inputs", "[Persistence][FunctionJSON]") {
		RealFunction sinFunction([](const Real x) { return std::sin(x); });
		std::stringstream stream;

		REQUIRE_FALSE(Persistence::SaveSampledFunction(stream, sinFunction, "bad", Persistence::SampleGrid{0.0, 1.0, 1}).success);
		REQUIRE_FALSE(Persistence::SaveSampledFunction(stream, sinFunction, "bad", Persistence::SampleGrid{1.0, 0.0, 4}).success);

		Persistence::SampledRealFunctionData duplicate;
		duplicate.nodes = Vector<Real>(std::vector<Real>{0.0, 0.0});
		duplicate.values = Vector<Real>(std::vector<Real>{1.0, 2.0});
		REQUIRE_FALSE(Persistence::SaveSampledFunction(stream, duplicate).success);

		Persistence::SampledRealFunctionData mismatched;
		mismatched.nodes = Vector<Real>(std::vector<Real>{0.0, 1.0});
		mismatched.values = Vector<Real>(std::vector<Real>{1.0});
		REQUIRE_FALSE(Persistence::SaveSampledFunction(stream, mismatched).success);
	}

	TEST_CASE("Persistence sampled function JSON validates schema scalar and data shape", "[Persistence][FunctionJSON]") {
		const std::string wrongScalar = R"json({
			"mml": {
				"format": "MML_JSON",
				"version": 1,
				"object": "SampledRealFunction",
				"object_version": 1,
				"scalar": "float",
				"scalar_bytes": 8,
				"scalar_encoding": "ieee754",
				"real_type": "double"
			},
			"domain": [0, 1],
			"nodes": [0, 1],
			"values": [0, 1]
		})json";

		std::stringstream wrongScalarStream(wrongScalar);
		Persistence::SampledRealFunctionData scalarRejected;
		IOResult scalarResult = Persistence::LoadSampledFunctionData(wrongScalarStream, scalarRejected);
		REQUIRE_FALSE(scalarResult.success);
		REQUIRE(scalarResult.error == IOError::UNSUPPORTED_SCALAR);

		const std::string wrongShape = WithRealMetadata(R"json({
			"mml": {
				"format": "MML_JSON",
				"version": 1,
				"object": "SampledRealFunction",
				"object_version": 1,
				"scalar": "Real",
				"scalar_bytes": 8,
				"scalar_encoding": "ieee754",
				"real_type": "double"
			},
			"domain": [0, 1],
			"nodes": [0, 1],
			"values": [0]
		})json");

		std::stringstream wrongShapeStream(wrongShape);
		Persistence::SampledRealFunctionData shapeRejected;
		IOResult shapeResult = Persistence::LoadSampledFunctionData(wrongShapeStream, shapeRejected);
		REQUIRE_FALSE(shapeResult.success);
		REQUIRE(shapeResult.error == IOError::SCHEMA_MISMATCH);

		const std::string duplicateNodes = WithRealMetadata(R"json({
			"mml": {
				"format": "MML_JSON",
				"version": 1,
				"object": "SampledRealFunction",
				"object_version": 1,
				"scalar": "Real",
				"scalar_bytes": 8,
				"scalar_encoding": "ieee754",
				"real_type": "double"
			},
			"domain": [0, 1],
			"nodes": [0, 0],
			"values": [0, 1]
		})json");

		std::stringstream duplicateStream(duplicateNodes);
		Persistence::SampledRealFunctionData duplicateRejected;
		IOResult duplicateResult = Persistence::LoadSampledFunctionData(duplicateStream, duplicateRejected);
		REQUIRE_FALSE(duplicateResult.success);
		REQUIRE(duplicateResult.error == IOError::INVALID_PARAMETERS);
	}

	TEST_CASE("Persistence sampled function loads into selected interpolation types", "[Persistence][FunctionJSON]") {
		Persistence::SampledRealFunctionData data;
		data.nodes = Vector<Real>(std::vector<Real>{0.0, 1.0, 2.0, 3.0});
		data.values = Vector<Real>(std::vector<Real>{0.0, 1.0, 4.0, 9.0});
		data.label = "quadratic samples";

		std::stringstream stream;
		REQUIRE(Persistence::SaveSampledFunction(stream, data).success);

		auto linear = Persistence::LoadLinearFunction(stream);
		REQUIRE(linear.result.success);
		REQUIRE(std::abs(linear.function(1.5) - 2.5) < 1e-12);

		std::stringstream polyStream(stream.str());
		auto polynomial = Persistence::LoadPolynomialFunction(polyStream, Persistence::PolynomialInterpolationOptions{3});
		REQUIRE(polynomial.result.success);
		REQUIRE(std::abs(polynomial.function(1.5) - 2.25) < 1e-12);

		std::stringstream splineStream(stream.str());
		auto spline = Persistence::LoadSplineFunction(splineStream);
		REQUIRE(spline.result.success);
		REQUIRE(std::abs(spline.function(1.0) - 1.0) < 1e-12);
		REQUIRE(std::abs(spline.function(2.0) - 4.0) < 1e-12);
	}

	TEST_CASE("Persistence sampled function rejects bad interpolation options", "[Persistence][FunctionJSON]") {
		Persistence::SampledRealFunctionData data;
		data.nodes = Vector<Real>(std::vector<Real>{0.0, 1.0, 2.0});
		data.values = Vector<Real>(std::vector<Real>{0.0, 1.0, 4.0});

		std::stringstream stream;
		REQUIRE(Persistence::SaveSampledFunction(stream, data).success);

		auto tooSmall = Persistence::LoadPolynomialFunction(stream, Persistence::PolynomialInterpolationOptions{1});
		REQUIRE_FALSE(tooSmall.result.success);
		REQUIRE(tooSmall.result.error == IOError::INVALID_PARAMETERS);

		std::stringstream tooLargeStream(stream.str());
		auto tooLarge = Persistence::LoadPolynomialFunction(tooLargeStream, Persistence::PolynomialInterpolationOptions{4});
		REQUIRE_FALSE(tooLarge.result.success);
		REQUIRE(tooLarge.result.error == IOError::INVALID_PARAMETERS);
	}

	// Test function for Real -> Real
	class SinFunction : public IRealFunction {
	public:
		Real operator()(Real x) const override { return std::sin(x); }
	};

	class CosFunction : public IRealFunction {
	public:
		Real operator()(Real x) const override { return std::cos(x); }
	};

	// Test function for scalar field f(x,y) = x*y
	class ScalarProd2D : public IScalarFunction<2> {
	public:
		Real operator()(const VectorN<Real, 2>& p) const override {
			return p[0] * p[1];
		}
	};

	// Test function for scalar field f(x,y,z) = x + y + z
	class ScalarSum3D : public IScalarFunction<3> {
	public:
		Real operator()(const VectorN<Real, 3>& p) const override {
			return p[0] + p[1] + p[2];
		}
	};

	// Test 2D vector field: F(x,y) = (x, y)
	class IdentityVectorField2D : public IVectorFunction<2> {
	public:
		VectorN<Real, 2> operator()(const VectorN<Real, 2>& p) const override {
			return p;
		}
	};

	// Test 3D vector field: F(x,y,z) = (x, y, z)
	class IdentityVectorField3D : public IVectorFunction<3> {
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& p) const override {
			return p;
		}
	};

	// Test parametric curve: circle
	class CircleCurve2D : public IRealToVectorFunction<2> {
	public:
		VectorN<Real, 2> operator()(Real t) const override {
			return VectorN<Real, 2>{std::cos(t), std::sin(t)};
		}
	};

	// Test parametric curve: helix
	class HelixCurve3D : public IRealToVectorFunction<3> {
	public:
		VectorN<Real, 3> operator()(Real t) const override {
			return VectorN<Real, 3>{std::cos(t), std::sin(t), static_cast<Real>(t / (2 * Constants::PI))};
		}
	};

	// Test IVectorFunctionNM<2,3> for parametric surface
	class SphereSurfaceNM : public IVectorFunctionNM<2, 3> {
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 2>& p) const override {
			Real u = p[0];
			Real w = p[1];
			return VectorN<Real, 3>{
				std::cos(u) * std::sin(w),
				std::sin(u) * std::sin(w),
				std::cos(w)
			};
		}
	};

	/////////////////////////////////////////////////////////////////////////////////////
	///                         REAL FUNCTION SERIALIZATION                           ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Persistence - SaveRealFunc equally spaced", "[serializer][realfunc]") {
		TEST_PRECISION_INFO();
		
		std::string testFile = GetTempFilePath("realfunc_equally_spaced");
		
		SECTION("Successful serialization with IRealFunction") {
			SinFunction sinFunc;
			auto result = Serializer::SaveRealFuncEquallySpaced(sinFunc, "sin(x)", 0.0, Constants::PI, 10, testFile);
			
			REQUIRE(result.success == true);
			REQUIRE(result.error == IOError::OK);

			// Verify file exists and has content
			std::string content = ReadFileContents(testFile);
			REQUIRE(content.find("MML_REAL_FUNCTION") != std::string::npos);
			REQUIRE(content.find("sin(x)") != std::string::npos);
			REQUIRE(content.find("x1:") != std::string::npos);
			REQUIRE(content.find("NumPoints: 10") != std::string::npos);
			
			CleanupTempFile(testFile);
		}
	}

	TEST_CASE("Serializer - SaveRealFunc specified points", "[serializer][realfunc]") {
		TEST_PRECISION_INFO();
		
		std::string testFile = GetTempFilePath("realfunc_specified");
		
		SECTION("Successful serialization with custom points") {
			SinFunction sinFunc;
			Vector<Real> points(std::vector<Real>{0.0, 0.5, 1.0, 1.5, 2.0});
			auto result = Serializer::SaveRealFunc(sinFunc, "sin(x) custom", points, testFile);
			
			REQUIRE(result.success == true);
			
			std::string content = ReadFileContents(testFile);
			REQUIRE(content.find("MML_REAL_FUNCTION") != std::string::npos);
			// Should have 5 data points
			std::istringstream iss(content);
			std::string line;
			int dataLines = 0;
			while (std::getline(iss, line)) {
				// Count lines that look like data (start with number)
				if (!line.empty() && (std::isdigit(line[0]) || line[0] == '-')) {
					dataLines++;
				}
			}
			REQUIRE(dataLines == 5);
			
			CleanupTempFile(testFile);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         MULTI-FUNCTION SERIALIZATION                          ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Serializer - SaveRealMultiFunc", "[serializer][multifunc]") {
		TEST_PRECISION_INFO();
		
		std::string testFile = GetTempFilePath("multifunc");
		
		SECTION("Multiple IRealFunction pointers") {
			SinFunction sinFunc;
			CosFunction cosFunc;
			std::vector<IRealFunction*> funcs = {&sinFunc, &cosFunc};
			std::vector<std::string> legend = {"sin(x)", "cos(x)"};
			
			auto result = Serializer::SaveRealMultiFunc(funcs, "Trig Functions", legend, 
			                                             0.0, Constants::PI, 20, testFile);
			
			REQUIRE(result.success == true);
			
			std::string content = ReadFileContents(testFile);
			REQUIRE(content.find("MML_MULTI_REAL_FUNCTION") != std::string::npos);
			REQUIRE(content.find("Trig Functions") != std::string::npos);
			REQUIRE(content.find("sin(x)") != std::string::npos);
			REQUIRE(content.find("cos(x)") != std::string::npos);
			
			CleanupTempFile(testFile);
		}

		SECTION("Linear interpolation functions") {
			// Create simple linear interpolations
			Vector<Real> x(std::vector<Real>{0.0, 1.0, 2.0});
			Vector<Real> y1(std::vector<Real>{0.0, 1.0, 0.0});
			Vector<Real> y2(std::vector<Real>{1.0, 0.0, 1.0});
			
			LinearInterpRealFunc interp1(x, y1);
			LinearInterpRealFunc interp2(x, y2);
			
			std::vector<LinearInterpRealFunc> funcs = {interp1, interp2};
			std::vector<std::string> legend = {"Triangle1", "Triangle2"};
			
			auto result = Serializer::SaveRealMultiFunc(funcs, "Linear Interps", legend,
			                                             0.0, 2.0, 10, testFile);
			
			REQUIRE(result.success == true);
			
			std::string content = ReadFileContents(testFile);
			REQUIRE(content.find("MML_MULTI_REAL_FUNCTION") != std::string::npos);
			REQUIRE(content.find("Triangle1") != std::string::npos);
			
			CleanupTempFile(testFile);
		}

		SECTION("Polynomial interpolation functions") {
			Vector<Real> x(std::vector<Real>{0.0, 1.0, 2.0});
			Vector<Real> y1(std::vector<Real>{0.0, 1.0, 4.0});
			Vector<Real> y2(std::vector<Real>{1.0, 2.0, 5.0});

			PolynomInterpRealFunc interp1(x, y1, 3);
			PolynomInterpRealFunc interp2(x, y2, 3);

			std::vector<PolynomInterpRealFunc> funcs = {interp1, interp2};
			std::vector<std::string> legend = {"Poly1", "Poly2"};

			auto result = Serializer::SaveRealMultiFunc(funcs, "Polynomial Interps", legend,
			                                             0.0, 2.0, 10, testFile);

			REQUIRE(result.success == true);

			std::string content = ReadFileContents(testFile);
			REQUIRE(content.find("MML_MULTI_REAL_FUNCTION") != std::string::npos);
			REQUIRE(content.find("Poly1") != std::string::npos);

			CleanupTempFile(testFile);
		}

		SECTION("Spline interpolation functions") {
			Vector<Real> x(std::vector<Real>{0.0, 1.0, 2.0, 3.0});
			Vector<Real> y1(std::vector<Real>{0.0, 1.0, 0.0, -1.0});
			Vector<Real> y2(std::vector<Real>{1.0, 0.0, -1.0, 0.0});

			SplineInterpRealFunc interp1(x, y1);
			SplineInterpRealFunc interp2(x, y2);

			std::vector<SplineInterpRealFunc> funcs = {interp1, interp2};
			std::vector<std::string> legend = {"Spline1", "Spline2"};

			auto result = Serializer::SaveRealMultiFunc(funcs, "Spline Interps", legend,
			                                             0.0, 3.0, 10, testFile);

			REQUIRE(result.success == true);

			std::string content = ReadFileContents(testFile);
			REQUIRE(content.find("MML_MULTI_REAL_FUNCTION") != std::string::npos);
			REQUIRE(content.find("Spline1") != std::string::npos);

			CleanupTempFile(testFile);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         PARAMETRIC CURVE SERIALIZATION                        ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Serializer - SaveAsParamCurve2D", "[serializer][paramcurve]") {
		TEST_PRECISION_INFO();
		
		std::string testFile = GetTempFilePath("paramcurve2d");
		
		SECTION("From vectors") {
			Vector<Real> x(std::vector<Real>{0.0, 1.0, 2.0, 3.0, 4.0});
			Vector<Real> y(std::vector<Real>{0.0, 1.0, 0.0, -1.0, 0.0});
			
			auto result = Serializer::SaveAsParamCurve2D(x, y, "Wave Curve", testFile, 0.0, 4.0);
			
			REQUIRE(result.success == true);
			
			std::string content = ReadFileContents(testFile);
			REQUIRE(content.find("MML_PARAMETRIC_CURVE_CARTESIAN_2D") != std::string::npos);
			REQUIRE(content.find("Wave Curve") != std::string::npos);
			REQUIRE(content.find("t1:") != std::string::npos);
			REQUIRE(content.find("t2:") != std::string::npos);
			
			CleanupTempFile(testFile);
		}
	}

	TEST_CASE("Serializer - SaveParamCurveCartesian2D", "[serializer][paramcurve]") {
		TEST_PRECISION_INFO();
		
		std::string testFile = GetTempFilePath("paramcurve2d_func");
		
		SECTION("From IRealToVectorFunction<2>") {
			CircleCurve2D circle;
			
			bool result = Serializer::SaveParamCurveCartesian2DResult(circle, "Circle", 
			                                                     0.0, 2 * Constants::PI, 36, testFile).success;
			
			REQUIRE(result == true);
			
			std::string content = ReadFileContents(testFile);
			REQUIRE(content.find("MML_PARAMETRIC_CURVE_CARTESIAN_2D") != std::string::npos);
			
			CleanupTempFile(testFile);
		}
	}

	TEST_CASE("Serializer - SaveParamCurveCartesian3D", "[serializer][paramcurve]") {
		TEST_PRECISION_INFO();
		
		std::string testFile = GetTempFilePath("paramcurve3d_func");
		
		SECTION("From IRealToVectorFunction<3>") {
			HelixCurve3D helix;
			
			bool result = Serializer::SaveParamCurveCartesian3DResult(helix, "Helix",
			                                                     0.0, 4 * Constants::PI, 100, testFile).success;
			
			REQUIRE(result == true);
			
			std::string content = ReadFileContents(testFile);
			REQUIRE(content.find("MML_PARAMETRIC_CURVE_CARTESIAN_3D") != std::string::npos);
			
			CleanupTempFile(testFile);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         PARAMETRIC SURFACE SERIALIZATION                      ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Serializer - SaveParametricSurface", "[serializer][paramsurface]") {
		TEST_PRECISION_INFO();
		
		std::string testFile = GetTempFilePath("paramsurface");
		
		SECTION("Sphere surface via IVectorFunctionNM") {
			SphereSurfaceNM sphere;
			
			auto result = Serializer::SaveParametricSurface(sphere, "Unit Sphere",
			                                                 0.0, 2 * Constants::PI, 18,
			                                                 0.0, Constants::PI, 9,
			                                                 testFile);
			
			REQUIRE(result.success == true);
			
			std::string content = ReadFileContents(testFile);
			REQUIRE(content.find("MML_PARAMETRIC_SURFACE_CARTESIAN") != std::string::npos);
			REQUIRE(content.find("Unit Sphere") != std::string::npos);
			REQUIRE(content.find("u1:") != std::string::npos);
			REQUIRE(content.find("w1:") != std::string::npos);
			REQUIRE(content.find("NumPointsU:") != std::string::npos);
			REQUIRE(content.find("NumPointsW:") != std::string::npos);
			
			CleanupTempFile(testFile);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         SCALAR FUNCTION SERIALIZATION                         ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Serializer - SaveScalarFunc2DCartesian", "[serializer][scalarfunc]") {
		TEST_PRECISION_INFO();
		
		std::string testFile = GetTempFilePath("scalar2d");
		
		SECTION("2D scalar field") {
			ScalarProd2D scalarFunc;
			
			auto result = Serializer::SaveScalarFunc2DCartesian(scalarFunc, "z=x*y",
			                                                     0.0, 5.0, 10,
			                                                     0.0, 5.0, 10,
			                                                     testFile);
			
			REQUIRE(result.success == true);
			
			std::string content = ReadFileContents(testFile);
			REQUIRE(content.find("MML_SCALAR_FUNCTION_CARTESIAN_2D") != std::string::npos);
			REQUIRE(content.find("z=x*y") != std::string::npos);
			REQUIRE(content.find("NumPointsX:") != std::string::npos);
			REQUIRE(content.find("NumPointsY:") != std::string::npos);
			
			CleanupTempFile(testFile);
		}
	}

	TEST_CASE("Serializer - SaveScalarFunc3DCartesian", "[serializer][scalarfunc]") {
		TEST_PRECISION_INFO();
		
		std::string testFile = GetTempFilePath("scalar3d");
		
		SECTION("3D scalar field") {
			ScalarSum3D scalarFunc;
			
			auto result = Serializer::SaveScalarFunc3DCartesian(scalarFunc, "w=x+y+z",
			                                                     0.0, 2.0, 3,
			                                                     0.0, 2.0, 3,
			                                                     0.0, 2.0, 3,
			                                                     testFile);
			
			REQUIRE(result.success == true);
			
			std::string content = ReadFileContents(testFile);
			REQUIRE(content.find("MML_SCALAR_FUNCTION_CARTESIAN_3D") != std::string::npos);
			REQUIRE(content.find("NumPointsZ:") != std::string::npos);
			
			CleanupTempFile(testFile);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         VECTOR FUNCTION SERIALIZATION                         ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Serializer - SaveVectorFunc2DCartesian", "[serializer][vectorfunc]") {
		TEST_PRECISION_INFO();
		
		std::string testFile = GetTempFilePath("vector2d");
		
		SECTION("2D vector field") {
			IdentityVectorField2D vecFunc;
			
			auto result = Serializer::SaveVectorFunc2DCartesian(vecFunc, "Identity Field",
			                                                     -5.0, 5.0, 5,
			                                                     -5.0, 5.0, 5,
			                                                     testFile);
			
			REQUIRE(result.success == true);
			
			std::string content = ReadFileContents(testFile);
			REQUIRE(content.find("MML_VECTOR_FIELD_2D_CARTESIAN") != std::string::npos);
			REQUIRE(content.find("Identity Field") != std::string::npos);
			
			CleanupTempFile(testFile);
		}

		SECTION("2D vector field with threshold") {
			IdentityVectorField2D vecFunc;
			
			auto result = Serializer::SaveVectorFunc2DCartesian(vecFunc, "Thresholded Field",
			                                                     -5.0, 5.0, 10,
			                                                     -5.0, 5.0, 10,
			                                                     testFile, 3.0);  // threshold = 3.0
			
			REQUIRE(result.success == true);
			
			// File should have fewer points due to threshold
			std::string content = ReadFileContents(testFile);
			REQUIRE(content.find("MML_VECTOR_FIELD_2D_CARTESIAN") != std::string::npos);
			
			CleanupTempFile(testFile);
		}
	}

	TEST_CASE("Serializer - SaveVectorFunc3DCartesian", "[serializer][vectorfunc]") {
		TEST_PRECISION_INFO();
		
		std::string testFile = GetTempFilePath("vector3d");
		
		SECTION("3D vector field") {
			IdentityVectorField3D vecFunc;
			
			auto result = Serializer::SaveVectorFunc3DCartesian(vecFunc, "3D Identity",
			                                                     -2.0, 2.0, 3,
			                                                     -2.0, 2.0, 3,
			                                                     -2.0, 2.0, 3,
			                                                     testFile);
			
			REQUIRE(result.success == true);
			
			std::string content = ReadFileContents(testFile);
			REQUIRE(content.find("MML_VECTOR_FIELD_3D_CARTESIAN") != std::string::npos);
			
			CleanupTempFile(testFile);
		}

		SECTION("3D vector field with threshold") {
			IdentityVectorField3D vecFunc;
			
			auto result = Serializer::SaveVectorFunc3DCartesian(vecFunc, "3D Thresholded",
			                                                     -2.0, 2.0, 4,
			                                                     -2.0, 2.0, 4,
			                                                     -2.0, 2.0, 4,
			                                                     testFile, 2.0);
			
			REQUIRE(result.success == true);
			
			CleanupTempFile(testFile);
		}
	}

	TEST_CASE("Serializer - SaveVectorFuncSpherical", "[serializer][vectorfunc]") {
		TEST_PRECISION_INFO();
		
		std::string testFile = GetTempFilePath("vectorspherical");
		
		SECTION("Spherical coordinate vector field") {
			IdentityVectorField3D vecFunc;
			
			auto result = Serializer::SaveVectorFuncSpherical(vecFunc, "Spherical Field",
			                                                   1.0, 5.0, 3,     // r
			                                                   0.0, Constants::PI, 3,  // theta
			                                                   0.0, 2 * Constants::PI, 4,  // phi
			                                                   testFile);
			
			REQUIRE(result.success == true);
			
			std::string content = ReadFileContents(testFile);
			REQUIRE(content.find("MML_VECTOR_FIELD_SPHERICAL") != std::string::npos);
			
			CleanupTempFile(testFile);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         PARTICLE SIMULATION SERIALIZATION                     ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Serializer - SaveParticleSimulation2D", "[serializer][particle]") {
		TEST_PRECISION_INFO();
		
		std::string testFile = GetTempFilePath("particle2d");
		
		SECTION("2D particle simulation") {
			// Create simple 2-ball simulation
			int numBalls = 2;
			Real width = 100.0, height = 100.0;
			Real dT = 0.01;
			
			std::vector<std::vector<Pnt2Cart>> positions(numBalls);
			for (int b = 0; b < numBalls; ++b) {
				for (int step = 0; step < 10; ++step) {
					positions[b].push_back(Pnt2Cart(b * 10.0 + step, b * 5.0 + step));
				}
			}
			
			std::vector<std::string> colors = {"red", "blue"};
			std::vector<Real> radii = {1.0, 2.0};
			
			auto result = Serializer::SaveParticleSimulation2D(testFile, numBalls, 
			                                                    width, height,
			                                                    positions, colors, radii,
			                                                    dT, 1);
			
			REQUIRE(result.success == true);
			
			std::string content = ReadFileContents(testFile);
			REQUIRE(content.find("MML_PARTICLE_SIMULATION_DATA_2D") != std::string::npos);
			REQUIRE(content.find("Width:") != std::string::npos);
			REQUIRE(content.find("Height:") != std::string::npos);
			REQUIRE(content.find("NumBalls:") != std::string::npos);
			REQUIRE(content.find("Ball_1") != std::string::npos);
			REQUIRE(content.find("red") != std::string::npos);
			REQUIRE(content.find("Step 0") != std::string::npos);
			
			CleanupTempFile(testFile);
		}

		SECTION("2D particle with saveEveryNSteps") {
			int numBalls = 1;
			std::vector<std::vector<Pnt2Cart>> positions(numBalls);
			for (int step = 0; step < 20; ++step) {
				positions[0].push_back(Pnt2Cart(static_cast<Real>(step), static_cast<Real>(step)));
			}
			
			std::vector<std::string> colors = {"green"};
			std::vector<Real> radii = {1.0};
			
			// Save every 5th step
			auto result = Serializer::SaveParticleSimulation2D(testFile, numBalls,
			                                                    100.0, 100.0,
			                                                    positions, colors, radii,
			                                                    0.01, 5);
			
			REQUIRE(result.success == true);
			
			std::string content = ReadFileContents(testFile);
			// Should have 4 steps (0, 5, 10, 15)
			REQUIRE(content.find("NumSteps: 4") != std::string::npos);
			
			CleanupTempFile(testFile);
		}
	}

	TEST_CASE("Serializer - SaveParticleSimulation3D", "[serializer][particle]") {
		TEST_PRECISION_INFO();
		
		std::string testFile = GetTempFilePath("particle3d");
		
		SECTION("3D particle simulation") {
			int numBalls = 2;
			std::vector<std::vector<Pnt3Cart>> positions(numBalls);
			for (int b = 0; b < numBalls; ++b) {
				for (int step = 0; step < 5; ++step) {
					positions[b].push_back(Pnt3Cart(b * 10.0, step * 1.0, static_cast<Real>(b + step)));
				}
			}
			
			std::vector<std::string> colors = {"white", "black"};
			std::vector<Real> radii = {0.5, 0.75};
			
			auto result = Serializer::SaveParticleSimulation3D(testFile, numBalls,
			                                                    50.0, 50.0, 50.0,
			                                                    positions, colors, radii,
			                                                    0.01, 1);
			
			REQUIRE(result.success == true);
			
			std::string content = ReadFileContents(testFile);
			REQUIRE(content.find("MML_PARTICLE_SIMULATION_DATA_3D") != std::string::npos);
			REQUIRE(content.find("Depth:") != std::string::npos);
			
			CleanupTempFile(testFile);
		}

		SECTION("3D particle with saveEveryNSteps") {
			int numBalls = 1;
			std::vector<std::vector<Pnt3Cart>> positions(numBalls);
			for (int step = 0; step < 20; ++step) {
				positions[0].push_back(Pnt3Cart(static_cast<Real>(step), static_cast<Real>(step), static_cast<Real>(step)));
			}

			std::vector<std::string> colors = {"green"};
			std::vector<Real> radii = {1.0};

			// Save every 5th step
			auto result = Serializer::SaveParticleSimulation3D(testFile, numBalls,
			                                                    100.0, 100.0, 100.0,
			                                                    positions, colors, radii,
			                                                    0.01, 5);

			REQUIRE(result.success == true);

			std::string content = ReadFileContents(testFile);
			// Should have 4 steps (0, 5, 10, 15)
			REQUIRE(content.find("NumSteps: 4") != std::string::npos);

			CleanupTempFile(testFile);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         ERROR HANDLING                                        ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Serializer - Error handling", "[serializer][error]") {
		TEST_PRECISION_INFO();

		SECTION("Empty filename") {
			SinFunction sinFunc;
			auto result = Serializer::SaveRealFuncEquallySpaced(sinFunc, "test", 0, 1, 10, "");
			
			REQUIRE(result.success == false);
			REQUIRE(result.error == IOError::INVALID_PARAMETERS);
		}

		SECTION("Invalid range (x1 >= x2)") {
			SinFunction sinFunc;
			std::string testFile = GetTempFilePath("error_range");
			
			auto result = Serializer::SaveRealFuncEquallySpaced(sinFunc, "test", 5.0, 1.0, 10, testFile);
			
			REQUIRE(result.success == false);
			REQUIRE(result.error == IOError::INVALID_PARAMETERS);
		}

		SECTION("Invalid numPoints (< 2)") {
			SinFunction sinFunc;
			std::string testFile = GetTempFilePath("error_numpoints");
			
			auto result = Serializer::SaveRealFuncEquallySpaced(sinFunc, "test", 0.0, 1.0, 1, testFile);
			
			REQUIRE(result.success == false);
			REQUIRE(result.error == IOError::INVALID_PARAMETERS);
		}

		SECTION("Mismatched funcs and legend sizes") {
			SinFunction sinFunc;
			std::vector<IRealFunction*> funcs = {&sinFunc};
			std::vector<std::string> legend = {"func1", "func2"};  // Too many legends
			std::string testFile = GetTempFilePath("error_mismatch");
			
			auto result = Serializer::SaveRealMultiFunc(funcs, "test", legend, 0.0, 1.0, 10, testFile);
			
			REQUIRE(result.success == false);
			REQUIRE(result.error == IOError::INVALID_PARAMETERS);
		}

		SECTION("Empty funcs vector") {
			std::vector<IRealFunction*> funcs;
			std::vector<std::string> legend;
			std::string testFile = GetTempFilePath("error_empty");
			
			auto result = Serializer::SaveRealMultiFunc(funcs, "test", legend, 0.0, 1.0, 10, testFile);
			
			REQUIRE(result.success == false);
			REQUIRE(result.error == IOError::INVALID_PARAMETERS);
		}

		SECTION("Invalid particle simulation parameters") {
			std::string testFile = GetTempFilePath("error_particle");
			
			// numBalls <= 0
			auto result = Serializer::SaveParticleSimulation2D(testFile, 0, 100, 100,
			                                                    {}, {}, {}, 0.01, 1);
			
			REQUIRE(result.success == false);
			REQUIRE(result.error == IOError::INVALID_PARAMETERS);
		}

		SECTION("File cannot be opened (invalid path)") {
			SinFunction sinFunc;
			// Try to write to a directory that doesn't exist
			std::string invalidPath = "/nonexistent_dir_12345/test.mml";
			
			auto result = Serializer::SaveRealFuncEquallySpaced(sinFunc, "test", 0.0, 1.0, 10, invalidPath);
			
			REQUIRE(result.success == false);
			REQUIRE(result.error == IOError::FILE_NOT_OPENED);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                         FORMAT CONSISTENCY                                    ///
	/////////////////////////////////////////////////////////////////////////////////////

	TEST_CASE("Serializer - Format consistency", "[serializer][format]") {
		TEST_PRECISION_INFO();

		SECTION("REAL_FUNCTION format structure") {
			std::string testFile = GetTempFilePath("format_real");
			SinFunction sinFunc;

			REQUIRE(Serializer::SaveRealFuncEquallySpaced(sinFunc, "test_title", 0.0, 1.0, 5, testFile).success);

			std::string content = ReadFileContents(testFile);
			std::istringstream iss(content);
			std::string line;
			
			// Line 1: Type
			std::getline(iss, line);
			REQUIRE(line == "MML_REAL_FUNCTION_EQUALLY_SPACED");
			
			// Line 2: Version
			std::getline(iss, line);
			REQUIRE(line == "VERSION: 1");
			
			// Line 3: Title
			std::getline(iss, line);
			REQUIRE(line == "test_title");
			
			// Line 4: x1
			std::getline(iss, line);
			REQUIRE(line.find("x1:") == 0);
			
			// Line 5: x2
			std::getline(iss, line);
			REQUIRE(line.find("x2:") == 0);
			
			// Line 6: NumPoints
			std::getline(iss, line);
			REQUIRE(line.find("NumPoints:") == 0);
			
			CleanupTempFile(testFile);
		}

		SECTION("MULTI_REAL_FUNCTION format structure") {
			std::string testFile = GetTempFilePath("format_multi");
			SinFunction sinFunc;
			CosFunction cosFunc;
			std::vector<IRealFunction*> funcs = {&sinFunc, &cosFunc};
			std::vector<std::string> legend = {"sin", "cos"};

			REQUIRE(Serializer::SaveRealMultiFunc(funcs, "multi_test", legend, 0.0, 1.0, 3, testFile).success);

			std::string content = ReadFileContents(testFile);
			std::istringstream iss(content);
			std::string line;
			
			// Line 1: Type
			std::getline(iss, line);
			REQUIRE(line == "MML_MULTI_REAL_FUNCTION");
			
			// Line 2: Version
			std::getline(iss, line);
			REQUIRE(line == "VERSION: 1");
			
			// Line 3: Title
			std::getline(iss, line);
			REQUIRE(line == "multi_test");
			
			// Line 4: Number of functions
			std::getline(iss, line);
			REQUIRE(line == "2");
			
			// Lines 5-6: Legend entries
			std::getline(iss, line);
			REQUIRE(line == "sin");
			std::getline(iss, line);
			REQUIRE(line == "cos");
			
			CleanupTempFile(testFile);
		}

		SECTION("Data rows have correct column count") {
			std::string testFile = GetTempFilePath("format_columns");
			SinFunction sinFunc;
			CosFunction cosFunc;
			std::vector<IRealFunction*> funcs = {&sinFunc, &cosFunc};
			std::vector<std::string> legend = {"sin", "cos"};

			REQUIRE(Serializer::SaveRealMultiFunc(funcs, "test", legend, 0.0, 1.0, 5, testFile).success);

			std::string content = ReadFileContents(testFile);
			std::istringstream iss(content);
			std::string line;
			
			// Skip header (9 lines: type, version, title, numFuncs, 2 legends, x1, x2, NumPoints)
			for (int i = 0; i < 9; ++i) {
				std::getline(iss, line);
			}
			
			// Data rows should have 3 columns: x, sin(x), cos(x)
			int dataRows = 0;
			while (std::getline(iss, line) && !line.empty()) {
				std::istringstream lineStream(line);
				Real val;
				int columns = 0;
				while (lineStream >> val) {
					columns++;
				}
				REQUIRE(columns == 3);  // x + 2 functions
				dataRows++;
			}
			REQUIRE(dataRows == 5);
			
			CleanupTempFile(testFile);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	///                    GOLDEN FILE CONTENT VALIDATION TESTS                       ///
	/////////////////////////////////////////////////////////////////////////////////////
	
	// These tests compare generated output against reference files to catch
	// both format changes AND numerical bugs. Each format has simple and complex
	// variants to exercise different code paths.
	
	TEST_CASE("Golden - REAL_FUNCTION_EQUALLY_SPACED", "[serializer][golden][realfunc]") {
		TEST_PRECISION_INFO();
		
		SinFunction sinFunc;
		
		SECTION("Simple - 5 points over half period") {
			std::string testFile = GetTempFilePath("golden_realfunc_es_simple");
			auto result = Serializer::SaveRealFuncEquallySpaced(sinFunc, "sin(x)",
			    0.0, Constants::PI, 5, testFile);
			REQUIRE(result.success);
			
			if constexpr (!std::is_same_v<Real, float>) {
				auto cmp = CompareFilesWithTolerance(testFile, 
				    GetReferenceFilePath("realfunc_equally_spaced_simple.mml"));
				INFO("Line " << cmp.line_number << ": " << cmp.error_message);
				REQUIRE(cmp.success);
			}
			
			CleanupTempFile(testFile);
		}
		
		SECTION("Complex - 13 points over full period") {
			std::string testFile = GetTempFilePath("golden_realfunc_es_complex");
			auto result = Serializer::SaveRealFuncEquallySpaced(sinFunc, "sin(x)",
			    0.0, 2*Constants::PI, 13, testFile);
			REQUIRE(result.success);
			
			if constexpr (!std::is_same_v<Real, float>) {
				auto cmp = CompareFilesWithTolerance(testFile,
				    GetReferenceFilePath("realfunc_equally_spaced_complex.mml"));
				INFO("Line " << cmp.line_number << ": " << cmp.error_message);
				REQUIRE(cmp.success);
			}
			
			CleanupTempFile(testFile);
		}
	}

	TEST_CASE("Golden - REAL_FUNCTION specified points", "[serializer][golden][realfunc]") {
		TEST_PRECISION_INFO();
		
		CosFunction cosFunc;
		
		SECTION("Simple - 4 uniform points") {
			std::string testFile = GetTempFilePath("golden_realfunc_sp_simple");
			Vector<Real> pts(std::vector<Real>{0.0, 1.0, 2.0, 3.0});
			auto result = Serializer::SaveRealFunc(cosFunc, "cos(x) custom", pts, testFile);
			REQUIRE(result.success);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("realfunc_specified_simple.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
		
		SECTION("Complex - 8 non-uniform points") {
			std::string testFile = GetTempFilePath("golden_realfunc_sp_complex");
			Vector<Real> pts(std::vector<Real>{0.0, 0.5, 1.0, 1.5, 2.5, 3.5, 5.0, 6.28});
			auto result = Serializer::SaveRealFunc(cosFunc, "cos(x) custom", pts, testFile);
			REQUIRE(result.success);
			
			if constexpr (!std::is_same_v<Real, float>) {
				auto cmp = CompareFilesWithTolerance(testFile,
				    GetReferenceFilePath("realfunc_specified_complex.mml"));
				INFO("Line " << cmp.line_number << ": " << cmp.error_message);
				REQUIRE(cmp.success);
			}
			
			CleanupTempFile(testFile);
		}
	}

	TEST_CASE("Golden - MULTI_REAL_FUNCTION", "[serializer][golden][multifunc]") {
		TEST_PRECISION_INFO();
		
		SinFunction sinFunc;
		CosFunction cosFunc;
		std::vector<IRealFunction*> funcs = {&sinFunc, &cosFunc};
		std::vector<std::string> legend = {"sin(x)", "cos(x)"};
		
		SECTION("Simple - 5 points") {
			std::string testFile = GetTempFilePath("golden_multi_simple");
			auto result = Serializer::SaveRealMultiFunc(funcs, "Trig Functions", legend,
			    0.0, Constants::PI, 5, testFile);
			REQUIRE(result.success);
			
			if constexpr (!std::is_same_v<Real, float>) {
				auto cmp = CompareFilesWithTolerance(testFile,
				    GetReferenceFilePath("multi_realfunc_simple.mml"));
				INFO("Line " << cmp.line_number << ": " << cmp.error_message);
				REQUIRE(cmp.success);
			}
			
			CleanupTempFile(testFile);
		}
		
		SECTION("Complex - 11 points over full period") {
			std::string testFile = GetTempFilePath("golden_multi_complex");
			auto result = Serializer::SaveRealMultiFunc(funcs, "Trig Functions", legend,
			    0.0, 2*Constants::PI, 11, testFile);
			REQUIRE(result.success);
			
			if constexpr (!std::is_same_v<Real, float>) {
				auto cmp = CompareFilesWithTolerance(testFile,
				    GetReferenceFilePath("multi_realfunc_complex.mml"));
				INFO("Line " << cmp.line_number << ": " << cmp.error_message);
				REQUIRE(cmp.success);
			}
			
			CleanupTempFile(testFile);
		}
	}

	TEST_CASE("Golden - PARAMETRIC_CURVE_2D", "[serializer][golden][paramcurve]") {
		TEST_PRECISION_INFO();
		
		CircleCurve2D circle;
		
		SECTION("Simple - 9 points") {
			std::string testFile = GetTempFilePath("golden_pc2d_simple");
			bool result = Serializer::SaveParamCurveCartesian2DResult(circle, "Unit Circle",
			    0.0, 2*Constants::PI, 9, testFile).success;
			REQUIRE(result);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("paramcurve2d_simple.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
		
		SECTION("Complex - 25 points") {
			std::string testFile = GetTempFilePath("golden_pc2d_complex");
			bool result = Serializer::SaveParamCurveCartesian2DResult(circle, "Unit Circle",
			    0.0, 2*Constants::PI, 25, testFile).success;
			REQUIRE(result);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("paramcurve2d_complex.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
	}

	TEST_CASE("Golden - PARAMETRIC_CURVE_3D", "[serializer][golden][paramcurve]") {
		TEST_PRECISION_INFO();
		
		HelixCurve3D helix;
		
		SECTION("Simple - 9 points, one turn") {
			std::string testFile = GetTempFilePath("golden_pc3d_simple");
			bool result = Serializer::SaveParamCurveCartesian3DResult(helix, "Helix",
			    0.0, 2*Constants::PI, 9, testFile).success;
			REQUIRE(result);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("paramcurve3d_simple.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
		
		SECTION("Complex - 25 points, two turns") {
			std::string testFile = GetTempFilePath("golden_pc3d_complex");
			bool result = Serializer::SaveParamCurveCartesian3DResult(helix, "Helix",
			    0.0, 4*Constants::PI, 25, testFile).success;
			REQUIRE(result);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("paramcurve3d_complex.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
	}

	TEST_CASE("Golden - PARAMETRIC_SURFACE", "[serializer][golden][paramsurface]") {
		TEST_PRECISION_INFO();
		
		SphereSurfaceNM sphere;
		
		SECTION("Simple - 5x5 grid") {
			std::string testFile = GetTempFilePath("golden_ps_simple");
			auto result = Serializer::SaveParametricSurface(sphere, "Unit Sphere",
			    0.0, 2*Constants::PI, 5,
			    0.0, Constants::PI, 5,
			    testFile);
			REQUIRE(result.success);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("paramsurface_simple.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
		
		SECTION("Complex - 9x9 grid") {
			std::string testFile = GetTempFilePath("golden_ps_complex");
			auto result = Serializer::SaveParametricSurface(sphere, "Unit Sphere",
			    0.0, 2*Constants::PI, 9,
			    0.0, Constants::PI, 9,
			    testFile);
			REQUIRE(result.success);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("paramsurface_complex.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
	}

	TEST_CASE("Golden - SCALAR_FUNCTION_2D", "[serializer][golden][scalarfunc]") {
		TEST_PRECISION_INFO();
		
		ScalarProd2D scalarFunc;
		
		SECTION("Simple - 4x4 grid") {
			std::string testFile = GetTempFilePath("golden_sf2d_simple");
			auto result = Serializer::SaveScalarFunc2DCartesian(scalarFunc, "z=x*y",
			    0.0, 3.0, 4,
			    0.0, 3.0, 4,
			    testFile);
			REQUIRE(result.success);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("scalar2d_simple.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
		
		SECTION("Complex - 7x7 grid with negative values") {
			std::string testFile = GetTempFilePath("golden_sf2d_complex");
			auto result = Serializer::SaveScalarFunc2DCartesian(scalarFunc, "z=x*y",
			    -2.0, 4.0, 7,
			    -2.0, 4.0, 7,
			    testFile);
			REQUIRE(result.success);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("scalar2d_complex.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
	}

	TEST_CASE("Golden - SCALAR_FUNCTION_3D", "[serializer][golden][scalarfunc]") {
		TEST_PRECISION_INFO();
		
		ScalarSum3D scalarFunc;
		
		SECTION("Simple - 3x3x3 grid") {
			std::string testFile = GetTempFilePath("golden_sf3d_simple");
			auto result = Serializer::SaveScalarFunc3DCartesian(scalarFunc, "w=x+y+z",
			    0.0, 2.0, 3,
			    0.0, 2.0, 3,
			    0.0, 2.0, 3,
			    testFile);
			REQUIRE(result.success);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("scalar3d_simple.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
		
		SECTION("Complex - 4x4x4 grid with negative values") {
			std::string testFile = GetTempFilePath("golden_sf3d_complex");
			auto result = Serializer::SaveScalarFunc3DCartesian(scalarFunc, "w=x+y+z",
			    -1.0, 2.0, 4,
			    -1.0, 2.0, 4,
			    -1.0, 2.0, 4,
			    testFile);
			REQUIRE(result.success);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("scalar3d_complex.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
	}

	TEST_CASE("Golden - VECTOR_FIELD_2D", "[serializer][golden][vectorfunc]") {
		TEST_PRECISION_INFO();
		
		IdentityVectorField2D vecFunc;
		
		SECTION("Simple - 3x3 grid") {
			std::string testFile = GetTempFilePath("golden_vf2d_simple");
			auto result = Serializer::SaveVectorFunc2DCartesian(vecFunc, "Identity 2D",
			    -1.0, 1.0, 3,
			    -1.0, 1.0, 3,
			    testFile);
			REQUIRE(result.success);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("vector2d_simple.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
		
		SECTION("Complex - 5x5 grid") {
			std::string testFile = GetTempFilePath("golden_vf2d_complex");
			auto result = Serializer::SaveVectorFunc2DCartesian(vecFunc, "Identity 2D",
			    -2.0, 2.0, 5,
			    -2.0, 2.0, 5,
			    testFile);
			REQUIRE(result.success);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("vector2d_complex.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
	}

	TEST_CASE("Golden - VECTOR_FIELD_3D", "[serializer][golden][vectorfunc]") {
		TEST_PRECISION_INFO();
		
		IdentityVectorField3D vecFunc;
		
		SECTION("Simple - 2x2x2 grid") {
			std::string testFile = GetTempFilePath("golden_vf3d_simple");
			auto result = Serializer::SaveVectorFunc3DCartesian(vecFunc, "Identity 3D",
			    -1.0, 1.0, 2,
			    -1.0, 1.0, 2,
			    -1.0, 1.0, 2,
			    testFile);
			REQUIRE(result.success);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("vector3d_simple.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
		
		SECTION("Complex - 3x3x3 grid") {
			std::string testFile = GetTempFilePath("golden_vf3d_complex");
			auto result = Serializer::SaveVectorFunc3DCartesian(vecFunc, "Identity 3D",
			    -2.0, 2.0, 3,
			    -2.0, 2.0, 3,
			    -2.0, 2.0, 3,
			    testFile);
			REQUIRE(result.success);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("vector3d_complex.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
	}

	TEST_CASE("Golden - VECTOR_FIELD_SPHERICAL", "[serializer][golden][vectorfunc]") {
		TEST_PRECISION_INFO();
		
		IdentityVectorField3D vecFunc;
		
		SECTION("Simple - 2x3x4 grid") {
			std::string testFile = GetTempFilePath("golden_vfs_simple");
			auto result = Serializer::SaveVectorFuncSpherical(vecFunc, "Spherical Field",
			    1.0, 2.0, 2,
			    0.0, Constants::PI, 3,
			    0.0, 2*Constants::PI, 4,
			    testFile);
			REQUIRE(result.success);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("vectorspherical_simple.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
		
		SECTION("Complex - 3x5x6 grid") {
			std::string testFile = GetTempFilePath("golden_vfs_complex");
			auto result = Serializer::SaveVectorFuncSpherical(vecFunc, "Spherical Field",
			    0.5, 2.5, 3,
			    0.0, Constants::PI, 5,
			    0.0, 2*Constants::PI, 6,
			    testFile);
			REQUIRE(result.success);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("vectorspherical_complex.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
	}

	TEST_CASE("Golden - PARTICLE_SIMULATION_2D", "[serializer][golden][particle]") {
		TEST_PRECISION_INFO();
		
		SECTION("Simple - 2 balls, 5 steps") {
			std::string testFile = GetTempFilePath("golden_p2d_simple");
			
			int numBalls = 2;
			std::vector<std::vector<Pnt2Cart>> positions(numBalls);
			for (int b = 0; b < numBalls; ++b) {
				for (int step = 0; step < 5; ++step) {
					positions[b].push_back(Pnt2Cart(10.0 + b * 20.0 + step * 2.0,
					                                 50.0 + b * 10.0 + step * 3.0));
				}
			}
			std::vector<std::string> colors = {"red", "blue"};
			std::vector<Real> radii = {5.0, 7.5};
			
			auto result = Serializer::SaveParticleSimulation2D(testFile, numBalls,
			    100.0, 100.0, positions, colors, radii, 0.01, 1);
			REQUIRE(result.success);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("particle2d_simple.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
		
		SECTION("Complex - 3 balls, 10 steps, circular motion") {
			std::string testFile = GetTempFilePath("golden_p2d_complex");
			
			int numBalls = 3;
			std::vector<std::vector<Pnt2Cart>> positions(numBalls);
			for (int b = 0; b < numBalls; ++b) {
				for (int step = 0; step < 10; ++step) {
					Real angle = step * 0.3 + b * Constants::PI / 3;
					positions[b].push_back(Pnt2Cart(50.0 + 20.0 * std::cos(angle),
					                                 50.0 + 20.0 * std::sin(angle)));
				}
			}
			std::vector<std::string> colors = {"green", "yellow", "purple"};
			std::vector<Real> radii = {3.0, 4.0, 5.0};
			
			auto result = Serializer::SaveParticleSimulation2D(testFile, numBalls,
			    100.0, 100.0, positions, colors, radii, 0.02, 1);
			REQUIRE(result.success);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("particle2d_complex.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
	}

	TEST_CASE("Golden - PARTICLE_SIMULATION_3D", "[serializer][golden][particle]") {
		TEST_PRECISION_INFO();
		
		SECTION("Simple - 2 balls, 4 steps") {
			std::string testFile = GetTempFilePath("golden_p3d_simple");
			
			int numBalls = 2;
			std::vector<std::vector<Pnt3Cart>> positions(numBalls);
			for (int b = 0; b < numBalls; ++b) {
				for (int step = 0; step < 4; ++step) {
					positions[b].push_back(Pnt3Cart(10.0 + b * 15.0,
					                                 20.0 + step * 5.0,
					                                 30.0 + b * step));
				}
			}
			std::vector<std::string> colors = {"white", "black"};
			std::vector<Real> radii = {2.0, 3.0};
			
			auto result = Serializer::SaveParticleSimulation3D(testFile, numBalls,
			    50.0, 50.0, 50.0, positions, colors, radii, 0.01, 1);
			REQUIRE(result.success);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("particle3d_simple.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
		
		SECTION("Complex - 3 balls, 8 steps, helical motion") {
			std::string testFile = GetTempFilePath("golden_p3d_complex");
			
			int numBalls = 3;
			std::vector<std::vector<Pnt3Cart>> positions(numBalls);
			for (int b = 0; b < numBalls; ++b) {
				for (int step = 0; step < 8; ++step) {
					Real t = step * 0.5;
					positions[b].push_back(Pnt3Cart(
					    25.0 + 10.0 * std::cos(t + b * 2.0),
					    25.0 + 10.0 * std::sin(t + b * 2.0),
					    25.0 + 5.0 * t
					));
				}
			}
			std::vector<std::string> colors = {"cyan", "magenta", "orange"};
			std::vector<Real> radii = {1.5, 2.0, 2.5};
			
			auto result = Serializer::SaveParticleSimulation3D(testFile, numBalls,
			    50.0, 50.0, 75.0, positions, colors, radii, 0.02, 1);
			REQUIRE(result.success);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("particle3d_complex.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////
	// Locale independence and honest failure reporting (MinimalMathLibrary-wewi.1.4)
	/////////////////////////////////////////////////////////////////////////////////////

	namespace {
		// A locale whose decimal point is ',' and which groups thousands with '.', built from a
		// facet so the test does not depend on any system locale being installed.
		struct CommaDecimalPunct : std::numpunct<char> {
			char do_decimal_point() const override { return ','; }
			char do_thousands_sep() const override { return '.'; }
			std::string do_grouping() const override { return "\3"; }
		};

		struct GlobalLocaleGuard {
			std::locale previous;
			explicit GlobalLocaleGuard(std::locale selected = std::locale(std::locale::classic(), new CommaDecimalPunct))
					: previous(std::locale::global(selected)) {}
			~GlobalLocaleGuard() { std::locale::global(previous); }
		};

		class LocaleProbeFunction : public IRealFunction {
		public:
			Real operator()(Real x) const override { return static_cast<MML::Real>(1234.5) * x + static_cast<MML::Real>(0.1); }
		};

		template<int N>
		class LocaleProbeCurve : public IRealToVectorFunction<N> {
		public:
			Real value;
			explicit LocaleProbeCurve(Real v)
					: value(v) {}
			VectorN<Real, N> operator()(Real) const override {
				VectorN<Real, N> result;
				result[0] = value;
				return result;
			}
		};

		class LocaleProbeSurface : public IParametricSurface<3>, public IParametricSurfaceRect<3> {
		public:
			Real value;
			explicit LocaleProbeSurface(Real v)
					: value(v) {}
			Real getMinU() const override { return 0; }
			Real getMaxU() const override { return 1; }
			Real getMinW() const override { return 0; }
			Real getMaxW() const override { return 1; }
			Real getMinW(Real) const override { return 0; }
			Real getMaxW(Real) const override { return 1; }
			VectorN<Real, 3> operator()(Real, Real) const override { return {value, 0, 0}; }
			VectorN<Real, 3> operator()(const VectorN<Real, 2>&) const override { return {value, 0, 0}; }
		};

		struct LocaleProbeGrid {
			Real value;
			const LocaleProbeGrid& grid() const { return *this; }
			const LocaleProbeGrid& domain() const { return *this; }
			int numNodesX() const { return 2; }
			int numNodesY() const { return 2; }
			Real xMin() const { return 0; }
			Real xMax() const { return 1; }
			Real yMin() const { return 0; }
			Real yMax() const { return 1; }
			Real x(int i) const { return Real(i); }
			Real y(int i) const { return Real(i); }
			Real operator()(int, int) const { return value; }
		};
	} // namespace

	TEST_CASE("Every visualization serializer preserves Real payloads under a decimal comma locale",
						"[serializer][locale][round-trip][wewi.2.5]") {
		const Real probe = std::nextafter(Real(1234.5), std::numeric_limits<Real>::infinity());
		const int precision = std::numeric_limits<Real>::max_digits10;
		RealFunctionFromStdFunc function([&](Real) { return probe; });
		ScalarFunctionFromStdFunc<2> scalar2([&](const VectorN<Real, 2>&) { return probe; });
		ScalarFunctionFromStdFunc<3> scalar3([&](const VectorN<Real, 3>&) { return probe; });
		VectorFunctionFromStdFunc<2> vector2([&](const VectorN<Real, 2>&) { return VectorN<Real, 2>{probe, 0}; });
		VectorFunctionFromStdFunc<3> vector3([&](const VectorN<Real, 3>&) { return VectorN<Real, 3>{probe, 0, 0}; });
		LocaleProbeCurve<2> curve2(probe);
		LocaleProbeCurve<3> curve3(probe);
		LocaleProbeSurface surface(probe);
		LocaleProbeGrid grid{probe};
		const Vector<Real> parameters{0, 1};
		const std::vector<VectorN<Real, 2>> points2{{probe, 0}, {probe, 1}};
		const std::vector<VectorN<Real, 3>> points3{{probe, 0, 0}, {probe, 1, 0}};
		const std::vector<IRealFunction*> functions{&function};
		const std::vector<LinearInterpRealFunc> interpolants{LinearInterpRealFunc(parameters, Vector<Real>{probe, probe})};
		ODESystemSolution solution(0, 1, 3, 2);
		solution.fillValues(0, 0, Vector<Real>{probe, 0, 0});
		solution.fillValues(1, 1, Vector<Real>{probe, 1, 0});
		struct Writer {
			std::string name;
			int values;
			std::function<IOResult(const std::string&)> write;
		};
		const std::vector<Writer> writers{
			{"real", 2, [&](const auto& p) { return Serializer::SaveRealFunc(function, "locale", 0, 1, 2, p, precision); }},
			{"real-points", 2, [&](const auto& p) { return Serializer::SaveRealFunc(function, "locale", parameters, p, precision); }},
			{"equally-spaced", 2,
			 [&](const auto& p) { return Serializer::SaveRealFuncEquallySpaced(function, "locale", 0, 1, 2, p, precision); }},
			{"multi", 2, [&](const auto& p) { return Serializer::SaveRealMultiFunc(functions, "locale", {"f"}, 0, 1, 2, p, precision); }},
			{"interpolants", 2, [&](const auto& p) { return Serializer::SaveRealMultiFunc(interpolants, "locale", {"f"}, 0, 1, 2, p); }},
			{"curve2", 2, [&](const auto& p) { return Serializer::SaveParamCurveCartesian2DResult(curve2, "locale", 0, 1, 2, p); }},
			{"curve3", 2, [&](const auto& p) { return Serializer::SaveParamCurveCartesian3DResult(curve3, "locale", 0, 1, 2, p); }},
			{"curve-points", 2, [&](const auto& p) { return Serializer::SaveParamCurve<3>(curve3, "CURVE", "locale", parameters, p); }},
			{"raw-curve", 2, [&](const auto& p) { return Serializer::SaveAsParamCurve<3>(points3, "CURVE", "locale", 0, 1, 2, p); }},
			{"raw-curve-points", 2, [&](const auto& p) { return Serializer::SaveAsParamCurve<3>(points3, "CURVE", "locale", parameters, p); }},
			{"xy-curve", 2, [&](const auto& p) { return Serializer::SaveAsParamCurve2D(Vector<Real>{probe, probe}, parameters, "locale", p); }},
			{"surface", 4,
			 [&](const auto& p) {
				 return Serializer::SaveParametricSurface(static_cast<const IParametricSurface<3>&>(surface), "locale", 0, 1, 2, 0, 1, 2, p);
			 }},
			{"surface-vector", 4,
			 [&](const auto& p) {
				 return Serializer::SaveParametricSurface(
					 static_cast<const IVectorFunctionNM<2, 3>&>(static_cast<const IParametricSurface<3>&>(surface)), "locale", 0, 1, 2, 0, 1, 2, p);
			 }},
			{"surface-rect", 4,
			 [&](const auto& p) {
				 return Serializer::SaveParametricSurface(static_cast<const IParametricSurfaceRect<3>&>(surface), "locale", 2, 2, p);
			 }},
			{"scalar2", 4, [&](const auto& p) { return Serializer::SaveScalarFunc2DCartesian(scalar2, "locale", 0, 1, 2, 0, 1, 2, p); }},
			{"scalar3", 8, [&](const auto& p) { return Serializer::SaveScalarFunc3DCartesian(scalar3, "locale", 0, 1, 2, 0, 1, 2, 0, 1, 2, p); }},
			{"grid", 4, [&](const auto& p) { return Serializer::SaveGridFunction2D(grid, "locale", p); }},
			{"vector2", 4, [&](const auto& p) { return Serializer::SaveVectorFunc2DCartesian(vector2, "locale", 0, 1, 2, 0, 1, 2, p); }},
			{"vector2-filter", 4,
			 [&](const auto& p) { return Serializer::SaveVectorFunc2DCartesian(vector2, "locale", 0, 1, 2, 0, 1, 2, p, 2 * probe); }},
			{"vector3", 8, [&](const auto& p) { return Serializer::SaveVectorFunc3DCartesian(vector3, "locale", 0, 1, 2, 0, 1, 2, 0, 1, 2, p); }},
			{"vector3-filter", 8,
			 [&](const auto& p) { return Serializer::SaveVectorFunc3DCartesian(vector3, "locale", 0, 1, 2, 0, 1, 2, 0, 1, 2, p, 2 * probe); }},
			{"spherical-vector", 8,
			 [&](const auto& p) { return Serializer::SaveVectorFuncSpherical(vector3, "locale", 1, 2, 2, 0, 1, 2, 0, 1, 2, p); }},
			{"spherical-vector-filter", 8,
			 [&](const auto& p) { return Serializer::SaveVectorFuncSpherical(vector3, "locale", 1, 2, 2, 0, 1, 2, 0, 1, 2, p, 2 * probe); }},
			{"ode-component", 2, [&](const auto& p) { return Serializer::SaveODESolutionComponentAsFunc(solution, 0, "locale", p); }},
			{"ode-multi", 2, [&](const auto& p) { return Serializer::SaveODESolutionAsMultiFunc(solution, "locale", {"x", "y", "z"}, p); }},
			{"ode-curve2", 2, [&](const auto& p) { return Serializer::SaveODESolAsParametricCurve2D(solution, p, 0, 1, "locale"); }},
			{"ode-curve3", 2, [&](const auto& p) { return Serializer::SaveODESolAsParametricCurve3D(solution, p, 0, 1, 2, "locale"); }},
			{"particles2", 2,
			 [&](const auto& p) {
				 return Serializer::SaveParticleSimulation2D(p, 1, 10, 10, {{Pnt2Cart(probe, 0), Pnt2Cart(probe, 1)}}, {"red"}, {1}, Real(0.1));
			 }},
			{"particles3", 2,
			 [&](const auto& p) {
				 return Serializer::SaveParticleSimulation3D(p, 1, 10, 10, 10, {{Pnt3Cart(probe, 0, 0), Pnt3Cart(probe, 1, 0)}}, {"red"}, {1},
																										 Real(0.1));
			 }},
			{"field-lines2", 2, [&](const auto& p) { return Serializer::SaveFieldLines2D({points2}, "locale", p); }},
			{"field-lines3", 2, [&](const auto& p) { return Serializer::SaveFieldLines3D({points3}, "locale", p); }},
		};
		for (const auto& writer : writers) {
			DYNAMIC_SECTION(writer.name) {
				const auto path = GetTempFilePath("locale_all_" + writer.name);
				std::string reference;
				{
					GlobalLocaleGuard classic(std::locale::classic());
					const auto result = writer.write(path);
					INFO(result.message);
					REQUIRE(result.success);
					reference = ReadFileContents(path);
					CleanupTempFile(path);
				}
				GlobalLocaleGuard comma;
				const auto result = writer.write(path);
				INFO(result.message);
				REQUIRE(result.success);
				const auto text = ReadFileContents(path);
				CleanupTempFile(path);
				REQUIRE(text == reference);
				REQUIRE(text.find(',') == std::string::npos);
				std::istringstream tokens(text);
				tokens.imbue(std::locale::classic());
				std::string token;
				int readBack = 0;
				while (tokens >> token) {
					std::istringstream number(token);
					number.imbue(std::locale::classic());
					Real value = 0;
					if (number >> value && number.peek() == std::char_traits<char>::eof() && value > 1234 && value < 1235) {
						INFO(token);
						REQUIRE(value == probe);
						++readBack;
					}
				}
				REQUIRE(readBack == writer.values);
			}
		}
	}

	TEST_CASE("Persistence Real corpus round trips exactly with decimal comma and full precision",
						"[persistence][locale][round-trip][wewi.2.5]") {
		GlobalLocaleGuard comma;
		const Vector<Real> corpus{Real(0.1),
															Real(1) / 3,
															std::nextafter(Real(1), Real(2)),
															std::nextafter(Real(1), Real(0)),
															Real(-0.1),
															Real(0),
															-Real(0),
															std::numeric_limits<Real>::min(),
															std::numeric_limits<Real>::max(),
															std::sqrt(std::numeric_limits<Real>::min()),
															-std::sqrt(std::numeric_limits<Real>::max())};
		Matrix<Real> original(1, corpus.size());
		for (int i = 0; i < corpus.size(); ++i)
			original(0, i) = corpus[i];
		for (int format = 0; format < 3; ++format) {
			CAPTURE(format);
			const auto path = GetTempFilePath("exact_corpus_" + std::to_string(format));
			const auto saved = format == 0	 ? Persistence::SaveMatrixToFile(original, path)
												 : format == 1 ? Persistence::SaveMatrixToCSV(original, path)
																			 : Persistence::SaveJson(original, path);
			REQUIRE(saved.success);
			Matrix<Real> restored;
			const auto loaded = format == 0		? Persistence::LoadMatrixFromFile(path, restored)
													: format == 1 ? Persistence::LoadMatrixFromCSV(path, restored)
																				: Persistence::LoadJson(path, restored);
			CleanupTempFile(path);
			INFO(loaded.message);
			REQUIRE(loaded.success);
			REQUIRE(restored.rows() == 1);
			REQUIRE(restored.cols() == corpus.size());
			for (int i = 0; i < corpus.size(); ++i) {
				CAPTURE(i);
				REQUIRE(restored(0, i) == corpus[i]);
				REQUIRE(std::signbit(restored(0, i)) == std::signbit(corpus[i]));
			}
		}
	}

	TEST_CASE("Serializer text output is locale independent and written at round-trip precision", "[serializer][locale][wewi.1.4]") {
		GlobalLocaleGuard commaLocale;

		// Sanity: the guard really changed how a default stream prints
		{ std::ostringstream probe; probe << static_cast<MML::Real>(1234.5); REQUIRE(probe.str() == "1.234,5"); }

		LocaleProbeFunction f;
		const std::string path = GetTempFilePath("locale_realfunc");
		auto result = Serializer::SaveRealFunc(f, "locale", static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 3, path, std::numeric_limits<Real>::max_digits10);
		REQUIRE(result.success);

		const std::string text = ReadFileContents(path);
		CleanupTempFile(path);
		INFO(text);
		REQUIRE(text.find(',') == std::string::npos);            // no comma decimal points, no grouping
		// Round-trip precision: the last data line must parse back to exactly f(1)
		{
			std::istringstream lines(text);
			std::string line, lastData;
			while (std::getline(lines, line))
				if (!line.empty() && (std::isdigit(static_cast<unsigned char>(line[0])) || line[0] == '-'))
					lastData = line;
			std::istringstream parse(lastData);
			parse.imbue(std::locale::classic());
			Real x = 0, y = 0;
			parse >> x >> y;
			REQUIRE_FALSE(parse.fail());
			REQUIRE(x == static_cast<MML::Real>(1.0));
			REQUIRE(y == f(static_cast<MML::Real>(1.0)));
		}

		// Stream version: the caller's stream is put into classic format as well
		std::ostringstream out;
		auto streamResult = Serializer::SaveRealFunc(out, f, "locale", static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 3);
		REQUIRE(streamResult.success);
		REQUIRE(out.str().find(',') == std::string::npos);
	}

	TEST_CASE("Persistence text/CSV/JSON round trips are locale independent", "[persistence][locale][wewi.1.4]") {
		GlobalLocaleGuard commaLocale;

		Matrix<Real> m(2, 2);
		m(0, 0) = static_cast<MML::Real>(1234.5); m(0, 1) = static_cast<MML::Real>(0.1);
		m(1, 0) = -static_cast<MML::Real>(2.5e-7); m(1, 1) = static_cast<MML::Real>(1.0) / static_cast<MML::Real>(3.0);

		SECTION("text format") {
			const std::string path = (std::filesystem::temp_directory_path() / "mml_locale_matrix.txt").string();
			REQUIRE(Persistence::SaveMatrixToFile(m, path).success);
			REQUIRE(ReadFileContents(path).find(',') == std::string::npos);
			Matrix<Real> back;
			REQUIRE(Persistence::LoadMatrixFromFile(path, back).success);
			std::filesystem::remove(path);
			for (int i = 0; i < 2; i++)
				for (int j = 0; j < 2; j++)
					REQUIRE(back(i, j) == m(i, j));   // bit-exact: max_digits10 round trip
		}
		SECTION("CSV format") {
			const std::string path = (std::filesystem::temp_directory_path() / "mml_locale_matrix.csv").string();
			REQUIRE(Persistence::SaveMatrixToCSV(m, path).success);
			Matrix<Real> back;
			REQUIRE(Persistence::LoadMatrixFromCSV(path, back).success);
			std::filesystem::remove(path);
			REQUIRE(back.rows() == 2);
			REQUIRE(back.cols() == 2);
			for (int i = 0; i < 2; i++)
				for (int j = 0; j < 2; j++)
					REQUIRE(back(i, j) == m(i, j));
		}
		SECTION("CSV reader reports unparsable cells instead of storing garbage") {
			const std::string path = (std::filesystem::temp_directory_path() / "mml_bad_matrix.csv").string();
			{ std::ofstream bad(path); bad << "1,2\n3,abc\n"; }
			Matrix<Real> back;
			auto r = Persistence::LoadMatrixFromCSV(path, back);
			std::filesystem::remove(path);
			REQUIRE_FALSE(r.success);
			REQUIRE(r.error == IOError::INVALID_FORMAT);
		}
		SECTION("JSON parse and write") {
			auto parsed = Persistence::ParseJson("{\"v\": 1234.5, \"e\": -2.5e-7}");
			REQUIRE(parsed.result.success);
			REQUIRE(parsed.value.object_value.at("v").number_value == 1234.5);
			REQUIRE(parsed.value.object_value.at("e").number_value == -2.5e-7);

			const std::string json = Persistence::ToJsonString(parsed.value);
			REQUIRE(json.find(',') != std::string::npos);        // the separator between members...
			REQUIRE(json.find("1234.5") != std::string::npos);   // ...but the number uses '.'
			REQUIRE(json.find("1.234,5") == std::string::npos);
		}
	}

	TEST_CASE("Vector persistence is locale independent and round-trips exactly", "[Vector][locale][wewi.1.4]") {
		GlobalLocaleGuard commaLocale;
		Vector<Real> v({static_cast<MML::Real>(1234.5), static_cast<MML::Real>(0.1), static_cast<MML::Real>(1.0) / static_cast<MML::Real>(3.0), -static_cast<MML::Real>(2.5e-7)});
		const std::string path = (std::filesystem::temp_directory_path() / "mml_locale_vector.txt").string();
		REQUIRE(Persistence::Save(v, path).success);
		REQUIRE(ReadFileContents(path).find(',') == std::string::npos);
		Vector<Real> back;
		REQUIRE(Persistence::Load(path, back).success);
		std::filesystem::remove(path);
		REQUIRE(back.size() == 4);
		for (int i = 0; i < 4; i++)
			REQUIRE(back[i] == v[i]);
		REQUIRE_FALSE(Persistence::Load((std::filesystem::temp_directory_path() / "mml_does_not_exist.txt").string(), back).success);
	}

	TEST_CASE("Serializer reports write failures instead of Success", "[serializer][errors][wewi.1.4]") {
		// A stream in a failed state: before the fix the try/catch never fired (streams do not
		// throw) and the header result was discarded, so this returned {true, OK, "Success"}.
		std::ostringstream broken;
		broken.setstate(std::ios::badbit);
		LocaleProbeFunction f;

		auto result = Serializer::SaveRealFunc(broken, f, "broken", static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 3);
		REQUIRE_FALSE(result.success);
		REQUIRE(result.error == IOError::WRITE_FAILED);

		std::ostringstream brokenMulti;
		brokenMulti.setstate(std::ios::failbit);
		std::vector<IRealFunction*> funcs{&f};
		auto multi = Serializer::SaveRealMultiFunc(brokenMulti, funcs, "broken", {"f"}, static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), 3);
		REQUIRE_FALSE(multi.success);
		REQUIRE(multi.error == IOError::WRITE_FAILED);

		// Header writers themselves report the stream state
		std::ostringstream brokenHeader;
		brokenHeader.setstate(std::ios::badbit);
		REQUIRE_FALSE(Serializer::WriteVectorFieldHeader(brokenHeader, "T", "title").success);
	}

} // namespace MML::Tests::Tools::SerializerTests
