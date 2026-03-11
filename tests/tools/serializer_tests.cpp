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
#include <fstream>
#include <sstream>
#include <cmath>
#include <filesystem>
#include <vector>
#include <string>
#include <regex>

#include "MMLBase.h"
#include "tools/Serializer.h"
#include "interfaces/IFunction.h"
#include "base/InterpolatedFunction.h"

using namespace MML;

// Alias for cleaner error type access (now at MML namespace level, not nested in class)
using MML::SerializeError;

#define TEST_PRECISION_INFO() \
	INFO("Test precision: " << Constants::Eps << " (" << (sizeof(Real) == 8 ? "double" : "float") << ")")

namespace MML::Tests::Tools::SerializerTests {

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
	                                         Real tolerance = 1e-10) {
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

	TEST_CASE("Serializer - SaveRealFunc equally spaced", "[serializer][realfunc]") {
		TEST_PRECISION_INFO();
		
		std::string testFile = GetTempFilePath("realfunc_equally_spaced");
		
		SECTION("Successful serialization with IRealFunction") {
			SinFunction sinFunc;
			auto result = Serializer::SaveRealFuncEquallySpaced(sinFunc, "sin(x)", 0.0, Constants::PI, 10, testFile);
			
			REQUIRE(result.success == true);
			REQUIRE(result.error == SerializeError::OK);
			
			// Verify file exists and has content
			std::string content = ReadFileContents(testFile);
			REQUIRE(content.find("REAL_FUNCTION") != std::string::npos);
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
			REQUIRE(content.find("REAL_FUNCTION") != std::string::npos);
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
			REQUIRE(content.find("MULTI_REAL_FUNCTION") != std::string::npos);
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
			REQUIRE(content.find("MULTI_REAL_FUNCTION") != std::string::npos);
			REQUIRE(content.find("Triangle1") != std::string::npos);
			
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
			REQUIRE(content.find("PARAMETRIC_CURVE_CARTESIAN_2D") != std::string::npos);
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
			REQUIRE(content.find("PARAMETRIC_CURVE_CARTESIAN_2D") != std::string::npos);
			
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
			REQUIRE(content.find("PARAMETRIC_CURVE_CARTESIAN_3D") != std::string::npos);
			
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
			REQUIRE(content.find("PARAMETRIC_SURFACE_CARTESIAN") != std::string::npos);
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
			REQUIRE(content.find("SCALAR_FUNCTION_CARTESIAN_2D") != std::string::npos);
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
			REQUIRE(content.find("SCALAR_FUNCTION_CARTESIAN_3D") != std::string::npos);
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
			REQUIRE(content.find("VECTOR_FIELD_2D_CARTESIAN") != std::string::npos);
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
			REQUIRE(content.find("VECTOR_FIELD_2D_CARTESIAN") != std::string::npos);
			
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
			REQUIRE(content.find("VECTOR_FIELD_3D_CARTESIAN") != std::string::npos);
			
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
			REQUIRE(content.find("VECTOR_FIELD_SPHERICAL") != std::string::npos);
			
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
			REQUIRE(content.find("PARTICLE_SIMULATION_DATA_2D") != std::string::npos);
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
			REQUIRE(content.find("PARTICLE_SIMULATION_DATA_3D") != std::string::npos);
			REQUIRE(content.find("Depth:") != std::string::npos);
			
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
			REQUIRE(result.error == SerializeError::INVALID_PARAMETERS);
		}

		SECTION("Invalid range (x1 >= x2)") {
			SinFunction sinFunc;
			std::string testFile = GetTempFilePath("error_range");
			
			auto result = Serializer::SaveRealFuncEquallySpaced(sinFunc, "test", 5.0, 1.0, 10, testFile);
			
			REQUIRE(result.success == false);
			REQUIRE(result.error == SerializeError::INVALID_PARAMETERS);
		}

		SECTION("Invalid numPoints (< 2)") {
			SinFunction sinFunc;
			std::string testFile = GetTempFilePath("error_numpoints");
			
			auto result = Serializer::SaveRealFuncEquallySpaced(sinFunc, "test", 0.0, 1.0, 1, testFile);
			
			REQUIRE(result.success == false);
			REQUIRE(result.error == SerializeError::INVALID_PARAMETERS);
		}

		SECTION("Mismatched funcs and legend sizes") {
			SinFunction sinFunc;
			std::vector<IRealFunction*> funcs = {&sinFunc};
			std::vector<std::string> legend = {"func1", "func2"};  // Too many legends
			std::string testFile = GetTempFilePath("error_mismatch");
			
			auto result = Serializer::SaveRealMultiFunc(funcs, "test", legend, 0.0, 1.0, 10, testFile);
			
			REQUIRE(result.success == false);
			REQUIRE(result.error == SerializeError::INVALID_PARAMETERS);
		}

		SECTION("Empty funcs vector") {
			std::vector<IRealFunction*> funcs;
			std::vector<std::string> legend;
			std::string testFile = GetTempFilePath("error_empty");
			
			auto result = Serializer::SaveRealMultiFunc(funcs, "test", legend, 0.0, 1.0, 10, testFile);
			
			REQUIRE(result.success == false);
			REQUIRE(result.error == SerializeError::INVALID_PARAMETERS);
		}

		SECTION("Invalid particle simulation parameters") {
			std::string testFile = GetTempFilePath("error_particle");
			
			// numBalls <= 0
			auto result = Serializer::SaveParticleSimulation2D(testFile, 0, 100, 100,
			                                                    {}, {}, {}, 0.01, 1);
			
			REQUIRE(result.success == false);
			REQUIRE(result.error == SerializeError::INVALID_PARAMETERS);
		}

		SECTION("File cannot be opened (invalid path)") {
			SinFunction sinFunc;
			// Try to write to a directory that doesn't exist
			std::string invalidPath = "/nonexistent_dir_12345/test.mml";
			
			auto result = Serializer::SaveRealFuncEquallySpaced(sinFunc, "test", 0.0, 1.0, 10, invalidPath);
			
			REQUIRE(result.success == false);
			REQUIRE(result.error == SerializeError::FILE_NOT_OPENED);
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
			
			Serializer::SaveRealFuncEquallySpaced(sinFunc, "test_title", 0.0, 1.0, 5, testFile);
			
			std::string content = ReadFileContents(testFile);
			std::istringstream iss(content);
			std::string line;
			
			// Line 1: Type
			std::getline(iss, line);
			REQUIRE(line == "REAL_FUNCTION_EQUALLY_SPACED");
			
			// Line 2: Title
			std::getline(iss, line);
			REQUIRE(line == "test_title");
			
			// Line 3: x1
			std::getline(iss, line);
			REQUIRE(line.find("x1:") == 0);
			
			// Line 4: x2
			std::getline(iss, line);
			REQUIRE(line.find("x2:") == 0);
			
			// Line 5: NumPoints
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
			
			Serializer::SaveRealMultiFunc(funcs, "multi_test", legend, 0.0, 1.0, 3, testFile);
			
			std::string content = ReadFileContents(testFile);
			std::istringstream iss(content);
			std::string line;
			
			// Line 1: Type
			std::getline(iss, line);
			REQUIRE(line == "MULTI_REAL_FUNCTION");
			
			// Line 2: Title
			std::getline(iss, line);
			REQUIRE(line == "multi_test");
			
			// Line 3: Number of functions
			std::getline(iss, line);
			REQUIRE(line == "2");
			
			// Lines 4-5: Legend entries
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
			
			Serializer::SaveRealMultiFunc(funcs, "test", legend, 0.0, 1.0, 5, testFile);
			
			std::string content = ReadFileContents(testFile);
			std::istringstream iss(content);
			std::string line;
			
			// Skip header (8 lines: type, title, numFuncs, 2 legends, x1, x2, NumPoints)
			for (int i = 0; i < 8; ++i) {
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
			
			auto cmp = CompareFilesWithTolerance(testFile, 
			    GetReferenceFilePath("realfunc_equally_spaced_simple.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
		
		SECTION("Complex - 13 points over full period") {
			std::string testFile = GetTempFilePath("golden_realfunc_es_complex");
			auto result = Serializer::SaveRealFuncEquallySpaced(sinFunc, "sin(x)",
			    0.0, 2*Constants::PI, 13, testFile);
			REQUIRE(result.success);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("realfunc_equally_spaced_complex.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
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
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("realfunc_specified_complex.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
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
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("multi_realfunc_simple.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
			CleanupTempFile(testFile);
		}
		
		SECTION("Complex - 11 points over full period") {
			std::string testFile = GetTempFilePath("golden_multi_complex");
			auto result = Serializer::SaveRealMultiFunc(funcs, "Trig Functions", legend,
			    0.0, 2*Constants::PI, 11, testFile);
			REQUIRE(result.success);
			
			auto cmp = CompareFilesWithTolerance(testFile,
			    GetReferenceFilePath("multi_realfunc_complex.mml"));
			INFO("Line " << cmp.line_number << ": " << cmp.error_message);
			REQUIRE(cmp.success);
			
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

} // namespace MML::Tests::Tools::SerializerTests

