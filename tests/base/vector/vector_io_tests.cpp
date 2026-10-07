///////////////////////////////////////////////////////////////////////////////////////////
// Vector I/O Tests - Round-trip tests for Vector file serialization
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include "../../TestPrecision.h"
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mml/base/Vector/Vector.h>

#include <filesystem>
#include <fstream>

using namespace MML;
using Catch::Matchers::WithinRel;

namespace {
    // Helper to create a temporary file path
    std::string TempFilePath(const std::string& name) {
        auto temp = std::filesystem::temp_directory_path() / ("mml_test_" + name);
        return temp.string();
    }
    
    // Helper to clean up test files
    struct TempFileGuard {
        std::string path;
        TempFileGuard(const std::string& p) : path(p) {}
        ~TempFileGuard() { 
            std::filesystem::remove(path); 
        }
    };
    
    // Helper to fill vector with deterministic values
    void FillVector(Vector<Real>& v, Real seed = 1.0) {
        for (int i = 0; i < v.size(); ++i) {
            v[i] = seed * i + 0.123456789012345;
        }
    }
    
    // Helper to compare vectors with relative tolerance
    bool VectorsEqual(const Vector<Real>& a, const Vector<Real>& b, Real relTol = TOL(1e-14, 1e-5)) {
        if (a.size() != b.size())
            return false;
        for (int i = 0; i < a.size(); ++i) {
            Real va = a[i];
            Real vb = b[i];
            Real diff = std::abs(va - vb);
            Real maxAbs = std::max(std::abs(va), std::abs(vb));
            // Use relative tolerance, but handle near-zero values with absolute tolerance
            Real tolerance = std::max(relTol * maxAbs, relTol);
            if (diff > tolerance)
                return false;
        }
        return true;
    }
    
}

///////////////////////////////////////////////////////////////////////////////////////////
//                           TEXT FORMAT TESTS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Vector I/O - Text format round-trip small", "[Vector][IO][Text]")
{
    std::string path = TempFilePath("vector_small.txt");
    TempFileGuard guard(path);
    
    Vector<Real> original(5);
    FillVector(original);
    
    REQUIRE(Vector<Real>::SaveToFile(original, path));
    
    Vector<Real> loaded;
    REQUIRE(Vector<Real>::LoadFromFile(path, loaded));
    
    REQUIRE(loaded.size() == original.size());
    // Text format uses default stream precision (~6 digits)
    REQUIRE(VectorsEqual(original, loaded, 1e-5));
}
TEST_CASE("Vector I/O - Text format round-trip large", "[Vector][IO][Text]")
{
    std::string path = TempFilePath("vector_large.txt");
    TempFileGuard guard(path);
    
    Vector<Real> original(1000);
    FillVector(original);
    
    REQUIRE(Vector<Real>::SaveToFile(original, path));
    
    Vector<Real> loaded;
    REQUIRE(Vector<Real>::LoadFromFile(path, loaded));
    
    REQUIRE(loaded.size() == 1000);
    // Text format uses default stream precision (~6 digits)
    REQUIRE(VectorsEqual(original, loaded, 1e-5));
}
TEST_CASE("Vector I/O - Text format single element", "[Vector][IO][Text]")
{
    std::string path = TempFilePath("vector_single.txt");
    TempFileGuard guard(path);
    
    Vector<Real> original(1);
    original[0] = 3.14159265358979323846;
    
    REQUIRE(Vector<Real>::SaveToFile(original, path));
    
    Vector<Real> loaded;
    REQUIRE(Vector<Real>::LoadFromFile(path, loaded));
    
    REQUIRE(loaded.size() == 1);
    // Text format uses default stream precision (~6 digits)
    REQUIRE_THAT(loaded[0], WithinRel(original[0], REAL(1e-5)));
}

TEST_CASE("Vector I/O - Text format special values", "[Vector][IO][Text]")
{
    std::string path = TempFilePath("vector_special.txt");
    TempFileGuard guard(path);
    
    // Values like 1e-300 and 1e+300 overflow/underflow for float
    if constexpr (!std::is_same_v<Real, float>) {
    Vector<Real> original(6);
    original[0] = 0.0;
    original[1] = -0.0;  // Negative zero
    original[2] = 1e-300;  // Very small
    original[3] = 1e+300;  // Very large
    original[4] = -1.23456789012345e-100;  // Negative scientific
    original[5] = 1.0 / 3.0;  // Repeating decimal
    
    REQUIRE(Vector<Real>::SaveToFile(original, path));
    
    Vector<Real> loaded;
    REQUIRE(Vector<Real>::LoadFromFile(path, loaded));
    
    // Text format loses precision - use relative tolerance
    REQUIRE(VectorsEqual(original, loaded, 1e-5));
    }
}

TEST_CASE("Vector I/O - Text format empty vector", "[Vector][IO][Text]")
{
    std::string path = TempFilePath("vector_empty.txt");
    TempFileGuard guard(path);
    
    Vector<Real> original(0);  // Empty vector
    
    REQUIRE(Vector<Real>::SaveToFile(original, path));
    
    Vector<Real> loaded;
    REQUIRE(Vector<Real>::LoadFromFile(path, loaded));
    
    REQUIRE(loaded.size() == 0);
}

///////////////////////////////////////////////////////////////////////////////////////////
//                           ERROR HANDLING TESTS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Vector I/O - Load from non-existent file", "[Vector][IO][Error]")
{
    Vector<Real> loaded;
    REQUIRE_FALSE(Vector<Real>::LoadFromFile("nonexistent_vector_file_12345.txt", loaded));
}

///////////////////////////////////////////////////////////////////////////////////////////
//                           INSTANCE METHOD ALIASES
///////////////////////////////////////////////////////////////////////////////////////////

// Note: Vector I/O methods are static-only. These tests verify the static API works correctly.

TEST_CASE("Vector I/O - Static SaveToFile API", "[Vector][IO][Static]")
{
    std::string path = TempFilePath("static_save.txt");
    TempFileGuard guard(path);
    
    Vector<Real> original(10);
    FillVector(original);
    
    // Test static method
    REQUIRE(Vector<Real>::SaveToFile(original, path));
    
    Vector<Real> loaded;
    REQUIRE(Vector<Real>::LoadFromFile(path, loaded));
    // Text format uses default stream precision (~6 digits)
    REQUIRE(VectorsEqual(original, loaded, 1e-5));
}

