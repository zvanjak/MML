///////////////////////////////////////////////////////////////////////////////////////////
// Vector I/O Tests - Round-trip tests for Vector file serialization
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mml/base/Vector/Vector.h>

#include <filesystem>
#include <fstream>
#include <cstring>

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
    bool VectorsEqual(const Vector<Real>& a, const Vector<Real>& b, Real relTol = 1e-14) {
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
    
    // Helper to check exact binary equality (no floating-point tolerance)
    bool VectorsExactlyEqual(const Vector<Real>& a, const Vector<Real>& b) {
        if (a.size() != b.size())
            return false;
        for (int i = 0; i < a.size(); ++i) {
            // Bit-exact comparison using memcmp
            Real va = a[i], vb = b[i];
            if (std::memcmp(&va, &vb, sizeof(Real)) != 0)
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
    REQUIRE_THAT(loaded[0], WithinRel(original[0], 1e-5));
}

TEST_CASE("Vector I/O - Text format special values", "[Vector][IO][Text]")
{
    std::string path = TempFilePath("vector_special.txt");
    TempFileGuard guard(path);
    
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
//                           BINARY FORMAT TESTS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Vector I/O - Binary format round-trip exact", "[Vector][IO][Binary]")
{
    std::string path = TempFilePath("vector.mmlv");
    TempFileGuard guard(path);
    
    Vector<Real> original(10);
    FillVector(original);
    
    REQUIRE(Vector<Real>::SaveToBinary(original, path));
    
    Vector<Real> loaded;
    REQUIRE(Vector<Real>::LoadFromBinary(path, loaded));
    
    REQUIRE(loaded.size() == original.size());
    
    // Binary format should preserve EXACT values (bit-for-bit)
    REQUIRE(VectorsExactlyEqual(original, loaded));
}

TEST_CASE("Vector I/O - Binary format large vector", "[Vector][IO][Binary]")
{
    std::string path = TempFilePath("vector_large.mmlv");
    TempFileGuard guard(path);
    
    Vector<Real> original(10000);
    FillVector(original);
    
    REQUIRE(Vector<Real>::SaveToBinary(original, path));
    
    Vector<Real> loaded;
    REQUIRE(Vector<Real>::LoadFromBinary(path, loaded));
    
    REQUIRE(loaded.size() == 10000);
    REQUIRE(VectorsExactlyEqual(original, loaded));
}

TEST_CASE("Vector I/O - Binary format special values", "[Vector][IO][Binary]")
{
    std::string path = TempFilePath("vector_special.mmlv");
    TempFileGuard guard(path);
    
    Vector<Real> original(9);
    original[0] = 0.0;
    original[1] = -0.0;
    original[2] = std::numeric_limits<Real>::epsilon();
    original[3] = std::numeric_limits<Real>::min();
    original[4] = std::numeric_limits<Real>::max();
    original[5] = std::numeric_limits<Real>::denorm_min();
    original[6] = 1e-300;
    original[7] = 1e+300;
    original[8] = -1.7976931348623157e+308;
    
    REQUIRE(Vector<Real>::SaveToBinary(original, path));
    
    Vector<Real> loaded;
    REQUIRE(Vector<Real>::LoadFromBinary(path, loaded));
    
    // Binary should be bit-exact for all these special values
    REQUIRE(VectorsExactlyEqual(original, loaded));
}

TEST_CASE("Vector I/O - Binary format single element", "[Vector][IO][Binary]")
{
    std::string path = TempFilePath("vector_1.mmlv");
    TempFileGuard guard(path);
    
    Vector<Real> original(1);
    original[0] = 3.14159265358979323846;
    
    REQUIRE(Vector<Real>::SaveToBinary(original, path));
    
    Vector<Real> loaded;
    REQUIRE(Vector<Real>::LoadFromBinary(path, loaded));
    
    REQUIRE(VectorsExactlyEqual(original, loaded));
}

TEST_CASE("Vector I/O - Binary format empty vector", "[Vector][IO][Binary]")
{
    std::string path = TempFilePath("vector_empty.mmlv");
    TempFileGuard guard(path);
    
    Vector<Real> original(0);
    
    REQUIRE(Vector<Real>::SaveToBinary(original, path));
    
    Vector<Real> loaded;
    REQUIRE(Vector<Real>::LoadFromBinary(path, loaded));
    
    REQUIRE(loaded.size() == 0);
}

TEST_CASE("Vector I/O - Binary format header verification", "[Vector][IO][Binary]")
{
    std::string path = TempFilePath("vector_header.mmlv");
    TempFileGuard guard(path);
    
    Vector<Real> original(42);
    FillVector(original);
    
    REQUIRE(Vector<Real>::SaveToBinary(original, path));
    
    // Manually verify header
    std::ifstream file(path, std::ios::binary);
    REQUIRE(file.is_open());
    
    uint32_t magic, version, size, elemSize;
    file.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    file.read(reinterpret_cast<char*>(&version), sizeof(version));
    file.read(reinterpret_cast<char*>(&size), sizeof(size));
    file.read(reinterpret_cast<char*>(&elemSize), sizeof(elemSize));
    
    REQUIRE(magic == MML::BinaryFormat::MAGIC_VECTOR);
    REQUIRE(version == MML::BinaryFormat::VERSION_VECTOR);
    REQUIRE(size == 42);
    REQUIRE(elemSize == sizeof(Real));
}

///////////////////////////////////////////////////////////////////////////////////////////
//                           ERROR HANDLING TESTS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Vector I/O - Load from non-existent file", "[Vector][IO][Error]")
{
    Vector<Real> loaded;
    REQUIRE_FALSE(Vector<Real>::LoadFromFile("nonexistent_vector_file_12345.txt", loaded));
    REQUIRE_FALSE(Vector<Real>::LoadFromBinary("nonexistent_vector_file_12345.mmlv", loaded));
}

TEST_CASE("Vector I/O - Binary wrong magic number", "[Vector][IO][Error]")
{
    std::string path = TempFilePath("bad_magic.mmlv");
    TempFileGuard guard(path);
    
    // Write file with wrong magic number
    std::ofstream file(path, std::ios::binary);
    uint32_t badMagic = 0x12345678;
    file.write(reinterpret_cast<const char*>(&badMagic), sizeof(badMagic));
    file.close();
    
    Vector<Real> loaded;
    REQUIRE_FALSE(Vector<Real>::LoadFromBinary(path, loaded));
}

TEST_CASE("Vector I/O - Binary truncated file", "[Vector][IO][Error]")
{
    std::string path = TempFilePath("truncated.mmlv");
    TempFileGuard guard(path);
    
    // Write only part of header
    std::ofstream file(path, std::ios::binary);
    uint32_t magic = MML::BinaryFormat::MAGIC_VECTOR;
    file.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
    // Don't write the rest
    file.close();
    
    Vector<Real> loaded;
    // Should fail or handle gracefully (not crash)
    bool result = Vector<Real>::LoadFromBinary(path, loaded);
    // Either fails or loads empty - should not crash
    (void)result;
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

TEST_CASE("Vector I/O - Static SaveToBinary API", "[Vector][IO][Static]")
{
    std::string path = TempFilePath("static_binary.mmlv");
    TempFileGuard guard(path);
    
    Vector<Real> original(10);
    FillVector(original);
    
    // Test static method
    REQUIRE(Vector<Real>::SaveToBinary(original, path));
    
    Vector<Real> loaded;
    REQUIRE(Vector<Real>::LoadFromBinary(path, loaded));
    REQUIRE(VectorsExactlyEqual(original, loaded));
}

///////////////////////////////////////////////////////////////////////////////////////////
//                           PRECISION PRESERVATION TESTS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Vector I/O - Binary preserves full precision", "[Vector][IO][Precision]")
{
    std::string path = TempFilePath("precision.mmlv");
    TempFileGuard guard(path);
    
    // Values that are difficult to preserve in text format
    Vector<Real> original(5);
    original[0] = 0.1;  // Cannot be represented exactly in binary floating-point
    original[1] = 0.2;
    original[2] = 0.1 + 0.2;  // Not exactly 0.3!
    original[3] = 1.0 / 7.0;  // Repeating decimal
    original[4] = std::sqrt(2.0);  // Irrational
    
    REQUIRE(Vector<Real>::SaveToBinary(original, path));
    
    Vector<Real> loaded;
    REQUIRE(Vector<Real>::LoadFromBinary(path, loaded));
    
    // Binary format must preserve exact bit patterns
    REQUIRE(VectorsExactlyEqual(original, loaded));
    
    // Verify the tricky 0.1 + 0.2 != 0.3 case is preserved exactly
    REQUIRE(loaded[2] == original[2]);
    REQUIRE(loaded[2] != 0.3);  // This should still be true!
}

TEST_CASE("Vector I/O - Text vs Binary precision comparison", "[Vector][IO][Precision]")
{
    std::string textPath = TempFilePath("vec_precision_compare.txt");
    std::string binPath = TempFilePath("vec_precision_compare.mmlv");
    TempFileGuard guardText(textPath);
    TempFileGuard guardBin(binPath);
    
    Vector<Real> original(3);
    original[0] = 1.0 / 3.0;  // 0.333333...
    original[1] = std::sqrt(2.0);
    original[2] = std::exp(1.0);
    
    REQUIRE(Vector<Real>::SaveToFile(original, textPath));
    REQUIRE(Vector<Real>::SaveToBinary(original, binPath));
    
    Vector<Real> fromText, fromBin;
    REQUIRE(Vector<Real>::LoadFromFile(textPath, fromText));
    REQUIRE(Vector<Real>::LoadFromBinary(binPath, fromBin));
    
    // Binary should be exact
    REQUIRE(VectorsExactlyEqual(original, fromBin));
    
    // Text should be close but may not be exact (limited to ~6 digits)
    REQUIRE(VectorsEqual(original, fromText, 1e-5));
}
