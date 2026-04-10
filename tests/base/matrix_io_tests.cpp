///////////////////////////////////////////////////////////////////////////////////////////
// Matrix I/O Tests - Round-trip tests for Matrix file serialization
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_test_macros.hpp>
#include "../TestPrecision.h"
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <mml/base/Matrix/Matrix.h>
#include <mml/tools/serializer/MatrixIO.h>

#include <filesystem>
#include <fstream>
#include <cstring>

using namespace MML;
using namespace MML::Serializer;
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
    
    // Helper to fill matrix with deterministic values
    void FillMatrix(Matrix<Real>& m, Real seed = 1.0) {
        for (int i = 0; i < m.rows(); ++i) {
            for (int j = 0; j < m.cols(); ++j) {
                m(i, j) = seed * (i * 100 + j) + 0.123456789012345;
            }
        }
    }
    
    // Helper to compare matrices with relative tolerance
    bool MatricesEqual(const Matrix<Real>& a, const Matrix<Real>& b, Real relTol = TOL(1e-14, 1e-5)) {
        if (a.rows() != b.rows() || a.cols() != b.cols())
            return false;
        for (int i = 0; i < a.rows(); ++i) {
            for (int j = 0; j < a.cols(); ++j) {
                Real va = a(i, j);
                Real vb = b(i, j);
                Real diff = std::abs(va - vb);
                Real maxAbs = std::max(std::abs(va), std::abs(vb));
                // Use relative tolerance, but handle near-zero values with absolute tolerance
                Real tolerance = std::max(relTol * maxAbs, relTol);
                if (diff > tolerance)
                    return false;
            }
        }
        return true;
    }
    
    // Helper to check exact binary equality (no floating-point tolerance)
    bool MatricesExactlyEqual(const Matrix<Real>& a, const Matrix<Real>& b) {
        if (a.rows() != b.rows() || a.cols() != b.cols())
            return false;
        for (int i = 0; i < a.rows(); ++i) {
            for (int j = 0; j < a.cols(); ++j) {
                Real va = a(i, j), vb = b(i, j);
                // Both NaN → considered equal for round-trip purposes
                if (std::isnan(va) && std::isnan(vb)) continue;
                // Value comparison (avoids comparing padding bytes of long double)
                if (va != vb) return false;
                // Distinguish +0 from -0
                if (va == Real(0) && std::signbit(va) != std::signbit(vb)) return false;
            }
        }
        return true;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
//                           TEXT FORMAT TESTS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Matrix I/O - Text format round-trip small", "[Matrix][IO][Text]")
{
    std::string path = TempFilePath("matrix_small.txt");
    TempFileGuard guard(path);
    
    Matrix<Real> original(3, 4);
    FillMatrix(original);
    
    REQUIRE(SaveMatrixToFile(original, path));
    
    Matrix<Real> loaded;
    REQUIRE(LoadMatrixFromFile(path, loaded));
    
    REQUIRE(loaded.rows() == original.rows());
    REQUIRE(loaded.cols() == original.cols());
    // Text format uses default stream precision (~6 digits)
    REQUIRE(MatricesEqual(original, loaded, 1e-5));
}

TEST_CASE("Matrix I/O - Text format round-trip square", "[Matrix][IO][Text]")
{
    std::string path = TempFilePath("matrix_square.txt");
    TempFileGuard guard(path);
    
    Matrix<Real> original(5, 5);
    FillMatrix(original, 2.5);
    
    REQUIRE(SaveMatrixToFile(original, path));
    
    Matrix<Real> loaded;
    REQUIRE(LoadMatrixFromFile(path, loaded));
    
    // Text format uses default stream precision (~6 digits)
    REQUIRE(MatricesEqual(original, loaded, 1e-5));
}

TEST_CASE("Matrix I/O - Text format round-trip large", "[Matrix][IO][Text]")
{
    std::string path = TempFilePath("matrix_large.txt");
    TempFileGuard guard(path);
    
    Matrix<Real> original(50, 50);
    FillMatrix(original);
    
    REQUIRE(SaveMatrixToFile(original, path));
    
    Matrix<Real> loaded;
    REQUIRE(LoadMatrixFromFile(path, loaded));
    
    REQUIRE(loaded.rows() == 50);
    REQUIRE(loaded.cols() == 50);
    // Text format uses default stream precision (~6 digits)
    REQUIRE(MatricesEqual(original, loaded, 1e-5));
}

TEST_CASE("Matrix I/O - Text format single element", "[Matrix][IO][Text]")
{
    std::string path = TempFilePath("matrix_single.txt");
    TempFileGuard guard(path);
    
    Matrix<Real> original(1, 1);
    original(0, 0) = 3.14159265358979323846;
    
    REQUIRE(SaveMatrixToFile(original, path));
    
    Matrix<Real> loaded;
    REQUIRE(LoadMatrixFromFile(path, loaded));
    
    REQUIRE(loaded.rows() == 1);
    REQUIRE(loaded.cols() == 1);
    // Text format uses default stream precision (~6 digits)
    REQUIRE_THAT(loaded(0, 0), WithinRel(original(0, 0), REAL(1e-5)));
}

TEST_CASE("Matrix I/O - Text format row vector", "[Matrix][IO][Text]")
{
    std::string path = TempFilePath("matrix_row.txt");
    TempFileGuard guard(path);
    
    Matrix<Real> original(1, 10);
    for (int j = 0; j < 10; ++j)
        original(0, j) = j * 1.1;
    
    REQUIRE(SaveMatrixToFile(original, path));
    
    Matrix<Real> loaded;
    REQUIRE(LoadMatrixFromFile(path, loaded));
    
    REQUIRE(loaded.rows() == 1);
    REQUIRE(loaded.cols() == 10);
    // Text format uses default stream precision (~6 digits)
    REQUIRE(MatricesEqual(original, loaded, 1e-5));
}

TEST_CASE("Matrix I/O - Text format column vector", "[Matrix][IO][Text]")
{
    std::string path = TempFilePath("matrix_col.txt");
    TempFileGuard guard(path);
    
    Matrix<Real> original(10, 1);
    for (int i = 0; i < 10; ++i)
        original(i, 0) = i * 2.2;
    
    REQUIRE(SaveMatrixToFile(original, path));
    
    Matrix<Real> loaded;
    REQUIRE(LoadMatrixFromFile(path, loaded));
    
    REQUIRE(loaded.rows() == 10);
    REQUIRE(loaded.cols() == 1);
    // Text format uses default stream precision (~6 digits)
    REQUIRE(MatricesEqual(original, loaded, 1e-5));
}

TEST_CASE("Matrix I/O - Text format special values", "[Matrix][IO][Text]")
{
    std::string path = TempFilePath("matrix_special.txt");
    TempFileGuard guard(path);
    
    // Values like 1e-300 and 1e+300 overflow/underflow for float
    if constexpr (!std::is_same_v<Real, float>) {
    Matrix<Real> original(2, 3);
    original(0, 0) = 0.0;
    original(0, 1) = -0.0;  // Negative zero
    original(0, 2) = 1e-300;  // Very small
    original(1, 0) = 1e+300;  // Very large
    original(1, 1) = -1.23456789012345e-100;  // Negative scientific
    original(1, 2) = 1.0 / 3.0;  // Repeating decimal
    
    REQUIRE(SaveMatrixToFile(original, path));
    
    Matrix<Real> loaded;
    REQUIRE(LoadMatrixFromFile(path, loaded));
    
    // Text format loses precision - use relative tolerance for large/small values
    REQUIRE(MatricesEqual(original, loaded, 1e-5));
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
//                           CSV FORMAT TESTS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Matrix I/O - CSV format round-trip", "[Matrix][IO][CSV]")
{
    std::string path = TempFilePath("matrix.csv");
    TempFileGuard guard(path);
    
    Matrix<Real> original(4, 5);
    FillMatrix(original);
    
    REQUIRE(SaveMatrixToCSV(original, path));
    
    Matrix<Real> loaded;
    REQUIRE(LoadMatrixFromCSV(path, loaded));
    
    REQUIRE(loaded.rows() == original.rows());
    REQUIRE(loaded.cols() == original.cols());
    // CSV format uses default stream precision (~6 digits)
    REQUIRE(MatricesEqual(original, loaded, 1e-5));
}

TEST_CASE("Matrix I/O - CSV format high precision values", "[Matrix][IO][CSV]")
{
    std::string path = TempFilePath("matrix_prec.csv");
    TempFileGuard guard(path);
    
    Matrix<Real> original(2, 2);
    original(0, 0) = 1.123456789012345;
    original(0, 1) = 2.234567890123456;
    original(1, 0) = 3.345678901234567;
    original(1, 1) = 4.456789012345678;
    
    REQUIRE(SaveMatrixToCSV(original, path));
    
    Matrix<Real> loaded;
    REQUIRE(LoadMatrixFromCSV(path, loaded));
    
    // CSV uses default precision (~6 digits), check with appropriate tolerance
    REQUIRE(MatricesEqual(original, loaded, 1e-5));
}

TEST_CASE("Matrix I/O - CSV format large matrix", "[Matrix][IO][CSV]")
{
    std::string path = TempFilePath("matrix_large.csv");
    TempFileGuard guard(path);
    
    Matrix<Real> original(100, 100);
    FillMatrix(original);
    
    REQUIRE(SaveMatrixToCSV(original, path));
    
    Matrix<Real> loaded;
    REQUIRE(LoadMatrixFromCSV(path, loaded));
    
    REQUIRE(loaded.rows() == 100);
    REQUIRE(loaded.cols() == 100);
    // CSV uses default precision (~6 digits)
    REQUIRE(MatricesEqual(original, loaded, 1e-5));
}

///////////////////////////////////////////////////////////////////////////////////////////
//                           BINARY FORMAT TESTS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Matrix I/O - Binary format round-trip exact", "[Matrix][IO][Binary]")
{
    std::string path = TempFilePath("matrix.mmlm");
    TempFileGuard guard(path);
    
    Matrix<Real> original(4, 5);
    FillMatrix(original);
    
    REQUIRE(SaveMatrixToBinary(original, path));
    
    Matrix<Real> loaded;
    REQUIRE(LoadMatrixFromBinary(path, loaded));
    
    REQUIRE(loaded.rows() == original.rows());
    REQUIRE(loaded.cols() == original.cols());
    
    // Binary format should preserve EXACT values (bit-for-bit)
    REQUIRE(MatricesExactlyEqual(original, loaded));
}

TEST_CASE("Matrix I/O - Binary format large matrix", "[Matrix][IO][Binary]")
{
    std::string path = TempFilePath("matrix_large.mmlm");
    TempFileGuard guard(path);
    
    Matrix<Real> original(200, 200);
    FillMatrix(original);
    
    REQUIRE(SaveMatrixToBinary(original, path));
    
    Matrix<Real> loaded;
    REQUIRE(LoadMatrixFromBinary(path, loaded));
    
    REQUIRE(loaded.rows() == 200);
    REQUIRE(loaded.cols() == 200);
    REQUIRE(MatricesExactlyEqual(original, loaded));
}

TEST_CASE("Matrix I/O - Binary format special values", "[Matrix][IO][Binary]")
{
    std::string path = TempFilePath("matrix_special.mmlm");
    TempFileGuard guard(path);
    
    Matrix<Real> original(3, 3);
    original(0, 0) = 0.0;
    original(0, 1) = -0.0;
    original(0, 2) = std::numeric_limits<Real>::epsilon();
    original(1, 0) = std::numeric_limits<Real>::min();
    original(1, 1) = std::numeric_limits<Real>::max();
    original(1, 2) = std::numeric_limits<Real>::denorm_min();
    original(2, 0) = 1e-300;
    original(2, 1) = 1e+300;
    original(2, 2) = -1.7976931348623157e+308;
    
    REQUIRE(SaveMatrixToBinary(original, path));
    
    Matrix<Real> loaded;
    REQUIRE(LoadMatrixFromBinary(path, loaded));
    
    // Binary should be bit-exact for all these special values
    REQUIRE(MatricesExactlyEqual(original, loaded));
}

TEST_CASE("Matrix I/O - Binary format single element", "[Matrix][IO][Binary]")
{
    std::string path = TempFilePath("matrix_1x1.mmlm");
    TempFileGuard guard(path);
    
    Matrix<Real> original(1, 1);
    original(0, 0) = 3.14159265358979323846;
    
    REQUIRE(SaveMatrixToBinary(original, path));
    
    Matrix<Real> loaded;
    REQUIRE(LoadMatrixFromBinary(path, loaded));
    
    REQUIRE(MatricesExactlyEqual(original, loaded));
}

TEST_CASE("Matrix I/O - Binary format header verification", "[Matrix][IO][Binary]")
{
    std::string path = TempFilePath("matrix_header.mmlm");
    TempFileGuard guard(path);
    
    Matrix<Real> original(7, 11);
    FillMatrix(original);
    
    REQUIRE(SaveMatrixToBinary(original, path));
    
    // Manually verify header
    std::ifstream file(path, std::ios::binary);
    REQUIRE(file.is_open());
    
    uint32_t magic, version, rows, cols, elemSize, reserved;
    file.read(reinterpret_cast<char*>(&magic), sizeof(magic));
    file.read(reinterpret_cast<char*>(&version), sizeof(version));
    file.read(reinterpret_cast<char*>(&rows), sizeof(rows));
    file.read(reinterpret_cast<char*>(&cols), sizeof(cols));
    file.read(reinterpret_cast<char*>(&elemSize), sizeof(elemSize));
    file.read(reinterpret_cast<char*>(&reserved), sizeof(reserved));
    
    REQUIRE(magic == MML::BinaryFormat::MAGIC_MATRIX);
    REQUIRE(version == MML::BinaryFormat::VERSION_MATRIX);
    REQUIRE(rows == 7);
    REQUIRE(cols == 11);
    REQUIRE(elemSize == sizeof(Real));
    REQUIRE(reserved == 0);
}

///////////////////////////////////////////////////////////////////////////////////////////
//                           ERROR HANDLING TESTS
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Matrix I/O - Load from non-existent file", "[Matrix][IO][Error]")
{
    Matrix<Real> loaded;
    REQUIRE_FALSE(LoadMatrixFromFile("nonexistent_matrix_file_12345.txt", loaded));
    REQUIRE_FALSE(LoadMatrixFromCSV("nonexistent_matrix_file_12345.csv", loaded));
    REQUIRE_FALSE(LoadMatrixFromBinary("nonexistent_matrix_file_12345.mmlm", loaded));
}

TEST_CASE("Matrix I/O - Binary wrong magic number", "[Matrix][IO][Error]")
{
    std::string path = TempFilePath("bad_magic.mmlm");
    TempFileGuard guard(path);
    
    // Write file with wrong magic number
    std::ofstream file(path, std::ios::binary);
    uint32_t badMagic = 0x12345678;
    file.write(reinterpret_cast<const char*>(&badMagic), sizeof(badMagic));
    file.close();
    
    Matrix<Real> loaded;
    REQUIRE_FALSE(LoadMatrixFromBinary(path, loaded));
}

TEST_CASE("Matrix I/O - Binary truncated file", "[Matrix][IO][Error]")
{
    std::string path = TempFilePath("truncated.mmlm");
    TempFileGuard guard(path);
    
    // Write only part of header
    std::ofstream file(path, std::ios::binary);
    uint32_t magic = MML::BinaryFormat::MAGIC_MATRIX;
    file.write(reinterpret_cast<const char*>(&magic), sizeof(magic));
    // Don't write the rest
    file.close();
    
    Matrix<Real> loaded;
    // Should fail or handle gracefully (not crash)
    bool result = LoadMatrixFromBinary(path, loaded);
    // Either fails or loads empty - should not crash
    (void)result;
}
