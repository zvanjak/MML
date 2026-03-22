#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Matrix/Matrix3D.h"
#endif

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Base::Matrix3DTests
{
	// ============================================================================
	// Construction Tests
	// ============================================================================

	TEST_CASE("Matrix3D::DefaultConstruction", "[Matrix3D]") {
		Matrix3D<double> m;
		
		REQUIRE(m.dim1() == 0);
		REQUIRE(m.dim2() == 0);
		REQUIRE(m.dim3() == 0);
		REQUIRE(m.size() == 0);
		REQUIRE(m.isEmpty());
	}

	TEST_CASE("Matrix3D::DimensionConstruction", "[Matrix3D]") {
		Matrix3D<double> m(3, 4, 5);
		
		REQUIRE(m.dim1() == 3);
		REQUIRE(m.dim2() == 4);
		REQUIRE(m.dim3() == 5);
		REQUIRE(m.size() == 60);
		REQUIRE_FALSE(m.isEmpty());
	}

	TEST_CASE("Matrix3D::InitValueConstruction", "[Matrix3D]") {
		Matrix3D<double> m(2, 3, 4, 7.5);
		
		REQUIRE(m.dim1() == 2);
		REQUIRE(m.size() == 24);
		
		// Check all values are initialized
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 3; ++j)
				for (int k = 0; k < 4; ++k)
					REQUIRE(m(i, j, k) == 7.5);
	}

	TEST_CASE("Matrix3D::InvalidDimensions", "[Matrix3D]") {
		REQUIRE_THROWS_AS(Matrix3D<double>(0, 5, 5), std::invalid_argument);
		REQUIRE_THROWS_AS(Matrix3D<double>(5, 0, 5), std::invalid_argument);
		REQUIRE_THROWS_AS(Matrix3D<double>(5, 5, 0), std::invalid_argument);
		REQUIRE_THROWS_AS(Matrix3D<double>(-1, 5, 5), std::invalid_argument);
	}

	// ============================================================================
	// Copy and Move Semantics Tests
	// ============================================================================

	TEST_CASE("Matrix3D::CopyConstruction", "[Matrix3D]") {
		Matrix3D<double> original(2, 3, 4, 1.5);
		original(0, 0, 0) = 99.0;
		original(1, 2, 3) = 77.0;
		
		Matrix3D<double> copy(original);
		
		REQUIRE(copy.dim1() == 2);
		REQUIRE(copy.dim2() == 3);
		REQUIRE(copy.dim3() == 4);
		REQUIRE(copy(0, 0, 0) == 99.0);
		REQUIRE(copy(1, 2, 3) == 77.0);
		REQUIRE(copy(0, 1, 1) == 1.5);
		
		// Verify deep copy - modifying copy doesn't affect original
		copy(0, 0, 0) = 123.0;
		REQUIRE(original(0, 0, 0) == 99.0);
	}

	TEST_CASE("Matrix3D::CopyAssignment", "[Matrix3D]") {
		Matrix3D<double> original(2, 2, 2, 3.14);
		Matrix3D<double> target(5, 5, 5, 0.0);
		
		target = original;
		
		REQUIRE(target.dim1() == 2);
		REQUIRE(target.dim2() == 2);
		REQUIRE(target.dim3() == 2);
		REQUIRE(target(0, 0, 0) == 3.14);
		
		// Self-assignment should be safe
		target = target;
		REQUIRE(target.dim1() == 2);
	}

	TEST_CASE("Matrix3D::MoveConstruction", "[Matrix3D]") {
		Matrix3D<double> original(3, 3, 3, 2.5);
		original(1, 1, 1) = 42.0;
		
		Matrix3D<double> moved(std::move(original));
		
		REQUIRE(moved.dim1() == 3);
		REQUIRE(moved(1, 1, 1) == 42.0);
		REQUIRE(original.isEmpty());  // Original should be empty after move
	}

	TEST_CASE("Matrix3D::MoveAssignment", "[Matrix3D]") {
		Matrix3D<double> original(2, 2, 2, 1.0);
		Matrix3D<double> target(5, 5, 5, 0.0);
		
		target = std::move(original);
		
		REQUIRE(target.dim1() == 2);
		REQUIRE(target(0, 0, 0) == 1.0);
		REQUIRE(original.isEmpty());
	}

	// ============================================================================
	// Element Access Tests
	// ============================================================================

	TEST_CASE("Matrix3D::BracketAccess", "[Matrix3D]") {
		Matrix3D<int> m(3, 4, 5, 0);
		
		// Write using brackets
		m[0][0][0] = 1;
		m[1][2][3] = 42;
		m[2][3][4] = 100;
		
		// Read using brackets
		REQUIRE(m[0][0][0] == 1);
		REQUIRE(m[1][2][3] == 42);
		REQUIRE(m[2][3][4] == 100);
	}

	TEST_CASE("Matrix3D::ParenthesisAccess", "[Matrix3D]") {
		Matrix3D<int> m(3, 4, 5, 0);
		
		// Write using parenthesis
		m(0, 0, 0) = 1;
		m(1, 2, 3) = 42;
		m(2, 3, 4) = 100;
		
		// Read using parenthesis
		REQUIRE(m(0, 0, 0) == 1);
		REQUIRE(m(1, 2, 3) == 42);
		REQUIRE(m(2, 3, 4) == 100);
	}

	TEST_CASE("Matrix3D::AtAccess", "[Matrix3D]") {
		Matrix3D<int> m(2, 3, 4, 5);
		
		REQUIRE(m.At(0, 0, 0) == 5);
		m.At(1, 2, 3) = 99;
		REQUIRE(m.At(1, 2, 3) == 99);
		
		// At should throw on out-of-bounds
		REQUIRE_THROWS_AS(m.At(-1, 0, 0), std::out_of_range);
		REQUIRE_THROWS_AS(m.At(0, 3, 0), std::out_of_range);
		REQUIRE_THROWS_AS(m.At(0, 0, 4), std::out_of_range);
		REQUIRE_THROWS_AS(m.At(2, 0, 0), std::out_of_range);
	}

	TEST_CASE("Matrix3D::DataAccess", "[Matrix3D]") {
		Matrix3D<double> m(2, 2, 2, 0.0);
		m(0, 0, 0) = 1.0;
		m(0, 0, 1) = 2.0;
		m(1, 1, 1) = 8.0;
		
		double* data = m.data();
		REQUIRE(data != nullptr);
		
		// Data should be contiguous
		REQUIRE(data[0] == 1.0);
		REQUIRE(data[1] == 2.0);
		REQUIRE(data[7] == 8.0);  // Element at (1,1,1) in row-major order
	}

	// ============================================================================
	// Utility Methods Tests
	// ============================================================================

	TEST_CASE("Matrix3D::Fill", "[Matrix3D]") {
		Matrix3D<double> m(3, 3, 3);
		m.Fill(42.0);
		
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				for (int k = 0; k < 3; ++k)
					REQUIRE(m(i, j, k) == 42.0);
	}

	TEST_CASE("Matrix3D::Apply", "[Matrix3D]") {
		Matrix3D<double> m(2, 2, 2, 4.0);
		
		m.Apply([](double x) { return x * 2.0; });
		
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 2; ++j)
				for (int k = 0; k < 2; ++k)
					REQUIRE(m(i, j, k) == 8.0);
	}

	TEST_CASE("Matrix3D::Resize", "[Matrix3D]") {
		Matrix3D<double> m(2, 2, 2, 1.0);
		
		m.Resize(3, 4, 5);
		
		REQUIRE(m.dim1() == 3);
		REQUIRE(m.dim2() == 4);
		REQUIRE(m.dim3() == 5);
		REQUIRE(m.size() == 60);
		
		// Resize with init value
		m.Resize(2, 2, 2, 99.0);
		REQUIRE(m(0, 0, 0) == 99.0);
		REQUIRE(m(1, 1, 1) == 99.0);
	}

	TEST_CASE("Matrix3D::ResizeSameDimensions", "[Matrix3D]") {
		Matrix3D<double> m(3, 3, 3, 5.0);
		double* originalData = m.data();
		
		// Resizing to same dimensions should be a no-op
		m.Resize(3, 3, 3);
		
		REQUIRE(m.data() == originalData);
		REQUIRE(m(0, 0, 0) == 5.0);
	}

	// ============================================================================
	// Comparison Tests
	// ============================================================================

	TEST_CASE("Matrix3D::Equality", "[Matrix3D]") {
		Matrix3D<int> m1(2, 2, 2, 5);
		Matrix3D<int> m2(2, 2, 2, 5);
		Matrix3D<int> m3(2, 2, 2, 6);
		Matrix3D<int> m4(3, 2, 2, 5);
		
		REQUIRE(m1 == m2);
		REQUIRE_FALSE(m1 == m3);  // Different values
		REQUIRE_FALSE(m1 == m4);  // Different dimensions
		
		REQUIRE(m1 != m3);
		REQUIRE_FALSE(m1 != m2);
	}

	// ============================================================================
	// Scalar Operations Tests
	// ============================================================================

	TEST_CASE("Matrix3D::ScalarMultiplication", "[Matrix3D]") {
		Matrix3D<double> m(2, 2, 2, 3.0);
		
		m *= 2.0;
		REQUIRE(m(0, 0, 0) == 6.0);
		
		auto m2 = m * 0.5;
		REQUIRE(m2(0, 0, 0) == 3.0);
		REQUIRE(m(0, 0, 0) == 6.0);  // Original unchanged
	}

	TEST_CASE("Matrix3D::ScalarDivision", "[Matrix3D]") {
		Matrix3D<double> m(2, 2, 2, 10.0);
		
		m /= 2.0;
		REQUIRE(m(0, 0, 0) == 5.0);
		
		auto m2 = m / 5.0;
		REQUIRE(m2(0, 0, 0) == 1.0);
	}

	TEST_CASE("Matrix3D::ScalarAddition", "[Matrix3D]") {
		Matrix3D<double> m(2, 2, 2, 5.0);
		
		m += 3.0;
		REQUIRE(m(0, 0, 0) == 8.0);
		
		auto m2 = m + 2.0;
		REQUIRE(m2(0, 0, 0) == 10.0);
	}

	// ============================================================================
	// Element-wise Operations Tests
	// ============================================================================

	TEST_CASE("Matrix3D::ElementWiseAddition", "[Matrix3D]") {
		Matrix3D<double> m1(2, 2, 2, 3.0);
		Matrix3D<double> m2(2, 2, 2, 4.0);
		
		auto m3 = m1 + m2;
		REQUIRE(m3(0, 0, 0) == 7.0);
		REQUIRE(m3(1, 1, 1) == 7.0);
		
		m1 += m2;
		REQUIRE(m1(0, 0, 0) == 7.0);
	}

	TEST_CASE("Matrix3D::ElementWiseSubtraction", "[Matrix3D]") {
		Matrix3D<double> m1(2, 2, 2, 10.0);
		Matrix3D<double> m2(2, 2, 2, 3.0);
		
		auto m3 = m1 - m2;
		REQUIRE(m3(0, 0, 0) == 7.0);
		
		m1 -= m2;
		REQUIRE(m1(0, 0, 0) == 7.0);
	}

	TEST_CASE("Matrix3D::DimensionMismatch", "[Matrix3D]") {
		Matrix3D<double> m1(2, 2, 2, 1.0);
		Matrix3D<double> m2(3, 3, 3, 1.0);
		
		REQUIRE_THROWS_AS(m1 + m2, std::invalid_argument);
		REQUIRE_THROWS_AS(m1 - m2, std::invalid_argument);
		REQUIRE_THROWS_AS(m1 += m2, std::invalid_argument);
	}

	// ============================================================================
	// String Output Tests
	// ============================================================================

	TEST_CASE("Matrix3D::ToString", "[Matrix3D]") {
		Matrix3D<double> m(3, 4, 5);
		
		std::string str = m.ToString();
		REQUIRE(str.find("3") != std::string::npos);
		REQUIRE(str.find("4") != std::string::npos);
		REQUIRE(str.find("5") != std::string::npos);
	}

	TEST_CASE("Matrix3D::StreamOutput", "[Matrix3D]") {
		Matrix3D<int> m(2, 3, 4);
		
		std::ostringstream oss;
		oss << m;
		
		REQUIRE(oss.str().find("Matrix3D") != std::string::npos);
	}

	// ============================================================================
	// Type Alias Tests
	// ============================================================================

	TEST_CASE("Matrix3D::TypeAliases", "[Matrix3D]") {
		Matrix3Dd md(2, 2, 2, 1.5);
		Matrix3Df mf(2, 2, 2, 1.5f);
		Matrix3Di mi(2, 2, 2, 42);
		Matrix3D<Real> mr(2, 2, 2, REAL(2.5));
		
		REQUIRE(md(0, 0, 0) == 1.5);
		REQUIRE(mf(0, 0, 0) == 1.5f);
		REQUIRE(mi(0, 0, 0) == 42);
		REQUIRE_THAT(mr(0, 0, 0), RealApprox(REAL(2.5)));
	}

	// ============================================================================
	// Memory Layout Tests
	// ============================================================================

	TEST_CASE("Matrix3D::ContiguousMemory", "[Matrix3D]") {
		Matrix3D<int> m(2, 3, 4, 0);
		
		// Fill with sequential values
		int val = 0;
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 3; ++j)
				for (int k = 0; k < 4; ++k)
					m(i, j, k) = val++;
		
		// Verify data is contiguous in memory
		int* data = m.data();
		for (int i = 0; i < 24; ++i)
			REQUIRE(data[i] == i);
	}
}
