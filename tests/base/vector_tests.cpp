#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector.h"
#include "base/BaseUtils.h"
#endif

using namespace MML;
using namespace MML::Testing;
using namespace MML::Testing::Matchers;

namespace MML::Tests::Base::VectorTests
{
	TEST_CASE("Vector::Default_ctor_has_zero_size", "[Vector]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> b, c(0);
		REQUIRE(0 == b.size());
		REQUIRE(0 == c.size());
	}

	TEST_CASE("Vector::Ctor_init_elems_to_zero", "[Vector]") {
		TEST_PRECISION_INFO();
		Vector<Real> a(3);

		REQUIRE(3 == a.size());
		REQUIRE(0.0 == a[0]);
		REQUIRE(0.0 == a[1]);
		REQUIRE(0.0 == a[2]);

		Vector<Real> big(10000);

		REQUIRE(10000 == big.size());
		REQUIRE(0.0 == big[0]);
		REQUIRE(0.0 == big[5000]);
		REQUIRE(0.0 == big[7258]);
	}

	TEST_CASE("Vector::Ctor_init_to_value", "[Vector]") {
		TEST_PRECISION_INFO();
		Vector<Real> a(3, REAL(1.5));

		REQUIRE(3 == a.size());

		REQUIRE(REAL(1.5) == a[0]);
		REQUIRE(REAL(1.5) == a[1]);
		REQUIRE(REAL(1.5) == a[2]);
	}

	TEST_CASE("Vector::Ctor_init_from_real_array", "[simple]") {
		TEST_PRECISION_INFO();
		Real x[] = { REAL(1.0), REAL(2.0), REAL(3.0) };
		Vector<Real> a(3, x);

		REQUIRE(3 == a.size());

		REQUIRE(REAL(1.0) == a[0]);
		REQUIRE(REAL(2.0) == a[1]);
		REQUIRE(REAL(3.0) == a[2]);
	}

	TEST_CASE("Vector::Ctor_negative_size_throws", "[simple]")
	{
		TEST_PRECISION_INFO();
		REQUIRE_THROWS_AS(Vector<Real>(-5), VectorInitializationError);
		REQUIRE_THROWS_AS(Vector<Real>(-10, REAL(0.5)), VectorInitializationError);

		Real x[] = { REAL(1.0), REAL(2.0), REAL(3.0) };
		REQUIRE_THROWS_AS(Vector<Real>(-3, x), VectorInitializationError);
	}

	TEST_CASE("Vector::Ctor_init_from_std vector", "[simple]") {
		TEST_PRECISION_INFO();
		std::vector<Real> x{ REAL(1.0), REAL(2.0), REAL(3.0) };
		Vector<Real> a(x);

		REQUIRE(3 == a.size());

		REQUIRE(REAL(1.0) == a[0]);
		REQUIRE(REAL(2.0) == a[1]);
		REQUIRE(REAL(3.0) == a[2]);
	}

	TEST_CASE("Vector::Ctor_initializer_list", "[simple]") {
		TEST_PRECISION_INFO();
		Vector<Real> a{ REAL(1.0), REAL(2.0), REAL(3.0) };

		REQUIRE(3 == a.size());

		REQUIRE(REAL(1.0) == a[0]);
		REQUIRE(REAL(2.0) == a[1]);
		REQUIRE(REAL(3.0) == a[2]);
	}

	TEST_CASE("Vector::IsEqual", "[simple]") {
		TEST_PRECISION_INFO();
		Vector<Real> a({ REAL(1.0), REAL(2.0) });
		// Difference of 1e-5 is detectable across all precisions after scaling
		Vector<Real> b({ REAL(1.0), REAL(2.00001) });

		// Use precision-aware tolerance
		// 1e-4 scales to 1e-2 (float), 1e-8 scales to 1e-6 (float)
		REQUIRE(true == a.IsEqualTo(b, ScaleTolerance(REAL(1e-4))));
		REQUIRE(false == a.IsEqualTo(b, ScaleTolerance(REAL(1e-8))));
	}

	TEST_CASE("Utils::AreEqual", "[simple]") {
		TEST_PRECISION_INFO();
		Vector<Real> a({ REAL(1.0), REAL(2.0) });
		// Difference of 1e-5 is detectable across all precisions after scaling
		Vector<Real> b({ REAL(1.0), REAL(2.00001) });

		// Use precision-aware tolerance
		// 1e-4 scales to 1e-2 (float), 1e-8 scales to 1e-6 (float)
		REQUIRE(true == Utils::AreEqual(a, b, ScaleTolerance(REAL(1e-4))));
		REQUIRE(false == Utils::AreEqual(a, b, ScaleTolerance(REAL(1e-8))));
	}

	TEST_CASE("Vector::access_operators", "[simple]") {
		TEST_PRECISION_INFO();
		Vector<Real> a({ REAL(1.0), REAL(2.0), REAL(3.0) });

		REQUIRE(REAL(1.0) == a[0]);
		REQUIRE(REAL(2.0) == a[1]);
		REQUIRE(REAL(3.0) == a[2]);

		REQUIRE(REAL(1.0) == a.at(0));
		REQUIRE(REAL(2.0) == a.at(1));
		REQUIRE(REAL(3.0) == a.at(2));
	}

	TEST_CASE("Vector::template_type_deduction", "[simple]") {
		TEST_PRECISION_INFO();
		Vector vec_double({ 1.0, 2.0, 3.0 });

		REQUIRE(1.0 == vec_double[0]);
		REQUIRE(2.0 == vec_double[1]);
		REQUIRE(3.0 == vec_double[2]);

		Vector vec_float({ 1.0F, 2.0F, 3.0F });

		REQUIRE(1.0F == vec_float[0]);
		REQUIRE(2.0F == vec_float[1]);
		REQUIRE(3.0F == vec_float[2]);

		Vector b_int({ 1, 2, 3 });

		REQUIRE(1 == b_int[0]);
		REQUIRE(2 == b_int[1]);
		REQUIRE(3 == b_int[2]);

		Vector vec_cmplx({ Complex(1,1), Complex(0,2) });

		REQUIRE(Complex(1, 1) == vec_cmplx[0]);
		REQUIRE(Complex(0, 2) == vec_cmplx[1]);

		int x = 5;
		Vector<Real> vec_dbl_1({ REAL(x * 1.0), REAL(sin(2.0)), REAL(3.0) });
		// CAN'T DO THIS!
		// Vector vec_dbl_2({ 1 , sin(2.0), 3.0 });
		// needs explicit template type argument
		Vector<Real> vec_dbl_2({ REAL(1) , REAL(sin(2.0)), REAL(3.0) });

		Vector<Complex> vec_cmplx_1({ REAL(3.0) * Complex(1,1), Complex(0,2) - REAL(sin(2.5)) });
		Vector<Complex> vec_cmplx_2({ REAL(1), REAL(3.5678), Complex(1, 1), REAL(sin(2.5)) });
	}

	TEST_CASE("Vector::Op+-", "[simple]") {
		TEST_PRECISION_INFO();
		Vector<Real> a({ REAL(1.0), REAL(2.0) });
		Vector<Real> b({ REAL(1.0), REAL(2.0) });

		auto c = a + b;
		auto d = a - b;

		REQUIRE(REAL(2.0) == c[0]);
		REQUIRE(REAL(4.0) == c[1]);

		REQUIRE(REAL(0.0) == d[0]);
		REQUIRE(REAL(0.0) == d[1]);
	}

	TEST_CASE("Vector::mul_double", "[simple]") {
		TEST_PRECISION_INFO();
		Vector a({ REAL(1.0), REAL(100.0) });

		auto b = a * REAL(2.0);
		auto c = REAL(2.0) * a;

		REQUIRE(REAL(2.0) == b[0]);
		REQUIRE(REAL(2.0) == c[0]);

		REQUIRE(REAL(200.0) == b[1]);
		REQUIRE(REAL(200.0) == c[1]);
	}

	TEST_CASE("Vector::div_double", "[simple]") {
		TEST_PRECISION_INFO();
		Vector a({ REAL(4.0), REAL(400.0) });

		auto b = a / REAL(2.0);

		REQUIRE(REAL(2.0) == b[0]);
		REQUIRE(REAL(200.0) == b[1]);
	}

	TEST_CASE("Utils::ScalarProduct", "[simple]") {
		TEST_PRECISION_INFO();
		Vector a({ REAL(1.0), REAL(2.0) });
		Vector b({ REAL(1.0), REAL(2.0) });

		REQUIRE(REAL(5.0) == Utils::ScalarProduct(a, b));
	}

	TEST_CASE("Vector::NormL2", "[simple]") {
		TEST_PRECISION_INFO();
		Vector a({ REAL(2.0), REAL(2.0) });

		REQUIRE_THAT(a.NormL2(), RealApprox(std::sqrt(REAL(8.0))));
	}

	//TEST_CASE("Vector::AngleToVector", "[simple]") {
	//	Vector a({ 0.0, 0.0, 5.0 });
	//	Vector b({ 1.0, 0.0, 0.0 });

	//	REQUIRE(Constants::PI / 2 == a.AngleToVector(b));

	//	Vector c({ 0.0, 1.0 });
	//	Vector d({ 1.0, 1.0 });

	//	REQUIRE_THAT(Constants::PI / 4, Catch::Matchers::WithinAbs(c.AngleToVector(d), REAL(1e-16)));
	//}

	TEST_CASE("Vector::to_string", "[simple]") {
		TEST_PRECISION_INFO();
		Vector a({ REAL(2.0), REAL(2.0) });

		REQUIRE("[    2,     2]" == a.to_string(5, 3));

		Vector b({ REAL(123.0), REAL(1.0), REAL(10.0), REAL(-8.0) });

		REQUIRE("[    123,       1,      10,      -8]" == b.to_string(7, 3));

		// Use float-safe values (powers of 2) for precision-independent testing
		Vector c({ REAL(128.0), REAL(2.0), REAL(16.0) });

		REQUIRE("[       128,          2,         16]" == c.to_string(10, 5));
	REQUIRE("[            128,               2,              16]" == c.to_string(15, 9));
	}

	TEST_CASE("Vector::exceptions", "[simple]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> vec_dim_5(5);
		Vector<Real> vec_dim_4({ REAL(1.0), REAL(0.0), REAL(0.0), REAL(1.0) });

		REQUIRE_THROWS_AS(Vector<Real>::GetUnitVector(-1, 5), VectorDimensionError);
		REQUIRE_THROWS_AS(Vector<Real>::GetUnitVector(5, 7), VectorDimensionError);

		REQUIRE_THROWS_AS(vec_dim_4.IsEqualTo(vec_dim_5), VectorDimensionError);

		REQUIRE_THROWS_AS(vec_dim_4.at(5), VectorDimensionError);
		REQUIRE_THROWS_AS(vec_dim_4.at(-1), VectorDimensionError);
		REQUIRE_THROWS_AS(vec_dim_4.at(4) = REAL(5.0), VectorDimensionError);
		REQUIRE_THROWS_AS(vec_dim_4.at(-1) = REAL(5.0), VectorDimensionError);

		REQUIRE_THROWS_AS(vec_dim_4 + vec_dim_5, VectorDimensionError);
		REQUIRE_THROWS_AS(vec_dim_4 - vec_dim_5, VectorDimensionError);
		REQUIRE_THROWS_AS(vec_dim_4 == vec_dim_5, VectorDimensionError);

		REQUIRE_THROWS_AS(Utils::ScalarProduct(vec_dim_4, vec_dim_5), VectorDimensionError);
		REQUIRE_THROWS_AS(Utils::VectorsAngle(vec_dim_4, vec_dim_5), VectorDimensionError);
	}

	TEST_CASE("Vector::NonNumericType", "[type_safety]")
	{
		TEST_PRECISION_INFO();
		// Test that Vector works with non-numeric types (e.g., for future statistics use)
		struct TestStruct {
			int id;
			std::string name;
			TestStruct() : id(0), name("") {}
			TestStruct(int i, std::string n) : id(i), name(n) {}
			bool operator==(const TestStruct& other) const {
				return id == other.id && name == other.name;
			}
		};

		// Default construction should work
		Vector<TestStruct> vec1(5);
		REQUIRE(vec1.size() == 5);
		REQUIRE(vec1[0].id == 0);
		REQUIRE(vec1[0].name == "");

		// Initialization with value should work
		TestStruct val{42, "test"};
		Vector<TestStruct> vec2(3, val);
		REQUIRE(vec2.size() == 3);
		REQUIRE(vec2[0].id == 42);
		REQUIRE(vec2[1].name == "test");
		REQUIRE(vec2[2] == val);

		// String vector for statistics
		Vector<std::string> str_vec(3);
		REQUIRE(str_vec.size() == 3);
		REQUIRE(str_vec[0] == "");
		
		Vector<std::string> str_vec2(2, "hello");
		REQUIRE(str_vec2.size() == 2);
		REQUIRE(str_vec2[0] == "hello");
		REQUIRE(str_vec2[1] == "hello");

		// Initializer list
		Vector<std::string> str_vec3{"alpha", "beta", "gamma"};
		REQUIRE(str_vec3.size() == 3);
		REQUIRE(str_vec3[0] == "alpha");
		REQUIRE(str_vec3[2] == "gamma");
	}

	// ========== EDGE CASE TESTS ==========

	TEST_CASE("Vector::EdgeCase_empty_vector_operations", "[Vector][edge_cases]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> empty1, empty2;
		
		// Empty vectors should be equal
		REQUIRE(empty1.IsEqualTo(empty2));
		
		// Arithmetic on empty vectors
		Vector<Real> sum = empty1 + empty2;
		REQUIRE(sum.size() == 0);
		
		Vector<Real> diff = empty1 - empty2;
		REQUIRE(diff.size() == 0);
		
		// Scalar operations on empty vector
		Vector<Real> scaled = empty1 * 5.0;
		REQUIRE(scaled.size() == 0);
		
		// NormL2 of empty vector should be 0
		REQUIRE(empty1.NormL2() == 0.0);
	}

	TEST_CASE("Vector::EdgeCase_single_element", "[Vector][edge_cases]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> single1({REAL(3.0)});
		Vector<Real> single2({REAL(4.0)});
		
		REQUIRE(single1.size() == 1);
		REQUIRE(single1[0] == REAL(3.0));
		
		// Arithmetic
		Vector<Real> sum = single1 + single2;
		REQUIRE(sum.size() == 1);
		REQUIRE(sum[0] == REAL(7.0));
		
		// Norm
		REQUIRE(single1.NormL2() == REAL(3.0));
		
		// Scalar product of single-element vectors
		Real dot = Utils::ScalarProduct(single1, single2);
		REQUIRE(dot == REAL(12.0));
	}

	TEST_CASE("Vector::EdgeCase_large_vector", "[Vector][edge_cases]")
	{
		TEST_PRECISION_INFO();
		const int N = 100000;
		Vector<Real> large(N, REAL(1.0));
		
		REQUIRE(large.size() == N);
		REQUIRE(large[0] == REAL(1.0));
		REQUIRE(large[N-1] == REAL(1.0));
		REQUIRE(large[N/2] == REAL(1.0));
		
		// NormL2 of vector of ones should be sqrt(N)
		Real norm = large.NormL2();
		REQUIRE_THAT(norm, RealApprox(std::sqrt(Real(N))));
	}

} // namespace MML::VectorTests