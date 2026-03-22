#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Vector/Vector.h"
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

		REQUIRE_THROWS_AS(Vector<Real>::UnitVector(-1, 5), VectorDimensionError);
		REQUIRE_THROWS_AS(Vector<Real>::UnitVector(5, 7), VectorDimensionError);

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

	///////////////////////////////////////////////////////////////////////////
	//                    Complex Vector Tests                               //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("Vector<Complex>::NormL2_single_element", "[Vector][complex]")
	{
		TEST_PRECISION_INFO();
		// v = [3+4i], |3+4i|² = 9+16 = 25, ||v|| = 5
		Vector<Complex> v({Complex(REAL(3.0), REAL(4.0))});
		
		REQUIRE(v.size() == 1);
		REQUIRE_THAT(v.NormL2(), RealApprox(REAL(5.0)));
	}

	TEST_CASE("Vector<Complex>::NormL2_multiple_elements", "[Vector][complex]")
	{
		TEST_PRECISION_INFO();
		// v = [3+4i, 0], |3+4i|² = 25, |0|² = 0, ||v|| = 5
		Vector<Complex> v1({Complex(REAL(3.0), REAL(4.0)), Complex(REAL(0.0), REAL(0.0))});
		REQUIRE_THAT(v1.NormL2(), RealApprox(REAL(5.0)));
		
		// v = [1+1i, 1+1i], |1+1i|² = 2, ||v||² = 4, ||v|| = 2
		Vector<Complex> v2({Complex(REAL(1.0), REAL(1.0)), Complex(REAL(1.0), REAL(1.0))});
		REQUIRE_THAT(v2.NormL2(), RealApprox(REAL(2.0)));
		
		// v = [1, i, 1+i], |1|² = 1, |i|² = 1, |1+i|² = 2, ||v||² = 4, ||v|| = 2
		Vector<Complex> v3({Complex(REAL(1.0), REAL(0.0)), 
		                    Complex(REAL(0.0), REAL(1.0)), 
		                    Complex(REAL(1.0), REAL(1.0))});
		REQUIRE_THAT(v3.NormL2(), RealApprox(REAL(2.0)));
	}

	TEST_CASE("Vector<Complex>::NormL2_pure_real", "[Vector][complex]")
	{
		TEST_PRECISION_INFO();
		// Complex vector with zero imaginary parts should match real vector norm
		Vector<Complex> vc({Complex(REAL(3.0), REAL(0.0)), Complex(REAL(4.0), REAL(0.0))});
		Vector<Real> vr({REAL(3.0), REAL(4.0)});
		
		// Both should give ||v|| = 5
		REQUIRE_THAT(vc.NormL2(), RealApprox(REAL(5.0)));
		REQUIRE_THAT(vr.NormL2(), RealApprox(REAL(5.0)));
		REQUIRE_THAT(vc.NormL2(), RealApprox(vr.NormL2()));
	}

	TEST_CASE("Vector<Complex>::NormL2_pure_imaginary", "[Vector][complex]")
	{
		TEST_PRECISION_INFO();
		// v = [3i, 4i], |3i|² = 9, |4i|² = 16, ||v|| = 5
		Vector<Complex> v({Complex(REAL(0.0), REAL(3.0)), Complex(REAL(0.0), REAL(4.0))});
		REQUIRE_THAT(v.NormL2(), RealApprox(REAL(5.0)));
	}

	TEST_CASE("Vector<Complex>::NormL2_unit_circle", "[Vector][complex]")
	{
		TEST_PRECISION_INFO();
		// All elements on unit circle: |e^(iθ)| = 1
		// v = [1, i, -1, -i] (4 elements on unit circle), ||v|| = 2
		Vector<Complex> v({Complex(REAL(1.0), REAL(0.0)),   // 1
		                   Complex(REAL(0.0), REAL(1.0)),   // i
		                   Complex(REAL(-1.0), REAL(0.0)),  // -1
		                   Complex(REAL(0.0), REAL(-1.0))}); // -i
		REQUIRE_THAT(v.NormL2(), RealApprox(REAL(2.0)));  // sqrt(4) = 2
	}

	TEST_CASE("Vector<Complex>::NormL2_empty", "[Vector][complex][edge_cases]")
	{
		TEST_PRECISION_INFO();
		Vector<Complex> empty;
		REQUIRE(empty.size() == 0);
		REQUIRE(empty.NormL2() == REAL(0.0));
	}

	///////////////////////////////////////////////////////////////////////////
	//                    Copy/Move Semantics Tests                          //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("Vector::copy_constructor", "[Vector][copy_move]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> original({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> copy(original);
		
		REQUIRE(copy.size() == 3);
		REQUIRE(copy[0] == REAL(1.0));
		REQUIRE(copy[1] == REAL(2.0));
		REQUIRE(copy[2] == REAL(3.0));
		
		// Modify copy, original should be unchanged
		copy[0] = REAL(99.0);
		REQUIRE(original[0] == REAL(1.0));
	}

	TEST_CASE("Vector::copy_assignment", "[Vector][copy_move]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> original({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> copy;
		copy = original;
		
		REQUIRE(copy.size() == 3);
		REQUIRE(copy[0] == REAL(1.0));
		
		// Modify copy, original should be unchanged
		copy[1] = REAL(99.0);
		REQUIRE(original[1] == REAL(2.0));
	}

	TEST_CASE("Vector::move_constructor", "[Vector][copy_move]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> original({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> moved(std::move(original));
		
		REQUIRE(moved.size() == 3);
		REQUIRE(moved[0] == REAL(1.0));
		REQUIRE(moved[1] == REAL(2.0));
		REQUIRE(moved[2] == REAL(3.0));
		// original is in valid but unspecified state after move
	}

	TEST_CASE("Vector::move_assignment", "[Vector][copy_move]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> original({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> moved;
		moved = std::move(original);
		
		REQUIRE(moved.size() == 3);
		REQUIRE(moved[0] == REAL(1.0));
	}

	///////////////////////////////////////////////////////////////////////////
	//                    STL-like Interface Tests                           //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("Vector::isEmpty", "[Vector][stl_interface]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> empty;
		Vector<Real> nonempty({REAL(1.0)});
		
		REQUIRE(empty.empty() == true);
		REQUIRE(nonempty.empty() == false);
	}

	TEST_CASE("Vector::front_back", "[Vector][stl_interface]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> v({REAL(1.0), REAL(2.0), REAL(3.0)});
		
		REQUIRE(v.front() == REAL(1.0));
		REQUIRE(v.back() == REAL(3.0));
		
		// Modify through front/back
		v.front() = REAL(10.0);
		v.back() = REAL(30.0);
		REQUIRE(v[0] == REAL(10.0));
		REQUIRE(v[2] == REAL(30.0));
		
		// Const version
		const Vector<Real> cv({REAL(5.0), REAL(6.0), REAL(7.0)});
		REQUIRE(cv.front() == REAL(5.0));
		REQUIRE(cv.back() == REAL(7.0));
	}

	TEST_CASE("Vector::iterators", "[Vector][stl_interface]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> v({REAL(1.0), REAL(2.0), REAL(3.0)});
		
		// Range-based for (uses begin/end)
		Real sum = REAL(0.0);
		for (const auto& x : v) {
			sum += x;
		}
		REQUIRE(sum == REAL(6.0));
		
		// Modify through iterator
		for (auto& x : v) {
			x *= REAL(2.0);
		}
		REQUIRE(v[0] == REAL(2.0));
		REQUIRE(v[1] == REAL(4.0));
		REQUIRE(v[2] == REAL(6.0));
		
		// cbegin/cend
		const Vector<Real> cv({REAL(1.0), REAL(2.0)});
		Real csum = REAL(0.0);
		for (auto it = cv.cbegin(); it != cv.cend(); ++it) {
			csum += *it;
		}
		REQUIRE(csum == REAL(3.0));
	}

	TEST_CASE("Vector::push_back", "[Vector][stl_interface]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> v;
		
		v.push_back(REAL(1.0));
		REQUIRE(v.size() == 1);
		REQUIRE(v[0] == REAL(1.0));
		
		v.push_back(REAL(2.0));
		REQUIRE(v.size() == 2);
		REQUIRE(v[1] == REAL(2.0));
		
		// Push_back with move
		Real val = REAL(3.0);
		v.push_back(std::move(val));
		REQUIRE(v.size() == 3);
		REQUIRE(v[2] == REAL(3.0));
	}

	TEST_CASE("Vector::insert", "[Vector][stl_interface]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> v({REAL(1.0), REAL(3.0)});
		
		v.insert(1, REAL(2.0));  // Insert at position 1
		REQUIRE(v.size() == 3);
		REQUIRE(v[0] == REAL(1.0));
		REQUIRE(v[1] == REAL(2.0));
		REQUIRE(v[2] == REAL(3.0));
		
		v.insert(0, REAL(0.0));  // Insert at beginning
		REQUIRE(v.size() == 4);
		REQUIRE(v[0] == REAL(0.0));
		REQUIRE(v[1] == REAL(1.0));
	}

	TEST_CASE("Vector::erase", "[Vector][stl_interface]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> v({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0)});
		
		// Erase single element
		v.erase(1);  // Erase element at position 1
		REQUIRE(v.size() == 3);
		REQUIRE(v[0] == REAL(1.0));
		REQUIRE(v[1] == REAL(3.0));
		REQUIRE(v[2] == REAL(4.0));
		
		// Erase range
		Vector<Real> v2({REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0)});
		v2.erase(1, 3);  // Erase elements at positions 1,2
		REQUIRE(v2.size() == 3);
		REQUIRE(v2[0] == REAL(1.0));
		REQUIRE(v2[1] == REAL(4.0));
		REQUIRE(v2[2] == REAL(5.0));
		
		// Erase by value
		Vector<Real> v3({REAL(1.0), REAL(2.0), REAL(2.0), REAL(3.0)});
		v3.erase(REAL(2.0));  // Erase all occurrences of 2.0
		REQUIRE(v3.size() == 2);
		REQUIRE(v3[0] == REAL(1.0));
		REQUIRE(v3[1] == REAL(3.0));
	}

	TEST_CASE("Vector::Clear", "[Vector][stl_interface]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> v({REAL(1.0), REAL(2.0), REAL(3.0)});
		
		v.Clear();
		REQUIRE(v.size() == 0);
		REQUIRE(v.empty() == true);
	}

	TEST_CASE("Vector::Resize", "[Vector][stl_interface]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> v({REAL(1.0), REAL(2.0), REAL(3.0)});
		
		// Resize smaller without preserving
		v.Resize(2);
		REQUIRE(v.size() == 2);
		
		// Resize larger
		v.Resize(5);
		REQUIRE(v.size() == 5);
		
		// Resize with preserve
		Vector<Real> v2({REAL(1.0), REAL(2.0), REAL(3.0)});
		v2.Resize(5, true);  // Preserve existing elements
		REQUIRE(v2.size() == 5);
		REQUIRE(v2[0] == REAL(1.0));
		REQUIRE(v2[1] == REAL(2.0));
		REQUIRE(v2[2] == REAL(3.0));
		
		// Resize smaller with preserve
		Vector<Real> v3({REAL(1.0), REAL(2.0), REAL(3.0)});
		v3.Resize(2, true);
		REQUIRE(v3.size() == 2);
		REQUIRE(v3[0] == REAL(1.0));
		REQUIRE(v3[1] == REAL(2.0));
	}

	///////////////////////////////////////////////////////////////////////////
	//                    Compound Assignment Operators                      //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("Vector::operator_plus_equals", "[Vector][operators]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> a({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> b({REAL(4.0), REAL(5.0), REAL(6.0)});
		
		a += b;
		REQUIRE(a[0] == REAL(5.0));
		REQUIRE(a[1] == REAL(7.0));
		REQUIRE(a[2] == REAL(9.0));
		
		// b should be unchanged
		REQUIRE(b[0] == REAL(4.0));
	}

	TEST_CASE("Vector::operator_minus_equals", "[Vector][operators]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> a({REAL(5.0), REAL(7.0), REAL(9.0)});
		Vector<Real> b({REAL(1.0), REAL(2.0), REAL(3.0)});
		
		a -= b;
		REQUIRE(a[0] == REAL(4.0));
		REQUIRE(a[1] == REAL(5.0));
		REQUIRE(a[2] == REAL(6.0));
	}

	TEST_CASE("Vector::operator_times_equals", "[Vector][operators]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> a({REAL(1.0), REAL(2.0), REAL(3.0)});
		
		a *= REAL(2.0);
		REQUIRE(a[0] == REAL(2.0));
		REQUIRE(a[1] == REAL(4.0));
		REQUIRE(a[2] == REAL(6.0));
	}

	TEST_CASE("Vector::operator_divide_equals", "[Vector][operators]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> a({REAL(2.0), REAL(4.0), REAL(6.0)});
		
		a /= REAL(2.0);
		REQUIRE(a[0] == REAL(1.0));
		REQUIRE(a[1] == REAL(2.0));
		REQUIRE(a[2] == REAL(3.0));
	}

	TEST_CASE("Vector::unary_negation", "[Vector][operators]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> a({REAL(1.0), REAL(-2.0), REAL(3.0)});
		
		Vector<Real> neg = -a;
		REQUIRE(neg[0] == REAL(-1.0));
		REQUIRE(neg[1] == REAL(2.0));
		REQUIRE(neg[2] == REAL(-3.0));
		
		// Original unchanged
		REQUIRE(a[0] == REAL(1.0));
	}

	TEST_CASE("Vector::operator_not_equals", "[Vector][operators]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> a({REAL(1.0), REAL(2.0)});
		Vector<Real> b({REAL(1.0), REAL(2.0)});
		Vector<Real> c({REAL(1.0), REAL(3.0)});
		
		REQUIRE((a != b) == false);
		REQUIRE((a != c) == true);
	}

	///////////////////////////////////////////////////////////////////////////
	//                    Norm and Equality Tests                            //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("Vector::NormL1", "[Vector][norms]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> v({REAL(1.0), REAL(-2.0), REAL(3.0)});
		
		REQUIRE(v.NormL1() == REAL(6.0));  // |1| + |-2| + |3| = 6
		
		Vector<Real> zeros(3, REAL(0.0));
		REQUIRE(zeros.NormL1() == REAL(0.0));
	}

	TEST_CASE("Vector::NormLInf", "[Vector][norms]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> v({REAL(1.0), REAL(-5.0), REAL(3.0)});
		
		REQUIRE(v.NormLInf() == REAL(5.0));  // max(|1|, |-5|, |3|) = 5
		
		Vector<Real> zeros(3, REAL(0.0));
		REQUIRE(zeros.NormLInf() == REAL(0.0));
	}

	TEST_CASE("Vector::IsNullVec", "[Vector][norms]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> zeros(3, REAL(0.0));
		Vector<Real> nonzero({REAL(0.0), REAL(1.0), REAL(0.0)});
		
		REQUIRE(zeros.isZero() == true);
		REQUIRE(nonzero.isZero() == false);
		
		Vector<Real> empty;
		REQUIRE(empty.isZero() == true);  // Empty vector is null
	}

	///////////////////////////////////////////////////////////////////////////
	//                    Print Methods Tests                                //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("Vector::Print", "[Vector][output]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> v({REAL(1.0), REAL(2.0), REAL(3.0)});
		
		std::ostringstream oss;
		v.Print(oss, 5, 2);
		std::string result = oss.str();
		
		REQUIRE(result.find('[') != std::string::npos);
		REQUIRE(result.find(']') != std::string::npos);
		REQUIRE(result.find('1') != std::string::npos);
		REQUIRE(result.find('2') != std::string::npos);
		REQUIRE(result.find('3') != std::string::npos);
	}

	TEST_CASE("Vector::Print_zeroThreshold", "[Vector][output]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> v({REAL(1.0), REAL(1e-10), REAL(3.0)});
		
		std::ostringstream oss;
		v.Print(oss, 5, 2, REAL(1e-5));  // Values below 1e-5 printed as 0
		std::string result = oss.str();
		
		// The small value should be printed as 0
		REQUIRE(result.find('0') != std::string::npos);
	}

	TEST_CASE("Vector::PrintLine", "[Vector][output]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> v({REAL(1.0), REAL(2.0)});
		
		std::ostringstream oss;
		v.PrintLine(oss, "Values: ", 5, 2);
		std::string result = oss.str();
		
		REQUIRE(result.find("Values:") != std::string::npos);
		REQUIRE(result.find('\n') != std::string::npos);
	}

	TEST_CASE("Vector::PrintCol", "[Vector][output]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> v({REAL(1.0), REAL(2.0), REAL(3.0)});
		
		std::ostringstream oss;
		v.PrintCol(oss, 5, 2);
		std::string result = oss.str();
		
		// Should have multiple newlines (one per element)
		size_t newline_count = std::count(result.begin(), result.end(), '\n');
		REQUIRE(newline_count == 3);
	}

	TEST_CASE("Vector::stream_output", "[Vector][output]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> v({REAL(1.0), REAL(2.0)});
		
		std::ostringstream oss;
		oss << v;
		std::string result = oss.str();
		
		REQUIRE(result.find('[') != std::string::npos);
		REQUIRE(result.find(']') != std::string::npos);
	}

	///////////////////////////////////////////////////////////////////////////
	//                    Type Alias Tests                                   //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("Vector::type_aliases", "[Vector][types]")
	{
		TEST_PRECISION_INFO();
		// Verify type aliases compile and work correctly
		VectorInt vi({1, 2, 3});
		REQUIRE(vi.size() == 3);
		REQUIRE(vi[0] == 1);
		
		VectorFlt vf({1.0f, 2.0f, 3.0f});
		REQUIRE(vf.size() == 3);
		REQUIRE(vf[0] == 1.0f);
		
		VectorDbl vd({1.0, 2.0, 3.0});
		REQUIRE(vd.size() == 3);
		REQUIRE(vd[0] == 1.0);
		
		VectorComplex vc({Complex(1.0, 0.0), Complex(0.0, 1.0)});
		REQUIRE(vc.size() == 2);
		
		// Short aliases
		VecI vi2({4, 5, 6});
		REQUIRE(vi2.size() == 3);
		
		VecF vf2({4.0f, 5.0f});
		REQUIRE(vf2.size() == 2);
		
		VecD vd2({4.0, 5.0});
		REQUIRE(vd2.size() == 2);
		
		VecC vc2({Complex(1.0, 1.0)});
		REQUIRE(vc2.size() == 1);
	}

	///////////////////////////////////////////////////////////////////////////
	//                    UnitVector Additional Tests                     //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("Vector::UnitVector", "[Vector][factory]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> e0 = Vector<Real>::UnitVector(3, 0);
		REQUIRE(e0[0] == REAL(1.0));
		REQUIRE(e0[1] == REAL(0.0));
		REQUIRE(e0[2] == REAL(0.0));
		
		Vector<Real> e1 = Vector<Real>::UnitVector(3, 1);
		REQUIRE(e1[0] == REAL(0.0));
		REQUIRE(e1[1] == REAL(1.0));
		REQUIRE(e1[2] == REAL(0.0));
		
		Vector<Real> e2 = Vector<Real>::UnitVector(3, 2);
		REQUIRE(e2[0] == REAL(0.0));
		REQUIRE(e2[1] == REAL(0.0));
		REQUIRE(e2[2] == REAL(1.0));
	}

	///////////////////////////////////////////////////////////////////////////
	//                    Compound Assignment Dimension Mismatch             //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("Vector::compound_assignment_dimension_mismatch", "[Vector][exceptions]")
	{
		TEST_PRECISION_INFO();
		Vector<Real> a({REAL(1.0), REAL(2.0), REAL(3.0)});
		Vector<Real> b({REAL(1.0), REAL(2.0)});
		
		REQUIRE_THROWS_AS(a += b, VectorDimensionError);
		REQUIRE_THROWS_AS(a -= b, VectorDimensionError);
	}

} // namespace MML::VectorTests
