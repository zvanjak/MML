#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Vector/VectorN.h>
#include <mml/base/BaseUtils.h>
#endif

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Base::VectorNTests
{
	TEST_CASE("VectorN::default_ctor_init_to_zero", "[VectorN]") {
			TEST_PRECISION_INFO();
		VectorN<Real, 3> a;

		REQUIRE(static_cast<MML::Real>(0.0) == a[0]);
		REQUIRE(static_cast<MML::Real>(0.0) == a[1]);
		REQUIRE(static_cast<MML::Real>(0.0) == a[2]);
	}

	TEST_CASE("VectorN::initializer_list_ctor", "[VectorN]") {
			TEST_PRECISION_INFO();
		VectorN<Real, 3> a{ static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0) };

		REQUIRE(static_cast<MML::Real>(1.0) == a[0]);
		REQUIRE(static_cast<MML::Real>(2.0) == a[1]);
		REQUIRE(static_cast<MML::Real>(3.0) == a[2]);
	}

	TEST_CASE("VectorN::init_to_value", "[VectorN]") {
			TEST_PRECISION_INFO();
		VectorN<Real, 3> a(static_cast<MML::Real>(1.5));

		REQUIRE(3 == a.size());

		REQUIRE(static_cast<MML::Real>(1.5) == a[0]);
		REQUIRE(static_cast<MML::Real>(1.5) == a[1]);
		REQUIRE(static_cast<MML::Real>(1.5) == a[2]);
	}

	TEST_CASE("VectorN::init_from_ std vector", "[VectorN]") {
			TEST_PRECISION_INFO();
		std::vector<Real> x{ static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0) };
		VectorN<Real, 3> a(x);

		REQUIRE(3 == a.size());

		REQUIRE(static_cast<MML::Real>(1.0) == a[0]);
		REQUIRE(static_cast<MML::Real>(2.0) == a[1]);
		REQUIRE(static_cast<MML::Real>(3.0) == a[2]);
	}

	TEST_CASE("VectorN::init_from_double_array", "[VectorN]") {
			TEST_PRECISION_INFO();
		Real x[] = { static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0) };
		VectorN<Real, 3> a(x);

		REQUIRE(3 == a.size());

		REQUIRE(static_cast<MML::Real>(1.0) == a[0]);
		REQUIRE(static_cast<MML::Real>(2.0) == a[1]);
		REQUIRE(static_cast<MML::Real>(3.0) == a[2]);
	}

	TEST_CASE("VectorN::IsEqual", "[VectorN]") {
			TEST_PRECISION_INFO();
		Vec2Dbl a({ static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0) });
		// Difference of 1e-5 is detectable across all precisions after scaling
		Vec2Dbl b({ static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.00001) });

		// 1e-4 scales to 1e-2 (float), TOL(1e-8, 1e-4) scales to 1e-6 (float)
		REQUIRE(true == a.IsEqualTo(b, TOL3(1e-4, 1e-2, 1e-4)));
		REQUIRE(false == a.IsEqualTo(b, ScaleTolerance(TOL(1e-8, 1e-8))));
	}

	TEST_CASE("VectorN::Op+-", "[VectorN]") {
			TEST_PRECISION_INFO();
		Vec2Dbl a({ static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0) });
		Vec2Dbl b({ static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0) });

		auto c = a + b;
		auto d = a - b;

		REQUIRE(static_cast<MML::Real>(2.0) == c[0]);
		REQUIRE(static_cast<MML::Real>(4.0) == c[1]);

		REQUIRE(static_cast<MML::Real>(0.0) == d[0]);
		REQUIRE(static_cast<MML::Real>(0.0) == d[1]);
	}

	TEST_CASE("VectorN::mul_double", "[VectorN]") {
			TEST_PRECISION_INFO();
		Vec2Dbl a({ static_cast<MML::Real>(1.0), static_cast<MML::Real>(100.0) });

		auto b = a * static_cast<MML::Real>(2.0);
		auto c = static_cast<MML::Real>(2.0) * a;

		REQUIRE(static_cast<MML::Real>(2.0) == b[0]);
		REQUIRE(static_cast<MML::Real>(2.0) == c[0]);

		REQUIRE(static_cast<MML::Real>(200.0) == b[1]);
		REQUIRE(static_cast<MML::Real>(200.0) == c[1]);
	}

	TEST_CASE("VectorN::div_double", "[VectorN]") {
			TEST_PRECISION_INFO();
		Vec2Dbl a({ static_cast<MML::Real>(4.0), static_cast<MML::Real>(400.0) });

		auto b = a / static_cast<MML::Real>(2.0);

		REQUIRE(static_cast<MML::Real>(2.0) == b[0]);
		REQUIRE(static_cast<MML::Real>(200.0) == b[1]);
	}

	TEST_CASE("VectorN::NormL2", "[VectorN]") {
			TEST_PRECISION_INFO();
		Vec2Dbl a({ static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0) });

		REQUIRE_THAT(a.NormL2(), RealApprox(std::sqrt(static_cast<MML::Real>(8.0))));
	}

	TEST_CASE("VectorN::to_string", "[VectorN]") {
			TEST_PRECISION_INFO();
		Vec2Dbl a({ static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0) });

		REQUIRE("[    2,     2]" == a.to_string(5, 3));

		Vec4Dbl b({ static_cast<MML::Real>(123.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(10.0), -static_cast<MML::Real>(8.0) });

		REQUIRE("[    123,       1,      10,      -8]" == b.to_string(7, 3));

		Vec3Dbl c({ static_cast<MML::Real>(123.123), static_cast<MML::Real>(1.9876543), static_cast<MML::Real>(10.0) });

	REQUIRE("[    123.12,     1.9877,         10]" == c.to_string(10, 5));
	if constexpr (std::is_same_v<Real, float>) {
		// Float has ~7 significant digits; 9-digit formatting shows precision artifacts
	} else {
		REQUIRE("[        123.123,       1.9876543,              10]" == c.to_string(15, 9));
	}
	}

	TEST_CASE("VectorN::exceptions", "[VectorN]")
	{
			TEST_PRECISION_INFO();
		Vec2Dbl vec2({ static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0) });
		REQUIRE_THROWS_AS(vec2.at(5), VectorAccessBoundsError);
		REQUIRE_THROWS_AS(vec2.at(-1), VectorAccessBoundsError);

		REQUIRE_THROWS_AS(vec2.at(4) = static_cast<MML::Real>(5.0), VectorAccessBoundsError);
		REQUIRE_THROWS_AS(vec2.at(-1) = static_cast<MML::Real>(5.0), VectorAccessBoundsError);
	}

	///////////////////////////////////////////////////////////////////////////
	//                    Factory Methods Tests                              //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("VectorN::UnitVector", "[VectorN][factory]")
	{
		TEST_PRECISION_INFO();
		auto e0 = VectorN<Real, 3>::UnitVector(0);
		REQUIRE(e0[0] == static_cast<MML::Real>(1.0));
		REQUIRE(e0[1] == static_cast<MML::Real>(0.0));
		REQUIRE(e0[2] == static_cast<MML::Real>(0.0));
		
		auto e1 = VectorN<Real, 3>::UnitVector(1);
		REQUIRE(e1[0] == static_cast<MML::Real>(0.0));
		REQUIRE(e1[1] == static_cast<MML::Real>(1.0));
		REQUIRE(e1[2] == static_cast<MML::Real>(0.0));
		
		auto e2 = VectorN<Real, 3>::UnitVector(2);
		REQUIRE(e2[0] == static_cast<MML::Real>(0.0));
		REQUIRE(e2[1] == static_cast<MML::Real>(0.0));
		REQUIRE(e2[2] == static_cast<MML::Real>(1.0));
	}

	TEST_CASE("VectorN::Normalized", "[VectorN][factory]")
	{
		TEST_PRECISION_INFO();
		using namespace MML::Testing::Matchers;
		
		Vec3 v({static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(0.0)});  // length = 5
		Vec3 normalized = v.Normalized();
		
		REQUIRE_THAT(normalized[0], RealApprox(static_cast<MML::Real>(0.6)));
		REQUIRE_THAT(normalized[1], RealApprox(static_cast<MML::Real>(0.8)));
		REQUIRE_THAT(normalized[2], RealApprox(static_cast<MML::Real>(0.0)));
		REQUIRE_THAT(normalized.NormL2(), RealApprox(static_cast<MML::Real>(1.0)));
		
		// Zero vector throws
		Vec3 zero;
		REQUIRE_THROWS_AS(zero.Normalized(), VectorDimensionError);
	}

	///////////////////////////////////////////////////////////////////////////
	//                    Size and Clear Tests                               //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("VectorN::size_and_clear", "[VectorN][basic]")
	{
		TEST_PRECISION_INFO();
		Vec2 v2;
		Vec3 v3;
		Vec4 v4;
		
		REQUIRE(v2.size() == 2);
		REQUIRE(v3.size() == 3);
		REQUIRE(v4.size() == 4);
		
		VectorN<Real, 10> v10({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)});
		REQUIRE(v10.size() == 10);
		REQUIRE(v10[0] == static_cast<MML::Real>(1.0));
		
		v10.SetZero();
		REQUIRE(v10[0] == static_cast<MML::Real>(0.0));
		REQUIRE(v10[1] == static_cast<MML::Real>(0.0));
		REQUIRE(v10[2] == static_cast<MML::Real>(0.0));
	}

	///////////////////////////////////////////////////////////////////////////
	//                    Equality Operators Tests                           //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("VectorN::operator_equals", "[VectorN][operators]")
	{
		TEST_PRECISION_INFO();
		Vec3 a({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)});
		Vec3 b({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)});
		Vec3 c({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(4.0)});
		
		REQUIRE((a == b) == true);
		REQUIRE((a == c) == false);
	}

	TEST_CASE("VectorN::operator_not_equals", "[VectorN][operators]")
	{
		TEST_PRECISION_INFO();
		Vec3 a({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)});
		Vec3 b({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)});
		Vec3 c({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(4.0)});
		
		REQUIRE((a != b) == false);
		REQUIRE((a != c) == true);
	}

	TEST_CASE("VectorN::AreEqual_static", "[VectorN][operators]")
	{
		TEST_PRECISION_INFO();
		Vec3 a({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)});
		Vec3 b({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.00001), static_cast<MML::Real>(3.0)});
		
		REQUIRE(Vec3::AreEqual(a, b, TOL3(1e-4, 1e-2, 1e-4)) == true);
		REQUIRE(Vec3::AreEqual(a, b, ScaleTolerance(TOL(1e-8, 1e-8))) == false);
	}

	TEST_CASE("VectorN::IsNullVec", "[VectorN][operators]")
	{
		TEST_PRECISION_INFO();
		Vec3 zero;
		Vec3 nonzero({static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0)});
		
		REQUIRE(zero.IsZero() == true);
		REQUIRE(nonzero.IsZero() == false);
	}

	///////////////////////////////////////////////////////////////////////////
	//                    Compound Assignment Operators                      //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("VectorN::operator_plus_equals", "[VectorN][operators]")
	{
		TEST_PRECISION_INFO();
		Vec3 a({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)});
		Vec3 b({static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0), static_cast<MML::Real>(6.0)});
		
		a += b;
		REQUIRE(a[0] == static_cast<MML::Real>(5.0));
		REQUIRE(a[1] == static_cast<MML::Real>(7.0));
		REQUIRE(a[2] == static_cast<MML::Real>(9.0));
	}

	TEST_CASE("VectorN::operator_minus_equals", "[VectorN][operators]")
	{
		TEST_PRECISION_INFO();
		Vec3 a({static_cast<MML::Real>(5.0), static_cast<MML::Real>(7.0), static_cast<MML::Real>(9.0)});
		Vec3 b({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)});
		
		a -= b;
		REQUIRE(a[0] == static_cast<MML::Real>(4.0));
		REQUIRE(a[1] == static_cast<MML::Real>(5.0));
		REQUIRE(a[2] == static_cast<MML::Real>(6.0));
	}

	TEST_CASE("VectorN::operator_times_equals", "[VectorN][operators]")
	{
		TEST_PRECISION_INFO();
		Vec3 a({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)});
		
		a *= static_cast<MML::Real>(2.0);
		REQUIRE(a[0] == static_cast<MML::Real>(2.0));
		REQUIRE(a[1] == static_cast<MML::Real>(4.0));
		REQUIRE(a[2] == static_cast<MML::Real>(6.0));
	}

	TEST_CASE("VectorN::operator_divide_equals", "[VectorN][operators]")
	{
		TEST_PRECISION_INFO();
		Vec3 a({static_cast<MML::Real>(2.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(6.0)});
		
		a /= static_cast<MML::Real>(2.0);
		REQUIRE(a[0] == static_cast<MML::Real>(1.0));
		REQUIRE(a[1] == static_cast<MML::Real>(2.0));
		REQUIRE(a[2] == static_cast<MML::Real>(3.0));
	}

	TEST_CASE("VectorN::unary_negation", "[VectorN][operators]")
	{
		TEST_PRECISION_INFO();
		Vec3 a({static_cast<MML::Real>(1.0), static_cast<MML::Real>(-2.0), static_cast<MML::Real>(3.0)});
		
		Vec3 neg = -a;
		REQUIRE(neg[0] == static_cast<MML::Real>(-1.0));
		REQUIRE(neg[1] == static_cast<MML::Real>(2.0));
		REQUIRE(neg[2] == static_cast<MML::Real>(-3.0));
		
		// Original unchanged
		REQUIRE(a[0] == static_cast<MML::Real>(1.0));
	}

	///////////////////////////////////////////////////////////////////////////
	//                    Norm Tests                                         //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("VectorN::NormL1", "[VectorN][norms]")
	{
		TEST_PRECISION_INFO();
		Vec3 v({static_cast<MML::Real>(1.0), static_cast<MML::Real>(-2.0), static_cast<MML::Real>(3.0)});
		
		REQUIRE(v.NormL1() == static_cast<MML::Real>(6.0));  // |1| + |-2| + |3| = 6
	}

	TEST_CASE("VectorN::NormLInf", "[VectorN][norms]")
	{
		TEST_PRECISION_INFO();
		Vec3 v({static_cast<MML::Real>(1.0), static_cast<MML::Real>(-5.0), static_cast<MML::Real>(3.0)});
		
		REQUIRE(v.NormLInf() == static_cast<MML::Real>(5.0));  // max(|1|, |-5|, |3|) = 5
	}

	///////////////////////////////////////////////////////////////////////////
	//                    Output Tests                                       //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("VectorN::Print", "[VectorN][output]")
	{
		TEST_PRECISION_INFO();
		Vec3 v({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)});
		
		std::ostringstream oss;
		v.Print(oss, 5, 2);
		std::string result = oss.str();
		
		REQUIRE(result.find('[') != std::string::npos);
		REQUIRE(result.find(']') != std::string::npos);
	}

	TEST_CASE("VectorN::Print_zeroThreshold", "[VectorN][output]")
	{
		TEST_PRECISION_INFO();
		Vec3 v({static_cast<MML::Real>(1.0), TOL(1e-10, 1e-5), static_cast<MML::Real>(3.0)});
		
		std::ostringstream oss;
		v.Print(oss, 5, 2, static_cast<MML::Real>(1e-5));
		std::string result = oss.str();
		
		// The small value should be printed as 0
		REQUIRE(result.find('0') != std::string::npos);
	}

	TEST_CASE("VectorN::PrintLine", "[VectorN][output]")
	{
		TEST_PRECISION_INFO();
		Vec3 v({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)});
		
		std::ostringstream oss;
		v.PrintLine(oss, "Values: ", 5, 2);
		std::string result = oss.str();
		
		REQUIRE(result.find("Values:") != std::string::npos);
		REQUIRE(result.find('\n') != std::string::npos);
	}

	TEST_CASE("VectorN::stream_output", "[VectorN][output]")
	{
		TEST_PRECISION_INFO();
		Vec3 v({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)});
		
		std::ostringstream oss;
		oss << v;
		std::string result = oss.str();
		
		REQUIRE(result.find('[') != std::string::npos);
		REQUIRE(result.find(']') != std::string::npos);
	}

	///////////////////////////////////////////////////////////////////////////
	//                    Type Alias Tests                                   //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("VectorN::type_aliases", "[VectorN][types]")
	{
		TEST_PRECISION_INFO();
		// Real-based aliases
		Vec2 v2({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)});
		REQUIRE(v2.size() == 2);
		
		Vec3 v3({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0)});
		REQUIRE(v3.size() == 3);
		
		Vec4 v4({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0)});
		REQUIRE(v4.size() == 4);
		
		// Float aliases
		Vec2Flt vf2({1.0f, 2.0f});
		REQUIRE(vf2.size() == 2);
		REQUIRE(vf2[0] == 1.0f);
		
		Vec3Flt vf3({1.0f, 2.0f, 3.0f});
		REQUIRE(vf3.size() == 3);
		
		Vec4Flt vf4({1.0f, 2.0f, 3.0f, 4.0f});
		REQUIRE(vf4.size() == 4);
		
		// Double aliases
		Vec2Dbl vd2({1.0, 2.0});
		REQUIRE(vd2.size() == 2);
		
		Vec3Dbl vd3({1.0, 2.0, 3.0});
		REQUIRE(vd3.size() == 3);
		
		Vec4Dbl vd4({1.0, 2.0, 3.0, 4.0});
		REQUIRE(vd4.size() == 4);
		
		// Short aliases
		Vec2F v2f({1.0f, 2.0f});
		REQUIRE(v2f.size() == 2);
		
		Vec3D v3d({1.0, 2.0, 3.0});
		REQUIRE(v3d.size() == 3);
		
		// Complex aliases
		Vec2Complex vc2({Complex(1.0, 0.0), Complex(0.0, 1.0)});
		REQUIRE(vc2.size() == 2);
		
		Vec3C vc3({Complex(1.0, 1.0), Complex(2.0, 2.0), Complex(3.0, 3.0)});
		REQUIRE(vc3.size() == 3);
	}

	///////////////////////////////////////////////////////////////////////////
	//                    Complex Vector Tests                               //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("VectorN<Complex>::NormL2", "[VectorN][complex]")
	{
		TEST_PRECISION_INFO();
		using namespace MML::Testing::Matchers;
		
		// v = [3+4i, 0], |3+4i|² = 25, ||v|| = 5
		Vec2C v({Complex(static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0)), Complex(static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0))});
		REQUIRE_THAT(v.NormL2(), RealApprox(static_cast<MML::Real>(5.0)));
		
		// v = [1+1i, 1+1i], |1+1i|² = 2, ||v||² = 4, ||v|| = 2
		Vec2C v2({Complex(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)), Complex(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0))});
		REQUIRE_THAT(v2.NormL2(), RealApprox(static_cast<MML::Real>(2.0)));
	}

	TEST_CASE("VectorN<Complex>::operations", "[VectorN][complex]")
	{
		TEST_PRECISION_INFO();
		Vec2C a({Complex(static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0)), Complex(static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0))});
		Vec2C b({Complex(static_cast<MML::Real>(3.0), static_cast<MML::Real>(3.0)), Complex(static_cast<MML::Real>(4.0), static_cast<MML::Real>(4.0))});
		
		// Addition
		Vec2C sum = a + b;
		REQUIRE(sum[0] == Complex(static_cast<MML::Real>(4.0), static_cast<MML::Real>(4.0)));
		REQUIRE(sum[1] == Complex(static_cast<MML::Real>(6.0), static_cast<MML::Real>(6.0)));
		
		// Subtraction
		Vec2C diff = b - a;
		REQUIRE(diff[0] == Complex(static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0)));
		REQUIRE(diff[1] == Complex(static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0)));
		
		// Scalar multiplication
		Vec2C scaled = a * Complex(static_cast<MML::Real>(2.0), static_cast<MML::Real>(0.0));
		REQUIRE(scaled[0] == Complex(static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0)));
	}

	///////////////////////////////////////////////////////////////////////////
	//                    Edge Cases                                         //
	///////////////////////////////////////////////////////////////////////////

	TEST_CASE("VectorN::initializer_list_truncation", "[VectorN][edge_cases]")
	{
		TEST_PRECISION_INFO();
		// More elements than N - should throw VectorDimensionError
		REQUIRE_THROWS_AS(Vec2({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0)}), VectorDimensionError);
	}

	TEST_CASE("VectorN::initializer_list_fewer", "[VectorN][edge_cases]")
	{
		TEST_PRECISION_INFO();
		// Fewer elements than N - remaining should stay at default (0)
		Vec4 v({static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)});
		REQUIRE(v.size() == 4);
		REQUIRE(v[0] == static_cast<MML::Real>(1.0));
		REQUIRE(v[1] == static_cast<MML::Real>(2.0));
		REQUIRE(v[2] == static_cast<MML::Real>(0.0));
		REQUIRE(v[3] == static_cast<MML::Real>(0.0));
	}

} // namespace MML::Tests::VectorNTests
