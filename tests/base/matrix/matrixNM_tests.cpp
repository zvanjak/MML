#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"
#include "../../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Matrix/MatrixNM.h"
#endif

using namespace MML;
using namespace MML::Testing;

namespace MML::Tests::Base::MatrixNMTests
{
	TEST_CASE("MatrixNM::default_ctor_init_to_zero", "[simple]") {
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a;

		REQUIRE(REAL(0.0) == a(0, 0));
		REQUIRE(REAL(0.0) == a(0, 1));
		REQUIRE(REAL(0.0) == a(1, 0));
		REQUIRE(REAL(0.0) == a(1, 1));
	}

	TEST_CASE("MatrixNM::initializer_list_ctor", "[simple]") {
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		REQUIRE(2 == a.rows());
		REQUIRE(2 == a.cols());

		REQUIRE(REAL(1.0) == a(0, 0));
		REQUIRE(REAL(2.0) == a(0, 1));
		REQUIRE(REAL(3.0) == a(1, 0));
		REQUIRE(REAL(4.0) == a(1, 1));
	}

	TEST_CASE("MatrixNM::MakeUnitMatrix", "[simple]") {
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a;

		a.MakeUnitMatrix();

		REQUIRE(REAL(1.0) == a(0, 0));
		REQUIRE(REAL(0.0) == a(0, 1));
		REQUIRE(REAL(0.0) == a(1, 0));
		REQUIRE(REAL(1.0) == a(1, 1));
	}

	TEST_CASE("MatrixNM::Identity", "[simple]") {
			TEST_PRECISION_INFO();
		auto a = MatrixNM<Real, 2, 2>::Identity();

		REQUIRE(REAL(1.0) == a[0][0]);
		REQUIRE(REAL(0.0) == a[0][1]);
		REQUIRE(REAL(0.0) == a[1][0]);
		REQUIRE(REAL(1.0) == a[1][1]);
	}

	TEST_CASE("MatrixNM::IsEqual", "[simple]") {
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		MatrixNM<Real, 2, 2> b({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		REQUIRE(true == a.IsEqualTo(b));
	}
	TEST_CASE("MatrixNM::IsEqual2", "[simple]") {
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		MatrixNM<Real, 2, 2> b({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(5.0) });

		REQUIRE(false == a.IsEqualTo(b));
	}
	TEST_CASE("MatrixNM::IsEqual3", "[simple]") {
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		MatrixNM<Real, 2, 2> b({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0001) });

		// Difference is 1e-4, use 1e-3 (yes) and 1e-7 (no)
		// 1e-3 scales to 1e-1 (float), 1e-7 scales to 1e-5 (float)
		REQUIRE(true == a.IsEqualTo(b, TOL3(1e-3, 1e-1, 1e-3)));
		REQUIRE(false == a.IsEqualTo(b, ScaleTolerance(REAL(1e-7))));
	}

	TEST_CASE("MatrixNM::Op+-", "[simple]") {
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		MatrixNM<Real, 2, 2> b({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		auto c = a + b;
		auto d = a - b;

		REQUIRE(REAL(2.0) == c(0, 0));
		REQUIRE(REAL(4.0) == c(0, 1));

		REQUIRE(REAL(0.0) == d(0, 0));
		REQUIRE(REAL(0.0) == d(0, 1));
	}

	TEST_CASE("MatrixNM::Op*", "[simple]") {
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		MatrixNM<Real, 2, 2> b({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		auto c = a * b;

		REQUIRE(REAL(7.0) == c(0, 0));
		REQUIRE(REAL(10.0) == c(0, 1));
		REQUIRE(REAL(15.0) == c(1, 0));
		REQUIRE(REAL(22.0) == c(1, 1));
	}

	TEST_CASE("MatrixNM::mul_double", "[simple]") {
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a({ REAL(1.0), REAL(100.0), REAL(50.0), REAL(100.0) });

		auto b = a * REAL(2.0);
		auto c = REAL(2.0) * a;

		REQUIRE(REAL(2.0) == b(0, 0));
		REQUIRE(REAL(2.0) == c(0, 0));

		REQUIRE(REAL(200.0) == b(0, 1));
		REQUIRE(REAL(200.0) == c(0, 1));
	}

	TEST_CASE("MatrixNM::div_double", "[simple]") {
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a({ REAL(4.0), REAL(400.0), REAL(1.0), 1 });

		auto b = a / REAL(2.0);

		REQUIRE(REAL(2.0) == b(0, 0));
		REQUIRE(REAL(200.0) == b(0, 1));
	}

	TEST_CASE("MatrixNM::mul_Vector", "[simple]") {
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a({ REAL(1.0), REAL(10.0),
																 REAL(5.0), REAL(2.0) });
		VectorN<Real, 2> b({ REAL(1.0), REAL(2.0) });

		auto c = a * b;
		auto d = b * a;

		REQUIRE(REAL(21.0) == c[0]);
		REQUIRE(REAL(9.0) == c[1]);

		REQUIRE(REAL(11.0) == d[0]);
		REQUIRE(REAL(14.0) == d[1]);
	}

	TEST_CASE("MatrixNM::Transpose", "[simple]")
	{
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> mat({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		MatrixNM<Real, 2, 2> matTransp({ REAL(1.0), REAL(3.0), REAL(2.0), REAL(4.0) });

		mat.Transpose();

		REQUIRE(mat.IsEqualTo(matTransp));
	}

	TEST_CASE("MatrixNM::transpose", "[simple]")
	{
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> mat({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		MatrixNM<Real, 2, 2> matTransp({ REAL(1.0), REAL(3.0), REAL(2.0), REAL(4.0) });

		auto trans = mat.transpose();

		REQUIRE(trans.IsEqualTo(matTransp));
	}

	TEST_CASE("MatrixNM::transpose_nonrectangular", "[simple]")
	{
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 3> mat({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0) });
		MatrixNM<Real, 3, 2> matTransp({ REAL(1.0), REAL(4.0), REAL(2.0), REAL(5.0), REAL(3.0), REAL(6.0) });

		auto trans = mat.transpose();

		REQUIRE(trans.IsEqualTo(matTransp));
	}

	TEST_CASE("MatrixNM::GetInverse", "[simple]")
	{
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> mat({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		auto b = mat.GetInverse();

		auto c = mat * b;

		MatrixNM<Real, 2, 2> d;
		d.MakeUnitMatrix();

		REQUIRE(c.IsEqualTo(d));
	}

	TEST_CASE("MatrixNM::exceptions", "[simple]")
	{
			TEST_PRECISION_INFO();
		MatrixNM<Real, 3, 4> mat_3;

		REQUIRE_THROWS_AS(mat_3.Invert(), MatrixDimensionError);
		REQUIRE_THROWS_AS(mat_3.GetInverse(), MatrixDimensionError);
		REQUIRE_THROWS_AS(mat_3.Transpose(), MatrixDimensionError);
	}

	///////////////////////          Copy and Move Semantics          //////////////////////
	TEST_CASE("MatrixNM::copy_constructor", "[MatrixNM][copy_move]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 3, 3> original({ REAL(1.0), REAL(2.0), REAL(3.0),
		                                REAL(4.0), REAL(5.0), REAL(6.0),
		                                REAL(7.0), REAL(8.0), REAL(9.0) });
		
		MatrixNM<Real, 3, 3> copy(original);
		
		REQUIRE(copy.IsEqualTo(original));
		
		// Modify copy, original should be unchanged
		copy(0, 0) = REAL(999.0);
		REQUIRE(original(0, 0) == REAL(1.0));
	}

	TEST_CASE("MatrixNM::copy_assignment", "[MatrixNM][copy_move]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		MatrixNM<Real, 2, 2> b;
		
		b = a;
		
		REQUIRE(b.IsEqualTo(a));
		
		// Modify b, a unchanged
		b(0, 0) = REAL(100.0);
		REQUIRE(a(0, 0) == REAL(1.0));
	}

	TEST_CASE("MatrixNM::scalar_assignment", "[MatrixNM][copy_move]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a;
		
		a = REAL(5.0);
		
		REQUIRE(a(0, 0) == REAL(5.0));
		REQUIRE(a(0, 1) == REAL(5.0));
		REQUIRE(a(1, 0) == REAL(5.0));
		REQUIRE(a(1, 1) == REAL(5.0));
	}

	TEST_CASE("MatrixNM::diagonal_constructor", "[MatrixNM][constructors]")
	{
		TEST_PRECISION_INFO();
		// Constructor from scalar creates diagonal matrix
		MatrixNM<Real, 3, 3> diag(REAL(2.0));
		
		REQUIRE(diag(0, 0) == REAL(2.0));
		REQUIRE(diag(1, 1) == REAL(2.0));
		REQUIRE(diag(2, 2) == REAL(2.0));
		REQUIRE(diag(0, 1) == REAL(0.0));
		REQUIRE(diag(1, 0) == REAL(0.0));
	}

	TEST_CASE("MatrixNM::nested_initializer_list_constructor", "[MatrixNM][constructors]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 3> m({
			{ REAL(1.0), REAL(2.0), REAL(3.0) },
			{ REAL(4.0), REAL(5.0), REAL(6.0) }
		});
		
		REQUIRE(m(0, 0) == REAL(1.0));
		REQUIRE(m(0, 2) == REAL(3.0));
		REQUIRE(m(1, 0) == REAL(4.0));
		REQUIRE(m(1, 2) == REAL(6.0));
	}

	TEST_CASE("MatrixNM::array_constructor", "[MatrixNM][constructors]")
	{
		TEST_PRECISION_INFO();
		Real arr[] = { REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) };
		MatrixNM<Real, 2, 2> m(arr, 4);
		
		REQUIRE(m(0, 0) == REAL(1.0));
		REQUIRE(m(0, 1) == REAL(2.0));
		REQUIRE(m(1, 0) == REAL(3.0));
		REQUIRE(m(1, 1) == REAL(4.0));
	}

	///////////////////////          Element Access                    //////////////////////
	TEST_CASE("MatrixNM::at_bounds_checked", "[MatrixNM][access]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 3> m({ REAL(1.0), REAL(2.0), REAL(3.0),
		                         REAL(4.0), REAL(5.0), REAL(6.0) });
		
		REQUIRE(m.at(0, 0) == REAL(1.0));
		REQUIRE(m.at(1, 2) == REAL(6.0));
		
		// Mutable access
		m.at(0, 1) = REAL(99.0);
		REQUIRE(m(0, 1) == REAL(99.0));
		
		// Out of bounds throws
		REQUIRE_THROWS_AS(m.at(-1, 0), MatrixAccessBoundsError);
		REQUIRE_THROWS_AS(m.at(0, 5), MatrixAccessBoundsError);
		REQUIRE_THROWS_AS(m.at(10, 0), MatrixAccessBoundsError);
	}

	///////////////////////          Triangular Extraction             //////////////////////
	TEST_CASE("MatrixNM::GetLower", "[MatrixNM][extraction]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 3, 3> m({ REAL(1.0), REAL(2.0), REAL(3.0),
		                         REAL(4.0), REAL(5.0), REAL(6.0),
		                         REAL(7.0), REAL(8.0), REAL(9.0) });
		
		auto lower = m.GetLower(true);  // With diagonal
		
		REQUIRE(lower(0, 0) == REAL(1.0));
		REQUIRE(lower(1, 0) == REAL(4.0));
		REQUIRE(lower(1, 1) == REAL(5.0));
		REQUIRE(lower(2, 0) == REAL(7.0));
		REQUIRE(lower(0, 1) == REAL(0.0));  // Upper part is zero
		REQUIRE(lower(0, 2) == REAL(0.0));
		
		auto lowerStrict = m.GetLower(false);  // Without diagonal
		REQUIRE(lowerStrict(0, 0) == REAL(0.0));
		REQUIRE(lowerStrict(1, 1) == REAL(0.0));
		REQUIRE(lowerStrict(1, 0) == REAL(4.0));
	}

	TEST_CASE("MatrixNM::GetUpper", "[MatrixNM][extraction]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 3, 3> m({ REAL(1.0), REAL(2.0), REAL(3.0),
		                         REAL(4.0), REAL(5.0), REAL(6.0),
		                         REAL(7.0), REAL(8.0), REAL(9.0) });
		
		auto upper = m.GetUpper(true);  // With diagonal
		
		REQUIRE(upper(0, 0) == REAL(1.0));
		REQUIRE(upper(0, 1) == REAL(2.0));
		REQUIRE(upper(0, 2) == REAL(3.0));
		REQUIRE(upper(1, 1) == REAL(5.0));
		REQUIRE(upper(1, 0) == REAL(0.0));  // Lower part is zero
		REQUIRE(upper(2, 0) == REAL(0.0));
		
		auto upperStrict = m.GetUpper(false);  // Without diagonal
		REQUIRE(upperStrict(0, 0) == REAL(0.0));
		REQUIRE(upperStrict(1, 1) == REAL(0.0));
		REQUIRE(upperStrict(0, 1) == REAL(2.0));
	}

	TEST_CASE("MatrixNM::GetLower_GetUpper_non_square_throws", "[MatrixNM][extraction]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 3> m;
		
		REQUIRE_THROWS_AS(m.GetLower(), MatrixDimensionError);
		REQUIRE_THROWS_AS(m.GetUpper(), MatrixDimensionError);
	}

	///////////////////////          Vector Extraction                 //////////////////////
	TEST_CASE("MatrixNM::VectorFromRow", "[MatrixNM][extraction]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 3> m({ REAL(1.0), REAL(2.0), REAL(3.0),
		                         REAL(4.0), REAL(5.0), REAL(6.0) });
		
		auto row0 = m.VectorFromRow(0);
		auto row1 = m.VectorFromRow(1);
		
		REQUIRE(row0[0] == REAL(1.0));
		REQUIRE(row0[1] == REAL(2.0));
		REQUIRE(row0[2] == REAL(3.0));
		REQUIRE(row1[0] == REAL(4.0));
		REQUIRE(row1[2] == REAL(6.0));
	}

	TEST_CASE("MatrixNM::VectorFromColumn", "[MatrixNM][extraction]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 3, 2> m({ REAL(1.0), REAL(2.0),
		                         REAL(3.0), REAL(4.0),
		                         REAL(5.0), REAL(6.0) });
		
		auto col0 = m.VectorFromColumn(0);
		auto col1 = m.VectorFromColumn(1);
		
		REQUIRE(col0[0] == REAL(1.0));
		REQUIRE(col0[1] == REAL(3.0));
		REQUIRE(col0[2] == REAL(5.0));
		REQUIRE(col1[0] == REAL(2.0));
		REQUIRE(col1[2] == REAL(6.0));
	}

	TEST_CASE("MatrixNM::VectorFromDiagonal", "[MatrixNM][extraction]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 3, 3> m({ REAL(1.0), REAL(2.0), REAL(3.0),
		                         REAL(4.0), REAL(5.0), REAL(6.0),
		                         REAL(7.0), REAL(8.0), REAL(9.0) });
		
		auto diag = m.VectorFromDiagonal();
		
		REQUIRE(diag[0] == REAL(1.0));
		REQUIRE(diag[1] == REAL(5.0));
		REQUIRE(diag[2] == REAL(9.0));
	}

	///////////////////////          Arithmetic Operators              //////////////////////
	TEST_CASE("MatrixNM::unary_negation", "[MatrixNM][arithmetic]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a({ REAL(1.0), REAL(-2.0), REAL(3.0), REAL(-4.0) });
		
		auto neg = -a;
		
		REQUIRE(neg(0, 0) == -REAL(1.0));
		REQUIRE(neg(0, 1) == REAL(2.0));
		REQUIRE(neg(1, 0) == -REAL(3.0));
		REQUIRE(neg(1, 1) == REAL(4.0));
	}

	TEST_CASE("MatrixNM::compound_multiply_scalar", "[MatrixNM][arithmetic]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		
		a *= REAL(2.0);
		
		REQUIRE(a(0, 0) == REAL(2.0));
		REQUIRE(a(0, 1) == REAL(4.0));
		REQUIRE(a(1, 0) == REAL(6.0));
		REQUIRE(a(1, 1) == REAL(8.0));
	}

	TEST_CASE("MatrixNM::rectangular_multiply", "[MatrixNM][arithmetic]")
	{
		TEST_PRECISION_INFO();
		// (2x3) * (3x4) = (2x4)
		MatrixNM<Real, 2, 3> a({ REAL(1.0), REAL(2.0), REAL(3.0),
		                         REAL(4.0), REAL(5.0), REAL(6.0) });
		MatrixNM<Real, 3, 4> b({ REAL(1.0), REAL(0.0), REAL(0.0), REAL(0.0),
		                         REAL(0.0), REAL(1.0), REAL(0.0), REAL(0.0),
		                         REAL(0.0), REAL(0.0), REAL(1.0), REAL(0.0) });
		
		auto c = a * b;
		
		REQUIRE(c.rows() == 2);
		REQUIRE(c.cols() == 4);
		REQUIRE(c(0, 0) == REAL(1.0));
		REQUIRE(c(0, 1) == REAL(2.0));
		REQUIRE(c(0, 2) == REAL(3.0));
		REQUIRE(c(0, 3) == REAL(0.0));
	}

	///////////////////////          Equality Operators                //////////////////////
	TEST_CASE("MatrixNM::operator_equals", "[MatrixNM][equality]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		MatrixNM<Real, 2, 2> b({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		MatrixNM<Real, 2, 2> c({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(5.0) });
		
		REQUIRE(a == b);
		REQUIRE(!(a == c));
	}

	TEST_CASE("MatrixNM::operator_not_equals", "[MatrixNM][equality]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		MatrixNM<Real, 2, 2> b({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(5.0) });
		
		REQUIRE(a != b);
	}

	TEST_CASE("MatrixNM::AreEqual_static", "[MatrixNM][equality]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		MatrixNM<Real, 2, 2> b({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0001) });
		
		REQUIRE(MatrixNM<Real, 2, 2>::AreEqual(a, b, REAL(0.001)) == true);
		REQUIRE(MatrixNM<Real, 2, 2>::AreEqual(a, b, REAL(0.00001)) == false);
	}

	///////////////////////          Trace                             //////////////////////
	TEST_CASE("MatrixNM::Trace", "[MatrixNM][trace]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 3, 3> m({ REAL(1.0), REAL(2.0), REAL(3.0),
		                         REAL(4.0), REAL(5.0), REAL(6.0),
		                         REAL(7.0), REAL(8.0), REAL(9.0) });
		
		REQUIRE_THAT(m.Trace(), RealWithinRel(REAL(15.0)));  // 1 + 5 + 9
	}

	TEST_CASE("MatrixNM::Trace_non_square_throws", "[MatrixNM][trace]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 3> m;
		
		REQUIRE_THROWS_AS(m.Trace(), MatrixDimensionError);
	}

	///////////////////////          In-place Invert                   //////////////////////
	TEST_CASE("MatrixNM::Invert_in_place", "[MatrixNM][invert]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> original({ REAL(4.0), REAL(7.0), REAL(2.0), REAL(6.0) });
		MatrixNM<Real, 2, 2> m(original);
		
		m.Invert();
		
		// m * original = I
		auto product = m * original;
		auto identity = MatrixNM<Real, 2, 2>::Identity();
		REQUIRE(product.IsEqualTo(identity, TOL(1e-10, 1e-5)));
	}

	///////////////////////          Type Aliases                      //////////////////////
	TEST_CASE("MatrixNM::type_aliases", "[MatrixNM][aliases]")
	{
		TEST_PRECISION_INFO();
		
		// Float aliases
		Mat22F mf22;
		mf22(0, 0) = 1.0f;
		REQUIRE(mf22(0, 0) == 1.0f);
		
		Mat33F mf33;
		REQUIRE(mf33.rows() == 3);
		REQUIRE(mf33.cols() == 3);
		
		Mat44F mf44;
		REQUIRE(mf44.rows() == 4);
		
		// Double aliases
		Mat22D md22({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		REQUIRE(md22(0, 0) == REAL(1.0));
		
		Mat33D md33;
		REQUIRE(md33.rows() == 3);
		
		Mat44D md44;
		REQUIRE(md44.rows() == 4);
	}

	TEST_CASE("MatrixNM::complex_type", "[MatrixNM][complex]")
	{
		TEST_PRECISION_INFO();
		
		Mat22C mc({ Complex(1.0, 2.0), Complex(3.0, 4.0),
		            Complex(5.0, 6.0), Complex(7.0, 8.0) });
		
		REQUIRE(mc(0, 0).real() == 1.0);
		REQUIRE(mc(0, 0).imag() == 2.0);
		REQUIRE(mc(1, 1).real() == 7.0);
	}

	///////////////////////          Output Methods                    //////////////////////
	TEST_CASE("MatrixNM::to_string", "[MatrixNM][output]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> m({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		
		std::string str = m.to_string(5, 2);
		
		REQUIRE(!str.empty());
		REQUIRE(str.find("Rows: 2") != std::string::npos);
		REQUIRE(str.find("Cols: 2") != std::string::npos);
	}

	TEST_CASE("MatrixNM::Print_stream", "[MatrixNM][output]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> m({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		
		std::stringstream ss;
		ss << m;
		
		std::string output = ss.str();
		REQUIRE(!output.empty());
	}

	TEST_CASE("MatrixNM::Print_with_format", "[MatrixNM][output]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> m({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		
		std::stringstream ss;
		m.Print(ss, MatrixPrintFormat::Scientific());
		
		std::string output = ss.str();
		REQUIRE(output.find("e") != std::string::npos);  // Scientific notation
	}

	///////////////////////          Edge Cases                        //////////////////////
	TEST_CASE("MatrixNM::identity_multiply", "[MatrixNM][edge_cases]")
	{
		TEST_PRECISION_INFO();
		auto I = MatrixNM<Real, 3, 3>::Identity();
		MatrixNM<Real, 3, 3> A({ REAL(1.0), REAL(2.0), REAL(3.0),
		                         REAL(4.0), REAL(5.0), REAL(6.0),
		                         REAL(7.0), REAL(8.0), REAL(10.0) });  // Non-singular
		
		// A * I = A
		auto AI = A * I;
		REQUIRE(A.IsEqualTo(AI, TOL(1e-10, 1e-5)));
		
		// I * A = A
		auto IA = I * A;
		REQUIRE(A.IsEqualTo(IA, TOL(1e-10, 1e-5)));
	}

	TEST_CASE("MatrixNM::1x1_matrix", "[MatrixNM][edge_cases]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 1, 1> m({ REAL(5.0) });
		
		REQUIRE(m.rows() == 1);
		REQUIRE(m.cols() == 1);
		REQUIRE(m(0, 0) == REAL(5.0));
		
		// Inverse of 1x1 [a] is [1/a]
		auto inv = m.GetInverse();
		REQUIRE_THAT(inv(0, 0), RealWithinRel(REAL(0.2)));
		
		// Trace of 1x1
		REQUIRE(m.Trace() == REAL(5.0));
	}

	TEST_CASE("MatrixNM::MakeUnitMatrix_non_square_throws", "[MatrixNM][edge_cases]")
	{
		TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 3> m;
		
		REQUIRE_THROWS_AS(m.MakeUnitMatrix(), MatrixDimensionError);
	}
}

