#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/MatrixNM.h"
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

		REQUIRE(2 == a.RowNum());
		REQUIRE(2 == a.ColNum());

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

	TEST_CASE("MatrixNM::GetUnitMatrix", "[simple]") {
			TEST_PRECISION_INFO();
		auto a = MatrixNM<Real, 2, 2>::GetUnitMatrix();

		REQUIRE(REAL(1.0) == a[0][0]);
		REQUIRE(REAL(0.0) == a[0][1]);
		REQUIRE(REAL(0.0) == a[1][0]);
		REQUIRE(REAL(1.0) == a[1][1]);
	}

	TEST_CASE("MatrixNM::IsEqual", "[simple]") {
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		MatrixNM<Real, 2, 2> b({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		REQUIRE(true == a.IsEqual(b));
	}
	TEST_CASE("MatrixNM::IsEqual2", "[simple]") {
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		MatrixNM<Real, 2, 2> b({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(5.0) });

		REQUIRE(false == a.IsEqual(b));
	}
	TEST_CASE("MatrixNM::IsEqual3", "[simple]") {
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> a({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		MatrixNM<Real, 2, 2> b({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0001) });

		// Difference is 1e-4, use 1e-3 (yes) and 1e-7 (no)
		// 1e-3 scales to 1e-1 (float), 1e-7 scales to 1e-5 (float)
		REQUIRE(true == a.IsEqual(b, ScaleTolerance(REAL(1e-3))));
		REQUIRE(false == a.IsEqual(b, ScaleTolerance(REAL(1e-7))));
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

		REQUIRE(mat.IsEqual(matTransp));
	}

	TEST_CASE("MatrixNM::GetTranspose", "[simple]")
	{
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> mat({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });
		MatrixNM<Real, 2, 2> matTransp({ REAL(1.0), REAL(3.0), REAL(2.0), REAL(4.0) });

		auto trans = mat.GetTranspose();

		REQUIRE(trans.IsEqual(matTransp));
	}

	TEST_CASE("MatrixNM::GetTranspose_nonrectangular", "[simple]")
	{
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 3> mat({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0), REAL(6.0) });
		MatrixNM<Real, 3, 2> matTransp({ REAL(1.0), REAL(4.0), REAL(2.0), REAL(5.0), REAL(3.0), REAL(6.0) });

		auto trans = mat.GetTranspose();

		REQUIRE(trans.IsEqual(matTransp));
	}

	TEST_CASE("MatrixNM::GetInverse", "[simple]")
	{
			TEST_PRECISION_INFO();
		MatrixNM<Real, 2, 2> mat({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0) });

		auto b = mat.GetInverse();

		auto c = mat * b;

		MatrixNM<Real, 2, 2> d;
		d.MakeUnitMatrix();

		REQUIRE(c.IsEqual(d));
	}

	TEST_CASE("MatrixNM::exceptions", "[simple]")
	{
			TEST_PRECISION_INFO();
		MatrixNM<Real, 3, 4> mat_3;

		REQUIRE_THROWS_AS(mat_3.Invert(), MatrixDimensionError);
		REQUIRE_THROWS_AS(mat_3.GetInverse(), MatrixDimensionError);
		REQUIRE_THROWS_AS(mat_3.Transpose(), MatrixDimensionError);
	}
}