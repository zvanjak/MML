#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/BaseUtils.h"
#endif

using namespace MML;

using namespace Catch::Matchers;

namespace MML::Tests::Base::BaseUtilsTests
{
	TEST_CASE("BaseUtils_Abs")
	{
		Complex a(1, -3);

		double b = Abs(a);

		REQUIRE(b == sqrt(10));
	}

	///////////////////                     Complex helpers                   ///////////////////
	TEST_CASE("BaseUtils_Complex_AreEqual", "[simple]")
	{
		Complex a(1, 2);
		Complex b(1, 2);

		REQUIRE(Utils::AreEqual(a, b));

		Complex c(1, 2.0001);
		REQUIRE(!Utils::AreEqual(a, c));
		REQUIRE(Utils::AreEqual(a, c, 1e-03));
		REQUIRE(!Utils::AreEqual(a, c, 1e-04));

		Complex d(1, 2.0000000001);
		REQUIRE(Utils::AreEqual(a, d, 1e-08));
		REQUIRE(Utils::AreEqual(a, d, 1e-09));
		REQUIRE(!Utils::AreEqual(a, d, 1e-10));
		REQUIRE(!Utils::AreEqual(a, d, 1e-11));
	}
	TEST_CASE("BaseUtils_Complex_AreEqualAbs", "[simple]")
	{
		Complex a(1, 2);
		Complex b(1, 2);

		REQUIRE(Utils::AreEqualAbs(a, b, 1e-08));

		Complex c(1, 2.0001);
		REQUIRE(Utils::AreEqualAbs(a, c, 1e-03));
		REQUIRE(!Utils::AreEqualAbs(a, c, 1e-04));

		Complex d(1, 2.0000000001);
		REQUIRE(Utils::AreEqualAbs(a, d, 1e-08));
		REQUIRE(Utils::AreEqualAbs(a, d, 1e-09));
		REQUIRE(!Utils::AreEqualAbs(a, d, 1e-10));
		REQUIRE(!Utils::AreEqualAbs(a, d, 1e-11));
	}
	TEST_CASE("BaseUtils_Vector<Complex>_AreEqual", "[simple]")
	{
		Vector<Complex> a({ Complex(1, 2), Complex(3, 4) });
		Vector<Complex> b({ Complex(1, 2), Complex(3, 4) });

		REQUIRE(Utils::AreEqual(a, b));

		Vector<Complex> c({ Complex(1, 2), Complex(3, 4.0001) });
		REQUIRE(!Utils::AreEqual(a, c));
		REQUIRE(!Utils::AreEqual(a, c, 1e-05));

		Vector<Complex> d({ Complex(1, 2), Complex(3, 4.0000000001) });
		REQUIRE(Utils::AreEqual(a, d, 1e-08));
		REQUIRE(Utils::AreEqual(a, d, 1e-09));
		REQUIRE(!Utils::AreEqual(a, d, 1e-10));
		REQUIRE(!Utils::AreEqual(a, d, 1e-11));
	}
	TEST_CASE("BaseUtils_Vector<Complex>_AreEqualAbs", "[simple]")
	{
		Vector<Complex> a({ Complex(1, 2), Complex(3, 4) });
		Vector<Complex> b({ Complex(1, 2), Complex(3, 4) });

		REQUIRE(Utils::AreEqualAbs(a, b));
		REQUIRE(Utils::AreEqualAbs(a, b, 1e-08));

		Vector<Complex> c({ Complex(1, 2), Complex(3, 4.0001) });
		REQUIRE(Utils::AreEqualAbs(a, c, 1e-03));
		REQUIRE(!Utils::AreEqualAbs(a, c, 1e-05));

		Vector<Complex> d({ Complex(1, 2), Complex(3, 4.0000000001) });
		REQUIRE(Utils::AreEqualAbs(a, d, 1e-08));
		REQUIRE(Utils::AreEqualAbs(a, d, 1e-09));
		REQUIRE(!Utils::AreEqualAbs(a, d, 1e-10));
		REQUIRE(!Utils::AreEqualAbs(a, d, 1e-11));
	}

	//////////////////                     Vector helpers                    ///////////////////
	TEST_CASE("BaseUtils_VectorProjectionParallelTo", "[simple]")
	{
		Vector<Real> a({ 1, 2 });
		Vector<Real> b({ 3, 4 });

		auto c = Utils::VectorProjectionParallelTo(a, b);

		REQUIRE(2 == c.size());
    REQUIRE_THAT(6.6, WithinAbs(c[0], 1e-15));
    REQUIRE_THAT(8.8, WithinAbs(c[1], 1e-15));
	}
	TEST_CASE("BaseUtils_VectorProjectionPerpendicularTo", "[simple]")
	{
		Vector<Real> a({ 1, 2 });
		Vector<Real> b({ 3, 4 });

		auto c = Utils::VectorProjectionPerpendicularTo(a, b);

		REQUIRE(2 == c.size());
    REQUIRE_THAT(-5.6, WithinAbs(c[0], 1e-15));
    REQUIRE_THAT(-6.8, WithinAbs(c[1], 1e-15));
	}
	TEST_CASE("BaseUtils_VectorScalarProductReal", "[simple]")
	{
		Vector<Real> a({ 1, 2 });
		Vector<Real> b({ 3, 4 });

		auto c = Utils::ScalarProduct(a, b);

		REQUIRE(11 == c);
	}
	TEST_CASE("BaseUtils_VectorScalarProductComplex", "[simple]")
	{
		Vector<Complex> a({ Complex(1, 2), Complex(3, 4) });
		Vector<Complex> b({ Complex(5, 6), Complex(7, 8) });

		auto c = Utils::ScalarProduct(a, b);

		REQUIRE(Complex(70.0, 8.0) == c);
	}
	TEST_CASE("BaseUtils_OuterProductReal", "[simple]")
	{
		Vector<Real> a({ 1, 2 });
		Vector<Real> b({ 3, 4 });
		Vector<Real> c({ 3, 4, 5 });

		auto a_x_b = Utils::OuterProduct(a, b);

		REQUIRE(2 == a_x_b.RowNum());
		REQUIRE(2 == a_x_b.ColNum());

		REQUIRE(3 == a_x_b(0, 0));
		REQUIRE(4 == a_x_b(0, 1));
		REQUIRE(6 == a_x_b(1, 0));
		REQUIRE(8 == a_x_b(1, 1));

		auto a_x_c = Utils::OuterProduct(a, c);

		REQUIRE(2 == a_x_c.RowNum());
		REQUIRE(3 == a_x_c.ColNum());

		REQUIRE(3 == a_x_c(0, 0));
		REQUIRE(4 == a_x_c(0, 1));
		REQUIRE(5 == a_x_c(0, 2));
		REQUIRE(6 == a_x_c(1, 0));
		REQUIRE(8 == a_x_c(1, 1));
		REQUIRE(10 == a_x_c(1, 2));
	}
	TEST_CASE("BaseUtils_OuterProductComplex", "[simple]")
	{
		Vector<Complex> a({ Complex(1, 2), Complex(3, 4) });
		Vector<Complex> b({ Complex(5, 6), Complex(7, 8) });

		auto c = Utils::OuterProduct(a, b);

		REQUIRE(2 == c.RowNum());
		REQUIRE(2 == c.ColNum());

		REQUIRE(Complex(-7, 16) == c(0, 0));
		REQUIRE(Complex(-9, 22) == c(0, 1));
		REQUIRE(Complex(-9, 38) == c(1, 0));
		REQUIRE(Complex(-11, 52) == c(1, 1));
	}

  ///////////////////       Vector<Complex> - Vector<Real> operations       ///////////////////
	TEST_CASE("BaseUtils_Utils_AddVec", "[simple]") {
		Vector<Complex> a({ Complex(1,0), Complex(0,1) });
		Vector<Real>    b({ 1,2 });

		auto c = Utils::AddVec(a, b);

		REQUIRE(2 == c.size());
		REQUIRE(Complex(2, 0) == c[0]);
		REQUIRE(Complex(2, 1) == c[1]);

		auto d = Utils::AddVec(b, a);
		REQUIRE(2 == d.size());
		REQUIRE(Complex(2, 0) == d[0]);
		REQUIRE(Complex(2, 1) == d[1]);
	}
	TEST_CASE("BaseUtils_Utils_SubVec", "[simple]") {
		Vector<Complex> a({ Complex(1,0), Complex(0,1) });
		Vector<Real>    b({ 1,2 });

		auto c = Utils::SubVec(a, b);

		REQUIRE(2 == c.size());
		REQUIRE(Complex(0, 0) == c[0]);
		REQUIRE(Complex(-2, 1) == c[1]);

		auto d = Utils::SubVec(b, a);
		REQUIRE(2 == d.size());
		REQUIRE(Complex(0, 0) == d[0]);
		REQUIRE(Complex(2, -1) == d[1]);
	}

	//////////////////             Creating Matrix from Vector              ///////////////////	
	TEST_CASE("MatrixUtils_RowMatrixFromVector", "[simple]") {
		Vector<Real> a{ 1.0, 2.0, 3.0 };

		auto b = MatrixUtils::RowMatrixFromVector<Real>(a);

		REQUIRE(1 == b.RowNum());
		REQUIRE(3 == b.ColNum());

		REQUIRE(1.0 == b[0][0]);
		REQUIRE(2.0 == b[0][1]);
		REQUIRE(3.0 == b[0][2]);
	}
	TEST_CASE("MatrixUtils_ColumnMatrixFromVector", "[simple]") {
		Vector<Real> a{ 1.0, 2.0, 3.0 };

		auto b = MatrixUtils::ColumnMatrixFromVector<Real>(a);

		REQUIRE(3 == b.RowNum());
		REQUIRE(1 == b.ColNum());

		REQUIRE(1.0 == b[0][0]);
		REQUIRE(2.0 == b[1][0]);
		REQUIRE(3.0 == b[2][0]);
	}

	///////////////////                   Matrix helpers                     ///////////////////
	TEST_CASE("MatrixUtils::Commutator", "[MatrixUtils]") 
	{
	}
	TEST_CASE("MatrixUtils::AntiCommutator", "[MatrixUtils]") 
	{
	}
	TEST_CASE("MatrixUtils::MatrixDecomposeToSymAntisym", "[MatrixUtils]") 
	{
	}
	///////////////////                  Matrix functions                    ///////////////////
	TEST_CASE("MatrixUtils::Exp", "[MatrixUtils]") 
	{
	}

	///////////////////                Real matrix helpers                   ///////////////////
	TEST_CASE("MatrixUtils::IsOrthogonal", "[MatrixUtils]") 
	{
	}

	///////////////////               Complex matrix helpers                 ///////////////////
	TEST_CASE("MatrixUtils::GetRealPart", "[MatrixUtils]") 
	{
	}
	TEST_CASE("MatrixUtils::GetImagPart", "[MatrixUtils]") 
	{
	}
	TEST_CASE("MatrixUtils::GetConjugateTranspose", "[MatrixUtils]") 
	{
	}
	TEST_CASE("MatrixUtils::CmplxMatFromRealMat", "[MatrixUtils]") 
	{
	}
	TEST_CASE("MatrixUtils::IsComplexMatReal", "[MatrixUtils]") 
	{
	}
	TEST_CASE("MatrixUtils::IsHermitian", "[MatrixUtils]") 
	{
	}
	TEST_CASE("MatrixUtils::IsUnitary", "[MatrixUtils]") 
	{
	}
	
	///////////////////       Matrix<Complex> - Matrix<Real>  operations     ///////////////////
	TEST_CASE("MatrixUtils::AddMat", "[MatrixUtils]") 
	{
	}
	TEST_CASE("MatrixUtils::SubMat", "[MatrixUtils]") 
	{
	}
	TEST_CASE("MatrixUtils::MulMat1", "[MatrixUtils]") 
	{
	}
	TEST_CASE("MatrixUtils::MulMat2", "[MatrixUtils]") 
	{
	}
	TEST_CASE("MatrixUtils::MulMatVec", "[MatrixUtils]") 
	{
	}
	TEST_CASE("MatrixUtils::MulVecMat", "[MatrixUtils]") 
	{
	}
}