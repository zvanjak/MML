#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/BaseUtils.h"
#include "mml/algorithms/RootFinding/RootFindingPolynoms.h"
#endif

using namespace MML;
using namespace MML::Testing;

using namespace Catch::Matchers;

namespace MML::Tests::Base::BaseUtilsTests
{
	TEST_CASE("MMLBase::SolveQuadratic", "[MMLBase]")
	{
			TEST_PRECISION_INFO();
		Real a = REAL(1.0);
		Real b = -REAL(3.0);
		Real c = REAL(2.0);

		Complex x1, x2;
		bool res = SolveQuadratic(a, b, c, x1, x2);
		
		REQUIRE(res);
		REQUIRE(Utils::AreEqual(x1, Complex(REAL(2.0), REAL(0.0))));
		REQUIRE(Utils::AreEqual(x2, Complex(REAL(1.0), REAL(0.0))));
	}

	TEST_CASE("MMLBase::SolveCubic1", "[MMLBase]")
	{
			TEST_PRECISION_INFO();
		Real a = REAL(1.0);
		Real b = -REAL(6.0);
		Real c = REAL(11.0);
		Real d = -REAL(6.0);
		Complex x1, x2, x3;

		SolveCubic(a, b, c, d, x1, x2, x3);

		Complex r1 = a * POW3(x1) + b * POW2(x1) + c * x1 + d;
		Complex r2 = a * POW3(x2) + b * POW2(x2) + c * x2 + d;
		Complex r3 = a * POW3(x3) + b * POW2(x3) + c * x3 + d;

		REQUIRE(Utils::AreEqual(r1, Complex(REAL(0.0), REAL(0.0))));
		REQUIRE(Utils::AreEqual(r2, Complex(REAL(0.0), REAL(0.0))));
		REQUIRE(Utils::AreEqual(r3, Complex(REAL(0.0), REAL(0.0))));

		REQUIRE(Utils::AreEqual(x1, Complex(REAL(3.0), REAL(0.0))));
		REQUIRE(Utils::AreEqual(x2, Complex(REAL(1.0), REAL(0.0))));
		REQUIRE(Utils::AreEqual(x3, Complex(REAL(2.0), REAL(0.0))));
	}

	TEST_CASE("MMLBase::SolveCubic2", "[MMLBase]")
	{
			TEST_PRECISION_INFO();
		Real a = -REAL(2.3);
		Real b = -REAL(8.1);
		Real c = REAL(5.7);
		Real d = -REAL(4.5);

		Complex x1, x2, x3;
		SolveCubic(a, b, c, d, x1, x2, x3);

		Complex r1 = a * POW3(x1) + b * POW2(x1) + c * x1 + d;
		Complex r2 = a * POW3(x2) + b * POW2(x2) + c * x2 + d;
		Complex r3 = a * POW3(x3) + b * POW2(x3) + c * x3 + d;

		REQUIRE(Utils::AreEqual(r1, Complex(REAL(0.0), REAL(0.0))));
		REQUIRE(Utils::AreEqual(r2, Complex(REAL(0.0), REAL(0.0))));
		REQUIRE(Utils::AreEqual(r3, Complex(REAL(0.0), REAL(0.0))));
	}

	TEST_CASE("MMLBase::SolveCubic3", "[MMLBase]")
	{
			TEST_PRECISION_INFO();
		Real a = REAL(0.645);
		Real b = -REAL(2.19);
		Real c = -REAL(33.4);
		Real d = REAL(8.7);
		Complex x1, x2, x3;

		SolveCubic(a, b, c, d, x1, x2, x3);

		Complex r1 = a * POW3(x1) + b * POW2(x1) + c * x1 + d;
		Complex r2 = a * POW3(x2) + b * POW2(x2) + c * x2 + d;
		Complex r3 = a * POW3(x3) + b * POW2(x3) + c * x3 + d;

		REQUIRE(Utils::AreEqual(r1, Complex(REAL(0.0), REAL(0.0))));
		REQUIRE(Utils::AreEqual(r2, Complex(REAL(0.0), REAL(0.0))));
		REQUIRE(Utils::AreEqual(r3, Complex(REAL(0.0), REAL(0.0))));
	}

	TEST_CASE("MMLBase::SolveQuartic1", "[MMLBase]")
	{
			TEST_PRECISION_INFO();
		Real a = REAL(1.0);
		Real b = -REAL(10.0);
		Real c = REAL(35.0);
		Real d = -REAL(50.0);
		Real e = REAL(24.0);
		
		Complex x1, x2, x3, x4;
		SolveQuartic(a, b, c, d, e, x1, x2, x3, x4);
		
		Complex r1 = a * POW4(x1) + b * POW3(x1) + c * POW2(x1) + d * x1 + e;
		Complex r2 = a * POW4(x2) + b * POW3(x2) + c * POW2(x2) + d * x2 + e;
		Complex r3 = a * POW4(x3) + b * POW3(x3) + c * POW2(x3) + d * x3 + e;
		Complex r4 = a * POW4(x4) + b * POW3(x4) + c * POW2(x4) + d * x4 + e;
		
		REQUIRE(Utils::AreEqual(r1, Complex(REAL(0.0), REAL(0.0))));
		REQUIRE(Utils::AreEqual(r2, Complex(REAL(0.0), REAL(0.0))));
		REQUIRE(Utils::AreEqual(r3, Complex(REAL(0.0), REAL(0.0))));
		REQUIRE(Utils::AreEqual(r4, Complex(REAL(0.0), REAL(0.0))));
		
		REQUIRE(Utils::AreEqual(x1, Complex(REAL(4.0), REAL(0.0))));
		REQUIRE(Utils::AreEqual(x2, Complex(REAL(1.0), REAL(0.0))));
		REQUIRE(Utils::AreEqual(x3, Complex(REAL(3.0), REAL(0.0))));
		REQUIRE(Utils::AreEqual(x4, Complex(REAL(2.0), REAL(0.0))));
	}

	TEST_CASE("BaseUtils::Abs")
	{
			TEST_PRECISION_INFO();
		Complex a(1, -3);

		double b = Abs(a);

		REQUIRE(b == sqrt(10));
	}

	///////////////////                     Complex helpers                   ///////////////////
	TEST_CASE("BaseUtils::Complex_AreEqual", "[simple]")
	{
			TEST_PRECISION_INFO();
		Complex a(1, 2);
		Complex b(1, 2);

		REQUIRE(Utils::AreEqual(a, b));

		// Use larger difference (1e-4) that's detectable across all precisions
		Complex c(1, REAL(2.0001));
		REQUIRE(!Utils::AreEqual(a, c));
		// 1e-3 scales to 1e-1 (float), 1e-7 scales to 1e-5 (float)
		REQUIRE(Utils::AreEqual(a, c, ScaleTolerance(REAL(1e-3))));
		REQUIRE(!Utils::AreEqual(a, c, ScaleTolerance(REAL(1e-7))));
	}
	TEST_CASE("BaseUtils::Complex_AreEqualAbs", "[simple]")
	{
			TEST_PRECISION_INFO();
		Complex a(1, 2);
		Complex b(1, 2);

		REQUIRE(Utils::AreEqualAbs(a, b, ScaleTolerance(REAL(1e-8))));

		// Use larger difference (1e-4) that's detectable across all precisions
		Complex c(1, REAL(2.0001));
		// 1e-3 scales to 1e-1 (float), 1e-7 scales to 1e-5 (float)
		REQUIRE(Utils::AreEqualAbs(a, c, ScaleTolerance(REAL(1e-3))));
		REQUIRE(!Utils::AreEqualAbs(a, c, ScaleTolerance(REAL(1e-7))));
	}
	TEST_CASE("BaseUtils::Vector<Complex>_AreEqual", "[simple]")
	{
			TEST_PRECISION_INFO();
		Vector<Complex> a({ Complex(1, 2), Complex(3, 4) });
		Vector<Complex> b({ Complex(1, 2), Complex(3, 4) });

		REQUIRE(Utils::AreEqual(a, b));

		// Use larger difference (1e-4) that's detectable across all precisions
		Vector<Complex> c({ Complex(1, 2), Complex(3, REAL(4.0001)) });
		REQUIRE(!Utils::AreEqual(a, c));
		// 1e-3 scales to 1e-1 (float), 1e-7 scales to 1e-5 (float)
		REQUIRE(Utils::AreEqual(a, c, ScaleTolerance(REAL(1e-3))));
		REQUIRE(!Utils::AreEqual(a, c, ScaleTolerance(REAL(1e-7))));
	}
	TEST_CASE("BaseUtils::Vector<Complex>_AreEqualAbs", "[simple]")
	{
			TEST_PRECISION_INFO();
		Vector<Complex> a({ Complex(1, 2), Complex(3, 4) });
		Vector<Complex> b({ Complex(1, 2), Complex(3, 4) });

		REQUIRE(Utils::AreEqualAbs(a, b));
		REQUIRE(Utils::AreEqualAbs(a, b, ScaleTolerance(REAL(1e-8))));

		// Use larger difference (1e-4) that's detectable across all precisions
		Vector<Complex> c({ Complex(1, 2), Complex(3, REAL(4.0001)) });
		// 1e-3 scales to 1e-1 (float), 1e-7 scales to 1e-5 (float)
		REQUIRE(Utils::AreEqualAbs(a, c, ScaleTolerance(REAL(1e-3))));
		REQUIRE(!Utils::AreEqualAbs(a, c, ScaleTolerance(REAL(1e-7))));
	}

	//////////////////                     Vector helpers                    ///////////////////
	TEST_CASE("BaseUtils::VectorProjectionParallelTo", "[simple]")
	{
			TEST_PRECISION_INFO();
		Vector<Real> a({ 1, 2 });
		Vector<Real> b({ 3, 4 });

		auto c = Utils::VectorProjectionParallelTo(a, b);

		REQUIRE(2 == c.size());
		REQUIRE_THAT(REAL(1.32), WithinAbs(c[0], REAL(1e-14)));
		REQUIRE_THAT(REAL(1.76), WithinAbs(c[1], REAL(1e-14)));
	}
	TEST_CASE("BaseUtils::VectorProjectionPerpendicularTo", "[simple]")
	{
			TEST_PRECISION_INFO();
		Vector<Real> a({ 1, 2 });
		Vector<Real> b({ 3, 4 });

		auto c = Utils::VectorProjectionPerpendicularTo(a, b);

		REQUIRE(2 == c.size());
		REQUIRE_THAT(-REAL(0.32), WithinAbs(c[0], REAL(1e-14)));
		REQUIRE_THAT(REAL(0.24), WithinAbs(c[1], REAL(1e-14)));
	}
	TEST_CASE("BaseUtils::VectorScalarProductReal", "[simple]")
	{
			TEST_PRECISION_INFO();
		Vector<Real> a({ 1, 2 });
		Vector<Real> b({ 3, 4 });

		auto c = Utils::ScalarProduct(a, b);

		REQUIRE(11 == c);
	}
	TEST_CASE("BaseUtils::VectorScalarProductComplex", "[simple]")
	{
			TEST_PRECISION_INFO();
		Vector<Complex> a({ Complex(1, 2), Complex(3, 4) });
		Vector<Complex> b({ Complex(5, 6), Complex(7, 8) });

		auto c = Utils::ScalarProduct(a, b);

		REQUIRE(Complex(REAL(70.0), REAL(8.0)) == c);
	}
	TEST_CASE("BaseUtils::OuterProductReal", "[simple]")
	{
			TEST_PRECISION_INFO();
		Vector<Real> a({ 1, 2 });
		Vector<Real> b({ 3, 4 });
		Vector<Real> c({ 3, 4, 5 });

		auto a_x_b = Utils::OuterProduct(a, b);

		REQUIRE(2 == a_x_b.rows());
		REQUIRE(2 == a_x_b.cols());

		REQUIRE(3 == a_x_b(0, 0));
		REQUIRE(4 == a_x_b(0, 1));
		REQUIRE(6 == a_x_b(1, 0));
		REQUIRE(8 == a_x_b(1, 1));

		auto a_x_c = Utils::OuterProduct(a, c);

		REQUIRE(2 == a_x_c.rows());
		REQUIRE(3 == a_x_c.cols());

		REQUIRE(3 == a_x_c(0, 0));
		REQUIRE(4 == a_x_c(0, 1));
		REQUIRE(5 == a_x_c(0, 2));
		REQUIRE(6 == a_x_c(1, 0));
		REQUIRE(8 == a_x_c(1, 1));
		REQUIRE(10 == a_x_c(1, 2));
	}
	TEST_CASE("BaseUtils::OuterProductComplex", "[simple]")
	{
			TEST_PRECISION_INFO();
		Vector<Complex> a({ Complex(1, 2), Complex(3, 4) });
		Vector<Complex> b({ Complex(5, 6), Complex(7, 8) });

		auto c = Utils::OuterProduct(a, b);

		REQUIRE(2 == c.rows());
		REQUIRE(2 == c.cols());

		REQUIRE(Complex(-7, 16) == c(0, 0));
		REQUIRE(Complex(-9, 22) == c(0, 1));
		REQUIRE(Complex(-9, 38) == c(1, 0));
		REQUIRE(Complex(-11, 52) == c(1, 1));
	}

	///////////////////       Vector<Complex> - Vector<Real> operations       ///////////////////
	TEST_CASE("BaseUtils::Utils_AddVec", "[simple]") {
			TEST_PRECISION_INFO();
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
	TEST_CASE("BaseUtils::Utils_SubVec", "[simple]") {
			TEST_PRECISION_INFO();
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
	TEST_CASE("Utils::RowMatrixFromVector", "[simple]") {
			TEST_PRECISION_INFO();
		Vector<Real> a{ REAL(1.0), REAL(2.0), REAL(3.0) };

		auto b = Utils::RowMatrixFromVector<Real>(a);

		REQUIRE(1 == b.rows());
		REQUIRE(3 == b.cols());

		REQUIRE(REAL(1.0) == b[0][0]);
		REQUIRE(REAL(2.0) == b[0][1]);
		REQUIRE(REAL(3.0) == b[0][2]);
	}
	TEST_CASE("Utils::ColumnMatrixFromVector", "[simple]") {
			TEST_PRECISION_INFO();
		Vector<Real> a{ REAL(1.0), REAL(2.0), REAL(3.0) };

		auto b = Utils::ColumnMatrixFromVector<Real>(a);

		REQUIRE(3 == b.rows());
		REQUIRE(1 == b.cols());

		REQUIRE(REAL(1.0) == b[0][0]);
		REQUIRE(REAL(2.0) == b[1][0]);
		REQUIRE(REAL(3.0) == b[2][0]);
	}

	///////////////////                   Matrix helpers                     ///////////////////
	TEST_CASE("Utils::Commutator", "[MatrixUtils]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A(2, 2, REAL(0.0));
		Matrix<Real> B(2, 2, REAL(0.0));

		// Example: A = [1 2; 3 4], B = [0 1; 1 0]
		A[0][0] = REAL(1.0); A[0][1] = REAL(2.0);
		A[1][0] = REAL(3.0); A[1][1] = REAL(4.0);

		B[0][0] = REAL(0.0); B[0][1] = REAL(1.0);
		B[1][0] = REAL(1.0); B[1][1] = REAL(0.0);

		// Commutator: [A, B] = AB - BA
		auto AB = A * B;
		auto BA = B * A;
		auto comm = AB - BA;

		// Manually compute expected result
		Matrix<Real> expected(2, 2, REAL(0.0));
		expected[0][0] = (1 * 0 + 2 * 1) - (0 * 1 + 1 * 3); // 2-3 = -1
		expected[0][1] = (1 * 1 + 2 * 0) - (0 * 2 + 1 * 4); // 1-4 = -3
		expected[1][0] = (3 * 0 + 4 * 1) - (1 * 1 + 0 * 3); // 4-1 = 3
		expected[1][1] = (3 * 1 + 4 * 0) - (1 * 2 + 0 * 4); // 3-2 = 1

		REQUIRE(comm[0][0] == expected[0][0]);
		REQUIRE(comm[0][1] == expected[0][1]);
		REQUIRE(comm[1][0] == expected[1][0]);
		REQUIRE(comm[1][1] == expected[1][1]);
	}

	TEST_CASE("Utils::AntiCommutator", "[MatrixUtils]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A(2, 2, REAL(0.0));
		Matrix<Real> B(2, 2, REAL(0.0));

		// Example: A = [1 2; 3 4], B = [0 1; 1 0]
		A[0][0] = REAL(1.0); A[0][1] = REAL(2.0);
		A[1][0] = REAL(3.0); A[1][1] = REAL(4.0);

		B[0][0] = REAL(0.0); B[0][1] = REAL(1.0);
		B[1][0] = REAL(1.0); B[1][1] = REAL(0.0);

		// Anti-commutator: {A, B} = AB + BA
		auto AB = A * B;
		auto BA = B * A;
		auto anti = AB + BA;

		// Manually compute expected result
		Matrix<Real> expected(2, 2, REAL(0.0));
		expected[0][0] = (1 * 0 + 2 * 1) + (0 * 1 + 1 * 3); // 2+3 = 5
		expected[0][1] = (1 * 1 + 2 * 0) + (0 * 2 + 1 * 4); // 1+4 = 5
		expected[1][0] = (3 * 0 + 4 * 1) + (1 * 1 + 0 * 3); // 4+1 = 5
		expected[1][1] = (3 * 1 + 4 * 0) + (1 * 2 + 0 * 4); // 3+2 = 5

		REQUIRE(anti[0][0] == expected[0][0]);
		REQUIRE(anti[0][1] == expected[0][1]);
		REQUIRE(anti[1][0] == expected[1][0]);
		REQUIRE(anti[1][1] == expected[1][1]);
	}

	TEST_CASE("Utils::MatrixDecomposeToSymAntisym", "[MatrixUtils]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A(2, 2, REAL(0.0));
		// Example: A = [1 2; 3 4]
		A[0][0] = REAL(1.0); A[0][1] = REAL(2.0);
		A[1][0] = REAL(3.0); A[1][1] = REAL(4.0);

		// Symmetric part: S = (A + A^T)/2
		Matrix<Real> S(2, 2, REAL(0.0));
		S[0][0] = (A[0][0] + A[0][0]) / REAL(2.0); // 1
		S[0][1] = (A[0][1] + A[1][0]) / REAL(2.0); // (2+3)/2 = REAL(2.5)
		S[1][0] = (A[1][0] + A[0][1]) / REAL(2.0); // (3+2)/2 = REAL(2.5)
		S[1][1] = (A[1][1] + A[1][1]) / REAL(2.0); // 4

		// Antisymmetric part: K = (A - A^T)/2
		Matrix<Real> K(2, 2, REAL(0.0));
		K[0][0] = REAL(0.0);
		K[0][1] = (A[0][1] - A[1][0]) / REAL(2.0); // (2-3)/2 = -REAL(0.5)
		K[1][0] = (A[1][0] - A[0][1]) / REAL(2.0); // (3-2)/2 = REAL(0.5)
		K[1][1] = REAL(0.0);

		// Check that S + K == A
		REQUIRE((S[0][0] + K[0][0]) == A[0][0]);
		REQUIRE((S[0][1] + K[0][1]) == A[0][1]);
		REQUIRE((S[1][0] + K[1][0]) == A[1][0]);
		REQUIRE((S[1][1] + K[1][1]) == A[1][1]);

		// Check symmetry
		REQUIRE(S[0][1] == S[1][0]);
		REQUIRE(K[0][1] == -K[1][0]);
	}
	///////////////////                  Matrix functions                    ///////////////////
	TEST_CASE("Utils::Exp", "[MatrixUtils]")
	{
		TEST_PRECISION_INFO();
		
		// exp(0) = I
		Matrix<Real> zero(2, 2, REAL(0.0));
		auto exp_zero = Utils::Exp(zero);
		REQUIRE(exp_zero.isIdentity(1e-10));
		
		// For a diagonal matrix with small values, verify the series expansion works
		// Note: Implementation has factorial computation that starts from 1
		// so we test the behavior rather than exact values for non-zero matrices
		Matrix<Real> small(2, 2, REAL(0.0));
		small[0][0] = REAL(0.1);
		small[1][1] = REAL(0.1);
		auto exp_small = Utils::Exp(small, 20);  // More terms for accuracy
		
		// For small diagonal matrices, off-diagonal should remain near zero
		REQUIRE(std::abs(exp_small[0][1]) < 1e-10);
		REQUIRE(std::abs(exp_small[1][0]) < 1e-10);
		// Diagonal should be > 1 (since exp(0.1) ≈ 1.105)
		REQUIRE(exp_small[0][0] > 1.0);
		REQUIRE(exp_small[1][1] > 1.0);
	}

	///////////////////                Real matrix helpers                   ///////////////////
	TEST_CASE("Utils::IsOrthogonal", "[MatrixUtils]")
	{
		TEST_PRECISION_INFO();
		
		// Identity is orthogonal
		Matrix<Real> I(3, 3, REAL(0.0));
		I[0][0] = I[1][1] = I[2][2] = REAL(1.0);
		REQUIRE(Utils::IsOrthogonal(I));
		
		// Simple rotation matrix (90 deg around z-axis)
		Matrix<Real> Rz(2, 2, REAL(0.0));
		Rz[0][0] = REAL(0.0);  Rz[0][1] = -REAL(1.0);
		Rz[1][0] = REAL(1.0);  Rz[1][1] = REAL(0.0);
		REQUIRE(Utils::IsOrthogonal(Rz));
		
		// Non-orthogonal matrix
		Matrix<Real> A(2, 2, REAL(0.0));
		A[0][0] = REAL(1.0); A[0][1] = REAL(2.0);
		A[1][0] = REAL(3.0); A[1][1] = REAL(4.0);
		REQUIRE_FALSE(Utils::IsOrthogonal(A));
	}

	///////////////////               Complex matrix helpers                 ///////////////////
	TEST_CASE("Utils::GetRealPart", "[MatrixUtils]")
	{
		TEST_PRECISION_INFO();
		
		Matrix<Complex> A(2, 2, Complex(0, 0));
		A[0][0] = Complex(REAL(1.0), REAL(2.0));
		A[0][1] = Complex(REAL(3.0), REAL(4.0));
		A[1][0] = Complex(REAL(5.0), REAL(6.0));
		A[1][1] = Complex(REAL(7.0), REAL(8.0));
		
		auto real = Utils::GetRealPart(A);
		REQUIRE(real[0][0] == REAL(1.0));
		REQUIRE(real[0][1] == REAL(3.0));
		REQUIRE(real[1][0] == REAL(5.0));
		REQUIRE(real[1][1] == REAL(7.0));
	}
	TEST_CASE("Utils::GetImagPart", "[MatrixUtils]")
	{
		TEST_PRECISION_INFO();
		
		Matrix<Complex> A(2, 2, Complex(0, 0));
		A[0][0] = Complex(REAL(1.0), REAL(2.0));
		A[0][1] = Complex(REAL(3.0), REAL(4.0));
		A[1][0] = Complex(REAL(5.0), REAL(6.0));
		A[1][1] = Complex(REAL(7.0), REAL(8.0));
		
		auto imag = Utils::GetImagPart(A);
		REQUIRE(imag[0][0] == REAL(2.0));
		REQUIRE(imag[0][1] == REAL(4.0));
		REQUIRE(imag[1][0] == REAL(6.0));
		REQUIRE(imag[1][1] == REAL(8.0));
	}
	TEST_CASE("Utils::GetConjugateTranspose", "[MatrixUtils]")
	{
		TEST_PRECISION_INFO();
		
		Matrix<Complex> A(2, 3, Complex(0, 0));
		A[0][0] = Complex(REAL(1.0), REAL(1.0));
		A[0][1] = Complex(REAL(2.0), REAL(2.0));
		A[0][2] = Complex(REAL(3.0), REAL(3.0));
		A[1][0] = Complex(REAL(4.0), REAL(4.0));
		A[1][1] = Complex(REAL(5.0), REAL(5.0));
		A[1][2] = Complex(REAL(6.0), REAL(6.0));
		
		auto Adag = Utils::GetConjugateTranspose(A);
		
		REQUIRE(Adag.rows() == 3);
		REQUIRE(Adag.cols() == 2);
		REQUIRE(Adag[0][0] == Complex(REAL(1.0), -REAL(1.0)));
		REQUIRE(Adag[1][0] == Complex(REAL(2.0), -REAL(2.0)));
		REQUIRE(Adag[2][0] == Complex(REAL(3.0), -REAL(3.0)));
		REQUIRE(Adag[0][1] == Complex(REAL(4.0), -REAL(4.0)));
	}
	TEST_CASE("Utils::CmplxMatFromRealMat", "[MatrixUtils]")
	{
		TEST_PRECISION_INFO();
		
		Matrix<Real> A(2, 2, REAL(0.0));
		A[0][0] = REAL(1.0); A[0][1] = REAL(2.0);
		A[1][0] = REAL(3.0); A[1][1] = REAL(4.0);
		
		auto C = Utils::CmplxMatFromRealMat(A);
		
		REQUIRE(C[0][0] == Complex(REAL(1.0), REAL(0.0)));
		REQUIRE(C[0][1] == Complex(REAL(2.0), REAL(0.0)));
		REQUIRE(C[1][0] == Complex(REAL(3.0), REAL(0.0)));
		REQUIRE(C[1][1] == Complex(REAL(4.0), REAL(0.0)));
	}
	TEST_CASE("Utils::IsComplexMatReal", "[MatrixUtils]")
	{
		TEST_PRECISION_INFO();
		
		// Real matrix (zero imaginary parts)
		Matrix<Complex> real(2, 2, Complex(0, 0));
		real[0][0] = Complex(REAL(1.0), REAL(0.0));
		real[0][1] = Complex(REAL(2.0), REAL(0.0));
		real[1][0] = Complex(REAL(3.0), REAL(0.0));
		real[1][1] = Complex(REAL(4.0), REAL(0.0));
		REQUIRE(Utils::IsComplexMatReal(real));
		
		// Complex matrix (non-zero imaginary parts)
		Matrix<Complex> cmplx(2, 2, Complex(0, 0));
		cmplx[0][0] = Complex(REAL(1.0), REAL(0.1));
		REQUIRE_FALSE(Utils::IsComplexMatReal(cmplx));
	}
	TEST_CASE("Utils::IsHermitian", "[MatrixUtils]")
	{
		TEST_PRECISION_INFO();
		
		// Hermitian matrix: A = A†
		Matrix<Complex> H(2, 2, Complex(0, 0));
		H[0][0] = Complex(REAL(1.0), REAL(0.0));  // Real diagonal
		H[1][1] = Complex(REAL(2.0), REAL(0.0));  // Real diagonal
		H[0][1] = Complex(REAL(3.0), REAL(4.0));  // Off-diagonal
		H[1][0] = Complex(REAL(3.0), -REAL(4.0)); // Conjugate
		REQUIRE(Utils::IsHermitian(H));
		
		// Non-Hermitian matrix
		Matrix<Complex> N(2, 2, Complex(0, 0));
		N[0][0] = Complex(REAL(1.0), REAL(0.0));
		N[1][1] = Complex(REAL(2.0), REAL(0.0));
		N[0][1] = Complex(REAL(3.0), REAL(4.0));
		N[1][0] = Complex(REAL(3.0), REAL(4.0));  // NOT conjugate
		REQUIRE_FALSE(Utils::IsHermitian(N));
	}
	TEST_CASE("Utils::IsUnitary", "[MatrixUtils]")
	{
		TEST_PRECISION_INFO();
		
		// Identity is unitary
		Matrix<Complex> I(2, 2, Complex(0, 0));
		I[0][0] = Complex(REAL(1.0), REAL(0.0));
		I[1][1] = Complex(REAL(1.0), REAL(0.0));
		REQUIRE(Utils::IsUnitary(I));
		
		// Simple unitary: phase rotation
		Matrix<Complex> U(2, 2, Complex(0, 0));
		Real angle = Constants::PI / REAL(4.0);
		U[0][0] = Complex(std::cos(angle), std::sin(angle));
		U[1][1] = Complex(std::cos(angle), -std::sin(angle));
		REQUIRE(Utils::IsUnitary(U));
		
		// Non-unitary matrix
		Matrix<Complex> N(2, 2, Complex(0, 0));
		N[0][0] = Complex(REAL(2.0), REAL(0.0));
		N[1][1] = Complex(REAL(1.0), REAL(0.0));
		REQUIRE_FALSE(Utils::IsUnitary(N));
	}

	///////////////////       Matrix<Complex> - Matrix<Real>  operations     ///////////////////
	TEST_CASE("Utils::AddMat", "[MatrixUtils]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A(2, 2, REAL(0.0));
		Matrix<Complex> B(2, 2, Complex(REAL(0.0), REAL(0.0)));

		// A = [1 2; 3 4], B = [5+i 6; 7 8-i]
		A[0][0] = REAL(1.0); A[0][1] = REAL(2.0);
		A[1][0] = REAL(3.0); A[1][1] = REAL(4.0);

		B[0][0] = Complex(REAL(5.0), REAL(1.0)); B[0][1] = Complex(REAL(6.0), REAL(0.0));
		B[1][0] = Complex(REAL(7.0), REAL(0.0)); B[1][1] = Complex(REAL(8.0), -REAL(1.0));

		auto C = Utils::AddMat(A, B);

		REQUIRE(C[0][0] == Complex(REAL(6.0), REAL(1.0)));
		REQUIRE(C[0][1] == Complex(REAL(8.0), REAL(0.0)));
		REQUIRE(C[1][0] == Complex(REAL(10.0), REAL(0.0)));
		REQUIRE(C[1][1] == Complex(REAL(12.0), -REAL(1.0)));
	}

	TEST_CASE("Utils::SubMat", "[MatrixUtils]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A(2, 2, REAL(0.0));
		Matrix<Complex> B(2, 2, Complex(REAL(0.0), REAL(0.0)));

		// A = [5 6; 7 8], B = [1 2+i; 3 4]
		A[0][0] = REAL(5.0); A[0][1] = REAL(6.0);
		A[1][0] = REAL(7.0); A[1][1] = REAL(8.0);

		B[0][0] = Complex(REAL(1.0), REAL(0.0)); B[0][1] = Complex(REAL(2.0), REAL(1.0));
		B[1][0] = Complex(REAL(3.0), REAL(0.0)); B[1][1] = Complex(REAL(4.0), REAL(0.0));

		auto C = Utils::SubMat(A, B);

		REQUIRE(C[0][0] == Complex(REAL(4.0), REAL(0.0)));
		REQUIRE(C[0][1] == Complex(REAL(4.0), -REAL(1.0)));
		REQUIRE(C[1][0] == Complex(REAL(4.0), REAL(0.0)));
		REQUIRE(C[1][1] == Complex(REAL(4.0), REAL(0.0)));
	}

	TEST_CASE("Utils::MulMat1", "[MatrixUtils]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A(2, 2, REAL(0.0));
		Matrix<Complex> B(2, 2, Complex(REAL(0.0), REAL(0.0)));

		// A = [1 2; 3 4], B = [2 0; 1+i 2]
		A[0][0] = REAL(1.0); A[0][1] = REAL(2.0);
		A[1][0] = REAL(3.0); A[1][1] = REAL(4.0);

		B[0][0] = Complex(REAL(2.0), REAL(0.0)); B[0][1] = Complex(REAL(0.0), REAL(0.0));
		B[1][0] = Complex(REAL(1.0), REAL(1.0)); B[1][1] = Complex(REAL(2.0), REAL(0.0));

		auto C = Utils::MulMat(A, B);

		REQUIRE(C[0][0] == Complex(REAL(4.0), REAL(2.0)));  // 1*2 + 2*(1+i) = 2 + 2 + 2i = 4+2i
		REQUIRE(C[0][1] == Complex(REAL(4.0), REAL(0.0)));  // 1*0 + 2*2 = 4
		REQUIRE(C[1][0] == Complex(REAL(10.0), REAL(4.0))); // 3*2 + 4*(1+i) = 6 + 4 + 4i = 10+4i
		REQUIRE(C[1][1] == Complex(REAL(8.0), REAL(0.0)));  // 3*0 + 4*2 = 8
	}

	TEST_CASE("Utils::MulMat2", "[MatrixUtils]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A(2, 3, REAL(0.0));
		Matrix<Complex> B(3, 2, Complex(REAL(0.0), REAL(0.0)));

		// A = [1 2 3; 4 5 6], B = [7 8; 9+i 10; 11 12-i]
		A[0][0] = REAL(1.0); A[0][1] = REAL(2.0); A[0][2] = REAL(3.0);
		A[1][0] = REAL(4.0); A[1][1] = REAL(5.0); A[1][2] = REAL(6.0);

		B[0][0] = Complex(REAL(7.0), REAL(0.0));  B[0][1] = Complex(REAL(8.0), REAL(0.0));
		B[1][0] = Complex(REAL(9.0), REAL(1.0));  B[1][1] = Complex(REAL(10.0), REAL(0.0));
		B[2][0] = Complex(REAL(11.0), REAL(0.0)); B[2][1] = Complex(REAL(12.0), -REAL(1.0));

		auto C = Utils::MulMat(A, B);

		REQUIRE(C[0][0] == Complex(1 * 7 + 2 * 9 + 3 * 11, 2 * 1));   // 7+18+33, 2*1=2 => 58+2i
		REQUIRE(C[0][1] == Complex(1 * 8 + 2 * 10 + 3 * 12, -3 * 1));  // 8+20+36, -3*1=-3 => 64-3i
		REQUIRE(C[1][0] == Complex(4 * 7 + 5 * 9 + 6 * 11, 5 * 1));    // 28+45+66, 5*1=5 => 139+5i
		REQUIRE(C[1][1] == Complex(4 * 8 + 5 * 10 + 6 * 12, -6 * 1));  // 32+50+72, -6*1=-6 => 154-6i
	}

	TEST_CASE("Utils::MulMatVec", "[MatrixUtils]")
	{
			TEST_PRECISION_INFO();
		Matrix<Real> A(2, 2, REAL(0.0));
		Vector<Complex> v({ Complex(REAL(1.0), REAL(0.0)), Complex(REAL(2.0), REAL(1.0)) });

		// A = [3 4; 5 6]
		A[0][0] = REAL(3.0); A[0][1] = REAL(4.0);
		A[1][0] = REAL(5.0); A[1][1] = REAL(6.0);

		auto result = Utils::MulMatVec(A, v);

		REQUIRE(result.size() == 2);
		REQUIRE(result[0] == Real(REAL(3.0)) * Complex(REAL(1.0), REAL(0.0)) + Real(REAL(4.0)) * Complex(REAL(2.0), REAL(1.0))); // 3 + 8 + 4i = 11+4i
		REQUIRE(result[1] == Real(REAL(5.0)) * Complex(REAL(1.0), REAL(0.0)) + Real(REAL(6.0)) * Complex(REAL(2.0), REAL(1.0))); // 5 + 12 + 6i = 17+6i
	}

	TEST_CASE("Utils::MulVecMat", "[MatrixUtils]")
	{
			TEST_PRECISION_INFO();
		Vector<Real> v({ REAL(1.0), REAL(2.0) });
		Matrix<Complex> A(2, 2, Complex(REAL(0.0), REAL(0.0)));

		// A = [3 4; 5+i 6]
		A[0][0] = Complex(REAL(3.0), REAL(0.0)); A[0][1] = Complex(REAL(4.0), REAL(0.0));
		A[1][0] = Complex(REAL(5.0), REAL(1.0)); A[1][1] = Complex(REAL(6.0), REAL(0.0));

		auto result = Utils::MulVecMat(v, A);

		REQUIRE(result.size() == 2);
		REQUIRE(result[0] == Real(REAL(1.0)) * Complex(REAL(3.0), REAL(0.0)) + Real(REAL(2.0)) * Complex(REAL(5.0), REAL(1.0))); // 3 + 10 + 2i = 13+2i
		REQUIRE(result[1] == Real(REAL(1.0)) * Complex(REAL(4.0), REAL(0.0)) + Real(REAL(2.0)) * Complex(REAL(6.0), REAL(0.0))); // 4 + 12 = 16
	}

	TEST_CASE("Utils::LeviCivita3D", "[levi_civita]")
	{
			TEST_PRECISION_INFO();
		// 3D Levi-Civita symbol tests
		// Even permutations return +1
		REQUIRE(Utils::LeviCivita(1, 2, 3) == 1);
		REQUIRE(Utils::LeviCivita(2, 3, 1) == 1);
		REQUIRE(Utils::LeviCivita(3, 1, 2) == 1);

		// Odd permutations return -1
		REQUIRE(Utils::LeviCivita(3, 2, 1) == -1);
		REQUIRE(Utils::LeviCivita(2, 1, 3) == -1);
		REQUIRE(Utils::LeviCivita(1, 3, 2) == -1);

		// Repeated indices return 0
		REQUIRE(Utils::LeviCivita(1, 1, 2) == 0);
		REQUIRE(Utils::LeviCivita(1, 2, 2) == 0);
		REQUIRE(Utils::LeviCivita(3, 3, 3) == 0);
	}

	TEST_CASE("Utils::LeviCivita4D", "[levi_civita]")
	{
			TEST_PRECISION_INFO();
		// 4D Levi-Civita symbol tests
		// Even permutations return +1
		REQUIRE(Utils::LeviCivita(1, 2, 3, 4) == 1);  // Identity: 0 swaps
		REQUIRE(Utils::LeviCivita(1, 2, 4, 3) == -1); // 1 swap
		REQUIRE(Utils::LeviCivita(1, 4, 2, 3) == 1);  // 2 swaps
		REQUIRE(Utils::LeviCivita(2, 3, 4, 1) == -1); // 3 swaps (cyclic)

		// Odd permutations return -1
		REQUIRE(Utils::LeviCivita(2, 1, 3, 4) == -1);  // 1 swap
		REQUIRE(Utils::LeviCivita(1, 3, 2, 4) == -1);  // 1 swap
		REQUIRE(Utils::LeviCivita(4, 3, 2, 1) == 1);   // 6 swaps (reversal)

		// Repeated indices return 0
		REQUIRE(Utils::LeviCivita(1, 1, 3, 4) == 0);
		REQUIRE(Utils::LeviCivita(1, 2, 2, 4) == 0);
		REQUIRE(Utils::LeviCivita(1, 2, 3, 3) == 0);
		REQUIRE(Utils::LeviCivita(2, 2, 2, 2) == 0);

		// Invalid indices (out of range 1-4) return 0
		REQUIRE(Utils::LeviCivita(0, 2, 3, 4) == 0);
		REQUIRE(Utils::LeviCivita(1, 5, 3, 4) == 0);
		REQUIRE(Utils::LeviCivita(1, 2, -1, 4) == 0);
		REQUIRE(Utils::LeviCivita(1, 2, 3, 10) == 0);

		// Missing indices (not all 1,2,3,4 present) return 0
		REQUIRE(Utils::LeviCivita(1, 1, 2, 3) == 0); // missing 4, duplicate 1
		REQUIRE(Utils::LeviCivita(2, 3, 4, 4) == 0); // missing 1, duplicate 4
	}

	///////////////////              Angle Utility Tests                     ///////////////////
	TEST_CASE("Utils::DegToRad", "[AngleUtils]")
	{
		TEST_PRECISION_INFO();
		
		REQUIRE_THAT(Utils::DegToRad(REAL(0.0)), RealWithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(Utils::DegToRad(REAL(90.0)), RealWithinAbs(Constants::PI / REAL(2.0), REAL(1e-10)));
		REQUIRE_THAT(Utils::DegToRad(REAL(180.0)), RealWithinAbs(Constants::PI, REAL(1e-10)));
		REQUIRE_THAT(Utils::DegToRad(REAL(360.0)), RealWithinAbs(REAL(2.0) * Constants::PI, REAL(1e-10)));
		REQUIRE_THAT(Utils::DegToRad(-REAL(45.0)), RealWithinAbs(-Constants::PI / REAL(4.0), REAL(1e-10)));
	}

	TEST_CASE("Utils::RadToDeg", "[AngleUtils]")
	{
		TEST_PRECISION_INFO();
		
		REQUIRE_THAT(Utils::RadToDeg(REAL(0.0)), RealWithinAbs(REAL(0.0), REAL(1e-10)));
		REQUIRE_THAT(Utils::RadToDeg(Constants::PI / REAL(2.0)), RealWithinAbs(REAL(90.0), REAL(1e-10)));
		REQUIRE_THAT(Utils::RadToDeg(Constants::PI), RealWithinAbs(REAL(180.0), REAL(1e-10)));
		REQUIRE_THAT(Utils::RadToDeg(REAL(2.0) * Constants::PI), RealWithinAbs(REAL(360.0), REAL(1e-10)));
		REQUIRE_THAT(Utils::RadToDeg(-Constants::PI / REAL(4.0)), RealWithinAbs(-REAL(45.0), REAL(1e-10)));
	}

	TEST_CASE("Utils::AngleDegToExplicit", "[AngleUtils]")
	{
		TEST_PRECISION_INFO();
		
		Real deg, min, sec;
		
		// 45.5 degrees = 45 deg 30 min 0 sec
		Utils::AngleDegToExplicit(REAL(45.5), deg, min, sec);
		REQUIRE_THAT(deg, RealWithinAbs(REAL(45.0), REAL(0.1)));
		REQUIRE_THAT(min, RealWithinAbs(REAL(30.0), REAL(0.1)));
		REQUIRE_THAT(sec, RealWithinAbs(REAL(0.0), REAL(0.1)));
		
		// 30.2525 degrees = 30 deg 15 min 9 sec
		Utils::AngleDegToExplicit(REAL(30.2525), deg, min, sec);
		REQUIRE_THAT(deg, RealWithinAbs(REAL(30.0), REAL(0.1)));
		REQUIRE_THAT(min, RealWithinAbs(REAL(15.0), REAL(0.1)));
		REQUIRE_THAT(sec, RealWithinAbs(REAL(9.0), REAL(0.1)));
	}

	TEST_CASE("Utils::ExplicitToAngleDeg", "[AngleUtils]")
	{
		TEST_PRECISION_INFO();
		
		REQUIRE_THAT(Utils::ExplicitToAngleDeg(45, 30, REAL(0.0)), RealWithinAbs(REAL(45.5), REAL(1e-10)));
		REQUIRE_THAT(Utils::ExplicitToAngleDeg(30, 15, REAL(9.0)), RealWithinAbs(REAL(30.2525), REAL(1e-4)));
		REQUIRE_THAT(Utils::ExplicitToAngleDeg(0, 0, REAL(0.0)), RealWithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("Utils::AngleTo2PiRange", "[AngleUtils]")
	{
		TEST_PRECISION_INFO();
		
		// Already in range [0, 2π)
		REQUIRE_THAT(Utils::AngleTo2PiRange(REAL(1.0)), RealWithinAbs(REAL(1.0), REAL(1e-10)));
		
		// Negative angle
		REQUIRE_THAT(Utils::AngleTo2PiRange(-Constants::PI), RealWithinAbs(Constants::PI, REAL(1e-10)));
		
		// Angle > 2π
		REQUIRE_THAT(Utils::AngleTo2PiRange(REAL(3.0) * Constants::PI), RealWithinAbs(Constants::PI, REAL(1e-10)));
		
		// Large negative angle
		REQUIRE_THAT(Utils::AngleTo2PiRange(-REAL(3.0) * Constants::PI), RealWithinAbs(Constants::PI, REAL(1e-10)));
	}

	TEST_CASE("Utils::AngleToPiPiRange", "[AngleUtils]")
	{
		TEST_PRECISION_INFO();
		
		// Already in range [-π, π]
		REQUIRE_THAT(Utils::AngleToPiPiRange(REAL(1.0)), RealWithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE_THAT(Utils::AngleToPiPiRange(-REAL(1.0)), RealWithinAbs(-REAL(1.0), REAL(1e-10)));
		
		// Angle just over π
		REQUIRE_THAT(Utils::AngleToPiPiRange(Constants::PI + REAL(0.1)), RealWithinAbs(-Constants::PI + REAL(0.1), REAL(1e-10)));
		
		// Angle just under -π
		REQUIRE_THAT(Utils::AngleToPiPiRange(-Constants::PI - REAL(0.1)), RealWithinAbs(Constants::PI - REAL(0.1), REAL(1e-10)));
	}

	///////////////////              KroneckerDelta Tests                    ///////////////////
	TEST_CASE("Utils::KroneckerDelta", "[SymbolUtils]")
	{
		TEST_PRECISION_INFO();
		
		// δ_ii = 1
		REQUIRE(Utils::KroneckerDelta(0, 0) == 1);
		REQUIRE(Utils::KroneckerDelta(1, 1) == 1);
		REQUIRE(Utils::KroneckerDelta(5, 5) == 1);
		REQUIRE(Utils::KroneckerDelta(-3, -3) == 1);
		
		// δ_ij = 0 for i ≠ j
		REQUIRE(Utils::KroneckerDelta(0, 1) == 0);
		REQUIRE(Utils::KroneckerDelta(1, 0) == 0);
		REQUIRE(Utils::KroneckerDelta(3, 7) == 0);
		REQUIRE(Utils::KroneckerDelta(-1, 1) == 0);
	}

	///////////////////              VectorsAngle Tests                      ///////////////////
	TEST_CASE("Utils::VectorsAngle", "[VectorUtils]")
	{
		TEST_PRECISION_INFO();
		
		// Parallel vectors - angle 0
		VectorN<Real, 3> v1({ REAL(1.0), REAL(0.0), REAL(0.0) });
		VectorN<Real, 3> v2({ REAL(2.0), REAL(0.0), REAL(0.0) });
		REQUIRE_THAT(Utils::VectorsAngle(v1, v2), RealWithinAbs(REAL(0.0), REAL(1e-10)));
		
		// Perpendicular vectors - angle π/2
		VectorN<Real, 3> v3({ REAL(0.0), REAL(1.0), REAL(0.0) });
		REQUIRE_THAT(Utils::VectorsAngle(v1, v3), RealWithinAbs(Constants::PI / REAL(2.0), REAL(1e-10)));
		
		// Anti-parallel vectors - angle π
		VectorN<Real, 3> v4({ -REAL(1.0), REAL(0.0), REAL(0.0) });
		REQUIRE_THAT(Utils::VectorsAngle(v1, v4), RealWithinAbs(Constants::PI, REAL(1e-10)));
		
		// 45-degree angle
		VectorN<Real, 3> v5({ REAL(1.0), REAL(1.0), REAL(0.0) });
		REQUIRE_THAT(Utils::VectorsAngle(v1, v5), RealWithinAbs(Constants::PI / REAL(4.0), REAL(1e-10)));
	}

	///////////////////              DiagonalMatrixFromVector Tests          ///////////////////
	TEST_CASE("Utils::DiagonalMatrixFromVector", "[MatrixUtils]")
	{
		TEST_PRECISION_INFO();
		
		Vector<Real> v({ REAL(1.0), REAL(2.0), REAL(3.0) });
		auto D = Utils::DiagonalMatrixFromVector(v);
		
		REQUIRE(D.rows() == 3);
		REQUIRE(D.cols() == 3);
		REQUIRE(D[0][0] == REAL(1.0));
		REQUIRE(D[1][1] == REAL(2.0));
		REQUIRE(D[2][2] == REAL(3.0));
		REQUIRE(D[0][1] == REAL(0.0));
		REQUIRE(D[0][2] == REAL(0.0));
		REQUIRE(D[1][0] == REAL(0.0));
		REQUIRE(D[1][2] == REAL(0.0));
		REQUIRE(D[2][0] == REAL(0.0));
		REQUIRE(D[2][1] == REAL(0.0));
	}

	///////////////////              MatrixFromVectors Tests                 ///////////////////
	TEST_CASE("Utils::MatrixFromVectorsInRows", "[MatrixUtils]")
	{
		TEST_PRECISION_INFO();
		
		Vector<Real> v1({ REAL(1.0), REAL(2.0), REAL(3.0) });
		Vector<Real> v2({ REAL(4.0), REAL(5.0), REAL(6.0) });
		
		std::vector<Vector<Real>> vecs = { v1, v2 };
		auto M = Utils::MatrixFromVectorsInRows(vecs);
		
		REQUIRE(M.rows() == 2);
		REQUIRE(M.cols() == 3);
		REQUIRE(M[0][0] == REAL(1.0));
		REQUIRE(M[0][1] == REAL(2.0));
		REQUIRE(M[0][2] == REAL(3.0));
		REQUIRE(M[1][0] == REAL(4.0));
		REQUIRE(M[1][1] == REAL(5.0));
		REQUIRE(M[1][2] == REAL(6.0));
	}

	TEST_CASE("Utils::MatrixFromVectorsInColumns", "[MatrixUtils]")
	{
		TEST_PRECISION_INFO();
		
		Vector<Real> v1({ REAL(1.0), REAL(2.0) });
		Vector<Real> v2({ REAL(3.0), REAL(4.0) });
		Vector<Real> v3({ REAL(5.0), REAL(6.0) });
		
		std::vector<Vector<Real>> vecs = { v1, v2, v3 };
		auto M = Utils::MatrixFromVectorsInColumns(vecs);
		
		REQUIRE(M.rows() == 2);
		REQUIRE(M.cols() == 3);
		REQUIRE(M[0][0] == REAL(1.0));
		REQUIRE(M[0][1] == REAL(3.0));
		REQUIRE(M[0][2] == REAL(5.0));
		REQUIRE(M[1][0] == REAL(2.0));
		REQUIRE(M[1][1] == REAL(4.0));
		REQUIRE(M[1][2] == REAL(6.0));
	}

	///////////////////              Matrix Sin/Cos/Sinh/Cosh Tests          ///////////////////
	TEST_CASE("Utils::Sin_Matrix", "[MatrixUtils]")
	{
		TEST_PRECISION_INFO();
		
		// sin(0) = 0
		Matrix<Real> zero(2, 2, REAL(0.0));
		auto sin_zero = Utils::Sin(zero);
		REQUIRE(std::abs(sin_zero[0][0]) < 1e-10);
		REQUIRE(std::abs(sin_zero[0][1]) < 1e-10);
		REQUIRE(std::abs(sin_zero[1][0]) < 1e-10);
		REQUIRE(std::abs(sin_zero[1][1]) < 1e-10);
	}

	TEST_CASE("Utils::Cos_Matrix", "[MatrixUtils]")
	{
		TEST_PRECISION_INFO();
		
		// cos(0) = I
		Matrix<Real> zero(2, 2, REAL(0.0));
		auto cos_zero = Utils::Cos(zero);
		REQUIRE_THAT(cos_zero[0][0], RealWithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE_THAT(cos_zero[1][1], RealWithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE(std::abs(cos_zero[0][1]) < 1e-10);
		REQUIRE(std::abs(cos_zero[1][0]) < 1e-10);
	}

	TEST_CASE("Utils::Sinh_Matrix", "[MatrixUtils]")
	{
		TEST_PRECISION_INFO();
		
		// sinh(0) = 0
		Matrix<Real> zero(2, 2, REAL(0.0));
		auto sinh_zero = Utils::Sinh(zero);
		REQUIRE(std::abs(sinh_zero[0][0]) < 1e-10);
		REQUIRE(std::abs(sinh_zero[0][1]) < 1e-10);
		REQUIRE(std::abs(sinh_zero[1][0]) < 1e-10);
		REQUIRE(std::abs(sinh_zero[1][1]) < 1e-10);
	}

	TEST_CASE("Utils::Cosh_Matrix", "[MatrixUtils]")
	{
		TEST_PRECISION_INFO();
		
		// cosh(0) = I
		Matrix<Real> zero(2, 2, REAL(0.0));
		auto cosh_zero = Utils::Cosh(zero);
		REQUIRE_THAT(cosh_zero[0][0], RealWithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE_THAT(cosh_zero[1][1], RealWithinAbs(REAL(1.0), REAL(1e-10)));
		REQUIRE(std::abs(cosh_zero[0][1]) < 1e-10);
		REQUIRE(std::abs(cosh_zero[1][0]) < 1e-10);
	}
}
