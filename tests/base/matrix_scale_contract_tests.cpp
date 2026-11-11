#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/Matrix/MatrixSym.h>
#include <mml/base/Matrix/MatrixTriDiag.h>
#include <mml/base/Matrix/MatrixBandDiag.h>
#include <mml/base/SparseMatrix/SparseMatrixCSR.h>
#endif

#include <array>
#include <complex>
#include <limits>

namespace MML::Tests::Base::MatrixScaleContracts
{
	TEMPLATE_TEST_CASE("Tridiagonal equality rejects nonfinite stored components",
		"[matrix-scale][tridiag][tridiag-equality]", Real, Complex)
	{
		const TridiagonalMatrix<TestType> zero(3, Vector<TestType>(3), Vector<TestType>(3), Vector<TestType>(3));
		for (const Real invalid : {std::numeric_limits<Real>::quiet_NaN(),
			std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity()}) {
			for (int band = 0; band < 3; ++band)
				for (int slot = 0; slot < 3; ++slot) {
					for (int component = 0; component < (MMLComplex<TestType> ? 2 : 1); ++component) {
						std::array<Vector<TestType>, 3> bands{Vector<TestType>(3), Vector<TestType>(3), Vector<TestType>(3)};
						if constexpr (MMLComplex<TestType>)
							bands[band][slot] = component == 0 ? TestType(invalid, 0) : TestType(0, invalid);
						else
							bands[band][slot] = invalid;
						const TridiagonalMatrix<TestType> bad(3, bands[0], bands[1], bands[2]);
						CHECK_FALSE(bad.IsEqualTo(zero));
						CHECK_FALSE(zero.IsEqualTo(bad));
						CHECK_FALSE(bad.IsEqualTo(bad));
						CHECK_FALSE(bad.IsEqualTo(zero, std::numeric_limits<Real>::max()));
					}
				}
		}
	}

	TEMPLATE_TEST_CASE("Tridiagonal equality validates tolerance and dimensions",
		"[matrix-scale][tridiag][tridiag-equality]", Real, Complex)
	{
		const TridiagonalMatrix<TestType> zero(1, {TestType{0}});
		const TridiagonalMatrix<TestType> differentSize(2, {TestType{0}, TestType{0}, TestType{0}, TestType{0}});
		CHECK(zero.IsEqualTo(zero, 0));
		CHECK_FALSE(zero.IsEqualTo(differentSize));
		CHECK_FALSE(differentSize.IsEqualTo(zero));
		for (const Real invalid : {Real{-1}, -std::numeric_limits<Real>::infinity(),
			std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()})
			CHECK_FALSE(zero.IsEqualTo(zero, invalid));
		const TridiagonalMatrix<TestType> closeValue(1, {TestType(Precision::MatrixIsEqualTolerance / Real{2})});
		CHECK(zero.IsEqualTo(closeValue));
		CHECK_FALSE(zero.IsEqualTo(closeValue, 0));
	}

	TEMPLATE_TEST_CASE("Tridiagonal equality retains inclusive absolute magnitude tolerance",
		"[matrix-scale][tridiag][tridiag-equality]", Real, Complex)
	{
		for (const int exponent : {std::max(-1000, std::numeric_limits<Real>::min_exponent / 2),
			std::max(-100, std::numeric_limits<Real>::min_exponent / 4), 0,
			std::min(100, std::numeric_limits<Real>::max_exponent / 4),
			std::min(1000, std::numeric_limits<Real>::max_exponent / 2)}) {
			const Real scale = std::scalbn(Real{1}, exponent);
			const TridiagonalMatrix<TestType> first(2, {TestType{0}, TestType(scale), TestType(-scale), TestType{0}});
			const TridiagonalMatrix<TestType> second(2, {TestType(scale), TestType(scale), TestType(-scale), TestType{0}});
			CHECK(first.IsEqualTo(second, scale));
			CHECK(second.IsEqualTo(first, scale));
			CHECK_FALSE(first.IsEqualTo(second, std::nextafter(scale, Real{0})));
			CHECK_FALSE(second.IsEqualTo(first, std::nextafter(scale, Real{0})));
			for (int band = 0; band < 3; ++band)
				for (int slot = 0; slot < 3; ++slot) {
					std::array<Vector<TestType>, 3> bands{Vector<TestType>(3), Vector<TestType>(3), Vector<TestType>(3)};
					const TridiagonalMatrix<TestType> zero(3, bands[0], bands[1], bands[2]);
					bands[band][slot] = TestType(-scale);
					const TridiagonalMatrix<TestType> value(3, bands[0], bands[1], bands[2]);
					CHECK(zero.IsEqualTo(value, scale));
					CHECK(value.IsEqualTo(zero, scale));
					CHECK_FALSE(zero.IsEqualTo(value, std::nextafter(scale, Real{0})));
				}
		}
		const Real tiny = std::numeric_limits<Real>::denorm_min();
		const Real huge = std::numeric_limits<Real>::max();
		const TridiagonalMatrix<TestType> zero(1, {TestType{0}}), subnormal(1, {TestType(tiny)});
		const TridiagonalMatrix<TestType> large(1, {TestType(huge)}), negativeLarge(1, {TestType(-huge)});
		CHECK(zero.IsEqualTo(subnormal, tiny));
		CHECK_FALSE(zero.IsEqualTo(subnormal, 0));
		CHECK(large.IsEqualTo(large, 0));
		CHECK_FALSE(large.IsEqualTo(negativeLarge, huge));
		CHECK_FALSE(negativeLarge.IsEqualTo(large, huge));
		if constexpr (MMLComplex<TestType>) {
			const TridiagonalMatrix<TestType> value(1, {TestType(3, 4)});
			CHECK(zero.IsEqualTo(value, 5));
			CHECK_FALSE(zero.IsEqualTo(value, 4));
		}
	}

	TEMPLATE_TEST_CASE("Tridiagonal integral equality avoids overflow wrapping and premature rounding",
		"[matrix-scale][tridiag][tridiag-equality]", int, unsigned, long long, unsigned long long)
	{
		const TestType low = std::numeric_limits<TestType>::lowest();
		const TestType high = std::numeric_limits<TestType>::max();
		const TridiagonalMatrix<TestType> first(1, {low}), last(1, {high}), adjacent(1, {TestType(high - 1)});
		CHECK(first.IsEqualTo(first, 0));
		CHECK_FALSE(first.IsEqualTo(last, 1));
		CHECK_FALSE(last.IsEqualTo(first, 1));
		CHECK_FALSE(last.IsEqualTo(adjacent, 0));
		CHECK(last.IsEqualTo(adjacent, 1));
		CHECK(adjacent.IsEqualTo(last, 1));
		CHECK(first.IsEqualTo(last, std::numeric_limits<Real>::max()));
		for (const Real invalid : {Real{-1}, std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()})
			CHECK_FALSE(first.IsEqualTo(first, invalid));
		if constexpr (std::numeric_limits<TestType>::digits >= 54) {
			const TestType exact = TestType{1} << 53;
			const TridiagonalMatrix<TestType> zero(1, {TestType{0}});
			const TridiagonalMatrix<TestType> boundary(1, {exact}), above(1, {TestType(exact + 1)});
			const Real tolerance = static_cast<Real>(exact);
			CHECK(zero.IsEqualTo(boundary, tolerance));
			CHECK_FALSE(zero.IsEqualTo(above, tolerance));
			CHECK_FALSE(above.IsEqualTo(zero, tolerance));
		}
	}

	template<class T>
	void CheckPredicatesReject(const Matrix<T>& matrix, Real tolerance)
	{
		CHECK_FALSE(matrix.IsZero(tolerance));
		CHECK_FALSE(matrix.IsIdentity(tolerance));
		CHECK_FALSE(matrix.IsSymmetric(tolerance));
		CHECK_FALSE(matrix.IsAntiSymmetric(tolerance));
		CHECK_FALSE(matrix.IsHermitian(tolerance));
		CHECK_FALSE(matrix.IsSkewHermitian(tolerance));
	}

	TEMPLATE_TEST_CASE("Matrix member predicates include diagonal constraints",
		"[matrix-scale][Matrix][A09]", Real, Complex, int, unsigned)
	{
		Matrix<TestType> one(1, 1, {TestType{2}});
		CHECK_FALSE(one.IsAntiSymmetric(0));
		CHECK_FALSE(one.IsSkewHermitian(0));
		CHECK(one.IsSymmetric(0));
		CHECK(one.IsHermitian(0));

		Matrix<TestType> zero(2, 2);
		CHECK(zero.IsZero(0));
		CHECK(zero.IsAntiSymmetric(0));
		CHECK(zero.IsSkewHermitian(0));
		zero.SetIdentity();
		CHECK(zero.IsIdentity(0));
		CHECK_FALSE(zero.IsAntiSymmetric(0));
		CHECK_FALSE(zero.IsSkewHermitian(0));

		if constexpr (!std::is_unsigned_v<TestType>) {
			Matrix<TestType> skew(2, 2, {TestType{0}, TestType{3}, TestType{-3}, TestType{0}});
			CHECK(skew.IsAntiSymmetric(0));
			CHECK(skew.IsSkewHermitian(0));
			CHECK_FALSE(skew.IsSymmetric(0));
			skew(1, 1) = TestType{1};
			CHECK_FALSE(skew.IsAntiSymmetric(0));
			CHECK_FALSE(skew.IsSkewHermitian(0));
		}
	}

	TEST_CASE("Complex antisymmetry and skew Hermiticity have different diagonals",
		"[matrix-scale][Matrix][complex][A09]")
	{
		Matrix<Complex> imaginary(2, 2, {Complex(0, 3), Complex(2, 5), Complex(-2, 5), Complex(0, -4)});
		CHECK(imaginary.IsSkewHermitian(0));
		CHECK_FALSE(imaginary.IsAntiSymmetric(0));
		CHECK_FALSE(imaginary.IsHermitian(0));
		const auto adjoint = imaginary.Transposed();
		for (int i = 0; i < imaginary.rows(); ++i)
			for (int j = 0; j < imaginary.cols(); ++j)
				CHECK(imaginary(i, j) == -std::conj(adjoint(i, j)));

		imaginary(0, 0) = Complex(Precision::EqualityTolerance, 3);
		CHECK(imaginary.IsSkewHermitian(Precision::EqualityTolerance));
		CHECK_FALSE(imaginary.IsSkewHermitian(0));

		Matrix<Complex> plainSkew(2, 2, {Complex(0), Complex(2, 5), Complex(-2, -5), Complex(0)});
		CHECK(plainSkew.IsAntiSymmetric(0));
		CHECK_FALSE(plainSkew.IsSkewHermitian(0));
	}

	TEMPLATE_TEST_CASE("Matrix approximate predicates fail closed for nonfinite entries",
		"[matrix-scale][Matrix][MatrixNM][MatrixSym][A18]", Real, Complex)
	{
		for (const Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
			for (int position = 0; position < 4; ++position) {
				Matrix<TestType> matrix(2, 2);
				matrix(position / 2, position % 2) = TestType(invalid);
				CheckPredicatesReject(matrix, Precision::EqualityTolerance);
				CHECK_FALSE(matrix.IsEqualTo(matrix));
				CHECK_FALSE(Matrix<TestType>::AreEqual(matrix, matrix));
			}
			MatrixNM<TestType, 1, 1> fixed{TestType(invalid)};
			MatrixSym<TestType> symmetric(1, {TestType(invalid)});
			CHECK_FALSE(fixed.IsEqualTo(fixed));
			CHECK_FALSE(symmetric.IsEqualTo(symmetric));
			if constexpr (MMLComplex<TestType>) {
				Matrix<TestType> imaginaryNaN(1, 1, {TestType(0, invalid)});
				CheckPredicatesReject(imaginaryNaN, Precision::EqualityTolerance);
			}
		}
	}

	TEMPLATE_TEST_CASE("Matrix comparison tolerance is inclusive and must be finite nonnegative",
		"[matrix-scale][Matrix][MatrixNM][MatrixSym][A18]", Real, Complex, int, unsigned)
	{
		Matrix<TestType> a(1, 1, {TestType{0}}), b(1, 1, {TestType{1}});
		MatrixNM<TestType, 1, 1> fixedA{TestType{0}}, fixedB{TestType{1}};
		MatrixSym<TestType> symA(1, {TestType{0}}), symB(1, {TestType{1}});
		CHECK(a.IsEqualTo(b, 1));
		CHECK(b.IsEqualTo(a, 1));
		CHECK(fixedA.IsEqualTo(fixedB, 1));
		CHECK(symA.IsEqualTo(symB, 1));
		CHECK_FALSE(a.IsEqualTo(b, 0));
		CHECK_FALSE(fixedA.IsEqualTo(fixedB, 0));
		CHECK_FALSE(symA.IsEqualTo(symB, 0));

		for (const Real invalid : {Real{-1}, std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
			CheckPredicatesReject(a, invalid);
			CheckPredicatesReject(Matrix<TestType>(), invalid);
			CHECK_FALSE(a.IsEqualTo(a, invalid));
			CHECK_FALSE(fixedA.IsEqualTo(fixedA, invalid));
			CHECK_FALSE(symA.IsEqualTo(symA, invalid));
			CHECK_FALSE(MatrixSym<TestType>().IsEqualTo(MatrixSym<TestType>(), invalid));
		}
	}

	TEMPLATE_TEST_CASE("Integral matrix member comparisons do not wrap or subtract in signed arithmetic",
		"[matrix-scale][Matrix][MatrixNM][MatrixSym][A10]", int, unsigned, long long, unsigned long long)
	{
		const TestType low = std::numeric_limits<TestType>::lowest();
		const TestType high = std::numeric_limits<TestType>::max();
		Matrix<TestType> a(1, 1, {low}), b(1, 1, {high});
		MatrixNM<TestType, 1, 1> fixedA{low}, fixedB{high};
		MatrixSym<TestType> symA(1, {low}), symB(1, {high});
		CHECK_FALSE(a.IsEqualTo(b, 1));
		CHECK_FALSE(b.IsEqualTo(a, 1));
		CHECK_FALSE(fixedA.IsEqualTo(fixedB, 1));
		CHECK_FALSE(fixedB.IsEqualTo(fixedA, 1));
		CHECK_FALSE(symA.IsEqualTo(symB, 1));
		CHECK_FALSE(symB.IsEqualTo(symA, 1));
		b(0, 0) = high - 1;
		a(0, 0) = high;
		CHECK_FALSE(a.IsEqualTo(b, 0));
		CHECK(a.IsEqualTo(b, 1));

		Matrix<TestType> asymmetric(2, 2, {TestType{0}, low, high, TestType{0}});
		CHECK_FALSE(asymmetric.IsSymmetric(1));
		if constexpr (std::is_unsigned_v<TestType>) {
			Matrix<TestType> wrapped(2, 2, {TestType{0}, TestType{1}, high, TestType{0}});
			CHECK_FALSE(wrapped.IsAntiSymmetric(0));
			CHECK_FALSE(wrapped.IsSkewHermitian(0));
		}
	}

	TEST_CASE("Structured and sparse Frobenius norms retain extreme finite scales",
		"[matrix-scale][MatrixSym][BandDiag][SparseMatrix][A14]")
	{
		for (const Real scale : {Real{1}, std::numeric_limits<Real>::max() / 64,
				std::numeric_limits<Real>::min() * 64, -std::numeric_limits<Real>::min() * 64}) {
			CAPTURE(scale);
			MatrixSym<Real> symmetric(2, {3 * scale, 4 * scale, Real{0}});
			CHECK(symmetric.NormFrobenius() / std::abs(scale) ==
				Catch::Approx(std::sqrt(Real{41})).epsilon(Precision::RoundoffZeroThreshold));
			MatrixSym<Complex> complexSymmetric(2, {Complex(3 * scale, 4 * scale), Complex(5 * scale, 12 * scale), Complex(0)});
			CHECK(complexSymmetric.NormFrobenius() / std::abs(scale) ==
				Catch::Approx(std::sqrt(Real{363})).epsilon(Precision::RoundoffZeroThreshold));

			BandDiagonalMatrix band(2, 1, 1, Matrix<Real>(2, 3));
			band(0, 0) = 3 * scale;
			band(0, 1) = 4 * scale;
			CHECK(band.NormFrobenius() / std::abs(scale) ==
				Catch::Approx(5).epsilon(Precision::RoundoffZeroThreshold));
			const BandDiagonalMatrix transposed = band.Transposed();
			CHECK(transposed.NormFrobenius() == band.NormFrobenius());

			SparseMatrix::SparseMatrixCSR<Real> sparse(2, 2,
				{3 * scale, 4 * scale}, {0, 1}, {0, 2, 2});
			CHECK(sparse.NormFrobenius() / std::abs(scale) ==
				Catch::Approx(5).epsilon(Precision::RoundoffZeroThreshold));
			SparseMatrix::SparseMatrixCSR<Complex> complexSparse(2, 2,
				{Complex(3 * scale, 4 * scale), Complex(0, 12 * scale)}, {0, 1}, {0, 2, 2});
			CHECK(complexSparse.NormFrobenius() / std::abs(scale) ==
				Catch::Approx(13).epsilon(Precision::RoundoffZeroThreshold));
		}
		CHECK(MatrixSym<Real>().NormFrobenius() == 0);
		CHECK(MatrixSym<int>(1, {std::numeric_limits<int>::lowest()}).NormFrobenius() ==
			-Real(std::numeric_limits<int>::lowest()));
		CHECK(BandDiagonalMatrix(2, 0, 0, Matrix<Real>(2, 1)).NormFrobenius() == 0);
		CHECK(SparseMatrix::SparseMatrixCSR<Real>(0, 0).NormFrobenius() == 0);
		const Real subnormal = std::numeric_limits<Real>::denorm_min();
		CHECK(MatrixSym<Real>(1, {subnormal}).NormFrobenius() == subnormal);
		CHECK(MatrixSym<Complex>(1, {Complex(3 * subnormal, 4 * subnormal)}).NormFrobenius() == 5 * subnormal);
		CHECK(SparseMatrix::SparseMatrixCSR<Complex>(1, 1,
			{Complex(3 * subnormal, 4 * subnormal)}, {0}, {0, 1}).NormFrobenius() == 5 * subnormal);
		CHECK(SparseMatrix::SparseMatrixCSR<int>(2, 2, {30000, 40000}, {0, 1}, {0, 2, 2}).NormFrobenius() == 50000);
	}

	TEST_CASE("Scaled Frobenius accumulation preserves infinity and NaN propagation",
		"[matrix-scale][MatrixSym][BandDiag][SparseMatrix][A14]")
	{
		const Real infinity = std::numeric_limits<Real>::infinity();
		const Real nan = std::numeric_limits<Real>::quiet_NaN();
		MatrixSym<Real> symmetric(2, {infinity, -infinity, infinity});
		CHECK(std::isinf(symmetric.NormFrobenius()));
		symmetric(1, 1) = nan;
		CHECK(std::isnan(symmetric.NormFrobenius()));
		CHECK(std::isnan(MatrixSym<Complex>(1, {Complex(infinity, nan)}).NormFrobenius()));
		Matrix<Real> diagonal(2, 1, {infinity, -infinity});
		CHECK(std::isinf(BandDiagonalMatrix(2, 0, 0, diagonal).NormFrobenius()));
		diagonal(1, 0) = nan;
		CHECK(std::isnan(BandDiagonalMatrix(2, 0, 0, diagonal).NormFrobenius()));
		CHECK(std::isinf(SparseMatrix::SparseMatrixCSR<Real>(1, 2,
			{infinity, -infinity}, {0, 1}, {0, 2}).NormFrobenius()));
		CHECK(std::isnan(SparseMatrix::SparseMatrixCSR<Complex>(1, 2,
			{Complex(infinity), Complex(infinity, nan)}, {0, 1}, {0, 2}).NormFrobenius()));
	}

	TEMPLATE_TEST_CASE("Tridiagonal solve is invariant under uniform finite scaling",
		"[matrix-scale][tridiag][solve][A15]", Real, Complex)
	{
		for (const Real scale : {Real{1}, Real(1e-20), std::numeric_limits<Real>::min() * 64,
				std::numeric_limits<Real>::min() / 64, std::numeric_limits<Real>::denorm_min() * 64,
				std::numeric_limits<Real>::max() / 128, -Real(1e-20)}) {
			CAPTURE(scale);
			TestType factor(scale);
			if constexpr (MMLComplex<TestType>) factor = TestType(scale, scale / 2);
			const std::array<TestType, 4> expected{TestType{1}, TestType{-2}, TestType{3}, TestType{-4}};
			TridiagonalMatrix<TestType> matrix(4,
				{TestType{4} * factor, TestType{2} * factor, -factor, TestType{4} * factor, TestType{2} * factor,
				 -factor, TestType{4} * factor, TestType{2} * factor, -factor, TestType{4} * factor});
			Vector<TestType> rhs{TestType{0} * factor, TestType{-3} * factor, TestType{6} * factor, TestType{-19} * factor};
			const auto solution = matrix.Solve(rhs);
			for (int row = 0; row < 4; ++row) {
				CHECK(Abs(solution[row] - expected[row]) <= Precision::RoundoffZeroThreshold);
				TestType product{};
				Real bound = Abs(rhs[row]);
				for (int column = 0; column < 4; ++column) {
					const TestType entry = static_cast<const TridiagonalMatrix<TestType>&>(matrix)(row, column);
					product += entry * solution[column];
					bound += Abs(entry) * Abs(solution[column]);
				}
				CHECK(Abs(product - rhs[row]) / bound <= Precision::RoundoffZeroThreshold);
			}
			Vector<TestType> aliased = rhs;
			matrix.Solve(aliased, aliased);
			for (int row = 0; row < 4; ++row)
				CHECK(Abs(aliased[row] - expected[row]) <= Precision::RoundoffZeroThreshold);
		}
	}

	TEST_CASE("Tridiagonal scalar pivots include subnormals but reject exact and elimination zeros",
		"[matrix-scale][tridiag][solve][A15]")
	{
		for (const Real scale : {Real(1e-20), std::numeric_limits<Real>::min(),
				std::numeric_limits<Real>::denorm_min(), std::numeric_limits<Real>::max()}) {
			CAPTURE(scale);
			TridiagonalMatrix<Real> scalar(1, {scale});
			CHECK(scalar.Solve(Vector<Real>{scale})[0] == 1);
			TridiagonalMatrix<Complex> complexScalar(1, {Complex(scale, scale)});
			CHECK(complexScalar.Solve(Vector<Complex>{Complex(scale, scale)})[0] == Complex(1));
			TridiagonalMatrix<Real> singular(2, {scale, scale, scale, scale});
			CHECK_THROWS_AS(singular.Solve(Vector<Real>{scale, scale}), SingularMatrixError);
		}
		TridiagonalMatrix<Real> zero(1, {0});
		Vector<Real> output{7};
		CHECK_THROWS_AS(zero.Solve(Vector<Real>{1}, output), SingularMatrixError);
		CHECK(output[0] == 7);
		CHECK_THROWS_AS(zero.Solve(Vector<Real>{}), VectorDimensionError);
		CHECK_THROWS_AS(zero.Solve(Vector<Real>{1}, output = Vector<Real>{}), VectorDimensionError);

		TridiagonalMatrix<Real> needsPivoting(2, {0, 1, 1, 0});
		CHECK_THROWS_AS(needsPivoting.Solve(Vector<Real>{1, 1}), SingularMatrixError);
		for (const Real scale : {Real(1e-20), Real{1}, std::numeric_limits<Real>::max() / 4}) {
			TridiagonalMatrix<Real> negligiblePivot(2,
				{scale, Real{0}, Real{0}, scale * (Precision::Eps / 2)});
			CHECK_THROWS_AS(negligiblePivot.Solve(Vector<Real>{scale, scale}), SingularMatrixError);
			TridiagonalMatrix<Real> boundaryPivot(2,
				{scale, Real{0}, Real{0}, scale * Precision::Eps});
			const auto boundarySolution = boundaryPivot.Solve(Vector<Real>{scale, scale * Precision::Eps});
			CHECK(boundarySolution[0] == 1);
			CHECK(boundarySolution[1] == 1);
		}
		for (const Real invalid : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
			TridiagonalMatrix<Real> nonfinite(1, {invalid});
			CHECK_THROWS_AS(nonfinite.Solve(Vector<Real>{1}), NumericInputError);
			TridiagonalMatrix<Real> finite(1, {1});
			CHECK_THROWS_AS(finite.Solve(Vector<Real>{invalid}), NumericInputError);
		}
	}

	TEST_CASE("Tridiagonal complex normalization preserves extreme diagonal solves and aliases",
		"[matrix-scale][tridiag][solve]")
	{
		for (Real scale : {std::numeric_limits<Real>::denorm_min(), std::numeric_limits<Real>::min(),
			std::numeric_limits<Real>::max()}) {
			CAPTURE(scale);
			const TridiagonalMatrix<Complex> matrix(2,
				{Complex(scale, scale), Complex(0), Complex(0), Complex(-scale, scale)});
			const Vector<Complex> rhs{Complex(scale, scale), Complex(-scale, -scale)};
			const auto solution = matrix.Solve(rhs);
			CHECK(solution[0] == Complex(1));
			CHECK(solution[1] == Complex(0, 1));
			auto alias = rhs;
			matrix.Solve(alias, alias);
			CHECK(alias == solution);
		}
	}

	TEST_CASE("Tridiagonal normalization preserves existing exact integer diagonal solves",
		"[matrix-scale][tridiag][solve][A10]")
	{
		const TridiagonalMatrix<int> diagonal(2,
			Vector<int>{0, 0}, Vector<int>{2, 4}, Vector<int>{0, 0});
		const auto solution = diagonal.Solve(Vector<int>{2, 4});
		CHECK(solution[0] == 1);
		CHECK(solution[1] == 1);
	}
}
