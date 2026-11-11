#ifndef MML_USE_SINGLE_HEADER
#include <mml/algorithms/MatrixAlg/MatrixInversion.h>
#endif
// Naming-convention regression tests for the matrix family (CONTRACT_POLICY section 11,
// task MinimalMathLibrary-wewi.5.2).  Every matrix type must expose the same vocabulary:
//   STL protocol lowercase (rows/cols/empty/resize), mathematical operations PascalCase
//   (Trace/Transposed/Inverse/Norm*), verb mutates / noun returns a copy (Invert/Inverse,
//   Transpose/Transposed), reset vocabulary SetZero/SetIdentity/Fill, accessors getX().
// The deprecated 2.0 spellings must keep compiling with identical semantics for one release.
#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/Matrix/MatrixSym.h>
#include <mml/base/Matrix/MatrixTriDiag.h>
#include <mml/base/Matrix/MatrixBandDiag.h>
#include <mml/base/Matrix/Matrix3D.h>
#include <mml/base/SparseMatrix/SparseMatrixCOO.h>
#include <mml/base/SparseMatrix/SparseMatrixCSR.h>
#include <mml/base/SparseMatrix/SparseMatrixCSC.h>
#include <mml/base/SparseMatrix/SparseMatrix.h>
#endif

#include <type_traits>

using namespace MML;

namespace MML::Tests::Base::MatrixNamingTests
{
	// ---- compile-time vocabulary checks -------------------------------------------------------
	template<typename M> using rows_t       = decltype(std::declval<const M&>().rows());
	template<typename M> using cols_t       = decltype(std::declval<const M&>().cols());
	template<typename M> using trace_t      = decltype(std::declval<const M&>().Trace());
	template<typename M> using transposed_t = decltype(std::declval<const M&>().Transposed());
	template<typename M> using inverse_t    = decltype(std::declval<const M&>().Inverse());
	template<typename M> using invert_t     = decltype(std::declval<M&>().Invert());
	template<typename M> using transpose_t  = decltype(std::declval<M&>().Transpose());
	template<typename M> using setzero_t    = decltype(std::declval<M&>().SetZero());
	template<typename M> using setid_t      = decltype(std::declval<M&>().SetIdentity());
	template<typename M> using fill_t       = decltype(std::declval<M&>().Fill(Real{1}));
	template<typename M> using tofull_t     = decltype(std::declval<const M&>().ToFullMatrix());
	template<typename M> using normfro_t    = decltype(std::declval<const M&>().NormFrobenius());
	template<typename M> using lower_t      = decltype(std::declval<const M&>().transpose());

	template<template<class> class Op, class M, class = void> struct has : std::false_type {};
	template<template<class> class Op, class M> struct has<Op, M, std::void_t<Op<M>>> : std::true_type {};

	using M  = Matrix<Real>;
	using M3 = MatrixNM<Real, 3, 3>;
	using MS = MatrixSym<Real>;
	using MT = TridiagonalMatrix<Real>;
	using MB = BandDiagonalMatrix;
	using CSR = SparseMatrix::SparseMatrixCSR<Real>;

	// rows()/cols() on every 2-D type
	static_assert(has<rows_t, M>::value && has<cols_t, M>::value);
	static_assert(has<rows_t, M3>::value && has<cols_t, M3>::value);
	static_assert(has<rows_t, MS>::value && has<cols_t, MS>::value);
	static_assert(has<rows_t, MT>::value && has<cols_t, MT>::value);
	static_assert(has<rows_t, MB>::value && has<cols_t, MB>::value);
	static_assert(has<rows_t, CSR>::value && has<cols_t, CSR>::value);

	// Trace() / Transposed() / Inverse() are the copy-returning spellings
	static_assert(has<trace_t, M>::value && has<trace_t, M3>::value && has<trace_t, MS>::value && has<trace_t, MB>::value);
	static_assert(has<transposed_t, M>::value && has<transposed_t, M3>::value && has<transposed_t, MT>::value && has<transposed_t, MB>::value);
	static_assert(!has<inverse_t, M>::value && !has<inverse_t, M3>::value && !has<inverse_t, MS>::value && !has<inverse_t, MT>::value);

	// the verb mutates and returns void
	static_assert(!has<invert_t, M>::value && !has<invert_t, M3>::value);
	static_assert(std::is_void_v<transpose_t<M>> && std::is_void_v<transpose_t<M3>>);

	// reset vocabulary
	static_assert(std::is_void_v<setzero_t<M>> && std::is_void_v<setzero_t<M3>> && std::is_void_v<setzero_t<MS>>);
	static_assert(std::is_void_v<setid_t<M>> && std::is_void_v<setid_t<M3>>);
	static_assert(std::is_void_v<fill_t<M>> && std::is_void_v<fill_t<M3>> && std::is_void_v<fill_t<MS>>);

	// conversions and norms
	static_assert(has<tofull_t, MS>::value && has<tofull_t, MT>::value && has<tofull_t, MB>::value);
	static_assert(has<normfro_t, MS>::value && has<normfro_t, MB>::value && has<normfro_t, CSR>::value);

	// the copy-returning spellings are [[nodiscard]]
	static_assert(std::is_same_v<transposed_t<M3>, MatrixNM<Real, 3, 3>>);

	TEST_CASE("MatrixNaming_Matrix_verb_mutates_noun_copies", "[Matrix][naming]")
	{
		Matrix<Real> a(2, 2, { 4.0, 7.0, 2.0, 6.0 });
		const Matrix<Real> original = a;

		Matrix<Real> t = a.Transposed();
		REQUIRE(a.IsEqualTo(original));
		REQUIRE(t(0, 1) == 2.0);
		a.Transpose();
		REQUIRE(a.IsEqualTo(t));

		Matrix<Real> inv = MatrixAlg::Inverse(original);
		REQUIRE(inv(0, 0) == Catch::Approx(0.6));
		REQUIRE(inv(0, 1) == Catch::Approx(-0.7));
		Matrix<Real> b = original;
		MatrixAlg::Invert(b);
		REQUIRE(b.IsEqualTo(inv, 1e-14));

		REQUIRE(original.Trace() == 10.0);
		REQUIRE(original.getDiagonal()[1] == 6.0);
		REQUIRE(original.IsSquare());
		REQUIRE_FALSE(original.empty());
		REQUIRE_FALSE(original.IsZero());
		REQUIRE_FALSE(original.IsIdentity());
		REQUIRE_FALSE(original.IsSymmetric());
		REQUIRE_FALSE(original.IsAntiSymmetric());

		Matrix<Real> s(2, 2, { 1.0, 2.0, 2.0, 1.0 });
		REQUIRE(s.IsSymmetric());
		REQUIRE_FALSE(Matrix<Real>(2, 3).IsSquare());

		Matrix<Real> r(2, 2);
		r.SetIdentity();
		REQUIRE(r.IsIdentity());
		r.Fill(3.0);
		REQUIRE(r(1, 0) == 3.0);
		r.SetZero();
		REQUIRE(r.IsZero(0.0));
		r.resize(3, 1);
		REQUIRE(r.rows() == 3);
		REQUIRE(r.cols() == 1);
		REQUIRE(Matrix<Real>().empty());
	}

	TEST_CASE("MatrixNaming_MatrixNM_vocabulary", "[MatrixNM][naming]")
	{
		MatrixNM<Real, 2, 3> a({ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 });
		MatrixNM<Real, 3, 2> t = a.Transposed();
		REQUIRE(t(2, 0) == 3.0);
		REQUIRE(t(0, 1) == 4.0);

		MatrixNM<Real, 2, 2> m({ 4.0, 7.0, 2.0, 6.0 });
		MatrixNM<Real, 2, 2> inv = MatrixAlg::Inverse(m);
		REQUIRE(inv(1, 1) == Catch::Approx(0.4));
		MatrixAlg::Invert(m);
		REQUIRE(m.IsEqualTo(inv, 1e-14));

		m.SetIdentity();
		REQUIRE(m(0, 0) == 1.0);
		REQUIRE(m(0, 1) == 0.0);
		m.Fill(2.5);
		REQUIRE(m(1, 0) == 2.5);
		m.SetZero();
		REQUIRE(m.Trace() == 0.0);
	}

	TEST_CASE("MatrixNaming_MatrixSym_vocabulary", "[MatrixSym][naming]")
	{
		MatrixSym<Real> s(2, { 2.0, 1.0, 3.0 });   // lower triangle: (0,0)=2, (1,0)=1, (1,1)=3
		REQUIRE(s.rows() == 2);
		REQUIRE(s.cols() == 2);

		Matrix<Real> full = s.ToFullMatrix();
		REQUIRE(full(0, 1) == full(1, 0));
		REQUIRE(full.IsSymmetric());

		Matrix<Real> inv = MatrixAlg::Inverse(s);
		Matrix<Real> prod = full * inv;
		REQUIRE(prod.IsIdentity(1e-12));

		MatrixSym<Real> d = MatrixSym<Real>::Diagonal(Vector<Real>{ 1.0, 2.0, 3.0 });
		REQUIRE(d.Trace() == 6.0);
		REQUIRE(d(0, 1) == 0.0);

		s.Fill(7.0);
		REQUIRE(s(1, 0) == 7.0);
		s.SetZero();
		REQUIRE(s.NormFrobenius() == 0.0);
		s.resize(4);
		REQUIRE(s.rows() == 4);
	}

	TEST_CASE("MatrixNaming_TriDiag_BandDiag_vocabulary", "[TridiagonalMatrix][BandDiag][naming]")
	{
		// band row-wise: {d0, a0,  b1, d1, a1,  b2, d2}
		TridiagonalMatrix<Real> tri(3, { 4.0, 7.0,  2.0, 5.0, 8.0,  1.0, 6.0 });
		TridiagonalMatrix<Real> tt = tri.Transposed();
		REQUIRE(tt(0, 1) == tri(1, 0));
		REQUIRE(tt(1, 2) == tri(2, 1));

		Matrix<Real> full = tri.ToFullMatrix();
		REQUIRE(full(0, 1) == 7.0);
		REQUIRE(full(1, 0) == 2.0);
		REQUIRE(full(0, 2) == 0.0);
		Matrix<Real> prod = full * MatrixAlg::Inverse(tri);
		REQUIRE(prod.IsIdentity(1e-12));

		Matrix<Real> data(3, 3, { 0.0, 1.0, 2.0,
		                          3.0, 4.0, 5.0,
		                          6.0, 7.0, 0.0 });
		BandDiagonalMatrix band(3, 1, 1, data);
		REQUIRE(band.rows() == 3);
		REQUIRE(band.getLowerBandwidth() == 1);
		REQUIRE(band.getUpperBandwidth() == 1);
		REQUIRE(band.getDiagonal()[1] == 4.0);
		REQUIRE(band.getDiagonal(1)[0] == 3.0);

		BandDiagonalMatrix bt = band.Transposed();
		REQUIRE(bt(1, 0) == band(0, 1));
		REQUIRE(bt.ToFullMatrix().IsEqualTo(band.ToFullMatrix().Transposed()));
		REQUIRE(band.Trace() == 12.0);
	}

	TEST_CASE("MatrixNaming_Sparse_vocabulary", "[Sparse][naming]")
	{
		auto L = SparseMatrix::Laplacian1D<Real>(4);
		REQUIRE(L.IsSquare());
		REQUIRE(L.IsSymmetric());
		REQUIRE(L.HasSymmetricPattern());
		REQUIRE(L.NormInf() == Catch::Approx(4.0));
		REQUIRE(L.Norm1() == Catch::Approx(4.0));
		REQUIRE(L.NormFrobenius() == Catch::Approx(std::sqrt(4 * 4.0 + 6 * 1.0)));

		std::vector<Real> x{ 1.0, 1.0, 1.0, 1.0 }, y;
		L.Multiply(x, y);
		REQUIRE(y[0] == 1.0);
		REQUIRE(y[1] == 0.0);
		L.MultiplyTranspose(x, y);
		REQUIRE(y[3] == 1.0);

		auto I = SparseMatrix::IdentityCSR<Real>(3);
		std::vector<Real> b{ 1.0, 2.0, 3.0 }, sol;
		I.SolveLower(b, sol);
		REQUIRE(sol == b);
		I.SolveUpper(b, sol);
		REQUIRE(sol == b);
		I.ZeroRow(1);
		REQUIRE(I.ToDense()[1][1] == 0.0);

		auto D = SparseMatrix::DiagonalCSR<Real>(std::vector<Real>{ 2.0, 3.0 });
		REQUIRE(D(1, 1) == 3.0);

		SparseMatrix::SparseMatrixCOO<Real> coo(2, 2);
		coo.addEntry(0, 1, 5.0);
		SparseMatrix::SparseMatrixCSC<Real> csc;
		csc.FromCOO(coo);
		REQUIRE(csc.ToDense()[0][1] == 5.0);
		csc.Multiply(std::vector<Real>{ 0.0, 1.0 }, y);
		REQUIRE(y[0] == 5.0);
	}

	// ---- deprecated 2.0 spellings must keep compiling with identical semantics ---------------
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4996)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
	TEST_CASE("MatrixNaming_deprecated_aliases_forward", "[Matrix][naming][deprecated]")
	{
		Matrix<Real> a(2, 2, { 4.0, 7.0, 2.0, 6.0 });
		REQUIRE(a.trace() == a.Trace());
		REQUIRE(a.transpose().IsEqualTo(a.Transposed()));
		REQUIRE(a.isSquare() == a.IsSquare());
		REQUIRE(a.isEmpty() == a.empty());
		REQUIRE(a.isZero() == a.IsZero());
		REQUIRE(a.isIdentity() == a.IsIdentity());
		REQUIRE(a.isAntiSymmetric() == a.IsAntiSymmetric());
		REQUIRE(a.diagonal()[0] == a.getDiagonal()[0]);
		Matrix<Real> u(2, 2);
		u.MakeUnitMatrix();
		REQUIRE(u.IsIdentity());
		u.Resize(3, 3);
		REQUIRE(u.rows() == 3);

		MatrixNM<Real, 2, 2> m({ 4.0, 7.0, 2.0, 6.0 });
		REQUIRE(m.transpose().IsEqualTo(m.Transposed()));
		m.MakeUnitMatrix();
		REQUIRE(m.Trace() == 2.0);

		MatrixSym<Real> s(2, { 2.0, 1.0, 3.0 });
		REQUIRE(s.Dim() == s.rows());
		REQUIRE(s.GetAsMatrix().IsEqualTo(s.ToFullMatrix()));
		REQUIRE(MatrixSym<Real>::GetDiagonalMatrix(Vector<Real>{ 1.0, 2.0 }).Trace() == 3.0);
		s.SetToValue(1.0);
		REQUIRE(s(1, 1) == 1.0);
		s.SetToZero();
		REQUIRE(s(1, 1) == 0.0);
		s.Resize(3);
		REQUIRE(s.rows() == 3);

		TridiagonalMatrix<Real> tri(2, { 2.0, 3.0,  1.0, 4.0 });   // {d0, a0, b1, d1}
		REQUIRE(tri.GetTranspose()(0, 1) == tri.Transposed()(0, 1));

		Matrix<Real> data(2, 2, { 0.0, 1.0, 2.0, 3.0 });
		BandDiagonalMatrix band(2, 1, 0, data);
		REQUIRE(band.GetDimension() == band.rows());
		REQUIRE(band.GetLowerBandwidth() == band.getLowerBandwidth());
		REQUIRE(band.GetUpperBandwidth() == band.getUpperBandwidth());
		REQUIRE(band.GetDiagonal()[1] == band.getDiagonal()[1]);
		REQUIRE(band.Transpose()(0, 1) == band.Transposed()(0, 1));   // copy semantics preserved

		Vector<Real> v(2);
		v.Resize(5);
		REQUIRE(v.size() == 5);
		Matrix3D<Real> m3(1, 1, 1);
		m3.Resize(2, 2, 2);
		REQUIRE(m3.dim1() == 2);

		auto L = SparseMatrix::laplacian1D<Real>(3);
		REQUIRE(L.isSymmetric() == L.IsSymmetric());
		REQUIRE(L.norm1() == L.Norm1());
		REQUIRE(L.normInf() == L.NormInf());
		REQUIRE(L.normFrobenius() == L.NormFrobenius());
		REQUIRE(L.isSquare());
		REQUIRE(L.hasSymmetricPattern());
		std::vector<Real> x{ 1.0, 2.0, 3.0 }, y1, y2;
		L.multiply(x, y1);
		L.Multiply(x, y2);
		REQUIRE(y1 == y2);
		REQUIRE(L.toDense() == L.ToDense());
		REQUIRE(SparseMatrix::identityCSR<Real>(2).nnz() == 2);
		REQUIRE(SparseMatrix::diagonalCSR<Real>(std::vector<Real>{ 1.0, 2.0 })(1, 1) == 2.0);
		REQUIRE(SparseMatrix::laplacian2D<Real>(2, 2).rows() == 4);
	}
#if defined(_MSC_VER)
#pragma warning(pop)
#elif defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif
}
