#include <catch2/catch_test_macros.hpp>

#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Matrix/MatrixSym.h>
#include <mml/base/Matrix/Matrix3D.h>
#include <mml/base/SparseMatrix/SparseMatrix.h>
#include <mml/base/Tensor.h>

#include <limits>
#include <new>
#include <utility>

namespace {
struct ThrowingScalar {
	int value = 0;
	static inline bool fail = false;

	ThrowingScalar() {
		if (fail) throw std::bad_alloc();
	}
	ThrowingScalar(int v) : value(v) {}
	ThrowingScalar(const ThrowingScalar& other) : value(other.value) {
		if (fail) throw std::bad_alloc();
	}
	ThrowingScalar& operator=(const ThrowingScalar& other) {
		if (fail) throw std::bad_alloc();
		value = other.value;
		return *this;
	}
	ThrowingScalar operator+(const ThrowingScalar& other) const { return value + other.value; }
	ThrowingScalar operator-(const ThrowingScalar& other) const { return value - other.value; }
	ThrowingScalar operator*(const ThrowingScalar& other) const { return value * other.value; }
	ThrowingScalar operator/(const ThrowingScalar& other) const { return value / other.value; }
};

struct FailCopies {
	FailCopies() { ThrowingScalar::fail = true; }
	~FailCopies() { ThrowingScalar::fail = false; }
};
}

TEST_CASE("Dense storage preserves shape and data on failed replacement", "[storage-regression][storage-safe-red]") {
	using namespace MML;
	SECTION("Matrix copy assignment") {
		Matrix<ThrowingScalar> destination(1, 1, ThrowingScalar{7});
		const Matrix<ThrowingScalar> source(2, 2, ThrowingScalar{9});
		{
			FailCopies failure;
			REQUIRE_THROWS_AS(destination = source, std::bad_alloc);
		}
		REQUIRE(destination.rows() == 1);
		REQUIRE(destination.cols() == 1);
		REQUIRE(destination(0, 0).value == 7);
	}
	SECTION("MatrixSym resize and copy assignment") {
		MatrixSym<ThrowingScalar> destination(1, ThrowingScalar{7});
		const MatrixSym<ThrowingScalar> source(2, ThrowingScalar{9});
		{
			FailCopies failure;
			REQUIRE_THROWS_AS(destination.resize(2), std::bad_alloc);
		}
		REQUIRE(destination.rows() == 1);
		REQUIRE(destination.size() == 1);
		REQUIRE(destination(0, 0).value == 7);
		{
			FailCopies failure;
			REQUIRE_THROWS_AS(destination = source, std::bad_alloc);
		}
		REQUIRE(destination.rows() == 1);
		REQUIRE(destination(0, 0).value == 7);
	}
	SECTION("Matrix3D resize and copy assignment") {
		Matrix3D<ThrowingScalar> destination(1, 1, 1, ThrowingScalar{7});
		const Matrix3D<ThrowingScalar> source(2, 2, 2, ThrowingScalar{9});
		{
			FailCopies failure;
			REQUIRE_THROWS_AS(destination.resize(2, 2, 2, ThrowingScalar{9}), std::bad_alloc);
		}
		REQUIRE(destination.dim1() == 1);
		REQUIRE(destination.dim2() == 1);
		REQUIRE(destination.dim3() == 1);
		REQUIRE(destination[0][0][0].value == 7);
		{
			FailCopies failure;
			REQUIRE_THROWS_AS(destination.resize(2, 2, 2), std::bad_alloc);
		}
		REQUIRE(destination.size() == 1);
		REQUIRE(destination[0][0][0].value == 7);
		{
			FailCopies failure;
			REQUIRE_THROWS_AS(destination = source, std::bad_alloc);
		}
		REQUIRE(destination.size() == 1);
		REQUIRE(destination[0][0][0].value == 7);
	}
}

TEST_CASE("Dense dimensions and block extents are checked before arithmetic", "[storage-regression]") {
	using namespace MML;
	const int maximum = std::numeric_limits<int>::max();
	Matrix<double> matrix(2, 2, 3.);
	REQUIRE_THROWS_AS(matrix.block(1, 0, maximum, 1), MatrixDimensionError);
	REQUIRE_THROWS_AS(std::as_const(matrix).block(0, 1, 1, maximum), MatrixDimensionError);
	REQUIRE_THROWS_AS(matrix.GetSubmatrix(1, 0, maximum, 1), MatrixDimensionError);
	REQUIRE_THROWS_AS(matrix.GetSubmatrix(0, 1, 1, maximum), MatrixDimensionError);
	REQUIRE(matrix.block(1, 1, 1, 1)(0, 0) == 3.);
	REQUIRE_THROWS_AS(Matrix3D<double>(65536, 65536, 1), MatrixAllocationError);
	Matrix3D<double> volume(1, 1, 1, 7.);
	REQUIRE_THROWS_AS(volume.resize(maximum, maximum, maximum), MatrixAllocationError);
	REQUIRE(volume.dim1() == 1);
	REQUIRE(volume[0][0][0] == 7.);
	REQUIRE_THROWS_AS(volume.Apply({}), NotImplementedError);
	volume.resize(0, 2, 3);
	REQUIRE(volume.empty());
	REQUIRE_THROWS_AS(volume.Apply({}), NotImplementedError);
}

TEST_CASE("Compressed sparse metadata is validated before segment access", "[storage-regression]") {
	using namespace MML;
	using namespace MML::SparseMatrix;
	for (const auto& pointers : {std::vector<int>{0, 2, 1}, {0, -1, 1}, {0, 1, 0}, {1, 1, 1}, {0, 1}}) {
		REQUIRE_THROWS_AS((SparseMatrixCSR<double>(2, 2, {1.}, {0}, pointers)), ArgumentError);
		REQUIRE_THROWS_AS((SparseMatrixCSC<double>(2, 2, {1.}, {0}, pointers)), ArgumentError);
	}
	REQUIRE_THROWS_AS((SparseMatrixCSR<double>(1, 2, {1., 2.}, {1, 0}, {0, 2})), ArgumentError);
	REQUIRE_THROWS_AS((SparseMatrixCSC<double>(2, 1, {1., 2.}, {1, 0}, {0, 2})), ArgumentError);
	REQUIRE(SparseMatrixCSR<double>().Norm1() == 0.);
	REQUIRE(SparseMatrixCSR<double>(2, 0).Norm1() == 0.);
	REQUIRE(SparseMatrixCSR<double>(0, 2).Norm1() == 0.);
	REQUIRE(SparseMatrixCSR<double>().validate());
	REQUIRE(SparseMatrixCSC<double>().validate());
}

TEST_CASE("Sparse vector products support input output aliasing", "[storage-regression][storage-safe-red]") {
	using namespace MML;
	using namespace MML::SparseMatrix;
	SparseMatrixCSR<double> csr(2, 2, {2., 3.}, {1, 0}, {0, 1, 2});
	SparseMatrixCSC<double> csc(csr);
	REQUIRE(csc.validate());
	auto check = [](const auto& matrix) {
		std::vector<double> standard{1., 2.};
		matrix.Multiply(standard, standard);
		CHECK(standard == std::vector<double>{4., 3.});
		standard = {1., 2.};
		matrix.MultiplyTranspose(standard, standard);
		CHECK(standard == std::vector<double>{6., 2.});
		Vector<double> vector{1., 2.};
		matrix.Multiply(vector, vector);
		CHECK(vector[0] == 4.);
		CHECK(vector[1] == 3.);
		vector = Vector<double>{1., 2.};
		matrix.MultiplyTranspose(vector, vector);
		CHECK(vector[0] == 6.);
		CHECK(vector[1] == 2.);
	};
	check(csr);
	check(csc);
	std::vector<double> standard{1., 2.};
	csr.gemv(2., standard, 3., standard);
	CHECK(standard == std::vector<double>{11., 12.});
	Vector<double> vector{1., 2.};
	csr.gemv(2., vector, 3., vector);
	CHECK(vector[0] == 11.);
	CHECK(vector[1] == 12.);

	const SparseMatrixCSR<double> rectangular(2, 3, {2., 3.}, {2, 0}, {0, 1, 2});
	std::vector<double> resized{1., 2., 4.};
	rectangular.Multiply(resized, resized);
	CHECK(resized == std::vector<double>{8., 3.});
	rectangular.MultiplyTranspose(resized, resized);
	CHECK(resized == std::vector<double>{9., 0., 16.});
}

TEST_CASE("Sparse helpers reject invalid dimensions grids and vector sizes", "[storage-regression]") {
	using namespace MML;
	using namespace MML::SparseMatrix;
	const std::vector<double> one{1.}, two{1., 2.};
	std::vector<double> output{7.};
	REQUIRE_THROWS_AS(dot(one, two), VectorDimensionError);
	REQUIRE_THROWS_AS(dot(two, one), VectorDimensionError);
	REQUIRE_THROWS_AS(axpy(1., two, output), VectorDimensionError);
	REQUIRE(output == std::vector<double>{7.});
	REQUIRE_THROWS_AS(residual(IdentityCSR<double>(2), two, one), VectorDimensionError);
	REQUIRE_THROWS_AS(IdentityCSR<double>(-1), MatrixDimensionError);
	REQUIRE_THROWS_AS(Laplacian1D<double>(-1), MatrixDimensionError);
	REQUIRE_THROWS_AS(Laplacian2D<double>(-2, -2), MatrixDimensionError);
	REQUIRE_THROWS_AS(Laplacian3D<double>(2, -2, -2), MatrixDimensionError);
	REQUIRE_THROWS_AS(Laplacian2D<double>(65536, 65536), MatrixAllocationError);
	REQUIRE_THROWS_AS(Laplacian3D<double>(65536, 65536, 65536), MatrixAllocationError);
	REQUIRE_THROWS_AS(Laplacian1D<double>(2, 0.), ArgumentError);
	REQUIRE_THROWS_AS(Laplacian2D<double>(2, 2, 1., -1.), ArgumentError);
	REQUIRE_THROWS_AS(Laplacian3D<double>(2, 2, 2, 1., 1., std::numeric_limits<double>::infinity()), ArgumentError);
	REQUIRE(Laplacian1D<double>(0).rows() == 0);
	const auto stencil = Laplacian2D<double>(2, 3, 2., 1.);
	CHECK(stencil(0, 0) == 2.5);
	CHECK(stencil(0, 1) == -.25);
	CHECK(stencil(0, 2) == -1.);
	SparseMatrixCOO<double> coo(1, 4);
	REQUIRE_THROWS_AS(coo.addLaplacian2D_5pt(0, 0, 0, 65536, 65536, 1.), ArgumentError);
	REQUIRE_THROWS_AS(coo.reserve(-1), ArgumentError);
}

TEST_CASE("High level fixed storage operations validate requested slots", "[storage-regression]") {
	using namespace MML;
	REQUIRE_THROWS_AS((VectorN<double, 3>::UnitVector(-1)), VectorAccessBoundsError);
	REQUIRE_THROWS_AS((VectorN<double, 3>::UnitVector(3)), VectorAccessBoundsError);
	CHECK(VectorN<double, 3>::UnitVector(2)[2] == 1.);
	MatrixNM<double, 2, 3> matrix;
	REQUIRE_THROWS_AS(matrix.VectorFromRow(-1), MatrixAccessBoundsError);
	REQUIRE_THROWS_AS(matrix.VectorFromRow(2), MatrixAccessBoundsError);
	REQUIRE_THROWS_AS(matrix.VectorFromColumn(-1), MatrixAccessBoundsError);
	REQUIRE_THROWS_AS(matrix.VectorFromColumn(3), MatrixAccessBoundsError);
	auto checkVariance = [](auto tensor, int rank) {
		const int before = tensor.numCovar();
		CHECK_THROWS_AS(tensor.setContravar(-1, true), TensorIndexError);
		CHECK_THROWS_AS(tensor.setContravar(rank, true), TensorIndexError);
		CHECK(tensor.numCovar() == before);
		tensor.setContravar(rank - 1, true);
		CHECK(tensor.numContravar() == 1);
		CHECK(tensor.numCovar() == rank - 1);
	};
	checkVariance(Tensor2<2>(2, 0), 2);
	checkVariance(Tensor3<2>(3, 0), 3);
	checkVariance(Tensor4<2>(4, 0), 4);
	checkVariance(Tensor5<2>(5, 0), 5);
}
