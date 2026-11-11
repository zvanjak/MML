#include <catch2/catch_test_macros.hpp>
#include <mml/algorithms/MatrixAlg/FieldMatrixAlgorithms.h>

TEST_CASE("Shared field-matrix algorithms remain available without discrete algebra", "[Algebra][MatrixAlg]")
{
	using namespace MML;
	SECTION("Exact pivoting, inverse and unchanged input") {
		const MatrixNM<Real, 2, 2> matrix{{0, 2}, {1, 3}};
		const auto inverse = Algebra::FieldMatrixInverse(matrix);
		REQUIRE(Algebra::FieldMatrixDeterminant(matrix) == Real(-2));
		REQUIRE(inverse == MatrixNM<Real, 2, 2>{{Real(-1.5), 1}, {Real(0.5), 0}});
		REQUIRE(matrix * inverse == MatrixNM<Real, 2, 2>::Identity());
		REQUIRE(matrix == MatrixNM<Real, 2, 2>{{0, 2}, {1, 3}});
	}
	SECTION("A singular field matrix returns zero determinant and throws on inversion") {
		const MatrixNM<Real, 2, 2> singular{{1, 2}, {2, 4}};
		REQUIRE(Algebra::FieldMatrixDeterminant(singular) == Real(0));
		REQUIRE_THROWS_AS(Algebra::FieldMatrixInverse(singular), DomainError);
	}
}
