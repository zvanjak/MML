#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/BaseUtils/MatrixOps.h>
#endif

using namespace MML;
using namespace MML::Testing;
using Catch::Matchers::WithinAbs;

namespace MML::Tests::Base::MatrixOpsTests
{
	TEST_CASE("MatrixOps computes maximum entry difference", "[MatrixOps][Comparison]")
	{
		const Matrix<Real> left{2, 2, {static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0)}};
		const Matrix<Real> right{2, 2, {static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(10.0)}};
		REQUIRE_THAT(Utils::MaxAbsDiff(left, left), WithinAbs(static_cast<MML::Real>(0.0), TOL(1e-14, 1e-5)));
		REQUIRE_THAT(Utils::MaxAbsDiff(left, right), WithinAbs(static_cast<MML::Real>(6.0), TOL(1e-14, 1e-5)));
	}

	TEST_CASE("MatrixOps applies real orthogonal similarity transforms", "[MatrixOps][Transform]")
	{
		const Matrix<Real> matrix{2, 2, {static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0)}};
		const Matrix<Real> identity = Matrix<Real>::Identity(2);
		REQUIRE(Utils::SimilarityTransform(identity, matrix).IsEqualTo(matrix, TOL(1e-12, 1e-5)));

		const Matrix<Real> permutation{2, 2, {static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0)}};
		const Matrix<Real> expected{2, 2, {static_cast<MML::Real>(4.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0)}};
		REQUIRE(Utils::SimilarityTransform(permutation, matrix).IsEqualTo(expected, TOL(1e-12, 1e-5)));
	}

	TEST_CASE("MatrixOps Gram Schmidt helpers produce orthonormal vectors", "[MatrixOps][GramSchmidt]")
	{
		const Matrix<Real> input{3, 2, {
			static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0),
			static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)
		}};
		const Matrix<Real> orthonormal = GramSchmidt(input);
		REQUIRE((orthonormal.Transposed() * orthonormal).IsIdentity(TOL(1e-12, 1e-5)));

		const std::vector<VectorN<Real, 3>> vectors{
			VectorN<Real, 3>{static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0)},
			VectorN<Real, 3>{static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0)}
		};
		const auto basis = GramSchmidtVectors(vectors);
		REQUIRE(basis.size() == 2);
		REQUIRE(IsOrthonormal(basis, TOL(1e-12, 1e-5)));
	}
}
