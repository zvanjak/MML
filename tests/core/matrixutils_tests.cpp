#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/CoreUtils.h"
#endif

#include "../test_data/linear_alg_eq_systems_test_bed.h"

using namespace MML;

// TODO 0.9 - HIGH - finish these tests
TEST_CASE("Test MatrixUtils::IsOrthogonal")
{
    Matrix<Real> A(3, 3, {1,0,0,0,1,0,0,0,1});

    REQUIRE(MatrixUtils::IsOrthogonal(A));

    Matrix<Real> B(3, 3, {1,0,0,0,1,0,0,0,2});

    REQUIRE_FALSE(MatrixUtils::IsOrthogonal(B));

    Matrix<Real> C(3, 3, {1,0,0,0,1,0,0,0,0});

    REQUIRE_FALSE(MatrixUtils::IsOrthogonal(C));

    Matrix<Real> D(3, 3, {1,0,0,0,1,0,0,0,0});

    REQUIRE_FALSE(MatrixUtils::IsOrthogonal(D));

    Matrix<Real> E(3, 3, {1,0,0,0,1,0,0,0,0});

    REQUIRE_FALSE(MatrixUtils::IsOrthogonal(E));

    Matrix<Real> F(3, 3, {1,0,0,0,1,0,0,0,0});

    REQUIRE_FALSE(MatrixUtils::IsOrthogonal(F));

    Matrix<Real> G(3, 3, {1,0,0,0,1,0,0,0,0});

    REQUIRE_FALSE(MatrixUtils::IsOrthogonal(G));

    Matrix<Real> H(3, 3, {1,0,0,0,1,0,0,0,0});

    REQUIRE_FALSE(MatrixUtils::IsOrthogonal(H));

    Matrix<Real> I(3, 3, {1,0,0,0,1,0,0,0,0});

    REQUIRE_FALSE(MatrixUtils::IsOrthogonal(I));

    Matrix<Real> J(3, 3, {1,0,0,0,1,0,0,0,0});

    REQUIRE_FALSE(MatrixUtils::IsOrthogonal(J));

    Matrix<Real> K(3, 3, {1,0,0,0,1,0,0,0,0});

    REQUIRE_FALSE(MatrixUtils::IsOrthogonal(K));

    Matrix<Real> L(3, 3, {1,0,0,0,1,0,0,0,0});

    REQUIRE_FALSE(MatrixUtils::IsOrthogonal(L));
}

TEST_CASE("Test MatrixUtils::Det")
{

}

TEST_CASE("Test Rank on TestBeds matrices")
{
    Matrix<Real> A(3, 3, {1,2,3,4,5,6,7,8,9});

    // interesting
    REQUIRE(2 == MatrixUtils::Rank(A));

    Matrix<Real> B(3, 3, {1,2,-1,4,2,5,-2,3,6});

    REQUIRE(3 == MatrixUtils::Rank(B));

    Matrix<Real> C(3, 3, {1,1,1,2,2,2,-2,3,6});

    REQUIRE(2 == MatrixUtils::Rank(C));

    Matrix<Real> D(3, 3, {1,1,1,2,2,2,3,3,3});

    REQUIRE(1 == MatrixUtils::Rank(D));

    REQUIRE(3 == MatrixUtils::Rank(TestBeds::mat_3x3));
    REQUIRE(3 == MatrixUtils::Rank(TestBeds::mat_3x3_1));
    REQUIRE(3 == MatrixUtils::Rank(TestBeds::mat_3x3_2));
    REQUIRE(3 == MatrixUtils::Rank(TestBeds::mat_3x3_3));
    REQUIRE(3 == MatrixUtils::Rank(TestBeds::mat_3x3_4));
    
    REQUIRE(5 == MatrixUtils::Rank(TestBeds::mat_5x5));
    REQUIRE(5 == MatrixUtils::Rank(TestBeds::mat_5x5_2));
    REQUIRE(5 == MatrixUtils::Rank(TestBeds::mat_5x5_3));
    REQUIRE(5 == MatrixUtils::Rank(TestBeds::mat_5x5_4));

    REQUIRE(8 == MatrixUtils::Rank(TestBeds::mat_8x8));

    REQUIRE(10 == MatrixUtils::Rank(TestBeds::mat_10x10));
    REQUIRE(10 == MatrixUtils::Rank(TestBeds::mat_10x10_1));
    REQUIRE(10 == MatrixUtils::Rank(TestBeds::mat_10x10_2));

    REQUIRE(20 == MatrixUtils::Rank(TestBeds::mat_20x20));
    REQUIRE(20 == MatrixUtils::Rank(TestBeds::mat_20x20_1));
    REQUIRE(20 == MatrixUtils::Rank(TestBeds::mat_20x20_2));

    REQUIRE(50 == MatrixUtils::Rank(TestBeds::mat_50x50));
    REQUIRE(50 == MatrixUtils::Rank(TestBeds::mat_50x50_1));
    REQUIRE(50 == MatrixUtils::Rank(TestBeds::mat_50x50_2));
}

TEST_CASE("Test MatrixUtils::Definitness")
{

}

TEST_CASE("Test MatrixUtils::IsHermitian")
{

}

TEST_CASE("Test MatrixUtils::IsUnitary")
{

}