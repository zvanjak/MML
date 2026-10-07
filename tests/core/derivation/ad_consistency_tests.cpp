#include <catch2/catch_all.hpp>

#include <mml/core/Derivation/ADJacobians.h>
#include <mml/core/Derivation/ForwardAD.h>
#include <mml/core/Derivation/ReverseAD.h>

#include <cmath>
#include <vector>

using namespace MML::AD;

namespace MML::Tests::Core::Derivation::ADConsistencyTests
{
    TEST_CASE("Automatic differentiation modes agree on scalar gradient", "[AutomaticDifferentiation][Derivation]")
    {
        const std::vector<Real> point{ 0.5, 1.5 };

        auto forwardFunc = [](const std::vector<Dual<Real>>& v) {
            return sin(v[0] * v[1]) + exp(v[0]);
        };

        auto reverseFunc = [](std::vector<ADVar>& v) {
            return sin(v[0] * v[1]) + exp(v[0]);
        };

        auto forwardGrad = gradient(forwardFunc, point);
        auto reverseGrad = gradientReverse(reverseFunc, point);

        REQUIRE(forwardGrad.size() == reverseGrad.size());
        REQUIRE(forwardGrad[0] == Catch::Approx(reverseGrad[0]).epsilon(std::numeric_limits<Real>::epsilon()));
        REQUIRE(forwardGrad[1] == Catch::Approx(reverseGrad[1]).epsilon(std::numeric_limits<Real>::epsilon()));
    }

    TEST_CASE("Automatic differentiation Jacobian helpers agree across modes", "[AutomaticDifferentiation][Derivation]")
    {
        const std::vector<Real> point{ 3.0, 4.0 };

        auto forwardFunc = [](const std::vector<Dual<Real>>& v) {
            return std::vector<Dual<Real>>{
                v[0] * v[0] + v[1],
                sin(v[0] * v[1])
            };
        };

        auto reverseFunc = [](std::vector<ADVar>& v) {
            return std::vector<ADVar>{
                v[0] * v[0] + v[1],
                sin(v[0] * v[1])
            };
        };

        auto forwardJac = jacobianForwardAD(forwardFunc, point);
        auto reverseJac = jacobianReverseAD(reverseFunc, point);

        REQUIRE(forwardJac.rows() == reverseJac.rows());
        REQUIRE(forwardJac.cols() == reverseJac.cols());

        for (int i = 0; i < forwardJac.rows(); ++i)
            for (int j = 0; j < forwardJac.cols(); ++j)
                REQUIRE(forwardJac(i, j) == Catch::Approx(reverseJac(i, j)).epsilon(std::numeric_limits<Real>::epsilon()));
    }
}