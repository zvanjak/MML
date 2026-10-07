#include <catch2/catch_all.hpp>

#include <mml/core/Derivation/ADJacobians.h>

#include <vector>

using namespace MML;
using namespace MML::AD;

namespace MML::Tests::Core::Derivation::ADJacobiansTests
{
    TEST_CASE("ADJacobians forward mode computes dynamic Jacobian", "[ADJacobians][Derivation]")
    {
        auto func = [](const std::vector<Dual<double>>& v) {
            return std::vector<Dual<double>>{
                v[0] * v[0] + v[1],
                v[0] * v[1]
            };
        };

        auto jac = jacobianForwardAD(func, std::vector<double>{ 3.0, 4.0 });

        REQUIRE(jac.rows() == 2);
        REQUIRE(jac.cols() == 2);
        REQUIRE(jac(0, 0) == Catch::Approx(6.0));
        REQUIRE(jac(0, 1) == Catch::Approx(1.0));
        REQUIRE(jac(1, 0) == Catch::Approx(4.0));
        REQUIRE(jac(1, 1) == Catch::Approx(3.0));
    }

    TEST_CASE("ADJacobians reverse mode computes dynamic Jacobian", "[ADJacobians][Derivation]")
    {
        auto func = [](std::vector<ADVar>& v) {
            return std::vector<ADVar>{
                v[0] * v[0] + v[1],
                v[0] * v[1]
            };
        };

        Vector<Real> point{ 3.0, 4.0 };
        auto jac = jacobianReverseAD(func, point);

        REQUIRE(jac.rows() == 2);
        REQUIRE(jac.cols() == 2);
        REQUIRE(jac(0, 0) == Catch::Approx(6.0));
        REQUIRE(jac(0, 1) == Catch::Approx(1.0));
        REQUIRE(jac(1, 0) == Catch::Approx(4.0));
        REQUIRE(jac(1, 1) == Catch::Approx(3.0));
    }

    TEST_CASE("ADJacobians forward-over-forward computes Hessian", "[ADJacobians][Derivation]")
    {
        using Dual2 = Dual<Dual<double>>;
        auto func = [](const std::vector<Dual2>& v) {
            return v[0] * v[0] + Dual2(2.0) * v[0] * v[1] + Dual2(3.0) * v[1] * v[1];
        };

        auto hessian = hessianForwardAD(func, std::vector<double>{ 1.0, 1.0 });

        REQUIRE(hessian.rows() == 2);
        REQUIRE(hessian.cols() == 2);
        REQUIRE(hessian(0, 0) == Catch::Approx(2.0));
        REQUIRE(hessian(0, 1) == Catch::Approx(2.0));
        REQUIRE(hessian(1, 0) == Catch::Approx(2.0));
        REQUIRE(hessian(1, 1) == Catch::Approx(6.0));
    }
}