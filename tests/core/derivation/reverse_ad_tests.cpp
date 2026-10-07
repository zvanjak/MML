#include <catch2/catch_all.hpp>

#include <mml/core/Derivation/ReverseAD.h>

#include <cmath>
#include <vector>

using namespace MML::AD;

namespace MML::Tests::Core::Derivation::ReverseADTests
{
    TEST_CASE("ReverseAD arithmetic propagates adjoints", "[ReverseAD][Derivation]")
    {
        Tape tape;
        ADVar x(&tape, 2.0);
        ADVar y(&tape, 3.0);

        ADVar z = x * y + x / y;
        z.backward();

        REQUIRE(z.value() == Catch::Approx(6.0 + 2.0 / 3.0));
        REQUIRE(x.adjoint() == Catch::Approx(3.0 + 1.0 / 3.0));
        REQUIRE(y.adjoint() == Catch::Approx(2.0 - 2.0 / 9.0));
    }

    TEST_CASE("ReverseAD math functions propagate analytic adjoints", "[ReverseAD][Derivation]")
    {
        Tape tape;
        ADVar x(&tape, 0.5);

        ADVar y = sin(x) + exp(x) + log(x);
        y.backward();

        REQUIRE(y.value() == Catch::Approx(std::sin(0.5) + std::exp(0.5) + std::log(0.5)));
        REQUIRE(x.adjoint() == Catch::Approx(std::cos(0.5) + std::exp(0.5) + 2.0));
    }

    TEST_CASE("ReverseAD scalar constants are tracked when combined with variables", "[ReverseAD][Derivation]")
    {
        Tape tape;
        ADVar x(&tape, 4.0);

        ADVar y = 2.0 * x + ADVar(3.0);
        y.backward();

        REQUIRE(y.value() == Catch::Approx(11.0));
        REQUIRE(x.adjoint() == Catch::Approx(2.0));
    }

    TEST_CASE("ReverseAD gradientReverse computes scalar gradients", "[ReverseAD][Derivation]")
    {
        auto quadratic = [](std::vector<ADVar>& v) {
            return v[0] * v[0] + v[0] * v[1] + v[1] * v[1];
        };

        auto grad = gradientReverse(quadratic, std::vector<Real>{ 2.0, 3.0 });

        REQUIRE(grad.size() == 2);
        REQUIRE(grad[0] == Catch::Approx(7.0));
        REQUIRE(grad[1] == Catch::Approx(8.0));
    }

    TEST_CASE("ReverseAD rejects unsafe tape and domain operations", "[ReverseAD][Derivation]")
    {
        SECTION("backward on untracked constant")
        {
            ADVar c(2.0);
            REQUIRE_THROWS_AS(c.backward(), MML::ArgumentError);
        }

        SECTION("mixed tapes")
        {
            Tape first;
            Tape second;
            ADVar x(&first, 1.0);
            ADVar y(&second, 2.0);
            REQUIRE_THROWS_AS(x + y, MML::ArgumentError);
        }

        SECTION("division by zero")
        {
            Tape tape;
            ADVar x(&tape, 1.0);
            REQUIRE_THROWS_AS(x / 0.0, MML::DivisionByZeroError);
        }

        SECTION("log and sqrt domains")
        {
            REQUIRE_THROWS_AS(log(ADVar(0.0)), MML::DomainError);
            REQUIRE_THROWS_AS(sqrt(ADVar(-1.0)), MML::DomainError);
        }
    }
}