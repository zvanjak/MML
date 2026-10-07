#include <catch2/catch_all.hpp>

#include <mml/core/Derivation/ForwardAD.h>

#include <cmath>
#include <vector>

using namespace MML::AD;

namespace MML::Tests::Core::Derivation::ForwardADTests
{
    TEST_CASE("ForwardAD dual arithmetic propagates derivatives", "[ForwardAD][Derivation]")
    {
        Dual<double> x(2.0, 1.0);
        Dual<double> y(3.0, 0.0);

        auto sum = x + y;
        REQUIRE(sum.value == Catch::Approx(5.0));
        REQUIRE(sum.deriv == Catch::Approx(1.0));

        auto product = x * y;
        REQUIRE(product.value == Catch::Approx(6.0));
        REQUIRE(product.deriv == Catch::Approx(3.0));

        auto quotient = x / y;
        REQUIRE(quotient.value == Catch::Approx(2.0 / 3.0));
        REQUIRE(quotient.deriv == Catch::Approx(1.0 / 3.0));
    }

    TEST_CASE("ForwardAD math functions propagate analytic derivatives", "[ForwardAD][Derivation]")
    {
        Dual<double> x(0.5, 1.0);

        auto sinValue = sin(x);
        REQUIRE(sinValue.value == Catch::Approx(std::sin(0.5)));
        REQUIRE(sinValue.deriv == Catch::Approx(std::cos(0.5)));

        auto expValue = exp(x);
        REQUIRE(expValue.value == Catch::Approx(std::exp(0.5)));
        REQUIRE(expValue.deriv == Catch::Approx(std::exp(0.5)));

        auto sqrtValue = sqrt(Dual<double>(4.0, 1.0));
        REQUIRE(sqrtValue.value == Catch::Approx(2.0));
        REQUIRE(sqrtValue.deriv == Catch::Approx(0.25));
    }

    TEST_CASE("ForwardAD derivative utility returns value and derivative", "[ForwardAD][Derivation]")
    {
        auto cubic = [](auto x) { return x * x * x - 2.0 * x + 1.0; };

        auto [value, deriv] = derivative(cubic, 2.0);

        REQUIRE(value == Catch::Approx(5.0));
        REQUIRE(deriv == Catch::Approx(10.0));
    }

    TEST_CASE("ForwardAD gradient utility computes scalar gradients", "[ForwardAD][Derivation]")
    {
        auto quadratic = [](const std::vector<Dual<double>>& v) {
            return v[0] * v[0] + v[0] * v[1] + v[1] * v[1];
        };

        auto grad = gradient(quadratic, std::vector<double>{2.0, 3.0});

        REQUIRE(grad.size() == 2);
        REQUIRE(grad[0] == Catch::Approx(7.0));
        REQUIRE(grad[1] == Catch::Approx(8.0));
    }
}