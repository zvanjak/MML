#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/base/Vector/Vector.h>
#endif
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <limits>
#include <stdexcept>

TEST_CASE("Vector unsigned sizes preserve int overloads and reject narrowing", "[Vector][unsigned-size]")
{
    using MML::Vector;
    const std::size_t n = 3;
    Vector<double> values(n, 4.0);
    REQUIRE(values.size() == 3);
    values[std::size_t{1}] = 7.0;
    const auto& constant = values;
    CHECK(constant[std::uint64_t{1}] == 7.0);
    CHECK(&values[1] == &values[std::size_t{1}]);
    double& (Vector<double>::*legacyIndex)(int) noexcept = &Vector<double>::operator[];
    CHECK(&(values.*legacyIndex)(1) == &values[std::size_t{1}]);

    double raw[] = {1.0, 2.0, 3.0};
    const Vector<double> copied(n, raw);
    CHECK(copied[2] == 3.0);
    const Vector<double> small(std::uint8_t{3});
    CHECK(small.size() == 3);
    values.resize(std::size_t{5}, true);
    CHECK(values.size() == 5);
    CHECK(values[1] == 7.0);
    values.resize(std::size_t{0});
    CHECK(values.empty());
    CHECK_THROWS_AS(Vector<double>(-1), MML::VectorInitializationError);

    for (const std::uint64_t count : {
        static_cast<std::uint64_t>((std::numeric_limits<int>::max)()) + 1,
        std::uint64_t{1} << 32, (std::numeric_limits<std::uint64_t>::max)()})
    {
        INFO(count);
        CHECK_THROWS_AS(Vector<double>(count), std::length_error);
        CHECK_THROWS_AS(Vector<double>(count, 1.0), std::length_error);
        CHECK_THROWS_AS(Vector<double>(count, raw), std::length_error);
        Vector<double> unchanged{1.0, 2.0};
        CHECK_THROWS_AS(unchanged.resize(count), std::length_error);
        CHECK(unchanged.size() == 2);
        CHECK(unchanged[0] == 1.0);
        CHECK(unchanged[1] == 2.0);
    }
}
