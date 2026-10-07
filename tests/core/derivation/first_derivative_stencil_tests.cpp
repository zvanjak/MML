#include <catch2/catch_all.hpp>

#include <mml/core/Derivation/FirstDerivativeStencil.h>

#include <vector>

using namespace MML;

namespace MML::Tests::Core::FirstDerivativeStencilTests
{
	using Derivation::Detail::EvaluateFirstDerivativeStencil;
	using Derivation::Detail::FirstDerivativeOrder;
	using Derivation::Detail::FirstDerivativeStencil;

	TEST_CASE("First derivative stencil metadata describes supported formulas", "[derivation][stencil-engine]")
	{
		STATIC_REQUIRE(FirstDerivativeStencil<FirstDerivativeOrder::One>::value_coefficients ==
		               std::array<Real, 2>{REAL(1.0), REAL(-1.0)});
		STATIC_REQUIRE(FirstDerivativeStencil<FirstDerivativeOrder::Two>::value_coefficients ==
		               std::array<Real, 2>{REAL(0.5), REAL(-0.5)});
		STATIC_REQUIRE(FirstDerivativeStencil<FirstDerivativeOrder::Four>::value_coefficients ==
		               std::array<Real, 4>{REAL(2.0 / 3.0), REAL(-2.0 / 3.0), REAL(-1.0 / 12.0), REAL(1.0 / 12.0)});
		STATIC_REQUIRE(FirstDerivativeStencil<FirstDerivativeOrder::Six>::value_coefficients ==
		               std::array<Real, 6>{REAL(3.0 / 4.0), REAL(-3.0 / 4.0), REAL(3.0 / 20.0), REAL(-3.0 / 20.0), REAL(1.0 / 60.0), REAL(-1.0 / 60.0)});
		STATIC_REQUIRE(FirstDerivativeStencil<FirstDerivativeOrder::Eight>::value_coefficients ==
		               std::array<Real, 8>{REAL(4.0 / 5.0), REAL(-4.0 / 5.0), REAL(1.0 / 5.0), REAL(-1.0 / 5.0), REAL(4.0 / 105.0), REAL(-4.0 / 105.0), REAL(1.0 / 280.0), REAL(-1.0 / 280.0)});
		STATIC_REQUIRE(FirstDerivativeStencil<FirstDerivativeOrder::One>::truncation_error_coefficients[1] == REAL(-1.0));
		STATIC_REQUIRE(FirstDerivativeStencil<FirstDerivativeOrder::Two>::truncation_error_coefficients[2] == REAL(1.0 / 12.0));
		STATIC_REQUIRE(FirstDerivativeStencil<FirstDerivativeOrder::Four>::truncation_error_coefficients[4] == REAL(1.0 / 60.0));
		STATIC_REQUIRE(FirstDerivativeStencil<FirstDerivativeOrder::Six>::truncation_error_coefficients[6] == REAL(1.0 / 280.0));
		STATIC_REQUIRE(FirstDerivativeStencil<FirstDerivativeOrder::Eight>::truncation_error_coefficients[8] == REAL(1.0 / 1260.0));
		STATIC_REQUIRE(FirstDerivativeStencil<FirstDerivativeOrder::Four>::max_reach == 3);
		STATIC_REQUIRE(FirstDerivativeStencil<FirstDerivativeOrder::Eight>::max_reach == 5);
	}

	TEST_CASE("First derivative stencil caches samples and preserves evaluation order", "[derivation][stencil-engine]")
	{
		constexpr Real h = REAL(0.125);
		auto characterize = [&]<FirstDerivativeOrder Order>(std::initializer_list<int> value_offsets,
		                                                    std::initializer_list<int> error_offsets) {
			std::vector<int> offsets;
			auto sample = [&](int offset) {
				offsets.push_back(offset);
				Real x = REAL(1.25) + h * offset;
				return x * x * x;
			};

			auto value_only = EvaluateFirstDerivativeStencil<Order>(sample,
				[](Real value) { return std::abs(value); }, h, false);
			REQUIRE(offsets == std::vector<int>(value_offsets));
			REQUIRE(value_only.function_evaluations == static_cast<int>(value_offsets.size()));
			REQUIRE(value_only.error == REAL(0.0));

			offsets.clear();
			auto with_error = EvaluateFirstDerivativeStencil<Order>(sample,
				[](Real value) { return std::abs(value); }, h, true);
			REQUIRE(offsets == std::vector<int>(error_offsets));
			REQUIRE(with_error.function_evaluations == static_cast<int>(error_offsets.size()));
			REQUIRE(with_error.error >= REAL(0.0));
		};

		characterize.template operator()<FirstDerivativeOrder::One>({1, 0}, {1, 0, -1});
		characterize.template operator()<FirstDerivativeOrder::Two>({1, -1}, {1, -1, 2, -2});
		characterize.template operator()<FirstDerivativeOrder::Four>({1, -1, 2, -2}, {1, -1, 2, -2, 3, -3});
		characterize.template operator()<FirstDerivativeOrder::Six>({1, -1, -2, 2, 3, -3}, {1, -1, -2, 2, 3, -3, 4, -4});
		characterize.template operator()<FirstDerivativeOrder::Eight>({1, -1, -2, 2, 3, -3, -4, 4}, {1, -1, -2, 2, 3, -3, -4, 4, 5, -5});
	}

	TEST_CASE("First derivative stencil supports scalar complex and vector outputs", "[derivation][stencil-engine]")
	{
		constexpr Real h = REAL(0.125);
		auto scalar = EvaluateFirstDerivativeStencil<FirstDerivativeOrder::Four>(
			[&](int offset) { Real x = REAL(1.25) + h * offset; return x * x * x; },
			[](Real value) { return std::abs(value); }, h, false);
		auto complex = EvaluateFirstDerivativeStencil<FirstDerivativeOrder::Four>(
			[&](int offset) { Complex x(REAL(1.25) + h * offset, REAL(0.0)); return x * x * x; },
			[](const Complex& value) { return std::abs(value); }, h, true);
		auto vector = EvaluateFirstDerivativeStencil<FirstDerivativeOrder::Four>(
			[&](int offset) { Real x = REAL(1.25) + h * offset; return VectorN<Real, 2>{x * x * x, 2 * x}; },
			[](const VectorN<Real, 2>& value) { return value.NormL2(); }, h, true);

		REQUIRE(scalar.value == Catch::Approx(REAL(4.6875)));
		REQUIRE(scalar.function_evaluations == 4);
		REQUIRE(scalar.error == REAL(0.0));
		REQUIRE(complex.value.real() == Catch::Approx(REAL(4.6875)));
		REQUIRE(complex.function_evaluations == 6);
		REQUIRE(vector.value[0] == Catch::Approx(REAL(4.6875)));
		REQUIRE(vector.value[1] == Catch::Approx(REAL(2.0)));
		REQUIRE(vector.function_evaluations == 6);
	}
}