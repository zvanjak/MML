#include <catch2/catch_all.hpp>

#include "../../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/core/Derivation.h>
#endif

#include <initializer_list>
#include <vector>

using namespace MML;

namespace MML::Tests::Core::DerivationCharacterizationTests
{
	class RecordingRealFunction : public IRealFunction
	{
	public:
		mutable std::vector<Real> arguments;

		Real operator()(Real x) const override
		{
			arguments.push_back(x);
			return x * x * x;
		}

		void Clear() const { arguments.clear(); }
	};

	void RequireOffsets(const std::vector<Real>& arguments, Real center, Real step,
	                    std::initializer_list<int> expected_offsets)
	{
		REQUIRE(arguments.size() == expected_offsets.size());

		auto argument = arguments.begin();
		for (int offset : expected_offsets)
		{
			CAPTURE(offset);
			REQUIRE(*argument == Catch::Approx(center + offset * step));
			++argument;
		}
	}

	template<typename Invoke>
	void RequireAdapterOffsets(std::vector<Real>& arguments, Invoke&& invoke,
	                           std::initializer_list<int> value_offsets,
	                           std::initializer_list<int> error_offsets)
	{
		arguments.clear();
		invoke(nullptr);
		RequireOffsets(arguments, REAL(1.25), REAL(0.125), value_offsets);

		Real error = REAL(0.0);
		arguments.clear();
		invoke(&error);
		RequireOffsets(arguments, REAL(1.25), REAL(0.125), error_offsets);
		REQUIRE(std::isfinite(error));
		REQUIRE(error >= REAL(0.0));
	}

	class RecordingComplexFunction : public IComplexFunction
	{
	public:
		mutable std::vector<Real> arguments;

		Complex operator()(Complex z) const override
		{
			arguments.push_back(z.real());
			return z * z * z;
		}
	};

	class RecordingScalarFunction : public IScalarFunction<2>
	{
	public:
		mutable std::vector<Real> arguments;

		Real operator()(const VectorN<Real, 2>& point) const override
		{
			arguments.push_back(point[0]);
			return point[0] * point[0] * point[0] + point[1];
		}
	};

	class RecordingVectorFunction : public IVectorFunction<2>
	{
	public:
		mutable std::vector<Real> arguments;

		VectorN<Real, 2> operator()(const VectorN<Real, 2>& point) const override
		{
			arguments.push_back(point[0]);
			return {point[0] * point[0] * point[0] + point[1], point[0] - point[1]};
		}
	};

	class RecordingCurve : public IParametricCurve<2>
	{
	public:
		mutable std::vector<Real> arguments;

		VectorN<Real, 2> operator()(Real t) const override
		{
			arguments.push_back(t);
			return {t * t * t, REAL(2.0) * t};
		}

		Real getMinT() const override { return REAL(-10.0); }
		Real getMaxT() const override { return REAL(10.0); }
	};

	class RecordingSurface : public IParametricSurfaceRect<2>
	{
	public:
		mutable std::vector<Real> u_arguments;
		mutable std::vector<Real> w_arguments;

		VectorN<Real, 2> operator()(Real u, Real w) const override
		{
			u_arguments.push_back(u);
			w_arguments.push_back(w);
			return {u * u * u + w, w * w * w + u};
		}

		Real getMinU() const override { return REAL(-10.0); }
		Real getMaxU() const override { return REAL(10.0); }
		Real getMinW() const override { return REAL(-10.0); }
		Real getMaxW() const override { return REAL(10.0); }
	};

	class RecordingTensorField2 : public ITensorField2<2>
	{
	public:
		mutable std::vector<Real> arguments;

		RecordingTensorField2() : ITensorField2<2>(0, 2) {}

		Real Component(int, int, const VectorN<Real, 2>& point) const override
		{
			arguments.push_back(point[0]);
			return point[0] * point[0] * point[0] + point[1];
		}

		Tensor2<2> operator()(const VectorN<Real, 2>&) const override { return Tensor2<2>(2, 0); }
	};

	class RecordingTensorField3 : public ITensorField3<2>
	{
	public:
		mutable std::vector<Real> arguments;

		Real Component(int, int, int, const VectorN<Real, 2>& point) const override
		{
			arguments.push_back(point[0]);
			return point[0] * point[0] * point[0] + point[1];
		}

		Tensor3<2> operator()(const VectorN<Real, 2>&) const override { return Tensor3<2>(3, 0); }
	};

	class RecordingTensorField4 : public ITensorField4<2>
	{
	public:
		mutable std::vector<Real> arguments;

		Real Component(int, int, int, int, const VectorN<Real, 2>& point) const override
		{
			arguments.push_back(point[0]);
			return point[0] * point[0] * point[0] + point[1];
		}

		Tensor4<2> operator()(const VectorN<Real, 2>&) const override { return Tensor4<2>(4, 0); }
	};

	TEST_CASE("Finite-difference real stencils preserve samples, counts, values, and errors",
	          "[derivation][characterization][stencil]")
	{
		constexpr Real x = REAL(1.25);
		constexpr Real h = REAL(0.125);
		RecordingRealFunction function;
		auto cubic = [](Real value) { return value * value * value; };

		auto characterize = [&](auto derivative,
		                        std::initializer_list<int> value_offsets,
		                        std::initializer_list<int> error_offsets,
		                        Real expected_value,
		                        Real expected_error)
		{
			DerivativeConfig config;
			config.estimate_error = false;

			function.Clear();
			auto value_only = derivative(function, x, h, config);
			RequireOffsets(function.arguments, x, h, value_offsets);
			REQUIRE(value_only.function_evaluations == static_cast<int>(value_offsets.size()));
			REQUIRE(value_only.value == Catch::Approx(expected_value));
			REQUIRE(value_only.error == REAL(0.0));

			config.estimate_error = true;
			function.Clear();
			auto with_error = derivative(function, x, h, config);
			RequireOffsets(function.arguments, x, h, error_offsets);
			REQUIRE(with_error.function_evaluations == static_cast<int>(error_offsets.size()));
			REQUIRE(with_error.value == Catch::Approx(expected_value));
			REQUIRE(with_error.error == Catch::Approx(expected_error).epsilon(MML::Testing::Tol(REAL(1e-12), REAL(2e-4))));
		};

		SECTION("NDer1")
		{
			characterize(
				[](const IRealFunction& f, Real at, Real step, const DerivativeConfig& config) {
					return Derivation::NDer1Detailed(f, at, step, config);
				},
				{1, 0}, {1, 0, -1}, REAL(5.171875), REAL(0.46875));
		}

		SECTION("NDer2")
		{
			characterize(
				[](const IRealFunction& f, Real at, Real step, const DerivativeConfig& config) {
					return Derivation::NDer2Detailed(f, at, step, config);
				},
				{1, -1}, {1, -1, 2, -2}, REAL(4.703125), REAL(0.015625));
		}

		SECTION("NDer4")
		{
			Real expected_error = Constants::Eps *
				(std::abs(cubic(x + 2 * h)) + std::abs(cubic(x - 2 * h)) +
				 8 * (std::abs(cubic(x - h)) + std::abs(cubic(x + h)))) / (12 * h);
			characterize(
				[](const IRealFunction& f, Real at, Real step, const DerivativeConfig& config) {
					return Derivation::NDer4Detailed(f, at, step, config);
				},
				{1, -1, 2, -2}, {1, -1, 2, -2, 3, -3}, REAL(4.6875), expected_error);
		}

		SECTION("NDer6")
		{
			Real expected_error = 5 *
				(std::abs(cubic(x + h)) + std::abs(cubic(x - h))) * Constants::Eps / h;
			characterize(
				[](const IRealFunction& f, Real at, Real step, const DerivativeConfig& config) {
					return Derivation::NDer6Detailed(f, at, step, config);
				},
				{1, -1, -2, 2, 3, -3}, {1, -1, -2, 2, 3, -3, 4, -4}, REAL(4.6875), expected_error);
		}

		SECTION("NDer8")
		{
			Real expected_error = 7 *
				(std::abs(cubic(x + h)) + std::abs(cubic(x - h))) * Constants::Eps / h;
			characterize(
				[](const IRealFunction& f, Real at, Real step, const DerivativeConfig& config) {
					return Derivation::NDer8Detailed(f, at, step, config);
				},
				{1, -1, -2, 2, 3, -3, -4, 4},
				{1, -1, -2, 2, 3, -3, -4, 4, 5, -5}, REAL(4.6875), expected_error);
		}
	}

	TEST_CASE("Finite-difference real wrappers preserve automatic step scaling",
	          "[derivation][characterization][stencil]")
	{
		constexpr Real x = REAL(2.0);
		RecordingRealFunction function;

		auto characterize = [&](auto invoke_automatic, auto invoke_explicit, Real default_step) {
			Real automatic_error = REAL(0.0);
			function.Clear();
			Real automatic_value = invoke_automatic(&automatic_error);
			auto automatic_arguments = function.arguments;

			Real explicit_error = REAL(0.0);
			function.Clear();
			Real explicit_value = invoke_explicit(Derivation::ScaleStep(default_step, x), &explicit_error);

			REQUIRE(function.arguments.size() == automatic_arguments.size());
			for (std::size_t i = 0; i < automatic_arguments.size(); ++i)
				REQUIRE(function.arguments[i] == Catch::Approx(automatic_arguments[i]));
			REQUIRE(automatic_value == Catch::Approx(explicit_value));
			REQUIRE(automatic_error == Catch::Approx(explicit_error));
		};

		characterize([&](Real* error) { return Derivation::NDer1(function, x, error); },
		             [&](Real h, Real* error) { return Derivation::NDer1(function, x, h, error); }, Derivation::NDer1_h);
		characterize([&](Real* error) { return Derivation::NDer2(function, x, error); },
		             [&](Real h, Real* error) { return Derivation::NDer2(function, x, h, error); }, Derivation::NDer2_h);
		characterize([&](Real* error) { return Derivation::NDer4(function, x, error); },
		             [&](Real h, Real* error) { return Derivation::NDer4(function, x, h, error); }, Derivation::NDer4_h);
		characterize([&](Real* error) { return Derivation::NDer6(function, x, error); },
		             [&](Real h, Real* error) { return Derivation::NDer6(function, x, h, error); }, Derivation::NDer6_h);
		characterize([&](Real* error) { return Derivation::NDer8(function, x, error); },
		             [&](Real h, Real* error) { return Derivation::NDer8(function, x, h, error); }, Derivation::NDer8_h);
	}

	TEST_CASE("Finite-difference real one-sided wrappers preserve sample reach",
	          "[derivation][characterization][stencil]")
	{
		constexpr Real x = REAL(1.25);
		constexpr Real h = REAL(0.125);
		RecordingRealFunction function;

		auto characterize = [&](auto invoke,
		                        std::initializer_list<int> value_offsets,
		                        std::initializer_list<int> error_offsets) {
			RequireAdapterOffsets(function.arguments, invoke, value_offsets, error_offsets);
		};

		characterize([&](Real* error) { return Derivation::NDer1Left(function, x, h, error); },
		             {0, -1}, {0, -1, -2});
		characterize([&](Real* error) { return Derivation::NDer1Right(function, x, h, error); },
		             {1, 0}, {1, 0, -1});
		characterize([&](Real* error) { return Derivation::NDer2Left(function, x, h, error); },
		             {-2, -4}, {-2, -4, -1, -5});
		characterize([&](Real* error) { return Derivation::NDer2Right(function, x, h, error); },
		             {4, 2}, {4, 2, 5, 1});
		characterize([&](Real* error) { return Derivation::NDer4Left(function, x, h, error); },
		             {-3, -5, -2, -6}, {-3, -5, -2, -6, -1, -7});
		characterize([&](Real* error) { return Derivation::NDer4Right(function, x, h, error); },
		             {5, 3, 6, 2}, {5, 3, 6, 2, 7, 1});
		characterize([&](Real* error) { return Derivation::NDer6Left(function, x, h, error); },
		             {-4, -6, -7, -3, -2, -8}, {-4, -6, -7, -3, -2, -8, -1, -9});
		characterize([&](Real* error) { return Derivation::NDer6Right(function, x, h, error); },
		             {6, 4, 3, 7, 8, 2}, {6, 4, 3, 7, 8, 2, 9, 1});
		characterize([&](Real* error) { return Derivation::NDer8Left(function, x, h, error); },
		             {-5, -7, -8, -4, -3, -9, -10, -2}, {-5, -7, -8, -4, -3, -9, -10, -2, -1, -11});
		characterize([&](Real* error) { return Derivation::NDer8Right(function, x, h, error); },
		             {7, 5, 4, 8, 9, 3, 2, 10}, {7, 5, 4, 8, 9, 3, 2, 10, 11, 1});
	}

	TEST_CASE("Finite-difference complex adapters preserve supported stencil reach",
	          "[derivation][characterization][stencil]")
	{
		constexpr Real x = REAL(1.25);
		constexpr Real h = REAL(0.125);
		RecordingComplexFunction function;

		auto characterize = [&](auto invoke, std::initializer_list<int> value_offsets,
		                        std::initializer_list<int> error_offsets)
		{
			DerivativeConfig config;
			config.estimate_error = false;
			function.arguments.clear();
			auto value_only = invoke(config);
			RequireOffsets(function.arguments, x, h, value_offsets);
			REQUIRE(value_only.function_evaluations == static_cast<int>(value_offsets.size()));

			config.estimate_error = true;
			function.arguments.clear();
			auto with_error = invoke(config);
			RequireOffsets(function.arguments, x, h, error_offsets);
			REQUIRE(with_error.function_evaluations == static_cast<int>(error_offsets.size()));
			REQUIRE(with_error.error >= REAL(0.0));
		};

		characterize([&](const DerivativeConfig& config) { return Derivation::NDer1ComplexDetailed(function, Complex(x), h, config); }, {1, 0}, {1, 0, -1});
		characterize([&](const DerivativeConfig& config) { return Derivation::NDer2ComplexDetailed(function, Complex(x), h, config); }, {1, -1}, {1, -1, 2, -2});
		characterize([&](const DerivativeConfig& config) { return Derivation::NDer4ComplexDetailed(function, Complex(x), h, config); }, {1, -1, 2, -2}, {1, -1, 2, -2, 3, -3});
		characterize([&](const DerivativeConfig& config) { return Derivation::NDer6ComplexDetailed(function, Complex(x), h, config); }, {1, -1, -2, 2, 3, -3}, {1, -1, -2, 2, 3, -3, 4, -4});
	}

	TEST_CASE("Finite-difference scalar and vector partial adapters preserve stencil reach",
	          "[derivation][characterization][stencil]")
	{
		const VectorN<Real, 2> point{REAL(1.25), REAL(0.5)};
		constexpr Real h = REAL(0.125);
		RecordingScalarFunction scalar;
		RecordingVectorFunction vector;

		auto characterize_scalar = [&](auto invoke, auto value_offsets, auto error_offsets) {
			RequireAdapterOffsets(scalar.arguments,
				invoke,
				value_offsets, error_offsets);
		};
		auto characterize_vector = [&](auto invoke, auto value_offsets, auto error_offsets) {
			RequireAdapterOffsets(vector.arguments,
				invoke,
				value_offsets, error_offsets);
		};

		#define CHARACTERIZE_PARTIAL(order, value_offsets, error_offsets) \
			characterize_scalar([&](Real* error) { return Derivation::NDer##order##Partial(scalar, 0, point, h, error); }, value_offsets, error_offsets); \
			characterize_vector([&](Real* error) { return Derivation::NDer##order##Partial(vector, 0, 0, point, h, error); }, value_offsets, error_offsets)
		CHARACTERIZE_PARTIAL(1, (std::initializer_list<int>{0, 1}), (std::initializer_list<int>{0, 1, -1}));
		CHARACTERIZE_PARTIAL(2, (std::initializer_list<int>{1, -1}), (std::initializer_list<int>{1, -1, 2, -2}));
		CHARACTERIZE_PARTIAL(4, (std::initializer_list<int>{1, -1, 2, -2}), (std::initializer_list<int>{1, -1, 2, -2, 3, -3}));
		CHARACTERIZE_PARTIAL(6, (std::initializer_list<int>{1, -1, 2, -2, 3, -3}), (std::initializer_list<int>{1, -1, 2, -2, 3, -3, 4, -4}));
		CHARACTERIZE_PARTIAL(8, (std::initializer_list<int>{1, -1, 2, -2, 3, -3, 4, -4}), (std::initializer_list<int>{1, -1, 2, -2, 3, -3, 4, -4, 5, -5}));
		#undef CHARACTERIZE_PARTIAL
	}

	TEST_CASE("Finite-difference curve and surface adapters preserve supported stencil reach",
	          "[derivation][characterization][stencil]")
	{
		constexpr Real x = REAL(1.25);
		constexpr Real h = REAL(0.125);
		RecordingCurve curve;
		RecordingSurface surface;

		auto characterize_curve = [&](auto invoke, std::initializer_list<int> value_offsets,
		                              std::initializer_list<int> error_offsets) {
			RequireAdapterOffsets(curve.arguments,
				invoke,
				value_offsets, error_offsets);
		};

		characterize_curve([&](Real* error) { return Derivation::NDer1(curve, x, h, error); }, {1, 0}, {1, 0, -1});
		characterize_curve([&](Real* error) { return Derivation::NDer2(curve, x, h, error); }, {1, -1}, {1, -1, 2, -2});
		characterize_curve([&](Real* error) { return Derivation::NDer4(curve, x, h, error); }, {1, -1, 2, -2}, {1, -1, 2, -2, 3, -3});
		characterize_curve([&](Real* error) { return Derivation::NDer6(curve, x, h, error); }, {1, -1, -2, 2, 3, -3}, {1, -1, -2, 2, 3, -3, 4, -4});
		characterize_curve([&](Real* error) { return Derivation::NDer8(curve, x, h, error); }, {1, -1, -2, 2, 3, -3, -4, 4}, {1, -1, -2, 2, 3, -3, -4, 4, 5, -5});

		RequireAdapterOffsets(surface.u_arguments,
			[&](Real* error) { surface.w_arguments.clear(); return Derivation::NDer1_u(surface, x, REAL(0.5), h, error); },
			{1, 0}, {1, 0, -1});
		RequireAdapterOffsets(surface.w_arguments,
			[&](Real* error) { surface.u_arguments.clear(); return Derivation::NDer1_w(surface, REAL(0.5), x, h, error); },
			{1, 0}, {1, 0, -1});
		RequireAdapterOffsets(surface.u_arguments,
			[&](Real* error) { surface.w_arguments.clear(); return Derivation::NDer2_u(surface, x, REAL(0.5), h, error); },
			{1, -1}, {1, -1, 2, -2});
		RequireAdapterOffsets(surface.w_arguments,
			[&](Real* error) { surface.u_arguments.clear(); return Derivation::NDer2_w(surface, REAL(0.5), x, h, error); },
			{1, -1}, {1, -1, 2, -2});
	}

	TEST_CASE("Finite-difference tensor adapters preserve rank and supported stencil reach",
	          "[derivation][characterization][stencil]")
	{
		const VectorN<Real, 2> point{REAL(1.25), REAL(0.5)};
		constexpr Real h = REAL(0.125);
		RecordingTensorField2 rank2;
		RecordingTensorField3 rank3;
		RecordingTensorField4 rank4;

		auto characterize = [&](auto& field, auto invoke, auto value_offsets, auto error_offsets) {
			RequireAdapterOffsets(field.arguments, invoke, value_offsets, error_offsets);
		};
		auto characterize_order = [&](auto rank2_invoke, auto rank3_invoke, auto rank4_invoke,
		                              std::initializer_list<int> value_offsets,
		                              std::initializer_list<int> error_offsets) {
			characterize(rank2, rank2_invoke, value_offsets, error_offsets);
			characterize(rank3, rank3_invoke, value_offsets, error_offsets);
			characterize(rank4, rank4_invoke, value_offsets, error_offsets);
		};

		characterize_order(
			[&](Real* error) { return Derivation::NDer1Partial(rank2, 0, 0, 0, point, h, error); },
			[&](Real* error) { return Derivation::NDer1Partial(rank3, 0, 0, 0, 0, point, h, error); },
			[&](Real* error) { return Derivation::NDer1Partial(rank4, 0, 0, 0, 0, 0, point, h, error); },
			{0, 1}, {0, 1, -1});
		characterize_order(
			[&](Real* error) { return Derivation::NDer2Partial(rank2, 0, 0, 0, point, h, error); },
			[&](Real* error) { return Derivation::NDer2Partial(rank3, 0, 0, 0, 0, point, h, error); },
			[&](Real* error) { return Derivation::NDer2Partial(rank4, 0, 0, 0, 0, 0, point, h, error); },
			{1, -1}, {1, -1, 2, -2});
		characterize_order(
			[&](Real* error) { return Derivation::NDer4Partial(rank2, 0, 0, 0, point, h, error); },
			[&](Real* error) { return Derivation::NDer4Partial(rank3, 0, 0, 0, 0, point, h, error); },
			[&](Real* error) { return Derivation::NDer4Partial(rank4, 0, 0, 0, 0, 0, point, h, error); },
			{1, -1, 2, -2}, {1, -1, 2, -2, 3, -3});
	}
}