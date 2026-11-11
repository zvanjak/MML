#include "../../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/algorithms/Analyzers/FieldAnalyzers.h>
#endif
#include <catch2/catch_test_macros.hpp>

#include <limits>
#include <stdexcept>
#include <type_traits>

using namespace MML;

namespace {
	class ScalarProbe : public IScalarFunction<3> {
	public:
		mutable int calls = 0;
		Real value = 0;
		bool forbid = false;
		Real operator()(const VectorN<Real, 3>&) const override {
			++calls;
			if (forbid) throw std::logic_error("scalar callback failure");
			return value;
		}
	};

	class VectorProbe : public IVectorFunction<3> {
	public:
		mutable int calls = 0;
		Real value = 0;
		bool forbid = false;
		VectorN<Real, 3> operator()(const VectorN<Real, 3>&) const override {
			++calls;
			if (forbid) throw std::logic_error("vector callback failure");
			return {value, value, value};
		}
	};

	class ScalarProvider : public ScalarProbe, public IScalarDerivativeProvider<3> {
	public:
		Real partial = 0;
		mutable std::vector<Vec3Cart> positions;
		Real EvaluatePartial(const VectorN<Real, 3>& pos, int) const override {
			positions.emplace_back(pos);
			return partial;
		}
		Real EvaluateSecondPartial(const VectorN<Real, 3>& pos, int, int) const override {
			positions.emplace_back(pos);
			return partial;
		}
	};

	class VectorProvider : public VectorProbe, public IVectorDerivativeProvider<3> {
	public:
		Real partial = 0;
		bool overflowingCurl = false;
		mutable int derivativeCalls = 0;
		Real EvaluateComponentPartial(const VectorN<Real, 3>&, int component, int axis) const override {
			++derivativeCalls;
			if (overflowingCurl)
				return component == (axis + 1) % 3 ? std::numeric_limits<Real>::max() : Real(0);
			return partial;
		}
	};

	template<class Operation>
	void EveryRegionMethod(const ScalarFieldAnalyzer& scalar, const VectorFieldAnalyzer& vector,
		const Vec3Cart& lower, const Vec3Cart& upper, int count, Operation operation) {
		operation([&] { return scalar.FindCriticalPoints(lower, upper, count); });
		operation([&] { return scalar.IsHarmonic(lower, upper, count); });
		operation([&] { return vector.IsSolenoidal(lower, upper, count); });
		operation([&] { return vector.IsIrrotational(lower, upper, count); });
		operation([&] { return vector.IsConservative(lower, upper, count); });
		operation([&] { return vector.MaxDivergence(lower, upper, count); });
		operation([&] { return vector.MaxCurlMagnitude(lower, upper, count); });
	}
}

TEST_CASE("Field analyzers reject invalid grid controls before evaluating fields", "[field_analyzer][field-analyzer-contract]") {
	ScalarProbe scalarField;
	VectorProbe vectorField;
	scalarField.forbid = vectorField.forbid = true;
	ScalarFieldAnalyzer scalar(scalarField);
	VectorFieldAnalyzer vector(vectorField);
	const Vec3Cart lower{-1, -1, -1}, upper{1, 1, 1};
	for (int count : {-1, 0}) {
		CAPTURE(count);
		EveryRegionMethod(scalar, vector, lower, upper, count, [](auto operation) {
			CHECK_THROWS_AS(operation(), ArgumentError);
		});
	}
	for (Real bad : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		for (int axis = 0; axis < 3; ++axis) {
			Vec3Cart invalidLower = lower, invalidUpper = upper;
			invalidLower[axis] = invalidUpper[axis] = bad;
			EveryRegionMethod(scalar, vector, invalidLower, upper, 1, [](auto operation) {
				CHECK_THROWS_AS(operation(), NumericInputError);
			});
			EveryRegionMethod(scalar, vector, lower, invalidUpper, 1, [](auto operation) {
				CHECK_THROWS_AS(operation(), NumericInputError);
			});
		}
		EveryRegionMethod(scalar, vector, lower, upper, std::numeric_limits<int>::max(), [](auto operation) {
			CHECK_THROWS_AS(operation(), ResourceLimitError);
		});
	}
	for (int axis = 0; axis < 3; ++axis) {
		Vec3Cart reversed = lower;
		reversed[axis] = upper[axis] + 1;
		EveryRegionMethod(scalar, vector, reversed, upper, 1, [](auto operation) {
			CHECK_THROWS_AS(operation(), ArgumentError);
		});
	}
	CHECK(scalarField.calls == 0);
	CHECK(vectorField.calls == 0);
}

TEST_CASE("Field analyzers reject invalid tolerances before point or region evaluation", "[field_analyzer][field-analyzer-contract]") {
	ScalarProbe scalarField;
	VectorProbe vectorField;
	scalarField.forbid = vectorField.forbid = true;
	ScalarFieldAnalyzer scalar(scalarField);
	VectorFieldAnalyzer vector(vectorField);
	const Vec3Cart lower{-1, -1, -1}, upper{1, 1, 1};
	for (Real tolerance : {Real(-1), Real(0), std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		CAPTURE(tolerance);
		CHECK_THROWS_AS(scalar.IsCriticalPoint(lower, tolerance), NumericInputError);
		CHECK_THROWS_AS(scalar.FindCriticalPoints(lower, upper, 1, tolerance), NumericInputError);
		CHECK_THROWS_AS(scalar.IsHarmonic(lower, upper, 1, tolerance), NumericInputError);
		CHECK_THROWS_AS(vector.IsSolenoidalAt(lower, tolerance), NumericInputError);
		CHECK_THROWS_AS(vector.IsIrrotationalAt(lower, tolerance), NumericInputError);
		CHECK_THROWS_AS(vector.IsSolenoidal(lower, upper, 1, tolerance), NumericInputError);
		CHECK_THROWS_AS(vector.IsIrrotational(lower, upper, 1, tolerance), NumericInputError);
		CHECK_THROWS_AS(vector.IsConservative(lower, upper, 1, tolerance), NumericInputError);
	}
	CHECK(scalarField.calls == 0);
	CHECK(vectorField.calls == 0);
}

TEST_CASE("Field analyzers never turn nonfinite fields into harmonic flags or zero maxima", "[field_analyzer][field-analyzer-contract]") {
	ScalarProbe scalarField;
	VectorProbe vectorField;
	ScalarFieldAnalyzer scalar(scalarField);
	VectorFieldAnalyzer vector(vectorField);
	const Vec3Cart lower{-1, -1, -1}, upper{1, 1, 1};
	for (Real bad : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		scalarField.value = vectorField.value = bad;
		CHECK_THROWS_AS(scalar.GradientAt(lower), NumericalMethodError);
		CHECK_THROWS_AS(scalar.LaplacianAt(lower), NumericalMethodError);
		CHECK_THROWS_AS(scalar.IsCriticalPoint(lower), NumericalMethodError);
		CHECK_THROWS_AS(vector.DivergenceAt(lower), NumericalMethodError);
		CHECK_THROWS_AS(vector.CurlAt(lower), NumericalMethodError);
		CHECK_THROWS_AS(vector.IsSolenoidalAt(lower), NumericalMethodError);
		CHECK_THROWS_AS(vector.IsIrrotationalAt(lower), NumericalMethodError);
		EveryRegionMethod(scalar, vector, lower, upper, 1, [](auto operation) {
			CHECK_THROWS_AS(operation(), NumericalMethodError);
		});
		CHECK_THROWS_AS(vector.CirculationAroundCircle(Vec3Cart{0, 0, 0}, Vec3Cart{0, 0, 1}, 1, 4), NumericalMethodError);
	}
}

TEST_CASE("Field analyzers check provider derivatives and overflowing gradient magnitudes", "[field_analyzer][field-analyzer-contract]") {
	ScalarProvider scalarField;
	VectorProvider vectorField;
	scalarField.forbid = vectorField.forbid = true;
	ScalarFieldAnalyzer scalar(scalarField);
	VectorFieldAnalyzer vector(vectorField);
	const Vec3Cart origin{0, 0, 0};
	for (Real bad : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		scalarField.partial = vectorField.partial = bad;
		CHECK_THROWS_AS(scalar.GradientAt(origin), NumericalMethodError);
		CHECK_THROWS_AS(scalar.LaplacianAt(origin), NumericalMethodError);
		CHECK_THROWS_AS(vector.DivergenceAt(origin), NumericalMethodError);
		CHECK_THROWS_AS(vector.CurlAt(origin), NumericalMethodError);
		EveryRegionMethod(scalar, vector, origin, origin, 1, [](auto operation) {
			CHECK_THROWS_AS(operation(), NumericalMethodError);
		});
	}
	scalarField.partial = std::numeric_limits<Real>::max();
	CHECK_THROWS_AS(scalar.IsCriticalPoint(origin), NumericalMethodError);
	vectorField.overflowingCurl = true;
	CHECK_THROWS_AS(vector.IsIrrotationalAt(origin), NumericalMethodError);
	CHECK_THROWS_AS(vector.IsIrrotational(origin, origin, 1), NumericalMethodError);
	CHECK_THROWS_AS(vector.IsConservative(origin, origin, 1), NumericalMethodError);
	CHECK_THROWS_AS(vector.MaxCurlMagnitude(origin, origin, 1), NumericalMethodError);
	CHECK(scalarField.calls == 0);
	CHECK(vectorField.calls == 0);
}

TEST_CASE("Field analyzer point operations reject nonfinite coordinates before providers", "[field_analyzer][field-analyzer-contract]") {
	ScalarProvider scalarField;
	VectorProvider vectorField;
	ScalarFieldAnalyzer scalar(scalarField);
	VectorFieldAnalyzer vector(vectorField);
	for (Real bad : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		for (int axis = 0; axis < 3; ++axis) {
			Vec3Cart pos{0, 0, 0};
			pos[axis] = bad;
			CHECK_THROWS_AS(scalar.GradientAt(pos), NumericInputError);
			CHECK_THROWS_AS(scalar.LaplacianAt(pos), NumericInputError);
			CHECK_THROWS_AS(scalar.IsCriticalPoint(pos), NumericInputError);
			CHECK_THROWS_AS(vector.DivergenceAt(pos), NumericInputError);
			CHECK_THROWS_AS(vector.CurlAt(pos), NumericInputError);
			CHECK_THROWS_AS(vector.IsSolenoidalAt(pos), NumericInputError);
			CHECK_THROWS_AS(vector.IsIrrotationalAt(pos), NumericInputError);
		}
	}
	CHECK(scalarField.positions.empty());
	CHECK(vectorField.derivativeCalls == 0);
	CHECK(scalarField.calls == 0);
	CHECK(vectorField.calls == 0);
}

TEST_CASE("Field analyzer grids include finite endpoints without overflowing the span", "[field_analyzer][field-analyzer-contract]") {
	ScalarProvider field;
	field.forbid = true;
	ScalarFieldAnalyzer analyzer(field);
	const Real largest = std::numeric_limits<Real>::max();
	const auto points = analyzer.FindCriticalPoints(Vec3Cart{-largest, 0, 0}, Vec3Cart{largest, 0, 0}, 2);
	REQUIRE(points.size() == 27);
	CHECK(points.front().X() == -largest);
	CHECK(points.back().X() == largest);
	CHECK(points[13].X() == 0);
	for (const auto& point : points) {
		CHECK(std::isfinite(point.X()));
		CHECK(point.Y() == 0);
		CHECK(point.Z() == 0);
	}
	CHECK(field.calls == 0);
	VectorProvider vectorField;
	vectorField.forbid = true;
	VectorFieldAnalyzer vector(vectorField);
	CHECK(vector.MaxDivergence(Vec3Cart{0, 0, 0}, Vec3Cart{0, 0, 0}, 2) == 0);
	CHECK(vectorField.derivativeCalls == 27 * 3);
}

TEST_CASE("Field analyzer circulation validates controls and propagates callback failures", "[field_analyzer][field-analyzer-contract]") {
	VectorProbe field;
	VectorFieldAnalyzer analyzer(field);
	const Vec3Cart origin{0, 0, 0}, normal{0, 0, 1};
	for (int count : {-1, 0}) CHECK_THROWS_AS(analyzer.CirculationAroundCircle(origin, normal, 1, count), ArgumentError);
	CHECK_THROWS_AS(analyzer.CirculationAroundCircle(origin, origin, 1, 4), ArgumentError);
	CHECK_THROWS_AS(analyzer.CirculationAroundCircle(origin, normal, -1, 4), ArgumentError);
	for (Real bad : {std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::infinity()}) {
		CHECK_THROWS_AS(analyzer.CirculationAroundCircle(Vec3Cart{bad, 0, 0}, normal, 1, 4), NumericInputError);
		CHECK_THROWS_AS(analyzer.CirculationAroundCircle(origin, Vec3Cart{0, 0, bad}, 1, 4), NumericInputError);
		CHECK_THROWS_AS(analyzer.CirculationAroundCircle(origin, normal, bad, 4), NumericInputError);
	}
	CHECK(field.calls == 0);
	CHECK(analyzer.CirculationAroundCircle(origin, normal, 0, 4) == 0);
	CHECK(analyzer.CirculationAroundCircle(origin, Vec3Cart{0, 0, std::numeric_limits<Real>::max()}, 1, 4) == 0);
	CHECK(analyzer.CirculationAroundCircle(origin, Vec3Cart{0, 0, std::numeric_limits<Real>::min()}, 1, 4) == 0);
	CHECK(analyzer.CirculationAroundCircle(Vec3Cart{std::numeric_limits<Real>::max(), 0, 0}, normal, 0, 4) == 0);
	field.value = std::numeric_limits<Real>::max();
	CHECK_THROWS_AS(analyzer.CirculationAroundCircle(origin, normal, std::numeric_limits<Real>::max() / 4, 4), NumericalMethodError);
	field.forbid = true;
	CHECK_THROWS_AS(analyzer.DivergenceAt(origin), std::logic_error);
	CHECK_THROWS_AS(analyzer.CirculationAroundCircle(origin, normal, 1, 4), std::logic_error);
	ScalarProbe scalarField;
	scalarField.forbid = true;
	ScalarFieldAnalyzer scalar(scalarField);
	CHECK_THROWS_AS(scalar.IsHarmonic(origin, normal, 1), std::logic_error);
}

TEST_CASE("Field analyzers borrow only live lvalue fields", "[field_analyzer][field-analyzer-contract]") {
	CHECK((std::is_constructible_v<ScalarFieldAnalyzer, ScalarProbe&>));
	CHECK((std::is_constructible_v<ScalarFieldAnalyzer, const ScalarProbe&>));
	CHECK_FALSE((std::is_constructible_v<ScalarFieldAnalyzer, ScalarProbe&&>));
	CHECK_FALSE((std::is_constructible_v<ScalarFieldAnalyzer, const ScalarProbe&&>));
	CHECK((std::is_constructible_v<VectorFieldAnalyzer, VectorProbe&>));
	CHECK((std::is_constructible_v<VectorFieldAnalyzer, const VectorProbe&>));
	CHECK_FALSE((std::is_constructible_v<VectorFieldAnalyzer, VectorProbe&&>));
	CHECK_FALSE((std::is_constructible_v<VectorFieldAnalyzer, const VectorProbe&&>));
}
