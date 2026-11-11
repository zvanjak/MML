#include <catch2/catch_all.hpp>

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/base/Function.h>
#include <mml/tools/Serializer.h>
#include <mml/tools/TextIOFormat.h>
#endif

#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace MML;

namespace {
	struct SamplingFile {
		std::filesystem::path path = std::filesystem::temp_directory_path() /
			("mml_sampling_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".mml");
		~SamplingFile() {
			std::error_code error;
			std::filesystem::remove(path, error);
			if (error) std::cerr << "Cannot remove sampling test file: " << error.message() << '\n';
		}
		std::string Read() const {
			std::ifstream in(path);
			REQUIRE(in.is_open());
			std::ostringstream out;
			out << in.rdbuf();
			REQUIRE_FALSE(in.bad());
			return out.str();
		}
	};

	using Rows = std::vector<std::vector<Real>>;

	Rows ReadRows(const std::string& text, int headerLines, int columns) {
		std::istringstream input(text);
		std::string line;
		for (int i = 0; i < headerLines; ++i)
			REQUIRE(bool(std::getline(input, line)));
		Rows rows;
		while (std::getline(input, line)) {
			std::istringstream row(line);
			std::vector<Real> values;
			for (int i = 0; i < columns; ++i) {
				std::string token;
				REQUIRE(bool(row >> token));
				INFO(token);
				Real value = 0;
				REQUIRE(TextIO::ParseReal(token, value));
				CHECK(std::isfinite(value));
				values.push_back(value);
			}
			std::string extra;
			CHECK_FALSE(bool(row >> extra));
			rows.push_back(values);
		}
		return rows;
	}

	struct Axis {
		Real first, last;
		std::vector<Real> expected;
		Real interiorTolerance = 0;
	};

	std::vector<Axis> Axes() {
		const Real big = std::ldexp(Real(1), std::numeric_limits<Real>::max_exponent - 1);
		const Real largest = std::numeric_limits<Real>::max();
		const Real tiny = std::numeric_limits<Real>::denorm_min();
		return {
			{-big, big, {-big, 0, big}},
			{0, big, {0, big / 2, big}},
			{-big, 0, {-big, -big / 2, 0}},
			{-Real(0), big, {-Real(0), big / 2, big}},
			{-big, -Real(0), {-big, -big / 2, -Real(0)}},
			{-largest, largest, {-largest, 0, largest}},
			{0, tiny, {0, 0, tiny}},
			{-tiny, tiny, {-tiny, 0, tiny}},
			{4 * tiny, 8 * tiny, {4 * tiny, 6 * tiny, 8 * tiny}},
			{Real(0.1), Real(0.3), {Real(0.1), (Real(0.1) + Real(0.3)) / Real(2), Real(0.3)},
				std::numeric_limits<Real>::epsilon() * Real(0.3)}
		};
	}

	bool MatchesCoordinate(Real actual, const Axis& axis, std::size_t index) {
		if (index == 0 || index + 1 == axis.expected.size())
			return actual == axis.expected[index];
		return std::abs(actual - axis.expected[index]) <= axis.interiorTolerance;
	}

	void CheckCoordinates(const std::vector<Real>& actual, const Axis& axis) {
		REQUIRE(actual.size() == axis.expected.size());
		for (std::size_t i = 0; i < actual.size(); ++i) {
			CAPTURE(i, axis.first, axis.last);
			CHECK(std::isfinite(actual[i]));
			CHECK(MatchesCoordinate(actual[i], axis, i));
		}
		CHECK(actual.front() == axis.first);
		CHECK(actual.back() == axis.last);
		CHECK(std::signbit(actual.front()) == std::signbit(axis.first));
		CHECK(std::signbit(actual.back()) == std::signbit(axis.last));
	}

	void CheckColumn(const Rows& rows, std::size_t column, const Axis& axis) {
		std::vector<Real> coordinates;
		for (const auto& row : rows) coordinates.push_back(row[column]);
		CheckCoordinates(coordinates, axis);
	}

	template<int N>
	struct SamplingCurve : IRealToVectorFunction<N> {
		mutable std::vector<Real> parameters;
		VectorN<Real, N> operator()(Real t) const override {
			parameters.push_back(t);
			return VectorN<Real, N>{};
		}
	};

	struct SamplingSurface : IParametricSurface<3>, IParametricSurfaceRect<3> {
		Real first, last;
		mutable std::vector<VectorN<Real, 2>> parameters;
		SamplingSurface(Real a, Real b) : first(a), last(b) {}
		Real getMinU() const override { return first; }
		Real getMaxU() const override { return last; }
		Real getMinW() const override { return first; }
		Real getMaxW() const override { return last; }
		Real getMinW(Real) const override { return first; }
		Real getMaxW(Real) const override { return last; }
		VectorN<Real, 3> operator()(Real u, Real w) const override {
			parameters.push_back({u, w});
			return {0, 0, 0};
		}
		VectorN<Real, 3> operator()(const VectorN<Real, 2>& p) const override {
			return (*this)(p[0], p[1]);
		}
	};

	template<int N>
	void CheckGrid(const std::vector<VectorN<Real, N>>& parameters, const Axis& axis) {
		std::size_t count = 1;
		for (int d = 0; d < N; ++d) count *= axis.expected.size();
		REQUIRE(parameters.size() == count);
		for (std::size_t i = 0; i < count; ++i) {
			std::size_t index = i;
			for (int d = N - 1; d >= 0; --d) {
				CHECK(MatchesCoordinate(parameters[i][d], axis, index % axis.expected.size()));
				index /= axis.expected.size();
			}
		}
	}

	template<int N>
	void CheckGridRows(const Rows& rows, const Axis& axis) {
		std::vector<VectorN<Real, N>> coordinates;
		for (const auto& row : rows) {
			VectorN<Real, N> point;
			for (int d = 0; d < N; ++d) point[d] = row[d];
			coordinates.push_back(point);
		}
		CheckGrid(coordinates, axis);
	}
}

TEST_CASE("Real function serializers sample finite represented intervals", "[serializer][serializer-sampling]")
{
	for (const auto& axis : Axes()) {
		for (bool valuesOnly : {false, true}) {
			std::vector<Real> parameters;
			RealFunctionFromStdFunc function([&](Real x) { parameters.push_back(x); return Real(1); });
			std::ostringstream out;
			const int precision = std::numeric_limits<Real>::max_digits10;
			const auto result = valuesOnly
				? Serializer::SaveRealFuncEquallySpaced(out, function, "constant", axis.first, axis.last, 3, precision)
				: Serializer::SaveRealFunc(out, function, "constant", axis.first, axis.last, 3, precision);
			REQUIRE(result.success);
			CheckCoordinates(parameters, axis);
			const auto rows = ReadRows(out.str(), 6, valuesOnly ? 1 : 2);
			REQUIRE(rows.size() == 3);
			if (!valuesOnly) CheckColumn(rows, 0, axis);
			for (const auto& row : rows) CHECK(row.back() == 1);
			parameters.clear();
			SamplingFile file;
			const auto saved = valuesOnly
				? Serializer::SaveRealFuncEquallySpaced(function, "constant", axis.first, axis.last, 3, file.path.string(), precision)
				: Serializer::SaveRealFunc(function, "constant", axis.first, axis.last, 3, file.path.string(), precision);
			REQUIRE(saved.success);
			CHECK(file.Read() == out.str());
			CheckCoordinates(parameters, axis);
		}
	}
}

TEST_CASE("Multi function overloads share finite sample coordinates", "[serializer][serializer-sampling]")
{
	for (const auto& axis : Axes()) {
		std::vector<Real> parameters;
		RealFunctionFromStdFunc function([&](Real x) { parameters.push_back(x); return Real(1); });
		const std::vector<IRealFunction*> functions{&function};
		std::ostringstream out;
		REQUIRE(Serializer::SaveRealMultiFunc(out, functions, "constant", {"one"}, axis.first, axis.last, 3,
			std::numeric_limits<Real>::max_digits10).success);
		CheckCoordinates(parameters, axis);
		CheckColumn(ReadRows(out.str(), 8, 2), 0, axis);
		parameters.clear();
		SamplingFile file;
		REQUIRE(Serializer::SaveRealMultiFunc(functions, "constant", {"one"}, axis.first, axis.last, 3,
			file.path.string(), std::numeric_limits<Real>::max_digits10).success);
		CHECK(file.Read() == out.str());
		CheckCoordinates(parameters, axis);
		parameters.clear();
		const std::vector<RealFunctionFromStdFunc> values{function};
		REQUIRE(Serializer::SaveRealMultiFunc(values, "constant", {"one"}, axis.first, axis.last, 3, file.path.string()).success);
		CheckColumn(ReadRows(file.Read(), 8, 2), 0, axis);
		CheckCoordinates(parameters, axis);
	}
}

TEST_CASE("Sample grids preserve exact endpoints and dyadic unit changes", "[serializer][serializer-sampling]")
{
	const int minimumExponent = std::numeric_limits<Real>::min_exponent - std::numeric_limits<Real>::digits;
	const int maximumExponent = std::numeric_limits<Real>::max_exponent;
	for (int exponent : {minimumExponent + 3, minimumExponent * 2 / 3, 0,
		maximumExponent * 2 / 3, maximumExponent - 3}) {
		const Real unit = std::ldexp(Real(1), exponent);
		for (int count : {3, 5, 9}) {
			for (int offset : {-4, -2}) {
				Axis axis{offset * unit, (offset + 8) * unit, {}};
				for (int i = 0; i < count; ++i)
					axis.expected.push_back((offset + i * (8 / (count - 1))) * unit);
				std::vector<Real> parameters;
				RealFunctionFromStdFunc function([&](Real x) { parameters.push_back(x); return Real(0); });
				std::ostringstream out;
				REQUIRE(Serializer::SaveRealFunc(out, function, "lattice", axis.first, axis.last, count,
					std::numeric_limits<Real>::max_digits10).success);
				CheckCoordinates(parameters, axis);
				CheckColumn(ReadRows(out.str(), 6, 2), 0, axis);
			}
		}
	}
	for (const auto& axis : Axes()) {
		for (int count : {2, 4, 8, 97}) {
			std::vector<Real> parameters;
			RealFunctionFromStdFunc function([&](Real x) { parameters.push_back(x); return Real(0); });
			std::ostringstream out;
			REQUIRE(Serializer::SaveRealFunc(out, function, "rounded", axis.first, axis.last, count).success);
			REQUIRE(parameters.size() == static_cast<std::size_t>(count));
			CHECK(parameters.front() == axis.first);
			CHECK(parameters.back() == axis.last);
			CHECK(std::signbit(parameters.front()) == std::signbit(axis.first));
			CHECK(std::signbit(parameters.back()) == std::signbit(axis.last));
			for (std::size_t i = 1; i < parameters.size(); ++i) {
				CHECK(std::isfinite(parameters[i]));
				CHECK(parameters[i] >= parameters[i - 1]);
				CHECK(parameters[i] <= axis.last);
			}
		}
	}
}

TEST_CASE("Curve exporters preserve finite generated parameter coordinates", "[serializer][serializer-sampling]")
{
	for (const auto& axis : Axes()) {
		SamplingFile file;
		SamplingCurve<2> curve2;
		REQUIRE(Serializer::SaveParamCurveCartesian2DResult(curve2, "curve", axis.first, axis.last, 3, file.path.string()).success);
		CheckCoordinates(curve2.parameters, axis);
		CheckColumn(ReadRows(file.Read(), 6, 3), 0, axis);
		SamplingCurve<3> curve3;
		REQUIRE(Serializer::SaveParamCurveCartesian3DResult(curve3, "curve", axis.first, axis.last, 3, file.path.string()).success);
		CheckCoordinates(curve3.parameters, axis);
		CheckColumn(ReadRows(file.Read(), 6, 4), 0, axis);
		REQUIRE(Serializer::SaveAsParamCurve<2>(std::vector<VectorN<Real, 2>>(3),
			SerializeFormatType::PARAMETRIC_CURVE_CARTESIAN_2D, "stored", axis.first, axis.last, 3, file.path.string()).success);
		CheckColumn(ReadRows(file.Read(), 6, 3), 0, axis);
		REQUIRE(Serializer::SaveAsParamCurve<3>(std::vector<VectorN<Real, 3>>(3),
			SerializeFormatType::PARAMETRIC_CURVE_CARTESIAN_3D, "stored", axis.first, axis.last, 3, file.path.string()).success);
		CheckColumn(ReadRows(file.Read(), 6, 4), 0, axis);
		REQUIRE(Serializer::SaveAsParamCurve2D(Vector<Real>{1, 2, 3}, Vector<Real>{4, 5, 6}, "components",
			file.path.string(), axis.first, axis.last).success);
		const auto rows = ReadRows(file.Read(), 6, 3);
		CheckColumn(rows, 0, axis);
		for (std::size_t i = 0; i < rows.size(); ++i) {
			CHECK(rows[i][1] == Real(i + 1));
			CHECK(rows[i][2] == Real(i + 4));
		}
	}
}

TEST_CASE("Surface exporter overloads sample both parameter axes safely", "[serializer][serializer-sampling]")
{
	for (const auto& axis : Axes()) {
		for (int overload = 0; overload < 3; ++overload) {
			SamplingFile file;
			SamplingSurface surface(axis.first, axis.last);
			IOResult result;
			if (overload == 0)
				result = Serializer::SaveParametricSurface(static_cast<const IParametricSurface<3>&>(surface),
					"surface", axis.first, axis.last, 3, axis.first, axis.last, 3, file.path.string());
			else if (overload == 1)
				result = Serializer::SaveParametricSurface(
					static_cast<const IVectorFunctionNM<2, 3>&>(static_cast<const IParametricSurfaceRect<3>&>(surface)),
					"surface", axis.first, axis.last, 3, axis.first, axis.last, 3, file.path.string());
			else
				result = Serializer::SaveParametricSurface(static_cast<const IParametricSurfaceRect<3>&>(surface),
					"surface", 3, 3, file.path.string());
			REQUIRE(result.success);
			CheckGrid(surface.parameters, axis);
			CheckGridRows<2>(ReadRows(file.Read(), 9, 5), axis);
		}
	}
}

TEST_CASE("Scalar field exporters sample every grid axis safely", "[serializer][serializer-sampling]")
{
	for (const auto& axis : Axes()) {
		SamplingFile file;
		std::vector<VectorN<Real, 2>> parameters2;
		ScalarFunctionFromStdFunc<2> scalar2([&](const VectorN<Real, 2>& p) { parameters2.push_back(p); return Real(1); });
		REQUIRE(Serializer::SaveScalarFunc2DCartesian(scalar2, "scalar", axis.first, axis.last, 3,
			axis.first, axis.last, 3, file.path.string()).success);
		CheckGrid(parameters2, axis);
		CheckGridRows<2>(ReadRows(file.Read(), 9, 3), axis);
		std::vector<VectorN<Real, 3>> parameters3;
		ScalarFunctionFromStdFunc<3> scalar3([&](const VectorN<Real, 3>& p) { parameters3.push_back(p); return Real(1); });
		REQUIRE(Serializer::SaveScalarFunc3DCartesian(scalar3, "scalar", axis.first, axis.last, 3,
			axis.first, axis.last, 3, axis.first, axis.last, 3, file.path.string()).success);
		CheckGrid(parameters3, axis);
		CheckGridRows<3>(ReadRows(file.Read(), 12, 4), axis);
	}
}

TEST_CASE("Vector exporters share finite coordinates with filtered and spherical paths", "[serializer][serializer-sampling]")
{
	for (const auto& axis : Axes()) {
		for (bool filtered : {false, true}) {
			SamplingFile file;
			std::vector<VectorN<Real, 2>> parameters2;
			VectorFunctionFromStdFunc<2> field2([&](const VectorN<Real, 2>& p) {
				parameters2.push_back(p); return VectorN<Real, 2>{0, 0};
			});
			const auto saved2 = filtered
				? Serializer::SaveVectorFunc2DCartesian(field2, "field", axis.first, axis.last, 3,
					axis.first, axis.last, 3, file.path.string(), 1)
				: Serializer::SaveVectorFunc2DCartesian(field2, "field", axis.first, axis.last, 3,
					axis.first, axis.last, 3, file.path.string());
			REQUIRE(saved2.success);
			CheckGrid(parameters2, axis);
			CheckGridRows<2>(ReadRows(file.Read(), 3, 4), axis);
			for (bool spherical : {false, true}) {
				std::vector<VectorN<Real, 3>> parameters3;
				VectorFunctionFromStdFunc<3> field3([&](const VectorN<Real, 3>& p) {
					parameters3.push_back(p); return VectorN<Real, 3>{0, 0, 0};
				});
				IOResult saved3;
				if (spherical)
					saved3 = filtered
						? Serializer::SaveVectorFuncSpherical(field3, "field", axis.first, axis.last, 3,
							axis.first, axis.last, 3, axis.first, axis.last, 3, file.path.string(), 1)
						: Serializer::SaveVectorFuncSpherical(field3, "field", axis.first, axis.last, 3,
							axis.first, axis.last, 3, axis.first, axis.last, 3, file.path.string());
				else
					saved3 = filtered
						? Serializer::SaveVectorFunc3DCartesian(field3, "field", axis.first, axis.last, 3,
							axis.first, axis.last, 3, axis.first, axis.last, 3, file.path.string(), 1)
						: Serializer::SaveVectorFunc3DCartesian(field3, "field", axis.first, axis.last, 3,
							axis.first, axis.last, 3, axis.first, axis.last, 3, file.path.string());
				REQUIRE(saved3.success);
				CheckGrid(parameters3, axis);
				CheckGridRows<3>(ReadRows(file.Read(), 3, 6), axis);
			}
		}

	}
}

TEST_CASE("Sampling grid axes retain independent units and sample counts", "[serializer][serializer-sampling]")
{
	const Real largest = std::numeric_limits<Real>::max();
	const Real largeUnit = std::ldexp(Real(1), std::numeric_limits<Real>::max_exponent - 3);
	const Real smallUnit = std::ldexp(Real(1),
		std::numeric_limits<Real>::min_exponent - std::numeric_limits<Real>::digits + 14);
	std::vector<VectorN<Real, 3>> parameters;
	VectorFunctionFromStdFunc<3> field([&](const VectorN<Real, 3>& p) {
		parameters.push_back(p); return VectorN<Real, 3>{1, 0, 0};
	});
	SamplingFile file;
	REQUIRE(Serializer::SaveVectorFunc3DCartesian(field, "mixed units", -largest, largest, 2,
		0, 4 * largeUnit, 5, -4 * smallUnit, 4 * smallUnit, 9, file.path.string(), 2).success);
	const auto rows = ReadRows(file.Read(), 3, 6);
	REQUIRE(parameters.size() == 2 * 5 * 9);
	REQUIRE(rows.size() == parameters.size());
	std::size_t row = 0;
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 5; ++j)
			for (int k = 0; k < 9; ++k, ++row) {
				const VectorN<Real, 3> expected{i == 0 ? -largest : largest, j * largeUnit, (k - 4) * smallUnit};
				for (int d = 0; d < 3; ++d) {
					CHECK(parameters[row][d] == expected[d]);
					CHECK(rows[row][d] == expected[d]);
				}
				CHECK(rows[row][3] == 1);
				CHECK(rows[row][4] == 0);
				CHECK(rows[row][5] == 0);
			}
}
