#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/tools/Persistence.h>
#endif
#include <catch2/catch_all.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <sstream>
#include <string>

using namespace MML;
using namespace MML::Persistence;

namespace {
	template<class T>
	IOResult SaveObject(std::ostream& out, const T& value) { return SaveJson(out, value); }
	IOResult SaveObject(std::ostream& out, const SampledRealFunctionData& value) { return SaveSampledFunction(out, value); }
	template<class T>
	IOResult LoadObject(std::istream& in, T& value, const LoadOptions& options) { return LoadJson(in, value, options); }
	IOResult LoadObject(std::istream& in, SampledRealFunctionData& value, const LoadOptions& options) {
		return LoadSampledFunctionData(in, value, options);
	}
	template<class T>
	std::string Encode(const T& value) {
		std::ostringstream out;
		REQUIRE(SaveObject(out, value).success);
		return out.str();
	}
	template<class T>
	JsonValue Root(const T& value) {
		const auto parsed = ParseJson(Encode(value));
		REQUIRE(parsed.result.success);
		return parsed.value;
	}
	template<class T>
	void CheckRejectedText(const T& before, const std::string& text, IOError expected, const LoadOptions& options = {}) {
		T destination = before;
		std::istringstream in(text);
		const auto result = LoadObject(in, destination, options);
		CHECK_FALSE(result.success);
		CHECK(result.error == expected);
		CHECK_FALSE(result.message.empty());
		CHECK(Encode(destination) == Encode(before));
	}
	template<class T>
	void CheckRejected(const T& before, const JsonValue& root, IOError expected, const LoadOptions& options = {}) {
		std::ostringstream out;
		REQUIRE(WriteJson(out, root).success);
		CheckRejectedText(before, out.str(), expected, options);
	}
	template<class T>
	void CheckHeaders(const T& value) {
		for (const auto& field : {"version", "object_version", "scalar_bytes"}) {
			auto root = Root(value);
			auto& number = root.object_value.at("mml").object_value.at(field);
			number = JsonValue::Number(number.number_value + 0.9);
			CheckRejected(value, root, std::string(field) == "scalar_bytes"
				? IOError::UNSUPPORTED_SCALAR : IOError::UNSUPPORTED_VERSION);
		}
	}
	template<class T>
	void CheckRealMetadata(const T& value) {
		for (const auto& invalid : {JsonValue::String("not-a-real-type"), JsonValue::Bool(true), JsonValue::Null()}) {
			auto root = Root(value);
			root.object_value.at("mml").object_value["real_type"] = invalid;
			CheckRejected(value, root, IOError::UNSUPPORTED_SCALAR);
		}
		auto root = Root(value);
		root.object_value.at("mml").object_value.erase("real_type");
		CheckRejected(value, root, IOError::UNSUPPORTED_SCALAR);
	}
	template<class Check>
	void EachObject(Check check) {
		check(Vector<Real>{41});
		check(VectorN<Real, 1>{41});
		Matrix<Real> matrix(1, 1);
		matrix(0, 0) = 41;
		check(matrix);
		MatrixNM<Real, 1, 1> fixed;
		fixed(0, 0) = 41;
		check(fixed);
		check(SampledRealFunctionData{Vector<Real>{0, 1}, Vector<Real>{41, 42}, "unchanged", "test"});
	}
	std::string WithRawNumber(const JsonValue& root, const std::string& token) {
		std::ostringstream out;
		REQUIRE(WriteJson(out, root).success);
		std::string text = out.str();
		const auto position = text.find("123.5");
		REQUIRE(position != std::string::npos);
		text.replace(position, 5, token);
		return text;
	}
	template<class Check>
	void EachDenseObject(Check check) {
		check(Vector<Real>{41});
		check(VectorN<Real, 1>{41});
		Matrix<Real> matrix(1, 1);
		matrix(0, 0) = 41;
		check(matrix);
		MatrixNM<Real, 1, 1> fixed;
		fixed(0, 0) = 41;
		check(fixed);
	}
	template<class Check>
	void EachFloatObject(Check check) {
		check(Vector<float>{11, 22});
		check(VectorN<float, 2>{11, 22});
		Matrix<float> matrix(1, 2);
		matrix(0, 0) = 11;
		matrix(0, 1) = 22;
		check(matrix);
		MatrixNM<float, 1, 2> fixed;
		fixed(0, 0) = 11;
		fixed(0, 1) = 22;
		check(fixed);
	}
	template<class T>
	void CheckAccepted(const T& before, const std::string& text, const LoadOptions& options = {}) {
		T destination = before;
		std::istringstream in(text);
		REQUIRE(LoadObject(in, destination, options).success);
		CHECK(Encode(destination) == Encode(before));
	}
}

TEST_CASE("JSON container object and scalar-size schema fields require integers", "[Persistence][JsonConversion]") {
	EachObject([](const auto& value) { CheckHeaders(value); });
}

TEST_CASE("Matrix JSON rejects fractional and negative dimensions without changing its target", "[Persistence][JsonConversion]") {
	Matrix<Real> matrix(1, 1);
	matrix(0, 0) = 41;
	for (double invalid : {1.9, -0.1}) {
		for (int dimension : {0, 1}) {
			auto root = Root(matrix);
			root.object_value.at("shape").array_value[dimension] = JsonValue::Number(invalid);
			if (invalid < 0)
				root.object_value.at("data") = JsonValue::ArrayValue();
			CheckRejected(matrix, root, IOError::SCHEMA_MISMATCH);
		}
	}
}

TEST_CASE("Logical Real JSON requires its emitted concrete backing-type metadata", "[Persistence][JsonConversion]") {
	EachObject([](const auto& value) { CheckRealMetadata(value); });
}

TEST_CASE("JSON schema rejects fractions hidden by double rounding and accepts integral decimal notation", "[Persistence][JsonConversion]") {
	EachObject([](const auto& value) {
		for (const char* token : {"1.00000000000000000001", "0.99999999999999999999", "10000000000000000001e-19"}) {
			auto root = Root(value);
			root.object_value.at("mml").object_value.at("version") = JsonValue::Number(123.5);
			CheckRejectedText(value, WithRawNumber(root, token), IOError::UNSUPPORTED_VERSION);
		}
		for (const char* token : {"1.0", "10e-1", "0.1e1", "1000e-3", "1e+0000000000", "1.0e-00000"}) {
			auto root = Root(value);
			root.object_value.at("mml").object_value.at("version") = JsonValue::Number(123.5);
			CheckAccepted(value, WithRawNumber(root, token));
		}
	});
	EachDenseObject([](const auto& value) {
		auto root = Root(value);
		root.object_value.at("shape").array_value[0] = JsonValue::Number(123.5);
		CheckRejectedText(value, WithRawNumber(root, "1.00000000000000000001"), IOError::SCHEMA_MISMATCH);
		CheckAccepted(value, WithRawNumber(root, "10e-1"));
	});
}

TEST_CASE("JSON dimension checks reject unsafe values before casting or allocating", "[Persistence][JsonConversion]") {
	for (const char* token : {"-0.1", "0.1", "2147483647.00000000001", "2147483648", "1e100"}) {
		const auto parsed = ParseJson(token);
		REQUIRE(parsed.result.success);
		int dimension = 77;
		const auto result = Persistence::Detail::ReadJsonDimension(parsed.value, dimension);
		CHECK_FALSE(result.success);
		CHECK(result.error == ((std::string(token) == "2147483648" || std::string(token) == "1e100")
			? IOError::ALLOCATION_LIMIT_EXCEEDED : IOError::SCHEMA_MISMATCH));
		CHECK(dimension == 77);
	}
	for (double invalid : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
		-std::numeric_limits<double>::infinity()}) {
		int dimension = 77;
		const auto number = JsonValue::Number(invalid);
		CHECK(Persistence::Detail::ReadJsonDimension(number, dimension).error == IOError::SCHEMA_MISMATCH);
		CHECK(dimension == 77);
		CHECK_FALSE(Persistence::Detail::JsonIntegerEquals(&number, 1));
	}
	for (const char* token : {"0", "-0", "0.000e-10", "2147483647", "21474836470e-1"}) {
		const auto parsed = ParseJson(token);
		REQUIRE(parsed.result.success);
		int dimension = -1;
		CHECK(Persistence::Detail::ReadJsonDimension(parsed.value, dimension).success);
		CHECK(dimension == (token[0] == '2' ? std::numeric_limits<int>::max() : 0));
	}
	EachDenseObject([](const auto& value) {
		auto root = Root(value);
		root.object_value.at("shape").array_value[0] = JsonValue::Number(1e100);
		CheckRejected(value, root, IOError::ALLOCATION_LIMIT_EXCEEDED);
	});
}

TEST_CASE("JSON metadata rejects large numeric fields and mismatched concrete Real under every option", "[Persistence][JsonConversion]") {
	EachObject([](const auto& value) {
		for (const char* field : {"version", "object_version", "scalar_bytes"}) {
			auto root = Root(value);
			root.object_value.at("mml").object_value.at(field) = JsonValue::Number(1e100);
			CheckRejected(value, root, std::string(field) == "scalar_bytes"
				? IOError::UNSUPPORTED_SCALAR : IOError::UNSUPPORTED_VERSION);
		}
		for (bool convert : {false, true}) {
			LoadOptions options;
			options.allow_scalar_conversion = convert;
			options.strict_schema = false;
			auto root = Root(value);
			const auto expected = root.object_value.at("mml").object_value.at("real_type").string_value;
			for (const char* type : {"float", "double", "long_double"}) {
				if (expected == type) continue;
				root.object_value.at("mml").object_value.at("real_type") = JsonValue::String(type);
				CheckRejected(value, root, IOError::UNSUPPORTED_SCALAR, options);
			}
			root = Root(value);
			root.object_value.at("mml").object_value.at("scalar") = JsonValue::String("not-the-target-type");
			CheckRejected(value, root, IOError::UNSUPPORTED_SCALAR, options);
		}
	});
}

TEST_CASE("JSON float payload ranges are validated before any destination changes", "[Persistence][JsonConversion]") {
	const double maximum = std::numeric_limits<float>::max();
	EachFloatObject([&](const auto& value) {
		for (double invalid : {1e100, -1e100, std::nextafter(maximum, std::numeric_limits<double>::infinity()),
			std::nextafter(-maximum, -std::numeric_limits<double>::infinity())}) {
			auto root = Root(value);
			root.object_value.at("data").array_value[0] = JsonValue::Number(0.5);
			root.object_value.at("data").array_value[1] = JsonValue::Number(invalid);
			CheckRejected(value, root, IOError::UNSUPPORTED_SCALAR);
			LoadOptions options;
			options.allow_scalar_conversion = true;
			CheckRejected(value, root, IOError::UNSUPPORTED_SCALAR, options);
		}
	});
	for (double invalid : {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()})
		CHECK_FALSE(Persistence::Detail::JsonScalarFits<float>(JsonValue::Number(invalid)));
}

TEST_CASE("JSON float payloads preserve rounding endpoints and user-approved underflow", "[Persistence][JsonConversion]") {
	for (double number : {0.1, double(std::numeric_limits<float>::max()), -double(std::numeric_limits<float>::max()),
		double(std::numeric_limits<float>::denorm_min()), 1e-50, -1e-50, 0.0, -0.0}) {
		Vector<float> before{11, 22};
		auto root = Root(before);
		root.object_value.at("data").array_value[0] = JsonValue::Number(number);
		std::ostringstream out;
		REQUIRE(WriteJson(out, root).success);
		std::istringstream in(out.str());
		REQUIRE(LoadJson(in, before).success);
		CHECK(before[0] == static_cast<float>(number));
		CHECK(before[1] == 22);
		if (before[0] == 0) CHECK(std::signbit(before[0]) == std::signbit(number));
	}
	Vector<float> value{0.123456f, -0.765432f};
	SaveOptions options;
	options.precision = 2;
	std::stringstream out;
	REQUIRE(SaveJson(out, value, options).success);
	Vector<float> loaded;
	REQUIRE(LoadJson(out, loaded).success);
	CHECK(loaded[0] == float(0.12));
	CHECK(loaded[1] == float(-0.77));
}

TEST_CASE("Sampled-function JSON validates every scalar-bearing field", "[Persistence][JsonConversion]") {
	const SampledRealFunctionData original{Vector<Real>{0, 1}, Vector<Real>{41, 42}, "unchanged", "test"};
	for (const char* field : {"domain", "nodes", "values"}) {
		auto root = Root(original);
		root.object_value.at(field).array_value[0] = JsonValue::Number(std::numeric_limits<double>::infinity());
		const auto result = Persistence::Detail::ValidateSampledFunctionJsonRoot(root, {});
		CHECK(result.error == IOError::UNSUPPORTED_SCALAR);
	}
	LoadOptions options;
	options.max_allocation_bytes = 0;
	CheckRejected(original, Root(original), IOError::ALLOCATION_LIMIT_EXCEEDED, options);
}

TEST_CASE("JSON allocation limits and fixed shapes remain enforced transactionally", "[Persistence][JsonConversion]") {
	EachDenseObject([](const auto& value) {
		LoadOptions options;
		options.max_allocation_bytes = 0;
		CheckRejected(value, Root(value), IOError::ALLOCATION_LIMIT_EXCEEDED, options);
		auto root = Root(value);
		root.object_value.at("shape").array_value[0] = JsonValue::Number(2);
		CheckRejected(value, root, IOError::SCHEMA_MISMATCH);
	});
}

TEST_CASE("JSON decimal integer classification agrees with an independent rational oracle", "[Persistence][JsonConversion]") {
	for (std::uint64_t coefficient : {0ull, 1ull, 10ull, 12ull, 100ull, 1200ull, 123450ull, 9999999ull}) {
		for (int fractional = 0; fractional <= 6; ++fractional) {
			std::string mantissa = std::to_string(coefficient);
			if (fractional > 0) {
				if (mantissa.size() <= static_cast<std::size_t>(fractional))
					mantissa.insert(0, fractional + 1 - mantissa.size(), '0');
				mantissa.insert(mantissa.size() - fractional, ".");
			}
			for (int exponent = -8; exponent <= 8; ++exponent) {
				std::uint64_t denominator = 1;
				for (int scale = 0; scale < fractional - exponent; ++scale) denominator *= 10;
				const bool integral = coefficient % denominator == 0;
				for (bool negative : {false, true}) {
					const std::string token = (negative ? "-" : "") + mantissa + "e" + std::to_string(exponent);
					CAPTURE(token);
					const auto parsed = ParseJson(token);
					REQUIRE(parsed.result.success);
					CHECK(Persistence::Detail::IsJsonInteger(parsed.value) == integral);
				}
			}
		}
	}
}
