#include <catch2/catch_all.hpp>

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/tools/Persistence.h>
#endif

#include <chrono>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

using namespace MML;

namespace {
	void SetHeaderSize(std::string& bytes, uint16_t size) {
		bytes[10] = static_cast<char>(size & 0xff);
		bytes[11] = static_cast<char>((size >> 8) & 0xff);
	}

	void SetSchemaVersion(std::string& bytes, uint32_t version) {
		for (int i = 0; i < 4; ++i) bytes[16 + i] = static_cast<char>((version >> (8 * i)) & 0xff);
	}

	template<class Object>
	std::string Encode(const Object& object) {
		std::ostringstream stream(std::ios::out | std::ios::binary);
		const auto result = Persistence::SaveBinary(stream, object);
		INFO(result.message);
		REQUIRE(result.success);
		return stream.str();
	}

	template<class Object>
	void CheckSame(const Object& actual, const Object& expected) {
		CHECK(actual == expected);
	}

	template<class Type>
	void CheckSame(const Matrix<Type>& actual, const Matrix<Type>& expected) {
		REQUIRE(actual.rows() == expected.rows());
		REQUIRE(actual.cols() == expected.cols());
		for (int row = 0; row < actual.rows(); ++row)
			for (int col = 0; col < actual.cols(); ++col)
				CHECK(actual(row, col) == expected(row, col));
	}

	template<class Object>
	void CheckVersions(const Object& source, const Object& sentinel) {
		for (uint32_t version : {uint32_t(0), uint32_t(2), std::numeric_limits<uint32_t>::max()})
			for (bool strict : {false, true}) {
				CAPTURE(version, strict);
				auto bytes = Encode(source);
				SetSchemaVersion(bytes, version);
				std::istringstream input(bytes, std::ios::in | std::ios::binary);
				auto destination = sentinel;
				Persistence::LoadOptions options;
				options.strict_schema = strict;
				options.allow_scalar_conversion = true;
				const auto result = Persistence::LoadBinary(input, destination, options);
				CHECK_FALSE(result.success);
				CHECK(result.error == (version == 0 ? IOError::SCHEMA_MISMATCH : IOError::UNSUPPORTED_VERSION));
				CHECK_FALSE(result.message.empty());
				CheckSame(destination, sentinel);
				CHECK(input.tellg() == static_cast<std::streamoff>(Persistence::BinaryEnvelope::HEADER_SIZE));
			}
	}

	template<class Object>
	void CheckExtents(const Object& source, const Object& sentinel) {
		for (uint16_t extent : {uint16_t(0), uint16_t(47), uint16_t(49), uint16_t(64), std::numeric_limits<uint16_t>::max()})
			for (bool extended : {false, true}) {
				CAPTURE(extent, extended);
				auto bytes = Encode(source);
				SetHeaderSize(bytes, extent);
				if (extended && extent > Persistence::BinaryEnvelope::HEADER_SIZE)
					bytes.insert(Persistence::BinaryEnvelope::HEADER_SIZE,
						extent - Persistence::BinaryEnvelope::HEADER_SIZE, '\0');
				for (bool strict : {false, true}) {
					Persistence::LoadOptions options;
					options.strict_schema = strict;
					std::istringstream input(bytes, std::ios::in | std::ios::binary);
					auto destination = sentinel;
					const auto result = Persistence::LoadBinary(input, destination, options);
					CHECK_FALSE(result.success);
					CHECK(result.error == IOError::SCHEMA_MISMATCH);
					CHECK_FALSE(result.message.empty());
					CheckSame(destination, sentinel);
					CHECK(input.tellg() == static_cast<std::streamoff>(Persistence::BinaryEnvelope::HEADER_SIZE));
				}
			}
	}

	struct TemporaryBinaryFile {
		std::filesystem::path path = std::filesystem::temp_directory_path() /
			("mml_binary_schema_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".mmlb");
		~TemporaryBinaryFile() {
			std::error_code error;
			std::filesystem::remove(path, error);
			if (error) std::cerr << "Cannot remove binary schema fixture: " << error.message() << '\n';
		}
	};
}

TEMPLATE_TEST_CASE("Dense binary codecs reject unknown schemas and noncanonical header extents before mutation",
	"[Persistence][Binary][binary-schema]", float, double, long double, (std::complex<float>), (std::complex<double>)) {
	const Vector<TestType> vector{TestType(1), TestType(-2)};
	const Vector<TestType> vectorSentinel{TestType(9)};
	const Matrix<TestType> matrix(2, 2, {TestType(1), TestType(-2), TestType(3), TestType(4)});
	const Matrix<TestType> matrixSentinel(1, 1, {TestType(9)});
	SECTION("Vector schema") { CheckVersions(vector, vectorSentinel); }
	SECTION("Matrix schema") { CheckVersions(matrix, matrixSentinel); }
	SECTION("Vector header extent") { CheckExtents(vector, vectorSentinel); }
	SECTION("Matrix header extent") { CheckExtents(matrix, matrixSentinel); }
	if constexpr (is_std_complex<TestType>::value) {
		const std::vector<TestType> complex{TestType(1, -2), TestType(3, 4)};
		const std::vector<TestType> sentinel{TestType(9, 8)};
		SECTION("std complex Vector schema") { CheckVersions(complex, sentinel); }
		SECTION("std complex Vector header extent") { CheckExtents(complex, sentinel); }
	}
}

TEST_CASE("Binary envelopes retain future object schemas but require the v1 header extent",
	"[Persistence][Binary][binary-schema]") {
	Persistence::BinaryEnvelopeHeader header;
	header.object_schema_version = 2;
	std::ostringstream output(std::ios::out | std::ios::binary);
	REQUIRE(Persistence::WriteBinaryEnvelopeHeader(output, header).success);
	REQUIRE(output.str().size() == 48);
	std::istringstream input(output.str(), std::ios::in | std::ios::binary);
	const auto decoded = Persistence::ReadBinaryEnvelopeHeader(input);
	REQUIRE(decoded.result.success);
	CHECK(decoded.header.object_schema_version == 2);
	CHECK(input.tellg() == 48);
	for (uint16_t extent : {uint16_t(47), uint16_t(49), std::numeric_limits<uint16_t>::max()}) {
		header.header_size = extent;
		CHECK(Persistence::ValidateBinaryEnvelopeHeader(header).error == IOError::SCHEMA_MISMATCH);
		auto bytes = output.str();
		SetHeaderSize(bytes, extent);
		std::istringstream invalid(bytes, std::ios::in | std::ios::binary);
		CHECK(Persistence::ReadBinaryEnvelopeHeader(invalid).result.error == IOError::SCHEMA_MISMATCH);
	}
	std::istringstream truncated(output.str().substr(0, 47), std::ios::in | std::ios::binary);
	CHECK(Persistence::ReadBinaryEnvelopeHeader(truncated).result.error == IOError::TRUNCATED_INPUT);
}

TEST_CASE("Dense binary path and format facades preserve schema failures and filenames",
	"[Persistence][Binary][binary-schema]") {
	TemporaryBinaryFile fixture;
	const auto filename = fixture.path.string();
	auto check = [&](const auto& source, const auto& sentinel) {
		for (bool version : {false, true}) {
			auto bytes = Encode(source);
			if (version) SetSchemaVersion(bytes, 2);
			else SetHeaderSize(bytes, 49);
			{
				std::ofstream output(fixture.path, std::ios::binary);
				REQUIRE(output.is_open());
				output.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
				output.close();
				REQUIRE_FALSE(output.fail());
			}
			const auto expected = version ? IOError::UNSUPPORTED_VERSION : IOError::SCHEMA_MISMATCH;
			auto destination = sentinel;
			auto result = Persistence::LoadBinary(filename, destination);
			CHECK_FALSE(result.success);
			CHECK(result.error == expected);
			CHECK(result.filename == filename);
			CheckSame(destination, sentinel);
			result = Persistence::Load(filename, destination);
			CHECK_FALSE(result.success);
			CHECK(result.error == expected);
			CHECK(result.filename == filename);
			CheckSame(destination, sentinel);
			std::istringstream input(bytes, std::ios::in | std::ios::binary);
			result = Persistence::Load(input, destination, Persistence::Format::Binary);
			CHECK_FALSE(result.success);
			CHECK(result.error == expected);
			CheckSame(destination, sentinel);
		}
	};
	check(Vector<double>{1, 2}, Vector<double>{9});
	check(Matrix<double>(1, 2, {1, 2}), Matrix<double>(1, 1, {9}));
}

TEMPLATE_TEST_CASE("Canonical binary v1 empty and populated objects round-trip and compose in streams",
	"[Persistence][Binary][binary-schema]", float, double, long double, (std::complex<float>), (std::complex<double>)) {
	auto check = [](const auto& source) {
		using Object = std::decay_t<decltype(source)>;
		const auto bytes = Encode(source);
		CHECK(bytes[16] == 1);
		CHECK(bytes[17] == 0);
		CHECK(bytes[18] == 0);
		CHECK(bytes[19] == 0);
		std::istringstream headerInput(bytes, std::ios::in | std::ios::binary);
		const auto envelope = Persistence::ReadBinaryEnvelopeHeader(headerInput);
		REQUIRE(envelope.result.success);
		CHECK(envelope.header.header_size == 48);
		CHECK(bytes.size() == 48 + envelope.header.payload_byte_count);
		std::istringstream input(bytes + bytes, std::ios::in | std::ios::binary);
		Object first;
		auto result = Persistence::LoadBinary(input, first);
		INFO(result.message);
		REQUIRE(result.success);
		CheckSame(first, source);
		CHECK(input.tellg() == static_cast<std::streamoff>(bytes.size()));
		Object second;
		result = Persistence::LoadBinary(input, second);
		REQUIRE(result.success);
		CheckSame(second, source);
		CHECK(input.tellg() == static_cast<std::streamoff>(2 * bytes.size()));
	};
	check(Vector<TestType>{});
	check(Vector<TestType>{TestType(1), TestType(-2)});
	check(Matrix<TestType>{});
	check(Matrix<TestType>(1, 2, {TestType(1), TestType(-2)}));
	if constexpr (is_std_complex<TestType>::value) {
		const std::vector<TestType> empty;
		check(empty);
		check(std::vector<TestType>{TestType(1, -2), TestType(3, 4)});
	}
}

TEST_CASE("Binary sidecar and trailing-byte policies remain permissive for streams and files",
	"[Persistence][Binary][binary-schema]") {
	TemporaryBinaryFile fixture;
	auto check = [&](const auto& source) {
		using Object = std::decay_t<decltype(source)>;
		for (auto policy : {Persistence::BinarySidecarPolicy::None, Persistence::BinarySidecarPolicy::Optional,
			Persistence::BinarySidecarPolicy::Required}) {
			auto bytes = Encode(source);
			bytes[44] = static_cast<char>(policy);
			const auto payloadEnd = static_cast<std::streamoff>(bytes.size());
			bytes += "TRAILING";
			{
				std::ofstream output(fixture.path, std::ios::binary);
				REQUIRE(output.is_open());
				output.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
				output.close();
				REQUIRE_FALSE(output.fail());
			}
			for (bool strict : {false, true}) {
				CAPTURE(policy, strict);
				Persistence::LoadOptions options;
				options.strict_schema = strict;
				std::istringstream input(bytes, std::ios::in | std::ios::binary);
				Object loaded;
				REQUIRE(Persistence::LoadBinary(input, loaded, options).success);
				CheckSame(loaded, source);
				CHECK(input.tellg() == payloadEnd);
				std::string tail;
				REQUIRE(input >> tail);
				CHECK(tail == "TRAILING");
				REQUIRE(Persistence::LoadBinary(fixture.path.string(), loaded, options).success);
				CheckSame(loaded, source);
			}
		}
	};
	check(Vector<double>{1, -2});
	check(Matrix<double>(1, 2, {1, -2}));
}
