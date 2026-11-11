#include <catch2/catch_all.hpp>

#include <mml/tools/Persistence.h>

#include <cmath>
#include <chrono>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <type_traits>

using namespace MML;

namespace {
	struct CommaDecimal : std::numpunct<char> {
		char do_decimal_point() const override { return ','; }
		char do_thousands_sep() const override { return '.'; }
		std::string do_grouping() const override { return "\3"; }
	};

	template<class Type> struct Component { using type = Type; };
	template<class Type> struct Component<std::complex<Type>> { using type = Type; };

	template<class Type>
	Matrix<Type> PrecisionMatrix() {
		using Scalar = typename Component<Type>::type;
		const Scalar next = std::nextafter(Scalar(1), Scalar(2));
		Matrix<Type> matrix(1, 3);
		if constexpr (is_std_complex<Type>::value) {
			matrix(0, 0) = Type(next, -next);
			matrix(0, 1) = Type(std::sqrt(Scalar(2)), Scalar(1) / Scalar(7));
			matrix(0, 2) = Type(std::numeric_limits<Scalar>::min() * Scalar(8), -Scalar(1));
		} else {
			matrix(0, 0) = next;
			matrix(0, 1) = std::sqrt(Scalar(2));
			matrix(0, 2) = std::numeric_limits<Scalar>::min() * Scalar(8);
		}
		return matrix;
	}

	struct TemporaryDirectory {
		std::filesystem::path path = std::filesystem::temp_directory_path() /
			("mml_portability_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
		TemporaryDirectory() { std::filesystem::create_directory(path); }
		~TemporaryDirectory() {
			std::error_code error;
			std::filesystem::remove_all(path, error);
			if (error) std::cerr << "Failed to clean persistence test files: " << error.message() << '\n';
		}
	};

	std::string LongDoubleRecord(uint16_t precision, int32_t exponent, uint64_t low, uint64_t high = 0) {
		std::ostringstream record(std::ios::out | std::ios::binary);
		record.put(0);
		record.put(0);
		Persistence::WriteUInt16LE(record, precision);
		Persistence::WriteUInt32LE(record, static_cast<uint32_t>(exponent));
		Persistence::WriteUInt64LE(record, low);
		Persistence::WriteUInt64LE(record, high);
		return record.str();
	}
}

TEMPLATE_TEST_CASE("Matrix text defaults preserve the actual scalar precision and caller format",
	"[Persistence][Portability][MatrixIO]", float, double, long double, (std::complex<float>), (std::complex<double>)) {
	const auto matrix = PrecisionMatrix<TestType>();
	std::stringstream stream;
	const std::locale locale(std::locale::classic(), new CommaDecimal);
	stream.imbue(locale);
	stream << std::hexfloat << std::showpos << std::uppercase << std::setprecision(2) << std::setfill('#') << std::setw(9);
	const auto flags = stream.flags();
	REQUIRE(Persistence::Save(stream, matrix, Persistence::Format::Text).success);
	REQUIRE(stream.getloc() == locale);
	REQUIRE(stream.flags() == flags);
	REQUIRE(stream.precision() == 2);
	REQUIRE(stream.fill() == '#');
	REQUIRE(stream.width() == 9);

	Matrix<TestType> loaded;
	const auto result = Persistence::Load(stream, loaded, Persistence::Format::Text);
	INFO(result.message);
	REQUIRE(result.success);
	REQUIRE(stream.getloc() == locale);
	REQUIRE(stream.flags() == flags);
	REQUIRE(stream.width() == 9);
	REQUIRE(loaded.rows() == 1);
	REQUIRE(loaded.cols() == 3);
	for (int col = 0; col < 3; ++col)
		REQUIRE(loaded(0, col) == matrix(0, col));
}

TEMPLATE_TEST_CASE("Matrix CSV defaults round-trip precise scalar values", "[Persistence][Portability][MatrixIO]",
	float, double, long double, (std::complex<float>), (std::complex<double>)) {
	const auto matrix = PrecisionMatrix<TestType>();
	std::stringstream stream;
	REQUIRE(Persistence::Save(stream, matrix, Persistence::Format::CSV).success);
	Matrix<TestType> loaded;
	const auto result = Persistence::Load(stream, loaded, Persistence::Format::CSV);
	INFO(result.message);
	REQUIRE(result.success);
	REQUIRE(loaded.rows() == 1);
	REQUIRE(loaded.cols() == 3);
	for (int col = 0; col < 3; ++col)
		REQUIRE(loaded(0, col) == matrix(0, col));
}

TEST_CASE("Matrix save precision is automatic by default and explicit when requested", "[Persistence][Portability][MatrixIO]") {
	Persistence::SaveOptions options;
	REQUIRE(options.precision == 0);
	const auto matrix = PrecisionMatrix<double>();
	options.precision = 6;
	std::stringstream stream;
	stream.precision(3);
	REQUIRE(Persistence::Save(stream, matrix, Persistence::Format::Text, options).success);
	REQUIRE(stream.precision() == 3);
	Matrix<double> loaded;
	REQUIRE(Persistence::Load(stream, loaded, Persistence::Format::Text).success);
	REQUIRE(loaded(0, 0) != matrix(0, 0));

	options.precision = -1;
	std::ostringstream invalid;
	const auto result = Persistence::Save(invalid, matrix, Persistence::Format::Text, options);
	REQUIRE_FALSE(result.success);
	REQUIRE(result.error == IOError::INVALID_PARAMETERS);
	REQUIRE(invalid.str().empty());
}

TEST_CASE("Matrix binary load forwards caller allocation limits", "[Persistence][Portability][Binary]") {
	const auto matrix = PrecisionMatrix<double>();
	std::stringstream stream(std::ios::in | std::ios::out | std::ios::binary);
	REQUIRE(Persistence::Save(stream, matrix, Persistence::Format::Binary).success);
	Persistence::LoadOptions options;
	options.max_allocation_bytes = 1;
	Matrix<double> loaded;
	const auto result = Persistence::Load(stream, loaded, Persistence::Format::Binary, options);
	REQUIRE_FALSE(result.success);
	REQUIRE(result.error == IOError::ALLOCATION_LIMIT_EXCEEDED);

	TemporaryDirectory directory;
	const auto filename = (directory.path / "matrix.mmlb").string();
	REQUIRE(Persistence::Save(matrix, filename).success);
	const auto fileResult = Persistence::Load(filename, loaded, options);
	REQUIRE_FALSE(fileResult.success);
	REQUIRE(fileResult.error == IOError::ALLOCATION_LIMIT_EXCEEDED);
}

TEMPLATE_TEST_CASE("Matrix file helpers and option overloads use round-trip defaults", "[Persistence][Portability][MatrixIO]",
	float, double, long double, (std::complex<float>), (std::complex<double>)) {
	TemporaryDirectory directory;
	const auto matrix = PrecisionMatrix<TestType>();
	for (const auto* extension : {".txt", ".csv"}) {
		const auto filename = (directory.path / (std::string("matrix") + extension)).string();
		REQUIRE(Persistence::Save(matrix, filename).success);
		Matrix<TestType> loaded;
		REQUIRE(Persistence::Load(filename, loaded).success);
		for (int col = 0; col < 3; ++col)
			REQUIRE(loaded(0, col) == matrix(0, col));
		REQUIRE(Persistence::Save(matrix, filename, Persistence::SaveOptions{}).success);
		REQUIRE(Persistence::Load(filename, loaded, Persistence::LoadOptions{}).success);
		for (int col = 0; col < 3; ++col)
			REQUIRE(loaded(0, col) == matrix(0, col));
	}
}

TEST_CASE("Binary scalar payloads have explicit little-endian IEEE widths", "[Persistence][Portability][Binary]") {
	std::ostringstream floatBytes(std::ios::out | std::ios::binary);
	REQUIRE(Persistence::Detail::WriteBinaryScalarLE(floatBytes, std::complex<float>(1, -2)).success);
	REQUIRE(floatBytes.str() == std::string("\x00\x00\x80\x3f\x00\x00\x00\xc0", 8));
	std::ostringstream doubleBytes(std::ios::out | std::ios::binary);
	REQUIRE(Persistence::Detail::WriteBinaryScalarLE(doubleBytes, std::complex<double>(1, -2)).success);
	REQUIRE(doubleBytes.str() == std::string("\x00\x00\x00\x00\x00\x00\xf0\x3f\x00\x00\x00\x00\x00\x00\x00\xc0", 16));

	Vector<long double> values{1};
	std::stringstream stream(std::ios::in | std::ios::out | std::ios::binary);
	REQUIRE(Persistence::SaveBinary(stream, values).success);
	REQUIRE(stream.str().size() == 48 + 8 + 24);
	const auto header = Persistence::ReadBinaryEnvelopeHeader(stream);
	REQUIRE(header.result.success);
	REQUIRE(header.header.scalar_byte_size == 24);
	REQUIRE(header.header.payload_byte_count == 8 + 24);
	auto bytes = stream.str();
	bytes[24] = 8; // A native-width record is not this canonical wire encoding.
	std::istringstream wrongWidth(bytes, std::ios::in | std::ios::binary);
	Vector<long double> loaded{9};
	const auto result = Persistence::LoadBinary(wrongWidth, loaded);
	REQUIRE_FALSE(result.success);
	REQUIRE(result.error == IOError::UNSUPPORTED_SCALAR);
	REQUIRE(loaded.size() == 1);
	REQUIRE(loaded[0] == 9);
}

TEST_CASE("Canonical long double reads foreign precision and rejects malformed records", "[Persistence][Portability][Binary]") {
	const auto one = LongDoubleRecord(24, 1, uint64_t{1} << 23);
	std::istringstream foreign(one, std::ios::in | std::ios::binary);
	long double value = 9;
	REQUIRE(Persistence::Detail::ReadFloatingScalarLE(foreign, value).success);
	REQUIRE(value == 1);
	const auto rejects = [&](std::string bytes, IOError error) {
		std::istringstream stream(bytes, std::ios::in | std::ios::binary);
		value = 9;
		const auto result = Persistence::Detail::ReadFloatingScalarLE(stream, value);
		INFO(result.message);
		REQUIRE_FALSE(result.success);
		REQUIRE(result.error == error);
		REQUIRE(value == 9);
	};
	auto badClass = one;
	badClass[0] = 3;
	rejects(badClass, IOError::INVALID_FORMAT);
	auto badSign = one;
	badSign[1] = 2;
	rejects(badSign, IOError::INVALID_FORMAT);
	rejects(LongDoubleRecord(0, 0, 0), IOError::INVALID_FORMAT);
	rejects(LongDoubleRecord(129, 0, 0), IOError::INVALID_FORMAT);
	rejects(LongDoubleRecord(24, 1, 0), IOError::INVALID_FORMAT);
	rejects(LongDoubleRecord(24, 1, 1), IOError::INVALID_FORMAT);
	rejects(LongDoubleRecord(24, 1, uint64_t{1} << 24), IOError::INVALID_FORMAT);
	auto nonfinitePayload = one;
	nonfinitePayload[0] = 1;
	rejects(nonfinitePayload, IOError::INVALID_FORMAT);
	rejects(LongDoubleRecord(24, std::numeric_limits<int32_t>::max(), uint64_t{1} << 23), IOError::UNSUPPORTED_SCALAR);
	rejects(LongDoubleRecord(24, std::numeric_limits<int32_t>::min(), uint64_t{1} << 23), IOError::UNSUPPORTED_SCALAR);
	const int minimumExponent = std::numeric_limits<long double>::min_exponent - std::numeric_limits<long double>::digits + 1;
	rejects(LongDoubleRecord(24, minimumExponent, (uint64_t{1} << 23) | (uint64_t{1} << 22)), IOError::UNSUPPORTED_SCALAR);
	if constexpr (std::numeric_limits<long double>::digits < 128)
		rejects(LongDoubleRecord(128, 1, 0, uint64_t{1} << 63), IOError::UNSUPPORTED_SCALAR);
	rejects(one.substr(0, 23), IOError::TRUNCATED_INPUT);
}

TEST_CASE("Canonical formatting restores caller state on write failure and exception", "[Persistence][Portability][MatrixIO]") {
	struct FailingBuffer : std::streambuf {
		int_type overflow(int_type) override { return traits_type::eof(); }
	} buffer;
	std::ostream stream(&buffer);
	const std::locale locale(std::locale::classic(), new CommaDecimal);
	stream.imbue(locale);
	stream << std::fixed << std::setprecision(2) << std::setfill('#') << std::setw(11);
	const auto flags = stream.flags();
	const auto matrix = PrecisionMatrix<double>();
	REQUIRE_FALSE(Persistence::Save(stream, matrix, Persistence::Format::Text).success);
	REQUIRE(stream.flags() == flags);
	REQUIRE(stream.getloc() == locale);
	REQUIRE(stream.precision() == 2);
	REQUIRE(stream.fill() == '#');
	REQUIRE(stream.width() == 11);
	stream.clear();
	stream.exceptions(std::ios::badbit | std::ios::failbit);
	REQUIRE_THROWS_AS(Persistence::Save(stream, matrix, Persistence::Format::Text), std::ios_base::failure);
	REQUIRE(stream.flags() == flags);
	REQUIRE(stream.getloc() == locale);
	REQUIRE(stream.precision() == 2);
	REQUIRE(stream.width() == 11);
}
