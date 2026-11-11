#include <catch2/catch_all.hpp>

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/tools/Persistence.h>
#endif

#include <chrono>
#include <complex>
#include <filesystem>
#include <fstream>
#include <locale>
#include <sstream>
#include <streambuf>
#include <string>

using namespace MML;
using namespace MML::Persistence;

namespace {
	struct TempFile {
		std::filesystem::path path;
		explicit TempFile(const std::string& extension = ".txt")
			: path(std::filesystem::temp_directory_path() /
			       ("mml_vector_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + extension)) {}
		~TempFile() {
			std::error_code error;
			std::filesystem::remove(path, error);
			if (error) std::cerr << "Cannot remove vector test file: " << error.message() << '\n';
		}
	};

	template<class Type>
	void CheckValue(const Type& actual, const Type& expected) {
		if constexpr (is_std_complex<Type>::value) {
			CheckValue(actual.real(), expected.real());
			CheckValue(actual.imag(), expected.imag());
		} else if constexpr (std::is_floating_point_v<Type>) {
			if (std::isnan(expected))
				REQUIRE(std::isnan(actual));
			else {
				REQUIRE(actual == expected);
				if (expected == 0 || std::isinf(expected))
					REQUIRE(std::signbit(actual) == std::signbit(expected));
			}
		} else
			REQUIRE(actual == expected);
	}

	template<class Type>
	void CheckVector(const Vector<Type>& actual, const Vector<Type>& expected) {
		REQUIRE(actual.size() == expected.size());
		for (int i = 0; i < actual.size(); ++i)
			CheckValue(actual[i], expected[i]);
	}

	template<class Type>
	Vector<Type> BoundaryValues() {
		if constexpr (is_std_complex<Type>::value) {
			using Part = typename Type::value_type;
			return Vector<Type>{Type(std::nextafter(Part(1), Part(2)), -Part(0)),
			        Type(std::numeric_limits<Part>::denorm_min(), std::numeric_limits<Part>::max()),
			        Type(std::numeric_limits<Part>::infinity(), -std::numeric_limits<Part>::infinity()),
			        Type(std::numeric_limits<Part>::quiet_NaN(), Part(1) / Part(3))};
		} else if constexpr (std::is_floating_point_v<Type>)
			return Vector<Type>{Type(0), -Type(0), Type(1) / Type(3), std::nextafter(Type(1), Type(2)),
			        std::numeric_limits<Type>::min(), std::numeric_limits<Type>::denorm_min(),
			        std::numeric_limits<Type>::max(), std::numeric_limits<Type>::lowest(),
			        std::numeric_limits<Type>::infinity(), -std::numeric_limits<Type>::infinity(),
			        std::numeric_limits<Type>::quiet_NaN()};
		else
			return Vector<Type>{std::numeric_limits<Type>::min(), Type(0), Type(1), std::numeric_limits<Type>::max()};
	}

	struct CommaDecimal : std::numpunct<char> {
		char do_decimal_point() const override { return ','; }
	};

	struct ThrowingOutput : std::streambuf {
		std::streamsize xsputn(const char*, std::streamsize) override { throw std::runtime_error("output failed"); }
		int_type overflow(int_type) override { throw std::runtime_error("output failed"); }
	};
}

TEMPLATE_TEST_CASE("Vector text preserves exact scalar values in files and streams",
                   "[VectorPersistence][Text]", float, double, long double,
                   std::complex<float>, std::complex<double>, std::complex<long double>,
                   int, unsigned int, long long, unsigned long long, signed char, unsigned char)
{
	const auto original = BoundaryValues<TestType>();
	std::stringstream stream;
	REQUIRE(Save(stream, original, Format::Text).success);
	Vector<TestType> restored;
	REQUIRE(Load(stream, restored, Format::Text).success);
	CheckVector(restored, original);
	TempFile file;
	REQUIRE(Save(original, file.path.string()).success);
	Vector<TestType> fromFile;
	REQUIRE(Load(file.path.string(), fromFile).success);
	CheckVector(fromFile, original);
	std::ifstream disk(file.path);
	std::ostringstream contents;
	contents << disk.rdbuf();
	// Windows text mode normalizes CRLF on input.
	REQUIRE(contents.str() == stream.str());
}

TEST_CASE("Vector text retains count-prefixed layout, empty vectors and whitespace parsing", "[VectorPersistence][Text]")
{
	std::ostringstream out;
	REQUIRE(Save(out, Vector<int>{1, -2, 3}, Format::Text).success);
	REQUIRE(out.str() == "3\n1\n-2\n3\n");
	std::ostringstream empty;
	REQUIRE(Save(empty, Vector<double>{}, Format::Text).success);
	REQUIRE(empty.str() == "0\n");
	TempFile emptyFile;
	REQUIRE(Save(Vector<double>{}, emptyFile.path.string()).success);
	Vector<double> fileDestination{9};
	REQUIRE(Load(emptyFile.path.string(), fileDestination).success);
	REQUIRE(fileDestination.empty());
	for (const char* encoded : {"0", "0\n \t", "3 1 -2 +3", "+3\n1\n-2\n3\n"}) {
		std::istringstream in(encoded);
		Vector<int> destination{99};
		REQUIRE(Load(in, destination, Format::Text).success);
		REQUIRE(destination == (encoded[0] == '0' ? Vector<int>{} : Vector<int>{1, -2, 3}));
	}
	std::istringstream complex("2\n( 1 , -0 )\n( inf , nan )");
	Vector<Complex> values;
	REQUIRE(Load(complex, values, Format::Text).success);
	REQUIRE(values[0].real() == 1);
	REQUIRE(std::signbit(values[0].imag()));
	REQUIRE(std::isinf(values[1].real()));
	REQUIRE(std::isnan(values[1].imag()));
}

TEST_CASE("Vector text rejects malformed or truncated data without changing destination", "[VectorPersistence][Errors]")
{
	for (const char* encoded : {"", "-1\n", "1.5\n", "hello\n", "1\n", "2\n1\n", "1\nabc\n",
	                            "1\n1e99999\n", "1\n4x\n", "1\n4\nextra", "0\nextra", "1\n+-1", "1\n++1"}) {
		INFO(encoded);
		std::istringstream in(encoded);
		Vector<double> destination{11, 22};
		const auto result = Load(in, destination, Format::Text);
		REQUIRE_FALSE(result.success);
		REQUIRE(result.offset_kind == IOOffsetKind::Line);
		REQUIRE(result.offset >= 1);
		REQUIRE_FALSE(result.message.empty());
		REQUIRE(destination == Vector<double>{11, 22});
	}
	std::istringstream invalid("2\n1\nwrong\n");
	Vector<double> values{7};
	const auto result = Load(invalid, values, Format::Text);
	REQUIRE(result.error == IOError::MALFORMED_INPUT);
	REQUIRE(result.offset == 3);
	std::istringstream incomplete("2\n1\n");
	const auto truncated = Load(incomplete, values, Format::Text);
	REQUIRE(truncated.error == IOError::TRUNCATED_INPUT);
	REQUIRE(truncated.offset == 3);
	std::istringstream negativeUnsigned("1\n-1\n");
	Vector<unsigned int> unsignedValues{9};
	REQUIRE(Load(negativeUnsigned, unsignedValues, Format::Text).error == IOError::MALFORMED_INPUT);
	CheckVector(unsignedValues, Vector<unsigned int>{9});
	std::istringstream permissive("1\n4\nextra");
	LoadOptions options;
	options.strict_schema = false;
	REQUIRE(Load(permissive, values, Format::Text, options).success);
	REQUIRE(values == Vector<double>{4});
}

TEST_CASE("Complex vector text validates delimiters, components and completeness", "[VectorPersistence][Errors]")
{
	for (const char* encoded : {"1\n1,2\n", "1\n(1)\n", "1\n(1,)\n", "1\n(1,2\n", "1\n(1,wrong)\n",
	                            "1\n(1;2)\n", "1\n(1,2)extra"}) {
		std::istringstream in(encoded);
		Vector<Complex> destination{Complex(3, 4)};
		REQUIRE_FALSE(Load(in, destination, Format::Text).success);
		REQUIRE(destination == Vector<Complex>{Complex(3, 4)});
	}
}

TEST_CASE("Vector text validates counts and configured payload allocation before resizing", "[VectorPersistence][Allocation]")
{
	for (const char* count : {"18446744073709551616", "2147483648", "99999999999999999999999999999999"}) {
		std::istringstream in(count);
		Vector<double> destination{9};
		const auto result = Load(in, destination, Format::Text);
		REQUIRE(result.error == IOError::ALLOCATION_LIMIT_EXCEEDED);
		REQUIRE(destination == Vector<double>{9});
	}
	LoadOptions options;
	options.max_allocation_bytes = sizeof(double);
	std::istringstream rejected("2\n1\n2\n");
	Vector<double> destination{9};
	REQUIRE(Load(rejected, destination, Format::Text, options).error == IOError::ALLOCATION_LIMIT_EXCEEDED);
	REQUIRE(destination == Vector<double>{9});
	std::istringstream accepted("1\n2\n");
	REQUIRE(Load(accepted, destination, Format::Text, options).success);
	options.max_allocation_bytes = 0;
	std::istringstream empty("0");
	REQUIRE(Load(empty, destination, Format::Text, options).success);
	REQUIRE(destination.empty());
}

TEST_CASE("Vector text preserves stream formatting and locale on success, failure and exceptions", "[VectorPersistence][StreamState]")
{
	std::stringstream stream;
	stream.imbue(std::locale(std::locale::classic(), new CommaDecimal));
	stream.setf(std::ios::scientific | std::ios::showpos | std::ios::uppercase);
	stream.precision(2);
	stream.width(17);
	stream.fill('#');
	const auto flags = stream.flags();
	const auto locale = stream.getloc();
	Vector<double> original{std::nextafter(1.0, 2.0), 1.0 / 3};
	REQUIRE(Save(stream, original, Format::Text).success);
	REQUIRE(stream.flags() == flags);
	REQUIRE(stream.getloc() == locale);
	REQUIRE(stream.precision() == 2);
	REQUIRE(stream.width() == 17);
	REQUIRE(stream.fill() == '#');
	Vector<double> restored;
	REQUIRE(Load(stream, restored, Format::Text).success);
	CheckVector(restored, original);
	REQUIRE(stream.flags() == flags);
	REQUIRE(stream.getloc() == locale);
	REQUIRE(stream.width() == 17);
	REQUIRE(stream.fill() == '#');
	std::istringstream malformed("1\nbad\n");
	malformed.imbue(locale);
	malformed.precision(3);
	REQUIRE_FALSE(Load(malformed, restored, Format::Text).success);
	REQUIRE(malformed.getloc() == locale);
	REQUIRE(malformed.precision() == 3);
	ThrowingOutput buffer;
	std::ostream failing(&buffer);
	failing.precision(3);
	failing.fill('@');
	failing.exceptions(std::ios::badbit);
	REQUIRE_THROWS(Save(failing, original, Format::Text));
	REQUIRE(failing.precision() == 3);
	REQUIRE(failing.fill() == '@');
	std::istringstream throwing("1\n");
	throwing.precision(3);
	throwing.exceptions(std::ios::eofbit | std::ios::failbit | std::ios::badbit);
	const auto before = restored;
	REQUIRE_THROWS(Load(throwing, restored, Format::Text));
	CheckVector(restored, before);
	REQUIRE(throwing.precision() == 3);
	std::ostringstream broken;
	broken.setstate(std::ios::badbit);
	REQUIRE(Save(broken, original, Format::Text).error == IOError::WRITE_FAILED);
	std::istringstream unreadable("1\n1\n");
	unreadable.setstate(std::ios::badbit);
	REQUIRE(Load(unreadable, restored, Format::Text).error == IOError::READ_FAILED);
	CheckVector(restored, before);
}

TEST_CASE("Vector options affect text and JSON and invalid requests do not clobber files", "[VectorPersistence][Options]")
{
	Vector<double> original{1.23456789};
	SaveOptions options;
	options.precision = 3;
	std::ostringstream shortText;
	REQUIRE(Save(shortText, original, Format::Text, options).success);
	REQUIRE(shortText.str() == "1\n1.23\n");
	std::ostringstream json;
	REQUIRE(Save(json, original, Format::JSON, options).success);
	const auto parsed = ParseJson(json.str());
	REQUIRE(parsed.result.success);
	REQUIRE(parsed.value.object_value.at("data").array_value[0].number_value == 1.23);
	options.precision = -1;
	std::ostringstream invalid;
	REQUIRE(Save(invalid, original, Format::Text, options).error == IOError::INVALID_PARAMETERS);
	REQUIRE(invalid.str().empty());
	TempFile file;
	{ std::ofstream sentinel(file.path); sentinel << "unchanged"; }
	REQUIRE_FALSE(Save(original, file.path.string(), options).success);
	REQUIRE(Save(Vector<int>{1}, file.path.string(), Format::JSON).error == IOError::UNSUPPORTED_SCALAR);
	std::ifstream contents(file.path);
	std::string sentinel;
	contents >> sentinel;
	REQUIRE(sentinel == "unchanged");
	std::ostringstream csv;
	REQUIRE(Save(csv, original, Format::CSV).error == IOError::INVALID_FORMAT);
	REQUIRE(Save(csv, original, Format::Auto).error == IOError::INVALID_FORMAT);
	REQUIRE(Save(csv, Vector<std::string>{"a"}, Format::Text).error == IOError::UNSUPPORTED_SCALAR);
	REQUIRE(Save(csv, Vector<int>{1}, Format::Binary).error == IOError::UNSUPPORTED_SCALAR);
	REQUIRE(Save(csv, Vector<Complex>{Complex(1, 2)}, Format::JSON).error == IOError::UNSUPPORTED_SCALAR);
	std::istringstream unsupported("not consumed");
	Vector<int> integers{9};
	REQUIRE(Load(unsupported, integers, Format::JSON).error == IOError::UNSUPPORTED_SCALAR);
	REQUIRE(integers == Vector<int>{9});
	REQUIRE(unsupported.tellg() == 0);
}

TEST_CASE("Generic vector loading is transactional for JSON and binary and forwards options", "[VectorPersistence][Binary][JSON]")
{
	const Vector<double> original{1, 2, 3};
	for (Format format : {Format::JSON, Format::Binary}) {
		std::stringstream encoded(std::ios::in | std::ios::out | std::ios::binary);
		REQUIRE(Save(encoded, original, format).success);
		Vector<double> destination{9, 8};
		LoadOptions limit;
		limit.max_allocation_bytes = sizeof(double);
		const auto rejected = Load(encoded, destination, format, limit);
		REQUIRE(rejected.error == IOError::ALLOCATION_LIMIT_EXCEEDED);
		REQUIRE(destination == Vector<double>{9, 8});
		encoded.clear();
		encoded.seekg(0);
		REQUIRE(Load(encoded, destination, format).success);
		CheckVector(destination, original);
		TempFile file(format == Format::JSON ? ".mmlj" : ".mmlb");
		SaveOptions options;
		options.metadata.title = "vector";
		REQUIRE(Save(original, file.path.string(), options).success);
		Vector<double> unchanged{9};
		REQUIRE(Load(file.path.string(), unchanged, limit).error == IOError::ALLOCATION_LIMIT_EXCEEDED);
		REQUIRE(unchanged == Vector<double>{9});
		REQUIRE(Load(file.path.string(), unchanged).success);
		CheckVector(unchanged, original);
	}
	std::stringstream binary(std::ios::in | std::ios::out | std::ios::binary);
	REQUIRE(Save(binary, original, Format::Binary).success);
	const auto bytes = binary.str();
	std::istringstream truncated(bytes.substr(0, bytes.size() - 1), std::ios::binary);
	Vector<double> destination{9};
	REQUIRE_FALSE(Load(truncated, destination, Format::Binary).success);
	REQUIRE(destination == Vector<double>{9});
	std::stringstream encoded;
	REQUIRE(Save(encoded, original, Format::JSON).success);
	std::string invalid = encoded.str();
	const auto shape = invalid.find("\"shape\"");
	const auto opening = invalid.find('[', shape);
	const auto closing = invalid.find(']', opening);
	for (const char* count : {"1.5", "-1", "1e100"}) {
		std::string malformed = invalid;
		malformed.replace(opening + 1, closing - opening - 1, count);
		std::istringstream input(malformed);
		REQUIRE_FALSE(Load(input, destination, Format::JSON).success);
		REQUIRE(destination == Vector<double>{9});
	}
}

TEST_CASE("Vector files support auto extensions and explicit formats without hidden fallbacks", "[VectorPersistence][Files]")
{
	for (const char* extension : {".txt", ".dat", ".mml", ".TXT"}) {
		TempFile file(extension);
		REQUIRE(Save(Vector<int>{1, 2}, file.path.string()).success);
		Vector<int> destination;
		REQUIRE(Load(file.path.string(), destination).success);
		REQUIRE(destination == Vector<int>{1, 2});
		LoadOptions limit;
		limit.max_allocation_bytes = sizeof(int);
		REQUIRE(Load(file.path.string(), destination, limit).error == IOError::ALLOCATION_LIMIT_EXCEEDED);
		REQUIRE(destination == Vector<int>{1, 2});
	}
	TempFile unknown(".unknown");
	REQUIRE(Save(Vector<int>{1}, unknown.path.string()).error == IOError::INVALID_FORMAT);
	REQUIRE_FALSE(std::filesystem::exists(unknown.path));
	REQUIRE(Save(Vector<int>{1}, unknown.path.string(), Format::Text).success);
	Vector<int> restored{9};
	REQUIRE(Load(unknown.path.string(), restored, Format::Text).success);
	REQUIRE(restored == Vector<int>{1});
	REQUIRE_FALSE(Load((unknown.path.string() + ".missing"), restored, Format::Text).success);
	REQUIRE(restored == Vector<int>{1});
	TempFile malformed;
	{ std::ofstream contents(malformed.path); contents << "2\n1\nbad\n"; }
	const auto failure = Load(malformed.path.string(), restored);
	REQUIRE_FALSE(failure.success);
	REQUIRE(failure.filename == malformed.path.string());
	REQUIRE(failure.offset == 3);
	REQUIRE(restored == Vector<int>{1});
	REQUIRE(Save(Vector<int>{1}, std::string("bad\0file.txt", 12)).error == IOError::INVALID_PARAMETERS);
}

TEST_CASE("Large vector text transport preserves every value exactly", "[VectorPersistence][Text]")
{
	Vector<Real> original(1000);
	for (int i = 0; i < original.size(); ++i)
		original[i] = Real(i) + Real(1) / 7;
	std::stringstream stream;
	REQUIRE(Save(stream, original, Format::Text).success);
	Vector<Real> restored;
	REQUIRE(Load(stream, restored, Format::Text).success);
	CheckVector(restored, original);
}
