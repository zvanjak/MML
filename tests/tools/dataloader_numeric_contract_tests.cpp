#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/tools/DataLoader.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <sstream>
#include <string>

using namespace MML;
using namespace MML::Data;
using Catch::Matchers::ContainsSubstring;

namespace {
	struct NumericInput {
		std::filesystem::path directory = std::filesystem::temp_directory_path() /
			("mml_numeric_loader_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
		std::string filename;

		NumericInput(const std::string& content, DataFormat format) {
			std::filesystem::create_directory(directory);
			filename = (directory / (format == DataFormat::JSON ? "data.json" :
				format == DataFormat::TSV ? "data.tsv" : "data.csv")).string();
			std::ofstream file(filename, std::ios::binary);
			file << content;
			file.close();
			REQUIRE(file.good());
		}
		~NumericInput() {
			std::error_code error;
			std::filesystem::remove_all(directory, error);
			if (error) std::cerr << "Failed to clean numeric loader input: " << error.message() << '\n';
		}
	};

	std::string MajorityInput(const std::string& tail, DataFormat format) {
		std::string content = format == DataFormat::JSON ? "[" : "value\n";
		for (int i = 1; i <= 9; ++i) {
			content += format == DataFormat::JSON ? "{\"value\":" + std::to_string(i) + "}," :
				std::to_string(i) + "\n";
		}
		return content + (format == DataFormat::JSON ? "{\"value\":" + tail + "}]" : tail + "\n");
	}

	std::string ExactIntegerBoundary() {
		std::ostringstream text;
		TextIO::ApplyClassicFormat(text);
		text << std::fixed << std::setprecision(0) << std::ldexp(Real(1),
			std::max(std::numeric_limits<Real>::digits, std::numeric_limits<int>::digits));
		return text.str();
	}

	std::string NextInteger(std::string value) {
		for (auto i = value.size(); i > 0; --i) {
			if (value[i - 1] != '9') {
				++value[i - 1];
				return value;
			}
			value[i - 1] = '0';
		}
		return "1" + value;
	}

	void CheckUnchanged(const Dataset& data) {
		CHECK(data.name == "original");
		REQUIRE(data.NumColumns() == 1);
		CHECK(data.NumRows() == 1);
		CHECK(data[0].name == "sentinel");
		REQUIRE(data[0].type == ColumnType::INT);
		CHECK(data[0].intData[0] == 42);
		CHECK_FALSE(data[0].IsMissing(0));
	}

	Dataset OriginalOutput() {
		Dataset data = LoadFromCSVString("sentinel\n42\n");
		data.name = "original";
		return data;
	}
}

TEST_CASE("Numeric loader promotes fractional and wide integer minorities", "[dataloader][numeric-contract]") {
	std::vector<std::string> values{"1", "2", "3", "4", "5", "6", "7", "8", "9"};
	for (const std::string tail : {"10.75", "1e2", "2147483648"}) {
		auto column = values;
		column.push_back(tail);
		CAPTURE(tail);
		CHECK(InferColumnType(column) == ColumnType::REAL);
	}
	CHECK(InferColumnType(values) == ColumnType::INT);
	CHECK(InferColumnType({"2147483648"}) == ColumnType::REAL);
}

TEST_CASE("Numeric loader integer parsing consumes the entire token transactionally", "[dataloader][numeric-contract]") {
	Real real = 7;
	int integer = 17;
	bool boolean = false;
	std::string text;
	for (const std::string value : {"10.75", "10oops", "1e2", "++1", "+-1", "2147483648", "-2147483649"}) {
		CAPTURE(value);
		CHECK(ParseValue(value, ColumnType::INT, real, integer, boolean, text) == ParseResult::Invalid);
		CHECK(integer == 17);
		integer = 17;
	}
	CHECK(ParseValue(" +00042 ", ColumnType::INT, real, integer, boolean, text) == ParseResult::Parsed);
	CHECK(integer == 42);
	CHECK(ParseValue(std::to_string(std::numeric_limits<int>::min()), ColumnType::INT,
		real, integer, boolean, text) == ParseResult::Parsed);
	CHECK(integer == std::numeric_limits<int>::min());
	CHECK(ParseValue(std::to_string(std::numeric_limits<int>::max()), ColumnType::INT,
		real, integer, boolean, text) == ParseResult::Parsed);
	CHECK(integer == std::numeric_limits<int>::max());
}

TEST_CASE("Numeric loader rejects inexact integer tokens in Real storage", "[dataloader][numeric-contract]") {
	Real real = 7;
	int integer = 0;
	bool boolean = false;
	std::string text;
	const auto exact = ExactIntegerBoundary();
	const auto inexact = NextInteger(exact);
	for (const auto& value : {inexact, "-" + inexact, "+000" + inexact}) {
		CAPTURE(value);
		CHECK(ParseValue(value, ColumnType::REAL, real, integer, boolean, text) == ParseResult::Invalid);
		CHECK(real == Real(7));
		real = 7;
	}
	for (const auto& value : {exact, "+000" + exact, "-" + exact, std::string("-000")}) {
		CAPTURE(value);
		CHECK(ParseValue(value, ColumnType::REAL, real, integer, boolean, text) == ParseResult::Parsed);
	}
	CHECK(real == Real(0));
}

TEST_CASE("Numeric loader preserves whole-token Real parsing and rejects range loss", "[dataloader][numeric-contract]") {
	Real real = 7;
	int integer = 0;
	bool boolean = false;
	std::string text;
	for (const std::string value : {"12.75oops", "1e99999", "1e-99999", "-1e-99999", ".01e-99999"}) {
		CAPTURE(value);
		CHECK(ParseValue(value, ColumnType::REAL, real, integer, boolean, text) == ParseResult::Invalid);
		CHECK(real == Real(7));
		real = 7;
	}
	for (const std::string value : {"0.1", "10.75", "1e2", "-0.0", "0e-20", "0.000e-99999"}) {
		CAPTURE(value);
		CHECK(ParseValue(value, ColumnType::REAL, real, integer, boolean, text) == ParseResult::Parsed);
	}
	CHECK(real == Real(0));
	CHECK(ParseValue("-0e-99999", ColumnType::REAL, real, integer, boolean, text) == ParseResult::Parsed);
	CHECK(std::signbit(real));
}

TEST_CASE("Numeric loader string entrypoints preserve fractional minorities and wide integers", "[dataloader][numeric-contract]") {
	for (auto format : {DataFormat::CSV, DataFormat::TSV, DataFormat::JSON}) {
		for (const std::string tail : {"10.75", "2147483648"}) {
			CAPTURE(format, tail);
			Dataset data = LoadFromString(MajorityInput(tail, format), format);
			REQUIRE(data["value"].type == ColumnType::REAL);
			CHECK(data.NumRows() == 10);
			CHECK(data["value"].realData[9] == (tail == "10.75" ? Real(10.75) : Real(2147483648.0)));
			CHECK(data["value"].CountMissing() == 0);
		}
	}
}

TEST_CASE("Numeric loader file entrypoints preserve fractional minorities and wide integers", "[dataloader][numeric-contract]") {
	for (auto format : {DataFormat::CSV, DataFormat::TSV, DataFormat::JSON}) {
		for (const std::string tail : {"10.75", "2147483648"}) {
			CAPTURE(format, tail);
			NumericInput input(MajorityInput(tail, format), format);
			Dataset data;
			const auto result = LoadData(input.filename, data);
			INFO(result.message);
			REQUIRE(result.success);
			REQUIRE(data["value"].type == ColumnType::REAL);
			CHECK(data["value"].realData[9] == (tail == "10.75" ? Real(10.75) : Real(2147483648.0)));
			CHECK(data["value"].CountMissing() == 0);
		}
	}
}

TEST_CASE("Numeric loader string entrypoints reject junk with data row and column diagnostics", "[dataloader][numeric-contract]") {
	for (auto format : {DataFormat::CSV, DataFormat::TSV, DataFormat::JSON}) {
		CAPTURE(format);
		const auto content = MajorityInput(format == DataFormat::JSON ? "\"10oops\"" : "10oops", format);
		CHECK_THROWS_WITH(LoadFromString(content, format),
			ContainsSubstring("10oops") && ContainsSubstring("data row 10") && ContainsSubstring("column 1"));
		CHECK_THROWS_AS(LoadFromString(content, format), DataError);
	}
}

TEST_CASE("Numeric loader result entrypoints reject junk without changing output", "[dataloader][numeric-contract]") {
	for (auto format : {DataFormat::CSV, DataFormat::TSV, DataFormat::JSON}) {
		CAPTURE(format);
		NumericInput input(MajorityInput(format == DataFormat::JSON ? "\"10oops\"" : "10oops", format), format);
		Dataset data = OriginalOutput();
		const auto result = format == DataFormat::JSON ? LoadJSONSafe(input.filename, data) :
			format == DataFormat::TSV ? LoadTSVSafe(input.filename, data) : LoadCSVSafe(input.filename, data);
		CHECK_FALSE(result.success);
		CHECK(result.error == IOError::MALFORMED_INPUT);
		CHECK(result.filename == input.filename);
		CHECK_THAT(result.message, ContainsSubstring("data row 10") && ContainsSubstring("column 1") &&
			ContainsSubstring("value") && ContainsSubstring("10oops"));
		CheckUnchanged(data);
		CHECK_THROWS_AS(LoadData(input.filename), DataError);
		CHECK_FALSE(LoadData(input.filename, data).success);
		CheckUnchanged(data);
	}
}

TEST_CASE("Numeric loader CSV diagnostics retain physical lines and column identity", "[dataloader][numeric-contract]") {
	std::string content = "\n\nid,value\n";
	for (int i = 1; i <= 9; ++i)
		content += std::to_string(i) + "," + std::to_string(i) + "\n";
	content += "10,10oops\n";
	NumericInput input(content, DataFormat::CSV);
	Dataset data = OriginalOutput();
	const auto result = LoadCSV(input.filename, data);
	CHECK_FALSE(result.success);
	CHECK(result.offset_kind == IOOffsetKind::Line);
	CHECK(result.offset == 13);
	CHECK_THAT(result.message, ContainsSubstring("data row 10") && ContainsSubstring("column 2") &&
		ContainsSubstring("value"));
	CheckUnchanged(data);
}

TEST_CASE("Numeric loader rejects inexact integer tokens through every data format", "[dataloader][numeric-contract]") {
	const auto token = NextInteger(ExactIntegerBoundary());
	for (auto format : {DataFormat::CSV, DataFormat::TSV, DataFormat::JSON}) {
		CAPTURE(format, token);
		const auto content = MajorityInput(token, format);
		CHECK_THROWS_WITH(LoadFromString(content, format), ContainsSubstring(token));
		NumericInput input(content, format);
		Dataset data = OriginalOutput();
		const auto result = LoadData(input.filename, data);
		CHECK_FALSE(result.success);
		CHECK(result.error == IOError::MALFORMED_INPUT);
		CHECK_THAT(result.message, ContainsSubstring(token));
		CheckUnchanged(data);
	}
}

TEST_CASE("Numeric loader keeps recognized missing cells distinct from invalid values", "[dataloader][numeric-contract]") {
	for (auto format : {DataFormat::CSV, DataFormat::TSV, DataFormat::JSON}) {
		Dataset data = LoadFromString(MajorityInput(format == DataFormat::JSON ? "null" : "NA", format), format);
		REQUIRE(data["value"].type == ColumnType::INT);
		CHECK(data["value"].CountMissing() == 1);
		CHECK(data["value"].IsMissing(9));
		CHECK(data["value"].intData[8] == 9);
	}
}

TEST_CASE("Numeric loader explicit string storage preserves otherwise invalid numeric tokens", "[dataloader][numeric-contract]") {
	NumericInput input(MajorityInput("10oops", DataFormat::CSV), DataFormat::CSV);
	CSVLoadOptions options;
	options.inferTypes = false;
	Dataset data;
	REQUIRE(LoadCSV(input.filename, data, options).success);
	REQUIRE(data["value"].type == ColumnType::STRING);
	CHECK(data["value"].stringData[9] == "10oops");
	CHECK(data["value"].CountMissing() == 0);
}
