#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/tools/DataLoader.h>
#endif

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace MML;
using namespace MML::Data;
using Catch::Matchers::ContainsSubstring;

namespace {
	struct CSVRecordInput {
		std::filesystem::path directory = std::filesystem::temp_directory_path() /
			("mml_csv_records_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
		std::string filename;

		CSVRecordInput(const std::string& content, char delimiter) {
			std::filesystem::create_directory(directory);
			filename = (directory / (delimiter == '\t' ? "data.tsv" : "data.csv")).string();
			std::ofstream file(filename, std::ios::binary);
			file << content;
			file.close();
			REQUIRE(file.good());
		}
		~CSVRecordInput() {
			std::error_code error;
			std::filesystem::remove_all(directory, error);
			if (error) std::cerr << "Failed to clean CSV record input: " << error.message() << '\n';
		}
	};

	std::string MultilineInput(char delimiter, const std::string& newline) {
		return "\xEF\xBB\xBF\"label\"" + std::string(1, delimiter) + "value" + delimiter + "empty" + newline +
			newline + "\"hello" + newline + newline + "world\"" + delimiter + "1" + delimiter + newline +
			"\"a" + delimiter + "b \"\"quoted\"\"\"" + delimiter + "2" + delimiter;
	}

	void CheckMultiline(const Dataset& data, char delimiter, const std::string& newline) {
		REQUIRE(data.NumRows() == 2);
		REQUIRE(data.NumColumns() == 3);
		REQUIRE(data["label"].type == ColumnType::STRING);
		CHECK(data["label"].stringData[0] == "hello" + newline + newline + "world");
		CHECK(data["label"].stringData[1] == "a" + std::string(1, delimiter) + "b \"quoted\"");
		REQUIRE(data["value"].type == ColumnType::INT);
		CHECK(data["value"].intData[0] == 1);
		CHECK(data["value"].intData[1] == 2);
		CHECK(data["empty"].IsMissing(0));
		CHECK(data["empty"].IsMissing(1));
	}

	Dataset OriginalCSVOutput() {
		Dataset data = LoadFromCSVString("sentinel\n42\n");
		data.name = "original";
		return data;
	}

	void CheckCSVUnchanged(const Dataset& data) {
		CHECK(data.name == "original");
		REQUIRE(data.NumColumns() == 1);
		REQUIRE(data.NumRows() == 1);
		CHECK(data[0].name == "sentinel");
		REQUIRE(data[0].type == ColumnType::INT);
		CHECK(data[0].intData[0] == 42);
		CHECK_FALSE(data[0].IsMissing(0));
	}

	struct MalformedRecord {
		std::string text;
		std::int64_t line;
	};

	const std::vector<MalformedRecord> malformedRecords{
		{"\"unterminated,1", 3},
		{"bad\"quote,1\n", 3},
		{"\"closed\"x,1\n", 3},
		{" \"padded\",1\n", 3},
		{"\"closed\" ,1\n", 3},
		{"\"first\nsecond\"x,1\n", 4},
		{"\"first\nsecond,1", 3}
	};

	class FailingCSVBuffer : public std::stringbuf {
	public:
		explicit FailingCSVBuffer(const std::string& prefix) : std::stringbuf(prefix, std::ios::in) {}
	protected:
		int_type underflow() override {
			throw std::ios_base::failure("Injected CSV stream failure");
		}
	};
}

TEST_CASE("CSV record string loaders preserve multiline fields and newline bytes", "[dataloader][csv-record]") {
	for (char delimiter : {',', '\t'}) {
		for (const std::string newline : {"\n", "\r\n"}) {
			CAPTURE(delimiter, newline);
			const auto content = MultilineInput(delimiter, newline);
			CheckMultiline(LoadFromCSVString(content, delimiter), delimiter, newline);
			CheckMultiline(LoadFromString(content, delimiter == '\t' ? DataFormat::TSV : DataFormat::CSV),
				delimiter, newline);
		}
	}
}

TEST_CASE("CSV record file loaders share multiline parsing with string loaders", "[dataloader][csv-record]") {
	for (char delimiter : {',', '\t'}) {
		for (const std::string newline : {"\n", "\r\n"}) {
			CAPTURE(delimiter, newline);
			CSVRecordInput input(MultilineInput(delimiter, newline), delimiter);
			CSVLoadOptions opts;
			opts.delimiter = delimiter;
			opts.strictSchema = true;
			Dataset output;
			REQUIRE(LoadCSV(input.filename, output, opts));
			CheckMultiline(output, delimiter, newline);
			REQUIRE(LoadCSVSafe(input.filename, output, true, delimiter));
			CheckMultiline(output, delimiter, newline);
			REQUIRE(LoadData(input.filename, output));
			CheckMultiline(output, delimiter, newline);
			CheckMultiline(LoadData(input.filename), delimiter, newline);
			CheckMultiline(LoadCSV(input.filename, opts), delimiter, newline);
			if (delimiter == '\t') {
				REQUIRE(LoadTSVSafe(input.filename, output));
				CheckMultiline(output, delimiter, newline);
				CheckMultiline(LoadTSV(input.filename), delimiter, newline);
			}
		}
	}
}

TEST_CASE("CSV records roundtrip escaped fields including final doubled quotes", "[dataloader][csv-record]") {
	for (char delimiter : {',', '\t', ';'}) {
		for (const std::string value : {"a\nb", "a\r\nb", "a\rb", "a,b", "a\tb", "a;b", "\"", "ends\""}) {
			CAPTURE(delimiter, value);
			const std::string content = "text\n" + CsvUtils::Escape(value, delimiter);
			auto data = LoadFromCSVString(content, delimiter);
			REQUIRE(data.NumRows() == 1);
			REQUIRE(data["text"].type == ColumnType::STRING);
			CHECK(data["text"].stringData[0] == value);
			CSVRecordInput input(content, delimiter);
			auto fromFile = LoadCSV(input.filename, true, delimiter);
			REQUIRE(fromFile.NumRows() == 1);
			CHECK(fromFile["text"].stringData[0] == value);
		}
	}
}

TEST_CASE("CSV records keep headerless multiline rows and final empty fields", "[dataloader][csv-record]") {
	for (char delimiter : {',', '\t', ';'}) {
		const std::string content = "\n \t \n\"first\nsecond\"" + std::string(1, delimiter) + "42" + delimiter;
		auto data = LoadFromCSVString(content, delimiter, false);
		REQUIRE(data.NumRows() == 1);
		REQUIRE(data.NumColumns() == 3);
		CHECK(data["Column0"].stringData[0] == "first\nsecond");
		CHECK(data["Column1"].intData[0] == 42);
		CHECK(data["Column2"].IsMissing(0));
		CSVRecordInput input(content, delimiter);
		auto fromFile = LoadCSV(input.filename, false, delimiter);
		REQUIRE(fromFile.NumRows() == 1);
		CHECK(fromFile["Column0"].stringData[0] == data["Column0"].stringData[0]);
		CHECK(fromFile["Column2"].IsMissing(0));
	}
}

TEST_CASE("CSV records preserve existing blank-record trimming and missing-marker rules", "[dataloader][csv-record]") {
	const std::string content = "\xEF\xBB\xBF\n \t \ntext\r\n\"  keep  \"\r\n\"NA\"\r\n\"\"\r\n \t \r\n";
	auto data = LoadFromCSVString(content);
	REQUIRE(data.NumRows() == 3);
	CHECK(data["text"].stringData[0] == "keep");
	CHECK(data["text"].IsMissing(1));
	CHECK(data["text"].IsMissing(2));
	CSVRecordInput input(content, ',');
	auto fromFile = LoadCSV(input.filename);
	REQUIRE(fromFile.NumRows() == 3);
	CHECK(fromFile["text"].stringData[0] == "keep");
	CHECK(fromFile["text"].IsMissing(1));
	CHECK(fromFile["text"].IsMissing(2));
	CHECK_THROWS_AS(LoadFromCSVString("\xEF\xBB\xBF\r\n \t \r\n"), DataError);
}

TEST_CASE("CSV record string loaders reject malformed quote syntax with physical lines", "[dataloader][csv-record]") {
	for (const auto& malformed : malformedRecords) {
		for (char delimiter : {',', '\t'}) {
			std::string record = malformed.text;
			std::replace(record.begin(), record.end(), ',', delimiter);
			const std::string content = "\ntext" + std::string(1, delimiter) + "value\n" + record;
			CAPTURE(record, delimiter);
			CHECK_THROWS_WITH(LoadFromCSVString(content, delimiter),
				ContainsSubstring("line " + std::to_string(malformed.line)));
			CHECK_THROWS_AS(LoadFromString(content, delimiter == '\t' ? DataFormat::TSV : DataFormat::CSV), DataError);
		}
	}
}

TEST_CASE("CSV record result APIs reject malformed quotes without replacing output", "[dataloader][csv-record]") {
	for (const auto& malformed : malformedRecords) {
		for (char delimiter : {',', '\t'}) {
			std::string record = malformed.text;
			std::replace(record.begin(), record.end(), ',', delimiter);
			CSVRecordInput input("\ntext" + std::string(1, delimiter) + "value\n" + record, delimiter);
			CAPTURE(record, delimiter);
			for (int route = 0; route < 3; ++route) {
				Dataset output = OriginalCSVOutput();
				CSVLoadOptions opts;
				opts.delimiter = delimiter;
				auto result = route == 0 ? LoadCSV(input.filename, output, opts) :
					route == 1 ? LoadCSVSafe(input.filename, output, true, delimiter) :
					LoadData(input.filename, output);
				CHECK_FALSE(result);
				CHECK(result.error == IOError::MALFORMED_INPUT);
				CHECK(result.filename == input.filename);
				CHECK(result.offset_kind == IOOffsetKind::Line);
				CHECK(result.offset == malformed.line);
				CHECK_THAT(result.message, ContainsSubstring("line " + std::to_string(malformed.line)));
				CheckCSVUnchanged(output);
			}
			CHECK_THROWS_AS(LoadCSV(input.filename, true, delimiter), DataError);
			CHECK_THROWS_AS(LoadData(input.filename), DataError);
		}
	}
}

TEST_CASE("CSV schema failures retain physical lines and logical warning rows after multiline records", "[dataloader][csv-record]") {
	for (char delimiter : {',', '\t'}) {
		const std::string content = "\n\ntext" + std::string(1, delimiter) + "value\n\"first\nsecond\"" +
			delimiter + "1\nwrong\n";
		CSVRecordInput input(content, delimiter);
		CSVLoadOptions opts;
		opts.delimiter = delimiter;
		opts.strictSchema = true;
		std::vector<std::string> warnings;
		opts.warnings = &warnings;
		Dataset output = OriginalCSVOutput();
		const auto result = LoadCSV(input.filename, output, opts);
		CHECK_FALSE(result);
		CHECK(result.error == IOError::SCHEMA_MISMATCH);
		CHECK(result.offset_kind == IOOffsetKind::Line);
		CHECK(result.offset == 6);
		CHECK_THAT(result.message, ContainsSubstring("Row 3 has 1 fields, expected 2"));
		CHECK(warnings.empty());
		CheckCSVUnchanged(output);
		opts.strictSchema = false;
		REQUIRE(LoadCSV(input.filename, output, opts));
		REQUIRE(warnings.size() == 1);
		CHECK(warnings[0] == "Row 3 has 1 fields, expected 2");
		REQUIRE(output.NumRows() == 2);
		CHECK(output["text"].stringData[0] == "first\nsecond");
		CHECK(output["value"].IsMissing(1));
	}
}

TEST_CASE("CSV numeric conversion failures retain multiline source locations", "[dataloader][csv-record]") {
	std::string content = "\ntext,value\n\"first\nsecond\",1\n";
	for (int i = 2; i <= 9; ++i) content += "plain," + std::to_string(i) + "\n";
	content += "bad,10oops\n";
	CSVRecordInput input(content, ',');
	Dataset output = OriginalCSVOutput();
	const auto result = LoadCSVSafe(input.filename, output);
	CHECK_FALSE(result);
	CHECK(result.error == IOError::MALFORMED_INPUT);
	CHECK(result.offset_kind == IOOffsetKind::Line);
	CHECK(result.offset == 13);
	CHECK_THAT(result.message, ContainsSubstring("data row 10, column 2 ('value')"));
	CheckCSVUnchanged(output);
	CHECK_THROWS_WITH(LoadFromCSVString(content),
		ContainsSubstring("data row 10, column 2 ('value')"));
}

TEST_CASE("CSV multiline parsing respects explicit string storage", "[dataloader][csv-record]") {
	CSVRecordInput input("text,value\n\"first\r\nsecond\",10oops\n", ',');
	CSVLoadOptions opts;
	opts.inferTypes = false;
	Dataset data;
	REQUIRE(LoadCSV(input.filename, data, opts));
	REQUIRE(data.NumRows() == 1);
	CHECK(data["text"].stringData[0] == "first\r\nsecond");
	REQUIRE(data["value"].type == ColumnType::STRING);
	CHECK(data["value"].stringData[0] == "10oops");
}

TEST_CASE("Shared CSV record reader publishes records and locations only on success", "[dataloader][csv-record]") {
	std::vector<std::vector<std::string>> records{{"original"}};
	std::vector<std::int64_t> lines{91};
	for (const auto& malformed : malformedRecords) {
		std::istringstream input("\ntext,value\n" + malformed.text);
		const auto result = CsvUtils::ReadRecords(input, records, lines);
		CHECK(result.error == IOError::MALFORMED_INPUT);
		CHECK(result.offset == malformed.line);
		CHECK(records == std::vector<std::vector<std::string>>{{"original"}});
		CHECK(lines == std::vector<std::int64_t>{91});
	}
	std::istringstream input("\n\"head\nname\",value\n\"first\nsecond\",1\nlast,2");
	REQUIRE(CsvUtils::ReadRecords(input, records, lines));
	CHECK(lines == std::vector<std::int64_t>{2, 4, 6});
	REQUIRE(records.size() == 3);
	CHECK(records[0] == std::vector<std::string>{"head\nname", "value"});
	CHECK(records[1] == std::vector<std::string>{"first\nsecond", "1"});
	CHECK(records[2] == std::vector<std::string>{"last", "2"});
}

TEST_CASE("Shared CSV record reader distinguishes I/O failure from EOF and malformed quotes", "[dataloader][csv-record]") {
	std::vector<std::vector<std::string>> records{{"original"}};
	std::vector<std::int64_t> lines{91};
	for (const auto flags : {std::ios::badbit, std::ios::failbit}) {
		std::istringstream input("text\nok\n");
		input.setstate(flags);
		const auto result = CsvUtils::ReadRecords(input, records, lines);
		CHECK(result.error == IOError::READ_FAILED);
		CHECK(result.offset_kind == IOOffsetKind::Line);
		CHECK(result.offset == 1);
	}
	for (const std::string prefix : {"text\nok\n", "text\nok\n\"unfinished\n"}) {
		FailingCSVBuffer buffer(prefix);
		std::istream input(&buffer);
		const auto result = CsvUtils::ReadRecords(input, records, lines);
		CHECK(result.error == IOError::READ_FAILED);
		CHECK(result.offset_kind == IOOffsetKind::Line);
		CHECK(result.offset == (prefix.back() == '\n' && prefix.find("unfinished") != std::string::npos ? 4 : 3));
		CHECK(records == std::vector<std::vector<std::string>>{{"original"}});
		CHECK(lines == std::vector<std::int64_t>{91});
	}
}

TEST_CASE("CSV record loaders reject ambiguous delimiters and preserve empty-input errors", "[dataloader][csv-record]") {
	CSVRecordInput input("text\nok\n", ',');
	for (char delimiter : {'"', '\n', '\r'}) {
		CSVLoadOptions opts;
		opts.delimiter = delimiter;
		Dataset output = OriginalCSVOutput();
		const auto result = LoadCSV(input.filename, output, opts);
		CHECK(result.error == IOError::INVALID_PARAMETERS);
		CheckCSVUnchanged(output);
		CHECK_THROWS_AS(LoadFromCSVString("text\nok\n", delimiter), DataError);
	}
	for (const std::string content : {"", "\n \t \r\n", "\xEF\xBB\xBF"}) {
		CSVRecordInput empty(content, ',');
		Dataset output = OriginalCSVOutput();
		const auto result = LoadCSVSafe(empty.filename, output);
		CHECK(result.error == IOError::INVALID_FORMAT);
		CHECK(result.offset == 1);
		CheckCSVUnchanged(output);
		CHECK_THROWS_AS(LoadFromCSVString(content), DataError);
	}
	auto headerOnly = LoadFromCSVString("text,value");
	CHECK(headerOnly.NumRows() == 0);
	CHECK(headerOnly.NumColumns() == 2);
}
