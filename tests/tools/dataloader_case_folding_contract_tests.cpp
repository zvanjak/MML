#include <catch2/catch_all.hpp>

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/tools/DataLoader.h>
#endif

#include <clocale>
#include <chrono>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using namespace MML;
using namespace MML::Data;

namespace {
	std::string OracleFold(std::string value) {
		const std::string uppercase = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
		const std::string lowercase = "abcdefghijklmnopqrstuvwxyz";
		for (char& byte : value) {
			const auto position = uppercase.find(byte);
			if (position != std::string::npos)
				byte = lowercase[position];
		}
		return value;
	}

	std::vector<std::string> CaseVariants(const std::string& token) {
		std::vector<std::string> variants{token};
		for (size_t index = 0; index < token.size(); ++index) {
			if (token[index] < 'a' || token[index] > 'z')
				continue;
			const size_t count = variants.size();
			for (size_t variant = 0; variant < count; ++variant) {
				std::string upper = variants[variant];
				upper[index] = static_cast<char>(token[index] - 'a' + 'A');
				variants.push_back(upper);
			}
		}
		return variants;
	}

	void CheckVocabulary() {
		for (const std::string sentinel : {"na", "n/a", "nan", "null", "none", ".", "?", "-", ""})
			for (const auto& token : CaseVariants(sentinel)) {
				CAPTURE(token);
				CHECK(IsMissingValue(token));
				CHECK(InferColumnType({token, "NaN"}) == ColumnType::STRING);
				Real real = 17;
				int integer = 19;
				bool boolean = true;
				std::string text = "original";
				CHECK(ParseValue(token, ColumnType::BOOL, real, integer, boolean, text) == ParseResult::Missing);
				CHECK(boolean);
				CHECK(text == "original");
			}
		for (const std::string spelling : {"true", "false", "yes", "no", "t", "f", "y", "n", "1", "0"}) {
			const bool expected = spelling == "true" || spelling == "yes" ||
				spelling == "t" || spelling == "y" || spelling == "1";
			for (const auto& token : CaseVariants(spelling)) {
				CAPTURE(token);
				Real real = 17;
				int integer = 19;
				bool boolean = !expected;
				std::string text = "original";
				CHECK_FALSE(IsMissingValue(token));
				CHECK(ParseValue(" \t" + token + "\r\n", ColumnType::BOOL,
					real, integer, boolean, text) == ParseResult::Parsed);
				CHECK(boolean == expected);
				CHECK(InferColumnType({token}) ==
					(spelling == "1" || spelling == "0" ? ColumnType::INT : ColumnType::BOOL));
			}
		}
		for (const std::string extension : {"csv", "tsv", "tab", "json"})
			for (const auto& token : CaseVariants(extension)) {
				const auto expected = extension == "json" ? DataFormat::JSON :
					extension == "csv" ? DataFormat::CSV : DataFormat::TSV;
				CHECK(DetectFormat("\xCE\x94" ".data." + token) == expected);
			}
		CHECK(DetectFormat("no_extension") == DataFormat::CSV);
		CHECK(DetectFormat("file.unknown") == DataFormat::CSV);
	}

	class ScopedCTypeLocale {
		std::string _original;
	public:
		ScopedCTypeLocale() {
			const char* name = std::setlocale(LC_CTYPE, nullptr);
			if (name == nullptr)
				throw DataError("Cannot query C character locale for regression test");
			_original = name;
		}
		~ScopedCTypeLocale() {
			if (std::setlocale(LC_CTYPE, _original.c_str()) == nullptr) {
				std::cerr << "Cannot restore C character locale after regression test\n";
				std::terminate();
			}
		}
	};

	struct InputFile {
		std::filesystem::path directory = std::filesystem::temp_directory_path() /
			("mml_case_folding_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
		std::string filename;
		InputFile(const std::string& content, DataFormat format) {
			std::filesystem::create_directory(directory);
			filename = (directory / (format == DataFormat::JSON ? "data.JsOn" :
				format == DataFormat::TSV ? "data.TsV" : "data.CsV")).string();
			std::ofstream stream(filename, std::ios::binary);
			stream << content;
			stream.close();
			REQUIRE(stream.good());
		}
		~InputFile() {
			std::error_code error;
			std::filesystem::remove_all(directory, error);
			if (error)
				std::cerr << "Cannot clean case-folding test input: " << error.message() << '\n';
		}
	};
}

TEST_CASE("DataLoader preserves every byte while folding only ASCII protocol vocabulary",
	"[dataloader][case-folding-contract]") {
	for (unsigned byte = 0; byte <= std::numeric_limits<unsigned char>::max(); ++byte) {
		CAPTURE(byte);
		const std::string token = "X" + std::string(1, static_cast<char>(byte)) + "Y";
		CHECK(static_cast<unsigned char>(token[1]) == byte);
		const auto folded = TextIO::ToLowerAscii(token[1]);
		CHECK(std::string(1, folded) == OracleFold(std::string(1, token[1])));
		CHECK(Persistence::ToLowerAscii(token[1]) == folded);
		CHECK_FALSE(IsMissingValue(token));
		CHECK(InferColumnType({token}) == ColumnType::STRING);
		Real real = 17;
		int integer = 19;
		bool boolean = true;
		std::string text;
		CHECK(ParseValue(token, ColumnType::STRING, real, integer, boolean, text) == ParseResult::Parsed);
		CHECK(text == token);
		CHECK(ParseValue(token, ColumnType::BOOL, real, integer, boolean, text) == ParseResult::Invalid);
		CHECK(boolean);
		const std::string extension = "ta" + std::string(1, static_cast<char>(byte));
		CHECK(DetectFormat("file." + extension) ==
			(OracleFold(extension) == "tab" ? DataFormat::TSV : DataFormat::CSV));
	}
}

TEST_CASE("DataLoader ASCII missing, boolean and extension vocabularies remain case insensitive",
	"[dataloader][case-folding-contract]") {
	CheckVocabulary();
}

TEST_CASE("DataLoader token recognition is independent of available C character locales",
	"[dataloader][case-folding-contract][locale]") {
	ScopedCTypeLocale restore;
	int available = 0;
	for (const char* name : {"C", ".UTF-8", "Turkish_Turkey.1254", "German_Germany.1252",
		"C.UTF-8", "tr_TR.UTF-8", "de_DE.UTF-8"}) {
		if (std::setlocale(LC_CTYPE, name) == nullptr)
			continue;
		++available;
		std::cout << "DataLoader case-folding locale: " << name << '\n';
		CAPTURE(name);
		CheckVocabulary();
		const std::string nonAscii = "N\xC4\xB0" "L";
		CHECK_FALSE(IsMissingValue(nonAscii));
		CHECK(InferColumnType({nonAscii}) == ColumnType::STRING);
	}
	REQUIRE(available >= 1);
}

TEST_CASE("DataLoader Unicode text is not Unicode case folded or treated as an ASCII sentinel",
	"[dataloader][case-folding-contract]") {
	const std::vector<std::string> values{
		"Caf\xC3\xA9", "\xCE\xA3", "\xCF\x83", "\xCF\x82", "\xE6\xBC\xA2\xE5\xAD\x97",
		"N\xC3\x84" "N", "N\xC4\xB0" "L", "\xF0\x9F\x9A\x80"};
	for (const auto& value : values) {
		CAPTURE(value);
		CHECK_FALSE(IsMissingValue(value));
		CHECK(InferColumnType({value}) == ColumnType::STRING);
		Real real = 17;
		int integer = 19;
		bool boolean = true;
		std::string text;
		CHECK(ParseValue(value, ColumnType::STRING, real, integer, boolean, text) == ParseResult::Parsed);
		CHECK(text == value);
		CHECK(ParseValue(value, ColumnType::BOOL, real, integer, boolean, text) == ParseResult::Invalid);
		CHECK(boolean);
		for (const auto& surrounding : {value + "TRUE", "NULL" + value, value + "NaN"})
			CHECK_FALSE(IsMissingValue(surrounding));
	}
}

TEST_CASE("DataLoader CSV TSV and JSON preserve decoded Unicode strings and headers",
	"[dataloader][case-folding-contract]") {
	const std::string header = "\xCE\x94";
	const std::vector<std::string> values{"Caf\xC3\xA9", "\xCE\xA3", "\xCF\x83",
		"N\xC3\x84" "N", "\xF0\x9F\x9A\x80"};
	for (auto format : {DataFormat::CSV, DataFormat::TSV, DataFormat::JSON}) {
		std::string content = format == DataFormat::JSON ? "[" : header + "\n";
		for (size_t row = 0; row < values.size(); ++row) {
			if (format == DataFormat::JSON) {
				if (row != 0) content += ",";
				content += "{\"" + header + "\":\"" + values[row] + "\"}";
			}
			else
				content += values[row] + "\n";
		}
		if (format == DataFormat::JSON) content += "]";
		const auto data = LoadFromString(content, format);
		REQUIRE(data.NumRows() == values.size());
		REQUIRE(data.NumColumns() == 1);
		CHECK(data[0].name == header);
		CHECK(data[0].type == ColumnType::STRING);
		CHECK(data[0].stringData == values);
		CHECK(data[0].CountMissing() == 0);
	}
	const auto escaped = LoadFromJSONString(
		R"([{"\u0394":"Caf\u00e9"},{"\u0394":"\u03a3"},{"\u0394":"\u03c3"},{"\u0394":"N\u00c4N"},{"\u0394":"\ud83d\ude80"}])");
	CHECK(escaped[0].name == header);
	CHECK(escaped[0].stringData == values);
}

TEST_CASE("DataLoader file dispatch preserves Unicode and rejects Unicode boolean minorities transactionally",
	"[dataloader][case-folding-contract]") {
	const std::string value = "N\xC3\x84" "N";
	for (auto format : {DataFormat::CSV, DataFormat::TSV, DataFormat::JSON}) {
		CAPTURE(format);
		const std::string valid = format == DataFormat::JSON ?
			"[{\"value\":\"" + value + "\"}]" : "value\n" + value + "\n";
		InputFile validInput(valid, format);
		Dataset loaded;
		const auto success = LoadData(validInput.filename, loaded);
		INFO(success.message);
		REQUIRE(success.success);
		CHECK(loaded["value"].stringData[0] == value);
		CHECK_FALSE(loaded["value"].IsMissing(0));
		std::string invalid = format == DataFormat::JSON ? "[" : "value\n";
		for (int row = 0; row < 9; ++row)
			invalid += format == DataFormat::JSON ? "{\"value\":\"TrUe\"}," : "TrUe\n";
		invalid += format == DataFormat::JSON ? "{\"value\":\"" + value + "\"}]" : value + "\n";
		InputFile invalidInput(invalid, format);
		Dataset original = LoadFromCSVString("sentinel\n42\n");
		original.name = "original";
		const auto rejected = LoadData(invalidInput.filename, original);
		CHECK_FALSE(rejected.success);
		CHECK(rejected.error == IOError::MALFORMED_INPUT);
		CHECK(rejected.message.find(value) != std::string::npos);
		CHECK(rejected.message.find("data row 10") != std::string::npos);
		CHECK(original.name == "original");
		CHECK(original.NumRows() == 1);
		CHECK(original.NumColumns() == 1);
		CHECK(original[0].name == "sentinel");
		CHECK(original[0].intData[0] == 42);
		CHECK_THROWS_AS(LoadData(invalidInput.filename), DataError);
	}
}
