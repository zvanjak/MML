///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        CsvUtils.h                                                          ///
///  Description: Single RFC 4180 CSV utility module - escape, unescape, and          ///
///               quote-honoring field splitting shared by ConsolePrinter,            ///
///               DataLoader, and MatrixIO                                            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CSV_UTILS_H
#define MML_CSV_UTILS_H

#include <mml/tools/IOResult.h>

#include <cstdint>
#include <istream>
#include <string>
#include <utility>
#include <vector>

namespace MML::CsvUtils
{
	/// @brief Check if a field needs quoting per RFC 4180
	/// (contains the delimiter, a double quote, or a newline)
	inline bool NeedsQuoting(const std::string& str, char delimiter = ',')
	{
		return str.find_first_of(std::string{ delimiter, '"', '\n', '\r' }) != std::string::npos;
	}

	/// @brief Escape a field for RFC 4180 compliant CSV output
	/// (quote if necessary, double any embedded quotes)
	inline std::string Escape(const std::string& str, char delimiter = ',')
	{
		if (!NeedsQuoting(str, delimiter))
			return str;

		std::string result;
		result.reserve(str.size() + 2);
		result += '"';
		for (char c : str) {
			if (c == '"')
				result += "\"\"";
			else
				result += c;
		}
		result += '"';
		return result;
	}

	/// @brief Reverse of Escape: strip surrounding quotes, un-double embedded quotes
	inline std::string Unescape(const std::string& str)
	{
		if (str.empty())
			return str;
		if (str.front() != '"' || str.back() != '"' || str.size() < 2)
			return str;   // not quoted

		std::string result;
		result.reserve(str.size() - 2);
		for (size_t i = 1; i < str.size() - 1; ++i) {
			if (str[i] == '"' && i + 1 < str.size() - 1 && str[i + 1] == '"') {
				result += '"';
				++i;
			}
			else
				result += str[i];
		}
		return result;
	}

	/// @brief Legacy permissive field splitter, honoring quoted fields and doubled quotes.
	/// @note Does not validate quote syntax or read multiline records. Use ReadRecords
	///       for validated CSV/TSV input.
	inline std::vector<std::string> SplitLine(const std::string& line, char delimiter = ',')
	{
		std::vector<std::string> result;
		std::string field;
		bool inQuotes = false;

		for (size_t i = 0; i < line.size(); ++i) {
			char c = line[i];

			if (c == '"') {
				if (inQuotes && i + 1 < line.size() && line[i + 1] == '"') {
					field += '"';   // doubled quote inside quoted field
					++i;
				}
				else
					inQuotes = !inQuotes;
			}
			else if (c == delimiter && !inQuotes) {
				result.push_back(field);
				field.clear();
			}
			else
				field += c;
		}
		result.push_back(field);

		return result;
	}

	/// @brief Read strict CSV/TSV records and their one-based starting physical lines.
	/// @details Accepts an initial UTF-8 BOM and skips blank physical lines outside
	/// quoted fields. Preserves embedded LF/CRLF bytes; file streams must be binary.
	/// Both outputs remain unchanged on syntax or stream failure.
	inline IOResult ReadRecords(std::istream& input, std::vector<std::vector<std::string>>& output,
		std::vector<std::int64_t>& sourceLines, char delimiter = ',')
	{
		if (delimiter == '"' || delimiter == '\r' || delimiter == '\n')
			return SerializeFailure(IOError::INVALID_PARAMETERS, "CSV delimiter cannot be a quote or newline");

		enum class State { FieldStart, Unquoted, Quoted, Closed };
		State state = State::FieldStart;
		std::vector<std::vector<std::string>> records;
		std::vector<std::int64_t> lines;
		std::vector<std::string> fields;
		std::string field, line;
		std::int64_t lineNumber = 0, recordLine = 1, quoteLine = 1;
		const auto malformed = [&](const std::string& message, std::int64_t location) {
			return SerializeFailure(IOError::MALFORMED_INPUT, "CSV line " + std::to_string(location) +
				", field " + std::to_string(fields.size() + 1) + ": " + message)
				.At(location, IOOffsetKind::Line);
		};

		while (std::getline(input, line)) {
			++lineNumber;
			if (lineNumber == 1 && line.compare(0, 3, "\xEF\xBB\xBF") == 0)
				line.erase(0, 3);
			if (state != State::Quoted) {
				if (line.find_first_not_of(" \t\r") == std::string::npos)
					continue;
				recordLine = lineNumber;
			}

			const bool hasLF = !input.eof();
			const bool hasCR = !line.empty() && line.back() == '\r';
			const size_t length = line.size() - (hasCR ? 1 : 0);
			for (size_t i = 0; i < length; ++i) {
				const char c = line[i];
				if (c == delimiter && state != State::Quoted) {
					fields.push_back(std::move(field));
					field.clear();
					state = State::FieldStart;
				}
				else if (state == State::Quoted) {
					if (c == '"')
						state = State::Closed;
					else
						field += c;
				}
				else if (state == State::Closed) {
					if (c != '"')
						return malformed("Unexpected character after closing quote", lineNumber);
					field += '"';
					state = State::Quoted;
				}
				else if (c == '"') {
					if (state != State::FieldStart)
						return malformed("Quote in unquoted field", lineNumber);
					quoteLine = lineNumber;
					state = State::Quoted;
				}
				else {
					field += c;
					state = State::Unquoted;
				}
			}

			if (state == State::Quoted) {
				if (hasCR) field += '\r';
				if (hasLF) field += '\n';
				continue;
			}
			fields.push_back(std::move(field));
			field.clear();
			records.push_back(std::move(fields));
			fields.clear();
			lines.push_back(recordLine);
			state = State::FieldStart;
		}
		if (input.bad() || !input.eof())
			return SerializeFailure(IOError::READ_FAILED, "CSV stream read failed")
				.At(lineNumber + 1, IOOffsetKind::Line);
		if (state == State::Quoted)
			return malformed("Unclosed quoted field (record starts at line " +
				std::to_string(recordLine) + ")", quoteLine);

		output = std::move(records);
		sourceLines = std::move(lines);
		return IOResult::Success();
	}
}

#endif // MML_CSV_UTILS_H
