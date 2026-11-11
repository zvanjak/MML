///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DataLoaderParsing.h                                                 ///
///  Description: Type inference and value parsing utilities for data loading        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DATA_LOADER_PARSING_H
#define MML_DATA_LOADER_PARSING_H

#include <mml/tools/data_loader/DataLoaderTypes.h>
#include <mml/tools/CsvUtils.h>
#include <mml/tools/TextIOFormat.h>

#include <algorithm>
#include <charconv>
#include <iomanip>
#include <regex>
#include <sstream>
#include <string>
#include <system_error>
#include <vector>

namespace MML {
	namespace Data {
		namespace Detail {
			inline bool TryParseInteger(const std::string& text, int& output) {
				const char* first = text.data();
				const char* last = first + text.size();
				if (first != last && *first == '+') {
					++first;
					if (first != last && *first == '-') return false;
				}
				int value = 0;
				const auto result = std::from_chars(first, last, value);
				if (result.ec != std::errc{} || result.ptr != last)
					return false;
				output = value;
				return true;
			}

			inline bool IsIntegerToken(const std::string& text) {
				size_t first = !text.empty() && (text[0] == '+' || text[0] == '-') ? 1 : 0;
				return first < text.size() && std::all_of(text.begin() + first, text.end(),
					[](char c) { return c >= '0' && c <= '9'; });
			}

			inline std::string NormalizeIntegerToken(const std::string& text) {
				const size_t first = text[0] == '+' || text[0] == '-' ? 1 : 0;
				const size_t significant = text.find_first_not_of('0', first);
				if (significant == std::string::npos) return "0";
				return (text[0] == '-' ? "-" : "") + text.substr(significant);
			}

			inline bool IntegerMatchesReal(const std::string& text, Real value) {
				std::ostringstream integer;
				TextIO::ApplyClassicFormat(integer);
				// Fixed, zero-decimal formatting exposes the stored integer, not a rounded significant-digit spelling.
				integer << std::fixed << std::setprecision(0) << value;
				return integer.good() && NormalizeIntegerToken(text) == NormalizeIntegerToken(integer.str());
			}

			inline IOResult InvalidValueResult(const std::string& value, const DataColumn& column,
				size_t rowIndex, size_t columnIndex) {
				return SerializeFailure(IOError::MALFORMED_INPUT,
					"Invalid " + ColumnTypeToString(column.type) + " value '" + value +
					"' at data row " + std::to_string(rowIndex + 1) + ", column " +
					std::to_string(columnIndex + 1) + " ('" + column.name + "')");
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////
		///                              TYPE INFERENCE                                    ///
		/////////////////////////////////////////////////////////////////////////////////////
		/// @brief Infer column type from string values
		/// @param values Sample of string values from the column
		/// @return Most appropriate ColumnType for the data
		/// @details Numeric majorities use REAL if any numeric token is fractional,
		/// exponent-spelled, or outside int range. The existing 90% threshold is retained.
		inline ColumnType InferColumnType(const std::vector<std::string>& values) {
			if (values.empty()) return ColumnType::STRING;

			// Count type matches
			size_t intCount = 0, realCount = 0, boolCount = 0;
			size_t dateCount = 0, timeCount = 0, dateTimeCount = 0;
			size_t missingCount = 0;
			bool requiresReal = false;

			// Patterns
			std::regex realPattern(R"(^[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$)");
			std::regex datePattern(R"(^\d{4}[-/]\d{1,2}[-/]\d{1,2}$)");
			std::regex timePattern(R"(^\d{1,2}:\d{2}(:\d{2})?(\.\d+)?$)");
			std::regex dateTimePattern(R"(^\d{4}[-/]\d{1,2}[-/]\d{1,2}[T ]\d{1,2}:\d{2}(:\d{2})?$)");

			for (const auto& val : values) {
				std::string trimmed = val;
				// Trim whitespace
				trimmed.erase(0, trimmed.find_first_not_of(" \t\r\n"));
				trimmed.erase(trimmed.find_last_not_of(" \t\r\n") + 1);

				if (IsMissingValue(trimmed)) {
					++missingCount;
					continue;
				}

				// Check numeric patterns FIRST (before bool, since "1"/"0" could be numeric)
				if (Detail::IsIntegerToken(trimmed)) {
					++intCount;
					int integer = 0;
					if (!Detail::TryParseInteger(trimmed, integer))
						requiresReal = true;
					continue;
				}
				if (std::regex_match(trimmed, realPattern)) {
					++realCount;
					continue;
				}

				// Check bool patterns (text-based only, not 1/0 which are handled as int)
				std::string lower = trimmed;
				std::transform(lower.begin(), lower.end(), lower.begin(), TextIO::ToLowerAscii);
				if (lower == "true" || lower == "false" || lower == "yes" || lower == "no" ||
				    lower == "t" || lower == "f" || lower == "y" || lower == "n") {
					++boolCount;
					continue;
				}

				// Check date/time patterns
				if (std::regex_match(trimmed, dateTimePattern)) {
					++dateTimeCount;
					continue;
				}
				if (std::regex_match(trimmed, datePattern)) {
					++dateCount;
					continue;
				}
				if (std::regex_match(trimmed, timePattern)) {
					++timeCount;
					continue;
				}
			}

			size_t totalNonMissing = values.size() - missingCount;
			if (totalNonMissing == 0) return ColumnType::STRING;

			// Determine type based on majority
			double threshold = 0.9;  // 90% match required

			if (dateTimeCount >= threshold * static_cast<double>(totalNonMissing)) return ColumnType::DATETIME;
			if (dateCount >= threshold * static_cast<double>(totalNonMissing)) return ColumnType::DATE;
			if (timeCount >= threshold * static_cast<double>(totalNonMissing)) return ColumnType::TIME;
			if (boolCount >= threshold * static_cast<double>(totalNonMissing)) return ColumnType::BOOL;
			if ((intCount + realCount) >= threshold * static_cast<double>(totalNonMissing))
				return realCount != 0 || requiresReal ? ColumnType::REAL : ColumnType::INT;

			return ColumnType::STRING;
		}

		/////////////////////////////////////////////////////////////////////////////////////
		///                              VALUE PARSING                                     ///
		/////////////////////////////////////////////////////////////////////////////////////

		/// @brief Result of parsing a string value to a typed value
		enum class ParseResult {
			Parsed,   ///< Successfully parsed to target type
			Missing,  ///< Value was empty or a recognized missing-value sentinel
			Invalid   ///< Value was present but could not be parsed as the target type
		};

		/// @brief Parse a string value to the target type
		/// @param value String value to parse
		/// @param targetType Target column type
		/// @param[out] realOut Real output (if REAL type)
		/// @param[out] intOut Integer output (if INT type)
		/// @param[out] boolOut Boolean output (if BOOL type)
		/// @param[out] stringOut String output (if STRING/DATE/TIME type)
		/// @return ParseResult indicating success, missing value, or parse failure
		/// @note Integer tokens require whole-token conversion and exact representation
		/// in the selected numeric storage. Decimal/exponent tokens use normal Real rounding.
		inline ParseResult ParseValue(const std::string& value, ColumnType targetType,
		                       Real& realOut, int& intOut, bool& boolOut, std::string& stringOut) {
			std::string trimmed = value;
			trimmed.erase(0, trimmed.find_first_not_of(" \t\r\n"));
			if (!trimmed.empty())
				trimmed.erase(trimmed.find_last_not_of(" \t\r\n") + 1);

			if (IsMissingValue(trimmed))
				return ParseResult::Missing;

			switch (targetType) {
				case ColumnType::REAL: {
					Real value = 0;
					if (!TextIO::ParseReal(trimmed, value) || !std::isfinite(value))
						return ParseResult::Invalid;
					if (Detail::IsIntegerToken(trimmed) && !Detail::IntegerMatchesReal(trimmed, value))
						return ParseResult::Invalid;
					realOut = value;
					return ParseResult::Parsed;
				}
				case ColumnType::INT:
					return Detail::TryParseInteger(trimmed, intOut) ? ParseResult::Parsed : ParseResult::Invalid;
				case ColumnType::BOOL: {
					std::string lower = trimmed;
					std::transform(lower.begin(), lower.end(), lower.begin(), TextIO::ToLowerAscii);
					if (lower == "true" || lower == "yes" || lower == "t" || lower == "y" || lower == "1") {
						boolOut = true;
						return ParseResult::Parsed;
					}
					if (lower == "false" || lower == "no" || lower == "f" || lower == "n" || lower == "0") {
						boolOut = false;
						return ParseResult::Parsed;
					}
					return ParseResult::Invalid;
				}
				case ColumnType::STRING:
				case ColumnType::DATE:
				case ColumnType::TIME:
				case ColumnType::DATETIME:
					stringOut = trimmed;
					return ParseResult::Parsed;
				default:
					return ParseResult::Invalid;
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////
		///                              STRING UTILITIES                                  ///
		/////////////////////////////////////////////////////////////////////////////////////
		/// @brief Split a string by delimiter, respecting quoted fields
		/// (delegates to the shared RFC 4180 module in tools/CsvUtils.h)
		inline std::vector<std::string> SplitLine(const std::string& line, char delimiter) {
			return CsvUtils::SplitLine(line, delimiter);
		}

		/// @brief Trim whitespace from string
		inline std::string Trim(const std::string& s) {
			size_t start = s.find_first_not_of(" \t\r\n");
			if (start == std::string::npos) return "";
			size_t end = s.find_last_not_of(" \t\r\n");
			return s.substr(start, end - start + 1);
		}

		/// @brief Remove UTF-8 BOM if present
		inline std::string RemoveBOM(const std::string& s) {
			if (s.size() >= 3 &&
			    static_cast<unsigned char>(s[0]) == 0xEF &&
			    static_cast<unsigned char>(s[1]) == 0xBB &&
			    static_cast<unsigned char>(s[2]) == 0xBF) {
				return s.substr(3);
			}
			return s;
		}

	}  // namespace Data
}  // namespace MML

#endif  // MML_DATA_LOADER_PARSING_H
