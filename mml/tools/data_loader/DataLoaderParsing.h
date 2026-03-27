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

#include "tools/data_loader/DataLoaderTypes.h"

#include <algorithm>
#include <cctype>
#include <regex>
#include <string>
#include <vector>

namespace MML {
	namespace Data {

		/////////////////////////////////////////////////////////////////////////////////////
		///                              TYPE INFERENCE                                    ///
		/////////////////////////////////////////////////////////////////////////////////////
		/// @brief Infer column type from string values
		/// @param values Sample of string values from the column
		/// @return Most appropriate ColumnType for the data
		inline ColumnType InferColumnType(const std::vector<std::string>& values) {
			if (values.empty()) return ColumnType::STRING;

			// Count type matches
			size_t intCount = 0, realCount = 0, boolCount = 0;
			size_t dateCount = 0, timeCount = 0, dateTimeCount = 0;
			size_t missingCount = 0;

			// Patterns
			std::regex intPattern(R"(^[-+]?\d+$)");
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
				if (std::regex_match(trimmed, intPattern)) {
					++intCount;
					continue;
				}
				if (std::regex_match(trimmed, realPattern)) {
					++realCount;
					continue;
				}

				// Check bool patterns (text-based only, not 1/0 which are handled as int)
				std::string lower = trimmed;
				std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
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
			if (intCount >= threshold * static_cast<double>(totalNonMissing)) return ColumnType::INT;
			if ((intCount + realCount) >= threshold * static_cast<double>(totalNonMissing)) return ColumnType::REAL;

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
		inline ParseResult ParseValue(const std::string& value, ColumnType targetType,
		                       Real& realOut, int& intOut, bool& boolOut, std::string& stringOut) {
			std::string trimmed = value;
			trimmed.erase(0, trimmed.find_first_not_of(" \t\r\n"));
			if (!trimmed.empty())
				trimmed.erase(trimmed.find_last_not_of(" \t\r\n") + 1);

			if (IsMissingValue(trimmed))
				return ParseResult::Missing;

			try {
				switch (targetType) {
					case ColumnType::REAL: {
						realOut = std::stod(trimmed);
						return ParseResult::Parsed;
					}
					case ColumnType::INT: {
						intOut = std::stoi(trimmed);
						return ParseResult::Parsed;
					}
					case ColumnType::BOOL: {
						std::string lower = trimmed;
						std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
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
			catch (...) {
				return ParseResult::Invalid;  // Parse error
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////
		///                              STRING UTILITIES                                  ///
		/////////////////////////////////////////////////////////////////////////////////////
		/// @brief Split a string by delimiter, respecting quoted fields
		inline std::vector<std::string> SplitLine(const std::string& line, char delimiter) {
			std::vector<std::string> result;
			std::string field;
			bool inQuotes = false;
			bool prevWasQuote = false;

			for (size_t i = 0; i < line.size(); ++i) {
				char c = line[i];

				if (c == '"') {
					if (inQuotes && i + 1 < line.size() && line[i + 1] == '"') {
						// Escaped quote
						field += '"';
						++i;
					}
					else {
						inQuotes = !inQuotes;
					}
					prevWasQuote = true;
				}
				else if (c == delimiter && !inQuotes) {
					result.push_back(field);
					field.clear();
					prevWasQuote = false;
				}
				else {
					field += c;
					prevWasQuote = false;
				}
			}
			result.push_back(field);

			return result;
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
