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

#include <string>
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

	/// @brief Split one CSV record into fields, honoring quoted fields and doubled quotes.
	/// @note Limitation: operates on a single line - quoted fields containing embedded
	///       newlines must be joined by the caller before splitting (record-oriented
	///       reading is not provided here).
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
}

#endif // MML_CSV_UTILS_H
