///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DataLoader.h                                                        ///
///  Description: Multi-format dataset loading infrastructure for CSV, TSV, JSON     ///
///               with automatic type inference and pandas-style data access         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DATA_LOADER_H
#define MML_DATA_LOADER_H

#include "MMLBase.h"
#include "MMLExceptions.h"
#include "base/Vector.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <numeric>
#include <regex>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

namespace MML {
	namespace Data {

		/////////////////////////////////////////////////////////////////////////////////////
		///                              DATA FORMAT ENUM                                 ///
		/////////////////////////////////////////////////////////////////////////////////////

		/// @brief Supported data file formats
		enum class DataFormat {
			CSV,      ///< Comma-separated values
			TSV,      ///< Tab-separated values
			JSON,     ///< JSON format (array of objects)
			Auto      ///< Auto-detect from file extension
		};

		/////////////////////////////////////////////////////////////////////////////////////
		///                              COLUMN TYPE ENUM                                  ///
		/////////////////////////////////////////////////////////////////////////////////////

		/// @brief Data types for columns
		enum class ColumnType {
			REAL,     ///< Floating-point numbers
			INT,      ///< Integer numbers
			BOOL,     ///< Boolean values (true/false, yes/no, 1/0)
			STRING,   ///< Text strings
			DATE,     ///< Date (YYYY-MM-DD)
			TIME,     ///< Time (HH:MM:SS)
			DATETIME, ///< Combined date and time
			MISSING   ///< Missing/null value marker
		};

		/// @brief Convert ColumnType to string for display
		inline std::string ColumnTypeToString(ColumnType type) {
			switch (type) {
				case ColumnType::REAL:     return "Real";
				case ColumnType::INT:      return "Int";
				case ColumnType::BOOL:     return "Bool";
				case ColumnType::STRING:   return "String";
				case ColumnType::DATE:     return "Date";
				case ColumnType::TIME:     return "Time";
				case ColumnType::DATETIME: return "DateTime";
				case ColumnType::MISSING:  return "Missing";
				default:                   return "Unknown";
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////
		///                              MISSING VALUE                                     ///
		/////////////////////////////////////////////////////////////////////////////////////

		/// @brief Marker for missing values
		struct MissingValue {
			bool operator==(const MissingValue&) const { return true; }
		};

		/// @brief Check if a string represents a missing value
		inline bool IsMissingValue(const std::string& s) {
			if (s.empty()) return true;
			std::string lower = s;
			std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
			return lower == "na" || lower == "n/a" || lower == "nan" ||
			       lower == "null" || lower == "none" || lower == "." ||
			       lower == "?" || lower == "-" || lower == "";
		}

		/////////////////////////////////////////////////////////////////////////////////////
		///                              DATA COLUMN STRUCT                                ///
		/////////////////////////////////////////////////////////////////////////////////////

		/// @brief A single column of data with typed storage
		struct DataColumn {
			std::string name;           ///< Column name/header
			ColumnType type;            ///< Inferred or specified type
			
			// Typed storage vectors - only one is populated based on type
			Vector<Real> realData;
			Vector<int> intData;
			std::vector<bool> boolData;           // Use std::vector for bool (MML Vector issues)
			std::vector<std::string> stringData;
			std::vector<std::string> dateData;    // Stored as YYYY-MM-DD strings
			std::vector<std::string> timeData;    // Stored as HH:MM:SS strings
			std::vector<bool> missingMask;        // true where value is missing

			/// @brief Get number of rows
			size_t Size() const {
				switch (type) {
					case ColumnType::REAL:     return realData.size();
					case ColumnType::INT:      return intData.size();
					case ColumnType::BOOL:     return boolData.size();
					case ColumnType::STRING:   return stringData.size();
					case ColumnType::DATE:     return dateData.size();
					case ColumnType::TIME:     return timeData.size();
					case ColumnType::DATETIME: return dateData.size();
					default:                   return missingMask.size();
				}
			}

			/// @brief Check if value at index is missing
			bool IsMissing(size_t index) const {
				if (index >= missingMask.size()) return false;
				return missingMask[index];
			}

			/// @brief Count missing values
			size_t CountMissing() const {
				size_t count = 0;
				for (bool m : missingMask) {
					if (m) ++count;
				}
				return count;
			}

			/// @brief Get value as string (for any type)
			std::string GetAsString(size_t index) const {
				if (index >= Size()) return "";
				if (IsMissing(index)) return "NA";

				switch (type) {
					case ColumnType::REAL: {
						std::ostringstream oss;
						oss << std::setprecision(6) << realData[index];
						return oss.str();
					}
					case ColumnType::INT:
						return std::to_string(intData[index]);
					case ColumnType::BOOL:
						return boolData[index] ? "true" : "false";
					case ColumnType::STRING:
						return stringData[index];
					case ColumnType::DATE:
						return dateData[index];
					case ColumnType::TIME:
						return timeData[index];
					case ColumnType::DATETIME:
						return dateData[index] + " " + timeData[index];
					default:
						return "NA";
				}
			}

			/// @brief Get real value (converts if possible)
			Real GetReal(size_t index) const {
				if (index >= Size() || IsMissing(index))
					return std::numeric_limits<Real>::quiet_NaN();

				switch (type) {
					case ColumnType::REAL:
						return realData[index];
					case ColumnType::INT:
						return static_cast<Real>(intData[index]);
					case ColumnType::BOOL:
						return boolData[index] ? 1.0 : 0.0;
					default:
						return std::numeric_limits<Real>::quiet_NaN();
				}
			}
		};

		/////////////////////////////////////////////////////////////////////////////////////
		///                              DATASET STRUCT                                    ///
		/////////////////////////////////////////////////////////////////////////////////////

		/// @brief A complete dataset with multiple columns
		struct Dataset {
			std::string name;                    ///< Dataset name
			std::vector<DataColumn> columns;     ///< Column data
			size_t rowCount;                     ///< Number of rows

			/// @brief Default constructor
			Dataset() : rowCount(0) {}

			/// @brief Get number of columns
			size_t NumColumns() const { return columns.size(); }

			/// @brief Get number of rows
			size_t NumRows() const { return rowCount; }

			/// @brief Access column by index
			DataColumn& operator[](size_t index) {
				if (index >= columns.size())
					throw DataError("Dataset: Column index " + std::to_string(index) + " out of range");
				return columns[index];
			}

			const DataColumn& operator[](size_t index) const {
				if (index >= columns.size())
					throw DataError("Dataset: Column index " + std::to_string(index) + " out of range");
				return columns[index];
			}

			/// @brief Access column by name
			DataColumn& operator[](const std::string& colName) {
				for (auto& col : columns) {
					if (col.name == colName) return col;
				}
				throw DataError("Dataset: Column '" + colName + "' not found");
			}

			const DataColumn& operator[](const std::string& colName) const {
				for (const auto& col : columns) {
					if (col.name == colName) return col;
				}
				throw DataError("Dataset: Column '" + colName + "' not found");
			}

			/// @brief Check if column exists
			bool HasColumn(const std::string& colName) const {
				for (const auto& col : columns) {
					if (col.name == colName) return true;
				}
				return false;
			}

			/// @brief Get column names
			std::vector<std::string> GetColumnNames() const {
				std::vector<std::string> names;
				names.reserve(columns.size());
				for (const auto& col : columns) {
					names.push_back(col.name);
				}
				return names;
			}

			/// @brief Get real-valued column as Vector<Real>
			/// @throws DataError if column not found or not convertible to Real
			Vector<Real> GetRealColumn(const std::string& colName) const {
				const DataColumn& col = (*this)[colName];
				if (col.type != ColumnType::REAL && col.type != ColumnType::INT)
					throw DataError("Dataset: Column '" + colName + "' is not numeric");

				if (col.type == ColumnType::REAL)
					return col.realData;

				// Convert int to real
				Vector<Real> result(col.intData.size());
				for (size_t i = 0; i < col.intData.size(); ++i)
					result[i] = static_cast<Real>(col.intData[i]);
				return result;
			}

			/// @brief Get integer column
			Vector<int> GetIntColumn(const std::string& colName) const {
				const DataColumn& col = (*this)[colName];
				if (col.type != ColumnType::INT)
					throw DataError("Dataset: Column '" + colName + "' is not integer type");
				return col.intData;
			}

			/// @brief Get string column
			std::vector<std::string> GetStringColumn(const std::string& colName) const {
				const DataColumn& col = (*this)[colName];
				// Can convert any column to string
				std::vector<std::string> result;
				result.reserve(rowCount);
				for (size_t i = 0; i < rowCount; ++i)
					result.push_back(col.GetAsString(i));
				return result;
			}

			/// @brief Get bool column
			std::vector<bool> GetBoolColumn(const std::string& colName) const {
				const DataColumn& col = (*this)[colName];
				if (col.type != ColumnType::BOOL)
					throw DataError("Dataset: Column '" + colName + "' is not boolean type");
				return col.boolData;
			}

			/// @brief Print dataset summary (pandas describe() style)
			std::string PrintSummary() const {
				std::ostringstream oss;
				oss << "Dataset: " << (name.empty() ? "(unnamed)" : name) << "\n";
				oss << "Shape: " << rowCount << " rows x " << columns.size() << " columns\n\n";

				// Column info
				oss << std::left << std::setw(20) << "Column"
				    << std::setw(12) << "Type"
				    << std::setw(12) << "Non-Null"
				    << std::setw(12) << "Missing"
				    << "\n";
				oss << std::string(56, '-') << "\n";

				for (const auto& col : columns) {
					size_t missing = col.CountMissing();
					oss << std::left << std::setw(20) << col.name
					    << std::setw(12) << ColumnTypeToString(col.type)
					    << std::setw(12) << (rowCount - missing)
					    << std::setw(12) << missing
					    << "\n";
				}

				// Numeric column statistics
				oss << "\nNumeric Column Statistics:\n";
				oss << std::string(80, '-') << "\n";
				oss << std::left << std::setw(15) << "Column"
				    << std::right << std::setw(12) << "Min"
				    << std::setw(12) << "Max"
				    << std::setw(12) << "Mean"
				    << std::setw(12) << "Std"
				    << "\n";

				for (const auto& col : columns) {
					if (col.type == ColumnType::REAL || col.type == ColumnType::INT) {
						Vector<Real> data = (col.type == ColumnType::REAL) ? col.realData :
							[&col]() {
								Vector<Real> r(col.intData.size());
								for (size_t i = 0; i < col.intData.size(); ++i)
									r[i] = static_cast<Real>(col.intData[i]);
								return r;
							}();

						// Filter out NaN/missing
						std::vector<Real> valid;
						for (size_t i = 0; i < data.size(); ++i) {
							if (!col.IsMissing(i) && !std::isnan(data[i]))
								valid.push_back(data[i]);
						}

						if (!valid.empty()) {
							Real minVal = *std::min_element(valid.begin(), valid.end());
							Real maxVal = *std::max_element(valid.begin(), valid.end());
							Real sum = std::accumulate(valid.begin(), valid.end(), 0.0);
							Real mean = sum / valid.size();
							Real sqSum = 0.0;
							for (Real v : valid) sqSum += (v - mean) * (v - mean);
							Real std = std::sqrt(sqSum / valid.size());

							oss << std::left << std::setw(15) << col.name
							    << std::right << std::setw(12) << std::setprecision(4) << minVal
							    << std::setw(12) << maxVal
							    << std::setw(12) << mean
							    << std::setw(12) << std
							    << "\n";
						}
					}
				}

				return oss.str();
			}

			/// @brief Get first n rows as string representation
			std::string Head(size_t n = 5) const {
				return RowsToString(0, std::min(n, rowCount));
			}

			/// @brief Get last n rows as string representation
			std::string Tail(size_t n = 5) const {
				size_t start = rowCount > n ? rowCount - n : 0;
				return RowsToString(start, rowCount);
			}

			/// @brief Get random sample of n rows
			std::string Sample(size_t n = 5, unsigned int seed = 42) const {
				if (n >= rowCount)
					return RowsToString(0, rowCount);

				// Simple LCG random selection
				std::vector<size_t> indices;
				indices.reserve(rowCount);
				for (size_t i = 0; i < rowCount; ++i)
					indices.push_back(i);

				// Fisher-Yates shuffle first n elements
				unsigned int state = seed;
				for (size_t i = 0; i < n; ++i) {
					state = state * 1103515245u + 12345u;
					size_t j = i + (state % (rowCount - i));
					std::swap(indices[i], indices[j]);
				}

				std::ostringstream oss;
				// Header
				for (size_t c = 0; c < columns.size(); ++c) {
					if (c > 0) oss << "\t";
					oss << columns[c].name;
				}
				oss << "\n";

				// Selected rows
				for (size_t i = 0; i < n; ++i) {
					size_t row = indices[i];
					for (size_t c = 0; c < columns.size(); ++c) {
						if (c > 0) oss << "\t";
						oss << columns[c].GetAsString(row);
					}
					oss << "\n";
				}

				return oss.str();
			}

			/// @brief Filter rows by predicate on a column
			/// @param colName Column to filter on
			/// @param predicate Function returning true for rows to keep
			Dataset FilterRows(const std::string& colName,
			                   std::function<bool(const std::string&)> predicate) const {
				const DataColumn& filterCol = (*this)[colName];

				// Find indices to keep
				std::vector<size_t> keepIndices;
				for (size_t i = 0; i < rowCount; ++i) {
					if (predicate(filterCol.GetAsString(i)))
						keepIndices.push_back(i);
				}

				// Build new dataset
				Dataset result;
				result.name = name + "_filtered";
				result.rowCount = keepIndices.size();
				result.columns.reserve(columns.size());

				for (const auto& srcCol : columns) {
					DataColumn newCol;
					newCol.name = srcCol.name;
					newCol.type = srcCol.type;
					newCol.missingMask.reserve(result.rowCount);

					for (size_t idx : keepIndices) {
						newCol.missingMask.push_back(srcCol.IsMissing(idx));
					}

					switch (srcCol.type) {
						case ColumnType::REAL:
							newCol.realData = Vector<Real>(result.rowCount);
							for (size_t i = 0; i < keepIndices.size(); ++i)
								newCol.realData[i] = srcCol.realData[keepIndices[i]];
							break;
						case ColumnType::INT:
							newCol.intData = Vector<int>(result.rowCount);
							for (size_t i = 0; i < keepIndices.size(); ++i)
								newCol.intData[i] = srcCol.intData[keepIndices[i]];
							break;
						case ColumnType::BOOL:
							newCol.boolData.resize(result.rowCount);
							for (size_t i = 0; i < keepIndices.size(); ++i)
								newCol.boolData[i] = srcCol.boolData[keepIndices[i]];
							break;
						case ColumnType::STRING:
							newCol.stringData.reserve(result.rowCount);
							for (size_t idx : keepIndices)
								newCol.stringData.push_back(srcCol.stringData[idx]);
							break;
						case ColumnType::DATE:
						case ColumnType::DATETIME:
							newCol.dateData.reserve(result.rowCount);
							for (size_t idx : keepIndices)
								newCol.dateData.push_back(srcCol.dateData[idx]);
							if (srcCol.type == ColumnType::DATETIME) {
								newCol.timeData.reserve(result.rowCount);
								for (size_t idx : keepIndices)
									newCol.timeData.push_back(srcCol.timeData[idx]);
							}
							break;
						case ColumnType::TIME:
							newCol.timeData.reserve(result.rowCount);
							for (size_t idx : keepIndices)
								newCol.timeData.push_back(srcCol.timeData[idx]);
							break;
						default:
							break;
					}

					result.columns.push_back(std::move(newCol));
				}

				return result;
			}

			/// @brief Select specific columns
			Dataset SelectColumns(const std::vector<std::string>& colNames) const {
				Dataset result;
				result.name = name + "_selected";
				result.rowCount = rowCount;
				result.columns.reserve(colNames.size());

				for (const auto& colName : colNames) {
					result.columns.push_back((*this)[colName]);
				}

				return result;
			}

		private:
			/// @brief Convert range of rows to string
			std::string RowsToString(size_t startRow, size_t endRow) const {
				std::ostringstream oss;

				// Header
				for (size_t c = 0; c < columns.size(); ++c) {
					if (c > 0) oss << "\t";
					oss << columns[c].name;
				}
				oss << "\n";

				// Data rows
				for (size_t r = startRow; r < endRow; ++r) {
					for (size_t c = 0; c < columns.size(); ++c) {
						if (c > 0) oss << "\t";
						oss << columns[c].GetAsString(r);
					}
					oss << "\n";
				}

				return oss.str();
			}
		};

		/////////////////////////////////////////////////////////////////////////////////////
		///                              TYPE INFERENCE                                    ///
		/////////////////////////////////////////////////////////////////////////////////////

		/// @brief Infer column type from string values
		/// @param values Sample of string values from the column
		/// @return Most appropriate ColumnType for the data
		inline ColumnType InferColumnType(const std::vector<std::string>& values) {
			if (values.empty()) return ColumnType::STRING;

			// Count type matches
			int intCount = 0, realCount = 0, boolCount = 0;
			int dateCount = 0, timeCount = 0, dateTimeCount = 0;
			int missingCount = 0;

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

			int totalNonMissing = static_cast<int>(values.size()) - missingCount;
			if (totalNonMissing == 0) return ColumnType::STRING;

			// Determine type based on majority
			double threshold = 0.9;  // 90% match required

			if (dateTimeCount >= threshold * totalNonMissing) return ColumnType::DATETIME;
			if (dateCount >= threshold * totalNonMissing) return ColumnType::DATE;
			if (timeCount >= threshold * totalNonMissing) return ColumnType::TIME;
			if (boolCount >= threshold * totalNonMissing) return ColumnType::BOOL;
			if (intCount >= threshold * totalNonMissing) return ColumnType::INT;
			if ((intCount + realCount) >= threshold * totalNonMissing) return ColumnType::REAL;

			return ColumnType::STRING;
		}

		/////////////////////////////////////////////////////////////////////////////////////
		///                              VALUE PARSING                                     ///
		/////////////////////////////////////////////////////////////////////////////////////

		/// @brief Parse a string value to the target type
		/// @param value String value to parse
		/// @param targetType Target column type
		/// @param[out] realOut Real output (if REAL type)
		/// @param[out] intOut Integer output (if INT type)
		/// @param[out] boolOut Boolean output (if BOOL type)
		/// @param[out] stringOut String output (if STRING/DATE/TIME type)
		/// @return true if parsing succeeded, false if value should be marked as missing
		inline bool ParseValue(const std::string& value, ColumnType targetType,
		                       Real& realOut, int& intOut, bool& boolOut, std::string& stringOut) {
			std::string trimmed = value;
			trimmed.erase(0, trimmed.find_first_not_of(" \t\r\n"));
			if (!trimmed.empty())
				trimmed.erase(trimmed.find_last_not_of(" \t\r\n") + 1);

			if (IsMissingValue(trimmed))
				return false;

			try {
				switch (targetType) {
					case ColumnType::REAL: {
						realOut = std::stod(trimmed);
						return true;
					}
					case ColumnType::INT: {
						intOut = std::stoi(trimmed);
						return true;
					}
					case ColumnType::BOOL: {
						std::string lower = trimmed;
						std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
						if (lower == "true" || lower == "yes" || lower == "t" || lower == "y" || lower == "1") {
							boolOut = true;
							return true;
						}
						if (lower == "false" || lower == "no" || lower == "f" || lower == "n" || lower == "0") {
							boolOut = false;
							return true;
						}
						return false;
					}
					case ColumnType::STRING:
					case ColumnType::DATE:
					case ColumnType::TIME:
					case ColumnType::DATETIME:
						stringOut = trimmed;
						return true;
					default:
						return false;
				}
			}
			catch (...) {
				return false;  // Parse error -> missing
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

		/////////////////////////////////////////////////////////////////////////////////////
		///                              CSV/TSV LOADING                                   ///
		/////////////////////////////////////////////////////////////////////////////////////

		/// @brief Load dataset from CSV file
		/// @param filename Path to CSV file
		/// @param hasHeader If true, first row is column names
		/// @param delimiter Field delimiter (default ',')
		/// @param inferTypes If true, automatically infer column types
		/// @return Loaded dataset
		/// @throws DataError on file/parse errors
		inline Dataset LoadCSV(const std::string& filename, bool hasHeader = true,
		                       char delimiter = ',', bool inferTypes = true) {
			std::ifstream file(filename);
			if (!file.is_open())
				throw DataError("LoadCSV: Cannot open file '" + filename + "'");

			Dataset dataset;
			dataset.name = filename;

			std::vector<std::vector<std::string>> allData;
			std::string line;
			bool firstLine = true;

			while (std::getline(file, line)) {
				// Handle BOM on first line
				if (firstLine) {
					line = RemoveBOM(line);
					firstLine = false;
				}

				// Skip empty lines
				if (Trim(line).empty()) continue;

				// Handle Windows line endings
				if (!line.empty() && line.back() == '\r')
					line.pop_back();

				std::vector<std::string> fields = SplitLine(line, delimiter);
				allData.push_back(fields);
			}

			if (allData.empty())
				throw DataError("LoadCSV: File is empty");

			// Determine column count
			size_t numCols = allData[0].size();
			size_t dataStartRow = hasHeader ? 1 : 0;

			// Setup columns
			dataset.columns.resize(numCols);

			for (size_t c = 0; c < numCols; ++c) {
				if (hasHeader) {
					dataset.columns[c].name = Trim(allData[0][c]);
					if (dataset.columns[c].name.empty())
						dataset.columns[c].name = "Column" + std::to_string(c);
				}
				else {
					dataset.columns[c].name = "Column" + std::to_string(c);
				}
			}

			// Collect string values for type inference
			if (inferTypes) {
				for (size_t c = 0; c < numCols; ++c) {
					std::vector<std::string> colValues;
					for (size_t r = dataStartRow; r < allData.size(); ++r) {
						if (c < allData[r].size())
							colValues.push_back(allData[r][c]);
						else
							colValues.push_back("");
					}
					dataset.columns[c].type = InferColumnType(colValues);
				}
			}
			else {
				// Default to STRING
				for (size_t c = 0; c < numCols; ++c)
					dataset.columns[c].type = ColumnType::STRING;
			}

			// Parse data
			dataset.rowCount = allData.size() - dataStartRow;

			for (size_t c = 0; c < numCols; ++c) {
				DataColumn& col = dataset.columns[c];
				col.missingMask.resize(dataset.rowCount, false);

				// Pre-allocate
				switch (col.type) {
					case ColumnType::REAL:
						col.realData = Vector<Real>(dataset.rowCount);
						break;
					case ColumnType::INT:
						col.intData = Vector<int>(dataset.rowCount);
						break;
					case ColumnType::BOOL:
						col.boolData.resize(dataset.rowCount);
						break;
					case ColumnType::STRING:
					case ColumnType::DATE:
					case ColumnType::TIME:
						col.stringData.resize(dataset.rowCount);
						break;
					case ColumnType::DATETIME:
						col.dateData.resize(dataset.rowCount);
						col.timeData.resize(dataset.rowCount);
						break;
					default:
						break;
				}
			}

			for (size_t r = dataStartRow; r < allData.size(); ++r) {
				size_t rowIdx = r - dataStartRow;
				const auto& row = allData[r];

				for (size_t c = 0; c < numCols; ++c) {
					DataColumn& col = dataset.columns[c];
					std::string value = (c < row.size()) ? row[c] : "";

					Real realVal = 0.0;
					int intVal = 0;
					bool boolVal = false;
					std::string strVal;

					bool parsed = ParseValue(value, col.type, realVal, intVal, boolVal, strVal);
					col.missingMask[rowIdx] = !parsed;

					switch (col.type) {
						case ColumnType::REAL:
							col.realData[rowIdx] = parsed ? realVal : std::numeric_limits<Real>::quiet_NaN();
							break;
						case ColumnType::INT:
							col.intData[rowIdx] = parsed ? intVal : 0;
							break;
						case ColumnType::BOOL:
							col.boolData[rowIdx] = parsed ? boolVal : false;
							break;
						case ColumnType::STRING:
							col.stringData[rowIdx] = parsed ? strVal : "";
							break;
						case ColumnType::DATE:
							col.dateData.resize(dataset.rowCount);
							col.dateData[rowIdx] = parsed ? strVal : "";
							break;
						case ColumnType::TIME:
							col.timeData.resize(dataset.rowCount);
							col.timeData[rowIdx] = parsed ? strVal : "";
							break;
						case ColumnType::DATETIME:
							if (parsed && !strVal.empty()) {
								// Split datetime into date and time
								size_t sep = strVal.find_first_of("T ");
								if (sep != std::string::npos) {
									col.dateData[rowIdx] = strVal.substr(0, sep);
									col.timeData[rowIdx] = strVal.substr(sep + 1);
								}
								else {
									col.dateData[rowIdx] = strVal;
									col.timeData[rowIdx] = "";
								}
							}
							break;
						default:
							break;
					}
				}
			}

			return dataset;
		}

		/// @brief Load dataset from TSV file
		inline Dataset LoadTSV(const std::string& filename, bool hasHeader = true, bool inferTypes = true) {
			return LoadCSV(filename, hasHeader, '\t', inferTypes);
		}

		/////////////////////////////////////////////////////////////////////////////////////
		///                              JSON LOADING                                      ///
		/////////////////////////////////////////////////////////////////////////////////////

		/// @brief Simple JSON value types for parsing
		enum class JsonType { Null, Bool, Number, String, Array, Object };

		/// @brief Minimal JSON parser state
		struct JsonParser {
			const std::string& json;
			size_t pos;

			JsonParser(const std::string& s) : json(s), pos(0) {}

			void SkipWhitespace() {
				while (pos < json.size() && std::isspace(json[pos])) ++pos;
			}

			char Peek() {
				SkipWhitespace();
				return pos < json.size() ? json[pos] : '\0';
			}

			char Get() {
				SkipWhitespace();
				return pos < json.size() ? json[pos++] : '\0';
			}

			std::string ParseString() {
				if (Get() != '"')
					throw DataError("JSON: Expected '\"'");

				std::string result;
				while (pos < json.size()) {
					char c = json[pos++];
					if (c == '"') return result;
					if (c == '\\' && pos < json.size()) {
						char esc = json[pos++];
						switch (esc) {
							case 'n': result += '\n'; break;
							case 't': result += '\t'; break;
							case 'r': result += '\r'; break;
							case '"': result += '"'; break;
							case '\\': result += '\\'; break;
							case '/': result += '/'; break;
							default: result += esc; break;
						}
					}
					else {
						result += c;
					}
				}
				throw DataError("JSON: Unterminated string");
			}

			std::string ParseNumber() {
				SkipWhitespace();
				size_t start = pos;
				if (json[pos] == '-') ++pos;
				while (pos < json.size() && (std::isdigit(json[pos]) || json[pos] == '.' ||
				       json[pos] == 'e' || json[pos] == 'E' || json[pos] == '+' || json[pos] == '-')) {
					++pos;
				}
				return json.substr(start, pos - start);
			}

			bool ParseBool() {
				SkipWhitespace();
				if (json.substr(pos, 4) == "true") {
					pos += 4;
					return true;
				}
				if (json.substr(pos, 5) == "false") {
					pos += 5;
					return false;
				}
				throw DataError("JSON: Expected boolean");
			}

			void ParseNull() {
				SkipWhitespace();
				if (json.substr(pos, 4) == "null") {
					pos += 4;
					return;
				}
				throw DataError("JSON: Expected null");
			}

			JsonType PeekType() {
				char c = Peek();
				if (c == '"') return JsonType::String;
				if (c == '[') return JsonType::Array;
				if (c == '{') return JsonType::Object;
				if (c == 't' || c == 'f') return JsonType::Bool;
				if (c == 'n') return JsonType::Null;
				if (c == '-' || std::isdigit(c)) return JsonType::Number;
				throw DataError("JSON: Unexpected character");
			}

			/// @brief Parse array of objects into Dataset
			Dataset ParseArrayOfObjects() {
				Dataset dataset;

				if (Get() != '[')
					throw DataError("JSON: Expected '[' at root");

				std::map<std::string, std::vector<std::string>> columnData;
				std::vector<std::string> columnOrder;

				// Parse each object
				while (Peek() != ']') {
					if (Peek() == ',') Get();  // Skip comma

					if (Get() != '{')
						throw DataError("JSON: Expected '{' for object");

					// Parse object key-value pairs
					while (Peek() != '}') {
						if (Peek() == ',') Get();

						std::string key = ParseString();
						if (Get() != ':')
							throw DataError("JSON: Expected ':' after key");

						// Track column order
						if (columnData.find(key) == columnData.end()) {
							columnOrder.push_back(key);
						}

						// Parse value as string
						std::string value;
						switch (PeekType()) {
							case JsonType::String:
								value = ParseString();
								break;
							case JsonType::Number:
								value = ParseNumber();
								break;
							case JsonType::Bool:
								value = ParseBool() ? "true" : "false";
								break;
							case JsonType::Null:
								ParseNull();
								value = "";
								break;
							default:
								throw DataError("JSON: Unsupported value type (nested objects/arrays not supported)");
						}

						columnData[key].push_back(value);
					}

					if (Get() != '}')
						throw DataError("JSON: Expected '}'");

					// Fill missing columns with empty for this row
					size_t maxRows = 0;
					for (const auto& col : columnOrder) {
						maxRows = std::max(maxRows, columnData[col].size());
					}
					for (const auto& col : columnOrder) {
						while (columnData[col].size() < maxRows)
							columnData[col].push_back("");
					}
				}

				if (Get() != ']')
					throw DataError("JSON: Expected ']'");

				// Build dataset
				if (!columnOrder.empty()) {
					dataset.rowCount = columnData[columnOrder[0]].size();
					dataset.columns.reserve(columnOrder.size());

					for (const auto& colName : columnOrder) {
						DataColumn col;
						col.name = colName;
						col.type = InferColumnType(columnData[colName]);
						col.missingMask.resize(dataset.rowCount, false);

						// Parse values
						const auto& values = columnData[colName];
						for (size_t i = 0; i < values.size(); ++i) {
							Real realVal = 0.0;
							int intVal = 0;
							bool boolVal = false;
							std::string strVal;

							bool parsed = ParseValue(values[i], col.type, realVal, intVal, boolVal, strVal);
							col.missingMask[i] = !parsed;

							switch (col.type) {
								case ColumnType::REAL:
									if (i == 0) col.realData = Vector<Real>(values.size());
									col.realData[i] = parsed ? realVal : std::numeric_limits<Real>::quiet_NaN();
									break;
								case ColumnType::INT:
									if (i == 0) col.intData = Vector<int>(values.size());
									col.intData[i] = parsed ? intVal : 0;
									break;
								case ColumnType::BOOL:
									if (i == 0) col.boolData.resize(values.size());
									col.boolData[i] = parsed ? boolVal : false;
									break;
								default:
									if (i == 0) col.stringData.resize(values.size());
									col.stringData[i] = parsed ? strVal : "";
									break;
							}
						}

						dataset.columns.push_back(std::move(col));
					}
				}

				return dataset;
			}
		};

		/// @brief Load dataset from JSON file (array of objects format)
		/// @param filename Path to JSON file
		/// @return Loaded dataset
		/// @throws DataError on file/parse errors
		inline Dataset LoadJSON(const std::string& filename) {
			std::ifstream file(filename);
			if (!file.is_open())
				throw DataError("LoadJSON: Cannot open file '" + filename + "'");

			std::stringstream buffer;
			buffer << file.rdbuf();
			std::string content = buffer.str();

			// Remove BOM if present
			content = RemoveBOM(content);

			JsonParser parser(content);
			Dataset dataset = parser.ParseArrayOfObjects();
			dataset.name = filename;

			return dataset;
		}

		/////////////////////////////////////////////////////////////////////////////////////
		///                              STRING LOADING                                    ///
		/////////////////////////////////////////////////////////////////////////////////////

		/// @brief Load dataset from string content
		/// @param content String containing data
		/// @param format Data format
		/// @param hasHeader If true, first row is column names
		/// @return Loaded dataset
		inline Dataset LoadFromString(const std::string& content, DataFormat format,
		                              bool hasHeader = true) {
			// Create a temporary stringstream
			std::string cleanContent = RemoveBOM(content);

			if (format == DataFormat::JSON) {
				JsonParser parser(cleanContent);
				return parser.ParseArrayOfObjects();
			}

			// CSV/TSV parsing
			char delimiter = (format == DataFormat::TSV) ? '\t' : ',';

			Dataset dataset;
			std::vector<std::vector<std::string>> allData;
			std::istringstream iss(cleanContent);
			std::string line;

			while (std::getline(iss, line)) {
				if (Trim(line).empty()) continue;
				if (!line.empty() && line.back() == '\r')
					line.pop_back();
				allData.push_back(SplitLine(line, delimiter));
			}

			if (allData.empty())
				throw DataError("LoadFromString: Content is empty");

			size_t numCols = allData[0].size();
			size_t dataStartRow = hasHeader ? 1 : 0;

			dataset.columns.resize(numCols);

			for (size_t c = 0; c < numCols; ++c) {
				if (hasHeader) {
					dataset.columns[c].name = Trim(allData[0][c]);
					if (dataset.columns[c].name.empty())
						dataset.columns[c].name = "Column" + std::to_string(c);
				}
				else {
					dataset.columns[c].name = "Column" + std::to_string(c);
				}

				// Type inference
				std::vector<std::string> colValues;
				for (size_t r = dataStartRow; r < allData.size(); ++r) {
					if (c < allData[r].size())
						colValues.push_back(allData[r][c]);
				}
				dataset.columns[c].type = InferColumnType(colValues);
			}

			dataset.rowCount = allData.size() - dataStartRow;

			// Parse data (same as LoadCSV)
			for (size_t c = 0; c < numCols; ++c) {
				DataColumn& col = dataset.columns[c];
				col.missingMask.resize(dataset.rowCount, false);

				switch (col.type) {
					case ColumnType::REAL:
						col.realData = Vector<Real>(dataset.rowCount);
						break;
					case ColumnType::INT:
						col.intData = Vector<int>(dataset.rowCount);
						break;
					case ColumnType::BOOL:
						col.boolData.resize(dataset.rowCount);
						break;
					default:
						col.stringData.resize(dataset.rowCount);
						break;
				}
			}

			for (size_t r = dataStartRow; r < allData.size(); ++r) {
				size_t rowIdx = r - dataStartRow;
				const auto& row = allData[r];

				for (size_t c = 0; c < numCols; ++c) {
					DataColumn& col = dataset.columns[c];
					std::string value = (c < row.size()) ? row[c] : "";

					Real realVal = 0.0;
					int intVal = 0;
					bool boolVal = false;
					std::string strVal;

					bool parsed = ParseValue(value, col.type, realVal, intVal, boolVal, strVal);
					col.missingMask[rowIdx] = !parsed;

					switch (col.type) {
						case ColumnType::REAL:
							col.realData[rowIdx] = parsed ? realVal : std::numeric_limits<Real>::quiet_NaN();
							break;
						case ColumnType::INT:
							col.intData[rowIdx] = parsed ? intVal : 0;
							break;
						case ColumnType::BOOL:
							col.boolData[rowIdx] = parsed ? boolVal : false;
							break;
						default:
							col.stringData[rowIdx] = parsed ? strVal : "";
							break;
					}
				}
			}

			return dataset;
		}

		/////////////////////////////////////////////////////////////////////////////////////
		///                              AUTO-FORMAT LOADING                               ///
		/////////////////////////////////////////////////////////////////////////////////////

		/// @brief Detect format from file extension
		inline DataFormat DetectFormat(const std::string& filename) {
			size_t dotPos = filename.rfind('.');
			if (dotPos == std::string::npos)
				return DataFormat::CSV;  // Default

			std::string ext = filename.substr(dotPos + 1);
			std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);

			if (ext == "csv") return DataFormat::CSV;
			if (ext == "tsv" || ext == "tab") return DataFormat::TSV;
			if (ext == "json") return DataFormat::JSON;

			return DataFormat::CSV;  // Default
		}

		/// @brief Load dataset with auto-detected format
		/// @param filename Path to file
		/// @param hasHeader If true, first row is column names (for CSV/TSV)
		/// @return Loaded dataset
		inline Dataset LoadData(const std::string& filename, bool hasHeader = true) {
			DataFormat format = DetectFormat(filename);

			switch (format) {
				case DataFormat::JSON:
					return LoadJSON(filename);
				case DataFormat::TSV:
					return LoadTSV(filename, hasHeader);
				case DataFormat::CSV:
				default:
					return LoadCSV(filename, hasHeader);
			}
		}

	}  // namespace Data
}  // namespace MML

#endif  // MML_DATA_LOADER_H
