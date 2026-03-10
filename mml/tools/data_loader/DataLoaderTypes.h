///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DataLoaderTypes.h                                                   ///
///  Description: Core types for data loading: DataFormat, ColumnType, DataColumn,   ///
///               Dataset, LoadResult                                                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DATA_LOADER_TYPES_H
#define MML_DATA_LOADER_TYPES_H

#include "MMLBase.h"
#include "MMLExceptions.h"
#include "base/Vector/Vector.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <string>
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
		///                              LOAD RESULT                                       ///
		/////////////////////////////////////////////////////////////////////////////////////
		/// @brief Result type for loading operations
		/// @details Provides success/error information and loaded dataset.
		/// Follows the pattern established by SerializeResult and VisualizerResult
		/// for consistent error handling across the tools layer.
		struct LoadResult {
			bool success;           ///< true if loading succeeded
			std::string errorMessage;  ///< Error description (empty on success)
			Dataset data;           ///< Loaded dataset (empty if failed)

			/// @brief Create success result with loaded dataset
			static LoadResult Success(const Dataset& dataset) {
				return {true, "", dataset};
			}

			/// @brief Create failure result with error message
			static LoadResult Failure(const std::string& message) {
				return {false, message, Dataset{}};
			}

			/// @brief Implicit conversion to bool for convenient checking
			/// @return true if loading succeeded
			operator bool() const { return success; }
		};

	}  // namespace Data
}  // namespace MML

#endif  // MML_DATA_LOADER_TYPES_H
