///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DataLoaderCSV.h                                                     ///
///  Description: CSV and TSV file loading functions                                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DATA_LOADER_CSV_H
#define MML_DATA_LOADER_CSV_H

#include "tools/data_loader/DataLoaderTypes.h"
#include "tools/data_loader/DataLoaderParsing.h"

#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace MML {
	namespace Data {

		/// @brief Options for CSV/TSV loading with schema validation
		struct CSVLoadOptions {
			bool hasHeader = true;          ///< First row is column names
			char delimiter = ',';           ///< Field delimiter
			bool inferTypes = true;         ///< Auto-detect column types
			bool strictSchema = false;      ///< Throw on row width mismatch
			std::vector<std::string>* warnings = nullptr;  ///< Collect non-fatal issues (e.g., width mismatches)
		};

		/////////////////////////////////////////////////////////////////////////////////////
		///                              CSV/TSV LOADING                                   ///
		/////////////////////////////////////////////////////////////////////////////////////
		/// @brief Load dataset from CSV file with full options
		/// @param filename Path to CSV file
		/// @param opts Loading options (schema validation, type inference, etc.)
		/// @return Loaded dataset
		/// @throws DataError on file/parse errors or schema violations (if strictSchema)
		inline Dataset LoadCSV(const std::string& filename, const CSVLoadOptions& opts) {
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

				std::vector<std::string> fields = SplitLine(line, opts.delimiter);
				allData.push_back(fields);
			}

			if (allData.empty())
				throw DataError("LoadCSV: File is empty");

			// Determine column count
			size_t numCols = allData[0].size();
			size_t dataStartRow = opts.hasHeader ? 1 : 0;

			// Setup columns
			dataset.columns.resize(numCols);

			for (size_t c = 0; c < numCols; ++c) {
				if (opts.hasHeader) {
					dataset.columns[c].name = Trim(allData[0][c]);
					if (dataset.columns[c].name.empty())
						dataset.columns[c].name = "Column" + std::to_string(c);
				}
				else {
					dataset.columns[c].name = "Column" + std::to_string(c);
				}
			}

			// Schema validation: check row widths
			if (opts.strictSchema || opts.warnings) {
				for (size_t r = dataStartRow; r < allData.size(); ++r) {
					if (allData[r].size() != numCols) {
						std::string msg = "Row " + std::to_string(r + 1) + " has " +
						                  std::to_string(allData[r].size()) + " fields, expected " +
						                  std::to_string(numCols);
						if (opts.strictSchema)
							throw DataError("LoadCSV: Schema mismatch - " + msg);
						if (opts.warnings)
							opts.warnings->push_back(msg);
					}
				}
			}

			// Collect string values for type inference
			if (opts.inferTypes) {
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

					auto result = ParseValue(value, col.type, realVal, intVal, boolVal, strVal);
					bool parsed = (result == ParseResult::Parsed);
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

		/// @brief Load dataset from CSV file (convenience overload)
		/// @param filename Path to CSV file
		/// @param hasHeader If true, first row is column names
		/// @param delimiter Field delimiter (default ',')
		/// @param inferTypes If true, automatically infer column types
		/// @return Loaded dataset
		/// @throws DataError on file/parse errors
		inline Dataset LoadCSV(const std::string& filename, bool hasHeader = true,
		                       char delimiter = ',', bool inferTypes = true) {
			CSVLoadOptions opts;
			opts.hasHeader = hasHeader;
			opts.delimiter = delimiter;
			opts.inferTypes = inferTypes;
			return LoadCSV(filename, opts);
		}

		/// @brief Load dataset from TSV file
		inline Dataset LoadTSV(const std::string& filename, bool hasHeader = true, bool inferTypes = true) {
			return LoadCSV(filename, hasHeader, '\t', inferTypes);
		}

		/////////////////////////////////////////////////////////////////////////////////////
		///                      CSV/TSV LOADING (SAFE VARIANTS)                          ///
		/////////////////////////////////////////////////////////////////////////////////////
		/// @brief Load dataset from CSV file with result-based error handling
		/// @param filename Path to CSV file
		/// @param hasHeader If true, first row is column names
		/// @param delimiter Field delimiter (default ',')
		/// @param inferTypes If true, automatically infer column types
		/// @return LoadResult with success status, error message, and loaded dataset
		/// @details This function returns errors instead of throwing exceptions.
		/// Use this for better error composition and when exceptions are not desired.
		/// For backward compatibility, use LoadCSV() which throws on error.
		inline LoadResult LoadCSVSafe(const std::string& filename, bool hasHeader = true,
		                              char delimiter = ',', bool inferTypes = true) {
			try {
				Dataset dataset = LoadCSV(filename, hasHeader, delimiter, inferTypes);
				return LoadResult::Success(dataset);
			}
			catch (const DataError& e) {
				return LoadResult::Failure(std::string("CSV load failed: ") + e.what());
			}
			catch (const std::exception& e) {
				return LoadResult::Failure(std::string("Unexpected error: ") + e.what());
			}
			catch (...) {
				return LoadResult::Failure("Unknown error occurred during CSV loading");
			}
		}

		/// @brief Load dataset from TSV file with result-based error handling
		/// @param filename Path to TSV file
		/// @param hasHeader If true, first row is column names
		/// @param inferTypes If true, automatically infer column types
		/// @return LoadResult with success status, error message, and loaded dataset
		/// @details This function returns errors instead of throwing exceptions.
		inline LoadResult LoadTSVSafe(const std::string& filename, bool hasHeader = true, bool inferTypes = true) {
			return LoadCSVSafe(filename, hasHeader, '\t', inferTypes);
		}

		/////////////////////////////////////////////////////////////////////////////////////
		///                              STRING LOADING (CSV/TSV)                          ///
		/////////////////////////////////////////////////////////////////////////////////////
		/// @brief Load dataset from CSV/TSV string content
		/// @param content String containing data
		/// @param delimiter Field delimiter (',' for CSV, '\t' for TSV)
		/// @param hasHeader If true, first row is column names
		/// @return Loaded dataset
		inline Dataset LoadFromCSVString(const std::string& content, char delimiter = ',',
		                                  bool hasHeader = true) {
			std::string cleanContent = RemoveBOM(content);

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
				throw DataError("LoadFromCSVString: Content is empty");

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

			// Parse data
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
					case ColumnType::STRING:
						col.stringData.resize(dataset.rowCount);
						break;
					case ColumnType::DATE:
					case ColumnType::TIME:
						col.stringData.resize(dataset.rowCount);
						break;
					case ColumnType::DATETIME:
						col.dateData.resize(dataset.rowCount);
						col.timeData.resize(dataset.rowCount);
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

					auto result = ParseValue(value, col.type, realVal, intVal, boolVal, strVal);
					bool parsed = (result == ParseResult::Parsed);
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
							col.stringData[rowIdx] = parsed ? strVal : "";
							break;
					}
				}
			}

			return dataset;
		}

	}  // namespace Data
}  // namespace MML

#endif  // MML_DATA_LOADER_CSV_H
