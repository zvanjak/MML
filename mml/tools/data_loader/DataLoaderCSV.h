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

#include <mml/tools/data_loader/DataLoaderTypes.h>
#include <mml/tools/data_loader/DataLoaderParsing.h>

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
		inline IOResult LoadCSV(const std::string& filename, Dataset& output, const CSVLoadOptions& opts = {}) {
			std::ifstream file(filename);
			if (!file.is_open())
				return SerializeFailure(IOError::FILE_NOT_OPENED, "Cannot open file").WithFile(filename);

			Dataset dataset;
			dataset.name = filename;

			std::vector<std::vector<std::string>> allData;
			std::string line;
			bool firstLine = true;
			std::int64_t lineNumber = 0;
			std::vector<std::int64_t> sourceLines;

			while (std::getline(file, line)) {
				++lineNumber;
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

				try {
					allData.push_back(SplitLine(line, opts.delimiter));
				} catch (const DataError& error) {
					return SerializeFailure(IOError::MALFORMED_INPUT, error.what()).WithFile(filename).At(lineNumber, IOOffsetKind::Line);
				}
				sourceLines.push_back(lineNumber);
			}
			if (file.bad() || !file.eof())
				return SerializeFailure(IOError::READ_FAILED, "CSV stream read failed").WithFile(filename).At(lineNumber + 1, IOOffsetKind::Line);

			if (allData.empty())
				return SerializeFailure(IOError::INVALID_FORMAT, "CSV file is empty").WithFile(filename).At(1, IOOffsetKind::Line);

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
							return SerializeFailure(IOError::SCHEMA_MISMATCH, msg).WithFile(filename).At(sourceLines[r], IOOffsetKind::Line);
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
						col.stringData.resize(dataset.rowCount);
						break;
					case ColumnType::DATE:
						col.dateData.resize(dataset.rowCount);
						break;
					case ColumnType::TIME:
						col.timeData.resize(dataset.rowCount);
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
							col.dateData[rowIdx] = parsed ? strVal : "";
							break;
						case ColumnType::TIME:
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

			output = std::move(dataset);
			return IOResult::Success(filename);
		}

		inline Dataset LoadCSV(const std::string& filename, const CSVLoadOptions& opts) {
			Dataset dataset;
			const auto result = LoadCSV(filename, dataset, opts);
			if (!result)
				throw DataError(result.message);
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
		/// Use LoadCSV() when exception-based error handling is preferred.
		inline IOResult LoadCSVSafe(const std::string& filename, Dataset& dataset, bool hasHeader = true, char delimiter = ',',
																bool inferTypes = true) {
			CSVLoadOptions opts;
			opts.hasHeader = hasHeader;
			opts.delimiter = delimiter;
			opts.inferTypes = inferTypes;
			return LoadCSV(filename, dataset, opts);
		}

		[[deprecated("use the IOResult overload with Dataset&")]]
		inline LegacyDatasetResult LoadCSVSafe(const std::string& filename, bool hasHeader = true, char delimiter = ',',
																					 bool inferTypes = true) {
			Dataset dataset;
			auto result = LoadCSVSafe(filename, dataset, hasHeader, delimiter, inferTypes);
			return LegacyDatasetResult::From(std::move(result), std::move(dataset));
		}

		/// @brief Load dataset from TSV file with result-based error handling
		/// @param filename Path to TSV file
		/// @param hasHeader If true, first row is column names
		/// @param inferTypes If true, automatically infer column types
		/// @return LoadResult with success status, error message, and loaded dataset
		/// @details This function returns errors instead of throwing exceptions.
		inline IOResult LoadTSVSafe(const std::string& filename, Dataset& dataset, bool hasHeader = true, bool inferTypes = true) {
			return LoadCSVSafe(filename, dataset, hasHeader, '\t', inferTypes);
		}

		[[deprecated("use the IOResult overload with Dataset&")]]
		inline LegacyDatasetResult LoadTSVSafe(const std::string& filename, bool hasHeader = true, bool inferTypes = true) {
			Dataset dataset;
			auto result = LoadTSVSafe(filename, dataset, hasHeader, inferTypes);
			return LegacyDatasetResult::From(std::move(result), std::move(dataset));
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
						col.dateData.resize(dataset.rowCount);
						break;
					case ColumnType::TIME:
						col.timeData.resize(dataset.rowCount);
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
							col.dateData[rowIdx] = parsed ? strVal : "";
							break;
						case ColumnType::TIME:
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
