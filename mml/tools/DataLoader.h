///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DataLoader.h                                                        ///
///  Description: Umbrella header for multi-format dataset loading infrastructure    ///
///               Supports CSV, TSV, JSON with automatic type inference              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
///
/// This is an umbrella header that includes all data loading functionality.
/// For finer-grained includes, use the individual headers:
///
///   - DataLoaderTypes.h   : DataFormat, ColumnType, DataColumn, Dataset, LoadResult
///   - DataLoaderParsing.h : InferColumnType, ParseValue, SplitLine, Trim, RemoveBOM
///   - DataLoaderCSV.h     : LoadCSV, LoadTSV, LoadCSVSafe, LoadTSVSafe
///   - DataLoaderJSON.h    : LoadJSON, LoadJSONSafe, JsonParser
///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DATA_LOADER_H
#define MML_DATA_LOADER_H

// Core types: DataFormat, ColumnType, DataColumn, Dataset, LoadResult
#include "tools/data_loader/DataLoaderTypes.h"

// Parsing utilities: type inference, value parsing, string utilities
#include "tools/data_loader/DataLoaderParsing.h"

// CSV/TSV loading
#include "tools/data_loader/DataLoaderCSV.h"

// JSON loading
#include "tools/data_loader/DataLoaderJSON.h"

#include <algorithm>

namespace MML {
	namespace Data {

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

		/// @brief Load dataset from string content
		/// @param content String containing data
		/// @param format Data format
		/// @param hasHeader If true, first row is column names
		/// @return Loaded dataset
		inline Dataset LoadFromString(const std::string& content, DataFormat format,
		                              bool hasHeader = true) {
			switch (format) {
				case DataFormat::JSON:
					return LoadFromJSONString(content);
				case DataFormat::TSV:
					return LoadFromCSVString(content, '\t', hasHeader);
				case DataFormat::CSV:
				default:
					return LoadFromCSVString(content, ',', hasHeader);
			}
		}

	}  // namespace Data
}  // namespace MML

#endif  // MML_DATA_LOADER_H
