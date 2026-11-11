///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DataLoaderJSON.h                                                    ///
///  Description: JSON dataset loading through the shared persistence parser         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DATA_LOADER_JSON_H
#define MML_DATA_LOADER_JSON_H

#include <mml/tools/data_loader/DataLoaderTypes.h>
#include <mml/tools/data_loader/DataLoaderParsing.h>
#include <mml/tools/persistence/JSON.h>

#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace MML {
	namespace Data {

		/////////////////////////////////////////////////////////////////////////////////////
		///                              JSON LOADING                                      ///
		/////////////////////////////////////////////////////////////////////////////////////
		namespace Detail {
			inline std::string JsonScalarToString(const Persistence::JsonValue& value) {
				switch (value.type) {
					case Persistence::JsonValueType::Null:
						return "";
					case Persistence::JsonValueType::Bool:
						return value.bool_value ? "true" : "false";
					case Persistence::JsonValueType::Number: {
						if (!value.number_text.empty())
							return value.number_text;
						std::ostringstream text;
						TextIO::ApplyClassicFormat(text, std::numeric_limits<double>::max_digits10);
						text << value.number_value;
						return text.str();
					}
					case Persistence::JsonValueType::String:
						return value.string_value;
					default:
						throw DataError("JSON: Unsupported value type (nested objects/arrays not supported)");
				}
			}

			inline IOResult DatasetFromJson(const Persistence::JsonValue& root, Dataset& output) {
				if (root.type != Persistence::JsonValueType::Array)
					throw DataError("JSON: Expected array at root");

				std::map<std::string, std::vector<std::string>> columnData;
				std::vector<std::string> columnOrder;
				for (std::size_t rowIndex = 0; rowIndex < root.array_value.size(); ++rowIndex) {
					const Persistence::JsonValue& row = root.array_value[rowIndex];
					if (row.type != Persistence::JsonValueType::Object)
						throw DataError("JSON: Expected object for each array element");

					for (auto& column : columnData)
						column.second.push_back("");

					auto addValue = [&](const std::string& name, const Persistence::JsonValue& value) {
						auto [column, inserted] = columnData.try_emplace(name, rowIndex, "");
						if (inserted) {
							columnOrder.push_back(name);
							column->second.push_back(JsonScalarToString(value));
						}
						else
							column->second.back() = JsonScalarToString(value);
					};

					if (!row.object_key_order.empty()) {
						for (const std::string& name : row.object_key_order)
							addValue(name, row.object_value.at(name));
					}
					else {
						for (const auto& [name, value] : row.object_value)
							addValue(name, value);
					}
				}

				Dataset dataset;
				dataset.rowCount = root.array_value.size();
				dataset.columns.reserve(columnData.size());
				for (const std::string& colName : columnOrder) {
					const std::vector<std::string>& values = columnData.at(colName);
					DataColumn col;
					col.name = colName;
					col.type = InferColumnType(values);
					col.missingMask.resize(dataset.rowCount, false);

					// Parse values
					for (size_t i = 0; i < values.size(); ++i) {
						Real realVal = 0.0;
						int intVal = 0;
						bool boolVal = false;
						std::string strVal;

						auto result = ParseValue(values[i], col.type, realVal, intVal, boolVal, strVal);
						if (result == ParseResult::Invalid)
							return InvalidValueResult(values[i], col, i, dataset.columns.size());
						bool parsed = (result == ParseResult::Parsed);
						col.missingMask[i] = result == ParseResult::Missing;

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
								case ColumnType::STRING:
									if (i == 0) col.stringData.resize(values.size());
									col.stringData[i] = parsed ? strVal : "";
									break;
								case ColumnType::DATE:
									if (i == 0) col.dateData.resize(values.size());
									col.dateData[i] = parsed ? strVal : "";
									break;
								case ColumnType::TIME:
									if (i == 0) col.timeData.resize(values.size());
									col.timeData[i] = parsed ? strVal : "";
									break;
								case ColumnType::DATETIME:
									if (i == 0) {
										col.dateData.resize(values.size());
										col.timeData.resize(values.size());
									}
									if (parsed && !strVal.empty()) {
										size_t sep = strVal.find_first_of("T ");
										if (sep != std::string::npos) {
											col.dateData[i] = strVal.substr(0, sep);
											col.timeData[i] = strVal.substr(sep + 1);
										}
										else {
											col.dateData[i] = strVal;
											col.timeData[i] = "";
										}
									}
									break;
								default:
									if (i == 0) col.stringData.resize(values.size());
									col.stringData[i] = parsed ? strVal : "";
									break;
						}
					}

					dataset.columns.push_back(std::move(col));
				}
				output = std::move(dataset);
				return IOResult::Success();
			}

			inline Dataset ParseDatasetJson(const std::string& content) {
				Persistence::JsonParseResult parsed = Persistence::ParseJson(content);
				if (!parsed.result.success)
					throw DataError("JSON: " + parsed.result.message);
				Dataset dataset;
				const auto result = DatasetFromJson(parsed.value, dataset);
				if (!result)
					throw DataError(result.message);
				return dataset;
			}
		} // namespace Detail

		/// @brief Load dataset from JSON file (array of objects format)
		/// @param filename Path to JSON file
		/// @return Loaded dataset
		/// @throws DataError on file/parse errors
		inline IOResult LoadJSONSafe(const std::string& filename, Dataset& output);
		inline Dataset LoadJSON(const std::string& filename) {
			Dataset dataset;
			const auto result = LoadJSONSafe(filename, dataset);
			if (!result)
				throw DataError(result.message);
			return dataset;
		}

		/// @brief Load dataset from JSON file with result-based error handling
		/// @param filename Path to JSON file
		/// @return LoadResult with success status, error message, and loaded dataset
		/// @details This function returns errors instead of throwing exceptions.
		/// Use this for better error composition and when exceptions are not desired.
		/// Use LoadJSON() when exception-based error handling is preferred.
		inline IOResult LoadJSONSafe(const std::string& filename, Dataset& output) {
			std::ifstream file(filename, std::ios::binary);
			if (!file.is_open())
				return SerializeFailure(IOError::FILE_NOT_OPENED, "Cannot open JSON dataset").WithFile(filename);
			std::stringstream buffer;
			buffer << file.rdbuf();
			if (file.bad() || buffer.bad())
				return InputResult(file, SerializeFailure(IOError::READ_FAILED, "JSON dataset read failed"), filename);
			const std::string content = buffer.str();
			const std::string cleanContent = RemoveBOM(content);
			auto parsed = Persistence::ParseJson(cleanContent);
			if (!parsed.result) {
				if (parsed.result.offset >= 0)
					parsed.result.offset += static_cast<std::int64_t>(content.size() - cleanContent.size());
				return parsed.result.WithFile(filename);
			}
			try {
				Dataset dataset;
				const auto result = Detail::DatasetFromJson(parsed.value, dataset);
				if (!result)
					return result.WithFile(filename);
				dataset.name = filename;
				output = std::move(dataset);
				return IOResult::Success(filename);
			} catch (const DataError& error) {
				return SerializeFailure(IOError::SCHEMA_MISMATCH, error.what()).WithFile(filename);
			}
		}

		inline IOResult LoadJSON(const std::string& filename, Dataset& output) { return LoadJSONSafe(filename, output); }

		[[deprecated("use the IOResult overload with Dataset&")]]
		inline LegacyDatasetResult LoadJSONSafe(const std::string& filename) {
			Dataset dataset;
			auto result = LoadJSONSafe(filename, dataset);
			return LegacyDatasetResult::From(std::move(result), std::move(dataset));
		}

		/// @brief Load dataset from JSON string content
		/// @param content String containing JSON data
		/// @return Loaded dataset
		inline Dataset LoadFromJSONString(const std::string& content) {
			std::string cleanContent = RemoveBOM(content);
			return Detail::ParseDatasetJson(cleanContent);
		}

	}  // namespace Data
}  // namespace MML

#endif  // MML_DATA_LOADER_JSON_H
