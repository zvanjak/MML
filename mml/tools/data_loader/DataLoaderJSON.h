///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DataLoaderJSON.h                                                    ///
///  Description: JSON file loading with minimal JSON parser                          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DATA_LOADER_JSON_H
#define MML_DATA_LOADER_JSON_H

#include "tools/data_loader/DataLoaderTypes.h"
#include "tools/data_loader/DataLoaderParsing.h"

#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace MML {
	namespace Data {

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

							auto result = ParseValue(values[i], col.type, realVal, intVal, boolVal, strVal);
							bool parsed = (result == ParseResult::Parsed);
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

		/// @brief Load dataset from JSON file with result-based error handling
		/// @param filename Path to JSON file
		/// @return LoadResult with success status, error message, and loaded dataset
		/// @details This function returns errors instead of throwing exceptions.
		/// Use this for better error composition and when exceptions are not desired.
		/// For backward compatibility, use LoadJSON() which throws on error.
		inline LoadResult LoadJSONSafe(const std::string& filename) {
			try {
				Dataset dataset = LoadJSON(filename);
				return LoadResult::Success(dataset);
			}
			catch (const DataError& e) {
				return LoadResult::Failure(std::string("JSON load failed: ") + e.what());
			}
			catch (const std::exception& e) {
				return LoadResult::Failure(std::string("Unexpected error: ") + e.what());
			}
			catch (...) {
				return LoadResult::Failure("Unknown error occurred during JSON loading");
			}
		}

		/// @brief Load dataset from JSON string content
		/// @param content String containing JSON data
		/// @return Loaded dataset
		inline Dataset LoadFromJSONString(const std::string& content) {
			std::string cleanContent = RemoveBOM(content);
			JsonParser parser(cleanContent);
			return parser.ParseArrayOfObjects();
		}

	}  // namespace Data
}  // namespace MML

#endif  // MML_DATA_LOADER_JSON_H
