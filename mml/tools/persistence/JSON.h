///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        persistence/JSON.h                                                  ///
///  Description: Minimal JSON value, parser, and writer for serializer schemas       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_PERSISTENCE_JSON_H
#define MML_PERSISTENCE_JSON_H

#include <mml/tools/persistence/PersistenceBase.h>

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace MML
{
namespace Persistence
{
	inline constexpr std::size_t MaxJsonNestingDepth = 256;

	enum class JsonValueType
	{
		Null,
		Bool,
		Number,
		String,
		Array,
		Object
	};

	struct JsonValue
	{
		using Array = std::vector<JsonValue>;
		using Object = std::map<std::string, JsonValue>;

		JsonValueType type = JsonValueType::Null;
		bool bool_value = false;
		double number_value = 0.0;
		std::string string_value;
		Array array_value;
		Object object_value;
		std::string number_text;
		std::vector<std::string> object_key_order;

		static JsonValue Null()
		{
			return JsonValue{};
		}

		static JsonValue Bool(bool value)
		{
			JsonValue json;
			json.type = JsonValueType::Bool;
			json.bool_value = value;
			return json;
		}

		static JsonValue Number(double value, std::string text = {})
		{
			JsonValue json;
			json.type = JsonValueType::Number;
			json.number_value = value;
			json.number_text = std::move(text);
			return json;
		}

		static JsonValue String(std::string value)
		{
			JsonValue json;
			json.type = JsonValueType::String;
			json.string_value = std::move(value);
			return json;
		}

		static JsonValue ArrayValue(Array value = {})
		{
			JsonValue json;
			json.type = JsonValueType::Array;
			json.array_value = std::move(value);
			return json;
		}

		static JsonValue ObjectValue(Object value = {})
		{
			JsonValue json;
			json.type = JsonValueType::Object;
			json.object_value = std::move(value);
			return json;
		}
	};

	struct JsonParseResult
	{
		SerializeResult result;
		JsonValue value;
	};

	namespace Detail
	{
		class JsonParser
		{
		public:
			JsonParser(const std::string& input) : _input(input) {}

			JsonParseResult Parse()
			{
				JsonParseResult parsed;
				parsed.value = ParseValue(parsed.result);
				if (!parsed.result.success)
					return parsed;

				SkipWhitespace();
				if (_pos != _input.size())
				{
					parsed.result = Fail("Unexpected trailing input");
					return parsed;
				}

				parsed.result = SerializeSuccess();
				return parsed;
			}

		private:
			const std::string& _input;
			std::size_t _pos = 0;

			SerializeResult Fail(const std::string& message) const
			{
				return SerializeFailure(SerializeError::MALFORMED_INPUT,
				                        message + " at byte " + std::to_string(_pos));
			}

			void SkipWhitespace()
			{
				while (_pos < _input.size())
				{
					const unsigned char c = static_cast<unsigned char>(_input[_pos]);
					if (c != ' ' && c != '\n' && c != '\r' && c != '\t')
						break;
					++_pos;
				}
			}

			bool Consume(char expected)
			{
				SkipWhitespace();
				if (_pos < _input.size() && _input[_pos] == expected)
				{
					++_pos;
					return true;
				}
				return false;
			}

			JsonValue ParseValue(SerializeResult& result, std::size_t depth = 0)
			{
				SkipWhitespace();
				if (_pos >= _input.size())
				{
					result = Fail("Expected JSON value");
					return JsonValue::Null();
				}

				const char c = _input[_pos];
				if (c == 'n') return ParseLiteral("null", JsonValue::Null(), result);
				if (c == 't') return ParseLiteral("true", JsonValue::Bool(true), result);
				if (c == 'f') return ParseLiteral("false", JsonValue::Bool(false), result);
				if (c == '"') return ParseString(result);
				if ((c == '[' || c == '{') && depth >= MaxJsonNestingDepth)
				{
					result = Fail("Maximum JSON nesting depth exceeded");
					return JsonValue::Null();
				}
				if (c == '[') return ParseArray(result, depth);
				if (c == '{') return ParseObject(result, depth);
				if (c == '-' || (c >= '0' && c <= '9')) return ParseNumber(result);

				result = Fail("Unexpected character while reading JSON value");
				return JsonValue::Null();
			}

			JsonValue ParseLiteral(const char* literal, JsonValue value, SerializeResult& result)
			{
				const std::size_t start = _pos;
				for (const char* p = literal; *p != '\0'; ++p)
				{
					if (_pos >= _input.size() || _input[_pos] != *p)
					{
						_pos = start;
						result = Fail(std::string("Expected literal ") + literal);
						return JsonValue::Null();
					}
					++_pos;
				}
				result = SerializeSuccess();
				return value;
			}

			JsonValue ParseString(SerializeResult& result)
			{
				if (!Consume('"'))
				{
					result = Fail("Expected string opening quote");
					return JsonValue::Null();
				}

				std::string value;
				while (_pos < _input.size())
				{
					char c = _input[_pos++];
					if (c == '"')
					{
						result = SerializeSuccess();
						return JsonValue::String(value);
					}
					if (static_cast<unsigned char>(c) < 0x20)
					{
						result = Fail("Unescaped control character in string");
						return JsonValue::Null();
					}
					if (c != '\\')
					{
						value += c;
						continue;
					}

					if (_pos >= _input.size())
					{
						result = Fail("Incomplete escape sequence");
						return JsonValue::Null();
					}

					const char esc = _input[_pos++];
					switch (esc)
					{
					case '"': value += '"'; break;
					case '\\': value += '\\'; break;
					case '/': value += '/'; break;
					case 'b': value += '\b'; break;
					case 'f': value += '\f'; break;
					case 'n': value += '\n'; break;
					case 'r': value += '\r'; break;
					case 't': value += '\t'; break;
					case 'u':
						if (!AppendUnicodeEscape(value, result))
							return JsonValue::Null();
						break;
					default:
						result = Fail("Invalid escape sequence");
						return JsonValue::Null();
					}
				}

				result = Fail("Unterminated string");
				return JsonValue::Null();
			}

			bool AppendUnicodeEscape(std::string& value, SerializeResult& result)
			{
				unsigned codePoint = 0;
				if (!ParseUnicodeCodeUnit(codePoint, result))
					return false;

				if (codePoint >= 0xD800 && codePoint <= 0xDBFF)
				{
					if (_pos + 2 > _input.size() || _input[_pos] != '\\' || _input[_pos + 1] != 'u')
					{
						result = Fail("Expected low surrogate after high surrogate");
						return false;
					}
					_pos += 2;
					unsigned lowSurrogate = 0;
					if (!ParseUnicodeCodeUnit(lowSurrogate, result))
						return false;
					if (lowSurrogate < 0xDC00 || lowSurrogate > 0xDFFF)
					{
						result = Fail("Invalid low surrogate");
						return false;
					}
					codePoint = 0x10000 + ((codePoint - 0xD800) << 10) + (lowSurrogate - 0xDC00);
				}
				else if (codePoint >= 0xDC00 && codePoint <= 0xDFFF)
				{
					result = Fail("Unexpected low surrogate");
					return false;
				}

				AppendUtf8(value, codePoint);
				result = SerializeSuccess();
				return true;
			}

			bool ParseUnicodeCodeUnit(unsigned& codeUnit, SerializeResult& result)
			{
				if (_pos + 4 > _input.size())
				{
					result = Fail("Incomplete unicode escape");
					return false;
				}

				codeUnit = 0;
				for (int i = 0; i < 4; ++i)
				{
					const char c = _input[_pos++];
					codeUnit <<= 4;
					if (c >= '0' && c <= '9') codeUnit += static_cast<unsigned>(c - '0');
					else if (c >= 'a' && c <= 'f') codeUnit += static_cast<unsigned>(c - 'a' + 10);
					else if (c >= 'A' && c <= 'F') codeUnit += static_cast<unsigned>(c - 'A' + 10);
					else
					{
						result = Fail("Invalid unicode escape");
						return false;
					}
				}
				return true;
			}

			static void AppendUtf8(std::string& value, unsigned codePoint)
			{
				if (codePoint <= 0x7F)
				{
					value += static_cast<char>(codePoint);
				}
				else if (codePoint <= 0x7FF)
				{
					value += static_cast<char>(0xC0 | (codePoint >> 6));
					value += static_cast<char>(0x80 | (codePoint & 0x3F));
				}
				else if (codePoint <= 0xFFFF)
				{
					value += static_cast<char>(0xE0 | (codePoint >> 12));
					value += static_cast<char>(0x80 | ((codePoint >> 6) & 0x3F));
					value += static_cast<char>(0x80 | (codePoint & 0x3F));
				}
				else
				{
					value += static_cast<char>(0xF0 | (codePoint >> 18));
					value += static_cast<char>(0x80 | ((codePoint >> 12) & 0x3F));
					value += static_cast<char>(0x80 | ((codePoint >> 6) & 0x3F));
					value += static_cast<char>(0x80 | (codePoint & 0x3F));
				}
			}

			JsonValue ParseNumber(SerializeResult& result)
			{
				const std::size_t start = _pos;
				if (_input[_pos] == '-') ++_pos;

				if (_pos >= _input.size())
				{
					result = Fail("Incomplete number");
					return JsonValue::Null();
				}

				if (_input[_pos] == '0')
					++_pos;
				else if (_input[_pos] >= '1' && _input[_pos] <= '9')
					while (_pos < _input.size() && _input[_pos] >= '0' && _input[_pos] <= '9') ++_pos;
				else
				{
					result = Fail("Invalid number");
					return JsonValue::Null();
				}

				if (_pos < _input.size() && _input[_pos] == '.')
				{
					++_pos;
					if (_pos >= _input.size() || _input[_pos] < '0' || _input[_pos] > '9')
					{
						result = Fail("Expected digit after decimal point");
						return JsonValue::Null();
					}
					while (_pos < _input.size() && _input[_pos] >= '0' && _input[_pos] <= '9') ++_pos;
				}

				if (_pos < _input.size() && (_input[_pos] == 'e' || _input[_pos] == 'E'))
				{
					++_pos;
					if (_pos < _input.size() && (_input[_pos] == '+' || _input[_pos] == '-')) ++_pos;
					if (_pos >= _input.size() || _input[_pos] < '0' || _input[_pos] > '9')
					{
						result = Fail("Expected exponent digits");
						return JsonValue::Null();
					}
					while (_pos < _input.size() && _input[_pos] >= '0' && _input[_pos] <= '9') ++_pos;
				}

				const std::string text = _input.substr(start, _pos - start);
				errno = 0;
				char* endPtr = nullptr;
				const double value = std::strtod(text.c_str(), &endPtr);
				if (errno == ERANGE || endPtr == text.c_str() || *endPtr != '\0')
				{
					result = Fail("Invalid or out-of-range number");
					return JsonValue::Null();
				}

				result = SerializeSuccess();
				return JsonValue::Number(value, text);
			}

			JsonValue ParseArray(SerializeResult& result, std::size_t depth)
			{
				Consume('[');
				JsonValue::Array values;
				SkipWhitespace();
				if (Consume(']'))
				{
					result = SerializeSuccess();
					return JsonValue::ArrayValue(values);
				}

				while (true)
				{
					values.push_back(ParseValue(result, depth + 1));
					if (!result.success) return JsonValue::Null();

					if (Consume(']'))
					{
						result = SerializeSuccess();
						return JsonValue::ArrayValue(values);
					}
					if (!Consume(','))
					{
						result = Fail("Expected ',' or ']' in array");
						return JsonValue::Null();
					}
				}
			}

			JsonValue ParseObject(SerializeResult& result, std::size_t depth)
			{
				Consume('{');
				JsonValue::Object object;
				std::vector<std::string> keyOrder;
				SkipWhitespace();
				if (Consume('}'))
				{
					result = SerializeSuccess();
					return JsonValue::ObjectValue(object);
				}

				while (true)
				{
					SkipWhitespace();
					if (_pos >= _input.size() || _input[_pos] != '"')
					{
						result = Fail("Expected object key");
						return JsonValue::Null();
					}

					JsonValue key = ParseString(result);
					if (!result.success) return JsonValue::Null();
					if (!Consume(':'))
					{
						result = Fail("Expected ':' after object key");
						return JsonValue::Null();
					}

					const bool isNewKey = object.find(key.string_value) == object.end();
					object[key.string_value] = ParseValue(result, depth + 1);
					if (!result.success) return JsonValue::Null();
					if (isNewKey)
						keyOrder.push_back(key.string_value);

					if (Consume('}'))
					{
						result = SerializeSuccess();
						JsonValue value = JsonValue::ObjectValue(object);
						value.object_key_order = std::move(keyOrder);
						return value;
					}
					if (!Consume(','))
					{
						result = Fail("Expected ',' or '}' in object");
						return JsonValue::Null();
					}
				}
			}
		};

		inline void WriteIndent(std::ostream& out, int indent)
		{
			for (int i = 0; i < indent; ++i) out << ' ';
		}

		inline SerializeResult ValidateJsonForWrite(const JsonValue& value)
		{
			switch (value.type)
			{
			case JsonValueType::Number:
				if (!std::isfinite(value.number_value))
					return SerializeFailure(SerializeError::UNSUPPORTED_SCALAR,
					                        "JSON writer cannot emit NaN or infinity as standard JSON numbers");
				break;
			case JsonValueType::Array:
				for (const JsonValue& item : value.array_value)
				{
					SerializeResult result = ValidateJsonForWrite(item);
					if (!result.success) return result;
				}
				break;
			case JsonValueType::Object:
				for (const auto& entry : value.object_value)
				{
					SerializeResult result = ValidateJsonForWrite(entry.second);
					if (!result.success) return result;
				}
				break;
			default:
				break;
			}
			return SerializeSuccess();
		}
	}

	inline JsonParseResult ParseJson(const std::string& input)
	{
		return Detail::JsonParser(input).Parse();
	}

	inline std::string EscapeJsonString(const std::string& input)
	{
		std::ostringstream out;
		for (char c : input)
		{
			switch (c)
			{
			case '"': out << "\\\""; break;
			case '\\': out << "\\\\"; break;
			case '\b': out << "\\b"; break;
			case '\f': out << "\\f"; break;
			case '\n': out << "\\n"; break;
			case '\r': out << "\\r"; break;
			case '\t': out << "\\t"; break;
			default:
				if (static_cast<unsigned char>(c) < 0x20)
				{
					out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
					    << static_cast<int>(static_cast<unsigned char>(c))
					    << std::dec << std::setfill(' ');
				}
				else
				{
					out << c;
				}
			}
		}
		return out.str();
	}

	inline SerializeResult WriteJson(std::ostream& out, const JsonValue& value, const SaveOptions& options = {})
	{
		SerializeResult validation = Detail::ValidateJsonForWrite(value);
		if (!validation.success)
			return validation;

		const int indentStep = options.pretty_json ? 2 : 0;
		std::function<void(const JsonValue&, int)> writeValue;
		writeValue = [&](const JsonValue& json, int indent)
		{
			switch (json.type)
			{
			case JsonValueType::Null:
				out << "null";
				break;
			case JsonValueType::Bool:
				out << (json.bool_value ? "true" : "false");
				break;
			case JsonValueType::Number:
				out << std::setprecision(options.precision) << json.number_value;
				break;
			case JsonValueType::String:
				out << '"' << EscapeJsonString(json.string_value) << '"';
				break;
			case JsonValueType::Array:
				out << '[';
				for (std::size_t i = 0; i < json.array_value.size(); ++i)
				{
					if (i > 0) out << ',';
					if (options.pretty_json) { out << '\n'; Detail::WriteIndent(out, indent + indentStep); }
					writeValue(json.array_value[i], indent + indentStep);
				}
				if (options.pretty_json && !json.array_value.empty()) { out << '\n'; Detail::WriteIndent(out, indent); }
				out << ']';
				break;
			case JsonValueType::Object:
				out << '{';
				{
					std::size_t i = 0;
					for (const auto& entry : json.object_value)
					{
						if (i++ > 0) out << ',';
						if (options.pretty_json) { out << '\n'; Detail::WriteIndent(out, indent + indentStep); }
						out << '"' << EscapeJsonString(entry.first) << '"' << ':';
						if (options.pretty_json) out << ' ';
						writeValue(entry.second, indent + indentStep);
					}
				}
				if (options.pretty_json && !json.object_value.empty()) { out << '\n'; Detail::WriteIndent(out, indent); }
				out << '}';
				break;
			}
		};

		writeValue(value, 0);
		if (out.fail())
			return SerializeFailure(SerializeError::WRITE_FAILED, "Failed to write JSON value");
		return SerializeSuccess();
	}

	inline std::string ToJsonString(const JsonValue& value, const SaveOptions& options = {})
	{
		std::ostringstream out;
		WriteJson(out, value, options);
		return out.str();
	}

} // namespace Persistence
} // namespace MML

#endif // MML_PERSISTENCE_JSON_H