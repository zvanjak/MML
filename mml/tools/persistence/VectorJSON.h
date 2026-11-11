///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        persistence/VectorJSON.h                                            ///
///  Description: JSON persistence for Vector and VectorN                             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_PERSISTENCE_VECTOR_JSON_H
#define MML_PERSISTENCE_VECTOR_JSON_H

#include <mml/base/Vector/Vector.h>
#include <mml/base/Vector/VectorN.h>
#include <mml/tools/persistence/JSON.h>
#include <mml/tools/persistence/VectorBinary.h>

#include <fstream>
#include <cmath>
#include <limits>
#include <sstream>
#include <type_traits>

namespace MML
{
namespace Persistence
{
	namespace Detail
	{
		template<class Type>
		constexpr const char* JsonScalarName()
		{
			if constexpr (std::is_same_v<Type, Real>) return "Real";
			else if constexpr (std::is_same_v<Type, float>) return "float";
			else if constexpr (std::is_same_v<Type, double>) return "double";
			else if constexpr (std::is_same_v<Type, long double>) return "long_double";
			else return "unsupported";
		}

		template<class Type>
		constexpr const char* JsonConcreteScalarName()
		{
			if constexpr (std::is_same_v<Type, float>) return "float";
			else if constexpr (std::is_same_v<Type, double>) return "double";
			else if constexpr (std::is_same_v<Type, long double>) return "long_double";
			else return "unsupported";
		}

		template<class Type>
		JsonValue BuildVectorJsonHeader(const char* objectName, int objectVersion)
		{
			JsonValue::Object mml;
			mml["format"] = JsonValue::String("MML_JSON");
			mml["version"] = JsonValue::Number(1);
			mml["object"] = JsonValue::String(objectName);
			mml["object_version"] = JsonValue::Number(objectVersion);
			mml["scalar"] = JsonValue::String(JsonScalarName<Type>());
			mml["scalar_bytes"] = JsonValue::Number(static_cast<double>(sizeof(Type)));
			mml["scalar_encoding"] = JsonValue::String("ieee754");
			if constexpr (std::is_same_v<Type, Real>)
				mml["real_type"] = JsonValue::String(JsonConcreteScalarName<Real>());
			return JsonValue::ObjectValue(mml);
		}

		inline const JsonValue* FindObjectField(const JsonValue& object, const std::string& key)
		{
			if (object.type != JsonValueType::Object)
				return nullptr;
			auto it = object.object_value.find(key);
			return it == object.object_value.end() ? nullptr : &it->second;
		}

		template<class Type>
		IOResult ValidateVectorJsonRoot(const JsonValue& root, const char* objectName, int expectedSize, const LoadOptions& options) {
			if (root.type != JsonValueType::Object)
				return SerializeFailure(IOError::SCHEMA_MISMATCH, "Vector JSON root must be an object");

			const JsonValue* mml = FindObjectField(root, "mml");
			if (mml == nullptr || mml->type != JsonValueType::Object)
				return SerializeFailure(IOError::SCHEMA_MISMATCH, "Vector JSON missing mml header object");

			const JsonValue* format = FindObjectField(*mml, "format");
			const JsonValue* version = FindObjectField(*mml, "version");
			const JsonValue* object = FindObjectField(*mml, "object");
			const JsonValue* objectVersion = FindObjectField(*mml, "object_version");
			const JsonValue* scalar = FindObjectField(*mml, "scalar");
			const JsonValue* scalarBytes = FindObjectField(*mml, "scalar_bytes");
			const JsonValue* scalarEncoding = FindObjectField(*mml, "scalar_encoding");

			if (format == nullptr || format->type != JsonValueType::String || format->string_value != "MML_JSON")
				return SerializeFailure(IOError::INVALID_FORMAT, "Vector JSON format must be MML_JSON");
			if (version == nullptr || version->type != JsonValueType::Number || static_cast<int>(version->number_value) != 1)
				return SerializeFailure(IOError::UNSUPPORTED_VERSION, "Unsupported Vector JSON container version");
			if (object == nullptr || object->type != JsonValueType::String || object->string_value != objectName)
				return SerializeFailure(IOError::TYPE_MISMATCH, "Vector JSON object kind does not match target type");
			if (objectVersion == nullptr || objectVersion->type != JsonValueType::Number || static_cast<int>(objectVersion->number_value) != 1)
				return SerializeFailure(IOError::UNSUPPORTED_VERSION, "Unsupported Vector JSON object schema version");
			if (scalar == nullptr || scalar->type != JsonValueType::String || scalar->string_value != JsonScalarName<Type>())
				return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Vector JSON scalar type does not match target type");
			if (scalarBytes == nullptr || scalarBytes->type != JsonValueType::Number || static_cast<int>(scalarBytes->number_value) != static_cast<int>(sizeof(Type)))
				return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Vector JSON scalar byte size does not match target type");
			if (scalarEncoding == nullptr || scalarEncoding->type != JsonValueType::String || scalarEncoding->string_value != "ieee754")
				return SerializeFailure(IOError::UNSUPPORTED_SCALAR, "Vector JSON scalar encoding is unsupported");

			const JsonValue* shape = FindObjectField(root, "shape");
			const JsonValue* data = FindObjectField(root, "data");
			if (shape == nullptr || shape->type != JsonValueType::Array || shape->array_value.size() != 1 || shape->array_value[0].type != JsonValueType::Number)
				return SerializeFailure(IOError::SCHEMA_MISMATCH, "Vector JSON shape must be a one-element array");
			const double encodedCount = shape->array_value[0].number_value;
			if (!std::isfinite(encodedCount) || encodedCount < 0 || encodedCount != std::trunc(encodedCount))
				return SerializeFailure(IOError::SCHEMA_MISMATCH, "Vector JSON shape count must be a nonnegative integer");
			if (encodedCount > std::numeric_limits<int>::max())
				return SerializeFailure(IOError::ALLOCATION_LIMIT_EXCEEDED, "Vector JSON count exceeds supported dimension range");
			const int count = static_cast<int>(encodedCount);
			if (expectedSize >= 0 && count != expectedSize)
				return SerializeFailure(IOError::SCHEMA_MISMATCH, "Vector JSON shape does not match fixed target dimension");
			if (static_cast<std::size_t>(count) > std::numeric_limits<std::size_t>::max() / sizeof(Type))
				return SerializeFailure(IOError::ALLOCATION_LIMIT_EXCEEDED, "Vector JSON byte count overflows size_t");
			if (static_cast<std::size_t>(count) * sizeof(Type) > options.max_allocation_bytes)
				return SerializeFailure(IOError::ALLOCATION_LIMIT_EXCEEDED, "Vector JSON data exceeds configured allocation limit");
			if (data == nullptr || data->type != JsonValueType::Array || data->array_value.size() != static_cast<std::size_t>(count))
				return SerializeFailure(IOError::SCHEMA_MISMATCH, "Vector JSON data length does not match shape");

			for (const JsonValue& value : data->array_value)
				if (value.type != JsonValueType::Number)
					return SerializeFailure(IOError::MALFORMED_INPUT, "Vector JSON data values must be numbers");

			return SerializeSuccess();
		}
	}

	template<class Type>
	IOResult SaveJson(std::ostream& out, const Vector<Type>& vec, const SaveOptions& options = {}) {
		static_assert(std::is_floating_point_v<Type>, "Vector JSON persistence currently supports floating-point scalar types");
		JsonValue::Object root;
		root["mml"] = Detail::BuildVectorJsonHeader<Type>("Vector", 1);
		root["shape"] = JsonValue::ArrayValue({ JsonValue::Number(static_cast<double>(vec.size())) });

		JsonValue::Array data;
		data.reserve(static_cast<std::size_t>(vec.size()));
		for (int i = 0; i < vec.size(); ++i)
			data.push_back(JsonValue::Number(static_cast<double>(vec[i])));
		root["data"] = JsonValue::ArrayValue(std::move(data));

		return WriteJson(out, JsonValue::ObjectValue(std::move(root)), options);
	}

	template<class Type>
	IOResult LoadJson(std::istream& in, Vector<Type>& vec, const LoadOptions& options = {}) {
		static_assert(std::is_floating_point_v<Type>, "Vector JSON persistence currently supports floating-point scalar types");
		std::stringstream buffer;
		buffer << in.rdbuf();
		JsonParseResult parsed = ParseJson(buffer.str());
		if (!parsed.result.success)
			return parsed.result;

		IOResult validation = Detail::ValidateVectorJsonRoot<Type>(parsed.value, "Vector", -1, options);
		if (!validation.success)
			return validation;

		const JsonValue* shape = Detail::FindObjectField(parsed.value, "shape");
		const JsonValue* data = Detail::FindObjectField(parsed.value, "data");
		const int count = static_cast<int>(shape->array_value[0].number_value);
		vec.resize(count);
		for (int i = 0; i < count; ++i)
			vec[i] = static_cast<Type>(data->array_value[static_cast<std::size_t>(i)].number_value);
		return SerializeSuccess();
	}

	template<class Type, int N>
	IOResult SaveJson(std::ostream& out, const VectorN<Type, N>& vec, const SaveOptions& options = {}) {
		static_assert(std::is_floating_point_v<Type>, "VectorN JSON persistence currently supports floating-point scalar types");
		JsonValue::Object root;
		root["mml"] = Detail::BuildVectorJsonHeader<Type>("VectorN", 1);
		root["shape"] = JsonValue::ArrayValue({ JsonValue::Number(static_cast<double>(N)) });

		JsonValue::Array data;
		data.reserve(static_cast<std::size_t>(N));
		for (int i = 0; i < N; ++i)
			data.push_back(JsonValue::Number(static_cast<double>(vec[i])));
		root["data"] = JsonValue::ArrayValue(std::move(data));

		return WriteJson(out, JsonValue::ObjectValue(std::move(root)), options);
	}

	template<class Type, int N>
	IOResult LoadJson(std::istream& in, VectorN<Type, N>& vec, const LoadOptions& options = {}) {
		static_assert(std::is_floating_point_v<Type>, "VectorN JSON persistence currently supports floating-point scalar types");
		std::stringstream buffer;
		buffer << in.rdbuf();
		JsonParseResult parsed = ParseJson(buffer.str());
		if (!parsed.result.success)
			return parsed.result;

		IOResult validation = Detail::ValidateVectorJsonRoot<Type>(parsed.value, "VectorN", N, options);
		if (!validation.success)
			return validation;

		const JsonValue* data = Detail::FindObjectField(parsed.value, "data");
		for (int i = 0; i < N; ++i)
			vec[i] = static_cast<Type>(data->array_value[static_cast<std::size_t>(i)].number_value);
		return SerializeSuccess();
	}

	template<class Type>
	IOResult SaveJson(const Vector<Type>& vec, const std::string& filename, const SaveOptions& options = {}) {
		std::ofstream file(filename);
		if (!file.is_open())
			return SerializeFailure(IOError::FILE_NOT_OPENED, "Could not create JSON file " + filename).WithFile(filename);
		return FinishOutputFile(file, filename, SaveJson(file, vec, options));
	}

	template<class Type>
	IOResult LoadJson(const std::string& filename, Vector<Type>& vec, const LoadOptions& options = {}) {
		std::ifstream file(filename);
		if (!file.is_open())
			return SerializeFailure(IOError::FILE_NOT_OPENED, "Could not open JSON file " + filename).WithFile(filename);
		return InputResult(file, LoadJson(file, vec, options), filename);
	}

	template<class Type, int N>
	IOResult SaveJson(const VectorN<Type, N>& vec, const std::string& filename, const SaveOptions& options = {}) {
		std::ofstream file(filename);
		if (!file.is_open())
			return SerializeFailure(IOError::FILE_NOT_OPENED, "Could not create JSON file " + filename).WithFile(filename);
		return FinishOutputFile(file, filename, SaveJson(file, vec, options));
	}

	template<class Type, int N>
	IOResult LoadJson(const std::string& filename, VectorN<Type, N>& vec, const LoadOptions& options = {}) {
		std::ifstream file(filename);
		if (!file.is_open())
			return SerializeFailure(IOError::FILE_NOT_OPENED, "Could not open JSON file " + filename).WithFile(filename);
		return InputResult(file, LoadJson(file, vec, options), filename);
	}

} // namespace Persistence
} // namespace MML

#endif // MML_PERSISTENCE_VECTOR_JSON_H