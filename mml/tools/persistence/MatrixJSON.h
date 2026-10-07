///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        persistence/MatrixJSON.h                                            ///
///  Description: JSON persistence for Matrix and MatrixNM                            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_PERSISTENCE_MATRIX_JSON_H
#define MML_PERSISTENCE_MATRIX_JSON_H

#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/tools/persistence/JSON.h>

#include <fstream>
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
		constexpr const char* MatrixJsonScalarName()
		{
			if constexpr (std::is_same_v<Type, Real>) return "Real";
			else if constexpr (std::is_same_v<Type, float>) return "float";
			else if constexpr (std::is_same_v<Type, double>) return "double";
			else if constexpr (std::is_same_v<Type, long double>) return "long_double";
			else return "unsupported";
		}

		template<class Type>
		constexpr const char* MatrixJsonConcreteScalarName()
		{
			if constexpr (std::is_same_v<Type, float>) return "float";
			else if constexpr (std::is_same_v<Type, double>) return "double";
			else if constexpr (std::is_same_v<Type, long double>) return "long_double";
			else return "unsupported";
		}

		template<class Type>
		JsonValue BuildMatrixJsonHeader(const char* objectName, int objectVersion)
		{
			JsonValue::Object mml;
			mml["format"] = JsonValue::String("MML_JSON");
			mml["version"] = JsonValue::Number(1);
			mml["object"] = JsonValue::String(objectName);
			mml["object_version"] = JsonValue::Number(objectVersion);
			mml["scalar"] = JsonValue::String(MatrixJsonScalarName<Type>());
			mml["scalar_bytes"] = JsonValue::Number(static_cast<double>(sizeof(Type)));
			mml["scalar_encoding"] = JsonValue::String("ieee754");
			mml["layout"] = JsonValue::String("row_major");
			if constexpr (std::is_same_v<Type, Real>)
				mml["real_type"] = JsonValue::String(MatrixJsonConcreteScalarName<Real>());
			return JsonValue::ObjectValue(mml);
		}

		inline const JsonValue* FindMatrixJsonField(const JsonValue& object, const std::string& key)
		{
			if (object.type != JsonValueType::Object)
				return nullptr;
			auto it = object.object_value.find(key);
			return it == object.object_value.end() ? nullptr : &it->second;
		}

		template<class Type>
		SerializeResult ValidateMatrixJsonRoot(const JsonValue& root, const char* objectName, int expectedRows, int expectedCols, const LoadOptions& options)
		{
			if (root.type != JsonValueType::Object)
				return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Matrix JSON root must be an object");

			const JsonValue* mml = FindMatrixJsonField(root, "mml");
			if (mml == nullptr || mml->type != JsonValueType::Object)
				return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Matrix JSON missing mml header object");

			const JsonValue* format = FindMatrixJsonField(*mml, "format");
			const JsonValue* version = FindMatrixJsonField(*mml, "version");
			const JsonValue* object = FindMatrixJsonField(*mml, "object");
			const JsonValue* objectVersion = FindMatrixJsonField(*mml, "object_version");
			const JsonValue* scalar = FindMatrixJsonField(*mml, "scalar");
			const JsonValue* scalarBytes = FindMatrixJsonField(*mml, "scalar_bytes");
			const JsonValue* scalarEncoding = FindMatrixJsonField(*mml, "scalar_encoding");
			const JsonValue* layout = FindMatrixJsonField(*mml, "layout");

			if (format == nullptr || format->type != JsonValueType::String || format->string_value != "MML_JSON")
				return SerializeFailure(SerializeError::INVALID_FORMAT, "Matrix JSON format must be MML_JSON");
			if (version == nullptr || version->type != JsonValueType::Number || static_cast<int>(version->number_value) != 1)
				return SerializeFailure(SerializeError::UNSUPPORTED_VERSION, "Unsupported Matrix JSON container version");
			if (object == nullptr || object->type != JsonValueType::String || object->string_value != objectName)
				return SerializeFailure(SerializeError::TYPE_MISMATCH, "Matrix JSON object kind does not match target type");
			if (objectVersion == nullptr || objectVersion->type != JsonValueType::Number || static_cast<int>(objectVersion->number_value) != 1)
				return SerializeFailure(SerializeError::UNSUPPORTED_VERSION, "Unsupported Matrix JSON object schema version");
			if (scalar == nullptr || scalar->type != JsonValueType::String || scalar->string_value != MatrixJsonScalarName<Type>())
				return SerializeFailure(SerializeError::UNSUPPORTED_SCALAR, "Matrix JSON scalar type does not match target type");
			if (scalarBytes == nullptr || scalarBytes->type != JsonValueType::Number || static_cast<int>(scalarBytes->number_value) != static_cast<int>(sizeof(Type)))
				return SerializeFailure(SerializeError::UNSUPPORTED_SCALAR, "Matrix JSON scalar byte size does not match target type");
			if (scalarEncoding == nullptr || scalarEncoding->type != JsonValueType::String || scalarEncoding->string_value != "ieee754")
				return SerializeFailure(SerializeError::UNSUPPORTED_SCALAR, "Matrix JSON scalar encoding is unsupported");
			if (layout == nullptr || layout->type != JsonValueType::String || layout->string_value != "row_major")
				return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Matrix JSON layout must be row_major");

			const JsonValue* shape = FindMatrixJsonField(root, "shape");
			const JsonValue* data = FindMatrixJsonField(root, "data");
			if (shape == nullptr || shape->type != JsonValueType::Array || shape->array_value.size() != 2 ||
				shape->array_value[0].type != JsonValueType::Number || shape->array_value[1].type != JsonValueType::Number)
				return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Matrix JSON shape must be a two-element array");

			const int rows = static_cast<int>(shape->array_value[0].number_value);
			const int cols = static_cast<int>(shape->array_value[1].number_value);
			if (rows < 0 || cols < 0)
				return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Matrix JSON shape dimensions cannot be negative");
			if (expectedRows >= 0 && rows != expectedRows)
				return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Matrix JSON row count does not match fixed target dimension");
			if (expectedCols >= 0 && cols != expectedCols)
				return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Matrix JSON column count does not match fixed target dimension");

			const std::size_t count = static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols);
			if (cols != 0 && count / static_cast<std::size_t>(cols) != static_cast<std::size_t>(rows))
				return SerializeFailure(SerializeError::ALLOCATION_LIMIT_EXCEEDED, "Matrix JSON shape overflows element count");
			if (count > std::numeric_limits<std::size_t>::max() / sizeof(Type))
				return SerializeFailure(SerializeError::ALLOCATION_LIMIT_EXCEEDED, "Matrix JSON byte count overflows size_t");
			if (count * sizeof(Type) > options.max_allocation_bytes)
				return SerializeFailure(SerializeError::ALLOCATION_LIMIT_EXCEEDED, "Matrix JSON data exceeds configured allocation limit");
			if (data == nullptr || data->type != JsonValueType::Array || data->array_value.size() != count)
				return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Matrix JSON data length does not match shape");

			for (const JsonValue& value : data->array_value)
				if (value.type != JsonValueType::Number)
					return SerializeFailure(SerializeError::MALFORMED_INPUT, "Matrix JSON data values must be numbers");

			return SerializeSuccess();
		}
	}

	template<class Type>
	SerializeResult SaveJson(std::ostream& out, const Matrix<Type>& mat, const SaveOptions& options = {})
	{
		static_assert(std::is_floating_point_v<Type>, "Matrix JSON persistence currently supports floating-point scalar types");
		JsonValue::Object root;
		root["mml"] = Detail::BuildMatrixJsonHeader<Type>("Matrix", 1);
		root["shape"] = JsonValue::ArrayValue({ JsonValue::Number(static_cast<double>(mat.rows())), JsonValue::Number(static_cast<double>(mat.cols())) });

		JsonValue::Array data;
		data.reserve(static_cast<std::size_t>(mat.rows()) * static_cast<std::size_t>(mat.cols()));
		for (int i = 0; i < mat.rows(); ++i)
			for (int j = 0; j < mat.cols(); ++j)
				data.push_back(JsonValue::Number(static_cast<double>(mat(i, j))));
		root["data"] = JsonValue::ArrayValue(std::move(data));

		return WriteJson(out, JsonValue::ObjectValue(std::move(root)), options);
	}

	template<class Type>
	SerializeResult LoadJson(std::istream& in, Matrix<Type>& mat, const LoadOptions& options = {})
	{
		static_assert(std::is_floating_point_v<Type>, "Matrix JSON persistence currently supports floating-point scalar types");
		std::stringstream buffer;
		buffer << in.rdbuf();
		JsonParseResult parsed = ParseJson(buffer.str());
		if (!parsed.result.success)
			return parsed.result;

		SerializeResult validation = Detail::ValidateMatrixJsonRoot<Type>(parsed.value, "Matrix", -1, -1, options);
		if (!validation.success)
			return validation;

		const JsonValue* shape = Detail::FindMatrixJsonField(parsed.value, "shape");
		const JsonValue* data = Detail::FindMatrixJsonField(parsed.value, "data");
		const int rows = static_cast<int>(shape->array_value[0].number_value);
		const int cols = static_cast<int>(shape->array_value[1].number_value);
		mat.Resize(rows, cols);
		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				mat(i, j) = static_cast<Type>(data->array_value[static_cast<std::size_t>(i) * cols + j].number_value);
		return SerializeSuccess();
	}

	template<class Type, int N, int M>
	SerializeResult SaveJson(std::ostream& out, const MatrixNM<Type, N, M>& mat, const SaveOptions& options = {})
	{
		static_assert(std::is_floating_point_v<Type>, "MatrixNM JSON persistence currently supports floating-point scalar types");
		JsonValue::Object root;
		root["mml"] = Detail::BuildMatrixJsonHeader<Type>("MatrixNM", 1);
		root["shape"] = JsonValue::ArrayValue({ JsonValue::Number(static_cast<double>(N)), JsonValue::Number(static_cast<double>(M)) });

		JsonValue::Array data;
		data.reserve(static_cast<std::size_t>(N) * static_cast<std::size_t>(M));
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < M; ++j)
				data.push_back(JsonValue::Number(static_cast<double>(mat(i, j))));
		root["data"] = JsonValue::ArrayValue(std::move(data));

		return WriteJson(out, JsonValue::ObjectValue(std::move(root)), options);
	}

	template<class Type, int N, int M>
	SerializeResult LoadJson(std::istream& in, MatrixNM<Type, N, M>& mat, const LoadOptions& options = {})
	{
		static_assert(std::is_floating_point_v<Type>, "MatrixNM JSON persistence currently supports floating-point scalar types");
		std::stringstream buffer;
		buffer << in.rdbuf();
		JsonParseResult parsed = ParseJson(buffer.str());
		if (!parsed.result.success)
			return parsed.result;

		SerializeResult validation = Detail::ValidateMatrixJsonRoot<Type>(parsed.value, "MatrixNM", N, M, options);
		if (!validation.success)
			return validation;

		const JsonValue* data = Detail::FindMatrixJsonField(parsed.value, "data");
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < M; ++j)
				mat(i, j) = static_cast<Type>(data->array_value[static_cast<std::size_t>(i) * M + j].number_value);
		return SerializeSuccess();
	}

	template<class Type>
	SerializeResult SaveJson(const Matrix<Type>& mat, const std::string& filename, const SaveOptions& options = {})
	{
		std::ofstream file(filename);
		if (!file.is_open())
			return SerializeFailure(SerializeError::FILE_NOT_OPENED, "Could not create JSON file " + filename);
		return SaveJson(file, mat, options);
	}

	template<class Type>
	SerializeResult LoadJson(const std::string& filename, Matrix<Type>& mat, const LoadOptions& options = {})
	{
		std::ifstream file(filename);
		if (!file.is_open())
			return SerializeFailure(SerializeError::FILE_NOT_OPENED, "Could not open JSON file " + filename);
		return LoadJson(file, mat, options);
	}

	template<class Type, int N, int M>
	SerializeResult SaveJson(const MatrixNM<Type, N, M>& mat, const std::string& filename, const SaveOptions& options = {})
	{
		std::ofstream file(filename);
		if (!file.is_open())
			return SerializeFailure(SerializeError::FILE_NOT_OPENED, "Could not create JSON file " + filename);
		return SaveJson(file, mat, options);
	}

	template<class Type, int N, int M>
	SerializeResult LoadJson(const std::string& filename, MatrixNM<Type, N, M>& mat, const LoadOptions& options = {})
	{
		std::ifstream file(filename);
		if (!file.is_open())
			return SerializeFailure(SerializeError::FILE_NOT_OPENED, "Could not open JSON file " + filename);
		return LoadJson(file, mat, options);
	}

} // namespace Persistence
} // namespace MML

#endif // MML_PERSISTENCE_MATRIX_JSON_H