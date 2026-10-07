///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        persistence/FunctionJSON.h                                          ///
///  Description: JSON persistence for sampled real-function data                     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_PERSISTENCE_FUNCTION_JSON_H
#define MML_PERSISTENCE_FUNCTION_JSON_H

#include <mml/base/InterpolatedFunctions/InterpolatedRealFunctionLinear.h>
#include <mml/base/InterpolatedFunctions/InterpolatedRealFunctionPolynomial.h>
#include <mml/base/InterpolatedFunctions/InterpolatedFunctionSpline.h>
#include <mml/base/Vector/Vector.h>
#include <mml/interfaces/IFunction.h>
#include <mml/tools/persistence/JSON.h>

#include <fstream>
#include <limits>
#include <sstream>
#include <string>

namespace MML
{
namespace Persistence
{
	struct SampleGrid
	{
		Real x_min = 0.0;
		Real x_max = 1.0;
		int num_points = 2;
	};

	struct SampledRealFunctionData
	{
		Vector<Real> nodes;
		Vector<Real> values;
		std::string label;
		std::string source;
	};

	enum class FunctionInterpolation
	{
		Linear,
		Polynomial,
		Spline
	};

	struct PolynomialInterpolationOptions
	{
		int order = 3;
	};

	struct SplineInterpolationOptions
	{
		Real first_derivative = SplineInterpRealFunc::NaturalBoundaryDerivative;
		Real last_derivative = SplineInterpRealFunc::NaturalBoundaryDerivative;
	};

	template<class InterpFunc>
	struct LoadFunctionResult
	{
		SerializeResult result;
		InterpFunc function;
	};

	namespace Detail
	{
		inline JsonValue BuildSampledFunctionJsonHeader()
		{
			JsonValue::Object mml;
			mml["format"] = JsonValue::String("MML_JSON");
			mml["version"] = JsonValue::Number(1);
			mml["object"] = JsonValue::String("SampledRealFunction");
			mml["object_version"] = JsonValue::Number(1);
			mml["scalar"] = JsonValue::String("Real");
			mml["scalar_bytes"] = JsonValue::Number(static_cast<double>(sizeof(Real)));
			mml["scalar_encoding"] = JsonValue::String("ieee754");
			if constexpr (std::is_same_v<Real, float>)
				mml["real_type"] = JsonValue::String("float");
			else if constexpr (std::is_same_v<Real, double>)
				mml["real_type"] = JsonValue::String("double");
			else if constexpr (std::is_same_v<Real, long double>)
				mml["real_type"] = JsonValue::String("long_double");
			else
				mml["real_type"] = JsonValue::String("custom");
			return JsonValue::ObjectValue(mml);
		}

		inline const JsonValue* FindFunctionJsonField(const JsonValue& object, const std::string& key)
		{
			if (object.type != JsonValueType::Object)
				return nullptr;
			auto it = object.object_value.find(key);
			return it == object.object_value.end() ? nullptr : &it->second;
		}

		inline SerializeResult ValidateSampledFunctionData(const SampledRealFunctionData& data)
		{
			if (data.nodes.size() < 2)
				return SerializeFailure(SerializeError::INVALID_PARAMETERS, "Sampled function requires at least two nodes");
			if (data.nodes.size() != data.values.size())
				return SerializeFailure(SerializeError::INVALID_PARAMETERS, "Sampled function nodes and values must have equal length");

			const Real firstDelta = data.nodes[1] - data.nodes[0];
			if (firstDelta == 0.0)
				return SerializeFailure(SerializeError::INVALID_PARAMETERS, "Sampled function nodes must be distinct");
			const bool increasing = firstDelta > 0.0;
			for (int i = 1; i < data.nodes.size(); ++i)
			{
				const Real delta = data.nodes[i] - data.nodes[i - 1];
				if (delta == 0.0)
					return SerializeFailure(SerializeError::INVALID_PARAMETERS, "Sampled function nodes must be distinct");
				if ((delta > 0.0) != increasing)
					return SerializeFailure(SerializeError::INVALID_PARAMETERS, "Sampled function nodes must be strictly monotonic");
			}

			return SerializeSuccess();
		}

		inline SerializeResult ValidateSampledFunctionJsonRoot(const JsonValue& root, const LoadOptions& options)
		{
			if (root.type != JsonValueType::Object)
				return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Sampled function JSON root must be an object");

			const JsonValue* mml = FindFunctionJsonField(root, "mml");
			if (mml == nullptr || mml->type != JsonValueType::Object)
				return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Sampled function JSON missing mml header object");

			const JsonValue* format = FindFunctionJsonField(*mml, "format");
			const JsonValue* version = FindFunctionJsonField(*mml, "version");
			const JsonValue* object = FindFunctionJsonField(*mml, "object");
			const JsonValue* objectVersion = FindFunctionJsonField(*mml, "object_version");
			const JsonValue* scalar = FindFunctionJsonField(*mml, "scalar");
			const JsonValue* scalarBytes = FindFunctionJsonField(*mml, "scalar_bytes");
			const JsonValue* scalarEncoding = FindFunctionJsonField(*mml, "scalar_encoding");

			if (format == nullptr || format->type != JsonValueType::String || format->string_value != "MML_JSON")
				return SerializeFailure(SerializeError::INVALID_FORMAT, "Sampled function JSON format must be MML_JSON");
			if (version == nullptr || version->type != JsonValueType::Number || static_cast<int>(version->number_value) != 1)
				return SerializeFailure(SerializeError::UNSUPPORTED_VERSION, "Unsupported sampled function JSON container version");
			if (object == nullptr || object->type != JsonValueType::String || object->string_value != "SampledRealFunction")
				return SerializeFailure(SerializeError::TYPE_MISMATCH, "Sampled function JSON object kind does not match target type");
			if (objectVersion == nullptr || objectVersion->type != JsonValueType::Number || static_cast<int>(objectVersion->number_value) != 1)
				return SerializeFailure(SerializeError::UNSUPPORTED_VERSION, "Unsupported sampled function JSON object schema version");
			if (scalar == nullptr || scalar->type != JsonValueType::String || scalar->string_value != "Real")
				return SerializeFailure(SerializeError::UNSUPPORTED_SCALAR, "Sampled function JSON scalar type must be Real");
			if (scalarBytes == nullptr || scalarBytes->type != JsonValueType::Number || static_cast<int>(scalarBytes->number_value) != static_cast<int>(sizeof(Real)))
				return SerializeFailure(SerializeError::UNSUPPORTED_SCALAR, "Sampled function JSON scalar byte size does not match Real");
			if (scalarEncoding == nullptr || scalarEncoding->type != JsonValueType::String || scalarEncoding->string_value != "ieee754")
				return SerializeFailure(SerializeError::UNSUPPORTED_SCALAR, "Sampled function JSON scalar encoding is unsupported");

			const JsonValue* domain = FindFunctionJsonField(root, "domain");
			const JsonValue* nodes = FindFunctionJsonField(root, "nodes");
			const JsonValue* values = FindFunctionJsonField(root, "values");

			if (domain == nullptr || domain->type != JsonValueType::Array || domain->array_value.size() != 2 ||
				domain->array_value[0].type != JsonValueType::Number || domain->array_value[1].type != JsonValueType::Number)
				return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Sampled function JSON domain must be a two-element numeric array");
			if (nodes == nullptr || nodes->type != JsonValueType::Array)
				return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Sampled function JSON nodes must be an array");
			if (values == nullptr || values->type != JsonValueType::Array)
				return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Sampled function JSON values must be an array");
			if (nodes->array_value.size() < 2)
				return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Sampled function JSON requires at least two nodes");
			if (nodes->array_value.size() != values->array_value.size())
				return SerializeFailure(SerializeError::SCHEMA_MISMATCH, "Sampled function JSON nodes and values lengths must match");
			if (nodes->array_value.size() > options.max_allocation_bytes / (2 * sizeof(Real)))
				return SerializeFailure(SerializeError::ALLOCATION_LIMIT_EXCEEDED, "Sampled function JSON data exceeds configured allocation limit");

			for (const JsonValue& node : nodes->array_value)
				if (node.type != JsonValueType::Number)
					return SerializeFailure(SerializeError::MALFORMED_INPUT, "Sampled function nodes must be numeric");
			for (const JsonValue& value : values->array_value)
				if (value.type != JsonValueType::Number)
					return SerializeFailure(SerializeError::MALFORMED_INPUT, "Sampled function values must be numeric");

			return SerializeSuccess();
		}
	}

	inline SerializeResult SaveSampledFunction(std::ostream& out, const SampledRealFunctionData& data, const SaveOptions& options = {})
	{
		SerializeResult validation = Detail::ValidateSampledFunctionData(data);
		if (!validation.success)
			return validation;

		JsonValue::Object root;
		root["mml"] = Detail::BuildSampledFunctionJsonHeader();
		root["domain"] = JsonValue::ArrayValue({
			JsonValue::Number(static_cast<double>(std::min(data.nodes[0], data.nodes[data.nodes.size() - 1]))),
			JsonValue::Number(static_cast<double>(std::max(data.nodes[0], data.nodes[data.nodes.size() - 1])))
		});

		JsonValue::Array nodes;
		JsonValue::Array values;
		nodes.reserve(static_cast<std::size_t>(data.nodes.size()));
		values.reserve(static_cast<std::size_t>(data.values.size()));
		for (int i = 0; i < data.nodes.size(); ++i)
		{
			nodes.push_back(JsonValue::Number(static_cast<double>(data.nodes[i])));
			values.push_back(JsonValue::Number(static_cast<double>(data.values[i])));
		}
		root["nodes"] = JsonValue::ArrayValue(std::move(nodes));
		root["values"] = JsonValue::ArrayValue(std::move(values));

		if (options.include_metadata && (!data.label.empty() || !data.source.empty()))
		{
			JsonValue::Object metadata;
			if (!data.label.empty()) metadata["label"] = JsonValue::String(data.label);
			if (!data.source.empty()) metadata["source"] = JsonValue::String(data.source);
			root["metadata"] = JsonValue::ObjectValue(std::move(metadata));
		}

		return WriteJson(out, JsonValue::ObjectValue(std::move(root)), options);
	}

	inline SerializeResult SaveSampledFunction(const SampledRealFunctionData& data, const std::string& filename, const SaveOptions& options = {})
	{
		std::ofstream file(filename);
		if (!file.is_open())
			return SerializeFailure(SerializeError::FILE_NOT_OPENED, "Could not create sampled function JSON file " + filename);
		return SaveSampledFunction(file, data, options);
	}

	inline SerializeResult SaveSampledFunction(std::ostream& out, const IRealFunction& function, std::string label, SampleGrid grid, const SaveOptions& options = {})
	{
		if (grid.num_points < 2)
			return SerializeFailure(SerializeError::INVALID_PARAMETERS, "SampleGrid num_points must be at least 2");
		if (grid.x_min >= grid.x_max)
			return SerializeFailure(SerializeError::INVALID_PARAMETERS, "SampleGrid x_min must be less than x_max");

		SampledRealFunctionData data;
		data.nodes.Resize(grid.num_points);
		data.values.Resize(grid.num_points);
		data.label = std::move(label);
		data.source = "sampled from IRealFunction";
		const Real step = (grid.x_max - grid.x_min) / (grid.num_points - 1);
		for (int i = 0; i < grid.num_points; ++i)
		{
			const Real x = grid.x_min + i * step;
			data.nodes[i] = x;
			data.values[i] = function(x);
		}
		return SaveSampledFunction(out, data, options);
	}

	inline SerializeResult SaveSampledFunction(const IRealFunction& function, std::string label, SampleGrid grid, const std::string& filename, const SaveOptions& options = {})
	{
		std::ofstream file(filename);
		if (!file.is_open())
			return SerializeFailure(SerializeError::FILE_NOT_OPENED, "Could not create sampled function JSON file " + filename);
		return SaveSampledFunction(file, function, std::move(label), grid, options);
	}

	inline SerializeResult SaveInterpolatedFunction(std::ostream& out, const RealFunctionInterpolated& function, std::string label = {}, const SaveOptions& options = {})
	{
		SampledRealFunctionData data;
		const int count = function.getNumPoints();
		data.nodes.Resize(count);
		data.values.Resize(count);
		data.label = std::move(label);
		data.source = function.InterpolationMethodName();
		for (int i = 0; i < count; ++i)
		{
			data.nodes[i] = function.X(i);
			data.values[i] = function.Y(i);
		}
		return SaveSampledFunction(out, data, options);
	}

	inline SerializeResult SaveInterpolatedFunction(const RealFunctionInterpolated& function, const std::string& filename, std::string label = {}, const SaveOptions& options = {})
	{
		std::ofstream file(filename);
		if (!file.is_open())
			return SerializeFailure(SerializeError::FILE_NOT_OPENED, "Could not create sampled function JSON file " + filename);
		return SaveInterpolatedFunction(file, function, std::move(label), options);
	}

	inline SerializeResult LoadSampledFunctionData(std::istream& in, SampledRealFunctionData& data, const LoadOptions& options = {})
	{
		std::stringstream buffer;
		buffer << in.rdbuf();
		JsonParseResult parsed = ParseJson(buffer.str());
		if (!parsed.result.success)
			return parsed.result;

		SerializeResult validation = Detail::ValidateSampledFunctionJsonRoot(parsed.value, options);
		if (!validation.success)
			return validation;

		const JsonValue* nodes = Detail::FindFunctionJsonField(parsed.value, "nodes");
		const JsonValue* values = Detail::FindFunctionJsonField(parsed.value, "values");
		const int count = static_cast<int>(nodes->array_value.size());
		SampledRealFunctionData loaded;
		loaded.nodes.Resize(count);
		loaded.values.Resize(count);
		for (int i = 0; i < count; ++i)
		{
			loaded.nodes[i] = static_cast<Real>(nodes->array_value[static_cast<std::size_t>(i)].number_value);
			loaded.values[i] = static_cast<Real>(values->array_value[static_cast<std::size_t>(i)].number_value);
		}

		const JsonValue* metadata = Detail::FindFunctionJsonField(parsed.value, "metadata");
		if (metadata != nullptr && metadata->type == JsonValueType::Object)
		{
			const JsonValue* label = Detail::FindFunctionJsonField(*metadata, "label");
			const JsonValue* source = Detail::FindFunctionJsonField(*metadata, "source");
			if (label != nullptr && label->type == JsonValueType::String) loaded.label = label->string_value;
			if (source != nullptr && source->type == JsonValueType::String) loaded.source = source->string_value;
		}

		validation = Detail::ValidateSampledFunctionData(loaded);
		if (!validation.success)
			return validation;

		data = std::move(loaded);
		return SerializeSuccess();
	}

	inline SerializeResult LoadSampledFunctionData(const std::string& filename, SampledRealFunctionData& data, const LoadOptions& options = {})
	{
		std::ifstream file(filename);
		if (!file.is_open())
			return SerializeFailure(SerializeError::FILE_NOT_OPENED, "Could not open sampled function JSON file " + filename);
		return LoadSampledFunctionData(file, data, options);
	}

	inline LoadFunctionResult<LinearInterpRealFunc> LoadLinearFunction(std::istream& in, bool extrapolateOutsideOfRange = false, const LoadOptions& options = {})
	{
		SampledRealFunctionData data;
		SerializeResult result = LoadSampledFunctionData(in, data, options);
		if (!result.success)
			return { result, LinearInterpRealFunc(Vector<Real>{0.0, 1.0}, Vector<Real>{0.0, 0.0}) };
		return { SerializeSuccess(), LinearInterpRealFunc(data.nodes, data.values, extrapolateOutsideOfRange) };
	}

	inline LoadFunctionResult<LinearInterpRealFunc> LoadLinearFunction(const std::string& filename, bool extrapolateOutsideOfRange = false, const LoadOptions& options = {})
	{
		std::ifstream file(filename);
		if (!file.is_open())
			return { SerializeFailure(SerializeError::FILE_NOT_OPENED, "Could not open sampled function JSON file " + filename), LinearInterpRealFunc(Vector<Real>{0.0, 1.0}, Vector<Real>{0.0, 0.0}) };
		return LoadLinearFunction(file, extrapolateOutsideOfRange, options);
	}

	inline LoadFunctionResult<PolynomInterpRealFunc> LoadPolynomialFunction(std::istream& in, PolynomialInterpolationOptions interpolationOptions = {}, const LoadOptions& options = {})
	{
		SampledRealFunctionData data;
		SerializeResult result = LoadSampledFunctionData(in, data, options);
		if (!result.success)
			return { result, PolynomInterpRealFunc(Vector<Real>{0.0, 1.0}, Vector<Real>{0.0, 0.0}, 2) };
		if (interpolationOptions.order < 2 || interpolationOptions.order > data.nodes.size())
			return { SerializeFailure(SerializeError::INVALID_PARAMETERS, "Polynomial interpolation order must be in [2, number of samples]"), PolynomInterpRealFunc(Vector<Real>{0.0, 1.0}, Vector<Real>{0.0, 0.0}, 2) };
		return { SerializeSuccess(), PolynomInterpRealFunc(data.nodes, data.values, interpolationOptions.order) };
	}

	inline LoadFunctionResult<PolynomInterpRealFunc> LoadPolynomialFunction(const std::string& filename, PolynomialInterpolationOptions interpolationOptions = {}, const LoadOptions& options = {})
	{
		std::ifstream file(filename);
		if (!file.is_open())
			return { SerializeFailure(SerializeError::FILE_NOT_OPENED, "Could not open sampled function JSON file " + filename), PolynomInterpRealFunc(Vector<Real>{0.0, 1.0}, Vector<Real>{0.0, 0.0}, 2) };
		return LoadPolynomialFunction(file, interpolationOptions, options);
	}

	inline LoadFunctionResult<SplineInterpRealFunc> LoadSplineFunction(std::istream& in, SplineInterpolationOptions interpolationOptions = {}, const LoadOptions& options = {})
	{
		SampledRealFunctionData data;
		SerializeResult result = LoadSampledFunctionData(in, data, options);
		if (!result.success)
			return { result, SplineInterpRealFunc(Vector<Real>{0.0, 1.0}, Vector<Real>{0.0, 0.0}) };
		return { SerializeSuccess(), SplineInterpRealFunc(data.nodes, data.values, interpolationOptions.first_derivative, interpolationOptions.last_derivative) };
	}

	inline LoadFunctionResult<SplineInterpRealFunc> LoadSplineFunction(const std::string& filename, SplineInterpolationOptions interpolationOptions = {}, const LoadOptions& options = {})
	{
		std::ifstream file(filename);
		if (!file.is_open())
			return { SerializeFailure(SerializeError::FILE_NOT_OPENED, "Could not open sampled function JSON file " + filename), SplineInterpRealFunc(Vector<Real>{0.0, 1.0}, Vector<Real>{0.0, 0.0}) };
		return LoadSplineFunction(file, interpolationOptions, options);
	}

} // namespace Persistence
} // namespace MML

#endif // MML_PERSISTENCE_FUNCTION_JSON_H