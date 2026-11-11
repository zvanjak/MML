///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        serializer/SerializerFunctions.h                                    ///
///  Description: Real function serialization utilities                               ///
///               Save single and multiple real functions to files/streams            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_SERIALIZER_FUNCTIONS_H
#define MML_SERIALIZER_FUNCTIONS_H

#include <mml/tools/serializer/SerializerBase.h>
#include <mml/interfaces/IFunction.h>
#include <mml/base/InterpolatedFunction.h>

namespace MML
{
	namespace Serializer
	{
		namespace detail
		{
			inline IOResult ValidateRealFunctionRange(Real x1, Real x2, int numPoints) {
				if (numPoints < 2)
					return SerializeFailure(IOError::INVALID_PARAMETERS, "numPoints must be >= 2");
				if (x1 >= x2)
					return SerializeFailure(IOError::INVALID_PARAMETERS, "x1 must be less than x2");
				return SerializeSuccess();
			}

			inline IOResult ValidateRealFunctionPoints(const Vector<Real>& points) {
				if (points.size() < 2)
					return SerializeFailure(IOError::INVALID_PARAMETERS, "points vector must have at least 2 elements");
				return SerializeSuccess();
			}
		}

		//===================================================================================
		// Real function serialization
		//===================================================================================

		/// @brief Serialize a real function to a stream at equally-spaced points
		/// @details Core stream version - writes to any std::ostream (stringstream, cout, etc.)
		/// @param out Output stream to write to
		/// @param f The real function to serialize
		/// @param title Display title for the data
		/// @param x1 Start of parameter range
		/// @param x2 End of parameter range
		/// @param numPoints Number of evaluation points (must be >= 2)
		/// @param precision Decimal places for output (default: 15)
		/// @return IOResult with success flag and error details
		inline IOResult SaveRealFunc(std::ostream& out, const IRealFunction& f, std::string title, Real x1, Real x2, int numPoints,
																 int precision = 15) {
			const auto validation = detail::ValidateRealFunctionRange(x1, x2, numPoints);
			if (!validation.success) return validation;

			try
			{
				Serializer::detail::ApplyTextFormat(out, precision);
				if (auto headerResult = WriteRealFuncHeader(out, SerializeFormatType::REAL_FUNCTION, title, x1, x2, numPoints); !headerResult.success)
					return headerResult;

				for (int i = 0; i < numPoints; i++)
				{
					Real x = detail::SampleCoordinate(x1, x2, i, numPoints);
					out << x << " " << f(x) << std::endl;
				}
				return Serializer::detail::StreamResult(out);
			}
			catch (const std::exception& e)
			{
				return {false, IOError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		/// @brief Serialize a real function to file at equally-spaced points (file version)
		/// @details Convenience wrapper that opens file and delegates to stream version
		inline IOResult SaveRealFunc(const IRealFunction& f, std::string title, Real x1, Real x2, int numPoints, std::string fileName,
																 int precision = 15) {
			const auto pathValidation = detail::ValidateFileName(fileName);
			if (!pathValidation.success) return pathValidation;
			const auto validation = detail::ValidateRealFunctionRange(x1, x2, numPoints);
			if (!validation.success) return validation.WithFile(fileName);

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			auto result = SaveRealFunc(file, f, title, x1, x2, numPoints, precision);
			return Serializer::detail::FinishFile(file, fileName, result);
		}

		/// @brief Serialize a real function to a stream at specified points
		/// @details Core stream version - writes to any std::ostream (stringstream, cout, etc.)
		inline IOResult SaveRealFunc(std::ostream& out, const IRealFunction& f, std::string title, Vector<Real> points, int precision = 15) {
			const auto validation = detail::ValidateRealFunctionPoints(points);
			if (!validation.success) return validation;

			try
			{
				Serializer::detail::ApplyTextFormat(out, precision);
				if (auto headerResult = WriteRealFuncHeader(out, SerializeFormatType::REAL_FUNCTION, title, points[0], points[points.size() - 1], static_cast<int>(points.size())); !headerResult.success)
					return headerResult;

				for (int i = 0; i < points.size(); i++)
				{
					Real x = points[i];
					out << x << " " << f(x) << std::endl;
				}
				return Serializer::detail::StreamResult(out);
			}
			catch (const std::exception& e)
			{
				return {false, IOError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		/// @brief Serialize a real function to file at specified points (file version)
		inline IOResult SaveRealFunc(const IRealFunction& f, std::string title, Vector<Real> points, std::string fileName, int precision = 15) {
			const auto pathValidation = detail::ValidateFileName(fileName);
			if (!pathValidation.success) return pathValidation;
			const auto validation = detail::ValidateRealFunctionPoints(points);
			if (!validation.success) return validation.WithFile(fileName);

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			auto result = SaveRealFunc(file, f, title, points, precision);
			return Serializer::detail::FinishFile(file, fileName, result);
		}

		/// @brief Serialize a real function to stream saving only function values (space-optimized)
		inline IOResult SaveRealFuncEquallySpaced(std::ostream& out, const IRealFunction& f, std::string title, Real x1, Real x2, int numPoints,
																							int precision = 15) {
			const auto validation = detail::ValidateRealFunctionRange(x1, x2, numPoints);
			if (!validation.success) return validation;

			try
			{
				Serializer::detail::ApplyTextFormat(out, precision);
				if (auto headerResult = WriteRealFuncHeader(out, SerializeFormatType::REAL_FUNCTION_EQUALLY_SPACED, title, x1, x2, numPoints); !headerResult.success)
					return headerResult;

				for (int i = 0; i < numPoints; i++)
				{
					Real x = detail::SampleCoordinate(x1, x2, i, numPoints);
					out << f(x) << std::endl;
				}
				return Serializer::detail::StreamResult(out);
			}
			catch (const std::exception& e)
			{
				return {false, IOError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		/// @brief Serialize a real function saving only function values to file (space-optimized)
		inline IOResult SaveRealFuncEquallySpaced(const IRealFunction& f, std::string title, Real x1, Real x2, int numPoints,
																							std::string fileName, int precision = 15) {
			const auto pathValidation = detail::ValidateFileName(fileName);
			if (!pathValidation.success) return pathValidation;
			const auto validation = detail::ValidateRealFunctionRange(x1, x2, numPoints);
			if (!validation.success) return validation.WithFile(fileName);

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			auto result = SaveRealFuncEquallySpaced(file, f, title, x1, x2, numPoints, precision);
			return Serializer::detail::FinishFile(file, fileName, result);
		}

		//===================================================================================
		// Multi-function serialization
		//===================================================================================

		namespace detail
		{
			template<class Type>
			IOResult ValidateMultiFunctionRequest(const std::vector<Type>& funcs, const std::vector<std::string>& legend,
												 Real x1, Real x2, int numPoints) {
				const auto validation = ValidateLabeledCollection(funcs, legend, "funcs");
				if (!validation.success) return validation;
				return ValidateRealFunctionRange(x1, x2, numPoints);
			}
		}

		/// @brief Serialize multiple real functions to a stream
		inline IOResult SaveRealMultiFunc(std::ostream& out, const std::vector<IRealFunction*>& funcs, std::string title,
																			std::vector<std::string> legend, Real x1, Real x2, int numPoints, int precision = 15) {
			const auto validation = detail::ValidateMultiFunctionRequest(funcs, legend, x1, x2, numPoints);
			if (!validation.success) return validation;

			try
			{
				Serializer::detail::ApplyTextFormat(out, precision);
				if (auto headerResult = WriteRealMultiFuncHeader(out, title, static_cast<int>(funcs.size()), legend, x1, x2, numPoints); !headerResult.success)
					return headerResult;

				for (int i = 0; i < numPoints; i++)
				{
					Real x = detail::SampleCoordinate(x1, x2, i, numPoints);
					out << x << " ";

					for (size_t j = 0; j < funcs.size(); j++)
						out << (*(funcs[j]))(x) << " ";

					out << std::endl;
				}
				return Serializer::detail::StreamResult(out);
			}
			catch (const std::exception& e)
			{
				return {false, IOError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		/// @brief Serialize multiple real functions to a single file
		inline IOResult SaveRealMultiFunc(const std::vector<IRealFunction*>& funcs, std::string title, std::vector<std::string> legend, Real x1,
																			Real x2, int numPoints, std::string fileName, int precision = 15) {
			const auto validation = detail::ValidateMultiFunctionRequest(funcs, legend, x1, x2, numPoints);
			if (!validation.success) return validation.WithFile(fileName);
			const auto pathValidation = detail::ValidateFileName(fileName);
			if (!pathValidation.success) return pathValidation;

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			auto result = SaveRealMultiFunc(file, funcs, title, legend, x1, x2, numPoints, precision);
			return Serializer::detail::FinishFile(file, fileName, result);
		}

		/// @brief Serialize multiple interpolation functions to file
		template<class InterpFunc>
		IOResult SaveRealMultiFunc(const std::vector<InterpFunc>& funcs, std::string title, std::vector<std::string> legend, Real x1, Real x2,
															 int numPoints, std::string fileName) {
			const auto validation = detail::ValidateMultiFunctionRequest(funcs, legend, x1, x2, numPoints);
			if (!validation.success) return validation.WithFile(fileName);
			const auto pathValidation = detail::ValidateFileName(fileName);
			if (!pathValidation.success) return pathValidation;

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			try
			{
				if (auto headerResult = WriteRealMultiFuncHeader(file, title, static_cast<int>(funcs.size()), legend, x1, x2, numPoints); !headerResult.success)
					return (headerResult).WithFile(fileName);

				for (int i = 0; i < numPoints; i++)
				{
					Real x = detail::SampleCoordinate(x1, x2, i, numPoints);
					file << x << " ";

					for (size_t j = 0; j < funcs.size(); j++)
						file << funcs[j](x) << " ";

					file << std::endl;
				}
				return Serializer::detail::FinishFile(file, fileName, "Success");
			}
			catch (const std::exception& e)
			{
				return IOResult{false, IOError::WRITE_FAILED, std::string("Write error: ") + e.what()}.WithFile(fileName);
			}
		}

	} // namespace Serializer
} // namespace MML

#endif // MML_SERIALIZER_FUNCTIONS_H
