///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        serializer/SerializerBase.h                                         ///
///  Description: Base types and utilities for serialization                          ///
///               Error codes, result struct, and common header writers               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_SERIALIZER_BASE_H
#define MML_SERIALIZER_BASE_H

#include <mml/MMLBase.h>
#include <mml/tools/IOResult.h>
#include <mml/tools/TextIOFormat.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <ios>
#include <limits>
#include <locale>
#include <cstddef>
#include <map>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace MML
{
	//===================================================================================
	// Text format type constants (single source of truth for all serializers)
	//===================================================================================
	namespace SerializeFormatType
	{
		constexpr const char* REAL_FUNCTION                  = "MML_REAL_FUNCTION";
		constexpr const char* REAL_FUNCTION_EQUALLY_SPACED   = "MML_REAL_FUNCTION_EQUALLY_SPACED";
		constexpr const char* MULTI_REAL_FUNCTION            = "MML_MULTI_REAL_FUNCTION";
		constexpr const char* PARAMETRIC_CURVE_CARTESIAN_2D  = "MML_PARAMETRIC_CURVE_CARTESIAN_2D";
		constexpr const char* PARAMETRIC_CURVE_CARTESIAN_3D  = "MML_PARAMETRIC_CURVE_CARTESIAN_3D";
		constexpr const char* PARAMETRIC_SURFACE_CARTESIAN   = "MML_PARAMETRIC_SURFACE_CARTESIAN";
		constexpr const char* SCALAR_FUNCTION_CARTESIAN_2D   = "MML_SCALAR_FUNCTION_CARTESIAN_2D";
		constexpr const char* SCALAR_FUNCTION_CARTESIAN_3D   = "MML_SCALAR_FUNCTION_CARTESIAN_3D";
		constexpr const char* VECTOR_FIELD_2D_CARTESIAN      = "MML_VECTOR_FIELD_2D_CARTESIAN";
		constexpr const char* VECTOR_FIELD_3D_CARTESIAN      = "MML_VECTOR_FIELD_3D_CARTESIAN";
		constexpr const char* VECTOR_FIELD_SPHERICAL         = "MML_VECTOR_FIELD_SPHERICAL";
		constexpr const char* FIELD_LINES_2D                 = "MML_FIELD_LINES_2D";
		constexpr const char* FIELD_LINES_3D                 = "MML_FIELD_LINES_3D";
		constexpr const char* PARTICLE_SIMULATION_DATA_2D    = "MML_PARTICLE_SIMULATION_DATA_2D";
		constexpr const char* PARTICLE_SIMULATION_DATA_3D    = "MML_PARTICLE_SIMULATION_DATA_3D";

		constexpr int CURRENT_VERSION = 1;
	}

	/// @brief Serializer namespace - data serialization utilities for MML objects
	/// @details Provides functions for saving mathematical objects to files and streams.
	/// All functions are free functions within the Serializer namespace.
	namespace Serializer
	{
		//===================================================================================
		// Stream hygiene helpers used by every serializer
		//===================================================================================
		namespace detail
		{
			// Interpolate the fraction first: a finite interval need not have a representable span or step.
			inline Real SampleCoordinate(Real first, Real last, int index, int count) {
				if (index == 0) return first;
				if (index == count - 1) return last;
				return std::lerp(first, last, Real(index) / Real(count - 1));
			}

			inline IOResult ValidateFileName(const std::string& fileName) {
				if (fileName.empty())
					return IOResult{false, IOError::INVALID_PARAMETERS, "fileName cannot be empty"}.WithFile(fileName);
				if (fileName.find('\0') != std::string::npos)
					return IOResult{false, IOError::INVALID_PARAMETERS, "fileName cannot contain embedded NUL"}.WithFile(fileName);
				return SerializeSuccess();
			}

			inline IOResult ValidateCollectionCount(std::size_t count, const char* name) {
				if (count == 0 || count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
					return SerializeFailure(IOError::INVALID_PARAMETERS, std::string(name) + " count must be in [1, INT_MAX]");
				return SerializeSuccess();
			}

			template<class Type>
			IOResult ValidateCollection(const std::vector<Type>& values, const char* name) {
				const auto result = ValidateCollectionCount(values.size(), name);
				if (!result.success) return result;
				if constexpr (std::is_pointer_v<Type>)
					for (std::size_t i = 0; i < values.size(); ++i)
						if (values[i] == nullptr)
							return SerializeFailure(IOError::INVALID_PARAMETERS,
								std::string(name) + " contains a null pointer at index " + std::to_string(i));
				return SerializeSuccess();
			}

			template<class Type>
			IOResult ValidateLabeledCollection(const std::vector<Type>& values, const std::vector<std::string>& legend, const char* name) {
				const auto result = ValidateCollection(values, name);
				if (!result.success) return result;
				if (values.size() != legend.size())
					return SerializeFailure(IOError::INVALID_PARAMETERS, std::string(name) + " and legend sizes must match");
				return SerializeSuccess();
			}

			inline IOResult ValidateHeaderCollectionCount(int count, std::size_t size, const char* name) {
				if (count < 0 || static_cast<std::size_t>(count) != size)
					return SerializeFailure(IOError::INVALID_PARAMETERS, std::string(name) + " size must match the nonnegative header count");
				return SerializeSuccess();
			}

			inline IOResult ValidateParticleCatalog(int numBalls, const std::vector<std::string>& colors, const std::vector<Real>& radii) {
				auto validation = ValidateHeaderCollectionCount(numBalls, colors.size(), "ballColors");
				if (!validation.success) return validation;
				validation = ValidateHeaderCollectionCount(numBalls, radii.size(), "ballRadius");
				if (!validation.success) return validation;
				for (std::size_t i = 0; i < colors.size(); ++i) {
					const auto& color = colors[i];
					if (color.empty() || std::any_of(color.begin(), color.end(),
						[](unsigned char character) { return character <= ' ' || character == '\x7f'; }))
						return SerializeFailure(IOError::INVALID_PARAMETERS, "ballColors[" + std::to_string(i) +
							"] must be a nonempty token without whitespace or control characters");
				}
				return SerializeSuccess();
			}

			/// @brief Put a text stream into the one format MML files are defined in:
			/// classic "C" locale (a '.' decimal point and no thousands grouping regardless of the
			/// process locale) and round-trip precision for Real by default.
			/// Without this, a program running under e.g. a German locale writes "3,14" and every
			/// reader in MML and the visualizers rejects or misparses the file (analysis v4, TOOLS T1).
			inline void ApplyTextFormat(std::ios& stream, int precision = std::numeric_limits<Real>::max_digits10)
			{
				TextIO::ApplyClassicFormat(stream, precision);
			}

			/// @brief Success result unless the stream has failed. Standard streams do not throw by
			/// default, so a try/catch around `<<` never fires; the state bits are the real signal.
			inline IOResult StreamResult(const std::ostream& out, std::string successMessage = "Success") {
				if (!out)
					return {false, IOError::WRITE_FAILED, "Stream write failed (failbit/badbit set)"};
				return {true, IOError::OK, std::move(successMessage)};
			}

			/// @brief Flush and close an output file, reporting any failure that the writes left
			/// behind or that close() itself produced (disk full, removed media, closed handle).
			inline IOResult FinishFile(std::ofstream& file, const std::string& fileName, std::string successMessage = "Success") {
				return FinishOutputFile(file, fileName, SerializeSuccess(std::move(successMessage)));
			}

			/// @brief FinishFile for the "open file, delegate to stream version" wrappers:
			/// a failed delegate result is returned as-is, otherwise the file is flushed/closed and checked.
			inline IOResult FinishFile(std::ofstream& file, const std::string& fileName, const IOResult& delegateResult) {
				return FinishOutputFile(file, fileName, delegateResult);
			}
		} // namespace detail

		//===================================================================================
		// Header writers - common utilities used by other serializer modules
		//===================================================================================

		/// @brief Sanitize a user-provided string for a line-based text header.
		/// Newlines/carriage returns would silently corrupt the format, so they become spaces.
		inline std::string SanitizeHeaderLine(std::string text)
		{
			for (char& character : text)
				if (character == '\n' || character == '\r')
					character = ' ';
			return text;
		}

		/// @brief Write header for real function data file
		/// @return IOResult with success flag and error details
		inline IOResult WriteRealFuncHeader(std::ostream& out, std::string type, std::string title, Real x1, Real x2, int numPoints) {
			try
			{
				out << type << std::endl;
				out << "VERSION: " << SerializeFormatType::CURRENT_VERSION << std::endl;
				out << SanitizeHeaderLine(title) << std::endl;
				out << "x1: " << x1 << std::endl; 
				out << "x2: " << x2 << std::endl;
				out << "NumPoints: " << numPoints << std::endl;
				return detail::StreamResult(out);
			}
			catch (const std::exception& e)
			{
				return {false, IOError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

		/// @brief Write header for multi-function data file
		/// @return IOResult with success flag and error details
		inline IOResult WriteRealMultiFuncHeader(std::ostream& out, std::string title, int numFuncs, std::vector<std::string> legend, Real x1,
																						 Real x2, int numPoints) {
			const auto validation = detail::ValidateHeaderCollectionCount(numFuncs, legend.size(), "legend");
			if (!validation.success) return validation;
			try
			{
				out << SerializeFormatType::MULTI_REAL_FUNCTION << std::endl;
				out << "VERSION: " << SerializeFormatType::CURRENT_VERSION << std::endl;
				out << SanitizeHeaderLine(title) << std::endl;
				out << numFuncs << std::endl;
				for (int i = 0; i < numFuncs; i++)
					out << SanitizeHeaderLine(legend[i]) << std::endl;
				out << "x1: " << x1 << std::endl;
				out << "x2: " << x2 << std::endl;
				out << "NumPoints: " << numPoints << std::endl;
				return detail::StreamResult(out);
			}
			catch (const std::exception& e)
			{
				return {false, IOError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

		/// @brief Write header for parametric curve data file
		/// @return IOResult with success flag and error details
		inline IOResult WriteParamCurveHeader(std::ostream& out, std::string type, std::string title, Real t1, Real t2, int numPoints) {
			try
			{
				out << type << std::endl;
				out << "VERSION: " << SerializeFormatType::CURRENT_VERSION << std::endl;
				out << SanitizeHeaderLine(title) << std::endl;
				out << "t1: " << t1 << std::endl;
				out << "t2: " << t2 << std::endl;
				out << "NumPoints: " << numPoints << std::endl;
				return detail::StreamResult(out);
			}
			catch (const std::exception& e)
			{
				return {false, IOError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

		/// @brief Write header for vector field data file
		/// @return IOResult with success flag and error details
		inline IOResult WriteVectorFieldHeader(std::ostream& out, const std::string& type, const std::string& title) {
			try
			{
				out << type << std::endl;
				out << "VERSION: " << SerializeFormatType::CURRENT_VERSION << std::endl;
				out << SanitizeHeaderLine(title) << std::endl;
				return detail::StreamResult(out);
			}
			catch (const std::exception& e)
			{
				return {false, IOError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

		/// @brief Write header for field lines data file
		inline IOResult WriteFieldLinesHeader(std::ostream& out, const std::string& type, const std::string& title, int numLines) {
			try
			{
				out << type << std::endl;
				out << "VERSION: " << SerializeFormatType::CURRENT_VERSION << std::endl;
				out << "Title: " << SanitizeHeaderLine(title) << std::endl;
				out << "NUM_LINES: " << numLines << std::endl;
				return detail::StreamResult(out);
			}
			catch (const std::exception& e)
			{
				return {false, IOError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

		/// @brief Write header for parametric surface data file
		inline IOResult WriteParametricSurfaceHeader(std::ostream& out, const std::string& title, Real u1, Real u2, int numPointsU, Real w1,
																								 Real w2, int numPointsW) {
			try
			{
				out << SerializeFormatType::PARAMETRIC_SURFACE_CARTESIAN << std::endl;
				out << "VERSION: " << SerializeFormatType::CURRENT_VERSION << std::endl;
				out << SanitizeHeaderLine(title) << std::endl;
				out << "u1: " << u1 << std::endl;
				out << "u2: " << u2 << std::endl;
				out << "NumPointsU: " << numPointsU << std::endl;
				out << "w1: " << w1 << std::endl;
				out << "w2: " << w2 << std::endl;
				out << "NumPointsW: " << numPointsW << std::endl;
				return detail::StreamResult(out);
			}
			catch (const std::exception& e)
			{
				return {false, IOError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

		/// @brief Write header for a 2D scalar function data file
		inline IOResult WriteScalarFunc2DHeader(std::ostream& out, const std::string& title, Real x1, Real x2, int numPointsX, Real y1, Real y2,
																						int numPointsY) {
			try
			{
				out << SerializeFormatType::SCALAR_FUNCTION_CARTESIAN_2D << std::endl;
				out << "VERSION: " << SerializeFormatType::CURRENT_VERSION << std::endl;
				out << SanitizeHeaderLine(title) << std::endl;
				out << "x1: " << x1 << std::endl;
				out << "x2: " << x2 << std::endl;
				out << "NumPointsX: " << numPointsX << std::endl;
				out << "y1: " << y1 << std::endl;
				out << "y2: " << y2 << std::endl;
				out << "NumPointsY: " << numPointsY << std::endl;
				return detail::StreamResult(out);
			}
			catch (const std::exception& e)
			{
				return {false, IOError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

		/// @brief Write header for a 3D scalar function data file
		inline IOResult WriteScalarFunc3DHeader(std::ostream& out, const std::string& title, Real x1, Real x2, int numPointsX, Real y1, Real y2,
																						int numPointsY, Real z1, Real z2, int numPointsZ) {
			try
			{
				out << SerializeFormatType::SCALAR_FUNCTION_CARTESIAN_3D << std::endl;
				out << "VERSION: " << SerializeFormatType::CURRENT_VERSION << std::endl;
				out << SanitizeHeaderLine(title) << std::endl;
				out << "x1: " << x1 << std::endl;
				out << "x2: " << x2 << std::endl;
				out << "NumPointsX: " << numPointsX << std::endl;
				out << "y1: " << y1 << std::endl;
				out << "y2: " << y2 << std::endl;
				out << "NumPointsY: " << numPointsY << std::endl;
				out << "z1: " << z1 << std::endl;
				out << "z2: " << z2 << std::endl;
				out << "NumPointsZ: " << numPointsZ << std::endl;
				return detail::StreamResult(out);
			}
			catch (const std::exception& e)
			{
				return {false, IOError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

		/// @brief Write header for a 2D particle simulation data file
		inline IOResult WriteParticleSimulation2DHeader(std::ostream& out, int numBalls, Real width, Real height,
																										const std::vector<std::string>& ballColors, const std::vector<Real>& ballRadius,
																										int numSteps) {
			const auto validation = detail::ValidateParticleCatalog(numBalls, ballColors, ballRadius);
			if (!validation.success) return validation;
			try
			{
				out << SerializeFormatType::PARTICLE_SIMULATION_DATA_2D << std::endl;
				out << "VERSION: " << SerializeFormatType::CURRENT_VERSION << std::endl;
				out << "Width: " << width << std::endl;
				out << "Height: " << height << std::endl;
				out << "NumBalls: " << numBalls << std::endl;
				for (int i = 0; i < numBalls; i++)
					out << "Ball_" << i + 1 << " " << ballColors[i] << " " << ballRadius[i] << std::endl;
				out << "NumSteps: " << numSteps << std::endl;
				return detail::StreamResult(out);
			}
			catch (const std::exception& e)
			{
				return {false, IOError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

		/// @brief Write header for a 3D particle simulation data file
		inline IOResult WriteParticleSimulation3DHeader(std::ostream& out, int numBalls, Real width, Real height, Real depth,
																										const std::vector<std::string>& ballColors, const std::vector<Real>& ballRadius,
																										int numSteps) {
			const auto validation = detail::ValidateParticleCatalog(numBalls, ballColors, ballRadius);
			if (!validation.success) return validation;
			try
			{
				out << SerializeFormatType::PARTICLE_SIMULATION_DATA_3D << std::endl;
				out << "VERSION: " << SerializeFormatType::CURRENT_VERSION << std::endl;
				out << "Width: " << width << std::endl;
				out << "Height: " << height << std::endl;
				out << "Depth: " << depth << std::endl;
				out << "NumBalls: " << numBalls << std::endl;
				for (int i = 0; i < numBalls; i++)
					out << "Ball_" << i + 1 << " " << ballColors[i] << " " << ballRadius[i] << std::endl;
				out << "NumSteps: " << numSteps << std::endl;
				return detail::StreamResult(out);
			}
			catch (const std::exception& e)
			{
				return {false, IOError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

	} // namespace Serializer
} // namespace MML

#endif // MML_SERIALIZER_BASE_H
