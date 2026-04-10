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

#include "mml/MMLBase.h"

#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace MML
{
	/// Error codes for serialization operations
	enum class SerializeError
	{
		OK,                      ///< Operation succeeded
		FILE_NOT_OPENED,         ///< Failed to open/create file
		INVALID_PARAMETERS,      ///< Invalid input parameters
		WRITE_FAILED,            ///< Write operation failed
		READ_FAILED,             ///< Read operation failed
		INVALID_FORMAT           ///< Invalid file format (bad magic, version, etc.)
	};

	/// Result of a serialization operation
	struct SerializeResult
	{
		bool success;                    ///< true if successful, false otherwise
		SerializeError error;            ///< Error code (if any)
		std::string message;             ///< Human-readable error message

		/// Implicit bool conversion for backward compatibility
		operator bool() const { return success; }
	};

	//===================================================================================
	// Text format type constants (single source of truth for all serializers)
	//===================================================================================
	namespace FormatType
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
		// Header writers - common utilities used by other serializer modules
		//===================================================================================

		/// @brief Write header for real function data file
		/// @return SerializeResult with success flag and error details
		inline SerializeResult WriteRealFuncHeader(std::ostream& out, std::string type, std::string title,
												   Real x1, Real x2, int numPoints)
		{
			try
			{
				out << type << std::endl;
				out << "VERSION: " << FormatType::CURRENT_VERSION << std::endl;
				out << title << std::endl;
				out << "x1: " << x1 << std::endl; 
				out << "x2: " << x2 << std::endl;
				out << "NumPoints: " << numPoints << std::endl;
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

		/// @brief Write header for multi-function data file
		/// @return SerializeResult with success flag and error details
		inline SerializeResult WriteRealMultiFuncHeader(std::ostream& out, std::string title, int numFuncs,
														std::vector<std::string> legend, Real x1, Real x2, int numPoints)
		{
			try
			{
				out << FormatType::MULTI_REAL_FUNCTION << std::endl;
				out << "VERSION: " << FormatType::CURRENT_VERSION << std::endl;
				out << title << std::endl;
				out << numFuncs << std::endl;
				for (int i = 0; i < numFuncs; i++)
					out << legend[i] << std::endl;
				out << "x1: " << x1 << std::endl;
				out << "x2: " << x2 << std::endl;
				out << "NumPoints: " << numPoints << std::endl;
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

		/// @brief Write header for parametric curve data file
		/// @return SerializeResult with success flag and error details
		inline SerializeResult WriteParamCurveHeader(std::ostream& out, std::string type, std::string title,
													 Real t1, Real t2, int numPoints)
		{
			try
			{
				out << type << std::endl;
				out << "VERSION: " << FormatType::CURRENT_VERSION << std::endl;
				if (!title.empty())
					out << title << std::endl;
				out << "t1: " << t1 << std::endl;
				out << "t2: " << t2 << std::endl;
				out << "NumPoints: " << numPoints << std::endl;
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

		/// @brief Write header for vector field data file
		/// @return SerializeResult with success flag and error details
		inline SerializeResult WriteVectorFieldHeader(std::ostream& out, const std::string& type, const std::string& title)
		{
			try
			{
				out << type << std::endl;
				out << "VERSION: " << FormatType::CURRENT_VERSION << std::endl;
				out << title << std::endl;
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

	} // namespace Serializer
} // namespace MML

#endif // MML_SERIALIZER_BASE_H
