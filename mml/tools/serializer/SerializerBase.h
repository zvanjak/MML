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
		WRITE_FAILED             ///< Write operation failed
	};

	/// Result of a serialization operation
	struct SerializeResult
	{
		bool success;                    ///< true if successful, false otherwise
		SerializeError error;            ///< Error code (if any)
		std::string message;             ///< Human-readable error message
	};

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
				out << "MULTI_REAL_FUNCTION" << std::endl;
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
		inline void WriteVectorFieldHeader(std::ofstream& file, const std::string& type, const std::string& title)
		{
			file << type << std::endl;
			file << title << std::endl;
		}

	} // namespace Serializer
} // namespace MML

#endif // MML_SERIALIZER_BASE_H
