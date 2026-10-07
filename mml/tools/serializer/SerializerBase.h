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

#include <fstream>
#include <iomanip>
#include <cstddef>
#include <map>
#include <string>
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
		/// @return SerializeResult with success flag and error details
		inline SerializeResult WriteRealFuncHeader(std::ostream& out, std::string type, std::string title,
												   Real x1, Real x2, int numPoints)
		{
			try
			{
				out << type << std::endl;
				out << "VERSION: " << SerializeFormatType::CURRENT_VERSION << std::endl;
				out << SanitizeHeaderLine(title) << std::endl;
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
				out << SerializeFormatType::MULTI_REAL_FUNCTION << std::endl;
				out << "VERSION: " << SerializeFormatType::CURRENT_VERSION << std::endl;
				out << SanitizeHeaderLine(title) << std::endl;
				out << numFuncs << std::endl;
				for (int i = 0; i < numFuncs; i++)
					out << SanitizeHeaderLine(legend[i]) << std::endl;
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
				out << "VERSION: " << SerializeFormatType::CURRENT_VERSION << std::endl;
				if (!title.empty())
					out << SanitizeHeaderLine(title) << std::endl;
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
				out << "VERSION: " << SerializeFormatType::CURRENT_VERSION << std::endl;
				out << SanitizeHeaderLine(title) << std::endl;
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

		/// @brief Write header for field lines data file
		inline SerializeResult WriteFieldLinesHeader(std::ostream& out, const std::string& type,
		                                              const std::string& title, int numLines)
		{
			try
			{
				out << type << std::endl;
				out << "VERSION: " << SerializeFormatType::CURRENT_VERSION << std::endl;
				out << "Title: " << SanitizeHeaderLine(title) << std::endl;
				out << "NUM_LINES: " << numLines << std::endl;
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

		/// @brief Write header for parametric surface data file
		inline SerializeResult WriteParametricSurfaceHeader(std::ostream& out, const std::string& title,
		                                                     Real u1, Real u2, int numPointsU,
		                                                     Real w1, Real w2, int numPointsW)
		{
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
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

		/// @brief Write header for a 2D scalar function data file
		inline SerializeResult WriteScalarFunc2DHeader(std::ostream& out, const std::string& title,
		                                               Real x1, Real x2, int numPointsX,
		                                               Real y1, Real y2, int numPointsY)
		{
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
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

		/// @brief Write header for a 3D scalar function data file
		inline SerializeResult WriteScalarFunc3DHeader(std::ostream& out, const std::string& title,
		                                               Real x1, Real x2, int numPointsX,
		                                               Real y1, Real y2, int numPointsY,
		                                               Real z1, Real z2, int numPointsZ)
		{
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
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

		/// @brief Write header for a 2D particle simulation data file
		inline SerializeResult WriteParticleSimulation2DHeader(std::ostream& out, int numBalls,
		                                                       Real width, Real height,
		                                                       const std::vector<std::string>& ballColors,
		                                                       const std::vector<Real>& ballRadius,
		                                                       int numSteps)
		{
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
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

		/// @brief Write header for a 3D particle simulation data file
		inline SerializeResult WriteParticleSimulation3DHeader(std::ostream& out, int numBalls,
		                                                       Real width, Real height, Real depth,
		                                                       const std::vector<std::string>& ballColors,
		                                                       const std::vector<Real>& ballRadius,
		                                                       int numSteps)
		{
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
