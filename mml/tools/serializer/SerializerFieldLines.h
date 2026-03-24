///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        serializer/SerializerFieldLines.h                                   ///
///  Description: Field line serialization utilities                                  ///
///               Save 2D/3D field lines to files for visualization                   ///
///                                                                                   ///
///  File Format: FIELD_LINES_2D / FIELD_LINES_3D                                     ///
///    - Header with title and line count                                             ///
///    - Each line: point count followed by coordinates                               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_SERIALIZER_FIELD_LINES_H
#define MML_SERIALIZER_FIELD_LINES_H

#include "mml/tools/serializer/SerializerBase.h"
#include "mml/algorithms/FieldLineTracer.h"

#include <vector>

namespace MML
{
	namespace Serializer
	{
		//===================================================================================
		// Field lines header writer
		//===================================================================================

		/// @brief Write header for field lines data file
		inline SerializeResult WriteFieldLinesHeader(std::ostream& out, const std::string& type, 
		                                              const std::string& title, int numLines)
		{
			try
			{
				out << type << std::endl;
				out << "Title: " << title << std::endl;
				out << "NUM_LINES: " << numLines << std::endl;
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Header write error: ") + e.what()};
			}
		}

		//===================================================================================
		// 2D field lines serialization
		//===================================================================================

		/// @brief Save 2D field lines to file
		/// @param lines Vector of field lines to save
		/// @param title Descriptive title for the visualization
		/// @param fileName Output file path
		/// @return SerializeResult with success/error information
		inline SerializeResult SaveFieldLines2D(const std::vector<FieldLine2D>& lines,
		                                         const std::string& title,
		                                         const std::string& fileName)
		{
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};
			if (lines.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "lines vector cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				// Write header
				auto headerResult = WriteFieldLinesHeader(file, "FIELD_LINES_2D", title, static_cast<int>(lines.size()));
				if (!headerResult.success)
					return headerResult;

				// Write each line
				for (size_t i = 0; i < lines.size(); ++i)
				{
					const auto& line = lines[i];
					file << "LINE " << i << std::endl;
					file << line.size() << std::endl;
					
					for (const auto& point : line.points)
					{
						file << std::setprecision(std::numeric_limits<Real>::max_digits10) << point[0] << " " << point[1] << std::endl;
					}
				}

				file.close();
				return {true, SerializeError::OK, "Successfully saved " + std::to_string(lines.size()) + " field lines"};
			}
			catch (const std::exception& e)
			{
				file.close();
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		/// @brief Save 2D field lines from raw point vectors
		/// @param lines Vector of point vectors (each inner vector is one field line)
		/// @param title Descriptive title for the visualization
		/// @param fileName Output file path
		/// @return SerializeResult with success/error information
		inline SerializeResult SaveFieldLines2D(const std::vector<std::vector<VectorN<Real, 2>>>& lines,
		                                         const std::string& title,
		                                         const std::string& fileName)
		{
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};
			if (lines.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "lines vector cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				// Write header
				auto headerResult = WriteFieldLinesHeader(file, "FIELD_LINES_2D", title, static_cast<int>(lines.size()));
				if (!headerResult.success)
					return headerResult;

				// Write each line
				for (size_t i = 0; i < lines.size(); ++i)
				{
					const auto& line = lines[i];
					file << "LINE " << i << std::endl;
					file << line.size() << std::endl;
					
					for (const auto& point : line)
					{
						file << std::setprecision(std::numeric_limits<Real>::max_digits10) << point[0] << " " << point[1] << std::endl;
					}
				}

				file.close();
				return {true, SerializeError::OK, "Successfully saved " + std::to_string(lines.size()) + " field lines"};
			}
			catch (const std::exception& e)
			{
				file.close();
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		//===================================================================================
		// 3D field lines serialization
		//===================================================================================

		/// @brief Save 3D field lines to file
		/// @param lines Vector of field lines to save
		/// @param title Descriptive title for the visualization
		/// @param fileName Output file path
		/// @return SerializeResult with success/error information
		inline SerializeResult SaveFieldLines3D(const std::vector<FieldLine3D>& lines,
		                                         const std::string& title,
		                                         const std::string& fileName)
		{
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};
			if (lines.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "lines vector cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				// Write header
				auto headerResult = WriteFieldLinesHeader(file, "FIELD_LINES_3D", title, static_cast<int>(lines.size()));
				if (!headerResult.success)
					return headerResult;

				// Write each line
				for (size_t i = 0; i < lines.size(); ++i)
				{
					const auto& line = lines[i];
					file << "LINE " << i << std::endl;
					file << line.size() << std::endl;
					
					for (const auto& point : line.points)
					{
						file << std::setprecision(std::numeric_limits<Real>::max_digits10) << point[0] << " " << point[1] << " " << point[2] << std::endl;
					}
				}

				file.close();
				return {true, SerializeError::OK, "Successfully saved " + std::to_string(lines.size()) + " field lines"};
			}
			catch (const std::exception& e)
			{
				file.close();
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		/// @brief Save 3D field lines from raw point vectors
		inline SerializeResult SaveFieldLines3D(const std::vector<std::vector<VectorN<Real, 3>>>& lines,
		                                         const std::string& title,
		                                         const std::string& fileName)
		{
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};
			if (lines.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "lines vector cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				// Write header
				auto headerResult = WriteFieldLinesHeader(file, "FIELD_LINES_3D", title, static_cast<int>(lines.size()));
				if (!headerResult.success)
					return headerResult;

				// Write each line
				for (size_t i = 0; i < lines.size(); ++i)
				{
					const auto& line = lines[i];
					file << "LINE " << i << std::endl;
					file << line.size() << std::endl;
					
					for (const auto& point : line)
					{
						file << std::setprecision(std::numeric_limits<Real>::max_digits10) << point[0] << " " << point[1] << " " << point[2] << std::endl;
					}
				}

				file.close();
				return {true, SerializeError::OK, "Successfully saved " + std::to_string(lines.size()) + " field lines"};
			}
			catch (const std::exception& e)
			{
				file.close();
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

	} // namespace Serializer
} // namespace MML

#endif // MML_SERIALIZER_FIELD_LINES_H
