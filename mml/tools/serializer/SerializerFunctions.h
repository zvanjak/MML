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

#include "mml/tools/serializer/SerializerBase.h"
#include "mml/interfaces/IFunction.h"
#include "mml/base/InterpolatedFunction.h"

namespace MML
{
	namespace Serializer
	{
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
		/// @return SerializeResult with success flag and error details
		inline SerializeResult SaveRealFunc(std::ostream& out, const IRealFunction& f, std::string title,
											Real x1, Real x2, int numPoints, int precision = 15)
		{
			if (numPoints < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "numPoints must be >= 2"};
			if (x1 >= x2)
				return {false, SerializeError::INVALID_PARAMETERS, "x1 must be less than x2"};

			try
			{
				out.precision(precision);
				WriteRealFuncHeader(out, FormatType::REAL_FUNCTION, title, x1, x2, numPoints);

				Real step = (x2 - x1) / (numPoints - 1);
				for (int i = 0; i < numPoints; i++)
				{
					Real x = x1 + i * step;
					out << x << " " << f(x) << std::endl;
				}
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		/// @brief Serialize a real function to file at equally-spaced points (file version)
		/// @details Convenience wrapper that opens file and delegates to stream version
		inline SerializeResult SaveRealFunc(const IRealFunction& f, std::string title,
											Real x1, Real x2, int numPoints, std::string fileName, int precision = 15)
		{
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			auto result = SaveRealFunc(file, f, title, x1, x2, numPoints, precision);
			file.close();
			return result;
		}

		/// @brief Serialize a real function to a stream at specified points
		/// @details Core stream version - writes to any std::ostream (stringstream, cout, etc.)
		inline SerializeResult SaveRealFunc(std::ostream& out, const IRealFunction& f, std::string title,
											Vector<Real> points, int precision = 15)
		{
			if (points.size() < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "points vector must have at least 2 elements"};

			try
			{
				out.precision(precision);
				WriteRealFuncHeader(out, FormatType::REAL_FUNCTION, title, points[0], points[points.size() - 1], static_cast<int>(points.size()));

				for (int i = 0; i < points.size(); i++)
				{
					Real x = points[i];
					out << x << " " << f(x) << std::endl;
				}
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		/// @brief Serialize a real function to file at specified points (file version)
		inline SerializeResult SaveRealFunc(const IRealFunction& f, std::string title,
											Vector<Real> points, std::string fileName, int precision = 15)
		{
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			auto result = SaveRealFunc(file, f, title, points, precision);
			file.close();
			return result;
		}

		/// @brief Serialize a real function to stream saving only function values (space-optimized)
		inline SerializeResult SaveRealFuncEquallySpaced(std::ostream& out, const IRealFunction& f, std::string title,
														 Real x1, Real x2, int numPoints, int precision = 15)
		{
			if (numPoints < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "numPoints must be >= 2"};
			if (x1 >= x2)
				return {false, SerializeError::INVALID_PARAMETERS, "x1 must be less than x2"};

			try
			{
				out.precision(precision);
				WriteRealFuncHeader(out, FormatType::REAL_FUNCTION_EQUALLY_SPACED, title, x1, x2, numPoints);

				Real step = (x2 - x1) / (numPoints - 1);
				for (int i = 0; i < numPoints; i++)
				{
					Real x = x1 + i * step;
					out << f(x) << std::endl;
				}
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		/// @brief Serialize a real function saving only function values to file (space-optimized)
		inline SerializeResult SaveRealFuncEquallySpaced(const IRealFunction& f, std::string title,
														 Real x1, Real x2, int numPoints, std::string fileName, int precision = 15)
		{
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			auto result = SaveRealFuncEquallySpaced(file, f, title, x1, x2, numPoints, precision);
			file.close();
			return result;
		}

		//===================================================================================
		// Multi-function serialization
		//===================================================================================

		/// @brief Serialize multiple real functions to a stream
		inline SerializeResult SaveRealMultiFunc(std::ostream& out, const std::vector<IRealFunction*> &funcs, std::string title,
												 std::vector<std::string> legend, 
												 Real x1, Real x2, int numPoints, int precision = 15)
		{
			if (funcs.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "funcs vector cannot be empty"};
			if (numPoints < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "numPoints must be >= 2"};
			if (x1 >= x2)
				return {false, SerializeError::INVALID_PARAMETERS, "x1 must be less than x2"};
			if (funcs.size() != legend.size())
				return {false, SerializeError::INVALID_PARAMETERS, "funcs and legend sizes must match"};

			try
			{
				out.precision(precision);
				WriteRealMultiFuncHeader(out, title, funcs.size(), legend, x1, x2, numPoints);

				for (int i = 0; i < numPoints; i++)
				{
					Real x = x1 + (x2 - x1) * i / (numPoints - 1);
					out << x << " ";

					for (size_t j = 0; j < funcs.size(); j++)
						out << (*(funcs[j]))(x) << " ";

					out << std::endl;
				}
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		/// @brief Serialize multiple real functions to a single file
		inline SerializeResult SaveRealMultiFunc(const std::vector<IRealFunction*> &funcs, std::string title,
												 std::vector<std::string> legend, 
												 Real x1, Real x2, int numPoints, std::string fileName, int precision = 15)
		{
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			auto result = SaveRealMultiFunc(file, funcs, title, legend, x1, x2, numPoints, precision);
			file.close();
			return result;
		}

		/// @brief Serialize multiple linear interpolation functions to file
		inline SerializeResult SaveRealMultiFunc(const std::vector<LinearInterpRealFunc>& funcs, std::string title,
												 std::vector<std::string> legend,
												 Real x1, Real x2, int numPoints, std::string fileName)
		{
			if (funcs.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "funcs vector cannot be empty"};
			if (numPoints < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "numPoints must be >= 2"};
			if (x1 >= x2)
				return {false, SerializeError::INVALID_PARAMETERS, "x1 must be less than x2"};
			if (funcs.size() != legend.size())
				return {false, SerializeError::INVALID_PARAMETERS, "funcs and legend sizes must match"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				WriteRealMultiFuncHeader(file, title, funcs.size(), legend, x1, x2, numPoints);

				for (int i = 0; i < numPoints; i++)
				{
					Real x = x1 + (x2 - x1) * i / (numPoints - 1);
					file << x << " ";

					for (int j = 0; j < funcs.size(); j++)
						file << funcs[j](x) << " ";

					file << std::endl;
				}
				file.close();
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		/// @brief Serialize multiple polynomial interpolation functions to file
		inline SerializeResult SaveRealMultiFunc(const std::vector<PolynomInterpRealFunc>& funcs, std::string title,
												 std::vector<std::string> legend,
												 Real x1, Real x2, int numPoints, std::string fileName)
		{
			if (funcs.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "funcs vector cannot be empty"};
			if (numPoints < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "numPoints must be >= 2"};
			if (x1 >= x2)
				return {false, SerializeError::INVALID_PARAMETERS, "x1 must be less than x2"};
			if (funcs.size() != legend.size())
				return {false, SerializeError::INVALID_PARAMETERS, "funcs and legend sizes must match"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				WriteRealMultiFuncHeader(file, title, funcs.size(), legend, x1, x2, numPoints);

				for (int i = 0; i < numPoints; i++)
				{
					Real x = x1 + (x2 - x1) * i / (numPoints - 1);
					file << x << " ";

					for (int j = 0; j < funcs.size(); j++)
						file << funcs[j](x) << " ";

					file << std::endl;
				}
				file.close();
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		/// @brief Serialize multiple spline interpolation functions to file
		inline SerializeResult SaveRealMultiFunc(const std::vector<SplineInterpRealFunc>& funcs, std::string title,
												 std::vector<std::string> legend,
												 Real x1, Real x2, int numPoints, std::string fileName)
		{
			if (funcs.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "funcs vector cannot be empty"};
			if (numPoints < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "numPoints must be >= 2"};
			if (x1 >= x2)
				return {false, SerializeError::INVALID_PARAMETERS, "x1 must be less than x2"};
			if (funcs.size() != legend.size())
				return {false, SerializeError::INVALID_PARAMETERS, "funcs and legend sizes must match"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				WriteRealMultiFuncHeader(file, title, funcs.size(), legend, x1, x2, numPoints);

				for (int i = 0; i < numPoints; i++)
				{
					Real x = x1 + (x2 - x1) * i / (numPoints - 1);
					file << x << " ";

					for (int j = 0; j < funcs.size(); j++)
						file << funcs[j](x) << " ";

					file << std::endl;
				}
				file.close();
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

	} // namespace Serializer
} // namespace MML

#endif // MML_SERIALIZER_FUNCTIONS_H
