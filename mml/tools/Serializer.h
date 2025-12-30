///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Serializer.h                                                        ///
///  Description: Data serialization utilities for MML objects                        ///
///               JSON, CSV, and binary export/import for vectors, matrices, etc.     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_SERIALIZER_H
#define MML_SERIALIZER_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/InterpolatedFunction.h"
#include "base/ODESystemSolution.h"

#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace MML
{
	class Serializer
	{
	public:
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

		/**
		 * @brief Serialize a real function to file at equally-spaced points
		 * @details Saves function values at equally-spaced points in an interval to file with
		 *          format: [x, f(x)] pairs, one per line. Useful for plotting and data analysis.
		 * @param f The real function to serialize
		 * @param title Display title for the data (written to file header)
		 * @param x1 Start of parameter range
		 * @param x2 End of parameter range
		 * @param numPoints Number of evaluation points (must be >= 2)
		 * @param fileName Output file path where data will be written
		 * @param precision Decimal places for output (default: 15)
		 * @return SerializeResult with success flag and error details
		 * @pre x1 < x2, numPoints >= 2, !fileName.empty()
		 * @post File is closed and flushed; existing file is overwritten
		 * @example
		 *   IRealFunction& f = ...;
		 *   auto result = Serializer::SaveRealFunc(f, "sin(x)", 0, 2*M_PI, 100, "output.txt", 8);
		 *   if (!result.success) std::cerr << result.message << std::endl;
		 * @see SaveRealFunc(const IRealFunction&, std::string, Vector<Real>, std::string, int)
		 */
		// Real function serialization 
		// serializing values at equally spaced points in given interval
		static SerializeResult SaveRealFunc(const IRealFunction& f, std::string title,
										 Real x1, Real x2, int numPoints, std::string fileName, int precision = 15)
		{
			// Validate parameters
			if (numPoints < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "numPoints must be >= 2"};
			if (x1 >= x2)
				return {false, SerializeError::INVALID_PARAMETERS, "x1 must be less than x2"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				file.precision(precision);
				WriteRealFuncHeader(file, "REAL_FUNCTION", title, x1, x2, numPoints);

				Real step = (x2 - x1) / (numPoints - 1);
				for (int i = 0; i < numPoints; i++)
				{
					Real x = x1 + i * step;
					file << x << " " << f(x) << std::endl;
				}
				file.close();
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		/**
		 * @brief Serialize a real function to file at specified points
		 * @details Saves function values at arbitrary (non-equally-spaced) points to file.
		 *          Format: [x, f(x)] pairs, one per line, with bounds from points vector.
		 * @param f The real function to serialize
		 * @param title Display title for the data (written to file header)
		 * @param points Vector of x-coordinates where function is evaluated (must be sorted)
		 * @param fileName Output file path where data will be written
		 * @param precision Decimal places for output (default: 15)
		 * @return SerializeResult with success flag and error details
		 * @pre points.size() >= 2, !fileName.empty()
		 * @post File is closed and flushed; existing file is overwritten
		 * @example
		 *   Vector<Real> points = {0, 1.5, 3.14, 6.28};
		 *   auto result = Serializer::SaveRealFunc(f, "sine", points, "sine_points.txt", 6);
		 * @see SaveRealFunc(const IRealFunction&, std::string, Real, Real, int, std::string, int)
		 */
		// serializing values at given list of points
		static SerializeResult SaveRealFunc(const IRealFunction& f, std::string title,
										 Vector<Real> points, std::string fileName, int precision = 15)
		{
			// Validate parameters
			if (points.size() < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "points vector must have at least 2 elements"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				file.precision(precision);
				WriteRealFuncHeader(file, "REAL_FUNCTION", title, points[0], points[points.size() - 1], static_cast<int>(points.size()));

				for (int i = 0; i < points.size(); i++)
				{
					Real x = points[i];
					file << x << " " << f(x) << std::endl;
				}
				file.close();
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		/**
		 * @brief Serialize a real function saving only function values (space-optimized)
		 * @details Like SaveRealFunc but omits explicit x-coordinates from file since they can be
		 *          reconstructed from x1, x2, and numPoints. File contains only f(x) values.
		 * @param f The real function to serialize
		 * @param title Display title for the data
		 * @param x1 Start of parameter range
		 * @param x2 End of parameter range
		 * @param numPoints Number of equally-spaced evaluation points (must be >= 2)
		 * @param fileName Output file path
		 * @param precision Decimal places for output (default: 15)
		 * @return SerializeResult with success flag and error details
		 * @pre x1 < x2, numPoints >= 2, !fileName.empty()
		 * @post File is closed and flushed; reduced file size vs SaveRealFunc
		 * @note Format: values only (one per line); x-values must be reconstructed by reader
		 * @see SaveRealFunc(const IRealFunction&, std::string, Real, Real, int, std::string, int)
		 */
		// same as SaveRealFunc, but points are not explicitly written in file
		// (as they can be calculated from x1, x2 and numPoints)
		static SerializeResult SaveRealFuncEquallySpaced(const IRealFunction& f, std::string title,
																				Real x1, Real x2, int numPoints, std::string fileName, int precision = 15)
		{
			// Validate parameters
			if (numPoints < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "numPoints must be >= 2"};
			if (x1 >= x2)
				return {false, SerializeError::INVALID_PARAMETERS, "x1 must be less than x2"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				file.precision(precision);
				WriteRealFuncHeader(file, "REAL_FUNCTION_EQUALLY_SPACED", title, x1, x2, numPoints);

				Real step = (x2 - x1) / (numPoints - 1);
				for (int i = 0; i < numPoints; i++)
				{
					Real x = x1 + i * step;
					file << f(x) << std::endl;
				}
				file.close();
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		// Helper function for writing real function headers
		static bool WriteRealFuncHeader(std::ofstream& file, std::string type, std::string title,
																Real x1, Real x2, int numPoints)
		{
			if (!file.is_open())
				return false;

			file << type << std::endl;
			file << title << std::endl;
			file << "x1: " << x1 << std::endl;
			file << "x2: " << x2 << std::endl;
			file << "NumPoints: " << numPoints << std::endl;

			return true;
		}

		// serializing multiple functions in a single files
		static bool WriteRealMultiFuncHeader(std::ofstream& file, std::string title, int numFuncs,
																				 std::vector<std::string> legend, Real x1, Real x2, int numPoints)
		{
			if (!file.is_open())
				return false;

			file << "MULTI_REAL_FUNCTION" << std::endl;

			file << title << std::endl;
			file << numFuncs << std::endl;
			for (int i = 0; i < numFuncs; i++)
				file << legend[i] << std::endl;

			file << "x1: " << x1 << std::endl;
			file << "x2: " << x2 << std::endl;
			file << "NumPoints: " << numPoints << std::endl;

			return true;
		}

		/**
		 * @brief Serialize multiple real functions to a single file
		 * @details Saves multiple function values at the same equally-spaced x-coordinates to file.
		 *          Each row contains: [x, f1(x), f2(x), ..., fn(x)]. Supports function pointers,
		 *          linear interpolation, polynomial interpolation, and spline interpolation functions.
		 * @param funcs Vector of real function pointers to serialize
		 * @param title Display title for the data set
		 * @param legend Vector of labels for each function (must match funcs.size())
		 * @param x1 Start of parameter range
		 * @param x2 End of parameter range
		 * @param numPoints Number of equally-spaced evaluation points (must be >= 2)
		 * @param fileName Output file path
		 * @param precision Decimal places for output (default: 15)
		 * @return SerializeResult with success flag and error details
		 * @pre x1 < x2, numPoints >= 2, funcs.size() == legend.size(), all funcs non-null, !fileName.empty()
		 * @post File is closed and flushed with header + data in columnar format
		 * @example
		 *   std::vector<IRealFunction*> funcs = {&f1, &f2};
		 *   std::vector<std::string> legend = {"sin(x)", "cos(x)"};
		 *   auto result = Serializer::SaveRealMultiFunc(funcs, "Trig", legend, 0, M_PI, 50, "trig.txt", 6);
		 * @see SaveRealMultiFunc(const std::vector<LinearInterpRealFunc>&, std::string, ...)
		 */
		static SerializeResult SaveRealMultiFunc(const std::vector<IRealFunction*> &funcs, std::string title,
												std::vector<std::string> legend, 
												Real x1, Real x2, int numPoints, std::string fileName, int precision = 15)
		{
			// Validate parameters
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
				file.precision(precision);
				WriteRealMultiFuncHeader(file, title, funcs.size(), legend, x1, x2, numPoints);

				for (int i = 0; i < numPoints; i++)
				{
					Real x = x1 + (x2 - x1) * i / (numPoints - 1);
					file << x << " ";

					for (int j = 0; j < funcs.size(); j++)
						file << (*(funcs[j]))(x) << " ";

					file << std::endl;
				}
				file.close();
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}		static SerializeResult SaveRealMultiFunc(const std::vector<LinearInterpRealFunc>& funcs, std::string title,
																	std::vector<std::string> legend,
																	Real x1, Real x2, int numPoints, std::string fileName)
		{
			// Validate parameters
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
		static SerializeResult SaveRealMultiFunc(const std::vector<PolynomInterpRealFunc>& funcs, std::string title,
																	std::vector<std::string> legend,
																	Real x1, Real x2, int numPoints, std::string fileName)
		{
			// Validate parameters
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
		static SerializeResult SaveRealMultiFunc(const std::vector<SplineInterpRealFunc>& funcs, std::string title,
																	std::vector<std::string> legend,
																	Real x1, Real x2, int numPoints, std::string fileName)
		{
			// Validate parameters
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

		// Helper function for writing parametric curve headers
		static bool WriteParamCurveHeader(std::ofstream& file, std::string type, std::string title,
																			Real t1, Real t2, int numPoints)
		{
			if (!file.is_open())
				return false;

			file << type << std::endl;
			if (!title.empty())
				file << title << std::endl;
			file << "t1: " << t1 << std::endl;
			file << "t2: " << t2 << std::endl;
			file << "NumPoints: " << numPoints << std::endl;

			return true;
		}

		// Parametric curve serialization
		template<int N>
		static bool SaveParamCurve(const IRealToVectorFunction<N>& f, std::string inType, std::string title, 
															 Real t1, Real t2, int numPoints, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			WriteParamCurveHeader(file, inType, title, t1, t2, numPoints);

			Real delta = (t2 - t1) / (numPoints - 1);
			for (Real t = t1; t <= t2; t += delta)
			{
				file << t << " ";
				for (int i = 0; i < N; i++)
					file << f(t)[i] << " ";
				file << std::endl;
			}
			file.close();
			return true;
		}

		template<int N>
		static bool SaveParamCurve(const IRealToVectorFunction<N>& f, std::string inType, std::string title, 
															 Vector<Real> points, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			WriteParamCurveHeader(file, inType, title, points[0], points[points.size() - 1], static_cast<int>(points.size()));

			for (int i = 0; i < points.size(); i++)
			{
				Real t = points[i];
				file << t << " ";
				for (int i = 0; i < N; i++)
					file << f(t)[i] << " ";
				file << std::endl;
			}

			file.close();
			return true;
		}

		template<int N>
		static bool SaveAsParamCurve(std::vector<VectorN<Real, N>> vals, std::string inType, std::string title, 
																 Real t1, Real t2, int numPoints, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			WriteParamCurveHeader(file, inType, title, t1, t2, numPoints);

			Real delta = (t2 - t1) / (numPoints - 1);
			for (int i = 0; i < numPoints; i++)
			{
				Real t = t1 + i * delta;
				file << t << " ";
				for (int j = 0; j < N; j++)
					file << vals[i][j] << " ";
				file << std::endl;
			}
			file.close();
			return true;
		}

		template<int N>
		static bool SaveAsParamCurve(std::vector<VectorN<Real, N>> vals, std::string inType, std::string title,
																 Vector<Real> points, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			WriteParamCurveHeader(file, inType, title, points[0], points[points.size() - 1], static_cast<int>(points.size()));

			for (int i = 0; i < points.size(); i++)
			{
				Real t = points[i];
				file << t << " ";
				for (int j = 0; j < N; j++)
					file << vals[i][j] << " ";
				file << std::endl;
			}
			file.close();
			return true;
		}

		/**
		 * @brief Serialize a parametric 2D curve from component vectors
		 * @details Saves a 2D parametric curve (x(t), y(t)) from explicit vector pairs to file.
		 *          Useful for plotting trajectories and curves defined by separate x and y vectors.
		 *          Format: [t, x(t), y(t)] triplets, one per line.
		 * @param vec_x Vector of x-coordinates (must be same size as vec_y)
		 * @param vec_y Vector of y-coordinates (must be same size as vec_x)
		 * @param fileName Output file path for the parametric curve
		 * @param t1 Start of parameter range (default: 0.0, used to label parameter axis)
		 * @param t2 End of parameter range (default: 1.0, used to label parameter axis)
		 * @return SerializeResult with success flag and error details
		 * @pre vec_x.size() >= 2, vec_x.size() == vec_y.size(), t1 < t2, !fileName.empty()
		 * @post File is closed and flushed; parameter range [t1, t2] mapped across vector indices
		 * @note Parameter values are computed linearly from t1 to t2 across indices
		 * @example
		 *   Vector<Real> xs = {0, 1, 2, 1, 0};
		 *   Vector<Real> ys = {0, 1.732, 0, -1.732, 0};
		 *   auto result = Serializer::SaveAsParamCurve2D(xs, ys, "Pentagon", "curve.txt", 0, 2*M_PI);
		 */
		// save parametric curve in 2D as a list of points
		static SerializeResult SaveAsParamCurve2D(const Vector<Real>& vec_x, const Vector<Real>& vec_y, 
													std::string title, std::string fileName, Real t1 = 0.0, Real t2 = 1.0)
		{
			// Validate parameters
			if (fileName.empty()) {
				return { false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty" };
			}
			if (vec_x.size() < 2) {
				return { false, SerializeError::INVALID_PARAMETERS, "vec_x must have at least 2 elements" };
			}
			if (vec_x.size() != vec_y.size()) {
				return { false, SerializeError::INVALID_PARAMETERS, "vec_x and vec_y must have same size" };
			}
			if (t1 >= t2) {
				return { false, SerializeError::INVALID_PARAMETERS, "t1 must be < t2" };
			}

			std::ofstream file(fileName);
			if (!file.is_open()) {
				return { false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName };
			}

			try {
				WriteParamCurveHeader(file, "PARAMETRIC_CURVE_CARTESIAN_2D", title, t1, t2, static_cast<int>(vec_x.size()));

				for (int i = 0; i < vec_x.size(); i++)
				{
					Real t = t1 + (t2 - t1) * i / (vec_x.size() - 1);
					file << t << " " << vec_x[i] << " " << vec_y[i] << std::endl;
				}
				file.close();
				return { true, SerializeError::OK, "Successfully saved parametric curve 2D" };
			}
			catch (const std::exception& e) {
				file.close();
				return { false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what() };
			}
		}

		// Helper/forwarding functions
		static bool SaveParamCurveCartesian2D(const IRealToVectorFunction<2>& f, std::string title, 
																					Real t1, Real t2, int numPoints, std::string fileName)
		{
			return SaveParamCurve<2>(f, "PARAMETRIC_CURVE_CARTESIAN_2D", title, t1, t2, numPoints, fileName);
		}
		
		static bool SaveParamCurveCartesian3D(const IRealToVectorFunction<3>& f, std::string title, 
																					Real t1, Real t2, int numPoints, std::string fileName)
		{
			return SaveParamCurve<3>(f, "PARAMETRIC_CURVE_CARTESIAN_3D", title, t1, t2, numPoints, fileName);
		}

		/**
		 * @brief Serialize a 2D scalar function (z=f(x,y)) to a grid file
		 * @details Saves a 2D scalar function evaluated on a rectangular grid to file.
		 *          Format: Header with grid parameters, then value rows (one per y-line, varying x).
		 *          Ideal for surface plots and heatmaps.
		 * @param f The 2D scalar function to serialize
		 * @param title Display title for the surface data
		 * @param x1 Start of x parameter range
		 * @param x2 End of x parameter range
		 * @param numPointsX Number of grid points along x (must be >= 2)
		 * @param y1 Start of y parameter range
		 * @param y2 End of y parameter range
		 * @param numPointsY Number of grid points along y (must be >= 2)
		 * @param fileName Output file path
		 * @return SerializeResult with success flag and error details
		 * @pre x1 < x2, y1 < y2, numPointsX >= 2, numPointsY >= 2, !fileName.empty()
		 * @post File contains grid layout with f(x,y) values; suitable for surface visualization
		 * @example
		 *   auto result = Serializer::SaveScalarFunc2DCartesian(f, "z=x*y", 0, 10, 20, 0, 10, 20, "surface.txt");
		 * @see SaveScalarFunc3DCartesian
		 */
		// Scalar function serialization
		static SerializeResult SaveScalarFunc2DCartesian(const IScalarFunction<2>& f, std::string title,
																				Real x1, Real x2, int numPointsX, 
																				Real y1, Real y2, int numPointsY, std::string fileName)
		{
			// Validate parameters
			if (numPointsX < 2 || numPointsY < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "numPointsX and numPointsY must be >= 2"};
			if (x1 >= x2 || y1 >= y2)
				return {false, SerializeError::INVALID_PARAMETERS, "x1 < x2 and y1 < y2 required"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				file << "SCALAR_FUNCTION_CARTESIAN_2D" << std::endl;
				file << title << std::endl;
				file << "x1: " << x1 << std::endl;
				file << "x2: " << x2 << std::endl;
				file << "NumPointsX: " << numPointsX << std::endl;
				file << "y1: " << y1 << std::endl;
				file << "y2: " << y2 << std::endl;
				file << "NumPointsY: " << numPointsY << std::endl;

				Real stepX = (x2 - x1) / (numPointsX - 1);
				Real stepY = (y2 - y1) / (numPointsY - 1);
				for (int i = 0; i < numPointsX; i++)
				{
					for (int j = 0; j < numPointsY; j++)
					{
						Real x = x1 + i * stepX;
						Real y = y1 + j * stepY;
						file << x << " " << y << " " << f(VectorN<Real, 2>{x, y}) << std::endl;
					}
				}
				file.close();
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		/**
		 * @brief Serialize a 3D scalar function (w=f(x,y,z)) to a grid file
		 * @details Saves a 3D scalar function evaluated on a rectangular 3D grid to file.
		 *          Format: Header with grid parameters, then value rows (one per point in 3D grid).
		 *          Useful for volumetric data and 3D scalar field visualization.
		 * @param f The 3D scalar function to serialize
		 * @param title Display title for the volumetric data
		 * @param x1 Start of x parameter range
		 * @param x2 End of x parameter range
		 * @param numPointsX Number of grid points along x (must be >= 2)
		 * @param y1 Start of y parameter range
		 * @param y2 End of y parameter range
		 * @param numPointsY Number of grid points along y (must be >= 2)
		 * @param z1 Start of z parameter range
		 * @param z2 End of z parameter range
		 * @param numPointsZ Number of grid points along z (must be >= 2)
		 * @param fileName Output file path
		 * @return SerializeResult with success flag and error details
		 * @pre All bounds ordered, all numPoints >= 2, !fileName.empty()
		 * @post File contains volumetric grid with f(x,y,z) values; large for dense grids
		 * @note Consider using coarse grids initially to manage file size
		 * @see SaveScalarFunc2DCartesian
		 */
		static SerializeResult SaveScalarFunc3DCartesian(const IScalarFunction<3>& f, std::string title, 
																				Real x1, Real x2, int numPointsX, 
																				Real y1, Real y2, int numPointsY, 
																				Real z1, Real z2, int numPointsZ, std::string fileName)
		{
			// Validate parameters
			if (numPointsX < 2 || numPointsY < 2 || numPointsZ < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "numPointsX, numPointsY, and numPointsZ must be >= 2"};
			if (x1 >= x2 || y1 >= y2 || z1 >= z2)
				return {false, SerializeError::INVALID_PARAMETERS, "x1 < x2, y1 < y2, and z1 < z2 required"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				file << "SCALAR_FUNCTION_CARTESIAN_3D" << std::endl;
				file << title << std::endl;
				file << "x1: " << x1 << std::endl;
				file << "x2: " << x2 << std::endl;
				file << "NumPointsX: " << numPointsX << std::endl;
				file << "y1: " << y1 << std::endl;
				file << "y2: " << y2 << std::endl;
				file << "NumPointsY: " << numPointsY << std::endl;
				file << "z1: " << z1 << std::endl;
				file << "z2: " << z2 << std::endl;
				file << "NumPointsZ: " << numPointsZ << std::endl;

				Real stepX = (x2 - x1) / (numPointsX - 1);
				Real stepY = (y2 - y1) / (numPointsY - 1);
				Real stepZ = (z2 - z1) / (numPointsZ - 1);
				for (int i = 0; i < numPointsX; i++)
				{
					for (int j = 0; j < numPointsY; j++)
					{
						for (int k = 0; k < numPointsZ; k++)
						{
							Real x = x1 + i * stepX;
							Real y = y1 + j * stepY;
							Real z = z1 + k * stepZ;
							file << x << " " << y << " " << z << " " << f(VectorN<Real, 3>{x, y, z}) << std::endl;
						}
					}
				}
				file.close();
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////
		// Vector field header writer - centralizes header generation for all vector field types
		///////////////////////////////////////////////////////////////////////////////////////
		static void WriteVectorFieldHeader(std::ofstream& file, const std::string& type, const std::string& title)
		{
			file << type << std::endl;
			file << title << std::endl;
		}

		// 2D vector function serialization
		static SerializeResult SaveVectorFunc2D( const IVectorFunction<2>& f, std::string inType, std::string title,
																	Real x1_start, Real x1_end, int numPointsX1,
																	Real x2_start, Real x2_end, int numPointsX2, std::string fileName)
		{
			// Validate parameters
			if (numPointsX1 < 1 || numPointsX2 < 1)
				return {false, SerializeError::INVALID_PARAMETERS, "numPointsX1 and numPointsX2 must be >= 1"};
			if (x1_start >= x1_end || x2_start >= x2_end)
				return {false, SerializeError::INVALID_PARAMETERS, "x1_start < x1_end and x2_start < x2_end required"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				WriteVectorFieldHeader(file, inType, title);

				Real stepX = (x1_end - x1_start) / (numPointsX1 - 1);
				Real stepY = (x2_end - x2_start) / (numPointsX2 - 1);
				for (int i = 0; i < numPointsX1; i++)
					for (int j = 0; j < numPointsX2; j++)
					{
						Real x = x1_start + i * stepX;
						Real y = x2_start + j * stepY;
						auto val = f(VectorN<Real, 2>{x, y});
						file << x << " " << y << " " << val[0] << " " << val[1] << std::endl;
					}
				file.close();
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		static SerializeResult SaveVectorFunc2D( const IVectorFunction<2>& f, std::string inType, std::string title,
																	Real x1_start, Real x1_end, int numPointsX1,
																	Real x2_start, Real x2_end, int numPointsX2,
																	std::string fileName, Real upper_threshold)
		{
			// Validate parameters
			if (numPointsX1 < 1 || numPointsX2 < 1)
				return {false, SerializeError::INVALID_PARAMETERS, "numPointsX1 and numPointsX2 must be >= 1"};
			if (x1_start >= x1_end || x2_start >= x2_end)
				return {false, SerializeError::INVALID_PARAMETERS, "x1_start < x1_end and x2_start < x2_end required"};
			if (upper_threshold <= 0)
				return {false, SerializeError::INVALID_PARAMETERS, "upper_threshold must be positive"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				WriteVectorFieldHeader(file, inType, title);

				Real stepX = (x1_end - x1_start) / (numPointsX1 - 1);
				Real stepY = (x2_end - x2_start) / (numPointsX2 - 1);
				for (int i = 0; i < numPointsX1; i++)
					for (int j = 0; j < numPointsX2; j++)
					{
						Real x = x1_start + i * stepX;
						Real y = x2_start + j * stepY;
						auto val = f(VectorN<Real, 2>{x, y});
						if (val.NormL2() < upper_threshold)
							file << x << " " << y << " " << val[0] << " " << val[1] << std::endl;
					}
				file.close();
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		static SerializeResult SaveVectorFunc2DCartesian(const IVectorFunction<2>& f, std::string title,
																				Real x1, Real x2, int numPointsX,
																				Real y1, Real y2, int numPointsY, std::string fileName)
		{
			return SaveVectorFunc2D(f, "VECTOR_FIELD_2D_CARTESIAN", title, x1, x2, numPointsX, y1, y2, numPointsY, fileName);
		}

		static SerializeResult SaveVectorFunc2DCartesian(const IVectorFunction<2>& f, std::string title,
																				Real x1, Real x2, int numPointsX,
																				Real y1, Real y2, int numPointsY, std::string fileName, 
																				Real upper_threshold)
		{
			return SaveVectorFunc2D(f, "VECTOR_FIELD_2D_CARTESIAN", title, x1, x2, numPointsX, y1, y2, numPointsY, fileName, upper_threshold);
		}

		// 3D vector function serialization
		static SerializeResult SaveVectorFunc3D(const IVectorFunction<3>& f, std::string inType, std::string title,
																 Real x1_start, Real x1_end, int numPointsX1, 
																 Real x2_start, Real x2_end, int numPointsX2, 
																 Real x3_start, Real x3_end, int numPointsX3, std::string fileName)
		{
			// Validate parameters
			if (numPointsX1 < 1 || numPointsX2 < 1 || numPointsX3 < 1)
				return {false, SerializeError::INVALID_PARAMETERS, "all numPoints must be >= 1"};
			if (x1_start >= x1_end || x2_start >= x2_end || x3_start >= x3_end)
				return {false, SerializeError::INVALID_PARAMETERS, "all start bounds must be less than end bounds"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				WriteVectorFieldHeader(file, inType, title);

				Real stepX = (x1_end - x1_start) / (numPointsX1 - 1);
				Real stepY = (x2_end - x2_start) / (numPointsX2 - 1);
				Real stepZ = (x3_end - x3_start) / (numPointsX3 - 1);
				for (int i = 0; i < numPointsX1; i++)
					for (int j = 0; j < numPointsX2; j++)
						for (int k = 0; k < numPointsX3; k++)
						{
							Real x = x1_start + i * stepX;
							Real y = x2_start + j * stepY;
							Real z = x3_start + k * stepZ;
							auto val = f(VectorN<Real, 3>{x, y, z});
							file << x << " " << y << " " << z << " " << val[0] << " " << val[1] << " " << val[2] << std::endl;
						}

				file.close();
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		static SerializeResult SaveVectorFunc3D(const IVectorFunction<3>& f, std::string inType, std::string title, 
																 Real x1_start, Real x1_end, int numPointsX1, 
																 Real x2_start, Real x2_end, int numPointsX2, 
																 Real x3_start, Real x3_end, int numPointsX3, 
																 std::string fileName, Real upper_threshold)
		{
			// Validate parameters
			if (numPointsX1 < 1 || numPointsX2 < 1 || numPointsX3 < 1)
				return {false, SerializeError::INVALID_PARAMETERS, "all numPoints must be >= 1"};
			if (x1_start >= x1_end || x2_start >= x2_end || x3_start >= x3_end)
				return {false, SerializeError::INVALID_PARAMETERS, "all start bounds must be less than end bounds"};
			if (upper_threshold <= 0)
				return {false, SerializeError::INVALID_PARAMETERS, "upper_threshold must be positive"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				WriteVectorFieldHeader(file, inType, title);

				Real stepX = (x1_end - x1_start) / (numPointsX1 - 1);
				Real stepY = (x2_end - x2_start) / (numPointsX2 - 1);
				Real stepZ = (x3_end - x3_start) / (numPointsX3 - 1);
				for (int i = 0; i < numPointsX1; i++)
					for (int j = 0; j < numPointsX2; j++)
						for (int k = 0; k < numPointsX3; k++)
						{
							Real x = x1_start + i * stepX;
							Real y = x2_start + j * stepY;
							Real z = x3_start + k * stepZ;
							auto val = f(VectorN<Real, 3>{x, y, z});

							if (val.NormL2() < upper_threshold)
								file << x << " " << y << " " << z << " " << val[0] << " " << val[1] << " " << val[2] << std::endl;
						}

				file.close();
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		static SerializeResult SaveVectorFunc3DCartesian(const IVectorFunction<3>& f, std::string title, 
																					Real x1, Real x2, int numPointsX, 
																					Real y1, Real y2, int numPointsY, 
																					Real z1, Real z2, int numPointsZ, std::string fileName)
		{
			return SaveVectorFunc3D(f, "VECTOR_FIELD_3D_CARTESIAN", title, x1, x2, numPointsX, y1, y2, numPointsY, z1, z2, numPointsZ, fileName);
		}
		static SerializeResult SaveVectorFunc3DCartesian(const IVectorFunction<3>& f, std::string title, 
																					Real x1, Real x2, int numPointsX, 
																					Real y1, Real y2, int numPointsY, 
																					Real z1, Real z2, int numPointsZ, std::string fileName, 
																					Real upper_threshold)
		{
			return SaveVectorFunc3D(f, "VECTOR_FIELD_3D_CARTESIAN", title, x1, x2, numPointsX, y1, y2, numPointsY, z1, z2, numPointsZ, fileName, upper_threshold);
		}
		
		static SerializeResult SaveVectorFuncSpherical(const IVectorFunction<3>& f, std::string title, Real r1, Real r2, int numPointsR, Real theta1, Real theta2, int numPointsTheta, Real phi1, Real phi2, int numPointsPhi, std::string fileName)
		{
			return SaveVectorFunc3D(f, "VECTOR_FIELD_SPHERICAL", title, r1, r2, numPointsR, theta1, theta2, numPointsTheta, phi1, phi2, numPointsPhi, fileName);
		}
		static SerializeResult SaveVectorFuncSpherical(const IVectorFunction<3>& f, std::string title, Real r1, Real r2, int numPointsR, Real theta1, Real theta2, int numPointsTheta, Real phi1, Real phi2, int numPointsPhi, std::string fileName, Real upper_threshold)
		{
			return SaveVectorFunc3D(f, "VECTOR_FIELD_SPHERICAL", title, r1, r2, numPointsR, theta1, theta2, numPointsTheta, phi1, phi2, numPointsPhi, fileName, upper_threshold);
		}

		// ODESolution serialization
		static SerializeResult SaveODESolutionComponentAsFunc(const ODESystemSolution& sol, int compInd, 
																							 std::string title, std::string fileName)
		{
			// Validate parameters
			if (sol.getSysDim() <= 0)
				return {false, SerializeError::INVALID_PARAMETERS, "ODE solution dimension must be positive"};
			if (compInd < 0 || compInd >= sol.getSysDim())
				return {false, SerializeError::INVALID_PARAMETERS, "component index out of range [0, " + std::to_string(sol.getSysDim() - 1) + "]"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				WriteRealFuncHeader(file, "REAL_FUNCTION", title, sol.getT1(), sol.getT2(), sol.getTotalSavedSteps());
				
				for (int i = 0; i < sol.getTotalSavedSteps(); i++)
				{
					file << sol.getTValues()[i] << " " << sol.getXValues()[compInd][i] << std::endl;
				}
				
				file.close();
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		static SerializeResult SaveODESolutionAsMultiFunc(const ODESystemSolution& sol, std::string title, std::vector<std::string> legend, std::string fileName)
		{
			// Validate parameters
			if (legend.size() != sol.getSysDim())
				return {false, SerializeError::INVALID_PARAMETERS, "legend size (" + std::to_string(legend.size()) + ") must match ODE dimension (" + std::to_string(sol.getSysDim()) + ")"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				WriteRealMultiFuncHeader(file, title, sol.getSysDim(), legend, sol.getT1(), sol.getT2(), sol.getTotalSavedSteps());

				for (int i = 0; i < sol.getTotalSavedSteps(); i++)
				{
					file << sol.getTValues()[i] << " ";
					for (int j = 0; j < sol.getSysDim(); j++)
					{
						file << sol.getXValues()[j][i] << " ";
					}
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

		static SerializeResult SaveODESolAsParametricCurve2D(const ODESystemSolution& sol, std::string fileName, 
																							int ind1, int ind2, std::string title)
		{
			// Validate parameters
			if (sol.getSysDim() <= 0)
				return {false, SerializeError::INVALID_PARAMETERS, "ODE solution dimension must be positive"};
			if (ind1 < 0 || ind1 >= sol.getSysDim() || ind2 < 0 || ind2 >= sol.getSysDim())
				return {false, SerializeError::INVALID_PARAMETERS, "component indices out of range [0, " + std::to_string(sol.getSysDim() - 1) + "]"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				WriteParamCurveHeader(file, "PARAMETRIC_CURVE_CARTESIAN_2D", title, sol.getT1(), sol.getT2(), sol.getTotalSavedSteps());
				
				for (int i = 0; i < sol.getTotalSavedSteps(); i++)
				{
					file << sol.getTValues()[i] << " ";
					file << sol.getXValues()[ind1][i] << " ";
					file << sol.getXValues()[ind2][i] << " ";
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
		
		static SerializeResult SaveODESolAsParametricCurve3D(const ODESystemSolution& sol, std::string fileName, 
																							int ind1, int ind2, int ind3, std::string title)
		{
			// Validate parameters
			if (sol.getSysDim() <= 0)
				return {false, SerializeError::INVALID_PARAMETERS, "ODE solution dimension must be positive"};
			if (ind1 < 0 || ind1 >= sol.getSysDim() || ind2 < 0 || ind2 >= sol.getSysDim() || ind3 < 0 || ind3 >= sol.getSysDim())
				return {false, SerializeError::INVALID_PARAMETERS, "component indices out of range [0, " + std::to_string(sol.getSysDim() - 1) + "]"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				WriteParamCurveHeader(file, "PARAMETRIC_CURVE_CARTESIAN_3D", title, sol.getT1(), sol.getT2(), sol.getTotalSavedSteps());

				for (int i = 0; i < sol.getTotalSavedSteps(); i++)
				{
					file << sol.getTValues()[i] << " ";
					file << sol.getXValues()[ind1][i] << " ";
					file << sol.getXValues()[ind2][i] << " ";
					file << sol.getXValues()[ind3][i] << " ";
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

		// particle simulation serialization
		static SerializeResult SaveParticleSimulation2D(std::string fileName, int numBalls, Real width, Real height,
																			 std::vector<std::vector<Pnt2Cart>> ballPositions, 
																			 std::vector<std::string> ballColors, std::vector<Real> ballRadius,
																			 Real dT, int saveEveryNSteps = 1)
		{
			// Validate parameters
			if (numBalls <= 0)
				return {false, SerializeError::INVALID_PARAMETERS, "numBalls must be positive"};
			if (width <= 0 || height <= 0)
				return {false, SerializeError::INVALID_PARAMETERS, "width and height must be positive"};
			if (dT <= 0)
				return {false, SerializeError::INVALID_PARAMETERS, "dT must be positive"};
			if (ballPositions.size() != numBalls || ballColors.size() != numBalls || ballRadius.size() != numBalls)
				return {false, SerializeError::INVALID_PARAMETERS, "ballPositions, ballColors, and ballRadius sizes must equal numBalls"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				std::ostringstream buffer;
				buffer << "PARTICLE_SIMULATION_DATA_2D\n";
				buffer << "Width: " << width << "\n";
				buffer << "Height: " << height << "\n";
				buffer << "NumBalls: " << numBalls << "\n";

				for (int i=0; i<numBalls; i++)
				{
					buffer << "Ball_" << i + 1 << " " << ballColors[i] << " " << ballRadius[i] << std::endl;
				}

				int numSteps = ballPositions[0].size() ;
				buffer << "NumSteps: " << numSteps / saveEveryNSteps << std::endl;

				int realStep = 0;
				for (int i = 0; i < numSteps; i+=saveEveryNSteps, realStep++)
				{
					buffer << "Step " << realStep << " " << i * dT << std::endl;
					for (int j = 0; j < numBalls; j++)
					{
						buffer << j << " " << ballPositions[j][i].X() << " " << ballPositions[j][i].Y() << "\n";
					}
				}
				file << buffer.str();
				file.close();
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		static SerializeResult SaveParticleSimulation3D(std::string fileName, int numBalls, Real width, Real height, Real depth,
																				 std::vector<std::vector<Pnt3Cart>> ballPositions,
																				 std::vector<std::string> ballColors, std::vector<Real> ballRadius,
																				 Real dT, int saveEveryNSteps = 1)
		{
			// Validate parameters
			if (numBalls <= 0)
				return {false, SerializeError::INVALID_PARAMETERS, "numBalls must be positive"};
			if (width <= 0 || height <= 0 || depth <= 0)
				return {false, SerializeError::INVALID_PARAMETERS, "width, height, and depth must be positive"};
			if (dT <= 0)
				return {false, SerializeError::INVALID_PARAMETERS, "dT must be positive"};
			if (ballPositions.size() != numBalls || ballColors.size() != numBalls || ballRadius.size() != numBalls)
				return {false, SerializeError::INVALID_PARAMETERS, "ballPositions, ballColors, and ballRadius sizes must equal numBalls"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				file << "PARTICLE_SIMULATION_DATA_3D" << std::endl;
				file << "Width: "    << width << std::endl;
				file << "Height: "   << height << std::endl;
				file << "Depth: "    << depth << std::endl;
				file << "NumBalls: " << numBalls << std::endl;

				for (int i = 0; i < numBalls; i++)
				{
					file << "Ball_" << i+1 << " " << ballColors[i] << " " << ballRadius[i] << std::endl;
				}

				int numSteps = ballPositions[0].size();
				file << "NumSteps: " << numSteps << std::endl;

				for (int i = 0; i < numSteps; i++)
				{
					file << "Step " << i << " " << i * dT << std::endl;
					for (int j = 0; j < numBalls; j++)
					{
						file << j << " " << ballPositions[j][i].X() << " " << ballPositions[j][i].Y() << " " << ballPositions[j][i].Z() << "\n";
					}
				}
				file.close();
				return {true, SerializeError::OK, "Success"};
			}
			catch (const std::exception& e)
			{
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}
	};
}
#endif 
