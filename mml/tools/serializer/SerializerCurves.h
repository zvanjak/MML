///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        serializer/SerializerCurves.h                                       ///
///  Description: Parametric curve serialization utilities                            ///
///               Save 2D/3D parametric curves to files                               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_SERIALIZER_CURVES_H
#define MML_SERIALIZER_CURVES_H

#include "mml/tools/serializer/SerializerBase.h"
#include "mml/interfaces/IFunction.h"
#include "mml/base/Vector/VectorN.h"

namespace MML
{
	namespace Serializer
	{
		//===================================================================================
		// Parametric curve serialization
		//===================================================================================

		/// @brief Serialize a parametric curve to file at equally-spaced parameter values
		/// @return SerializeResult with success flag and error details
		template<int N>
		inline SerializeResult SaveParamCurve(const IRealToVectorFunction<N>& f, std::string inType, std::string title, 
		                                      Real t1, Real t2, int numPoints, std::string fileName)
		{
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};
			if (numPoints < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "numPoints must be >= 2"};
			if (t1 >= t2)
				return {false, SerializeError::INVALID_PARAMETERS, "t1 must be < t2"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				auto headerResult = WriteParamCurveHeader(file, inType, title, t1, t2, numPoints);
				if (!headerResult.success)
					return headerResult;

				Real delta = (t2 - t1) / (numPoints - 1);
				for (int n = 0; n < numPoints; n++)
				{
					Real t = t1 + n * delta;
					file << t << " ";
					auto val = f(t);
					for (int i = 0; i < N; i++)
						file << val[i] << " ";
					file << std::endl;
				}
				file.close();
				return {true, SerializeError::OK, "Successfully saved parametric curve"};
			}
			catch (const std::exception& e)
			{
				file.close();
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		/// @brief Serialize a parametric curve to file at specified parameter values
		/// @return SerializeResult with success flag and error details
		template<int N>
		inline SerializeResult SaveParamCurve(const IRealToVectorFunction<N>& f, std::string inType, std::string title, 
		                                      Vector<Real> points, std::string fileName)
		{
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};
			if (points.size() < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "points vector must have at least 2 elements"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				auto headerResult = WriteParamCurveHeader(file, inType, title, points[0], points[points.size() - 1], static_cast<int>(points.size()));
				if (!headerResult.success)
					return headerResult;

				for (int i = 0; i < points.size(); i++)
				{
					Real t = points[i];
					file << t << " ";
					auto val = f(t);
					for (int j = 0; j < N; j++)
						file << val[j] << " ";
					file << std::endl;
				}

				file.close();
				return {true, SerializeError::OK, "Successfully saved parametric curve"};
			}
			catch (const std::exception& e)
			{
				file.close();
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		/// @brief Serialize pre-computed curve points to file at equally-spaced parameter values
		/// @return SerializeResult with success flag and error details
		template<int N>
		inline SerializeResult SaveAsParamCurve(std::vector<VectorN<Real, N>> vals, std::string inType, std::string title, 
		                                        Real t1, Real t2, int numPoints, std::string fileName)
		{
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};
			if (numPoints < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "numPoints must be >= 2"};
			if (t1 >= t2)
				return {false, SerializeError::INVALID_PARAMETERS, "t1 must be < t2"};
			if (static_cast<int>(vals.size()) < numPoints)
				return {false, SerializeError::INVALID_PARAMETERS, "vals vector size must be >= numPoints"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				auto headerResult = WriteParamCurveHeader(file, inType, title, t1, t2, numPoints);
				if (!headerResult.success)
					return headerResult;

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
				return {true, SerializeError::OK, "Successfully saved parametric curve"};
			}
			catch (const std::exception& e)
			{
				file.close();
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		/// @brief Serialize pre-computed curve points to file at specified parameter values
		/// @return SerializeResult with success flag and error details
		template<int N>
		inline SerializeResult SaveAsParamCurve(std::vector<VectorN<Real, N>> vals, std::string inType, std::string title,
		                                        Vector<Real> points, std::string fileName)
		{
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};
			if (points.size() < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "points vector must have at least 2 elements"};
			if (vals.size() < points.size())
				return {false, SerializeError::INVALID_PARAMETERS, "vals vector size must be >= points size"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				auto headerResult = WriteParamCurveHeader(file, inType, title, points[0], points[points.size() - 1], static_cast<int>(points.size()));
				if (!headerResult.success)
					return headerResult;

				for (int i = 0; i < points.size(); i++)
				{
					Real t = points[i];
					file << t << " ";
					for (int j = 0; j < N; j++)
						file << vals[i][j] << " ";
					file << std::endl;
				}
				file.close();
				return {true, SerializeError::OK, "Successfully saved parametric curve"};
			}
			catch (const std::exception& e)
			{
				file.close();
				return {false, SerializeError::WRITE_FAILED, std::string("Write error: ") + e.what()};
			}
		}

		/// @brief Serialize a parametric 2D curve from component vectors
		/// @details Saves a 2D parametric curve (x(t), y(t)) from explicit vector pairs to file.
		/// Useful for plotting trajectories and curves defined by separate x and y vectors.
		/// Format: [t, x(t), y(t)] triplets, one per line.
		/// @param vec_x Vector of x-coordinates (must be same size as vec_y)
		/// @param vec_y Vector of y-coordinates (must be same size as vec_x)
		/// @param fileName Output file path for the parametric curve
		/// @param t1 Start of parameter range (default: 0.0, used to label parameter axis)
		/// @param t2 End of parameter range (default: 1.0, used to label parameter axis)
		/// @return SerializeResult with success flag and error details
		/// @pre vec_x.size() >= 2, vec_x.size() == vec_y.size(), t1 < t2, !fileName.empty()
		/// @post File is closed and flushed; parameter range [t1, t2] mapped across vector indices
		/// @note Parameter values are computed linearly from t1 to t2 across indices
		inline SerializeResult SaveAsParamCurve2D(const Vector<Real>& vec_x, const Vector<Real>& vec_y, 
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

		//===================================================================================
		// Helper/forwarding functions - SerializeResult versions
		//===================================================================================

		/// @brief Serialize a 2D Cartesian parametric curve to file
		/// @return SerializeResult with success flag and error details
		inline SerializeResult SaveParamCurveCartesian2DResult(const IRealToVectorFunction<2>& f, std::string title, 
		                                                       Real t1, Real t2, int numPoints, std::string fileName)
		{
			return SaveParamCurve<2>(f, "PARAMETRIC_CURVE_CARTESIAN_2D", title, t1, t2, numPoints, fileName);
		}
		
		/// @brief Serialize a 3D Cartesian parametric curve to file
		/// @return SerializeResult with success flag and error details
		inline SerializeResult SaveParamCurveCartesian3DResult(const IRealToVectorFunction<3>& f, std::string title, 
		                                                       Real t1, Real t2, int numPoints, std::string fileName)
		{
			return SaveParamCurve<3>(f, "PARAMETRIC_CURVE_CARTESIAN_3D", title, t1, t2, numPoints, fileName);
		}

	} // namespace Serializer
} // namespace MML

#endif // MML_SERIALIZER_CURVES_H
