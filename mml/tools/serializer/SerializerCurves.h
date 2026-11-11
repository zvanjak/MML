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

#include <mml/tools/serializer/SerializerBase.h>
#include <mml/interfaces/IFunction.h>
#include <mml/base/Vector/VectorN.h>

namespace MML
{
	namespace Serializer
	{
		//===================================================================================
		// Parametric curve serialization
		//===================================================================================

		/// @brief Serialize a parametric curve to file at equally-spaced parameter values
		/// @return IOResult with success flag and error details
		template<int N>
		inline IOResult SaveParamCurve(const IRealToVectorFunction<N>& f, std::string inType, std::string title, Real t1, Real t2,
																	 int numPoints, std::string fileName) {
			const auto pathValidation = detail::ValidateFileName(fileName);
			if (!pathValidation.success) return pathValidation;
			if (numPoints < 2)
				return IOResult{false, IOError::INVALID_PARAMETERS, "numPoints must be >= 2"}.WithFile(fileName);
			if (t1 >= t2)
				return IOResult{false, IOError::INVALID_PARAMETERS, "t1 must be < t2"}.WithFile(fileName);

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			try
			{
				auto headerResult = WriteParamCurveHeader(file, inType, title, t1, t2, numPoints);
				if (!headerResult.success)
					return (headerResult).WithFile(fileName);

				for (int n = 0; n < numPoints; n++)
				{
					Real t = detail::SampleCoordinate(t1, t2, n, numPoints);
					file << t << " ";
					auto val = f(t);
					for (int i = 0; i < N; i++)
						file << val[i] << " ";
					file << std::endl;
				}
				return Serializer::detail::FinishFile(file, fileName, "Successfully saved parametric curve");
			}
			catch (const std::exception& e)
			{
				file.close();
				return IOResult{false, IOError::WRITE_FAILED, std::string("Write error: ") + e.what()}.WithFile(fileName);
			}
		}

		/// @brief Serialize a parametric curve to file at specified parameter values
		/// @return IOResult with success flag and error details
		template<int N>
		inline IOResult SaveParamCurve(const IRealToVectorFunction<N>& f, std::string inType, std::string title, Vector<Real> points,
																	 std::string fileName) {
			const auto pathValidation = detail::ValidateFileName(fileName);
			if (!pathValidation.success) return pathValidation;
			if (points.size() < 2)
				return IOResult{false, IOError::INVALID_PARAMETERS, "points vector must have at least 2 elements"}.WithFile(fileName);

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			try
			{
				auto headerResult = WriteParamCurveHeader(file, inType, title, points[0], points[points.size() - 1], static_cast<int>(points.size()));
				if (!headerResult.success)
					return (headerResult).WithFile(fileName);

				for (int i = 0; i < points.size(); i++)
				{
					Real t = points[i];
					file << t << " ";
					auto val = f(t);
					for (int j = 0; j < N; j++)
						file << val[j] << " ";
					file << std::endl;
				}

				return Serializer::detail::FinishFile(file, fileName, "Successfully saved parametric curve");
			}
			catch (const std::exception& e)
			{
				file.close();
				return IOResult{false, IOError::WRITE_FAILED, std::string("Write error: ") + e.what()}.WithFile(fileName);
			}
		}

		/// @brief Serialize pre-computed curve points to file at equally-spaced parameter values
		/// @return IOResult with success flag and error details
		template<int N>
		inline IOResult SaveAsParamCurve(std::vector<VectorN<Real, N>> vals, std::string inType, std::string title, Real t1, Real t2,
																		 int numPoints, std::string fileName) {
			const auto pathValidation = detail::ValidateFileName(fileName);
			if (!pathValidation.success) return pathValidation;
			if (numPoints < 2)
				return IOResult{false, IOError::INVALID_PARAMETERS, "numPoints must be >= 2"}.WithFile(fileName);
			if (t1 >= t2)
				return IOResult{false, IOError::INVALID_PARAMETERS, "t1 must be < t2"}.WithFile(fileName);
			if (vals.size() < static_cast<std::size_t>(numPoints))
				return IOResult{false, IOError::INVALID_PARAMETERS, "vals vector size must be >= numPoints"}.WithFile(fileName);

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			try
			{
				auto headerResult = WriteParamCurveHeader(file, inType, title, t1, t2, numPoints);
				if (!headerResult.success)
					return (headerResult).WithFile(fileName);

				for (int i = 0; i < numPoints; i++)
				{
					Real t = detail::SampleCoordinate(t1, t2, i, numPoints);
					file << t << " ";
					for (int j = 0; j < N; j++)
						file << vals[i][j] << " ";
					file << std::endl;
				}
				return Serializer::detail::FinishFile(file, fileName, "Successfully saved parametric curve");
			}
			catch (const std::exception& e)
			{
				file.close();
				return IOResult{false, IOError::WRITE_FAILED, std::string("Write error: ") + e.what()}.WithFile(fileName);
			}
		}

		/// @brief Serialize pre-computed curve points to file at specified parameter values
		/// @return IOResult with success flag and error details
		template<int N>
		inline IOResult SaveAsParamCurve(std::vector<VectorN<Real, N>> vals, std::string inType, std::string title, Vector<Real> points,
																		 std::string fileName) {
			const auto pathValidation = detail::ValidateFileName(fileName);
			if (!pathValidation.success) return pathValidation;
			if (points.size() < 2)
				return IOResult{false, IOError::INVALID_PARAMETERS, "points vector must have at least 2 elements"}.WithFile(fileName);
			if (vals.size() < static_cast<size_t>(points.size()))
				return IOResult{false, IOError::INVALID_PARAMETERS, "vals vector size must be >= points size"}.WithFile(fileName);

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			try
			{
				auto headerResult = WriteParamCurveHeader(file, inType, title, points[0], points[points.size() - 1], static_cast<int>(points.size()));
				if (!headerResult.success)
					return (headerResult).WithFile(fileName);

				for (int i = 0; i < points.size(); i++)
				{
					Real t = points[i];
					file << t << " ";
					for (int j = 0; j < N; j++)
						file << vals[i][j] << " ";
					file << std::endl;
				}
				return Serializer::detail::FinishFile(file, fileName, "Successfully saved parametric curve");
			}
			catch (const std::exception& e)
			{
				file.close();
				return IOResult{false, IOError::WRITE_FAILED, std::string("Write error: ") + e.what()}.WithFile(fileName);
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
		/// @return IOResult with success flag and error details
		/// @pre vec_x.size() >= 2, vec_x.size() == vec_y.size(), t1 < t2, nonempty filename without NUL
		/// @post File is closed and flushed; parameter range [t1, t2] mapped across vector indices
		/// @note Parameter values are computed linearly from t1 to t2 across indices
		inline IOResult SaveAsParamCurve2D(const Vector<Real>& vec_x, const Vector<Real>& vec_y, std::string title, std::string fileName,
																			 Real t1 = 0.0, Real t2 = 1.0) {
			// Validate parameters
			const auto pathValidation = detail::ValidateFileName(fileName);
			if (!pathValidation.success) return pathValidation;
			if (vec_x.size() < 2) {
				return IOResult{false, IOError::INVALID_PARAMETERS, "vec_x must have at least 2 elements"}.WithFile(fileName);
			}
			if (vec_x.size() != vec_y.size()) {
				return IOResult{false, IOError::INVALID_PARAMETERS, "vec_x and vec_y must have same size"}.WithFile(fileName);
			}
			if (t1 >= t2) {
				return IOResult{false, IOError::INVALID_PARAMETERS, "t1 must be < t2"}.WithFile(fileName);
			}

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open()) {
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);
			}

			try {
				if (auto headerResult = WriteParamCurveHeader(file, SerializeFormatType::PARAMETRIC_CURVE_CARTESIAN_2D, title, t1, t2, static_cast<int>(vec_x.size())); !headerResult.success)
					return (headerResult).WithFile(fileName);

				for (int i = 0; i < vec_x.size(); i++)
				{
					Real t = detail::SampleCoordinate(t1, t2, i, static_cast<int>(vec_x.size()));
					file << t << " " << vec_x[i] << " " << vec_y[i] << std::endl;
				}
				return Serializer::detail::FinishFile(file, fileName, "Successfully saved parametric curve 2D");
			}
			catch (const std::exception& e) {
				file.close();
				return IOResult{false, IOError::WRITE_FAILED, std::string("Write error: ") + e.what()}.WithFile(fileName);
			}
		}

		//===================================================================================
		// Helper/forwarding functions - IOResult versions
		//===================================================================================

		/// @brief Serialize a 2D Cartesian parametric curve to file
		/// @return IOResult with success flag and error details
		inline IOResult SaveParamCurveCartesian2DResult(const IRealToVectorFunction<2>& f, std::string title, Real t1, Real t2, int numPoints,
																										std::string fileName) {
			return (SaveParamCurve<2>(f, SerializeFormatType::PARAMETRIC_CURVE_CARTESIAN_2D, title, t1, t2, numPoints, fileName))
				.WithFile(fileName);
		}

		/// @brief Serialize a 3D Cartesian parametric curve to file
		/// @return IOResult with success flag and error details
		inline IOResult SaveParamCurveCartesian3DResult(const IRealToVectorFunction<3>& f, std::string title, Real t1, Real t2, int numPoints,
																										std::string fileName) {
			return (SaveParamCurve<3>(f, SerializeFormatType::PARAMETRIC_CURVE_CARTESIAN_3D, title, t1, t2, numPoints, fileName))
				.WithFile(fileName);
		}

	} // namespace Serializer
} // namespace MML

#endif // MML_SERIALIZER_CURVES_H
