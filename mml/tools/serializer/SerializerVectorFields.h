///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        serializer/SerializerVectorFields.h                                 ///
///  Description: Vector field serialization utilities                               ///
///               Save 2D/3D vector fields to files                                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_SERIALIZER_VECTOR_FIELDS_H
#define MML_SERIALIZER_VECTOR_FIELDS_H

#include <mml/tools/serializer/SerializerBase.h>
#include <mml/interfaces/IFunction.h>
#include <mml/base/Vector/VectorN.h>

namespace MML
{
	namespace Serializer
	{
		//===================================================================================
		// 2D vector function serialization
		//===================================================================================

		inline IOResult SaveVectorFunc2D(const IVectorFunction<2>& f, std::string inType, std::string title, Real x1_start, Real x1_end,
																		 int numPointsX1, Real x2_start, Real x2_end, int numPointsX2, std::string fileName) {
			// Validate parameters
			if (numPointsX1 < 2 || numPointsX2 < 2)
				return IOResult{false, IOError::INVALID_PARAMETERS, "numPointsX1 and numPointsX2 must be >= 2"}.WithFile(fileName);
			if (x1_start >= x1_end || x2_start >= x2_end)
				return IOResult{false, IOError::INVALID_PARAMETERS, "x1_start < x1_end and x2_start < x2_end required"}.WithFile(fileName);
			const auto pathValidation = detail::ValidateFileName(fileName);
			if (!pathValidation.success) return pathValidation;

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			try
			{
				if (auto headerResult = WriteVectorFieldHeader(file, inType, title); !headerResult.success)
					return (headerResult).WithFile(fileName);

				for (int i = 0; i < numPointsX1; i++)
					for (int j = 0; j < numPointsX2; j++)
					{
						Real x = detail::SampleCoordinate(x1_start, x1_end, i, numPointsX1);
						Real y = detail::SampleCoordinate(x2_start, x2_end, j, numPointsX2);
						auto val = f(VectorN<Real, 2>{x, y});
						file << x << " " << y << " " << val[0] << " " << val[1] << std::endl;
					}
				return Serializer::detail::FinishFile(file, fileName, "Success");
			}
			catch (const std::exception& e)
			{
				return IOResult{false, IOError::WRITE_FAILED, std::string("Write error: ") + e.what()}.WithFile(fileName);
			}
		}

		inline IOResult SaveVectorFunc2D(const IVectorFunction<2>& f, std::string inType, std::string title, Real x1_start, Real x1_end,
																		 int numPointsX1, Real x2_start, Real x2_end, int numPointsX2, std::string fileName,
																		 Real upper_threshold) {
			// Validate parameters
			if (numPointsX1 < 2 || numPointsX2 < 2)
				return IOResult{false, IOError::INVALID_PARAMETERS, "numPointsX1 and numPointsX2 must be >= 2"}.WithFile(fileName);
			if (x1_start >= x1_end || x2_start >= x2_end)
				return IOResult{false, IOError::INVALID_PARAMETERS, "x1_start < x1_end and x2_start < x2_end required"}.WithFile(fileName);
			if (!(upper_threshold > 0))
				return IOResult{false, IOError::INVALID_PARAMETERS, "upper_threshold must be positive"}.WithFile(fileName);
			const auto pathValidation = detail::ValidateFileName(fileName);
			if (!pathValidation.success) return pathValidation;

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			try
			{
				if (auto headerResult = WriteVectorFieldHeader(file, inType, title); !headerResult.success)
					return (headerResult).WithFile(fileName);

				for (int i = 0; i < numPointsX1; i++)
					for (int j = 0; j < numPointsX2; j++)
					{
						Real x = detail::SampleCoordinate(x1_start, x1_end, i, numPointsX1);
						Real y = detail::SampleCoordinate(x2_start, x2_end, j, numPointsX2);
						auto val = f(VectorN<Real, 2>{x, y});
						if (val.NormL2() < upper_threshold)
							file << x << " " << y << " " << val[0] << " " << val[1] << std::endl;
					}
				return Serializer::detail::FinishFile(file, fileName, "Success");
			}
			catch (const std::exception& e)
			{
				return IOResult{false, IOError::WRITE_FAILED, std::string("Write error: ") + e.what()}.WithFile(fileName);
			}
		}

		inline IOResult SaveVectorFunc2DCartesian(const IVectorFunction<2>& f, std::string title, Real x1, Real x2, int numPointsX, Real y1,
																							Real y2, int numPointsY, std::string fileName) {
			return (SaveVectorFunc2D(f, SerializeFormatType::VECTOR_FIELD_2D_CARTESIAN, title, x1, x2, numPointsX, y1, y2, numPointsY, fileName))
				.WithFile(fileName);
		}

		inline IOResult SaveVectorFunc2DCartesian(const IVectorFunction<2>& f, std::string title, Real x1, Real x2, int numPointsX, Real y1,
																							Real y2, int numPointsY, std::string fileName, Real upper_threshold) {
			return (SaveVectorFunc2D(f, SerializeFormatType::VECTOR_FIELD_2D_CARTESIAN, title, x1, x2, numPointsX, y1, y2, numPointsY, fileName,
															 upper_threshold))
				.WithFile(fileName);
		}

		//===================================================================================
		// 3D vector function serialization
		//===================================================================================

		inline IOResult SaveVectorFunc3D(const IVectorFunction<3>& f, std::string inType, std::string title, Real x1_start, Real x1_end,
																		 int numPointsX1, Real x2_start, Real x2_end, int numPointsX2, Real x3_start, Real x3_end,
																		 int numPointsX3, std::string fileName) {
			// Validate parameters
			if (numPointsX1 < 2 || numPointsX2 < 2 || numPointsX3 < 2)
				return IOResult{false, IOError::INVALID_PARAMETERS, "all numPoints must be >= 2"}.WithFile(fileName);
			if (x1_start >= x1_end || x2_start >= x2_end || x3_start >= x3_end)
				return IOResult{false, IOError::INVALID_PARAMETERS, "all start bounds must be less than end bounds"}.WithFile(fileName);
			const auto pathValidation = detail::ValidateFileName(fileName);
			if (!pathValidation.success) return pathValidation;

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			try
			{
				if (auto headerResult = WriteVectorFieldHeader(file, inType, title); !headerResult.success)
					return (headerResult).WithFile(fileName);

				for (int i = 0; i < numPointsX1; i++)
					for (int j = 0; j < numPointsX2; j++)
						for (int k = 0; k < numPointsX3; k++)
						{
							Real x = detail::SampleCoordinate(x1_start, x1_end, i, numPointsX1);
							Real y = detail::SampleCoordinate(x2_start, x2_end, j, numPointsX2);
							Real z = detail::SampleCoordinate(x3_start, x3_end, k, numPointsX3);
							auto val = f(VectorN<Real, 3>{x, y, z});
							file << x << " " << y << " " << z << " " << val[0] << " " << val[1] << " " << val[2] << std::endl;
						}

				return Serializer::detail::FinishFile(file, fileName, "Success");
			}
			catch (const std::exception& e)
			{
				return IOResult{false, IOError::WRITE_FAILED, std::string("Write error: ") + e.what()}.WithFile(fileName);
			}
		}

		inline IOResult SaveVectorFunc3D(const IVectorFunction<3>& f, std::string inType, std::string title, Real x1_start, Real x1_end,
																		 int numPointsX1, Real x2_start, Real x2_end, int numPointsX2, Real x3_start, Real x3_end,
																		 int numPointsX3, std::string fileName, Real upper_threshold) {
			// Validate parameters
			if (numPointsX1 < 2 || numPointsX2 < 2 || numPointsX3 < 2)
				return IOResult{false, IOError::INVALID_PARAMETERS, "all numPoints must be >= 2"}.WithFile(fileName);
			if (x1_start >= x1_end || x2_start >= x2_end || x3_start >= x3_end)
				return IOResult{false, IOError::INVALID_PARAMETERS, "all start bounds must be less than end bounds"}.WithFile(fileName);
			if (!(upper_threshold > 0))
				return IOResult{false, IOError::INVALID_PARAMETERS, "upper_threshold must be positive"}.WithFile(fileName);
			const auto pathValidation = detail::ValidateFileName(fileName);
			if (!pathValidation.success) return pathValidation;

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			try
			{
				if (auto headerResult = WriteVectorFieldHeader(file, inType, title); !headerResult.success)
					return (headerResult).WithFile(fileName);

				for (int i = 0; i < numPointsX1; i++)
					for (int j = 0; j < numPointsX2; j++)
						for (int k = 0; k < numPointsX3; k++)
						{
							Real x = detail::SampleCoordinate(x1_start, x1_end, i, numPointsX1);
							Real y = detail::SampleCoordinate(x2_start, x2_end, j, numPointsX2);
							Real z = detail::SampleCoordinate(x3_start, x3_end, k, numPointsX3);
							auto val = f(VectorN<Real, 3>{x, y, z});

							if (val.NormL2() < upper_threshold)
								file << x << " " << y << " " << z << " " << val[0] << " " << val[1] << " " << val[2] << std::endl;
						}

				return Serializer::detail::FinishFile(file, fileName, "Success");
			}
			catch (const std::exception& e)
			{
				return IOResult{false, IOError::WRITE_FAILED, std::string("Write error: ") + e.what()}.WithFile(fileName);
			}
		}

		inline IOResult SaveVectorFunc3DCartesian(const IVectorFunction<3>& f, std::string title, Real x1, Real x2, int numPointsX, Real y1,
																							Real y2, int numPointsY, Real z1, Real z2, int numPointsZ, std::string fileName) {
			return (SaveVectorFunc3D(f, SerializeFormatType::VECTOR_FIELD_3D_CARTESIAN, title, x1, x2, numPointsX, y1, y2, numPointsY, z1, z2,
															 numPointsZ, fileName))
				.WithFile(fileName);
		}

		inline IOResult SaveVectorFunc3DCartesian(const IVectorFunction<3>& f, std::string title, Real x1, Real x2, int numPointsX, Real y1,
																							Real y2, int numPointsY, Real z1, Real z2, int numPointsZ, std::string fileName,
																							Real upper_threshold) {
			return (SaveVectorFunc3D(f, SerializeFormatType::VECTOR_FIELD_3D_CARTESIAN, title, x1, x2, numPointsX, y1, y2, numPointsY, z1, z2,
															 numPointsZ, fileName, upper_threshold))
				.WithFile(fileName);
		}

		//===================================================================================
		// Spherical coordinate vector field serialization
		//===================================================================================

		inline IOResult SaveVectorFuncSpherical(const IVectorFunction<3>& f, std::string title, Real r1, Real r2, int numPointsR, Real theta1,
																						Real theta2, int numPointsTheta, Real phi1, Real phi2, int numPointsPhi, std::string fileName) {
			return (SaveVectorFunc3D(f, SerializeFormatType::VECTOR_FIELD_SPHERICAL, title, r1, r2, numPointsR, theta1, theta2, numPointsTheta,
															 phi1, phi2, numPointsPhi, fileName))
				.WithFile(fileName);
		}

		inline IOResult SaveVectorFuncSpherical(const IVectorFunction<3>& f, std::string title, Real r1, Real r2, int numPointsR, Real theta1,
																						Real theta2, int numPointsTheta, Real phi1, Real phi2, int numPointsPhi, std::string fileName,
																						Real upper_threshold) {
			return (SaveVectorFunc3D(f, SerializeFormatType::VECTOR_FIELD_SPHERICAL, title, r1, r2, numPointsR, theta1, theta2, numPointsTheta,
															 phi1, phi2, numPointsPhi, fileName, upper_threshold))
				.WithFile(fileName);
		}

	} // namespace Serializer
} // namespace MML

#endif // MML_SERIALIZER_VECTOR_FIELDS_H
