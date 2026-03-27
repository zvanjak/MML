///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        serializer/SerializerVectors.h                                      ///
///  Description: Vector field serialization utilities                               ///
///               Save 2D/3D vector fields to files                                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_SERIALIZER_VECTORS_H
#define MML_SERIALIZER_VECTORS_H

#include "mml/tools/serializer/SerializerBase.h"
#include "mml/interfaces/IFunction.h"
#include "mml/base/Vector/VectorN.h"

namespace MML
{
	namespace Serializer
	{
		//===================================================================================
		// 2D vector function serialization
		//===================================================================================

		inline SerializeResult SaveVectorFunc2D(const IVectorFunction<2>& f, std::string inType, std::string title,
		                                        Real x1_start, Real x1_end, int numPointsX1,
		                                        Real x2_start, Real x2_end, int numPointsX2, std::string fileName)
		{
			// Validate parameters
			if (numPointsX1 < 2 || numPointsX2 < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "numPointsX1 and numPointsX2 must be >= 2"};
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

		inline SerializeResult SaveVectorFunc2D(const IVectorFunction<2>& f, std::string inType, std::string title,
		                                        Real x1_start, Real x1_end, int numPointsX1,
		                                        Real x2_start, Real x2_end, int numPointsX2,
		                                        std::string fileName, Real upper_threshold)
		{
			// Validate parameters
			if (numPointsX1 < 2 || numPointsX2 < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "numPointsX1 and numPointsX2 must be >= 2"};
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

		inline SerializeResult SaveVectorFunc2DCartesian(const IVectorFunction<2>& f, std::string title,
		                                                 Real x1, Real x2, int numPointsX,
		                                                 Real y1, Real y2, int numPointsY, std::string fileName)
		{
			return SaveVectorFunc2D(f, "VECTOR_FIELD_2D_CARTESIAN", title, x1, x2, numPointsX, y1, y2, numPointsY, fileName);
		}

		inline SerializeResult SaveVectorFunc2DCartesian(const IVectorFunction<2>& f, std::string title,
		                                                 Real x1, Real x2, int numPointsX,
		                                                 Real y1, Real y2, int numPointsY, std::string fileName, 
		                                                 Real upper_threshold)
		{
			return SaveVectorFunc2D(f, "VECTOR_FIELD_2D_CARTESIAN", title, x1, x2, numPointsX, y1, y2, numPointsY, fileName, upper_threshold);
		}

		//===================================================================================
		// 3D vector function serialization
		//===================================================================================

		inline SerializeResult SaveVectorFunc3D(const IVectorFunction<3>& f, std::string inType, std::string title,
		                                        Real x1_start, Real x1_end, int numPointsX1, 
		                                        Real x2_start, Real x2_end, int numPointsX2, 
		                                        Real x3_start, Real x3_end, int numPointsX3, std::string fileName)
		{
			// Validate parameters
			if (numPointsX1 < 2 || numPointsX2 < 2 || numPointsX3 < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "all numPoints must be >= 2"};
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

		inline SerializeResult SaveVectorFunc3D(const IVectorFunction<3>& f, std::string inType, std::string title, 
		                                        Real x1_start, Real x1_end, int numPointsX1, 
		                                        Real x2_start, Real x2_end, int numPointsX2, 
		                                        Real x3_start, Real x3_end, int numPointsX3, 
		                                        std::string fileName, Real upper_threshold)
		{
			// Validate parameters
			if (numPointsX1 < 2 || numPointsX2 < 2 || numPointsX3 < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "all numPoints must be >= 2"};
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

		inline SerializeResult SaveVectorFunc3DCartesian(const IVectorFunction<3>& f, std::string title, 
		                                                 Real x1, Real x2, int numPointsX, 
		                                                 Real y1, Real y2, int numPointsY, 
		                                                 Real z1, Real z2, int numPointsZ, std::string fileName)
		{
			return SaveVectorFunc3D(f, "VECTOR_FIELD_3D_CARTESIAN", title, x1, x2, numPointsX, y1, y2, numPointsY, z1, z2, numPointsZ, fileName);
		}

		inline SerializeResult SaveVectorFunc3DCartesian(const IVectorFunction<3>& f, std::string title, 
		                                                 Real x1, Real x2, int numPointsX, 
		                                                 Real y1, Real y2, int numPointsY, 
		                                                 Real z1, Real z2, int numPointsZ, std::string fileName, 
		                                                 Real upper_threshold)
		{
			return SaveVectorFunc3D(f, "VECTOR_FIELD_3D_CARTESIAN", title, x1, x2, numPointsX, y1, y2, numPointsY, z1, z2, numPointsZ, fileName, upper_threshold);
		}

		//===================================================================================
		// Spherical coordinate vector field serialization
		//===================================================================================
		
		inline SerializeResult SaveVectorFuncSpherical(const IVectorFunction<3>& f, std::string title, 
		                                               Real r1, Real r2, int numPointsR, 
		                                               Real theta1, Real theta2, int numPointsTheta, 
		                                               Real phi1, Real phi2, int numPointsPhi, std::string fileName)
		{
			return SaveVectorFunc3D(f, "VECTOR_FIELD_SPHERICAL", title, r1, r2, numPointsR, theta1, theta2, numPointsTheta, phi1, phi2, numPointsPhi, fileName);
		}

		inline SerializeResult SaveVectorFuncSpherical(const IVectorFunction<3>& f, std::string title, 
		                                               Real r1, Real r2, int numPointsR, 
		                                               Real theta1, Real theta2, int numPointsTheta, 
		                                               Real phi1, Real phi2, int numPointsPhi, 
		                                               std::string fileName, Real upper_threshold)
		{
			return SaveVectorFunc3D(f, "VECTOR_FIELD_SPHERICAL", title, r1, r2, numPointsR, theta1, theta2, numPointsTheta, phi1, phi2, numPointsPhi, fileName, upper_threshold);
		}

	} // namespace Serializer
} // namespace MML

#endif // MML_SERIALIZER_VECTORS_H
