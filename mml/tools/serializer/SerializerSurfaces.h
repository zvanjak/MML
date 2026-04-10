///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        serializer/SerializerSurfaces.h                                     ///
///  Description: Surface and scalar function serialization utilities                 ///
///               Save parametric surfaces and 2D/3D scalar functions to files        ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_SERIALIZER_SURFACES_H
#define MML_SERIALIZER_SURFACES_H

#include "mml/tools/serializer/SerializerBase.h"
#include "mml/interfaces/IFunction.h"
#include "mml/base/Vector/VectorN.h"

namespace MML
{
	namespace Serializer
	{
		//===================================================================================
		// Parametric Surface Serialization
		//===================================================================================

		/// @brief Serialize a 3D parametric surface to file
		/// @details Saves a parametric surface r(u,w) = (x(u,w), y(u,w), z(u,w)) evaluated on a 
		/// rectangular parameter grid to file. Each point contains parameter values and 3D position.
		/// Format: Header with grid parameters, then rows of [u, w, x, y, z].
		/// Ideal for 3D surface visualization of spheres, tori, and other parametric shapes.
		/// @param surface The parametric surface function (R² → R³)
		/// @param title Display title for the surface
		/// @param u1 Start of u parameter range
		/// @param u2 End of u parameter range
		/// @param numPointsU Number of grid points along u (must be >= 2)
		/// @param w1 Start of w parameter range
		/// @param w2 End of w parameter range
		/// @param numPointsW Number of grid points along w (must be >= 2)
		/// @param fileName Output file path
		/// @return SerializeResult with success flag and error details
		inline SerializeResult SaveParametricSurface(const IParametricSurface<3>& surface, std::string title,
		                                             Real u1, Real u2, int numPointsU,
		                                             Real w1, Real w2, int numPointsW,
		                                             std::string fileName)
		{
			// Validate parameters
			if (numPointsU < 2 || numPointsW < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "numPointsU and numPointsW must be >= 2"};
			if (u1 >= u2 || w1 >= w2)
				return {false, SerializeError::INVALID_PARAMETERS, "u1 < u2 and w1 < w2 required"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				// Write header
				file << FormatType::PARAMETRIC_SURFACE_CARTESIAN << std::endl;
				file << "VERSION: " << FormatType::CURRENT_VERSION << std::endl;
				file << title << std::endl;
				file << "u1: " << u1 << std::endl;
				file << "u2: " << u2 << std::endl;
				file << "NumPointsU: " << numPointsU << std::endl;
				file << "w1: " << w1 << std::endl;
				file << "w2: " << w2 << std::endl;
				file << "NumPointsW: " << numPointsW << std::endl;

				// Write grid data
				Real stepU = (u2 - u1) / (numPointsU - 1);
				Real stepW = (w2 - w1) / (numPointsW - 1);
				
				for (int i = 0; i < numPointsU; i++)
				{
					Real u = u1 + i * stepU;
					for (int j = 0; j < numPointsW; j++)
					{
						Real w = w1 + j * stepW;
						VectorN<Real, 3> pos = surface(u, w);
						file << u << " " << w << " " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
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

		/// @brief Serialize a general R² → R³ function as a parametric surface
		/// @details Lower-level function that works with any IVectorFunctionNM<2,3>.
		/// Useful when you have a lambda or custom function that isn't wrapped in IParametricSurface.
		inline SerializeResult SaveParametricSurface(const IVectorFunctionNM<2, 3>& f, std::string title,
		                                             Real u1, Real u2, int numPointsU,
		                                             Real w1, Real w2, int numPointsW,
		                                             std::string fileName)
		{
			// Validate parameters
			if (numPointsU < 2 || numPointsW < 2)
				return {false, SerializeError::INVALID_PARAMETERS, "numPointsU and numPointsW must be >= 2"};
			if (u1 >= u2 || w1 >= w2)
				return {false, SerializeError::INVALID_PARAMETERS, "u1 < u2 and w1 < w2 required"};
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				// Write header
				file << FormatType::PARAMETRIC_SURFACE_CARTESIAN << std::endl;
				file << "VERSION: " << FormatType::CURRENT_VERSION << std::endl;
				file << title << std::endl;
				file << "u1: " << u1 << std::endl;
				file << "u2: " << u2 << std::endl;
				file << "NumPointsU: " << numPointsU << std::endl;
				file << "w1: " << w1 << std::endl;
				file << "w2: " << w2 << std::endl;
				file << "NumPointsW: " << numPointsW << std::endl;

				// Write grid data
				Real stepU = (u2 - u1) / (numPointsU - 1);
				Real stepW = (w2 - w1) / (numPointsW - 1);
				
				for (int i = 0; i < numPointsU; i++)
				{
					Real u = u1 + i * stepU;
					for (int j = 0; j < numPointsW; j++)
					{
						Real w = w1 + j * stepW;
						VectorN<Real, 3> pos = f(VectorN<Real, 2>{u, w});
						file << u << " " << w << " " << pos[0] << " " << pos[1] << " " << pos[2] << std::endl;
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

		/// @brief Serialize a 3D parametric surface using its intrinsic bounds
		/// @details Convenience overload that uses the surface's getMinU/MaxU/MinW/MaxW methods
		/// to determine the parameter range automatically.
		inline SerializeResult SaveParametricSurface(const IParametricSurfaceRect<3>& surface, std::string title,
		                                             int numPointsU, int numPointsW,
		                                             std::string fileName)
		{
			return SaveParametricSurface(
				static_cast<const IVectorFunctionNM<2, 3>&>(surface), title,
				surface.getMinU(), surface.getMaxU(), numPointsU,
				surface.getMinW(), surface.getMaxW(), numPointsW,
				fileName);
		}

		//===================================================================================
		// Scalar function serialization
		//===================================================================================

		/// @brief Serialize a 2D scalar function (z=f(x,y)) to a grid file
		/// @details Saves a 2D scalar function evaluated on a rectangular grid to file.
		/// Format: Header with grid parameters, then value rows (one per y-line, varying x).
		/// Ideal for surface plots and heatmaps.
		inline SerializeResult SaveScalarFunc2DCartesian(const IScalarFunction<2>& f, std::string title,
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
				file << FormatType::SCALAR_FUNCTION_CARTESIAN_2D << std::endl;
				file << "VERSION: " << FormatType::CURRENT_VERSION << std::endl;
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

		/// @brief Serialize a 3D scalar function (w=f(x,y,z)) to a grid file
		/// @details Saves a 3D scalar function evaluated on a rectangular 3D grid to file.
		/// Format: Header with grid parameters, then value rows (one per point in 3D grid).
		/// Useful for volumetric data and 3D scalar field visualization.
		inline SerializeResult SaveScalarFunc3DCartesian(const IScalarFunction<3>& f, std::string title, 
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
				file << FormatType::SCALAR_FUNCTION_CARTESIAN_3D << std::endl;
				file << "VERSION: " << FormatType::CURRENT_VERSION << std::endl;
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

		//===================================================================================
		// Grid function serialization
		//===================================================================================

		/// @brief Serialize a 2D grid function (e.g., PDE solution) directly to surface visualization format
		/// @details Saves grid data without the overhead of wrapping in IScalarFunction<2>.
		/// Output format is identical to SaveScalarFunc2DCartesian, allowing direct
		/// visualization with the surface visualizer.
		///
		/// @tparam GridFunc Type with grid() method returning a grid with domain(), numNodesX/Y(), x(i), y(j)
		/// and operator()(i,j) for value access
		/// @param gridFunc The grid function to serialize (e.g., MML::PDE::GridFunction2D)
		/// @param title Display title for the surface plot
		/// @param fileName Output file path
		/// @param scaleXY Scale factor for x,y coordinates (default 1.0)
		/// @param scaleValue Scale factor for values (default 1.0)
		/// @return SerializeResult with success flag and error details
		template<typename GridFunc>
		inline SerializeResult SaveGridFunction2D(const GridFunc& gridFunc, std::string title,
		                                          std::string fileName,
		                                          Real scaleXY = 1.0, Real scaleValue = 1.0)
		{
			if (fileName.empty())
				return {false, SerializeError::INVALID_PARAMETERS, "fileName cannot be empty"};

			std::ofstream file(fileName);
			if (!file.is_open())
				return {false, SerializeError::FILE_NOT_OPENED, "Cannot open file: " + fileName};

			try
			{
				const auto& grid = gridFunc.grid();
				const auto& domain = grid.domain();
				
				int numPointsX = grid.numNodesX();
				int numPointsY = grid.numNodesY();
				
				// Apply scaling to domain bounds
				Real x1 = domain.xMin() * scaleXY;
				Real x2 = domain.xMax() * scaleXY;
				Real y1 = domain.yMin() * scaleXY;
				Real y2 = domain.yMax() * scaleXY;
				
				file << FormatType::SCALAR_FUNCTION_CARTESIAN_2D << std::endl;
				file << "VERSION: " << FormatType::CURRENT_VERSION << std::endl;
				file << title << std::endl;
				file << "x1: " << x1 << std::endl;
				file << "x2: " << x2 << std::endl;
				file << "NumPointsX: " << numPointsX << std::endl;
				file << "y1: " << y1 << std::endl;
				file << "y2: " << y2 << std::endl;
				file << "NumPointsY: " << numPointsY << std::endl;

				// Output grid values directly (no interpolation needed)
				for (int i = 0; i < numPointsX; i++)
				{
					for (int j = 0; j < numPointsY; j++)
					{
						Real x = grid.x(i) * scaleXY;
						Real y = grid.y(j) * scaleXY;
						Real value = gridFunc(i, j) * scaleValue;
						file << x << " " << y << " " << value << std::endl;
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

	} // namespace Serializer
} // namespace MML

#endif // MML_SERIALIZER_SURFACES_H
