///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        serializer/SerializerODE.h                                          ///
///  Description: ODE solution serialization utilities                                ///
///               Save ODE solutions as functions and parametric curves               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_SERIALIZER_ODE_H
#define MML_SERIALIZER_ODE_H

#include "mml/tools/serializer/SerializerBase.h"
#include "mml/base/ODESystemSolution.h"

namespace MML
{
	namespace Serializer
	{
		//===================================================================================
		// ODE solution serialization
		//===================================================================================

		/// @brief Serialize a single component of an ODE solution as a real function
		/// @param sol The ODE system solution
		/// @param compInd Component index to serialize (0-based)
		/// @param title Display title for the function
		/// @param fileName Output file path
		/// @return SerializeResult with success flag and error details
		inline SerializeResult SaveODESolutionComponentAsFunc(const ODESystemSolution& sol, int compInd, 
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

		/// @brief Serialize all components of an ODE solution as a multi-function
		/// @param sol The ODE system solution
		/// @param title Display title for the multi-function
		/// @param legend Labels for each component
		/// @param fileName Output file path
		/// @return SerializeResult with success flag and error details
		inline SerializeResult SaveODESolutionAsMultiFunc(const ODESystemSolution& sol, std::string title, 
		                                                  std::vector<std::string> legend, std::string fileName)
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

		/// @brief Serialize two components of an ODE solution as a 2D parametric curve
		/// @param sol The ODE system solution
		/// @param fileName Output file path
		/// @param ind1 Index of first component (x-axis)
		/// @param ind2 Index of second component (y-axis)
		/// @param title Display title for the curve
		/// @return SerializeResult with success flag and error details
		inline SerializeResult SaveODESolAsParametricCurve2D(const ODESystemSolution& sol, std::string fileName, 
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
		
		/// @brief Serialize three components of an ODE solution as a 3D parametric curve
		/// @param sol The ODE system solution
		/// @param fileName Output file path
		/// @param ind1 Index of first component (x-axis)
		/// @param ind2 Index of second component (y-axis)
		/// @param ind3 Index of third component (z-axis)
		/// @param title Display title for the curve
		/// @return SerializeResult with success flag and error details
		inline SerializeResult SaveODESolAsParametricCurve3D(const ODESystemSolution& sol, std::string fileName, 
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

	} // namespace Serializer
} // namespace MML

#endif // MML_SERIALIZER_ODE_H
