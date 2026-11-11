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

#include <mml/tools/serializer/SerializerBase.h>
#include <mml/base/ODESystemSolution.h>

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
		/// @return IOResult with success flag and error details
		inline IOResult SaveODESolutionComponentAsFunc(const ODESystemSolution& sol, int compInd, std::string title, std::string fileName) {
			// Validate parameters
			if (sol.getSysDim() <= 0)
				return IOResult{false, IOError::INVALID_PARAMETERS, "ODE solution dimension must be positive"}.WithFile(fileName);
			if (compInd < 0 || compInd >= sol.getSysDim())
				return IOResult{false, IOError::INVALID_PARAMETERS, "component index out of range [0, " + std::to_string(sol.getSysDim() - 1) + "]"}
					.WithFile(fileName);
			if (fileName.empty())
				return IOResult{false, IOError::INVALID_PARAMETERS, "fileName cannot be empty"}.WithFile(fileName);

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			try
			{
				if (auto headerResult = WriteRealFuncHeader(file, SerializeFormatType::REAL_FUNCTION, title, sol.getT1(), sol.getT2(), sol.size()); !headerResult.success)
					return (headerResult).WithFile(fileName);

				for (int i = 0; i < sol.size(); i++)
				{
					file << sol.getTValue(i) << " " << sol.getXValue(i, compInd) << std::endl;
				}
				
				return Serializer::detail::FinishFile(file, fileName, "Success");
			}
			catch (const std::exception& e)
			{
				return IOResult{false, IOError::WRITE_FAILED, std::string("Write error: ") + e.what()}.WithFile(fileName);
			}
		}

		/// @brief Serialize all components of an ODE solution as a multi-function
		/// @param sol The ODE system solution
		/// @param title Display title for the multi-function
		/// @param legend Labels for each component
		/// @param fileName Output file path
		/// @return IOResult with success flag and error details
		inline IOResult SaveODESolutionAsMultiFunc(const ODESystemSolution& sol, std::string title, std::vector<std::string> legend,
																							 std::string fileName) {
			// Validate parameters
			if (legend.size() != sol.getSysDim())
				return IOResult{false, IOError::INVALID_PARAMETERS,
												"legend size (" + std::to_string(legend.size()) + ") must match ODE dimension (" + std::to_string(sol.getSysDim()) +
													")"}
					.WithFile(fileName);
			if (fileName.empty())
				return IOResult{false, IOError::INVALID_PARAMETERS, "fileName cannot be empty"}.WithFile(fileName);

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			try
			{
				if (auto headerResult = WriteRealMultiFuncHeader(file, title, sol.getSysDim(), legend, sol.getT1(), sol.getT2(), sol.size()); !headerResult.success)
					return (headerResult).WithFile(fileName);

				for (int i = 0; i < sol.size(); i++)
				{
					file << sol.getTValue(i) << " ";
					for (int j = 0; j < sol.getSysDim(); j++)
					{
						file << sol.getXValue(i, j) << " ";
					}
					file << std::endl;
				}
				return Serializer::detail::FinishFile(file, fileName, "Success");
			}
			catch (const std::exception& e)
			{
				return IOResult{false, IOError::WRITE_FAILED, std::string("Write error: ") + e.what()}.WithFile(fileName);
			}
		}

		/// @brief Serialize two components of an ODE solution as a 2D parametric curve
		/// @param sol The ODE system solution
		/// @param fileName Output file path
		/// @param ind1 Index of first component (x-axis)
		/// @param ind2 Index of second component (y-axis)
		/// @param title Display title for the curve
		/// @return IOResult with success flag and error details
		inline IOResult SaveODESolAsParametricCurve2D(const ODESystemSolution& sol, std::string fileName, int ind1, int ind2,
																									std::string title) {
			// Validate parameters
			if (sol.getSysDim() <= 0)
				return IOResult{false, IOError::INVALID_PARAMETERS, "ODE solution dimension must be positive"}.WithFile(fileName);
			if (ind1 < 0 || ind1 >= sol.getSysDim() || ind2 < 0 || ind2 >= sol.getSysDim())
				return IOResult{false, IOError::INVALID_PARAMETERS,
												"component indices out of range [0, " + std::to_string(sol.getSysDim() - 1) + "]"}
					.WithFile(fileName);
			if (fileName.empty())
				return IOResult{false, IOError::INVALID_PARAMETERS, "fileName cannot be empty"}.WithFile(fileName);

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			try
			{
				if (auto headerResult = WriteParamCurveHeader(file, SerializeFormatType::PARAMETRIC_CURVE_CARTESIAN_2D, title, sol.getT1(), sol.getT2(), sol.size()); !headerResult.success)
					return (headerResult).WithFile(fileName);

				for (int i = 0; i < sol.size(); i++)
				{
					file << sol.getTValue(i) << " ";
					file << sol.getXValue(i, ind1) << " ";
					file << sol.getXValue(i, ind2) << " ";
					file << std::endl;
				}
				
				return Serializer::detail::FinishFile(file, fileName, "Success");
			}
			catch (const std::exception& e)
			{
				return IOResult{false, IOError::WRITE_FAILED, std::string("Write error: ") + e.what()}.WithFile(fileName);
			}
		}

		/// @brief Serialize three components of an ODE solution as a 3D parametric curve
		/// @param sol The ODE system solution
		/// @param fileName Output file path
		/// @param ind1 Index of first component (x-axis)
		/// @param ind2 Index of second component (y-axis)
		/// @param ind3 Index of third component (z-axis)
		/// @param title Display title for the curve
		/// @return IOResult with success flag and error details
		inline IOResult SaveODESolAsParametricCurve3D(const ODESystemSolution& sol, std::string fileName, int ind1, int ind2, int ind3,
																									std::string title) {
			// Validate parameters
			if (sol.getSysDim() <= 0)
				return IOResult{false, IOError::INVALID_PARAMETERS, "ODE solution dimension must be positive"}.WithFile(fileName);
			if (ind1 < 0 || ind1 >= sol.getSysDim() || ind2 < 0 || ind2 >= sol.getSysDim() || ind3 < 0 || ind3 >= sol.getSysDim())
				return IOResult{false, IOError::INVALID_PARAMETERS,
												"component indices out of range [0, " + std::to_string(sol.getSysDim() - 1) + "]"}
					.WithFile(fileName);
			if (fileName.empty())
				return IOResult{false, IOError::INVALID_PARAMETERS, "fileName cannot be empty"}.WithFile(fileName);

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			try
			{
				if (auto headerResult = WriteParamCurveHeader(file, SerializeFormatType::PARAMETRIC_CURVE_CARTESIAN_3D, title, sol.getT1(), sol.getT2(), sol.size()); !headerResult.success)
					return (headerResult).WithFile(fileName);

				for (int i = 0; i < sol.size(); i++)
				{
					file << sol.getTValue(i) << " ";
					file << sol.getXValue(i, ind1) << " ";
					file << sol.getXValue(i, ind2) << " ";
					file << sol.getXValue(i, ind3) << " ";
					file << std::endl;
				}

				return Serializer::detail::FinishFile(file, fileName, "Success");
			}
			catch (const std::exception& e)
			{
				return IOResult{false, IOError::WRITE_FAILED, std::string("Write error: ") + e.what()}.WithFile(fileName);
			}
		}

	} // namespace Serializer
} // namespace MML

#endif // MML_SERIALIZER_ODE_H
