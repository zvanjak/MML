///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        serializer/SerializerSimulation.h                                   ///
///  Description: Particle simulation serialization utilities                         ///
///               Save 2D/3D particle simulations to files                            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_SERIALIZER_SIMULATION_H
#define MML_SERIALIZER_SIMULATION_H

#include <mml/tools/serializer/SerializerBase.h>
#include <mml/base/Geometry/Geometry.h>

namespace MML
{
	namespace Serializer
	{
		//===================================================================================
		// Particle simulation serialization
		//===================================================================================

		/// @brief Serialize a 2D particle simulation to file
		/// @param fileName Output file path
		/// @param numBalls Number of particles in the simulation
		/// @param width Width of the simulation domain
		/// @param height Height of the simulation domain
		/// @param ballPositions Position history for each particle [numBalls][numSteps]
		/// @param ballColors Color string for each particle
		/// @param ballRadius Radius for each particle
		/// @param dT Time step between samples
		/// @param saveEveryNSteps Save every Nth step (default: 1 = save all)
		/// @return IOResult with success flag and error details
		inline IOResult SaveParticleSimulation2D(std::string fileName, int numBalls, Real width, Real height,
																						 std::vector<std::vector<Pnt2Cart>> ballPositions, std::vector<std::string> ballColors,
																						 std::vector<Real> ballRadius, Real dT, int saveEveryNSteps = 1) {
			// Validate parameters
			if (numBalls <= 0)
				return IOResult{false, IOError::INVALID_PARAMETERS, "numBalls must be positive"}.WithFile(fileName);
			if (width <= 0 || height <= 0)
				return IOResult{false, IOError::INVALID_PARAMETERS, "width and height must be positive"}.WithFile(fileName);
			if (dT <= 0)
				return IOResult{false, IOError::INVALID_PARAMETERS, "dT must be positive"}.WithFile(fileName);
			if (ballPositions.size() != numBalls || ballColors.size() != numBalls || ballRadius.size() != numBalls)
				return IOResult{false, IOError::INVALID_PARAMETERS, "ballPositions, ballColors, and ballRadius sizes must equal numBalls"}.WithFile(
					fileName);
			if (fileName.empty())
				return IOResult{false, IOError::INVALID_PARAMETERS, "fileName cannot be empty"}.WithFile(fileName);
			if (saveEveryNSteps < 1)
				return IOResult{false, IOError::INVALID_PARAMETERS, "saveEveryNSteps must be >= 1"}.WithFile(fileName);

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			try
			{
				std::ostringstream buffer;
				Serializer::detail::ApplyTextFormat(buffer);
				int numSteps = ballPositions[0].size() ;
				int actualFrames = (numSteps + saveEveryNSteps - 1) / saveEveryNSteps;
				auto headerResult = WriteParticleSimulation2DHeader(buffer, numBalls, width, height, ballColors, ballRadius, actualFrames);
				if (!headerResult.success)
					return (headerResult).WithFile(fileName);

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
				return Serializer::detail::FinishFile(file, fileName, "Success");
			}
			catch (const std::exception& e)
			{
				return IOResult{false, IOError::WRITE_FAILED, std::string("Write error: ") + e.what()}.WithFile(fileName);
			}
		}

		/// @brief Serialize a 3D particle simulation to file
		/// @param fileName Output file path
		/// @param numBalls Number of particles in the simulation
		/// @param width Width of the simulation domain
		/// @param height Height of the simulation domain
		/// @param depth Depth of the simulation domain
		/// @param ballPositions Position history for each particle [numBalls][numSteps]
		/// @param ballColors Color string for each particle
		/// @param ballRadius Radius for each particle
		/// @param dT Time step between samples
		/// @param saveEveryNSteps Save every Nth step (default: 1 = save all)
		/// @return IOResult with success flag and error details
		inline IOResult SaveParticleSimulation3D(std::string fileName, int numBalls, Real width, Real height, Real depth,
																						 std::vector<std::vector<Pnt3Cart>> ballPositions, std::vector<std::string> ballColors,
																						 std::vector<Real> ballRadius, Real dT, int saveEveryNSteps = 1) {
			// Validate parameters
			if (numBalls <= 0)
				return IOResult{false, IOError::INVALID_PARAMETERS, "numBalls must be positive"}.WithFile(fileName);
			if (width <= 0 || height <= 0 || depth <= 0)
				return IOResult{false, IOError::INVALID_PARAMETERS, "width, height, and depth must be positive"}.WithFile(fileName);
			if (dT <= 0)
				return IOResult{false, IOError::INVALID_PARAMETERS, "dT must be positive"}.WithFile(fileName);
			if (ballPositions.size() != numBalls || ballColors.size() != numBalls || ballRadius.size() != numBalls)
				return IOResult{false, IOError::INVALID_PARAMETERS, "ballPositions, ballColors, and ballRadius sizes must equal numBalls"}.WithFile(
					fileName);
			if (fileName.empty())
				return IOResult{false, IOError::INVALID_PARAMETERS, "fileName cannot be empty"}.WithFile(fileName);
			if (saveEveryNSteps < 1)
				return IOResult{false, IOError::INVALID_PARAMETERS, "saveEveryNSteps must be >= 1"}.WithFile(fileName);

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			try
			{
				int numSteps = ballPositions[0].size();
				int actualFrames = (numSteps + saveEveryNSteps - 1) / saveEveryNSteps;
				auto headerResult = WriteParticleSimulation3DHeader(file, numBalls, width, height, depth, ballColors, ballRadius, actualFrames);
				if (!headerResult.success)
					return (headerResult).WithFile(fileName);

				int realStep = 0;
				for (int i = 0; i < numSteps; i += saveEveryNSteps, realStep++)
				{
					file << "Step " << realStep << " " << i * dT << std::endl;
					for (int j = 0; j < numBalls; j++)
					{
						file << j << " " << ballPositions[j][i].X() << " " << ballPositions[j][i].Y() << " " << ballPositions[j][i].Z() << "\n";
					}
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

#endif // MML_SERIALIZER_SIMULATION_H
