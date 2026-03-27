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

#include "mml/tools/serializer/SerializerBase.h"
#include "mml/base/Geometry/Geometry.h"

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
		/// @return SerializeResult with success flag and error details
		inline SerializeResult SaveParticleSimulation2D(std::string fileName, int numBalls, Real width, Real height,
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
			if (saveEveryNSteps < 1)
				return {false, SerializeError::INVALID_PARAMETERS, "saveEveryNSteps must be >= 1"};

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
				int actualFrames = (numSteps + saveEveryNSteps - 1) / saveEveryNSteps;
				buffer << "NumSteps: " << actualFrames << std::endl;

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
		/// @return SerializeResult with success flag and error details
		inline SerializeResult SaveParticleSimulation3D(std::string fileName, int numBalls, Real width, Real height, Real depth,
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
			if (saveEveryNSteps < 1)
				return {false, SerializeError::INVALID_PARAMETERS, "saveEveryNSteps must be >= 1"};

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
				int actualFrames = (numSteps + saveEveryNSteps - 1) / saveEveryNSteps;
				file << "NumSteps: " << actualFrames << std::endl;

				int realStep = 0;
				for (int i = 0; i < numSteps; i += saveEveryNSteps, realStep++)
				{
					file << "Step " << realStep << " " << i * dT << std::endl;
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

	} // namespace Serializer
} // namespace MML

#endif // MML_SERIALIZER_SIMULATION_H
