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
		namespace detail
		{
			inline IOResult ParticleSimulationFrameCount(std::size_t numSteps, int stride, int& frames) {
				if (stride < 1)
					return SerializeFailure(IOError::INVALID_PARAMETERS, "saveEveryNSteps must be >= 1");
				if (numSteps > static_cast<std::size_t>(std::numeric_limits<int>::max()))
					return SerializeFailure(IOError::INVALID_PARAMETERS, "Particle history length exceeds the supported int range");
				const auto stepStride = static_cast<std::size_t>(stride);
				frames = static_cast<int>(numSteps / stepStride + (numSteps % stepStride != 0));
				return SerializeSuccess();
			}

			template<class Point>
			IOResult ValidateParticleHistories(const std::vector<std::vector<Point>>& histories, int stride, int& frames) {
				if (histories.empty())
					return SerializeFailure(IOError::INVALID_PARAMETERS, "Particle histories must contain at least one particle");
				const auto numSteps = histories.front().size();
				for (std::size_t particle = 1; particle < histories.size(); ++particle)
					if (histories[particle].size() != numSteps)
						return SerializeFailure(IOError::INVALID_PARAMETERS,
							"Particle history " + std::to_string(particle) + " length must match the first particle history");
				return ParticleSimulationFrameCount(numSteps, stride, frames);
			}
		}

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
			if (ballPositions.size() != static_cast<size_t>(numBalls) || ballColors.size() != static_cast<size_t>(numBalls) ||
					ballRadius.size() != static_cast<size_t>(numBalls))
				return IOResult{false, IOError::INVALID_PARAMETERS, "ballPositions, ballColors, and ballRadius sizes must equal numBalls"}.WithFile(
					fileName);
			const auto catalogValidation = detail::ValidateParticleCatalog(numBalls, ballColors, ballRadius);
			if (!catalogValidation.success) return catalogValidation.WithFile(fileName);
			const auto pathValidation = detail::ValidateFileName(fileName);
			if (!pathValidation.success) return pathValidation;
			int actualFrames = 0;
			const auto validation = detail::ValidateParticleHistories(ballPositions, saveEveryNSteps, actualFrames);
			if (!validation.success)
				return validation.WithFile(fileName);

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			try
			{
				std::ostringstream buffer;
				Serializer::detail::ApplyTextFormat(buffer);
				auto headerResult = WriteParticleSimulation2DHeader(buffer, numBalls, width, height, ballColors, ballRadius, actualFrames);
				if (!headerResult.success)
					return (headerResult).WithFile(fileName);

				for (int realStep = 0; realStep < actualFrames; ++realStep)
				{
					const auto i = static_cast<std::size_t>(realStep) * static_cast<std::size_t>(saveEveryNSteps);
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
			if (ballPositions.size() != static_cast<size_t>(numBalls) || ballColors.size() != static_cast<size_t>(numBalls) ||
					ballRadius.size() != static_cast<size_t>(numBalls))
				return IOResult{false, IOError::INVALID_PARAMETERS, "ballPositions, ballColors, and ballRadius sizes must equal numBalls"}.WithFile(
					fileName);
			const auto catalogValidation = detail::ValidateParticleCatalog(numBalls, ballColors, ballRadius);
			if (!catalogValidation.success) return catalogValidation.WithFile(fileName);
			const auto pathValidation = detail::ValidateFileName(fileName);
			if (!pathValidation.success) return pathValidation;
			int actualFrames = 0;
			const auto validation = detail::ValidateParticleHistories(ballPositions, saveEveryNSteps, actualFrames);
			if (!validation.success)
				return validation.WithFile(fileName);

			std::ofstream file(fileName);
			Serializer::detail::ApplyTextFormat(file);
			if (!file.is_open())
				return IOResult{false, IOError::FILE_NOT_OPENED, "Cannot open file: " + fileName, fileName}.WithFile(fileName);

			try
			{
				auto headerResult = WriteParticleSimulation3DHeader(file, numBalls, width, height, depth, ballColors, ballRadius, actualFrames);
				if (!headerResult.success)
					return (headerResult).WithFile(fileName);

				for (int realStep = 0; realStep < actualFrames; ++realStep)
				{
					const auto i = static_cast<std::size_t>(realStep) * static_cast<std::size_t>(saveEveryNSteps);
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
