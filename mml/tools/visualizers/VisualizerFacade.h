///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///  File:        tools/visualizers/VisualizerFacade.h
///  Description: Public visualization workflows over shared configuration, files and processes
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                         ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_VISUALIZER_FACADE_H
#define MML_VISUALIZER_FACADE_H

#include <mml/MMLBase.h>
#include <mml/interfaces/IFunction.h>
#include <mml/base/InterpolatedFunction.h>
#include <mml/base/ODESystemSolution.h>
#include <mml/tools/Serializer.h>
#include <mml/tools/visualizers/VisualizerConfig.h>
#include <mml/tools/visualizers/VisualizerFiles.h>
#include <mml/tools/visualizers/VisualizerProcess.h>

#include <string>
#include <vector>

namespace MML {
	class Visualizer {
		using Kind = detail::visualizer::Kind;
		static inline std::string _pathResultFiles{GetResultFilesPath()};

		static IOResult Launch(Kind kind, const std::vector<std::string>& paths) {
			if (paths.empty())
				return IOResult::Failure("No visualizer data files specified", {}, -1, IOError::INVALID_PARAMETERS);
			return detail::visualizer::ExecuteVisualizer(detail::visualizer::ExecutablePath(kind), paths);
		}

		template<class Save>
		static IOResult SaveAndLaunch(Kind kind, const std::string& filename, Save save) {
			const auto output = detail::visualizer::MakeOutputPath(filename, _pathResultFiles);
			if (!output)
				return output;
			const auto saved = save(output.filename);
			if (!saved)
				return saved.WithFile(output.filename);
			return Launch(kind, {output.filename});
		}

		static IOResult LaunchFromFiles(Kind kind, const std::vector<std::string>& filenames) {
			std::vector<std::string> paths;
			const auto result = detail::visualizer::ResolveInputFiles(filenames, _pathResultFiles, paths);
			if (!result)
				return result;
			return Launch(kind, paths);
		}

	public:
		/// @brief Visualize a real function over an interval
		static IOResult VisualizeRealFunction(const IRealFunction& f, std::string title, Real x1, Real x2, int numPoints,
																					std::string fileName) {
			return SaveAndLaunch(Kind::RealFunction, fileName, [&](const std::string& name) {
				return Serializer::SaveRealFunc(f, title, x1, x2, numPoints, name);
			});
		}
		/// @brief Visualize a real function at specified points
		static IOResult VisualizeRealFunction(const IRealFunction& f, std::string title, Vector<Real> points, std::string fileName) {
			return SaveAndLaunch(Kind::RealFunction, fileName, [&](const std::string& name) {
				return Serializer::SaveRealFunc(f, title, points, name);
			});
		}
		/// @brief Visualize multiple real functions on the same plot
		static IOResult VisualizeMultiRealFunction(std::vector<IRealFunction*> funcs, std::string title, std::vector<std::string> func_legend,
																							 Real x1, Real x2, int numPoints, std::string fileName) {
			return SaveAndLaunch(Kind::RealFunction, fileName, [&](const std::string& name) {
				return Serializer::SaveRealMultiFunc(funcs, title, func_legend, x1, x2, numPoints, name);
			});
		}
		/// @brief Visualize multiple linear interpolation functions on the same plot
		static IOResult VisualizeMultiRealFunction(std::vector<LinearInterpRealFunc> funcs, std::string title,
																							 std::vector<std::string> func_legend, Real x1, Real x2, int numPoints,
																							 std::string fileName) {
			return SaveAndLaunch(Kind::RealFunction, fileName, [&](const std::string& name) {
				return Serializer::SaveRealMultiFunc(funcs, title, func_legend, x1, x2, numPoints, name);
			});
		}
		/// @brief Visualize multiple functions as separate plots
		static IOResult VisualizeMultiRealFunctionSeparately(std::vector<LinearInterpRealFunc> funcs, std::string title,
																												 std::vector<std::string> func_legend, Real x1, Real x2, int numPoints,
																												 std::string fileName) {
			// for each function, serialize data to file, with name generated from fileName
			std::vector<std::string> fileNames;
			int i = 0;
			for (auto& func : funcs) {
				std::string suffix = "_" + std::to_string(i + 1) + ".mml";
				const auto output = detail::visualizer::MakeOutputPath(fileName, _pathResultFiles, suffix);
				if (!output)
					return output;
				const std::string& name = output.filename;

				Real xLow = x1 < funcs[i].MinX() ? funcs[i].MinX() : x1;
				Real xHigh = x2 > funcs[i].MaxX() ? funcs[i].MaxX() : x2;

				auto saveResult = Serializer::SaveRealFunc(funcs[i], func_legend[i], xLow, xHigh, numPoints, name);
				if (!saveResult.success) {
					return saveResult.WithFile(name);
				}
				fileNames.push_back(name);
				i++;
			}
			return Launch(Kind::RealFunction, fileNames);
		}
		/// @brief Visualize multiple polynomial interpolation functions on the same plot
		static IOResult VisualizeMultiRealFunction(std::vector<PolynomInterpRealFunc> funcs, std::string title,
																							 std::vector<std::string> func_legend, Real x1, Real x2, int numPoints,
																							 std::string fileName) {
			return SaveAndLaunch(Kind::RealFunction, fileName, [&](const std::string& name) {
				return Serializer::SaveRealMultiFunc(funcs, title, func_legend, x1, x2, numPoints, name);
			});
		}
		/// @brief Visualize multiple spline interpolation functions on the same plot
		static IOResult VisualizeMultiRealFunction(std::vector<SplineInterpRealFunc> funcs, std::string title,
																							 std::vector<std::string> func_legend, Real x1, Real x2, int numPoints,
																							 std::string fileName) {
			return SaveAndLaunch(Kind::RealFunction, fileName, [&](const std::string& name) {
				return Serializer::SaveRealMultiFunc(funcs, title, func_legend, x1, x2, numPoints, name);
			});
		}
		/// @brief Visualize a pre-saved real function data file
		/// @param fileNameOrPath Full path or filename (looked up in results folder as fallback)
		/// @return IOResult indicating success or failure
		static IOResult VisualizeRealFunctionFromFile(std::string fileNameOrPath) {
			return LaunchFromFiles(Kind::RealFunction, {fileNameOrPath});
		}

		/// @brief Visualize a 2D scalar function as a surface plot
		static IOResult VisualizeScalarFunc2DCartesian(const IScalarFunction<2>& func, std::string title, Real x1, Real x2, int numPointsX,
																									 Real y1, Real y2, int numPointsY, std::string fileName) {
			return SaveAndLaunch(Kind::ScalarFunction2D, fileName, [&](const std::string& name) {
				return Serializer::SaveScalarFunc2DCartesian(func, title, x1, x2, numPointsX, y1, y2, numPointsY, name);
			});
		}

		// Visualization of 2D grid functions (e.g., PDE solutions) - direct serialization without IScalarFunction wrapper
		/// @brief Visualize a 2D grid function (e.g., PDE solution) as a surface plot
		/// @details Directly serializes grid data without wrapping in IScalarFunction<2>.
		/// Ideal for visualizing PDE solutions from GridFunction2D.
		///
		/// @tparam GridFunc Type with grid() and operator()(i,j) interface
		/// @param gridFunc The grid function to visualize
		/// @param title Display title for the surface plot
		/// @param fileName Output file name (without path)
		/// @param scaleXY Scale factor for x,y coordinates (e.g., 100.0 to scale [0,1] to [0,100])
		/// @param scaleValue Scale factor for values (e.g., 100.0 for better visibility)
		/// @return IOResult with success/failure info
		///
		/// @example
		/// auto solution = poissonSolver.solve();
		/// Visualizer::VisualizeGridFunction2D(solution, "Temperature Distribution", "temp.mml", 100.0, 100.0);
		template<typename GridFunc>
		static IOResult VisualizeGridFunction2D(const GridFunc& gridFunc, std::string title, std::string fileName, Real scaleXY = 1.0,
																						Real scaleValue = 1.0) {
			return SaveAndLaunch(Kind::ScalarFunction2D, fileName, [&](const std::string& name) {
				return Serializer::SaveGridFunction2D(gridFunc, title, name, scaleXY, scaleValue);
			});
		}

		/// @brief Visualize a pre-saved 2D scalar function data file
		/// @param fileNameOrPath Full path or filename (looked up in results folder as fallback)
		/// @return IOResult indicating success or failure
		static IOResult VisualizeScalarFunc2DFromFile(std::string fileNameOrPath) {
			return LaunchFromFiles(Kind::ScalarFunction2D, {fileNameOrPath});
		}

		// visualizations of Vector fields
		static IOResult VisualizeVectorField2DCartesian(const IVectorFunction<2>& func, std::string title, Real x1, Real x2, int numPointsX,
																										Real y1, Real y2, int numPointsY, std::string fileName) {
			return SaveAndLaunch(Kind::VectorField2D, fileName, [&](const std::string& name) {
				return Serializer::SaveVectorFunc2DCartesian(func, title, x1, x2, numPointsX, y1, y2, numPointsY, name);
			});
		}

		static IOResult VisualizeVectorField3DCartesian(const IVectorFunction<3>& func, std::string title, Real x1, Real x2, int numPointsX,
																										Real y1, Real y2, int numPointsY, Real z1, Real z2, int numPointsZ,
																										std::string fileName) {
			return SaveAndLaunch(Kind::VectorField3D, fileName, [&](const std::string& name) {
				return Serializer::SaveVectorFunc3DCartesian(func, title, x1, x2, numPointsX, y1, y2, numPointsY, z1, z2, numPointsZ, name, 3.0);
			});
		}

		/// @brief Visualize a pre-saved 2D vector field data file
		/// @param fileNameOrPath Full path or filename (looked up in results folder as fallback)
		/// @return IOResult indicating success or failure
		static IOResult VisualizeVectorField2DFromFile(std::string fileNameOrPath) {
			return LaunchFromFiles(Kind::VectorField2D, {fileNameOrPath});
		}

		/// @brief Visualize a pre-saved 3D vector field data file
		/// @param fileNameOrPath Full path or filename (looked up in results folder as fallback)
		/// @return IOResult indicating success or failure
		static IOResult VisualizeVectorField3DFromFile(std::string fileNameOrPath) {
			return LaunchFromFiles(Kind::VectorField3D, {fileNameOrPath});
		}

		/// @brief Visualize a pre-saved MML world scene file
		/// @param fileNameOrPath Full path or filename (looked up in results folder as fallback)
		/// @return IOResult indicating success or failure
		static IOResult VisualizeWorldSceneFromFile(std::string fileNameOrPath) {
			return LaunchFromFiles(Kind::World, {fileNameOrPath});
		}

		// visualizations of Parametric curves
		static IOResult VisualizeParamCurve2D(const IRealToVectorFunction<2>& f, std::string title, Real t1, Real t2, int numPoints,
																					std::string fileName) {
			return SaveAndLaunch(Kind::Curve2D, fileName, [&](const std::string& name) {
				return Serializer::SaveParamCurveCartesian2DResult(f, title, t1, t2, numPoints, name);
			});
		}

		static IOResult VisualizeMultiParamCurve2D(std::vector<IRealToVectorFunction<2>*> curves, std::string title, Real t1, Real t2,
																							 int numPoints, std::string fileName) {
			// for each function, serialize data to file, with name generated from fileName
			std::vector<std::string> fileNames;
			int i = 0;
			for (auto& curve : curves) {
				std::string suffix = "_" + std::to_string(i + 1) + ".mml";
				const auto output = detail::visualizer::MakeOutputPath(fileName, _pathResultFiles, suffix);
				if (!output)
					return output;
				const std::string& name = output.filename;
				auto saveResult = Serializer::SaveParamCurveCartesian2DResult(*curve, title, t1, t2, numPoints, name);
				if (!saveResult) {
					return saveResult.WithFile(name);
				}
				fileNames.push_back(name);
				i++;
			}
			return Launch(Kind::Curve2D, fileNames);
		}

		static IOResult VisualizeMultiParamCurve2D(std::vector<IRealToVectorFunction<2>*> curves, std::vector<std::string> legend, Real t1,
																							 Real t2, int numPoints, std::string fileName) {
			// for each function, serialize data to file, with name generated from fileName
			std::vector<std::string> fileNames;
			int i = 0;
			for (auto& curve : curves) {
				std::string suffix = "_" + std::to_string(i + 1) + ".mml";
				const auto output = detail::visualizer::MakeOutputPath(fileName, _pathResultFiles, suffix);
				if (!output)
					return output;
				const std::string& name = output.filename;
				auto saveResult = Serializer::SaveParamCurveCartesian2DResult(*curve, legend[i], t1, t2, numPoints, name);
				if (!saveResult) {
					return saveResult.WithFile(name);
				}
				fileNames.push_back(name);
				i++;
			}
			return Launch(Kind::Curve2D, fileNames);
		}

		static IOResult VisualizeMultiParamCurve2D(std::vector<std::string> fileNames) {
			return LaunchFromFiles(Kind::Curve2D, fileNames);
		}

		static IOResult VisualizeParamCurve3D(const IRealToVectorFunction<3>& f, std::string title, Real t1, Real t2, int numPoints,
																					std::string fileName) {
			return SaveAndLaunch(Kind::Curve3D, fileName, [&](const std::string& name) {
				return Serializer::SaveParamCurveCartesian3DResult(f, title, t1, t2, numPoints, name);
			});
		}

		static IOResult VisualizeMultiParamCurve3D(std::vector<IRealToVectorFunction<3>*> curves, std::string title, Real t1, Real t2,
																							 int numPoints, std::string fileName) {
			// for each function, serialize data to file, with name generated from fileName
			std::vector<std::string> fileNames;
			int i = 0;
			for (auto& curve : curves) {
				std::string suffix = "_" + std::to_string(i + 1) + ".mml";
				const auto output = detail::visualizer::MakeOutputPath(fileName, _pathResultFiles, suffix);
				if (!output)
					return output;
				const std::string& name = output.filename;
				auto saveResult = Serializer::SaveParamCurveCartesian3DResult(*curve, title, t1, t2, numPoints, name);
				if (!saveResult) {
					return saveResult.WithFile(name);
				}
				fileNames.push_back(name);
				i++;
			}
			return Launch(Kind::Curve3D, fileNames);
		}

		static IOResult VisualizeMultiParamCurve3D(std::vector<std::string> fileNames) {
			return LaunchFromFiles(Kind::Curve3D, fileNames);
		}

		// ODE Solution visualizations
		// Visualizing single variable of ODE system solution as Real function
		static IOResult VisualizeODESysSolCompAsFunc(const ODESystemSolution& sol, int compInd, std::string title, std::string fileName) {
			return SaveAndLaunch(Kind::RealFunction, fileName, [&](const std::string& name) {
				return Serializer::SaveODESolutionComponentAsFunc(sol, compInd, title, name);
			});
		}

		// Visualizing ODE system solution as a multi-function (all variables)
		static IOResult VisualizeODESysSolAsMultiFunc(const ODESystemSolution& sol, std::string title, std::vector<std::string> legend,
																									std::string fileName) {
			return SaveAndLaunch(Kind::RealFunction, fileName, [&](const std::string& name) {
				return Serializer::SaveODESolutionAsMultiFunc(sol, title, legend, name);
			});
		}

		// Visualizing two variables of ODE system solution as a parametric curve in 2D
		static IOResult VisualizeODESysSolAsParamCurve2(const ODESystemSolution& sol, int ind1, int ind2, std::string title,
																										std::string fileName) {
			return SaveAndLaunch(Kind::Curve2D, fileName, [&](const std::string& name) {
				return Serializer::SaveODESolAsParametricCurve2D(sol, name, ind1, ind2, title);
			});
		}

		// Visualizing three variables of ODE system solution as a parametric curve in 3D
		static IOResult VisualizeODESysSolAsParamCurve3(const ODESystemSolution& sol, int ind1, int ind2, int ind3, std::string title,
																										std::string fileName) {
			return SaveAndLaunch(Kind::Curve3D, fileName, [&](const std::string& name) {
				return Serializer::SaveODESolAsParametricCurve3D(sol, name, ind1, ind2, ind3, title);
			});
		}

		// Particle simulation visualizations
		static IOResult VisualizeParticleSimulation2D(std::string fileNameOrPath) {
			return LaunchFromFiles(Kind::Particles2D, {fileNameOrPath});
		}

		static IOResult VisualizeParticleSimulation3D(std::string fileNameOrPath) {
			return LaunchFromFiles(Kind::Particles3D, {fileNameOrPath});
		}

		// Parametric surface visualizations
		/// @brief Visualize a 3D parametric surface over a rectangular parameter domain
		/// @details Samples the surface on a uniform grid in (u, w) parameters and launches
		/// the parametric surface visualizer to display the result.
		/// @param surface The parametric surface to visualize
		/// @param title Title for the visualization
		/// @param u1 Minimum u parameter value
		/// @param u2 Maximum u parameter value
		/// @param numPointsU Number of sample points along u
		/// @param w1 Minimum w parameter value
		/// @param w2 Maximum w parameter value
		/// @param numPointsW Number of sample points along w
		/// @param fileName Output filename (without path)
		/// @return IOResult indicating success or failure
		static IOResult VisualizeParametricSurface(const IParametricSurface<3>& surface, std::string title, Real u1, Real u2, int numPointsU,
																							 Real w1, Real w2, int numPointsW, std::string fileName) {
			return SaveAndLaunch(Kind::ParametricSurface, fileName, [&](const std::string& name) {
				return Serializer::SaveParametricSurface(surface, title, u1, u2, numPointsU, w1, w2, numPointsW, name);
			});
		}

		/// @brief Visualize a 3D parametric surface using its intrinsic bounds
		/// @details For surfaces with rectangular domains (IParametricSurfaceRect), uses
		/// the surface's getMinU/getMaxU/getMinW/getMaxW bounds automatically.
		/// @param surface The parametric surface with rectangular domain
		/// @param title Title for the visualization
		/// @param numPointsU Number of sample points along u
		/// @param numPointsW Number of sample points along w
		/// @param fileName Output filename (without path)
		/// @return IOResult indicating success or failure
		static IOResult VisualizeParametricSurface(const IParametricSurfaceRect<3>& surface, std::string title, int numPointsU, int numPointsW,
																							 std::string fileName) {
			return SaveAndLaunch(Kind::ParametricSurface, fileName, [&](const std::string& name) {
				return Serializer::SaveParametricSurface(surface, title, numPointsU, numPointsW, name);
			});
		}

		/// @brief Visualize a general R² → R³ vector function as a parametric surface
		/// @param f The vector function mapping (u, w) to (x, y, z)
		/// @param title Title for the visualization
		/// @param u1 Minimum u parameter value
		/// @param u2 Maximum u parameter value
		/// @param numPointsU Number of sample points along u
		/// @param w1 Minimum w parameter value
		/// @param w2 Maximum w parameter value
		/// @param numPointsW Number of sample points along w
		/// @param fileName Output filename (without path)
		/// @return IOResult indicating success or failure
		static IOResult VisualizeParametricSurface(const IVectorFunctionNM<2, 3>& f, std::string title, Real u1, Real u2, int numPointsU,
																							 Real w1, Real w2, int numPointsW, std::string fileName) {
			return SaveAndLaunch(Kind::ParametricSurface, fileName, [&](const std::string& name) {
				return Serializer::SaveParametricSurface(f, title, u1, u2, numPointsU, w1, w2, numPointsW, name);
			});
		}

		/// @brief Visualize a pre-saved parametric surface data file
		/// @param fileName Name of the data file (in results folder)
		/// @return IOResult indicating success or failure
		static IOResult VisualizeParametricSurfaceFromFile(std::string fileNameOrPath) {
			return LaunchFromFiles(Kind::ParametricSurface, {fileNameOrPath});
		}

		// 3D Scalar function visualizations (volumetric data)
		/// @brief Visualize a 3D scalar function as volumetric data
		/// @details Samples the function on a uniform 3D grid and launches the
		/// 3D scalar function visualizer for isosurface/volume rendering.
		/// @param func The 3D scalar function f(x, y, z) to visualize
		/// @param title Title for the visualization
		/// @param x1 Minimum x value
		/// @param x2 Maximum x value
		/// @param numPointsX Number of sample points along x axis
		/// @param y1 Minimum y value
		/// @param y2 Maximum y value
		/// @param numPointsY Number of sample points along y axis
		/// @param z1 Minimum z value
		/// @param z2 Maximum z value
		/// @param numPointsZ Number of sample points along z axis
		/// @param fileName Output filename (without path)
		/// @return IOResult indicating success or failure
		/// @note Large grids (e.g., 50x50x50 = 125,000 points) can create large data files
		static IOResult VisualizeScalarFunc3DCartesian(const IScalarFunction<3>& func, std::string title, Real x1, Real x2, int numPointsX,
																									 Real y1, Real y2, int numPointsY, Real z1, Real z2, int numPointsZ,
																									 std::string fileName) {
			return SaveAndLaunch(Kind::ScalarFunction3D, fileName, [&](const std::string& name) {
				return Serializer::SaveScalarFunc3DCartesian(func, title, x1, x2, numPointsX, y1, y2, numPointsY, z1, z2, numPointsZ, name);
			});
		}

		/// @brief Visualize a pre-saved 3D scalar function data file
		/// @param fileName Name of the data file (in results folder)
		/// @return IOResult indicating success or failure
		static IOResult VisualizeScalarFunc3DFromFile(std::string fileNameOrPath) {
			return LaunchFromFiles(Kind::ScalarFunction3D, {fileNameOrPath});
		}

		// Rigid body trajectory visualizations
		/// @brief Visualize a single rigid body trajectory file
		/// @details Launches the rigid body visualizer with the specified data file.
		/// The file should be in RIGID_BODY_TRAJECTORY_3D format (single body trajectory).
		/// Use this with files like rigid_body_box1.mml, rigid_body_box2.mml.
		/// @param fileName Name or path of the data file (.mml format)
		/// @return IOResult indicating success or failure
		static IOResult VisualizeRigidBodyTrajectory(std::string fileName) {
			return LaunchFromFiles(Kind::RigidBody, {fileName});
		}

		/// @brief Visualize multiple rigid body trajectory files together
		/// @details Useful for comparing trajectories side by side.
		/// Each file should be in RIGID_BODY_TRAJECTORY_3D format.
		/// @param fileNames Vector of file paths to visualize together
		/// @return IOResult indicating success or failure
		static IOResult VisualizeRigidBodyTrajectory(std::vector<std::string> fileNames) {
			return LaunchFromFiles(Kind::RigidBody, fileNames);
		}
	};
}

#endif // MML_VISUALIZER_FACADE_H
