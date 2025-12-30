///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Visualizer.h                                                        ///
///  Description: Visualization tools for functions, curves, and data                 ///
///               Gnuplot integration for 2D/3D plotting                              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_VIZUALIZER_H
#define MML_VIZUALIZER_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/InterpolatedFunction.h"
#include "base/ODESystemSolution.h"

#include "tools/Serializer.h"

#include <filesystem>
#include <algorithm>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#else
#include <unistd.h>
#include <sys/wait.h>
#endif


namespace MML
{
	// Result struct for visualizer operations - provides error information
	struct VisualizerResult
	{
		bool success;
		int exitCode;
		std::string errorMessage;
		std::string dataFilePath;  // Path to the data file that was created
		
		static VisualizerResult Success(const std::string& dataPath = "") {
			return VisualizerResult{true, 0, "", dataPath};
		}
		
		static VisualizerResult Failure(const std::string& message, const std::string& dataPath = "", int code = -1) {
			return VisualizerResult{false, code, message, dataPath};
		}
		
		// Allow implicit conversion to bool for backward compatibility
		operator bool() const { return success; }
	};

	class Visualizer
	{
		// Use getter functions for cross-platform path resolution
		static inline std::string _pathResultFiles{ GetResultFilesPath() };

		static inline std::string _pathRealFuncViz{ GetRealFuncVisualizerPath() };
		static inline std::string _pathSurfaceViz{ GetSurfaceVisualizerPath() };
		static inline std::string _pathParametricCurve3DViz{ GetParametricCurve3DVisualizerPath() };
		static inline std::string _pathParametricCurve2DViz{ GetParametricCurve2DVisualizerPath() };
		static inline std::string _pathVectorField2DViz{ GetVectorField2DVisualizerPath() };
		static inline std::string _pathVectorField3DViz{ GetVectorField3DVisualizerPath() };

		static inline std::string _pathParticle2DViz{ GetParticle2DVisualizerPath() };
		static inline std::string _pathParticle3DViz{ GetParticle3DVisualizerPath() };

		// Helper: validate filename contains only safe characters
		static bool IsFilenameSafe(const std::string& filename) {
			for (char c : filename) {
				// Allow alphanumeric, underscore, hyphen, dot, path separators, and colon (for drive letters)
				if (!std::isalnum(static_cast<unsigned char>(c)) && 
					c != '_' && c != '-' && c != '.' && c != '/' && c != '\\' && c != ':') {
					return false;
				}
			}
			return true;
		}
		
		// Helper: sanitize filename by replacing unsafe characters
		static std::string SanitizeFilename(const std::string& filename) {
			std::string result = filename;
			for (char& c : result) {
				if (!std::isalnum(static_cast<unsigned char>(c)) && 
					c != '_' && c != '-' && c != '.' && c != '/' && c != '\\' && c != ':') {
					c = '_';
				}
			}
			return result;
		}

		// Helper: normalize path separators for the current platform
		static std::string NormalizePath(const std::string& path) {
			std::string result = path;
#ifdef _WIN32
			for (char& c : result) {
				if (c == '/') c = '\\';
			}
#else
			for (char& c : result) {
				if (c == '\\') c = '/';
			}
#endif
			return result;
		}

		// Safe process execution without shell - Windows implementation
#ifdef _WIN32
		static VisualizerResult ExecuteVisualizerWindows(const std::string& executable, 
														const std::vector<std::string>& args) {
			// Build command line: "exe" "arg1" "arg2" ...
			std::string cmdLine = "\"" + NormalizePath(executable) + "\"";
			for (const auto& arg : args) {
				cmdLine += " \"" + NormalizePath(arg) + "\"";
			}
			
			// Get working directory (directory containing the executable)
			std::string exePath = NormalizePath(executable);
			std::string workDir;
			auto lastSep = exePath.find_last_of("\\");
			if (lastSep != std::string::npos) {
				workDir = exePath.substr(0, lastSep);
			}
			
			STARTUPINFOA si = {};
			si.cb = sizeof(si);
			PROCESS_INFORMATION pi = {};
			
			// CreateProcess needs a modifiable string
			std::vector<char> cmdLineBuffer(cmdLine.begin(), cmdLine.end());
			cmdLineBuffer.push_back('\0');
			
			BOOL success = CreateProcessA(
				nullptr,                    // Application name (use command line)
				cmdLineBuffer.data(),       // Command line
				nullptr,                    // Process security attributes
				nullptr,                    // Thread security attributes
				FALSE,                      // Inherit handles
				0,                          // Creation flags
				nullptr,                    // Environment
				workDir.empty() ? nullptr : workDir.c_str(),  // Working directory
				&si,                        // Startup info
				&pi                         // Process info
			);
			
			if (!success) {
				DWORD error = GetLastError();
				return VisualizerResult::Failure(
					"Failed to start visualizer process. Error code: " + std::to_string(error),
					args.empty() ? "" : args[0],
					static_cast<int>(error)
				);
			}
			
			// Wait for GUI visualizer to close (use 0xFFFFFFFF for INFINITE timeout, compatible with all include orders)
			WaitForSingleObject(pi.hProcess, 0xFFFFFFFF);
			
			DWORD dwExitCode = 0;
			GetExitCodeProcess(pi.hProcess, &dwExitCode);
			
			CloseHandle(pi.hProcess);
			CloseHandle(pi.hThread);
			
			return VisualizerResult::Success(args.empty() ? "" : args[0]);
		}
#else
		// Safe process execution without shell - POSIX implementation
		static VisualizerResult ExecuteVisualizerPosix(const std::string& executable,
													  const std::vector<std::string>& args) {
			pid_t pid = fork();
			
			if (pid < 0) {
				return VisualizerResult::Failure("Failed to fork process",
					args.empty() ? "" : args[0], -1);
			}
			
			if (pid == 0) {
				// Child process
				std::vector<char*> argv;
				std::string exePath = NormalizePath(executable);
				argv.push_back(const_cast<char*>(exePath.c_str()));
				
				std::vector<std::string> normalizedArgs;
				for (const auto& arg : args) {
					normalizedArgs.push_back(NormalizePath(arg));
				}
				for (auto& arg : normalizedArgs) {
					argv.push_back(const_cast<char*>(arg.c_str()));
				}
				argv.push_back(nullptr);
				
				execv(exePath.c_str(), argv.data());
				
				// If execv returns, it failed
				_exit(127);
			}
			
			// Parent process - don't wait for GUI visualizers
			// Let the visualizer run independently
			return VisualizerResult::Success(args.empty() ? "" : args[0]);
		}
#endif

		// Main execution dispatcher
		static VisualizerResult ExecuteVisualizer(const std::string& executable,
												 const std::vector<std::string>& args) {
			// Validate executable path
			if (executable.empty()) {
				return VisualizerResult::Failure("Visualizer executable path is empty");
			}
			
			// Validate all arguments (filenames)
			for (const auto& arg : args) {
				if (!IsFilenameSafe(arg)) {
					return VisualizerResult::Failure(
						"Invalid characters in filename: " + arg,
						arg
					);
				}
			}

#ifdef _WIN32
			return ExecuteVisualizerWindows(executable, args);
#else
			return ExecuteVisualizerPosix(executable, args);
#endif
		}

	public:
		// visualizations of Real function
		static VisualizerResult VisualizeRealFunction(const IRealFunction& f, std::string title,
																			Real x1, Real x2, int numPoints, std::string fileName)
		{
			std::string name = _pathResultFiles + SanitizeFilename(fileName);
			auto saveResult = Serializer::SaveRealFunc(f, title, x1, x2, numPoints, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(_pathRealFuncViz, {name});
		}

		static VisualizerResult VisualizeRealFunction(const IRealFunction& f, std::string title,
																			Vector<Real> points, std::string fileName)
		{
			std::string name = _pathResultFiles + SanitizeFilename(fileName);
			auto saveResult = Serializer::SaveRealFunc(f, title, points, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(_pathRealFuncViz, {name});
		}

		// visualizations of multiple Real functions
		static VisualizerResult VisualizeMultiRealFunction(std::vector<IRealFunction*> funcs, std::string title,
																		 std::vector<std::string> func_legend,
																		 Real x1, Real x2, int numPoints, std::string fileName)
		{
			std::string name = _pathResultFiles + SanitizeFilename(fileName);
			auto saveResult = Serializer::SaveRealMultiFunc(funcs, title, func_legend, x1, x2, numPoints, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(_pathRealFuncViz, {name});
		}

		static VisualizerResult VisualizeMultiRealFunction(std::vector<LinearInterpRealFunc> funcs, std::string title,
																		 std::vector<std::string> func_legend,
																		 Real x1, Real x2, int numPoints, std::string fileName)
		{
			std::string name = _pathResultFiles + SanitizeFilename(fileName);
			auto saveResult = Serializer::SaveRealMultiFunc(funcs, title, func_legend, x1, x2, numPoints, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(_pathRealFuncViz, {name});
		}

		// visualization of multiple real functions where each function is saved (and visualized) separately
		static VisualizerResult VisualizeMultiRealFunctionSeparately(std::vector<LinearInterpRealFunc> funcs, std::string title,
																										 std::vector<std::string> func_legend,
																										 Real x1, Real x2, int numPoints, std::string fileName)
		{
			// for each function, serialize data to file, with name generated from fileName
			std::vector<std::string> fileNames;
			int i = 0;
			for (auto& func : funcs)
			{
				std::string name = _pathResultFiles + SanitizeFilename(fileName) + "_" + std::to_string(i + 1) + ".txt";

				Real xLow  = x1 < funcs[i].MinX() ? funcs[i].MinX() : x1;
				Real xHigh = x2 > funcs[i].MaxX() ? funcs[i].MaxX() : x2;

				auto saveResult = Serializer::SaveRealFunc(funcs[i], func_legend[i], xLow, xHigh, numPoints, name);
				if (!saveResult.success) {
					return VisualizerResult::Failure("Failed to save function " + std::to_string(i) + ": " + saveResult.message, name);
				}
				fileNames.push_back(name);
				i++;
			}
			return ExecuteVisualizer(_pathRealFuncViz, fileNames);
		}

		static VisualizerResult VisualizeMultiRealFunction(std::vector<PolynomInterpRealFunc> funcs, std::string title,
																		 std::vector<std::string> func_legend,
																		 Real x1, Real x2, int numPoints, std::string fileName)
		{
			std::string name = _pathResultFiles + SanitizeFilename(fileName);
			auto saveResult = Serializer::SaveRealMultiFunc(funcs, title, func_legend, x1, x2, numPoints, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(_pathRealFuncViz, {name});
		}

		static VisualizerResult VisualizeMultiRealFunction(std::vector<SplineInterpRealFunc> funcs, std::string title,
																		 std::vector<std::string> func_legend,
																		 Real x1, Real x2, int numPoints, std::string fileName)
		{
			std::string name = _pathResultFiles + SanitizeFilename(fileName);
			auto saveResult = Serializer::SaveRealMultiFunc(funcs, title, func_legend, x1, x2, numPoints, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(_pathRealFuncViz, {name});
		}

		// visualizations of Scalar function in 2D
		static VisualizerResult VisualizeScalarFunc2DCartesian(const IScalarFunction<2>& func, std::string title,
																							 Real x1, Real x2, int numPointsX,
																							 Real y1, Real y2, int numPointsY, std::string fileName)
		{
			std::string name = _pathResultFiles + SanitizeFilename(fileName);
			auto saveResult = Serializer::SaveScalarFunc2DCartesian(func, title, x1, x2, numPointsX, y1, y2, numPointsY, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(_pathSurfaceViz, {name});
		}

		// visualizations of Vector fields
		static VisualizerResult VisualizeVectorField2DCartesian(const IVectorFunction<2>& func, std::string title,
																								Real x1, Real x2, int numPointsX,
																								Real y1, Real y2, int numPointsY, std::string fileName)
		{
			std::string name = _pathResultFiles + SanitizeFilename(fileName);
			auto saveResult = Serializer::SaveVectorFunc2DCartesian(func, title, x1, x2, numPointsX, y1, y2, numPointsY, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(_pathVectorField2DViz, {name});
		}

		static VisualizerResult VisualizeVectorField3DCartesian(const IVectorFunction<3>& func, std::string title,
																								Real x1, Real x2, int numPointsX,
																								Real y1, Real y2, int numPointsY,
																								Real z1, Real z2, int numPointsZ, std::string fileName)
		{
			std::string name = _pathResultFiles + SanitizeFilename(fileName);
			auto saveResult = Serializer::SaveVectorFunc3DCartesian(func, title, x1, x2, numPointsX, y1, y2, numPointsY, z1, z2, numPointsZ, name, 3.0);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(_pathVectorField3DViz, {name});
		}

		// visualizations of Parametric curves
		static VisualizerResult VisualizeParamCurve2D(const IRealToVectorFunction<2>& f, std::string title, 
																			Real t1, Real t2, int numPoints, 
																			std::string fileName)
		{
			std::string name = _pathResultFiles + SanitizeFilename(fileName);
			bool saved = Serializer::SaveParamCurveCartesian2D(f, title, t1, t2, numPoints, name);
			if (!saved) {
				return VisualizerResult::Failure("Failed to save parametric curve data", name);
			}
			return ExecuteVisualizer(_pathParametricCurve2DViz, {name});
		}

		static VisualizerResult VisualizeMultiParamCurve2D(std::vector<IRealToVectorFunction<2>*> curves, 
																					 std::string title,
																					 Real t1, Real t2,
																					 int numPoints, std::string fileName)
		{
			// for each function, serialize data to file, with name generated from fileName
			std::vector<std::string> fileNames;
			int i = 0;
			for (auto& curve : curves)
			{
				std::string name = _pathResultFiles + SanitizeFilename(fileName) + "_" + std::to_string(i + 1) + ".txt";
				bool saved = Serializer::SaveParamCurveCartesian2D(*curve, title, t1, t2, numPoints, name);
				if (!saved) {
					return VisualizerResult::Failure("Failed to save curve " + std::to_string(i), name);
				}
				fileNames.push_back(name);
				i++;
			}
			return ExecuteVisualizer(_pathParametricCurve2DViz, fileNames);
		}

		static VisualizerResult VisualizeMultiParamCurve2D(std::vector<IRealToVectorFunction<2>*> curves,
																						std::vector<std::string> legend,
																						Real t1, Real t2,
																						int numPoints, std::string fileName)
		{
			// for each function, serialize data to file, with name generated from fileName
			std::vector<std::string> fileNames;
			int i = 0;
			for (auto& curve : curves)
			{
				std::string name = _pathResultFiles + SanitizeFilename(fileName) + "_" + std::to_string(i + 1) + ".txt";
				bool saved = Serializer::SaveParamCurveCartesian2D(*curve, legend[i], t1, t2, numPoints, name);
				if (!saved) {
					return VisualizerResult::Failure("Failed to save curve " + std::to_string(i), name);
				}
				fileNames.push_back(name);
				i++;
			}
			return ExecuteVisualizer(_pathParametricCurve2DViz, fileNames);
		}

		static VisualizerResult VisualizeMultiParamCurve2D(std::vector<std::string> fileNames)
		{
			std::vector<std::string> fullPaths;
			for (auto& name : fileNames)
			{
				std::string fullPath = _pathResultFiles + SanitizeFilename(name);
				if (!std::filesystem::exists(fullPath)) {
					return VisualizerResult::Failure("Data file not found: " + fullPath, fullPath);
				}
				fullPaths.push_back(fullPath);
			}
			return ExecuteVisualizer(_pathParametricCurve2DViz, fullPaths);
		}

		static VisualizerResult VisualizeParamCurve3D(const IRealToVectorFunction<3>& f, std::string title, 
																			Real t1, Real t2, int numPoints, 
																			std::string fileName)
		{
			std::string name = _pathResultFiles + SanitizeFilename(fileName);
			bool saved = Serializer::SaveParamCurveCartesian3D(f, title, t1, t2, numPoints, name);
			if (!saved) {
				return VisualizerResult::Failure("Failed to save parametric curve data", name);
			}
			return ExecuteVisualizer(_pathParametricCurve3DViz, {name});
		}

		static VisualizerResult VisualizeMultiParamCurve3D(std::vector<IRealToVectorFunction<3>*> curves, 
																					 std::string title,
																					 Real t1, Real t2, int numPoints, 
																					 std::string fileName)
		{
			// for each function, serialize data to file, with name generated from fileName
			std::vector<std::string> fileNames;
			int i = 0;
			for (auto& curve : curves)
			{
				std::string name = _pathResultFiles + SanitizeFilename(fileName) + "_" + std::to_string(i + 1) + ".txt";
				bool saved = Serializer::SaveParamCurveCartesian3D(*curve, title, t1, t2, numPoints, name);
				if (!saved) {
					return VisualizerResult::Failure("Failed to save curve " + std::to_string(i), name);
				}
				fileNames.push_back(name);
				i++;
			}
			return ExecuteVisualizer(_pathParametricCurve3DViz, fileNames);
		}

		static VisualizerResult VisualizeMultiParamCurve3D(std::vector<std::string> fileNames)
		{
			std::vector<std::string> fullPaths;
			for (auto& name : fileNames)
			{
				std::string fullPath = _pathResultFiles + SanitizeFilename(name);
				if (!std::filesystem::exists(fullPath)) {
					return VisualizerResult::Failure("Data file not found: " + fullPath, fullPath);
				}
				fullPaths.push_back(fullPath);
			}
			return ExecuteVisualizer(_pathParametricCurve3DViz, fullPaths);
		}

		// ODE Solution visualizations
		// Visualizing single variable of ODE system solution as Real function
		static VisualizerResult VisualizeODESysSolCompAsFunc(const ODESystemSolution& sol, int compInd,
																				 std::string title, std::string fileName)
		{
			std::string name = _pathResultFiles + SanitizeFilename(fileName);
			auto saveResult = Serializer::SaveODESolutionComponentAsFunc(sol, compInd, title, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(_pathRealFuncViz, {name});
		}

		// Visualizing ODE system solution as a multi-function (all variables)
		static VisualizerResult VisualizeODESysSolAsMultiFunc(const ODESystemSolution& sol,
																				std::string title, std::vector<std::string> legend, std::string fileName)
		{
			std::string name = _pathResultFiles + SanitizeFilename(fileName);
			auto saveResult = Serializer::SaveODESolutionAsMultiFunc(sol, title, legend, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(_pathRealFuncViz, {name});
		}

		// Visualizing two variables of ODE system solution as a parametric curve in 2D
		static VisualizerResult VisualizeODESysSolAsParamCurve2(const ODESystemSolution& sol,
																								int ind1, int ind2,
																								std::string title, std::string fileName)
		{
			std::string name = _pathResultFiles + SanitizeFilename(fileName);
			auto saveResult = Serializer::SaveODESolAsParametricCurve2D(sol, name, ind1, ind2, title);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(_pathParametricCurve2DViz, {name});
		}

		// Visualizing three variables of ODE system solution as a parametric curve in 3D
		static VisualizerResult VisualizeODESysSolAsParamCurve3(const ODESystemSolution& sol,
																								int ind1, int ind2, int ind3,
																								std::string title, std::string fileName)
		{
			std::string name = _pathResultFiles + SanitizeFilename(fileName);
			auto saveResult = Serializer::SaveODESolAsParametricCurve3D(sol, name, ind1, ind2, ind3, title);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(_pathParametricCurve3DViz, {name});
		}

		// Particle simulation visualizations
		static VisualizerResult VisualizeParticleSimulation2D(std::string fileName)
		{
			std::string name = _pathResultFiles + SanitizeFilename(fileName);
			if (!std::filesystem::exists(name)) {
				return VisualizerResult::Failure("Data file not found: " + name, name);
			}
			return ExecuteVisualizer(_pathParticle2DViz, {name});
		}

		static VisualizerResult VisualizeParticleSimulation3D(std::string fileName)
		{
			std::string name = _pathResultFiles + SanitizeFilename(fileName);
			if (!std::filesystem::exists(name)) {
				return VisualizerResult::Failure("Data file not found: " + name, name);
			}
			return ExecuteVisualizer(_pathParticle3DViz, {name});
		}
	};
}
#endif
