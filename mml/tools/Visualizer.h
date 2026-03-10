///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Visualizer.h                                                        ///
///  Description: Visualization tools for functions, curves, and data                 ///
///               Launches external visualizer apps for 2D/3D plotting                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_VIZUALIZER_H
#define MML_VIZUALIZER_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/InterpolatedFunction.h"
#include "base/ODESystemSolution.h"

#include "algorithms/FieldLineTracer.h"
#include "tools/Serializer.h"

#include <filesystem>
#include <algorithm>
#include <string>
#include <vector>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#else
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>
#endif


namespace MML {
	/// @brief Special timeout value indicating no timeout (wait indefinitely)
	static constexpr int VISUALIZER_NO_TIMEOUT = -1;
	
	/// @brief Default timeout for visualizer operations (30 seconds)
	static constexpr int VISUALIZER_DEFAULT_TIMEOUT_MS = 30000;

	// Result struct for visualizer operations - provides error information
	struct VisualizerResult {
		bool success;
		bool timedOut;  // True if operation timed out
		int exitCode;
		std::string errorMessage;
		std::string dataFilePath; // Path to the data file that was created

		static VisualizerResult Success(const std::string& dataPath = "") { 
			return VisualizerResult{true, false, 0, "", dataPath}; 
		}

		static VisualizerResult Failure(const std::string& message, const std::string& dataPath = "", int code = -1) {
			return VisualizerResult{false, false, code, message, dataPath};
		}
		
		static VisualizerResult Timeout(const std::string& message, const std::string& dataPath = "") {
			return VisualizerResult{false, true, -1, message, dataPath};
		}

		// Allow implicit conversion to bool for backward compatibility
		operator bool() const { return success; }
	};

	class Visualizer {
		// Result files path is static (doesn't depend on backend)
		static inline std::string _pathResultFiles{GetResultFilesPath()};

		// Dynamic path accessors - call getter each time to respect runtime backend changes
		static std::string pathRealFuncViz() { return GetRealFuncVisualizerPath(); }
		static std::string pathSurfaceViz() { return GetSurfaceVisualizerPath(); }
		static std::string pathParametricCurve3DViz() { return GetParametricCurve3DVisualizerPath(); }
		static std::string pathParametricCurve2DViz() { return GetParametricCurve2DVisualizerPath(); }
		static std::string pathVectorField2DViz() { return GetVectorField2DVisualizerPath(); }
		static std::string pathVectorField3DViz() { return GetVectorField3DVisualizerPath(); }
		static std::string pathParticle2DViz() { return GetParticle2DVisualizerPath(); }
		static std::string pathParticle3DViz() { return GetParticle3DVisualizerPath(); }
		static std::string pathParametricSurfaceViz() { return GetParametricSurfaceVisualizerPath(); }
		static std::string pathScalarFunc3DViz() { return GetScalarFunction3DVisualizerPath(); }
		static std::string pathRigidBodyViz() { return GetRigidBodyVisualizerPath(); }

		// Helper: resolve a file path for "FromFile" methods.
		// Tries the path as-is first (absolute or relative), then falls back to the results folder.
		// Returns the resolved path, or empty string if the file is not found anywhere.
		static std::string ResolveFilePath(const std::string& fileNameOrPath) {
			// 1. Try as-is (absolute path or relative to cwd)
			if (std::filesystem::exists(fileNameOrPath))
				return std::filesystem::weakly_canonical(fileNameOrPath).string();

			// 2. Try in results folder (legacy: bare filename)
			std::string inResults = MakeSafeOutputPath(fileNameOrPath);
			if (!inResults.empty() && std::filesystem::exists(inResults))
				return inResults;

			return {};  // not found
		}

		// Helper: resolve multiple file paths (for multi-file visualizers)
		static std::pair<bool, std::vector<std::string>> ResolveFilePaths(const std::vector<std::string>& fileNames) {
			std::vector<std::string> resolved;
			for (const auto& f : fileNames) {
				std::string path = ResolveFilePath(f);
				if (path.empty())
					return {false, {f}};  // store failed name in vector for error reporting
				resolved.push_back(path);
			}
			return {true, resolved};
		}

		// Helper: validate filename contains only safe basename characters (no path separators)
		static bool IsFilenameSafe(const std::string& filename) {
			for (char c : filename) {
				// Allow alphanumeric, underscore, hyphen, dot ONLY - no path separators
				if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_' && c != '-' && c != '.') {
					return false;
				}
			}
			return !filename.empty();
		}

		// Helper: sanitize filename to safe basename (no path components)
		// Extracts basename, replaces unsafe chars with underscore
		static std::string SanitizeFilename(const std::string& filename) {
			// Extract basename only - remove any directory components
			std::filesystem::path p(filename);
			std::string basename = p.filename().string();

			// If empty after extraction, use a default name
			if (basename.empty()) {
				basename = "output";
			}

			// Replace any remaining unsafe characters
			std::string result;
			result.reserve(basename.size());
			for (char c : basename) {
				if (std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-' || c == '.') {
					result += c;
				} else {
					result += '_';
				}
			}

			// Ensure result is not empty and doesn't start with dots (hidden files / traversal)
			while (!result.empty() && result[0] == '.') {
				result = result.substr(1);
			}
			if (result.empty()) {
				result = "output";
			}

			return result;
		}

		// Helper: build safe output path within results folder
		// Returns empty string if path escapes the results folder
		static std::string MakeSafeOutputPath(const std::string& fileName) {
			std::string safeBasename = SanitizeFilename(fileName);
			std::filesystem::path resultsDir(_pathResultFiles);
			std::filesystem::path outputPath = resultsDir / safeBasename;

			// Normalize the path (resolve . and ..)
			std::filesystem::path canonical;
			try {
				// Create parent directories if needed
				if (!std::filesystem::exists(resultsDir)) {
					std::filesystem::create_directories(resultsDir);
				}
				// Use weakly_canonical since output file may not exist yet
				canonical = std::filesystem::weakly_canonical(outputPath);
				std::filesystem::path canonicalResults = std::filesystem::weakly_canonical(resultsDir);

				// Verify containment: canonical path must start with results directory
				auto mismatch = std::mismatch(canonicalResults.begin(), canonicalResults.end(), canonical.begin(), canonical.end());
				if (mismatch.first != canonicalResults.end()) {
					// Path escapes the results folder - return empty (will cause failure)
					return "";
				}
			} catch (const std::filesystem::filesystem_error&) {
				// Path resolution failed - return empty
				return "";
			}

			return canonical.string();
		}

		// Overload for appending suffix (e.g., "_1.mml")
		static std::string MakeSafeOutputPath(const std::string& fileName, const std::string& suffix) {
			std::string safeBasename = SanitizeFilename(fileName) + suffix;
			std::filesystem::path resultsDir(_pathResultFiles);
			std::filesystem::path outputPath = resultsDir / safeBasename;

			std::filesystem::path canonical;
			try {
				if (!std::filesystem::exists(resultsDir)) {
					std::filesystem::create_directories(resultsDir);
				}
				canonical = std::filesystem::weakly_canonical(outputPath);
				std::filesystem::path canonicalResults = std::filesystem::weakly_canonical(resultsDir);

				auto mismatch = std::mismatch(canonicalResults.begin(), canonicalResults.end(), canonical.begin(), canonical.end());
				if (mismatch.first != canonicalResults.end()) {
					return "";
				}
			} catch (const std::filesystem::filesystem_error&) {
				return "";
			}

			return canonical.string();
		}

		// Helper: normalize path separators for the current platform
		static std::string NormalizePath(const std::string& path) {
			std::string result = path;
#ifdef _WIN32
			for (char& c : result) {
				if (c == '/')
					c = '\\';
			}
#else
			for (char& c : result) {
				if (c == '\\')
					c = '/';
			}
#endif
			return result;
		}

		// Safe process execution without shell - Windows implementation
#ifdef _WIN32
		/// @brief Execute visualizer on Windows with optional timeout
		/// @param executable Path to the visualizer executable
		/// @param args Command line arguments (file paths)
		/// @param timeout_ms Timeout in milliseconds. Use VISUALIZER_NO_TIMEOUT (-1) for infinite wait.
		/// @return VisualizerResult with success/failure/timeout status
		static VisualizerResult ExecuteVisualizerWindows(const std::string& executable, 
																										 const std::vector<std::string>& args,
																										 int timeout_ms = VISUALIZER_NO_TIMEOUT) {
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

			BOOL success = CreateProcessA(nullptr,                         // Application name (use command line)
																		cmdLineBuffer.data(),                // Command line
																		nullptr,                             // Process security attributes
																		nullptr,                             // Thread security attributes
																		FALSE,                               // Inherit handles
																		0,                                   // Creation flags
																		nullptr,                             // Environment
																		workDir.empty() ? nullptr : workDir.c_str(), // Working directory
																		&si,                                 // Startup info
																		&pi                                  // Process info
			);

			if (!success) {
				DWORD error = GetLastError();
				return VisualizerResult::Failure("Failed to start visualizer process. Error code: " + std::to_string(error),
																				 args.empty() ? "" : args[0], static_cast<int>(error));
			}

			// Convert timeout: -1 means INFINITE (0xFFFFFFFF)
			DWORD waitTimeout = (timeout_ms < 0) ? 0xFFFFFFFF : static_cast<DWORD>(timeout_ms);
			
			// Wait for GUI visualizer to close with specified timeout
			DWORD waitResult = WaitForSingleObject(pi.hProcess, waitTimeout);
			
			if (waitResult == WAIT_TIMEOUT) {
				// Process is still running - terminate it and report timeout
				TerminateProcess(pi.hProcess, 1);
				CloseHandle(pi.hProcess);
				CloseHandle(pi.hThread);
				return VisualizerResult::Timeout(
					"Visualizer process timed out after " + std::to_string(timeout_ms) + " ms",
					args.empty() ? "" : args[0]);
			}

			DWORD dwExitCode = 0;
			GetExitCodeProcess(pi.hProcess, &dwExitCode);

			CloseHandle(pi.hProcess);
			CloseHandle(pi.hThread);

			return VisualizerResult::Success(args.empty() ? "" : args[0]);
		}
#else
		// Safe process execution without shell - POSIX implementation
		/// @brief Execute visualizer on POSIX with optional timeout
		/// @param executable Path to the visualizer executable
		/// @param args Command line arguments (file paths)
		/// @param timeout_ms Timeout in milliseconds. Use VISUALIZER_NO_TIMEOUT (-1) for fire-and-forget.
		/// @return VisualizerResult with success/failure/timeout status
		static VisualizerResult ExecuteVisualizerPosix(const std::string& executable, 
																											   const std::vector<std::string>& args,
																											   int timeout_ms = VISUALIZER_NO_TIMEOUT) {
#ifdef __APPLE__
			// On macOS, detect .app bundles and launch via 'open' command.
			// Direct execv of a GUI binary doesn't register with the window server,
			// causing the app window to bounce in the Dock without receiving focus.
			std::string exePath = NormalizePath(executable);
			auto appPos = exePath.find(".app/Contents/MacOS/");
			if (appPos != std::string::npos) {
				// Extract the .app bundle path
				std::string appBundlePath = exePath.substr(0, appPos + 4); // includes ".app"
				
				// Build argv: open -W -a <bundle> --args <file1> <file2> ...
				// -W = wait for the app to quit before returning
				// -a = specify application by path
				std::vector<std::string> openArgs;
				openArgs.push_back("/usr/bin/open");
				openArgs.push_back("-W");
				openArgs.push_back("-a");
				openArgs.push_back(appBundlePath);
				if (!args.empty()) {
					openArgs.push_back("--args");
					for (const auto& arg : args) {
						openArgs.push_back(NormalizePath(arg));
					}
				}
				
				pid_t pid = fork();
				if (pid < 0) {
					return VisualizerResult::Failure("Failed to fork process", args.empty() ? "" : args[0], -1);
				}
				
				if (pid == 0) {
					// Child process
					std::vector<char*> argv;
					for (auto& a : openArgs) {
						argv.push_back(const_cast<char*>(a.c_str()));
					}
					argv.push_back(nullptr);
					
					execv("/usr/bin/open", argv.data());
					_exit(127);
				}
				
				// Parent: wait for 'open -W' to return (it waits for the app to quit)
				int status;
				waitpid(pid, &status, 0);
				if (WIFEXITED(status) && WEXITSTATUS(status) == 127) {
					return VisualizerResult::Failure("Failed to launch app bundle via open", 
																				  args.empty() ? "" : args[0], 127);
				}
				return VisualizerResult::Success(args.empty() ? "" : args[0]);
			}
			// Fall through for non-.app executables (e.g., FLTK plain binaries)
#endif

			pid_t pid = fork();

			if (pid < 0) {
				return VisualizerResult::Failure("Failed to fork process", args.empty() ? "" : args[0], -1);
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

#ifdef __linux__
				// Set Qt environment variables for the child process (Linux only)
				// Qt visualizers need LD_LIBRARY_PATH, QT_PLUGIN_PATH, and
				// QT_QPA_PLATFORM_PLUGIN_PATH to find bundled Qt libs and plugins
				std::string exeDir = exePath.substr(0, exePath.find_last_of('/'));
				std::string libPath = exeDir + "/lib";
				std::string pluginPath = exeDir + "/plugins";
				std::string platformPath = exeDir + "/plugins/platforms";
				
				// Prepend to existing LD_LIBRARY_PATH if it exists
				const char* existingLdPath = getenv("LD_LIBRARY_PATH");
				std::string fullLdPath = existingLdPath 
					? libPath + ":" + std::string(existingLdPath) 
					: libPath;
				setenv("LD_LIBRARY_PATH", fullLdPath.c_str(), 1);
				setenv("QT_PLUGIN_PATH", pluginPath.c_str(), 1);
				setenv("QT_QPA_PLATFORM_PLUGIN_PATH", platformPath.c_str(), 1);
#endif

				execv(exePath.c_str(), argv.data());

				// If execv returns, it failed
				_exit(127);
			}

			// Parent process
			if (timeout_ms < 0) {
				// Infinite wait mode (matches Windows INFINITE behavior)
				// Wait for the visualizer to close before returning
				int status;
				pid_t result = waitpid(pid, &status, 0);  // Blocking wait
				if (result == pid && WIFEXITED(status)) {
					int exitCode = WEXITSTATUS(status);
					if (exitCode == 127) {
						return VisualizerResult::Failure("Failed to execute visualizer", 
																					  args.empty() ? "" : args[0], 127);
					}
				}
				return VisualizerResult::Success(args.empty() ? "" : args[0]);
			}
			
			// Wait with timeout using polling (portable approach)
			const int poll_interval_ms = 50;  // Check every 50ms
			int elapsed_ms = 0;
			
			while (elapsed_ms < timeout_ms) {
				int status;
				pid_t result = waitpid(pid, &status, WNOHANG);
				
				if (result == pid) {
					// Child has exited
					if (WIFEXITED(status)) {
						int exitCode = WEXITSTATUS(status);
						if (exitCode == 127) {
							return VisualizerResult::Failure("Failed to execute visualizer", 
																						  args.empty() ? "" : args[0], 127);
						}
						return VisualizerResult::Success(args.empty() ? "" : args[0]);
					}
					return VisualizerResult::Success(args.empty() ? "" : args[0]);
				}
				else if (result < 0) {
					// Error occurred
					return VisualizerResult::Failure("Error waiting for process", 
																		  args.empty() ? "" : args[0], -1);
				}
				
				// Child still running - sleep and try again
				usleep(poll_interval_ms * 1000);  // usleep takes microseconds
				elapsed_ms += poll_interval_ms;
			}
			
			// Timeout reached - kill the process
			kill(pid, SIGTERM);
			usleep(100000);  // Give it 100ms to terminate gracefully
			kill(pid, SIGKILL);  // Force kill if still running
			waitpid(pid, nullptr, 0);  // Reap the zombie
			
			return VisualizerResult::Timeout(
				"Visualizer process timed out after " + std::to_string(timeout_ms) + " ms",
				args.empty() ? "" : args[0]);
		}
#endif

		// Main execution dispatcher
		/// @brief Execute visualizer with optional timeout
		/// @param executable Path to the visualizer executable
		/// @param args Command line arguments (file paths)
		/// @param timeout_ms Timeout in milliseconds. Default is VISUALIZER_NO_TIMEOUT (-1) for 
		///                   backward compatibility (Windows: infinite wait, POSIX: fire-and-forget).
		///                   Use VISUALIZER_DEFAULT_TIMEOUT_MS for 30-second timeout.
		/// @return VisualizerResult with success/failure/timeout status
		static VisualizerResult ExecuteVisualizer(const std::string& executable, 
																						   const std::vector<std::string>& args,
																						   int timeout_ms = VISUALIZER_NO_TIMEOUT) {
			// Validate executable path
			if (executable.empty()) {
				return VisualizerResult::Failure("Visualizer executable path is empty");
			}

			// Validate all arguments (file paths) - check for shell injection characters
			// Allow path separators (/, \), drive letters (:), and safe filename chars
			for (const auto& arg : args) {
				for (char c : arg) {
					// Reject shell metacharacters that could enable command injection
					if (c == '|' || c == '&' || c == ';' || c == '`' || 
					    c == '$' || c == '(' || c == ')' || c == '<' || c == '>' ||
					    c == '\'' || c == '"' || c == '\n' || c == '\r') {
						return VisualizerResult::Failure("Invalid characters in file path: " + arg, arg);
					}
				}
			}

#ifdef _WIN32
			return ExecuteVisualizerWindows(executable, args, timeout_ms);
#else
			return ExecuteVisualizerPosix(executable, args, timeout_ms);
#endif
		}

	public:
		/// @brief Visualize a real function over an interval
		static VisualizerResult VisualizeRealFunction(const IRealFunction& f, std::string title, Real x1, Real x2, int numPoints,
																									std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult = Serializer::SaveRealFunc(f, title, x1, x2, numPoints, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(pathRealFuncViz(), {name});
		}
		/// @brief Visualize a real function at specified points
		static VisualizerResult VisualizeRealFunction(const IRealFunction& f, std::string title, Vector<Real> points, std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult = Serializer::SaveRealFunc(f, title, points, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(pathRealFuncViz(), {name});
		}
		/// @brief Visualize multiple real functions on the same plot
		static VisualizerResult VisualizeMultiRealFunction(std::vector<IRealFunction*> funcs, std::string title,
																											 std::vector<std::string> func_legend, Real x1, Real x2, int numPoints,
																											 std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult = Serializer::SaveRealMultiFunc(funcs, title, func_legend, x1, x2, numPoints, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(pathRealFuncViz(), {name});
		}
		/// @brief Visualize multiple linear interpolation functions on the same plot
		static VisualizerResult VisualizeMultiRealFunction(std::vector<LinearInterpRealFunc> funcs, std::string title,
																											 std::vector<std::string> func_legend, Real x1, Real x2, int numPoints,
																											 std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult = Serializer::SaveRealMultiFunc(funcs, title, func_legend, x1, x2, numPoints, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(pathRealFuncViz(), {name});
		}
		/// @brief Visualize multiple functions as separate plots
		static VisualizerResult VisualizeMultiRealFunctionSeparately(std::vector<LinearInterpRealFunc> funcs, std::string title,
																																 std::vector<std::string> func_legend, Real x1, Real x2, int numPoints,
																																 std::string fileName) {
			// for each function, serialize data to file, with name generated from fileName
			std::vector<std::string> fileNames;
			int i = 0;
			for (auto& func : funcs) {
				std::string suffix = "_" + std::to_string(i + 1) + ".mml";
				std::string name = MakeSafeOutputPath(fileName, suffix);
				if (name.empty()) {
					return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
				}

				Real xLow = x1 < funcs[i].MinX() ? funcs[i].MinX() : x1;
				Real xHigh = x2 > funcs[i].MaxX() ? funcs[i].MaxX() : x2;

				auto saveResult = Serializer::SaveRealFunc(funcs[i], func_legend[i], xLow, xHigh, numPoints, name);
				if (!saveResult.success) {
					return VisualizerResult::Failure("Failed to save function " + std::to_string(i) + ": " + saveResult.message, name);
				}
				fileNames.push_back(name);
				i++;
			}
			return ExecuteVisualizer(pathRealFuncViz(), fileNames);
		}
		/// @brief Visualize multiple polynomial interpolation functions on the same plot
		static VisualizerResult VisualizeMultiRealFunction(std::vector<PolynomInterpRealFunc> funcs, std::string title,
																											 std::vector<std::string> func_legend, Real x1, Real x2, int numPoints,
																											 std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult = Serializer::SaveRealMultiFunc(funcs, title, func_legend, x1, x2, numPoints, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(pathRealFuncViz(), {name});
		}
		/// @brief Visualize multiple spline interpolation functions on the same plot
		static VisualizerResult VisualizeMultiRealFunction(std::vector<SplineInterpRealFunc> funcs, std::string title,
																											 std::vector<std::string> func_legend, Real x1, Real x2, int numPoints,
																											 std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult = Serializer::SaveRealMultiFunc(funcs, title, func_legend, x1, x2, numPoints, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(pathRealFuncViz(), {name});
		}
		/// @brief Visualize a pre-saved real function data file
		/// @param fileNameOrPath Full path or filename (looked up in results folder as fallback)
		/// @return VisualizerResult indicating success or failure
		static VisualizerResult VisualizeRealFunctionFromFile(std::string fileNameOrPath) {
			std::string resolved = ResolveFilePath(fileNameOrPath);
			if (resolved.empty())
				return VisualizerResult::Failure("Data file not found: " + fileNameOrPath, fileNameOrPath);
			return ExecuteVisualizer(pathRealFuncViz(), {resolved});
		}

		/// @brief Visualize a 2D scalar function as a surface plot
		static VisualizerResult VisualizeScalarFunc2DCartesian(const IScalarFunction<2>& func, std::string title, Real x1, Real x2,
																													 int numPointsX, Real y1, Real y2, int numPointsY, std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult = Serializer::SaveScalarFunc2DCartesian(func, title, x1, x2, numPointsX, y1, y2, numPointsY, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(pathSurfaceViz(), {name});
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
		/// @return VisualizerResult with success/failure info
		///
		/// @example
		/// auto solution = poissonSolver.solve();
		/// Visualizer::VisualizeGridFunction2D(solution, "Temperature Distribution", "temp.mml", 100.0, 100.0);
		template<typename GridFunc>
		static VisualizerResult VisualizeGridFunction2D(const GridFunc& gridFunc, std::string title, std::string fileName, Real scaleXY = 1.0,
																										Real scaleValue = 1.0) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult = Serializer::SaveGridFunction2D(gridFunc, title, name, scaleXY, scaleValue);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save grid data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(pathSurfaceViz(), {name});
		}

		/// @brief Visualize a pre-saved 2D scalar function data file
		/// @param fileNameOrPath Full path or filename (looked up in results folder as fallback)
		/// @return VisualizerResult indicating success or failure
		static VisualizerResult VisualizeScalarFunc2DFromFile(std::string fileNameOrPath) {
			std::string resolved = ResolveFilePath(fileNameOrPath);
			if (resolved.empty())
				return VisualizerResult::Failure("Data file not found: " + fileNameOrPath, fileNameOrPath);
			return ExecuteVisualizer(pathSurfaceViz(), {resolved});
		}

		// visualizations of Vector fields
		static VisualizerResult VisualizeVectorField2DCartesian(const IVectorFunction<2>& func, std::string title, Real x1, Real x2,
																														int numPointsX, Real y1, Real y2, int numPointsY, std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult = Serializer::SaveVectorFunc2DCartesian(func, title, x1, x2, numPointsX, y1, y2, numPointsY, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(pathVectorField2DViz(), {name});
		}

		static VisualizerResult VisualizeVectorField3DCartesian(const IVectorFunction<3>& func, std::string title, Real x1, Real x2,
																														int numPointsX, Real y1, Real y2, int numPointsY, Real z1, Real z2,
																														int numPointsZ, std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult =
				Serializer::SaveVectorFunc3DCartesian(func, title, x1, x2, numPointsX, y1, y2, numPointsY, z1, z2, numPointsZ, name, 3.0);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(pathVectorField3DViz(), {name});
		}

		/// @brief Visualize a pre-saved 2D vector field data file
		/// @param fileNameOrPath Full path or filename (looked up in results folder as fallback)
		/// @return VisualizerResult indicating success or failure
		static VisualizerResult VisualizeVectorField2DFromFile(std::string fileNameOrPath) {
			std::string resolved = ResolveFilePath(fileNameOrPath);
			if (resolved.empty())
				return VisualizerResult::Failure("Data file not found: " + fileNameOrPath, fileNameOrPath);
			return ExecuteVisualizer(pathVectorField2DViz(), {resolved});
		}

		/// @brief Visualize a pre-saved 3D vector field data file
		/// @param fileNameOrPath Full path or filename (looked up in results folder as fallback)
		/// @return VisualizerResult indicating success or failure
		static VisualizerResult VisualizeVectorField3DFromFile(std::string fileNameOrPath) {
			std::string resolved = ResolveFilePath(fileNameOrPath);
			if (resolved.empty())
				return VisualizerResult::Failure("Data file not found: " + fileNameOrPath, fileNameOrPath);
			return ExecuteVisualizer(pathVectorField3DViz(), {resolved});
		}

		// visualizations of Parametric curves
		static VisualizerResult VisualizeParamCurve2D(const IRealToVectorFunction<2>& f, std::string title, Real t1, Real t2, int numPoints,
																									std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			bool saved = Serializer::SaveParamCurveCartesian2DResult(f, title, t1, t2, numPoints, name).success;
			if (!saved) {
				return VisualizerResult::Failure("Failed to save parametric curve data", name);
			}
			return ExecuteVisualizer(pathParametricCurve2DViz(), {name});
		}

		static VisualizerResult VisualizeMultiParamCurve2D(std::vector<IRealToVectorFunction<2>*> curves, std::string title, Real t1, Real t2,
																											 int numPoints, std::string fileName) {
			// for each function, serialize data to file, with name generated from fileName
			std::vector<std::string> fileNames;
			int i = 0;
			for (auto& curve : curves) {
				std::string suffix = "_" + std::to_string(i + 1) + ".mml";
				std::string name = MakeSafeOutputPath(fileName, suffix);
				if (name.empty()) {
					return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
				}
				bool saved = Serializer::SaveParamCurveCartesian2DResult(*curve, title, t1, t2, numPoints, name).success;
				if (!saved) {
					return VisualizerResult::Failure("Failed to save curve " + std::to_string(i), name);
				}
				fileNames.push_back(name);
				i++;
			}
			return ExecuteVisualizer(pathParametricCurve2DViz(), fileNames);
		}

		static VisualizerResult VisualizeMultiParamCurve2D(std::vector<IRealToVectorFunction<2>*> curves, std::vector<std::string> legend,
																											 Real t1, Real t2, int numPoints, std::string fileName) {
			// for each function, serialize data to file, with name generated from fileName
			std::vector<std::string> fileNames;
			int i = 0;
			for (auto& curve : curves) {
				std::string suffix = "_" + std::to_string(i + 1) + ".mml";
				std::string name = MakeSafeOutputPath(fileName, suffix);
				if (name.empty()) {
					return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
				}
				bool saved = Serializer::SaveParamCurveCartesian2DResult(*curve, legend[i], t1, t2, numPoints, name).success;
				if (!saved) {
					return VisualizerResult::Failure("Failed to save curve " + std::to_string(i), name);
				}
				fileNames.push_back(name);
				i++;
			}
			return ExecuteVisualizer(pathParametricCurve2DViz(), fileNames);
		}

		static VisualizerResult VisualizeMultiParamCurve2D(std::vector<std::string> fileNames) {
			auto [ok, resolved] = ResolveFilePaths(fileNames);
			if (!ok)
				return VisualizerResult::Failure("Data file not found: " + resolved[0], resolved[0]);
			return ExecuteVisualizer(pathParametricCurve2DViz(), resolved);
		}

		static VisualizerResult VisualizeParamCurve3D(const IRealToVectorFunction<3>& f, std::string title, Real t1, Real t2, int numPoints,
																									std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			bool saved = Serializer::SaveParamCurveCartesian3DResult(f, title, t1, t2, numPoints, name).success;
			if (!saved) {
				return VisualizerResult::Failure("Failed to save parametric curve data", name);
			}
			return ExecuteVisualizer(pathParametricCurve3DViz(), {name});
		}

		static VisualizerResult VisualizeMultiParamCurve3D(std::vector<IRealToVectorFunction<3>*> curves, std::string title, Real t1, Real t2,
																											 int numPoints, std::string fileName) {
			// for each function, serialize data to file, with name generated from fileName
			std::vector<std::string> fileNames;
			int i = 0;
			for (auto& curve : curves) {
				std::string suffix = "_" + std::to_string(i + 1) + ".mml";
				std::string name = MakeSafeOutputPath(fileName, suffix);
				if (name.empty()) {
					return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
				}
				bool saved = Serializer::SaveParamCurveCartesian3DResult(*curve, title, t1, t2, numPoints, name).success;
				if (!saved) {
					return VisualizerResult::Failure("Failed to save curve " + std::to_string(i), name);
				}
				fileNames.push_back(name);
				i++;
			}
			return ExecuteVisualizer(pathParametricCurve3DViz(), fileNames);
		}

		static VisualizerResult VisualizeMultiParamCurve3D(std::vector<std::string> fileNames) {
			auto [ok, resolved] = ResolveFilePaths(fileNames);
			if (!ok)
				return VisualizerResult::Failure("Data file not found: " + resolved[0], resolved[0]);
			return ExecuteVisualizer(pathParametricCurve3DViz(), resolved);
		}

		// ODE Solution visualizations
		// Visualizing single variable of ODE system solution as Real function
		static VisualizerResult VisualizeODESysSolCompAsFunc(const ODESystemSolution& sol, int compInd, std::string title,
																												 std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult = Serializer::SaveODESolutionComponentAsFunc(sol, compInd, title, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(pathRealFuncViz(), {name});
		}

		// Visualizing ODE system solution as a multi-function (all variables)
		static VisualizerResult VisualizeODESysSolAsMultiFunc(const ODESystemSolution& sol, std::string title, std::vector<std::string> legend,
																													std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult = Serializer::SaveODESolutionAsMultiFunc(sol, title, legend, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(pathRealFuncViz(), {name});
		}

		// Visualizing two variables of ODE system solution as a parametric curve in 2D
		static VisualizerResult VisualizeODESysSolAsParamCurve2(const ODESystemSolution& sol, int ind1, int ind2, std::string title,
																														std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult = Serializer::SaveODESolAsParametricCurve2D(sol, name, ind1, ind2, title);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(pathParametricCurve2DViz(), {name});
		}

		// Visualizing three variables of ODE system solution as a parametric curve in 3D
		static VisualizerResult VisualizeODESysSolAsParamCurve3(const ODESystemSolution& sol, int ind1, int ind2, int ind3, std::string title,
																														std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult = Serializer::SaveODESolAsParametricCurve3D(sol, name, ind1, ind2, ind3, title);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(pathParametricCurve3DViz(), {name});
		}

		// Particle simulation visualizations
		static VisualizerResult VisualizeParticleSimulation2D(std::string fileNameOrPath) {
			std::string resolved = ResolveFilePath(fileNameOrPath);
			if (resolved.empty())
				return VisualizerResult::Failure("Data file not found: " + fileNameOrPath, fileNameOrPath);
			return ExecuteVisualizer(pathParticle2DViz(), {resolved});
		}

		static VisualizerResult VisualizeParticleSimulation3D(std::string fileNameOrPath) {
			std::string resolved = ResolveFilePath(fileNameOrPath);
			if (resolved.empty())
				return VisualizerResult::Failure("Data file not found: " + fileNameOrPath, fileNameOrPath);
			return ExecuteVisualizer(pathParticle3DViz(), {resolved});
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
		/// @return VisualizerResult indicating success or failure
		static VisualizerResult VisualizeParametricSurface(const IParametricSurface<3>& surface, std::string title, Real u1, Real u2,
																											 int numPointsU, Real w1, Real w2, int numPointsW, std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult = Serializer::SaveParametricSurface(surface, title, u1, u2, numPointsU, w1, w2, numPointsW, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save parametric surface data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(pathParametricSurfaceViz(), {name});
		}

		/// @brief Visualize a 3D parametric surface using its intrinsic bounds
		/// @details For surfaces with rectangular domains (IParametricSurfaceRect), uses
		/// the surface's getMinU/getMaxU/getMinW/getMaxW bounds automatically.
		/// @param surface The parametric surface with rectangular domain
		/// @param title Title for the visualization
		/// @param numPointsU Number of sample points along u
		/// @param numPointsW Number of sample points along w
		/// @param fileName Output filename (without path)
		/// @return VisualizerResult indicating success or failure
		static VisualizerResult VisualizeParametricSurface(const IParametricSurfaceRect<3>& surface, std::string title, int numPointsU,
																											 int numPointsW, std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult = Serializer::SaveParametricSurface(surface, title, numPointsU, numPointsW, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save parametric surface data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(pathParametricSurfaceViz(), {name});
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
		/// @return VisualizerResult indicating success or failure
		static VisualizerResult VisualizeParametricSurface(const IVectorFunctionNM<2, 3>& f, std::string title, Real u1, Real u2,
																											 int numPointsU, Real w1, Real w2, int numPointsW, std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult = Serializer::SaveParametricSurface(f, title, u1, u2, numPointsU, w1, w2, numPointsW, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save parametric surface data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(pathParametricSurfaceViz(), {name});
		}

		/// @brief Visualize a pre-saved parametric surface data file
		/// @param fileName Name of the data file (in results folder)
		/// @return VisualizerResult indicating success or failure
		static VisualizerResult VisualizeParametricSurfaceFromFile(std::string fileNameOrPath) {
			std::string resolved = ResolveFilePath(fileNameOrPath);
			if (resolved.empty())
				return VisualizerResult::Failure("Data file not found: " + fileNameOrPath, fileNameOrPath);
			return ExecuteVisualizer(pathParametricSurfaceViz(), {resolved});
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
		/// @return VisualizerResult indicating success or failure
		/// @note Large grids (e.g., 50x50x50 = 125,000 points) can create large data files
		static VisualizerResult VisualizeScalarFunc3DCartesian(const IScalarFunction<3>& func, std::string title, Real x1, Real x2,
																													 int numPointsX, Real y1, Real y2, int numPointsY, Real z1, Real z2,
																													 int numPointsZ, std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult =
				Serializer::SaveScalarFunc3DCartesian(func, title, x1, x2, numPointsX, y1, y2, numPointsY, z1, z2, numPointsZ, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save 3D scalar function data: " + saveResult.message, name);
			}
			return ExecuteVisualizer(pathScalarFunc3DViz(), {name});
		}

		/// @brief Visualize a pre-saved 3D scalar function data file
		/// @param fileName Name of the data file (in results folder)
		/// @return VisualizerResult indicating success or failure
		static VisualizerResult VisualizeScalarFunc3DFromFile(std::string fileNameOrPath) {
			std::string resolved = ResolveFilePath(fileNameOrPath);
			if (resolved.empty())
				return VisualizerResult::Failure("Data file not found: " + fileNameOrPath, fileNameOrPath);
			return ExecuteVisualizer(pathScalarFunc3DViz(), {resolved});
		}

		// Rigid body trajectory visualizations
		/// @brief Visualize a single rigid body trajectory file
		/// @details Launches the rigid body visualizer with the specified data file.
		/// The file should be in RIGID_BODY_TRAJECTORY_3D format (single body trajectory).
		/// Use this with files like rigid_body_box1.mml, rigid_body_box2.mml.
		/// @param fileName Name or path of the data file (.mml format)
		/// @return VisualizerResult indicating success or failure
		static VisualizerResult VisualizeRigidBodyTrajectory(std::string fileName) {
			// First try as-is (for full paths)
			if (std::filesystem::exists(fileName)) {
				return ExecuteVisualizer(pathRigidBodyViz(), {fileName});
			}
			// Try in results folder
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			if (!std::filesystem::exists(name)) {
				return VisualizerResult::Failure("Data file not found: " + name, name);
			}
			return ExecuteVisualizer(pathRigidBodyViz(), {name});
		}

		/// @brief Visualize multiple rigid body trajectory files together
		/// @details Useful for comparing trajectories side by side.
		/// Each file should be in RIGID_BODY_TRAJECTORY_3D format.
		/// @param fileNames Vector of file paths to visualize together
		/// @return VisualizerResult indicating success or failure
		static VisualizerResult VisualizeRigidBodyTrajectory(std::vector<std::string> fileNames) {
			std::vector<std::string> validPaths;
			for (const auto& fileName : fileNames) {
				// First try as-is (for full paths)
				if (std::filesystem::exists(fileName)) {
					validPaths.push_back(fileName);
					continue;
				}
				// Try in results folder
				std::string name = MakeSafeOutputPath(fileName);
				if (name.empty()) {
					return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
				}
				if (!std::filesystem::exists(name)) {
					return VisualizerResult::Failure("Data file not found: " + name, name);
				}
				validPaths.push_back(name);
			}
			return ExecuteVisualizer(pathRigidBodyViz(), validPaths);
		}

		//===================================================================================
		// Field Lines Visualization
		//===================================================================================

		/// @brief Visualize 2D field lines
		/// @details Saves field lines to file and launches the field lines visualizer.
		/// Currently outputs the file for external visualization (visualizer app coming soon).
		/// @param lines Vector of 2D field lines to visualize
		/// @param title Descriptive title for the visualization
		/// @param fileName Output file name (will be placed in results folder)
		/// @return VisualizerResult with success/failure and file path
		static VisualizerResult VisualizeFieldLines2D(const std::vector<FieldLine2D>& lines,
		                                               std::string title,
		                                               std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult = Serializer::SaveFieldLines2D(lines, title, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save field lines: " + saveResult.message, name);
			}
			// TODO: Launch MML_FieldLines2D_Visualizer when available
			// For now, just return success with the file path
			return VisualizerResult::Success(name);
		}

		/// @brief Visualize 3D field lines
		/// @details Saves field lines to file and launches the field lines visualizer.
		/// Currently outputs the file for external visualization (visualizer app coming soon).
		/// @param lines Vector of 3D field lines to visualize
		/// @param title Descriptive title for the visualization
		/// @param fileName Output file name (will be placed in results folder)
		/// @return VisualizerResult with success/failure and file path
		static VisualizerResult VisualizeFieldLines3D(const std::vector<FieldLine3D>& lines,
		                                               std::string title,
		                                               std::string fileName) {
			std::string name = MakeSafeOutputPath(fileName);
			if (name.empty()) {
				return VisualizerResult::Failure("Invalid filename: path escapes results folder", fileName);
			}
			auto saveResult = Serializer::SaveFieldLines3D(lines, title, name);
			if (!saveResult.success) {
				return VisualizerResult::Failure("Failed to save field lines: " + saveResult.message, name);
			}
			// TODO: Launch MML_FieldLines3D_Visualizer when available
			// For now, just return success with the file path
			return VisualizerResult::Success(name);
		}
	};
} // namespace MML
#endif
