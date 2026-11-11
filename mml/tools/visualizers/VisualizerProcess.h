///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///  File:        tools/visualizers/VisualizerProcess.h
///  Description: Shell-free platform process execution and checked outcomes
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                         ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_VISUALIZER_PROCESS_H
#define MML_VISUALIZER_PROCESS_H

#include <mml/tools/IOResult.h>

#include <cerrno>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>
#endif

namespace MML {
	/// @brief Special timeout value indicating no timeout (wait indefinitely)
	inline constexpr int VISUALIZER_NO_TIMEOUT = -1;

	/// @brief Default timeout for visualizer operations (30 seconds)
	inline constexpr int VISUALIZER_DEFAULT_TIMEOUT_MS = 30000;

}

namespace MML::detail::visualizer {
#ifdef _WIN32
	inline std::string QuoteWindowsArgument(const std::string& arg) {
			std::string quoted = "\"";
			std::size_t backslashes = 0;
			for (char c : arg) {
				if (c == '\\') {
					++backslashes;
					continue;
				}
				quoted.append(c == '"' ? 2 * backslashes + 1 : backslashes, '\\');
				quoted += c;
				backslashes = 0;
			}
			quoted.append(2 * backslashes, '\\');
			return quoted + '"';
		}
#else
	inline IOResult PosixExitResult(int status, const std::vector<std::string>& args);
#endif

	inline std::string NormalizePath(const std::string& path) {
		std::string result = path;
#ifdef _WIN32
		for (char& c : result) {
			if (c == '/')
				c = '\\';
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
	/// @return IOResult with success/failure/timeout status
	inline IOResult ExecuteVisualizerWindows(const std::string& executable, const std::vector<std::string>& args,
																					 int timeout_ms = VISUALIZER_NO_TIMEOUT) {
		// Build command line: "exe" "arg1" "arg2" ...
		const std::string exePath = NormalizePath(executable);
		std::string cmdLine = QuoteWindowsArgument(exePath);
		for (const auto& arg : args) {
			cmdLine += " " + QuoteWindowsArgument(arg);
		}

		// Get working directory (directory containing the executable)
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

		BOOL success = CreateProcessA(exePath.c_str(),                  // Explicit executable, never a shell
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
			return IOResult::Failure("Failed to start visualizer process. Error code: " + std::to_string(error), args.empty() ? "" : args[0],
															 static_cast<int>(error));
		}

		// Convert timeout: -1 means INFINITE (0xFFFFFFFF)
		DWORD waitTimeout = (timeout_ms < 0) ? std::numeric_limits<DWORD>::max() : static_cast<DWORD>(timeout_ms);

		// Wait for GUI visualizer to close with specified timeout
		DWORD waitResult = WaitForSingleObject(pi.hProcess, waitTimeout);

		if (waitResult == WAIT_TIMEOUT) {
			// Process is still running - terminate it and report timeout
			if (!TerminateProcess(pi.hProcess, 1)) {
				const DWORD error = GetLastError();
				CloseHandle(pi.hProcess);
				CloseHandle(pi.hThread);
				return IOResult::Failure("Failed to terminate timed-out visualizer. Error code: " + std::to_string(error),
					args.empty() ? "" : args[0], static_cast<int>(error));
			}
			const DWORD terminationWait = WaitForSingleObject(pi.hProcess, std::numeric_limits<DWORD>::max());
			const DWORD terminationError = terminationWait == WAIT_FAILED ? GetLastError() : 0;
			CloseHandle(pi.hProcess);
			CloseHandle(pi.hThread);
			if (terminationWait != WAIT_OBJECT_0)
				return IOResult::Failure("Failed to wait for terminated visualizer. Error code: " + std::to_string(terminationError),
					args.empty() ? "" : args[0], static_cast<int>(terminationError));
			return IOResult::Timeout("Visualizer process timed out after " + std::to_string(timeout_ms) + " ms", args.empty() ? "" : args[0]);
		}

		if (waitResult != WAIT_OBJECT_0) {
			const DWORD error = GetLastError();
			CloseHandle(pi.hProcess);
			CloseHandle(pi.hThread);
			return IOResult::Failure("Failed to wait for visualizer process. Error code: " + std::to_string(error),
				args.empty() ? "" : args[0], static_cast<int>(error));
		}

		DWORD dwExitCode = 0;
		const BOOL gotExitCode = GetExitCodeProcess(pi.hProcess, &dwExitCode);
		const DWORD exitCodeError = gotExitCode ? 0 : GetLastError();

		CloseHandle(pi.hProcess);
		CloseHandle(pi.hThread);

		if (!gotExitCode)
			return IOResult::Failure("Failed to read visualizer exit code. Error code: " + std::to_string(exitCodeError),
				args.empty() ? "" : args[0], static_cast<int>(exitCodeError));
		if (dwExitCode != 0)
			return IOResult::Failure("Visualizer exited with code " + std::to_string(dwExitCode),
				args.empty() ? "" : args[0], static_cast<int>(dwExitCode));
		return IOResult::Success(args.empty() ? "" : args[0]);
	}
#else
	// Safe process execution without shell - POSIX implementation
	/// @brief Execute visualizer on POSIX with optional timeout
	/// @param executable Path to the visualizer executable
	/// @param args Command line arguments (file paths)
	/// @param timeout_ms Timeout in milliseconds. Use VISUALIZER_NO_TIMEOUT (-1) for infinite wait.
	/// @return IOResult with success/failure/timeout status
	inline IOResult ExecuteVisualizerPosix(const std::string& executable, const std::vector<std::string>& args,
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
					openArgs.push_back(arg);
				}
			}

			pid_t pid = fork();
			if (pid < 0) {
				return IOResult::Failure("Failed to fork process", args.empty() ? "" : args[0], -1);
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
			int status = 0;
			pid_t result;
			do { result = waitpid(pid, &status, 0); } while (result < 0 && errno == EINTR);
			if (result < 0)
				return IOResult::Failure("Failed to wait for app bundle: " + std::string(std::strerror(errno)),
					args.empty() ? "" : args[0], errno);
			return PosixExitResult(status, args);
		}
		// Fall through for non-.app executables (e.g., FLTK plain binaries)
#endif

		pid_t pid = fork();

		if (pid < 0) {
			return IOResult::Failure("Failed to fork process", args.empty() ? "" : args[0], -1);
		}

		if (pid == 0) {
			// Child process
			std::vector<char*> argv;
			std::string exePath = NormalizePath(executable);
			argv.push_back(const_cast<char*>(exePath.c_str()));

			for (const auto& arg : args) {
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
			int status = 0;
			pid_t result;
			do { result = waitpid(pid, &status, 0); } while (result < 0 && errno == EINTR);
			if (result < 0) {
				return IOResult::Failure("waitpid failed: " + std::string(std::strerror(errno)), args.empty() ? "" : args[0], errno);
			}
			return PosixExitResult(status, args);
		}

		// Wait with timeout using polling (portable approach)
		const int poll_interval_ms = 50;  // Check every 50ms
		int elapsed_ms = 0;

		while (elapsed_ms < timeout_ms) {
			int status;
			pid_t result = waitpid(pid, &status, WNOHANG);

			if (result == pid) {
				return PosixExitResult(status, args);
			}
			else if (result < 0) {
				if (errno == EINTR)
					continue;
				// Error occurred
				return IOResult::Failure("Error waiting for process: " + std::string(std::strerror(errno)),
					args.empty() ? "" : args[0], errno);
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

		return IOResult::Timeout("Visualizer process timed out after " + std::to_string(timeout_ms) + " ms", args.empty() ? "" : args[0]);
	}

	inline IOResult PosixExitResult(int status, const std::vector<std::string>& args) {
		const std::string filename = args.empty() ? "" : args[0];
		if (WIFEXITED(status)) {
			const int exitCode = WEXITSTATUS(status);
			if (exitCode == 0)
				return IOResult::Success(filename);
			return IOResult::Failure("Visualizer exited with code " + std::to_string(exitCode), filename, exitCode);
		}
		if (WIFSIGNALED(status))
			return IOResult::Failure("Visualizer terminated by signal " + std::to_string(WTERMSIG(status)),
				filename, -WTERMSIG(status));
		return IOResult::Failure("Visualizer did not exit normally", filename);
	}
#endif

	// Main execution dispatcher
	/// @brief Execute visualizer with optional timeout
	/// @param executable Path to the visualizer executable
	/// @param args Command line arguments (file paths)
	/// @param timeout_ms Timeout in milliseconds. Default is VISUALIZER_NO_TIMEOUT (-1) for
	///                   an infinite wait on every platform.
	///                   Use VISUALIZER_DEFAULT_TIMEOUT_MS for 30-second timeout.
	/// @return IOResult with success/failure/timeout status
	inline IOResult ExecuteVisualizer(const std::string& executable, const std::vector<std::string>& args,
																		int timeout_ms = VISUALIZER_NO_TIMEOUT) {
		// Validate executable path
		if (executable.empty()) {
			return IOResult::Failure("Visualizer executable path is empty");
		}

		if (executable.find('\0') != std::string::npos)
			return IOResult::Failure("Visualizer executable path contains an embedded null", executable, -1, IOError::INVALID_PARAMETERS);
		if (timeout_ms < VISUALIZER_NO_TIMEOUT)
			return IOResult::Failure("Visualizer timeout must be -1 or nonnegative", executable, -1, IOError::INVALID_PARAMETERS);

		for (const auto& arg : args) {
			if (arg.find('\0') != std::string::npos)
				return IOResult::Failure("Visualizer argument contains an embedded null", arg, -1, IOError::INVALID_PARAMETERS);
		}
		if (!std::filesystem::is_regular_file(executable))
			return IOResult::Failure("Visualizer executable not found: " + executable, args.empty() ? "" : args[0],
				-1, IOError::FILE_NOT_OPENED);

#ifdef _WIN32
		return ExecuteVisualizerWindows(executable, args, timeout_ms);
#else
		return ExecuteVisualizerPosix(executable, args, timeout_ms);
#endif
	}
}

#endif // MML_VISUALIZER_PROCESS_H
