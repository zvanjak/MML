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
#include <algorithm>
#include <chrono>
#include <thread>
#include <spawn.h>
#include <unistd.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>
#ifdef __APPLE__
#include <crt_externs.h>
#else
extern char** environ;
#endif
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
	inline IOResult PosixExitResult(int status, const std::vector<std::string>& args, bool appBundle = false);
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

		// Loader dialogs can block CreateProcess itself, before the process wait timeout applies.
		DWORD previousErrorMode = 0;
		if (!SetThreadErrorMode(GetThreadErrorMode() | SEM_FAILCRITICALERRORS, &previousErrorMode)) {
			const DWORD error = GetLastError();
			return IOResult::Failure("Failed to suppress visualizer loader dialogs. Error code: " + std::to_string(error),
				args.empty() ? "" : args[0], static_cast<int>(error));
		}
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
		const DWORD startError = success ? ERROR_SUCCESS : GetLastError();
		const BOOL restoredErrorMode = SetThreadErrorMode(previousErrorMode, nullptr);
		const DWORD restoreError = restoredErrorMode ? ERROR_SUCCESS : GetLastError();
		const auto finish = [restoreError](IOResult result) {
			if (restoreError != ERROR_SUCCESS) {
				const auto message = "Failed to restore visualizer caller error mode. Error code: " + std::to_string(restoreError);
				if (result.success)
					return IOResult::Failure(message, result.filename, static_cast<int>(restoreError));
				result.message += "; " + message;
			}
			return result;
		};

		if (!success) {
			return finish(IOResult::Failure("Failed to start visualizer process. Error code: " + std::to_string(startError),
				args.empty() ? "" : args[0], static_cast<int>(startError)));
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
				return finish(IOResult::Failure("Failed to terminate timed-out visualizer. Error code: " + std::to_string(error),
					args.empty() ? "" : args[0], static_cast<int>(error)));
			}
			const DWORD terminationWait = WaitForSingleObject(pi.hProcess, std::numeric_limits<DWORD>::max());
			const DWORD terminationError = terminationWait == WAIT_FAILED ? GetLastError() : 0;
			CloseHandle(pi.hProcess);
			CloseHandle(pi.hThread);
			if (terminationWait != WAIT_OBJECT_0)
				return finish(IOResult::Failure("Failed to wait for terminated visualizer. Error code: " + std::to_string(terminationError),
					args.empty() ? "" : args[0], static_cast<int>(terminationError)));
			return finish(IOResult::Timeout("Visualizer process timed out after " + std::to_string(timeout_ms) + " ms", args.empty() ? "" : args[0]));
		}

		if (waitResult != WAIT_OBJECT_0) {
			const DWORD error = GetLastError();
			CloseHandle(pi.hProcess);
			CloseHandle(pi.hThread);
			return finish(IOResult::Failure("Failed to wait for visualizer process. Error code: " + std::to_string(error),
				args.empty() ? "" : args[0], static_cast<int>(error)));
		}

		DWORD dwExitCode = 0;
		const BOOL gotExitCode = GetExitCodeProcess(pi.hProcess, &dwExitCode);
		const DWORD exitCodeError = gotExitCode ? 0 : GetLastError();

		CloseHandle(pi.hProcess);
		CloseHandle(pi.hThread);

		if (!gotExitCode)
			return finish(IOResult::Failure("Failed to read visualizer exit code. Error code: " + std::to_string(exitCodeError),
				args.empty() ? "" : args[0], static_cast<int>(exitCodeError)));
		if (dwExitCode != 0)
			return finish(IOResult::Failure("Visualizer exited with code " + std::to_string(dwExitCode),
				args.empty() ? "" : args[0], static_cast<int>(dwExitCode)));
		return finish(IOResult::Success(args.empty() ? "" : args[0]));
	}
#else
	// Safe process execution without shell - POSIX implementation
	inline std::vector<std::string> AppBundleCommand(const std::string& executable, const std::vector<std::string>& args) {
		const auto appPos = executable.find(".app/Contents/MacOS/");
		if (appPos == std::string::npos)
			return {};
		std::vector<std::string> command{"/usr/bin/open", "-W", "-a", executable.substr(0, appPos + 4)};
		if (!args.empty()) {
			command.push_back("--args");
			command.insert(command.end(), args.begin(), args.end());
		}
		return command;
	}

	inline std::vector<std::string> PosixChildEnvironment(const std::string& executable) {
		std::vector<std::string> environment;
#ifdef __APPLE__
		char** inherited = *_NSGetEnviron();
#else
		char** inherited = ::environ;
#endif
		if (inherited) {
			for (char** entry = inherited; *entry; ++entry)
				environment.emplace_back(*entry);
		}
#ifdef __linux__
		auto directory = std::filesystem::path(executable).parent_path().string();
		if (directory.empty())
			directory = ".";
		std::string libraryPath = directory + "/lib";
		for (const auto& entry : environment) {
			if (entry.compare(0, 16, "LD_LIBRARY_PATH=") == 0) {
				libraryPath += ":" + entry.substr(16);
				break;
			}
		}
		const auto replace = [&environment](const std::string& key, const std::string& value) {
			const auto prefix = key + "=";
			environment.erase(std::remove_if(environment.begin(), environment.end(), [&prefix](const std::string& entry) {
				return entry.compare(0, prefix.size(), prefix) == 0;
			}), environment.end());
			environment.push_back(prefix + value);
		};
		replace("LD_LIBRARY_PATH", libraryPath);
		replace("QT_PLUGIN_PATH", directory + "/plugins");
		replace("QT_QPA_PLATFORM_PLUGIN_PATH", directory + "/plugins/platforms");
#else
		(void)executable;
#endif
		return environment;
	}

	inline pid_t WaitPosixChildUntil(pid_t pid, int& status, std::chrono::steady_clock::time_point deadline) {
		for (;;) {
			const pid_t result = waitpid(pid, &status, WNOHANG);
			if (result == pid || (result < 0 && errno != EINTR))
				return result;
			const auto now = std::chrono::steady_clock::now();
			if (now >= deadline)
				return 0;
			if (result == 0)
				std::this_thread::sleep_until(std::min(deadline, now + std::chrono::milliseconds(10)));
		}
	}

	inline IOResult ExecutePosixCommand(const std::vector<std::string>& command, const std::vector<std::string>& args,
																	 int timeout_ms, bool appBundle = false) {
		const std::string filename = args.empty() ? "" : args[0];
		const std::string subject = appBundle ? "App-bundle launcher (open -W)" : "Visualizer process";
		const auto failure = [&filename, &subject](const std::string& operation, int error) {
			return IOResult::Failure(subject + ": " + operation + ": " + std::string(std::strerror(error)), filename, error);
		};
		if (command.empty())
			return IOResult::Failure(subject + ": empty command", filename, -1, IOError::INVALID_PARAMETERS);

		// Prepare all storage in the parent; no application code runs between fork and exec.
		std::vector<char*> argv;
		for (const auto& argument : command)
			argv.push_back(const_cast<char*>(argument.c_str()));
		argv.push_back(nullptr);
		auto environment = PosixChildEnvironment(command[0]);
		std::vector<char*> envp;
		for (auto& entry : environment)
			envp.push_back(entry.data());
		envp.push_back(nullptr);
		const auto deadline = std::chrono::steady_clock::now() + std::chrono::milliseconds(std::max(timeout_ms, 0));
		pid_t pid = 0;
		const int spawnError = posix_spawn(&pid, command[0].c_str(), nullptr, nullptr, argv.data(), envp.data());
		if (spawnError != 0)
			return failure("posix_spawn failed", spawnError);

		int status = 0;
		pid_t result;
		if (timeout_ms < 0) {
			do { result = waitpid(pid, &status, 0); } while (result < 0 && errno == EINTR);
		} else {
			result = WaitPosixChildUntil(pid, status, deadline);
		}
		if (result == pid)
			return PosixExitResult(status, args, appBundle);
		if (result < 0)
			return failure("waitpid failed; child ownership or status unavailable", errno);

		int terminationError = 0;
		if (kill(pid, SIGTERM) != 0 && errno != ESRCH)
			terminationError = errno;
		result = WaitPosixChildUntil(pid, status, std::chrono::steady_clock::now() + std::chrono::milliseconds(100));
		if (result < 0)
			return failure("waitpid during timeout cleanup failed", errno);
		if (result == 0) {
			if (kill(pid, SIGKILL) != 0 && errno != ESRCH) {
				const int killError = errno;
				// Do not block indefinitely if the OS refused forced termination.
				do { result = waitpid(pid, &status, WNOHANG); } while (result < 0 && errno == EINTR);
				if (result < 0) {
					const int reapError = errno;
					return failure("SIGKILL failed (" + std::string(std::strerror(killError)) + "); final waitpid failed", reapError);
				}
				return failure(result == pid ? "SIGKILL failed; child subsequently reaped" :
					"SIGKILL failed; child could not be reaped and may still be running", killError);
			}
			do { result = waitpid(pid, &status, 0); } while (result < 0 && errno == EINTR);
			if (result < 0)
				return failure("final waitpid during timeout cleanup failed", errno);
		}
		if (terminationError != 0)
			return failure("SIGTERM failed; child reaped after timeout", terminationError);
		return IOResult::Timeout(subject + " timed out after " + std::to_string(timeout_ms) + " ms" +
			(appBundle ? "; launcher terminated and reaped, but the GUI may still be running; GUI exit code is unavailable" : ""),
			filename);
	}

	/// @brief Execute a POSIX visualizer, timing out and reaping only the owned child.
	inline IOResult ExecuteVisualizerPosix(const std::string& executable, const std::vector<std::string>& args,
																			 int timeout_ms = VISUALIZER_NO_TIMEOUT) {
		std::vector<std::string> command;
		bool appBundle = false;
#ifdef __APPLE__
		// Preserve LaunchServices activation; open's status is not the GUI's exit code.
		command = AppBundleCommand(executable, args);
		appBundle = !command.empty();
#endif
		if (!appBundle) {
			command.push_back(executable);
			command.insert(command.end(), args.begin(), args.end());
		}
		return ExecutePosixCommand(command, args, timeout_ms, appBundle);
	}

	inline IOResult PosixExitResult(int status, const std::vector<std::string>& args, bool appBundle) {
		const std::string filename = args.empty() ? "" : args[0];
		const std::string subject = appBundle ? "App-bundle launcher (open -W)" : "Visualizer";
		if (WIFEXITED(status)) {
			const int exitCode = WEXITSTATUS(status);
			if (exitCode == 0) {
				auto result = IOResult::Success(filename);
				if (appBundle)
					result.message = subject + " completed; GUI exit code is unavailable";
				return result;
			}
			return IOResult::Failure(subject + " exited with code " + std::to_string(exitCode), filename, exitCode);
		}
		if (WIFSIGNALED(status))
			return IOResult::Failure(subject + " terminated by signal " + std::to_string(WTERMSIG(status)),
				filename, -WTERMSIG(status));
		return IOResult::Failure(subject + " did not exit normally", filename);
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
