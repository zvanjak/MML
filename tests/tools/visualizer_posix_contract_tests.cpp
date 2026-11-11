#include <catch2/catch_all.hpp>

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/tools/visualizers/VisualizerProcess.h>
#endif

#ifndef _WIN32
#include <atomic>
#include <chrono>
#include <csignal>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <pthread.h>
#include <string>
#include <stdexcept>
#include <thread>
#include <vector>

namespace {
	using namespace MML;
	using namespace MML::detail::visualizer;

	struct PosixProcessFiles {
		std::filesystem::path directory = std::filesystem::temp_directory_path() /
			("mml_posix_process_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
		PosixProcessFiles() { std::filesystem::create_directory(directory); }
		~PosixProcessFiles() {
			std::error_code error;
			std::filesystem::remove_all(directory, error);
			if (error) std::cerr << "Cannot remove POSIX process fixture files: " << error.message() << '\n';
		}
		std::string executable() const { return MML_VISUALIZER_PROCESS_FIXTURE; }
	};

	void NoOpSignal(int) {}

	struct ParentAllocationActivity {
		std::atomic<bool> stop{false};
		std::thread allocator;
		ParentAllocationActivity() : allocator([this] {
			while (!stop.load()) {
				std::vector<std::string> storage(32, std::string(4096, 'x'));
				storage.back() += "allocation activity";
			}
		}) {}
		~ParentAllocationActivity() {
			stop = true;
			allocator.join();
		}
	};

	struct InterruptWait {
		struct sigaction previous {};
		std::atomic<bool> stop{false};
		std::thread sender;
		InterruptWait() {
			struct sigaction action {};
			action.sa_handler = NoOpSignal;
			sigemptyset(&action.sa_mask);
			if (sigaction(SIGUSR1, &action, &previous) != 0)
				throw std::runtime_error("Cannot install POSIX wait interruption fixture");
			const auto waitingThread = pthread_self();
			sender = std::thread([this, waitingThread] {
				while (!stop.load()) {
					if (pthread_kill(waitingThread, SIGUSR1) != 0) {
						std::cerr << "Cannot interrupt POSIX wait fixture\n";
						return;
					}
					std::this_thread::sleep_for(std::chrono::milliseconds(2));
				}
			});
		}
		~InterruptWait() {
			stop = true;
			sender.join();
			if (sigaction(SIGUSR1, &previous, nullptr) != 0)
				std::cerr << "Cannot restore POSIX signal fixture\n";
		}
	};

	void CheckReaped(const std::filesystem::path& report) {
		std::ifstream in(report);
		pid_t child = 0;
		REQUIRE(static_cast<bool>(in >> child));
		REQUIRE(child > 0);
		int status = 0;
		errno = 0;
		CHECK(waitpid(child, &status, WNOHANG) == -1);
		CHECK(errno == ECHILD);
		errno = 0;
		CHECK(kill(child, 0) == -1);
		CHECK(errno == ESRCH);
	}
}

TEST_CASE("POSIX visualizer spawn errors retain operating system diagnostics", "[Visualizer][PosixProcess]")
{
	PosixProcessFiles files;
	const auto executable = files.directory / "not executable";
	{
		std::ofstream out(executable);
		out << "not a runnable file\n";
		REQUIRE(static_cast<bool>(out));
	}
	std::filesystem::permissions(executable, std::filesystem::perms::owner_read | std::filesystem::perms::owner_write);
	const auto result = ExecuteVisualizer(executable.string(), {"diagnostic input"});
	CHECK_FALSE(result.success);
	CHECK_FALSE(result.timedOut);
	CHECK(result.error == IOError::PROCESS_FAILED);
	CHECK(result.exitCode == EACCES);
	CHECK(result.filename == "diagnostic input");
	CHECK(result.message.find("Permission denied") != std::string::npos);
}

TEST_CASE("POSIX finite visualizer waits use elapsed deadlines despite interrupted waits", "[Visualizer][PosixProcess]")
{
	PosixProcessFiles files;
	const auto report = files.directory / "sleep.pid";
	InterruptWait interrupts;
	const auto start = std::chrono::steady_clock::now();
	const auto result = ExecuteVisualizer(files.executable(), {"--sleep-report", report.string()}, 250);
	const auto elapsed = std::chrono::steady_clock::now() - start;
	CHECK_FALSE(result.success);
	CHECK(result.timedOut);
	CHECK(result.error == IOError::TIMEOUT);
	CHECK(elapsed >= std::chrono::milliseconds(250));
	CHECK(elapsed < std::chrono::seconds(2));
	CheckReaped(report);
}

TEST_CASE("POSIX visualizer child environment is isolated from its parent", "[Visualizer][PosixProcess]")
{
	PosixProcessFiles files;
	const auto report = files.directory / "environment.txt";
	std::vector<std::optional<std::string>> before;
	for (const auto* key : {"LD_LIBRARY_PATH", "QT_PLUGIN_PATH", "QT_QPA_PLATFORM_PLUGIN_PATH", "PATH"}) {
		const auto* value = std::getenv(key);
		before.emplace_back(value ? std::optional<std::string>(value) : std::nullopt);
	}
	const auto result = ExecuteVisualizer(files.executable(), {"--environment", report.string()});
	REQUIRE(result.success);
	std::vector<std::string> expected;
#ifdef __linux__
	const auto directory = std::filesystem::path(files.executable()).parent_path().string();
	expected = {directory + "/lib" + (before[0] ? ":" + *before[0] : ""),
		directory + "/plugins", directory + "/plugins/platforms", before[3].value_or("<unset>")};
#else
	for (const auto& value : before) expected.push_back(value.value_or("<unset>"));
#endif
	std::ifstream in(report);
	for (const auto& value : expected) {
		std::string line;
		REQUIRE(static_cast<bool>(std::getline(in, line)));
		CHECK(line == value);
	}
	std::size_t i = 0;
	for (const auto* key : {"LD_LIBRARY_PATH", "QT_PLUGIN_PATH", "QT_QPA_PLATFORM_PLUGIN_PATH", "PATH"}) {
		const auto* value = std::getenv(key);
		CHECK((value ? std::optional<std::string>(value) : std::nullopt) == before[i++]);
	}
}

TEST_CASE("POSIX visualizer waits preserve exits and reap graceful timeout children", "[Visualizer][PosixProcess]")
{
	PosixProcessFiles files;
	const auto report = files.directory / "child.pid";
	InterruptWait interrupts;
	SECTION("blocking wait survives signals") {
		const auto result = ExecuteVisualizer(files.executable(), {"--exit-after", report.string(), "100"});
		CHECK(result.success);
		CHECK_FALSE(result.timedOut);
		CheckReaped(report);
	}
	SECTION("finite wait completes normally") {
		const auto result = ExecuteVisualizer(files.executable(), {"--exit-after", report.string(), "100"}, 2000);
		CHECK(result.success);
		CHECK_FALSE(result.timedOut);
		CheckReaped(report);
	}
	SECTION("SIGTERM is sufficient") {
		const auto result = ExecuteVisualizer(files.executable(), {"--terminate-report", report.string()}, 250);
		CHECK_FALSE(result.success);
		CHECK(result.timedOut);
		CheckReaped(report);
	}
	SECTION("signal exit is not success") {
		const auto result = ExecuteVisualizer(files.executable(), {"--signal"});
		CHECK_FALSE(result.success);
		CHECK_FALSE(result.timedOut);
		CHECK(result.exitCode == -SIGTERM);
	}
	SECTION("zero timeout is bounded") {
		const auto start = std::chrono::steady_clock::now();
		const auto result = ExecuteVisualizer(files.executable(), {"--sleep"}, 0);
		CHECK_FALSE(result.success);
		CHECK(result.timedOut);
		CHECK(std::chrono::steady_clock::now() - start < std::chrono::seconds(2));
	}
}

TEST_CASE("App bundle command preserves literal arguments and launcher identity", "[Visualizer][PosixProcess]")
{
	PosixProcessFiles files;
	const std::string bundle = (files.directory / "Demo ;$literal.app").string();
	const std::string executable = bundle + "/Contents/MacOS/Demo";
	const std::vector<std::string> arguments{"file with spaces", "", "\"quote\\", "$(not a shell)"};
	auto command = AppBundleCommand(executable, arguments);
	REQUIRE(command == std::vector<std::string>{"/usr/bin/open", "-W", "-a", bundle, "--args",
		arguments[0], arguments[1], arguments[2], arguments[3]});
	CHECK(AppBundleCommand(executable, {}) == std::vector<std::string>{"/usr/bin/open", "-W", "-a", bundle});
	CHECK(AppBundleCommand(files.executable(), arguments).empty());
}

TEST_CASE("POSIX visualizer spawn runs alongside parent allocation activity", "[Visualizer][PosixProcess]")
{
	PosixProcessFiles files;
	const auto input = files.directory / "input.mml";
	{ std::ofstream out(input); out << 0; REQUIRE(static_cast<bool>(out)); }
	ParentAllocationActivity activity;
	for (int i = 0; i < 16; ++i) {
		const auto result = ExecuteVisualizer(files.executable(), {input.string()}, 2000);
		CHECK(result.success);
		CHECK_FALSE(result.timedOut);
	}
}

TEST_CASE("App bundle launcher shares checked timeout and status policy", "[Visualizer][PosixProcess]")
{
	PosixProcessFiles files;
	const auto report = files.directory / "launcher.pid";
	const auto executable = (files.directory / "Demo.app/Contents/MacOS/Demo").string();
	SECTION("finite timeout terminates and reaps only the owned launcher") {
		const std::vector<std::string> arguments{"--sleep-report", report.string()};
		auto command = AppBundleCommand(executable, arguments);
		command[0] = files.executable();
		InterruptWait interrupts;
		const auto start = std::chrono::steady_clock::now();
		const auto result = ExecutePosixCommand(command, arguments, 250, true);
		CHECK_FALSE(result.success);
		CHECK(result.timedOut);
		CHECK(result.error == IOError::TIMEOUT);
		CHECK(result.message.find("GUI may still be running") != std::string::npos);
		CHECK(std::chrono::steady_clock::now() - start >= std::chrono::milliseconds(250));
		CHECK(std::chrono::steady_clock::now() - start < std::chrono::seconds(2));
		CheckReaped(report);
	}
	SECTION("blocking success is launcher success") {
		const std::vector<std::string> arguments{"--exit-after", report.string(), "100"};
		auto command = AppBundleCommand(executable, arguments);
		command[0] = files.executable();
		InterruptWait interrupts;
		const auto result = ExecutePosixCommand(command, arguments, VISUALIZER_NO_TIMEOUT, true);
		CHECK(result.success);
		CHECK(result.message.find("GUI exit code is unavailable") != std::string::npos);
		CheckReaped(report);
	}
	SECTION("nonzero launcher status is retained") {
		const auto input = files.directory / "literal & O'Brien.mml";
		{ std::ofstream out(input); out << 7; }
		const std::vector<std::string> arguments{input.string(), "", "\"quote\\", "$(not a shell)"};
		auto command = AppBundleCommand(executable, arguments);
		command[0] = files.executable();
		const auto result = ExecutePosixCommand(command, arguments, 2000, true);
		CHECK_FALSE(result.success);
		CHECK_FALSE(result.timedOut);
		CHECK(result.exitCode == 7);
		CHECK(result.filename == input.string());
		CHECK(result.message.find("App-bundle launcher (open -W)") != std::string::npos);
		std::ifstream in(input.string() + ".argv");
		std::string line;
		REQUIRE(static_cast<bool>(std::getline(in, line)));
		CHECK(line == "4");
		REQUIRE(static_cast<bool>(std::getline(in, line)));
		for (const auto& argument : arguments) {
			REQUIRE(static_cast<bool>(std::getline(in, line)));
			CHECK(line == argument);
		}
	}
}
#endif
