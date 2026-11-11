#include <catch2/catch_all.hpp>

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/tools/Timer.h>
#include <mml/tools/visualizers/VisualizerProcess.h>
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>

using namespace MML;

namespace {
	struct CaptureOutput {
		std::ostringstream text;
		std::streambuf* original = std::cout.rdbuf(text.rdbuf());
		~CaptureOutput() { std::cout.rdbuf(original); }
	};

	struct DiagnosticFile {
		std::filesystem::path path = std::filesystem::temp_directory_path() /
			("mml_scoped_timer_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".log");
		~DiagnosticFile() {
			std::error_code error;
			std::filesystem::remove(path, error);
			if (error) std::cerr << "Cannot remove timer fixture diagnostic: " << error.message() << '\n';
		}
	};
}

TEST_CASE("ScopedTimer preserves silent result and elapsed controls", "[timer][scoped-timer]") {
	static_assert(std::is_nothrow_destructible_v<ScopedTimer>);
	static_assert(!std::is_copy_constructible_v<ScopedTimer> && !std::is_move_constructible_v<ScopedTimer>);
	CaptureOutput output;
	const auto before = std::chrono::steady_clock::now();
	double result = std::numeric_limits<double>::quiet_NaN();
	{
		ScopedTimer timer("silent", result);
		CHECK(timer.name() == "silent");
		const double first = timer.elapsed();
		const double second = timer.elapsed();
		CHECK(std::isfinite(first));
		CHECK(first >= 0);
		CHECK(second >= first);
	}
	const double total = std::chrono::duration<double>(std::chrono::steady_clock::now() - before).count();
	CHECK(std::isfinite(result));
	CHECK(result >= 0);
	CHECK(result <= total);
	CHECK(output.text.str().empty());
}

TEST_CASE("ScopedTimer preserves callbacks stdout format and nested reporting order", "[timer][scoped-timer]") {
	CaptureOutput output;
	int calls = 0;
	double observed = -1;
	bool callbackBeforePrint = false;
	std::string callbackName;
	{
		ScopedTimer timer("timed", [&](const std::string& name, double elapsed) {
			++calls;
			callbackName = name;
			observed = elapsed;
			callbackBeforePrint = output.text.str().empty();
		}, false);
	}
	CHECK(calls == 1);
	CHECK(callbackName == "timed");
	CHECK(callbackBeforePrint);
	CHECK(std::isfinite(observed));
	CHECK(observed >= 0);
	std::istringstream printed(output.text.str());
	std::string name, units;
	double value = -1;
	REQUIRE(printed >> name >> value >> units);
	CHECK(name == "timed:");
	CHECK(value >= 0);
	CHECK(units == "seconds");
	CHECK(output.text.str().back() == '\n');
	output.text.str("");
	{
		ScopedTimer outer("outer");
		{ ScopedTimer inner("inner"); }
	}
	CHECK(output.text.str().find("inner: ") == 0);
	CHECK(output.text.str().find("outer: ") > 0);
	output.text.str("");
	{ ScopedTimer empty("empty callback", ScopedTimer::Callback{}); }
	CHECK(output.text.str().empty());
}

TEST_CASE("ScopedTimer output and callback failures are contained and diagnosed in bounded children",
	"[timer][scoped-timer]") {
	struct Scenario { const char* mode; bool stdoutFailure; bool callbackFailure; };
	for (const auto& scenario : {
		Scenario{"throw-write", true, false}, Scenario{"short-write", true, false},
		Scenario{"throw-flush", true, false}, Scenario{"flush-state", true, false},
		Scenario{"throw-unknown", true, false}, Scenario{"bad-state", true, false},
		Scenario{"unwind-output", true, false}, Scenario{"callback-std", false, true},
		Scenario{"callback-unknown", false, true}, Scenario{"callback-print", false, true},
		Scenario{"callback-writes", false, true}, Scenario{"both", true, true},
		Scenario{"unwind-both", true, true}, Scenario{"silent-bad", false, false}}) {
		CAPTURE(scenario.mode);
		DiagnosticFile diagnostic;
		const auto result = detail::visualizer::ExecuteVisualizer(MML_SCOPED_TIMER_FIXTURE,
			{scenario.mode, diagnostic.path.string()}, 5000);
		INFO(result.message);
		CHECK_FALSE(result.timedOut);
		CHECK(result.exitCode != 86);
		CHECK(result.success);
		std::ifstream file(diagnostic.path, std::ios::binary);
		REQUIRE(file.is_open());
		std::ostringstream content;
		content << file.rdbuf();
		REQUIRE_FALSE(file.bad());
		const auto text = content.str();
		CHECK((text.find("stdout reporting failed") != std::string::npos) == scenario.stdoutFailure);
		CHECK((text.find("callback failed") != std::string::npos) == scenario.callbackFailure);
		const int expected = static_cast<int>(scenario.stdoutFailure) + static_cast<int>(scenario.callbackFailure);
		CHECK(std::count(text.begin(), text.end(), '\n') == expected);
		if (expected > 0) CHECK(text.find("ScopedTimer 'fixture': ") == 0);
		if (std::string(scenario.mode) == "callback-unknown")
			CHECK(text.find("non-standard exception") != std::string::npos);
	}
}
