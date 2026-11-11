#include <catch2/catch_all.hpp>
#include <src/visualization_examples/VisualizationResult.h>
#include <mml/tools/Visualizer.h>

#include <sstream>
#include <stdexcept>
#include <streambuf>

using namespace MML;

TEST_CASE("Visualization example reporting does not emit errors for success or data-only results",
          "[Visualizer][ExampleResults]")
{
	std::ostringstream errors;
	Examples::ReportVisualizationResult(IOResult::Success("plot.mml"), "plot", errors);
	Examples::ReportVisualizationResult(IOResult::DataSaved("plot.mml"), "export", errors);
	REQUIRE(errors.str().empty());
}

TEST_CASE("Visualization examples report current IOResult fields and preserve caller context",
          "[Visualizer][ExampleResults]")
{
	const auto failure = IOResult::Failure("child exited", "plot (one) & O'Brien.mml", 7);
	std::ostringstream errors;
	Examples::ReportVisualizationResult(failure, "curve export and launch", errors);
	REQUIRE(errors.str() ==
	        "curve export and launch failed [PROCESS_FAILED]: child exited\n"
	        "  File: plot (one) & O'Brien.mml\n"
	        "  Exit code: 7\n");
	REQUIRE(failure.filename == "plot (one) & O'Brien.mml");
	REQUIRE(failure.exitCode == 7);
	REQUIRE_FALSE(failure.success);
}

TEST_CASE("Visualization example reporting covers timeouts and missing-file errors",
          "[Visualizer][ExampleResults]")
{
	std::ostringstream timeout;
	Examples::ReportVisualizationResult(IOResult::Timeout("deadline reached", "curve.mml"), "curve", timeout);
	REQUIRE(timeout.str().find("[TIMEOUT]") != std::string::npos);
	REQUIRE(timeout.str().find("Timed out: yes") != std::string::npos);
	std::ostringstream missing;
	Examples::ReportVisualizationResult(
		IOResult::Failure("missing executable", "curve.mml", 0, IOError::FILE_NOT_OPENED), "curve", missing);
	REQUIRE(missing.str().find("[FILE_NOT_OPENED]") != std::string::npos);
	REQUIRE(missing.str().find("File: curve.mml") != std::string::npos);
	REQUIRE(missing.str().find("Exit code:") == std::string::npos);
	REQUIRE(missing.str().find("Timed out:") == std::string::npos);
}

TEST_CASE("Visualization example reporting handles an actual rejected request without a GUI",
          "[Visualizer][ExampleResults]")
{
	const auto invalid = Visualizer::VisualizeRealFunctionFromFile("");
	REQUIRE_FALSE(invalid.success);
	std::ostringstream errors;
	Examples::ReportVisualizationResult(invalid, "real function", errors);
	REQUIRE(errors.str().find("[INVALID_PARAMETERS]") != std::string::npos);
	REQUIRE(errors.str().find(invalid.message) != std::string::npos);
	REQUIRE(errors.str().find("File:") == std::string::npos);
}

TEST_CASE("Visualization example diagnostic stream failures propagate",
          "[Visualizer][ExampleResults]")
{
	class ThrowingBuffer : public std::streambuf {
		std::streamsize xsputn(const char*, std::streamsize) override {
			throw std::runtime_error("diagnostic output failed");
		}
		int_type overflow(int_type) override {
			throw std::runtime_error("diagnostic output failed");
		}
	} buffer;
	std::ostream errors(&buffer);
	errors.exceptions(std::ios::badbit);
	REQUIRE_THROWS(Examples::ReportVisualizationResult(IOResult::Failure("failed"), "plot", errors));
}
