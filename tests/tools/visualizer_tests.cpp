/**
 * @file visualizer_tests.cpp
 * @brief Infrastructure tests for Visualizer.h
 * 
 * Verifies result contracts, bundle discovery, and actual shell-free launching
 * using a small non-GUI child process.
 */

#include <catch2/catch_all.hpp>
#include <string>
#include <chrono>
#include <fstream>
#include <optional>
#include <stdexcept>
#include <array>
#include <sstream>
#include <limits>

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/tools/visualizers/VisualizerConfig.h>
#include <mml/MMLVisualizators.h>
#include <mml/base/Function.h>
#include <mml/tools/Visualizer.h>
#endif

using namespace MML;

namespace MML::Tests::Tools::VisualizerTests {

namespace {
#ifdef _WIN32
	struct ThreadErrorModeFixture {
		DWORD previous = 0;

		explicit ThreadErrorModeFixture(DWORD mode) {
			if (!SetThreadErrorMode(mode, &previous))
				throw std::runtime_error("Failed to set test thread error mode: " + std::to_string(GetLastError()));
		}
		~ThreadErrorModeFixture() { CHECK(SetThreadErrorMode(previous, nullptr) != FALSE); }
	};
#endif

	void SetTestEnvironment(const char* name, const std::optional<std::string>& value) {
#ifdef _WIN32
		const int error = _putenv_s(name, value ? value->c_str() : "");
#else
		const int error = value ? setenv(name, value->c_str(), 1) : unsetenv(name);
#endif
		if (error != 0)
			throw std::runtime_error("Failed to update test environment");
	}

	struct LauncherFixture {
		std::filesystem::path directory = std::filesystem::temp_directory_path() /
			("mml launcher (test) & O'Brien " + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
		std::filesystem::path previousRoot = detail::GetVisualizerRootState();
		VisualizerBackend previousBackend = detail::GetUnifiedBackendState();
		std::optional<std::string> previousEnvRoot;
		std::optional<std::string> previousEnvBackend;
		std::optional<std::string> previousProject;

		LauncherFixture() {
			if (const auto* value = detail::visualizer::GetEnv("MML_VISUALIZER_ROOT")) previousEnvRoot = value;
			if (const auto* value = detail::visualizer::GetEnv("MML_VISUALIZER_BACKEND")) previousEnvBackend = value;
			if (const auto* value = detail::visualizer::GetEnv("MML_PROJECT_PATH")) previousProject = value;
			SetTestEnvironment("MML_VISUALIZER_ROOT", std::nullopt);
			SetTestEnvironment("MML_VISUALIZER_BACKEND", std::nullopt);
			SetTestEnvironment("MML_PROJECT_PATH", std::nullopt);
			std::filesystem::create_directories(directory);
			SetVisualizerRoot(directory);
			detail::GetUnifiedBackendState() = VisualizerBackend::FLTK;
		}

		~LauncherFixture() {
			SetVisualizerRoot(previousRoot);
			detail::GetUnifiedBackendState() = previousBackend;
			SetTestEnvironment("MML_VISUALIZER_ROOT", previousEnvRoot);
			SetTestEnvironment("MML_VISUALIZER_BACKEND", previousEnvBackend);
			SetTestEnvironment("MML_PROJECT_PATH", previousProject);
			std::error_code error;
			std::filesystem::remove_all(directory, error);
			if (error)
				std::cerr << "Failed to clean launcher test directory: " << error.message() << '\n';
		}

		std::filesystem::path Install(const std::string& name = "MML_RealFunctionVisualizer") {
			const auto executable = std::filesystem::path(GetVisualizerPath(name, VisualizerBackend::FLTK));
			std::filesystem::create_directories(executable.parent_path());
			std::filesystem::copy_file(MML_VISUALIZER_PROCESS_FIXTURE, executable);
			return executable;
		}

		std::filesystem::path Data(const std::string& name, int exitCode = 0) {
			const auto path = directory / name;
			std::ofstream file(path);
			file << exitCode;
			if (!file)
				throw std::runtime_error("Failed to write launcher test data");
			return path;
		}
	};

	struct GeneratedOutput {
		std::filesystem::path path = std::filesystem::path(GetResultFilesPath()) /
			("visualizer_workflow_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".mml");

		GeneratedOutput() {
			std::filesystem::create_directories(path.parent_path());
		}

		~GeneratedOutput() {
			for (const auto& file : {path, std::filesystem::path(path.string() + ".argv")}) {
				std::error_code error;
				std::filesystem::remove(file, error);
				if (error)
					std::cerr << "Failed to remove visualizer test output: " << error.message() << '\n';
			}
		}
	};

	std::string ReadGeneratedOutput(const std::filesystem::path& path) {
		std::ifstream in(path);
		REQUIRE(in.is_open());
		std::ostringstream out;
		out << in.rdbuf();
		REQUIRE_FALSE(in.bad());
		return out.str();
	}

	void CheckVectorRows(const std::string& text, int dimension, int expectedRows, Real magnitude) {
		std::istringstream input(text);
		std::string line;
		for (int header = 0; header < 3; ++header) REQUIRE(static_cast<bool>(std::getline(input, line)));
		int rows = 0;
		while (std::getline(input, line)) {
			std::istringstream values(line);
			Real value;
			for (int coordinate = 0; coordinate < dimension; ++coordinate) REQUIRE(static_cast<bool>(values >> value));
			for (int component = 0; component < dimension; ++component) {
				REQUIRE(static_cast<bool>(values >> value));
				CHECK(value == (component == 0 ? magnitude : Real(0)));
			}
			std::string extra;
			CHECK_FALSE(static_cast<bool>(values >> extra));
			++rows;
		}
		CHECK(rows == expectedRows);
	}
}

TEST_CASE("Ordinary 3D vector visualization retains all sampled magnitudes", "[Visualizer][VectorFiltering]") {
	LauncherFixture fixture;
	fixture.Install("MML_VectorField3D_Visualizer");
	for (Real magnitude : {Real(0), Real(3), Real(4), Real(40), std::numeric_limits<Real>::max() / 4}) {
		CAPTURE(magnitude);
		GeneratedOutput output, reference;
		int evaluations = 0;
		VectorFunctionFromStdFunc<3> field([&](const VectorN<Real, 3>&) {
			++evaluations;
			return VectorN<Real, 3>{magnitude, 0, 0};
		});
		const auto result = Visualizer::VisualizeVectorField3DCartesian(
			field, "constant", 0, 1, 2, 0, 1, 2, 0, 1, 2, output.path.filename().string());
		INFO(result.message);
		REQUIRE(result.success);
		CHECK(evaluations == 8);
		CHECK(std::filesystem::path(result.filename) == std::filesystem::weakly_canonical(output.path));
		CHECK(std::filesystem::exists(output.path.string() + ".argv"));
		REQUIRE(Serializer::SaveVectorFunc3DCartesian(field, "constant", 0, 1, 2, 0, 1, 2, 0, 1, 2, reference.path.string()).success);
		const auto text = ReadGeneratedOutput(output.path);
		CHECK(text == ReadGeneratedOutput(reference.path));
		CheckVectorRows(text, 3, 8, magnitude);
	}
}

TEST_CASE("Ordinary 2D vector visualization remains unfiltered", "[Visualizer][VectorFiltering]") {
	LauncherFixture fixture;
	fixture.Install("MML_VectorField2D_Visualizer");
	GeneratedOutput output, reference;
	VectorFunctionFromStdFunc<2> field([](const VectorN<Real, 2>&) { return VectorN<Real, 2>{4, 0}; });
	const auto result = Visualizer::VisualizeVectorField2DCartesian(
		field, "constant", 0, 1, 2, 0, 1, 2, output.path.filename().string());
	INFO(result.message);
	REQUIRE(result.success);
	CHECK(std::filesystem::exists(output.path.string() + ".argv"));
	REQUIRE(Serializer::SaveVectorFunc2DCartesian(field, "constant", 0, 1, 2, 0, 1, 2, reference.path.string()).success);
	const auto text = ReadGeneratedOutput(output.path);
	CHECK(text == ReadGeneratedOutput(reference.path));
	CheckVectorRows(text, 2, 4, 4);
}

TEST_CASE("Explicit 3D vector cutoff preserves exclusive serializer semantics", "[Visualizer][VectorFiltering]") {
	LauncherFixture fixture;
	fixture.Install("MML_VectorField3D_Visualizer");
	VectorFunctionFromStdFunc<3> field([](const VectorN<Real, 3>&) { return VectorN<Real, 3>{4, 0, 0}; });
	for (Real cutoff : {Real(3), Real(4), Real(5), std::numeric_limits<Real>::infinity()}) {
		CAPTURE(cutoff);
		GeneratedOutput output, reference;
		const auto result = Visualizer::VisualizeVectorField3DCartesian(
			field, "constant", 0, 1, 2, 0, 1, 2, 0, 1, 2, output.path.filename().string(), cutoff);
		INFO(result.message);
		REQUIRE(result.success);
		CHECK(std::filesystem::exists(output.path.string() + ".argv"));
		REQUIRE(Serializer::SaveVectorFunc3DCartesian(
			field, "constant", 0, 1, 2, 0, 1, 2, 0, 1, 2, reference.path.string(), cutoff).success);
		const auto text = ReadGeneratedOutput(output.path);
		CHECK(text == ReadGeneratedOutput(reference.path));
		CheckVectorRows(text, 3, cutoff > 4 ? 8 : 0, 4);
	}
}

TEST_CASE("Explicit 2D vector cutoff matches the serializer", "[Visualizer][VectorFiltering]") {
	LauncherFixture fixture;
	fixture.Install("MML_VectorField2D_Visualizer");
	VectorFunctionFromStdFunc<2> field([](const VectorN<Real, 2>&) { return VectorN<Real, 2>{4, 0}; });
	for (Real cutoff : {Real(3), Real(4), Real(5), std::numeric_limits<Real>::infinity()}) {
		CAPTURE(cutoff);
		GeneratedOutput output, reference;
		const auto result = Visualizer::VisualizeVectorField2DCartesian(
			field, "constant", 0, 1, 2, 0, 1, 2, output.path.filename().string(), cutoff);
		INFO(result.message);
		REQUIRE(result.success);
		CHECK(std::filesystem::exists(output.path.string() + ".argv"));
		REQUIRE(Serializer::SaveVectorFunc2DCartesian(field, "constant", 0, 1, 2, 0, 1, 2, reference.path.string(), cutoff).success);
		const auto text = ReadGeneratedOutput(output.path);
		CHECK(text == ReadGeneratedOutput(reference.path));
		CheckVectorRows(text, 2, cutoff > 4 ? 4 : 0, 4);
	}
}

TEST_CASE("Explicit vector cutoff selectively retains samples below its bound", "[Visualizer][VectorFiltering]") {
	LauncherFixture fixture;
	fixture.Install("MML_VectorField2D_Visualizer");
	fixture.Install("MML_VectorField3D_Visualizer");
	VectorFunctionFromStdFunc<2> field2([](const VectorN<Real, 2>& p) { return VectorN<Real, 2>{2 + 2 * p[0], 0}; });
	VectorFunctionFromStdFunc<3> field3([](const VectorN<Real, 3>& p) { return VectorN<Real, 3>{2 + 2 * p[0], 0, 0}; });
	GeneratedOutput output2, output3;
	REQUIRE(Visualizer::VisualizeVectorField2DCartesian(
		field2, "mixed", 0, 1, 2, 0, 1, 2, output2.path.filename().string(), 3).success);
	REQUIRE(Visualizer::VisualizeVectorField3DCartesian(
		field3, "mixed", 0, 1, 2, 0, 1, 2, 0, 1, 2, output3.path.filename().string(), 3).success);
	CheckVectorRows(ReadGeneratedOutput(output2.path), 2, 2, 2);
	CheckVectorRows(ReadGeneratedOutput(output3.path), 3, 4, 2);
}

TEST_CASE("Invalid vector cutoffs reject before evaluation file writes or launching", "[Visualizer][VectorFiltering]") {
	LauncherFixture fixture;
	fixture.Install("MML_VectorField2D_Visualizer");
	fixture.Install("MML_VectorField3D_Visualizer");
	int evaluations = 0;
	VectorFunctionFromStdFunc<2> field2([&](const VectorN<Real, 2>&) { ++evaluations; return VectorN<Real, 2>{4, 0}; });
	VectorFunctionFromStdFunc<3> field3([&](const VectorN<Real, 3>&) { ++evaluations; return VectorN<Real, 3>{4, 0, 0}; });
	for (int dimension : {2, 3}) {
		for (Real cutoff : {Real(0), Real(-1), -std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
			for (bool facade : {false, true}) {
				for (bool seed : {false, true}) {
					CAPTURE(dimension, cutoff, facade, seed);
					GeneratedOutput output;
					if (seed) {
						std::ofstream out(output.path);
						out << "SENTINEL";
						out.close();
						REQUIRE_FALSE(out.fail());
					}
					IOResult result;
					if (dimension == 2) {
						result = facade ? Visualizer::VisualizeVectorField2DCartesian(
							field2, "invalid", 0, 1, 2, 0, 1, 2, output.path.filename().string(), cutoff) :
							Serializer::SaveVectorFunc2DCartesian(field2, "invalid", 0, 1, 2, 0, 1, 2, output.path.string(), cutoff);
					} else {
						result = facade ? Visualizer::VisualizeVectorField3DCartesian(
							field3, "invalid", 0, 1, 2, 0, 1, 2, 0, 1, 2, output.path.filename().string(), cutoff) :
							Serializer::SaveVectorFunc3DCartesian(field3, "invalid", 0, 1, 2, 0, 1, 2, 0, 1, 2, output.path.string(), cutoff);
					}
					CHECK_FALSE(result.success);
					CHECK(result.error == IOError::INVALID_PARAMETERS);
					CHECK_FALSE(result.message.empty());
					CHECK(std::filesystem::path(result.filename) == std::filesystem::weakly_canonical(output.path));
					CHECK(evaluations == 0);
					CHECK_FALSE(std::filesystem::exists(output.path.string() + ".argv"));
					if (seed) CHECK(ReadGeneratedOutput(output.path) == "SENTINEL");
					else CHECK_FALSE(std::filesystem::exists(output.path));
				}
			}
		}
	}
}

TEST_CASE("Filtered vector facade propagates invalid sampling requests", "[Visualizer][VectorFiltering]") {
	LauncherFixture fixture;
	fixture.Install("MML_VectorField2D_Visualizer");
	fixture.Install("MML_VectorField3D_Visualizer");
	int evaluations = 0;
	VectorFunctionFromStdFunc<2> field2([&](const VectorN<Real, 2>&) { ++evaluations; return VectorN<Real, 2>{4, 0}; });
	VectorFunctionFromStdFunc<3> field3([&](const VectorN<Real, 3>&) { ++evaluations; return VectorN<Real, 3>{4, 0, 0}; });
	GeneratedOutput output2, output3;
	const auto result2 = Visualizer::VisualizeVectorField2DCartesian(
		field2, "invalid", 0, 1, 1, 0, 1, 2, output2.path.filename().string(), 5);
	const auto result3 = Visualizer::VisualizeVectorField3DCartesian(
		field3, "invalid", 0, 1, 2, 0, 1, 2, 0, 1, 1, output3.path.filename().string(), 5);
	CHECK(result2.error == IOError::INVALID_PARAMETERS);
	CHECK(result3.error == IOError::INVALID_PARAMETERS);
	CHECK(evaluations == 0);
	CHECK_FALSE(std::filesystem::exists(output2.path));
	CHECK_FALSE(std::filesystem::exists(output3.path));
	CHECK_FALSE(std::filesystem::exists(output2.path.string() + ".argv"));
	CHECK_FALSE(std::filesystem::exists(output3.path.string() + ".argv"));
}

TEST_CASE("Visualizer catalog drives all legacy path getters", "[Visualizer][Launcher][Catalog]") {
	LauncherFixture fixture;
	using Kind = detail::visualizer::Kind;
	struct Entry { Kind kind; const char* name; std::string (*getter)(); };
	const std::array<Entry, 12> expected{{
		{Kind::RealFunction, "MML_RealFunctionVisualizer", GetRealFuncVisualizerPath},
		{Kind::ScalarFunction2D, "MML_ScalarFunction2D_Visualizer", GetSurfaceVisualizerPath},
		{Kind::Curve2D, "MML_ParametricCurve2D_Visualizer", GetParametricCurve2DVisualizerPath},
		{Kind::Curve3D, "MML_ParametricCurve3D_Visualizer", GetParametricCurve3DVisualizerPath},
		{Kind::VectorField2D, "MML_VectorField2D_Visualizer", GetVectorField2DVisualizerPath},
		{Kind::VectorField3D, "MML_VectorField3D_Visualizer", GetVectorField3DVisualizerPath},
		{Kind::Particles2D, "MML_ParticleVisualizer2D", GetParticle2DVisualizerPath},
		{Kind::Particles3D, "MML_ParticleVisualizer3D", GetParticle3DVisualizerPath},
		{Kind::ParametricSurface, "MML_ParametricSurface_Visualizer", GetParametricSurfaceVisualizerPath},
		{Kind::ScalarFunction3D, "MML_ScalarFunction3D_Visualizer", GetScalarFunction3DVisualizerPath},
		{Kind::RigidBody, "MML_RigidBodyMovement_Visualizer", GetRigidBodyVisualizerPath},
		{Kind::World, "MML_WorldVisualizer", GetWorldVisualizerPath}
	}};
	REQUIRE(detail::visualizer::ExecutableNames.size() == expected.size());
	for (const auto& entry : expected) {
		INFO(entry.name);
		REQUIRE(std::string(detail::visualizer::ExecutableName(entry.kind)) == entry.name);
		REQUIRE(entry.getter() == GetVisualizerPath(entry.name));
		REQUIRE(detail::visualizer::ExecutablePath(entry.kind) == entry.getter());
	}
	REQUIRE_THROWS_AS(detail::visualizer::ExecutableName(static_cast<Kind>(-1)), VisualizerError);
	REQUIRE_THROWS_AS(detail::visualizer::ExecutableName(static_cast<Kind>(expected.size())), VisualizerError);
}

TEST_CASE("Every visualizer file facade dispatches to its catalog executable", "[Visualizer][Launcher][Catalog]") {
	struct Entry { const char* name; IOResult (*launch)(const std::string&); };
	const std::array<Entry, 13> entries{{
		{"MML_RealFunctionVisualizer", [](const std::string& f) { return Visualizer::VisualizeRealFunctionFromFile(f); }},
		{"MML_ScalarFunction2D_Visualizer", [](const std::string& f) { return Visualizer::VisualizeScalarFunc2DFromFile(f); }},
		{"MML_ParametricCurve2D_Visualizer", [](const std::string& f) { return Visualizer::VisualizeMultiParamCurve2D(std::vector<std::string>{f}); }},
		{"MML_ParametricCurve3D_Visualizer", [](const std::string& f) { return Visualizer::VisualizeMultiParamCurve3D(std::vector<std::string>{f}); }},
		{"MML_VectorField2D_Visualizer", [](const std::string& f) { return Visualizer::VisualizeVectorField2DFromFile(f); }},
		{"MML_VectorField3D_Visualizer", [](const std::string& f) { return Visualizer::VisualizeVectorField3DFromFile(f); }},
		{"MML_ParticleVisualizer2D", [](const std::string& f) { return Visualizer::VisualizeParticleSimulation2D(f); }},
		{"MML_ParticleVisualizer3D", [](const std::string& f) { return Visualizer::VisualizeParticleSimulation3D(f); }},
		{"MML_ParametricSurface_Visualizer", [](const std::string& f) { return Visualizer::VisualizeParametricSurfaceFromFile(f); }},
		{"MML_ScalarFunction3D_Visualizer", [](const std::string& f) { return Visualizer::VisualizeScalarFunc3DFromFile(f); }},
		{"MML_RigidBodyMovement_Visualizer", [](const std::string& f) { return Visualizer::VisualizeRigidBodyTrajectory(f); }},
		{"MML_RigidBodyMovement_Visualizer", [](const std::string& f) { return Visualizer::VisualizeRigidBodyTrajectory(std::vector<std::string>{f}); }},
		{"MML_WorldVisualizer", [](const std::string& f) { return Visualizer::VisualizeWorldSceneFromFile(f); }}
	}};
	for (const auto& entry : entries) {
		INFO(entry.name);
		LauncherFixture fixture;
		fixture.Install(entry.name);
		const auto input = fixture.Data("dispatch (one).mml");
		const auto result = entry.launch(input.string());
		INFO(result.message);
		REQUIRE(result.success);
		REQUIRE(std::filesystem::path(result.filename) == std::filesystem::weakly_canonical(input));
		REQUIRE(std::filesystem::is_regular_file(input.string() + ".argv"));
	}
}

TEST_CASE("Visualizer file helpers preserve containment and filesystem diagnostics", "[Visualizer][Files]") {
	LauncherFixture fixture;
	const auto root = fixture.directory / "results";
	const auto output = detail::visualizer::MakeOutputPath(
		(std::filesystem::path("nested") / ".plot (one).mml").string(), root.string());
	REQUIRE(output.success);
	REQUIRE(std::filesystem::path(output.filename).parent_path() == std::filesystem::weakly_canonical(root));
	REQUIRE(std::filesystem::path(output.filename).filename() == "plot__one_.mml");
	const auto suffix = detail::visualizer::MakeOutputPath("curve", root.string(), "_2.mml");
	REQUIRE(std::filesystem::path(suffix.filename).filename() == "curve_2.mml");
	REQUIRE(std::filesystem::path(detail::visualizer::MakeOutputPath("...", root.string()).filename).filename() == "output");
	const auto traversal = detail::visualizer::MakeOutputPath(
		"curve", root.string(), std::string(1, std::filesystem::path::preferred_separator) +
		(std::filesystem::path("..") / ".." / "escaped.mml").string());
	REQUIRE_FALSE(traversal.success);
	REQUIRE(traversal.error == IOError::INVALID_PARAMETERS);
	const auto nul = detail::visualizer::MakeOutputPath(std::string("bad\0path", 8), root.string());
	REQUIRE_FALSE(nul.success);
	REQUIRE(nul.error == IOError::INVALID_PARAMETERS);
	const auto blocked = fixture.Data("blocking");
	const auto failure = detail::visualizer::MakeOutputPath("plot.mml", (blocked / "results").string());
	REQUIRE_FALSE(failure.success);
	REQUIRE(failure.error == IOError::FILE_NOT_OPENED);
	REQUIRE(failure.exitCode != 0);
	REQUIRE(failure.message.find("results directory") != std::string::npos);
}

TEST_CASE("Visualizer input resolution preserves precedence and atomic multi-file results", "[Visualizer][Files]") {
	LauncherFixture fixture;
	const auto root = fixture.directory / "results";
	std::filesystem::create_directory(root);
	const std::string name = "lookup_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".mml";
	const auto fallback = fixture.Data((std::filesystem::path("results") / name).string());
	const auto found = detail::visualizer::ResolveInputFile(name, root.string());
	REQUIRE(found.success);
	REQUIRE(std::filesystem::path(found.filename) == std::filesystem::weakly_canonical(fallback));
	const auto direct = fixture.Data(name);
	const auto preferred = detail::visualizer::ResolveInputFile(direct.string(), root.string());
	REQUIRE(preferred.success);
	REQUIRE(std::filesystem::path(preferred.filename) == std::filesystem::weakly_canonical(direct));
	const auto directory = detail::visualizer::ResolveInputFile(root.string(), root.string());
	REQUIRE_FALSE(directory.success);
	REQUIRE(directory.message.find("not a regular file") != std::string::npos);
	std::vector<std::string> resolved{"unchanged"};
	const auto missing = detail::visualizer::ResolveInputFiles({direct.string(), "missing.mml"}, root.string(), resolved);
	REQUIRE_FALSE(missing.success);
	REQUIRE(missing.filename == "missing.mml");
	REQUIRE(resolved == std::vector<std::string>{"unchanged"});
	const auto empty = detail::visualizer::ResolveInputFiles({}, root.string(), resolved);
	REQUIRE_FALSE(empty.success);
	REQUIRE(empty.error == IOError::INVALID_PARAMETERS);
	REQUIRE_FALSE(detail::visualizer::ResolveInputFile("", root.string()).success);
	REQUIRE_FALSE(detail::visualizer::ResolveInputFile(std::string("bad\0path", 8), root.string()).success);
	REQUIRE_FALSE(Visualizer::VisualizeMultiParamCurve2D(std::vector<std::string>{}).success);
	const auto noCurves = Visualizer::VisualizeMultiParamCurve2D(
		std::vector<IRealToVectorFunction<2>*>{}, "empty", Real(0), Real(1), 4, "empty.mml");
	REQUIRE_FALSE(noCurves.success);
	REQUIRE(noCurves.error == IOError::INVALID_PARAMETERS);
}

TEST_CASE("Visualizer save-and-launch propagates serialization failure before starting a process", "[Visualizer][Launcher]") {
	LauncherFixture fixture;
	fixture.Install();
	GeneratedOutput output;
	int evaluations = 0;
	RealFunctionFromStdFunc function([&](Real x) { ++evaluations; return x * x; });
	const auto invalid = Visualizer::VisualizeRealFunction(function, "invalid", Real(0), Real(1), 0, output.path.filename().string());
	REQUIRE_FALSE(invalid.success);
	REQUIRE(invalid.error == IOError::INVALID_PARAMETERS);
	REQUIRE(evaluations == 0);
	REQUIRE_FALSE(std::filesystem::exists(output.path.string() + ".argv"));
	const auto nul = Visualizer::VisualizeRealFunction(function, "invalid path", Real(0), Real(1), 4, std::string("bad\0path", 8));
	REQUIRE_FALSE(nul.success);
	REQUIRE(nul.error == IOError::INVALID_PARAMETERS);
	REQUIRE(evaluations == 0);
	const auto valid = Visualizer::VisualizeRealFunction(function, "square", Real(0), Real(1), 4, output.path.filename().string());
	INFO(valid.message);
	REQUIRE(valid.success);
	REQUIRE(evaluations == 4);
	REQUIRE(std::filesystem::path(valid.filename) == std::filesystem::weakly_canonical(output.path));
	REQUIRE(std::filesystem::is_regular_file(output.path.string() + ".argv"));
}

TEST_CASE("Rigid-body visualizers resolve relative inputs before changing child working directory", "[Visualizer][Launcher][Files]") {
	LauncherFixture fixture;
	fixture.Install("MML_RigidBodyMovement_Visualizer");
	const auto input = fixture.Data("relative.mml");
	struct RestoreDirectory {
		std::filesystem::path previous = std::filesystem::current_path();
		~RestoreDirectory() { std::filesystem::current_path(previous); }
	} restore;
	std::filesystem::current_path(fixture.directory);
	for (bool multiple : {false, true}) {
		const auto result = multiple
			? Visualizer::VisualizeRigidBodyTrajectory(std::vector<std::string>{"relative.mml"})
			: Visualizer::VisualizeRigidBodyTrajectory("relative.mml");
		INFO(result.message);
		REQUIRE(result.success);
		REQUIRE(std::filesystem::path(result.filename) == std::filesystem::weakly_canonical(input));
	}
}

TEST_CASE("Visualizer process layer preserves quoting, failures and timeout cleanup", "[Visualizer][Launcher][Process]") {
	LauncherFixture fixture;
	const auto executable = fixture.Install();
	const auto input = fixture.Data("quoting.mml");
	const std::vector<std::string> arguments{input.string(), "", "quote\"inside", "trailing\\", "double\\\\\"quote", "&;$`( )'"};
	const auto result = detail::visualizer::ExecuteVisualizer(executable.string(), arguments);
	INFO(result.message);
	REQUIRE(result.success);
	std::ifstream report(input.string() + ".argv");
	std::string line;
	std::getline(report, line);
	REQUIRE(line == std::to_string(arguments.size()));
	std::getline(report, line);
	for (const auto& argument : arguments) {
		std::getline(report, line);
		REQUIRE(line == argument);
	}
	const auto invalidTimeout = detail::visualizer::ExecuteVisualizer(executable.string(), {input.string()}, -2);
	REQUIRE_FALSE(invalidTimeout.success);
	REQUIRE(invalidTimeout.error == IOError::INVALID_PARAMETERS);
	const auto nul = detail::visualizer::ExecuteVisualizer(executable.string(), {std::string("bad\0arg", 7)});
	REQUIRE_FALSE(nul.success);
	REQUIRE(nul.error == IOError::INVALID_PARAMETERS);
	const auto timeout = detail::visualizer::ExecuteVisualizer(executable.string(), {"--sleep"}, 25);
	REQUIRE_FALSE(timeout.success);
	REQUIRE(timeout.timedOut);
	REQUIRE(timeout.error == IOError::TIMEOUT);
}

TEST_CASE("Visualizer root precedence and deployment layout", "[Visualizer][Launcher]") {
	LauncherFixture fixture;
	const auto environmentRoot = fixture.directory / "environment";
	SetTestEnvironment("MML_VISUALIZER_ROOT", environmentRoot.string());
	REQUIRE(std::filesystem::path(GetVisualizerRoot()) == fixture.directory);
	SetVisualizerRoot({});
	REQUIRE(std::filesystem::path(GetVisualizerRoot()) == environmentRoot);
	SetTestEnvironment("MML_VISUALIZER_ROOT", std::nullopt);
	SetTestEnvironment("MML_PROJECT_PATH", fixture.directory.string());
	REQUIRE(std::filesystem::path(GetVisualizerRoot()) == std::filesystem::current_path() / "visualizers");

	SetVisualizerRoot("relative bundle");
	REQUIRE(std::filesystem::path(GetVisualizerRoot()) == std::filesystem::current_path() / "relative bundle");
	REQUIRE(GetVisualizerPath("MML_RealFunctionVisualizer", VisualizerBackend::FLTK).find(GetVisualizerRoot()) == 0);
	const auto savedRoot = GetVisualizerRoot();
	REQUIRE_THROWS_AS(SetVisualizerRoot(std::string("invalid\0root", 12)), VisualizerError);
	REQUIRE(GetVisualizerRoot() == savedRoot);
}

TEST_CASE("Visualizer backend paths use the current bundle without cached discovery", "[Visualizer][Launcher]") {
	LauncherFixture fixture;
	const auto first = GetVisualizerPath("MML_RealFunctionVisualizer", VisualizerBackend::Qt);
	SetVisualizerRoot(fixture.directory / "other bundle");
	const auto second = GetVisualizerPath("MML_RealFunctionVisualizer", VisualizerBackend::Qt);
	REQUIRE(first != second);
#ifdef _WIN32
	REQUIRE(std::filesystem::path(second) == fixture.directory / "other bundle" / "win" / "Qt" / "MML_RealFunctionVisualizer.exe");
	REQUIRE(std::filesystem::path(GetWPFVisualizerPath("MML_RealFunctionVisualizer", "MML_RealFunctionVisualizer")) ==
		fixture.directory / "other bundle" / "win" / "WPF" / "MML_RealFunctionVisualizer" / "MML_RealFunctionVisualizer.exe");
#elif defined(__APPLE__)
	REQUIRE(std::filesystem::path(second) == fixture.directory / "other bundle" / "mac" / "Qt" /
		"MML_RealFunctionVisualizer.app" / "Contents" / "MacOS" / "MML_RealFunctionVisualizer");
#else
	REQUIRE(std::filesystem::path(second) == fixture.directory / "other bundle" / "linux" / "Qt" / "MML_RealFunctionVisualizer");
#endif
}

TEST_CASE("Visualizer availability requires a supported regular executable", "[Visualizer][Launcher]") {
	LauncherFixture fixture;
	REQUIRE_FALSE(IsBackendAvailable(VisualizerBackend::Auto));
	const auto executable = std::filesystem::path(GetVisualizerPath("MML_RealFunctionVisualizer", VisualizerBackend::FLTK));
	std::filesystem::create_directories(executable);
	REQUIRE_FALSE(IsBackendAvailable(VisualizerBackend::FLTK));
	std::filesystem::remove(executable);
	fixture.Install();
	REQUIRE(IsBackendAvailable(VisualizerBackend::FLTK));
	REQUIRE(IsBackendAvailable(VisualizerBackend::Auto));
	REQUIRE_FALSE(IsBackendAvailable(VisualizerBackend::Qt));
	SetTestEnvironment("MML_VISUALIZER_BACKEND", "AUTO");
	REQUIRE_FALSE(IsBackendAvailable(VisualizerBackend::Auto));
#ifdef _WIN32
	REQUIRE(GetVisualizerBackend() == VisualizerBackend::WPF);
#else
	REQUIRE(GetVisualizerBackend() == VisualizerBackend::Qt);
#endif
	SetTestEnvironment("MML_VISUALIZER_BACKEND", "FLTK");
	REQUIRE(IsBackendAvailable(VisualizerBackend::Auto));
#ifndef _WIN32
	REQUIRE_FALSE(IsBackendAvailable(VisualizerBackend::WPF));
#endif
}

TEST_CASE("Visualizer launches literal metacharacter paths and reports child exit status", "[Visualizer][Launcher]") {
	LauncherFixture fixture;
	const auto executable = fixture.Install();
	const auto input = fixture.Data("data (curve) & O'Brien;$`plot.mml");
	const auto result = Visualizer::VisualizeRealFunctionFromFile(input.string());
	INFO(result.message);
	REQUIRE(result.success);
	std::ifstream report(input.string() + ".argv");
	std::string count, workDir, argument;
	std::getline(report, count);
	std::getline(report, workDir);
	std::getline(report, argument);
	REQUIRE(count == "1");
#ifdef _WIN32
	REQUIRE(std::filesystem::path(workDir) == executable.parent_path());
#endif
	REQUIRE(std::filesystem::path(argument) == std::filesystem::weakly_canonical(input));

	const auto failingInput = fixture.Data("failure (7).mml", 7);
	const auto failure = Visualizer::VisualizeRealFunctionFromFile(failingInput.string());
	REQUIRE_FALSE(failure.success);
	REQUIRE(failure.exitCode == 7);
	REQUIRE(failure.message.find("7") != std::string::npos);
	REQUIRE(std::filesystem::path(failure.filename) == std::filesystem::weakly_canonical(failingInput));
}

TEST_CASE("Visualizer preserves separate arguments and checks the exact executable", "[Visualizer][Launcher]") {
	LauncherFixture fixture;
	fixture.Install("MML_ParametricCurve3D_Visualizer");
	const auto first = fixture.Data("first (curve) & one.mml");
	const auto second = fixture.Data("second O'Brien;$`curve.mml");
	const auto result = Visualizer::VisualizeMultiParamCurve3D(std::vector<std::string>{first.string(), second.string()});
	INFO(result.message);
	REQUIRE(result.success);
	std::ifstream report(first.string() + ".argv");
	std::string line;
	std::getline(report, line);
	REQUIRE(line == "2");
	std::getline(report, line);
	std::getline(report, line);
	REQUIRE(std::filesystem::path(line) == std::filesystem::weakly_canonical(first));
	std::getline(report, line);
	REQUIRE(std::filesystem::path(line) == std::filesystem::weakly_canonical(second));

	REQUIRE(IsBackendAvailable(VisualizerBackend::FLTK));
	const auto missing = Visualizer::VisualizeRealFunctionFromFile(first.string());
	REQUIRE_FALSE(missing.success);
	REQUIRE(missing.error == IOError::FILE_NOT_OPENED);
	REQUIRE(missing.message.find("executable not found") != std::string::npos);
}

TEST_CASE("Visualizer surfaces process-start failures for an invalid executable", "[Visualizer][Launcher]") {
	LauncherFixture fixture;
	const auto executable = std::filesystem::path(GetVisualizerPath("MML_RealFunctionVisualizer", VisualizerBackend::FLTK));
	std::filesystem::create_directories(executable.parent_path());
	std::ofstream invalid(executable);
	invalid << "not an executable";
	invalid.close();
	REQUIRE(IsBackendAvailable(VisualizerBackend::FLTK));
	const auto input = fixture.Data("start failure.mml");
#ifdef _WIN32
	const DWORD processMode = GetErrorMode();
	DWORD expectedError = ERROR_SUCCESS;
	{
		ThreadErrorModeFixture nativeMode(GetThreadErrorMode() | SEM_FAILCRITICALERRORS);
		STARTUPINFOA startup = {};
		startup.cb = sizeof(startup);
		PROCESS_INFORMATION process = {};
		const BOOL created = CreateProcessA(executable.string().c_str(), nullptr, nullptr, nullptr,
			FALSE, 0, nullptr, nullptr, &startup, &process);
		expectedError = created ? ERROR_SUCCESS : GetLastError();
		if (created) {
			CloseHandle(process.hThread);
			CloseHandle(process.hProcess);
		}
		REQUIRE_FALSE(created);
	}
	REQUIRE(expectedError != ERROR_SUCCESS);
	for (DWORD mode : {DWORD(0), DWORD(SEM_NOOPENFILEERRORBOX), DWORD(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX)}) {
		CAPTURE(mode);
		ThreadErrorModeFixture errorMode(mode);
#endif
	const auto result = Visualizer::VisualizeRealFunctionFromFile(input.string());
	REQUIRE_FALSE(result.success);
	REQUIRE(result.exitCode != 0);
	REQUIRE_FALSE(result.message.empty());
	REQUIRE(std::filesystem::path(result.filename) == std::filesystem::weakly_canonical(input));
#ifdef _WIN32
		CHECK(result.error == IOError::PROCESS_FAILED);
		CHECK_FALSE(result.timedOut);
		CHECK(result.exitCode == static_cast<int>(expectedError));
		CHECK(result.message.find(std::to_string(expectedError)) != std::string::npos);
		CHECK(GetThreadErrorMode() == mode);
		CHECK(GetErrorMode() == processMode);
	}
#endif
}

#ifdef _WIN32
TEST_CASE("Windows visualizer launches preserve caller error modes on success failure and timeout", "[Visualizer][Launcher][Process]") {
	LauncherFixture fixture;
	const auto executable = fixture.Install();
	const auto input = fixture.Data("error mode success.mml");
	const auto failingInput = fixture.Data("error mode exit7.mml", 7);
	const DWORD processMode = GetErrorMode();
	for (DWORD mode : {DWORD(0), DWORD(SEM_NOOPENFILEERRORBOX), DWORD(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX)}) {
		CAPTURE(mode);
		ThreadErrorModeFixture errorMode(mode);
		const auto success = detail::visualizer::ExecuteVisualizer(executable.string(), {input.string()});
		INFO(success.message);
		CHECK(success.success);
		CHECK(GetThreadErrorMode() == mode);
		CHECK(GetErrorMode() == processMode);
		const auto failure = detail::visualizer::ExecuteVisualizer(executable.string(), {failingInput.string()});
		CHECK_FALSE(failure.success);
		CHECK(failure.exitCode == 7);
		CHECK(GetThreadErrorMode() == mode);
		CHECK(GetErrorMode() == processMode);
		const auto timeout = detail::visualizer::ExecuteVisualizer(executable.string(), {"--sleep"}, 25);
		CHECK_FALSE(timeout.success);
		CHECK(timeout.timedOut);
		CHECK(timeout.error == IOError::TIMEOUT);
		CHECK(GetThreadErrorMode() == mode);
		CHECK(GetErrorMode() == processMode);
	}
}
#endif

///////////////////////////////////////////////////////////////////////////////
//                        VISUALIZER RESULT - SUCCESS                        //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("IOResult::Success - Default", "[Visualizer][Result]") {
	auto result = IOResult::Success();

	REQUIRE(result.success == true);
	REQUIRE(result.exitCode == 0);
	REQUIRE(result.message.empty());
	REQUIRE(result.filename.empty());
}

TEST_CASE("IOResult::Success - With data path", "[Visualizer][Result]") {
	auto result = IOResult::Success("/path/to/data.mml");

	REQUIRE(result.success == true);
	REQUIRE(result.exitCode == 0);
	REQUIRE(result.filename == "/path/to/data.mml");
	REQUIRE(result.message.empty());
}

TEST_CASE("IOResult::Success - With Windows path", "[Visualizer][Result]") {
	auto result = IOResult::Success("C:\\Users\\data\\file.mml");

	REQUIRE(result.success == true);
	REQUIRE(result.filename == "C:\\Users\\data\\file.mml");
}

TEST_CASE("IOResult::Success - With relative path", "[Visualizer][Result]") {
	auto result = IOResult::Success("./results/output.dat");

	REQUIRE(result.success == true);
	REQUIRE(result.filename == "./results/output.dat");
}

///////////////////////////////////////////////////////////////////////////////
//                        VISUALIZER RESULT - FAILURE                        //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("IOResult::Failure - Message only", "[Visualizer][Result]") {
	auto result = IOResult::Failure("Something went wrong");

	REQUIRE(result.success == false);
	REQUIRE(result.exitCode == -1); // Default error code
	REQUIRE(result.message == "Something went wrong");
	REQUIRE(result.filename.empty());
}

TEST_CASE("IOResult::Failure - With data path", "[Visualizer][Result]") {
	auto result = IOResult::Failure("Error occurred", "/path/to/partial.mml");

	REQUIRE(result.success == false);
	REQUIRE(result.exitCode == -1);
	REQUIRE(result.message == "Error occurred");
	REQUIRE(result.filename == "/path/to/partial.mml");
}

TEST_CASE("IOResult::Failure - With exit code", "[Visualizer][Result]") {
	auto result = IOResult::Failure("Process failed", "/path/to/file.mml", 42);

	REQUIRE(result.success == false);
	REQUIRE(result.exitCode == 42);
	REQUIRE(result.message == "Process failed");
	REQUIRE(result.filename == "/path/to/file.mml");
}

TEST_CASE("IOResult::Failure - With negative exit code", "[Visualizer][Result]") {
	auto result = IOResult::Failure("Crashed", "", -9);

	REQUIRE(result.success == false);
	REQUIRE(result.exitCode == -9);
}

TEST_CASE("IOResult::Failure - Empty message", "[Visualizer][Result]") {
	auto result = IOResult::Failure("");

	REQUIRE(result.success == false);
	REQUIRE(result.message.empty());
}

///////////////////////////////////////////////////////////////////////////////
//                      EXPLICIT SUCCESS CHECKS                              //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("IOResult - Explicit success checks", "[Visualizer][Result]") {
	auto success = IOResult::Success();
	auto failure = IOResult::Failure("fail");

	REQUIRE(success.success);
	REQUIRE_FALSE(failure.success);
}

///////////////////////////////////////////////////////////////////////////////
//                       STRUCT MEMBER ACCESS                                //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("IOResult - Aggregate initialization", "[Visualizer][Result]") {
	IOResult result{true, IOError::OK, "", "/data/file.mml"};

	REQUIRE(result.success == true);
	REQUIRE(result.exitCode == 0);
	REQUIRE(result.message == "");
	REQUIRE(result.filename == "/data/file.mml");
}

TEST_CASE("IOResult - Failure aggregate initialization", "[Visualizer][Result]") {
	IOResult result = IOResult::Failure("Command not found", "/usr/bin/viz", 127);

	REQUIRE(result.success == false);
	REQUIRE(result.exitCode == 127);
	REQUIRE(result.message == "Command not found");
	REQUIRE(result.filename == "/usr/bin/viz");
}

TEST_CASE("IOResult - Members are modifiable", "[Visualizer][Result]") {
	IOResult result = IOResult::Success();

	result.success = false;
	result.exitCode = 99;
	result.message = "Modified";
	result.filename = "/new/path.mml";

	REQUIRE(result.success == false);
	REQUIRE(result.exitCode == 99);
	REQUIRE(result.message == "Modified");
	REQUIRE(result.filename == "/new/path.mml");
}

///////////////////////////////////////////////////////////////////////////////
//                       STRING CONTENT TESTS                                //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("IOResult - Long error message", "[Visualizer][Result]") {
	std::string longError = "This is a very long error message that describes "
													"in great detail what went wrong during the visualization "
													"process, including potentially stack traces, file paths, "
													"and other diagnostic information that would be helpful "
													"for debugging purposes.";

	auto result = IOResult::Failure(longError);

	REQUIRE(result.message == longError);
	REQUIRE(result.message.length() > 200);
}

TEST_CASE("IOResult - Unicode in paths", "[Visualizer][Result]") {
	std::string unicodePath = "/path/to/données.mml";
	auto result = IOResult::Success(unicodePath);

	REQUIRE(result.filename == unicodePath);
}

TEST_CASE("IOResult - Special characters in error message", "[Visualizer][Result]") {
	std::string specialMsg = "Error: File \"test.txt\" not found!\n\tCheck path & retry.";
	auto result = IOResult::Failure(specialMsg);

	REQUIRE(result.message == specialMsg);
	REQUIRE(result.message.find('\n') != std::string::npos);
	REQUIRE(result.message.find('\t') != std::string::npos);
}

///////////////////////////////////////////////////////////////////////////////
//                       COPY AND MOVE SEMANTICS                             //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("IOResult - Copy construction", "[Visualizer][Result]") {
	auto original = IOResult::Failure("original error", "/original/path.mml", 42);
	IOResult copy = original;

	REQUIRE(copy.success == original.success);
	REQUIRE(copy.exitCode == original.exitCode);
	REQUIRE(copy.message == original.message);
	REQUIRE(copy.filename == original.filename);

	// Modify copy, original should be unchanged
	copy.message = "modified";
	REQUIRE(original.message == "original error");
}

TEST_CASE("IOResult - Copy assignment", "[Visualizer][Result]") {
	auto original = IOResult::Success("/path/to/file.mml");
	IOResult copy = IOResult::Failure("dummy");

	copy = original;

	REQUIRE(copy.success == true);
	REQUIRE(copy.filename == "/path/to/file.mml");
}

TEST_CASE("IOResult - Move construction", "[Visualizer][Result]") {
	auto original = IOResult::Success("/path/to/file.mml");
	IOResult moved = std::move(original);

	REQUIRE(moved.success == true);
	REQUIRE(moved.filename == "/path/to/file.mml");
}

TEST_CASE("IOResult - Move assignment", "[Visualizer][Result]") {
	auto original = IOResult::Failure("error", "/path.mml", 10);
	IOResult target = IOResult::Success();

	target = std::move(original);

	REQUIRE(target.success == false);
	REQUIRE(target.exitCode == 10);
}

///////////////////////////////////////////////////////////////////////////////
//                       EDGE CASES                                          //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("IOResult - Default construction", "[Visualizer][Result]") {
	IOResult result{};

	REQUIRE(result.success == false); // bool default-initializes to false
	REQUIRE(result.exitCode == 0);		// int default-initializes to 0
	REQUIRE(result.message.empty());
	REQUIRE(result.filename.empty());
}

TEST_CASE("IOResult - Exit code zero with failure", "[Visualizer][Result]") {
	// Rare case: process exits with 0 but we mark as failure
	IOResult result = IOResult::Failure("Logic error", "", 0);

	REQUIRE(result.success == false);
	REQUIRE(result.exitCode == 0);
	REQUIRE_FALSE(result.success);
}

TEST_CASE("IOResult - Non-zero exit with success", "[Visualizer][Result]") {
	// Rare case: marking as success despite non-zero exit (shouldn't happen normally)
	IOResult result = IOResult::Success();
	result.exitCode = 1;

	REQUIRE(result.success == true);
	REQUIRE(result.exitCode == 1);
	REQUIRE(result.success);
}

} // namespace MML::Tests::Tools::VisualizerTests
