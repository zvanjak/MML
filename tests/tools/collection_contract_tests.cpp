#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/Function.h>
#include <mml/tools/Visualizer.h>
#endif
#include <catch2/catch_all.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

using namespace MML;

namespace {
	struct CollectionFile {
		std::filesystem::path path = std::filesystem::temp_directory_path() /
			("mml_collection_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".mml");
		~CollectionFile() {
			std::error_code error;
			std::filesystem::remove(path, error);
			if (error) std::cerr << "Cannot remove collection test file: " << error.message() << '\n';
		}
		void Seed() const {
			std::ofstream out(path);
			out << "SENTINEL";
			out.close();
			REQUIRE_FALSE(out.fail());
		}
		std::string Read() const {
			std::ifstream in(path);
			REQUIRE(in.is_open());
			std::ostringstream out;
			out << in.rdbuf();
			REQUIRE_FALSE(in.bad());
			return out.str();
		}
	};
	template<class Write>
	void RejectStream(Write write) {
		std::ostringstream out;
		out.setf(std::ios::scientific, std::ios::floatfield);
		out.precision(3);
		out << "SENTINEL";
		const auto flags = out.flags();
		const auto result = write(out);
		CHECK_FALSE(result.success);
		CHECK(result.error == IOError::INVALID_PARAMETERS);
		CHECK_FALSE(result.message.empty());
		CHECK(out.str() == "SENTINEL");
		CHECK(out.flags() == flags);
		CHECK(out.precision() == 3);
	}
	struct CountingFunction : IRealFunction {
		mutable int calls = 0;
		Real operator()(Real x) const override { ++calls; return x; }
	};
	template<int N>
	struct CountingCurve : IRealToVectorFunction<N> {
		mutable int calls = 0;
		VectorN<Real, N> operator()(Real x) const override {
			++calls;
			VectorN<Real, N> result;
			for (int i = 0; i < N; ++i) result[i] = x + i;
			return result;
		}
	};
	struct VisualizerFiles {
		std::filesystem::path previousRoot = detail::GetVisualizerRootState();
		std::string base = "mml_collection_visualizer_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count());
		std::vector<std::filesystem::path> paths;
		VisualizerFiles() {
			SetVisualizerRoot(std::filesystem::temp_directory_path() / (base + "_missing_bundle"));
			const auto directory = std::filesystem::path(GetResultFilesPath());
			paths.push_back(directory / base);
			for (int i = 1; i <= 3; ++i) paths.push_back(directory / (base + "_" + std::to_string(i) + ".mml"));
		}
		~VisualizerFiles() {
			SetVisualizerRoot(previousRoot);
			for (const auto& path : paths) {
				std::error_code error;
				std::filesystem::remove(path, error);
				if (error) std::cerr << "Cannot remove visualizer collection test file: " << error.message() << '\n';
			}
		}
		void Seed() const {
			std::filesystem::create_directories(paths.front().parent_path());
			for (const auto& path : paths) {
				std::ofstream out(path);
				out << "SENTINEL";
				out.close();
				REQUIRE_FALSE(out.fail());
			}
		}
		std::string Read(std::size_t index) const {
			std::ifstream in(paths[index]);
			REQUIRE(in.is_open());
			std::ostringstream out;
			out << in.rdbuf();
			return out.str();
		}
		template<class Call>
		void Reject(Call call, bool seed) {
			if (seed) Seed();
			const auto result = call(base);
			CHECK_FALSE(result.success);
			CHECK(result.error == IOError::INVALID_PARAMETERS);
			CHECK_FALSE(result.message.empty());
			CHECK(result.filename == base);
			for (std::size_t i = 0; i < paths.size(); ++i) {
				if (seed) CHECK(Read(i) == "SENTINEL");
				else CHECK_FALSE(std::filesystem::exists(paths[i]));
			}
		}
	};
}

TEST_CASE("Collection headers reject negative counts and ignored extra entries", "[serializer][Collections]") {
	for (int count : {-1, 1}) {
		const std::vector<std::string> labels = count < 0 ? std::vector<std::string>{} : std::vector<std::string>{"a", "b"};
		const std::vector<Real> radii(labels.size(), Real(0.5));
		RejectStream([&](auto& out) { return Serializer::WriteRealMultiFuncHeader(out, "test", count, labels, 0, 1, 2); });
		RejectStream([&](auto& out) { return Serializer::WriteParticleSimulation2DHeader(out, count, 10, 20, labels, radii, 2); });
		RejectStream([&](auto& out) { return Serializer::WriteParticleSimulation3DHeader(out, count, 10, 20, 30, labels, radii, 2); });
	}
}

TEST_CASE("Invalid multi-function collections do not truncate existing files", "[serializer][Collections]") {
	RealFunction function([](Real x) { return x; });
	for (const auto& labels : {std::vector<std::string>{}, std::vector<std::string>{"a", "b"}}) {
		CollectionFile file;
		file.Seed();
		const auto result = Serializer::SaveRealMultiFunc(std::vector<IRealFunction*>{&function}, "test", labels, 0, 1, 2, file.path.string());
		CHECK(result.error == IOError::INVALID_PARAMETERS);
		CHECK(result.filename == file.path.string());
		CHECK(file.Read() == "SENTINEL");
	}
}

TEST_CASE("Null multi-function entries reject the entire request before evaluation or writes", "[serializer][Collections]") {
		CountingFunction function;
		for (const auto& funcs : {std::vector<IRealFunction*>{nullptr}, std::vector<IRealFunction*>{&function, nullptr},
			std::vector<IRealFunction*>{nullptr, &function}}) {
			const std::vector<std::string> labels(funcs.size(), "test");
			RejectStream([&](auto& out) { return Serializer::SaveRealMultiFunc(out, funcs, "test", labels, 0, 1, 3); });
			CollectionFile file;
			file.Seed();
			const auto result = Serializer::SaveRealMultiFunc(funcs, "test", labels, 0, 1, 3, file.path.string());
			CHECK(result.error == IOError::INVALID_PARAMETERS);
			CHECK(result.filename == file.path.string());
			CHECK(file.Read() == "SENTINEL");
			CHECK(function.calls == 0);
			VisualizerFiles visualizer;
			visualizer.Reject([&](const auto& name) {
				return Visualizer::VisualizeMultiRealFunction(funcs, "test", labels, 0, 1, 3, name);
			}, true);
		}
	}

	TEST_CASE("Collection headers reject short parallel arrays before changing streams", "[serializer][Collections]") {
		RejectStream([](auto& out) { return Serializer::WriteRealMultiFuncHeader(out, "test", 2, {"one"}, 0, 1, 2); });
		for (const auto& colors : {std::vector<std::string>{}, std::vector<std::string>{"a"}, std::vector<std::string>{"a", "b"}}) {
			for (const auto& radii : {std::vector<Real>{}, std::vector<Real>{1}, std::vector<Real>{1, 2}}) {
				if (colors.size() == 1 && radii.size() == 1) continue;
				RejectStream([&](auto& out) { return Serializer::WriteParticleSimulation2DHeader(out, 1, 10, 20, colors, radii, 2); });
				RejectStream([&](auto& out) { return Serializer::WriteParticleSimulation3DHeader(out, 1, 10, 20, 30, colors, radii, 2); });
			}
		}
	}

	TEST_CASE("Zero-count low-level collection headers remain valid with empty collections", "[serializer][Collections]") {
		std::ostringstream multi, particles2D, particles3D;
		REQUIRE(Serializer::WriteRealMultiFuncHeader(multi, "zero", 0, {}, 0, 1, 2).success);
		CHECK(multi.str() == "MML_MULTI_REAL_FUNCTION\nVERSION: 1\nzero\n0\nx1: 0\nx2: 1\nNumPoints: 2\n");
		REQUIRE(Serializer::WriteParticleSimulation2DHeader(particles2D, 0, 10, 20, {}, {}, 0).success);
		CHECK(particles2D.str() == "MML_PARTICLE_SIMULATION_DATA_2D\nVERSION: 1\nWidth: 10\nHeight: 20\nNumBalls: 0\nNumSteps: 0\n");
		REQUIRE(Serializer::WriteParticleSimulation3DHeader(particles3D, 0, 10, 20, 30, {}, {}, 0).success);
		CHECK(particles3D.str() == "MML_PARTICLE_SIMULATION_DATA_3D\nVERSION: 1\nWidth: 10\nHeight: 20\nDepth: 30\nNumBalls: 0\nNumSteps: 0\n");
		RejectStream([](auto& out) { return Serializer::WriteRealMultiFuncHeader(out, "zero", 0, {"extra"}, 0, 1, 2); });
	}

	TEST_CASE("Collection count validation covers representability limits without allocating collections", "[serializer][Collections]") {
		const auto maximum = static_cast<std::size_t>(std::numeric_limits<int>::max());
		CHECK(Serializer::detail::ValidateCollectionCount(maximum, "test").success);
		for (std::size_t size : {std::size_t(0), maximum + 1, std::numeric_limits<std::size_t>::max()}) {
			const auto result = Serializer::detail::ValidateCollectionCount(size, "test");
			CHECK(result.error == IOError::INVALID_PARAMETERS);
			CHECK_FALSE(result.message.empty());
		}
		CHECK(Serializer::detail::ValidateHeaderCollectionCount(std::numeric_limits<int>::max(), maximum, "test").success);
		CHECK(Serializer::detail::ValidateHeaderCollectionCount(std::numeric_limits<int>::max(), maximum + 1, "test").error == IOError::INVALID_PARAMETERS);
	}

	TEST_CASE("Separate-function visualization validates all legends before creating any output", "[Visualizer][Collections]") {
		const LinearInterpRealFunc linear(Vector<Real>{0, 1}, Vector<Real>{0, 1});
		for (const auto& labels : {std::vector<std::string>{}, std::vector<std::string>{"one"},
			std::vector<std::string>{"one", "two", "extra"}}) {
			for (bool seed : {false, true}) {
				VisualizerFiles files;
				files.Reject([&](const auto& name) {
					return Visualizer::VisualizeMultiRealFunctionSeparately({linear, linear}, "test", labels, 0, 1, 3, name);
				}, seed);
			}
		}
	}

	TEMPLATE_TEST_CASE("Interpolant collections share exact legend checks in serializer and visualizer",
		"[serializer][Visualizer][Collections]", LinearInterpRealFunc, PolynomInterpRealFunc, SplineInterpRealFunc) {
		const TestType function = [] {
			const Vector<Real> nodes{0, Real(0.5), 1};
			if constexpr (std::is_same_v<TestType, PolynomInterpRealFunc>)
				return TestType(nodes, nodes, 2);
			else
				return TestType(nodes, nodes);
		}();
		for (const auto& labels : {std::vector<std::string>{}, std::vector<std::string>{"one", "extra"}}) {
			CollectionFile file;
			file.Seed();
			const auto result = Serializer::SaveRealMultiFunc(std::vector<TestType>{function}, "test", labels, 0, 1, 3, file.path.string());
			CHECK(result.error == IOError::INVALID_PARAMETERS);
			CHECK(file.Read() == "SENTINEL");
			VisualizerFiles files;
			files.Reject([&](const auto& name) {
				return Visualizer::VisualizeMultiRealFunction(std::vector<TestType>{function}, "test", labels, 0, 1, 3, name);
			}, false);
		}
		CollectionFile valid;
		REQUIRE(Serializer::SaveRealMultiFunc(std::vector<TestType>{function}, "valid", {"label\nline"}, 0, 1, 3, valid.path.string()).success);
		CHECK(valid.Read().find("valid\n1\nlabel line\n") != std::string::npos);
		VisualizerFiles files;
		CHECK(Visualizer::VisualizeMultiRealFunction(std::vector<TestType>{function}, "valid", {"label"}, 0, 1, 3, files.base).error == IOError::FILE_NOT_OPENED);
		CHECK(files.Read(0).find("valid\n1\nlabel\n") != std::string::npos);
	}

	TEST_CASE("Curve collections reject every null entry before evaluating earlier curves or creating files", "[Visualizer][Collections]") {
		CountingCurve<2> curve2D;
		CountingCurve<3> curve3D;
		for (bool seed : {false, true}) {
			VisualizerFiles first, second, third;
			first.Reject([&](const auto& name) {
				return Visualizer::VisualizeMultiParamCurve2D({&curve2D, nullptr}, std::string("test"), 0, 1, 3, name);
			}, seed);
			second.Reject([&](const auto& name) {
				return Visualizer::VisualizeMultiParamCurve2D({&curve2D, nullptr}, std::vector<std::string>{"a", "b"}, 0, 1, 3, name);
			}, seed);
			third.Reject([&](const auto& name) {
				return Visualizer::VisualizeMultiParamCurve3D({&curve3D, nullptr}, "test", 0, 1, 3, name);
			}, seed);
		}
		CHECK(curve2D.calls == 0);
		CHECK(curve3D.calls == 0);
	}

	TEST_CASE("Curve legends must match the complete collection rather than a serialized prefix", "[Visualizer][Collections]") {
		CountingCurve<2> curve;
		for (const auto& labels : {std::vector<std::string>{}, std::vector<std::string>{"one"},
			std::vector<std::string>{"one", "two", "extra"}}) {
			VisualizerFiles files;
			files.Reject([&](const auto& name) {
				return Visualizer::VisualizeMultiParamCurve2D({&curve, &curve}, labels, 0, 1, 3, name);
			}, true);
		}
		CHECK(curve.calls == 0);
	}

	TEST_CASE("High-level empty collections are rejected while low-level zero headers stay supported", "[Visualizer][serializer][Collections]") {
		RejectStream([](auto& out) { return Serializer::SaveRealMultiFunc(out, {}, "empty", {}, 0, 1, 3); });
		VisualizerFiles functions, separate, curves2D, curves3D;
		functions.Reject([](const auto& name) {
			return Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{}, "test", {}, 0, 1, 3, name);
		}, false);
		separate.Reject([](const auto& name) {
			return Visualizer::VisualizeMultiRealFunctionSeparately({}, "test", {}, 0, 1, 3, name);
		}, false);
		curves2D.Reject([](const auto& name) {
			return Visualizer::VisualizeMultiParamCurve2D(std::vector<IRealToVectorFunction<2>*>{}, std::string("test"), 0, 1, 3, name);
		}, false);
		curves3D.Reject([](const auto& name) {
			return Visualizer::VisualizeMultiParamCurve3D(std::vector<IRealToVectorFunction<3>*>{}, "test", 0, 1, 3, name);
		}, false);
	}

	TEST_CASE("Valid multi-function encoding and callback-error behavior are preserved", "[serializer][Collections]") {
		CountingFunction first, second;
		std::ostringstream out;
		REQUIRE(Serializer::SaveRealMultiFunc(out, {&first, &second}, "test", {"a", "b"}, 0, 1, 3).success);
		CHECK(out.str() == "MML_MULTI_REAL_FUNCTION\nVERSION: 1\ntest\n2\na\nb\nx1: 0\nx2: 1\nNumPoints: 3\n0 0 0 \n0.5 0.5 0.5 \n1 1 1 \n");
		CHECK(first.calls == 3);
		CHECK(second.calls == 3);
		RealFunction broken([](Real) -> Real { throw std::runtime_error("callback failed"); });
		std::ostringstream failed;
		CHECK(Serializer::SaveRealMultiFunc(failed, {&broken}, "test", {"a"}, 0, 1, 3).error == IOError::WRITE_FAILED);
	}

	TEST_CASE("Valid visualizer collections serialize all outputs before a safely missing executable is reported", "[Visualizer][Collections]") {
		CountingCurve<2> curve2D;
		CountingCurve<3> curve3D;
		{
			VisualizerFiles files;
			const auto result = Visualizer::VisualizeMultiParamCurve2D({&curve2D, &curve2D}, std::vector<std::string>{"first", "second"}, 0, 1, 3, files.base);
			CHECK(result.error == IOError::FILE_NOT_OPENED);
			CHECK(files.Read(1).find("first\n") != std::string::npos);
			CHECK(files.Read(2).find("second\n") != std::string::npos);
			CHECK(curve2D.calls == 6);
		}
		{
			VisualizerFiles files;
			const auto result = Visualizer::VisualizeMultiParamCurve3D({&curve3D}, "third", 0, 1, 3, files.base);
			CHECK(result.error == IOError::FILE_NOT_OPENED);
			CHECK(files.Read(1).find("third\n") != std::string::npos);
			CHECK(curve3D.calls == 3);
		}
		{
			VisualizerFiles files;
			const LinearInterpRealFunc function(Vector<Real>{0, 1}, Vector<Real>{0, 1});
			const auto result = Visualizer::VisualizeMultiRealFunctionSeparately({function, function}, "unused", {"first", "second"}, -1, 2, 3, files.base);
			CHECK(result.error == IOError::FILE_NOT_OPENED);
			CHECK(files.Read(1).find("first\n") != std::string::npos);
			CHECK(files.Read(2).find("second\n") != std::string::npos);
			CHECK(files.Read(1).find("x1: 0\nx2: 1\nNumPoints: 3\n") != std::string::npos);
		}
		{
			VisualizerFiles files;
			CHECK(Visualizer::VisualizeMultiParamCurve2D({&curve2D}, std::string("common"), 0, 1, 3, files.base).error == IOError::FILE_NOT_OPENED);
			CHECK(files.Read(1).find("common\n") != std::string::npos);
			CHECK(curve2D.calls == 9);
		}
		{
			VisualizerFiles files;
			CountingFunction function;
			CHECK(Visualizer::VisualizeMultiRealFunction(std::vector<IRealFunction*>{&function}, "pointer", {"label"}, 0, 1, 3, files.base).error == IOError::FILE_NOT_OPENED);
			CHECK(files.Read(0).find("pointer\n1\nlabel\n") != std::string::npos);
			CHECK(function.calls == 3);
		}
	}

TEST_CASE("Invalid complete multi-function requests leave both new and existing destinations untouched", "[serializer][Collections]") {
	CountingFunction function;
	for (int scenario = 0; scenario < 4; ++scenario) {
		const std::vector<IRealFunction*> funcs = scenario == 0 ? std::vector<IRealFunction*>{} : std::vector<IRealFunction*>{&function};
		const std::vector<std::string> labels(funcs.size(), "label");
		const Real end = scenario == 1 ? Real(0) : Real(1);
		const int points = scenario == 2 ? 1 : 3;
		const auto actualLabels = scenario == 3 ? std::vector<std::string>{"extra", "extra"} : labels;
		RejectStream([&](auto& out) { return Serializer::SaveRealMultiFunc(out, funcs, "test", actualLabels, 0, end, points); });
		for (bool seed : {false, true}) {
			CollectionFile file;
			if (seed) file.Seed();
			const auto result = Serializer::SaveRealMultiFunc(funcs, "test", actualLabels, 0, end, points, file.path.string());
			CHECK(result.error == IOError::INVALID_PARAMETERS);
			CHECK_FALSE(result.message.empty());
			CHECK(result.filename == file.path.string());
			if (seed) CHECK(file.Read() == "SENTINEL");
			else CHECK_FALSE(std::filesystem::exists(file.path));
		}
	}
	CHECK(function.calls == 0);
}

TEST_CASE("Precomputed curves retain intentional sample-prefix serialization", "[serializer][Collections]") {
	const std::vector<VectorN<Real, 2>> points{{1, 2}, {3, 4}, {999, 999}};
	CollectionFile file;
	REQUIRE(Serializer::SaveAsParamCurve(points, "MML_PARAMETRIC_CURVE_CARTESIAN_2D", "prefix", 0, 1, 2, file.path.string()).success);
	const auto bytes = file.Read();
	CHECK(bytes.find("0 1 2 \n1 3 4 \n") != std::string::npos);
	CHECK(bytes.find("999") == std::string::npos);
	const auto rejected = Serializer::SaveAsParamCurve(std::vector<VectorN<Real, 2>>{{1, 2}}, "unused", "short", 0, 1, 2, file.path.string());
	CHECK(rejected.error == IOError::INVALID_PARAMETERS);
	CHECK(file.Read() == bytes);
}
