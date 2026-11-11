#include <catch2/catch_all.hpp>

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/base/Function.h>
#include <mml/tools/Serializer.h>
#include <mml/tools/serializer/SerializerFieldLines.h>
#endif

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace MML;

namespace {
	struct OutputFile {
		std::filesystem::path path = std::filesystem::temp_directory_path() /
			("mml_serializer_preflight_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".mml");
		~OutputFile() {
			std::error_code error;
			std::filesystem::remove(path, error);
			if (error) std::cerr << "Cannot remove preflight test file: " << error.message() << '\n';
		}
		void Seed() const {
			std::ofstream out(path, std::ios::binary);
			out << "SENTINEL";
			out.close();
			REQUIRE_FALSE(out.fail());
		}
		std::string Read() const {
			std::ifstream in(path, std::ios::binary);
			REQUIRE(in.is_open());
			std::ostringstream out;
			out << in.rdbuf();
			REQUIRE_FALSE(in.bad());
			return out.str();
		}
	};

	void CheckInvalid(const IOResult& result, const std::string& path) {
		CHECK_FALSE(result.success);
		CHECK(result.error == IOError::INVALID_PARAMETERS);
		CHECK(result.filename == path);
		CHECK_FALSE(result.message.empty());
	}

	template<int N>
	struct Curve : IRealToVectorFunction<N> {
		VectorN<Real, N> operator()(Real t) const override {
			VectorN<Real, N> result;
			result[0] = t;
			return result;
		}
	};

	struct Surface : IParametricSurface<3>, IParametricSurfaceRect<3> {
		Real getMinU() const override { return 0; }
		Real getMaxU() const override { return 1; }
		Real getMinW() const override { return 0; }
		Real getMaxW() const override { return 1; }
		Real getMinW(Real) const override { return 0; }
		Real getMaxW(Real) const override { return 1; }
		VectorN<Real, 3> operator()(Real u, Real w) const override { return {u, w, 0}; }
		VectorN<Real, 3> operator()(const VectorN<Real, 2>& p) const override { return {p[0], p[1], 0}; }
	};

	struct Grid {
		const Grid& grid() const { return *this; }
		const Grid& domain() const { return *this; }
		int numNodesX() const { return 2; }
		int numNodesY() const { return 2; }
		Real xMin() const { return 0; }
		Real xMax() const { return 1; }
		Real yMin() const { return 0; }
		Real yMax() const { return 1; }
		Real x(int i) const { return Real(i); }
		Real y(int i) const { return Real(i); }
		Real operator()(int i, int j) const { return Real(i + j); }
	};
}

TEST_CASE("Invalid real function requests never open output files", "[serializer][serializer-preflight]")
{
	int calls = 0;
	RealFunctionFromStdFunc function([&](Real x) { ++calls; return x; });
	using Writer = std::function<IOResult(const std::string&)>;
	std::vector<Writer> invalid;
	for (int count : {-1, 0, 1}) {
		invalid.push_back([&, count](const auto& p) { return Serializer::SaveRealFunc(function, "test", 0, 1, count, p); });
		invalid.push_back([&, count](const auto& p) { return Serializer::SaveRealFuncEquallySpaced(function, "test", 0, 1, count, p); });
	}
	for (Real end : {Real(0), Real(-1)}) {
		invalid.push_back([&, end](const auto& p) { return Serializer::SaveRealFunc(function, "test", 0, end, 2, p); });
		invalid.push_back([&, end](const auto& p) { return Serializer::SaveRealFuncEquallySpaced(function, "test", 0, end, 2, p); });
	}
	for (const Vector<Real>& points : {Vector<Real>{}, Vector<Real>{0}}) {
		invalid.push_back([&, points](const auto& p) { return Serializer::SaveRealFunc(function, "test", points, p); });
	}
	for (bool seed : {false, true}) {
		for (const auto& write : invalid) {
			OutputFile file;
			if (seed) file.Seed();
			CheckInvalid(write(file.path.string()), file.path.string());
			CHECK(calls == 0);
			if (seed) CHECK(file.Read() == "SENTINEL");
			else CHECK_FALSE(std::filesystem::exists(file.path));
		}
	}
}

TEST_CASE("Function stream preflight rejects before formatting or evaluation", "[serializer][serializer-preflight]")
{
	int calls = 0;
	RealFunctionFromStdFunc function([&](Real x) { ++calls; return x; });
	using Writer = std::function<IOResult(std::ostream&)>;
	const std::vector<Writer> invalid{
		[&](auto& out) { return Serializer::SaveRealFunc(out, function, "test", 0, 1, 1); },
		[&](auto& out) { return Serializer::SaveRealFunc(out, function, "test", 1, 0, 2); },
		[&](auto& out) { return Serializer::SaveRealFunc(out, function, "test", Vector<Real>{0}); },
		[&](auto& out) { return Serializer::SaveRealFuncEquallySpaced(out, function, "test", 0, 1, 1); },
		[&](auto& out) { return Serializer::SaveRealFuncEquallySpaced(out, function, "test", 1, 0, 2); }};
	for (const auto& write : invalid) {
		std::ostringstream out;
		out << "SENTINEL";
		out << std::scientific;
		out.precision(3);
		const auto flags = out.flags();
		CheckInvalid(write(out), "");
		CHECK(out.str() == "SENTINEL");
		CHECK(out.flags() == flags);
		CHECK(out.precision() == 3);
		CHECK(calls == 0);
	}
}

TEST_CASE("Every serializer file writer rejects the full embedded NUL filename", "[serializer][serializer-preflight]")
{
	int evaluations = 0;
	RealFunctionFromStdFunc function([&](Real x) { ++evaluations; return x; });
	ScalarFunctionFromStdFunc<2> scalar2([](const VectorN<Real, 2>& p) { return p[0]; });
	ScalarFunctionFromStdFunc<3> scalar3([](const VectorN<Real, 3>& p) { return p[0]; });
	VectorFunctionFromStdFunc<2> vector2([](const VectorN<Real, 2>& p) { return p; });
	VectorFunctionFromStdFunc<3> vector3([](const VectorN<Real, 3>& p) { return p; });
	const Curve<2> curve2;
	const Curve<3> curve3;
	const Surface surface;
	const Grid grid;
	const Vector<Real> parameters{0, 1};
#ifndef MML_SINGLE_HEADER_INCLUDE
	const std::vector<VectorN<Real, 2>> points2{{0, 0}, {1, 1}};
#endif
	const std::vector<VectorN<Real, 3>> points3{{0, 0, 0}, {1, 1, 1}};
	const std::vector<IRealFunction*> functions{&function};
	const std::vector<LinearInterpRealFunc> interpolants{LinearInterpRealFunc(parameters, parameters)};
	ODESystemSolution solution(0, 1, 3, 2);
	solution.fillValues(0, 0, Vector<Real>{0, 0, 0});
	solution.fillValues(1, 1, Vector<Real>{1, 1, 1});
	struct Writer {
		const char* name;
		std::function<IOResult(const std::string&)> write;
	};
	const std::vector<Writer> writers{
		{"real", [&](const auto& p) { return Serializer::SaveRealFunc(function, "test", 0, 1, 2, p); }},
		{"real-points", [&](const auto& p) { return Serializer::SaveRealFunc(function, "test", parameters, p); }},
		{"equally-spaced", [&](const auto& p) { return Serializer::SaveRealFuncEquallySpaced(function, "test", 0, 1, 2, p); }},
		{"multi", [&](const auto& p) { return Serializer::SaveRealMultiFunc(functions, "test", {"f"}, 0, 1, 2, p); }},
		{"interpolants", [&](const auto& p) { return Serializer::SaveRealMultiFunc(interpolants, "test", {"f"}, 0, 1, 2, p); }},
		{"curve2", [&](const auto& p) { return Serializer::SaveParamCurveCartesian2DResult(curve2, "test", 0, 1, 2, p); }},
		{"curve3", [&](const auto& p) { return Serializer::SaveParamCurveCartesian3DResult(curve3, "test", 0, 1, 2, p); }},
		{"curve-points", [&](const auto& p) { return Serializer::SaveParamCurve<3>(curve3, "CURVE", "test", parameters, p); }},
		{"raw-curve", [&](const auto& p) { return Serializer::SaveAsParamCurve<3>(points3, "CURVE", "test", 0, 1, 2, p); }},
		{"raw-curve-points", [&](const auto& p) { return Serializer::SaveAsParamCurve<3>(points3, "CURVE", "test", parameters, p); }},
		{"xy-curve", [&](const auto& p) { return Serializer::SaveAsParamCurve2D(parameters, parameters, "test", p); }},
		{"surface", [&](const auto& p) { return Serializer::SaveParametricSurface(
			static_cast<const IParametricSurface<3>&>(surface), "test", 0, 1, 2, 0, 1, 2, p); }},
		{"surface-vector", [&](const auto& p) { return Serializer::SaveParametricSurface(
			static_cast<const IVectorFunctionNM<2, 3>&>(static_cast<const IParametricSurface<3>&>(surface)), "test", 0, 1, 2, 0, 1, 2, p); }},
		{"surface-rect", [&](const auto& p) { return Serializer::SaveParametricSurface(
			static_cast<const IParametricSurfaceRect<3>&>(surface), "test", 2, 2, p); }},
		{"scalar2", [&](const auto& p) { return Serializer::SaveScalarFunc2DCartesian(scalar2, "test", 0, 1, 2, 0, 1, 2, p); }},
		{"scalar3", [&](const auto& p) { return Serializer::SaveScalarFunc3DCartesian(scalar3, "test", 0, 1, 2, 0, 1, 2, 0, 1, 2, p); }},
		{"grid", [&](const auto& p) { return Serializer::SaveGridFunction2D(grid, "test", p); }},
		{"vector2", [&](const auto& p) { return Serializer::SaveVectorFunc2DCartesian(vector2, "test", 0, 1, 2, 0, 1, 2, p); }},
		{"vector2-filter", [&](const auto& p) { return Serializer::SaveVectorFunc2DCartesian(vector2, "test", 0, 1, 2, 0, 1, 2, p, 10); }},
		{"vector3", [&](const auto& p) { return Serializer::SaveVectorFunc3DCartesian(vector3, "test", 0, 1, 2, 0, 1, 2, 0, 1, 2, p); }},
		{"vector3-filter", [&](const auto& p) { return Serializer::SaveVectorFunc3DCartesian(vector3, "test", 0, 1, 2, 0, 1, 2, 0, 1, 2, p, 10); }},
		{"spherical", [&](const auto& p) { return Serializer::SaveVectorFuncSpherical(vector3, "test", 1, 2, 2, 0, 1, 2, 0, 1, 2, p); }},
		{"spherical-filter", [&](const auto& p) { return Serializer::SaveVectorFuncSpherical(vector3, "test", 1, 2, 2, 0, 1, 2, 0, 1, 2, p, 10); }},
		{"ode-component", [&](const auto& p) { return Serializer::SaveODESolutionComponentAsFunc(solution, 0, "test", p); }},
		{"ode-multi", [&](const auto& p) { return Serializer::SaveODESolutionAsMultiFunc(solution, "test", {"x", "y", "z"}, p); }},
		{"ode-curve2", [&](const auto& p) { return Serializer::SaveODESolAsParametricCurve2D(solution, p, 0, 1, "test"); }},
		{"ode-curve3", [&](const auto& p) { return Serializer::SaveODESolAsParametricCurve3D(solution, p, 0, 1, 2, "test"); }},
		{"particles2", [&](const auto& p) { return Serializer::SaveParticleSimulation2D(
			p, 1, 10, 10, {{{0, 0}, {1, 1}}}, {"red"}, {1}, 1); }},
		{"particles3", [&](const auto& p) { return Serializer::SaveParticleSimulation3D(
			p, 1, 10, 10, 10, {{{0, 0, 0}, {1, 1, 1}}}, {"red"}, {1}, 1); }},
#ifndef MML_SINGLE_HEADER_INCLUDE
		{"field-lines2", [&](const auto& p) { return Serializer::SaveFieldLines2D({points2}, "test", p); }},
		{"field-lines3", [&](const auto& p) { return Serializer::SaveFieldLines3D({points3}, "test", p); }},
#endif
	};
	for (const auto& writer : writers) {
		CAPTURE(writer.name);
		CheckInvalid(writer.write(""), "");
		for (bool seed : {false, true}) {
			OutputFile file;
			if (seed) file.Seed();
			const std::string path = file.path.string() + std::string(1, '\0') + "_ignored_suffix.mml";
			CheckInvalid(writer.write(path), path);
			if (seed) CHECK(file.Read() == "SENTINEL");
			else CHECK_FALSE(std::filesystem::exists(file.path));
		}
	}
	CHECK(evaluations == 0);
}

TEST_CASE("Valid function file exports retain stream bytes and explicit point order", "[serializer][serializer-preflight]")
{
	RealFunctionFromStdFunc function([](Real x) { return x * x; });
	const Vector<Real> points{1, 0, 1};
	for (int kind = 0; kind < 3; ++kind) {
		OutputFile file;
		file.Seed();
		std::ostringstream expected;
		IOResult saved, streamed;
		if (kind == 0) {
			streamed = Serializer::SaveRealFunc(expected, function, "valid", 0, 1, 3);
			saved = Serializer::SaveRealFunc(function, "valid", 0, 1, 3, file.path.string());
		} else if (kind == 1) {
			streamed = Serializer::SaveRealFuncEquallySpaced(expected, function, "valid", 0, 1, 3);
			saved = Serializer::SaveRealFuncEquallySpaced(function, "valid", 0, 1, 3, file.path.string());
		} else {
			streamed = Serializer::SaveRealFunc(expected, function, "valid", points);
			saved = Serializer::SaveRealFunc(function, "valid", points, file.path.string());
		}
		REQUIRE(streamed.success);
		REQUIRE(saved.success);
		CHECK(saved.filename == file.path.string());
		std::string actual = file.Read();
		actual.erase(std::remove(actual.begin(), actual.end(), '\r'), actual.end());
		CHECK(actual == expected.str());
	}
}

TEST_CASE("Valid serializer requests still report open and write finalization failures", "[serializer][serializer-preflight]")
{
	RealFunctionFromStdFunc throwing([](Real) -> Real { throw std::runtime_error("evaluation failed"); });
	OutputFile file;
	file.Seed();
	const auto failed = Serializer::SaveRealFunc(throwing, "valid request", 0, 1, 2, file.path.string());
	CHECK_FALSE(failed.success);
	CHECK(failed.error == IOError::WRITE_FAILED);
	CHECK(failed.filename == file.path.string());
	CHECK_FALSE(failed.message.empty());
	std::ofstream closed(file.path);
	closed.close();
	const auto finish = Serializer::detail::FinishFile(closed, file.path.string());
	CHECK_FALSE(finish.success);
	CHECK(finish.error == IOError::WRITE_FAILED);
	RealFunctionFromStdFunc function([](Real x) { return x; });
	const auto missingParent = file.path / "output.mml";
	const auto unavailable = Serializer::SaveRealFunc(function, "test", 0, 1, 2, missingParent.string());
	CHECK(unavailable.error == IOError::FILE_NOT_OPENED);
	CHECK(unavailable.filename == missingParent.string());
}
