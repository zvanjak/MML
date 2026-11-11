#include <catch2/catch_all.hpp>

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/tools/Serializer.h>
#endif

#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace MML;

namespace {
	struct MetadataFile {
		std::filesystem::path path = std::filesystem::temp_directory_path() /
			("mml_serializer_metadata_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".mml");
		~MetadataFile() {
			std::error_code error;
			std::filesystem::remove(path, error);
			if (error) std::cerr << "Cannot remove metadata test file: " << error.message() << '\n';
		}
		void Seed() const {
			std::ofstream out(path, std::ios::binary);
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

	std::vector<std::string> Lines(const std::string& text) {
		std::istringstream in(text);
		std::vector<std::string> result;
		for (std::string line; std::getline(in, line);)
			result.push_back(line);
		return result;
	}

	template<int N>
	struct MetadataCurve : IRealToVectorFunction<N> {
		VectorN<Real, N> operator()(Real t) const override {
			VectorN<Real, N> result;
			for (int i = 0; i < N; ++i) result[i] = t;
			return result;
		}
	};

	// Qt and WPF consume a mandatory title slot before the three curve header fields.
	void CheckCurveReaderContract(const std::string& text, int dimension, const std::string& title) {
		const auto lines = Lines(text);
		REQUIRE(lines.size() == 8);
		CHECK(lines[0] == (dimension == 2 ? SerializeFormatType::PARAMETRIC_CURVE_CARTESIAN_2D :
			SerializeFormatType::PARAMETRIC_CURVE_CARTESIAN_3D));
		CHECK(lines[1] == "VERSION: 1");
		CHECK(lines[2] == title);
		CHECK(lines[3] == "t1: 0");
		CHECK(lines[4] == "t2: 1");
		CHECK(lines[5] == "NumPoints: 2");
		for (int row = 0; row < 2; ++row) {
			std::istringstream data(lines[6 + row]);
			for (int column = 0; column <= dimension; ++column) {
				Real value;
				REQUIRE(static_cast<bool>(data >> value));
				CHECK(value == Real(row));
			}
			std::string extra;
			CHECK_FALSE(static_cast<bool>(data >> extra));
		}
	}

	const std::vector<std::string> invalidColors{
		"", " ", "red blue", " red", "red ", "red\tblue", "red\nNumSteps: 999",
		"red\rNumSteps: 999", "red\r\nNumSteps: 999", "red\vblue", "red\fblue",
		std::string("red\0blue", 8), std::string("red") + '\x1b', std::string("red") + '\x7f'};

	IOResult ParticleHeader(std::ostream& out, int dimension, const std::vector<std::string>& colors) {
		const std::vector<Real> radii(colors.size(), Real(1));
		const int count = static_cast<int>(colors.size());
		if (dimension == 2)
			return Serializer::WriteParticleSimulation2DHeader(out, count, 10, 20, colors, radii, 1);
		return Serializer::WriteParticleSimulation3DHeader(out, count, 10, 20, 30, colors, radii, 1);
	}

	IOResult ParticleFile(const std::string& path, int dimension, const std::vector<std::string>& colors) {
		if (dimension == 2)
			return Serializer::SaveParticleSimulation2D(path, 2, 10, 20, {{{0, 0}}, {{1, 1}}}, colors, {1, 2}, 1);
		return Serializer::SaveParticleSimulation3D(path, 2, 10, 20, 30, {{{0, 0, 0}}, {{1, 1, 1}}}, colors, {1, 2}, 1);
	}
}

TEST_CASE("Curve headers retain the mandatory empty title slot", "[serializer][serializer-metadata]")
{
	for (const auto* type : {SerializeFormatType::PARAMETRIC_CURVE_CARTESIAN_2D,
		SerializeFormatType::PARAMETRIC_CURVE_CARTESIAN_3D}) {
		for (const auto& title : {std::string{}, std::string("ordinary title"), std::string("first\r\nsecond")}) {
			std::ostringstream out;
			REQUIRE(Serializer::WriteParamCurveHeader(out, type, title, 0, 1, 2).success);
			const std::string expectedTitle = title == "first\r\nsecond" ? "first  second" : title;
			CHECK(out.str() == std::string(type) + "\nVERSION: 1\n" + expectedTitle + "\nt1: 0\nt2: 1\nNumPoints: 2\n");
		}
	}
}

TEST_CASE("Curve and ODE exporters preserve fixed reader header positions", "[serializer][serializer-metadata]")
{
	const MetadataCurve<2> curve2;
	const MetadataCurve<3> curve3;
	const Vector<Real> parameters{0, 1};
	const std::vector<VectorN<Real, 2>> points2{{0, 0}, {1, 1}};
	const std::vector<VectorN<Real, 3>> points3{{0, 0, 0}, {1, 1, 1}};
	ODESystemSolution solution(0, 1, 3, 2);
	solution.fillValues(0, 0, Vector<Real>{0, 0, 0});
	solution.fillValues(1, 1, Vector<Real>{1, 1, 1});
	struct Writer {
		const char* name;
		int dimension;
		std::function<IOResult(const std::string&, const std::string&)> write;
	};
	const std::vector<Writer> writers{
		{"curve2", 2, [&](const auto& p, const auto& title) { return Serializer::SaveParamCurveCartesian2DResult(curve2, title, 0, 1, 2, p); }},
		{"curve3", 3, [&](const auto& p, const auto& title) { return Serializer::SaveParamCurveCartesian3DResult(curve3, title, 0, 1, 2, p); }},
		{"curve2-points", 2, [&](const auto& p, const auto& title) {
			return Serializer::SaveParamCurve<2>(curve2, SerializeFormatType::PARAMETRIC_CURVE_CARTESIAN_2D, title, parameters, p); }},
		{"curve3-points", 3, [&](const auto& p, const auto& title) {
			return Serializer::SaveParamCurve<3>(curve3, SerializeFormatType::PARAMETRIC_CURVE_CARTESIAN_3D, title, parameters, p); }},
		{"raw2", 2, [&](const auto& p, const auto& title) {
			return Serializer::SaveAsParamCurve<2>(points2, SerializeFormatType::PARAMETRIC_CURVE_CARTESIAN_2D, title, 0, 1, 2, p); }},
		{"raw3", 3, [&](const auto& p, const auto& title) {
			return Serializer::SaveAsParamCurve<3>(points3, SerializeFormatType::PARAMETRIC_CURVE_CARTESIAN_3D, title, 0, 1, 2, p); }},
		{"raw2-points", 2, [&](const auto& p, const auto& title) {
			return Serializer::SaveAsParamCurve<2>(points2, SerializeFormatType::PARAMETRIC_CURVE_CARTESIAN_2D, title, parameters, p); }},
		{"raw3-points", 3, [&](const auto& p, const auto& title) {
			return Serializer::SaveAsParamCurve<3>(points3, SerializeFormatType::PARAMETRIC_CURVE_CARTESIAN_3D, title, parameters, p); }},
		{"xy", 2, [&](const auto& p, const auto& title) { return Serializer::SaveAsParamCurve2D(parameters, parameters, title, p); }},
		{"ode2", 2, [&](const auto& p, const auto& title) { return Serializer::SaveODESolAsParametricCurve2D(solution, p, 0, 1, title); }},
		{"ode3", 3, [&](const auto& p, const auto& title) { return Serializer::SaveODESolAsParametricCurve3D(solution, p, 0, 1, 2, title); }}};
	for (const auto& writer : writers) {
		CAPTURE(writer.name);
		for (const auto& title : {std::string{}, std::string("nonempty control")}) {
			MetadataFile file;
			const auto result = writer.write(file.path.string(), title);
			REQUIRE(result.success);
			CHECK(result.filename == file.path.string());
			CheckCurveReaderContract(file.Read(), writer.dimension, title);
		}
	}
}

TEST_CASE("Invalid particle color tokens reject the entire header without output", "[serializer][serializer-metadata]")
{
	for (int dimension : {2, 3}) {
		for (const auto& color : invalidColors) {
			CAPTURE(dimension, color);
			std::ostringstream out;
			out << "SENTINEL" << std::scientific;
			out.precision(3);
			const auto flags = out.flags();
			const auto result = ParticleHeader(out, dimension, {"red", color});
			CHECK_FALSE(result.success);
			CHECK(result.error == IOError::INVALID_PARAMETERS);
			CHECK_FALSE(result.message.empty());
			CHECK(out.str() == "SENTINEL");
			CHECK(out.flags() == flags);
			CHECK(out.precision() == 3);
		}
	}
}

TEST_CASE("Invalid particle catalogs are rejected before creating or truncating files", "[serializer][serializer-metadata]")
{
	for (int dimension : {2, 3}) {
		for (const auto& color : invalidColors) {
			for (bool seed : {false, true}) {
				CAPTURE(dimension, color, seed);
				MetadataFile file;
				if (seed) file.Seed();
				const auto result = ParticleFile(file.path.string(), dimension, {"red", color});
				CHECK_FALSE(result.success);
				CHECK(result.error == IOError::INVALID_PARAMETERS);
				CHECK_FALSE(result.message.empty());
				CHECK(result.filename == file.path.string());
				if (seed) CHECK(file.Read() == "SENTINEL");
				else CHECK_FALSE(std::filesystem::exists(file.path));
			}
		}
	}
}

TEST_CASE("Valid particle catalog tokens retain golden bytes and reader fields", "[serializer][serializer-metadata]")
{
	for (int dimension : {2, 3}) {
		std::ostringstream header;
		REQUIRE(ParticleHeader(header, dimension, {"red", "blue"}).success);
		const std::string expected = std::string(dimension == 2 ? SerializeFormatType::PARTICLE_SIMULATION_DATA_2D :
			SerializeFormatType::PARTICLE_SIMULATION_DATA_3D) + "\nVERSION: 1\nWidth: 10\nHeight: 20\n" +
			(dimension == 3 ? "Depth: 30\n" : "") + "NumBalls: 2\nBall_1 red 1\nBall_2 blue 1\nNumSteps: 1\n";
		CHECK(header.str() == expected);
		for (const auto& color : {std::string("RebeccaPurple"), std::string("#ff00cc"), std::string("custom_token")}) {
			MetadataFile file;
			REQUIRE(ParticleFile(file.path.string(), dimension, {"red", color}).success);
			const auto lines = Lines(file.Read());
			const std::size_t catalog = dimension == 2 ? 5 : 6;
			REQUIRE(lines.size() == catalog + 6);
			for (std::size_t ball = 0; ball < 2; ++ball) {
				std::istringstream record(lines[catalog + ball]);
				std::string name, decoded, extra;
				Real radius;
				REQUIRE(static_cast<bool>(record >> name >> decoded >> radius));
				CHECK(name == "Ball_" + std::to_string(ball + 1));
				CHECK(decoded == (ball == 0 ? "red" : color));
				CHECK(radius == Real(ball + 1));
				CHECK_FALSE(static_cast<bool>(record >> extra));
			}
			CHECK(lines[catalog + 2] == "NumSteps: 1");
			CHECK(lines[catalog + 3] == "Step 0 0");
		}
	}
}
