#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/tools/serializer/SerializerSimulation.h>
#endif
#include <catch2/catch_all.hpp>
#include <chrono>
#include <cstdint>
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
	struct ParticleFile {
		std::filesystem::path path = std::filesystem::temp_directory_path() /
			("mml_particle_contract_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".mml");
		~ParticleFile() {
			std::error_code error;
			std::filesystem::remove(path, error);
			if (error) std::cerr << "Cannot remove particle-contract file: " << error.message() << '\n';
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
	template<class Point>
	Point Position(int particle, int step) {
		if constexpr (std::is_same_v<Point, Pnt2Cart>) return Point(Real(step), Real(particle));
		else return Point(Real(step), Real(particle), Real(particle + step));
	}
	template<class Point>
	std::vector<std::vector<Point>> History(const std::vector<int>& lengths) {
		std::vector<std::vector<Point>> histories(lengths.size());
		for (std::size_t particle = 0; particle < lengths.size(); ++particle)
			for (int step = 0; step < lengths[particle]; ++step)
				histories[particle].push_back(Position<Point>(static_cast<int>(particle), step));
		return histories;
	}
	template<class Point>
	IOResult SaveHistory(const std::string& filename, const std::vector<std::vector<Point>>& histories,
						 int stride = 1, Real timeStep = Real(0.25)) {
		const int particles = static_cast<int>(histories.size());
		const std::vector<std::string> colors(histories.size(), "red");
		const std::vector<Real> radii(histories.size(), Real(0.5));
		if constexpr (std::is_same_v<Point, Pnt2Cart>)
			return Serializer::SaveParticleSimulation2D(filename, particles, 10, 20, histories, colors, radii, timeStep, stride);
		else
			return Serializer::SaveParticleSimulation3D(filename, particles, 10, 20, 30, histories, colors, radii, timeStep, stride);
	}
	template<class Point>
	void RejectHistory(const std::vector<std::vector<Point>>& history, bool seed, int stride = 1) {
		ParticleFile file;
		if (seed) file.Seed();
		const auto result = SaveHistory(file.path.string(), history, stride);
		CHECK_FALSE(result.success);
		CHECK(result.error == IOError::INVALID_PARAMETERS);
		CHECK_FALSE(result.message.empty());
		CHECK(result.filename == file.path.string());
		if (seed) CHECK(file.Read() == "SENTINEL");
		else CHECK_FALSE(std::filesystem::exists(file.path));
	}
	template<class Point>
	void CheckExport(int steps, int stride, Real timeStep) {
		ParticleFile file;
		const auto histories = History<Point>({steps, steps});
		const auto result = SaveHistory(file.path.string(), histories, stride, timeStep);
		REQUIRE(result.success);
		CHECK(result.filename == file.path.string());
		std::istringstream in(file.Read());
		std::string line;
		REQUIRE(static_cast<bool>(std::getline(in, line)));
		CHECK(line == (std::is_same_v<Point, Pnt2Cart> ? "MML_PARTICLE_SIMULATION_DATA_2D" : "MML_PARTICLE_SIMULATION_DATA_3D"));
		for (const char* header : {"VERSION: 1", "Width: 10", "Height: 20"}) {
			REQUIRE(static_cast<bool>(std::getline(in, line)));
			CHECK(line == header);
		}
		if constexpr (std::is_same_v<Point, Pnt3Cart>) {
			REQUIRE(static_cast<bool>(std::getline(in, line)));
			CHECK(line == "Depth: 30");
		}
		for (const char* header : {"NumBalls: 2", "Ball_1 red 0.5", "Ball_2 red 0.5"}) {
			REQUIRE(static_cast<bool>(std::getline(in, line)));
			CHECK(line == header);
		}
		const auto expectedFrames = (static_cast<std::uint64_t>(steps) + static_cast<std::uint64_t>(stride) - 1) /
			static_cast<std::uint64_t>(stride);
		REQUIRE(static_cast<bool>(std::getline(in, line)));
		CHECK(line == "NumSteps: " + std::to_string(expectedFrames));
		int frame = 0;
		for (std::size_t step = 0; step < static_cast<std::size_t>(steps); step += static_cast<std::size_t>(stride), ++frame) {
			REQUIRE(static_cast<bool>(std::getline(in, line)));
			std::istringstream stepLine(line);
			std::string label;
			int emittedFrame = -1;
			Real timestamp = -1;
			REQUIRE(static_cast<bool>(stepLine >> label >> emittedFrame >> timestamp));
			CHECK(label == "Step");
			CHECK(emittedFrame == frame);
			CHECK(timestamp == static_cast<Real>(step) * timeStep);
			stepLine >> std::ws;
			CHECK(stepLine.eof());
			for (int particle = 0; particle < 2; ++particle) {
				REQUIRE(static_cast<bool>(std::getline(in, line)));
				std::istringstream positionLine(line);
				int emittedParticle = -1;
				Real x = -1, y = -1, z = -1;
				REQUIRE(static_cast<bool>(positionLine >> emittedParticle >> x >> y));
				CHECK(emittedParticle == particle);
				CHECK(x == static_cast<Real>(step));
				CHECK(y == particle);
				if constexpr (std::is_same_v<Point, Pnt3Cart>) {
					REQUIRE(static_cast<bool>(positionLine >> z));
					CHECK(z == static_cast<Real>(particle) + static_cast<Real>(step));
				}
				positionLine >> std::ws;
				CHECK(positionLine.eof());
			}
		}
		CHECK(static_cast<std::uint64_t>(frame) == expectedFrames);
		in >> std::ws;
		CHECK(in.eof());
	}
}

TEMPLATE_TEST_CASE("Particle export rejects longer later histories without truncating files",
	"[serializer][ParticleHistory]", Pnt2Cart, Pnt3Cart) {
	for (const auto& lengths : {std::vector<int>{1, 2}, std::vector<int>{0, 1}, std::vector<int>{1, 1, 2}}) {
		CAPTURE(lengths);
		RejectHistory(History<TestType>(lengths), true);
		RejectHistory(History<TestType>(lengths), false);
		RejectHistory(History<TestType>(lengths), true, std::numeric_limits<int>::max());
		RejectHistory(History<TestType>(lengths), false, std::numeric_limits<int>::max());
	}
}

TEMPLATE_TEST_CASE("Particle export rejects shorter and empty later histories before indexing",
	"[serializer][ParticleHistory]", Pnt2Cart, Pnt3Cart) {
	for (const auto& lengths : {std::vector<int>{2, 1}, std::vector<int>{1, 0}, std::vector<int>{2, 2, 1}}) {
		CAPTURE(lengths);
		RejectHistory(History<TestType>(lengths), true);
		RejectHistory(History<TestType>(lengths), false);
		RejectHistory(History<TestType>(lengths), true, std::numeric_limits<int>::max());
		RejectHistory(History<TestType>(lengths), false, std::numeric_limits<int>::max());
	}
}

TEMPLATE_TEST_CASE("Particle export retains rectangular layouts frame ceilings and original timestamps",
	"[serializer][ParticleHistory]", Pnt2Cart, Pnt3Cart) {
	for (int steps : {0, 1, 2, 5, 6, 7}) {
		for (int stride : {1, 2, 3, 7, std::numeric_limits<int>::max() - 1, std::numeric_limits<int>::max()}) {
			CAPTURE(steps, stride);
			CheckExport<TestType>(steps, stride, Real(0.25));
		}
	}
	CheckExport<TestType>(7, 2, Real(0.1));
}

TEMPLATE_TEST_CASE("Particle export rejects nonpositive strides without changing files",
	"[serializer][ParticleHistory]", Pnt2Cart, Pnt3Cart) {
	const auto history = History<TestType>({2, 2});
	for (int stride : {0, -1, std::numeric_limits<int>::min()}) {
		CAPTURE(stride);
		RejectHistory(history, true, stride);
		RejectHistory(history, false, stride);
	}
}

TEST_CASE("Particle frame counts cover int limits without allocating histories", "[serializer][ParticleHistory]") {
	const auto maximum = static_cast<std::size_t>(std::numeric_limits<int>::max());
	for (std::size_t steps : {std::size_t(0), std::size_t(1), std::size_t(2), maximum - 1, maximum}) {
		for (int stride : {1, 2, 3, std::numeric_limits<int>::max() / 2, std::numeric_limits<int>::max()}) {
			CAPTURE(steps, stride);
			int frames = -1;
			REQUIRE(Serializer::detail::ParticleSimulationFrameCount(steps, stride, frames).success);
			const auto oracle = (static_cast<std::uint64_t>(steps) + static_cast<std::uint64_t>(stride) - 1) /
				static_cast<std::uint64_t>(stride);
			CHECK(static_cast<std::uint64_t>(frames) == oracle);
			if (frames > 0) {
				const auto lastIndex = static_cast<std::size_t>(frames - 1) * static_cast<std::size_t>(stride);
				CHECK(lastIndex < steps);
				CHECK(steps - lastIndex <= static_cast<std::size_t>(stride));
			}
		}
	}
	for (std::size_t steps : {maximum + 1, std::numeric_limits<std::size_t>::max()}) {
		int frames = 77;
		const auto result = Serializer::detail::ParticleSimulationFrameCount(steps, 1, frames);
		CHECK(result.error == IOError::INVALID_PARAMETERS);
		CHECK_FALSE(result.message.empty());
		CHECK(frames == 77);
	}
	for (int stride : {0, -1, std::numeric_limits<int>::min()}) {
		int frames = 77;
		CHECK(Serializer::detail::ParticleSimulationFrameCount(2, stride, frames).error == IOError::INVALID_PARAMETERS);
		CHECK(frames == 77);
	}
}

TEST_CASE("Particle history validation preserves zero-frame rectangles and validates empty outer histories", "[serializer][ParticleHistory]") {
	int frames = 77;
	CHECK(Serializer::detail::ValidateParticleHistories(std::vector<std::vector<Pnt2Cart>>{}, 1, frames).error == IOError::INVALID_PARAMETERS);
	CHECK(frames == 77);
	CHECK(Serializer::detail::ValidateParticleHistories(History<Pnt2Cart>({0, 0}), 1, frames).success);
	CHECK(frames == 0);
}
