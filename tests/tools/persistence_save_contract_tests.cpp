#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/Function.h>
#include <mml/tools/Persistence.h>
#endif
#include <catch2/catch_all.hpp>
#include <chrono>
#include <complex>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <vector>

using namespace MML;
using namespace MML::Persistence;

namespace {
	struct SavePath {
		std::filesystem::path path = std::filesystem::temp_directory_path() /
			("mml_save_contract_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".mmlj");
		~SavePath() {
			std::error_code error;
			std::filesystem::remove(path, error);
			if (error) std::cerr << "Cannot remove save-contract file: " << error.message() << '\n';
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
	template<class SaveFunction>
	void RejectWithoutTruncation(SaveFunction save, IOError expected, bool nul = false) {
		SavePath file;
		file.Seed();
		std::string filename = file.path.string();
		if (nul) filename += std::string("\0suffix", 7);
		const auto result = save(filename);
		CHECK_FALSE(result.success);
		CHECK(result.error == expected);
		CHECK_FALSE(result.message.empty());
		CHECK(result.filename == filename);
		CHECK(file.Read() == "SENTINEL");
	}
	template<class Check>
	void DenseObjects(Check check) {
		check(Vector<Real>{1, 2});
		check(VectorN<Real, 2>{1, 2});
		Matrix<Real> matrix(1, 2);
		matrix(0, 0) = 1;
		matrix(0, 1) = 2;
		check(matrix);
		MatrixNM<Real, 1, 2> fixed;
		fixed(0, 0) = 1;
		fixed(0, 1) = 2;
		check(fixed);
	}
	void Nonfinite(Vector<Real>& value, Real invalid) { value[1] = invalid; }
	void Nonfinite(VectorN<Real, 2>& value, Real invalid) { value[1] = invalid; }
	void Nonfinite(Matrix<Real>& value, Real invalid) { value(0, 1) = invalid; }
	void Nonfinite(MatrixNM<Real, 1, 2>& value, Real invalid) { value(0, 1) = invalid; }
}

TEST_CASE("JSON file saves prevalidate dense scalar data and precision", "[Persistence][SaveContracts]") {
	DenseObjects([](auto value) {
		SaveOptions invalid;
		invalid.precision = -1;
		RejectWithoutTruncation([&](const auto& path) { return SaveJson(value, path, invalid); }, IOError::INVALID_PARAMETERS);
		for (Real number : {std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
			Nonfinite(value, number);
			RejectWithoutTruncation([&](const auto& path) { return SaveJson(value, path); }, IOError::UNSUPPORTED_SCALAR);
		}
	});
	Vector<Real> vector{1, std::numeric_limits<Real>::infinity()};
	RejectWithoutTruncation([&](const auto& path) { return Save(vector, path); }, IOError::UNSUPPORTED_SCALAR);
	Matrix<Real> matrix(1, 1);
	matrix(0, 0) = std::numeric_limits<Real>::infinity();
	RejectWithoutTruncation([&](const auto& path) { return Save(matrix, path); }, IOError::UNSUPPORTED_SCALAR);
}

TEST_CASE("Matrix file dispatch prevalidates textual options and unsupported formats", "[Persistence][SaveContracts]") {
	Matrix<Real> matrix(1, 1);
	matrix(0, 0) = 1;
	SaveOptions options;
	options.precision = -1;
	for (Format format : {Format::Text, Format::CSV}) {
		RejectWithoutTruncation([&](const auto& path) { return Save(matrix, path, options, format); }, IOError::INVALID_PARAMETERS);
	}
	RejectWithoutTruncation([&](const auto& path) { return Save(matrix, path, SaveOptions{}, static_cast<Format>(99)); }, IOError::INVALID_FORMAT);
}

TEST_CASE("Sampled-function file saves reject invalid data before opening", "[Persistence][SaveContracts]") {
	const SampledRealFunctionData valid{Vector<Real>{0, 1}, Vector<Real>{2, 3}, "test", "fixture"};
	auto data = valid;
	data.values.resize(1);
	RejectWithoutTruncation([&](const auto& path) { return SaveSampledFunction(data, path); }, IOError::INVALID_PARAMETERS);
	data = valid;
	data.values[1] = std::numeric_limits<Real>::infinity();
	RejectWithoutTruncation([&](const auto& path) { return SaveSampledFunction(data, path); }, IOError::UNSUPPORTED_SCALAR);
	SaveOptions options;
	options.precision = -1;
	RejectWithoutTruncation([&](const auto& path) { return SaveSampledFunction(valid, path, options); }, IOError::INVALID_PARAMETERS);
}

TEST_CASE("Function file saves prevalidate grids and options without truncating", "[Persistence][SaveContracts]") {
	RealFunction function([](Real x) { return x * x; });
	RejectWithoutTruncation([&](const auto& path) {
		return SaveSampledFunction(function, "test", SampleGrid{0, 1, 1}, path);
	}, IOError::INVALID_PARAMETERS);
	SaveOptions options;
	options.precision = -1;
	RejectWithoutTruncation([&](const auto& path) {
		return SaveSampledFunction(function, "test", SampleGrid{0, 1, 3}, path, options);
	}, IOError::INVALID_PARAMETERS);
	LinearInterpRealFunc interpolated(Vector<Real>{0, 1}, Vector<Real>{1, 2});
	RejectWithoutTruncation([&](const auto& path) {
		return SaveInterpolatedFunction(interpolated, path, "test", options);
	}, IOError::INVALID_PARAMETERS);
}

TEST_CASE("Persistence save filenames reject embedded NUL before opening prefix paths", "[Persistence][SaveContracts]") {
	DenseObjects([](const auto& value) {
		RejectWithoutTruncation([&](const auto& path) { return SaveJson(value, path); }, IOError::INVALID_PARAMETERS, true);
	});
	Vector<Real> vector{1};
	Matrix<Real> matrix(1, 1);
	matrix(0, 0) = 1;
	for (Format format : {Format::Text, Format::CSV, Format::Binary, Format::JSON}) {
		RejectWithoutTruncation([&](const auto& path) { return Save(matrix, path, format); }, IOError::INVALID_PARAMETERS, true);
		RejectWithoutTruncation([&](const auto& path) { return Save(matrix, path, SaveOptions{}, format); }, IOError::INVALID_PARAMETERS, true);
	}
	RejectWithoutTruncation([&](const auto& path) { return SaveBinary(vector, path); }, IOError::INVALID_PARAMETERS, true);
	RejectWithoutTruncation([&](const auto& path) { return SaveMatrixToFile(matrix, path); }, IOError::INVALID_PARAMETERS, true);
	RejectWithoutTruncation([&](const auto& path) { return SaveMatrixToCSV(matrix, path); }, IOError::INVALID_PARAMETERS, true);
	const SampledRealFunctionData data{Vector<Real>{0, 1}, Vector<Real>{1, 2}, "", ""};
	RejectWithoutTruncation([&](const auto& path) { return SaveSampledFunction(data, path); }, IOError::INVALID_PARAMETERS, true);
	const std::vector<std::complex<double>> complex{{1, 2}};
	RejectWithoutTruncation([&](const auto& path) { return SaveComplexVector(complex, path); }, IOError::INVALID_PARAMETERS, true);
	RealFunction function([](Real x) { return x; });
	RejectWithoutTruncation([&](const auto& path) {
		return SaveSampledFunction(function, "test", SampleGrid{}, path);
	}, IOError::INVALID_PARAMETERS, true);
	LinearInterpRealFunc interpolated(Vector<Real>{0, 1}, Vector<Real>{1, 2});
	RejectWithoutTruncation([&](const auto& path) {
		return SaveInterpolatedFunction(interpolated, path);
	}, IOError::INVALID_PARAMETERS, true);
}

TEST_CASE("Invalid persistence saves do not create output files", "[Persistence][SaveContracts]") {
	SavePath file;
	REQUIRE_FALSE(std::filesystem::exists(file.path));
	Vector<Real> invalid{std::numeric_limits<Real>::infinity()};
	CHECK(SaveJson(invalid, file.path.string()).error == IOError::UNSUPPORTED_SCALAR);
	CHECK_FALSE(std::filesystem::exists(file.path));
	Matrix<Real> matrix(1, 1);
	SaveOptions options;
	options.precision = -1;
	CHECK(Save(matrix, file.path.string(), options, Format::Text).error == IOError::INVALID_PARAMETERS);
	CHECK_FALSE(std::filesystem::exists(file.path));
	const SampledRealFunctionData data{Vector<Real>{0, 0}, Vector<Real>{1, 2}, "", ""};
	CHECK(SaveSampledFunction(data, file.path.string()).error == IOError::INVALID_PARAMETERS);
	CHECK_FALSE(std::filesystem::exists(file.path));
}

TEST_CASE("Function sampling completes exactly once before the output is opened", "[Persistence][SaveContracts]") {
	SavePath file;
	file.Seed();
	int calls = 0;
	RealFunctionFromStdFunc function([&](Real x) {
		CHECK(file.Read() == "SENTINEL");
		++calls;
		return x * x;
	});
	const auto result = SaveSampledFunction(function, "quadratic", SampleGrid{0, 1, 5}, file.path.string());
	REQUIRE(result.success);
	CHECK(result.filename == file.path.string());
	CHECK(calls == 5);
	SampledRealFunctionData restored;
	REQUIRE(LoadSampledFunctionData(file.path.string(), restored).success);
	REQUIRE(restored.nodes.size() == 5);
	CHECK(restored.label == "quadratic");
	CHECK(restored.source == "sampled from IRealFunction");
	for (int i = 0; i < 5; ++i) {
		CHECK(restored.nodes[i] == Real(i) / 4);
		CHECK(restored.values[i] == restored.nodes[i] * restored.nodes[i]);
	}
}

TEST_CASE("Invalid function requests do not evaluate callbacks or alter output", "[Persistence][SaveContracts]") {
	int calls = 0;
	RealFunctionFromStdFunc function([&](Real x) { ++calls; return x; });
	SaveOptions options;
	options.precision = -1;
	RejectWithoutTruncation([&](const auto& path) {
		return SaveSampledFunction(function, "", SampleGrid{}, path, options);
	}, IOError::INVALID_PARAMETERS);
	CHECK(calls == 0);
	for (SampleGrid grid : {SampleGrid{0, 1, 1}, SampleGrid{1, 0, 2},
		SampleGrid{std::numeric_limits<Real>::quiet_NaN(), 1, 2},
		SampleGrid{0, std::numeric_limits<Real>::infinity(), 2},
		SampleGrid{-std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max(), 2},
		SampleGrid{1, std::nextafter(Real(1), Real(2)), 3}}) {
		RejectWithoutTruncation([&](const auto& path) { return SaveSampledFunction(function, "", grid, path); }, IOError::INVALID_PARAMETERS);
		CHECK(calls == 0);
	}
	RejectWithoutTruncation([&](const auto& path) {
		return SaveSampledFunction(function, "", SampleGrid{}, path);
	}, IOError::INVALID_PARAMETERS, true);
	CHECK(calls == 0);
}

TEST_CASE("Nonfinite results and callback exceptions preserve existing function files", "[Persistence][SaveContracts]") {
	RealFunction invalid([](Real) { return std::numeric_limits<Real>::quiet_NaN(); });
	RejectWithoutTruncation([&](const auto& path) {
		return SaveSampledFunction(invalid, "", SampleGrid{}, path);
	}, IOError::UNSUPPORTED_SCALAR);
	SavePath file;
	file.Seed();
	RealFunction throwing([](Real) -> Real { throw std::runtime_error("callback failed"); });
	CHECK_THROWS_AS(SaveSampledFunction(throwing, "", SampleGrid{}, file.path.string()), std::runtime_error);
	CHECK(file.Read() == "SENTINEL");
	for (Real invalidNode : {std::numeric_limits<Real>::infinity(), std::numeric_limits<Real>::quiet_NaN()}) {
		const SampledRealFunctionData data{Vector<Real>{0, invalidNode}, Vector<Real>{1, 2}, "", ""};
		RejectWithoutTruncation([&](const auto& path) { return SaveSampledFunction(data, path); }, IOError::UNSUPPORTED_SCALAR);
	}
}

TEST_CASE("Valid JSON file helpers preserve stream encodings and finish every function route", "[Persistence][SaveContracts]") {
	DenseObjects([](const auto& value) {
		for (bool pretty : {false, true}) {
			SaveOptions options;
			options.pretty_json = pretty;
			options.precision = 3;
			std::ostringstream expected;
			REQUIRE(SaveJson(expected, value, options).success);
			SavePath file;
			const auto result = SaveJson(value, file.path.string(), options);
			REQUIRE(result.success);
			CHECK(result.filename == file.path.string());
			std::ifstream in(file.path);
			std::ostringstream actual;
			actual << in.rdbuf();
			CHECK(actual.str() == expected.str());
		}
	});
	const SampledRealFunctionData data{Vector<Real>{1, 0}, Vector<Real>{2, 3}, "descending", "fixture"};
	SavePath direct;
	REQUIRE(SaveSampledFunction(data, direct.path.string()).success);
	SampledRealFunctionData restored;
	REQUIRE(LoadSampledFunctionData(direct.path.string(), restored).success);
	CHECK(restored.nodes == data.nodes);
	CHECK(restored.values == data.values);
	CHECK(restored.label == data.label);
	LinearInterpRealFunc function(Vector<Real>{0, 1}, Vector<Real>{1, 2});
	SavePath interpolated;
	const auto result = SaveInterpolatedFunction(function, interpolated.path.string(), "linear");
	REQUIRE(result.success);
	CHECK(result.filename == interpolated.path.string());
	REQUIRE(LoadSampledFunctionData(interpolated.path.string(), restored).success);
	CHECK(restored.label == "linear");
	CHECK(restored.source == function.InterpolationMethodName());
}

TEST_CASE("All persistence file routes reject empty names with explicit context", "[Persistence][SaveContracts]") {
	DenseObjects([](const auto& value) { CHECK(SaveJson(value, "").error == IOError::INVALID_PARAMETERS); });
	Vector<Real> vector{1};
	Matrix<Real> matrix(1, 1);
	for (Format format : {Format::Text, Format::CSV, Format::JSON, Format::Binary}) {
		CHECK(Save(matrix, "", format).error == IOError::INVALID_PARAMETERS);
		CHECK(Save(matrix, "", SaveOptions{}, format).error == IOError::INVALID_PARAMETERS);
	}
	CHECK(Save(vector, "", Format::Text).error == IOError::INVALID_PARAMETERS);
	CHECK(SaveBinary(vector, "").error == IOError::INVALID_PARAMETERS);
	const SampledRealFunctionData data{Vector<Real>{0, 1}, Vector<Real>{1, 2}, "", ""};
	CHECK(SaveSampledFunction(data, "").error == IOError::INVALID_PARAMETERS);
}

TEST_CASE("Binary save preflight retains existing envelope limits without allocating payloads", "[Persistence][SaveContracts]") {
	const auto limit = static_cast<uint64_t>(LoadOptions{}.max_allocation_bytes);
	BinaryEnvelopeHeader header;
	const uint64_t vectorCount = (limit - 8) / sizeof(double);
	CHECK(Persistence::Detail::BuildVectorBinaryHeader<double>(vectorCount, header).success);
	CHECK(header.payload_byte_count == limit);
	CHECK(Persistence::Detail::BuildVectorBinaryHeader<double>(vectorCount + 1, header).error == IOError::ALLOCATION_LIMIT_EXCEEDED);
	CHECK(Persistence::Detail::BuildVectorBinaryHeader<double>(std::numeric_limits<uint64_t>::max(), header).error == IOError::ALLOCATION_LIMIT_EXCEEDED);
	const uint64_t matrixCols = (limit - 16) / sizeof(double);
	CHECK(Persistence::Detail::BuildMatrixBinaryHeader<double>(1, matrixCols, header).success);
	CHECK(header.payload_byte_count == limit);
	CHECK(Persistence::Detail::BuildMatrixBinaryHeader<double>(1, matrixCols + 1, header).error == IOError::ALLOCATION_LIMIT_EXCEEDED);
	CHECK(Persistence::Detail::BuildMatrixBinaryHeader<double>(std::numeric_limits<uint64_t>::max(), 2, header).error == IOError::ALLOCATION_LIMIT_EXCEEDED);
	SaveOptions options;
	options.precision = -1; // Binary intentionally ignores text precision.
	SavePath file;
	Vector<Real> vector{1, 2};
	REQUIRE(Save(vector, file.path.string(), options, Format::Binary).success);
	Vector<Real> restored;
	REQUIRE(LoadBinary(file.path.string(), restored).success);
	CHECK(restored == vector);
	Matrix<Real> matrix(1, 1);
	matrix(0, 0) = 3;
	REQUIRE(Save(matrix, file.path.string(), options, Format::Binary).success);
	Matrix<Real> matrixRestored;
	REQUIRE(LoadBinary(file.path.string(), matrixRestored).success);
	CHECK(matrixRestored(0, 0) == 3);
}

TEST_CASE("Persistence output finalization detects close failure and preserves original error context", "[Persistence][SaveContracts]") {
	SavePath file;
	std::ofstream out(file.path);
	out << "written";
	REQUIRE(out.good());
	REQUIRE(out.rdbuf()->close() != nullptr);
	const auto result = FinishOutputFile(out, file.path.string(), SerializeSuccess());
	CHECK_FALSE(result.success);
	CHECK(result.error == IOError::WRITE_FAILED);
	CHECK(result.filename == file.path.string());
	CHECK_FALSE(result.message.empty());
	std::ofstream failed(file.path);
	failed.setstate(std::ios::badbit);
	const auto original = SerializeFailure(IOError::INVALID_PARAMETERS, "original validation failure").At(7);
	const auto preserved = FinishOutputFile(failed, file.path.string(), original);
	CHECK(preserved.error == original.error);
	CHECK(preserved.message == original.message);
	CHECK(preserved.offset == 7);
	CHECK(preserved.filename == file.path.string());
}

TEST_CASE("Persistence output finalization reports a delayed flush failure after successful serialization", "[Persistence][SaveContracts]") {
	struct FailingSyncBuffer : std::stringbuf {
		int sync() override { return -1; }
	} buffer;
	SavePath path;
	std::ofstream file(path.path);
	file << "written";
	REQUIRE(file.good());
	std::ostream& stream = file;
	auto* original = stream.rdbuf(&buffer);
	REQUIRE(file.good());
	const auto result = FinishOutputFile(file, path.path.string(), SerializeSuccess());
	CHECK(result.error == IOError::WRITE_FAILED);
	CHECK_FALSE(result.success);
	CHECK(result.filename == path.path.string());
	CHECK(result.offset == 7);
	CHECK(result.offset_kind == IOOffsetKind::Byte);
	CHECK_FALSE(file.is_open());
	stream.rdbuf(original);
}

TEST_CASE("Valid generic persistence saves retain bytes for explicit formats and default helpers", "[Persistence][SaveContracts]") {
	Matrix<Real> matrix(1, 2);
	matrix(0, 0) = 1;
	matrix(0, 1) = 2;
	for (Format format : {Format::Text, Format::CSV, Format::JSON, Format::Binary}) {
		std::ostringstream expected;
		REQUIRE(Save(expected, matrix, format).success);
		SavePath file;
		REQUIRE(Save(matrix, file.path.string(), format).success);
		std::ifstream disk(file.path, format == Format::Binary ? std::ios::binary : std::ios::in);
		std::ostringstream actual;
		actual << disk.rdbuf();
		CHECK(actual.str() == expected.str());
		disk.close();
		REQUIRE(Save(matrix, file.path.string(), SaveOptions{}, format).success);
	}
	Vector<Real> vector{1, 2};
	for (Format format : {Format::Text, Format::JSON, Format::Binary}) {
		std::ostringstream expected;
		REQUIRE(Save(expected, vector, format).success);
		SavePath file;
		REQUIRE(Save(vector, file.path.string(), format).success);
		std::ifstream disk(file.path, format == Format::Binary ? std::ios::binary : std::ios::in);
		std::ostringstream actual;
		actual << disk.rdbuf();
		CHECK(actual.str() == expected.str());
	}
}

TEST_CASE("Persistence stream write failures remain explicit and invalid JSON does not touch streams", "[Persistence][SaveContracts]") {
	struct RefusingBuffer : std::streambuf {
		std::streamsize xsputn(const char*, std::streamsize) override { return 0; }
		int_type overflow(int_type) override { return traits_type::eof(); }
	} buffer;
	std::ostream out(&buffer);
	const SampledRealFunctionData data{Vector<Real>{0, 1}, Vector<Real>{1, 2}, "", ""};
	CHECK(SaveSampledFunction(out, data).error == IOError::WRITE_FAILED);
	DenseObjects([](auto value) {
		Nonfinite(value, std::numeric_limits<Real>::infinity());
		std::ostringstream stream;
		stream << "SENTINEL";
		CHECK(SaveJson(stream, value).error == IOError::UNSUPPORTED_SCALAR);
		CHECK(stream.str() == "SENTINEL");
	});
	CHECK(Persistence::Detail::JsonNumberIsWritable(std::numeric_limits<double>::max()));
	CHECK(Persistence::Detail::JsonNumberIsWritable(std::numeric_limits<double>::lowest()));
	CHECK(Persistence::Detail::JsonNumberIsWritable(std::numeric_limits<long double>::denorm_min()));
	if constexpr (std::numeric_limits<long double>::max_exponent > std::numeric_limits<double>::max_exponent)
		CHECK_FALSE(Persistence::Detail::JsonNumberIsWritable(std::numeric_limits<long double>::max()));
}
