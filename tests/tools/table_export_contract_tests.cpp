#include <catch2/catch_all.hpp>

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/tools/ConsolePrinter.h>
#endif

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <streambuf>
#include <string>

using namespace MML;

namespace {
	struct ExportPath {
		std::filesystem::path directory = std::filesystem::temp_directory_path() /
			("mml_table_export_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
		std::filesystem::path file;
		ExportPath() {
			std::filesystem::create_directory(directory);
			file = directory / "table.txt";
		}
		~ExportPath() {
			std::error_code error;
			std::filesystem::remove_all(directory, error);
			if (error)
				std::cerr << "Cannot clean table export fixture: " << error.message() << '\n';
		}
	};

	enum class FailureMode { Write, Flush, Close };

	struct ControlledBuffer : std::stringbuf {
		FailureMode mode;
		int syncCalls = 0;
		explicit ControlledBuffer(FailureMode failure) : mode(failure) {}
		std::streamsize xsputn(const char* data, std::streamsize count) override {
			return mode == FailureMode::Write ? 0 : std::stringbuf::xsputn(data, count);
		}
		int_type overflow(int_type character = traits_type::eof()) override {
			return mode == FailureMode::Write ? traits_type::eof() : std::stringbuf::overflow(character);
		}
		int sync() override {
			++syncCalls;
			return mode == FailureMode::Flush ? -1 : 0;
		}
	};

	struct ControlledTag {
		ControlledBuffer* buffer;
		bool* injected;
		friend std::ostream& operator<<(std::ostream& output, const ControlledTag& tag) {
			auto* file = dynamic_cast<std::ofstream*>(&output);
			if (file != nullptr) {
				REQUIRE(file->is_open());
				REQUIRE(file->good());
				if (tag.buffer->mode == FailureMode::Close)
					REQUIRE(file->rdbuf()->close() != nullptr);
				output.rdbuf(tag.buffer);
				*tag.injected = true;
			}
			return output << "controlled tag";
		}
	};

	struct ThrowingTag {
		friend std::ostream& operator<<(std::ostream&, const ThrowingTag&) {
			throw ArgumentError("custom table value cannot be rendered");
		}
	};

	std::string TextFileBytes(const std::string& text) {
#ifdef _WIN32
		std::string result;
		for (char character : text) {
			if (character == '\n') result += '\r';
			result += character;
		}
		return result;
#else
		return text;
#endif
	}

	std::string ReadBytes(const std::filesystem::path& path) {
		std::ifstream input(path, std::ios::binary);
		REQUIRE(input.is_open());
		return std::string(std::istreambuf_iterator<char>(input), {});
	}

	template<class Write>
	void CheckFileError(Write write, const std::string& filename) {
		bool caught = false;
		try {
			write();
		} catch (const FileIOError& error) {
			caught = true;
			CHECK(error.filename() == filename);
			CHECK(std::string(error.what()).find(filename) != std::string::npos);
			CHECK(std::string(error.what()).find("WRITE_FAILED") != std::string::npos);
		}
		CHECK(caught);
	}
}

TEST_CASE("Table file exports report controlled write flush and close failures",
	"[ConsolePrinter][TableExportContract]") {
	for (auto mode : {FailureMode::Write, FailureMode::Flush, FailureMode::Close}) {
		CAPTURE(mode);
		ExportPath path;
		ControlledBuffer buffer(mode);
		bool injected = false;
		TablePrinter<ControlledTag, int> table("tag", {"value"});
		table.addRow({&buffer, &injected}, {42});
		const auto filename = path.file.string();
		CheckFileError([&] { table.exportToFile(filename, ExportFormat::TSV); }, filename);
		CHECK(injected);
		if (mode == FailureMode::Flush || mode == FailureMode::Close)
			CHECK(buffer.syncCalls == 1);
		if (mode == FailureMode::Close)
			CHECK(buffer.str() == "controlled tag\t42\n");
	}
}

TEST_CASE("Valid owned table files preserve every existing export encoding and close on return",
	"[ConsolePrinter][TableExportContract]") {
	TablePrinter<std::string, std::string> table("Tag", {"Value"});
	table.addRow("\xCE\xB1", {"A,\"B\"\nC<&|_"});
	table.style(TableStyle().border(BorderStyle::Rounded));
	for (auto format : {ExportFormat::Console, ExportFormat::CSV, ExportFormat::TSV,
		ExportFormat::Markdown, ExportFormat::LaTeX, ExportFormat::HTML}) {
		CAPTURE(format);
		ExportPath path;
		std::ostringstream expected;
		table.exportTo(expected, format);
		table.exportToFile(path.file.string(), format);
		CHECK(ReadBytes(path.file) == TextFileBytes(expected.str()));
		const auto renamed = path.directory / "renamed.txt";
		REQUIRE_NOTHROW(std::filesystem::rename(path.file, renamed));
		CHECK(ReadBytes(renamed) == TextFileBytes(expected.str()));
	}
	ExportPath path;
	table.printToFile(path.file.string());
	CHECK(ReadBytes(path.file) == TextFileBytes(table.toString()));
}

TEST_CASE("Table file helpers preserve explicit open-error filename diagnostics",
	"[ConsolePrinter][TableExportContract]") {
	ExportPath path;
	const auto filename = (path.directory / "missing" / "table.txt").string();
	const TablePrinter<int, int> table("tag", {"value"});
	for (bool console : {false, true}) {
		bool caught = false;
		try {
			if (console) table.printToFile(filename);
			else table.exportToFile(filename, ExportFormat::CSV);
		} catch (const FileIOError& error) {
			caught = true;
			CHECK(error.filename() == filename);
			CHECK(std::string(error.what()).find(filename) != std::string::npos);
		}
		CHECK(caught);
	}
	CHECK_FALSE(std::filesystem::exists(path.directory / "missing"));
}

TEST_CASE("Caller-owned table streams retain nonthrowing failure state without being finalized",
	"[ConsolePrinter][TableExportContract]") {
	TablePrinter<int, int> table("tag", {"value"});
	table.addRow(1, {2});
	for (auto format : {ExportFormat::Console, ExportFormat::CSV, ExportFormat::TSV,
		ExportFormat::Markdown, ExportFormat::LaTeX, ExportFormat::HTML}) {
		ControlledBuffer buffer(FailureMode::Write);
		std::ostream output(&buffer);
		output << std::scientific << std::showpos;
		output.precision(8);
		output.width(9);
		output.fill('@');
		const auto flags = output.flags();
		CHECK_NOTHROW(table.exportTo(output, format));
		CHECK(output.fail());
		CHECK(buffer.syncCalls == 0);
		CHECK(output.flags() == flags);
		CHECK(output.precision() == 8);
		CHECK(output.width() == 9);
		CHECK(output.fill() == '@');
	}
	ControlledBuffer buffer(FailureMode::Flush);
	std::ostream output(&buffer);
	table.print(output);
	CHECK(output.good());
	CHECK(buffer.syncCalls == 0);
	CHECK(buffer.str().find("2") != std::string::npos);
}

TEST_CASE("Caller-owned table streams retain their requested stream exception behavior",
	"[ConsolePrinter][TableExportContract]") {
	TablePrinter<int, int> table("tag", {"value"});
	table.addRow(1, {2});
	for (auto format : {ExportFormat::Console, ExportFormat::CSV, ExportFormat::TSV,
		ExportFormat::Markdown, ExportFormat::LaTeX, ExportFormat::HTML}) {
		ControlledBuffer buffer(FailureMode::Write);
		std::ostream output(&buffer);
		output.precision(8);
		output.width(9);
		output.exceptions(std::ios::badbit | std::ios::failbit);
		CHECK_THROWS_AS(table.exportTo(output, format), std::ios_base::failure);
		CHECK(output.precision() == 8);
		CHECK(output.width() == 9);
		CHECK(output.exceptions() == (std::ios::badbit | std::ios::failbit));
		CHECK(buffer.syncCalls == 0);
	}
}

TEST_CASE("Empty tables finish successfully without manufacturing a write failure",
	"[ConsolePrinter][TableExportContract]") {
	const TablePrinter<int, int> table("tag", {});
	for (auto format : {ExportFormat::Console, ExportFormat::CSV, ExportFormat::TSV,
		ExportFormat::Markdown, ExportFormat::LaTeX, ExportFormat::HTML}) {
		ExportPath path;
		std::ostringstream expected;
		table.exportTo(expected, format);
		CHECK_NOTHROW(table.exportToFile(path.file.string(), format));
		CHECK(ReadBytes(path.file) == TextFileBytes(expected.str()));
	}
	TablePrinter<int, int> blank("tag", {});
	blank.style(TableStyle().border(BorderStyle::None).header(false));
	ExportPath path;
	CHECK_NOTHROW(blank.printToFile(path.file.string()));
	CHECK(ReadBytes(path.file).empty());
}

TEST_CASE("File export preserves exceptions raised by custom table value rendering",
	"[ConsolePrinter][TableExportContract]") {
	TablePrinter<ThrowingTag, int> table("tag", {"value"});
	table.addRow({}, {42});
	for (bool console : {false, true}) {
		ExportPath path;
		bool caught = false;
		try {
			if (console) table.printToFile(path.file.string());
			else table.exportToFile(path.file.string(), ExportFormat::CSV);
		} catch (const ArgumentError& error) {
			caught = true;
			CHECK(std::string(error.what()) == "custom table value cannot be rendered");
		}
		CHECK(caught);
		const auto renamed = path.directory / "after_exception.txt";
		CHECK_NOTHROW(std::filesystem::rename(path.file, renamed));
	}
}
