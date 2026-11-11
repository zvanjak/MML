#include <catch2/catch_all.hpp>

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/tools/ConsolePrinter.h>
#endif

#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

using namespace MML;

TEST_CASE("Column widths reject invalid and border-overflowing values before changing state",
	"[ConsolePrinter][ColumnFormatContract]") {
	ColumnFormat column("value");
	column.width(7);
	for (int width : {-2, std::numeric_limits<int>::min(),
		std::numeric_limits<int>::max() - 1, std::numeric_limits<int>::max()}) {
		CAPTURE(width);
		CHECK_THROWS_AS(column.width(width), ArgumentError);
		CHECK(column.width() == 7);
		CHECK(column.rawWidth() == 7);
		CHECK_THROWS_AS(column.setCalculatedWidth(width), ArgumentError);
		CHECK(column.width() == 7);
	}
	CHECK_THROWS_AS(ColumnFormat("value").width(-2).formatAligned("x"), ArgumentError);
	CHECK_NOTHROW(column.width(std::numeric_limits<int>::max() - 2));
	CHECK(column.width() == std::numeric_limits<int>::max() - 2);
	CHECK_NOTHROW(column.setCalculatedWidth(std::numeric_limits<int>::max() - 2));
}

TEST_CASE("Column automatic width is safe alone and does not shadow later explicit widths",
	"[ConsolePrinter][ColumnFormatContract]") {
	ColumnFormat column("value");
	column.autoWidth();
	CHECK(column.formatAligned("x") == "x");
	CHECK(column.formatAligned("") == "");
	column.setCalculatedWidth(6);
	CHECK(column.formatAligned("x") == "     x");
	column.width(3);
	CHECK(column.width() == 3);
	CHECK(column.formatAligned("x") == "  x");
	column.setCalculatedWidth(8);
	column.autoWidth();
	CHECK(column.width() == ColumnFormat::AUTO_WIDTH);
	CHECK(column.formatAligned("x") == "x");
	column.width(0);
	CHECK(column.formatAligned("long content") == "long content");
	column.setCalculatedWidth(5);
	column.setCalculatedWidth(ColumnFormat::AUTO_WIDTH);
	CHECK(column.width() == 0);
}

TEST_CASE("Column integer formatting rounds finite floating scalars halfway away from zero",
	"[ConsolePrinter][ColumnFormatContract][PrintFormatting]") {
	const auto column = ColumnFormat("value").format(FormatType::Integer).precision(9);
	for (const auto& sample : {std::pair<Real, std::string>{Real(1.25), "1"},
		{Real(1.5), "2"}, {Real(2.5), "3"}, {Real(-1.25), "-1"},
		{Real(-1.5), "-2"}, {Real(-2.5), "-3"}, {Real(0.25), "0"},
		{Real(-0.25), "-0"}, {Real(0), "0"}, {Real(-0.0), "-0"}}) {
		CAPTURE(sample.first);
		CHECK(column.formatValue(sample.first) == sample.second);
		CHECK(PrintFormatting::Text(sample.first, 2, PrintFormatting::Notation::Integer) == sample.second);
	}
	CHECK(column.formatValue(1.5f) == "2");
	CHECK(column.formatValue(-2.5L) == "-3");
	for (int value = -128; value <= 128; ++value) {
		CHECK(column.formatValue(Real(value)) == std::to_string(value));
		CHECK(column.formatValue(Real(value) + Real(0.25)) == std::to_string(value));
		if (value != 0)
			CHECK(column.formatValue(Real(value) - Real(0.25)) == std::to_string(value));
	}
}

TEST_CASE("Integer formatting does not narrow floating values to an integer storage type",
	"[ConsolePrinter][ColumnFormatContract][PrintFormatting]") {
	for (Real value : {std::ldexp(Real(1), 80), -std::ldexp(Real(1), 80),
		std::numeric_limits<Real>::max(), std::numeric_limits<Real>::lowest()}) {
		std::ostringstream expected;
		expected << std::fixed << std::setprecision(0) << value;
		const auto text = PrintFormatting::Text(value, 3, PrintFormatting::Notation::Integer);
		CHECK(text == expected.str());
		CHECK(text.find('.') == std::string::npos);
		CHECK(text.find_first_of("eE") == std::string::npos);
	}
	for (Real value : {std::numeric_limits<Real>::infinity(), -std::numeric_limits<Real>::infinity(),
		std::numeric_limits<Real>::quiet_NaN()}) {
		std::ostringstream expected;
		expected << std::fixed << std::setprecision(0) << value;
		CHECK(PrintFormatting::Text(value, 3, PrintFormatting::Notation::Integer) == expected.str());
	}
}

TEST_CASE("Special numeric formats reject unsupported types while preserving other format semantics",
	"[ConsolePrinter][ColumnFormatContract][PrintFormatting]") {
	const auto integer = ColumnFormat("value").format(FormatType::Integer);
	const auto hexadecimal = ColumnFormat("value").format(FormatType::Hexadecimal);
	CHECK_THROWS_AS(integer.formatValue(std::string("12")), ArgumentError);
	CHECK_THROWS_AS(integer.formatValue(Complex(1, 2)), ArgumentError);
	CHECK_THROWS_AS(hexadecimal.formatValue(Real(1.25)), ArgumentError);
	CHECK_THROWS_AS(hexadecimal.formatValue(Complex(1, 2)), ArgumentError);
	CHECK_THROWS_AS(hexadecimal.formatValue(std::string("ff")), ArgumentError);
	CHECK(integer.formatValue(std::numeric_limits<std::int64_t>::min()) ==
		std::to_string(std::numeric_limits<std::int64_t>::min()));
	CHECK(integer.formatValue(std::numeric_limits<std::uint64_t>::max()) ==
		std::to_string(std::numeric_limits<std::uint64_t>::max()));
	CHECK(hexadecimal.formatValue(255) == "ff");
	CHECK(integer.formatValue(true) == "1");
	CHECK(hexadecimal.formatValue(false) == "0");
	CHECK(integer.formatValue('A') == "65");
	CHECK(hexadecimal.formatValue('A') == "41");
	CHECK(integer.formatValue(static_cast<unsigned char>(255)) == "255");
	CHECK(hexadecimal.formatValue(static_cast<unsigned char>(255)) == "ff");
	CHECK(ColumnFormat("value").format(FormatType::Fixed).precision(2).formatValue(Real(1.25)) == "1.25");
	CHECK(ColumnFormat("value").format(FormatType::Scientific).precision(2).formatValue(Real(1.25)) == "1.25e+00");
	CHECK(ColumnFormat("value").formatValue(std::string("text")) == "text");
}

TEST_CASE("Generic and borrowed vector tables share rounded values auto widths and every border style",
	"[ConsolePrinter][ColumnFormatContract]") {
	const Vector<Real> values{Real(1.5), Real(-2.5), Real(100.25)};
	for (const auto border : {BorderStyle::None, BorderStyle::Simple, BorderStyle::Markdown,
		BorderStyle::Rounded, BorderStyle::Double, BorderStyle::Bold})
		for (const auto align : {Alignment::Left, Alignment::Right, Alignment::Center}) {
			CAPTURE(border, align);
			const auto format = ColumnFormat("v").autoWidth().format(FormatType::Integer).align(align);
			VectorTablePrinter vectorTable;
			vectorTable.addVector(format, values);
			TablePrinter<Real, Real> generic(format, {});
			for (Real value : values) generic.addRow(value, {});
			const auto style = TableStyle().border(border);
			vectorTable.style(style);
			generic.style(style);
			const auto text = vectorTable.toString();
			CHECK(text == generic.toString());
			CHECK(text.find("2") != std::string::npos);
			CHECK(text.find("-3") != std::string::npos);
			CHECK(text.find("100") != std::string::npos);
			CHECK(text.find('.') == std::string::npos);
		}
}

TEST_CASE("Alignment preserves Unicode content and computes exact nonnegative padding",
	"[ConsolePrinter][ColumnFormatContract]") {
	const std::string content = "\xCE\xB1" "b";
	for (int width = 0; width <= 20; ++width)
		for (const auto alignment : {Alignment::Left, Alignment::Right, Alignment::Center}) {
			const auto column = ColumnFormat("value").width(width).align(alignment);
			const int padding = std::max(0, width - 2);
			const int left = alignment == Alignment::Left ? 0 :
				alignment == Alignment::Right ? padding : padding / 2;
			const int right = padding - left;
			CHECK(column.formatAligned(content) ==
				std::string(static_cast<size_t>(left), ' ') + content +
				std::string(static_cast<size_t>(right), ' '));
		}
}

TEST_CASE("Border repetition rejects impossible sizes without overflowing or allocating",
	"[ConsolePrinter][ColumnFormatContract]") {
	const auto maximum = std::numeric_limits<size_t>::max();
	CHECK_THROWS_AS(Utf8::repeatString("ab", maximum), ResourceLimitError);
	CHECK_THROWS_AS(Utf8::repeatString("\xE2\x94\x80", maximum / 2), ResourceLimitError);
	CHECK_THROWS_AS(Utf8::repeatString("-", std::string{}.max_size() + 1), ResourceLimitError);
	CHECK(Utf8::repeatString("", maximum).empty());
	CHECK(Utf8::repeatString("abc", 0).empty());
	CHECK(Utf8::repeatString("\xE2\x94\x80", 2) == "\xE2\x94\x80\xE2\x94\x80");
}

TEST_CASE("Integer rounding distinguishes halfway values from adjacent representable floats",
	"[ConsolePrinter][ColumnFormatContract][PrintFormatting]") {
	for (Real tie : {Real(0.5), Real(1.5), Real(2.5), Real(100.5)}) {
		const auto below = std::nextafter(tie, Real(0));
		const auto above = std::nextafter(tie, std::numeric_limits<Real>::infinity());
		const auto smaller = std::to_string(static_cast<int>(std::floor(tie)));
		const auto larger = std::to_string(static_cast<int>(std::ceil(tie)));
		CHECK(PrintFormatting::Text(below, 0, PrintFormatting::Notation::Integer) == smaller);
		CHECK(PrintFormatting::Text(tie, 0, PrintFormatting::Notation::Integer) == larger);
		CHECK(PrintFormatting::Text(above, 0, PrintFormatting::Notation::Integer) == larger);
		const auto negativeSmaller = smaller == "0" ? "-0" : "-" + smaller;
		CHECK(PrintFormatting::Text(-below, 0, PrintFormatting::Notation::Integer) == negativeSmaller);
		CHECK(PrintFormatting::Text(-tie, 0, PrintFormatting::Notation::Integer) == "-" + larger);
		CHECK(PrintFormatting::Text(-above, 0, PrintFormatting::Notation::Integer) == "-" + larger);
	}
}

TEST_CASE("Only console rendering applies column integer formats and export formats retain source values",
	"[ConsolePrinter][ColumnFormatContract]") {
	const auto column = ColumnFormat("v").width(0).format(FormatType::Integer);
	TablePrinter<Real, Real> table(column, {});
	table.addRow(Real(1.25), {});
	table.style(TableStyle().border(BorderStyle::None).header(false));
	CHECK(table.toString() == "1\n");
	for (auto format : {ExportFormat::CSV, ExportFormat::TSV, ExportFormat::Markdown,
		ExportFormat::LaTeX, ExportFormat::HTML}) {
		std::ostringstream output;
		table.exportTo(output, format);
		CHECK(output.str().find("1.25") != std::string::npos);
	}
}
