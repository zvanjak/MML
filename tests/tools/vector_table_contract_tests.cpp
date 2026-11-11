#include <catch2/catch_all.hpp>

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#else
#include <mml/tools/ConsolePrinter.h>
#endif

#include <sstream>
#include <limits>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace MML;

namespace {
	template<class Label, class Value, class = void>
	struct CanBorrowVector : std::false_type {};

	template<class Label, class Value>
	struct CanBorrowVector<Label, Value, std::void_t<decltype(
		std::declval<VectorTablePrinter&>().addVector(std::declval<Label>(), std::declval<Value>()))>> : std::true_type {};
}

TEST_CASE("Vector tables accept live lvalues but reject temporary vector borrowing", "[ConsolePrinter][VectorTableContract]")
{
	CHECK((CanBorrowVector<const std::string&, Vector<Real>&>::value));
	CHECK((CanBorrowVector<const std::string&, const Vector<Real>&>::value));
	CHECK((CanBorrowVector<const ColumnFormat&, Vector<Real>&>::value));
	CHECK((CanBorrowVector<const ColumnFormat&, const Vector<Real>&>::value));
	CHECK_FALSE((CanBorrowVector<const std::string&, Vector<Real>&&>::value));
	CHECK_FALSE((CanBorrowVector<const std::string&, const Vector<Real>&&>::value));
	CHECK_FALSE((CanBorrowVector<const ColumnFormat&, Vector<Real>&&>::value));
	CHECK_FALSE((CanBorrowVector<const ColumnFormat&, const Vector<Real>&&>::value));
}

TEST_CASE("Borrowed vector tables recompute live row counts on every print", "[ConsolePrinter][VectorTableContract]")
{
	Vector<Real> first{1}, second{9};
	VectorTablePrinter table;
	table.addVector(ColumnFormat("a").width(0), first);
	table.addVector(ColumnFormat("b").width(0), second);
	table.style(TableStyle().border(BorderStyle::None).header(false));
	CHECK(table.toString() == "1  9\n");
	first.push_back(2);
	first.push_back(3);
	second.resize(0);
	CHECK(table.toString() == "1  \n2  \n3  \n");
	first.resize(0);
	second = Vector<Real>{4, 5};
	CHECK(table.toString() == "  4\n  5\n");
	first.resize(0);
	second.resize(0);
	CHECK(table.toString().empty());
	second.push_back(8);
	const auto copied = table;
	VectorTablePrinter assigned;
	assigned = table;
	CHECK(copied.toString() == "  8\n");
	CHECK(assigned.toString() == "  8\n");
	second[0] = 7;
	CHECK(table.toString() == "  7\n");
	CHECK(copied.toString() == table.toString());
	CHECK(assigned.toString() == table.toString());
}

TEST_CASE("Vector automatic widths match generic table formatting and borders", "[ConsolePrinter][VectorTableContract]")
{
	const Vector<Real> first{Real(1.25), Real(-123.5)}, second{Real(3.5), Real(16)};
	for (const auto border : {BorderStyle::None, BorderStyle::Simple, BorderStyle::Markdown,
		BorderStyle::Rounded, BorderStyle::Double, BorderStyle::Bold}) {
		for (const auto alignment : {Alignment::Left, Alignment::Center, Alignment::Right}) {
			for (bool header : {false, true}) {
				CAPTURE(border, alignment, header);
				const auto a = ColumnFormat("a").autoWidth().precision(2).format(FormatType::Fixed).align(alignment);
				const auto b = ColumnFormat("long heading").autoWidth().precision(3).format(FormatType::Scientific).align(alignment);
				VectorTablePrinter vectors;
				vectors.addVector(a, first);
				vectors.addVector(b, second);
				TablePrinter<Real, Real> generic(a, {b});
				for (int i = 0; i < first.size(); ++i) generic.addRow(first[i], {second[i]});
				const auto style = TableStyle().border(border).header(header);
				vectors.style(style);
				generic.style(style);
				std::ostringstream expected;
				generic.print(expected);
				std::string actual;
				REQUIRE_NOTHROW(actual = vectors.toString());
				CHECK(actual == expected.str());
			}
		}
	}
}

TEST_CASE("Automatic vector widths track changing data without retaining stale calculations", "[ConsolePrinter][VectorTableContract]")
{
	Vector<Real> values{1};
	VectorTablePrinter table;
	auto format = ColumnFormat("v").autoWidth().format(FormatType::Fixed).precision(2);
	format.setCalculatedWidth(50);
	table.addVector(format, values);
	std::string text;
	REQUIRE_NOTHROW(text = table.toString());
	CHECK(text == "+------+\n|    v |\n+------+\n| 1.00 |\n+------+\n");
	values = Vector<Real>{Real(-1234.5)};
	CHECK(table.toString() == "+----------+\n|        v |\n+----------+\n| -1234.50 |\n+----------+\n");
	values = Vector<Real>{2};
	CHECK(table.toString() == "+------+\n|    v |\n+------+\n| 2.00 |\n+------+\n");
	values.resize(0);
	CHECK(table.toString() == "+---+\n| v |\n+---+\n+---+\n");
	values.push_back(3);
	CHECK(table.toString() == "+------+\n|    v |\n+------+\n| 3.00 |\n+------+\n");
}

TEST_CASE("Automatic vector widths handle empty and unequal columns and UTF8 headings", "[ConsolePrinter][VectorTableContract]")
{
	const Vector<Real> empty;
	const Vector<Real> values{1, 2};
	VectorTablePrinter table;
	table.addVector(ColumnFormat("\xCE\xB1").autoWidth().align(Alignment::Left), empty);
	table.addVector(ColumnFormat("b").autoWidth(), values);
	std::string text;
	REQUIRE_NOTHROW(text = table.toString());
	CHECK(text == "+---+---+\n| \xCE\xB1 | b |\n+---+---+\n|   | 1 |\n|   | 2 |\n+---+---+\n");
	const Vector<Real> noValues;
	VectorTablePrinter nameless;
	nameless.addVector(ColumnFormat("").autoWidth(), noValues);
	CHECK(nameless.toString() == "+--+\n|  |\n+--+\n+--+\n");
}

TEST_CASE("Pointer vector table constructors reject nulls and preserve default format behavior", "[ConsolePrinter][VectorTableContract]")
{
	Vector<Real> first{1}, second{2};
	CHECK_THROWS_AS((VectorTablePrinter({}, {nullptr})), std::invalid_argument);
	CHECK_THROWS_AS((VectorTablePrinter({}, {&first, nullptr})), std::invalid_argument);
	CHECK_THROWS_AS((VectorTablePrinter({ColumnFormat("a"), ColumnFormat("b")}, {&first, nullptr})), std::invalid_argument);
	VectorTablePrinter constructed({ColumnFormat("a").autoWidth()}, {&first, &second});
	VectorTablePrinter expected;
	expected.addVector(ColumnFormat("a").autoWidth(), first);
	expected.addVector("Column1", second);
	CHECK(constructed.toString() == expected.toString());
	second.push_back(3);
	CHECK(constructed.toString() == expected.toString());
	VectorTablePrinter extra({ColumnFormat("a").autoWidth(), ColumnFormat("unused")}, {&first});
	VectorTablePrinter single;
	single.addVector(ColumnFormat("a").autoWidth(), first);
	CHECK(extra.toString() == single.toString());
	VectorTablePrinter empty({}, {});
	CHECK(empty.toString() == VectorTablePrinter().toString());
}

TEST_CASE("Invalid vector column widths fail explicitly before changing the table", "[ConsolePrinter][VectorTableContract]")
{
	const Vector<Real> values{1};
	VectorTablePrinter table;
	table.addVector(ColumnFormat("valid").width(0), values);
	const auto original = table.toString();
	for (int width : {-2, std::numeric_limits<int>::min()}) {
		CHECK_THROWS_AS(table.addVector(ColumnFormat("invalid").width(width), values), std::invalid_argument);
		CHECK(table.toString() == original);
	}
	Vector<Real> mutableValues{2};
	CHECK_THROWS_AS((VectorTablePrinter({ColumnFormat("invalid").width(-2)}, {&mutableValues})), std::invalid_argument);
}

TEST_CASE("Vector table automatic widths preserve caller state on output and stream failure", "[ConsolePrinter][VectorTableContract]")
{
	class ThrowingBuffer : public std::streambuf {
		std::streamsize xsputn(const char*, std::streamsize) override { throw std::runtime_error("write failed"); }
		int_type overflow(int_type) override { throw std::runtime_error("write failed"); }
	};
	const Vector<Real> values{Real(1.25), 100};
	VectorTablePrinter printer;
	printer.addVector(ColumnFormat("live").autoWidth().format(FormatType::Fixed).precision(2), values);
	const auto& view = printer;
	for (bool fail : {false, true}) {
		ThrowingBuffer buffer;
		std::ostringstream successful;
		std::ostream failing(&buffer);
		auto& out = fail ? failing : successful;
		out << std::scientific << std::showpos;
		out.precision(3);
		out.width(7);
		out.fill('@');
		out.exceptions(std::ios::badbit);
		const auto flags = out.flags();
		const auto locale = out.getloc();
		if (fail) CHECK_THROWS(view.print(out));
		else CHECK_NOTHROW(view.print(out));
		CHECK(out.flags() == flags);
		CHECK(out.precision() == 3);
		CHECK(out.width() == 7);
		CHECK(out.fill() == '@');
		CHECK(out.getloc() == locale);
	}
}
