/**
 * @file console_printer_tests.cpp
 * @brief Comprehensive tests for ConsolePrinter.h
 * 
 * Tests cover:
 * - Utf8 namespace helpers (repeatString, isSingleDisplayChar)
 * - Enumerations (FormatType, Alignment, BorderStyle, ExportFormat)
 * - StreamStateGuard RAII class
 * - ColumnFormat class with builder pattern
 * - TableStyle class with builder pattern
 * - TablePrinter template class
 * - VectorTablePrinter class
 */

#include <catch2/catch_all.hpp>
#include <sstream>
#include <iomanip>
#include <string>

#include "tools/ConsolePrinter.h"
#include "base/Vector/Vector.h"

using namespace MML;

namespace MML::Tests::Tools::ConsolePrinterTests {

///////////////////////////////////////////////////////////////////////////////
//                           UTF-8 HELPERS                                   //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("Utf8::repeatString - Empty string", "[ConsolePrinter][Utf8]") {
    REQUIRE(Utf8::repeatString("", 5) == "");
    REQUIRE(Utf8::repeatString("abc", 0) == "");
}

TEST_CASE("Utf8::repeatString - ASCII characters", "[ConsolePrinter][Utf8]") {
    REQUIRE(Utf8::repeatString("-", 5) == "-----");
    REQUIRE(Utf8::repeatString("ab", 3) == "ababab");
    REQUIRE(Utf8::repeatString("x", 1) == "x");
}

TEST_CASE("Utf8::repeatString - UTF-8 multi-byte characters", "[ConsolePrinter][Utf8]") {
    // Box-drawing characters are 3 bytes each
    std::string dash = "─";  // U+2500
    std::string result = Utf8::repeatString(dash, 3);
    REQUIRE(result == "───");
    REQUIRE(result.size() == 9);  // 3 chars × 3 bytes each
}

TEST_CASE("Utf8::isSingleDisplayChar - ASCII", "[ConsolePrinter][Utf8]") {
    REQUIRE(Utf8::isSingleDisplayChar("a") == true);
    REQUIRE(Utf8::isSingleDisplayChar("X") == true);
    REQUIRE(Utf8::isSingleDisplayChar(" ") == true);
    REQUIRE(Utf8::isSingleDisplayChar("ab") == false);  // Two chars
}

TEST_CASE("Utf8::isSingleDisplayChar - UTF-8", "[ConsolePrinter][Utf8]") {
    REQUIRE(Utf8::isSingleDisplayChar("─") == true);   // 3-byte
    REQUIRE(Utf8::isSingleDisplayChar("═") == true);   // 3-byte
    REQUIRE(Utf8::isSingleDisplayChar("╔") == true);   // 3-byte
    REQUIRE(Utf8::isSingleDisplayChar("") == false);   // Empty
}

///////////////////////////////////////////////////////////////////////////////
//                           ENUMERATIONS                                    //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("FormatType - All values exist", "[ConsolePrinter][Enums]") {
    REQUIRE(static_cast<int>(FormatType::General) >= 0);
    REQUIRE(static_cast<int>(FormatType::Fixed) >= 0);
    REQUIRE(static_cast<int>(FormatType::Scientific) >= 0);
    REQUIRE(static_cast<int>(FormatType::Hexadecimal) >= 0);
    REQUIRE(static_cast<int>(FormatType::Integer) >= 0);
}

TEST_CASE("Alignment - All values exist", "[ConsolePrinter][Enums]") {
    REQUIRE(static_cast<int>(Alignment::Left) >= 0);
    REQUIRE(static_cast<int>(Alignment::Right) >= 0);
    REQUIRE(static_cast<int>(Alignment::Center) >= 0);
}

TEST_CASE("BorderStyle - All values exist", "[ConsolePrinter][Enums]") {
    REQUIRE(static_cast<int>(BorderStyle::None) >= 0);
    REQUIRE(static_cast<int>(BorderStyle::Simple) >= 0);
    REQUIRE(static_cast<int>(BorderStyle::Markdown) >= 0);
    REQUIRE(static_cast<int>(BorderStyle::Rounded) >= 0);
    REQUIRE(static_cast<int>(BorderStyle::Double) >= 0);
    REQUIRE(static_cast<int>(BorderStyle::Bold) >= 0);
}

TEST_CASE("ExportFormat - All values exist", "[ConsolePrinter][Enums]") {
    REQUIRE(static_cast<int>(ExportFormat::Console) >= 0);
    REQUIRE(static_cast<int>(ExportFormat::CSV) >= 0);
    REQUIRE(static_cast<int>(ExportFormat::TSV) >= 0);
    REQUIRE(static_cast<int>(ExportFormat::Markdown) >= 0);
    REQUIRE(static_cast<int>(ExportFormat::LaTeX) >= 0);
    REQUIRE(static_cast<int>(ExportFormat::HTML) >= 0);
}

///////////////////////////////////////////////////////////////////////////////
//                        STREAM STATE GUARD                                 //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("StreamStateGuard - Preserves precision", "[ConsolePrinter][StreamStateGuard]") {
    std::ostringstream oss;
    oss << std::setprecision(3);
    
    {
        StreamStateGuard guard(oss);
        oss << std::setprecision(10);
        REQUIRE(oss.precision() == 10);
    }
    
    REQUIRE(oss.precision() == 3);  // Restored
}

TEST_CASE("StreamStateGuard - Preserves flags", "[ConsolePrinter][StreamStateGuard]") {
    std::ostringstream oss;
    oss << std::fixed;
    auto originalFlags = oss.flags();
    
    {
        StreamStateGuard guard(oss);
        oss << std::scientific;
        REQUIRE(oss.flags() != originalFlags);
    }
    
    REQUIRE(oss.flags() == originalFlags);  // Restored
}

TEST_CASE("StreamStateGuard - Preserves width", "[ConsolePrinter][StreamStateGuard]") {
    std::ostringstream oss;
    oss << std::setw(15);
    
    {
        StreamStateGuard guard(oss);
        oss << std::setw(5);
    }
    
    // Width is consumed after first output, but guard saves initial state
    REQUIRE(oss.width() == 15);
}

///////////////////////////////////////////////////////////////////////////////
//                          COLUMN FORMAT                                    //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("ColumnFormat - Default constructor", "[ConsolePrinter][ColumnFormat]") {
    ColumnFormat col("TestColumn");
    
    REQUIRE(col.name() == "TestColumn");
    REQUIRE(col.width() == 12);  // Default width
    REQUIRE(col.precision() == 6);  // Default precision
    REQUIRE(col.formatType() == FormatType::General);
    REQUIRE(col.alignment() == Alignment::Right);
}

TEST_CASE("ColumnFormat - Builder pattern width", "[ConsolePrinter][ColumnFormat]") {
    ColumnFormat col("Col");
    col.width(20);
    REQUIRE(col.width() == 20);
}

TEST_CASE("ColumnFormat - Builder pattern precision", "[ConsolePrinter][ColumnFormat]") {
    ColumnFormat col("Col");
    col.precision(10);
    REQUIRE(col.precision() == 10);
}

TEST_CASE("ColumnFormat - Builder pattern format", "[ConsolePrinter][ColumnFormat]") {
    ColumnFormat col("Col");
    col.format(FormatType::Scientific);
    REQUIRE(col.formatType() == FormatType::Scientific);
}

TEST_CASE("ColumnFormat - Builder pattern alignment", "[ConsolePrinter][ColumnFormat]") {
    ColumnFormat col("Col");
    col.align(Alignment::Left);
    REQUIRE(col.alignment() == Alignment::Left);
}

TEST_CASE("ColumnFormat - Builder pattern chaining", "[ConsolePrinter][ColumnFormat]") {
    ColumnFormat col("Value");
    col.width(15).precision(4).format(FormatType::Fixed).align(Alignment::Center);
    
    REQUIRE(col.width() == 15);
    REQUIRE(col.precision() == 4);
    REQUIRE(col.formatType() == FormatType::Fixed);
    REQUIRE(col.alignment() == Alignment::Center);
}

TEST_CASE("ColumnFormat - Auto width", "[ConsolePrinter][ColumnFormat]") {
    ColumnFormat col("Col");
    col.autoWidth();
    REQUIRE(col.rawWidth() == ColumnFormat::AUTO_WIDTH);
}

TEST_CASE("ColumnFormat - formatValue fixed", "[ConsolePrinter][ColumnFormat]") {
    ColumnFormat col("Val");
    col.format(FormatType::Fixed).precision(2);
    
    std::string result = col.formatValue(3.14159);
    REQUIRE(result == "3.14");
}

TEST_CASE("ColumnFormat - formatValue scientific", "[ConsolePrinter][ColumnFormat]") {
    ColumnFormat col("Val");
    col.format(FormatType::Scientific).precision(2);
    
    std::string result = col.formatValue(12345.0);
    REQUIRE(result.find('e') != std::string::npos);  // Contains exponent
}

TEST_CASE("ColumnFormat - formatAligned left", "[ConsolePrinter][ColumnFormat]") {
    ColumnFormat col("Col");
    col.width(10).align(Alignment::Left);
    
    std::string result = col.formatAligned("Hi");
    REQUIRE(result == "Hi        ");
    REQUIRE(result.length() == 10);
}

TEST_CASE("ColumnFormat - formatAligned right", "[ConsolePrinter][ColumnFormat]") {
    ColumnFormat col("Col");
    col.width(10).align(Alignment::Right);
    
    std::string result = col.formatAligned("Hi");
    REQUIRE(result == "        Hi");
    REQUIRE(result.length() == 10);
}

TEST_CASE("ColumnFormat - formatAligned center", "[ConsolePrinter][ColumnFormat]") {
    ColumnFormat col("Col");
    col.width(10).align(Alignment::Center);
    
    std::string result = col.formatAligned("Hi");
    REQUIRE(result == "    Hi    ");
    REQUIRE(result.length() == 10);
}

TEST_CASE("ColumnFormat - formatAligned content exceeds width", "[ConsolePrinter][ColumnFormat]") {
    ColumnFormat col("Col");
    col.width(5);
    
    std::string result = col.formatAligned("LongContent");
    REQUIRE(result == "LongContent");  // Not truncated
}

TEST_CASE("ColumnFormat - Legacy constructor", "[ConsolePrinter][ColumnFormat]") {
    ColumnFormat col("Name", 15, 4, 'F');
    
    REQUIRE(col.name() == "Name");
    REQUIRE(col.width() == 15);
    REQUIRE(col.precision() == 4);
    REQUIRE(col.formatType() == FormatType::Fixed);
}

///////////////////////////////////////////////////////////////////////////////
//                           TABLE STYLE                                     //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("TableStyle - Default values", "[ConsolePrinter][TableStyle]") {
    TableStyle style;
    
    REQUIRE(style.borderStyle() == BorderStyle::Simple);
    REQUIRE(style.showHeader() == true);
    REQUIRE(style.showRowSeparators() == false);
    REQUIRE(style.compactMode() == false);
}

TEST_CASE("TableStyle - Builder pattern", "[ConsolePrinter][TableStyle]") {
    TableStyle style;
    style.border(BorderStyle::Double)
         .header(false)
         .rowSeparators(true)
         .compact(true);
    
    REQUIRE(style.borderStyle() == BorderStyle::Double);
    REQUIRE(style.showHeader() == false);
    REQUIRE(style.showRowSeparators() == true);
    REQUIRE(style.compactMode() == true);
}

TEST_CASE("TableStyle - getChars None", "[ConsolePrinter][TableStyle]") {
    TableStyle style;
    style.border(BorderStyle::None);
    auto chars = style.getChars();
    
    REQUIRE(chars.topLeft == "");
    REQUIRE(chars.horizontal == "");
    REQUIRE(chars.vertical == "");
}

TEST_CASE("TableStyle - getChars Simple", "[ConsolePrinter][TableStyle]") {
    TableStyle style;
    style.border(BorderStyle::Simple);
    auto chars = style.getChars();
    
    REQUIRE(chars.topLeft == "+");
    REQUIRE(chars.horizontal == "-");
    REQUIRE(chars.vertical == "|");
}

TEST_CASE("TableStyle - getChars Rounded", "[ConsolePrinter][TableStyle]") {
    TableStyle style;
    style.border(BorderStyle::Rounded);
    auto chars = style.getChars();
    
    REQUIRE(chars.topLeft == "╭");
    REQUIRE(chars.topRight == "╮");
    REQUIRE(chars.horizontal == "─");
}

TEST_CASE("TableStyle - getChars Double", "[ConsolePrinter][TableStyle]") {
    TableStyle style;
    style.border(BorderStyle::Double);
    auto chars = style.getChars();
    
    REQUIRE(chars.topLeft == "╔");
    REQUIRE(chars.horizontal == "═");
    REQUIRE(chars.vertical == "║");
}

TEST_CASE("TableStyle - getChars Bold", "[ConsolePrinter][TableStyle]") {
    TableStyle style;
    style.border(BorderStyle::Bold);
    auto chars = style.getChars();
    
    REQUIRE(chars.topLeft == "┏");
    REQUIRE(chars.horizontal == "━");
    REQUIRE(chars.vertical == "┃");
}

///////////////////////////////////////////////////////////////////////////////
//                          TABLE PRINTER                                    //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("TablePrinter - Construction with column names", "[ConsolePrinter][TablePrinter]") {
    TablePrinter<std::string, double> table("Index", {"Value1", "Value2"});
    
    REQUIRE(table.rowCount() == 0);
    REQUIRE(table.columnCount() == 2);
}

TEST_CASE("TablePrinter - Add row", "[ConsolePrinter][TablePrinter]") {
    TablePrinter<std::string, double> table("X", {"Y", "Z"});
    table.addRow("A", {1.0, 2.0});
    table.addRow("B", {3.0, 4.0});
    
    REQUIRE(table.rowCount() == 2);
}

TEST_CASE("TablePrinter - Add row wrong size throws", "[ConsolePrinter][TablePrinter]") {
    TablePrinter<std::string, double> table("X", {"Y", "Z"});
    
    REQUIRE_THROWS_AS(table.addRow("A", {1.0}), VectorDimensionError);
}

TEST_CASE("TablePrinter - Clear", "[ConsolePrinter][TablePrinter]") {
    TablePrinter<std::string, double> table("X", {"Y"});
    table.addRow("A", {1.0});
    table.addRow("B", {2.0});
    table.clear();
    
    REQUIRE(table.rowCount() == 0);
}

TEST_CASE("TablePrinter - toString produces output", "[ConsolePrinter][TablePrinter]") {
    TablePrinter<std::string, double> table("Item", {"Price"});
    table.addRow("Apple", {1.50});
    
    std::string output = table.toString();
    REQUIRE(output.find("Item") != std::string::npos);
    REQUIRE(output.find("Price") != std::string::npos);
    REQUIRE(output.find("Apple") != std::string::npos);
}

TEST_CASE("TablePrinter - Print to stream", "[ConsolePrinter][TablePrinter]") {
    TablePrinter<int, double> table("N", {"Square"});
    table.addRow(1, {1.0});
    table.addRow(2, {4.0});
    table.addRow(3, {9.0});
    
    std::ostringstream oss;
    table.print(oss);
    
    std::string output = oss.str();
    REQUIRE(output.length() > 0);
    REQUIRE(output.find("Square") != std::string::npos);
}

TEST_CASE("TablePrinter - Style configuration", "[ConsolePrinter][TablePrinter]") {
    TablePrinter<std::string, double> table("X", {"Y"});
    table.style(TableStyle().border(BorderStyle::None).header(false));
    table.addRow("A", {1.0});
    
    std::string output = table.toString();
    // Without borders and headers, output should be minimal
    REQUIRE(output.find("+") == std::string::npos);  // No border chars
}

TEST_CASE("TablePrinter - Export CSV", "[ConsolePrinter][TablePrinter]") {
    TablePrinter<std::string, double> table("Name", {"Value"});
    table.addRow("Test", {42.5});
    
    std::ostringstream oss;
    table.exportTo(oss, ExportFormat::CSV);
    
    std::string output = oss.str();
    REQUIRE(output.find("Name,Value") != std::string::npos);
    REQUIRE(output.find("Test,42.5") != std::string::npos);
}

TEST_CASE("TablePrinter - Export TSV", "[ConsolePrinter][TablePrinter]") {
    TablePrinter<std::string, int> table("Key", {"Count"});
    table.addRow("Items", {10});
    
    std::ostringstream oss;
    table.exportTo(oss, ExportFormat::TSV);
    
    std::string output = oss.str();
    REQUIRE(output.find("Key\tCount") != std::string::npos);
    REQUIRE(output.find("Items\t10") != std::string::npos);
}

TEST_CASE("TablePrinter - Export Markdown", "[ConsolePrinter][TablePrinter]") {
    TablePrinter<std::string, double> table("Col1", {"Col2"});
    table.addRow("Row1", {1.0});
    
    std::ostringstream oss;
    table.exportTo(oss, ExportFormat::Markdown);
    
    std::string output = oss.str();
    REQUIRE(output.find("| Col1 | Col2 |") != std::string::npos);
    REQUIRE(output.find("|---") != std::string::npos);  // Separator line
}

TEST_CASE("TablePrinter - Export LaTeX", "[ConsolePrinter][TablePrinter]") {
    TablePrinter<std::string, int> table("A", {"B"});
    table.addRow("x", {1});
    
    std::ostringstream oss;
    table.exportTo(oss, ExportFormat::LaTeX);
    
    std::string output = oss.str();
    REQUIRE(output.find("\\begin{tabular}") != std::string::npos);
    REQUIRE(output.find("\\end{tabular}") != std::string::npos);
    REQUIRE(output.find("\\hline") != std::string::npos);
    REQUIRE(output.find(" & ") != std::string::npos);  // Column separator
}

TEST_CASE("TablePrinter - Export HTML", "[ConsolePrinter][TablePrinter]") {
    TablePrinter<std::string, int> table("Header", {"Data"});
    table.addRow("row", {42});
    
    std::ostringstream oss;
    table.exportTo(oss, ExportFormat::HTML);
    
    std::string output = oss.str();
    REQUIRE(output.find("<table") != std::string::npos);
    REQUIRE(output.find("</table>") != std::string::npos);
    REQUIRE(output.find("<th>") != std::string::npos);
    REQUIRE(output.find("<td>") != std::string::npos);
}

TEST_CASE("TablePrinter - Reserve does not affect row count", "[ConsolePrinter][TablePrinter]") {
    TablePrinter<int, double> table("N", {"Val"});
    table.reserve(100);
    
    REQUIRE(table.rowCount() == 0);  // Reserve doesn't add rows
}

///////////////////////////////////////////////////////////////////////////////
//                       VECTOR TABLE PRINTER                                //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("VectorTablePrinter - Default construction", "[ConsolePrinter][VectorTablePrinter]") {
    VectorTablePrinter printer;
    std::string output = printer.toString();
    // Empty printer should still produce valid (possibly minimal) output
    REQUIRE_NOTHROW(printer.toString());
}

TEST_CASE("VectorTablePrinter - Add vectors", "[ConsolePrinter][VectorTablePrinter]") {
    Vector<Real> v1({1.0, 2.0, 3.0});
    Vector<Real> v2({4.0, 5.0, 6.0});
    
    VectorTablePrinter printer;
    printer.addVector("First", v1);
    printer.addVector("Second", v2);
    
    std::string output = printer.toString();
    REQUIRE(output.find("First") != std::string::npos);
    REQUIRE(output.find("Second") != std::string::npos);
}

TEST_CASE("VectorTablePrinter - Different length vectors", "[ConsolePrinter][VectorTablePrinter]") {
    Vector<Real> v1({1.0, 2.0});
    Vector<Real> v2({4.0, 5.0, 6.0, 7.0});
    
    VectorTablePrinter printer;
    printer.addVector("Short", v1);
    printer.addVector("Long", v2);
    
    // Should handle mismatched lengths gracefully
    std::string output = printer.toString();
    REQUIRE(output.length() > 0);
}

TEST_CASE("VectorTablePrinter - With ColumnFormat", "[ConsolePrinter][VectorTablePrinter]") {
    Vector<Real> v({1.23456, 2.34567});
    
    VectorTablePrinter printer;
    ColumnFormat fmt("Values");
    fmt.width(15).precision(3).format(FormatType::Fixed);
    printer.addVector(fmt, v);
    
    std::string output = printer.toString();
    REQUIRE(output.find("Values") != std::string::npos);
    REQUIRE(output.find("1.235") != std::string::npos);  // Rounded to 3 decimals
}

TEST_CASE("VectorTablePrinter - Style configuration", "[ConsolePrinter][VectorTablePrinter]") {
    Vector<Real> v({1.0, 2.0});
    
    VectorTablePrinter printer;
    printer.addVector("Col", v);
    printer.style(TableStyle().border(BorderStyle::Double));
    
    std::string output = printer.toString();
    REQUIRE(output.find("═") != std::string::npos);  // Double border char
}

TEST_CASE("VectorTablePrinter - Print to stream", "[ConsolePrinter][VectorTablePrinter]") {
    Vector<Real> v({10.0, 20.0, 30.0});
    
    VectorTablePrinter printer;
    printer.addVector("Numbers", v);
    
    std::ostringstream oss;
    printer.print(oss);
    
    std::string output = oss.str();
    REQUIRE(output.find("Numbers") != std::string::npos);
    REQUIRE(output.find("10") != std::string::npos);
}

///////////////////////////////////////////////////////////////////////////////
//                          INTEGRATION TESTS                                //
///////////////////////////////////////////////////////////////////////////////

TEST_CASE("TablePrinter - Complete workflow", "[ConsolePrinter][Integration]") {
    // Create table with custom formatting
    TablePrinter<double, double> table("x", {"sin(x)", "cos(x)"});
    
    table.style(TableStyle().border(BorderStyle::Rounded));
    table.tagFormat(ColumnFormat("x").width(8).precision(4).format(FormatType::Fixed));
    table.columnFormat(0, ColumnFormat("sin(x)").width(12).precision(6).format(FormatType::Fixed));
    table.columnFormat(1, ColumnFormat("cos(x)").width(12).precision(6).format(FormatType::Fixed));
    
    // Add data
    for (double x = 0.0; x <= 1.0; x += 0.25) {
        table.addRow(x, {std::sin(x), std::cos(x)});
    }
    
    // Export to different formats
    std::ostringstream console, csv, markdown;
    table.exportTo(console, ExportFormat::Console);
    table.exportTo(csv, ExportFormat::CSV);
    table.exportTo(markdown, ExportFormat::Markdown);
    
    REQUIRE(console.str().find("sin(x)") != std::string::npos);
    REQUIRE(csv.str().find(",") != std::string::npos);
    REQUIRE(markdown.str().find("|") != std::string::npos);
}

TEST_CASE("TablePrinter - Numeric types", "[ConsolePrinter][Integration]") {
    // Test with different numeric types
    TablePrinter<int, double> intTable("N", {"Value"});
    intTable.addRow(1, {1.5});
    intTable.addRow(2, {2.5});
    REQUIRE(intTable.rowCount() == 2);
    
    TablePrinter<double, int> doubleTable("X", {"Count"});
    doubleTable.addRow(1.5, {10});
    doubleTable.addRow(2.5, {20});
    REQUIRE(doubleTable.rowCount() == 2);
}

TEST_CASE("All border styles produce valid output", "[ConsolePrinter][Integration]") {
    TablePrinter<std::string, int> table("A", {"B"});
    table.addRow("x", {1});
    
    std::vector<BorderStyle> styles = {
        BorderStyle::None,
        BorderStyle::Simple,
        BorderStyle::Markdown,
        BorderStyle::Rounded,
        BorderStyle::Double,
        BorderStyle::Bold
    };
    
    for (auto style : styles) {
        table.style(TableStyle().border(style));
        std::string output = table.toString();
        REQUIRE(output.length() > 0);
    }
}

TEST_CASE("All export formats produce valid output", "[ConsolePrinter][Integration]") {
    TablePrinter<std::string, double> table("Key", {"Value"});
    table.addRow("test", {42.0});
    
    std::vector<ExportFormat> formats = {
        ExportFormat::Console,
        ExportFormat::CSV,
        ExportFormat::TSV,
        ExportFormat::Markdown,
        ExportFormat::LaTeX,
        ExportFormat::HTML
    };
    
    for (auto format : formats) {
        std::ostringstream oss;
        table.exportTo(oss, format);
        REQUIRE(oss.str().length() > 0);
    }
}

///////////////////////////       RFC 4180 CSV TESTS       ///////////////////////////

TEST_CASE("Csv::needsQuoting - Plain text", "[ConsolePrinter][Csv]") {
    REQUIRE_FALSE(Csv::needsQuoting(""));
    REQUIRE_FALSE(Csv::needsQuoting("hello"));
    REQUIRE_FALSE(Csv::needsQuoting("123"));
    REQUIRE_FALSE(Csv::needsQuoting("hello world"));
}

TEST_CASE("Csv::needsQuoting - Contains comma", "[ConsolePrinter][Csv]") {
    REQUIRE(Csv::needsQuoting(","));
    REQUIRE(Csv::needsQuoting("a,b"));
    REQUIRE(Csv::needsQuoting("one, two, three"));
}

TEST_CASE("Csv::needsQuoting - Contains double quote", "[ConsolePrinter][Csv]") {
    REQUIRE(Csv::needsQuoting("\""));
    REQUIRE(Csv::needsQuoting("He said \"hi\""));
    REQUIRE(Csv::needsQuoting("quote\"midway"));
}

TEST_CASE("Csv::needsQuoting - Contains newlines", "[ConsolePrinter][Csv]") {
    REQUIRE(Csv::needsQuoting("\n"));
    REQUIRE(Csv::needsQuoting("line1\nline2"));
    REQUIRE(Csv::needsQuoting("\r"));
    REQUIRE(Csv::needsQuoting("line1\r\nline2"));
}

TEST_CASE("Csv::escape - Plain text unchanged", "[ConsolePrinter][Csv]") {
    REQUIRE(Csv::escape("") == "");
    REQUIRE(Csv::escape("hello") == "hello");
    REQUIRE(Csv::escape("123.456") == "123.456");
    REQUIRE(Csv::escape("hello world") == "hello world");
}

TEST_CASE("Csv::escape - Comma triggers quoting", "[ConsolePrinter][Csv]") {
    REQUIRE(Csv::escape("a,b") == "\"a,b\"");
    REQUIRE(Csv::escape(",") == "\",\"");
    REQUIRE(Csv::escape("one, two") == "\"one, two\"");
}

TEST_CASE("Csv::escape - Double quotes are doubled", "[ConsolePrinter][Csv]") {
    REQUIRE(Csv::escape("\"") == "\"\"\"\"");  // " becomes ""
    REQUIRE(Csv::escape("He said \"Hi\"") == "\"He said \"\"Hi\"\"\"");
    REQUIRE(Csv::escape("a\"b") == "\"a\"\"b\"");
}

TEST_CASE("Csv::escape - Newlines trigger quoting", "[ConsolePrinter][Csv]") {
    REQUIRE(Csv::escape("a\nb") == "\"a\nb\"");
    REQUIRE(Csv::escape("line1\r\nline2") == "\"line1\r\nline2\"");
}

TEST_CASE("Csv::escape - Combined special characters", "[ConsolePrinter][Csv]") {
    // Comma + quote + newline
    REQUIRE(Csv::escape("a,\"b\"\nc") == "\"a,\"\"b\"\"\nc\"");
}

TEST_CASE("Csv::unescape - Plain text unchanged", "[ConsolePrinter][Csv]") {
    REQUIRE(Csv::unescape("") == "");
    REQUIRE(Csv::unescape("hello") == "hello");
    REQUIRE(Csv::unescape("123.456") == "123.456");
}

TEST_CASE("Csv::unescape - Removes quotes", "[ConsolePrinter][Csv]") {
    REQUIRE(Csv::unescape("\"a,b\"") == "a,b");
    REQUIRE(Csv::unescape("\"hello\"") == "hello");
}

TEST_CASE("Csv::unescape - Undoubles double quotes", "[ConsolePrinter][Csv]") {
    REQUIRE(Csv::unescape("\"a\"\"b\"") == "a\"b");
    REQUIRE(Csv::unescape("\"He said \"\"Hi\"\"\"") == "He said \"Hi\"");
}

TEST_CASE("Csv::escape and unescape are inverse operations", "[ConsolePrinter][Csv]") {
    std::vector<std::string> testCases = {
        "",
        "hello",
        "hello world",
        "123.456",
        "a,b",
        "He said \"Hi\"",
        "line1\nline2",
        "complex, \"quoted\", and\nmulti-line"
    };
    
    for (const auto& original : testCases) {
        std::string escaped = Csv::escape(original);
        std::string unescaped = Csv::unescape(escaped);
        REQUIRE(unescaped == original);
    }
}

TEST_CASE("TablePrinter - CSV export with special characters", "[ConsolePrinter][Csv][Integration]") {
    TablePrinter<std::string, double> table("Name", {"Value", "Description"});
    table.addRow("Normal", {1.0, 0.0});
    table.addRow("With, Comma", {2.0, 0.0});  // Row tag with comma
    
    std::ostringstream oss;
    table.exportTo(oss, ExportFormat::CSV);
    
    std::string output = oss.str();
    // Should properly escape the row tag with comma
    REQUIRE(output.find("\"With, Comma\"") != std::string::npos);
}

TEST_CASE("TablePrinter - CSV export with quoted column names", "[ConsolePrinter][Csv][Integration]") {
    // Column name with special characters
    TablePrinter<std::string, double> table("Item", {"Value (x,y)", "Note"});
    table.addRow("Test", {42.5, 0.0});
    
    std::ostringstream oss;
    table.exportTo(oss, ExportFormat::CSV);
    
    std::string output = oss.str();
    // Should properly escape column header with comma
    REQUIRE(output.find("\"Value (x,y)\"") != std::string::npos);
}

} // namespace MML::Tests::Tools::ConsolePrinterTests
