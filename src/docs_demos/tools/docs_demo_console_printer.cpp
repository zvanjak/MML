///////////////////////////////////////////////////////////////////////////////////////////
/// MML Documentation Demo: ConsolePrinter - Formatted Console Output
///////////////////////////////////////////////////////////////////////////////////////////

#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#include <mml/tools/ConsolePrinter.h>

#include <iostream>
#include <iomanip>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>

using namespace MML;

namespace MML::docs_demos::console_printer
{
    ///////////////////////////////////////////////////////////////////////////////////////
    ///                              BORDER STYLES                                       ///
    ///////////////////////////////////////////////////////////////////////////////////////
    // BorderStyle::None     - No borders
    // BorderStyle::Simple   - ASCII: + - |
    // BorderStyle::Markdown - | - |
    // BorderStyle::Rounded  - Unicode: ╭─╮│╰─╯
    // BorderStyle::Double   - Unicode: ╔═╗║╚═╝
    // BorderStyle::Bold     - Unicode: ┏━┓┃┗━┛

    void demo_table_printer()
    {
        std::cout << "\n=== TablePrinter Demo ===\n";

        // TablePrinter<RowTag, CellValue> - templated on tag and cell types
        // Constructor takes tag column name and vector of column names
        TablePrinter<std::string, double> table("Method", {"x", "f(x)", "Error"});

        // Add rows: tag first, then values
        table.addRow("Newton", {1.0, 0.841471, 1.2e-10});
        table.addRow("Bisection", {2.0, 0.909297, 3.5e-8});
        table.addRow("Secant", {3.0, 0.141120, 8.7e-12});

        // Print with different styles
        std::cout << "\nSimple border style:\n";
        table.style(TableStyle().border(BorderStyle::Simple)).print();

        std::cout << "\nDouble border style:\n";
        table.style(TableStyle().border(BorderStyle::Double)).print();

        std::cout << "\nRounded border style:\n";
        table.style(TableStyle().border(BorderStyle::Rounded)).print();
    }

    void demo_column_format()
    {
        std::cout << "\n=== ColumnFormat Demo ===\n";

        // ColumnFormat uses builder pattern
        ColumnFormat col1("Value");
        col1.width(12).precision(6).format(FormatType::Fixed).align(Alignment::Right);

        ColumnFormat col2("Scientific");
        col2.width(15).precision(4).format(FormatType::Scientific).align(Alignment::Right);

        // Format values
        double val = 12345.6789;
        std::cout << "Value: " << val << "\n";
        std::cout << "Fixed format:      " << col1.formatValue(val) << "\n";
        std::cout << "Scientific format: " << col2.formatValue(val) << "\n";
    }

    void demo_table_style()
    {
        std::cout << "\n=== TableStyle Demo ===\n";

        // TableStyle also uses builder pattern
        TableStyle style;
        style.border(BorderStyle::Bold)
             .header(true)
             .rowSeparators(true)
             .compact(false);

        // Create table with this style - constructor takes tag name and column names
        TablePrinter<int, double> table("n", {"n^2", "sqrt(n)"});

        for (int n = 1; n <= 5; ++n) {
            table.addRow(n, {static_cast<double>(n * n), std::sqrt(static_cast<double>(n))});
        }

        table.style(style).print();
    }

    void demo_vector_table_printer()
    {
        std::cout << "\n=== VectorTablePrinter Demo ===\n";

        // Create sample vectors
        Vector<Real> x({0.0, 0.5, 1.0, 1.5, 2.0});
        Vector<Real> sinx({0.0, 0.479, 0.841, 0.997, 0.909});
        Vector<Real> cosx({1.0, 0.878, 0.540, 0.071, -0.416});

        // Print vectors side by side
        VectorTablePrinter vtp;
        vtp.addVector(ColumnFormat("x").width(10).precision(2), x);
        vtp.addVector(ColumnFormat("sin(x)").width(12).precision(4), sinx);
        vtp.addVector(ColumnFormat("cos(x)").width(12).precision(4), cosx);

        vtp.style(TableStyle().border(BorderStyle::Rounded)).print();
    }

    void demo_utf8_helpers()
    {
        std::cout << "\n=== UTF-8 Helpers ===\n";

        // String repetition (handles multi-byte UTF-8)
        std::string boxChar = "═";
        std::cout << "Box line: " << Utf8::repeatString(boxChar, 20) << "\n";

        // Check single display character
        std::cout << "Is '═' single char: " << (Utf8::isSingleDisplayChar("═") ? "yes" : "no") << "\n";
        std::cout << "Is 'abc' single char: " << (Utf8::isSingleDisplayChar("abc") ? "yes" : "no") << "\n";
    }

    void demo_stream_state_guard()
    {
        std::cout << "\n=== StreamStateGuard Demo ===\n";

        // StreamStateGuard preserves stream formatting across operations
        std::cout << "Before guard - default precision\n";
        std::cout << "Pi = " << 3.14159265358979 << "\n";

        {
            StreamStateGuard guard(std::cout);
            std::cout << std::fixed << std::setprecision(2);
            std::cout << "Inside guard - fixed with 2 decimals\n";
            std::cout << "Pi = " << 3.14159265358979 << "\n";
        } // guard restores state here

        std::cout << "After guard - restored default\n";
        std::cout << "Pi = " << 3.14159265358979 << "\n";
    }

    void Run()
    {
        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "       MML Documentation Demo: ConsolePrinter\n";
        std::cout << std::string(70, '=') << "\n";

        demo_column_format();
        demo_table_printer();
        demo_table_style();
        demo_vector_table_printer();
        demo_utf8_helpers();
        demo_stream_state_guard();

        std::cout << "\n" << std::string(70, '=') << "\n";
    }
} // namespace

void Docs_Demo_ConsolePrinter() {
    MML::docs_demos::console_printer::Run();
}

void Docs_Demo_TableExportContracts()
{
	const auto directory = std::filesystem::temp_directory_path() /
		("mml_table_export_demo_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
	std::filesystem::create_directory(directory);
	struct Cleanup {
		std::filesystem::path directory;
		~Cleanup() {
			std::error_code error;
			std::filesystem::remove_all(directory, error);
			if (error)
				std::cerr << "Cannot clean table export demo: " << error.message() << '\n';
		}
	} cleanup{directory};
	const std::string filename = (directory / "table.csv").string();
	TablePrinter<std::string, Real> table("Tag", {"Value"});
	table.addRow("row", {Real(1.25)});
	table.exportToFile(filename, ExportFormat::CSV);
	{
		std::ifstream input(filename);
		if (!input)
			throw FileIOError("Cannot read table export demo", filename);
		const std::string content(std::istreambuf_iterator<char>(input), {});
		if (content != "Tag,Value\nrow,1.25\n")
			throw FileIOError("Unexpected table export content", filename);
	}
	const auto badFilename = (directory / "missing" / "table.csv").string();
	bool rejected = false;
	try {
		table.printToFile(badFilename);
	} catch (const FileIOError& error) {
		rejected = true;
		if (error.filename() != badFilename)
			throw ArgumentError("Table export error lost the destination filename");
		std::cout << "Open failure reported: " << error.what() << '\n';
	}
	if (!rejected)
		throw ArgumentError("Table export accepted an unavailable destination");
	std::cout << "Table CSV written, flushed, closed, and read back successfully\n";
}

void Docs_Demo_ColumnFormatContracts()
{
	const auto integer = ColumnFormat("Rounded").width(5).format(FormatType::Integer);
	if (integer.formatValue(Real(1.5)) != "2" || integer.formatValue(Real(-2.5)) != "-3")
		throw ArgumentError("Integer display rounding contract failed");
	std::cout << integer.formatAligned(integer.formatValue(Real(-2.5))) << '\n';
	bool invalidWidth = false;
	try {
		ColumnFormat("Invalid").width(-2);
	} catch (const ArgumentError& error) {
		invalidWidth = true;
		std::cout << error.what() << '\n';
	}
	if (!invalidWidth)
		throw ArgumentError("Negative column width was accepted");
	bool invalidType = false;
	try {
		ColumnFormat("Hex").format(FormatType::Hexadecimal).formatValue(Real(1.25));
	} catch (const ArgumentError& error) {
		invalidType = true;
		std::cout << error.what() << '\n';
	}
	if (!invalidType)
		throw ArgumentError("Floating hexadecimal column was accepted");
}

void Docs_Demo_BorrowedVectorTables()
{
	Vector<Real> values{1};
	VectorTablePrinter table;
	table.addVector(ColumnFormat("live").autoWidth().precision(2).format(FormatType::Fixed), values);
	const auto initial = table.toString();
	values.push_back(Real(1234.5));
	values.push_back(Real(-2));
	const auto grown = table.toString();
	if (grown.find("1234.50") == std::string::npos || grown.find("-2.00") == std::string::npos || grown == initial)
		throw std::runtime_error("Borrowed vector table did not observe its live size and automatic width");
	std::cout << "Live vector after growth:\n" << grown;
	values = Vector<Real>{2};
	const auto shrunk = table.toString();
	if (shrunk != "+------+\n| live |\n+------+\n| 2.00 |\n+------+\n")
		throw std::runtime_error("Borrowed vector table retained stale rows or automatic width after shrinking");
	std::cout << "Live vector after shrinking:\n" << shrunk;
	const Vector<Real> constant{3};
	VectorTablePrinter constBorrow;
	constBorrow.addVector("const lvalue", constant);
	if (constBorrow.toString().find("const lvalue") == std::string::npos)
		throw std::runtime_error("Const lvalue borrowing failed");
	bool rejectedNull = false;
	try {
		VectorTablePrinter invalid({}, {nullptr});
	} catch (const std::invalid_argument& error) {
		rejectedNull = true;
		std::cout << "Null input rejected: " << error.what() << '\n';
	}
	if (!rejectedNull) throw std::runtime_error("Null vector pointer was accepted");
	std::ostringstream out;
	out << std::scientific;
	out.precision(3);
	const auto flags = out.flags();
	table.print(out);
	if (out.flags() != flags || out.precision() != 3)
		throw std::runtime_error("Vector table modified caller stream formatting");
}
