///////////////////////////////////////////////////////////////////////////////////////////
/// MML Documentation Demo: ConsolePrinter - Formatted Console Output
///////////////////////////////////////////////////////////////////////////////////////////

#include "MMLBase.h"
#include "base/Vector/Vector.h"
#include "tools/ConsolePrinter.h"

#include <iostream>
#include <iomanip>
#include <cmath>

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
