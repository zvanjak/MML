# ConsolePrinter - Professional Table Formatting for C++

**Version:** 2.0 (Modernized December 2025)  
**Header:** `mml/tools/ConsolePrinter.h`  
**Namespace:** `MML`

---

## Overview

ConsolePrinter is a modern C++ library for creating beautifully formatted tables in console applications. It provides a fluent builder-pattern API for easy configuration, multiple export formats (Console, CSV, Markdown, LaTeX, HTML), customizable border styles, and professional-quality output.

### Key Features

✅ **Multiple Export Formats** - Console, CSV, TSV, Markdown, LaTeX, HTML  
✅ **Beautiful Borders** - Simple, Rounded, Double, Bold, Markdown styles  
✅ **Fluent Builder API** - Intuitive, chainable configuration  
✅ **Type-Safe Formatting** - Enum classes instead of cryptic characters  
✅ **Stream Abstraction** - Print to any `std::ostream`  
✅ **RAII Stream Safety** - Automatic stream state restoration  
✅ **Auto-Width Calculation** - Automatic column width sizing  
✅ **Backward Compatible** - Legacy API still supported  

---

## Quick Start

### Basic Table (Modern API)

```cpp
#include "mml/tools/ConsolePrinter.h"
using namespace MML;

TablePrinter<double, double> table("x", {"sin(x)", "cos(x)"});
table.addRow(0.0, {0.0, 1.0});
table.addRow(1.57, {1.0, 0.0});
table.addRow(3.14, {0.0, -1.0});
table.print();
```

**Output:**
```
+--------------+--------------+--------------+
|            x |       sin(x) |       cos(x) |
+--------------+--------------+--------------+
|            0 |            0 |            1 |
|         1.57 |            1 |            0 |
|         3.14 |            0 |           -1 |
+--------------+--------------+--------------+
```

### Custom Formatting

```cpp
TablePrinter<double, double> table(
    ColumnFormat("x").width(10).precision(3),
    {
        ColumnFormat("Result").width(15).precision(6).format(FormatType::Scientific),
        ColumnFormat("Error").width(12).precision(3).format(FormatType::Scientific)
    }
);

table.addRow(1.0, {0.123456789, 1.23e-6});
table.addRow(2.0, {0.987654321, 4.56e-5});
table.print();
```

**Output:**
```
+------------+-----------------+--------------+
|          x |          Result |        Error |
+------------+-----------------+--------------+
|          1 |    1.234568e-01 |    1.230e-06 |
|          2 |    9.876543e-01 |    4.560e-05 |
+------------+-----------------+--------------+
```

### Beautiful Borders

```cpp
TablePrinter<double, double> table("x", {"y", "z"});
table.style(TableStyle().border(BorderStyle::Rounded));
table.addRow(1.0, {2.0, 3.0});
table.addRow(4.0, {5.0, 6.0});
table.print();
```

**Output (with Unicode box-drawing):**
```
╭──────────────╮──────────────╮──────────────╮
│            x │            y │            z │
├──────────────┼──────────────┼──────────────┤
│            1 │            2 │            3 │
│            4 │            5 │            6 │
╰──────────────┴──────────────┴──────────────╯
```

---

## Usage Examples

### Example 1: Export to Multiple Formats

```cpp
TablePrinter<double, double> results("Time", {"Position", "Velocity"});
results.addRow(0.0, {0.0, 1.0});
results.addRow(1.0, {0.5, 0.87});
results.addRow(2.0, {0.87, 0.5});

// Console output
results.print();

// CSV export
results.exportToFile("results.csv", ExportFormat::CSV);

// Markdown for documentation
results.exportToFile("results.md", ExportFormat::Markdown);

// LaTeX for papers
results.exportToFile("results.tex", ExportFormat::LaTeX);

// HTML for web
results.exportToFile("results.html", ExportFormat::HTML);
```

**CSV Output:**
```csv
Time,Position,Velocity
0,0,1
1,0.5,0.87
2,0.87,0.5
```

**Markdown Output:**
```markdown
| Time | Position | Velocity |
|------|----------|----------|
| 0 | 0 | 1 |
| 1 | 0.5 | 0.87 |
| 2 | 0.87 | 0.5 |
```

### Example 2: Numerical Analysis Results

```cpp
TablePrinter<double, double> analysis("x", 8, 3,
    {"Exact", "Approx", "Abs Error", "Rel Error"},
    {{12, 7, 'F'}, {12, 7, 'F'}, {13, 6, 'S'}, {10, 3, 'S'}}
);

for (double x = 0.0; x <= 2.0; x += 0.5) {
    double exact = std::sin(x);
    double approx = x - x*x*x/6.0;  // Taylor approximation
    double abs_err = std::abs(exact - approx);
    double rel_err = abs_err / exact;
    
    analysis.addRow(x, {exact, approx, abs_err, rel_err});
}
analysis.Print();
```

### Example 3: Vector Data

```cpp
Vector<Real> time_vec = {0.0, 0.1, 0.2, 0.3, 0.4};
Vector<Real> pos_vec = {0.0, 0.05, 0.2, 0.45, 0.8};
Vector<Real> vel_vec = {0.0, 1.0, 2.0, 3.0, 4.0};

VectorTablePrinter vp;
vp.addVector(ColumnFormat("t").width(8).precision(2), time_vec);
vp.addVector(ColumnFormat("x").width(10).precision(4), pos_vec);
vp.addVector(ColumnFormat("v").width(10).precision(4), vel_vec);
vp.style(TableStyle().border(BorderStyle::Simple));
vp.print();
```

**Output:**
```
+----------+------------+------------+
|        t |          x |          v |
+----------+------------+------------+
|     0.00 |     0.0000 |     0.0000 |
|     0.10 |     0.0500 |     1.0000 |
|     0.20 |     0.2000 |     2.0000 |
|     0.30 |     0.4500 |     3.0000 |
|     0.40 |     0.8000 |     4.0000 |
+----------+------------+------------+
```

### Example 4: Print to String or File Stream

```cpp
TablePrinter<int, std::string> roster("ID", {"Name", "Role"});
roster.addRow(1, {"Alice", "Developer"});
roster.addRow(2, {"Bob", "Tester"});

// Print to std::cout (default)
roster.print();

// Print to std::cerr
roster.print(std::cerr);

// Print to string
std::string table_str = roster.toString();

// Print to file
std::ofstream outfile("roster.txt");
roster.print(outfile);
outfile.close();

// Or use convenience method
roster.printToFile("roster.txt");
```

---

## API Reference

### Enumerations

#### FormatType
Controls numeric formatting style.

```cpp
enum class FormatType {
    General,      // Default formatting (no flags set)
    Fixed,        // Fixed-point notation (std::fixed)
    Scientific,   // Scientific notation (std::scientific)
    Hexadecimal,  // Hexadecimal output (std::hex)
    Integer       // Integer formatting (no decimal point)
};
```

#### Alignment
Controls column text alignment.

```cpp
enum class Alignment {
    Left,    // Left-aligned text
    Right,   // Right-aligned text (default for numbers)
    Center   // Center-aligned text
};
```

#### BorderStyle
Controls table border appearance.

```cpp
enum class BorderStyle {
    None,      // No borders (plain text columns)
    Simple,    // ASCII borders: + - |
    Markdown,  // Markdown table format: | - |
    Rounded,   // Unicode rounded corners: ╭─╮│╰─╯
    Double,    // Unicode double lines: ╔═╗║╚═╝
    Bold       // Unicode bold lines: ┏━┓┃┗━┛
};
```

#### ExportFormat
Controls output format.

```cpp
enum class ExportFormat {
    Console,   // Pretty-printed console output with borders
    CSV,       // Comma-separated values
    TSV,       // Tab-separated values
    Markdown,  // Markdown table format
    LaTeX,     // LaTeX tabular environment
    HTML       // HTML table element
};
```

---

### ColumnFormat Class

Specifies formatting for a single column using builder pattern.

#### Constructors

```cpp
// Simple constructor with column name
explicit ColumnFormat(std::string name);

// Legacy constructor (backward compatibility)
ColumnFormat(std::string name, int width, int precision, char format);
```

#### Builder Methods

```cpp
ColumnFormat& width(int w);              // Set column width
ColumnFormat& precision(int p);          // Set numeric precision
ColumnFormat& format(FormatType fmt);    // Set format type
ColumnFormat& align(Alignment a);        // Set text alignment
ColumnFormat& autoWidth();               // Enable auto-width calculation
```

#### Accessors

```cpp
const std::string& name() const;         // Get column name
int width() const;                       // Get effective width
int precision() const;                   // Get numeric precision
FormatType formatType() const;           // Get format type
Alignment alignment() const;             // Get alignment
```

#### Example

```cpp
auto fmt = ColumnFormat("Temperature")
    .width(12)
    .precision(2)
    .format(FormatType::Fixed)
    .align(Alignment::Right);
```

---

### TableStyle Class

Controls table visual styling.

#### Constructor

```cpp
TableStyle();  // Default: Simple border, show header, no row separators
```

#### Builder Methods

```cpp
TableStyle& border(BorderStyle style);   // Set border style
TableStyle& header(bool show);           // Show/hide header row
TableStyle& rowSeparators(bool show);    // Show/hide separators between rows
TableStyle& compact(bool enable);        // Enable compact mode (future)
```

#### Accessors

```cpp
BorderStyle borderStyle() const;
bool showHeader() const;
bool showRowSeparators() const;
bool compactMode() const;
```

#### Example

```cpp
auto style = TableStyle()
    .border(BorderStyle::Double)
    .header(true)
    .rowSeparators(true);
```

---

### TablePrinter Template Class

Main class for creating and printing tables.

```cpp
template<typename RowTag, typename CellValue>
class TablePrinter;
```

#### Type Parameters

- **RowTag**: Type of row identifiers (e.g., `double` for x-axis values, `int` for IDs)
- **CellValue**: Type of cell values (e.g., `double` for numeric data, `std::string` for text)

#### Constructors

```cpp
// Modern API: simple column names
TablePrinter(std::string tagName, std::vector<std::string> columnNames);

// Modern API: explicit column formats
TablePrinter(ColumnFormat tagFormat, std::vector<ColumnFormat> columnFormats);

// Legacy API (backward compatibility)
TablePrinter(std::string tagName, int tagWidth, int tagPrecision,
             std::vector<std::string> valueNames,
             std::vector<std::tuple<int, int, char>> formatSpecs);
```

#### Configuration Methods

```cpp
TablePrinter& style(const TableStyle& s);
TablePrinter& tagFormat(const ColumnFormat& fmt);
TablePrinter& columnFormat(size_t col, const ColumnFormat& fmt);
```

#### Data Methods

```cpp
void addRow(RowTag tag, std::vector<CellValue> values);
void clear();
void reserve(size_t rows);
```

#### Output Methods

```cpp
// Print to stream
void print(std::ostream& os = std::cout) const;
void Print();  // Legacy method (capital P)

// Export to file
void printToFile(const std::string& filename) const;
void exportToFile(const std::string& filename, ExportFormat format) const;

// Export to stream
void exportTo(std::ostream& os, ExportFormat format) const;

// Convert to string
std::string toString() const;
```

#### Query Methods

```cpp
size_t rowCount() const;
size_t columnCount() const;
```

---

### VectorTablePrinter Class

Specialized printer for Vector<Real> data (common in numerical computing).

#### Methods

```cpp
void addVector(const std::string& name, const Vector<Real>& vec);
void addVector(const ColumnFormat& format, const Vector<Real>& vec);
VectorTablePrinter& style(const TableStyle& s);
void print(std::ostream& os = std::cout) const;
std::string toString() const;
```

#### Example

```cpp
Vector<Real> x_data = {1.0, 2.0, 3.0};
Vector<Real> y_data = {1.414, 1.732, 2.0};

VectorTablePrinter vp;
vp.addVector("x", x_data);
vp.addVector(ColumnFormat("√x").width(10).precision(3), y_data);
vp.print();
```

---

### StreamStateGuard Class (RAII Helper)

Automatically saves and restores stream formatting state.

```cpp
class StreamStateGuard {
public:
    explicit StreamStateGuard(std::ostream& os);
    ~StreamStateGuard();  // Restores state automatically
};
```

Used internally by ConsolePrinter to ensure stream state is never corrupted.

---

## Legacy API (Deprecated)

The legacy API is still supported for backward compatibility but deprecated. Consider migrating to the modern API.

### ColDesc (Deprecated)

```cpp
using ColDesc [[deprecated("Use ColumnFormat instead")]] = ColumnFormat;
```

### VerticalVectorPrinter (Deprecated)

```cpp
class [[deprecated("Use VectorTablePrinter instead")]] VerticalVectorPrinter;
```

**Old usage:**
```cpp
std::vector<ColDesc> names{
    ColDesc("t", 10, 2, 'F'),
    ColDesc("x", 12, 5, 'F')
};
std::vector<Vector<Real>*> vecs{&time_vec, &pos_vec};
VerticalVectorPrinter vvp(names, vecs);
vvp.Print();
```

**Modern equivalent:**
```cpp
VectorTablePrinter vp;
vp.addVector(ColumnFormat("t").width(10).precision(2), time_vec);
vp.addVector(ColumnFormat("x").width(12).precision(5), pos_vec);
vp.print();
```

---

## Advanced Topics

### Custom Output Streams

Print to any C++ output stream:

```cpp
TablePrinter<int, double> table("ID", {"Value"});
table.addRow(1, {3.14});

// Standard streams
table.print(std::cout);
table.print(std::cerr);
table.print(std::clog);

// File stream
std::ofstream file("output.txt");
table.print(file);

// String stream
std::ostringstream buffer;
table.print(buffer);
std::string result = buffer.str();
```

### Export Format Details

#### CSV Format
- Standard comma-separated values
- Header row with column names
- No borders or formatting
- Compatible with Excel, Python pandas, R

#### TSV Format
- Tab-separated values
- Header row with column names
- Useful when data contains commas

#### Markdown Format
- GitHub-flavored markdown tables
- Perfect for README files
- Header with separator row
- Pipe-delimited

#### LaTeX Format
- Standard `tabular` environment
- All columns centered with borders
- `\hline` separators
- Ready for `\begin{table}` wrapper

#### HTML Format
- Semantic `<table>` with `<thead>` and `<tbody>`
- Border attribute set to "1"
- Ready for web embedding

### Auto-Width Calculation

Set column width to `AUTO_WIDTH` to automatically calculate based on content:

```cpp
TablePrinter<std::string, std::string> table(
    ColumnFormat("Category").autoWidth(),
    {
        ColumnFormat("Description").autoWidth(),
        ColumnFormat("Status").width(10)  // Fixed width
    }
);

table.addRow("A", {"Short", "OK"});
table.addRow("Very Long Category Name", {"This is a much longer description", "PASS"});
table.print();
// Columns 0 and 1 will auto-size to fit content
```

### Error Handling

```cpp
try {
    TablePrinter<double, double> table("x", {"y"});
    table.addRow(1.0, {2.0, 3.0});  // Wrong number of values
} catch (const std::invalid_argument& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    // "Number of values does not match number of columns"
}
```

---

## Performance Tips

1. **Reserve rows** if you know the size: `table.reserve(1000);`
2. **Use fixed width** instead of auto-width for large tables
3. **Minimize precision** for faster formatting
4. **Use CSV/TSV** for very large datasets (no border overhead)
5. **Reuse table objects** - call `clear()` instead of creating new instances

---

## Migration Guide

### From Old API to New API

**Before:**
```cpp
TablePrinter<double, double> table("x", 8, 3,
    {"y1", "y2"},
    {{10, 5, 'F'}, {12, 6, 'S'}}
);
```

**After:**
```cpp
TablePrinter<double, double> table(
    ColumnFormat("x").width(8).precision(3),
    {
        ColumnFormat("y1").width(10).precision(5).format(FormatType::Fixed),
        ColumnFormat("y2").width(12).precision(6).format(FormatType::Scientific)
    }
);
```

**Benefits:**
- ✅ More readable and self-documenting
- ✅ Type-safe enums instead of chars
- ✅ Easier to modify individual columns
- ✅ Chainable builder pattern

---

## Examples from Real Code

### Example from readme7_derivation.cpp

```cpp
TablePrinter<double, double> print_data("x", 8, 3,
    {"Exact der.", "Nder1", "Nder1 err.", "Nder2", "Nder2 err.", 
     "Nder4", "Nder4 err.", "Nder6", "Nder6 err.", "Nder8", "Nder8 err."},
    {{12,7,'F'}, {13,7,'F'}, {15,6,'S'}, {13,7,'F'}, {15,6,'S'}, 
     {13,7,'F'}, {15,6,'S'}, {13,7,'F'}, {15,6,'S'}, {13,7,'F'}, {15,6,'S'}}
);

for (int i = 0; i < 20; i++) {
    double x = -7.0 + 14.0 * i / 19.0;
    double exact_der = f_der(x);
    double num_der1 = MML::Derivation::NDer1(f, x);
    // ... more calculations ...
    print_data.addRow(x, {exact_der, num_der1, err1, num_der2, err2, 
                          num_der4, err4, num_der6, err6, num_der8, err8});
}
print_data.Print();
```

### Example from example8_gravity.cpp

```cpp
std::vector<ColDesc> vecNames{
    ColDesc("t", 11, 2, 'F'),
    ColDesc("x1", 10, 5, 'F'), ColDesc("y1", 10, 5, 'F'),
    ColDesc("x2", 10, 5, 'F'), ColDesc("y2", 10, 5, 'F'),
    ColDesc("v1_x", 10, 5, 'F'), ColDesc("v1_y", 10, 5, 'F'),
    ColDesc("v2_x", 10, 5, 'F'), ColDesc("v2_y", 10, 5, 'F')
};

std::vector<Vector<Real>*> vecVals{
    &t_vals, &x1_vals, &y1_vals, &x2_vals, &y2_vals,
    &v1_x_vals, &v1_y_vals, &v2_x_vals, &v2_y_vals
};

VerticalVectorPrinter vvp(vecNames, vecVals);
vvp.Print();
```

---

## Troubleshooting

### Issue: Borders look wrong on Windows

**Cause:** Windows console may not support Unicode box-drawing characters.

**Solution:** Use `BorderStyle::Simple` (ASCII) instead of Unicode styles:
```cpp
table.style(TableStyle().border(BorderStyle::Simple));
```

### Issue: Precision not applied correctly

**Cause:** Forgot to set format type to Fixed or Scientific.

**Solution:** Explicitly set format type:
```cpp
ColumnFormat("value").precision(3).format(FormatType::Fixed)
```

### Issue: Columns too narrow/too wide

**Solution:** Either set explicit width or use auto-width:
```cpp
// Explicit width
ColumnFormat("data").width(15)

// Auto-width
ColumnFormat("data").autoWidth()
```

### Issue: Stream formatting changed after printing

**Cause:** Legacy code or custom stream manipulation.

**Solution:** The StreamStateGuard automatically handles this in modern API. For manual control:
```cpp
{
    StreamStateGuard guard(std::cout);
    // Your printing code
}  // State restored automatically
```

---

## See Also

- **Serializer** - Save/load data structures to files
- **Visualizers** - Graphical plotting and visualization
- **Vector Class** - Mathematical vector container

---

## Version History

**v2.0 (December 2025)**
- Complete modernization with builder pattern
- Multiple export formats (CSV, TSV, Markdown, LaTeX, HTML)
- Border styles (Simple, Rounded, Double, Bold, Markdown)
- RAII stream state management
- Auto-width calculation
- Backward-compatible legacy API

**v1.0 (Original)**
- Basic table printing
- Simple column formatting
- VerticalVectorPrinter for vectors
