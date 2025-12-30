# ConsolePrinter - Professional Table Formatting for C++

**Version:** 2.0 (Modernized December 2025)  
**Header:** `mml/tools/ConsolePrinter.h`  
**Namespace:** `MML`

---

## Overview

ConsolePrinter is a modern C++ library for creating beautifully formatted tables in console applications. It provides a fluent builder-pattern API for easy configuration, multiple export formats (Console, CSV, Markdown, LaTeX, HTML), customizable border styles, and professional-quality output.

### Key Features

✅ **Multiple Export Formats** - Console, CSV  
✅ **Beautiful Borders** - Simple, Rounded, Double, Bold, Markdown styles  
✅ **Fluent Builder API** - Intuitive, chainable configuration  
✅ **Type-Safe Formatting** - Enum classes instead of cryptic characters  
✅ **Stream Abstraction** - Print to any `std::ostream`  
✅ **RAII Stream Safety** - Automatic stream state restoration  
✅ **Auto-Width Calculation** - Automatic column width sizing  
✅ **Backward Compatible** - Legacy API still supported  

## Highlights
- Professional, multi-format table output for scientific results.
- Fluent builder API with type-safe enums and RAII stream safety.
- Auto width computation, border styles, and export to Markdown/LaTeX/HTML.
- Interoperates with Tools: Serializer and Visualizers for reporting.
- Backward-compatible legacy API maintained alongside modern interfaces.

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

**Double borders:**
```
╔════════╗
║ Header ║
╠════════╣
║ Data   ║
╚════════╝
```

**Bold borders:**
```
┏━━━━━━━━┓
┃ Header ┃
┡━━━━━━━━┩
│ Data   │
└────────┘
```

---

## Core Components

### 1. ColumnFormat

Builder-pattern class for column configuration.

```cpp
ColumnFormat(std::string name);  // Basic constructor

// Builder methods
ColumnFormat& width(int w);
ColumnFormat& precision(int p);
ColumnFormat& format(FormatType fmt);
ColumnFormat& align(Alignment a);
ColumnFormat& autoWidth();
```

**Example:**
```cpp
auto fmt = ColumnFormat("Temperature")
    .width(12)
    .precision(2)
    .format(FormatType::Fixed)
    .align(Alignment::Right);
```

### 2. TableStyle

Configuration for table appearance.

```cpp
TableStyle();  // Default: Simple border, show header, no row separators

// Builder methods
TableStyle& border(BorderStyle style);
TableStyle& header(bool show);
TableStyle& rowSeparators(bool show);
TableStyle& compact(bool enable);
```

**Example:**
```cpp
auto style = TableStyle()
    .border(BorderStyle::Double)
    .header(true)
    .rowSeparators(true);
```

### 3. TablePrinter

Main table printing class with templated row tags and cell values.

```cpp
template<typename RowTag, typename CellValue>
class TablePrinter {
    // Modern API: simple column names
    TablePrinter(std::string tagName, std::vector<std::string> columnNames);
    
    // Modern API: explicit column formats
    TablePrinter(ColumnFormat tagFormat, std::vector<ColumnFormat> columnFormats);
    
    void addRow(RowTag tag, std::vector<CellValue> values);
    void print(std::ostream& os = std::cout) const;
    void exportTo(std::ostream& os, ExportFormat format) const;
    void exportToFile(const std::string& filename, ExportFormat format) const;
};
```

### 4. VectorTablePrinter

Specialized printer for displaying multiple vectors side-by-side.

```cpp
class VectorTablePrinter {
    void addVector(const std::string& name, const Vector<Real>& vec);
    void addVector(const ColumnFormat& format, const Vector<Real>& vec);
    VectorTablePrinter& style(const TableStyle& s);
    void print(std::ostream& os = std::cout) const;
    std::string toString() const;
};
```

### 5. StreamStateGuard

RAII class preserving stream formatting state.

```cpp
class StreamStateGuard {
public:
    explicit StreamStateGuard(std::ostream& os);
    ~StreamStateGuard();  // Restores state automatically
};
```

Used internally by ConsolePrinter to ensure stream state is never corrupted.

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

---

### TableStyle Class

Controls table visual styling.

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

---

### TablePrinter Template Class

#### Type Parameters

- **RowTag**: Type of row identifiers (e.g., `double` for x-axis values, `int` for IDs)
- **CellValue**: Type of cell values (e.g., `double` for numeric data, `std::string` for text)

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

**LaTeX Output:**
```latex
\begin{tabular}{|c|c|c|}
\hline
Time & Position & Velocity \\
\hline
0 & 0 & 1 \\
1 & 0.5 & 0.87 \\
2 & 0.87 & 0.5 \\
\hline
\end{tabular}
```

**HTML Output:**
```html
<table border="1">
  <thead>
    <tr><th>Time</th><th>Position</th><th>Velocity</th></tr>
  </thead>
  <tbody>
    <tr><td>0</td><td>0</td><td>1</td></tr>
    <tr><td>1</td><td>0.5</td><td>0.87</td></tr>
  </tbody>
</table>
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

### Example 5: Convergence Analysis Table

```cpp
TablePrinter<int, double> convTable("N", {"h", "Error", "Ratio"});

ColumnFormat nCol("N").width(8).align(Alignment::Right);
ColumnFormat hCol("h").width(12).precision(6).format(FormatType::Scientific);
ColumnFormat errCol("Error").width(12).precision(3).format(FormatType::Scientific);
ColumnFormat ratioCol("Ratio").width(8).precision(2);

TablePrinter<int, double> table(nCol, {hCol, errCol, ratioCol});

double prevError = 0;
for (int n = 10; n <= 10000; n *= 10) {
    double h = 1.0 / n;
    double error = computeError(n);
    double ratio = (prevError > 0) ? prevError / error : 0;
    
    table.addRow(n, {h, error, ratio});
    prevError = error;
}

table.style(TableStyle().border(BorderStyle::Bold));
table.print();
```

**Output:**
```
┏━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━━━━━━━━━━┳━━━━━━━━┓
┃      N ┃ h            ┃ Error        ┃  Ratio ┃
┡━━━━━━━━╇━━━━━━━━━━━━━━╇━━━━━━━━━━━━━━╇━━━━━━━━┩
│     10 │ 1.000000e-01 │ 3.142e-03    │   0.00 │
│    100 │ 1.000000e-02 │ 3.142e-05    │ 100.00 │
│   1000 │ 1.000000e-03 │ 3.142e-07    │ 100.00 │
│  10000 │ 1.000000e-04 │ 3.142e-09    │ 100.00 │
└────────┴──────────────┴──────────────┴────────┘
```

### Example 6: String Tables

```cpp
TablePrinter<std::string, std::string> fileTable("File", {"Size", "Modified"});

fileTable.addRow("data.txt", {"1.2 MB", "2025-01-15"});
fileTable.addRow("results.csv", {"3.5 MB", "2025-01-16"});
fileTable.addRow("plot.png", {"512 KB", "2025-01-14"});

TableStyle markdown;
markdown.border(BorderStyle::Markdown);
fileTable.style(markdown).print();
```

**Output:**
```
| File        | Size   | Modified   |
|-------------|--------|------------|
| data.txt    | 1.2 MB | 2025-01-15 |
| results.csv | 3.5 MB | 2025-01-16 |
| plot.png    | 512 KB | 2025-01-14 |
```

### Example 7: Builder Pattern Chaining

```cpp
TablePrinter<double, double> table("x", {"f(x)", "g(x)"});

table.style(TableStyle()
            .border(BorderStyle::Rounded)
            .header(true)
            .rowSeparators(false))
     .columnFormat(0, ColumnFormat("f(x)")
                       .width(12)
                       .precision(4)
                       .format(FormatType::Scientific))
     .columnFormat(1, ColumnFormat("g(x)")
                       .width(10)
                       .precision(2));

for (double x : xValues) {
    table.addRow(x, {f(x), g(x)});
}

table.print();
```

### Example 8: Performance Results Table

```cpp
#include "mml/tools/Timer.h"
#include "mml/tools/ConsolePrinter.h"

Timer timer;
std::vector<std::string> algorithms = {"Euler", "RK4", "RKCK"};
std::vector<double> times;

for (auto& algo : algorithms) {
    timer.Start();
    runAlgorithm(algo);
    timer.MarkTime(algo);
    times.push_back(timer.GetIntervalTime(times.size()));
}

TablePrinter<std::string, double> table("Algorithm", {"Time (s)", "Speedup"});

double baseline = times[0];
for (size_t i = 0; i < algorithms.size(); ++i) {
    table.addRow(algorithms[i], {times[i], baseline / times[i]});
}

TableStyle style;
style.border(BorderStyle::Bold).header(true);
table.style(style).print();
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

**Output:**
```
+--------------------------+-----------------------------------+------------+
| Category                 | Description                       | Status     |
+--------------------------+-----------------------------------+------------+
| A                        | Short                             | OK         |
| Very Long Category Name  | This is a much longer description | PASS       |
+--------------------------+-----------------------------------+------------+
```

### Stream State Preservation

```cpp
void printFormattedData() {
    StreamStateGuard guard(std::cout);  // Save state
    
    std::cout << std::scientific << std::setprecision(3);
    std::cout << someValue << "\n";
    
    // Original formatting restored automatically on guard destruction
}

std::cout << nextValue;  // Uses original formatting
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

try {
    table.exportToFile("output.csv", ExportFormat::CSV);
} catch (const std::runtime_error& e) {
    std::cerr << "Export failed: " << e.what() << '\n';
}
```

---

## Best Practices

### 1. Use Appropriate Precision

Choose precision based on your data and readability needs.

```cpp
// Too much precision - cluttered
ColumnFormat col("Value").precision(15);  // 3.14159265358979

// Good - readable
ColumnFormat col("Value").precision(3);   // 3.142
```

### 2. Choose Format Type Based on Data

```cpp
// Small numbers - Fixed
ColumnFormat("Temperature").format(FormatType::Fixed).precision(1);

// Very large/small - Scientific
ColumnFormat("Population").format(FormatType::Scientific).precision(2);

// Whole numbers - Integer
ColumnFormat("Count").format(FormatType::Integer);
```

### 3. Use Auto-Width for Variable Data

```cpp
// Fixed width
ColumnFormat("Name").width(20);  // May truncate or waste space

// Auto-width
ColumnFormat("Name").autoWidth();  // Adapts to content
```

### 4. Match Border Style to Output Medium

```cpp
// Console - Simple or Rounded
TableStyle().border(BorderStyle::Rounded);

// Markdown files - Markdown
TableStyle().border(BorderStyle::Markdown);

// Papers/reports - None (export to LaTeX)
TableStyle().border(BorderStyle::None);
```

### 5. Reserve Capacity for Large Tables

Pre-allocate memory for better performance with large datasets.

```cpp
TablePrinter<int, double> table("Index", {"Value"});
table.reserve(10000);  // Pre-allocate for performance

for (int i = 0; i < 10000; ++i) {
    table.addRow(i, {compute(i)});
}
```

### 6. Error Handling

Always handle potential errors when adding rows or exporting data.

```cpp
try {
    table.addRow(x, {value1, value2});
} catch (const std::invalid_argument& e) {
    std::cerr << "Column count mismatch: " << e.what() << '\n';
}
```

---

## Performance Considerations

- **Auto-width calculation**: O(rows × columns) - done once before printing
- **Memory**: ~40 bytes per cell + data size
- **Large tables**: 
  - Use `reserve()` to avoid reallocations
  - Example: 10,000 rows × 5 columns ≈ 2 MB
  - Export to CSV (~50 KB) for large datasets instead of console output (~500 KB with borders)
- **Export overhead**: Minimal, I/O-bound
- **Stream state**: RAII guard has negligible overhead
- **File size estimates**:
  - Console format: ~100 bytes per row (with borders)
  - CSV/TSV: ~50 bytes per row
  - Markdown: ~80 bytes per row
  - LaTeX: ~120 bytes per row
  - HTML: ~150 bytes per row

**Performance Tips:**

1. **Reserve rows** if you know the size: `table.reserve(1000);`
2. **Use fixed width** instead of auto-width for large tables
3. **Minimize precision** for faster formatting
4. **Use CSV/TSV** for very large datasets (no border overhead)
5. **Reuse table objects** - call `clear()` instead of creating new instances

---

## Common Issues

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

### Issue: Column count mismatch

**Cause:** Number of values in `addRow()` doesn't match column count.

**Solution:** Ensure the correct number of values:
```cpp
TablePrinter<int, double> table("ID", {"Value1", "Value2"});  // 2 columns
table.addRow(1, {3.14, 2.71});  // ✓ Correct: 2 values
table.addRow(2, {1.41});        // ✗ Error: only 1 value
```

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

### Format Character Migration

| Old Char | New Enum |
|----------|----------|
| `'F'` | `FormatType::Fixed` |
| `'S'` | `FormatType::Scientific` |
| `'G'` | `FormatType::General` |
| `'X'` | `FormatType::Hexadecimal` |
| `'I'` | `FormatType::Integer` |

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

## Summary

**ConsolePrinter** provides:
- ✅ Modern builder-pattern API
- ✅ 6 border styles (Simple, Markdown, Rounded, Double, Bold, None)
- ✅ 6 export formats (Console, CSV, TSV, Markdown, LaTeX, HTML)
- ✅ Auto-width column sizing
- ✅ Flexible formatting (Fixed, Scientific, Integer, Hex, General)
- ✅ Left/Right/Center alignment
- ✅ RAII stream state management
- ✅ Specialized vector table printer
- ✅ Backward compatible legacy API

**Key takeaway:** Professional-quality formatted output for console, files, and documents with minimal code.

---

## See Also

- **[Serializer](Serializer.md)** - Save/load data structures to files
- **[Visualizers](Visualizers.md)** - Graphical plotting and visualization
- **[Vectors](../base/Vectors.md)** - Mathematical vector containers
- **[Timer & ThreadPool](Timer_ThreadPool.md)** - Performance measurement for table generation

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
