# ConsolePrinter - Formatted Console Output

## Overview

Modern C++ API for formatted table output to console and file with multiple export formats. Provides builder-pattern configuration, multiple border styles, alignment options, and export to CSV/TSV/Markdown/LaTeX/HTML.

**Key Features:**
- Builder pattern for column and table configuration
- Multiple border styles (Simple, Markdown, Rounded, Double, Bold, None)
- Export formats: Console, CSV, TSV, Markdown, LaTeX, HTML
- Auto-width calculation
- RAII stream state management
- Specialized vector table printer

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

**Enums:**
```cpp
enum class FormatType {
    General,      // Default formatting
    Fixed,        // std::fixed
    Scientific,   // std::scientific
    Hexadecimal,  // std::hex
    Integer       // No decimal point
};

enum class Alignment {
    Left, Right, Center
};
```

**Example:**
```cpp
ColumnFormat col("Temperature")
    .width(12)
    .precision(3)
    .format(FormatType::Fixed)
    .align(Alignment::Right);
```

### 2. TableStyle

Configuration for table appearance.

```cpp
enum class BorderStyle {
    None, Simple, Markdown, Rounded, Double, Bold
};

TableStyle style;
style.border(BorderStyle::Rounded)
     .header(true)
     .rowSeparators(false)
     .compact(false);
```

### 3. TablePrinter

Main table printing class with templated row tags and cell values.

```cpp
template<typename RowTag, typename CellValue>
class TablePrinter {
    TablePrinter(std::string tagName, std::vector<std::string> columnNames);
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
    void print(std::ostream& os = std::cout) const;
};
```

### 5. StreamStateGuard

RAII class preserving stream formatting state.

```cpp
StreamStateGuard guard(std::cout);  // Saves state
// Modify stream formatting
// Destructor restores original state
```

---

## Usage Examples

### Basic Table

```cpp
#include "tools/ConsolePrinter.h"

TablePrinter<double, double> table("x", {"sin(x)", "cos(x)", "tan(x)"});

for (double x = 0; x <= M_PI; x += M_PI/4) {
    table.addRow(x, {std::sin(x), std::cos(x), std::tan(x)});
}

table.print();
```

**Output:**
```
+--------+--------+--------+--------+
| x      | sin(x) | cos(x) | tan(x) |
+--------+--------+--------+--------+
| 0      | 0      | 1      | 0      |
| 0.7854 | 0.7071 | 0.7071 | 1      |
| 1.5708 | 1      | 0      | inf    |
| 2.3562 | 0.7071 | -0.7071| -1     |
| 3.1416 | 0      | -1     | -0     |
+--------+--------+--------+--------+
```

### Custom Column Formatting

```cpp
ColumnFormat xCol("x")
    .width(10)
    .precision(4)
    .format(FormatType::Fixed)
    .align(Alignment::Right);

ColumnFormat sinCol("sin(x)")
    .width(12)
    .precision(6)
    .format(FormatType::Scientific);

ColumnFormat cosCol("cos(x)")
    .width(12)
    .precision(6)
    .format(FormatType::Fixed);

TablePrinter<double, double> table(xCol, {sinCol, cosCol});
table.addRow(0.5, {std::sin(0.5), std::cos(0.5)});
table.print();
```

### Border Styles

```cpp
TableStyle rounded;
rounded.border(BorderStyle::Rounded)
       .header(true)
       .rowSeparators(true);

table.style(rounded).print();
```

**Rounded borders:**
```
╭────────╮
│ Header │
├────────┤
│ Data   │
├────────┤
│ More   │
╰────────╯
```

**Double borders:**
```
╔════════╗
║ Header ║
╠════════╣
║ Data   ║
╚════════╝
```

**Markdown style:**
```
| Column1 | Column2 |
|---------|---------|
| Data1   | Data2   |
```

### Export Formats

```cpp
// Export to CSV
table.exportToFile("data.csv", ExportFormat::CSV);

// Export to Markdown
table.exportToFile("table.md", ExportFormat::Markdown);

// Export to LaTeX
table.exportToFile("table.tex", ExportFormat::LaTeX);

// Export to HTML
table.exportToFile("table.html", ExportFormat::HTML);

// Export to TSV
table.exportToFile("data.tsv", ExportFormat::TSV);
```

**CSV output:**
```csv
x,sin(x),cos(x)
0,0,1
0.7854,0.7071,0.7071
1.5708,1,0
```

**Markdown output:**
```markdown
| x | sin(x) | cos(x) |
|-----------|---------|---------|
| 0 | 0 | 1 |
| 0.7854 | 0.7071 | 0.7071 |
```

**LaTeX output:**
```latex
\begin{tabular}{|c|c|c|}
\hline
x & sin(x) & cos(x) \\
\hline
0 & 0 & 1 \\
0.7854 & 0.7071 & 0.7071 \\
\hline
\end{tabular}
```

### Auto-Width Columns

```cpp
ColumnFormat col1("Name").autoWidth();
ColumnFormat col2("Value").autoWidth().precision(3);

TablePrinter<std::string, double> table(col1, {col2});
table.addRow("pi", M_PI);
table.addRow("e", M_E);
table.addRow("sqrt(2)", std::sqrt(2));
table.print();
```

**Output automatically sized:**
```
+---------+-----------+
| Name    | Value     |
+---------+-----------+
| pi      | 3.142     |
| e       | 2.718     |
| sqrt(2) | 1.414     |
+---------+-----------+
```

### Vector Tables

```cpp
VectorTablePrinter vecTable;

Vector<Real> x = {0, 0.5, 1.0, 1.5, 2.0};
Vector<Real> y = {1, 2, 4, 8, 16};
Vector<Real> z = {0, 0.25, 1.0, 2.25, 4.0};

vecTable.addVector("x", x);
vecTable.addVector("y", y);
vecTable.addVector("z", z);

TableStyle style;
style.border(BorderStyle::Double).header(true);
vecTable.style(style).print();
```

**Output:**
```
╔════╦══════╦══════╗
║ x  ║ y    ║ z    ║
╠════╬══════╬══════╣
║ 0  ║ 1    ║ 0    ║
║0.5 ║ 2    ║ 0.25 ║
║1.0 ║ 4    ║ 1.0  ║
║1.5 ║ 8    ║ 2.25 ║
║2.0 ║ 16   ║ 4.0  ║
╚════╩══════╩══════╝
```

### Performance Results Table

```cpp
#include "tools/Timer.h"
#include "tools/ConsolePrinter.h"

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

### Convergence Analysis Table

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

### String Tables

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

---

## Advanced Features

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

### Mixed Data Types

```cpp
TablePrinter<int, std::variant<double, std::string>> mixedTable(
    "ID", {"Value", "Status"}
);

mixedTable.addRow(1, {3.14159, "OK"});
mixedTable.addRow(2, {2.71828, "OK"});
mixedTable.addRow(3, {1.41421, "Warning"});
```

### Custom Row Tags

```cpp
struct DataPoint {
    double x, y;
};

std::ostream& operator<<(std::ostream& os, const DataPoint& p) {
    return os << "(" << p.x << ", " << p.y << ")";
}

TablePrinter<DataPoint, double> table("Point", {"f(x,y)"});
table.addRow({1.0, 2.0}, {func(1.0, 2.0)});
```

### Builder Pattern Chaining

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

---

## Best Practices

### 1. Use Appropriate Precision

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

```cpp
TablePrinter<int, double> table("Index", {"Value"});
table.reserve(10000);  // Pre-allocate for performance

for (int i = 0; i < 10000; ++i) {
    table.addRow(i, {compute(i)});
}
```

### 6. Error Handling

```cpp
try {
    table.addRow(x, {value1, value2});  // Wrong count throws
} catch (const std::invalid_argument& e) {
    std::cerr << "Column count mismatch: " << e.what() << '\n';
}

try {
    table.exportToFile("output.csv", ExportFormat::CSV);
} catch (const std::runtime_error& e) {
    std::cerr << "Export failed: " << e.what() << '\n';
}
```

---

## Export Format Details

### CSV (Comma-Separated Values)

```
header1,header2,header3
value1,value2,value3
```

**Use:** Excel, data analysis tools

### TSV (Tab-Separated Values)

```
header1	header2	header3
value1	value2	value3
```

**Use:** Databases, text processing

### Markdown

```
| Header1 | Header2 |
|---------|---------|
| Value1  | Value2  |
```

**Use:** GitHub README, documentation

### LaTeX

```latex
\begin{tabular}{|c|c|}
\hline
Header1 & Header2 \\
\hline
Value1 & Value2 \\
\hline
\end{tabular}
```

**Use:** Academic papers, reports

### HTML

```html
<table border="1">
  <thead><tr><th>Header1</th></tr></thead>
  <tbody><tr><td>Value1</td></tr></tbody>
</table>
```

**Use:** Web pages, reports

---

## Performance Considerations

- **Auto-width calculation**: O(rows × columns) - done once before printing
- **Memory**: ~40 bytes per cell + data size
- **Large tables**: Use `reserve()` to avoid reallocations
- **Export overhead**: Minimal, I/O-bound
- **Stream state**: RAII guard has negligible overhead

---

## Migration from Legacy API

### Old (Deprecated)

```cpp
ColDesc col("Name", 10, 3, 'F');  // Old style
VerticalVectorPrinter printer;   // Old class
```

### New (Recommended)

```cpp
ColumnFormat col("Name")
    .width(10)
    .precision(3)
    .format(FormatType::Fixed);

VectorTablePrinter printer;
```

---

## See Also

- **[Timer](Timer_ThreadPool.md#timer-class)**: Performance measurement for table generation
- **[Visualizer](Visualizer.md)**: Graphical data visualization
- **[Serializer](Serializer.md)**: Data export to files

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

**Key takeaway:** Professional-quality formatted output for console, files, and documents with minimal code.
