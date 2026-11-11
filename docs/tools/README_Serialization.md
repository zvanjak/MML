# MML Serialization & Data I/O Guide

Current entry points for object persistence, tabular data loading, and visualizer
exports. Detailed persistence and export schemas live in the linked guides;
the DataLoader examples below have checked, runnable counterparts.

## Table of Contents

1. [Overview](#overview)
2. [Quick Reference](#quick-reference)
3. [Object Persistence](#object-persistence)
4. [Matrix I/O](#matrix-io)
5. [Vector I/O](#vector-io)
6. [DataLoader](#dataloader)
7. [Serializer](#serializer)
8. [RigidBody Simulation Serializer](#rigidbody-simulation-serializer)
9. [File Format Specifications](#file-format-specifications)
10. [Best Practices](#best-practices)
11. [Runnable Examples](#examples)

---

## Overview

| Subsystem | Purpose | Direction | Header |
|-----------|---------|-----------|--------|
| **Persistence** | Durable vector/matrix object save/load | Bidirectional | `<mml/tools/Persistence.h>` |
| **DataLoader** | Column-oriented CSV, TSV and JSON datasets | Read-only | `<mml/tools/DataLoader.h>` |
| **Serializer** | Sample mathematical objects for visualizers | Write-only presentation export | `<mml/tools/Serializer.h>` |

Matrix and vector file I/O belongs to the free `Persistence` APIs, not members of
`Matrix` or `Vector`. Presentation exporters are free functions in
`MML::Serializer`, with names such as `SaveRealFunc`, not `SaveRealFunction`.

Visualizer `.mml` files are **presentation exports**: they do not reconstruct
arbitrary functions, fields or solver instances. Use explicit `.mmlj` JSON and
`.mmlb` binary formats for durable object persistence. Count-prefixed object text
can also use `.mml`; that layout is not a visualizer presentation format.

### Supported File Formats

| Format | Extension | Purpose |
|--------|-----------|---------|
| MML JSON object | `.mmlj` | Readable object persistence |
| MML binary object | `.mmlb` | Canonical vector/matrix binary persistence |
| Object text | `.txt`, `.dat` | Count-prefixed vector or dimension-prefixed matrix data |
| CSV | `.csv` | Tabular datasets or numeric matrix I/O |
| TSV | `.tsv`, `.tab` | Tab-separated datasets |
| JSON | `.json` | DataLoader arrays of flat objects |
| Visualizer text | `.mml` | Sampled functions, curves, fields and simulation data |

An ordinary DataLoader JSON array is not an MML `.mmlj` object schema.

---

## Quick Reference

| Goal | Current API |
|------|-------------|
| Save/load a vector or matrix | `Persistence::Save`, `Persistence::Load` |
| Load a file by extension | `Data::LoadData(filename)` |
| Load a file with structured errors | `Data::LoadData(filename, output)` |
| Parse in-memory CSV/TSV | `Data::LoadFromCSVString(content, delimiter, hasHeader)` |
| Parse in-memory JSON | `Data::LoadFromJSONString(content)` |
| Parse a string in an explicit format | `Data::LoadFromString(content, format, hasHeader)` |
| Configure CSV schema/type inference | `Data::CSVLoadOptions` and `Data::LoadCSV` |
| Read a column | `dataset["x"]`, `dataset.GetRealColumn("x")` |
| Export a sampled real function | `Serializer::SaveRealFunc` |
| Export an ODE solution | `Serializer::SaveODESolutionAsMultiFunc` |

---

## Object Persistence

Include `<mml/tools/Persistence.h>` and use `MML::Persistence::Save` / `Load`.
These APIs return `IOResult`; always check it before using the output. Filename
overloads can infer the format from the extension, while stream overloads require
an explicit `Persistence::Format`.

Dynamic vectors and matrices support text, JSON and binary storage. Matrix CSV is
also supported; fixed-size `VectorN` and `MatrixNM` support JSON. Scalar support,
conversion rules, allocation limits, and transactional-load guarantees depend on
the selected object and format.

See [SerializationPersistence.md](SerializationPersistence.md) for exact
signatures, schemas, options, failure contracts and runnable examples.

## Matrix I/O

Use the free persistence APIs with a `Matrix<T>`:

- `.txt` / `Persistence::Format::Text`: dimensions followed by matrix entries.
- `.csv` / `Persistence::Format::CSV`: matrix CSV, without named dataset columns.
- `.mmlj` / `Persistence::Format::JSON`: structured object schema.
- `.mmlb` / `Persistence::Format::Binary`: canonical binary object storage.

The old matrix `SaveToFile`, `SaveToCSV`, `LoadFromFile` and `LoadFromCSV` members
are not the current API. For named or mixed-type CSV columns, use DataLoader
instead of the matrix reader. DataLoader treats the first record as a header by
default; matrix CSV has no such named-column header.

Default text/CSV precision is the scalar component's `max_digits10`, not a fixed
15 decimal places. Positive `Persistence::SaveOptions::precision` requests a
specific number of significant digits and can lose information.

The checked matrix text example is in
[docs_demo_serializer.cpp](../../src/docs_demos/tools/docs_demo_serializer.cpp);
run `MML_DocsApp --persistence-portability`.

## Vector I/O

Use the same free persistence APIs with a `Vector<T>`. Text stores the element
count followed by values; JSON and binary use their explicit object schemas.
Vector CSV is not implemented. The removed vector `SaveToFile` / `LoadFromFile`
members have no compatibility aliases.

Run `MML_DocsApp --vector-persistence` for checked file/stream round trips,
formatting restoration, malformed-input diagnostics and payload budgets.

---

## DataLoader

### Header

```cpp
#include <mml/MMLBase.h>
#include <mml/tools/DataLoader.h>
#include <iostream>
using namespace MML;
using namespace MML::Data;
```

Run **`MML_DocsApp --dataloader`** for the examples in this section, including
Unicode, numeric inference and CSV-record checks. Their source is
[docs_demo_data_loader.cpp](../../src/docs_demos/tools/docs_demo_data_loader.cpp).
The file-loading examples create, validate and clean up their own temporary
fixtures; they do not depend on a working-directory `sample_data.csv`.

### Core Types and Column Access

`Dataset` is **column-oriented**. Its `columns` member is a
`std::vector<DataColumn>`, not a list of column names; there is no `Dataset.rows`,
`DataRow` or `DataValue` API.

| Type/member | Meaning |
|-------------|---------|
| `Dataset::name`, `rowCount`, `columns` | Dataset metadata and typed columns |
| `NumRows()`, `NumColumns()`, `GetColumnNames()` | Shape and column names |
| `dataset[index]`, `dataset[name]` | A `DataColumn`; invalid lookup throws `DataError` |
| `HasColumn(name)` | Check a column name before lookup |
| `GetRealColumn(name)` | Copy a `REAL` column or convert an `INT` column to `Vector<Real>` |
| `GetIntColumn`, `GetBoolColumn`, `GetStringColumn` | Typed column extraction; string extraction uses display conversion |
| `DataColumn::name`, `type` | Name and `ColumnType` tag |
| `realData`, `intData`, `boolData`, `stringData`, `dateData`, `timeData` | Storage selected by the column type |
| `Size()`, `IsMissing(row)`, `CountMissing()` | Column length and missing-value mask |
| `GetReal(row)`, `GetAsString(row)` | Numeric conversion or display access to a cell |

`ColumnType` includes `REAL`, `INT`, `BOOL`, `STRING`, `DATE`, `TIME`, `DATETIME`
and `MISSING`. Check the tag before reading typed storage. `GetRealColumn`
accepts only `REAL`/`INT`; `GetReal` converts `REAL`/`INT`/`BOOL` cells and returns
NaN for missing values or other types. `GetAsString` is display formatting, not
lossless numeric persistence.

```cpp
std::string csv = "x,y,label\n1.0,2.0,A\n2.0,4.0,B\n3.0,6.0,C\n4.0,8.0,D\n5.0,10.0,E\n";
Dataset data = LoadFromCSVString(csv, ',', true);
const DataColumn& xColumn = data["x"];
const DataColumn& labelColumn = data["label"];
for (size_t row = 0; row < data.NumRows(); ++row) {
    if (xColumn.IsMissing(row) || xColumn.GetReal(row) != static_cast<Real>(row + 1))
        throw DataError("Column access demo lost a numeric value");
    std::cout << labelColumn.GetAsString(row) << ": " << xColumn.GetReal(row) << '\n';
}
Vector<Real> xVec = data.GetRealColumn("x");
Vector<Real> yVec = data.GetRealColumn("y");
```

### File Loading Versus String Parsing

| Input | Throwing dataset-returning API | Result-based API |
|-------|--------------------------------|------------------|
| CSV filename | `LoadCSV(filename)` or `LoadCSV(filename, options)` | `LoadCSVSafe(filename, output, hasHeader, delimiter, inferTypes)` or `LoadCSV(filename, output, options)` |
| TSV filename | `LoadTSV(filename, hasHeader, inferTypes)` | `LoadTSVSafe(filename, output, hasHeader, inferTypes)` |
| JSON filename | `LoadJSON(filename)` | `LoadJSONSafe(filename, output)` |
| Filename with detected format | `LoadData(filename, hasHeader)` | `LoadData(filename, output, hasHeader)` |

`LoadData` recognizes `.csv`, `.tsv`, `.tab` and `.json`, using ASCII-only
case-insensitive extension matching; other/missing extensions default to CSV.
The header flag applies to CSV/TSV, not JSON.

For **in-memory content**, use `LoadFromCSVString(content, delimiter, hasHeader)`,
`LoadFromJSONString(content)`, or `LoadFromString(content, DataFormat, hasHeader)`.
TSV strings use the CSV parser with `'\t'`. These functions throw `DataError` on
invalid input.

```cpp
const std::string goodCSV = "a,b,c\n1,2,3\n4,5,6\n";
const Dataset parsed = LoadFromCSVString(goodCSV);
```

There is **no `fromString` flag** on a safe file loader. Do not pass CSV content
as its filename. In particular, the string parser's argument order
`(content, delimiter, hasHeader)` differs from the file loader's options.
Use the `Dataset&` output overload rather than deprecated result-owning
`LoadCSVSafe` / `LoadTSVSafe` / `LoadJSONSafe` overloads.

### Safe File Loading and Error Reporting

`IOResult` reports `success`, `error`, `message`, `filename`, and optional
`offset` / `offset_kind`; it does not own the dataset. The demo creates
`csvFilename` before executing this fragment:

```cpp
Dataset loaded;
const IOResult result = LoadCSVSafe(csvFilename, loaded, true, ',');
if (!result)
    throw FileIOError(result.message, result.filename);
std::cout << "CSV file load: " << IOErrorName(result.error) << '\n';
```

The checked counterpart verifies two rows, three integer columns, all six values,
the success code and filename. It also checks CSV/TSV/JSON convenience loaders,
both `LoadData` overloads, and explicit-format string parsing.

A failed result-based load leaves the output dataset unchanged. In the demo,
`missingFilename` names an absent file in its own newly created directory:

```cpp
const IOResult missing = LoadCSVSafe(missingFilename, loaded);
if (missing || missing.error != IOError::FILE_NOT_OPENED ||
    missing.filename != missingFilename || missing.message.empty())
    throw DataError("Missing CSV did not report FILE_NOT_OPENED");
std::cout << "Expected missing-file failure: "
          << IOErrorName(missing.error) << ": " << missing.message << '\n';
```

The demo also verifies `MALFORMED_INPUT` for malformed quotes, the physical line
diagnostic, and preservation of the previously loaded dataset on both failures.
Unexpected outcomes throw rather than silently printing a success-shaped result.

### CSVLoadOptions and Headerless Data

The actual options type is `CSVLoadOptions`, with these fields:

| Field | Default | Purpose |
|-------|---------|---------|
| `hasHeader` | `true` | First logical record supplies column names |
| `delimiter` | `','` | CSV/TSV field separator |
| `inferTypes` | `true` | Infer typed columns; `false` retains strings |
| `strictSchema` | `false` | Reject row-width mismatches |
| `warnings` | `nullptr` | Optional `std::vector<std::string>*` for nonfatal width diagnostics |

There is no DataLoader `LoadOptions` or `columnNames` option. Without a header,
names are `Column0`, `Column1`, etc.; rename the resulting `DataColumn::name`
fields if required. The demo creates `noHeaderFilename` containing
`1.25;alpha` and `2.5;beta`, then executes:

```cpp
CSVLoadOptions options;
options.hasHeader = false;
options.delimiter = ';';
options.inferTypes = false;
options.strictSchema = true;
Dataset strings = LoadCSV(noHeaderFilename, options);
strings.columns[0].name = "x";
strings.columns[1].name = "label";
```

Its checks verify both rows, generated names before renaming, string storage,
and the result-based `LoadCSV(filename, output, options)` overload.
`Persistence::LoadOptions` is a different type for object persistence.

### ASCII Tokens and Unicode Text

File extensions, missing-value markers, and boolean aliases are matched
case-insensitively using **ASCII-only** folding, independent of the C character
locale. Non-ASCII bytes are never passed to character case-conversion functions
or changed by token normalization. This does not perform Unicode case folding,
normalization, or validation: Greek capital/small sigma, for example, remain
different text. Existing whitespace trimming and CSV/JSON decoding still apply.

UTF-8 string values and column names retain their bytes. JSON Unicode escapes
are decoded to UTF-8 before type inference, so a string such as `N\u00c4N`
is text, not the ASCII missing marker `NaN`.

```cpp
const Dataset data = LoadFromJSONString(
    R"([{"name":"Caf\u00e9","active":"TrUe"},{"name":"N\u00c4N","active":"FaLsE"}])");
std::cout << data["name"].stringData[0] << '\n';
```

This snippet and explicit byte/boolean checks also run independently through
`MML_DocsApp --dataloader-unicode`.

### Numeric Inference and Conversion Errors

The DataLoader keeps its 90% majority inference threshold, excluding recognized
missing values. Within a numeric column, any decimal/exponent-spelled token or
integer outside `int` range promotes the column to `REAL`, even when most rows
contain small integers. For example:

```cpp
const std::string csv = "value\n1\n2\n3\n4\n5\n6\n7\n8\n9\n10.75\n";
Dataset data = LoadFromCSVString(csv);
std::cout << data["value"].realData[9];  // 10.75, not 10
```

This example is exercised by `demo_numeric_import_contract` in the
[DataLoader documentation demo](../../src/docs_demos/tools/docs_demo_data_loader.cpp).

Numeric conversions must consume the entire trimmed token and fit the selected
storage. Integer-spelled tokens (an optional sign followed only by digits) must
be represented exactly: a value outside `int` range may use `Real`, but an
integer that `Real` would round is rejected. Decimal/exponent-spelled tokens
retain normal floating-point rounding and the existing range checks.

Invalid present values are **not missing cells**. Once a numeric/bool type has
been inferred, an invalid conversion fails the entire load with a diagnostic
containing the value, one-based data-row number, column number and column name.
Recognized missing markers such as `NA` and JSON `null` remain missing.

The `IOResult` file overloads report `MALFORMED_INPUT`, attach the filename and
leave the caller's `Dataset` unchanged; CSV/TSV failures also carry the physical
line in `offset` with `IOOffsetKind::Line`. JSON conversion diagnostics identify
the data row in the message, not a fabricated byte/line offset. Dataset-returning
file overloads and string loaders throw `DataError`. Setting
`CSVLoadOptions::inferTypes` to `false` explicitly retains string storage.

This conversion policy does not change CSV record/quote parsing.

### CSV/TSV Format Requirements

CSV and TSV loaders share a strict, record-oriented quote parser. A quoted field
may contain delimiters, doubled quotes (`""` represents one quote), and physical
newlines. Quotes may open only at the beginning of a field; after the closing
quote only a delimiter or record end is allowed. Spaces/tabs outside a quoted
field are not accepted as padding.

Embedded LF and CRLF bytes are preserved, including blank lines inside quoted
fields. File loaders use binary input for the same results as string loaders.
An initial UTF-8 BOM is accepted, blank physical lines outside quoted fields are
skipped, and trailing empty fields are retained. Existing value trimming and
missing-marker rules still apply after CSV decoding, even to quoted values.

```cpp
const std::string csv = "name,value\r\n\"hello\r\nworld\",1\r\n";
Dataset data = LoadFromCSVString(csv);
std::cout << data.NumRows() << '\n';              // 1 logical data row
std::cout << data["name"].stringData[0] << '\n';  // contains the original CRLF
```

This example is checked by `demo_csv_record_contract` in the
[DataLoader documentation demo](../../src/docs_demos/tools/docs_demo_data_loader.cpp).

Malformed quote syntax fails the entire load: result-based file APIs return
`MALFORMED_INPUT` without replacing the destination, and dataset-returning APIs
throw `DataError`. Syntax diagnostics identify the physical line and field;
file results also carry the line in `offset` with `IOOffsetKind::Line`. An
unclosed field is located at its opening quote. Delimiters that are themselves
quotes or newlines are rejected with `INVALID_PARAMETERS`.

`CSVLoadOptions::strictSchema` still rejects width mismatches with
`SCHEMA_MISMATCH`; its offset is the logical record's starting physical line.
Non-strict warning messages retain logical `Row N` numbering (including the
header, if present), not physical line numbering.

### JSON Format Requirements

JSON input must be an **array of flat objects**. Numbers, strings (including
escapes), booleans and `null` are supported. Nested objects, nested arrays and
top-level objects are not supported. See `demo_json_loading` for the runnable
array-of-objects example.

---

## Serializer

Include `<mml/tools/Serializer.h>` and call functions in namespace `MML::Serializer`.
These are write-only presentation exports returning `IOResult`, not object
round-trip persistence.

| Object | Current function |
|--------|----------------|
| Real function | `Serializer::SaveRealFunc` |
| Multiple functions | `Serializer::SaveRealMultiFunc` |
| Parametric curve | `Serializer::SaveParamCurveCartesian2DResult`, `SaveParamCurveCartesian3DResult` |
| Scalar field | `Serializer::SaveScalarFunc2DCartesian`, `SaveScalarFunc3DCartesian` |
| Vector field | `Serializer::SaveVectorFunc2DCartesian`, `SaveVectorFunc3DCartesian` |
| ODE solution | `Serializer::SaveODESolutionAsMultiFunc` |
| Particle simulation | `Serializer::SaveParticleSimulation2D`, `SaveParticleSimulation3D` |

See the modules included by [Serializer.h](../../mml/tools/Serializer.h) for exact
argument order, and [SerializerFormats.md](SerializerFormats.md) for the actual
`MML_`-prefixed identifiers and versioned headers. The old `# Format: ...` layouts
are not the current exporter header specification.

## RigidBody Simulation Serializer

Rigid-body simulation and its specialized serializer belong to **MPL**, not MML
core. Use the MPL project's documentation and demos for that API; MML's particle
export methods do not substitute for a rigid-body simulator.

## File Format Specifications

- [SerializationPersistence.md](SerializationPersistence.md): object text/JSON,
  binary envelopes, supported scalars and conversion rules.
- [SerializerFormats.md](SerializerFormats.md): visualizer presentation formats.
- [Visualizers.md](Visualizers.md): consuming visualizer exports.

Do not treat DataLoader JSON arrays, MML JSON object schemas, and sampled
visualizer files as interchangeable representations.

## Best Practices

1. **Check the result.** Persistence, presentation exporters and result-based
   DataLoader file APIs return `IOResult`; report its error code, message and
   filename. Check offset kind before interpreting an offset as a line or byte.
2. **Choose the correct input API.** File loaders take paths, string parsers take
   content, and CSV options belong to `CSVLoadOptions`.
3. **Validate the schema and missing values.** Use `HasColumn`, dimensions,
   `ColumnType` and `IsMissing` before numerical processing. Enable `strictSchema`
   when ragged records must be rejected.
4. **Preserve precision deliberately.** Use persistence round-trip defaults for
   stored vectors/matrices, rather than display conversions or arbitrary
   decimal-place assumptions.
5. **Keep presentation separate from persistence.** Sampled `.mml` output is for
   visualization; use the object schemas for durable reconstruction.

## Examples

| Command | Checked counterpart |
|---------|---------------------|
| `MML_DocsApp --dataloader` | File/string loading, automatic detection, column access, CSV options, success/failure diagnostics, destination preservation, numeric/CSV/Unicode contracts |
| `MML_DocsApp --dataloader-unicode` | Focused ASCII-token and UTF-8 checks |
| `MML_DocsApp --vector-persistence` | Vector file/stream round trips and failure contracts |
| `MML_DocsApp --persistence-portability` | Matrix text precision and persistence portability |
| `MML_DocsApp --serializer-sampling` | Finite sample coordinates and exported coordinate round trips |

Demo sources:
[DataLoader](../../src/docs_demos/tools/docs_demo_data_loader.cpp) and
[persistence/serializer](../../src/docs_demos/tools/docs_demo_serializer.cpp).
These selectors do not launch a graphical visualizer.
