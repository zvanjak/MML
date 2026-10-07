# MML Object Persistence

This guide documents the MML object-persistence framework implemented under
`mml/tools/persistence/`. It covers the durable `.mmlj` JSON and `.mmlb` binary
formats used by generic `Persistence::Save` / `Persistence::Load` workflows.

For legacy visualizer exports such as sampled functions, curves, fields, ODE
solutions, and particle simulations, use `.mml` presentation files documented in
[SerializerFormats.md](SerializerFormats.md). Those files are stable visualizer
inputs, but they are not general object round-trip persistence files.

---

## Format Overview

| Format | Extension | Purpose | Human readable | Current object support |
|--------|-----------|---------|----------------|------------------------|
| Visualizer text | `.mml` | Presentation export for plotting and visualizers | Yes | Export-only sampled data |
| MML JSON object | `.mmlj` | Structured object persistence | Yes | `Vector`, `VectorN`, `Matrix`, `MatrixNM` |
| MML binary object | `.mmlb` | Compact object persistence | No | real/complex dynamic vectors and matrices |
| CSV/text helpers | `.csv`, `.txt`, `.dat` | Simple matrix I/O and inspection | Yes | existing matrix helper APIs |

Use `.mmlj` when you want readable object files, schema metadata, golden-file
fixtures, or easy debugging. Use `.mmlb` when you want compact, exact binary
payloads for real or complex vectors and dense matrices.

---

## Headers

Use the persistence umbrella header for the public surface:

```cpp
#include <mml/tools/Persistence.h>
```

The concrete implementation headers are:

```cpp
#include <mml/tools/persistence/PersistenceBase.h>
#include <mml/tools/persistence/JSON.h>
#include <mml/tools/persistence/BinaryBase.h>
#include <mml/tools/persistence/VectorBinary.h>
#include <mml/tools/persistence/MatrixBinary.h>
#include <mml/tools/persistence/VectorJSON.h>
#include <mml/tools/persistence/MatrixJSON.h>
```

---

## Basic Examples

### Vector JSON Round Trip

```cpp
#include <mml/tools/Persistence.h>

using namespace MML;

Vector<Real> v{1.25, -2.5, 3.75};

SerializeResult saved = Persistence::Save(v, "vector.mmlj");
if (!saved.success) {
    std::cerr << saved.message << '\n';
}

Vector<Real> loaded;
SerializeResult loadedResult = Persistence::Load("vector.mmlj", loaded);
```

The same operation can be performed explicitly:

```cpp
Persistence::SaveJson(v, "vector.mmlj");
Persistence::LoadJson("vector.mmlj", loaded);
```

### Vector Binary Round Trip

```cpp
Vector<double> v{1.0, 2.0, 3.0};

Persistence::Save(v, "vector.mmlb");

Vector<double> loaded;
Persistence::Load("vector.mmlb", loaded);
```

Binary dynamic vector files use the envelope magic `MML_VECT`.

Vectors and matrices support `float`, `double`, `long double`,
`std::complex<float>`, and `std::complex<double>`. Complex values encode each
real/imaginary component in little-endian order. The scalar metadata identifies
the exact stored type.

### Matrix JSON Round Trip

```cpp
Matrix<Real> A(2, 3);
A(0, 0) = 1.0;
A(0, 1) = 2.0;
A(0, 2) = 3.0;
A(1, 0) = 4.0;
A(1, 1) = 5.0;
A(1, 2) = 6.0;

Persistence::Save(A, "matrix.mmlj");

Matrix<Real> B;
Persistence::Load("matrix.mmlj", B);
```

Matrix JSON records `layout: row_major` and stores the matrix data as a flat
row-major array.

### Matrix Binary Round Trip

```cpp
Matrix<double> A(2, 2);
A(0, 0) = 1.0;
A(0, 1) = 2.0;
A(1, 0) = 3.0;
A(1, 1) = 4.0;

Persistence::Save(A, "matrix.mmlb");

Matrix<double> B;
Persistence::Load("matrix.mmlb", B);
```

Binary dynamic matrix files use the envelope magic `MML_MATX`.

### Fixed-Size JSON Types

`VectorN` and `MatrixNM` currently support JSON persistence:

```cpp
VectorN<Real, 3> direction{1.0, 0.0, 0.0};
Persistence::SaveJson(direction, "direction.mmlj");

VectorN<Real, 3> loadedDirection;
Persistence::LoadJson("direction.mmlj", loadedDirection);

MatrixNM<Real, 2, 3> T{1.0, 2.0, 3.0,
                       4.0, 5.0, 6.0};
Persistence::SaveJson(T, "transform.mmlj");
```

A fixed-size load rejects files whose `shape` does not match the compile-time
size.

### Presentation Export Example

A sampled ODE solution or function export is still a `.mml` presentation file:

```cpp
std::vector<std::string> legend = {"x", "v"};
Persistence::SaveODESolutionAsMultiFunc(solution, "Oscillator", legend,
                                       "oscillator.mml");
```

That file is intended for visualizer workflows. It is not the same contract as
`Persistence::Save(solution, "solution.mmlj")`, which is not implemented yet for
ODE solutions.

### Planned Sampled Function Persistence

Executable function bodies cannot be persisted. MML can persist sampled function
data: nodes, values, labels, domain metadata, scalar metadata, and optional source
annotations. Loading that data should reconstruct an interpolated function chosen
explicitly by the caller.

The intended API shape is:

```cpp
RealFunction f([](Real x) { return std::sin(x); });

Persistence::SaveSampledFunction(
    f,
    "sin.mmlj",
    Persistence::SampleGrid{0.0, Constants::PI, 64},
    Persistence::SaveOptions{}
);

auto linear = Persistence::LoadLinearFunction("sin.mmlj");
if (linear.result.success) {
  Real y = linear.function(Constants::PI / 4.0);
}

auto spline = Persistence::LoadSplineFunction("sin.mmlj");

auto polynomial = Persistence::LoadPolynomialFunction(
  "sin.mmlj",
  Persistence::PolynomialInterpolationOptions{.order = 5}
);
```

The important contract is that the load APIs do not recover executable code.
They load persisted sample data into a concrete interpolation model selected by
the user. The first implementation supports:

| Selection | Reconstructed type | Extra options |
|-----------|--------------------|---------------|
| `Linear` | `LinearInterpRealFunc` | optional extrapolation policy |
| `Polynomial` | `PolynomInterpRealFunc` | interpolation order / point count |
| `Spline` | `SplineInterpRealFunc` | optional endpoint derivatives |

The schema should store enough information to recreate the data independently of
the interpolation chosen at load time:

```json
{
  "mml": {
    "format": "MML_JSON",
    "version": 1,
    "object": "SampledRealFunction",
    "object_version": 1,
    "scalar": "Real",
    "scalar_bytes": 8,
    "scalar_encoding": "ieee754",
    "real_type": "double"
  },
  "domain": [0.0, 3.141592653589793],
  "nodes": [0.0, 0.5, 1.0],
  "values": [0.0, 0.4794255386, 0.8414709848],
  "metadata": {
    "label": "sin(x)",
    "source": "sampled from IRealFunction"
  }
}
```

Readers must validate matching node/value lengths, strictly valid node data for
the requested interpolator, scalar metadata, schema version, and interpolation
options before constructing the target interpolation object.

Implemented load helpers return `LoadFunctionResult<T>`, which contains both a
`SerializeResult` and the reconstructed interpolation object. Check
`result.success` before using `function`:

```cpp
auto loaded = Persistence::LoadSplineFunction("sin.mmlj");
if (loaded.result.success) {
  Real y = loaded.function(0.25);
}
```

Current sampled-function validation uses the same `SerializeError` vocabulary as
other persistence APIs:

| Condition | Error |
|-----------|-------|
| invalid sampling range or point count | `INVALID_PARAMETERS` |
| duplicate or non-monotonic nodes | `INVALID_PARAMETERS` |
| mismatched node/value lengths | `SCHEMA_MISMATCH` on load, `INVALID_PARAMETERS` on explicit save |
| wrong scalar metadata | `UNSUPPORTED_SCALAR` |
| unsupported schema version | `UNSUPPORTED_VERSION` |
| bad polynomial order | `INVALID_PARAMETERS` |

Existing `.mml` function visualizer exports remain unchanged and export-only.
Use `.mmlj` sampled-function persistence when you need to load the samples back
into an interpolation object.

---

## Options

### SaveOptions

```cpp
Persistence::SaveOptions options;
options.precision = 17;
options.pretty_json = true;
options.include_metadata = true;
options.allow_lossy = false;
options.metadata.title = "test fixture";
options.metadata.tags["source"] = "unit-test";
```

`precision` controls JSON numeric output. `pretty_json` controls indentation and
newlines. User metadata is available for object-persistence formats, but the
current dense vector/matrix serializers only emit required schema metadata.

### LoadOptions

```cpp
Persistence::LoadOptions options;
options.strict_schema = true;
options.allow_scalar_conversion = false;
options.max_allocation_bytes = 1ull << 30;

Matrix<Real> A;
Persistence::LoadJson("matrix.mmlj", A, options);
```

`max_allocation_bytes` is enforced before allocating vector or matrix storage.
This protects readers from unreasonable or corrupted file dimensions.

---

## JSON Schema

All MML JSON object-persistence files use a top-level object with an `mml`
header plus object-specific fields.

### Vector

```json
{
  "mml": {
    "format": "MML_JSON",
    "version": 1,
    "object": "Vector",
    "object_version": 1,
    "scalar": "Real",
    "scalar_bytes": 8,
    "scalar_encoding": "ieee754",
    "real_type": "double"
  },
  "shape": [3],
  "data": [1.25, -2.5, 3.75]
}
```

### VectorN

`VectorN<T,N>` uses the same payload shape as `Vector`, but the object kind is
`VectorN`. The loader checks that `shape[0] == N`.

```json
{
  "mml": {
    "format": "MML_JSON",
    "version": 1,
    "object": "VectorN",
    "object_version": 1,
    "scalar": "Real",
    "scalar_bytes": 8,
    "scalar_encoding": "ieee754",
    "real_type": "double"
  },
  "shape": [3],
  "data": [1, 0, 0]
}
```

### Matrix

```json
{
  "mml": {
    "format": "MML_JSON",
    "version": 1,
    "object": "Matrix",
    "object_version": 1,
    "scalar": "Real",
    "scalar_bytes": 8,
    "scalar_encoding": "ieee754",
    "layout": "row_major",
    "real_type": "double"
  },
  "shape": [2, 3],
  "data": [1, 2, 3, 4, 5, 6]
}
```

### MatrixNM

`MatrixNM<T,N,M>` uses the same schema as `Matrix`, but the object kind is
`MatrixNM`. The loader checks that `shape == [N, M]`.

### Required JSON Fields

| Field | Meaning |
|-------|---------|
| `mml.format` | Must be `MML_JSON` |
| `mml.version` | JSON container version, currently `1` |
| `mml.object` | `Vector`, `VectorN`, `Matrix`, or `MatrixNM` |
| `mml.object_version` | Object schema version, currently `1` |
| `mml.scalar` | Logical scalar name, such as `Real`, `float`, or `double` |
| `mml.scalar_bytes` | Stored scalar size in bytes |
| `mml.scalar_encoding` | Currently `ieee754` |
| `mml.real_type` | Concrete `Real` backing type when `scalar == Real` |
| `mml.layout` | Required for matrix types, currently `row_major` |
| `shape` | `[count]` for vectors, `[rows, cols]` for matrices |
| `data` | Flat numeric array matching `shape` |

The current JSON serializers support floating-point scalar types. Exact numeric
types such as rational and big integer are intentionally not part of this first
dense implementation.

---

## Binary Layout

All generic `.mmlb` files use a 48-byte envelope header followed by an
object-specific payload. Integer header and payload shape fields are encoded as
little-endian values. Floating-point payload values are stored as their IEEE 754
bit representation in little-endian byte order.

### Envelope Header

| Offset | Size | Field |
|--------|------|-------|
| 0 | 8 | magic, `MML_` plus object code, such as `MML_VECT` or `MML_MATX` |
| 8 | 2 | container version, currently `1` |
| 10 | 2 | header size, currently `48` |
| 12 | 4 | object kind enum |
| 16 | 4 | object schema version, currently `1` for dense vectors/matrices |
| 20 | 4 | scalar type enum |
| 24 | 4 | scalar byte size |
| 28 | 4 | endian marker `0x01020304` |
| 32 | 8 | payload byte count |
| 40 | 4 | flags, currently `0` |
| 44 | 4 | sidecar metadata policy enum |

### Object Codes

| Magic | Object |
|-------|--------|
| `MML_VECT` | dynamic real or complex `Vector<T>` |
| `MML_MATX` | dynamic real or complex `Matrix<T>` |

The envelope helpers also reserve codes for future `VectorN`, `MatrixNM`, sparse
matrices, tensors, ODE solutions, and sampled data. Those object payloads are
not implemented yet.

### Scalar Type Enum

| Enum | Meaning | Byte size |
|------|---------|-----------|
| `Float32` | `float` | 4 |
| `Float64` | `double` / current `Real` | 8 |
| `LongDouble` | `long double` canonical binary value | 24 |
| `ComplexFloat32` | `std::complex<float>` | 8 |
| `ComplexFloat64` | `std::complex<double>` | 16 |

The `LongDouble` payload does not copy ABI bytes. It stores value class, sign,
source precision, binary exponent, and a 128-bit significand in a fixed 24-byte
record. This avoids platform-specific padding and x87 representation details.
A reader rejects a finite value whose source precision exceeds the target
platform's `long double` precision rather than silently rounding it. Signed
zero, infinity, and NaN classification are preserved; NaN payload bits are not.

### Vector Payload

After the 48-byte envelope:

| Offset from payload start | Size | Field |
|---------------------------|------|-------|
| 0 | 8 | element count, `uint64` |
| 8 | `count * scalar_bytes` | scalar payload |

### Matrix Payload

After the 48-byte envelope:

| Offset from payload start | Size | Field |
|---------------------------|------|-------|
| 0 | 8 | row count, `uint64` |
| 8 | 8 | column count, `uint64` |
| 16 | `rows * cols * scalar_bytes` | row-major scalar payload |

---

## Validation And Troubleshooting

All new object-persistence APIs return `SerializeResult` and do not throw for
normal validation failures.

```cpp
auto result = Persistence::Load("matrix.mmlj", A);
if (!result.success) {
    std::cerr << SerializeErrorName(result.error) << ": "
              << result.message << '\n';
}
```

| Error | Common cause | Typical remedy |
|-------|--------------|----------------|
| `FILE_NOT_OPENED` | Path is wrong or output directory does not exist | Check path and permissions |
| `INVALID_FORMAT` | Wrong extension or invalid envelope magic | Use the correct file format |
| `UNSUPPORTED_VERSION` | Future container or object schema version | Upgrade MML or migrate the file explicitly |
| `TYPE_MISMATCH` | Loading `Matrix` data into `Vector`, or magic object code disagrees with header kind | Check target type and file identity |
| `ENDIAN_MISMATCH` | Binary endian marker is not `0x01020304` | Regenerate the file or add a future byte-swap path |
| `ALLOCATION_LIMIT_EXCEEDED` | Shape or payload size exceeds `LoadOptions::max_allocation_bytes`, overflows byte count, or exceeds supported dimensions | Raise the limit only for trusted files or fix corrupted metadata |
| `SCHEMA_MISMATCH` | Missing `mml` header, wrong `shape`, wrong `layout`, data length mismatch, nonzero binary flags | Inspect schema metadata and regenerate the file |
| `MALFORMED_INPUT` | Invalid JSON syntax or non-numeric dense data | Fix file syntax or data array values |
| `UNSUPPORTED_SCALAR` | Scalar metadata does not match target type or uses unsupported scalar encoding | Load into matching type or regenerate with supported scalar |
| `TRUNCATED_INPUT` | Binary file ends before the declared header, shape, or payload data | Recopy or regenerate the file |

### Allocation Limits

For untrusted files, set a smaller limit:

```cpp
Persistence::LoadOptions options;
options.max_allocation_bytes = 64 * 1024 * 1024;

Matrix<Real> A;
auto result = Persistence::LoadJson("matrix.mmlj", A, options);
```

The dense JSON and binary loaders validate shape and byte counts before resizing
`Vector` or `Matrix` storage.

### Scalar Mismatch

Scalar metadata is strict. A file with `scalar: "float"` is rejected when loading
into `Vector<Real>` or `Matrix<Real>`. Scalar conversion is intentionally not
performed in the first implementation.

### Version Policy

The current JSON container version is `1`. The current binary envelope version is
`1`. Readers reject unsupported future versions with `UNSUPPORTED_VERSION`.

---

## Testing Commands

Focused serializer tests:

```powershell
& .\build\tests\Release\MML_Tests.exe '[Serializer][VectorJSON],[Serializer][MatrixJSON],[Serializer][DenseBinary],[Serializer][DenseGuards]'
```

Single-header smoke test after regenerating `MML.h`:

```powershell
python tools\create_MML_single_header\create_header.py
cmake --build build --config Release --target MML_SingleHeader_Tests --parallel
& .\build\tests\Release\MML_SingleHeader_Tests.exe
```

Full suite:

```powershell
& .\build\tests\Release\MML_Tests.exe
```
