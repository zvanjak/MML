# Fourier Transforms

MML provides complex DFT and FFT transforms, real-input FFT support, cosine and sine
transforms, window functions, spectrum utilities, and convolution algorithms through
`mml/algorithms/Fourier/Fourier.h`.

## Spectrum ownership

`mml/algorithms/Fourier/FourierSpectrum.h` owns operations on transform output:
`Magnitude`, `MagnitudeSquared`, `Phase`, `UnwrapPhase`, `FrequencyAxis` and its
layout enum, `FFTShift`, and `IFFTShift`. `MagnitudeSquared` is raw `|X[k]|^2`,
not power per Hz or a normalized one-sided spectrum. The redundant `Power`
spelling is removed without a forwarding alias.
These core spectrum utilities are also included in the generated single-header
distribution.

```cpp
#include <mml/algorithms/Fourier/FourierSpectrum.h>

Vector<Complex> coefficients{Complex(3, 4), Complex(0, 2)};
const auto squared = MML::Fourier::MagnitudeSquared(coefficients);
const auto frequencies = MML::Fourier::FrequencyAxis(
    8, Real(80), MML::Fourier::FrequencyAxisLayout::OneSided);
```

Runnable counterpart: `Docs_Demo_SpectrumUtilities()` in
[docs_demo_fourier.cpp](../../src/docs_demos/algorithms/docs_demo_fourier.cpp);
run `MML_DocsApp --spectrum-utilities`.

Native `fourier/Spectrum.h` in MML-Packages-Private owns windowed periodograms,
Welch averaging, STFT/power spectrograms, frame-start timing, and dB presentation
policy. It consumes these core utilities rather than supplying another frequency
axis implementation. This separation introduces no package dependency into core.

## FFT sizes

`FFT::Forward`, `FFT::Inverse`, and `FFT::Transform` accept every non-empty input size.
Power-of-two lengths use the iterative radix-2 Cooley-Tukey kernel. Other lengths use
Bluestein's chirp-z algorithm, reducing the transform to a power-of-two convolution.
Both paths have $O(N \log N)$ time complexity.

Plain transforms and default Detailed calls preserve the input length. Zero-padding
is not required and should only be used when the application specifically wants a
different frequency grid.

`FFT::ForwardDetailed` and `FFT::InverseDetailed` accept explicitly gated padding:
set `allow_padding=true`. A positive `target_size` must be at least the input length;
zero selects the next power of two. Arbitrary larger targets are supported through
Bluestein. Samples/bins are extended with zeros at the end, and normalization uses
the resulting length. Padding an inverse spectrum does **not** reorder bins, split
Nyquist coefficients, or perform signal resampling. An equal-size target does not
count as applied padding. Empty inputs, truncation, negative targets, and a target
without `allow_padding` are rejected. DFT, DCT and DST Detailed APIs reject all padding
requests rather than silently ignoring them.

```cpp
Vector<Complex> signal{Complex(1, 0), Complex(2, 0), Complex(3, 0)};
MML::Fourier::FourierConfig config;
config.allow_padding = true;
config.target_size = 8;
config.normalization = MML::Fourier::TransformNormalization::Orthonormal;
auto spectrum = MML::Fourier::FFT::ForwardDetailed(signal, config);
```

Runnable counterpart: `Docs_Demo_FourierDetailedContracts()` in
[docs_demo_fourier.cpp](../../src/docs_demos/algorithms/docs_demo_fourier.cpp);
run `MML_DocsApp --fourier-contracts`.

`DFT` remains the direct $O(N^2)$ reference implementation. It is useful for small
inputs and numerical verification.

## Normalization

All complex DFT and FFT entry points support `TransformNormalization`:

| Mode | Forward scale | Inverse scale |
|---|---:|---:|
| `Legacy` | $1$ | $1/N$ |
| `None` | $1$ | $1$ |
| `Forward` | $1/N$ | $1$ |
| `Inverse` | $1$ | $1/N$ |
| `Orthonormal` | $1/\sqrt{N}$ | $1/\sqrt{N}$ |

Matching forward and inverse calls round-trip for every mode except `None`, whose
round trip is scaled by $N$.

### Large finite complex transforms

Complex `DFT` and `FFT` operations, their Detailed/real-input convenience
overloads, and `RealFFT::Forward` / `Inverse` protect large finite inputs with
power-of-two headroom before summations, butterflies and real-spectrum
reconstruction. Bluestein's bound includes its convolution FFTs and products,
not only the requested transform length. Normalization is combined with
restoring the input exponent, so a representable normalized result does not
require an overflowing unnormalized intermediate or reciprocal input scale.

```cpp
const Real huge = std::numeric_limits<Real>::max() * Real(0.75);
const Vector<Complex> signal(2, Complex(huge, 0));
const auto spectrum = FFT::Forward(signal, TransformNormalization::Forward);
// spectrum[0] == Complex(huge, 0), spectrum[1] == Complex(0, 0)
```

Runnable counterpart: `Docs_Demo_FourierRange()` in
[docs_demo_fourier.cpp](../../src/docs_demos/algorithms/docs_demo_fourier.cpp);
run `MML_DocsApp --fourier-range`. It also checks a three-point Bluestein
transform, the default normalized inverse DFT, and normalized real-input FFT.

The safeguard changes neither the normalization table nor Bluestein's
convolution normalization. `FFT::Transform` remains unnormalized in both
directions. Detailed padding uses the effective padded length for both scaling
and work metadata. Ordinary and tiny inputs are not predivided by the transform
length, avoiding unnecessary subnormal underflow.

This is finite-range protection, not extended-precision arithmetic: cancellation,
roundoff and underflow of very small components in extreme mixed-scale inputs
remain possible. It does not extend the DCT/DST range contract or promise to
represent genuinely overflowing coefficients. With the default finite checks,
Detailed APIs still report `NumericalInstability` (or propagate
`NumericalMethodError`) for nonfinite outputs and reject nonfinite inputs.
Disabling `check_finite` retains the existing unchecked behavior.

## DCT dispatch

`DCT::ForwardII` and `DCT::InverseII` use direct formulas for small inputs and the
FFT-based $O(N \log N)$ implementation for larger inputs. The fast path supports both
power-of-two and arbitrary lengths.

The plain DCT-II/III and DST-I entry points also accept a normalization argument;
their Detailed counterparts use `FourierConfig::normalization`. Omitting it preserves
the existing `Legacy` coefficients.

### Window and verification edge contracts

`Windows::Kaiser(1, beta)` and `Windows::Gaussian(1, sigma)` return a single
unit weight, consistent with the other singleton windows. Parameters are still
validated first: beta must be finite and nonnegative; sigma must be finite and
positive. Empty windows and invalid parameters raise `FourierError`.

`DCT::VerifyRoundTrip` and `DCT::VerifyFastMatchesReference` require nonempty,
finite data and a finite, nonnegative **absolute** tolerance. Empty input and
invalid tolerance raise `FourierError`; nonfinite samples raise
`NumericInputError`. A numerical discrepancy or nonfinite computed output
returns false. Zero tolerance requests an exact comparison.
Fast/reference verification actually executes both implementations at every
supported nonempty length, including singleton, odd, prime and non-power-of-two
lengths; there is no skipped-validation success path. The transform work-size
and resource checks still apply and may throw. These helpers verify the legacy
DCT-II/III normalization, not every normalization convention.

### General phase unwrapping

`UnwrapPhase` accepts arbitrary finite radians, not only the principal
`[-pi, pi]` output of `Phase`. It preserves the first sample. For each adjacent
input jump strictly above the positive `discontinuity` threshold, it removes
the smallest positive number of whole `2*pi` turns needed to bring that jump
at or below the threshold; a jump below the negative threshold is corrected
symmetrically. Exact threshold ties are unchanged. With the default threshold
`pi`, corrected jumps lie in `[-pi, pi]`; positive/negative odd-half-turn ties
retain their sign. Thus `[0, 5*pi]` becomes `[0, pi]`, not `[0, 3*pi]`.

Custom thresholds retain their previous meaning. In particular, a threshold
below `pi` can produce an opposite-sign jump larger than that threshold
(there may be no congruent jump inside such a narrow interval). A threshold
above `pi` need not select the nearest principal jump. The operation uses the
floating-point period `2 * Constants::PI`; very large phases are reduced
before subtraction to avoid overflow and loss of whole turns. Returned
values remain subject to ordinary `Real` rounding.

Input must be nonempty and the threshold finite and positive (`FourierError`);
nonfinite phase samples raise `NumericInputError`. An unrepresentable
unwrapped result raises `NumericalMethodError`. `FrequencyAxis` likewise
requires a finite positive sample rate through the shared scalar validator.

```cpp
const auto kaiser = Windows::Kaiser(1, Real{8.6});   // {1}
const auto gaussian = Windows::Gaussian(1, Real{0.4}); // {1}
const Vector<Real> samples{1, 2, 3};
const bool roundTrip = DCT::VerifyRoundTrip(samples, Precision::EqualityTolerance);
const bool fastMatches = DCT::VerifyFastMatchesReference(samples, Precision::EqualityTolerance);
const auto unwrapped = UnwrapPhase(Vector<Real>{0, Real{5} * Constants::PI}); // {0, pi}
```

Runnable counterpart: `Docs_Demo_FourierDetailedContracts()` in
[docs_demo_fourier.cpp](../../src/docs_demos/algorithms/docs_demo_fourier.cpp);
run `MML_DocsApp --fourier-contracts`.

Define the raw cosine pair as

$$C_k=2\sum_{j=0}^{N-1}x_j\cos\frac{\pi k(j+1/2)}N,\qquad
I_j=c_0+2\sum_{k=1}^{N-1}c_k\cos\frac{\pi k(j+1/2)}N.$$

The raw sine transform is

$$S_k=2\sum_{j=0}^{N-1}x_j\sin\frac{\pi(k+1)(j+1)}{N+1}.$$

Their logical normalization lengths are $L=2N$ for DCT and $L=2(N+1)$ for DST,
matching FFTW/SciPy's unnormalized conventions:

| Mode | Forward | Inverse |
|---|---|---|
| `None` | raw | raw |
| `Forward` | raw / $L$ | raw |
| `Inverse` | raw | raw / $L$ |
| `Orthonormal` | orthonormal basis | transposed orthonormal basis |
| `Legacy` DCT | $C/N$ | $I/2$ |
| `Legacy` DST | $S/\sqrt L$ | $S/\sqrt L$ |

Orthonormal DCT uses $C_0/(2\sqrt N)$ for DC and $C_k/\sqrt{2N}$ otherwise.
Its inverse is $c_0/\sqrt N+\sqrt{2/N}\sum_{k>0}c_k\cos(\ldots)$.
Orthonormal DST is $S/\sqrt L$ in both directions. Thus DCT of eight ones has
DC coefficient $\sqrt8$ in `Orthonormal` mode, but retains DC coefficient $2$
in `Legacy` mode. Matching modes round-trip unchanged except `None`, whose
round trip is multiplied by $L$. Unknown normalization enum values are errors.

## Detailed validation and work metadata

All Detailed transforms check every real and imaginary input component, and all
returned output components, when `check_finite=true` (the default). Under
`ConvertToStatus`, nonfinite input and invalid options report `InvalidInput`;
nonfinite computed output reports `NumericalInstability`; unsupported workspace
sizes and allocation failures report `ResourceLimitExceeded`. The default
`Propagate` policy throws the corresponding exception instead. Setting
`check_finite=false` explicitly opts out of both input and output checks and may
return nonfinite values. Finite checking diagnoses overflow; it does not guarantee
that every mathematically representable normalized transform avoids intermediate
overflow. No error bound is exposed, and the inherited `estimate_error` flag does
not change this.

`input_size` is the original length; `output_size` is the produced length.
`normalization` records the requested convention, including `Legacy`, and
`padding_applied` indicates actual zero extension. Failure results retain available
metadata and elapsed time; check `IsSuccess()` before using their payload.

`algorithm_name` keeps the stable operation name. `kernel_name` distinguishes
`Direct`, `Radix2`, `Bluestein`, `FFT/Radix2`, and `FFT/Bluestein`.
`work_size` is the direct length or largest FFT length, including Bluestein's
convolution workspace. DCT uses direct sums below 32 samples and a length-$2N$
FFT otherwise, regardless of whether $N$ itself is a power of two.

`work_units` is a deterministic **kernel-work metric**, not measured callbacks,
time or FLOPs:
- Direct DFT/DCT/DST: $N^2$ summation terms.
- Radix-2 FFT: $(N/2)\log_2N$ butterflies.
- Bluestein: three radix-2 FFTs of length $M$, plus $M$ convolution products and
  $2N$ chirp products, where $M$ is the next power of two at least $2N-1$.
- Fast DCT: the length-$2N$ FFT's metric plus $N$ twiddle products.

Scaling, copying and finite-validation scans are excluded. The metric uses a
64-bit unsigned counter, and doubled/convolution work sizes are checked before
allocation. `function_evaluations` is now zero: these transforms invoke no user
function. Consumers of the former complexity surrogate in that field should use
`work_units` instead.
