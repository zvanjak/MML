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

The result always has the same length as the input. Zero-padding is therefore not
required and should only be used when the application specifically wants a different
frequency grid.

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

## DCT dispatch

`DCT::ForwardII` and `DCT::InverseII` use direct formulas for small inputs and the
FFT-based $O(N \log N)$ implementation for larger inputs. The fast path supports both
power-of-two and arbitrary lengths.
