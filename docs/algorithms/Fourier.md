# Fourier Transforms

MML provides complex DFT and FFT transforms, real-input FFT support, cosine and sine
transforms, window functions, spectrum utilities, and convolution algorithms through
`mml/algorithms/Fourier/Fourier.h`.

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
