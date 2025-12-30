# Fourier Analysis and Transforms

Comprehensive guide to Fourier transforms, spectral analysis, convolution, and frequency-domain signal processing.

## Overview

Fourier analysis decomposes signals into frequency components, enabling:
- **Signal Processing**: Filtering, noise reduction, feature extraction
- **Spectral Analysis**: Frequency content, power spectra, time-frequency analysis
- **Convolution/Correlation**: Fast algorithms via FFT, signal matching
- **Function Approximation**: Chebyshev via DCT, spectral methods for PDEs
- **Image/Audio Compression**: JPEG (DCT), MP3 (MDCT), video codecs

The library provides **8 algorithm classes** covering discrete transforms (DFT/FFT/DCT), applications (convolution, correlation, spectra), and continuous expansions (Fourier series).

## Quick Reference

| Class | Purpose | Key Methods | Performance | Use Case |
|-------|---------|-------------|-------------|----------|
| **DFT** | Naive discrete Fourier transform | Forward, Inverse | O(N²) | **Testing, reference** |
| **FFT** | Fast Fourier transform | Forward, Inverse, Transform | O(N log N) | **Production FFT** |
| **RealFFT** | Optimized real-valued FFT | Forward, Inverse | 2× faster | **Real signals** |
| **DCT** | Discrete Cosine Transform | ForwardII, InverseII, ForwardDST | O(N²) naive | **JPEG, Chebyshev** |
| **Convolution** | Fast convolution via FFT | Linear, Circular | O(N log N) | **Filtering, systems** |
| **Correlation** | Cross/auto-correlation | Cross, Auto, Normalized | O(N log N) | **Signal matching** |
| **PowerSpectrum** | Spectral analysis | Compute, Welch | O(N log N) | **Frequency analysis** |
| **Spectrogram** | Time-frequency analysis | STFT, PowerSpectrogram | O(M·N log N) | **Chirps, time-varying** |
| **FourierSeries** | Continuous Fourier series | Derivative, Integral, Energy | Function space | **Periodic functions** |
| **FourierBasis** | Orthogonal basis | Evaluate, Normalization | Integration | **Spectral methods** |
| **Windows** | Window functions | Hann, Hamming, Blackman, etc. | O(N) | **Spectral leakage** |
| **FFTNdim** | Multi-dimensional FFT | Forward2D, Inverse2D | O(N² log N) | **Images, 2D signals** |

**Key Insight**: Use FFT/RealFFT for speed (O(N log N)), DFT for reference/testing. Apply windowing to reduce spectral leakage.

## Mathematical Background

### Fourier Transform Pairs

**Continuous Fourier Transform**:
```
F(ω) = ∫_{-∞}^{∞} f(t)·e^{-iωt} dt
f(t) = (1/2π) ∫_{-∞}^{∞} F(ω)·e^{iωt} dω
```

**Discrete Fourier Transform (DFT)**:
```
X[k] = Σ_{n=0}^{N-1} x[n]·e^{-2πikn/N}     (Forward)
x[n] = (1/N) Σ_{k=0}^{N-1} X[k]·e^{2πikn/N} (Inverse)
```

**Key Properties**:
- **Linearity**: DFT(ax + by) = a·DFT(x) + b·DFT(y)
- **Shift**: DFT(x[n-m]) = e^{-2πikm/N}·X[k]
- **Convolution Theorem**: DFT(x * y) = DFT(x) · DFT(y)
- **Parseval's Theorem**: Σ|x[n]|² = (1/N)Σ|X[k]|² (energy conservation)

### Fast Fourier Transform (FFT)

**Cooley-Tukey Algorithm**: Recursively splits DFT into even/odd indices.

**Complexity**: O(N log N) vs O(N²) for naive DFT
- N=1024: FFT needs ~10K ops vs DFT needs ~1M ops (100× faster)
- **Requirement**: N must be power of 2 (use zero-padding if not)

**Algorithm**:
1. **Bit-reversal permutation**: Reorder input
2. **Butterfly operations**: Combine pairs with twiddle factors W_N^k = e^{-2πik/N}
3. **Log₂(N) stages**: Each halves problem size

### Real FFT Optimization

Real signals have **conjugate-symmetric spectrum**: X[k] = X*[N-k]

**Key Insight**: Only need to store N/2+1 complex values instead of N.
- **Memory**: 50% reduction
- **Speed**: ~2× faster than complex FFT
- **Output**: DC (k=0), positive frequencies (k=1..N/2), Nyquist (k=N/2)

### Discrete Cosine Transform (DCT)

**DCT-II** (most common, used in JPEG):
```
X[k] = (2/N) Σ_{n=0}^{N-1} x[n]·cos(π·k·(n+0.5)/N)
```

**Properties**:
- **Real-to-real** transform (no complex arithmetic)
- **Energy compaction**: Most energy in low frequencies (better than DFT for images)
- **Connection to Chebyshev**: DCT-II at Chebyshev nodes computes polynomial coefficients

**Parseval's Theorem for DCT-II**:
```
Σ x[n]² = (N/2)·[0.5·X[0]² + Σ_{k=1}^{N-1} X[k]²]
```

### Windowing

**Problem**: Finite signal length causes **spectral leakage** (energy spreads across frequencies).

**Solution**: Multiply signal by window function that tapers to zero at edges.

**Common Windows**:
- **Rectangular**: w[n] = 1 (no windowing, worst leakage)
- **Hann**: w[n] = 0.5(1 - cos(2πn/N)) (good general purpose)
- **Hamming**: w[n] = 0.54 - 0.46·cos(2πn/N) (slightly better sidelobe)
- **Blackman**: w[n] = 0.42 - 0.5·cos(2πn/N) + 0.08·cos(4πn/N) (best sidelobe, wider main lobe)

**Trade-off**: Better frequency resolution (narrower main lobe) vs better sidelobe suppression.

## Core Transforms

### DFT - Reference Implementation

**Purpose**: Naive O(N²) implementation for testing and non-power-of-2 sizes.

```cpp
#include "mml/algorithms/Fourier/DFT.h"

Vector<Complex> signal(128);
// ... fill signal ...

// Forward DFT
auto spectrum = DFT::Forward(signal);

// Inverse DFT
auto recovered = DFT::Inverse(spectrum);
```

**Use Cases**:
- Testing FFT correctness
- Small signals (N < 32)
- Non-power-of-2 sizes
- Educational/reference

### FFT - Fast Fourier Transform

**Purpose**: Production O(N log N) implementation using Cooley-Tukey algorithm.

```cpp
#include "mml/algorithms/Fourier/FFT.h"

// Method 1: In-place transform
Vector<Complex> data(1024);
// ... fill data ...
FFT::Transform(data, 1);  // 1 = forward, -1 = inverse

// Method 2: Copy-based
auto spectrum = FFT::Forward(data);
auto recovered = FFT::Inverse(spectrum);

// Check if size is valid
if (!FFT::IsPowerOfTwo(n)) {
    int padded = FFT::NextPowerOfTwo(n);
    // Zero-pad to next power of 2
}
```

**Key Features**:
- **Requirement**: N must be power of 2
- **In-place option**: Transform() modifies input (memory efficient)
- **Helpers**: IsPowerOfTwo(), NextPowerOfTwo()
- **Normalization**: Inverse divides by N

### RealFFT - Optimized for Real Signals

**Purpose**: 2× faster FFT for real-valued input exploiting conjugate symmetry.

```cpp
#include "mml/algorithms/Fourier/RealFFT.h"

Vector<Real> signal(1024);
// ... fill signal ...

// Forward: N real → N/2+1 complex
auto spectrum = RealFFT::Forward(signal);  // Length 513

// Inverse: N/2+1 complex → N real
auto recovered = RealFFT::Inverse(spectrum);  // Length 1024
```

**Output Structure**:
- `spectrum[0]` = DC component (k=0)
- `spectrum[1..N/2-1]` = Positive frequencies
- `spectrum[N/2]` = Nyquist frequency

**When to Use**: Always prefer RealFFT over FFT for real signals (audio, sensor data, etc.)

### DCT - Discrete Cosine Transform

**Purpose**: Real-to-real transform used in compression (JPEG) and Chebyshev approximation.

```cpp
#include "mml/algorithms/Fourier/DCT.h"

Vector<Real> signal(64);
// ... fill signal ...

// DCT-II (forward transform)
auto coeffs = DCT::ForwardII(signal);

// DCT-III (inverse transform)
auto recovered = DCT::InverseII(coeffs);

// Verify round-trip
bool ok = DCT::VerifyRoundTrip(signal);  // Should be true

// DST-I for boundary value problems
auto dstCoeffs = DCT::ForwardDST(signal);
```

**DCT Variants**:
- **DCT-II**: Standard "DCT", used in JPEG (8×8 blocks)
- **DCT-III**: Inverse of DCT-II
- **DST-I**: Sine transform for Dirichlet boundary conditions

**Connection to Chebyshev**: Sampling function at Chebyshev nodes and applying DCT-II computes Chebyshev polynomial coefficients.

## Applications

### Convolution - Fast Convolution Theorem

**Purpose**: Compute linear/circular convolution in O(N log N) time.

```cpp
#include "mml/algorithms/Fourier/Convolution.h"

Vector<Real> signal(1000);
Vector<Real> kernel(50);
// ... fill signal and kernel ...

// Linear convolution (length N+M-1)
auto result = Convolution::Linear(signal, kernel);

// Circular convolution (same length as input)
auto circular = Convolution::Circular(signal, kernel);

// Complex convolution
Vector<Complex> x(512), y(512);
auto conv = Convolution::Linear(x, y);
```

**Convolution Theorem**: 
```
x * y = IFFT(FFT(x) · FFT(y))
```

**Applications**:
- **Filtering**: Apply FIR filter kernel
- **Blur**: Convolve image with Gaussian kernel
- **Differentiation**: Convolve with derivative kernel
- **System response**: Input signal * impulse response

### Correlation - Signal Matching

**Purpose**: Measure similarity between signals, detect time delays.

```cpp
#include "mml/algorithms/Fourier/Correlation.h"

Vector<Real> signal1(256);
Vector<Real> signal2(256);
// ... fill signals ...

// Cross-correlation (find time delay)
auto crossCorr = Correlation::Cross(signal1, signal2);

// Auto-correlation (periodicity detection)
auto autoCorr = Correlation::Auto(signal1);

// Normalized correlation (range [-1, 1])
auto normalized = Correlation::CrossNormalized(signal1, signal2);
```

**Output Interpretation**:
- Cross-correlation peak location indicates time delay
- Auto-correlation peak at zero lag = signal energy
- Auto-correlation periodic peaks indicate periodicity

**Applications**:
- **Template matching**: Find pattern in signal
- **Echo detection**: Locate time-delayed copies
- **Pitch detection**: Find fundamental frequency (music/speech)
- **Radar/sonar**: Measure distance via time delay

### PowerSpectrum - Spectral Analysis

**Purpose**: Estimate frequency content and power distribution.

```cpp
#include "mml/algorithms/Fourier/Spectrum.h"
#include "mml/algorithms/Fourier/Windowing.h"

Vector<Real> signal(2048);
// ... fill signal ...

// Simple periodogram
auto spectrum = PowerSpectrum::Compute(signal);

// With custom window
auto window = Windows::Blackman(2048);
auto windowed = PowerSpectrum::Compute(signal, window);

// Welch method (reduced variance)
auto welch = PowerSpectrum::Welch(signal, 1024, 512);  // 1024-pt segments, 50% overlap

// Frequency axis
Real sampleRate = 1000.0;  // Hz
auto freqs = PowerSpectrum::FrequencyAxis(spectrum.size(), sampleRate);
```

**Methods**:
- **Compute**: Simple periodogram (high variance)
- **Welch**: Averaged periodogram (reduced variance, lower resolution)

**Applications**:
- **Peak detection**: Find dominant frequencies
- **Energy distribution**: Where is signal energy concentrated?
- **Noise analysis**: Characterize noise spectrum
- **Quality metrics**: Signal-to-noise ratio (SNR)

### Spectrogram - Time-Frequency Analysis

**Purpose**: Visualize how frequency content changes over time.

```cpp
#include "mml/algorithms/Fourier/Spectrum.h"

Vector<Real> signal(10000);
// ... fill signal (e.g., chirp) ...

int windowSize = 256;
int hopSize = 128;  // 50% overlap

// STFT (Short-Time Fourier Transform)
auto stft = Spectrogram::STFT(signal, windowSize, hopSize);
// Result: Matrix<Complex> of size [numWindows × (windowSize/2+1)]

// Power spectrogram
auto powerSpec = Spectrogram::PowerSpectrogram(signal, windowSize, hopSize);
// Result: Matrix<Real> for visualization
```

**Parameters**:
- **windowSize**: Frequency resolution (larger = better freq resolution, worse time resolution)
- **hopSize**: Time resolution (smaller = better time resolution, more computation)
- **Trade-off**: Time resolution ↔ Frequency resolution (uncertainty principle)

**Applications**:
- **Chirp analysis**: Frequency sweeps over time
- **Speech analysis**: Formant tracking, phoneme segmentation
- **Music analysis**: Note onset detection, tempo tracking
- **Non-stationary signals**: Any signal with time-varying frequency content

### Windows - Spectral Leakage Reduction

**Purpose**: Taper signals to reduce frequency leakage in FFT.

```cpp
#include "mml/algorithms/Fourier/Windowing.h"

int N = 512;

// Common windows
auto rectangular = Windows::Rectangular(N);  // No windowing
auto hann = Windows::Hann(N);                // General purpose
auto hamming = Windows::Hamming(N);          // Slightly better sidelobes
auto blackman = Windows::Blackman(N);        // Best sidelobe suppression
auto bartlett = Windows::Bartlett(N);        // Triangular
auto welch = Windows::Welch(N);              // Parabolic
auto kaiser = Windows::Kaiser(N, 8.6);       // Adjustable (beta parameter)
auto gaussian = Windows::Gaussian(N, 0.4);   // Smooth Gaussian

// Apply window
Vector<Real> signal(N);
for (int i = 0; i < N; i++)
    signal[i] *= hann[i];
```

**Selection Guide**:
- **Hann**: Good default, excellent general-purpose window
- **Hamming**: Better stopband attenuation than Hann
- **Blackman**: Best sidelobe suppression, wider main lobe
- **Kaiser**: Adjustable trade-off via beta parameter
- **Rectangular**: No windowing (use only for complete periods)

## Advanced Features

### FourierSeries - Continuous Fourier Series

**Purpose**: Represent periodic functions as sum of sines/cosines.

```cpp
#include "mml/algorithms/Fourier/FourierSeries.h"

// Define periodic function
auto func = [](Real x) { return x * x; };  // f(x) = x² on [-π, π]
Real L = M_PI;  // Half-period
int n = 20;     // Number of terms

// Create Fourier series approximation
FourierSeries fs(func, L, n);

// Evaluate
Real value = fs(0.5);

// Calculus operations
FourierSeries derivative = fs.Derivative();
FourierSeries integral = fs.Integral();

// Energy (Parseval's theorem)
Real energy = fs.Energy();

// Get complex coefficients
auto complexCoeffs = fs.ComplexCoefficients();
```

**Applications**:
- **PDE solving**: Separation of variables (heat equation, wave equation)
- **Function approximation**: Represent periodic functions
- **Signal synthesis**: Build complex waveforms from harmonics

### FourierBasis - Orthogonal Basis Integration

**Purpose**: Fourier basis for spectral expansion framework.

```cpp
#include "mml/algorithms/Fourier/FourierBasis.h"

// Real Fourier basis (cos/sin)
FourierBasis basis(M_PI);  // Period 2π

// Evaluate basis function n at x
Real value = basis.Evaluate(n, x);

// Complex exponential basis
ComplexFourierBasis complexBasis(M_PI);
Complex cValue = complexBasis.Evaluate(n, x);
```

**Integration**: Implements `IOrthogonalBasis` interface for use with spectral expansion framework.

### FFTNdim - Multi-Dimensional FFT

**Purpose**: 2D FFT for images and 2D signals.

```cpp
#include "mml/algorithms/Fourier/FourierNdim.h"

Matrix<Complex> image(512, 512);
// ... fill image ...

// 2D FFT (separable: FFT rows, then columns)
auto spectrum = FFTNdim::Forward2D(image);

// 2D Inverse FFT
auto recovered = FFTNdim::Inverse2D(spectrum);

// Real 2D FFT (for real-valued images)
Matrix<Real> realImage(512, 512);
auto realSpectrum = FFTNdim::ForwardReal2D(realImage);
```

**Applications**:
- **Image processing**: Filtering, compression, feature extraction
- **Pattern matching**: 2D correlation
- **Deconvolution**: Image restoration
- **Texture analysis**: Frequency content of textures

## Usage Examples

### Example 1: Signal Filtering (Low-Pass)

```cpp
#include "mml/algorithms/Fourier/RealFFT.h"

// Create noisy signal
Vector<Real> signal(1024);
for (int i = 0; i < 1024; i++) {
    signal[i] = sin(2*M_PI*10*i/1024) +      // 10 Hz signal
                0.5*sin(2*M_PI*100*i/1024);  // 100 Hz noise
}

// Transform to frequency domain
auto spectrum = RealFFT::Forward(signal);

// Low-pass filter: Keep frequencies below 50 Hz
Real sampleRate = 1024.0;  // Hz
Real cutoff = 50.0;        // Hz
int cutoffBin = (int)(cutoff * spectrum.size() / (sampleRate/2));

for (int k = cutoffBin; k < spectrum.size(); k++) {
    spectrum[k] = Complex(0, 0);  // Zero out high frequencies
}

// Transform back to time domain
auto filtered = RealFFT::Inverse(spectrum);
// Result: Clean 10 Hz signal without 100 Hz noise
```

### Example 2: Spectral Analysis

```cpp
#include "mml/algorithms/Fourier/Spectrum.h"
#include "mml/algorithms/Fourier/Windowing.h"

Vector<Real> signal(2048);
// ... load signal from sensor ...

// Apply Hann window to reduce leakage
auto window = Windows::Hann(2048);
auto spectrum = PowerSpectrum::Compute(signal, window);

// Find peak frequency
Real maxPower = 0.0;
int peakBin = 0;
for (int k = 0; k < spectrum.size(); k++) {
    if (spectrum[k] > maxPower) {
        maxPower = spectrum[k];
        peakBin = k;
    }
}

Real sampleRate = 1000.0;  // Hz
Real peakFreq = peakBin * (sampleRate / 2) / (spectrum.size() - 1);
std::cout << "Dominant frequency: " << peakFreq << " Hz\n";
```

### Example 3: Convolution (System Response)

```cpp
#include "mml/algorithms/Fourier/Convolution.h"

// Input signal: unit impulse train
Vector<Real> input(1000);
for (int i = 0; i < 1000; i += 100) {
    input[i] = 1.0;  // Impulse every 100 samples
}

// Impulse response: exponential decay
Vector<Real> impulseResponse(50);
for (int i = 0; i < 50; i++) {
    impulseResponse[i] = exp(-0.1 * i);
}

// System output = input convolved with impulse response
auto output = Convolution::Linear(input, impulseResponse);

// Result: Decaying pulses at each impulse location
```

### Example 4: Time-Delay Detection

```cpp
#include "mml/algorithms/Fourier/Correlation.h"

// Original signal
Vector<Real> signal1(512);
for (int i = 0; i < 512; i++) {
    signal1[i] = sin(2*M_PI*5*i/512);
}

// Delayed copy (delay = 20 samples)
Vector<Real> signal2(512, 0.0);
for (int i = 20; i < 512; i++) {
    signal2[i] = signal1[i - 20];
}

// Cross-correlation
auto corr = Correlation::Cross(signal1, signal2);

// Find peak
Real maxCorr = 0.0;
int peakLag = 0;
int center = signal1.size() - 1;
for (int i = 0; i < corr.size(); i++) {
    if (abs(corr[i]) > abs(maxCorr)) {
        maxCorr = corr[i];
        peakLag = i - center;
    }
}

std::cout << "Detected delay: " << abs(peakLag) << " samples\n";
// Output: Detected delay: 20 samples
```

## Performance Notes

### Complexity Summary

| Operation | Complexity | Notes |
|-----------|-----------|-------|
| DFT | O(N²) | Reference only |
| FFT | O(N log N) | N must be power of 2 |
| RealFFT | O(N log N) | ~2× faster than FFT |
| DCT/DST | O(N²) naive | Fast DCT via FFT possible (future) |
| Convolution | O(N log N) | Via FFT, much faster than O(N²) direct |
| Correlation | O(N log N) | Via FFT |
| Power Spectrum | O(N log N) | FFT-based |
| Spectrogram | O(M·N log N) | M windows of size N |

### Optimization Tips

1. **Power-of-2 sizes**: Always use power-of-2 for FFT
   - Zero-pad if necessary: `FFT::NextPowerOfTwo(n)`
   
2. **Real signals**: Use RealFFT instead of FFT
   - 2× faster, 50% less memory

3. **In-place transforms**: Use `FFT::Transform()` for memory efficiency
   - Modifies input but avoids allocation

4. **Welch method**: Reduce spectral variance
   - Trade-off: Reduced frequency resolution

5. **Window selection**: Balance leakage vs resolution
   - Hann/Hamming: Good general purpose
   - Blackman: Best sidelobe suppression

### Memory Usage

| Transform | Input Size | Output Size | Memory |
|-----------|-----------|-------------|--------|
| FFT | N complex | N complex | O(N) |
| RealFFT | N real | N/2+1 complex | O(N) |
| DCT-II | N real | N real | O(N) |
| Convolution | N + M | N + M - 1 | O(N + M) |
| Spectrogram | N real | [W × N/2] complex | O(W·N) |

## Integration with MML

### Connection to ChebyshevApproximation

Chebyshev polynomial approximation uses DCT-II internally:

```cpp
#include "mml/algorithms/ChebyshevApproximation.h"

// Sample function at Chebyshev nodes
auto func = [](Real x) { return sin(x); };
ChebyshevApproximation approx(func, -1, 1, 50);

// Internally uses DCT::ForwardII() to compute coefficients
// This is the discrete orthogonality relation for Chebyshev polynomials
```

**Why it works**: At Chebyshev nodes x_k = cos(π(k+0.5)/N), the continuous orthogonality integral for Chebyshev polynomials becomes exactly the DCT-II formula.

### Function Spaces

FourierBasis implements IOrthogonalBasis:

```cpp
#include "mml/core/OrthogonalBasis.h"
#include "mml/algorithms/Fourier/FourierBasis.h"

// Can be used with spectral expansion framework
FourierBasis basis(M_PI);
// SpectralExpansion<FourierBasis> expansion(...);  // Future feature
```

## References

### Books
- **Numerical Recipes, 3rd Edition** - Chapters on FFT, convolution, correlation, spectral analysis
- **Digital Signal Processing** by Oppenheim & Schafer - Comprehensive DSP text
- **The Scientist and Engineer's Guide to Digital Signal Processing** by Steven W. Smith - Free online

### Papers
- Cooley, J. W., & Tukey, J. W. (1965). "An algorithm for the machine calculation of complex Fourier series"
- Welch, P. (1967). "The use of fast Fourier transform for the estimation of power spectra"

### Online Resources
- [FFTW](http://www.fftw.org/) - Fastest FFT in the West (reference implementation)
- [Wikipedia: FFT](https://en.wikipedia.org/wiki/Fast_Fourier_transform)
- [Wikipedia: Window function](https://en.wikipedia.org/wiki/Window_function)
- [DSP Guide Online](http://www.dspguide.com/) - Free comprehensive guide

## Future Enhancements

1. **Fast DCT via FFT**: O(N log N) DCT implementation
2. **2D DCT**: For image compression (JPEG)
3. **Chirp Z-transform**: Non-uniform frequency sampling
4. **Arbitrary-length FFT**: Bluestein algorithm for non-power-of-2
5. **Fast convolution**: Overlap-add, overlap-save methods for streaming
6. **Wavelet transforms**: Multi-resolution analysis
7. **Hilbert transform**: Analytic signal, instantaneous frequency

## See Also

- [Function_analyzer.md](Function_analyzer.md) - Function analysis and feature extraction
- [Statistics.md](Statistics.md) - Statistical analysis of signals
- [README_Algorithms.md](README_Algorithms.md) - Algorithm overview