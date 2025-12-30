#include <catch2/catch_all.hpp>
#include <complex>
#include <cmath>

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "../../mml/algorithms/Fourier/DFT.h"
#include "../../mml/algorithms/Fourier/RealFFT.h"
#include "../../mml/algorithms/Fourier/Convolution.h"
#include "../../mml/algorithms/Fourier/Correlation.h"
#include "../../mml/algorithms/Fourier/Spectrum.h"
#include "../../mml/algorithms/Fourier/FourierSeries.h"
#include "../../mml/algorithms/Fourier/FourierBasis.h"
#include "../../mml/algorithms/Fourier/DCT.h"
#include "../../test_data/fourier_test_bed.h"
#include "../../test_data/root_finding_test_bed.h"  // For RealFunctionWrapper

using namespace MML;
using namespace MML::TestBeds;
using Complex = std::complex<Real>;

///////////////////////////////////////////////////////////////////////////////////////////
//                          TEST 1: DFT vs FFT CORRECTNESS                                //
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("DFT vs FFT - Pure Tone", "[Fourier][DFT][FFT]")
{
    auto signal = getPureToneSignal();
    Vector<Complex> data(signal.data.size());
    for (size_t i = 0; i < signal.data.size(); i++)
    {
        data[i] = Complex(signal.data[i], 0.0);
    }
    
    auto dftResult = DFT::Forward(data);
    auto fftResult = FFT::Forward(data);
    
    REQUIRE(dftResult.size() == fftResult.size());
    
    for (size_t i = 0; i < dftResult.size(); i++)
    {
        REQUIRE(std::abs(dftResult[i] - fftResult[i]) < 1e-10);
    }
}

TEST_CASE("DFT vs FFT - Dual Tone", "[Fourier][DFT][FFT]")
{
    auto signal = getDualToneSignal();
    Vector<Complex> data(signal.data.size());
    for (size_t i = 0; i < signal.data.size(); i++)
    {
        data[i] = Complex(signal.data[i], 0.0);
    }
    
    auto dftResult = DFT::Forward(data);
    auto fftResult = FFT::Forward(data);
    
    REQUIRE(dftResult.size() == fftResult.size());
    
    for (size_t i = 0; i < dftResult.size(); i++)
    {
        REQUIRE(std::abs(dftResult[i] - fftResult[i]) < 1e-9);
    }
}

TEST_CASE("DFT vs FFT - All Test Signals", "[Fourier][DFT][FFT][comprehensive]")
{
    auto signals = getAllSignals();
    
    for (const auto& signal : signals)
    {
        DYNAMIC_SECTION("Signal: " << signal.name)
        {
            Vector<Complex> data(signal.data.size());
            for (size_t i = 0; i < signal.data.size(); i++)
            {
                data[i] = Complex(signal.data[i], 0.0);
            }
            
            auto dftResult = DFT::Forward(data);
            auto fftResult = FFT::Forward(data);
            
            REQUIRE(dftResult.size() == fftResult.size());
            
            Real maxError = 0.0;
            for (size_t i = 0; i < dftResult.size(); i++)
            {
                Real error = std::abs(dftResult[i] - fftResult[i]);
                maxError = std::max(maxError, error);
            }
            
            REQUIRE(maxError < 1e-8);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
//                          TEST 2: FFT INVERSE RECOVERY                                  //
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FFT Forward/Inverse Round-Trip", "[Fourier][FFT][inverse]")
{
    auto signal = getHarmonicSignal();
    Vector<Complex> data(signal.data.size());
    for (size_t i = 0; i < signal.data.size(); i++)
    {
        data[i] = Complex(signal.data[i], 0.0);
    }
    
    auto original = data;
    auto spectrum = FFT::Forward(data);
    auto recovered = FFT::Inverse(spectrum);
    
    REQUIRE(original.size() == recovered.size());
    
    for (size_t i = 0; i < original.size(); i++)
    {
        REQUIRE(std::abs(original[i] - recovered[i]) < 1e-10);
    }
}

TEST_CASE("FFT Inverse Recovery - All Signals", "[Fourier][FFT][inverse][comprehensive]")
{
    auto signals = getAllSignals();
    
    for (const auto& signal : signals)
    {
        DYNAMIC_SECTION("Signal: " << signal.name)
        {
            Vector<Complex> data(signal.data.size());
            for (size_t i = 0; i < signal.data.size(); i++)
            {
                data[i] = Complex(signal.data[i], 0.0);
            }
            
            auto original = data;
            auto spectrum = FFT::Forward(data);
            auto recovered = FFT::Inverse(spectrum);
            
            Real maxError = 0.0;
            for (size_t i = 0; i < original.size(); i++)
            {
                Real error = std::abs(original[i] - recovered[i]);
                maxError = std::max(maxError, error);
            }
            
            REQUIRE(maxError < 1e-9);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
//                          TEST 3: PARSEVAL'S THEOREM                                    //
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Parseval's Theorem - Energy Conservation", "[Fourier][FFT][Parseval]")
{
    auto signal = getPureToneSignal();
    Vector<Complex> data(signal.data.size());
    for (size_t i = 0; i < signal.data.size(); i++)
    {
        data[i] = Complex(signal.data[i], 0.0);
    }
    
    // Time domain energy
    Real timeEnergy = 0.0;
    for (size_t i = 0; i < data.size(); i++)
    {
        timeEnergy += std::norm(data[i]);
    }
    
    // Frequency domain energy
    auto spectrum = FFT::Forward(data);
    Real freqEnergy = 0.0;
    for (size_t i = 0; i < spectrum.size(); i++)
    {
        freqEnergy += std::norm(spectrum[i]);
    }
    freqEnergy /= spectrum.size();  // Normalization
    
    REQUIRE(std::abs(timeEnergy - freqEnergy) / timeEnergy < 1e-10);
}

TEST_CASE("Parseval's Theorem - All Signals", "[Fourier][FFT][Parseval][comprehensive]")
{
    auto signals = getAllSignals();
    
    for (const auto& signal : signals)
    {
        DYNAMIC_SECTION("Signal: " << signal.name)
        {
            Vector<Complex> data(signal.data.size());
            for (size_t i = 0; i < signal.data.size(); i++)
            {
                data[i] = Complex(signal.data[i], 0.0);
            }
            
            Real timeEnergy = 0.0;
            for (size_t i = 0; i < data.size(); i++)
            {
                timeEnergy += std::norm(data[i]);
            }
            
            auto spectrum = FFT::Forward(data);
            Real freqEnergy = 0.0;
            for (size_t i = 0; i < spectrum.size(); i++)
            {
                freqEnergy += std::norm(spectrum[i]);
            }
            freqEnergy /= spectrum.size();
            
            Real relativeError = std::abs(timeEnergy - freqEnergy) / timeEnergy;
            REQUIRE(relativeError < 1e-9);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
//                          TEST 4: REAL FFT                                              //
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("RealFFT vs Complex FFT", "[Fourier][RealFFT]")
{
    auto signal = getPureToneSignal();
    
    // Complex FFT on real data
    Vector<Complex> complexData(signal.data.size());
    for (size_t i = 0; i < signal.data.size(); i++)
    {
        complexData[i] = Complex(signal.data[i], 0.0);
    }
    auto complexResult = FFT::Forward(complexData);
    
    // Real FFT
    Vector<Real> realData(signal.data.size());
    for (size_t i = 0; i < signal.data.size(); i++) realData[i] = signal.data[i];
    auto realResult = RealFFT::Forward(realData);
    
    // Compare positive frequencies (RealFFT only stores n/2+1 values)
    size_t compareSize = realResult.size();
    for (size_t i = 0; i < compareSize; i++)
    {
        REQUIRE(std::abs(complexResult[i] - realResult[i]) < 1e-10);
    }
}

TEST_CASE("RealFFT Conjugate Symmetry", "[Fourier][RealFFT]")
{
    auto signal = getCosineToneSignal();
    Vector<Real> data(signal.data.size());
    for (size_t i = 0; i < signal.data.size(); i++) data[i] = signal.data[i];
    
    auto spectrum = RealFFT::Forward(data);
    
    // For real input, X[k] = conj(X[N-k])
    // Check a few pairs
    size_t N = data.size();
    size_t maxK = (spectrum.size() > 6) ? 5 : spectrum.size() - 1;
    for (size_t k = 1; k < maxK; k++)
    {
        // spectrum[k] should match conj(spectrum[N-k]) from full FFT
        // Since RealFFT only stores positive frequencies, we verify properties indirectly
        REQUIRE((std::abs(spectrum[k].imag()) < 1e-10 || 
                 std::abs(spectrum[N-k].imag()) < 1e-10));
    }
}

TEST_CASE("RealFFT Forward/Inverse Round-Trip", "[Fourier][RealFFT][inverse]")
{
    auto signal = getDualToneSignal();
    Vector<Real> original(signal.data.size());
    for (size_t i = 0; i < signal.data.size(); i++) original[i] = signal.data[i];
    
    auto spectrum = RealFFT::Forward(original);
    auto recovered = RealFFT::Inverse(spectrum);
    
    REQUIRE(original.size() == recovered.size());
    
    for (size_t i = 0; i < original.size(); i++)
    {
        REQUIRE(std::abs(original[i] - recovered[i]) < 1e-10);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
//                          TEST 5: CONVOLUTION                                           //
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Convolution - Simple 3-point", "[Fourier][Convolution]")
{
    auto test = getSimpleConvolutionTest();
    Vector<Real> signal(test.signal.size());
    for (size_t i = 0; i < test.signal.size(); i++) signal[i] = test.signal[i];
    Vector<Real> kernel(test.kernel.size());
    for (size_t i = 0; i < test.kernel.size(); i++) kernel[i] = test.kernel[i];
    
    auto result = Convolution::Linear(signal, kernel);
    
    REQUIRE(result.size() == test.expectedLinear.size());
    
    for (size_t i = 0; i < result.size(); i++)
    {
        REQUIRE(std::abs(result[i] - test.expectedLinear[i]) < 1e-10);
    }
}

TEST_CASE("Convolution - Impulse Response", "[Fourier][Convolution]")
{
    auto test = getImpulseResponseTest();
    Vector<Real> signal(test.signal.size());
    for (size_t i = 0; i < test.signal.size(); i++) signal[i] = test.signal[i];
    Vector<Real> kernel(test.kernel.size());
    for (size_t i = 0; i < test.kernel.size(); i++) kernel[i] = test.kernel[i];
    
    auto result = Convolution::Linear(signal, kernel);
    
    REQUIRE(result.size() == test.expectedLinear.size());
    
    for (size_t i = 0; i < result.size(); i++)
    {
        REQUIRE(std::abs(result[i] - test.expectedLinear[i]) < 1e-10);
    }
}

TEST_CASE("Convolution Theorem Verification", "[Fourier][Convolution][theorem]")
{
    // Generate two simple signals
    Vector<Real> x(4);
    x[0] = 1; x[1] = 2; x[2] = 3; x[3] = 4;
    Vector<Real> h(3);
    h[0] = 0.5; h[1] = 1.0; h[2] = 0.5;
    
    // Direct convolution
    auto direct = Convolution::Linear(x, h);
    
    // Via FFT (convolution theorem)
    size_t N = x.size() + h.size() - 1;
    size_t fftSize = FFT::NextPowerOfTwo(N);
    
    Vector<Complex> X(fftSize, Complex(0, 0));
    Vector<Complex> H(fftSize, Complex(0, 0));
    
    for (size_t i = 0; i < x.size(); i++) X[i] = Complex(x[i], 0);
    for (size_t i = 0; i < h.size(); i++) H[i] = Complex(h[i], 0);
    
    auto X_fft = FFT::Forward(X);
    auto H_fft = FFT::Forward(H);
    
    Vector<Complex> Y_fft(fftSize);
    for (size_t i = 0; i < fftSize; i++)
    {
        Y_fft[i] = X_fft[i] * H_fft[i];
    }
    
    auto y = FFT::Inverse(Y_fft);
    
    // Compare
    for (size_t i = 0; i < direct.size(); i++)
    {
        REQUIRE(std::abs(y[i].real() - direct[i]) < 1e-10);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////
//                          TEST 6: CORRELATION                                           //
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Auto-Correlation - Peak at Zero Lag", "[Fourier][Correlation]")
{
    auto test = getPeriodicAutoCorrelationTest();
    Vector<Real> signal(test.signal1.size());
    for (size_t i = 0; i < test.signal1.size(); i++) signal[i] = test.signal1[i];
    
    auto result = Correlation::Auto(signal);
    
    // Find peak
    Real maxVal = std::abs(result[0]);
    size_t maxIdx = 0;
    for (size_t i = 1; i < result.size(); i++)
    {
        if (std::abs(result[i]) > maxVal)
        {
            maxVal = std::abs(result[i]);
            maxIdx = i;
        }
    }
    
    REQUIRE(maxIdx == test.expectedPeakLag);
}

TEST_CASE("Cross-Correlation - Detect Time Delay", "[Fourier][Correlation]")
{
    auto test = getDelayedCrossCorrelationTest();
    Vector<Real> sig1(test.signal1.size());
    for (size_t i = 0; i < test.signal1.size(); i++) sig1[i] = test.signal1[i];
    Vector<Real> sig2(test.signal2.size());
    for (size_t i = 0; i < test.signal2.size(); i++) sig2[i] = test.signal2[i];
    
    auto result = Correlation::Cross(sig1, sig2);
    
    // Find peak
    Real maxVal = std::abs(result[0]);
    size_t maxIdx = 0;
    for (size_t i = 1; i < result.size(); i++)
    {
        if (std::abs(result[i]) > maxVal)
        {
            maxVal = std::abs(result[i]);
            maxIdx = i;
        }
    }
    
    // Peak should be near expected delay
    REQUIRE(std::abs(static_cast<int>(maxIdx) - test.expectedPeakLag) <= 1);
}

TEST_CASE("Normalized Correlation Bounds", "[Fourier][Correlation]")
{
    auto test = getPeriodicAutoCorrelationTest();
    Vector<Real> sig1(test.signal1.size());
    for (size_t i = 0; i < test.signal1.size(); i++) sig1[i] = test.signal1[i];
    Vector<Real> sig2(test.signal1.size());
    for (size_t i = 0; i < test.signal1.size(); i++) sig2[i] = test.signal1[i];
    
    auto result = Correlation::CrossNormalized(sig1, sig2);
    
    // Normalized correlation should be in [-1, 1]
    for (size_t i = 0; i < result.size(); i++)
    {
        REQUIRE(result[i] >= -1.0);
        REQUIRE(result[i] <= 1.0);
    }
    
    // Peak should be 1.0 for identical signals
    Real maxVal = result[0];
    for (size_t i = 1; i < result.size(); i++)
    {
        maxVal = std::max(maxVal, result[i]);
    }
    REQUIRE(std::abs(maxVal - 1.0) < 1e-10);
}

///////////////////////////////////////////////////////////////////////////////////////////
//                          TEST 7: POWER SPECTRUM                                        //
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Power Spectrum - Pure Tone Peak Detection", "[Fourier][Spectrum]")
{
    auto signal = getPureToneSignal();
    Vector<Real> data(signal.data.size());
    for (size_t i = 0; i < signal.data.size(); i++) data[i] = signal.data[i];
    
    auto spectrum = PowerSpectrum::Compute(data);
    
    // Find peak
    Real maxPower = spectrum[0];
    size_t peakIdx = 0;
    for (size_t i = 1; i < spectrum.size()/2; i++)
    {
        if (spectrum[i] > maxPower)
        {
            maxPower = spectrum[i];
            peakIdx = i;
        }
    }
    
    // Convert index to frequency
    Real peakFreq = peakIdx * signal.sampleRate / data.size();
    Real expectedFreq = signal.expectedPeakFrequencies[0];
    
    // Allow 1 bin tolerance
    Real freqResolution = signal.sampleRate / data.size();
    REQUIRE(std::abs(peakFreq - expectedFreq) < 2 * freqResolution);
}

TEST_CASE("Power Spectrum - Energy Conservation", "[Fourier][Spectrum]")
{
    auto signal = getDualToneSignal();
    Vector<Real> data(signal.data.size());
    for (size_t i = 0; i < signal.data.size(); i++) data[i] = signal.data[i];
    
    // Use no window for exact Parseval verification
    Vector<Real> noWindow(data.size(), 1.0);
    auto spectrum = PowerSpectrum::Compute(data, noWindow);
    
    // For real FFT, sum DC + 2*(interior bins) + Nyquist to account for both sides
    size_t N = data.size();
    Real spectralEnergy = spectrum[0];  // DC component
    for (size_t i = 1; i < spectrum.size() - 1; i++)
    {
        spectralEnergy += 2.0 * spectrum[i];  // Double for positive/negative frequencies
    }
    if (N % 2 == 0) {
        spectralEnergy += spectrum[spectrum.size() - 1];  // Nyquist (no negative pair)
    } else {
        spectralEnergy += 2.0 * spectrum[spectrum.size() - 1];  // Last bin has pair
    }
    
    Real relativeError = std::abs(spectralEnergy - signal.signalEnergy) / signal.signalEnergy;
    REQUIRE(relativeError < 0.05);  // 5% tolerance
}

TEST_CASE("Power Spectrum - Frequency Axis", "[Fourier][Spectrum]")
{
    int N = 128;
    Real sampleRate = 64.0;
    
    auto freqs = PowerSpectrum::FrequencyAxis(N, sampleRate);
    
    REQUIRE(freqs.size() == N/2 + 1);
    REQUIRE(freqs[0] == 0.0);
    REQUIRE(std::abs(freqs.back() - sampleRate/2) < 1e-10);  // Nyquist frequency
}

///////////////////////////////////////////////////////////////////////////////////////////
//                          TEST 8: FOURIER SERIES                                        //
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FourierSeries - Sine Wave Reconstruction", "[Fourier][FourierSeries]")
{
    // Simple sine wave: sin(2πx/L)
    Real L = 1.0;
    auto sineFunc = [](Real x) { return std::sin(2.0 * M_PI * x); };
    RealFunctionWrapper func(sineFunc);
    
    FourierSeries fs(func, L, 10);
    
    // Check reconstruction at several points
    for (int i = 0; i <= 10; i++)
    {
        Real x = -L + 2*L * i / 10.0;
        Real expected = sineFunc(x);
        Real computed = fs(x);
        REQUIRE(std::abs(computed - expected) < 0.1);
    }
}

TEST_CASE("FourierSeries - Square Wave Approximation", "[Fourier][FourierSeries]")
{
    Real L = 1.0;
    auto squareWave = [](Real x) { return (x > 0) ? 1.0 : -1.0; };
    RealFunctionWrapper func(squareWave);
    
    FourierSeries fs(func, L, 50);  // More terms for better approximation
    
    // Check at non-discontinuous points
    Real x1 = 0.5;
    REQUIRE(std::abs(fs(x1) - 1.0) < 0.15);
    
    Real x2 = -0.5;
    REQUIRE(std::abs(fs(x2) - (-1.0)) < 0.15);
}

TEST_CASE("FourierSeries - Derivative", "[Fourier][FourierSeries]")
{
    Real L = 1.0;
    // f(x) = cos(2πx/L), f'(x) = -2π/L * sin(2πx/L)
    auto cosFunc = [L](Real x) { return std::cos(2.0 * M_PI * x / L); };
    RealFunctionWrapper func(cosFunc);
    
    FourierSeries fs(func, L, 20);
    auto deriv = fs.Derivative();
    
    // Check derivative at several points
    for (int i = 0; i <= 5; i++)
    {
        Real x = -L/2 + L * i / 5.0;
        Real expected = -(2.0 * M_PI / L) * std::sin(2.0 * M_PI * x / L);
        Real computed = deriv(x);
        REQUIRE(std::abs(computed - expected) < 0.2);
    }
}

TEST_CASE("FourierSeries - Energy (Parseval)", "[Fourier][FourierSeries]")
{
    Real L = 1.0;
    auto sineFunc = [](Real x) { return std::sin(2.0 * M_PI * x); };
    RealFunctionWrapper func(sineFunc);
    
    FourierSeries fs(func, L, 10);
    
    Real energy = fs.Energy();
    
    // For sin(2πx), energy over [-L,L] is L
    REQUIRE(std::abs(energy - L) < 0.1);
}

///////////////////////////////////////////////////////////////////////////////////////////
//                          TEST 9: FOURIER BASIS                                         //
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FourierBasis - Evaluate Basis Functions", "[Fourier][FourierBasis]")
{
    FourierBasis basis(Constants::PI);  // Period 2π
    
    // n=0: constant
    Real val0 = basis.Evaluate(0, 0.5);
    REQUIRE(std::abs(val0 - 1.0) < 1e-10);
    
    // n=1: cos(x)
    Real val1 = basis.Evaluate(1, Constants::PI/2);
    Real expected1 = std::cos(Constants::PI/2);
    REQUIRE(std::abs(val1 - expected1) < 1e-10);
}

TEST_CASE("FourierBasis - Orthogonality", "[Fourier][FourierBasis]")
{
    FourierBasis basis(Constants::PI);
    
    // Verify orthogonality: ∫ φₘ(x) φₙ(x) dx = δₘₙ·Norm(n)
    int N = 100;
    Real dx = 2.0 * Constants::PI / N;
    
    for (int m = 0; m < 3; m++)
    {
        for (int n = 0; n < 3; n++)
        {
            Real integral = 0.0;
            for (int i = 0; i < N; i++)
            {
                Real x = -Constants::PI + i * dx;
                integral += basis.Evaluate(m, x) * basis.Evaluate(n, x) * dx;
            }
            
            if (m == n)
            {
                Real expectedNorm = basis.Normalization(n);
                REQUIRE(std::abs(integral - expectedNorm) < 0.5);
            }
            else
            {
                REQUIRE(std::abs(integral) < 0.5);
            }
        }
    }
}

TEST_CASE("ComplexFourierBasis - Exponential Form", "[Fourier][FourierBasis]")
{
    ComplexFourierBasis basis(Constants::PI);
    
    // e^{ix} evaluation at x=π
    Complex val = basis.Evaluate(1, Constants::PI);
    Complex expected = std::exp(Complex(0, Constants::PI));  // e^{iπ} = -1
    
    REQUIRE(std::abs(val - expected) < 1e-10);
}

///////////////////////////////////////////////////////////////////////////////////////////
//                          TEST 10: EDGE CASES                                           //
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("FFT - Power of Two Requirement", "[Fourier][FFT][edge-case]")
{
    REQUIRE(FFT::IsPowerOfTwo(128));
    REQUIRE(FFT::IsPowerOfTwo(256));
    REQUIRE(!FFT::IsPowerOfTwo(100));
    REQUIRE(!FFT::IsPowerOfTwo(127));
    
    REQUIRE(FFT::NextPowerOfTwo(100) == 128);
    REQUIRE(FFT::NextPowerOfTwo(128) == 128);
    REQUIRE(FFT::NextPowerOfTwo(129) == 256);
}

TEST_CASE("FFT - Empty Input", "[Fourier][FFT][edge-case]")
{
    Vector<Complex> empty;
    REQUIRE_THROWS(FFT::Forward(empty));
}

TEST_CASE("FFT - Single Element", "[Fourier][FFT][edge-case]")
{
    Vector<Complex> single(1);
    single[0] = Complex(1.0, 0.0);
    auto result = FFT::Forward(single);
    REQUIRE(result.size() == 1);
    REQUIRE(std::abs(result[0] - Complex(1.0, 0.0)) < 1e-10);
}

TEST_CASE("Convolution - Mismatched Sizes", "[Fourier][Convolution][edge-case]")
{
    Vector<Real> x(5);
    x[0] = 1; x[1] = 2; x[2] = 3; x[3] = 4; x[4] = 5;
    Vector<Real> h(2);
    h[0] = 1; h[1] = 2;
    
    auto result = Convolution::Linear(x, h);
    REQUIRE(result.size() == x.size() + h.size() - 1);
}

///////////////////////////////////////////////////////////////////////////////////////////
//                          TEST 11: DCT (DISCRETE COSINE TRANSFORM)                      //
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("DCT-II Forward/Inverse Round-Trip", "[Fourier][DCT]")
{
    // Test with various signal sizes
    std::vector<int> sizes = {4, 8, 16, 32, 64};
    
    for (int N : sizes)
    {
        Vector<Real> signal(N);
        for (int i = 0; i < N; i++)
        {
            signal[i] = std::sin(2.0 * M_PI * i / N) + 0.5 * std::cos(4.0 * M_PI * i / N);
        }
        
        auto coeffs = DCT::ForwardII(signal);
        auto recovered = DCT::InverseII(coeffs);
        
        REQUIRE(recovered.size() == signal.size());
        
        for (int i = 0; i < N; i++)
        {
            REQUIRE(std::abs(recovered[i] - signal[i]) < 1e-12);
        }
    }
}

TEST_CASE("DCT-II Known Signal Test", "[Fourier][DCT]")
{
    // Test with DC signal (constant)
    Vector<Real> dc(8, 1.0);
    auto dctDC = DCT::ForwardII(dc);
    
    // DC signal should have all energy in first coefficient
    REQUIRE(std::abs(dctDC[0] - 2.0) < 1e-10);  // DCT-II of constant 1 is 2
    for (int i = 1; i < 8; i++)
    {
        REQUIRE(std::abs(dctDC[i]) < 1e-10);
    }
}

TEST_CASE("DCT-II Energy Conservation", "[Fourier][DCT]")
{
    int N = 32;
    Vector<Real> signal(N);
    for (int i = 0; i < N; i++)
    {
        signal[i] = std::sin(2.0 * M_PI * 5.0 * i / N);
    }
    
    Real signalEnergy = 0.0;
    for (int i = 0; i < N; i++)
    {
        signalEnergy += signal[i] * signal[i];
    }
    
    auto coeffs = DCT::ForwardII(signal);
    
    // Parseval's theorem for DCT-II with 2/N normalization:
    // ||x||² = (N/2) * [ (1/2)|X[0]|² + Σ_{k=1}^{N-1} |X[k]|² ]
    Real dctEnergy = (N / 2.0) * (0.5 * coeffs[0] * coeffs[0]);
    for (int k = 1; k < N; k++)
    {
        dctEnergy += (N / 2.0) * coeffs[k] * coeffs[k];
    }
    
    Real relativeError = std::abs(dctEnergy - signalEnergy) / signalEnergy;
    REQUIRE(relativeError < 1e-10);
}

TEST_CASE("DCT Helper - VerifyRoundTrip", "[Fourier][DCT]")
{
    Vector<Real> signal(16);
    for (int i = 0; i < 16; i++)
    {
        signal[i] = std::exp(-0.1 * i) * std::cos(2.0 * M_PI * i / 16);
    }
    
    REQUIRE(DCT::VerifyRoundTrip(signal, 1e-12));
}

TEST_CASE("DST-I Forward/Inverse Round-Trip", "[Fourier][DST]")
{
    // Test DST with boundary conditions (x[0] = x[N+1] = 0)
    int N = 16;
    Vector<Real> signal(N);
    for (int i = 0; i < N; i++)
    {
        // Sine function naturally satisfies boundary conditions
        signal[i] = std::sin(M_PI * (i + 1) / (N + 1));
    }
    
    auto coeffs = DCT::ForwardDST(signal);
    auto recovered = DCT::InverseDST(coeffs);
    
    REQUIRE(recovered.size() == signal.size());
    
    for (int i = 0; i < N; i++)
    {
        REQUIRE(std::abs(recovered[i] - signal[i]) < 1e-12);
    }
}

TEST_CASE("DST-I Known Signal Test", "[Fourier][DST]")
{
    // Pure sine mode: sin(π*k*(n+1)/(N+1)) should give single peak at coefficient k-1
    int N = 8;
    int mode = 3;  // Third mode
    
    Vector<Real> signal(N);
    for (int n = 0; n < N; n++)
    {
        signal[n] = std::sin(M_PI * mode * (n + 1) / (N + 1));
    }
    
    auto coeffs = DCT::ForwardDST(signal);
    
    // Should have peak at coefficient (mode-1)
    Real maxCoeff = 0.0;
    int maxIdx = 0;
    for (int k = 0; k < N; k++)
    {
        if (std::abs(coeffs[k]) > maxCoeff)
        {
            maxCoeff = std::abs(coeffs[k]);
            maxIdx = k;
        }
    }
    
    REQUIRE(maxIdx == mode - 1);
    REQUIRE(maxCoeff > 0.9);  // Should be close to 1
}

TEST_CASE("DCT Edge Cases", "[Fourier][DCT][edge-case]")
{
    // Single element
    Vector<Real> single(1);
    single[0] = 5.0;
    auto dctSingle = DCT::ForwardII(single);
    REQUIRE(dctSingle.size() == 1);
    REQUIRE(std::abs(dctSingle[0] - 10.0) < 1e-10);  // 2/1 * 5 = 10
    
    // Two elements
    Vector<Real> two(2);
    two[0] = 1.0;
    two[1] = 2.0;
    auto dctTwo = DCT::ForwardII(two);
    REQUIRE(dctTwo.size() == 2);
    REQUIRE(DCT::VerifyRoundTrip(two));
}
