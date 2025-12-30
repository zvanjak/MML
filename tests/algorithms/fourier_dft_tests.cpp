#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "../../mml/algorithms/Fourier/DFT.h"

using namespace MML;

TEST_CASE("DFT - Basic functionality", "[fourier][dft]")
{
    SECTION("Forward DFT of impulse at t=0")
    {
        // Impulse at t=0 should give constant spectrum (all 1s)
        Vector<Complex> data(4);
        data[0] = Complex(1.0, 0.0);
        data[1] = Complex(0.0, 0.0);
        data[2] = Complex(0.0, 0.0);
        data[3] = Complex(0.0, 0.0);

        Vector<Complex> spectrum = DFT::Forward(data);

        // All frequency bins should be 1.0
        for (int k = 0; k < 4; k++) {
            REQUIRE(std::abs(spectrum[k] - Complex(1.0, 0.0)) < 1e-10);
        }
    }

    SECTION("Forward DFT of constant signal")
    {
        // Constant signal should give impulse at DC (k=0)
        Vector<Complex> data(4, Complex(1.0, 0.0));

        Vector<Complex> spectrum = DFT::Forward(data);

        // DC component should be N
        REQUIRE(std::abs(spectrum[0] - Complex(4.0, 0.0)) < 1e-10);
        
        // All other bins should be zero
        for (int k = 1; k < 4; k++) {
            REQUIRE(std::abs(spectrum[k]) < 1e-10);
        }
    }

    SECTION("Inverse DFT recovers original signal")
    {
        Vector<Complex> original(8);
        for (int i = 0; i < 8; i++) {
            original[i] = Complex(std::sin(i * 0.5), std::cos(i * 0.3));
        }

        Vector<Complex> spectrum = DFT::Forward(original);
        Vector<Complex> recovered = DFT::Inverse(spectrum);

        for (int i = 0; i < 8; i++) {
            REQUIRE(std::abs(recovered[i] - original[i]) < 1e-10);
        }
    }

    SECTION("Real-valued input version")
    {
        Vector<Real> real_data{1.0, 2.0, 3.0, 4.0};
        
        Vector<Complex> spectrum = DFT::Forward(real_data);
        Vector<Real> recovered = DFT::InverseReal(spectrum);

        for (int i = 0; i < 4; i++) {
            REQUIRE(std::abs(recovered[i] - real_data[i]) < 1e-10);
        }
    }
}

TEST_CASE("DFT - Pure sinusoid analysis", "[fourier][dft]")
{
    SECTION("Single frequency sinusoid")
    {
        // Create cos(2π*k0*n/N) signal - should have peaks at k=k0 and k=N-k0
        int N = 16;
        int k0 = 3;  // Frequency bin
        Vector<Complex> data(N);
        
        for (int n = 0; n < N; n++) {
            Real angle = 2.0 * Constants::PI * k0 * n / N;
            data[n] = Complex(std::cos(angle), 0.0);
        }

        Vector<Complex> spectrum = DFT::Forward(data);

        // Should have peaks at k=3 and k=13 (N-k0 = 16-3 = 13)
        // Each peak should have magnitude N/2
        REQUIRE(std::abs(spectrum[k0]) > N/2 - 1e-6);
        REQUIRE(std::abs(spectrum[N - k0]) > N/2 - 1e-6);

        // Other bins should be near zero
        for (int k = 0; k < N; k++) {
            if (k != k0 && k != (N - k0)) {
                REQUIRE(std::abs(spectrum[k]) < 1e-10);
            }
        }
    }

    SECTION("Sine wave - imaginary spectrum")
    {
        // Create sin(2π*k0*n/N) signal
        int N = 16;
        int k0 = 4;
        Vector<Complex> data(N);
        
        for (int n = 0; n < N; n++) {
            Real angle = 2.0 * Constants::PI * k0 * n / N;
            data[n] = Complex(std::sin(angle), 0.0);
        }

        Vector<Complex> spectrum = DFT::Forward(data);

        // Sine wave produces imaginary spectrum
        // Peak at k=k0 should be -i*N/2 (negative imaginary)
        // Peak at k=N-k0 should be +i*N/2 (positive imaginary)
        
        // Check magnitude
        REQUIRE(std::abs(spectrum[k0]) > N/2 - 1e-6);
        REQUIRE(std::abs(spectrum[N - k0]) > N/2 - 1e-6);
        
        // Check it's purely imaginary at peak
        REQUIRE(std::abs(spectrum[k0].real()) < 1e-10);
        REQUIRE(std::abs(spectrum[N - k0].real()) < 1e-10);
    }
}

TEST_CASE("DFT - Mathematical properties", "[fourier][dft]")
{
    SECTION("Parseval's theorem - energy conservation")
    {
        // Sum of |x[n]|² = (1/N) * Sum of |X[k]|²
        Vector<Complex> signal(16);
        for (int i = 0; i < 16; i++) {
            signal[i] = Complex(std::sin(i * 0.3), std::cos(i * 0.5));
        }

        Vector<Complex> spectrum = DFT::Forward(signal);

        Real time_energy = 0.0;
        for (int i = 0; i < 16; i++) {
            time_energy += std::norm(signal[i]);  // |x[n]|²
        }

        Real freq_energy = 0.0;
        for (int k = 0; k < 16; k++) {
            freq_energy += std::norm(spectrum[k]);  // |X[k]|²
        }
        freq_energy /= 16;  // Normalize by N

        REQUIRE(std::abs(time_energy - freq_energy) < 1e-8);
    }

    SECTION("Linearity - DFT(a*x + b*y) = a*DFT(x) + b*DFT(y)")
    {
        Vector<Complex> x(8), y(8);
        for (int i = 0; i < 8; i++) {
            x[i] = Complex(i * 0.5, 0.0);
            y[i] = Complex(0.0, i * 0.3);
        }

        Complex a(2.0, 1.0);
        Complex b(1.0, -1.0);

        // Compute DFT(a*x + b*y) directly
        Vector<Complex> combined(8);
        for (int i = 0; i < 8; i++) {
            combined[i] = a * x[i] + b * y[i];
        }
        Vector<Complex> spectrum_combined = DFT::Forward(combined);

        // Compute a*DFT(x) + b*DFT(y)
        Vector<Complex> spectrum_x = DFT::Forward(x);
        Vector<Complex> spectrum_y = DFT::Forward(y);
        Vector<Complex> spectrum_linear(8);
        for (int k = 0; k < 8; k++) {
            spectrum_linear[k] = a * spectrum_x[k] + b * spectrum_y[k];
        }

        // Should be equal
        for (int k = 0; k < 8; k++) {
            REQUIRE(std::abs(spectrum_combined[k] - spectrum_linear[k]) < 1e-10);
        }
    }

    SECTION("Conjugate symmetry for real signals")
    {
        // Real signal -> X[k] = X*[N-k] (conjugate symmetry)
        Vector<Real> real_signal{1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 1.0, 0.0};
        Vector<Complex> spectrum = DFT::Forward(real_signal);

        int N = spectrum.size();
        for (int k = 1; k < N/2; k++) {
            Complex conj_diff = spectrum[k] - std::conj(spectrum[N - k]);
            REQUIRE(std::abs(conj_diff) < 1e-10);
        }

        // DC and Nyquist bins should be real
        REQUIRE(std::abs(spectrum[0].imag()) < 1e-10);
        if (N % 2 == 0) {
            REQUIRE(std::abs(spectrum[N/2].imag()) < 1e-10);
        }
    }
}

TEST_CASE("DFT - Edge cases", "[fourier][dft]")
{
    SECTION("Single element")
    {
        Vector<Complex> data{Complex(42.0, 17.0)};
        Vector<Complex> spectrum = DFT::Forward(data);
        
        REQUIRE(spectrum.size() == 1);
        REQUIRE(std::abs(spectrum[0] - Complex(42.0, 17.0)) < 1e-10);

        Vector<Complex> recovered = DFT::Inverse(spectrum);
        REQUIRE(std::abs(recovered[0] - Complex(42.0, 17.0)) < 1e-10);
    }

    SECTION("Two elements")
    {
        Vector<Complex> data{Complex(1.0, 0.0), Complex(-1.0, 0.0)};
        Vector<Complex> spectrum = DFT::Forward(data);
        Vector<Complex> recovered = DFT::Inverse(spectrum);

        for (int i = 0; i < 2; i++) {
            REQUIRE(std::abs(recovered[i] - data[i]) < 1e-10);
        }
    }

    SECTION("Empty input throws")
    {
        Vector<Complex> empty;
        REQUIRE_THROWS_AS(DFT::Forward(empty), std::invalid_argument);
        REQUIRE_THROWS_AS(DFT::Inverse(empty), std::invalid_argument);
    }
}

TEST_CASE("DFT - Utility functions", "[fourier][dft]")
{
    SECTION("IsReasonableSize")
    {
        REQUIRE(DFT::IsReasonableSize(1024));
        REQUIRE(DFT::IsReasonableSize(4096));
        REQUIRE_FALSE(DFT::IsReasonableSize(8192));
    }

    SECTION("ComputeBin - single frequency bin")
    {
        Vector<Complex> data(8);
        for (int i = 0; i < 8; i++) {
            data[i] = Complex(i * 0.5, 0.0);
        }

        Vector<Complex> full_spectrum = DFT::Forward(data);

        for (int k = 0; k < 8; k++) {
            Complex bin = DFT::ComputeBin(data, k);
            REQUIRE(std::abs(bin - full_spectrum[k]) < 1e-10);
        }
    }

    SECTION("ComputeBin - invalid index throws")
    {
        Vector<Complex> data(8, Complex(1.0, 0.0));
        
        REQUIRE_THROWS_AS(DFT::ComputeBin(data, -1), std::invalid_argument);
        REQUIRE_THROWS_AS(DFT::ComputeBin(data, 8), std::invalid_argument);
    }
}
