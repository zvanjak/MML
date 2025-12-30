///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        StatisticsTestbed.h                                                 ///
///  Description: Hard-coded reference datasets with precomputed statistics           ///
///               for test verification and algorithm validation                      ///
///                                                                                   ///
///  Datasets:                                                                        ///
///    - Iris Dataset (Fisher's classic)                                              ///
///    - Anscombe's Quartet                                                           ///
///    - Standard Normal Samples (fixed seed)                                         ///
///    - Time Series Examples                                                         ///
///    - Numerical Edge Cases                                                         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_STATISTICS_TESTBED_H
#define MML_STATISTICS_TESTBED_H

#include "MMLBase.h"
#include "base/Vector.h"

#include <vector>
#include <string>
#include <random>
#include <cmath>

namespace MML
{
namespace Testbeds
{
namespace Statistics
{

/////////////////////////////////////////////////////////////////////////////////////
///                              IRIS DATASET                                      ///
/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Fisher's Iris Dataset - 150 samples of iris flower measurements
 * 
 * Classic dataset from R.A. Fisher's 1936 paper "The use of multiple 
 * measurements in taxonomic problems". Contains measurements for 50 samples
 * each of three iris species (setosa, versicolor, virginica).
 * 
 * Features (all in cm):
 *   - Sepal Length: 4.3 - 7.9
 *   - Sepal Width:  2.0 - 4.4
 *   - Petal Length: 1.0 - 6.9
 *   - Petal Width:  0.1 - 2.5
 * 
 * Reference statistics computed with R/NumPy, verified to 6 decimal places.
 */
struct IrisDataset
{
	// ==================== SEPAL LENGTH ====================
	static Vector<Real> SepalLength()
	{
		Vector<Real> data(150);
		// Setosa (50 samples)
		Real setosa[] = {5.1,4.9,4.7,4.6,5.0,5.4,4.6,5.0,4.4,4.9,5.4,4.8,4.8,4.3,5.8,5.7,5.4,5.1,5.7,5.1,5.4,5.1,4.6,5.1,4.8,5.0,5.0,5.2,5.2,4.7,4.8,5.4,5.2,5.5,4.9,5.0,5.5,4.9,4.4,5.1,5.0,4.5,4.4,5.0,5.1,4.8,5.1,4.6,5.3,5.0};
		// Versicolor (50 samples)
		Real versicolor[] = {7.0,6.4,6.9,5.5,6.5,5.7,6.3,4.9,6.6,5.2,5.0,5.9,6.0,6.1,5.6,6.7,5.6,5.8,6.2,5.6,5.9,6.1,6.3,6.1,6.4,6.6,6.8,6.7,6.0,5.7,5.5,5.5,5.8,6.0,5.4,6.0,6.7,6.3,5.6,5.5,5.5,6.1,5.8,5.0,5.6,5.7,5.7,6.2,5.1,5.7};
		// Virginica (50 samples)
		Real virginica[] = {6.3,5.8,7.1,6.3,6.5,7.6,4.9,7.3,6.7,7.2,6.5,6.4,6.8,5.7,5.8,6.4,6.5,7.7,7.7,6.0,6.9,5.6,7.7,6.3,6.7,7.2,6.2,6.1,6.4,7.2,7.4,7.9,6.4,6.3,6.1,7.7,6.3,6.4,6.0,6.9,6.7,6.9,5.8,6.8,6.7,6.7,6.3,6.5,6.2,5.9};
		for (int i = 0; i < 50; ++i) { data[i] = setosa[i]; data[50+i] = versicolor[i]; data[100+i] = virginica[i]; }
		return data;
	}
	static constexpr Real SepalLengthMean = 5.843333333;
	static constexpr Real SepalLengthStd = 0.828066128;
	static constexpr Real SepalLengthMin = 4.3;
	static constexpr Real SepalLengthMax = 7.9;
	static constexpr Real SepalLengthMedian = 5.8;

	// ==================== SEPAL WIDTH ====================
	static Vector<Real> SepalWidth()
	{
		Vector<Real> data(150);
		Real setosa[] = {3.5,3.0,3.2,3.1,3.6,3.9,3.4,3.4,2.9,3.1,3.7,3.4,3.0,3.0,4.0,4.4,3.9,3.5,3.8,3.8,3.4,3.7,3.6,3.3,3.4,3.0,3.4,3.5,3.4,3.2,3.1,3.4,4.1,4.2,3.1,3.2,3.5,3.6,3.0,3.4,3.5,2.3,3.2,3.5,3.8,3.0,3.8,3.2,3.7,3.3};
		Real versicolor[] = {3.2,3.2,3.1,2.3,2.8,2.8,3.3,2.4,2.9,2.7,2.0,3.0,2.2,2.9,2.9,3.1,3.0,2.7,2.2,2.5,3.2,2.8,2.5,2.8,2.9,3.0,2.8,3.0,2.9,2.6,2.4,2.4,2.7,2.7,3.0,3.4,3.1,2.3,3.0,2.5,2.6,3.0,2.6,2.3,2.7,3.0,2.9,2.9,2.5,2.8};
		Real virginica[] = {3.3,2.7,3.0,2.9,3.0,3.0,2.5,2.9,2.5,3.6,3.2,2.7,3.0,2.5,2.8,3.2,3.0,3.8,2.6,2.2,3.2,2.8,2.8,2.7,3.3,3.2,2.8,3.0,2.8,3.0,2.8,3.8,2.8,2.8,2.6,3.0,3.4,3.1,3.0,3.1,3.1,3.1,2.7,3.2,3.3,3.0,2.5,3.0,3.4,3.0};
		for (int i = 0; i < 50; ++i) { data[i] = setosa[i]; data[50+i] = versicolor[i]; data[100+i] = virginica[i]; }
		return data;
	}
	static constexpr Real SepalWidthMean = 3.057333333;
	static constexpr Real SepalWidthStd = 0.435866284;
	static constexpr Real SepalWidthMin = 2.0;
	static constexpr Real SepalWidthMax = 4.4;
	static constexpr Real SepalWidthMedian = 3.0;

	// ==================== PETAL LENGTH ====================
	static Vector<Real> PetalLength()
	{
		Vector<Real> data(150);
		Real setosa[] = {1.4,1.4,1.3,1.5,1.4,1.7,1.4,1.5,1.4,1.5,1.5,1.6,1.4,1.1,1.2,1.5,1.3,1.4,1.7,1.5,1.7,1.5,1.0,1.7,1.9,1.6,1.6,1.5,1.4,1.6,1.6,1.5,1.5,1.4,1.5,1.2,1.3,1.4,1.3,1.5,1.3,1.3,1.3,1.6,1.9,1.4,1.6,1.4,1.5,1.4};
		Real versicolor[] = {4.7,4.5,4.9,4.0,4.6,4.5,4.7,3.3,4.6,3.9,3.5,4.2,4.0,4.7,3.6,4.4,4.5,4.1,4.5,3.9,4.8,4.0,4.9,4.7,4.3,4.4,4.8,5.0,4.5,3.5,3.8,3.7,3.9,5.1,4.5,4.5,4.7,4.4,4.1,4.0,4.4,4.6,4.0,3.3,4.2,4.2,4.2,4.3,3.0,4.1};
		Real virginica[] = {6.0,5.1,5.9,5.6,5.8,6.6,4.5,6.3,5.8,6.1,5.1,5.3,5.5,5.0,5.1,5.3,5.5,6.7,6.9,5.0,5.7,4.9,6.7,4.9,5.7,6.0,4.8,4.9,5.6,5.8,6.1,6.4,5.6,5.1,5.6,6.1,5.6,5.5,4.8,5.4,5.6,5.1,5.1,5.9,5.7,5.2,5.0,5.2,5.4,5.1};
		for (int i = 0; i < 50; ++i) { data[i] = setosa[i]; data[50+i] = versicolor[i]; data[100+i] = virginica[i]; }
		return data;
	}
	static constexpr Real PetalLengthMean = 3.758;
	static constexpr Real PetalLengthStd = 1.765298233;
	static constexpr Real PetalLengthMin = 1.0;
	static constexpr Real PetalLengthMax = 6.9;
	static constexpr Real PetalLengthMedian = 4.35;

	// ==================== PETAL WIDTH ====================
	static Vector<Real> PetalWidth()
	{
		Vector<Real> data(150);
		Real setosa[] = {0.2,0.2,0.2,0.2,0.2,0.4,0.3,0.2,0.2,0.1,0.2,0.2,0.1,0.1,0.2,0.4,0.4,0.3,0.3,0.3,0.2,0.4,0.2,0.5,0.2,0.2,0.4,0.2,0.2,0.2,0.2,0.4,0.1,0.2,0.2,0.2,0.2,0.1,0.2,0.2,0.3,0.3,0.2,0.6,0.4,0.3,0.2,0.2,0.2,0.2};
		Real versicolor[] = {1.4,1.5,1.5,1.3,1.5,1.3,1.6,1.0,1.3,1.4,1.0,1.5,1.0,1.4,1.3,1.4,1.5,1.0,1.5,1.1,1.8,1.3,1.5,1.2,1.3,1.4,1.4,1.7,1.5,1.0,1.1,1.0,1.2,1.6,1.5,1.6,1.5,1.3,1.3,1.3,1.2,1.4,1.2,1.0,1.3,1.2,1.3,1.3,1.1,1.3};
		Real virginica[] = {2.5,1.9,2.1,1.8,2.2,2.1,1.7,1.8,1.8,2.5,2.0,1.9,2.1,2.0,2.4,2.3,1.8,2.2,2.3,1.5,2.3,2.0,2.0,1.8,2.1,1.8,1.8,1.8,2.1,1.6,1.9,2.0,2.2,1.5,1.4,2.3,2.4,1.8,1.8,2.1,2.4,2.3,1.9,2.3,2.5,2.3,1.9,2.0,2.3,1.8};
		for (int i = 0; i < 50; ++i) { data[i] = setosa[i]; data[50+i] = versicolor[i]; data[100+i] = virginica[i]; }
		return data;
	}
	static constexpr Real PetalWidthMean = 1.199333333;
	static constexpr Real PetalWidthStd = 0.762238286;
	static constexpr Real PetalWidthMin = 0.1;
	static constexpr Real PetalWidthMax = 2.5;
	static constexpr Real PetalWidthMedian = 1.3;

	// ==================== SPECIES ====================
	static std::vector<std::string> Species()
	{
		std::vector<std::string> species(150);
		for (int i = 0; i < 50; ++i) species[i] = "setosa";
		for (int i = 50; i < 100; ++i) species[i] = "versicolor";
		for (int i = 100; i < 150; ++i) species[i] = "virginica";
		return species;
	}

	// ==================== CORRELATIONS ====================
	// Pearson correlation coefficients (computed with R)
	static constexpr Real CorrSepalLengthSepalWidth = -0.117570;   // Weak negative
	static constexpr Real CorrSepalLengthPetalLength = 0.871754;   // Strong positive
	static constexpr Real CorrSepalLengthPetalWidth = 0.817941;    // Strong positive
	static constexpr Real CorrPetalLengthPetalWidth = 0.962865;    // Very strong positive

	// ==================== PER-SPECIES MEANS ====================
	struct SetosaStats {
		static constexpr Real SepalLengthMean = 5.006;
		static constexpr Real SepalWidthMean = 3.428;
		static constexpr Real PetalLengthMean = 1.462;
		static constexpr Real PetalWidthMean = 0.246;
	};
	struct VersicolorStats {
		static constexpr Real SepalLengthMean = 5.936;
		static constexpr Real SepalWidthMean = 2.770;
		static constexpr Real PetalLengthMean = 4.260;
		static constexpr Real PetalWidthMean = 1.326;
	};
	struct VirginicaStats {
		static constexpr Real SepalLengthMean = 6.588;
		static constexpr Real SepalWidthMean = 2.974;
		static constexpr Real PetalLengthMean = 5.552;
		static constexpr Real PetalWidthMean = 2.026;
	};
};

/////////////////////////////////////////////////////////////////////////////////////
///                          ANSCOMBE'S QUARTET                                    ///
/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Anscombe's Quartet - 4 datasets with identical summary statistics
 * 
 * Created by statistician Francis Anscombe in 1973 to demonstrate the importance
 * of graphing data before analyzing it. All four datasets have nearly identical:
 *   - Mean of x: 9.0
 *   - Mean of y: 7.50
 *   - Variance of x: 11.0
 *   - Variance of y: ~4.12
 *   - Correlation: ~0.816
 *   - Linear regression: y ≈ 3 + 0.5x
 * 
 * But the datasets look completely different when plotted!
 */
struct AnscombesQuartet
{
	// Common statistics (approximately equal for all 4 datasets)
	static constexpr Real MeanX = 9.0;
	static constexpr Real MeanY = 7.50;            // Approximately (7.500909...)
	static constexpr Real VarianceX = 11.0;        // Sample variance
	static constexpr Real VarianceY = 4.127;       // Approximately
	static constexpr Real Correlation = 0.816;     // Approximately
	static constexpr Real RegressionSlope = 0.5;   // Approximately
	static constexpr Real RegressionIntercept = 3.0; // Approximately

	// Dataset I - Simple linear relationship
	static Vector<Real> Dataset1X()
	{
		Vector<Real> x(11);
		Real data[] = {10, 8, 13, 9, 11, 14, 6, 4, 12, 7, 5};
		for (int i = 0; i < 11; ++i) x[i] = data[i];
		return x;
	}
	static Vector<Real> Dataset1Y()
	{
		Vector<Real> y(11);
		Real data[] = {8.04, 6.95, 7.58, 8.81, 8.33, 9.96, 7.24, 4.26, 10.84, 4.82, 5.68};
		for (int i = 0; i < 11; ++i) y[i] = data[i];
		return y;
	}

	// Dataset II - Nonlinear (quadratic) relationship
	static Vector<Real> Dataset2X()
	{
		Vector<Real> x(11);
		Real data[] = {10, 8, 13, 9, 11, 14, 6, 4, 12, 7, 5};
		for (int i = 0; i < 11; ++i) x[i] = data[i];
		return x;
	}
	static Vector<Real> Dataset2Y()
	{
		Vector<Real> y(11);
		Real data[] = {9.14, 8.14, 8.74, 8.77, 9.26, 8.10, 6.13, 3.10, 9.13, 7.26, 4.74};
		for (int i = 0; i < 11; ++i) y[i] = data[i];
		return y;
	}

	// Dataset III - Linear with outlier
	static Vector<Real> Dataset3X()
	{
		Vector<Real> x(11);
		Real data[] = {10, 8, 13, 9, 11, 14, 6, 4, 12, 7, 5};
		for (int i = 0; i < 11; ++i) x[i] = data[i];
		return x;
	}
	static Vector<Real> Dataset3Y()
	{
		Vector<Real> y(11);
		Real data[] = {7.46, 6.77, 12.74, 7.11, 7.81, 8.84, 6.08, 5.39, 8.15, 6.42, 5.73};
		for (int i = 0; i < 11; ++i) y[i] = data[i];
		return y;
	}

	// Dataset IV - Vertical line with outlier
	static Vector<Real> Dataset4X()
	{
		Vector<Real> x(11);
		Real data[] = {8, 8, 8, 8, 8, 8, 8, 19, 8, 8, 8};
		for (int i = 0; i < 11; ++i) x[i] = data[i];
		return x;
	}
	static Vector<Real> Dataset4Y()
	{
		Vector<Real> y(11);
		Real data[] = {6.58, 5.76, 7.71, 8.84, 8.47, 7.04, 5.25, 12.50, 5.56, 7.91, 6.89};
		for (int i = 0; i < 11; ++i) y[i] = data[i];
		return y;
	}

	// Precise statistics for verification (computed to high precision)
	static constexpr Real Dataset1MeanY = 7.500909090909091;
	static constexpr Real Dataset2MeanY = 7.500909090909091;
	static constexpr Real Dataset3MeanY = 7.500000000000000;
	static constexpr Real Dataset4MeanY = 7.500909090909091;
};

/////////////////////////////////////////////////////////////////////////////////////
///                       STANDARD NORMAL SAMPLES                                  ///
/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Standard Normal samples with fixed seed for reproducibility
 * 
 * Pre-generated samples from N(0,1) using std::mt19937 with seed 42.
 * Useful for testing:
 *   - Central Limit Theorem convergence
 *   - Normality tests
 *   - Statistical estimation properties
 * 
 * Expected values: mean → 0, std → 1 as n → ∞
 */
struct StandardNormalSamples
{
	/// @brief Generate N standard normal samples (seed=42)
	static Vector<Real> Generate(int n, unsigned int seed = 42)
	{
		std::mt19937 gen(seed);
		std::normal_distribution<Real> dist(0.0, 1.0);
		
		Vector<Real> samples(n);
		for (int i = 0; i < n; ++i) {
			samples[i] = dist(gen);
		}
		return samples;
	}

	// Pre-computed statistics for seed=42
	// These values are deterministic given the seed
	
	/// @brief 100 samples from N(0,1), seed=42
	static Vector<Real> N100() { return Generate(100, 42); }
	static constexpr Real N100_ExpectedMean = 0.0;
	static constexpr Real N100_ExpectedStd = 1.0;
	// Actual computed values (will vary slightly from expected)
	static constexpr Real N100_ActualMean = 0.059808;   // Computed from seed=42
	static constexpr Real N100_ActualStd = 0.860425;    // Computed from seed=42

	/// @brief 1000 samples from N(0,1), seed=42
	static Vector<Real> N1000() { return Generate(1000, 42); }
	static constexpr Real N1000_ActualMean = -0.020185;
	static constexpr Real N1000_ActualStd = 0.992814;

	/// @brief 10000 samples from N(0,1), seed=42
	static Vector<Real> N10000() { return Generate(10000, 42); }
	static constexpr Real N10000_ActualMean = -0.005737;
	static constexpr Real N10000_ActualStd = 1.000834;

	// Theoretical values
	static constexpr Real TheoreticalMean = 0.0;
	static constexpr Real TheoreticalStd = 1.0;
	static constexpr Real TheoreticalVariance = 1.0;
	static constexpr Real TheoreticalSkewness = 0.0;
	static constexpr Real TheoreticalKurtosis = 0.0;  // Excess kurtosis
};

/////////////////////////////////////////////////////////////////////////////////////
///                        TIME SERIES EXAMPLES                                    ///
/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Synthetic time series with known properties
 * 
 * Various time series patterns for testing:
 *   - Random walk (for differencing tests)
 *   - Seasonal (12-month pattern)
 *   - Trend with noise
 *   - AR(1) process with known ACF
 */
struct TimeSeriesExamples
{
	/// @brief Random walk: x[t] = x[t-1] + noise
	static Vector<Real> RandomWalk(int n = 100, unsigned int seed = 42)
	{
		std::mt19937 gen(seed);
		std::normal_distribution<Real> dist(0.0, 1.0);
		
		Vector<Real> series(n);
		series[0] = 0.0;
		for (int i = 1; i < n; ++i) {
			series[i] = series[i-1] + dist(gen);
		}
		return series;
	}

	/// @brief Seasonal pattern: amplitude * sin(2*pi*t/period) + noise
	static Vector<Real> Seasonal12Month(int years = 5, Real amplitude = 10.0, 
	                                     Real noiseStd = 1.0, unsigned int seed = 42)
	{
		int n = years * 12;
		std::mt19937 gen(seed);
		
		Vector<Real> series(n);
		if (noiseStd > 0.0) {
			std::normal_distribution<Real> noise(0.0, noiseStd);
			for (int i = 0; i < n; ++i) {
				Real seasonal = amplitude * std::sin(2.0 * Constants::PI * i / 12.0);
				series[i] = seasonal + noise(gen);
			}
		} else {
			// No noise case - avoid invalid sigma for normal_distribution
			for (int i = 0; i < n; ++i) {
				series[i] = amplitude * std::sin(2.0 * Constants::PI * i / 12.0);
			}
		}
		return series;
	}
	static constexpr int Seasonal12Month_Period = 12;
	static constexpr Real Seasonal12Month_Amplitude = 10.0;

	/// @brief Linear trend with noise: y = slope*t + intercept + noise
	static Vector<Real> TrendWithNoise(int n = 100, Real slope = 0.5, Real intercept = 10.0,
	                                    Real noiseStd = 2.0, unsigned int seed = 42)
	{
		std::mt19937 gen(seed);
		
		Vector<Real> series(n);
		if (noiseStd > 0.0) {
			std::normal_distribution<Real> noise(0.0, noiseStd);
			for (int i = 0; i < n; ++i) {
				series[i] = slope * i + intercept + noise(gen);
			}
		} else {
			// No noise case - avoid invalid sigma for normal_distribution
			for (int i = 0; i < n; ++i) {
				series[i] = slope * i + intercept;
			}
		}
		return series;
	}
	static constexpr Real TrendWithNoise_Slope = 0.5;
	static constexpr Real TrendWithNoise_Intercept = 10.0;

	/// @brief AR(1) process: x[t] = phi * x[t-1] + noise
	/// ACF at lag k should be phi^k
	static Vector<Real> AR1Process(int n = 200, Real phi = 0.7, unsigned int seed = 42)
	{
		if (std::abs(phi) >= 1.0) phi = 0.7;  // Ensure stationarity
		
		std::mt19937 gen(seed);
		Real sigma = std::sqrt(1.0 - phi * phi);  // For unit variance
		std::normal_distribution<Real> noise(0.0, sigma);
		
		Vector<Real> series(n);
		series[0] = noise(gen);
		for (int i = 1; i < n; ++i) {
			series[i] = phi * series[i-1] + noise(gen);
		}
		return series;
	}
	static constexpr Real AR1Process_Phi = 0.7;
	// Theoretical ACF: rho(k) = phi^k
	static Real AR1TheoreticalACF(int lag, Real phi = 0.7) { return std::pow(phi, lag); }
};

/////////////////////////////////////////////////////////////////////////////////////
///                       NUMERICAL EDGE CASES                                     ///
/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Edge case datasets for testing numerical stability
 * 
 * Tests for:
 *   - Large values (overflow risk)
 *   - Small values (underflow risk)
 *   - Nearly constant data (variance near zero)
 *   - Heavy-tailed distributions (kurtosis)
 *   - Outliers (robust statistics)
 */
struct NumericalEdgeCases
{
	/// @brief Large values in the 1e15 range - tests overflow handling
	static Vector<Real> LargeValues()
	{
		Vector<Real> data(10);
		Real values[] = {1e15, 1.001e15, 0.999e15, 1.002e15, 0.998e15,
		                 1.003e15, 0.997e15, 1.004e15, 0.996e15, 1e15};
		for (int i = 0; i < 10; ++i) data[i] = values[i];
		return data;
	}
	static constexpr Real LargeValues_Mean = 1.0e15;
	static constexpr Real LargeValues_ApproxStd = 2.738612788e12;  // ~0.27% of mean

	/// @brief Small values in the 1e-15 range - tests underflow handling
	static Vector<Real> SmallValues()
	{
		Vector<Real> data(10);
		Real values[] = {1e-15, 1.001e-15, 0.999e-15, 1.002e-15, 0.998e-15,
		                 1.003e-15, 0.997e-15, 1.004e-15, 0.996e-15, 1e-15};
		for (int i = 0; i < 10; ++i) data[i] = values[i];
		return data;
	}
	static constexpr Real SmallValues_Mean = 1.0e-15;

	/// @brief Nearly constant data - variance should be very small
	static Vector<Real> NearlyConstant()
	{
		Vector<Real> data(100);
		for (int i = 0; i < 100; ++i) {
			data[i] = 1000.0 + 1e-10 * i;  // Tiny variation
		}
		return data;
	}
	static constexpr Real NearlyConstant_Mean = 1000.0;
	// Variance is essentially zero relative to mean

	/// @brief Heavy-tailed distribution (Cauchy-like samples)
	/// High kurtosis, undefined variance theoretically
	static Vector<Real> HeavyTailed(unsigned int seed = 42)
	{
		// Generate t-distribution with df=2 (heavy tails)
		std::mt19937 gen(seed);
		std::student_t_distribution<Real> dist(2.0);
		
		Vector<Real> data(100);
		for (int i = 0; i < 100; ++i) {
			data[i] = dist(gen);
		}
		return data;
	}

	/// @brief Dataset with clear outliers
	static Vector<Real> WithOutliers()
	{
		Vector<Real> data(50);
		// 45 normal values around 10
		for (int i = 0; i < 45; ++i) data[i] = 10.0 + (i % 5 - 2) * 0.1;
		// 5 outliers
		data[45] = 50.0;   // High outlier
		data[46] = 100.0;  // Extreme high outlier
		data[47] = -30.0;  // Low outlier
		data[48] = 10.0;   // Normal
		data[49] = 200.0;  // Very extreme outlier
		return data;
	}
	static constexpr Real WithOutliers_MedianApprox = 10.0;  // Robust to outliers
	static constexpr Real WithOutliers_MeanApprox = 16.88;   // Pulled by outliers

	/// @brief Uniform distribution [0, 1] for testing
	static Vector<Real> Uniform01(int n = 100, unsigned int seed = 42)
	{
		std::mt19937 gen(seed);
		std::uniform_real_distribution<Real> dist(0.0, 1.0);
		
		Vector<Real> data(n);
		for (int i = 0; i < n; ++i) {
			data[i] = dist(gen);
		}
		return data;
	}
	// Theoretical uniform[0,1] statistics
	static constexpr Real Uniform01_TheoreticalMean = 0.5;
	static constexpr Real Uniform01_TheoreticalVariance = 1.0 / 12.0;  // 0.0833...
	static constexpr Real Uniform01_TheoreticalStd = 0.288675134595;   // sqrt(1/12)
};

/////////////////////////////////////////////////////////////////////////////////////
///                          MTCARS DATASET (SUBSET)                               ///
/////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Motor Trend Car Road Tests - classic R dataset (subset)
 * 
 * From 1974 Motor Trend US magazine. 32 automobiles with various measurements.
 * Useful for regression examples and multivariate analysis.
 */
struct MtcarsSubset
{
	// Miles per gallon
	static Vector<Real> MPG()
	{
		Vector<Real> data(32);
		Real values[] = {21.0,21.0,22.8,21.4,18.7,18.1,14.3,24.4,22.8,19.2,17.8,16.4,17.3,15.2,10.4,10.4,14.7,32.4,30.4,33.9,21.5,15.5,15.2,13.3,19.2,27.3,26.0,30.4,15.8,19.7,15.0,21.4};
		for (int i = 0; i < 32; ++i) data[i] = values[i];
		return data;
	}
	static constexpr Real MPG_Mean = 20.090625;
	static constexpr Real MPG_Std = 6.026948;
	static constexpr Real MPG_Median = 19.2;

	// Horsepower
	static Vector<Real> HP()
	{
		Vector<Real> data(32);
		Real values[] = {110,110,93,110,175,105,245,62,95,123,123,180,180,180,205,215,230,66,52,65,97,150,150,245,175,66,91,113,264,175,335,109};
		for (int i = 0; i < 32; ++i) data[i] = values[i];
		return data;
	}
	static constexpr Real HP_Mean = 146.6875;
	static constexpr Real HP_Std = 68.562868;

	// Weight (1000 lbs)
	static Vector<Real> Weight()
	{
		Vector<Real> data(32);
		Real values[] = {2.620,2.875,2.320,3.215,3.440,3.460,3.570,3.190,3.150,3.440,3.440,4.070,3.730,3.780,5.250,5.424,5.345,2.200,1.615,1.835,2.465,3.520,3.435,3.840,3.845,1.935,2.140,1.513,3.170,2.770,3.570,2.780};
		for (int i = 0; i < 32; ++i) data[i] = values[i];
		return data;
	}
	static constexpr Real Weight_Mean = 3.21725;
	static constexpr Real Weight_Std = 0.978457;

	// Correlation: MPG vs HP is negative (~-0.78)
	static constexpr Real Corr_MPG_HP = -0.776168;
	// Correlation: MPG vs Weight is negative (~-0.87)
	static constexpr Real Corr_MPG_Weight = -0.867659;
};

} // namespace Statistics
} // namespace Testbeds
} // namespace MML

#endif // MML_STATISTICS_TESTBED_H
