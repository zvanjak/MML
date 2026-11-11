#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>

#include <mml/base/Vector/Vector.h>
#include <mml/algorithms/Statistics.h>
#include <mml/algorithms/Statistics/Distributions.h>
#include <mml/algorithms/Statistics/DiscreteDistributions.h>
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                    BASIC DESCRIPTIVE STATISTICS                                     ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Statistics_Basic()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Basic Descriptive Statistics\n";
	std::cout << "==========================================================================\n";
	
	// Sample data
	Vector<Real> data({2.5, 3.1, 4.7, 5.2, 3.8, 4.1, 5.9, 2.9, 4.3, 3.6});
	
	std::cout << "\nSample data (n=10): ";
	for (int i = 0; i < data.size(); i++)
		std::cout << data[i] << " ";
	std::cout << std::endl;
	
	// Mean, variance, standard deviation
	Real avg, var;
	Statistics::AvgVar(data, avg, var);
	Real stdDev = Statistics::StdDev(data);
	
	std::cout << "\n--- Central Tendency & Dispersion ---\n";
	std::cout << "Mean:               " << avg << std::endl;
	std::cout << "Sample Variance:    " << var << std::endl;
	std::cout << "Sample Std Dev:     " << stdDev << std::endl;
	std::cout << "Population Variance:" << Statistics::PopulationVariance(data) << std::endl;
	std::cout << "Population Std Dev: " << Statistics::PopulationStdDev(data) << std::endl;
	
	// Median and range
	Real median = Statistics::Median(data);
	Real range = Statistics::Range(data);
	
	std::cout << "Median:             " << median << std::endl;
	std::cout << "Range:              " << range << std::endl;
}

void Docs_Demo_Statistics_Moments()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Statistical Moments\n";
	std::cout << "==========================================================================\n";
	
	// Create dataset with some skew
	Vector<Real> data({1.0, 2.0, 2.5, 3.0, 3.2, 3.5, 4.0, 4.5, 5.0, 8.0});
	
	std::cout << "\nData with positive skew: ";
	for (int i = 0; i < data.size(); i++)
		std::cout << data[i] << " ";
	std::cout << std::endl;
	
	Real ave, adev, sdev, var, skew, curt;
	Statistics::Moments(data, ave, adev, sdev, var, skew, curt);
	
	std::cout << "\n--- Statistical Moments ---\n";
	std::cout << "Mean:               " << ave << std::endl;
	std::cout << "Average Deviation:  " << adev << std::endl;
	std::cout << "Standard Deviation: " << sdev << std::endl;
	std::cout << "Variance:           " << var << std::endl;
	std::cout << "Skewness:           " << skew << " (>0 = right-skewed)" << std::endl;
	std::cout << "Kurtosis:           " << curt << " (>0 = heavy-tailed)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    ORDER STATISTICS                                                 ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Statistics_Order()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Order Statistics\n";
	std::cout << "==========================================================================\n";
	
	Vector<Real> data({12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45});
	
	std::cout << "\nSorted data (n=12): ";
	for (int i = 0; i < data.size(); i++)
		std::cout << data[i] << " ";
	std::cout << std::endl;
	
	// Quartiles
	Real q1, median, q3;
	Statistics::Quartiles(data, q1, median, q3);
	
	std::cout << "\n--- Quartiles ---\n";
	std::cout << "Q1 (25th percentile): " << q1 << std::endl;
	std::cout << "Q2 (Median):          " << median << std::endl;
	std::cout << "Q3 (75th percentile): " << q3 << std::endl;
	std::cout << "IQR (Q3-Q1):          " << Statistics::IQR(data) << std::endl;
	
	// Percentiles
	std::cout << "\n--- Percentiles ---\n";
	std::cout << "10th percentile: " << Statistics::Percentile(data, 10) << std::endl;
	std::cout << "90th percentile: " << Statistics::Percentile(data, 90) << std::endl;
	
	// Min/Max
	Real minVal, maxVal;
	Statistics::MinMax(data, minVal, maxVal);
	std::cout << "\n--- Range ---\n";
	std::cout << "Min: " << minVal << ", Max: " << maxVal << std::endl;
	std::cout << "Range: " << Statistics::Range(data) << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    ROBUST STATISTICS                                                ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Statistics_Robust()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Robust Statistics\n";
	std::cout << "==========================================================================\n";
	
	// Data with outliers
	Vector<Real> data({10, 12, 11, 13, 12, 11, 14, 100});  // 100 is outlier
	
	std::cout << "\nData with outlier (100): ";
	for (int i = 0; i < data.size(); i++)
		std::cout << data[i] << " ";
	std::cout << std::endl;
	
	std::cout << "\n--- Comparison: Standard vs Robust ---\n";
	std::cout << "Mean:         " << Statistics::Mean(data) << " (sensitive to outlier)" << std::endl;
	std::cout << "Median:       " << Statistics::Median(data) << " (robust)" << std::endl;
	std::cout << "Trimmed Mean (10%): " << Statistics::TrimmedMean(data, 10.0) << std::endl;
	
	std::cout << "\n--- Robust Dispersion ---\n";
	std::cout << "StdDev:       " << Statistics::StdDev(data) << " (inflated by outlier)" << std::endl;
	std::cout << "MAD:          " << Statistics::MAD(data) << " (robust)" << std::endl;
	std::cout << "IQR:          " << Statistics::IQR(data) << " (robust)" << std::endl;
}

void Docs_Demo_Statistics_Means()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Different Types of Means\n";
	std::cout << "==========================================================================\n";
	
	// Growth rates (multiplicative data)
	Vector<Real> growthRates({1.05, 1.10, 0.95, 1.08, 1.12});  // 5%, 10%, -5%, 8%, 12%
	
	std::cout << "\nGrowth factors: ";
	for (int i = 0; i < growthRates.size(); i++)
		std::cout << growthRates[i] << " ";
	std::cout << std::endl;
	
	std::cout << "\n--- Types of Averages ---\n";
	std::cout << "Arithmetic Mean: " << Statistics::Mean(growthRates) << std::endl;
	std::cout << "Geometric Mean:  " << Statistics::GeometricMean(growthRates) 
	          << " (appropriate for growth rates)" << std::endl;
	std::cout << "Harmonic Mean:   " << Statistics::HarmonicMean(growthRates) << std::endl;
	
	// Weighted mean example
	Vector<Real> values({85, 90, 78, 92});
	Vector<Real> weights({2, 3, 1, 4});  // credits/importance
	
	std::cout << "\n--- Weighted Statistics (GPA example) ---\n";
	std::cout << "Grades: 85, 90, 78, 92\n";
	std::cout << "Credits: 2, 3, 1, 4\n";
	std::cout << "Simple average:   " << Statistics::Mean(values) << std::endl;
	std::cout << "Weighted average: " << Statistics::WeightedMean(values, weights) << std::endl;
	std::cout << "Weighted std dev: " << Statistics::WeightedStdDev(values, weights) << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                    CORRELATION AND COVARIANCE                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Statistics_Correlation()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Correlation and Covariance\n";
	std::cout << "==========================================================================\n";
	
	// Positively correlated data (study hours vs test score)
	Vector<Real> studyHours({2, 3, 4, 5, 6, 7, 8, 9});
	Vector<Real> testScore({65, 70, 72, 80, 85, 88, 92, 95});
	
	std::cout << "\n--- Positive Correlation (Study Hours vs Test Score) ---\n";
	std::cout << "Study Hours: ";
	for (int i = 0; i < studyHours.size(); i++) std::cout << studyHours[i] << " ";
	std::cout << "\nTest Scores: ";
	for (int i = 0; i < testScore.size(); i++) std::cout << testScore[i] << " ";
	std::cout << std::endl;
	
	Real cov = Statistics::Covariance(studyHours, testScore);
	Real r = Statistics::PearsonCorrelation(studyHours, testScore);
	Real r2 = Statistics::RSquared(studyHours, testScore);
	Vector<Real> frequencyWeights({1, 1, 2, 2, 3, 3, 1, 1});
	Real weightedCov = Statistics::WeightedCovariance(studyHours, testScore, frequencyWeights);
	Real weightedR = Statistics::WeightedPearsonCorrelation(studyHours, testScore, frequencyWeights);
	
	std::cout << "\nCovariance:    " << cov << std::endl;
	std::cout << "Correlation r: " << r << " (strong positive)" << std::endl;
	std::cout << "R-squared:     " << r2 << " (" << r2*100 << "% variance explained)" << std::endl;
	std::cout << "Weighted covariance: " << weightedCov << std::endl;
	std::cout << "Weighted correlation: " << weightedR << std::endl;
	
	// Significance test
	auto result = Statistics::PearsonCorrelationWithTest(studyHours, testScore);
	std::cout << "\nSignificance Test:\n";
	std::cout << "  t-statistic: " << result.tStatistic << std::endl;
	std::cout << "  df: " << result.degreesOfFreedom << std::endl;
	std::cout << "  (Compare |t| to t-critical for significance)\n";
}

void Docs_Demo_Statistics_CovMatrix()
{
	std::cout << "\n==========================================================================\n";
	std::cout << "Demo: Covariance Matrix\n";
	std::cout << "==========================================================================\n";
	
	// Multivariate data: 5 observations of 3 variables
	Matrix<Real> data(5, 3);
	// Variable 1: Height (cm)
	data(0,0) = 170; data(1,0) = 175; data(2,0) = 165; data(3,0) = 180; data(4,0) = 168;
	// Variable 2: Weight (kg)
	data(0,1) = 70;  data(1,1) = 78;  data(2,1) = 65;  data(3,1) = 85;  data(4,1) = 72;
	// Variable 3: Age
	data(0,2) = 25;  data(1,2) = 30;  data(2,2) = 22;  data(3,2) = 35;  data(4,2) = 28;
	
	std::cout << "\nMultivariate data (5 obs x 3 vars):\n";
	std::cout << "Height(cm)  Weight(kg)  Age\n";
	for (int i = 0; i < 5; i++) {
		std::cout << data(i,0) << "        " << data(i,1) << "          " << data(i,2) << std::endl;
	}
	
	Matrix<Real> covMat = Statistics::CovarianceMatrix(data);
	Matrix<Real> corrMat = Statistics::CorrelationMatrix(data);
	
	std::cout << "\n--- Covariance Matrix ---\n";
	std::cout << "           Height   Weight    Age\n";
	std::cout << "Height    " << covMat(0,0) << "     " << covMat(0,1) << "      " << covMat(0,2) << std::endl;
	std::cout << "Weight    " << covMat(1,0) << "     " << covMat(1,1) << "      " << covMat(1,2) << std::endl;
	std::cout << "Age       " << covMat(2,0) << "     " << covMat(2,1) << "       " << covMat(2,2) << std::endl;
	
	std::cout << "\nDiagonal = variances, off-diagonal = covariances\n";
	std::cout << "\n--- Correlation Matrix ---\n";
	std::cout << "           Height   Weight    Age\n";
	std::cout << "Height    " << corrMat(0,0) << "     " << corrMat(0,1) << "      " << corrMat(0,2) << std::endl;
	std::cout << "Weight    " << corrMat(1,0) << "     " << corrMat(1,1) << "      " << corrMat(1,2) << std::endl;
	std::cout << "Age       " << corrMat(2,0) << "     " << corrMat(2,1) << "       " << corrMat(2,2) << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         MAIN DEMO FUNCTION                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Statistics_Stability()
{
	const Real huge = std::pow(Real(10), std::numeric_limits<Real>::max_exponent10);
	const Real unit = std::pow(Real(10), std::min(100, std::numeric_limits<Real>::max_exponent10 / 3));
	Vector<Real> pair({huge, huge});
	Vector<Real> large({-unit, 0, unit});
	Vector<Real> small({-Real(1) / unit, 0, Real(1) / unit});
	Real mean = Statistics::Avg(pair);
	Real median = Statistics::Median(pair);
	Real largeCorrelation = Statistics::PearsonCorrelation(large, large);
	Real smallCorrelation = Statistics::PearsonCorrelation(small, small);
	std::cout << "\n--- Unit-stable summaries (double: 1e308, 1e100, 1e-100) ---\n"
	          << "Mean: " << mean << ", median: " << median
	          << ", correlations: " << largeCorrelation << ", " << smallCorrelation << '\n';
	if (!std::isfinite(mean) || !std::isfinite(median)
	    || !IsWithinRel(mean, huge, Precision::EqualityTolerance)
	    || !IsWithinRel(median, huge, Precision::EqualityTolerance)
	    || !IsWithinRel(largeCorrelation, Real(1), Precision::EqualityTolerance)
	    || !IsWithinRel(smallCorrelation, Real(1), Precision::EqualityTolerance))
		throw std::runtime_error("Statistical stability demo failed");

	Vector<Real> extreme({-huge, 0, huge});
	Real stddev = Statistics::SampleStdDev(extreme);
	bool varianceRejected = false;
	try {
		Statistics::SampleVariance(extreme);
	}
	catch (const StatisticsError&) {
		varianceRejected = true;
	}
	if (!varianceRejected || !std::isfinite(stddev) || !IsWithinRel(stddev, huge, Precision::EqualityTolerance))
		throw std::runtime_error("Statistical output range demo failed");
	std::cout << "Finite sample stddev: " << stddev << "; overflowing variance rejected\n";
}

void Docs_Demo_Statistics_ProbabilityTails()
{
	Real lower = Statistics::StandardNormal::cdf(-10);
	Real upper = Statistics::StandardNormal::sf(10);
	Real quantile = Statistics::InverseStandardNormalCdf(Real(1e-20));
	Real poisson = Statistics::PoissonDistribution(40).cdf(0);
	Real gammaUpper = SpecialFunctions::RegularizedGammaQ(1, 40);
	Real inverseErfc = SpecialFunctions::ErfcInv(Real(1e-20));
	std::cout << "\n--- Representable probability tails ---\n"
	          << "Normal lower/upper: " << lower << ", " << upper
	          << "; probit(1e-20): " << quantile << '\n'
	          << "Poisson CDF(0): " << poisson << "; GammaQ(1,40): " << gammaUpper
	          << "; ErfcInv(1e-20): " << inverseErfc << '\n';
	Real normalReference = Real(0.5) * std::erfc(Real(10) / std::sqrt(Real(2)));
	Real poissonReference = std::exp(Real(-40));
	if (!(lower > 0 && upper > 0 && poisson > 0 && gammaUpper > 0)
	    || !IsWithinRel(lower, normalReference, Precision::EqualityTolerance)
	    || !IsWithinRel(upper, normalReference, Precision::EqualityTolerance)
	    || !IsWithinRel(poisson, poissonReference, Precision::EqualityTolerance)
	    || !IsWithinRel(gammaUpper, poissonReference, Precision::EqualityTolerance)
	    || !std::isfinite(quantile) || !std::isfinite(inverseErfc)
	    || !IsWithinRel(quantile, Real(-9.2623400897984075737), Precision::EqualityTolerance)
	    || !IsWithinRel(inverseErfc, Real(6.6015806223551425615), Precision::EqualityTolerance))
		throw std::runtime_error("Probability tail demo failed");
}

void Docs_Demo_Statistics_Quantiles()
{
	Real fCritical = Statistics::FDistribution(1, 1).criticalValue(Real(0.05));
	Real betaQuantile = Statistics::BetaDistribution(Real(0.5), 1).inverseCdf(Real(1e-6));
	Real tQuantile = Statistics::TDistribution(2).inverseCdf(Real(0.999999));
	Real chiCritical = Statistics::ChiSquareDistribution(2).criticalValue(Real(1e-20));
	Real gammaQuantile = Statistics::GammaDistribution(1, 2).inverseCdf(Real(0.9));
	std::cout << "\n--- Safeguarded continuous quantiles ---\n"
	          << "F critical: " << fCritical << "; Beta quantile: " << betaQuantile
	          << "; Student-t quantile: " << tQuantile << '\n'
	          << "Chi-square critical: " << chiCritical << "; Gamma quantile: " << gammaQuantile << '\n';
	Real p = Real(0.999999);
	if (!IsWithinRel(fCritical, 1 / std::pow(std::tan(Constants::PI * Real(0.05) / 2), 2), Precision::EqualityTolerance)
	    || !IsWithinRel(betaQuantile, Real(1e-12), Precision::EqualityTolerance)
	    || !IsWithinRel(tQuantile, (2 * p - 1) / std::sqrt(2 * p * (1 - p)), Precision::EqualityTolerance)
	    || !IsWithinRel(chiCritical, -2 * std::log(Real(1e-20)), Precision::EqualityTolerance)
	    || !IsWithinRel(gammaQuantile, -2 * std::log1p(-Real(0.9)), Precision::EqualityTolerance))
		throw std::runtime_error("Continuous quantile demo failed");
}

void Docs_Demo_Statistics_DiscreteBoundaries()
{
	Statistics::BinomialDistribution binomial(20, Real(0.4));
	Real probability = binomial.pmf(8);
	Real cumulative = binomial.cdf(8);
	int quantile = binomial.inverseCdf(Real(0.95));
	Statistics::PoissonDistribution poisson(3);
	Statistics::HypergeometricDistribution hypergeometric(52, 13, 5);
	int geometricLower = Statistics::GeometricDistribution(Real(0.5)).inverseCdf(Real(1e-20));
	int certainGeometric = Statistics::GeometricDistribution(1).inverseCdf(1);
	int certainNegative = Statistics::NegativeBinomialDistribution(3, 1).inverseCdf(1);
	int zeroBinomial = Statistics::BinomialDistribution(5, 0).inverseCdf(1);
	int poissonQuantile = Statistics::PoissonDistribution(40).inverseCdf(Real(0.99));
	std::cout << "\n--- Discrete quantile support and checked boundaries ---\n"
	          << "Binomial PMF/CDF: " << probability << ", " << cumulative << "; 95% quantile: " << quantile << '\n'
	          << "Geometric lower/degenerate: " << geometricLower << ", " << certainGeometric
	          << "; negative-binomial degenerate: " << certainNegative << "; binomial zero mass: " << zeroBinomial << '\n'
	          << "Poisson(40) 99% quantile: " << poissonQuantile << '\n';
	if (geometricLower != 0 || certainGeometric != 0 || certainNegative != 0 || zeroBinomial != 0 || poissonQuantile != 55
	    || !(probability > 0 && cumulative > 0 && cumulative <= 1)
	    || binomial.cdf(quantile) < Real(0.95) || binomial.cdf(quantile - 1) >= Real(0.95)
	    || poisson.cdf(0) <= 0 || hypergeometric.pmf(2) <= 0)
		throw std::runtime_error("Discrete quantile support demo failed");
	bool rejected = false;
	try {
		Statistics::GeometricDistribution(Real(1e-20)).inverseCdf(Real(0.5));
	} catch (const ConvergenceError&) {
		rejected = true;
	}
	if (!rejected)
		throw std::runtime_error("Unrepresentable discrete quantile was not rejected");
}

void Docs_Demo_Statistics_ParameterValidation()
{
	Statistics::GammaDistribution gamma(2, 1);
	gamma.shape = std::numeric_limits<Real>::quiet_NaN();
	std::mt19937 gen(17019);
	const auto before = gen;
	bool rejected = false;
	try {
		gamma.sample(gen);
	} catch (const StatisticsError&) {
		rejected = true;
	}
	if (!rejected || gen != before)
		throw std::runtime_error("Invalid gamma state was not rejected before RNG use");
	std::cout << "\n--- Distribution parameter validation ---\n"
	          << "Invalid mutable gamma state rejected without advancing the generator.\n";
}

void Docs_Demo_Statistics_DiscreteSampling()
{
	std::mt19937 gen(34023);
	int largeCount = Statistics::PoissonDistribution(1000000).sample(gen);
	const auto before = gen;
	bool rejected = false;
	try {
		Statistics::PoissonDistribution(Real(std::numeric_limits<int>::max()) + 1).sample(gen);
	} catch (const NumericalMethodError&) {
		rejected = true;
	}
	if (largeCount < 0 || !rejected || gen != before)
		throw std::runtime_error("Checked Poisson sampling demo failed");
	std::cout << "\n--- Checked discrete sampling ---\n"
	          << "Poisson(1000000) sample: " << largeCount
	          << "; unsupported sampling rate rejected without advancing the generator.\n";
}

void Docs_Demo_Statistics()
{
	std::cout << "\n##########################################################################\n";
	std::cout << "#                    STATISTICS DEMOS                                     #\n";
	std::cout << "##########################################################################\n";
	
	Docs_Demo_Statistics_Basic();
	Docs_Demo_Statistics_Moments();
	Docs_Demo_Statistics_Order();
	Docs_Demo_Statistics_Robust();
	Docs_Demo_Statistics_Means();
	Docs_Demo_Statistics_Correlation();
	Docs_Demo_Statistics_CovMatrix();
	Docs_Demo_Statistics_Stability();
	Docs_Demo_Statistics_ProbabilityTails();
	Docs_Demo_Statistics_Quantiles();
	Docs_Demo_Statistics_DiscreteBoundaries();
	Docs_Demo_Statistics_ParameterValidation();
	Docs_Demo_Statistics_DiscreteSampling();
	
	std::cout << "\n=== All Statistics Demos Complete ===\n";
}
