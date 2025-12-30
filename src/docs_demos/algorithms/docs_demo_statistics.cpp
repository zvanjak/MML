#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Vector.h"
#include "algorithms/Statistics.h"
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
	std::cout << "Variance:           " << var << std::endl;
	std::cout << "Standard Deviation: " << stdDev << std::endl;
	
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
	
	std::cout << "\n--- Weighted Mean (GPA example) ---\n";
	std::cout << "Grades: 85, 90, 78, 92\n";
	std::cout << "Credits: 2, 3, 1, 4\n";
	std::cout << "Simple average:   " << Statistics::Mean(values) << std::endl;
	std::cout << "Weighted average: " << Statistics::WeightedMean(values, weights) << std::endl;
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
	
	std::cout << "\nCovariance:    " << cov << std::endl;
	std::cout << "Correlation r: " << r << " (strong positive)" << std::endl;
	std::cout << "R-squared:     " << r2 << " (" << r2*100 << "% variance explained)" << std::endl;
	
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
	
	std::cout << "\n--- Covariance Matrix ---\n";
	std::cout << "           Height   Weight    Age\n";
	std::cout << "Height    " << covMat(0,0) << "     " << covMat(0,1) << "      " << covMat(0,2) << std::endl;
	std::cout << "Weight    " << covMat(1,0) << "     " << covMat(1,1) << "      " << covMat(1,2) << std::endl;
	std::cout << "Age       " << covMat(2,0) << "     " << covMat(2,1) << "       " << covMat(2,2) << std::endl;
	
	std::cout << "\nDiagonal = variances, off-diagonal = covariances\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                         MAIN DEMO FUNCTION                                          ///
///////////////////////////////////////////////////////////////////////////////////////////

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
	
	std::cout << "\n=== All Statistics Demos Complete ===\n";
}
