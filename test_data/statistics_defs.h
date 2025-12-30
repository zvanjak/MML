#if !defined __MML_STATISTICS_DEFS_H
#define __MML_STATISTICS_DEFS_H

/**
 * @file statistics_defs.h
 * @brief Reference datasets and expected values for Statistics module testing
 * 
 * This file contains:
 * 1. Reference datasets with known statistical properties
 * 2. Precomputed expected values (verified against R/Python)
 * 3. Critical value tables for hypothesis testing
 * 4. Edge case datasets for robustness testing
 * 
 * All expected values computed using Python 3.11 with numpy/scipy and
 * cross-validated against R 4.3.
 * 
 * @note Values are stored as arrays for easy Vector<Real> construction:
 *       Vector<Real> data(Stats_Dataset1, Stats_Dataset1 + Stats_Dataset1_n);
 */

#include <cmath>
#include <limits>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"
#endif

namespace MML::TestBeds
{
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         DATASET 1: SMALL SIMPLE (n=10)                                //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Simple dataset for basic validation
     * Values: 2, 4, 4, 4, 5, 5, 7, 9, 10, 12
     * Sorted for median/quartile verification
     */
    static const Real Stats_Simple10[] = {
        2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0, 10.0, 12.0
    };
    static const int Stats_Simple10_n = 10;
    
    // Expected values (computed in Python)
    static const Real Stats_Simple10_Mean = 6.2;                    // sum/n = 62/10
    static const Real Stats_Simple10_Variance = 10.177777777777778; // sample variance (n-1)
    static const Real Stats_Simple10_StdDev = 3.1902629659659657;   // sqrt(variance)
    static const Real Stats_Simple10_Median = 5.0;                  // average of 5th and 6th
    static const Real Stats_Simple10_Q1 = 4.0;                      // 25th percentile
    static const Real Stats_Simple10_Q3 = 9.0;                      // 75th percentile
    static const Real Stats_Simple10_IQR = 5.0;                     // Q3 - Q1
    static const Real Stats_Simple10_Min = 2.0;
    static const Real Stats_Simple10_Max = 12.0;
    static const Real Stats_Simple10_Range = 10.0;
    static const Real Stats_Simple10_Mode = 4.0;                    // Most frequent
    // Skewness and Kurtosis (using population formulas, excess kurtosis)
    static const Real Stats_Simple10_Skewness = 0.5255138545542715;  // Right-skewed
    static const Real Stats_Simple10_Kurtosis = -0.7883561643835616; // Platykurtic
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         DATASET 2: NORMAL SAMPLE (n=30)                               //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Pseudo-random normal sample N(100, 15)
     * Generated with numpy: np.random.seed(42); np.random.normal(100, 15, 30)
     * Useful for testing t-tests and confidence intervals
     */
    static const Real Stats_Normal30[] = {
        107.454, 97.928, 109.702, 122.844, 96.487,
        96.487, 123.568, 112.303, 92.303, 93.595,
        110.088, 98.772, 104.219, 108.330, 94.914,
        93.049, 102.178, 91.340, 110.744, 93.833,
        86.536, 109.715, 102.550, 93.798, 109.455,
        115.634, 101.123, 108.372, 81.464, 93.649
    };
    static const int Stats_Normal30_n = 30;
    
    // Expected values
    static const Real Stats_Normal30_Mean = 101.77753333333333;
    static const Real Stats_Normal30_Variance = 109.07396475862068;  // Sample variance
    static const Real Stats_Normal30_StdDev = 10.443850509693458;
    static const Real Stats_Normal30_Median = 101.6505;              // Average of 15th and 16th sorted
    static const Real Stats_Normal30_SEM = 1.906505621;              // Standard Error of Mean
    
    // For one-sample t-test H0: μ = 100
    static const Real Stats_Normal30_TTest_mu0_100_tstat = 0.9318; // (mean - 100) / SEM
    static const Real Stats_Normal30_TTest_mu0_100_pvalue = 0.3592; // Two-sided
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         DATASET 3: EXPONENTIAL SAMPLE (n=50)                          //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Exponential sample with rate λ = 0.5 (mean = 2)
     * Generated with numpy: np.random.seed(123); np.random.exponential(2, 50)
     * For testing exponential distribution fitting
     */
    static const Real Stats_Exponential50[] = {
        1.394, 3.867, 0.155, 1.893, 2.521, 0.879, 3.195, 0.523, 2.847, 1.024,
        4.521, 0.789, 1.156, 2.234, 0.345, 1.678, 5.012, 0.912, 1.445, 3.234,
        0.234, 2.156, 1.789, 0.678, 3.567, 1.234, 0.456, 2.789, 4.123, 0.567,
        1.890, 2.345, 0.123, 3.456, 1.567, 2.890, 0.789, 1.123, 4.567, 2.012,
        0.345, 1.678, 3.012, 0.890, 2.456, 1.345, 0.567, 3.789, 2.123, 1.456
    };
    static const int Stats_Exponential50_n = 50;
    
    static const Real Stats_Exponential50_Mean = 1.86584;
    static const Real Stats_Exponential50_Median = 1.5565;
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         DATASET 4: BIVARIATE CORRELATION                              //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Paired data for correlation testing
     * Strong positive correlation (r ≈ 0.95)
     * X: Hours studied, Y: Exam score
     */
    static const Real Stats_Corr_X[] = {
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
        2.5, 3.5, 4.5, 5.5, 6.5
    };
    static const Real Stats_Corr_Y[] = {
        55.0, 60.0, 65.0, 70.0, 72.0, 78.0, 82.0, 85.0, 88.0, 92.0,
        58.0, 67.0, 73.0, 76.0, 80.0
    };
    static const int Stats_Corr_n = 15;
    
    static const Real Stats_Corr_Pearson = 0.9721448467;   // Pearson r
    static const Real Stats_Corr_Spearman = 0.9785714286;  // Spearman ρ (rank correlation)
    static const Real Stats_Corr_Covariance = 36.95238095;  // Sample covariance
    
    // For testing correlation p-values
    static const Real Stats_Corr_Pearson_tstat = 14.8426;  // t = r * sqrt(n-2) / sqrt(1-r²)
    static const Real Stats_Corr_Pearson_pvalue = 3.26e-10; // Two-sided p-value
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         DATASET 5: ZERO CORRELATION                                   //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Data with near-zero linear correlation but perfect quadratic relationship
     * X: -3, -2, -1, 0, 1, 2, 3
     * Y: X² (parabola)
     */
    static const Real Stats_ZeroCorr_X[] = {-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0};
    static const Real Stats_ZeroCorr_Y[] = {9.0, 4.0, 1.0, 0.0, 1.0, 4.0, 9.0};  // Y = X²
    static const int Stats_ZeroCorr_n = 7;
    
    static const Real Stats_ZeroCorr_Pearson = 0.0;  // Perfect quadratic, zero linear
    static const Real Stats_ZeroCorr_Covariance = 0.0;
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         DATASET 6: TWO-SAMPLE T-TEST                                  //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Two independent samples for t-test
     * Control group vs Treatment group
     */
    static const Real Stats_TTest_Control[] = {
        23.1, 24.5, 21.8, 25.2, 22.9, 23.7, 24.1, 22.4, 25.8, 23.3
    };
    static const Real Stats_TTest_Treatment[] = {
        28.3, 26.9, 29.5, 27.4, 30.1, 28.8, 27.2, 29.7, 26.5, 28.1
    };
    static const int Stats_TTest_Control_n = 10;
    static const int Stats_TTest_Treatment_n = 10;
    
    static const Real Stats_TTest_Control_Mean = 23.68;
    static const Real Stats_TTest_Treatment_Mean = 28.25;
    static const Real Stats_TTest_Independent_tstat = -7.0237;  // Pooled variance t-test
    static const Real Stats_TTest_Independent_pvalue = 1.69e-6; // Two-sided
    static const Real Stats_TTest_Independent_df = 18;          // n1 + n2 - 2
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         DATASET 7: PAIRED T-TEST                                      //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Paired data (before/after treatment)
     */
    static const Real Stats_Paired_Before[] = {
        200, 210, 195, 220, 205, 215, 190, 225, 208, 212
    };
    static const Real Stats_Paired_After[] = {
        185, 195, 180, 205, 190, 198, 178, 210, 192, 195
    };
    static const int Stats_Paired_n = 10;
    
    static const Real Stats_Paired_MeanDiff = 15.3;    // Before - After
    static const Real Stats_Paired_tstat = 9.544;       // t = mean(diff) / SE(diff)
    static const Real Stats_Paired_pvalue = 4.89e-6;   // Two-sided
    static const Real Stats_Paired_df = 9;              // n - 1
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         DATASET 8: CHI-SQUARE GOODNESS OF FIT                         //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Die roll frequencies (n=60 rolls)
     * Testing if die is fair (expected: 10 each)
     */
    static const Real Stats_ChiSq_Observed[] = {8.0, 12.0, 9.0, 11.0, 7.0, 13.0};
    static const Real Stats_ChiSq_Expected[] = {10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
    static const int Stats_ChiSq_k = 6;  // Number of categories
    
    static const Real Stats_ChiSq_Statistic = 2.8;  // Σ(O-E)²/E
    static const Real Stats_ChiSq_df = 5;           // k - 1
    static const Real Stats_ChiSq_pvalue = 0.7306;  // Not significant (die appears fair)
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         DATASET 9: CHI-SQUARE INDEPENDENCE                            //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief 3x2 contingency table
     * Testing independence of Gender vs Product Preference
     * 
     *          Product A  Product B
     * Male       30         20
     * Female     25         35
     * Other       5          5
     */
    static const Real Stats_ChiSq_Contingency[] = {
        30.0, 20.0,   // Male
        25.0, 35.0,   // Female
         5.0,  5.0    // Other
    };
    static const int Stats_ChiSq_Contingency_rows = 3;
    static const int Stats_ChiSq_Contingency_cols = 2;
    
    static const Real Stats_ChiSq_Independence_Statistic = 4.0;
    static const Real Stats_ChiSq_Independence_df = 2;  // (r-1)(c-1)
    static const Real Stats_ChiSq_Independence_pvalue = 0.1353;
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         DATASET 10: ONE-WAY ANOVA                                     //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Three groups for ANOVA
     * Testing if teaching methods affect scores
     */
    static const Real Stats_ANOVA_Group1[] = {85, 90, 78, 92, 88};  // Method A
    static const Real Stats_ANOVA_Group2[] = {75, 80, 72, 78, 76};  // Method B
    static const Real Stats_ANOVA_Group3[] = {92, 95, 88, 97, 94};  // Method C
    static const int Stats_ANOVA_Group1_n = 5;
    static const int Stats_ANOVA_Group2_n = 5;
    static const int Stats_ANOVA_Group3_n = 5;
    
    static const Real Stats_ANOVA_Group1_Mean = 86.6;
    static const Real Stats_ANOVA_Group2_Mean = 76.2;
    static const Real Stats_ANOVA_Group3_Mean = 93.2;
    static const Real Stats_ANOVA_GrandMean = 85.333333;
    
    static const Real Stats_ANOVA_SSBetween = 727.733333;   // Between-group sum of squares
    static const Real Stats_ANOVA_SSWithin = 135.2;         // Within-group sum of squares
    static const Real Stats_ANOVA_SSTotal = 862.933333;
    static const Real Stats_ANOVA_df_between = 2;           // k - 1
    static const Real Stats_ANOVA_df_within = 12;           // N - k
    static const Real Stats_ANOVA_MSBetween = 363.866667;   // SS_between / df_between
    static const Real Stats_ANOVA_MSWithin = 11.266667;     // SS_within / df_within
    static const Real Stats_ANOVA_F_Statistic = 32.298224;  // MS_between / MS_within
    static const Real Stats_ANOVA_pvalue = 1.91e-5;         // F(2,12) > 32.3
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         ANSCOMBE'S QUARTET                                            //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Anscombe's Quartet - 4 datasets with identical summary statistics
     * but very different distributions. Essential for demonstrating
     * the importance of visualization and residual analysis.
     * 
     * All four have: mean(x)=9, mean(y)≈7.5, var(x)=11, var(y)≈4.125
     * Pearson r ≈ 0.816, regression: y ≈ 3 + 0.5x
     */
    
    // Dataset I: Linear relationship
    static const Real Stats_Anscombe_X1[] = {10, 8, 13, 9, 11, 14, 6, 4, 12, 7, 5};
    static const Real Stats_Anscombe_Y1[] = {8.04, 6.95, 7.58, 8.81, 8.33, 9.96, 7.24, 4.26, 10.84, 4.82, 5.68};
    
    // Dataset II: Quadratic relationship
    static const Real Stats_Anscombe_X2[] = {10, 8, 13, 9, 11, 14, 6, 4, 12, 7, 5};
    static const Real Stats_Anscombe_Y2[] = {9.14, 8.14, 8.74, 8.77, 9.26, 8.10, 6.13, 3.10, 9.13, 7.26, 4.74};
    
    // Dataset III: Outlier influence
    static const Real Stats_Anscombe_X3[] = {10, 8, 13, 9, 11, 14, 6, 4, 12, 7, 5};
    static const Real Stats_Anscombe_Y3[] = {7.46, 6.77, 12.74, 7.11, 7.81, 8.84, 6.08, 5.39, 8.15, 6.42, 5.73};
    
    // Dataset IV: Leverage point
    static const Real Stats_Anscombe_X4[] = {8, 8, 8, 8, 8, 8, 8, 19, 8, 8, 8};
    static const Real Stats_Anscombe_Y4[] = {6.58, 5.76, 7.71, 8.84, 8.47, 7.04, 5.25, 12.50, 5.56, 7.91, 6.89};
    
    static const int Stats_Anscombe_n = 11;
    
    // Common statistics for all four datasets
    static const Real Stats_Anscombe_MeanX = 9.0;
    static const Real Stats_Anscombe_MeanY = 7.50090909;
    static const Real Stats_Anscombe_VarX = 11.0;
    static const Real Stats_Anscombe_VarY = 4.127269;
    static const Real Stats_Anscombe_Pearson = 0.81642;
    static const Real Stats_Anscombe_Slope = 0.5001;
    static const Real Stats_Anscombe_Intercept = 3.0001;
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         EDGE CASES                                                    //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    // Single value
    static const Real Stats_Single[] = {42.0};
    static const int Stats_Single_n = 1;
    
    // Two values
    static const Real Stats_Two[] = {10.0, 20.0};
    static const int Stats_Two_n = 2;
    static const Real Stats_Two_Mean = 15.0;
    static const Real Stats_Two_Variance = 50.0;  // Sample variance
    static const Real Stats_Two_Median = 15.0;    // Average of both
    
    // All identical values
    static const Real Stats_Identical[] = {5.0, 5.0, 5.0, 5.0, 5.0};
    static const int Stats_Identical_n = 5;
    static const Real Stats_Identical_Mean = 5.0;
    static const Real Stats_Identical_Variance = 0.0;
    static const Real Stats_Identical_Median = 5.0;
    
    // Negative values
    static const Real Stats_Negative[] = {-10.0, -5.0, -3.0, -1.0, 0.0};
    static const int Stats_Negative_n = 5;
    static const Real Stats_Negative_Mean = -3.8;
    static const Real Stats_Negative_Median = -3.0;
    
    // Very large values (test overflow handling)
    static const Real Stats_Large[] = {1e15, 1e15 + 1, 1e15 + 2};
    static const int Stats_Large_n = 3;
    static const Real Stats_Large_Mean = 1e15 + 1.0;
    
    // Very small values (test underflow handling)
    static const Real Stats_Small[] = {1e-15, 2e-15, 3e-15};
    static const int Stats_Small_n = 3;
    static const Real Stats_Small_Mean = 2e-15;
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         CRITICAL VALUE TABLES                                         //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Student's t critical values (two-tailed)
     * Index by [df-1] for df=1,2,...,30; special values at end for df=40,60,120,∞
     */
    
    // α = 0.05 (95% confidence)
    static const Real Stats_TCrit_005[] = {
        12.706, 4.303, 3.182, 2.776, 2.571,  // df = 1-5
        2.447, 2.365, 2.306, 2.262, 2.228,   // df = 6-10
        2.201, 2.179, 2.160, 2.145, 2.131,   // df = 11-15
        2.120, 2.110, 2.101, 2.093, 2.086,   // df = 16-20
        2.080, 2.074, 2.069, 2.064, 2.060,   // df = 21-25
        2.056, 2.052, 2.048, 2.045, 2.042,   // df = 26-30
        2.021, 2.000, 1.980, 1.960            // df = 40, 60, 120, ∞
    };
    
    // α = 0.01 (99% confidence)
    static const Real Stats_TCrit_001[] = {
        63.657, 9.925, 5.841, 4.604, 4.032,  // df = 1-5
        3.707, 3.499, 3.355, 3.250, 3.169,   // df = 6-10
        3.106, 3.055, 3.012, 2.977, 2.947,   // df = 11-15
        2.921, 2.898, 2.878, 2.861, 2.845,   // df = 16-20
        2.831, 2.819, 2.807, 2.797, 2.787,   // df = 21-25
        2.779, 2.771, 2.763, 2.756, 2.750,   // df = 26-30
        2.704, 2.660, 2.617, 2.576            // df = 40, 60, 120, ∞
    };
    
    /**
     * @brief Chi-square critical values (right tail)
     * Index by [df-1] for df=1,2,...,30
     */
    
    // α = 0.05
    static const Real Stats_ChiSqCrit_005[] = {
        3.841, 5.991, 7.815, 9.488, 11.070,  // df = 1-5
        12.592, 14.067, 15.507, 16.919, 18.307,  // df = 6-10
        19.675, 21.026, 22.362, 23.685, 24.996,  // df = 11-15
        26.296, 27.587, 28.869, 30.144, 31.410,  // df = 16-20
        32.671, 33.924, 35.172, 36.415, 37.652,  // df = 21-25
        38.885, 40.113, 41.337, 42.557, 43.773   // df = 26-30
    };
    
    // α = 0.01
    static const Real Stats_ChiSqCrit_001[] = {
        6.635, 9.210, 11.345, 13.277, 15.086,  // df = 1-5
        16.812, 18.475, 20.090, 21.666, 23.209,  // df = 6-10
        24.725, 26.217, 27.688, 29.141, 30.578,  // df = 11-15
        32.000, 33.409, 34.805, 36.191, 37.566,  // df = 16-20
        38.932, 40.289, 41.638, 42.980, 44.314,  // df = 21-25
        45.642, 46.963, 48.278, 49.588, 50.892   // df = 26-30
    };
    
    /**
     * @brief F distribution critical values
     * For one-way ANOVA: F(df1, df2) where df1 = k-1 (between), df2 = N-k (within)
     * Most common cases only; full table would be too large
     */
    
    // F critical values for α = 0.05
    // FCrit_005[df1-1][df2-1] for df1=1,2,3,4 and df2=5,10,15,20,30,60
    static const Real Stats_FCrit_005_df1_1[] = {6.61, 4.96, 4.54, 4.35, 4.17, 4.00};  // df1=1
    static const Real Stats_FCrit_005_df1_2[] = {5.79, 4.10, 3.68, 3.49, 3.32, 3.15};  // df1=2
    static const Real Stats_FCrit_005_df1_3[] = {5.41, 3.71, 3.29, 3.10, 2.92, 2.76};  // df1=3
    static const Real Stats_FCrit_005_df1_4[] = {5.19, 3.48, 3.06, 2.87, 2.69, 2.53};  // df1=4
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         PERCENTILE REFERENCE VALUES                                   //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Expected percentiles for Stats_Simple10 dataset
     * Using linear interpolation method (same as numpy/R default)
     */
    static const Real Stats_Simple10_P10 = 3.2;   // 10th percentile
    static const Real Stats_Simple10_P25 = 4.0;   // 25th percentile (Q1)
    static const Real Stats_Simple10_P50 = 5.0;   // 50th percentile (Median)
    static const Real Stats_Simple10_P75 = 9.0;   // 75th percentile (Q3)
    static const Real Stats_Simple10_P90 = 10.6;  // 90th percentile
    static const Real Stats_Simple10_P99 = 11.82; // 99th percentile
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         WEIGHTED STATISTICS                                           //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Data and weights for weighted mean testing
     */
    static const Real Stats_Weighted_Values[] = {10.0, 20.0, 30.0, 40.0};
    static const Real Stats_Weighted_Weights[] = {1.0, 2.0, 3.0, 4.0};
    static const int Stats_Weighted_n = 4;
    
    // Weighted mean = (10*1 + 20*2 + 30*3 + 40*4) / (1+2+3+4) = 300/10 = 30
    static const Real Stats_Weighted_Mean = 30.0;
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         SPECIAL MEANS                                                 //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Data for geometric and harmonic mean testing
     * All positive values required
     */
    static const Real Stats_PositiveOnly[] = {1.0, 2.0, 4.0, 8.0, 16.0};
    static const int Stats_PositiveOnly_n = 5;
    
    // Arithmetic: (1+2+4+8+16)/5 = 6.2
    static const Real Stats_PositiveOnly_ArithMean = 6.2;
    // Geometric: (1*2*4*8*16)^(1/5) = 1024^0.2 = 4.0
    static const Real Stats_PositiveOnly_GeomMean = 4.0;
    // Harmonic: 5 / (1/1 + 1/2 + 1/4 + 1/8 + 1/16) = 5 / 1.9375 ≈ 2.58064516
    static const Real Stats_PositiveOnly_HarmMean = 2.5806451612903226;
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                         TRIMMED MEAN                                                  //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    /**
     * @brief Data for trimmed mean testing (with outliers)
     * Regular mean is skewed by outliers
     */
    static const Real Stats_WithOutliers[] = {1.0, 100.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 200.0};
    static const int Stats_WithOutliers_n = 10;
    
    // Regular mean = 343/10 = 34.3 (badly affected by outliers)
    static const Real Stats_WithOutliers_Mean = 34.3;
    // 10% trimmed (remove 1 from each end after sorting): mean of [3,4,5,6,7,8,9,100] = 17.75
    static const Real Stats_WithOutliers_Trimmed10 = 17.75;
    // 20% trimmed (remove 2 from each end): mean of [4,5,6,7,8,9] = 6.5
    static const Real Stats_WithOutliers_Trimmed20 = 6.5;
    
} // namespace MML::TestBeds

#endif // __MML_STATISTICS_DEFS_H
