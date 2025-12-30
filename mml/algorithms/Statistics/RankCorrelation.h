///////////////////////////////////////////////////////////////////////////////////////////
// RankCorrelation.h
// 
// Rank-based correlation methods for MML (MinimalMathLibrary)
// Extracted from Statistics.h as part of refactoring to improve organization
//
// Contains: Spearman's rho, Kendall's tau, and helper functions
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_RANK_CORRELATION_H
#define MML_RANK_CORRELATION_H

#include "../Statistics.h"

namespace MML
{
	namespace Statistics
	{
		/**
		 * @brief Compute ranks for a dataset
		 * 
		 * Converts raw data values to ranks (1-based). Tied values receive
		 * the average of the ranks they would have received.
		 * 
		 * Example: [5, 2, 2, 8] → [3, 1.5, 1.5, 4]
		 * 
		 * @param data Vector of values to rank
		 * @return Vector of ranks (1-based, with average ranks for ties)
		 * 
		 * @note Uses average ranking for ties: if values at positions 2,3,4
		 *       are tied, they all receive rank (2+3+4)/3 = 3
		 * 
		 * Complexity: O(n log n) for sorting
		 */
		static Vector<Real> ComputeRanks(const Vector<Real>& data)
		{
			int n = data.size();
			if (n == 0)
				return Vector<Real>();

			// Create index array and sort by values
			std::vector<int> indices(n);
			for (int i = 0; i < n; i++)
				indices[i] = i;

			std::sort(indices.begin(), indices.end(), [&data](int a, int b) {
				return data[a] < data[b];
			});

			// Assign ranks with tie handling (average ranks for ties)
			Vector<Real> ranks(n);
			int i = 0;
			while (i < n) {
				int j = i;
				// Find extent of tied values
				while (j < n - 1 && std::abs(data[indices[j + 1]] - data[indices[i]]) < 1e-14)
					j++;

				// Average rank for tied values
				Real avgRank = (i + 1 + j + 1) / 2.0;  // 1-based ranks
				for (int k = i; k <= j; k++)
					ranks[indices[k]] = avgRank;

				i = j + 1;
			}

			return ranks;
		}

		/**
		 * @brief Compute Spearman's rank correlation coefficient
		 * 
		 * ρ (rho) = 1 - 6·Σd²ᵢ / (n(n²-1))  [no ties formula]
		 * 
		 * For data with ties, computes Pearson correlation on ranks.
		 * 
		 * Spearman correlation measures monotonic relationship (not just linear).
		 * Values range from -1 (perfect negative monotonic) to +1 (perfect positive).
		 * 
		 * @param x First variable (vector of values)
		 * @param y Second variable (must have same size as x)
		 * @return Spearman correlation coefficient ρ in [-1, 1]
		 * @throws StatisticsError if vectors are empty, have different sizes,
		 *         or have fewer than 2 elements
		 * 
		 * Complexity: O(n log n) due to ranking
		 */
		static Real SpearmanCorrelation(const Vector<Real>& x, const Vector<Real>& y)
		{
			int n = x.size();
			if (n < 2)
				throw StatisticsError("Vector size must be at least 2 in SpearmanCorrelation");
			if (y.size() != n)
				throw StatisticsError("Vectors must have the same size in SpearmanCorrelation");

			// Compute ranks
			Vector<Real> rankX = ComputeRanks(x);
			Vector<Real> rankY = ComputeRanks(y);

			// Use Pearson correlation on ranks (handles ties correctly)
			return PearsonCorrelation(rankX, rankY);
		}

		/**
		 * @brief Result structure for rank correlation with significance test
		 */
		struct RankCorrelationResult
		{
			Real rho;              ///< Correlation coefficient (Spearman or Kendall)
			Real zScore;           ///< z-score for large sample approximation
			int n;                 ///< Sample size

			/**
			 * @brief Check if correlation is significant at given alpha level
			 * 
			 * Uses normal approximation (valid for n > 10)
			 * Critical z-values: α=0.05 → z=1.96, α=0.01 → z=2.576
			 */
			bool IsSignificant(Real alpha = 0.05) const
			{
				Real criticalZ = (alpha <= 0.01) ? 2.576 : 1.96;
				return std::abs(zScore) > criticalZ;
			}
		};

		/**
		 * @brief Compute Spearman correlation with significance test
		 * 
		 * Returns correlation coefficient with z-score for significance testing.
		 * Uses large-sample approximation: z = ρ·√(n-1)
		 * 
		 * @param x First variable
		 * @param y Second variable
		 * @return RankCorrelationResult with rho, zScore, and n
		 * @throws StatisticsError if vectors are invalid
		 */
		static RankCorrelationResult SpearmanCorrelationWithTest(const Vector<Real>& x, const Vector<Real>& y)
		{
			int n = x.size();
			if (n < 3)
				throw StatisticsError("Need at least 3 observations for significance test in SpearmanCorrelationWithTest");

			Real rho = SpearmanCorrelation(x, y);

			// Large sample approximation: z = rho * sqrt(n - 1)
			Real zScore = rho * std::sqrt(static_cast<Real>(n - 1));

			return RankCorrelationResult{rho, zScore, n};
		}

		/**
		 * @brief Compute Kendall's tau-b rank correlation coefficient
		 * 
		 * τ_b = (C - D) / √((C + D + Tx)(C + D + Ty))
		 * 
		 * where:
		 * - C = number of concordant pairs
		 * - D = number of discordant pairs
		 * - Tx = pairs tied only in x
		 * - Ty = pairs tied only in y
		 * 
		 * Kendall's tau is more robust than Spearman for small samples
		 * and has a more intuitive interpretation (probability difference).
		 * 
		 * @param x First variable (vector of values)
		 * @param y Second variable (must have same size as x)
		 * @return Kendall tau-b coefficient in [-1, 1]
		 * @throws StatisticsError if vectors are empty, have different sizes,
		 *         or have fewer than 2 elements
		 * 
		 * Complexity: O(n²) - naive implementation
		 */
		static Real KendallCorrelation(const Vector<Real>& x, const Vector<Real>& y)
		{
			int n = x.size();
			if (n < 2)
				throw StatisticsError("Vector size must be at least 2 in KendallCorrelation");
			if (y.size() != n)
				throw StatisticsError("Vectors must have the same size in KendallCorrelation");

			long long concordant = 0;
			long long discordant = 0;
			long long tiedX = 0;
			long long tiedY = 0;

			const Real eps = 1e-14;

			// Count all pairs
			for (int i = 0; i < n - 1; i++) {
				for (int j = i + 1; j < n; j++) {
					Real dx = x[j] - x[i];
					Real dy = y[j] - y[i];

					bool xTied = std::abs(dx) < eps;
					bool yTied = std::abs(dy) < eps;

					if (xTied && yTied) {
						// Both tied - don't count
						continue;
					} else if (xTied) {
						tiedX++;
					} else if (yTied) {
						tiedY++;
					} else if ((dx > 0 && dy > 0) || (dx < 0 && dy < 0)) {
						concordant++;
					} else {
						discordant++;
					}
				}
			}

			// tau-b formula with tie correction
			Real numerator = static_cast<Real>(concordant - discordant);
			Real denom1 = static_cast<Real>(concordant + discordant + tiedX);
			Real denom2 = static_cast<Real>(concordant + discordant + tiedY);

			if (denom1 <= 0.0 || denom2 <= 0.0)
				return 0.0;  // All pairs tied

			return numerator / std::sqrt(denom1 * denom2);
		}

		/**
		 * @brief Compute Kendall correlation with significance test
		 * 
		 * Returns correlation coefficient with z-score for significance testing.
		 * Uses asymptotic normal approximation:
		 * z = τ / √(2(2n+5) / 9n(n-1))
		 * 
		 * @param x First variable
		 * @param y Second variable
		 * @return RankCorrelationResult with tau (as rho), zScore, and n
		 * @throws StatisticsError if vectors are invalid
		 */
		static RankCorrelationResult KendallCorrelationWithTest(const Vector<Real>& x, const Vector<Real>& y)
		{
			int n = x.size();
			if (n < 3)
				throw StatisticsError("Need at least 3 observations for significance test in KendallCorrelationWithTest");

			Real tau = KendallCorrelation(x, y);

			// Standard error for Kendall's tau (no ties):
			// SE = sqrt(2(2n+5) / 9n(n-1))
			Real variance = (2.0 * (2.0 * n + 5.0)) / (9.0 * n * (n - 1));
			Real zScore = tau / std::sqrt(variance);

			return RankCorrelationResult{tau, zScore, n};
		}
	}
}

#endif // MML_RANK_CORRELATION_H
