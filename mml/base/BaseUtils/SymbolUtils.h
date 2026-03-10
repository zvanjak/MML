///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        SymbolUtils.h                                                       ///
///  Description: Mathematical symbol utilities (Levi-Civita, Kronecker delta, etc.)  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_SYMBOL_UTILS_H
#define MML_SYMBOL_UTILS_H

namespace MML
{
	namespace Utils
	{
		/// @brief Computes the 3D Levi-Civita symbol (permutation symbol).
		/// @details Returns +1 for even permutations of (1,2,3), -1 for odd permutations,
		///          and 0 if any indices are equal or out of range.
		/// @param i First index (1-3)
		/// @param j Second index (1-3)
		/// @param k Third index (1-3)
		/// @return +1, -1, or 0 depending on permutation parity
		static constexpr int LeviCivita(int i, int j, int k) noexcept
		{
			if (i == j || j == k || i == k)
				return 0;

			if (i == 1 && j == 2 && k == 3)
				return 1;
			if (i == 2 && j == 3 && k == 1)
				return 1;
			if (i == 3 && j == 1 && k == 2)
				return 1;
			if (i == 3 && j == 2 && k == 1)
				return -1;
			if (i == 2 && j == 1 && k == 3)
				return -1;
			if (i == 1 && j == 3 && k == 2)
				return -1;

			return 0;
		}

		/// @brief Computes the 4D Levi-Civita symbol (permutation symbol).
		/// @details Returns +1 for even permutations of (1,2,3,4), -1 for odd permutations,
		///          and 0 if any indices are repeated or out of range.
		/// @param i First index (1-4)
		/// @param j Second index (1-4)
		/// @param k Third index (1-4)
		/// @param l Fourth index (1-4)
		/// @return +1, -1, or 0 depending on permutation parity
		static constexpr int LeviCivita(int i, int j, int k, int l) noexcept
		{
			// Validate input: must be permutation of {1, 2, 3, 4}
			bool present[5] = { false, false, false, false, false }; // index 0 unused
			int indices[4] = { i, j, k, l };
			
			for (int idx = 0; idx < 4; idx++)
			{
				if (indices[idx] < 1 || indices[idx] > 4)
					return 0; // Invalid index
				if (present[indices[idx]])
					return 0; // Duplicate index
				present[indices[idx]] = true;
			}
			
			// All indices 1,2,3,4 present - compute sign from permutation
			int a[4] = { i, j, k, l };
			int ret = 1;

			for (int i = 0; i < 4; i++)
				for (int j = i + 1; j < 4; j++)
					if (a[i] > a[j])
					{
						int tmp = a[i];
						a[i] = a[j];
						a[j] = tmp;
						ret *= -1;
					}
			return ret;
		}

		/// @brief Kronecker delta function.
		/// @details Returns 1 if indices are equal, 0 otherwise.
		/// @param i First index
		/// @param j Second index
		/// @return 1 if i == j, 0 otherwise
		static constexpr int KroneckerDelta(int i, int j) noexcept
		{
			return (i == j) ? 1 : 0;
		}

	} // namespace Utils
} // namespace MML

#endif // MML_SYMBOL_UTILS_H
