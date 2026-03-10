///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ChebyshevPolynom.h                                                  ///
///  Description: Chebyshev polynomial evaluation (first and second kind)             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CHEBYSHEV_POLYNOM_H
#define MML_CHEBYSHEV_POLYNOM_H

#include "../MMLBase.h"
#include "../MMLExceptions.h"

namespace MML
{
	///////////////////////////////////////////////////////////////////////////////////////////
	///                         CHEBYSHEV POLYNOMIAL EVALUATION                             ///
	///////////////////////////////////////////////////////////////////////////////////////////

	/// @brief Evaluate Chebyshev polynomial of the first kind T_n(x)
	/// @details Uses the recurrence relation:
	///          - T_0(x) = 1
	///          - T_1(x) = x
	///          - T_{n+1}(x) = 2x · T_n(x) - T_{n-1}(x)
	///
	///          Key properties:
	///          - T_n(1) = 1 for all n
	///          - T_n(-1) = (-1)^n
	///          - T_n(cos θ) = cos(nθ)  (trigonometric definition)
	///          - |T_n(x)| ≤ 1 for x ∈ [-1, 1]
	/// @param n Polynomial degree (must be non-negative)
	/// @param x Evaluation point
	/// @return Value of T_n(x)
	/// @throws ArgumentError if n < 0
	inline Real ChebyshevT(int n, Real x)
	{
		if (n < 0)
			throw ArgumentError("ChebyshevT: n must be non-negative");
		
		if (n == 0) return 1.0;
		if (n == 1) return x;
		
		Real T_prev2 = 1.0;   // T_0
		Real T_prev1 = x;     // T_1
		Real T_curr = 0.0;
		
		for (int k = 2; k <= n; k++)
		{
			T_curr = 2.0 * x * T_prev1 - T_prev2;
			T_prev2 = T_prev1;
			T_prev1 = T_curr;
		}
		
		return T_curr;
	}

	/// @brief Evaluate Chebyshev polynomial of the second kind U_n(x)
	/// @details Uses the recurrence relation:
	///          - U_0(x) = 1
	///          - U_1(x) = 2x
	///          - U_{n+1}(x) = 2x · U_n(x) - U_{n-1}(x)
	///
	///          Key properties:
	///          - U_n(1) = n + 1
	///          - U_n(-1) = (-1)^n · (n + 1)
	///          - U_n(cos θ) = sin((n+1)θ) / sin θ
	///
	///          Relationship to T_n:
	///          - dT_n/dx = n · U_{n-1}(x)
	/// @param n Polynomial degree (must be non-negative)
	/// @param x Evaluation point
	/// @return Value of U_n(x)
	/// @throws ArgumentError if n < 0
	inline Real ChebyshevU(int n, Real x)
	{
		if (n < 0)
			throw ArgumentError("ChebyshevU: n must be non-negative");
		
		if (n == 0) return 1.0;
		if (n == 1) return 2.0 * x;
		
		Real U_prev2 = 1.0;       // U_0
		Real U_prev1 = 2.0 * x;   // U_1
		Real U_curr = 0.0;
		
		for (int k = 2; k <= n; k++)
		{
			U_curr = 2.0 * x * U_prev1 - U_prev2;
			U_prev2 = U_prev1;
			U_prev1 = U_curr;
		}
		
		return U_curr;
	}

}   // namespace MML

#endif // MML_CHEBYSHEV_POLYNOM_H
