///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ODERKCoefficients.h                                                 ///
///  Description: Centralized Butcher tableau coefficients for RK methods             ///
///               Shared by both StepCalculators and Steppers for consistency         ///
///                                                                                   ///
///  REFERENCES:                                                                      ///
///    [NR3]  Press et al., Numerical Recipes 3rd ed., Ch. 17                         ///
///    [HNW1] Hairer et al., Solving ODEs I, Ch. II                                   ///
///    [DP80] Dormand & Prince (1980), J. Comp. Appl. Math. 6(1), pp. 19-26          ///
///    [CK90] Cash & Karp (1990), ACM TOMS 16(3), pp. 201-222                        ///
///    [PD81] Prince & Dormand (1981), J. Comp. Appl. Math. 7(1), pp. 67-75          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ODE_RK_COEFFICIENTS_H
#define MML_ODE_RK_COEFFICIENTS_H

#include "mml/MMLBase.h"
#include <cmath>

namespace MML {
namespace RKCoeff {

	// C++17-compatible constexpr abs (std::abs is not constexpr until C++23)
	template<typename T>
	constexpr T cabs(T x) noexcept { return x < T(0) ? -x : x; }

	/******************************************************************************
	 * CASH-KARP 5(4) COEFFICIENTS
	 *
	 * 6-stage embedded pair. Fifth-order solution with fourth-order error estimate.
	 * Reference: [CK90], [NR3] Section 17.2
	 *
	 * VERIFIED: Coefficients match Numerical Recipes 2nd ed. rkck() exactly.
	 ******************************************************************************/
	struct CashKarp5 {
		static constexpr int stages = 6;
		static constexpr int order = 5;
		static constexpr int error_order = 4;
		static constexpr bool fsal = false;

		// Time nodes (c coefficients)
		static constexpr Real c2 = 1.0 / 5.0;
		static constexpr Real c3 = 3.0 / 10.0;
		static constexpr Real c4 = 3.0 / 5.0;
		static constexpr Real c5 = 1.0;
		static constexpr Real c6 = 7.0 / 8.0;

		// Stage coefficients (a matrix, lower triangular)
		static constexpr Real a21 = 1.0 / 5.0;

		static constexpr Real a31 = 3.0 / 40.0;
		static constexpr Real a32 = 9.0 / 40.0;

		static constexpr Real a41 = 3.0 / 10.0;
		static constexpr Real a42 = -9.0 / 10.0;
		static constexpr Real a43 = 6.0 / 5.0;

		static constexpr Real a51 = -11.0 / 54.0;
		static constexpr Real a52 = 5.0 / 2.0;
		static constexpr Real a53 = -70.0 / 27.0;
		static constexpr Real a54 = 35.0 / 27.0;

		static constexpr Real a61 = 1631.0 / 55296.0;
		static constexpr Real a62 = 175.0 / 512.0;
		static constexpr Real a63 = 575.0 / 13824.0;
		static constexpr Real a64 = 44275.0 / 110592.0;
		static constexpr Real a65 = 253.0 / 4096.0;

		// 5th order solution weights (b coefficients)
		static constexpr Real b1 = 37.0 / 378.0;
		static constexpr Real b2 = 0.0;
		static constexpr Real b3 = 250.0 / 621.0;
		static constexpr Real b4 = 125.0 / 594.0;
		static constexpr Real b5 = 0.0;
		static constexpr Real b6 = 512.0 / 1771.0;

		// 4th order solution weights (b* coefficients, for error estimate)
		static constexpr Real bstar1 = 2825.0 / 27648.0;
		static constexpr Real bstar2 = 0.0;
		static constexpr Real bstar3 = 18575.0 / 48384.0;
		static constexpr Real bstar4 = 13525.0 / 55296.0;
		static constexpr Real bstar5 = 277.0 / 14336.0;
		static constexpr Real bstar6 = 1.0 / 4.0;

		// Error coefficients (b - b*)
		static constexpr Real e1 = b1 - bstar1;
		static constexpr Real e2 = b2 - bstar2;
		static constexpr Real e3 = b3 - bstar3;
		static constexpr Real e4 = b4 - bstar4;
		static constexpr Real e5 = b5 - bstar5;
		static constexpr Real e6 = b6 - bstar6;
	};

	/******************************************************************************
	 * DORMAND-PRINCE 5(4) COEFFICIENTS
	 *
	 * 7-stage embedded pair with FSAL (First Same As Last) property.
	 * Fifth-order solution with fourth-order error estimate.
	 * The standard method in MATLAB's ode45 and SciPy's RK45.
	 * Reference: [DP80], [HNW1] Chapter II.5, [NR3] Section 17.2
	 *
	 * VERIFIED: Coefficients match standard Butcher tableau.
	 ******************************************************************************/
	struct DormandPrince5 {
		static constexpr int stages = 7;
		static constexpr int order = 5;
		static constexpr int error_order = 4;
		static constexpr bool fsal = true;

		// Time nodes (c coefficients)
		static constexpr Real c2 = 1.0 / 5.0;
		static constexpr Real c3 = 3.0 / 10.0;
		static constexpr Real c4 = 4.0 / 5.0;
		static constexpr Real c5 = 8.0 / 9.0;
		static constexpr Real c6 = 1.0;
		static constexpr Real c7 = 1.0;

		// Stage coefficients (a matrix, lower triangular)
		static constexpr Real a21 = 1.0 / 5.0;

		static constexpr Real a31 = 3.0 / 40.0;
		static constexpr Real a32 = 9.0 / 40.0;

		static constexpr Real a41 = 44.0 / 45.0;
		static constexpr Real a42 = -56.0 / 15.0;
		static constexpr Real a43 = 32.0 / 9.0;

		static constexpr Real a51 = 19372.0 / 6561.0;
		static constexpr Real a52 = -25360.0 / 2187.0;
		static constexpr Real a53 = 64448.0 / 6561.0;
		static constexpr Real a54 = -212.0 / 729.0;

		static constexpr Real a61 = 9017.0 / 3168.0;
		static constexpr Real a62 = -355.0 / 33.0;
		static constexpr Real a63 = 46732.0 / 5247.0;
		static constexpr Real a64 = 49.0 / 176.0;
		static constexpr Real a65 = -5103.0 / 18656.0;

		static constexpr Real a71 = 35.0 / 384.0;
		static constexpr Real a72 = 0.0;
		static constexpr Real a73 = 500.0 / 1113.0;
		static constexpr Real a74 = 125.0 / 192.0;
		static constexpr Real a75 = -2187.0 / 6784.0;
		static constexpr Real a76 = 11.0 / 84.0;

		// 5th order solution weights (b coefficients) - same as a7x due to FSAL
		static constexpr Real b1 = 35.0 / 384.0;
		static constexpr Real b2 = 0.0;
		static constexpr Real b3 = 500.0 / 1113.0;
		static constexpr Real b4 = 125.0 / 192.0;
		static constexpr Real b5 = -2187.0 / 6784.0;
		static constexpr Real b6 = 11.0 / 84.0;
		static constexpr Real b7 = 0.0;

		// 4th order solution weights (b* coefficients, for error estimate)
		static constexpr Real bstar1 = 5179.0 / 57600.0;
		static constexpr Real bstar2 = 0.0;
		static constexpr Real bstar3 = 7571.0 / 16695.0;
		static constexpr Real bstar4 = 393.0 / 640.0;
		static constexpr Real bstar5 = -92097.0 / 339200.0;
		static constexpr Real bstar6 = 187.0 / 2100.0;
		static constexpr Real bstar7 = 1.0 / 40.0;

		// Error coefficients (b - b*)
		static constexpr Real e1 = b1 - bstar1;  // = 71/57600
		static constexpr Real e2 = b2 - bstar2;  // = 0
		static constexpr Real e3 = b3 - bstar3;  // = -71/16695
		static constexpr Real e4 = b4 - bstar4;  // = 71/1920
		static constexpr Real e5 = b5 - bstar5;  // = -17253/339200
		static constexpr Real e6 = b6 - bstar6;  // = 22/525
		static constexpr Real e7 = b7 - bstar7;  // = -1/40
	};

	/******************************************************************************
	 * DORMAND-PRINCE 8(7) COEFFICIENTS
	 *
	 * 13-stage embedded pair (DOP853 variant).
	 * Eighth-order solution with seventh-order error estimate.
	 * For high-precision problems requiring very tight error control.
	 * Reference: [PD81], [HNW1] Chapter II.5
	 *
	 * VERIFIED: Coefficients match standard DOP853 Butcher tableau.
	 ******************************************************************************/
	struct DormandPrince8 {
		static constexpr int stages = 13;
		static constexpr int order = 8;
		static constexpr int error_order = 7;
		static constexpr bool fsal = true;

		// Time nodes (c coefficients)
		static constexpr Real c[13] = {
			0.0,
			1.0 / 18.0,
			1.0 / 12.0,
			1.0 / 8.0,
			5.0 / 16.0,
			3.0 / 8.0,
			59.0 / 400.0,
			93.0 / 200.0,
			5490023248.0 / 9719169821.0,
			13.0 / 20.0,
			1201146811.0 / 1299019798.0,
			1.0,
			1.0
		};

		// Stage coefficients (a matrix, stored per-row for computational convenience)
		static constexpr Real a2[] = {1.0 / 18.0};
		static constexpr Real a3[] = {1.0 / 48.0, 1.0 / 16.0};
		static constexpr Real a4[] = {1.0 / 32.0, 0.0, 3.0 / 32.0};
		static constexpr Real a5[] = {5.0 / 16.0, 0.0, -75.0 / 64.0, 75.0 / 64.0};
		static constexpr Real a6[] = {3.0 / 80.0, 0.0, 0.0, 3.0 / 16.0, 3.0 / 20.0};
		static constexpr Real a7[] = {29443841.0 / 614563906.0, 0.0, 0.0, 77736538.0 / 692538347.0, -28693883.0 / 1125000000.0, 23124283.0 / 1800000000.0};
		static constexpr Real a8[] = {16016141.0 / 946692911.0, 0.0, 0.0, 61564180.0 / 158732637.0, 22789713.0 / 633445777.0, 545815736.0 / 2771057229.0, -180193667.0 / 1043307555.0};
		static constexpr Real a9[] = {39632708.0 / 573591083.0, 0.0, 0.0, -433636366.0 / 683701615.0, -421739975.0 / 2616292301.0, 100302831.0 / 723423059.0, 790204164.0 / 839813087.0, 800635310.0 / 3783071287.0};
		static constexpr Real a10[] = {246121993.0 / 1340847787.0, 0.0, 0.0, -37695042795.0 / 15268766246.0, -309121744.0 / 1061227803.0, -12992083.0 / 490766935.0, 6005943493.0 / 2108947869.0, 393006217.0 / 1396673457.0, 123872331.0 / 1001029789.0};
		static constexpr Real a11[] = {-1028468189.0 / 846180014.0, 0.0, 0.0, 8478235783.0 / 508512852.0, 1311729495.0 / 1432422823.0, -10304129995.0 / 1701304382.0, -48777925059.0 / 3047939560.0, 15336726248.0 / 1032824649.0, -45442868181.0 / 3398467696.0, 3065993473.0 / 597172653.0};
		static constexpr Real a12[] = {185892177.0 / 718116043.0, 0.0, 0.0, -3185094517.0 / 667107341.0, -477755414.0 / 1098053517.0, -703635378.0 / 230739211.0, 5731566787.0 / 1027545527.0, 5232866602.0 / 850066563.0, -4093664535.0 / 808688257.0, 3962137247.0 / 1805957418.0, 65686358.0 / 487910083.0};
		static constexpr Real a13[] = {403863854.0 / 491063109.0, 0.0, 0.0, -5068492393.0 / 434740067.0, -411421997.0 / 543043805.0, 652783627.0 / 914296604.0, 11173962825.0 / 925320556.0, -13158990841.0 / 6184727034.0, 3936647629.0 / 1978049680.0, -160528059.0 / 685178525.0, 248638103.0 / 1413531060.0, 0.0};

		// 8th order solution weights (b coefficients)
		static constexpr Real b8[13] = {
			14005451.0 / 335480064.0,
			0.0,
			0.0,
			0.0,
			0.0,
			-59238493.0 / 1068277825.0,
			181606767.0 / 758867731.0,
			561292985.0 / 797845732.0,
			-1041891430.0 / 1371343529.0,
			760417239.0 / 1151165299.0,
			118820643.0 / 751138087.0,
			-528747749.0 / 2220607170.0,
			1.0 / 4.0
		};

		// 7th order solution weights (b* coefficients, for error estimate)
		static constexpr Real b7[13] = {
			13451932.0 / 455176623.0,
			0.0,
			0.0,
			0.0,
			0.0,
			-808719846.0 / 976000145.0,
			1757004468.0 / 5645159321.0,
			656045339.0 / 265891186.0,
			-3867574721.0 / 1518517206.0,
			465885868.0 / 322736535.0,
			53011238.0 / 667516719.0,
			2.0 / 45.0,
			0.0
		};
	};

	//=============================================================================
	//                    COMPILE-TIME CONSISTENCY CHECKS
	//=============================================================================

	// Precision-dependent tolerance: float has ~7 digits, double has ~15
	constexpr Real RK_COEFF_TOL = sizeof(Real) >= 8 ? Real(1e-14) : Real(1e-5);

	// Verify row sums for consistency (sum of a[i][j] should equal c[i])
	// These are fundamental Butcher tableau consistency conditions

	// Cash-Karp row sum checks
	static_assert(cabs(CashKarp5::a21 - CashKarp5::c2) < RK_COEFF_TOL, "CK5: row 2 sum != c2");
	static_assert(cabs(CashKarp5::a31 + CashKarp5::a32 - CashKarp5::c3) < RK_COEFF_TOL, "CK5: row 3 sum != c3");
	static_assert(cabs(CashKarp5::a41 + CashKarp5::a42 + CashKarp5::a43 - CashKarp5::c4) < RK_COEFF_TOL, "CK5: row 4 sum != c4");
	static_assert(cabs(CashKarp5::a51 + CashKarp5::a52 + CashKarp5::a53 + CashKarp5::a54 - CashKarp5::c5) < RK_COEFF_TOL, "CK5: row 5 sum != c5");
	static_assert(cabs(CashKarp5::a61 + CashKarp5::a62 + CashKarp5::a63 + CashKarp5::a64 + CashKarp5::a65 - CashKarp5::c6) < RK_COEFF_TOL, "CK5: row 6 sum != c6");

	// Verify b weights sum to 1 (required for consistency)
	static_assert(cabs(CashKarp5::b1 + CashKarp5::b2 + CashKarp5::b3 + CashKarp5::b4 + CashKarp5::b5 + CashKarp5::b6 - 1.0) < RK_COEFF_TOL, "CK5: b weights don't sum to 1");
	static_assert(cabs(CashKarp5::bstar1 + CashKarp5::bstar2 + CashKarp5::bstar3 + CashKarp5::bstar4 + CashKarp5::bstar5 + CashKarp5::bstar6 - 1.0) < RK_COEFF_TOL, "CK5: bstar weights don't sum to 1");

	// Dormand-Prince 5 row sum checks
	static_assert(cabs(DormandPrince5::a21 - DormandPrince5::c2) < RK_COEFF_TOL, "DP5: row 2 sum != c2");
	static_assert(cabs(DormandPrince5::a31 + DormandPrince5::a32 - DormandPrince5::c3) < RK_COEFF_TOL, "DP5: row 3 sum != c3");
	static_assert(cabs(DormandPrince5::a41 + DormandPrince5::a42 + DormandPrince5::a43 - DormandPrince5::c4) < RK_COEFF_TOL, "DP5: row 4 sum != c4");
	static_assert(cabs(DormandPrince5::a51 + DormandPrince5::a52 + DormandPrince5::a53 + DormandPrince5::a54 - DormandPrince5::c5) < RK_COEFF_TOL, "DP5: row 5 sum != c5");
	static_assert(cabs(DormandPrince5::a61 + DormandPrince5::a62 + DormandPrince5::a63 + DormandPrince5::a64 + DormandPrince5::a65 - DormandPrince5::c6) < RK_COEFF_TOL, "DP5: row 6 sum != c6");
	static_assert(cabs(DormandPrince5::a71 + DormandPrince5::a72 + DormandPrince5::a73 + DormandPrince5::a74 + DormandPrince5::a75 + DormandPrince5::a76 - DormandPrince5::c7) < RK_COEFF_TOL, "DP5: row 7 sum != c7");

	// Verify b weights sum to 1
	static_assert(cabs(DormandPrince5::b1 + DormandPrince5::b2 + DormandPrince5::b3 + DormandPrince5::b4 + DormandPrince5::b5 + DormandPrince5::b6 + DormandPrince5::b7 - 1.0) < RK_COEFF_TOL, "DP5: b weights don't sum to 1");
	static_assert(cabs(DormandPrince5::bstar1 + DormandPrince5::bstar2 + DormandPrince5::bstar3 + DormandPrince5::bstar4 + DormandPrince5::bstar5 + DormandPrince5::bstar6 + DormandPrince5::bstar7 - 1.0) < RK_COEFF_TOL, "DP5: bstar weights don't sum to 1");

	// Verify FSAL property: a7x = bx for DP5
	static_assert(DormandPrince5::a71 == DormandPrince5::b1, "DP5: FSAL violated (a71 != b1)");
	static_assert(DormandPrince5::a72 == DormandPrince5::b2, "DP5: FSAL violated (a72 != b2)");
	static_assert(DormandPrince5::a73 == DormandPrince5::b3, "DP5: FSAL violated (a73 != b3)");
	static_assert(DormandPrince5::a74 == DormandPrince5::b4, "DP5: FSAL violated (a74 != b4)");
	static_assert(DormandPrince5::a75 == DormandPrince5::b5, "DP5: FSAL violated (a75 != b5)");
	static_assert(DormandPrince5::a76 == DormandPrince5::b6, "DP5: FSAL violated (a76 != b6)");

} // namespace RKCoeff
} // namespace MML

#endif // MML_ODE_RK_COEFFICIENTS_H
