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

#include <mml/MMLBase.h>
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

		// Uniform tableau view consumed by ExplicitRKStageEvaluator. Named coefficients
		// above remain available for source compatibility and coefficient verification.
		static constexpr Real nodes[stages] = {0.0, c2, c3, c4, c5, c6};
		static constexpr Real stage_coefficients[stages][stages] = {
			{0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
			{a21, 0.0, 0.0, 0.0, 0.0, 0.0},
			{a31, a32, 0.0, 0.0, 0.0, 0.0},
			{a41, a42, a43, 0.0, 0.0, 0.0},
			{a51, a52, a53, a54, 0.0, 0.0},
			{a61, a62, a63, a64, a65, 0.0}
		};
		static constexpr Real solution_weights[stages] = {b1, b2, b3, b4, b5, b6};
		static constexpr Real error_weights[stages] = {e1, e2, e3, e4, e5, e6};

		static constexpr Real node(int stage) { return nodes[stage]; }
		static constexpr Real stageCoefficient(int stage, int previousStage) {
			return stage_coefficients[stage][previousStage];
		}
		static constexpr Real solutionWeight(int stage) { return solution_weights[stage]; }
		static constexpr Real errorWeight(int stage) { return error_weights[stage]; }
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

		// Uniform tableau view consumed by ExplicitRKStageEvaluator. Named coefficients
		// above remain available for source compatibility and coefficient verification.
		static constexpr Real nodes[stages] = {0.0, c2, c3, c4, c5, c6, c7};
		static constexpr Real stage_coefficients[stages][stages] = {
			{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
			{a21, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
			{a31, a32, 0.0, 0.0, 0.0, 0.0, 0.0},
			{a41, a42, a43, 0.0, 0.0, 0.0, 0.0},
			{a51, a52, a53, a54, 0.0, 0.0, 0.0},
			{a61, a62, a63, a64, a65, 0.0, 0.0},
			{a71, a72, a73, a74, a75, a76, 0.0}
		};
		static constexpr Real solution_weights[stages] = {b1, b2, b3, b4, b5, b6, b7};
		static constexpr Real error_weights[stages] = {e1, e2, e3, e4, e5, e6, e7};

		static constexpr Real node(int stage) { return nodes[stage]; }
		static constexpr Real stageCoefficient(int stage, int previousStage) {
			return stage_coefficients[stage][previousStage];
		}
		static constexpr Real solutionWeight(int stage) { return solution_weights[stage]; }
		static constexpr Real errorWeight(int stage) { return error_weights[stage]; }
	};

	/******************************************************************************
	 * DOP853 8(5,3) COEFFICIENTS
	 *
	 * Hairer's 12-stage eighth-order propagation formula with embedded fifth-
	 * and third-order estimators. An endpoint derivative and three additional
	 * stages support the method's seventh-order dense output.
	 * Reference: Hairer, Norsett & Wanner, Solving ODEs I, Section II.5.
	 ******************************************************************************/
	struct DormandPrince8 {
		static constexpr int stages = 12;
		static constexpr int extended_stages = 16;
		static constexpr int order = 8;
		static constexpr int error_order = 7;
		static constexpr bool fsal = false;

		static constexpr Real c[extended_stages] = {
			0.0, 0.526001519587677318785587544488e-1, 0.789002279381515978178381316732e-1,
			0.118350341907227396726757197510, 0.281649658092772603273242802490,
			1.0 / 3.0, 0.25, 4.0 / 13.0, 0.651282051282051282051282051282,
			0.6, 6.0 / 7.0, 1.0, 1.0, 0.1, 0.2, 7.0 / 9.0
		};

		static constexpr Real a[extended_stages][extended_stages] = {
			{},
			{5.26001519587677318785587544488e-2},
			{1.97250569845378994544595329183e-2, 5.91751709536136983633785987549e-2},
			{2.95875854768068491816892993775e-2, 0, 8.87627564304205475450678981324e-2},
			{2.41365134159266685502369798665e-1, 0, -8.84549479328286085344864962717e-1, 9.24834003261792003115737966543e-1},
			{3.7037037037037037037037037037e-2, 0, 0, 1.70828608729473871279604482173e-1, 1.25467687566822425016691814123e-1},
			{3.7109375e-2, 0, 0, 1.70252211019544039314978060272e-1, 6.02165389804559606850219397283e-2, -1.7578125e-2},
			{3.70920001185047927108779319836e-2, 0, 0, 1.70383925712239993810214054705e-1, 1.07262030446373284651809199168e-1, -1.53194377486244017527936158236e-2, 8.27378916381402288758473766002e-3},
			{6.24110958716075717114429577812e-1, 0, 0, -3.36089262944694129406857109825, -8.68219346841726006818189891453e-1, 2.75920996994467083049415600797e1, 2.01540675504778934086186788979e1, -4.34898841810699588477366255144e1},
			{4.77662536438264365890433908527e-1, 0, 0, -2.48811461997166764192642586468, -5.90290826836842996371446475743e-1, 2.12300514481811942347288949897e1, 1.52792336328824235832596922938e1, -3.32882109689848629194453265587e1, -2.03312017085086261358222928593e-2},
			{-9.3714243008598732571704021658e-1, 0, 0, 5.18637242884406370830023853209, 1.09143734899672957818500254654, -8.14978701074692612513997267357, -1.85200656599969598641566180701e1, 2.27394870993505042818970056734e1, 2.49360555267965238987089396762, -3.0467644718982195003823669022},
			{2.27331014751653820792359768449, 0, 0, -1.05344954667372501984066689879e1, -2.00087205822486249909675718444, -1.79589318631187989172765950534e1, 2.79488845294199600508499808837e1, -2.85899827713502369474065508674, -8.87285693353062954433549289258, 1.23605671757943030647266201528e1, 6.43392746015763530355970484046e-1},
			{5.42937341165687622380535766363e-2, 0, 0, 0, 0, 4.45031289275240888144113950566, 1.89151789931450038304281599044, -5.8012039600105847814672114227, 3.1116436695781989440891606237e-1, -1.52160949662516078556178806805e-1, 2.01365400804030348374776537501e-1, 4.47106157277725905176885569043e-2},
			{5.61675022830479523392909219681e-2, 0, 0, 0, 0, 0, 2.53500210216624811088794765333e-1, -2.46239037470802489917441475441e-1, -1.24191423263816360469010140626e-1, 1.5329179827876569731206322685e-1, 8.20105229563468988491666602057e-3, 7.56789766054569976138603589584e-3, -8.298e-3},
			{3.18346481635021405060768473261e-2, 0, 0, 0, 0, 2.83009096723667755288322961402e-2, 5.35419883074385676223797384372e-2, -5.49237485713909884646569340306e-2, 0, 0, -1.08347328697249322858509316994e-4, 3.82571090835658412954920192323e-4, -3.40465008687404560802977114492e-4, 1.41312443674632500278074618366e-1},
			{-4.28896301583791923408573538692e-1, 0, 0, 0, 0, -4.69762141536116384314449447206, 7.68342119606259904184240953878, 4.06898981839711007970213554331, 3.56727187455281109270669543021e-1, 0, 0, 0, -1.39902416515901462129418009734e-3, 2.9475147891527723389556272149, -9.15095847217987001081870187138}
		};

		static constexpr Real e3[13] = {
			5.42937341165687622380535766363e-2 - 0.244094488188976377952755905512,
			0, 0, 0, 0, 4.45031289275240888144113950566, 1.89151789931450038304281599044,
			-5.8012039600105847814672114227,
			3.1116436695781989440891606237e-1 - 0.733846688281611857341361741547,
			-1.52160949662516078556178806805e-1, 2.01365400804030348374776537501e-1,
			4.47106157277725905176885569043e-2 - 0.220588235294117647058823529412e-1, 0
		};
		static constexpr Real e5[13] = {
			0.1312004499419488073250102996e-1, 0, 0, 0, 0,
			-0.1225156446376204440720569753e1, -0.4957589496572501915214079952,
			0.1664377182454986536961530415e1, -0.3503288487499736816886487290,
			0.3341791187130174790297318841, 0.8192320648511571246570742613e-1,
			-0.2235530786388629525884427845e-1, 0
		};
		static constexpr Real d[4][extended_stages] = {
			{-0.84289382761090128651353491142e1, 0, 0, 0, 0, 0.56671495351937776962531783590, -0.30689499459498916912797304727e1, 0.23846676565120698287728149680e1, 0.21170345824450282767155149946e1, -0.87139158377797299206789907490, 0.22404374302607882758541771650e1, 0.63157877876946881815570249290, -0.88990336451333310820698117400e-1, 0.18148505520854727256656404962e2, -0.91946323924783554000451984436e1, -0.44360363875948939664310572000e1},
			{0.10427508642579134603413151009e2, 0, 0, 0, 0, 0.24228349177525818288430175319e3, 0.16520045171727028198505394887e3, -0.37454675472269020279518312152e3, -0.22113666853125306036270938578e2, 0.77334326684722638389603898808e1, -0.30674084731089398182061213626e2, -0.93321305264302278729567221706e1, 0.15697238121770843886131091075e2, -0.31139403219565177677282850411e2, -0.93529243588444783865713862664e1, 0.35816841486394083752465898540e2},
			{0.19985053242002433820987653617e2, 0, 0, 0, 0, -0.38703730874935176555105901742e3, -0.18917813819516756882830838328e3, 0.52780815920542364900561016686e3, -0.11573902539959630126141871134e2, 0.68812326946963000169666922661e1, -0.10006050966910838403183860980e1, 0.77771377980534432092869265740, -0.27782057523535084065932004339e1, -0.60196695231264120758267380846e2, 0.84320405506677161018159903784e2, 0.11992291136182789328035130030e2},
			{-0.25693933462703749003312586129e2, 0, 0, 0, 0, -0.15418974869023643374053993627e3, -0.23152937917604549567536039109e3, 0.35763911791061412378285349910e3, 0.93405324183624310003907691704e2, -0.37458323136451633156875139351e2, 0.10409964950896230045147246184e3, 0.29840293426660503123344363579e2, -0.43533456590011143754432175058e2, 0.96324553959188282948394950600e2, -0.39177261675615439165231486172e2, -0.14972683625798562581422125276e3}
		};

		static constexpr Real node(int stage) { return c[stage]; }
		static constexpr Real stageCoefficient(int stage, int previousStage) { return a[stage][previousStage]; }
		static constexpr Real solutionWeight(int stage) { return a[12][stage]; }
		static constexpr Real errorWeight(int stage) { return e5[stage]; }
		static constexpr Real error3Weight(int stage) { return e3[stage]; }
		static constexpr Real error5Weight(int stage) { return e5[stage]; }
		static constexpr Real denseWeight(int coefficient, int stage) { return d[coefficient][stage]; }
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
