///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        GaussKronrod.h                                                      ///
///  Description: Gauss-Kronrod adaptive quadrature (G7K15, G10K21, G15K31)           ///
///               Industry standard for general-purpose numerical integration         ///
///               Used in QUADPACK, GSL, scipy, MATLAB, etc.                          ///
///                                                                                   ///
///  REFERENCES:                                                                      ///
///    [K65]     Kronrod, A.S. (1965). Nodes and Weights of Quadrature Formulas       ///
///    [PS73]    Patterson, T.N.L. (1968). Math. Comp. 22(104), pp. 847-856          ///
///    [PSWDK83] Piessens et al. (1983). QUADPACK: Automatic Integration. Springer   ///
///    [NR3]     Press et al., Numerical Recipes 3rd ed., Ch. 4.5                    ///
///                                                                                   ///
///  See references/book_references.md and references/paperes_references.md          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_GAUSS_KRONROD_H
#define MML_GAUSS_KRONROD_H

#include <vector>
#include <cmath>
#include <limits>
#include <functional>

#include "MMLBase.h"
#include "core/Integration/IntegrationBase.h"

namespace MML {
	namespace Integration 
	{
		/// @brief Result of Gauss-Kronrod integration
		struct GKResult {
			Real value;					 ///< Integral estimate (from Kronrod rule)
			Real error_estimate; ///< Error estimate |Gauss - Kronrod|
			int function_evals;	 ///< Number of function evaluations
			bool converged;			 ///< True if error < tolerance

			GKResult(Real v = 0.0, Real e = 0.0, int n = 0, bool c = false)
					: value(v)
					, error_estimate(e)
					, function_evals(n)
					, converged(c) {}
		};

		/*************************************************************************/
		/*****              G7K15 GAUSS-KRONROD RULE (15 points)             *****/
		/*****  7-point Gauss embedded in 15-point Kronrod extension         *****/
		/*****  Most commonly used rule for adaptive integration             *****/
		/*************************************************************************/

		/// @brief 15-point Kronrod nodes (symmetric about 0)
		/// @details Only positive nodes listed; negative nodes are -xGK15[i]
		///          Node 0 is at origin (weight wGK15_0)
		static constexpr int GK15_N = 8; // Number of unique positive nodes + origin
		static constexpr Real xGK15[GK15_N] = {
			0.0,																	 // x0 (origin)
			0.20778495500789846760068940377324491, // x1
			0.40584515137739716690660641207696146, // x2
			0.58608723546769113029414483825872959, // x3
			0.74153118559939443986386477328078841, // x4
			0.86486442335976907278971278864092620, // x5
			0.94910791234275852452618968404785126, // x6
			0.99145537112081263920685469752632851	 // x7
		};

		/// @brief 15-point Kronrod weights (for all 15 points including symmetric)
		static constexpr Real wGK15[GK15_N] = {
			0.20948214108472782801299917489171426, // w0 (origin)
			0.20443294007529889241416199923464908, // w1
			0.19035057806478540991325640242101368, // w2
			0.16900472663926790282658342659855028, // w3
			0.14065325971552591874518959051023792, // w4
			0.10479001032225018383987632254151801, // w5
			0.06309209262997855329070066318920429, // w6
			0.02293532201052922496373200805896959	 // w7
		};

		/// @brief 7-point Gauss weights (for odd-indexed Kronrod nodes only)
		/// @details Gauss nodes are xGK15[1,3,5,7] (every other Kronrod node starting from x1)
		static constexpr int G7_N = 4; // Half of 7 (symmetric) + different handling
		static constexpr Real wG7[G7_N] = {
			0.41795918367346938775510204081632653, // w0 (at x=0, but stored at index 0)
			0.38183005050511894495036977548897513, // w1 (at xGK15[2])
			0.27970539148927666790146777142377958, // w2 (at xGK15[4])
			0.12948496616886969327061143267908202	 // w3 (at xGK15[6])
		};

		/// @brief Apply G7K15 rule to integrate f over [a,b].
		/// @details Evaluates function at 15 points to compute both 7-point Gauss
		///          and 15-point Kronrod estimates. Error = |Gauss - Kronrod|.
		/// @tparam Func Callable with signature Real(Real)
		/// @param f Function to integrate
		/// @param a Lower limit
		/// @param b Upper limit
		/// @return GKResult with integral estimate and error
		template<typename Func>
		static GKResult IntegrateGK15(Func f, Real a, Real b) {
			const Real center = 0.5 * (a + b);
			const Real half_length = 0.5 * (b - a);

			// Function value at center (used by both Gauss and Kronrod)
			Real f_center = f(center);

			// Compute sums
			Real result_gauss = wG7[0] * f_center;		 // G7 weight at origin
			Real result_kronrod = wGK15[0] * f_center; // K15 weight at origin

			// Symmetric pairs
			for (int i = 1; i < GK15_N; i++) {
				Real x_scaled = half_length * xGK15[i];
				Real f_plus = f(center + x_scaled);
				Real f_minus = f(center - x_scaled);
				Real f_sum = f_plus + f_minus;

				// All points contribute to Kronrod
				result_kronrod += wGK15[i] * f_sum;

				// Only odd-indexed points (i = 2, 4, 6) contribute to Gauss
				// i=2 uses wG7[1], i=4 uses wG7[2], i=6 uses wG7[3]
				if (i == 2)
					result_gauss += wG7[1] * f_sum;
				else if (i == 4)
					result_gauss += wG7[2] * f_sum;
				else if (i == 6)
					result_gauss += wG7[3] * f_sum;
			}

			// Scale by half-length for change of variables
			result_gauss *= half_length;
			result_kronrod *= half_length;

			Real error = std::abs(result_kronrod - result_gauss);

			return GKResult(result_kronrod, error, 15, true);
		}

		/*************************************************************************/
		/*****              G10K21 GAUSS-KRONROD RULE (21 points)            *****/
		/*****  10-point Gauss embedded in 21-point Kronrod extension        *****/
		/*************************************************************************/

		/// @brief Apply G10K21 rule to integrate f over [a,b].
		/// @details Uses tabulated nodes and weights from GSL (QUADPACK).
		///          21 Kronrod points = 11 unique abscissae (mirrored about center)
		///          10 Gauss points at odd Kronrod indices (1,3,5,7,9)
		template<typename Func>
		static GKResult IntegrateGK21(Func f, Real a, Real b) {
			// 21-point Kronrod abscissae (11 unique, from GSL qk21.c)
			// xgk[1], xgk[3], ... are the 10-point Gauss abscissae
			// xgk[0], xgk[2], ... are added by Kronrod extension
			static constexpr Real xgk[11] = {
				0.995657163025808080735527280689003, // 0 - Kronrod only
				0.973906528517171720077964012084452, // 1 - Gauss
				0.930157491355708226001207180059508, // 2 - Kronrod only
				0.865063366688984510732096688423493, // 3 - Gauss
				0.780817726586416897063717578345042, // 4 - Kronrod only
				0.679409568299024406234327365114874, // 5 - Gauss
				0.562757134668604683339000099272694, // 6 - Kronrod only
				0.433395394129247190799265943165784, // 7 - Gauss
				0.294392862701460198131126603103866, // 8 - Kronrod only
				0.148874338981631210884826001129720, // 9 - Gauss
				0.000000000000000000000000000000000	 // 10 - Center (Kronrod only)
			};

			// 10-point Gauss weights (for odd-indexed Kronrod nodes)
			static constexpr Real wg[5] = {
				0.066671344308688137593568809893332, // for xgk[1]
				0.149451349150580593145776339657697, // for xgk[3]
				0.219086362515982043995534934228163, // for xgk[5]
				0.269266719309996355091226921569469, // for xgk[7]
				0.295524224714752870173892994651338	 // for xgk[9]
			};

			// 21-point Kronrod weights (for all nodes)
			static constexpr Real wgk[11] = {
				0.011694638867371874278064396062192, // for xgk[0]
				0.032558162307964727478818972459390, // for xgk[1]
				0.054755896574351996031381300244580, // for xgk[2]
				0.075039674810919952767043140916190, // for xgk[3]
				0.093125454583697605535065465083366, // for xgk[4]
				0.109387158802297641899210590325805, // for xgk[5]
				0.123491976262065851077958109831074, // for xgk[6]
				0.134709217311473325928054001771707, // for xgk[7]
				0.142775938577060080797094273138717, // for xgk[8]
				0.147739104901338491374841515972068, // for xgk[9]
				0.149445554002916905664936468389821	 // for xgk[10] (center)
			};

			const Real center = 0.5 * (a + b);
			const Real half_length = 0.5 * (b - a);

			Real f_center = f(center);

			// Center point (index 10) - Kronrod only (n=11 is odd, so no Gauss at center)
			Real result_kronrod = wgk[10] * f_center;
			Real result_gauss = 0.0;

			// Process symmetric pairs (indices 0-9)
			for (int j = 0; j < 5; j++) {
				// Odd indices (1,3,5,7,9) contribute to both Gauss and Kronrod
				int jtw = j * 2 + 1; // 1, 3, 5, 7, 9
				Real x_scaled = half_length * xgk[jtw];
				Real f_plus = f(center + x_scaled);
				Real f_minus = f(center - x_scaled);
				Real f_sum = f_plus + f_minus;

				result_gauss += wg[j] * f_sum;
				result_kronrod += wgk[jtw] * f_sum;
			}

			for (int j = 0; j < 5; j++) {
				// Even indices (0,2,4,6,8) contribute to Kronrod only
				int jtwm1 = j * 2; // 0, 2, 4, 6, 8
				Real x_scaled = half_length * xgk[jtwm1];
				Real f_plus = f(center + x_scaled);
				Real f_minus = f(center - x_scaled);
				Real f_sum = f_plus + f_minus;

				result_kronrod += wgk[jtwm1] * f_sum;
			}

			result_gauss *= half_length;
			result_kronrod *= half_length;

			return GKResult(result_kronrod, std::abs(result_kronrod - result_gauss), 21, true);
		}

		/*************************************************************************/
		/*****              G15K31 GAUSS-KRONROD RULE (31 points)            *****/
		/*****  15-point Gauss embedded in 31-point Kronrod extension        *****/
		/*************************************************************************/

		/// @brief Apply G15K31 rule to integrate f over [a,b].
		/// @details Uses tabulated nodes and weights from GSL (QUADPACK).
		///          31 Kronrod points = 16 unique abscissae (mirrored about center)
		///          15 Gauss points at odd Kronrod indices (1,3,5,7,9,11,13,15)
		///          Note: n=16 is even, so center IS a Gauss point
		template<typename Func>
		static GKResult IntegrateGK31(Func f, Real a, Real b) {
			// 31-point Kronrod abscissae (16 unique, from GSL qk31.c)
			// xgk[1], xgk[3], ... are the 15-point Gauss abscissae
			// xgk[0], xgk[2], ... are added by Kronrod extension
			static constexpr Real xgk[16] = {
				0.998002298693397060285172840152271, // 0 - Kronrod only
				0.987992518020485428489565718586613, // 1 - Gauss
				0.967739075679139134257347978784337, // 2 - Kronrod only
				0.937273392400705904307758947710209, // 3 - Gauss
				0.897264532344081900882509656454496, // 4 - Kronrod only
				0.848206583410427216200648320774217, // 5 - Gauss
				0.790418501442465932967649294817947, // 6 - Kronrod only
				0.724417731360170047416186054613938, // 7 - Gauss
				0.650996741297416970533735895313275, // 8 - Kronrod only
				0.570972172608538847537226737253911, // 9 - Gauss
				0.485081863640239680693655740232351, // 10 - Kronrod only
				0.394151347077563369897207370981045, // 11 - Gauss
				0.299180007153168812166780024266389, // 12 - Kronrod only
				0.201194093997434522300628303394596, // 13 - Gauss
				0.101142066918717499027074231447392, // 14 - Kronrod only
				0.000000000000000000000000000000000	 // 15 - Center (Gauss + Kronrod, n=16 even)
			};

			// 15-point Gauss weights (for odd-indexed nodes plus center)
			// wg[0..6] for xgk[1,3,5,7,9,11,13], wg[7] for center xgk[15]
			static constexpr Real wg[8] = {
				0.030753241996117268354628393577204, // for xgk[1]
				0.070366047488108124709267416450667, // for xgk[3]
				0.107159220467171935011869546685869, // for xgk[5]
				0.139570677926154314447804794511028, // for xgk[7]
				0.166269205816993933553200860481209, // for xgk[9]
				0.186161000015562211026800561866423, // for xgk[11]
				0.198431485327111576456118326443839, // for xgk[13]
				0.202578241925561272880620199967519	 // for xgk[15] (center)
			};

			// 31-point Kronrod weights (for all nodes)
			static constexpr Real wgk[16] = {
				0.005377479872923348987792051430128, // for xgk[0]
				0.015007947329316122538374763075807, // for xgk[1]
				0.025460847326715320186874001019653, // for xgk[2]
				0.035346360791375846222037948478360, // for xgk[3]
				0.044589751324764876608227299373280, // for xgk[4]
				0.053481524690928087265343147239430, // for xgk[5]
				0.062009567800670640285139230960803, // for xgk[6]
				0.069854121318728258709520077099147, // for xgk[7]
				0.076849680757720378894432777482659, // for xgk[8]
				0.083080502823133021038289247286104, // for xgk[9]
				0.088564443056211770647275443693774, // for xgk[10]
				0.093126598170825321225486872747346, // for xgk[11]
				0.096642726983623678505179907627589, // for xgk[12]
				0.099173598721791959332393173484603, // for xgk[13]
				0.100769845523875595044946662617570, // for xgk[14]
				0.101330007014791549017374792767493	 // for xgk[15] (center)
			};

			const Real center = 0.5 * (a + b);
			const Real half_length = 0.5 * (b - a);

			Real f_center = f(center);

			// Center point (index 15) - Both Gauss and Kronrod (n=16 is even)
			Real result_kronrod = wgk[15] * f_center;
			Real result_gauss = wg[7] * f_center; // wg[n/2-1] = wg[7]

			// Process symmetric pairs (indices 0-14)
			for (int j = 0; j < 7; j++) {
				// Odd indices (1,3,5,7,9,11,13) contribute to both Gauss and Kronrod
				int jtw = j * 2 + 1; // 1, 3, 5, 7, 9, 11, 13
				Real x_scaled = half_length * xgk[jtw];
				Real f_plus = f(center + x_scaled);
				Real f_minus = f(center - x_scaled);
				Real f_sum = f_plus + f_minus;

				result_gauss += wg[j] * f_sum;
				result_kronrod += wgk[jtw] * f_sum;
			}

			for (int j = 0; j < 8; j++) {
				// Even indices (0,2,4,6,8,10,12,14) contribute to Kronrod only
				int jtwm1 = j * 2; // 0, 2, 4, 6, 8, 10, 12, 14
				Real x_scaled = half_length * xgk[jtwm1];
				Real f_plus = f(center + x_scaled);
				Real f_minus = f(center - x_scaled);
				Real f_sum = f_plus + f_minus;

				result_kronrod += wgk[jtwm1] * f_sum;
			}

			result_gauss *= half_length;
			result_kronrod *= half_length;

			return GKResult(result_kronrod, std::abs(result_kronrod - result_gauss), 31, true);
		}

		/*************************************************************************/
		/*****              ADAPTIVE GAUSS-KRONROD INTEGRATION               *****/
		/*****  Recursively subdivides until error tolerance met             *****/
		/*************************************************************************/

		/// @brief Rule selection for adaptive Gauss-Kronrod
		enum class GKRule { GK15, GK21, GK31 };

		/// @brief Adaptive Gauss-Kronrod integration with automatic subdivision.
		/// @details Uses recursive bisection when local error exceeds tolerance.
		///          This is the standard approach used by QUADPACK/GSL.
		/// @tparam Func Callable with signature Real(Real)
		/// @param f Function to integrate
		/// @param a Lower limit
		/// @param b Upper limit
		/// @param tol_abs Absolute error tolerance (default 1e-10)
		/// @param tol_rel Relative error tolerance (default 1e-10)
		/// @param max_depth Maximum recursion depth (default 50)
		/// @param rule Which Gauss-Kronrod rule to use (default GK15)
		/// @return GKResult with integral and accumulated error estimate
		template<typename Func>
		static GKResult IntegrateGKAdaptive(Func f, Real a, Real b, Real tol_abs = 1e-10, Real tol_rel = 1e-10, int max_depth = 50,
																				GKRule rule = GKRule::GK15) {
			// Apply basic rule first
			GKResult local = (rule == GKRule::GK15)		? IntegrateGK15(f, a, b)
											 : (rule == GKRule::GK21) ? IntegrateGK21(f, a, b)
																								: IntegrateGK31(f, a, b);

			Real tolerance = std::max(tol_abs, tol_rel * std::abs(local.value));

			// Check if we meet tolerance or hit max depth
			if (local.error_estimate <= tolerance || max_depth <= 0) {
				local.converged = (local.error_estimate <= tolerance);
				return local;
			}

			// Subdivide and recurse
			Real mid = 0.5 * (a + b);
			GKResult left = IntegrateGKAdaptive(f, a, mid, tol_abs / 2.0, tol_rel, max_depth - 1, rule);
			GKResult right = IntegrateGKAdaptive(f, mid, b, tol_abs / 2.0, tol_rel, max_depth - 1, rule);

			return GKResult(left.value + right.value, left.error_estimate + right.error_estimate, left.function_evals + right.function_evals,
											left.converged && right.converged);
		}

		/*************************************************************************/
		/*****                    CONVENIENCE FUNCTIONS                      *****/
		/*************************************************************************/

		/// @brief Integrate using G7K15 rule (non-adaptive, single application).
		/// @param f Function object (IRealFunction or callable)
		/// @param a Lower limit
		/// @param b Upper limit
		/// @return Integral value
		template<typename Func>
		static Real IntegrateGaussKronrod15(Func f, Real a, Real b) {
			return IntegrateGK15(f, a, b).value;
		}

		/// @brief Integrate using G10K21 rule (non-adaptive, single application).
		template<typename Func>
		static Real IntegrateGaussKronrod21(Func f, Real a, Real b) {
			return IntegrateGK21(f, a, b).value;
		}

		/// @brief Integrate using G15K31 rule (non-adaptive, single application).
		template<typename Func>
		static Real IntegrateGaussKronrod31(Func f, Real a, Real b) {
			return IntegrateGK31(f, a, b).value;
		}

		/// @brief Adaptive Gauss-Kronrod integration (recommended general use).
		/// @param f Function to integrate
		/// @param a Lower limit
		/// @param b Upper limit
		/// @param tol Absolute tolerance (default 1e-10)
		/// @return Integral value
		template<typename Func>
		static Real IntegrateGaussKronrod(Func f, Real a, Real b, Real tol = 1e-10) {
			return IntegrateGKAdaptive(f, a, b, tol, tol, 50, GKRule::GK15).value;
		}

		/// @brief Adaptive Gauss-Kronrod with error estimate.
		/// @param f Function to integrate
		/// @param a Lower limit
		/// @param b Upper limit
		/// @param error Output parameter for error estimate
		/// @param tol Tolerance
		/// @return Integral value
		template<typename Func>
		static Real IntegrateGaussKronrod(Func f, Real a, Real b, Real* error, Real tol = 1e-10) {
			auto result = IntegrateGKAdaptive(f, a, b, tol, tol, 50, GKRule::GK15);
			if (error)
				*error = result.error_estimate;
			return result.value;
		}

		/*****************************************************************/
		/*****        Detailed API - Gauss-Kronrod                  *****/
		/*****************************************************************/

		/// Helper to populate IntegrationDetailedResult from GKResult
		namespace GKDetail
		{
			inline void PopulateFromGK(IntegrationDetailedResult& out,
			                           const GKResult& gk,
			                           const char* algorithm_name)
			{
				out.value = gk.value;
				out.error_estimate = gk.error_estimate;
				out.function_evaluations = gk.function_evals;
				out.converged = gk.converged;
				if (!gk.converged) {
					out.status = AlgorithmStatus::MaxIterationsExceeded;
					out.error_message = std::string(algorithm_name) +
						" did not converge (error=" + std::to_string(gk.error_estimate) + ")";
				}
			}
		} // namespace GKDetail

		/// GK15 integration with full diagnostics
		template<typename Func>
		static IntegrationDetailedResult IntegrateGK15Detailed(
			Func f, Real a, Real b,
			const IntegrationConfig& config = {})
		{
			return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
				"IntegrateGK15", config,
				[&](IntegrationDetailedResult& result) {
					auto gk = IntegrateGK15(f, a, b);
					GKDetail::PopulateFromGK(result, gk, "IntegrateGK15");
				});
		}

		/// GK21 integration with full diagnostics
		template<typename Func>
		static IntegrationDetailedResult IntegrateGK21Detailed(
			Func f, Real a, Real b,
			const IntegrationConfig& config = {})
		{
			return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
				"IntegrateGK21", config,
				[&](IntegrationDetailedResult& result) {
					auto gk = IntegrateGK21(f, a, b);
					GKDetail::PopulateFromGK(result, gk, "IntegrateGK21");
				});
		}

		/// GK31 integration with full diagnostics
		template<typename Func>
		static IntegrationDetailedResult IntegrateGK31Detailed(
			Func f, Real a, Real b,
			const IntegrationConfig& config = {})
		{
			return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
				"IntegrateGK31", config,
				[&](IntegrationDetailedResult& result) {
					auto gk = IntegrateGK31(f, a, b);
					GKDetail::PopulateFromGK(result, gk, "IntegrateGK31");
				});
		}

		/// Adaptive Gauss-Kronrod with full diagnostics
		template<typename Func>
		static IntegrationDetailedResult IntegrateGKAdaptiveDetailed(
			Func f, Real a, Real b,
			const IntegrationConfig& config = {},
			Real tol_abs = 1e-10, Real tol_rel = 1e-10,
			int max_depth = 50, GKRule rule = GKRule::GK15)
		{
			return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
				"IntegrateGKAdaptive", config,
				[&](IntegrationDetailedResult& result) {
					auto gk = IntegrateGKAdaptive(f, a, b, tol_abs, tol_rel, max_depth, rule);
					GKDetail::PopulateFromGK(result, gk, "IntegrateGKAdaptive");
				});
		}

	} // namespace Integration
} // namespace MML

#endif // MML_GAUSS_KRONROD_H
