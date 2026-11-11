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
#include <algorithm>
#include <cmath>
#include <limits>
#include <functional>
#include <array>

#include <mml/MMLBase.h>
#include <mml/core/Integration/IntegrationBase.h>

namespace MML {
	namespace Integration 
	{
		/// @brief Result of Gauss-Kronrod integration: the common IntegrationResult with the GK
		///        constructor order (value, error, function evaluations, converged).
		struct GKResult : public IntegrationResult {
			Real roundoff_floor = 0; ///< Arithmetic floor used by adaptive control, separate from the fixed embedded estimate.
			GKResult(Real v = 0.0, Real e = 0.0, int n = 0, bool c = false)
				: IntegrationResult(v, e, 1, c, /*error_estimated*/ true)
			{
				function_evaluations = n;
				algorithm_name = "GaussKronrod";
			}
			GKResult(const IntegrationResult& r) : IntegrationResult(r) {}
			static GKResult Fixed(Real value, Real error, int evaluations) {
				IntegrationDetail::RequireFinite(value, "Gauss-Kronrod");
				IntegrationDetail::RequireFinite(error, "Gauss-Kronrod error");
				GKResult result(value, error, evaluations, false);
				result.status = AlgorithmStatus::Completed;
				return result;
			}

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
			MML::Real(0.20778495500789846760068940377324491L), // x1
			MML::Real(0.40584515137739716690660641207696146L), // x2
			MML::Real(0.58608723546769113029414483825872959L), // x3
			MML::Real(0.74153118559939443986386477328078841L), // x4
			MML::Real(0.86486442335976907278971278864092620L), // x5
			MML::Real(0.94910791234275852452618968404785126L), // x6
			MML::Real(0.99145537112081263920685469752632851L)	 // x7
		};

		/// @brief 15-point Kronrod weights (for all 15 points including symmetric)
		static constexpr Real wGK15[GK15_N] = {
			MML::Real(0.20948214108472782801299917489171426L), // w0 (origin)
			MML::Real(0.20443294007529889241416199923464908L), // w1
			MML::Real(0.19035057806478540991325640242101368L), // w2
			MML::Real(0.16900472663926790282658342659855028L), // w3
			MML::Real(0.14065325971552591874518959051023792L), // w4
			MML::Real(0.10479001032225018383987632254151801L), // w5
			MML::Real(0.06309209262997855329070066318920429L), // w6
			MML::Real(0.02293532201052922496373200805896959L)	 // w7
		};

		/// @brief 7-point Gauss weights (for odd-indexed Kronrod nodes only)
		/// @details Gauss nodes are xGK15[1,3,5,7] (every other Kronrod node starting from x1)
		static constexpr int G7_N = 4; // Half of 7 (symmetric) + different handling
		static constexpr Real wG7[G7_N] = {
			MML::Real(0.41795918367346938775510204081632653L), // w0 (at x=0, but stored at index 0)
			MML::Real(0.38183005050511894495036977548897513L), // w1 (at xGK15[2])
			MML::Real(0.27970539148927666790146777142377958L), // w2 (at xGK15[4])
			MML::Real(0.12948496616886969327061143267908202L)	 // w3 (at xGK15[6])
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
		static GKResult IntegrateGK15Impl(Func f, Real a, Real b, IntegrationResult* progress, bool checkFinite) {
			const Real half_length = IntegrationDetail::HalfWidth(a, b);
			const Real center = std::midpoint(a, b);
			int evaluations = 0;
			auto counted = [&](Real x) {
				if (progress) progress->function_evaluations = evaluations;
				return f(x);
			};
			auto sample = [&](Real x) { return IntegrationDetail::Sample(counted, evaluations, checkFinite, x); };
			auto weighted = [&](Real weight, Real value) { return IntegrationDetail::Product({half_length, weight, value}); };

			// Function value at center (used by both Gauss and Kronrod)
			Real f_center = sample(center);

			// Compute sums
			Real result_gauss = weighted(wG7[0], f_center);
			Real result_kronrod = weighted(wGK15[0], f_center);
			Real roundoffFloor = Precision::IntegrationRoundoffFloorScale * std::abs(result_kronrod);

			// Symmetric pairs
			for (int i = 1; i < GK15_N; i++) {
				Real x_scaled = half_length * xGK15[i];
				Real f_plus = sample(center + x_scaled);
				Real f_minus = sample(center - x_scaled);

				// All points contribute to Kronrod
				result_kronrod += weighted(wGK15[i], f_plus) + weighted(wGK15[i], f_minus);
				roundoffFloor += Precision::IntegrationRoundoffFloorScale * std::abs(weighted(wGK15[i], f_plus)) +
					Precision::IntegrationRoundoffFloorScale * std::abs(weighted(wGK15[i], f_minus));

				// Only odd-indexed points (i = 2, 4, 6) contribute to Gauss
				// i=2 uses wG7[1], i=4 uses wG7[2], i=6 uses wG7[3]
				if (i == 2)
					result_gauss += weighted(wG7[1], f_plus) + weighted(wG7[1], f_minus);
				else if (i == 4)
					result_gauss += weighted(wG7[2], f_plus) + weighted(wG7[2], f_minus);
				else if (i == 6)
					result_gauss += weighted(wG7[3], f_plus) + weighted(wG7[3], f_minus);
			}

			// Scale by half-length for change of variables
			auto result = GKResult::Fixed(result_kronrod, std::abs(result_kronrod - result_gauss), evaluations);
			result.roundoff_floor = roundoffFloor;
			return result;
		}

		template<typename Func>
		static GKResult IntegrateGK15Core(Func f, Real a, Real b) {
			return IntegrateGK15Impl(f, a, b, nullptr, true);
		}

		template<typename Func>
		static Real IntegrateGK15(Func f, Real a, Real b)
		{
			return IntegrateGK15Core(f, a, b).ValueOrThrow("IntegrateGK15");
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
		static GKResult IntegrateGK21Impl(Func f, Real a, Real b, IntegrationResult* progress, bool checkFinite) {
			// 21-point Kronrod abscissae (11 unique, from GSL qk21.c)
			// xgk[1], xgk[3], ... are the 10-point Gauss abscissae
			// xgk[0], xgk[2], ... are added by Kronrod extension
			static constexpr Real xgk[11] = {
				MML::Real(0.995657163025808080735527280689003L), // 0 - Kronrod only
				MML::Real(0.973906528517171720077964012084452L), // 1 - Gauss
				MML::Real(0.930157491355708226001207180059508L), // 2 - Kronrod only
				MML::Real(0.865063366688984510732096688423493L), // 3 - Gauss
				MML::Real(0.780817726586416897063717578345042L), // 4 - Kronrod only
				MML::Real(0.679409568299024406234327365114874L), // 5 - Gauss
				MML::Real(0.562757134668604683339000099272694L), // 6 - Kronrod only
				MML::Real(0.433395394129247190799265943165784L), // 7 - Gauss
				MML::Real(0.294392862701460198131126603103866L), // 8 - Kronrod only
				MML::Real(0.148874338981631210884826001129720L), // 9 - Gauss
				0.000000000000000000000000000000000	 // 10 - Center (Kronrod only)
			};

			// 10-point Gauss weights (for odd-indexed Kronrod nodes)
			static constexpr Real wg[5] = {
				MML::Real(0.066671344308688137593568809893332L), // for xgk[1]
				MML::Real(0.149451349150580593145776339657697L), // for xgk[3]
				MML::Real(0.219086362515982043995534934228163L), // for xgk[5]
				MML::Real(0.269266719309996355091226921569469L), // for xgk[7]
				MML::Real(0.295524224714752870173892994651338L)	 // for xgk[9]
			};

			// 21-point Kronrod weights (for all nodes)
			static constexpr Real wgk[11] = {
				MML::Real(0.011694638867371874278064396062192L), // for xgk[0]
				MML::Real(0.032558162307964727478818972459390L), // for xgk[1]
				MML::Real(0.054755896574351996031381300244580L), // for xgk[2]
				MML::Real(0.075039674810919952767043140916190L), // for xgk[3]
				MML::Real(0.093125454583697605535065465083366L), // for xgk[4]
				MML::Real(0.109387158802297641899210590325805L), // for xgk[5]
				MML::Real(0.123491976262065851077958109831074L), // for xgk[6]
				MML::Real(0.134709217311473325928054001771707L), // for xgk[7]
				MML::Real(0.142775938577060080797094273138717L), // for xgk[8]
				MML::Real(0.147739104901338491374841515972068L), // for xgk[9]
				MML::Real(0.149445554002916905664936468389821L)	 // for xgk[10] (center)
			};

			const Real half_length = IntegrationDetail::HalfWidth(a, b);
			const Real center = std::midpoint(a, b);
			int evaluations = 0;
			auto counted = [&](Real x) {
				if (progress) progress->function_evaluations = evaluations;
				return f(x);
			};
			auto sample = [&](Real x) { return IntegrationDetail::Sample(counted, evaluations, checkFinite, x); };
			auto weighted = [&](Real weight, Real value) { return IntegrationDetail::Product({half_length, weight, value}); };

			Real f_center = sample(center);

			// Center point (index 10) - Kronrod only (n=11 is odd, so no Gauss at center)
			Real result_kronrod = weighted(wgk[10], f_center);
			Real result_gauss = 0.0;
			Real roundoffFloor = Precision::IntegrationRoundoffFloorScale * std::abs(result_kronrod);

			// Process symmetric pairs (indices 0-9)
			for (int j = 0; j < 5; j++) {
				// Odd indices (1,3,5,7,9) contribute to both Gauss and Kronrod
				int jtw = j * 2 + 1; // 1, 3, 5, 7, 9
				Real x_scaled = half_length * xgk[jtw];
				Real f_plus = sample(center + x_scaled);
				Real f_minus = sample(center - x_scaled);

				result_gauss += weighted(wg[j], f_plus) + weighted(wg[j], f_minus);
				result_kronrod += weighted(wgk[jtw], f_plus) + weighted(wgk[jtw], f_minus);
				roundoffFloor += Precision::IntegrationRoundoffFloorScale * std::abs(weighted(wgk[jtw], f_plus)) +
					Precision::IntegrationRoundoffFloorScale * std::abs(weighted(wgk[jtw], f_minus));
			}

			for (int j = 0; j < 5; j++) {
				// Even indices (0,2,4,6,8) contribute to Kronrod only
				int jtwm1 = j * 2; // 0, 2, 4, 6, 8
				Real x_scaled = half_length * xgk[jtwm1];
				Real f_plus = sample(center + x_scaled);
				Real f_minus = sample(center - x_scaled);

				result_kronrod += weighted(wgk[jtwm1], f_plus) + weighted(wgk[jtwm1], f_minus);
				roundoffFloor += Precision::IntegrationRoundoffFloorScale * std::abs(weighted(wgk[jtwm1], f_plus)) +
					Precision::IntegrationRoundoffFloorScale * std::abs(weighted(wgk[jtwm1], f_minus));
			}

			auto result = GKResult::Fixed(result_kronrod, std::abs(result_kronrod - result_gauss), evaluations);
			result.roundoff_floor = roundoffFloor;
			return result;
		}

		template<typename Func>
		static GKResult IntegrateGK21Core(Func f, Real a, Real b) {
			return IntegrateGK21Impl(f, a, b, nullptr, true);
		}

		template<typename Func>
		static Real IntegrateGK21(Func f, Real a, Real b)
		{
			return IntegrateGK21Core(f, a, b).ValueOrThrow("IntegrateGK21");
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
		static GKResult IntegrateGK31Impl(Func f, Real a, Real b, IntegrationResult* progress, bool checkFinite) {
			// 31-point Kronrod abscissae (16 unique, from GSL qk31.c)
			// xgk[1], xgk[3], ... are the 15-point Gauss abscissae
			// xgk[0], xgk[2], ... are added by Kronrod extension
			static constexpr Real xgk[16] = {
				MML::Real(0.998002298693397060285172840152271L), // 0 - Kronrod only
				MML::Real(0.987992518020485428489565718586613L), // 1 - Gauss
				MML::Real(0.967739075679139134257347978784337L), // 2 - Kronrod only
				MML::Real(0.937273392400705904307758947710209L), // 3 - Gauss
				MML::Real(0.897264532344081900882509656454496L), // 4 - Kronrod only
				MML::Real(0.848206583410427216200648320774217L), // 5 - Gauss
				MML::Real(0.790418501442465932967649294817947L), // 6 - Kronrod only
				MML::Real(0.724417731360170047416186054613938L), // 7 - Gauss
				MML::Real(0.650996741297416970533735895313275L), // 8 - Kronrod only
				MML::Real(0.570972172608538847537226737253911L), // 9 - Gauss
				MML::Real(0.485081863640239680693655740232351L), // 10 - Kronrod only
				MML::Real(0.394151347077563369897207370981045L), // 11 - Gauss
				MML::Real(0.299180007153168812166780024266389L), // 12 - Kronrod only
				MML::Real(0.201194093997434522300628303394596L), // 13 - Gauss
				MML::Real(0.101142066918717499027074231447392L), // 14 - Kronrod only
				0.000000000000000000000000000000000	 // 15 - Center (Gauss + Kronrod, n=16 even)
			};

			// 15-point Gauss weights (for odd-indexed nodes plus center)
			// wg[0..6] for xgk[1,3,5,7,9,11,13], wg[7] for center xgk[15]
			static constexpr Real wg[8] = {
				MML::Real(0.030753241996117268354628393577204L), // for xgk[1]
				MML::Real(0.070366047488108124709267416450667L), // for xgk[3]
				MML::Real(0.107159220467171935011869546685869L), // for xgk[5]
				MML::Real(0.139570677926154314447804794511028L), // for xgk[7]
				MML::Real(0.166269205816993933553200860481209L), // for xgk[9]
				MML::Real(0.186161000015562211026800561866423L), // for xgk[11]
				MML::Real(0.198431485327111576456118326443839L), // for xgk[13]
				MML::Real(0.202578241925561272880620199967519L)	 // for xgk[15] (center)
			};

			// 31-point Kronrod weights (for all nodes)
			static constexpr Real wgk[16] = {
				MML::Real(0.005377479872923348987792051430128L), // for xgk[0]
				MML::Real(0.015007947329316122538374763075807L), // for xgk[1]
				MML::Real(0.025460847326715320186874001019653L), // for xgk[2]
				MML::Real(0.035346360791375846222037948478360L), // for xgk[3]
				MML::Real(0.044589751324764876608227299373280L), // for xgk[4]
				MML::Real(0.053481524690928087265343147239430L), // for xgk[5]
				MML::Real(0.062009567800670640285139230960803L), // for xgk[6]
				MML::Real(0.069854121318728258709520077099147L), // for xgk[7]
				MML::Real(0.076849680757720378894432777482659L), // for xgk[8]
				MML::Real(0.083080502823133021038289247286104L), // for xgk[9]
				MML::Real(0.088564443056211770647275443693774L), // for xgk[10]
				MML::Real(0.093126598170825321225486872747346L), // for xgk[11]
				MML::Real(0.096642726983623678505179907627589L), // for xgk[12]
				MML::Real(0.099173598721791959332393173484603L), // for xgk[13]
				MML::Real(0.100769845523875595044946662617570L), // for xgk[14]
				MML::Real(0.101330007014791549017374792767493L)	 // for xgk[15] (center)
			};

			const Real half_length = IntegrationDetail::HalfWidth(a, b);
			const Real center = std::midpoint(a, b);
			int evaluations = 0;
			auto counted = [&](Real x) {
				if (progress) progress->function_evaluations = evaluations;
				return f(x);
			};
			auto sample = [&](Real x) { return IntegrationDetail::Sample(counted, evaluations, checkFinite, x); };
			auto weighted = [&](Real weight, Real value) { return IntegrationDetail::Product({half_length, weight, value}); };

			Real f_center = sample(center);

			// Center point (index 15) - Both Gauss and Kronrod (n=16 is even)
			Real result_kronrod = weighted(wgk[15], f_center);
			Real result_gauss = weighted(wg[7], f_center);
			Real roundoffFloor = Precision::IntegrationRoundoffFloorScale * std::abs(result_kronrod);

			// Process symmetric pairs (indices 0-14)
			for (int j = 0; j < 7; j++) {
				// Odd indices (1,3,5,7,9,11,13) contribute to both Gauss and Kronrod
				int jtw = j * 2 + 1; // 1, 3, 5, 7, 9, 11, 13
				Real x_scaled = half_length * xgk[jtw];
				Real f_plus = sample(center + x_scaled);
				Real f_minus = sample(center - x_scaled);

				result_gauss += weighted(wg[j], f_plus) + weighted(wg[j], f_minus);
				result_kronrod += weighted(wgk[jtw], f_plus) + weighted(wgk[jtw], f_minus);
				roundoffFloor += Precision::IntegrationRoundoffFloorScale * std::abs(weighted(wgk[jtw], f_plus)) +
					Precision::IntegrationRoundoffFloorScale * std::abs(weighted(wgk[jtw], f_minus));
			}

			for (int j = 0; j < 8; j++) {
				// Even indices (0,2,4,6,8,10,12,14) contribute to Kronrod only
				int jtwm1 = j * 2; // 0, 2, 4, 6, 8, 10, 12, 14
				Real x_scaled = half_length * xgk[jtwm1];
				Real f_plus = sample(center + x_scaled);
				Real f_minus = sample(center - x_scaled);

				result_kronrod += weighted(wgk[jtwm1], f_plus) + weighted(wgk[jtwm1], f_minus);
				roundoffFloor += Precision::IntegrationRoundoffFloorScale * std::abs(weighted(wgk[jtwm1], f_plus)) +
					Precision::IntegrationRoundoffFloorScale * std::abs(weighted(wgk[jtwm1], f_minus));
			}

			auto result = GKResult::Fixed(result_kronrod, std::abs(result_kronrod - result_gauss), evaluations);
			result.roundoff_floor = roundoffFloor;
			return result;
		}

		template<typename Func>
		static GKResult IntegrateGK31Core(Func f, Real a, Real b) {
			return IntegrateGK31Impl(f, a, b, nullptr, true);
		}

		template<typename Func>
		static Real IntegrateGK31(Func f, Real a, Real b)
		{
			return IntegrateGK31Core(f, a, b).ValueOrThrow("IntegrateGK31");
		}

		/*************************************************************************/
		/*****              ADAPTIVE GAUSS-KRONROD INTEGRATION               *****/
		/*****  Recursively subdivides until error tolerance met             *****/
		/*************************************************************************/

		/// @brief Rule selection for adaptive Gauss-Kronrod
		enum class GKRule { GK15, GK21, GK31 };

		namespace GKDetail {
			template<class Function>
			void Adaptive(Function& function, Real a, Real b, Real absolute, Real relative, int depth,
				GKRule rule, IntegrationResult& result, bool checkFinite = true) {
				IntegrationDetail::HalfWidth(a, b);
				if (rule != GKRule::GK15 && rule != GKRule::GK21 && rule != GKRule::GK31)
					throw ArgumentError("Adaptive Gauss-Kronrod: invalid rule");
				using Interval = std::array<Real, 2>;
				auto sample = [&](Real x) {
					return IntegrationDetail::Sample(function, result.function_evaluations, checkFinite, x);
				};
				auto evaluate = [&](const Interval& interval) -> IntegrationResult {
					GKResult local;
					switch (rule) {
					case GKRule::GK15: local = IntegrateGK15Impl(sample, interval[0], interval[1], nullptr, false); break;
					case GKRule::GK21: local = IntegrateGK21Impl(sample, interval[0], interval[1], nullptr, false); break;
					case GKRule::GK31: local = IntegrateGK31Impl(sample, interval[0], interval[1], nullptr, false); break;
					}
					local.error_estimate = std::max(local.error_estimate, local.roundoff_floor);
					return local;
				};
				auto split = [](const Interval& interval) {
					const Real midpoint = std::midpoint(interval[0], interval[1]);
					if (midpoint == interval[0] || midpoint == interval[1]) return std::vector<Interval>{};
					return std::vector<Interval>{{interval[0], midpoint}, {midpoint, interval[1]}};
				};
				const int cost = rule == GKRule::GK15 ? 15 : rule == GKRule::GK21 ? 21 : 31;
				IntegrationDetail::RefineCells(result, std::vector<Interval>{{a, b}}, evaluate, split,
					absolute, relative, depth, std::numeric_limits<int>::max(), cost);
			}
		}

		/// @brief Adaptive Gauss-Kronrod integration with automatic subdivision.
		/// @details Bisects the largest refinable error and checks the aggregate frontier.
		///          Includes parent work, arithmetic floors and stagnation diagnostics.
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
		static GKResult IntegrateGKAdaptiveCore(Func f, Real a, Real b, Real tol_abs = PrecisionValues<Real>::IntegrationTolerance, Real tol_rel = PrecisionValues<Real>::IntegrationTolerance, int max_depth = 50,
																				GKRule rule = GKRule::GK15) {
			GKResult result;
			GKDetail::Adaptive(f, a, b, tol_abs, tol_rel, max_depth, rule, result);
			return result;
		}

		/// Result-returning variant with the plain argument list (RESULT_POLICY §4).
		template<typename Func>
		static auto IntegrateGKAdaptiveDetailed(Func f, Real a, Real b, Real tol_abs, Real tol_rel = PrecisionValues<Real>::IntegrationTolerance, int max_depth = 50, GKRule rule = GKRule::GK15)
		{
			auto r = IntegrateGKAdaptiveCore(f, a, b, tol_abs, tol_rel, max_depth, rule);
			r.algorithm_name = "IntegrateGKAdaptive";
			return r;
		}

		template<typename Func>
		static Real IntegrateGKAdaptive(Func f, Real a, Real b, Real tol_abs = PrecisionValues<Real>::IntegrationTolerance, Real tol_rel = PrecisionValues<Real>::IntegrationTolerance, int max_depth = 50,
																				GKRule rule = GKRule::GK15)
		{
			return IntegrateGKAdaptiveCore(f, a, b, tol_abs, tol_rel, max_depth, rule).ValueOrThrow("IntegrateGKAdaptive");
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
			return IntegrateGK15Core(f, a, b).value;
		}

		/// @brief Integrate using G10K21 rule (non-adaptive, single application).
		template<typename Func>
		static Real IntegrateGaussKronrod21(Func f, Real a, Real b) {
			return IntegrateGK21Core(f, a, b).value;
		}

		/// @brief Integrate using G15K31 rule (non-adaptive, single application).
		template<typename Func>
		static Real IntegrateGaussKronrod31(Func f, Real a, Real b) {
			return IntegrateGK31Core(f, a, b).value;
		}

		/// @brief Adaptive Gauss-Kronrod integration (recommended general use).
		/// @param f Function to integrate
		/// @param a Lower limit
		/// @param b Upper limit
		/// @param tol Absolute tolerance (default 1e-10)
		/// @return Integral value
		template<typename Func>
		static Real IntegrateGaussKronrod(Func f, Real a, Real b, Real tol = Precision::IntegrationTolerance) {
			return IntegrateGKAdaptiveCore(f, a, b, tol, tol, 50, GKRule::GK15).ValueOrThrow("IntegrateGaussKronrod");
		}

		/// @brief Adaptive Gauss-Kronrod with error estimate.
		/// @param f Function to integrate
		/// @param a Lower limit
		/// @param b Upper limit
		/// @param error Output parameter for error estimate
		/// @param tol Tolerance
		/// @return Integral value
		template<typename Func>
		static Real IntegrateGaussKronrod(Func f, Real a, Real b, Real* error, Real tol = Precision::IntegrationTolerance) {
			auto result = IntegrateGKAdaptiveCore(f, a, b, tol, tol, 50, GKRule::GK15);
			if (error)
				*error = result.error_estimate;
			return result.ValueOrThrow("IntegrateGaussKronrod");
		}

		/*****************************************************************/
		/*****        Detailed API - Gauss-Kronrod                  *****/
		/*****************************************************************/

		/// Helper to populate IntegrationDetailedResult from GKResult
		namespace GKDetail
		{
			inline void PopulateFromGK(IntegrationResult& out, const GKResult& gk, const char* algorithm_name)
			{
				IntegrationDetail::Adopt(out, gk, algorithm_name);
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
					auto gk = IntegrateGK15Impl(f, a, b, &result, config.check_finite);
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
					auto gk = IntegrateGK21Impl(f, a, b, &result, config.check_finite);
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
					auto gk = IntegrateGK31Impl(f, a, b, &result, config.check_finite);
					GKDetail::PopulateFromGK(result, gk, "IntegrateGK31");
				});
		}

		/// Adaptive Gauss-Kronrod with full diagnostics
		template<typename Func>
		static IntegrationDetailedResult IntegrateGKAdaptiveDetailed(
			Func f, Real a, Real b,
			const IntegrationConfig& config = {},
			Real tol_abs = PrecisionValues<Real>::IntegrationTolerance, Real tol_rel = PrecisionValues<Real>::IntegrationTolerance,
			int max_depth = 50, GKRule rule = GKRule::GK15)
		{
			return IntegrationDetail::ExecuteIntegrationDetailed<IntegrationDetailedResult>(
				"IntegrateGKAdaptive", config,
				[&](IntegrationDetailedResult& result) {
					GKDetail::Adaptive(f, a, b, tol_abs, tol_rel, max_depth, rule, result, config.check_finite);
				});
		}

	} // namespace Integration
} // namespace MML

#endif // MML_GAUSS_KRONROD_H
