#if !defined MML_LINEAR_ALG_EQ_SYSTEMS_SYM_DEFS_H
#define MML_LINEAR_ALG_EQ_SYSTEMS_SYM_DEFS_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Matrix/MatrixSym.h"
#endif

// SIOF-SAFE: All Matrix/Vector types use construct-on-first-use idiom.
// Scalar primitives (Real, int) use inline constexpr - safe by design.

namespace MML::TestBeds
{
	/***********************************************************************************************/
	/**************                 3 x 3 symmetric test matrices                           ********/
	/***********************************************************************************************/
	inline const MatrixSym<Real>& symm_mat_3x3() {
		static const MatrixSym<Real> m{ 3, { 7.9, -9.3,  1.7,
		                                               -5.9,  9.9,
		                                                     -5.5 } };
		return m;
	}

	inline const Vector<Real>& symm_mat_3x3_rhs0() {
		static const Vector<Real> v{ -2.7, 1.3, 5.4 };
		return v;
	}
	inline const Vector<Real>& symm_mat_3x3_rhs0_sol() {
		static const Vector<Real> v{ -0.45139279169212, 0.103906536346976, -0.310565058032987 };
		return v;
	}

	inline const Vector<Real>& symm_mat_3x3_eigen_val() {
		static const Vector<Real> v{ 19.2792340213883, -2.727695435303, -12.4515385860853 };
		return v;
	}
	inline const std::vector<Vector<Real>>& symm_mat_3x3_eigen_vecs() {
		static const std::vector<Vector<Real>> vecs
		{
			Vector<Real>{ -0.6953039682,  0.5949327881,  0.4032460409 },
			Vector<Real>{ -0.7178924215, -0.5480450913, -0.4292750273 },
			Vector<Real>{ 0.03439277556,  0.5879639066, -0.8081556666 },
		};
		return vecs;
	}

	// Matrix properties for symm_mat_3x3 - scalar primitives are SIOF-safe
	inline constexpr Real symm_mat_3x3_det = 658.7220000000001;
	inline const Vector<Real>& symm_mat_3x3_singular_values() {
		static const Vector<Real> v{19.279234021388303, 12.451538586085298, 2.727695435303002};
		return v;
	}
	inline constexpr Real symm_mat_3x3_cond_2 = 7.068434587696746;
	inline constexpr int  symm_mat_3x3_rank = 3;

	inline const MatrixSym<Real>& symm_mat_3x3_1() {
		static const MatrixSym<Real> m{ 3, { 2, -1,  0,
		                                              2, -1,
		                                                  2 } };
		return m;
	}

	inline const Vector<Real>& symm_mat_3x3_1_rhs0() {
		static const Vector<Real> v{ -2.7, 1.3, 5.4 };
		return v;
	}
	inline const Vector<Real>& symm_mat_3x3_1_rhs0_sol() {
		static const Vector<Real> v{ -0.45139279169212, 0.103906536346976, -0.310565058032987 };
		return v;
	}

	// Matrix properties for symm_mat_3x3_1 - scalar primitives are SIOF-safe
	inline constexpr Real symm_mat_3x3_1_det = 4.0;
	inline const Vector<Real>& symm_mat_3x3_1_singular_values() {
		static const Vector<Real> v{3.414213562373095, 2.0, 0.5857864376269049};
		return v;
	}
	inline constexpr Real symm_mat_3x3_1_cond_2 = 5.828427124746191;
	inline constexpr int  symm_mat_3x3_1_rank = 3;

	/***********************************************************************************************/
	/**************                 5 x 5 symmetric test matrices                           ********/
	/***********************************************************************************************/
	inline const MatrixSym<Real>& symm_mat_5x5() {
		static const MatrixSym<Real> m{ 5,  {1.4, 2.1, 2.1, 7.4, 9.6,
		                                              1.5, 1.1, 0.7, 5.0,
		                                                   9.6, 5.4, 8.8,
		                                                        0.4, 8.0,
		                                                             7.7 } };
		return m;
	}

	inline const Vector<Real>& symm_mat_5x5_rhs0() {
		static const Vector<Real> v{ -2.7, 1.3, 5.4, 1.5, -2.7 };
		return v;
	}
	inline const Vector<Real>& symm_mat_5x5_rhs0_sol() {
		static const Vector<Real> v{ 6.33596564609689, -4.17485303537051, -0.47300554791769, -0.152420675302137, 0.160137642067717 };
		return v;
	}

	inline const Vector<Real>& symm_mat_5x5_eigen_val() {
		static const Vector<Real> v{ 24.5884918308427, 8.00842513434593, 3.5416658121729, -0.570094077835055, -13.2684886995265 };
		return v;
	}
	inline const std::vector<Vector<Real>>& symm_mat_5x5_eigen_vecs() {
		static const std::vector<Vector<Real>> vecs
		{
			Vector<Real>{ 0.3207830408,  0.4359308891,  0.3974911379,  0.4784849759,  0.5657874369 },
			Vector<Real>{ 0.323313998,   0.4160284348,  0.4267540611, -0.7051973277, -0.2072826582 },
			Vector<Real>{ 0.1398659018, -0.2065443391,  0.5233038154,  0.4572026692, -0.6744596367 },
			Vector<Real>{ 0.8045643778, -0.5452049336, -0.1559364946, -0.09146643767, 0.1508159234 },
			Vector<Real>{ 0.3545137387,  0.5449546933, -0.6014305388,  0.2373783201, -0.3990955168 },
		};
		return vecs;
	}

	// Matrix properties for symm_mat_5x5 - scalar primitives are SIOF-safe
	inline constexpr Real symm_mat_5x5_det = 3958.6609039999995;
	inline const Vector<Real>& symm_mat_5x5_singular_values() {
		static const Vector<Real> v{24.588491830842698, 13.268488699526523, 8.00842513434593, 3.5416658121729, 0.5700940778350548};
		return v;
	}
	inline constexpr Real symm_mat_5x5_cond_2 = 43.13301067461779;
	inline constexpr int  symm_mat_5x5_rank = 5;

	/***********************************************************************************************/
	/*************              10 x 10 symmetric test matrices                         ************/
	/***********************************************************************************************/
	inline const MatrixSym<Real>& symm_mat_10x10() {
		static const MatrixSym<Real> m{ 10, {  -0.91, -0.88, -0.32, -0.96, -0.27,  0.5,  -0.81,  0.65, -0.87, -0.37,
		                                             -0.73, -0.072, 0.92,  0.11, -0.25,  1.0,   0.3,  -0.24,  0.71,
		                                                    -0.85, -0.18, -0.78,  0.29, -0.05, -0.47,  0.57,  0.64,
		                                                           -0.27, -0.15, -0.85, -0.41,  0.29,  0.073, 0.43,
		                                                                  -0.63, -0.41,  0.37,  0.36,  0.58,  0.71,
		                                                                          0.52,  0.84,  0.66, -0.14,  0.29,
		                                                                                 0.82,  0.27,  0.14,  0.79,
		                                                                                        0.47, -0.89,  0.11,
		                                                                                               0.42,  0.81,
		                                                                                                     -0.32 } };
		return m;
	}

	inline const Vector<Real>& symm_mat_10x10_rhs0() {
		static const Vector<Real> v{ -0.57, 0.92, 0.53, -0.8, -1.4, 0.8, 0.72, -1.5, 1.2, -0.9 };
		return v;
	}
	inline const Vector<Real>& symm_mat_10x10_rhs0_sol() {
		static const Vector<Real> v{ 3.52975177499083, 1.14755244118994, -3.37526055097626, -3.58478581007823, 0.0630707984729596, -1.85597545512844, 4.29069683741684, -5.26485523582046, 0.0269309930785835, 2.38666042289625 };
		return v;
	}

	inline const Vector<Real>& symm_mat_10x10_eigen_val() {
		static const Vector<Real> v{ 2.46552276858711, 2.07894422432798, 1.12637210737141, 0.877705460150441, 0.069595393281048, -0.289946145524018, -1.31538396615402, -1.68322193756041, -1.77811740978995, -3.79147049468958 };
		return v;
	}
	inline const std::vector<Vector<Real>>& symm_mat_10x10_eigen_vecs() {
		static const std::vector<Vector<Real>> vecs
		{
			Vector<Real>{  0.3660146,     -0.1394605012, -0.5662781991,   0.04407610321, -0.4997210282,  0.2820569235, -0.2654623887,  0.185069201,    -0.2576078803,  -0.1536569441 },
			Vector<Real>{ -0.0786040476,  -0.3145579783,  0.1993343822,  -0.4630688613,   0.0137507947, -0.4303463398, -0.2145106699, -0.003346098308, -0.588877082,   -0.2500368327 },
			Vector<Real>{  0.186581028,   -0.2850284346, -0.03055385745,  0.2087237244,   0.214203072,  -0.3358183881, -0.3951749565,  0.41632425,      0.09657341182,  0.5847860147 },
			Vector<Real>{ -0.442592973,    0.461657522,  -0.3760287681,   0.4244650158,   0.0340738307, -0.3407530939, -0.0955367936, -0.054844722,    -0.3734622731,   0.02307855419 },
			Vector<Real>{ -0.4032869848,  -0.1886901214,  0.2446369709,   0.3029474893,   0.1245023272,  0.3438925557, -0.1216928937,  0.5841845899,   -0.03900630761, -0.3984545273 },
			Vector<Real>{ -0.04862557386, -0.295455427,  -0.4920414413,  -0.1418043838,   0.3507408272, -0.04655910028, 0.6679025287,  0.2548089709,   -0.1053825257,   0.02854167303 },
			Vector<Real>{ -0.106837945,    0.3104465485,  0.2774032194,  -0.2263131742,  -0.4500895918,  0.1000492797,  0.3064679175,  0.4209284903,   -0.2705969647,   0.4551063616 },
			Vector<Real>{  0.2961310056,   0.531108144,  -0.06145183558, -0.3250429634,   0.5651514285,  0.2740815745, -0.2360767942,  0.1976298839,   -0.1718546758,  -0.04427878569 },
			Vector<Real>{  0.1887690916,  -0.2204005835,  0.2333845627,   0.3990940422,   0.2008414292,  0.3745029919,  0.1232095207, -0.3413452318,   -0.5660476865,   0.2633410279 },
			Vector<Real>{ -0.5739400396,  -0.1931526915, -0.2451823337,  -0.3631199308,   0.04702991421, 0.4006823686, -0.3041151033, -0.2246111993,    0.04762110897,  0.3651743111 },
		};
		return vecs;
	}

	// Matrix properties for symm_mat_10x10 - scalar primitives are SIOF-safe
	inline constexpr Real symm_mat_10x10_det = -0.31016234912654906;
	inline const Vector<Real>& symm_mat_10x10_singular_values() {
		static const Vector<Real> v{3.79147049468958, 2.46552276858711, 2.07894422432798, 1.7781174097899468, 1.6832219375604082, 1.3153839661540222, 1.1263721073714085, 0.8777054601504407, 0.28994614552401766, 0.06959539328104801};
		return v;
	}
	inline constexpr Real symm_mat_10x10_cond_2 = 54.47879912099627;
	inline constexpr int  symm_mat_10x10_rank = 10;
}

#endif // MML_LINEAR_ALG_EQ_SYSTEMS_SYM_DEFS_H

