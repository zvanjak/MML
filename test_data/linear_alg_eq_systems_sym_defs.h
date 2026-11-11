#if !defined MML_LINEAR_ALG_EQ_SYSTEMS_SYM_DEFS_H
#define MML_LINEAR_ALG_EQ_SYSTEMS_SYM_DEFS_H

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Matrix/MatrixSym.h>
#endif

// SIOF-SAFE: All Matrix/Vector types use construct-on-first-use idiom.
// Scalar primitives (Real, int) use inline constexpr - safe by design.

namespace MML::TestBeds
{
	/***********************************************************************************************/
	/**************                 3 x 3 symmetric test matrices                           ********/
	/***********************************************************************************************/
	inline const MatrixSym<Real>& symm_mat_3x3() {
		static const MatrixSym<Real> m{ 3, { MML::Real(7.9), -MML::Real(9.3),  MML::Real(1.7),
		                                               -MML::Real(5.9),  MML::Real(9.9),
		                                                     -5.5 } };
		return m;
	}

	inline const Vector<Real>& symm_mat_3x3_rhs0() {
		static const Vector<Real> v{ -MML::Real(2.7), MML::Real(1.3), MML::Real(5.4) };
		return v;
	}
	inline const Vector<Real>& symm_mat_3x3_rhs0_sol() {
		static const Vector<Real> v{ -MML::Real(0.45139279169212), MML::Real(0.103906536346976), -MML::Real(0.310565058032987) };
		return v;
	}

	inline const Vector<Real>& symm_mat_3x3_eigen_val() {
		static const Vector<Real> v{ MML::Real(19.2792340213883), -MML::Real(2.727695435303), -MML::Real(12.4515385860853) };
		return v;
	}
	inline const std::vector<Vector<Real>>& symm_mat_3x3_eigen_vecs() {
		static const std::vector<Vector<Real>> vecs
		{
			Vector<Real>{ -MML::Real(0.6953039682),  MML::Real(0.5949327881),  MML::Real(0.4032460409) },
			Vector<Real>{ -MML::Real(0.7178924215), -MML::Real(0.5480450913), -MML::Real(0.4292750273) },
			Vector<Real>{ MML::Real(0.03439277556),  MML::Real(0.5879639066), -MML::Real(0.8081556666) },
		};
		return vecs;
	}

	// Matrix properties for symm_mat_3x3 - scalar primitives are SIOF-safe
	inline constexpr Real symm_mat_3x3_det = MML::Real(658.7220000000001);
	inline const Vector<Real>& symm_mat_3x3_singular_values() {
		static const Vector<Real> v{MML::Real(19.279234021388303), MML::Real(12.451538586085298), MML::Real(2.727695435303002)};
		return v;
	}
	inline constexpr Real symm_mat_3x3_cond_2 = MML::Real(7.068434587696746);
	inline constexpr int  symm_mat_3x3_rank = 3;

	inline const MatrixSym<Real>& symm_mat_3x3_1() {
		static const MatrixSym<Real> m{ 3, { 2, -1,  0,
		                                              2, -1,
		                                                  2 } };
		return m;
	}

	inline const Vector<Real>& symm_mat_3x3_1_rhs0() {
		static const Vector<Real> v{ -MML::Real(2.7), MML::Real(1.3), MML::Real(5.4) };
		return v;
	}
	inline const Vector<Real>& symm_mat_3x3_1_rhs0_sol() {
		static const Vector<Real> v{ -MML::Real(0.45139279169212), MML::Real(0.103906536346976), -MML::Real(0.310565058032987) };
		return v;
	}

	// Matrix properties for symm_mat_3x3_1 - scalar primitives are SIOF-safe
	inline constexpr Real symm_mat_3x3_1_det = 4.0;
	inline const Vector<Real>& symm_mat_3x3_1_singular_values() {
		static const Vector<Real> v{MML::Real(3.414213562373095), MML::Real(2.0), MML::Real(0.5857864376269049)};
		return v;
	}
	inline constexpr Real symm_mat_3x3_1_cond_2 = MML::Real(5.828427124746191);
	inline constexpr int  symm_mat_3x3_1_rank = 3;

	/***********************************************************************************************/
	/**************                 5 x 5 symmetric test matrices                           ********/
	/***********************************************************************************************/
	inline const MatrixSym<Real>& symm_mat_5x5() {
		static const MatrixSym<Real> m{ 5,  {MML::Real(1.4), MML::Real(2.1), MML::Real(2.1), MML::Real(7.4), MML::Real(9.6),
		                                              MML::Real(1.5), MML::Real(1.1), MML::Real(0.7), MML::Real(5.0),
		                                                   MML::Real(9.6), MML::Real(5.4), MML::Real(8.8),
		                                                        MML::Real(0.4), MML::Real(8.0),
		                                                             MML::Real(7.7) } };
		return m;
	}

	inline const Vector<Real>& symm_mat_5x5_rhs0() {
		static const Vector<Real> v{ -MML::Real(2.7), MML::Real(1.3), MML::Real(5.4), MML::Real(1.5), -MML::Real(2.7) };
		return v;
	}
	inline const Vector<Real>& symm_mat_5x5_rhs0_sol() {
		static const Vector<Real> v{ MML::Real(6.33596564609689), -MML::Real(4.17485303537051), -MML::Real(0.47300554791769), -MML::Real(0.152420675302137), MML::Real(0.160137642067717) };
		return v;
	}

	inline const Vector<Real>& symm_mat_5x5_eigen_val() {
		static const Vector<Real> v{ MML::Real(24.5884918308427), MML::Real(8.00842513434593), MML::Real(3.5416658121729), -MML::Real(0.570094077835055), -MML::Real(13.2684886995265) };
		return v;
	}
	inline const std::vector<Vector<Real>>& symm_mat_5x5_eigen_vecs() {
		static const std::vector<Vector<Real>> vecs
		{
			Vector<Real>{ MML::Real(0.3207830408),  MML::Real(0.4359308891),  MML::Real(0.3974911379),  MML::Real(0.4784849759),  MML::Real(0.5657874369) },
			Vector<Real>{ MML::Real(0.323313998),   MML::Real(0.4160284348),  MML::Real(0.4267540611), -MML::Real(0.7051973277), -MML::Real(0.2072826582) },
			Vector<Real>{ MML::Real(0.1398659018), -MML::Real(0.2065443391),  MML::Real(0.5233038154),  MML::Real(0.4572026692), -MML::Real(0.6744596367) },
			Vector<Real>{ MML::Real(0.8045643778), -MML::Real(0.5452049336), -MML::Real(0.1559364946), -MML::Real(0.09146643767), MML::Real(0.1508159234) },
			Vector<Real>{ MML::Real(0.3545137387),  MML::Real(0.5449546933), -MML::Real(0.6014305388),  MML::Real(0.2373783201), -MML::Real(0.3990955168) },
		};
		return vecs;
	}

	// Matrix properties for symm_mat_5x5 - scalar primitives are SIOF-safe
	inline constexpr Real symm_mat_5x5_det = MML::Real(3958.6609039999995);
	inline const Vector<Real>& symm_mat_5x5_singular_values() {
		static const Vector<Real> v{MML::Real(24.588491830842698), MML::Real(13.268488699526523), MML::Real(8.00842513434593), MML::Real(3.5416658121729), MML::Real(0.5700940778350548)};
		return v;
	}
	inline constexpr Real symm_mat_5x5_cond_2 = MML::Real(43.13301067461779);
	inline constexpr int  symm_mat_5x5_rank = 5;

	/***********************************************************************************************/
	/*************              10 x 10 symmetric test matrices                         ************/
	/***********************************************************************************************/
	inline const MatrixSym<Real>& symm_mat_10x10() {
		static const MatrixSym<Real> m{ 10, {  -MML::Real(0.91), -MML::Real(0.88), -MML::Real(0.32), -MML::Real(0.96), -MML::Real(0.27),  MML::Real(0.5),  -MML::Real(0.81),  MML::Real(0.65), -MML::Real(0.87), -MML::Real(0.37),
		                                             -MML::Real(0.73), -MML::Real(0.072), MML::Real(0.92),  MML::Real(0.11), -MML::Real(0.25),  MML::Real(1.0),   MML::Real(0.3),  -MML::Real(0.24),  MML::Real(0.71),
		                                                    -MML::Real(0.85), -MML::Real(0.18), -MML::Real(0.78),  MML::Real(0.29), -MML::Real(0.05), -MML::Real(0.47),  MML::Real(0.57),  MML::Real(0.64),
		                                                           -MML::Real(0.27), -MML::Real(0.15), -MML::Real(0.85), -MML::Real(0.41),  MML::Real(0.29),  MML::Real(0.073), MML::Real(0.43),
		                                                                  -MML::Real(0.63), -MML::Real(0.41),  MML::Real(0.37),  MML::Real(0.36),  MML::Real(0.58),  MML::Real(0.71),
		                                                                          MML::Real(0.52),  MML::Real(0.84),  MML::Real(0.66), -MML::Real(0.14),  MML::Real(0.29),
		                                                                                 MML::Real(0.82),  MML::Real(0.27),  MML::Real(0.14),  MML::Real(0.79),
		                                                                                        MML::Real(0.47), -MML::Real(0.89),  MML::Real(0.11),
		                                                                                               MML::Real(0.42),  MML::Real(0.81),
		                                                                                                     -MML::Real(0.32) } };
		return m;
	}

	inline const Vector<Real>& symm_mat_10x10_rhs0() {
		static const Vector<Real> v{ -MML::Real(0.57), MML::Real(0.92), MML::Real(0.53), -MML::Real(0.8), -MML::Real(1.4), MML::Real(0.8), MML::Real(0.72), -MML::Real(1.5), MML::Real(1.2), -MML::Real(0.9) };
		return v;
	}
	inline const Vector<Real>& symm_mat_10x10_rhs0_sol() {
		static const Vector<Real> v{ MML::Real(3.52975177499083), MML::Real(1.14755244118994), -MML::Real(3.37526055097626), -MML::Real(3.58478581007823), MML::Real(0.0630707984729596), -MML::Real(1.85597545512844), MML::Real(4.29069683741684), -MML::Real(5.26485523582046), MML::Real(0.0269309930785835), MML::Real(2.38666042289625) };
		return v;
	}

	inline const Vector<Real>& symm_mat_10x10_eigen_val() {
		static const Vector<Real> v{ MML::Real(2.46552276858711), MML::Real(2.07894422432798), MML::Real(1.12637210737141), MML::Real(0.877705460150441), MML::Real(0.069595393281048), -MML::Real(0.289946145524018), -MML::Real(1.31538396615402), -MML::Real(1.68322193756041), -MML::Real(1.77811740978995), -MML::Real(3.79147049468958) };
		return v;
	}
	inline const std::vector<Vector<Real>>& symm_mat_10x10_eigen_vecs() {
		static const std::vector<Vector<Real>> vecs
		{
			Vector<Real>{  MML::Real(0.3660146),     -MML::Real(0.1394605012), -MML::Real(0.5662781991),   MML::Real(0.04407610321), -MML::Real(0.4997210282),  MML::Real(0.2820569235), -MML::Real(0.2654623887),  MML::Real(0.185069201),    -MML::Real(0.2576078803),  -MML::Real(0.1536569441) },
			Vector<Real>{ -MML::Real(0.0786040476),  -MML::Real(0.3145579783),  MML::Real(0.1993343822),  -MML::Real(0.4630688613),   MML::Real(0.0137507947), -MML::Real(0.4303463398), -MML::Real(0.2145106699), -MML::Real(0.003346098308), -MML::Real(0.588877082),   -MML::Real(0.2500368327) },
			Vector<Real>{  MML::Real(0.186581028),   -MML::Real(0.2850284346), -MML::Real(0.03055385745),  MML::Real(0.2087237244),   MML::Real(0.214203072),  -MML::Real(0.3358183881), -MML::Real(0.3951749565),  MML::Real(0.41632425),      MML::Real(0.09657341182),  MML::Real(0.5847860147) },
			Vector<Real>{ -MML::Real(0.442592973),    MML::Real(0.461657522),  -MML::Real(0.3760287681),   MML::Real(0.4244650158),   MML::Real(0.0340738307), -MML::Real(0.3407530939), -MML::Real(0.0955367936), -MML::Real(0.054844722),    -MML::Real(0.3734622731),   MML::Real(0.02307855419) },
			Vector<Real>{ -MML::Real(0.4032869848),  -MML::Real(0.1886901214),  MML::Real(0.2446369709),   MML::Real(0.3029474893),   MML::Real(0.1245023272),  MML::Real(0.3438925557), -MML::Real(0.1216928937),  MML::Real(0.5841845899),   -MML::Real(0.03900630761), -MML::Real(0.3984545273) },
			Vector<Real>{ -MML::Real(0.04862557386), -MML::Real(0.295455427),  -MML::Real(0.4920414413),  -MML::Real(0.1418043838),   MML::Real(0.3507408272), -MML::Real(0.04655910028), MML::Real(0.6679025287),  MML::Real(0.2548089709),   -MML::Real(0.1053825257),   MML::Real(0.02854167303) },
			Vector<Real>{ -MML::Real(0.106837945),    MML::Real(0.3104465485),  MML::Real(0.2774032194),  -MML::Real(0.2263131742),  -MML::Real(0.4500895918),  MML::Real(0.1000492797),  MML::Real(0.3064679175),  MML::Real(0.4209284903),   -MML::Real(0.2705969647),   MML::Real(0.4551063616) },
			Vector<Real>{  MML::Real(0.2961310056),   MML::Real(0.531108144),  -MML::Real(0.06145183558), -MML::Real(0.3250429634),   MML::Real(0.5651514285),  MML::Real(0.2740815745), -MML::Real(0.2360767942),  MML::Real(0.1976298839),   -MML::Real(0.1718546758),  -MML::Real(0.04427878569) },
			Vector<Real>{  MML::Real(0.1887690916),  -MML::Real(0.2204005835),  MML::Real(0.2333845627),   MML::Real(0.3990940422),   MML::Real(0.2008414292),  MML::Real(0.3745029919),  MML::Real(0.1232095207), -MML::Real(0.3413452318),   -MML::Real(0.5660476865),   MML::Real(0.2633410279) },
			Vector<Real>{ -MML::Real(0.5739400396),  -MML::Real(0.1931526915), -MML::Real(0.2451823337),  -MML::Real(0.3631199308),   MML::Real(0.04702991421), MML::Real(0.4006823686), -MML::Real(0.3041151033), -MML::Real(0.2246111993),    MML::Real(0.04762110897),  MML::Real(0.3651743111) },
		};
		return vecs;
	}

	// Matrix properties for symm_mat_10x10 - scalar primitives are SIOF-safe
	inline constexpr Real symm_mat_10x10_det = -MML::Real(0.31016234912654906);
	inline const Vector<Real>& symm_mat_10x10_singular_values() {
		static const Vector<Real> v{MML::Real(3.79147049468958), MML::Real(2.46552276858711), MML::Real(2.07894422432798), MML::Real(1.7781174097899468), MML::Real(1.6832219375604082), MML::Real(1.3153839661540222), MML::Real(1.1263721073714085), MML::Real(0.8777054601504407), MML::Real(0.28994614552401766), MML::Real(0.06959539328104801)};
		return v;
	}
	inline constexpr Real symm_mat_10x10_cond_2 = MML::Real(54.47879912099627);
	inline constexpr int  symm_mat_10x10_rank = 10;
}

#endif // MML_LINEAR_ALG_EQ_SYSTEMS_SYM_DEFS_H

