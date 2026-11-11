#if !defined MML_LINEAR_ALG_EQ_SYSTEMS_COMPLEX_DEFS_H
#define MML_LINEAR_ALG_EQ_SYSTEMS_COMPLEX_DEFS_H

// SIOF-SAFE PATTERN: All test data uses inline functions with static local variables
// to avoid Static Initialization Order Fiasco issues across translation units.

#ifdef MML_SINGLE_HEADER_INCLUDE
#include MML_SINGLE_HEADER_INCLUDE
#elif defined(MML_USE_SINGLE_HEADER)
#include <MML.h>
#else
#include <mml/base/Matrix/Matrix.h>
#endif

namespace MML::TestBeds
{
    /***********************************************************************************************/
    /**************               Structured 2 x 2 complex matrices                    **************/
    /***********************************************************************************************/
    inline const Matrix<Complex>& mat_cmplx_symmetric_2x2() {
        static const Matrix<Complex> m{2, 2, {
            Complex(2.0, 1.0), Complex(1.0, -2.0),
            Complex(1.0, -2.0), Complex(3.0, 4.0)
        }};
        return m;
    }

    inline const Vector<Complex>& mat_cmplx_symmetric_2x2_rhs0() {
        static const Vector<Complex> v{Complex(5.0, -5.0), Complex(5.0, 5.0)};
        return v;
    }
    inline const Vector<Complex>& mat_cmplx_symmetric_2x2_rhs0_sol() {
        static const Vector<Complex> v{Complex(1.0, -1.0), Complex(2.0, 0.0)};
        return v;
    }

    inline Real mat_cmplx_symmetric_2x2_det_abs() { return std::sqrt(static_cast<MML::Real>(250.0)); }
    inline const Vector<Real>& mat_cmplx_symmetric_2x2_singular_values() {
        static const Vector<Real> v{MML::Real(5.678683713116302), MML::Real(2.7843403749343954)};
        return v;
    }
    inline Real mat_cmplx_symmetric_2x2_cond_2() { return MML::Real(2.039507733308834); }
    inline int mat_cmplx_symmetric_2x2_rank() { return 2; }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    inline const Matrix<Complex>& mat_cmplx_hermitian_2x2() {
        static const Matrix<Complex> m{2, 2, {
            Complex(4.0, 0.0), Complex(1.0, 2.0),
            Complex(1.0, -2.0), Complex(3.0, 0.0)
        }};
        return m;
    }

    inline const Vector<Complex>& mat_cmplx_hermitian_2x2_rhs0() {
        static const Vector<Complex> v{Complex(6.0, 0.0), Complex(5.0, -3.0)};
        return v;
    }
    inline const Vector<Complex>& mat_cmplx_hermitian_2x2_rhs0_sol() {
        static const Vector<Complex> v{Complex(1.0, -1.0), Complex(2.0, 0.0)};
        return v;
    }

    inline Real mat_cmplx_hermitian_2x2_det_abs() { return 7.0; }
    inline const Vector<Real>& mat_cmplx_hermitian_2x2_singular_values() {
        static const Vector<Real> v{MML::Real(5.79128784747792), MML::Real(1.2087121525220803)};
        return v;
    }
    inline Real mat_cmplx_hermitian_2x2_cond_2() { return MML::Real(4.7912878474779195); }
    inline int mat_cmplx_hermitian_2x2_rank() { return 2; }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    inline const Matrix<Complex>& mat_cmplx_unitary_2x2() {
        static const Real invSqrt2 = static_cast<MML::Real>(1.0) / std::sqrt(static_cast<MML::Real>(2.0));
        static const Matrix<Complex> m{2, 2, {
            Complex(invSqrt2, 0.0), Complex(0.0, invSqrt2),
            Complex(0.0, invSqrt2), Complex(invSqrt2, 0.0)
        }};
        return m;
    }

    inline const Vector<Complex>& mat_cmplx_unitary_2x2_rhs0() {
        static const Real invSqrt2 = static_cast<MML::Real>(1.0) / std::sqrt(static_cast<MML::Real>(2.0));
        static const Vector<Complex> v{Complex(invSqrt2, invSqrt2), Complex(MML::Real(3.0) * invSqrt2, invSqrt2)};
        return v;
    }
    inline const Vector<Complex>& mat_cmplx_unitary_2x2_rhs0_sol() {
        static const Vector<Complex> v{Complex(1.0, -1.0), Complex(2.0, 0.0)};
        return v;
    }

    inline Real mat_cmplx_unitary_2x2_det_abs() { return 1.0; }
    inline const Vector<Real>& mat_cmplx_unitary_2x2_singular_values() {
        static const Vector<Real> v{1.0, 1.0};
        return v;
    }
    inline Real mat_cmplx_unitary_2x2_cond_2() { return 1.0; }
    inline int mat_cmplx_unitary_2x2_rank() { return 2; }

    /***********************************************************************************************/
    /**************                     3 x 3 test matrices                           **************/
    /***********************************************************************************************/
    inline const Matrix<Complex>& mat_cmplx_3x3() {
        static const Matrix<Complex> m{3, 3, { 1.0, 2.0, -1.0, 
                                                    -1.0, 5.0, 6.0, 
                                                     3.0, 1.0, 1.0 }};
        return m;
    }

    inline const Vector<Complex>& mat_cmplx_3x3_rhs0() {
        static const Vector<Complex> v{1.0, 2.0, 1.0};
        return v;
    }
    inline const Vector<Complex>& mat_cmplx_3x3_rhs0_sol() {
        static const Vector<Complex> v{MML::Real(0.18867924528301885), MML::Real(0.41509433962264153), MML::Real(0.018867924528301921)};
        return v;
    }

    // Matrix properties for mat_cmplx_3x3 (same as real mat_3x3)
    inline Real mat_cmplx_3x3_det_abs() { return 53.0; }  // |det(A)|
    inline const Vector<Real>& mat_cmplx_3x3_singular_values() {
        static const Vector<Real> v{MML::Real(7.965662078134766), MML::Real(3.4347838553531793), MML::Real(1.9371568814853632)};
        return v;
    }
    inline Real mat_cmplx_3x3_cond_2() { return MML::Real(4.112039181691377); }
    inline int mat_cmplx_3x3_rank() { return 3; }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    inline const Matrix<Complex>& mat_cmplx_1_3x3() {
        static const Matrix<Complex> m{3, 3, { Complex(MML::Real(1.0), MML::Real(2.0)), Complex(MML::Real(2.0), -MML::Real(0.5)), Complex(-MML::Real(1.0), MML::Real(2.3)),
                                                         Complex(-2.0,1.0), Complex(4.0,-2.0),  Complex(6.0,5.0), 
                                                         Complex(2.5,1.0),  Complex(1.0,-1.0),  Complex(-2.0, -4.0) }};
        return m;
    }

    inline const Vector<Complex>& mat_cmplx_1_3x3_rhs0() {
        static const Vector<Complex> v{ Complex(1.5,0.5), Complex(-2.0,1.5), Complex(-1.0,-2.0)};
        return v;
    }
    inline const Vector<Complex>& mat_cmplx_1_3x3_rhs0_sol() {
        static const Vector<Complex> v{ Complex(MML::Real(0.693423783306577),-MML::Real(1.25000647849031)),
                                              Complex(-MML::Real(1.09549140996106),-MML::Real(0.596939956473329)),
                                              Complex(MML::Real(0.242505265145942),-MML::Real(0.451531010007621))};
        return v;
    }

    // Matrix properties for mat_cmplx_1_3x3
    inline Real mat_cmplx_1_3x3_det_abs() { return MML::Real(144.36547296686918); }  // |det(A)|
    inline const Vector<Real>& mat_cmplx_1_3x3_singular_values() {
        static const Vector<Real> v{MML::Real(10.20912316755949), MML::Real(5.028959556700652), MML::Real(2.815033298093671)};
        return v;
    }
    inline Real mat_cmplx_1_3x3_cond_2() { return MML::Real(3.6263579631689193); }
    inline int mat_cmplx_1_3x3_rank() { return 3; }

    /////////////////////////////////////////////////////////////////////////////////////////////////
    inline const Matrix<Complex>& mat_cmplx_1_5x5() {
        static const Matrix<Complex> m{5, 5, { Complex(MML::Real(1.0),MML::Real(2.0)),  Complex(MML::Real(2.0),-MML::Real(0.5)), Complex(-MML::Real(1.0), MML::Real(2.3)), Complex(MML::Real(1.8),-MML::Real(3.0)),  Complex(-MML::Real(2.5),MML::Real(3.7)),
                                                         Complex(-2.0,1.0), Complex(4.0,-2.0), Complex(6.0,5.0),   Complex(-1.0,-2.0), Complex(-1.0,-2.0), 
                                                         Complex(MML::Real(3.0),-MML::Real(3.0)), Complex(MML::Real(1.0),MML::Real(1.4)),  Complex(MML::Real(2.0),-MML::Real(1.0)),  Complex(-MML::Real(1.7),-MML::Real(2.6)), Complex(MML::Real(7.0),-MML::Real(3.0)),
                                                         Complex(-4.0,1.5), Complex(2.0,-3.0), Complex(-2.0,3.0),  Complex(1.0,2.0),   Complex(-5.0,6.0), 
                                                         Complex(2.5,1.0),  Complex(1.0,-1.0), Complex(-2.0,-4.0), Complex(-1.0,2.0),  Complex(2.0,3.0) }};
        return m;
    }

    inline const Vector<Complex>& mat_cmplx_1_5x5_rhs0() {
        static const Vector<Complex> v{ Complex(MML::Real(1.5),MML::Real(0.5)), Complex(-MML::Real(2.0),MML::Real(1.5)), Complex(-MML::Real(1.0),-MML::Real(2.0)), Complex(MML::Real(1.0),-MML::Real(2.5)), Complex(-MML::Real(1.8),MML::Real(2.1))};
        return v;
    }
    inline const Vector<Complex>& mat_cmplx_1_5x5_rhs0_sol() {
        static const Vector<Complex> v{ Complex(MML::Real(0.901430270043644),MML::Real(1.58987437993965)),
                                              Complex(-MML::Real(0.351975602660937),MML::Real(0.144870096280711)),
                                              Complex(MML::Real(0.558130434051583),-MML::Real(0.123040787906678)),
                                              Complex(MML::Real(0.894179352753342),-MML::Real(0.202080868658524)),
                                              Complex(-MML::Real(0.780553254235851),-MML::Real(0.467650784949547))};
        return v;
    }

    // Matrix properties for mat_cmplx_1_5x5
    inline Real mat_cmplx_1_5x5_det_abs() { return MML::Real(10168.952903809948); }  // |det(A)|
    inline const Vector<Real>& mat_cmplx_1_5x5_singular_values() {
        static const Vector<Real> v{MML::Real(15.64174282618698), MML::Real(9.00912568656689), MML::Real(5.696377816847538), MML::Real(4.137858398839877), MML::Real(2.380802820168614)};
        return v;
    }
    inline Real mat_cmplx_1_5x5_cond_2() { return MML::Real(6.570844779227217); }
    inline int mat_cmplx_1_5x5_rank() { return 5; }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    inline const Matrix<Complex>& mat_cmplx_2_5x5() {
        static const Matrix<Complex> m{5, 5, {
            Complex(-MML::Real(4.7),-MML::Real(1.0)), Complex(-MML::Real(0.1),-MML::Real(3.4)), Complex(MML::Real(3.5),MML::Real(0.9)), Complex(-MML::Real(0.9),-MML::Real(1.2)), Complex(MML::Real(1.3),-MML::Real(2.1)),
            Complex(-MML::Real(2.1),-MML::Real(4.8)), Complex(MML::Real(0.4),-MML::Real(0.3)), Complex(-MML::Real(2.7),MML::Real(4.1)), Complex(-MML::Real(2.6),-MML::Real(4.0)), Complex(MML::Real(4.5),-MML::Real(0.7)),
            Complex(MML::Real(1.5),-MML::Real(3.7)), Complex(-MML::Real(4.8),-MML::Real(3.9)), Complex(MML::Real(5.0),MML::Real(4.3)), Complex(MML::Real(2.4),MML::Real(1.2)), Complex(MML::Real(2.2),-MML::Real(2.3)),
            Complex(-MML::Real(0.9),MML::Real(2.3)), Complex(-MML::Real(3.0),-MML::Real(3.1)), Complex(MML::Real(4.6),-MML::Real(2.2)), Complex(-MML::Real(2.8),-MML::Real(2.1)), Complex(-MML::Real(4.5),MML::Real(0.2)),
            Complex(MML::Real(3.8),-MML::Real(1.1)), Complex(MML::Real(3.9),-MML::Real(2.7)), Complex(MML::Real(4.6),MML::Real(3.1)), Complex(MML::Real(1.1),MML::Real(2.1)), Complex(MML::Real(4.9),MML::Real(3.5))
        }};
        return m;
    }

    inline const Vector<Complex>& mat_cmplx_2_5x5_rhs0() {
        static const Vector<Complex> v{
            Complex(MML::Real(1.1),MML::Real(2.7)),
            Complex(MML::Real(0.1),MML::Real(2.9)),
            Complex(MML::Real(0.7),MML::Real(2.1)),
            Complex(-MML::Real(0.1),MML::Real(0.5)),
            Complex(-MML::Real(0.2),-MML::Real(1.8))
        };
        return v;
    }
    inline const Vector<Complex>& mat_cmplx_2_5x5_rhs0_sol() {
        static const Vector<Complex> v{
            Complex(-MML::Real(0.364324249141192),-MML::Real(0.683814875129271)),
            Complex(-MML::Real(0.422202239716293),-MML::Real(0.369635385572325)),
            Complex(MML::Real(0.189016352402793),-MML::Real(0.296292202554985)),
            Complex(MML::Real(0.145166170094375),MML::Real(0.115685130841116)),
            Complex(MML::Real(0.489514812714704),-MML::Real(0.136538611915636))
        };
        return v;
    }

    // Matrix properties for mat_cmplx_2_5x5
    inline Real mat_cmplx_2_5x5_det_abs() { return MML::Real(13745.1507590657); }
    inline const Vector<Real>& mat_cmplx_2_5x5_singular_values() {
        static const Vector<Real> v{MML::Real(16.1973119257448), MML::Real(10.7879205362089), MML::Real(6.40784516509668), MML::Real(6.0423260356214), MML::Real(2.03166761650746)};
        return v;
    }
    inline Real mat_cmplx_2_5x5_cond_2() { return MML::Real(7.97242215908764); }
    inline int mat_cmplx_2_5x5_rank() { return 5; }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    inline const Matrix<Complex>& mat_cmplx_3_5x5() {
        static const Matrix<Complex> m{5, 5, {
            Complex(-MML::Real(3.9),-MML::Real(1.8)), Complex(MML::Real(1.8),-MML::Real(4.3)), Complex(MML::Real(0.3),-MML::Real(4.5)), Complex(-MML::Real(1.3),-MML::Real(0.6)), Complex(-MML::Real(0.9),MML::Real(2.0)),
            Complex(MML::Real(0.9),MML::Real(2.9)), Complex(MML::Real(2.2),-MML::Real(2.6)), Complex(-MML::Real(3.3),MML::Real(1.5)), Complex(-MML::Real(4.3),MML::Real(3.7)), Complex(MML::Real(0.0),-MML::Real(2.9)),
            Complex(-MML::Real(0.7),-MML::Real(3.7)), Complex(MML::Real(4.5),MML::Real(3.1)), Complex(-MML::Real(0.8),-MML::Real(4.6)), Complex(-MML::Real(2.5),-MML::Real(1.4)), Complex(-MML::Real(1.0),MML::Real(1.8)),
            Complex(MML::Real(1.3),MML::Real(3.3)), Complex(MML::Real(0.8),-MML::Real(0.6)), Complex(MML::Real(2.8),-MML::Real(4.2)), Complex(MML::Real(1.2),-MML::Real(4.1)), Complex(MML::Real(1.3),MML::Real(1.9)),
            Complex(-MML::Real(1.3),MML::Real(2.9)), Complex(-MML::Real(3.3),-MML::Real(1.8)), Complex(MML::Real(2.5),-MML::Real(1.2)), Complex(-MML::Real(1.0),MML::Real(2.3)), Complex(MML::Real(1.6),-MML::Real(2.2))
        }};
        return m;
    }

    inline const Vector<Complex>& mat_cmplx_3_5x5_rhs0() {
        static const Vector<Complex> v{
            Complex(-MML::Real(2.2),MML::Real(2.5)),
            Complex(-MML::Real(1.1),-MML::Real(0.7)),
            Complex(MML::Real(0.3),-MML::Real(2.2)),
            Complex(-MML::Real(2.2),MML::Real(3.0)),
            Complex(-MML::Real(0.0),MML::Real(1.6))
        };
        return v;
    }
    inline const Vector<Complex>& mat_cmplx_3_5x5_rhs0_sol() {
        static const Vector<Complex> v{
            Complex(MML::Real(0.061637980467496),MML::Real(0.231305484458262)),
            Complex(-MML::Real(0.608924639903936),-MML::Real(0.383281231895997)),
            Complex(-MML::Real(0.283005443992466),-MML::Real(0.007473368712111)),
            Complex(-MML::Real(0.291587621773846),MML::Real(0.038722221408967)),
            Complex(MML::Real(0.062649872788515),-MML::Real(0.070415668769379))
        };
        return v;
    }

    // Matrix properties for mat_cmplx_3_5x5
    inline Real mat_cmplx_3_5x5_det_abs() { return MML::Real(6319.66771674941); }
    inline const Vector<Real>& mat_cmplx_3_5x5_singular_values() {
        static const Vector<Real> v{MML::Real(13.3205460279412), MML::Real(9.39155667384619), MML::Real(6.02046598723545), MML::Real(4.85473027591569), MML::Real(1.72838074627172)};
        return v;
    }
    inline Real mat_cmplx_3_5x5_cond_2() { return MML::Real(7.70695117767013); }
    inline int mat_cmplx_3_5x5_rank() { return 5; }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    inline const Matrix<Complex>& mat_cmplx_4_5x5() {
        static const Matrix<Complex> m{5, 5, {
            Complex(-MML::Real(2.2),-MML::Real(4.9)), Complex(MML::Real(1.4),-MML::Real(0.6)), Complex(MML::Real(1.2),-MML::Real(4.5)), Complex(MML::Real(2.8),MML::Real(1.7)), Complex(-MML::Real(0.2),-MML::Real(4.7)),
            Complex(MML::Real(4.6),-MML::Real(1.2)), Complex(MML::Real(5.0),MML::Real(3.2)), Complex(MML::Real(3.3),MML::Real(4.0)), Complex(-MML::Real(3.0),-MML::Real(1.7)), Complex(-MML::Real(3.7),-MML::Real(4.0)),
            Complex(MML::Real(1.2),-MML::Real(3.1)), Complex(MML::Real(0.9),-MML::Real(3.4)), Complex(MML::Real(1.4),MML::Real(4.4)), Complex(MML::Real(4.9),MML::Real(4.9)), Complex(-MML::Real(1.6),-MML::Real(4.1)),
            Complex(MML::Real(3.1),-MML::Real(2.4)), Complex(MML::Real(1.7),-MML::Real(3.7)), Complex(MML::Real(0.4),MML::Real(2.2)), Complex(-MML::Real(4.3),-MML::Real(3.1)), Complex(-MML::Real(4.1),-MML::Real(4.8)),
            Complex(MML::Real(3.7),-MML::Real(3.5)), Complex(MML::Real(0.9),MML::Real(4.0)), Complex(-MML::Real(4.2),MML::Real(1.1)), Complex(-MML::Real(1.6),-MML::Real(4.0)), Complex(-MML::Real(1.3),-MML::Real(3.3))
        }};
        return m;
    }

    inline const Vector<Complex>& mat_cmplx_4_5x5_rhs0() {
        static const Vector<Complex> v{
            Complex(MML::Real(1.5),-MML::Real(1.1)),
            Complex(-MML::Real(2.3),-MML::Real(2.4)),
            Complex(-3.0,0.5),
            Complex(-MML::Real(1.2),-MML::Real(1.8)),
            Complex(-3.0,2.5)
        };
        return v;
    }
    inline const Vector<Complex>& mat_cmplx_4_5x5_rhs0_sol() {
        static const Vector<Complex> v{
            Complex(-MML::Real(0.926407716133575),MML::Real(0.033047493912236)),
            Complex(-MML::Real(0.055819119213253),MML::Real(0.016480814215369)),
            Complex(MML::Real(0.411015404552370),-MML::Real(0.188896066902875)),
            Complex(-MML::Real(0.400191856408234),MML::Real(0.311772178213589)),
            Complex(MML::Real(0.781799961064350),MML::Real(0.325947296047793))
        };
        return v;
    }

    // Matrix properties for mat_cmplx_4_5x5
    inline Real mat_cmplx_4_5x5_det_abs() { return MML::Real(35140.2960834746); }
    inline const Vector<Real>& mat_cmplx_4_5x5_singular_values() {
        static const Vector<Real> v{MML::Real(16.4968668117334), MML::Real(11.4058088404484), MML::Real(8.20284126627676), MML::Real(6.24295130117635), MML::Real(3.64689790129684)};
        return v;
    }
    inline Real mat_cmplx_4_5x5_cond_2() { return MML::Real(4.52353404406171); }
    inline int mat_cmplx_4_5x5_rank() { return 5; }

    /***********************************************************************************************/
    /**************                     8 x 8 test matrices                           **************/
    /***********************************************************************************************/
    inline const Matrix<Complex>& mat_cmplx_1_8x8() {
        static const Matrix<Complex> m{8, 8, {
            Complex(MML::Real(4.5),MML::Real(4.8)), Complex(MML::Real(2.4),-MML::Real(0.5)), Complex(-MML::Real(4.6),MML::Real(4.0)), Complex(-MML::Real(1.7),MML::Real(1.1)), Complex(-MML::Real(1.8),MML::Real(1.4)), Complex(-MML::Real(0.9),MML::Real(3.7)), Complex(-MML::Real(4.3),MML::Real(4.4)), Complex(-MML::Real(1.7),-MML::Real(1.2)),
            Complex(MML::Real(2.8),MML::Real(3.0)), Complex(-MML::Real(3.3),MML::Real(4.6)), Complex(MML::Real(0.3),-MML::Real(1.3)), Complex(-MML::Real(2.5),-MML::Real(0.8)), Complex(-MML::Real(4.4),-MML::Real(1.1)), Complex(-MML::Real(3.1),MML::Real(0.4)), Complex(-MML::Real(0.6),-MML::Real(2.7)), Complex(MML::Real(1.5),MML::Real(2.7)),
            Complex(-MML::Real(1.0),-MML::Real(1.8)), Complex(MML::Real(4.1),MML::Real(3.2)), Complex(MML::Real(1.6),-MML::Real(3.9)), Complex(MML::Real(2.0),MML::Real(0.8)), Complex(MML::Real(0.4),-MML::Real(0.5)), Complex(-MML::Real(0.8),MML::Real(0.5)), Complex(-MML::Real(4.0),-MML::Real(0.1)), Complex(-MML::Real(1.2),MML::Real(2.0)),
            Complex(-MML::Real(1.9),-MML::Real(0.9)), Complex(MML::Real(4.9),-MML::Real(2.0)), Complex(MML::Real(4.2),-MML::Real(1.2)), Complex(-MML::Real(2.0),-MML::Real(1.6)), Complex(-MML::Real(3.0),MML::Real(4.1)), Complex(-MML::Real(2.2),MML::Real(4.1)), Complex(-MML::Real(4.4),MML::Real(0.9)), Complex(MML::Real(1.7),-MML::Real(2.5)),
            Complex(-MML::Real(2.1),MML::Real(3.3)), Complex(MML::Real(3.2),-MML::Real(3.6)), Complex(MML::Real(3.3),MML::Real(3.3)), Complex(-MML::Real(1.5),-MML::Real(3.9)), Complex(MML::Real(0.3),MML::Real(4.8)), Complex(-MML::Real(0.3),MML::Real(2.2)), Complex(-MML::Real(2.1),-MML::Real(4.2)), Complex(MML::Real(0.7),-MML::Real(1.6)),
            Complex(-MML::Real(1.4),MML::Real(4.2)), Complex(-MML::Real(1.4),MML::Real(1.0)), Complex(-MML::Real(4.2),MML::Real(2.6)), Complex(-MML::Real(3.8),MML::Real(3.9)), Complex(MML::Real(2.2),-MML::Real(3.4)), Complex(-MML::Real(4.5),MML::Real(2.7)), Complex(-MML::Real(2.6),MML::Real(0.3)), Complex(-MML::Real(3.5),MML::Real(1.3)),
            Complex(-MML::Real(2.8),MML::Real(3.0)), Complex(MML::Real(2.7),MML::Real(3.2)), Complex(-MML::Real(3.7),MML::Real(4.0)), Complex(MML::Real(1.6),-MML::Real(0.7)), Complex(-MML::Real(4.8),MML::Real(0.8)), Complex(-MML::Real(3.4),-MML::Real(1.4)), Complex(MML::Real(0.3),-MML::Real(3.0)), Complex(-MML::Real(1.7),MML::Real(4.2)),
            Complex(MML::Real(1.8),MML::Real(2.8)), Complex(MML::Real(4.9),-MML::Real(4.7)), Complex(MML::Real(3.0),-MML::Real(4.8)), Complex(-MML::Real(4.8),-MML::Real(1.6)), Complex(MML::Real(4.8),MML::Real(4.3)), Complex(-MML::Real(1.1),MML::Real(0.4)), Complex(-MML::Real(0.1),-MML::Real(3.8)), Complex(-MML::Real(2.9),MML::Real(1.5))
        }};
        return m;
    }

    inline const Vector<Complex>& mat_cmplx_1_8x8_rhs0() {
        static const Vector<Complex> v{
            Complex(-MML::Real(2.4),-MML::Real(0.6)),
            Complex(-MML::Real(1.2),MML::Real(1.5)),
            Complex(MML::Real(1.7),-MML::Real(2.0)),
            Complex(MML::Real(0.0),MML::Real(2.7)),
            Complex(-MML::Real(2.3),-MML::Real(2.7)),
            Complex(-MML::Real(2.4),MML::Real(2.2)),
            Complex(MML::Real(2.2),MML::Real(0.4)),
            Complex(MML::Real(2.2),MML::Real(2.3))
        };
        return v;
    }
    inline const Vector<Complex>& mat_cmplx_1_8x8_rhs0_sol() {
        static const Vector<Complex> v{
            Complex(MML::Real(0.279983758960108),-MML::Real(0.280318189874449)),
            Complex(-MML::Real(0.185720534820703),MML::Real(0.661437184412004)),
            Complex(MML::Real(0.080737838282574),MML::Real(0.657637228492826)),
            Complex(MML::Real(1.476088586253125),-MML::Real(1.410794553549068)),
            Complex(-MML::Real(0.290319020357946),-MML::Real(0.864416476092015)),
            Complex(-MML::Real(1.032934559655212),MML::Real(0.267746969846193)),
            Complex(MML::Real(0.293628372974839),MML::Real(0.398814739983594)),
            Complex(-MML::Real(0.066786432249243),MML::Real(0.250909440756107))
        };
        return v;
    }

    // Matrix properties for mat_cmplx_1_8x8
    inline Real mat_cmplx_1_8x8_det_abs() { return MML::Real(8254361.91150796); }
    inline const Vector<Real>& mat_cmplx_1_8x8_singular_values() {
        static const Vector<Real> v{MML::Real(20.7329056394107), MML::Real(17.1020495055299), MML::Real(12.3211712343142), MML::Real(9.14017764485158), MML::Real(7.61605888733823), MML::Real(5.18445288908348), MML::Real(3.62879622218723), MML::Real(1.44268906430648)};
        return v;
    }
    inline Real mat_cmplx_1_8x8_cond_2() { return MML::Real(14.3710146228753); }
    inline int mat_cmplx_1_8x8_rank() { return 8; }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    inline const Matrix<Complex>& mat_cmplx_2_8x8() {
        static const Matrix<Complex> m{8, 8, {
            Complex(-MML::Real(1.8),MML::Real(0.4)), Complex(MML::Real(2.6),MML::Real(3.1)), Complex(MML::Real(3.8),MML::Real(0.8)), Complex(MML::Real(1.9),-MML::Real(0.6)), Complex(MML::Real(2.8),MML::Real(1.7)), Complex(MML::Real(4.1),-MML::Real(1.1)), Complex(MML::Real(2.3),MML::Real(4.8)), Complex(-MML::Real(2.0),MML::Real(1.1)),
            Complex(-MML::Real(2.2),-MML::Real(1.7)), Complex(MML::Real(2.1),-MML::Real(1.8)), Complex(MML::Real(0.4),-MML::Real(0.3)), Complex(MML::Real(5.0),MML::Real(4.7)), Complex(MML::Real(4.7),MML::Real(0.1)), Complex(MML::Real(3.5),MML::Real(3.2)), Complex(-MML::Real(1.9),MML::Real(0.1)), Complex(MML::Real(1.8),-MML::Real(1.4)),
            Complex(-MML::Real(4.4),MML::Real(3.7)), Complex(MML::Real(0.2),-MML::Real(2.9)), Complex(-MML::Real(3.1),MML::Real(2.0)), Complex(MML::Real(4.3),MML::Real(0.7)), Complex(MML::Real(5.0),MML::Real(4.3)), Complex(-MML::Real(4.0),MML::Real(4.1)), Complex(-MML::Real(1.2),MML::Real(4.7)), Complex(MML::Real(3.7),-MML::Real(3.4)),
            Complex(MML::Real(3.9),MML::Real(1.5)), Complex(-MML::Real(3.7),MML::Real(2.9)), Complex(MML::Real(0.9),MML::Real(3.6)), Complex(MML::Real(0.4),MML::Real(4.2)), Complex(MML::Real(4.2),-MML::Real(0.1)), Complex(MML::Real(0.1),-MML::Real(4.9)), Complex(-MML::Real(2.6),-MML::Real(1.6)), Complex(MML::Real(0.7),MML::Real(1.6)),
            Complex(MML::Real(1.8),MML::Real(4.2)), Complex(-MML::Real(0.5),-MML::Real(1.6)), Complex(-MML::Real(2.6),-MML::Real(4.9)), Complex(MML::Real(1.0),-MML::Real(3.0)), Complex(-MML::Real(0.0),-MML::Real(2.4)), Complex(-MML::Real(0.6),-MML::Real(3.0)), Complex(-MML::Real(1.3),-MML::Real(0.4)), Complex(-MML::Real(1.4),MML::Real(0.4)),
            Complex(-MML::Real(0.8),MML::Real(4.0)), Complex(MML::Real(0.6),MML::Real(4.9)), Complex(MML::Real(3.8),MML::Real(2.3)), Complex(MML::Real(1.9),-MML::Real(4.1)), Complex(-MML::Real(1.9),MML::Real(3.3)), Complex(-MML::Real(3.6),MML::Real(1.0)), Complex(-MML::Real(2.8),-MML::Real(4.0)), Complex(MML::Real(3.0),MML::Real(1.7)),
            Complex(MML::Real(1.4),MML::Real(3.0)), Complex(-MML::Real(4.7),-MML::Real(1.3)), Complex(-MML::Real(4.6),MML::Real(1.3)), Complex(-MML::Real(3.2),-MML::Real(4.2)), Complex(MML::Real(0.4),MML::Real(4.9)), Complex(MML::Real(1.9),MML::Real(1.8)), Complex(MML::Real(0.0),-MML::Real(1.2)), Complex(-MML::Real(1.8),-MML::Real(1.9)),
            Complex(-MML::Real(1.9),MML::Real(3.7)), Complex(-MML::Real(1.2),-MML::Real(0.2)), Complex(-MML::Real(1.5),MML::Real(3.1)), Complex(-MML::Real(3.4),-MML::Real(2.5)), Complex(-MML::Real(1.0),MML::Real(2.9)), Complex(MML::Real(1.3),MML::Real(1.3)), Complex(-MML::Real(1.7),-MML::Real(4.3)), Complex(MML::Real(4.5),-MML::Real(3.7))
        }};
        return m;
    }

    inline const Vector<Complex>& mat_cmplx_2_8x8_rhs0() {
        static const Vector<Complex> v{
            Complex(-MML::Real(0.8),MML::Real(0.5)),
            Complex(-MML::Real(0.8),-MML::Real(0.5)),
            Complex(MML::Real(0.8),MML::Real(1.8)),
            Complex(-MML::Real(1.6),-MML::Real(0.6)),
            Complex(-1.0,1.0),
            Complex(MML::Real(1.1),MML::Real(2.7)),
            Complex(-MML::Real(1.2),MML::Real(0.9)),
            Complex(MML::Real(2.6),MML::Real(1.9))
        };
        return v;
    }
    inline const Vector<Complex>& mat_cmplx_2_8x8_rhs0_sol() {
        static const Vector<Complex> v{
            Complex(-MML::Real(0.071230196831102),-MML::Real(0.301989480770654)),
            Complex(MML::Real(0.579940424169923),MML::Real(0.349107697938878)),
            Complex(MML::Real(0.082682783590385),-MML::Real(0.579849142434740)),
            Complex(-MML::Real(0.357408952376676),MML::Real(0.047989986784227)),
            Complex(-MML::Real(0.090219498196085),MML::Real(0.319245310867363)),
            Complex(-MML::Real(0.095237471765996),-MML::Real(0.142075932285791)),
            Complex(-MML::Real(0.236103821389095),-MML::Real(0.005801913163186)),
            Complex(-MML::Real(0.098903720964927),-MML::Real(0.051871533889682))
        };
        return v;
    }

    // Matrix properties for mat_cmplx_2_8x8
    inline Real mat_cmplx_2_8x8_det_abs() { return MML::Real(9681017.30839477); }
    inline const Vector<Real>& mat_cmplx_2_8x8_singular_values() {
        static const Vector<Real> v{MML::Real(17.9783905280538), MML::Real(17.0468894260995), MML::Real(13.0089171954949), MML::Real(10.7846443438304), MML::Real(8.43566113581774), MML::Real(5.55392678710442), MML::Real(2.69870097150435), MML::Real(1.78075439196711)};
        return v;
    }
    inline Real mat_cmplx_2_8x8_cond_2() { return MML::Real(10.0959405795394); }
    inline int mat_cmplx_2_8x8_rank() { return 8; }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    inline const Matrix<Complex>& mat_cmplx_3_8x8() {
        static const Matrix<Complex> m{8, 8, {
            Complex(-MML::Real(3.5),MML::Real(2.4)), Complex(-MML::Real(3.8),-MML::Real(3.1)), Complex(-MML::Real(1.8),-MML::Real(4.5)), Complex(MML::Real(1.8),MML::Real(1.6)), Complex(MML::Real(1.1),-MML::Real(2.6)), Complex(MML::Real(3.8),MML::Real(3.4)), Complex(-MML::Real(1.9),-MML::Real(1.3)), Complex(-MML::Real(1.8),MML::Real(3.9)),
            Complex(-MML::Real(2.5),MML::Real(4.4)), Complex(-MML::Real(0.4),-MML::Real(0.9)), Complex(MML::Real(1.9),-MML::Real(0.8)), Complex(MML::Real(2.2),-MML::Real(2.2)), Complex(-MML::Real(3.5),-MML::Real(1.3)), Complex(MML::Real(1.8),-MML::Real(3.1)), Complex(-MML::Real(3.5),-MML::Real(2.2)), Complex(MML::Real(2.8),MML::Real(4.0)),
            Complex(-MML::Real(4.2),-MML::Real(3.2)), Complex(-MML::Real(4.1),-MML::Real(2.8)), Complex(-MML::Real(3.4),-MML::Real(2.2)), Complex(-MML::Real(4.3),-MML::Real(4.8)), Complex(MML::Real(0.2),-MML::Real(4.1)), Complex(MML::Real(4.2),-MML::Real(0.2)), Complex(MML::Real(0.6),-MML::Real(2.0)), Complex(MML::Real(0.0),MML::Real(0.0)),
            Complex(-MML::Real(2.3),MML::Real(3.9)), Complex(MML::Real(4.0),MML::Real(1.1)), Complex(-MML::Real(4.1),-MML::Real(4.5)), Complex(-MML::Real(0.2),MML::Real(1.5)), Complex(-MML::Real(0.3),-MML::Real(3.3)), Complex(MML::Real(2.1),-MML::Real(0.3)), Complex(MML::Real(4.1),MML::Real(4.1)), Complex(-MML::Real(1.0),MML::Real(1.1)),
            Complex(MML::Real(2.5),MML::Real(2.9)), Complex(-MML::Real(1.3),-MML::Real(3.1)), Complex(MML::Real(0.1),-MML::Real(3.7)), Complex(MML::Real(2.0),MML::Real(3.0)), Complex(-MML::Real(2.1),MML::Real(0.9)), Complex(MML::Real(2.9),MML::Real(0.0)), Complex(-MML::Real(3.5),-MML::Real(1.8)), Complex(MML::Real(1.8),MML::Real(4.4)),
            Complex(MML::Real(0.6),-MML::Real(3.0)), Complex(MML::Real(4.7),-MML::Real(1.2)), Complex(MML::Real(2.5),MML::Real(4.0)), Complex(-MML::Real(4.4),-MML::Real(0.9)), Complex(-MML::Real(1.9),MML::Real(4.8)), Complex(MML::Real(0.7),MML::Real(0.2)), Complex(-MML::Real(3.6),MML::Real(1.1)), Complex(-MML::Real(0.7),MML::Real(0.5)),
            Complex(MML::Real(2.0),MML::Real(1.8)), Complex(MML::Real(2.2),MML::Real(0.0)), Complex(MML::Real(4.4),MML::Real(1.7)), Complex(MML::Real(4.6),MML::Real(1.9)), Complex(-MML::Real(4.7),-MML::Real(1.4)), Complex(-MML::Real(4.1),MML::Real(1.9)), Complex(-MML::Real(1.5),MML::Real(2.6)), Complex(-MML::Real(1.6),MML::Real(1.5)),
            Complex(MML::Real(0.4),MML::Real(3.1)), Complex(-MML::Real(3.7),-MML::Real(0.8)), Complex(-MML::Real(1.9),-MML::Real(3.4)), Complex(MML::Real(0.6),MML::Real(0.3)), Complex(MML::Real(3.7),MML::Real(1.0)), Complex(MML::Real(3.6),MML::Real(0.0)), Complex(-MML::Real(0.6),MML::Real(0.0)), Complex(MML::Real(1.3),MML::Real(4.2))
        }};
        return m;
    }

    inline const Vector<Complex>& mat_cmplx_3_8x8_rhs0() {
        static const Vector<Complex> v{
            Complex(MML::Real(0.8),MML::Real(1.1)),
            Complex(MML::Real(2.3),MML::Real(2.7)),
            Complex(-MML::Real(1.4),MML::Real(0.0)),
            Complex(-MML::Real(1.1),MML::Real(0.8)),
            Complex(MML::Real(2.6),MML::Real(2.2)),
            Complex(MML::Real(1.5),MML::Real(2.1)),
            Complex(-MML::Real(3.0),MML::Real(1.4)),
            Complex(-MML::Real(2.9),MML::Real(0.9))
        };
        return v;
    }
    inline const Vector<Complex>& mat_cmplx_3_8x8_rhs0_sol() {
        static const Vector<Complex> v{
            Complex(-MML::Real(0.639415913735030),-MML::Real(0.641550528548533)),
            Complex(MML::Real(0.468937846956392),MML::Real(0.485780618653904)),
            Complex(MML::Real(0.162250528053269),-MML::Real(0.458584525975030)),
            Complex(MML::Real(0.403585883119101),-MML::Real(0.615326804902230)),
            Complex(-MML::Real(0.225178977849940),-MML::Real(0.170686259188212)),
            Complex(MML::Real(1.311877929063035),-MML::Real(1.146147631413710)),
            Complex(-MML::Real(0.351915453593622),MML::Real(0.227370481021649)),
            Complex(MML::Real(1.839401613755739),MML::Real(2.025842111839458))
        };
        return v;
    }

    // Matrix properties for mat_cmplx_3_8x8
    inline Real mat_cmplx_3_8x8_det_abs() { return MML::Real(5323054.71961043); }
    inline const Vector<Real>& mat_cmplx_3_8x8_singular_values() {
        static const Vector<Real> v{MML::Real(19.1813429563283), MML::Real(15.3590286176189), MML::Real(12.4259062472838), MML::Real(9.69670848034898), MML::Real(6.7671725337737), MML::Real(5.44102025627927), MML::Real(3.49898640472805), MML::Real(1.16395425345961)};
        return v;
    }
    inline Real mat_cmplx_3_8x8_cond_2() { return MML::Real(16.4794646347274); }
    inline int mat_cmplx_3_8x8_rank() { return 8; }

    /***********************************************************************************************/
    /**************                    10 x 10 test matrices                          **************/
    /***********************************************************************************************/
    inline const Matrix<Complex>& mat_cmplx_1_10x10() {
        static const Matrix<Complex> m{10, 10, {
            Complex(-MML::Real(4.9),MML::Real(0.3)), Complex(-MML::Real(0.7),MML::Real(3.7)), Complex(-MML::Real(2.5),MML::Real(0.5)), Complex(MML::Real(0.1),-MML::Real(0.7)), Complex(MML::Real(2.5),MML::Real(2.5)), Complex(-MML::Real(1.5),MML::Real(3.7)), Complex(-MML::Real(2.2),MML::Real(3.8)), Complex(MML::Real(4.0),-MML::Real(3.3)), Complex(MML::Real(2.1),-MML::Real(3.0)), Complex(MML::Real(4.0),MML::Real(1.0)),
            Complex(-MML::Real(0.5),MML::Real(1.7)), Complex(-MML::Real(3.7),-MML::Real(2.4)), Complex(-MML::Real(4.9),MML::Real(3.5)), Complex(MML::Real(4.0),MML::Real(5.0)), Complex(MML::Real(3.7),MML::Real(1.2)), Complex(-MML::Real(4.7),MML::Real(0.2)), Complex(-MML::Real(1.3),MML::Real(2.9)), Complex(-MML::Real(0.7),MML::Real(0.8)), Complex(-MML::Real(1.5),MML::Real(0.5)), Complex(MML::Real(1.0),MML::Real(1.2)),
            Complex(-MML::Real(2.6),MML::Real(3.1)), Complex(-MML::Real(3.4),MML::Real(3.7)), Complex(-MML::Real(0.9),-MML::Real(4.4)), Complex(-MML::Real(1.8),MML::Real(2.6)), Complex(-MML::Real(4.3),-MML::Real(2.0)), Complex(-MML::Real(4.8),-MML::Real(4.2)), Complex(MML::Real(1.0),MML::Real(1.2)), Complex(MML::Real(1.5),-MML::Real(4.4)), Complex(-MML::Real(0.2),-MML::Real(3.9)), Complex(MML::Real(1.5),MML::Real(0.3)),
            Complex(MML::Real(1.2),-MML::Real(4.9)), Complex(MML::Real(2.8),MML::Real(0.7)), Complex(-MML::Real(2.3),MML::Real(2.8)), Complex(MML::Real(3.3),MML::Real(1.8)), Complex(MML::Real(4.4),-MML::Real(0.8)), Complex(-MML::Real(2.3),-MML::Real(0.1)), Complex(MML::Real(2.0),-MML::Real(2.0)), Complex(MML::Real(1.1),MML::Real(4.4)), Complex(MML::Real(3.2),MML::Real(1.2)), Complex(-MML::Real(0.2),-MML::Real(4.8)),
            Complex(-MML::Real(2.8),-MML::Real(1.8)), Complex(MML::Real(3.9),MML::Real(4.1)), Complex(-MML::Real(0.7),-MML::Real(2.7)), Complex(-MML::Real(0.9),MML::Real(0.0)), Complex(MML::Real(3.5),-MML::Real(5.0)), Complex(-MML::Real(0.0),-MML::Real(3.6)), Complex(MML::Real(3.0),MML::Real(4.3)), Complex(-MML::Real(3.0),-MML::Real(4.4)), Complex(MML::Real(4.5),MML::Real(4.0)), Complex(-MML::Real(1.0),MML::Real(3.0)),
            Complex(MML::Real(3.2),-MML::Real(3.2)), Complex(MML::Real(3.0),-MML::Real(3.6)), Complex(-MML::Real(2.9),MML::Real(3.0)), Complex(-MML::Real(3.2),-MML::Real(3.9)), Complex(-MML::Real(1.8),MML::Real(2.9)), Complex(-MML::Real(4.5),MML::Real(3.2)), Complex(-MML::Real(0.6),MML::Real(2.7)), Complex(MML::Real(3.7),MML::Real(4.5)), Complex(-MML::Real(3.9),MML::Real(2.3)), Complex(MML::Real(0.3),MML::Real(1.2)),
            Complex(-MML::Real(4.5),MML::Real(1.6)), Complex(-MML::Real(4.9),-MML::Real(4.1)), Complex(MML::Real(4.9),MML::Real(1.3)), Complex(MML::Real(2.4),-MML::Real(2.9)), Complex(MML::Real(4.8),MML::Real(4.7)), Complex(MML::Real(4.8),MML::Real(4.4)), Complex(MML::Real(3.2),-MML::Real(1.7)), Complex(-MML::Real(0.5),-MML::Real(4.3)), Complex(-MML::Real(3.1),-MML::Real(4.8)), Complex(-MML::Real(2.4),-MML::Real(0.6)),
            Complex(MML::Real(2.1),MML::Real(3.4)), Complex(MML::Real(3.8),-MML::Real(4.0)), Complex(MML::Real(0.5),MML::Real(1.9)), Complex(-MML::Real(4.2),MML::Real(4.3)), Complex(-MML::Real(1.0),-MML::Real(4.2)), Complex(MML::Real(2.8),MML::Real(0.1)), Complex(-MML::Real(1.7),MML::Real(0.5)), Complex(-MML::Real(3.0),MML::Real(4.3)), Complex(MML::Real(4.4),MML::Real(0.0)), Complex(MML::Real(0.4),MML::Real(3.6)),
            Complex(MML::Real(4.4),MML::Real(0.0)), Complex(MML::Real(1.2),MML::Real(3.7)), Complex(-MML::Real(0.3),-MML::Real(2.1)), Complex(MML::Real(1.9),MML::Real(4.1)), Complex(MML::Real(1.2),MML::Real(0.2)), Complex(MML::Real(0.0),MML::Real(0.4)), Complex(-MML::Real(3.5),MML::Real(3.0)), Complex(MML::Real(1.2),MML::Real(1.2)), Complex(MML::Real(2.9),-MML::Real(2.3)), Complex(MML::Real(4.4),-MML::Real(4.7)),
            Complex(MML::Real(1.8),-MML::Real(3.0)), Complex(MML::Real(2.6),MML::Real(3.8)), Complex(MML::Real(2.9),MML::Real(3.5)), Complex(-MML::Real(4.9),MML::Real(0.0)), Complex(MML::Real(2.2),-MML::Real(3.7)), Complex(-MML::Real(1.6),-MML::Real(2.2)), Complex(MML::Real(3.9),MML::Real(0.6)), Complex(MML::Real(0.6),MML::Real(1.8)), Complex(MML::Real(4.4),MML::Real(2.0)), Complex(-MML::Real(3.2),MML::Real(0.1))
        }};
        return m;
    }

    inline const Vector<Complex>& mat_cmplx_1_10x10_rhs0() {
        static const Vector<Complex> v{
            Complex(MML::Real(1.9),MML::Real(0.2)),
            Complex(MML::Real(1.0),-MML::Real(2.8)),
            Complex(MML::Real(2.8),MML::Real(0.4)),
            Complex(MML::Real(0.7),-MML::Real(1.2)),
            Complex(MML::Real(1.7),-MML::Real(2.3)),
            Complex(-MML::Real(1.1),MML::Real(2.9)),
            Complex(-MML::Real(2.6),MML::Real(0.8)),
            Complex(MML::Real(1.0),MML::Real(1.3)),
            Complex(MML::Real(1.8),MML::Real(1.4)),
            Complex(MML::Real(2.4),-MML::Real(2.3))
        };
        return v;
    }
    inline const Vector<Complex>& mat_cmplx_1_10x10_rhs0_sol() {
        static const Vector<Complex> v{
            Complex(MML::Real(1.424455431769854),-MML::Real(0.098032182842441)),
            Complex(-MML::Real(0.250718193209870),MML::Real(0.054061631825379)),
            Complex(-MML::Real(1.313751354614997),-MML::Real(0.054192017868504)),
            Complex(-MML::Real(0.219616343264554),-MML::Real(0.020504031689833)),
            Complex(MML::Real(0.341075078386610),MML::Real(0.201100460621572)),
            Complex(MML::Real(0.415000822479465),MML::Real(0.064838930474745)),
            Complex(MML::Real(0.320142292054508),MML::Real(0.568435781963075)),
            Complex(MML::Real(0.911687848983700),MML::Real(1.087865845194611)),
            Complex(MML::Real(0.626424756320989),MML::Real(0.415393058674079)),
            Complex(-MML::Real(0.338535794500345),-MML::Real(0.532388399034952))
        };
        return v;
    }

    // Matrix properties for mat_cmplx_1_10x10
    inline Real mat_cmplx_1_10x10_det_abs() { return MML::Real(1042200039.51399); }
    inline const Vector<Real>& mat_cmplx_1_10x10_singular_values() {
        static const Vector<Real> v{MML::Real(22.7185465473606), MML::Real(19.8152614906404), MML::Real(17.208921749824), MML::Real(14.7435183297016), MML::Real(11.5502006454006), MML::Real(10.2323923940024), MML::Real(7.627380102872), MML::Real(5.79022085629406), MML::Real(2.68282848064847), MML::Real(0.651607085870877)};
        return v;
    }
    inline Real mat_cmplx_1_10x10_cond_2() { return MML::Real(34.8654074517885); }
    inline int mat_cmplx_1_10x10_rank() { return 10; }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    inline const Matrix<Complex>& mat_cmplx_2_10x10() {
        static const Matrix<Complex> m{10, 10, {
            Complex(-MML::Real(1.1),-MML::Real(0.7)), Complex(MML::Real(2.2),MML::Real(0.3)), Complex(MML::Real(0.6),MML::Real(4.8)), Complex(MML::Real(2.7),-MML::Real(2.0)), Complex(MML::Real(4.4),-MML::Real(1.5)), Complex(-MML::Real(2.7),MML::Real(4.5)), Complex(-MML::Real(4.3),-MML::Real(3.4)), Complex(-MML::Real(3.0),-MML::Real(3.8)), Complex(MML::Real(2.9),-MML::Real(3.5)), Complex(-MML::Real(3.4),-MML::Real(2.4)),
            Complex(-MML::Real(0.7),MML::Real(3.5)), Complex(-MML::Real(0.6),MML::Real(3.7)), Complex(-MML::Real(0.9),MML::Real(3.7)), Complex(-MML::Real(1.5),-MML::Real(0.4)), Complex(MML::Real(3.0),MML::Real(1.0)), Complex(MML::Real(5.0),MML::Real(2.2)), Complex(MML::Real(4.5),-MML::Real(1.2)), Complex(-MML::Real(0.9),MML::Real(3.6)), Complex(MML::Real(3.2),MML::Real(0.0)), Complex(MML::Real(0.9),MML::Real(0.9)),
            Complex(MML::Real(4.3),MML::Real(1.6)), Complex(-MML::Real(2.7),-MML::Real(1.1)), Complex(-MML::Real(4.1),MML::Real(3.3)), Complex(MML::Real(2.1),-MML::Real(3.7)), Complex(-MML::Real(0.4),-MML::Real(2.2)), Complex(MML::Real(1.0),-MML::Real(0.2)), Complex(-MML::Real(4.5),-MML::Real(3.4)), Complex(-MML::Real(4.0),MML::Real(3.8)), Complex(-MML::Real(4.7),-MML::Real(0.6)), Complex(MML::Real(3.4),-MML::Real(4.2)),
            Complex(-MML::Real(2.8),MML::Real(2.8)), Complex(-MML::Real(4.1),-MML::Real(4.1)), Complex(MML::Real(2.0),MML::Real(0.7)), Complex(MML::Real(4.4),MML::Real(3.0)), Complex(MML::Real(0.5),MML::Real(1.9)), Complex(MML::Real(2.8),-MML::Real(4.7)), Complex(-MML::Real(4.0),MML::Real(0.2)), Complex(-MML::Real(0.8),-MML::Real(0.2)), Complex(MML::Real(3.7),-MML::Real(2.6)), Complex(MML::Real(2.0),-MML::Real(3.2)),
            Complex(MML::Real(0.1),MML::Real(1.8)), Complex(MML::Real(4.7),-MML::Real(2.7)), Complex(-MML::Real(4.1),MML::Real(4.1)), Complex(-MML::Real(0.3),MML::Real(2.4)), Complex(-MML::Real(4.2),-MML::Real(3.7)), Complex(-MML::Real(2.6),-MML::Real(1.8)), Complex(MML::Real(1.6),MML::Real(2.9)), Complex(-MML::Real(2.6),-MML::Real(2.3)), Complex(-MML::Real(1.7),MML::Real(3.6)), Complex(MML::Real(4.9),-MML::Real(0.9)),
            Complex(MML::Real(4.4),MML::Real(0.3)), Complex(MML::Real(4.7),MML::Real(1.2)), Complex(-MML::Real(1.4),MML::Real(4.3)), Complex(-MML::Real(2.5),MML::Real(4.1)), Complex(MML::Real(3.7),MML::Real(4.2)), Complex(MML::Real(4.9),MML::Real(0.7)), Complex(MML::Real(1.9),MML::Real(0.0)), Complex(-MML::Real(4.3),-MML::Real(2.7)), Complex(MML::Real(4.8),MML::Real(2.0)), Complex(MML::Real(1.9),-MML::Real(4.3)),
            Complex(MML::Real(2.0),MML::Real(4.0)), Complex(-MML::Real(0.9),-MML::Real(3.7)), Complex(-MML::Real(0.2),-MML::Real(0.9)), Complex(MML::Real(4.7),MML::Real(0.8)), Complex(MML::Real(1.3),-MML::Real(4.8)), Complex(-MML::Real(4.6),MML::Real(4.7)), Complex(MML::Real(0.4),-MML::Real(2.7)), Complex(-MML::Real(0.8),-MML::Real(0.5)), Complex(-MML::Real(0.5),MML::Real(1.6)), Complex(-MML::Real(2.6),-MML::Real(3.6)),
            Complex(MML::Real(0.9),-MML::Real(3.2)), Complex(-MML::Real(1.5),-MML::Real(3.8)), Complex(MML::Real(3.2),MML::Real(1.7)), Complex(MML::Real(4.9),-MML::Real(4.8)), Complex(-MML::Real(0.1),-MML::Real(2.2)), Complex(MML::Real(4.0),MML::Real(1.4)), Complex(-MML::Real(2.8),MML::Real(3.8)), Complex(MML::Real(1.8),MML::Real(1.9)), Complex(-MML::Real(1.3),MML::Real(0.4)), Complex(MML::Real(2.6),-MML::Real(0.9)),
            Complex(-MML::Real(2.8),MML::Real(3.5)), Complex(MML::Real(0.8),-MML::Real(2.2)), Complex(MML::Real(0.3),-MML::Real(4.5)), Complex(MML::Real(3.1),-MML::Real(2.5)), Complex(-MML::Real(4.1),MML::Real(4.5)), Complex(MML::Real(4.2),-MML::Real(2.8)), Complex(MML::Real(2.2),-MML::Real(2.8)), Complex(MML::Real(2.9),-MML::Real(3.4)), Complex(-MML::Real(3.2),MML::Real(0.1)), Complex(-MML::Real(4.8),MML::Real(4.5)),
            Complex(-MML::Real(2.6),MML::Real(3.9)), Complex(MML::Real(1.4),MML::Real(0.7)), Complex(MML::Real(2.3),MML::Real(3.8)), Complex(MML::Real(0.5),-MML::Real(3.5)), Complex(-MML::Real(3.1),-MML::Real(4.7)), Complex(MML::Real(3.9),-MML::Real(3.3)), Complex(MML::Real(4.1),-MML::Real(4.0)), Complex(-MML::Real(0.2),MML::Real(2.1)), Complex(-MML::Real(0.8),-MML::Real(1.3)), Complex(-MML::Real(3.7),MML::Real(0.8))
        }};
        return m;
    }

    inline const Vector<Complex>& mat_cmplx_2_10x10_rhs0() {
        static const Vector<Complex> v{
            Complex(MML::Real(1.3),-MML::Real(2.2)),
            Complex(-1.5,2.5),
            Complex(3.0,0.0),
            Complex(-MML::Real(0.6),-MML::Real(0.6)),
            Complex(MML::Real(1.3),MML::Real(1.1)),
            Complex(MML::Real(0.7),-MML::Real(0.2)),
            Complex(-MML::Real(2.9),MML::Real(0.5)),
            Complex(-MML::Real(2.5),MML::Real(2.6)),
            Complex(MML::Real(0.6),MML::Real(0.3)),
            Complex(-MML::Real(2.7),-MML::Real(0.7))
        };
        return v;
    }
    inline const Vector<Complex>& mat_cmplx_2_10x10_rhs0_sol() {
        static const Vector<Complex> v{
            Complex(MML::Real(0.706007116481964),MML::Real(0.013421144138452)),
            Complex(-MML::Real(0.668876690042133),MML::Real(1.976565088318023)),
            Complex(-MML::Real(1.550612209444719),-MML::Real(1.287669181025966)),
            Complex(-MML::Real(1.704881833429372),-MML::Real(0.054239536509271)),
            Complex(MML::Real(0.105269869191405),-MML::Real(0.404764634517515)),
            Complex(-MML::Real(0.009580978242681),-MML::Real(0.232977323759457)),
            Complex(MML::Real(0.500618952349040),-MML::Real(0.372903391127123)),
            Complex(MML::Real(1.651960217269608),MML::Real(0.563573074178541)),
            Complex(-MML::Real(0.457765142442384),MML::Real(2.177776204331938)),
            Complex(-MML::Real(0.596210857708472),MML::Real(0.239850315657696))
        };
        return v;
    }

    // Matrix properties for mat_cmplx_2_10x10
    inline Real mat_cmplx_2_10x10_det_abs() { return MML::Real(3907151884.60909); }
    inline const Vector<Real>& mat_cmplx_2_10x10_singular_values() {
        static const Vector<Real> v{MML::Real(22.3486893729876), MML::Real(18.4057371396619), MML::Real(16.0143867894748), MML::Real(14.9306536661108), MML::Real(13.7285234336414), MML::Real(9.71270306087155), MML::Real(8.32484844148929), MML::Real(7.14073463579936), MML::Real(6.0902048638322), MML::Real(0.82290710430045)};
        return v;
    }
    inline Real mat_cmplx_2_10x10_cond_2() { return MML::Real(27.1582165911499); }
    inline int mat_cmplx_2_10x10_rank() { return 10; }

    //////////////////////////////////////////////////////////////////////////////////////////////////
    inline const Matrix<Complex>& mat_cmplx_3_10x10() {
        static const Matrix<Complex> m{10, 10, {
            Complex(-MML::Real(2.5),MML::Real(1.0)), Complex(MML::Real(3.8),-MML::Real(3.9)), Complex(MML::Real(4.4),-MML::Real(1.8)), Complex(-MML::Real(1.9),-MML::Real(2.3)), Complex(-MML::Real(0.3),-MML::Real(1.7)), Complex(-MML::Real(3.3),-MML::Real(3.9)), Complex(-MML::Real(2.8),MML::Real(0.9)), Complex(-MML::Real(2.2),MML::Real(4.9)), Complex(MML::Real(1.2),MML::Real(1.4)), Complex(-MML::Real(3.0),MML::Real(4.7)),
            Complex(MML::Real(3.6),-MML::Real(1.3)), Complex(-MML::Real(1.0),MML::Real(1.9)), Complex(MML::Real(2.0),MML::Real(3.1)), Complex(-MML::Real(0.3),-MML::Real(4.8)), Complex(-MML::Real(3.9),MML::Real(0.7)), Complex(MML::Real(0.2),-MML::Real(3.5)), Complex(MML::Real(3.9),-MML::Real(0.6)), Complex(-MML::Real(2.1),-MML::Real(2.4)), Complex(MML::Real(3.1),-MML::Real(1.9)), Complex(-MML::Real(4.8),MML::Real(1.5)),
            Complex(MML::Real(2.7),-MML::Real(3.6)), Complex(MML::Real(1.4),MML::Real(3.7)), Complex(MML::Real(1.0),MML::Real(0.3)), Complex(-MML::Real(0.1),-MML::Real(3.4)), Complex(MML::Real(2.6),-MML::Real(2.8)), Complex(-MML::Real(2.6),-MML::Real(4.5)), Complex(MML::Real(4.7),-MML::Real(4.9)), Complex(MML::Real(2.6),MML::Real(0.3)), Complex(-MML::Real(4.5),-MML::Real(1.6)), Complex(-MML::Real(4.3),MML::Real(0.4)),
            Complex(-MML::Real(1.7),-MML::Real(3.7)), Complex(-MML::Real(1.8),MML::Real(5.0)), Complex(-MML::Real(3.2),MML::Real(3.6)), Complex(-MML::Real(0.4),-MML::Real(2.8)), Complex(MML::Real(2.3),MML::Real(3.0)), Complex(MML::Real(4.2),MML::Real(0.4)), Complex(-MML::Real(4.9),-MML::Real(3.6)), Complex(-MML::Real(4.8),-MML::Real(0.2)), Complex(-MML::Real(2.5),-MML::Real(3.7)), Complex(-MML::Real(0.2),MML::Real(3.7)),
            Complex(-MML::Real(4.0),MML::Real(2.5)), Complex(MML::Real(0.9),-MML::Real(0.3)), Complex(-MML::Real(2.2),-MML::Real(0.5)), Complex(MML::Real(1.5),MML::Real(3.2)), Complex(-MML::Real(4.0),MML::Real(4.5)), Complex(MML::Real(1.9),MML::Real(0.7)), Complex(-MML::Real(4.4),MML::Real(3.5)), Complex(MML::Real(2.7),-MML::Real(0.4)), Complex(-MML::Real(0.5),MML::Real(3.5)), Complex(-MML::Real(1.8),MML::Real(1.3)),
            Complex(MML::Real(1.7),MML::Real(0.7)), Complex(MML::Real(2.1),MML::Real(0.6)), Complex(MML::Real(0.1),-MML::Real(3.8)), Complex(-MML::Real(3.0),-MML::Real(3.9)), Complex(MML::Real(4.5),-MML::Real(2.9)), Complex(-MML::Real(1.4),MML::Real(0.5)), Complex(MML::Real(0.9),-MML::Real(3.6)), Complex(MML::Real(4.7),-MML::Real(1.5)), Complex(MML::Real(2.8),MML::Real(0.2)), Complex(MML::Real(5.0),MML::Real(2.4)),
            Complex(-MML::Real(3.7),MML::Real(3.9)), Complex(-MML::Real(0.8),-MML::Real(1.6)), Complex(MML::Real(3.3),MML::Real(0.6)), Complex(-MML::Real(0.7),-MML::Real(2.3)), Complex(MML::Real(3.3),MML::Real(4.1)), Complex(MML::Real(1.0),-MML::Real(1.0)), Complex(MML::Real(0.5),-MML::Real(2.2)), Complex(MML::Real(1.0),MML::Real(2.5)), Complex(MML::Real(2.3),MML::Real(4.0)), Complex(MML::Real(1.1),-MML::Real(4.4)),
            Complex(MML::Real(0.5),MML::Real(4.4)), Complex(MML::Real(4.8),-MML::Real(3.8)), Complex(MML::Real(2.0),-MML::Real(3.5)), Complex(-MML::Real(2.5),-MML::Real(2.4)), Complex(MML::Real(1.7),MML::Real(1.1)), Complex(-MML::Real(3.6),MML::Real(4.0)), Complex(MML::Real(1.4),MML::Real(2.5)), Complex(MML::Real(2.4),-MML::Real(0.9)), Complex(-MML::Real(0.9),-MML::Real(4.9)), Complex(-MML::Real(3.6),MML::Real(4.6)),
            Complex(-MML::Real(2.6),MML::Real(2.3)), Complex(MML::Real(4.4),-MML::Real(1.5)), Complex(MML::Real(4.4),-MML::Real(4.5)), Complex(-MML::Real(0.7),MML::Real(3.6)), Complex(-MML::Real(4.7),MML::Real(1.0)), Complex(MML::Real(0.7),-MML::Real(2.7)), Complex(-MML::Real(4.6),MML::Real(5.0)), Complex(MML::Real(0.9),-MML::Real(4.3)), Complex(MML::Real(3.9),-MML::Real(4.1)), Complex(-MML::Real(0.1),-MML::Real(4.8)),
            Complex(MML::Real(4.4),MML::Real(3.1)), Complex(-MML::Real(0.7),MML::Real(1.6)), Complex(-MML::Real(4.4),MML::Real(1.6)), Complex(MML::Real(3.4),-MML::Real(3.3)), Complex(MML::Real(2.6),-MML::Real(1.6)), Complex(MML::Real(0.3),MML::Real(3.8)), Complex(-MML::Real(3.1),MML::Real(1.8)), Complex(MML::Real(4.9),-MML::Real(2.5)), Complex(-MML::Real(2.9),-MML::Real(0.5)), Complex(MML::Real(0.0),MML::Real(1.5))
        }};
        return m;
    }

    inline const Vector<Complex>& mat_cmplx_3_10x10_rhs0() {
        static const Vector<Complex> v{
            Complex(-MML::Real(0.0),MML::Real(0.6)),
            Complex(-MML::Real(1.9),MML::Real(1.8)),
            Complex(MML::Real(2.2),MML::Real(1.8)),
            Complex(-MML::Real(2.9),MML::Real(0.4)),
            Complex(MML::Real(1.5),-MML::Real(1.8)),
            Complex(MML::Real(2.6),MML::Real(2.1)),
            Complex(MML::Real(0.8),MML::Real(2.6)),
            Complex(1.0,-0.5),
            Complex(-MML::Real(1.3),-MML::Real(0.2)),
            Complex(-MML::Real(2.4),-MML::Real(2.3))
        };
        return v;
    }
    inline const Vector<Complex>& mat_cmplx_3_10x10_rhs0_sol() {
        static const Vector<Complex> v{
            Complex(MML::Real(0.205631403861014),-MML::Real(0.069120662300960)),
            Complex(MML::Real(0.607552380277906),-MML::Real(0.133279742702262)),
            Complex(-MML::Real(0.627691172799875),MML::Real(0.093459414320121)),
            Complex(-MML::Real(0.588592081392438),-MML::Real(0.390493395743841)),
            Complex(MML::Real(0.051256231959498),-MML::Real(0.186957099902922)),
            Complex(-MML::Real(0.059106269104728),MML::Real(0.076716221074585)),
            Complex(MML::Real(0.226099826745400),MML::Real(0.256171743798648)),
            Complex(MML::Real(0.167444522224238),-MML::Real(0.330821725681724)),
            Complex(MML::Real(0.211147906212214),-MML::Real(0.251467714908403)),
            Complex(-MML::Real(0.216706429129029),MML::Real(0.165134497920343))
        };
        return v;
    }

    // Matrix properties for mat_cmplx_3_10x10
    inline Real mat_cmplx_3_10x10_det_abs() { return MML::Real(2412748639.16755); }
    inline const Vector<Real>& mat_cmplx_3_10x10_singular_values() {
        static const Vector<Real> v{MML::Real(22.5163950951509), MML::Real(20.275670042204), MML::Real(17.0822720041155), MML::Real(13.923744195654), MML::Real(12.0162310557434), MML::Real(9.04808201323421), MML::Real(6.67171858689903), MML::Real(4.89006586995963), MML::Real(3.45262581668275), MML::Real(1.81430205907735)};
        return v;
    }
    inline Real mat_cmplx_3_10x10_cond_2() { return MML::Real(12.410499664318); }
    inline int mat_cmplx_3_10x10_rank() { return 10; }
}
#endif
