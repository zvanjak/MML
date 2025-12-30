#if !defined MML_LINEAR_ALG_EQ_SYSTEMS_COMPLEX_DEFS_H
#define MML_LINEAR_ALG_EQ_SYSTEMS_COMPLEX_DEFS_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Matrix.h"
#endif

namespace MML::TestBeds
{
    /***********************************************************************************************/
    /**************                     3 x 3 test matrices                           **************/
    /***********************************************************************************************/
    const static inline MML::Matrix<Complex>     mat_cmplx_3x3{3, 3, { 1.0, 2.0, -1.0, 
                                                                      -1.0, 5.0, 6.0, 
                                                                       3.0, 1.0, 1.0 }};

    const static inline MML::Vector<Complex>     mat_cmplx_3x3_rhs0{1.0, 2.0, 1.0};
    const static inline MML::Vector<Complex>     mat_cmplx_3x3_rhs0_sol{0.18867924528301885, 0.41509433962264153, 0.018867924528301921};

    // Matrix properties for mat_cmplx_3x3 (same as real mat_3x3)
    const static inline Real                     mat_cmplx_3x3_det_abs = 53.0;  // |det(A)|
    const static inline MML::Vector<Real>        mat_cmplx_3x3_singular_values{7.965662078134766, 3.4347838553531793, 1.9371568814853632};
    const static inline Real                     mat_cmplx_3x3_cond_2 = 4.112039181691377;
    const static inline int                      mat_cmplx_3x3_rank = 3;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Complex>     mat_cmplx_1_3x3{3, 3, { Complex(1.0, 2.0), Complex(2.0, -0.5), Complex(-1.0, 2.3), 
                                                                         Complex(-2.0,1.0), Complex(4.0,-2.0),  Complex(6.0,5.0), 
                                                                         Complex(2.5,1.0),  Complex(1.0,-1.0),  Complex(-2.0, -4.0) }};

    const static inline MML::Vector<Complex>     mat_cmplx_1_3x3_rhs0{ Complex(1.5,0.5), Complex(-2.0,1.5), Complex(-1.0,-2.0)};
    const static inline MML::Vector<Complex>     mat_cmplx_1_3x3_rhs0_sol{ Complex(0.693423783306577,-1.25000647849031),
                                                                            Complex(-1.09549140996106,-0.596939956473329),
                                                                            Complex(0.242505265145942,-0.451531010007621)};

    // Matrix properties for mat_cmplx_1_3x3
    const static inline Real                     mat_cmplx_1_3x3_det_abs = 144.36547296686918;  // |det(A)|
    const static inline MML::Vector<Real>        mat_cmplx_1_3x3_singular_values{10.20912316755949, 5.028959556700652, 2.815033298093671};
    const static inline Real                     mat_cmplx_1_3x3_cond_2 = 3.6263579631689193;
    const static inline int                      mat_cmplx_1_3x3_rank = 3;

    /////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Complex>     mat_cmplx_1_5x5{5, 5, { Complex(1.0,2.0),  Complex(2.0,-0.5), Complex(-1.0, 2.3), Complex(1.8,-3.0),  Complex(-2.5,3.7), 
                                                                         Complex(-2.0,1.0), Complex(4.0,-2.0), Complex(6.0,5.0),   Complex(-1.0,-2.0), Complex(-1.0,-2.0), 
                                                                         Complex(3.0,-3.0), Complex(1.0,1.4),  Complex(2.0,-1.0),  Complex(-1.7,-2.6), Complex(7.0,-3.0), 
                                                                         Complex(-4.0,1.5), Complex(2.0,-3.0), Complex(-2.0,3.0),  Complex(1.0,2.0),   Complex(-5.0,6.0), 
                                                                         Complex(2.5,1.0),  Complex(1.0,-1.0), Complex(-2.0,-4.0), Complex(-1.0,2.0),  Complex(2.0,3.0) }};

    const static inline MML::Vector<Complex>     mat_cmplx_1_5x5_rhs0{ Complex(1.5,0.5), Complex(-2.0,1.5), Complex(-1.0,-2.0), Complex(1.0,-2.5), Complex(-1.8,2.1)};
    const static inline MML::Vector<Complex>     mat_cmplx_1_5x5_rhs0_sol{ Complex(0.901430270043644,1.58987437993965),
                                                                            Complex(-0.351975602660937,0.144870096280711),
                                                                            Complex(0.558130434051583,-0.123040787906678),
                                                                            Complex(0.894179352753342,-0.202080868658524),
                                                                            Complex(-0.780553254235851,-0.467650784949547)};

    // Matrix properties for mat_cmplx_1_5x5
    const static inline Real                     mat_cmplx_1_5x5_det_abs = 10168.952903809948;  // |det(A)|
    const static inline MML::Vector<Real>        mat_cmplx_1_5x5_singular_values{15.64174282618698, 9.00912568656689, 5.696377816847538, 4.137858398839877, 2.380802820168614};
    const static inline Real                     mat_cmplx_1_5x5_cond_2 = 6.570844779227217;
    const static inline int                      mat_cmplx_1_5x5_rank = 5;

    //////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Complex>     mat_cmplx_2_5x5{5, 5, {
            Complex(-4.7,-1.0), Complex(-0.1,-3.4), Complex(3.5,0.9), Complex(-0.9,-1.2), Complex(1.3,-2.1),
            Complex(-2.1,-4.8), Complex(0.4,-0.3), Complex(-2.7,4.1), Complex(-2.6,-4.0), Complex(4.5,-0.7),
            Complex(1.5,-3.7), Complex(-4.8,-3.9), Complex(5.0,4.3), Complex(2.4,1.2), Complex(2.2,-2.3),
            Complex(-0.9,2.3), Complex(-3.0,-3.1), Complex(4.6,-2.2), Complex(-2.8,-2.1), Complex(-4.5,0.2),
            Complex(3.8,-1.1), Complex(3.9,-2.7), Complex(4.6,3.1), Complex(1.1,2.1), Complex(4.9,3.5)
        }};

    const static inline MML::Vector<Complex>     mat_cmplx_2_5x5_rhs0{
            Complex(1.1,2.7),
            Complex(0.1,2.9),
            Complex(0.7,2.1),
            Complex(-0.1,0.5),
            Complex(-0.2,-1.8)
        };
    const static inline MML::Vector<Complex>     mat_cmplx_2_5x5_rhs0_sol{
            Complex(-0.364324249141192,-0.683814875129271),
            Complex(-0.422202239716293,-0.369635385572325),
            Complex(0.189016352402793,-0.296292202554985),
            Complex(0.145166170094375,0.115685130841116),
            Complex(0.489514812714704,-0.136538611915636)
        };

    // Matrix properties for mat_cmplx_2_5x5
    const static inline Real                     mat_cmplx_2_5x5_det_abs = 13745.1507590657;
    const static inline MML::Vector<Real>        mat_cmplx_2_5x5_singular_values{16.1973119257448, 10.7879205362089, 6.40784516509668, 6.0423260356214, 2.03166761650746};
    const static inline Real                     mat_cmplx_2_5x5_cond_2 = 7.97242215908764;
    const static inline int                      mat_cmplx_2_5x5_rank = 5;

    //////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Complex>     mat_cmplx_3_5x5{5, 5, {
            Complex(-3.9,-1.8), Complex(1.8,-4.3), Complex(0.3,-4.5), Complex(-1.3,-0.6), Complex(-0.9,2.0),
            Complex(0.9,2.9), Complex(2.2,-2.6), Complex(-3.3,1.5), Complex(-4.3,3.7), Complex(0.0,-2.9),
            Complex(-0.7,-3.7), Complex(4.5,3.1), Complex(-0.8,-4.6), Complex(-2.5,-1.4), Complex(-1.0,1.8),
            Complex(1.3,3.3), Complex(0.8,-0.6), Complex(2.8,-4.2), Complex(1.2,-4.1), Complex(1.3,1.9),
            Complex(-1.3,2.9), Complex(-3.3,-1.8), Complex(2.5,-1.2), Complex(-1.0,2.3), Complex(1.6,-2.2)
        }};

    const static inline MML::Vector<Complex>     mat_cmplx_3_5x5_rhs0{
            Complex(-2.2,2.5),
            Complex(-1.1,-0.7),
            Complex(0.3,-2.2),
            Complex(-2.2,3.0),
            Complex(-0.0,1.6)
        };
    const static inline MML::Vector<Complex>     mat_cmplx_3_5x5_rhs0_sol{
            Complex(0.061637980467496,0.231305484458262),
            Complex(-0.608924639903936,-0.383281231895997),
            Complex(-0.283005443992466,-0.007473368712111),
            Complex(-0.291587621773846,0.038722221408967),
            Complex(0.062649872788515,-0.070415668769379)
        };

    // Matrix properties for mat_cmplx_3_5x5
    const static inline Real                     mat_cmplx_3_5x5_det_abs = 6319.66771674941;
    const static inline MML::Vector<Real>        mat_cmplx_3_5x5_singular_values{13.3205460279412, 9.39155667384619, 6.02046598723545, 4.85473027591569, 1.72838074627172};
    const static inline Real                     mat_cmplx_3_5x5_cond_2 = 7.70695117767013;
    const static inline int                      mat_cmplx_3_5x5_rank = 5;

    //////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Complex>     mat_cmplx_4_5x5{5, 5, {
            Complex(-2.2,-4.9), Complex(1.4,-0.6), Complex(1.2,-4.5), Complex(2.8,1.7), Complex(-0.2,-4.7),
            Complex(4.6,-1.2), Complex(5.0,3.2), Complex(3.3,4.0), Complex(-3.0,-1.7), Complex(-3.7,-4.0),
            Complex(1.2,-3.1), Complex(0.9,-3.4), Complex(1.4,4.4), Complex(4.9,4.9), Complex(-1.6,-4.1),
            Complex(3.1,-2.4), Complex(1.7,-3.7), Complex(0.4,2.2), Complex(-4.3,-3.1), Complex(-4.1,-4.8),
            Complex(3.7,-3.5), Complex(0.9,4.0), Complex(-4.2,1.1), Complex(-1.6,-4.0), Complex(-1.3,-3.3)
        }};

    const static inline MML::Vector<Complex>     mat_cmplx_4_5x5_rhs0{
            Complex(1.5,-1.1),
            Complex(-2.3,-2.4),
            Complex(-3.0,0.5),
            Complex(-1.2,-1.8),
            Complex(-3.0,2.5)
        };
    const static inline MML::Vector<Complex>     mat_cmplx_4_5x5_rhs0_sol{
            Complex(-0.926407716133575,0.033047493912236),
            Complex(-0.055819119213253,0.016480814215369),
            Complex(0.411015404552370,-0.188896066902875),
            Complex(-0.400191856408234,0.311772178213589),
            Complex(0.781799961064350,0.325947296047793)
        };

    // Matrix properties for mat_cmplx_4_5x5
    const static inline Real                     mat_cmplx_4_5x5_det_abs = 35140.2960834746;
    const static inline MML::Vector<Real>        mat_cmplx_4_5x5_singular_values{16.4968668117334, 11.4058088404484, 8.20284126627676, 6.24295130117635, 3.64689790129684};
    const static inline Real                     mat_cmplx_4_5x5_cond_2 = 4.52353404406171;
    const static inline int                      mat_cmplx_4_5x5_rank = 5;

    /***********************************************************************************************/
    /**************                     8 x 8 test matrices                           **************/
    /***********************************************************************************************/
    const static inline MML::Matrix<Complex>     mat_cmplx_1_8x8{8, 8, {
            Complex(4.5,4.8), Complex(2.4,-0.5), Complex(-4.6,4.0), Complex(-1.7,1.1), Complex(-1.8,1.4), Complex(-0.9,3.7), Complex(-4.3,4.4), Complex(-1.7,-1.2),
            Complex(2.8,3.0), Complex(-3.3,4.6), Complex(0.3,-1.3), Complex(-2.5,-0.8), Complex(-4.4,-1.1), Complex(-3.1,0.4), Complex(-0.6,-2.7), Complex(1.5,2.7),
            Complex(-1.0,-1.8), Complex(4.1,3.2), Complex(1.6,-3.9), Complex(2.0,0.8), Complex(0.4,-0.5), Complex(-0.8,0.5), Complex(-4.0,-0.1), Complex(-1.2,2.0),
            Complex(-1.9,-0.9), Complex(4.9,-2.0), Complex(4.2,-1.2), Complex(-2.0,-1.6), Complex(-3.0,4.1), Complex(-2.2,4.1), Complex(-4.4,0.9), Complex(1.7,-2.5),
            Complex(-2.1,3.3), Complex(3.2,-3.6), Complex(3.3,3.3), Complex(-1.5,-3.9), Complex(0.3,4.8), Complex(-0.3,2.2), Complex(-2.1,-4.2), Complex(0.7,-1.6),
            Complex(-1.4,4.2), Complex(-1.4,1.0), Complex(-4.2,2.6), Complex(-3.8,3.9), Complex(2.2,-3.4), Complex(-4.5,2.7), Complex(-2.6,0.3), Complex(-3.5,1.3),
            Complex(-2.8,3.0), Complex(2.7,3.2), Complex(-3.7,4.0), Complex(1.6,-0.7), Complex(-4.8,0.8), Complex(-3.4,-1.4), Complex(0.3,-3.0), Complex(-1.7,4.2),
            Complex(1.8,2.8), Complex(4.9,-4.7), Complex(3.0,-4.8), Complex(-4.8,-1.6), Complex(4.8,4.3), Complex(-1.1,0.4), Complex(-0.1,-3.8), Complex(-2.9,1.5)
        }};

    const static inline MML::Vector<Complex>     mat_cmplx_1_8x8_rhs0{
            Complex(-2.4,-0.6),
            Complex(-1.2,1.5),
            Complex(1.7,-2.0),
            Complex(0.0,2.7),
            Complex(-2.3,-2.7),
            Complex(-2.4,2.2),
            Complex(2.2,0.4),
            Complex(2.2,2.3)
        };
    const static inline MML::Vector<Complex>     mat_cmplx_1_8x8_rhs0_sol{
            Complex(0.279983758960108,-0.280318189874449),
            Complex(-0.185720534820703,0.661437184412004),
            Complex(0.080737838282574,0.657637228492826),
            Complex(1.476088586253125,-1.410794553549068),
            Complex(-0.290319020357946,-0.864416476092015),
            Complex(-1.032934559655212,0.267746969846193),
            Complex(0.293628372974839,0.398814739983594),
            Complex(-0.066786432249243,0.250909440756107)
        };

    // Matrix properties for mat_cmplx_1_8x8
    const static inline Real                     mat_cmplx_1_8x8_det_abs = 8254361.91150796;
    const static inline MML::Vector<Real>        mat_cmplx_1_8x8_singular_values{20.7329056394107, 17.1020495055299, 12.3211712343142, 9.14017764485158, 7.61605888733823, 5.18445288908348, 3.62879622218723, 1.44268906430648};
    const static inline Real                     mat_cmplx_1_8x8_cond_2 = 14.3710146228753;
    const static inline int                      mat_cmplx_1_8x8_rank = 8;

    //////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Complex>     mat_cmplx_2_8x8{8, 8, {
            Complex(-1.8,0.4), Complex(2.6,3.1), Complex(3.8,0.8), Complex(1.9,-0.6), Complex(2.8,1.7), Complex(4.1,-1.1), Complex(2.3,4.8), Complex(-2.0,1.1),
            Complex(-2.2,-1.7), Complex(2.1,-1.8), Complex(0.4,-0.3), Complex(5.0,4.7), Complex(4.7,0.1), Complex(3.5,3.2), Complex(-1.9,0.1), Complex(1.8,-1.4),
            Complex(-4.4,3.7), Complex(0.2,-2.9), Complex(-3.1,2.0), Complex(4.3,0.7), Complex(5.0,4.3), Complex(-4.0,4.1), Complex(-1.2,4.7), Complex(3.7,-3.4),
            Complex(3.9,1.5), Complex(-3.7,2.9), Complex(0.9,3.6), Complex(0.4,4.2), Complex(4.2,-0.1), Complex(0.1,-4.9), Complex(-2.6,-1.6), Complex(0.7,1.6),
            Complex(1.8,4.2), Complex(-0.5,-1.6), Complex(-2.6,-4.9), Complex(1.0,-3.0), Complex(-0.0,-2.4), Complex(-0.6,-3.0), Complex(-1.3,-0.4), Complex(-1.4,0.4),
            Complex(-0.8,4.0), Complex(0.6,4.9), Complex(3.8,2.3), Complex(1.9,-4.1), Complex(-1.9,3.3), Complex(-3.6,1.0), Complex(-2.8,-4.0), Complex(3.0,1.7),
            Complex(1.4,3.0), Complex(-4.7,-1.3), Complex(-4.6,1.3), Complex(-3.2,-4.2), Complex(0.4,4.9), Complex(1.9,1.8), Complex(0.0,-1.2), Complex(-1.8,-1.9),
            Complex(-1.9,3.7), Complex(-1.2,-0.2), Complex(-1.5,3.1), Complex(-3.4,-2.5), Complex(-1.0,2.9), Complex(1.3,1.3), Complex(-1.7,-4.3), Complex(4.5,-3.7)
        }};

    const static inline MML::Vector<Complex>     mat_cmplx_2_8x8_rhs0{
            Complex(-0.8,0.5),
            Complex(-0.8,-0.5),
            Complex(0.8,1.8),
            Complex(-1.6,-0.6),
            Complex(-1.0,1.0),
            Complex(1.1,2.7),
            Complex(-1.2,0.9),
            Complex(2.6,1.9)
        };
    const static inline MML::Vector<Complex>     mat_cmplx_2_8x8_rhs0_sol{
            Complex(-0.071230196831102,-0.301989480770654),
            Complex(0.579940424169923,0.349107697938878),
            Complex(0.082682783590385,-0.579849142434740),
            Complex(-0.357408952376676,0.047989986784227),
            Complex(-0.090219498196085,0.319245310867363),
            Complex(-0.095237471765996,-0.142075932285791),
            Complex(-0.236103821389095,-0.005801913163186),
            Complex(-0.098903720964927,-0.051871533889682)
        };

    // Matrix properties for mat_cmplx_2_8x8
    const static inline Real                     mat_cmplx_2_8x8_det_abs = 9681017.30839477;
    const static inline MML::Vector<Real>        mat_cmplx_2_8x8_singular_values{17.9783905280538, 17.0468894260995, 13.0089171954949, 10.7846443438304, 8.43566113581774, 5.55392678710442, 2.69870097150435, 1.78075439196711};
    const static inline Real                     mat_cmplx_2_8x8_cond_2 = 10.0959405795394;
    const static inline int                      mat_cmplx_2_8x8_rank = 8;

    //////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Complex>     mat_cmplx_3_8x8{8, 8, {
            Complex(-3.5,2.4), Complex(-3.8,-3.1), Complex(-1.8,-4.5), Complex(1.8,1.6), Complex(1.1,-2.6), Complex(3.8,3.4), Complex(-1.9,-1.3), Complex(-1.8,3.9),
            Complex(-2.5,4.4), Complex(-0.4,-0.9), Complex(1.9,-0.8), Complex(2.2,-2.2), Complex(-3.5,-1.3), Complex(1.8,-3.1), Complex(-3.5,-2.2), Complex(2.8,4.0),
            Complex(-4.2,-3.2), Complex(-4.1,-2.8), Complex(-3.4,-2.2), Complex(-4.3,-4.8), Complex(0.2,-4.1), Complex(4.2,-0.2), Complex(0.6,-2.0), Complex(0.0,0.0),
            Complex(-2.3,3.9), Complex(4.0,1.1), Complex(-4.1,-4.5), Complex(-0.2,1.5), Complex(-0.3,-3.3), Complex(2.1,-0.3), Complex(4.1,4.1), Complex(-1.0,1.1),
            Complex(2.5,2.9), Complex(-1.3,-3.1), Complex(0.1,-3.7), Complex(2.0,3.0), Complex(-2.1,0.9), Complex(2.9,0.0), Complex(-3.5,-1.8), Complex(1.8,4.4),
            Complex(0.6,-3.0), Complex(4.7,-1.2), Complex(2.5,4.0), Complex(-4.4,-0.9), Complex(-1.9,4.8), Complex(0.7,0.2), Complex(-3.6,1.1), Complex(-0.7,0.5),
            Complex(2.0,1.8), Complex(2.2,0.0), Complex(4.4,1.7), Complex(4.6,1.9), Complex(-4.7,-1.4), Complex(-4.1,1.9), Complex(-1.5,2.6), Complex(-1.6,1.5),
            Complex(0.4,3.1), Complex(-3.7,-0.8), Complex(-1.9,-3.4), Complex(0.6,0.3), Complex(3.7,1.0), Complex(3.6,0.0), Complex(-0.6,0.0), Complex(1.3,4.2)
        }};

    const static inline MML::Vector<Complex>     mat_cmplx_3_8x8_rhs0{
            Complex(0.8,1.1),
            Complex(2.3,2.7),
            Complex(-1.4,0.0),
            Complex(-1.1,0.8),
            Complex(2.6,2.2),
            Complex(1.5,2.1),
            Complex(-3.0,1.4),
            Complex(-2.9,0.9)
        };
    const static inline MML::Vector<Complex>     mat_cmplx_3_8x8_rhs0_sol{
            Complex(-0.639415913735030,-0.641550528548533),
            Complex(0.468937846956392,0.485780618653904),
            Complex(0.162250528053269,-0.458584525975030),
            Complex(0.403585883119101,-0.615326804902230),
            Complex(-0.225178977849940,-0.170686259188212),
            Complex(1.311877929063035,-1.146147631413710),
            Complex(-0.351915453593622,0.227370481021649),
            Complex(1.839401613755739,2.025842111839458)
        };

    // Matrix properties for mat_cmplx_3_8x8
    const static inline Real                     mat_cmplx_3_8x8_det_abs = 5323054.71961043;
    const static inline MML::Vector<Real>        mat_cmplx_3_8x8_singular_values{19.1813429563283, 15.3590286176189, 12.4259062472838, 9.69670848034898, 6.7671725337737, 5.44102025627927, 3.49898640472805, 1.16395425345961};
    const static inline Real                     mat_cmplx_3_8x8_cond_2 = 16.4794646347274;
    const static inline int                      mat_cmplx_3_8x8_rank = 8;

    /***********************************************************************************************/
    /**************                    10 x 10 test matrices                          **************/
    /***********************************************************************************************/
    const static inline MML::Matrix<Complex>     mat_cmplx_1_10x10{10, 10, {
            Complex(-4.9,0.3), Complex(-0.7,3.7), Complex(-2.5,0.5), Complex(0.1,-0.7), Complex(2.5,2.5), Complex(-1.5,3.7), Complex(-2.2,3.8), Complex(4.0,-3.3), Complex(2.1,-3.0), Complex(4.0,1.0),
            Complex(-0.5,1.7), Complex(-3.7,-2.4), Complex(-4.9,3.5), Complex(4.0,5.0), Complex(3.7,1.2), Complex(-4.7,0.2), Complex(-1.3,2.9), Complex(-0.7,0.8), Complex(-1.5,0.5), Complex(1.0,1.2),
            Complex(-2.6,3.1), Complex(-3.4,3.7), Complex(-0.9,-4.4), Complex(-1.8,2.6), Complex(-4.3,-2.0), Complex(-4.8,-4.2), Complex(1.0,1.2), Complex(1.5,-4.4), Complex(-0.2,-3.9), Complex(1.5,0.3),
            Complex(1.2,-4.9), Complex(2.8,0.7), Complex(-2.3,2.8), Complex(3.3,1.8), Complex(4.4,-0.8), Complex(-2.3,-0.1), Complex(2.0,-2.0), Complex(1.1,4.4), Complex(3.2,1.2), Complex(-0.2,-4.8),
            Complex(-2.8,-1.8), Complex(3.9,4.1), Complex(-0.7,-2.7), Complex(-0.9,0.0), Complex(3.5,-5.0), Complex(-0.0,-3.6), Complex(3.0,4.3), Complex(-3.0,-4.4), Complex(4.5,4.0), Complex(-1.0,3.0),
            Complex(3.2,-3.2), Complex(3.0,-3.6), Complex(-2.9,3.0), Complex(-3.2,-3.9), Complex(-1.8,2.9), Complex(-4.5,3.2), Complex(-0.6,2.7), Complex(3.7,4.5), Complex(-3.9,2.3), Complex(0.3,1.2),
            Complex(-4.5,1.6), Complex(-4.9,-4.1), Complex(4.9,1.3), Complex(2.4,-2.9), Complex(4.8,4.7), Complex(4.8,4.4), Complex(3.2,-1.7), Complex(-0.5,-4.3), Complex(-3.1,-4.8), Complex(-2.4,-0.6),
            Complex(2.1,3.4), Complex(3.8,-4.0), Complex(0.5,1.9), Complex(-4.2,4.3), Complex(-1.0,-4.2), Complex(2.8,0.1), Complex(-1.7,0.5), Complex(-3.0,4.3), Complex(4.4,0.0), Complex(0.4,3.6),
            Complex(4.4,0.0), Complex(1.2,3.7), Complex(-0.3,-2.1), Complex(1.9,4.1), Complex(1.2,0.2), Complex(0.0,0.4), Complex(-3.5,3.0), Complex(1.2,1.2), Complex(2.9,-2.3), Complex(4.4,-4.7),
            Complex(1.8,-3.0), Complex(2.6,3.8), Complex(2.9,3.5), Complex(-4.9,0.0), Complex(2.2,-3.7), Complex(-1.6,-2.2), Complex(3.9,0.6), Complex(0.6,1.8), Complex(4.4,2.0), Complex(-3.2,0.1)
        }};

    const static inline MML::Vector<Complex>     mat_cmplx_1_10x10_rhs0{
            Complex(1.9,0.2),
            Complex(1.0,-2.8),
            Complex(2.8,0.4),
            Complex(0.7,-1.2),
            Complex(1.7,-2.3),
            Complex(-1.1,2.9),
            Complex(-2.6,0.8),
            Complex(1.0,1.3),
            Complex(1.8,1.4),
            Complex(2.4,-2.3)
        };
    const static inline MML::Vector<Complex>     mat_cmplx_1_10x10_rhs0_sol{
            Complex(1.424455431769854,-0.098032182842441),
            Complex(-0.250718193209870,0.054061631825379),
            Complex(-1.313751354614997,-0.054192017868504),
            Complex(-0.219616343264554,-0.020504031689833),
            Complex(0.341075078386610,0.201100460621572),
            Complex(0.415000822479465,0.064838930474745),
            Complex(0.320142292054508,0.568435781963075),
            Complex(0.911687848983700,1.087865845194611),
            Complex(0.626424756320989,0.415393058674079),
            Complex(-0.338535794500345,-0.532388399034952)
        };

    // Matrix properties for mat_cmplx_1_10x10
    const static inline Real                     mat_cmplx_1_10x10_det_abs = 1042200039.51399;
    const static inline MML::Vector<Real>        mat_cmplx_1_10x10_singular_values{22.7185465473606, 19.8152614906404, 17.208921749824, 14.7435183297016, 11.5502006454006, 10.2323923940024, 7.627380102872, 5.79022085629406, 2.68282848064847, 0.651607085870877};
    const static inline Real                     mat_cmplx_1_10x10_cond_2 = 34.8654074517885;
    const static inline int                      mat_cmplx_1_10x10_rank = 10;

    //////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Complex>     mat_cmplx_2_10x10{10, 10, {
            Complex(-1.1,-0.7), Complex(2.2,0.3), Complex(0.6,4.8), Complex(2.7,-2.0), Complex(4.4,-1.5), Complex(-2.7,4.5), Complex(-4.3,-3.4), Complex(-3.0,-3.8), Complex(2.9,-3.5), Complex(-3.4,-2.4),
            Complex(-0.7,3.5), Complex(-0.6,3.7), Complex(-0.9,3.7), Complex(-1.5,-0.4), Complex(3.0,1.0), Complex(5.0,2.2), Complex(4.5,-1.2), Complex(-0.9,3.6), Complex(3.2,0.0), Complex(0.9,0.9),
            Complex(4.3,1.6), Complex(-2.7,-1.1), Complex(-4.1,3.3), Complex(2.1,-3.7), Complex(-0.4,-2.2), Complex(1.0,-0.2), Complex(-4.5,-3.4), Complex(-4.0,3.8), Complex(-4.7,-0.6), Complex(3.4,-4.2),
            Complex(-2.8,2.8), Complex(-4.1,-4.1), Complex(2.0,0.7), Complex(4.4,3.0), Complex(0.5,1.9), Complex(2.8,-4.7), Complex(-4.0,0.2), Complex(-0.8,-0.2), Complex(3.7,-2.6), Complex(2.0,-3.2),
            Complex(0.1,1.8), Complex(4.7,-2.7), Complex(-4.1,4.1), Complex(-0.3,2.4), Complex(-4.2,-3.7), Complex(-2.6,-1.8), Complex(1.6,2.9), Complex(-2.6,-2.3), Complex(-1.7,3.6), Complex(4.9,-0.9),
            Complex(4.4,0.3), Complex(4.7,1.2), Complex(-1.4,4.3), Complex(-2.5,4.1), Complex(3.7,4.2), Complex(4.9,0.7), Complex(1.9,0.0), Complex(-4.3,-2.7), Complex(4.8,2.0), Complex(1.9,-4.3),
            Complex(2.0,4.0), Complex(-0.9,-3.7), Complex(-0.2,-0.9), Complex(4.7,0.8), Complex(1.3,-4.8), Complex(-4.6,4.7), Complex(0.4,-2.7), Complex(-0.8,-0.5), Complex(-0.5,1.6), Complex(-2.6,-3.6),
            Complex(0.9,-3.2), Complex(-1.5,-3.8), Complex(3.2,1.7), Complex(4.9,-4.8), Complex(-0.1,-2.2), Complex(4.0,1.4), Complex(-2.8,3.8), Complex(1.8,1.9), Complex(-1.3,0.4), Complex(2.6,-0.9),
            Complex(-2.8,3.5), Complex(0.8,-2.2), Complex(0.3,-4.5), Complex(3.1,-2.5), Complex(-4.1,4.5), Complex(4.2,-2.8), Complex(2.2,-2.8), Complex(2.9,-3.4), Complex(-3.2,0.1), Complex(-4.8,4.5),
            Complex(-2.6,3.9), Complex(1.4,0.7), Complex(2.3,3.8), Complex(0.5,-3.5), Complex(-3.1,-4.7), Complex(3.9,-3.3), Complex(4.1,-4.0), Complex(-0.2,2.1), Complex(-0.8,-1.3), Complex(-3.7,0.8)
        }};

    const static inline MML::Vector<Complex>     mat_cmplx_2_10x10_rhs0{
            Complex(1.3,-2.2),
            Complex(-1.5,2.5),
            Complex(3.0,0.0),
            Complex(-0.6,-0.6),
            Complex(1.3,1.1),
            Complex(0.7,-0.2),
            Complex(-2.9,0.5),
            Complex(-2.5,2.6),
            Complex(0.6,0.3),
            Complex(-2.7,-0.7)
        };
    const static inline MML::Vector<Complex>     mat_cmplx_2_10x10_rhs0_sol{
            Complex(0.706007116481964,0.013421144138452),
            Complex(-0.668876690042133,1.976565088318023),
            Complex(-1.550612209444719,-1.287669181025966),
            Complex(-1.704881833429372,-0.054239536509271),
            Complex(0.105269869191405,-0.404764634517515),
            Complex(-0.009580978242681,-0.232977323759457),
            Complex(0.500618952349040,-0.372903391127123),
            Complex(1.651960217269608,0.563573074178541),
            Complex(-0.457765142442384,2.177776204331938),
            Complex(-0.596210857708472,0.239850315657696)
        };

    // Matrix properties for mat_cmplx_2_10x10
    const static inline Real                     mat_cmplx_2_10x10_det_abs = 3907151884.60909;
    const static inline MML::Vector<Real>        mat_cmplx_2_10x10_singular_values{22.3486893729876, 18.4057371396619, 16.0143867894748, 14.9306536661108, 13.7285234336414, 9.71270306087155, 8.32484844148929, 7.14073463579936, 6.0902048638322, 0.82290710430045};
    const static inline Real                     mat_cmplx_2_10x10_cond_2 = 27.1582165911499;
    const static inline int                      mat_cmplx_2_10x10_rank = 10;

    //////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Complex>     mat_cmplx_3_10x10{10, 10, {
            Complex(-2.5,1.0), Complex(3.8,-3.9), Complex(4.4,-1.8), Complex(-1.9,-2.3), Complex(-0.3,-1.7), Complex(-3.3,-3.9), Complex(-2.8,0.9), Complex(-2.2,4.9), Complex(1.2,1.4), Complex(-3.0,4.7),
            Complex(3.6,-1.3), Complex(-1.0,1.9), Complex(2.0,3.1), Complex(-0.3,-4.8), Complex(-3.9,0.7), Complex(0.2,-3.5), Complex(3.9,-0.6), Complex(-2.1,-2.4), Complex(3.1,-1.9), Complex(-4.8,1.5),
            Complex(2.7,-3.6), Complex(1.4,3.7), Complex(1.0,0.3), Complex(-0.1,-3.4), Complex(2.6,-2.8), Complex(-2.6,-4.5), Complex(4.7,-4.9), Complex(2.6,0.3), Complex(-4.5,-1.6), Complex(-4.3,0.4),
            Complex(-1.7,-3.7), Complex(-1.8,5.0), Complex(-3.2,3.6), Complex(-0.4,-2.8), Complex(2.3,3.0), Complex(4.2,0.4), Complex(-4.9,-3.6), Complex(-4.8,-0.2), Complex(-2.5,-3.7), Complex(-0.2,3.7),
            Complex(-4.0,2.5), Complex(0.9,-0.3), Complex(-2.2,-0.5), Complex(1.5,3.2), Complex(-4.0,4.5), Complex(1.9,0.7), Complex(-4.4,3.5), Complex(2.7,-0.4), Complex(-0.5,3.5), Complex(-1.8,1.3),
            Complex(1.7,0.7), Complex(2.1,0.6), Complex(0.1,-3.8), Complex(-3.0,-3.9), Complex(4.5,-2.9), Complex(-1.4,0.5), Complex(0.9,-3.6), Complex(4.7,-1.5), Complex(2.8,0.2), Complex(5.0,2.4),
            Complex(-3.7,3.9), Complex(-0.8,-1.6), Complex(3.3,0.6), Complex(-0.7,-2.3), Complex(3.3,4.1), Complex(1.0,-1.0), Complex(0.5,-2.2), Complex(1.0,2.5), Complex(2.3,4.0), Complex(1.1,-4.4),
            Complex(0.5,4.4), Complex(4.8,-3.8), Complex(2.0,-3.5), Complex(-2.5,-2.4), Complex(1.7,1.1), Complex(-3.6,4.0), Complex(1.4,2.5), Complex(2.4,-0.9), Complex(-0.9,-4.9), Complex(-3.6,4.6),
            Complex(-2.6,2.3), Complex(4.4,-1.5), Complex(4.4,-4.5), Complex(-0.7,3.6), Complex(-4.7,1.0), Complex(0.7,-2.7), Complex(-4.6,5.0), Complex(0.9,-4.3), Complex(3.9,-4.1), Complex(-0.1,-4.8),
            Complex(4.4,3.1), Complex(-0.7,1.6), Complex(-4.4,1.6), Complex(3.4,-3.3), Complex(2.6,-1.6), Complex(0.3,3.8), Complex(-3.1,1.8), Complex(4.9,-2.5), Complex(-2.9,-0.5), Complex(0.0,1.5)
        }};

    const static inline MML::Vector<Complex>     mat_cmplx_3_10x10_rhs0{
            Complex(-0.0,0.6),
            Complex(-1.9,1.8),
            Complex(2.2,1.8),
            Complex(-2.9,0.4),
            Complex(1.5,-1.8),
            Complex(2.6,2.1),
            Complex(0.8,2.6),
            Complex(1.0,-0.5),
            Complex(-1.3,-0.2),
            Complex(-2.4,-2.3)
        };
    const static inline MML::Vector<Complex>     mat_cmplx_3_10x10_rhs0_sol{
            Complex(0.205631403861014,-0.069120662300960),
            Complex(0.607552380277906,-0.133279742702262),
            Complex(-0.627691172799875,0.093459414320121),
            Complex(-0.588592081392438,-0.390493395743841),
            Complex(0.051256231959498,-0.186957099902922),
            Complex(-0.059106269104728,0.076716221074585),
            Complex(0.226099826745400,0.256171743798648),
            Complex(0.167444522224238,-0.330821725681724),
            Complex(0.211147906212214,-0.251467714908403),
            Complex(-0.216706429129029,0.165134497920343)
        };

    // Matrix properties for mat_cmplx_3_10x10
    const static inline Real                     mat_cmplx_3_10x10_det_abs = 2412748639.16755;
    const static inline MML::Vector<Real>        mat_cmplx_3_10x10_singular_values{22.5163950951509, 20.275670042204, 17.0822720041155, 13.923744195654, 12.0162310557434, 9.04808201323421, 6.67171858689903, 4.89006586995963, 3.45262581668275, 1.81430205907735};
    const static inline Real                     mat_cmplx_3_10x10_cond_2 = 12.410499664318;
    const static inline int                      mat_cmplx_3_10x10_rank = 10;
}
#endif