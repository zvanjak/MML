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

    /////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Complex>     mat_cmplx_1_3x3{3, 3, { Complex(1.0, 2.0), Complex(2.0, -0.5), Complex(-1.0, 2.3), 
                                                                         Complex(-2.0,1.0), Complex(4.0,-2.0),  Complex(6.0,5.0), 
                                                                         Complex(2.5,1.0),  Complex(1.0,-1.0),  Complex(-2.0, -4.0) }};

    const static inline MML::Vector<Complex>     mat_cmplx_1_3x3_rhs0{ Complex(1.5,0.5), Complex(-2.0,1.5), Complex(-1.0,-2.0)};
    const static inline MML::Vector<Complex>     mat_cmplx_1_3x3_rhs0_sol{ Complex(0.693423783306577,-1.25000647849031),
                                                                            Complex(-1.09549140996106,-0.596939956473329),
                                                                            Complex(0.242505265145942,-0.451531010007621)};

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

    // TODO - LOW, dodati 10x10 i 20x20 complex test
}
#endif