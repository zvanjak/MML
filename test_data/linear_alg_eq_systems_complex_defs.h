#if !defined MML_LINEAR_ALG_EQ_SYSTEMS_COMPLEX_DEFS_H
#define MML_LINEAR_ALG_EQ_SYSTEMS_COMPLEX_DEFS_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Matrix.h"
#include "core/MatrixSym.h"
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

    const static inline VectorComplex   mat_cmplx_3x3_eigen_val{ Complex(6.6894939101402411,  0), 
                                                           Complex(0.15525304492987929, 2.8104746742964433), 
                                                           Complex(0.15525304492987929, -2.8104746742964433) };
    const static inline std::vector<VectorComplex> mat_cmplx_3x3_eigen_vecs
    {
        VectorComplex{ Complex( 0.2559569335, 0),            Complex( 0.8722704144, 0),           Complex(0.2882754144, 0) },
        VectorComplex{ Complex(-0.6589654354, 0.4832617401), Complex(-0.189800524,-0.8080623097), Complex(0.4219347887, 0.6441149315) },
        VectorComplex{ Complex(-0.6589654354,-0.4832617401), Complex(-0.189800524, 0.8080623097), Complex(0.4219347887,-0.6441149315) }
    };

    /////////////////////////////////////////////////////////////////////////////////////////////////
    const static inline MML::Matrix<Complex>     mat_cmplx_1_3x3{3, 3, { Complex(1.0, 2.0), Complex(2.0, -0.5), Complex(-1.0, 2.3), 
                                                                         Complex(-2.0,1.0), Complex(4.0,-2.0),  Complex(6.0,5.0), 
                                                                         Complex(2.5,1.0),  Complex(1.0,-1.0),  Complex(-2.0, -4.0) }};

    const static inline MML::Vector<Complex>     mat_cmplx_1_3x3_rhs0{ Complex(1.5,0.5), Complex(-2.0,1.5), Complex(-1.0,-2.0)};
    const static inline MML::Vector<Complex>     mat_cmplx_1_3x3_rhs0_sol{  Complex( 0.69342378330657772, -1.2500064784903142),
                                                                            Complex(-1.0954914099610633,  -0.59693995647332976),
                                                                            Complex( 0.24250526514594178, -0.45153101000762030) 
                                                                         };

    const static inline VectorComplex   mat_cmplx_1_3x3_eigen_val{ Complex(6.6894939101402411,  0), 
                                                           Complex(0.15525304492987929, 2.8104746742964433), 
                                                           Complex(0.15525304492987929, -2.8104746742964433) };
    const static inline std::vector<VectorComplex> mat_cmplx_1_3x3_eigen_vecs
    {
        VectorComplex{ Complex( 0.2559569335, 0),            Complex( 0.8722704144, 0),           Complex(0.2882754144, 0) },
        VectorComplex{ Complex(-0.6589654354, 0.4832617401), Complex(-0.189800524,-0.8080623097), Complex(0.4219347887, 0.6441149315) },
        VectorComplex{ Complex(-0.6589654354,-0.4832617401), Complex(-0.189800524, 0.8080623097), Complex(0.4219347887,-0.6441149315) }
    };    
}
#endif