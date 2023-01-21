#if !defined __MML_LINEAR_ALG_EQ_SOLVERS_TEST_BED_H
#define __MML_LINEAR_ALG_EQ_SOLVERS_TEST_BED_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/Matrix.h"
#endif

namespace MML::Tests
{
    // known LinAlgSys
    // - matrica, rhs, solution, i eigen
    class LinearAlgEqTestBed
    {
    public:
        const static inline MML::Matrix     mat0{3, 3, {1.0, 2.0, -1.0, 
                                                        -1.0, 5.0, 6.0, 
                                                        3.0, 1.0, 1.0 }};

        const static inline MML::Vector     mat0_rhs0{1.0, 2.0, 1.0};
        const static inline MML::Vector     mat0_rhs0_sol{0.18867924528301885, 0.41509433962264153, 0.018867924528301921};

        const static inline MML::Matrix     mat4{5, 5, {1.4, 2.1, 2.1, 7.4, 9.6,
                                                        1.6, 1.5, 1.1, 0.7, 5.0,
                                                        3.8, 8.0, 9.6, 5.4, 8.8,
                                                        4.6, 8.2, 8.4, 0.4, 8.0,
                                                        2.6, 2.9, 0.1, 9.6, 7.7}};

        const static inline MML::Matrix     mat4_rhs4{5, 2, {1.1, 1.6, 
                                                            4.7, 9.1, 
                                                            0.1, 4.0, 
                                                            9.3, 8.4, 
                                                            0.4, 4.1}};  

        const static inline MML::Matrix     mat4_rhs0_sol{5, 2, {-3.9032710424808688,  15.643114174796667, 
                                                                  5.2353433160849479, -11.587503332831671, 
                                                                 -3.2920957702478550,   4.4111268480786325, 
                                                                 -1.7183300108528281,   0.21432757972725644, 
                                                                  1.5832710097423177,  -0.70999930382454612}}; 

    };
}

#endif