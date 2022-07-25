#if !defined __MML_LINEAR_ALG_EQ_SOLVERS_TEST_BED_H
#define __MML_LINEAR_ALG_EQ_SOLVERS_TEST_BED_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/Matrix.h"
#endif

namespace MML::Tests
{
    class LinearAlgEqTestBed
    {
    public:
        const static inline MML::Matrix     mat0{3, 3, {1.0, 2.0, -1.0, 
                                                        -1.0, 5.0, 6.0, 
                                                        3.0, 1.0, 1.0 }};

        const static inline MML::Vector     mat0_rhs0{1.0, 2.0, 1.0};
        const static inline MML::Vector     mat0_rhs0_sol{0.18867924528301885, 0.41509433962264153, 0.018867924528301921};

    };
}

#endif