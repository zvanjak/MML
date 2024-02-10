#if !defined MML_LINEAR_ALG_EQ_SYSTEMS_SINGULAR_DEFS_H
#define MML_LINEAR_ALG_EQ_SYSTEMS_SINGULAR_DEFS_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Matrix.h"
#include "base/MatrixSym.h"
#endif

namespace MML::TestBeds
{

    const static inline MML::Matrix<Real>     mat_5x5_singular1{5, 5, { 1.4, 2.1, 2.1, 7.4, 9.6,
                                                                  1.6, 1.5, 1.1, 0.7, 5.0,
                                                                  1.6, 1.5, 1.1, 0.7, 5.0,
                                                                  4.6, 8.2, 8.4, 0.4, 8.0,
                                                                  2.6, 2.9, 0.1, 9.6, 7.7 } };

    const static inline MML::Matrix<Real>     mat_5x5_almost_singular1{5, 5, { 1.4, 2.1, 2.1, 7.4, 9.6,
                                                                         1.6, 1.5, 1.1, 0.7, 5.0,
                                                                         1.6, 1.5, 1.1, 0.7, 5.000001,
                                                                         4.6, 8.2, 8.4, 0.4, 8.0,
                                                                         2.6, 2.9, 0.1, 9.6, 7.7 } };

}

#endif 