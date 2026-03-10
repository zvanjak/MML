#if !defined MML_LINEAR_ALG_EQ_SYSTEMS_SINGULAR_DEFS_H
#define MML_LINEAR_ALG_EQ_SYSTEMS_SINGULAR_DEFS_H

// Singular and near-singular test matrices
// NOTE: All Matrix/Vector objects use construct-on-first-use idiom to avoid
// Static Initialization Order Fiasco (SIOF) across translation units.

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Matrix/Matrix.h"
#include "base/Matrix/MatrixSym.h"
#endif

namespace MML::TestBeds
{
    // 5x5 Singular matrix (rows 2 and 3 are identical)
    inline const Matrix<Real>& mat_5x5_singular1() {
        static const Matrix<Real> m{5, 5, { 
            1.4, 2.1, 2.1, 7.4, 9.6,
            1.6, 1.5, 1.1, 0.7, 5.0,
            1.6, 1.5, 1.1, 0.7, 5.0,
            4.6, 8.2, 8.4, 0.4, 8.0,
            2.6, 2.9, 0.1, 9.6, 7.7 
        }};
        return m;
    }

    // 5x5 Almost singular matrix (row 3 differs from row 2 by 1e-6 in last element)
    inline const Matrix<Real>& mat_5x5_almost_singular1() {
        static const Matrix<Real> m{5, 5, { 
            1.4, 2.1, 2.1, 7.4, 9.6,
            1.6, 1.5, 1.1, 0.7, 5.0,
            1.6, 1.5, 1.1, 0.7, 5.000001,
            4.6, 8.2, 8.4, 0.4, 8.0,
            2.6, 2.9, 0.1, 9.6, 7.7 
        }};
        return m;
    }

}

#endif 