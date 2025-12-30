///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        LinAlgEqSystem.h                                                    ///
///  Description: Linear algebraic equation system class                              ///
///               High-level interface for solving Ax=b systems                       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_LIN_ALG_EQ_SYSTEM_H
#define MML_LIN_ALG_EQ_SYSTEM_H

#include "MMLBase.h"

namespace MML
{
	class LinAlgEqSysSolver
	{
		// dva ctora - Matrix, Matrix i Matrix, Vector

		// solve by GJ

		// solve by LU - dvije verzije - in place, i verzija koja vraca dekompoziciju
		// perform LU
		//
		// QR, Cholesky

		// SVD decomposition

		// find eigenvalues

		// ima i Verify - das mu solution i onda vidis kolika je norma razlike u odnosu na zadani rhs

	};
}
#endif


