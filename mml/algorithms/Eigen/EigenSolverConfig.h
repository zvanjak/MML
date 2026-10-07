///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        EigenSolverConfig.h                                                 ///
///  Description: Shared configuration for eigenvalue/eigenvector solvers             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined(MML_EIGEN_SOLVER_CONFIG_H)
#define MML_EIGEN_SOLVER_CONFIG_H

#include <mml/MMLBase.h>

namespace MML {

	/// Configuration for eigenvalue/eigenvector solvers.
	///
	/// Provides user control over convergence criteria, iteration limits,
	/// and algorithm behavior. All fields have sensible defaults.
	struct EigenSolverConfig {
		/// Convergence tolerance for off-diagonal elements.
		Real tolerance = PrecisionValues<Real>::EigenSolverConvergenceTolerance;

		/// Maximum number of iterations/sweeps. Set to 0 for the solver default.
		int max_iterations = 100;

		/// Sort eigenvalues in ascending order.
		bool sort_eigenvalues = true;

		/// Compute eigenvectors when supported by the solver.
		bool compute_eigenvectors = true;

		/// Enable verbose diagnostic output.
		bool verbose = false;

		static EigenSolverConfig HighPrecision() {
			EigenSolverConfig config;
			config.tolerance = 1e-14;
			config.max_iterations = 500;
			return config;
		}

		static EigenSolverConfig Fast() {
			EigenSolverConfig config;
			config.tolerance = 1e-6;
			config.max_iterations = 50;
			return config;
		}
	};

} // namespace MML

#endif // MML_EIGEN_SOLVER_CONFIG_H