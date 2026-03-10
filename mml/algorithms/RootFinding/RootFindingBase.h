///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        RootFindingBase.h                                                   ///
///  Description: Configuration and result types for root-finding algorithms         ///
///               Configuration structs, result types, and forward declarations       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ROOTFINDING_BASE_H
#define MML_ROOTFINDING_BASE_H

#include "mml/MMLBase.h"
#include "mml/core/AlgorithmTypes.h"

namespace MML {
	namespace RootFinding {

		/*********************************************************************/
		/*****           Configuration and Result Types                  *****/
		/*********************************************************************/

		/// Configuration parameters for root-finding algorithms.
		/// 
		/// Provides user control over convergence criteria, iteration limits,
		/// and diagnostic output. All fields have sensible defaults.
		/// 
		/// @example
		/// RootFindingConfig config;
		/// config.tolerance = 1e-12;        // Higher precision
		/// config.max_iterations = 200;     // Allow more iterations
		/// auto result = FindRootBrent(f, 0, 1, config);
		struct RootFindingConfig {
			/// Absolute accuracy tolerance for root (default: 1e-10)
			Real tolerance = 1e-10;
			
			/// Maximum number of iterations (default: 100)
			/// Set to 0 to use algorithm-specific default from Defaults namespace
			int max_iterations = 100;
			
			/// Initial step size for derivative-free methods (default: 0.01)
			/// Used by some algorithms for initial search or derivative approximation
			Real initial_step = 0.01;
			
			/// Enable verbose output for debugging (default: false)
			/// When true, prints iteration details to std::cout
			bool verbose = false;
		};

		/// Result of a root-finding operation.
		/// 
		/// Contains the root value along with diagnostic information about
		/// the convergence process. Always check `converged` before using `root`.
		/// 
		/// @example
		/// auto result = FindRootBrent(f, 0, 1, config);
		/// if (result.converged) {
		///     std::cout << "Root: " << result.root << " found in " 
		///               << result.iterations_used << " iterations\n";
		/// } else {
		///     std::cerr << "Failed to converge after " << result.iterations_used << " iterations\n";
		/// }
		struct RootFindingResult {
			/// The computed root value
			Real root = 0.0;
			
			/// Function value at root: f(root), should be near zero if converged
			Real function_value = 0.0;
			
			/// Number of iterations actually used
			int iterations_used = 0;
			
			/// True if algorithm converged within tolerance and iteration limits
			bool converged = false;
			
			/// Actual achieved tolerance (may be better than requested)
			Real achieved_tolerance = 0.0;
			
			/// Algorithm termination status (Phase 8 standardization)
			AlgorithmStatus status = AlgorithmStatus::Success;
			
			/// Error message if not converged (empty string if successful)
			std::string error_message;
			
			/// Name of the algorithm used (e.g., "Brent", "Bisection", "Newton")
			std::string algorithm_name;
			
			/// Elapsed wall-clock time in milliseconds
			double elapsed_time_ms = 0.0;
		};

	} // namespace RootFinding
} // namespace MML

#endif // MML_ROOTFINDING_BASE_H
