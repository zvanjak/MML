///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DynamicalSystemAnalyzer.h                                           ///
///  Description: Unified facade for dynamical systems analysis                       ///
///               One class to access ALL dynamical system analysis capabilities      ///
///               Fixed points, Lyapunov exponents, bifurcations, phase space         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DYNAMICAL_SYSTEM_ANALYZER_H
#define MML_DYNAMICAL_SYSTEM_ANALYZER_H

#include "MMLBase.h"
#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"
#include "interfaces/IDynamicalSystem.h"
#include "systems/DynamicalSystemTypes.h"
#include "systems/DynamicalSystemAnalyzers.h"

#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <optional>
#include <cmath>

namespace MML::Systems 
{
	//=============================================================================
	// COMPREHENSIVE ANALYSIS RESULT
	//=============================================================================

	/// @brief Complete analysis report for a dynamical system
	template<typename Type = Real>
	struct DynamicalSystemReport {
		// System info
		int dimension;
		int numParameters;
		std::vector<std::string> stateNames;
		std::vector<std::string> paramNames;
		Vector<Type> currentParams;

		// Fixed points
		std::vector<FixedPoint<Type>> fixedPoints;
		int numStableFixedPoints;
		int numUnstableFixedPoints;
		int numSaddlePoints;

		// Lyapunov analysis
		LyapunovResult<Type> lyapunov;
		bool isChaotic;
		bool isDissipative;	  ///< Sum of Lyapunov exponents < 0
		bool isConservative;  ///< Sum ≈ 0 (Hamiltonian-like)

		// Attractor properties
		Type kaplanYorkeDimension;
		Type lyapunovTime;  ///< 1/λ_max - predictability time scale

		// Trajectory info (if computed)
		std::vector<Vector<Type>> trajectory;
		std::vector<Vector<Type>> poincareSection;

		// Text report
		std::string summary;
	};

	//=============================================================================
	// DYNAMICAL SYSTEM ANALYZER
	//=============================================================================

	/// @class DynamicalSystemAnalyzer
	///
	/// @brief Unified facade for all MML dynamical systems analysis capabilities
	///
	/// DynamicalSystemAnalyzer provides a single, comprehensive interface to:
	/// - Find and classify fixed points
	/// - Compute Lyapunov exponents and assess chaos
	/// - Generate bifurcation diagrams
	/// - Compute Poincaré sections
	/// - Integrate trajectories
	/// - Generate comprehensive analysis reports
	///
	/// @par Example: Basic Analysis
	/// @code
	/// LorenzSystem lorenz(10, 28, 8.0/3.0);
	/// DynamicalSystemAnalyzer analyzer(lorenz);
	///
	/// // Quick chaos check
	/// bool chaotic = analyzer.IsChaotic({1, 1, 1});
	///
	/// // Full analysis
	/// auto report = analyzer.Analyze({1, 1, 1});
	/// std::cout << report.summary;
	/// @endcode
	///
	/// @par Example: Specific Analysis
	/// @code
	/// DynamicalSystemAnalyzer analyzer(rossler);
	///
	/// // Find fixed points
	/// auto fps = analyzer.FindFixedPoints({{0,0,0}, {1,1,1}, {-1,-1,1}});
	///
	/// // Lyapunov exponents only
	/// auto lyap = analyzer.ComputeLyapunov({0.1, 0.1, 0.1}, 1000);
	///
	/// // Bifurcation diagram
	/// auto bif = analyzer.ComputeBifurcation(0, 20.0, 30.0, 100, {0.1, 0.1, 0.1});
	/// @endcode
	template<typename Type = Real>
	class DynamicalSystemAnalyzer {
	public:
		//=========================================================================
		// CONSTRUCTOR
		//=========================================================================

		/// @brief Construct analyzer for a dynamical system
		/// @param sys Reference to the dynamical system to analyze
		explicit DynamicalSystemAnalyzer(IDynamicalSystem& sys)
				: _sys(sys) {}

		//=========================================================================
		// QUICK QUERIES
		//=========================================================================

		/// @brief Check if system exhibits chaotic behavior from given initial condition
		/// @param x0 Initial condition
		/// @param integrationTime Time to integrate (longer = more accurate)
		/// @return True if largest Lyapunov exponent > 0
		bool IsChaotic(const Vector<Type>& x0, Type integrationTime = 500) const {
			auto result = LyapunovAnalyzer::Compute(_sys, x0, integrationTime);
			return result.isChaotic;
		}

		/// @brief Check if system is dissipative (contracting phase space volume)
		/// @param x0 Initial condition for Lyapunov computation
		/// @return True if sum of Lyapunov exponents < 0
		bool IsDissipative(const Vector<Type>& x0, Type integrationTime = 500) const {
			auto result = LyapunovAnalyzer::Compute(_sys, x0, integrationTime);
			return result.sum < -1e-6;
		}

		/// @brief Get the Lyapunov time (predictability horizon)
		/// @param x0 Initial condition
		/// @return Time scale 1/λ_max (smaller = less predictable)
		Type GetLyapunovTime(const Vector<Type>& x0, Type integrationTime = 500) const {
			auto result = LyapunovAnalyzer::Compute(_sys, x0, integrationTime);
			if (result.maxExponent > 1e-10)
				return 1.0 / result.maxExponent;
			else
				return std::numeric_limits<Type>::infinity();
		}

		/// @brief Get the fractal (Kaplan-Yorke) dimension of the attractor
		Type GetFractalDimension(const Vector<Type>& x0, Type integrationTime = 500) const {
			auto result = LyapunovAnalyzer::Compute(_sys, x0, integrationTime);
			return result.kaplanYorkeDimension;
		}

		//=========================================================================
		// FIXED POINT ANALYSIS
		//=========================================================================

		/// @brief Find a single fixed point near initial guess
		/// @param initialGuess Starting point for Newton iteration
		/// @param tol Convergence tolerance
		/// @return Fixed point result (check convergenceResidual for success)
		FixedPoint<Type> FindFixedPoint(const Vector<Type>& initialGuess, Type tol = 1e-10) const {
			return FixedPointFinder::Find(_sys, initialGuess, tol);
		}

		/// @brief Find multiple fixed points from several initial guesses
		/// @param initialGuesses Vector of starting points
		/// @param tol Convergence tolerance
		/// @param uniqueTol Distance below which two points are considered same
		/// @return Vector of unique fixed points found
		std::vector<FixedPoint<Type>> FindFixedPoints(const std::vector<Vector<Type>>& initialGuesses,
																									Type tol = 1e-10, Type uniqueTol = 1e-6) const {
			return FixedPointFinder::FindMultiple(_sys, initialGuesses, tol, uniqueTol);
		}

		/// @brief Generate a grid of initial guesses and find all fixed points
		/// @param boxMin Lower corner of search box
		/// @param boxMax Upper corner of search box
		/// @param gridPoints Number of grid points per dimension
		/// @return Vector of unique fixed points found
		std::vector<FixedPoint<Type>> FindFixedPointsInBox(const Vector<Type>& boxMin, const Vector<Type>& boxMax,
																											 int gridPoints = 5) const {
			int n = _sys.getDim();
			std::vector<Vector<Type>> guesses;

			// Generate grid
			std::vector<int> indices(n, 0);
			bool done = false;

			while (!done) {
				Vector<Type> guess(n);
				for (int i = 0; i < n; ++i) {
					Type t = (gridPoints > 1) ? static_cast<Type>(indices[i]) / (gridPoints - 1) : 0.5;
					guess[i] = boxMin[i] + t * (boxMax[i] - boxMin[i]);
				}
				guesses.push_back(guess);

				// Increment indices (like counting in base gridPoints)
				int carry = 1;
				for (int i = 0; i < n && carry; ++i) {
					indices[i] += carry;
					if (indices[i] >= gridPoints) {
						indices[i] = 0;
						carry = 1;
					} else {
						carry = 0;
					}
				}
				if (carry) done = true;
			}

			return FindFixedPoints(guesses);
		}

		//=========================================================================
		// LYAPUNOV ANALYSIS
		//=========================================================================

		/// @brief Compute all Lyapunov exponents
		/// @param x0 Initial condition
		/// @param totalTime Total integration time
		/// @param orthTime Time between orthonormalizations
		/// @param stepSize Integration step size
		/// @return LyapunovResult with exponents and derived quantities
		LyapunovResult<Type> ComputeLyapunov(const Vector<Type>& x0, Type totalTime = 1000,
																				 Type orthTime = 1.0, Type stepSize = 0.01) const {
			return LyapunovAnalyzer::Compute(_sys, x0, totalTime, orthTime, stepSize);
		}

		/// @brief Compute just the largest Lyapunov exponent (faster)
		/// @param x0 Initial condition
		/// @param totalTime Integration time
		/// @return Largest Lyapunov exponent λ_1
		Type ComputeMaxLyapunov(const Vector<Type>& x0, Type totalTime = 1000) const {
			auto result = LyapunovAnalyzer::Compute(_sys, x0, totalTime);
			return result.maxExponent;
		}

		//=========================================================================
		// BIFURCATION ANALYSIS
		//=========================================================================

		/// @brief Compute bifurcation diagram by sweeping a parameter
		/// @param paramIndex Which parameter to sweep
		/// @param paramMin Starting parameter value
		/// @param paramMax Ending parameter value
		/// @param numSteps Number of parameter values to sample
		/// @param x0 Initial condition
		/// @param component Which state component to record (default: 0)
		/// @param tTransient Transient integration time
		/// @param tRecord Recording time after transient
		/// @return BifurcationDiagram with parameter values and attractor points
		BifurcationDiagram<Type> ComputeBifurcation(int paramIndex, Type paramMin, Type paramMax, int numSteps,
																								const Vector<Type>& x0, int component = 0,
																								Type tTransient = 100, Type tRecord = 50) const {
			return BifurcationAnalyzer::Sweep(_sys, paramIndex, paramMin, paramMax, numSteps,
																				x0, component, tTransient, tRecord);
		}

		//=========================================================================
		// PHASE SPACE ANALYSIS
		//=========================================================================

		/// @brief Compute Poincaré section
		/// @param x0 Initial condition
		/// @param section Definition of the Poincaré section
		/// @param numIntersections Number of intersections to collect
		/// @param stepSize Integration step size
		/// @return Vector of intersection points
		std::vector<Vector<Type>> ComputePoincareSection(const Vector<Type>& x0, const PoincareSection<Type>& section,
																										 int numIntersections = 1000, Type stepSize = 0.01) const {
			return PhaseSpaceAnalyzer::ComputePoincareSection(_sys, x0, section, numIntersections, stepSize);
		}

		/// @brief Integrate trajectory from initial condition
		/// @param x0 Initial condition
		/// @param totalTime Total integration time
		/// @param outputInterval Time between output points
		/// @param stepSize Integration step size
		/// @return Vector of state points along trajectory
		std::vector<Vector<Type>> IntegrateTrajectory(const Vector<Type>& x0, Type totalTime,
																									Type outputInterval = 0.1, Type stepSize = 0.01) const {
			return PhaseSpaceAnalyzer::IntegrateTrajectory(_sys, x0, totalTime, outputInterval, stepSize);
		}

		//=========================================================================
		// COMPREHENSIVE ANALYSIS
		//=========================================================================

		/// @brief Perform comprehensive analysis of the system
		/// @param x0 Initial condition for trajectory-based analyses
		/// @param fixedPointGuesses Initial guesses for fixed point search (optional)
		/// @param lyapunovTime Integration time for Lyapunov computation
		/// @return Complete analysis report
		DynamicalSystemReport<Type> Analyze(const Vector<Type>& x0,
																				const std::vector<Vector<Type>>& fixedPointGuesses = {},
																				Type lyapunovTime = 500) const {
			DynamicalSystemReport<Type> report;
			int n = _sys.getDim();
			int p = _sys.getNumParam();

			// Basic info
			report.dimension = n;
			report.numParameters = p;
			report.stateNames.resize(n);
			report.paramNames.resize(p);
			report.currentParams.Resize(p);

			for (int i = 0; i < n; ++i)
				report.stateNames[i] = _sys.getStateName(i);
			for (int i = 0; i < p; ++i) {
				report.paramNames[i] = _sys.getParamName(i);
				report.currentParams[i] = _sys.getParam(i);
			}

			// Fixed points
			if (!fixedPointGuesses.empty()) {
				report.fixedPoints = FindFixedPoints(fixedPointGuesses);
			} else {
				// Use origin and a few random-ish points
				std::vector<Vector<Type>> defaultGuesses;
				Vector<Type> origin(n, 0.0);
				defaultGuesses.push_back(origin);

				// Add some points near x0
				defaultGuesses.push_back(x0);
				Vector<Type> scaled = x0 * 0.5;
				defaultGuesses.push_back(scaled);
				scaled = x0 * (-1.0);
				defaultGuesses.push_back(scaled);

				report.fixedPoints = FindFixedPoints(defaultGuesses);
			}

			// Count fixed point types
			report.numStableFixedPoints = 0;
			report.numUnstableFixedPoints = 0;
			report.numSaddlePoints = 0;
			for (const auto& fp : report.fixedPoints) {
				if (fp.isStable)
					report.numStableFixedPoints++;
				else if (fp.type == FixedPointType::Saddle)
					report.numSaddlePoints++;
				else
					report.numUnstableFixedPoints++;
			}

			// Lyapunov analysis
			report.lyapunov = ComputeLyapunov(x0, lyapunovTime);
			report.isChaotic = report.lyapunov.isChaotic;
			report.isDissipative = report.lyapunov.sum < -1e-6;
			report.isConservative = std::abs(report.lyapunov.sum) < 1e-4;
			report.kaplanYorkeDimension = report.lyapunov.kaplanYorkeDimension;

			if (report.lyapunov.maxExponent > 1e-10)
				report.lyapunovTime = 1.0 / report.lyapunov.maxExponent;
			else
				report.lyapunovTime = std::numeric_limits<Type>::infinity();

			// Generate summary
			report.summary = GenerateSummary(report);

			return report;
		}

		//=========================================================================
		// PARAMETER SENSITIVITY
		//=========================================================================

		/// @brief Check how chaos depends on a parameter
		/// @param paramIndex Which parameter to vary
		/// @param paramMin Starting value
		/// @param paramMax Ending value
		/// @param numSteps Number of parameter values
		/// @param x0 Initial condition
		/// @return Vector of pairs (param_value, max_lyapunov)
		std::vector<std::pair<Type, Type>> ScanChaosVsParameter(int paramIndex, Type paramMin, Type paramMax,
																														int numSteps, const Vector<Type>& x0,
																														Type lyapunovTime = 200) const {
			std::vector<std::pair<Type, Type>> results;
			Type dParam = (paramMax - paramMin) / (numSteps - 1);
			Type originalParam = _sys.getParam(paramIndex);

			for (int i = 0; i < numSteps; ++i) {
				Type param = paramMin + i * dParam;
				_sys.setParam(paramIndex, param);

				auto lyap = ComputeLyapunov(x0, lyapunovTime);
				results.push_back({param, lyap.maxExponent});
			}

			// Restore original parameter
			_sys.setParam(paramIndex, originalParam);
			return results;
		}

	private:
		IDynamicalSystem& _sys;

		/// @brief Generate human-readable summary
		std::string GenerateSummary(const DynamicalSystemReport<Type>& report) const {
			std::ostringstream ss;
			ss << std::fixed << std::setprecision(4);

			ss << "═══════════════════════════════════════════════════════════════\n";
			ss << "                DYNAMICAL SYSTEM ANALYSIS REPORT                \n";
			ss << "═══════════════════════════════════════════════════════════════\n\n";

			// System overview
			ss << "SYSTEM OVERVIEW\n";
			ss << "───────────────────────────────────────────────────────────────\n";
			ss << "  Dimension: " << report.dimension << "\n";
			ss << "  State variables: ";
			for (int i = 0; i < report.dimension; ++i) {
				if (i > 0) ss << ", ";
				ss << report.stateNames[i];
			}
			ss << "\n";

			ss << "  Parameters (" << report.numParameters << "): ";
			for (int i = 0; i < report.numParameters; ++i) {
				if (i > 0) ss << ", ";
				ss << report.paramNames[i] << "=" << report.currentParams[i];
			}
			ss << "\n\n";

			// Fixed points
			ss << "FIXED POINTS\n";
			ss << "───────────────────────────────────────────────────────────────\n";
			ss << "  Found: " << report.fixedPoints.size() << " fixed point(s)\n";
			ss << "    Stable:   " << report.numStableFixedPoints << "\n";
			ss << "    Saddles:  " << report.numSaddlePoints << "\n";
			ss << "    Unstable: " << report.numUnstableFixedPoints << "\n";

			for (size_t i = 0; i < report.fixedPoints.size(); ++i) {
				const auto& fp = report.fixedPoints[i];
				ss << "\n  Fixed Point #" << (i+1) << ": " << ToString(fp.type) << "\n";
				ss << "    Location: (";
				for (int j = 0; j < fp.location.size(); ++j) {
					if (j > 0) ss << ", ";
					ss << fp.location[j];
				}
				ss << ")\n";
				ss << "    Eigenvalues: ";
				for (size_t j = 0; j < fp.eigenvalues.size(); ++j) {
					if (j > 0) ss << ", ";
					ss << fp.eigenvalues[j].real();
					if (std::abs(fp.eigenvalues[j].imag()) > 1e-10) {
						ss << (fp.eigenvalues[j].imag() > 0 ? "+" : "") << fp.eigenvalues[j].imag() << "i";
					}
				}
				ss << "\n";
			}
			ss << "\n";

			// Lyapunov analysis
			ss << "LYAPUNOV ANALYSIS\n";
			ss << "───────────────────────────────────────────────────────────────\n";
			ss << "  Lyapunov exponents: (";
			for (int i = 0; i < report.lyapunov.exponents.size(); ++i) {
				if (i > 0) ss << ", ";
				ss << report.lyapunov.exponents[i];
			}
			ss << ")\n";
			ss << "  Maximum exponent (λ₁): " << report.lyapunov.maxExponent << "\n";
			ss << "  Sum of exponents: " << report.lyapunov.sum << "\n";
			ss << "  Kaplan-Yorke dimension: " << report.kaplanYorkeDimension << "\n";
			if (report.lyapunovTime < 1e10)
				ss << "  Lyapunov time (1/λ₁): " << report.lyapunovTime << "\n";
			ss << "\n";

			// Classification
			ss << "DYNAMICAL CLASSIFICATION\n";
			ss << "───────────────────────────────────────────────────────────────\n";
			if (report.isChaotic) {
				ss << "  ⚡ CHAOTIC DYNAMICS (λ₁ > 0)\n";
				ss << "     Sensitive dependence on initial conditions.\n";
				ss << "     Predictability horizon ~ " << report.lyapunovTime << " time units.\n";
			} else {
				ss << "  ✓ REGULAR DYNAMICS (λ₁ ≤ 0)\n";
				ss << "     Trajectories converge or remain bounded.\n";
			}

			if (report.isDissipative) {
				ss << "  📉 DISSIPATIVE (Σλ < 0)\n";
				ss << "     Phase space volume contracts.\n";
				ss << "     Attractor has dimension < state space.\n";
			} else if (report.isConservative) {
				ss << "  ⚖ CONSERVATIVE (Σλ ≈ 0)\n";
				ss << "     Phase space volume preserved (Hamiltonian-like).\n";
			}

			ss << "\n═══════════════════════════════════════════════════════════════\n";

			return ss.str();
		}
	};

} // namespace MML::Systems
#endif // MML_DYNAMICAL_SYSTEM_ANALYZER_H
