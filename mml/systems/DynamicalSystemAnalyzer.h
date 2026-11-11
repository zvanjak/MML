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

#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/interfaces/IDynamicalSystem.h>
#include <mml/systems/DynamicalSystem/DynamicalSystemTypes.h>
#include <mml/systems/DynamicalSystem/FixedPointAnalysis.h>
#include <mml/systems/DynamicalSystem/LyapunovAnalysis.h>
#include <mml/systems/DynamicalSystem/BifurcationAnalysis.h>
#include <mml/systems/DynamicalSystem/PhaseSpaceAnalysis.h>

#include <vector>
#include <string>
#include <type_traits>
#include <sstream>
#include <iomanip>
#include <optional>
#include <cmath>

namespace MML::Systems 
{
	//=============================================================================
	// COMPREHENSIVE ANALYSIS RESULT
	//=============================================================================

	/// @brief Local fixed-point searches and finite-time trajectory indicators, not a global classification
	template<typename Type = Real>
	struct DynamicalSystemReport {
		// System info
		int dimension = 0;
		int numParameters = 0;
		std::vector<std::string> stateNames;
		std::vector<std::string> paramNames;
		Vector<Type> currentParams;

		// Fixed points
		std::vector<FixedPoint<Type>> fixedPoints;
		int numStableFixedPoints = 0; ///< Strictly decaying numerical linearizations among the found points
		int numUnstableFixedPoints = 0;
		int numSaddlePoints = 0;

		// Lyapunov analysis
		LyapunovResult<Type> lyapunov;
		bool isChaotic = false; ///< Legacy positive-growth heuristic; not a proof of chaos
		bool isDissipative = false; ///< Legacy finite-time sum < -Precision::LooseTolerance indicator
		bool isConservative = false; ///< Legacy near-zero finite-time sum indicator; not a Hamiltonian proof

		// Attractor properties
		Type kaplanYorkeDimension = std::numeric_limits<Type>::quiet_NaN(); ///< Spectrum-based estimate, unavailable until analyzed
		Type lyapunovTime = std::numeric_limits<Type>::quiet_NaN(); ///< Reciprocal detected growth rate; not a certified predictability horizon

		// Trajectory info (if computed)
		std::vector<Vector<Type>> trajectory;
		std::vector<Vector<Type>> poincareSection;

		// Text report
		std::string summary;
		int numNeutralFixedPoints = 0; ///< Linear center labels, not a nonlinear stability guarantee
		int numUnclassifiedFixedPoints = 0;
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
	/// - Estimate finite-time Lyapunov exponents and tangent-growth indicators
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
	/// // Quick positive-growth heuristic (not a proof of chaos)
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
	/// @note The template parameter is a placeholder: every engine the facade delegates to
	///       (FixedPointFinder, LyapunovAnalyzer, BifurcationAnalyzer, PhaseSpaceAnalyzer) is
	///       implemented for `Real` only, so `Type` must be `Real`. It is kept so that existing
	///       spellings `DynamicalSystemAnalyzer<>` / `DynamicalSystemAnalyzer<Real>` continue to
	///       compile; a different `Type` is rejected at compile time instead of silently
	///       converting (analysis v4, SYSTEMS S5).
	template<typename Type = Real>
	class DynamicalSystemAnalyzer {
		static_assert(std::is_same_v<Type, Real>,
			"DynamicalSystemAnalyzer is implemented for Real only; the Type parameter is a placeholder");
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

		/// @brief Compute the legacy finite-time positive-growth heuristic, not a proof of chaos
		/// @param x0 Initial condition
		/// @param integrationTime Finite integration duration; convergence must be checked separately
		/// @return True if the largest estimated exponent exceeds Precision::LyapunovChaosThreshold
		bool IsChaotic(const Vector<Type>& x0, Type integrationTime = 500) const {
			auto result = LyapunovAnalyzer::Compute(_sys, x0, integrationTime);
			return result.isChaotic;
		}

		/// @brief Compute a finite-time exponent-sum contraction indicator along this trajectory
		/// @param x0 Initial condition for Lyapunov computation
		/// @return True if the estimated sum is below -Precision::LooseTolerance; not a global property
		bool IsDissipative(const Vector<Type>& x0, Type integrationTime = 500) const {
			auto result = LyapunovAnalyzer::Compute(_sys, x0, integrationTime);
			return result.sum < -Precision::LooseTolerance;
		}

		/// @brief Get the reciprocal detected positive-growth rate, not a certified predictability horizon
		/// @param x0 Initial condition
		/// @return 1/λ_max above the growth threshold; otherwise infinity means no reported growth time scale
		Type GetLyapunovTime(const Vector<Type>& x0, Type integrationTime = 500) const {
			auto result = LyapunovAnalyzer::Compute(_sys, x0, integrationTime);
			if (result.maxExponent > Precision::LyapunovChaosThreshold)
				return MML::Real(1.0) / result.maxExponent;
			else
				return std::numeric_limits<Type>::infinity();
		}

		/// @brief Get the spectrum-based Kaplan-Yorke estimate, not a measured attractor/fractal dimension
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
		/// @return Fixed point result (check status; false isStable does not imply instability)
		FixedPoint<Type> FindFixedPoint(const Vector<Type>& initialGuess, Type tol = PrecisionValues<Type>::FixedPointTolerance) const {
			return FixedPointFinder::Find(_sys, initialGuess, tol);
		}

		/// @brief Find multiple fixed points from several initial guesses
		/// @param initialGuesses Vector of starting points
		/// @param tol Convergence tolerance
		/// @param uniqueTol Distance below which two points are considered same
		/// @return Vector of unique fixed points found
		std::vector<FixedPoint<Type>> FindFixedPoints(const std::vector<Vector<Type>>& initialGuesses,
																									Type tol = PrecisionValues<Type>::FixedPointTolerance, Type uniqueTol = PrecisionValues<Type>::FixedPointUniquenessTolerance) const {
			return FixedPointFinder::FindMultiple(_sys, initialGuesses, tol, uniqueTol);
		}

		/// @brief Search for unique fixed points from a bounded grid of initial guesses
		/// @param boxMin Lower corner of search box
		/// @param boxMax Upper corner of search box
		/// @param gridPoints Number of grid points per dimension
		/// @param maxGridStarts Maximum number of starting guesses, checked before allocation
		/// A one-point grid uses the box midpoint; this search does not prove completeness.
		/// @return Vector of unique fixed points found
		std::vector<FixedPoint<Type>> FindFixedPointsInBox(const Vector<Type>& boxMin, const Vector<Type>& boxMax,
																											 int gridPoints = 5,
																											 int maxGridStarts = PrecisionValues<Type>::FixedPointMaxGridStarts) const {
			int n = _sys.getDim();
			DynSysDetail::CheckInitialState(_sys, boxMin);
			DynSysDetail::CheckInitialState(_sys, boxMax);
			if (!_sys.isAutonomous())
				throw ArgumentError("Fixed-point box searches require an autonomous system");
			if (gridPoints < 1 || maxGridStarts < 1)
				throw ArgumentError("Box grid size and starting-guess budget must be positive");
			int total = 1;
			for (int i = 0; i < n; ++i) {
				if (boxMin[i] > boxMax[i])
					throw ArgumentError("Box lower bounds must not exceed upper bounds");
				if (total > maxGridStarts / gridPoints)
					throw ResourceLimitError("Fixed-point box grid exceeds the starting-guess budget");
				total *= gridPoints;
			}
			std::vector<Vector<Type>> guesses;
			guesses.reserve(total);

			// Generate grid
			std::vector<int> indices(n, 0);
			bool done = false;

			while (!done) {
				Vector<Type> guess(n);
				for (int i = 0; i < n; ++i) {
					Type t = (gridPoints > 1) ? static_cast<Type>(indices[i]) / (gridPoints - 1) : Type(0.5);
					guess[i] = std::lerp(boxMin[i], boxMax[i], t);
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
		/// State and physical time continue across phases and samples, starting at t=0.
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
		/// @param stepSize Initial adaptive integration step size
		/// @param maxSearchTime Integration-time budget for the crossing search (throws ConvergenceError when exhausted)
		/// @param transientTime Time integrated and discarded before crossings are recorded
		/// @return Vector of intersection points
		std::vector<Vector<Type>> ComputePoincareSection(const Vector<Type>& x0, const PoincareSection<Type>& section,
																										 int numIntersections = 1000, Type stepSize = 0.01,
																										 Type maxSearchTime = DefaultPoincareMaxSearchTime,
																										 Type transientTime = DefaultTransientTime) const {
			return PhaseSpaceAnalyzer::ComputePoincareSection(_sys, x0, section, numIntersections, stepSize, maxSearchTime, transientTime);
		}

		/// @brief Integrate trajectory from initial condition
		/// @param x0 Initial condition
		/// @param totalTime Total integration time
		/// @param outputInterval Time between output points
		/// @param stepSize Initial adaptive integration step size
		/// @return Vector of state points along trajectory
		std::vector<Vector<Type>> IntegrateTrajectory(const Vector<Type>& x0, Type totalTime,
																									Type outputInterval = 0.1, Type stepSize = 0.01) const {
			return PhaseSpaceAnalyzer::IntegrateTrajectory(_sys, x0, totalTime, outputInterval, stepSize);
		}

		//=========================================================================
		// COMPREHENSIVE ANALYSIS
		//=========================================================================

		/// @brief Combine local fixed-point searches with finite-time Lyapunov indicators
		/// @param x0 Initial condition for trajectory-based analyses
		/// @param fixedPointGuesses Initial guesses for fixed point search (optional)
		/// @param lyapunovTime Integration time for Lyapunov computation
		/// @return Non-exhaustive fixed points and finite-time estimates; no boundedness or global stability proof
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
			report.currentParams.resize(p);

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
				else if (fp.type == FixedPointType::Saddle || fp.type == FixedPointType::SaddleFocus)
					report.numSaddlePoints++;
				else if (fp.type == FixedPointType::UnstableNode || fp.type == FixedPointType::UnstableFocus
					|| fp.type == FixedPointType::UnstableStarNode)
					report.numUnstableFixedPoints++;
				else if (fp.type == FixedPointType::Center)
					report.numNeutralFixedPoints++;
				else
					report.numUnclassifiedFixedPoints++;
			}

			// Lyapunov analysis
			report.lyapunov = ComputeLyapunov(x0, lyapunovTime);
			report.isChaotic = report.lyapunov.isChaotic;
			report.isDissipative = report.lyapunov.sum < -Precision::LooseTolerance;
			report.isConservative = IsNearlyZero(report.lyapunov.sum, Precision::CoarseTolerance);
			report.kaplanYorkeDimension = report.lyapunov.kaplanYorkeDimension;

			if (report.lyapunov.maxExponent > Precision::LyapunovChaosThreshold)
				report.lyapunovTime = MML::Real(1.0) / report.lyapunov.maxExponent;
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
		/// Each sample restarts from x0 (unlike bifurcation continuation).
		/// One sample evaluates paramMin; descending ranges are supported.
		/// @param paramIndex Which parameter to vary
		/// @param paramMin Starting value
		/// @param paramMax Ending value
		/// @param numSteps Number of parameter values
		/// @param x0 Initial condition
		/// @return Vector of pairs (param_value, max_lyapunov)
		std::vector<std::pair<Type, Type>> ScanChaosVsParameter(int paramIndex, Type paramMin, Type paramMax,
																														int numSteps, const Vector<Type>& x0,
																														Type lyapunovTime = 200) const {
			DynSysDetail::CheckVariationalDimension(_sys.getDim());
			DynSysDetail::CheckInitialState(_sys, x0);
			DynSysDetail::CheckParameterScan(_sys, paramIndex, paramMin, paramMax, numSteps);
			DynSysDetail::RequirePositive(lyapunovTime, "Chaos scan Lyapunov time");
			std::vector<std::pair<Type, Type>> results;
			results.reserve(numSteps);
			return DynSysDetail::WithRestoredParameter(_sys, paramIndex, [&] {
				for (int i = 0; i < numSteps; ++i) {
					const Type fraction = numSteps > 1 ? Type(i) / (numSteps - 1) : Type(0);
					const Type param = std::lerp(paramMin, paramMax, fraction);
					_sys.setParam(paramIndex, param);
					auto lyap = ComputeLyapunov(x0, lyapunovTime);
					results.push_back({param, lyap.maxExponent});
				}
				return std::move(results);
			});
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
			ss << "  Search is non-exhaustive; no points found does not prove none exist.\n";
			ss << "  Stability counts describe local numerical linearizations, not global stability.\n";
			ss << "    Stable:   " << report.numStableFixedPoints << "\n";
			ss << "    Saddles:  " << report.numSaddlePoints << "\n";
			ss << "    Unstable: " << report.numUnstableFixedPoints << "\n";
			ss << "    Neutral (linear): " << report.numNeutralFixedPoints << "\n";
			ss << "    Unclassified: " << report.numUnclassifiedFixedPoints << "\n";

			for (size_t i = 0; i < report.fixedPoints.size(); ++i) {
				const auto& fp = report.fixedPoints[i];
				ss << "\n  Fixed Point #" << (i+1) << ": " << to_string(fp.type) << "\n";
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
					if (std::abs(fp.eigenvalues[j].imag()) > Precision::EqualityTolerance) {
						ss << (fp.eigenvalues[j].imag() > 0 ? "+" : "") << fp.eigenvalues[j].imag() << "i";
					}
				}
				ss << "\n";
			}
			ss << "\n";

			// Lyapunov analysis
			ss << "LYAPUNOV ANALYSIS\n";
			ss << "───────────────────────────────────────────────────────────────\n";
			ss << "  Finite-time duration: " << report.lyapunov.totalTime << "\n";
			ss << "  Lyapunov exponents: (";
			for (int i = 0; i < report.lyapunov.exponents.size(); ++i) {
				if (i > 0) ss << ", ";
				ss << report.lyapunov.exponents[i];
			}
			ss << ")\n";
			ss << "  Maximum exponent (λ₁): " << report.lyapunov.maxExponent << "\n";
			ss << "  Sum of exponents: " << report.lyapunov.sum << "\n";
			ss << "  Kaplan-Yorke dimension estimate: " << report.kaplanYorkeDimension << "\n";
			ss << "    Spectrum-based; invariant-set and asymptotic-spectrum assumptions are not verified.\n";
			ss << "    An attractor or measured fractal dimension is not established.\n";
			if (std::isfinite(report.lyapunovTime))
				ss << "  Reciprocal positive-growth rate (1/λ₁): " << report.lyapunovTime << "\n";
			else
				ss << "  No positive-growth time scale reported (threshold not exceeded).\n";
			ss << "    This is not a certified predictability horizon.\n";
			ss << "\n";

			// Classification
			ss << "FINITE-TIME INDICATORS\n";
			ss << "───────────────────────────────────────────────────────────────\n";
			ss << "  Finite-time positive growth: " << (report.isChaotic ? "detected" : "not detected")
				<< " (threshold " << Precision::LyapunovChaosThreshold << ")\n";
			ss << "    Tangent expansion is not a proof of chaos or a bounded chaotic attractor.\n";
			ss << "    No boundedness or convergence conclusion follows from these indicators.\n";
			ss << "  Volume contraction indicator: " << (report.isDissipative ? "detected" : "not detected")
				<< " (sum < -" << Precision::LooseTolerance << ")\n";
			ss << "    This concerns average volume growth along the finite-time trajectory, not global phase space.\n";
			ss << "  Near-zero-sum indicator: " << (report.isConservative ? "detected" : "not detected")
				<< " (absolute sum < " << Precision::CoarseTolerance << ")\n";
			ss << "    A near-zero sum is not proof of Hamiltonian structure or exact volume preservation.\n";
			ss << "    The two volume indicators use independent tolerances and can overlap.\n";

			ss << "\n═══════════════════════════════════════════════════════════════\n";

			return ss.str();
		}
	};

} // namespace MML::Systems
#endif // MML_DYNAMICAL_SYSTEM_ANALYZER_H
