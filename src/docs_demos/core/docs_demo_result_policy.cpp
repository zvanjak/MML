///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        docs_demo_result_policy.cpp                                         ///
///  Description: Runnable counterpart of docs/RESULT_POLICY.md (section 6 examples). ///
///               Every snippet in that page must compile and run here.               ///
///                                                                                   ///
///  MinimalMathLibrary-wewi.3.1                                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>

#include <mml/base/AlgorithmTypes.h>
#include <mml/base/Function.h>
#include <mml/core/Integration.h>
#include <mml/algorithms/RootFinding.h>
#include <mml/tools/Persistence.h>
#include <mml/tools/DataLoader.h>
#include <mml/MMLExceptions.h>
#endif

#include <cassert>
#include <iostream>
#include <type_traits>
#include <sstream>

using namespace MML;

void Docs_Demo_ExceptionOwnership()
{
	const FourierError fourier("Fourier validation");
	const StatisticsError statistics("Statistics validation");
	const CurveFittingError fitting("Fitting validation");
	for (const MMLException* error : {static_cast<const MMLException*>(&fourier),
	                                static_cast<const MMLException*>(&statistics),
	                                static_cast<const MMLException*>(&fitting)})
		std::cout << "Centralized shared exception: " << error->message() << '\n';
	try {
		throw StatisticsError("same type in core and packages");
	} catch (const MMLException& error) {
		std::cout << "Caught through MMLException: " << error.message() << '\n';
	}
}

void Docs_Demo_ResultPolicy()
{
    std::cout << "=== RESULT_POLICY.md - worked examples ===\n";

    RealFunctionFromStdFunc f([](Real x) { return std::cos(x) - x; });
    RealFunctionFromStdFunc g([](Real x) { return std::exp(-x * x); });
    const Real a = 0.0, b = 1.0;

    // --- Iterative: Success only after the check -------------------------------------
    // (RootFinding spells its Detailed form as the config overload)
    RootFinding::RootFindingConfig cfg;
    cfg.tolerance = 1e-12;
    auto r = RootFinding::FindRootBrentDetailed(f, 0.0, 2.0, cfg);
    if (!r) {
        // r.status, r.error_message, r.root is the best bracket midpoint
        std::cout << "  Brent failed: " << to_string(r.status) << " - " << r.error_message << "\n";
    } else {
        std::cout << "  Brent: root = " << r.root << ", status = " << to_string(r.status)
                  << ", converged = " << std::boolalpha << r.converged
                  << ", iterations = " << r.iterations_used << "\n";
    }
    assert(r.converged == (r.status == AlgorithmStatus::Success));

    // --- One-shot fixed rule: Completed, no estimate ---------------------------------
    auto gauss = IntegrateGauss10Detailed(g, a, b);
    assert(gauss.status == AlgorithmStatus::Completed);
    assert(!gauss.converged && !HasEstimate(gauss.error_estimate));
    assert(gauss.IsSuccess());          // the value is usable ...
    std::cout << "  Gauss10: value = " << gauss.value << ", status = " << to_string(gauss.status)
              << ", error estimate present = " << HasEstimate(gauss.error_estimate) << "\n";

    // --- Adaptive rule: Success with a real estimate ---------------------------------
    auto t = IntegrateTrapDetailed(g, a, b);
    if (t.converged)
        std::cout << "  Trap:    value = " << t.value << " +/- " << t.error_estimate
                  << " after " << t.iterations << " refinements\n";
    assert(HasEstimate(t.error_estimate));

    // --- Plain call returns the value; result conversion stays explicit ------------------
    Real I1 = IntegrateSimpson(g, a, b);                              // value, or ConvergenceError
    Real I2 = IntegrateSimpsonDetailed(g, a, b).value;                // explicit opt-out from a result
    Real I3 = static_cast<Real>(IntegrateSimpsonDetailed(g, a, b));  // likewise
    static_assert(!std::is_convertible_v<IntegrationResult, Real>,
                  "RESULT_POLICY section 4: no implicit conversion from a result");
    std::cout << "  Simpson: " << I1 << " (plain) == " << I2 << " (.value) == " << I3 << " (static_cast)\n";

    // --- Status vocabulary -----------------------------------------------------------
    std::cout << "  IsUsable(Success) = " << IsUsable(AlgorithmStatus::Success)
              << ", IsUsable(Completed) = " << IsUsable(AlgorithmStatus::Completed)
              << ", IsUsable(MaxIterationsExceeded) = " << IsUsable(AlgorithmStatus::MaxIterationsExceeded) << "\n";

    // --- Building results with the helpers -------------------------------------------
    struct DemoResult : IterativeResultBase { Real payload = 0.0; };
    auto ok = MakeSuccessResult<DemoResult>(7, 1e-13, "Demo");
    auto done = MakeCompletedResult<DemoResult>("DemoDirect");
    auto bad = MakeFailureResult<DemoResult>(AlgorithmStatus::Stalled, "no progress after 12 iterations", 12, "Demo");
    std::cout << "  helpers: " << to_string(ok.status) << " (converged " << ok.converged << "), "
              << to_string(done.status) << " (converged " << done.converged
              << ", achieved_tolerance has estimate " << HasEstimate(done.achieved_tolerance) << "), "
              << to_string(bad.status) << " (usable " << static_cast<bool>(bad) << ")\n";

		std::stringstream storage;
		IOResult saved = Persistence::SaveJson(storage, Vector<Real>{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)});
		Vector<Real> restored;
		IOResult loaded = Persistence::LoadJson(storage, restored);
		if (!saved || !loaded || restored != Vector<Real>{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)})
			throw FileIOError("IOResult stream round-trip failed");
		static_assert(!std::is_convertible_v<IOResult, bool>);
		std::cout << "  IOResult: " << restored.size() << " values restored; explicit success check\n";

		Data::Dataset dataset;
		const IOResult datasetResult = Data::LoadCSVSafe("", dataset);
		if (datasetResult || datasetResult.error != IOError::FILE_NOT_OPENED)
			throw FileIOError("Dataset IOResult failed-open contract was not preserved");
}
