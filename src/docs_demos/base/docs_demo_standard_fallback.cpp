#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/StandardFunctions.h>
#endif

#include <iostream>

using namespace MML;

void Docs_Demo_StandardFallbackParity()
{
    std::cout << "Standard special-function backend: "
              << (MML_HAS_STD_SPECIAL_FUNCTIONS ? "native" : "portable subset") << '\n';
    std::cout << "Complete elliptic E uses the shared portable AGM implementation\n";
    const Real associated = Functions::AssocLegendre(1, 1, 0);
    const Real spherical = Functions::SphLegendre(1, 1, Constants::PI / 2);
    const Real elliptic = Functions::Comp_ellint_2(Real(0.5));
    std::cout << associated << '\n';  // 1: no Condon-Shortley phase
    std::cout << spherical << '\n';   // approximately -0.345494
    std::cout << elliptic << '\n';    // approximately 1.467462209339427

    if (associated != 1 ||
        !IsNearlyEqual(spherical, -std::sqrt(3 / (8 * Constants::PI)), Precision::RoundoffZeroThreshold) ||
        !IsNearlyEqual(elliptic, Real(1.467462209339427155459795266990916L), Precision::RoundoffZeroThreshold) ||
        !IsNearlyEqual(Functions::Comp_ellint_2(Real(0.9L)),
                       Real(1.171697052781614141185913957957410L), Precision::RoundoffZeroThreshold) ||
        !IsNearlyEqual(Functions::Comp_ellint_2(Real(0.9375L)),
                       Real(1.1215931702323738421144018082563278L), Precision::RoundoffZeroThreshold) ||
        Functions::Comp_ellint_2(Real(1)) != 1 || Functions::Comp_ellint_2(Real(-1)) != 1)
        throw NumericalMethodError("Standard function fallback documentation example failed");
}
