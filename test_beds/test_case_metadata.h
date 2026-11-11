#ifndef MML_TEST_CASE_METADATA_H
#define MML_TEST_CASE_METADATA_H

#include <mml/MMLTypeDefs.h>

#include <optional>
#include <string>

namespace MML::TestBeds
{
    enum class VerificationMethod
    {
        Unspecified,
        Analytical,
        ExternalReference,
        IndependentNumerical
    };

    struct TestCaseMetadata
    {
        std::string source;
        std::string category;
        std::optional<int> difficulty;  // Existing catalog scale: 1 (easy) through 4 (expert).
        std::optional<Real> conditionNumber;
        std::optional<Real> tolerance;
        VerificationMethod verificationMethod = VerificationMethod::Unspecified;
    };
}

#endif
