#include <test_beds/linear_alg_eq_systems_test_bed.h>

#include <cmath>
#include <iostream>
#include <stdexcept>

void Docs_Demo_TestBedMetadata()
{
    using namespace MML::TestBeds;
    const auto info = LinearAlgEqTestBed::getMetadata("symm_mat_3x3()");
    if (!info.conditionNumber || info.verificationMethod != VerificationMethod::Unspecified)
        throw std::runtime_error("Unexpected symmetric fixture metadata");
    std::cout << info.source << ": condition number = " << *info.conditionNumber << '\n';

    const auto unknown = LinearAlgEqTestBed::getMetadata("mat_20x20");
    if (unknown.conditionNumber || unknown.difficulty || unknown.tolerance)
        throw std::runtime_error("Unavailable references must remain unknown");
    std::cout << "mat_20x20: condition number unknown\n";

    const auto singular = LinearAlgEqTestBed::getMetadata("spd_6x6_graph_laplacian");
    if (!singular.conditionNumber || !std::isinf(*singular.conditionNumber))
        throw std::runtime_error("Singular reference must retain its infinite condition number");
    std::cout << "spd_6x6_graph_laplacian: infinite condition number\n";
}
