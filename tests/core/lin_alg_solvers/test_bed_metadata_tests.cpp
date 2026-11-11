#ifdef MML_SINGLE_HEADER_INCLUDE
#define MML_USE_SINGLE_HEADER
#endif
#include <test_beds/linear_alg_eq_systems_test_bed.h>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <string>
#include <set>
#include <type_traits>

namespace
{
    using namespace MML;
    using namespace MML::TestBeds;

    template<class System>
    void CheckSpectralReferences(const System& system, const Vector<Real>& singularValues,
                                 Real condition, int rank)
    {
        CHECK(system._cond_2 == condition);
        CHECK(system._rank == rank);
        CHECK(system._singular_values.size() == singularValues.size());
        if (system._singular_values.size() == singularValues.size())
            for (int i = 0; i < singularValues.size(); ++i)
                CHECK(system._singular_values[i] == singularValues[i]);
    }
}

TEST_CASE("Symmetric catalog retains existing fixture properties", "[TestBed][metadata]")
{
    const auto a = LinearAlgEqTestBed::getLinAlgEqSystemSymmetric("symm_mat_3x3()");
    CheckSpectralReferences(a, symm_mat_3x3_singular_values(), symm_mat_3x3_cond_2, symm_mat_3x3_rank);
    CHECK(a._determinant == symm_mat_3x3_det);
    CHECK(a._eigen_vectors.size() == symm_mat_3x3_eigen_vecs().size());

    const auto b = LinearAlgEqTestBed::getLinAlgEqSystemSymmetric("symm_mat_5x5()");
    CheckSpectralReferences(b, symm_mat_5x5_singular_values(), symm_mat_5x5_cond_2, symm_mat_5x5_rank);
    CHECK(b._determinant == symm_mat_5x5_det);
    CHECK(b._eigen_vectors.size() == symm_mat_5x5_eigen_vecs().size());

    const auto c = LinearAlgEqTestBed::getLinAlgEqSystemSymmetric("symm_mat_10x10()");
    CheckSpectralReferences(c, symm_mat_10x10_singular_values(), symm_mat_10x10_cond_2, symm_mat_10x10_rank);
    CHECK(c._determinant == symm_mat_10x10_det);
    CHECK(c._eigen_vectors.size() == symm_mat_10x10_eigen_vecs().size());
}

TEST_CASE("SPD catalog and compatibility helpers retain existing fixture properties", "[TestBed][metadata]")
{
    const auto check = [](const std::string& name, const TestLinearSystemSymmetric& helper,
                          const Vector<Real>& values, Real condition, int rank, Real determinant)
    {
        INFO(name);
        const auto system = LinearAlgEqTestBed::getLinAlgEqSystemSymmetric(name);
        CheckSpectralReferences(system, values, condition, rank);
        CheckSpectralReferences(helper, values, condition, rank);
        CHECK(system._determinant == determinant);
        CHECK(helper._determinant == determinant);
    };
    check("spd_3x3", spd_3x3(), spd_3x3_singular_values, spd_3x3_cond_2, spd_3x3_rank, spd_3x3_det);
    check("spd_4x4_correlation", spd_4x4_correlation(), spd_4x4_correlation_singular_values,
          spd_4x4_correlation_cond_2, spd_4x4_correlation_rank, spd_4x4_correlation_det);
    check("spd_5x5_mass_matrix", spd_5x5_mass_matrix(), spd_5x5_mass_matrix_singular_values,
          spd_5x5_mass_matrix_cond_2, spd_5x5_mass_matrix_rank, spd_5x5_mass_matrix_det);
    check("spd_6x6_graph_laplacian", spd_6x6_graph_laplacian(), spd_6x6_graph_laplacian_singular_values,
          spd_6x6_graph_laplacian_cond_2, spd_6x6_graph_laplacian_rank, spd_6x6_graph_laplacian_det);
}

TEST_CASE("Real compatibility helpers retain existing fixture properties", "[TestBed][metadata]")
{
    const auto a = overdetermined_4x3();
    CheckSpectralReferences(a, overdetermined_4x3_singular_values, overdetermined_4x3_cond_2, overdetermined_4x3_rank);
    CHECK(a._mat.rows() == 4);
    CHECK(a._mat.cols() == 3);
    CheckSpectralReferences(overdetermined_5x3_regression(), overdetermined_5x3_regression_singular_values,
                            overdetermined_5x3_regression_cond_2, overdetermined_5x3_regression_rank);
    CheckSpectralReferences(overdetermined_6x4_polynomial(), overdetermined_6x4_polynomial_singular_values,
                            overdetermined_6x4_polynomial_cond_2, overdetermined_6x4_polynomial_rank);

    const auto b = diag_dominant_4x4();
    CheckSpectralReferences(b, diag_dominant_4x4_singular_values, diag_dominant_4x4_cond_2, diag_dominant_4x4_rank);
    CHECK(b._determinant == diag_dominant_4x4_det);
    const auto c = diag_dominant_5x5_tridiag();
    CheckSpectralReferences(c, diag_dominant_5x5_tridiag_singular_values,
                            diag_dominant_5x5_tridiag_cond_2, diag_dominant_5x5_tridiag_rank);
    CHECK(c._determinant == diag_dominant_5x5_tridiag_det);
    const auto d = diag_dominant_6x6_poisson2d();
    CheckSpectralReferences(d, diag_dominant_6x6_poisson2d_singular_values,
                            diag_dominant_6x6_poisson2d_cond_2, diag_dominant_6x6_poisson2d_rank);
    CHECK(d._determinant == diag_dominant_6x6_poisson2d_det);
}

TEST_CASE("Common metadata distinguishes unknown values from stored zero", "[TestBed][metadata]")
{
    static_assert(std::is_aggregate_v<TestCaseMetadata>);
    const TestCaseMetadata unknown;
    CHECK(unknown.source.empty());
    CHECK(unknown.category.empty());
    CHECK_FALSE(unknown.difficulty.has_value());
    CHECK_FALSE(unknown.conditionNumber.has_value());
    CHECK_FALSE(unknown.tolerance.has_value());
    CHECK(unknown.verificationMethod == VerificationMethod::Unspecified);

    const TestCaseMetadata supplied{"example", "example", 1, Real(1), Real(0), VerificationMethod::Analytical};
    REQUIRE(supplied.tolerance.has_value());
    CHECK(*supplied.tolerance == 0.0);
    CHECK(supplied.difficulty == 1);
    CHECK(supplied.conditionNumber == 1.0);
    CHECK(supplied.verificationMethod == VerificationMethod::Analytical);
}

TEST_CASE("Every linear-system catalog entry has honest metadata", "[TestBed][metadata]")
{
    std::set<std::string> seen;
    const auto check = [&](const std::string& name, Real condition, const std::string& source, const std::string& category)
    {
        INFO(name);
        CHECK(seen.insert(name).second);
        auto info = LinearAlgEqTestBed::getMetadata(name);
        CHECK(info.source == source);
        CHECK(info.category == category);
        CHECK_FALSE(info.difficulty.has_value());
        CHECK_FALSE(info.tolerance.has_value());
        CHECK(info.verificationMethod == VerificationMethod::Unspecified);
        if (condition > 0)
            CHECK(info.conditionNumber == condition);
        else
            CHECK_FALSE(info.conditionNumber.has_value());
        info.source = "changed snapshot";
        CHECK(LinearAlgEqTestBed::getMetadata(name).source == source);
    };

    const auto realGroup = [&](const auto& entries, const std::string& source, const std::string& category)
    {
        for (const auto& entry : entries)
        {
            check(entry.first, entry.second->_cond_2, source, category);
            CHECK(entry.second == &LinearAlgEqTestBed::getLinAlgEqSystem(entry.first));
        }
    };
    realGroup(LinearAlgEqTestBed::getRandomMatrices(), "test_data/linear_alg_eq_systems_real_defs.h", "random");
    const std::string classic = "test_data/linear_alg_eq_systems_classic_defs.h";
    realGroup(LinearAlgEqTestBed::getHilbertMatrices(), classic, "hilbert");
    realGroup(LinearAlgEqTestBed::getPascalMatrices(), classic, "pascal");
    realGroup(LinearAlgEqTestBed::getVandermondeMatrices(), classic, "vandermonde");
    realGroup(LinearAlgEqTestBed::getFrankMatrices(), classic, "frank");
    realGroup(LinearAlgEqTestBed::getKahanMatrices(), classic, "kahan");
    realGroup(LinearAlgEqTestBed::getDiagDominantMatrices(),
              "test_data/linear_alg_eq_systems_solverspecific_defs.h", "diagonally-dominant");
    CHECK(seen.size() == static_cast<std::size_t>(LinearAlgEqTestBed::numLinAlgEqSystems()));

    for (int i = 0; i < LinearAlgEqTestBed::numLinAlgEqSystemsSymmetric(); ++i)
    {
        const auto entries = LinearAlgEqTestBed::getAllSymmetricSystems();
        const auto& entry = entries[i];
        const bool spd = entry.first.find("spd_") == 0;
        check(entry.first, entry.second->_cond_2,
              spd ? "test_data/linear_alg_eq_systems_solverspecific_defs.h" : "test_data/linear_alg_eq_systems_sym_defs.h",
              spd ? "spd" : "symmetric");
        CHECK(entry.second == &LinearAlgEqTestBed::getLinAlgEqSystemSymmetric(i));
    }
    for (const auto& entry : LinearAlgEqTestBed::getAllComplexSystems())
        check(entry.first, entry.second->_cond_2, "test_data/linear_alg_eq_systems_complex_defs.h", "complex");
    for (const auto& entry : LinearAlgEqTestBed::getAllOverdeterminedSystems())
        check(entry.first, entry.second->_cond_2, "test_data/linear_alg_eq_systems_solverspecific_defs.h", "overdetermined");
    check("mat_5x5_multi_rhs1", LinearAlgEqTestBed::getLinAlgEqSystemMultiRHS(0)._cond_2,
          "test_beds/linear_alg_eq_systems_test_bed.h", "multi-rhs");

    CHECK(seen.size() == static_cast<std::size_t>(LinearAlgEqTestBed::numLinAlgEqSystems()
          + LinearAlgEqTestBed::numLinAlgEqSystemsSymmetric() + LinearAlgEqTestBed::numComplexSystems()
          + LinearAlgEqTestBed::numOverdeterminedSystems() + LinearAlgEqTestBed::numLinAlgEqSystemsMultiRHS()));
    CHECK_THROWS_AS(LinearAlgEqTestBed::getMetadata("absent"), std::runtime_error);
    CHECK_THROWS_AS(LinearAlgEqTestBed::getMetadata(""), std::runtime_error);
    CHECK_THROWS_AS(LinearAlgEqTestBed::getMetadata("symm_mat_3x3"), std::runtime_error);
}

TEST_CASE("Missing references are not confused with singular reference cases", "[TestBed][metadata]")
{
    for (const std::string name : {"mat_8x8", "mat_10x10", "mat_20x20"})
    {
        INFO(name);
        CHECK_FALSE(LinearAlgEqTestBed::getMetadata(name).conditionNumber.has_value());
        const auto system = LinearAlgEqTestBed::getLinAlgEqSystem(name);
        CHECK(system._cond_2 == 0.0);
        CHECK(system._singular_values.size() == 0);
        CHECK(system._rank == system._n);  // Legacy assumption, not a reference rank.
        for (const auto& group : {LinearAlgEqTestBed::getWellConditioned(),
                                 LinearAlgEqTestBed::getModeratelyConditioned(),
                                 LinearAlgEqTestBed::getIllConditioned()})
            for (const auto& entry : group)
                CHECK(entry.first != name);
    }
    const auto singular = LinearAlgEqTestBed::getLinAlgEqSystemSymmetric("spd_6x6_graph_laplacian");
    REQUIRE(LinearAlgEqTestBed::getMetadata("spd_6x6_graph_laplacian").conditionNumber.has_value());
    CHECK(std::isinf(*LinearAlgEqTestBed::getMetadata("spd_6x6_graph_laplacian").conditionNumber));
    CHECK(singular._rank == 5);
    CHECK(singular._determinant == 0.0);
    REQUIRE(singular._singular_values.size() == 6);
    CHECK(singular._singular_values[5] == 0.0);
}

TEST_CASE("Multi-RHS case reuses properties of its existing coefficient matrix", "[TestBed][metadata]")
{
    const auto* system = &LinearAlgEqTestBed::getLinAlgEqSystemMultiRHS("mat_5x5_multi_rhs1");
    CheckSpectralReferences(*system, mat_5x5_singular_values(), mat_5x5_cond_2, mat_5x5_rank);
    CHECK(system->_determinant == mat_5x5_det);
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j)
            CHECK(system->_mat(i, j) == mat_5x5()(i, j));
    CHECK(system->_rhs.cols() == 2);
    CHECK(system->_sol.cols() == 2);
    CHECK(system == &LinearAlgEqTestBed::getLinAlgEqSystemMultiRHS(0));
}

TEST_CASE("Legacy linear-system constructors and catalog filters remain compatible", "[TestBed][metadata]")
{
    const TestLinearSystem real(3, mat_3x3(), mat_3x3_rhs0(), mat_3x3_rhs0_sol(), mat_3x3_eigen_val());
    const TestLinearSystem realWithVectors(3, mat_3x3(), mat_3x3_rhs0(), mat_3x3_rhs0_sol(),
                                           mat_3x3_eigen_val(), mat_3x3_eigen_vecs());
    const TestLinearSystemSymmetric symmetric(3, symm_mat_3x3(), symm_mat_3x3_rhs0(),
                                              symm_mat_3x3_rhs0_sol(), symm_mat_3x3_eigen_val());
    const TestLinearSystemSymmetric symmetricWithVectors(3, symm_mat_3x3(), symm_mat_3x3_rhs0(),
        symm_mat_3x3_rhs0_sol(), symm_mat_3x3_eigen_val(), symm_mat_3x3_eigen_vecs());
    const auto& complexCase = LinearAlgEqTestBed::getComplexSystem(0);
    const TestLinearSystemComplex complex(2, complexCase._mat, complexCase._rhs, complexCase._sol);
    const TestLinearSystemOverdetermined rectangular(4, 3, overdetermined_4x3_mat,
                                                     overdetermined_4x3_rhs, overdetermined_4x3_sol);
    const TestLinearSystemMultiRHS multi(5, mat_5x5(), mat_5x5_rhs_multi(), mat_5x5_rhs_multi());
    CHECK(real._rank == 3);
    CHECK(realWithVectors._rank == 3);
    CHECK(symmetric._rank == 3);
    CHECK(symmetricWithVectors._rank == 3);
    CHECK(complex._rank == 2);
    CHECK(rectangular._rank == 3);
    CHECK(multi._cond_2 == 0.0);
    CHECK(realWithVectors._eigen_vectors.size() == mat_3x3_eigen_vecs().size());
    CHECK(symmetricWithVectors._eigen_vectors.size() == symm_mat_3x3_eigen_vecs().size());

    std::set<std::string> classified;
    for (const auto& group : {LinearAlgEqTestBed::getWellConditioned(),
                             LinearAlgEqTestBed::getModeratelyConditioned(), LinearAlgEqTestBed::getIllConditioned()})
        for (const auto& entry : group)
        {
            CHECK(classified.insert(entry.first).second);
            CHECK(entry.second == &LinearAlgEqTestBed::getLinAlgEqSystem(entry.first));
        }
    CHECK(classified.size() + 3 == LinearAlgEqTestBed::getAllSystems().size());
    for (int i = 0; i < LinearAlgEqTestBed::numLinAlgEqSystems(); ++i)
        CHECK(&LinearAlgEqTestBed::getLinAlgEqSystem(i) == LinearAlgEqTestBed::getAllSystems()[i].second);
}
