// Test that MML single headers work standalone.
// This file explicitly includes ONLY the selected single header - no other MML headers.

#include <catch2/catch_all.hpp>
#include <type_traits>

// Force single header mode and include the selected variant directly.
#ifndef MML_SINGLE_HEADER_INCLUDE
#define MML_SINGLE_HEADER_INCLUDE <MML.h>
#endif
#define MML_USE_SINGLE_HEADER
#include MML_SINGLE_HEADER_INCLUDE

#if defined(MML_SINGLE_HEADER_EXPECT_FLOAT)
static_assert(std::is_same_v<Real, float>);
#elif defined(MML_SINGLE_HEADER_EXPECT_LONG_DOUBLE)
static_assert(std::is_same_v<Real, long double>);
#else
static_assert(std::is_same_v<Real, double>);
#endif

using namespace MML;

namespace MML::Tests::SingleHeader
{
    using namespace MML::Systems;

    constexpr Real single_header_tolerance()
    {
        if constexpr (std::is_same_v<Real, float>)
            return REAL(1e-4);
        else
            return REAL(1e-8);
    }

    template<typename Scalar>
    constexpr auto single_header_matrix_tolerance()
    {
        using Magnitude = MatrixAlg::MatrixMagnitude<Scalar>;
        const auto tolerance = static_cast<Magnitude>(single_header_tolerance());
        return MatrixAlg::MatrixComparisonTolerance<Magnitude>{tolerance, tolerance};
    }

    template<typename MatrixType>
    concept HasLegacyMatrixIsDiagonal = requires(const MatrixType& matrix) { matrix.isDiagonal(); };

    template<typename MatrixType>
    concept HasLegacyMatrixIsDiagonallyDominant = requires(const MatrixType& matrix) { matrix.isDiagonallyDominant(); };

    template<typename MatrixType>
    concept HasLegacyMatrixIsSymmetric = requires(const MatrixType& matrix) { matrix.isSymmetric(); };

    template<typename MatrixType>
    concept HasLegacyMatrixNormL1 = requires(const MatrixType& matrix) { matrix.NormL1(); };

    template<typename MatrixType>
    concept HasLegacyMatrixNormL2 = requires(const MatrixType& matrix) { matrix.NormL2(); };

    template<typename MatrixType>
    concept HasLegacyMatrixNormLInf = requires(const MatrixType& matrix) { matrix.NormLInf(); };

    template<typename MatrixType>
    concept HasLegacyIsOrthonormalColumns = requires(const MatrixType& matrix) { IsOrthonormalColumns(matrix); };

    template<typename MatrixType>
    concept HasLegacyAnalyzeMatrix = requires(const MatrixType& matrix) { AnalyzeMatrix(matrix); };

    template<typename Analyzer>
    concept HasLegacyAnalyzerGetEigen = requires(const Analyzer& analyzer) { analyzer.GetEigen(); };

    template<typename Analyzer>
    concept HasLegacyAnalyzerEigenvaluesSymmetric = requires(const Analyzer& analyzer) { analyzer.EigenvaluesSymmetric(); };

    template<typename System>
    concept HasLegacySystemGetEigen = requires(const System& system) { system.GetEigen(); };

    template<typename System>
    concept HasLegacySystemEigenvaluesSymmetric = requires(const System& system) { system.EigenvaluesSymmetric(); };

    template<typename System>
    concept HasLegacySystemIsSymmetric = requires(const System& system) { system.isSymmetric(); };

    struct SingleHeaderGF4Modulus
    {
        static Algebra::Polynomial<Algebra::PrimeFieldElement<2>> modulus()
        {
            return {1, 1, 1};
        }
    };

    class SingleHeaderCyclicGroup3
    {
    public:
        using element_type = int;

        const std::vector<int>& elements() const
        {
            static const std::vector<int> values{0, 1, 2};
            return values;
        }

        int identity() const { return 0; }
        int compose(int left, int right) const { return (left + right) % 3; }
        int inverse(int element) const { return (3 - element) % 3; }
    };

    class SingleHeaderRealFunction : public IRealFunction
    {
    public:
        Real operator()(Real x) const override { return x * x; }
    };

    class SingleHeaderComplexFunction : public IComplexFunction
    {
    public:
        Complex operator()(Complex z) const override { return z * z; }
    };

    class SingleHeaderScalarFunction : public IScalarFunction<2>
    {
    public:
        Real operator()(const VectorN<Real, 2>& point) const override { return point[0] * point[0] + point[1]; }
    };

    class SingleHeaderVectorFunction : public IVectorFunction<2>
    {
    public:
        VectorN<Real, 2> operator()(const VectorN<Real, 2>& point) const override { return {point[0] * point[0], point[0] + point[1]}; }
    };

    class SingleHeaderCurve : public IParametricCurve<2>
    {
    public:
        VectorN<Real, 2> operator()(Real t) const override { return {t * t, t}; }
        Real getMinT() const override { return REAL(-10.0); }
        Real getMaxT() const override { return REAL(10.0); }
    };

    class SingleHeaderSurface : public IParametricSurfaceRect<2>
    {
    public:
        VectorN<Real, 2> operator()(Real u, Real w) const override { return {u * u, u + w}; }
        Real getMinU() const override { return REAL(-10.0); }
        Real getMaxU() const override { return REAL(10.0); }
        Real getMinW() const override { return REAL(-10.0); }
        Real getMaxW() const override { return REAL(10.0); }
    };

    class SingleHeaderTensorField : public ITensorField2<2>
    {
    public:
        SingleHeaderTensorField() : ITensorField2<2>(0, 2) {}
        Real Component(int, int, const VectorN<Real, 2>& point) const override { return point[0] * point[0] + point[1]; }
        Tensor2<2> operator()(const VectorN<Real, 2>&) const override { return Tensor2<2>(2, 0); }
    };

    TEST_CASE("SingleHeader_FirstDerivativeAdapters", "[single_header][derivation]")
    {
        const Real h = REAL(1e-3);
        const VectorN<Real, 2> point{REAL(2.0), REAL(1.0)};
        SingleHeaderRealFunction real_function;
        SingleHeaderComplexFunction complex_function;
        SingleHeaderScalarFunction scalar_function;
        SingleHeaderVectorFunction vector_function;
        SingleHeaderCurve curve;
        SingleHeaderSurface surface;
        SingleHeaderTensorField tensor_field;

        const auto expected_derivative = Catch::Approx(REAL(4.0)).epsilon(static_cast<double>(single_header_tolerance()));
        REQUIRE(Derivation::NDer4(real_function, REAL(2.0), h) == expected_derivative);
        REQUIRE(Derivation::NDer4Complex(complex_function, Complex(REAL(2.0)), h).real() == expected_derivative);
        REQUIRE(Derivation::NDer4Partial(scalar_function, 0, point, h) == expected_derivative);
        REQUIRE(Derivation::NDer4Partial(vector_function, 0, 0, point, h) == expected_derivative);
        REQUIRE(Derivation::NDer4(curve, REAL(2.0), h)[0] == expected_derivative);
        REQUIRE(Derivation::NDer2_u(surface, REAL(2.0), REAL(1.0), h)[0] == expected_derivative);
        REQUIRE(Derivation::NDer4Partial(tensor_field, 0, 0, 0, point, h) == expected_derivative);
    }

    TEST_CASE("SingleHeader_Algebra_GroupLaws", "[single_header][Algebra]")
    {
        REQUIRE(Algebra::CheckGroupLaws(SingleHeaderCyclicGroup3{}));
    }

    TEST_CASE("SingleHeader_Algebra_Permutation", "[single_header][Algebra]")
    {
        const auto rotation = Algebra::Permutation<4>::from_cycles({{0, 1, 2, 3}});
        const std::array<int, 4> vertices{10, 20, 30, 40};

        REQUIRE(rotation.order() == 4);
        REQUIRE(rotation.compose(rotation.inverse()) == Algebra::Permutation<4>::identity());
        REQUIRE(rotation.apply(vertices) == std::array<int, 4>{40, 10, 20, 30});
    }

    TEST_CASE("SingleHeader_Algebra_FiniteGroups", "[single_header][Algebra]")
    {
        Algebra::CyclicGroup cyclic(5);
        Algebra::DihedralGroup dihedral(4);
        const auto rotation = dihedral.generator();
        const auto reflection = dihedral.reflection();

        REQUIRE(Algebra::element_order(cyclic, cyclic.generator()) == 5);
        REQUIRE(dihedral.compose(dihedral.compose(reflection, rotation), reflection) ==
            dihedral.inverse(rotation));
        REQUIRE(Algebra::CheckGroupLaws(dihedral));

        const auto graph = Algebra::MakeCayleyGraphData(
            dihedral, std::vector<Algebra::DihedralElement>{rotation, reflection});
        REQUIRE(graph.elements.size() == 8);
        REQUIRE(graph.edges.size() == 16);
    }

    TEST_CASE("SingleHeader_Algebra_GroupActions", "[single_header][Algebra]")
    {
        Algebra::CyclicGroup group(3);
        auto action = Algebra::MakeGroupAction<Algebra::CyclicElement, int>(
            [&group](const Algebra::CyclicElement& element, const int& vertex) {
                return group.permutation(element).apply(vertex);
            });

        REQUIRE(Algebra::Orbit(group, action, 0).size() == 3);
        REQUIRE(Algebra::Stabilizer(group, action, 0).size() == 1);
        REQUIRE(Algebra::Orbit(group, action, 0).size() *
            Algebra::Stabilizer(group, action, 0).size() == group.order());
    }

    TEST_CASE("SingleHeader_Algebra_PrimeFields", "[single_header][Algebra]")
    {
        using F5 = Algebra::PrimeFieldElement<5>;
        Algebra::PrimeField<5> field;
        MatrixNM<F5, 2, 2> matrix{{1, 2}, {3, 4}};
        const auto inverse = Algebra::FieldMatrixInverse(matrix);
        const auto identity = matrix * inverse;

        REQUIRE(Algebra::CheckFieldLaws(field));
        REQUIRE(F5(3).inverse() == F5(2));
        REQUIRE(Algebra::PrimitiveRoot<5>() == F5(2));
        REQUIRE(identity(0, 0) == F5(1));
        REQUIRE(identity(0, 1) == F5(0));
        REQUIRE(identity(1, 0) == F5(0));
        REQUIRE(identity(1, 1) == F5(1));
    }

    TEST_CASE("SingleHeader_Algebra_ExtensionFields", "[single_header][Algebra]")
    {
        using GF4 = Algebra::FiniteFieldElement<2, 2, SingleHeaderGF4Modulus>;
        const GF4 alpha{0, 1};
        Algebra::ExtensionField<2, 2, SingleHeaderGF4Modulus> field;

        REQUIRE(alpha * alpha == alpha + GF4(1));
        REQUIRE(alpha.inverse() == alpha + GF4(1));
        REQUIRE(Algebra::IsIrreducible(SingleHeaderGF4Modulus::modulus()));
        REQUIRE(Algebra::CheckFieldLaws(field));
    }

    TEST_CASE("SingleHeader_Algebra_Representations", "[single_header][Algebra]")
    {
        Algebra::CyclicGroup group(2);
        auto representation = Algebra::MakeRepresentation<Algebra::CyclicElement, Real, 2>(
            [](const Algebra::CyclicElement& element) {
                return element.exponent == 0
                    ? MatrixNM<Real, 2, 2>::Identity()
                    : MatrixNM<Real, 2, 2>{{1, 0}, {0, -1}};
            });
        const auto projection = Algebra::InvariantProjection(group, representation);
        const auto vector = Algebra::SymmetrizeVector(group, representation, VectorN<Real, 2>{2, 3});

        REQUIRE(Algebra::VerifyRepresentation(group, representation));
        REQUIRE(projection * projection == projection);
        REQUIRE(vector == VectorN<Real, 2>{2, 0});
        REQUIRE(representation.apply(group.generator(), vector) == vector);
    }

    TEST_CASE("SingleHeader_Algebra_LieGroups", "[single_header][Algebra]")
    {
        const Algebra::SO2 planar(Constants::PI / REAL(2.0));
        const auto spatial = Algebra::SO3::Exp({REAL(0.2), REAL(-0.1), REAL(0.3)});
        MatrixNM<Real, 3, 3> drifted = spatial.matrix();
        drifted(0, 0) += REAL(1e-5);

        REQUIRE(planar.apply({1, 0}).IsEqualTo({0, 1}, single_header_tolerance()));
        REQUIRE(spatial.log().IsEqualTo({REAL(0.2), REAL(-0.1), REAL(0.3)}, single_header_tolerance()));
        REQUIRE(Algebra::SO3::IsRotationMatrix(spatial.matrix()));
        REQUIRE(Algebra::SO3::IsRotationMatrix(Algebra::SO3::Project(drifted).matrix()));
    }

    TEST_CASE("SingleHeader_Algebra_RigidMotions", "[single_header][Algebra]")
    {
        const Algebra::SE3 transform(
            Algebra::SO3::FromAxisAngle({0, 0, 1}, Constants::PI / REAL(2.0)), {2, 3, 4});
        const auto restored = Algebra::SE3::FromHomogeneousMatrix(transform.homogeneous_matrix());

        REQUIRE(transform.apply_point({1, 0, 0}).IsEqualTo({2, 4, 4}, single_header_tolerance()));
        REQUIRE(transform.apply_vector({1, 0, 0}).IsEqualTo({0, 1, 0}, single_header_tolerance()));
        REQUIRE(restored.homogeneous_matrix().IsEqualTo(transform.homogeneous_matrix(), single_header_tolerance()));
        REQUIRE(transform.inverse().apply_point(transform.apply_point({1, 2, 3})).IsEqualTo({1, 2, 3}, single_header_tolerance()));
    }

    TEST_CASE("SingleHeader_Vector_Basic", "[single_header]")
    {
        Vector<Real> v1({ REAL(1.0), REAL(2.0), REAL(3.0) });
        Vector<Real> v2({ REAL(4.0), REAL(5.0), REAL(6.0) });
        
        // Vector addition
        auto v3 = v1 + v2;
        REQUIRE(v3.size() == 3);
        REQUIRE(v3[0] == REAL(5.0));
        REQUIRE(v3[1] == REAL(7.0));
        REQUIRE(v3[2] == REAL(9.0));
        
        // Scalar product (free function in VectorUtils)
        Real dot = Utils::ScalarProduct(v1, v2);
        REQUIRE(dot == REAL(32.0)); // 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
    }
    
    TEST_CASE("SingleHeader_Matrix_Basic", "[single_header]")
    {
        Matrix<Real> A(2, 3, {
            REAL(1.0), REAL(2.0), REAL(3.0),
            REAL(4.0), REAL(5.0), REAL(6.0)
        });
        
        Matrix<Real> B(3, 2, {
            REAL(7.0), REAL(8.0),
            REAL(9.0), REAL(10.0),
            REAL(11.0), REAL(12.0)
        });
        
        // Matrix multiplication: (2x3) * (3x2) = (2x2)
        auto C = A * B;
        
        REQUIRE(C.rows() == 2);
        REQUIRE(C.cols() == 2);
        
        // C[0][0] = 1*7 + 2*9 + 3*11 = 7 + 18 + 33 = 58
        REQUIRE(C[0][0] == REAL(58.0));
        // C[0][1] = 1*8 + 2*10 + 3*12 = 8 + 20 + 36 = 64
        REQUIRE(C[0][1] == REAL(64.0));
        // C[1][0] = 4*7 + 5*9 + 6*11 = 28 + 45 + 66 = 139
        REQUIRE(C[1][0] == REAL(139.0));
        // C[1][1] = 4*8 + 5*10 + 6*12 = 32 + 50 + 72 = 154
        REQUIRE(C[1][1] == REAL(154.0));
    }
    
    TEST_CASE("SingleHeader_MatrixVector_Multiply", "[single_header]")
    {
        Matrix<Real> A(2, 3, {
            REAL(1.0), REAL(2.0), REAL(3.0),
            REAL(4.0), REAL(5.0), REAL(6.0)
        });
        
        Vector<Real> v({ REAL(1.0), REAL(2.0), REAL(3.0) });
        
        // Matrix-vector multiplication: (2x3) * (3x1) = (2x1)
        auto result = A * v;
        
        REQUIRE(result.size() == 2);
        // result[0] = 1*1 + 2*2 + 3*3 = 1 + 4 + 9 = 14
        REQUIRE(result[0] == REAL(14.0));
        // result[1] = 4*1 + 5*2 + 6*3 = 4 + 10 + 18 = 32
        REQUIRE(result[1] == REAL(32.0));
    }

    TEST_CASE("SingleHeader_ComplexQRAndHermitianEigen", "[single_header][complex_linalg]")
    {
        Matrix<Complex> matrix(2, 2, {
            {2.0, 0.0}, {0.0, 1.0},
            {0.0, -1.0}, {3.0, 0.0}
        });
        Vector<Complex> expected{ {1.0, 1.0}, {-2.0, 0.5} };

        QRSolver<Complex> qr(matrix);
        REQUIRE(qr.Solve(matrix * expected).IsEqualTo(expected, single_header_tolerance()));

        auto eigen = HermitianMatEigenSolverJacobi::Solve(matrix);
        REQUIRE(eigen.converged);
        REQUIRE(eigen.eigenvalues.size() == 2);
        REQUIRE(MatrixAlg::IsUnitary(eigen.eigenvectors, single_header_matrix_tolerance<Complex>()));
    }

    TEST_CASE("SingleHeader_MatrixAnalysisResultTypes", "[single_header][matrix_analysis][MML2]")
    {
        STATIC_REQUIRE_FALSE(HasLegacyMatrixIsDiagonal<Matrix<Real>>);
        STATIC_REQUIRE_FALSE(HasLegacyMatrixIsDiagonallyDominant<Matrix<Real>>);
        STATIC_REQUIRE_FALSE(HasLegacyMatrixIsSymmetric<Matrix<Real>>);
        STATIC_REQUIRE_FALSE(HasLegacyMatrixNormL1<Matrix<Real>>);
        STATIC_REQUIRE_FALSE(HasLegacyMatrixNormL2<Matrix<Real>>);
        STATIC_REQUIRE_FALSE(HasLegacyMatrixNormLInf<Matrix<Real>>);
        STATIC_REQUIRE_FALSE(HasLegacyIsOrthonormalColumns<Matrix<Real>>);
        STATIC_REQUIRE_FALSE(HasLegacyAnalyzeMatrix<Matrix<Real>>);
        STATIC_REQUIRE_FALSE(HasLegacyAnalyzerGetEigen<MatrixAnalyzer<Real>>);
        STATIC_REQUIRE_FALSE(HasLegacyAnalyzerEigenvaluesSymmetric<MatrixAnalyzer<Real>>);
        STATIC_REQUIRE_FALSE(HasLegacySystemGetEigen<Systems::LinearSystem<Real>>);
        STATIC_REQUIRE_FALSE(HasLegacySystemEigenvaluesSymmetric<Systems::LinearSystem<Real>>);
        STATIC_REQUIRE_FALSE(HasLegacySystemIsSymmetric<Systems::LinearSystem<Real>>);

        STATIC_REQUIRE(std::same_as<MatrixAlg::MatrixMagnitude<Complex>, Real>);
        STATIC_REQUIRE(std::same_as<MatrixAlg::MatrixComplexScalar<Real>, Complex>);
        STATIC_REQUIRE(std::same_as<
            decltype(std::declval<MatrixAlg::SVDDecomposition<Complex>>().singularValues),
            Vector<Real>>);
        STATIC_REQUIRE(std::same_as<
            decltype(std::declval<MatrixAlg::EigensystemResult<Real>>().eigenvectors),
            Matrix<Complex>>);
        STATIC_REQUIRE(std::same_as<
            decltype(std::declval<MatrixAlg::MatrixAnalysis<Complex>>().determinant),
            std::optional<Complex>>);

        const MatrixAlg::MatrixAnalysis<Real> analysis;
        REQUIRE(analysis.isSquare);
        REQUIRE(analysis.stability == MatrixAlg::MatrixStability::Singular);
    }

    TEST_CASE("SingleHeader_MatrixAlgStructuralAndScalarOperations", "[single_header][matrix_analysis][MML2]")
    {
        const Matrix<Complex> hermitian{2, 2, {
            Complex{2.0, 0.0}, Complex{0.0, 1.0},
            Complex{0.0, -1.0}, Complex{3.0, 0.0}
        }};
        REQUIRE(MatrixAlg::IsHermitian(hermitian));
        REQUIRE_FALSE(MatrixAlg::IsSymmetric(hermitian));
        REQUIRE(MatrixAlg::Trace(hermitian) == Complex{5.0, 0.0});

        const Matrix<Real> diagonal{2, 2, {REAL(1e20), REAL(0.0), REAL(0.0), REAL(1.0)}};
        REQUIRE(std::abs(MatrixAlg::Determinant(diagonal) / REAL(1e20) - REAL(1.0)) < REAL(1e-10));
        REQUIRE((diagonal * MatrixAlg::Inverse(diagonal)).IsEqualTo(Matrix<Real>::Identity(2), REAL(1e-10)));
    }

    TEST_CASE("SingleHeader_MatrixAlgSVDAndSubspaces", "[single_header][matrix_analysis][MML2]")
    {
        const Matrix<Real> matrix{2, 3, {
            REAL(1.0), REAL(2.0), REAL(3.0),
            REAL(2.0), REAL(4.0), REAL(6.0)
        }};
        const auto svd = MatrixAlg::SVDDecompose(matrix);
        REQUIRE(svd.U.rows() == 2);
        REQUIRE(svd.U.cols() == 2);
        REQUIRE(svd.V.rows() == 3);
        REQUIRE(svd.V.cols() == 3);
        REQUIRE(svd.rank == 1);

        const Matrix<Real> pseudoinverse = MatrixAlg::PseudoInverse(matrix);
        REQUIRE((matrix * pseudoinverse * matrix).IsEqualTo(matrix, single_header_tolerance()));

        const auto spaces = MatrixAlg::FundamentalSubspacesOf(Matrix<Real>::Identity(3));
        REQUIRE(spaces.nullSpace.rows() == 3);
        REQUIRE(spaces.nullSpace.cols() == 0);
    }

    TEST_CASE("SingleHeader_MatrixAlgComplexSVD", "[single_header][matrix_analysis][complex_svd][MML2]")
    {
        const Matrix<Complex> matrix{2, 2, {
            Complex{2.0, 1.0}, Complex{1.0, -2.0},
            Complex{-1.0, 0.5}, Complex{3.0, 1.0}
        }};
        const auto svd = MatrixAlg::SVDDecompose(matrix);
        Matrix<Complex> sigma(2, 2);
        sigma(0, 0) = Complex{svd.singularValues[0], REAL(0.0)};
        sigma(1, 1) = Complex{svd.singularValues[1], REAL(0.0)};
        Matrix<Complex> adjointV(2, 2);
        for (int row = 0; row < 2; ++row)
            for (int col = 0; col < 2; ++col)
                adjointV(row, col) = std::conj(svd.V(col, row));

        REQUIRE(svd.rank == 2);
        REQUIRE((svd.U * sigma * adjointV).IsEqualTo(matrix, single_header_tolerance()));

        const auto pseudoinverse = MatrixAlg::PseudoInverse(matrix);
        REQUIRE((matrix * pseudoinverse * matrix).IsEqualTo(matrix, single_header_tolerance()));
    }

    TEST_CASE("SingleHeader_HermitianCholeskyAndDefiniteness", "[single_header][matrix_analysis][hermitian_cholesky][MML2]")
    {
        const Matrix<Complex> matrix{2, 2, {
            Complex{4.0, 0.0}, Complex{1.0, 2.0},
            Complex{1.0, -2.0}, Complex{3.0, 0.0}
        }};
        const auto cholesky = MatrixAlg::CholeskyDecompose(matrix);
        Matrix<Complex> adjointLower(2, 2);
        for (int row = 0; row < 2; ++row)
            for (int col = 0; col < 2; ++col)
                adjointLower(row, col) = std::conj(cholesky.L(col, row));

        REQUIRE((cholesky.L * adjointLower).IsEqualTo(matrix, REAL(1e-8)));
        REQUIRE(MatrixAlg::ClassifyDefiniteness(matrix) == MatrixAlg::Definiteness::PositiveDefinite);
        REQUIRE(MatrixAlg::IsPositiveDefinite(matrix));
    }

    TEST_CASE("SingleHeader_GeneralComplexEigenanalysis", "[single_header][matrix_analysis][complex_eigen][MML2]")
    {
        const Matrix<Complex> matrix{2, 2, {
            Complex{1.0, 2.0}, Complex{3.0, -1.0},
            Complex{}, Complex{-2.0, 0.5}
        }};
        const auto result = MatrixAlg::Eigensystem(matrix);

        REQUIRE(result.converged);
        REQUIRE(result.eigenvalues.size() == 2);
        REQUIRE(result.eigenvectors.rows() == 2);
        REQUIRE(result.eigenvectors.cols() == 2);
        REQUIRE(result.maxResidual < single_header_tolerance());
        REQUIRE(std::abs(MatrixAlg::SpectralRadius(matrix) - std::abs(Complex{1.0, 2.0})) < single_header_tolerance());
    }

    TEST_CASE("SingleHeader_MatrixAnalyzer", "[single_header][matrix_analyzer][MML2]")
    {
        MatrixAnalyzer<Real> analyzer(Matrix<Real>{2, 2, {
            REAL(4.0), REAL(1.0),
            REAL(1.0), REAL(3.0)
        }});
        const auto& svd = analyzer.SVDDecompose();
        REQUIRE(&svd == &analyzer.SVDDecompose());
        REQUIRE(analyzer.IsPositiveDefinite());

        const auto analysis = analyzer.Analyze();
        REQUIRE(analysis.isSquare);
        REQUIRE(analysis.rank == 2);
        REQUIRE(analysis.definiteness == MatrixAlg::Definiteness::PositiveDefinite);
    }

    TEST_CASE("SingleHeader_ComplexMatrixAnalyzer", "[single_header][matrix_analyzer][complex][MML2]")
    {
        const Matrix<Complex> matrix{2, 2, {
            Complex{4.0, 0.0}, Complex{1.0, 2.0},
            Complex{1.0, -2.0}, Complex{3.0, 0.0}
        }};
        MatrixAnalyzer<Complex> analyzer(matrix);

        REQUIRE(&analyzer.SVDDecompose() == &analyzer.SVDDecompose());
        REQUIRE(&analyzer.Eigensystem() == &analyzer.Eigensystem());
        REQUIRE(analyzer.Eigensystem().algorithmName == "HermitianJacobi");
        REQUIRE(analyzer.HermitianEigenvalues().size() == 2);
        REQUIRE(analyzer.Analyze().definiteness == MatrixAlg::Definiteness::PositiveDefinite);
    }

    TEST_CASE("SingleHeader_LinearSystemComposedAnalysis", "[single_header][linear_system][MML2]")
    {
        const Matrix<Real> coefficient{2, 2, {
            REAL(1.0), REAL(1.0),
            REAL(2.0), REAL(2.0)
        }};
        Systems::LinearSystem<Real> system(coefficient, Vector<Real>{REAL(1.0), REAL(3.0)});
        const auto analysis = system.Analyze();

        REQUIRE(analysis.matrix.rank == 1);
        REQUIRE(analysis.solutionStatuses == std::vector<Systems::SolutionStatus>{Systems::SolutionStatus::Inconsistent});
        REQUIRE(analysis.recommendedSolver == Systems::LinearSolverRecommendation::SVD);
    }

    TEST_CASE("SingleHeader_SparseSolvers_CG", "[single_header][sparse]")
    {
        auto A = SparseMatrix::laplacian1D<Real>(3);
        std::vector<Real> b{ REAL(1.0), REAL(2.0), REAL(1.0) };
        std::vector<Real> x(3, REAL(0.0));

        SparseSolvers::SolverConfig<Real> config;
        config.setTolerance(REAL(1e-12)).setMaxIterations(50);

        auto result = SparseSolvers::solveCG(A, b, x, config);

        REQUIRE(result.converged());
        REQUIRE(x[0] == Catch::Approx(REAL(2.0)).epsilon(1e-12));
        REQUIRE(x[1] == Catch::Approx(REAL(3.0)).epsilon(1e-12));
        REQUIRE(x[2] == Catch::Approx(REAL(2.0)).epsilon(1e-12));
    }

    TEST_CASE("SingleHeader_AutomaticDifferentiation", "[single_header][AutomaticDifferentiation]")
    {
        auto [value, deriv] = AD::derivative([](auto x) { return x * x * x; }, REAL(2.0));
        REQUIRE(value == Catch::Approx(REAL(8.0)));
        REQUIRE(deriv == Catch::Approx(REAL(12.0)));

        auto grad = AD::gradientReverse([](std::vector<AD::ADVar>& vars) {
            return vars[0] * vars[0] + vars[0] * vars[1] + vars[1] * vars[1];
        }, std::vector<Real>{ REAL(2.0), REAL(3.0) });

        REQUIRE(grad[0] == Catch::Approx(REAL(7.0)));
        REQUIRE(grad[1] == Catch::Approx(REAL(8.0)));

        auto jac = AD::jacobianForwardAD([](const std::vector<AD::Dual<Real>>& vars) {
            return std::vector<AD::Dual<Real>>{ vars[0] * vars[0], vars[0] * vars[1] };
        }, std::vector<Real>{ REAL(3.0), REAL(4.0) });

        REQUIRE(jac.rows() == 2);
        REQUIRE(jac.cols() == 2);
        REQUIRE(jac(0, 0) == Catch::Approx(REAL(6.0)));
        REQUIRE(jac(1, 1) == Catch::Approx(REAL(3.0)));
    }

    TEST_CASE("SingleHeader_Statistics_Descriptive", "[single_header][statistics]")
    {
        Vector<Real> data({ REAL(1.0), REAL(2.0), REAL(3.0), REAL(4.0), REAL(5.0) });
        Vector<Real> weights({ REAL(1.0), REAL(1.0), REAL(2.0), REAL(2.0), REAL(4.0) });

        REQUIRE(Statistics::Mean(data) == Catch::Approx(REAL(3.0)));
        REQUIRE(Statistics::Median(data) == Catch::Approx(REAL(3.0)));
        REQUIRE(Statistics::SampleVariance(data) == Catch::Approx(REAL(2.5)));
        REQUIRE(Statistics::WeightedMean(data, weights) == Catch::Approx(REAL(3.7)));
    }

    TEST_CASE("SingleHeader_Statistics_HistogramAndECDF", "[single_header][statistics]")
    {
        using namespace Statistics::Histogram;

        Vector<Real> sample({ REAL(1.0), REAL(2.0), REAL(2.0), REAL(3.0), REAL(4.0) });

        auto histogram = ComputeHistogram(sample, 3);
        auto table = FrequencyTable(sample);
        auto ecdf = EmpiricalCDF(sample);

        REQUIRE(histogram.totalCount == 5);
        REQUIRE(histogram.numBins == 3);
        REQUIRE(table.totalCount == 5);
        REQUIRE(table.uniqueCount == 4);
        REQUIRE(table.counts.at(REAL(2.0)) == 2);
        REQUIRE(ecdf.n == 5);
        REQUIRE(EvaluateECDF(sample, REAL(2.0)) == Catch::Approx(REAL(3.0) / REAL(5.0)));
        REQUIRE(Quantile(sample, REAL(0.5)) == Catch::Approx(REAL(2.0)));
    }

    TEST_CASE("SingleHeader_Statistics_DiscreteDistributions", "[single_header][statistics]")
    {
        Statistics::BinomialDistribution binomial(10, REAL(0.5));
        Statistics::PoissonDistribution poisson(REAL(3.0));
        Statistics::HypergeometricDistribution hypergeometric(52, 13, 5);

        REQUIRE(binomial.pmf(5) == Catch::Approx(REAL(0.24609375)));
        REQUIRE(binomial.cdf(4) < binomial.cdf(5));
        REQUIRE(poisson.mean() == Catch::Approx(REAL(3.0)));
        REQUIRE(poisson.variance() == Catch::Approx(REAL(3.0)));
        REQUIRE(hypergeometric.pmf(0) > REAL(0.0));
        REQUIRE(hypergeometric.inverseCdf(REAL(0.95)) >= 0);
    }
}
