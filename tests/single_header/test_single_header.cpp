// Test that MML single headers work standalone.
// This file explicitly includes ONLY the selected single header - no other MML headers.

#include <cmath>

#ifndef __STDCPP_WANT_MATH_SPEC_FUNCS__
#define MML_TEST_EXPECT_NO_MATH_OPT_IN
#endif

// Force single header mode and include the selected variant directly.
#ifndef MML_SINGLE_HEADER_INCLUDE
#define MML_SINGLE_HEADER_INCLUDE <MML.h>
#endif
#define MML_USE_SINGLE_HEADER
#include MML_SINGLE_HEADER_INCLUDE

namespace MML::FunctionSpaces {
	class ChebyshevCollocationSpace1D;
	class TrialSpace1D;
	class FunctionSpace1D;
}

template<class Type>
concept CompleteFunctionSpaceType = requires { sizeof(Type); };

static_assert(!CompleteFunctionSpaceType<MML::FunctionSpaces::ChebyshevCollocationSpace1D>);
static_assert(!CompleteFunctionSpaceType<MML::FunctionSpaces::TrialSpace1D>);
static_assert(!CompleteFunctionSpaceType<MML::FunctionSpaces::FunctionSpace1D>);

namespace MML::Curves {
	class LemniscateCurve;
	class DeltoidCurve;
	class AstroidCurve;
	class EpitrochoidCurve;
	class ButterflyCurve;
	class TwistedCubicCurve;
	class ToroidalSpiralCurve;
}
namespace MML::Surfaces {
	class MonkeySaddle;
	class MobiusStrip;
	class Helicoid;
	class EnneperSurface;
	class KleinBottle;
	class DiniSurface;
}

template<class Type>
concept CompleteShowcaseType = requires { sizeof(Type); };

static_assert(!CompleteShowcaseType<MML::Curves::LemniscateCurve>);
static_assert(!CompleteShowcaseType<MML::Curves::DeltoidCurve>);
static_assert(!CompleteShowcaseType<MML::Curves::AstroidCurve>);
static_assert(!CompleteShowcaseType<MML::Curves::EpitrochoidCurve>);
static_assert(!CompleteShowcaseType<MML::Curves::ButterflyCurve>);
static_assert(!CompleteShowcaseType<MML::Curves::TwistedCubicCurve>);
static_assert(!CompleteShowcaseType<MML::Curves::ToroidalSpiralCurve>);
static_assert(!CompleteShowcaseType<MML::Surfaces::MonkeySaddle>);
static_assert(!CompleteShowcaseType<MML::Surfaces::MobiusStrip>);
static_assert(!CompleteShowcaseType<MML::Surfaces::Helicoid>);
static_assert(!CompleteShowcaseType<MML::Surfaces::EnneperSurface>);
static_assert(!CompleteShowcaseType<MML::Surfaces::KleinBottle>);
static_assert(!CompleteShowcaseType<MML::Surfaces::DiniSurface>);

#if defined(MML_TEST_EXPECT_NO_MATH_OPT_IN) && defined(__STDCPP_WANT_MATH_SPEC_FUNCS__)
#error MML must not inject a special-math opt-in after cmath has already been included.
#endif

#include <catch2/catch_all.hpp>
#include <type_traits>

#if defined(MML_SINGLE_HEADER_EXPECT_FLOAT)
static_assert(std::is_same_v<MML::Real, float>);
#elif defined(MML_SINGLE_HEADER_EXPECT_LONG_DOUBLE)
static_assert(std::is_same_v<MML::Real, long double>);
#else
static_assert(std::is_same_v<MML::Real, double>);
#endif

using namespace MML;

TEST_CASE("Single header rank significance uses scalar-typed critical values", "[single-header][statistics]")
{
	Statistics::RankCorrelationDetailedResult result;
	result.zScore = 2;
	CHECK(result.IsSignificant());
	CHECK_FALSE(result.IsSignificant(Real(0.01)));
	result.zScore = -2;
	CHECK(result.IsSignificant());
	CHECK_FALSE(result.IsSignificant(Real(0.01)));
}

TEST_CASE("Single header retains plane coefficients and protected Dirac parameter after member cleanup",
          "[single-header][reserved-members]")
{
	Plane3D plane(Pnt3Cart(0, 0, 2), Vec3Cart(0, 0, 1));
	REQUIRE(plane.C() == 1);
	REQUIRE(plane.D() == -2);
	plane.D() = -3;
	REQUIRE(plane.D() == -3);
	class CustomDelta : public DiracFunction {
	public:
		using DiracFunction::DiracFunction;
		Real operator()(Real) const override { return static_cast<Real>(_sharpness); }
	};
	REQUIRE(CustomDelta(7)(0) == 7);
}

template<class Type>
concept SingleHeaderContainerInverse = requires(const Type& matrix) { matrix.Inverse(); };
static_assert(!SingleHeaderContainerInverse<Matrix<Real>>);
static_assert(!SingleHeaderContainerInverse<MatrixNM<Real, 2, 2>>);
static_assert(!SingleHeaderContainerInverse<MatrixSym<Real>>);
static_assert(!SingleHeaderContainerInverse<TridiagonalMatrix<Real>>);

TEST_CASE("Single header keeps inversion algorithm-only and prints matrix options consistently",
          "[single-header][MatrixInversion][PrintFormatting]")
{
	Matrix<Real> matrix(2, 2, {4, 7, 2, 6});
	const auto inverse = MatrixAlg::Inverse(matrix);
	REQUIRE((matrix * inverse).IsIdentity(Precision::LooseTolerance));
	MatrixNM<Real, 2, 2> fixed{2, 0, 0, 4};
	const auto metric = Metric<2, Cartesian2>(fixed);
	REQUIRE(IsNearlyEqual(metric.contravariantComponents()(1, 1), Real(1) / 4));
	MatrixPrintFormat format = MatrixPrintFormat::Compact();
	format.width = 0;
	std::ostringstream output;
	matrix.Print(output, format);
	REQUIRE(output.str() == "[4.00, 7.00; 2.00, 6.00]\n");
}

TEST_CASE("Single header includes each centralized exception exactly once", "[single-header][exceptions]")
{
	STATIC_REQUIRE(std::is_base_of_v<std::invalid_argument, FourierError>);
	STATIC_REQUIRE(std::is_base_of_v<std::runtime_error, StatisticsError>);
	STATIC_REQUIRE(std::is_base_of_v<std::invalid_argument, CurveFittingError>);
	const FourierError fourier("fourier");
	const StatisticsError statistics("statistics");
	const CurveFittingError fitting("fitting");
	REQUIRE(std::string(fourier.message()) == "fourier");
	REQUIRE(std::string(statistics.message()) == "statistics");
	REQUIRE(std::string(fitting.message()) == "fitting");
}

TEST_CASE("Single header supports unsigned vector printing and norms", "[single-header][UnsignedVector]")
{
	const auto maximum = std::numeric_limits<unsigned long long>::max();
	Vector<unsigned long long> values{maximum, 0};
	std::ostringstream out;
	values.Print(out, 0, 6, Real(2));
	REQUIRE(out.str() == "[" + std::to_string(maximum) + ", 0]");
	REQUIRE(Abs(maximum) == static_cast<Real>(maximum));
	Vector<unsigned int> small{3, 4, 0};
	REQUIRE(IsNearlyEqual(small.NormL2(), Real(5)));
	VectorN<unsigned int, 3> fixed{3, 4, 0};
	REQUIRE(IsNearlyEqual(fixed.NormL1(), Real(7)));
}

TEST_CASE("Single header vector persistence has exact text and transactional failure", "[single-header][VectorPersistence]")
{
	Vector<double> original{std::nextafter(1.0, 2.0), -0.0};
	std::stringstream text;
	REQUIRE(Persistence::Save(text, original, Persistence::Format::Text).success);
	Vector<double> restored;
	REQUIRE(Persistence::Load(text, restored, Persistence::Format::Text).success);
	REQUIRE(restored[0] == original[0]);
	REQUIRE(std::signbit(restored[1]));
	std::istringstream invalid("2\n1\nbad\n");
	REQUIRE_FALSE(Persistence::Load(invalid, restored, Persistence::Format::Text).success);
	REQUIRE(restored[0] == original[0]);
	std::stringstream integers;
	REQUIRE(Persistence::Save(integers, Vector<int>{1, -2}, Persistence::Format::Text).success);
	Vector<int> decoded;
	REQUIRE(Persistence::Load(integers, decoded, Persistence::Format::Text).success);
	REQUIRE(decoded == Vector<int>{1, -2});
}

TEST_CASE("Single header separates numerical analyzers from optional report rendering", "[single-header][AnalysisReports]")
{
	RealFunction function(+[](Real x) { return x * x; });
	const RealFunctionAnalyzer analyzer(function);
	const auto interval = analyzer.AnalyzeInterval(Real(-1), Real(1), 16);
	const auto continuity = analyzer.AnalyzeContinuity(Real(-1), Real(1), 16);
	REQUIRE(interval.minValue == 0);
	std::ostringstream reports;
	AnalysisReports::PrintPointAnalysis(reports, analyzer.AnalyzePoint(Real(0)));
	AnalysisReports::PrintIntervalAnalysis(reports, interval);
	AnalysisReports::PrintDetailedIntervalAnalysis(reports, analyzer.AnalyzePoints(Real(-1), Real(1), 16));
	AnalysisReports::PrintContinuityAnalysis(reports, continuity);
	REQUIRE(reports.str().find("Sampled minimum: 0") != std::string::npos);
	REQUIRE(reports.str().find("No discontinuity candidates") != std::string::npos);
	REQUIRE_FALSE(analyzer.IsMonotonicOnSamples(Real(-1), Real(1), 16));
	REQUIRE(analyzer.AppearsContinuous(Real(-1), Real(1), 16));
	const auto matrix = MatrixAnalyzer<Real>(Matrix<Real>{2, 2, {2, 0, 0, 3}}).Analyze();
	AnalysisReports::PrintMatrixAnalysis(reports, matrix);
	REQUIRE(reports.str().find("Rank: 2") != std::string::npos);
}

TEST_CASE("Single header persistence uses exact automatic text precision and canonical widths", "[single-header][Persistence]")
{
	Matrix<double> original(1, 1);
	original(0, 0) = std::nextafter(1.0, 2.0);
	std::stringstream text;
	text << std::fixed << std::setprecision(2);
	REQUIRE(Persistence::SaveOptions{}.precision == 0);
	REQUIRE(Persistence::Save(text, original, Persistence::Format::Text).success);
	REQUIRE(text.precision() == 2);
	Matrix<double> decoded;
	REQUIRE(Persistence::Load(text, decoded, Persistence::Format::Text).success);
	REQUIRE(decoded(0, 0) == original(0, 0));
	REQUIRE(Persistence::ExpectedScalarByteSize(Persistence::BinaryScalarType::LongDouble) == 24);
}

TEST_CASE("Single header supports explicit visualizer bundle configuration", "[single-header][Visualizer]")
{
	struct RestoreRoot {
		std::filesystem::path previous = detail::GetVisualizerRootState();
		~RestoreRoot() { SetVisualizerRoot(previous); }
	} restore;
	const auto root = std::filesystem::temp_directory_path() / "mml standalone bundle";
	SetVisualizerRoot(root);
	REQUIRE(std::filesystem::path(GetVisualizerRoot()) == root);
	REQUIRE(GetVisualizerPath("MML_RealFunctionVisualizer", VisualizerBackend::Qt).find(root.string()) == 0);
	REQUIRE_THROWS_AS(SetVisualizerRoot(std::string("bad\0root", 8)), VisualizerError);
	REQUIRE(std::filesystem::path(GetVisualizerRoot()) == root);
	REQUIRE(GetRealFuncVisualizerPath() ==
		detail::visualizer::ExecutablePath(detail::visualizer::Kind::RealFunction));
	REQUIRE(detail::visualizer::ExecutableNames.size() == 12);
	const auto empty = Visualizer::VisualizeRealFunctionFromFile("");
	REQUIRE_FALSE(empty.success);
	REQUIRE(empty.error == IOError::INVALID_PARAMETERS);
	const auto invalidTimeout = detail::visualizer::ExecuteVisualizer(GetRealFuncVisualizerPath(), {}, -2);
	REQUIRE_FALSE(invalidTimeout.success);
	REQUIRE(invalidTimeout.error == IOError::INVALID_PARAMETERS);
}

namespace MML::Tests::SingleHeader
{
    using namespace MML::Systems;

    constexpr Real single_header_tolerance()
    {
        if constexpr (std::is_same_v<Real, float>)
            return static_cast<MML::Real>(1e-4);
        else
            return static_cast<MML::Real>(1e-8);
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
        Real getMinT() const override { return static_cast<MML::Real>(-10.0); }
        Real getMaxT() const override { return static_cast<MML::Real>(10.0); }
    };

    class SingleHeaderSurface : public IParametricSurfaceRect<2>
    {
    public:
        VectorN<Real, 2> operator()(Real u, Real w) const override { return {u * u, u + w}; }
        Real getMinU() const override { return static_cast<MML::Real>(-10.0); }
        Real getMaxU() const override { return static_cast<MML::Real>(10.0); }
        Real getMinW() const override { return static_cast<MML::Real>(-10.0); }
        Real getMaxW() const override { return static_cast<MML::Real>(10.0); }
    };

    class SingleHeaderTensorField : public ITensorField2<2>
    {
    public:
        SingleHeaderTensorField() : ITensorField2<2>(0, 2) {}
        Real component(int, int, const VectorN<Real, 2>& point) const override { return point[0] * point[0] + point[1]; }
        Tensor2<2> operator()(const VectorN<Real, 2>&) const override { return Tensor2<2>(2, 0); }
    };

    TEST_CASE("SingleHeader_FirstDerivativeAdapters", "[single_header][derivation]")
    {
        const Real h = static_cast<MML::Real>(1e-3);
        const VectorN<Real, 2> point{static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0)};
        SingleHeaderRealFunction real_function;
        SingleHeaderComplexFunction complex_function;
        SingleHeaderScalarFunction scalar_function;
        SingleHeaderVectorFunction vector_function;
        SingleHeaderCurve curve;
        SingleHeaderSurface surface;
        SingleHeaderTensorField tensor_field;

        const auto expected_derivative = Catch::Approx(static_cast<MML::Real>(4.0)).epsilon(static_cast<double>(single_header_tolerance()));
        REQUIRE(Derivation::NDer4(real_function, static_cast<MML::Real>(2.0), h) == expected_derivative);
        REQUIRE(Derivation::NDer4Complex(complex_function, Complex(static_cast<MML::Real>(2.0)), h).real() == expected_derivative);
        REQUIRE(Derivation::NDer4Partial(scalar_function, 0, point, h) == expected_derivative);
        REQUIRE(Derivation::NDer4Partial(vector_function, 0, 0, point, h) == expected_derivative);
        REQUIRE(Derivation::NDer4(curve, static_cast<MML::Real>(2.0), h)[0] == expected_derivative);
        REQUIRE(Derivation::NDer2_u(surface, static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0), h)[0] == expected_derivative);
        REQUIRE(Derivation::NDer4Partial(tensor_field, 0, 0, 0, point, h) == expected_derivative);
    }

    TEST_CASE("SingleHeader_RODAS4", "[single_header][RODAS4]")
    {
        class DecayDAE : public IODESystemDAEWithJacobian {
        public:
            int getDiffDim() const override { return 1; }
            int getAlgDim() const override { return 1; }
            void diffEqs(Real, const Vector<Real>& x, const Vector<Real>&, Vector<Real>& f) const override { f[0] = -x[0]; }
            void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override { g[0] = y[0] - x[0]; }
            void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = -1; }
            void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 0; }
            void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = -1; }
            void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& J) const override { J(0, 0) = 1; }
        } system;
        DAESolverConfig config;
        config.abs_tolerance = config.rel_tolerance = single_header_tolerance();
        auto result = SolveDAERODAS4(system, 0, Vector<Real>{1}, Vector<Real>{1}, 1, config);
        INFO(result.error_message);
        REQUIRE(result.status == AlgorithmStatus::Success);
        REQUIRE(result.solution.getTValue(result.solution.size() - 1) == 1);
        REQUIRE(IsWithinAbs(result.solution.getXValuesAtEnd()[0], std::exp(Real(-1)),
            single_header_tolerance() * 10));
        REQUIRE(result.linear_solves == 6 * result.lu_factorizations);
        REQUIRE(result.newton_iterations == 0);
    }















    TEST_CASE("SingleHeader_Algebra_LieGroups", "[single_header][Algebra]")
    {
        const Algebra::SO2 planar(Constants::PI / static_cast<MML::Real>(2.0));
        const auto spatial = Algebra::SO3::Exp({static_cast<MML::Real>(0.2), static_cast<MML::Real>(-0.1), static_cast<MML::Real>(0.3)});
        MatrixNM<Real, 3, 3> drifted = spatial.matrix();
        drifted(0, 0) += static_cast<MML::Real>(1e-5);

        REQUIRE(planar.apply({1, 0}).IsEqualTo({0, 1}, single_header_tolerance()));
        REQUIRE(spatial.log().IsEqualTo({static_cast<MML::Real>(0.2), static_cast<MML::Real>(-0.1), static_cast<MML::Real>(0.3)}, single_header_tolerance()));
        REQUIRE(Algebra::SO3::IsRotationMatrix(spatial.matrix()));
        REQUIRE(Algebra::SO3::IsRotationMatrix(Algebra::SO3::Project(drifted).matrix()));
    }

    TEST_CASE("SingleHeader_Algebra_RigidMotions", "[single_header][Algebra]")
    {
        const Algebra::SE3 transform(
            Algebra::SO3::FromAxisAngle({0, 0, 1}, Constants::PI / static_cast<MML::Real>(2.0)), {2, 3, 4});
        const auto restored = Algebra::SE3::FromHomogeneousMatrix(transform.homogeneous_matrix());

        REQUIRE(transform.apply_point({1, 0, 0}).IsEqualTo({2, 4, 4}, single_header_tolerance()));
        REQUIRE(transform.apply_vector({1, 0, 0}).IsEqualTo({0, 1, 0}, single_header_tolerance()));
        REQUIRE(restored.homogeneous_matrix().IsEqualTo(transform.homogeneous_matrix(), single_header_tolerance()));
        REQUIRE(transform.inverse().apply_point(transform.apply_point({1, 2, 3})).IsEqualTo({1, 2, 3}, single_header_tolerance()));
    }

    TEST_CASE("SingleHeader_Vector_Basic", "[single_header]")
    {
        Vector<Real> v1({ static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0) });
        Vector<Real> v2({ static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0), static_cast<MML::Real>(6.0) });
        
        // Vector addition
        auto v3 = v1 + v2;
        REQUIRE(v3.size() == 3);
        REQUIRE(v3[0] == static_cast<MML::Real>(5.0));
        REQUIRE(v3[1] == static_cast<MML::Real>(7.0));
        REQUIRE(v3[2] == static_cast<MML::Real>(9.0));
        
        // Scalar product (free function in VectorUtils)
        Real dot = Utils::ScalarProduct(v1, v2);
        REQUIRE(dot == static_cast<MML::Real>(32.0)); // 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
    }
    
    TEST_CASE("SingleHeader_Matrix_Basic", "[single_header]")
    {
        Matrix<Real> A(2, 3, {
            static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0),
            static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0), static_cast<MML::Real>(6.0)
        });
        
        Matrix<Real> B(3, 2, {
            static_cast<MML::Real>(7.0), static_cast<MML::Real>(8.0),
            static_cast<MML::Real>(9.0), static_cast<MML::Real>(10.0),
            static_cast<MML::Real>(11.0), static_cast<MML::Real>(12.0)
        });
        
        // Matrix multiplication: (2x3) * (3x2) = (2x2)
        auto C = A * B;
        
        REQUIRE(C.rows() == 2);
        REQUIRE(C.cols() == 2);
        
        // C[0][0] = 1*7 + 2*9 + 3*11 = 7 + 18 + 33 = 58
        REQUIRE(C[0][0] == static_cast<MML::Real>(58.0));
        // C[0][1] = 1*8 + 2*10 + 3*12 = 8 + 20 + 36 = 64
        REQUIRE(C[0][1] == static_cast<MML::Real>(64.0));
        // C[1][0] = 4*7 + 5*9 + 6*11 = 28 + 45 + 66 = 139
        REQUIRE(C[1][0] == static_cast<MML::Real>(139.0));
        // C[1][1] = 4*8 + 5*10 + 6*12 = 32 + 50 + 72 = 154
        REQUIRE(C[1][1] == static_cast<MML::Real>(154.0));
    }
    
    TEST_CASE("SingleHeader_MatrixVector_Multiply", "[single_header]")
    {
        Matrix<Real> A(2, 3, {
            static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0),
            static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0), static_cast<MML::Real>(6.0)
        });
        
        Vector<Real> v({ static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0) });
        
        // Matrix-vector multiplication: (2x3) * (3x1) = (2x1)
        auto result = A * v;
        
        REQUIRE(result.size() == 2);
        // result[0] = 1*1 + 2*2 + 3*3 = 1 + 4 + 9 = 14
        REQUIRE(result[0] == static_cast<MML::Real>(14.0));
        // result[1] = 4*1 + 5*2 + 6*3 = 4 + 10 + 18 = 32
        REQUIRE(result[1] == static_cast<MML::Real>(32.0));
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

        const Matrix<Real> diagonal{2, 2, {static_cast<MML::Real>(1e20), static_cast<MML::Real>(0.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(1.0)}};
        REQUIRE(std::abs(MatrixAlg::Determinant(diagonal) / static_cast<MML::Real>(1e20) - static_cast<MML::Real>(1.0)) < static_cast<MML::Real>(1e-10));
        REQUIRE((diagonal * MatrixAlg::Inverse(diagonal)).IsEqualTo(Matrix<Real>::Identity(2), static_cast<MML::Real>(1e-10)));
    }

    TEST_CASE("SingleHeader_MatrixAlgSVDAndSubspaces", "[single_header][matrix_analysis][MML2]")
    {
        const Matrix<Real> matrix{2, 3, {
            static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0),
            static_cast<MML::Real>(2.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(6.0)
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
        sigma(0, 0) = Complex{svd.singularValues[0], static_cast<MML::Real>(0.0)};
        sigma(1, 1) = Complex{svd.singularValues[1], static_cast<MML::Real>(0.0)};
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

        REQUIRE((cholesky.L * adjointLower).IsEqualTo(matrix, static_cast<MML::Real>(1e-8)));
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
            static_cast<MML::Real>(4.0), static_cast<MML::Real>(1.0),
            static_cast<MML::Real>(1.0), static_cast<MML::Real>(3.0)
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
            static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0),
            static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0)
        }};
        Systems::LinearSystem<Real> system(coefficient, Vector<Real>{static_cast<MML::Real>(1.0), static_cast<MML::Real>(3.0)});
        const auto analysis = system.Analyze();

        REQUIRE(analysis.matrix.rank == 1);
        REQUIRE(analysis.solutionStatuses == std::vector<Systems::SolutionStatus>{Systems::SolutionStatus::Inconsistent});
        REQUIRE(analysis.recommendedSolver == Systems::LinearSolverRecommendation::SVD);
    }

    TEST_CASE("SingleHeader_SparseSolvers_CG", "[single_header][sparse]")
    {
        auto A = SparseMatrix::Laplacian1D<Real>(3);
        Vector<Real> b{ static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0) };
        Vector<Real> x(3, static_cast<MML::Real>(0.0));

        SparseSolvers::SolverConfig<Real> config;
        config.setTolerance(static_cast<MML::Real>(1e-12)).setMaxIterations(50);

        auto result = SparseSolvers::solveCG(A, b, x, config);

        REQUIRE(result.converged);
        REQUIRE(x[0] == Catch::Approx(static_cast<MML::Real>(2.0)).epsilon(1e-12));
        REQUIRE(x[1] == Catch::Approx(static_cast<MML::Real>(3.0)).epsilon(1e-12));
        REQUIRE(x[2] == Catch::Approx(static_cast<MML::Real>(2.0)).epsilon(1e-12));
    }

    TEST_CASE("SingleHeader_AutomaticDifferentiation", "[single_header][AutomaticDifferentiation]")
    {
        auto [value, deriv] = AD::CalcDerivativeForwardAD([](auto x) { return x * x * x; }, static_cast<MML::Real>(2.0));
        REQUIRE(value == Catch::Approx(static_cast<MML::Real>(8.0)));
        REQUIRE(deriv == Catch::Approx(static_cast<MML::Real>(12.0)));

        auto grad = AD::CalcGradientReverseAD([](std::vector<AD::ADVar>& vars) {
            return vars[0] * vars[0] + vars[0] * vars[1] + vars[1] * vars[1];
        }, std::vector<Real>{ static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0) });

        REQUIRE(grad[0] == Catch::Approx(static_cast<MML::Real>(7.0)));
        REQUIRE(grad[1] == Catch::Approx(static_cast<MML::Real>(8.0)));

        auto jac = AD::CalcJacobianForwardAD([](const std::vector<AD::Dual<Real>>& vars) {
            return std::vector<AD::Dual<Real>>{ vars[0] * vars[0], vars[0] * vars[1] };
        }, std::vector<Real>{ static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0) });

        REQUIRE(jac.rows() == 2);
        REQUIRE(jac.cols() == 2);
        REQUIRE(jac(0, 0) == Catch::Approx(static_cast<MML::Real>(6.0)));
        REQUIRE(jac(1, 1) == Catch::Approx(static_cast<MML::Real>(3.0)));
    }

    TEST_CASE("SingleHeader_Statistics_Descriptive", "[single_header][statistics]")
    {
        Vector<Real> data({ static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0), static_cast<MML::Real>(5.0) });
        Vector<Real> weights({ static_cast<MML::Real>(1.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(4.0) });

        REQUIRE(Statistics::Mean(data) == Catch::Approx(static_cast<MML::Real>(3.0)));
        REQUIRE(Statistics::Median(data) == Catch::Approx(static_cast<MML::Real>(3.0)));
        REQUIRE(Statistics::SampleVariance(data) == Catch::Approx(static_cast<MML::Real>(2.5)));
        REQUIRE(Statistics::WeightedMean(data, weights) == Catch::Approx(static_cast<MML::Real>(3.7)));
    }

    TEST_CASE("SingleHeader_Statistics_HistogramAndECDF", "[single_header][statistics]")
    {
        using namespace Statistics::Histogram;

        Vector<Real> sample({ static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(2.0), static_cast<MML::Real>(3.0), static_cast<MML::Real>(4.0) });

        auto histogram = ComputeHistogram(sample, 3);
        auto table = FrequencyTable(sample);
        auto ecdf = EmpiricalCDF(sample);

        REQUIRE(histogram.totalCount == 5);
        REQUIRE(histogram.numBins == 3);
        REQUIRE(table.totalCount == 5);
        REQUIRE(table.uniqueCount == 4);
        REQUIRE(table.counts.at(static_cast<MML::Real>(2.0)) == 2);
        REQUIRE(ecdf.n == 5);
        REQUIRE(EvaluateECDF(sample, static_cast<MML::Real>(2.0)) == Catch::Approx(static_cast<MML::Real>(3.0) / static_cast<MML::Real>(5.0)));
        REQUIRE(Quantile(sample, static_cast<MML::Real>(0.5)) == Catch::Approx(static_cast<MML::Real>(2.0)));
    }

    TEST_CASE("SingleHeader_UnifiedResults", "[single_header][result-policy]")
    {
        STATIC_REQUIRE(std::is_base_of_v<IterativeResultBase, Optimization::LinearProgrammingResult>);
        STATIC_REQUIRE(std::is_base_of_v<IterativeResultBase, SparseSolvers::KrylovResult<Real>>);
        STATIC_REQUIRE_FALSE(std::is_convertible_v<IOResult, bool>);
        std::stringstream stream;
        Vector<Real> values{static_cast<MML::Real>(1.0), static_cast<MML::Real>(2.0)};
        REQUIRE(Persistence::SaveJson(stream, values));
        Vector<Real> restored;
        REQUIRE(Persistence::LoadJson(stream, restored));
        REQUIRE(restored == values);
        const auto malformed = Persistence::ParseJson("{\"x\":}");
        REQUIRE_FALSE(malformed.result);
        REQUIRE(malformed.result.offset == 5);
        REQUIRE(malformed.result.offset_kind == IOOffsetKind::Byte);
    }

		TEST_CASE("SingleHeader_FieldComponentConversions", "[single_header][fields][wewi.8.1]") {
			MetricTensorCylindrical metric;
			const Vec3Cyl point{2, Real(0.5), 1};
			const VectorN<Real, 3> coordinate{1, 3, -2};
			const auto physical = VectorFieldOperations::ToPhysical(coordinate, point, metric);
			REQUIRE(physical == VectorN<Real, 3>{1, 6, -2});
			REQUIRE(VectorFieldOperations::ToContravariant(physical, point, metric).IsEqualTo(coordinate));
			REQUIRE_THROWS_AS(VectorFieldOperations::ToPhysical(coordinate, Vec3Cyl{0, 0, 1}, metric), DomainError);
		}

		TEST_CASE("SingleHeader_AnalyticalChristoffels", "[single_header][metrics][wewi.8.2]") {
			MetricTensorSpherical spherical;
			const MetricTensorField<3>& base = spherical;
			const VectorN<Real, 3> pos{2, Real(0.7), Real(0.3)};
			REQUIRE(base.GetChristoffelSymbolSecondKind(0, 1, 1, pos) == -pos[0]);
			REQUIRE(base.GetChristoffelSymbolSecondKind(1, 0, 1, pos) == 1 / pos[0]);
			REQUIRE_THAT(base.GetChristoffelSymbolSecondKindNumerical(0, 1, 1, pos),
									 Catch::Matchers::WithinAbs(-pos[0], single_header_tolerance()));
			REQUIRE_THROWS_AS(base.GetChristoffelSymbolSecondKind(0, 0, 0, VectorN<Real, 3>{0, 1, 0}), DomainError);
			MetricTensorSphericalContravar inverse;
			REQUIRE_THROWS_AS(inverse.GetChristoffelSymbolSecondKind(0, 1, 1, pos), DomainError);
		}

		TEST_CASE("SingleHeader_MatrixSingularityPolicy", "[single_header][metrics][wewi.8.3]") {
			const VectorN<Real, 3> origin{0, 1, 0};
			MetricTensorSpherical metric;
			REQUIRE_THROWS_AS(metric.GetContravariantMetric(origin), SingularMatrixError);
			REQUIRE(std::isnan(metric.GetContravariantMetric(origin, SingularityPolicy::ReturnNaN)(0, 0)));
			REQUIRE(std::isnan(CoordTransfSpherToCart.inverseJacobian(VectorN<Real, 3>{0, 0, 0}, SingularityPolicy::ReturnNaN)(0, 0)));
			VectorFunction<3> field([](const VectorN<Real, 3>& p) { return p; });
			FieldOperationConfig config;
			config.exception_policy = EvaluationExceptionPolicy::ConvertToStatus;
			const auto result = VectorFieldOperations::DivergenceViaDetDetailed(field, origin, metric, config);
			REQUIRE_FALSE(result.IsSuccess());
			REQUIRE(result.status == AlgorithmStatus::InvalidInput);
			REQUIRE_FALSE(result.error_message.empty());
		}

		TEST_CASE("SingleHeader_ExceptionClassification", "[single_header][exceptions]")
    {
        Vector<Real> values(2);
        REQUIRE_THROWS_AS(values.at(-1), IndexError);
        const VectorN<Real, 2> fixed;
        REQUIRE_THROWS_AS(fixed.at(2), std::out_of_range);
        Matrix<Real> matrix(2, 2, {static_cast<MML::Real>(2.0), static_cast<MML::Real>(1.0), static_cast<MML::Real>(0.0), static_cast<MML::Real>(2.0)});
        REQUIRE_THROWS_AS(CholeskySolver<Real>(matrix), DomainError);
        STATIC_REQUIRE_FALSE(std::is_convertible_v<std::string, DataError>);
    }

		TEST_CASE("SingleHeader_ADGeometryAndFields", "[single_header][autodiff][wewi.8.4]") {
			auto curve = AD::MakeADCurve([](auto t) {
				using std::sin;
				using std::cos;
				using S = decltype(t);
				return std::array<S, 3>{2 * cos(t), 2 * sin(t), t};
			});
			REQUIRE_THAT(curve.getCurvature(Real(0.7)), Catch::Matchers::WithinAbs(Real(0.4), single_header_tolerance()));
			REQUIRE_THAT(curve.getTorsion(Real(0.7)), Catch::Matchers::WithinAbs(Real(0.2), single_header_tolerance()));
			auto surface = AD::MakeADSurface([](auto u, auto v) {
				using std::sin;
				using std::cos;
				using S = decltype(u);
				return std::array<S, 3>{sin(u) * cos(v), sin(u) * sin(v), cos(u)};
			});
			REQUIRE_THAT(surface.GaussianCurvature(Real(0.8), Real(0.4)), Catch::Matchers::WithinAbs(Real(1), single_header_tolerance()));
			REQUIRE_THAT(DifferentialGeometry::GaussianCurvatureIntrinsic(surface, Real(0.8), Real(0.4)),
									 Catch::Matchers::WithinAbs(Real(1), single_header_tolerance()));
			auto potential = AD::MakeADScalarField<2>([](const auto& p) { return p[0] * p[0] + p[1] * p[1]; });
			REQUIRE(ScalarFieldOperations::GradientCart(potential, VectorN<Real, 2>{1, 2})[1] == 4);
			REQUIRE(ScalarFieldOperations::LaplacianCart(potential, VectorN<Real, 2>{1, 2}) == 4);
			REQUIRE(potential.Hessian({1, 2})(0, 0) == 2);
			REQUIRE(std::isnan(ScalarFieldOperations::GradientCartDetailed(potential, VectorN<Real, 2>{1, 2}).error[0]));
		}

		TEST_CASE("SingleHeader_Statistics_DiscreteDistributions", "[single_header][statistics]")
    {
        Statistics::BinomialDistribution binomial(10, static_cast<MML::Real>(0.5));
        Statistics::PoissonDistribution poisson(static_cast<MML::Real>(3.0));
        Statistics::HypergeometricDistribution hypergeometric(52, 13, 5);

        REQUIRE(binomial.pmf(5) == Catch::Approx(static_cast<MML::Real>(0.24609375)));
        REQUIRE(binomial.cdf(4) < binomial.cdf(5));
        REQUIRE(poisson.mean() == Catch::Approx(static_cast<MML::Real>(3.0)));
        REQUIRE(poisson.variance() == Catch::Approx(static_cast<MML::Real>(3.0)));
        REQUIRE(hypergeometric.pmf(0) > static_cast<MML::Real>(0.0));
        REQUIRE(hypergeometric.inverseCdf(static_cast<MML::Real>(0.95)) >= 0);
    }
}
