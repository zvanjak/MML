#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Function.h"
#include "core/MetricTensor.h"
#include "core/CoordTransf.h"
#include "core/CoordTransf/CoordTransfSpherical.h"
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////
///                    PREDEFINED METRIC TENSORS
///////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_MetricTensorCartesian()
{
    std::cout << "=== METRIC TENSOR: CARTESIAN (EUCLIDEAN) ===" << std::endl;
    
    MetricTensorCartesian3D metricCart;
    
    Vector3Cartesian pos{1.0, 2.0, 3.0};
    
    // Get covariant metric (gᵢⱼ)
    auto g = metricCart.GetCovariantMetric(pos);
    
    std::cout << "\nCovariant metric g_ij at (" << pos[0] << ", " << pos[1] << ", " << pos[2] << "):" << std::endl;
    for (int i = 0; i < 3; i++) {
        std::cout << "  [";
        for (int j = 0; j < 3; j++)
            std::cout << " " << g[i][j];
        std::cout << " ]" << std::endl;
    }
    std::cout << "Expected: Identity matrix (Euclidean space)" << std::endl;
    
    // Get contravariant metric (gⁱʲ)
    auto g_inv = metricCart.GetContravariantMetric(pos);
    
    std::cout << "\nContravariant metric g^ij:" << std::endl;
    for (int i = 0; i < 3; i++) {
        std::cout << "  [";
        for (int j = 0; j < 3; j++)
            std::cout << " " << g_inv[i][j];
        std::cout << " ]" << std::endl;
    }
    std::cout << "Expected: Identity matrix (inverse of identity)" << std::endl;
    
    // Compute distance using metric
    Vector3Cartesian dr{0.1, 0.2, 0.3};
    Real ds_squared = 0.0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            ds_squared += g[i][j] * dr[i] * dr[j];
    
    std::cout << "\nLine element ds² for dr = (0.1, 0.2, 0.3):" << std::endl;
    std::cout << "  ds² = " << ds_squared << std::endl;
    std::cout << "  Expected: 0.01 + 0.04 + 0.09 = 0.14" << std::endl;
    
    // Christoffel symbols (all zero in Cartesian)
    std::cout << "\nChristoffel symbols (all should be ~0 in flat space):" << std::endl;
    Real Gamma_0_11 = metricCart.GetChristoffelSymbolSecondKind(0, 1, 1, pos);
    Real Gamma_1_00 = metricCart.GetChristoffelSymbolSecondKind(1, 0, 0, pos);
    std::cout << "  Γ⁰₁₁ = " << Gamma_0_11 << std::endl;
    std::cout << "  Γ¹₀₀ = " << Gamma_1_00 << std::endl;
}

void Docs_Demo_MetricTensorSpherical()
{
    std::cout << "\n=== METRIC TENSOR: SPHERICAL ===" << std::endl;
    
    MetricTensorSpherical metricSpher;
    
    // Position: r=5, θ=π/4, φ=π/3
    Vector3Spherical pos{5.0, Constants::PI/4, Constants::PI/3};
    
    std::cout << "\nPosition: r=" << pos[0] << ", θ=" << pos[1] << " (π/4), φ=" << pos[2] << " (π/3)" << std::endl;
    
    // Covariant metric
    auto g = metricSpher.GetCovariantMetric(pos);
    
    std::cout << "\nCovariant metric g_ij:" << std::endl;
    std::cout << "  g_rr = " << g[0][0] << " (expected: 1)" << std::endl;
    std::cout << "  g_θθ = " << g[1][1] << " (expected: r² = 25)" << std::endl;
    std::cout << "  g_φφ = " << g[2][2] << " (expected: r²sin²θ = 25 * 0.5 = 12.5)" << std::endl;
    
    // Contravariant metric
    auto g_inv = metricSpher.GetContravariantMetric(pos);
    
    std::cout << "\nContravariant metric g^ij:" << std::endl;
    std::cout << "  g^rr = " << g_inv[0][0] << " (expected: 1)" << std::endl;
    std::cout << "  g^θθ = " << g_inv[1][1] << " (expected: 1/r² = 0.04)" << std::endl;
    std::cout << "  g^φφ = " << g_inv[2][2] << " (expected: 1/(r²sin²θ) = 0.08)" << std::endl;
    
    // Christoffel symbols
    std::cout << "\nChristoffel symbols (second kind):" << std::endl;
    
    // Γʳθθ = -r
    Real Gamma_r_theta_theta = metricSpher.GetChristoffelSymbolSecondKind(0, 1, 1, pos);
    std::cout << "  Γʳ_θθ = " << Gamma_r_theta_theta << " (expected: -r = -5)" << std::endl;
    
    // Γʳφφ = -r*sin²(θ)
    Real Gamma_r_phi_phi = metricSpher.GetChristoffelSymbolSecondKind(0, 2, 2, pos);
    std::cout << "  Γʳ_φφ = " << Gamma_r_phi_phi << " (expected: -r*sin²θ ≈ -2.5)" << std::endl;
    
    // Γθrθ = 1/r
    Real Gamma_theta_r_theta = metricSpher.GetChristoffelSymbolSecondKind(1, 0, 1, pos);
    std::cout << "  Γθ_rθ = " << Gamma_theta_r_theta << " (expected: 1/r = 0.2)" << std::endl;
    
    // Γθφφ = -sin(θ)*cos(θ)
    Real Gamma_theta_phi_phi = metricSpher.GetChristoffelSymbolSecondKind(1, 2, 2, pos);
    std::cout << "  Γθ_φφ = " << Gamma_theta_phi_phi << " (expected: -sinθ*cosθ ≈ -0.5)" << std::endl;
    
    // Γφrφ = 1/r
    Real Gamma_phi_r_phi = metricSpher.GetChristoffelSymbolSecondKind(2, 0, 2, pos);
    std::cout << "  Γφ_rφ = " << Gamma_phi_r_phi << " (expected: 1/r = 0.2)" << std::endl;
    
    // Γφθφ = cot(θ)
    Real Gamma_phi_theta_phi = metricSpher.GetChristoffelSymbolSecondKind(2, 1, 2, pos);
    std::cout << "  Γφ_θφ = " << Gamma_phi_theta_phi << " (expected: cot(θ) = 1)" << std::endl;
}

void Docs_Demo_MetricTensorCylindrical()
{
    std::cout << "\n=== METRIC TENSOR: CYLINDRICAL ===" << std::endl;
    
    MetricTensorCylindrical metricCyl;
    
    // Position: ρ=3, φ=π/6, z=2
    Vector3Cylindrical pos{3.0, Constants::PI/6, 2.0};
    
    std::cout << "\nPosition: ρ=" << pos[0] << ", φ=" << pos[1] << " (π/6), z=" << pos[2] << std::endl;
    
    // Covariant metric
    auto g = metricCyl.GetCovariantMetric(pos);
    
    std::cout << "\nCovariant metric g_ij:" << std::endl;
    std::cout << "  g_ρρ = " << g[0][0] << " (expected: 1)" << std::endl;
    std::cout << "  g_φφ = " << g[1][1] << " (expected: ρ² = 9)" << std::endl;
    std::cout << "  g_zz = " << g[2][2] << " (expected: 1)" << std::endl;
    
    // Contravariant metric
    auto g_inv = metricCyl.GetContravariantMetric(pos);
    
    std::cout << "\nContravariant metric g^ij:" << std::endl;
    std::cout << "  g^ρρ = " << g_inv[0][0] << " (expected: 1)" << std::endl;
    std::cout << "  g^φφ = " << g_inv[1][1] << " (expected: 1/ρ² ≈ 0.111)" << std::endl;
    std::cout << "  g^zz = " << g_inv[2][2] << " (expected: 1)" << std::endl;
    
    // Arc length along circle at constant ρ and z
    Vector3Cylindrical dr{0.0, 0.1, 0.0};  // dρ=0, dφ=0.1 rad, dz=0
    Real ds_squared = 0.0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            ds_squared += g[i][j] * dr[i] * dr[j];
    Real ds = sqrt(ds_squared);
    
    std::cout << "\nArc length for dφ = 0.1 rad at ρ = 3:" << std::endl;
    std::cout << "  ds = " << ds << " (expected: ρ·dφ = 3 * 0.1 = 0.3)" << std::endl;
}

void Docs_Demo_MetricTensorMinkowski()
{
    std::cout << "\n=== METRIC TENSOR: MINKOWSKI (SPECIAL RELATIVITY) ===" << std::endl;
    
    MetricTensorMinkowski metricMinkowski;
    
    // Event in spacetime: ct=5, x=3, y=0, z=0
    VectorN<Real, 4> event{5.0, 3.0, 0.0, 0.0};
    
    std::cout << "\nSpacetime event: (ct=" << event[0] << ", x=" << event[1] 
              << ", y=" << event[2] << ", z=" << event[3] << ")" << std::endl;
    
    // Get metric (η_μν)
    auto eta = metricMinkowski.GetCovariantMetric(event);
    
    std::cout << "\nMinkowski metric η_μν (signature -,+,+,+):" << std::endl;
    std::cout << "  η_00 = " << eta[0][0] << " (time-time, expected: -1)" << std::endl;
    std::cout << "  η_11 = " << eta[1][1] << " (space-space, expected: +1)" << std::endl;
    std::cout << "  η_22 = " << eta[2][2] << " (space-space, expected: +1)" << std::endl;
    std::cout << "  η_33 = " << eta[3][3] << " (space-space, expected: +1)" << std::endl;
    
    // Compute spacetime interval
    VectorN<Real, 4> displacement{2.0, 1.5, 0.0, 0.0};  // Δct=2, Δx=1.5
    Real interval_squared = 0.0;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            interval_squared += eta[i][j] * displacement[i] * displacement[j];
    
    std::cout << "\nSpacetime interval for Δ(ct, x, y, z) = (2, 1.5, 0, 0):" << std::endl;
    std::cout << "  ds² = " << interval_squared << std::endl;
    std::cout << "  Expected: -(2)² + (1.5)² = -4 + 2.25 = -1.75" << std::endl;
    
    if (interval_squared < 0) {
        std::cout << "  Interpretation: TIMELIKE (causal connection possible)" << std::endl;
        Real proper_time = sqrt(-interval_squared);
        std::cout << "  Proper time: τ = √(-ds²) = " << proper_time << std::endl;
    } else if (interval_squared > 0) {
        std::cout << "  Interpretation: SPACELIKE (causally disconnected)" << std::endl;
        Real proper_distance = sqrt(interval_squared);
        std::cout << "  Proper distance: d = √(ds²) = " << proper_distance << std::endl;
    } else {
        std::cout << "  Interpretation: LIGHTLIKE (photon path)" << std::endl;
    }
    
    // Example 2: Lightlike interval
    std::cout << "\nLightlike example: Δ(ct, x, y, z) = (5, 5, 0, 0) (photon path)" << std::endl;
    VectorN<Real, 4> lightPath{5.0, 5.0, 0.0, 0.0};
    Real lightInterval = 0.0;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            lightInterval += eta[i][j] * lightPath[i] * lightPath[j];
    std::cout << "  ds² = " << lightInterval << " (expected: 0 for light)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////
///                    METRIC FROM COORDINATE TRANSFORMATION
///////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_MetricFromCoordTransf()
{
    std::cout << "\n=== METRIC FROM COORDINATE TRANSFORMATION ===" << std::endl;
    
    // Build metric automatically from Jacobian: g_ij = Σ_k (∂x^k/∂q^i)(∂x^k/∂q^j)
    // Use static transformation instance
    MetricTensorFromCoordTransf<Vector3Spherical, Vector3Cartesian, 3> 
        metricFromTransf(CoordTransfSpherToCart);
    
    // Also get predefined spherical metric for comparison
    MetricTensorSpherical metricPredefined;
    
    Vector3Spherical pos{4.0, Constants::PI/3, Constants::PI/4};  // r=4, θ=60°, φ=45°
    
    std::cout << "\nPosition: r=" << pos[0] << ", θ=" << pos[1] << " (π/3), φ=" << pos[2] << " (π/4)" << std::endl;
    
    // Compare metrics
    auto g_transf = metricFromTransf.GetCovariantMetric(pos);
    auto g_pred = metricPredefined.GetCovariantMetric(pos);
    
    std::cout << "\nMetric from transformation vs Predefined:" << std::endl;
    std::cout << "  g_rr: " << g_transf[0][0] << " vs " << g_pred[0][0] << " (expected: 1)" << std::endl;
    std::cout << "  g_θθ: " << g_transf[1][1] << " vs " << g_pred[1][1] << " (expected: r² = 16)" << std::endl;
    std::cout << "  g_φφ: " << g_transf[2][2] << " vs " << g_pred[2][2] << " (expected: r²sin²θ ≈ 12)" << std::endl;
    
    // Verify they match
    Real diff_rr = std::abs(g_transf[0][0] - g_pred[0][0]);
    Real diff_tt = std::abs(g_transf[1][1] - g_pred[1][1]);
    Real diff_pp = std::abs(g_transf[2][2] - g_pred[2][2]);
    
    std::cout << "\nDifferences: " << diff_rr << ", " << diff_tt << ", " << diff_pp << std::endl;
    std::cout << "Match: " << ((diff_rr < 0.01 && diff_tt < 0.01 && diff_pp < 0.01) ? "YES" : "NO") << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////
///                    COVARIANT DERIVATIVES
///////////////////////////////////////////////////////////////////////////////////////

// Define a velocity field in spherical coordinates
class RotatingFlowSpher : public IVectorFunction<3>
{
public:
    VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
    {
        Real r = pos[0];
        // Solid body rotation: only φ-component
        VectorN<Real, 3> v;
        v[0] = 0.0;              // v^r = 0
        v[1] = 0.0;              // v^θ = 0
        v[2] = 1.0 / r;          // v^φ = 1/r (angular velocity)
        return v;
    }
};

void Docs_Demo_CovariantDerivative()
{
    std::cout << "\n=== COVARIANT DERIVATIVES ===" << std::endl;
    
    MetricTensorSpherical metric;
    RotatingFlowSpher velocity;
    
    Vector3Spherical pos{3.0, Constants::PI/2, 0.0};  // r=3, θ=90°, φ=0
    
    std::cout << "\nVelocity field: v = (0, 0, 1/r) in spherical coordinates" << std::endl;
    std::cout << "Position: r=" << pos[0] << ", θ=π/2, φ=0" << std::endl;
    
    // Compute covariant derivative ∇_j v^i for j=1 (∂/∂θ direction)
    VectorN<Real, 3> covarDeriv_theta = 
        metric.CovariantDerivativeContravar(velocity, 1, pos);
    
    std::cout << "\nCovariant derivative ∇_θ v (derivative in θ direction):" << std::endl;
    std::cout << "  ∇_θ v^r = " << covarDeriv_theta[0] << std::endl;
    std::cout << "  ∇_θ v^θ = " << covarDeriv_theta[1] << std::endl;
    std::cout << "  ∇_θ v^φ = " << covarDeriv_theta[2] << std::endl;
    
    // Covariant derivative in r direction
    VectorN<Real, 3> covarDeriv_r = 
        metric.CovariantDerivativeContravar(velocity, 0, pos);
    
    std::cout << "\nCovariant derivative ∇_r v (derivative in r direction):" << std::endl;
    std::cout << "  ∇_r v^r = " << covarDeriv_r[0] << std::endl;
    std::cout << "  ∇_r v^θ = " << covarDeriv_r[1] << std::endl;
    std::cout << "  ∇_r v^φ = " << covarDeriv_r[2] << std::endl;
    
    // Single component
    Real comp = metric.CovariantDerivativeContravarComp(velocity, 2, 0, pos);
    std::cout << "\nSingle component ∇_r v^φ = " << comp << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////
///                    FIRST KIND CHRISTOFFEL SYMBOLS
///////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_ChristoffelFirstKind()
{
    std::cout << "\n=== CHRISTOFFEL SYMBOLS (FIRST KIND) ===" << std::endl;
    
    MetricTensorSpherical metric;
    
    Vector3Spherical pos{5.0, Constants::PI/4, Constants::PI/3};
    
    std::cout << "\nPosition: r=" << pos[0] << ", θ=π/4, φ=π/3" << std::endl;
    std::cout << "\nFirst kind Christoffel symbols Γ_ijk:" << std::endl;
    
    // Γ_rθθ = g_rr * Γʳ_θθ = 1 * (-r) = -r
    Real Gamma_r_theta_theta_1st = metric.GetChristoffelSymbolFirstKind(0, 1, 1, pos);
    std::cout << "  Γ_rθθ = " << Gamma_r_theta_theta_1st << " (expected: -r = -5)" << std::endl;
    
    // Second kind for comparison
    Real Gamma_r_theta_theta_2nd = metric.GetChristoffelSymbolSecondKind(0, 1, 1, pos);
    std::cout << "  Γʳ_θθ = " << Gamma_r_theta_theta_2nd << " (second kind)" << std::endl;
    
    // Relation: Γ_ijk = g_il * Γˡ_jk
    std::cout << "\nRelation: Γ_rθθ = g_rr * Γʳ_θθ = 1 * " << Gamma_r_theta_theta_2nd 
              << " = " << Gamma_r_theta_theta_2nd << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////
///                    GEOMETRIC APPLICATIONS
///////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_LineElement()
{
    std::cout << "\n=== GEOMETRIC APPLICATION: LINE ELEMENT ===" << std::endl;
    
    MetricTensorSpherical metric;
    
    // Different positions on a sphere of radius r=10
    Real r = 10.0;
    
    // At equator (θ = π/2)
    Vector3Spherical posEquator{r, Constants::PI/2, 0.0};
    auto g_equator = metric.GetCovariantMetric(posEquator);
    
    // At pole (θ = π/8, near pole)
    Vector3Spherical posPole{r, Constants::PI/8, 0.0};
    auto g_pole = metric.GetCovariantMetric(posPole);
    
    // Small displacement in φ (longitude)
    Vector3Spherical dPhi{0.0, 0.0, 0.01};  // dφ = 0.01 rad
    
    // Line element at equator
    Real ds_equator = sqrt(g_equator[2][2]) * dPhi[2];
    
    // Line element near pole
    Real ds_pole = sqrt(g_pole[2][2]) * dPhi[2];
    
    std::cout << "\nFor same dφ = 0.01 rad on sphere r=" << r << ":" << std::endl;
    std::cout << "  At equator (θ=π/2): ds = " << ds_equator 
              << " (r*sin(90°)*dφ = " << r << "*1*0.01)" << std::endl;
    std::cout << "  Near pole (θ=π/8):  ds = " << ds_pole 
              << " (r*sin(22.5°)*dφ = " << r*sin(Constants::PI/8) << "*0.01)" << std::endl;
    std::cout << "\nLatitude lines are shorter near poles!" << std::endl;
}

void Docs_Demo_VolumeElement()
{
    std::cout << "\n=== GEOMETRIC APPLICATION: VOLUME ELEMENT ===" << std::endl;
    
    MetricTensorSpherical metric;
    
    Vector3Spherical pos{5.0, Constants::PI/4, 0.0};
    auto g = metric.GetCovariantMetric(pos);
    
    // Volume element: dV = √det(g) dr dθ dφ
    // For spherical: det(g) = r⁴ sin²(θ), so √det(g) = r² sin(θ)
    
    // Compute determinant
    Real det_g = g[0][0] * g[1][1] * g[2][2];  // Diagonal matrix
    Real sqrt_det_g = sqrt(det_g);
    
    std::cout << "\nAt position r=" << pos[0] << ", θ=" << pos[1] << ":" << std::endl;
    std::cout << "  det(g) = " << det_g << std::endl;
    std::cout << "  √det(g) = " << sqrt_det_g << std::endl;
    std::cout << "  Expected: r²sin(θ) = " << pos[0]*pos[0]*sin(pos[1]) << std::endl;
    std::cout << "\nVolume element: dV = " << sqrt_det_g << " dr dθ dφ" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////
///                           MAIN ENTRY POINT
///////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Metric_Tensor()
{
    std::cout << "***********************************************************" << std::endl;
    std::cout << "*****              METRIC TENSOR FRAMEWORK            *****" << std::endl;
    std::cout << "***********************************************************" << std::endl;

    // Predefined metric tensors
    Docs_Demo_MetricTensorCartesian();
    Docs_Demo_MetricTensorSpherical();
    Docs_Demo_MetricTensorCylindrical();
    Docs_Demo_MetricTensorMinkowski();
    
    // Building metric from transformation
    Docs_Demo_MetricFromCoordTransf();
    
    // Christoffel symbols
    Docs_Demo_ChristoffelFirstKind();
    
    // Covariant derivatives
    Docs_Demo_CovariantDerivative();
    
    // Geometric applications
    std::cout << "\n***********************************************************" << std::endl;
    std::cout << "*****          GEOMETRIC APPLICATIONS                 *****" << std::endl;
    std::cout << "***********************************************************" << std::endl;
    
    Docs_Demo_LineElement();
    Docs_Demo_VolumeElement();
    
    std::cout << "\n***********************************************************" << std::endl;
    std::cout << "*****        METRIC TENSOR DEMO COMPLETE              *****" << std::endl;
    std::cout << "***********************************************************" << std::endl;
}
