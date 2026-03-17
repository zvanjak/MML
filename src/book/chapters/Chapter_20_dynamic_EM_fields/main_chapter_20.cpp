#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Vector/VectorN.h"
#include "mml/base/Vector/VectorTypes.h"
#include "mml/base/Tensor.h"

#include "mml/core/Derivation.h"
#include "mml/algorithms/FieldLineTracer.h"
#include "mml/tools/Serializer.h"

#include "mpl/Electromagnetism/EMTensor.h"
#include "mpl/Electromagnetism/LienardWiechertPotential.h"

#include "mpl/SpecialRelativity/LorentzTransformation.h"
#endif

#include <iomanip>
#include <cmath>
#include <filesystem>

using namespace MML;
using namespace MPL;

//////////////////////////////////////////////////////////////////////////////
/// SCENARIO 1: Stationary Charge - Verify EM Tensor Transformation
/// 
/// A charge at rest in frame S has only an electric field (Coulomb's law).
/// An observer in frame S' moving at velocity v sees:
///   - A compressed ("pancaked") electric field
///   - A magnetic field (moving charge = current!)
///
/// We verify that the EM tensor transforms correctly by comparing:
///   1) Tensor transformation of F^μν
///   2) Analytic field transformation formulas E'_⊥ = γE_⊥, B' = v × E
///
/// Note: The Heaviside-Feynman formula gives the field pattern for a charge
/// that IS moving. That's verified separately in Demo1c.
//////////////////////////////////////////////////////////////////////////////

void Demo1_StationaryCharge_EMTensorVerification()
{
	std::cout << "\n╔═══════════════════════════════════════════════════════════════════╗\n";
	std::cout << "║   SCENARIO 1: EM Tensor Transformation Verification               ║\n";
	std::cout << "║   Tensor F^μν transforms as F'=ΛFΛᵀ                                ║\n";
	std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
	std::cout << "║   Compare tensor transformation with analytic E,B formulas        ║\n";
	std::cout << "║   E' = γE_⊥, B'_y = -γβE_z, B'_z = +γβE_y (for B=0)              ║\n";
	std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n\n";

	Real beta = 0.6;  // v/c
	Real gamma = 1.0 / std::sqrt(1.0 - beta * beta);
	Real charge = 1.0;  // Unit charge (q/4πε₀ = 1)

	std::cout << "Boost velocity: β = " << beta << ", γ = " << std::fixed << std::setprecision(6) << gamma << "\n\n";

	// Test points (field position)
	std::vector<Vector3Cartesian> testPoints = {
		Vector3Cartesian(1.0, 0.0, 0.0),   // Along boost direction - no perpendicular component
		Vector3Cartesian(0.0, 1.0, 0.0),   // Perpendicular to boost
		Vector3Cartesian(0.0, 0.0, 1.0),   // Perpendicular (z)
		Vector3Cartesian(1.0, 1.0, 0.0),   // Mixed
		Vector3Cartesian(3.0, 4.0, 0.0),   // General point
	};

	std::cout << std::setw(20) << "Field Point" 
	          << " │ " << std::setw(10) << "Method"
	          << " │ " << std::setw(30) << "E' field"
	          << " │ " << std::setw(30) << "B' field"
	          << " │ Match?\n";
	std::cout << std::string(110, '-') << "\n";

	bool allPassed = true;
	CoordTransfLorentzXAxis lorentzTransf(beta);  // Boost by +β along x-axis

	for (const auto& point : testPoints)
	{
		// E field in rest frame (Coulomb's law)
		Real r = point.NormL2();
		Vector3Cartesian E_rest = point * (charge / (r * r * r));
		Vector3Cartesian B_rest(0.0, 0.0, 0.0);

		// ══════════════════════════════════════════════════════════════════
		// METHOD 1: Full tensor transformation
		// ══════════════════════════════════════════════════════════════════
		Tensor2<4> F_rest = GetEMTensorContravariant(E_rest, B_rest);
		Vector4Minkowski point4({0.0, point.X(), point.Y(), point.Z()});
		Tensor2<4> F_boosted = lorentzTransf.transfTensor2(F_rest, point4);
		
		Vector3Cartesian E_tensor, B_tensor;
		GetEandBFromEMTensorContravariant(F_boosted, E_tensor, B_tensor);

		// ══════════════════════════════════════════════════════════════════
		// METHOD 2: Analytic field transformation formulas
		// For boost along +x with velocity β (matching CoordTransfLorentzXAxis):
		// E'_x = E_x, E'_y = γE_y, E'_z = γE_z (when B=0)
		// B'_x = 0, B'_y = -γβE_z, B'_z = +γβE_y (when B=0)
		// Note: Sign convention matches our Lorentz transform implementation
		// ══════════════════════════════════════════════════════════════════
		Vector3Cartesian E_analytic(
			E_rest.X(),
			gamma * E_rest.Y(),
			gamma * E_rest.Z()
		);
		Vector3Cartesian B_analytic(
			0.0,
			-gamma * beta * E_rest.Z(),
			+gamma * beta * E_rest.Y()
		);

		// Compare
		Real E_diff = (E_tensor - E_analytic).NormL2();
		Real B_diff = (B_tensor - B_analytic).NormL2();
		bool match = (E_diff < 1e-10 && B_diff < 1e-10);
		if (!match) allPassed = false;

		// Print
		std::ostringstream pointStr;
		pointStr << "(" << point.X() << ", " << point.Y() << ", " << point.Z() << ")";
		
		std::cout << std::setw(20) << pointStr.str()
		          << " │ " << std::setw(10) << "Tensor"
		          << " │ " << std::fixed << std::setprecision(6)
		          << "(" << std::setw(8) << E_tensor.X() << ", " << std::setw(8) << E_tensor.Y() << ", " << std::setw(8) << E_tensor.Z() << ")"
		          << " │ (" << std::setw(8) << B_tensor.X() << ", " << std::setw(8) << B_tensor.Y() << ", " << std::setw(8) << B_tensor.Z() << ")"
		          << " │\n";
		std::cout << std::setw(20) << ""
		          << " │ " << std::setw(10) << "Analytic"
		          << " │ (" << std::setw(8) << E_analytic.X() << ", " << std::setw(8) << E_analytic.Y() << ", " << std::setw(8) << E_analytic.Z() << ")"
		          << " │ (" << std::setw(8) << B_analytic.X() << ", " << std::setw(8) << B_analytic.Y() << ", " << std::setw(8) << B_analytic.Z() << ")"
		          << " │ " << (match ? "✓ PASS" : "✗ FAIL") << "\n";
		std::cout << std::string(110, '-') << "\n";
	}

	std::cout << "\n╔═══════════════════════════════════════════════════════════════════╗\n";
	if (allPassed)
		std::cout << "║   ✓ ALL TESTS PASSED - EM Tensor is a true Lorentz tensor!       ║\n";
	else
		std::cout << "║   ✗ SOME TESTS FAILED - Check tensor transformation!             ║\n";
	std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n";
}

//////////////////////////////////////////////////////////////////////////////
/// SCENARIO 2: Infinite Current-Carrying Wire
/// 
/// This is the classic demonstration that MAGNETISM IS A RELATIVISTIC EFFECT!
///
/// In the wire rest frame S:
///   - Positive ions at rest, negative electrons moving at drift velocity v_d
///   - Net charge density = 0 (wire is electrically neutral)
///   - Current I creates magnetic field B = μ₀I/(2πr) circling the wire
///   - A test charge at rest feels no force (no E field)
///
/// In frame S' moving with the electrons (at velocity v_d):
///   - Electrons appear at rest, positive ions moving at -v_d  
///   - Due to LENGTH CONTRACTION:
///     * Positive ions are closer together → positive charge density increases
///     * Electrons are spread apart → negative charge density decreases
///   - Net POSITIVE charge appears! → Electric field appears!
///   - The magnetic force in S becomes an electric force in S'
///
/// The EM tensor transformation shows this automatically!
//////////////////////////////////////////////////////////////////////////////

void Demo2_InfiniteWire_MagnetismIsRelativistic()
{
	std::cout << "\n╔═══════════════════════════════════════════════════════════════════╗\n";
	std::cout << "║   SCENARIO 2: Infinite Current-Carrying Wire                      ║\n";
	std::cout << "║   Why Magnetism is a Relativistic Effect                          ║\n";
	std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
	std::cout << "║   Wire frame S:  Neutral wire with current I along z-axis         ║\n";
	std::cout << "║                  B field circles wire, no E field                 ║\n";
	std::cout << "║   Moving frame S': Observer moves with electrons at v_drift       ║\n";
	std::cout << "║                    E field appears due to length contraction!     ║\n";
	std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n\n";

	// Wire along z-axis, current flows in +z direction
	// Electron drift velocity (typically ~mm/s, but we use relativistic for dramatic effect)
	Real beta_drift = 0.6;  // Drift velocity as fraction of c (exaggerated for demo)
	Real gamma = 1.0 / std::sqrt(1.0 - beta_drift * beta_drift);
	
	// Current creates B field: B = μ₀I/(2πr) in φ direction
	// We use normalized units: B₀ = μ₀I/(2π) = 1
	Real B0 = 1.0;  // Magnetic field strength parameter
	
	std::cout << "Electron drift velocity: β = " << beta_drift << " (exaggerated for demo)\n";
	std::cout << "Lorentz factor:          γ = " << std::fixed << std::setprecision(6) << gamma << "\n\n";

	// Test points at different distances from wire (wire along z-axis)
	std::vector<Vector3Cartesian> testPoints = {
		Vector3Cartesian(1.0, 0.0, 0.0),   // On +x axis, distance r=1
		Vector3Cartesian(0.0, 1.0, 0.0),   // On +y axis, distance r=1
		Vector3Cartesian(2.0, 0.0, 0.0),   // Farther out
		Vector3Cartesian(1.0, 1.0, 0.0),   // At 45°, distance r=√2
	};

	std::cout << "═══════════════════════════════════════════════════════════════════\n";
	std::cout << "Wire Rest Frame (S): B field only, no E field\n";
	std::cout << "═══════════════════════════════════════════════════════════════════\n\n";

	std::cout << std::setw(20) << "Point" 
	          << " │ " << std::setw(30) << "E (should be 0)"
	          << " │ " << std::setw(30) << "B (circles wire)\n";
	std::cout << std::string(90, '-') << "\n";

	bool allPassed = true;

	for (const auto& point : testPoints)
	{
		// Distance from z-axis
		Real r = std::sqrt(point.X() * point.X() + point.Y() * point.Y());
		
		// In wire rest frame: E = 0, B = B₀/r in φ̂ direction
		// φ̂ = (-y/r, x/r, 0) for current in +z direction (right-hand rule)
		Vector3Cartesian E_wire(0.0, 0.0, 0.0);
		Vector3Cartesian B_wire(
			-point.Y() / (r * r) * B0,  // Bx = -y/r² × B₀
			 point.X() / (r * r) * B0,  // By = +x/r² × B₀  
			 0.0                         // Bz = 0
		);

		std::ostringstream pointStr, EStr, BStr;
		pointStr << "(" << point.X() << ", " << point.Y() << ", " << point.Z() << ")";
		EStr << std::fixed << std::setprecision(6) 
		     << "(" << E_wire.X() << ", " << E_wire.Y() << ", " << E_wire.Z() << ")";
		BStr << std::fixed << std::setprecision(6)
		     << "(" << B_wire.X() << ", " << B_wire.Y() << ", " << B_wire.Z() << ")";

		std::cout << std::setw(20) << pointStr.str()
		          << " │ " << std::setw(30) << EStr.str()
		          << " │ " << std::setw(30) << BStr.str() << "\n";

		// ══════════════════════════════════════════════════════════════════
		// Transform to electron rest frame (moving at +β in z direction)
		// ══════════════════════════════════════════════════════════════════
		
		// We need a Lorentz boost along z-axis
		// For boost along z: E'_x = γ(E_x + βB_y), E'_y = γ(E_y - βB_x)
		//                    B'_x = γ(B_x - βE_y), B'_y = γ(B_y + βE_x)
		// But let's use tensor transformation for generality
	}

	std::cout << "\n═══════════════════════════════════════════════════════════════════\n";
	std::cout << "Electron Rest Frame (S'): E field appears!\n";
	std::cout << "Observer moves with electrons at β = " << beta_drift << " along z-axis\n";
	std::cout << "═══════════════════════════════════════════════════════════════════\n\n";

	std::cout << std::setw(20) << "Point" 
	          << " │ " << std::setw(30) << "E' (radial - appears!)"
	          << " │ " << std::setw(30) << "B' (still circles)\n";
	std::cout << std::string(90, '-') << "\n";

	for (const auto& point : testPoints)
	{
		Real r = std::sqrt(point.X() * point.X() + point.Y() * point.Y());
		
		// Wire frame fields
		Vector3Cartesian E_wire(0.0, 0.0, 0.0);
		Vector3Cartesian B_wire(
			-point.Y() / (r * r) * B0,
			 point.X() / (r * r) * B0,
			 0.0
		);

		// Form EM tensor
		Tensor2<4> F_wire = GetEMTensorContravariant(E_wire, B_wire);

		// Lorentz boost along z-axis
		// For boost along z with velocity β:
		// Λ = | γ      0  0  -γβ |
		//     | 0      1  0   0  |
		//     | 0      0  1   0  |
		//     |-γβ     0  0   γ  |
		//
		// We need to use CoordTransfLorentzGeneral for z-axis boost
		Vec3Cart zDir(0.0, 0.0, 1.0);
		CoordTransfLorentzGeneral lorentzZ(beta_drift, zDir);
		
		Vector4Minkowski point4({0.0, point.X(), point.Y(), point.Z()});
		Tensor2<4> F_electron = lorentzZ.transfTensor2(F_wire, point4);

		Vector3Cartesian E_electron, B_electron;
		GetEandBFromEMTensorContravariant(F_electron, E_electron, B_electron);

		// ══════════════════════════════════════════════════════════════════
		// Verify with analytic formulas for z-boost
		// E'_x = γ(E_x + βB_y), E'_y = γ(E_y - βB_x), E'_z = E_z
		// B'_x = γ(B_x - βE_y), B'_y = γ(B_y + βE_x), B'_z = B_z
		// ══════════════════════════════════════════════════════════════════
		Vector3Cartesian E_analytic(
			gamma * (E_wire.X() + beta_drift * B_wire.Y()),
			gamma * (E_wire.Y() - beta_drift * B_wire.X()),
			E_wire.Z()
		);
		Vector3Cartesian B_analytic(
			gamma * (B_wire.X() - beta_drift * E_wire.Y()),
			gamma * (B_wire.Y() + beta_drift * E_wire.X()),
			B_wire.Z()
		);

		// Check match
		Real E_diff = (E_electron - E_analytic).NormL2();
		Real B_diff = (B_electron - B_analytic).NormL2();
		bool match = (E_diff < 1e-9 && B_diff < 1e-9);
		if (!match) allPassed = false;

		std::ostringstream pointStr, EStr, BStr;
		pointStr << "(" << point.X() << ", " << point.Y() << ", " << point.Z() << ")";
		EStr << std::fixed << std::setprecision(6)
		     << "(" << E_electron.X() << ", " << E_electron.Y() << ", " << E_electron.Z() << ")";
		BStr << std::fixed << std::setprecision(6)
		     << "(" << B_electron.X() << ", " << B_electron.Y() << ", " << B_electron.Z() << ")";

		std::cout << std::setw(20) << pointStr.str()
		          << " │ " << std::setw(30) << EStr.str()
		          << " │ " << std::setw(30) << BStr.str()
		          << " │ " << (match ? "✓" : "✗") << "\n";
	}

	std::cout << "\n╔═══════════════════════════════════════════════════════════════════╗\n";
	std::cout << "║   Physical Interpretation:                                        ║\n";
	std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
	std::cout << "║   In wire frame: Neutral wire, B field from current               ║\n";
	std::cout << "║   In electron frame: Positive charge density appears!             ║\n";
	std::cout << "║   → Length contraction makes positive ions denser                 ║\n";
	std::cout << "║   → The magnetic force IS the electric force seen differently!    ║\n";
	std::cout << "╠═══════════════════════════════════════════════════════════════════╣\n";
	if (allPassed)
		std::cout << "║   ✓ Tensor transformation correctly predicts E field appearance! ║\n";
	else
		std::cout << "║   ✗ Check transformation - some tests failed                     ║\n";
	std::cout << "╚═══════════════════════════════════════════════════════════════════╝\n";
}

//////////////////////////////////////////////////////////////////////////////
/// SCENARIO 3: Dual Validation - Moving Charge Field Formulas
///
/// Compare two independent methods for computing the field of a moving charge:
///
/// Method 1: HEAVISIDE-FEYNMAN FORMULA  
///   - Direct formula for field of uniformly moving charge
///   - E = q(1-β²)/[r²(1-β²sin²θ)^(3/2)] · r̂
///   - B = v × E
///
/// Method 2: LIÉNARD-WIECHERT FIELDS
///   - General retarded potential theory
///   - For uniform velocity, simplifies to the same result as H-F
///
/// Both methods must give identical results for uniform velocity!
//////////////////////////////////////////////////////////////////////////////
 
void Demo3_TripleValidation_MovingCharge()
{
	std::cout << "\n╔══════════════════════════════════════════════════════════════════════════════╗\n";
	std::cout << "║   SCENARIO 3: Moving Charge Field Validation                                 ║\n";
	std::cout << "║   Heaviside-Feynman vs Liénard-Wiechert                                      ║\n";
	std::cout << "╠══════════════════════════════════════════════════════════════════════════════╣\n";
	std::cout << "║   Method 1: Heaviside-Feynman      E ∝ (1-β²)/(1-β²sin²θ)^(3/2)              ║\n";
	std::cout << "║   Method 2: Liénard-Wiechert       Same result for uniform velocity          ║\n";
	std::cout << "╚══════════════════════════════════════════════════════════════════════════════╝\n";

	// Charge at origin, velocity along x-axis
	Real charge = 1.0;
	Real beta = 0.8;  // v = 0.8c - highly relativistic
	Real gamma = 1.0 / std::sqrt(1.0 - beta * beta);
	
	std::cout << "\nCharge at origin, moving at β = " << beta << " along x-axis\n";
	std::cout << "Lorentz factor: γ = " << std::fixed << std::setprecision(6) << gamma << "\n\n";

	// Test points at various angles around the charge
	std::vector<Vector3Cartesian> testPoints = {
		{1.0, 0.0, 0.0},   // Along velocity (θ = 0)
		{-1.0, 0.0, 0.0},  // Opposite velocity (θ = π)
		{0.0, 1.0, 0.0},   // Perpendicular (θ = π/2)
		{0.0, 0.0, 1.0},   // Perpendicular in z
		{1.0, 1.0, 0.0},   // 45 degrees
		{2.0, 1.0, 0.0},   // Different angle
		{0.5, 0.5, 0.5},   // 3D diagonal
	};

	// ═══════════════════════════════════════════════════════════════════════
	// METHOD 1: Tensor Transformation
	// We have a stationary charge in frame S'. The observer (us) moves at -β.
	// Equivalently: charge moves at +β in our frame S.
	// We boost from charge rest frame to lab frame.
	// ═══════════════════════════════════════════════════════════════════════

	std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
	std::cout << " Point         │ Method      │ E field                    │ B field                    \n";
	std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";

	bool allPassed = true;
	Real tolerance = 1e-10;

	for (const auto& point : testPoints) {
		// ─────────────────────────────────────────────────────────────────────
		// METHOD 1: Heaviside-Feynman Formula
		// ─────────────────────────────────────────────────────────────────────
		// Direct formula for field of uniformly moving charge:
		// E = q(1-β²) / [r²(1-β²sin²θ)^(3/2)] · r̂
		// B = v × E
		Vector3Cartesian E_HF, B_HF;
		Vector3Cartesian chargePos(0, 0, 0);
		Vector3Cartesian velocity(beta, 0, 0);
		GetFieldsOfMovingCharge(chargePos, point, velocity, charge, E_HF, B_HF);

		// ─────────────────────────────────────────────────────────────────────
		// METHOD 2: Liénard-Wiechert (for uniform velocity)
		// ─────────────────────────────────────────────────────────────────────
		// For uniform velocity, L-W reduces to the same formula as H-F
		Vector3Cartesian E_LW, B_LW;
		GetLienardWiechertFieldsUniformVelocity(chargePos, point, velocity, charge, E_LW, B_LW);

		// ─────────────────────────────────────────────────────────────────────
		// Compare: H-F vs L-W (they must be identical for uniform velocity)
		// ─────────────────────────────────────────────────────────────────────
		Real diff_HFLW = (E_HF - E_LW).NormL2() + (B_HF - B_LW).NormL2();

		bool methodsAgree = (diff_HFLW < tolerance);
		if (!methodsAgree) allPassed = false;

		// Format output
		auto formatVec = [](const Vector3Cartesian& v) {
			std::ostringstream ss;
			ss << std::fixed << std::setprecision(4) 
			   << "(" << std::setw(8) << v.X() << ", " << std::setw(8) << v.Y() << ", " << std::setw(8) << v.Z() << ")";
			return ss.str();
		};

		auto formatPoint = [](const Vector3Cartesian& p) {
			std::ostringstream ss;
			ss << "(" << p.X() << "," << p.Y() << "," << p.Z() << ")";
			return ss.str();
		};

		// Show only H-F and L-W (tensor shown for reference but not validated)
		std::cout << std::setw(14) << formatPoint(point) << " │ Heaviside   │ " 
		          << std::setw(27) << formatVec(E_HF) << " │ " << std::setw(27) << formatVec(B_HF) << "\n";
		std::cout << "               │ Liénard-W   │ " 
		          << std::setw(27) << formatVec(E_LW) << " │ " << std::setw(27) << formatVec(B_LW) 
		          << " │ " << (methodsAgree ? "✓ AGREE" : "✗ DIFFER") << "\n";
		std::cout << "───────────────┼─────────────┼────────────────────────────┼────────────────────────────\n";
	}

	std::cout << "\n╔══════════════════════════════════════════════════════════════════════════════╗\n";
	if (allPassed) {
		std::cout << "║   ✓ HEAVISIDE-FEYNMAN AND LIÉNARD-WIECHERT AGREE PERFECTLY!                 ║\n";
		std::cout << "╠══════════════════════════════════════════════════════════════════════════════╣\n";
		std::cout << "║   For uniformly moving charges, both formulas give identical results.       ║\n";
		std::cout << "║   This confirms our moving charge field implementations are correct.        ║\n";
	} else {
		std::cout << "║   ✗ METHODS DISAGREE - Check Implementation!                                 ║\n";
	}
	std::cout << "╚══════════════════════════════════════════════════════════════════════════════╝\n";
}

//////////////////////////////////////////////////////////////////////////////
/// SCENARIO 4: Field Line Visualization - The Pancake Effect!
///
/// As a charge moves faster, its electric field gets "squished" perpendicular
/// to the direction of motion. At v → c, the field becomes a thin disk!
///
/// Physics:
/// - Field along velocity (θ=0°): E ∝ (1-β²) → diminishes as β→1
/// - Field perpendicular (θ=90°): E ∝ γ → INCREASES as β→1
/// - Ratio: E_perp / E_parallel = γ³
///
/// At v = 0.99c, γ ≈ 7.1, so γ³ ≈ 358 - field is 358× stronger perpendicular!
//////////////////////////////////////////////////////////////////////////////

void Demo4_FieldVisualization_PancakeEffect()
{
	std::cout << "\n╔══════════════════════════════════════════════════════════════════════════════╗\n";
	std::cout << "║   SCENARIO 4: The Pancake Effect - Field Compression Visualization          ║\n";
	std::cout << "╠══════════════════════════════════════════════════════════════════════════════╣\n";
	std::cout << "║   As velocity → c, the E field gets compressed into a thin disk!            ║\n";
	std::cout << "║   Charge moves along x-axis, field shown in x-y plane                       ║\n";
	std::cout << "╚══════════════════════════════════════════════════════════════════════════════╝\n\n";

	Real charge = 1.0;
	
	// Different velocities to show compression
	std::vector<Real> betas = {0.0, 0.5, 0.9, 0.99};
	
	// Show field magnitude at key angles
	std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
	std::cout << " E-field magnitude at distance r=1 from charge (normalized to static case)\n";
	std::cout << "═══════════════════════════════════════════════════════════════════════════════\n\n";
	
	std::cout << "      β (v/c) │      γ │ E(θ=0°) │ E(θ=90°) │  Ratio  │ Field Shape\n";
	std::cout << "──────────────┼────────┼─────────┼──────────┼─────────┼────────────────────────\n";
	
	for (Real beta : betas) {
		Real gamma = (beta < 1e-10) ? 1.0 : 1.0 / std::sqrt(1.0 - beta * beta);
		
		// Field along velocity (θ = 0, point on x-axis)
		Vector3Cartesian E_parallel, B_parallel;
		GetFieldsOfMovingCharge({0,0,0}, {1,0,0}, {beta,0,0}, charge, E_parallel, B_parallel);
		Real E_par_mag = E_parallel.NormL2();
		
		// Field perpendicular (θ = 90°, point on y-axis)
		Vector3Cartesian E_perp, B_perp;
		GetFieldsOfMovingCharge({0,0,0}, {0,1,0}, {beta,0,0}, charge, E_perp, B_perp);
		Real E_perp_mag = E_perp.NormL2();
		
		// Ratio shows how "pancaked" the field is
		Real ratio = (E_par_mag > 1e-15) ? E_perp_mag / E_par_mag : 0;
		
		// Visual indicator of pancaking
		std::string shape;
		if (beta < 0.01) shape = "○ Spherical";
		else if (beta < 0.6) shape = "◯ Slightly squished";
		else if (beta < 0.95) shape = "⬭ Elliptical";
		else shape = "━ PANCAKE! 🥞";
		
		std::cout << std::fixed << std::setprecision(2)
		          << std::setw(13) << beta << " │"
		          << std::setw(7) << gamma << " │"
		          << std::setw(8) << std::setprecision(4) << E_par_mag << " │"
		          << std::setw(9) << std::setprecision(4) << E_perp_mag << " │"
		          << std::setw(8) << std::setprecision(1) << ratio << " │ "
		          << shape << "\n";
	}
	
	std::cout << "\n";
	
	// Now show ASCII art visualization of the field pattern!
	std::cout << "═══════════════════════════════════════════════════════════════════════════════\n";
	std::cout << " Field pattern visualization (charge at center, velocity → right)\n";
	std::cout << " Dots show relative field strength at each angle\n";
	std::cout << "═══════════════════════════════════════════════════════════════════════════════\n\n";
	
	for (Real beta : betas) {
		Real gamma = (beta < 1e-10) ? 1.0 : 1.0 / std::sqrt(1.0 - beta * beta);
		
		std::cout << "β = " << std::fixed << std::setprecision(2) << beta 
		          << " (γ = " << std::setprecision(2) << gamma << ")";
		if (beta > 0.98) std::cout << " 🥞 PANCAKE!";
		std::cout << "\n";
		
		// Create ASCII visualization of field pattern
		// Grid is 31x15, charge at center (15, 7)
		const int width = 41;
		const int height = 21;
		const int cx = width / 2;
		const int cy = height / 2;
		
		char grid[height][width + 1];
		
		// Initialize grid with spaces
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				grid[y][x] = ' ';
			}
			grid[y][width] = '\0';
		}
		
		// Draw border
		for (int x = 0; x < width; x++) {
			grid[0][x] = '-';
			grid[height-1][x] = '-';
		}
		for (int y = 0; y < height; y++) {
			grid[y][0] = '|';
			grid[y][width-1] = '|';
		}
		grid[0][0] = '+'; grid[0][width-1] = '+';
		grid[height-1][0] = '+'; grid[height-1][width-1] = '+';
		
		// Mark center (charge position)
		grid[cy][cx] = 'Q';
		
		// Draw velocity arrow
		grid[cy][cx+1] = '-';
		grid[cy][cx+2] = '-';
		grid[cy][cx+3] = '>';
		
		// Sample field at various angles and radii
		// Normalize to max field strength for this velocity
		Real maxE = 0;
		std::vector<std::pair<int, int>> points;
		std::vector<Real> fieldStrengths;
		
		for (int angle = 0; angle < 360; angle += 10) {
			Real theta = angle * 3.14159265 / 180.0;
			Real r = 8.0;  // Radius in grid units
			
			// Position in grid coordinates
			int gx = cx + static_cast<int>(r * std::cos(theta));
			int gy = cy - static_cast<int>(r * 0.5 * std::sin(theta));  // 0.5 for aspect ratio
			
			if (gx < 1 || gx >= width-1 || gy < 1 || gy >= height-1) continue;
			
			// Calculate field at this point
			Real fx = std::cos(theta);
			Real fy = std::sin(theta);
			Vector3Cartesian E, B;
			GetFieldsOfMovingCharge({0,0,0}, {fx, fy, 0}, {beta,0,0}, charge, E, B);
			Real Emag = E.NormL2();
			
			points.push_back({gx, gy});
			fieldStrengths.push_back(Emag);
			if (Emag > maxE) maxE = Emag;
		}
		
		// Place markers based on relative field strength
		for (size_t i = 0; i < points.size(); i++) {
			Real relStrength = fieldStrengths[i] / maxE;
			char marker;
			if (relStrength > 0.9) marker = '#';
			else if (relStrength > 0.7) marker = '*';
			else if (relStrength > 0.5) marker = 'o';
			else if (relStrength > 0.3) marker = '.';
			else marker = ' ';
			
			grid[points[i].second][points[i].first] = marker;
		}
		
		// Print the grid
		for (int y = 0; y < height; y++) {
			std::cout << "    " << grid[y] << "\n";
		}
		std::cout << "\n";
	}
	
	std::cout << "╔══════════════════════════════════════════════════════════════════════════════╗\n";
	std::cout << "║   Key: # = strongest, * = strong, o = medium, . = weak                      ║\n";
	std::cout << "║   Q = charge position, --> = velocity direction                             ║\n";
	std::cout << "╠══════════════════════════════════════════════════════════════════════════════╣\n";
	std::cout << "║   Notice how the field concentrates perpendicular to velocity!              ║\n";
	std::cout << "║   At β = 0.99, almost ALL the field is in a thin disk - the PANCAKE! 🥞    ║\n";
	std::cout << "╚══════════════════════════════════════════════════════════════════════════════╝\n";
}

//////////////////////////////////////////////////////////////////////////////
/// SCENARIO 5: Field Line Tracing - Export for Proper Visualization!
///
/// Instead of ASCII art, we trace the actual field lines using the MML
/// FieldLineTracer and export them to .mml files for visualization.
///
/// This demonstrates the FieldLineTracer algorithm with a 2D projection
/// of the electric field from a moving charge (the pancake effect).
//////////////////////////////////////////////////////////////////////////////

// Wrapper class to make the E-field callable by FieldLineTracer
class MovingChargeEField2D : public IVectorFunction<2>
{
private:
	Real m_beta;    // v/c
	Real m_charge;
	
public:
	MovingChargeEField2D(Real beta, Real charge = 1.0) 
		: m_beta(beta), m_charge(charge) {}
	
	VectorN<Real, 2> operator()(const VectorN<Real, 2>& pos) const override
	{
		// Convert 2D position to 3D (in x-y plane)
		Vector3Cartesian point3D(pos[0], pos[1], 0.0);
		Vector3Cartesian chargePos(0.0, 0.0, 0.0);
		Vector3Cartesian velocity(m_beta, 0.0, 0.0);
		
		Vector3Cartesian E, B;
		GetFieldsOfMovingCharge(chargePos, point3D, velocity, m_charge, E, B);
		
		// Return 2D projection (x,y components)
		return VectorN<Real, 2>({E.X(), E.Y()});
	}
};

void Demo5_FieldLineTracing_Export()
{
	std::cout << "\n╔══════════════════════════════════════════════════════════════════════════════╗\n";
	std::cout << "║   SCENARIO 5: Field Line Tracing with FieldLineTracer                       ║\n";
	std::cout << "╠══════════════════════════════════════════════════════════════════════════════╣\n";
	std::cout << "║   Tracing E-field lines and exporting to .mml files for visualization       ║\n";
	std::cout << "║   Watch the pancake form as β increases!                                    ║\n";
	std::cout << "╚══════════════════════════════════════════════════════════════════════════════╝\n\n";

	// Create output directory if it doesn't exist
	std::string outputDir = "field_lines_output";
	std::filesystem::create_directories(outputDir);

	std::vector<Real> betas = {0.0, 0.5, 0.9, 0.99};
	Real charge = 1.0;
	
	// FieldLineTracer configuration
	FieldLineTracer<2>::Config config;
	config.stepSize = 0.05;
	config.maxLength = 5.0;       // Stop after this arc length
	config.maxPoints = 500;
	config.minFieldMagnitude = 1e-6;
	config.traceForward = true;   // Away from positive charge
	config.traceBackward = false;
	
	// Bounding box for tracing (avoid going to infinity)
	BoundingBox<2> bounds;
	bounds.minCorner = VectorN<Real, 2>({-5.0, -5.0});
	bounds.maxCorner = VectorN<Real, 2>({5.0, 5.0});
	
	// Seed points: ring around the charge
	int numLines = 24;
	Real seedRadius = 0.1;  // Small radius near charge
	std::vector<VectorN<Real, 2>> seedPoints;
	
	for (int i = 0; i < numLines; ++i) {
		Real angle = 2.0 * 3.14159265358979 * i / numLines;
		seedPoints.push_back(VectorN<Real, 2>({
			seedRadius * std::cos(angle),
			seedRadius * std::sin(angle)
		}));
	}
	
	std::cout << "Tracing " << numLines << " field lines for each velocity...\n\n";
	std::cout << "      β      │  γ factor  │ Lines Traced │ Output File\n";
	std::cout << "─────────────┼────────────┼──────────────┼───────────────────────────────\n";
	
	for (Real beta : betas) {
		Real gamma = (beta < 1e-10) ? 1.0 : 1.0 / std::sqrt(1.0 - beta * beta);
		
		// Create field function for this velocity
		MovingChargeEField2D field(beta, charge);
		
		// Create tracer with config
		FieldLineTracer<2> tracer(config);
		
		// Trace from seed points
		std::vector<FieldLine<2>> fieldLines = tracer.TraceFromSeedPoints(field, seedPoints, bounds);
		
		// Generate filename
		std::string filename = outputDir + "/Pancake_beta_" + 
			std::to_string(static_cast<int>(beta * 100)) + ".mml";
		
		// Save to file (order: lines, title, fileName)
		SerializeResult result = Serializer::SaveFieldLines2D(fieldLines, 
			"Electric field lines of moving charge, beta = " + std::to_string(beta),
			filename);
		
		std::cout << std::fixed << std::setprecision(2)
		          << std::setw(12) << beta << " │"
		          << std::setw(11) << std::setprecision(4) << gamma << " │"
		          << std::setw(13) << fieldLines.size() << " │ "
		          << (result.success ? filename : "FAILED") << "\n";
	}
	
	std::cout << "\n╔══════════════════════════════════════════════════════════════════════════════╗\n";
	std::cout << "║   Field lines exported! Open .mml files with MML visualizer to see:         ║\n";
	std::cout << "║   • β=0.00: Perfect radial pattern (Coulomb field)                          ║\n";
	std::cout << "║   • β=0.50: Slight compression perpendicular to velocity                    ║\n";
	std::cout << "║   • β=0.90: Significant pancaking visible!                                  ║\n";
	std::cout << "║   • β=0.99: Extreme pancake - field lines nearly in a disk! 🥞              ║\n";
	std::cout << "╚══════════════════════════════════════════════════════════════════════════════╝\n";
}

void Chapter20_EM_field_investigations()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****        CHAPTER 20 - Dynamic EM Fields & Tensor Verification    ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Demo1_StationaryCharge_EMTensorVerification();
	Demo2_InfiniteWire_MagnetismIsRelativistic();
	Demo3_TripleValidation_MovingCharge();
	Demo4_FieldVisualization_PancakeEffect();
	Demo5_FieldLineTracing_Export();
}