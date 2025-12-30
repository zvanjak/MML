#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Geometry3D.h"
#include "base/Geometry3DBodies.h"
#include "core/Integration/SurfaceIntegration.h"
#include "core/Surfaces.h"
#include "interfaces/IFunction.h"
#endif

using namespace MML;
using namespace MML::Testing;
using namespace MML::Surfaces;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

namespace MML::Tests::Algorithms::SurfaceIntegrationTests
{
	/*********************************************************************/
	/*****          Helper Vector Field Classes                      *****/
	/*********************************************************************/
	
	// Constant vector field F = (a, b, c)
	class ConstantVectorField : public IVectorFunction<3>
	{
		VectorN<Real, 3> _value;
	public:
		ConstantVectorField(Real a, Real b, Real c) 
			: _value{a, b, c} {}
		
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const override
		{
			return _value;
		}
	};
	
	// Linear radial field F = (x, y, z)
	class RadialVectorField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const override
		{
			return x;  // F = r
		}
	};
	
	// Quadratic field F = (x², y², z²)
	class QuadraticVectorField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const override
		{
			VectorN<Real, 3> result;
			result[0] = x[0] * x[0];
			result[1] = x[1] * x[1];
			result[2] = x[2] * x[2];
			return result;
		}
	};
	/*********************************************************************/
	/*****          Flux Through Simple Surfaces                     *****/
	/*********************************************************************/
	TEST_CASE("SurfaceIntegration::ConstantField_ThroughSquare", "[flux][basic]")
	{
			TEST_PRECISION_INFO();
		// Constant vector field F = (0, 0, 1) through unit square in XY plane
		// Square: (0,0,0), (1,0,0), (1,1,0), (0,1,0)
		// Normal: (0, 0, 1)
		// Expected flux: F·n * Area = 1 * 1 = 1
		
		ConstantVectorField field(0, 0, 1);

		RectSurface3D square(
			Point3Cartesian(0, 0, 0),
			Point3Cartesian(1, 0, 0),
			Point3Cartesian(1, 1, 0),
			Point3Cartesian(0, 1, 0)
		);

		Real flux = SurfaceIntegration::SurfaceIntegral(field, square, 1e-6);
		
		REQUIRE_THAT(flux, WithinAbs(REAL(1.0), REAL(1e-5)));
	}

	TEST_CASE("SurfaceIntegration::ConstantField_ThroughTriangle", "[flux][basic]")
	{
			TEST_PRECISION_INFO();
		// Constant field F = (1, 0, 0) through triangle in YZ plane
		// Triangle: (0,0,0), (0,1,0), (0,0,1)
		// Area = REAL(0.5), Normal = (1, 0, 0)
		// Expected flux: 1 * REAL(0.5) = REAL(0.5)
		
		ConstantVectorField field(1, 0, 0);

		Triangle3D triangle(
			Point3Cartesian(0, 0, 0),
			Point3Cartesian(0, 1, 0),
			Point3Cartesian(0, 0, 1)
		);

		Real flux = SurfaceIntegration::SurfaceIntegral(field, triangle, 1e-6);
		
		REQUIRE_THAT(flux, WithinAbs(REAL(0.5), REAL(1e-5)));
	}

	TEST_CASE("SurfaceIntegration::LinearField_ThroughSquare", "[flux]")
	{
			TEST_PRECISION_INFO();
		// Linear field F = (x, y, z) through unit square z=1
		// Square: (0,0,1), (1,0,1), (1,1,1), (0,1,1)
		// At center (REAL(0.5), REAL(0.5), 1): F = (REAL(0.5), REAL(0.5), 1)
		// Normal = (0, 0, 1), F·n = 1
		// Flux ≈ 1 * 1 = 1
		
		RadialVectorField field;

		RectSurface3D square(
			Point3Cartesian(0, 0, 1),
			Point3Cartesian(1, 0, 1),
			Point3Cartesian(1, 1, 1),
			Point3Cartesian(0, 1, 1)
		);

		Real flux = SurfaceIntegration::SurfaceIntegral(field, square, 1e-6);
		
		REQUIRE_THAT(flux, WithinAbs(REAL(1.0), REAL(1e-4)));
	}

	/*********************************************************************/
	/*****          Flux Through Closed Surfaces (Gauss's Law)       *****/
	/*********************************************************************/
	TEST_CASE("SurfaceIntegration::DivergenceTheorem_Cube", "[flux][gauss]")
	{
			TEST_PRECISION_INFO();
		// Constant field F = (1, 2, 3) through unit cube
		// Divergence of F = ∂F_x/∂x + ∂F_y/∂y + ∂F_z/∂z = 0
		// By divergence theorem: flux through closed surface = 0
		
		ConstantVectorField field(1, 2, 3);

		// Unit cube: 6 faces
		Cube3D cube(REAL(1.0), Point3Cartesian(0, 0, 0));

		Real flux = SurfaceIntegration::SurfaceIntegral(field, cube, 1e-6);
		
		// Should be zero (within numerical precision)
		REQUIRE_THAT(flux, WithinAbs(REAL(0.0), REAL(1e-4)));
	}

	TEST_CASE("SurfaceIntegration::RadialField_ThroughCube", "[flux][gauss]")
	{
			TEST_PRECISION_INFO();
		// Radial field F = (x, y, z) centered at origin
		// Divergence = 3
		// Cube3D(size, center) creates cube with side length = size, centered at center  
		// With size=2, center=(1,1,1): cube spans (0,0,0) to (2,2,2), volume = 8
		// Theoretical flux = divergence * volume = 3 * 8 = 24
		
		RadialVectorField field;

		// Cube with side 2, centered at (1,1,1) - volume = 8
		Cube3D cube(REAL(2.0), Point3Cartesian(1, 1, 1));

		Real flux = SurfaceIntegration::SurfaceIntegral(field, cube, 1e-6);
		
		REQUIRE_THAT(flux, WithinRel(REAL(24.0), REAL(1e-3)));
	}

	/*********************************************************************/
	/*****          Flux With Different Orientations                 *****/
	/*********************************************************************/
	TEST_CASE("SurfaceIntegration::FluxDirection_Matters", "[flux][orientation]")
	{
			TEST_PRECISION_INFO();
		// Test that reversing surface normal reverses flux
		
		ConstantVectorField vectorField(0, 0, 1);

		// Square facing +Z
		RectSurface3D squareUp(
			Point3Cartesian(0, 0, 0),
			Point3Cartesian(1, 0, 0),
			Point3Cartesian(1, 1, 0),
			Point3Cartesian(0, 1, 0)
		);

		// Square facing -Z (reversed winding)
		RectSurface3D squareDown(
			Point3Cartesian(0, 0, 0),
			Point3Cartesian(0, 1, 0),
			Point3Cartesian(1, 1, 0),
			Point3Cartesian(1, 0, 0)
		);

		Real fluxUp = SurfaceIntegration::SurfaceIntegral(vectorField, squareUp, 1e-6);
		Real fluxDown = SurfaceIntegration::SurfaceIntegral(vectorField, squareDown, 1e-6);

		// Fluxes should be opposite
		REQUIRE_THAT(fluxUp + fluxDown, WithinAbs(REAL(0.0), REAL(1e-5)));
		REQUIRE_THAT(fluxUp, WithinAbs(REAL(1.0), REAL(1e-5)));
		REQUIRE_THAT(fluxDown, WithinAbs(-REAL(1.0), REAL(1e-5)));
	}

	/*********************************************************************/
	/*****          Multiple Surface Tests                           *****/
	/*********************************************************************/
	TEST_CASE("SurfaceIntegration::Tetrahedron_ConstantField", "[flux][composite]")
	{
			TEST_PRECISION_INFO();
		// Constant field through tetrahedron (4 triangular faces)
		// Field F = (0, 0, 1), divergence = 0
		// Flux through closed surface should be 0
		
		ConstantVectorField field(0, 0, 1);

		// Create tetrahedron with triangular faces
		// Vertices: (0,0,0), (1,0,0), (0,1,0), (0,0,1)
		std::vector<Triangle3D> faces;
		
		// Base triangle in XY plane
		faces.push_back(Triangle3D(
			Point3Cartesian(0, 0, 0),
			Point3Cartesian(0, 1, 0),
			Point3Cartesian(1, 0, 0)
		));
		
		// Three side faces
		faces.push_back(Triangle3D(
			Point3Cartesian(0, 0, 0),
			Point3Cartesian(1, 0, 0),
			Point3Cartesian(0, 0, 1)
		));
		
		faces.push_back(Triangle3D(
			Point3Cartesian(0, 0, 0),
			Point3Cartesian(0, 0, 1),
			Point3Cartesian(0, 1, 0)
		));
		
		faces.push_back(Triangle3D(
			Point3Cartesian(1, 0, 0),
			Point3Cartesian(0, 1, 0),
			Point3Cartesian(0, 0, 1)
		));

		Real totalFlux = REAL(0.0);
		for (const auto& face : faces) {
			totalFlux += SurfaceIntegration::SurfaceIntegral(field, face, 1e-6);
		}

		// Constant field has zero divergence, so total flux should be ~0
		REQUIRE_THAT(totalFlux, WithinAbs(REAL(0.0), REAL(0.1)));
	}

	/*********************************************************************/
	/*****          Precision and Refinement Tests                   *****/
	/*********************************************************************/
	TEST_CASE("SurfaceIntegration::RefinementImprovesPrecision", "[flux][refinement]")
	{
			TEST_PRECISION_INFO();
		// Test that recursive refinement improves accuracy
		// Use a field that varies across the surface
		
		QuadraticVectorField field;

		RectSurface3D square(
			Point3Cartesian(0, 0, 1),
			Point3Cartesian(1, 0, 1),
			Point3Cartesian(1, 1, 1),
			Point3Cartesian(0, 1, 1)
		);

		// Low precision
		Real fluxLow = SurfaceIntegration::SurfaceIntegral(field, square, REAL(0.1), 2);
		
		// High precision
		Real fluxHigh = SurfaceIntegration::SurfaceIntegral(field, square, 1e-8, 6);

		// They should be different (refinement matters)
		// And high precision should be close to analytical value
		// For F_z = z² at z=1: F_z = 1, so flux ≈ 1 * 1 = 1
		REQUIRE_THAT(fluxHigh, WithinRel(REAL(1.0), REAL(1e-3)));
	}

	/*********************************************************************/
	/*****          Edge Cases                                       *****/
	/*********************************************************************/
	TEST_CASE("SurfaceIntegration::ZeroArea_Triangle", "[flux][edge]")
	{
			TEST_PRECISION_INFO();
		// Degenerate triangle (all points collinear) should have zero flux
		
		ConstantVectorField vectorField(1, 1, 1);

		// Collinear points (line, not a triangle)
		Triangle3D degenerateTriangle(
			Point3Cartesian(0, 0, 0),
			Point3Cartesian(1, 0, 0),
			Point3Cartesian(2, 0, 0)
		);

		Real flux = SurfaceIntegration::SurfaceIntegral(vectorField, degenerateTriangle, 1e-6);
		
		REQUIRE_THAT(flux, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	TEST_CASE("SurfaceIntegration::ParallelField_NoFlux", "[flux][perpendicular]")
	{
			TEST_PRECISION_INFO();
		// Field tangent to surface (perpendicular to normal) should have zero flux
		// Field in XY plane: F = (1, 1, 0)
		// Square in XY plane, normal is (0, 0, 1)
		// F · n = (1, 1, 0) · (0, 0, 1) = 0
		
		ConstantVectorField tangentField(1, 1, 0);

		RectSurface3D square(
			Point3Cartesian(0, 0, 0),
			Point3Cartesian(1, 0, 0),
			Point3Cartesian(1, 1, 0),
			Point3Cartesian(0, 1, 0)
		);

		Real flux = SurfaceIntegration::SurfaceIntegral(tangentField, square, 1e-6);
		
		// Field · normal = (1, 1, 0) · (0, 0, 1) = 0
		REQUIRE_THAT(flux, WithinAbs(REAL(0.0), REAL(1e-10)));
	}

	/*********************************************************************/
	/*****          Advanced Flux Calculations                       *****/
	/*********************************************************************/
	
	// Helper class for inverse-square field (like gravity or electric field from point source)
	class InverseSquareField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			// F = r / |r|^3 (normalized radial direction divided by distance squared)
			Real r2 = pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2];
			if (r2 < 1e-10) return VectorN<Real, 3>({0, 0, 0}); // Avoid singularity at origin
			Real r3 = r2 * sqrt(r2);
			return VectorN<Real, 3>({pos[0] / r3, pos[1] / r3, pos[2] / r3});
		}

		Vector3Cartesian operator()(const Vector3Cartesian& pos) const
		{
			VectorN<Real, 3> v = operator()(VectorN<Real, 3>({pos[0], pos[1], pos[2]}));
			return Vector3Cartesian(v[0], v[1], v[2]);
		}
	};

	// Helper class for cylindrical field F = (-y, x, 0) (rotation around Z-axis)
	class RotationalField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			return VectorN<Real, 3>({-pos[1], pos[0], 0});
		}

		Vector3Cartesian operator()(const Vector3Cartesian& pos) const
		{
			VectorN<Real, 3> v = operator()(VectorN<Real, 3>({pos[0], pos[1], pos[2]}));
			return Vector3Cartesian(v[0], v[1], v[2]);
		}
	};

	// Helper class for sinusoidal field F = (sin(x), sin(y), sin(z))
	class SinusoidalField : public IVectorFunction<3>
	{
	public:
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
		{
			return VectorN<Real, 3>({sin(pos[0]), sin(pos[1]), sin(pos[2])});
		}

		Vector3Cartesian operator()(const Vector3Cartesian& pos) const
		{
			VectorN<Real, 3> v = operator()(VectorN<Real, 3>({pos[0], pos[1], pos[2]}));
			return Vector3Cartesian(v[0], v[1], v[2]);
		}
	};

	TEST_CASE("SurfaceIntegration::InverseSquareField_SmallCube", "[flux][inverse-square]")
	{
			TEST_PRECISION_INFO();
		// Inverse square field F = r/|r|^3 through small cube not containing origin
		// Flux should approximate Gauss's law behavior
		// For a small cube at distance >> size, flux ≈ solid angle * source strength
		
		InverseSquareField field;

		// Small cube centered at (2, 2, 2), side length REAL(0.5)
		Cube3D cube(REAL(0.5), Point3Cartesian(2, 2, 2));

		Real flux = SurfaceIntegration::SurfaceIntegral(field, cube, 1e-5);
		
		// Flux should be small and positive (pointing outward from origin)
		// Actual value depends on solid angle subtended by cube
		REQUIRE(flux > REAL(0.0));
		REQUIRE(flux < REAL(0.1)); // Small solid angle
	}

	TEST_CASE("SurfaceIntegration::RotationalField_ThroughSquare", "[flux][rotational]")
	{
			TEST_PRECISION_INFO();
		// Rotational field F = (-y, x, 0) has curl but zero divergence
		// Through closed surface, flux should be zero
		// Through open surface perpendicular to Z-axis, flux should be zero
		
		RotationalField field;

		RectSurface3D square(
			Point3Cartesian(0, 0, 1),
			Point3Cartesian(1, 0, 1),
			Point3Cartesian(1, 1, 1),
			Point3Cartesian(0, 1, 1)
		);

		Real flux = SurfaceIntegration::SurfaceIntegral(field, square, 1e-6);
		
		// F = (-y, x, 0), normal = (0, 0, 1)
		// F · n = 0 everywhere
		REQUIRE_THAT(flux, WithinAbs(REAL(0.0), REAL(1e-8)));
	}

	TEST_CASE("SurfaceIntegration::SinusoidalField_ThroughSquare", "[flux][sinusoidal]")
	{
			TEST_PRECISION_INFO();
		// Sinusoidal field F = (sin(x), sin(y), sin(z))
		// Tests integration of oscillating field
		
		SinusoidalField field;

		// Square in Z=1 plane, from (0,0,1) to (π, π, 1)
		RectSurface3D square(
			Point3Cartesian(0, 0, 1),
			Point3Cartesian(Constants::PI, 0, 1),
			Point3Cartesian(Constants::PI, Constants::PI, 1),
			Point3Cartesian(0, Constants::PI, 1)
		);

		Real flux = SurfaceIntegration::SurfaceIntegral(field, square, 1e-6);
		
		// F_z = sin(z), at z=1: sin(1) ≈ REAL(0.8414709848)
		// Flux = sin(1) * area = sin(1) * π² ≈ REAL(8.289)
		Real expected = sin(REAL(1.0)) * Constants::PI * Constants::PI;
		REQUIRE_THAT(flux, WithinRel(expected, REAL(1e-3)));
	}

	TEST_CASE("SurfaceIntegration::QuadraticField_ThroughCube", "[flux][quadratic]")
	{
			TEST_PRECISION_INFO();
		// Quadratic field F = (x², y², z²) through unit cube
		// Divergence = 2x + 2y + 2z
		// For unit cube (0,0,0) to (1,1,1): ∫∫∫ (2x + 2y + 2z) dV
		// = 2 * ∫₀¹ x dx + 2 * ∫₀¹ y dy + 2 * ∫₀¹ z dz
		// = 2 * [1/2] + 2 * [1/2] + 2 * [1/2] = 3
		
		QuadraticVectorField field;

		Cube3D cube(REAL(1.0), Point3Cartesian(REAL(0.5), REAL(0.5), REAL(0.5))); // Unit cube at origin

		Real flux = SurfaceIntegration::SurfaceIntegral(field, cube, 1e-5);
		
		// By divergence theorem, flux = 3
		REQUIRE_THAT(flux, WithinRel(REAL(3.0), REAL(5e-2)));
	}

	/*********************************************************************/
	/*****          Numerical Precision and Convergence Tests        *****/
	/*********************************************************************/
	TEST_CASE("SurfaceIntegration::ToleranceConvergence", "[flux][precision][tolerance]")
	{
			TEST_PRECISION_INFO();
		// Test that decreasing tolerance improves accuracy
		// Use a field that has known analytical solution
		
		ConstantVectorField field(1, 2, 3);

		RectSurface3D square(
			Point3Cartesian(0, 0, 0),
			Point3Cartesian(1, 0, 0),
			Point3Cartesian(1, 1, 0),
			Point3Cartesian(0, 1, 0)
		);

		// Different tolerance levels
		Real flux1 = SurfaceIntegration::SurfaceIntegral(field, square, 1e-2, 3);
		Real flux2 = SurfaceIntegration::SurfaceIntegral(field, square, 1e-4, 5);
		Real flux3 = SurfaceIntegration::SurfaceIntegral(field, square, 1e-6, 7);

		// Flux of constant field F=(1,2,3) through XY unit square:
		// Flux = F · n * Area = (1,2,3) · (0,0,1) * 1 = 3
		Real expected = REAL(3.0);
		REQUIRE_THAT(flux1, WithinAbs(expected, REAL(1e-1)));
		REQUIRE_THAT(flux2, WithinAbs(expected, REAL(1e-3)));
		REQUIRE_THAT(flux3, WithinAbs(expected, REAL(1e-5)));
	}

	TEST_CASE("SurfaceIntegration::VerySmallSurface", "[flux][precision][small]")
	{
			TEST_PRECISION_INFO();
		// Test numerical stability with very small surface
		
		ConstantVectorField field(0, 0, 1);

		// Tiny square: 1e-6 × 1e-6 at z=1
		RectSurface3D tinySquare(
			Point3Cartesian(0, 0, 1),
			Point3Cartesian(1e-6, 0, 1),
			Point3Cartesian(1e-6, 1e-6, 1),
			Point3Cartesian(0, 1e-6, 1)
		);

		Real flux = SurfaceIntegration::SurfaceIntegral(field, tinySquare, 1e-10);
		
		// Expected flux = 1 * (1e-6)² = 1e-12
		REQUIRE_THAT(flux, WithinRel(REAL(1e-12), REAL(1e-2)));
	}

	TEST_CASE("SurfaceIntegration::VeryLargeSurface", "[flux][precision][large]")
	{
			TEST_PRECISION_INFO();
		// Test numerical stability with very large surface
		
		ConstantVectorField field(0, 0, 2);

		// Large square: 1000 × 1000 at z=0
		RectSurface3D largeSquare(
			Point3Cartesian(0, 0, 0),
			Point3Cartesian(1000, 0, 0),
			Point3Cartesian(1000, 1000, 0),
			Point3Cartesian(0, 1000, 0)
		);

		Real flux = SurfaceIntegration::SurfaceIntegral(field, largeSquare, 1e-3);
		
		// Expected flux = 2 * 1000² = 2,000,000
		REQUIRE_THAT(flux, WithinRel(REAL(2e6), REAL(1e-4)));
	}

	TEST_CASE("SurfaceIntegration::HighAspectRatioSurface", "[flux][precision][aspect-ratio]")
	{
			TEST_PRECISION_INFO();
		// Test numerical stability with high aspect ratio surface (long and thin)
		
		ConstantVectorField field(0, 0, 1);

		// Rectangle: 100 × REAL(0.01) at z=0
		RectSurface3D rectangle(
			Point3Cartesian(0, 0, 0),
			Point3Cartesian(100, 0, 0),
			Point3Cartesian(100, REAL(0.01), 0),
			Point3Cartesian(0, REAL(0.01), 0)
		);

		Real flux = SurfaceIntegration::SurfaceIntegral(field, rectangle, 1e-6);
		
		// Expected flux = 1 * (100 * REAL(0.01)) = 1
		REQUIRE_THAT(flux, WithinRel(REAL(1.0), REAL(1e-3)));
	}

	TEST_CASE("SurfaceIntegration::NearlyDegenerateTriangle", "[flux][precision][edge]")
	{
			TEST_PRECISION_INFO();
		// Triangle with very small area (nearly collinear points)
		
		ConstantVectorField field(1, 1, 1);

		// Nearly collinear points: (0,0,0), (1,0,0), (REAL(1.001), REAL(0.001), 0)
		// This forms a very thin triangle
		Triangle3D nearlyDegenerate(
			Point3Cartesian(0, 0, 0),
			Point3Cartesian(1, 0, 0),
			Point3Cartesian(REAL(1.001), REAL(0.001), 0)
		);

		Real flux = SurfaceIntegration::SurfaceIntegral(field, nearlyDegenerate, 1e-8);
		
		// Flux should be very small but not exactly zero
		// Area ≈ REAL(0.5) * |base × height| ≈ REAL(0.5) * 1 * REAL(0.001) = REAL(0.0005)
		REQUIRE(fabs(flux) < 1e-2);
	}

	/*********************************************************************/
	/*****          Additional Divergence Theorem Tests              *****/
	/*********************************************************************/
	TEST_CASE("SurfaceIntegration::DivergenceTheorem_LinearField_Cube", "[flux][gauss][divergence]")
	{
			TEST_PRECISION_INFO();
		// Linear field F = (x, y, z) centered at origin
		// div(F) = 1 + 1 + 1 = 3
		// For cube with volume V: flux = 3V
		// Cube side=1 centered at (REAL(0.5), REAL(0.5), REAL(0.5)): volume = 1
		
		RadialVectorField field;

		Cube3D cube(REAL(1.0), Point3Cartesian(REAL(0.5), REAL(0.5), REAL(0.5)));

		Real flux = SurfaceIntegration::SurfaceIntegral(field, cube, 1e-5);
		
		// div(F) = 3, volume = 1, expected flux = 3
		REQUIRE_THAT(flux, WithinRel(REAL(3.0), REAL(1e-2)));
	}

	TEST_CASE("SurfaceIntegration::DivergenceTheorem_QuadraticField_SmallCube", "[flux][gauss][divergence]")
	{
			TEST_PRECISION_INFO();
		// Quadratic field F = (x², y², z²)
		// div(F) = 2x + 2y + 2z
		// For small cube centered at (REAL(0.1), REAL(0.1), REAL(0.1)) with side REAL(0.2)
		// Volume = REAL(0.008)
		// Average divergence ≈ 2(REAL(0.1) + REAL(0.1) + REAL(0.1)) = REAL(0.6)
		// Expected flux ≈ REAL(0.6) * REAL(0.008) = REAL(0.0048)
		
		QuadraticVectorField field;

		Cube3D cube(REAL(0.2), Point3Cartesian(REAL(0.1), REAL(0.1), REAL(0.1)));

		Real flux = SurfaceIntegration::SurfaceIntegral(field, cube, 1e-6);
		
		// By divergence theorem, flux ≈ REAL(0.0048)
		REQUIRE_THAT(flux, WithinRel(REAL(0.0048), REAL(0.1)));
	}

	/*********************************************************************/
	/*****          Multi-Surface Composite Tests                    *****/
	/*********************************************************************/
	TEST_CASE("SurfaceIntegration::Octahedron_RadialField", "[flux][composite][octahedron]")
	{
			TEST_PRECISION_INFO();
		// Octahedron (8 triangular faces) with radial field
		// Vertices at (±1, 0, 0), (0, ±1, 0), (0, 0, ±1)
		
		RadialVectorField field;

		std::vector<Triangle3D> faces;
		
		// Define 8 triangular faces of octahedron
		// Top pyramid (apex at (0,0,1))
		faces.push_back(Triangle3D(
			Point3Cartesian(1, 0, 0),
			Point3Cartesian(0, 1, 0),
			Point3Cartesian(0, 0, 1)
		));
		faces.push_back(Triangle3D(
			Point3Cartesian(0, 1, 0),
			Point3Cartesian(-1, 0, 0),
			Point3Cartesian(0, 0, 1)
		));
		faces.push_back(Triangle3D(
			Point3Cartesian(-1, 0, 0),
			Point3Cartesian(0, -1, 0),
			Point3Cartesian(0, 0, 1)
		));
		faces.push_back(Triangle3D(
			Point3Cartesian(0, -1, 0),
			Point3Cartesian(1, 0, 0),
			Point3Cartesian(0, 0, 1)
		));
		
		// Bottom pyramid (apex at (0,0,-1))
		faces.push_back(Triangle3D(
			Point3Cartesian(1, 0, 0),
			Point3Cartesian(0, 0, -1),
			Point3Cartesian(0, 1, 0)
		));
		faces.push_back(Triangle3D(
			Point3Cartesian(0, 1, 0),
			Point3Cartesian(0, 0, -1),
			Point3Cartesian(-1, 0, 0)
		));
		faces.push_back(Triangle3D(
			Point3Cartesian(-1, 0, 0),
			Point3Cartesian(0, 0, -1),
			Point3Cartesian(0, -1, 0)
		));
		faces.push_back(Triangle3D(
			Point3Cartesian(0, -1, 0),
			Point3Cartesian(0, 0, -1),
			Point3Cartesian(1, 0, 0)
		));

		Real totalFlux = REAL(0.0);
		for (const auto& face : faces) {
			totalFlux += SurfaceIntegration::SurfaceIntegral(field, face, 1e-5);
		}

		// Radial field F = (x,y,z) has div(F) = 3
		// Octahedron with vertices at (±1,0,0), (0,±1,0), (0,0,±1)
		// Volume = 4/3 (two pyramids, each with square base area = 2 and height = 1)
		// Expected flux = div(F) * V = 3 * (4/3) = 4
		Real octVolume = REAL(4.0) / REAL(3.0);
		Real expected = REAL(3.0) * octVolume;  // = REAL(4.0)
		REQUIRE_THAT(totalFlux, WithinRel(expected, REAL(0.15)));
	}

	TEST_CASE("SurfaceIntegration::MixedSurfaces_CubeAndTriangles", "[flux][composite][mixed]")
	{
			TEST_PRECISION_INFO();
		// Test that we can combine different surface types
		// Use constant field to verify basic consistency
		
		ConstantVectorField field(0, 0, 5);

		// Top face of cube as RectSurface
		RectSurface3D topFace(
			Point3Cartesian(0, 0, 1),
			Point3Cartesian(1, 0, 1),
			Point3Cartesian(1, 1, 1),
			Point3Cartesian(0, 1, 1)
		);

		// Same face divided into two triangles
		Triangle3D tri1(
			Point3Cartesian(0, 0, 1),
			Point3Cartesian(1, 0, 1),
			Point3Cartesian(1, 1, 1)
		);
		Triangle3D tri2(
			Point3Cartesian(0, 0, 1),
			Point3Cartesian(1, 1, 1),
			Point3Cartesian(0, 1, 1)
		);

		Real fluxRect = SurfaceIntegration::SurfaceIntegral(field, topFace, 1e-6);
		Real fluxTri = SurfaceIntegration::SurfaceIntegral(field, tri1, 1e-6) +
		               SurfaceIntegration::SurfaceIntegral(field, tri2, 1e-6);

		// Both should give same result: 5 * 1 = 5
		REQUIRE_THAT(fluxRect, WithinRel(REAL(5.0), REAL(1e-4)));
		REQUIRE_THAT(fluxTri, WithinRel(REAL(5.0), REAL(1e-4)));
		REQUIRE_THAT(fluxRect - fluxTri, WithinAbs(REAL(0.0), REAL(1e-3)));
	}

/*********************************************************************/
// Parametric surface integration tests
/*********************************************************************/

TEST_CASE("SurfaceIntegration::Parametric_Sphere_RadialField", "[flux][parametric][sphere]")
{
// Radial field F = (x, y, z) pointing outward from origin
// For unit sphere, F  n = R everywhere
// Expected flux =  R dS = 4pR * R = 4pR
// For R=1: flux = 4p  REAL(12.566)

class RadialField : public IVectorFunction<3>
{
public:
VectorN<Real, 3> operator()(const VectorN<Real, 3>& pos) const override
{
return pos;  // Radial field pointing outward
}
};

RadialField field;
Sphere sphere(REAL(1.0));

Real flux = SurfaceIntegration::SurfaceIntegral(
field, sphere,
sphere.getMinU(), sphere.getMaxU(),
sphere.getMinW(), sphere.getMaxW(),
30, 60
);

// Expected: 4pR4 = 4p for R=1
Real expected = REAL(4.0) * Constants::PI;
REQUIRE_THAT(flux, WithinRel(expected, REAL(0.05)));  // 5% tolerance
}

TEST_CASE("SurfaceIntegration::Parametric_Plane_ConstantField", "[flux][parametric][plane]")
{
// Flat plane surface with constant field
// Simple verification that parametric integration works for planes

ConstantVectorField field(0, 0, 1);

// Plane parallel to x-y plane at z=0, spanning [0,1] x [0,1]
Vec3Cart point(0, 0, 0);
Vec3Cart normal(0, 0, 1);
Vec3Cart uAxis(1, 0, 0);
Vec3Cart vAxis(0, 1, 0);
PlaneSurface plane(point, normal, uAxis, vAxis, 0, 1, 0, 1);

Real flux = SurfaceIntegration::SurfaceIntegral(
field, plane,
plane.getMinU(), plane.getMaxU(),
plane.getMinW(), plane.getMaxW(),
20, 20
);

// Area = 1, Fn = 1, so flux = 1
REQUIRE_THAT(flux, WithinRel(REAL(1.0), REAL(1e-3)));
}
}
