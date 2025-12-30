#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Curves.h"
#include "core/Surfaces.h"
#include "core/Integration/PathIntegration.h"
#endif

#include "../test_data/parametric_curves_test_bed.h"

using namespace MML;
using namespace MML::Curves;
using namespace MML::Surfaces;

///////////////////////////////////////////////////////////////////////////////////////////
/////                              2D CURVES DEMOS                                    /////
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Curves_2D_Basic()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****          2D Curves - Basic Predefined Curves                 ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// Circle2D
	std::cout << "\n--- Circle2DCurve ---" << std::endl;
	Circle2DCurve circle(5.0, Pnt2Cart(2.0, 3.0));  // Radius 5 at center (2, 3)
	
	Vec2Cart pt_top = circle(Constants::PI / 2);      // Top point
	Vec2Cart pt_right = circle(0.0);                  // Right point
	Vec2Cart tangent = circle.getTangent(0.0);        // Tangent at t=0
	Vec2Cart normal = circle.getNormalUnit(0.0);      // Normal at t=0
	
	std::cout << "Circle: radius=5, center=(2, 3)" << std::endl;
	std::cout << "  Point at t=PI/2: "; pt_top.Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "  Point at t=0:    "; pt_right.Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "  Tangent at t=0:  "; tangent.Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "  Normal at t=0:   "; normal.Print(std::cout, 8, 4); std::cout << std::endl;

	// Logarithmic Spiral
	std::cout << "\n--- LogSpiralCurve ---" << std::endl;
	LogSpiralCurve spiral(-0.1, 2.0);  // Slow inward spiral
	
	std::cout << "Logarithmic Spiral: lambda=-0.1, c=2" << std::endl;
	std::cout << "  Points along spiral:" << std::endl;
	for (int i = 0; i <= 4; i++) {
		Real t = i * 2.5;
		Vec2Cart pt = spiral(t);
		std::cout << "    t=" << std::setw(4) << t << ": "; pt.Print(std::cout, 10, 4); std::cout << std::endl;
	}

	// Lemniscate
	std::cout << "\n--- LemniscateCurve (Figure-8) ---" << std::endl;
	LemniscateCurve lemniscate;
	
	std::cout << "  Points along lemniscate:" << std::endl;
	for (int i = 0; i <= 4; i++) {
		Real t = i * Constants::PI / 4;
		Vec2Cart pt = lemniscate(t);
		std::cout << "    t=" << std::setw(6) << std::fixed << std::setprecision(4) << t << ": "; 
		pt.Print(std::cout, 10, 4); std::cout << std::endl;
	}

	// Deltoid
	std::cout << "\n--- DeltoidCurve ---" << std::endl;
	DeltoidCurve deltoid(2);  // Scaled deltoid
	
	std::cout << "Deltoid: n=2" << std::endl;
	Vec2Cart deltoid_pt = deltoid(Constants::PI / 3);
	std::cout << "  Point at t=PI/3: "; deltoid_pt.Print(std::cout, 10, 4); std::cout << std::endl;

	// Astroid
	std::cout << "\n--- AstroidCurve ---" << std::endl;
	AstroidCurve astroid(3.0);  // c=3
	
	std::cout << "Astroid: c=3" << std::endl;
	for (int i = 0; i <= 4; i++) {
		Real t = i * Constants::PI / 4;
		Vec2Cart pt = astroid(t);
		std::cout << "  t=" << std::setw(6) << std::fixed << std::setprecision(4) << t << ": "; 
		pt.Print(std::cout, 10, 4); std::cout << std::endl;
	}

	// Archimedean Spiral
	std::cout << "\n--- ArchimedeanSpiralCurve ---" << std::endl;
	ArchimedeanSpiralCurve archspiral(0.5);
	
	std::cout << "Archimedean Spiral: a=0.5" << std::endl;
	for (int i = 0; i <= 3; i++) {
		Real t = i * Constants::PI;
		Vec2Cart pt = archspiral(t);
		std::cout << "  t=" << std::setw(6) << std::fixed << std::setprecision(4) << t << ": "; 
		pt.Print(std::cout, 10, 4); std::cout << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
/////                              3D CURVES DEMOS                                    /////
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Curves_3D_Basic()
{
	std::cout << "\n***********************************************************************" << std::endl;
	std::cout << "****          3D Curves - Basic Predefined Curves                 ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// Line curve
	std::cout << "\n--- LineCurve ---" << std::endl;
	LineCurve line(0.0, 10.0, Point3Cartesian(0, 0, 0), Vector3Cartesian(1, 1, 1));
	
	std::cout << "Line: from (0,0,0) in direction (1,1,1), t in [0, 10]" << std::endl;
	Vec3Cart line_pt5 = line(5.0);
	std::cout << "  Point at t=5: "; line_pt5.Print(std::cout, 8, 4); std::cout << std::endl;
	Vec3Cart line_tang = line.getTangentUnit(5.0);
	std::cout << "  Unit tangent: "; line_tang.Print(std::cout, 8, 4); std::cout << std::endl;

	// Circles in coordinate planes
	std::cout << "\n--- Circle3DXY, Circle3DXZ, Circle3DYZ ---" << std::endl;
	Circle3DXY circleXY(3.0);
	Circle3DXZ circleXZ(2.0);
	Circle3DYZ circleYZ(4.0);
	
	std::cout << "CircleXY (R=3) at t=PI: "; circleXY(Constants::PI).Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "CircleXZ (R=2) at t=PI: "; circleXZ(Constants::PI).Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "CircleYZ (R=4) at t=PI: "; circleYZ(Constants::PI).Print(std::cout, 8, 4); std::cout << std::endl;

	// General Circle in arbitrary plane
	std::cout << "\n--- Circle (General 3D) ---" << std::endl;
	Vec3Cart normal(1, 1, 1);  // Will be normalized
	Pnt3Cart center(0, 0, 5);
	Circle circle3D(2.0, normal, center);
	
	std::cout << "Circle: R=2, normal=(1,1,1), center=(0,0,5)" << std::endl;
	Vec3Cart pt_circle3D = circle3D(Constants::PI / 4);
	std::cout << "  Point at t=PI/4: "; pt_circle3D.Print(std::cout, 8, 4); std::cout << std::endl;

	// Helix curve
	std::cout << "\n--- HelixCurve ---" << std::endl;
	HelixCurve helix(2.0, 0.5);  // Radius 2, pitch parameter 0.5
	
	std::cout << "Helix: R=2, b=0.5" << std::endl;
	std::cout << "  Point at t=0: "; helix(0.0).Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "  Point at t=PI: "; helix(Constants::PI).Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "  Point at t=2PI: "; helix(2*Constants::PI).Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "  Curvature (constant): " << helix.getCurvature(0.0) << std::endl;
	std::cout << "  Torsion (constant):   " << helix.getTorsion(0.0) << std::endl;

	// Twisted cubic
	std::cout << "\n--- TwistedCubicCurve ---" << std::endl;
	TwistedCubicCurve twisted;
	
	std::cout << "Twisted Cubic: (t, t^2, t^3)" << std::endl;
	for (int i = -2; i <= 2; i++) {
		Real t = static_cast<Real>(i);
		Vec3Cart pt = twisted(t);
		std::cout << "  t=" << std::setw(2) << i << ": "; pt.Print(std::cout, 8, 4); std::cout << std::endl;
	}

	// Toroidal spiral
	std::cout << "\n--- ToroidalSpiralCurve ---" << std::endl;
	ToroidalSpiralCurve toroidal(3, 1.0);  // 3 windings
	
	std::cout << "Toroidal Spiral: n=3, scale=1.0" << std::endl;
	for (int i = 0; i <= 4; i++) {
		Real t = i * Constants::PI / 2;
		Vec3Cart pt = toroidal(t);
		std::cout << "  t=" << std::setw(6) << std::fixed << std::setprecision(4) << t << ": "; 
		pt.Print(std::cout, 8, 4); std::cout << std::endl;
	}
}

void Docs_Demo_Curves_Frenet_Frame()
{
	std::cout << "\n***********************************************************************" << std::endl;
	std::cout << "****          Frenet Frame - Tangent, Normal, Binormal            ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// Using Helix curve from TestData
	const CurveCartesian3D& helix = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;
	
	std::cout << "\nHelix curve - Frenet frame along the curve:" << std::endl;
	std::cout << "t         tangent                  unit_tangent              normal                   unit_normal                binormal" << std::endl;
	
	for (int i = -5; i <= 5; ++i)
	{
		Real t = i / 5.0;

		auto tangent = helix.getTangent(t);
		auto unit_tang = helix.getTangentUnit(t);
		auto normal = helix.getNormal(t);
		auto unit_norm = helix.getNormalUnit(t);
		auto binormal = helix.getBinormal(t);

		std::cout << "t=" << std::setw(4) << t << " : ";
		tangent.Print(std::cout, 9, 4, 1e-10); std::cout << "  ";
		unit_tang.Print(std::cout, 9, 4, 1e-10); std::cout << "  ";
		normal.Print(std::cout, 9, 4, 1e-10); std::cout << "  ";
		unit_norm.Print(std::cout, 9, 4, 1e-10); std::cout << "  ";
		binormal.Print(std::cout, 9, 4, 1e-10);
		std::cout << std::endl;
	}

	// Moving trihedron demonstration
	std::cout << "\n--- Moving Trihedron at specific point ---" << std::endl;
	HelixCurve myHelix(1.0, 1.0);
	Vector3Cartesian T, N, B;
	myHelix.getMovingTrihedron(Constants::PI/4, T, N, B);
	
	std::cout << "Helix at t=PI/4:" << std::endl;
	std::cout << "  Tangent (T):  "; T.Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "  Normal  (N):  "; N.Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "  Binormal (B): "; B.Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "  T·N = " << ScalarProduct(T, N) << " (should be ~0)" << std::endl;
	std::cout << "  T·B = " << ScalarProduct(T, B) << " (should be ~0)" << std::endl;
	std::cout << "  N·B = " << ScalarProduct(N, B) << " (should be ~0)" << std::endl;
}

void Docs_Demo_Curves_Curvature_Torsion()
{
	std::cout << "\n***********************************************************************" << std::endl;
	std::cout << "****          Curvature and Torsion                               ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// Helix - constant curvature and torsion
	std::cout << "\n--- Helix (constant curvature and torsion) ---" << std::endl;
	HelixCurve helix(2.0, 1.0);
	
	std::cout << "Helix: R=2, b=1" << std::endl;
	std::cout << "  Analytical curvature: R/(R² + b²) = 2/5 = " << 2.0/(4.0+1.0) << std::endl;
	std::cout << "  Analytical torsion:   b/(R² + b²) = 1/5 = " << 1.0/(4.0+1.0) << std::endl;
	std::cout << "  Computed curvature at t=0:   " << helix.getCurvature(0.0) << std::endl;
	std::cout << "  Computed torsion at t=0:     " << helix.getTorsion(0.0) << std::endl;
	std::cout << "  Computed curvature at t=PI:  " << helix.getCurvature(Constants::PI) << std::endl;
	std::cout << "  Computed torsion at t=PI:    " << helix.getTorsion(Constants::PI) << std::endl;

	// Circle - zero torsion, constant curvature
	std::cout << "\n--- Circle (zero torsion, constant curvature) ---" << std::endl;
	Circle3DXY circle(3.0);  // Radius 3
	
	std::cout << "Circle: R=3" << std::endl;
	std::cout << "  Expected curvature: 1/R = " << 1.0/3.0 << std::endl;
	std::cout << "  Expected torsion: 0" << std::endl;
	std::cout << "  Computed curvature at t=0:   " << circle.getCurvature(0.0) << std::endl;
	std::cout << "  Computed curvature at t=PI:  " << circle.getCurvature(Constants::PI) << std::endl;
	std::cout << "  Computed torsion at t=0:     " << circle.getTorsion(0.0) << std::endl;

	// Twisted cubic - varying curvature and torsion
	std::cout << "\n--- Twisted Cubic (varying curvature and torsion) ---" << std::endl;
	TwistedCubicCurve twisted;
	
	std::cout << "Twisted Cubic: (t, t², t³)" << std::endl;
	for (int i = 1; i <= 5; i++) {
		Real t = i * 0.5;
		std::cout << "  t=" << std::setw(3) << t 
				  << ": curvature=" << std::setw(10) << twisted.getCurvature(t)
				  << ", torsion=" << std::setw(10) << twisted.getTorsion(t) << std::endl;
	}

	// Curvature vector
	std::cout << "\n--- Curvature Vector ---" << std::endl;
	const CurveCartesian3D& testHelix = TestBeds::ParametricCurvesTestBed::getTestCurve("Helix")._curve;
	Vec3Cart curv_vec = testHelix.getCurvatureVector(0.0);
	std::cout << "Helix curvature vector at t=0: "; curv_vec.Print(std::cout, 10, 6); std::cout << std::endl;
}

void Docs_Demo_Curves_Arc_Length()
{
	std::cout << "\n***********************************************************************" << std::endl;
	std::cout << "****          Arc Length Calculation                              ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// Circle arc length
	std::cout << "\n--- Circle Arc Length ---" << std::endl;
	Circle3DXY circle(2.0);  // Radius 2
	
	Real circle_full = PathIntegration::ParametricCurveLength(circle, 0.0, 2*Constants::PI);
	Real circle_half = PathIntegration::ParametricCurveLength(circle, 0.0, Constants::PI);
	
	std::cout << "Circle R=2:" << std::endl;
	std::cout << "  Full circumference (expected 4*PI=" << 4*Constants::PI << "): " << circle_full << std::endl;
	std::cout << "  Half circumference (expected 2*PI=" << 2*Constants::PI << "): " << circle_half << std::endl;

	// Helix arc length
	std::cout << "\n--- Helix Arc Length ---" << std::endl;
	HelixCurve helix(1.0, 1.0);  // R=1, b=1
	
	// One full turn: length = 2*PI*sqrt(R² + b²) = 2*PI*sqrt(2)
	Real helix_turn = PathIntegration::ParametricCurveLength(helix, 0.0, 2*Constants::PI);
	Real expected = 2*Constants::PI*std::sqrt(2.0);
	
	std::cout << "Helix R=1, b=1 (one full turn):" << std::endl;
	std::cout << "  Expected: 2*PI*sqrt(2) = " << expected << std::endl;
	std::cout << "  Computed: " << helix_turn << std::endl;

	// Line segment arc length
	std::cout << "\n--- Line Segment Arc Length ---" << std::endl;
	LineCurve line(0.0, 10.0, Point3Cartesian(0, 0, 0), Vector3Cartesian(3, 4, 0));
	
	Real line_len = PathIntegration::ParametricCurveLength(line, 0.0, 2.0);
	std::cout << "Line from (0,0,0) direction (3,4,0), t in [0, 2]:" << std::endl;
	std::cout << "  Expected: 2*||(3,4,0)|| = 2*5 = 10" << std::endl;
	std::cout << "  Computed: " << line_len << std::endl;
}

void Docs_Demo_Curves_Planes()
{
	std::cout << "\n***********************************************************************" << std::endl;
	std::cout << "****          Osculating, Normal, Rectifying Planes               ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	HelixCurve helix(2.0, 0.5);
	Real t = Constants::PI / 4;
	
	std::cout << "Helix R=2, b=0.5 at t=PI/4:" << std::endl;
	
	// Get the point on the curve
	Vec3Cart pt = helix(t);
	std::cout << "\n  Point: "; pt.Print(std::cout, 8, 4); std::cout << std::endl;
	
	// Osculating plane (contains T and N)
	Plane3D oscPlane = helix.getOsculationPlane(t);
	std::cout << "\n  Osculating plane (contains curve locally):" << std::endl;
	std::cout << "    Normal: "; oscPlane.Normal().Print(std::cout, 8, 4); std::cout << std::endl;
	
	// Normal plane (perpendicular to T)
	Plane3D normPlane = helix.getNormalPlane(t);
	std::cout << "\n  Normal plane (perpendicular to tangent):" << std::endl;
	std::cout << "    Normal: "; normPlane.Normal().Print(std::cout, 8, 4); std::cout << std::endl;
	
	// Rectifying plane (contains T and B)
	Plane3D rectPlane = helix.getRectifyingPlane(t);
	std::cout << "\n  Rectifying plane (contains tangent and binormal):" << std::endl;
	std::cout << "    Normal: "; rectPlane.Normal().Print(std::cout, 8, 4); std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
/////                              SURFACES DEMOS                                     /////
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Surfaces_Basic()
{
	std::cout << "\n***********************************************************************" << std::endl;
	std::cout << "****          Basic Surface Properties                            ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// Sphere
	std::cout << "\n--- Sphere ---" << std::endl;
	Sphere sphere(5.0);  // Radius 5
	
	Real u = Constants::PI / 2;  // Equator
	Real w = 0.0;
	
	std::cout << "Sphere R=5:" << std::endl;
	Vec3Cart pt = sphere(u, w);
	std::cout << "  Point at (u=PI/2, w=0): "; pt.Print(std::cout, 8, 4); std::cout << std::endl;
	
	Vec3Cart n = sphere.Normal(u, w);
	std::cout << "  Normal at equator: "; n.Print(std::cout, 8, 4); std::cout << std::endl;
	
	std::cout << "  Gaussian curvature K = 1/R² = " << sphere.GaussianCurvature(u, w) 
			  << " (expected " << 1.0/25.0 << ")" << std::endl;
	std::cout << "  Mean curvature H = 1/R = " << sphere.MeanCurvature(u, w) 
			  << " (expected " << 1.0/5.0 << ")" << std::endl;

	// Cylinder
	std::cout << "\n--- Cylinder ---" << std::endl;
	Cylinder cylinder(2.0, 5.0);  // R=2, H=5
	
	u = Constants::PI / 4;
	w = 2.5;  // Middle height
	
	std::cout << "Cylinder R=2, H=5:" << std::endl;
	pt = cylinder(u, w);
	std::cout << "  Point at (u=PI/4, w=2.5): "; pt.Print(std::cout, 8, 4); std::cout << std::endl;
	
	n = cylinder.Normal(u, w);
	std::cout << "  Normal: "; n.Print(std::cout, 8, 4); std::cout << std::endl;
	
	std::cout << "  Gaussian curvature K = 0: " << cylinder.GaussianCurvature(u, w) << std::endl;
	std::cout << "  Mean curvature H = 1/(2R) = " << cylinder.MeanCurvature(u, w) 
			  << " (expected " << 1.0/4.0 << ")" << std::endl;

	// Torus
	std::cout << "\n--- Torus ---" << std::endl;
	Torus torus(3.0, 1.0);  // Major R=3, minor r=1
	
	std::cout << "Torus R=3, r=1:" << std::endl;
	
	// Outer point (u=0, w=0)
	pt = torus(0.0, 0.0);
	std::cout << "  Outer point (u=0, w=0): "; pt.Print(std::cout, 8, 4); std::cout << std::endl;
	
	// Inner point (u=0, w=PI)
	pt = torus(0.0, Constants::PI);
	std::cout << "  Inner point (u=0, w=PI): "; pt.Print(std::cout, 8, 4); std::cout << std::endl;
	
	// Curvatures at outer rim (positive K)
	u = 0.0; w = 0.0;
	std::cout << "  Gaussian curvature at outer rim: " << torus.GaussianCurvature(u, w) 
			  << " (positive - elliptic)" << std::endl;
	
	// Curvatures at inner rim (negative K)
	u = 0.0; w = Constants::PI;
	std::cout << "  Gaussian curvature at inner rim: " << torus.GaussianCurvature(u, w) 
			  << " (negative - hyperbolic)" << std::endl;
}

void Docs_Demo_Surfaces_Curvature()
{
	std::cout << "\n***********************************************************************" << std::endl;
	std::cout << "****          Surface Curvatures - Gaussian and Mean              ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// Ellipsoid - varying curvature
	std::cout << "\n--- Ellipsoid (varying curvature) ---" << std::endl;
	Ellipsoid ellipsoid(3.0, 2.0, 1.0);  // Flattened
	
	std::cout << "Ellipsoid a=3, b=2, c=1:" << std::endl;
	
	// At pole (u=0)
	std::cout << "  At pole (u=0):" << std::endl;
	std::cout << "    Gaussian K = " << ellipsoid.GaussianCurvature(0.01, 0.0) << std::endl;
	std::cout << "    Mean H = " << ellipsoid.MeanCurvature(0.01, 0.0) << std::endl;
	
	// At equator (u=PI/2)
	std::cout << "  At equator (u=PI/2, w=0):" << std::endl;
	std::cout << "    Point: "; ellipsoid(Constants::PI/2, 0.0).Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "    Gaussian K = " << ellipsoid.GaussianCurvature(Constants::PI/2, 0.0) << std::endl;
	std::cout << "    Mean H = " << ellipsoid.MeanCurvature(Constants::PI/2, 0.0) << std::endl;

	// Hyperboloid - negative Gaussian curvature
	std::cout << "\n--- Hyperboloid (negative Gaussian curvature) ---" << std::endl;
	Hyperboloid hyperboloid(1.0, 1.0, 1.0);
	
	std::cout << "Hyperboloid a=b=c=1:" << std::endl;
	std::cout << "  At u=0.5, w=0:" << std::endl;
	std::cout << "    Point: "; hyperboloid(0.5, 0.0).Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "    Gaussian K = " << hyperboloid.GaussianCurvature(0.5, 0.0) << " (negative - saddle)" << std::endl;

	// Monkey Saddle - saddle point with K=0
	std::cout << "\n--- MonkeySaddle ---" << std::endl;
	MonkeySaddle saddle;
	
	std::cout << "Monkey Saddle z = u(u² - 3w²):" << std::endl;
	std::cout << "  At origin (u=0, w=0):" << std::endl;
	std::cout << "    Point: "; saddle(0.0, 0.0).Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "    Gaussian K = " << saddle.GaussianCurvature(0.001, 0.001) << std::endl;
	std::cout << "  Away from origin (u=2, w=1):" << std::endl;
	std::cout << "    Point: "; saddle(2.0, 1.0).Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "    Gaussian K = " << saddle.GaussianCurvature(2.0, 1.0) << std::endl;
}

void Docs_Demo_Surfaces_Principal_Curvatures()
{
	std::cout << "\n***********************************************************************" << std::endl;
	std::cout << "****          Principal Curvatures and Directions                 ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// Sphere - all points umbilical (k1 = k2)
	std::cout << "\n--- Sphere (umbilical - k1 = k2) ---" << std::endl;
	Sphere sphere(4.0);
	
	Real k1, k2;
	sphere.PrincipalCurvatures(Constants::PI/2, 0.0, k1, k2);
	
	std::cout << "Sphere R=4 at equator:" << std::endl;
	std::cout << "  Principal curvatures: k1=" << k1 << ", k2=" << k2 << std::endl;
	std::cout << "  Expected: both = 1/R = " << 1.0/4.0 << std::endl;

	// Torus - different principal curvatures
	std::cout << "\n--- Torus (varying principal curvatures) ---" << std::endl;
	Torus torus(3.0, 1.0);
	
	// Outer equator
	torus.PrincipalCurvatures(0.0, 0.0, k1, k2);
	std::cout << "Torus R=3, r=1 at outer equator:" << std::endl;
	std::cout << "  Principal curvatures: k1=" << k1 << ", k2=" << k2 << std::endl;
	std::cout << "  K = k1*k2 = " << k1*k2 << std::endl;
	std::cout << "  H = (k1+k2)/2 = " << (k1+k2)/2 << std::endl;

	// Principal directions
	std::cout << "\n--- Principal Directions ---" << std::endl;
	Ellipsoid ellipsoid(3.0, 2.0, 1.0);
	
	Vec3Cart dir1, dir2;
	ellipsoid.PrincipalDirections(Constants::PI/4, Constants::PI/4, dir1, dir2);
	ellipsoid.PrincipalCurvatures(Constants::PI/4, Constants::PI/4, k1, k2);
	
	std::cout << "Ellipsoid at (PI/4, PI/4):" << std::endl;
	std::cout << "  Principal curvatures: k1=" << k1 << ", k2=" << k2 << std::endl;
	std::cout << "  Direction 1: "; dir1.Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "  Direction 2: "; dir2.Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "  dir1·dir2 = " << ScalarProduct(dir1, dir2) << " (should be ~0)" << std::endl;
}

void Docs_Demo_Surfaces_Fundamental_Forms()
{
	std::cout << "\n***********************************************************************" << std::endl;
	std::cout << "****          First and Second Fundamental Forms                  ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// First fundamental form (metric)
	std::cout << "\n--- First Fundamental Form (E, F, G) ---" << std::endl;
	
	Sphere sphere(3.0);
	Real E, F, G;
	sphere.GetFirstNormalFormCoefficients(Constants::PI/2, 0.0, E, F, G);
	
	std::cout << "Sphere R=3 at equator (u=PI/2, w=0):" << std::endl;
	std::cout << "  E = " << E << " (r_u · r_u)" << std::endl;
	std::cout << "  F = " << F << " (r_u · r_w, should be 0 for orthogonal params)" << std::endl;
	std::cout << "  G = " << G << " (r_w · r_w)" << std::endl;
	std::cout << "  Area element dA = sqrt(EG - F²) = " << std::sqrt(E*G - F*F) << std::endl;

	// Torus
	Torus torus(3.0, 1.0);
	torus.GetFirstNormalFormCoefficients(Constants::PI/4, Constants::PI/4, E, F, G);
	
	std::cout << "\nTorus R=3, r=1 at (PI/4, PI/4):" << std::endl;
	std::cout << "  E = " << E << std::endl;
	std::cout << "  F = " << F << std::endl;
	std::cout << "  G = " << G << std::endl;

	// Second fundamental form (shape operator)
	std::cout << "\n--- Second Fundamental Form (L, M, N) ---" << std::endl;
	
	Real L, M, N;
	sphere.GetSecondNormalFormCoefficients(Constants::PI/2, 0.0, L, M, N);
	
	std::cout << "Sphere R=3 at equator:" << std::endl;
	std::cout << "  L = " << L << " (r_uu · n)" << std::endl;
	std::cout << "  M = " << M << " (r_uw · n)" << std::endl;
	std::cout << "  N = " << N << " (r_ww · n)" << std::endl;
	
	std::cout << "\n  From fundamental forms:" << std::endl;
	std::cout << "    Gaussian K = (LN - M²)/(EG - F²) = " 
			  << (L*N - M*M)/(E*G - F*F) << std::endl;
	std::cout << "    Mean H = (EN + GL - 2FM)/(2(EG - F²)) = " 
			  << (E*N + G*L - 2*F*M)/(2*(E*G - F*F)) << std::endl;
}

void Docs_Demo_Surfaces_Special()
{
	std::cout << "\n***********************************************************************" << std::endl;
	std::cout << "****          Special Surfaces                                    ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// Möbius Strip
	std::cout << "\n--- MobiusStrip (non-orientable) ---" << std::endl;
	MobiusStrip mobius;
	
	std::cout << "Möbius Strip:" << std::endl;
	std::cout << "  Parameterization: u in [0, 2PI], w in [-1, 1]" << std::endl;
	
	// Points along center line (w=0)
	std::cout << "  Center line (w=0):" << std::endl;
	for (int i = 0; i <= 4; i++) {
		Real u = i * Constants::PI / 2;
		Vec3Cart pt = mobius(u, 0.0);
		std::cout << "    u=" << std::setw(6) << std::fixed << std::setprecision(4) << u << ": ";
		pt.Print(std::cout, 8, 4); std::cout << std::endl;
	}
	
	// Normal reversal check
	Vec3Cart n1 = mobius.Normal(0.0, 0.0);
	Vec3Cart n2 = mobius.Normal(2*Constants::PI - 0.01, 0.0);  // Just before completing loop
	std::cout << "\n  Normal at u=0: "; n1.Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "  Normal at u≈2PI: "; n2.Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "  (Normals should point opposite - non-orientable!)" << std::endl;

	// Paraboloid
	std::cout << "\n--- Paraboloid ---" << std::endl;
	Paraboloid paraboloid(1.0, 1.0);
	
	std::cout << "Paraboloid a=1, h=1:" << std::endl;
	Vec3Cart pt = paraboloid(4.0, 0.0);
	std::cout << "  Point at (u=4, w=0): "; pt.Print(std::cout, 8, 4); std::cout << std::endl;
	
	// Plane surface
	std::cout << "\n--- PlaneSurface ---" << std::endl;
	PlaneSurface plane(Vec3Cart(0, 0, 0), Vec3Cart(0, 0, 1));  // XY plane through origin
	
	std::cout << "Plane through origin with normal (0, 0, 1):" << std::endl;
	pt = plane(3.0, 2.0);
	std::cout << "  Point at (u=3, w=2): "; pt.Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "  Gaussian curvature K = " << plane.GaussianCurvature(0.0, 0.0) << " (should be 0)" << std::endl;
}

void Docs_Demo_Surfaces_Tangents_Normals()
{
	std::cout << "\n***********************************************************************" << std::endl;
	std::cout << "****          Surface Tangent Vectors and Normals                 ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	Sphere sphere(2.0);
	Real u = Constants::PI / 4;
	Real w = Constants::PI / 3;
	
	std::cout << "Sphere R=2 at (u=PI/4, w=PI/3):" << std::endl;
	
	Vec3Cart pt = sphere(u, w);
	std::cout << "  Point: "; pt.Print(std::cout, 8, 4); std::cout << std::endl;
	
	Vec3Cart tU, tW;
	sphere.Tangents(u, w, tU, tW);
	std::cout << "  Tangent in u direction: "; tU.Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "  Tangent in w direction: "; tW.Print(std::cout, 8, 4); std::cout << std::endl;
	
	Vec3Cart n = sphere.Normal(u, w);
	std::cout << "  Unit normal: "; n.Print(std::cout, 8, 4); std::cout << std::endl;
	
	// Verify orthogonality
	std::cout << "\n  Orthogonality check:" << std::endl;
	std::cout << "    tU · n = " << ScalarProduct(tU, n) << " (should be ~0)" << std::endl;
	std::cout << "    tW · n = " << ScalarProduct(tW, n) << " (should be ~0)" << std::endl;
	
	// For sphere, normal should be radial
	std::cout << "\n  Sphere normal should be radial (point/R):" << std::endl;
	Vec3Cart expected_n = pt / pt.NormL2();
	std::cout << "    Expected: "; expected_n.Print(std::cout, 8, 4); std::cout << std::endl;
	std::cout << "    Computed: "; n.Print(std::cout, 8, 4); std::cout << std::endl;
}

void Docs_Demo_Surfaces_Classification()
{
	std::cout << "\n***********************************************************************" << std::endl;
	std::cout << "****          Surface Point Classification                        ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	std::cout << "\nPoint classification based on Gaussian curvature K:" << std::endl;
	std::cout << "  K > 0: Elliptic (dome-like)" << std::endl;
	std::cout << "  K < 0: Hyperbolic (saddle-like)" << std::endl;
	std::cout << "  K = 0: Parabolic (cylinder-like) or flat" << std::endl;

	// Sphere - elliptic everywhere
	std::cout << "\n--- Sphere (elliptic everywhere) ---" << std::endl;
	Sphere sphere(3.0);
	Real K = sphere.GaussianCurvature(Constants::PI/2, 0.0);
	std::cout << "  K = " << K << " > 0: elliptic (dome)" << std::endl;

	// Cylinder - parabolic everywhere
	std::cout << "\n--- Cylinder (parabolic everywhere) ---" << std::endl;
	Cylinder cylinder(2.0, 3.0);
	K = cylinder.GaussianCurvature(Constants::PI/4, 1.0);
	std::cout << "  K = " << K << " = 0: parabolic" << std::endl;
	std::cout << "  isFlat: " << (cylinder.isFlat(Constants::PI/4, 1.0) ? "yes" : "no") << std::endl;

	// Hyperboloid - hyperbolic everywhere
	std::cout << "\n--- Hyperboloid (hyperbolic everywhere) ---" << std::endl;
	Hyperboloid hyperboloid(1.0, 1.0, 1.0);
	K = hyperboloid.GaussianCurvature(0.5, 0.0);
	std::cout << "  K = " << K << " < 0: hyperbolic (saddle)" << std::endl;

	// Torus - varies
	std::cout << "\n--- Torus (varies by location) ---" << std::endl;
	Torus torus(3.0, 1.0);
	
	std::cout << "  Outer rim (u=0, w=0):" << std::endl;
	K = torus.GaussianCurvature(0.0, 0.0);
	std::cout << "    K = " << K << " > 0: elliptic" << std::endl;
	
	std::cout << "  Inner rim (u=0, w=PI):" << std::endl;
	K = torus.GaussianCurvature(0.0, Constants::PI);
	std::cout << "    K = " << K << " < 0: hyperbolic" << std::endl;
	
	std::cout << "  Top/bottom (u=PI/2, w=PI/2):" << std::endl;
	K = torus.GaussianCurvature(Constants::PI/2, Constants::PI/2);
	std::cout << "    K = " << K << " (parabolic circles)" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
/////                              MAIN DEMO FUNCTIONS                                /////
///////////////////////////////////////////////////////////////////////////////////////////

void Docs_Demo_Diff_geometry_curves()
{
	Docs_Demo_Curves_2D_Basic();
	Docs_Demo_Curves_3D_Basic();
	Docs_Demo_Curves_Frenet_Frame();
	Docs_Demo_Curves_Curvature_Torsion();
	Docs_Demo_Curves_Arc_Length();
	Docs_Demo_Curves_Planes();
}

void Docs_Demo_Diff_geometry_surfaces()
{
	Docs_Demo_Surfaces_Basic();
	Docs_Demo_Surfaces_Curvature();
	Docs_Demo_Surfaces_Principal_Curvatures();
	Docs_Demo_Surfaces_Fundamental_Forms();
	Docs_Demo_Surfaces_Special();
	Docs_Demo_Surfaces_Tangents_Normals();
	Docs_Demo_Surfaces_Classification();
}

void Docs_Demo_Diff_geometry()
{
	Docs_Demo_Diff_geometry_curves();
	Docs_Demo_Diff_geometry_surfaces();
}
