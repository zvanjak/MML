#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/core/Surfaces.h"
#endif

#include <iostream>

using namespace MML;
using namespace MML::Surfaces;

void Chapter20_Diff_geometry_surfaces()
{
	// testing surfaces
	PlaneSurface	plane({ 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 });
	Cylinder			cylinder(2.0, 10.0);
	Torus					torus(10.0, 2.0);
	Sphere				sphere(5.0);
	MonkeySaddle	monkeySaddle;
	MobiusStrip		mobiusStrip;
	Ellipsoid			ellipsoid(3.0, 2.0, 1.0);
	Hyperboloid		hyperboloid(2.0, 1.0, 1.0);
	Paraboloid		paraboloid(1.0, 1.0);

	ISurfaceCartesian& surface = plane; // or monkeySaddle, depending on what you want to test

	Real u = 1.0;
	Real w = 2.0;

	std::cout << "Surface at (" << u << ", " << w << "): " << surface(u, w) << std::endl;

	// is flat?
	std::cout << "Is surface flat at (" << u << ", " << w << "): " << (surface.isFlat(u, w) ? "Yes" : "No") << std::endl;

	// calculate first normal form coefficients at (u, w) on the given surface
	Real E, F, G;
	surface.GetFirstNormalFormCoefficients(u, w, E, F, G);
	std::cout << "First normal form coefficients at (" << u << ", " << w << "): E=" << E << ", F=" << F << ", G=" << G << std::endl;

	// calculate second fundamental form coefficients at (u, w) on the given surface
	Real L, M, N;
	surface.GetSecondNormalFormCoefficients(u, w, L, M, N);
	std::cout << "Second normal form coefficients at (" << u << ", " << w << "): L=" << L << ", M=" << M << ", N=" << N << std::endl;

	std::cout << "Gaussian curvature at (" << u << ", " << w << "): " << surface.GaussianCurvature(u, w) << std::endl;
	std::cout << "Mean curvature at (" << u << ", " << w << "): " << surface.MeanCurvature(u, w) << std::endl;
}
