///////////////////////////////////////////////////////////////////////////////////////////
// Chapter 02 - Parametric Surface Visualization Demos
// 
// This file demonstrates MML's parametric surface visualization capabilities.
// Surfaces are 2-parameter mappings (u,v) -> (x,y,z) in 3D space.
//
// Demos:
//   1. Torus - Classic donut shape with major/minor radii
//   2. Möbius Strip - Non-orientable surface with one side
//   3. Enneper Surface - Elegant minimal surface from differential geometry
//   4. Klein Bottle - Non-orientable closed surface (figure-8 immersion)
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Vector/VectorN.h"
#include "mml/base/Function.h"
#include "mml/core/Surfaces.h"
#include "mml/tools/Visualizer.h"
#endif

using namespace MML;
using namespace MML::Surfaces;

///////////////////////////////////////////////////////////////////////////////////////////
//                           PARAMETRIC SURFACE DEMOS
///////////////////////////////////////////////////////////////////////////////////////////

void Chapter02_Demo_ParametricSurfaces()
{
	std::cout << "\n=== Chapter 02: Parametric Surface Visualization ===" << std::endl;

	// Demo 1: Torus - the classic donut
	// R = major radius (center of tube to center of torus)
	// r = minor radius (radius of the tube)
	{
		std::cout << "  Visualizing Torus (R=3, r=1)..." << std::endl;
		Surfaces::Torus torus(3.0, 1.0);
		auto result = Visualizer::VisualizeParametricSurface(
			torus, 
			"Torus (R=3, r=1)", 
			60, 30,  // Resolution: 60 points around major circle, 30 around minor
			"ch02_surface_torus.mml");
		if (!result.success) {
			std::cerr << "    Failed: " << result.errorMessage << std::endl;
		}
	}

	// Demo 2: Möbius Strip - non-orientable surface with one side and one edge
	// Famous example: walk along it and you end up on the "other side"
	{
		std::cout << "  Visualizing Mobius Strip..." << std::endl;
		Surfaces::MobiusStrip mobius(2.0);  // Scale factor for visibility
		auto result = Visualizer::VisualizeParametricSurface(
			mobius, 
			"Mobius Strip", 
			80, 20,  // More points around the strip, fewer across width
			"ch02_surface_mobius.mml");
		if (!result.success) {
			std::cerr << "    Failed: " << result.errorMessage << std::endl;
		}
	}

	// Demo 3: Enneper Surface - minimal surface with elegant saddle shape
	// Self-intersects but has zero mean curvature everywhere
	{
		std::cout << "  Visualizing Enneper Surface..." << std::endl;
		Surfaces::EnneperSurface enneper(1.5);  // Scale for visibility
		auto result = Visualizer::VisualizeParametricSurface(
			enneper, 
			"Enneper Minimal Surface", 
			50, 50,  // Equal resolution in both directions
			"ch02_surface_enneper.mml");
		if (!result.success) {
			std::cerr << "    Failed: " << result.errorMessage << std::endl;
		}
	}

	// Demo 4: Klein Bottle - non-orientable closed surface
	// Like a Möbius strip but closed - has no boundary
	// Uses figure-8 immersion to display in 3D (self-intersects)
	{
		std::cout << "  Visualizing Klein Bottle..." << std::endl;
		Surfaces::KleinBottle klein(0.5);  // Scale factor
		auto result = Visualizer::VisualizeParametricSurface(
			klein, 
			"Klein Bottle (Figure-8 Immersion)", 
			80, 40,  // High resolution for the complex shape
			"ch02_surface_klein.mml");
		if (!result.success) {
			std::cerr << "    Failed: " << result.errorMessage << std::endl;
		}
	}

	std::cout << "  Parametric surface demos complete!" << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////
//                           MAIN FUNCTION FOR THIS MODULE
///////////////////////////////////////////////////////////////////////////////////////////

void Chapter02_ParametricSurfaceVisualization()
{
	std::cout << "\n========================================" << std::endl;
	std::cout << "Chapter 02: Parametric Surface Visualization" << std::endl;
	std::cout << "========================================" << std::endl;

	Chapter02_Demo_ParametricSurfaces();

	std::cout << "\nAll parametric surface visualizations complete!" << std::endl;
}
