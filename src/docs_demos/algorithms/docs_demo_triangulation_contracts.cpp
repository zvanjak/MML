#include <mml/algorithms/CompGeometry/Triangulation.h>
#include <iostream>
#include <stdexcept>

void Docs_Demo_TriangulationContracts() {
	using namespace MML;
	using namespace MML::CompGeometry;

	const auto skinny = Triangulation::ComputeDelaunay({{0, 0}, {1, 0}, {Real(0.5), Real(1e-6)}});
	if (skinny.NumTriangles() != 1)
		throw std::runtime_error("A noncollinear three-point set must have one face");

	const Real side = Real(1e-6);
	const auto tiny = Triangulation::ComputeDelaunay({{0, 0}, {side, 0}, {0, side}});
	if (tiny.Interpolate({side / 4, side / 4}, {1, 1, 1}, -77) != 1 ||
		tiny.Interpolate({2 * side, 2 * side}, {1, 1, 1}, -77) != -77)
		throw std::runtime_error("Interpolation must distinguish tiny interior points from outside points");

	const Real offset = Real(8192);
	const Polygon2D clockwise({{offset, offset}, {offset, offset + 1},
		{offset + 1, offset + 1}, {offset + 1, offset}});
	const auto clipped = Triangulation::EarClipTriangulationDetailed(clockwise);
	if (!clipped.IsSuccess() || clipped.used_fan_fallback || clipped.triangles.size() != 2)
		throw std::runtime_error("Translation must not change polygon winding or trigger fallback");

	std::cout << "Planar triangulation topology, winding and interpolation contracts verified.\n";
}
