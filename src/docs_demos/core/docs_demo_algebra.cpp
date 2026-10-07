#include <mml/base/Algebra_base.h>

#include <iostream>
#include <vector>

using namespace MML;
using namespace MML::Algebra;

void Docs_Demo_Algebra()
{
	CyclicGroup cyclic(5);
	DihedralGroup dihedral(5);
	auto vertexAction = MakeGroupAction<DihedralElement, int>(
		[&dihedral](const DihedralElement& element, const int& vertex) {
			return dihedral.permutation(element).apply(vertex);
		});
	auto orbit = Orbit(dihedral, vertexAction, 0);
	auto stabilizer = Stabilizer(dihedral, vertexAction, 0);

	std::vector<std::vector<int>> colorings;
	for (int mask = 0; mask < 16; mask++) {
		std::vector<int> coloring(4);
		for (int vertex = 0; vertex < 4; vertex++) coloring[vertex] = (mask >> vertex) & 1;
		colorings.push_back(coloring);
	}
	CyclicGroup rotations(4);
	auto coloringAction = MakeGroupAction<CyclicElement, std::vector<int>>(
		[&rotations](const CyclicElement& element, const std::vector<int>& coloring) {
			return rotations.permutation(element).apply(coloring);
		});

	auto spatial = SO3::FromAxisAngle({0, 0, 1}, Constants::PI / REAL(2.0));
	auto rotated = spatial.apply({1, 0, 0});

	std::cout << "C5 order: " << cyclic.order() << ", D5 order: " << dihedral.order() << '\n';
	std::cout << "vertex orbit/stabilizer: " << orbit.size() << "/" << stabilizer.size() << '\n';
	std::cout << "binary necklaces of length 4: " << BurnsideCount(rotations, coloringAction, colorings) << '\n';
	std::cout << "SO3 quarter-turn of e1: (" << rotated[0] << ", " << rotated[1] << ", " << rotated[2] << ")\n";
}