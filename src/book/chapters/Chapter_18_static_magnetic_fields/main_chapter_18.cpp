#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Geometry3D.h"

#include "mml/base/Function.h"
#include "mml/core/Derivation.h"
#include "mml/tools/Serializer.h"
#include "mml/tools/Visualizer.h"

#include "mpl/Electromagnetism/MagneticFields.h"
#endif


using namespace MML;
using namespace MPL;

void Chapter18_Infinite_line_magnetic_field()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****          Example 18 - infinite line magnetic field            ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	InfiniteLineCurrent_Field_B   EM_field;
	EM_field.AddLine(300.0, Line3D(Point3Cartesian(120, 50, -50), Vector3Cartesian(0, 1, 1)));
	EM_field.AddLine(200.0, Line3D(Point3Cartesian(-150, 100, 0), Vector3Cartesian(0, 0, 1)));
	EM_field.AddLine(200.0, Line3D(Point3Cartesian(20, -100, 00), Vector3Cartesian(1, 0, 0)));

	// Visualizer::VisualizeVectorField3DCartesian(EM_field, "EM_field", -300.0, 300.0, 30, -300.0, 300.0, 30, -300.0, 300.0, 30, "EM_field.txt");
}

