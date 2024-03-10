#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Geometry3D.h"

#include "core/Function.h"
#include "core/Derivation.h"
#include "core/Serializer.h"
#include "core/Visualizer.h"
#endif


using namespace MML;

struct LineCurrent {
	double _currentI;
	Line3D _line;
};
class InfiniteLineCurrent_Field_B : public IVectorFunction<3>
{
	std::vector<LineCurrent> _lines;
public:
	void AddLine(double currentI, const Line3D& line) { _lines.push_back({ currentI, line }); }

	VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const override {
		VectorN<Real, 3>  ret;
		Point3Cartesian   pos(x[0], x[1], x[2]);

		for (int i = 0; i < _lines.size(); i++) {
			Point3Cartesian  nearest_point = _lines[i]._line.NearestPoint(pos);
			Vector3Cartesian vec_to_pos(nearest_point, pos);
			Vector3Cartesian fieldDirection = VectorProd(_lines[i]._line.Direction(), vec_to_pos);

			double B_magnitude = _lines[i]._currentI / (2 * Constants::PI * pos.Dist(nearest_point));

			ret = ret + B_magnitude * fieldDirection.GetAsUnitVector();
		}
		return ret;
	}
};

void Example_Main()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                         MAIN EXAMPLE                          ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	InfiniteLineCurrent_Field_B   EM_field;
	EM_field.AddLine(300.0, Line3D(Point3Cartesian(120, 50, -50), Vector3Cartesian(0, 1, 1)));
	EM_field.AddLine(200.0, Line3D(Point3Cartesian(-150, 100, 0), Vector3Cartesian(0, 0, 1)));
	EM_field.AddLine(200.0, Line3D(Point3Cartesian(20, -100, 00), Vector3Cartesian(1, 0, 0)));

	Visualizer::VisualizeVectorField3DCartesian(EM_field, "EM_field", -300.0, 300.0, 30, -300.0, 300.0, 30, -300.0, 300.0, 30, "..\\..\\results\\EM_field.txt");
}