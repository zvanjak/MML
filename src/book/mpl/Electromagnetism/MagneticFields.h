#if !defined MPL_MAGNETICFIELDS_H
#define MPL_MAGNETICFIELDS_H

#include "MMLBase.h"

#include "base/Geometry3D.h"

#include "base/Function.h"

using namespace MML;

namespace MPL
{
	struct LineCurrent {
		double _currentI;
		Line3D _line;
	};

	class SingleInfiniteLineCurrent_Field_B : public IVectorFunction<3>
	{
		LineCurrent _lineCurrent;
	
	public:
		SingleInfiniteLineCurrent_Field_B(double currentI, const Line3D& line) : _lineCurrent{ currentI, line } {}

		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const override 
		{
			VectorN<Real, 3>  ret;
			Point3Cartesian   pos(x[0], x[1], x[2]);
			
			// find nearest point on line to given point
			Point3Cartesian		nearest_point = _lineCurrent._line.NearestPointOnLine(pos);
			Vector3Cartesian	vec_to_pos(nearest_point, pos);
			
			Vector3Cartesian	fieldDirection = VectorProduct(_lineCurrent._line.Direction(), vec_to_pos);
			
			double B_magnitude = _lineCurrent._currentI / (2 * Constants::PI * pos.Dist(nearest_point));
			ret = B_magnitude * fieldDirection.GetAsUnitVector();
			
			return ret;
		}
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
				Point3Cartesian  nearest_point = _lines[i]._line.NearestPointOnLine(pos);
				Vector3Cartesian vec_to_pos(nearest_point, pos);
				Vector3Cartesian fieldDirection = VectorProduct(_lines[i]._line.Direction(), vec_to_pos);

				double B_magnitude = _lines[i]._currentI / (2 * Constants::PI * pos.Dist(nearest_point));

				ret = ret + B_magnitude * fieldDirection.GetAsUnitVector();
			}
			return ret;
		}
	};

}

#endif // MPL_MAGNETICFIELDS_H