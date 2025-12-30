#if !defined MPL_MOMENT_OF_INERTIA_CALCULATOR_H
#define MPL_MOMENT_OF_INERTIA_CALCULATOR_H


#include "MMLBase.h"

#include "base/VectorN.h"
#include "base/Geometry3DBodies.h"

#include "core/Derivation.h"
#include "core/Integration.h"
#include "core/CoordTransf.h"
#include "core/CoordTransf/CoordTransf3D.h"

using namespace MML;

namespace MPL
{
	struct DiscreteMass {
		Vector3Cartesian _position;
		double _mass;

		DiscreteMass(const Vector3Cartesian& position, const double& mass)
			: _position(position), _mass(mass) {
		}
	};
	struct DiscreteMassesConfig {
		std::vector<DiscreteMass> _masses;

		DiscreteMassesConfig(const std::vector<DiscreteMass>& masses)
			: _masses(masses) {
		}
	};

	class DiscreteMassMomentOfInertiaTensorCalculator
	{
		DiscreteMassesConfig _massesConfig;
	public:
		DiscreteMassMomentOfInertiaTensorCalculator(const DiscreteMassesConfig& massesConfig)
			: _massesConfig(massesConfig) {
		}

		Tensor2<3> calculate()
		{
			Tensor2<3> tensor(2, 0);  // can be (0,2) or (1,1) as well (it is a Cartesian tensor)
			for (const auto& mass : _massesConfig._masses)
			{
				Vector3Cartesian pos = mass._position;
				tensor(0, 0) += mass._mass * (pos.Y() * pos.Y() + pos.Z() * pos.Z());
				tensor(1, 1) += mass._mass * (pos.X() * pos.X() + pos.Z() * pos.Z());
				tensor(2, 2) += mass._mass * (pos.X() * pos.X() + pos.Y() * pos.Y());

				tensor(0, 1) -= mass._mass * pos.X() * pos.Y();
				tensor(0, 2) -= mass._mass * pos.X() * pos.Z();
				tensor(1, 2) -= mass._mass * pos.Y() * pos.Z();
			}
			tensor(1, 0) = tensor(0, 1);
			tensor(2, 0) = tensor(0, 2);
			tensor(2, 1) = tensor(1, 2);

			return tensor;
		}
	};

	// continuous mass with defined boundaries (needs density function)

	class ContMassScalarFuncBase : public IScalarFunction<3>
	{
	protected:
		ISolidBodyWithBoundary& _mass;
	public:
		ContMassScalarFuncBase(ISolidBodyWithBoundary& mass) : _mass(mass) {}
	};

	struct Func11 : public ContMassScalarFuncBase
	{
		Func11(ISolidBodyWithBoundary& mass) : ContMassScalarFuncBase(mass) {}
		Real operator()(const VectorN<Real, 3>& x) const
		{
			return _mass.getDensity(x) * (x[1] * x[1] + x[2] * x[2]);
		}
	};
	struct Func22 : public ContMassScalarFuncBase
	{
		Func22(ISolidBodyWithBoundary& mass) : ContMassScalarFuncBase(mass) {}
		Real operator()(const VectorN<Real, 3>& x) const
		{
			return _mass.getDensity(x) * (x[0] * x[0] + x[2] * x[2]);
		}
	};
	struct Func33 : public ContMassScalarFuncBase
	{
		Func33(ISolidBodyWithBoundary& mass) : ContMassScalarFuncBase(mass) {}
		Real operator()(const VectorN<Real, 3>& x) const
		{
			return _mass.getDensity(x) * (x[0] * x[0] + x[1] * x[1]);
		}
	};
	struct Func12 : public ContMassScalarFuncBase
	{
		Func12(ISolidBodyWithBoundary& mass) : ContMassScalarFuncBase(mass) {}
		Real operator()(const VectorN<Real, 3>& x) const
		{
			return -_mass.getDensity(x) * x[0] * x[1];
		}
	};
	struct Func13 : public ContMassScalarFuncBase
	{
		Func13(ISolidBodyWithBoundary& mass) : ContMassScalarFuncBase(mass) {}
		Real operator()(const VectorN<Real, 3>& x) const
		{
			return -_mass.getDensity(x) * x[0] * x[2];
		}
	};
	struct Func23 : public ContMassScalarFuncBase
	{
		Func23(ISolidBodyWithBoundary& mass) : ContMassScalarFuncBase(mass) {}
		Real operator()(const VectorN<Real, 3>& x) const
		{
			return -_mass.getDensity(x) * x[1] * x[2];
		}
	};

	struct HelperIntegrator : public IScalarFunction<3>
	{
	protected:
		ISolidBodyWithBoundary& _mass;
		Real(*_func)(const VectorN<Real, 3>& x);

	public:
		HelperIntegrator(ISolidBodyWithBoundary& mass, Real(*func)(const VectorN<Real, 3>& x)) : _mass(mass), _func(func) {}
		
		Real operator()(const VectorN<Real, 3>& x) const
		{
			return _mass.getDensity(x) * _func(x);
		}
	};

	class ContinuousMassMomentOfInertiaTensorCalculator
	{
		ISolidBodyWithBoundary& _mass;
	public:
		ContinuousMassMomentOfInertiaTensorCalculator(ISolidBodyWithBoundary& mass)
			: _mass(mass)
		{	}

		Tensor2<3> calculate()
		{
			Tensor2<3> tensor(2, 0);

			HelperIntegrator _f11(_mass, [](const Vec3 &x) { return x[1] * x[1] + x[2] * x[2]; } );
			//Func11 _f11(_mass);
			Func22 _f22(_mass);
			Func33 _f33(_mass);
			Func12 _f12(_mass);
			Func13 _f13(_mass);
			Func23 _f23(_mass);
			tensor(0, 0) = Integrate3D(_f11, _mass._x1, _mass._x2, _mass._y1, _mass._y2, _mass._z1, _mass._z2);
			tensor(1, 1) = Integrate3D(_f22, _mass._x1, _mass._x2, _mass._y1, _mass._y2, _mass._z1, _mass._z2);
			tensor(2, 2) = Integrate3D(_f33, _mass._x1, _mass._x2, _mass._y1, _mass._y2, _mass._z1, _mass._z2);
			tensor(0, 1) = Integrate3D(_f12, _mass._x1, _mass._x2, _mass._y1, _mass._y2, _mass._z1, _mass._z2);
			tensor(0, 2) = Integrate3D(_f13, _mass._x1, _mass._x2, _mass._y1, _mass._y2, _mass._z1, _mass._z2);
			tensor(1, 2) = Integrate3D(_f23, _mass._x1, _mass._x2, _mass._y1, _mass._y2, _mass._z1, _mass._z2);
			tensor(1, 0) = tensor(0, 1);
			tensor(2, 0) = tensor(0, 2);
			tensor(2, 1) = tensor(1, 2);

			ScalarFunctionFromStdFunc<3> fDensity(std::function<Real(const VectorN<Real, 3>&)>{ std::bind(&ISolidBodyWithBoundary::getDensity, &_mass, std::placeholders::_1) });

			Real vol = Integrate3D(fDensity,
														_mass._x1, _mass._x2,
														_mass._y1, _mass._y2,
														_mass._z1, _mass._z2);

			return tensor;
		}
	};
} // namespace MPL

#endif // MPL_MOMENT_OF_INERTIA_CALCULATOR_H