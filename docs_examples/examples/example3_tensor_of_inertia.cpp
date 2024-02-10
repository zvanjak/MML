#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/VectorTypes.h"

#include "core/Derivation.h"
#include "core/CoordTransf.h"
#endif


using namespace MML;

// TODO 0.8

struct DiscreteMass
{
    Vector3Cartesian _position;
    double _mass;

    DiscreteMass(const Vector3Cartesian &position, const double& mass)
        : _position(position), _mass(mass)
    { }
};

struct DiscreteMassesConfig
{
    std::vector<DiscreteMass> _masses;

    DiscreteMassesConfig(const std::vector<DiscreteMass>& masses)
        : _masses(masses)
    { }
};

class DiscreteMassMomentOfInertiaTensorCalculator
{
    DiscreteMassesConfig _massesConfig;
    // treba zadati os kalkulacije
public:
    DiscreteMassMomentOfInertiaTensorCalculator(const DiscreteMassesConfig& massesConfig)
        : _massesConfig(massesConfig)
    { }

    Tensor2<3> calculate()
    {
        Tensor2<3> tensor(2,0);      // can be (0,2) or (1,1) as well (it is Cartesian tensor)
        for (const auto& mass : _massesConfig._masses)
        {
            Vector3Cartesian pos = mass._position;
            tensor.Component(0,0) += mass._mass * (pos.Y() * pos.Y() + pos.Z() * pos.Z());
            tensor.Component(1,1) += mass._mass * (pos.X() * pos.X() + pos.Z() * pos.Z());
            tensor.Component(2,2) += mass._mass * (pos.X() * pos.X() + pos.Y() * pos.Y());
            
            tensor.Component(0,1) -= mass._mass * pos.X() * pos.Y();
            tensor.Component(0,2) -= mass._mass * pos.X() * pos.Z();
            tensor.Component(1,2) -= mass._mass * pos.Y() * pos.Z();

            tensor.Component(1,0) = tensor.Component(0,1);
            tensor.Component(2,0) = tensor.Component(0,2);
            tensor.Component(2,1) = tensor.Component(1,2);
        }
        return tensor;
    }
};

// continuous mass
// ima kljucno definiranu funkciju gustoce na prostoru
// treba zadati os kalkulacije
class ContinuousMassMomentOfInertiaTensorCalculator
{
    // treba zadati os kalkulacije
public:
    ContinuousMassMomentOfInertiaTensorCalculator()
    {
    }

    Tensor2<3> calculate()
    {
        Tensor2<3> tensor(2,0);
        return tensor;
    }
};

void Example3_tensor_of_inertia()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                  EXAMPLE 3 - tensor of inertia                ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // define set of discrete masses
    double a = 1;
    Vector3Cartesian pos1(a, a, 0);
    Vector3Cartesian pos2(-a, a, 0);
    Vector3Cartesian pos3(-a, -a, 0);
    Vector3Cartesian pos4(a, -a, 0);
    Vector3Cartesian pos5(0, 0, 4*a);
    double m1 = 2;
    double m2 = 1;
    double m3 = 4;
    double m4 = 1;
    double m5 = 1;

    DiscreteMass mass1(pos1, m1);
    DiscreteMass mass2(pos2, m2);
    DiscreteMass mass3(pos3, m3);
    DiscreteMass mass4(pos4, m4);
    DiscreteMass mass5(pos5, m5);

    std::vector<DiscreteMass> masses = {mass1, mass2, mass3, mass4, mass5};
    DiscreteMassesConfig massesConfig(masses);

    DiscreteMassMomentOfInertiaTensorCalculator calculator(massesConfig);
    Tensor2<3> tensor = calculator.calculate();

    std::cout << "Tensor of inertia: " << std::endl;
    std::cout << tensor << std::endl;

    CoordTransfCart3DRotationXAxis transf(30.0 * Constants::PI / 180.0);

    
    Tensor2<3> tensor3 = transf.transfTensor2(tensor, Vector3Cartesian(1, 1, 1));

    std::cout << "Tensor of inertia transf: " << std::endl;
    std::cout << tensor3 << std::endl;



    DiscreteMassesConfig massesTransConfig(masses);
    for (auto& mass : massesTransConfig._masses)
	{
		mass._position = transf.transf(mass._position);
	}

    DiscreteMassMomentOfInertiaTensorCalculator calculator2(massesTransConfig);
    Tensor2<3> tensor2 = calculator2.calculate();

    std::cout << "Tensor of inertia rotated maasses: " << std::endl;
    std::cout << tensor2 << std::endl;
}