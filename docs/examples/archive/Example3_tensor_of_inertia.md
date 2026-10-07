# Calculating tensor of inertia

## Calculation for a set of discrete masses

Classes for modeling set of discrete masses and calculating tensor of inertia.
~~~c++
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
~~~

Calculating tensor of inertia:

~~~c++
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
Tensor2<3> tensor_orig = calculator.calculate();

std::cout << "Tensor of inertia: " << std::endl;
std::cout << tensor_orig << std::endl;

// investigating what happens if we change coord.system, int two cases:
// 1. using coord.system transform we calculate TRANSFORMED (original) tensor
// 2. using coord.system transform we calculate NEW set of masses and then calculate tensor

// new coord.system is rotated around x axis for 30 degrees
CoordTransfCart3DRotationXAxis coord_transf(30.0 * Constants::PI / 180.0);

// 1) - calculated tensor transformation
Tensor2<3> tensor_transf = coord_transf.transfTensor2(tensor_orig, Vector3Cartesian(1, 1, 1));

std::cout << "Tensor of inertia transformed: " << std::endl;
std::cout << tensor_transf << std::endl;

// 2) - change masses position and calculate new tensor 
DiscreteMassesConfig massesTransConfig(masses);
for (auto& mass : massesTransConfig._masses)
    mass._position = coord_transf.transf(mass._position);

DiscreteMassMomentOfInertiaTensorCalculator calculator2(massesTransConfig);
Tensor2<3> tensor_changed = calculator2.calculate();

std::cout << "Tensor of inertia rotated maasses: " << std::endl;
std::cout << tensor_changed << std::endl;
~~~

## Calculation for a continuous mass distribution (solid body)


[C++ file with implementation](/docs_examples/examples/example3_tensor_of_inertia.cpp)

