#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/VectorTypes.h"

#include "core/Derivation.h"
#endif


using namespace MML;

// TODO 0.9
class Body2D
{
    double _mass;
    Vector2Cartesian _position;
    Vector2Cartesian _velocity;
};

class Body3D
{
    double _mass;
    Vector3Cartesian _position;
    Vector3Cartesian _velocity;

};

// 2D i 3D verzija
class CollisionSimulator2D
{
public:
    CollisionSimulator2D(const std::vector<Body2D>& bodies)
        : _bodies(bodies)
    {
    }

    std::vector<Body2D> _bodies;

    // HasBodiesCollided
    // SimulateOneStep(dT)
};

class CollisionSimulator3D
{
public:
    CollisionSimulator3D(const std::vector<Body3D>& bodies)
        : _bodies(bodies)
    {
    }

    std::vector<Body3D> _bodies;

    // HasBodiesCollided
    // SimulateOneStep(dT)
};

void Example2_collision_calculator()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                 EXAMPLE 2 - collision calculator              ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // napraviti 2D prikaz u konzoli s rubovima i kuglama kao * koje se micu!!!
}