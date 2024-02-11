#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "interfaces/IODESystem.h"

#include "base/Geometry3D.h"

#include "core/Derivation.h"
#include "core/Fields.h"
#endif


using namespace MML;

// TODO 0.8
class GravityMass
{
public:
    Real _mass;
    Vector3Cartesian _position;
    Vector3Cartesian _velocity;

    GravityMass(const Real& mass, const Vector3Cartesian &position) : _mass(mass), _position(position) {}
    GravityMass(const Real& mass, const Vector3Cartesian &position, const Vector3Cartesian &velocity) : _mass(mass), _position(position), _velocity(velocity) {}
};

class GravityMultibodyConfigCart
{
    std::vector<GravityMass> _masses;
public:
    int NumBodies() const { return (int) _masses.size(); }
    
    void AddBody(Real mass, Vector3Cartesian position, Vector3Cartesian velocity)
    {
        _masses.push_back(GravityMass(mass, position, velocity));
    }

    Real Mass(int i) const { return _masses[i]._mass; }
    Vector3Cartesian Position(int i) const { return _masses[i]._position; }
    Vector3Cartesian Velocity(int i) const { return _masses[i]._velocity; }
};

class MultibodyGravityPotentialFieldCart : public IScalarFunction<3>
{
private:
    Real _G;
    GravityMultibodyConfigCart _config;
public:
    MultibodyGravityPotentialFieldCart(GravityMultibodyConfigCart config) : _G(6.67430e-11), _config(config) {}
    MultibodyGravityPotentialFieldCart(Real G, GravityMultibodyConfigCart config) : _G(G), _config(config) {}

    Real operator()(const VectorN<Real, 3> &x) const  
    { 
        Real pot = 0.0;
        for(int i = 0; i < _config.NumBodies(); i++)
        {
            pot += InverseRadialPotentialFieldCart(_G * _config.Mass(i), _config.Position(i) - x);
        }
        return pot;
    }
};     
class MultibodyGravityForceFieldCart : public IVectorFunction<3>
{
private:
    Real _G;
    GravityMultibodyConfigCart _config;
public:
    MultibodyGravityForceFieldCart(GravityMultibodyConfigCart config) : _G(6.67430e-11), _config(config) {}
    MultibodyGravityForceFieldCart(Real G, GravityMultibodyConfigCart config) : _G(G), _config(config) {}

    VectorN<Real, 3> operator()(const VectorN<Real, 3> &x) const  
    { 
        VectorN<Real, 3> force(0.0);
        for(int i = 0; i < _config.NumBodies(); i++)
        {
            force = force + InverseRadialPotentialForceFieldCart(_G * _config.Mass(i), _config.Position(i) - x);
        }
        return force;
    }
};

// set up ODESystem za simuliranje sustava tijela u gravitacijskom polju
// a ima i druga opcija - u DANOM sustavu tijela, simulirati neko trece tijelo!
class GravitySystemMotionSolverSimple 
{
protected:
    GravityMultibodyConfigCart _config;

public:
    GravitySystemMotionSolverSimple(GravityMultibodyConfigCart inConfig) : _config(inConfig) { }
    
    int getDim() const { return 3 * _config.NumBodies(); }
    // simulate one step by calculating force on each body and updating velocity and position

    void simulateOneStep(const Real dt) const 
    {
        for(int i = 0; i < _config.NumBodies(); i++)
        {
            // calculate force on body i
            Vector3Cartesian force(0.0);
            for(int j = 0; j < _config.NumBodies(); j++)
            {
                if(i != j)
                {
                    force = force + InverseRadialPotentialForceFieldCart(6.67430e-11 * _config.Mass(j), _config.Position(j) - _config.Position(i));
                }
            }
            // update velocity and position
            // eeee, nemozeee, poremetit ce s eostatak petlje!!!
            // treba radit s kopijom!
            _config.Velocity(i) = _config.Velocity(i) + force * dt / _config.Mass(i);
            _config.Position(i) = _config.Position(i) + _config.Velocity(i) * dt;
        }
    }
};

class GravitySystemMotionSolver : public IODESystem
{
protected:
    GravityMultibodyConfigCart _initialConfig;

public:
    GravitySystemMotionSolver(GravityMultibodyConfigCart inConfig) : _initialConfig(inConfig) { }
    
    int getDim() const { return 3 * _initialConfig.NumBodies(); }
    void derivs(const Real t, const Vector<Real> &x, Vector<Real> &dxdt) const 
    {
        // moramo popuniti dxdt vektor s N * 6 derivacija nasih varijabli (x, y, z, vx, vy, vz)
        for(int i = 0; i < _initialConfig.NumBodies(); i++)
        {
            Vector3Cartesian force(0,0,0);
            // calculating force on body i
            for(int j = 0; j < _initialConfig.NumBodies(); j++)
            {
                if(i != j)      // self-force would be infinite :)
                {
                    Vector3Cartesian vec_dist(dxdt[6 * j] - dxdt[6 * i], dxdt[6 * j + 2] - dxdt[6 * i + 2], dxdt[6 * j + 4] - dxdt[6 * i + 4]);
                    Vector3Cartesian f = InverseRadialPotentialForceFieldCart(6.67430e-11 * _initialConfig.Mass(i) * _initialConfig.Mass(j), vec_dist);
                    force = force + f;
                }
            }

            // x coord
            dxdt[6 * i]     = x[6 * i + 1];
            dxdt[6 * i + 1] = 1 / _initialConfig.Mass(i) * force.X();
            // y coord
            dxdt[6 * i + 2] = x[6 * i + 3];
            dxdt[6 * i + 3] = 1 / _initialConfig.Mass(i) * force.Y();
            // z coord
            dxdt[6 * i + 4] = x[6 * i + 5];
            dxdt[6 * i + 5] = 1 / _initialConfig.Mass(i) * force.Z();
        }
    }
};

void Example4_Gravity_field_visualization()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****            EXAMPLE 4 - gravity field investigations           ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;

    // prvo, simulacija 5 tijela

    // drugo, Voyager kroz solarni sustav
}