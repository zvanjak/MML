/**
 * @file generate_parametric_surface_data.cpp
 * @brief Generate sample parametric surface data for visualization testing
 * 
 * Creates serialized data files for various interesting parametric surfaces
 * that can be used to test and develop the ParametricSurface visualizer.
 */

#include "MMLBase.h"
#include "base/Vector/Vector.h"
#include "base/Vector/VectorN.h"
#include "interfaces/IFunction.h"
#include "tools/Serializer.h"

#include <iostream>
#include <cmath>

using namespace MML;

///////////////////////////////////////////////////////////////////////////////
// 1. Sphere - the classic parametric surface
// x = r * cos(u) * sin(w)
// y = r * sin(u) * sin(w)  
// z = r * cos(w)
// u ∈ [0, 2π], w ∈ [0, π]
///////////////////////////////////////////////////////////////////////////////
class Sphere : public IParametricSurface<3>
{
    Real _radius;
public:
    Sphere(Real radius = 1.0) : _radius(radius) {}
    
    VectorN<Real, 3> operator()(Real u, Real w) const override
    {
        return VectorN<Real, 3>{
            _radius * std::cos(u) * std::sin(w),
            _radius * std::sin(u) * std::sin(w),
            _radius * std::cos(w)
        };
    }
    
    Real getMinU() const override { return 0.0; }
    Real getMaxU() const override { return 2.0 * Constants::PI; }
    Real getMinW(Real) const override { return 0.0; }
    Real getMaxW(Real) const override { return Constants::PI; }
};

///////////////////////////////////////////////////////////////////////////////
// 2. Torus - donut shape
// x = (R + r*cos(w)) * cos(u)
// y = (R + r*cos(w)) * sin(u)
// z = r * sin(w)
// u, w ∈ [0, 2π]
///////////////////////////////////////////////////////////////////////////////
class Torus : public IParametricSurface<3>
{
    Real _R;  // Major radius (distance from center to tube center)
    Real _r;  // Minor radius (tube radius)
public:
    Torus(Real majorRadius = 2.0, Real minorRadius = 0.5) 
        : _R(majorRadius), _r(minorRadius) {}
    
    VectorN<Real, 3> operator()(Real u, Real w) const override
    {
        return VectorN<Real, 3>{
            (_R + _r * std::cos(w)) * std::cos(u),
            (_R + _r * std::cos(w)) * std::sin(u),
            _r * std::sin(w)
        };
    }
    
    Real getMinU() const override { return 0.0; }
    Real getMaxU() const override { return 2.0 * Constants::PI; }
    Real getMinW(Real) const override { return 0.0; }
    Real getMaxW(Real) const override { return 2.0 * Constants::PI; }
};

///////////////////////////////////////////////////////////////////////////////
// 3. Möbius Strip - non-orientable surface
// x = (1 + (w/2)*cos(u/2)) * cos(u)
// y = (1 + (w/2)*cos(u/2)) * sin(u)
// z = (w/2) * sin(u/2)
// u ∈ [0, 2π], w ∈ [-1, 1]
///////////////////////////////////////////////////////////////////////////////
class MobiusStrip : public IParametricSurface<3>
{
    Real _width;
public:
    MobiusStrip(Real width = 0.5) : _width(width) {}
    
    VectorN<Real, 3> operator()(Real u, Real w) const override
    {
        Real halfU = u / 2.0;
        return VectorN<Real, 3>{
            (1.0 + _width * w * std::cos(halfU)) * std::cos(u),
            (1.0 + _width * w * std::cos(halfU)) * std::sin(u),
            _width * w * std::sin(halfU)
        };
    }
    
    Real getMinU() const override { return 0.0; }
    Real getMaxU() const override { return 2.0 * Constants::PI; }
    Real getMinW(Real) const override { return -1.0; }
    Real getMaxW(Real) const override { return 1.0; }
};

///////////////////////////////////////////////////////////////////////////////
// 4. Klein Bottle (immersion in R³) - non-orientable, self-intersecting
// A classic topological object!
///////////////////////////////////////////////////////////////////////////////
class KleinBottle : public IParametricSurface<3>
{
public:
    VectorN<Real, 3> operator()(Real u, Real v) const override
    {
        Real x, y, z;
        
        if (u < Constants::PI)
        {
            x = 3.0 * std::cos(u) * (1.0 + std::sin(u)) 
                + 2.0 * (1.0 - std::cos(u) / 2.0) * std::cos(u) * std::cos(v);
            y = 8.0 * std::sin(u) 
                + 2.0 * (1.0 - std::cos(u) / 2.0) * std::sin(u) * std::cos(v);
        }
        else
        {
            x = 3.0 * std::cos(u) * (1.0 + std::sin(u)) 
                + 2.0 * (1.0 - std::cos(u) / 2.0) * std::cos(v + Constants::PI);
            y = 8.0 * std::sin(u);
        }
        z = 2.0 * (1.0 - std::cos(u) / 2.0) * std::sin(v);
        
        return VectorN<Real, 3>{x, y, z};
    }
    
    Real getMinU() const override { return 0.0; }
    Real getMaxU() const override { return 2.0 * Constants::PI; }
    Real getMinW(Real) const override { return 0.0; }
    Real getMaxW(Real) const override { return 2.0 * Constants::PI; }
};

///////////////////////////////////////////////////////////////////////////////
// 5. Helicoid - a ruled surface (minimal surface like a spiral ramp)
// x = u * cos(w)
// y = u * sin(w)
// z = c * w
///////////////////////////////////////////////////////////////////////////////
class Helicoid : public IParametricSurface<3>
{
    Real _pitch;  // How much z changes per revolution
public:
    Helicoid(Real pitch = 0.5) : _pitch(pitch) {}
    
    VectorN<Real, 3> operator()(Real u, Real w) const override
    {
        return VectorN<Real, 3>{
            u * std::cos(w),
            u * std::sin(w),
            _pitch * w
        };
    }
    
    Real getMinU() const override { return -1.0; }
    Real getMaxU() const override { return 1.0; }
    Real getMinW(Real) const override { return 0.0; }
    Real getMaxW(Real) const override { return 4.0 * Constants::PI; }  // Two full turns
};


int main()
{
    std::cout << "Generating parametric surface data for visualization...\n\n";
    
    const std::string outputDir = "src/visualization_examples/visualization_data/";
    
    // Grid resolution - balance between detail and file size
    const int gridU = 50;  // Points along u
    const int gridW = 50;  // Points along w
    
    //=========================================================================
    // 1. Sphere
    //=========================================================================
    {
        Sphere sphere(1.5);
        std::string filename = outputDir + "parametric_sphere.mml";
        
        std::cout << "Generating Sphere... ";
        auto result = Serializer::SaveParametricSurface(
            sphere, "Unit Sphere (r=1.5)",
            0.0, 2.0 * Constants::PI, gridU,   // u range [0, 2π]
            0.0, Constants::PI, gridW,          // w range [0, π]
            filename
        );
        
        if (result.success)
            std::cout << "OK -> " << filename << "\n";
        else
            std::cout << "FAILED: " << result.message << "\n";
    }
    
    //=========================================================================
    // 2. Torus
    //=========================================================================
    {
        Torus torus(2.0, 0.7);  // Major radius 2, minor radius 0.7
        std::string filename = outputDir + "parametric_torus.mml";
        
        std::cout << "Generating Torus... ";
        auto result = Serializer::SaveParametricSurface(
            torus, "Torus (R=2.0, r=0.7)",
            0.0, 2.0 * Constants::PI, gridU,
            0.0, 2.0 * Constants::PI, gridW,
            filename
        );
        
        if (result.success)
            std::cout << "OK -> " << filename << "\n";
        else
            std::cout << "FAILED: " << result.message << "\n";
    }
    
    //=========================================================================
    // 3. Möbius Strip
    //=========================================================================
    {
        MobiusStrip mobius(0.4);
        std::string filename = outputDir + "parametric_mobius.mml";
        
        std::cout << "Generating Mobius Strip... ";
        auto result = Serializer::SaveParametricSurface(
            mobius, "Mobius Strip",
            0.0, 2.0 * Constants::PI, gridU,
            -1.0, 1.0, gridW / 2,  // Fewer points across width
            filename
        );
        
        if (result.success)
            std::cout << "OK -> " << filename << "\n";
        else
            std::cout << "FAILED: " << result.message << "\n";
    }
    
    //=========================================================================
    // 4. Klein Bottle (this will look wild!)
    //=========================================================================
    {
        KleinBottle klein;
        std::string filename = outputDir + "parametric_klein_bottle.mml";
        
        std::cout << "Generating Klein Bottle... ";
        auto result = Serializer::SaveParametricSurface(
            klein, "Klein Bottle (immersion)",
            0.0, 2.0 * Constants::PI, gridU,
            0.0, 2.0 * Constants::PI, gridW,
            filename
        );
        
        if (result.success)
            std::cout << "OK -> " << filename << "\n";
        else
            std::cout << "FAILED: " << result.message << "\n";
    }
    
    //=========================================================================
    // 5. Helicoid
    //=========================================================================
    {
        Helicoid helicoid(0.3);
        std::string filename = outputDir + "parametric_helicoid.mml";
        
        std::cout << "Generating Helicoid... ";
        auto result = Serializer::SaveParametricSurface(
            helicoid, "Helicoid (spiral ramp)",
            -1.5, 1.5, gridW / 2,              // u controls radial extent
            0.0, 4.0 * Constants::PI, gridU,   // w controls spiral angle
            filename
        );
        
        if (result.success)
            std::cout << "OK -> " << filename << "\n";
        else
            std::cout << "FAILED: " << result.message << "\n";
    }
    
    std::cout << "\nDone! Generated " << gridU << " x " << gridW << " = " 
              << gridU * gridW << " points per surface.\n";
    
    return 0;
}
