/**
 * @file generate_scalar3d_data.cpp
 * @brief Generate sample 3D scalar function data for visualization testing
 * 
 * Creates serialized data files for various interesting 3D scalar functions
 * that can be used to test and develop the ScalarFunction3D visualizer.
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
// 1. 3D Gaussian Blob - smooth, symmetric, peaks at origin
///////////////////////////////////////////////////////////////////////////////
class GaussianBlob3D : public IScalarFunction<3>
{
    Real _sigma;
public:
    GaussianBlob3D(Real sigma = 1.0) : _sigma(sigma) {}
    
    Real operator()(const VectorN<Real, 3>& x) const override
    {
        Real r2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
        return std::exp(-r2 / (2.0 * _sigma * _sigma));
    }
};

///////////////////////////////////////////////////////////////////////////////
// 2. Triply Periodic Sinusoidal - interesting wave interference pattern
///////////////////////////////////////////////////////////////////////////////
class SinusoidalWave3D : public IScalarFunction<3>
{
    Real _freq;
public:
    SinusoidalWave3D(Real frequency = 1.0) : _freq(frequency) {}
    
    Real operator()(const VectorN<Real, 3>& x) const override
    {
        return std::sin(_freq * x[0]) * std::sin(_freq * x[1]) * std::sin(_freq * x[2]);
    }
};

///////////////////////////////////////////////////////////////////////////////
// 3. Gyroid - beautiful triply periodic minimal surface (implicit)
// The surface is where f(x,y,z) = 0
///////////////////////////////////////////////////////////////////////////////
class Gyroid3D : public IScalarFunction<3>
{
    Real _scale;
public:
    Gyroid3D(Real scale = 1.0) : _scale(scale) {}
    
    Real operator()(const VectorN<Real, 3>& x) const override
    {
        Real sx = _scale * x[0];
        Real sy = _scale * x[1];
        Real sz = _scale * x[2];
        return std::sin(sx) * std::cos(sy) + 
               std::sin(sy) * std::cos(sz) + 
               std::sin(sz) * std::cos(sx);
    }
};

///////////////////////////////////////////////////////////////////////////////
// Bonus: Sphere Distance Function - useful for isosurface testing
///////////////////////////////////////////////////////////////////////////////
class SphereDistance3D : public IScalarFunction<3>
{
    Real _radius;
public:
    SphereDistance3D(Real radius = 1.0) : _radius(radius) {}
    
    Real operator()(const VectorN<Real, 3>& x) const override
    {
        Real r = std::sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
        return r - _radius;  // Negative inside, zero on surface, positive outside
    }
};


int main()
{
    std::cout << "Generating 3D scalar function data for visualization...\n\n";
    
    const std::string outputDir = "src/visualization_examples/visualization_data/";
    
    // Grid resolution - balance between detail and file size
    const int gridPoints = 30;  // 30^3 = 27,000 points
    
    //=========================================================================
    // 1. Gaussian Blob - centered at origin, extends from -3 to 3
    //=========================================================================
    {
        GaussianBlob3D gaussian(1.0);
        std::string filename = outputDir + "gaussian_blob_3d.mml";
        
        std::cout << "Generating Gaussian Blob 3D... ";
        auto result = Serializer::SaveScalarFunc3DCartesian(
            gaussian, "3D Gaussian Blob (sigma=1.0)",
            -3.0, 3.0, gridPoints,   // x range
            -3.0, 3.0, gridPoints,   // y range
            -3.0, 3.0, gridPoints,   // z range
            filename
        );
        
        if (result.success)
            std::cout << "OK -> " << filename << "\n";
        else
            std::cout << "FAILED: " << result.message << "\n";
    }
    
    //=========================================================================
    // 2. Sinusoidal Wave - triply periodic pattern
    //=========================================================================
    {
        SinusoidalWave3D wave(2.0);  // frequency = 2
        std::string filename = outputDir + "sinusoidal_wave_3d.mml";
        
        std::cout << "Generating Sinusoidal Wave 3D... ";
        auto result = Serializer::SaveScalarFunc3DCartesian(
            wave, "Triply Periodic Sinusoidal Wave (freq=2.0)",
            0.0, Constants::PI * 2, gridPoints,   // x range [0, 2π]
            0.0, Constants::PI * 2, gridPoints,   // y range
            0.0, Constants::PI * 2, gridPoints,   // z range
            filename
        );
        
        if (result.success)
            std::cout << "OK -> " << filename << "\n";
        else
            std::cout << "FAILED: " << result.message << "\n";
    }
    
    //=========================================================================
    // 3. Gyroid - triply periodic minimal surface
    //=========================================================================
    {
        Gyroid3D gyroid(1.0);
        std::string filename = outputDir + "gyroid_3d.mml";
        
        std::cout << "Generating Gyroid 3D... ";
        auto result = Serializer::SaveScalarFunc3DCartesian(
            gyroid, "Gyroid Minimal Surface",
            -Constants::PI * 2, Constants::PI * 2, gridPoints,   // x range
            -Constants::PI * 2, Constants::PI * 2, gridPoints,   // y range
            -Constants::PI * 2, Constants::PI * 2, gridPoints,   // z range
            filename
        );
        
        if (result.success)
            std::cout << "OK -> " << filename << "\n";
        else
            std::cout << "FAILED: " << result.message << "\n";
    }
    
    //=========================================================================
    // Bonus: Sphere distance function (good for isosurface extraction testing)
    //=========================================================================
    {
        SphereDistance3D sphere(1.5);
        std::string filename = outputDir + "sphere_distance_3d.mml";
        
        std::cout << "Generating Sphere Distance 3D... ";
        auto result = Serializer::SaveScalarFunc3DCartesian(
            sphere, "Sphere Distance Function (r=1.5)",
            -3.0, 3.0, gridPoints,
            -3.0, 3.0, gridPoints,
            -3.0, 3.0, gridPoints,
            filename
        );
        
        if (result.success)
            std::cout << "OK -> " << filename << "\n";
        else
            std::cout << "FAILED: " << result.message << "\n";
    }
    
    std::cout << "\nDone! Generated " << gridPoints << "^3 = " 
              << gridPoints * gridPoints * gridPoints << " points per file.\n";
    
    return 0;
}
