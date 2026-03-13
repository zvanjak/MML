///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        generate_serializer_reference_files.cpp                             ///
///  Description: Generates reference files for serializer golden file tests         ///
///                                                                                   ///
///  Run this once to generate all reference files in:                               ///
///    tests/tools/serialization_test_data/                                           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                    ///
///////////////////////////////////////////////////////////////////////////////////////////

#include "MMLBase.h"
#include "tools/Serializer.h"
#include "interfaces/IFunction.h"
#include "base/InterpolatedFunction.h"

#include <iostream>
#include <filesystem>
#include <cmath>

using namespace MML;

namespace MML::Tests::Tools::GenerateSerializerReferenceFiles {

//////////////////////////////////////////////////////////////////////////////////
// Test Functions - Same as in serializer_tests.cpp for consistency
//////////////////////////////////////////////////////////////////////////////////

class SinFunction : public IRealFunction {
public:
    Real operator()(Real x) const override { return std::sin(x); }
};

class CosFunction : public IRealFunction {
public:
    Real operator()(Real x) const override { return std::cos(x); }
};

class SquareFunction : public IRealFunction {
public:
    Real operator()(Real x) const override { return x * x; }
};

class ScalarProd2D : public IScalarFunction<2> {
public:
    Real operator()(const VectorN<Real, 2>& p) const override {
        return p[0] * p[1];
    }
};

class ScalarSum3D : public IScalarFunction<3> {
public:
    Real operator()(const VectorN<Real, 3>& p) const override {
        return p[0] + p[1] + p[2];
    }
};

class IdentityVectorField2D : public IVectorFunction<2> {
public:
    VectorN<Real, 2> operator()(const VectorN<Real, 2>& p) const override {
        return p;
    }
};

class IdentityVectorField3D : public IVectorFunction<3> {
public:
    VectorN<Real, 3> operator()(const VectorN<Real, 3>& p) const override {
        return p;
    }
};

class CircleCurve2D : public IRealToVectorFunction<2> {
public:
    VectorN<Real, 2> operator()(Real t) const override {
        return VectorN<Real, 2>{std::cos(t), std::sin(t)};
    }
};

class HelixCurve3D : public IRealToVectorFunction<3> {
public:
    VectorN<Real, 3> operator()(Real t) const override {
        return VectorN<Real, 3>{std::cos(t), std::sin(t), static_cast<Real>(t / (2 * Constants::PI))};
    }
};

class SphereSurfaceNM : public IVectorFunctionNM<2, 3> {
public:
    VectorN<Real, 3> operator()(const VectorN<Real, 2>& p) const override {
        Real u = p[0];
        Real w = p[1];
        return VectorN<Real, 3>{
            std::cos(u) * std::sin(w),
            std::sin(u) * std::sin(w),
            std::cos(w)
        };
    }
};

} // namespace MML::Tests::Tools::GenerateSerializerReferenceFiles

//////////////////////////////////////////////////////////////////////////////////
// Main - Generate all reference files
//////////////////////////////////////////////////////////////////////////////////

using namespace MML::Tests::Tools::GenerateSerializerReferenceFiles;

int main() {
    std::string baseDir = "tests/tools/serialization_test_data/";
    
    // Ensure directory exists
    std::filesystem::create_directories(baseDir);
    
    std::cout << "Generating serializer reference files in: " << baseDir << "\n\n";

    //============================================================================
    // 1. REAL_FUNCTION_EQUALLY_SPACED
    //============================================================================
    {
        SinFunction sinFunc;
        
        // Simple: 5 points
        auto r1 = Serializer::SaveRealFuncEquallySpaced(sinFunc, "sin(x)", 
            0.0, Constants::PI, 5, baseDir + "realfunc_equally_spaced_simple.mml");
        std::cout << "realfunc_equally_spaced_simple: " << (r1.success ? "OK" : "FAILED") << "\n";
        
        // Complex: 13 points (full cycle)
        auto r2 = Serializer::SaveRealFuncEquallySpaced(sinFunc, "sin(x)", 
            0.0, 2*Constants::PI, 13, baseDir + "realfunc_equally_spaced_complex.mml");
        std::cout << "realfunc_equally_spaced_complex: " << (r2.success ? "OK" : "FAILED") << "\n";
    }

    //============================================================================
    // 2. REAL_FUNCTION (specified points)
    //============================================================================
    {
        CosFunction cosFunc;
        
        // Simple: 4 points
        Vector<Real> pts_simple(std::vector<Real>{0.0, 1.0, 2.0, 3.0});
        auto r1 = Serializer::SaveRealFunc(cosFunc, "cos(x) custom", pts_simple, 
            baseDir + "realfunc_specified_simple.mml");
        std::cout << "realfunc_specified_simple: " << (r1.success ? "OK" : "FAILED") << "\n";
        
        // Complex: 8 points (non-uniform spacing)
        Vector<Real> pts_complex(std::vector<Real>{0.0, 0.5, 1.0, 1.5, 2.5, 3.5, 5.0, 6.28});
        auto r2 = Serializer::SaveRealFunc(cosFunc, "cos(x) custom", pts_complex, 
            baseDir + "realfunc_specified_complex.mml");
        std::cout << "realfunc_specified_complex: " << (r2.success ? "OK" : "FAILED") << "\n";
    }

    //============================================================================
    // 3. MULTI_REAL_FUNCTION
    //============================================================================
    {
        SinFunction sinFunc;
        CosFunction cosFunc;
        std::vector<IRealFunction*> funcs = {&sinFunc, &cosFunc};
        std::vector<std::string> legend = {"sin(x)", "cos(x)"};
        
        // Simple: 5 points
        auto r1 = Serializer::SaveRealMultiFunc(funcs, "Trig Functions", legend,
            0.0, Constants::PI, 5, baseDir + "multi_realfunc_simple.mml");
        std::cout << "multi_realfunc_simple: " << (r1.success ? "OK" : "FAILED") << "\n";
        
        // Complex: 11 points
        auto r2 = Serializer::SaveRealMultiFunc(funcs, "Trig Functions", legend,
            0.0, 2*Constants::PI, 11, baseDir + "multi_realfunc_complex.mml");
        std::cout << "multi_realfunc_complex: " << (r2.success ? "OK" : "FAILED") << "\n";
    }

    //============================================================================
    // 4. PARAMETRIC_CURVE_CARTESIAN_2D
    //============================================================================
    {
        CircleCurve2D circle;
        
        // Simple: 9 points (quarter turns)
        bool r1 = Serializer::SaveParamCurveCartesian2DResult(circle, "Unit Circle",
            0.0, 2*Constants::PI, 9, baseDir + "paramcurve2d_simple.mml").success;
        std::cout << "paramcurve2d_simple: " << (r1 ? "OK" : "FAILED") << "\n";
        
        // Complex: 25 points
        bool r2 = Serializer::SaveParamCurveCartesian2DResult(circle, "Unit Circle",
            0.0, 2*Constants::PI, 25, baseDir + "paramcurve2d_complex.mml").success;
        std::cout << "paramcurve2d_complex: " << (r2 ? "OK" : "FAILED") << "\n";
    }

    //============================================================================
    // 5. PARAMETRIC_CURVE_CARTESIAN_3D
    //============================================================================
    {
        HelixCurve3D helix;
        
        // Simple: 9 points (one turn)
        bool r1 = Serializer::SaveParamCurveCartesian3DResult(helix, "Helix",
            0.0, 2*Constants::PI, 9, baseDir + "paramcurve3d_simple.mml").success;
        std::cout << "paramcurve3d_simple: " << (r1 ? "OK" : "FAILED") << "\n";
        
        // Complex: 25 points (two turns)
        bool r2 = Serializer::SaveParamCurveCartesian3DResult(helix, "Helix",
            0.0, 4*Constants::PI, 25, baseDir + "paramcurve3d_complex.mml").success;
        std::cout << "paramcurve3d_complex: " << (r2 ? "OK" : "FAILED") << "\n";
    }

    //============================================================================
    // 6. PARAMETRIC_SURFACE_CARTESIAN
    //============================================================================
    {
        SphereSurfaceNM sphere;
        
        // Simple: 5x5 grid
        auto r1 = Serializer::SaveParametricSurface(sphere, "Unit Sphere",
            0.0, 2*Constants::PI, 5,
            0.0, Constants::PI, 5,
            baseDir + "paramsurface_simple.mml");
        std::cout << "paramsurface_simple: " << (r1.success ? "OK" : "FAILED") << "\n";
        
        // Complex: 9x9 grid
        auto r2 = Serializer::SaveParametricSurface(sphere, "Unit Sphere",
            0.0, 2*Constants::PI, 9,
            0.0, Constants::PI, 9,
            baseDir + "paramsurface_complex.mml");
        std::cout << "paramsurface_complex: " << (r2.success ? "OK" : "FAILED") << "\n";
    }

    //============================================================================
    // 7. SCALAR_FUNCTION_CARTESIAN_2D
    //============================================================================
    {
        ScalarProd2D scalarFunc;
        
        // Simple: 4x4 grid
        auto r1 = Serializer::SaveScalarFunc2DCartesian(scalarFunc, "z=x*y",
            0.0, 3.0, 4,
            0.0, 3.0, 4,
            baseDir + "scalar2d_simple.mml");
        std::cout << "scalar2d_simple: " << (r1.success ? "OK" : "FAILED") << "\n";
        
        // Complex: 7x7 grid
        auto r2 = Serializer::SaveScalarFunc2DCartesian(scalarFunc, "z=x*y",
            -2.0, 4.0, 7,
            -2.0, 4.0, 7,
            baseDir + "scalar2d_complex.mml");
        std::cout << "scalar2d_complex: " << (r2.success ? "OK" : "FAILED") << "\n";
    }

    //============================================================================
    // 8. SCALAR_FUNCTION_CARTESIAN_3D
    //============================================================================
    {
        ScalarSum3D scalarFunc;
        
        // Simple: 3x3x3 grid
        auto r1 = Serializer::SaveScalarFunc3DCartesian(scalarFunc, "w=x+y+z",
            0.0, 2.0, 3,
            0.0, 2.0, 3,
            0.0, 2.0, 3,
            baseDir + "scalar3d_simple.mml");
        std::cout << "scalar3d_simple: " << (r1.success ? "OK" : "FAILED") << "\n";
        
        // Complex: 4x4x4 grid
        auto r2 = Serializer::SaveScalarFunc3DCartesian(scalarFunc, "w=x+y+z",
            -1.0, 2.0, 4,
            -1.0, 2.0, 4,
            -1.0, 2.0, 4,
            baseDir + "scalar3d_complex.mml");
        std::cout << "scalar3d_complex: " << (r2.success ? "OK" : "FAILED") << "\n";
    }

    //============================================================================
    // 9. VECTOR_FIELD_2D_CARTESIAN
    //============================================================================
    {
        IdentityVectorField2D vecFunc;
        
        // Simple: 3x3 grid
        auto r1 = Serializer::SaveVectorFunc2DCartesian(vecFunc, "Identity 2D",
            -1.0, 1.0, 3,
            -1.0, 1.0, 3,
            baseDir + "vector2d_simple.mml");
        std::cout << "vector2d_simple: " << (r1.success ? "OK" : "FAILED") << "\n";
        
        // Complex: 5x5 grid
        auto r2 = Serializer::SaveVectorFunc2DCartesian(vecFunc, "Identity 2D",
            -2.0, 2.0, 5,
            -2.0, 2.0, 5,
            baseDir + "vector2d_complex.mml");
        std::cout << "vector2d_complex: " << (r2.success ? "OK" : "FAILED") << "\n";
    }

    //============================================================================
    // 10. VECTOR_FIELD_3D_CARTESIAN
    //============================================================================
    {
        IdentityVectorField3D vecFunc;
        
        // Simple: 2x2x2 grid
        auto r1 = Serializer::SaveVectorFunc3DCartesian(vecFunc, "Identity 3D",
            -1.0, 1.0, 2,
            -1.0, 1.0, 2,
            -1.0, 1.0, 2,
            baseDir + "vector3d_simple.mml");
        std::cout << "vector3d_simple: " << (r1.success ? "OK" : "FAILED") << "\n";
        
        // Complex: 3x3x3 grid
        auto r2 = Serializer::SaveVectorFunc3DCartesian(vecFunc, "Identity 3D",
            -2.0, 2.0, 3,
            -2.0, 2.0, 3,
            -2.0, 2.0, 3,
            baseDir + "vector3d_complex.mml");
        std::cout << "vector3d_complex: " << (r2.success ? "OK" : "FAILED") << "\n";
    }

    //============================================================================
    // 11. VECTOR_FIELD_SPHERICAL
    //============================================================================
    {
        IdentityVectorField3D vecFunc;
        
        // Simple: 2x3x4 grid
        auto r1 = Serializer::SaveVectorFuncSpherical(vecFunc, "Spherical Field",
            1.0, 2.0, 2,
            0.0, Constants::PI, 3,
            0.0, 2*Constants::PI, 4,
            baseDir + "vectorspherical_simple.mml");
        std::cout << "vectorspherical_simple: " << (r1.success ? "OK" : "FAILED") << "\n";
        
        // Complex: 3x5x6 grid
        auto r2 = Serializer::SaveVectorFuncSpherical(vecFunc, "Spherical Field",
            0.5, 2.5, 3,
            0.0, Constants::PI, 5,
            0.0, 2*Constants::PI, 6,
            baseDir + "vectorspherical_complex.mml");
        std::cout << "vectorspherical_complex: " << (r2.success ? "OK" : "FAILED") << "\n";
    }

    //============================================================================
    // 12. PARTICLE_SIMULATION_DATA_2D
    //============================================================================
    {
        // Simple: 2 balls, 5 steps
        int numBalls = 2;
        std::vector<std::vector<Pnt2Cart>> positions(numBalls);
        for (int b = 0; b < numBalls; ++b) {
            for (int step = 0; step < 5; ++step) {
                positions[b].push_back(Pnt2Cart(10.0 + b * 20.0 + step * 2.0, 
                                                 50.0 + b * 10.0 + step * 3.0));
            }
        }
        std::vector<std::string> colors = {"red", "blue"};
        std::vector<Real> radii = {5.0, 7.5};
        
        auto r1 = Serializer::SaveParticleSimulation2D(baseDir + "particle2d_simple.mml", 
            numBalls, 100.0, 100.0, positions, colors, radii, 0.01, 1);
        std::cout << "particle2d_simple: " << (r1.success ? "OK" : "FAILED") << "\n";
        
        // Complex: 3 balls, 10 steps
        numBalls = 3;
        positions.clear();
        positions.resize(numBalls);
        for (int b = 0; b < numBalls; ++b) {
            for (int step = 0; step < 10; ++step) {
                Real angle = step * 0.3 + b * Constants::PI / 3;
                positions[b].push_back(Pnt2Cart(50.0 + 20.0 * std::cos(angle), 
                                                 50.0 + 20.0 * std::sin(angle)));
            }
        }
        colors = {"green", "yellow", "purple"};
        radii = {3.0, 4.0, 5.0};
        
        auto r2 = Serializer::SaveParticleSimulation2D(baseDir + "particle2d_complex.mml",
            numBalls, 100.0, 100.0, positions, colors, radii, 0.02, 1);
        std::cout << "particle2d_complex: " << (r2.success ? "OK" : "FAILED") << "\n";
    }

    //============================================================================
    // 13. PARTICLE_SIMULATION_DATA_3D
    //============================================================================
    {
        // Simple: 2 balls, 4 steps
        int numBalls = 2;
        std::vector<std::vector<Pnt3Cart>> positions(numBalls);
        for (int b = 0; b < numBalls; ++b) {
            for (int step = 0; step < 4; ++step) {
                positions[b].push_back(Pnt3Cart(10.0 + b * 15.0, 
                                                 20.0 + step * 5.0,
                                                 30.0 + b * step));
            }
        }
        std::vector<std::string> colors = {"white", "black"};
        std::vector<Real> radii = {2.0, 3.0};
        
        auto r1 = Serializer::SaveParticleSimulation3D(baseDir + "particle3d_simple.mml",
            numBalls, 50.0, 50.0, 50.0, positions, colors, radii, 0.01, 1);
        std::cout << "particle3d_simple: " << (r1.success ? "OK" : "FAILED") << "\n";
        
        // Complex: 3 balls, 8 steps
        numBalls = 3;
        positions.clear();
        positions.resize(numBalls);
        for (int b = 0; b < numBalls; ++b) {
            for (int step = 0; step < 8; ++step) {
                Real t = step * 0.5;
                positions[b].push_back(Pnt3Cart(
                    25.0 + 10.0 * std::cos(t + b * 2.0),
                    25.0 + 10.0 * std::sin(t + b * 2.0),
                    25.0 + 5.0 * t
                ));
            }
        }
        colors = {"cyan", "magenta", "orange"};
        radii = {1.5, 2.0, 2.5};
        
        auto r2 = Serializer::SaveParticleSimulation3D(baseDir + "particle3d_complex.mml",
            numBalls, 50.0, 50.0, 75.0, positions, colors, radii, 0.02, 1);
        std::cout << "particle3d_complex: " << (r2.success ? "OK" : "FAILED") << "\n";
    }

    std::cout << "\n=== Generation Complete ===\n";
    std::cout << "Generated 26 reference files in " << baseDir << "\n";
    
    return 0;
}
