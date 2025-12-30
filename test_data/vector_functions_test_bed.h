#if !defined __MML_VECTOR_FUNCTIONS_TEST_BED_H
#define __MML_VECTOR_FUNCTIONS_TEST_BED_H

#include <string>
#include <cmath>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "base/Function.h"
#endif

namespace MML::TestBeds
{
    template<int N>
    struct TestFunctionVector
    {
        std::string _funcName;

        VectorFunction<N> _func;
        VectorN<Real, N> (*_funcDerived)(const MML::VectorN<Real, N> &, int ind);

        std::string _funcExpr;
        std::string _funcDerivedExpr;
        
        // Expected vector calculus properties for validation
        Real _expectedDiv;      // Expected divergence (NaN if position-dependent)
        Real _expectedCurlMag;  // Expected curl magnitude (NaN if position-dependent)
        bool _isSolenoidal;     // div F = 0 (incompressible)
        bool _isIrrotational;   // curl F = 0 (conservative)

        TestFunctionVector( std::string funcName,
                            VectorN<Real, N> (*f1)(const VectorN<Real, N> &), std::string funcExpr, 
                            VectorN<Real, N> (*f2)(const VectorN<Real, N> &, int ind), std::string funcDerivedExpr,
                            Real expectedDiv = std::numeric_limits<Real>::quiet_NaN(),
                            Real expectedCurlMag = std::numeric_limits<Real>::quiet_NaN(),
                            bool isSolenoidal = false, bool isIrrotational = false
                            ) : _funcName(funcName),
                                _func(f1), _funcDerived(f2), 
                                _funcExpr(funcExpr), _funcDerivedExpr(funcDerivedExpr),
                                _expectedDiv(expectedDiv), _expectedCurlMag(expectedCurlMag),
                                _isSolenoidal(isSolenoidal), _isIrrotational(isIrrotational)
        {}
    };    

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    FUNCTION 1: IDENTITY FIELD (Baseline)                              //
    ///////////////////////////////////////////////////////////////////////////////////////////
    // F(x,y,z) = (x, y, z)
    // Jacobian = Identity matrix
    // div F = 3 (constant), curl F = 0
    // Use: Simplest possible test, sanity check for Jacobian computation
    
    static MML::VectorN<Real, 3> TestVectorFunc_Identity(const VectorN<Real, 3> &v) 
    {
        return VectorN<Real, 3>{v[0], v[1], v[2]};
    }
    static MML::VectorN<Real, 3> TestVectorFunc_Identity_derived(const VectorN<Real, 3> &v, int ind) 
    { 
        // Jacobian row 'ind': d(F_i)/dx_ind for all i
        // J = I, so row 0 = (1,0,0), row 1 = (0,1,0), row 2 = (0,0,1)
        if (ind == 0) return VectorN<Real, 3>{1, 0, 0};
        else if (ind == 1) return VectorN<Real, 3>{0, 1, 0};
        else return VectorN<Real, 3>{0, 0, 1};
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    FUNCTION 2: VORTEX/ROTATION FIELD                                  //
    ///////////////////////////////////////////////////////////////////////////////////////////
    // F(x,y,z) = (-y, x, 0)
    // Rigid body rotation around z-axis with angular velocity ω=1
    // Jacobian = [[0, -1, 0], [1, 0, 0], [0, 0, 0]]
    // div F = 0 (solenoidal/incompressible), curl F = (0, 0, 2) (constant vorticity)
    // Use: Tests divergence-free fields, constant non-zero curl
    
    static MML::VectorN<Real, 3> TestVectorFunc_Vortex(const VectorN<Real, 3> &v) 
    {
        return VectorN<Real, 3>{-v[1], v[0], 0};
    }
    static MML::VectorN<Real, 3> TestVectorFunc_Vortex_derived(const VectorN<Real, 3> &v, int ind) 
    { 
        // dF/dx = (0, 1, 0), dF/dy = (-1, 0, 0), dF/dz = (0, 0, 0)
        if (ind == 0) return VectorN<Real, 3>{0, 1, 0};
        else if (ind == 1) return VectorN<Real, 3>{-1, 0, 0};
        else return VectorN<Real, 3>{0, 0, 0};
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    FUNCTION 3: RADIAL INVERSE-SQUARE FIELD                            //
    ///////////////////////////////////////////////////////////////////////////////////////////
    // F(x,y,z) = (x, y, z) / r³  where r = sqrt(x²+y²+z²)
    // This is the unit radial direction scaled by 1/r² (gravitational/Coulomb field direction)
    // div F = 0 (away from origin), curl F = 0 (irrotational)
    // Jacobian_ij = δ_ij/r³ - 3*x_i*x_j/r⁵
    // Use: Tests Laplacian fields, both solenoidal and irrotational (harmonic)
    
    static MML::VectorN<Real, 3> TestVectorFunc_RadialInvSq(const VectorN<Real, 3> &v) 
    {
        Real r2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
        if (r2 < 1e-20) return VectorN<Real, 3>{0, 0, 0};  // Regularize at origin
        Real r3_inv = 1.0 / (r2 * std::sqrt(r2));
        return VectorN<Real, 3>{v[0] * r3_inv, v[1] * r3_inv, v[2] * r3_inv};
    }
    static MML::VectorN<Real, 3> TestVectorFunc_RadialInvSq_derived(const VectorN<Real, 3> &v, int ind) 
    { 
        // dF_i/dx_j = δ_ij/r³ - 3*x_i*x_j/r⁵
        Real r2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
        if (r2 < 1e-20) return VectorN<Real, 3>{0, 0, 0};
        Real r = std::sqrt(r2);
        Real r3_inv = 1.0 / (r2 * r);
        Real r5_inv = r3_inv / r2;
        
        // Row ind of Jacobian: [dF_0/dx_ind, dF_1/dx_ind, dF_2/dx_ind]
        Real x_ind = v[ind];
        return VectorN<Real, 3>{
            (ind == 0 ? r3_inv : 0) - 3 * v[0] * x_ind * r5_inv,
            (ind == 1 ? r3_inv : 0) - 3 * v[1] * x_ind * r5_inv,
            (ind == 2 ? r3_inv : 0) - 3 * v[2] * x_ind * r5_inv
        };
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    FUNCTION 4: GRADIENT OF PRODUCT (Conservative)                     //
    ///////////////////////////////////////////////////////////////////////////////////////////
    // F(x,y,z) = ∇(xyz) = (yz, xz, xy)
    // Gradient of scalar field f=xyz, so F is conservative (curl F = 0)
    // Jacobian = [[0, z, y], [z, 0, x], [y, x, 0]] (symmetric!)
    // div F = 0, curl F = 0
    // Use: Tests conservative fields, symmetric Jacobian
    
    static MML::VectorN<Real, 3> TestVectorFunc_GradProduct(const VectorN<Real, 3> &v) 
    {
        return VectorN<Real, 3>{v[1]*v[2], v[0]*v[2], v[0]*v[1]};
    }
    static MML::VectorN<Real, 3> TestVectorFunc_GradProduct_derived(const VectorN<Real, 3> &v, int ind) 
    { 
        // dF/dx = (0, z, y), dF/dy = (z, 0, x), dF/dz = (y, x, 0)
        if (ind == 0) return VectorN<Real, 3>{0, v[2], v[1]};
        else if (ind == 1) return VectorN<Real, 3>{v[2], 0, v[0]};
        else return VectorN<Real, 3>{v[1], v[0], 0};
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    //                    FUNCTION 5: COMPLEX MIXED (Original)                               //
    ///////////////////////////////////////////////////////////////////////////////////////////
    // F(x,y,z) = (x*cos(y)*z², sin(x)*(y²+z²), exp(xy/(z²+1)))
    // Complex function with trigonometric, polynomial, and exponential terms
    // Non-trivial Jacobian, tests full numerical differentiation accuracy
    // Use: Thorough test of general vector function handling
    
    static MML::VectorN<Real, 3> TestVectorFunc_Complex(const VectorN<Real, 3> &v) 
    {
        Real x = v[0], y = v[1], z = v[2];
        return VectorN<Real, 3>{
            x * std::cos(y) * z * z,
            std::sin(x) * (y*y + z*z),
            std::exp(x*y / (z*z + 1))
        };
    }
    static MML::VectorN<Real, 3> TestVectorFunc_Complex_derived(const VectorN<Real, 3> &v, int ind) 
    { 
        Real x = v[0], y = v[1], z = v[2];
        Real zz1 = z*z + 1;
        Real expTerm = std::exp(x*y / zz1);
        
        // dF/dx: (z²cos(y), cos(x)(y²+z²), y*exp(xy/(z²+1))/(z²+1))
        // dF/dy: (-xz²sin(y), 2y*sin(x), x*exp(xy/(z²+1))/(z²+1))
        // dF/dz: (2xz*cos(y), 2z*sin(x), -2xyz*exp(xy/(z²+1))/(z²+1)²)
        
        if (ind == 0) 
            return VectorN<Real, 3>{
                z*z * std::cos(y),
                std::cos(x) * (y*y + z*z),
                y * expTerm / zz1
            };
        else if (ind == 1) 
            return VectorN<Real, 3>{
                -x * z*z * std::sin(y),
                2 * y * std::sin(x),
                x * expTerm / zz1
            };
        else 
            return VectorN<Real, 3>{
                2 * x * z * std::cos(y),
                2 * z * std::sin(x),
                -2 * x * y * z * expTerm / (zz1 * zz1)
            };
    }
    
    ///////////////////////////////////////////////////////////////////////////////////////////
    //                             TEST BED CLASS                                            //
    ///////////////////////////////////////////////////////////////////////////////////////////
    
    class VectorFunctionsTestBed
    {
    public:
        static int getNumTestFunctionVector() { return 5; }

        const static TestFunctionVector<3>& getTestFunctionVector(int i)  { return _listFuncVector3[i]; }

        const static TestFunctionVector<3>& getTestFunctionVector(const std::string &funcName)
        {
            for (int i = 0; i < getNumTestFunctionVector(); i++)
            {
                if (_listFuncVector3[i]._funcName == funcName)
                    return _listFuncVector3[i];
            }
            throw std::runtime_error("TestFunctionVector " + funcName + " not found!");
        }
    private:
        const static inline TestFunctionVector<3> _listFuncVector3[] = { 
            // Function 1: Identity - simplest baseline
            { "Identity", 
              TestVectorFunc_Identity, "F = (x, y, z)", 
              TestVectorFunc_Identity_derived, "J = I",
              3.0, 0.0, false, true },  // div=3, curl=0, not solenoidal, is irrotational
              
            // Function 2: Vortex - solenoidal with constant curl
            { "Vortex", 
              TestVectorFunc_Vortex, "F = (-y, x, 0)", 
              TestVectorFunc_Vortex_derived, "J = [[0,-1,0],[1,0,0],[0,0,0]]",
              0.0, 2.0, true, false },  // div=0, |curl|=2, is solenoidal, not irrotational
              
            // Function 3: Radial inverse-square - harmonic (both solenoidal and irrotational)
            { "Radial Inverse Square", 
              TestVectorFunc_RadialInvSq, "F = r/|r|^3", 
              TestVectorFunc_RadialInvSq_derived, "J_ij = delta_ij/r^3 - 3*x_i*x_j/r^5",
              0.0, 0.0, true, true },  // div=0, curl=0 (away from origin)
              
            // Function 4: Gradient of xyz - conservative with symmetric Jacobian
            { "Gradient of xyz", 
              TestVectorFunc_GradProduct, "F = grad(xyz) = (yz, xz, xy)", 
              TestVectorFunc_GradProduct_derived, "J = [[0,z,y],[z,0,x],[y,x,0]]",
              0.0, 0.0, true, true },  // div=0, curl=0
              
            // Function 5: Complex mixed - thorough testing
            { "Complex Mixed", 
              TestVectorFunc_Complex, "F = (x*cos(y)*z^2, sin(x)*(y^2+z^2), exp(xy/(z^2+1)))", 
              TestVectorFunc_Complex_derived, "complex Jacobian",
              std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::quiet_NaN(), 
              false, false }  // position-dependent div/curl
        };
    };

}

#endif