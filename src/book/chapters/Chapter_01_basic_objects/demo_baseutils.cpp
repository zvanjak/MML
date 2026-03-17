///////////////////////////////////////////////////////////////////////////////////////////
///  File:        demo_baseutils.cpp
///  Description: Comprehensive BaseUtils demonstrations for Chapter 01
///               Shows utility functions for angles, vectors, and numerical helpers
///////////////////////////////////////////////////////////////////////////////////////////

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Vector/Vector.h"
#include "mml/base/Vector/VectorN.h"
#include "mml/base/Matrix/Matrix.h"
#include "mml/base/BaseUtils.h"
#endif

using namespace MML;

///////////////////////////////////////////////////////////////////////////////////////////
///                        Constants Overview                                           ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Constants()
{
    std::cout << "\n=== MML Constants ===\n\n";

    std::cout << "--- Mathematical Constants ---\n";
    std::cout << "PI        = " << Constants::PI << "\n";
    std::cout << "E (Euler) = " << Constants::E << "\n";
    std::cout << "SQRT2     = " << Constants::SQRT2 << "\n";
    std::cout << "SQRT3     = " << Constants::SQRT3 << "\n";
    std::cout << "LN2       = " << Constants::LN2 << "\n";
    std::cout << "LN10      = " << Constants::LN10 << "\n";

    std::cout << "\n--- Derived Values ---\n";
    std::cout << "2*PI      = " << 2 * Constants::PI << " (full circle)\n";
    std::cout << "PI/2      = " << Constants::PI / 2 << " (right angle)\n";
    std::cout << "PI/4      = " << Constants::PI / 4 << " (45 degrees)\n";
    std::cout << "180/PI    = " << 180.0 / Constants::PI << " (rad to deg factor)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Angle Utilities                                              ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Angle_Utilities()
{
    std::cout << "\n=== Angle Utilities ===\n\n";

    // Degree-Radian conversion
    std::cout << "--- Degree-Radian Conversion ---\n";
    Real deg90 = 90.0;
    Real rad90 = Utils::DegToRad(deg90);
    std::cout << "DegToRad(90) = " << rad90 << " rad (PI/2 = " << Constants::PI/2 << ")\n";
    
    Real radPI = Constants::PI;
    Real degPI = Utils::RadToDeg(radPI);
    std::cout << "RadToDeg(PI) = " << degPI << " deg\n";

    // Common conversions
    std::cout << "\n--- Common Angle Conversions ---\n";
    std::cout << "  30 deg = " << Utils::DegToRad(30.0) << " rad\n";
    std::cout << "  45 deg = " << Utils::DegToRad(45.0) << " rad\n";
    std::cout << "  60 deg = " << Utils::DegToRad(60.0) << " rad\n";
    std::cout << "  PI/6   = " << Utils::RadToDeg(Constants::PI/6) << " deg\n";
    std::cout << "  PI/4   = " << Utils::RadToDeg(Constants::PI/4) << " deg\n";
    std::cout << "  PI/3   = " << Utils::RadToDeg(Constants::PI/3) << " deg\n";

    // Degree-Minute-Second format
    std::cout << "\n--- Degree-Minute-Second Format ---\n";
    Real latitude = 40.7128;  // NYC latitude
    Real d, m, s;
    Utils::AngleDegToExplicit(latitude, d, m, s);
    std::cout << "NYC latitude: " << latitude << " deg\n";
    std::cout << "  = " << d << " deg " << m << "' " << s << "\"\n";

    // Convert back
    Real back = Utils::ExplicitToAngleDeg(40, 42, 46.08);
    std::cout << "40 deg 42' 46.08\" = " << back << " deg\n";

    // Angle normalization
    std::cout << "\n--- Angle Normalization ---\n";
    Real bigAngle = 5.0 * Constants::PI;  // 900 degrees
    Real normalized = Utils::AngleTo2PiRange(bigAngle);
    std::cout << "AngleTo2PiRange(5*PI) = " << normalized << " rad (in [0, 2*PI])\n";

    Real negAngle = -Constants::PI / 4;  // -45 degrees
    Real normPiPi = Utils::AngleToPiPiRange(negAngle);
    std::cout << "AngleToPiPiRange(-PI/4) = " << normPiPi << " rad (in [-PI, PI])\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Vector Utilities                                             ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Vector_Utilities()
{
    std::cout << "\n=== Vector Utilities ===\n\n";

    // Scalar product (dot product)
    std::cout << "--- Scalar Product ---\n";
    Vector<Real> a({1, 2, 3});
    Vector<Real> b({4, 5, 6});
    Real dot = Utils::ScalarProduct(a, b);
    std::cout << "a = " << a << "\n";
    std::cout << "b = " << b << "\n";
    std::cout << "ScalarProduct(a, b) = " << dot << " (1*4 + 2*5 + 3*6 = 32)\n";

    // Angle between vectors
    std::cout << "\n--- Angle Between Vectors ---\n";
    Vector<Real> x_axis({1, 0, 0});
    Vector<Real> y_axis({0, 1, 0});
    Vector<Real> diagonal({1, 1, 0});
    
    Real angle_xy = Utils::VectorsAngle(x_axis, y_axis);
    Real angle_xd = Utils::VectorsAngle(x_axis, diagonal);
    
    std::cout << "Angle(X-axis, Y-axis) = " << Utils::RadToDeg(angle_xy) << " deg\n";
    std::cout << "Angle(X-axis, diagonal) = " << Utils::RadToDeg(angle_xd) << " deg\n";

    // Vector projections
    std::cout << "\n--- Vector Projections ---\n";
    Vector<Real> v({3, 4, 0});
    Vector<Real> onto({1, 0, 0});  // X-axis
    
    Vector<Real> parallel = Utils::VectorProjectionParallelTo(v, onto);
    Vector<Real> perp = Utils::VectorProjectionPerpendicularTo(v, onto);
    
    std::cout << "v = " << v << "\n";
    std::cout << "Projection onto X-axis:\n";
    std::cout << "  Parallel component: " << parallel << "\n";
    std::cout << "  Perpendicular component: " << perp << "\n";
    std::cout << "  Sum (should = v): " << (parallel + perp) << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Matrix Construction Utilities                                ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_MatrixFromVector_Construction()
{
    std::cout << "\n=== Matrix Construction Utilities ===\n\n";

    Vector<Real> v({1, 2, 3});
    std::cout << "Vector v = " << v << "\n\n";

    // Row matrix from vector
    std::cout << "--- Row Matrix from Vector ---\n";
    Matrix<Real> rowMat = Utils::RowMatrixFromVector(v);
    std::cout << "RowMatrixFromVector(v):\n" << rowMat << "\n";

    // Column matrix from vector
    std::cout << "--- Column Matrix from Vector ---\n";
    Matrix<Real> colMat = Utils::ColumnMatrixFromVector(v);
    std::cout << "ColumnMatrixFromVector(v):\n" << colMat << "\n";

    // Diagonal matrix from vector
    std::cout << "--- Diagonal Matrix from Vector ---\n";
    Matrix<Real> diagMat = Utils::DiagonalMatrixFromVector(v);
    std::cout << "DiagonalMatrixFromVector(v):\n" << diagMat << "\n";

    // Outer product (tensor product)
    std::cout << "--- Outer Product ---\n";
    Vector<Real> u({1, 2});
    Vector<Real> w({3, 4, 5});
    Matrix<Real> outer = Utils::OuterProduct(u, w);
    std::cout << "u = " << u << "\n";
    std::cout << "w = " << w << "\n";
    std::cout << "OuterProduct(u, w):\n" << outer << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Equality and Comparison                                      ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Equality_Comparison()
{
    std::cout << "\n=== Numerical Equality ===\n\n";

    std::cout << "Floating-point comparison requires tolerance!\n\n";

    // Real number comparison
    std::cout << "--- Real Number Comparison ---\n";
    Real a = 1.0;
    Real b = 1.0 + 1e-10;
    Real c = 1.0 + 1e-6;
    
    std::cout << "a = 1.0\n";
    std::cout << "b = 1.0 + 1e-10\n";
    std::cout << "c = 1.0 + 1e-6\n\n";
    
    std::cout << "a == b (exact): " << (a == b ? "true" : "false") << "\n";
    std::cout << "AreEqual(a, b, 1e-9): " << (Utils::AreEqual(a, b, 1e-9) ? "true" : "false") << "\n";
    std::cout << "AreEqual(a, c, 1e-9): " << (Utils::AreEqual(a, c, 1e-9) ? "true" : "false") << "\n";
    std::cout << "AreEqual(a, c, 1e-5): " << (Utils::AreEqual(a, c, 1e-5) ? "true" : "false") << "\n";

    // Complex number comparison
    std::cout << "\n--- Complex Number Comparison ---\n";
    Complex z1(1.0, 2.0);
    Complex z2(1.0 + 1e-10, 2.0 - 1e-10);
    
    std::cout << "z1 = " << z1 << "\n";
    std::cout << "z2 = " << z2 << " (slightly different)\n";
    std::cout << "AreEqual(z1, z2): " << (Utils::AreEqual(z1, z2) ? "true" : "false") << "\n";

    // Vector comparison
    std::cout << "\n--- Vector Comparison ---\n";
    Vector<Real> v1({1.0, 2.0, 3.0});
    Vector<Real> v2({1.0 + 1e-10, 2.0, 3.0 - 1e-10});
    
    std::cout << "v1 = " << v1 << "\n";
    std::cout << "v2 = " << v2 << "\n";
    std::cout << "AreEqual(v1, v2): " << (Utils::AreEqual(v1, v2) ? "true" : "false") << "\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Practical Applications                                       ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_Utils_Applications()
{
    std::cout << "\n=== Practical Applications ===\n\n";

    // Navigation: GPS coordinates
    std::cout << "--- Navigation: GPS Coordinate Conversion ---\n";
    // New York City
    Real nyc_lat_deg = 40.7128;
    Real nyc_lon_deg = -74.0060;
    
    Real lat_rad = Utils::DegToRad(nyc_lat_deg);
    Real lon_rad = Utils::DegToRad(nyc_lon_deg);
    
    std::cout << "NYC: (" << nyc_lat_deg << " N, " << -nyc_lon_deg << " W)\n";
    std::cout << "  In radians: (" << lat_rad << ", " << lon_rad << ")\n";

    // Physics: Force decomposition
    std::cout << "\n--- Physics: Force Decomposition ---\n";
    Vector<Real> gravity({0, -9.81, 0});  // Downward
    Vector<Real> incline({1, 0.5, 0});    // 30 degree incline direction
    incline = incline * (1.0 / incline.NormL2());  // Normalize
    
    Vector<Real> F_parallel = Utils::VectorProjectionParallelTo(gravity, incline);
    Vector<Real> F_normal = Utils::VectorProjectionPerpendicularTo(gravity, incline);
    
    std::cout << "Gravity: " << gravity << " m/s^2\n";
    std::cout << "Incline direction: " << incline << "\n";
    std::cout << "Force along incline: " << F_parallel << " (|F| = " << F_parallel.NormL2() << ")\n";
    std::cout << "Normal force: " << F_normal << " (|F| = " << F_normal.NormL2() << ")\n";

    // Graphics: Lighting calculation
    std::cout << "\n--- Graphics: Diffuse Lighting (Lambert) ---\n";
    Vector<Real> light_dir({1, 1, 1});
    light_dir = light_dir * (1.0 / light_dir.NormL2());  // Normalize
    
    Vector<Real> surface_normal({0, 1, 0});  // Pointing up
    
    Real cos_angle = Utils::ScalarProduct(light_dir, surface_normal);
    Real intensity = std::max(0.0, cos_angle);  // Clamp negative
    
    std::cout << "Light direction: " << light_dir << "\n";
    std::cout << "Surface normal: " << surface_normal << "\n";
    std::cout << "cos(angle) = " << cos_angle << "\n";
    std::cout << "Diffuse intensity = " << intensity << " (" << intensity * 100 << "%)\n";
}

///////////////////////////////////////////////////////////////////////////////////////////
///                        Main Demo Entry Point                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

void Demo_BaseUtils()
{
    std::cout << "\n";
    std::cout << "***********************************************************************\n";
    std::cout << "****                    BASE UTILITIES IN MML                      ****\n";
    std::cout << "****         Constants, Angles, Vector Ops, and Helpers            ****\n";
    std::cout << "***********************************************************************\n";

    Demo_Constants();
    Demo_Angle_Utilities();
    Demo_Vector_Utilities();
    Demo_MatrixFromVector_Construction();
    Demo_Equality_Comparison();
    Demo_Utils_Applications();
}
