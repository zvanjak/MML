/**
 * @file show_spherical_coord.cpp
 * @brief Visualization of the spherical coordinate system and local basis vectors
 *
 * Draws a wireframe sphere, the three Cartesian axes, a selected point on
 * the sphere, and the three unit basis vectors of the local spherical frame
 * at that point:
 *
 *   e_r     (radial, outward)        – orange
 *   e_theta (toward south pole)      – cyan
 *   e_phi   (azimuthal, eastward)    – magenta
 */

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>
#include <mml/core/CoordTransf/CoordTransfSpherical.h>
#include <mml/tools/Visualizer.h>
#endif

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace MML;

namespace {

    constexpr double kPI = 3.14159265358979323846;

    static CoordTransfSphericalToCartesian gSpherToCart;
    static CoordTransfCartesianToSpherical gCartToSpher;

    std::string ScenePath(const char* fileName)
    {
        std::filesystem::path results = GetResultFilesPath();
        std::filesystem::create_directories(results);
        return (results / fileName).string();
    }

    void WriteLine(std::ostream& out,
                   double x0, double y0, double z0,
                   double x1, double y1, double z1,
                   double radius, const char* color)
    {
        out << "LINE "
            << x0 << " " << y0 << " " << z0 << "  "
            << x1 << " " << y1 << " " << z1 << "  "
            << radius << " " << color << "\n";
    }

    void WriteVector(std::ostream& out,
                     double px, double py, double pz,
                     double vx, double vy, double vz,
                     double radius, const char* color)
    {
        out << "VECTOR "
            << px << " " << py << " " << pz << "  "
            << vx << " " << vy << " " << vz << "  "
            << radius << " " << color << "\n";
    }

    // Arrow that STARTS at (px,py,pz) and points toward (px+vx, py+vy, pz+vz)
    void WriteVectorAtPoint(std::ostream& out,
                            double px, double py, double pz,
                            double vx, double vy, double vz,
                            double radius, const char* color)
    {
        out << "VECTOR_AT_POINT "
            << px << " " << py << " " << pz << "  "
            << vx << " " << vy << " " << vz << "  "
            << radius << " " << color << "\n";
    }

    void WritePoint(std::ostream& out,
                    double px, double py, double pz,
                    double radius, const char* color)
    {
        out << "POINT " << px << " " << py << " " << pz << "  " << radius << " " << color << "\n";
    }

    // Latitude circle: constant theta, all phi in [0, 2pi)
    void WriteLatitudeCircle(std::ostream& out, double R, double theta,
                             int N, const char* color, double lineRadius = 0.8)
    {
        double z   = R * std::cos(theta);
        double rxy = R * std::sin(theta);
        for (int i = 0; i < N; ++i) {
            double phi0 = 2.0 * kPI * i       / N;
            double phi1 = 2.0 * kPI * (i + 1) / N;
            WriteLine(out,
                rxy * std::cos(phi0), rxy * std::sin(phi0), z,
                rxy * std::cos(phi1), rxy * std::sin(phi1), z,
                lineRadius, color);
        }
    }

    // Longitude arc: constant phi, theta from 0 to pi
    void WriteLongitudeArc(std::ostream& out, double R, double phi,
                           int N, const char* color, double lineRadius = 0.8)
    {
        for (int i = 0; i < N; ++i) {
            double t0 = kPI * i       / N;
            double t1 = kPI * (i + 1) / N;
            WriteLine(out,
                R * std::sin(t0) * std::cos(phi),
                R * std::sin(t0) * std::sin(phi),
                R * std::cos(t0),
                R * std::sin(t1) * std::cos(phi),
                R * std::sin(t1) * std::sin(phi),
                R * std::cos(t1),
                lineRadius, color);
        }
    }

    // Full sphere wireframe via single SPHERE_WIREFRAME command.
    // nLat/nLon = number of latitude/longitude lines;
    // segsLat/segsLon = segments per latitude circle / longitude arc.
    void WriteSphereWireframe(std::ostream& out,
                              double cx, double cy, double cz,
                              double R, int nLat, int nLon,
                              int segsLat, int segsLon,
                              double lineRadius, const char* color)
    {
        out << "SPHERE_WIREFRAME "
            << cx << " " << cy << " " << cz << "  "
            << R << "  " << nLat << " " << nLon << "  "
            << segsLat << " " << segsLon << "  "
            << lineRadius << " " << color << "\n";
    }

    // Full circle: centre (cx,cy,cz), radius R, in plane with given normal
    void WriteCircle(std::ostream& out,
                     double cx, double cy, double cz,
                     double R,
                     double nx, double ny, double nz,
                     double lineRadius, const char* color)
    {
        out << "CIRCLE "
            << cx << " " << cy << " " << cz << "  "
            << R << "  "
            << nx << " " << ny << " " << nz << "  "
            << lineRadius << " " << color << "\n";
    }

    // Three Cartesian coordinate axes: arrows FROM origin in positive direction only
    void WriteAxes(std::ostream& out, double len)
    {
        WriteVectorAtPoint(out, 0, 0, 0,  len, 0, 0,   0.375, "#FF4444");  // x – red
        WriteVectorAtPoint(out, 0, 0, 0,  0, len, 0,   0.375, "#44FF44");  // y – green
        WriteVectorAtPoint(out, 0, 0, 0,  0, 0, len,   0.375, "#4488FF");  // z – blue
    }

} // anonymous namespace

void Show_Spherical_Coord_Visualization()
{
    std::cout << "\n=== Spherical Coordinate System – Local Basis Visualization ===\n\n";

    const Real R     = 100.0;
    const Real theta = kPI / 3.0;   // 60 deg from z-axis
    const Real phi   = kPI / 4.0;   // 45 deg in xy-plane

    // -------------------------------------------------------------------
    // Cartesian position of the selected point
    // -------------------------------------------------------------------
    Vector3Spherical sphPos{ R, theta, phi };
    Vector3Cartesian cartPos = gSpherToCart.transf(sphPos);

    const double px = cartPos[0];
    const double py = cartPos[1];
    const double pz = cartPos[2];

    // -------------------------------------------------------------------
    // Local unit basis vectors (e_r, e_theta, e_phi) in Cartesian coords
    // -------------------------------------------------------------------
    Vector3Cartesian er     = gSpherToCart.getUnitBasisVec(0, sphPos);  // radial
    Vector3Cartesian etheta = gSpherToCart.getUnitBasisVec(1, sphPos);  // polar (toward south)
    Vector3Cartesian ephi   = gSpherToCart.getUnitBasisVec(2, sphPos);  // azimuthal (eastward)

    const double arrowScale  = 48.0;
    const double arrowRadius = 0.7;

    // -------------------------------------------------------------------
    // Velocity (Cartesian) and decomposition in local spherical basis
    // -------------------------------------------------------------------
    const Real vCartX = 5.0, vCartY = 35.0, vCartZ = -20.0;
    Vector3Cartesian vCart{ vCartX, vCartY, vCartZ };

    auto dot = [](const Vector3Cartesian& a, const Vector3Cartesian& b) {
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    };

    // Physical unit-basis components via dot product (valid since basis is orthonormal)
    const double vr_dot     = dot(vCart, er);
    const double vtheta_dot = dot(vCart, etheta);
    const double vphi_dot   = dot(vCart, ephi);

    // MML verification: transfInverseVecContravariant gives (dr/dt, dtheta/dt, dphi/dt)
    // Physical component = contravariant * scale factor: h_r=1, h_theta=R, h_phi=R*sin(theta)
    Vector3Spherical sphContravar = gSpherToCart.transfInverseVecContravariant(vCart, cartPos);
    const double vr_mml     = sphContravar[0];
    const double vtheta_mml = R * sphContravar[1];
    const double vphi_mml   = R * std::sin(theta) * sphContravar[2];

    // MML verification 2: gCartToSpher.transfVecContravariant – the natural forward direction
    // Source=Cartesian, Target=Spherical: directly gives (dr/dt, dtheta/dt, dphi/dt)
    Vector3Spherical sphContravar2 = gCartToSpher.transfVecContravariant(vCart, cartPos);
    const double vr_mml2     = sphContravar2[0];
    const double vtheta_mml2 = R * sphContravar2[1];
    const double vphi_mml2   = R * std::sin(theta) * sphContravar2[2];

    // Component vectors in Cartesian (for drawing and reconstruction check)
    const double cr_x = vr_dot * er[0],         cr_y = vr_dot * er[1],         cr_z = vr_dot * er[2];
    const double ct_x = vtheta_dot * etheta[0],  ct_y = vtheta_dot * etheta[1],  ct_z = vtheta_dot * etheta[2];
    const double cp_x = vphi_dot * ephi[0],      cp_y = vphi_dot * ephi[1],      cp_z = vphi_dot * ephi[2];

    // -------------------------------------------------------------------
    // Write scene file
    // -------------------------------------------------------------------
    std::string path = ScenePath("spherical_coord.mmlworld");
    std::ofstream scene(path);
    scene << std::fixed << std::setprecision(5);

    scene << "MML_WORLD_SCENE 1\n";
    scene << "TITLE Spherical coord – local basis at (r=100, theta=60 deg, phi=45 deg)\n";
    scene << "CAMERA 360 200 260\n\n";

    // Axes
    scene << "# Cartesian axes  (x=red, y=green, z=blue)  arrows from origin\n";
    WriteAxes(scene, 140.0);
    scene << "\n";

    // Sphere wireframe – dim grey
    scene << "# Sphere wireframe  R=" << R << "\n";
    WriteSphereWireframe(scene, 0, 0, 0, R, 10, 12, 60, 60, 0.175, "#334455");
    scene << "\n";

    // Highlight the meridian (full great circle at phi=45 deg) and latitude through the chosen point
    scene << "# Selected meridian great circle (phi=45 deg) and latitude (theta=60 deg)\n";
    // Great circle normal = direction perpendicular to the meridian plane: (-sin(phi), cos(phi), 0)
    WriteCircle(scene, 0, 0, 0, R,
                -std::sin(phi), std::cos(phi), 0.0,
                0.325, "#778800");
    // Latitude circle: centre on Z-axis at R*cos(theta), radius = R*sin(theta), normal = Z
    WriteCircle(scene, 0, 0, R * std::cos(theta), R * std::sin(theta),
                0, 0, 1,
                0.325, "#778800");
    scene << "\n";

    // Radial line from origin to the point
    scene << "# Radial line O -> P\n";
    WriteLine(scene, 0, 0, 0, px, py, pz, 0.25, "#AAAAAA");
    scene << "\n";

    // Point on sphere
    scene << "# Selected point on sphere\n";
    WritePoint(scene, px, py, pz, 1.25, "#FFFFFF");
    scene << "\n";

    // Local basis vectors
    scene << "# e_r     (radial, outward)       - orange\n";
    WriteVectorAtPoint(scene, px, py, pz,
                er[0] * arrowScale, er[1] * arrowScale, er[2] * arrowScale,
                arrowRadius * 0.25, "#FF8C00");

    scene << "# e_theta (toward south pole)     - cyan\n";
    WriteVectorAtPoint(scene, px, py, pz,
                etheta[0] * arrowScale, etheta[1] * arrowScale, etheta[2] * arrowScale,
                arrowRadius * 0.25, "#00D4D4");

    scene << "# e_phi   (azimuthal, eastward)   - magenta\n";
    WriteVectorAtPoint(scene, px, py, pz,
                ephi[0] * arrowScale, ephi[1] * arrowScale, ephi[2] * arrowScale,
                arrowRadius * 0.25, "#DD44FF");

    // -------------------------------------------------------------------
    // Velocity vector and spherical component arrows
    // -------------------------------------------------------------------
    scene << "\n# Velocity vector v = (5, 35, -20)  – yellow\n";
    WriteVectorAtPoint(scene, px, py, pz,
                vCartX, vCartY, vCartZ,
                arrowRadius * 1.4, "#FFEE00");

    scene << "\n# v_r component   (v · e_r)      – light orange\n";
    WriteVectorAtPoint(scene, px, py, pz, cr_x, cr_y, cr_z, arrowRadius * 0.85, "#FFC860");
    scene << "# v_theta component (v · e_theta) – light cyan\n";
    WriteVectorAtPoint(scene, px, py, pz, ct_x, ct_y, ct_z, arrowRadius * 0.85, "#44FFFF");
    scene << "# v_phi component  (v · e_phi)    – light magenta\n";
    WriteVectorAtPoint(scene, px, py, pz, cp_x, cp_y, cp_z, arrowRadius * 0.85, "#FF88FF");

    // Parallelepiped wireframe showing velocity = v_r + v_theta + v_phi
    // 8 corners: P + i*a + j*b + k*c  (a=v_r, b=v_theta, c=v_phi, i/j/k in {0,1})
    scene << "\n# Parallelepiped: P + i*v_r + j*v_theta + k*v_phi  (12 edges)\n";

    // 4 edges along v_r
    WriteLine(scene, px,           py,           pz,            px+cr_x,              py+cr_y,              pz+cr_z,              0.15, "#88CCFF");
    WriteLine(scene, px+ct_x,      py+ct_y,      pz+ct_z,       px+ct_x+cr_x,         py+ct_y+cr_y,         pz+ct_z+cr_z,         0.15, "#88CCFF");
    WriteLine(scene, px+cp_x,      py+cp_y,      pz+cp_z,       px+cp_x+cr_x,         py+cp_y+cr_y,         pz+cp_z+cr_z,         0.15, "#88CCFF");
    WriteLine(scene, px+ct_x+cp_x, py+ct_y+cp_y, pz+ct_z+cp_z,  px+ct_x+cp_x+cr_x,    py+ct_y+cp_y+cr_y,    pz+ct_z+cp_z+cr_z,    0.15, "#88CCFF");

    // 4 edges along v_theta
    WriteLine(scene, px,           py,           pz,            px+ct_x,              py+ct_y,              pz+ct_z,              0.15, "#88CCFF");
    WriteLine(scene, px+cr_x,      py+cr_y,      pz+cr_z,       px+cr_x+ct_x,         py+cr_y+ct_y,         pz+cr_z+ct_z,         0.15, "#88CCFF");
    WriteLine(scene, px+cp_x,      py+cp_y,      pz+cp_z,       px+cp_x+ct_x,         py+cp_y+ct_y,         pz+cp_z+ct_z,         0.15, "#88CCFF");
    WriteLine(scene, px+cr_x+cp_x, py+cr_y+cp_y, pz+cr_z+cp_z,  px+cr_x+cp_x+ct_x,    py+cr_y+cp_y+ct_y,    pz+cr_z+cp_z+ct_z,    0.15, "#88CCFF");

    // 4 edges along v_phi
    WriteLine(scene, px,           py,           pz,            px+cp_x,              py+cp_y,              pz+cp_z,              0.15, "#88CCFF");
    WriteLine(scene, px+cr_x,      py+cr_y,      pz+cr_z,       px+cr_x+cp_x,         py+cr_y+cp_y,         pz+cr_z+cp_z,         0.15, "#88CCFF");
    WriteLine(scene, px+ct_x,      py+ct_y,      pz+ct_z,       px+ct_x+cp_x,         py+ct_y+cp_y,         pz+ct_z+cp_z,         0.15, "#88CCFF");
    WriteLine(scene, px+cr_x+ct_x, py+cr_y+ct_y, pz+cr_z+ct_z,  px+cr_x+ct_x+cp_x,    py+cr_y+ct_y+cp_y,    pz+cr_z+ct_z+cp_z,    0.15, "#88CCFF");

    scene.close();

    // -------------------------------------------------------------------
    // Console summary
    // -------------------------------------------------------------------
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  theta = " << theta * 180.0 / kPI << " deg   phi = " << phi * 180.0 / kPI << " deg\n\n";
    std::cout << "  P      = (" << px << ", " << py << ", " << pz << ")\n";
    std::cout << "  e_r    = (" << er[0]     << ", " << er[1]     << ", " << er[2]     << ")\n";
    std::cout << "  e_theta= (" << etheta[0] << ", " << etheta[1] << ", " << etheta[2] << ")\n";
    std::cout << "  e_phi  = (" << ephi[0]   << ", " << ephi[1]   << ", " << ephi[2]   << ")\n\n";
    std::cout << "  v      = (" << vCartX << ", " << vCartY << ", " << vCartZ << ")  [Cartesian]\n\n";

    std::cout << "  Orthogonality check:\n";
    std::cout << "    e_r · e_theta = " << dot(er, etheta) << "\n";
    std::cout << "    e_r · e_phi   = " << dot(er, ephi)   << "\n";
    std::cout << "    e_theta·e_phi = " << dot(etheta, ephi) << "\n\n";

    std::cout << "  -- Velocity decomposition in local spherical basis --\n\n";
    std::cout << "  Method 1 - dot products with orthonormal unit basis:\n";
    std::cout << "    v_r     = v . e_r     = " << vr_dot     << "\n";
    std::cout << "    v_theta = v . e_theta = " << vtheta_dot << "\n";
    std::cout << "    v_phi   = v . e_phi   = " << vphi_dot   << "\n\n";

    std::cout << "  Method 2 - MML gSpherToCart.transfInverseVecContravariant (scale-factor corrected):\n";
    std::cout << "    Contravariant sph. = ("
              << sphContravar[0] << ", " << sphContravar[1] << ", " << sphContravar[2] << ")\n";
    std::cout << "      [i.e. (dr/dt, dtheta/dt, dphi/dt)]\n";
    std::cout << "    v_r     = 1            * " << sphContravar[0] << " = " << vr_mml     << "\n";
    std::cout << "    v_theta = R            * " << sphContravar[1] << " = " << vtheta_mml << "\n";
    std::cout << "    v_phi   = R*sin(theta) * " << sphContravar[2] << " = " << vphi_mml   << "\n\n";

    std::cout << "  Method 3 - MML gCartToSpher.transfVecContravariant (natural forward direction):\n";
    std::cout << "    Contravariant sph. = ("
              << sphContravar2[0] << ", " << sphContravar2[1] << ", " << sphContravar2[2] << ")\n";
    std::cout << "      [i.e. (dr/dt, dtheta/dt, dphi/dt)]\n";
    std::cout << "    v_r     = 1            * " << sphContravar2[0] << " = " << vr_mml2     << "\n";
    std::cout << "    v_theta = R            * " << sphContravar2[1] << " = " << vtheta_mml2 << "\n";
    std::cout << "    v_phi   = R*sin(theta) * " << sphContravar2[2] << " = " << vphi_mml2   << "\n\n";

    const double recX = cr_x + ct_x + cp_x;
    const double recY = cr_y + ct_y + cp_y;
    const double recZ = cr_z + ct_z + cp_z;
    std::cout << "  Reconstruction  v_r*e_r + v_theta*e_theta + v_phi*e_phi:\n";
    std::cout << "    = (" << recX << ", " << recY << ", " << recZ << ")\n";
    std::cout << "    Original v = (" << vCartX << ", " << vCartY << ", " << vCartZ << ")\n";
    std::cout << "    Error      = (" << (recX-vCartX) << ", " << (recY-vCartY) << ", " << (recZ-vCartZ) << ")\n\n";

    std::cout << "  Scene: " << path << "\n";

    if (std::getenv("MML_VISUALIZER_SAVE_ONLY")) {
        std::cout << "MML_VISUALIZER_SAVE_ONLY is set; not launching MML_WorldVisualizer.\n";
        return;
    }

    Visualizer::VisualizeWorldSceneFromFile(path);
}
