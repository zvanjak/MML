#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/Function.h"

#include "core/Fields.h"
#endif


using namespace MML;

void Demo4_vector_field_gravity()
{

// konstruirati multibod konfiguraciju on spot

    // GravityMultibodyConfigCart config;
    // config.AddBody(100, Vector3Cartesian{-100.0, 0.0, 0.0});
    // config.AddBody(200, Vector3Cartesian{100.0, 0.0, 0.0});
    // config.AddBody(100, Vector3Cartesian{0.0, 0.0, 100.0});

    // MultibodyGravityForceFieldCart gravity_force_field(20.0, config);

    // int numPnt = 10;
    // gravity_force_field.Serialize3DCartesian(-150.0, 150.0, numPnt, -150.0, 150.0, numPnt, -150.0, 150.0, numPnt, "multibody_gravity_field.txt");
    // auto ret = std::system("..\\tools\\visualizers\\vector_field_visualizer\\MML_VectorFieldVisualizer.exe multibody_gravity_field.txt");

}

