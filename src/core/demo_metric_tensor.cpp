// TODO - BIG!!! - Demo metric tensor 
#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/Tensor.h"
#include "core/CoordTransf.h"
#include "core/MetricTensor.h"

#include "base/Geometry3D.h"
#endif

using namespace MML;

void Demo_Metric_Tensors_comparison()
{
    MetricTensorCartesian<3> metricCart;
    MetricTensorSpherical metricSpher;
    MetricTensorCylindrical metricCyl;

    // using instanced object    
    CoordTransfSphericalToCartesian coordTransfSpherToCart;
    MetricTensorFromCoordTransf<Vector3Spherical, Vector3Cartesian, 3> metricSpherFromCart(coordTransfSpherToCart);
    // using static 
    MetricTensorFromCoordTransf<Vector3Spherical, Vector3Cartesian, 3>   metricSpherFromTransf(CoordTransfSpherToCart);
    MetricTensorFromCoordTransf<Vector3Cylindrical, Vector3Cartesian, 3> metricCylFromTransf(CoordTransfCylToCart);

    Vector3Cartesian   pos(1.0, 2.0, -1.0);
    Vector3Spherical   posSpher = CoordTransfSpherToCart.transf(pos);
    Vector3Cylindrical posCyl   = CoordTransfCylToCart.transf(pos);

    auto cart_metric = metricCart(pos);

    auto spher_metric = metricSpher(posSpher);

    auto cyl_metric = metricCyl(posCyl);

    auto spher_metric2 = metricSpherFromTransf(posSpher);

    auto cyl_metric2 = metricCylFromTransf(posCyl);

    std::cout << "Point (Cart): " << pos << std::endl;
    std::cout << "Cartesian Metric Tensor:" << cart_metric << std::endl;
    std::cout << "Spherical Metric Tensor:" << spher_metric << std::endl;
    std::cout << "Cylindrical Metric Tensor:" << cyl_metric << std::endl;
    std::cout << "Spherical Metric Tensor (from transf):" << spher_metric2 << std::endl;
    std::cout << "Cylindrical Metric Tensor (from transf):" << cyl_metric2 << std::endl;
}

void Demo_Metric_Tensors()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                       METRIC TENSORS                          ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    
    Tensor2<3> t2(2,0);
    Tensor3<3> t3(1,2);
    Tensor4<3> t4(3,1);

    Demo_Metric_Tensors_comparison();
}