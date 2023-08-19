#include <iostream>
#include <iomanip>
#include <cmath>

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "core/Tensor.h"
#include "basic_types/MetricTensor.h"
#endif

using namespace MML;

void Demo_Tensors()
{
    std::cout << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                           TENSORS                             ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    
    Tensor2<3> t2(1,1);
    Tensor3<3> t3(1,1);
    Tensor4<3> t4(1,1);

    MetricTensorCartesian<3> mtc;
    MetricTensorSpherical mts;
    MetricTensorCylindrical mtcyl;
}