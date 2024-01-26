#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/VectorN.h"
#include "base/VectorSpace.h"

#endif

using namespace MML;

void Demo_VectorSpaces()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                     VECTOR SPACES                             ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    
    VectorSpace<Real, VectorN<Real,3>> vecSpace;

    RealVectorSpaceN<3> vecSpaceReal3;
    ComplexVectorSpaceN<3> vecSpaceComplex3;
}