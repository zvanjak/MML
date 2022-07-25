#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/Vector.h"
#include "basic_types/Matrix.h"
#endif

TEST_CASE("Test_Matrix", "[simple]") {
    MML::Matrix mat(2,2, {1.0, 0.0, 0.0, 1.0} );
    
}