#include "../catch/catch.hpp"

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "basic_types/MatrixNM.h"
#endif

TEST_CASE("Test_MatrixNM", "[simple]") {
    MML::MatrixNM<2,2> mat({1.0, 0.0, 0.0, 1.0} );

}