#include <catch2/catch_all.hpp>
#include "../../TestPrecision.h"

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/core/CoordTransf/CoordTransfBase.h>
#include <mml/core/CoordTransf/CoordTransfSpherical.h>
#include <mml/core/Fields/FieldOperations.h>
#include <mml/core/Fields/Fields.h>
#include <mml/core/Fields/TensorFieldAdapters.h>
#endif


using namespace MML;

namespace MML::Tests::Core::CoordTransfTensorTests
{
	class ConnectionLikeRank3Field : public ITensorField3<3>
	{
	public:
		ConnectionLikeRank3Field() : ITensorField3<3>(1, 2) {}

		bool TransformsAsTensor() const override { return false; }

		Tensor3<3> operator()(const VectorN<Real, 3>& pos) const override
		{
			Tensor3<3> ret(2, 1);
			ret(0, 1, 1) = REAL(1.0);
			return ret;
		}

		Real Component(int i, int j, int k, const VectorN<Real, 3>& pos) const override
		{
			return (*this)(pos)(i, j, k);
		}
	};

	TEST_CASE("TransformedTensorField3 rejects connection-like non-tensor fields", "[CoordTransf][TensorField][Christoffel]")
	{
		ConnectionLikeRank3Field connection;
		CoordTransfCartesianToSpherical coordTransf;

		REQUIRE_FALSE(connection.TransformsAsTensor());
		REQUIRE_THROWS_AS((TransformedTensorField3<Vector3Cartesian, Vector3Spherical, 3>(connection, coordTransf)), std::invalid_argument);
	}
} // namespace MML::Tests::Core::CoordTransfTensorTests
