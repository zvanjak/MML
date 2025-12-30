#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"
#include "../TestMatchers.h"

#include "../../mml/core/Derivation/DerivationTensorField.h"
#include "../../mml/core/Derivation/Jacobians.h"
#include "../../mml/interfaces/ITensorField.h"
#include "../../mml/interfaces/IFunction.h"
#include "../../mml/base/VectorN.h"
#include "../../mml/base/MatrixNM.h"

using namespace MML;
using namespace MML::Testing;

/**************************************************************************
 * TEST TENSOR FIELDS WITH KNOWN ANALYTICAL DERIVATIVES
 **************************************************************************/

// Simple 2D tensor field: T_ij(x,y) = x^i * y^j
// Derivatives: dT_ij/dx = i*x^(i-1)*y^j, dT_ij/dy = j*x^i*y^(j-1)
class TestTensorField2_Polynomial : public ITensorField2<2>
{
public:
    TestTensorField2_Polynomial() : ITensorField2<2>(2, 0) {}

    Real Component(int i, int j, const VectorN<Real, 2>& pos) const override
    {
        Real x = pos[0];
        Real y = pos[1];
        return std::pow(x, i) * std::pow(y, j);
    }

    Tensor2<2> operator()(const VectorN<Real, 2>& pos) const override
    {
        Tensor2<2> ret(2, 0);
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                ret(i, j) = Component(i, j, pos);
        return ret;
    }

    Real ComponentDerivX(int i, int j, const VectorN<Real, 2>& pos) const
    {
        if (i == 0) return REAL(0.0);
        Real x = pos[0];
        Real y = pos[1];
        return i * std::pow(x, i - 1) * std::pow(y, j);
    }

    Real ComponentDerivY(int i, int j, const VectorN<Real, 2>& pos) const
    {
        if (j == 0) return REAL(0.0);
        Real x = pos[0];
        Real y = pos[1];
        return j * std::pow(x, i) * std::pow(y, j - 1);
    }
};

// 3D tensor field: T_ij(x,y,z) = sin(x*i) * cos(y*j) * z
class TestTensorField2_Trigonometric : public ITensorField2<3>
{
public:
    TestTensorField2_Trigonometric() : ITensorField2<3>(2, 0) {}

    Real Component(int i, int j, const VectorN<Real, 3>& pos) const override
    {
        Real x = pos[0];
        Real y = pos[1];
        Real z = pos[2];
        return std::sin(x * i) * std::cos(y * j) * z;
    }

    Tensor2<3> operator()(const VectorN<Real, 3>& pos) const override
    {
        Tensor2<3> ret(2, 0);
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                ret(i, j) = Component(i, j, pos);
        return ret;
    }

    Real ComponentDerivX(int i, int j, const VectorN<Real, 3>& pos) const
    {
        if (i == 0) return REAL(0.0);
        Real x = pos[0];
        Real y = pos[1];
        Real z = pos[2];
        return i * std::cos(x * i) * std::cos(y * j) * z;
    }

    Real ComponentDerivY(int i, int j, const VectorN<Real, 3>& pos) const
    {
        if (j == 0) return REAL(0.0);
        Real x = pos[0];
        Real y = pos[1];
        Real z = pos[2];
        return -j * std::sin(x * i) * std::sin(y * j) * z;
    }

    Real ComponentDerivZ(int i, int j, const VectorN<Real, 3>& pos) const
    {
        Real x = pos[0];
        Real y = pos[1];
        return std::sin(x * i) * std::cos(y * j);
    }
};

/**************************************************************************
 * TEST VECTOR FUNCTIONS FOR JACOBIAN TESTS
 **************************************************************************/

// 2D rotation: [x', y'] = [x*cos(theta) - y*sin(theta), x*sin(theta) + y*cos(theta)]
class TestRotationVectorFunction : public IVectorFunction<2>
{
    Real theta = REAL(0.5);
public:
    VectorN<Real, 2> operator()(const VectorN<Real, 2>& x) const override
    {
        VectorN<Real, 2> result;
        result[0] = x[0] * std::cos(theta) - x[1] * std::sin(theta);
        result[1] = x[0] * std::sin(theta) + x[1] * std::cos(theta);
        return result;
    }

    MatrixNM<Real, 2, 2> AnalyticalJacobian() const
    {
        MatrixNM<Real, 2, 2> jac;
        jac(0, 0) = std::cos(theta);
        jac(0, 1) = -std::sin(theta);
        jac(1, 0) = std::sin(theta);
        jac(1, 1) = std::cos(theta);
        return jac;
    }
};

// 3D quadratic: [x^2, y^2, z^2]
class TestQuadraticVectorFunction : public IVectorFunction<3>
{
public:
    VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const override
    {
        VectorN<Real, 3> result;
        result[0] = x[0] * x[0];
        result[1] = x[1] * x[1];
        result[2] = x[2] * x[2];
        return result;
    }

    MatrixNM<Real, 3, 3> AnalyticalJacobian(const VectorN<Real, 3>& pos) const
    {
        MatrixNM<Real, 3, 3> jac;
        jac(0, 0) = REAL(2.0) * pos[0];
        jac(0, 1) = REAL(0.0);
        jac(0, 2) = REAL(0.0);
        jac(1, 0) = REAL(0.0);
        jac(1, 1) = REAL(2.0) * pos[1];
        jac(1, 2) = REAL(0.0);
        jac(2, 0) = REAL(0.0);
        jac(2, 1) = REAL(0.0);
        jac(2, 2) = REAL(2.0) * pos[2];
        return jac;
    }
};

/**************************************************************************
 * TEST CASES FOR TENSOR FIELD DERIVATIVES
 **************************************************************************/

TEST_CASE("TensorField2 - NDer1Partial with polynomial tensor", "[tensor_field][NDer1]")
{
    TestTensorField2_Polynomial tensor;
    VectorN<Real, 2> pos{REAL(2.0), REAL(3.0)};

    SECTION("Derivatives with respect to x (index 0)")
    {
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                Real numerical = Derivation::NDer1Partial(tensor, i, j, 0, pos);
                Real analytical = tensor.ComponentDerivX(i, j, pos);
                REQUIRE(std::abs(numerical - analytical) < 1e-4);
            }
        }
    }

    SECTION("Derivatives with respect to y (index 1)")
    {
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                Real numerical = Derivation::NDer1Partial(tensor, i, j, 1, pos);
                Real analytical = tensor.ComponentDerivY(i, j, pos);
                REQUIRE(std::abs(numerical - analytical) < 1e-4);
            }
        }
    }
}

TEST_CASE("TensorField2 - NDer2Partial with trigonometric tensor", "[tensor_field][NDer2]")
{
    TestTensorField2_Trigonometric tensor;
    VectorN<Real, 3> pos{REAL(0.5), REAL(1.0), REAL(2.0)};

    SECTION("Derivatives with respect to x")
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                Real numerical = Derivation::NDer2Partial(tensor, i, j, 0, pos);
                Real analytical = tensor.ComponentDerivX(i, j, pos);
                REQUIRE(std::abs(numerical - analytical) < 1e-7);
            }
        }
    }

    SECTION("Derivatives with respect to y")
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                Real numerical = Derivation::NDer2Partial(tensor, i, j, 1, pos);
                Real analytical = tensor.ComponentDerivY(i, j, pos);
                REQUIRE(std::abs(numerical - analytical) < 1e-7);
            }
        }
    }

    SECTION("Derivatives with respect to z")
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                Real numerical = Derivation::NDer2Partial(tensor, i, j, 2, pos);
                Real analytical = tensor.ComponentDerivZ(i, j, pos);
                REQUIRE(std::abs(numerical - analytical) < 1e-7);
            }
        }
    }
}

TEST_CASE("TensorField2 - NDer4Partial with polynomial tensor", "[tensor_field][NDer4]")
{
    TestTensorField2_Polynomial tensor;
    VectorN<Real, 2> pos{REAL(1.5), REAL(2.5)};

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            Real numerical_x = Derivation::NDer4Partial(tensor, i, j, 0, pos);
            Real analytical_x = tensor.ComponentDerivX(i, j, pos);
            REQUIRE(std::abs(numerical_x - analytical_x) < 1e-9);

            Real numerical_y = Derivation::NDer4Partial(tensor, i, j, 1, pos);
            Real analytical_y = tensor.ComponentDerivY(i, j, pos);
            REQUIRE(std::abs(numerical_y - analytical_y) < 1e-9);
        }
    }
}

TEST_CASE("TensorField - High accuracy verification", "[tensor_field][accuracy]")
{
    TestTensorField2_Polynomial tensor;
    VectorN<Real, 2> pos{REAL(2.0), REAL(3.0)};

    SECTION("NDer1Partial accuracy")
    {
        Real analytical = tensor.ComponentDerivX(1, 1, pos);
        Real numerical = Derivation::NDer1Partial(tensor, 1, 1, 0, pos);
        REQUIRE(std::abs(numerical - analytical) < 1e-3);
    }
    
    SECTION("NDer2Partial accuracy")
    {
        Real analytical = tensor.ComponentDerivX(1, 1, pos);
        Real numerical = Derivation::NDer2Partial(tensor, 1, 1, 0, pos);
        REQUIRE(std::abs(numerical - analytical) < 1e-6);
    }
    
    SECTION("NDer4Partial accuracy")
    {
        Real analytical = tensor.ComponentDerivX(1, 1, pos);
        Real numerical = Derivation::NDer4Partial(tensor, 1, 1, 0, pos);
        REQUIRE(std::abs(numerical - analytical) < 1e-9);
    }
}

/**************************************************************************
 * TEST CASES FOR JACOBIAN CALCULATIONS
 **************************************************************************/

TEST_CASE("Jacobian - 2D rotation matrix", "[jacobian][2d][rotation]")
{
    TestRotationVectorFunction func;
    VectorN<Real, 2> pos{REAL(1.0), REAL(2.0)};

    auto numerical_jac = Derivation::calcJacobian(func, pos);
    auto analytical_jac = func.AnalyticalJacobian();

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            REQUIRE(std::abs(numerical_jac(i, j) - analytical_jac(i, j)) < 1e-9);
        }
    }
}

TEST_CASE("Jacobian - 3D diagonal quadratic", "[jacobian][3d][quadratic]")
{
    TestQuadraticVectorFunction func;
    VectorN<Real, 3> pos{REAL(1.0), REAL(2.0), REAL(3.0)};

    auto numerical_jac = Derivation::calcJacobian(func, pos);
    auto analytical_jac = func.AnalyticalJacobian(pos);

    SECTION("Diagonal elements are 2x, 2y, 2z")
    {
        REQUIRE(std::abs(numerical_jac(0, 0) - REAL(2.0) * pos[0]) < 1e-9);
        REQUIRE(std::abs(numerical_jac(1, 1) - REAL(2.0) * pos[1]) < 1e-9);
        REQUIRE(std::abs(numerical_jac(2, 2) - REAL(2.0) * pos[2]) < 1e-9);
    }

    SECTION("Off-diagonal elements are zero")
    {
        REQUIRE(std::abs(numerical_jac(0, 1)) < 1e-9);
        REQUIRE(std::abs(numerical_jac(0, 2)) < 1e-9);
        REQUIRE(std::abs(numerical_jac(1, 0)) < 1e-9);
        REQUIRE(std::abs(numerical_jac(1, 2)) < 1e-9);
        REQUIRE(std::abs(numerical_jac(2, 0)) < 1e-9);
        REQUIRE(std::abs(numerical_jac(2, 1)) < 1e-9);
    }

    SECTION("Full comparison with analytical")
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                REQUIRE(std::abs(numerical_jac(i, j) - analytical_jac(i, j)) < 1e-9);
            }
        }
    }
}

TEST_CASE("Jacobian - Multiple evaluation points", "[jacobian][robustness]")
{
    TestQuadraticVectorFunction func;
    
    VectorN<Real, 3> test_points[] = {
        {REAL(0.5), REAL(1.0), REAL(1.5)},
        {REAL(1.0), REAL(2.0), REAL(3.0)},
        {REAL(2.0), REAL(2.0), REAL(2.0)},
        {-REAL(1.0), REAL(1.0), REAL(2.0)}
    };

    for (const auto& pos : test_points)
    {
        auto numerical_jac = Derivation::calcJacobian(func, pos);
        auto analytical_jac = func.AnalyticalJacobian(pos);

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                REQUIRE(std::abs(numerical_jac(i, j) - analytical_jac(i, j)) < 1e-9);
            }
        }
    }
}

TEST_CASE("Jacobian - Identity function has identity Jacobian", "[jacobian][identity]")
{
    class IdentityFunction : public IVectorFunction<3>
    {
    public:
        VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const override
        {
            return x;
        }
    };

    IdentityFunction func;
    VectorN<Real, 3> pos{REAL(1.0), REAL(2.0), REAL(3.0)};

    auto jac = Derivation::calcJacobian(func, pos);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            if (i == j)
                REQUIRE(std::abs(jac(i, j) - REAL(1.0)) < 1e-10);
            else
                REQUIRE(std::abs(jac(i, j)) < 1e-10);
        }
    }
}

TEST_CASE("Jacobian - Linear transformation", "[jacobian][linear]")
{
    // F(x,y) = [2x + 3y, 4x - y]
    class LinearFunction : public IVectorFunction<2>
    {
    public:
        VectorN<Real, 2> operator()(const VectorN<Real, 2>& x) const override
        {
            VectorN<Real, 2> result;
            result[0] = REAL(2.0) * x[0] + REAL(3.0) * x[1];
            result[1] = REAL(4.0) * x[0] - x[1];
            return result;
        }
    };

    LinearFunction func;
    VectorN<Real, 2> pos{REAL(1.0), REAL(2.0)};

    auto jac = Derivation::calcJacobian(func, pos);

    REQUIRE(std::abs(jac(0, 0) - REAL(2.0)) < 1e-10);
    REQUIRE(std::abs(jac(0, 1) - REAL(3.0)) < 1e-10);
    REQUIRE(std::abs(jac(1, 0) - REAL(4.0)) < 1e-10);
    REQUIRE(std::abs(jac(1, 1) - (-REAL(1.0))) < 1e-10);
}
