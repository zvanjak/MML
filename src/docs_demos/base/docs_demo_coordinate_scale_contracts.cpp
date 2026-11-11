#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/base/Vector/Vector.h>
#include <mml/base/Vector/VectorTypes2D.h>
#include <mml/base/Vector/VectorTypes3D.h>
#include <mml/base/Vector/VectorTypes4D.h>
#endif

#include <iostream>
#include <stdexcept>

void Docs_Demo_CoordinateScaleContracts()
{
    using namespace MML;
    Vector<int> values{1, 2, 1};
    values.erase(0);
    values.eraseValue(1);
    if (values != Vector<int>{2})
        throw std::runtime_error("Integer vector erase contract failed");

    Vector2Polar polar(1, 0);
    polar += Vector2Polar(1, Constants::PI / Real{2});
    const auto difference = Vector2Polar(1, 0) - Vector2Polar(1, Real{1e-10});
    if (!IsNearlyEqual(polar.R(), std::sqrt(Real{2})) || !(difference.R() > 0))
        throw std::runtime_error("Polar arithmetic contract failed");
    Vector3Cartesian huge(Real{1e200}, 0, 0), diagonal(Real{1e200}, Real{1e200}, 0);
    if (huge.IsParallelTo(diagonal) || !IsNearlyEqual(huge.AngleToVector(diagonal), Constants::PI / Real{4}))
        throw std::runtime_error("Scaled direction contract failed");

    const Vector4Minkowski time{Real{5e200}, Real{3e200}, 0, 0};
    const Vector4Minkowski mixed{Real{1e200}, Real{1e200}, 1, 0};
    if (!time.IsTimelike() || !IsNearlyEqual(time.Norm() / Real{1e200}, Real{4}) ||
        !mixed.IsSpacelike() || mixed.Norm() != -1)
        throw std::runtime_error("Scaled Minkowski contract failed");
    std::cout << "Polar sum radius: " << polar.R() << "; close-angle difference: " << difference.R()
              << "\nHuge-vector angle: " << huge.AngleToVector(diagonal)
              << "\nMinkowski norm: " << time.Norm() << "; mixed-scale norm: " << mixed.Norm() << '\n';
}
