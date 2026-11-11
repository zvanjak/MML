#include <mml/base/Quaternions.h>

#include <cmath>
#include <iostream>
#include <limits>

using namespace MML;

void Docs_Demo_Quaternion_Scale_Contracts()
{
	const Real large = std::numeric_limits<Real>::max();
	const Quaternion q(large, -large, large, -large);
	const Quaternion unit = q.Normalized();
	const Quaternion identity = q * q.Inverse();
	if (!unit.IsUnit(Precision::EqualityTolerance) ||
		!identity.IsIdentity(Precision::EqualityTolerance))
		throw QuaternionError("Scaled quaternion normalization or inversion demo failed");
	std::cout << "Unrepresentable norm: " << q.Norm() << ", normalized: " << unit << '\n';

	const Real small = std::sqrt(PrecisionValues<Real>::QuaternionZeroThreshold) / Real(10);
	const Quaternion smallQ(small, 0, 0, 0);
	if (!(smallQ * smallQ.Inverse()).IsIdentity(Precision::EqualityTolerance))
		throw QuaternionError("Small quaternion inverse demo failed");
	std::cout << "Small quaternion inverse: " << smallQ.Inverse() << '\n';

	const Real tiny = std::numeric_limits<Real>::denorm_min() * Real(1024);
	Quaternion quotient(tiny, 0, 0, 0);
	quotient /= tiny;
	if (!quotient.IsIdentity())
		throw QuaternionError("Subnormal quaternion division demo failed");
	std::cout << "Subnormal divided by itself: " << quotient << '\n';

	bool zeroRejected = false;
	try {
		quotient /= Real(0);
	} catch (const DivisionByZeroError&) {
		zeroRejected = true;
	}
	if (!zeroRejected || !quotient.IsIdentity())
		throw QuaternionError("Quaternion division guard demo failed");

	const Quaternion rotation = Quaternion::FromEulerZYZ(Real(0.3), Constants::PI, Real(-0.7));
	for (Real sign : {Real(1), Real(-1)}) {
		const Vec3Cart angles = (rotation * sign).ToEulerZYZ();
		const Quaternion reconstructed = Quaternion::FromEulerZYZ(angles[0], angles[1], angles[2]);
		for (int axis = 0; axis < 3; ++axis) {
			Vec3Cart basis(0, 0, 0);
			basis[axis] = 1;
			const Vec3Cart expected = rotation.Rotate(basis);
			const Vec3Cart actual = reconstructed.Rotate(basis);
			for (int component = 0; component < 3; ++component)
				if (!IsNearlyEqual(actual[component], expected[component], Precision::EqualityTolerance))
					throw QuaternionError("ZYZ pole reconstruction demo failed");
		}
		std::cout << "ZYZ beta=pi, quaternion sign " << sign << ": " << angles << '\n';
	}
}
