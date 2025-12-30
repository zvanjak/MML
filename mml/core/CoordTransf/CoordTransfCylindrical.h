///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        CoordTransfCylindrical.h                                            ///
///  Description: Cylindrical coordinate transformations                              ///
///               Cylindrical <-> Cartesian conversions with Jacobians                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COORD_TRANSF_CYLINDRICAL_H
#define MML_COORD_TRANSF_CYLINDRICAL_H

#include "MMLBase.h"

#include "core/CoordTransf.h"


namespace MML
{

	class CoordTransfCylindricalToCartesian : public CoordTransfWithInverse<Vector3Cylindrical, Vector3Cartesian, 3>
	{
	private:
		// q1 = r   - distance from symmetry axis
		// q2 = phi - angle to symmetry axis
		// q3 = z   - z
		static Real x(const VectorN<Real, 3>& q) { return q[0] * cos(q[1]); }
		static Real y(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]); }
		static Real z(const VectorN<Real, 3>& q) { return q[2]; }

		// q[0] = x
		// q[1] = y
		// q[2] = z
		static Real r(const VectorN<Real, 3>& q) { return sqrt(q[0] * q[0] + q[1] * q[1]); }
		static Real phi(const VectorN<Real, 3>& q) { return atan2(q[1], q[0]); }
		// z-coordinate is the same, so we'll use for inverse the same function as for forward transformation

		inline static ScalarFunction<3> _func[3] = { ScalarFunction<3>{x},
																								 ScalarFunction<3>{y},
																								 ScalarFunction<3>{z}
		};

		inline static ScalarFunction<3> _funcInverse[3] = { ScalarFunction<3>{r},
																												ScalarFunction<3>{phi},
																												ScalarFunction<3>{z}
		};
	public:
		Vector3Cartesian     transf(const Vector3Cylindrical& q)      const { return Vector3Cartesian{ x(q), y(q), z(q) }; }
		Vector3Cylindrical   transfInverse(const Vector3Cartesian& q) const { return Vector3Cylindrical{ r(q), phi(q), z(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i)        const { return _func[i]; }
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
	};

	class CoordTransfCartesianToCylindrical : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cylindrical, 3>
	{
	private:
		// q[0] = x
		// q[1] = y
		// q[2] = z
		static Real r(const VectorN<Real, 3>& q) { return sqrt(q[0] * q[0] + q[1] * q[1]); }
		static Real phi(const VectorN<Real, 3>& q) { return atan2(q[1], q[0]); }
		static Real z(const VectorN<Real, 3>& q) { return q[2]; }

		// q1 = r   - distance from symmetry axis
		// q2 = phi - angle to symmetry axis
		// q3 = z   - z
		static Real x(const VectorN<Real, 3>& q) { return q[0] * cos(q[1]); }
		static Real y(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]); }

		inline static ScalarFunction<3> _func[3] = { ScalarFunction<3>{r},
																								 ScalarFunction<3>{phi},
																								 ScalarFunction<3>{z}
		};

		inline static ScalarFunction<3> _funcInverse[3] = { ScalarFunction<3>{x},
																												ScalarFunction<3>{y},
																												ScalarFunction<3>{z}
		};
	public:
		Vector3Cylindrical transf(const Vector3Cartesian& q)          const { return Vector3Cylindrical{ r(q), phi(q), z(q) }; }
		Vector3Cartesian   transfInverse(const Vector3Cylindrical& q) const { return Vector3Cartesian{ x(q), y(q), z(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i)        const { return _func[i]; }
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
	};

	static CoordTransfCartesianToCylindrical    CoordTransfCartToCyl;
	static CoordTransfCylindricalToCartesian    CoordTransfCylToCart;

}

#endif
