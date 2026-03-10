///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        CoordTransfCylindrical.h                                            ///
///  Description: Cylindrical coordinate transformations                              ///
///               Cylindrical <-> Cartesian conversions with Jacobians                ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COORD_TRANSF_CYLINDRICAL_H
#define MML_COORD_TRANSF_CYLINDRICAL_H

#include "MMLBase.h"

#include "core/CoordTransf.h"


namespace MML
{
	///////////////////////////////////////////////////////////////////////////////
	// CYLINDRICAL COORDINATE SYSTEM CONVENTIONS
	//
	// This implementation uses the standard (r, φ, z) ordering:
	//
	//   Coordinate ordering: (r, φ, z)
	//   - r (radius):      distance from the z-axis (symmetry axis), r ≥ 0
	//   - φ (phi):         azimuthal angle in xy-plane from positive x-axis
	//                      φ ∈ (-π, π]  (range of atan2)
	//   - z (height):      same as Cartesian z, z ∈ (-∞, +∞)
	//
	//   Transformation formulas:
	//     x = r cos(φ)
	//     y = r sin(φ)
	//     z = z
	//
	//   Inverse:
	//     r = √(x² + y²)
	//     φ = atan2(y, x)
	//     z = z
	//
	//   Angle units:   Radians throughout
	//   Handedness:    Right-handed (r, φ, z) — consistent with (x, y, z)
	//   Singularity:   r = 0 (z-axis) — φ is undefined
	//
	// See also: CoordTransfSpherical.h for spherical conventions
	///////////////////////////////////////////////////////////////////////////////

	/// @brief Cylindrical to Cartesian coordinate transformation
	/// @details Transforms cylindrical coordinates (r, φ, z) to Cartesian (x, y, z)
	///          x = r·cos(φ), y = r·sin(φ), z = z
	class CoordTransfCylindricalToCartesian : public CoordTransfWithInverse<Vector3Cylindrical, Vector3Cartesian, 3>
	{
	private:
		/// q1 = r   - distance from symmetry axis
		/// q2 = phi - angle to symmetry axis
		/// q3 = z   - z
		/// @brief Convert cylindrical to Cartesian x-coordinate
		static Real x(const VectorN<Real, 3>& q) { return q[0] * cos(q[1]); }
		static Real y(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]); }
		/// @brief Convert cylindrical to Cartesian z-coordinate (unchanged)
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
		/// @brief Transform from Cartesian to cylindrical coordinates (inverse)
		Vector3Cylindrical   transfInverse(const Vector3Cartesian& q) const { return Vector3Cylindrical{ r(q), phi(q), z(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i)        const { return _func[i]; }
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
	};

	/// @brief Cartesian to cylindrical coordinate transformation
	/// @details Transforms Cartesian coordinates (x, y, z) to cylindrical (r, φ, z)
	///          r = √(x²+y²), φ = atan2(y,x), z = z
	class CoordTransfCartesianToCylindrical : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cylindrical, 3>
	{
	private:
		/// q[0] = x
		/// q[1] = y
		/// q[2] = z
		/// @brief Convert Cartesian to cylindrical radial distance
		static Real r(const VectorN<Real, 3>& q) { return sqrt(q[0] * q[0] + q[1] * q[1]); }
		/// @brief Convert Cartesian to cylindrical azimuthal angle
		static Real phi(const VectorN<Real, 3>& q) { return atan2(q[1], q[0]); }
		/// @brief Convert Cartesian to cylindrical z-coordinate (unchanged)
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
		/// @brief Transform from cylindrical to Cartesian coordinates (inverse)
		Vector3Cartesian   transfInverse(const Vector3Cylindrical& q) const { return Vector3Cartesian{ x(q), y(q), z(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i)        const { return _func[i]; }
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
	};

	/// @brief Global instance for Cartesian to cylindrical transformation
	static CoordTransfCartesianToCylindrical    CoordTransfCartToCyl;
	/// @brief Global instance for cylindrical to Cartesian transformation
	static CoordTransfCylindricalToCartesian    CoordTransfCylToCart;

}

#endif
