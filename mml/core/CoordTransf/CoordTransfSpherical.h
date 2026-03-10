///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        CoordTransfSpherical.h                                              ///
///  Description: Spherical coordinate transformations                                ///
///               Spherical <-> Cartesian conversions with Jacobians                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COORD_TRANSF_SPHERICAL_H
#define MML_COORD_TRANSF_SPHERICAL_H

#include "MMLBase.h"

#include "core/CoordTransf.h"


namespace MML
{
	///////////////////////////////////////////////////////////////////////////////
	// SPHERICAL COORDINATE SYSTEM CONVENTIONS
	//
	// This implementation follows the Mathematics/ISO 31-11 convention:
	//
	//   Coordinate ordering: (r, θ, φ)
	//   - r (radius):        radial distance from origin, r ≥ 0
	//   - θ (theta):         polar angle (inclination) from positive z-axis, θ ∈ [0, π]
	//   - φ (phi):           azimuthal angle in xy-plane from positive x-axis, φ ∈ [0, 2π) or (-π, π]
	//
	//   Transformation formulas:
	//     x = r sin(θ) cos(φ)
	//     y = r sin(θ) sin(φ)
	//     z = r cos(θ)
	//
	//   Inverse:
	//     r     = √(x² + y² + z²)
	//     θ     = arccos(z / r)
	//     φ     = atan2(y, x)
	//
	// ALTERNATIVE CONVENTIONS (not used here):
	//
	//   Physics convention (common in US physics textbooks):
	//     Often uses (r, φ, θ) ordering or defines θ as azimuthal and φ as polar.
	//
	//   Geodesy/Geography:
	//     Uses (latitude, longitude, altitude) where latitude = 90° - θ (colatitude).
	//
	// RATIONALE:
	//   The Math/ISO convention is used here because:
	//   - It's standard in differential geometry and tensor calculus
	//   - θ ranges match typical polar angle definitions
	//   - It aligns with spherical harmonics conventions (Ylm)
	//   - Basis vectors are naturally ordered (∂r, ∂θ, ∂φ)
	//
	///////////////////////////////////////////////////////////////////////////////

	/// @brief Spherical to Cartesian coordinate transformation (Math/ISO 31-11 convention)
	/// @details Transforms spherical (r, θ, φ) to Cartesian (x, y, z)
	///          x = r·sin(θ)·cos(φ), y = r·sin(θ)·sin(φ), z = r·cos(θ)
	///          θ is polar angle from z-axis, φ is azimuthal angle in xy-plane
	class CoordTransfSphericalToCartesian : public CoordTransfWithInverse<Vector3Spherical, Vector3Cartesian, 3>
	{
	private:
		/// q[0] = r     - radial distance
		/// q[1] = theta - inclination
		/// q[2] = phi   - azimuthal angle
		/// @brief Convert spherical to Cartesian x-coordinate
		static Real x(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]) * cos(q[2]); }
		/// @brief Convert spherical to Cartesian y-coordinate
		static Real y(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]) * sin(q[2]); }
		/// @brief Convert spherical to Cartesian z-coordinate
		static Real z(const VectorN<Real, 3>& q) { return q[0] * cos(q[1]); }

		/// q[0] = x
		/// q[1] = y
		/// q[2] = z
		/// @brief Convert Cartesian to spherical radial distance
		static Real r(const VectorN<Real, 3>& q)		 { return sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]); }
		/// @brief Convert Cartesian to spherical polar angle (inclination from z-axis)
		static Real theta(const VectorN<Real, 3>& q) { return acos(q[2] / sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2])); }
		/// @brief Convert Cartesian to spherical azimuthal angle
		static Real phi(const VectorN<Real, 3>& q)	 { return atan2(q[1], q[0]); }

		inline static ScalarFunction<3> _func[3] = { ScalarFunction<3>{x},
																								 ScalarFunction<3>{y},
																								 ScalarFunction<3>{z}
		};

		inline static ScalarFunction<3> _funcInverse[3] = { ScalarFunction<3>{r},
																												ScalarFunction<3>{theta},
																																ScalarFunction<3>{phi}
		};
	public:
		/// @brief Transform from spherical to Cartesian coordinates
		Vector3Cartesian     transf(const Vector3Spherical& q)				const override { return Vector3Cartesian{ x(q), y(q), z(q) }; }
		/// @brief Transform from Cartesian to spherical coordinates (inverse)
		Vector3Spherical     transfInverse(const Vector3Cartesian& q) const override { return Vector3Spherical{ r(q), theta(q), phi(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i)				const override { return _func[i]; }
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const override { return _funcInverse[i]; }

		/// @brief Get covariant basis vector at position in spherical coordinates
		/// @param ind Basis vector index (0=∂r, 1=∂θ, 2=∂φ)
		/// @param pos Position in spherical coordinates
		virtual Vector3Cartesian getBasisVec(int ind, const Vector3Spherical& pos) override
		{
			const Real r = pos[0];
			const Real theta = pos[1];
			const Real phi = pos[2];
			switch (ind)
			{
			case 0: return Vector3Cartesian{ sin(theta) * cos(phi),     sin(theta) * sin(phi),      cos(theta) };
			case 1: return Vector3Cartesian{ r * cos(theta) * cos(phi), r * cos(theta) * sin(phi), -r * sin(theta) };
			case 2: return Vector3Cartesian{ -r * sin(theta) * sin(phi), r * sin(theta) * cos(phi),						  REAL(0.0) };
			default:
			return Vector3Cartesian{ REAL(0.0), REAL(0.0), REAL(0.0) };
			}
		}

		/// @brief Get unit (normalized) basis vector at position
		/// @param ind Basis vector index (0=e_r, 1=e_θ, 2=e_φ)
		/// @param pos Position in spherical coordinates
		Vector3Cartesian getUnitBasisVec(int ind, const Vector3Spherical& pos)
		{
			const Real r = pos[0];
			const Real theta = pos[1];
			const Real phi = pos[2];
			switch (ind)
			{
			case 0: return Vector3Cartesian{ sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta) };
			case 1: return Vector3Cartesian{ cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta) };
			case 2: return Vector3Cartesian{ -sin(phi), cos(phi), REAL(0.0) };
			default:
			return Vector3Cartesian{ REAL(0.0), REAL(0.0), REAL(0.0) };
			}
		}

		/// @brief Get contravariant (dual) basis vector at position
		/// @param ind Basis vector index (0, 1, 2)
		/// @param pos Position in spherical coordinates
		Vector3Spherical getInverseBasisVec(int ind, const Vector3Spherical& pos) override
		{
			const Real r = pos[0];
			const Real theta = pos[1];
			const Real phi = pos[2];
			switch(ind)
			{
			case 0: return Vector3Spherical{ sin(theta) * cos(phi), r * cos(theta) * cos(phi), -r * sin(theta) * sin(phi) };
			case 1: return Vector3Spherical{ sin(theta) * sin(phi), r * cos(theta) * sin(phi),  r * sin(theta) * cos(phi) };
			case 2: return Vector3Spherical{ cos(theta)           ,-r * sin(theta)           ,                        REAL(0.0) };
			default: 
			return Vector3Spherical{ REAL(0.0), REAL(0.0), REAL(0.0) };
			}
		}
		/// @brief Get unit contravariant basis vector at position
		/// @param ind Basis vector index
		/// @param pos Position in spherical coordinates
		Vector3Spherical getInverseUnitBasisVec(int ind, const Vector3Spherical& pos)
		{
			const Real r = pos[0];
			const Real theta = pos[1];
			const Real phi = pos[2];
			switch(ind)
			{
			case 0: return Vector3Spherical{ sin(theta) * cos(phi),  cos(theta) * cos(phi), -sin(phi) };
			case 1: return Vector3Spherical{ sin(theta) * sin(phi),  cos(theta) * sin(phi),  cos(phi) };
			case 2: return Vector3Spherical{ cos(theta)           , -sin(theta)           ,  REAL(0.0) };
			default: 
			return Vector3Spherical{ REAL(0.0), REAL(0.0), REAL(0.0) };
			}
		}
	};

	/// @brief Cartesian to spherical coordinate transformation (Math/ISO 31-11 convention)
	/// @details Transforms Cartesian (x, y, z) to spherical (r, θ, φ)
	///          r = √(x²+y²+z²), θ = arccos(z/r), φ = atan2(y,x)
	class CoordTransfCartesianToSpherical : public CoordTransfWithInverse<Vector3Cartesian, Vector3Spherical, 3>
	{
	private:
		/// q[0] = x
		/// q[1] = y
		/// q[2] = z
		/// @brief Convert Cartesian to spherical radial distance
		static Real r(const VectorN<Real, 3>& q) { return sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]); }
		/// @brief Convert Cartesian to spherical polar angle
		static Real theta(const VectorN<Real, 3>& q) { return acos(q[2] / sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2])); }
		/// @brief Convert Cartesian to spherical azimuthal angle
		static Real phi(const VectorN<Real, 3>& q) { return atan2(q[1], q[0]); }

		// q[0] = r     - radial distance
		// q[1] = theta - inclination
		// q[2] = phi   - azimuthal angle
		static Real x(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]) * cos(q[2]); }
		static Real y(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]) * sin(q[2]); }
		static Real z(const VectorN<Real, 3>& q) { return q[0] * cos(q[1]); }

		inline static ScalarFunction<3> _func[3] = { ScalarFunction<3>{r},
																								 ScalarFunction<3>{theta},
																								 ScalarFunction<3>{phi}
		};

		inline static ScalarFunction<3> _funcInverse[3] = { ScalarFunction<3>{x},
																												ScalarFunction<3>{y},
																												ScalarFunction<3>{z}
		};
	public:
		/// @brief Transform from Cartesian to spherical coordinates
		Vector3Spherical     transf(const Vector3Cartesian& q) const { return Vector3Spherical{ r(q), theta(q), phi(q) }; }
		/// @brief Transform from spherical to Cartesian coordinates (inverse)
		Vector3Cartesian     transfInverse(const Vector3Spherical& q) const { return Vector3Cartesian{ x(q), y(q), z(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const { return _func[i]; }
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
	};

	/// @brief Global instance for spherical to Cartesian transformation
	static CoordTransfSphericalToCartesian      CoordTransfSpherToCart;
	/// @brief Global instance for Cartesian to spherical transformation
	static CoordTransfCartesianToSpherical      CoordTransfCartToSpher;
}

#endif
