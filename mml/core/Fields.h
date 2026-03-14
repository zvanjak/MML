///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Fields.h                                                            ///
///  Description: Scalar, vector, and tensor field classes                            ///
///               Field representations for physical simulations                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FIELDS_H
#define MML_FIELDS_H


#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "interfaces/IFunction.h"

#include "base/Vector/VectorN.h"
#include "base/Vector/VectorTypes.h"
#include "base/Geometry/Geometry3D.h"

#include "core/FieldOperations.h"

#endif


namespace MML::Fields
{
	///////////////////////////////////////////////////////////////////////////////
	// FIELD CONVENTIONS
	//
	//   Sign convention:   Potentials are defined as Φ(r) = C/r where
	//                      C > 0 for repulsive (Coulomb like-charge)
	//                      C < 0 for attractive (gravity: C = -GM)
	//                      Force field:  F = -∇Φ  (negative gradient of potential)
	//
	//   Coordinate suffixes:
	//     *Cart   → Cartesian (x, y, z):     position = (x[0], x[1], x[2])
	//     *Spher  → Spherical (r, θ, φ):      uses Math/ISO 31-11 convention
	//     *Cyl    → Cylindrical (r, φ, z):    standard cylindrical ordering
	//     See CoordTransfSpherical.h, CoordTransfCylindrical.h for details.
	//
	//   Inverse-square law: F ∝ C·r̂/r² (Gravity, Coulomb)
	//                       Potential:  Φ ∝ C/r
	//
	//   Units:             All quantities are dimensionless — the caller is
	//                      responsible for consistent unit systems.
	//
	// See also: FieldOperations.h, FieldAnalyzers.h
	///////////////////////////////////////////////////////////////////////////////

	////////////////////             INVERSE RADIAL FIELD                /////////////////
	/// @brief Inverse radial potential field functions (Φ ∝ 1/r)
	/// @details Used for gravitational and Coulomb potentials: Φ(r) = -GM/r, Φ(r) = kQ/r
	
	/// @brief Cartesian inverse radial potential Φ = 1/|r|
	static Real InverseRadialPotentialFieldCart(const VectorN<Real, 3>& x) { return 1.0 / x.NormL2(); }
	/// @brief Cartesian inverse radial potential with constant Φ = C/|r|
	static Real InverseRadialPotentialFieldCart(Real constant, const VectorN<Real, 3>& x) { return constant / x.NormL2(); }
	/// @brief Spherical inverse radial potential Φ = 1/r
	static Real InverseRadialPotentialFieldSpher(const VectorN<Real, 3>& x) { return 1.0 / x[0]; }
	/// @brief Spherical inverse radial potential with constant Φ = C/r
	static Real InverseRadialPotentialFieldSpher(Real constant, const VectorN<Real, 3>& x) { return constant / x[0]; }
	/// @brief Cylindrical inverse radial potential Φ = 1/√(r²+z²)
	static Real InverseRadialPotentialFieldCyl(const VectorN<Real, 3>& x) { return 1.0 / sqrt(x[0] * x[0] + x[2] * x[2]); }
	/// @brief Cylindrical inverse radial potential with constant Φ = C/√(r²+z²)
	static Real InverseRadialPotentialFieldCyl(Real constant, const VectorN<Real, 3>& x) { return constant / sqrt(x[0] * x[0] + x[2] * x[2]); }

	/// @brief Force field functions (F = -∇Φ for inverse square law)
	/// @details F(r) ∝ r/|r|³, used for gravity and Coulomb force
	
	/// @brief Cartesian force field F = -r/|r|³
	/// @note Returns VectorN<Real,3> for VectorFunction<3> interface compatibility
	static VectorN<Real, 3> InverseRadialPotentialForceFieldCart(const VectorN<Real, 3>& x)
	{
		return (-1) * x / std::pow((x).NormL2(), 3);
	}
	/// @brief Cartesian force field with constant F = -C·r/|r|³
	static VectorN<Real, 3> InverseRadialPotentialForceFieldCart(Real constant, const VectorN<Real, 3>& x)
	{
		return -constant * x / std::pow((x).NormL2(), 3);
	}
	/// @brief Spherical force field F = (-1/r²)e_r
	static VectorN<Real, 3> InverseRadialPotentialForceFieldSph(const VectorN<Real, 3>& x)
	{
		return VectorN<Real, 3>{-1 / (x[0] * x[0]), 0.0, 0.0};
	}
	/// @brief Spherical force field with constant F = (-C/r²)e_r
	static VectorN<Real, 3> InverseRadialPotentialForceFieldSph(Real constant, const VectorN<Real, 3>& x)
	{
		return VectorN<Real, 3>{-constant / (x[0] * x[0]), 0.0, 0.0};
	}

	/// @brief Inverse radial scalar field in Cartesian coordinates
	/// @details Represents gravitational or Coulomb potential Φ(r) = C/|r|
	class InverseRadialFieldCart : public IScalarFunction<3>
	{
	protected:
		Real _constant;
	public:
		InverseRadialFieldCart() : _constant(-1.0) {}
		InverseRadialFieldCart(Real constant) : _constant(constant) {}

		Real operator()(const VectorN<Real, 3>& x) const { return _constant * Fields::InverseRadialPotentialFieldCart(x); }
	};
	/// @brief Inverse radial scalar field in spherical coordinates
	/// @details Represents gravitational or Coulomb potential Φ(r) = C/r
	class InverseRadialFieldSpher : public IScalarFunction<3>
	{
	protected:
		Real _constant;
	public:
		InverseRadialFieldSpher() : _constant(-1.0) {}
		InverseRadialFieldSpher(Real constant) : _constant(constant) {}

		Real operator()(const VectorN<Real, 3>& x) const { return _constant * InverseRadialPotentialFieldSpher(x); }
	};

	/// @brief Inverse radial force field in Cartesian coordinates
	/// @details Represents gravitational or Coulomb force F(r) = C·r/|r|³
	class InverseRadialForceFieldCart : public IVectorFunction<3>
	{
	private:
		Real _constant;
	public:
		InverseRadialForceFieldCart() : _constant(-1.0) {}
		InverseRadialForceFieldCart(Real constant) : _constant(constant) {}

		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const { return _constant * Fields::InverseRadialPotentialForceFieldCart(x); }
	};
	/// @brief Inverse radial force field in spherical coordinates
	/// @details Represents gravitational or Coulomb force F(r) = (C/r²)e_r
	class InverseRadialForceFieldSpher : public IVectorFunction<3>
	{
	private:
		Real _constant;
	public:
		InverseRadialForceFieldSpher() : _constant(-1.0) {}
		InverseRadialForceFieldSpher(Real constant) : _constant(constant) {}

		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const { return _constant * InverseRadialPotentialForceFieldSph(x); }
	};

	////////////////////             MAGNETIC FIELDS                     /////////////////

	/// @brief Magnetic field from an infinite straight current-carrying wire (Biot-Savart)
	/// @details B = (μ₀I)/(2πr) in azimuthal direction (right-hand rule)
	/// Wire along z-axis, field circles around it in x-y plane
	/// @note Uses Cartesian coordinates; constant = μ₀I/(2π)
	class InfiniteLineCurrentField : public IVectorFunction<3>
	{
	private:
		Real _constant;  ///< μ₀I/(2π) - magnetic constant times current over 2π
	public:
		/// @brief Default constructor with unit constant
		InfiniteLineCurrentField() : _constant(1.0) {}
		/// @brief Constructor with specified constant (μ₀I/(2π))
		/// @param constant The value μ₀I/(2π) where I is the current
		InfiniteLineCurrentField(Real constant) : _constant(constant) {}

		/// @brief Compute magnetic field at position x
		/// @param x Position in Cartesian coordinates (wire along z-axis)
		/// @return Magnetic field vector B = (constant/r²)(-y, x, 0)
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const
		{
			Real r_sq = x[0] * x[0] + x[1] * x[1];  // Distance squared from z-axis
			if (r_sq < 1e-15) return VectorN<Real, 3>{0.0, 0.0, 0.0};  // Avoid singularity
			Real factor = _constant / r_sq;
			// B direction is tangential: (-y, x, 0) / r, magnitude is constant/r
			return VectorN<Real, 3>{-factor * x[1], factor * x[0], 0.0};
		}
	};

	/// @brief Magnetic field from an infinite wire along an arbitrary Line3D
	/// @details Generalized Biot-Savart: B = (μ₀I)/(2πr) perpendicular to both wire and radial
	/// @note constant = μ₀I/(2π) where I is the current
	class InfiniteLineCurrentFieldGeneral : public IVectorFunction<3>
	{
	private:
		Real _constant;    ///< μ₀I/(2π)
		Line3D _line;      ///< The wire geometry
	public:
		/// @brief Constructor with Line3D and optional constant
		/// @param line The line representing the wire (point + direction)
		/// @param constant The value μ₀I/(2π), defaults to 1.0
		InfiniteLineCurrentFieldGeneral(const Line3D& line, Real constant = 1.0)
			: _constant(constant), _line(line) {}

		/// @brief Compute magnetic field at position x
		/// @param x Position in Cartesian coordinates
		/// @return Magnetic field vector B = (constant/r) * (d × r̂)
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const
		{
			// Convert position to Point3D
			Pnt3Cart pos(x[0], x[1], x[2]);
			
			// Find closest point on wire to field point
			Pnt3Cart closest = _line.NearestPointOnLine(pos);
			
			// Radial vector from wire to field point
			Vec3Cart radial(closest, pos);
			Real r_sq = radial.NormL2() * radial.NormL2();
			
			if (r_sq < 1e-15) return VectorN<Real, 3>{0.0, 0.0, 0.0};  // On the wire
			
			// B direction: d × r̂ (right-hand rule), magnitude: constant/r
			Vec3Cart B_dir = VectorProduct(_line.Direction(), radial.GetAsUnitVector());
			Real factor = _constant / std::sqrt(r_sq);
			
			return VectorN<Real, 3>{factor * B_dir.X(), factor * B_dir.Y(), factor * B_dir.Z()};
		}
	};

	/// @brief Magnetic field from an infinite current-carrying plane
	/// @details B = (μ₀K)/2, uniform field parallel to plane (perpendicular to current)
	/// Plane at z=0 with surface current density K in x-direction
	/// Field is ±(μ₀K/2) in y-direction (sign depends on which side of plane)
	/// @note constant = μ₀K/2 where K is surface current density
	class InfinitePlaneCurrentField : public IVectorFunction<3>
	{
	private:
		Real _constant;        ///< μ₀K/2 - field magnitude
		VectorN<Real, 3> _dir; ///< Field direction (perpendicular to current, parallel to plane)
	public:
		/// @brief Default constructor: plane at z=0, current in x, field in ±y
		InfinitePlaneCurrentField() : _constant(1.0), _dir{0.0, 1.0, 0.0} {}
		/// @brief Constructor with specified constant
		/// @param constant The value μ₀K/2 where K is surface current density
		InfinitePlaneCurrentField(Real constant) : _constant(constant), _dir{0.0, 1.0, 0.0} {}
		/// @brief Full constructor with custom field direction
		/// @param constant The value μ₀K/2
		/// @param fieldDirection Unit vector in field direction (should be normalized)
		InfinitePlaneCurrentField(Real constant, const VectorN<Real, 3>& fieldDirection)
			: _constant(constant), _dir(fieldDirection) {}

		/// @brief Compute magnetic field at position x
		/// @param x Position in Cartesian coordinates
		/// @return Uniform magnetic field, sign flips across z=0 plane
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const
		{
			// Field direction flips across the plane (z=0)
			Real sign = (x[2] >= 0.0) ? 1.0 : -1.0;
			return sign * _constant * _dir;
		}
	};

	////////////////////             MATHEMATICAL VECTOR FIELDS          /////////////////

	/// @brief Rigid body rotation (vortex) field: v = ω × r
	/// @details Represents solid body rotation about an axis. Has constant curl = 2ω, zero divergence.
	/// Default: rotation about z-axis with angular velocity ω
	class VortexField : public IVectorFunction<3>
	{
	private:
		VectorN<Real, 3> _omega;  ///< Angular velocity vector (axis and magnitude)
	public:
		/// @brief Default constructor: unit rotation about z-axis
		VortexField() : _omega{0.0, 0.0, 1.0} {}
		/// @brief Constructor with angular velocity magnitude (rotation about z-axis)
		VortexField(Real omega) : _omega{0.0, 0.0, omega} {}
		/// @brief Full constructor with arbitrary rotation axis
		VortexField(const VectorN<Real, 3>& omega) : _omega(omega) {}

		/// @brief Compute velocity at position x: v = ω × r
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const
		{
			// Cross product: ω × r
			return VectorN<Real, 3>{
				_omega[1] * x[2] - _omega[2] * x[1],
				_omega[2] * x[0] - _omega[0] * x[2],
				_omega[0] * x[1] - _omega[1] * x[0]
			};
		}
	};

	/// @brief Potential (irrotational) vortex field: v = (Γ/2πr)ê_θ
	/// @details Irrotational everywhere except at origin (singularity). Zero curl, zero divergence.
	/// Classic example of multiply-connected domain with non-zero circulation.
	/// @note constant = Γ/(2π) where Γ is circulation
	class PotentialVortexField : public IVectorFunction<3>
	{
	private:
		Real _constant;           ///< Γ/(2π) - circulation over 2π
		VectorN<Real, 3> _axis;   ///< Vortex axis (default z)
	public:
		/// @brief Default constructor: unit circulation about z-axis
		PotentialVortexField() : _constant(1.0), _axis{0.0, 0.0, 1.0} {}
		/// @brief Constructor with circulation constant (about z-axis)
		PotentialVortexField(Real constant) : _constant(constant), _axis{0.0, 0.0, 1.0} {}
		/// @brief Full constructor with arbitrary axis
		PotentialVortexField(Real constant, const VectorN<Real, 3>& axis)
			: _constant(constant), _axis(axis / std::sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2])) {}

		/// @brief Compute velocity at position x
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const
		{
			// For z-axis: v = (Γ/2πr²)(-y, x, 0)
			// General: project x onto plane perpendicular to axis, then tangent direction
			Real ax = _axis[0], ay = _axis[1], az = _axis[2];
			
			// Component of x along axis
			Real proj = x[0]*ax + x[1]*ay + x[2]*az;
			
			// Perpendicular component (radial in the rotation plane)
			VectorN<Real, 3> r_perp{x[0] - proj*ax, x[1] - proj*ay, x[2] - proj*az};
			Real r_sq = r_perp[0]*r_perp[0] + r_perp[1]*r_perp[1] + r_perp[2]*r_perp[2];
			
			if (r_sq < 1e-15) return VectorN<Real, 3>{0.0, 0.0, 0.0};  // On the axis
			
			// Tangent direction: axis × r_perp (normalized by r²)
			Real factor = _constant / r_sq;
			return VectorN<Real, 3>{
				factor * (ay * r_perp[2] - az * r_perp[1]),
				factor * (az * r_perp[0] - ax * r_perp[2]),
				factor * (ax * r_perp[1] - ay * r_perp[0])
			};
		}
	};

	/// @brief Saddle point (hyperbolic) flow field: v = (x, -y, 0)
	/// @details Classic saddle flow with stagnation point at origin.
	/// Properties: ∇·v = 0 (incompressible), ∇×v = 0 (irrotational)
	/// Eigenvalues: +1, -1, 0 (hyperbolic fixed point)
	class SaddleField : public IVectorFunction<3>
	{
	private:
		Real _strength;  ///< Scaling factor for the field
	public:
		/// @brief Default constructor with unit strength
		SaddleField() : _strength(1.0) {}
		/// @brief Constructor with specified strength
		SaddleField(Real strength) : _strength(strength) {}

		/// @brief Compute velocity at position x
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const
		{
			return VectorN<Real, 3>{_strength * x[0], -_strength * x[1], 0.0};
		}
	};

	/// @brief Source/sink field: v = Q·r/|r|³
	/// @details Radial flow from (source, Q>0) or toward (sink, Q<0) origin.
	/// Properties: ∇·v = 0 (except at origin), ∇×v = 0
	/// @note Q is the volumetric flow rate (source strength)
	class SourceSinkField : public IVectorFunction<3>
	{
	private:
		Real _strength;              ///< Q - source strength (negative for sink)
		VectorN<Real, 3> _center;    ///< Location of source/sink
	public:
		/// @brief Default constructor: unit source at origin
		SourceSinkField() : _strength(1.0), _center{0.0, 0.0, 0.0} {}
		/// @brief Constructor with specified strength at origin
		SourceSinkField(Real strength) : _strength(strength), _center{0.0, 0.0, 0.0} {}
		/// @brief Full constructor with arbitrary location
		SourceSinkField(Real strength, const VectorN<Real, 3>& center)
			: _strength(strength), _center(center) {}

		/// @brief Compute velocity at position x
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const
		{
			VectorN<Real, 3> r{x[0] - _center[0], x[1] - _center[1], x[2] - _center[2]};
			Real r_norm = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
			if (r_norm < 1e-15) return VectorN<Real, 3>{0.0, 0.0, 0.0};
			Real factor = _strength / (r_norm * r_norm * r_norm);
			return VectorN<Real, 3>{factor * r[0], factor * r[1], factor * r[2]};
		}
	};

	/// @brief Dipole field: superposition of source and sink
	/// @details Classic dipole with moment vector p. Far field: v ∝ (3(p·r̂)r̂ - p)/r³
	/// Used in electrostatics, magnetostatics, and potential flow.
	class DipoleField : public IVectorFunction<3>
	{
	private:
		VectorN<Real, 3> _moment;    ///< Dipole moment vector (direction and strength)
		VectorN<Real, 3> _center;    ///< Location of dipole
	public:
		/// @brief Default constructor: unit dipole along z at origin
		DipoleField() : _moment{0.0, 0.0, 1.0}, _center{0.0, 0.0, 0.0} {}
		/// @brief Constructor with specified moment at origin
		DipoleField(const VectorN<Real, 3>& moment) : _moment(moment), _center{0.0, 0.0, 0.0} {}
		/// @brief Full constructor with arbitrary location
		DipoleField(const VectorN<Real, 3>& moment, const VectorN<Real, 3>& center)
			: _moment(moment), _center(center) {}

		/// @brief Compute field at position x: (3(p·r̂)r̂ - p)/r³
		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const
		{
			VectorN<Real, 3> r{x[0] - _center[0], x[1] - _center[1], x[2] - _center[2]};
			Real r_sq = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
			if (r_sq < 1e-15) return VectorN<Real, 3>{0.0, 0.0, 0.0};
			
			Real r_norm = std::sqrt(r_sq);
			Real r5 = r_sq * r_sq * r_norm;  // r^5
			Real p_dot_r = _moment[0]*r[0] + _moment[1]*r[1] + _moment[2]*r[2];
			
			// v = (3(p·r)r - r²p) / r^5
			return VectorN<Real, 3>{
				static_cast<Real>((3.0 * p_dot_r * r[0] - r_sq * _moment[0]) / r5),
				static_cast<Real>((3.0 * p_dot_r * r[1] - r_sq * _moment[1]) / r5),
				static_cast<Real>((3.0 * p_dot_r * r[2] - r_sq * _moment[2]) / r5)
			};
		}
	};

	////////////////////             MATHEMATICAL SCALAR FIELDS          /////////////////

	/// @brief Gaussian bump scalar field: Φ = A·exp(-|r-c|²/σ²)
	/// @details Smooth, localized perturbation. Useful for testing, initial conditions.
	/// Properties: max at center, decays to zero, infinitely differentiable.
	class GaussianBumpField : public IScalarFunction<3>
	{
	private:
		Real _amplitude;            ///< Peak amplitude A
		Real _sigma;                ///< Width parameter σ
		VectorN<Real, 3> _center;   ///< Center of the bump
	public:
		/// @brief Default constructor: unit bump at origin with σ=1
		GaussianBumpField() : _amplitude(1.0), _sigma(1.0), _center{0.0, 0.0, 0.0} {}
		/// @brief Constructor with amplitude and width
		GaussianBumpField(Real amplitude, Real sigma)
			: _amplitude(amplitude), _sigma(sigma), _center{0.0, 0.0, 0.0} {}
		/// @brief Full constructor
		GaussianBumpField(Real amplitude, Real sigma, const VectorN<Real, 3>& center)
			: _amplitude(amplitude), _sigma(sigma), _center(center) {}

		/// @brief Compute scalar value at position x
		Real operator()(const VectorN<Real, 3>& x) const
		{
			Real dx = x[0] - _center[0];
			Real dy = x[1] - _center[1];
			Real dz = x[2] - _center[2];
			Real r_sq = dx*dx + dy*dy + dz*dz;
			return _amplitude * std::exp(-r_sq / (_sigma * _sigma));
		}
	};

	/// @brief Mexican hat (Ricker wavelet) scalar field: Φ = A(1 - r²/σ²)exp(-r²/2σ²)
	/// @details Laplacian of Gaussian. Zero mean, good for edge detection, wavelets.
	/// Properties: positive center, negative ring, zero at r = σ.
	class MexicanHatField : public IScalarFunction<3>
	{
	private:
		Real _amplitude;            ///< Peak amplitude A
		Real _sigma;                ///< Width parameter σ
		VectorN<Real, 3> _center;   ///< Center of the field
	public:
		/// @brief Default constructor
		MexicanHatField() : _amplitude(1.0), _sigma(1.0), _center{0.0, 0.0, 0.0} {}
		/// @brief Constructor with amplitude and width
		MexicanHatField(Real amplitude, Real sigma)
			: _amplitude(amplitude), _sigma(sigma), _center{0.0, 0.0, 0.0} {}
		/// @brief Full constructor
		MexicanHatField(Real amplitude, Real sigma, const VectorN<Real, 3>& center)
			: _amplitude(amplitude), _sigma(sigma), _center(center) {}

		/// @brief Compute scalar value at position x
		Real operator()(const VectorN<Real, 3>& x) const
		{
			Real dx = x[0] - _center[0];
			Real dy = x[1] - _center[1];
			Real dz = x[2] - _center[2];
			Real r_sq = dx*dx + dy*dy + dz*dz;
			Real sigma_sq = _sigma * _sigma;
			return _amplitude * (1.0 - r_sq / sigma_sq) * std::exp(-r_sq / (2.0 * sigma_sq));
		}
	};
}

#endif
