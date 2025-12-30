///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        Fields.h                                                            ///
///  Description: Scalar, vector, and tensor field classes                            ///
///               Field representations for physical simulations                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FIELDS_H
#define MML_FIELDS_H


#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "interfaces/IFunction.h"

#include "base/VectorN.h"
#include "base/VectorTypes.h"

#include "core/FieldOperations.h"

#endif


namespace MML::Fields
{
	////////////////////             INVERSE RADIAL FIELD                /////////////////
	// Potential fields
	static Real InverseRadialPotentialFieldCart(const VectorN<Real, 3>& x) { return 1.0 / x.NormL2(); }
	static Real InverseRadialPotentialFieldCart(Real constant, const VectorN<Real, 3>& x) { return constant / x.NormL2(); }
	static Real InverseRadialPotentialFieldSpher(const VectorN<Real, 3>& x) { return 1.0 / x[0]; }
	static Real InverseRadialPotentialFieldSpher(Real constant, const VectorN<Real, 3>& x) { return constant / x[0]; }
	static Real InverseRadialPotentialFieldCyl(const VectorN<Real, 3>& x) { return 1.0 / sqrt(x[0] * x[0] + x[2] * x[2]); }
	static Real InverseRadialPotentialFieldCyl(Real constant, const VectorN<Real, 3>& x) { return constant / sqrt(x[0] * x[0] + x[2] * x[2]); }

	// Force fields
	// NOTE: Cannot change to Vector3Cartesian - VectorFunction<3> interface requires VectorN<Real,3> function pointers
	static VectorN<Real, 3> InverseRadialPotentialForceFieldCart(const VectorN<Real, 3>& x)
	{
		return (-1) * x / std::pow((x).NormL2(), 3);
	}
	static VectorN<Real, 3> InverseRadialPotentialForceFieldCart(Real constant, const VectorN<Real, 3>& x)
	{
		return -constant * x / std::pow((x).NormL2(), 3);
	}
	// can't be Vector3Spherical, because in that case we can't form a ScalarFunction<3> out of it :()
	static VectorN<Real, 3> InverseRadialPotentialForceFieldSph(const VectorN<Real, 3>& x)
	{
		return VectorN<Real, 3>{-1 / (x[0] * x[0]), 0.0, 0.0};
	}
	static VectorN<Real, 3> InverseRadialPotentialForceFieldSph(Real constant, const VectorN<Real, 3>& x)
	{
		return VectorN<Real, 3>{-constant / (x[0] * x[0]), 0.0, 0.0};
	}

	class InverseRadialFieldCart : public IScalarFunction<3>
	{
	protected:
		Real _constant;
	public:
		InverseRadialFieldCart() : _constant(-1.0) {}
		InverseRadialFieldCart(Real constant) : _constant(constant) {}

		Real operator()(const VectorN<Real, 3>& x) const { return _constant * Fields::InverseRadialPotentialFieldCart(x); }
	};
	class InverseRadialFieldSpher : public IScalarFunction<3>
	{
	protected:
		Real _constant;
	public:
		InverseRadialFieldSpher() : _constant(-1.0) {}
		InverseRadialFieldSpher(Real constant) : _constant(constant) {}

		Real operator()(const VectorN<Real, 3>& x) const { return _constant * InverseRadialPotentialFieldSpher(x); }
	};

	class InverseRadialForceFieldCart : public IVectorFunction<3>
	{
	private:
		Real _constant;
	public:
		InverseRadialForceFieldCart() : _constant(-1.0) {}
		InverseRadialForceFieldCart(Real constant) : _constant(constant) {}

		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const { return _constant * Fields::InverseRadialPotentialForceFieldCart(x); }
	};
	class InverseRadialForceFieldSpher : public IVectorFunction<3>
	{
	private:
		Real _constant;
	public:
		InverseRadialForceFieldSpher() : _constant(-1.0) {}
		InverseRadialForceFieldSpher(Real constant) : _constant(constant) {}

		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const { return _constant * InverseRadialPotentialForceFieldSph(x); }
	};
	// TODO - dodati multibody gravity
	// TODO - dodati InfiniteLineCurrentField
	// TODO - dodati InfinitePlaneCurrentField

}

#endif
