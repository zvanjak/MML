///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        PathIntegration.h                                                   ///
///  Description: Line/path integrals along parametric curves                         ///
///               Scalar and vector field integration along paths                     ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_PATH_INTEGRATION_H
#define MML_PATH_INTEGRATION_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/Vector.h"
#include "base/BaseUtils.h"

#include "core/Derivation.h"
#include "core/Integration.h"

#include "core/FieldOperations.h"

namespace MML
{
	class PathIntegration
	{
		template<int N>
		class HelperCurveLen : public IRealFunction
		{
			const IParametricCurve<N>& _curve;
		public:
			HelperCurveLen(const IParametricCurve<N>& curve) : _curve(curve) {}

			Real operator()(Real t) const
			{
				auto tangent_vec = Derivation::DeriveCurve<N>(_curve, t, nullptr);
				return tangent_vec.NormL2();
			}
		};

		template<int N>
		class HelperCurveMass : public IRealFunction
		{
			const IParametricCurve<N>& _curve;
			const IRealFunction& _density;	
		public:
			HelperCurveMass(const IParametricCurve<N>& curve, const IRealFunction &density) : _curve(curve), _density(density) {}

			Real operator()(Real t) const
			{
				auto tangent_vec = Derivation::DeriveCurve<N>(_curve, t, nullptr);
				return tangent_vec.NormL2() * _density(t);
			}
		};

		template<int N>
		class HelperLineIntegralScalarFunc : public IRealFunction
		{
			const IScalarFunction<N>& _scalar_field;
			const IParametricCurve<N>& _curve;
		public:
			HelperLineIntegralScalarFunc(const IScalarFunction<N>& scalarField, const IParametricCurve<N>& curve) 
					: _scalar_field(scalarField), _curve(curve) {}

			Real operator()(Real t) const
			{
				auto tangent_vec = Derivation::DeriveCurve<N>(_curve, t, nullptr);

				auto field_val = _scalar_field(_curve(t));
				auto ret = field_val * tangent_vec.NormL2();

				return ret;
			}
		};

		template<int N>
		class HelperLineIntegralVectorFunc : public IRealFunction
		{
			const IVectorFunction<N>& _vector_field;
			const IParametricCurve<N>& _curve;
		public:
			HelperLineIntegralVectorFunc(const IVectorFunction<N>& vectorField, const IParametricCurve<N>& curve) 
					: _vector_field(vectorField), _curve(curve) {}

			Real operator()(Real t) const
			{
				auto tangent_vec = Derivation::DeriveCurve<N>(_curve, t, nullptr);
				auto field_vec = _vector_field(_curve(t));

				return Utils::ScalarProduct<N>(tangent_vec, field_vec);
			}
		};

	public:
		template<int N>
		static Real ParametricCurveLength(const IParametricCurve<N>& curve, const Real a, const Real b)
		{
			HelperCurveLen helper(curve);

			return IntegrateTrap(helper, a, b);
		}
		
		template<int N>
		static Real ParametricCurveMass(const IParametricCurve<N>& curve, const IRealFunction &density, 
																		const Real a, const Real b)
		{
      HelperCurveMass<N> helper(curve, density);
      
      return IntegrateTrap(helper, a, b);
		}

		static Real LineIntegral(const IScalarFunction<3>& scalarField, const IParametricCurve<3>& curve, 
														 const Real t1, const Real t2, const Real eps = Defaults::WorkIntegralPrecision)
		{
			HelperLineIntegralScalarFunc helper(scalarField, curve);

			return IntegrateTrap(helper, t1, t2, nullptr, nullptr, eps);
		}

		static Real LineIntegral(const IVectorFunction<3>& vectorField, const IParametricCurve<3>& curve, 
														 const Real t1, const Real t2, const Real eps = Defaults::LineIntegralPrecision)
		{
			HelperLineIntegralVectorFunc helper(vectorField, curve);

			return IntegrateTrap(helper, t1, t2, nullptr, nullptr, eps);
		}
	};
} // end namespace

#endif // MML_PATH_INTEGRATION_H
