///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        PathIntegration.h                                                   ///
///  Description: Line/path integrals along parametric curves                         ///
///               Scalar and vector field integration along paths                     ///
///                                                                                   ///
///  Features:    - Parametric curve arc length computation                           ///
///               - Curve mass with density function                                  ///
///               - Line integrals of scalar fields: ∫_C f(r) ds                     ///
///               - Line integrals of vector fields: ∫_C F·dr (work integrals)       ///
///                                                                                   ///
///  Usage:                                                                           ///
///    // Arc length of helix from t=0 to t=2π                                        ///
///    Real len = PathIntegration::ParametricCurveLength<3>(helix, 0, 2*PI);          ///
///                                                                                   ///
///    // Work done by force field F along path C                                     ///
///    Real work = PathIntegration::LineIntegral(F, curve, t1, t2);                   ///
///                                                                                   ///
///  Mathematical Background:                                                         ///
///    Arc length:     L = ∫_a^b |r'(t)| dt                                          ///
///    Scalar line:    ∫_C f ds = ∫_a^b f(r(t)) |r'(t)| dt                           ///
///    Vector line:    ∫_C F·dr = ∫_a^b F(r(t))·r'(t) dt                             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#ifndef MML_PATH_INTEGRATION_H
#define MML_PATH_INTEGRATION_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/Vector/Vector.h"
#include "base/BaseUtils.h"

#include "core/Derivation.h"
#include "core/Integration.h"

#include "core/FieldOperations.h"

namespace MML
{
	///////////////////////////////////////////////////////////////////////////
	///                         PathIntegration                             ///
	///////////////////////////////////////////////////////////////////////////
	/// @brief Static class for computing line/path integrals along curves
	/// 
	/// Provides methods for:
	/// - Computing arc length of parametric curves
	/// - Computing mass along curves with varying density
	/// - Line integrals of scalar fields (∫_C f ds)
	/// - Line integrals of vector fields (∫_C F·dr) - work integrals
	///
	/// @note All methods use trapezoidal integration internally
	/// @see IParametricCurve, IScalarFunction, IVectorFunction
	///////////////////////////////////////////////////////////////////////////
	class PathIntegration
	{
	private:
		//////////////////////////////////////////////////////////////////
		// Helper functors for converting path integrals to 1D integrals
		//////////////////////////////////////////////////////////////////

		/// @brief Helper for arc length: |r'(t)|
		/// @details Computes the speed (magnitude of tangent) at parameter t
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

		/// @brief Helper for curve mass: ρ(t) * |r'(t)|
		/// @details Integrand for mass computation with density function
		template<int N>
		class HelperCurveMass : public IRealFunction
		{
			const IParametricCurve<N>& _curve;
			const IRealFunction& _density;	
		public:
			HelperCurveMass(const IParametricCurve<N>& curve, const IRealFunction &density) 
				: _curve(curve), _density(density) {}

			Real operator()(Real t) const
			{
				auto tangent_vec = Derivation::DeriveCurve<N>(_curve, t, nullptr);
				return tangent_vec.NormL2() * _density(t);
			}
		};

		/// @brief Helper for scalar field line integral: f(r(t)) * |r'(t)|
		/// @details Converts ∫_C f ds to ∫_a^b f(r(t))|r'(t)| dt
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
				return field_val * tangent_vec.NormL2();
			}
		};

		/// @brief Helper for vector field line integral: F(r(t)) · r'(t)
		/// @details Converts ∫_C F·dr to ∫_a^b F(r(t))·r'(t) dt
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
		///////////////////////////////////////////////////////////////////////
		///                     Arc Length Computation                      ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Compute the arc length of a parametric curve
		/// 
		/// Computes L = ∫_a^b |r'(t)| dt using trapezoidal integration.
		///
		/// @tparam N Dimension of the curve (2 or 3 typically)
		/// @param curve The parametric curve r(t)
		/// @param a Starting parameter value
		/// @param b Ending parameter value
		/// @return IntegrationResult with arc length as value (implicitly converts to Real)
		///
		/// @par Example:
		/// @code
		/// // Length of a helix: r(t) = (cos(t), sin(t), t/2π), t ∈ [0, 2π]
		/// ParametricCurve<3> helix([](Real t) {
		///     return VectorN<3>({cos(t), sin(t), t/(2*PI)});
		/// });
		/// Real length = PathIntegration::ParametricCurveLength<3>(helix, 0, 2*PI);
		/// // length ≈ 2π√(1 + 1/(4π²)) for one turn
		/// @endcode
		template<int N>
		static IntegrationResult ParametricCurveLength(const IParametricCurve<N>& curve, const Real a, const Real b)
		{
			HelperCurveLen<N> helper(curve);
			return IntegrateTrap(helper, a, b);
		}
		
		///////////////////////////////////////////////////////////////////////
		///                     Curve Mass Computation                      ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Compute mass of a wire with variable density
		/// 
		/// Computes M = ∫_a^b ρ(t) |r'(t)| dt where ρ is the linear density.
		///
		/// @tparam N Dimension of the curve
		/// @param curve The parametric curve r(t) representing the wire
		/// @param density Linear density function ρ(t) [mass per unit length]
		/// @param a Starting parameter value
		/// @param b Ending parameter value
		/// @return IntegrationResult with total mass as value (implicitly converts to Real)
		///
		/// @par Physical interpretation:
		/// For a wire bent along curve C with density ρ, the mass element
		/// is dm = ρ ds where ds = |r'(t)| dt is the arc length element.
		///
		/// @par Example:
		/// @code
		/// // Wire along unit circle with density ρ(t) = 1 + sin(t)
		/// ParametricCurve<2> circle([](Real t) {
		///     return VectorN<2>({cos(t), sin(t)});
		/// });
		/// RealFunction density([](Real t) { return 1 + sin(t); });
		/// Real mass = PathIntegration::ParametricCurveMass<2>(circle, density, 0, 2*PI);
		/// @endcode
		template<int N>
		static IntegrationResult ParametricCurveMass(const IParametricCurve<N>& curve, const IRealFunction &density, 
															const Real a, const Real b)
		{
      HelperCurveMass<N> helper(curve, density);
      return IntegrateTrap(helper, a, b);
		}

		///////////////////////////////////////////////////////////////////////
		///                   Scalar Field Line Integral                    ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Line integral of a scalar field over a curve
		/// 
		/// Computes ∫_C f ds = ∫_{t1}^{t2} f(r(t)) |r'(t)| dt
		///
		/// @param scalarField Scalar field f: ℝ³ → ℝ to integrate
		/// @param curve Parametric curve C: [t1,t2] → ℝ³
		/// @param t1 Starting parameter value
		/// @param t2 Ending parameter value
		/// @param eps Desired relative precision (default: WorkIntegralPrecision)
		/// @return IntegrationResult with line integral value (implicitly converts to Real)
		///
		/// @par Physical interpretation:
		/// If f represents a scalar quantity (temperature, density, etc.),
		/// this integral gives the "accumulated" value along the path,
		/// weighted by path length.
		///
		/// @par Example:
		/// @code
		/// // Integrate temperature field T(x,y,z) = x² + y² along helix
		/// ScalarFunction<3> temperature([](const VectorN<3>& p) {
		///     return p[0]*p[0] + p[1]*p[1];
		/// });
		/// Real total = PathIntegration::LineIntegral(temperature, helix, 0, 2*PI);
		/// @endcode
		static IntegrationResult LineIntegral(const IScalarFunction<3>& scalarField, const IParametricCurve<3>& curve, 
														 const Real t1, const Real t2, const Real eps = Defaults::WorkIntegralPrecision)
		{
			HelperLineIntegralScalarFunc<3> helper(scalarField, curve);
			return IntegrateTrap(helper, t1, t2, eps);
		}

		///////////////////////////////////////////////////////////////////////
		///                   Vector Field Line Integral                    ///
		///////////////////////////////////////////////////////////////////////

		/// @brief Line integral of a vector field over a curve (work integral)
		/// 
		/// Computes ∫_C F·dr = ∫_{t1}^{t2} F(r(t))·r'(t) dt
		///
		/// @param vectorField Vector field F: ℝ³ → ℝ³ to integrate
		/// @param curve Parametric curve C: [t1,t2] → ℝ³
		/// @param t1 Starting parameter value
		/// @param t2 Ending parameter value
		/// @param eps Desired relative precision (default: LineIntegralPrecision)
		/// @return IntegrationResult with work integral value (implicitly converts to Real)
		///
		/// @par Physical interpretation:
		/// This computes the work done by force field F moving a particle
		/// along curve C. Also known as circulation when curve is closed.
		///
		/// @par Conservative fields:
		/// If F = ∇φ (gradient of potential), then ∫_C F·dr = φ(B) - φ(A)
		/// (path-independent, depends only on endpoints).
		///
		/// @par Example:
		/// @code
		/// // Work done by gravitational field F = -mg k̂ along parabolic path
		/// VectorFunction<3> gravity([](const VectorN<3>& p) {
		///     return VectorN<3>({0, 0, -9.81});  // F = -g in z-direction
		/// });
		/// // Path: projectile motion r(t) = (v₀t, 0, -½gt²)
		/// ParametricCurve<3> path([](Real t) {
		///     return VectorN<3>({10*t, 0, -0.5*9.81*t*t});
		/// });
		/// Real work = PathIntegration::LineIntegral(gravity, path, 0, 2);
		/// @endcode
		///
		/// @note For circulation of vector fields, ensure curve is closed:
		///       r(t1) = r(t2). For conservative fields, result will be ≈ 0.
		static IntegrationResult LineIntegral(const IVectorFunction<3>& vectorField, const IParametricCurve<3>& curve, 
														 const Real t1, const Real t2, const Real eps = Defaults::LineIntegralPrecision)
		{
			HelperLineIntegralVectorFunc<3> helper(vectorField, curve);
			return IntegrateTrap(helper, t1, t2, eps);
		}
	};

} // end namespace MML

#endif // MML_PATH_INTEGRATION_H
