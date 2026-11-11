///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DerivationScalarFunction.h                                          ///
///  Description: Partial derivatives of scalar functions f:R^n->R                    ///
///               Gradient, Hessian, directional derivatives                          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DERIVATION_SCALAR_FUNCTION_H
#define MML_DERIVATION_SCALAR_FUNCTION_H

#include <mml/MMLBase.h>
#include <algorithm>

#include "DerivationBase.h"
#include "FirstDerivativeStencil.h"

#include <mml/base/Vector/VectorN.h>

namespace MML
{
	namespace Derivation
	{
		/********************************************************************************************************************/
		/********                               Numerical derivatives of FIRST order                                 ********/
		/********************************************************************************************************************/		/// @brief First-order partial derivative ∂f/∂xᵢ using forward difference (O(h) accuracy)
		/// @tparam N Dimension of domain R^N
		/// @param f Scalar function f:R^N→R
		/// @param deriv_index Index i of variable (0-based)
		/// @param point Evaluation point
		/// @param h Step size
		/// @param error Output: error estimate (optional)
		/// @return ∂f/∂xᵢ(point) ≈ [f(x+h·eᵢ) - f(x)] / h
		template <int N>
		static Real NDer1Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, 
														 Real h, Real* error = nullptr)
		{
			auto result = Detail::EvaluateScalarPartialFirstDerivativeStencil<Detail::FirstDerivativeOrder::One>(
				[&](int offset) { auto x = point; x[deriv_index] += offset * h; return f(x); }, h, error != nullptr);
			if (error) *error = result.error;
			return result.value;
		}
		template <int N>
		static Real NDer1Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, 
														 Real* error = nullptr)
		{
			return NDer1Partial(f, deriv_index, point, ScaleStep(NDer1_h, point[deriv_index]), error);
		}
		/// @brief Compute all N first-order partial derivatives (gradient components) at once
		/// @return Vector [∂f/∂x₁, ∂f/∂x₂, ..., ∂f/∂xₙ] (gradient ∇f)
		template <int N>
		static VectorN<Real, N> NDer1PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, 
																							Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer1Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer1Partial(f, i, point, h);
			}

			return ret;
		}
		template <int N>
		static VectorN<Real, N> NDer1PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, 
																							VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;
			for (int i = 0; i < N; i++) {
				if (error)
					ret[i] = NDer1Partial(f, i, point, &(*error)[i]);
				else
					ret[i] = NDer1Partial(f, i, point);
			}
			return ret;
		}

		/********************************************************************************************************************/
		/// @brief Second-order partial derivative ∂f/∂xᵢ using central difference (O(h²) accuracy)
		/// @note More accurate than 1st order: [f(x+h·eᵢ) - f(x-h·eᵢ)] / (2h)
		template <int N>
		static Real NDer2Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			auto result = Detail::EvaluateScalarPartialFirstDerivativeStencil<Detail::FirstDerivativeOrder::Two>(
				[&](int offset) { auto x = point; x[deriv_index] += offset * h; return f(x); }, h, error != nullptr);
			if (error) *error = result.error;
			return result.value;
		}
		template <int N>
		static Real NDer2Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, deriv_index, point, ScaleStep(NDer2_h, point[deriv_index]), error);
		}
		
		template <int N>
		static VectorN<Real, N> NDer2PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer2Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer2Partial(f, i, point, h);
			}

			return ret;
		}
		template <int N>
		static VectorN<Real, N> NDer2PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;
			for (int i = 0; i < N; i++) {
				if (error)
					ret[i] = NDer2Partial(f, i, point, &(*error)[i]);
				else
					ret[i] = NDer2Partial(f, i, point);
			}
			return ret;
		}
		
		/********************************************************************************************************************/
		/********                               Numerical derivatives of FOURTH order                                ********/
		/********************************************************************************************************************/
		/// @brief Fourth-order partial derivative ∂f/∂xᵢ using 5-point stencil (O(h⁴) accuracy)
		/// @note High accuracy: [f(x-2h) - 8f(x-h) + 8f(x+h) - f(x+2h)] / (12h)
		template <int N>
		static Real NDer4Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, 
														 Real h, Real* error = nullptr)
		{
			auto result = Detail::EvaluateScalarPartialFirstDerivativeStencil<Detail::FirstDerivativeOrder::Four>(
				[&](int offset) { auto x = point; x[deriv_index] += offset * h; return f(x); }, h, error != nullptr);
			if (error) *error = result.error;
			return result.value;
		}
		template <int N>
		static Real NDer4Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer4Partial(f, deriv_index, point, ScaleStep(NDer4_h, point[deriv_index]), error);
		}

		template <int N>
		static VectorN<Real, N> NDer4PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer4Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer4Partial(f, i, point, h);
			}

			return ret;
		}
		template <int N>
		static VectorN<Real, N> NDer4PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;
			for (int i = 0; i < N; i++) {
				if (error)
					ret[i] = NDer4Partial(f, i, point, &(*error)[i]);
				else
					ret[i] = NDer4Partial(f, i, point);
			}
			return ret;
		}
		
		/********************************************************************************************************************/
		/********                               Numerical derivatives of SIXTH order                                 ********/
		/********************************************************************************************************************/
		/// @brief Sixth-order partial derivative ∂f/∂xᵢ using 7-point stencil (O(h⁶) accuracy)
		/// @note Very high accuracy for smooth functions
		template <int N>
		static Real NDer6Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			auto result = Detail::EvaluateScalarPartialFirstDerivativeStencil<Detail::FirstDerivativeOrder::Six>(
				[&](int offset) { auto x = point; x[deriv_index] += offset * h; return f(x); }, h, error != nullptr);
			if (error) *error = result.error;
			return result.value;
		}
		template <int N>
		static Real NDer6Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer6Partial(f, deriv_index, point, ScaleStep(NDer6_h, point[deriv_index]), error);
		}

		template <int N>
		static VectorN<Real, N> NDer6PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer6Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer6Partial(f, i, point, h);
			}

			return ret;
		}
		template <int N>
		static VectorN<Real, N> NDer6PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;
			for (int i = 0; i < N; i++) {
				if (error)
					ret[i] = NDer6Partial(f, i, point, &(*error)[i]);
				else
					ret[i] = NDer6Partial(f, i, point);
			}
			return ret;
		}
		
		/********************************************************************************************************************/
		/********                               Numerical derivatives of EIGHTH order                                ********/
		/********************************************************************************************************************/
		/// @brief Eighth-order partial derivative ∂f/∂xᵢ using 9-point stencil (O(h⁸) accuracy)
		/// @note Highest accuracy, 9 function evaluations
		template <int N>
		static Real NDer8Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			auto result = Detail::EvaluateScalarPartialFirstDerivativeStencil<Detail::FirstDerivativeOrder::Eight>(
				[&](int offset) { auto x = point; x[deriv_index] += offset * h; return f(x); }, h, error != nullptr);
			if (error) *error = result.error;
			return result.value;
		}
		template <int N>
		static Real NDer8Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer8Partial(f, deriv_index, point, ScaleStep(NDer8_h, point[deriv_index]), error);
		}

		template <int N>
		static VectorN<Real, N> NDer8PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer8Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer8Partial(f, i, point, h);
			}

			return ret;
		}
		template <int N>
		static VectorN<Real, N> NDer8PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;
			for (int i = 0; i < N; i++) {
				if (error)
					ret[i] = NDer8Partial(f, i, point, &(*error)[i]);
				else
					ret[i] = NDer8Partial(f, i, point);
			}
			return ret;
		}
		
		/********************************************************************************************************************/
		/********                                      SECOND PARTIAL DERIVATIVES                                    ********/
		/********  NOTE: Direct finite difference formulas optimized for minimal function evaluations                ********/
		/********        For mixed partials (der_ind1 != der_ind2): ∂²f/∂x∂y                                        ********/
		/********        For pure partials (der_ind1 == der_ind2): ∂²f/∂x²                                          ********/
		/********                                                                                                    ********/
		/********        NSecDer2Partial: 3/4 evals (pure/mixed), 5/8 with error                                     ********/
		/********        NSecDer4Partial: 5/8 evals (pure/mixed), 7/12 with error                                    ********/
		/********************************************************************************************************************/
		
		// Second-order accurate (O(h²)) second partial derivative
		// For pure second partial (∂²f/∂xᵢ²): [f(x-h) - 2f(x) + f(x+h)] / h² - 3 function evaluations
		// For mixed partial (∂²f/∂xᵢ∂xⱼ): [f(x+h_i+h_j) - f(x+h_i-h_j) - f(x-h_i+h_j) + f(x-h_i-h_j)] / (4h²) - 4 evaluations
		template <int N>
		static Real NSecDer2Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, 
		                            const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			ValidateStepSize(h, "NSecDer2Partial");
			Real h2 = h * h;
			
			if (der_ind1 == der_ind2)
			{
				// Pure second partial ∂²f/∂xᵢ² - use standard 3-point formula
				auto x_plus = point;
				auto x_minus = point;
				
				x_plus[der_ind1] = point[der_ind1] + h;
				x_minus[der_ind1] = point[der_ind1] - h;
				
				Real f0 = f(point);
				Real fh = f(x_plus);
				Real fmh = f(x_minus);
				
				Real result = (fmh - MML::Real(2.0)*f0 + fh) / h2;
				
				if (error)
				{
					// Error estimate using 4th derivative
					auto x_2plus = point;
					auto x_2minus = point;
					x_2plus[der_ind1] = point[der_ind1] + 2*h;
					x_2minus[der_ind1] = point[der_ind1] - 2*h;
					
					const Real fm2h = f(x_2minus), f2h = f(x_2plus);
					*error = Detail::EstimateSecondTruncation2(fm2h, fmh, f0, fh, f2h, h,
						[](Real value) { return std::abs(value); }) +
					         Constants::Eps * (std::abs(fmh) + 2*std::abs(f0) + std::abs(fh)) / h2;
				}
				
				return result;
			}
			else
			{
				// Mixed partial ∂²f/∂xᵢ∂xⱼ - use 4-point cross-difference formula
				const Real halfI = error ? Detail::RefineMixedDerivativeStep(point[der_ind1], h) : Real(0);
				const Real halfJ = error ? Detail::RefineMixedDerivativeStep(point[der_ind2], h) : Real(0);
				auto x_pp = point;  // +h in both directions
				auto x_pm = point;  // +h in i, -h in j
				auto x_mp = point;  // -h in i, +h in j
				auto x_mm = point;  // -h in both directions
				
				x_pp[der_ind1] = point[der_ind1] + h;
				x_pp[der_ind2] = point[der_ind2] + h;
				
				x_pm[der_ind1] = point[der_ind1] + h;
				x_pm[der_ind2] = point[der_ind2] - h;
				
				x_mp[der_ind1] = point[der_ind1] - h;
				x_mp[der_ind2] = point[der_ind2] + h;
				
				x_mm[der_ind1] = point[der_ind1] - h;
				x_mm[der_ind2] = point[der_ind2] - h;
				
				Real fpp = f(x_pp);
				Real fpm = f(x_pm);
				Real fmp = f(x_mp);
				Real fmm = f(x_mm);
				
				Real result = (fpp - fpm - fmp + fmm) / (MML::Real(4.0) * h2);
				
				if (error)
				{
					auto fineSample = [&](Real signI, Real signJ) {
						auto x = point;
						x[der_ind1] += signI * halfI;
						x[der_ind2] += signJ * halfJ;
						return f(x);
					};
					const Real gpp = fineSample(1, 1), gpm = fineSample(1, -1);
					const Real gmp = fineSample(-1, 1), gmm = fineSample(-1, -1);
					const Real fine = (gpp - gpm - gmp + gmm) / h2;
					const Real coarseRoundoff = Constants::Eps *
						(std::abs(fpp) + std::abs(fpm) + std::abs(fmp) + std::abs(fmm)) / (4 * h2);
					const Real fineRoundoff = Constants::Eps *
						(std::abs(gpp) + std::abs(gpm) + std::abs(gmp) + std::abs(gmm)) / h2;
					*error = Detail::EstimateMixedDerivativeError<2>(result, fine, coarseRoundoff, fineRoundoff,
						[](Real value) { return std::abs(value); });
				}
				
				return result;
			}
		}
		template <int N>
		static Real NSecDer2Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, 
		                            const VectorN<Real, N>& point, Real* error = nullptr)
		{
			Real scale = std::max(std::abs(point[der_ind1]), std::abs(point[der_ind2]));
			return NSecDer2Partial(f, der_ind1, der_ind2, point, Precision::FourthRootEps * std::max(Real(1), scale), error);
		}

		/// @brief Fourth-order accurate (O(h⁴)) second partial derivative
		/// @note Pure: 5 evals (7 with error); mixed: 8 evals (12 with error).
		// Fourth-order accurate (O(h⁴)) second partial derivative
		// For pure second partial (∂²f/∂xᵢ²): [-f(x-2h) + 16f(x-h) - 30f(x) + 16f(x+h) - f(x+2h)] / (12h²) - 5 evaluations
		// For mixed partial (∂²f/∂xᵢ∂xⱼ): Richardson-extrapolated cross-differences - 8 evaluations
		template <int N>
		static Real NSecDer4Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, 
		                            const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			ValidateStepSize(h, "NSecDer4Partial");
			Real h2 = h * h;
			
			if (der_ind1 == der_ind2)
			{
				// Pure second partial ∂²f/∂xᵢ² - use 5-point formula
				auto x = point;
				Real x_orig = point[der_ind1];
				
				Real f0 = f(point);
				
				x[der_ind1] = x_orig + h;
				Real fh = f(x);
				
				x[der_ind1] = x_orig - h;
				Real fmh = f(x);
				
				x[der_ind1] = x_orig + 2*h;
				Real f2h = f(x);
				
				x[der_ind1] = x_orig - 2*h;
				Real fm2h = f(x);
				
				Real result = (-fm2h + MML::Real(16.0)*fmh - MML::Real(30.0)*f0 + MML::Real(16.0)*fh - f2h) / (MML::Real(12.0) * h2);
				
				if (error)
				{
					// Error estimate using 6th derivative
					x[der_ind1] = x_orig + 3*h;
					Real f3h = f(x);
					x[der_ind1] = x_orig - 3*h;
					Real fm3h = f(x);
					
					*error = Detail::EstimateSecondTruncation4(fm3h, fm2h, fmh, f0, fh, f2h, f3h, h,
						[](Real value) { return std::abs(value); }) +
					         Constants::Eps * (std::abs(fm2h) + 16*std::abs(fmh) + 30*std::abs(f0) + 
					                           16*std::abs(fh) + std::abs(f2h)) / (Real(12) * h2);
				}
				
				return result;
			}
			else
			{
				// Mixed partial ∂²f/∂xᵢ∂xⱼ using Richardson extrapolation for O(h⁴) accuracy
				// 
				// Standard 4-point cross-difference formula (O(h²)):
				//   D_h = [f(x+h,y+h) - f(x+h,y-h) - f(x-h,y+h) + f(x-h,y-h)] / (4h²)
				//
				// Richardson extrapolation with step h and 2h:
				//   D_h  = exact + c·h² + O(h⁴)
				//   D_2h = exact + 4c·h² + O(h⁴)
				//   (4·D_h - D_2h) / 3 = exact + O(h⁴)
				//
				// Simplifying: (16·d1 - d2) / (48·h²) where:
				//   d1 = f(±h,±h) cross-difference
				//   d2 = f(±2h,±2h) cross-difference
				
				Real x_i = point[der_ind1];
				Real x_j = point[der_ind2];
				const Real halfI = error ? Detail::RefineMixedDerivativeStep(x_i, h) : Real(0);
				const Real halfJ = error ? Detail::RefineMixedDerivativeStep(x_j, h) : Real(0);
				
				auto eval = [&](Real di, Real dj) {
					auto x = point;
					x[der_ind1] = x_i + di;
					x[der_ind2] = x_j + dj;
					return f(x);
				};
				
				// Cross-difference at step h
				Real f_p1_p1 = eval( h,  h);
				Real f_p1_m1 = eval( h, -h);
				Real f_m1_p1 = eval(-h,  h);
				Real f_m1_m1 = eval(-h, -h);
				Real d1 = (f_p1_p1 - f_p1_m1 - f_m1_p1 + f_m1_m1);
				
				// Cross-difference at step 2h  
				Real f_p2_p2 = eval( 2*h,  2*h);
				Real f_p2_m2 = eval( 2*h, -2*h);
				Real f_m2_p2 = eval(-2*h,  2*h);
				Real f_m2_m2 = eval(-2*h, -2*h);
				Real d2 = (f_p2_p2 - f_p2_m2 - f_m2_p2 + f_m2_m2);
				
				// Richardson extrapolation: (4·D_h - D_2h) / 3 = (16·d1 - d2) / (48·h²)
				Real result = (MML::Real(16.0) * d1 - d2) / (MML::Real(48.0) * h2);
				
				if (error)
				{
					const Real gpp = eval(halfI, halfJ), gpm = eval(halfI, -halfJ);
					const Real gmp = eval(-halfI, halfJ), gmm = eval(-halfI, -halfJ);
					const Real dHalf = gpp - gpm - gmp + gmm;
					const Real fine = (16 * dHalf - d1) / (12 * h2);
					const Real sum1 = std::abs(f_p1_p1) + std::abs(f_p1_m1) + std::abs(f_m1_p1) + std::abs(f_m1_m1);
					const Real sum2 = std::abs(f_p2_p2) + std::abs(f_p2_m2) + std::abs(f_m2_p2) + std::abs(f_m2_m2);
					const Real sumHalf = std::abs(gpp) + std::abs(gpm) + std::abs(gmp) + std::abs(gmm);
					const Real coarseRoundoff = Constants::Eps * (16 * sum1 + sum2) / (48 * h2);
					const Real fineRoundoff = Constants::Eps * (16 * sumHalf + sum1) / (12 * h2);
					*error = Detail::EstimateMixedDerivativeError<4>(result, fine, coarseRoundoff, fineRoundoff,
						[](Real value) { return std::abs(value); });
				}
				
				return result;
			}
		}
		template <int N>
		static Real NSecDer4Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, 
		                            const VectorN<Real, N>& point, Real* error = nullptr)
		{
			Real scale = std::max(std::abs(point[der_ind1]), std::abs(point[der_ind2]));
			return NSecDer4Partial(f, der_ind1, der_ind2, point, Precision::SixthRootEps * std::max(Real(1), scale), error);
		}

		/********************************************************************************************************************/
		/********                            Definitions of default derivation functions                             ********/
		/********************************************************************************************************************/
		template<int N>
		static inline Real(*DerivePartial)(const IScalarFunction<N>& f, int deriv_index, 
																			 const VectorN<Real, N>& point, Real* error) = Derivation::NDer4Partial;
		template<int N>
		static inline Real(*DeriveSecPartial)(const IScalarFunction<N>& f, int der_ind1, int der_ind2, 
																					const VectorN<Real, N>& point, Real* error) = Derivation::NSecDer4Partial;
		template<int N>
		static inline VectorN<Real, N>(*DerivePartialAll)(const IScalarFunction<N>& f, const VectorN<Real, N>& point, 
																											VectorN<Real, N>* error) = Derivation::NDer4PartialByAll;
	}
}

#endif // MML_DERIVATION_SCALAR_FUNCTION_H
