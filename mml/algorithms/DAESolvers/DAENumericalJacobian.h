///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DAENumericalJacobian.h                                              ///
///  Description: Adapter that provides numerical Jacobians for DAE systems           ///
///               Wraps any IODESystemDAE to satisfy IODESystemDAEWithJacobian        ///
///               using 4th-order central finite differences                          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DAE_NUMERICAL_JACOBIAN_H
#define MML_DAE_NUMERICAL_JACOBIAN_H

#include <mml/interfaces/IODESystemDAE.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>

#include <cmath>
#include <limits>

namespace MML {

	/// @brief Adapter that wraps any IODESystemDAE and provides numerical Jacobians.
	///
	/// Computes the four Jacobian matrices (∂f/∂x, ∂f/∂y, ∂g/∂x, ∂g/∂y)
	/// via 4th-order central finite differences, allowing DAE systems without
	/// analytical Jacobians to be used with all implicit DAE solvers.
	///
	/// Uses the stencil: f'(x) ≈ (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / (12h)
	///
	/// Each Jacobian column requires only 4 function evaluations regardless of
	/// the output dimension, making this efficient for moderate-dimensional systems.
	///
	/// @example
	/// @code
	/// // System without analytical Jacobians
	/// class MyDAE : public IODESystemDAE {
	///     int getDiffDim() const override { return 2; }
	///     int getAlgDim() const override { return 1; }
	///     void diffEqs(...) const override { /* ... */ }
	///     void algConstraints(...) const override { /* ... */ }
	/// };
	///
	/// MyDAE system;
	/// DAESystemNumericalJacobian wrapped(system);
	/// auto result = SolveDAEBackwardEuler(wrapped, t0, x0, y0, tEnd, config);
	/// @endcode
	class DAESystemNumericalJacobian : public IODESystemDAEWithJacobian
	{
		const IODESystemDAE& _system;
		Real _h;

	public:
		/// @brief Construct adapter wrapping a DAE system.
		/// @param system The DAE system to wrap (must outlive this adapter)
		/// @param h Step size for finite differences (0 = auto-select optimal h ≈ ε^(1/5))
		explicit DAESystemNumericalJacobian(const IODESystemDAE& system, Real h = 0.0)
			: _system(system), _h(h) {}

		//=========================================================================
		//                     Forwarded base interface
		//=========================================================================

		int getDiffDim() const override { return _system.getDiffDim(); }
		int getAlgDim() const override { return _system.getAlgDim(); }

		void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>& y,
		             Vector<Real>& dxdt) const override
		{
			_system.diffEqs(t, x, y, dxdt);
		}

		void algConstraints(Real t, const Vector<Real>& x, const Vector<Real>& y,
		                    Vector<Real>& g) const override
		{
			_system.algConstraints(t, x, y, g);
		}

		std::string getDiffVarName(int i) const override { return _system.getDiffVarName(i); }
		std::string getAlgVarName(int i) const override { return _system.getAlgVarName(i); }

		//=========================================================================
		//         Numerical Jacobians via 4th-order central differences
		//=========================================================================

		/// @brief Compute ∂f/∂x numerically (nDiff × nDiff matrix).
		void jacobian_fx(Real t, const Vector<Real>& x, const Vector<Real>& y,
		                 Matrix<Real>& df_dx) const override
		{
			const int nDiff = getDiffDim();
			const Real h = stepSize();
			Vector<Real> x_pert = x;
			Vector<Real> f_p2h(nDiff), f_ph(nDiff), f_mh(nDiff), f_m2h(nDiff);

			for (int j = 0; j < nDiff; ++j)
			{
				const Real x_orig = x_pert[j];

				x_pert[j] = x_orig + 2 * h;
				_system.diffEqs(t, x_pert, y, f_p2h);

				x_pert[j] = x_orig + h;
				_system.diffEqs(t, x_pert, y, f_ph);

				x_pert[j] = x_orig - h;
				_system.diffEqs(t, x_pert, y, f_mh);

				x_pert[j] = x_orig - 2 * h;
				_system.diffEqs(t, x_pert, y, f_m2h);

				x_pert[j] = x_orig;

				for (int i = 0; i < nDiff; ++i)
					df_dx(i, j) = (-f_p2h[i] + 8 * f_ph[i] - 8 * f_mh[i] + f_m2h[i]) / (12 * h);
			}
		}

		/// @brief Compute ∂f/∂y numerically (nDiff × nAlg matrix).
		void jacobian_fy(Real t, const Vector<Real>& x, const Vector<Real>& y,
		                 Matrix<Real>& df_dy) const override
		{
			const int nDiff = getDiffDim();
			const int nAlg = getAlgDim();
			const Real h = stepSize();
			Vector<Real> y_pert = y;
			Vector<Real> f_p2h(nDiff), f_ph(nDiff), f_mh(nDiff), f_m2h(nDiff);

			for (int j = 0; j < nAlg; ++j)
			{
				const Real y_orig = y_pert[j];

				y_pert[j] = y_orig + 2 * h;
				_system.diffEqs(t, x, y_pert, f_p2h);

				y_pert[j] = y_orig + h;
				_system.diffEqs(t, x, y_pert, f_ph);

				y_pert[j] = y_orig - h;
				_system.diffEqs(t, x, y_pert, f_mh);

				y_pert[j] = y_orig - 2 * h;
				_system.diffEqs(t, x, y_pert, f_m2h);

				y_pert[j] = y_orig;

				for (int i = 0; i < nDiff; ++i)
					df_dy(i, j) = (-f_p2h[i] + 8 * f_ph[i] - 8 * f_mh[i] + f_m2h[i]) / (12 * h);
			}
		}

		/// @brief Compute ∂g/∂x numerically (nAlg × nDiff matrix).
		void jacobian_gx(Real t, const Vector<Real>& x, const Vector<Real>& y,
		                 Matrix<Real>& dg_dx) const override
		{
			const int nDiff = getDiffDim();
			const int nAlg = getAlgDim();
			const Real h = stepSize();
			Vector<Real> x_pert = x;
			Vector<Real> g_p2h(nAlg), g_ph(nAlg), g_mh(nAlg), g_m2h(nAlg);

			for (int j = 0; j < nDiff; ++j)
			{
				const Real x_orig = x_pert[j];

				x_pert[j] = x_orig + 2 * h;
				_system.algConstraints(t, x_pert, y, g_p2h);

				x_pert[j] = x_orig + h;
				_system.algConstraints(t, x_pert, y, g_ph);

				x_pert[j] = x_orig - h;
				_system.algConstraints(t, x_pert, y, g_mh);

				x_pert[j] = x_orig - 2 * h;
				_system.algConstraints(t, x_pert, y, g_m2h);

				x_pert[j] = x_orig;

				for (int i = 0; i < nAlg; ++i)
					dg_dx(i, j) = (-g_p2h[i] + 8 * g_ph[i] - 8 * g_mh[i] + g_m2h[i]) / (12 * h);
			}
		}

		/// @brief Compute ∂g/∂y numerically (nAlg × nAlg matrix).
		void jacobian_gy(Real t, const Vector<Real>& x, const Vector<Real>& y,
		                 Matrix<Real>& dg_dy) const override
		{
			const int nAlg = getAlgDim();
			const Real h = stepSize();
			Vector<Real> y_pert = y;
			Vector<Real> g_p2h(nAlg), g_ph(nAlg), g_mh(nAlg), g_m2h(nAlg);

			for (int j = 0; j < nAlg; ++j)
			{
				const Real y_orig = y_pert[j];

				y_pert[j] = y_orig + 2 * h;
				_system.algConstraints(t, x, y_pert, g_p2h);

				y_pert[j] = y_orig + h;
				_system.algConstraints(t, x, y_pert, g_ph);

				y_pert[j] = y_orig - h;
				_system.algConstraints(t, x, y_pert, g_mh);

				y_pert[j] = y_orig - 2 * h;
				_system.algConstraints(t, x, y_pert, g_m2h);

				y_pert[j] = y_orig;

				for (int i = 0; i < nAlg; ++i)
					dg_dy(i, j) = (-g_p2h[i] + 8 * g_ph[i] - 8 * g_mh[i] + g_m2h[i]) / (12 * h);
			}
		}

		/// @brief Compute all four Jacobians efficiently.
		///
		/// Optimized: computes ∂f/∂x and ∂g/∂x together (shared x-perturbation),
		/// then ∂f/∂y and ∂g/∂y together (shared y-perturbation).
		void allJacobians(Real t, const Vector<Real>& x, const Vector<Real>& y,
		                  Matrix<Real>& df_dx, Matrix<Real>& df_dy,
		                  Matrix<Real>& dg_dx, Matrix<Real>& dg_dy) const override
		{
			const int nDiff = getDiffDim();
			const int nAlg = getAlgDim();
			const Real h = stepSize();

			// Compute ∂f/∂x and ∂g/∂x together (perturb x, evaluate both f and g)
			{
				Vector<Real> x_pert = x;
				Vector<Real> f_p2h(nDiff), f_ph(nDiff), f_mh(nDiff), f_m2h(nDiff);
				Vector<Real> g_p2h(nAlg), g_ph(nAlg), g_mh(nAlg), g_m2h(nAlg);

				for (int j = 0; j < nDiff; ++j)
				{
					const Real x_orig = x_pert[j];

					x_pert[j] = x_orig + 2 * h;
					_system.diffEqs(t, x_pert, y, f_p2h);
					_system.algConstraints(t, x_pert, y, g_p2h);

					x_pert[j] = x_orig + h;
					_system.diffEqs(t, x_pert, y, f_ph);
					_system.algConstraints(t, x_pert, y, g_ph);

					x_pert[j] = x_orig - h;
					_system.diffEqs(t, x_pert, y, f_mh);
					_system.algConstraints(t, x_pert, y, g_mh);

					x_pert[j] = x_orig - 2 * h;
					_system.diffEqs(t, x_pert, y, f_m2h);
					_system.algConstraints(t, x_pert, y, g_m2h);

					x_pert[j] = x_orig;

					for (int i = 0; i < nDiff; ++i)
						df_dx(i, j) = (-f_p2h[i] + 8 * f_ph[i] - 8 * f_mh[i] + f_m2h[i]) / (12 * h);

					for (int i = 0; i < nAlg; ++i)
						dg_dx(i, j) = (-g_p2h[i] + 8 * g_ph[i] - 8 * g_mh[i] + g_m2h[i]) / (12 * h);
				}
			}

			// Compute ∂f/∂y and ∂g/∂y together (perturb y, evaluate both f and g)
			{
				Vector<Real> y_pert = y;
				Vector<Real> f_p2h(nDiff), f_ph(nDiff), f_mh(nDiff), f_m2h(nDiff);
				Vector<Real> g_p2h(nAlg), g_ph(nAlg), g_mh(nAlg), g_m2h(nAlg);

				for (int j = 0; j < nAlg; ++j)
				{
					const Real y_orig = y_pert[j];

					y_pert[j] = y_orig + 2 * h;
					_system.diffEqs(t, x, y_pert, f_p2h);
					_system.algConstraints(t, x, y_pert, g_p2h);

					y_pert[j] = y_orig + h;
					_system.diffEqs(t, x, y_pert, f_ph);
					_system.algConstraints(t, x, y_pert, g_ph);

					y_pert[j] = y_orig - h;
					_system.diffEqs(t, x, y_pert, f_mh);
					_system.algConstraints(t, x, y_pert, g_mh);

					y_pert[j] = y_orig - 2 * h;
					_system.diffEqs(t, x, y_pert, f_m2h);
					_system.algConstraints(t, x, y_pert, g_m2h);

					y_pert[j] = y_orig;

					for (int i = 0; i < nDiff; ++i)
						df_dy(i, j) = (-f_p2h[i] + 8 * f_ph[i] - 8 * f_mh[i] + f_m2h[i]) / (12 * h);

					for (int i = 0; i < nAlg; ++i)
						dg_dy(i, j) = (-g_p2h[i] + 8 * g_ph[i] - 8 * g_mh[i] + g_m2h[i]) / (12 * h);
				}
			}
		}

	private:
		Real stepSize() const
		{
			if (_h > 0.0) return _h;
			return std::pow(std::numeric_limits<Real>::epsilon(), 0.2);
		}
	};

} // namespace MML

#endif // MML_DAE_NUMERICAL_JACOBIAN_H
