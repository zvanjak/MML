///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DAESystem.h                                                         ///
///  Description: Differential-algebraic equation system representations              ///
///               Initial value problem definitions for DAE solvers                   ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DAE_SYSTEM_H
#define MML_DAE_SYSTEM_H

#include <mml/MMLBase.h>
#include <algorithm>
#include <mml/MMLExceptions.h>

#include <mml/interfaces/IODESystemDAE.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/InterpolatedFunction.h>
#include <mml/base/TrajectoryStorage.h>

namespace MML
{
	///////////////////////////////////////////////////////////////////////////
	//                         DAE System Solution
	///////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Solution container for DAE system integration.
	 * 
	 * Stores time points, differential state vectors (x), and algebraic 
	 * state vectors (y) from DAE integration. Automatically grows storage 
	 * as needed using exponential growth strategy.
	 * 
	 * Provides interpolation methods for querying both differential and
	 * algebraic variables between saved points.
	 * 
	 * @note Unlike ODESystemSolution which stores only x(t), this class
	 *       stores both x(t) (differential) and y(t) (algebraic) trajectories.
	 */
	class DAESolution
	{
		int _numStepsOK, _numStepsBad;

		int _diffDim;    ///< Dimension of differential variables
		int _algDim;     ///< Dimension of algebraic variables
		int _totalSavedSteps;
		int _count;      ///< Actual number of data points stored (vs capacity)

		Real _t1, _t2;
		Vector<Real> _tval;    ///< Time values
		Matrix<Real> _xval;    ///< Differential variables (diffDim x numPoints)
		Matrix<Real> _yval;    ///< Algebraic variables (algDim x numPoints)

		/// @brief Grow storage (1.5x geometric growth) until it can hold requiredSize points.
		/// @note Loops rather than growing once (analysis v4, BASE D2).
		void EnsureCapacity(std::int64_t requiredSize)
		{
			if (requiredSize <= _totalSavedSteps)
				return;
			const int capacity = detail::TrajectoryCapacity(_totalSavedSteps, requiredSize, max_size(), true);
			detail::ResizeTrajectoryStorage(capacity, _tval, {&_xval, &_yval});
			_totalSavedSteps = capacity;
		}

		/// @brief Copy of the first `size()` columns of a (rows x capacity) storage matrix
		Matrix<Real> TrimmedColumns(const Matrix<Real>& storage, int rows) const
		{
			Matrix<Real> out(rows, _count);
			for (int r = 0; r < rows; r++)
				for (int i = 0; i < _count; i++)
					out[r][i] = storage[r][i];
			return out;
		}

	public:
		/// @brief Construct solution container with specified capacity
		/// @param t1 Initial time
		/// @param t2 Final time
		/// @param diffDim Dimension of differential variables
		/// @param algDim Dimension of algebraic variables
		/// @param maxSteps Initial storage capacity
		DAESolution(Real t1, Real t2, int diffDim, int algDim, int maxSteps)
			: _numStepsOK(0), _numStepsBad(0), _diffDim(diffDim), _algDim(algDim),
			  _totalSavedSteps(0), _count(0), _t1(t1), _t2(t2)
		{
			const int capacity = detail::TrajectoryCapacity(0, static_cast<std::int64_t>(maxSteps) + 1,
				max_size(), false);
			_xval.resize(diffDim, 0);
			_yval.resize(algDim, 0);
			detail::ResizeTrajectoryStorage(capacity, _tval, {&_xval, &_yval});
			_totalSavedSteps = capacity;
		}

		/// @brief Construct solution container with default capacity (1000 steps)
		DAESolution(Real t1, Real t2, int diffDim, int algDim)
			: DAESolution(t1, t2, diffDim, algDim, 1000) {}

		//=========================================================================
		//                           Statistics
		//=========================================================================

		/// @brief Increment count of successful integration steps
		void incrementSuccessfulSteps() { _numStepsOK++; }
		/// @brief Increment count of rejected integration steps
		void incrementRejectedSteps() { _numStepsBad++; }

		/// @brief Get count of successful integration steps
		int getNumStepsOK() const { return _numStepsOK; }
		/// @brief Get count of rejected integration steps
		int getNumStepsBad() const { return _numStepsBad; }
		/// @brief Get total integration steps attempted (OK + rejected)
		int getTotalNumSteps() const { return _numStepsOK + _numStepsBad; }

		//=========================================================================
		//                           Dimensions
		//=========================================================================

		/// @brief Get dimension of differential variables
		int getDiffDim() const { return _diffDim; }
		/// @brief Get dimension of algebraic variables
		int getAlgDim() const { return _algDim; }
		/// @brief Get total dimension (differential + algebraic)
		int getTotalDim() const { return _diffDim + _algDim; }

		/// @brief Get initial time
		Real getT1() const { return _t1; }
		/// @brief Get final time
		Real getT2() const { return _t2; }

		/// @brief Get current storage capacity
		int getTotalSavedSteps() const { return _totalSavedSteps; }
		/// @brief Get number of saved solution intervals
		int getNumSteps() const { return _count > 0 ? _count - 1 : 0; }
		/// @brief Check if solution is empty (no points saved)
		bool isEmpty() const { return _count == 0; }
		/// @brief Get number of saved solution points
		int size() const { return _count; }
		/// @brief Get current storage capacity
		int capacity() const { return _totalSavedSteps; }
		int max_size() const { return detail::TrajectoryMaxPoints({_diffDim, _algDim}); }

		/// Reserve a prospective point without changing the saved size or step counters.
		void reserveForPoint(int index) {
			if (index < 0) throw IndexError("Negative trajectory point index");
			EnsureCapacity(static_cast<std::int64_t>(index) + 1);
		}

		//=========================================================================
		//                           Data Access
		//=========================================================================

		/// @brief Get all time points as vector (only the `size()` stored points, not the capacity)
		Vector<Real> getTValues() const
		{
			Vector<Real> t(_count);
			for (int i = 0; i < _count; i++)
				t[i] = _tval[i];
			return t;
		}
		/// @brief Get all differential state vectors as matrix (diffDim x size())
		Matrix<Real> getXValues() const { return TrimmedColumns(_xval, _diffDim); }
		/// @brief Get all algebraic state vectors as matrix (algDim x size())
		Matrix<Real> getYValues() const { return TrimmedColumns(_yval, _algDim); }

		/// @brief Get t-value at a specific saved point
		/// @param ind Index of saved point
		/// @throws IndexError if index is invalid
		Real getTValue(int ind) const
		{
			if (ind < 0 || ind >= _count)
				throw IndexError("Index out of range in DAESolution::getTValue");
			return _tval[ind];
		}

		/// @brief Get differential variable value at a specific saved point
		/// @param ind Index of saved point
		/// @param component Differential component index (0 to diffDim-1)
		/// @throws IndexError if indices are invalid
		Real getXValue(int ind, int component) const
		{
			if (ind < 0 || ind >= _count)
				throw IndexError("Index out of range in DAESolution::getXValue(ind)");
			if (component < 0 || component >= _diffDim)
				throw IndexError("Component out of range in DAESolution::getXValue(component)");
			return _xval[component][ind];
		}

		/// @brief Get algebraic variable value at a specific saved point
		/// @param ind Index of saved point
		/// @param component Algebraic component index (0 to algDim-1)
		/// @throws IndexError if indices are invalid
		Real getYValue(int ind, int component) const
		{
			if (ind < 0 || ind >= _count)
				throw IndexError("Index out of range in DAESolution::getYValue(ind)");
			if (component < 0 || component >= _algDim)
				throw IndexError("Component out of range in DAESolution::getYValue(component)");
			return _yval[component][ind];
		}

		/// @brief Get all values for a specific differential component as vector
		/// @param component Differential component index (0 to diffDim-1)
		/// @throws IndexError if component index is invalid
		Vector<Real> getXValues(int component) const
		{
			if (component < 0 || component >= _diffDim)
				throw IndexError("Component out of range in DAESolution::getXValues(component)");
			Vector<Real> values(_count);
			for (int i = 0; i < _count; ++i) values[i] = _xval[component][i];
			return values;
		}

		/// @brief Get all values for a specific algebraic component as vector
		/// @param component Algebraic component index (0 to algDim-1)
		/// @throws IndexError if component index is invalid
		Vector<Real> getYValues(int component) const
		{
			if (component < 0 || component >= _algDim)
				throw IndexError("Component out of range in DAESolution::getYValues(component)");
			Vector<Real> values(_count);
			for (int i = 0; i < _count; ++i) values[i] = _yval[component][i];
			return values;
		}

		/// @brief Get differential state vector at final saved point
		/// @return Vector of all differential components at last time point
		Vector<Real> getXValuesAtEnd() const
		{
			if (_count == 0)
				throw IndexError("No data in DAESolution::getXValuesAtEnd");
			Vector<Real> res(_diffDim);
			for (int i = 0; i < _diffDim; i++)
				res[i] = _xval[i][_count - 1];
			return res;
		}

		/// @brief Get algebraic state vector at final saved point
		/// @return Vector of all algebraic components at last time point
		Vector<Real> getYValuesAtEnd() const
		{
			if (_count == 0)
				throw IndexError("No data in DAESolution::getYValuesAtEnd");
			Vector<Real> res(_algDim);
			for (int i = 0; i < _algDim; i++)
				res[i] = _yval[i][_count - 1];
			return res;
		}

		//=========================================================================
		//                           Data Storage
		//=========================================================================

		/// @brief Store time and state vectors at specified index
		/// @param ind Index where to store values
		/// @param t Time value
		/// @param x Differential state vector (must match diffDim)
		/// @param y Algebraic state vector (must match algDim)
		/// @throws IndexError if index is negative
		/// @throws ArgumentError if vector sizes don't match dimensions
		void fillValues(int ind, Real t, const Vector<Real>& x, const Vector<Real>& y)
		{
			if (ind < 0)
				throw IndexError("Index must be non-negative in DAESolution::fillValues");
			if (x.size() != _diffDim)
				throw ArgumentError("Differential vector size mismatch in DAESolution::fillValues");
			if (y.size() != _algDim)
				throw ArgumentError("Algebraic vector size mismatch in DAESolution::fillValues");

			reserveForPoint(ind);

			_tval[ind] = t;
			for (int i = 0; i < _diffDim; i++)
				_xval[i][ind] = x[i];
			for (int i = 0; i < _algDim; i++)
				_yval[i][ind] = y[i];
			if (ind + 1 > _count) _count = ind + 1;
		}

		/// @brief Set time value at specified index
		/// @param ind Index where to store time value
		/// @param t Time value to store
		void setTVal(int ind, Real t)
		{
			if (ind < 0)
				throw IndexError("Index must be non-negative in DAESolution::setTVal");
			reserveForPoint(ind);
			_tval[ind] = t;
			if (ind + 1 > _count) _count = ind + 1;
		}

		/// @brief Set differential variable value at specified index
		/// @param ind Index where to store value
		/// @param component Differential component index (0 to diffDim-1)
		/// @param val Value to store
		void setXVal(int ind, int component, Real val)
		{
			if (ind < 0)
				throw IndexError("Index must be non-negative in DAESolution::setXVal");
			if (component < 0 || component >= _diffDim)
				throw IndexError("Component index out of range in DAESolution::setXVal");
			reserveForPoint(ind);
			_xval[component][ind] = val;
			if (ind + 1 > _count) _count = ind + 1;
		}

		/// @brief Set algebraic variable value at specified index
		/// @param ind Index where to store value
		/// @param component Algebraic component index (0 to algDim-1)
		/// @param val Value to store
		void setYVal(int ind, int component, Real val)
		{
			if (ind < 0)
				throw IndexError("Index must be non-negative in DAESolution::setYVal");
			if (component < 0 || component >= _algDim)
				throw IndexError("Component index out of range in DAESolution::setYVal");
			reserveForPoint(ind);
			_yval[component][ind] = val;
			if (ind + 1 > _count) _count = ind + 1;
		}

		/// @brief Trim storage to actual number of saved steps
		/// @param numDoneSteps Number of steps actually completed
		/// @note Call this after integration completes to free unused memory
		void setFinalSize(int numDoneSteps)
		{
			const int capacity = detail::TrajectoryCapacity(_totalSavedSteps,
				static_cast<std::int64_t>(numDoneSteps) + 1, max_size(), false);
			detail::ResizeTrajectoryStorage(capacity, _tval, {&_xval, &_yval});
			_totalSavedSteps = _count = capacity;
		}

		//=========================================================================
		//                           Interpolation
		//=========================================================================

		/// @brief Create linear interpolator for specified differential component
		/// @param component Differential component index (0 to diffDim-1)
		/// @return Linear interpolation function
		LinearInterpRealFunc getXAsLinInterp(int component) const
		{
			Vector<Real> xsave = getTValues();
			Vector<Real> ysave = getXValues(component);
			return LinearInterpRealFunc(xsave, ysave);
		}

		/// @brief Create linear interpolator for specified algebraic component
		/// @param component Algebraic component index (0 to algDim-1)
		/// @return Linear interpolation function
		LinearInterpRealFunc getYAsLinInterp(int component) const
		{
			Vector<Real> xsave = getTValues();
			Vector<Real> ysave = getYValues(component);
			return LinearInterpRealFunc(xsave, ysave);
		}

		/// @brief Create cubic spline interpolator for specified differential component
		/// @param component Differential component index (0 to diffDim-1)
		/// @return Cubic spline interpolation function
		SplineInterpRealFunc getXAsSplineInterp(int component) const
		{
			Vector<Real> xsave = getTValues();
			Vector<Real> ysave = getXValues(component);
			return SplineInterpRealFunc(xsave, ysave);
		}

		/// @brief Create cubic spline interpolator for specified algebraic component
		/// @param component Algebraic component index (0 to algDim-1)
		/// @return Cubic spline interpolation function
		SplineInterpRealFunc getYAsSplineInterp(int component) const
		{
			Vector<Real> xsave = getTValues();
			Vector<Real> ysave = getYValues(component);
			return SplineInterpRealFunc(xsave, ysave);
		}

		/// @brief Create 2D parametric curve from two differential components
		/// @param ind1 First differential component index (x-coordinate)
		/// @param ind2 Second differential component index (y-coordinate)
		/// @return Spline-interpolated parametric curve in 2D
		SplineInterpParametricCurve<2> getXAsParamCurve2D(int ind1, int ind2) const
		{
			Matrix<Real> curve_points(_count, 2);
			for (int i = 0; i < _count; i++)
			{
				int source = _t1 > _t2 ? _count - 1 - i : i;
				curve_points(i, 0) = _xval[ind1][source];
				curve_points(i, 1) = _xval[ind2][source];
			}
			return SplineInterpParametricCurve<2>(std::min(_t1, _t2), std::max(_t1, _t2), curve_points);
		}

		/// @brief Create 3D parametric curve from three differential components
		/// @param ind1 First differential component index (x-coordinate)
		/// @param ind2 Second differential component index (y-coordinate)
		/// @param ind3 Third differential component index (z-coordinate)
		/// @return Spline-interpolated parametric curve in 3D
		SplineInterpParametricCurve<3> getXAsParamCurve3D(int ind1, int ind2, int ind3) const
		{
			Matrix<Real> curve_points(_count, 3);
			for (int i = 0; i < _count; i++)
			{
				int source = _t1 > _t2 ? _count - 1 - i : i;
				curve_points(i, 0) = _xval[ind1][source];
				curve_points(i, 1) = _xval[ind2][source];
				curve_points(i, 2) = _xval[ind3][source];
			}
			return SplineInterpParametricCurve<3>(std::min(_t1, _t2), std::max(_t1, _t2), curve_points);
		}
	};

	///////////////////////////////////////////////////////////////////////////
	//                         DAE System Classes
	///////////////////////////////////////////////////////////////////////////

	/**
	 * @brief Basic DAE system representation using function pointers.
	 * 
	 * Wraps C-style function pointers for use with DAE solvers.
	 * The functions compute differential equations and algebraic constraints.
	 */
	class DAESystem : public IODESystemDAE
	{
	protected:
		int _diffDim;
		int _algDim;

		/// Differential equations: f(t, x, y, dxdt)
		void (*_diffFunc)(Real, const Vector<Real>&, const Vector<Real>&, Vector<Real>&);
		
		/// Algebraic constraints: g(t, x, y, residual)
		void (*_algFunc)(Real, const Vector<Real>&, const Vector<Real>&, Vector<Real>&);

	public:
		/// @brief Default constructor (creates empty system)
		DAESystem()
			: _diffDim(0), _algDim(0), _diffFunc(nullptr), _algFunc(nullptr) {}

		/// @brief Construct DAE system with function pointers
		/// @param diffDim Dimension of differential variables
		/// @param algDim Dimension of algebraic variables
		/// @param diffFunc Function computing differential equations: f(t, x, y, dxdt)
		/// @param algFunc Function computing algebraic constraints: g(t, x, y, residual)
		DAESystem(int diffDim, int algDim,
		          void (*diffFunc)(Real, const Vector<Real>&, const Vector<Real>&, Vector<Real>&),
		          void (*algFunc)(Real, const Vector<Real>&, const Vector<Real>&, Vector<Real>&))
			: _diffDim(diffDim), _algDim(algDim), _diffFunc(diffFunc), _algFunc(algFunc)
		{
			if (diffFunc == nullptr || algFunc == nullptr)
				throw ArgumentError("DAESystem: function pointers cannot be null");
		}

		/// @brief Virtual destructor for proper cleanup in derived classes
		virtual ~DAESystem() = default;

		/// @brief Get dimension of differential variables
		int getDiffDim() const override { return _diffDim; }
		
		/// @brief Get dimension of algebraic variables
		int getAlgDim() const override { return _algDim; }

		/// @brief Compute differential equations dx/dt = f(t, x, y)
		/// @throws NotImplementedError if no differential function was provided
		void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>& y,
		             Vector<Real>& dxdt) const override
		{
			if (_diffFunc == nullptr)
				throw NotImplementedError("DAESystem::diffEqs - no differential function provided");
			_diffFunc(t, x, y, dxdt);
		}

		/// @brief Compute algebraic constraint residuals g(t, x, y)
		/// @throws NotImplementedError if no algebraic function was provided
		void algConstraints(Real t, const Vector<Real>& x, const Vector<Real>& y,
		                    Vector<Real>& g) const override
		{
			if (_algFunc == nullptr)
				throw NotImplementedError("DAESystem::algConstraints - no algebraic function provided");
			_algFunc(t, x, y, g);
		}
	};

	/**
	 * @brief DAE system with analytic Jacobians for implicit solvers.
	 * 
	 * Implements the Jacobian-capable interface through one unambiguous base chain.
	 * Reuses DAESystem internally for differential and algebraic callbacks.
	 */
	class DAESystemWithJacobian : public IODESystemDAEWithJacobian
	{
	private:
		DAESystem _system;
		/// Jacobian ∂f/∂x
		void (*_jacobian_fx)(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&);
		/// Jacobian ∂f/∂y
		void (*_jacobian_fy)(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&);
		/// Jacobian ∂g/∂x
		void (*_jacobian_gx)(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&);
		/// Jacobian ∂g/∂y
		void (*_jacobian_gy)(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&);

	public:
		/// @brief Default constructor (creates empty system)
		DAESystemWithJacobian()
			: _jacobian_fx(nullptr), _jacobian_fy(nullptr),
			  _jacobian_gx(nullptr), _jacobian_gy(nullptr) {}

		/// @brief Construct DAE system with all function pointers
		/// @param diffDim Dimension of differential variables
		/// @param algDim Dimension of algebraic variables
		/// @param diffFunc Function computing differential equations
		/// @param algFunc Function computing algebraic constraints
		/// @param jacobian_fx Function computing ∂f/∂x
		/// @param jacobian_fy Function computing ∂f/∂y
		/// @param jacobian_gx Function computing ∂g/∂x
		/// @param jacobian_gy Function computing ∂g/∂y
		/// @throws ArgumentError if any equation or Jacobian callback is null
		DAESystemWithJacobian(
			int diffDim, int algDim,
			void (*diffFunc)(Real, const Vector<Real>&, const Vector<Real>&, Vector<Real>&),
			void (*algFunc)(Real, const Vector<Real>&, const Vector<Real>&, Vector<Real>&),
			void (*jacobian_fx)(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&),
			void (*jacobian_fy)(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&),
			void (*jacobian_gx)(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&),
			void (*jacobian_gy)(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>&))
			: _system(diffDim, algDim, diffFunc, algFunc),
			  _jacobian_fx(jacobian_fx), _jacobian_fy(jacobian_fy),
			  _jacobian_gx(jacobian_gx), _jacobian_gy(jacobian_gy)
		{
			if (jacobian_fx == nullptr || jacobian_fy == nullptr || jacobian_gx == nullptr || jacobian_gy == nullptr)
				throw ArgumentError("DAESystemWithJacobian: all Jacobian function pointers are required");
		}

		int getDiffDim() const override { return _system.getDiffDim(); }
		int getAlgDim() const override { return _system.getAlgDim(); }

		void diffEqs(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& dxdt) const override
		{
			_system.diffEqs(t, x, y, dxdt);
		}

		void algConstraints(Real t, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& residual) const override
		{
			_system.algConstraints(t, x, y, residual);
		}

		/// @brief Compute Jacobian ∂f/∂x
		void jacobian_fx(Real t, const Vector<Real>& x, const Vector<Real>& y,
		                 Matrix<Real>& df_dx) const override
		{
			if (_jacobian_fx == nullptr)
				throw NotImplementedError("DAESystemWithJacobian::jacobian_fx - no Jacobian function provided");
			_jacobian_fx(t, x, y, df_dx);
		}

		/// @brief Compute Jacobian ∂f/∂y
		void jacobian_fy(Real t, const Vector<Real>& x, const Vector<Real>& y,
		                 Matrix<Real>& df_dy) const override
		{
			if (_jacobian_fy == nullptr)
				throw NotImplementedError("DAESystemWithJacobian::jacobian_fy - no Jacobian function provided");
			_jacobian_fy(t, x, y, df_dy);
		}

		/// @brief Compute Jacobian ∂g/∂x
		void jacobian_gx(Real t, const Vector<Real>& x, const Vector<Real>& y,
		                 Matrix<Real>& dg_dx) const override
		{
			if (_jacobian_gx == nullptr)
				throw NotImplementedError("DAESystemWithJacobian::jacobian_gx - no Jacobian function provided");
			_jacobian_gx(t, x, y, dg_dx);
		}

		/// @brief Compute Jacobian ∂g/∂y
		void jacobian_gy(Real t, const Vector<Real>& x, const Vector<Real>& y,
		                 Matrix<Real>& dg_dy) const override
		{
			if (_jacobian_gy == nullptr)
				throw NotImplementedError("DAESystemWithJacobian::jacobian_gy - no Jacobian function provided");
			_jacobian_gy(t, x, y, dg_dy);
		}
	};

} // namespace MML

#endif // MML_DAE_SYSTEM_H
