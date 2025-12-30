///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ODESystemSolution.h                                                 ///
///  Description: Solution container for ODE integration results                      ///
///               Storage and interpolation of computed trajectories                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_ODE_SYSTEM_SOLUTON_H
#define MML_ODE_SYSTEM_SOLUTON_H

#include "MMLBase.h"
#include "MMLExceptions.h"

#include "base/InterpolatedFunction.h"
//#include "base/ODESystem.h"


namespace MML
{
	/// @brief Solution container for ODE system integration
	/// @details Stores time points and state vectors from ODE integration.
	///          Automatically grows storage as needed using exponential growth strategy.
	///          Provides interpolation methods for querying solution between saved points.
	class ODESystemSolution
	{
		int  _numStepsOK, _numStepsBad;

		int  _sysDim;
		int  _totalSavedSteps;
	
		Real _t1, _t2;
		Vector<Real> _tval;
		Matrix<Real> _xval;
	
		/// @brief Extend storage capacity using exponential growth (1.5x)
		void ExtendSavedSteps()
		{
			_totalSavedSteps = static_cast<int>(_totalSavedSteps * 1.5) + 1;
			
			// resizing the storage while preserving old values
			_tval.Resize(_totalSavedSteps, true);
			_xval.Resize(_sysDim, _totalSavedSteps, true);
		}

	public:
		/// @brief Construct solution container with specified capacity
		/// @param x1 Initial time
		/// @param x2 Final time
		/// @param dim System dimension (number of equations)
		/// @param maxSteps Initial storage capacity
		ODESystemSolution(Real x1, Real x2, int dim, int maxSteps)
			: _t1(x1), _t2(x2), _sysDim(dim),
				_numStepsOK(0), _numStepsBad(0), _totalSavedSteps(maxSteps + 1)
		{
			_tval.Resize(_totalSavedSteps);
			_xval.Resize(dim, _totalSavedSteps);
		}
		
		/// @brief Construct solution container with default capacity (1000 steps)
		ODESystemSolution(Real x1, Real x2, int dim) : ODESystemSolution(x1, x2, dim, 1000)	{ }

		/// @brief Increment count of successful integration steps
		void incrementSuccessfulSteps()	{ _numStepsOK++; }
		/// @brief Increment count of rejected integration steps
		void incrementRejectedSteps() { _numStepsBad++; }

		/// @brief Get system dimension (number of equations)
		int getSysDim() const { return _sysDim; }
		
		/// @brief Get initial time
		Real getT1() const { return _t1; }
		/// @brief Get final time
		Real getT2() const { return _t2; }

		/// @brief Get count of successful integration steps
		int getNumStepsOK() const { return _numStepsOK; }
		/// @brief Get count of rejected integration steps
		int getNumStepsBad() const { return _numStepsBad; }
		/// @brief Get total integration steps attempted (OK + rejected)
		int getTotalNumSteps() const { return _numStepsOK + _numStepsBad; }
		/// @brief Get current storage capacity
		int getTotalSavedSteps() const { return _totalSavedSteps; }
		/// @brief Get number of saved solution intervals
		int getNumSteps() const { return _totalSavedSteps - 1; }  // Number of intervals
		
		/// @brief Check if solution is empty (no points saved)
		bool isEmpty() const { return _totalSavedSteps == 0; }
		/// @brief Get number of saved solution points
		int size() const { return _totalSavedSteps; }
		/// @brief Get current storage capacity
		int capacity() const { return _totalSavedSteps; }
		
		/// @brief Pre-allocate storage for specified number of points
		/// @param n Number of points to reserve space for
		void reserve(int n) {
			if (n > _totalSavedSteps) {
				_totalSavedSteps = n;
				_tval.Resize(_totalSavedSteps, true);
				_xval.Resize(_sysDim, _totalSavedSteps, true);
			}
		}

		/// @brief Get all time points as vector
		Vector<Real> getTValues() const { return _tval; }
		/// @brief Get all state vectors as matrix (dim x numPoints)
		Matrix<Real> getXValues() const { return _xval; }
		
		/// @brief Get t-value at a specific saved point
		/// @param ind Index of saved point
		/// @throws std::out_of_range if index is invalid
		Real getTValue(int ind) const { 
			if (ind < 0 || ind >= _totalSavedSteps)
				throw std::out_of_range("Index out of range in ODESystemSolution::getTValue");
			return _tval[ind]; 
		}
		
		/// @brief Get x-value for a component at a specific saved point
		/// @param ind Index of saved point
		/// @param component Component index (0 to dim-1)
		/// @throws std::out_of_range if indices are invalid
		Real getXValue(int ind, int component) const { 
			if (ind < 0 || ind >= _totalSavedSteps)
				throw std::out_of_range("Index out of range in ODESystemSolution::getXValue(ind)");
			if (component < 0 || component >= _sysDim)
				throw std::out_of_range("Component out of range in ODESystemSolution::getXValue(component)");
			return _xval[component][ind]; 
		}
		
		/// @brief Get all values for a specific component as vector
		/// @param component Component index (0 to dim-1)
		/// @throws std::out_of_range if component index is invalid
		Vector<Real> getXValues(int component) const 
		{ 
			if (component < 0 || component >= _sysDim)
				throw std::out_of_range("Component out of range in ODESystemSolution::getXValues(component)");
			
			return _xval.VectorFromRow(component); 
		}

		/// @brief Get state vector at final saved point
		/// @return Vector of all components at last time point
		Vector<Real> getXValuesAtEnd() const
		{
			Vector<Real> res(_sysDim);
			for (int i = 0; i < _sysDim; i++)
				res[i] = _xval[i][_totalSavedSteps - 1];
			return res;
		}

		/// @brief Set time value at specified index
		/// @param ind Index where to store time value
		/// @param x Time value to store
		/// @throws std::out_of_range if index is negative
		void setTVal(int ind, Real x)
		{
			if (ind < 0 )
				throw std::out_of_range("Index must be non-negative in ODESystemSolution::setTVal");
			
			if (ind >= _totalSavedSteps)
				ExtendSavedSteps();

			_tval[ind] = x;
		}
		/// @brief Set state value for specific component at specified index
		/// @param ind Index where to store value
		/// @param component Component index (0 to dim-1)
		/// @param x Value to store
		/// @throws std::out_of_range if indices are invalid
		void setXVal(int ind, int component, Real x)
		{
			if (ind < 0)
				throw std::out_of_range("Index must be non-negative in ODESystemSolution::setXVal(ind)");

			if (component < 0 || component >= _sysDim)
				throw std::out_of_range("Component index out of range in ODESystemSolution::setXVal(component)");

			if (ind >= _totalSavedSteps)
				ExtendSavedSteps();

			_xval[component][ind] = x;
		}

		/// @brief Store time and state vector at specified index
		/// @param ind Index where to store values
		/// @param x Time value
		/// @param y State vector (must match system dimension)
		/// @throws std::out_of_range if index is negative
		/// @throws std::invalid_argument if vector size doesn't match system dimension
		void fillValues(int ind, Real x, const Vector<Real>& y)
		{
			if (ind < 0 )
				throw std::out_of_range("Index must be non-negative in ODESystemSolution::fillValues");

			if( y.size() != _sysDim )
				throw std::invalid_argument("Vector size mismatch in ODESystemSolution::fillValues");

			if(ind >= _totalSavedSteps)
				ExtendSavedSteps();

			_tval[ind] = x;
			for (int i = 0; i < _sysDim; i++)
				_xval[i][ind] = y[i];
		}
		/// @brief Trim storage to actual number of saved steps
		/// @param numDoneSteps Number of steps actually completed
		/// @note Call this after integration completes to free unused memory
		void setFinalSize(int numDoneSteps)
		{
			_totalSavedSteps = numDoneSteps + 1;

			_tval.Resize(_totalSavedSteps, true);
			_xval.Resize(_sysDim, _totalSavedSteps, true);
		}

		/// @brief Create linear interpolator for specified component
		/// @param component Component index (0 to dim-1)
		/// @return Linear interpolation function
		LinearInterpRealFunc  getSolAsLinInterp(int component) const
		{
			Vector<Real> xsave = _tval;
			Vector<Real> ysave = _xval.VectorFromRow(component);

			return LinearInterpRealFunc(xsave, ysave);
		}
		
		/// @brief Create polynomial interpolator for specified component
		/// @param component Component index (0 to dim-1)
		/// @param polyOrder Polynomial order for interpolation
		/// @return Polynomial interpolation function
		PolynomInterpRealFunc getSolAsPolyInterp(int component, int polyOrder) const
		{
			Vector<Real> xsave = _tval;
			Vector<Real> ysave = _xval.VectorFromRow(component);
			return PolynomInterpRealFunc(xsave, ysave, polyOrder + 1);
		}
		
		/// @brief Create cubic spline interpolator for specified component
		/// @param component Component index (0 to dim-1)
		/// @return Cubic spline interpolation function
		SplineInterpRealFunc  getSolAsSplineInterp(int component) const
		{
			Vector<Real> xsave = _tval;
			Vector<Real> ysave = _xval.VectorFromRow(component);
			return SplineInterpRealFunc(xsave, ysave);
		}

		/// @brief Create 2D parametric curve from two components
		/// @param ind1 First component index (x-coordinate)
		/// @param ind2 Second component index (y-coordinate)
		/// @return Spline-interpolated parametric curve in 2D
		SplineInterpParametricCurve<2> getSolAsParamCurve2D(int ind1, int ind2) const
		{
			Matrix<Real> curve_points(_totalSavedSteps, 2);
			for (int i = 0; i < _totalSavedSteps; i++)
			{
				curve_points(i, 0) = _xval[ind1][i];
				curve_points(i, 1) = _xval[ind2][i];
			}
			return SplineInterpParametricCurve<2>(_t1, _t2, curve_points);
		}

		/// @brief Create 3D parametric curve from three components
		/// @param ind1 First component index (x-coordinate)
		/// @param ind2 Second component index (y-coordinate)
		/// @param ind3 Third component index (z-coordinate)
		/// @return Spline-interpolated parametric curve in 3D
		SplineInterpParametricCurve<3> getSolAsParamCurve3D(int ind1, int ind2, int ind3) const
		{
			Matrix<Real> curve_points(_totalSavedSteps, 3);
			for (int i = 0; i < _totalSavedSteps; i++)
			{
				curve_points(i, 0) = _xval[ind1][i];
				curve_points(i, 1) = _xval[ind2][i];
				curve_points(i, 2) = _xval[ind3][i];
			}
			return SplineInterpParametricCurve<3>(_t1, _t2, curve_points);
		}
	};
}

#endif