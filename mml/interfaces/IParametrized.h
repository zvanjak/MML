///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        IParametrized.h                                                     ///
///  Description: Shared interface for objects with adjustable parameters              ///
///               Eliminates duplication across function and ODE interfaces            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/**
 * @file IParametrized.h
 * @brief Shared base interface for parametrized mathematical objects.
 *
 * Defines the common parameter access API used across:
 * - Parametrized functions (IRealFunctionParametrized, IScalarFunctionParametrized, etc.)
 * - Parametrized ODE systems (IODESystemParametrized, IODESystemDAEParametrized)
 * - Dynamical systems (IDynamicalSystem via IODESystemParametrized)
 *
 * Provides three pure virtual methods that derived classes must implement:
 * - getNumParam() — number of parameters
 * - getParam(i) — read a single parameter
 * - setParam(i, val) — write a single parameter
 *
 * And two virtual methods with default implementations:
 * - getParams() — bulk read (calls getParam in a loop)
 * - setParams(v) — bulk write (calls setParam in a loop)
 *
 * @see IRealFunctionParametrized, IODESystemParametrized, IDynamicalSystem
 */

#if !defined MML_IPARAMETRIZED_H
#define MML_IPARAMETRIZED_H

#include "MMLBase.h"
#include "base/Vector/Vector.h"

#include <algorithm>

namespace MML
{
	/**
	 * @brief Interface for objects with adjustable real-valued parameters.
	 *
	 * This is the shared base for all parametrized mathematical objects in MML.
	 * Concrete classes only need to implement getNumParam(), getParam(), and
	 * setParam(); bulk accessors getParams()/setParams() have sensible defaults.
	 *
	 * @par Design Notes
	 * - Pure interface (no data members, no constructor).
	 * - Default getParams()/setParams() are virtual so they can be overridden
	 *   for efficiency (e.g., a single memcpy when parameters are contiguous).
	 * - Works with any class hierarchy via regular (non-virtual) inheritance;
	 *   virtual inheritance is not required because no diamond occurs in MML.
	 */
	class IParametrized
	{
	public:
		virtual ~IParametrized() = default;

		/** @brief Get the number of adjustable parameters. */
		virtual int  getNumParam() const = 0;

		/** @brief Get the i-th parameter value. */
		virtual Real getParam(int i) const = 0;

		/** @brief Set the i-th parameter to val. */
		virtual void setParam(int i, Real val) = 0;

		/** @brief Get all parameters as a vector.
		 *  Default implementation loops over getParam().
		 */
		virtual Vector<Real> getParams() const
		{
			Vector<Real> p(getNumParam());
			for (int i = 0; i < getNumParam(); ++i)
				p[i] = getParam(i);
			return p;
		}

		/** @brief Set all parameters from a vector.
		 *  Default implementation loops over setParam().
		 */
		virtual void setParams(const Vector<Real>& params)
		{
			for (int i = 0; i < std::min(getNumParam(), static_cast<int>(params.size())); ++i)
				setParam(i, params[i]);
		}
	};

} // namespace MML

#endif // MML_IPARAMETRIZED_H
