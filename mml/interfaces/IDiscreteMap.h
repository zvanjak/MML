///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        IDiscreteMap.h                                                      ///
///  Description: Interface for discrete maps                                         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_IDISCRETE_MAP_H
#define MML_IDISCRETE_MAP_H

#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>

#include <cmath>
#include <vector>
#include <string>

namespace MML::Systems 
{
	//=============================================================================
	// DISCRETE MAP INTERFACE
	//=============================================================================

	/// @brief Base interface for discrete maps x_{n+1} = f(x_n)
	/// @tparam N State dimension
	template <int N>
	class IDiscreteMap {
	public:
		virtual ~IDiscreteMap() = default;

		/// @brief Apply map once: x_{n+1} = f(x_n)
		virtual void map(const Vector<Real>& x, Vector<Real>& xNext) const = 0;

		/// @brief Apply map once (convenience version returning result)
		Vector<Real> iterate(const Vector<Real>& x) const {
			Vector<Real> xNext(N);
			map(x, xNext);
			return xNext;
		}

		/// @brief Jacobian of map Df(x)
		virtual void jacobian(const Vector<Real>& x, Matrix<Real>& J) const = 0;

		/// @brief State dimension
		int getDim() const { return N; }

		/// @brief Get parameter
		virtual Real getParam(int index) const = 0;

		/// @brief Set parameter
		virtual void setParam(int index, Real value) = 0;

		/// @brief Get parameter name
		virtual std::string getParamName(int index) const = 0;

		/// @brief Get state name
		virtual std::string getStateName(int index) const = 0;
	};
}
#endif // MML_IDISCRETE_MAP_H