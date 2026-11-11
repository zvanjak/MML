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
#include <mml/interfaces/IFunction.h>
#include <mml/interfaces/IParametrized.h>

#include <cmath>
#include <vector>
#include <string>

namespace MML
{
	//=============================================================================
	// DISCRETE MAP INTERFACE
	//=============================================================================

	/// @brief Base interface for discrete maps x_{n+1} = f(x_n)
	///
	/// Member of the interface family (MinimalMathLibrary-wewi.6.2): a discrete map is an
	/// `IFunction<Vector<Real>, const Vector<Real>&>` (operator() == iterate) with adjustable
	/// parameters through the shared `IParametrized` mixin (getNumParam / getParam / setParam /
	/// getParams / setParams), mirroring `IODESystemParametrized` for continuous systems.
	/// Implementations live in `MML::Systems` (DiscreteMaps.h); `MML::Systems::IDiscreteMap`
	/// remains available as an alias.
	///
	/// @tparam N State dimension
	template <int N>
	class IDiscreteMap : public IFunction<Vector<Real>, const Vector<Real>&>, public IParametrized {
		static_assert(N > 0, "IDiscreteMap requires a positive state dimension");
	public:
		static constexpr int Dimension = N;

		virtual ~IDiscreteMap() = default;

		/// @brief Apply map once: x_{n+1} = f(x_n). xNext must have size N.
		virtual void map(const Vector<Real>& x, Vector<Real>& xNext) const = 0;

		/// @brief Apply map once (convenience version returning result)
		/// @throws VectorDimensionError if x.size() != N
		Vector<Real> iterate(const Vector<Real>& x) const {
			if (x.size() != N)
				throw VectorDimensionError("IDiscreteMap::iterate - state size must equal the map dimension", x.size(), N);
			Vector<Real> xNext(N);
			map(x, xNext);
			return xNext;
		}

		/// @brief Function-call form: operator()(x) == iterate(x) (IFunction contract)
		Vector<Real> operator()(const Vector<Real>& x) const override { return iterate(x); }

		/// @brief Jacobian of map Df(x)
		virtual void jacobian(const Vector<Real>& x, Matrix<Real>& J) const = 0;

		/// @brief State dimension
		int getDim() const { return N; }

		// Parameters: getNumParam() / getParam(i) / setParam(i, v) are inherited pure virtuals
		// of IParametrized; getParams() / setParams() come with checked defaults.

		/// @brief Get parameter name (default "p0", "p1", ...)
		virtual std::string getParamName(int index) const {
			if (index < 0 || index >= getNumParam())
				throw IndexError("IDiscreteMap::getParamName - parameter index out of range");
			return "p" + std::to_string(index);
		}

		/// @brief Get state name (default "x0", "x1", ...)
		virtual std::string getStateName(int index) const {
			if (index < 0 || index >= N)
				throw IndexError("IDiscreteMap::getStateName - state index out of range");
			return "x" + std::to_string(index);
		}
	};
}

namespace MML::Systems
{
	/// @deprecated IDiscreteMap moved to namespace MML (the interface family); this alias keeps
	///             existing `MML::Systems::IDiscreteMap<N>` spellings compiling.
	template <int N>
	using IDiscreteMap = MML::IDiscreteMap<N>;
}
#endif // MML_IDISCRETE_MAP_H