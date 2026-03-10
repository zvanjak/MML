///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DynamicalSystemBase.h                                               ///
///  Description: Base class template for dynamical systems                           ///
///               Provides common infrastructure for N-dimensional systems            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DYNAMICAL_SYSTEM_BASE_H
#define MML_DYNAMICAL_SYSTEM_BASE_H

#include "MMLBase.h"
#include "base/Vector/Vector.h"
#include "interfaces/IDynamicalSystem.h"

#include <vector>
#include <string>
#include <utility>

namespace MML::Systems 
{
	//=============================================================================
	// DYNAMICAL SYSTEM BASE CLASS
	//=============================================================================

	/// @brief Template base class for dynamical systems
	///
	/// Provides common infrastructure for continuous dynamical systems:
	/// - Parameter storage and access
	/// - State and parameter naming
	/// - Parameter range specifications
	///
	/// @tparam N State space dimension
	/// @tparam P Number of parameters
	template<int N, int P>
	class DynamicalSystemBase : public IDynamicalSystem {
	protected:
		std::vector<Real> _params;
		std::vector<std::string> _stateNames;
		std::vector<std::string> _paramNames;
		std::vector<std::pair<Real, Real>> _paramRanges;

	public:
		DynamicalSystemBase()
				: _params(P, 0.0)
				, _stateNames(N)
				, _paramNames(P)
				, _paramRanges(P, std::make_pair(-1e10, 1e10)) {
			for (int i = 0; i < N; ++i)
				_stateNames[i] = "x" + std::to_string(i);
			for (int i = 0; i < P; ++i)
				_paramNames[i] = "p" + std::to_string(i);
		}

		int getDim() const override { return N; }
		int getNumParam() const override { return P; }

		Real getParam(int i) const override { return _params[i]; }
		void setParam(int i, Real val) override { _params[i] = val; }

		Vector<Real> getParams() const override {
			Vector<Real> p(P);
			for (int i = 0; i < P; ++i)
				p[i] = _params[i];
			return p;
		}

		void setParams(const Vector<Real>& p) override {
			for (int i = 0; i < P; ++i)
				_params[i] = p[i];
		}

		std::string getStateName(int i) const override { return _stateNames[i]; }
		std::string getParamName(int i) const override { return _paramNames[i]; }
		std::pair<Real, Real> getParamRange(int i) const override { return _paramRanges[i]; }
	};

} // namespace MML::Systems
#endif // MML_DYNAMICAL_SYSTEM_BASE_H
