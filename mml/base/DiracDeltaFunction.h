///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DiracDeltaFunction.h                                                ///
///  Description: Dirac delta function approximations for discrete sampling           ///
///               Gaussian and sinc approximations for numerical integration          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DELTA_FUNCTION_H
#define MML_DELTA_FUNCTION_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

namespace MML
{
	class DiracFunction : public IRealFunction
	{
	protected:
		int _N;
	public:
		DiracFunction(int N) : _N(N) {}
	};

	class DiracStep : public DiracFunction
	{
	public:
		DiracStep(int N) : DiracFunction(N) {}

		Real operator()(const Real x) const
		{
			if (x < -1.0 / (2 * _N) || x > 1.0 / (2 * _N))
				return 0.0;
			else
				return _N;
		}
	};
	class DiracExp : public DiracFunction
	{
	public:
		DiracExp(int N) : DiracFunction(N) {}

		Real operator()(const Real x) const { return _N / sqrt(2 * Constants::PI) * exp(-x * x * _N * _N); }
	};
	class DiracSqr : public DiracFunction
	{
	public:
		DiracSqr(int N) : DiracFunction(N) {}

		Real operator()(const Real x) const { return _N / Constants::PI / (1 + _N * _N * x * x); }
	};
	class DiracSin : public DiracFunction
	{
	public:
		DiracSin(int N) : DiracFunction(N) {}

		Real operator()(const Real x) const { return sin(_N * x) / (Constants::PI * x); }
	};
}

#endif
