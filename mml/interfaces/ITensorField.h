///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ITensorField.h                                                      ///
///  Description: Interface for tensor-valued fields (scalar, vector, tensor fields)  ///
///               Base classes for field operations on manifolds                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_ITENSOR_FIELD_H
#define MML_ITENSOR_FIELD_H

#include "MMLBase.h"

#include "IFunction.h"

#include "base/Vector.h"
#include "base/VectorN.h"
#include "base/Tensor.h"

namespace MML
{
	//////////////////////////////////////////////////////////////////////
	template<int N>
	class ITensorField2 : public IFunction<Tensor2<N>, const VectorN<Real, N>& >
	{
		int _numContravar;
		int _numCovar;
	public:
		ITensorField2(int numContra, int numCo) : _numContravar(numContra), _numCovar(numCo) {}

		int getNumContravar() const { return _numContravar; }
		int getNumCovar()			const { return _numCovar; }

		// concrete implementations need to provide this 
		virtual Real Component(int i, int j, const VectorN<Real, N>& pos) const = 0;

		virtual ~ITensorField2() {}
	};

	template<int N>
	class ITensorField3 : public IFunction<Tensor3<N>, const VectorN<Real, N>& >
	{
		int _numContravar;
		int _numCovar;
	public:
		int getNumContravar() const { return _numContravar; }
		int getNumCovar() const { return _numCovar; }

		virtual Real    Component(int i, int j, int k, const VectorN<Real, N>& pos) const = 0;

		virtual ~ITensorField3() {}
	};

	template<int N>
	class ITensorField4 : public IFunction<Tensor4<N>, const VectorN<Real, N>& >
	{
		int _numContravar;
		int _numCovar;
	public:
		int getNumContravar() const { return _numContravar; }
		int getNumCovar() const { return _numCovar; }

		virtual Real    Component(int i, int j, int k, int l, const VectorN<Real, N>& pos) const = 0;

		virtual ~ITensorField4() {}
	};

	template<int N>
	class ITensorField5 : public IFunction<Tensor5<N>, const VectorN<Real, N>& >
	{
		int _numContravar;
		int _numCovar;
	public:
		int getNumContravar() const { return _numContravar; }
		int getNumCovar() const { return _numCovar; }

		virtual Real    Component(int i, int j, int k, int l, int m, const VectorN<Real, N>& pos) const = 0;
		virtual ~ITensorField5() {}
	};
}
#endif