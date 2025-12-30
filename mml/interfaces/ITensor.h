///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        ITensor.h                                                           ///
///  Description: Interface for tensor objects (contravariant/covariant indices)      ///
///               Base class for tensor algebra operations                            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_ITENSOR_H
#define MML_ITENSOR_H

#include "MMLBase.h"

namespace MML
{
	enum TensorIndexType { CONTRAVARIANT, COVARIANT };
    
	template<int N>
	class ITensor2
	{
	public:
		virtual int   NumContravar() const = 0;
		virtual int   NumCovar() const = 0;

		virtual Real  operator()(int i, int j) const = 0;
		virtual Real& operator()(int i, int j) = 0;
	};

	template<int N>
	class ITensor3
	{
	public:
		virtual int   NumContravar() const = 0;
		virtual int   NumCovar() const = 0;

		virtual Real  operator()(int i, int j, int k) const = 0;
		virtual Real& operator()(int i, int j, int k) = 0;
	};

	template<int N>
	class ITensor4
	{
	public:
		virtual int   NumContravar() const = 0;
		virtual int   NumCovar() const = 0;

		virtual Real  operator()(int i, int j, int k, int l) const = 0;
		virtual Real& operator()(int i, int j, int k, int l) = 0;
	};

	template<int N>
	class ITensor5
	{
	public:
		virtual int   NumContravar() const = 0;
		virtual int   NumCovar() const = 0;

		virtual Real  operator()(int i, int j, int k, int l, int m) const = 0;
		virtual Real& operator()(int i, int j, int k, int l, int m) = 0;
	};
}
#endif


