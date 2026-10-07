///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        QuasiRandom.h                                                       ///
///  Description: Low-discrepancy (quasi-random) sequences: Halton and Sobol          ///
///               For quasi-Monte Carlo integration and sampling                       ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

#if !defined MML_QUASI_RANDOM_H
#define MML_QUASI_RANDOM_H

#include <mml/MMLBase.h>
#include <mml/base/Vector/VectorN.h>

namespace MML
{
	///////////////////////////////////////////////////////////////////////////
	///                          HaltonSequence                             ///
	///////////////////////////////////////////////////////////////////////////

	/// @brief Halton low-discrepancy sequence (radical inverse in coprime bases).
	/// @details Dimension d uses the (d+1)-th prime as its base. Supports up to 40
	///          dimensions; higher dimensions of plain Halton correlate, so prefer
	///          Sobol for small dimensions.
	class HaltonSequence
	{
		static constexpr int kMaxDim = 40;
		int  _dim;
		long _index;

		static int prime(int d)
		{
			static const int primes[kMaxDim] = {
				2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
				73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173 };
			return primes[d];
		}
	public:
		static constexpr int MaxDim = kMaxDim;

		/// @brief Construct a Halton generator for `dim` dimensions (index starts at 1 by default).
		explicit HaltonSequence(int dim, long startIndex = 1) : _dim(dim), _index(startIndex)
		{
			if (dim < 1 || dim > kMaxDim)
				throw ArgumentError("HaltonSequence: dim must be in [1, 40]");
			if (startIndex < 0)
				throw ArgumentError("HaltonSequence: startIndex must be >= 0");
		}

		int dim() const { return _dim; }
		void Reset(long startIndex = 1) { _index = startIndex; }

		/// @brief Van der Corput radical inverse of `i` in the given base.
		static Real RadicalInverse(long i, int base)
		{
			Real f = REAL(1.0) / base;
			Real r = REAL(0.0);
			while (i > 0) { r += f * (i % base); i /= base; f /= base; }
			return r;
		}

		/// @brief Fill out[0.._dim-1] with the next point in [0,1)^dim.
		void Next(Real* out)
		{
			for (int d = 0; d < _dim; ++d) out[d] = RadicalInverse(_index, prime(d));
			++_index;
		}

		/// @brief Return the next point as a fixed-size vector (N must equal dim()).
		template<int N>
		VectorN<Real, N> Next()
		{
			static_assert(N >= 1, "HaltonSequence::Next<N> requires N >= 1");
			VectorN<Real, N> v;
			Next(&v[0]);
			return v;
		}
	};

	///////////////////////////////////////////////////////////////////////////
	///                           SobolSequence                             ///
	///////////////////////////////////////////////////////////////////////////

	/// @brief Sobol low-discrepancy sequence (Gray-code construction, up to 6 dimensions).
	/// @details Direction numbers follow the classic six-dimensional set (Numerical Recipes).
	///          Each instance keeps its own state, so multiple generators are independent.
	class SobolSequence
	{
		static const int kMaxBit = 30;
		static const int kMaxDim = 6;

		int      _dim;
		unsigned _in;
		unsigned _ix[kMaxDim];
		unsigned _iv[kMaxDim * kMaxBit];
		Real     _fac;

		void initialize()
		{
			static const int      mdeg[kMaxDim] = { 1, 2, 3, 3, 4, 4 };
			static const unsigned ip[kMaxDim]   = { 0, 1, 1, 2, 1, 4 };
			// Initial direction numbers, laid out row-major as _iv[row * kMaxDim + dim].
			static const unsigned iv0[24] = {
				1, 1, 1, 1, 1, 1,
				3, 1, 3, 3, 1, 1,
				5, 7, 7, 3, 3, 5,
				15, 11, 5, 15, 13, 9 };

			for (int i = 0; i < kMaxDim * kMaxBit; ++i) _iv[i] = 0;
			for (int i = 0; i < 24; ++i) _iv[i] = iv0[i];
			for (int k = 0; k < kMaxDim; ++k) _ix[k] = 0;
			_in = 0;
			_fac = REAL(1.0) / (unsigned(1) << kMaxBit);

			// _iv[j*kMaxDim + k] is direction number for bit row j, dimension k.
			for (int k = 0; k < kMaxDim; ++k)
			{
				for (int j = 0; j < mdeg[k]; ++j)
					_iv[j * kMaxDim + k] <<= (kMaxBit - 1 - j);

				for (int j = mdeg[k]; j < kMaxBit; ++j)
				{
					unsigned ipp = ip[k];
					unsigned i = _iv[(j - mdeg[k]) * kMaxDim + k];
					i ^= (i >> mdeg[k]);
					for (int l = mdeg[k] - 1; l >= 1; --l)
					{
						if (ipp & 1u) i ^= _iv[(j - l) * kMaxDim + k];
						ipp >>= 1;
					}
					_iv[j * kMaxDim + k] = i;
				}
			}
		}
	public:
		static constexpr int MaxDim = kMaxDim;

		/// @brief Construct a Sobol generator for `dim` dimensions (1..6).
		explicit SobolSequence(int dim) : _dim(dim)
		{
			if (dim < 1 || dim > kMaxDim)
				throw ArgumentError("SobolSequence: dim must be in [1, 6]");
			initialize();
		}

		int dim() const { return _dim; }
		void Reset() { initialize(); }

		/// @brief Fill out[0.._dim-1] with the next point in [0,1)^dim.
		void Next(Real* out)
		{
			unsigned im = _in++;
			int j = 0;
			for (; j < kMaxBit; ++j) { if (!(im & 1u)) break; im >>= 1; }
			if (j >= kMaxBit)
				throw ArgumentError("SobolSequence: sequence exhausted (2^30 points)");

			im = static_cast<unsigned>(j) * kMaxDim;
			for (int k = 0; k < _dim; ++k)
			{
				_ix[k] ^= _iv[im + k];
				out[k] = _ix[k] * _fac;
			}
		}

		/// @brief Return the next point as a fixed-size vector (1 <= N <= 6).
		template<int N>
		VectorN<Real, N> Next()
		{
			static_assert(N >= 1 && N <= kMaxDim, "SobolSequence::Next<N> requires 1 <= N <= 6");
			VectorN<Real, N> v;
			Next(&v[0]);
			return v;
		}
	};
}
#endif
