///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        CoordTransf.h                                                       ///
///  Description: Coordinate transformation implementations                           ///
///               Spherical-Cartesian, Cylindrical-Cartesian conversions              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COORD_TRANSF_H
#define MML_COORD_TRANSF_H

#include "MMLBase.h"

#include "interfaces/ICoordTransf.h"

#include "base/Vector.h"
#include "base/VectorTypes.h"
#include "base/Matrix.h"
#include "base/MatrixNM.h"
#include "base/Tensor.h"
#include "base/Function.h"

#include "core/Derivation.h"

namespace MML
{
	template<typename VectorFrom, typename VectorTo, int N>
	class CoordTransf : public virtual ICoordTransf<VectorFrom, VectorTo, N>
	{
	public:
		// inherited from IVectorFunction
		VectorN<Real, N> operator()(const VectorN<Real, N>& x) const
    {
      VectorN<Real, N> ret;
      for (int i = 0; i < N; i++)
        ret[i] = this->coordTransfFunc(i)(x);
      return ret;
    }

		virtual VectorTo   getBasisVec(int ind, const VectorFrom& pos)
		{
			VectorTo ret;

			for (int i = 0; i < N; i++)
				ret[i] = Derivation::NDer4Partial(this->coordTransfFunc(i), ind, pos);

			return ret;
		}
		virtual VectorFrom getInverseBasisVec(int ind, const VectorFrom& pos)
		{
			VectorFrom ret;

			for (int i = 0; i < N; i++)
				ret[i] = Derivation::NDer4Partial(this->coordTransfFunc(ind), i, pos);

			return ret;
		}

		MatrixNM<Real, N, N> jacobian(const VectorN<Real, N>& pos)
		{
			MatrixNM<Real, N, N> jac;

			for (int i = 0; i < N; ++i)
				for (int j = 0; j < N; ++j)
				{
					jac(i, j) = Derivation::NDer4Partial(this->coordTransfFunc(i), j, pos);
				}

			return jac;
		}

		VectorTo   transfVecContravariant(const VectorFrom& vec, const VectorFrom& pos)
		{
			VectorFrom ret;
			for (int i = 0; i < N; i++) {
				ret[i] = 0;
				for (int j = 0; j < N; j++)
					ret[i] += Derivation::NDer4Partial(this->coordTransfFunc(i), j, pos) * vec[j];
			}
			return ret;
		}
		VectorFrom transfInverseVecCovariant(const VectorTo& vec, const VectorFrom& pos)
		{
			VectorFrom ret;
			for (int i = 0; i < N; i++)
			{
				ret[i] = 0;
				for (int j = 0; j < N; j++)
					ret[i] += Derivation::NDer4Partial(this->coordTransfFunc(j), i, pos) * vec[j];
			}
			return ret;
		}
	};
 
	template<typename VectorFrom, typename VectorTo, int N>
	class CoordTransfWithInverse : public virtual CoordTransf<VectorFrom, VectorTo, N>,
																 public virtual ICoordTransfWithInverse<VectorFrom, VectorTo, N>
	{
	public:
		virtual VectorFrom getContravarBasisVec(int ind, const VectorTo& pos)
		{
			VectorFrom ret;

			for (int i = 0; i < N; i++)
				ret[i] = Derivation::NDer4Partial(this->inverseCoordTransfFunc(ind), i, pos);

			return ret;
		}
		virtual VectorTo   getInverseContravarBasisVec(int ind, const VectorTo& pos)
		{
			VectorFrom ret;

			for (int i = 0; i < N; i++)
				ret[i] = Derivation::NDer4Partial(this->inverseCoordTransfFunc(i), ind, pos);

			return ret;
		}

		VectorTo   transfVecCovariant(const VectorFrom& vec, const VectorTo& pos)
		{
			VectorTo ret;
			for (int i = 0; i < N; i++)
			{
				ret[i] = 0;
				for (int j = 0; j < N; j++)
					ret[i] += Derivation::NDer4Partial(this->inverseCoordTransfFunc(j), i, pos) * vec[j];
			}

			return ret;
		}
		VectorFrom transfInverseVecContravariant(const VectorTo& vec, const VectorTo& pos)
		{
			VectorFrom ret;
			for (int i = 0; i < N; i++)
			{
				ret[i] = 0;
				for (int j = 0; j < N; j++)
					ret[i] += Derivation::NDer4Partial(this->inverseCoordTransfFunc(i), j, pos, PrecisionValues<Real>::DerivativeStepSize) * vec[j];
			}

			return ret;
		}

		Tensor2<N> transfTensor2(const Tensor2<N>& tensor, const VectorFrom& pos)
		{
			Tensor2<N> ret(tensor.NumContravar(), tensor.NumCovar());

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					ret(i, j) = 0;
					for (int k = 0; k < N; k++)
						for (int l = 0; l < N; l++)
						{
							double coef1, coef2;
							if (tensor._isContravar[0])
								coef1 = Derivation::NDer1Partial(this->coordTransfFunc(i), k, pos);
							else
								coef1 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(k), i, pos);

							if (tensor._isContravar[1])
								coef2 = Derivation::NDer1Partial(this->coordTransfFunc(j), l, pos);
							else
								coef2 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(l), j, pos);

							ret(i, j) += coef1 * coef2 * tensor(k, l);
						}
				}

			return ret;
		}
		Tensor3<N> transfTensor3(const Tensor3<N>& tensor, const VectorFrom& pos)
		{
			Tensor3<N> ret(tensor.NumContravar(), tensor.NumCovar());

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
					{
						ret.Component(i, j, k) = 0;
						for (int l = 0; l < N; l++)
							for (int m = 0; m < N; m++)
								for (int n = 0; n < N; n++)
								{
									double coef1, coef2, coef3;
									if (tensor._isContravar[0])
										coef1 = Derivation::NDer1Partial(this->coordTransfFunc(i), l, pos);
									else
										coef1 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(l), i, pos);

									if (tensor._isContravar[1])
										coef2 = Derivation::NDer1Partial(this->coordTransfFunc(j), m, pos);
									else
										coef2 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(m), j, pos);

									if (tensor._isContravar[2])
										coef3 = Derivation::NDer1Partial(this->coordTransfFunc(k), n, pos);
									else
										coef3 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(n), k, pos);

									ret.Component(i, j, k) += coef1 * coef2 * coef3 * tensor.Component(l, m, n);
								}
					}

			return ret;
		}
		Tensor4<N> transfTensor4(const Tensor4<N>& tensor, const VectorFrom& pos)
		{
			Tensor4<N> ret(tensor.NumContravar(), tensor.NumCovar());

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
						for (int l = 0; l < N; l++)
						{
							ret[i][j][k][l] = 0;
							for (int m = 0; m < N; m++)
								for (int n = 0; n < N; n++)
									for (int o = 0; o < N; o++)
										for (int p = 0; p < N; p++)
										{
											double coef1, coef2, coef3, coef4;
											if (tensor._isContravar[0])
												coef1 = Derivation::NDer1Partial(this->coordTransfFunc(i), m, pos);
											else
												coef1 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(m), i, pos);

											if (tensor._isContravar[1])
												coef2 = Derivation::NDer1Partial(this->coordTransfFunc(j), n, pos);
											else
												coef2 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(n), j, pos);

											if (tensor._isContravar[2])
												coef3 = Derivation::NDer1Partial(this->coordTransfFunc(k), o, pos);
											else
												coef3 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(o), k, pos);

											if (tensor._isContravar[3])
												coef4 = Derivation::NDer1Partial(this->coordTransfFunc(l), p, pos);
											else
												coef4 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(p), l, pos);

											ret[i][j][k][l] += coef1 * coef2 * coef3 * coef4 * tensor[m][n][o][p];
										}
						}

			return ret;
		}
		Tensor5<N> transfTensor5(const Tensor5<N>& tensor, const VectorFrom& pos)
		{
			Tensor5<N> ret(tensor.NumContravar(), tensor.NumCovar());

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
						for (int l = 0; l < N; l++)
							for (int m = 0; m < N; m++)
							{
								ret[i][j][k][l][m] = 0;
								for (int n = 0; n < N; n++)
									for (int o = 0; o < N; o++)
										for (int p = 0; p < N; p++)
											for (int q = 0; q < N; q++)
												for (int r = 0; r < N; r++)
												{
													double coef1, coef2, coef3, coef4, coef5;
													if (tensor._isContravar[0])
														coef1 = Derivation::NDer1Partial(this->coordTransfFunc(i), n, pos);
													else
														coef1 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(n), i, pos);

													if (tensor._isContravar[1])
														coef2 = Derivation::NDer1Partial(this->coordTransfFunc(j), o, pos);
													else
														coef2 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(o), j, pos);

													if (tensor._isContravar[2])
														coef3 = Derivation::NDer1Partial(this->coordTransfFunc(k), p, pos);
													else
														coef3 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(p), k, pos);

													if (tensor._isContravar[3])
														coef4 = Derivation::NDer1Partial(this->coordTransfFunc(l), q, pos);
													else
														coef4 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(q), l, pos);

													if (tensor._isContravar[4])
														coef4 = Derivation::NDer1Partial(this->coordTransfFunc(m), r, pos);
													else
														coef4 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(r), m, pos);

													ret[i][j][k][l][m] += coef1 * coef2 * coef3 * coef4 * coef5 * tensor[n][o][p][q][r];
												}
							}

			return ret;
		}
	};
}

#endif // MML_COORD_TRANSF_H
