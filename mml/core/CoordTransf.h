///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        CoordTransf.h                                                       ///
///  Description: Coordinate transformation implementations                           ///
///               Spherical-Cartesian, Cylindrical-Cartesian conversions              ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COORD_TRANSF_H
#define MML_COORD_TRANSF_H

#include "MMLBase.h"

#include "interfaces/ICoordTransf.h"

#include "base/Vector/Vector.h"
#include "base/Vector/VectorTypes.h"
#include "base/Matrix/Matrix.h"
#include "base/Matrix/MatrixNM.h"
#include "base/Tensor.h"
#include "base/Function.h"

#include "core/Derivation.h"

namespace MML
{
	///////////////////////////////////////////////////////////////////////////////
	// COORDINATE TRANSFORMATION CONVENTIONS
	//
	//   Jacobian:       J_ij = ∂x'_i / ∂x_j  where x' = target, x = source
	//                   Matrix layout: J[row=target_dim][col=source_dim]
	//
	//   Basis vectors:  Covariant  eᵢ = ∂r/∂qⁱ  (tangent vectors, unnormalized)
	//                   Contravariant eⁱ = ∇qⁱ   (gradient 1-forms, dual basis)
	//                   Contravariant = rows of inverse Jacobian
	//
	//   Vector transforms:
	//     Contravariant components (velocity-like):  v'ⁱ = J_ij · vʲ
	//     Covariant components (gradient-like):       v'ᵢ = (J⁻¹)ʲᵢ · vⱼ
	//
	//   Indexing:       0-based throughout (q[0], q[1], q[2])
	//
	//   Numerical diff: Jacobian computed via NDer4 (4th-order central differences,
	//                   error ~h⁴) unless analytic Jacobian is provided.
	//
	//   Handedness:     Right-handed orientation assumed for all
	//                   standard coordinate systems.
	//
	// Coordinate conventions per system:
	//   Spherical:    (r, θ, φ)  — see CoordTransfSpherical.h for full details
	//   Cylindrical:  (r, φ, z)  — see CoordTransfCylindrical.h for full details
	//
	// See also: MetricTensor.h, FieldOperations.h
	///////////////////////////////////////////////////////////////////////////////

	/// @brief Base coordinate transformation class
	/// 
	/// Implements coordinate transformations between different coordinate systems.
	/// Provides basis vectors, Jacobians, and vector transformations.
	/// @tparam VectorFrom Source coordinate system vector type
	/// @tparam VectorTo Target coordinate system vector type
	/// @tparam N Dimension of coordinate space
	template<typename VectorFrom, typename VectorTo, int N>
	class CoordTransf : public virtual ICoordTransf<VectorFrom, VectorTo, N>
	{
	public:
		/// @brief Transform coordinates from source to target system
		/// @param x Position in source coordinate system
		/// @return Position in target coordinate system
		VectorN<Real, N> operator()(const VectorN<Real, N>& x) const
    {
      VectorN<Real, N> ret;
      for (int i = 0; i < N; i++)
        ret[i] = this->coordTransfFunc(i)(x);
      return ret;
    }

		/// @brief Get basis vector e_i at position (covariant basis)
		/// @param ind Index of basis vector (0 to N-1)
		/// @param pos Position in source coordinates
		/// @return Covariant basis vector in target coordinates
		virtual VectorTo   getBasisVec(int ind, const VectorFrom& pos)
		{
			VectorTo ret;

			for (int i = 0; i < N; i++)
				ret[i] = Derivation::NDer4Partial(this->coordTransfFunc(i), ind, pos);

			return ret;
		}
		/// @brief Get inverse basis vector (contravariant dual basis)
		/// @param ind Index of basis vector (0 to N-1)
		/// @param pos Position in source coordinates
		/// @return Contravariant basis vector
		virtual VectorFrom getInverseBasisVec(int ind, const VectorFrom& pos)
		{
			VectorFrom ret;

			for (int i = 0; i < N; i++)
				ret[i] = Derivation::NDer4Partial(this->coordTransfFunc(ind), i, pos);

			return ret;
		}

		/// @brief Compute Jacobian matrix of transformation at position
		/// @param pos Position in source coordinates
		/// @return Jacobian matrix J_ij = ∂x'_i/∂x_j
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

		/// @brief Transform contravariant vector (tangent vector)
		/// @param vec Vector in source coordinates
		/// @param pos Position where transformation is evaluated
		/// @return Transformed contravariant vector
		VectorTo   transfVecContravariant(const VectorFrom& vec, const VectorFrom& pos)
		{
			VectorTo ret;
			for (int i = 0; i < N; i++) {
				ret[i] = 0;
				for (int j = 0; j < N; j++)
					ret[i] += Derivation::NDer4Partial(this->coordTransfFunc(i), j, pos) * vec[j];
			}
			return ret;
		}
		/// @brief Transform covariant vector inversely (gradient-like vector)
		/// @param vec Vector in target coordinates
		/// @param pos Position in source coordinates
		/// @return Transformed covariant vector
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
 
	/// @brief Coordinate transformation with explicit inverse
	/// 
	/// Extends CoordTransf with inverse transformation capabilities.
	/// Allows transformations in both directions and tensor transformations.
	/// @tparam VectorFrom Source coordinate system vector type
	/// @tparam VectorTo Target coordinate system vector type
	/// @tparam N Dimension of coordinate space
	template<typename VectorFrom, typename VectorTo, int N>
	class CoordTransfWithInverse : public virtual CoordTransf<VectorFrom, VectorTo, N>,
																 public virtual ICoordTransfWithInverse<VectorFrom, VectorTo, N>
	{
	public:
		/// @brief Get contravariant basis vector using inverse transformation
		/// @param ind Index of basis vector (0 to N-1)
		/// @param pos Position in target coordinates
		/// @return Contravariant basis vector in source coordinates
		virtual VectorFrom getContravarBasisVec(int ind, const VectorTo& pos)
		{
			VectorFrom ret;

			for (int i = 0; i < N; i++)
				ret[i] = Derivation::NDer4Partial(this->inverseCoordTransfFunc(ind), i, pos);

			return ret;
		}
		/// @brief Get inverse contravariant basis vector
		/// @param ind Index of basis vector (0 to N-1)
		/// @param pos Position in target coordinates
		/// @return Contravariant basis vector in target coordinates
		virtual VectorTo   getInverseContravarBasisVec(int ind, const VectorTo& pos)
		{
			VectorTo ret;

			for (int i = 0; i < N; i++)
				ret[i] = Derivation::NDer4Partial(this->inverseCoordTransfFunc(i), ind, pos);

			return ret;
		}

		/// @brief Transform covariant vector (one-form, gradient-like)
		/// @param vec Vector in source coordinates
		/// @param pos Position in target coordinates where transformation is evaluated
		/// @return Transformed covariant vector
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
		/// @brief Transform contravariant vector inversely
		/// @param vec Vector in target coordinates
		/// @param pos Position in target coordinates
		/// @return Transformed contravariant vector in source coordinates
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

		/// @brief Transform rank-2 tensor (matrix-like object)
		/// @param tensor Input tensor with mixed contravariant/covariant indices
		/// @param pos Position in source coordinates
		/// @return Transformed tensor in target coordinates
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
		/// @brief Transform rank-3 tensor
		/// @param tensor Input tensor with mixed contravariant/covariant indices
		/// @param pos Position in source coordinates
		/// @return Transformed tensor in target coordinates
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
		/// @brief Transform rank-4 tensor (e.g., Riemann curvature tensor)
		/// @param tensor Input tensor with mixed contravariant/covariant indices
		/// @param pos Position in source coordinates
		/// @return Transformed tensor in target coordinates
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
		/// @brief Transform rank-5 tensor
		/// @param tensor Input tensor with mixed contravariant/covariant indices
		/// @param pos Position in source coordinates
		/// @return Transformed tensor in target coordinates
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
