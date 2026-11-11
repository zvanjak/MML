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

#include <mml/algorithms/MatrixSingularity.h>
#include <mml/MMLBase.h>
#include <mml/MMLSingularityHandling.h>

#include <mml/interfaces/ICoordTransf.h>

#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/Tensor.h>
#include <mml/base/Function.h>

#include <mml/core/Derivation.h>

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
	template<typename VectorFrom, typename VectorTo, int N = Detail::VectorDimension<VectorFrom>>
	class CoordTransf : public virtual ICoordTransf<VectorFrom, VectorTo, N>
	{
	public:
		using Base = ICoordTransf<VectorFrom, VectorTo, N>;
		using Base::FromDim;   ///< Dimension of the source space (deduced from VectorFrom)
		using Base::ToDim;     ///< Dimension of the target space (deduced from VectorTo)
		/// Jacobian J_ij = dx'_i/dx_j of the forward map is ToDim x FromDim
		using JacobianType        = MatrixNM<Real, ToDim, FromDim>;
		/// Jacobian K_ij = dx_i/dx'_j of the inverse map is FromDim x ToDim
		using InverseJacobianType = MatrixNM<Real, FromDim, ToDim>;

		/// @brief Transform coordinates component-wise through coordTransfFunc(i)
		/// @param x Position in source coordinate system
		/// @return Position in target coordinate system
		/// @note Must be consistent with transf() — both use coordTransfFunc(i)
		VectorTo operator()(const VectorFrom& x) const override
		{
			VectorTo ret;
			for (int i = 0; i < ToDim; i++)
				ret[i] = this->coordTransfFunc(i)(x);
			return ret;
		}

		/// @brief Get basis vector e_i at position (covariant basis)
		/// @param ind Index of basis vector (0 to N-1)
		/// @param pos Position in source coordinates
		/// @return Covariant basis vector in target coordinates
		virtual VectorTo   getBasisVec(int ind, const VectorFrom& pos) const
		{
			const auto jac = jacobian(pos);
			VectorTo ret;

			for (int i = 0; i < ToDim; i++)
				ret[i] = jac(i, ind);

			return ret;
		}
		/// @brief Get inverse basis vector (contravariant dual basis)
		/// @param ind Index of basis vector (0 to N-1)
		/// @param pos Position in source coordinates
		/// @return Contravariant basis vector
		virtual VectorFrom getInverseBasisVec(int ind, const VectorFrom& pos) const
		{
			const auto jac = jacobian(pos);
			VectorFrom ret;

			for (int i = 0; i < FromDim; i++)
				ret[i] = jac(ind, i);

			return ret;
		}

		/// @brief Compute Jacobian matrix of transformation at position
		/// @param pos Position in source coordinates
		/// @return Jacobian matrix J_ij = ∂x'_i/∂x_j (ToDim x FromDim)
		virtual JacobianType jacobian(const VectorN<Real, FromDim>& pos) const
		{
			JacobianType jac;

			for (int i = 0; i < ToDim; ++i)
				for (int j = 0; j < FromDim; ++j)
				{
					jac(i, j) = Derivation::NDer4Partial(this->coordTransfFunc(i), j, pos);
				}

			return jac;
		}

		/// @brief Transform contravariant vector (tangent vector)
		/// @param vec Vector in source coordinates
		/// @param pos Position where transformation is evaluated
		/// @return Transformed contravariant vector
		VectorTo   transfVecContravariant(const VectorFrom& vec, const VectorFrom& pos) const
		{
			const auto jac = jacobian(pos);
			VectorTo ret;
			for (int i = 0; i < ToDim; i++) {
				ret[i] = 0;
				for (int j = 0; j < FromDim; j++)
					ret[i] += jac(i, j) * vec[j];
			}
			return ret;
		}
		/// @brief Transform covariant vector inversely (gradient-like vector)
		/// @param vec Vector in target coordinates
		/// @param pos Position in source coordinates
		/// @return Transformed covariant vector
		VectorFrom transfInverseVecCovariant(const VectorTo& vec, const VectorFrom& pos) const
		{
			const auto jac = jacobian(pos);
			VectorFrom ret;
			for (int i = 0; i < FromDim; i++)
			{
				ret[i] = 0;
				for (int j = 0; j < ToDim; j++)
					ret[i] += jac(j, i) * vec[j];
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
	template<typename VectorFrom, typename VectorTo, int N = Detail::VectorDimension<VectorFrom>>
	class CoordTransfWithInverse : public virtual CoordTransf<VectorFrom, VectorTo, N>,
																 public virtual ICoordTransfWithInverse<VectorFrom, VectorTo, N>
	{
	public:
		using CoordTransf<VectorFrom, VectorTo, N>::FromDim;
		using CoordTransf<VectorFrom, VectorTo, N>::ToDim;
		using typename CoordTransf<VectorFrom, VectorTo, N>::InverseJacobianType;

		/// @brief Compute Jacobian matrix of the inverse transformation at a target position
		/// @param pos Position in target coordinates
		/// @return Inverse Jacobian K_ij = ∂x_i/∂x'_j (FromDim x ToDim)
		virtual InverseJacobianType inverseJacobian(const VectorN<Real, ToDim>& pos) const
		{
			Singularity::RequireFiniteCoordinates<ToDim>(pos, "Inverse Jacobian");
			InverseJacobianType jac;

			for (int i = 0; i < FromDim; ++i)
				for (int j = 0; j < ToDim; ++j)
					jac(i, j) = Derivation::NDer4Partial(this->inverseCoordTransfFunc(i), j, pos);

			Singularity::RequireFiniteMatrix(jac, "Inverse Jacobian");
			if constexpr (FromDim == ToDim)
				(void)Singularity::InverseMatrix(jac);
			return jac;
		}

		/// Explicit singular-matrix policy; retains dispatch to one-argument user overrides.
		InverseJacobianType inverseJacobian(const VectorN<Real, ToDim>& pos, SingularityPolicy policy) const {
			Singularity::RequireMatrixPolicy(policy);
			Singularity::RequireFiniteCoordinates<ToDim>(pos, "Inverse Jacobian");
			try {
				const auto jac = inverseJacobian(pos);
				Singularity::RequireFiniteMatrix(jac, "Inverse Jacobian");
				return jac;
			} catch (const SingularMatrixError&) {
				if (policy == SingularityPolicy::Throw)
					throw;
				return Singularity::MatrixAtSingularity<InverseJacobianType>(policy, "Inverse Jacobian: singular parametrization");
			}
		}

		/// @brief Get contravariant basis vector using inverse transformation
		/// @param ind Index of basis vector (0 to N-1)
		/// @param pos Position in target coordinates
		/// @return Contravariant basis vector in source coordinates
		/// @note Row ind of K has ToDim entries but is returned as a VectorFrom; only meaningful
		///       for square transformations (the original, unchanged semantics).
		VectorFrom getContravarBasisVec(int ind, const VectorTo& pos) const
			requires (FromDim == ToDim)
		{
			const auto jac = inverseJacobian(pos);
			VectorFrom ret;

			for (int i = 0; i < FromDim; i++)
				ret[i] = jac(ind, i);

			return ret;
		}
		/// @brief Get inverse contravariant basis vector
		/// @param ind Index of basis vector (0 to N-1)
		/// @param pos Position in target coordinates
		/// @return Contravariant basis vector in target coordinates
		/// @note Column ind of K has FromDim entries but is returned as a VectorTo; only meaningful
		///       for square transformations (the original, unchanged semantics).
		VectorTo   getInverseContravarBasisVec(int ind, const VectorTo& pos) const
			requires (FromDim == ToDim)
		{
			const auto jac = inverseJacobian(pos);
			VectorTo ret;

			for (int i = 0; i < ToDim; i++)
				ret[i] = jac(i, ind);

			return ret;
		}

		/// @brief Transform covariant vector (one-form, gradient-like)
		/// @param vec Vector in source coordinates
		/// @param pos Position in target coordinates where transformation is evaluated
		/// @return Transformed covariant vector
		VectorTo   transfVecCovariant(const VectorFrom& vec, const VectorTo& pos) const
		{
			const auto jac = inverseJacobian(pos);
			VectorTo ret;
			for (int i = 0; i < ToDim; i++)
			{
				ret[i] = 0;
				for (int j = 0; j < FromDim; j++)
					ret[i] += jac(j, i) * vec[j];
			}

			return ret;
		}
		/// @brief Transform contravariant vector inversely
		/// @param vec Vector in target coordinates
		/// @param pos Position in target coordinates
		/// @return Transformed contravariant vector in source coordinates
		VectorFrom transfInverseVecContravariant(const VectorTo& vec, const VectorTo& pos) const
		{
			const auto jac = inverseJacobian(pos);
			VectorFrom ret;
			for (int i = 0; i < FromDim; i++)
			{
				ret[i] = 0;
				for (int j = 0; j < ToDim; j++)
					ret[i] += jac(i, j) * vec[j];
			}

			return ret;
		}

		/// @brief Transform rank-2 tensor (matrix-like object)
		/// @param tensor Input tensor with mixed contravariant/covariant indices
		/// @param pos Position in source coordinates
		/// @return Transformed tensor in target coordinates
		Tensor2<N> transfTensor2(const Tensor2<N>& tensor, const VectorFrom& pos) const
			requires (FromDim == ToDim)
		{
			const auto forwardJac = this->jacobian(pos);
			const auto inverseJac = inverseJacobian(this->transf(pos));
			// Tensor ctor signature is (nCovar, nContravar); preserve per-slot variance of the input
			Tensor2<N> ret(tensor.numCovar(), tensor.numContravar());
			for (int s = 0; s < 2; s++)
				ret.setContravar(s, tensor.isContravar(s));

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					ret(i, j) = 0;
					for (int k = 0; k < N; k++)
						for (int l = 0; l < N; l++)
						{
              Real coef1, coef2;
							if (tensor.isContravar(0))
								coef1 = forwardJac(i, k);
							else
								coef1 = inverseJac(k, i);

							if (tensor.isContravar(1))
								coef2 = forwardJac(j, l);
							else
								coef2 = inverseJac(l, j);

							ret(i, j) += coef1 * coef2 * tensor(k, l);
						}
				}

			return ret;
		}
		/// @brief Transform rank-3 tensor
		/// @param tensor Input tensor with mixed contravariant/covariant indices
		/// @param pos Position in source coordinates
		/// @return Transformed tensor in target coordinates
		Tensor3<N> transfTensor3(const Tensor3<N>& tensor, const VectorFrom& pos) const
			requires (FromDim == ToDim)
		{
			const auto forwardJac = this->jacobian(pos);
			const auto inverseJac = inverseJacobian(this->transf(pos));
			// Tensor ctor signature is (nCovar, nContravar); preserve per-slot variance of the input
			Tensor3<N> ret(tensor.numCovar(), tensor.numContravar());
			for (int s = 0; s < 3; s++)
				ret.setContravar(s, tensor.isContravar(s));

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
					{
						ret(i, j, k) = 0;
						for (int l = 0; l < N; l++)
							for (int m = 0; m < N; m++)
								for (int n = 0; n < N; n++)
								{
									Real coef1, coef2, coef3;
									if (tensor.isContravar(0))
										coef1 = forwardJac(i, l);
									else
										coef1 = inverseJac(l, i);

									if (tensor.isContravar(1))
										coef2 = forwardJac(j, m);
									else
										coef2 = inverseJac(m, j);

									if (tensor.isContravar(2))
										coef3 = forwardJac(k, n);
									else
										coef3 = inverseJac(n, k);

									ret(i, j, k) += coef1 * coef2 * coef3 * tensor(l, m, n);
								}
					}

			return ret;
		}
		/// @brief Transform rank-4 tensor (e.g., Riemann curvature tensor)
		/// @param tensor Input tensor with mixed contravariant/covariant indices
		/// @param pos Position in source coordinates
		/// @return Transformed tensor in target coordinates
		Tensor4<N> transfTensor4(const Tensor4<N>& tensor, const VectorFrom& pos) const
			requires (FromDim == ToDim)
		{
			const auto forwardJac = this->jacobian(pos);
			const auto inverseJac = inverseJacobian(this->transf(pos));
			// Tensor ctor signature is (nCovar, nContravar); preserve per-slot variance of the input
			Tensor4<N> ret(tensor.numCovar(), tensor.numContravar());
			for (int s = 0; s < 4; s++)
				ret.setContravar(s, tensor.isContravar(s));

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
											Real coef1, coef2, coef3, coef4;
											if (tensor.isContravar(0))
												coef1 = forwardJac(i, m);
											else
												coef1 = inverseJac(m, i);

											if (tensor.isContravar(1))
												coef2 = forwardJac(j, n);
											else
												coef2 = inverseJac(n, j);

											if (tensor.isContravar(2))
												coef3 = forwardJac(k, o);
											else
												coef3 = inverseJac(o, k);

											if (tensor.isContravar(3))
												coef4 = forwardJac(l, p);
											else
												coef4 = inverseJac(p, l);

											ret[i][j][k][l] += coef1 * coef2 * coef3 * coef4 * tensor[m][n][o][p];
										}
						}

			return ret;
		}
		/// @brief Transform rank-5 tensor
		/// @param tensor Input tensor with mixed contravariant/covariant indices
		/// @param pos Position in source coordinates
		/// @return Transformed tensor in target coordinates
		Tensor5<N> transfTensor5(const Tensor5<N>& tensor, const VectorFrom& pos) const
			requires (FromDim == ToDim)
		{
			const auto forwardJac = this->jacobian(pos);
			const auto inverseJac = inverseJacobian(this->transf(pos));
			// Tensor ctor signature is (nCovar, nContravar); preserve per-slot variance of the input
			Tensor5<N> ret(tensor.numCovar(), tensor.numContravar());
			for (int s = 0; s < 5; s++)
				ret.setContravar(s, tensor.isContravar(s));

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
														Real coef1, coef2, coef3, coef4, coef5;
													if (tensor.isContravar(0))
													coef1 = forwardJac(i, n);
													else
													coef1 = inverseJac(n, i);

													if (tensor.isContravar(1))
													coef2 = forwardJac(j, o);
													else
													coef2 = inverseJac(o, j);

													if (tensor.isContravar(2))
													coef3 = forwardJac(k, p);
													else
													coef3 = inverseJac(p, k);

													if (tensor.isContravar(3))
													coef4 = forwardJac(l, q);
													else
													coef4 = inverseJac(q, l);

													if (tensor.isContravar(4))
												    coef5 = forwardJac(m, r);
											    else
												    coef5 = inverseJac(r, m);
													
                          ret[i][j][k][l][m] += coef1 * coef2 * coef3 * coef4 * coef5 * tensor[n][o][p][q][r];
												}
							}

			return ret;
		}
	};
}

#endif // MML_COORD_TRANSF_H
