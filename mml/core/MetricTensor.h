///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        MetricTensor.h                                                      ///
///  Description: Metric tensor calculations for curvilinear coordinates              ///
///               Jacobians, Christoffel symbols, geodesics                           ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_METRIC_TENSOR_H
#define MML_METRIC_TENSOR_H

#include "interfaces/IFunction.h"
#include "interfaces/ICoordTransf.h"

#include "base/Vector/VectorN.h"
#include "base/Matrix/MatrixNM.h"

#include "core/Derivation.h"
#include "core/Derivation/DerivationTensorField.h"
#include "core/SingularityHandling.h"

namespace MML
{
	///////////////////////////////////////////////////////////////////////////////
	// DIFFERENTIAL GEOMETRY CONVENTIONS
	//
	//   Index convention:  0-based indexing for all tensor components
	//                      (i, j, k ∈ {0, 1, ..., N-1})
	//
	//   Metric tensor:     gᵢⱼ = eᵢ · eⱼ  (covariant metric, symmetric)
	//                      gⁱʲ = inverse of gᵢⱼ  (contravariant metric)
	//                      Used to raise/lower indices.
	//
	//   Metric signature:  Positive-definite (Riemannian) assumed.
	//                      For Lorentzian metrics, see CoordTransfLorentz.h.
	//
	//   Christoffel symbols:
	//     First kind:   Γᵢⱼₖ = ½(∂gᵢₖ/∂qʲ + ∂gⱼₖ/∂qⁱ - ∂gᵢⱼ/∂qᵏ)
	//     Second kind:  Γⁱⱼₖ = gⁱˡ Γₗⱼₖ  (connection coefficients)
	//                   Symmetric in lower indices: Γⁱⱼₖ = Γⁱₖⱼ
	//
	//   Covariant derivative:  ∇ⱼvⁱ = ∂vⁱ/∂qʲ + Γⁱⱼₖ vᵏ
	//
	//   Numerical diff:  All derivatives computed via NDer4 (4th-order central
	//                    differences) unless analytic forms are provided.
	//
	// See also: CoordTransf.h, FieldOperations.h
	///////////////////////////////////////////////////////////////////////////////

	/// @brief Metric tensor field for curvilinear coordinates in differential geometry
	/// @tparam N Dimension of the manifold
	/// @note Provides covariant/contravariant metrics gᵢⱼ, Christoffel symbols Γᵢⱼₖ, covariant derivatives
	/// @note Used for geodesic equations, parallel transport, curvature calculations
	template<int N>
	class MetricTensorField : public ITensorField2<N>
	{
	public:
		/// @brief Default constructor (0 contravariant indices, 2 covariant)
		MetricTensorField() : ITensorField2<N>(0, 2) { }
		/// @brief Constructor with custom index configuration
		/// @param numContra Number of contravariant indices
		/// @param numCo Number of covariant indices
		MetricTensorField(int numContra, int numCo) : ITensorField2<N>(numContra, numCo) { }

		// implementing operator() required by IFunction interface
		virtual Tensor2<N>   operator()(const VectorN<Real, N>& pos) const override
		{
			Tensor2<N> ret(this->getNumCovar(), this->getNumContravar());

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					ret(i, j) = this->Component(i, j, pos);

			return ret;
		}

		/// @brief Get covariant metric tensor components gᵢⱼ at a point
		/// @param pos Position in coordinate space
		/// @return N×N matrix of covariant metric components
		/// @note Measures squared arc length: ds² = gᵢⱼ dxⁱ dxʲ
		// Get the covariant metric tensor components at a point (gᵢⱼ)
		MatrixNM<Real, N, N> GetCovariantMetric(const VectorN<Real, N>& pos) const
		{
			MatrixNM<Real, N, N> g_covar;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					g_covar[i][j] = this->Component(i, j, pos);
			return g_covar;
		}

		/// @brief Get contravariant metric tensor components gⁱʲ (inverse of gᵢⱼ)
		/// @param pos Position in coordinate space
		/// @return N×N matrix of contravariant metric components
		/// @note Used for raising indices: vⁱ = gⁱʲ vⱼ
		// Get the contravariant metric tensor components at a point (gⁱʲ = inverse of gᵢⱼ)
		MatrixNM<Real, N, N> GetContravariantMetric(const VectorN<Real, N>& pos) const
		{
			MatrixNM<Real, N, N> g_covar = GetCovariantMetric(pos);
			return g_covar.GetInverse();
		}

		/// @brief Raise a covariant vector index: vⁱ = gⁱʲ vⱼ
		/// @param v_covar Covariant (lower-index) vector
		/// @param pos Position where metric is evaluated
		/// @return Contravariant (upper-index) vector
		VectorN<Real, N> RaiseIndex(const VectorN<Real, N>& v_covar, const VectorN<Real, N>& pos) const
		{
			MatrixNM<Real, N, N> g_inv = GetContravariantMetric(pos);
			VectorN<Real, N> result;
			for (int i = 0; i < N; i++) {
				result[i] = 0.0;
				for (int j = 0; j < N; j++)
					result[i] += g_inv[i][j] * v_covar[j];
			}
			return result;
		}

		/// @brief Lower a contravariant vector index: vᵢ = gᵢⱼ vʲ
		/// @param v_contra Contravariant (upper-index) vector
		/// @param pos Position where metric is evaluated
		/// @return Covariant (lower-index) vector
		VectorN<Real, N> LowerIndex(const VectorN<Real, N>& v_contra, const VectorN<Real, N>& pos) const
		{
			MatrixNM<Real, N, N> g = GetCovariantMetric(pos);
			VectorN<Real, N> result;
			for (int i = 0; i < N; i++) {
				result[i] = 0.0;
				for (int j = 0; j < N; j++)
					result[i] += g[i][j] * v_contra[j];
			}
			return result;
		}

		/// @brief Raise one index of a rank-2 tensor: T^i_j = gⁱᵏ T_kj
		/// @param t Rank-2 covariant tensor T_ij (must have numContravar == 0)
		/// @param index_to_raise Which index to raise (0 or 1)
		/// @param pos Position where metric is evaluated
		/// @return Mixed tensor with one raised index
		Tensor2<N> RaiseIndex(const Tensor2<N>& t, int index_to_raise, const VectorN<Real, N>& pos) const
		{
			MatrixNM<Real, N, N> g_inv = GetContravariantMetric(pos);
			Tensor2<N> result(1, 1);  // one covariant, one contravariant

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++) {
					result(i, j) = 0.0;
					for (int k = 0; k < N; k++) {
						if (index_to_raise == 0)
							result(i, j) += g_inv[i][k] * t(k, j);
						else
							result(i, j) += g_inv[j][k] * t(i, k);
					}
				}

			return result;
		}

		/// @brief Lower one index of a rank-2 tensor: T_ij = g_ik T^k_j
		/// @param t Rank-2 tensor with at least one contravariant index
		/// @param index_to_lower Which index to lower (0 or 1)
		/// @param pos Position where metric is evaluated
		/// @return Tensor with one lowered index
		Tensor2<N> LowerIndex(const Tensor2<N>& t, int index_to_lower, const VectorN<Real, N>& pos) const
		{
			MatrixNM<Real, N, N> g = GetCovariantMetric(pos);
			Tensor2<N> result(1, 1);  // one covariant, one contravariant

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++) {
					result(i, j) = 0.0;
					for (int k = 0; k < N; k++) {
						if (index_to_lower == 0)
							result(i, j) += g[i][k] * t(k, j);
						else
							result(i, j) += g[j][k] * t(i, k);
					}
				}

			return result;
		}

		/// @brief Get Christoffel symbol of the first kind Γᵢⱼₖ (all indices lowered)
		/// @param i,j,k Indices
		/// @param pos Position in coordinate space
		/// @return Γᵢⱼₖ = gₘₖ Γᵐᵢⱼ (related to second kind via metric)
		Real GetChristoffelSymbolFirstKind(int i, int j, int k, const VectorN<Real, N>& pos) const
		{
			const MetricTensorField<N>& g = *this;

			Real gamma_ijk = 0.0;
			for (int m = 0; m < N; m++)
			{
				gamma_ijk += g.Component(m, k, pos) * GetChristoffelSymbolSecondKind(m, i, j, pos);
			}
			return gamma_ijk;
		}
		/// @brief Get Christoffel symbol of the second kind Γᵐᵢⱼ (one index raised)
		/// @param i,j,k Indices (i=contravariant, j,k=covariant)
		/// @param pos Position in coordinate space
		/// @return Γᵐᵢⱼ = ½ gᵐˡ (∂ⱼ gˡₖ + ∂ₖ gˡⱼ - ∂ˡ gⱼₖ) (connection coefficients)
		/// @note Used in geodesic equation and covariant derivatives
		Real GetChristoffelSymbolSecondKind(int i, int j, int k, const VectorN<Real, N>& pos) const
		{
			const MetricTensorField<N>& g = *this;

			// Γⁱⱼₖ = ½ gⁱˡ (∂ⱼgₗₖ + ∂ₖgₗⱼ - ∂ₗgⱼₖ)
			// Need contravariant metric gⁱˡ to raise the index
			MatrixNM<Real, N, N> g_contravar = GetContravariantMetric(pos);

			Real gamma_ijk = 0.0;
			for (int l = 0; l < N; l++)
			{
				Real coef1 = Derivation::NDer4Partial<N>(g, l, k, j, pos, nullptr);  // ∂ⱼgₗₖ
				Real coef2 = Derivation::NDer4Partial<N>(g, l, j, k, pos, nullptr);  // ∂ₖgₗⱼ
				Real coef3 = Derivation::NDer4Partial<N>(g, j, k, l, pos, nullptr);  // ∂ₗgⱼₖ

				gamma_ijk += 0.5 * g_contravar[i][l] * (coef1 + coef2 - coef3);
			}
			return gamma_ijk;
		}

		/// @brief Covariant derivative of contravariant vector: ∇ⱼ vⁱ = ∂ⱼ vⁱ + Γⁱₖⱼ vₖ
		/// @param func Vector field
		/// @param j Derivative direction
		/// @param pos Position
		/// @return Vector of covariant derivatives
		VectorN<Real, N> CovariantDerivativeContravar(const IVectorFunction<N>& func, int j, 
																									const VectorN<Real, N>& pos) const
		{
			VectorN<Real, N> ret;
			VectorN<Real, N> vec_val = func(pos);

			for (int i = 0; i < N; i++) {
				Real comp_val = Derivation::DeriveVecPartial<N>(func, i, j, pos, nullptr);

				for (int k = 0; k < N; k++)
					comp_val += GetChristoffelSymbolSecondKind(i, k, j, pos) * vec_val[k];

				ret[i] = comp_val;
			}
			return ret;
		}
		/// @brief Single component of covariant derivative (contravariant)
		Real CovariantDerivativeContravarComp(const IVectorFunction<N>& func, int i, int j, 
																					const VectorN<Real, N>& pos) const
		{
			Real ret = Derivation::DeriveVecPartial<N>(func, i, j, pos, nullptr);

			for (int k = 0; k < N; k++)
				ret += GetChristoffelSymbolSecondKind(i, k, j, pos) * func(pos)[k];

			return ret;
		}

		/// @brief Covariant derivative of covariant vector: ∇ⱼ vᵢ = ∂ⱼ vᵢ - Γₖᵢⱼ vₖ
		VectorN<Real, N> CovariantDerivativeCovar(const IVectorFunction<N>& func, int j, 
																							const VectorN<Real, N>& pos) const
		{
			VectorN<Real, N> ret;
			VectorN<Real, N> vec_val = func(pos);

			for (int i = 0; i < N; i++) {
				Real comp_val = Derivation::DeriveVecPartial<N>(func, i, j, pos, nullptr);

				for (int k = 0; k < N; k++)
					comp_val -= GetChristoffelSymbolSecondKind(k, i, j, pos) * vec_val[k];

				ret[i] = comp_val;
			}
			return ret;
		}
		/// @brief Single component of covariant derivative (covariant)
		Real CovariantDerivativeCovarComp(const IVectorFunction<N>& func, int i, int j, 
																			const VectorN<Real, N>& pos) const
		{
			Real comp_val = Derivation::DeriveVecPartial<N>(func, i, j, pos, nullptr);

			for (int k = 0; k < N; k++)
				comp_val -= GetChristoffelSymbolSecondKind(k, i, j, pos) * func(pos)[k];

			return comp_val;
		}
	};

	/// @brief Flat metric tensor for Cartesian 3D (gᵢⱼ = δᵢⱼ, all Christoffel symbols vanish)
	class MetricTensorCartesian3D : public MetricTensorField<3>
	{
	public:
		MetricTensorCartesian3D() : MetricTensorField<3>(0, 2) { }

		Real Component(int i, int j, const VectorN<Real, 3>& pos) const
		{
			if (i == j)
				return 1.0;
			else
				return 0.0;
		}
	};

	/// @brief Metric for spherical coords (r,θ,φ): ds²=dr²+r²dθ²+r²sin²θdφ²
	/// @note Diagonal: g=diag(1, r², r²sin²θ)
	/// @brief Metric for spherical coords (r,θ,φ): ds²=dr²+r²dθ²+r²sin²θdφ²
	/// @note Diagonal: g=diag(1, r², r²sin²θ)
	class MetricTensorSpherical : public MetricTensorField<3>
	{
	public:
		MetricTensorSpherical() : MetricTensorField<3>(0, 2) { }

		virtual  Real Component(int i, int j, const VectorN<Real, 3>& pos) const override
		{
			if (i == 0 && j == 0)
				return 1.0;
			else if (i == 1 && j == 1)
				return POW2(pos[0]);
			else if (i == 2 && j == 2)
				return pos[0] * pos[0] * sin(pos[1]) * sin(pos[1]);
			else
				return 0.0;
		}
	};
	/// @brief Contravariant spherical metric: gⁱʲ=diag(1, 1/r², 1/(r²sin²θ))
	class MetricTensorSphericalContravar : public MetricTensorField<3>
	{
	public:
		MetricTensorSphericalContravar() : MetricTensorField<3>(2, 0) { }

		virtual Real Component(int i, int j, const VectorN<Real, 3>& pos) const
		{
			if (i == 0 && j == 0)
				return 1.0;
			else if (i == 1 && j == 1)
			{
				Real r2 = pos[0] * pos[0];
				return Singularity::SafeDivide(1.0, r2, SingularityPolicy::Throw,
					"MetricTensorSphericalContravar g^11: r=0 singularity");
			}
			else if (i == 2 && j == 2)
			{
				Real r2sin2 = pos[0] * pos[0] * std::sin(pos[1]) * std::sin(pos[1]);
				return Singularity::SafeDivide(1.0, r2sin2, SingularityPolicy::Throw,
					"MetricTensorSphericalContravar g^22: r=0 or theta=0/pi singularity");
			}
			else
				return 0.0;
		}
	};
	/// @brief Metric for cylindrical coords (ρ,φ,z): ds²=dρ²+ρ²dφ²+dz², g=diag(1,ρ²,1)
	class MetricTensorCylindrical : public MetricTensorField<3>
	{
	public:
		MetricTensorCylindrical() : MetricTensorField<3>(0, 2) { }

		virtual Real Component(int i, int j, const VectorN<Real, 3>& pos) const
		{
			if (i == 0 && j == 0)
				return 1.0;
			else if (i == 1 && j == 1)
				return pos[0] * pos[0];
			else if (i == 2 && j == 2)
				return 1.0;
			else
				return 0.0;
		}
	};

	/// @brief Compute metric from coord transformation Jacobian: gᵢⱼ=∂xₖ/∂ξⁱ ∂xₖ/∂ξʲ
	template<typename VectorFrom, typename VectorTo, int N>
	class MetricTensorFromCoordTransf : public MetricTensorField<N>
	{
		const ICoordTransfWithInverse<VectorFrom, VectorTo, N>& _coordTransf;

	public:
		MetricTensorFromCoordTransf(ICoordTransfWithInverse<VectorFrom, VectorTo, N>& inTransf) : _coordTransf(inTransf)
		{ }

		virtual Real Component(int i, int j, const VectorN<Real, N>& pos) const
		{
			Real g_ij = 0.0;
			for (int k = 0; k < N; k++)
			{
				auto der_k_by_i = Derivation::DerivePartial<N>(_coordTransf.coordTransfFunc(k), i, pos, nullptr);
				auto der_k_by_j = Derivation::DerivePartial<N>(_coordTransf.coordTransfFunc(k), j, pos, nullptr);

				g_ij += der_k_by_i * der_k_by_j;
			}
			return g_ij;
		}
	};

	/// @brief Minkowski metric for special relativity: η=diag(-1,1,1,1), signature (−,+,+,+)
	/// @note Flat spacetime with coords (ct,x,y,z), ds²=-c²dt²+dx²+dy²+dz²
	class MetricTensorMinkowski : public MetricTensorField<4>
	{
	public:
		MetricTensorMinkowski() : MetricTensorField<4>(0, 2) {}

		virtual Real Component(int i, int j, const VectorN<Real, 4>& pos) const override
		{
			if (i == 0 && j == 0)
				return -1.0;
			else if (i == 1 && j == 1)
				return 1.0;
			else if (i == 2 && j == 2)
				return 1.0;
			else if (i == 3 && j == 3)
				return 1.0;
			else
				return 0.0;
		}
	};
}
#endif