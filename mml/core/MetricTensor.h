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

#include <mml/MMLSingularityHandling.h>

#include <mml/interfaces/IFunction.h>

#include <mml/base/Vector/VectorN.h>
#include <mml/base/Matrix/MatrixNM.h>

#include <mml/core/Derivation.h>
#include <mml/core/Derivation/DerivationTensorField.h>
#include <mml/core/CoordTransf/CoordTransfBase.h>

#include <cmath>

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
	//   Metric signature:  Positive-definite (Riemannian) by default.
	//                      Lorentzian metrics inherit LorentzianMetric<N> and use
	//                      signature (-,+,+,...) for interval classification.
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
		enum class SignatureType
		{
			Riemannian,
			Lorentzian
		};

		/// @brief Default constructor (0 contravariant indices, 2 covariant)
		MetricTensorField() : ITensorField2<N>(0, 2) { }
		/// @brief Constructor with custom index configuration
		/// @param numContra Number of contravariant indices
		/// @param numCo Number of covariant indices
		MetricTensorField(int numContra, int numCo) : ITensorField2<N>(numContra, numCo) { }

		virtual SignatureType Signature() const { return SignatureType::Riemannian; }
		bool IsRiemannian() const { return Signature() == SignatureType::Riemannian; }
		bool IsLorentzian() const { return Signature() == SignatureType::Lorentzian; }

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
		virtual MatrixNM<Real, N, N> GetCovariantMetric(const VectorN<Real, N>& pos) const
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

		/// @brief Get Riemann curvature tensor component R^rho_sigma_mu_nu.
		/// @param rho Contravariant index
		/// @param sigma,mu,nu Covariant indices
		/// @param pos Position in coordinate space
		/// @return R^rho_sigma_mu_nu = partial_mu Gamma^rho_sigma_nu - partial_nu Gamma^rho_sigma_mu + Gamma^rho_lambda_mu Gamma^lambda_sigma_nu - Gamma^rho_lambda_nu Gamma^lambda_sigma_mu
		Real GetRiemannCurvatureTensor(int rho, int sigma, int mu, int nu, const VectorN<Real, N>& pos) const
		{
			Real value = DeriveChristoffelSymbolSecondKind(rho, sigma, nu, mu, pos)
				- DeriveChristoffelSymbolSecondKind(rho, sigma, mu, nu, pos);

			for (int lambda = 0; lambda < N; lambda++)
			{
				value += GetChristoffelSymbolSecondKind(rho, lambda, mu, pos) * GetChristoffelSymbolSecondKind(lambda, sigma, nu, pos);
				value -= GetChristoffelSymbolSecondKind(rho, lambda, nu, pos) * GetChristoffelSymbolSecondKind(lambda, sigma, mu, pos);
			}

			return value;
		}

		/// @brief Get Riemann curvature tensor R^rho_sigma_mu_nu at a point.
		/// @param pos Position in coordinate space
		/// @return Rank-4 tensor with index variance (1 contravariant, 3 covariant)
		Tensor4<N> GetRiemannCurvatureTensor(const VectorN<Real, N>& pos) const
		{
			Tensor4<N> riemann(3, 1);

			for (int rho = 0; rho < N; rho++)
				for (int sigma = 0; sigma < N; sigma++)
					for (int mu = 0; mu < N; mu++)
						for (int nu = 0; nu < N; nu++)
							riemann(rho, sigma, mu, nu) = GetRiemannCurvatureTensor(rho, sigma, mu, nu, pos);

			return riemann;
		}

		/// @brief Get Ricci tensor component R_sigma_nu = R^rho_sigma_rho_nu.
		/// @param sigma,nu Covariant indices
		/// @param pos Position in coordinate space
		/// @return Ricci tensor component R_sigma_nu
		Real GetRicciTensor(int sigma, int nu, const VectorN<Real, N>& pos) const
		{
			Real value = REAL(0.0);
			for (int rho = 0; rho < N; rho++)
				value += GetRiemannCurvatureTensor(rho, sigma, rho, nu, pos);
			return value;
		}

		/// @brief Get Ricci tensor R_sigma_nu at a point.
		/// @param pos Position in coordinate space
		/// @return Rank-2 covariant Ricci tensor
		Tensor2<N> GetRicciTensor(const VectorN<Real, N>& pos) const
		{
			Tensor2<N> ricci(2, 0);

			for (int sigma = 0; sigma < N; sigma++)
				for (int nu = 0; nu < N; nu++)
					ricci(sigma, nu) = GetRicciTensor(sigma, nu, pos);

			return ricci;
		}

		/// @brief Get scalar curvature R = g^sigma_nu R_sigma_nu.
		/// @param pos Position in coordinate space
		/// @return Ricci scalar curvature
		Real GetRicciScalar(const VectorN<Real, N>& pos) const
		{
			MatrixNM<Real, N, N> gContravar = GetContravariantMetric(pos);
			Tensor2<N> ricci = GetRicciTensor(pos);
			Real value = REAL(0.0);

			for (int sigma = 0; sigma < N; sigma++)
				for (int nu = 0; nu < N; nu++)
					value += gContravar(sigma, nu) * ricci(sigma, nu);

			return value;
		}

		/// @brief Get Einstein tensor component G_mu_nu = R_mu_nu - 1/2 g_mu_nu R.
		/// @param mu,nu Covariant indices
		/// @param pos Position in coordinate space
		/// @return Einstein tensor component G_mu_nu
		Real GetEinsteinTensor(int mu, int nu, const VectorN<Real, N>& pos) const
		{
			return GetRicciTensor(mu, nu, pos) - REAL(0.5) * this->Component(mu, nu, pos) * GetRicciScalar(pos);
		}

		/// @brief Get Einstein tensor G_mu_nu at a point.
		/// @param pos Position in coordinate space
		/// @return Rank-2 covariant Einstein tensor
		Tensor2<N> GetEinsteinTensor(const VectorN<Real, N>& pos) const
		{
			Tensor2<N> einstein(2, 0);

			for (int mu = 0; mu < N; mu++)
				for (int nu = 0; nu < N; nu++)
					einstein(mu, nu) = GetEinsteinTensor(mu, nu, pos);

			return einstein;
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

	private:
		Real DeriveChristoffelSymbolSecondKind(int i, int j, int k, int derivIndex, const VectorN<Real, N>& pos) const
		{
			VectorN<Real, N> x = pos;
			Real original = pos[derivIndex];
			Real h = Derivation::ScaleStep(Derivation::NDer4_h, original);

			x[derivIndex] = original + h;
			Real yh = GetChristoffelSymbolSecondKind(i, j, k, x);

			x[derivIndex] = original - h;
			Real ymh = GetChristoffelSymbolSecondKind(i, j, k, x);

			x[derivIndex] = original + 2 * h;
			Real y2h = GetChristoffelSymbolSecondKind(i, j, k, x);

			x[derivIndex] = original - 2 * h;
			Real ym2h = GetChristoffelSymbolSecondKind(i, j, k, x);

			return (ym2h - y2h + 8 * (yh - ymh)) / (12 * h);
		}
	};

	/// @brief Base class for Lorentzian metric fields using signature (-,+,+,...).
	/// @details Provides causal classification from ds² = g_ij dx^i dx^j.
	template<int N>
	class LorentzianMetric : public MetricTensorField<N>
	{
	public:
		enum class IntervalType
		{
			Timelike,
			Spacelike,
			Null
		};

		LorentzianMetric() : MetricTensorField<N>(0, 2) { }
		LorentzianMetric(int numContra, int numCo) : MetricTensorField<N>(numContra, numCo) { }

		virtual typename MetricTensorField<N>::SignatureType Signature() const override
		{
			return MetricTensorField<N>::SignatureType::Lorentzian;
		}

		Real IntervalSquared(const VectorN<Real, N>& displacement, const VectorN<Real, N>& pos) const
		{
			MatrixNM<Real, N, N> g = this->GetCovariantMetric(pos);
			Real ds2 = REAL(0.0);
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					ds2 += g(i, j) * displacement[i] * displacement[j];
			return ds2;
		}

		IntervalType ClassifyInterval(const VectorN<Real, N>& displacement, const VectorN<Real, N>& pos,
			Real tolerance = REAL(1e-12)) const
		{
			Real ds2 = IntervalSquared(displacement, pos);
			if (std::abs(ds2) <= tolerance)
				return IntervalType::Null;
			return ds2 < REAL(0.0) ? IntervalType::Timelike : IntervalType::Spacelike;
		}

		bool IsTimelike(const VectorN<Real, N>& displacement, const VectorN<Real, N>& pos,
			Real tolerance = REAL(1e-12)) const
		{
			return ClassifyInterval(displacement, pos, tolerance) == IntervalType::Timelike;
		}

		bool IsSpacelike(const VectorN<Real, N>& displacement, const VectorN<Real, N>& pos,
			Real tolerance = REAL(1e-12)) const
		{
			return ClassifyInterval(displacement, pos, tolerance) == IntervalType::Spacelike;
		}

		bool IsNull(const VectorN<Real, N>& displacement, const VectorN<Real, N>& pos,
			Real tolerance = REAL(1e-12)) const
		{
			return ClassifyInterval(displacement, pos, tolerance) == IntervalType::Null;
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
		const CoordTransf<VectorFrom, VectorTo, N>& _coordTransf;

	public:
		explicit MetricTensorFromCoordTransf(const CoordTransf<VectorFrom, VectorTo, N>& inTransf) : _coordTransf(inTransf)
		{ }

		Real Component(int i, int j, const VectorN<Real, N>& pos) const override
		{
			const auto jac = _coordTransf.jacobian(pos);
			Real g_ij = 0.0;
			for (int k = 0; k < N; k++)
				g_ij += jac(k, i) * jac(k, j);
			return g_ij;
		}

		MatrixNM<Real, N, N> GetCovariantMetric(const VectorN<Real, N>& pos) const override
		{
			const auto jac = _coordTransf.jacobian(pos);
			MatrixNM<Real, N, N> metric;

			for (int i = 0; i < N; ++i)
				for (int j = 0; j < N; ++j)
				{
					metric(i, j) = REAL(0.0);
					for (int k = 0; k < N; ++k)
						metric(i, j) += jac(k, i) * jac(k, j);
				}

			return metric;
		}
	};

	/// @brief Minkowski metric for special relativity: η=diag(-1,1,1,1), signature (−,+,+,+)
	/// @note Flat spacetime with coords (ct,x,y,z), ds²=-c²dt²+dx²+dy²+dz²
	class MetricTensorMinkowski : public LorentzianMetric<4>
	{
	public:
		MetricTensorMinkowski() : LorentzianMetric<4>(0, 2) {}

		virtual Real Component(int i, int j, const VectorN<Real, 4>& pos) const override
		{
			if (i != j)
				return REAL(0.0);

			return i == 0 ? -REAL(1.0) : REAL(1.0);
		}
	};
}
#endif