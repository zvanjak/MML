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

#include <mml/algorithms/MatrixSingularity.h>
#include <mml/MMLSingularityHandling.h>

#include <mml/interfaces/IFunction.h>

#include <mml/base/Vector/VectorN.h>
#include <mml/base/Matrix/MatrixNM.h>

#include <mml/core/Derivation.h>
#include <mml/core/Derivation/DerivativeProviders.h>
#include <mml/core/Derivation/DerivationTensorField.h>
#include <mml/core/CoordTransf/CoordTransfBase.h>

#include <cmath>
#include <typeinfo>

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
		MatrixNM<Real, N, N> ReadStoredMetric(const VectorN<Real, N>& pos) const {
			if (!((this->getNumCovar() == 2 && this->getNumContravar() == 0) ||
			      (this->getNumCovar() == 0 && this->getNumContravar() == 2)))
				throw DomainError("Metric tensor requires fully covariant or fully contravariant storage");
			Singularity::RequireFiniteCoordinates<N>(pos, "Metric tensor");
			MatrixNM<Real, N, N> g;
			for (int i = 0; i < N; ++i)
				for (int j = 0; j < N; ++j)
					g(i,j) = this->component(i,j,pos);
			Singularity::RequireSymmetricMatrix(g, "Metric tensor");
			return g;
		}

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
					ret(i, j) = this->component(i, j, pos);

			return ret;
		}

		/// @brief Get covariant metric tensor components gᵢⱼ at a point
		/// @param pos Position in coordinate space
		/// @return N×N matrix of covariant metric components
		/// @note Measures squared arc length: ds² = gᵢⱼ dxⁱ dxʲ
		// Get the covariant metric tensor components at a point (gᵢⱼ)
		virtual MatrixNM<Real, N, N> GetCovariantMetric(const VectorN<Real, N>& pos) const
		{
			const auto stored = ReadStoredMetric(pos);
			return this->getNumContravar() == 2
				? Singularity::InverseMatrix(stored, Singularity::DEFAULT_POLICY) : stored;
		}

		/// @brief Get contravariant metric tensor components gⁱʲ (inverse of gᵢⱼ)
		/// @param pos Position in coordinate space
		/// @return N×N matrix of contravariant metric components
		/// @note Used for raising indices: vⁱ = gⁱʲ vⱼ
		// Get the contravariant metric tensor components at a point (gⁱʲ = inverse of gᵢⱼ)
		MatrixNM<Real, N, N> GetContravariantMetric(const VectorN<Real, N>& pos) const
		{
			return GetContravariantMetric(pos, Singularity::DEFAULT_POLICY);
		}

		/// Explicit sentinel policy for a singular metric; Clamp does not regularize it.
		MatrixNM<Real, N, N> GetContravariantMetric(const VectorN<Real, N>& pos, SingularityPolicy policy) const {
			Singularity::RequireMatrixPolicy(policy);
			Singularity::RequireFiniteCoordinates<N>(pos, "Metric inversion");
			if (this->getNumContravar() == 2)
				return ReadStoredMetric(pos);
			return Singularity::InverseMatrix(GetCovariantMetric(pos), policy);
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
		/// @param t Rank-2 tensor whose selected slot is covariant
		/// @param index_to_raise Which index to raise (0 or 1)
		/// @param pos Position where metric is evaluated
		/// @return Tensor with only the selected slot raised
		Tensor2<N> RaiseIndex(const Tensor2<N>& t, int index_to_raise, const VectorN<Real, N>& pos) const
		{
			if (index_to_raise < 0 || index_to_raise > 1)
				throw ArgumentError("MetricTensorField::RaiseIndex - slot must be 0 or 1");
			const auto components = GetContravariantMetric(pos);
			Tensor2<N> inverse(0, 2);
			for (int i = 0; i < N; ++i)
				for (int j = 0; j < N; ++j)
					inverse(i,j) = components(i,j);
			return index_to_raise == 0 ? t.RaiseFirstIndex(inverse) : t.RaiseSecondIndex(inverse);
		}

		/// @brief Lower one index of a rank-2 tensor: T_ij = g_ik T^k_j
		/// @param t Rank-2 tensor with at least one contravariant index
		/// @param index_to_lower Which index to lower (0 or 1)
		/// @param pos Position where metric is evaluated
		/// @return Tensor with one lowered index
		Tensor2<N> LowerIndex(const Tensor2<N>& t, int index_to_lower, const VectorN<Real, N>& pos) const
		{
			if (index_to_lower < 0 || index_to_lower > 1)
				throw ArgumentError("MetricTensorField::LowerIndex - slot must be 0 or 1");
			const auto components = GetCovariantMetric(pos);
			Tensor2<N> metric(2, 0);
			for (int i = 0; i < N; ++i)
				for (int j = 0; j < N; ++j)
					metric(i,j) = components(i,j);
			return index_to_lower == 0 ? t.LowerFirstIndex(metric) : t.LowerSecondIndex(metric);
		}

		/// @brief Get Christoffel symbol of the first kind Γᵢⱼₖ (all indices lowered)
		/// @param i,j,k Indices
		/// @param pos Position in coordinate space
		/// @return Γᵢⱼₖ = gₘₖ Γᵐᵢⱼ (related to second kind via metric)
		Real GetChristoffelSymbolFirstKind(int i, int j, int k, const VectorN<Real, N>& pos) const
		{
			ValidateChristoffelArguments(i, j, k, pos);
			const MetricTensorField<N>& g = *this;

			Real gamma_ijk = 0.0;
			for (int m = 0; m < N; m++)
			{
				gamma_ijk += g.component(m, k, pos) * GetChristoffelSymbolSecondKind(m, i, j, pos);
			}
			return gamma_ijk;
		}
		/// @brief Get Christoffel symbol of the second kind Γᵐᵢⱼ (one index raised)
		/// @param i,j,k Indices (i=contravariant, j,k=covariant)
		/// @param pos Position in coordinate space
		/// @return Γᵐᵢⱼ = ½ gᵐˡ (∂ⱼ gˡₖ + ∂ₖ gˡⱼ - ∂ˡ gⱼₖ) (connection coefficients)
		/// @note Used in geodesic equation and covariant derivatives
		virtual Real GetChristoffelSymbolSecondKind(int i, int j, int k, const VectorN<Real, N>& pos) const {
			ValidateChristoffelArguments(i, j, k, pos);
			if (const auto* provider = dynamic_cast<const IMetricJetProvider<N>*>(this))
				return ConnectionFromJet(provider->EvaluateMetricJet(pos), i, j, k);
			return GetChristoffelSymbolSecondKindNumerical(i, j, k, pos);
		}

		/// @brief Fourth-order numerical connection reference, bypassing analytical overrides.
		/// @note Requires covariant metric storage. Singular matrices fail during inversion.
		Real GetChristoffelSymbolSecondKindNumerical(int i, int j, int k, const VectorN<Real, N>& pos) const {
			ValidateChristoffelArguments(i, j, k, pos);
			const MetricTensorField<N>& g = *this;

			// Γⁱⱼₖ = ½ gⁱˡ (∂ⱼgₗₖ + ∂ₖgₗⱼ - ∂ₗgⱼₖ)
			// Need contravariant metric gⁱˡ to raise the index
			const auto g_contravar = GetContravariantMetric(pos);

			Real gamma_ijk = 0.0;
			for (int l = 0; l < N; l++)
			{
				Real coef1 = Derivation::NDer4Partial<N>(g, l, k, j, pos, nullptr);  // ∂ⱼgₗₖ
				Real coef2 = Derivation::NDer4Partial<N>(g, l, j, k, pos, nullptr);  // ∂ₖgₗⱼ
				Real coef3 = Derivation::NDer4Partial<N>(g, j, k, l, pos, nullptr);  // ∂ₗgⱼₖ

				gamma_ijk += MML::Real(0.5) * g_contravar[i][l] * (coef1 + coef2 - coef3);
			}
			if (!std::isfinite(gamma_ijk))
				throw DomainError("Numerical Christoffel calculation produced a non-finite coefficient");
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
			Tensor4<N> riemann(4, 0);
			riemann.setContravar(0, true);

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
			Real value = static_cast<MML::Real>(0.0);
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
			Real value = static_cast<MML::Real>(0.0);

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
			return GetRicciTensor(mu, nu, pos) - static_cast<MML::Real>(0.5) * this->component(mu, nu, pos) * GetRicciScalar(pos);
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
				Real comp_val = Differential::VectorPartial<N>(func, i, j, pos);

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
			Real ret = Differential::VectorPartial<N>(func, i, j, pos);

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
				Real comp_val = Differential::VectorPartial<N>(func, i, j, pos);

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
			Real comp_val = Differential::VectorPartial<N>(func, i, j, pos);

			for (int k = 0; k < N; k++)
				comp_val -= GetChristoffelSymbolSecondKind(k, i, j, pos) * func(pos)[k];

			return comp_val;
		}

	protected:
		void ValidateChristoffelArguments(int i, int j, int k, const VectorN<Real, N>& pos) const {
			if (i < 0 || i >= N || j < 0 || j >= N || k < 0 || k >= N)
				throw IndexError("Christoffel index outside metric dimension");
			if (this->getNumContravar() != 0 || this->getNumCovar() != 2)
				throw DomainError("Christoffel symbols require covariant metric storage; use the covariant metric field");
			for (int axis = 0; axis < N; ++axis)
				if (!std::isfinite(pos[axis]))
					throw DomainError("Christoffel symbols require finite coordinates");
		}

		// d(g^-1) = -g^-1 (dg) g^-1; no differentiation through matrix inversion.
		static Real ConnectionFromJet(const MetricJet<N>& jet, int i, int j, int k, int axis = -1) {
			if (i < 0 || i >= N || j < 0 || j >= N || k < 0 || k >= N || axis < -1 || axis >= N)
				throw IndexError("Metric jet connection index outside dimension");
			Singularity::RequireFiniteMatrix(jet.value, "Metric derivative provider");
			for (int a = 0; a < N; ++a) {
				Singularity::RequireFiniteMatrix(jet.first[a], "Metric derivative provider");
				for (int b = 0; b < N; ++b)
					Singularity::RequireFiniteMatrix(jet.second[a][b], "Metric derivative provider");
			}
			const auto inv = Singularity::InverseMatrix(jet.value, Singularity::DEFAULT_POLICY);
			Real result = 0;
			for (int l = 0; l < N; ++l) {
				const Real coefficient = jet.first[j](l, k) + jet.first[k](l, j) - jet.first[l](j, k);
				if (axis < 0) {
					result += Real(0.5) * inv(i, l) * coefficient;
				} else {
					Real inverseDerivative = 0;
					for (int a = 0; a < N; ++a)
						for (int b = 0; b < N; ++b)
							inverseDerivative -= inv(i, a) * jet.first[axis](a, b) * inv(b, l);
					result += Real(0.5) * (inverseDerivative * coefficient +
																 inv(i, l) * (jet.second[axis][j](l, k) + jet.second[axis][k](l, j) - jet.second[axis][l](j, k)));
				}
			}
			if (!std::isfinite(result))
				throw NumericalMethodError("Metric derivative provider produced a non-finite connection derivative");
			return result;
		}

		virtual Real DeriveChristoffelSymbolSecondKind(int i, int j, int k, int derivIndex, const VectorN<Real, N>& pos) const {
			ValidateChristoffelArguments(i, j, k, pos);
			if (derivIndex < 0 || derivIndex >= N)
				throw IndexError("Connection derivative index outside metric dimension");
			if (const auto* provider = dynamic_cast<const IMetricJetProvider<N>*>(this))
				return ConnectionFromJet(provider->EvaluateMetricJet(pos), i, j, k, derivIndex);
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
			Singularity::RequireFiniteCoordinates<N>(displacement, "Lorentzian interval displacement");
			MatrixNM<Real, N, N> g = this->GetCovariantMetric(pos);
			Real ds2 = static_cast<MML::Real>(0.0);
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					ds2 += g(i, j) * displacement[i] * displacement[j];
			if (!std::isfinite(ds2))
				throw NumericalMethodError("Lorentzian interval is not representable");
			return ds2;
		}

		IntervalType ClassifyInterval(const VectorN<Real, N>& displacement, const VectorN<Real, N>& pos,
			Real tolerance = Precision::ZeroTolerance) const
		{
			ValidateNonNegativeTolerance(tolerance, "Lorentzian interval tolerance");
			Real ds2 = IntervalSquared(displacement, pos);
			if (std::abs(ds2) <= tolerance)
				return IntervalType::Null;
			return ds2 < static_cast<MML::Real>(0.0) ? IntervalType::Timelike : IntervalType::Spacelike;
		}

		bool IsTimelike(const VectorN<Real, N>& displacement, const VectorN<Real, N>& pos,
			Real tolerance = Precision::ZeroTolerance) const
		{
			return ClassifyInterval(displacement, pos, tolerance) == IntervalType::Timelike;
		}

		bool IsSpacelike(const VectorN<Real, N>& displacement, const VectorN<Real, N>& pos,
			Real tolerance = Precision::ZeroTolerance) const
		{
			return ClassifyInterval(displacement, pos, tolerance) == IntervalType::Spacelike;
		}

		bool IsNull(const VectorN<Real, N>& displacement, const VectorN<Real, N>& pos,
			Real tolerance = Precision::ZeroTolerance) const
		{
			return ClassifyInterval(displacement, pos, tolerance) == IntervalType::Null;
		}
	};

	/// @brief Flat metric tensor for Cartesian 3D (gᵢⱼ = δᵢⱼ, all Christoffel symbols vanish)
	class MetricTensorCartesian3D : public MetricTensorField<3>
	{
	public:
		MetricTensorCartesian3D() : MetricTensorField<3>(0, 2) { }

		Real GetChristoffelSymbolSecondKind(int i, int j, int k, const VectorN<Real, 3>& pos) const override {
			// A subclass may change the metric without supplying its own connection.
			if (typeid(*this) != typeid(MetricTensorCartesian3D))
				return GetChristoffelSymbolSecondKindNumerical(i, j, k, pos);
			ValidateChristoffelArguments(i, j, k, pos);
			return 0.0;
		}

		Real component(int i, int j, const VectorN<Real, 3>& pos) const override
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

		Real GetChristoffelSymbolSecondKind(int i, int j, int k, const VectorN<Real, 3>& pos) const override {
			if (typeid(*this) != typeid(MetricTensorSpherical))
				return GetChristoffelSymbolSecondKindNumerical(i, j, k, pos);
			ValidateChristoffelArguments(i, j, k, pos);
			const Real r = pos[0], s = std::sin(pos[1]), c = std::cos(pos[1]);
			const Real inverseR = Singularity::SafeDivide(1.0, r, SingularityPolicy::Throw, "Spherical Christoffel symbols: origin");
			const Real cotTheta = Singularity::SafeDivide(c, s, SingularityPolicy::Throw, "Spherical Christoffel symbols: pole");
			if (i == 0 && j == 1 && k == 1)
				return -r;
			if (i == 0 && j == 2 && k == 2)
				return -r * s * s;
			if (i == 1 && ((j == 0 && k == 1) || (j == 1 && k == 0)))
				return inverseR;
			if (i == 1 && j == 2 && k == 2)
				return -s * c;
			if (i == 2 && ((j == 0 && k == 2) || (j == 2 && k == 0)))
				return inverseR;
			if (i == 2 && ((j == 1 && k == 2) || (j == 2 && k == 1)))
				return cotTheta;
			return 0.0;
		}

		virtual  Real component(int i, int j, const VectorN<Real, 3>& pos) const override
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

		virtual Real component(int i, int j, const VectorN<Real, 3>& pos) const override
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

		Real GetChristoffelSymbolSecondKind(int i, int j, int k, const VectorN<Real, 3>& pos) const override {
			if (typeid(*this) != typeid(MetricTensorCylindrical))
				return GetChristoffelSymbolSecondKindNumerical(i, j, k, pos);
			ValidateChristoffelArguments(i, j, k, pos);
			const Real inverseRho = Singularity::SafeDivide(1.0, pos[0], SingularityPolicy::Throw, "Cylindrical Christoffel symbols: axis");
			if (i == 0 && j == 1 && k == 1)
				return -pos[0];
			if (i == 1 && ((j == 0 && k == 1) || (j == 1 && k == 0)))
				return inverseRho;
			return 0.0;
		}

		virtual Real component(int i, int j, const VectorN<Real, 3>& pos) const override
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
		MetricTensorFromCoordTransf(CoordTransf<VectorFrom, VectorTo, N>&&) = delete;
		MetricTensorFromCoordTransf(const CoordTransf<VectorFrom, VectorTo, N>&&) = delete;

		Real component(int i, int j, const VectorN<Real, N>& pos) const override
		{
			if (i < 0 || i >= N || j < 0 || j >= N)
				throw ArgumentError("Induced metric - component index out of range");
			Singularity::RequireFiniteCoordinates<N>(pos, "Induced metric");
			const auto jac = _coordTransf.jacobian(pos);
			Singularity::RequireFiniteMatrix(jac, "Induced metric Jacobian");
			Real g_ij = 0.0;
			for (int k = 0; k < CoordTransf<VectorFrom, VectorTo, N>::ToDim; k++)
				g_ij += jac(k, i) * jac(k, j);
			if (!std::isfinite(g_ij))
				throw NumericalMethodError("Induced metric component is not representable");
			return g_ij;
		}

		MatrixNM<Real, N, N> GetCovariantMetric(const VectorN<Real, N>& pos) const override
		{
			Singularity::RequireFiniteCoordinates<N>(pos, "Induced metric");
			const auto jac = _coordTransf.jacobian(pos);
			Singularity::RequireFiniteMatrix(jac, "Induced metric Jacobian");
			MatrixNM<Real, N, N> metric;

			for (int i = 0; i < N; ++i)
				for (int j = 0; j < N; ++j)
				{
					metric(i, j) = static_cast<MML::Real>(0.0);
					for (int k = 0; k < CoordTransf<VectorFrom, VectorTo, N>::ToDim; ++k)
						metric(i, j) += jac(k, i) * jac(k, j);
				}

			Singularity::RequireFiniteMatrix(metric, "Induced metric");
			return metric;
		}
	};

	/// @brief Minkowski metric for special relativity: η=diag(-1,1,1,1), signature (−,+,+,+)
	/// @note Flat spacetime with coords (ct,x,y,z), ds²=-c²dt²+dx²+dy²+dz²
	class MetricTensorMinkowski : public LorentzianMetric<4>
	{
	public:
		MetricTensorMinkowski() : LorentzianMetric<4>(0, 2) {}

		Real GetChristoffelSymbolSecondKind(int i, int j, int k, const VectorN<Real, 4>& pos) const override {
			if (typeid(*this) != typeid(MetricTensorMinkowski))
				return GetChristoffelSymbolSecondKindNumerical(i, j, k, pos);
			ValidateChristoffelArguments(i, j, k, pos);
			return static_cast<MML::Real>(0.0);
		}

		virtual Real component(int i, int j, const VectorN<Real, 4>& pos) const override
		{
			if (i != j)
				return static_cast<MML::Real>(0.0);

			return i == 0 ? -static_cast<MML::Real>(1.0) : static_cast<MML::Real>(1.0);
		}
	};
}
#endif
