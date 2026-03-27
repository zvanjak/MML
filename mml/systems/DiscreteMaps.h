///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DiscreteMaps.h                                                      ///
///  Description: Discrete dynamical systems (iterated maps)                          ///
///               Logistic, Hénon, Standard, Tent maps + Lyapunov analysis            ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                        ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_DISCRETE_MAPS_H
#define MML_DISCRETE_MAPS_H

#include "MMLBase.h"
#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"

#include <cmath>
#include <vector>
#include <string>

namespace MML::Systems 
{
	//=============================================================================
	// DISCRETE MAP INTERFACE
	//=============================================================================

	/// @brief Base interface for discrete maps x_{n+1} = f(x_n)
	/// @tparam N State dimension
	template <int N>
	class IDiscreteMap {
	public:
		virtual ~IDiscreteMap() = default;

		/// @brief Apply map once: x_{n+1} = f(x_n)
		virtual void map(const Vector<Real>& x, Vector<Real>& xNext) const = 0;

		/// @brief Apply map once (convenience version returning result)
		Vector<Real> iterate(const Vector<Real>& x) const {
			Vector<Real> xNext(N);
			map(x, xNext);
			return xNext;
		}

		/// @brief Jacobian of map Df(x)
		virtual void jacobian(const Vector<Real>& x, Matrix<Real>& J) const = 0;

		/// @brief State dimension
		int getDim() const { return N; }

		/// @brief Get parameter
		virtual Real getParam(int index) const = 0;

		/// @brief Set parameter
		virtual void setParam(int index, Real value) = 0;

		/// @brief Get parameter name
		virtual std::string getParamName(int index) const = 0;

		/// @brief Get state name
		virtual std::string getStateName(int index) const = 0;
	};

	//=============================================================================
	// LOGISTIC MAP
	//=============================================================================

	/// @brief Logistic map: x_{n+1} = r * x_n * (1 - x_n)
	///
	/// Classic 1D map exhibiting period-doubling route to chaos.
	/// Chaotic for r ≈ 3.57 and above.
	class LogisticMap : public IDiscreteMap<1> {
	public:
		LogisticMap(Real r = 3.9) : _r(r) {}

		void map(const Vector<Real>& x, Vector<Real>& xNext) const override { xNext[0] = _r * x[0] * (1 - x[0]); }

		void jacobian(const Vector<Real>& x, Matrix<Real>& J) const override {
			J.Resize(1, 1);
			J(0, 0) = _r * (1 - 2 * x[0]);
		}

		Real getParam(int index) const override { return (index == 0) ? _r : 0; }
		void setParam(int index, Real value) override {
			if (index == 0)
				_r = value;
		}
		std::string getParamName(int /*index*/) const override { return "r"; }
		std::string getStateName(int /*index*/) const override { return "x"; }

		/// @brief Analytical Lyapunov exponent for r=4: λ = ln(2)
		Real analyticalLyapunov() const { return (_r == 4.0) ? std::log(2.0) : 0.0; }

		/// @brief Generate orbit of n iterations
		std::vector<Vector<Real>> orbit(const Vector<Real>& x0, int n) const {
			std::vector<Vector<Real>> result;
			Vector<Real> x = x0;
			for (int i = 0; i < n; ++i) {
				x = iterate(x);
				result.push_back(x);
			}
			return result;
		}

	private:
		Real _r;
	};

	//=============================================================================
	// HÉNON MAP
	//=============================================================================

	/// @brief Hénon map: x_{n+1} = 1 - a*x_n² + y_n
	///                   y_{n+1} = b*x_n
	///
	/// Classic 2D map with a strange attractor.
	/// Standard parameters: a = 1.4, b = 0.3
	class HenonMap : public IDiscreteMap<2> {
	public:
		HenonMap(Real a = 1.4, Real b = 0.3) : _a(a), _b(b) {}

		void map(const Vector<Real>& x, Vector<Real>& xNext) const override {
			xNext[0] = 1 - _a * x[0] * x[0] + x[1];
			xNext[1] = _b * x[0];
		}

		void jacobian(const Vector<Real>& x, Matrix<Real>& J) const override {
			J.Resize(2, 2);
			J(0, 0) = -2 * _a * x[0];
			J(0, 1) = 1;
			J(1, 0) = _b;
			J(1, 1) = 0;
		}

		Real getParam(int index) const override { return (index == 0) ? _a : _b; }
		void setParam(int index, Real value) override {
			if (index == 0)
				_a = value;
			else if (index == 1)
				_b = value;
		}
		std::string getParamName(int index) const override { return (index == 0) ? "a" : "b"; }
		std::string getStateName(int index) const override { return (index == 0) ? "x" : "y"; }

		/// @brief Get parameter a
		Real getA() const { return _a; }

		/// @brief Get parameter b
		Real getB() const { return _b; }

		/// @brief Jacobian determinant is constant: det(J) = -b
		Real jacobianDeterminant() const { return -_b; }

		/// @brief Generate orbit of n iterations
		std::vector<Vector<Real>> orbit(const Vector<Real>& x0, int n) const {
			std::vector<Vector<Real>> result;
			Vector<Real> x = x0;
			for (int i = 0; i < n; ++i) {
				x = iterate(x);
				result.push_back(x);
			}
			return result;
		}

	private:
		Real _a, _b;
	};

	//=============================================================================
	// STANDARD (CHIRIKOV) MAP
	//=============================================================================

	/// @brief Standard (Chirikov-Taylor) map
	///
	/// p_{n+1} = p_n + K*sin(θ_n)    (mod 2π)
	/// θ_{n+1} = θ_n + p_{n+1}       (mod 2π)
	///
	/// Canonical model for Hamiltonian chaos.
	class StandardMap : public IDiscreteMap<2> {
	public:
		StandardMap(Real K = 0.971635) : _K(K) {}

		void map(const Vector<Real>& x, Vector<Real>& xNext) const override {
			Real p = x[0] + _K * std::sin(x[1]);
			Real theta = x[1] + p;

			// Wrap to [0, 2π)
			const Real twopi = 2 * Constants::PI;
			p = p - twopi * std::floor(p / twopi);
			theta = theta - twopi * std::floor(theta / twopi);

			xNext[0] = p;
			xNext[1] = theta;
		}

		void jacobian(const Vector<Real>& x, Matrix<Real>& J) const override {
			J.Resize(2, 2);
			Real c = _K * std::cos(x[1]);
			J(0, 0) = 1;
			J(0, 1) = c;
			J(1, 0) = 1;
			J(1, 1) = 1 + c;
		}

		Real getParam(int /*index*/) const override { return _K; }
		void setParam(int /*index*/, Real value) override { _K = value; }
		std::string getParamName(int /*index*/) const override { return "K"; }
		std::string getStateName(int index) const override { return (index == 0) ? "p" : "theta"; }

		/// @brief Standard map preserves area (det(J) = 1)
		bool isAreaPreserving() const { return true; }

	private:
		Real _K;
	};

	//=============================================================================
	// TENT MAP
	//=============================================================================

	/// @brief Tent map
	///
	/// x_{n+1} = μ * min(x_n, 1 - x_n)
	///
	/// Piecewise linear, exhibits chaos for μ = 2.
	class TentMap : public IDiscreteMap<1> {
	public:
		TentMap(Real mu = 2.0) : _mu(mu) {}

		void map(const Vector<Real>& x, Vector<Real>& xNext) const override { xNext[0] = _mu * std::min(x[0], 1 - x[0]); }

		void jacobian(const Vector<Real>& x, Matrix<Real>& J) const override {
			J.Resize(1, 1);
			// Derivative is ±μ depending on which side of peak
			J(0, 0) = (x[0] < 0.5) ? _mu : -_mu;
		}

		Real getParam(int /*index*/) const override { return _mu; }
		void setParam(int /*index*/, Real value) override { _mu = value; }
		std::string getParamName(int /*index*/) const override { return "mu"; }
		std::string getStateName(int /*index*/) const override { return "x"; }

		/// @brief Analytical Lyapunov exponent for μ=2: λ = ln(2)
		Real analyticalLyapunov() const { return (_mu == 2.0) ? std::log(2.0) : std::log(_mu); }

	private:
		Real _mu;
	};

	//=============================================================================
	// DISCRETE MAP LYAPUNOV RESULT
	//=============================================================================

	/// @brief Result of discrete map Lyapunov computation
	template<typename Type = Real>
	struct DiscreteMapLyapunovResult {
		Vector<Type> exponents;  ///< Lyapunov exponents (descending order)
		Type maxExponent;        ///< Largest exponent (λ₁)
		bool isChaotic;          ///< True if λ₁ > 0
	};

	//=============================================================================
	// DISCRETE MAP LYAPUNOV ANALYZER
	//=============================================================================

	/// @brief Compute Lyapunov exponent for discrete maps
	/// @tparam N State dimension
	template <int N>
	class DiscreteMapLyapunov {
	public:
		/// @brief Compute Lyapunov exponents
		/// @param map The discrete map
		/// @param x0 Initial condition
		/// @param numIterations Total iterations
		/// @param skipTransient Iterations to skip
		/// @return DiscreteMapLyapunovResult with exponents and chaos indicator
		static DiscreteMapLyapunovResult<Real> Compute(IDiscreteMap<N>& map, const Vector<Real>& x0, int numIterations = 10000,
																int skipTransient = 1000) {
			Vector<Real> x = x0;
			Vector<Real> xNext(N);
			Matrix<Real> J(N, N);
			Matrix<Real> Q = Matrix<Real>::Identity(N);
			Vector<Real> lyapunovSums(N, 0.0);

			// Skip transient
			for (int i = 0; i < skipTransient; ++i) {
				map.map(x, xNext);
				x = xNext;
			}

			// Accumulate stretching
			for (int iter = 0; iter < numIterations; ++iter) {
				// Get Jacobian
				map.jacobian(x, J);

				// Multiply Q by J
				Matrix<Real> JQ(N, N);
				for (int i = 0; i < N; ++i) {
					for (int j = 0; j < N; ++j) {
						Real sum = 0;
						for (int k = 0; k < N; ++k)
							sum += J(i, k) * Q(k, j);
						JQ(i, j) = sum;
					}
				}
				Q = JQ;

				// Gram-Schmidt orthonormalization
				for (int j = 0; j < N; ++j) {
					// Get column j
					Vector<Real> v(N);
					for (int i = 0; i < N; ++i)
						v[i] = Q(i, j);

					// Orthogonalize
					for (int k = 0; k < j; ++k) {
						Real dot = 0;
						for (int i = 0; i < N; ++i)
							dot += v[i] * Q(i, k);
						for (int i = 0; i < N; ++i)
							v[i] -= dot * Q(i, k);
					}

					// Normalize
					Real norm = 0;
					for (int i = 0; i < N; ++i)
						norm += v[i] * v[i];
					norm = std::sqrt(norm);

					if (norm > Precision::DivisionSafetyThreshold) {
						lyapunovSums[j] += std::log(norm);
						for (int i = 0; i < N; ++i)
							Q(i, j) = v[i] / norm;
					}
				}

				// Iterate map
				map.map(x, xNext);
				x = xNext;
			}

			// Build result
			DiscreteMapLyapunovResult<Real> result;
			result.exponents.Resize(N);
			for (int i = 0; i < N; ++i)
				result.exponents[i] = lyapunovSums[i] / numIterations;

			result.maxExponent = result.exponents[0];
			result.isChaotic = result.maxExponent > 0.01;  // Small threshold

			return result;
		}

};

} // namespace MML::Systems
#endif // MML_DISCRETE_MAPS_H
