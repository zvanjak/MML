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

#include <mml/MMLBase.h>
#include <algorithm>
#include <functional>
#include <limits>
#include <mml/interfaces/IDiscreteMap.h>
#include <mml/systems/DynamicalSystem/LyapunovQR.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>

#include <cmath>
#include <vector>
#include <string>

namespace MML::Systems 
{
	//=============================================================================
	// LOGISTIC MAP
	//=============================================================================

	/// @brief Logistic map: x_{n+1} = r * x_n * (1 - x_n)
	///
	/// Classic 1D map exhibiting period-doubling route to chaos.
	/// Chaotic for r ≈ 3.57 and above.
	class LogisticMap : public IDiscreteMap<1> {
	public:
		LogisticMap(Real r = 3.9) : _r(r) {
			DynSysDetail::RequireFinite(r, "LogisticMap parameter");
		}

		void map(const Vector<Real>& x, Vector<Real>& xNext) const override { xNext[0] = _r * x[0] * (1 - x[0]); }

		void jacobian(const Vector<Real>& x, Matrix<Real>& J) const override {
			J.resize(1, 1);
			J(0, 0) = _r * (1 - 2 * x[0]);
		}

		int  getNumParam() const override { return 1; }
		Real getParam(int index) const override {
			if (index != 0) throw IndexError("LogisticMap::getParam - parameter index out of range");
			return _r;
		}
		void setParam(int index, Real value) override {
			if (index != 0) throw IndexError("LogisticMap::setParam - parameter index out of range");
			DynSysDetail::RequireFinite(value, "LogisticMap parameter");
			_r = value;
		}
		std::string getParamName(int index) const override {
			if (index != 0) throw IndexError("LogisticMap::getParamName - parameter index out of range");
			return "r";
		}
		std::string getStateName(int index) const override {
			if (index != 0) throw IndexError("LogisticMap::getStateName - state index out of range");
			return "x";
		}

		/// @brief Analytical Lyapunov exponent for r=4: λ = ln(2)
		Real analyticalLyapunov() const { return (_r == MML::Real(4.0)) ? std::log(MML::Real(2.0)) : MML::Real(0.0); }

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
		HenonMap(Real a = 1.4, Real b = 0.3) : _a(a), _b(b) {
			DynSysDetail::RequireFinite(a, "HenonMap parameter a");
			DynSysDetail::RequireFinite(b, "HenonMap parameter b");
		}

		void map(const Vector<Real>& x, Vector<Real>& xNext) const override {
			xNext[0] = 1 - _a * x[0] * x[0] + x[1];
			xNext[1] = _b * x[0];
		}

		void jacobian(const Vector<Real>& x, Matrix<Real>& J) const override {
			J.resize(2, 2);
			J(0, 0) = -2 * _a * x[0];
			J(0, 1) = 1;
			J(1, 0) = _b;
			J(1, 1) = 0;
		}

		int  getNumParam() const override { return 2; }
		Real getParam(int index) const override {
			if (index < 0 || index > 1) throw IndexError("HenonMap::getParam - parameter index out of range");
			return (index == 0) ? _a : _b;
		}
		void setParam(int index, Real value) override {
			if (index < 0 || index > 1) throw IndexError("HenonMap::setParam - parameter index out of range");
			DynSysDetail::RequireFinite(value, "HenonMap parameter");
			(index == 0 ? _a : _b) = value;
		}
		void setParams(const Vector<Real>& values) override {
			if (values.size() != 2)
				throw VectorDimensionError("HenonMap::setParams - parameter vector size must equal 2", values.size(), 2);
			DynSysDetail::RequireFinite(values[0], "HenonMap parameter a");
			DynSysDetail::RequireFinite(values[1], "HenonMap parameter b");
			_a = values[0];
			_b = values[1];
		}
		std::string getParamName(int index) const override {
			if (index < 0 || index > 1) throw IndexError("HenonMap::getParamName - parameter index out of range");
			return (index == 0) ? "a" : "b";
		}
		std::string getStateName(int index) const override {
			if (index < 0 || index > 1) throw IndexError("HenonMap::getStateName - state index out of range");
			return (index == 0) ? "x" : "y";
		}

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
		StandardMap(Real K = 0.971635) : _kickStrength(K) {
			DynSysDetail::RequireFinite(K, "StandardMap parameter");
		}

		void map(const Vector<Real>& x, Vector<Real>& xNext) const override {
			Real p = x[0] + _kickStrength * std::sin(x[1]);
			Real theta = x[1] + p;

			// Wrap to [0, 2π)
			const Real twopi = 2 * Constants::PI;
			p = p - twopi * std::floor(p / twopi);
			theta = theta - twopi * std::floor(theta / twopi);

			xNext[0] = p;
			xNext[1] = theta;
		}

		void jacobian(const Vector<Real>& x, Matrix<Real>& J) const override {
			J.resize(2, 2);
			Real c = _kickStrength * std::cos(x[1]);
			J(0, 0) = 1;
			J(0, 1) = c;
			J(1, 0) = 1;
			J(1, 1) = 1 + c;
		}

		int  getNumParam() const override { return 1; }
		Real getParam(int index) const override {
			if (index != 0) throw IndexError("StandardMap::getParam - parameter index out of range");
			return _kickStrength;
		}
		void setParam(int index, Real value) override {
			if (index != 0) throw IndexError("StandardMap::setParam - parameter index out of range");
			DynSysDetail::RequireFinite(value, "StandardMap parameter");
			_kickStrength = value;
		}
		std::string getParamName(int index) const override {
			if (index != 0) throw IndexError("StandardMap::getParamName - parameter index out of range");
			return "K";
		}
		std::string getStateName(int index) const override {
			if (index < 0 || index > 1) throw IndexError("StandardMap::getStateName - state index out of range");
			return (index == 0) ? "p" : "theta";
		}

		/// @brief Standard map preserves area (det(J) = 1)
		bool isAreaPreserving() const { return true; }

	private:
		Real _kickStrength;
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
		TentMap(Real mu = 2.0) : _mu(mu) {
			DynSysDetail::RequireFinite(mu, "TentMap parameter");
		}

		void map(const Vector<Real>& x, Vector<Real>& xNext) const override { xNext[0] = _mu * std::min(x[0], 1 - x[0]); }

		void jacobian(const Vector<Real>& x, Matrix<Real>& J) const override {
			J.resize(1, 1);
			// Derivative is ±μ depending on which side of peak
			J(0, 0) = (x[0] < 0.5) ? _mu : -_mu;
		}

		int  getNumParam() const override { return 1; }
		Real getParam(int index) const override {
			if (index != 0) throw IndexError("TentMap::getParam - parameter index out of range");
			return _mu;
		}
		void setParam(int index, Real value) override {
			if (index != 0) throw IndexError("TentMap::setParam - parameter index out of range");
			DynSysDetail::RequireFinite(value, "TentMap parameter");
			_mu = value;
		}
		std::string getParamName(int index) const override {
			if (index != 0) throw IndexError("TentMap::getParamName - parameter index out of range");
			return "mu";
		}
		std::string getStateName(int index) const override {
			if (index != 0) throw IndexError("TentMap::getStateName - state index out of range");
			return "x";
		}

		/// @brief Analytical Lyapunov exponent for μ=2: λ = ln(2)
		Real analyticalLyapunov() const { return (_mu == MML::Real(2.0)) ? std::log(MML::Real(2.0)) : std::log(_mu); }

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
		bool isChaotic;          ///< Positive-growth heuristic using Precision::DiscreteLyapunovChaosThreshold
	};

	//=============================================================================
	// DISCRETE MAP LYAPUNOV ANALYZER
	//=============================================================================

	/// @brief Compute Lyapunov exponent for discrete maps
	/// @tparam N State dimension
	template <int N>
	class DiscreteMapLyapunov {
		static_assert(N > 0, "Discrete Lyapunov analysis requires a positive state dimension");
	public:
		/// @brief Compute Lyapunov exponents
		/// @param map The discrete map
		/// @param x0 Initial condition
		/// @param numIterations Total iterations
		/// @param skipTransient Iterations to skip
		/// @return DiscreteMapLyapunovResult with exponents and chaos indicator
		/// @throws SingularMatrixError for zero or numerically dependent tangent columns
		/// @throws NumericalMethodError for non-finite callback outputs or tangent frames
		static DiscreteMapLyapunovResult<Real> Compute(IDiscreteMap<N>& map, const Vector<Real>& x0, int numIterations = 10000,
																int skipTransient = 1000) {
			if (numIterations <= 0 || skipTransient < 0)
				throw ArgumentError("Discrete Lyapunov: require positive measurement count and nonnegative transient count");
			DynSysDetail::CheckState<ArgumentError>(x0, N, "Discrete Lyapunov initial state");
			Vector<Real> x = x0;
			Vector<Real> xNext(N);
			Matrix<Real> J(N, N);
			Matrix<Real> Q = Matrix<Real>::Identity(N);
			Vector<Real> lyapunovSums(N, 0.0);

			// Skip transient
			for (int i = 0; i < skipTransient; ++i) {
				xNext = Vector<Real>(N, std::numeric_limits<Real>::quiet_NaN());
				map.map(x, xNext);
				DynSysDetail::CheckState(xNext, N, "Discrete Lyapunov map output");
				x = xNext;
			}

			// Accumulate stretching
			for (int iter = 0; iter < numIterations; ++iter) {
				// Get Jacobian
				J = Matrix<Real>(N, N, std::numeric_limits<Real>::quiet_NaN());
				map.jacobian(x, J);
				DynSysDetail::CheckJacobian(J, N);

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

				const auto logarithms = LyapunovDetail::Orthonormalize(Q);
				for (int j = 0; j < N; ++j)
					lyapunovSums[j] += logarithms[j];

				// Iterate map
				xNext = Vector<Real>(N, std::numeric_limits<Real>::quiet_NaN());
				map.map(x, xNext);
				DynSysDetail::CheckState(xNext, N, "Discrete Lyapunov map output");
				x = xNext;
			}

			// Build result
			DiscreteMapLyapunovResult<Real> result;
			result.exponents.resize(N);
			for (int i = 0; i < N; ++i)
				result.exponents[i] = lyapunovSums[i] / numIterations;
			std::sort(result.exponents.begin(), result.exponents.end(), std::greater<Real>());

			result.maxExponent = result.exponents[0];
			result.isChaotic = result.maxExponent > Precision::DiscreteLyapunovChaosThreshold;

			return result;
		}

};

} // namespace MML::Systems
#endif // MML_DISCRETE_MAPS_H
