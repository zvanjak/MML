///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        CurveFitting.h                                                      ///
///  Description: Curve fitting algorithms (linear, polynomial, exponential, spline)  ///
///               Least squares fitting and regression utilities                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CURVE_FITTING_H
#define MML_CURVE_FITTING_H

#include <mml/MMLBase.h>
#include <mml/MMLExceptions.h>
#include <algorithm>

#include <mml/interfaces/IFunction.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/AlgorithmTypes.h>
#include <mml/core/LinAlgEqSolvers.h>
#include <mml/core/OrthogonalBasis/ChebyshevBasis.h>
#include <mml/core/OrthogonalBasis/LegendreBasis.h>
#include <mml/core/OrthogonalBasis/HermiteBasis.h>
#include <mml/core/OrthogonalBasis/LaguerreBasis.h>

#include <functional>
#include <type_traits>

namespace MML {

	///////////////////////////////////////////////////////////////////////////////////////////
	// CurveFittingConfig - Configuration for curve fitting detailed APIs
	///////////////////////////////////////////////////////////////////////////////////////////
	struct CurveFittingConfig : public EvaluationConfigBase {
		// Inherits: estimate_error, check_finite, exception_policy
		// estimate_error: when true, goodness-of-fit statistics (R², MSE) are populated
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// CurveFittingResult - Result type for curve fitting detailed APIs
	///////////////////////////////////////////////////////////////////////////////////////////
	template<typename Scalar>
	struct CurveFittingResult : public EvaluationResultBase {
		/// Fitted coefficients (for linear: [a, b]; for general: [c0, c1, ...])
		Vector<Scalar> coefficients{};

		/// Euclidean residual norm; GLS uses the covariance-weighted Mahalanobis norm.
		Scalar residual_norm = 0.0;

		/// Coefficient of determination (may be negative); NaN when not requested.
		Scalar r_squared = std::numeric_limits<Scalar>::quiet_NaN();

		/// Mean squared error (average squared residual)
		Scalar mean_squared_error = std::numeric_limits<Scalar>::quiet_NaN();

		/// Adjusted R² (adjusted for number of parameters); NaN if not applicable
		Scalar adjusted_r_squared = std::numeric_limits<Scalar>::quiet_NaN();

		/// Effective rank of the design matrix (from SVD); -1 if not applicable
		int effective_rank = -1;

		/// Condition number of the design matrix; NaN if not applicable
		Scalar condition_number = std::numeric_limits<Scalar>::quiet_NaN();
		/// Condition number restricted to the retained singular subspace.
		Scalar effective_condition_number = std::numeric_limits<Scalar>::quiet_NaN();

		/// Sum of observation weights; zero for unweighted fits
		Scalar weight_sum = 0.0;

		Scalar regularization_parameter = 0.0;
		Scalar coefficient_norm = 0.0;
	};

	///////////////////////////////////////////////////////////////////////////////////////////
	// CurveFittingDetail - Internal helpers for Detailed API execution
	///////////////////////////////////////////////////////////////////////////////////////////
	namespace CurveFittingDetail
	{
		template<typename Scalar>
		void CheckFinite(Scalar value) {
			if (!std::isfinite(value))
				throw MatrixNumericalError("Curve fitting: nonfinite arithmetic result");
		}

		template<typename Scalar>
		void ValidateData(const Vector<Scalar>& x, const Vector<Scalar>& y) {
			if (x.empty() || x.size() != y.size())
				throw CurveFittingError("Curve fitting: nonempty matching data required");
			for (int i = 0; i < x.size(); ++i)
				if (!std::isfinite(x[i]) || !std::isfinite(y[i]))
					throw CurveFittingError("Curve fitting: observations must be finite");
		}

		template<typename Scalar>
		Scalar Norm(const Vector<Scalar>& values) {
			Scalar norm = 0;
			for (int i = 0; i < values.size(); ++i) norm = std::hypot(norm, values[i]);
			CheckFinite(norm);
			return norm;
		}

		template<typename Scalar>
		Scalar DataScale(const Vector<Scalar>& values) {
			Scalar scale = 0;
			for (int i = 0; i < values.size(); ++i) scale = std::max(scale, std::abs(values[i]));
			return scale == 0 ? Scalar(1) : scale;
		}

		template<typename Scalar>
		Vector<Scalar> CenterData(const Vector<Scalar>& values, Scalar scale, Scalar& mean,
			const Vector<Scalar>* weights = nullptr) {
			Vector<Scalar> centered(values.size());
			Scalar offsetMean = 0, total = 0, compensation = 0;
			for (int i = 0; i < values.size(); ++i) {
				// Subtract the anchor before normalizing nearby values; opposite signs need scaling first.
				centered[i] = std::signbit(values[i]) == std::signbit(values[0])
					? (values[i] - values[0]) / scale : values[i] / scale - values[0] / scale;
				const Scalar weight = weights ? (*weights)[i] : Scalar(1);
				total += weight;
				CheckFinite(total);
				const Scalar increment = (centered[i] - offsetMean) * (weight / total) - compensation;
				const Scalar nextMean = offsetMean + increment;
				compensation = (nextMean - offsetMean) - increment;
				offsetMean = nextMean;
			}
			mean = std::fma(offsetMean, scale, values[0]);
			CheckFinite(mean);
			for (int i = 0; i < values.size(); ++i) centered[i] -= offsetMean;
			return centered;
		}

		template<typename Scalar>
		struct FitStatistics {
			Scalar residual_norm = 0;
			Scalar mean_squared_error = std::numeric_limits<Scalar>::quiet_NaN();
			Scalar r_squared = std::numeric_limits<Scalar>::quiet_NaN();
			Scalar adjusted_r_squared = std::numeric_limits<Scalar>::quiet_NaN();
		};

		template<typename Scalar>
		FitStatistics<Scalar> Statistics(const Vector<Scalar>& y, const Vector<Scalar>& residual,
			int rank, bool statistics, const Vector<Scalar>* weights = nullptr) {
			FitStatistics<Scalar> result;
			const Scalar rootWeightScale = weights ? std::sqrt(DataScale(*weights)) : Scalar(1);
			Scalar totalWeight = 0;
			Vector<Scalar> weightedResidual(y.size());
			for (int i = 0; i < y.size(); ++i) {
				const Scalar rootWeight = weights ? std::sqrt((*weights)[i]) / rootWeightScale : Scalar(1);
				totalWeight += rootWeight * rootWeight;
				CheckFinite(totalWeight);
				weightedResidual[i] = residual[i] * rootWeight;
			}
			const Scalar normalizedResidualNorm = Norm(weightedResidual);
			result.residual_norm = normalizedResidualNorm * rootWeightScale;
			CheckFinite(result.residual_norm);
			if (!statistics) return result;
			const Scalar rms = normalizedResidualNorm / std::sqrt(totalWeight);
			result.mean_squared_error = rms * rms;
			CheckFinite(result.mean_squared_error);
			const Scalar scale = DataScale(y);
			if (std::all_of(y.begin(), y.end(), [&](Scalar value) { return value == y[0]; })) {
				// Preserve the legacy perfect/failed constant-response convention, in response units.
				result.r_squared = rms / scale <= PrecisionValues<Scalar>::RoundoffZeroThreshold ? Scalar(1) : Scalar(0);
			} else {
				Scalar mean = 0;
				auto centered = CenterData(y, scale, mean, weights);
				for (int i = 0; i < y.size(); ++i)
					centered[i] *= weights ? std::sqrt((*weights)[i]) / rootWeightScale : Scalar(1);
				const Scalar totalNorm = Norm(centered);
				if (totalNorm == 0)
					throw MatrixNumericalError("Curve fitting: response variation is unrepresentable");
				const Scalar ratio = (rms / scale) / (totalNorm / std::sqrt(totalWeight));
				result.r_squared = 1 - ratio * ratio;
				CheckFinite(result.r_squared);
			}
			if (y.size() > rank) {
				result.adjusted_r_squared = 1 - (1 - result.r_squared) * Scalar(y.size() - 1) / Scalar(y.size() - rank);
				CheckFinite(result.adjusted_r_squared);
			}
			return result;
		}

		template<typename Scalar>
		Vector<std::function<Scalar(Scalar)>> CallableBasis(const Vector<const IRealFunction*>& basis) {
			Vector<std::function<Scalar(Scalar)>> functions(basis.size());
			for (int j = 0; j < basis.size(); ++j) {
				const auto* function = basis[j];
				if (!function)
					throw CurveFittingError("Curve fitting: basis pointers must be non-null");
				functions[j] = [function](Scalar x) { return static_cast<Scalar>((*function)(x)); };
			}
			return functions;
		}

		/// Execute a curve fitting Detailed operation with timing and exception handling.
		template<typename ResultType, typename ComputeFn>
		ResultType ExecuteCurveFittingDetailed(const char* algorithm_name,
		                                      const CurveFittingConfig& config,
		                                      ComputeFn&& compute)
		{
			auto execute = [&]() {
				AlgorithmTimer timer;

				ResultType result = MakeEvaluationSuccessResult<ResultType>(algorithm_name);

				compute(result);

				result.elapsed_time_ms = timer.elapsed_ms();
				return result;
			};

			if (config.exception_policy == EvaluationExceptionPolicy::Propagate)
				return execute();

			try {
				return execute();
			}
			catch (const SingularMatrixError& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::SingularMatrix, ex.what(), algorithm_name);
			}
			catch (const MatrixNumericalError& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::NumericalInstability, ex.what(), algorithm_name);
			}
			catch (const CurveFittingError& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::InvalidInput, ex.what(), algorithm_name);
			}
			catch (const std::exception& ex) {
				return MakeEvaluationFailureResult<ResultType>(
					AlgorithmStatus::AlgorithmSpecificFailure, ex.what(), algorithm_name);
			}
		}
	} // namespace CurveFittingDetail

	template<typename Scalar>
	struct LinearFitResult {
		Scalar a = 0;
		Scalar b = 0;
		Scalar residual_norm = 0;
		Scalar r_squared = std::numeric_limits<Scalar>::quiet_NaN();
		Scalar mean_squared_error = std::numeric_limits<Scalar>::quiet_NaN();
	};

	namespace CurveFittingDetail {
		template<typename Scalar>
		LinearFitResult<Scalar> FitLine(const Vector<Scalar>& x, const Vector<Scalar>& y, bool statistics) {
			ValidateData(x, y);
			const int n = x.size();
			const Scalar xScale = DataScale(x), yScale = DataScale(y);
			Scalar meanX = 0, meanY = 0;
			const auto centeredX = CenterData(x, xScale, meanX), centeredY = CenterData(y, yScale, meanY);
			LinearFitResult<Scalar> result;
			const bool horizontal = n == 1 || (n == 2 &&
				std::abs(x[1] / xScale - x[0] / xScale) < PrecisionValues<Scalar>::RoundoffZeroThreshold);
			if (horizontal) {
				result.b = meanY;
			} else {
				if (std::all_of(x.begin(), x.end(), [&](Scalar value) { return value == x[0]; }))
					throw SingularMatrixError("LinearLeastSquares: all x values are identical", 0);
				Scalar xx = 0, xy = 0;
				for (int i = 0; i < n; ++i) {
					const Scalar dx = centeredX[i], dy = centeredY[i];
					xx = std::fma(dx, dx, xx);
					xy = std::fma(dx, dy, xy);
				}
				if (xx == 0)
					throw SingularMatrixError("LinearLeastSquares: zero centered design scale", 0);
				const Scalar slope = xy / xx;
				int xExponent = 0, yExponent = 0;
				const Scalar xFraction = std::frexp(xScale, &xExponent), yFraction = std::frexp(yScale, &yExponent);
				result.a = std::scalbn(slope * yFraction / xFraction, yExponent - xExponent);
				result.b = std::fma(-result.a, meanX, meanY);
			}
			CheckFinite(result.a);
			CheckFinite(result.b);
			Vector<Scalar> residual(n);
			for (int i = 0; i < n; ++i) {
				const Scalar predicted = std::fma(result.a, x[i], result.b);
				CheckFinite(predicted);
				residual[i] = y[i] - predicted;
			}
			const auto stats = Statistics(y, residual, horizontal ? 1 : 2, statistics);
			result.residual_norm = stats.residual_norm;
			result.r_squared = stats.r_squared;
			result.mean_squared_error = stats.mean_squared_error;
			return result;
		}
	}
	///////////////////////////           LINEAR LEAST SQUARES           ///////////////////////////

	// Fits a linear function y = a*x + b to a set of data points using least squares method
	// Complexity: O(n), using centered, unit-scaled moments.
	//
	// Parameters:
	//   x_data - Vector of x-coordinates (independent variable)
	//   y_data - Vector of y-coordinates (dependent variable)
	//   a      - Output: slope coefficient
	//   b      - Output: intercept coefficient
	//
	// Returns:
	//   residual_norm - The norm of the residual vector (measure of fit quality)
	//
	// Mathematical formulation:
	//   Minimize: sum_i (y_i - (a*x_i + b))^2
	//
	//   a = sum((x_i-mean_x)*(y_i-mean_y)) / sum((x_i-mean_x)^2)
	//   b = mean_y - a*mean_x
	//
	// Throws:
	//   std::invalid_argument if x_data and y_data have different sizes or are empty
	//   std::runtime_error if the normal equations system is singular (e.g., all x values identical)
	template<typename Scalar>
	Scalar LinearLeastSquares(const Vector<Scalar>& x_data, const Vector<Scalar>& y_data, Scalar& a, Scalar& b) {
		const auto fit = CurveFittingDetail::FitLine(x_data, y_data, false);
		a = fit.a;
		b = fit.b;
		return fit.residual_norm;
	}


	// Fits a linear function y = a*x + b and returns detailed statistics
	template<typename Scalar>
	LinearFitResult<Scalar> LinearLeastSquaresDetailed(const Vector<Scalar>& x_data, const Vector<Scalar>& y_data) {
		return CurveFittingDetail::FitLine(x_data, y_data, true);
	}


	///////////////////////////    GENERAL LINEAR LEAST SQUARES    ///////////////////////////

	// Result structure for general linear least squares fitting
	// Contains fitted coefficients and comprehensive statistics
	template<typename Scalar>
	struct GeneralLinearFitResult {
		Vector<Scalar> coefficients; // Fitted coefficients for each basis function
		Scalar residual_norm;		   // Euclidean norm, or covariance-weighted norm for GLS
		Scalar r_squared;			   // Coefficient of determination; may be negative
		Scalar mean_squared_error;   // Average squared residual
		Scalar adjusted_r_squared;   // R² adjusted for number of parameters
		int num_data_points;	   // Number of data points (m)
		int num_basis_functions;   // Number of basis functions (n)
		int effective_rank;		   // Effective rank of the design matrix (from SVD)
		Scalar condition_number;	   // Condition number of the design matrix (ratio of largest/smallest singular values)
		Scalar effective_condition_number; // Ratio on the retained singular subspace
		Scalar weight_sum;          // Sum of observation weights; zero for unweighted fits
		Scalar regularization_parameter;
		Scalar coefficient_norm;

		GeneralLinearFitResult()
			: residual_norm(0)
			, r_squared(std::numeric_limits<Scalar>::quiet_NaN())
			, mean_squared_error(std::numeric_limits<Scalar>::quiet_NaN())
			, adjusted_r_squared(std::numeric_limits<Scalar>::quiet_NaN())
			, num_data_points(0)
			, num_basis_functions(0)
			, effective_rank(0)
			, condition_number(std::numeric_limits<Scalar>::quiet_NaN())
			, effective_condition_number(std::numeric_limits<Scalar>::quiet_NaN())
			, weight_sum(0)
			, regularization_parameter(0)
			, coefficient_norm(0) {}

		// Evaluate the fitted function at a point x using the basis functions
		Scalar evaluate(Scalar x, const Vector<const IRealFunction*>& basis_functions) const {
			return evaluate(x, CurveFittingDetail::CallableBasis<Scalar>(basis_functions));
		}

		// Evaluate using a function wrapper (std::function version)
		Scalar evaluate(Scalar x, const Vector<std::function<Scalar(Scalar)>>& basis_functions) const {
			if (coefficients.size() != basis_functions.size())
				throw CurveFittingError("GeneralLinearFitResult::evaluate: mismatched basis functions count");
			if (!std::isfinite(x))
				throw CurveFittingError("GeneralLinearFitResult::evaluate: coordinate must be finite");
			for (int i = 0; i < basis_functions.size(); ++i)
				if (!basis_functions[i])
					throw CurveFittingError("GeneralLinearFitResult::evaluate: basis functions must be callable");
			Scalar result = 0.0;
			for (int i = 0; i < coefficients.size(); i++) {
				const Scalar value = basis_functions[i](x);
				if (!std::isfinite(value))
					throw CurveFittingError("GeneralLinearFitResult::evaluate: basis values must be finite");
				CurveFittingDetail::CheckFinite(coefficients[i]);
				result = std::fma(coefficients[i], value, result);
				CurveFittingDetail::CheckFinite(result);
			}
			return result;
		}
	};

	namespace CurveFittingDetail {
		template<typename Scalar>
		void SetConditionNumbers(SVDecompositionSolver<Scalar>& svd, GeneralLinearFitResult<Scalar>& result) {
			result.effective_rank = svd.Rank();
			const auto singular = svd.getW();
			const int rank = result.effective_rank;
			result.effective_condition_number = rank > 0 ? singular[0] / singular[rank - 1]
				: std::numeric_limits<Scalar>::infinity();
			result.condition_number = rank == singular.size() ? result.effective_condition_number
				: std::numeric_limits<Scalar>::infinity();
		}

		template<typename Scalar>
		GeneralLinearFitResult<Scalar> FitLinearModel(const Vector<Scalar>& x, const Vector<Scalar>& y,
			const Vector<std::function<Scalar(Scalar)>>& basis, bool statistics,
			const Vector<Scalar>* weights = nullptr, const Vector<Scalar>* regularization = nullptr) {
			ValidateData(x, y);
			const int m = x.size(), n = basis.size();
			if (n == 0 || m < n)
				throw CurveFittingError("Curve fitting: 1..m basis functions required");
			for (int j = 0; j < n; ++j)
				if (!basis[j])
					throw CurveFittingError("Curve fitting: basis functions must be callable");
			GeneralLinearFitResult<Scalar> result;
			if (weights) {
				if (weights->size() != m)
					throw CurveFittingError("Curve fitting: weight count must match observations");
				for (int i = 0; i < m; ++i) {
					if (!std::isfinite((*weights)[i]) || (*weights)[i] <= 0)
						throw CurveFittingError("Curve fitting: weights must be finite and positive");
					result.weight_sum += (*weights)[i];
					CheckFinite(result.weight_sum);
				}
			}
			if (regularization) {
				if (regularization->size() != n || n > std::numeric_limits<int>::max() - m)
					throw CurveFittingError("Curve fitting: invalid regularization dimensions");
				for (int j = 0; j < n; ++j)
					if (!std::isfinite((*regularization)[j]) || (*regularization)[j] < 0)
						throw CurveFittingError("Curve fitting: regularization must be finite and nonnegative");
			}
			Matrix<Scalar> design(m, n), system(m + (regularization ? n : 0), n);
			Vector<Scalar> rhs(system.rows());
			const Scalar rootWeightScale = weights ? std::sqrt(DataScale(*weights)) : Scalar(1);
			for (int i = 0; i < m; ++i) {
				const Scalar scale = weights ? std::sqrt((*weights)[i]) / rootWeightScale : Scalar(1);
				rhs[i] = scale * y[i];
				CheckFinite(rhs[i]);
				for (int j = 0; j < n; ++j) {
					design(i, j) = basis[j](x[i]);
					if (!std::isfinite(design(i, j)))
						throw CurveFittingError("Curve fitting: basis values must be finite");
					system(i, j) = scale * design(i, j);
					CheckFinite(system(i, j));
				}
			}
			if (regularization)
				for (int j = 0; j < n; ++j) system(m + j, j) = (*regularization)[j];
			SVDecompositionSolver<Scalar> svd(system);
			result.coefficients = svd.Solve(rhs);
			for (int j = 0; j < n; ++j) CheckFinite(result.coefficients[j]);
			result.num_data_points = m;
			result.num_basis_functions = n;
			SetConditionNumbers(svd, result);
			Vector<Scalar> residual(m);
			for (int i = 0; i < m; ++i) {
				Scalar prediction = 0;
				for (int j = 0; j < n; ++j) {
					prediction = std::fma(result.coefficients[j], design(i, j), prediction);
					CheckFinite(prediction);
				}
				residual[i] = y[i] - prediction;
			}
			const auto stats = Statistics(y, residual, result.effective_rank, statistics, weights);
			result.residual_norm = stats.residual_norm;
			result.mean_squared_error = stats.mean_squared_error;
			result.r_squared = stats.r_squared;
			result.adjusted_r_squared = stats.adjusted_r_squared;
			if (regularization) {
				result.coefficient_norm = Norm(result.coefficients);
				Scalar maximum = 0;
				for (int j = 0; j < n; ++j) maximum = std::max(maximum, (*regularization)[j]);
				result.regularization_parameter = maximum * maximum;
				CheckFinite(result.regularization_parameter);
			}
			return result;
		}

		template<typename Scalar>
		Vector<std::function<Scalar(Scalar)>> PolynomialBasis(int degree, int observations) {
			if (degree < 0 || degree >= observations)
				throw CurveFittingError("PolynomialFit: degree must be nonnegative and smaller than the observation count");
			Vector<std::function<Scalar(Scalar)>> basis(degree + 1);
			for (int power = 0; power <= degree; ++power)
				basis[power] = [power](Scalar x) {
					return power == 0 ? Scalar(1) : static_cast<Scalar>(std::pow(x, power));
				};
			return basis;
		}
	}

	// Fits a linear combination of arbitrary basis functions to data using least squares
	//
	// Mathematical formulation:
	//   Find coefficients c_0, c_1, ..., c_{n-1} that minimize:
	//     sum_i (y_i - sum_j c_j * f_j(x_i))^2
	//
	//   This is equivalent to solving the overdetermined system:
	//     A * c = y
	//   where A[i][j] = f_j(x_i) is the design matrix
	//
	//   The solution uses SVD: A = U * W * V^T
	//   Least squares solution: c = V * W^{-1} * U^T * y (pseudoinverse)
	//
	// Parameters:
	//   x_data          - Vector of x-coordinates (independent variable)
	//   y_data          - Vector of y-coordinates (dependent variable)
	//   basis_functions - Vector of pointers to basis functions f_j(x)
	//
	// Returns:
	//   GeneralLinearFitResult containing coefficients and statistics
	//
	// Example basis functions:
	//   - Polynomial: {1, x, x^2, x^3, ...}
	//   - Fourier:    {1, sin(x), cos(x), sin(2x), cos(2x), ...}
	//   - Custom:     {exp(-x^2), log(x+1), sqrt(x), ...}
	//
	// Throws:
	//   std::invalid_argument if vectors have mismatched sizes, are empty,
	//                         or basis_functions is empty
	template<typename Scalar>
	GeneralLinearFitResult<Scalar> GeneralLinearLeastSquares(const Vector<Scalar>& x_data, const Vector<Scalar>& y_data,
														   const Vector<const IRealFunction*>& basis_functions) {
		return CurveFittingDetail::FitLinearModel(x_data, y_data,
			CurveFittingDetail::CallableBasis<Scalar>(basis_functions), true);
	}

	// Weighted least squares minimizes sum_i weights[i] * residual[i]^2.
	// The SVD solve uses the equivalent row-scaled system sqrt(W) A c = sqrt(W) y.
	template<typename Scalar>
	GeneralLinearFitResult<Scalar> WeightedGeneralLinearLeastSquares(
		const Vector<Scalar>& x_data, const Vector<Scalar>& y_data,
		const Vector<Scalar>& weights,
		const Vector<std::function<Scalar(Scalar)>>& basis_functions) {
		return CurveFittingDetail::FitLinearModel(x_data, y_data, basis_functions, true, &weights);
	}

	namespace CurveFittingDetail {
		template<typename Scalar>
		GeneralLinearFitResult<Scalar> FitGeneralizedLinear(
			const Vector<Scalar>& x, const Vector<Scalar>& y, const Matrix<Scalar>& covariance,
			const Vector<std::function<Scalar(Scalar)>>& basis, bool statistics)
		{
			static_assert(std::is_floating_point_v<Scalar>, "GLS requires a real floating-point scalar");
			const int m = x.size(), n = basis.size();
			if (m == 0 || n == 0 || y.size() != m || m < n)
				throw CurveFittingError("GeneralizedLinearLeastSquares: nonempty matching data and 1..m basis functions required");
			if (covariance.rows() != m || covariance.cols() != m)
				throw CurveFittingError("GeneralizedLinearLeastSquares: covariance must be m by m");
			for (int j = 0; j < n; ++j)
				if (!basis[j])
					throw CurveFittingError("GeneralizedLinearLeastSquares: basis functions must be callable");
			Scalar scale = 0;
			for (int i = 0; i < m; ++i) {
				if (!std::isfinite(x[i]) || !std::isfinite(y[i]))
					throw CurveFittingError("GeneralizedLinearLeastSquares: observations must be finite");
				for (int j = 0; j < m; ++j) {
					if (!std::isfinite(covariance(i, j)))
						throw CurveFittingError("GeneralizedLinearLeastSquares: covariance must be finite");
					scale = std::max(scale, std::abs(covariance(i, j)));
				}
			}
			if (scale == 0)
				throw CurveFittingError("GeneralizedLinearLeastSquares: covariance must be positive definite");
			// Normalize units before Cholesky's numerical pivot checks; do not invert covariance.
			Matrix<Scalar> normalized(m, m);
			for (int i = 0; i < m; ++i)
				for (int j = 0; j <= i; ++j) {
					const Scalar a = covariance(i, j) / scale, b = covariance(j, i) / scale;
					if (std::abs(a - b) > PrecisionValues<Scalar>::IsMatrixSymmetricTolerance)
						throw CurveFittingError("GeneralizedLinearLeastSquares: covariance must be symmetric");
					normalized(i, j) = normalized(j, i) = a / 2 + b / 2;
				}
			const auto factor = [&]() {
				try {
					return CholeskySolver<Scalar>(normalized);
				}
				catch (const SingularMatrixError& error) {
					throw CurveFittingError(std::string("GeneralizedLinearLeastSquares: covariance must be numerically positive definite: ") + error.what());
				}
			}();
			const auto& lower = factor.L();
			auto whiten = [&](const Vector<Scalar>& values) {
				Vector<Scalar> transformed(m);
				for (int i = 0; i < m; ++i) {
					Scalar value = values[i];
					for (int k = 0; k < i; ++k) value -= lower(i, k) * transformed[k];
					transformed[i] = value / lower(i, i);
					if (!std::isfinite(transformed[i]))
						throw MatrixNumericalError("GeneralizedLinearLeastSquares: non-finite whitened value");
				}
				return transformed;
			};
			Matrix<Scalar> design(m, n), whitenedDesign(m, n);
			for (int j = 0; j < n; ++j) {
				Vector<Scalar> column(m);
				for (int i = 0; i < m; ++i) {
					column[i] = basis[j](x[i]);
					if (!std::isfinite(column[i]))
						throw CurveFittingError("GeneralizedLinearLeastSquares: basis values must be finite");
					design(i, j) = column[i];
				}
				const auto transformed = whiten(column);
				for (int i = 0; i < m; ++i) whitenedDesign(i, j) = transformed[i];
			}
			const auto whitenedY = whiten(y);
			SVDecompositionSolver<Scalar> svd(whitenedDesign);
			GeneralLinearFitResult<Scalar> result;
			result.coefficients = svd.Solve(whitenedY);
			result.num_data_points = m;
			result.num_basis_functions = n;
			SetConditionNumbers(svd, result);
			Vector<Scalar> residual(m);
			for (int i = 0; i < m; ++i) {
				Scalar predicted = 0;
				for (int j = 0; j < n; ++j) {
					if (!std::isfinite(result.coefficients[j]))
						throw MatrixNumericalError("GeneralizedLinearLeastSquares: non-finite fitted coefficient");
					predicted += design(i, j) * result.coefficients[j];
				}
				residual[i] = y[i] - predicted;
			}
			const Scalar normalizedResidualNorm = Norm(whiten(residual));
			result.residual_norm = normalizedResidualNorm / std::sqrt(scale);
			if (!std::isfinite(result.residual_norm))
				throw MatrixNumericalError("GeneralizedLinearLeastSquares: non-finite covariance-weighted residual norm");
			result.r_squared = result.mean_squared_error = result.adjusted_r_squared =
				std::numeric_limits<Scalar>::quiet_NaN();
			if (statistics) {
				result.mean_squared_error = (result.residual_norm / std::sqrt(Scalar(m))) *
					(result.residual_norm / std::sqrt(Scalar(m)));
				if (!std::isfinite(result.mean_squared_error))
					throw MatrixNumericalError("GeneralizedLinearLeastSquares: non-finite covariance-weighted MSE");
				if (std::all_of(y.begin(), y.end(), [&](Scalar value) { return value == y[0]; }))
					return result;
				const auto whitenedOnes = whiten(Vector<Scalar>(m, Scalar(1)));
				Scalar numerator = 0, denominator = 0;
				for (int i = 0; i < m; ++i) {
					numerator += whitenedOnes[i] * whitenedY[i];
					denominator += whitenedOnes[i] * whitenedOnes[i];
				}
				if (!std::isfinite(numerator) || !std::isfinite(denominator) || denominator <= 0)
					throw MatrixNumericalError("GeneralizedLinearLeastSquares: non-finite generalized-mean weights");
				const Scalar mean = numerator / denominator;
				if (!std::isfinite(mean))
					throw MatrixNumericalError("GeneralizedLinearLeastSquares: non-finite generalized mean");
				Vector<Scalar> centered(m);
				for (int i = 0; i < m; ++i) centered[i] = y[i] - mean;
				const Scalar totalNorm = Norm(whiten(centered));
				if (!std::isfinite(totalNorm))
					throw MatrixNumericalError("GeneralizedLinearLeastSquares: non-finite centered norm");
				if (totalNorm > 0) {
					const Scalar ratio = normalizedResidualNorm / totalNorm;
					result.r_squared = 1 - ratio * ratio;
					if (!std::isfinite(result.r_squared))
						throw MatrixNumericalError("GeneralizedLinearLeastSquares: non-finite R-squared");
					const int rank = result.effective_rank;
					if (m > rank) {
						result.adjusted_r_squared = 1 - (1 - result.r_squared) * Scalar(m - 1) / Scalar(m - rank);
						if (!std::isfinite(result.adjusted_r_squared))
							throw MatrixNumericalError("GeneralizedLinearLeastSquares: non-finite adjusted R-squared");
					}
				}
				// Constant response has no centered variation: R-squared is undefined.
			}
			return result;
		}

	}

	/// Minimize (y-Ac)^T covariance^{-1} (y-Ac) for a finite symmetric positive-definite
	/// observation covariance. Cholesky whitening + SVD also handles rank-deficient designs.
	/// residual_norm is the Mahalanobis norm; MSE is its square divided by observations.
	template<typename Scalar>
	GeneralLinearFitResult<Scalar> GeneralizedLinearLeastSquares(
		const Vector<Scalar>& x_data, const Vector<Scalar>& y_data, const Matrix<Scalar>& covariance,
		const Vector<std::function<Scalar(Scalar)>>& basis_functions)
	{
		return CurveFittingDetail::FitGeneralizedLinear(x_data, y_data, covariance, basis_functions, true);
	}

	template<typename Scalar>
	GeneralLinearFitResult<Scalar> GeneralizedLinearLeastSquares(
		const Vector<Scalar>& x_data, const Vector<Scalar>& y_data, const Matrix<Scalar>& covariance,
		const Vector<const IRealFunction*>& basis_functions)
	{
		return GeneralizedLinearLeastSquares(x_data, y_data, covariance,
			CurveFittingDetail::CallableBasis<Scalar>(basis_functions));
	}

	// Overload using std::function for more flexible basis function specification
	// This allows lambda functions, function pointers, and functors
	template<typename Scalar>
	GeneralLinearFitResult<Scalar> GeneralLinearLeastSquares(const Vector<Scalar>& x_data, const Vector<Scalar>& y_data,
														   const Vector<std::function<Scalar(Scalar)>>& basis_functions) {
		return CurveFittingDetail::FitLinearModel(x_data, y_data, basis_functions, true);
	}

	template<typename Scalar>
	GeneralLinearFitResult<Scalar> TikhonovGeneralLinearLeastSquares(
		const Vector<Scalar>& x_data, const Vector<Scalar>& y_data,
		const Vector<std::function<Scalar(Scalar)>>& basis_functions,
		const Vector<Scalar>& regularization_diagonal)
	{
		return CurveFittingDetail::FitLinearModel<Scalar>(x_data, y_data, basis_functions, true,
			nullptr, &regularization_diagonal);
	}

	template<typename Scalar>
	GeneralLinearFitResult<Scalar> RidgeGeneralLinearLeastSquares(
		const Vector<Scalar>& x_data, const Vector<Scalar>& y_data,
		const Vector<std::function<Scalar(Scalar)>>& basis_functions,
		Scalar lambda, bool regularize_constant = false)
	{
		if (!std::isfinite(lambda) || lambda < 0)
			throw CurveFittingError("RidgeGeneralLinearLeastSquares: lambda must be finite and non-negative");
		Vector<Scalar> diagonal(static_cast<int>(basis_functions.size()), std::sqrt(lambda));
		if (!regularize_constant && diagonal.size() > 0) diagonal[0] = 0.0;
		return TikhonovGeneralLinearLeastSquares(x_data, y_data, basis_functions, diagonal);
	}

	template<typename Basis, typename Scalar = Real>
	Vector<std::function<Scalar(Scalar)>> MakeOrthogonalPolynomialBasis(int max_degree)
	{
		if (max_degree < 0)
			throw CurveFittingError("MakeOrthogonalPolynomialBasis: max_degree must be non-negative");
		Vector<std::function<Scalar(Scalar)>> functions(max_degree + 1);
		for (int degree = 0; degree <= max_degree; ++degree)
			functions[degree] = [degree](Scalar x) { return static_cast<Scalar>(Basis().Evaluate(degree, x)); };
		return functions;
	}

	template<typename Scalar = Real>
	auto MakeChebyshevFitBasis(int max_degree) { return MakeOrthogonalPolynomialBasis<ChebyshevBasis, Scalar>(max_degree); }
	template<typename Scalar = Real>
	auto MakeLegendreFitBasis(int max_degree) { return MakeOrthogonalPolynomialBasis<LegendreBasis, Scalar>(max_degree); }
	template<typename Scalar = Real>
	auto MakeHermiteFitBasis(int max_degree) { return MakeOrthogonalPolynomialBasis<HermiteBasis, Scalar>(max_degree); }
	template<typename Scalar = Real>
	auto MakeLaguerreFitBasis(int max_degree) { return MakeOrthogonalPolynomialBasis<LaguerreBasis, Scalar>(max_degree); }


	///////////////////////////    POLYNOMIAL FITTING    ///////////////////////////

	// Fits a polynomial of specified degree to data: y = c_0 + c_1*x + c_2*x^2 + ... + c_n*x^n
	// Complexity: O(n*m² + m³) where n=data size, m=degree+1 (dominated by solving normal equations)
	//
	// This is a convenience wrapper around GeneralLinearLeastSquares with polynomial basis functions
	//
	// Parameters:
	//   x_data - Vector of x-coordinates
	//   y_data - Vector of y-coordinates
	//   degree - Degree of the polynomial (0 = constant, 1 = linear, 2 = quadratic, etc.)
	//
	// Returns:
	//   GeneralLinearFitResult where coefficients[i] is the coefficient of x^i
	//
	// Example: degree=2 fits y = c[0] + c[1]*x + c[2]*x^2
	template<typename Scalar>
	GeneralLinearFitResult<Scalar> PolynomialFit(const Vector<Scalar>& x_data, const Vector<Scalar>& y_data, int degree) {
		return GeneralLinearLeastSquares(x_data, y_data,
			CurveFittingDetail::PolynomialBasis<Scalar>(degree, x_data.size()));
	}

	template<typename Scalar>
	GeneralLinearFitResult<Scalar> WeightedPolynomialFit(
		const Vector<Scalar>& x_data, const Vector<Scalar>& y_data,
		const Vector<Scalar>& weights, int degree) {
		return WeightedGeneralLinearLeastSquares(x_data, y_data, weights,
			CurveFittingDetail::PolynomialBasis<Scalar>(degree, x_data.size()));
	}

	// Evaluate a fitted polynomial at a point
	template<typename Scalar>
	Scalar EvaluatePolynomial(Scalar x, const Vector<Scalar>& coefficients) {
		// Use Horner's method for numerical stability: c_0 + x*(c_1 + x*(c_2 + ...))
		int n = coefficients.size();
		if (n == 0)
			return 0.0;

		Scalar result = coefficients[n - 1];
		for (int i = n - 2; i >= 0; i--) {
			result = coefficients[i] + x * result;
		}
		return result;
	}


	///////////////////////////    DETAILED API    ///////////////////////////

	/// Linear least squares fit with full instrumentation.
	/// Returns CurveFittingResult with coefficients=[a, b], plus R², MSE, timing, and AlgorithmStatus.
	template<typename Scalar>
	CurveFittingResult<Scalar> LinearFitDetailed(
		const Vector<Scalar>& x_data,
		const Vector<Scalar>& y_data,
		const CurveFittingConfig& config = {})
	{
		return CurveFittingDetail::ExecuteCurveFittingDetailed<CurveFittingResult<Scalar>>(
			"LinearLeastSquares", config,
			[&](CurveFittingResult<Scalar>& result) {
				auto fit = CurveFittingDetail::FitLine(x_data, y_data, config.estimate_error);

				result.coefficients = Vector<Scalar>(2);
				result.coefficients[0] = fit.a;
				result.coefficients[1] = fit.b;
				result.residual_norm   = fit.residual_norm;
				result.r_squared       = fit.r_squared;
				result.mean_squared_error = fit.mean_squared_error;
			});
	}

	/// General linear least squares fit with full instrumentation.
	/// Accepts std::function basis functions.
	template<typename Scalar>
	CurveFittingResult<Scalar> GeneralLinearFitDetailed(
		const Vector<Scalar>& x_data,
		const Vector<Scalar>& y_data,
		const Vector<std::function<Scalar(Scalar)>>& basis_functions,
		const CurveFittingConfig& config = {})
	{
		return CurveFittingDetail::ExecuteCurveFittingDetailed<CurveFittingResult<Scalar>>(
			"GeneralLinearLeastSquares", config,
			[&](CurveFittingResult<Scalar>& result) {
				auto fit = CurveFittingDetail::FitLinearModel(x_data, y_data, basis_functions, config.estimate_error);

				result.coefficients       = fit.coefficients;
				result.residual_norm      = fit.residual_norm;
				result.r_squared          = fit.r_squared;
				result.mean_squared_error = fit.mean_squared_error;
				result.adjusted_r_squared = fit.adjusted_r_squared;
				result.effective_rank     = fit.effective_rank;
				result.condition_number   = fit.condition_number;
				result.effective_condition_number = fit.effective_condition_number;
			});
	}

	/// Polynomial fit with full instrumentation.
	template<typename Scalar>
	CurveFittingResult<Scalar> PolynomialFitDetailed(
		const Vector<Scalar>& x_data,
		const Vector<Scalar>& y_data,
		int degree,
		const CurveFittingConfig& config = {})
	{
		return CurveFittingDetail::ExecuteCurveFittingDetailed<CurveFittingResult<Scalar>>(
			"PolynomialFit", config,
			[&](CurveFittingResult<Scalar>& result) {
				auto fit = CurveFittingDetail::FitLinearModel(x_data, y_data,
					CurveFittingDetail::PolynomialBasis<Scalar>(degree, x_data.size()), config.estimate_error);

				result.coefficients       = fit.coefficients;
				result.residual_norm      = fit.residual_norm;
				result.r_squared          = fit.r_squared;
				result.mean_squared_error = fit.mean_squared_error;
				result.adjusted_r_squared = fit.adjusted_r_squared;
				result.effective_rank     = fit.effective_rank;
				result.condition_number   = fit.condition_number;
				result.effective_condition_number = fit.effective_condition_number;
			});
	}

	template<typename Scalar>
	CurveFittingResult<Scalar> WeightedGeneralLinearFitDetailed(
		const Vector<Scalar>& x_data, const Vector<Scalar>& y_data,
		const Vector<Scalar>& weights,
		const Vector<std::function<Scalar(Scalar)>>& basis_functions,
		const CurveFittingConfig& config = {}) {
		return CurveFittingDetail::ExecuteCurveFittingDetailed<CurveFittingResult<Scalar>>(
			"WeightedGeneralLinearLeastSquares", config,
			[&](CurveFittingResult<Scalar>& result) {
				auto fit = CurveFittingDetail::FitLinearModel(x_data, y_data, basis_functions, config.estimate_error, &weights);
				result.coefficients = fit.coefficients;
				result.residual_norm = fit.residual_norm;
				result.r_squared = fit.r_squared;
				result.mean_squared_error = fit.mean_squared_error;
				result.adjusted_r_squared = fit.adjusted_r_squared;
				result.effective_rank = fit.effective_rank;
				result.condition_number = fit.condition_number;
				result.effective_condition_number = fit.effective_condition_number;
				result.weight_sum = fit.weight_sum;
				result.function_evaluations = x_data.size() * basis_functions.size();
			});
	}

	/// GLS with structured status/timing. Completed denotes a direct fit, not a checked
	/// accuracy target; estimate_error=false leaves optional goodness-of-fit fields NaN.
	template<typename Scalar>
	CurveFittingResult<Scalar> GeneralizedLinearFitDetailed(
		const Vector<Scalar>& x_data, const Vector<Scalar>& y_data, const Matrix<Scalar>& covariance,
		const Vector<std::function<Scalar(Scalar)>>& basis_functions,
		const CurveFittingConfig& config = {})
	{
		return CurveFittingDetail::ExecuteCurveFittingDetailed<CurveFittingResult<Scalar>>(
			"GeneralizedLinearLeastSquares", config, [&](CurveFittingResult<Scalar>& result) {
				const auto fit = CurveFittingDetail::FitGeneralizedLinear(
					x_data, y_data, covariance, basis_functions, config.estimate_error);
				result.coefficients = fit.coefficients;
				result.residual_norm = fit.residual_norm;
				result.mean_squared_error = fit.mean_squared_error;
				result.r_squared = fit.r_squared;
				result.adjusted_r_squared = fit.adjusted_r_squared;
				result.effective_rank = fit.effective_rank;
				result.condition_number = fit.condition_number;
				result.effective_condition_number = fit.effective_condition_number;
				result.function_evaluations = x_data.size() * basis_functions.size();
				result.status = AlgorithmStatus::Completed;
			});
	}

	template<typename Scalar>
	CurveFittingResult<Scalar> GeneralizedLinearFitDetailed(
		const Vector<Scalar>& x_data, const Vector<Scalar>& y_data, const Matrix<Scalar>& covariance,
		const Vector<const IRealFunction*>& basis_functions,
		const CurveFittingConfig& config = {})
	{
		// Keep pointer validation inside the common exception-policy boundary.
		return CurveFittingDetail::ExecuteCurveFittingDetailed<CurveFittingResult<Scalar>>(
			"GeneralizedLinearLeastSquares", config, [&](CurveFittingResult<Scalar>& result) {
				result = GeneralizedLinearFitDetailed(x_data, y_data, covariance,
					CurveFittingDetail::CallableBasis<Scalar>(basis_functions), config);
			});
	}

} // namespace MML

#endif // MML_CURVE_FITTING_H
