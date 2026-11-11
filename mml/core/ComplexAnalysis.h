///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///  File:        ComplexAnalysis.h                                                   ///
///  Description: Contours, Cauchy formulas, residues and argument principle           ///
///  References: Ahlfors, Complex Analysis (1979); Numerical Recipes, 3rd ed.           ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                         ///
///  License:     MIT License (see LICENSE.md)                                         ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COMPLEX_ANALYSIS_H
#define MML_COMPLEX_ANALYSIS_H

#include <mml/MMLBase.h>
#include <mml/interfaces/IComplexFunction.h>
#include <mml/base/ComplexFunction.h>
#include <mml/core/Derivation/DerivationComplex.h>
#include <mml/core/Integration/GaussKronrod.h>

#include <array>
#include <cmath>
#include <functional>
#include <limits>

namespace MML
{
	namespace ComplexAnalysis
	{
		using ContourIntegrationResult = BasicIntegrationResult<Complex>;

		struct ContourIntegrationConfig : IntegrationConfig {
			Real tolerance = Precision::ComplexAnalysisTolerance;
			int max_depth = Precision::ComplexAnalysisMaxDepth;
			int max_evaluations = Precision::ComplexAnalysisMaxEvaluations;
			int initial_panels = Precision::ComplexAnalysisInitialPanels;
		};

		/// The complex payload is the unrounded estimate. Certification is numerical, not a proof of analyticity.
		struct [[nodiscard]] IntegerContourResult : ContourIntegrationResult {
			int integer_value = 0;
			bool integer_certified = false;
			Real derivative_error_estimate = 0;

			int IntegerOrThrow(const char* context) const {
				ValueOrThrow(context);
				if (!integer_certified)
					throw ConvergenceError(std::string(context) + ": integer estimate is not certified",
						iterations, error_estimate);
				return integer_value;
			}
		};

		namespace Detail
		{
			inline void ValidatePoint(Complex point) {
				if (!std::isfinite(point.real()) || !std::isfinite(point.imag()))
					throw ArgumentError("Complex analysis requires finite points");
			}

			inline void ValidateCircle(Complex center, Real radius) {
				ValidatePoint(center);
				if (!std::isfinite(radius) || radius <= 0)
					throw ArgumentError("Complex contour radius must be finite and positive");
				for (Real coordinate : {center.real(), center.imag()})
					if (!std::isfinite(coordinate + radius) || !std::isfinite(coordinate - radius) ||
						coordinate + radius == coordinate || coordinate - radius == coordinate)
						throw ArgumentError("Complex contour is not representable at its center and radius");
			}

			inline void ValidateConfig(const ContourIntegrationConfig& config) {
				IntegrationDetail::ValidateAdaptive(config.tolerance, 0, config.max_depth, config.max_evaluations);
				if (config.initial_panels <= 0)
					throw ArgumentError("Complex contour integration requires positive initial_panels");
			}

			inline ContourIntegrationConfig WithTolerance(Real tolerance) {
				ContourIntegrationConfig config;
				config.tolerance = tolerance;
				return config;
			}

			class CountedFunction : public IComplexFunction {
				const IComplexFunction& _function;
				ContourIntegrationResult& _result;
				const ContourIntegrationConfig& _config;
			public:
				CountedFunction(const IComplexFunction& function, ContourIntegrationResult& result,
					const ContourIntegrationConfig& config) : _function(function), _result(result), _config(config) {}
				Complex operator()(Complex point) const override {
					ValidatePoint(point);
					if (_result.function_evaluations >= _config.max_evaluations)
						throw ResourceLimitError("Complex contour function-evaluation budget exhausted");
					++_result.function_evaluations;
					const Complex value = _function(point);
					if (_config.check_finite) IntegrationDetail::RequireFinite(value, "Complex integrand");
					return value;
				}
			};
		}

		class CircleContour : public IRealToComplexFunction {
			Complex _center;
			Real _radius;
		public:
			CircleContour(Complex center, Real radius) : _center(center), _radius(radius) {
				Detail::ValidateCircle(center, radius);
			}
			Complex operator()(Real t) const override {
				return _center + _radius * Complex(std::cos(t), std::sin(t));
			}
			Complex derivative(Real t) const { return _radius * Complex(-std::sin(t), std::cos(t)); }
			Complex center() const { return _center; }
			Real radius() const { return _radius; }
			Real t_start() const { return 0; }
			Real t_end() const { return 2 * Constants::PI; }
		};

		class LineSegmentContour : public IRealToComplexFunction {
			Complex _z1, _z2;
		public:
			LineSegmentContour(Complex z1, Complex z2) : _z1(z1), _z2(z2) {
				Detail::ValidatePoint(z1);
				Detail::ValidatePoint(z2);
				Detail::ValidatePoint(z2 - z1);
			}
			Complex operator()(Real t) const override {
				return {std::lerp(_z1.real(), _z2.real(), t), std::lerp(_z1.imag(), _z2.imag(), t)};
			}
			Complex derivative(Real) const { return _z2 - _z1; }
			Real t_start() const { return 0; }
			Real t_end() const { return 1; }
		};

		class ArcContour : public IRealToComplexFunction {
			Complex _center;
			Real _radius, _t1, _t2;
		public:
			ArcContour(Complex center, Real radius, Real t1, Real t2)
				: _center(center), _radius(radius), _t1(t1), _t2(t2) {
				Detail::ValidateCircle(center, radius);
				if (!std::isfinite(t1) || !std::isfinite(t2))
					throw ArgumentError("Complex contour parameter endpoints must be finite");
			}
			Complex operator()(Real t) const override {
				return _center + _radius * Complex(std::cos(t), std::sin(t));
			}
			Complex derivative(Real t) const { return _radius * Complex(-std::sin(t), std::cos(t)); }
			Real t_start() const { return _t1; }
			Real t_end() const { return _t2; }
		};

		namespace Detail
		{
			// G7/K15 uses nonuniform nodes; the five-node dyadic Simpson test aliases harmonic four.
			template<class Function>
			ContourIntegrationResult Kronrod15(Function& function, Real a, Real b) {
				const Real half = IntegrationDetail::HalfWidth(a, b);
				const Real center = std::midpoint(a, b);
				IntegrationDetail::ComplexSum kronrod, gauss;
				IntegrationDetail::Sum floor;
				auto add = [&](Real t, int node) {
					const Complex value = function(t);
					IntegrationDetail::RequireFinite(value, "Complex quadrature sample");
					auto weighted = [&](Real weight) {
						return Complex(IntegrationDetail::Product({half, weight, value.real()}),
							IntegrationDetail::Product({half, weight, value.imag()}));
					};
					const Complex contribution = weighted(Integration::wGK15[node]);
					kronrod.Add(contribution);
					floor.Add(Precision::IntegrationRoundoffFloorScale * std::abs(contribution));
					if (node % 2 == 0) gauss.Add(weighted(Integration::wG7[node / 2]));
				};
				add(center, 0);
				for (int node = 1; node < Integration::GK15_N; ++node) {
					const Real offset = half * Integration::xGK15[node];
					add(center + offset, node);
					add(center - offset, node);
				}
				return {kronrod.Value(), std::max(std::abs(kronrod.Value() - gauss.Value()), floor.Value()), 0, false};
			}

			template<class Function>
			void Integrate(Function& function, Real a, Real b, const ContourIntegrationConfig& config,
				ContourIntegrationResult& result) {
				ValidateConfig(config);
				IntegrationDetail::HalfWidth(a, b);
				if (a == b) {
					result.value = 0;
					result.error_estimate = 0;
					result.error_estimated = result.converged = true;
					return;
				}
				if (config.max_evaluations < 2 || config.initial_panels > (config.max_evaluations - 2) / 15)
					throw ResourceLimitError("Budget cannot cover the initial complex quadrature panels");
				using Interval = std::array<Real, 2>;
				std::vector<Interval> intervals;
				for (int i = 0; i < config.initial_panels; ++i) {
					const Real left = std::lerp(a, b, Real(i) / config.initial_panels);
					const Real right = std::lerp(a, b, Real(i + 1) / config.initial_panels);
					if (left == right) {
						result.status = AlgorithmStatus::ToleranceUnachievable;
						result.error_message = "Complex quadrature initial panels cannot make floating-point progress";
						return;
					}
					intervals.push_back({left, right});
				}
				IntegrationDetail::RequireFinite(function(a), "Complex contour endpoint");
				IntegrationDetail::RequireFinite(function(b), "Complex contour endpoint");
				auto evaluate = [&](const Interval& interval) { return Kronrod15(function, interval[0], interval[1]); };
				auto split = [](const Interval& interval) {
					const Real mid = std::midpoint(interval[0], interval[1]);
					if (mid == interval[0] || mid == interval[1]) return std::vector<Interval>{};
					return std::vector<Interval>{{interval[0], mid}, {mid, interval[1]}};
				};
				IntegrationDetail::RefineCells(result, intervals, evaluate, split, config.tolerance, 0,
					config.max_depth, config.max_evaluations, 15);
			}

			template<class Function, class Gamma, class Tangent>
			void ContourCore(Function& function, const Gamma& gamma, const Tangent& tangent, Real a, Real b,
				const ContourIntegrationConfig& config, ContourIntegrationResult& result) {
				auto integrand = [&](Real t) {
					const Complex point = gamma(t), derivative = tangent(t);
					ValidatePoint(point);
					IntegrationDetail::RequireFinite(derivative, "Complex contour tangent");
					return function(point) * derivative;
				};
				Integrate(integrand, a, b, config, result);
			}

			inline void Scale(ContourIntegrationResult& result, Real factor) noexcept {
				result.value = Complex(result.value.imag() * factor, -result.value.real() * factor);
				if (result.error_estimated) result.error_estimate *= std::abs(factor);
			}

			// Preserve output units and sampled derivative uncertainty even when a callback interrupts refinement.
			template<class Result>
			struct ScaledResult {
				Result& result;
				Real factor;
				~ScaledResult() {
					Scale(result, factor);
					if constexpr (requires { result.derivative_error_estimate; })
						if (result.error_estimated) result.error_estimate += result.derivative_error_estimate;
				}
			};

			inline ContourIntegrationConfig ScaledConfig(const ContourIntegrationConfig& config, Real factor) {
				ValidateConfig(config);
				auto scaled = config;
				scaled.tolerance = config.tolerance / std::abs(factor);
				ValidateConfig(scaled);
				return scaled;
			}

			inline void ValidateInterior(const CircleContour& contour, Complex point) {
				ValidatePoint(point);
				if (!(std::abs(point - contour.center()) < contour.radius()))
					throw ArgumentError("Cauchy's formula requires a point strictly inside the circle");
			}

			inline void CertifyInteger(IntegerContourResult& result, Real tolerance) {
				if (!result.IsUsable()) return;
				if (!std::isfinite(result.value.real()) || !std::isfinite(result.value.imag()) ||
					(result.error_estimated && (!std::isfinite(result.error_estimate) || result.error_estimate < 0))) {
					result.status = AlgorithmStatus::NumericalInstability;
					result.error_message = "Contour count or its uncertainty is non-finite or invalid";
					result.converged = false;
					return;
				}
				const Real nearest = std::round(result.value.real());
				if (!std::isfinite(nearest) || static_cast<long double>(nearest) < std::numeric_limits<int>::min() || // tol-ok: integer conversion bound
					static_cast<long double>(nearest) > std::numeric_limits<int>::max()) {
					result.status = AlgorithmStatus::NumericalInstability;
					result.error_message = "Contour count is outside the representable integer range";
					result.converged = false;
					return;
				}
				const Real discrepancy = std::abs(result.value - Complex(nearest, 0));
				if (!result.converged || !result.error_estimated ||
					discrepancy + result.error_estimate > std::min(tolerance, Precision::ComplexAnalysisIntegerRadius)) {
					result.status = AlgorithmStatus::ToleranceUnachievable;
					result.error_message = "Contour estimate, imaginary part and uncertainty do not certify a unique integer";
					result.converged = false;
					return;
				}
				result.integer_value = static_cast<int>(nearest);
				result.integer_certified = true;
			}
		}

		template<class Contour>
		static ContourIntegrationResult ContourIntegralDetailed(const IComplexFunction& f, const Contour& contour,
			const ContourIntegrationConfig& config = {}) {
			return IntegrationDetail::ExecuteIntegrationDetailed<ContourIntegrationResult>("ContourIntegral", config,
				[&](auto& result) {
					Detail::CountedFunction counted(f, result, config);
					auto tangent = [&](Real t) { return contour.derivative(t); };
					Detail::ContourCore(counted, contour, tangent, contour.t_start(), contour.t_end(), config, result);
				});
		}

		static inline ContourIntegrationResult ContourIntegralDetailed(const IComplexFunction& f,
			const std::function<Complex(Real)>& gamma, const std::function<Complex(Real)>& gamma_deriv,
			Real t_start, Real t_end, const ContourIntegrationConfig& config = {}) {
			return IntegrationDetail::ExecuteIntegrationDetailed<ContourIntegrationResult>("ContourIntegral", config,
				[&](auto& result) {
					if (!gamma || !gamma_deriv) throw ArgumentError("Contour parameterization and derivative must be callable");
					Detail::CountedFunction counted(f, result, config);
					Detail::ContourCore(counted, gamma, gamma_deriv, t_start, t_end, config, result);
				});
		}

		template<class Contour>
		static ContourIntegrationResult ContourIntegralDetailed(const IComplexFunction& f, const Contour& contour, Real tol) {
			return ContourIntegralDetailed(f, contour, Detail::WithTolerance(tol));
		}

		static inline ContourIntegrationResult ContourIntegralDetailed(const IComplexFunction& f,
			const std::function<Complex(Real)>& gamma, const std::function<Complex(Real)>& gamma_deriv,
			Real a, Real b, Real tol) {
			return ContourIntegralDetailed(f, gamma, gamma_deriv, a, b, Detail::WithTolerance(tol));
		}

		/// Compatibility API: historically returns a result, not a plain Complex.
		template<class Contour>
		static ContourIntegrationResult ContourIntegral(const IComplexFunction& f, const Contour& contour,
			Real tol = Precision::ComplexAnalysisTolerance) {
			return ContourIntegralDetailed(f, contour, tol);
		}

		static inline ContourIntegrationResult ContourIntegral(const IComplexFunction& f,
			const std::function<Complex(Real)>& gamma, const std::function<Complex(Real)>& gamma_deriv,
			Real a, Real b, Real tol = Precision::ComplexAnalysisTolerance) {
			return ContourIntegralDetailed(f, gamma, gamma_deriv, a, b, tol);
		}

		static inline ContourIntegrationResult CauchyDerivativeDetailed(const IComplexFunction& f,
			const CircleContour& contour, Complex z0, int n, const ContourIntegrationConfig& config = {}) {
			return IntegrationDetail::ExecuteIntegrationDetailed<ContourIntegrationResult>("CauchyDerivative", config,
				[&](auto& result) {
					Detail::ValidateInterior(contour, z0);
					if (n < 0) throw ArgumentError("Cauchy derivative order must be nonnegative");
					Real factorial = 1;
					for (int k = 2; k <= n; ++k) {
						if (factorial > std::numeric_limits<Real>::max() / k)
							throw ArgumentError("Cauchy derivative factorial exceeds the scalar range");
						factorial *= k;
					}
					const Real factor = factorial / (2 * Constants::PI);
					auto scaled = Detail::ScaledConfig(config, factor);
					scaled.initial_panels = std::max(scaled.initial_panels, n + 1);
					Detail::CountedFunction counted(f, result, config);
					auto function = [&](Complex z) {
						Complex denominator(1, 0);
						for (int k = 0; k <= n; ++k) denominator *= z - z0;
						IntegrationDetail::RequireFinite(denominator, "Cauchy denominator");
						if (denominator == Complex(0, 0)) throw NumericalMethodError("Cauchy denominator underflowed");
						return counted(z) / denominator;
					};
					auto tangent = [&](Real t) { return contour.derivative(t); };
					Detail::ScaledResult scale{result, factor};
					Detail::ContourCore(function, contour, tangent, contour.t_start(), contour.t_end(), scaled, result);
				});
		}

		static inline ContourIntegrationResult CauchyDerivativeDetailed(const IComplexFunction& f,
			const CircleContour& contour, Complex z0, int n, Real tol) {
			return CauchyDerivativeDetailed(f, contour, z0, n, Detail::WithTolerance(tol));
		}

		static inline Complex CauchyDerivative(const IComplexFunction& f, const CircleContour& contour, Complex z0,
			int n, const ContourIntegrationConfig& config) {
			return CauchyDerivativeDetailed(f, contour, z0, n, config).ValueOrThrow("CauchyDerivative");
		}

		static inline Complex CauchyDerivative(const IComplexFunction& f, const CircleContour& contour, Complex z0,
			int n, Real tol = Precision::ComplexAnalysisTolerance) {
			return CauchyDerivative(f, contour, z0, n, Detail::WithTolerance(tol));
		}

		static inline ContourIntegrationResult CauchyIntegralFormulaDetailed(const IComplexFunction& f,
			const CircleContour& contour, Complex z0, const ContourIntegrationConfig& config = {}) {
			auto result = CauchyDerivativeDetailed(f, contour, z0, 0, config);
			result.algorithm_name = "CauchyIntegralFormula";
			return result;
		}

		static inline ContourIntegrationResult CauchyIntegralFormulaDetailed(const IComplexFunction& f,
			const CircleContour& contour, Complex z0, Real tol) {
			return CauchyIntegralFormulaDetailed(f, contour, z0, Detail::WithTolerance(tol));
		}

		static inline Complex CauchyIntegralFormula(const IComplexFunction& f, const CircleContour& contour, Complex z0,
			const ContourIntegrationConfig& config) {
			return CauchyIntegralFormulaDetailed(f, contour, z0, config).ValueOrThrow("CauchyIntegralFormula");
		}

		static inline Complex CauchyIntegralFormula(const IComplexFunction& f, const CircleContour& contour, Complex z0,
			Real tol = Precision::ComplexAnalysisTolerance) {
			return CauchyIntegralFormula(f, contour, z0, Detail::WithTolerance(tol));
		}

		static inline ContourIntegrationResult ResidueDetailed(const IComplexFunction& f, Complex z0, Real radius,
			const ContourIntegrationConfig& config) {
			return IntegrationDetail::ExecuteIntegrationDetailed<ContourIntegrationResult>("Residue", config,
				[&](auto& result) {
					CircleContour contour(z0, radius);
					const Real factor = 1 / (2 * Constants::PI);
					const auto scaled = Detail::ScaledConfig(config, factor);
					Detail::CountedFunction counted(f, result, config);
					auto tangent = [&](Real t) { return contour.derivative(t); };
					Detail::ScaledResult scale{result, factor};
					Detail::ContourCore(counted, contour, tangent, contour.t_start(), contour.t_end(), scaled, result);
				});
		}

		static inline ContourIntegrationResult ResidueDetailed(const IComplexFunction& f, Complex z0,
			Real radius = Real(0.1), Real tol = Precision::DefaultToleranceStrict) {
			return ResidueDetailed(f, z0, radius, Detail::WithTolerance(tol));
		}

		static inline Complex Residue(const IComplexFunction& f, Complex z0, Real radius,
			const ContourIntegrationConfig& config) {
			return ResidueDetailed(f, z0, radius, config).ValueOrThrow("Residue");
		}

		static inline Complex Residue(const IComplexFunction& f, Complex z0, Real radius = Real(0.1),
			Real tol = Precision::DefaultToleranceStrict) {
			return Residue(f, z0, radius, Detail::WithTolerance(tol));
		}

		static inline ContourIntegrationResult ResidueSimplePoleDetailed(const IComplexFunction& f, Complex z0,
			Real h = Precision::DefaultToleranceStrict, const ContourIntegrationConfig& config = {}) {
			return IntegrationDetail::ExecuteIntegrationDetailed<ContourIntegrationResult>("ResidueSimplePole", config,
				[&](auto& result) {
					Detail::ValidateConfig(config);
					Detail::ValidatePoint(z0);
					if (!std::isfinite(h) || h <= 0) throw ArgumentError("Simple-pole step must be finite and positive");
					const std::array<Complex, 4> points{z0 + Complex(h, 0), z0 - Complex(h, 0),
						z0 + Complex(0, h), z0 - Complex(0, h)};
					for (Complex point : points) {
						Detail::ValidatePoint(point);
						if (point == z0) throw ArgumentError("Simple-pole step cannot make floating-point progress");
					}
					Detail::CountedFunction counted(f, result, config);
					IntegrationDetail::ComplexSum sum;
					for (int axis = 0; axis < 2; ++axis) {
						const Complex plus = points[2 * axis] - z0, minus = points[2 * axis + 1] - z0;
						// Actual displacements may differ at a binade boundary: interpolate each opposite pair to zero.
						const Real denominator = std::abs(plus) / 2 + std::abs(minus) / 2;
						const Real weight = (std::abs(minus) / 2) / denominator;
						sum.Add((plus * counted(points[2 * axis])) * (weight / 2));
						sum.Add((minus * counted(points[2 * axis + 1])) * ((1 - weight) / 2));
					}
					result.value = sum.Value();
					result.MarkFixedRule(result.function_evaluations);
				});
		}

		static inline Complex ResidueSimplePole(const IComplexFunction& f, Complex z0, Real h = Precision::DefaultToleranceStrict) {
			return ResidueSimplePoleDetailed(f, z0, h).ValueOrThrow("ResidueSimplePole");
		}

		static inline IntegerContourResult WindingNumberDetailed(const CircleContour& contour, Complex z0,
			const ContourIntegrationConfig& config = {}) {
			return IntegrationDetail::ExecuteIntegrationDetailed<IntegerContourResult>("WindingNumber", config,
				[&](auto& result) {
					Detail::ValidatePoint(z0);
					if (std::abs(z0 - contour.center()) == contour.radius())
						throw ArgumentError("Winding number is undefined for a point on the contour");
					const Real factor = 1 / (2 * Constants::PI);
					const auto scaled = Detail::ScaledConfig(config, factor);
					ComplexFunctionFromStdFunc reciprocal([&](Complex z) {
						if (z == z0) throw NumericalMethodError("Winding-number contour intersects its pole");
						return Real(1) / (z - z0);
					});
					Detail::CountedFunction counted(reciprocal, result, config);
					auto tangent = [&](Real t) { return contour.derivative(t); };
					{
						Detail::ScaledResult scale{result, factor};
						Detail::ContourCore(counted, contour, tangent, contour.t_start(), contour.t_end(), scaled, result);
					}
					Detail::CertifyInteger(result, config.tolerance);
				});
		}

		static inline IntegerContourResult WindingNumberDetailed(const CircleContour& contour, Complex z0, Real tol) {
			return WindingNumberDetailed(contour, z0, Detail::WithTolerance(tol));
		}

		static inline int WindingNumber(const CircleContour& contour, Complex z0, const ContourIntegrationConfig& config) {
			return WindingNumberDetailed(contour, z0, config).IntegerOrThrow("WindingNumber");
		}

		static inline int WindingNumber(const CircleContour& contour, Complex z0, Real tol = Precision::ComplexAnalysisTolerance) {
			return WindingNumber(contour, z0, Detail::WithTolerance(tol));
		}

		static inline IntegerContourResult ArgumentPrincipleDetailed(const IComplexFunction& f, const CircleContour& contour,
			const ContourIntegrationConfig& config) {
			return IntegrationDetail::ExecuteIntegrationDetailed<IntegerContourResult>("ArgumentPrinciple", config,
				[&](auto& result) {
					const Real factor = 1 / (2 * Constants::PI);
					auto scaled = Detail::ScaledConfig(config, factor);
					scaled.tolerance /= 2;
					Detail::CountedFunction counted(f, result, config);
					auto logarithmicDerivative = [&](Complex z) {
						const Complex value = counted(z);
						if (value == Complex(0, 0)) throw NumericalMethodError("Argument-principle function is zero on the contour");
						DerivativeConfig derivativeConfig;
						derivativeConfig.estimate_error = true;
						Real h = Derivation::Detail::ScaleStepComplex(Derivation::NDer4_h, z);
						for (int refinement = 0; refinement <= Precision::ComplexAnalysisDerivativeRefinements; ++refinement) {
							const auto derivative = Derivation::NDer4ComplexDetailed(counted, z, h, derivativeConfig);
							if (!MML::IsUsable(derivative.status)) throw NumericalMethodError(derivative.error_message);
							const Real uncertainty = IntegrationDetail::Product({derivative.error, contour.radius()}, std::abs(value));
							if (uncertainty <= config.tolerance / 4) {
								result.derivative_error_estimate = std::max(result.derivative_error_estimate, uncertainty);
								return derivative.value / value;
							}
							h /= 2;
						}
						throw NumericalMethodError("Argument-principle numerical derivative cannot resolve the requested count tolerance");
					};
					auto tangent = [&](Real t) { return contour.derivative(t); };
					{
						Detail::ScaledResult scale{result, factor};
						Detail::ContourCore(logarithmicDerivative, contour, tangent, contour.t_start(), contour.t_end(), scaled, result);
					}
					Detail::CertifyInteger(result, config.tolerance);
				});
		}

		static inline IntegerContourResult ArgumentPrincipleDetailed(const IComplexFunction& f, const CircleContour& contour,
			Real tol = Precision::DefaultToleranceStrict) {
			return ArgumentPrincipleDetailed(f, contour, Detail::WithTolerance(tol));
		}

		static inline int ArgumentPrinciple(const IComplexFunction& f, const CircleContour& contour,
			const ContourIntegrationConfig& config) {
			return ArgumentPrincipleDetailed(f, contour, config).IntegerOrThrow("ArgumentPrinciple");
		}

		static inline int ArgumentPrinciple(const IComplexFunction& f, const CircleContour& contour,
			Real tol = Precision::DefaultToleranceStrict) {
			return ArgumentPrinciple(f, contour, Detail::WithTolerance(tol));
		}

		static inline IntegerContourResult CountZerosDetailed(const IComplexFunction& f, const CircleContour& contour,
			const ContourIntegrationConfig& config) {
			auto result = ArgumentPrincipleDetailed(f, contour, config);
			result.algorithm_name = "CountZeros";
			return result;
		}

		static inline IntegerContourResult CountZerosDetailed(const IComplexFunction& f, const CircleContour& contour,
			Real tol = Precision::DefaultToleranceStrict) {
			return CountZerosDetailed(f, contour, Detail::WithTolerance(tol));
		}

		static inline int CountZeros(const IComplexFunction& f, const CircleContour& contour,
			const ContourIntegrationConfig& config) {
			return CountZerosDetailed(f, contour, config).IntegerOrThrow("CountZeros");
		}

		static inline int CountZeros(const IComplexFunction& f, const CircleContour& contour,
			Real tol = Precision::DefaultToleranceStrict) {
			return CountZeros(f, contour, Detail::WithTolerance(tol));
		}
	}
}

#endif // MML_COMPLEX_ANALYSIS_H
