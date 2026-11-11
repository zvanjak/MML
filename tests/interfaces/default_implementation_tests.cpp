#include <catch2/catch_all.hpp>
#include "../TestPrecision.h"

#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/base/Intervals.h>
#include <mml/interfaces/IFunction.h>
#include <mml/systems/DynamicalSystem.h>
#endif

using namespace MML;
using namespace MML::Systems;

namespace MML::Tests::Interfaces::DefaultImplementationTests {
	namespace {
		class NumericalSystem : public DynamicalSystemBase<3, 0> {
		public:
			void derivs(Real t, const Vector<Real>& x, Vector<Real>& f) const override {
				f[0] = x[0] * x[1] + t * x[2];
				f[1] = std::sin(x[0]) - 3 * x[1];
				f[2] = x[2] * x[2] + x[0] - t;
			}
		};

		class TraceCurve : public IParametricCurve<3> {
		public:
			mutable std::vector<Real> evaluated;
			Real getMinT() const override { return -2; }
			Real getMaxT() const override { return 4; }
			VectorN<Real, 3> operator()(Real t) const override {
				evaluated.push_back(t);
				return {t, t * t, 2 * t + 1};
			}
		};
	} // namespace

	TEST_CASE("Interface defaults - numerical Jacobian retains state and differentiates every column",
						"[interfaces][adversarial][wewi.2.7]") {
		NumericalSystem system;
		REQUIRE_FALSE(system.hasAnalyticalJacobian());
		const Real time = Real(0.75);
		for (const auto& x : {Vector<Real>{0, 0, 0}, Vector<Real>{Real(0.25), -2, 3}, Vector<Real>{-8, Real(0.001), -4}}) {
			const auto original = x;
			Matrix<Real> jacobian(1, 1);
			system.jacobian(time, x, jacobian);
			REQUIRE(x == original);
			REQUIRE(jacobian.rows() == 3);
			REQUIRE(jacobian.cols() == 3);
			const Matrix<Real> expected(3, 3, {x[1], x[0], time, std::cos(x[0]), -3, 0, 1, 0, 2 * x[2]});
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j) {
					CAPTURE(x, i, j);
					REQUIRE(IsWithinAbs(jacobian(i, j), expected(i, j), Precision::CoarseTolerance));
				}
		}
	}

	TEST_CASE("Interface defaults - zero parameters round trip without per-index access",
						"[interfaces][parametrized][adversarial][wewi.2.7]") {
		class NoParameters : public IParametrized {
		public:
			int getNumParam() const override { return 0; }
			Real getParam(int) const override {
				FAIL("getParam must not be called");
				return 0;
			}
			void setParam(int, Real) override { FAIL("setParam must not be called"); }
		} parameters;
		IParametrized& view = parameters;
		REQUIRE(view.getParams().size() == 0);
		REQUIRE_NOTHROW(view.setParams(Vector<Real>{}));
		REQUIRE_THROWS_AS(view.setParams(Vector<Real>{1}), VectorDimensionError);
	}

	TEST_CASE("Interface defaults - curve traces handle empty singleton and endpoint samples", "[interfaces][adversarial][wewi.2.7]") {
		for (int count : {-1, 0, 1, 2, 5}) {
			for (const auto& limits : {std::pair<Real, Real>{Real(-2), Real(4)}, {Real(4), Real(-2)}, {Real(1), Real(1)}}) {
				CAPTURE(count, limits.first, limits.second);
				TraceCurve curve;
				const IParametricCurve<3>& view = curve;
				auto trace = view.getTrace(limits.first, limits.second, count);
				REQUIRE(trace.size() == static_cast<size_t>(std::max(0, count)));
				REQUIRE(curve.evaluated.size() == trace.size());
				for (size_t i = 0; i < trace.size(); ++i) {
					const Real t = count == 1 ? limits.first : limits.first + Real(i) * (limits.second - limits.first) / Real(count - 1);
					REQUIRE(trace[i][0] == t);
					REQUIRE(trace[i][1] == t * t);
					REQUIRE(trace[i][2] == 2 * t + 1);
				}
			}
		}
	}

	TEST_CASE("Interface defaults - interval covering edge counts replace old output and respect endpoints",
						"[interfaces][intervals][adversarial][wewi.2.7]") {
		ClosedInterval closed(-2, 4);
		OpenInterval open(-2, 4);
		CompleteRInterval unbounded;
		for (const IInterval* interval :
				 {static_cast<const IInterval*>(&closed), static_cast<const IInterval*>(&open), static_cast<const IInterval*>(&unbounded)}) {
			for (int count : {-1, 0, 1, 2, 5}) {
				CAPTURE(count, interval->getLowerBound(), interval->getUpperBound());
				std::vector<Real> points{99, 100};
				interval->getEquidistantCovering(count, points);
				REQUIRE(points.size() == static_cast<size_t>(std::max(0, count)));
				for (size_t i = 0; i < points.size(); ++i) {
					REQUIRE(std::isfinite(points[i]));
					REQUIRE(interval->contains(points[i]));
					if (i > 0)
						REQUIRE(points[i] > points[i - 1]);
				}
				if (count == 1)
					REQUIRE(points[0] == (interval == &unbounded ? 0 : 1));
				if (count >= 2 && interval == &closed) {
					REQUIRE(points.front() == -2);
					REQUIRE(points.back() == 4);
				}
			}
		}
		ClosedInterval point(3, 3);
		std::vector<Real> points;
		point.getEquidistantCovering(2, points);
		REQUIRE(points == std::vector<Real>{3});
	}
} // namespace MML::Tests::Interfaces::DefaultImplementationTests
