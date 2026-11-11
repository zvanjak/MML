#include "../TestPrecision.h"
#ifndef MML_SINGLE_HEADER_INCLUDE
#include <mml/core/Surfaces.h>
#endif
#include <catch2/catch_all.hpp>
#include <algorithm>
#include <array>
#include <vector>

using namespace MML;

namespace MML::Tests::Core::SurfaceStepContracts {
	using Point = std::pair<Real, Real>;

	class RecordedSurface : public IParametricSurfaceRect<3> {
		bool bounded;
	public:
		mutable std::vector<Point> samples;
		explicit RecordedSurface(bool bounded = false) : bounded(bounded) {}
		Real getMinU() const override { return bounded ? Real(0) : Constants::NegInf; }
		Real getMaxU() const override { return bounded ? Real(1) : Constants::PosInf; }
		Real getMinW() const override { return getMinU(); }
		Real getMaxW() const override { return getMaxU(); }
		VectorN<Real, 3> operator()(Real u, Real w) const override {
			samples.emplace_back(u, w);
			if (bounded && (u < 0 || u > 1 || w < 0 || w > 1))
				throw DomainError("Recorded surface sampled outside its chart");
			return {std::cos(u) * std::sin(w), std::sin(u) * std::cos(w), std::exp(u + w)};
		}
	};

	void CheckSamples(const RecordedSurface& surface, std::vector<Point> expected) {
		auto actual = surface.samples;
		REQUIRE(actual.size() == expected.size());
		std::sort(actual.begin(), actual.end());
		std::sort(expected.begin(), expected.end());
		for (size_t i = 0; i < actual.size(); ++i) {
			CAPTURE(i);
			CHECK(actual[i].first == expected[i].first);
			CHECK(actual[i].second == expected[i].second);
		}
	}

	TEST_CASE("Surface second partial defaults use the existing fourth-root scale", "[surface-steps]") {
		RecordedSurface surface;
		const Real u = 2, w = -3;
		const Real hu = Derivation::ScaleStep(Precision::FourthRootEps, u);
		const Real hw = Derivation::ScaleStep(Precision::FourthRootEps, w);
		Real error;
		Derivation::NDer2_uu(surface, u, w, &error);
		CheckSamples(surface, {{u + hu, w}, {u - hu, w}, {u, w}, {u + 2 * hu, w}, {u - 2 * hu, w}});
		surface.samples.clear();
		Derivation::NDer2_ww(surface, u, w, &error);
		CheckSamples(surface, {{u, w + hw}, {u, w - hw}, {u, w}, {u, w + 2 * hw}, {u, w - 2 * hw}});
		surface.samples.clear();
		const Real h = Derivation::ScaleStep(Precision::FourthRootEps, std::max(std::abs(u), std::abs(w)));
		Derivation::NDer2_uw(surface, u, w, &error);
		CheckSamples(surface, {{u + h, w + h}, {u + h, w - h}, {u - h, w + h}, {u - h, w - h},
			{u + h / 2, w + h / 2}, {u + h / 2, w - h / 2}, {u - h / 2, w + h / 2}, {u - h / 2, w - h / 2}});
	}

	TEST_CASE("Surface first partial defaults retain their original scales", "[surface-steps]") {
		RecordedSurface surface;
		const Real u = Real(0.3), w = Real(1.7);
		const Real forward = 2 * std::sqrt(Constants::Eps);
		const Real centered = std::pow(3 * Constants::Eps, Real(1) / Real(3));
		CHECK(Derivation::NDer1_h == forward);
		CHECK(Derivation::NDer2_h == centered);
		CHECK((Derivation::NDer1_u(surface, u, w) -
			Derivation::NDer1_u(surface, u, w, Derivation::ScaleStep(forward, u))).NormL2() == 0);
		CHECK((Derivation::NDer1_w(surface, u, w) -
			Derivation::NDer1_w(surface, u, w, Derivation::ScaleStep(forward, w))).NormL2() == 0);
		CHECK((Derivation::NDer2_u(surface, u, w) -
			Derivation::NDer2_u(surface, u, w, Derivation::ScaleStep(centered, u))).NormL2() == 0);
		CHECK((Derivation::NDer2_w(surface, u, w) -
			Derivation::NDer2_w(surface, u, w, Derivation::ScaleStep(centered, w))).NormL2() == 0);
	}

	TEST_CASE("Default surface partials and curvatures meet a square-root-epsilon budget", "[surface-steps]") {
		const Real budget = 16 * Precision::SqrtEps;
		for (Real radius : {Real(1e-8), Real(1e-4), Real(1), Real(1e4), Real(1e8)})
			for (Real u : {Real(0.3), Real(0.9), Real(1.7), Real(2.5)})
				for (Real w : {Real(-0.9), Real(0.2), Real(2.3)}) {
					CAPTURE(radius, u, w);
					Surfaces::Sphere sphere(radius);
					const Real su = std::sin(u), cu = std::cos(u), sw = std::sin(w), cw = std::cos(w);
					const std::array<VectorN<Real, 3>, 3> expected = {{
						{-su * cw, -su * sw, -cu}, {-cu * sw, cu * cw, 0}, {-su * cw, -su * sw, 0}
					}};
					std::array<Real, 3> estimates;
					const std::array<VectorN<Real, 3>, 3> actual = {{
						Derivation::NDer2_uu(sphere, u, w, &estimates[0]),
						Derivation::NDer2_uw(sphere, u, w, &estimates[1]),
						Derivation::NDer2_ww(sphere, u, w, &estimates[2])
					}};
					for (size_t partial = 0; partial < actual.size(); ++partial) {
						CAPTURE(partial);
						for (int component = 0; component < 3; ++component)
							CHECK(std::abs(actual[partial][component] / radius - expected[partial][component]) <= budget);
						CHECK(std::isfinite(estimates[partial]));
						CHECK(estimates[partial] >= 0);
						CHECK(estimates[partial] / radius <= budget);
					}
					CHECK(std::abs(sphere.GaussianCurvature(u, w) * radius * radius - 1) <= budget);
					CHECK(std::abs(sphere.MeanCurvature(u, w) * radius + 1) <= budget);
				}
	}

	TEST_CASE("Explicit surface steps retain second-order convergence and error controls", "[surface-steps]") {
		RecordedSurface surface;
		const Real u = Real(0.3), w = Real(0.2);
		const Real x = std::cos(u) * std::sin(w), y = std::sin(u) * std::cos(w), z = std::exp(u + w);
		const std::array<VectorN<Real, 3>, 3> expected = {{
			{-x, -y, z}, {-std::sin(u) * std::cos(w), -std::cos(u) * std::sin(w), z}, {-x, -y, z}
		}};
		std::array<Real, 3> previous{};
		for (int refinement = 0; refinement < 3; ++refinement) {
			const Real h = std::ldexp(Real(1) / 4, -refinement);
			CAPTURE(h);
			std::array<Real, 3> estimates;
			const std::array<VectorN<Real, 3>, 3> estimated = {{
				Derivation::NDer2_uu(surface, u, w, h, &estimates[0]),
				Derivation::NDer2_uw(surface, u, w, h, &estimates[1]),
				Derivation::NDer2_ww(surface, u, w, h, &estimates[2])
			}};
			const std::array<VectorN<Real, 3>, 3> plain = {{
				Derivation::NDer2_uu(surface, u, w, h),
				Derivation::NDer2_uw(surface, u, w, h),
				Derivation::NDer2_ww(surface, u, w, h)
			}};
			for (size_t partial = 0; partial < expected.size(); ++partial) {
				CAPTURE(partial);
				const Real actualError = (estimated[partial] - expected[partial]).NormL2();
				CHECK((estimated[partial] - plain[partial]).NormL2() == 0);
				CHECK(actualError > 0);
				CHECK(estimates[partial] >= actualError / 2);
				CHECK(estimates[partial] <= 2 * actualError);
				if (refinement != 0) {
					CHECK(previous[partial] / actualError > 3);
					CHECK(previous[partial] / actualError < 5);
				}
				previous[partial] = actualError;
			}
		}
	}

	TEST_CASE("Surface defaults preserve centered chart-boundary contracts", "[surface-steps]") {
		RecordedSurface surface(true);
		const Real h = Precision::FourthRootEps;
		CHECK_THROWS_AS(Derivation::NDer2_uu(surface, Real(0), Real(0.5)), DomainError);
		CHECK_THROWS_AS(Derivation::NDer2_ww(surface, Real(0.5), Real(0)), DomainError);
		CHECK_THROWS_AS(Derivation::NDer2_uw(surface, Real(1), Real(1)), DomainError);
		CHECK_NOTHROW(Derivation::NDer2_uu(surface, Real(1.5) * h, Real(0.5)));
		Real error;
		CHECK_THROWS_AS(Derivation::NDer2_uu(surface, Real(1.5) * h, Real(0.5), &error), DomainError);
		surface.samples.clear();
		CHECK_NOTHROW(Derivation::NDer2_uw(surface, h, h, &error));
		CHECK(surface.samples.size() == 8);
	}

	TEST_CASE("Anisotropic explicit surface steps keep their independent sample spans", "[surface-steps]") {
		RecordedSurface surface;
		const Real u = Real(0.3), w = Real(0.2), hu = Real(1) / 8, hw = Real(1) / 32;
		Real error;
		const auto estimated = Derivation::NDer2_uw(surface, u, w, hu, hw, &error);
		CheckSamples(surface, {{u + hu, w + hw}, {u + hu, w - hw}, {u - hu, w + hw}, {u - hu, w - hw},
			{u + hu / 2, w + hw / 2}, {u + hu / 2, w - hw / 2}, {u - hu / 2, w + hw / 2}, {u - hu / 2, w - hw / 2}});
		surface.samples.clear();
		const auto plain = Derivation::NDer2_uw(surface, u, w, hu, hw);
		CHECK(surface.samples.size() == 4);
		CHECK((estimated - plain).NormL2() == 0);
		CHECK(std::isfinite(error));
		CHECK(error > 0);
	}
}
