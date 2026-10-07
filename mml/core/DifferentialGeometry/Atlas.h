///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        DifferentialGeometry/Atlas.h                                        ///
///  Description: Minimal atlas support for homogeneous chart collections             ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CORE_DIFFERENTIAL_GEOMETRY_ATLAS_H
#define MML_CORE_DIFFERENTIAL_GEOMETRY_ATLAS_H

#include <mml/MMLBase.h>

#include <mml/core/DifferentialGeometry/Chart.h>

#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace MML::DifferentialGeometry
{
	template<int DomainN, int AmbientN, class DomainFrame, class AmbientFrame>
	class Atlas
	{
	public:
		using ChartType = Chart<DomainN, AmbientN, DomainFrame, AmbientFrame>;
		using DomainPoint = Point<Real, DomainN, DomainFrame>;
		using AmbientPoint = Point<Real, AmbientN, AmbientFrame>;
		using DomainPredicate = std::function<bool(const DomainPoint&)>;
		using AmbientPredicate = std::function<bool(const AmbientPoint&)>;
		using TransitionFunction = std::function<DomainPoint(const DomainPoint&)>;

	private:
		struct ChartEntry
		{
			std::string name;
			std::shared_ptr<const ChartType> chart;
			DomainPredicate containsDomain;
			AmbientPredicate containsAmbient;
		};

		struct TransitionEntry
		{
			int fromChart = -1;
			int toChart = -1;
			TransitionFunction transition;
		};

		std::vector<ChartEntry> _charts;
		std::vector<TransitionEntry> _transitions;

		void ValidateChartIndex(int chartIndex, const char* functionName) const
		{
			if (chartIndex < 0 || chartIndex >= static_cast<int>(_charts.size()))
				throw IndexError(std::string(functionName) + " - chart index out of range");
		}

	public:
		int AddChart(std::string name,
			std::shared_ptr<const ChartType> chart,
			DomainPredicate containsDomain = {},
			AmbientPredicate containsAmbient = {})
		{
			if (!chart)
				throw ArgumentError("Atlas::AddChart - chart cannot be null");

			_charts.push_back({ std::move(name), std::move(chart), std::move(containsDomain), std::move(containsAmbient) });
			return static_cast<int>(_charts.size()) - 1;
		}

		void AddTransition(int fromChart, int toChart, TransitionFunction transition)
		{
			ValidateChartIndex(fromChart, "Atlas::AddTransition");
			ValidateChartIndex(toChart, "Atlas::AddTransition");
			if (!transition)
				throw ArgumentError("Atlas::AddTransition - transition cannot be empty");

			_transitions.push_back({ fromChart, toChart, std::move(transition) });
		}

		int chartCount() const noexcept { return static_cast<int>(_charts.size()); }

		const std::string& chartName(int chartIndex) const
		{
			ValidateChartIndex(chartIndex, "Atlas::chartName");
			return _charts[chartIndex].name;
		}

		const ChartType& chart(int chartIndex) const
		{
			ValidateChartIndex(chartIndex, "Atlas::chart");
			return *_charts[chartIndex].chart;
		}

		bool contains(int chartIndex, const DomainPoint& point) const
		{
			ValidateChartIndex(chartIndex, "Atlas::contains");
			if (!_charts[chartIndex].containsDomain)
				return true;
			return _charts[chartIndex].containsDomain(point);
		}

		bool containsAmbient(int chartIndex, const AmbientPoint& point) const
		{
			ValidateChartIndex(chartIndex, "Atlas::containsAmbient");
			if (!_charts[chartIndex].containsAmbient)
				return false;
			return _charts[chartIndex].containsAmbient(point);
		}

		int selectChartForAmbientPoint(const AmbientPoint& point) const
		{
			for (int i = 0; i < static_cast<int>(_charts.size()); i++)
				if (containsAmbient(i, point))
					return i;
			return -1;
		}

		DomainPoint transition(int fromChart, int toChart, const DomainPoint& point) const
		{
			ValidateChartIndex(fromChart, "Atlas::transition");
			ValidateChartIndex(toChart, "Atlas::transition");
			if (fromChart == toChart)
				return point;

			for (const auto& entry : _transitions) {
				if (entry.fromChart == fromChart && entry.toChart == toChart)
					return entry.transition(point);
			}

			throw ArgumentError("Atlas::transition - transition map is not registered");
		}
	};

	class UnitSphereStereographicChart : public Chart<2, 3, Parameter2, Cartesian3>
	{
	public:
		enum class Pole
		{
			North,
			South
		};

	private:
		Pole _pole;

	public:
		explicit UnitSphereStereographicChart(Pole pole) : _pole(pole) { }

		Pole pole() const noexcept { return _pole; }

		Point<Real, 3, Cartesian3> map_point(const Point<Real, 2, Parameter2>& point) const override
		{
			Real u = point[0];
			Real v = point[1];
			Real radiusSquared = u * u + v * v;
			Real denominator = REAL(1.0) + radiusSquared;
			Real z = _pole == Pole::North
				? (REAL(1.0) - radiusSquared) / denominator
				: (-REAL(1.0) + radiusSquared) / denominator;

			return Point<Real, 3, Cartesian3>{ REAL(2.0) * u / denominator, REAL(2.0) * v / denominator, z };
		}

		MatrixNM<Real, 3, 2> jacobian(const Point<Real, 2, Parameter2>& point) const override
		{
			Real u = point[0];
			Real v = point[1];
			Real radiusSquared = u * u + v * v;
			Real denominator = REAL(1.0) + radiusSquared;
			Real denominatorSquared = denominator * denominator;
			Real zSign = _pole == Pole::North ? -REAL(4.0) : REAL(4.0);

			MatrixNM<Real, 3, 2> jac;
			jac(0, 0) = REAL(2.0) * (REAL(1.0) - u * u + v * v) / denominatorSquared;
			jac(0, 1) = -REAL(4.0) * u * v / denominatorSquared;
			jac(1, 0) = -REAL(4.0) * u * v / denominatorSquared;
			jac(1, 1) = REAL(2.0) * (REAL(1.0) + u * u - v * v) / denominatorSquared;
			jac(2, 0) = zSign * u / denominatorSquared;
			jac(2, 1) = zSign * v / denominatorSquared;
			return jac;
		}
	};

	inline bool IsFiniteStereographicPoint(const Point<Real, 2, Parameter2>& point)
	{
		return std::isfinite(point[0]) && std::isfinite(point[1]);
	}

	inline Point<Real, 2, Parameter2> StereographicSphereTransition(const Point<Real, 2, Parameter2>& point)
	{
		Real radiusSquared = point[0] * point[0] + point[1] * point[1];
		if (radiusSquared <= PrecisionValues<Real>::NumericalZeroThreshold)
			throw ArgumentError("StereographicSphereTransition - pole point is not in chart overlap");

		return Point<Real, 2, Parameter2>{ point[0] / radiusSquared, point[1] / radiusSquared };
	}

	inline Atlas<2, 3, Parameter2, Cartesian3> MakeUnitSphereStereographicAtlas()
	{
		Atlas<2, 3, Parameter2, Cartesian3> atlas;
		auto north = std::make_shared<UnitSphereStereographicChart>(UnitSphereStereographicChart::Pole::North);
		auto south = std::make_shared<UnitSphereStereographicChart>(UnitSphereStereographicChart::Pole::South);

		int northId = atlas.AddChart("north", north, IsFiniteStereographicPoint,
			[](const Point<Real, 3, Cartesian3>& point) { return point[2] >= REAL(0.0); });
		int southId = atlas.AddChart("south", south, IsFiniteStereographicPoint,
			[](const Point<Real, 3, Cartesian3>& point) { return point[2] < REAL(0.0); });

		atlas.AddTransition(northId, southId, StereographicSphereTransition);
		atlas.AddTransition(southId, northId, StereographicSphereTransition);
		return atlas;
	}
} // namespace MML::DifferentialGeometry

#endif // MML_CORE_DIFFERENTIAL_GEOMETRY_ATLAS_H