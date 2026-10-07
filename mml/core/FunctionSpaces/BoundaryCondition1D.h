///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionSpaces/BoundaryCondition1D.h                                ///
///  Description: One-dimensional boundary condition descriptors                      ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_BOUNDARY_CONDITION_1D_H
#define MML_BOUNDARY_CONDITION_1D_H

#include <mml/core/FunctionSpaces/FunctionSpace1D.h>
#include <mml/core/FunctionSpaces/TrialSpace1D.h>

#include <mml/base/Vector/Vector.h>

#include <cmath>
#include <vector>

namespace MML::FunctionSpaces
{
	enum class BoundaryConditionKind
	{
		Dirichlet,
		Neumann,
		Robin,
		Periodic
	};

	struct BoundaryCondition1D
	{
		BoundaryConditionKind kind = BoundaryConditionKind::Dirichlet;
		Real location = REAL(0.0);
		Real alpha = REAL(1.0);
		Real beta = REAL(0.0);
		Real value = REAL(0.0);

		static BoundaryCondition1D Dirichlet(Real location, Real value)
		{
			return BoundaryCondition1D{ BoundaryConditionKind::Dirichlet, location, REAL(1.0), REAL(0.0), value };
		}

		static BoundaryCondition1D Neumann(Real location, Real value)
		{
			return BoundaryCondition1D{ BoundaryConditionKind::Neumann, location, REAL(0.0), REAL(1.0), value };
		}

		static BoundaryCondition1D Robin(Real location, Real alpha, Real beta, Real value)
		{
			return BoundaryCondition1D{ BoundaryConditionKind::Robin, location, alpha, beta, value };
		}

		static BoundaryCondition1D Periodic()
		{
			return BoundaryCondition1D{ BoundaryConditionKind::Periodic, REAL(0.0), REAL(1.0), REAL(0.0), REAL(0.0) };
		}

		void validate(const FunctionSpace1D& space, Real tolerance = Defaults::VectorIsEqualTolerance) const
		{
			if (kind != BoundaryConditionKind::Periodic && !space.contains(location, tolerance))
				throw FunctionSpaceInputError("BoundaryCondition1D: boundary location must be inside function-space domain");
			if (!std::isfinite(value) || !std::isfinite(alpha) || !std::isfinite(beta))
				throw FunctionSpaceInputError("BoundaryCondition1D: alpha, beta, and value must be finite");
			if (kind == BoundaryConditionKind::Robin && std::abs(alpha) <= tolerance && std::abs(beta) <= tolerance)
				throw FunctionSpaceInputError("BoundaryCondition1D: Robin condition requires non-zero alpha or beta");
		}

		Vector<Real> row(const TrialSpace1D& trialSpace, Real tolerance = Defaults::VectorIsEqualTolerance) const
		{
			validate(trialSpace.functionSpace(), tolerance);
			if (kind == BoundaryConditionKind::Periodic)
				throw FunctionSpaceError("BoundaryCondition1D::row: periodic boundary rows require paired assembly");

			Vector<Real> result(trialSpace.dimension());
			for (int index = 0; index < trialSpace.dimension(); ++index) {
				Real basis = trialSpace.basisValue(index, location);
				Real derivative = REAL(0.0);
				if (kind == BoundaryConditionKind::Neumann || kind == BoundaryConditionKind::Robin) {
					if (!trialSpace.hasBasisDerivative(1))
						throw FunctionSpaceError("BoundaryCondition1D::row: trial space does not provide first derivatives");
					derivative = trialSpace.basisDerivative(index, 1, location);
				}
				result[index] = alpha * basis + beta * derivative;
			}
			return result;
		}
	};

	class BoundaryConditions1D
	{
		std::vector<BoundaryCondition1D> _conditions;

	public:
		BoundaryConditions1D() = default;
		BoundaryConditions1D(std::initializer_list<BoundaryCondition1D> conditions)
			: _conditions(conditions) { }

		int size() const noexcept { return static_cast<int>(_conditions.size()); }
		bool empty() const noexcept { return _conditions.empty(); }
		const std::vector<BoundaryCondition1D>& conditions() const noexcept { return _conditions; }

		const BoundaryCondition1D& operator[](int index) const { return _conditions.at(static_cast<size_t>(index)); }

		void Add(const BoundaryCondition1D& condition)
		{
			_conditions.push_back(condition);
		}

		void validate(const FunctionSpace1D& space, Real tolerance = Defaults::VectorIsEqualTolerance) const
		{
			for (const auto& condition : _conditions)
				condition.validate(space, tolerance);
		}
	};
}

#endif // MML_BOUNDARY_CONDITION_1D_H
