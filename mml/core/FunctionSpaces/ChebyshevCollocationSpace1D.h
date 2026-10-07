///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionSpaces/ChebyshevCollocationSpace1D.h                        ///
///  Description: Chebyshev-Lobatto collocation trial spaces in one dimension         ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_CHEBYSHEV_COLLOCATION_SPACE_1D_H
#define MML_CHEBYSHEV_COLLOCATION_SPACE_1D_H

#include <mml/core/FunctionSpaces/TrialSpace1D.h>

#include <mml/base/ChebyshevPolynom.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Vector/Vector.h>

#include <cmath>

namespace MML::FunctionSpaces
{
	class ChebyshevCollocationSpace1D : public TrialSpace1D
	{
		L2IntervalSpace _space;
		int _dimension;
		Vector<Real> _nodes;

		static Vector<Real> BuildLobattoNodes(Real min, Real max, int dimension)
		{
			if (dimension < 2)
				throw FunctionSpaceInputError("ChebyshevCollocationSpace1D: dimension must be at least 2");

			L2IntervalSpace space(min, max);
			Vector<Real> nodes(dimension);
			for (int index = 0; index < dimension; ++index) {
				Real theta = Constants::PI * static_cast<Real>(dimension - 1 - index) / static_cast<Real>(dimension - 1);
				Real xi = std::cos(theta);
				nodes[index] = (space.domainMin() + space.domainMax()) / REAL(2.0) + (space.domainLength() / REAL(2.0)) * xi;
			}
			return nodes;
		}

	public:
		ChebyshevCollocationSpace1D(Real min, Real max, int dimension)
			: _space(min, max), _dimension(dimension), _nodes(BuildLobattoNodes(min, max, dimension))
		{
		}

		const FunctionSpace1D& functionSpace() const noexcept override { return _space; }
		int dimension() const noexcept override { return _dimension; }
		bool hasNodes() const noexcept override { return true; }
		int nodeCount() const noexcept override { return _dimension; }

		Real node(int nodeIndex) const override
		{
			validateNodeIndex(nodeIndex, "ChebyshevCollocationSpace1D::node");
			return _nodes[nodeIndex];
		}

		const Vector<Real>& nodes() const noexcept { return _nodes; }

		Real basisValue(int basisIndex, Real x) const override
		{
			validateBasisIndex(basisIndex, "ChebyshevCollocationSpace1D::basisValue");
			Real xi = (REAL(2.0) * x - functionSpace().domainMin() - functionSpace().domainMax()) / functionSpace().domainLength();
			return ChebyshevT(basisIndex, xi);
		}

		Matrix<Real> firstDerivativeMatrix() const
		{
			Matrix<Real> result(_dimension, _dimension);
			const int n = _dimension - 1;
			const Real scale = REAL(2.0) / functionSpace().domainLength();

			for (int row = 0; row < _dimension; ++row) {
				int descendingRow = n - row;
				Real xRow = std::cos(Constants::PI * static_cast<Real>(descendingRow) / static_cast<Real>(n));
				Real cRow = (descendingRow == 0 || descendingRow == n) ? REAL(2.0) : REAL(1.0);
				if (descendingRow % 2 != 0)
					cRow = -cRow;

				for (int column = 0; column < _dimension; ++column) {
					int descendingColumn = n - column;
					Real xColumn = std::cos(Constants::PI * static_cast<Real>(descendingColumn) / static_cast<Real>(n));
					Real cColumn = (descendingColumn == 0 || descendingColumn == n) ? REAL(2.0) : REAL(1.0);
					if (descendingColumn % 2 != 0)
						cColumn = -cColumn;

					if (row != column)
						result(row, column) = scale * (cRow / cColumn) / (xRow - xColumn);
					else if (descendingRow == 0)
						result(row, column) = scale * (REAL(2.0) * n * n + REAL(1.0)) / REAL(6.0);
					else if (descendingRow == n)
						result(row, column) = -scale * (REAL(2.0) * n * n + REAL(1.0)) / REAL(6.0);
					else
						result(row, column) = scale * (-xRow / (REAL(2.0) * (REAL(1.0) - xRow * xRow)));
				}
			}
			return result;
		}

		Matrix<Real> secondDerivativeMatrix() const
		{
			Matrix<Real> first = firstDerivativeMatrix();
			Matrix<Real> second(_dimension, _dimension);
			for (int row = 0; row < _dimension; ++row) {
				for (int column = 0; column < _dimension; ++column) {
					Real value = REAL(0.0);
					for (int inner = 0; inner < _dimension; ++inner)
						value += first(row, inner) * first(inner, column);
					second(row, column) = value;
				}
			}
			return second;
		}
	};
}

#endif // MML_CHEBYSHEV_COLLOCATION_SPACE_1D_H
