///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        InterpolatedRealFunctionAdvanced.h                                  ///
///  Description: Hermite and Akima interpolation                                    ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

/// @file InterpolatedRealFunctionAdvanced.h
/// @brief Hermite and Akima interpolation.
/// @ingroup Interpolation

#if !defined MML_INTERPOLATED_REAL_FUNCTION_ADVANCED_H
#define MML_INTERPOLATED_REAL_FUNCTION_ADVANCED_H

#include <algorithm>
#include <cmath>

#include <mml/base/InterpolatedFunctions/InterpolatedRealFunctionLinear.h>

namespace MML {
	/////////////////////////////////////////////////////////////////////////////////////
	///                  USER-DERIVATIVE HERMITE INTERPOLATION                       ///
	/////////////////////////////////////////////////////////////////////////////////////

	class HermiteInterpRealFunc : public RealFunctionInterpolated {
		Vector<Real> _derivatives;

	public:
		HermiteInterpRealFunc(const Vector<Real>& xv, const Vector<Real>& yv,
			const Vector<Real>& derivatives)
			: RealFunctionInterpolated(xv, yv, 2), _derivatives(derivatives) {
			if (derivatives.size() != xv.size())
				throw RealFuncInterpInitError("HermiteInterpRealFunc: derivative vector size mismatch");
			validateOrderedNodes();
		}

		Real calcInterpValue(int startInd, Real x) const override {
			return detail::EvaluateCubicHermiteSegment(x, X(startInd), X(startInd + 1),
				Y(startInd), Y(startInd + 1), _derivatives[startInd], _derivatives[startInd + 1]);
		}

		const char* InterpolationMethodName() const override { return "CubicHermite"; }

		Real Derivative(Real x) const {
			int interval = locate(x);
			return detail::DifferentiateCubicHermiteSegment(x, X(interval), X(interval + 1),
				Y(interval), Y(interval + 1), _derivatives[interval], _derivatives[interval + 1]);
		}

		Real Integrate(Real from, Real to) const {
			return integratePiecewise(from, to);
		}

		Real GetDerivative(int index) const {
			if (index < 0 || index >= _derivatives.size())
				throw IndexError("HermiteInterpRealFunc derivative index out of range");
			return _derivatives[index];
		}

	private:
		Real integratePiecewise(Real from, Real to) const {
			if (from > to) return -integratePiecewise(to, from);
			if (from < MinX() || to > MaxX())
				throw DomainError("Hermite interpolation: integration bounds exceed data domain");
			Real result = 0.0;
			for (int i = 0; i < getNumPoints() - 1; ++i) {
				Real segmentMin = std::min(X(i), X(i + 1));
				Real segmentMax = std::max(X(i), X(i + 1));
				Real left = std::max(from, segmentMin);
				Real right = std::min(to, segmentMax);
				if (left < right)
					result += detail::IntegrateCubicHermiteSegment(left, right, X(i), X(i + 1),
						Y(i), Y(i + 1), _derivatives[i], _derivatives[i + 1]);
			}
			return result;
		}

		void validateOrderedNodes() const {
			const bool ascending = X(getNumPoints() - 1) > X(0);
			for (int i = 1; i < getNumPoints(); ++i)
				if ((X(i) > X(i - 1)) != ascending)
					throw RealFuncInterpInitError("HermiteInterpRealFunc: x-values must be strictly monotonic");
		}
	};


	/////////////////////////////////////////////////////////////////////////////////////
	///                           AKIMA INTERPOLATION                                  ///
	/////////////////////////////////////////////////////////////////////////////////////

	class AkimaInterpRealFunc : public RealFunctionInterpolated {
		Vector<Real> _derivatives;

	public:
		AkimaInterpRealFunc(const Vector<Real>& xv, const Vector<Real>& yv)
			: RealFunctionInterpolated(xv, yv, 2), _derivatives(xv.size()) {
			if (xv.size() < 5)
				throw RealFuncInterpInitError("AkimaInterpRealFunc requires at least five points");
			validateOrderedNodes();
			initializeDerivatives();
		}

		Real calcInterpValue(int startInd, Real x) const override {
			return detail::EvaluateCubicHermiteSegment(x, X(startInd), X(startInd + 1),
				Y(startInd), Y(startInd + 1), _derivatives[startInd], _derivatives[startInd + 1]);
		}

		const char* InterpolationMethodName() const override { return "Akima"; }

		Real GetDerivative(int index) const {
			if (index < 0 || index >= _derivatives.size())
				throw IndexError("AkimaInterpRealFunc derivative index out of range");
			return _derivatives[index];
		}

	private:
		void validateOrderedNodes() const {
			const bool ascending = X(getNumPoints() - 1) > X(0);
			for (int i = 1; i < getNumPoints(); ++i)
				if ((X(i) > X(i - 1)) != ascending)
					throw RealFuncInterpInitError("AkimaInterpRealFunc: x-values must be strictly monotonic");
		}

		void initializeDerivatives() {
			const int n = getNumPoints();
			Vector<Real> slopes(n + 3);
			for (int i = 0; i < n - 1; ++i)
				slopes[i + 2] = (Y(i + 1) - Y(i)) / (X(i + 1) - X(i));

			slopes[1] = REAL(2.0) * slopes[2] - slopes[3];
			slopes[0] = REAL(2.0) * slopes[1] - slopes[2];
			slopes[n + 1] = REAL(2.0) * slopes[n] - slopes[n - 1];
			slopes[n + 2] = REAL(2.0) * slopes[n + 1] - slopes[n];

			for (int i = 0; i < n; ++i) {
				Real weightLeft = std::abs(slopes[i + 3] - slopes[i + 2]);
				Real weightRight = std::abs(slopes[i + 1] - slopes[i]);
				Real weightSum = weightLeft + weightRight;
				_derivatives[i] = weightSum > PrecisionValues<Real>::DivisionSafetyThreshold
					? (weightLeft * slopes[i + 1] + weightRight * slopes[i + 2]) / weightSum
					: REAL(0.5) * (slopes[i + 1] + slopes[i + 2]);
			}
		}
	};

} // namespace MML

#endif // MML_INTERPOLATED_REAL_FUNCTION_ADVANCED_H
