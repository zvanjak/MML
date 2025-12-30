///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        FunctionsAnalyzer.h                                                 ///
///  Description: Function analysis tools (zeros, extrema, inflection points)         ///
///               Asymptotes, periodicity, series expansion analysis                  ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_FUNCTION_ANALYZER_H
#define MML_FUNCTION_ANALYZER_H

#include "MMLBase.h"

#include "interfaces/IFunction.h"

#include "base/VectorN.h"

#include "core/Integration.h"
#include "core/FunctionHelpers.h"
#include "core/FieldOperations.h"

#include "algorithms/RootFinding.h"
#include "algorithms/Statistics.h"

// Undefine Windows.h macros that conflict with our enum values
#ifdef _WIN32
#ifdef INFINITE
#undef INFINITE
#endif
#ifdef DOMAIN
#undef DOMAIN
#endif
#endif


namespace MML
{
	// Critical point classification
	enum class CriticalPointType
	{
		LOCAL_MINIMUM,
		LOCAL_MAXIMUM,
		SADDLE_POINT    // For 1D functions, this would be a point where f'=0 but f'' changes sign (rare)
	};

	struct CriticalPoint
	{
		Real x;
		CriticalPointType type;
		Real value;  // f(x)
	};

	// Discontinuity classification
	enum class DiscontinuityType
	{
		JUMP,           // Left and right limits exist but differ (e.g., step function)
		REMOVABLE,      // Limits exist and agree, but f(x) is different or undefined
		INFINITE,       // At least one limit is infinite (e.g., 1/x at x=0)
		OSCILLATORY,    // Limit does not exist due to oscillation (e.g., sin(1/x) at x=0)
		UNKNOWN         // Unable to classify
	};

	struct DiscontinuityPoint
	{
		Real x;
		DiscontinuityType type;
		Real leftLimit;       // lim(x→a⁻) f(x)
		Real rightLimit;      // lim(x→a⁺) f(x)
		Real valueAtPoint;    // f(a), may be NaN if undefined
		Real jumpSize;        // |rightLimit - leftLimit| for jump discontinuities
	};
	class RealFunctionAnalyzer
	{
		const IRealFunction& _f;
		std::string _funcName;
	public:
		RealFunctionAnalyzer(const IRealFunction& f) : _f(f) {}
		RealFunctionAnalyzer(const IRealFunction& f, std::string inName) : _f(f), _funcName(inName) {}

		void PrintPointAnalysis(Real x, Real eps = 1e-6)
		{
			std::cout << "Function analysis at point: " << x << ":" << std::endl;
			std::cout << "  Defined at point:      " << (isDefinedAtPoint(x) ? "yes" : "no") << std::endl;
			std::cout << "  Continuous at point:   " << (isContinuousAtPoint(x, eps) ? "yes" : "no") << std::endl;
			std::cout << "  Derivative defined:    " << (isDerivativeDefinedAtPoint(x, eps) ? "yes" : "no") << std::endl;
			std::cout << "  Is inflection point:   " << (isInflectionPoint(x, eps) ? "yes" : "no") << std::endl;
		}

		void PrintIntervalAnalysis(Real x1, Real x2, int numPoints, Real eps = 1e-6)
		{
			if (_funcName != "")
				std::cout << std::fixed << "f(x) = " << _funcName << " - Function analysis in interval [" << x1 << ", " << x2 << "] with " << numPoints << " points:" << std::endl;
			else
				std::cout << std::fixed << "Function analysis in interval [" << x1 << ", " << x2 << "] with " << numPoints << " points:" << std::endl;

			bool isDef = true;
			bool isCont = true;
			std::vector<Real> _notDefinedPoints;
			std::vector<Real> _notContinuousPoints;

			for (int i = 0; i < numPoints; i++)
			{
				Real x = x1 + i * (x2 - x1) / numPoints;
				if (!isDefinedAtPoint(x))
				{
					isDef = false;
					_notDefinedPoints.push_back(x);
				}
				else if (!isContinuousAtPoint(x, eps))
				{
					isCont = false;
					_notContinuousPoints.push_back(x);
				}
			}

			std::cout << "  Defined    : " << (isDef ? "yes" : "no");
			if (!_notDefinedPoints.empty())
			{
				std::cout << "  Not defined at points: ";
				for (int i = 0; i < _notDefinedPoints.size(); i++)
					std::cout << _notDefinedPoints[i] << " ";
				std::cout << std::endl;
			}
			else
				std::cout << std::endl;
			std::cout << "  Continuous : " << (isCont ? "yes" : "no");
			if (!_notContinuousPoints.empty())
			{
				std::cout << "  Not continuous at points: ";
				for (int i = 0; i < _notContinuousPoints.size(); i++)
					std::cout << _notContinuousPoints[i] << " ";
				std::cout << std::endl;
			}
			else
				std::cout << std::endl;
			std::cout << "  Monotonic  : " << (isMonotonic(x1, x2, numPoints) ? "yes" : "no") << std::endl;
			std::cout << "  Min        : " << MinInNPoints(x1, x2, numPoints) << std::endl;
			std::cout << "  Max        : " << MaxInNPoints(x1, x2, numPoints) << std::endl;
		}
		void PrintDetailedIntervalAnalysis(Real x1, Real x2, int numPoints, Real eps = 1e-6)
		{
			std::cout << std::fixed << "Function analysis in interval [" << x1 << ", " << x2 << "] with " << numPoints << " points:" << std::endl;
			std::cout << " Point " << "           Value     " << "         First.der.     " << "     Sec.der.     " << "Defined " << " Continuous " << " Deriv.Def. " << " Inflection " << std::endl;
			for (int i = 0; i < numPoints; i++)
			{
				Real x = x1 + i * (x2 - x1) / numPoints;

				std::cout << std::setw(6) << x << " :  ";
				std::cout << std::setw(16) << _f(x) << " :  ";
				std::cout << std::setw(16) << Derivation::NDer1(_f, x) << " :  ";
				std::cout << std::setw(16) << Derivation::NSecDer2(_f, x) << " :  ";
				std::cout << (isDefinedAtPoint(x) ? "   yes   " : "    no   ");
				std::cout << (isContinuousAtPoint(x, eps) ? "   yes   " : "    no   ");
				std::cout << (isDerivativeDefinedAtPoint(x, eps) ? "   yes   " : "    no   ");
				std::cout << (isInflectionPoint(x, eps) ? "     yes   " : "      no   ") << std::endl;
			}
		}

		std::vector<Real> GetRoots(Real x1, Real x2, Real eps)
		{
			std::vector<Real> roots;
			Real step = (x2 - x1) / 1000;
			
			// Check if x1 is a root (within tolerance)
			Real val_x1 = _f(x1);
			if (std::abs(val_x1) < eps)
				roots.push_back(x1);
			
			Real prev = val_x1;
			for (int i = 1; i <= 1000; i++)
			{
				Real x_curr = x1 + i * step;
				Real curr = _f(x_curr);
				
				// Check if current point is a root
				if (std::abs(curr) < eps && (roots.empty() || std::abs(x_curr - roots.back()) > eps))
				{
					roots.push_back(x_curr);
				}
				// Check for sign change (bracket around root)
				else if (prev * curr < 0)
				{
					Real root = RootFinding::FindRootBisection(_f, x1 + (i - 1) * step, x_curr, eps);
					// Avoid duplicate roots
					if (roots.empty() || std::abs(root - roots.back()) > eps)
						roots.push_back(root);
				}
				prev = curr;
			}
			return roots;
		}
		Vector<Real> GetLocalOptimums(Real x1, Real x2, Real eps = 1e-6)
		{
			Vector<Real> optimums;

			// Create a wrapper for f'(x) to find critical points where f'(x) = 0
			std::function<Real(Real)> derivative = [this](Real x) -> Real {
				return Derivation::NDer4(_f, x);
			};
			RealFunctionFromStdFunc derivativeFunc(derivative);
			RealFunctionAnalyzer derivAnalyzer(derivativeFunc);

			// Find roots of f'(x) - these are potential critical points
			std::vector<Real> criticalPoints = derivAnalyzer.GetRoots(x1, x2, eps);

			// For each critical point, verify it's a local optimum using second derivative test
			for (Real x : criticalPoints)
			{
				if (isLocalOptimum(x, eps))
				{
					optimums.push_back(x);
				}
			}

			return optimums;
		}
		std::vector<CriticalPoint> GetLocalOptimumsClassified(Real x1, Real x2, Real eps = 1e-6)
		{
			std::vector<CriticalPoint> classified;

			// Create a wrapper for f'(x) to find critical points where f'(x) = 0
			std::function<Real(Real)> derivative = [this](Real x) -> Real {
				return Derivation::NDer4(_f, x);
			};
			RealFunctionFromStdFunc derivativeFunc(derivative);
			RealFunctionAnalyzer derivAnalyzer(derivativeFunc);

			// Find roots of f'(x) - these are potential critical points
			std::vector<Real> criticalPoints = derivAnalyzer.GetRoots(x1, x2, eps);

			// For each critical point, verify it's a local optimum and classify
			for (Real x : criticalPoints)
			{
				if (isLocalOptimum(x, eps))
				{
					CriticalPoint cp;
					cp.x = x;
					cp.value = _f(x);
					
					// Classify using second derivative test
					Real secDeriv = Derivation::NSecDer4(_f, x);
					
					if (secDeriv > eps)
						cp.type = CriticalPointType::LOCAL_MINIMUM;
					else if (secDeriv < -eps)
						cp.type = CriticalPointType::LOCAL_MAXIMUM;
					else
						cp.type = CriticalPointType::SADDLE_POINT;  // Very rare in 1D
					
					classified.push_back(cp);
				}
			}

			return classified;
		}
		Vector<Real> GetInflectionPoints(Real x1, Real x2, Real eps = 1e-6)
		{
			Vector<Real> inflection_points;

			// Create a wrapper for f''(x) to find inflection points where f''(x) = 0
			std::function<Real(Real)> secondDerivative = [this](Real x) -> Real {
				return Derivation::NSecDer4(_f, x);
			};
			RealFunctionFromStdFunc secondDerivFunc(secondDerivative);
			RealFunctionAnalyzer secDerivAnalyzer(secondDerivFunc);

			// Find roots of f''(x) - these are potential inflection points
			std::vector<Real> candidates = secDerivAnalyzer.GetRoots(x1, x2, eps);

			// For each candidate, verify it's an inflection point (concavity changes)
			for (Real x : candidates)
			{
				if (isInflectionPoint(x, eps))
				{
					inflection_points.push_back(x);
				}
			}

			return inflection_points;
		}
		
		bool isDefinedAtPoint(Real x) const
		{
			Real y = _f(x);
			return !std::isnan(y) && !std::isinf(y);
		}
		bool isDerivativeDefinedAtPoint(Real x, Real eps = 1e-6)
		{
			// Derivative exists at x if left and right derivatives exist and are equal
			// Left derivative:  lim (h→0-) [f(x+h) - f(x)] / h
			// Right derivative: lim (h→0+) [f(x+h) - f(x)] / h
			
			if (!isDefinedAtPoint(x))
				return false;

			// Check if function is defined in a neighborhood
			Real h = eps * 10;
			if (!isDefinedAtPoint(x - h) || !isDefinedAtPoint(x + h))
				return false;

			// Compute left and right derivatives using NDer2
			Real leftDeriv = Derivation::NDer2Left(_f, x, eps);
			Real rightDeriv = Derivation::NDer2Right(_f, x, eps);

			// Check if both derivatives are defined (not NaN or Inf)
			if (std::isnan(leftDeriv) || std::isinf(leftDeriv) ||
				std::isnan(rightDeriv) || std::isinf(rightDeriv))
				return false;

			// Check if left and right derivatives match (within tolerance)
			Real deriv_diff = std::abs(leftDeriv - rightDeriv);
			Real avg_deriv = (std::abs(leftDeriv) + std::abs(rightDeriv)) / 2.0;
			
			// Use relative tolerance for large derivatives, absolute for small
			Real tolerance = eps * 100 + eps * 100 * avg_deriv;
			
			return deriv_diff < tolerance;
		}
		bool isContinuousAtPoint(Real x, Real eps)
		{
			// Function is continuous at x if:
			// lim(h→0) f(x+h) = f(x) and lim(h→0) f(x-h) = f(x)
			// 
			// We use epsilon-delta definition: for any eps > 0, there exists delta > 0
			// such that |x - y| < delta => |f(x) - f(y)| < eps
			
			if (!isDefinedAtPoint(x))
				return false;

			Real val = _f(x);
			
			// Adaptive tolerance: use absolute tolerance for small values, 
			// relative tolerance for large values
			Real abs_tol = eps * 100;  // Base absolute tolerance
			Real rel_tol = eps * 100;  // Base relative tolerance
			
			// Test with progressively smaller h values
			// If function is continuous, |f(x±h) - f(x)| should decrease as h decreases
			Real h = eps * 10;
			Real prev_max_diff = std::numeric_limits<Real>::max();
			
			for (int i = 0; i < 5; i++)
			{
				Real left = _f(x - h);
				Real right = _f(x + h);
				
				// Check if left and right values are defined
				if (!isDefinedAtPoint(x - h) || !isDefinedAtPoint(x + h))
					return false;
				
				// Calculate differences
				Real diff_left = std::abs(left - val);
				Real diff_right = std::abs(right - val);
				Real max_diff = std::max(diff_left, diff_right);
				
				// For continuity, the difference should be getting smaller
				// Use combined absolute and relative tolerance
				Real tolerance = abs_tol + rel_tol * std::abs(val);
				
				if (max_diff < tolerance)
					return true;  // Close enough
				
				// If difference is not decreasing, it's likely discontinuous
				if (i > 0 && max_diff >= prev_max_diff * 0.9)
					return false;
				
				prev_max_diff = max_diff;
				h /= 2.0;  // Try with smaller h
			}
			
			return false;
		}
		bool isLocalOptimum(Real x, Real eps)
		{
			// For a local optimum:
			// 1. First derivative must be (approximately) zero: f'(x) ≈ 0
			// 2. Second derivative must be non-zero with consistent sign: f''(x) ≠ 0
			//    - f''(x) > 0 => local minimum
			//    - f''(x) < 0 => local maximum
			
			// Check first derivative is near zero
			Real first_der = Derivation::NDer4(_f, x);
			if (std::abs(first_der) > eps)
				return false;  // Not a stationary point

			// Check second derivative has consistent sign (non-zero)
			Real left_sec_der = Derivation::NSecDer4(_f, x - 4 * eps);
			Real right_sec_der = Derivation::NSecDer4(_f, x + 4 * eps);

			// Both should have the same sign AND be non-zero
			return left_sec_der * right_sec_der > 0;
		}
		bool isInflectionPoint(Real x, Real eps)
		{
			// Inflection point: second derivative changes sign (concavity changes)
			// Note: First derivative does NOT need to be zero (that's for stationary points)
			Real left_sec_der = Derivation::NSecDer4(_f, x - 4 * eps);
			Real right_sec_der = Derivation::NSecDer4(_f, x + 4 * eps);

			// Check if second derivatives have opposite signs (sign change)
			return left_sec_der * right_sec_der < 0;
		}
		bool isContinuous(Real x1, Real x2, int numPoints, Real eps = 1e-6, std::vector<Real>* discontinuities = nullptr)
		{
			// Multi-stage approach to detect discontinuities:
			// 1) Sample at numPoints
			// 2) Check continuity at each point using isContinuousAtPoint()
			// 3) Adaptive refinement where suspicious behavior detected
			
			bool allContinuous = true;
			Real step = (x2 - x1) / numPoints;
			std::vector<Real> suspiciousPoints;
			
			// Stage 1: Initial sampling
			for (int i = 0; i <= numPoints; i++)
			{
				Real x = x1 + i * step;
				if (!isContinuousAtPoint(x, eps))
				{
					allContinuous = false;
					suspiciousPoints.push_back(x);
				}
			}
			
			// Stage 2: Check for sign changes in derivative (potential discontinuities)
			// Sample derivative at midpoints to catch jump discontinuities
			for (int i = 0; i < numPoints; i++)
			{
				Real x_left = x1 + i * step;
				Real x_right = x1 + (i + 1) * step;
				Real x_mid = (x_left + x_right) / 2.0;
				
				// Check for large derivative changes (potential discontinuity)
				Real der_left = Derivation::NDer2(_f, x_left);
				Real der_right = Derivation::NDer2(_f, x_right);
				
				// If derivative changes dramatically or function value jumps
				Real val_left = _f(x_left);
				Real val_right = _f(x_right);
				Real val_mid = _f(x_mid);
				
				// Check for non-smooth behavior (jump or cusp)
				bool potentialJump = std::abs(val_right - val_left) > 10 * eps * std::abs(val_left + val_right);
				bool derivativeJump = !std::isnan(der_left) && !std::isnan(der_right) && 
					                   std::abs(der_right - der_left) > 100 * eps;
				
				if (potentialJump || derivativeJump)
				{
					// Refine search in this interval
					if (!isContinuousAtPoint(x_mid, eps))
					{
						allContinuous = false;
						// Avoid duplicates
						if (std::find(suspiciousPoints.begin(), suspiciousPoints.end(), x_mid) == suspiciousPoints.end())
						{
							suspiciousPoints.push_back(x_mid);
						}
					}
				}
			}
			
			// Return discontinuity locations if requested
			if (discontinuities != nullptr && !suspiciousPoints.empty())
			{
				*discontinuities = suspiciousPoints;
			}
			
			return allContinuous;
		}
		bool isMonotonic(Real x1, Real x2, int numPoints)
		{
			Real step = (x2 - x1) / numPoints;
			Real prev = _f(x1 + step);
			if (_f(x1) < _f(x1 + step))
			{
				for (int i = 2; i < numPoints; i++)
				{
					Real curr = _f(x1 + i * step);
					if (curr < prev)
						return false;
					prev = curr;
				}
			}
			else
			{
				for (int i = 2; i < numPoints; i++)
				{
					Real curr = _f(x1 + i * step);
					if (curr > prev)
						return false;
					prev = curr;
				}
			}
			return true;
		}
		Real MinInNPoints(Real x1, Real x2, int numPoints)
		{
			Real step = (x2 - x1) / numPoints;
			Real min = _f(x1);

			for (int i = 1; i <= numPoints; i++)
			{
				Real curr = _f(x1 + i * step);
				if (curr < min)
					min = curr;
			}
			return min;
		}
		Real MaxInNPoints(Real x1, Real x2, int numPoints)
		{
			Real step = (x2 - x1) / numPoints;
			Real max = _f(x1);

			for (int i = 1; i < numPoints; i++)
			{
				Real curr = _f(x1 + i * step);
				if (curr > max)
					max = curr;
			}
			return max;
		}

		//////////////////////////////////////////////////////////////////////////
		// CONTINUITY ANALYSIS METHODS
		//////////////////////////////////////////////////////////////////////////

		// Compute left limit: lim(x→a⁻) f(x)
		Real ComputeLeftLimit(Real x, Real eps = 1e-6, int maxSteps = 20) const
		{
			// Approach from the left with decreasing step sizes
			Real h = eps * 100;
			Real prevValue = std::numeric_limits<Real>::quiet_NaN();
			int stableCount = 0;

			for (int i = 0; i < maxSteps; i++)
			{
				Real testPoint = x - h;
				if (!isDefinedAtPoint(testPoint))
					return std::numeric_limits<Real>::quiet_NaN();

				Real value = _f(testPoint);
				
				// Check for convergence
				if (!std::isnan(prevValue))
				{
					Real diff = std::abs(value - prevValue);
					Real tolerance = eps * (1.0 + std::abs(value));
					
					if (diff < tolerance)
					{
						stableCount++;
						if (stableCount >= 3)  // Require 3 stable iterations
							return value;
					}
					else
					{
						stableCount = 0;
					}
				}
				
				prevValue = value;
				h /= 2.0;
			}
			
			return prevValue;  // Best approximation
		}

		// Compute right limit: lim(x→a⁺) f(x)
		Real ComputeRightLimit(Real x, Real eps = 1e-6, int maxSteps = 20) const
		{
			// Approach from the right with decreasing step sizes
			Real h = eps * 100;
			Real prevValue = std::numeric_limits<Real>::quiet_NaN();
			int stableCount = 0;

			for (int i = 0; i < maxSteps; i++)
			{
				Real testPoint = x + h;
				if (!isDefinedAtPoint(testPoint))
					return std::numeric_limits<Real>::quiet_NaN();

				Real value = _f(testPoint);
				
				// Check for convergence
				if (!std::isnan(prevValue))
				{
					Real diff = std::abs(value - prevValue);
					Real tolerance = eps * (1.0 + std::abs(value));
					
					if (diff < tolerance)
					{
						stableCount++;
						if (stableCount >= 3)  // Require 3 stable iterations
							return value;
					}
					else
					{
						stableCount = 0;
					}
				}
				
				prevValue = value;
				h /= 2.0;
			}
			
			return prevValue;  // Best approximation
		}

		// Classify a discontinuity at point x
		DiscontinuityPoint ClassifyDiscontinuity(Real x, Real eps = 1e-6) const
		{
			DiscontinuityPoint disc;
			disc.x = x;
			disc.valueAtPoint = isDefinedAtPoint(x) ? _f(x) : std::numeric_limits<Real>::quiet_NaN();
			disc.leftLimit = ComputeLeftLimit(x, eps);
			disc.rightLimit = ComputeRightLimit(x, eps);
			disc.jumpSize = 0.0;

			// Check if limits exist (are finite)
			bool leftLimitExists = !std::isnan(disc.leftLimit) && !std::isinf(disc.leftLimit);
			bool rightLimitExists = !std::isnan(disc.rightLimit) && !std::isinf(disc.rightLimit);
			
			// Check if limits are very large (approaching infinity)
			const Real largeThreshold = 1e10;
			bool leftLimitVeryLarge = std::abs(disc.leftLimit) > largeThreshold;
			bool rightLimitVeryLarge = std::abs(disc.rightLimit) > largeThreshold;

			// Classify discontinuity type
			if (!leftLimitExists && !rightLimitExists)
			{
				disc.type = DiscontinuityType::OSCILLATORY;  // Or could be UNKNOWN
			}
			else if (!leftLimitExists || !rightLimitExists || leftLimitVeryLarge || rightLimitVeryLarge)
			{
				// At least one limit is infinite or doesn't exist
				disc.type = DiscontinuityType::INFINITE;
			}
			else if (std::abs(disc.leftLimit - disc.rightLimit) > eps * 1000 * (1.0 + std::abs(disc.leftLimit)))
			{
				// Limits exist but differ significantly - JUMP discontinuity
				disc.type = DiscontinuityType::JUMP;
				disc.jumpSize = std::abs(disc.rightLimit - disc.leftLimit);
			}
			else
			{
				// Limits exist and agree
				Real limitValue = (disc.leftLimit + disc.rightLimit) / 2.0;
				
				if (!isDefinedAtPoint(x))
				{
					disc.type = DiscontinuityType::REMOVABLE;  // Undefined at point
				}
				else if (std::abs(_f(x) - limitValue) > eps * 1000 * (1.0 + std::abs(limitValue)))
				{
					disc.type = DiscontinuityType::REMOVABLE;  // Value differs from limit
				}
				else
				{
					// This shouldn't happen - if we're here, it should be continuous
					disc.type = DiscontinuityType::UNKNOWN;
				}
			}

			return disc;
		}

		// Find all discontinuities in an interval
		std::vector<DiscontinuityPoint> FindDiscontinuities(Real x1, Real x2, int numPoints = 1000, Real eps = 1e-6)
		{
			std::vector<DiscontinuityPoint> discontinuities;
			Real step = (x2 - x1) / numPoints;

			// Stage 1: Check each sample point for continuity
			std::vector<Real> candidates;
			for (int i = 0; i <= numPoints; i++)
			{
				Real x = x1 + i * step;
				if (!isContinuousAtPoint(x, eps))
				{
					candidates.push_back(x);
				}
			}

			// Stage 2: Check between sample points for rapid changes (potential jumps)
			for (int i = 0; i < numPoints; i++)
			{
				Real x_left = x1 + i * step;
				Real x_right = x1 + (i + 1) * step;
				Real x_mid = (x_left + x_right) / 2.0;

				if (!isDefinedAtPoint(x_left) || !isDefinedAtPoint(x_right))
				{
					// If points aren't defined, check the midpoint
					if (!isContinuousAtPoint(x_mid, eps))
					{
						bool alreadyFound = false;
						for (const auto& candX : candidates)
						{
							if (std::abs(candX - x_mid) < step / 2.0)
							{
								alreadyFound = true;
								break;
							}
						}
						if (!alreadyFound)
							candidates.push_back(x_mid);
					}
					continue;
				}

				Real val_left = _f(x_left);
				Real val_right = _f(x_right);

				// Check for large jump between consecutive points
				// Use a more sensitive threshold for detecting jumps
				Real jumpThreshold = eps * 100;  // More sensitive
				
				// Also check relative changes for functions with large values
				Real maxVal = std::max(std::abs(val_left), std::abs(val_right));
				if (maxVal > 1.0)
					jumpThreshold = std::max(jumpThreshold, maxVal * eps * 10);

				if (std::abs(val_right - val_left) > jumpThreshold)
				{
					// Potential discontinuity - check multiple points in this interval
					for (int j = 1; j <= 4; j++)
					{
						Real x_test = x_left + (j / 5.0) * (x_right - x_left);
						if (!isContinuousAtPoint(x_test, eps))
						{
							// Avoid duplicates
							bool alreadyFound = false;
							for (const auto& candX : candidates)
							{
								if (std::abs(candX - x_test) < step / 5.0)
								{
									alreadyFound = true;
									break;
								}
							}
							if (!alreadyFound)
								candidates.push_back(x_test);
						}
					}
				}
			}

			// Stage 3: Classify each candidate discontinuity
			for (Real x : candidates)
			{
				DiscontinuityPoint disc = ClassifyDiscontinuity(x, eps);
				discontinuities.push_back(disc);
			}

			// Sort by x coordinate
			std::sort(discontinuities.begin(), discontinuities.end(),
				[](const DiscontinuityPoint& a, const DiscontinuityPoint& b) {
					return a.x < b.x;
				});

			return discontinuities;
		}

		// Print detailed continuity analysis for an interval
		void PrintContinuityAnalysis(Real x1, Real x2, int numPoints = 1000, Real eps = 1e-6)
		{
			std::cout << std::fixed << std::setprecision(6);
			std::cout << "\n=== CONTINUITY ANALYSIS ===" << std::endl;
			std::cout << "Interval: [" << x1 << ", " << x2 << "]" << std::endl;
			std::cout << "Sample points: " << numPoints << std::endl;
			std::cout << "Tolerance: " << eps << "\n" << std::endl;

			auto discontinuities = FindDiscontinuities(x1, x2, numPoints, eps);

			if (discontinuities.empty())
			{
				std::cout << "✓ Function is CONTINUOUS throughout the interval." << std::endl;
			}
			else
			{
				std::cout << "✗ Function has " << discontinuities.size() << " discontinuity point(s):\n" << std::endl;

				for (size_t i = 0; i < discontinuities.size(); i++)
				{
					const auto& disc = discontinuities[i];
					std::cout << "Discontinuity #" << (i + 1) << " at x = " << disc.x << std::endl;
					
					// Print type
					std::cout << "  Type: ";
					switch (disc.type)
					{
					case DiscontinuityType::JUMP:
						std::cout << "JUMP (jump size: " << disc.jumpSize << ")" << std::endl;
						break;
					case DiscontinuityType::REMOVABLE:
						std::cout << "REMOVABLE" << std::endl;
						break;
					case DiscontinuityType::INFINITE:
						std::cout << "INFINITE" << std::endl;
						break;
					case DiscontinuityType::OSCILLATORY:
						std::cout << "OSCILLATORY (limit does not exist)" << std::endl;
						break;
					case DiscontinuityType::UNKNOWN:
						std::cout << "UNKNOWN" << std::endl;
						break;
					}

					// Print limits
					std::cout << "  Left limit:  ";
					if (std::isnan(disc.leftLimit))
						std::cout << "undefined (NaN)" << std::endl;
					else if (std::isinf(disc.leftLimit))
						std::cout << (disc.leftLimit > 0 ? "+∞" : "-∞") << std::endl;
					else
						std::cout << disc.leftLimit << std::endl;

					std::cout << "  Right limit: ";
					if (std::isnan(disc.rightLimit))
						std::cout << "undefined (NaN)" << std::endl;
					else if (std::isinf(disc.rightLimit))
						std::cout << (disc.rightLimit > 0 ? "+∞" : "-∞") << std::endl;
					else
						std::cout << disc.rightLimit << std::endl;

					std::cout << "  f(" << disc.x << "):      ";
					if (std::isnan(disc.valueAtPoint))
						std::cout << "undefined" << std::endl;
					else if (std::isinf(disc.valueAtPoint))
						std::cout << (disc.valueAtPoint > 0 ? "+∞" : "-∞") << std::endl;
					else
						std::cout << disc.valueAtPoint << std::endl;

					std::cout << std::endl;
				}
			}
		}
	
		// If given a periodic function, calculates the average period of the roots
		// of the function in the given interval
		Real calcRootsPeriod(Real t1, Real t2, int numPoints)
		{
			Vector<Real> root_brack_x1(10), root_brack_x2(10);
			int	numFoundRoots = RootFinding::FindRootBrackets(_f, t1, t2, numPoints, root_brack_x1, root_brack_x2);

			if( numFoundRoots == 0 )
				return 0.0;

			Vector<Real> roots(numFoundRoots);
			Vector<Real> rootDiffs(numFoundRoots - 1);
			for (int i = 0; i < numFoundRoots; i++)
			{
				roots[i] = RootFinding::FindRootBisection(_f, root_brack_x1[i], root_brack_x2[i], 1e-7);

				if (i > 0)
					rootDiffs[i - 1] = roots[i] - roots[i - 1];
			}

			if( numFoundRoots == 1 )
				return 0.0; // only one root found, no period

			return Statistics::Avg(rootDiffs);
		}
	};

	class RealFunctionComparer
	{
		IRealFunction& _f1;
		IRealFunction& _f2;

	public:
		RealFunctionComparer(IRealFunction& f1, IRealFunction& f2) : _f1(f1), _f2(f2) {}

		Real getAbsDiffSum(Real a, Real b, int numPoints)
		{
			Real step = (b - a) / numPoints;
			Real sum = 0.0;

			for (int i = 0; i < numPoints; i++)
				sum += std::abs(_f1(a + i * step) - _f2(a + i * step));

			return sum;
		}
		Real getAbsDiffAvg(Real a, Real b, int numPoints)
		{
			return getAbsDiffSum(a, b, numPoints) / numPoints;
		}
		Real getAbsDiffMax(Real a, Real b, int numPoints)
		{
			Real step = (b - a) / numPoints;
			Real max = std::abs(_f1(a) - _f2(a));

			for (int i = 0; i < numPoints; i++)
			{
				Real diff = std::abs(_f1(a + i * step) - _f2(a + i * step));
				if (diff > max)
					max = diff;
			}
			return max;
		}
		Real getRelDiffSum(Real a, Real b, int numPoints)
		{
			Real step = (b - a) / numPoints;
			Real sum = 0.0;

			for (int i = 0; i < numPoints; i++)
				if (_f1(a + i * step) != 0.0)
					sum += std::abs(_f1(a + i * step) - _f2(a + i * step)) / std::abs(_f1(a + i * step));

			return sum;
		}
		Real getRelDiffAvg(Real a, Real b, int numPoints)
		{
			return getRelDiffSum(a, b, numPoints) / numPoints;
		}
		Real getRelDiffMax(Real a, Real b, int numPoints)
		{
			Real step = (b - a) / numPoints;
			Real max = 0.0;

			for (int i = 0; i < numPoints; i++)
			{
				if (_f1(a + i * step) != 0.0)
				{
					Real diff = std::abs(_f1(a + i * step) - _f2(a + i * step)) / std::abs(_f1(a + i * step));
					if (diff > max)
						max = diff;
				}
			}
			return max;
		}

		///////////                  Integration measures                /////////
		Real getIntegratedDiff(Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP)
		{
			return getIntegratedDiff(_f1, _f2, a, b, method);
		}
		static Real getIntegratedDiff(IRealFunction& f1, IRealFunction& f2, Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP)
		{
			RealFuncDiffHelper helper(f1, f2);

			switch (method)
			{
			case IntegrationMethod::SIMPSON:
				return IntegrateSimpson(helper, a, b);
			case IntegrationMethod::ROMBERG:
				return IntegrateRomberg(helper, a, b);
			default:
				return IntegrateTrap(helper, a, b);
			}
		}

		Real getIntegratedAbsDiff(Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP)
		{
			return getIntegratedAbsDiff(_f1, _f2, a, b, method);
		}
		static Real getIntegratedAbsDiff(IRealFunction& f1, IRealFunction& f2, Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP)
		{
			RealFuncAbsDiffHelper helper(f1, f2);

			switch (method)
			{
			case IntegrationMethod::SIMPSON:
				return IntegrateSimpson(helper, a, b);
			case IntegrationMethod::ROMBERG:
				return IntegrateRomberg(helper, a, b);
			default:
				return IntegrateTrap(helper, a, b);
			}
		}

		Real getIntegratedSqrDiff(Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP)
		{
			return getIntegratedSqrDiff(_f1, _f2, a, b, method);
		}
		static Real getIntegratedSqrDiff(IRealFunction& f1, IRealFunction& f2, Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP)
		{
			RealFuncDiffSqrHelper helper(f1, f2);

			switch (method)
			{
			case IntegrationMethod::SIMPSON:
				return IntegrateSimpson(helper, a, b);
			case IntegrationMethod::ROMBERG:
				return IntegrateRomberg(helper, a, b);
			default:
				return IntegrateTrap(helper, a, b);
			}
		}
	};
}
#endif