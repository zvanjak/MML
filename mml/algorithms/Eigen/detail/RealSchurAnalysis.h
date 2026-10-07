#ifndef MML_EIGEN_DETAIL_REAL_SCHUR_ANALYSIS_H
#define MML_EIGEN_DETAIL_REAL_SCHUR_ANALYSIS_H

#include <mml/MMLBase.h>
#include <mml/base/Vector/Vector.h>
#include <mml/base/Matrix/Matrix.h>
#include <mml/base/BaseUtils.h>
#include <algorithm>
#include <cmath>
#include <complex>

namespace MML::detail {
	// =========================================================================
	// BUILDING BLOCK 4: EIGENVALUE EXTRACTION FROM 2x2 BLOCK
	// =========================================================================

	/// @struct Eigenvalue2x2Result
	/// @brief Eigenvalues of a 2x2 matrix (may be complex conjugate pair)
	struct Eigenvalue2x2Result {
		Real real1, imag1; // First eigenvalue: real1 + i*imag1
		Real real2, imag2; // Second eigenvalue: real2 + i*imag2
		bool isComplex;	   // True if complex conjugate pair
	};

	/// Compute eigenvalues of 2x2 matrix [[a, b], [c, d]]
	///
	/// Characteristic polynomial: λ² - (a+d)λ + (ad-bc) = 0
	/// λ = (a+d)/2 ± sqrt((a+d)²/4 - (ad-bc))
	/// = (a+d)/2 ± sqrt((a-d)²/4 + bc)
	inline Eigenvalue2x2Result Eigenvalues2x2(Real a, Real b, Real c, Real d) {
		Eigenvalue2x2Result result;

		Real trace = a + d;
		Real det = a * d - b * c;

		// Discriminant = trace²/4 - det = (a-d)²/4 + bc
		Real p = (a - d) / 2.0;
		Real disc = p * p + b * c;

		result.real1 = trace / 2.0;
		result.real2 = trace / 2.0;

		if (disc >= 0.0) {
			// Two real eigenvalues
			Real sqrtDisc = std::sqrt(disc);
			result.real1 += sqrtDisc;
			result.real2 -= sqrtDisc;
			result.imag1 = 0.0;
			result.imag2 = 0.0;
			result.isComplex = false;
		} else {
			// Complex conjugate pair
			Real sqrtDisc = std::sqrt(-disc);
			result.imag1 = sqrtDisc;
			result.imag2 = -sqrtDisc;
			result.isComplex = true;
		}

		return result;
	}

	/// @struct DeflationResult
	/// @brief Result of checking for deflation in Hessenberg matrix
	struct DeflationResult {
		bool canDeflate;	// True if deflation is possible
		int deflationIndex; // Index where deflation occurs (H[idx, idx-1] ≈ 0)
		int blockSize;		// Size of deflated block (1 = real eigenvalue, 2 = complex pair)
	};

	/// Check if matrix can be deflated at any position.
	///
	/// Deflation occurs when |H[k, k-1]| < tol * (|H[k-1,k-1]| + |H[k,k]|)
	/// This means the matrix decouples into independent subproblems.
	///
	/// @param H Upper Hessenberg matrix
	/// @param iLo Start index of active region
	/// @param iHi End index of active region
	/// @param tol Tolerance for considering element as zero
	/// @return DeflationResult indicating if/where deflation is possible
	inline DeflationResult CheckDeflation(const Matrix<Real>& H, int iLo, int iHi, Real tol = 1e-10) {
		DeflationResult result;
		result.canDeflate = false;
		result.deflationIndex = -1;
		result.blockSize = 0;

		// Check from bottom up for deflation
		for (int k = iHi; k > iLo; k--) {
			Real off = std::abs(H(k, k - 1));
			Real diag = std::abs(H(k - 1, k - 1)) + std::abs(H(k, k));

			if (off < tol * std::max(diag, REAL(1.0)) || off < PrecisionValues<Real>::DivisionSafetyThreshold) {
				result.canDeflate = true;
				result.deflationIndex = k;

				// Determine block size: check if this isolates a 1x1 or 2x2 block
				if (k == iHi) {
					// Single eigenvalue at position iHi
					result.blockSize = 1;
				} else {
					// Check next subdiagonal
					Real off2 = std::abs(H(k + 1, k));
					Real diag2 = std::abs(H(k, k)) + std::abs(H(k + 1, k + 1));
					if (off2 < tol * std::max(diag2, REAL(1.0)) || off2 < PrecisionValues<Real>::DivisionSafetyThreshold)
						result.blockSize = 1;
					else
						result.blockSize = 2; // 2x2 block (complex pair)
				}
				return result;
			}
		}

		return result;
	}

	/// @struct ComplexEigenvalue
	/// @brief Represents a potentially complex eigenvalue
	struct ComplexEigenvalue {
		Real real;
		Real imag;
		bool isComplex;

		ComplexEigenvalue(Real r = 0.0, Real i = 0.0)
			: real(r)
			, imag(i)
			, isComplex(std::abs(i) > PrecisionValues<Real>::DivisionSafetyThreshold) {}
	};

	/// @struct EigenvalueExtractionResult
	/// @brief All eigenvalues extracted from quasi-upper-triangular form
	struct EigenvalueExtractionResult {
		std::vector<ComplexEigenvalue> eigenvalues;
		int realCount;	  // Number of real eigenvalues
		int complexPairs; // Number of complex conjugate pairs
	};

	/// Extract all eigenvalues from quasi-upper-triangular (real Schur) form.
	///
	/// The matrix should be the result of QR iteration:
	/// - 1x1 diagonal blocks contain real eigenvalues
	/// - 2x2 diagonal blocks contain complex conjugate pairs
	///
	/// A 2x2 block is detected when H[i+1, i] is non-negligible.
	///
	/// @param H Quasi-upper-triangular matrix
	/// @param tol Tolerance for detecting 2x2 blocks
	/// @return EigenvalueExtractionResult with all eigenvalues
	inline EigenvalueExtractionResult ExtractEigenvalues(const Matrix<Real>& H, Real tol = 1e-10) {
		int n = H.rows();
		EigenvalueExtractionResult result;
		result.realCount = 0;
		result.complexPairs = 0;

		int i = 0;
		while (i < n) {
			if (i == n - 1) {
				// Last element: 1x1 block (real eigenvalue)
				result.eigenvalues.push_back(ComplexEigenvalue(H(i, i), 0.0));
				result.realCount++;
				i++;
			} else {
				// Check if this is a 2x2 block
				Real subdiag = std::abs(H(i + 1, i));
				Real diagSum = std::abs(H(i, i)) + std::abs(H(i + 1, i + 1));

				if (subdiag < tol * std::max(diagSum, REAL(1.0)) || subdiag < PrecisionValues<Real>::DivisionSafetyThreshold) {
					// 1x1 block: real eigenvalue
					result.eigenvalues.push_back(ComplexEigenvalue(H(i, i), REAL(0.0)));
					result.realCount++;
					i++;
				} else {
					// 2x2 block: extract eigenvalues
					auto eig = Eigenvalues2x2(H(i, i), H(i, i + 1), H(i + 1, i), H(i + 1, i + 1));

					result.eigenvalues.push_back(ComplexEigenvalue(eig.real1, eig.imag1));
					result.eigenvalues.push_back(ComplexEigenvalue(eig.real2, eig.imag2));

					if (eig.isComplex) {
						result.complexPairs++;
					} else {
						result.realCount += 2;
					}
					i += 2;
				}
			}
		}

		return result;
	}

	/// Force deflation at a specific index by zeroing the subdiagonal element.
	/// Call this after CheckDeflation returns canDeflate = true.
	///
	/// @param H Upper Hessenberg matrix (modified in place)
	/// @param deflationIndex Index returned by CheckDeflation
	inline void ApplyDeflation(Matrix<Real>& H, int deflationIndex) {
		if (deflationIndex > 0 && deflationIndex < H.rows()) {
			H(deflationIndex, deflationIndex - 1) = 0.0;
		}
	}


} // namespace MML::detail

#endif // MML_EIGEN_DETAIL_REAL_SCHUR_ANALYSIS_H
