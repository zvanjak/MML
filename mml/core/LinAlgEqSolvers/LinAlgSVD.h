///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        LinAlgSVD.h                                                         ///
///  Description: Singular Value Decomposition solver for general linear systems,     ///
///               least-squares problems, and rank-deficient matrices                 ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined  MML_LINEAR_ALG_SVD_H
#define MML_LINEAR_ALG_SVD_H

#include "MMLBase.h"

#include "base/Vector/Vector.h"
#include "base/Matrix/Matrix.h"

namespace MML
{
	/// @brief Singular Value Decomposition solver for general linear systems and least-squares
	/// @tparam Type Numeric type (Real, float, long double, etc.)
	/// @note Decomposes A=UΣVᵀ where U (m×n) and V (n×n) are orthogonal, Σ is diagonal with singular values
	/// @note Works for any m×n matrix (square, overdetermined, underdetermined)
	/// @note Complexity: O(mn²+n³) for m>n
	/// @note Provides pseudoinverse solution for rank-deficient and least-squares problems
	template<class Type>
	class SVDecompositionSolver
	{
	private:
		int m, n;
		Matrix<Type> u, v;
		Vector<Type> w;
		Type eps, tsh;

		// Helper function: numerically stable sqrt(a^2 + b^2)
		Type pythag(const Type a, const Type b) const {
			Type absa = std::abs(a), absb = std::abs(b);
			return (absa > absb ? absa * std::sqrt(1.0 + (absb / absa) * (absb / absa)) :
				(absb == 0.0 ? 0.0 : absb * std::sqrt(1.0 + (absa / absb) * (absa / absb))));
		}

		// Core SVD decomposition using Householder reduction and QR iteration
		void decompose();

		// Reorder singular values in descending order
		void reorder();

	public:
		/// @brief Get singular values vector σᵢ (descending order)
		Vector<Type> getW() const { return w; }
		/// @brief Get left singular vectors matrix U (m×n)
		Matrix<Type> getU() const { return u; }
		/// @brief Get right singular vectors matrix V (n×n)
		Matrix<Type> getV() const { return v; }

	public:
		/// @brief Constructor - performs SVD decomposition A=UΣVᵀ
		/// @param a Input matrix m×n
		/// @note Uses Householder reduction to bidiagonal form, then QR iteration
		/// @note Singular values ordered in descending order on output
		SVDecompositionSolver(const Matrix<Type>& a) : m(a.rows()), n(a.cols()), u(a), v(n, n), w(n)
		{
			// Given a Matrix a[m][n], this routine computes its singular value decomposition, A = U·W·V^T
			// The Matrix U replaces a on output (stored in u)
			// The diagonal matrix of singular values W is output as a vector w[n]
			// The matrix V (not the transpose V^T) is output as v[n][n]
			
			eps = std::numeric_limits<Type>::epsilon();
			decompose();
			reorder();
			tsh = 0.5 * std::sqrt(m + n + 1.0) * w[0] * eps;
		}

		/// @brief Compute inverse condition number (ratio of smallest to largest singular value)
		/// @return σₘᵢₙ/σₘₐₓ (0 if any singular value ≤0)
		/// @note Measures numerical stability - smaller values indicate ill-conditioning
		Type inv_condition() {
			return (w[0] <= 0. || w[n - 1] <= 0.) ? 0. : w[n - 1] / w[0];
		}

		/// @brief Solve Ax=b using pseudoinverse from SVD (handles rank deficiency)
		/// @param b Right-hand side vector (must have size m)
		/// @param x Solution vector (output, will be resized to n)
		/// @param thresh Singular value threshold (values below are treated as zero)
		///               If negative, uses default: 0.5√(m+n+1)·σ₁·εₘₐᶜₕᵢₙₑ
		/// @throws VectorDimensionError if b.size() != m
		/// @note Computes least-squares solution for overdetermined systems
		void Solve(const Vector<Type>& b, Vector<Type>& x, Type thresh = -1.)
		{
			// Solve A·x = b for a vector x using the pseudoinverse of A as obtained by SVD. If positive,
			// thresh is the threshold value below which singular values are considered as zero. If thresh is
			// negative, a default based on expected roundoff error is used.
			
			// Dimension validation
			if (b.size() != m)
				throw VectorDimensionError("SVDecompositionSolver::Solve - vector b must have size m (number of rows)", m, b.size());
			
			// Resize output vector if needed
			if (x.size() != n)
				x.Resize(n);
			
			Type tsh = (thresh >= 0. ? thresh : 0.5 * std::sqrt(m + n + 1.0) * w[0] * eps);
			
			Vector<Type> tmp(n);
			// Calculate U^T · b
			for (int j = 0; j < n; j++)
			{
				Type s = 0.0;
				if (w[j] > tsh)  // Only include non-zero singular values
				{
					for (int i = 0; i < m; i++)
						s += u[i][j] * b[i];
					s /= w[j];  // Multiply by inverse singular value
				}
				tmp[j] = s;
			}
			
			// Calculate x = V · tmp
			for (int j = 0; j < n; j++)
			{
				Type s = 0.0;
				for (int jj = 0; jj < n; jj++)
					s += v[j][jj] * tmp[jj];
				x[j] = s;
			}
		}

		Vector<Type> Solve(const Vector<Type>& b, Type thresh = -1.)
		{
			Vector<Type> x(n);
			Solve(b, x, thresh);
			return x;
		}

		// Solves m sets of n equations A·X = B using the pseudoinverse of A. The right-hand sides are
		// input as b[m][p], while x[n][p] returns the solutions. thresh as above.
		void Solve(const Matrix<Type>& b, Matrix<Type>& x, Type thresh = -1.)
		{
			int p = b.cols();
			if (b.rows() != m || x.rows() != n || x.cols() != p)
				throw MatrixDimensionError("SVD::Solve - bad dimensions", m, n, b.rows(), x.rows());
			
			Vector<Type> bcol(m), xcol(n);
			for (int j = 0; j < p; j++)
			{
				// Extract column j from b
				for (int i = 0; i < m; i++)
					bcol[i] = b[i][j];
				
				// Solve for column j
				Solve(bcol, xcol, thresh);
				
				// Store result in column j of x
				for (int i = 0; i < n; i++)
					x[i][j] = xcol[i];
			}
		}

		// Return the rank of A, after zeroing any singular values smaller than thresh. If thresh is
		// negative, a default value based on estimated roundoff is used.        
		int Rank(Type thresh = -1.) {
			Type tsh = (thresh >= 0. ? thresh : 0.5 * std::sqrt(m + n + 1.0) * w[0] * eps);
			int rank = 0;
			for (int j = 0; j < n; j++)
				if (w[j] > tsh) rank++;
			return rank;
		}

		// Return the nullity of A, after zeroing any singular values smaller than thresh. Default value as above.
		int Nullity(Type thresh = -1.) {
			Type tsh = (thresh >= 0. ? thresh : 0.5 * std::sqrt(m + n + 1.0) * w[0] * eps);
			int nullity = 0;
			for (int j = 0; j < n; j++)
				if (w[j] <= tsh) nullity++;
			return nullity;
		}

		// Gives an orthonormal basis for the range of A as the columns of a returned matrix. thresh as above.
		Matrix<Type> Range(Type thresh = -1.) {
			Type tsh = (thresh >= 0. ? thresh : 0.5 * std::sqrt(m + n + 1.0) * w[0] * eps);
			int rank = Rank(tsh);
			
			Matrix<Type> range(m, rank);
			int col = 0;
			for (int j = 0; j < n; j++)
			{
				if (w[j] > tsh)
				{
					for (int i = 0; i < m; i++)
						range[i][col] = u[i][j];
					col++;
				}
			}
			return range;
		}

		// Gives an orthonormal basis for the nullspace of A as the columns of a returned matrix. thresh as above
		Matrix<Type> Nullspace(Type thresh = -1.) {
			Type tsh = (thresh >= 0. ? thresh : 0.5 * std::sqrt(m + n + 1.0) * w[0] * eps);
			int nullity = Nullity(tsh);
			
			Matrix<Type> nullspace(n, nullity);
			int col = 0;
			for (int j = 0; j < n; j++)
			{
				if (w[j] <= tsh)
				{
					for (int i = 0; i < n; i++)
						nullspace[i][col] = v[i][j];
					col++;
				}
			}
			return nullspace;
		}
	};

	///////////////////// SVDecompositionSolver Implementation /////////////////////

	template<class Type>
	inline void SVDecompositionSolver<Type>::decompose()
	{
		bool flag;
		int i, its, j, jj, k, l, nm;
		Type anorm, c, f, g, h, s, scale, x, y, z;
		Vector<Type> rv1(n);

		g = scale = anorm = 0.0;

		// Householder reduction to bidiagonal form
		for (i = 0; i < n; i++)
		{
			l = i + 2;
			rv1[i] = scale * g;
			g = s = scale = 0.0;

			if (i < m)
			{
				for (k = i; k < m; k++) scale += std::abs(u[k][i]);
				if (scale != 0.0)
				{
					for (k = i; k < m; k++)
					{
						u[k][i] /= scale;
						s += u[k][i] * u[k][i];
					}
					f = u[i][i];
					g = -std::copysign(std::sqrt(s), f);
					h = f * g - s;
					u[i][i] = f - g;
					for (j = l - 1; j < n; j++)
					{
						for (s = 0.0, k = i; k < m; k++) s += u[k][i] * u[k][j];
						f = s / h;
						for (k = i; k < m; k++) u[k][j] += f * u[k][i];
					}
					for (k = i; k < m; k++) u[k][i] *= scale;
				}
			}

			w[i] = scale * g;
			g = s = scale = 0.0;

			if (i + 1 <= m && i + 1 != n)
			{
				for (k = l - 1; k < n; k++) scale += std::abs(u[i][k]);
				if (scale != 0.0)
				{
					for (k = l - 1; k < n; k++)
					{
						u[i][k] /= scale;
						s += u[i][k] * u[i][k];
					}
					f = u[i][l - 1];
					g = -std::copysign(std::sqrt(s), f);
					h = f * g - s;
					u[i][l - 1] = f - g;
					for (k = l - 1; k < n; k++) rv1[k] = u[i][k] / h;
					for (j = l - 1; j < m; j++)
					{
						for (s = 0.0, k = l - 1; k < n; k++) s += u[j][k] * u[i][k];
						for (k = l - 1; k < n; k++) u[j][k] += s * rv1[k];
					}
					for (k = l - 1; k < n; k++) u[i][k] *= scale;
				}
			}
			anorm = std::max(anorm, (std::abs(w[i]) + std::abs(rv1[i])));
		}

		// Accumulation of right-hand transformations
		for (i = n - 1; i >= 0; i--)
		{
			if (i < n - 1)
			{
				if (g != 0.0)
				{
					for (j = l; j < n; j++)
						v[j][i] = (u[i][j] / u[i][l]) / g;
					for (j = l; j < n; j++)
					{
						for (s = 0.0, k = l; k < n; k++) s += u[i][k] * v[k][j];
						for (k = l; k < n; k++) v[k][j] += s * v[k][i];
					}
				}
				for (j = l; j < n; j++) v[i][j] = v[j][i] = 0.0;
			}
			v[i][i] = 1.0;
			g = rv1[i];
			l = i;
		}

		// Accumulation of left-hand transformations
		for (i = std::min(m, n) - 1; i >= 0; i--)
		{
			l = i + 1;
			g = w[i];
			for (j = l; j < n; j++) u[i][j] = 0.0;
			if (g != 0.0)
			{
				g = 1.0 / g;
				for (j = l; j < n; j++)
				{
					for (s = 0.0, k = l; k < m; k++) s += u[k][i] * u[k][j];
					f = (s / u[i][i]) * g;
					for (k = i; k < m; k++) u[k][j] += f * u[k][i];
				}
				for (j = i; j < m; j++) u[j][i] *= g;
			}
			else
				for (j = i; j < m; j++) u[j][i] = 0.0;
			++u[i][i];
		}

		// Diagonalization of the bidiagonal form via QR iteration
		for (k = n - 1; k >= 0; k--)
		{
			for (its = 0; its < 30; its++)
			{
				flag = true;
				for (l = k; l >= 0; l--)
				{
					nm = l - 1;
					if (l == 0 || std::abs(rv1[l]) <= eps * anorm)
					{
						flag = false;
						break;
					}
					if (std::abs(w[nm]) <= eps * anorm) break;
				}

				if (flag)
				{
					c = 0.0;
					s = 1.0;
					for (i = l; i < k + 1; i++)
					{
						f = s * rv1[i];
						rv1[i] = c * rv1[i];
						if (std::abs(f) <= eps * anorm) break;
						g = w[i];
						h = pythag(f, g);
						w[i] = h;
						h = 1.0 / h;
						c = g * h;
						s = -f * h;
						for (j = 0; j < m; j++)
						{
							y = u[j][nm];
							z = u[j][i];
							u[j][nm] = y * c + z * s;
							u[j][i] = z * c - y * s;
						}
					}
				}

				z = w[k];
				if (l == k)
				{
					if (z < 0.0)
					{
						w[k] = -z;
						for (j = 0; j < n; j++) v[j][k] = -v[j][k];
					}
					break;
				}

				if (its == 29)
						throw ConvergenceError("SVD: no convergence in 30 iterations", 30);
				x = w[l];
				nm = k - 1;
				y = w[nm];
				g = rv1[nm];
				h = rv1[k];
				f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
				g = pythag(f, 1.0);
				f = ((x - z) * (x + z) + h * ((y / (f + std::copysign(g, f))) - h)) / x;
				c = s = 1.0;

				for (j = l; j <= nm; j++)
				{
					i = j + 1;
					g = rv1[i];
					y = w[i];
					h = s * g;
					g = c * g;
					z = pythag(f, h);
					rv1[j] = z;
					c = f / z;
					s = h / z;
					f = x * c + g * s;
					g = g * c - x * s;
					h = y * s;
					y *= c;
					for (jj = 0; jj < n; jj++)
					{
						x = v[jj][j];
						z = v[jj][i];
						v[jj][j] = x * c + z * s;
						v[jj][i] = z * c - x * s;
					}
					z = pythag(f, h);
					w[j] = z;
					if (z)
					{
						z = 1.0 / z;
						c = f * z;
						s = h * z;
					}
					f = c * g + s * y;
					x = c * y - s * g;
					for (jj = 0; jj < m; jj++)
					{
						y = u[jj][j];
						z = u[jj][i];
						u[jj][j] = y * c + z * s;
						u[jj][i] = z * c - y * s;
					}
				}
				rv1[l] = 0.0;
				rv1[k] = f;
				w[k] = x;
			}
		}
	}

	template<class Type>
	inline void SVDecompositionSolver<Type>::reorder()
	{
		int i, j, k, s, inc = 1;
		Type sw;
		Vector<Type> su(m), sv(n);

		// Shell sort to order singular values in descending order
		do { inc *= 3; inc++; } while (inc <= n);
		do {
			inc /= 3;
			for (i = inc; i < n; i++)
			{
				sw = w[i];
				for (k = 0; k < m; k++) su[k] = u[k][i];
				for (k = 0; k < n; k++) sv[k] = v[k][i];
				j = i;
				while (w[j - inc] < sw)
				{
					w[j] = w[j - inc];
					for (k = 0; k < m; k++) u[k][j] = u[k][j - inc];
					for (k = 0; k < n; k++) v[k][j] = v[k][j - inc];
					j -= inc;
					if (j < inc) break;
				}
				w[j] = sw;
				for (k = 0; k < m; k++) u[k][j] = su[k];
				for (k = 0; k < n; k++) v[k][j] = sv[k];
			}
		} while (inc > 1);

		// Flip signs for consistent sign convention
		for (k = 0; k < n; k++)
		{
			s = 0;
			for (i = 0; i < m; i++) if (u[i][k] < 0.) s++;
			for (j = 0; j < n; j++) if (v[j][k] < 0.) s++;
			if (s > (m + n) / 2)
			{
				for (i = 0; i < m; i++) u[i][k] = -u[i][k];
				for (j = 0; j < n; j++) v[j][k] = -v[j][k];
			}
		}
	}

} // namespace MML

#endif // MML_LINEAR_ALG_SVD_H
