///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        krylov_solvers_tests.cpp                                            ///
///  Description: Convergence, breakdown and preconditioner tests for the Krylov      ///
///               solvers (CG, BiCGSTAB, GMRES) and the four preconditioners          ///
///               (Jacobi, SSOR, ILU0, BlockJacobi).                                  ///
///                                                                                   ///
///  Analysis v4 (CORE C3, MinimalMathLibrary-wewi.2.3): BiCGSTAB had zero tests,     ///
///  GMRES one wiring test, preconditioner convergence was untested.                  ///
///////////////////////////////////////////////////////////////////////////////////////////

#include <catch2/catch_all.hpp>

#include <mml/base/SparseMatrix/SparseMatrix.h>
#include <mml/core/SparseSolvers/IterativeSolvers.h>

#include <cmath>
#include <memory>
#include <random>
#include <vector>

namespace MML::Tests::Core::KrylovSolversTests
{
    using SparseMatrix::SparseMatrixCOO;
    using SparseMatrix::SparseMatrixCSR;
    using namespace SparseSolvers;

    namespace
    {
			double Norm2(const Vector<double>& v) {
				double s = 0.0;
				for (double x : v)
					s += x * x;
				return std::sqrt(s);
			}

			double ResidualNorm(const SparseMatrixCSR<double>& A, const Vector<double>& x, const Vector<double>& b) {
				Vector<double> Ax = A * x;
				double s = 0.0;
				for (int i = 0; i < b.size(); ++i)
					s += (b[i] - Ax[i]) * (b[i] - Ax[i]);
				return std::sqrt(s);
			}

			double MaxAbsDiff(const Vector<double>& a, const Vector<double>& b) {
				double m = 0.0;
				for (int i = 0; i < a.size(); ++i)
					m = std::max(m, std::abs(a[i] - b[i]));
				return m;
			}

				// 1D convection-diffusion  -u'' + c u' = f  with central differences: nonsymmetric, diagonally dominant for |c h| < 2
        SparseMatrixCSR<double> ConvectionDiffusion1D(int n, double peclet)
        {
            SparseMatrixCOO<double> coo(n, n);
            for (int i = 0; i < n; ++i) {
                coo.addEntry(i, i, 2.0);
                if (i > 0)     coo.addEntry(i, i - 1, -1.0 - peclet / 2.0);
                if (i + 1 < n) coo.addEntry(i, i + 1, -1.0 + peclet / 2.0);
            }
            return SparseMatrixCSR<double>(coo);
        }

        // Random sparse strictly diagonally dominant nonsymmetric matrix (always nonsingular)
        SparseMatrixCSR<double> RandomDominant(int n, int perRow, unsigned seed)
        {
            std::mt19937 gen(seed);
            std::uniform_real_distribution<double> u(-1.0, 1.0);
            std::uniform_int_distribution<int> col(0, n - 1);
            SparseMatrixCOO<double> coo(n, n);
            for (int i = 0; i < n; ++i) {
                double offSum = 0.0;
                for (int k = 0; k < perRow; ++k) {
                    int j = col(gen);
                    if (j == i) continue;
                    double v = u(gen);
                    coo.addEntry(i, j, v);
                    offSum += std::abs(v);
                }
                coo.addEntry(i, i, offSum + 1.0 + std::abs(u(gen)));
            }
            return SparseMatrixCSR<double>(coo);
        }

				Vector<double> KnownSolution(int n) {
					Vector<double> x(n);
					for (int i = 0; i < n; ++i)
						x[i] = std::sin(0.37 * i) + 0.5;
					return x;
				}

				SolverConfig<double> Tight(int maxIter = 2000)
        {
            SolverConfig<double> c;
            c.setTolerance(1e-12).setAbsTolerance(1e-14).setMaxIterations(maxIter);
            return c;
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    //                                   BiCGSTAB
    ///////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("BiCGSTAB converges on SPD, nonsymmetric and random dominant systems with a consistent residual report",
              "[SparseSolvers][BiCGSTAB][wewi.2.3]")
    {
        SECTION("2D Laplacian 12x12 (SPD): agrees with CG")
        {
            auto A = SparseMatrix::Laplacian2D<double>(12, 12);
            const int n = A.rows();
            auto xExact = KnownSolution(n);
            auto b = A * xExact;

						Vector<double> xB(n, 0.0), xC(n, 0.0);
						auto rb = solveBiCGSTAB(A, b, xB, Tight());
            auto rc = solveCG(A, b, xC, Tight());

						INFO(rb.error_message);
						REQUIRE(rb.converged);
						REQUIRE(rc.converged);
						REQUIRE(MaxAbsDiff(xB, xExact) < 1e-8);
            REQUIRE(MaxAbsDiff(xB, xC) < 1e-8);
            REQUIRE(rb.relativeResidual < 1e-12);
            // the reported residual must describe the returned x
            REQUIRE(ResidualNorm(A, xB, b) <= 10 * std::max(rb.residualNorm, 1e-14));
        }

        SECTION("1D convection-diffusion (nonsymmetric), n = 200, Peclet 1.0")
        {
            auto A = ConvectionDiffusion1D(200, 1.0);
            auto xExact = KnownSolution(200);
            auto b = A * xExact;
						Vector<double> x(200, 0.0);

						auto r = solveBiCGSTAB(A, b, x, Tight(5000));
						INFO(r.error_message << " iterations " << r.iterations_used);
						REQUIRE(r.converged);
						REQUIRE(MaxAbsDiff(x, xExact) < 1e-7);
            REQUIRE(ResidualNorm(A, x, b) < 1e-9 * Norm2(b));
        }

        SECTION("random strictly diagonally dominant nonsymmetric, n = 300")
        {
            auto A = RandomDominant(300, 4, 123);
            auto xExact = KnownSolution(300);
            auto b = A * xExact;
						Vector<double> x(300, 0.0);

						auto r = solveBiCGSTAB(A, b, x, Tight());
						INFO(r.error_message << " iterations " << r.iterations_used);
						REQUIRE(r.converged);
						REQUIRE(r.iterations_used < 300); // far fewer than n for a dominant matrix
						REQUIRE(MaxAbsDiff(x, xExact) < 1e-9);
        }

        SECTION("non-zero initial guess is used (fewer iterations from a near-solution)")
        {
            auto A = RandomDominant(200, 4, 7);
            auto xExact = KnownSolution(200);
            auto b = A * xExact;

						Vector<double> xZero(200, 0.0);
						auto rZero = solveBiCGSTAB(A, b, xZero, Tight());
						Vector<double> xNear = xExact;
						for (auto& v : xNear) v += 1e-6;
            auto rNear = solveBiCGSTAB(A, b, xNear, Tight());

						REQUIRE(rZero.converged);
						REQUIRE(rNear.converged);
						REQUIRE(rNear.iterations_used <= rZero.iterations_used);
						REQUIRE(MaxAbsDiff(xNear, xExact) < 1e-9);
        }
    }

    TEST_CASE("BiCGSTAB reports breakdown, max-iterations and invalid input honestly", "[SparseSolvers][BiCGSTAB][wewi.2.3]")
    {
        SECTION("(r0, A r0) = 0 at the first step is a breakdown, not a Success")
        {
            // A = [[0,1],[1,0]], b = (1,0): r0 = b, v = A r0 = (0,1), (r0, v) = 0
            SparseMatrixCOO<double> coo(2, 2);
            coo.addEntry(0, 1, 1.0);
            coo.addEntry(1, 0, 1.0);
            SparseMatrixCSR<double> A(coo);
						Vector<double> b{1.0, 0.0}, x(2, 0.0);

						auto r = solveBiCGSTAB(A, b, x, Tight());
						REQUIRE_FALSE(r.converged);
						REQUIRE(r.status == AlgorithmStatus::NumericalInstability);
						REQUIRE_FALSE(r.error_message.empty());
				}

        SECTION("maxIterations exhausted -> MaxIterations with the current (non-converged) residual")
        {
            auto A = SparseMatrix::Laplacian2D<double>(20, 20);
						Vector<double> b(A.rows(), 1.0), x(A.rows(), 0.0);
						SolverConfig<double> c;
            c.setTolerance(1e-14).setAbsTolerance(0.0).setMaxIterations(2);

            auto r = solveBiCGSTAB(A, b, x, c);
						REQUIRE_FALSE(r.converged);
						REQUIRE(r.status == AlgorithmStatus::MaxIterationsExceeded);
						REQUIRE(r.iterations_used == 2);
						REQUIRE(r.relativeResidual > 1e-14);
            REQUIRE(ResidualNorm(A, x, b) <= 10 * r.residualNorm);
        }

        SECTION("dimension mismatch -> InvalidInput; zero RHS -> zero solution")
        {
            auto A = SparseMatrix::Laplacian1D<double>(5);
						Vector<double> bBad(4, 1.0), x(5, 0.0);
						REQUIRE(solveBiCGSTAB(A, bBad, x, Tight()).status == AlgorithmStatus::InvalidInput);

						Vector<double> bZero(5, 0.0), x0(5, 3.0);
						auto r = solveBiCGSTAB(A, bZero, x0, Tight());
						REQUIRE(r.converged);
						REQUIRE(r.iterations_used == 0);
						for (double v : x0) REQUIRE(v == 0.0);
        }

        SECTION("progress callback is invoked once per iteration")
        {
            auto A = ConvectionDiffusion1D(50, 0.5);
						Vector<double> b(50, 1.0), x(50, 0.0);
						BiCGSTAB<double> solver;
            solver.setConfig(Tight());
            int calls = 0;
            double last = std::numeric_limits<double>::infinity();
            solver.setCallback([&](int, double res) { ++calls; last = res; });
            auto r = solver.solve(A, b, x);
						REQUIRE(r.converged);
						REQUIRE(calls >= r.iterations_used); // initial report + one per iteration
						REQUIRE(std::isfinite(last));
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    //                                     GMRES
    ///////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("GMRES: full and restarted variants converge on nonsymmetric systems", "[SparseSolvers][GMRES][wewi.2.3]")
    {
        SECTION("full GMRES (restart = n) terminates within n iterations on a 40x40 nonsymmetric system")
        {
            const int n = 40;
            auto A = RandomDominant(n, 3, 99);
            auto xExact = KnownSolution(n);
            auto b = A * xExact;
						Vector<double> x(n, 0.0);

						auto r = solveGMRES(A, b, x, n, Tight(n));
						INFO(r.error_message << " iterations " << r.iterations_used);
						REQUIRE(r.converged);
						REQUIRE(r.iterations_used <= n);
						REQUIRE(MaxAbsDiff(x, xExact) < 1e-9);
        }

        SECTION("restarted GMRES(5) and GMRES(50) both converge on convection-diffusion n = 150; full GMRES in <= n")
        {
            // Note: iteration counts are deliberately NOT compared between restart lengths.
            // Restarted GMRES is not monotone in m (Embree, "The tortoise and the hare restart
            // GMRES", SIAM Review 45, 2003); on this matrix GMRES(5) beats GMRES(50) (458 vs 952).
            auto A = ConvectionDiffusion1D(150, 0.8);
            auto xExact = KnownSolution(150);
            auto b = A * xExact;

						Vector<double> x5(150, 0.0), x50(150, 0.0), xFull(150, 0.0);
						auto r5 = solveGMRES(A, b, x5, 5, Tight(20000));
            auto r50 = solveGMRES(A, b, x50, 50, Tight(20000));
            auto rFull = solveGMRES(A, b, xFull, 150, Tight(20000));
						INFO("GMRES(5): " << r5.iterations_used << ", GMRES(50): " << r50.iterations_used << ", full: " << rFull.iterations_used);
						REQUIRE(r5.converged);
						REQUIRE(r50.converged);
						REQUIRE(rFull.converged);
						REQUIRE(MaxAbsDiff(x5, xExact) < 1e-7);
            REQUIRE(MaxAbsDiff(x50, xExact) < 1e-7);
            REQUIRE(MaxAbsDiff(xFull, xExact) < 1e-9);
						REQUIRE(rFull.iterations_used <= 150);
						REQUIRE(rFull.iterations_used <= r5.iterations_used);
						REQUIRE(rFull.iterations_used <= r50.iterations_used);
				}

        SECTION("residual norms reported through the callback are non-increasing within a cycle")
        {
            auto A = RandomDominant(60, 3, 5);
						Vector<double> b(60, 1.0), x(60, 0.0);
						GMRES<double> solver(60);
            solver.setConfig(Tight(60));
						Vector<double> history;
						solver.setCallback([&](int, double res) { history.push_back(res); });
            auto r = solver.solve(A, b, x);
						REQUIRE(r.converged);
						REQUIRE(history.size() >= 2);
            for (int i = 1; i < history.size(); ++i)
                REQUIRE(history[i] <= history[i - 1] * (1.0 + 1e-12));
        }

        SECTION("inconsistent singular system is not reported as Success")
        {
            // A = diag(1, 0), b = (1, 1): no solution
            auto A = SparseMatrix::DiagonalCSR(std::vector<double>{1.0, 0.0});
						Vector<double> b{1.0, 1.0}, x(2, 0.0);
						auto r = solveGMRES(A, b, x, 2, Tight(10));
						REQUIRE_FALSE(r.converged);
				}

        SECTION("dimension mismatch -> InvalidInput")
        {
            auto A = SparseMatrix::Laplacian1D<double>(5);
						Vector<double> bBad(6, 1.0), x(5, 0.0);
						REQUIRE(solveGMRES(A, bBad, x, 5, Tight()).status == AlgorithmStatus::InvalidInput);
				}
    }

    ///////////////////////////////////////////////////////////////////////////////////////
    //                                 Preconditioners
    ///////////////////////////////////////////////////////////////////////////////////////

    TEST_CASE("Preconditioners apply M^-1 correctly on cases with known inverses", "[SparseSolvers][Preconditioner][wewi.2.3]")
    {
			const Vector<double> r{1.0, -2.0, 3.0, 0.5, -1.5};
			Vector<double> z;

			SECTION("Jacobi: z = r / diag") {
				auto A = SparseMatrix::DiagonalCSR(std::vector<double>{2.0, 4.0, -8.0, 0.5, 10.0});
				JacobiPreconditioner<double> M;
				M.setup(A);
				M.apply(r, z);
				const Vector<double> expected{0.5, -0.5, -0.375, 1.0, -0.15};
				REQUIRE(MaxAbsDiff(z, expected) < 1e-15);
        }

        SECTION("SSOR on a diagonal matrix: z = omega r / diag (documented (2-omega)-scaled SSOR)")
        {
            auto A = SparseMatrix::DiagonalCSR(std::vector<double>{2.0, 4.0, -8.0, 0.5, 10.0});
            SSORPreconditioner<double> M(1.5);
            M.setup(A);
            M.apply(r, z);
						const Vector<double> expected{1.5 * 0.5, 1.5 * -0.5, 1.5 * -0.375, 1.5 * 1.0, 1.5 * -0.15};
						REQUIRE(MaxAbsDiff(z, expected) < 1e-15);
        }

        SECTION("SSOR is symmetric for symmetric A: (u, M^-1 v) = (v, M^-1 u)")
        {
            auto A = SparseMatrix::Laplacian2D<double>(4, 4);
            SSORPreconditioner<double> M(1.3);
            M.setup(A);
						Vector<double> u(16), v(16), Mu, Mv;
						for (int i = 0; i < 16; ++i) { u[i] = std::cos(0.7 * i); v[i] = std::sin(1.3 * i + 0.2); }
            M.apply(u, Mu);
            M.apply(v, Mv);
            double uMv = 0.0, vMu = 0.0;
            for (int i = 0; i < 16; ++i) { uMv += u[i] * Mv[i]; vMu += v[i] * Mu[i]; }
            REQUIRE(uMv == Catch::Approx(vMu).epsilon(1e-12));
        }

        SECTION("ILU(0) on a tridiagonal matrix is the exact LU: z = A^-1 r")
        {
            auto A = SparseMatrix::Laplacian1D<double>(5);
            ILU0Preconditioner<double> M;
            M.setup(A);
            M.apply(r, z);
            REQUIRE(ResidualNorm(A, z, r) < 1e-13);
        }

        SECTION("Block Jacobi with one block covering the matrix is the exact inverse")
        {
            auto A = ConvectionDiffusion1D(5, 0.6);
            BlockJacobiPreconditioner<double> M(5);
            M.setup(A);
            M.apply(r, z);
            REQUIRE(ResidualNorm(A, z, r) < 1e-13);
        }

        SECTION("Block Jacobi with blockSize 1 equals Jacobi")
        {
            auto A = ConvectionDiffusion1D(5, 0.6);
            BlockJacobiPreconditioner<double> MB(1);
            JacobiPreconditioner<double> MJ;
            MB.setup(A);
            MJ.setup(A);
						Vector<double> zb, zj;
						MB.apply(r, zb);
            MJ.apply(r, zj);
            REQUIRE(MaxAbsDiff(zb, zj) < 1e-15);
        }
    }

    TEST_CASE("Preconditioned Krylov solvers converge and the strong preconditioners reduce iteration counts",
              "[SparseSolvers][Preconditioner][CG][BiCGSTAB][GMRES][wewi.2.3]")
    {
        auto A = SparseMatrix::Laplacian2D<double>(25, 25);   // n = 625, condition number ~ 270
        const int n = A.rows();
        auto xExact = KnownSolution(n);
        auto b = A * xExact;
        auto cfg = Tight(5000);

				Vector<double> xPlain(n, 0.0);
				auto plainCG = solveCG(A, b, xPlain, cfg);
				REQUIRE(plainCG.converged);
				REQUIRE(MaxAbsDiff(xPlain, xExact) < 1e-8);

        SECTION("CG + Jacobi (constant diagonal: same iteration count, still correct)")
        {
            auto M = std::make_shared<JacobiPreconditioner<double>>(); M->setup(A);
						Vector<double> x(n, 0.0);
						auto r = solveCG(A, b, x, M, cfg);
						REQUIRE(r.converged);
						REQUIRE(MaxAbsDiff(x, xExact) < 1e-8);
						REQUIRE(std::abs(r.iterations_used - plainCG.iterations_used) <= 2);
				}

        SECTION("CG + SSOR converges in fewer iterations than plain CG")
        {
            auto M = std::make_shared<SSORPreconditioner<double>>(1.0); M->setup(A);
						Vector<double> x(n, 0.0);
						auto r = solveCG(A, b, x, M, cfg);
						INFO("plain " << plainCG.iterations_used << " SSOR " << r.iterations_used);
						REQUIRE(r.converged);
						REQUIRE(MaxAbsDiff(x, xExact) < 1e-8);
						REQUIRE(r.iterations_used < plainCG.iterations_used);
				}

        SECTION("CG + Block Jacobi (block 25 = one grid line) converges in fewer iterations than plain CG")
        {
            auto M = std::make_shared<BlockJacobiPreconditioner<double>>(25); M->setup(A);
						Vector<double> x(n, 0.0);
						auto r = solveCG(A, b, x, M, cfg);
						INFO("plain " << plainCG.iterations_used << " BlockJacobi " << r.iterations_used);
						REQUIRE(r.converged);
						REQUIRE(MaxAbsDiff(x, xExact) < 1e-8);
						REQUIRE(r.iterations_used < plainCG.iterations_used);
				}

        SECTION("BiCGSTAB + ILU0 converges in fewer iterations than plain BiCGSTAB")
        {
					Vector<double> xP(n, 0.0), xI(n, 0.0);
					auto plain = solveBiCGSTAB(A, b, xP, cfg);
					auto M = std::make_shared<ILU0Preconditioner<double>>();
					M->setup(A);
					auto ilu = solveBiCGSTAB(A, b, xI, M, cfg);
					INFO("plain " << plain.iterations_used << " ILU0 " << ilu.iterations_used);
					REQUIRE(plain.converged);
					REQUIRE(ilu.converged);
					REQUIRE(MaxAbsDiff(xI, xExact) < 1e-8);
					REQUIRE(ilu.iterations_used < plain.iterations_used);
				}

        SECTION("GMRES(20) + ILU0 converges in fewer iterations than plain GMRES(20)")
        {
					Vector<double> xP(n, 0.0), xI(n, 0.0);
					auto plain = solveGMRES(A, b, xP, 20, cfg);
					auto M = std::make_shared<ILU0Preconditioner<double>>();
					M->setup(A);
					auto ilu = solveGMRES(A, b, xI, M, 20, cfg);
					INFO("plain " << plain.iterations_used << " ILU0 " << ilu.iterations_used);
					REQUIRE(plain.converged);
					REQUIRE(ilu.converged);
					REQUIRE(MaxAbsDiff(xI, xExact) < 1e-7);
					REQUIRE(ilu.iterations_used < plain.iterations_used);
				}
    }

    TEST_CASE("CG and BiCGSTAB do not report false breakdown near convergence at tight tolerances (regression)",
              "[SparseSolvers][CG][BiCGSTAB][wewi.2.3]")
    {
			// Before the fix, (p,Ap), rho and (t,t) were compared with an absolute 100*eps: these
			// quantities are quadratic in the residual, so every solve with tolerance ~1e-12 on a
			// system larger than a handful of unknowns ended in AlgorithmStatus::NumericalInstability.
			auto A = SparseMatrix::Laplacian2D<double>(30, 30); // n = 900
			const int n = A.rows();
			Vector<double> b(n, 1e-3), x(n, 0.0); // small RHS makes the old thresholds bite early
			SolverConfig<double> c;
			c.setTolerance(1e-13).setAbsTolerance(0.0).setMaxIterations(10000);

			auto rc = solveCG(A, b, x, c);
			INFO("CG: " << rc.error_message);
			REQUIRE(rc.status == AlgorithmStatus::Success);
			REQUIRE(ResidualNorm(A, x, b) < 1e-13 * Norm2(b) * 1.01);

			std::fill(x.begin(), x.end(), 0.0);
			auto rb = solveBiCGSTAB(A, b, x, c);
			INFO("BiCGSTAB: " << rb.error_message);
			REQUIRE(rb.status == AlgorithmStatus::Success);
			REQUIRE(ResidualNorm(A, x, b) < 1e-13 * Norm2(b) * 10);
    }

    TEST_CASE("solvePoisson2D reproduces a manufactured solution", "[SparseSolvers][Poisson][wewi.2.3]")
    {
        // -Laplace u = f on the unit square, u = sin(pi x) sin(pi y): f = 2 pi^2 u
        const int nx = 33, ny = 33;
        const double h = 1.0 / (nx - 1);
				Vector<double> f, uExact;
				for (int j = 1; j < ny - 1; ++j)
            for (int i = 1; i < nx - 1; ++i) {
                const double x = i * h, y = j * h;
                const double u = std::sin(Constants::PI * x) * std::sin(Constants::PI * y);
                uExact.push_back(u);
                f.push_back(2.0 * Constants::PI * Constants::PI * u);
            }
        auto u = solvePoisson2D<double>(nx, ny, f, Tight(5000));
        REQUIRE(u.size() == uExact.size());
        // second-order discretisation error ~ h^2 pi^2 / 12 * max|u| ~ 8e-4
        REQUIRE(MaxAbsDiff(u, uExact) < 2e-3);
    }
}
