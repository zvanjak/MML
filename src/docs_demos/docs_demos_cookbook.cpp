///////////////////////////////////////////////////////////////////////////////////////////
// docs_demos_cookbook.cpp - runnable counterpart of docs/COOKBOOK.md
// ALL cookbook recipes live in this one file; every recipe snippet in COOKBOOK.md
// must compile and execute here (docs<->demos contract, AGENTS.md).
///////////////////////////////////////////////////////////////////////////////////////////
#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>
#include <mml/base/Function.h>
#include <mml/base/ComplexFunction.h>
#include <mml/base/Algebra_base.h>
#include <mml/base/ChebyshevApproximation.h>
#include <mml/base/Combinatorics.h>
#include <mml/base/InterpolatedFunction.h>
#include <mml/base/Polynom.h>
#include <mml/base/Geometry/Geometry3D.h>
#include <mml/base/Geometry/Geometry3DBodies.h>
#include <mml/base/Graph.h>
#include <mml/base/NumberTheory.h>
#include <mml/base/Quaternions.h>
#include <mml/base/ODESystem.h>
#include <mml/interfaces/IODESystemWithEvents.h>
#include <mml/core/LinAlgEqSolvers/LinAlgDirect.h>
#include <mml/core/ComplexAnalysis.h>
#include <mml/core/Derivation.h>
#include <mml/core/Curves.h>
#include <mml/core/Surfaces.h>
#include <mml/core/CoordTransf/CoordTransfBase.h>
#include <mml/core/CoordTransf/CoordTransfSpherical.h>
#include <mml/core/CoordTransf/CoordTransfCylindrical.h>
#include <mml/core/Fields/Fields.h>
#include <mml/core/Fields/ScalarFieldOperations.h>
#include <mml/core/Fields/VectorFieldOperations.h>
#include <mml/core/MetricTensor.h>
#include <mml/base/Tensor/Tensor2.h>
#include <mml/base/Tensor/Tensor4.h>
#include <mml/core/DifferentialGeometry/Frames.h>
#include <mml/core/Integration.h>
#include <mml/core/Integration/MonteCarloIntegration.h>
#include <mml/core/Integration/PathIntegration.h>
#include <mml/core/Integration/SurfaceIntegration.h>
#include <mml/systems/LinearSystem.h>
#include <mml/algorithms/EigenSystemSolvers.h>
#include <mml/algorithms/MatrixAlg.h>
#include <mml/algorithms/Analyzers/MatrixAnalyzer.h>
#include <mml/algorithms/Analyzers/FunctionsAnalyzer.h>
#include <mml/algorithms/Fourier/Fourier.h>
#include <mml/algorithms/Fourier/FourierConvolution.h>
#include <mml/algorithms/Fourier/FourierRealFFT.h>
#include <mml/algorithms/Fourier/FourierSpectrum.h>
#include <mml/algorithms/Fourier/FourierWindowing.h>
#include <mml/algorithms/GraphAlg.h>
#include <mml/algorithms/ComputationalGeometry.h>
#include <mml/algorithms/CompGeometry/KDTree.h>
#include <mml/algorithms/Statistics.h>
#include <mml/algorithms/Statistics/Distributions.h>
#include <mml/algorithms/Statistics/Histogram.h>
#include <mml/algorithms/RootFinding.h>
#include <mml/algorithms/RootFinding/RootFindingPolynoms.h>
#include <mml/algorithms/RootFinding/RootFindingComplex.h>
#include <mml/algorithms/DAESolvers.h>
#include <mml/algorithms/ODESolvers.h>
#include <mml/algorithms/ODESolvers/ODESolverEventDetection.h>
#include <mml/algorithms/CurveFitting.h>
#include <mml/algorithms/Optimization/Optimization.h>
#include <mml/algorithms/Optimization/OptimizationMultidim.h>
#include <mml/algorithms/Optimization/LinearProgramming.h>
#include <mml/algorithms/Geodesic.h>
#include <mml/tools/DataLoader.h>
#endif

#include <array>
#include <filesystem>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace MML;

namespace
{
	using F2 = Algebra::PrimeFieldElement<2>;

	struct Cookbook_GF4_Modulus
	{
		static Algebra::Polynomial<F2> modulus()
		{
			return Algebra::Polynomial<F2>({ F2(1), F2(1), F2(1) });
		}
	};

	VectorN<Real, 3> Cookbook_CustomCurvePoint(Real t)
	{
		return VectorN<Real, 3>{ t, std::sin(t), Real{0.25} * t * t };
	}

	class Cookbook_WavySheet : public Surfaces::ISurfaceCartesian
	{
	public:
		Real getMinU() const override { return -2.0; }
		Real getMaxU() const override { return 2.0; }
		Real getMinW() const override { return -2.0; }
		Real getMaxW() const override { return 2.0; }

		VectorN<Real, 3> operator()(Real u, Real w) const override
		{
			return VectorN<Real, 3>{ u, w, Real{0.20} * std::sin(u) * std::cos(w) };
		}
	};

	class Cookbook_LinearDAE : public IODESystemDAEWithJacobian
	{
	public:
		int getDiffDim() const override { return 1; }
		int getAlgDim() const override { return 1; }

		void diffEqs(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& dxdt) const override
		{
			dxdt[0] = -x[0] + y[0];
		}

		void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& g) const override
		{
			g[0] = x[0] + y[0] - REAL(1.0);
		}

		std::string getDiffVarName(int) const override { return "x"; }
		std::string getAlgVarName(int) const override { return "y"; }

		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& df_dx) const override
		{
			df_dx(0, 0) = REAL(-1.0);
		}

		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& df_dy) const override
		{
			df_dy(0, 0) = REAL(1.0);
		}

		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& dg_dx) const override
		{
			dg_dx(0, 0) = REAL(1.0);
		}

		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& dg_dy) const override
		{
			dg_dy(0, 0) = REAL(1.0);
		}
	};

	class Cookbook_TimedDAEEvent : public IODESystemDAEWithEvents
	{
		Real _eventTime;

	public:
		explicit Cookbook_TimedDAEEvent(Real eventTime)
			: _eventTime(eventTime) {}

		int getDiffDim() const override { return 1; }
		int getAlgDim() const override { return 1; }

		void diffEqs(Real, const Vector<Real>&, const Vector<Real>&, Vector<Real>& dxdt) const override
		{
			dxdt[0] = REAL(1.0);
		}

		void algConstraints(Real, const Vector<Real>& x, const Vector<Real>& y, Vector<Real>& constraints) const override
		{
			constraints[0] = x[0] + y[0] - REAL(1.0);
		}

		void jacobian_fx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& value) const override
		{
			value(0, 0) = REAL(0.0);
		}

		void jacobian_fy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& value) const override
		{
			value(0, 0) = REAL(0.0);
		}

		void jacobian_gx(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& value) const override
		{
			value(0, 0) = REAL(1.0);
		}

		void jacobian_gy(Real, const Vector<Real>&, const Vector<Real>&, Matrix<Real>& value) const override
		{
			value(0, 0) = REAL(1.0);
		}

		int getNumEvents() const override { return 1; }

		Real eventFunction(int, Real t, const Vector<Real>&, const Vector<Real>&) const override
		{
			return t - _eventTime;
		}

		EventDirection getEventDirection(int) const override { return EventDirection::Increasing; }
		EventAction getEventAction(int) const override { return EventAction::Restart; }

		void handleEvent(int, Real, Vector<Real>& x, Vector<Real>& y) const override
		{
			x[0] = REAL(0.25);
			y[0] = REAL(42.0);
		}
	};

	std::string Cookbook_FormatPoint(const Pnt3Cart& point)
	{
		std::ostringstream oss;
		oss << "(" << point.X() << "," << point.Y() << "," << point.Z() << ")";
		return oss.str();
	}

	std::string Cookbook_FormatLineIntersectionType(LineIntersectionType3D type)
	{
		switch (type) {
		case LineIntersectionType3D::Point: return "point";
		case LineIntersectionType3D::Parallel: return "parallel";
		case LineIntersectionType3D::Coincident: return "coincident";
		case LineIntersectionType3D::Skew: return "skew";
		}
		return "unknown";
	}

	template<class Container>
	std::string Cookbook_FormatIntContainer(const Container& values)
	{
		std::ostringstream text;
		text << "[";
		for (std::size_t index = 0; index < values.size(); ++index) {
			if (index > 0) text << ",";
			text << values[index];
		}
		text << "]";
		return text.str();
	}

	template<int Modulus>
	std::string Cookbook_FormatModInt(Algebra::ModInt<Modulus> value)
	{
		return std::to_string(value.value());
	}

	template<class FieldElement>
	std::string Cookbook_FormatGFElement(const FieldElement& value)
	{
		std::ostringstream text;
		text << "[";
		const auto& coefficients = value.coefficients();
		for (std::size_t index = 0; index < coefficients.size(); ++index) {
			if (index > 0) text << ",";
			text << coefficients[index].value();
		}
		text << "]";
		return text.str();
	}

	std::vector<std::vector<int>> Cookbook_BinaryColorings(int count)
	{
		std::vector<std::vector<int>> colorings;
		const int total = 1 << count;
		colorings.reserve(static_cast<std::size_t>(total));
		for (int mask = 0; mask < total; ++mask) {
			std::vector<int> coloring(static_cast<std::size_t>(count));
			for (int position = 0; position < count; ++position)
				coloring[static_cast<std::size_t>(position)] = (mask >> position) & 1;
			colorings.push_back(std::move(coloring));
		}
		return colorings;
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 1: Solving linear systems
	///////////////////////////////////////////////////////////////////////////

	// 1.1 Real 3x3 system - LU direct, automatic selection, LinearSystem facade
	void Cookbook_Recipe01_1_RealSystem(std::ostream& out)
	{
		// --- Setup: the classic 3x3 system, exact solution x = (2, 3, -1)
		Matrix<Real> A(3, 3, {  2,  1, -1,
		                       -3, -1,  2,
		                       -2,  1,  2 });
		Vector<Real> b({ 8, -11, -3 });

		// --- Solve: directly with an explicit solver (decompose once, reuse for many b)
		LUSolver<Real> lu(A);              // decomposition happens here; throws if singular
		Vector<Real> x = lu.Solve(b);

		// --- Solve: automatic solver selection, the one-liner
		Vector<Real> x2 = Systems::SolveLinearSystem(A, b);

		// --- Solve: the LinearSystem facade with diagnostics and control
		Systems::LinearSystem<Real> sys(A, b);

		Vector<Real> x3 = sys.Solve();       // auto-selects the best solver
		auto check = sys.Verify(x3);         // residual-based quality check

		Vector<Real> x4 = sys.SolveByQR();   // or force a specific method

		out << "  [1.1] LUSolver:           x = " << x << "\n";
		out << "  [1.1] SolveLinearSystem:  x = " << x2 << "\n";
		out << "  [1.1] LinearSystem+Verify: x = " << x3
		    << ", isAccurate = " << std::boolalpha << check.isAccurate << "\n";
		out << "  [1.1] SolveByQR:           x = " << x4 << "\n";
	}

	// 1.2 Complex 3x3 system - same recipe, Complex elements
	void Cookbook_Recipe01_2_ComplexSystem(std::ostream& out)
	{
		// --- Setup: manufacture b from a known solution x = (1, i, 1-i)
		Matrix<Complex> Ac(3, 3, { Complex(2, 1), Complex(1, 0),  Complex(0, -1),
		                           Complex(0, 2), Complex(3, -1), Complex(1, 0),
		                           Complex(1, 0), Complex(0, 1),  Complex(2, 2) });

		Vector<Complex> xExpected({ Complex(1, 0), Complex(0, 1), Complex(1, -1) });
		Vector<Complex> bc = Ac * xExpected;

		// --- Solve: identical API, complex element type
		LUSolver<Complex> luc(Ac);
		Vector<Complex> xc = luc.Solve(bc);

		out << "  [1.2] LUSolver<Complex>:  x = " << xc
		    << ", matches expected = " << std::boolalpha
		    << xc.IsEqualTo(xExpected, 1e-12) << "\n";
	}

	void Cookbook_Recipe01_SolvingLinearSystems(std::ostream& out)
	{
		Cookbook_Recipe01_1_RealSystem(out);
		Cookbook_Recipe01_2_ComplexSystem(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 2: Root finding
	///////////////////////////////////////////////////////////////////////////

	// 2.1 Classic: bracket the root first, then refine with Brent
	void Cookbook_Recipe02_1_BracketThenSolve(std::ostream& out)
	{
		// --- Setup: Kepler's equation E - e*sin(E) = M for a highly eccentric orbit
		const Real ecc = 0.9, M = 1.0;
		RealFunction kepler([](Real E) { return E - 0.9 * std::sin(E) - 1.0; });

		// --- Solve: expand an initial interval until it brackets a sign change...
		Real x1 = 0.0, x2 = 0.5;
		bool bracketed = RootFinding::BracketRoot(kepler, x1, x2);

		// ...then hand the bracket to a robust solver
		RootFinding::RootFindingConfig config;
		config.tolerance = 1e-14;
		auto result = RootFinding::FindRootBrent(kepler, x1, x2, config);

		out << "  [2.1] bracketed [" << x1 << ", " << x2 << "]: " << std::boolalpha << bracketed
		    << ", eccentric anomaly E = " << result.root
		    << " (" << result.iterations_used << " iters, |f| = "
		    << std::abs(result.function_value) << ")\n";
		(void)ecc; (void)M;
	}

	// 2.2 All roots in an interval
	void Cookbook_Recipe02_2_AllRoots(std::ostream& out)
	{
		// --- Setup: cantilever-beam frequency equation cos(x)*cosh(x) + 1 = 0
		RealFunction beam([](Real x) { return std::cos(x) * std::cosh(x) + 1.0; });

		// --- Solve: isolate + refine every root in [0, 12]
		auto all = RootFinding::FindAllRealRootsInInterval(beam, 0.0, 12.0);

		out << "  [2.2] found " << all.roots.size() << " beam frequencies:";
		for (Real r : all.roots)
			out << " " << r;
		out << "\n";
	}

	// 2.3 Polynomial roots
	void Cookbook_Recipe02_3_PolynomialRoots(std::ostream& out)
	{
		// --- Setup: p(x) = (x-1)(x-2)(x-3)(x-4)(x-5), ascending coefficients
		PolynomReal p({ -120, 274, -225, 85, -15, 1 });

		// --- Solve: Laguerre's method returns all (complex) roots at once
		Vector<Complex> roots = RootFinding::LaguerreRoots(p);

		out << "  [2.3] roots of Wilkinson-5:";
		for (int i = 0; i < roots.size(); i++)
			out << " " << roots[i].real();
		out << "\n";
	}

	// 2.4 Roots of a complex function
	void Cookbook_Recipe02_4_ComplexRoots(std::ostream& out)
	{
		// --- Setup: e^z = z has NO real solutions - its roots are genuinely complex
		ComplexFunction f([](Complex z) { return std::exp(z) - z; });

		// --- Solve: Muller's method from a single complex starting point
		auto result = RootFinding::FindRootMuller(f, Complex(0.5, 1.0));

		out << "  [2.4] root of e^z - z: " << result.root
		    << " (converged: " << std::boolalpha << result.converged
		    << ", |f| = " << std::abs(result.function_value) << ")\n";
	}

	// 2.5 Roots in multiple dimensions (nonlinear system)
	void Cookbook_Recipe02_5_NonlinearSystem(std::ostream& out)
	{
		// --- Setup: intersection of the circle x^2+y^2 = 4 with the curve e^x + y = 1
		RootFinding::DynamicSystemFunction F = [](const Vector<Real>& p) {
			return Vector<Real>({ p[0] * p[0] + p[1] * p[1] - 4.0,
			                      std::exp(p[0]) + p[1] - 1.0 });
		};

		// --- Solve: damped Newton with automatic numerical Jacobian
		auto result = RootFinding::SolveNonlinearSystemNewton(F, Vector<Real>({ 1.0, -1.7 }));

		out << "  [2.5] intersection point: " << result.solution
		    << " (residual norm: " << result.residual_norm
		    << ", " << result.iterations_used << " iters)\n";
	}

	void Cookbook_Recipe02_RootFinding(std::ostream& out)
	{
		Cookbook_Recipe02_1_BracketThenSolve(out);
		Cookbook_Recipe02_2_AllRoots(out);
		Cookbook_Recipe02_3_PolynomialRoots(out);
		Cookbook_Recipe02_4_ComplexRoots(out);
		Cookbook_Recipe02_5_NonlinearSystem(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 3: Calculating integrals
	///////////////////////////////////////////////////////////////////////////

	// 3.1 Real function - trapezoid, Romberg, adaptive Gauss-Kronrod
	void Cookbook_Recipe03_1_RealFunction(std::ostream& out)
	{
		// --- Setup: damped oscillation on [0, 2*pi]; exact antiderivative known
		RealFunction f([](Real x) { return std::exp(-x / 2) * std::cos(3 * x); });
		const Real a = 0.0, b = 2 * Constants::PI;
		const Real exact = (Real(0.5) - Real(0.5) * std::exp(-Constants::PI)) / Real(9.25);

		// --- Solve: three quadratures, increasing sophistication
		Real trap    = IntegrateTrap(f, a, b);        // IntegrationResult converts to Real
		Real romberg = IntegrateRomberg(f, a, b);

		Real gkError = 0.0;
		Real gk = Integration::IntegrateGaussKronrod(f, a, b, &gkError);   // adaptive, error estimate

		out << "  [3.1] exact = " << exact
		    << ", trap = " << trap
		    << ", romberg = " << romberg
		    << ", gauss-kronrod = " << gk << " (err est " << gkError << ")\n";
	}

	// 3.2 2D integration - simple (curved domain) and adaptive (rectangle)
	void Cookbook_Recipe03_2_Integrate2D(std::ostream& out)
	{
		// --- Simple: integral of x^2+y^2 over the unit disk = pi/2 (y-limits as lambdas)
		ScalarFunction<2> f2([](const VectorN<Real, 2>& p) { return p[0] * p[0] + p[1] * p[1]; });

		Real disk = Integrate2D(f2, IntegrationMethod::GAUSS10,
			-1, 1,
			[](Real x) { return -std::sqrt(1 - x * x); },
			[](Real x) { return  std::sqrt(1 - x * x); });

		// --- Adaptive: sin(x)cos(y) on [0,pi] x [0,pi/2] = 2, with error control
		auto g2 = [](Real x, Real y) { return std::sin(x) * std::cos(y); };

		auto adaptive = Integration::IntegrateAdaptive2D(g2, 0.0, Constants::PI, 0.0, Constants::PI / 2, 1e-10);

		out << "  [3.2] disk integral = " << disk << " (exact pi/2 = " << Constants::PI / 2 << ")"
		    << ", adaptive = " << adaptive.value << " (exact 2, err est " << adaptive.error_estimate
		    << ", " << adaptive.function_evaluations << " evals)\n";
	}

	// 3.3 3D integration - simple (sphere volume) and adaptive (box)
	void Cookbook_Recipe03_3_Integrate3D(std::ostream& out)
	{
		// --- Simple: volume of the unit sphere = 4/3 pi, limits as nested lambdas
		ScalarFunction<3> one([](const VectorN<Real, 3>&) { return Real(1); });

		Real vol = Integrate3D(one,
			-1, 1,
			[](Real x) { return -std::sqrt(1 - x * x); },
			[](Real x) { return  std::sqrt(1 - x * x); },
			[](Real x, Real y) { return -std::sqrt(std::max(Real(0), 1 - x * x - y * y)); },
			[](Real x, Real y) { return  std::sqrt(std::max(Real(0), 1 - x * x - y * y)); });

		// --- Adaptive: x*y*z over [0,1]^3 = 1/8
		auto g3 = [](Real x, Real y, Real z) { return x * y * z; };

		auto adaptive = Integration::IntegrateAdaptive3D(g3, 0, 1, 0, 1, 0, 1, 1e-10);

		out << "  [3.3] sphere volume = " << vol << " (exact " << 4.0 / 3.0 * Constants::PI << ")"
		    << ", adaptive xyz = " << adaptive.value << " (exact 0.125)\n";
	}

	// 3.4 Improper integrals - infinite ranges and an endpoint singularity
	void Cookbook_Recipe03_4_ImproperIntegrals(std::ostream& out)
	{
		// Gaussian tail: integral 0..inf of e^(-x^2) = sqrt(pi)/2
		RealFunction gaussian([](Real x) { return std::exp(-x * x); });
		Real tail = IntegrateUpperInf(gaussian, 0.0);

		// Lorentzian over the whole real line: integral -inf..inf of 1/(1+x^2) = pi
		RealFunction lorentzian([](Real x) { return 1.0 / (1.0 + x * x); });
		Real whole = IntegrateInf(lorentzian);

		// Integrable singularity at the lower endpoint: integral 0..1 of 1/sqrt(x) = 2
		RealFunction invSqrt([](Real x) { return 1.0 / std::sqrt(x); });
		Real singular = IntegrateLowerSingular(invSqrt, 0.0, 1.0);

		out << "  [3.4] gaussian tail = " << tail << " (exact " << std::sqrt(Constants::PI) / 2 << ")"
		    << ", lorentzian = " << whole << " (exact pi)"
		    << ", 1/sqrt(x) = " << singular << " (exact 2)\n";
	}

	// 3.5 Path integration - work done by a force field along a helix
	void Cookbook_Recipe03_5_WorkIntegral(std::ostream& out)
	{
		// --- Setup: rotational + axial force F = (-y, x, z), one helix turn
		VectorFunction<3> force([](const VectorN<Real, 3>& p) {
			return VectorN<Real, 3>{ -p[1], p[0], p[2] };
		});
		ParametricCurve<3> helix([](Real t) { return VectorN<Real, 3>{ std::cos(t), std::sin(t), t }; });

		// --- Solve: W = integral F . dr; exact = 2*pi + 2*pi^2
		Real work = PathIntegration::LineIntegral(force, helix, 0.0, 2 * Constants::PI, 1e-8);

		out << "  [3.5] work along helix = " << work
		    << " (exact " << 2 * Constants::PI + 2 * Constants::PI * Constants::PI << ")\n";
	}

	// 3.6 Surface integration - flux of a vector field through a closed surface
	void Cookbook_Recipe03_6_FluxIntegral(std::ostream& out)
	{
		// --- Setup: radial field F = (x, y, z) through a cube of side 2 centered at origin
		VectorFunction<3> field([](const VectorN<Real, 3>& p) { return p; });
		Cube3D cube(2.0);

		// --- Solve: flux = closed-surface integral F . dS; div F = 3, so exact = 3 * V = 24
		Real flux = SurfaceIntegration::SurfaceIntegral(field, cube);

		out << "  [3.6] flux through cube = " << flux << " (exact 3 * volume = 24)\n";
	}

	// 3.7 Monte Carlo integration - a smooth 3D integral with statistical error
	void Cookbook_Recipe03_7_MonteCarlo(std::ostream& out)
	{
		// --- Setup: e^(-|p|^2) over [-1,1]^3; exact = (sqrt(pi) * erf(1))^3
		ScalarFunction<3> gauss3([](const VectorN<Real, 3>& p) {
			return std::exp(-(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]));
		});
		const Real exact = std::pow(std::sqrt(Constants::PI) * std::erf(1.0), 3);

		// --- Solve: seeded sampler, 500k samples, statistical error estimate
		MonteCarloIntegrator<3> mc(42);
		MonteCarloConfig config;
		config.num_samples = 500000;

		auto result = mc.integrate(gauss3,
			VectorN<Real, 3>{ -1, -1, -1 }, VectorN<Real, 3>{ 1, 1, 1 }, config);

		out << "  [3.7] monte carlo = " << result.value << " (exact " << exact
		    << ", std error " << result.error_estimate
		    << ", " << result.samples_used << " samples)\n";
	}

	void Cookbook_Recipe03_Integrals(std::ostream& out)
	{
		Cookbook_Recipe03_1_RealFunction(out);
		Cookbook_Recipe03_2_Integrate2D(out);
		Cookbook_Recipe03_3_Integrate3D(out);
		Cookbook_Recipe03_4_ImproperIntegrals(out);
		Cookbook_Recipe03_5_WorkIntegral(out);
		Cookbook_Recipe03_6_FluxIntegral(out);
		Cookbook_Recipe03_7_MonteCarlo(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 4: Calculating derivatives
	///////////////////////////////////////////////////////////////////////////

	// 4.1 Real function - first, second, and third derivatives
	void Cookbook_Recipe04_1_RealFunction(std::ostream& out)
	{
		// --- Setup: f(x) = sin(x), evaluated at x = pi/4
		RealFunction f([](Real x) { return std::sin(x); });
		const Real x = Constants::PI / 4.0;

		// --- Solve: fourth-order formulas for derivative orders one, two, and three
		Real first  = Derivation::NDer4(f, x);
		Real second = Derivation::NSecDer4(f, x);
		Real third  = Derivation::NThirdDer4(f, x);

		out << "  [4.1] sin(x) at pi/4: f' = " << first
		    << ", f'' = " << second << ", f''' = " << third << "\n";
	}

	// 4.2 Scalar function - one partial derivative and the full gradient
	void Cookbook_Recipe04_2_ScalarFunction(std::ostream& out)
	{
		// --- Setup: f(x,y,z) = x^2*y + sin(z), at p = (1,2,0)
		ScalarFunction<3> f([](const VectorN<Real, 3>& p) {
			return p[0] * p[0] * p[1] + std::sin(p[2]);
		});
		VectorN<Real, 3> point{ 1.0, 2.0, 0.0 };

		// --- Solve: one component and all components of the gradient
		Real df_dx = Derivation::NDer4Partial(f, 0, point);
		VectorN<Real, 3> gradient = Derivation::NDer4PartialByAll(f, point);

		out << "  [4.2] scalar partial df/dx = " << df_dx
		    << ", gradient = " << gradient << "\n";
	}

	// 4.3 Vector function - one partial, one output gradient, and the full Jacobian
	void Cookbook_Recipe04_3_VectorFunction(std::ostream& out)
	{
		// --- Setup: F(x,y,z) = (xy, yz, zx), at p = (1,2,3)
		VectorFunction<3> field([](const VectorN<Real, 3>& p) {
			return VectorN<Real, 3>{ p[0] * p[1], p[1] * p[2], p[2] * p[0] };
		});
		VectorN<Real, 3> point{ 1.0, 2.0, 3.0 };

		// --- Solve: dF0/dy, grad(F0), and all dFi/dxj
		Real partial = Derivation::NDer4Partial(field, 0, 1, point);
		VectorN<Real, 3> firstComponentGradient =
			Derivation::NDer4PartialByAll(field, 0, point);
		MatrixNM<Real, 3, 3> allByAll =
			Derivation::NDer4PartialAllByAll(field, point);

		out << "  [4.3] vector partial dF0/dy = " << partial
		    << ", grad(F0) = " << firstComponentGradient
		    << ", all-by-all = " << allByAll << "\n";
	}

	// 4.4 Parametric curve - first, second, and third derivatives
	void Cookbook_Recipe04_4_ParametricCurve(std::ostream& out)
	{
		// --- Setup: circular helix r(t) = (cos(t), sin(t), t)
		ParametricCurve<3> helix([](Real t) {
			return VectorN<Real, 3>{ std::cos(t), std::sin(t), t };
		});

		// --- Solve: velocity/tangent, acceleration, and jerk at t = 0
		VectorN<Real, 3> first  = Derivation::NDer4(helix, 0.0);
		VectorN<Real, 3> second = Derivation::NSecDer4(helix, 0.0);
		VectorN<Real, 3> third  = Derivation::NThirdDer4(helix, 0.0);

		out << "  [4.4] helix at t=0: r' = " << first
		    << ", r'' = " << second << ", r''' = " << third << "\n";
	}

	// 4.5 Parametric surface - tangent plane, unit normal, and area scale
	void Cookbook_Recipe04_5_ParametricSurface(std::ostream& out)
	{
		// --- Setup: sphere of radius 2, parameterized by polar angle u and azimuth w
		ParametricSurfaceRect<3> sphere([](Real u, Real w) {
			return VectorN<Real, 3>{
				2.0 * std::sin(u) * std::cos(w),
				2.0 * std::sin(u) * std::sin(w),
				2.0 * std::cos(u)
			};
		});
		const Real u = Constants::PI / 4.0;
		const Real w = Constants::PI / 3.0;

		// --- Solve: the two partials span the tangent plane
		VectorN<Real, 3> tangentU = Derivation::NDer2_u(sphere, u, w);
		VectorN<Real, 3> tangentW = Derivation::NDer2_w(sphere, u, w);
		VectorN<Real, 3> normalCross{
			tangentU[1] * tangentW[2] - tangentU[2] * tangentW[1],
			tangentU[2] * tangentW[0] - tangentU[0] * tangentW[2],
			tangentU[0] * tangentW[1] - tangentU[1] * tangentW[0]
		};
		Real areaScale = normalCross.NormL2();
		VectorN<Real, 3> unitNormal = normalCross / areaScale;

		out << "  [4.5] sphere tangents: r_u = " << tangentU << ", r_w = " << tangentW
		    << ", unit normal = " << unitNormal << ", area scale = " << areaScale << "\n";
	}

	// 4.6 Jacobian and Hessian convenience functions
	void Cookbook_Recipe04_6_JacobianAndHessian(std::ostream& out)
	{
		// --- Setup: a vector map and a scalar quadratic, both at (1,2)
		VectorFunction<2> map([](const VectorN<Real, 2>& p) {
			return VectorN<Real, 2>{ p[0] * p[0] + p[1], p[0] * p[1] };
		});
		ScalarFunction<2> quadratic([](const VectorN<Real, 2>& p) {
			return p[0] * p[0] + p[0] * p[1] + 3.0 * p[1] * p[1];
		});
		VectorN<Real, 2> point{ 1.0, 2.0 };

		// --- Solve: J(i,j) = dFi/dxj; H(i,j) = d2f/dxi dxj
		MatrixNM<Real, 2, 2> jacobian = Derivation::calcJacobian(map, point);
		MatrixNM<Real, 2, 2> hessian = Derivation::calcHessian(quadratic, point);

		out << "  [4.6] Jacobian = " << jacobian << ", Hessian = " << hessian << "\n";
	}

	// 4.7 Complex-step differentiation - high precision without cancellation
	void Cookbook_Recipe04_7_ComplexStep(std::ostream& out)
	{
		// --- Setup: one analytic formula valid for both Real and Complex arguments
		auto f = [](auto x) { return std::exp(x) * std::sin(x); };
		const Real x = 1.0;
		const Real exact = std::exp(x) * (std::sin(x) + std::cos(x));

		// --- Solve: subtraction loses the tiny real step; complex-step does not subtract
		Real tinyFiniteDifference = Derivation::NDer2(f, x, 1e-20);
		Real complexStep = Derivation::ComplexStep(f, x); // default h is extremely small

		out << "  [4.7] exact f'(1) = " << exact
		    << ", finite difference (h=1e-20) = " << tinyFiniteDifference
		    << ", complex-step = " << complexStep << "\n";
	}

	void Cookbook_Recipe04_Derivatives(std::ostream& out)
	{
		Cookbook_Recipe04_1_RealFunction(out);
		Cookbook_Recipe04_2_ScalarFunction(out);
		Cookbook_Recipe04_3_VectorFunction(out);
		Cookbook_Recipe04_4_ParametricCurve(out);
		Cookbook_Recipe04_5_ParametricSurface(out);
		Cookbook_Recipe04_6_JacobianAndHessian(out);
		Cookbook_Recipe04_7_ComplexStep(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 5: Optimization recipes
	///////////////////////////////////////////////////////////////////////////

	// 5.1 Real function - golden section search and Brent's method
	void Cookbook_Recipe05_1_RealFunctionGoldenAndBrent(std::ostream& out)
	{
		// --- Setup: smooth one-dimensional objective with a minimum near x = 1.275
		RealFunction objective([](Real x) {
			return (x - 1.25) * (x - 1.25) + 0.15 * std::sin(3.0 * x);
		});

		// --- Solve: bracket once, then compare a robust golden search with Brent
		MinimumBracket bracket = Minimization::BracketMinimum(objective, 0.0, 1.0);
		MinimizationResult golden = Minimization::GoldenSectionSearch(objective, bracket, 1e-10);
		MinimizationResult brent = Minimization::BrentMinimize(objective, bracket, 1e-10);

		out << "  [5.1] one-dimensional minimum: golden x = " << golden.xmin
		    << ", f = " << golden.fmin << " (" << golden.iterations << " iters); Brent x = "
		    << brent.xmin << ", f = " << brent.fmin << " (" << brent.iterations << " iters)\n";
	}

	// 5.2 Multidimensional function - Powell, Nelder-Mead, and BFGS
	class Cookbook_OffsetBowlWithGradient : public Optimization::IDifferentiableScalarFunction<2> {
	public:
		Real operator()(const VectorN<Real, 2>& x) const override
		{
			Real dx = x[0] - 1.0;
			Real dy = x[1] + 2.0;
			return dx * dx + 2.0 * dy * dy;
		}

		void Gradient(const VectorN<Real, 2>& x, VectorN<Real, 2>& grad) const override
		{
			grad[0] = 2.0 * (x[0] - 1.0);
			grad[1] = 4.0 * (x[1] + 2.0);
		}
	};

	void Cookbook_Recipe05_2_MultidimFunction(std::ostream& out)
	{
		// --- Setup: same quadratic bowl, with gradient available for BFGS
		Cookbook_OffsetBowlWithGradient objective;
		VectorN<Real, 2> start{ -3.0, 4.0 };

		// --- Solve: compare two derivative-free methods and one gradient-based method
		auto powell = Optimization::PowellMinimize<2>(objective, start);
		auto nelderMead = Optimization::NelderMeadMinimize<2>(objective, start, 0.75, 1e-10);
		auto bfgs = Optimization::BFGSMinimize<2>(objective, start);

		out << "  [5.2] multidim minimum from " << start
		    << ": Powell x = " << powell.xmin << ", f = " << powell.fmin
		    << "; Nelder-Mead x = " << nelderMead.xmin << ", f = " << nelderMead.fmin
		    << "; BFGS x = " << bfgs.xmin << ", f = " << bfgs.fmin << "\n";
	}

	// 5.3 Linear programming - production planning model
	void Cookbook_Recipe05_3_LinearProgramming(std::ostream& out)
	{
		// --- Setup: maximize profit 40*x + 30*y subject to resource limits
		Optimization::LinearProgram lp(2, "ProductionPlan");
		lp.SetVariableNames({ "standard", "premium" });
		lp.SetObjective({ 40.0, 30.0 }, Optimization::LPObjective::Maximize);
		lp.AddConstraint({ 2.0, 1.0 }, Optimization::LPConstraintType::LessEqual, 100.0, "machine-hours");
		lp.AddConstraint({ 1.0, 1.0 }, Optimization::LPConstraintType::LessEqual, 80.0, "assembly-hours");
		lp.AddConstraint({ 1.0, 0.0 }, Optimization::LPConstraintType::LessEqual, 40.0, "standard-demand");

		// --- Solve: dense full-tableau simplex; variables are nonnegative by default
		Optimization::LPResult result = Optimization::SolveLP(lp);

		out << "  [5.3] LP status: " << result.statusMessage();
		if (result.IsOptimal()) {
			out << ", standard = " << result.x[0]
			    << ", premium = " << result.x[1]
			    << ", profit = " << result.objectiveValue;
		}
		out << "\n";
	}

	void Cookbook_Recipe05_Optimization(std::ostream& out)
	{
		Cookbook_Recipe05_1_RealFunctionGoldenAndBrent(out);
		Cookbook_Recipe05_2_MultidimFunction(out);
		Cookbook_Recipe05_3_LinearProgramming(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 6: Solving differential equations
	///////////////////////////////////////////////////////////////////////////

	// 5.1 Fixed-step RK4 on a simple problem
	void Cookbook_Recipe06_1_RK4(std::ostream& out)
	{
		// --- Setup: harmonic oscillator x'' = -x as a first-order system (x, v)
		ODESystem sho(2, [](Real t, const Vector<Real>& y, Vector<Real>& dydt) {
			dydt[0] = y[1];
			dydt[1] = -y[0];
		});
		Vector<Real> y0({ 1.0, 0.0 });

		// --- Solve: 100 RK4 steps over one period; exact solution returns to (1, 0)
		RungeKutta4_StepCalculator rk4;
		ODESystemFixedStepSolver solver(sho, rk4);
		auto sol = solver.integrate(y0, 0.0, 2 * Constants::PI, 100);

		out << "  [6.1] RK4 after one period: x = " << sol.getXValue(100, 0)
		    << ", v = " << sol.getXValue(100, 1) << " (exact 1, 0)\n";
	}

	// 5.2 Adaptive Dormand-Prince on a hard problem (chaotic Lorenz system)
	void Cookbook_Recipe06_2_AdaptiveDormandPrince(std::ostream& out)
	{
		// --- Setup: Lorenz attractor, sigma=10, rho=28, beta=8/3
		ODESystem lorenz(3, [](Real t, const Vector<Real>& y, Vector<Real>& dydt) {
			dydt[0] = 10.0 * (y[1] - y[0]);
			dydt[1] = y[0] * (28.0 - y[2]) - y[1];
			dydt[2] = y[0] * y[1] - 8.0 / 3.0 * y[2];
		});
		Vector<Real> y0({ 1.0, 1.0, 1.0 });

		// --- Solve: DP5 with error control; step size adapts to the dynamics
		ODEAdaptiveIntegrator<DormandPrince5_Stepper> integrator(lorenz);
		auto sol = integrator.integrate(y0, 0.0, 10.0, 0.1, 1e-8);

		out << "  [6.2] Lorenz to t=10: " << sol.size() << " saved points, final = "
		    << sol.getXValuesAtEnd() << "\n";
	}

	// 5.3 Event detection - artillery shot, integration stops when the shell lands
	class ArtilleryShot : public IODESystemWithEvents {
		Real _g, _drag;
	public:
		ArtilleryShot(Real g = 9.81, Real drag = 0.001) : _g(g), _drag(drag) {}

		int getDim() const override { return 4; }   // (x, y, vx, vy)

		void derivs(Real t, const Vector<Real>& s, Vector<Real>& dsdt) const override {
			Real speed = std::sqrt(s[2] * s[2] + s[3] * s[3]);
			dsdt[0] = s[2];
			dsdt[1] = s[3];
			dsdt[2] = -_drag * speed * s[2];          // quadratic air drag
			dsdt[3] = -_g - _drag * speed * s[3];
		}

		int getNumEvents() const override { return 1; }
		Real eventFunction(int, Real, const Vector<Real>& s) const override { return s[1]; } // y = 0
		EventDirection getEventDirection(int) const override { return EventDirection::Decreasing; }
		EventAction getEventAction(int) const override { return EventAction::Stop; }
	};

	void Cookbook_Recipe06_3_EventDetection(std::ostream& out)
	{
		// --- Setup: shell fired at 45 degrees, 100 m/s, with quadratic drag
		ArtilleryShot shot;
		Real v0 = 100.0, angle = Constants::PI / 4;
		Vector<Real> s0({ 0.0, 0.0, v0 * std::cos(angle), v0 * std::sin(angle) });

		// --- Solve: integrate until the ground-hit event terminates the run
		DormandPrince5EventIntegrator integrator(shot);
		auto result = integrator.integrateWithEvents(shot, s0, 0.0, 60.0, 0.05, 1e-10, 1e-12);

		const auto& impact = result.events.back();
		out << "  [6.3] shell lands at t = " << impact.time
		    << " s, range = " << impact.state[0]
		    << " m (terminatedByEvent: " << std::boolalpha << result.terminatedByEvent
		    << "; vacuum range would be " << v0 * v0 / 9.81 << " m)\n";
	}

	// 5.6 Stiff system - semi-implicit Rosenbrock where explicit RK4 would need h < 0.002
	// Stiff solvers need the Jacobian df/dy, so the system implements IODESystemWithJacobian
	class StiffRelaxation : public IODESystemWithJacobian {
	public:
		// y' = -1000 (y - sin t) + cos t, y(0) = 0; exact solution y = sin t
		int getDim() const override { return 1; }

		void derivs(Real t, const Vector<Real>& y, Vector<Real>& dydt) const override {
			dydt[0] = -1000.0 * (y[0] - std::sin(t)) + std::cos(t);
		}

		void jacobian(const Real t, const Vector<Real>& y, Vector<Real>& dydt, Matrix<Real>& J) const override {
			J(0, 0) = -1000.0;
		}
	};

	void Cookbook_Recipe06_6_StiffSystem(std::ostream& out)
	{
		StiffRelaxation stiff;
		Vector<Real> y0({ 0.0 });

		// --- Solve: Rosenbrock 2(3) takes large stable steps through the stiffness
		// (2nd-order method: keep tolerances moderate, 1e-6 here)
		Rosenbrock23Solver solver(stiff, 1e-6, 1e-6);
		auto sol = solver.Solve(0.0, y0, 1.0, 0.1);

		out << "  [6.6] stiff y(1) = " << sol.getXValuesAtEnd()[0]
		    << " (exact sin 1 = " << std::sin(1.0) << "), "
		    << sol.getNumStepsOK() << " accepted / " << sol.getNumStepsBad() << " rejected\n";
	}

	// 5.4 Boundary value problem via the shooting method
	void Cookbook_Recipe06_4_ShootingMethod(std::ostream& out)
	{
		// --- Setup: y'' = -y with y(0) = 0 and y(pi/2) = 1; exact y = sin t, so y'(0) = 1
		ODESystem sho(2, [](Real t, const Vector<Real>& y, Vector<Real>& dydt) {
			dydt[0] = y[1];
			dydt[1] = -y[0];
		});

		// --- Solve: shoot for the unknown initial velocity from two guesses
		BVPShootingSolver solver(sho);
		auto result = solver.solvePositionBVP(0.0, Constants::PI / 2,
		                                      0.0, 1.0,     // y(0) = 0, y(pi/2) = 1
		                                      0.5, 2.0);    // initial velocity guesses

		out << "  [6.4] shooting converged: " << std::boolalpha << result.converged
		    << ", found y'(0) = " << result.shootingParams[0]
		    << " (exact 1), residual = " << result.residual
		    << ", " << result.iterations << " iterations\n";
	}

	// 5.5 Where to aim? - the original shooting problem, with a free boundary
	// Muzzle speed is fixed by the gun; find the elevation that lands the shot on target.
	void Cookbook_Recipe06_5_AimTheGun(std::ostream& out)
	{
		const Real v0 = 100.0, targetRange = 400.0;
		const Real deg = Constants::PI / 180.0;

		// range(theta): fire with the ArtilleryShot system from 6.3, read the impact event
		auto rangeFor = [&](Real theta) {
			ArtilleryShot shot;
			Vector<Real> s0({ 0.0, 0.0, v0 * std::cos(theta), v0 * std::sin(theta) });
			DormandPrince5EventIntegrator integrator(shot);
			auto res = integrator.integrateWithEvents(shot, s0, 0.0, 60.0, 0.05, 1e-10, 1e-12);
			return res.events.back().state[0];
		};

		// --- Solve: root-find the elevation on the miss distance range(theta) - target.
		// With drag there are TWO solutions: direct fire (low) and mortar lob (high).
		auto miss = [&](Real theta) { return rangeFor(theta) - targetRange; };

		RootFinding::RootFindingConfig config;
		config.tolerance = 1e-10;
		auto lowSol  = RootFinding::FindRootBrent(miss,  5 * deg, 45 * deg, config);
		auto highSol = RootFinding::FindRootBrent(miss, 45 * deg, 85 * deg, config);

		out << "  [6.5] to land at " << targetRange << " m with v0 = " << v0 << " m/s:"
		    << " direct fire at " << lowSol.root / deg << " deg (lands " << rangeFor(lowSol.root)
		    << " m), mortar lob at " << highSol.root / deg << " deg (lands "
		    << rangeFor(highSol.root) << " m)\n";
	}

	void Cookbook_Recipe06_DifferentialEquations(std::ostream& out)
	{
		Cookbook_Recipe06_1_RK4(out);
		Cookbook_Recipe06_2_AdaptiveDormandPrince(out);
		Cookbook_Recipe06_3_EventDetection(out);
		Cookbook_Recipe06_4_ShootingMethod(out);
		Cookbook_Recipe06_5_AimTheGun(out);
		Cookbook_Recipe06_6_StiffSystem(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 7: Eigenvalues and eigenvectors
	///////////////////////////////////////////////////////////////////////////

	// 7.1 Eigenvalues and eigenvectors of a real symmetric matrix
	void Cookbook_Recipe07_1_RealSymmetric(std::ostream& out)
	{
		// --- Setup: symmetric matrix with eigenvalues 3, 3, and 6
		Matrix<Real> A(3, 3, { 4.0, 1.0, 1.0,
		                       1.0, 4.0, 1.0,
		                       1.0, 1.0, 4.0 });

		// --- Solve: eigenvectors are the columns corresponding to sorted eigenvalues
		auto result = SymmMatEigenSolverJacobi::Solve(A);

		Real maxEigenpairResidual = 0.0;
		for (int column = 0; column < A.cols(); column++) {
			Vector<Real> eigenvector = result.eigenvectors.VectorFromColumn(column);
			Vector<Real> residual = A * eigenvector - result.eigenvalues[column] * eigenvector;
			maxEigenpairResidual = std::max(maxEigenpairResidual, residual.NormL2());
		}

		out << "  [7.1] symmetric eigenvalues: " << result.eigenvalues
		    << " (converged: " << std::boolalpha << result.converged
		    << ", max ||Av-lambda*v||: " << maxEigenpairResidual << ")\n";
	}

	// 7.2 Eigenvalues of a general real matrix
	void Cookbook_Recipe07_2_GeneralReal(std::ostream& out)
	{
		// --- Setup: planar rotation plus independent scaling; eigenvalues are +/-i and 2
		Matrix<Real> A(3, 3, { 0.0, -1.0, 0.0,
		                       1.0,  0.0, 0.0,
		                       0.0,  0.0, 2.0 });

		// --- Solve: a real nonsymmetric matrix can have complex conjugate eigenvalues
		auto result = EigenSolver::Solve(A);

		out << "  [7.2] general real eigenvalues:";
		for (const auto& eigenvalue : result.eigenvalues)
			out << " " << eigenvalue;
		out << " (converged: " << std::boolalpha << result.converged << ")\n";
	}

	// 7.3 Eigenvalues and eigenvectors of a complex Hermitian matrix
	void Cookbook_Recipe07_3_ComplexHermitian(std::ostream& out)
	{
		// --- Setup: A = A*, so all eigenvalues are real and eigenvectors are unitary
		Matrix<Complex> A(3, 3, {
			{2.0, 0.0}, {1.0, 1.0},  {0.0, -1.0},
			{1.0, -1.0}, {3.0, 0.0}, {2.0, 0.5},
			{0.0, 1.0}, {2.0, -0.5}, {5.0, 0.0}
		});

		// --- Solve: eigenvalues are real; eigenvector columns are complex
		auto result = HermitianMatEigenSolverJacobi::Solve(A);

		Real maxEigenpairResidual = 0.0;
		for (int column = 0; column < A.cols(); column++) {
			Vector<Complex> eigenvector = result.eigenvectors.VectorFromColumn(column);
			Vector<Complex> residual = A * eigenvector
				- Complex(result.eigenvalues[column], 0.0) * eigenvector;
			maxEigenpairResidual = std::max(maxEigenpairResidual, residual.NormL2());
		}

		out << "  [7.3] Hermitian eigenvalues: " << result.eigenvalues
		    << " (converged: " << std::boolalpha << result.converged
		    << ", max ||Av-lambda*v||: " << maxEigenpairResidual << ")\n";
	}

	// 7.4 Eigenvalues of a general complex matrix
	void Cookbook_Recipe07_4_GeneralComplex(std::ostream& out)
	{
		// --- Setup: triangular but non-Hermitian; eigenvalues are its diagonal entries
		Matrix<Complex> A(3, 3, {
			{1.0, 1.0}, {2.0, 0.0},  {0.0, 0.0},
			{0.0, 0.0}, {-2.0, 0.5}, {1.0, -1.0},
			{0.0, 0.0}, {0.0, 0.0},  {3.0, -2.0}
		});

		// --- Solve: shifted complex QR computes all eigenvalues and right eigenvectors
		auto result = ComplexEigenSolver::Solve(A);

		out << "  [7.4] general complex eigenvalues: " << result.eigenvalues
		    << " (converged: " << std::boolalpha << result.converged
		    << ", max residual: " << result.maxResidual << ")\n";
	}

	void Cookbook_Recipe07_EigenvaluesAndEigenvectors(std::ostream& out)
	{
		Cookbook_Recipe07_1_RealSymmetric(out);
		Cookbook_Recipe07_2_GeneralReal(out);
		Cookbook_Recipe07_3_ComplexHermitian(out);
		Cookbook_Recipe07_4_GeneralComplex(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 8: Matrix decompositions
	///////////////////////////////////////////////////////////////////////////

	Matrix<Real> Cookbook_MakeSigma(const MatrixAlg::SVDDecomposition<Real>& svd, int rows, int cols, int keptValues)
	{
		Matrix<Real> sigma(rows, cols);
		const int count = std::min({ keptValues, svd.singularValues.size(), rows, cols });
		for (int index = 0; index < count; ++index)
			sigma(index, index) = svd.singularValues[index];
		return sigma;
	}

	Matrix<Real> Cookbook_PermuteRows(const Matrix<Real>& matrix, const std::vector<int>& permutation)
	{
		Matrix<Real> result(matrix.rows(), matrix.cols());
		for (int row = 0; row < matrix.rows(); ++row)
			for (int col = 0; col < matrix.cols(); ++col)
				result(row, col) = matrix(permutation[row], col);
		return result;
	}

	// 8.1 SVD: singular values, rank, and condition number
	void Cookbook_Recipe08_1_SVDAnatomy(std::ostream& out)
	{
		Matrix<Real> A(4, 3, { 1.0,  0.0,  2.0,
		                       0.0,  1.0, -1.0,
		                       2.0,  1.0,  0.0,
		                       1.0, -1.0,  1.0 });

		auto svd = MatrixAlg::SVDDecompose(A);
		Matrix<Real> reconstructed = svd.U * Cookbook_MakeSigma(svd, A.rows(), A.cols(), svd.singularValues.size()) * svd.V.transpose();
		Real residual = MatrixAlg::FrobeniusNorm(A - reconstructed);

		out << "  [8.1] SVD singular values: " << svd.singularValues
		    << " rank = " << svd.rank
		    << ", cond2 = " << MatrixAlg::ConditionNumber(A)
		    << ", reconstruction residual = " << residual << "\n";
	}

	// 8.2 SVD: low-rank approximation
	void Cookbook_Recipe08_2_SVDLowRankApproximation(std::ostream& out)
	{
		Matrix<Real> measurements(5, 4, { 10.0,  9.8,  6.0,  5.9,
		                               8.0,  7.9,  4.8,  4.7,
		                               6.0,  6.1,  3.5,  3.6,
		                               4.0,  4.1,  2.5,  2.4,
		                               2.0,  2.1,  1.2,  1.1 });

		auto svd = MatrixAlg::SVDDecompose(measurements);
		Matrix<Real> rank1 = svd.U * Cookbook_MakeSigma(svd, measurements.rows(), measurements.cols(), 1) * svd.V.transpose();
		Matrix<Real> rank2 = svd.U * Cookbook_MakeSigma(svd, measurements.rows(), measurements.cols(), 2) * svd.V.transpose();

		out << "  [8.2] SVD low-rank errors: rank-1 = "
		    << MatrixAlg::FrobeniusNorm(measurements - rank1)
		    << ", rank-2 = " << MatrixAlg::FrobeniusNorm(measurements - rank2) << "\n";
	}

	// 8.3 SVD: null space and rank deficiency
	void Cookbook_Recipe08_3_SVDNullSpace(std::ostream& out)
	{
		Matrix<Real> A(3, 4, { 1.0, 2.0, 0.0,  3.0,
		                       2.0, 4.0, 0.0,  6.0,
		                       0.0, 1.0, 1.0, -1.0 });

		auto spaces = MatrixAlg::FundamentalSubspacesOf(A);
		Real nullResidual = MatrixAlg::FrobeniusNorm(A * spaces.nullSpace);

		out << "  [8.3] rank-deficient matrix: rank = " << spaces.rank
		    << ", nullity = " << spaces.nullSpace.cols()
		    << ", ||A*N||F = " << nullResidual << "\n";
	}

	// 8.4 QR: orthogonal basis and reconstruction
	void Cookbook_Recipe08_4_QRBasis(std::ostream& out)
	{
		Matrix<Real> A(4, 2, { 1.0, 2.0,
		                       2.0, 0.0,
		                       0.0, 3.0,
		                       1.0, 1.0 });

		auto qr = MatrixAlg::QRDecompose(A);
		Real reconstruction = MatrixAlg::FrobeniusNorm(A - qr.Q * qr.R);
		Real orthogonality = MatrixAlg::FrobeniusNorm(qr.Q.transpose() * qr.Q - Matrix<Real>::Identity(qr.Q.cols()));

		out << "  [8.4] QR: ||A-QR||F = " << reconstruction
		    << ", ||Q^T Q-I||F = " << orthogonality << "\n";
	}

	// 8.5 LU: factors and determinant
	void Cookbook_Recipe08_5_LUFactors(std::ostream& out)
	{
		Matrix<Real> A(3, 3, { 0.0, 2.0, 1.0,
		                       1.0, 1.0, 0.0,
		                       2.0, 0.0, 1.0 });

		auto lu = MatrixAlg::LUDecompose(A);
		Matrix<Real> permuted = Cookbook_PermuteRows(A, lu.permutation);
		Real residual = MatrixAlg::FrobeniusNorm(lu.L * lu.U - permuted);

		out << "  [8.5] LU: det = " << lu.determinant
		    << ", ||LU-PA||F = " << residual
		    << ", first pivot row = " << lu.permutation.front() << "\n";
	}

	// 8.6 Cholesky: positive-definite decomposition
	void Cookbook_Recipe08_6_CholeskySPD(std::ostream& out)
	{
		Matrix<Real> covariance(3, 3, { 4.0, 2.0, 0.6,
		                              2.0, 3.0, 0.5,
		                              0.6, 0.5, 1.5 });

		auto chol = MatrixAlg::CholeskyDecompose(covariance);
		Real residual = MatrixAlg::FrobeniusNorm(chol.L * chol.L.transpose() - covariance);

		out << "  [8.6] Cholesky: L(0,0) = " << chol.L(0, 0)
		    << ", ||LL^T-A||F = " << residual << "\n";
	}

	void Cookbook_Recipe08_MatrixDecompositions(std::ostream& out)
	{
		Cookbook_Recipe08_1_SVDAnatomy(out);
		Cookbook_Recipe08_2_SVDLowRankApproximation(out);
		Cookbook_Recipe08_3_SVDNullSpace(out);
		Cookbook_Recipe08_4_QRBasis(out);
		Cookbook_Recipe08_5_LUFactors(out);
		Cookbook_Recipe08_6_CholeskySPD(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 9: Matrix properties and diagnostics
	///////////////////////////////////////////////////////////////////////////

	const char* Cookbook_StabilityName(MatrixAlg::MatrixStability stability)
	{
		switch (stability) {
		case MatrixAlg::MatrixStability::WellConditioned: return "well-conditioned";
		case MatrixAlg::MatrixStability::ModeratelyConditioned: return "moderately conditioned";
		case MatrixAlg::MatrixStability::IllConditioned: return "ill-conditioned";
		case MatrixAlg::MatrixStability::Singular: return "singular";
		}
		return "unknown";
	}

	class Cookbook_TimeDependentMatrix : public IFunction<Matrix<Real>, Real>
	{
	public:
		Matrix<Real> operator()(Real t) const override
		{
			return Matrix<Real>(2, 2, { 1.0 + t, 0.25,
			                           0.25, 2.0 - 0.5 * t });
		}
	};

	// 9.1 Inspect matrix structure
	void Cookbook_Recipe09_1_Structure(std::ostream& out)
	{
		Matrix<Real> A(3, 3, { 4.0, 1.0, 0.0,
		                       1.0, 3.0, 0.0,
		                       0.0, 0.0, 2.0 });
		MatrixAnalyzer<Real> analyzer(A);

		out << "  [9.1] structure: square = " << std::boolalpha << analyzer.IsSquare()
		    << ", symmetric = " << analyzer.IsSymmetric()
		    << ", diagonal = " << analyzer.IsDiagonal()
		    << ", diagonally dominant = " << analyzer.IsDiagonallyDominant() << "\n";
	}

	// 9.2 Norms and approximation error
	void Cookbook_Recipe09_2_NormsAndError(std::ostream& out)
	{
		Matrix<Real> A(3, 3, { 3.0, -1.0, 0.0,
		                       2.0,  4.0, 1.0,
		                       0.0, -2.0, 5.0 });
		Matrix<Real> approximation(3, 3, { 3.0, -1.0, 0.0,
		                                 2.0,  4.0, 0.9,
		                                 0.0, -2.0, 5.0 });

		out << "  [9.2] norms: ||A||F = " << MatrixAlg::FrobeniusNorm(A)
		    << ", ||A||1 = " << MatrixAlg::OneNorm(A)
		    << ", ||A||inf = " << MatrixAlg::InfinityNorm(A)
		    << ", approximation error = " << MatrixAlg::FrobeniusNorm(A - approximation) << "\n";
	}

	// 9.3 Trace, determinant, and rank
	void Cookbook_Recipe09_3_TraceDeterminantRank(std::ostream& out)
	{
		Matrix<Real> A(3, 3, { 2.0, 1.0, 0.0,
		                       1.0, 3.0, 1.0,
		                       0.0, 1.0, 2.0 });

		out << "  [9.3] invariants: trace = " << MatrixAlg::Trace(A)
		    << ", det = " << MatrixAlg::Determinant(A)
		    << ", rank = " << MatrixAlg::Rank(A)
		    << ", Gaussian rank = " << MatrixAlg::RankGaussian(A) << "\n";
	}

	// 9.4 Condition number and sensitivity
	void Cookbook_Recipe09_4_Conditioning(std::ostream& out)
	{
		Matrix<Real> A(2, 2, { 1.0, 0.999,
		                       0.999, 0.998 });
		MatrixAnalyzer<Real> analyzer(A);
		auto digitsLost = analyzer.ExpectedDigitsLost();

		out << "  [9.4] conditioning: cond2 = " << analyzer.ConditionNumber()
		    << ", stability = " << Cookbook_StabilityName(analyzer.AssessStability())
		    << ", digits lost ~= " << (digitsLost ? std::to_string(*digitsLost) : std::string("all")) << "\n";
	}

	// 9.5 Inverse and identity verification
	void Cookbook_Recipe09_5_InverseVerification(std::ostream& out)
	{
		Matrix<Real> A(3, 3, { 4.0, 1.0, 0.0,
		                       1.0, 3.0, 1.0,
		                       0.0, 1.0, 2.0 });
		Matrix<Real> inverse = MatrixAlg::Inverse(A);
		Real residual = MatrixAlg::FrobeniusNorm(A * inverse - Matrix<Real>::Identity(A.rows()));

		out << "  [9.5] inverse check: ||A*A^-1-I||F = " << residual
		    << ", inverse(0,0) = " << inverse(0, 0) << "\n";
	}

	// 9.6 Matrix slicing and transformations
	void Cookbook_Recipe09_6_SlicingTransformations(std::ostream& out)
	{
		Matrix<Real> samples(4, 4, { 1.0,  2.0,  3.0,  4.0,
		                          5.0,  6.0,  7.0,  8.0,
		                          9.0, 10.0, 11.0, 12.0,
		                         13.0, 14.0, 15.0, 16.0 });
		Matrix<Real> center(samples, 1, 1, 2, 2);
		Matrix<Real> transformed = center.transpose() * center;

		out << "  [9.6] slice transform: center trace = " << MatrixAlg::Trace(center)
		    << ", transformed(0,0) = " << transformed(0, 0)
		    << ", transformed sparsity = " << MatrixAlg::Sparsity(transformed) << "\n";
	}

	// 9.7 Matrix-valued functions
	void Cookbook_Recipe09_7_MatrixValuedFunctions(std::ostream& out)
	{
		Cookbook_TimeDependentMatrix A;
		Matrix<Real> atQuarter = A(0.25);
		Matrix<Real> atHalf = A(0.5);
		Matrix<Real> rate = (atHalf - atQuarter) * (1.0 / 0.25);

		out << "  [9.7] matrix-valued function: trace A(0.5) = " << MatrixAlg::Trace(atHalf)
		    << ", det A(0.5) = " << MatrixAlg::Determinant(atHalf)
		    << ", finite-difference rate(0,0) = " << rate(0, 0) << "\n";
	}

	void Cookbook_Recipe09_MatrixPropertiesAndDiagnostics(std::ostream& out)
	{
		Cookbook_Recipe09_1_Structure(out);
		Cookbook_Recipe09_2_NormsAndError(out);
		Cookbook_Recipe09_3_TraceDeterminantRank(out);
		Cookbook_Recipe09_4_Conditioning(out);
		Cookbook_Recipe09_5_InverseVerification(out);
		Cookbook_Recipe09_6_SlicingTransformations(out);
		Cookbook_Recipe09_7_MatrixValuedFunctions(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 10: Statistics
	///////////////////////////////////////////////////////////////////////////

	std::filesystem::path CookbookRepoPath(const std::filesystem::path& relativePath)
	{
		return std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() / relativePath;
	}

	// 10.1 Quick descriptive summary
	void Cookbook_Recipe10_1_DescriptiveSummary(std::ostream& out)
	{
		Vector<Real> data({ 2, 4, 4, 4, 5, 5, 7, 9 });

		Real mean, sampleStdDev, q1, median, q3, minimum, maximum;
		Statistics::AvgStdDev(data, mean, sampleStdDev);
		Statistics::Quartiles(data, q1, median, q3);
		Statistics::MinMax(data, minimum, maximum);

		out << "  [10.1] mean = " << mean << ", median = " << median
		    << ", sample sd = " << sampleStdDev
		    << ", population sd = " << Statistics::PopulationStdDev(data)
		    << ", quartiles = (" << q1 << ", " << median << ", " << q3 << ")"
		    << ", range = [" << minimum << ", " << maximum << "]\n";
	}

	// 10.2 Robust statistics in the presence of an outlier
	void Cookbook_Recipe10_2_RobustStatistics(std::ostream& out)
	{
		Vector<Real> responseMs({ 10, 12, 11, 13, 12, 11, 14, 100 });

		out << "  [10.2] response times: mean = " << Statistics::Mean(responseMs)
		    << ", median = " << Statistics::Median(responseMs)
		    << ", sd = " << Statistics::StdDev(responseMs)
		    << ", trimmed mean = " << Statistics::TrimmedMean(responseMs, 12.5)
		    << ", MAD = " << Statistics::MAD(responseMs)
		    << ", IQR = " << Statistics::IQR(responseMs) << "\n";
	}

	// 10.3 Weighted measurements
	void Cookbook_Recipe10_3_WeightedMeasurements(std::ostream& out)
	{
		Vector<Real> grades({ 85, 90, 78, 92 });
		Vector<Real> credits({ 2, 3, 1, 4 });
		Vector<Real> studyHours({ 6, 8, 4, 10 });

		out << "  [10.3] weighted grade = " << Statistics::WeightedMean(grades, credits)
		    << ", weighted sd = " << Statistics::WeightedStdDev(grades, credits)
		    << ", weighted corr(hours, grade) = "
		    << Statistics::WeightedPearsonCorrelation(studyHours, grades, credits) << "\n";
	}

	// 10.4 Relationships between variables
	void Cookbook_Recipe10_4_Relationships(std::ostream& out)
	{
		Vector<Real> studyHours({ 2, 3, 4, 5, 6, 7, 8, 9 });
		Vector<Real> testScores({ 65, 70, 72, 80, 85, 88, 92, 95 });
		Matrix<Real> observations(8, 3);
		for (int row = 0; row < 8; ++row) {
			observations(row, 0) = studyHours[row];
			observations(row, 1) = testScores[row];
			observations(row, 2) = 10.0 - 0.5 * studyHours[row]; // hours of sleep
		}

		Real correlation = Statistics::PearsonCorrelation(studyHours, testScores);
		Matrix<Real> correlationMatrix = Statistics::CorrelationMatrix(observations);

		out << "  [10.4] covariance = " << Statistics::Covariance(studyHours, testScores)
		    << ", Pearson r = " << correlation << ", R^2 = " << Statistics::RSquared(studyHours, testScores)
		    << ", corr(score,sleep) = " << correlationMatrix(1, 2) << "\n";
	}

	// 10.5 Histogram and empirical cumulative distribution
	void Cookbook_Recipe10_5_DistributionShape(std::ostream& out)
	{
		Vector<Real> deliveryMinutes({ 18, 19, 20, 21, 21, 22, 23, 24, 24, 25,
		                               26, 27, 29, 31, 35, 48 });
		auto histogram = Statistics::Histogram::ComputeHistogramAuto(
			deliveryMinutes, Statistics::Histogram::BinningMethod::FreedmanDiaconis);
		auto cumulativeCounts = histogram.GetCumulativeCounts();
		Real fractionWithinThirty = Statistics::Histogram::EvaluateECDF(deliveryMinutes, 30.0);

		out << "  [10.5] delivery histogram: " << histogram.numBins << " bins, counts = "
		    << histogram.counts << ", cumulative = " << cumulativeCounts
		    << ", P(time <= 30 min) = " << fractionWithinThirty << "\n";
	}

	// 10.6 Confidence interval for an unknown mean
	void Cookbook_Recipe10_6_MeanConfidenceInterval(std::ostream& out)
	{
		Vector<Real> fillWeights({ 500.2, 499.8, 500.5, 501.0, 499.4, 500.1, 500.7, 499.9, 500.3, 499.6 });
		Real mean, sampleStdDev;
		Statistics::AvgStdDev(fillWeights, mean, sampleStdDev);
		Statistics::TDistribution studentT(fillWeights.size() - 1);
		Real critical = studentT.inverseCdf(0.975);
		Real margin = critical * sampleStdDev / std::sqrt(static_cast<Real>(fillWeights.size()));

		out << "  [10.6] mean fill = " << mean << " g, 95% CI = ["
		    << mean - margin << ", " << mean + margin << "]\n";
	}

	// 10.7 Real dataset: Palmer Penguins
	void Cookbook_Recipe10_7_PalmerPenguins(std::ostream& out)
	{
		auto path = CookbookRepoPath("test_data/statistics/palmer_penguins.csv");
		Data::Dataset penguins = Data::LoadCSV(path.string());
		const auto species = penguins.GetStringColumn("species");
		const Data::DataColumn& massColumn = penguins["body_mass_g"];
		const Data::DataColumn& flipperColumn = penguins["flipper_length_mm"];

		Vector<Real> bodyMass;
		Vector<Real> flipperLength;
		std::map<std::string, Vector<Real>> massBySpecies;
		for (std::size_t row = 0; row < penguins.NumRows(); ++row) {
			Real mass = massColumn.GetReal(row);
			Real flipper = flipperColumn.GetReal(row);
			if (!std::isfinite(mass) || !std::isfinite(flipper))
				continue;
			bodyMass.push_back(mass);
			flipperLength.push_back(flipper);
			massBySpecies[species[row]].push_back(mass);
		}

		Real meanMass, massStdDev, q1, median, q3;
		Statistics::AvgStdDev(bodyMass, meanMass, massStdDev);
		Statistics::Quartiles(bodyMass, q1, median, q3);
		auto histogram = Statistics::Histogram::ComputeHistogramAuto(
			bodyMass, Statistics::Histogram::BinningMethod::FreedmanDiaconis);
		Statistics::TDistribution studentT(bodyMass.size() - 1);
		Real margin = studentT.inverseCdf(0.975) * massStdDev
			/ std::sqrt(static_cast<Real>(bodyMass.size()));

		out << "  [10.7] Palmer Penguins: " << penguins.NumRows() << " rows, "
		    << bodyMass.size() << " complete mass/flipper pairs; mass mean/median = "
		    << meanMass << "/" << median << " g, IQR = " << q3 - q1
		    << " g, MAD = " << Statistics::MAD(bodyMass)
		    << ", 95% mean CI = [" << meanMass - margin << ", " << meanMass + margin << "]"
		    << ", histogram bins = " << histogram.numBins
		    << ", P(mass <= 5000 g) = " << Statistics::Histogram::EvaluateECDF(bodyMass, 5000.0)
		    << ", corr(flipper,mass) = " << Statistics::PearsonCorrelation(flipperLength, bodyMass)
		    << ", species means:";
		for (const auto& [name, masses] : massBySpecies)
			out << " " << name << "=" << Statistics::Mean(masses) << "g";
		out << "\n";
	}

	struct MonzaLap
	{
		std::string driver;
		std::string compound;
		Real raceLap;
		Real tyreAge;
		Real lapTime;
		int stint;
	};

	struct TwoFactorTrend
	{
		Real rawTyreAgeSecondsPerLap = 0.0;
		Real tyreAgeSecondsPerLap = 0.0;
		Real raceLapSecondsPerLap = 0.0;
		int observations = 0;
		bool identifiable = false;
	};

	TwoFactorTrend FitMonzaTrend(const std::vector<MonzaLap>& laps, const std::string& compound)
	{
		std::map<std::string, Vector<Real>> timesByStint;
		for (const auto& lap : laps)
			if (lap.compound == compound)
				timesByStint[lap.driver + ":" + std::to_string(lap.stint)].push_back(lap.lapTime);

		std::map<std::string, Real> medianByStint;
		std::map<std::string, Real> madByStint;
		for (const auto& [stint, times] : timesByStint) {
			medianByStint[stint] = Statistics::Median(times);
			madByStint[stint] = Statistics::MAD(times);
		}

		std::map<std::string, std::vector<const MonzaLap*>> cleanByDriver;
		for (const auto& lap : laps) {
			if (lap.compound != compound)
				continue;
			std::string stint = lap.driver + ":" + std::to_string(lap.stint);
			Real mad = madByStint[stint];
			Real cutoff = std::max(Real(1.0), 4.0 * mad);
			if (std::abs(lap.lapTime - medianByStint[stint]) <= cutoff)
				cleanByDriver[lap.driver].push_back(&lap);
		}

		Real ageAge = 0.0, ageRace = 0.0, raceRace = 0.0;
		Real ageTime = 0.0, raceTime = 0.0;
		int observations = 0;
		for (const auto& [driver, driverLaps] : cleanByDriver) {
			if (driverLaps.size() < 5)
				continue;
			Real meanAge = 0.0, meanRace = 0.0, meanTime = 0.0;
			for (const MonzaLap* lap : driverLaps) {
				meanAge += lap->tyreAge;
				meanRace += lap->raceLap;
				meanTime += lap->lapTime;
			}
			meanAge /= driverLaps.size();
			meanRace /= driverLaps.size();
			meanTime /= driverLaps.size();
			for (const MonzaLap* lap : driverLaps) {
				Real age = lap->tyreAge - meanAge;
				Real race = lap->raceLap - meanRace;
				Real time = lap->lapTime - meanTime;
				ageAge += age * age;
				ageRace += age * race;
				raceRace += race * race;
				ageTime += age * time;
				raceTime += race * time;
				++observations;
			}
		}

		Real determinant = ageAge * raceRace - ageRace * ageRace;
		Real rawSlope = ageAge > 0.0 ? ageTime / ageAge : 0.0;
		if (observations == 0 || std::abs(determinant) < 1e-12)
			return { rawSlope, 0.0, 0.0, observations, false };
		return {
			rawSlope,
			(ageTime * raceRace - raceTime * ageRace) / determinant,
			(raceTime * ageAge - ageTime * ageRace) / determinant,
			observations,
			true
		};
	}

	// 10.8 Real dataset: tyre degradation at the 2024 Italian Grand Prix
	void Cookbook_Recipe10_8_MonzaTyreDegradation(std::ostream& out)
	{
		auto path = CookbookRepoPath(
			"src/examples/03_formula_1_sim/data/2024/16_italian_grand_prix/lap_times.csv");
		Data::Dataset data = Data::LoadCSV(path.string());
		std::vector<MonzaLap> completeLaps;
		struct StintState { std::string compound; Real lastTyreAge = 0.0; int stint = 0; };
		std::map<std::string, StintState> stintByDriver;
		for (std::size_t row = 0; row < data.NumRows(); ++row) {
			Real lap = data["lap_number"].GetReal(row);
			Real tyreAge = data["tyre_life"].GetReal(row);
			Real lapTime = data["lap_time_s"].GetReal(row);
			Real sector1 = data["sector1_s"].GetReal(row);
			Real sector2 = data["sector2_s"].GetReal(row);
			Real sector3 = data["sector3_s"].GetReal(row);
			if (!std::isfinite(lap) || !std::isfinite(tyreAge) || !std::isfinite(lapTime)
				|| !std::isfinite(sector1) || !std::isfinite(sector2) || !std::isfinite(sector3))
				continue;
			std::string driver = data["driver"].GetAsString(row);
			std::string compound = data["compound"].GetAsString(row);
			auto& state = stintByDriver[driver];
			if (state.stint == 0 || state.compound != compound || tyreAge <= state.lastTyreAge)
				++state.stint;
			state.compound = compound;
			state.lastTyreAge = tyreAge;
			completeLaps.push_back({ driver, compound, lap, tyreAge, lapTime, state.stint });
		}

		TwoFactorTrend medium = FitMonzaTrend(completeLaps, "MEDIUM");
		TwoFactorTrend hard = FitMonzaTrend(completeLaps, "HARD");

		out << "  [10.8] Monza 2024 (driver-centered, robust-cleaned): MEDIUM raw trend "
		    << medium.rawTyreAgeSecondsPerLap << " s/lap (n=" << medium.observations
		    << "; tyre age and race lap are not separately identifiable); HARD adjusted tyre effect "
		    << hard.tyreAgeSecondsPerLap << " s/tyre-lap, race-lap effect "
		    << hard.raceLapSecondsPerLap << " s/lap (n=" << hard.observations
		    << ", identifiable=" << std::boolalpha << hard.identifiable << ")\n";
	}

	void Cookbook_Recipe10_Statistics(std::ostream& out)
	{
		Cookbook_Recipe10_1_DescriptiveSummary(out);
		Cookbook_Recipe10_2_RobustStatistics(out);
		Cookbook_Recipe10_3_WeightedMeasurements(out);
		Cookbook_Recipe10_4_Relationships(out);
		Cookbook_Recipe10_5_DistributionShape(out);
		Cookbook_Recipe10_6_MeanConfidenceInterval(out);
		Cookbook_Recipe10_7_PalmerPenguins(out);
		Cookbook_Recipe10_8_MonzaTyreDegradation(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 11: 3D geometry: planes, lines, and bodies
	///////////////////////////////////////////////////////////////////////////

	// 11.1 Build points, vectors, lines, and segments
	void Cookbook_Recipe11_1_PointsVectorsLinesSegments(std::ostream& out)
	{
		Pnt3Cart stationA(0.0, 0.0, 0.0);
		Pnt3Cart stationB(3.0, 4.0, 0.0);
		Vec3Cart baseline(stationA, stationB);
		Line3D surveyLine(stationA, baseline);
		SegmentLine3D surveySegment(stationA, stationB);

		Pnt3Cart pointAfterTwoUnits = surveyLine(2.0);

		out << "  [11.1] baseline length = " << baseline.NormL2()
		    << ", unit direction = " << baseline.GetAsUnitVector()
		    << ", segment length = " << surveySegment.Length()
		    << ", line(2) = " << Cookbook_FormatPoint(pointAfterTwoUnits) << "\n";
	}

	// 11.2 Intersect a line with a plane
	void Cookbook_Recipe11_2_LinePlaneIntersection(std::ostream& out)
	{
		Plane3D ground = Plane3D::GetXYPlane();
		Line3D descendingPath(Pnt3Cart(2.0, -1.0, 5.0), Vec3Cart(0.0, 1.0, -2.0));
		Line3D parallelPath(Pnt3Cart(0.0, 0.0, 3.0), Vec3Cart(1.0, 0.0, 0.0));

		Pnt3Cart hitPoint;
		bool hitsGround = ground.IntersectionWithLine(descendingPath, hitPoint);
		bool parallelHitsGround = ground.IntersectionWithLine(parallelPath, hitPoint);

		out << "  [11.2] descending path hits ground = " << std::boolalpha << hitsGround;
		if (hitsGround)
			out << " at " << Cookbook_FormatPoint(hitPoint);
		out << ", parallel path hits ground = " << parallelHitsGround << "\n";
	}

	// 11.3 Project points and measure distances
	void Cookbook_Recipe11_3_ProjectionAndDistances(std::ostream& out)
	{
		Plane3D platform(Pnt3Cart(0.0, 0.0, 1.0), Vec3Cart(0.0, 0.0, 1.0));
		Pnt3Cart sensor(2.0, 3.0, 5.0);
		Pnt3Cart footprint = platform.ProjectionToPlane(sensor);
		Real heightAbovePlatform = platform.DistToPoint(sensor);

		Line3D cableAxis(Pnt3Cart(0.0, 0.0, 0.0), Vec3Cart(1.0, 1.0, 0.0));
		Pnt3Cart nearestCablePoint = cableAxis.NearestPointOnLine(sensor);
		Real distanceToCable = cableAxis.Dist(sensor);

		out << "  [11.3] sensor footprint = " << Cookbook_FormatPoint(footprint)
		    << ", plane distance = " << heightAbovePlatform
		    << ", nearest cable point = " << Cookbook_FormatPoint(nearestCablePoint)
		    << ", cable distance = " << distanceToCable << "\n";
	}

	// 11.4 Work with triangles and rectangular surfaces
	void Cookbook_Recipe11_4_TrianglesAndSurfaces(std::ostream& out)
	{
		Triangle3D brace(Pnt3Cart(0.0, 0.0, 0.0), Pnt3Cart(4.0, 0.0, 0.0), Pnt3Cart(0.0, 3.0, 0.0));
		Pnt3Cart centroid = brace.Centroid();
		bool loadPointInside = brace.IsPointInside(Pnt3Cart(1.0, 1.0, 0.0));

		RectSurface3D panel(Pnt3Cart(0.0, 0.0, 0.0), Pnt3Cart(2.0, 0.0, 0.0),
		                    Pnt3Cart(2.0, 3.0, 0.0), Pnt3Cart(0.0, 3.0, 0.0));
		VectorN<Real, 3> panelCenter = panel(0.0, 0.0);

		out << "  [11.4] triangle area = " << brace.Area()
		    << ", centroid = " << Cookbook_FormatPoint(centroid)
		    << ", point inside = " << std::boolalpha << loadPointInside
		    << ", panel area = " << panel.getArea()
		    << ", panel(0,0) = " << panelCenter << "\n";
	}

	// 11.5 Create and query 3D bodies
	void Cookbook_Recipe11_5_Bodies(std::ostream& out)
	{
		Cube3D equipmentBox(2.0, Pnt3Cart(1.0, 2.0, 3.0));
		Sphere3D safetyZone(1.5, Pnt3Cart(1.0, 2.0, 3.0));
		Cylinder3D tank(1.0, 3.0, Pnt3Cart(0.0, 0.0, 0.0));

		Pnt3Cart probe(1.5, 2.0, 3.0);
		bool probeInBox = equipmentBox.IsInside(probe);
		bool probeInZone = safetyZone.IsInside(probe);

		out << "  [11.5] cube volume = " << equipmentBox.Volume()
		    << ", cube area = " << equipmentBox.SurfaceArea()
		    << ", probe in cube = " << std::boolalpha << probeInBox
		    << ", sphere volume = " << safetyZone.Volume()
		    << ", probe in sphere = " << probeInZone
		    << ", tank center = " << Cookbook_FormatPoint(tank.GetCenter()) << "\n";
	}

	// 11.6 Use bounding volumes for spatial checks
	void Cookbook_Recipe11_6_BoundingVolumes(std::ostream& out)
	{
		Cube3D equipmentBox(2.0, Pnt3Cart(1.0, 2.0, 3.0));
		Box3D bounds = equipmentBox.GetBoundingBox();
		BoundingSphere3D sphereBounds = equipmentBox.GetBoundingSphere();
		BoundingSphere3D nearbyObject(Pnt3Cart(3.0, 2.0, 3.0), 0.75);

		Pnt3Cart probe(1.8, 2.0, 3.0);
		Box3D expandedBounds = bounds.Expanded(0.5);
		bool containsProbe = bounds.Contains(probe);
		bool sphereOverlap = sphereBounds.Intersects(nearbyObject);

		out << "  [11.6] bounds size = " << bounds.Size()
		    << ", contains probe = " << std::boolalpha << containsProbe
		    << ", expanded volume = " << expandedBounds.Volume()
		    << ", bounding spheres overlap = " << sphereOverlap << "\n";
	}

	void Cookbook_Recipe11_Geometry3D(std::ostream& out)
	{
		Cookbook_Recipe11_1_PointsVectorsLinesSegments(out);
		Cookbook_Recipe11_2_LinePlaneIntersection(out);
		Cookbook_Recipe11_3_ProjectionAndDistances(out);
		Cookbook_Recipe11_4_TrianglesAndSurfaces(out);
		Cookbook_Recipe11_5_Bodies(out);
		Cookbook_Recipe11_6_BoundingVolumes(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 12: 3D geometry: algorithms
	///////////////////////////////////////////////////////////////////////////

	// 12.1 Classify two 3D lines
	void Cookbook_Recipe12_1_LineClassification(std::ostream& out)
	{
		Line3D craneRail(Pnt3Cart(0.0, 0.0, 0.0), Vec3Cart(1.0, 0.0, 0.0));
		Line3D walkway(Pnt3Cart(2.0, -1.0, 0.0), Vec3Cart(0.0, 1.0, 0.0));
		Line3D parallelRail(Pnt3Cart(0.0, 2.0, 0.0), Vec3Cart(1.0, 0.0, 0.0));
		Line3D sameRail(Pnt3Cart(3.0, 0.0, 0.0), Vec3Cart(-1.0, 0.0, 0.0));
		Line3D skewPipe(Pnt3Cart(0.0, 1.0, 1.0), Vec3Cart(0.0, 1.0, 0.0));

		auto crossing = craneRail.Intersection(walkway);
		auto parallel = craneRail.Intersection(parallelRail);
		auto coincident = craneRail.Intersection(sameRail);
		auto skew = craneRail.Intersection(skewPipe);

		out << "  [12.1] line classes: crossing=" << Cookbook_FormatLineIntersectionType(crossing.type)
		    << " at " << Cookbook_FormatPoint(crossing.Point())
		    << ", parallel=" << Cookbook_FormatLineIntersectionType(parallel.type)
		    << ", coincident=" << Cookbook_FormatLineIntersectionType(coincident.type)
		    << ", skew=" << Cookbook_FormatLineIntersectionType(skew.type) << "\n";
	}

	// 12.2 Find closest approach between skew lines
	void Cookbook_Recipe12_2_SkewClosestApproach(std::ostream& out)
	{
		Line3D cameraRay(Pnt3Cart(0.0, 0.0, 1.0), Vec3Cart(1.0, 0.0, 0.0));
		Line3D robotAxis(Pnt3Cart(0.0, 2.0, 0.0), Vec3Cart(0.0, 0.0, 1.0));
		auto approach = cameraRay.Intersection(robotAxis);

		out << "  [12.2] skew closest distance = " << approach.distance
		    << ", ray point = " << Cookbook_FormatPoint(approach.point1)
		    << ", axis point = " << Cookbook_FormatPoint(approach.point2) << "\n";
	}

	// 12.3 Project points onto lines, segments, and planes
	void Cookbook_Recipe12_3_ProjectAndMeasure(std::ostream& out)
	{
		Pnt3Cart sensor(2.0, 3.0, 5.0);
		Line3D cable(Pnt3Cart(0.0, 0.0, 1.0), Vec3Cart(1.0, 1.0, 0.0));
		SegmentLine3D boom(Pnt3Cart(0.0, 0.0, 0.0), Pnt3Cart(4.0, 0.0, 0.0));
		Plane3D deck(Pnt3Cart(0.0, 0.0, 2.0), Vec3Cart(0.0, 0.0, 1.0));

		Pnt3Cart cableFoot = cable.NearestPointOnLine(sensor);
		Pnt3Cart deckFoot = deck.ProjectionToPlane(sensor);

		out << "  [12.3] cable foot = " << Cookbook_FormatPoint(cableFoot)
		    << ", line distance = " << cable.Dist(sensor)
		    << ", segment distance = " << boom.Dist(sensor)
		    << ", deck foot = " << Cookbook_FormatPoint(deckFoot)
		    << ", plane distance = " << deck.DistToPoint(sensor) << "\n";
	}

	// 12.4 Intersect lines with planes and planes with planes
	void Cookbook_Recipe12_4_LinePlaneAndPlanePlane(std::ostream& out)
	{
		Plane3D floor = Plane3D::GetXYPlane();
		Plane3D wall(Pnt3Cart(2.0, 0.0, 0.0), Vec3Cart(1.0, 0.0, 0.0));
		Line3D drillPath(Pnt3Cart(2.0, 1.0, 3.0), Vec3Cart(0.0, 0.0, -1.0));

		Pnt3Cart floorHit;
		bool lineHitsFloor = floor.IntersectionWithLine(drillPath, floorHit);
		Line3D floorWallEdge;
		bool planesIntersect = floor.IntersectionWithPlane(wall, floorWallEdge);

		out << "  [12.4] line hits floor=" << std::boolalpha << lineHitsFloor
		    << " at " << Cookbook_FormatPoint(floorHit)
		    << ", floor/wall edge exists=" << planesIntersect
		    << ", edge direction=" << floorWallEdge.Direction() << "\n";
	}

	// 12.5 Work with triangle geometry and point containment
	void Cookbook_Recipe12_5_TriangleGeometry(std::ostream& out)
	{
		Triangle3D facet(Pnt3Cart(0.0, 0.0, 0.0), Pnt3Cart(4.0, 0.0, 0.0), Pnt3Cart(0.0, 3.0, 0.0));
		Pnt3Cart loadPoint(1.0, 1.0, 0.0);
		Plane3D supportPlane = facet.getDefinedPlane();

		out << "  [12.5] triangle area=" << facet.Area()
		    << ", centroid=" << Cookbook_FormatPoint(facet.Centroid())
		    << ", contains load=" << std::boolalpha << facet.IsPointInside(loadPoint)
		    << ", normal=" << supportPlane.Normal() << "\n";
	}

	// 12.6 Ray-pick a triangle
	void Cookbook_Recipe12_6_RayTrianglePicking(std::ostream& out)
	{
		Triangle3D facet(Pnt3Cart(0.0, 0.0, 0.0), Pnt3Cart(4.0, 0.0, 0.0), Pnt3Cart(0.0, 3.0, 0.0));
		Pnt3Cart rayOrigin(1.0, 1.0, 5.0);
		Vec3Cart rayDirection(0.0, 0.0, -1.0);
		auto hit = CompGeometry::Intersections::IntersectRayTriangle(rayOrigin, rayDirection, facet);
		auto [u, v, w] = hit.GetBarycentricCoords();

		out << "  [12.6] ray hit=" << std::boolalpha << hit.hit
		    << ", t=" << hit.t
		    << ", point=" << Cookbook_FormatPoint(hit.point)
		    << ", barycentric(w,u,v)=(" << w << "," << u << "," << v << ")\n";
	}

	// 12.7 Query a 3D point cloud efficiently
	void Cookbook_Recipe12_7_KDTreePointCloud(std::ostream& out)
	{
		std::vector<Pnt3Cart> landmarks = {
			Pnt3Cart(0.0, 0.0, 0.0), Pnt3Cart(2.0, 0.0, 0.0), Pnt3Cart(0.0, 2.0, 0.0),
			Pnt3Cart(0.0, 0.0, 2.0), Pnt3Cart(2.0, 2.0, 2.0), Pnt3Cart(5.0, 5.0, 0.0)
		};
		KDTree3D tree;
		tree.build(landmarks);

		Pnt3Cart query(1.1, 0.2, 0.1);
		auto nearest = tree.findNearest(query);
		auto neighbors = tree.findKNearest(query, 3);
		auto inRadius = tree.findInRadius(query, 2.25);
		auto inBox = tree.findInBox(Pnt3Cart(-0.1, -0.1, -0.1), Pnt3Cart(2.1, 2.1, 0.1));

		out << "  [12.7] nearest index=" << nearest.index
		    << ", distance=" << nearest.distance
		    << ", kNN count=" << neighbors.indices.size()
		    << ", radius count=" << inRadius.indices.size()
		    << ", box count=" << inBox.size() << "\n";
	}

	// 12.8 Build a convex hull and test containment
	void Cookbook_Recipe12_8_ConvexHull3D(std::ostream& out)
	{
		std::vector<Pnt3Cart> points = {
			Pnt3Cart(0.0, 0.0, 0.0), Pnt3Cart(1.0, 0.0, 0.0), Pnt3Cart(1.0, 1.0, 0.0), Pnt3Cart(0.0, 1.0, 0.0),
			Pnt3Cart(0.0, 0.0, 1.0), Pnt3Cart(1.0, 0.0, 1.0), Pnt3Cart(1.0, 1.0, 1.0), Pnt3Cart(0.0, 1.0, 1.0),
			Pnt3Cart(0.4, 0.4, 0.4)
		};
		auto hull = CompGeometry::ConvexHull3DComputer::Compute(points);

		out << "  [12.8] hull vertices=" << hull.NumVertices()
		    << ", faces=" << hull.NumFaces()
		    << ", volume=" << hull.Volume()
		    << ", area=" << hull.SurfaceArea()
		    << ", contains center=" << std::boolalpha << hull.Contains(Pnt3Cart(0.5, 0.5, 0.5))
		    << ", contains outside=" << hull.Contains(Pnt3Cart(2.0, 0.5, 0.5)) << "\n";
	}

	void Cookbook_Recipe12_Geometry3DAlgorithms(std::ostream& out)
	{
		Cookbook_Recipe12_1_LineClassification(out);
		Cookbook_Recipe12_2_SkewClosestApproach(out);
		Cookbook_Recipe12_3_ProjectAndMeasure(out);
		Cookbook_Recipe12_4_LinePlaneAndPlanePlane(out);
		Cookbook_Recipe12_5_TriangleGeometry(out);
		Cookbook_Recipe12_6_RayTrianglePicking(out);
		Cookbook_Recipe12_7_KDTreePointCloud(out);
		Cookbook_Recipe12_8_ConvexHull3D(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 13: Interpolating functions
	///////////////////////////////////////////////////////////////////////////

	// 13.1 Real function interpolation with several methods
	void Cookbook_Recipe13_1_RealFunctionMethods(std::ostream& out)
	{
		Vector<Real> x({ 0.0, 0.35, 0.8, 1.2, 1.65, 2.1, 2.55, 3.0, 3.45, 3.9, 4.4 });
		Vector<Real> y({ 0.02, 0.28, 0.61, 0.86, 1.01, 0.92, 0.63, 0.21, -0.18, -0.52, -0.77 });

		LinearInterpRealFunc linear(x, y);
		PolynomInterpRealFunc polynomial(x, y, 5);
		SplineInterpRealFunc spline(x, y);
		BarycentricRationalInterp rational(x, y, 3);

		const Real query = 2.35;
		const Real linearValue = linear(query);
		const Real polynomialValue = polynomial(query);
		const Real polynomialErrorEstimate = polynomial.getLastErrorEst();
		const Real splineValue = spline(query);
		const Real rationalValue = rational(query);

		out << "  [13.1] tabulated calibration at x = " << query
		    << ": linear = " << linearValue
		    << ", polynomial = " << polynomialValue << " (est " << polynomialErrorEstimate << ")"
		    << ", spline = " << splineValue << " (f' " << spline.Derivative(query) << ")"
		    << ", barycentric rational = " << rationalValue << "\n";
	}

	// 13.2 2D grid interpolation with bilinear and bicubic spline methods
	void Cookbook_Recipe13_2_Grid2D(std::ostream& out)
	{
		Vector<Real> gridX({ 0.0, 1.0, 2.0, 3.0, 4.0 });
		Vector<Real> gridY({ 0.0, 0.8, 1.6, 2.4 });
		Matrix<Real> z(5, 4, {
			10.0, 10.5, 11.0, 11.4,
			10.8, 11.4, 11.8, 12.1,
			11.3, 12.0, 12.2, 12.4,
			11.1, 11.7, 12.0, 12.2,
			10.6, 11.1, 11.5, 11.9
		});

		BilinearInterp2D bilinear(gridX, gridY, z);
		BicubicSplineInterp2D spline2D(gridX, gridY, z);

		const Real qx = 1.7, qy = 1.1;
		const Real bilinearValue = bilinear(qx, qy);
		const Real splineValue = spline2D(qx, qy);
		Real zValue = 0.0, dzdx = 0.0, dzdy = 0.0;
		bilinear.interpWithDerivatives(qx, qy, zValue, dzdx, dzdy);

		out << "  [13.2] tabulated 2D surface at (" << qx << ", " << qy << ")"
		    << ": bilinear = " << bilinearValue
		    << " (dz/dx approx " << dzdx << ", dz/dy approx " << dzdy << ")"
		    << ", bicubic spline = " << splineValue << "\n";
	}

	// 13.3 Parametric curve interpolation from waypoints
	void Cookbook_Recipe13_3_ParametricCurve(std::ostream& out)
	{
		Matrix<Real> points(5, 2, {
			0.0, 0.0,
			1.0, 0.2,
			1.8, 1.0,
			2.8, 0.8,
			3.5, 1.6
		});

		LinInterpParametricCurve<2> polyline(points);
		SplineInterpParametricCurve<2> smooth(points);

		VectorN<Real, 2> polylineMid = polyline(0.5);
		VectorN<Real, 2> smoothMid = smooth(0.5);
		VectorN<Real, 2> smoothQuarter = smooth(0.25);

		out << "  [13.3] waypoint curve: polyline(0.5) = " << polylineMid
		    << ", spline(0.25) = " << smoothQuarter
		    << ", spline(0.5) = " << smoothMid << "\n";
	}

	void Cookbook_Recipe13_InterpolatingFunctions(std::ostream& out)
	{
		Cookbook_Recipe13_1_RealFunctionMethods(out);
		Cookbook_Recipe13_2_Grid2D(out);
		Cookbook_Recipe13_3_ParametricCurve(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 14: Curve fitting
	///////////////////////////////////////////////////////////////////////////

	// 14.1 Sensor calibration line with ordinary least squares
	void Cookbook_Recipe14_1_SensorCalibration(std::ostream& out)
	{
		Vector<Real> voltage({ 0.12, 0.38, 0.64, 0.91, 1.17, 1.44, 1.70, 1.96, 2.23, 2.49 });
		Vector<Real> temperatureC({ 3.1, 9.8, 16.4, 23.2, 29.9, 36.7, 43.1, 49.8, 56.6, 63.0 });

		auto fit = LinearLeastSquaresDetailed(voltage, temperatureC);
		Real predicted = fit.a * 1.85 + fit.b;

		out << "  [14.1] thermistor calibration: T = " << fit.a << "*V + " << fit.b
		    << ", R^2 = " << fit.r_squared
		    << ", T(1.85 V) = " << predicted << " C\n";
	}

	// 14.2 Weighted least squares when measurements have different uncertainty
	void Cookbook_Recipe14_2_WeightedCalibration(std::ostream& out)
	{
		Vector<Real> referenceKg({ 0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0 });
		Vector<Real> sensorMv({ 0.03, 1.02, 2.03, 3.05, 4.10, 5.02, 6.38, 6.95 });
		Vector<Real> weights({ 4.0, 4.0, 4.0, 4.0, 3.0, 3.0, 0.5, 2.0 });
		Vector<std::function<Real(Real)>> basis({
			[](Real) { return 1.0; },
			[](Real x) { return x; }
		});

		auto fit = WeightedGeneralLinearLeastSquares(referenceKg, sensorMv, weights, basis);
		Real predicted = fit.evaluate(9.0, basis);

		out << "  [14.2] weighted load-cell fit: mV = " << fit.coefficients[1]
		    << "*kg + " << fit.coefficients[0]
		    << ", weighted R^2 = " << fit.r_squared
		    << ", mV at 9 kg = " << predicted << "\n";
	}

	// 14.3 Polynomial performance curve
	void Cookbook_Recipe14_3_PumpPerformanceCurve(std::ostream& out)
	{
		Vector<Real> flowLps({ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 });
		Vector<Real> headM({ 39.9, 39.2, 38.1, 36.5, 34.7, 32.3, 29.5, 26.2, 22.7 });

		auto fit = PolynomialFit(flowLps, headM, 2);
		Real headAt55 = EvaluatePolynomial(5.5, fit.coefficients);

		out << "  [14.3] pump curve coefficients = " << fit.coefficients
		    << ", adjusted R^2 = " << fit.adjusted_r_squared
		    << ", head at 5.5 L/s = " << headAt55 << " m\n";
	}

	// 14.4 Seasonal model with Fourier basis functions
	void Cookbook_Recipe14_4_SeasonalDemand(std::ostream& out)
	{
		Vector<Real> monthAngle(12);
		for (int month = 0; month < 12; ++month)
			monthAngle[month] = 2.0 * Constants::PI * month / 12.0;
		Vector<Real> energyMWh({ 42.1, 39.8, 36.5, 32.4, 29.6, 28.3, 30.1, 33.7, 37.9, 41.8, 45.6, 46.9 });
		Vector<std::function<Real(Real)>> basis({
			[](Real) { return 1.0; },
			[](Real t) { return std::sin(t); },
			[](Real t) { return std::cos(t); },
			[](Real t) { return std::sin(2.0 * t); },
			[](Real t) { return std::cos(2.0 * t); }
		});

		auto fit = GeneralLinearLeastSquares(monthAngle, energyMWh, basis);
		Real seasonalAmplitude = std::sqrt(fit.coefficients[1] * fit.coefficients[1]
		                                + fit.coefficients[2] * fit.coefficients[2]);
		Real julyForecast = fit.evaluate(Constants::PI, basis);

		out << "  [14.4] seasonal demand mean = " << fit.coefficients[0]
		    << " MWh, annual amplitude = " << seasonalAmplitude
		    << ", July forecast = " << julyForecast
		    << ", R^2 = " << fit.r_squared << "\n";
	}

	// 14.5 Exponential response curve with fixed decay bases
	void Cookbook_Recipe14_5_ExponentialResponse(std::ostream& out)
	{
		Vector<Real> minutes({ 0.0, 2.0, 4.0, 7.0, 10.0, 14.0, 18.0, 24.0, 30.0, 40.0 });
		Vector<Real> temperatureC({ 92.0, 78.6, 68.8, 57.1, 49.6, 42.9, 38.5, 34.4, 31.9, 29.6 });
		Vector<std::function<Real(Real)>> basis({
			[](Real) { return 1.0; },
			[](Real t) { return std::exp(-t / 5.0); },
			[](Real t) { return std::exp(-t / 20.0); }
		});

		auto fit = GeneralLinearLeastSquares(minutes, temperatureC, basis);
		Real temperatureAt12 = fit.evaluate(12.0, basis);

		out << "  [14.5] cooling model coefficients = " << fit.coefficients
		    << ", T(12 min) = " << temperatureAt12
		    << " C, R^2 = " << fit.r_squared << "\n";
	}

	// 14.6 Fit a smooth 2D path from measured points
	void Cookbook_Recipe14_6_PathFromPoints(std::ostream& out)
	{
		Vector<Real> x({ 0.0, 1.2, 2.5, 3.9, 5.1, 6.4, 7.2, 8.0 });
		Vector<Real> y({ 0.0, 0.4, 1.1, 1.5, 1.2, 0.6, -0.1, -0.4 });
		Vector<Real> t(x.size());
		t[0] = -1.0;
		Real totalLength = 0.0;
		for (int i = 1; i < x.size(); ++i) {
			Real dx = x[i] - x[i - 1];
			Real dy = y[i] - y[i - 1];
			totalLength += std::sqrt(dx * dx + dy * dy);
			t[i] = totalLength;
		}
		for (int i = 1; i < t.size(); ++i)
			t[i] = -1.0 + 2.0 * t[i] / totalLength;

		auto basis = MakeLegendreFitBasis(3);
		auto fitX = RidgeGeneralLinearLeastSquares(t, x, basis, Real{1e-4});
		auto fitY = RidgeGeneralLinearLeastSquares(t, y, basis, Real{1e-4});

		Real midX = fitX.evaluate(0.0, basis);
		Real midY = fitY.evaluate(0.0, basis);

		out << "  [14.6] fitted path midpoint = (" << midX << ", " << midY << ")"
		    << ", x R^2 = " << fitX.r_squared
		    << ", y R^2 = " << fitY.r_squared
		    << ", path length from data = " << totalLength << "\n";
	}

	void Cookbook_Recipe14_CurveFitting(std::ostream& out)
	{
		Cookbook_Recipe14_1_SensorCalibration(out);
		Cookbook_Recipe14_2_WeightedCalibration(out);
		Cookbook_Recipe14_3_PumpPerformanceCurve(out);
		Cookbook_Recipe14_4_SeasonalDemand(out);
		Cookbook_Recipe14_5_ExponentialResponse(out);
		Cookbook_Recipe14_6_PathFromPoints(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 15: Curves and surfaces
	///////////////////////////////////////////////////////////////////////////

	// 15.1 Use built-in curves and surfaces
	void Cookbook_Recipe15_1_BuiltIns(std::ostream& out)
	{
		Curves::Circle2DCurve circle(2.0, Pnt2Cart(1.0, -1.0));
		Curves::HelixCurve helix(1.5, 0.4);
		Surfaces::Sphere sphere(3.0);
		Surfaces::Torus torus(3.0, 0.75);

		VectorN<Real, 2> circlePoint = circle(Constants::PI / 3.0);
		VectorN<Real, 3> helixPoint = helix(Constants::PI);
		VectorN<Real, 3> spherePoint = sphere(Constants::PI / 2.0, Constants::PI / 4.0);
		VectorN<Real, 3> torusPoint = torus(0.0, Constants::PI);

		out << "  [15.1] built-ins: circle(pi/3)=" << circlePoint
		    << ", helix(pi)=" << helixPoint
		    << ", sphere equator=" << spherePoint
		    << ", torus inner=" << torusPoint << "\n";
	}

	// 15.2 Define custom curves and surfaces
	void Cookbook_Recipe15_2_CustomDefinitions(std::ostream& out)
	{
		Curves::CurveCartesian3D trajectory(0.0, 4.0, Cookbook_CustomCurvePoint);
		Cookbook_WavySheet sheet;

		VectorN<Real, 3> point = trajectory(2.0);
		VectorN<Real, 3> tangent = trajectory.getTangent(2.0);
		VectorN<Real, 3> sheetPoint = sheet(1.0, 0.5);
		VectorN<Real, 3> sheetNormal = sheet.Normal(1.0, 0.5);

		out << "  [15.2] custom curve p(2)=" << point
		    << ", tangent=" << tangent
		    << ", sheet(1,0.5)=" << sheetPoint
		    << ", normal=" << sheetNormal << "\n";
	}

	// 15.3 Calculate basic curve properties
	void Cookbook_Recipe15_3_CurveProperties(std::ostream& out)
	{
		Curves::Circle2DCurve circle(2.0);
		Curves::HelixCurve helix(1.5, 0.4);
		Real circleCurvature = circle.getCurvature(Constants::PI / 4.0);
		Real helixLength = PathIntegration::ParametricCurveLength<3>(helix, 0.0, 2.0 * Constants::PI);
		Real helixCurvature = helix.getCurvature(1.0);
		Real helixTorsion = helix.getTorsion(1.0);

		out << "  [15.3] curve properties: circle curvature=" << circleCurvature
		    << ", helix one-turn length=" << helixLength
		    << ", helix curvature/torsion=" << helixCurvature << "/" << helixTorsion << "\n";
	}

	// 15.4 Compute a Frenet frame
	void Cookbook_Recipe15_4_FrenetFrame(std::ostream& out)
	{
		Curves::HelixCurve helix(1.5, 0.4);
		Real t = 1.0;
		Vector3Cartesian tangent, normal, binormal;
		helix.getMovingTrihedron(t, tangent, normal, binormal);
		auto frame = DifferentialGeometry::ComputeFrenetFrame<3>(helix, t);

		out << "  [15.4] Frenet frame at t=1: T=" << tangent
		    << ", N=" << normal
		    << ", B=" << binormal
		    << ", speed=" << frame.speed
		    << ", curvature=" << frame.curvature << "\n";
	}

	// 15.5 Inspect surface normals and curvatures
	void Cookbook_Recipe15_5_SurfaceProperties(std::ostream& out)
	{
		Surfaces::Sphere sphere(2.0);
		Surfaces::CylinderSurface cylinder(1.5, 4.0);
		Surfaces::Helicoid helicoid(0.4);
		Real sphereU = Constants::PI / 3.0;
		Real sphereW = Constants::PI / 5.0;
		Real k1, k2;
		sphere.PrincipalCurvatures(sphereU, sphereW, k1, k2);

		out << "  [15.5] surfaces: sphere normal=" << sphere.Normal(sphereU, sphereW)
		    << ", K/H=" << sphere.GaussianCurvature(sphereU, sphereW) << "/" << sphere.MeanCurvature(sphereU, sphereW)
		    << ", principal=" << k1 << "/" << k2
		    << ", cylinder K=" << cylinder.GaussianCurvature(0.7, 2.0)
		    << ", helicoid H=" << helicoid.MeanCurvature(1.0, 0.7) << "\n";
	}

	// 15.6 Calculate first and second fundamental forms
	void Cookbook_Recipe15_6_FundamentalForms(std::ostream& out)
	{
		Surfaces::Torus torus(3.0, 0.75);
		Real u = 0.8;
		Real w = 1.1;
		Real E, F, G;
		Real L, M, N;
		torus.GetFirstNormalFormCoefficients(u, w, E, F, G);
		torus.GetSecondNormalFormCoefficients(u, w, L, M, N);

		out << "  [15.6] torus forms: I=(" << E << ", " << F << ", " << G << ")"
		    << ", II=(" << L << ", " << M << ", " << N << ")"
		    << ", K=" << torus.GaussianCurvature(u, w)
		    << ", H=" << torus.MeanCurvature(u, w) << "\n";
	}

	void Cookbook_Recipe15_CurvesAndSurfaces(std::ostream& out)
	{
		Cookbook_Recipe15_1_BuiltIns(out);
		Cookbook_Recipe15_2_CustomDefinitions(out);
		Cookbook_Recipe15_3_CurveProperties(out);
		Cookbook_Recipe15_4_FrenetFrame(out);
		Cookbook_Recipe15_5_SurfaceProperties(out);
		Cookbook_Recipe15_6_FundamentalForms(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 16: Polynoms
	///////////////////////////////////////////////////////////////////////////

	// 16.1 Build and evaluate a polynomial model with MML's coefficient order
	void Cookbook_Recipe16_1_CoefficientsAndEvaluation(std::ostream& out)
	{
		Polynom<Real> calibration({ 1.20, -0.35, 0.08, -0.004 });
		const Real input = 8.0;
		Real corrected = calibration(input);

		out << "  [16.1] p(x)=a0+a1*x+... with a0=" << calibration[0]
		    << ", degree=" << calibration.degree()
		    << ", leading=" << calibration.leadingTerm()
		    << ", p(" << input << ")=" << corrected << "\n";
	}

	// 16.2 Compose polynomial arithmetic into a power model
	void Cookbook_Recipe16_2_ArithmeticModel(std::ostream& out)
	{
		Polynom<Real> voltage({ 12.0, -0.30 });
		Polynom<Real> current({ 1.50, 0.20, -0.01 });
		Polynom<Real> cableLoss({ 0.0, 0.0, 0.35 });
		Polynom<Real> netPower = voltage * current - cableLoss;
		const Real load = 4.0;

		out << "  [16.2] voltage degree=" << voltage.degree()
		    << ", current degree=" << current.degree()
		    << ", net power degree=" << netPower.degree()
		    << ", netPower(" << load << ")=" << netPower(load) << "\n";
	}

	// 16.3 Divide out a known factor and inspect the residual
	void Cookbook_Recipe16_3_DivisionAndResidual(std::ostream& out)
	{
		Polynom<Real> knownFactor({ -2.0, 1.0 });
		Polynom<Real> processTrend({ 3.0, 0.5, -0.1 });
		Polynom<Real> measured = knownFactor * processTrend + Polynom<Real>::Constant(0.25);
		Polynom<Real> quotient, remainder;
		Polynom<Real>::poldiv(measured, knownFactor, quotient, remainder);

		out << "  [16.3] measured degree=" << measured.degree()
		    << ", quotient degree=" << quotient.degree()
		    << ", remainder=" << remainder.constantTerm()
		    << ", reconstruction at x=5 -> " << (knownFactor * quotient + remainder)(5.0) << "\n";
	}

	// 16.4 Differentiate and integrate polynomial models exactly
	void Cookbook_Recipe16_4_Calculus(std::ostream& out)
	{
		Polynom<Real> height({ 120.0, 35.0, -4.9, 0.15 });
		Polynom<Real> velocity = height.derivative();
		Polynom<Real> acceleration = velocity.derivative();
		Polynom<Real> accumulatedHeight = height.integral();
		Vector<Real> jet(4);
		const Real t = 3.0;
		height.Derive(t, jet);
		Real area = accumulatedHeight(4.0) - accumulatedHeight(0.0);

		out << "  [16.4] h(" << t << ")=" << jet[0]
		    << ", v=" << velocity(t)
		    << ", a=" << acceleration(t)
		    << ", Derive-pack p'/p''=" << jet[1] << "/" << jet[2]
		    << ", integral[0,4]=" << area << "\n";
	}

	// 16.5 Construct an explicit polynomial from measurement points
	void Cookbook_Recipe16_5_FromValues(std::ostream& out)
	{
		std::vector<Real> dose({ 0.0, 1.0, 2.0, 3.0, 4.0 });
		std::vector<Real> response({ 2.0, 3.4, 4.0, 4.1, 4.4 });
		Polynom<Real> responseCurve = Polynom<Real>::FromValues(dose, response);
		Polynom<Real> sensitivity = responseCurve.derivative();

		out << "  [16.5] interpolating polynomial degree=" << responseCurve.degree()
		    << ", p(2.5)=" << responseCurve(2.5)
		    << ", sensitivity p'(2.5)=" << sensitivity(2.5)
		    << ", p(4)=" << responseCurve(4.0) << "\n";
	}

	// 16.6 Use Chebyshev approximation and convert to a power-basis polynomial
	void Cookbook_Recipe16_6_ChebyshevApproximation(std::ostream& out)
	{
		auto target = [](Real x) { return std::exp(x); };
		ChebyshevApproximation approximation(target, -1.0, 1.0, 12);
		Polynom<Real> powerSeries = approximation.ToPolynomial();
		const Real x = 0.35;
		Real exact = target(x);
		Real chebError = std::abs(approximation(x) - exact);
		Real powerError = std::abs(powerSeries(x) - exact);
		Real maxError = approximation.MaxError(target, 101);

		out << "  [16.6] exp(x) Chebyshev terms=" << approximation.NumTerms()
		    << ", power degree=" << powerSeries.degree()
		    << ", error at " << x << " cheb/power=" << chebError << "/" << powerError
		    << ", max sampled error=" << maxError << "\n";
	}

	void Cookbook_Recipe16_Polynoms(std::ostream& out)
	{
		Cookbook_Recipe16_1_CoefficientsAndEvaluation(out);
		Cookbook_Recipe16_2_ArithmeticModel(out);
		Cookbook_Recipe16_3_DivisionAndResidual(out);
		Cookbook_Recipe16_4_Calculus(out);
		Cookbook_Recipe16_5_FromValues(out);
		Cookbook_Recipe16_6_ChebyshevApproximation(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 17: Quaternions
	///////////////////////////////////////////////////////////////////////////

	// 17.1 Rotate a vector around an arbitrary axis
	void Cookbook_Recipe17_1_AxisAngleRotation(std::ostream& out)
	{
		Vec3Cart axis = Vec3Cart(1.0, 1.0, 0.0).Normalized();
		Quaternion rotation = Quaternion::FromAxisAngle(axis, Constants::PI / 2.0);
		Vec3Cart original(0.0, 0.0, 1.0);
		Vec3Cart rotated = rotation.Rotate(original);

		out << "  [17.1] rotate " << original << " by 90 deg around " << axis
		    << " -> " << rotated << ", length = " << rotated.NormL2() << "\n";
	}

	Quaternion RotationBetweenDirections(const Vec3Cart& fromDirection, const Vec3Cart& toDirection)
	{
		Vec3Cart from = fromDirection.Normalized();
		Vec3Cart to = toDirection.Normalized();
		Real dot = std::clamp(from[0] * to[0] + from[1] * to[1] + from[2] * to[2], -1.0, 1.0);
		if (dot > 1.0 - 1e-12)
			return Quaternion::Identity();
		if (dot < -1.0 + 1e-12) {
			Vec3Cart helper = std::abs(from[0]) < 0.9 ? Vec3Cart(1, 0, 0) : Vec3Cart(0, 1, 0);
			Vec3Cart axis = VectorProduct(from, helper).Normalized();
			return Quaternion::FromAxisAngle(axis, Constants::PI);
		}

		Vec3Cart axis = VectorProduct(from, to).Normalized();
		return Quaternion::FromAxisAngle(axis, std::acos(dot));
	}

	// 17.2 Point one direction toward another
	void Cookbook_Recipe17_2_AlignDirections(std::ostream& out)
	{
		Vec3Cart forward(1.0, 0.0, 0.0);
		Vec3Cart target = Vec3Cart(1.0, 2.0, 1.0).Normalized();
		Quaternion aim = RotationBetweenDirections(forward, target);
		Vec3Cart aligned = aim.Rotate(forward);
		Vec3Cart opposite = RotationBetweenDirections(forward, -forward).Rotate(forward);

		out << "  [17.2] align forward to target: " << aligned
		    << " (target " << target << "), opposite case -> " << opposite << "\n";
	}

	// 17.3 Compose rotations and observe order
	void Cookbook_Recipe17_3_ComposeRotations(std::ostream& out)
	{
		Quaternion yawZ = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), Constants::PI / 2.0);
		Quaternion rollX = Quaternion::FromAxisAngle(Vec3Cart(1, 0, 0), Constants::PI / 2.0);
		Vec3Cart vector(1.0, 0.0, 0.0);

		Vec3Cart yawThenRoll = (rollX * yawZ).Rotate(vector);
		Vec3Cart rollThenYaw = (yawZ * rollX).Rotate(vector);

		out << "  [17.3] yaw then roll = " << yawThenRoll
		    << ", roll then yaw = " << rollThenYaw << "\n";
	}

	// 17.4 Convert among Euler angles, quaternion, and rotation matrix
	void Cookbook_Recipe17_4_Conversions(std::ostream& out)
	{
		const Real deg = Constants::PI / 180.0;
		Quaternion orientation = Quaternion::FromEulerZYX(30.0 * deg, 20.0 * deg, 10.0 * deg);
		Vec3Cart recoveredEuler = orientation.ToEulerZYX();
		MatrixNM<Real, 3, 3> matrix = orientation.ToRotationMatrix();
		Quaternion reconstructed = Quaternion::FromRotationMatrix(matrix);
		Vec3Cart probe(1.0, 2.0, 3.0);

		out << "  [17.4] recovered yaw/pitch/roll = ("
		    << recoveredEuler[0] / deg << ", " << recoveredEuler[1] / deg << ", "
		    << recoveredEuler[2] / deg << ") deg, matrix round-trip probe = "
		    << reconstructed.Rotate(probe) << "\n";
	}

	// 17.5 Undo rotations and calculate relative orientation
	void Cookbook_Recipe17_5_InverseAndRelative(std::ostream& out)
	{
		const Real deg = Constants::PI / 180.0;
		Quaternion from = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), 30.0 * deg);
		Quaternion to = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), 100.0 * deg);
		Quaternion relative = to * from.Inverse();
		Vec3Cart vector(2.0, -1.0, 3.0);
		Vec3Cart restored = from.Inverse().Rotate(from.Rotate(vector));

		out << "  [17.5] inverse restores " << restored
		    << ", relative orientation = " << relative.GetRotationAngle() / deg << " deg around "
		    << relative.GetRotationAxis() << "\n";
	}

	// 17.6 Smooth interpolation with SLERP
	void Cookbook_Recipe17_6_Slerp(std::ostream& out)
	{
		const Real deg = Constants::PI / 180.0;
		Quaternion start = Quaternion::Identity();
		Quaternion finish = Quaternion::FromAxisAngle(Vec3Cart(0, 0, 1), 120.0 * deg);
		Quaternion quarter = Quaternion::Slerp(start, finish, 0.25);
		Quaternion halfway = Quaternion::Slerp(start, finish, 0.5);
		Quaternion lerpQuarter = Quaternion::Lerp(start, finish, 0.25).Normalized();

		out << "  [17.6] SLERP angles at t=0.25/0.5 = "
		    << quarter.GetRotationAngle() / deg << "/" << halfway.GetRotationAngle() / deg
		    << " deg; normalized LERP at t=0.25 = " << lerpQuarter.GetRotationAngle() / deg << " deg\n";
	}

	void Cookbook_Recipe17_Quaternions(std::ostream& out)
	{
		Cookbook_Recipe17_1_AxisAngleRotation(out);
		Cookbook_Recipe17_2_AlignDirections(out);
		Cookbook_Recipe17_3_ComposeRotations(out);
		Cookbook_Recipe17_4_Conversions(out);
		Cookbook_Recipe17_5_InverseAndRelative(out);
		Cookbook_Recipe17_6_Slerp(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 18: Analyzing functions
	///////////////////////////////////////////////////////////////////////////

	template<typename Container>
	std::string AnalyzerValues(const Container& values)
	{
		std::ostringstream text;
		text << "[";
		for (size_t i = 0; i < values.size(); ++i) {
			if (i > 0) text << ", ";
			text << values[i];
		}
		text << "]";
		return text.str();
	}

	const char* DiscontinuityName(DiscontinuityType type)
	{
		switch (type) {
		case DiscontinuityType::JUMP: return "jump";
		case DiscontinuityType::REMOVABLE: return "removable";
		case DiscontinuityType::INFINITE: return "infinite";
		case DiscontinuityType::OSCILLATORY: return "oscillatory";
		default: return "unknown";
		}
	}

	// 18.1 Produce point and interval health reports
	void Cookbook_Recipe18_1_HealthReport(std::ostream& out)
	{
		RealFunctionFromStdFunc logarithm([](Real x) { return std::log(x); });
		RealFunctionAnalyzer analyzer(logarithm, "log(x)");
		auto validPoint = analyzer.AnalyzePoint(1.0);
		auto invalidPoint = analyzer.AnalyzePoint(-1.0);
		auto interval = analyzer.AnalyzeInterval(0.25, 4.0, 200);
		Real minimum = analyzer.MinInNPoints(0.25, 4.0, 200);
		Real maximum = analyzer.MaxInNPoints(0.25, 4.0, 200);

		RealFunctionFromStdFunc exponential([](Real x) { return std::exp(x); });
		auto overflow = SafeEvaluate(exponential, 1000.0);

		out << "  [18.1] log health: x=1 defined/continuous/derivative="
		    << validPoint.isDefined << "/" << validPoint.isContinuous << "/"
		    << validPoint.isDerivativeDefined << ", x=-1 defined=" << invalidPoint.isDefined
		    << ", [0.25,4] monotonic=" << interval.isMonotonic << ", sampled range=["
		    << minimum << ", " << maximum << "], exp(1000) overflow=" << overflow.isOverflow << "\n";
	}

	// 18.2 Build a complete polynomial feature map
	void Cookbook_Recipe18_2_PolynomialFeatures(std::ostream& out)
	{
		RealFunctionFromStdFunc polynomial([](Real x) { return x*x*x - 3.0*x; });
		RealFunctionAnalyzer analyzer(polynomial, "x^3-3x");
		auto roots = analyzer.GetRoots(-2.0, 2.0, 1e-7);
		auto extrema = analyzer.GetLocalOptimumsClassified(-2.0, 2.0, 1e-4);
		auto inflections = analyzer.GetInflectionPoints(-2.0, 2.0, 1e-4);
		std::vector<Real> extremaLocations;
		for (const auto& point : extrema) extremaLocations.push_back(point.x);

		bool outerLeftMonotonic = analyzer.isMonotonic(-2.0, -1.0, 100);
		bool wholeIntervalMonotonic = analyzer.isMonotonic(-2.0, 2.0, 100);

		out << "  [18.2] x^3-3x roots=" << AnalyzerValues(roots)
		    << ", extrema=" << AnalyzerValues(extremaLocations)
		    << ", inflections=" << AnalyzerValues(inflections)
		    << ", monotonic outer/whole=" << outerLeftMonotonic << "/" << wholeIntervalMonotonic << "\n";
	}

	// 18.3 Detect and classify discontinuities
	void Cookbook_Recipe18_3_Discontinuities(std::ostream& out)
	{
		RealFunctionFromStdFunc step([](Real x) { return x < 0.0 ? 0.0 : 1.0; });
		RealFunctionFromStdFunc removable([](Real x) {
			if (x == 1.0) return std::numeric_limits<Real>::quiet_NaN();
			return (x*x - 1.0) / (x - 1.0);
		});
		RealFunctionFromStdFunc pole([](Real x) { return 1.0 / (x - 2.0); });
		RealFunctionAnalyzer stepAnalyzer(step);
		RealFunctionAnalyzer removableAnalyzer(removable);
		RealFunctionAnalyzer poleAnalyzer(pole);
		auto jump = stepAnalyzer.ClassifyDiscontinuity(0.0);
		auto hole = removableAnalyzer.ClassifyDiscontinuity(1.0);
		auto infinite = poleAnalyzer.ClassifyDiscontinuity(2.0, 1e-8);
		auto discovered = stepAnalyzer.FindDiscontinuities(-1.0, 1.0, 100);

		out << "  [18.3] discontinuities jump/removable/pole="
		    << DiscontinuityName(jump.type) << "/" << DiscontinuityName(hole.type) << "/"
		    << DiscontinuityName(infinite.type) << ", jump limits=" << jump.leftLimit << "/"
		    << jump.rightLimit << ", removable limit=" << hole.leftLimit
		    << ", interval scan found=" << discovered.size() << "\n";
	}

	// 18.4 Estimate oscillation period from zero crossings
	void Cookbook_Recipe18_4_PeriodFromZeros(std::ostream& out)
	{
		RealFunctionFromStdFunc damped([](Real t) {
			return std::exp(-0.05 * t) * std::sin(5.0 * t);
		});
		RealFunctionAnalyzer analyzer(damped, "exp(-0.05t)sin(5t)");
		auto roots = analyzer.GetRoots(0.2, 5.5, 1e-7);
		Real zeroSpacing = analyzer.calcRootsPeriod(0.2, 5.5, 1000);
		Real inferredPeriod = 2.0 * zeroSpacing;

		out << "  [18.4] damped oscillation roots=" << roots.size()
		    << ", average zero spacing=" << zeroSpacing
		    << ", inferred full period=" << inferredPeriod
		    << " (exact " << 2.0 * Constants::PI / 5.0 << ")\n";
	}

	// 18.5 Compare an approximation against a reference function
	void Cookbook_Recipe18_5_CompareApproximation(std::ostream& out)
	{
		RealFunctionFromStdFunc exact([](Real x) { return 9.0 - (x - 3.0) * (x - 3.0); });
		Vector<Real> nodes({0.0, 3.0, 6.0});
		Vector<Real> values({exact(0.0), exact(3.0), exact(6.0)});
		LinearInterpRealFunc approximation(nodes, values);
		RealFunctionComparer comparer(approximation, exact);

		Real averageAbsolute = comparer.getAbsDiffAvg(0.0, 6.0, 1000);
		Real maximumAbsolute = comparer.getAbsDiffMax(0.0, 6.0, 1000);
		Real signedArea = comparer.getIntegratedDiff(0.0, 6.0, IntegrationMethod::ROMBERG);
		Real absoluteArea = comparer.getIntegratedAbsDiff(0.0, 6.0, IntegrationMethod::ROMBERG);
		Real squaredError = comparer.getIntegratedSqrDiff(0.0, 6.0, IntegrationMethod::ROMBERG);

		out << "  [18.5] linear-vs-parabola error: avg=" << averageAbsolute
		    << ", max=" << maximumAbsolute << ", signed integral=" << signedArea
		    << ", absolute integral=" << absoluteArea << ", squared integral=" << squaredError << "\n";
	}

	void Cookbook_Recipe18_AnalyzingFunctions(std::ostream& out)
	{
		Cookbook_Recipe18_1_HealthReport(out);
		Cookbook_Recipe18_2_PolynomialFeatures(out);
		Cookbook_Recipe18_3_Discontinuities(out);
		Cookbook_Recipe18_4_PeriodFromZeros(out);
		Cookbook_Recipe18_5_CompareApproximation(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 23: Fourier algorithms
	///////////////////////////////////////////////////////////////////////////

	Real MaxComplexDifference(const Vector<Complex>& first, const Vector<Complex>& second)
	{
		Real maximum = 0.0;
		for (int i = 0; i < first.size(); ++i)
			maximum = std::max(maximum, static_cast<Real>(std::abs(first[i] - second[i])));
		return maximum;
	}

	Real RootMeanSquareError(const Vector<Real>& first, const Vector<Real>& second)
	{
		Real sum = 0.0;
		for (int i = 0; i < first.size(); ++i) {
			Real difference = first[i] - second[i];
			sum += difference * difference;
		}
		return std::sqrt(sum / first.size());
	}

	// 23.1 FFT round trip and normalization conventions
	void Cookbook_Recipe23_1_RoundTripAndNormalization(std::ostream& out)
	{
		Vector<Complex> signal(16);
		for (int i = 0; i < signal.size(); ++i)
			signal[i] = Complex(std::sin(0.3 * i), 0.25 * std::cos(0.7 * i));

		auto legacySpectrum = Fourier::FFT::Forward(signal);
		auto legacyRecovered = Fourier::FFT::Inverse(legacySpectrum);
		auto orthonormalSpectrum = Fourier::FFT::Forward(signal, Fourier::TransformNormalization::Orthonormal);
		auto orthonormalRecovered = Fourier::FFT::Inverse(
			orthonormalSpectrum, Fourier::TransformNormalization::Orthonormal);

		Vector<Complex> raw = signal;
		Fourier::FFT::Transform(raw, 1);
		Fourier::FFT::Transform(raw, -1);
		for (int i = 0; i < raw.size(); ++i) raw[i] /= raw.size();

		out << "  [23.1] FFT round-trip errors: legacy=" << MaxComplexDifference(signal, legacyRecovered)
		    << ", orthonormal=" << MaxComplexDifference(signal, orthonormalRecovered)
		    << ", raw in-place after /N=" << MaxComplexDifference(signal, raw) << "\n";
	}

	// 23.2 Detect frequencies and amplitudes in a real signal
	void Cookbook_Recipe23_2_RealSignalSpectrum(std::ostream& out)
	{
		const int sampleCount = 256;
		const Real sampleRate = 256.0;
		Vector<Real> signal(sampleCount);
		for (int n = 0; n < sampleCount; ++n) {
			Real time = n / sampleRate;
			signal[n] = std::sin(2.0 * Constants::PI * 20.0 * time)
				+ 0.4 * std::sin(2.0 * Constants::PI * 60.0 * time);
		}

		auto spectrum = Fourier::RealFFT::Forward(signal);
		auto frequencies = Fourier::FrequencyAxis(
			sampleCount, sampleRate, Fourier::FrequencyAxisLayout::OneSided);
		std::vector<std::pair<Real, int>> peaks;
		for (int bin = 1; bin < spectrum.size() - 1; ++bin)
			peaks.push_back({ 2.0 * std::abs(spectrum[bin]) / sampleCount, bin });
		std::partial_sort(peaks.begin(), peaks.begin() + 2, peaks.end(), std::greater<>());

		out << "  [23.2] dominant real-signal tones: " << frequencies[peaks[0].second]
		    << " Hz @ " << peaks[0].first << ", " << frequencies[peaks[1].second]
		    << " Hz @ " << peaks[1].first << "; stored bins=" << spectrum.size() << "\n";
	}

	// 23.3 Compare spectral leakage and window tradeoffs
	void Cookbook_Recipe23_3_Windowing(std::ostream& out)
	{
		const int sampleCount = 256;
		const Real frequency = 20.5;
		Vector<Real> signal(sampleCount);
		for (int n = 0; n < sampleCount; ++n)
			signal[n] = std::sin(2.0 * Constants::PI * frequency * n / sampleCount);

		auto leakageOutside = [&](const Vector<Real>& window, int halfWidth) {
			auto metrics = Fourier::Windows::Metrics(window);
			auto spectrum = Fourier::RealFFT::Forward(Fourier::Windows::ApplyWindow(signal, window));
			int peak = 1;
			for (int bin = 2; bin < spectrum.size(); ++bin)
				if (std::abs(spectrum[bin]) > std::abs(spectrum[peak])) peak = bin;
			Real total = 0.0, outside = 0.0;
			for (int bin = 0; bin < spectrum.size(); ++bin) {
				Real power = std::norm(spectrum[bin]);
				total += power;
				if (std::abs(bin - peak) > halfWidth) outside += power;
			}
			Real correctedAmplitude = 2.0 * std::abs(spectrum[peak])
				/ (sampleCount * metrics.coherent_gain);
			return std::pair<Real, Real>{ outside / total, correctedAmplitude };
		};

		auto rectangular = leakageOutside(Fourier::Windows::Rectangular(sampleCount), 2);
		auto hann = leakageOutside(Fourier::Windows::Hann(sampleCount), 2);
		auto flatTopMetrics = Fourier::Windows::Metrics(Fourier::Windows::FlatTop(sampleCount));

		out << "  [23.3] off-bin tone leakage outside +/-2 bins: rectangular=" << rectangular.first
		    << ", Hann=" << hann.first << "; corrected peak amplitudes=" << rectangular.second
		    << "/" << hann.second << ", flat-top ENBW=" << flatTopMetrics.enbw_bins << " bins\n";
	}

	// 23.4 Remove interference with frequency-domain filtering
	void Cookbook_Recipe23_4_FrequencyDomainFiltering(std::ostream& out)
	{
		const int sampleCount = 256;
		const Real sampleRate = 256.0;
		Vector<Real> clean(sampleCount), noisy(sampleCount);
		for (int n = 0; n < sampleCount; ++n) {
			Real time = n / sampleRate;
			clean[n] = std::sin(2.0 * Constants::PI * 10.0 * time);
			noisy[n] = clean[n] + 0.7 * std::sin(2.0 * Constants::PI * 40.0 * time)
				+ 0.2 * std::sin(2.0 * Constants::PI * 70.0 * time);
		}

		auto spectrum = Fourier::RealFFT::Forward(noisy);
		auto frequencies = Fourier::FrequencyAxis(
			sampleCount, sampleRate, Fourier::FrequencyAxisLayout::OneSided);
		for (int bin = 0; bin < spectrum.size(); ++bin)
			if (frequencies[bin] > 20.0) spectrum[bin] = Complex{};
		auto filtered = Fourier::RealFFT::Inverse(spectrum);

		out << "  [23.4] low-pass RMSE: noisy=" << RootMeanSquareError(noisy, clean)
		    << ", filtered=" << RootMeanSquareError(filtered, clean) << "\n";
	}

	// 23.5 Smooth and shape signals with convolution
	void Cookbook_Recipe23_5_Convolution(std::ostream& out)
	{
		Vector<Real> signal({ 1, 2, 3, 4, 5, 4, 3, 2 });
		Vector<Real> kernel({ 0.25, 0.5, 0.25 });
		auto full = Fourier::Convolve(signal, kernel, Fourier::ConvolutionMode::Full);
		auto sameDirect = Fourier::Convolve(signal, kernel, Fourier::ConvolutionMode::Same,
			Fourier::ConvolutionMethod::Direct);
		auto sameFFT = Fourier::Convolve(signal, kernel, Fourier::ConvolutionMode::Same,
			Fourier::ConvolutionMethod::FFT);
		auto valid = Fourier::Convolve(signal, kernel, Fourier::ConvolutionMode::Valid);

		Real methodDifference = 0.0;
		for (int i = 0; i < sameDirect.size(); ++i)
			methodDifference = std::max(methodDifference, std::abs(sameDirect[i] - sameFFT[i]));

		out << "  [23.5] convolution sizes full/same/valid=" << full.size() << "/"
		    << sameDirect.size() << "/" << valid.size() << ", smoothed=" << sameDirect
		    << ", direct-vs-FFT max diff=" << methodDifference << "\n";
	}

	// 23.6 Estimate a sample delay from phase slope
	void Cookbook_Recipe23_6_PhaseDelay(std::ostream& out)
	{
		const int sampleCount = 128;
		const int delay = 7;
		Vector<Complex> original(sampleCount, Complex{}), delayed(sampleCount, Complex{});
		original[0] = 1.0;
		delayed[delay] = 1.0;
		auto originalSpectrum = Fourier::FFT::Forward(original);
		auto delayedSpectrum = Fourier::FFT::Forward(delayed);
		Vector<Complex> phaseRatio(sampleCount / 2 + 1);
		for (int bin = 0; bin < phaseRatio.size(); ++bin)
			phaseRatio[bin] = delayedSpectrum[bin] * std::conj(originalSpectrum[bin]);
		auto unwrapped = Fourier::UnwrapPhase(Fourier::Phase(phaseRatio));

		Real sumXX = 0.0, sumXY = 0.0;
		for (int bin = 1; bin < unwrapped.size(); ++bin) {
			sumXX += bin * bin;
			sumXY += bin * unwrapped[bin];
		}
		Real slope = sumXY / sumXX;
		Real estimatedDelay = -slope * sampleCount / (2.0 * Constants::PI);

		out << "  [23.6] phase-slope delay estimate = " << estimatedDelay
		    << " samples (actual " << delay << ")\n";
	}

	// 23.7 Arbitrary-size FFT through Bluestein
	void Cookbook_Recipe23_7_Bluestein(std::ostream& out)
	{
		const int primeLength = 127;
		Vector<Complex> signal(primeLength);
		for (int n = 0; n < primeLength; ++n)
			signal[n] = Complex(std::sin(0.17 * n) + 0.3 * std::cos(0.41 * n), 0.1 * std::sin(0.07 * n));

		auto fft = Fourier::FFT::Forward(signal);
		auto dft = Fourier::DFT::Forward(signal);
		auto recovered = Fourier::FFT::Inverse(fft);

		out << "  [23.7] prime-length N=" << primeLength
		    << " Bluestein-vs-DFT max diff=" << MaxComplexDifference(fft, dft)
		    << ", round-trip error=" << MaxComplexDifference(signal, recovered) << "\n";
	}

	// 23.8 Verify Parseval energy conservation for a real FFT
	void Cookbook_Recipe23_8_Parseval(std::ostream& out)
	{
		const int sampleCount = 128;
		Vector<Real> signal(sampleCount);
		Real timeEnergy = 0.0;
		for (int n = 0; n < sampleCount; ++n) {
			signal[n] = std::sin(2.0 * Constants::PI * 7.0 * n / sampleCount)
				+ 0.5 * std::cos(2.0 * Constants::PI * 19.0 * n / sampleCount);
			timeEnergy += signal[n] * signal[n];
		}
		auto spectrum = Fourier::RealFFT::Forward(signal);
		Real frequencyEnergy = std::norm(spectrum[0]) + std::norm(spectrum[sampleCount / 2]);
		for (int bin = 1; bin < sampleCount / 2; ++bin)
			frequencyEnergy += 2.0 * std::norm(spectrum[bin]);
		frequencyEnergy /= sampleCount;

		out << "  [23.8] Parseval energy time/frequency=" << timeEnergy << "/"
		    << frequencyEnergy << ", difference=" << std::abs(timeEnergy - frequencyEnergy) << "\n";
	}

	// 23.9 Compress smooth data with the DCT
	void Cookbook_Recipe23_9_DCTCompression(std::ostream& out)
	{
		const int sampleCount = 64;
		const int retained = 8;
		Vector<Real> signal(sampleCount);
		for (int n = 0; n < sampleCount; ++n) {
			Real x = (n + 0.5) / sampleCount;
			signal[n] = std::exp(-2.0 * x) + 0.2 * std::cos(2.0 * Constants::PI * x);
		}
		auto coefficients = Fourier::DCT::ForwardII(signal);
		for (int i = retained; i < coefficients.size(); ++i) coefficients[i] = 0.0;
		auto reconstructed = Fourier::DCT::InverseII(coefficients);

		out << "  [23.9] DCT compression retained " << retained << "/" << sampleCount
		    << " coefficients, RMSE=" << RootMeanSquareError(signal, reconstructed) << "\n";
	}

	// 23.10 Analyze zero-boundary sine modes with DST-I
	void Cookbook_Recipe23_10_DSTModes(std::ostream& out)
	{
		const int sampleCount = 32;
		Vector<Real> signal(sampleCount);
		for (int n = 0; n < sampleCount; ++n) {
			signal[n] = std::sin(Constants::PI * 3.0 * (n + 1) / (sampleCount + 1))
				+ 0.3 * std::sin(Constants::PI * 7.0 * (n + 1) / (sampleCount + 1));
		}
		auto coefficients = Fourier::DCT::ForwardDST(signal);
		auto recovered = Fourier::DCT::InverseDST(coefficients);
		std::vector<std::pair<Real, int>> modes;
		for (int i = 0; i < coefficients.size(); ++i)
			modes.push_back({ std::abs(coefficients[i]), i + 1 });
		std::partial_sort(modes.begin(), modes.begin() + 2, modes.end(), std::greater<>());

		out << "  [23.10] dominant DST modes=" << modes[0].second << "/" << modes[1].second
		    << ", reconstruction RMSE=" << RootMeanSquareError(signal, recovered) << "\n";
	}

	void Cookbook_Recipe23_FourierAlgorithms(std::ostream& out)
	{
		Cookbook_Recipe23_1_RoundTripAndNormalization(out);
		Cookbook_Recipe23_2_RealSignalSpectrum(out);
		Cookbook_Recipe23_3_Windowing(out);
		Cookbook_Recipe23_4_FrequencyDomainFiltering(out);
		Cookbook_Recipe23_5_Convolution(out);
		Cookbook_Recipe23_6_PhaseDelay(out);
		Cookbook_Recipe23_7_Bluestein(out);
		Cookbook_Recipe23_8_Parseval(out);
		Cookbook_Recipe23_9_DCTCompression(out);
		Cookbook_Recipe23_10_DSTModes(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 24: Graph algorithms
	///////////////////////////////////////////////////////////////////////////

	template<typename V, typename E>
	std::string GraphPathNames(const Graph<V, E>& graph, const std::vector<size_t>& path)
	{
		std::ostringstream text;
		for (size_t i = 0; i < path.size(); ++i) {
			if (i > 0) text << " -> ";
			text << graph.vertexData(path[i]);
		}
		return text.str();
	}

	std::string GraphIndices(const std::vector<size_t>& indices)
	{
		std::ostringstream text;
		text << "[";
		for (size_t i = 0; i < indices.size(); ++i) {
			if (i > 0) text << ", ";
			text << indices[i];
		}
		text << "]";
		return text.str();
	}

	// 24.1 Build, traverse, and inspect a graph
	void Cookbook_Recipe24_1_TraversalAndComponents(std::ostream& out)
	{
		Graph<std::string> network(Graph<std::string>::Type::Undirected);
		for (const char* name : { "A", "B", "C", "D", "E", "F" })
			network.addVertex(name);
		network.addEdge(0, 1);
		network.addEdge(1, 2);
		network.addEdge(0, 2);
		network.addEdge(3, 4);

		auto bfs = BFS(network, 0);
		auto dfs = DFS(network, 0);
		auto components = ConnectedComponents(network);

		out << "  [24.1] BFS = " << GraphPathNames(network, bfs.visitOrder)
		    << ", DFS = " << GraphPathNames(network, dfs.visitOrder)
		    << ", connected components = " << components.numComponents << "\n";
	}

	// 24.2 Find best routes with BFS, Dijkstra, and A*
	void Cookbook_Recipe24_2_BestRoutes(std::ostream& out)
	{
		Graph<std::string> roads(Graph<std::string>::Type::Undirected);
		for (const char* name : { "Depot", "Museum", "Park", "Station", "Harbor" })
			roads.addVertex(name);
		roads.addEdge(0, 1, 2.0);
		roads.addEdge(0, 2, 4.0);
		roads.addEdge(1, 2, 1.0);
		roads.addEdge(1, 3, 5.0);
		roads.addEdge(2, 3, 1.5);
		roads.addEdge(3, 4, 2.0);
		roads.addEdge(2, 4, 6.0);

		std::vector<Vec3Cart> position{
			{0, 0, 0}, {2, 0, 0}, {2, 1, 0}, {3, 2, 0}, {5, 2, 0}
		};
		auto heuristic = [&](size_t vertex, size_t target) {
			return (position[target] - position[vertex]).NormL2();
		};

		auto hops = ShortestPathUnweighted(roads, 0, 4);
		auto dijkstra = DijkstraPath(roads, 0, 4);
		auto astar = AStarPath(roads, 0, 4, heuristic);

		out << "  [24.2] fewest hops: " << GraphPathNames(roads, hops.path)
		    << "; shortest weighted route: " << GraphPathNames(roads, dijkstra.path)
		    << " (" << dijkstra.totalWeight << "), A* explored " << astar.nodesExplored << " nodes\n";
	}

	// 24.3 Schedule dependent tasks and find the critical path
	void Cookbook_Recipe24_3_DAGScheduling(std::ostream& out)
	{
		Graph<std::string> project(Graph<std::string>::Type::Directed);
		for (const char* task : { "Start", "Design", "Procure", "Build", "Test", "Deploy" })
			project.addVertex(task);
		project.addEdge(0, 1, 3.0);
		project.addEdge(0, 2, 2.0);
		project.addEdge(1, 3, 5.0);
		project.addEdge(2, 3, 4.0);
		project.addEdge(3, 4, 2.0);
		project.addEdge(4, 5, 1.0);

		auto order = TopologicalSort(project);
		auto critical = LongestPathDAG(project, 0);
		Graph<std::string> invalid = project;
		invalid.addEdge(5, 1);
		auto cycle = FindCycle(invalid);

		out << "  [24.3] topological order: " << GraphPathNames(project, order.order)
		    << "; critical path: " << GraphPathNames(project, critical.pathTo(5))
		    << " (duration " << critical.distance[5] << "), added back-edge cycle = "
		    << std::boolalpha << cycle.hasCycle << "\n";
	}

	// 24.4 Design the cheapest connected network
	void Cookbook_Recipe24_4_MinimumSpanningTree(std::ostream& out)
	{
		Graph<std::string> sites(Graph<std::string>::Type::Undirected);
		for (const char* site : { "A", "B", "C", "D", "E" })
			sites.addVertex(site);
		sites.addEdge(0, 1, 4.0);
		sites.addEdge(0, 2, 2.0);
		sites.addEdge(1, 2, 1.0);
		sites.addEdge(1, 3, 5.0);
		sites.addEdge(2, 3, 8.0);
		sites.addEdge(2, 4, 10.0);
		sites.addEdge(3, 4, 2.0);

		auto kruskal = Kruskal(sites);
		auto prim = Prim(sites);

		out << "  [24.4] cheapest network = " << kruskal.totalWeight
		    << " using " << kruskal.edges.size() << " edges; Prim agrees = "
		    << (std::abs(prim.totalWeight - kruskal.totalWeight) < 1e-12) << "\n";
	}

	// 24.5 Find single points of failure
	void Cookbook_Recipe24_5_NetworkResilience(std::ostream& out)
	{
		Graph<std::string> routers(Graph<std::string>::Type::Undirected);
		for (const char* router : { "A", "B", "C", "D", "E", "F" })
			routers.addVertex(router);
		routers.addEdge(0, 1);
		routers.addEdge(1, 2);
		routers.addEdge(2, 0); // resilient triangle
		routers.addEdge(1, 3); // bridge to second region
		routers.addEdge(3, 4);
		routers.addEdge(4, 5);
		routers.addEdge(5, 3); // resilient triangle

		auto resilience = UndirectedConnectivity(routers);

		out << "  [24.5] articulation routers = " << GraphIndices(resilience.articulationPoints)
		    << ", bridges = " << resilience.bridges.size()
		    << ", biconnected regions = " << resilience.biconnectedComponents.size() << "\n";
	}

	// 24.6 Detect cycles and strongly connected subsystems
	void Cookbook_Recipe24_6_StrongComponents(std::ostream& out)
	{
		Graph<std::string> services(Graph<std::string>::Type::Directed);
		for (const char* service : { "API", "Auth", "Users", "Billing", "Email", "Audit" })
			services.addVertex(service);
		services.addEdge(0, 1);
		services.addEdge(1, 2);
		services.addEdge(2, 0); // mutually dependent subsystem
		services.addEdge(2, 3);
		services.addEdge(3, 4);
		services.addEdge(4, 3); // second SCC
		services.addEdge(4, 5);

		auto cycle = FindCycle(services);
		auto components = CondensationDAG(services);

		out << "  [24.6] cycle = " << GraphPathNames(services, cycle.cycle)
		    << ", SCCs = " << components.numComponents
		    << ", condensation edges = " << components.condensationEdges.size() << "\n";
	}

	// 24.7 Calculate maximum throughput and bottleneck cut
	void Cookbook_Recipe24_7_MaxFlowMinCut(std::ostream& out)
	{
		Graph<std::string> pipes(Graph<std::string>::Type::Directed);
		for (const char* node : { "Source", "A", "B", "C", "D", "Sink" })
			pipes.addVertex(node);
		pipes.addEdge(0, 1, 16);
		pipes.addEdge(0, 2, 13);
		pipes.addEdge(1, 2, 10);
		pipes.addEdge(2, 1, 4);
		pipes.addEdge(1, 3, 12);
		pipes.addEdge(2, 4, 14);
		pipes.addEdge(3, 2, 9);
		pipes.addEdge(4, 3, 7);
		pipes.addEdge(3, 5, 20);
		pipes.addEdge(4, 5, 4);

		auto flow = Dinic(pipes, 0, 5);

		out << "  [24.7] maximum flow = " << flow.maxFlow
		    << ", minimum cut has " << flow.cutEdges.size() << " edges; source side = "
		    << GraphIndices(flow.sourceSide) << "\n";
	}

	// 24.8 Assign jobs to workers with bipartite matching
	void Cookbook_Recipe24_8_BipartiteMatching(std::ostream& out)
	{
		Graph<std::string> assignments(Graph<std::string>::Type::Undirected);
		for (const char* name : { "Ana", "Bo", "Cy", "Weld", "Inspect", "Pack" })
			assignments.addVertex(name);
		assignments.addEdge(0, 3);
		assignments.addEdge(0, 4);
		assignments.addEdge(1, 3);
		assignments.addEdge(1, 5);
		assignments.addEdge(2, 4);
		assignments.addEdge(2, 5);

		auto matching = HopcroftKarp(assignments, std::vector<size_t>{ 0, 1, 2 });

		out << "  [24.8] assigned " << matching.cardinality << " workers:";
		for (const auto& [worker, job] : matching.matching)
			out << " " << assignments.vertexData(worker) << "->" << assignments.vertexData(job);
		out << "\n";
	}

	void Cookbook_Recipe24_GraphAlgorithms(std::ostream& out)
	{
		Cookbook_Recipe24_1_TraversalAndComponents(out);
		Cookbook_Recipe24_2_BestRoutes(out);
		Cookbook_Recipe24_3_DAGScheduling(out);
		Cookbook_Recipe24_4_MinimumSpanningTree(out);
		Cookbook_Recipe24_5_NetworkResilience(out);
		Cookbook_Recipe24_6_StrongComponents(out);
		Cookbook_Recipe24_7_MaxFlowMinCut(out);
		Cookbook_Recipe24_8_BipartiteMatching(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 25: DAE solvers
	///////////////////////////////////////////////////////////////////////////

	void Cookbook_Recipe25_1_Index1DAE(std::ostream& out)
	{
		Cookbook_LinearDAE system;
		Vector<Real> x0{ REAL(1.0) };
		Vector<Real> y0{ REAL(0.0) };
		DAESolverConfig config;
		config.step_size = REAL(0.05);
		config.constraint_tol = REAL(1e-10);

		auto result = SolveDAEBackwardEuler(system, REAL(0.0), x0, y0, REAL(1.0), config);
		Vector<Real> xEnd = result.solution.getXValuesAtEnd();
		Vector<Real> yEnd = result.solution.getYValuesAtEnd();

		out << "  [25.1] Backward Euler status=" << static_cast<int>(result.status)
		    << ", x(1)=" << xEnd[0]
		    << ", y(1)=" << yEnd[0]
		    << ", |g|=" << result.final_constraint_norm << "\n";
	}

	void Cookbook_Recipe25_2_ConsistentIC(std::ostream& out)
	{
		Cookbook_LinearDAE system;
		Vector<Real> x0{ REAL(0.8) };
		Vector<Real> yGuess{ REAL(0.0) };
		auto ic = ComputeConsistentICDetailed(system, REAL(0.0), x0, yGuess, 20, REAL(1e-12));
		bool verified = VerifyConsistentIC(system, REAL(0.0), x0, yGuess, REAL(1e-10));

		out << "  [25.2] consistent IC converged=" << std::boolalpha << ic.converged
		    << ", y0=" << yGuess[0]
		    << ", residual=" << ic.final_residual_norm
		    << ", verified=" << verified << "\n";
	}

	void Cookbook_Recipe25_3_SolverComparison(std::ostream& out)
	{
		Cookbook_LinearDAE system;
		Vector<Real> x0{ REAL(1.0) };
		Vector<Real> y0{ REAL(0.0) };
		DAESolverConfig config;
		config.step_size = REAL(0.05);
		config.constraint_tol = REAL(1e-10);

		auto backwardEuler = SolveDAEBackwardEuler(system, REAL(0.0), x0, y0, REAL(1.0), config);
		auto bdf2 = SolveDAEBDF2(system, REAL(0.0), x0, y0, REAL(1.0), config);
		auto bdf4 = SolveDAEBDF4(system, REAL(0.0), x0, y0, REAL(1.0), config);
		auto radau = SolveDAERadauIIA(system, REAL(0.0), x0, y0, REAL(1.0), config);

		out << "  [25.3] final x: BE=" << backwardEuler.solution.getXValuesAtEnd()[0]
		    << ", BDF2=" << bdf2.solution.getXValuesAtEnd()[0]
		    << ", BDF4=" << bdf4.solution.getXValuesAtEnd()[0]
		    << ", Radau=" << radau.solution.getXValuesAtEnd()[0]
		    << "; max |g| Radau=" << radau.max_constraint_violation << "\n";
	}

	void Cookbook_Recipe25_4_EventsAndRestart(std::ostream& out)
	{
		Cookbook_TimedDAEEvent system(REAL(0.5));
		Vector<Real> x0{ REAL(0.0) };
		Vector<Real> y0{ REAL(1.0) };
		DAESolverConfig config;
		config.step_size = REAL(0.2);
		config.constraint_tol = REAL(1e-10);

		auto result = SolveDAEBackwardEulerWithEvents(system, REAL(0.0), x0, y0, REAL(1.0), config);

		out << "  [25.4] events=" << result.events.size()
		    << ", event t=" << (result.events.empty() ? REAL(-1.0) : result.events[0].time)
		    << ", final x=" << result.final_differential_state[0]
		    << ", final y=" << result.final_algebraic_state[0]
		    << ", |g|=" << result.integration.final_constraint_norm << "\n";
	}

	void Cookbook_Recipe25_DAESolvers(std::ostream& out)
	{
		Cookbook_Recipe25_1_Index1DAE(out);
		Cookbook_Recipe25_2_ConsistentIC(out);
		Cookbook_Recipe25_3_SolverComparison(out);
		Cookbook_Recipe25_4_EventsAndRestart(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 26: Algebra
	///////////////////////////////////////////////////////////////////////////

	// 26.1 Reorder real data with permutations
	void Cookbook_Recipe26_1_DataPermutation(std::ostream& out)
	{
		std::array<std::string, 4> wiredChannels = { "pressure", "temperature", "humidity", "voltage" };
		Algebra::Permutation<4> wiringFix = Algebra::Permutation<4>::from_cycles({ {0, 2, 1} });
		auto correctedChannels = wiringFix.apply(wiredChannels);
		auto rawChannelsAgain = wiringFix.inverse().apply(correctedChannels);

		out << "  [26.1] corrected channel order = " << Cookbook_FormatIntContainer(wiringFix.images())
		    << ", corrected[0] = " << correctedChannels[0]
		    << ", inverse restores[0] = " << rawChannelsAgain[0]
		    << ", order = " << wiringFix.order()
		    << ", sign = " << wiringFix.sign() << "\n";
	}

	// 26.2 Model polygon symmetries with a dihedral group
	void Cookbook_Recipe26_2_DihedralSymmetry(std::ostream& out)
	{
		Algebra::DihedralGroup square(4);
		const auto rotation = square.rotation(1);
		const auto reflection = square.reflection();
		const auto relation = square.compose(square.compose(reflection, rotation), reflection);
		const auto rotationInverse = square.inverse(rotation);
		const auto reflectedVertices = square.permutation(reflection).images();
		const auto generated = Algebra::GeneratedSubgroup(square, { rotation, reflection });

		out << "  [26.2] D4 order = " << square.order()
		    << ", reflection maps vertices to " << Cookbook_FormatIntContainer(reflectedVertices)
		    << ", s*r*s = r^-1 -> " << std::boolalpha << (relation == rotationInverse)
		    << ", generated subgroup size = " << generated.size() << "\n";
	}

	// 26.3 Count distinct binary bracelets with Burnside's lemma
	void Cookbook_Recipe26_3_BurnsideBracelets(std::ostream& out)
	{
		Algebra::DihedralGroup braceletSymmetry(6);
		const auto colorings = Cookbook_BinaryColorings(6);
		auto action = Algebra::MakeGroupAction<Algebra::DihedralElement, std::vector<int>>(
			[&braceletSymmetry](const Algebra::DihedralElement& element, const std::vector<int>& coloring) {
				return braceletSymmetry.permutation(element).apply(coloring);
			});

		const std::size_t distinctBracelets = Algebra::BurnsideCount(braceletSymmetry, action, colorings);
		const std::size_t fixedByOneStepRotation = Algebra::FixedPointCount(
			braceletSymmetry, action, braceletSymmetry.rotation(1), colorings);
		const auto orbits = Algebra::OrbitPartition(braceletSymmetry, action, colorings);

		out << "  [26.3] binary bracelets with 6 beads: raw patterns = " << colorings.size()
		    << ", symmetry classes = " << distinctBracelets
		    << ", orbit partition size = " << orbits.size()
		    << ", fixed by one-step rotation = " << fixedByOneStepRotation << "\n";
	}

	// 26.4 Exact modular arithmetic and prime fields
	void Cookbook_Recipe26_4_ModularArithmetic(std::ostream& out)
	{
		using Z12 = Algebra::ModInt<12>;
		using F7 = Algebra::PrimeFieldElement<7>;
		Algebra::ModularRing<12> clockRing;
		Algebra::PrimeField<7> field7;

		Z12 clockValue = Z12(10) + Z12(5);
		bool sixHasInverse = true;
		try {
			(void)Z12(6).inverse();
		}
		catch (const DomainError&) {
			sixHasInverse = false;
		}

		F7 quotient = F7(3) / F7(2);
		F7 primitive = Algebra::PrimitiveRoot<7>();
		bool ringLaws = Algebra::CheckRingLaws(clockRing);
		bool fieldLaws = Algebra::CheckFieldLaws(field7);

		out << "  [26.4] 10+5 mod 12 = " << Cookbook_FormatModInt(clockValue)
		    << ", 6 inverse exists in Z12 = " << std::boolalpha << sixHasInverse
		    << ", 3/2 in F7 = " << Cookbook_FormatModInt(quotient)
		    << ", primitive root F7 = " << Cookbook_FormatModInt(primitive)
		    << ", ring/field laws = " << ringLaws << "/" << fieldLaws << "\n";
	}

	// 26.5 Solve an exact linear system over a finite field
	void Cookbook_Recipe26_5_FiniteFieldLinearSystem(std::ostream& out)
	{
		using F7 = Algebra::PrimeFieldElement<7>;
		MatrixNM<F7, 2, 2> system{{ F7(2), F7(3) }, { F7(4), F7(1) }};
		std::array<F7, 2> rightHandSide = { F7(1), F7(6) };

		const auto inverse = Algebra::FieldMatrixInverse(system);
		const F7 determinant = Algebra::FieldMatrixDeterminant(system);
		const std::array<F7, 2> solution = {
			inverse(0, 0) * rightHandSide[0] + inverse(0, 1) * rightHandSide[1],
			inverse(1, 0) * rightHandSide[0] + inverse(1, 1) * rightHandSide[1]
		};
		const F7 check0 = system(0, 0) * solution[0] + system(0, 1) * solution[1];
		const F7 check1 = system(1, 0) * solution[0] + system(1, 1) * solution[1];

		out << "  [26.5] exact F7 solve: det = " << Cookbook_FormatModInt(determinant)
		    << ", solution = [" << Cookbook_FormatModInt(solution[0]) << ","
		    << Cookbook_FormatModInt(solution[1]) << "]"
		    << ", A*x = [" << Cookbook_FormatModInt(check0) << ","
		    << Cookbook_FormatModInt(check1) << "]\n";
	}

	// 26.6 Build a tiny extension field from a polynomial modulus
	void Cookbook_Recipe26_6_ExtensionField(std::ostream& out)
	{
		using GF4 = Algebra::FiniteFieldElement<2, 2, Cookbook_GF4_Modulus>;
		const GF4 alpha{ F2(0), F2(1) };
		Algebra::ExtensionField<2, 2, Cookbook_GF4_Modulus> field;

		const GF4 alphaSquared = alpha * alpha;
		const GF4 alphaPlusOne = alpha + GF4(1);
		const GF4 inverse = alpha.inverse();
		const bool irreducible = Algebra::IsIrreducible(Cookbook_GF4_Modulus::modulus());
		const bool fieldLaws = Algebra::CheckFieldLaws(field);

		out << "  [26.6] GF(4): alpha^2 = " << Cookbook_FormatGFElement(alphaSquared)
		    << ", alpha+1 = " << Cookbook_FormatGFElement(alphaPlusOne)
		    << ", alpha^-1 = " << Cookbook_FormatGFElement(inverse)
		    << ", modulus irreducible = " << std::boolalpha << irreducible
		    << ", field laws = " << fieldLaws << "\n";
	}

	void Cookbook_Recipe26_Algebra(std::ostream& out)
	{
		Cookbook_Recipe26_1_DataPermutation(out);
		Cookbook_Recipe26_2_DihedralSymmetry(out);
		Cookbook_Recipe26_3_BurnsideBracelets(out);
		Cookbook_Recipe26_4_ModularArithmetic(out);
		Cookbook_Recipe26_5_FiniteFieldLinearSystem(out);
		Cookbook_Recipe26_6_ExtensionField(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 27: Complex analysis
	///////////////////////////////////////////////////////////////////////////

	// 27.1 Numerically differentiate an analytic complex function
	void Cookbook_Recipe27_1_ComplexDerivative(std::ostream& out)
	{
		ComplexFunctionFromStdFunc function([](Complex z) { return std::exp(z) * std::sin(z); });
		Complex point(0.7, 0.4);
		Complex exact = std::exp(point) * (std::sin(point) + std::cos(point));
		Complex secondOrder = Derivation::NDer2Complex(function, point);
		Complex fourthOrder = Derivation::NDer4Complex(function, point);
		DerivativeConfig config;
		config.estimate_error = true;
		auto detailed = Derivation::NDer4ComplexDetailed(function, point, config);

		out << "  [27.1] f'(0.7+0.4i): exact=" << exact
		    << ", NDer2 error=" << std::abs(secondOrder - exact)
		    << ", NDer4 error=" << std::abs(fourthOrder - exact)
		    << ", estimate=" << detailed.error << "\n";
	}

	// 27.2 Integrate along line, arc, circle, and custom contours
	void Cookbook_Recipe27_2_ContourIntegration(std::ostream& out)
	{
		ComplexFunctionFromStdFunc identity([](Complex z) { return z; });
		ComplexFunctionFromStdFunc one([](Complex) { return Complex(1.0, 0.0); });
		ComplexFunctionFromStdFunc square([](Complex z) { return z * z; });
		ComplexFunctionFromStdFunc reciprocal([](Complex z) { return Complex(1.0, 0.0) / z; });

		ComplexAnalysis::LineSegmentContour line(Complex(0, 0), Complex(1, 1));
		ComplexAnalysis::ArcContour quarterArc(Complex(0, 0), 1.0, 0.0, Constants::PI / 2.0);
		ComplexAnalysis::CircleContour unitCircle(Complex(0, 0), 1.0);
		auto lineIntegral = ComplexAnalysis::ContourIntegral(identity, line);
		auto arcIntegral = ComplexAnalysis::ContourIntegral(one, quarterArc);
		auto analyticCircle = ComplexAnalysis::ContourIntegral(square, unitCircle);
		auto singularCircle = ComplexAnalysis::ContourIntegral(reciprocal, unitCircle);

		auto ellipse = [](Real t) { return Complex(2.0 * std::cos(t), std::sin(t)); };
		auto ellipseDerivative = [](Real t) { return Complex(-2.0 * std::sin(t), std::cos(t)); };
		auto custom = ComplexAnalysis::ContourIntegral(
			identity, ellipse, ellipseDerivative, 0.0, 2.0 * Constants::PI);

		out << "  [27.2] contour integrals: line(z dz)=" << lineIntegral.value
		    << ", quarter arc(dz)=" << arcIntegral.value
		    << ", circle(z^2 dz)=" << analyticCircle.value
		    << ", circle(dz/z)=" << singularCircle.value
		    << ", custom ellipse(z dz)=" << custom.value << "\n";
	}

	// 27.3 Compute winding numbers
	void Cookbook_Recipe27_3_WindingNumber(std::ostream& out)
	{
		ComplexAnalysis::CircleContour shifted(Complex(2.0, 1.0), 3.0);
		int center = ComplexAnalysis::WindingNumber(shifted, Complex(2.0, 1.0));
		int offCenterInside = ComplexAnalysis::WindingNumber(shifted, Complex(0.0, 1.0));
		int outside = ComplexAnalysis::WindingNumber(shifted, Complex(6.0, 1.0));

		out << "  [27.3] shifted-circle winding numbers center/inside/outside = "
		    << center << "/" << offCenterInside << "/" << outside << "\n";
	}

	// 27.4 Recover values and derivatives with Cauchy's formula
	void Cookbook_Recipe27_4_CauchyFormula(std::ostream& out)
	{
		ComplexFunctionFromStdFunc exponential([](Complex z) { return std::exp(z); });
		Complex point(0.2, 0.1);
		ComplexAnalysis::CircleContour contour(point, 1.0);
		Complex value = ComplexAnalysis::CauchyIntegralFormula(exponential, contour, point);
		Complex first = ComplexAnalysis::CauchyDerivative(exponential, contour, point, 1);
		Complex second = ComplexAnalysis::CauchyDerivative(exponential, contour, point, 2);
		Complex exact = std::exp(point);

		out << "  [27.4] Cauchy exp(z0) value/first/second errors = "
		    << std::abs(value - exact) << "/" << std::abs(first - exact)
		    << "/" << std::abs(second - exact) << "\n";
	}

	// 27.5 Compute residues and verify the residue theorem
	void Cookbook_Recipe27_5_Residues(std::ostream& out)
	{
		ComplexFunctionFromStdFunc rational([](Complex z) {
			return Complex(1.0, 0.0) / (z * z + Complex(1.0, 0.0));
		});
		Complex poleI(0.0, 1.0);
		Complex poleMinusI(0.0, -1.0);
		Complex residueI = ComplexAnalysis::Residue(rational, poleI, 0.4);
		Complex simpleI = ComplexAnalysis::ResidueSimplePole(rational, poleI);
		Complex residueMinusI = ComplexAnalysis::Residue(rational, poleMinusI, 0.4);
		auto onePoleIntegral = ComplexAnalysis::ContourIntegral(
			rational, ComplexAnalysis::CircleContour(poleI, 0.5));
		auto bothPolesIntegral = ComplexAnalysis::ContourIntegral(
			rational, ComplexAnalysis::CircleContour(Complex(0, 0), 2.0));

		out << "  [27.5] residues at +/-i=" << residueI << "/" << residueMinusI
		    << ", simple-pole check=" << simpleI
		    << ", one-pole integral=" << onePoleIntegral.value
		    << ", both-poles integral=" << bothPolesIntegral.value << "\n";
	}

	// 27.6 Count zeros and poles with the argument principle
	void Cookbook_Recipe27_6_ArgumentPrinciple(std::ostream& out)
	{
		ComplexFunctionFromStdFunc polynomial([](Complex z) { return z * z + Complex(1.0, 0.0); });
		ComplexFunctionFromStdFunc meromorphic([](Complex z) {
			return (z * z * z - Complex(1.0, 0.0)) / (z - Complex(0.5, 0.0));
		});
		int bothRoots = ComplexAnalysis::CountZeros(
			polynomial, ComplexAnalysis::CircleContour(Complex(0, 0), 2.0));
		int upperRoot = ComplexAnalysis::CountZeros(
			polynomial, ComplexAnalysis::CircleContour(Complex(0, 1), 0.5));
		int zerosMinusPoles = ComplexAnalysis::ArgumentPrinciple(
			meromorphic, ComplexAnalysis::CircleContour(Complex(0, 0), 2.0));

		out << "  [27.6] zero counts z^2+1 (large/around i)=" << bothRoots << "/" << upperRoot
		    << ", (z^3-1)/(z-0.5) has N-P=" << zerosMinusPoles << "\n";
	}

	// 27.7 Count roots first, then locate them
	void Cookbook_Recipe27_7_CountAndLocateRoots(std::ostream& out)
	{
		ComplexFunctionFromStdFunc cubic([](Complex z) { return z * z * z - Complex(1.0, 0.0); });
		ComplexAnalysis::CircleContour contour(Complex(0, 0), 2.0);
		int count = ComplexAnalysis::CountZeros(cubic, contour);
		auto root1 = RootFinding::FindRootNewtonComplex(cubic, Complex(1.1, 0.1));
		auto root2 = RootFinding::FindRootNewtonComplex(cubic, Complex(-0.4, 0.8));
		auto root3 = RootFinding::FindRootMuller(cubic, Complex(-0.4, -0.8));

		out << "  [27.7] counted " << count << " roots of z^3-1; located "
		    << root1.root << ", " << root2.root << ", " << root3.root
		    << ", max |f|=" << std::max({ std::abs(root1.function_value),
		       std::abs(root2.function_value), std::abs(root3.function_value) }) << "\n";
	}

	void Cookbook_Recipe27_ComplexAnalysis(std::ostream& out)
	{
		Cookbook_Recipe27_1_ComplexDerivative(out);
		Cookbook_Recipe27_2_ContourIntegration(out);
		Cookbook_Recipe27_3_WindingNumber(out);
		Cookbook_Recipe27_4_CauchyFormula(out);
		Cookbook_Recipe27_5_Residues(out);
		Cookbook_Recipe27_6_ArgumentPrinciple(out);
		Cookbook_Recipe27_7_CountAndLocateRoots(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 30: Combinatorics and number theory
	///////////////////////////////////////////////////////////////////////////

	// 30.1 Select a committee, then assign distinct roles
	void Cookbook_Recipe30_1_CommitteesAndRoles(std::ostream& out)
	{
		const int people = 10;
		const int committeeSize = 4;
		long long committeeCount = Combinatorics::BinomialCoefficient(people, committeeSize);
		long long rolesPerCommittee = Combinatorics::FallingFactorial(committeeSize, 3);

		long long enumeratedCommittees = 0;
		long long enumeratedRoleAssignments = 0;
		std::vector<int> firstCommittee;
		Combinatorics::ForEachCombination(people, committeeSize, [&](const std::vector<int>& committee) {
			if (firstCommittee.empty()) firstCommittee = committee;
			++enumeratedCommittees;
			Combinatorics::ForEachPermutation(committeeSize, [&](const std::vector<int>&) {
				++enumeratedRoleAssignments;
			});
		});

		out << "  [30.1] committees = " << committeeCount << " (enumerated " << enumeratedCommittees
		    << "), 3-role assignments per committee = " << rolesPerCommittee
		    << ", total assignments = " << enumeratedRoleAssignments << "\n";
	}

	// 30.2 Enumerate constrained integer resource allocations
	void Cookbook_Recipe30_2_ResourcePartitions(std::ostream& out)
	{
		const int units = 12;
		long long allPartitions = Combinatorics::PartitionCount(units);
		long long feasible = 0;
		std::vector<int> firstFeasible;
		Combinatorics::ForEachPartition(units, [&](const std::vector<int>& allocation) {
			if (allocation.size() == 3 && allocation.front() <= 6 && allocation.back() >= 2) {
				if (firstFeasible.empty()) firstFeasible = allocation;
				++feasible;
			}
		});

		out << "  [30.2] partitions of 12 = " << allPartitions
		    << ", feasible 3-project allocations (2..6 each) = " << feasible
		    << ", first = " << GraphIndices(std::vector<size_t>(firstFeasible.begin(), firstFeasible.end())) << "\n";
	}

	// 30.3 Count groupings of labeled objects
	void Cookbook_Recipe30_3_LabeledGroups(std::ostream& out)
	{
		const int students = 8;
		long long exactlyThreeTeams = Combinatorics::StirlingSecond(students, 3);
		long long anyNumberOfTeams = Combinatorics::BellNumber(students);
		long long threeNamedTeams = exactlyThreeTeams * Combinatorics::FallingFactorial(3, 3);

		out << "  [30.3] group 8 students: exactly 3 unlabeled teams = " << exactlyThreeTeams
		    << ", 3 named teams = " << threeNamedTeams
		    << ", any number of nonempty teams = " << anyNumberOfTeams << "\n";
	}

	// 30.4 Analyze an integer through its prime structure
	void Cookbook_Recipe30_4_PrimeStructure(std::ostream& out)
	{
		const long long value = 360360;
		auto powers = NumberTheory::FactorizePowers(value);
		const std::uint64_t primeA = 1000000007ull;
		const std::uint64_t primeB = 1000000009ull;
		auto semiprimeFactors = NumberTheory::Factorize(primeA * primeB);

		out << "  [30.4] " << value << " factors:";
		for (const auto& [prime, exponent] : powers)
			out << " " << prime << "^" << exponent;
		out << "; divisors = " << NumberTheory::DivisorCount(value)
		    << ", divisor sum = " << NumberTheory::DivisorSum(value)
		    << ", phi = " << NumberTheory::EulerTotient(value)
		    << ", mu = " << NumberTheory::MobiusMu(value)
		    << "; large semiprime -> " << semiprimeFactors[0] << " * " << semiprimeFactors[1] << "\n";
	}

	// 30.5 Educational RSA-style encryption round trip
	void Cookbook_Recipe30_5_ToyRSA(std::ostream& out)
	{
		const std::uint64_t p = 61, q = 53;
		const std::uint64_t modulus = p * q;
		const long long phi = static_cast<long long>((p - 1) * (q - 1));
		const std::uint64_t publicExponent = 17;
		const std::uint64_t privateExponent = NumberTheory::ModInverse(publicExponent, phi);
		const std::uint64_t message = 65;
		const std::uint64_t encrypted = NumberTheory::ModPow(message, publicExponent, modulus);
		const std::uint64_t decrypted = NumberTheory::ModPow(encrypted, privateExponent, modulus);

		out << "  [30.5] toy RSA: n=" << modulus << ", phi=" << phi << ", e=" << publicExponent
		    << ", d=" << privateExponent << ", " << message << " -> " << encrypted
		    << " -> " << decrypted << " (educational only)\n";
	}

	// 30.6 Synchronize repeating schedules with CRT
	void Cookbook_Recipe30_6_ChineseRemainder(std::ostream& out)
	{
		auto coprime = NumberTheory::ChineseRemainder({ 2, 3, 2 }, { 3, 5, 7 });
		auto nonCoprime = NumberTheory::ChineseRemainder({ 1, 4 }, { 6, 9 });
		long long bezoutX, bezoutY;
		long long gcd = NumberTheory::ExtendedGcd(6, 9, bezoutX, bezoutY);
		bool inconsistentRejected = false;
		try {
			(void)NumberTheory::ChineseRemainder({ 0, 1 }, { 2, 4 });
		} catch (const DomainError&) {
			inconsistentRejected = true;
		}

		out << "  [30.6] CRT: first alignment = " << coprime.first << " mod " << coprime.second
		    << ", non-coprime alignment = " << nonCoprime.first << " mod " << nonCoprime.second
		    << ", gcd(6,9) = " << gcd << " = 6*(" << bezoutX << ") + 9*(" << bezoutY
		    << "), inconsistent rejected = " << std::boolalpha << inconsistentRejected << "\n";
	}

	void Cookbook_Recipe30_CombinatoricsAndNumberTheory(std::ostream& out)
	{
		Cookbook_Recipe30_1_CommitteesAndRoles(out);
		Cookbook_Recipe30_2_ResourcePartitions(out);
		Cookbook_Recipe30_3_LabeledGroups(out);
		Cookbook_Recipe30_4_PrimeStructure(out);
		Cookbook_Recipe30_5_ToyRSA(out);
		Cookbook_Recipe30_6_ChineseRemainder(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 19: Coordinate transformations  (story arc 19 -> 20 -> 21 -> 22)
	///////////////////////////////////////////////////////////////////////////

	// 19.1 Points: Cartesian <-> spherical <-> cylindrical round-trips
	void Cookbook_Recipe19_1_Points(std::ostream& out)
	{
		Vector3Cartesian p{ 1.0, 2.0, 2.0 };            // |p| = 3 exactly

		Vector3Spherical   pSpher = CoordTransfCartToSpher.transf(p);
		Vector3Cylindrical pCyl   = CoordTransfCartToCyl.transf(p);

		Vector3Cartesian backS = CoordTransfSpherToCart.transf(pSpher);
		Vector3Cartesian backC = CoordTransfCylToCart.transf(pCyl);

		out << "  [19.1] cart " << p << " -> spher " << pSpher
		    << " (r = " << pSpher[0] << ", exact 3)\n"
		    << "  [19.1] round-trips: spher " << backS << ", cyl " << backC << "\n";
	}

	// 19.2 Vectors are NOT points: a velocity transforms contravariantly, at a point
	void Cookbook_Recipe19_2_VelocityContravariant(std::ostream& out)
	{
		Vector3Cartesian x_cart{ 1.0, 2.0, 2.0 };
		Vector3Spherical x_spher = CoordTransfCartToSpher.transf(x_cart);

		Vec3Cart v_cart{ 1.0, 1.0, 0.0 };               // |v|^2 = 2

		Vector3Spherical v_spher = CoordTransfCartToSpher.transfVecContravariant(v_cart, x_cart);
		Vector3Cartesian v_back  = CoordTransfSpherToCart.transfVecContravariant(v_spher, x_spher);

		out << "  [19.2] v_cart " << v_cart << " -> v_spher " << v_spher
		    << " -> back " << v_back << "\n";
	}

	// 19.3 Gradients transform by the OTHER rule: covariantly
	void Cookbook_Recipe19_3_GradientCovariant(std::ostream& out)
	{
		Vector3Cartesian p_cart{ 1.0, 2.0, 2.0 };
		Vector3Spherical p_spher = CoordTransfCartToSpher.transf(p_cart);

		// gradient of the 1/r potential, computed natively in spherical coordinates
		ScalarFunction<3> potSpher(Fields::InverseRadialPotentialFieldSpher);
		Vector3Spherical grad_spher = ScalarFieldOperations::GradientSpher(potSpher, p_spher);

		Vector3Cartesian grad_cart = CoordTransfSpherToCart.transfVecCovariant(grad_spher, p_cart);
		Vector3Spherical grad_back = CoordTransfCartToSpher.transfVecCovariant(grad_cart, p_spher);

		out << "  [19.3] grad(spher) " << grad_spher << " -> covariant to cart " << grad_cart
		    << " -> back " << grad_back << "\n";
	}

	void Cookbook_Recipe19_CoordinateTransformations(std::ostream& out)
	{
		Cookbook_Recipe19_1_Points(out);
		Cookbook_Recipe19_2_VelocityContravariant(out);
		Cookbook_Recipe19_3_GradientCovariant(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 20: Fields  (story arc 19 -> 20 -> 21 -> 22)
	///////////////////////////////////////////////////////////////////////////

	// 20.1 One field, two coordinate systems, one physical answer
	void Cookbook_Recipe20_1_SamePhysics(std::ostream& out)
	{
		Vector3Cartesian p_cart{ 1.0, 1.0, 1.0 };
		Vector3Spherical p_spher = CoordTransfCartToSpher.transf(p_cart);

		ScalarFunction<3> potCart(Fields::InverseRadialPotentialFieldCart);
		ScalarFunction<3> potSpher(Fields::InverseRadialPotentialFieldSpher);

		Vector3Cartesian grad_cart  = ScalarFieldOperations::GradientCart<3>(potCart, p_cart);
		Vector3Spherical grad_spher = ScalarFieldOperations::GradientSpher(potSpher, p_spher);

		// transform the spherical gradient to Cartesian - must match grad_cart
		Vector3Cartesian grad_transf = CoordTransfSpherToCart.transfVecCovariant(grad_spher, p_cart);

		out << "  [20.1] grad in cart:  " << grad_cart << "\n"
		    << "  [20.1] grad via spher:" << grad_transf
		    << "  (match: " << std::boolalpha
		    << grad_cart.IsEqualTo(grad_transf, 1e-7) << ")\n";
	}

	// 20.2 The identities every field obeys: curl grad = 0, div curl = 0
	void Cookbook_Recipe20_2_FieldIdentities(std::ostream& out)
	{
		Vector3Cartesian p{ 1.2, -0.7, 0.4 };

		ScalarFunction<3> phi([](const VectorN<Real, 3>& x) {
			return x[0] * x[0] * x[1] + x[2] * x[2] * x[2];
		});
		VectorFunctionFromStdFunc<3> gradPhi(std::function<VectorN<Real, 3>(const VectorN<Real, 3>&)>(
			[&phi](const VectorN<Real, 3>& x) {
				return ScalarFieldOperations::GradientCart<3>(phi, x);
			}));
		Vec3Cart curlOfGrad = VectorFieldOperations::CurlCart(gradPhi, p);

		VectorFunction<3> A([](const VectorN<Real, 3>& x) {
			return VectorN<Real, 3>{ x[0] * x[0] * x[1], x[1] * x[1] * x[2], x[2] * x[2] * x[0] };
		});
		VectorFunctionFromStdFunc<3> curlA(std::function<VectorN<Real, 3>(const VectorN<Real, 3>&)>(
			[&A](const VectorN<Real, 3>& x) {
				return VectorN<Real, 3>(VectorFieldOperations::CurlCart(A, x));
			}));
		Real divOfCurl = VectorFieldOperations::DivCart<3>(curlA, p);

		out << "  [20.2] curl(grad phi) = " << curlOfGrad << " (exact 0)\n"
		    << "  [20.2] div(curl A)    = " << divOfCurl << " (exact 0)\n";
	}

	// 20.3 Empty-space gravity: Laplacian of 1/r vanishes, in BOTH coordinate systems
	void Cookbook_Recipe20_3_HarmonicPotential(std::ostream& out)
	{
		Vector3Cartesian p_cart{ 1.0, 1.0, 1.0 };
		Vector3Spherical p_spher = CoordTransfCartToSpher.transf(p_cart);

		ScalarFunction<3> potCart(Fields::InverseRadialPotentialFieldCart);
		ScalarFunction<3> potSpher(Fields::InverseRadialPotentialFieldSpher);

		Real lapCart  = ScalarFieldOperations::LaplacianCart<3>(potCart, p_cart);
		Real lapSpher = ScalarFieldOperations::LaplacianSpher(potSpher, p_spher);

		out << "  [20.3] Laplacian(1/r): cartesian = " << lapCart
		    << ", spherical = " << lapSpher << " (exact 0 away from origin)\n";
	}

	void Cookbook_Recipe20_Fields(std::ostream& out)
	{
		Cookbook_Recipe20_1_SamePhysics(out);
		Cookbook_Recipe20_2_FieldIdentities(out);
		Cookbook_Recipe20_3_HarmonicPotential(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 21: Tensors  (story arc 19 -> 20 -> 21 -> 22)
	///////////////////////////////////////////////////////////////////////////

	// 21.1 Tensor2 basics: index variance, contraction, evaluation
	void Cookbook_Recipe21_1_TensorBasics(std::ostream& out)
	{
		Tensor2<3> T(1, 1, { 2, 0, 1,
		                     0, 3, 0,
		                     0, 0, 4 });   // mixed tensor: 1 covariant, 1 contravariant index

		VectorN<Real, 3> u{ 1, 0, 0 }, w{ 0, 0, 1 };

		out << "  [21.1] covar indices: " << T.NumCovar() << ", contravar: " << T.NumContravar()
		    << ", trace (contraction) = " << T.Contract() << " (exact 9)"
		    << ", T(u,w) = " << T(u, w) << " (exact 1)\n";
	}

	// 21.2 THE tensor: the metric - what makes components physically meaningful
	void Cookbook_Recipe21_2_MetricMeasures(std::ostream& out)
	{
		Vector3Cartesian x_cart{ 1.0, 2.0, 2.0 };
		Vector3Spherical x_spher = CoordTransfCartToSpher.transf(x_cart);

		Vec3Cart v_cart{ 1.0, 1.0, 0.0 };               // |v|^2 = 2 in Cartesian
		Vector3Spherical v_spher = CoordTransfCartToSpher.transfVecContravariant(v_cart, x_cart);

		MetricTensorSpherical g;
		Tensor2<3> g_ij = g(x_spher);

		Real lenSq = 0.0;                                // |v|^2 = g_ij v^i v^j
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				lenSq += g_ij(i, j) * v_spher[i] * v_spher[j];

		out << "  [21.2] g_rr = " << g_ij(0, 0) << ", g_theta,theta = " << g_ij(1, 1)
		    << " (r^2 = 9), g_phi,phi = " << g_ij(2, 2) << "\n"
		    << "  [21.2] |v|^2 via spherical metric = " << lenSq << " (Cartesian value: 2)\n";
	}

	// 21.3 Raising an index: gradient (covariant) -> force (contravariant)
	void Cookbook_Recipe21_3_RaiseIndex(std::ostream& out)
	{
		Vector3Cartesian p_cart{ 1.0, 1.0, 1.0 };
		Vector3Spherical p_spher = CoordTransfCartToSpher.transf(p_cart);

		ScalarFunction<3> potSpher(Fields::InverseRadialPotentialFieldSpher);
		Vector3Spherical grad_spher = ScalarFieldOperations::GradientSpher(potSpher, p_spher); // covariant

		MetricTensorSpherical g;
		auto g_inv = g.GetContravariantMetric(p_spher);

		VectorN<Real, 3> force;                          // F^i = g^ij grad_j
		for (int i = 0; i < 3; i++) {
			force[i] = 0.0;
			for (int j = 0; j < 3; j++)
				force[i] += g_inv(i, j) * grad_spher[j];
		}

		// contravariant components transform with the velocity rule from 19.2 -
		// carried to Cartesian they must equal the Cartesian gradient (metric = identity there)
		Vector3Cartesian force_cart = CoordTransfSpherToCart.transfVecContravariant(
			Vector3Spherical(force), p_spher);

		ScalarFunction<3> potCart(Fields::InverseRadialPotentialFieldCart);
		Vector3Cartesian grad_cart = ScalarFieldOperations::GradientCart<3>(potCart, p_cart);

		out << "  [21.3] raised F^i (spher) = " << VectorN<Real, 3>(force)
		    << ", to cart = " << force_cart
		    << " (matches grad_cart: " << std::boolalpha
		    << force_cart.IsEqualTo(grad_cart, 1e-7) << ")\n";
	}

	// 21.4 The twist: curvy coordinates, FLAT space - Riemann tensor vanishes
	void Cookbook_Recipe21_4_FlatSpace(std::ostream& out)
	{
		MetricTensorSpherical g;
		Vector3Spherical pos{ 2.0, Constants::PI / 3, Constants::PI / 4 };

		Real Gamma_r_thth = g.GetChristoffelSymbolSecondKind(0, 1, 1, pos);   // = -r
		Real Gamma_th_rth = g.GetChristoffelSymbolSecondKind(1, 0, 1, pos);   // = 1/r

		Tensor4<3> riemann = g.GetRiemannCurvatureTensor(pos);
		Real maxAbs = 0.0;
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 3; j++)
				for (int k = 0; k < 3; k++)
					for (int l = 0; l < 3; l++)
						maxAbs = std::max(maxAbs, std::abs(riemann(i, j, k, l)));

		out << "  [21.4] Christoffels nonzero (Gamma^r_tt = " << Gamma_r_thth
		    << ", Gamma^t_rt = " << Gamma_th_rth << ") - coordinates are curvy...\n"
		    << "  [21.4] ...but max |Riemann| = " << maxAbs
		    << " - the SPACE is flat\n";
	}

	void Cookbook_Recipe21_Tensors(std::ostream& out)
	{
		Cookbook_Recipe21_1_TensorBasics(out);
		Cookbook_Recipe21_2_MetricMeasures(out);
		Cookbook_Recipe21_3_RaiseIndex(out);
		Cookbook_Recipe21_4_FlatSpace(out);
	}

	///////////////////////////////////////////////////////////////////////////
	// Recipe 22: Differential geometry & geodesics  (story arc finale)
	///////////////////////////////////////////////////////////////////////////

	// 22.1 Truly curved: the sphere's induced metric and Theorema Egregium
	void Cookbook_Recipe22_1_InducedMetric(std::ostream& out)
	{
		Surfaces::Sphere sphere(2.0);                   // R = 2, parametrized by (theta, phi)
		DifferentialGeometry::InducedMetric2D g(sphere);

		VectorN<Real, 2> pos{ Constants::PI / 3, 0.5 }; // theta = 60 deg
		Tensor2<2> g_ij = g(pos);

		Real K = DifferentialGeometry::GaussianCurvatureIntrinsic(g, pos);

		out << "  [22.1] induced metric: g_tt = " << g_ij(0, 0) << " (R^2 = 4), g_pp = "
		    << g_ij(1, 1) << " (R^2 sin^2 t = 3)\n"
		    << "  [22.1] Gaussian curvature K = " << K << " (exact 1/R^2 = 0.25)\n";
	}

	// 22.2 Geodesics on the sphere: great circles ('why planes fly over the pole')
	void Cookbook_Recipe22_2_Geodesics(std::ostream& out)
	{
		Surfaces::Sphere sphere(1.0);
		DifferentialGeometry::InducedMetric2D metric(sphere);

		// start on the equator heading north-east at unit speed
		VectorN<Real, 2> pos0{ Constants::PI / 2, 0.0 };
		VectorN<Real, 2> vel0{ -1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0) };

		auto sol = IntegrateGeodesicFixedStep<2>(metric, pos0, vel0, 0.0, 2 * Constants::PI, 2000);

		auto end = sol.getXValuesAtEnd();               // (theta, phi, vtheta, vphi)
		Real speedSq = end[2] * end[2] + std::sin(end[0]) * std::sin(end[0]) * end[3] * end[3];

		out << "  [22.2] after one circumference (lambda = 2 pi): theta = " << end[0]
		    << " (start pi/2 = " << Constants::PI / 2 << "), phi = " << end[1]
		    << " (start + 2 pi = " << 2 * Constants::PI << ")\n"
		    << "  [22.2] |v|^2 conserved: " << speedSq << " (exact 1)\n";
	}

	// 22.3 FINALE - holonomy: parallel transport around a latitude circle
	void Cookbook_Recipe22_3_Holonomy(std::ostream& out)
	{
		Surfaces::Sphere sphere(1.0);
		DifferentialGeometry::InducedMetric2D metric(sphere);

		const Real theta0 = Constants::PI / 3;          // 60 deg colatitude circle
		ParametricCurveFromStdFunc<2> latitude(std::function<VectorN<Real, 2>(Real)>(
			[theta0](Real t) { return VectorN<Real, 2>{ theta0, t }; }));

		VectorN<Real, 2> V0{ 1.0, 0.0 };                // unit vector pointing south
		auto sol = IntegrateParallelTransportFixedStep<2>(metric, latitude, V0, 0.0, 2 * Constants::PI, 2000);

		auto Vend = sol.getXValuesAtEnd();

		// angle between start and transported vector, measured BY THE METRIC
		VectorN<Real, 2> pos{ theta0, 0.0 };
		Tensor2<2> g_ij = metric(pos);
		Real dot = 0.0, n0 = 0.0, n1 = 0.0;
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++) {
				dot += g_ij(i, j) * V0[i] * Vend[j];
				n0  += g_ij(i, j) * V0[i] * V0[j];
				n1  += g_ij(i, j) * Vend[i] * Vend[j];
			}
		Real angle = std::acos(dot / std::sqrt(n0 * n1));

		out << "  [22.3] transported around the 60 deg circle: V = (" << Vend[0] << ", " << Vend[1]
		    << "), started as (1, 0)\n"
		    << "  [22.3] holonomy angle = " << angle << " (Gauss-Bonnet: K x cap area = 2 pi (1 - cos 60) = pi = "
		    << Constants::PI << ")\n";
	}

	void Cookbook_Recipe22_DifferentialGeometry(std::ostream& out)
	{
		Cookbook_Recipe22_1_InducedMetric(out);
		Cookbook_Recipe22_2_Geodesics(out);
		Cookbook_Recipe22_3_Holonomy(out);
	}
} // namespace

void Docs_Demo_Cookbook()
{
	std::ostream& out = std::cout;
	out << "=== Docs_Demo_Cookbook (docs/COOKBOOK.md) ===\n";
	Cookbook_Recipe01_SolvingLinearSystems(out);
	Cookbook_Recipe02_RootFinding(out);
	Cookbook_Recipe03_Integrals(out);
	Cookbook_Recipe04_Derivatives(out);
	Cookbook_Recipe05_Optimization(out);
	Cookbook_Recipe06_DifferentialEquations(out);
	Cookbook_Recipe07_EigenvaluesAndEigenvectors(out);
	Cookbook_Recipe08_MatrixDecompositions(out);
	Cookbook_Recipe09_MatrixPropertiesAndDiagnostics(out);
	Cookbook_Recipe10_Statistics(out);
	Cookbook_Recipe11_Geometry3D(out);
	Cookbook_Recipe12_Geometry3DAlgorithms(out);
	Cookbook_Recipe13_InterpolatingFunctions(out);
	Cookbook_Recipe14_CurveFitting(out);
	Cookbook_Recipe15_CurvesAndSurfaces(out);
	Cookbook_Recipe16_Polynoms(out);
	Cookbook_Recipe17_Quaternions(out);
	Cookbook_Recipe18_AnalyzingFunctions(out);
	Cookbook_Recipe19_CoordinateTransformations(out);
	Cookbook_Recipe20_Fields(out);
	Cookbook_Recipe21_Tensors(out);
	Cookbook_Recipe22_DifferentialGeometry(out);
	Cookbook_Recipe23_FourierAlgorithms(out);
	Cookbook_Recipe24_GraphAlgorithms(out);
	Cookbook_Recipe25_DAESolvers(out);
	Cookbook_Recipe26_Algebra(out);
	Cookbook_Recipe27_ComplexAnalysis(out);
	Cookbook_Recipe30_CombinatoricsAndNumberTheory(out);
	out << "=== Cookbook demo done ===\n";
}
