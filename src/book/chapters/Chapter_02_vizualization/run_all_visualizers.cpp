#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "mml/base/Vector/VectorN.h"

#include "mml/base/Function.h"
#include "mml/base/InterpolatedFunction.h"

#include "mml/core/FunctionHelpers.h"
#include "mml/core/Curves.h"

#include "mml/base/ODESystem.h"
#include "mml/algorithms/ODESolvers/ODEAdaptiveIntegrator.h"

#include "mml/tools/Serializer.h"
#include "mml/tools/Visualizer.h"
#endif

using namespace MML;
using namespace MML::Curves;

// Forward declarations - defined in separate files
void Chapter02_RealFunctionVisualization();
void Chapter02_ParametricCurveVisualization();
void Chapter02_ScalarFunctionVisualization();
void Chapter02_VectorFieldVisualization();
void Chapter02_ParametricSurfaceVisualization();
void Chapter02_ParticleVisualization();

void Example3_WPF_Scalar_function_2D_visualization() {
	// Monkey saddle surface
	ScalarFunction<2> testFunc1{[](const VectorN<Real, 2>& x) { return 3 * x[0] / 5 * (x[0] * x[0] / 25 - 3 * x[1] * x[1] / 25); }};

	Visualizer::VisualizeScalarFunc2DCartesian(testFunc1, "Monkey saddle", -10.0, 10.0, 20, -10.0, 10.0, 20, "example3_wpf_surface1.mml");

	ScalarFunction<2> testFunc2{
		[](const VectorN<Real, 2>& x) { return (std::abs(x[0]) - 10) * (std::abs(x[1]) - 10) * sin(x[0]) * cos(x[1]); }};

	Visualizer::VisualizeScalarFunc2DCartesian(testFunc2, "Example surface", -10.0, 10.0, 50, -10.0, 10.0, 50, "example3_wpf_surface2.mml");
}

void Example3_WPF_Vector_field_2D_visualization() {
	VectorFunction<2> rotational_field{[](const VectorN<Real, 2>& x) {
		Real norm = x.NormL2();
		return VectorN<Real, 2>{-x[1] / norm, x[0] / norm};
	}};

	Visualizer::VisualizeVectorField2DCartesian(rotational_field, "rotational_field", -10.0, 10.0, 20, -10.0, 10.0, 20,
																							"example3_wpf_vector_field_2d.mml");
}

void Example3_WPF_Vector_field_3D_visualization() {
	// let's define potential gravity field of two masses as ScalarFunction
	ScalarFunction<3> gravity_field{[](const VectorN<Real, 3>& x) {
		const VectorN<Real, 3> x1{100.0, 0.0, 0.0};
		const VectorN<Real, 3> x2{-100.0, 50.0, 0.0};
		const Real m1 = 100.0;
		const Real m2 = -20.0;
		const Real G = 10.0;
		return -G * m1 / (x - x1).NormL2() - G * m2 / (x - x2).NormL2();
	}};

	// defining force field of two masses as VectorFunction
	VectorFunction<3> gravity_force_field{[](const VectorN<Real, 3>& x) {
		const VectorN<Real, 3> x1{100, 0, 0}, x2{-100, 0, 0};
		const Real m1 = 100.0, m2 = 100.0;
		const Real G = 10;
		return -G * m1 * (x - x1) / std::pow((x - x1).NormL2(), 3) - G * m2 * (x - x2) / std::pow((x - x2).NormL2(), 3);
	}};

	Visualizer::VisualizeVectorField3DCartesian(gravity_force_field, "Gravity field", -200, 200, 30, -200, 200, 30, -200, 200, 30,
																							"example3_wpf_vector_field_3d.mml");
}

void Chapter02_RunAllVisualizers() {
	// Select visualization backend - uncomment the one you want to use
	// MML::SetVisualizerBackend(MML::VisualizerBackend::FLTK);
	// MML::SetVisualizerBackend(MML::VisualizerBackend::Qt);

	// Debug: Show which backend is active
	auto backend = MML::GetVisualizerBackend();
	std::cout << "Active backend: ";
	switch (backend) {
	case MML::VisualizerBackend::WPF:
		std::cout << "WPF\n";
		break;
	case MML::VisualizerBackend::Qt:
		std::cout << "Qt\n";
		break;
	case MML::VisualizerBackend::FLTK:
		std::cout << "FLTK\n";
		break;
	case MML::VisualizerBackend::Auto:
		std::cout << "Auto\n";
		break;
	}
	std::cout << "Visualizer path: " << MML::GetRealFuncVisualizerPath() << "\n\n";

	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                  CHAPTER 02 - Visualization                   ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// Chapter02_RealFunctionVisualization();
	// Chapter02_ParametricCurveVisualization();
	Chapter02_ScalarFunctionVisualization();
	// Chapter02_VectorFieldVisualization();
	// Chapter02_ParametricSurfaceVisualization();
	// Chapter02_ParticleVisualization();
}
