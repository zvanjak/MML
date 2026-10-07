/**
 * @file show_typed_forms_3d.cpp
 * @brief Visualizes typed differential forms with MML_WorldVisualizer
 */

#ifdef MML_USE_SINGLE_HEADER
#include <MML.h>
#else
#include <mml/MMLBase.h>
#include <mml/base/DifferentialGeometry/DifferentialForm.h>
#include <mml/base/DifferentialGeometry/Hodge.h>
#include <mml/base/DifferentialGeometry/Metric.h>
#include <mml/tools/Visualizer.h>
#endif

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace MML;

namespace {

	std::string ScenePath()
	{
		std::filesystem::path results = GetResultFilesPath();
		std::filesystem::create_directories(results);
		return (results / "typed_forms_3d.mmlworld").string();
	}

	std::string ScenePath(const char* fileName)
	{
		std::filesystem::path results = GetResultFilesPath();
		std::filesystem::create_directories(results);
		return (results / fileName).string();
	}

	Real Norm(const TangentVector<3, Cartesian3>& vector)
	{
		return std::sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
	}

	TangentVector<3, Cartesian3> Unit(const TangentVector<3, Cartesian3>& vector)
	{
		Real norm = Norm(vector);
		return TangentVector<3, Cartesian3>{ vector[0] / norm, vector[1] / norm, vector[2] / norm };
	}

	void WriteVector(std::ostream& out, double px, double py, double pz,
					 TangentVector<3, Cartesian3> vector, double scale, double radius, const char* color)
	{
		out << "VECTOR " << px << " " << py << " " << pz << " "
			<< vector[0] * scale << " " << vector[1] * scale << " " << vector[2] * scale << " "
			<< radius << " " << color << "\n";
	}

	void WriteLine(std::ostream& out, double x0, double y0, double z0,
				   double x1, double y1, double z1, double radius, const char* color)
	{
		out << "LINE " << x0 << " " << y0 << " " << z0 << " "
			<< x1 << " " << y1 << " " << z1 << " " << radius << " " << color << "\n";
	}

	void WritePlane(std::ostream& out, const TangentVector<3, Cartesian3>& center,
					const TangentVector<3, Cartesian3>& normal, double size, const char* color, double opacity)
	{
		out << "PLANE "
			<< center[0] << " " << center[1] << " " << center[2] << " "
			<< normal[0] << " " << normal[1] << " " << normal[2] << " "
			<< size << " " << color << " " << opacity << "\n";
	}

	void WritePatch(std::ostream& out, const TangentVector<3, Cartesian3>& center,
					const TangentVector<3, Cartesian3>& u, const TangentVector<3, Cartesian3>& v,
					const char* color, double opacity)
	{
		out << "PATCH "
			<< center[0] << " " << center[1] << " " << center[2] << " "
			<< u[0] << " " << u[1] << " " << u[2] << " "
			<< v[0] << " " << v[1] << " " << v[2] << " "
			<< color << " " << opacity << "\n";
	}

	TangentVector<3, Cartesian3> VortexTubeField(const TangentVector<3, Cartesian3>& point)
	{
		Real x = point[0];
		Real y = point[1];
		Real z = point[2];
		Real r2 = x * x + y * y;
		Real damping = std::exp(Real{-0.00016} * r2);
		Real swirl = Real{0.018} * damping;
		Real upward = Real{0.85} + Real{0.15} * std::cos(Real{0.03} * z);
		return TangentVector<3, Cartesian3>{ -swirl * y, swirl * x, upward };
	}

	TangentVector<3, Cartesian3> CrossEuclidean(const TangentVector<3, Cartesian3>& a,
										 const TangentVector<3, Cartesian3>& b)
	{
		return TangentVector<3, Cartesian3>{
			a[1] * b[2] - a[2] * b[1],
			a[2] * b[0] - a[0] * b[2],
			a[0] * b[1] - a[1] * b[0]
		};
	}

	TangentVector<3, Cartesian3> MagneticDipoleField(const TangentVector<3, Cartesian3>& point)
	{
		Real x = point[0] / Real{55};
		Real y = point[1] / Real{55};
		Real z = point[2] / Real{55};
		Real r2 = x * x + y * y + z * z + Real{0.08};
		Real r = std::sqrt(r2);
		Real r5 = r2 * r2 * r;
		Real factor = Real{1.6} / r5;

		return TangentVector<3, Cartesian3>{
			Real{3} * z * x * factor,
			Real{3} * z * y * factor,
			(Real{3} * z * z - r2) * factor
		};
	}

	void WriteCircle(std::ostream& scene, double radius, double z, const char* color)
	{
		constexpr int numSegments = 40;
		for (int i = 0; i < numSegments; i++) {
			double a0 = 2.0 * 3.14159265358979323846 * i / numSegments;
			double a1 = 2.0 * 3.14159265358979323846 * (i + 1) / numSegments;
			WriteLine(scene,
				radius * std::cos(a0), radius * std::sin(a0), z,
				radius * std::cos(a1), radius * std::sin(a1), z,
				0.75, color);
		}
	}

	void WriteLocalFormGlyph(std::ostream& scene, const TangentVector<3, Cartesian3>& point,
						 TangentVector<3, Cartesian3> vector, Real magnitudeScale,
						 const char* arrowColor = "#FFD23F", const char* planeColor = "#56D6C9",
						 const char* positivePatchColor = "#FF6B35", const char* negativePatchColor = "#7D5CFF")
	{
		Metric<3, Cartesian3> metric = Metric<3, Cartesian3>::Euclidean();
		Orientation orientation = Orientation::Positive;
		Form1<3, Cartesian3> flat = ToForm(Flat(metric, vector));
		Form2<3, Cartesian3> fluxForm = hodge_star(flat, metric, orientation);

		Real magnitude = Norm(vector);
		if (magnitude < 1e-10)
			return;

		TangentVector<3, Cartesian3> normal = Unit(vector);
		TangentVector<3, Cartesian3> reference{ REAL(0.0), REAL(0.0), REAL(1.0) };
		if (std::abs(normal[2]) > 0.88)
			reference = TangentVector<3, Cartesian3>{ REAL(1.0), REAL(0.0), REAL(0.0) };
		TangentVector<3, Cartesian3> tangentA = Unit(CrossEuclidean(normal, reference));
		TangentVector<3, Cartesian3> tangentB = Unit(CrossEuclidean(normal, tangentA));

		Real arrowScale = Real{42} * magnitudeScale;
		Real planeSpacing = Real{6};
		Real planeSize = Real{34} + Real{8} * magnitude;
		Real patchSize = Real{20} + Real{7} * magnitude;

		WriteVector(scene, point[0], point[1], point[2], vector, arrowScale, 2.2, arrowColor);
		for (int i = -1; i <= 1; i++) {
			TangentVector<3, Cartesian3> center{
				point[0] + normal[0] * planeSpacing * i,
				point[1] + normal[1] * planeSpacing * i,
				point[2] + normal[2] * planeSpacing * i
			};
			WritePlane(scene, center, normal, planeSize, planeColor, 0.34);
		}

		TangentVector<3, Cartesian3> patchU{ tangentA[0] * patchSize, tangentA[1] * patchSize, tangentA[2] * patchSize };
		TangentVector<3, Cartesian3> patchV{ tangentB[0] * patchSize, tangentB[1] * patchSize, tangentB[2] * patchSize };
		Real flux = fluxForm(patchU, patchV);
		WritePatch(scene, point, patchU, patchV, flux >= 0 ? positivePatchColor : negativePatchColor, 0.50);
		scene << "POINT " << point[0] << " " << point[1] << " " << point[2] << " 3.4 #FFFFFF\n";
	}

} // anonymous namespace

void Show_Typed_Forms_3D_Visualization()
{
	std::cout << "\n=== Typed Differential Forms 3D Visualization ===\n\n";

	Metric<3, Cartesian3> metric = Metric<3, Cartesian3>::Euclidean();
	Orientation orientation = Orientation::Positive;

	TangentVector<3, Cartesian3> velocity{ REAL(2.0), -REAL(1.0), REAL(3.0) };
	TangentVector<3, Cartesian3> ru{ REAL(2.4), REAL(0.0), REAL(0.0) };
	TangentVector<3, Cartesian3> rv{ REAL(0.0), REAL(1.7), REAL(0.0) };

	Form1<3, Cartesian3> velocityFlat = ToForm(Flat(metric, velocity));
	Form2<3, Cartesian3> fluxForm = hodge_star(velocityFlat, metric, orientation);
	TangentVector<3, Cartesian3> normal = cross(ru, rv, metric, orientation);

	Real fluxThroughPatch = fluxForm(ru, rv);
	TangentVector<3, Cartesian3> oneFormNormal = Unit(velocity);

	std::string path = ScenePath();
	std::ofstream scene(path);
	scene << std::fixed << std::setprecision(6);
	scene << "MML_WORLD_SCENE 1\n";
	scene << "TITLE Typed forms: vector, one-form planes, flux two-form\n";
	scene << "CAMERA 420 260 280\n";
	scene << "# VECTOR px py pz vx vy vz radius color\n";
	scene << "# PLANE cx cy cz nx ny nz size color opacity\n";
	scene << "# PATCH cx cy cz ux uy uz vx vy vz color opacity\n";
	scene << "# LINE x0 y0 z0 x1 y1 z1 radius color\n\n";

	WriteLine(scene, -160, 0, 0, 160, 0, 0, 1.2, "#777777");
	WriteLine(scene, 0, -160, 0, 0, 160, 0, 1.2, "#777777");
	WriteLine(scene, 0, 0, -80, 0, 0, 180, 1.2, "#777777");

	WriteVector(scene, 0, 0, 0, velocity, 35.0, 4.0, "#FFD23F");
	WriteVector(scene, 0, 0, 0, ru, 45.0, 2.5, "#4CB3FF");
	WriteVector(scene, 0, 0, 0, rv, 45.0, 2.5, "#4CB3FF");
	WriteVector(scene, 0, 0, 0, Unit(normal), 95.0, 3.0, "#27D17F");

	for (int i = -3; i <= 3; i++) {
		double offset = i * 22.0;
		scene << "PLANE "
			<< oneFormNormal[0] * offset << " " << oneFormNormal[1] * offset << " " << oneFormNormal[2] * offset << " "
			<< oneFormNormal[0] << " " << oneFormNormal[1] << " " << oneFormNormal[2] << " "
			<< 145.0 << " #56D6C9 0.22\n";
	}

	scene << "PATCH 0 0 0 "
		  << ru[0] * 45.0 << " " << ru[1] * 45.0 << " " << ru[2] * 45.0 << " "
		  << rv[0] * 45.0 << " " << rv[1] * 45.0 << " " << rv[2] * 45.0 << " "
		  << "#FF6B35 0.48\n";
	scene << "POINT 0 0 0 4.0 #FFFFFF\n";
	scene.close();

	std::cout << "Generated world scene: " << path << "\n";
	std::cout << "velocity^flat components: (" << velocityFlat[0] << ", " << velocityFlat[1] << ", " << velocityFlat[2] << ")\n";
	std::cout << "flux_form(ru, rv): " << fluxThroughPatch << "\n";
	std::cout << "\nLegend:\n";
	std::cout << "  yellow arrow  : tangent vector v\n";
	std::cout << "  cyan planes   : level planes of the one-form v^flat\n";
	std::cout << "  orange patch  : oriented surface for the two-form flux\n";
	std::cout << "  green arrow   : oriented normal from cross(ru, rv)\n\n";

	if (std::getenv("MML_VISUALIZER_SAVE_ONLY") != nullptr) {
		std::cout << "MML_VISUALIZER_SAVE_ONLY is set; not launching MML_WorldVisualizer.\n";
		return;
	}

	auto result = Visualizer::VisualizeWorldSceneFromFile(path);
	if (!result.success)
		std::cerr << "World visualizer error: " << result.errorMessage << "\n";
}

void Show_Typed_Forms_Vortex_Tube_Visualization()
{
	std::cout << "\n=== Typed Differential Forms: Vortex Tube Field ===\n\n";

	std::string path = ScenePath("typed_forms_vortex_tube.mmlworld");
	std::ofstream scene(path);
	scene << std::fixed << std::setprecision(6);
	scene << "MML_WORLD_SCENE 1\n";
	scene << "TITLE Typed forms sampled on a damped vortex tube field\n";
	scene << "CAMERA 460 340 310\n";
	scene << "# v(x,y,z) = exp(-0.00016(x^2+y^2))*(-0.018y, 0.018x, 0) + upward z-flow\n\n";

	WriteLine(scene, -150, 0, 0, 150, 0, 0, 1.0, "#666666");
	WriteLine(scene, 0, -150, 0, 0, 150, 0, 1.0, "#666666");
	WriteLine(scene, 0, 0, -80, 0, 0, 160, 1.0, "#666666");
	for (int i = -2; i <= 2; i++) {
		WriteLine(scene, -120, i * 40, -60, 120, i * 40, -60, 0.55, "#343434");
		WriteLine(scene, i * 40, -120, -60, i * 40, 120, -60, 0.55, "#343434");
	}

	Real maxMagnitude = Real{0};
	for (Real z : { Real{-45}, Real{0}, Real{45} }) {
		for (int i = 0; i < 8; i++) {
			Real angle = i * Constants::PI / Real{4} + (z > Real{0} ? Real{0.18} : (z < Real{0} ? Real{-0.18} : Real{0}));
			TangentVector<3, Cartesian3> point{ Real{74} * std::cos(angle), Real{74} * std::sin(angle), z };
			maxMagnitude = std::max(maxMagnitude, Norm(VortexTubeField(point)));
		}
	}

	for (Real z : { Real{-45}, Real{0}, Real{45} }) {
		for (int i = 0; i < 8; i++) {
			Real angle = i * Constants::PI / Real{4} + (z > Real{0} ? Real{0.18} : (z < Real{0} ? Real{-0.18} : Real{0}));
			TangentVector<3, Cartesian3> point{ Real{74} * std::cos(angle), Real{74} * std::sin(angle), z };
			TangentVector<3, Cartesian3> vector = VortexTubeField(point);
			WriteLocalFormGlyph(scene, point, vector, 1.0 / maxMagnitude);
		}
	}

	scene.close();

	std::cout << "Generated world scene: " << path << "\n";
	std::cout << "Sampled 24 points on a damped helical vortex tube field.\n";
	std::cout << "\nLegend:\n";
	std::cout << "  yellow arrows : field vectors v(p)\n";
	std::cout << "  cyan planes   : local one-form v_flat(p)\n";
	std::cout << "  orange patches: local flux two-form *v_flat(p)\n";
	std::cout << "  white points  : sample locations\n\n";

	if (std::getenv("MML_VISUALIZER_SAVE_ONLY") != nullptr) {
		std::cout << "MML_VISUALIZER_SAVE_ONLY is set; not launching MML_WorldVisualizer.\n";
		return;
	}

	auto result = Visualizer::VisualizeWorldSceneFromFile(path);
	if (!result.success)
		std::cerr << "World visualizer error: " << result.errorMessage << "\n";
}

void Show_Typed_Forms_EM_Dipole_Visualization()
{
	std::cout << "\n=== Typed Differential Forms: Electromagnetic Dipole Field ===\n\n";

	std::string path = ScenePath("typed_forms_em_dipole.mmlworld");
	std::ofstream scene(path);
	scene << std::fixed << std::setprecision(6);
	scene << "MML_WORLD_SCENE 1\n";
	scene << "TITLE Typed forms sampled on a magnetic dipole-style field\n";
	scene << "CAMERA 500 360 330\n";
	scene << "# B(r) = (3(m.r)r - m|r|^2)/|r|^5, regularized near the origin, m along +z\n\n";

	WriteLine(scene, -155, 0, 0, 155, 0, 0, 1.0, "#666666");
	WriteLine(scene, 0, -155, 0, 0, 155, 0, 1.0, "#666666");
	WriteLine(scene, 0, 0, -135, 0, 0, 150, 1.6, "#C7C7C7");
	WriteCircle(scene, 32.0, 0.0, "#D98A1F");
	WriteCircle(scene, 42.0, 0.0, "#D98A1F");
	scene << "POINT 0 0 42 5.0 #FF4D4D\n";
	scene << "POINT 0 0 -42 5.0 #4CB3FF\n";

	std::vector<TangentVector<3, Cartesian3>> samplePoints;
	for (Real z : { Real{-85}, Real{-40}, Real{40}, Real{85} }) {
		for (int i = 0; i < 8; i++) {
			Real angle = i * Constants::PI / Real{4} + (z > Real{0} ? Real{0.22} : Real{-0.22});
			Real radius = std::abs(z) > Real{60} ? Real{58} : Real{86};
			samplePoints.push_back(TangentVector<3, Cartesian3>{
				radius * std::cos(angle), radius * std::sin(angle), z
			});
		}
	}
	for (Real x : { Real{-105}, Real{-70}, Real{70}, Real{105} }) {
		samplePoints.push_back(TangentVector<3, Cartesian3>{ x, 0.0, 0.0 });
		samplePoints.push_back(TangentVector<3, Cartesian3>{ 0.0, x, 0.0 });
	}

	Real maxMagnitude = Real{0};
	for (const auto& point : samplePoints)
		maxMagnitude = std::max(maxMagnitude, Norm(MagneticDipoleField(point)));

	for (const auto& point : samplePoints) {
		TangentVector<3, Cartesian3> field = MagneticDipoleField(point);
		WriteLocalFormGlyph(scene, point, field, 1.0 / maxMagnitude,
			"#4CB3FF", "#9CE8FF", "#FFB000", "#7D5CFF");
	}

	scene.close();

	std::cout << "Generated world scene: " << path << "\n";
	std::cout << "Sampled " << samplePoints.size() << " points on a regularized magnetic dipole-style field.\n";
	std::cout << "\nLegend:\n";
	std::cout << "  blue arrows  : magnetic field vectors B(p)\n";
	std::cout << "  pale planes  : local one-form B_flat(p)\n";
	std::cout << "  amber patches: local magnetic flux two-form *B_flat(p)\n";
	std::cout << "  orange rings : stylized current loop / dipole source\n\n";

	if (std::getenv("MML_VISUALIZER_SAVE_ONLY") != nullptr) {
		std::cout << "MML_VISUALIZER_SAVE_ONLY is set; not launching MML_WorldVisualizer.\n";
		return;
	}

	auto result = Visualizer::VisualizeWorldSceneFromFile(path);
	if (!result.success)
		std::cerr << "World visualizer error: " << result.errorMessage << "\n";
}