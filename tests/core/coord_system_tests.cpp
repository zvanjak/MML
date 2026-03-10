///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        coord_system_tests.cpp                                              ///
///  Description: Tests for CoordSystem.h - Reference frame hierarchy and transformations
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#include <catch2/catch_all.hpp>

#include "MMLBase.h"
#include "core/CoordSystem.h"

using namespace MML;

namespace MML::Tests::Core::CoordSystemTests {

///////////////////////////////////////////////////////////////////////////////////////////
///                           ReferenceFrame3D - Base Class                             ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("ReferenceFrame3D - Default construction", "[CoordSystem][ReferenceFrame]")
{
	ReferenceFrame3D frame;
	
	REQUIRE(frame.GetParentFrame() == nullptr);
	REQUIRE(frame.GetChildFrames().empty());
}

TEST_CASE("ReferenceFrame3D - Default origin and position", "[CoordSystem][ReferenceFrame]")
{
	ReferenceFrame3D frame;
	
	// Default frame returns zero position
	Vector3Cartesian origin = frame.GetOriginPositionAtTime(0.0);
	REQUIRE(origin[0] == Catch::Approx(0.0));
	REQUIRE(origin[1] == Catch::Approx(0.0));
	REQUIRE(origin[2] == Catch::Approx(0.0));
	
	// Default frame returns position unchanged
	VectorN<Real, 3> localPos({1.0, 2.0, 3.0});
	Vector3Cartesian result = frame.GetLocalPosInParentFrameAtTime(localPos, 1.0);
	REQUIRE(result[0] == Catch::Approx(1.0));
	REQUIRE(result[1] == Catch::Approx(2.0));
	REQUIRE(result[2] == Catch::Approx(3.0));
}

TEST_CASE("ReferenceFrame3D - Parent-child relationship", "[CoordSystem][ReferenceFrame]")
{
	ReferenceFrame3D parent;
	ReferenceFrame3D child;
	
	REQUIRE(child.SetParentFrame(&parent));
	
	REQUIRE(child.GetParentFrame() == &parent);
	REQUIRE(parent.GetChildFrames().size() == 1);
	REQUIRE(parent.GetChildFrames()[0] == &child);
}

TEST_CASE("ReferenceFrame3D - Multiple children", "[CoordSystem][ReferenceFrame]")
{
	ReferenceFrame3D parent;
	ReferenceFrame3D child1, child2, child3;
	
	REQUIRE(parent.AddChildFrame(&child1));
	REQUIRE(parent.AddChildFrame(&child2));
	REQUIRE(parent.AddChildFrame(&child3));
	
	REQUIRE(parent.GetChildFrames().size() == 3);
	REQUIRE(child1.GetParentFrame() == &parent);
	REQUIRE(child2.GetParentFrame() == &parent);
	REQUIRE(child3.GetParentFrame() == &parent);
}

TEST_CASE("ReferenceFrame3D - Remove child", "[CoordSystem][ReferenceFrame]")
{
	ReferenceFrame3D parent;
	ReferenceFrame3D child1, child2;
	
	parent.AddChildFrame(&child1);
	parent.AddChildFrame(&child2);
	
	REQUIRE(parent.RemoveChildFrame(&child1));
	
	REQUIRE(parent.GetChildFrames().size() == 1);
	REQUIRE(child1.GetParentFrame() == nullptr);
	REQUIRE(child2.GetParentFrame() == &parent);
}

TEST_CASE("ReferenceFrame3D - Reparenting", "[CoordSystem][ReferenceFrame]")
{
	ReferenceFrame3D parent1, parent2;
	ReferenceFrame3D child;
	
	child.SetParentFrame(&parent1);
	REQUIRE(parent1.GetChildFrames().size() == 1);
	
	// Reparent to new parent
	child.SetParentFrame(&parent2);
	
	REQUIRE(parent1.GetChildFrames().empty());
	REQUIRE(parent2.GetChildFrames().size() == 1);
	REQUIRE(child.GetParentFrame() == &parent2);
}

TEST_CASE("ReferenceFrame3D - Cycle prevention (direct)", "[CoordSystem][ReferenceFrame]")
{
	ReferenceFrame3D parent;
	ReferenceFrame3D child;
	
	child.SetParentFrame(&parent);
	
	// Parent cannot become child of its child (would create cycle)
	REQUIRE_FALSE(parent.SetParentFrame(&child));
	
	// Original relationship unchanged
	REQUIRE(child.GetParentFrame() == &parent);
	REQUIRE(parent.GetParentFrame() == nullptr);
}

TEST_CASE("ReferenceFrame3D - Cycle prevention (indirect)", "[CoordSystem][ReferenceFrame]")
{
	ReferenceFrame3D grandparent;
	ReferenceFrame3D parent;
	ReferenceFrame3D child;
	
	parent.SetParentFrame(&grandparent);
	child.SetParentFrame(&parent);
	
	// Grandparent cannot become child of grandchild
	REQUIRE_FALSE(grandparent.SetParentFrame(&child));
	
	// Chain unchanged
	REQUIRE(child.GetParentFrame() == &parent);
	REQUIRE(parent.GetParentFrame() == &grandparent);
}

TEST_CASE("ReferenceFrame3D - Add null child", "[CoordSystem][ReferenceFrame][validation]")
{
	ReferenceFrame3D parent;
	
	REQUIRE_FALSE(parent.AddChildFrame(nullptr));
	REQUIRE(parent.GetChildFrames().empty());
}

TEST_CASE("ReferenceFrame3D - Add duplicate child", "[CoordSystem][ReferenceFrame]")
{
	ReferenceFrame3D parent;
	ReferenceFrame3D child;
	
	REQUIRE(parent.AddChildFrame(&child));
	REQUIRE(parent.AddChildFrame(&child));  // Should be no-op, return true
	
	REQUIRE(parent.GetChildFrames().size() == 1);
}

TEST_CASE("ReferenceFrame3D - Remove non-existent child", "[CoordSystem][ReferenceFrame]")
{
	ReferenceFrame3D parent;
	ReferenceFrame3D notAChild;
	
	REQUIRE_FALSE(parent.RemoveChildFrame(&notAChild));
}

TEST_CASE("ReferenceFrame3D - Set same parent (no-op)", "[CoordSystem][ReferenceFrame]")
{
	ReferenceFrame3D parent;
	ReferenceFrame3D child;
	
	child.SetParentFrame(&parent);
	REQUIRE(child.SetParentFrame(&parent));  // Should be no-op
	
	REQUIRE(parent.GetChildFrames().size() == 1);
}

TEST_CASE("ReferenceFrame3D - Clear parent", "[CoordSystem][ReferenceFrame]")
{
	ReferenceFrame3D parent;
	ReferenceFrame3D child;
	
	child.SetParentFrame(&parent);
	REQUIRE(child.SetParentFrame(nullptr));
	
	REQUIRE(child.GetParentFrame() == nullptr);
	REQUIRE(parent.GetChildFrames().empty());
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           InertialFrame3D                                           ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("InertialFrame3D - Default construction", "[CoordSystem][InertialFrame]")
{
	InertialFrame3D frame;
	
	REQUIRE(frame.GetParentFrame() == nullptr);
	
	// Default is stationary at origin
	Vector3Cartesian pos = frame.GetOriginPositionAtTime(0.0);
	REQUIRE(pos[0] == Catch::Approx(0.0));
	REQUIRE(pos[1] == Catch::Approx(0.0));
	REQUIRE(pos[2] == Catch::Approx(0.0));
}

TEST_CASE("InertialFrame3D - Constant velocity motion", "[CoordSystem][InertialFrame]")
{
	ReferenceFrame3D parent;
	Vector3Cartesian velocity({1.0, 2.0, 0.0});
	Vector3Cartesian initial_pos({0.0, 0.0, 0.0});
	
	InertialFrame3D frame(&parent, velocity, initial_pos);
	
	// Position at t=0
	Vector3Cartesian pos0 = frame.GetOriginPositionAtTime(0.0);
	REQUIRE(pos0[0] == Catch::Approx(0.0));
	REQUIRE(pos0[1] == Catch::Approx(0.0));
	
	// Position at t=1
	Vector3Cartesian pos1 = frame.GetOriginPositionAtTime(1.0);
	REQUIRE(pos1[0] == Catch::Approx(1.0));
	REQUIRE(pos1[1] == Catch::Approx(2.0));
	
	// Position at t=5
	Vector3Cartesian pos5 = frame.GetOriginPositionAtTime(5.0);
	REQUIRE(pos5[0] == Catch::Approx(5.0));
	REQUIRE(pos5[1] == Catch::Approx(10.0));
}

TEST_CASE("InertialFrame3D - With initial offset", "[CoordSystem][InertialFrame]")
{
	ReferenceFrame3D parent;
	Vector3Cartesian velocity({1.0, 0.0, 0.0});
	Vector3Cartesian initial_pos({10.0, 20.0, 30.0});
	
	InertialFrame3D frame(&parent, velocity, initial_pos);
	
	// Position at t=0 equals initial position
	Vector3Cartesian pos0 = frame.GetOriginPositionAtTime(0.0);
	REQUIRE(pos0[0] == Catch::Approx(10.0));
	REQUIRE(pos0[1] == Catch::Approx(20.0));
	REQUIRE(pos0[2] == Catch::Approx(30.0));
	
	// Position at t=5
	Vector3Cartesian pos5 = frame.GetOriginPositionAtTime(5.0);
	REQUIRE(pos5[0] == Catch::Approx(15.0));  // 10 + 5*1
	REQUIRE(pos5[1] == Catch::Approx(20.0));
	REQUIRE(pos5[2] == Catch::Approx(30.0));
}

TEST_CASE("InertialFrame3D - Local to parent transformation", "[CoordSystem][InertialFrame]")
{
	ReferenceFrame3D parent;
	Vector3Cartesian velocity({1.0, 0.0, 0.0});
	Vector3Cartesian initial_pos({0.0, 0.0, 0.0});
	
	InertialFrame3D frame(&parent, velocity, initial_pos);
	
	VectorN<Real, 3> localPos({5.0, 5.0, 0.0});
	
	// At t=0: frame at origin, so parent pos = local pos
	Vector3Cartesian result0 = frame.GetLocalPosInParentFrameAtTime(localPos, 0.0);
	REQUIRE(result0[0] == Catch::Approx(5.0));
	REQUIRE(result0[1] == Catch::Approx(5.0));
	
	// At t=10: frame at (10,0,0), so parent pos = (15,5,0)
	Vector3Cartesian result10 = frame.GetLocalPosInParentFrameAtTime(localPos, 10.0);
	REQUIRE(result10[0] == Catch::Approx(15.0));
	REQUIRE(result10[1] == Catch::Approx(5.0));
}

TEST_CASE("InertialFrame3D - Negative time", "[CoordSystem][InertialFrame]")
{
	ReferenceFrame3D parent;
	Vector3Cartesian velocity({2.0, 0.0, 0.0});
	Vector3Cartesian initial_pos({0.0, 0.0, 0.0});
	
	InertialFrame3D frame(&parent, velocity, initial_pos);
	
	// Negative time works (extrapolate backwards)
	Vector3Cartesian pos_neg = frame.GetOriginPositionAtTime(-5.0);
	REQUIRE(pos_neg[0] == Catch::Approx(-10.0));
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           CircleOrbitingFrame3DCartesian                            ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("CircleOrbitingFrame3DCartesian - Circular orbit position", "[CoordSystem][CircleOrbitingFrame]")
{
	ReferenceFrame3D parent;
	Real radius = 10.0;
	Real period = 2.0 * Constants::PI;  // Period = 2π seconds, so ω = 1 rad/s
	
	CircleOrbitingFrame3DCartesian orbit(&parent, radius, period);
	
	// At t=0: position at (radius, 0, 0)
	Vector3Cartesian pos0 = orbit.GetOriginPositionAtTime(0.0);
	REQUIRE(pos0[0] == Catch::Approx(10.0));
	REQUIRE(pos0[1] == Catch::Approx(0.0));
	REQUIRE(pos0[2] == Catch::Approx(0.0));
	
	// At t=π/2: quarter orbit, at (0, radius, 0)
	Vector3Cartesian pos_quarter = orbit.GetOriginPositionAtTime(Constants::PI / 2);
	REQUIRE(pos_quarter[0] == Catch::Approx(0.0).margin(1e-10));
	REQUIRE(pos_quarter[1] == Catch::Approx(10.0));
	
	// At t=π: half orbit, at (-radius, 0, 0)
	Vector3Cartesian pos_half = orbit.GetOriginPositionAtTime(Constants::PI);
	REQUIRE(pos_half[0] == Catch::Approx(-10.0));
	REQUIRE(pos_half[1] == Catch::Approx(0.0).margin(1e-10));
	
	// At t=2π: full orbit, back to start
	Vector3Cartesian pos_full = orbit.GetOriginPositionAtTime(2.0 * Constants::PI);
	REQUIRE(pos_full[0] == Catch::Approx(10.0));
	REQUIRE(pos_full[1] == Catch::Approx(0.0).margin(1e-10));
}

TEST_CASE("CircleOrbitingFrame3DCartesian - With initial angle", "[CoordSystem][CircleOrbitingFrame]")
{
	ReferenceFrame3D parent;
	Real radius = 5.0;
	Real period = 2.0 * Constants::PI;
	Real initial_angle = Constants::PI / 2;  // Start at 90 degrees
	
	CircleOrbitingFrame3DCartesian orbit(&parent, radius, period, initial_angle);
	
	// At t=0: position at (0, radius, 0) due to 90° initial angle
	Vector3Cartesian pos0 = orbit.GetOriginPositionAtTime(0.0);
	REQUIRE(pos0[0] == Catch::Approx(0.0).margin(1e-10));
	REQUIRE(pos0[1] == Catch::Approx(5.0));
}

TEST_CASE("CircleOrbitingFrame3DCartesian - Local to parent transformation", "[CoordSystem][CircleOrbitingFrame]")
{
	ReferenceFrame3D parent;
	Real radius = 10.0;
	Real period = 2.0 * Constants::PI;
	
	CircleOrbitingFrame3DCartesian orbit(&parent, radius, period);
	
	VectorN<Real, 3> localPos({1.0, 0.0, 0.0});  // 1 unit in local x
	
	// At t=0: orbit at (10, 0, 0), so parent pos = (11, 0, 0)
	Vector3Cartesian result0 = orbit.GetLocalPosInParentFrameAtTime(localPos, 0.0);
	REQUIRE(result0[0] == Catch::Approx(11.0));
	REQUIRE(result0[1] == Catch::Approx(0.0));
	
	// At t=π/2: orbit at (0, 10, 0), so parent pos = (1, 10, 0)
	// Note: frame orientation stays fixed (not rotating with orbit)
	Vector3Cartesian result_quarter = orbit.GetLocalPosInParentFrameAtTime(localPos, Constants::PI / 2);
	REQUIRE(result_quarter[0] == Catch::Approx(1.0).margin(1e-10));
	REQUIRE(result_quarter[1] == Catch::Approx(10.0));
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           RotatingFrame3D                                           ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("RotatingFrame3D - Origin stays fixed", "[CoordSystem][RotatingFrame]")
{
	ReferenceFrame3D parent;
	Real period = 1.0;
	VectorN<Real, 3> axis({0, 0, 1});  // Z-axis rotation
	
	RotatingFrame3D frame(&parent, period, axis);
	
	// Origin doesn't move regardless of time
	for (Real t : {0.0, 0.5, 1.0, 10.0}) {
		Vector3Cartesian origin = frame.GetOriginPositionAtTime(t);
		REQUIRE(origin[0] == Catch::Approx(0.0));
		REQUIRE(origin[1] == Catch::Approx(0.0));
		REQUIRE(origin[2] == Catch::Approx(0.0));
	}
}

TEST_CASE("RotatingFrame3D - Basic construction", "[CoordSystem][RotatingFrame]")
{
	ReferenceFrame3D parent;
	Real period = 24.0;  // 24 hour rotation
	VectorN<Real, 3> axis({0, 0, 1});
	
	RotatingFrame3D frame(&parent, period, axis);
	
	REQUIRE(frame._period == Catch::Approx(24.0));
	REQUIRE(frame._axis[2] == Catch::Approx(1.0));
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           HardSphereRotatingFrame                                   ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("HardSphereRotatingFrame - Construction", "[CoordSystem][HardSphereRotating]")
{
	ReferenceFrame3D parent;
	Real radius = 6371.0;  // Earth radius in km
	Real period = 24.0;    // 24 hour rotation
	VectorN<Real, 3> axis({0, 0, 1});
	
	HardSphereRotatingFrame earth(&parent, radius, period, axis);
	
	REQUIRE(earth._radius == Catch::Approx(6371.0));
	REQUIRE(earth._period == Catch::Approx(24.0));
}

TEST_CASE("HardSphereRotatingFrame - Origin stays fixed", "[CoordSystem][HardSphereRotating]")
{
	ReferenceFrame3D parent;
	HardSphereRotatingFrame sphere(&parent, 100.0, 10.0, VectorN<Real, 3>({0, 0, 1}));
	
	for (Real t : {0.0, 5.0, 10.0}) {
		Vector3Cartesian origin = sphere.GetOriginPositionAtTime(t);
		REQUIRE(origin[0] == Catch::Approx(0.0));
		REQUIRE(origin[1] == Catch::Approx(0.0));
		REQUIRE(origin[2] == Catch::Approx(0.0));
	}
}

TEST_CASE("HardSphereRotatingFrame - Equator point transformation", "[CoordSystem][HardSphereRotating]")
{
	ReferenceFrame3D parent;
	Real radius = 1.0;
	Real period = 2.0 * Constants::PI;
	VectorN<Real, 3> axis({0, 0, 1});
	
	HardSphereRotatingFrame sphere(&parent, radius, period, axis);
	
	// Point on equator (lat=0, long=0, alt=0)
	VectorN<Real, 3> localPos({0.0, 0.0, 0.0});  // (lat, long, altitude)
	
	Vector3Cartesian result = sphere.GetLocalPosInParentFrameAtTime(localPos, 0.0);
	
	// At equator, lat=0, long=0: transforms to spherical (r=1, θ=90°, φ=0°)
	// Which in Cartesian is (1, 0, 0)
	REQUIRE(result[0] == Catch::Approx(1.0));
	REQUIRE(result[1] == Catch::Approx(0.0).margin(1e-10));
	REQUIRE(result[2] == Catch::Approx(0.0).margin(1e-10));
}

TEST_CASE("HardSphereRotatingFrame - North pole transformation", "[CoordSystem][HardSphereRotating]")
{
	ReferenceFrame3D parent;
	Real radius = 1.0;
	HardSphereRotatingFrame sphere(&parent, radius, 10.0, VectorN<Real, 3>({0, 0, 1}));
	
	// North pole (lat=90, any long, alt=0)
	VectorN<Real, 3> northPole({90.0, 0.0, 0.0});
	
	Vector3Cartesian result = sphere.GetLocalPosInParentFrameAtTime(northPole, 0.0);
	
	// North pole: spherical (r=1, θ=0°, φ=any) -> Cartesian (0, 0, 1)
	REQUIRE(result[0] == Catch::Approx(0.0).margin(1e-10));
	REQUIRE(result[1] == Catch::Approx(0.0).margin(1e-10));
	REQUIRE(result[2] == Catch::Approx(1.0));
}

TEST_CASE("HardSphereRotatingFrame - With altitude", "[CoordSystem][HardSphereRotating]")
{
	ReferenceFrame3D parent;
	Real radius = 10.0;
	HardSphereRotatingFrame sphere(&parent, radius, 10.0, VectorN<Real, 3>({0, 0, 1}));
	
	// Point at equator with altitude
	VectorN<Real, 3> elevated({0.0, 0.0, 5.0});  // 5 units above surface
	
	Vector3Cartesian result = sphere.GetLocalPosInParentFrameAtTime(elevated, 0.0);
	
	// Should be at (radius + altitude, 0, 0) = (15, 0, 0)
	REQUIRE(result[0] == Catch::Approx(15.0));
	REQUIRE(result[1] == Catch::Approx(0.0).margin(1e-10));
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           RotatingSphereLocalCartesian                              ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("RotatingSphereLocalCartesian - Construction", "[CoordSystem][LocalCartesian]")
{
	ReferenceFrame3D grandparent;
	HardSphereRotatingFrame earth(&grandparent, 6371.0, 24.0, VectorN<Real, 3>({0, 0, 1}));
	
	// Create local frame at latitude 45°, longitude 0°
	RotatingSphereLocalCartesian local(earth, 45.0, 0.0);
	
	// Should compile and construct without issues
	REQUIRE(true);
}

TEST_CASE("RotatingSphereLocalCartesian - Origin position", "[CoordSystem][LocalCartesian]")
{
	ReferenceFrame3D grandparent;
	HardSphereRotatingFrame earth(&grandparent, 6371.0, 24.0, VectorN<Real, 3>({0, 0, 1}));
	
	RotatingSphereLocalCartesian local(earth, 0.0, 0.0);
	
	// Origin position (currently returns zero - TODO in implementation)
	Vector3Cartesian origin = local.GetOriginPositionAtTime(0.0);
	REQUIRE(origin[0] == Catch::Approx(0.0));  // Current implementation
}

TEST_CASE("RotatingSphereLocalCartesian - Local to geographic", "[CoordSystem][LocalCartesian]")
{
	ReferenceFrame3D grandparent;
	HardSphereRotatingFrame earth(&grandparent, 6371.0, 24.0, VectorN<Real, 3>({0, 0, 1}));
	
	// Local frame at equator, prime meridian
	RotatingSphereLocalCartesian local(earth, 0.0, 0.0);
	
	// Origin point (0,0,0) should map back to (0°, 0°, 0)
	VectorN<Real, 3> origin({0.0, 0.0, 0.0});
	Vector3Cartesian result = local.GetLocalPosInParentFrameAtTime(origin, 0.0);
	
	REQUIRE(result[0] == Catch::Approx(0.0));  // latitude
	REQUIRE(result[1] == Catch::Approx(0.0));  // longitude
	REQUIRE(result[2] == Catch::Approx(0.0));  // altitude
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           Frame Hierarchy                                           ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("Frame hierarchy - Three level chain", "[CoordSystem][Hierarchy]")
{
	ReferenceFrame3D sun;  // Fixed reference
	
	Vector3Cartesian earth_velocity({1.0, 0.0, 0.0});
	Vector3Cartesian earth_initial({100.0, 0.0, 0.0});
	InertialFrame3D earth(&sun, earth_velocity, earth_initial);
	
	Vector3Cartesian moon_velocity({0.0, 0.5, 0.0});
	Vector3Cartesian moon_initial({10.0, 0.0, 0.0});
	InertialFrame3D moon(&earth, moon_velocity, moon_initial);
	
	REQUIRE(sun.GetChildFrames().size() == 1);
	REQUIRE(earth.GetParentFrame() == &sun);
	REQUIRE(earth.GetChildFrames().size() == 1);
	REQUIRE(moon.GetParentFrame() == &earth);
}

TEST_CASE("Frame hierarchy - Destruction detaches children", "[CoordSystem][Hierarchy][lifecycle]")
{
	ReferenceFrame3D* parent = new ReferenceFrame3D();
	ReferenceFrame3D child;
	
	child.SetParentFrame(parent);
	REQUIRE(child.GetParentFrame() == parent);
	
	delete parent;
	
	// Child should be orphaned after parent destruction
	REQUIRE(child.GetParentFrame() == nullptr);
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           SphericalRotatingFrame                                    ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("SphericalRotatingFrame - Construction", "[CoordSystem][SphericalRotatingFrame]")
{
	ReferenceFrame3D parent;
	Real period = 24.0;
	VectorN<Real, 3> axis({0, 0, 1});
	
	SphericalRotatingFrame frame(&parent, period, axis);
	
	REQUIRE(frame._period == Catch::Approx(24.0));
	REQUIRE(frame._axis[2] == Catch::Approx(1.0));
}

TEST_CASE("SphericalRotatingFrame - GetPositionAtTime stub", "[CoordSystem][SphericalRotatingFrame]")
{
	ReferenceFrame3D parent;
	SphericalRotatingFrame frame(&parent, 10.0, VectorN<Real, 3>({0, 0, 1}));
	
	Vector3Spherical pos({1.0, Constants::PI/4, Constants::PI/2});
	VectorN<Real, 3> result = frame.GetPositionAtTime(pos, 0.0);
	
	// Current implementation returns zero (stub)
	REQUIRE(result[0] == Catch::Approx(0.0));
	REQUIRE(result[1] == Catch::Approx(0.0));
	REQUIRE(result[2] == Catch::Approx(0.0));
}

///////////////////////////////////////////////////////////////////////////////////////////
///                           NonInertialFrame3D                                        ///
///////////////////////////////////////////////////////////////////////////////////////////

TEST_CASE("NonInertialFrame3D - Construction", "[CoordSystem][NonInertialFrame]")
{
	ReferenceFrame3D parent;
	NonInertialFrame3D frame(&parent);
	
	REQUIRE(frame.GetParentFrame() == &parent);
}

TEST_CASE("NonInertialFrame3D - Default construction", "[CoordSystem][NonInertialFrame]")
{
	NonInertialFrame3D frame;
	
	REQUIRE(frame.GetParentFrame() == nullptr);
}

} // namespace MML::Tests::Core::CoordSystemTests
