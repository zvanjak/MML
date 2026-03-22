///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        CoordSystem.h                                                       ///
///  Description: Coordinate system definitions (Cartesian, Spherical, Cylindrical)   ///
///               Basis vectors, scale factors, and coordinate relationships          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COORD_SYSTEM_H
#define MML_COORD_SYSTEM_H

#include "MMLBase.h"

#include <vector>

#include "base/Vector/VectorN.h"
#include "base/Matrix/MatrixNM.h"
#include "base/BaseUtils.h"
#include "mml/base/Geometry/Geometry3D.h"

#include "core/CoordTransf.h"
#include "core/CoordTransf/CoordTransfSpherical.h"

namespace MML
{
	///////////////////////////////////////////////////////////////////////////////
	// REFERENCE FRAME CONVENTIONS
	//
	//   Default frame:    Right-handed Cartesian (x, y, z)
	//
	//   Frame hierarchy:  Parent-child tree using non-owning raw pointers.
	//                     Each child frame has position and velocity relative
	//                     to its parent. Root frame has no parent (nullptr).
	//
	//   Time parameter:   Frame positions/velocities are functions of time t.
	//                     Units of t are user-defined (typically seconds).
	//
	//   Kinematics:       Galilean (non-relativistic) — velocities add linearly.
	//                     Position in parent = child_origin + child_offset.
	//
	//   Coordinate systems: Coordinate systems (Cartesian, Spherical, Cylindrical)
	//                     follow the conventions defined in CoordTransfSpherical.h
	//                     and CoordTransfCylindrical.h.
	//
	// See also: CoordTransf.h, CoordTransfSpherical.h, CoordTransfCylindrical.h
	///////////////////////////////////////////////////////////////////////////////

	/// @brief Base class for 3D reference frames with hierarchical parent-child relationships
	/// 
	/// Provides foundation for coordinate transformations between reference frames.
	/// Supports frame hierarchies where child frames can have positions/velocities relative to parent frames.
	/// Default implementation represents a fixed Cartesian reference frame.
	class ReferenceFrame3D
	{ 
		ReferenceFrame3D* _parentFrame = nullptr;
		std::vector<ReferenceFrame3D*> _childFrames;

		/// @brief Check if 'frame' is an ancestor of this frame (cycle detection)
		/// @param frame Frame to check
		/// @return true if 'frame' is this frame or any of its ancestors
		bool IsAncestor(const ReferenceFrame3D* frame) const
		{
			const ReferenceFrame3D* current = this;
			while (current != nullptr)
			{
				if (current == frame)
					return true;
				current = current->_parentFrame;
			}
			return false;
		}

		/// @brief Remove a child from the children list (internal helper)
		/// @param childFrame Child to remove
		void RemoveChildInternal(ReferenceFrame3D* childFrame)
		{
			auto it = std::find(_childFrames.begin(), _childFrames.end(), childFrame);
			if (it != _childFrames.end())
				_childFrames.erase(it);
		}

	public:
		/// @brief Default constructor - creates independent reference frame
		ReferenceFrame3D() {}
		
		/// @brief Constructor with parent frame
		/// @param parentFrame Parent reference frame (can be nullptr)
		ReferenceFrame3D(ReferenceFrame3D* parentFrame)
		{
			SetParentFrame(parentFrame);
		}

		/// @brief Virtual destructor - detaches from parent and orphans children
		/// 
		/// Lifetime model: Non-owning raw pointers. Parent does not own children.
		/// When a frame is destroyed, children become orphaned (parent = nullptr).
		/// Users must ensure frames outlive their usage or manage lifetimes externally.
		virtual ~ReferenceFrame3D()
		{
			// Detach from parent
			if (_parentFrame != nullptr)
				_parentFrame->RemoveChildInternal(this);

			// Orphan all children (set their parent to nullptr)
			for (auto* child : _childFrames)
				child->_parentFrame = nullptr;
			
			_childFrames.clear();
		}

		/// @brief Set or change the parent reference frame
		/// @param parentFrame New parent frame (nullptr to remove parent)
		/// @return true if parent was set successfully, false if cycle would be created
		/// @note Automatically detaches from old parent before reparenting
		bool SetParentFrame(ReferenceFrame3D* parentFrame)
		{
			// No-op if already the parent
			if (parentFrame == _parentFrame)
				return true;

			// Prevent cycles: can't set parent if it would create a cycle
			if (parentFrame != nullptr && parentFrame->IsAncestor(this))
				return false;  // Would create a cycle

			// Detach from old parent first
			if (_parentFrame != nullptr)
				_parentFrame->RemoveChildInternal(this);

			_parentFrame = parentFrame;

			if (parentFrame != nullptr)
				parentFrame->AddChildFrame(this);

			return true;
		}

		/// @brief Add a child reference frame
		/// @param childFrame Child frame to add
		/// @return true if child was added, false if null, duplicate, or would create cycle
		bool AddChildFrame(ReferenceFrame3D* childFrame)
		{
			if (childFrame == nullptr)
				return false;

			// Prevent cycles: can't add child if this frame is already a descendant of it
			if (IsAncestor(childFrame))
				return false;

			// Prevent duplicates
			if (std::find(_childFrames.begin(), _childFrames.end(), childFrame) != _childFrames.end())
				return true;  // Already a child, no-op (return true since it's already there)

			// Detach child from its old parent
			if (childFrame->_parentFrame != nullptr)
				childFrame->_parentFrame->RemoveChildInternal(childFrame);

			_childFrames.push_back(childFrame);
			childFrame->_parentFrame = this;
			return true;
		}

		/// @brief Remove a child reference frame
		/// @param childFrame Child to remove
		/// @return true if child was removed, false if not found
		bool RemoveChildFrame(ReferenceFrame3D* childFrame)
		{
			auto it = std::find(_childFrames.begin(), _childFrames.end(), childFrame);
			if (it != _childFrames.end())
			{
				(*it)->_parentFrame = nullptr;
				_childFrames.erase(it);
				return true;
			}
			return false;
		}

		/// @brief Get the parent reference frame
		/// @return Parent frame or nullptr if this is a root frame
		ReferenceFrame3D* GetParentFrame() const { return _parentFrame; }

		/// @brief Get the child reference frames
		/// @return Read-only reference to children vector
		const std::vector<ReferenceFrame3D*>& GetChildFrames() const { return _childFrames; }

		// get child origin position at T, with reference to parent frame
		virtual Vector3Cartesian GetOriginPositionAtTime(Real t) const
		{
			return Vector3Cartesian({ 0,0,0 });
		}

		// get LocalPoint position at time T, in parent frame
		virtual Vector3Cartesian GetLocalPosInParentFrameAtTime(const VectorN<Real, 3> &pos, Real t) const 
		{
			return pos;
		}
	};

	/// @brief Inertial reference frame moving with constant velocity
	/// 
	/// Represents a reference frame that moves with constant velocity relative to its parent.
	/// Position evolves linearly with time: position(t) = initial_position + velocity * t
	class InertialFrame3D : public ReferenceFrame3D
	{
		Vector3Cartesian _velocity;
		Vector3Cartesian _pos_at_0;

	public:
		/// @brief Default constructor - creates stationary inertial frame
		InertialFrame3D() {}
		
		/// @brief Constructor with parent frame
		/// @param parentFrame Parent reference frame
		InertialFrame3D(ReferenceFrame3D* parentFrame)
		{
			SetParentFrame(parentFrame);
		}
		/// @brief Constructor with parent frame, velocity, and initial position
		/// @param parentFrame Parent reference frame
		/// @param velocity Constant velocity vector
		/// @param pos_at_0 Initial position at t=0
		InertialFrame3D(ReferenceFrame3D* parentFrame, Vector3Cartesian velocity, Vector3Cartesian pos_at_0) : _velocity(velocity), _pos_at_0(pos_at_0)
		{
			SetParentFrame(parentFrame);
		}

		/// @brief Get origin position at time t (linear motion with constant velocity)
		/// @param t Time parameter
		/// @return Position = initial_position + velocity * t
		virtual Vector3Cartesian GetOriginPositionAtTime(Real t) const override 
		{
			return _pos_at_0 + t * _velocity;
		}
		
		/// @brief Transform local position to parent frame (includes frame motion)
		/// @param pos Local position in this frame
		/// @param t Time parameter
		/// @return Position in parent frame = origin_at_t + local_pos
		virtual Vector3Cartesian GetLocalPosInParentFrameAtTime(const VectorN<Real, 3> &pos, Real t)  const override
		{
			return GetOriginPositionAtTime(t) + pos;
		}
	};

	/// @brief Non-inertial reference frame (accelerating or rotating)
	/// 
	/// Base class for reference frames that are accelerating or rotating.
	/// Derived classes implement specific non-inertial motion (circular orbits, rotation, etc.)
	class NonInertialFrame3D : public ReferenceFrame3D
	{
	public:
		/// @brief Default constructor
		NonInertialFrame3D() {}
		
		/// @brief Constructor with parent frame
		/// @param parentFrame Parent reference frame
		NonInertialFrame3D(ReferenceFrame3D* parentFrame)
		{
			SetParentFrame(parentFrame);
		}
		// getOriginPositionAtTime - vraca poziciju u odnosu na ReferentialFrame3D
		// getSpeedAtTime - vraca brzinu u odnosu na ReferentialFrame3D
	};

	/// @brief Reference frame in circular orbit around center of mass
	/// 
	/// Models a reference frame orbiting in a circular path around its parent's origin.
	/// Uses Cartesian coordinates for orbit parameters (radius, period).
	/// The frame itself maintains fixed orientation (no rotation relative to parent).
	class CircleOrbitingFrame3DCartesian : public NonInertialFrame3D
	{
	public:
		Real _radius;          ///< Orbital radius from parent origin
		Real _speed;           ///< Orbital speed (calculated from period)
		Real _period;          ///< Orbital period (time for full revolution)
		Real _angle_at_t0;     ///< Initial angle at t=0 (in radians)
		Vector3Cartesian _axis; ///< Normal to orbital plane (rotation axis)

		/// @brief Constructor with radius and period (orbit in XY plane)
		/// @param parentFrame Parent reference frame
		/// @param radius Orbital radius
		/// @param period Orbital period
		CircleOrbitingFrame3DCartesian(ReferenceFrame3D* parentFrame, Real radius, Real period) : NonInertialFrame3D(parentFrame)
		{
			_radius = radius;
			// _speed = speed; // izracunati
			_period = period;
			_angle_at_t0 = 0;
		}
		/// @brief Constructor with initial angle
		/// @param parentFrame Parent reference frame
		/// @param radius Orbital radius
		/// @param period Orbital period
		/// @param angle_at_t0 Initial angle at t=0
		CircleOrbitingFrame3DCartesian(ReferenceFrame3D* parentFrame, Real radius, Real period, Real angle_at_t0) : NonInertialFrame3D(parentFrame)
		{
			_radius = radius;
			_period = period;
			_angle_at_t0 = angle_at_t0;
		}
		/// @brief Constructor with orbital inclination
		/// @param parentFrame Parent reference frame
		/// @param radius Orbital radius
		/// @param period Orbital period
		/// @param orbit_inclination Inclination of orbital plane
		/// @param angle_at_t0 Initial angle at t=0
		CircleOrbitingFrame3DCartesian(ReferenceFrame3D* parentFrame, Real radius, Real period, Real orbit_inclination, Real angle_at_t0) : NonInertialFrame3D(parentFrame)
		{
			_radius = radius;
			_period = period;
			_angle_at_t0 = angle_at_t0;
			// calc axis based on inclination
		}

		/// @brief Get orbital position at time t
		/// @param t Time parameter
		/// @return Position on circular orbit (in XY plane by default)
		virtual Vector3Cartesian GetOriginPositionAtTime(Real t) const override
		{
			Real angle = _angle_at_t0 + 2 * Constants::PI * t / _period;
			// in z-plane!
			Vector3Cartesian CM_pos({ _radius * cos(angle), _radius * sin(angle), 0 });

			return CM_pos;
		}

		/// @brief Transform local position to parent frame (orbital motion + local offset)
		/// @param pos Local position in this frame
		/// @param t Time parameter
		/// @return Position in parent frame = orbital_position(t) + local_pos
		/// @note Frame orientation remains fixed (no rotation)
		virtual Vector3Cartesian GetLocalPosInParentFrameAtTime(const VectorN<Real, 3> &pos, Real t) const override
		{
			Vector3Cartesian CM_pos = GetOriginPositionAtTime(t);

			// add local coordinates to CM position
			// BITNA PRETPOSTAVKA - kako naš sustav rotira oko CM, njegova apsolutna orijentacije se ne mijenja
			// ie, Zemljina (lokalna) os rotacije je jednom nagnuta OD Sunca, a za sest mjeseci nagnuta PREMA Suncu
			return CM_pos + pos;
		}
	};

	/// @brief Rotating reference frame (rotation around arbitrary axis)
	/// 
	/// Models a reference frame rotating with constant angular velocity around an axis.
	/// Uses cylindrical coordinate system for representing rotating positions.
	/// Origin remains fixed; only orientation changes with time.
	class RotatingFrame3D : public NonInertialFrame3D
	{
	public:
		Real _period;              ///< Rotation period (time for full rotation)
		Real _angle_at_t0;         ///< Initial rotation angle at t=0
		VectorN<Real, 3> _axis;    ///< Rotation axis (default: z-axis)

		/// @brief Constructor for rotating frame
		/// @param parentFrame Parent reference frame
		/// @param period Rotation period
		/// @param axis Rotation axis vector
		RotatingFrame3D(ReferenceFrame3D* parentFrame, Real period, VectorN<Real, 3> axis) : NonInertialFrame3D(parentFrame)
		{
			_period = period;
			_axis = axis;
		}

		/// @brief Get origin position (fixed at parent origin)
		/// @param t Time parameter (unused for rotating frame)
		/// @return Zero vector (origin doesn't move)
		virtual Vector3Cartesian GetOriginPositionAtTime(Real t) const override
		{
			return Vector3Cartesian({ 0,0,0 });
		}

		/// @brief Transform rotating cylindrical position to parent Cartesian coordinates
		/// @param pos Position in cylindrical coordinates (r, φ, z)
		/// @param t Time parameter (determines rotation angle)
		/// @return Position in parent Cartesian frame
		virtual Vector3Cartesian GetLocalPosInParentFrameAtTime(const VectorN<Real, 3> &pos, Real t) const override
		{
			
			Real angle = fmod(2 * Constants::PI / _period * t, 2 * Constants::PI);

			return VectorN<Real, 3>({ _axis[0] * cos(angle), _axis[1] * sin(angle), pos[2] });
		}
	};

	/// @brief Rotating frame using spherical coordinates
	/// 
	/// Extension of RotatingFrame3D that uses spherical coordinates (r, θ, φ)
	/// instead of cylindrical. Useful for planetary rotation models.
	class SphericalRotatingFrame : public RotatingFrame3D
	{
	public:
		/// @brief Constructor for spherical rotating frame
		/// @param parentFrame Parent reference frame
		/// @param period Rotation period
		/// @param axis Rotation axis vector
		SphericalRotatingFrame(ReferenceFrame3D* parentFrame, Real period, VectorN<Real, 3> axis)
			: RotatingFrame3D(parentFrame, period, axis) {}

		/// @brief Get position at time t from spherical coordinates
		/// @param pos Position in spherical coordinates (r, θ, φ)
		/// @param t Time parameter
		/// @return Position vector (implementation needed)
		VectorN<Real, 3> GetPositionAtTime(Vector3Spherical pos, Real t)
		{
			return VectorN<Real, 3>({ 0,0,0 });
		}
	};

	/// @brief Rotating spherical body (e.g., rotating planet)
	/// 
	/// Models a rigid sphere rotating around its axis (e.g., Earth).
	/// Local coordinates are latitude, longitude, and altitude.
	/// Center of mass at sphere center. Coordinates aligned with rotation axis.
	class HardSphereRotatingFrame : public RotatingFrame3D
	{
	public:
		double _radius;  ///< Radius of the rotating sphere

		/// @brief Constructor for rotating sphere
		/// @param parentFrame Parent reference frame
		/// @param radius Sphere radius
		/// @param period Rotation period
		/// @param axis Rotation axis vector
		HardSphereRotatingFrame(ReferenceFrame3D* parentFrame, Real radius, Real period, VectorN<Real, 3> axis)
			: _radius(radius), RotatingFrame3D(parentFrame, period, axis) {}

		/// @brief Get origin position (sphere center is fixed)
		/// @param t Time parameter (unused)
		/// @return Zero vector (only rotating, not moving)
		virtual Vector3Cartesian GetOriginPositionAtTime(Real t) const override
		{
			return Vector3Cartesian({ 0,0,0 });
		}

		/// @brief Transform local (lat, long, altitude) to parent Cartesian coordinates
		/// @param localPos Local position as (latitude_deg, longitude_deg, altitude)
		/// @param t Time parameter (determines rotation)
		/// @return Position in parent Cartesian frame
		virtual Vector3Cartesian GetLocalPosInParentFrameAtTime(const VectorN<Real, 3> &localPos, Real t) const override
		{
			Real latitudeDeg = localPos[0];
			Real longitudeDeg = localPos[1];
			Real h = localPos[2];

			// formirati sferni vektor, i vidjeti koliko se zarotira za T
		Vector3Spherical spherePos({Real(_radius + h), Real(Utils::DegToRad(90 - latitudeDeg)), Real(Utils::DegToRad(longitudeDeg)) });
			// transf u Cartesian
			Vector3Cartesian cartPos = CoordTransfSpherToCart.transf(spherePos);

			return cartPos;
		}
	};

	/// @brief Local Cartesian frame on rotating sphere surface
	/// 
	/// Represents local Cartesian coordinates (x=east, y=north, z=up) at a specific
	/// point on a rotating sphere. Used for modeling trajectories on rotating planets
	/// (e.g., projectile motion on Earth accounting for rotation).
	class RotatingSphereLocalCartesian : public InertialFrame3D
	{
		const HardSphereRotatingFrame& _parentFrame;  ///< Reference to parent rotating sphere
		const Real _latitude;                          ///< Latitude of local frame origin (degrees)
		const Real _longitude;                         ///< Longitude of local frame origin (degrees)
	public:
		/// @brief Constructor for local Cartesian frame on sphere
		/// @param parent Parent rotating sphere frame
		/// @param latitude Latitude of origin point (degrees)
		/// @param longitude Longitude of origin point (degrees)
		RotatingSphereLocalCartesian(const HardSphereRotatingFrame &parent, Real latitude, Real longitude)  : _parentFrame(parent), _latitude(latitude), _longitude(longitude)
		{
		}
		// kosi hitac zadan u lokalnom kartezijevom, i izracunam
		// onda vidim gdje je taj lokalni kartezije u trenutku deltaT, i da li se 
		// slaze TRENUTNA tocka (x,y,z) di je sletio hitac, s onom kako sam izracunao

		/// @brief Get origin position (accounts for sphere rotation)
		/// @param t Time parameter
		/// @return Origin position in parent frame (TODO: implement rotation)
		virtual Vector3Cartesian GetOriginPositionAtTime(Real t) const override
		{
			return Vector3Cartesian({ 0,0,0 });
		}

		/// @brief Transform local Cartesian (x=east, y=north, z=up) to (lat, long, altitude)
		/// @param localPos Local Cartesian position (x, y, z) in meters
		/// @param t Time parameter (currently ignored, rotation handled by parent)
		/// @return Geographic coordinates (latitude_deg, longitude_deg, altitude)
		virtual Vector3Cartesian GetLocalPosInParentFrameAtTime(const VectorN<Real, 3> &localPos, Real t) const override
		{
			// VRIJEME IGNORIRAMO!!!
			// ukljuceno je u parent kalkulacije, i ovdje samo vracamo transf lokalnih (x, y, z) u (lat, long, h)
			Vector3Cartesian pos(localPos);

			double one_km_in_lat_deg = 1 / 111.32;
			double one_lon_deg_in_km = 2 * _parentFrame._radius * cos(Utils::DegToRad(_latitude)) * Constants::PI / 360;
			double one_km_in_long_deg = 1 / one_lon_deg_in_km;

			// local x axis is oriented towards east, y towards north
			// so dx is in longitude direction, dy in latitude direction
			double latitudeDeg = _latitude + pos[1] * one_km_in_lat_deg;
			double longitudeDeg = _longitude + pos[0] * one_km_in_long_deg;
			double h = localPos[2];								// h is equal to z locally


			Vector3Cartesian ret(latitudeDeg, longitudeDeg, h);

			return ret;
		}
	};

	// da li mi treba Local3D koji za parenta ima HardSphereRotatingFrameToSpherical?
	// lokalni sustav, baziran na TOCNO ODREDJENOJ TOCKI SFERE, s x, y i z
	// za njega NE TREBA davati lat, long i h jer vec ima, a x, y i z transformira lokalno

 
}

#endif
