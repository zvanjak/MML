///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        CoordSystem.h                                                       ///
///  Description: Coordinate system definitions (Cartesian, Spherical, Cylindrical)   ///
///               Basis vectors, scale factors, and coordinate relationships          ///
///                                                                                   ///
///  Copyright:   (c) 2024-2025 Zvonimir Vanjak                                       ///
///  License:     Licensed under MML dual-license (see LICENSE.md)                    ///
///               - Free for non-commercial use                                       ///
///               - Commercial license available                                      ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_COORD_SYSTEM_H
#define MML_COORD_SYSTEM_H

#include "MMLBase.h"

#include "base/VectorN.h"
#include "base/MatrixNM.h"
#include "base/BaseUtils.h"
#include "base/Geometry3D.h"

#include "core/CoordTransf.h"
#include "core/CoordTransf/CoordTransfSpherical.h"

namespace MML
{
	// FIXED referential frame, Cartesian local coordinates
	class ReferenceFrame3D
	{ 
		ReferenceFrame3D* _parentFrame = nullptr;
		std::vector<ReferenceFrame3D*> _childFrames;

	public:
		ReferenceFrame3D() {}
		ReferenceFrame3D(ReferenceFrame3D* parentFrame)
		{
			SetParentFrame(parentFrame);
		}

		void SetParentFrame(ReferenceFrame3D* parentFrame)
		{
			if (parentFrame != nullptr)
			{
				_parentFrame = parentFrame;
				parentFrame->AddChildFrame(this);
			}
			else
			{
				_parentFrame = nullptr;
			}
		}
		void AddChildFrame(ReferenceFrame3D* childFrame)
		{
			_childFrames.push_back(childFrame);
			childFrame->_parentFrame = this;
		}

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

	class InertialFrame3D : public ReferenceFrame3D
	{
		// u odnosu na drugi referential frame ima samo konstantu brzinu
		Vector3Cartesian _velocity;
		Vector3Cartesian _pos_at_0;

	public:
		InertialFrame3D() {}
		InertialFrame3D(ReferenceFrame3D* parentFrame)
		{
			SetParentFrame(parentFrame);
		}
		InertialFrame3D(ReferenceFrame3D* parentFrame, Vector3Cartesian velocity, Vector3Cartesian pos_at_0) : _velocity(velocity), _pos_at_0(pos_at_0)
		{
			SetParentFrame(parentFrame);
		}

		virtual Vector3Cartesian GetOriginPositionAtTime(Real t) const override 
		{
			return _pos_at_0 + t * _velocity;
		}
		virtual Vector3Cartesian GetLocalPosInParentFrameAtTime(const VectorN<Real, 3> &pos, Real t)  const override
		{
			return GetOriginPositionAtTime(t) + pos;
		}
	};

	class NonInertialFrame3D : public ReferenceFrame3D
	{
	public:
		NonInertialFrame3D() {}
		NonInertialFrame3D(ReferenceFrame3D* parentFrame)
		{
			SetParentFrame(parentFrame);
		}
		// getOriginPositionAtTime - vraca poziciju u odnosu na ReferentialFrame3D
		// getSpeedAtTime - vraca brzinu u odnosu na ReferentialFrame3D
	};

	// ovo je referentni frame koji se rotira oko nekog centra mase, i treba ga zamisliti kao kocku koja rotira oko CM
	class CircleOrbitingFrame3DCartesian : public NonInertialFrame3D
	{
		// koristimo Cartesian sustav - vraca pozicije u odnosu na CM oko kojeg orbitira U CARTESIAN KOORDINATAMA
		// što ukoliko parametre orbite ne zelim u Cartesian sustavu? - nova klasa CircleOrbitingFrame3DSpherical
	public:
		Real _radius;
		Real _speed;
		Real _period;
		Real _angle_at_t0;
		// normal to plane (axis of rotation), kad je ravnina orbite zakrenuta
		Vector3Cartesian _axis;

		CircleOrbitingFrame3DCartesian(ReferenceFrame3D* parentFrame, Real radius, Real period) : NonInertialFrame3D(parentFrame)
		{
			_radius = radius;
			// _speed = speed; // izracunati
			_period = period;
			_angle_at_t0 = 0;
		}
		CircleOrbitingFrame3DCartesian(ReferenceFrame3D* parentFrame, Real radius, Real period, Real angle_at_t0) : NonInertialFrame3D(parentFrame)
		{
			_radius = radius;
			_period = period;
			_angle_at_t0 = angle_at_t0;
		}
		CircleOrbitingFrame3DCartesian(ReferenceFrame3D* parentFrame, Real radius, Real period, Real orbit_inclination, Real angle_at_t0) : NonInertialFrame3D(parentFrame)
		{
			_radius = radius;
			_period = period;
			_angle_at_t0 = angle_at_t0;
			// calc axis based on inclination
		}

		virtual Vector3Cartesian GetOriginPositionAtTime(Real t) const override
		{
			// calculate rotational evolution of position of center of mass
			Real angle = _angle_at_t0 + 2 * Constants::PI * t / _period;
			// in z-plane!
			Vector3Cartesian CM_pos({ _radius * cos(angle), _radius * sin(angle), 0 });

			return CM_pos;
		}

		virtual Vector3Cartesian GetLocalPosInParentFrameAtTime(const VectorN<Real, 3> &pos, Real t) const override
		{
			// calculate evolution of position of center of mass
			Vector3Cartesian CM_pos = GetOriginPositionAtTime(t);

			// add local coordinates to CM position
			// BITNA PRETPOSTAVKA - kako naš sustav rotira oko CM, njegova apsolutna orijentacije se ne mijenja
			// ie, Zemljina (lokalna) os rotacije je jednom nagnuta OD Sunca, a za sest mjeseci nagnuta PREMA Suncu
			return CM_pos + pos;
		}
	};

	class RotatingFrame3D : public NonInertialFrame3D
	{
		// ovaj frame koristi cilindrični sustav (generalna rotacija oko osi)
	public:
		Real _period;
		Real _angle_at_t0;
		VectorN<Real, 3> _axis;     // pretpostavljamo z-axis za pocetak

		RotatingFrame3D(ReferenceFrame3D* parentFrame, Real period, VectorN<Real, 3> axis) : NonInertialFrame3D(parentFrame)
		{
			_period = period;
			_axis = axis;
		}

		// get child origin position at T, with reference to parent frame
		virtual Vector3Cartesian GetOriginPositionAtTime(Real t) const override
		{
			return Vector3Cartesian({ 0,0,0 });
		}

		// get LocalPoint position at time T, in parent frame
		virtual Vector3Cartesian GetLocalPosInParentFrameAtTime(const VectorN<Real, 3> &pos, Real t) const override
		{
			// pos is in cylindrical coordinates
			
			Real angle = fmod(2 * Constants::PI / _period * t, 2 * Constants::PI);

			return VectorN<Real, 3>({ _axis[0] * cos(angle), _axis[1] * sin(angle), pos[2] });
		}
	};

	class SphericalRotatingFrame : public RotatingFrame3D
	{
		// ovaj frame radi sa spherical koordinatama
	public:
		SphericalRotatingFrame(ReferenceFrame3D* parentFrame, Real period, VectorN<Real, 3> axis)
			: RotatingFrame3D(parentFrame, period, axis) {}

		VectorN<Real, 3> GetPositionAtTime(Vector3Spherical pos, Real t)
		{
			return VectorN<Real, 3>({ 0,0,0 });
		}
	};

	class HardSphereRotatingFrame : public RotatingFrame3D
	{
		// lokalne koordinate - lat, long, h 
		// ima svoj CENTAR MASE u sredini sfere, i u odnosu na njega vraća pozicije
		// koje su usuglasene s axisom rotacije (lat, long)
	public:
		// axis rotacije zadan base klasom
		double _radius;

		HardSphereRotatingFrame(ReferenceFrame3D* parentFrame, Real radius, Real period, VectorN<Real, 3> axis)
			: _radius(radius), RotatingFrame3D(parentFrame, period, axis) {}

		// get origin position at T, with reference to parent frame
		virtual Vector3Cartesian GetOriginPositionAtTime(Real t) const override
		{
			// it is NOT moving, only rotating
			return Vector3Cartesian({ 0,0,0 });
		}

		// get LocalPoint position at time T, in parent frame, meaning Cartesian coordinates
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

	class RotatingSphereLocalCartesian : public InertialFrame3D
	{
		// u ctor dobije ref na HardSphereRotatingFrame, I TOCNO ODREDJENU TOCKU NA SFERI!!!
		// ima smisla - gleda nakon deltaT gdje je pozicija tocke u jednom i drugom
		const HardSphereRotatingFrame& _parentFrame;
		const Real _latitude;
		const Real _longitude;
	public:
		RotatingSphereLocalCartesian(const HardSphereRotatingFrame &parent, Real latitude, Real longitude)  : _parentFrame(parent), _latitude(latitude), _longitude(longitude)
		{
		}
		// kosi hitac zadan u lokalnom kartezijevom, i izracunam
		// onda vidim gdje je taj lokalni kartezije u trenutku deltaT, i da li se 
		// slaze TRENUTNA tocka (x,y,z) di je sletio hitac, s onom kako sam izracunao

		// get origin position at T, with reference to parent frame
		// UZETI U OBZIR DA SE ORIGIN TOCKA ROTIRA!!!
		virtual Vector3Cartesian GetOriginPositionAtTime(Real t) const override
		{
			// depends on HardSphere rotating period
			return Vector3Cartesian({ 0,0,0 });
		}

		// return type NISU CARTESIAN!!! - to su lat, long, h
		virtual Vector3Cartesian GetLocalPosInParentFrameAtTime(const VectorN<Real, 3> &localPos, Real t) const override
		{
			// VRIJEME IGNORIRAMO!!!
			// ukljuceno je u parent kalkulacije, i ovdje samo vracamo transf lokalnih (x, y, z) u (lat, long, h)
			Vector3Cartesian pos(localPos);

			double one_km_in_lat_deg = 1 / 111.32;
			double one_lat_deg_in_km = 2 * _parentFrame._radius * cos(Utils::DegToRad(_latitude)) * Constants::PI / 360;
			double one_km_in_long_deg = 1 / one_lat_deg_in_km;

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
