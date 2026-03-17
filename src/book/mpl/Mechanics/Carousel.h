#if !defined MPL_CAROUSEL_H
#define MPL_CAROUSEL_H

#include "MMLBase.h"


using namespace MML;

namespace MPL
{
	// Local Cartesian coordinate system on the rotating carousel, describing positions with x and y coordinates.
	// As it rotates, all transformations to GlobalCartSystem need to take into account the rotation of the carousel
	// and have t as an additional parameter
	// At t=0, orientation of the local cylindrical coordinate system is such that
	// local x-axis is pointing in the direction of the global x-axis,
	// and local y-axis is pointing in the direction of the global y-axis
	class RotatingCarouselFrame_LocalCartSystem
	{
		Real _w;			// angular velocity of the carousel

	public:
		RotatingCarouselFrame_LocalCartSystem(Real w) : _w(w) {}

		// Transforming local Cartesian coordinates to global Cartesian coordinates
		void transfToGlobal(Real t,
												const Real& x_local, const Real& y_local, const Real& z_local,
												Real& x_global, Real& y_global, Real& z_global)
		{
			Real phi_global = _w * t; // global azimuthal angle considering carousel rotation
			x_global = x_local * cos(phi_global) - y_local * sin(phi_global);
			y_global = x_local * sin(phi_global) + y_local * cos(phi_global);
			z_global = z_local;
		}
		// Transforming global Cartesian coordinates to local Cartesian coordinates
		void transfFromGlobal(Real t,
													const Real& x_global, const Real& y_global, const Real& z_global,
													Real& x_local, Real& y_local, Real& z_local)
		{
			Real phi_global = _w * t; // global azimuthal angle considering carousel rotation
			x_local = x_global * cos(-phi_global) - y_global * sin(-phi_global);
			y_local = x_global * sin(-phi_global) + y_global * cos(-phi_global);
			z_local = z_global;
		}

		// transfroming local Cartesian coordinates to local cylindrical coordinates
		void transfToLocalCyl(const Real& x_local, const Real& y_local, const Real& z_local,
													Real& r_cyl, Real& phi_cyl, Real& z_cyl)
		{
			// r_cyl - radial distance from the center of carousel in local cylindrical coordinates
			// phi_cyl - azimuthal angle in local cylindrical coordinates
			r_cyl = std::sqrt(x_local * x_local + y_local * y_local);
			if (r_cyl > 0.0)
			{
				phi_cyl = std::atan2(y_local, x_local); // azimuthal angle in local cylindrical coordinates
			}
			else
			{
				phi_cyl = 0.0; // if r_local is zero, set phi_local to zero
			}
			z_cyl = z_local; // z-coordinate remains the same in both coordinate systems
		}
		// Transforming local cylindrical coordinates to local Cartesian coordinates
		void transfFromLocalCyl(const Real& r_cyl, const Real& phi_cyl, const Real& z_cyl,
														Real& x_local, Real& y_local, Real& z_local)
		{
			x_local = r_cyl * cos(phi_cyl);
			y_local = r_cyl * sin(phi_cyl);
			z_local = z_cyl;
		}
	};

	// Local cylindrical coordinate system on the rotating carousel, describing positions with R and phi coordinates.
	// As it rotates, all transformations to GlobalCartSystem need to take into account the rotation of the carousel
	// and have t as an additional parameter
	// At t=0, orientation of the local cylindrical coordinate system is such that
	// "meridian" with phi=0 is pointing in the direction of the global x-axis,
	// and "meridian" with phi=pi/2 is pointing in the direction of the global y-axis
	class RotatingCarouselFrame_LocalCylSystem
	{
		Real _w;			// angular velocity of the carousel

	public:
		RotatingCarouselFrame_LocalCylSystem(Real w) : _w(w) {}

		// Transforming local cylindrical coordinates to global Cartesian coordinates
		void transfToGlobal(Real t,
												const Real& r_local, const Real& phi_local, const Real& z_local,
												Real& x_global, Real& y_global, Real& z_global)
		{
			Real phi_global = phi_local + _w * t;

			x_global = r_local * cos(phi_global);
			y_global = r_local * sin(phi_global);
			z_global = z_local;
		}

		// Transforming global Cartesian coordinates to local cylindrical coordinates
		void transfFromGlobal(Real t,
													const Real& x_global, const Real& y_global, const Real& z_global,
													Real& r_local, Real& phi_local, Real& z_local)
		{
			Real phi_global = std::atan2(y_global, x_global); // global azimuthal angle in global Cartesian coordinates

			phi_local = phi_global - _w * t;
			r_local = std::sqrt(x_global * x_global + y_global * y_global);
			z_local = z_global;
		}
	};
	// make Focoult Pendulum
}

#endif