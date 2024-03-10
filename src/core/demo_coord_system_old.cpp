#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/CoordSystem.h"
#include "core/CoordTransf.h"
#endif

using namespace MML;

void Demo_Inertial()
{
	// zadati isti origin, neku brzinu, i vidjeti kako ide transformacija
}

/*
	SUSTINA
	- imam koordinate u jednom, i (lokalne) koordinate u drugom koordinatnom sistemu
	- kakva je veza?
	- i kako izracunati jedne na osnovu drugih?

	PlanarRotatingSystem disk_rotation(pocetni phi, brzina rotacije);
	- za dane dvije koord, lat i long, daje poziciju u odnosu na dani fiksni koord sustav
	LocalCartesian disk_surface(disk_rotation, lat, long);

	- što izracunati?
			- artiljerijski hitac s dane pozicije i po danoj paraboli
			- gdje ce pasti - koordinate u jednom i drugom sustavu

	- i onda još dodati vrtuljak na toj površini!

	MovingDynamicalSytem3D earth_around_sun(funkcija ovisnosti pozicije u odnosu na GLOBALNI KARTEZIJEV sustav);
	RotatingSystem3D earth_rotation(earth_around_sun);
	- za dane dvije koord, lat i long, daje poziciju u odnosu na dani koord sustav
	LocalCartesian3D earth_surface(earth_rotation, lat, long);

	LorentzInertialMovingFrame observer_moving_frame(vektor smjera, ovisnost pozicije o t); /// moze i (0,0,0,0) - stoji na mjestu
	LorentzInertialMovingFrame s1(vektor smjera, ovisnost pozicije o t);
	LorentzInertialMovingFrame s2(vektor smjera, ovisnost pozicije o t);

	LocalLorent s1;
	LorentzBoosted s2;
	LorentTranslated s3;
	*/

void Demo_CoordSystem_old()
{
	std::cout << std::endl;
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                         COORD SYSTEM                          ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	//InertialMovingFrame train(Vector3Cartesian{ 0,0,0 }, Vector3Cartesian{ 1,0,0 });
	//InertialMovingFrame table_moving_on_train(Vector3Cartesian{ 1,0,0 }, Vector3Cartesian{ 0.5,0,0 });
	//RotatingFrame disc_spinning_on_table;
	//InertialMovingFrame ball_thrown_up;

	//// point on disc - polar coordinates
	//Vector3Cartesian point_on_disc{ 0.5, 0.5, 0 };
	//auto transf0 = disc_spinning_on_table.transf(point_on_disc, 1.0);
	//auto transf1 = table_moving_on_train.transf(transf0, 1.0);
	//auto transf2 = train.transf(transf1, 1.0);
	/*
			// TODO - simulation of Solar system

			- sunce nam je u ishodistu i definira koord sustav
			- moving coord system earth_rotation_around_sun, jupier around sun, mars around sun
					- elipsa s centrom mase u zaristu
					- definira tocku centra planete u odnosu na ishodiste
			- rotating coord system earth_rotation_around_earth_center
					- sferni coord system
					- daje transformaciju iz lokalne (lat, long, height) coord u globalnu (x,y,z) coord
			- local observer coord system - "nebo"
					- za danu lat, long tocku na zemlji
					- definira oko te tocke sferni sustav, gdje je azimut poravnan sa smjerom sjevera kako treba
			- SUSTINA - zelimo transformaciju jupitera i marsa na taj sustav, da vidimo gdje se nalaze na nebo
	*/

	//RotatingFrame earthSolarSystemCM;    // opisuje kružnicu oko sunca - KLJUČNO JE PORAVNANJE!!!
	// T = 0 - Greenwhich je točno na x osi
	// MovingCoordSystem<3> earthRotational;       // kao koordinate ima sferne koordinate iznad površine zemlje
	// MovingCoordSystem<3> earthLocal;            // projekcija na 2D coord system, sa zadanim ishodištem

	// ČIM GIBANJE čestice po putanji odstupa od tangente, znači da je prisutna sila

	// najprije u 2D ?

	// imamo statičnu transf - sve zadano, pa da vidimo

	// imamo simulaciju - ovisnost o t, pa se mijenja

	// usporediti
	//  1. kompletnqa simulacija zemlje oko sunca i topovske granate ispaljene
	//  2. zemlja se zadano giva oko sunca, a imamo lokalni rotacijski sustav
}