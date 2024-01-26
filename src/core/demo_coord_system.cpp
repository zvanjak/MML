#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "core/CoordTransf.h"
#endif

using namespace MML;

// TODO 0.8 - 2D verzije!!!
// carousel - rotirajuci sustav
// 1. simple- carousel promjera 10 km, i na njemu ljudi
// 2. carousel promjera 10 km, i na njemu ljudi, i na njemu je jos jedan carousel promjera 1 km
// na njemu se odvija streljačko natjecanje
class ReferentialFrame3D
{ };
class InertialFrame3D : public ReferentialFrame3D
{ };
class MovingFrame3D : public ReferentialFrame3D
{ 
    // getOriginPositionAtTime - vraca poziciju u odnosu na ReferentialFrame3D
    // getSppedAtTime - vraca brzinu u odnosu na ReferentialFrame3D
};

// TODO - CircleOrbitingCartesianFrame3DToSpherical
// simulacija solarnog sustava u sfernim koordinatama oko centra mase
class CircleOrbitingCartesianFrame3DToCartesian : public MovingFrame3D
{
    // koristimo Cartesian sustav - vraca pozicije u odnosu na CM oko kojeg orbitira U CARTESIAN KOORDINATAMA
    // što ukoliko parametre orbite ne zelim u Cartesian sustavu? - nova klasa
public:
    double _radius;
    double _speed;
    double _period;
    double _angle_at_t0;
    // axis, za 3d slucaj
    ReferentialFrame3D _parentFrame;

    CircleOrbitingCartesianFrame3DToCartesian(ReferentialFrame3D parentFrame, double radius, double period)
    {
        _radius = radius;
        // _speed = speed; // izracunati
        _period = period;
        _parentFrame = parentFrame;
        _angle_at_t0 = 0;
    }
    CircleOrbitingCartesianFrame3DToCartesian(ReferentialFrame3D parentFrame, double radius, double period, double angle_at_t0)
    {
        _radius = radius;
        // _speed = speed;
        _period = period;
        _angle_at_t0 = angle_at_t0;
        _parentFrame = parentFrame;
    }

    Vector3Cartesian GetOriginPositionAtTime(double t)
    {
        // calculate evolution of position of center of mass
        double angle = _angle_at_t0 + 2 * Constants::PI * t / _period;
        Vector3Cartesian CM_pos({_radius * cos(angle), _radius * sin(angle), 0});

        return CM_pos;
    }
    // za danu tocku zadanu Cartesian koordinatama u lokalnom sustavu
    // vraca Cartesian poziciju (u LOKALNOM frameu u odnosu na CM), nakon vremena T
    Vector3Cartesian GetLocalPositionAtTime(Vector3Cartesian pos, double t)
    {
        // calculate evolution of position of center of mass
        Vector3Cartesian CM_pos = GetOriginPositionAtTime(t);

        // add local coordinates to CM position
        return CM_pos + pos;
    }
    Vector3Cartesian GetLocalPositionAtTime(Vector3Spherical pos, double t)
    {
        // calculate evolution of position of center of mass
        Vector3Cartesian CM_pos = GetOriginPositionAtTime(t);

        // transform given spherical coordinates to cartesian
        // add local coordinates to CM position        
        return Vector3Cartesian({0,0,0});
    }    
};

class RotatingFrame3D : public ReferentialFrame3D
{
    // ovaj frame koristi cilindrični sustav (generalna rotacija oko osi)
public:
    double _period;
    double _angle_at_t0;
    VectorN<Real, 3> _axis;     // pretpostavljamo z-axis za pocetak
    ReferentialFrame3D _parentFrame;

    RotatingFrame3D(ReferentialFrame3D parentFrame, double period, VectorN<Real, 3> axis)
    {
        _period = period;
        _axis = axis;
        _parentFrame = parentFrame;
    }
    VectorN<Real, 3> GetPositionAtTime(Vector3Cylindrical pos, double t)
    {
        double angle = 2 * Constants::PI * t / _period;
        return VectorN<Real, 3>({_axis[0] * cos(angle), _axis[1] * sin(angle), 0});
        // TODO
    }
};

class SphericalRotatingFrame : public RotatingFrame3D
{
    // ovaj frame radi sa spherical koordinatama
public:
    SphericalRotatingFrame(ReferentialFrame3D parentFrame, double period, VectorN<Real, 3> axis) 
        : RotatingFrame3D(parentFrame, period, axis) {}

    VectorN<Real, 3> GetPositionAtTime(Vector3Spherical pos, double t)
    {
        return VectorN<Real, 3>({0,0,0});
    }
};

class HardSphereRotatingFrameToSpherical : public RotatingFrame3D
{
    // lokalne koordinate - lat, long, h 
    // vraca - spherical
    // ima svoj CENTAR MASE u sredini sfere, i u odnosu na njega vraća pozicije
    // koje su usglasene s axisom rotacije (lat, long)
public:
    double _radius;
    
    HardSphereRotatingFrameToSpherical(ReferentialFrame3D parentFrame, double radius, double period, VectorN<Real, 3> axis) 
        : _radius(radius), RotatingFrame3D(parentFrame, period, axis) {}

    Vector3Cartesian GetOriginPositionAtTime(double t)
    {
        return Vector3Cartesian({0,0,0});       // ne mice se!!!!
    }

    // za danu tocku u lokalnom sustavu, vraca spherical poziciju (u LOKALNOM frameu u odnosu na CM), nakon vremena T
    Vector3Spherical GetPositionAtTime(double latitude, double longitude, double height, double t)
    {
        // taking into consideration rotation of the sphere
        // calculate position after time t

        // pretvoriti u pravi spherical vektor
        // 
        return Vector3Spherical({0,0,0});
    }
};

class HardSphereToLocalCartesian : public InertialFrame3D
{
    // u ctor dobije ref na HardSphereRotatingFrameToSpherical
    // ima smisla - gleda nakon deltaT gdje je pozicija tocke u jednom i drugom
    // kosi hitac zadan u lokalnom kartezijevom, i izracunam
    // onda vidim gdje je taj lokalni kartezije u trenutku deltaT, i da li se 
    // slaze TRENUTNA tocka (x,y,z) di je sletio hitac, s onom kako sam izracunao
};

// da li mi treba Local3D koji za parenta ima HardSphereRotatingFrameToSpherical?
// lokalni sustav, baziran na TOCNO ODREDJENOJ TOCKI SFERE, s x, y i z
// za njega NE TREBA davati lat, long i h jer vec ima, a x, y i z transformira lokalno


// TODO 0.8 - zanemari sunce, i samo lokalni rotacijski sustav zemlje

void Demo_Coord_system()
{
    std::cout << "***********************************************************************" << std::endl;
    std::cout << "****                         COORD SYSTEM                          ****" << std::endl;
    std::cout << "***********************************************************************" << std::endl;
    
    InertialFrame3D SolarSystemCMFrame;
    CircleOrbitingCartesianFrame3DToCartesian EarthCMFrame(SolarSystemCMFrame, 150000, 365*24*3600);    // radijus kruženja (jeinica 10^6 m), period kruženja (s)
    
    // za bilo koji (x,y,z) u EarthCM frameu, mogu dobiti poziciju u SolarSystem frameu za vrijeme T
    
    // za bilo koji (r,theta,phi) u EarthCM frameu, mogu dobiti poziciju u SolarSystem frameu za vrijeme T
    // to me zanima, jer ce mi tako sustav Zemlje davati svoje lokalne koordinate

    // sad kreiramo rotirajuci sustav Zemlje, s osi rotacije pod kutom 23 stupnja prema z-osi sustava CM
    // u trenuttku t=0, os rotacije je nagnuta u y-z ravnini
    VectorN<Real, 3> axis({0, sin(23 * Constants::PI / 180), cos(23 * Constants::PI / 180)}); 
    HardSphereRotatingFrameToSpherical EarthSurfaceFrame(EarthCMFrame, 6.4, 24 * 3600, axis); 

    // pozicija na površini zemlje - izračunati u EarthCM i SolarSystemCM frameu
    VectorN<Real, 3> pos({0, 0, 0});
    double latitude = 45;
    double longitude = 16;
    double height = 0;          // on the surface
    double T = 10;

    // pozicija tocke u lokalnom sustavu EarthCMFrame, ali je u SFERNIM koordinatama!!
    Vector3Spherical pos_EarthCM = EarthSurfaceFrame.GetPositionAtTime(latitude, longitude, height, T);
    // koristeci PARENT referentni sustav, transformiramo to u kartezian kord. u EarthCM frameu
    Vector3Cartesian pos_EarthCM_cart = EarthCMFrame.GetLocalPositionAtTime(pos_EarthCM, T);

    // is sun visible from this position?

}