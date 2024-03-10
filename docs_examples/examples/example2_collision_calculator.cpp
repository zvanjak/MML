#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/VectorTypes.h"
#endif

using namespace MML;


struct Body2D
{
	double _mass;
	double _radius;
	Point2Cartesian _position;
	Vector2Cartesian _velocity;

	Point2Cartesian  Position() const { return _position; }
	Point2Cartesian& Position() { return _position; }

	Vector2Cartesian  Velocity() const { return _velocity; }
	Vector2Cartesian& Velocity() { return _velocity; }
};

struct Container2D
{
	double _width;
	double _height;

	std::vector<Body2D> _bodies;

	void CheckAndHandleOutOfBounds(int ballIndex)
	{
		// ako je izvan boxa i ako ide "u stijenku"
		if (_bodies[ballIndex].Position().X() < _bodies[ballIndex]._radius && _bodies[ballIndex].Velocity().X() < 0)
		{
			_bodies[ballIndex].Position().X() = _bodies[ballIndex]._radius + (_bodies[ballIndex]._radius - _bodies[ballIndex].Position().X()); // vraćamo ga u box
			_bodies[ballIndex].Velocity().X() *= -1;
		}

		if (_bodies[ballIndex].Position().X() > _width - _bodies[ballIndex]._radius && _bodies[ballIndex].Velocity().X() > 0)
		{
			_bodies[ballIndex].Position().X() -= (_bodies[ballIndex].Position().X() + _bodies[ballIndex]._radius) - _width;
			_bodies[ballIndex].Velocity().X() *= -1;
		}

		if (_bodies[ballIndex].Position().Y() < _bodies[ballIndex]._radius && _bodies[ballIndex].Velocity().Y() < 0)
		{
			_bodies[ballIndex].Position().Y() = _bodies[ballIndex]._radius + (_bodies[ballIndex]._radius - _bodies[ballIndex].Position().Y());
			_bodies[ballIndex].Velocity().Y() *= -1;
		}

		if (_bodies[ballIndex].Position().Y() > _height - _bodies[ballIndex]._radius && _bodies[ballIndex].Velocity().Y() > 0)
		{
			_bodies[ballIndex].Position().Y() -= (_bodies[ballIndex].Position().Y() + _bodies[ballIndex]._radius) - _height;
			_bodies[ballIndex].Velocity().Y() *= -1;
		}
	}
};

class CollisionSimulator2D
{
	Container2D _box;

public:
	CollisionSimulator2D() { }

	double DistBalls(int i, int j)
	{
		return _box._bodies[i]._position.Dist(_box._bodies[j]._position);
	}

	bool HasBallsCollided(int i, int j)
	{
		// ako je udaljenost izmedju njihovih centara manja od zbroja radijusa
		if (DistBalls(i, j) < _box._bodies[i]._radius + _box._bodies[j]._radius)
			return true;
		else
			return false;
	}

	void SimulateOneStep(double dt) 
	{
		int NumBalls = _box._bodies.size();

		for (int i = 0; i < NumBalls; i++)
		{
			// update pozicije kugle
			_box._bodies[i].Position() = _box._bodies[i].Position() + _box._bodies[i].Velocity() * dt;

			_box.CheckAndHandleOutOfBounds(i);
		}

		// provjeriti za sve parove da li su se sudarili
		for (int m = 0; m < NumBalls - 1; m++)
		{
			for (int n = m + 1; n < NumBalls; n++)
			{
				if (HasBallsCollided(m, n))
				{
					// računamo točku u kojoj su bili prije sudara
					Point2Cartesian x10 = _box._bodies[m].Position() - _box._bodies[m].Velocity() * dt;
					Point2Cartesian x20 = _box._bodies[n].Position() - _box._bodies[n].Velocity() * dt;

					Vector2Cartesian dx0(x10, x20);

					Vector2Cartesian dv = _box._bodies[n].Velocity() - _box._bodies[m].Velocity();

					double A = dv.ScalarProductCartesian(dv);
					double B = 2 * dx0.ScalarProductCartesian(dv);
					double C = dx0.ScalarProductCartesian(dx0) - (_box._bodies[n]._radius + _box._bodies[m]._radius) * (_box._bodies[n]._radius + _box._bodies[m]._radius);

					double t1 = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
					double t2 = (-B - sqrt(B * B - 4 * A * C)) / (2 * A);

					double tCollision = 0.0;
					if (t1 < t2)
						tCollision = t1;
					else
						tCollision = t2;
					double tReverseBalls = tCollision - dt;    // with respect to CURRENT balls position

					// moving both balls to point of collision
					_box._bodies[m].Position() = _box._bodies[m].Position() + tReverseBalls * _box._bodies[m].Velocity();
					_box._bodies[n].Position() = _box._bodies[n].Position() + tReverseBalls * _box._bodies[n].Velocity();

					// ovo mora biti jednako zbroju radijusa kugli
					double r = DistBalls(m, n);

					// COLLISION - algoritam https://en.wikipedia.org/wiki/Elastic_collision
					// za dvije kugle, točno u "sudarnoj točki", izračunavamo nove brzine nakon sudara
					Vector2Cartesian v1_v2 = _box._bodies[m].Velocity() - _box._bodies[n].Velocity();
					Vector2Cartesian x1_x2(_box._bodies[n].Position(), _box._bodies[m].Position());

					double scalProd = v1_v2.ScalarProductCartesian(x1_x2);
					double denom = x1_x2.NormL2() * x1_x2.NormL2();

					double factor1 = 2 * _box._bodies[n]._mass / (_box._bodies[m]._mass + _box._bodies[n]._mass) * scalProd / denom;
					double factor2 = 2 * _box._bodies[m]._mass / (_box._bodies[m]._mass + _box._bodies[n]._mass) * scalProd / denom;

					Vector2Cartesian v1_new = _box._bodies[m].Velocity() - factor1 * Vector2Cartesian(_box._bodies[n].Position(), _box._bodies[m].Position());
					Vector2Cartesian v2_new = _box._bodies[n].Velocity() - factor2 * Vector2Cartesian(_box._bodies[m].Position(), _box._bodies[n].Position());

					_box._bodies[m].Velocity() = v1_new;
					_box._bodies[n].Velocity() = v2_new;

					// treba i obje kugle pomaknuti u novom smjeru!! za preostali dT od trenutka sudara
					_box._bodies[m].Position() = _box._bodies[m].Position() + _box._bodies[m].Velocity() * (dt - tCollision);
					_box._bodies[n].Position() = _box._bodies[n].Position() + _box._bodies[n].Velocity() * (dt - tCollision);
				}
			}
		}
	}
};


struct Body3D
{
	double _mass;
	double _radius;
	Point3Cartesian _position;
	Vector3Cartesian _velocity;
};

class CollisionSimulator3D
{
public:
	CollisionSimulator3D(const std::vector<Body3D>& bodies)
		: _bodies(bodies)
	{
	}

	std::vector<Body3D> _bodies;

	// HasBodiesCollided
	// SimulateOneStep(dT)
};

void Example2_collision_calculator()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                 EXAMPLE 2 - collision calculator              ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// napraviti 2D prikaz u konzoli s rubovima i kuglama kao * koje se micu!!!
}