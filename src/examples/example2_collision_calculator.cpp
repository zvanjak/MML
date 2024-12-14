#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "MMLBase.h"

#include "base/VectorTypes.h"
#endif

using namespace MML;

// Calculating deviation and changes for ever decreasing dT
//   naci onu tocku nakon koje smanjivanje dT ne vodi drugacijim rezultatima

namespace CollisionSimulator
{
	struct Body2D
	{
	private:
		double _mass;
		double _radius;
		Point2Cartesian _position;
		Vector2Cartesian _velocity;

	public:
		Body2D(double mass, double radius, const Point2Cartesian& position, const Vector2Cartesian& velocity)
			: _mass(mass), _radius(radius), _position(position), _velocity(velocity)
		{ }

		double  Mass() const { return _mass; }
		double& Mass() { return _mass; }

		double  Rad() const { return _radius; }
		double& Rad() { return _radius; }

		Point2Cartesian  Pos() const { return _position; }
		Point2Cartesian& Pos() { return _position; }

		Vector2Cartesian  V() const { return _velocity; }
		Vector2Cartesian& V() { return _velocity; }
	};

	struct Container2D
	{
		double _width;
		double _height;

		std::vector<Body2D> _bodies;

		void CheckAndHandleOutOfBounds(int ballIndex)
		{
			const Body2D& ball = ball;

			// ako je izvan boxa i ako ide "u stijenku"
			if (ball.Pos().X() < ball.Rad() && ball.V().X() < 0)
			{
				ball.Pos().X() = ball.Rad() + (ball.Rad() - ball.Pos().X()); // vraćamo ga u box
				ball.V().X() *= -1;
			}

			if (ball.Pos().X() > _width - ball.Rad() && ball.V().X() > 0)
			{
				ball.Pos().X() -= (ball.Pos().X() + ball.Rad()) - _width;
				ball.V().X() *= -1;
			}

			if (ball.Pos().Y() < ball.Rad() && ball.V().Y() < 0)
			{
				ball.Pos().Y() = ball.Rad() + (ball.Rad() - ball.Pos().Y());
				ball.V().Y() *= -1;
			}

			if (ball.Pos().Y() > _height - ball.Rad() && ball.V().Y() > 0)
			{
				ball.Pos().Y() -= (ball.Pos().Y() + ball.Rad()) - _height;
				ball.V().Y() *= -1;
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
			return _box._bodies[i].Pos().Dist(_box._bodies[j].Pos());
		}

		bool HasBallsCollided(int i, int j)
		{
			// ako je udaljenost izmedju njihovih centara manja od zbroja radijusa
			if (DistBalls(i, j) < _box._bodies[i].Rad() + _box._bodies[j].Rad())
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
				_box._bodies[i].Pos() = _box._bodies[i].Pos() + _box._bodies[i].V() * dt;

				_box.CheckAndHandleOutOfBounds(i);
			}

			// provjeriti za sve parove da li su se sudarili
			for (int m = 0; m < NumBalls - 1; m++) {
				for (int n = m + 1; n < NumBalls; n++)
				{
					if (HasBallsCollided(m, n))
					{
						Body2D& ball1 = _box._bodies[m];
						Body2D& ball2 = _box._bodies[n];

						// calculating point where they were before collision
						// (if there was collision with the box wall, then calc.pos. will be outside the box
						// but it doesn't matter, since we need only direction, ie. velocity, to calculate exact collision point)
						Pnt2Cart x10 = ball1.Pos() - ball1.V() * dt;
						Pnt2Cart x20 = ball2.Pos() - ball2.V() * dt;

						Vec2Cart dx0(x10, x20);
						Vec2Cart dv(ball2.V() - ball1.V());

						double A = dv * dv;
						double B = 2 * dx0 * dv;
						double C = dx0 * dx0 - POW2(ball2.Rad() + ball1.Rad());

						double t1 = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
						double t2 = (-B - sqrt(B * B - 4 * A * C)) / (2 * A);

						double tCollision = t1 < t2 ? t1 : t2;
						//double tReverseBalls = tCollision - dt;    // with respect to CURRENT balls position

						// calculating position of balls at the point of collision (moving them backwards)
						Pnt2Cart x1 = ball1.Pos() + (tCollision - dt) * ball1.V();
						Pnt2Cart x2 = ball2.Pos() + (tCollision - dt) * ball2.V();

						// https://en.wikipedia.org/wiki/Elastic_collision - calculating new velocities after collision
						double   m1 = ball1.Mass(), m2 = ball2.Mass();

						Vec2Cart v1(ball1.V()), v2(ball2.V());

						Vec2Cart v1_v2(v1 - v2);
						Vec2Cart x1_x2(x2, x1);

						Vec2Cart v1_new = v1 - 2 * m2 / (m1 + m2) * (v1_v2 * x1_x2) / POW2(x1_x2.NormL2()) * Vec2Cart(x2, x1);
						Vec2Cart v2_new = v2 - 2 * m1 / (m1 + m2) * (v1_v2 * x1_x2) / POW2(x1_x2.NormL2()) * Vec2Cart(x1, x2);

						ball1.V() = v1_new;
						ball2.V() = v2_new;

						// adjusting new ball positions
						ball1.Pos() = x1 + ball1.V() * (dt - tCollision);
						ball2.Pos() = x2 + ball2.V() * (dt - tCollision);
					}
				}
			}
		}
	};
}

void Example2_collision_calculator()
{
	std::cout << "***********************************************************************" << std::endl;
	std::cout << "****                 EXAMPLE 2 - collision calculator              ****" << std::endl;
	std::cout << "***********************************************************************" << std::endl;

	// napraviti 2D prikaz u konzoli s rubovima i kuglama kao * koje se micu!!!
}