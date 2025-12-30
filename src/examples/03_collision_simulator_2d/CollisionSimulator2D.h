/******************************************************************************
 * CollisionSimulator2D.h - Self-contained 2D Collision Physics Engine
 * ============================================================================
 * 
 * This header provides a complete 2D collision simulation system:
 * 
 *   - Ball2D: Particle with mass, radius, position, velocity, color
 *   - BoxContainer2D: Rectangular container with wall collision handling
 *   - CollisionSimulator2D: Physics engine with space subdivision & parallel
 *   - ContainerFactory2D: Pre-built configurations for common scenarios
 *   - SimResultsCollSim2D: Simulation results storage
 * 
 * Physics:
 *   - Elastic collisions (conserve momentum and kinetic energy)
 *   - Exact collision time calculation for sub-timestep accuracy
 *   - Space subdivision for O(N) average collision detection
 *   - Parallel execution using std::thread for large simulations
 * 
 * Performance:
 *   - RunTypeExact: O(N²) - accurate for small N
 *   - RunTypeFast: O(N) average - space subdivision
 *   - RunTypeFastMultithread: Parallel space subdivision
 * 
 *****************************************************************************/

#ifndef COLLISION_SIMULATOR_2D_SELF_CONTAINED_H
#define COLLISION_SIMULATOR_2D_SELF_CONTAINED_H

#include "MMLBase.h"
#include "base/Vector.h"
#include "base/VectorTypes.h"
#include "base/Matrix.h"
#include "base/Random.h"
#include "tools/Serializer.h"
#include "tools/Visualizer.h"

#include <thread>
#include <atomic>
#include <string>
#include <vector>
#include <map>
#include <cmath>

namespace Collision2D
{
    using namespace MML;

    /**************************************************************************
     * BALL2D - Particle with physical properties
     **************************************************************************/
    
    struct Ball2D
    {
    private:
        double _mass, _radius;
        std::string _color;
        Pnt2Cart _position;
        Vec2Cart _velocity;

    public:
        Ball2D(double mass, double radius, std::string color, 
               const Pnt2Cart& position, const Vec2Cart& velocity)
            : _mass(mass), _radius(radius), _color(color), 
              _position(position), _velocity(velocity) {}

        double  Mass() const { return _mass; }
        double& Mass() { return _mass; }
        double  Rad()  const { return _radius; }
        double& Rad() { return _radius; }

        Pnt2Cart  Pos() const { return _position; }
        Pnt2Cart& Pos() { return _position; }
        Vec2Cart  V() const { return _velocity; }
        Vec2Cart& V() { return _velocity; }

        std::string Color() const { return _color; }
        std::string& Color() { return _color; }
        
        double Speed() const { return _velocity.NormL2(); }
        double KineticEnergy() const { return 0.5 * _mass * (_velocity * _velocity); }
    };

    /**************************************************************************
     * BOX CONTAINER 2D - Rectangular container with wall reflections
     **************************************************************************/
    
    struct BoxContainer2D
    {
        double _width;
        double _height;
        std::vector<Ball2D> _balls;

        BoxContainer2D() : _width(1000), _height(1000) {}
        BoxContainer2D(double width, double height) : _width(width), _height(height) {}

        void AddBall(const Ball2D& ball) { _balls.push_back(ball); }
        Ball2D& Ball(int i) { return _balls[i]; }
        const Ball2D& Ball(int i) const { return _balls[i]; }
        size_t NumBalls() const { return _balls.size(); }

        // Check and handle wall collisions with elastic reflection
        void CheckAndHandleOutOfBounds(int ballIndex)
        {
            Ball2D& ball = _balls[ballIndex];

            // Left wall
            if (ball.Pos().X() < ball.Rad() && ball.V().X() < 0)
            {
                ball.Pos().X() = ball.Rad() + (ball.Rad() - ball.Pos().X());
                ball.V().X() *= -1;
            }
            // Right wall
            if (ball.Pos().X() > _width - ball.Rad() && ball.V().X() > 0)
            {
                ball.Pos().X() -= (ball.Pos().X() + ball.Rad()) - _width;
                ball.V().X() *= -1;
            }
            // Bottom wall
            if (ball.Pos().Y() < ball.Rad() && ball.V().Y() < 0)
            {
                ball.Pos().Y() = ball.Rad() + (ball.Rad() - ball.Pos().Y());
                ball.V().Y() *= -1;
            }
            // Top wall
            if (ball.Pos().Y() > _height - ball.Rad() && ball.V().Y() > 0)
            {
                ball.Pos().Y() -= (ball.Pos().Y() + ball.Rad()) - _height;
                ball.V().Y() *= -1;
            }
        }

        // Fill matrix with ball indices for each subdivided cell
        void SetNumBallsInSubdividedContainer(int nRows, int nCols, Matrix<Vector<int>>& M)
        {
            double cellWidth = _width / nCols;
            double cellHeight = _height / nRows;

            for (size_t i = 0; i < _balls.size(); i++)
            {
                Ball2D& ball = _balls[i];
                int row = static_cast<int>(ball.Pos().Y() / cellHeight);
                int col = static_cast<int>(ball.Pos().X() / cellWidth);

                if (row >= 0 && row < nRows && col >= 0 && col < nCols)
                    M(row, col).push_back(static_cast<int>(i));
            }
        }
    };

    /**************************************************************************
     * SIMULATION RESULTS
     **************************************************************************/
    
    class SimResultsCollSim2D
    {
        int _numBalls;
    public:
        std::vector<Ball2D> Balls;
        std::vector<std::vector<Pnt2Cart>> BallPosList;
        std::vector<std::vector<Vec2Cart>> BallVelList;

        SimResultsCollSim2D(int numBalls)
            : BallPosList(numBalls), BallVelList(numBalls), _numBalls(numBalls) {}

        SimResultsCollSim2D(const std::vector<Ball2D>& balls)
            : Balls(balls), BallPosList(balls.size()), BallVelList(balls.size()), 
              _numBalls(static_cast<int>(balls.size())) {}

        int NumBalls() const { return _numBalls; }

        const std::vector<Pnt2Cart>& getBallPositions(int indBall) const
        {
            return BallPosList[indBall];
        }

        void CalcAllBallsStatistic(int timeStep, double& minSpeed, double& maxSpeed, 
                                   double& avgSpeed, double& speedDev)
        {
            minSpeed = std::numeric_limits<double>::max();
            maxSpeed = std::numeric_limits<double>::lowest();
            double sumSpeed = 0.0, sumSpeedSq = 0.0;
            int count = 0;

            for (const auto& ballVelList : BallVelList)
            {
                if (timeStep < (int)ballVelList.size())
                {
                    double speed = ballVelList[timeStep].NormL2();
                    minSpeed = std::min(minSpeed, speed);
                    maxSpeed = std::max(maxSpeed, speed);
                    sumSpeed += speed;
                    sumSpeedSq += speed * speed;
                    count++;
                }
            }
            if (count > 0)
            {
                avgSpeed = sumSpeed / count;
                speedDev = sqrt(sumSpeedSq / count - avgSpeed * avgSpeed);
            }
        }
    };

    /**************************************************************************
     * COLLISION SIMULATOR RUN TYPES
     **************************************************************************/
    
    enum CollisionSimulatorRunType
    {
        RunTypeExact,              // O(N²) - check all pairs
        RunTypeFast,               // O(N) avg - space subdivision
        RunTypeFastMultithread     // Parallel space subdivision
    };

    /**************************************************************************
     * COLLISION SIMULATOR 2D - Main physics engine
     **************************************************************************/
    
    class CollisionSimulator2D
    {
        BoxContainer2D _box;
        Matrix<Vector<int>> M;  // Space subdivision grid

    public:
        CollisionSimulator2D() : M(2, 2) {}
        
        CollisionSimulator2D(const BoxContainer2D& box) 
            : _box(box), M(2, 2) {}
        
        CollisionSimulator2D(const BoxContainer2D& box, int nRows, int nCols) 
            : _box(box), M(nRows, nCols) {}

        BoxContainer2D& Container() { return _box; }
        const BoxContainer2D& Container() const { return _box; }

        double DistBalls(int i, int j)
        {
            return _box._balls[i].Pos().Dist(_box._balls[j].Pos());
        }

        bool HasBallsCollided(int i, int j)
        {
            return DistBalls(i, j) < _box._balls[i].Rad() + _box._balls[j].Rad();
        }

        // Elastic collision with exact collision time calculation
        void HandleCollision(int stepNum, int m, int n, double dt)
        {
            Ball2D& ball1 = _box._balls[m];
            Ball2D& ball2 = _box._balls[n];

            // Calculate positions before collision
            Pnt2Cart x10 = ball1.Pos() - ball1.V() * dt;
            Pnt2Cart x20 = ball2.Pos() - ball2.V() * dt;

            Vec2Cart dx0(x10, x20);
            Vec2Cart dv(ball2.V() - ball1.V());

            // Solve quadratic for exact collision time
            double A = dv * dv;
            double B = 2 * dx0 * dv;
            double C = dx0 * dx0 - POW2(ball2.Rad() + ball1.Rad());

            double discriminant = B * B - 4 * A * C;
            if (discriminant < 0 || A == 0) return;  // No collision (shouldn't happen)

            double t1 = (-B + sqrt(discriminant)) / (2 * A);
            double t2 = (-B - sqrt(discriminant)) / (2 * A);
            double tCollision = std::min(t1, t2);

            // Positions at collision
            Pnt2Cart x1 = ball1.Pos() + (tCollision - dt) * ball1.V();
            Pnt2Cart x2 = ball2.Pos() + (tCollision - dt) * ball2.V();

            // Elastic collision formula (Wikipedia)
            double m1 = ball1.Mass(), m2 = ball2.Mass();
            Vec2Cart v1 = ball1.V(), v2 = ball2.V();
            Vec2Cart v1_v2 = v1 - v2;
            Vec2Cart x1_x2(x2, x1);

            double normSq = POW2(x1_x2.NormL2());
            if (normSq == 0) return;

            Vec2Cart v1_new = v1 - 2 * m2 / (m1 + m2) * (v1_v2 * x1_x2) / normSq * Vec2Cart(x2, x1);
            Vec2Cart v2_new = v2 - 2 * m1 / (m1 + m2) * (v1_v2 * x1_x2) / normSq * Vec2Cart(x1, x2);

            ball1.V() = v1_new;
            ball2.V() = v2_new;

            // Update positions for remaining time
            ball1.Pos() = x1 + ball1.V() * (dt - tCollision);
            ball2.Pos() = x2 + ball2.V() * (dt - tCollision);
        }

        // O(N²) exact method
        void SimulateOneStepExact(double dt, int stepNum)
        {
            int NumBalls = static_cast<int>(_box._balls.size());

            // Update positions
            for (int i = 0; i < NumBalls; i++)
            {
                _box._balls[i].Pos() = _box._balls[i].Pos() + _box._balls[i].V() * dt;
                _box.CheckAndHandleOutOfBounds(i);
            }

            // Check all pairs
            for (int m = 0; m < NumBalls - 1; m++)
            {
                for (int n = m + 1; n < NumBalls; n++)
                {
                    if (HasBallsCollided(m, n))
                        HandleCollision(stepNum, m, n, dt);
                }
            }
        }

        // O(N) average with space subdivision
        void SimulateOneStepFast(double dt, int stepNum)
        {
            int NumBalls = static_cast<int>(_box._balls.size());

            // Update positions
            for (int i = 0; i < NumBalls; i++)
                _box._balls[i].Pos() = _box._balls[i].Pos() + _box._balls[i].V() * dt;

            // Build spatial grid
            _box.SetNumBallsInSubdividedContainer(M.RowNum(), M.ColNum(), M);

            // Handle wall collisions
            for (int i = 0; i < NumBalls; i++)
                _box.CheckAndHandleOutOfBounds(i);

            // Check collisions within each cell
            for (int i = 0; i < M.RowNum(); i++)
            {
                for (int j = 0; j < M.ColNum(); j++)
                {
                    if (M(i, j).size() > 1)
                    {
                        for (size_t m = 0; m < M(i, j).size() - 1; m++)
                        {
                            for (size_t n = m + 1; n < M(i, j).size(); n++)
                            {
                                int idx1 = M(i, j)[m];
                                int idx2 = M(i, j)[n];
                                if (HasBallsCollided(idx1, idx2))
                                    HandleCollision(stepNum, idx1, idx2, dt);
                            }
                        }
                        M(i, j).Clear();
                    }
                }
            }
        }

        // Parallel space subdivision
        void SimulateOneStepFastMultithread(double dt, int stepNum)
        {
            int NumBalls = static_cast<int>(_box._balls.size());

            // Update positions and handle walls
            for (int i = 0; i < NumBalls; i++)
            {
                _box._balls[i].Pos() = _box._balls[i].Pos() + _box._balls[i].V() * dt;
                _box.CheckAndHandleOutOfBounds(i);
            }

            // Build spatial grid
            _box.SetNumBallsInSubdividedContainer(M.RowNum(), M.ColNum(), M);

            // Parallel collision detection per cell
            std::vector<std::thread> threads;
            for (int i = 0; i < M.RowNum(); i++)
            {
                for (int j = 0; j < M.ColNum(); j++)
                {
                    if (M(i, j).size() > 1)
                    {
                        threads.emplace_back([this, i, j, stepNum, dt]()
                        {
                            for (size_t m = 0; m < M(i, j).size() - 1; m++)
                            {
                                for (size_t n = m + 1; n < M(i, j).size(); n++)
                                {
                                    int idx1 = M(i, j)[m];
                                    int idx2 = M(i, j)[n];
                                    if (HasBallsCollided(idx1, idx2))
                                        HandleCollision(stepNum, idx1, idx2, dt);
                                }
                            }
                            M(i, j).Clear();
                        });
                    }
                }
            }
            for (auto& t : threads) t.join();
        }

        // Main simulation loop
        SimResultsCollSim2D Simulate(int numSteps, double timeStep, 
                                     CollisionSimulatorRunType runType = RunTypeFastMultithread,
                                     bool verbose = false)
        {
            int numBalls = static_cast<int>(_box._balls.size());
            SimResultsCollSim2D results(_box._balls);

            for (int i = 0; i < numSteps; i++)
            {
                if (verbose && i % 10 == 0)
                    std::cout << "Step " << i << "/" << numSteps << "\r" << std::flush;

                // Record state
                for (int j = 0; j < numBalls; j++)
                {
                    results.BallPosList[j].push_back(_box._balls[j].Pos());
                    results.BallVelList[j].push_back(_box._balls[j].V());
                }

                // Advance physics
                switch (runType)
                {
                case RunTypeExact:
                    SimulateOneStepExact(timeStep, i);
                    break;
                case RunTypeFast:
                    SimulateOneStepFast(timeStep, i);
                    break;
                case RunTypeFastMultithread:
                    SimulateOneStepFastMultithread(timeStep, i);
                    break;
                }
            }

            if (verbose)
                std::cout << "\nSimulation complete!\n";

            return results;
        }

        // Serialize simulation for visualization
        bool Serialize(const std::string& fileName, const SimResultsCollSim2D& simResults, 
                       Real dT, int saveEveryNSteps = 1)
        {
            int numBalls = static_cast<int>(_box._balls.size());
            std::vector<std::string> ballColors;
            std::vector<Real> ballRadius;

            for (const auto& ball : _box._balls)
            {
                ballColors.push_back(ball.Color());
                ballRadius.push_back(ball.Rad());
            }

            auto result = Serializer::SaveParticleSimulation2D(
                fileName, numBalls, _box._width, _box._height,
                simResults.BallPosList, ballColors, ballRadius, dT, saveEveryNSteps);
            
            return result.success;
        }
    };

    /**************************************************************************
     * CONTAINER FACTORY - Pre-built configurations
     **************************************************************************/
    
    class ContainerFactory2D
    {
    public:
        // Two groups of balls, each in one half of the container
        // Perfect for mixing/diffusion visualization
        static BoxContainer2D CreateTwoHalves(
            double width, double height,
            const std::string& color1, int nBall1, double mass1, double rad1, double vel1,
            const std::string& color2, int nBall2, double mass2, double rad2, double vel2)
        {
            BoxContainer2D box(width, height);

            // Left half
            for (int i = 0; i < nBall1; i++)
            {
                double x = Random::UniformReal(rad1, width / 2 - rad1);
                double y = Random::UniformReal(rad1, height - rad1);
                Real vx, vy;
                Random::UniformVecDirection2(vx, vy, vel1);
                box.AddBall(Ball2D(mass1, rad1, color1, Pnt2Cart(x, y), Vec2Cart(vx, vy)));
            }

            // Right half
            for (int i = 0; i < nBall2; i++)
            {
                double x = Random::UniformReal(width / 2 + rad2, width - rad2);
                double y = Random::UniformReal(rad2, height - rad2);
                Real vx, vy;
                Random::UniformVecDirection2(vx, vy, vel2);
                box.AddBall(Ball2D(mass2, rad2, color2, Pnt2Cart(x, y), Vec2Cart(vx, vy)));
            }

            return box;
        }

        // Shock wave configuration: energetic core + cold surrounding gas
        static BoxContainer2D CreateShockWave(
            double width, double height,
            int numSurrounding, double mass1, double rad1, double vel1,
            int numEnergetic, double coreRadius, double mass2, double rad2, double vel2)
        {
            BoxContainer2D box(width, height);

            double cx = width / 2.0;
            double cy = height / 2.0;

            // Surrounding cold gas (avoid center)
            for (int i = 0; i < numSurrounding; i++)
            {
                Real vx, vy;
                Random::UniformVecDirection2(vx, vy, vel1);

                Pnt2Cart position(
                    Random::UniformReal(rad1, width - rad1),
                    Random::UniformReal(rad1, height - rad1)
                );

                // Keep away from energetic core
                while (Pnt2Cart(cx, cy).Dist(position) < coreRadius + rad1 + rad2)
                {
                    position = Pnt2Cart(
                        Random::UniformReal(rad1, width - rad1),
                        Random::UniformReal(rad1, height - rad1)
                    );
                }

                box.AddBall(Ball2D(mass1, rad1, "Blue", position, Vec2Cart(vx, vy)));
            }

            // Energetic core at center
            for (int i = 0; i < numEnergetic; i++)
            {
                double angle = Random::UniformReal(0, 2 * Constants::PI);
                double r = coreRadius * Random::UniformReal(0, 1);
                Pnt2Cart position(cx + r * cos(angle), cy + r * sin(angle));

                Real vx, vy;
                Random::UniformVecDirection2(vx, vy, vel2);

                box.AddBall(Ball2D(mass2, rad2, "Red", position, Vec2Cart(vx, vy)));
            }

            return box;
        }

        // N identical balls with random positions and velocities
        static BoxContainer2D CreateUniform(
            double width, double height, int N,
            double mass, double radius, double velocityRange, const std::string& color)
        {
            BoxContainer2D box(width, height);
            
            for (int i = 0; i < N; i++)
            {
                double x = Random::UniformReal(radius, width - radius);
                double y = Random::UniformReal(radius, height - radius);
                Real vx, vy;
                Random::UniformVecDirection2(vx, vy, velocityRange);
                box.AddBall(Ball2D(mass, radius, color, Pnt2Cart(x, y), Vec2Cart(vx, vy)));
            }
            
            return box;
        }

        // Brownian motion: one large particle surrounded by many small ones
        static BoxContainer2D CreateBrownianMotion(
            double width, double height, int numSmallBalls,
            double bigMass, double bigRadius,
            double smallMass, double smallRadius, double smallSpeed)
        {
            BoxContainer2D box(width, height);

            // Central large particle (stationary)
            box.AddBall(Ball2D(bigMass, bigRadius, "Red", 
                              Pnt2Cart(width / 2, height / 2), Vec2Cart(0, 0)));

            // Many small particles
            for (int i = 0; i < numSmallBalls; i++)
            {
                double x = Random::UniformReal(smallRadius, width - smallRadius);
                double y = Random::UniformReal(smallRadius, height - smallRadius);

                // Ensure not overlapping with big particle
                while (box.Ball(0).Pos().Dist(Pnt2Cart(x, y)) < bigRadius + smallRadius + 1)
                {
                    x = Random::UniformReal(smallRadius, width - smallRadius);
                    y = Random::UniformReal(smallRadius, height - smallRadius);
                }

                Real vx, vy;
                Random::UniformVecDirection2(vx, vy, smallSpeed);
                box.AddBall(Ball2D(smallMass, smallRadius, "Blue", 
                                  Pnt2Cart(x, y), Vec2Cart(vx, vy)));
            }

            return box;
        }

        // Grid of evenly spaced particles
        static BoxContainer2D CreateGrid(
            double width, double height, int N,
            double mass, double radius, double speed)
        {
            BoxContainer2D box(width, height);
            double stepX = width / N;
            double stepY = height / N;

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    double x = i * stepX + stepX / 2;
                    double y = j * stepY + stepY / 2;
                    Real vx, vy;
                    Random::UniformVecDirection2(vx, vy, speed);
                    box.AddBall(Ball2D(mass, radius, "Green", Pnt2Cart(x, y), Vec2Cart(vx, vy)));
                }
            }

            return box;
        }
    };

} // namespace Collision2D

#endif // COLLISION_SIMULATOR_2D_SELF_CONTAINED_H
