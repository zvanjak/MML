#if !defined MPL_SCHWARZSCHILD_METRIC_H
#define MPL_SCHWARZSCHILD_METRIC_H

#include "MMLBase.h"


using namespace MML;

namespace MPL
{
  class SchwarzschildMetric : public MetricTensorField<4>
  {
    // Parameters for Schwarzschild black hole
    Real _mass; // Mass of the black hole (in geometric units, G = c = 1)
    Real _G;    // Gravitational constant
    Real _c;    // Speed of light

  public:
    SchwarzschildMetric(Real mass = 1.0, Real G = 1.0, Real c = 1.0)
      : MetricTensorField<4>(2, 0), _mass(mass), _G(G), _c(c) {}

    // Schwarzschild radius
    Real schwarzschild_radius() const { return 2.0 * _G * _mass / (_c * _c); }

    Real getMassInGeometricUnits(Real massKg)
    {
      // Physical constants in SI units
      constexpr Real G = 6.67430e-11;           // m^3 kg^-1 s^-2
      constexpr Real c = 2.99792458e8;          // m/s

      // Conversion formula: M_geom = G * M_SI / c^2
      return G * massKg / (c * c);              // Result in meters      
    }
    // pos: [t, r, theta, phi]
    virtual Real Component(int i, int j, const VectorN<Real, 4>& pos) const override
    {
      const Real t = pos[0];
      const Real r = pos[1];
      const Real theta = pos[2];
      const Real phi = pos[3];

      const Real rs = schwarzschild_radius();

      // Metric signature (-,+,+,+)
      if (i == 0 && j == 0)
      {
        // g_tt = -(1 - rs/r) * c^2
        return -(1.0 - rs / r) * _c * _c;
      }
      else if (i == 1 && j == 1)
      {
        // g_rr = 1 / (1 - rs/r)
        return 1.0 / (1.0 - rs / r);
      }
      else if (i == 2 && j == 2)
      {
        // g_theta_theta = r^2
        return r * r;
      }
      else if (i == 3 && j == 3)
      {
        // g_phi_phi = r^2 * sin^2(theta)
        return r * r * std::sin(theta) * std::sin(theta);
      }
      else
      {
        // Off-diagonal terms are zero
        return 0.0;
      }
    }

    // Accessors for parameters
    Real getMass() const { return _mass; }
    void setMass(Real m) { _mass = m; }
    Real getG() const { return _G; }
    void setG(Real g) { _G = g; }
    Real getC() const { return _c; }
    void setC(Real c) { _c = c; }
  };
}

#endif // MPL_SCHWARZSCHILD_METRIC_H