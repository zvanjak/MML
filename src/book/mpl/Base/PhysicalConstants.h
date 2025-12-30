#ifndef MPL_PHYSICS_CONSTANTS_H
#define MPL_PHYSICS_CONSTANTS_H

#include "MMLBase.h"

namespace MPL
{
  // https://physics.nist.gov/cuu/Constants/

  // klasifikacija
  //  - exact
  //  - measured
  //  - computed - za njih definirati i kako su izracunate
  //      ove computed kao constexpr definirati!!!

  class PhyConst
  {
  public:
    // FUNDAMENTAL
    static constexpr double GravityConstant = 6.67430e-11;        // N m^2 / kg^2
    static inline constexpr double PlanckConstant = 6.62607015e-34;      // J s 
    static inline constexpr double SpeedOfLight = 299792458;             // m / s

    // permitivity, el. mag.
    static inline constexpr double VacuumElectricPermitivity = 8.8541878128e-12;   // C^2 / ( N m^2 )
    static inline constexpr double VaccumMagneticPermeability = 1.25663706212e-6;   // T m / A

    // alpha - fine structure
    static inline constexpr double FineStructureConstant = 7.2973525693e-3;
    static inline constexpr double InverseFineStructureConstant = 137.035999084;

    // basic charge
    static inline constexpr double CoulumbConstant = 8.9875517923e9;       // kg m^3 / ( C^2 s^2 )
    static inline constexpr double ElementaryCharge = 1.602176634e-19;      // C

    // electron rest mass, proton, neutron
    static inline constexpr double ElectronRestMass = 9.1093837015e-31;     // kg
    static inline constexpr double ProtonRestMass = 1.67262192369e-27;    // kg
    static inline constexpr double NeutronRestMass = 1.67492749804e-27;    // kg

    static inline constexpr double AvogardoConstant = 6.02214076e23;       // mol^-1

    static inline constexpr double BoltzmannConstant = 1.380649e-23;        // J / K

    static inline constexpr double MolarGasConstant = 8.314462618;         // J / ( mol K )

    static inline constexpr double StefanBoltzmannConstant = 5.670374419e-8;     // W / ( m^2 K^4 )

    // TODO - Raydberg
  };

} // namespace MPL

#endif // MPL_PHYSICS_CONSTANTS_H