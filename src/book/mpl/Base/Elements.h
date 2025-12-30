#ifndef MPL_ELEMENTS_H
#define MPL_ELEMENTS_H

#include <map>
#include <string>

namespace MPL
{
  enum ElementSymbol {
    Empty = 0,
    H = 1,
    He = 2,
    Li = 3,
    Be = 4,
    B = 5,
    C = 6,
    N = 7,
    O = 8,
    F = 9,
    Ne = 10,
    Na = 11,
    Mg = 12,
    Al = 13,
    Si = 14,
    P = 15,
    S = 16,
    Cl = 17,
    Ar = 18,
    K = 19,
    Ca = 20,
    Sc = 21,
    Ti = 22,
    V = 23,
    Cr = 24,
    Mn = 25
  };

  class ElementData
  {
  public:
    const int _atomicNumber;
    const ElementSymbol _symbol;
    const std::string _symbolStr;
    const std::string _name;
    const double _atomicMass;
    const int _neutronNumber;
    const std::string _electronConfigurationShort;
    const std::string _electronConfiguration;

    ElementData() : _atomicNumber(0), _symbol(ElementSymbol::Empty), _symbolStr(""), _name(""), _atomicMass(0.0), _neutronNumber(0), _electronConfigurationShort(""), _electronConfiguration("") {}
    ElementData(int atomicNumber, ElementSymbol symbol, std::string symbolStr, std::string name, double atomicMass, int neutronNumber, std::string electronConfigurationShort, std::string electronConfiguration) :
      _atomicNumber(atomicNumber), _symbol(symbol), _symbolStr(symbolStr), _name(name), _atomicMass(atomicMass), _neutronNumber(neutronNumber), _electronConfigurationShort(electronConfigurationShort), _electronConfiguration(electronConfiguration) {
    }
  };

  enum MatterState {
    Solid = 0,
    Liquid = 1,
    Gas = 2,
    Plasma = 3
	};
  enum BondingType {
    Ionic = 0,
    Covalent = 1,
    Metallic = 2,
    VanDerWaals = 3,
    HydrogenBonding = 4
  };
  enum LatticeType {
    SimpleCubic = 0,
    BodyCenteredCubic = 1,
    FaceCenteredCubic = 2,
    HexagonalClosePacked = 3
	};


	class ElementsDataExtended
  {
  public:
    const int _atomicNumber;
    const ElementSymbol _symbol;

    double _atomicRadius;
    double _covalentRadius;
    double _vanDerWaalsRadius;

    double _density_at_STP;

		MatterState _matterState_at_STP;
		BondingType _bondingType;
		LatticeType _latticeType;

    double _meltingPoint;
    double _boilingPoint;

    double _ionizationEnergy;

    double _heatCapacity;
    double _molarHeatCapacity;

    // TODO - izotopi

    //double _thermalConductivity;
    //double _coefficientOfLinearThermalExpansion;
    //double _electricalResistivity;
    //double _youngModulus;
    //double _shearModulus;
    //double _bulkModulus;
    //double _magneticSusceptibility;
    //double _electronAffinity;
    //double _oxidationStates;
    //double _electronicConfiguration;
  };

  class Elements {
  private:
    const static inline ElementDataExt _emptyElementData[] = {
        ElementDataExt(0, ElementSymbol::Empty, 0.0,  0.0,  0.0,  0.0, MatterState::Solid, BondingType::Ionic,      LatticeType::SimpleCubic, 0.0, 0.0, 0.0, 0.0),
        ElementDataExt(1, ElementSymbol::H,     0.53, 0.31, 1.2,  0.0, MatterState::Gas,   BondingType::Covalent,   LatticeType::SimpleCubic,  14.01, 20.28, 1312.0, 28.836), // Hydrogen
        ElementDataExt(2, ElementSymbol::He,    0.31, 0.28, 1.4,  0.0, MatterState::Gas,   BondingType::VanDerWaals,LatticeType::SimpleCubic,   0.95,  4.22, 2372.0, 20.786), // Helium
        ElementDataExt(3, ElementSymbol::Li,    1.52, 1.23, 1.82, 0.0, MatterState::Solid, BondingType::Ionic,      LatticeType::BodyCenteredCubic,     453.69, 1560.0, 520.2, 24.860), // Lithium
        ElementDataExt(4, ElementSymbol::Be,    1.12, 0.96, 1.53, 0.0, MatterState::Solid, BondingType::Covalent,   LatticeType::HexagonalClosePacked, 1560.0, 2742.0, 899.5, 16.443), // Beryllium
        ElementDataExt(5, ElementSymbol::B,     0.82, 0.82, 1.92, 0.0, MatterState::Solid, BondingType::Covalent,   LatticeType::HexagonalClosePacked, 2075.0, 4000.0, 800.6, 11.087) // Boron
    };

    const static inline ElementData _elementsData[] = {
        ElementData(0,  ElementSymbol::Empty, "", "", 0.0, 0, "", ""),
        ElementData(1,  ElementSymbol::H,  "H",  "Hydrogen",    1.008,      0,  "1s1",          "1s1"),
        ElementData(2,  ElementSymbol::He, "He", "Helium",      4.002602,   2,  "1s2",          "1s2"),
        ElementData(3,  ElementSymbol::Li, "Li", "Lithium",     6.94,       4,  "[He] 2s1",     "1s2, 2s1"),
        ElementData(4,  ElementSymbol::Be, "Be", "Beryllium",   9.0121831,  5,  "[He] 2s2",     "1s2 2s2"),
        ElementData(5,  ElementSymbol::B,  "B",  "Boron",       10.81,      6,  "[He] 2s2 2p1", "1s2 2s2 2p1"),
        ElementData(6,  ElementSymbol::C,  "C",  "Carbon",      12.011,     6,  "[He] 2s2 2p2", "1s2 2s2 2p2"),
        ElementData(7,  ElementSymbol::N,  "N",  "Nitrogen",    14.007,     7,  "[He] 2s2 2p3", "1s2 2s2 2p3"),
        ElementData(8,  ElementSymbol::O,  "O",  "Oxygen",      15.999,     8,  "[He] 2s2 2p4", "1s2 2s2 2p4"),
        ElementData(9,  ElementSymbol::F,  "F",  "Fluorine",    18.9984031, 10, "[He] 2s2 2p5", "1s2 2s2 2p5"),
        ElementData(10, ElementSymbol::Ne, "Ne", "Neon",        20.1797,    10, "[He] 2s2 2p6", "1s2 2s2 2p6"),
        ElementData(11, ElementSymbol::Na, "Na", "Sodium",      22.9897692, 12, "[Ne] 3s1",     "1s2 2s2 2p6 3s1"),
        ElementData(12, ElementSymbol::Mg, "Mg", "Magnesium",   24.305,     12, "[Ne] 3s2",     "1s2 2s2 2p6 3s2"),
        ElementData(13, ElementSymbol::Al, "Al", "Aluminium",   26.9815385, 14, "[Ne] 3s2 3p1", "1s2 2s2 2p6 3s2 3p1"),
        ElementData(14, ElementSymbol::Si, "Si", "Silicon",     28.085,     14, "[Ne] 3s2 3p2", "1s2 2s2 2p6 3s2 3p2"),
        ElementData(15, ElementSymbol::P,  "P",  "Phosphorus",  30.9737619, 16, "[Ne] 3s2 3p3", "1s2 2s2 2p6 3s2 3p3"),
        ElementData(16, ElementSymbol::S,  "S",  "Sulfur",      32.06,      16, "[Ne] 3s2 3p4", "1s2 2s2 2p6 3s2 3p4"),
        ElementData(17, ElementSymbol::Cl, "Cl", "Chlorine",    35.45,      18, "[Ne] 3s2 3p5", "1s2 2s2 2p6 3s2 3p5"),
        ElementData(18, ElementSymbol::Ar, "Ar", "Argon",       39.948,     22, "[Ne] 3s2 3p6", "1s2 2s2 2p6 3s2 3p6"),
        ElementData(19, ElementSymbol::K,  "K",  "Potassium",   39.0983,    20, "[Ar] 4s1",     "1s2 2s2 2p6 3s2 3p6 4s1"),
        ElementData(20, ElementSymbol::Ca, "Ca", "Calcium",     40.078,     20, "[Ar] 4s2",     "1s2 2s2 2p6 3s2 3p6 4s2"),
        ElementData(21, ElementSymbol::Sc, "Sc", "Scandium",    44.955908,  24, "[Ar] 3d1 4s2", "1s2 2s2 2p6 3s2 3p6 3d1 4s2"),
        ElementData(22, ElementSymbol::Ti, "Ti", "Titanium",    47.867,     26, "[Ar] 3d2 4s2", "1s2 2s2 2p6 3s2 3p6 3d2 4s2"),
        ElementData(23, ElementSymbol::V,  "V",  "Vanadium",    50.9415,    28, "[Ar] 3d3 4s2", "1s2 2s2 2p6 3s2 3p6 3d3 4s2"),
        ElementData(24, ElementSymbol::Cr, "Cr", "Chromium",    51.9961,    28, "[Ar] 3d5 4s1", "1s2 2s2 2p6 3s2 3p6 3d5 4s1"),
        ElementData(25, ElementSymbol::Mn, "Mn", "Manganese",   54.938044,  30, "[Ar] 3d5 4s2", "1s2 2s2 2p6 3s2 3p6 3d5 4s2")

    };

    const static inline std::map<ElementSymbol, ElementData> _elements = {
        {ElementSymbol::Empty, _elementsData[0]},
        {ElementSymbol::H, _elementsData[1]}
    };

  public:
    static const ElementData getElem(int atomicNumber) {
      return _elementsData[atomicNumber];
    }
    static const ElementData getElem(const ElementSymbol element) {
      return _elements.at(element);
    }
    static const ElementData getElem(const std::string& name) {
      for (auto& elem : _elements) {
        if (elem.second._name == name) {
          return elem.second;
        }
      }
      throw std::runtime_error("Element with symbol " + name + " not found!");
    }

    static std::string getSymbolStr(ElementSymbol element) {
      return _elements.at(element)._symbolStr;
    }
    static std::string getName(ElementSymbol element) {
      return _elements.at(element)._name;
    }
  };
} // namespace MPL

#endif // MPL_ELEMENTS_H