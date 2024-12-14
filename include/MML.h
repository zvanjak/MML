//  __ __ __ __ _
// |  |  |  |  | |    Minimal Math Library for Modern C++
// | | | | | | | |__  version 1.0
// |_| |_|_| |_|____| https://github.com/zvanjak/mml
//
// Copyright: 2023 - 2024, Zvonimir Vanjak 
//
// LICENSE:
// Code is given as it is, without any warranty. Use it at your own risk.
// STRICTLY NON-COMMERCIAL USE ONLY!
// Unfortunately, also unavailable for Open Source projects, due to restrictive Numerical Recipes license.
// So basically, it is for personal, educational and research use only.

#ifndef MML_SINGLE_HEADER
#define MML_SINGLE_HEADER

#include <exception>
#include <stdexcept>
#include <initializer_list>
#include <algorithm>
#include <memory>
#include <functional>

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <cmath>
#include <limits>
#include <complex>
#include <numbers>


// https://opensource.apple.com/source/CarbonHeaders/CarbonHeaders-18.1/TargetConditionals.h.auto.html
#ifdef __APPLE__
#  include <TargetConditionals.h>
#  if (defined(TARGET_OS_OSX) && TARGET_OS_OSX == 1) || (defined(TARGET_OS_MAC) && TARGET_OS_MAC == 1)
#    define MML_PLATFORM_MAC
#  elif (defined(TARGET_OS_IPHONE) && TARGET_OS_IPHONE == 1)
#    define MML_PLATFORM_IPHONE
#  endif

#elif defined(linux) || defined(__linux) || defined(__linux__)
#  define MML_PLATFORM_LINUX

#elif defined(WIN32) || defined(__WIN32__) || defined(_WIN32) || defined(_MSC_VER) || defined(__MINGW32__)
#  define MML_PLATFORM_WINDOWS
#endif

///////////////////////////   ./include/MMLBase.h   ///////////////////////////





// https://opensource.apple.com/source/CarbonHeaders/CarbonHeaders-18.1/TargetConditionals.h.auto.html



// Complex must have the same underlaying type as Real
typedef double               Real;      // default real type
typedef std::complex<double> Complex;   // default complex type

// Global paths for Visualizers
static const std::string GLOB_PATH_ResultFiles = "E:\\Projects\\MinimalMathLibrary\\results\\";
static const std::string GLOB_PATH_RealFuncViz = "E:\\Projects\\MinimalMathLibrary\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe";
static const std::string GLOB_PATH_SurfaceViz = "E:\\Projects\\MinimalMathLibrary\\tools\\visualizers\\scalar_function_2d_visualizer\\MML_ScalarFunction2Visualizer.exe";
static const std::string GLOB_PATH_ParametricCurveViz = "E:\\Projects\\MinimalMathLibrary\\tools\\visualizers\\parametric_curve_visualizer\\MML_ParametricCurveVisualizer.exe";
static const std::string GLOB_PATH_VectorFieldViz = "E:\\Projects\\MinimalMathLibrary\\tools\\visualizers\\vector_field_visualizer\\MML_VectorFieldVisualizer.exe";

namespace MML
{
	template<class Type>
	static Real Abs(const Type& a)
	{
		return std::abs(a);
	}
	template<class Type>
	static Real Abs(const std::complex<Type>& a)
	{
		return sqrt(a.real() * a.real() + a.imag() * a.imag());
	}

	template<class T> inline T POW2(const T &a) { return ((a) * (a)); }
	template<class T> inline T POW3(const T &a) { return ((a) * (a) * (a)); }
	template<class T> inline T POW4(const T &a) { return ((a) * (a) * (a) * (a)); }
	template<class T> inline T POW5(const T &a) { return ((a) * (a) * (a) * (a) * (a)); }
	template<class T> inline T POW6(const T &a) { return ((a) * (a) * (a) * (a) * (a) * (a)); }

	template<class T>
	inline T SIGN(const T& a, const T& b)
	{
		return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
	}

	////////////                  Constants                ////////////////
	namespace Constants
	{
		static inline const Real PI = std::numbers::pi;
		static inline const Real Epsilon = std::numeric_limits<Real>::epsilon();
		static inline const Real PositiveInf = std::numeric_limits<Real>::max();
		static inline const Real NegativeInf = -std::numeric_limits<Real>::max();
	}

	namespace Defaults
	{
    static int VectorPrintWidth = 15;
    static int VectorPrintPrecision = 10;

		//////////               Default precisions             ///////////
		// TODO - make dependent on Real type (ie. different values for float, double and long double)
		static inline const double ComplexEqualityPrecision = 1e-15;
		static inline const double ComplexAbsEqualityPrecision = 1e-15;
		static inline const double MatrixEqualityPrecision = 1e-15;
		static inline const double VectorEqualityPrecision = 1e-15;

		static inline const double IsMatrixSymmetricPrecision = 1e-15;
		static inline const double IsMatrixDiagonalPrecision = 1e-15;
		static inline const double IsMatrixUnitPrecision = 1e-15;
		static inline const double IsMatrixOrthogonalPrecision = 1e-15;

		static inline const double DerivationDefaultStep = 1e-6;

		static inline const int    IntegrateTrapMaxSteps = 20;
		static inline const double IntegrateTrapEPS = 1.0e-5;

		static inline const int    IntegrateSimpMaxSteps = 20;
		static inline const double IntegrateSimpEPS = 1.0e-5;

		static inline const int    IntegrateRombMaxSteps = 20;
		static inline const double IntegrateRombEPS = 1.0e-6;

		static inline const double WorkIntegralPrecision = 1e-05;
		static inline const double LineIntegralPrecision = 1e-05;
	}
}
///////////////////////////   ./include/MMLExceptions.h   ///////////////////////////

namespace MML
{
	//////////             Vector error exceptions            ///////////
	class VectorInitializationError : public std::invalid_argument
	{
	public:
		int _size1;
		VectorInitializationError(std::string inMessage, int size1) : std::invalid_argument(inMessage), _size1(size1)
		{ }
	};
	class VectorDimensionError : public std::invalid_argument
	{
	public:
		int _size1, _size2;
		VectorDimensionError(std::string inMessage, int size1, int size2) : std::invalid_argument(inMessage), _size1(size1), _size2(size2)
		{ }
	};
	class VectorAccessBoundsError : public std::out_of_range
	{
	public:
		int _i, _n;
		VectorAccessBoundsError(std::string inMessage, int i, int n) : std::out_of_range(inMessage), _i(i), _n(n)
		{ }
	};

	//////////             Matrix error exceptions            ///////////
	class MatrixAllocationError : public std::out_of_range
	{
	public:
		int _rows, _cols;
		MatrixAllocationError(std::string inMessage, int rows, int cols) : std::out_of_range(inMessage), _rows(rows), _cols(cols)
		{ }
	};
	class MatrixAccessBoundsError : public std::out_of_range
	{
	public:
		int _i, _j, _rows, _cols;
		MatrixAccessBoundsError(std::string inMessage, int i, int j, int rows, int cols) : std::out_of_range(inMessage), _i(i), _j(j), _rows(rows), _cols(cols)
		{ }
	};
	class MatrixDimensionError : public std::invalid_argument
	{
	public:
		int _rows1, _cols1, _rows2, _cols2;

		MatrixDimensionError(std::string inMessage, int r1, int c1, int r2, int c2) : std::invalid_argument(inMessage), _rows1(r1), _cols1(c1), _rows2(r2), _cols2(c2)
		{ }
	};
	class SingularMatrixError : public std::domain_error
	{
	public:
		SingularMatrixError(std::string inMessage) : std::domain_error(inMessage)
		{ }
	};

	//////////             Integration exceptions            ///////////
	class IntegrationTooManySteps : public std::domain_error
	{
	public:
		IntegrationTooManySteps(std::string inMessage) : std::domain_error(inMessage)
		{ }
	};

	////////////             Tensor exceptions             /////////////
	class TensorCovarContravarNumError : public std::invalid_argument
	{
	public:
		int _numContra, _numCo;
		TensorCovarContravarNumError(std::string inMessage, int size1, int size2) : std::invalid_argument(inMessage), _numContra(size1), _numCo(size2)
		{ }
	};

	class TensorCovarContravarAirthmeticError : public std::invalid_argument
	{
	public:
		int _numContra, _numCo;
        int _bContra, _bCo;

		TensorCovarContravarAirthmeticError(std::string inMessage, int contra, int co, int b_contra, int b_co) : std::invalid_argument(inMessage), _numContra(contra), _numCo(co), _bContra(b_contra), _bCo(b_co)
		{ }
	};

	class TensorIndexError : public std::invalid_argument
	{
	public:
		TensorIndexError(std::string inMessage) : std::invalid_argument(inMessage)
		{ }
	};    
}

///////////////////////////   ./include/interfaces/IInterval.h   ///////////////////////////

namespace MML
{
	// group
	class IInterval
	{
	public:
		virtual Real getLowerBound() const = 0;
		virtual Real getUpperBound() const = 0;
		virtual Real getLength() const = 0;

		virtual bool isContinuous() const = 0;
		virtual bool contains(Real x) const = 0;

		virtual void GetEquidistantCovering(int numPoints) = 0;

		virtual ~IInterval() {}
	};
}
///////////////////////////   ./include/base/Intervals.h   ///////////////////////////


namespace MML
{
	// TODO - finalize intervals properly (and use them in test beds)
	///////////////////////////////////////////////   Interfaces    ///////////////////////////////////////////
	enum class EndpointType
	{
		OPEN,
		CLOSED,
		NEG_INF,
		POS_INF
	};

	class BaseInterval : public IInterval
	{
	protected:
		Real _lower, _upper;
		EndpointType _lowerType, _upperType;

		BaseInterval(Real lower, EndpointType lowerType, Real upper, EndpointType upperType) : _lower(lower), _lowerType(lowerType), _upper(upper), _upperType(upperType) { }
	public:
		virtual ~BaseInterval() {}

		Real getLowerBound() const { return _lower; }
		Real getUpperBound() const { return _upper; }
		Real getLength()     const { return _upper - _lower; }

		virtual bool isContinuous()  const { return true; }     // we suppose continuous intervals by default

		void GetEquidistantCovering(int numPoints) 
		{ 
			// TODO 1.1 - implement equidistant covering
		}
	};

	class CompleteRInterval : public BaseInterval
	{
	public:
		CompleteRInterval() : BaseInterval(-std::numeric_limits<double>::max(), EndpointType::NEG_INF, std::numeric_limits<double>::max(), EndpointType::POS_INF) { }
		bool contains(Real x) const { return true; }
	};
	class CompleteRWithReccuringPointHoles : public BaseInterval
	{
		Real _hole0, _holeDelta;
	public:
		CompleteRWithReccuringPointHoles(Real hole0, Real holeDelta) : BaseInterval(-std::numeric_limits<double>::max(), EndpointType::NEG_INF, std::numeric_limits<double>::max(), EndpointType::POS_INF)
		{
			_hole0 = hole0;
			_holeDelta = holeDelta;
		}

		bool contains(Real x) const {


			if (x == _hole0)
				return false;

			Real diff = (x - _hole0) / _holeDelta;
			if (diff == (int)diff)
				return false;

			return true;
		}
		bool isContinuous() const { return false; }
	};
	class OpenInterval : public BaseInterval
	{
		// for equidistant covering
		double _lowerRealDif = 0.0000001;
	public:
		OpenInterval(Real lower, Real upper) : BaseInterval(lower, EndpointType::OPEN, upper, EndpointType::OPEN) { }
		bool contains(Real x) const { return (x > _lower) && (x < _upper); }
	};
	class OpenClosedInterval : public BaseInterval
	{
	public:
		OpenClosedInterval(Real lower, Real upper) : BaseInterval(lower, EndpointType::OPEN, upper, EndpointType::CLOSED) { }
		bool contains(Real x) const { return (x > _lower) && (x <= _upper); }
	};
	class ClosedInterval : public BaseInterval
	{
	public:
		ClosedInterval(Real lower, Real upper) : BaseInterval(lower, EndpointType::CLOSED, upper, EndpointType::CLOSED) { }

		bool contains(Real x) const {
			return (x >= _lower) && (x <= _upper);
		}
	};
	class ClosedIntervalWithReccuringPointHoles : public BaseInterval
	{
		Real _hole0, _holeDelta;
	public:
		ClosedIntervalWithReccuringPointHoles(Real lower, Real upper, Real hole0, Real holeDelta)
			: BaseInterval(lower, EndpointType::CLOSED, upper, EndpointType::CLOSED)
		{
			_hole0 = hole0;
			_holeDelta = holeDelta;
		}

		bool contains(Real x) const {
			// check for hole!
			if (x == _hole0)
				return false;

			Real diff = (x - _hole0) / _holeDelta;
			if (diff == (int)diff)
				return false;

			return (x >= _lower) && (x <= _upper);
		}
		bool isContinuous() const { return false; }
	};
	class ClosedOpenInterval : public BaseInterval
	{
	public:
		ClosedOpenInterval(Real lower, Real upper) : BaseInterval(lower, EndpointType::CLOSED, upper, EndpointType::OPEN) { }

		bool contains(Real x) const {
			return (x >= _lower) && (x < _upper);
		}
	};
	class NegInfToOpenInterval : public BaseInterval
	{
	public:
		NegInfToOpenInterval(Real upper) : BaseInterval(-std::numeric_limits<double>::max(), EndpointType::NEG_INF, upper, EndpointType::OPEN) { }

		bool contains(Real x) const {
			return x < _upper;
		}
	};
	class NegInfToClosedInterval : public BaseInterval
	{
	public:
		NegInfToClosedInterval(Real upper) : BaseInterval(-std::numeric_limits<double>::max(), EndpointType::NEG_INF, upper, EndpointType::CLOSED) { }

		bool contains(Real x) const {
			return x <= _upper;
		}
	};
	class OpenToInfInterval : public BaseInterval
	{
	public:
		OpenToInfInterval(Real lower) : BaseInterval(lower, EndpointType::OPEN, std::numeric_limits<double>::max(), EndpointType::POS_INF) { }

		bool contains(Real x) const {
			return x > _lower;
		}
	};
	class ClosedToInfInterval : public BaseInterval
	{
	public:
		ClosedToInfInterval(Real lower) : BaseInterval(lower, EndpointType::CLOSED, std::numeric_limits<double>::max(), EndpointType::POS_INF) { }

		bool contains(Real x) const {
			return x >= _lower;
		}
	};

	class Interval : public IInterval
	{
		Real _lower, _upper;
		std::vector<std::shared_ptr<BaseInterval>> _intervals;
	public:
		Interval() {}
		Interval(std::initializer_list<BaseInterval*> intervals)
		{
			for (BaseInterval* interval : intervals)
			{
				_intervals.emplace_back(std::shared_ptr<BaseInterval>(interval));
			}

			// sortirati po lower bound
			// i redom provjeriti presijecanja
		}

		template<class _IntervalType>
		Interval& AddInterval(const _IntervalType& interval)
		{
			_intervals.emplace_back(std::make_shared<_IntervalType>(interval));
			return *this;
		}

		static Interval Intersection(const IInterval& a, const IInterval& b)
		{
			Interval ret;
			// TODO - implement intersection
			return ret;
		}
		static Interval Difference(const IInterval& a, const IInterval& b)
		{
			Interval ret;
			// TODO - implement diff
			return ret;
		}
		static Interval Complement(const IInterval& a)
		{
			Interval ret;
			// TODO - implement complement
			return ret;
		}
		Real getLowerBound() const { return 0; }
		Real getUpperBound() const { return 0; }
		Real getLength()     const { return 0; }

		bool isContinuous()  const { return false; }
		bool contains(Real x) const
		{
			// check for each interval if it contains x
			for (auto& interval : _intervals)
			{
				if (interval->contains(x))
					return true;
			}
			return false;
		}
		// TODO - implement this
		// bool contains(const IInterval &other) const = 0;
		// bool intersects(const IInterval &other) const = 0;

		void GetEquidistantCovering(int numPoints) { }
	};
}

///////////////////////////   ./include/base/StandardFunctions.h   ///////////////////////////


namespace MML
{
    namespace Functions
    {
        // Functions of REAL domain
        // TODO - umjesto Real, naprosto staviti sve tri dostupne varijante?
        static inline Real Sin(Real x) { return sin(x); }
        static inline Real Cos(Real x) { return cos(x); }
        static inline Real Sec(Real x) { return 1.0 / cos(x); }
        static inline Real Csc(Real x) { return 1.0 / sin(x); }
        static inline Real Tan(Real x) { return tan(x); }
        static inline Real Ctg(Real x) { return 1.0 / tan(x); }
        
        static inline Real Exp(Real x) { return exp(x); }
        static inline Real Log(Real x) { return log(x); }
        static inline Real Log10(Real x){ return log10(x); }
        static inline Real Sqrt(Real x) { return sqrt(x); }
        static inline Real Pow(Real x, Real y) { return pow(x, y); }
        
        static inline Real Sinh(Real x) { return sinh(x); }
        static inline Real Cosh(Real x) { return cosh(x); }
        static inline Real Sech(Real x) { return 1.0 / cosh(x); }
        static inline Real Csch(Real x) { return 1.0 / sinh(x); }        
        static inline Real Tanh(Real x) { return tanh(x); }
        static inline Real Ctgh(Real x) { return 1.0 / tanh(x); }        
        
        static inline Real Asin(Real x) { return asin(x); }
        static inline Real Acos(Real x) { return acos(x); }
        static inline Real Atan(Real x) { return atan(x); }

        static inline Real Asinh(Real x) { return asinh(x); }
        static inline Real Acosh(Real x) { return acosh(x); }
        static inline Real Atanh(Real x) { return atanh(x); }

        static inline Real Erf(Real x)  { return std::erf(x); }
        static inline Real Erfc(Real x) { return std::erfc(x); }

        static inline Real TGamma(Real x) { return std::tgamma(x); }
        static inline Real LGamma(Real x) { return std::lgamma(x); }
        static inline Real RiemannZeta(Real x) { return std::riemann_zeta(x); }
        static inline Real Comp_ellint_1(Real x) { return std::comp_ellint_1(x); }
        static inline Real Comp_ellint_2(Real x) { return std::comp_ellint_2(x); }

        static inline Real Hermite(unsigned int n, Real x) { return std::hermite(n, x); }
        static inline Real Legendre(unsigned int n, Real x) { return std::legendre(n, x); }
        static inline Real Laguerre(unsigned int n, Real x) { return std::laguerre(n, x); }
        static inline Real SphBessel(unsigned int n, Real x) { return std::sph_bessel(n, x); }
        static inline Real SphLegendre(int n1, int n2, Real x) { return std::sph_legendre(n1, n2, x); }
        
        // Functions of COMPLEX domain
        static inline Complex Sin(Complex x) { return sin(x); }
        static inline Complex Cos(Complex x) { return cos(x); }
        static inline Complex Sec(Complex x) { return Real{1.0} / cos(x); }
        static inline Complex Csc(Complex x) { return Real{1.0} / sin(x); }
        static inline Complex Tan(Complex x) { return tan(x); }
        static inline Complex Ctg(Complex x) { return Real{1.0} / tan(x); }
        
        static inline Complex Exp(Complex x) { return exp(x); }
        static inline Complex Log(Complex x) { return log(x); }
        static inline Complex Sqrt(Complex x) { return sqrt(x); }
        static inline Complex Pow(Complex x, Complex y) { return pow(x, y); }
        
        static inline Complex Sinh(Complex x) { return sinh(x); }
        static inline Complex Cosh(Complex x) { return cosh(x); }
        static inline Complex Sech(Complex x) { return Complex(1.0) / cosh(x); }
        static inline Complex Csch(Complex x) { return Complex(1.0) / sinh(x); }        
        static inline Complex Tanh(Complex x) { return tanh(x); }
        static inline Complex Ctgh(Complex x) { return Complex(1.0) / tanh(x); } 

        static inline Complex Asin(Complex x) { return asin(x); }
        static inline Complex Acos(Complex x) { return acos(x); }
        static inline Complex Atan(Complex x) { return atan(x); }

        static inline Complex Asinh(Complex x) { return asinh(x); }
        static inline Complex Acosh(Complex x) { return acosh(x); }
        static inline Complex Atanh(Complex x) { return atanh(x); }    

        static inline Real Factorial(int n) {
            Real fact = 1.0;
            for( int i=2; i<=n; i++)
                fact *= i;
            return fact;
        }   
        static inline long long FactorialInt(int n) {
            long long fact = 1;
            for( int i=2; i<=n; i++)
                fact *= i;
            return fact;
        }           
    }
}

///////////////////////////   ./include/base/Geometry.h   ///////////////////////////


namespace MML
{
	class Point2Cartesian
	{
	private:
		Real _x, _y;

	public:
		Real  X() const { return _x; }
		Real& X() { return _x; }
		Real  Y() const { return _y; }
		Real& Y() { return _y; }

		Point2Cartesian() : _x(0), _y(0) {}
		Point2Cartesian(Real x, Real y) : _x(x), _y(y) {}

		Real Dist(const Point2Cartesian& b) const { return sqrt(POW2(b.X() - X()) + POW2(b.Y() - Y())); }

		bool	operator==(const Point2Cartesian& b) const { return (X() == b.X()) && (Y() == b.Y()); }
		bool	operator!=(const Point2Cartesian& b) const { return (X() != b.X()) || (Y() != b.Y()); }

		Point2Cartesian operator+(const Point2Cartesian& b) const { return Point2Cartesian(X() + b.X(), Y() + b.Y()); }

		friend Point2Cartesian operator*(const Point2Cartesian& a, Real b) { return Point2Cartesian(a.X() * b, a.Y() * b); }
		friend Point2Cartesian operator*(Real a, const Point2Cartesian& b) { return Point2Cartesian(a * b.X(), a * b.Y()); }
		friend Point2Cartesian operator/(const Point2Cartesian& a, Real b) { return Point2Cartesian(a.X() / b, a.Y() / b); }
	};	

	class Point2Polar
	{
	private:
		Real _r, _phi;

	public:
		Real  R() const { return _r; }
		Real& R() { return _r; }
		Real  Phi() const { return _phi; }
		Real& Phi() { return _phi; }

		Point2Polar() : _r(0), _phi(0) {}
		Point2Polar(Real r, Real phi) : _r(r), _phi(phi) {}
		Point2Polar(const Point2Cartesian& pnt) 
		{ 
			_r = sqrt(POW2(pnt.X()) + POW2(pnt.Y()));
			_phi = atan2(pnt.Y(), pnt.X());
		}

		static Point2Polar GetFromCartesian(const Point2Cartesian &pnt)	{ return Point2Polar(pnt); }
		Point2Cartesian TransfToCartesian() const {  return Point2Cartesian(R() * cos(Phi()), R() * sin(Phi())); }

		Real Dist(const Point2Polar& b) const { return sqrt(R() * R() + b.R() * b.R() - 2 * R() * b.R() * cos(b.Phi() - Phi())); }
	};

	class Point3Cartesian
	{
	private:
		Real _x, _y, _z;

	public:
		Real  X() const { return _x; }
		Real& X() { return _x; }
		Real  Y() const { return _y; }
		Real& Y() { return _y; }
		Real  Z() const { return _z; }
		Real& Z() { return _z; }

		Point3Cartesian() : _x(0), _y(0), _z(0) {}
		Point3Cartesian(Real x, Real y, Real z) : _x(x), _y(y), _z(z) {}

		Real Dist(const Point3Cartesian& b) const { return sqrt(POW2(b.X() - X()) + POW2(b.Y() - Y()) + POW2(b.Z() - Z())); }

		bool	operator==(const Point3Cartesian& b) const {return (X() == b.X()) && (Y() == b.Y()) && (Z() == b.Z()); }
		bool	operator!=(const Point3Cartesian& b) const {return (X() != b.X()) || (Y() != b.Y()) || (Z() != b.Z()); }

		Point3Cartesian operator+(const Point3Cartesian& b) const { return Point3Cartesian(X() + b.X(), Y() + b.Y(), Z() + b.Z()); }

		friend Point3Cartesian operator*(const Point3Cartesian& a, Real b) { return Point3Cartesian(a.X() * b, a.Y() * b, a.Z() * b); }
		friend Point3Cartesian operator*(Real a, const Point3Cartesian& b) { return Point3Cartesian(a * b.X(), a * b.Y(), a * b.Z()); }
		friend Point3Cartesian operator/(const Point3Cartesian& a, Real b) { return Point3Cartesian(a.X() / b, a.Y() / b, a.Z() / b); }
	};

	class Triangle
	{
	private:
		Real _a, _b, _c;

	public:
		Real  A() const { return _a; }
		Real& A() { return _a; }
		Real  B() const { return _b; }
		Real& B() { return _b; }
		Real  C() const { return _c; }
		Real& C() { return _c; }

		Triangle() : _a(0.0), _b(0.0), _c(0.0){}
		Triangle(Real a, Real b, Real c) : _a(a), _b(b), _c(c) {}

		Real Area() const
		{
			Real s = (_a + _b + _c) / 2.0;
			return sqrt(s * (s - _a) * (s - _b) * (s - _c));
		}
		bool IsRight() const
		{
			return ( POW2(_a) + POW2(_b) == POW2(_c)) || 
							(POW2(_a) + POW2(_c) == POW2(_b)) || 
							(POW2(_b) + POW2(_c) == POW2(_a) );
		}
		bool IsIsosceles() const
		{
			return (_a == _b) || (_a == _c) || (_b == _c);
		}
		bool IsEquilateral() const
		{
			return (_a == _b) && (_a == _c);
		}
	};

	typedef Point2Cartesian Pnt2Cart;
	typedef Point2Polar			Pnt2Pol;
	typedef Point3Cartesian Pnt3Cart;
}
///////////////////////////   ./include/base/Vector.h   ///////////////////////////

namespace MML
{
	template<class Type>
	class	Vector
	{
	private:
		std::vector<Type> _elems;

	public:
		///////////////////////          Constructors and destructor       //////////////////////
		Vector() {}
		explicit Vector(int n) {
			if(n < 0)
				throw VectorInitializationError("Vector::Vector - negative size", n);

			_elems.resize(n);
		}
		explicit Vector(int n, const Type &val) {
			if (n < 0)
				throw VectorInitializationError("Vector::Vector - negative size", n);

			_elems.resize(n, val);
		}
		explicit Vector(int n, Type* vals) 
		{
			if (n < 0)
				throw VectorInitializationError("Vector::Vector - negative size", n);

			_elems.resize(n);
			for (int i = 0; i < n; ++i)
				_elems[i] = vals[i];
		}
		explicit Vector(std::vector<Type> values) : _elems(values) {}
		explicit Vector(std::initializer_list<Type> list) : _elems(list) {}

		// not really needed, but let's be explicit
		Vector(const Vector& b) = default; 
		Vector(Vector&& b) = default;
		Vector& operator=(const Vector& b) = default; 
		Vector& operator=(Vector&& b) = default;

		////////////////            Standard std::vector stuff             ////////////////////
		int  size() const { return (int)_elems.size(); }
		bool empty() const { return _elems.empty(); }

		void clear() { _elems.clear(); }
		void resize(int newLen)		{ _elems.resize(newLen); }

		////////////////////////            Standard stuff             ////////////////////////
		static Vector GetUnitVector(int dimVec, int indUnit)
		{
			if (indUnit < 0 || indUnit >= dimVec)
				throw VectorDimensionError("Vector::GetUnitVector - wrong unit index", dimVec, indUnit);

			Vector ret(dimVec);
			ret[indUnit] = Type{ 1.0 };
			return ret;
		}
		
		static bool AreEqual(const Vector& a, const Vector& b, Type eps = Defaults::VectorEqualityPrecision)
		{
			return a.IsEqualTo(b, eps);
		}
		bool IsEqualTo(const Vector& b, Real eps = Defaults::VectorEqualityPrecision) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::IsEqual - vectors must be equal size", size(), b.size());

			for (int i = 0; i < size(); i++)
			{
				if (Abs((*this)[i] - b[i]) > eps)
					return false;
			}
			return true;
		}

		///////////////////////////            Operators             ///////////////////////////
		Type&       operator[](int n)       { return _elems[n]; }
		const Type& operator[](int n) const { return _elems[n]; }

		// checked access
		Type& at(int n)	{
			if(n < 0 || n >= size())
				throw VectorDimensionError("Vector::at - index out of bounds", size(), n);
			else
				return _elems[n];
		}
		Type  at(int n) const { 
			if(n < 0 || n >= size())
				throw VectorDimensionError("Vector::at - index out of bounds", size(), n);
			else
				return _elems[n];
		}

		Vector operator-() const         // unary minus
		{
			Vector ret(size());
			for (int i = 0; i < size(); i++)
				ret._elems[i] = Type{ -1 } * (*this)[i];
			return ret;
		}
		Vector operator+(const Vector& b) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::operator+() - vectors must be equal size", size(), b.size());

			Vector ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret._elems[i] = (*this)[i] + b._elems[i];
			return ret;
		}
		Vector operator-(const Vector& b) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::operator-() - vectors must be equal size", size(), b.size());

			Vector ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret._elems[i] = (*this)[i] - b._elems[i];
			return ret;
		}
		bool   operator==(const Vector& b) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::operator==() - vectors must be equal size", size(), b.size());

			for (int i = 0; i < size(); i++)
				if ((*this)[i] != b[i])
					return false;

			return true;
		}

		Vector operator*(Type b)
		{
			Vector ret(size());;
			for (int i = 0; i < size(); i++)
				ret._elems[i] = b * _elems[i];
			return ret;
		}
		Vector operator/(Type b)
		{
			Vector ret(size());
			for (int i = 0; i < size(); i++)
				ret._elems[i] = _elems[i] / b;
			return ret;
		}
		
		friend Vector operator*(Type a, const Vector& b)
		{
			Vector ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret._elems[i] = a * b._elems[i];
			return ret;
		}
		
		//////////////////////                 Operations                 ///////////////////////
		Type ScalarProductCartesian(const Vector& b) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::ScalarProductCartesian - vectors must be equal size", size(), b.size());

			Type product{ 0.0 };
			for (int i = 0; i < size(); i++)
				product += (*this)[i] * b[i];
			return product;
		}
		Type NormL2() const
		{
			Type norm{ 0.0 };
			for (int i = 0; i < size(); i++)
				norm += (*this)[i] * (*this)[i];
			return std::sqrt(norm);
		}
		Type AngleToVector(const Vector& b) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::AngleToVector - vectors must be equal size", size(), b.size());

			Type cosAngle = ScalarProductCartesian(b) / (NormL2() * b.NormL2());
			return std::acos(cosAngle);
		}

		///////////////////////////               I/O                 ///////////////////////////
		std::string to_string(int width, int precision) const
		{
			std::stringstream str;

			Print(str, width, precision);

			return str.str();
		}
		void Print(std::ostream& stream, int width, int precision) const
		{
			stream << "[";
			bool first = true;
			for (const Type& x : _elems)
			{
				if (!first)
					stream << ", ";
				else
					first = false;

				stream << std::setw(width) << std::setprecision(precision) << x;
			}
			stream << "]";
		}
    std::ostream& PrintLine(std::ostream& stream, std::string msg, int width, int precision) const
		{
			stream << msg;
			Print(stream, width, precision);
			stream << std::endl;

			return stream;
		}
		std::ostream& Print(std::ostream& stream, int width, int precision, Real zeroThreshold) const
		{
			stream << "[";
			bool first = true;
			for (const Type& x : _elems)
			{
				if (!first)
					stream << ", ";
				else
					first = false;

        if( Abs(x) > zeroThreshold )
				  stream << std::setw(width) << std::setprecision(precision) << x;
        else
          stream << std::setw(width) << std::setprecision(precision) << 0.0;
			}
			stream << "]";

			return stream;
		}      
		friend std::ostream& operator<<(std::ostream& stream, const Vector& a)
		{
			a.Print(stream, Defaults::VectorPrintWidth, Defaults::VectorPrintPrecision);

			return stream;
		}
	};

	typedef Vector<int>     VectorInt;
	typedef Vector<float>   VectorFlt;
	typedef Vector<double>  VectorDbl;
	typedef Vector<Complex> VectorComplex;

	typedef Vector<int>     VecI;
	typedef Vector<float>   VecF;
	typedef Vector<double>  VecD;
	typedef Vector<Complex> VecC;
}

///////////////////////////   ./include/base/VectorN.h   ///////////////////////////

namespace MML
{
	template<class Type, int N>
	class VectorN
	{
	protected:
		Type  _val[N] = { 0 };

	public:
		///////////////////////          Constructors and destructor       //////////////////////
		VectorN() {}
		explicit VectorN(const Type& init_val) {
			for (int i = 0; i < N; ++i)
				_val[i] = init_val;
		}
		explicit VectorN(std::initializer_list<Type> list)
		{
			int count = 0;
			for (auto element : list)
			{
				_val[count] = element;
				++count;

				if (count >= N)
					break;
			}
		}
		explicit VectorN(std::vector<Type> list)
		{
			int count{ 0 };
			for (auto element : list)
			{
				_val[count] = element;
				++count;

				if (count >= N)
					break;
			}
		}
		explicit VectorN(Type* vals)
		{
			for (int i = 0; i < N; ++i)
				_val[i] = vals[i];
		}
		static VectorN GetUnitVector(int indUnit)
		{
			VectorN ret;
			ret[indUnit] = 1.0;
			return ret;
		}

		////////////////////////            Standard stuff             ////////////////////////
		int size() const { return N; }
		void clear() {
			for (int i = 0; i < N; ++i)
				_val[i] = Type{ 0 };
		}

		VectorN GetAsUnitVector() const
		{
			return VectorN{ (*this) / NormL2() };
		}
		void MakeUnitVector()
		{
			(*this) = GetAsUnitVector();
		}

		bool IsEqualTo(const VectorN& b, Type eps = Defaults::VectorEqualityPrecision) const
		{
			for (int i = 0; i < N; i++)
			{
				if (Abs((*this)[i] - b[i]) > eps)
					return false;
			}
			return true;
		}
		static bool AreEqual(const VectorN& a, const VectorN& b, Type eps = Defaults::VectorEqualityPrecision)
		{
			return a.IsEqualTo(b, eps);
		}
		bool IsNullVec() const
		{
			for (int i = 0; i < N; i++)
				if (_val[i] != 0.0)
					return false;

			return true;
		}
		///////////////////////////            Operators             ///////////////////////////
		Type& operator[](int n) { return _val[n]; }
		Type  operator[](int n) const { return _val[n]; }

		// checked access
		Type& at(int n) {
			if (n < 0 || n >= N)
				throw VectorDimensionError("VectorN::at - index out of bounds", N, n);
			else
				return _val[n];
		}
		Type  at(int n) const {
			if (n < 0 || n >= N)
				throw VectorDimensionError("VectorN::at - index out of bounds", N, n);
			else
				return _val[n];
		}

		VectorN operator-() const        // unary minus
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = Type{ -1 } *_val[i];
			return ret;
		}
		VectorN operator+(const VectorN& b) const
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = _val[i] + b._val[i];
			return ret;
		}
		VectorN operator-(const VectorN& b) const
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = _val[i] - b._val[i];
			return ret;
		}
		bool operator==(const VectorN& b) const
		{
			for (int i = 0; i < size(); i++)
				if ((*this)[i] != b[i])
					return false;

			return true;
		}

		VectorN operator*(const Type &b) const
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = _val[i] * b;
			return ret;
		}
		VectorN operator/(const Type &b) const
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = _val[i] / b;
			return ret;
		}		
		friend VectorN operator*(Type a, const VectorN<Type, N>& b)
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = a * b[i];
			return ret;
		}

		//////////////////////                 Operations                 ///////////////////////
		Type ScalarProductCartesian(const VectorN& b) const
		{
			Type product{ 0.0 };
			for (int i = 0; i < N; i++)
				product += (*this)[i] * b[i];
			return product;
		}
		Type NormL2() const
		{
			Type norm{ 0.0 };
			for (int i = 0; i < N; i++)
				norm += (*this)[i] * (*this)[i];
			return std::sqrt(norm);
		}
		Type AngleToVector(const VectorN& b) const
		{
			if (size() != b.size())
				throw VectorDimensionError("VectorN::AngleToVector - vectors must be equal size", size(), b.size());

			Type cosAngle = ScalarProductCartesian(b) / (NormL2() * b.NormL2());
			return std::acos(cosAngle);
		}

		///////////////////////////               I/O                 ///////////////////////////
		std::string to_string(int width, int precision) const
		{
			std::stringstream str;

			Print(str, width, precision);

			return str.str();
		}
		std::ostream& Print(std::ostream& stream, int width, int precision) const
		{
			stream << "[";
			bool first = true;
			for (const Type& x : _val)
			{
				if (!first)
					stream << ", ";
				else
					first = false;

				stream << std::setw(width) << std::setprecision(precision) << x;
			}
			stream << "]";

			return stream;
		}
		std::ostream& PrintLine(std::ostream& stream, std::string msg, int width, int precision) const
		{
			stream << msg;
			Print(stream, width, precision);
			stream << std::endl;

			return stream;
		}

		std::ostream& Print(std::ostream& stream, int width, int precision, Real zeroThreshold) const
		{
			stream << "[";
			bool first = true;
			for (const Type& x : _val)
			{
				if (!first)
					stream << ", ";
				else
					first = false;

				if (Abs(x) > zeroThreshold)
					stream << std::setw(width) << std::setprecision(precision) << x;
				else
					stream << std::setw(width) << std::setprecision(precision) << 0.0;
			}
			stream << "]";

			return stream;
		}
		friend std::ostream& operator<<(std::ostream& stream, const VectorN<Type, N>& a)
		{
			a.Print(stream, 15, 10);

			return stream;
		}
	};


	typedef VectorN<float, 2> Vec2Flt;
	typedef VectorN<float, 3> Vec3Flt;
	typedef VectorN<float, 4> Vec4Flt;

	typedef VectorN<double, 2> Vec2Dbl;
	typedef VectorN<double, 3> Vec3Dbl;
	typedef VectorN<double, 4> Vec4Dbl;

	typedef VectorN<Complex, 2> Vec2Complex;
	typedef VectorN<Complex, 3> Vec3Complex;
	typedef VectorN<Complex, 4> Vec4Complex;

	typedef VectorN<float, 2> Vec2F;
	typedef VectorN<float, 3> Vec3F;
	typedef VectorN<float, 4> Vec4F;

	typedef VectorN<double, 2> Vec2D;
	typedef VectorN<double, 3> Vec3D;
	typedef VectorN<double, 4> Vec4D;

	typedef VectorN<Complex, 2> Vec2C;
	typedef VectorN<Complex, 3> Vec3C;
	typedef VectorN<Complex, 4> Vec4C;
}

///////////////////////////   ./include/base/Matrix.h   ///////////////////////////


namespace MML
{
	template<class Type>
	class Matrix
	{
	private:
		int  _rows;
		int  _cols;
		Type** _data;

		void Init(int rows, int cols)
		{
			if (rows <= 0 || cols < 0)
				throw MatrixDimensionError("Matrix::Init - rowNum and colNum must be positive", rows, cols, -1, -1);

			_rows = rows;
			_cols = cols;
			int numElem = rows * cols;

			_data = new Type * [rows];
			if (_data)
			{
				_data[0] = numElem > 0 ? new Type[numElem] : nullptr;
				
				if( numElem > 0 && _data[0] == nullptr)
					throw MatrixAllocationError("Matrix::Init - allocation error", rows, cols);

				for (int i = 1; i < rows; i++)
					_data[i] = _data[i - 1] + cols;
			}
			else
				throw MatrixAllocationError("Matrix::Init - allocation error", rows, cols);
		}

	public:
		typedef Type value_type;      // make Type available externally

		///////////////////////          Constructors and destructor       //////////////////////
		explicit Matrix() : _rows(0), _cols(0), _data{ nullptr } {}
		explicit Matrix(int rows, int cols) : _rows(rows), _cols(cols)
		{
			Init(rows, cols);
			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					_data[i][j] = 0;
		}
		explicit Matrix(int rows, int cols, const Type& val) : _rows(rows), _cols(cols)
		{
			Init(rows, cols);
			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					_data[i][j] = val;
		}
		// useful if you have a pointer to continuous 2D array (can be in row-, or column-wise memory layout)
		explicit Matrix(int rows, int cols, Type* val, bool isRowWise = true) : _rows(rows), _cols(cols)
		{
			Init(rows, cols);

			if( isRowWise)
				for (int i = 0; i < _rows; ++i)
					for (int j = 0; j < _cols; ++j)
						_data[i][j] = *val++;
			else
				for (int j = 0; j < _cols; ++j)
					for (int i = 0; i < _rows; ++i)
						_data[i][j] = *val++;
		}
		// in strict mode, you must supply ALL necessary values for complete matrix initialization
		explicit Matrix(int rows, int cols, std::initializer_list<Type> values, bool strictMode = true) : _rows(rows), _cols(cols)
		{
			Init(rows, cols);

			auto val = values.begin();
			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					if (val != values.end())
						_data[i][j] = *val++;
					else {
						if (strictMode)
							throw MatrixDimensionError("Matrix::Matrix - not enough values in initializer list", _rows, _cols, -1, -1);
						else
							_data[i][j] = Type{0};
					}
		}
		
		Matrix(const Matrix& m) : _rows(m._rows), _cols(m._cols)
		{
			Init(m._rows, m._cols);

			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					_data[i][j] = m._data[i][j];
		}
		// creating submatrix from given matrix 'm'
		Matrix(const Matrix& m, int ind_row, int ind_col, int row_num, int col_num)
		{
			if (ind_row < 0 || ind_row >= m._rows || ind_col < 0 || ind_col >= m._cols)
				throw MatrixDimensionError("Matrix::Matrix - invalid row or column index", m._rows, m._cols, ind_row, ind_col);

			if (row_num <= 0 || col_num <= 0)
				throw MatrixDimensionError("Matrix::Matrix - rowNum and colNum must be positive", row_num, col_num, -1, -1);

			if (ind_row + row_num > m._rows || ind_col + col_num > m._cols)
				throw MatrixDimensionError("Matrix::Matrix - submatrix out of bounds", m._rows, m._cols, ind_row, ind_col);

			_rows = row_num;
			_cols = col_num;

			Init(row_num, col_num);

			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					_data[i][j] = m._data[ind_row + i][ind_col + j];
		}
		Matrix(Matrix&& m)
		{
			_data = m._data;

			_rows = m._rows;
			_cols = m._cols;

			m._rows = 0;
			m._cols = 0;
			m._data = nullptr;
		}
		~Matrix()
		{
			if (_data != NULL) {
				delete[](_data[0]);
				delete[](_data);
			}
		}

		void Resize(int rows, int cols)
		{
			if (rows <= 0 || cols <= 0)
				throw MatrixDimensionError("Matrix::Resize - rowNum and colNum must be positive", rows, cols, -1, -1);

			if (rows == RowNum() && cols == ColNum())      // nice :)
				return;

			if (_data != NULL) {
				delete[](_data[0]);
				delete[](_data);
			}

			Init(rows, cols);

			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					_data[i][j] = 0;
		}

		void Expand(int newRows, int newCols, bool retainValues = true)
		{
		}
		void AddRow(bool retainValues = true)
		{
		}
		void AddCol(bool retainValues = true)
		{
		}

		///////////////////////              Standard stuff                //////////////////////
		int RowNum() const { return _rows; }
		int ColNum() const { return _cols; }

		static Matrix GetUnitMatrix(int dim)
		{
			if (dim <= 0)
				throw MatrixDimensionError("Matrix::GetUnitMatrix - dimension must be positive", dim, dim, -1, -1);

			Matrix unitMat(dim, dim);
			unitMat.MakeUnitMatrix();

			return unitMat;
		}
		void   MakeUnitMatrix(void)
		{
			if (_rows == _cols)
			{
				for (int i = 0; i < _rows; i++)
					for (int j = 0; j < _cols; j++)
						if (i == j)
							_data[i][j] = 1;
						else
							_data[i][j] = 0;
			}
			else
				throw MatrixDimensionError("Matrix::MakeUnitMatrix - must be square matrix", _rows, _cols, -1, -1);
		}

		Matrix GetLower(bool includeDiagonal = true) const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("Matrix::GetLower - must be square matrix", _rows, _cols, -1, -1);

			Matrix ret(RowNum(), ColNum());
			for (int i = 0; i < RowNum(); i++)
			{
				if (includeDiagonal)
					for (int j = 0; j <= i; j++)
						ret[i][j] = _data[i][j];
				else
					for (int j = 0; j < i; j++)
						ret[i][j] = _data[i][j];
			}

			return ret;
		}
		Matrix GetUpper(bool includeDiagonal = true) const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("Matrix::GetUpper - must be square matrix", _rows, _cols, -1, -1);

			Matrix ret(RowNum(), ColNum());
			for (int i = 0; i < RowNum(); i++)
			{
				if (includeDiagonal)
					for (int j = i; j < ColNum(); j++)
						ret[i][j] = _data[i][j];
				else
					for (int j = i + 1; j < ColNum(); j++)
						ret[i][j] = _data[i][j];
			}

			return ret;
		}

		///////////////////////          Matrix to Vector conversions      //////////////////////
		Vector<Type> VectorFromRow(int rowInd) const
		{
			if (rowInd < 0 || rowInd >= RowNum())
				throw MatrixAccessBoundsError("VectorFromRow - invalid row index", rowInd, 0, RowNum(), ColNum());

			Vector<Type> ret(ColNum());
			for (int i = 0; i < ColNum(); i++)
				ret[i] = (*this)(rowInd, i);

			return ret;
		}
		Vector<Type> VectorFromColumn(int colInd) const
		{
			if (colInd < 0 || colInd >= ColNum())
				throw MatrixAccessBoundsError("VectorFromColumn - invalid column index", 0, colInd, RowNum(), ColNum());

			Vector<Type> ret(RowNum());
			for (int i = 0; i < RowNum(); i++)
				ret[i] = (*this)(i, colInd);

			return ret;
		}
		Vector<Type> VectorFromDiagonal() const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("VectorFromDiagonal - must be square matrix", RowNum(), ColNum(), -1, -1);

			Vector<Type> ret(RowNum());
			for (int i = 0; i < RowNum(); i++)
				ret[i] = (*this)(i, i);

			return ret;
		}

		///////////////////////               Matrix properties            //////////////////////
		bool IsUnit(double eps = Defaults::IsMatrixUnitPrecision) const
		{
			for (int i = 0; i < RowNum(); i++)
				for (int j = 0; j < ColNum(); j++)
					if (i != j && Abs((*this)[i][j]) > eps)
						return false;

			for (int i = 0; i < RowNum(); i++)
				if (Abs((*this)[i][i] - Real{ 1.0 }) > eps)
					return false;

			return true;
		}
		bool IsDiagonal(double eps = Defaults::IsMatrixDiagonalPrecision) const
		{
			for (int i = 0; i < RowNum(); i++)
				for (int j = 0; j < ColNum(); j++)
					if (i != j && Abs((*this)[i][j]) > eps)
						return false;

			return true;
		}
		bool IsDiagDominant() const
		{
			for (int i = 0; i < RowNum(); i++)
			{
				Type sum = 0.0;
				for (int j = 0; j < ColNum(); j++)
					if (i != j)
						sum += Abs((*this)[i][j]);

				if (Abs((*this)[i][i]) < sum)
					return false;
			}
			return true;
		}
		bool IsSymmetric() const
		{
			if (RowNum() != ColNum())
				return false;

			for (int i = 0; i < RowNum(); i++)
				for (int j = i + 1; j < ColNum(); j++)
					if ((*this)[i][j] != (*this)[j][i])
						return false;

			return true;
		}
		bool IsAntiSymmetric() const
		{
			if (RowNum() != ColNum())
				return false;
			for (int i = 0; i < RowNum(); i++)
				for (int j = i + 1; j < ColNum(); j++)
					if ((*this)[i][j] != -(*this)[j][i])
						return false;
			return true;
		}

		///////////////////////             Assignment operators           //////////////////////
		Matrix& operator=(const Matrix& m)
		{
			if (this == &m)
				return *this;

			if (_rows != m._rows || _cols != m._cols)
			{
				delete[](_data[0]);
				delete[](_data);

				Init(m._rows, m._cols);
			}

			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					_data[i][j] = m._data[i][j];

			return *this;
		}
		Matrix& operator=(Matrix&& m) noexcept
		{
			if (this == &m)
				return *this;

			std::swap(_data, m._data);
			std::swap(_rows, m._rows);
			std::swap(_cols, m._cols);

			return *this;
		}

		///////////////////////               Access operators             //////////////////////
		Type*				operator[](int i)							{ return _data[i]; }
		const Type* operator[](const int i) const { return _data[i]; }

		Type  operator()(int i, int j) const { return _data[i][j]; }
		Type& operator()(int i, int j)			 { return _data[i][j]; }

		// version with checked access
		Type  at(int i, int j) const
		{
			if (i < 0 || i >= RowNum() || j < 0 || j >= ColNum())
				throw MatrixAccessBoundsError("Matrix::at", i, j, RowNum(), ColNum());

			return _data[i][j];
		}
		Type& at(int i, int j)
		{
			if (i < 0 || i >= RowNum() || j < 0 || j >= ColNum())
				throw MatrixAccessBoundsError("Matrix::at", i, j, RowNum(), ColNum());

			return _data[i][j];
		}

		///////////////////////              Equality operations           //////////////////////
		bool operator==(const Matrix& b) const
		{
			if (_rows != b._rows || _cols != b._cols)
				return false;

			for (int i = 0; i < _rows; i++)
				for (int j = 0; j < _cols; j++)
					if (_data[i][j] != b._data[i][j])
						return false;

			return true;
		}
		bool operator!=(const Matrix& b) const
		{
			return !(*this == b);
		}

		bool IsEqualTo(const Matrix<Type>& b, Type eps = Defaults::MatrixEqualityPrecision) const
		{
			if (RowNum() != b.RowNum() || ColNum() != b.ColNum())
				return false;

			for (int i = 0; i < RowNum(); i++)
				for (int j = 0; j < ColNum(); j++)
				{
					if (Abs(_data[i][j] - b._data[i][j]) > eps)
						return false;
				}

			return true;
		}
		static bool AreEqual(const Matrix& a, const Matrix& b, Type eps = Defaults::MatrixEqualityPrecision)
		{
			return a.IsEqualTo(b, eps);
		}

		///////////////////////              Arithmetic operators          //////////////////////
		Matrix operator-() const            // unary minus
		{
			Matrix temp(_rows, _cols);
			for (int i = 0; i < _rows; i++)
				for (int j = 0; j < _cols; j++)
					temp._data[i][j] = -_data[i][j];

			return temp;
		}
		Matrix operator+(const Matrix& b) const
		{
			if (_rows != b._rows || _cols != b._cols)
				throw MatrixDimensionError("Matrix::operator+() - must be same dim", _rows, _cols, b._rows, b._cols);

			Matrix temp(_rows, _cols);
			for (int i = 0; i < _rows; i++)
				for (int j = 0; j < _cols; j++)
					temp._data[i][j] = b._data[i][j] + _data[i][j];

			return temp;
		}
		Matrix operator-(const Matrix& b) const
		{
			if (_rows != b._rows || _cols != b._cols)
				throw MatrixDimensionError("Matrix::operator-() - must be same dim", _rows, _cols, b._rows, b._cols);

			Matrix temp(_rows, _cols);
			for (int i = 0; i < _rows; i++)
				for (int j = 0; j < _cols; j++)
					temp._data[i][j] = _data[i][j] - b._data[i][j];

			return temp;
		}
		Matrix operator*(const Matrix& b) const
		{
			if (ColNum() != b.RowNum())
				throw MatrixDimensionError("Matrix::operator*() - a.colNum must be equal to b.rowNum", _rows, _cols, b._rows, b._cols);

			Matrix	ret(RowNum(), b.ColNum());
			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++)
				{
					ret._data[i][j] = 0;
					for (int k = 0; k < ColNum(); k++)
						ret._data[i][j] += _data[i][k] * b._data[k][j];
				}

			return	ret;
		}

		Matrix operator*(const Type &b) const
		{
			int	i, j;
			Matrix	ret(RowNum(), ColNum());

			for (i = 0; i < RowNum(); i++)
				for (j = 0; j < ColNum(); j++)
					ret[i][j] = _data[i][j] * b;

			return ret;
		}
		Matrix operator/(const Type &b) const
		{
			int	i, j;
			Matrix	ret(RowNum(), ColNum());

			for (i = 0; i < RowNum(); i++)
				for (j = 0; j < ColNum(); j++)
					ret[i][j] = _data[i][j] / b;

			return ret;
		}
		Vector<Type> operator*(const Vector<Type>& b) const
		{
			if (ColNum() != b.size())
				throw MatrixDimensionError("operator*(Mat a, Vec b) - a.colNum must be equal to vector size", _rows, _cols, (int)b.size(), -1);

			Vector<Type>	ret(RowNum());
			for (int i = 0; i < RowNum(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < ColNum(); j++)
					ret[i] += _data[i][j] * b[j];
			}

			return ret;
		}

		friend Matrix operator*(const Type &a, const Matrix<Type>& b)
		{
			int	i, j;
			Matrix	ret(b.RowNum(), b.ColNum());

			for (i = 0; i < b.RowNum(); i++)
				for (j = 0; j < b.ColNum(); j++)
					ret[i][j] = a * b._data[i][j];

			return ret;
		}
		friend Vector<Type> operator*(const Vector<Type>& a, const Matrix<Type>& b)
		{
			if (a.size() != b.RowNum())
			{
				//std::string error = std::format("Hello {}!\n", "world");
				throw MatrixDimensionError("operator*(Vec a, Mat b) - vector size must be equal to b.rowNum", (int)a.size(), -1, b._rows, b._cols);
			}

			Vector<Type>	ret(b.ColNum());
			for (int i = 0; i < b.ColNum(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < b.RowNum(); j++)
					ret[i] += a[j] * b(j, i);
			}

			return ret;
		}

		///////////////////////            Trace, Inverse & Transpose      //////////////////////
		Type   Trace() const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("Matrix::Trace - must be square matrix", _rows, _cols, -1, -1);

			Type sum = 0;
			for (int i = 0; i < RowNum(); i++)
				sum += _data[i][i];

			return sum;
		}

		void   Invert()
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("Matrix::Invert - must be square matrix", _rows, _cols, -1, -1);

			Matrix& a = *this;
			Matrix  b(RowNum(), 1);      // dummy rhs

			b(0, 0) = 1.0;

			int i, icol, irow, j, k, l, ll;
			Type dum, pivinv;
			Real big;

			int n = a.RowNum();
			int m = b.ColNum();
			std::vector<int> indxc(n), indxr(n), ipiv(n);
			for (j = 0; j < n; j++) ipiv[j] = 0;
			for (i = 0; i < n; i++) {
				big = 0.0;
				for (j = 0; j < n; j++)
					if (ipiv[j] != 1)
						for (k = 0; k < n; k++) {
							if (ipiv[k] == 0) {
								if (Abs(a[j][k]) >= big) {
									big = Abs(a[j][k]);
									irow = j;
									icol = k;
								}
							}
						}
				++(ipiv[icol]);
				if (irow != icol) {
					for (l = 0; l < n; l++) std::swap(a[irow][l], a[icol][l]);
					for (l = 0; l < m; l++) std::swap(b[irow][l], b[icol][l]);
				}
				indxr[i] = irow;
				indxc[i] = icol;

				if (a[icol][icol] == 0.0)
					throw SingularMatrixError("Matrix::Invert, Singular Matrix");

				pivinv = 1.0 / a[icol][icol];
				a[icol][icol] = 1.0;
				for (l = 0; l < n; l++) a[icol][l] *= pivinv;
				for (l = 0; l < m; l++) b[icol][l] *= pivinv;
				for (ll = 0; ll < n; ll++)
					if (ll != icol) {
						dum = a[ll][icol];
						a[ll][icol] = 0.0;
						for (l = 0; l < n; l++) a[ll][l] -= a[icol][l] * dum;
						for (l = 0; l < m; l++) b[ll][l] -= b[icol][l] * dum;
					}
			}
			for (l = n - 1; l >= 0; l--) {
				if (indxr[l] != indxc[l])
					for (k = 0; k < n; k++)
						std::swap(a[k][indxr[l]], a[k][indxc[l]]);
			}
		}
		Matrix GetInverse() const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("Matrix::GetInverse - must be square matrix", _rows, _cols, -1, -1);

			Matrix a(*this);              // making a copy, where inverse will be stored at the end
			a.Invert();

			return a;
		}

		void   Transpose()
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("Matrix::Transpose - in-place Transpose possible only for square matrix", _rows, _cols, -1, -1);

			for (int i = 0; i < RowNum(); i++)
				for (int j = i + 1; j < ColNum(); j++)
					std::swap(_data[i][j], _data[j][i]);
		}
		Matrix GetTranspose() const
		{
			Matrix ret(ColNum(), RowNum());

			for (int i = 0; i < ColNum(); i++)
				for (int j = 0; j < RowNum(); j++)
					ret[i][j] = _data[j][i];

			return ret;
		}

		///////////////////////                    I/O                    //////////////////////
		std::string to_string(int width, int precision) const
		{
			std::stringstream str;

			Print(str, width, precision);

			return str.str();
		}
		void   Print(std::ostream& stream, int width, int precision) const
		{
			stream << "Rows: " << RowNum() << " Cols: " << ColNum();
      
      if( RowNum() == 0 || ColNum() == 0) {
        stream << " - Empty matrix" << std::endl;
        return;
      }
      else
        stream << std::endl;

			for (int i = 0; i < RowNum(); i++)
			{
				stream << "[ ";
				for (int j = 0; j < ColNum(); j++)
				{
					if( j == ColNum() - 1 )
						stream << std::setw(width) << std::setprecision(precision) << _data[i][j];
					else
						stream << std::setw(width) << std::setprecision(precision) << _data[i][j] << ", ";
				}
				if( i == RowNum() - 1 )
					stream << " ]";
				else
					stream << " ]" << std::endl;
			}
		}
		void   Print(std::ostream& stream, int width, int precision, Real zeroThreshold) const
		{
			stream << "Rows: " << RowNum() << " Cols: " << ColNum(); 
      
      if( RowNum() == 0 || ColNum() == 0) {
        stream << " - Empty matrix" << std::endl;
        return;
      }
      else
        stream << std::endl;

			for (int i = 0; i < RowNum(); i++)
			{
				stream << "[ ";
				for (int j = 0; j < ColNum(); j++)
				{
					Type value{0};
					if (Abs(_data[i][j]) > zeroThreshold)
						value = _data[i][j];

					if( j == ColNum() - 1 )
						stream << std::setw(width) << std::setprecision(precision) << value;
					else
						stream << std::setw(width) << std::setprecision(precision) << value << ", ";
				}
				if( i == RowNum() - 1 )
					stream << " ]";
				else
					stream << " ]" << std::endl;
			}
		}
		friend std::ostream& operator<<(std::ostream& stream, const Matrix& a)
		{
			a.Print(stream, 10, 3);

			return stream;
		}

		static bool LoadFromFile(std::string inFileName, Matrix& outMat)
		{
			std::ifstream file(inFileName);

			if (file.is_open())
			{
				int rows, cols;
				file >> rows >> cols;

				outMat.Resize(rows, cols);
				for (int  i = 0; i < outMat.RowNum(); i++)
					for (int j = 0; j < outMat.ColNum(); j++)
						file >> outMat[i][j];

				file.close();
			}
			else {
				std::cerr << "Error: could not open file " << inFileName << " for reading." << std::endl;
				return false;
			}

			return true;
		}
		static bool SaveToFile(const Matrix& mat, std::string inFileName)
		{
			std::ofstream file(inFileName);

			if (file.is_open())
			{
				file << mat.RowNum() << " " << mat.ColNum() << std::endl;
				for (int i = 0; i < mat.RowNum(); i++)
				{
					for (int j = 0; j < mat.ColNum(); j++)
						file << mat(i, j) << " ";
					file << std::endl;
				}
				file.close();
			}
			else {
				std::cerr << "Error: could not create file " << inFileName << " for writing." << std::endl;
				return false;
			}

			return true;
		}
	};

	//////////////////////               Default Matrix typdefs         ////////////////////
	typedef Matrix<int>     MatrixInt;
	typedef Matrix<float>   MatrixFlt;
	typedef Matrix<double>  MatrixDbl;
	typedef Matrix<Complex> MatrixComplex;

	typedef Matrix<int>     MatI;
	typedef Matrix<float>   MatF;
	typedef Matrix<double>  MatD;
	typedef Matrix<Complex> MatC;
}
///////////////////////////   ./include/base/MatrixNM.h   ///////////////////////////


namespace MML
{
	template <class Type, int N, int M>
	class MatrixNM
	{
	public:
		Type _vals[N][M] = { {0} };

	public:
		//////////////////////////             Constructors           /////////////////////////
		MatrixNM() {}
		MatrixNM(std::initializer_list<Type> values)
		{
			auto val = values.begin();
			for (size_t i = 0; i < RowNum(); ++i)
				for (size_t j = 0; j < ColNum(); ++j)
					if (val != values.end())
						_vals[i][j] = *val++;
					else
						_vals[i][j] = 0.0;
		}
		MatrixNM(const MatrixNM& m)
		{
			for (size_t i = 0; i < RowNum(); ++i)
				for (size_t j = 0; j < ColNum(); ++j)
					_vals[i][j] = m._vals[i][j];
		}
		MatrixNM(const Type& m)        // initialize as diagonal matrix
		{
			for (int i = 0; i < N; i++)
				_vals[i][i] = Type{ m };
		}

		typedef Type value_type;      // make T available externally

		////////////////////////            Standard stuff             ////////////////////////
		int RowNum() const { return N; }
		int ColNum() const { return M; }

		static MatrixNM GetUnitMatrix()
		{
			MatrixNM unitMat;

			for (int i = 0; i < N; i++)
				unitMat._vals[i][i] = 1.0;

			return unitMat;
		}
		void   MakeUnitMatrix(void)
		{
			if (RowNum() == ColNum())
			{
				for (int i = 0; i < RowNum(); i++)
					for (int j = 0; j < ColNum(); j++)
						if (i == j)
							_vals[i][j] = 1;
						else
							_vals[i][j] = 0;
			}
			else
				throw MatrixDimensionError("MatrixNM::MakeUnitMatrix - must be square matrix", N, M, -1, -1);
		}

		MatrixNM GetLower(bool includeDiagonal = true) const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("Matrix::GetLower - must be square matrix", N, M, -1, -1);

			MatrixNM ret;
			for (int i = 0; i < RowNum(); i++)
			{
				if (includeDiagonal)
					for (int j = 0; j < i + 1; j++)
						ret[i][j] = _vals[i][j];
				else
					for (int j = 0; j < i; j++)
						ret[i][j] = _vals[i][j];
			}

			return ret;
		}
		MatrixNM GetUpper(bool includeDiagonal = true) const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("Matrix::GetUpper - must be square matrix", N, M, -1, -1);

			MatrixNM ret;
			for (int i = 0; i < RowNum(); i++)
			{
				if (includeDiagonal)
					for (int j = i; j < ColNum(); j++)
						ret[i][j] = _vals[i][j];
				else
					for (int j = i + 1; j < ColNum(); j++)
						ret[i][j] = _vals[i][j];
			}

			return ret;
		}

		/////////////////////          Vector-Matrix conversion           /////////////////////
		static MatrixNM<Type, 1, M> RowMatrixFromVector(const VectorN<Type, M>& b)
		{
			MatrixNM<Type, 1, M>  ret;
			for (int j = 0; j < M; j++)
				ret._vals[0][j] = b[j];

			return ret;
		}
		static MatrixNM<Type, N, 1> ColumnMatrixFromVector(const VectorN<Type, N>& b)
		{
			MatrixNM<Type, N, 1>  ret;
			for (int i = 0; i < M; i++)
				ret._vals[i][0] = b[i];

			return ret;
		}
		static MatrixNM<Type, N, N> DiagonalMatrixFromVector(const VectorN<Type, N>& b)
		{
			MatrixNM<Type, N, N> ret;
			for (int i = 0; i < N; i++)
				ret[i][i] = b[i];

			return ret;
		}
		static VectorN<Type, M> VectorFromRow(const MatrixNM<Type, N, M>& a, int rowInd)
		{
			VectorN<Type, M> ret;
			for (int j = 0; j < M; j++)
				ret[j] = a._vals[rowInd][j];

			return ret;
		}
		static VectorN<Type, N> VectorFromColumn(const MatrixNM<Type, N, M>& a, int colInd)
		{
			VectorN<Type, N> ret;
			for (int i = 0; i < N; i++)
				ret[i] = a._vals[i][colInd];

			return ret;
		}
		static VectorN<Type, N> VectorFromDiagonal(const MatrixNM<Type, N, N>& a)
		{
			VectorN<Type, N> ret;
			for (int i = 0; i < N; i++)
				ret[i] = a._vals[i][i];

			return ret;
		}

		/////////////////////            Assignment operators             ////////////////////
		MatrixNM& operator=(const MatrixNM& m)
		{
			if (this == &m)
				return *this;

			for (size_t i = 0; i < RowNum(); ++i)
				for (size_t j = 0; j < ColNum(); ++j)
					_vals[i][j] = m._vals[i][j];

			return *this;
		}
		MatrixNM& operator=(const Type& m)
		{
			if (this == &m)
				return *this;

			for (size_t i = 0; i < RowNum(); ++i)
				for (size_t j = 0; j < ColNum(); ++j)
					_vals[i][j] = m;

			return *this;
		}

		////////////////////            Access operators             ///////////////////////
		Type* operator[](int i) { return _vals[i]; }
		const Type* operator[](const int i) const { return _vals[i]; }

		Type  operator()(int i, int j) const { return _vals[i][j]; }
		Type& operator()(int i, int j) { return _vals[i][j]; }

		// version with checking bounds
		Type  ElemAt(int i, int j) const
		{
			if (i < 0 || i >= RowNum() || j < 0 || j >= ColNum())
				throw MatrixAccessBoundsError("MatrixNM::ElemAt", i, j, RowNum(), ColNum());

			return _vals[i][j];
		}
		Type& ElemAt(int i, int j)
		{
			if (i < 0 || i >= RowNum() || j < 0 || j >= ColNum())
				throw MatrixAccessBoundsError("MatrixNM::ElemAt", i, j, RowNum(), ColNum());

			return _vals[i][j];
		}

		////////////////////            Arithmetic operators             ////////////////////
		MatrixNM operator-()            // unary minus
		{
			MatrixNM temp;
			for (size_t i = 0; i < RowNum(); i++)
				for (size_t j = 0; j < ColNum(); j++)
					temp._vals[i][j] = -_vals[i][j];
			return temp;
		}
		MatrixNM operator+(const MatrixNM& b) const
		{
			MatrixNM temp;
			for (size_t i = 0; i < RowNum(); i++)
				for (size_t j = 0; j < ColNum(); j++)
					temp._vals[i][j] = b._vals[i][j] + _vals[i][j];
			return temp;
		}
		MatrixNM operator-(const MatrixNM& b) const
		{
			MatrixNM temp;
			for (size_t i = 0; i < RowNum(); i++)
				for (size_t j = 0; j < ColNum(); j++)
					temp._vals[i][j] = _vals[i][j] - b._vals[i][j];
			return temp;
		}
		template<int K>
		MatrixNM<Type, N, K>  operator*(const MatrixNM<Type, M, K>& b) const
		{
			MatrixNM<Type, N, K>	ret;

			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++)
				{
					ret._vals[i][j] = 0;
					for (int k = 0; k < ColNum(); k++)
						ret._vals[i][j] += _vals[i][k] * b._vals[k][j];
				}

			return	ret;
		}

		MatrixNM operator*(Type b)
		{
			int	i, j;
			MatrixNM	ret(*this);

			for (i = 0; i < RowNum(); i++)
				for (j = 0; j < ColNum(); j++)
					ret._vals[i][j] *= b;

			return ret;
		}
		MatrixNM& operator*=(Type b)
		{
			int	i, j;

			for (i = 0; i < RowNum(); i++)
				for (j = 0; j < ColNum(); j++)
					_vals[i][j] *= b;

			return *this;
		}

		MatrixNM operator/(const Type& b) const
		{
			int	i, j;
			MatrixNM	ret(*this);

			for (i = 0; i < RowNum(); i++)
				for (j = 0; j < ColNum(); j++)
					ret._vals[i][j] /= b;

			return ret;
		}
		VectorN<Type, N> operator*(const VectorN<Type, M>& b) const
		{
			int	i, j;
			VectorN<Type, N>	ret;

			for (i = 0; i < N; i++)
			{
				ret[i] = 0;
				for (j = 0; j < M; j++)
					ret[i] += _vals[i][j] * b[j];
			}

			return ret;
		}

		friend MatrixNM operator*(Type a, const MatrixNM& b)
		{
			int	i, j;
			MatrixNM	ret;

			for (i = 0; i < b.RowNum(); i++)
				for (j = 0; j < b.ColNum(); j++)
					ret._vals[i][j] = a * b._vals[i][j];

			return ret;
		}
		friend VectorN<Type, M> operator*(const VectorN<Type, N>& a, const MatrixNM<Type, N, M>& b)
		{
			int	i, j;
			VectorN<Type, M>	ret;

			for (i = 0; i < M; i++)
			{
				ret[i] = 0;
				for (j = 0; j < N; j++)
					ret[i] += a[j] * b._vals[j][i];
			}

			return ret;
		}

		///////////////////////            Equality operations             //////////////////////
		bool operator==(const MatrixNM& b) const
		{
			for (int i = 0; i < N; i++)
				for (int j = 0; j < M; j++)
					if (_vals[i][j] != b._vals[i][j])
						return false;

			return true;
		}
		bool operator!=(const MatrixNM& b) const
		{
			return !(*this == b);
		}

		bool IsEqual(const MatrixNM& b, Type eps = Defaults::MatrixEqualityPrecision) const
		{
			for (int i = 0; i < RowNum(); i++)
				for (int j = 0; j < ColNum(); j++)
					if (std::abs(_vals[i][j] - b._vals[i][j]) > eps)
						return false;

			return true;
		}
		bool AreEqual(const MatrixNM& a, const MatrixNM& b, Type eps = Defaults::MatrixEqualityPrecision) const
		{
			return a.IsEqual(b, eps);
		}

		///////////////////            Trace, Inverse & Transpose             ///////////////////
		Type   Trace() const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("MatrixNM::Trace - must be square matrix", N, M, -1, -1);

			Type sum = 0;
			for (int i = 0; i < RowNum(); i++)
				sum += _vals[i][i];

			return sum;
		}

		void Invert()
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("MatrixNM::Invert - must be square matrix", N, M, -1, -1);

			MatrixNM& a = *this;
			MatrixNM<Type, N, 1>  b;      // dummy rhs

			b(0, 0) = 1.0;

			int i, icol, irow, j, k, l, ll;
			Type big, dum, pivinv;

			int n = RowNum();
			int m = b.ColNum();
			std::vector<int> indxc(n), indxr(n), ipiv(n);
			for (j = 0; j < n; j++) ipiv[j] = 0;
			for (i = 0; i < n; i++) {
				big = 0.0;
				for (j = 0; j < n; j++)
					if (ipiv[j] != 1)
						for (k = 0; k < n; k++) {
							if (ipiv[k] == 0) {
								if (std::abs(a._vals[j][k]) >= big) {
									big = std::abs(a._vals[j][k]);
									irow = j;
									icol = k;
								}
							}
						}
				++(ipiv[icol]);
				if (irow != icol) {
					for (l = 0; l < n; l++) std::swap(a._vals[irow][l], a._vals[icol][l]);
					for (l = 0; l < m; l++) std::swap(b._vals[irow][l], b._vals[icol][l]);
				}
				indxr[i] = irow;
				indxc[i] = icol;

				if (a._vals[icol][icol] == 0.0)
					throw SingularMatrixError("MatrixNM::Invert, gaussj: Singular Matrix");

				pivinv = 1.0 / a._vals[icol][icol];
				a._vals[icol][icol] = 1.0;
				for (l = 0; l < n; l++) a._vals[icol][l] *= pivinv;
				for (l = 0; l < m; l++) b._vals[icol][l] *= pivinv;
				for (ll = 0; ll < n; ll++)
					if (ll != icol) {
						dum = a._vals[ll][icol];
						a._vals[ll][icol] = 0.0;
						for (l = 0; l < n; l++) a._vals[ll][l] -= a._vals[icol][l] * dum;
						for (l = 0; l < m; l++) b._vals[ll][l] -= b._vals[icol][l] * dum;
					}
			}
			for (l = n - 1; l >= 0; l--) {
				if (indxr[l] != indxc[l])
					for (k = 0; k < n; k++)
						std::swap(a._vals[k][indxr[l]], a._vals[k][indxc[l]]);
			}
		}
		MatrixNM GetInverse() const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("MatrixNM::GetInverse - must be square matrix", N, M, -1, -1);

			MatrixNM a(*this);              // making a copy, where inverse will be stored at the end

			a.Invert();

			return a;
		}

		void Transpose()
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("MatrixNM::Transpose - inplace Transpose possible only for square  matrix", N, M, -1, -1);

			for (size_t i = 0; i < RowNum(); i++)
				for (size_t j = i + 1; j < ColNum(); j++)
					std::swap(_vals[i][j], _vals[j][i]);
		}
		MatrixNM<Type, M, N> GetTranspose() const
		{
			MatrixNM<Type, M, N> ret;

			for (size_t i = 0; i < ColNum(); i++)
				for (size_t j = 0; j < RowNum(); j++)
					ret._vals[i][j] = _vals[j][i];

			return ret;
		}

		///////////////////////////               I/O                 ///////////////////////////
		std::string to_string(int width, int precision) const
		{
			std::stringstream str;

			Print(str, width, precision);

			return str.str();
		}
		void Print(std::ostream& stream, int width, int precision) const
		{
			stream << "Rows: " << RowNum() << "  Cols: " << ColNum() << std::endl;

			/*	for (size_t i = 0; i < RowNum(); i++)
				{
					stream << "[ ";
					for (size_t j = 0; j < ColNum(); j++)
					{
						stream << std::setw(width) << std::setprecision(precision) << _vals[i][j] << ", ";
					}
					stream << " ]" << std::endl;
				}*/

			for (int i = 0; i < RowNum(); i++)
			{
				stream << "[ ";
				for (int j = 0; j < ColNum(); j++)
				{
					if (j == ColNum() - 1)
						stream << std::setw(width) << std::setprecision(precision) << _vals[i][j];
					else
						stream << std::setw(width) << std::setprecision(precision) << _vals[i][j] << ", ";
				}
				if (i == RowNum() - 1)
					stream << " ]";
				else
					stream << " ]" << std::endl;
			}
		}
		void   Print(std::ostream& stream, int width, int precision, Real zeroThreshold) const
		{
			stream << "Rows: " << RowNum() << " Cols: " << ColNum();

			for (int i = 0; i < RowNum(); i++)
			{
				stream << "[ ";
				for (int j = 0; j < ColNum(); j++)
				{
					Type value{ 0 };
					if (Abs(_vals[i][j]) > zeroThreshold)
						value = _vals[i][j];

					if (j == ColNum() - 1)
						stream << std::setw(width) << std::setprecision(precision) << value;
					else
						stream << std::setw(width) << std::setprecision(precision) << value << ", ";
				}
				if (i == RowNum() - 1)
					stream << " ]";
				else
					stream << " ]" << std::endl;
			}
		}
		friend std::ostream& operator<<(std::ostream& stream, const MatrixNM& a)
		{
			a.Print(stream, 10, 3);

			return stream;
		}
	};

	///////////////////               Default MatrixNM typdefs                 ///////////////////
	typedef MatrixNM<float, 2, 2> Matrix22Flt;
	typedef MatrixNM<float, 3, 3> Matrix33Flt;
	typedef MatrixNM<float, 4, 4> Matrix44Flt;

	typedef MatrixNM<Real, 2, 2> Matrix22Dbl;
	typedef MatrixNM<Real, 3, 3> Matrix33Dbl;
	typedef MatrixNM<Real, 4, 4> Matrix44Dbl;

	typedef MatrixNM<Complex, 2, 2> Matrix22Complex;
	typedef MatrixNM<Complex, 3, 3> Matrix33Complex;
	typedef MatrixNM<Complex, 4, 4> Matrix44Complex;

	typedef MatrixNM<float, 2, 2> Mat22F;
	typedef MatrixNM<float, 3, 3> Mat33F;
	typedef MatrixNM<float, 4, 4> Mat44F;

	typedef MatrixNM<Real, 2, 2> Mat22D;
	typedef MatrixNM<Real, 3, 3> Mat33D;
	typedef MatrixNM<Real, 4, 4> Mat44D;

	typedef MatrixNM<Complex, 2, 2> Mat22C;
	typedef MatrixNM<Complex, 3, 3> Mat33C;
	typedef MatrixNM<Complex, 4, 4> Mat44C;
}

///////////////////////////   ./include/base/MatrixSym.h   ///////////////////////////


namespace MML
{
	template<class Type>
	class MatrixSym
	{
	private:
		int  _dim;
		Type** _ptrData;

		void Init(int dim)
		{
			if (dim <= 0)
				throw MatrixDimensionError("MatrixSym::Init - dimension must be positive", dim, -1, -1, -1);

			_dim = dim;

			_ptrData = new Type * [dim];

			int numElem = (dim * dim + dim) / 2;
			if (_ptrData)
			{
				_ptrData[0] = numElem > 0 ? new Type[numElem] : nullptr;

				if (numElem > 0 && _ptrData[0] == nullptr)
					throw MatrixAllocationError("MatrixSym::Init - allocation error", dim, -1);

				for (int i = 1; i < dim; i++)
					_ptrData[i] = _ptrData[i - 1] + i;
			}
			else
				throw MatrixAllocationError("Matrix::Init - allocation error", dim, -1);
		}

	public:
		typedef Type value_type;      // make Type available externally

		///////////////////////          Constructors and destructor       //////////////////////
		MatrixSym() : _dim(0), _ptrData{ nullptr } {}
		MatrixSym(int dim) : _dim(dim)
		{
			Init(dim);
			for (int i = 0; i < _dim; ++i)
				for (int j = 0; j < _dim; ++j)
					_ptrData[i][j] = 0;
		}
		MatrixSym(int dim, Type val) : _dim(dim)
		{
			Init(dim);
			for (int i = 0; i < _dim; ++i)
				for (int j = 0; j <= i; ++j)
					_ptrData[i][j] = val;
		}
		MatrixSym(int dim, std::initializer_list<Type> values) : _dim(dim)
		{
			Init(dim);
			if (values.size() != (dim * dim + dim) / 2)
				throw MatrixDimensionError("Error in MatrixSym constructor: wrong dimensions", dim, -1, -1, -1);

			auto val = values.begin();
			for (int i = 0; i < _dim; ++i)
				for (int j = 0; j <= i; ++j)
					if (val != values.end())
						_ptrData[i][j] = *val++;
		}
		MatrixSym(const MatrixSym& m) : _dim(m._dim)
		{
			Init(m._dim);

			for (int i = 0; i < _dim; ++i)
				for (int j = 0; j <= i; ++j)
					_ptrData[i][j] = m._ptrData[i][j];
		}
		MatrixSym(MatrixSym&& m)
		{
			_ptrData = m._ptrData;

			_dim = m._dim;

			m._dim = 0;
			m._ptrData = nullptr;
		}
		~MatrixSym()
		{
			if (_ptrData != NULL) {
				delete[](_ptrData[0]);
				delete[](_ptrData);
			}
		}

		////////////////////////               Standard stuff             ////////////////////////
		int Dim() const { return (int)_dim; }
		int RowNum() const { return (int)_dim; }
		int ColNum() const { return (int)_dim; }

		bool IsEqual(const MatrixSym& b, Type eps = Defaults::MatrixEqualityPrecision) const
		{
			if (Dim() != b.Dim())
				return false;

			for (int i = 0; i < Dim(); i++)
				for (int j = 0; j <= i; j++)
				{
					if (Abs(_ptrData[i][j] - b._ptrData[i][j]) > eps)
						return false;
				}

			return true;
		}
		static bool AreEqual(const MatrixSym& a, const MatrixSym& b, Type eps = Defaults::MatrixEqualityPrecision)
		{
			return a.IsEqual(b, eps);
		}

		Matrix<Type> GetAsMatrix() const
		{
			Matrix<Type> ret(Dim(), Dim());

			for (int i = 0; i < Dim(); i++)
				for (int j = 0; j < Dim(); j++)
					ret[i][j] = (*this)(i, j);

			return ret;
		}

		/////////////////////          Init from regular Matrix           /////////////////////
		MatrixSym InitFromLower(Matrix<Type>& b)
		{
			if (b.RowNum() != b.ColNum())
				throw MatrixDimensionError("MatrixSym::InitFromLower - must be square matrix", b.RowNum(), b.ColNum(), -1, -1);

			MatrixSym ret(b.RowNum());
			for (int i = 0; i < b.RowNum(); i++)
				for (int j = 0; j <= i; j++)
					ret._ptrData[i][j] = b[i][j];

			return ret;
		}
		MatrixSym InitFromUpper(Matrix<Type>& b)
		{
			if (b.RowNum() != b.ColNum())
				throw MatrixDimensionError("MatrixSym::InitFromUpper - must be square matrix", b.RowNum(), b.ColNum(), -1, -1);

			MatrixSym ret(b.RowNum());
			for (int i = 0; i < b.RowNum(); i++)
				for (int j = i; j < b.RowNum(); j++)
					ret._ptrData[i][j] = b[i][j];

			return ret;
		}

		/////////////////////          Vector-Matrix conversion           /////////////////////
		Vector<Type> VectorFromRow(int rowInd)
		{
			if (rowInd >= RowNum())
				throw MatrixAccessBoundsError("VectorFromRow - row index must be less then a.RowNum()", rowInd, 0, RowNum(), ColNum());

			Vector<Type> ret(ColNum());
			for (int i = 0; i < ColNum(); i++)
				ret[i] = this->a(rowInd, i);

			return ret;
		}
		Vector<Type> VectorFromColumn(int colInd)
		{
			if (colInd >= ColNum())
				throw MatrixAccessBoundsError("VectorFromColumn - column index must be less then a.ColNum()", 0, colInd, RowNum(), ColNum());

			Vector<Type> ret(RowNum());
			for (int i = 0; i < RowNum(); i++)
				ret[i] = this->a(i, colInd);

			return ret;
		}
		Vector<Type> VectorFromDiagonal()
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("VectorFromDiagonal - must be square matrix", RowNum(), ColNum(), -1, -1);

			Vector<Type> ret(RowNum());
			for (int i = 0; i < RowNum(); i++)
				ret[i] = this->a(i, i);

			return ret;
		}

		///////////////////////////            Operators             ///////////////////////////
		MatrixSym& operator=(const MatrixSym& m)
		{
			if (this == &m)
				return *this;

			if (_dim != m._dim)
			{
				delete[](_ptrData[0]);
				delete[](_ptrData);

				Init(m.Dim());
			}

			for (size_t i = 0; i < m.Dim(); ++i)
				for (size_t j = 0; j <= i; ++j)
					_ptrData[i][j] = m._ptrData[i][j];

			return *this;
		}
		MatrixSym& operator=(MatrixSym&& m)
		{
			if (this == &m)
				return *this;

			std::swap(_ptrData, m._ptrData);
			std::swap(_dim, m._dim);

			return *this;
		}

		////////////////////////           Access operators             ///////////////////////
		Type  operator()(int i, int j) const {
			if (i < j)
				return _ptrData[j][i];
			else
				return _ptrData[i][j];
		}
		Type& operator()(int i, int j) {
			if (i < j)
				return _ptrData[j][i];
			else
				return _ptrData[i][j];
		}

		// version with checking bounds
		Type  ElemAt(int i, int j) const
		{
			if (i < 0 || i >= RowNum() || j < 0 || j >= ColNum())
				throw MatrixAccessBoundsError("MatrixSym::ElemAt", i, j, RowNum(), ColNum());

			return operator()(i, j);
		}
		Type& ElemAt(int i, int j)
		{
			if (i < 0 || i >= RowNum() || j < 0 || j >= ColNum())
				throw MatrixAccessBoundsError("MatrixSym::ElemAt", i, j, RowNum(), ColNum());

			return operator()(i, j);
		}

		///////////////////////           Arithmetic operators             //////////////////////
		MatrixSym operator+(const MatrixSym& b) const
		{
			if (_dim != b._dim)
				throw MatrixDimensionError("MatrixSym::operator+() - must be same dim", _dim, -1, b._dim, -1);

			MatrixSym temp(_dim);
			for (size_t i = 0; i < Dim(); i++)
				for (size_t j = 0; j <= i; j++)
					temp._ptrData[i][j] = b._ptrData[i][j] + _ptrData[i][j];

			return temp;
		}
		MatrixSym operator-()
		{
			MatrixSym temp(_dim);
			for (size_t i = 0; i < Dim(); i++)
				for (size_t j = 0; j <= i; j++)
					temp._ptrData[i][j] = -_ptrData[i][j];

			return temp;
		}
		MatrixSym operator-(const MatrixSym& b) const
		{
			if (_dim != b._dim)
				throw MatrixDimensionError("MatrixSym::operator-() - must be same dim", _dim, -1, b._dim, -1);

			MatrixSym temp(_dim);
			for (size_t i = 0; i < Dim(); i++)
				for (size_t j = 0; j <= i; j++)
					temp._ptrData[i][j] = b._ptrData[i][j] - _ptrData[i][j];

			return temp;
		}

		Matrix<Type> operator*(const MatrixSym& b) const
		{
			if (Dim() != b.Dim())
				throw MatrixDimensionError("MatrixSym::operator*(MatrixSym &) - a.colNum must be equal to b.rowNum", _dim, _dim, b._dim, b._dim);

			Matrix<Type>	ret(RowNum(), b.ColNum());
			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++)
				{
					ret[i][j] = 0;
					for (int k = 0; k < ColNum(); k++)
						ret[i][j] += (*this)(i, k) * b(k, j);
				}

			return	ret;
		}

		Matrix<Type> operator*(const Matrix<Type>& b) const
		{
			if (Dim() != b.RowNum())
				throw MatrixDimensionError("MatrixSym::operator*(Matrix &) - a.colNum must be equal to b.rowNum", _dim, _dim, b.RowNum(), b.ColNum());

			Matrix<Type>	ret(RowNum(), b.ColNum());
			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++)
				{
					ret[i][j] = 0;
					for (int k = 0; k < ColNum(); k++)
						ret[i][j] += (*this)(i, k) * b(k, j);
				}

			return	ret;
		}

		friend MatrixSym operator*(const MatrixSym& a, Type b)
		{
			int	i, j;
			MatrixSym	ret(a.Dim());

			for (i = 0; i < a.Dim(); i++)
				for (j = 0; j <= i; j++)
					ret[i][j] = a._ptrData[i][j] * b;

			return ret;
		}
		friend MatrixSym operator*(Type a, const MatrixSym& b)
		{
			int	i, j;
			MatrixSym	ret(a.Dim());

			for (i = 0; i < a.Dim(); i++)
				for (j = 0; j <= i; j++)
					ret[i][j] = a * b._ptrData[i][j];

			return ret;
		}
		friend MatrixSym operator/(const MatrixSym& a, Type b)
		{
			int	i, j;
			MatrixSym	ret(a.Dim());

			for (i = 0; i < a.Dim(); i++)
				for (j = 0; j <= i; j++)
					ret[i][j] = a._ptrData[i][j] / b;

			return ret;
		}

		friend Vector<Type> operator*(const MatrixSym& a, const Vector<Type>& b)
		{
			if (a.Dim() != b.size())
				throw MatrixDimensionError("operator*(MatSym a, Vec b) - a.Dim must be equal to vector size", a.Dim(), a.Dim(), (int)b.size(), -1);

			Vector<Type>	ret(a.RowNum());
			for (int i = 0; i < a.Dim(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < a.Dim(); j++)
					ret[i] += a(i, j) * b[j];
			}

			return ret;
		}
		friend Vector<Type> operator*(const Vector<Type>& a, const MatrixSym& b)
		{
			if (a.size() != b.Dim())
			{
				//std::string error = std::format("Hello {}!\n", "world");
				throw MatrixDimensionError("operator*(Vec a, MatSym b) - vector size must be equal to b.Dim", (int)a.size(), -1, b.Dim(), -1);
			}

			Vector<Type>	ret(b.Dim());
			for (int i = 0; i < b.Dim(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < b.Dim(); j++)
					ret[i] += a[j] * b(i, j);
			}

			return ret;
		}

		////////////////////////            Inverse              ///////////////////////
		Matrix<Type> GetInverse() const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("MatrixSym::GetInverse - must be square matrix", _dim, _dim, -1, -1);

			Matrix<Type> a = this->GetAsMatrix();              // making a copy, where inverse will be stored at the end

			a.Invert();

			return a;
		}

		///////////////////////////               I/O                 ///////////////////////////
		std::string to_string(int width, int precision) const
		{
			std::stringstream str;

			Print(str, width, precision);

			return str.str();
		}
		void   Print(std::ostream& stream, int width, int precision) const
		{
			stream << "Rows: " << Dim() << std::endl;

			for (int i = 0; i < Dim(); i++)
			{
				stream << "[ ";
				for (int j = 0; j < Dim(); j++)
				{
					stream << std::setw(width) << std::setprecision(precision) << (*this)(i, j) << ", ";
				}
				stream << " ]" << std::endl;
			}
		}
		friend std::ostream& operator<<(std::ostream& stream, const MatrixSym& a)
		{
			a.Print(stream, 10, 3);

			return stream;
		}
	};
}
///////////////////////////   ./include/base/MatrixTriDiag.h   ///////////////////////////

namespace MML
{
	template<class Type>
	class TridiagonalMatrix
	{
	private:
		int _dim;
		Vector<Type> _belowDiag;
		Vector<Type> _diag;
		Vector<Type> _aboveDiag;

	public:
		TridiagonalMatrix(int dim, const Vector<Type>& a, const Vector<Type>& diag, const Vector<Type>& c) : _belowDiag(a), _diag(diag), _aboveDiag(c), _dim(dim)
		{
			if( dim < 3)
				throw("Error in TridiagonalMatrix constructor: dim < 3");

			if (a.size() != dim || diag.size() != dim || c.size() != dim)
				throw("Error in TridiagonalMatrix constructor: wrong dimensions");
		}

		TridiagonalMatrix(int dim, std::initializer_list<Type> values) : _dim(dim), _belowDiag(dim), _diag(dim), _aboveDiag(dim)
		{
			if (values.size() != dim * 3 - 2)
				throw("Error in TridiagonalMatrix constructor: wrong dimensions");

			auto val = values.begin();
			_belowDiag[0] = 0.0;
			_diag[0] = *val++;
			_aboveDiag[0] = *val++;
			for (int i = 1; i < dim - 1; ++i)
			{
				_belowDiag[i] = *val++;
				_diag[i] = *val++;
				_aboveDiag[i] = *val++;
			}
			_belowDiag[dim - 1] = *val++;
			_diag[dim - 1] = *val++;
			_aboveDiag[dim - 1] = 0.0;
		}

		int RowNum() const { return _dim; }
		int ColNum() const { return _dim; }

		Type  operator()(int i, int j) const {
			if (i == j)
				return _diag[i];
			else if (i == j - 1)
				return _aboveDiag[i];
			else if (i == j + 1 && j < _dim - 1)
				return _belowDiag[i];
			else
				return 0.0;
		}
		Type& operator()(int i, int j) {
			if (i == j)
				return _diag[i];
			else if (i == j - 1)
				return _aboveDiag[i];
			else if (i == j + 1 && j < _dim - 1)
				return _belowDiag[i];
			else
				throw MatrixAccessBoundsError("TridiagonalMatrix::operator()", i, j, _dim, _dim);
		}

		// TODO 0.9 - dodati IsEqual
		// TODO 1.0 - imaju li smisla operacije s regularnim matricama?
		TridiagonalMatrix operator+(const TridiagonalMatrix& b) const
		{
			if (_dim != b._dim)
				throw MatrixDimensionError("TridiagonalMatrix::operator+() - must be same dim", _dim, _dim, b._dim, b._dim);

			TridiagonalMatrix temp(*this);
			for (int i = 0; i < _belowDiag.size(); i++)
				temp._belowDiag[i] += b._belowDiag[i];
			for (int i = 0; i < _diag.size(); i++)
				temp._diag[i] += b._diag[i];
			for (int i = 0; i < _aboveDiag.size(); i++)
				temp._aboveDiag[i] += b._aboveDiag[i];

			return temp;
		}
		TridiagonalMatrix operator-(const TridiagonalMatrix& b) const
		{
			if (_dim != b._dim)
				throw MatrixDimensionError("TridiagonalMatrix::operator+() - must be same dim", _dim, _dim, b._dim, b._dim);

			TridiagonalMatrix temp(*this);
			for (int i = 0; i < _belowDiag.size(); i++)
				temp._belowDiag[i] -= b._belowDiag[i];
			for (int i = 0; i < _diag.size(); i++)
				temp._diag[i] -= b._diag[i];
			for (int i = 0; i < _aboveDiag.size(); i++)
				temp._aboveDiag[i] -= b._aboveDiag[i];

			return temp;
		}

		friend TridiagonalMatrix operator*(const TridiagonalMatrix& a, Type b)
		{
			int	i, j;
			TridiagonalMatrix	ret(a);

			for (int i = 0; i < ret._a.size(); i++)
				ret._a[i] *= b;
			for (int i = 0; i < ret._diag.size(); i++)
				ret._diag[i] *= b;
			for (int i = 0; i < ret._c.size(); i++)
				ret._c[i] *= b;

			return ret;
		}
		friend TridiagonalMatrix operator*(Type a, const TridiagonalMatrix& b)
		{
			int	i, j;
			TridiagonalMatrix	ret(b);

			for (int i = 0; i < ret._a.size(); i++)
				ret._a[i] *= a;
			for (int i = 0; i < ret._diag.size(); i++)
				ret._diag[i] *= a;
			for (int i = 0; i < ret._c.size(); i++)
				ret._c[i] *= a;

			return ret;
		}
		friend TridiagonalMatrix operator/(const TridiagonalMatrix& a, Type b)
		{
			int	i, j;
			TridiagonalMatrix	ret(a);

			for (int i = 0; i < ret._a.size(); i++)
				ret._a[i] /= b;
			for (int i = 0; i < ret._diag.size(); i++)
				ret._diag[i] /= b;
			for (int i = 0; i < ret._c.size(); i++)
				ret._c[i] /= b;

			return ret;
		}

		// TODO - invert
		// TODO - 0.9 transpose

		void Solve(Vector<Type>& rhs, Vector<Type>& sol)
		{
			int j, n = _belowDiag.size();
			Real bet;
			Vector<Type> gam(n);

			if (_diag[0] == 0.0)
				throw("Error 1 in tridag");

			sol[0] = rhs[0] / (bet = _diag[0]);

			for (j = 1; j < n; j++) {
				gam[j] = _aboveDiag[j - 1] / bet;
				bet = _diag[j] - _belowDiag[j] * gam[j];

				if (bet == 0.0)
					throw("Error 2 in tridag");

				sol[j] = (rhs[j] - _belowDiag[j] * sol[j - 1]) / bet;
			}
			for (j = (n - 2); j >= 0; j--)
				sol[j] -= gam[j + 1] * sol[j + 1];
		}
		Vector<Type> Solve(Vector<Type>& rhs)
		{
			Vector<Type> sol(rhs.size());
			Solve(rhs, sol);
			return sol;
		}

		void   Print(std::ostream& stream, int width, int precision) const
		{
			stream << "Dim: " << _dim << std::endl;
			for (int i = 0; i < _dim; i++)
			{
				stream << "[ ";
				for (int j = 0; j < _dim; j++) {
					stream << std::setw(width) << std::setprecision(precision) << (*this)(i, j) << ", ";
				}
				stream << " ]" << std::endl;
			}
		}
	};
}

///////////////////////////   ./include/interfaces/ITensor.h   ///////////////////////////

namespace MML
{
	enum TensorIndexType { CONTRAVARIANT, COVARIANT };
    
	template<int N>
	class ITensor2
	{
	public:
		virtual int   NumContravar() const = 0;
		virtual int   NumCovar() const = 0;

		virtual Real  operator()(int i, int j) const = 0;
		virtual Real& operator()(int i, int j) = 0;
	};

	template<int N>
	class ITensor3
	{
	public:
		virtual int   NumContravar() const = 0;
		virtual int   NumCovar() const = 0;

		virtual Real  operator()(int i, int j, int k) const = 0;
		virtual Real& operator()(int i, int j, int k) = 0;
	};

	template<int N>
	class ITensor4
	{
	public:
		virtual int   NumContravar() const = 0;
		virtual int   NumCovar() const = 0;

		virtual Real  operator()(int i, int j, int k, int l) const = 0;
		virtual Real& operator()(int i, int j, int k, int l) = 0;
	};

	template<int N>
	class ITensor5
	{
	public:
		virtual int   NumContravar() const = 0;
		virtual int   NumCovar() const = 0;

		virtual Real  operator()(int i, int j, int k, int l, int m) const = 0;
		virtual Real& operator()(int i, int j, int k, int l, int m) = 0;
	};
}


///////////////////////////   ./include/base/Tensor.h   ///////////////////////////



namespace MML
{
	template <int N>
	class Tensor2 : public ITensor2<N>
	{
		MatrixNM<Real, N, N> _coeff = { 0 };
	public:
		int _numContravar = 0;
		int _numCovar = 0;
		bool _isContravar[2];

		Tensor2(int nCovar, int nContraVar) : _numCovar(nCovar), _numContravar(nContraVar)
		{
			if ( _numContravar < 0 || _numCovar  < 0 || _numContravar + _numCovar != 2 )
				throw TensorCovarContravarNumError("Tensor2 ctor, wrong number of contravariant and covariant indices", nCovar, nContraVar);

			for (int i = 0; i < _numCovar; i++)
				_isContravar[i] = false;

			for (int i = _numCovar; i < _numCovar + _numContravar; i++)
				_isContravar[i] = true;
		}
		Tensor2(int nCovar, int nContraVar, std::initializer_list<Real> values) : _numCovar(nCovar), _numContravar(nContraVar)
		{
			if ( _numContravar < 0 || _numCovar  < 0 || _numContravar + _numCovar != 2 )
				throw TensorCovarContravarNumError("Tensor2 ctor, wrong number of covariant and contravariant indices", nCovar, nContraVar);

			for (int i = 0; i < _numCovar; i++)
				_isContravar[i] = false;

			for (int i = _numCovar; i < _numCovar + _numContravar; i++)
				_isContravar[i] = true;

			auto val = values.begin();
			for (size_t i = 0; i < N; ++i)
				for (size_t j = 0; j < N; ++j)
					if (val != values.end())
						_coeff[i][j] = *val++;
					else
						_coeff[i][j] = 0.0;
		}

		int   NumContravar() const { return _numContravar; }
		int   NumCovar()     const { return _numCovar; }

		// add IsContravar() method
		bool IsContravar(int i) const { return _isContravar[i]; }

		// add IsCovar() method
		bool IsCovar(int i) const { return !_isContravar[i]; }

		Real  operator()(int i, int j) const { return _coeff[i][j]; }
		Real& operator()(int i, int j) { return _coeff[i][j]; }

		//Real  Component(int i, int j) const { return _coeff[i][j]; }
		//Real& Component(int i, int j) { return _coeff[i][j]; }

		Tensor2 operator+(const Tensor2& other) const
		{
			if (_numContravar != other._numContravar || _numCovar != other._numCovar)
				throw TensorCovarContravarAirthmeticError("Tensor2 operator+, wrong number of contravariant and covariant indices", _numContravar, _numCovar, other._numContravar, other._numCovar);

			Tensor2 result(_numContravar, _numCovar);

			result._coeff = _coeff + other._coeff;

			return result;
		}

		Tensor2 operator-(const Tensor2& other) const
		{
			if (_numContravar != other._numContravar || _numCovar != other._numCovar)
				throw TensorCovarContravarAirthmeticError("Tensor2 operator-, wrong number of contravariant and covariant indices", _numContravar, _numCovar, other._numContravar, other._numCovar);

			Tensor2 result(_numContravar, _numCovar);

			result._coeff = _coeff - other._coeff;

			return result;
		}

		Tensor2 operator*=(Real scalar) const
		{
			Tensor2 result(_numContravar, _numCovar);

			result._coeff = _coeff * scalar;

			return result;
		}

		Tensor2 operator/(Real scalar) const
		{
			Tensor2 result(_numContravar, _numCovar);

			result._coeff = _coeff / scalar;

			return result;
		}

		friend Tensor2 operator*(Real scalar, const Tensor2& b)
		{
			Tensor2 result(b.NumContravar(), b.NumCovar());

			result._coeff = b._coeff * scalar;

			return result;
		}

		Real Contract() const
		{
			if (_numContravar != 1 || _numCovar != 1)
				throw TensorCovarContravarNumError("Tensor2 Contract, wrong number of contravariant and covariant indices", _numContravar, _numCovar);

			Real result = 0.0;
			for (int i = 0; i < N; i++)
				result += _coeff[i][i];

			return result;
		}
		Real operator()(const VectorN<Real, N>& v1, const VectorN<Real, N>& v2) const
		{
			Real sum = 0.0;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					sum += _coeff[i][j] * v1[i] * v2[j];

			return sum;
		}

		void   Print(std::ostream& stream, int width, int precision) const
		{
			stream << std::fixed << "N = " << N << std::endl;

			for (size_t i = 0; i < N; i++)
			{
				stream << "[ ";
				for (size_t j = 0; j < N; j++)
					stream << std::setw(width) << std::setprecision(precision) << _coeff[i][j] << ", ";
				stream << " ]" << std::endl;
			}
		}
		friend std::ostream& operator<<(std::ostream& stream, const Tensor2& a)
		{
			a.Print(stream, 15, 10);

			return stream;
		}
	};

	template <int N>
	class Tensor3 : public ITensor3<N>
	{
		Real _coeff[N][N][N] = { 0 };
	public:
		int _numContravar;
		int _numCovar;
		bool _isContravar[3];

		Tensor3(int nCovar, int nContraVar) : _numCovar(nCovar), _numContravar(nContraVar)
		{
			if (_numContravar + _numCovar != 3)
				throw TensorCovarContravarNumError("Tensor3 ctor, wrong number of contravariant and covariant indices", nCovar, nContraVar);

			for (int i = 0; i < _numCovar; i++)
				_isContravar[i] = false;

			for (int i = _numCovar; i < _numCovar + _numContravar; i++)
				_isContravar[i] = true;
		}

		int   NumContravar() const { return _numContravar; }
		int   NumCovar()     const { return _numCovar; }

		Real  operator()(int i, int j, int k) const { return _coeff[i][j][k]; }
		Real& operator()(int i, int j, int k) { return _coeff[i][j][k]; }

		Real operator()(const VectorN<Real, N>& v1, const VectorN<Real, N>& v2, const VectorN<Real, N>& v3) const
		{
			Real sum = 0.0;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
						sum += _coeff[i][j][k] * v1[i] * v2[j] * v3[k];

			return sum;
		}

		VectorN<Real, N> Contract(int ind1, int ind2) const
		{
			VectorN<Real, N> result;

			if (ind1 < 0 || ind1 >= N || ind2 < 0 || ind2 >= N)
				throw TensorIndexError("Tensor3 Contract, wrong index");

			if (ind1 == ind2)
				throw TensorIndexError("Tensor3 Contract, indices are the same");

			if (ind1 == 0 && ind2 == 1)
			{
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						result[i] += _coeff[j][j][i];
			}
			else if (ind1 == 0 && ind2 == 2)
			{
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						result[i] += _coeff[j][i][j];
			}
			else if (ind1 == 1 && ind2 == 2)
			{
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						result[i] += _coeff[i][j][j];
			}
			else
				throw TensorIndexError("Tensor3 Contract, wrong indices");

			return result;
		}
	};

	template <int N>
	class Tensor4 : public ITensor4<N>
	{
		Real _coeff[N][N][N][N] = { 0 };
	public:
		int _numContravar;
		int _numCovar;
		bool _isContravar[4];

		Tensor4(int nCovar, int nContraVar) : _numCovar(nCovar), _numContravar(nContraVar) 
		{
			if (_numContravar + _numCovar != 4)
				throw TensorCovarContravarNumError("Tensor4 ctor, wrong number of contravariant and covariant indices", nCovar, nContraVar);

			for (int i = 0; i < _numCovar; i++)
				_isContravar[i] = false;

			for (int i = _numCovar; i < _numCovar + _numContravar; i++)
				_isContravar[i] = true;
		}

		int   NumContravar() const { return _numContravar; }
		int   NumCovar()     const { return _numCovar; }

		Real  operator()(int i, int j, int k, int l) const { return _coeff[i][j][k][l]; }
		Real& operator()(int i, int j, int k, int l) { return _coeff[i][j][k][l]; }

		Real operator()(const VectorN<Real, N>& v1, const VectorN<Real, N>& v2, const VectorN<Real, N>& v3, const VectorN<Real, N>& v4) const
		{
			Real sum = 0.0;
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
						for (int l = 0; l < N; l++)
							sum += _coeff[i][j][k][l] * v1[i] * v2[j] * v3[k] * v4[l];

			return sum;
		}

		Tensor2<N> Contract(int ind1, int ind2) const
		{
			MatrixNM<Real, N, N> result;
			Tensor2<N> ret;

			// check covar contravar!!!!
			if (ind1 < 0 || ind1 >= N || ind2 < 0 || ind2 >= N)
				throw TensorIndexError("Tensor4 Contract, wrong index");

			if (ind1 == ind2)
				throw TensorIndexError("Tensor4 Contract, indices are the same");

			if (ind1 == 0 && ind2 == 1)
			{
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						for(int k = 0; k < N; k++)
							result[i][j] += _coeff[k][k][i][j];
			}
			else if (ind1 == 0 && ind2 == 2)
			{
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						for(int k = 0; k < N; k++)
							result[i][j] += _coeff[k][i][k][j];
			}
			else if (ind1 == 0 && ind2 == 3)
			{
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						for(int k = 0; k < N; k++)
							result[i][j] += _coeff[k][i][j][k];
			}
			else if (ind1 == 1 && ind2 == 2)
			{
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						for(int k = 0; k < N; k++)
							result[i][j] += _coeff[i][k][k][j];
			}
			else if (ind1 == 1 && ind2 == 3)
			{
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						for(int k = 0; k < N; k++)
							result[i][j] += _coeff[i][k][j][k];
			}
			else if (ind1 == 2 && ind2 == 3)
			{
				for (int i = 0; i < N; i++)
					for (int j = 0; j < N; j++)
						for(int k = 0; k < N; k++)
							result[i][j] += _coeff[i][j][k][k];
			}
			else
				throw TensorIndexError("Tensor4 Contract, wrong indices");

			ret._coeff = result;

			return ret;
		}
	};

	template <int N>
	class Tensor5 : public ITensor5<N>
	{
		Real _coeff[N][N][N][N][N];
	public:
		int _numContravar;
		int _numCovar;
		bool _isContravar[5];

		Tensor5(int nCovar, int nContraVar) : _numCovar(nCovar), _numContravar(nContraVar) 
		{
			if (_numContravar + _numCovar != 5)
				throw TensorCovarContravarNumError("Tensor5 ctor, wrong number of contravariant and covariant indices", nCovar, nContraVar);

			for (int i = 0; i < _numCovar; i++)
				_isContravar[i] = false;

			for (int i = _numCovar; i < _numCovar + _numContravar; i++)
				_isContravar[i] = true;

		}

		int   NumContravar() const { return _numContravar; }
		int   NumCovar()     const { return _numCovar; }

		Real  operator()(int i, int j, int k, int l, int m) const { return _coeff[i][j][k][l][m]; }
		Real& operator()(int i, int j, int k, int l, int m) { return _coeff[i][j][k][l][m]; }
	};
}
///////////////////////////   ./include/interfaces/IFunction.h   ///////////////////////////


namespace MML
{
	//////////////////////////////////////////////////////////////////////
	template<typename _RetType, typename _ArgType>
	class IFunction
	{
	public:
		virtual _RetType operator()(_ArgType) const = 0;

		virtual ~IFunction() {}
	};

	//////////////////////////////////////////////////////////////////////
	class IRealFunction : public IFunction<Real, Real>
	{
	public:
		virtual Real operator()(Real) const = 0;
		
		virtual ~IRealFunction() {}

		void GetValues(Real x1, Real x2, int numPnt, Vector<Real>& outX, Vector<Real>& outY)
		{
			outX.resize(numPnt);
			outY.resize(numPnt);

			for (int i = 0; i < numPnt; i++) {
				outX[i] = x1 + i * (x2 - x1) / (numPnt - 1);
				outY[i] = (*this)(outX[i]);
			}
		}
		// ovo u console printer!
		void Print(Real x1, Real x2, int numPnt)
		{
			for (int i = 0; i < numPnt; i++) {
				Real x = x1 + i * (x2 - x1) / (numPnt - 1);
				std::cout << x << " " << (*this)(x) << std::endl;
			}
		}
	};

	// RFExt - zna derivaciju
	class IRealFunctionParametrized : public IRealFunction
	{
	public:
		virtual int		getNumParam() const = 0;
		virtual Real	getParam(int i) const = 0;
		virtual void	setParam(int i, Real val) = 0;

		virtual Vector<Real>	getParams() const = 0;
		virtual void					setParams(const Vector<Real>&) = 0;

	};

	//////////////////////////////////////////////////////////////////////
	template<int N>
	class IScalarFunction : public IFunction<Real, const VectorN<Real, N>&>
	{
	public:
		virtual Real operator()(const VectorN<Real, N>& x) const = 0;

		virtual ~IScalarFunction() {}
	};

	template<int N>
	class IScalarFunctionParametrized : public IScalarFunction<N>
	{
	public:
		virtual int		getNumParam() const = 0;
		virtual Real	getParam(int i) const = 0;
		virtual void	setParam(int i, Real val) = 0;

		virtual Vector<Real>	getParams() const = 0;
		virtual void					setParams(const Vector<Real>&) = 0;

	};

	//////////////////////////////////////////////////////////////////////
	template<int N>
	class IRealToVectorFunction : public IFunction<VectorN<Real, N>, Real>
	{
	public:
		virtual VectorN<Real, N> operator()(Real x) const = 0;

		virtual ~IRealToVectorFunction() {}
	};

	//////////////////////////////////////////////////////////////////////
	template<int N>
	class IVectorFunction : public IFunction<VectorN<Real, N>, const VectorN<Real, N>&>
	{
	public:
		virtual VectorN<Real, N> operator()(const VectorN<Real, N>& x) const = 0;

		virtual Real operator()(const VectorN<Real, N>& x, int component) const
		{
			VectorN<Real, N> val = (*this)(x);
			return val[component];
		}

		virtual ~IVectorFunction() {}
	};

	template<int N>
	class IVectorFunctionParametrized : public IVectorFunction<N>
	{
	public:
		virtual int		getNumParam() const = 0;
		virtual Real	getParam(int i) const = 0;
		virtual void	setParam(int i, Real val) = 0;

		virtual Vector<Real>	getParams() const = 0;
		virtual void					setParams(const Vector<Real>&) = 0;

	};

	//////////////////////////////////////////////////////////////////////
	template<int N, int M>
	class IVectorFunctionNM : public IFunction<VectorN<Real, M>, const VectorN<Real, N>&>
	{
	public:
		virtual VectorN<Real, M> operator()(const VectorN<Real, N>& x) const = 0;
		
		virtual Real operator()(const VectorN<Real, N>& x, int component) const
		{
			VectorN<Real, M> val = (*this)(x);
			return val[component];
		}
		
		virtual ~IVectorFunctionNM() {}
	};

	//////////////////////////////////////////////////////////////////////
	template<int N>
	class IParametricCurve : public IRealToVectorFunction<N>
	{
	public:
		virtual VectorN<Real, N> operator()(Real t) const = 0;

		virtual Real getMinT() const = 0;
		virtual Real getMaxT() const = 0;

		std::vector<VectorN<Real, N>> GetTrace(double t1, double t2, int numPoints) const
		{
			std::vector<VectorN<Real, N>> ret;
			double deltaT = (t2 - t1) / (numPoints - 1);
			for (Real t = t1; t <= t2; t += deltaT)
				ret.push_back((*this)(t));
			return ret;
		}

		virtual ~IParametricCurve() {}
	};

	template<int N>
	class IParametricCurveParametrized : public IParametricCurve<N>
	{
	public:
			virtual int		getNumParam() const = 0;
			virtual Real	getParam(int i) const = 0;
			virtual void	setParam(int i, Real val) = 0;

			// is this neccessary?
			//virtual Vector<Real>	getParams() const = 0;
			//virtual void					setParams(const Vector<Real>&) = 0;
	};

	//////////////////////////////////////////////////////////////////////
	// simple regular surface, defined on rectangular coordinate patch
	template<int N>
	class IParametricSurface : public IFunction<VectorN<Real, N>, const VectorN<Real, 2>&>
	{
	public:
		virtual VectorN<Real, N> operator()(Real u, Real w) const = 0;

		virtual Real getMinU() const = 0;
		virtual Real getMaxU() const = 0;
		virtual Real getMinW() const = 0;
		virtual Real getMaxW() const = 0;

		virtual VectorN<Real, N> operator()(const VectorN<Real, 2>& coord) const
		{
			return operator()(coord[0], coord[1]);
		}

		virtual ~IParametricSurface() {}
	};

	// complex surface, with fixed u limits, but variable w limits (dependent on u)
	template<int N>
	class IParametricSurfaceComplex : public IFunction<VectorN<Real, N>, const VectorN<Real, 2>&>
	{
	public:
		virtual VectorN<Real, N> operator()(Real u, Real w) const = 0;

		virtual Real getMinU() const = 0;
		virtual Real getMaxU() const = 0;
		virtual Real getMinW(Real u) const = 0;
		virtual Real getMaxW(Real u) const = 0;

		virtual VectorN<Real, N> operator()(const VectorN<Real, 2>& coord) const
		{
			return operator()(coord[0], coord[1]);
		}

		virtual ~IParametricSurfaceComplex() {}
	};

	//////////////////////////////////////////////////////////////////////
	template<int N>
	class ITensorField2 : public IFunction<Tensor2<N>, const VectorN<Real, N>& >
	{
		int _numContravar;
		int _numCovar;
	public:
		ITensorField2(int numContra, int numCo) : _numContravar(numContra), _numCovar(numCo) { }

		int getNumContravar() const { return _numContravar; }
		int getNumCovar() const { return _numCovar; }

		virtual Real    Component(int i, int j, const VectorN<Real, N>& pos) const = 0;
		virtual ~ITensorField2() {}
	};

	template<int N>
	class ITensorField3 : public IFunction<Tensor3<N>, const VectorN<Real, N>& >
	{
		int _numContravar;
		int _numCovar;
	public:
		int getNumContravar() const { return _numContravar; }
		int getNumCovar() const { return _numCovar; }

		virtual Real    Component(int i, int j, int k, const VectorN<Real, N>& pos) const = 0;
		virtual ~ITensorField3() {}
	};

	template<int N>
	class ITensorField4 : public IFunction<Tensor4<N>, const VectorN<Real, N>& >
	{
		int _numContravar;
		int _numCovar;
	public:
		int getNumContravar() const { return _numContravar; }
		int getNumCovar() const { return _numCovar; }

		virtual Real    Component(int i, int j, int k, int l, const VectorN<Real, N>& pos) const = 0;
		virtual ~ITensorField4() {}
	};

	template<int N>
	class ITensorField5 : public IFunction<Tensor5<N>, const VectorN<Real, N>& >
	{
		int _numContravar;
		int _numCovar;
	public:
		int getNumContravar() const { return _numContravar; }
		int getNumCovar() const { return _numCovar; }

		virtual Real    Component(int i, int j, int k, int l, int m, const VectorN<Real, N>& pos) const = 0;
		virtual ~ITensorField5() {}
	};
}
///////////////////////////   ./include/interfaces/ICoordTransf.h   ///////////////////////////



namespace MML
{
	template<typename VectorFrom, typename VectorTo, int N>
	class ICoordTransf
	{
	public:
		virtual       VectorTo            transf(const VectorFrom& in) const = 0;
		virtual const IScalarFunction<N>& coordTransfFunc(int i) const = 0;

		virtual ~ICoordTransf() {}
	};

	template<typename VectorFrom, typename VectorTo, int N>
	class ICoordTransfWithInverse : public virtual ICoordTransf<VectorFrom, VectorTo, N>
	{
	public:
		virtual       VectorFrom          transfInverse(const VectorTo& in) const = 0;
		virtual const IScalarFunction<N>& inverseCoordTransfFunc(int i) const = 0;

		virtual ~ICoordTransfWithInverse() {}
	};
}
///////////////////////////   ./include/interfaces/IODESystem.h   ///////////////////////////


namespace MML
{
	class IODESystem
	{
	public:
		virtual int     getDim() const = 0;
		virtual void    derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const = 0;
	};

	class IODESystemParametrized : public IODESystem
	{
	public:
		virtual int		getNumParam() const = 0;
		virtual Real	getParam(int i) const = 0;
		virtual void	setParam(int i, Real val) = 0;

		virtual Vector<Real>	getParams() const = 0;
		virtual void					setParams(const Vector<Real>&) = 0;

	};

}
///////////////////////////   ./include/base/Polynom.h   ///////////////////////////



namespace MML
{
	template <typename _Field>
	class Polynom
	{
	protected:
		// polynom coefficients - with _vecCoef[0] being the constant term
		std::vector<Real> _vecCoef;
	public:
		Polynom() {}
		Polynom(int n) { _vecCoef.resize(n + 1); }
		Polynom(const std::vector<Real>& vecCoef) : _vecCoef(vecCoef) {}
		Polynom(std::initializer_list<Real> list) : _vecCoef(list) {}
		Polynom(const Polynom& Copy) = default;

		Polynom& operator=(const Polynom& Copy) { _vecCoef = Copy._vecCoef; return *this; }

		int  GetDegree() const { return (int)_vecCoef.size() - 1; }
		void SetDegree(int newDeg) { _vecCoef.resize(newDeg + 1); }
		bool IsNullPolynom() const { return _vecCoef.size() == 0; }

		Real  operator[] (int i) const { return _vecCoef[i]; }
		Real& operator[] (int i) { return _vecCoef[i]; }

		_Field operator() (const _Field& x) const {
			int j = GetDegree();
			_Field p = _vecCoef[j] * _Field(1);

			while (j > 0)
				p = p * x + _vecCoef[--j];
			return p;
		}
		// Given the coefficients of a polynomial of degree nc as an array c[0..nc] of size nc+1 (with
		// c[0] being the constant term), and given a value x, this routine fills an output array pd of size
		// nd+1 with the value of the polynomial evaluated at x in pd[0], and the first nd derivatives at
		// x in pd[1..nd].
		void Derive(const Real x, Vector<Real>& pd)
		{
			int  nnd, j, i;
			int  nc = GetDegree();
			int  nd = pd.size() - 1;
			Real cnst = 1.0;

			pd[0] = (*this)[nc];
			for (j = 1; j < nd + 1; j++)
				pd[j] = 0.0;

			for (i = nc - 1; i >= 0; i--)
			{
				nnd = (nd < (nc - i) ? nd : nc - i);
				for (j = nnd; j > 0; j--)
					pd[j] = pd[j] * x + pd[j - 1];

				pd[0] = pd[0] * x + (*this)[i];
			}
			for (i = 2; i < nd + 1; i++) {
				cnst *= i;
				pd[i] *= cnst;
			}
		}

		bool IsEqual(const Polynom& b) const
		{
			if (_vecCoef.size() != b._vecCoef.size())
				return false;
			for (int i = 0; i < _vecCoef.size(); i++)
				if (_vecCoef[i] != b._vecCoef[i])
					return false;
			return true;
		}

		bool operator==(const Polynom& b) const
		{
			// TODO - nije nuzno - mogu biti isti, a razlicitog reda?
			// ako ne reduciramo red kad su highest coef 0!
			if (_vecCoef.size() != b._vecCoef.size())
				return false;
			for (int i = 0; i < _vecCoef.size(); i++)
				if (_vecCoef[i] != b._vecCoef[i])
					return false;
			return true;
		}

		Polynom operator+(const Polynom& b) const
		{
			Polynom result;
			int n = (int)std::max(_vecCoef.size(), b._vecCoef.size());
			result._vecCoef.resize(n);
			for (int i = 0; i < n; i++)
			{
				if (i < _vecCoef.size())
					result._vecCoef[i] += _vecCoef[i];
				if (i < b._vecCoef.size())
					result._vecCoef[i] += b._vecCoef[i];
			}
			return result;
		}

		Polynom operator-(const Polynom& b) const
		{
			Polynom result;
			int n = (int)std::max(_vecCoef.size(), b._vecCoef.size());
			result._vecCoef.resize(n);
			for (int i = 0; i < n; i++)
			{
				if (i < _vecCoef.size())
					result._vecCoef[i] += _vecCoef[i];
				if (i < b._vecCoef.size())
					result._vecCoef[i] -= b._vecCoef[i];
			}
			return result;
		}

		Polynom operator*(const Polynom& b) const
		{
			Polynom result;

			int n = (int)(_vecCoef.size() + b._vecCoef.size() - 1);
			result._vecCoef.resize(n);
			for (int i = 0; i < _vecCoef.size(); i++)
				for (int j = 0; j < b._vecCoef.size(); j++)
					result._vecCoef[i + j] += _vecCoef[i] * b._vecCoef[j];
			return result;
		}

		static void poldiv(const Polynom& u, const Polynom& v, Polynom& qout, Polynom& rout)
		{
			int k, j, n = u.GetDegree(), nv = v.GetDegree();

			// find real degree of v
			while (nv >= 0 && v._vecCoef[nv] == 0.)
				nv--;

			if (nv < 0)
				throw std::domain_error("poldiv divide by zero polynomial");

			Polynom r = u;
			Polynom q(u.GetDegree());
			for (k = n - nv; k >= 0; k--)
			{
				q[k] = r[nv + k] / v[nv];

				for (j = nv + k - 1; j >= k; j--)
					r[j] -= q[k] * v[j - k];
			}
			for (j = nv; j <= n; j++)
				r[j] = 0.0;


			int nq = q.GetDegree();
			while (nq >= 0 && q[nq] == 0.)
				nq--;

			// setting exact size for quotient
			qout.SetDegree(nq);
			for (j = 0; j <= nq; j++)
				qout[j] = q[j];

			// setting exact size for remainder
			rout.SetDegree(nv - 1);
			for (j = 0; j < nv; j++)
				rout[j] = r[j];
		}

		friend Polynom operator*(const Polynom& a, Real b)
		{
			Polynom ret;
			ret._vecCoef.resize(a.GetDegree() + 1);
			for (int i = 0; i <= a.GetDegree(); i++)
				ret._vecCoef[i] = a._vecCoef[i] * b;
			return ret;
		}

		friend Polynom operator*(Real a, const Polynom& b)
		{
			Polynom ret;
			ret._vecCoef.resize(b.GetDegree() + 1);
			for (int i = 0; i <= b.GetDegree(); i++)
				ret._vecCoef[i] = a * b._vecCoef[i];
			return ret;
		}

		friend Polynom operator/(const Polynom& a, Real b)
		{
			Polynom ret;
			ret._vecCoef.resize(a.GetDegree() + 1);
			for (int i = 0; i <= a.GetDegree(); i++)
				ret._vecCoef[i] = a._vecCoef[i] / b;
			return ret;
		}

		///////////////////////////               I/O                 ///////////////////////////
		std::string to_string(int width, int precision) const
		{
			std::stringstream str;

			Print(str, width, precision);

			return str.str();
		}

		std::ostream& Print(std::ostream& stream, int width, int precision) const
		{
			// first, save current stream state
			std::ios_base::fmtflags f(stream.flags());

			// change formatting
			stream << std::fixed << std::setw(width) << std::setprecision(precision);

			// print the polynomial
			Print(stream);

			// restore stream state
			stream.flags(f);

			return stream;
		}

		std::ostream& Print(std::ostream& stream) const
		{
			for (int i = (int)_vecCoef.size() - 1; i >= 0; i--)
			{
				if (std::abs(_vecCoef[i]) == 0.0)
					continue;

				// handling first term
				if (i == _vecCoef.size() - 1)
				{
					if (i == 0) // means we have only constant term
					{
						stream << _vecCoef[i];
						return stream;
					}

					if (_vecCoef[i] < 0.0)
					{
						if (_vecCoef[i] == -1.0)
							stream << "-x";
						else
							stream << _vecCoef[i] << "*x";;
					}
					else
					{
						if (_vecCoef[i] == 1.0)
							stream << "x";
						else
							stream << _vecCoef[i] << "*x";;
					}
					// adding x exponent
					if (i > 1)
						stream << "^" << i;
				}
				else // handling other terms
				{
					if (i == 0)
					{
						if (_vecCoef[i] > 0.0)
							stream << " + " << _vecCoef[i];
						else
							stream << " - " << std::abs(_vecCoef[i]);
						return stream;
					}

					if (_vecCoef[i] > 0.0)
					{
						if (_vecCoef[i] == 1.0)
							stream << " + x";
						else
							stream << " + " << _vecCoef[i] << "*x";;
					}
					else
					{
						if (_vecCoef[i] == -1.0)
							stream << " - x";
						else
							stream << " - " << std::abs(_vecCoef[i]) << "*x";;
					}

					if (i > 1)
						stream << "^" << i;
				}
			}

			return stream;
		}

		friend std::ostream& operator<<(std::ostream& stream, const Polynom<_Field>& a)
		{
			a.Print(stream);

			return stream;
		}
	};

	class PolynomFunc : public Polynom<Real>, public IRealFunction
	{
		// staviti samo daima const ref na Polynom<Real>
	public:
		PolynomFunc() {}
		PolynomFunc(int n) { _vecCoef.resize(n + 1); }
		PolynomFunc(const std::vector<Real>& vecCoef) : Polynom(vecCoef) {}
		PolynomFunc(std::initializer_list<Real> list) : Polynom(list) {}
		PolynomFunc(const Polynom& Copy) : Polynom(Copy) {}
		~PolynomFunc() {}

		Real operator()(Real x) const { return (*this)(x); }
	};

	typedef Polynom<Real>         RealPolynom;
	typedef Polynom<Complex>   ComplexPolynom;

	typedef Polynom<MatrixNM<Real, 2, 2>>       Matrix2Polynom;
	typedef Polynom<MatrixNM<Real, 3, 3>>       Matrix3Polynom;
	typedef Polynom<MatrixNM<Real, 4, 4>>       Matrix4Polynom;
}

///////////////////////////   ./include/base/VectorTypes.h   ///////////////////////////


namespace MML
{
	class Vector2Cartesian : public VectorN<Real, 2>
	{
	public:
		Vector2Cartesian() {}
		Vector2Cartesian(Real x, Real y)
		{
			_val[0] = x;
			_val[1] = y;
		}
		Vector2Cartesian(const VectorN<Real, 2>& b) : VectorN<Real, 2>{ b[0], b[1] } {}
		Vector2Cartesian(const Point2Cartesian& a, const Point2Cartesian& b)
		{
			_val[0] = b.X() - a.X();
			_val[1] = b.Y() - a.Y();
		}

		Real  X() const { return _val[0]; }
		Real& X() { return _val[0]; }
		Real  Y() const { return _val[1]; }
		Real& Y() { return _val[1]; }
			
		// For Cartesian vector, we will enable operator* to represent standard scalar product
		Real operator*(const Vector2Cartesian& b)
		{
			return this->ScalarProductCartesian(b);
		}

		friend Vector2Cartesian operator*(const Vector2Cartesian& a, Real b)
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = a[i] * b;
			return ret;
		}
		friend Vector2Cartesian operator*(Real a, const Vector2Cartesian& b)
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = a * b[i];
			return ret;
		}
		friend Vector2Cartesian operator/(const Vector2Cartesian& a, Real b)
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = a[i] / b;
			return ret;
		}

		Vector2Cartesian GetUnitVector() const
		{
			VectorN<Real, 2> res = (*this) / NormL2();

			return Vector2Cartesian(res[0], res[1]);
		}
		Vector2Cartesian GetAsUnitVectorAtPos(const Vector2Cartesian& pos) const
		{
			return Vector2Cartesian{ (*this) / NormL2() };
		}

		friend Point2Cartesian operator+(const Point2Cartesian& a, const Vector2Cartesian& b) { return Point2Cartesian(a.X() + b[0], a.Y() + b[1]); }
		friend Point2Cartesian operator-(const Point2Cartesian& a, const Vector2Cartesian& b) { return Point2Cartesian(a.X() - b[0], a.Y() - b[1]); }
	};

	class Vector2Polar : public VectorN<Real, 2>
	{
	public:
		Real  R() const { return _val[0]; }
		Real& R() { return _val[0]; }
		Real  Phi() const { return _val[1]; }
		Real& Phi() { return _val[1]; }

		Vector2Polar() {}
		Vector2Polar(Real r, Real phi)
		{
			_val[0] = r;
			_val[1] = phi;
		}
		Vector2Polar(const VectorN<Real, 2>& b) : VectorN<Real, 2>{ b[0], b[1] } {}

		Vector2Polar GetAsUnitVectorAtPos(const Vector2Polar& pos) const
		{
			// TODO 1.1 - BIG!!!
			return Vector2Polar{ (*this) / NormL2() };
		}
	};

	class Vector3Cartesian : public VectorN<Real, 3>
	{
	public:
		Real  X() const { return _val[0]; }
		Real& X() { return _val[0]; }
		Real  Y() const { return _val[1]; }
		Real& Y() { return _val[1]; }
		Real  Z() const { return _val[2]; }
		Real& Z() { return _val[2]; }

		Vector3Cartesian() : VectorN<Real, 3>{ 0.0, 0.0, 0.0 } {}
		Vector3Cartesian(const VectorN<Real, 3>& b) : VectorN<Real, 3>{ b } {}
		Vector3Cartesian(Real x, Real y, Real z) : VectorN<Real, 3>{ x, y, z } {}
		Vector3Cartesian(std::initializer_list<Real> list) : VectorN<Real, 3>(list) { }
		Vector3Cartesian(const Point3Cartesian& a, const Point3Cartesian& b)
		{
			_val[0] = b.X() - a.X();
			_val[1] = b.Y() - a.Y();
			_val[2] = b.Z() - a.Z();
		}

		// For Cartesian vector, we will enable operator* to represent standard scalar product
		Real operator*(const Vector3Cartesian& b)
		{
			return this->ScalarProductCartesian(b);
		}

		friend Vector3Cartesian operator*(const Vector3Cartesian& a, Real b)
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = a[i] * b;
			return ret;
		}
		friend Vector3Cartesian operator*(Real a, const Vector3Cartesian& b)
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = a * b[i];
			return ret;
		}
		friend Vector3Cartesian operator/(const Vector3Cartesian& a, Real b)
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = a[i] / b;
			return ret;
		}

		friend Vector3Cartesian operator-(const Point3Cartesian& a, const Point3Cartesian& b) { return  Vector3Cartesian{ a.X() - b.X(), a.Y() - b.Y(), a.Z() - b.Z() }; }

		friend Point3Cartesian operator+(const Point3Cartesian& a, const Vector3Cartesian& b) { return Point3Cartesian(a.X() + b[0], a.Y() + b[1], a.Z() + b[2]); }
		friend Point3Cartesian operator-(const Point3Cartesian& a, const Vector3Cartesian& b) { return Point3Cartesian(a.X() - b[0], a.Y() - b[1], a.Z() - b[2]); }

		Point3Cartesian getAsPoint()
		{
			return Point3Cartesian(_val[0], _val[1], _val[2]);
		}

		bool IsParallelTo(const Vector3Cartesian& b, Real eps = 1e-15) const
		{
			Real norm1 = NormL2();
			Real norm2 = b.NormL2();

			return std::abs(X() / norm1 - b.X() / norm2) < eps &&
				std::abs(Y() / norm1 - b.Y() / norm2) < eps &&
				std::abs(Z() / norm1 - b.Z() / norm2) < eps;
		}
		bool IsPerpendicularTo(const Vector3Cartesian& b, Real eps = 1e-15) const
		{
			if (std::abs(ScalarProd(*this, b)) < eps)
				return true;
			else
				return false;
		}
		Real AngleToVector(const Vector3Cartesian& b)
		{
			Real cos_phi = ScalarProd(*this, b) / (NormL2() * b.NormL2());

			return acos(cos_phi);
		}

		Vector3Cartesian GetAsUnitVector() const
		{
			return Vector3Cartesian{ (*this) / NormL2() };
		}
		Vector3Cartesian GetAsUnitVectorAtPos(const Vector3Cartesian& pos) const
		{
			return Vector3Cartesian{ (*this) / NormL2() };
		}

		friend Real ScalarProd(const Vector3Cartesian& a, const Vector3Cartesian& b)
		{
			return a.ScalarProductCartesian(b);
		}
		friend Vector3Cartesian VectorProd(const Vector3Cartesian& a, const Vector3Cartesian& b)
		{
			Vector3Cartesian ret;

			ret.X() = a.Y() * b.Z() - a.Z() * b.Y();
			ret.Y() = a.Z() * b.X() - a.X() * b.Z();
			ret.Z() = a.X() * b.Y() - a.Y() * b.X();

			return ret;
		}
	};

	class Vector3Spherical : public VectorN<Real, 3>
	{
	public:
		Real  R()     const { return _val[0]; }
		Real& R() { return _val[0]; }
		Real  Theta() const { return _val[1]; }
		Real& Theta() { return _val[1]; }
		Real  Phi()   const { return _val[2]; }
		Real& Phi() { return _val[2]; }

		Vector3Spherical() : VectorN<Real, 3>{ 0.0, 0.0, 0.0 } {}
		Vector3Spherical(const VectorN<Real, 3>& b) : VectorN<Real, 3>{ b[0], b[1], b[2] } {}
		Vector3Spherical(Real r, Real theta, Real phi) : VectorN<Real, 3>{ r, theta, phi } {}
		Vector3Spherical(std::initializer_list<Real> list) : VectorN<Real, 3>(list) { }

		// TODO - HIGH, HARD, osnovne operacije +, -
		// sve operaciju pretpostavlju da se odvijaju NA ISTOJ TOCKI U PROSTORU
		Vector3Spherical GetAsUnitVectorAtPos(const Vector3Spherical& pos) const
		{
			// TODO 1.1 - VERIFY this!!!
			return Vector3Spherical{ R(), Theta() / pos.R(), Phi() / (pos.R() * sin(pos.Theta())) };
		}

		std::ostream& PrintDeg(std::ostream& stream, int width, int precision) const
		{
			stream << "[ ";
			stream << std::fixed << std::setw(width) << std::setprecision(precision);
			stream << R();
			stream << ", " << Theta() * 180.0 / Constants::PI;
			stream << ", " << Phi() * 180.0 / Constants::PI << " ]" << std::endl;

			return stream;
		}
	};

	class Vector3Cylindrical : public VectorN<Real, 3>
	{
	public:
		Real    R()   const { return _val[0]; }
		Real& R() { return _val[0]; }
		Real    Phi() const { return _val[1]; }
		Real& Phi() { return _val[1]; }
		Real    Z()   const { return _val[2]; }
		Real& Z() { return _val[2]; }

		Vector3Cylindrical() : VectorN<Real, 3>{ 0.0, 0.0, 0.0 } {}
		Vector3Cylindrical(const VectorN<Real, 3>& b) : VectorN<Real, 3>{ b[0], b[1], b[2] } {}
		Vector3Cylindrical(Real r, Real phi, Real z) : VectorN<Real, 3>{ r, phi, z } {}
		Vector3Cylindrical(std::initializer_list<Real> list) : VectorN<Real, 3>(list) { }

		// TODO - MED, implement Vector3Cylindrical GetAsUniTVector()
		// TODO - MED, razmisliti generalno vektor i vektor at point
		Vector3Cylindrical GetAsUnitVectorAtPos(const Vector3Cylindrical& pos) const
		{
			return Vector3Cylindrical{ R(), Phi() / pos.R(), Z() };
		}
	};

	class Vector4Lorentz : public VectorN<Real, 4>
	{
	public:
		Real  T() const { return _val[0]; }
		Real& T() { return _val[0]; }
		Real  X() const { return _val[1]; }
		Real& X() { return _val[1]; }
		Real  Y() const { return _val[2]; }
		Real& Y() { return _val[2]; }
		Real  Z() const { return _val[3]; }
		Real& Z() { return _val[3]; }

		Vector4Lorentz() : VectorN<Real, 4>{ 0.0, 0.0, 0.0, 0.0 } {}
		Vector4Lorentz(std::initializer_list<Real> list) : VectorN<Real, 4>(list) { }

		Vector4Lorentz GetAsUnitVectorAtPos(const Vector4Lorentz& pos) const
		{
			return *this;
		}
	};

	typedef Vector2Cartesian    Vec2Cart;
	typedef Vector3Cartesian    Vec3Cart;
	typedef Vector3Spherical    Vec3Sph;
	typedef Vector3Cylindrical  Vec3Cyl;
}
///////////////////////////   ./include/base/Geometry2D.h   ///////////////////////////


namespace MML
{
	class Line2D
	{
	private:
		Point2Cartesian _point;
		Vector2Cartesian _direction; // unit vector in line direction

	public:
		Line2D(const Point2Cartesian& pnt, const Vector2Cartesian dir)
		{
			_point = pnt;
			_direction = dir.GetUnitVector();
		}

		Line2D(const Point2Cartesian& a, const Point2Cartesian& b)
		{
			Vector2Cartesian dir(a, b);
			_point = a;
			_direction = dir.GetUnitVector();
		}

		Point2Cartesian   StartPoint() const { return _point; }
		Point2Cartesian&  StartPoint() { return _point; }

		Vector2Cartesian  Direction() const { return _direction; }
		Vector2Cartesian& Direction() { return _direction; }

		Point2Cartesian PointOnLine(Real t)
		{
			Vector2Cartesian dist = t * _direction;
			Point2Cartesian ret = _point + dist;
			return ret;
		}
	};

	class SegmentLine2D
	{
	private:
		Point2Cartesian _point1;
		Point2Cartesian _point2;

	public:
		SegmentLine2D(Point2Cartesian pnt1, Point2Cartesian pnt2) : _point1(pnt1), _point2(pnt2)
		{ }

		SegmentLine2D(const Point2Cartesian& pnt1, const Vector2Cartesian& direction, Real t) : _point1(pnt1)
		{
			_point2 = pnt1 + direction * t;
		}

		Point2Cartesian  StartPoint() const { return _point1; }
		Point2Cartesian& StartPoint() { return _point1; }

		Point2Cartesian  EndPoint()  const { return _point2; }
		Point2Cartesian& EndPoint() { return _point2; }

		Point2Cartesian PointOnSegment(Real t)
		{
			if (t < 0.0 || t > 1.0)
				throw std::invalid_argument("SegmentLine2D::PointOnSegment t must be in [0,1]");

			Vector2Cartesian dist = t * Direction();
			Point2Cartesian ret = _point1 + dist;
			return ret;
		}

		Real                Length()    const { return _point1.Dist(_point2); }
		Vector2Cartesian    Direction() const { return Vector2Cartesian(_point1, _point2); }
	};

	class Triangle2D
	{
	private:
		Point2Cartesian _pnt1, _pnt2, _pnt3;

	public:
		Triangle2D(Point2Cartesian pnt1, Point2Cartesian pnt2, Point2Cartesian pnt3) : _pnt1(pnt1), _pnt2(pnt2), _pnt3(pnt3)
		{ }

		Point2Cartesian  Pnt1() const { return _pnt1; }
		Point2Cartesian& Pnt1() { return _pnt1; }
		Point2Cartesian  Pnt2() const { return _pnt2; }
		Point2Cartesian& Pnt2() { return _pnt2; }
		Point2Cartesian  Pnt3() const { return _pnt3; }
		Point2Cartesian& Pnt3() { return _pnt3; }

		Real Area() const
		{
			Real a = _pnt1.Dist(_pnt2);
			Real b = _pnt2.Dist(_pnt3);
			Real c = _pnt3.Dist(_pnt1);

			Real s = (a + b + c) / 2.0;

			return sqrt(s * (s - a) * (s - b) * (s - c));
		}
	};

	class Polygon2D
	{
	private:
		std::vector<Point2Cartesian> _points;
	public:
		Polygon2D() {}
		Polygon2D(std::vector<Point2Cartesian> points) : _points(points) {}
		Polygon2D(std::initializer_list<Point2Cartesian> list)
		{
			for (auto element : list)
				_points.push_back(element);
		}

		std::vector<Point2Cartesian>  Points() const { return _points; }
		std::vector<Point2Cartesian>& Points() { return _points; }

		Real Area() const
		{
			Real area = 0.0;
			int n = (int)_points.size();
			for (int i = 0; i < n; i++)
			{
				area += _points[i].X() * _points[(i + 1) % n].Y();
				area -= _points[i].Y() * _points[(i + 1) % n].X();
			}
			area /= 2.0;
			return area;
		}

		bool IsSimple() const
		{
			return false;
		}

		bool IsConvex() const
		{
			int n = (int)_points.size();
			if (n < 3)
				return false;
			bool sign = false;
			for (int i = 0; i < n; i++)
			{
				Vector2Cartesian v1(_points[(i + 1) % n], _points[i]);
				Vector2Cartesian v2(_points[(i + 2) % n], _points[(i + 1) % n]);
				Real cross = v1.X() * v2.Y() - v1.Y() * v2.X();
				if (i == 0)
					sign = cross > 0;
				else if (sign != (cross > 0))
					return false;
			}
			return true;
		}

		std::vector<Triangle2D> Triangularization() const
		{
			std::vector<Triangle2D> triangles;
			int n = (int)_points.size();
			if (n < 3)
				return triangles;
			std::vector<int> indices(n);
			for (int i = 0; i < n; i++)
				indices[i] = i;
			int i = 0;
			while (n > 3)
			{
				int i1 = indices[i];
				int i2 = indices[(i + 1) % n];
				int i3 = indices[(i + 2) % n];
				Triangle2D tri(_points[i1], _points[i2], _points[i3]);
				triangles.push_back(tri);
				indices.erase(indices.begin() + (i + 1) % n);
				n--;
				i = (i + 1) % n;
			}
			Triangle2D tri(_points[indices[0]], _points[indices[1]], _points[indices[2]]);
			triangles.push_back(tri);
			return triangles;
		}

		bool IsInside(Point2Cartesian pnt) const
		{
			int n = (int)_points.size();
			if (n < 3)
				return false;
		
			return false;
		}
	};
}
///////////////////////////   ./include/base/Geometry3D.h   ///////////////////////////



namespace MML
{
	class Line3D
	{
	private:
		Point3Cartesian  _point;
		Vector3Cartesian _direction;

	public:
		Line3D() {}
		// by default, direction vector is normalized to unit vector (but, it need not be such!)
		Line3D(const Point3Cartesian& pnt, const Vector3Cartesian dir)
		{
			// check for null vector as direction
			if (dir.X() == 0.0 && dir.Y() == 0.0 && dir.Z() == 0.0)
				throw std::runtime_error("Line3D ctor - null vector as direction");

			_point = pnt;
			_direction = dir.GetAsUnitVector();
		}

		Line3D(const Point3Cartesian& a, const Point3Cartesian& b)
		{
			// check for same points
			if (a == b)
				throw std::runtime_error("Line3D ctor - same points");

			Vector3Cartesian dir(a, b);
			_point = a;
			_direction = dir.GetAsUnitVector();
		}

		Point3Cartesian   StartPoint() const { return _point; }
		Point3Cartesian&  StartPoint() { return _point; }

		Vector3Cartesian  Direction() const { return _direction; }
		Vector3Cartesian& Direction() { return _direction; }

		Point3Cartesian PointOnLine(Real t) const { return _point + t * _direction; }

		bool IsPerpendicular(const Line3D& b) const
		{
			return ScalarProd(Direction(), b.Direction()) == 0.0f;
		}
		bool IsParallel(const Line3D& b) const
		{
			return Direction() == b.Direction();
		}

		Real Dist(const Point3Cartesian& pnt) const
		{
			// Bronshtein 3.394
			const Real a = pnt.X();
			const Real b = pnt.Y();
			const Real c = pnt.Z();

			const Real x1 = StartPoint().X();
			const Real y1 = StartPoint().Y();
			const Real z1 = StartPoint().Z();

			const Real l = Direction().X();
			const Real m = Direction().Y();
			const Real n = Direction().Z();

			Real numer = POW2((a - x1) * m - (b - y1) * l) + POW2((b - y1) * n - (c - z1) * m) + POW2((c - z1) * l - (a - x1) * n);
			Real denom = l * l + m * m + n * n;

			return sqrt(numer / denom);
		}

		Real Dist(const Line3D& line) const
		{
			// https://math.stackexchange.com/questions/2213165/distance-between-two-lines-in-3d-space
			// https://en.wikipedia.org/wiki/Skew_lines#Nearest_points
			Point3Cartesian  p1 = StartPoint();
			Vector3Cartesian d1 = Direction();
			Point3Cartesian  p2 = line.StartPoint();
			Vector3Cartesian d2 = line.Direction();

			Vector3Cartesian n = VectorProd(d1, d2);

			if (n.IsNullVec())
			{
				// parallel lines
				return Vector3Cartesian(p1, p2) * d1 / d2.NormL2();
			}
			else
			{
				// skew lines
				return Abs(n * Vector3Cartesian(p1, p2)) / n.NormL2();
			}
		}

		Real Dist(const Line3D& line, Point3Cartesian& out_line1_pnt, Point3Cartesian& out_line2_pnt) const
		{
			// https://math.stackexchange.com/questions/2213165/distance-between-two-lines-in-3d-space
			// https://en.wikipedia.org/wiki/Skew_lines#Nearest_points
			Point3Cartesian  p1 = StartPoint();
			Vector3Cartesian d1 = Direction();
			Point3Cartesian  p2 = line.StartPoint();
			Vector3Cartesian d2 = line.Direction();
			Real dist = 0.0;

			Vector3Cartesian n = VectorProd(d1, d2);

			if (n.IsNullVec())
			{
				// parallel lines
				dist = Vector3Cartesian(p1, p2) * d1 / d2.NormL2();
				// should we throw exception???
				// because we can't get out_line_pnt's
			}
			else
			{
				// skew lines
				dist = Abs(n * Vector3Cartesian(p1, p2)) / n.NormL2();

				Real t1 = (VectorProd(d2, n) * Vector3Cartesian(p1, p2)) / POW2(n.NormL2());
				Real t2 = (VectorProd(d1, n) * Vector3Cartesian(p1, p2)) / POW2(n.NormL2());

				out_line1_pnt = p1 + t1 * d1;
				out_line2_pnt = p2 + t2 * d2;
			}

			return dist;
		}

		Point3Cartesian NearestPointOnLine(const Point3Cartesian& pnt) const
		{
			// https://math.stackexchange.com/questions/1521128/given-a-line-and-a-point-in-3d-how-to-find-the-closest-point-on-the-line         
			Vector3Cartesian line_dir = this->Direction();
			Vector3Cartesian rad_vec_AP(StartPoint(), pnt);

			double t = rad_vec_AP * line_dir / POW2(line_dir.NormL2());

			return StartPoint() + t * line_dir;
		}

		// TODO - Intersection of two lines
		bool Intersection(const Line3D& line, Point3Cartesian& out_inter_pnt) const
		{
			// https://en.wikipedia.org/wiki/Skew_lines#Nearest_points
			// https://math.stackexchange.com/questions/2213165/distance-between-two-lines-in-3d-space

			return true;
		}

		// TODO - pravac koji prolazi kroz tocku i sijece zadani pravac okomito
		Line3D PerpendicularLineThroughPoint(const Point3Cartesian& pnt)
		{
			Line3D ret;

			return ret;
		}
	};

	class SegmentLine3D
	{
	private:
		Point3Cartesian _point1;
		Point3Cartesian _point2;

	public:
		SegmentLine3D(Point3Cartesian pnt1, Point3Cartesian pnt2) : _point1(pnt1), _point2(pnt2)
		{
		}

		SegmentLine3D(Point3Cartesian pnt1, Vector3Cartesian direction, Real t)
		{
			_point1 = pnt1;
			_point2 = pnt1 + t * direction;
		}

		Point3Cartesian   StartPoint() const { return _point1; }
		Point3Cartesian& StartPoint() { return _point1; }

		Point3Cartesian   EndPoint() const { return _point2; }
		Point3Cartesian& EndPoint() { return _point2; }

		Point3Cartesian PointOnSegment(Real t)
		{
			if (t < 0.0 || t > 1.0)
				throw std::runtime_error("SegmentLine3D::PointOnSegment - t not in [0,1]");

			return _point1 + t * Direction();
		}

		Real                Length()    const { return _point1.Dist(_point2); }
		Vector3Cartesian    Direction() const { return Vector3Cartesian(_point1, _point2); }
	};

	class Plane3D
	{
	private:
		Real _A, _B, _C, _D;

	public:
		Plane3D(const Point3Cartesian& a, const Vector3Cartesian& normal)
		{
			if (normal.IsNullVec())
				throw std::runtime_error("Plane3D ctor - normal is null vector");

			Vector3Cartesian unitNormal = normal.GetAsUnitVector();

			_A = unitNormal.X();
			_B = unitNormal.Y();
			_C = unitNormal.Z();
			_D = -(a.X() * unitNormal.X() + a.Y() * unitNormal.Y() + a.Z() * unitNormal.Z());
		}

		Plane3D(const Point3Cartesian& a, const Point3Cartesian& b, const Point3Cartesian& c)
			: Plane3D(a, VectorProd(Vector3Cartesian(a, b), Vector3Cartesian(a, c)))
		{
		}

		Plane3D(Real alpha, Real beta, Real gamma, Real d)      // Hesseov (normalni) oblik
		{
			_A = cos(alpha);
			_B = cos(beta);
			_C = cos(gamma);
			_D = -d;
		}

		// tri segmenta na koord osima ctor
		Plane3D(Real seg_x, Real seg_y, Real seg_z)
		{
			if (seg_x == 0 || seg_y == 0 || seg_z == 0)
				throw std::runtime_error("Plane3D ctor - zero segment");

			_A = 1 / seg_x;
			_B = 1 / seg_y;
			_C = 1 / seg_z;
			_D = -1;
		}

		static Plane3D GetXYPlane() { return Plane3D(Point3Cartesian(0, 0, 0), Vector3Cartesian(0, 0, 1)); }
		static Plane3D GetXZPlane() { return Plane3D(Point3Cartesian(0, 0, 0), Vector3Cartesian(0, 1, 0)); }
		static Plane3D GetYZPlane() { return Plane3D(Point3Cartesian(0, 0, 0), Vector3Cartesian(1, 0, 0)); }

		Real  A() const { return _A; }
		Real& A() { return _A; }
		Real  B() const { return _B; }
		Real& B() { return _B; }
		Real  C() const { return _C; }
		Real& C() { return _C; }
		Real  D() const { return _D; }
		Real& D() { return _D; }

		Point3Cartesian GetPointOnPlane() const {
			if (_A != 0.0)
				return Point3Cartesian(-_D / _A, 0, 0);
			else  if (_B != 0.0)
				return Point3Cartesian(0, -_D / _B, 0);
			else
				return Point3Cartesian(0, 0, -_D / _C);
		}

		Vector3Cartesian Normal() const { return Vector3Cartesian(_A, _B, _C); }

		void GetCoordAxisSegments(Real& outseg_x, Real& outseg_y, Real& outseg_z)
		{
			if (_A != 0.0)
				outseg_x = -_D / _A;
			else
				outseg_x = Constants::PositiveInf;

			if (_B != 0.0)
				outseg_y = -_D / _B;
			else
				outseg_y = Constants::PositiveInf;

			if (_C != 0.0)
				outseg_z = -_D / _C;
			else
				outseg_z = Constants::PositiveInf;
		}

		bool IsPointOnPlane(const Point3Cartesian& pnt, Real defEps = 1e-15) const
		{
			return std::abs(_A * pnt.X() + _B * pnt.Y() + _C * pnt.Z() + _D) < defEps;
		}
		Real DistToPoint(const Point3Cartesian& pnt) const
		{
			Real a = _A * pnt.X() + _B * pnt.Y() + _C * pnt.Z() + _D;
			Real b = sqrt(_A * _A + _B * _B + _C * _C);

			return std::abs(a / b);
		}

		Point3Cartesian ProjectionToPlane(const Point3Cartesian& pnt) const
		{
			if (IsPointOnPlane(pnt))
				return pnt;

			Line3D line(pnt, Normal());

			// find intersection of this line with plane
			Point3Cartesian ret;
			IntersectionWithLine(line, ret);

			return ret;
		}

		bool IsLineOnPlane(const Line3D& line) const
		{
			// get two points
			const Point3Cartesian pnt1 = line.StartPoint();
			const Point3Cartesian pnt2 = line.PointOnLine(1.0);

			if (IsPointOnPlane(pnt1) && IsPointOnPlane(pnt2))
				return true;
			else
				return false;
		}

		Real AngleToLine(const Line3D& line) const
		{
			// angle between line and normal to plane
			return Constants::PI / 2.0 - line.Direction().AngleToVector(this->Normal());
		}

		bool IntersectionWithLine(const Line3D& line, Point3Cartesian& out_inter_pnt) const
		{
			// Calculate the direction vector of the line
			Vector3Cartesian line_dir = line.Direction();

			// Calculate the normal vector of the plane
			Vector3Cartesian plane_normal = Normal();

			// Check if the line is parallel to the plane
			if (line_dir.ScalarProductCartesian(plane_normal) == 0) {
				return false;
			}

			// Calculate the distance between the point on the line and the plane
			double dist = Vector3Cartesian(GetPointOnPlane(), line.StartPoint()).ScalarProductCartesian(plane_normal) / line_dir.ScalarProductCartesian(plane_normal);

			// Calculate the point of intersection
			Point3Cartesian inter_pnt = line.StartPoint() + line_dir * dist;

			// Set the output parameter to the point of intersection
			out_inter_pnt = inter_pnt;

			return true;
		}

		bool IsParallelToPlane(const Plane3D& plane) const
		{
			Vector3Cartesian norm1(_A, _B, _C);

			Vector3Cartesian norm2(plane._A, plane._B, plane._C);

			return norm1.IsParallelTo(norm2);
		}
		bool IsPerpendicularToPlane(const Plane3D& plane) const
		{
			Vector3Cartesian norm1(_A, _B, _C);
			Vector3Cartesian norm2(plane._A, plane._B, plane._C);

			return norm1.IsPerpendicularTo(norm2);
		}
		Real AngleToPlane(const Plane3D& plane) const
		{
			// to je kut normala
			return this->Normal().AngleToVector(plane.Normal());
		}
		Real DistToPlane(const Plane3D& plane) const
		{
			// TODO finish
			// ili su paralelne, pa imamo neki broj, ili se sijeku pa je 0
			return 0.0;
		}
		// TODO - check implementation IntersectionWithPlane
		bool IntersectionWithPlane(const Plane3D& plane, Line3D& out_inter_line) const
		{
			Vector3Cartesian inter_dir = VectorProd(Normal(), plane.Normal());

			// Check if the planes are parallel
			if (inter_dir.NormL2() == 0.0) {
				return false;
			}

			// Calculate a point on the intersection line
			Point3Cartesian inter_pnt = GetPointOnPlane();

			// Calculate the distance between the intersection line and the two planes
			double dist1 = Vector3Cartesian(inter_pnt, GetPointOnPlane()).ScalarProductCartesian(plane.Normal());
			double dist2 = Vector3Cartesian(inter_pnt, plane.GetPointOnPlane()).ScalarProductCartesian(Normal());

			// Calculate the point of intersection
			inter_pnt = inter_pnt - inter_dir * (dist1 / inter_dir.ScalarProductCartesian(plane.Normal()));

			// Set the output parameter to the intersection line
			out_inter_line = Line3D(inter_pnt, inter_dir);

			return true;
			return false;
		}
	};

	class Triangle3D
	{
	private:
		Point3Cartesian _pnt1, _pnt2, _pnt3;
	public:
		Triangle3D(Point3Cartesian pnt1, Point3Cartesian pnt2, Point3Cartesian pnt3)
			: _pnt1(pnt1), _pnt2(pnt2), _pnt3(pnt3)
		{
		}

		Point3Cartesian& Pnt1() { return _pnt1; }
		Point3Cartesian& Pnt2() { return _pnt2; }
		Point3Cartesian& Pnt3() { return _pnt3; }
	};

	class TriangleSurface3D : public Triangle3D, IParametricSurface<3>
	{
	public:
		Real _minX, _maxX, _minY, _maxY;
		Point3Cartesian _center;
		Vector3Cartesian _localX, _localY;

		TriangleSurface3D(Point3Cartesian pnt1, Point3Cartesian pnt2, Point3Cartesian pnt3)
			: Triangle3D(pnt1, pnt2, pnt3)
		{
			// calculate min and max
			_minX = std::min(pnt1.X(), std::min(pnt2.X(), pnt3.X()));
			_maxX = std::max(pnt1.X(), std::max(pnt2.X(), pnt3.X()));
			_minY = std::min(pnt1.Y(), std::min(pnt2.Y(), pnt3.Y()));
			_maxY = std::max(pnt1.Y(), std::max(pnt2.Y(), pnt3.Y()));

			// calculate center
			_center = (pnt1 + pnt2 + pnt3) / 3.0;

			// calculate local coordinate system
			_localX = Vector3Cartesian(pnt1, pnt2).GetAsUnitVector();
			_localY = Vector3Cartesian(pnt1, pnt3).GetAsUnitVector();
		}

		virtual Real getMinX() const { return _minX; }
		virtual Real getMaxX() const { return _maxX; }
		virtual Real getMinY() const { return _minY; }
		virtual Real getMaxY() const { return _maxY; }
	};

	class RectSurface3D : public IParametricSurface<3>
	{
	public:
		// predstavlja PRAVOKUTNU povrsinu u 3D
		Point3Cartesian _pnt1, _pnt2, _pnt3, _pnt4;
		Real _minX, _maxX, _minY, _maxY;
		Point3Cartesian _center;
		Vector3Cartesian _localX, _localY;

		RectSurface3D() {}
		// zadati i centralnom tockom, vektorom normale, uz dodatni a i b!
		RectSurface3D(Point3Cartesian pnt1, Point3Cartesian pnt2, Point3Cartesian pnt3, Point3Cartesian pnt4)
			: _pnt1(pnt1), _pnt2(pnt2), _pnt3(pnt3), _pnt4(pnt4)
		{
			// KLJUCNO - provjeriti da li su sve u ravnini!

			// podesiti min i max
			Real lenX = _pnt1.Dist(_pnt2);
			Real lenY = _pnt1.Dist(_pnt4);
			_minX = -lenX / 2;
			_maxX = lenX / 2;
			_minY = -lenY / 2;
			_maxY = lenY / 2;

			// calculate center
			_center = (_pnt1 + _pnt2 + _pnt3 + _pnt4) / 4.0;

			// calculate local coordinate system
			_localX = Vector3Cartesian(_pnt1, _pnt2).GetAsUnitVector();
			_localY = Vector3Cartesian(_pnt1, _pnt4).GetAsUnitVector();
		}
		virtual Real getMinU() const { return _minX; }
		virtual Real getMaxU() const { return _maxX; }
		virtual Real getMinW() const { return _minY; }
		virtual Real getMaxW() const { return _maxY; }

		Vector3Cartesian getNormal() const {
			return VectorProd(Vector3Cartesian(_pnt1, _pnt2), Vector3Cartesian(_pnt1, _pnt4)).GetAsUnitVector();
		}
		Point3Cartesian getCenter() const { return _center; }

		Real getArea() const {
			return VectorProd(Vector3Cartesian(_pnt1, _pnt2), Vector3Cartesian(_pnt1, _pnt4)).NormL2();
		}

		// vraca dva Triangle3D - i orijentacija je parametar! (kako ce odabrati tocke)

		VectorN<Real, 3> operator()(Real u, Real w) const {
			Point3Cartesian ret = _center + u * _localX + w * _localY;
			return VectorN<Real, 3>({ ret.X(), ret.Y(), ret.Z() });
		}
	};

	class IBody
	{
	public:
		virtual bool isInside(const Point3Cartesian& pnt) const = 0;

		// must also have boundary defined

	};

	class BodyWithRectSurfaces : public IBody
	{

	};

	class BodyWithTriangleSurfaces : public IBody
	{

	};

	class BodyWithBoundary : public IBody
	{

	};

	// TODO - IntegrableSolid? koji ima i potrebne funkcije kojima definira granice tijela?
	class IntegrableVolume3D
	{
		// osigurava da se znaju funkcije koje definiraju granice tijela
		// da se moze obaviti volume integracija
	};

	class SolidSurfaces3D : public BodyWithRectSurfaces
	{
	public:
		// solid body in 3D defined by surfaces
		std::vector<RectSurface3D> _surfaces;

	public:
		bool isInside(const Point3Cartesian& pnt) const
		{
			// TODO - implement
			return false;
		}
		// isClosed() - iz centra mase (?) odasilje zrake u svim smje
		// rovima i gleda da li je pogodio iti jednu povrsinu
		// vraca listu povrsina, koje bi TREBALE omedjivati tijelo!

		// i po tome se moze obaviti surface integracija!!!
	};

	// represents solid that is composed of other solids
	// dobar nacin za provjeriti je li composed solid "tight" je izracunati fluks kroz sve povrsine
	class ComposedSolidSurfaces3D
	{
		bool IsInside(const Point3Cartesian& pnt) const
		{
			// TODO - implement
			return false;
		}
	};

	class Cube3D : public SolidSurfaces3D
	{
		// kocka
		Real _a;
		Vector3Cartesian _center;
	public:
		Cube3D(Real a) : _a(a)
		{
			Point3Cartesian pnt1(a / 2, -a / 2, -a / 2);
			Point3Cartesian pnt2(a / 2, a / 2, -a / 2);
			Point3Cartesian pnt3(-a / 2, a / 2, -a / 2);
			Point3Cartesian pnt4(-a / 2, -a / 2, -a / 2);
			Point3Cartesian pnt5(a / 2, -a / 2, a / 2);
			Point3Cartesian pnt6(a / 2, a / 2, a / 2);
			Point3Cartesian pnt7(-a / 2, a / 2, a / 2);
			Point3Cartesian pnt8(-a / 2, -a / 2, a / 2);

			// dodati svih 6 stranica u popis povrsina
			_surfaces.push_back(RectSurface3D(pnt1, pnt4, pnt3, pnt2));     // lower side in xy plane
			_surfaces.push_back(RectSurface3D(pnt5, pnt6, pnt7, pnt8));     // upper side in xy plane
			_surfaces.push_back(RectSurface3D(pnt1, pnt2, pnt6, pnt5));     // front side in yz plane
			_surfaces.push_back(RectSurface3D(pnt4, pnt8, pnt7, pnt3));     // back side in yz plane
			_surfaces.push_back(RectSurface3D(pnt1, pnt5, pnt8, pnt4));     // left side in xz plane
			_surfaces.push_back(RectSurface3D(pnt2, pnt3, pnt7, pnt6));     // right side in xz plane
		}
		Cube3D(Real a, const Vector3Cartesian& center) : Cube3D(a)
		{
			_center = center;
		}
	};
}

///////////////////////////   ./include/base/BaseUtils.h   ///////////////////////////


namespace MML
{
	namespace Utils
	{
		static int LeviCivita(int i, int j, int k)
		{
			if (i == j || j == k || i == k)
				return 0;

			if (i == 1 && j == 2 && k == 3)
				return 1;
			if (i == 2 && j == 3 && k == 1)
				return 1;
			if (i == 3 && j == 1 && k == 2)
				return 1;
			if (i == 3 && j == 2 && k == 1)
				return -1;
			if (i == 2 && j == 1 && k == 3)
				return -1;
			if (i == 1 && j == 3 && k == 2)
				return -1;

			return 0;
		}
		static int LeviCivita(int i, int j, int k, int l)
		{
			int a[4] = { i, j, k, l };
			int ret = 1;

			// TODO check if 1, 2, 3, 4 are present
			for (int i = 0; i < 4; i++)
				for (int j = i + 1; j < 4; j++)
					if (a[i] > a[j])
					{
						int tmp = a[i];
						a[i] = a[j];
						a[j] = tmp;
						ret *= -1;
					}
			return ret;
		}

		static Real DegToRad(Real angleDeg) { return angleDeg * Constants::PI / 180.0; }
		static Real RadToDeg(Real angleRad) { return angleRad * 180.0 / Constants::PI; }

		static void AngleDegToExplicit(Real angle, Real& deg, Real& min, Real& sec)
		{
			deg = floor(angle);
			min = floor((angle - deg) * 60.0);
			sec = (angle - deg - min / 60.0) * 3600.0;
		}
		static void AngleRadToExplicit(Real angleRad, Real& deg, Real& min, Real& sec) { AngleDegToExplicit(angleRad * 180.0 / Constants::PI, deg, min, sec); }
		static Real ExplicitToAngleDeg(Real deg, Real min, Real sec) { return deg + min / 60.0 + sec / 3600.0; }
		static Real ExplicitToAngleRad(Real deg, Real min, Real sec) { return ExplicitToAngleDeg(deg, min, sec) * Constants::PI / 180.0; }

		///////////////////                     Complex helpers                   ///////////////////
		static bool AreEqual(const Complex& a, const Complex& b, double eps = Defaults::ComplexEqualityPrecision)
		{
			if (std::abs(a.real() - b.real()) > eps || std::abs(a.imag() - b.imag()) > eps)
				return false;
			return true;
		}
		static bool AreEqualAbs(const Complex& a, const Complex& b, double eps = Defaults::ComplexAbsEqualityPrecision)
		{
			if (Abs(a - b) > eps)
				return false;
			return true;
		}
		static bool AreEqual(const Vector<Complex>& a, const Vector<Complex>& b, double eps = Defaults::ComplexEqualityPrecision)
		{
			if (a.size() != b.size())
				return false;

			for (int i = 0; i < a.size(); i++)
				if( !AreEqual(a[i], b[i], eps) )
					return false;

			return true;
		}
		static bool AreEqualAbs(const Vector<Complex>& a, const Vector<Complex>& b, double eps = Defaults::ComplexAbsEqualityPrecision)
		{
			if (a.size() != b.size())
				return false;

			for (int i = 0; i < a.size(); i++)
				if ( !AreEqualAbs(a[i], b[i], eps) )
					return false;

			return true;
		}

		///////////////////                     Vector helpers                    ///////////////////
		static Vector<Real> VectorProjectionParallelTo(const Vector<Real>& orig, const Vector<Real>& b)
		{
			return orig.ScalarProductCartesian(b) / b.NormL2() * b;
		}
		static Vector<Real> VectorProjectionPerpendicularTo(const Vector<Real>& orig, const Vector<Real>& b)
		{
			return orig - VectorProjectionParallelTo(orig, b);
		}

		static Real ScalarProduct(const Vector<Real>& a, const Vector<Real>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("ScalarProduct - must be same dim", a.size(), b.size());

			Real ret = 0;
			for (int i = 0; i < a.size(); i++)
				ret += a[i] * b[i];

			return ret;
		}
		static Complex ScalarProduct(const Vector<Complex>& a, const Vector<Complex>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("ScalarProduct - must be same dim", a.size(), b.size());

			Complex ret = 0;
			for (int i = 0; i < a.size(); i++)
				ret += a[i] * std::conj(b[i]);

			return ret;
		}

		template<class Type> static Matrix<Type> OuterProduct(const Vector<Type>& a, const Vector<Type>& b)
		{
			Matrix<Type> ret(a.size(), b.size());

			for (int i = 0; i < a.size(); i++)
				for (int j = 0; j < b.size(); j++)
					ret[i][j] = a[i] * b[j];

			return ret;
		}

		///////////////////       Vector<Complex> - Vector<Real> operations       ///////////////////
		static Vector<Complex> AddVec(const Vector<Complex>& a, const Vector<Real>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("AddVec(Complex, Real) - must be same dim", a.size(), b.size());

			Vector<Complex> ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret[i] = a[i] + b[i];
			return ret;
		}
		static Vector<Complex> AddVec(const Vector<Real>& a, const Vector<Complex>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("AddVec(Real, Complex) - must be same dim", a.size(), b.size());

			Vector<Complex> ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret[i] = a[i] + b[i];
			return ret;
		}

		static Vector<Complex> SubVec(const Vector<Complex>& a, const Vector<Real>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("SubVec(Complex, Real) - must be same dim", a.size(), b.size());

			Vector<Complex> ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret[i] = a[i] - b[i];
			return ret;
		}
		static Vector<Complex> SubVec(const Vector<Real>& a, const Vector<Complex>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("SubVec(Real, Complex) - must be same dim", a.size(), b.size());

			Vector<Complex> ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret[i] = a[i] - b[i];
			return ret;
		}
	};

	namespace MatrixUtils
	{
		const static inline MatrixNM<Complex, 2, 2> Pauli[] = {
				MatrixNM<Complex, 2, 2>{ 0, 1,
																 1, 0 },
				MatrixNM<Complex, 2, 2>{ 0,              Complex(0, -1),
																 Complex(0, 1),  0},
				MatrixNM<Complex, 2, 2>{ 1,  0,
																 0, -1 }
		};
		const static inline MatrixNM<Complex, 4, 4> DiracGamma[] = {
				MatrixNM<Complex, 4, 4>{ 1,  0,  0,  0,
																 0,  1,  0,  0,
																 0,  0, -1,  0,
																 0,  0,  0, -1 },
				MatrixNM<Complex, 4, 4>{ 0,  0,  0,  1,
																 0,  0,  1,  0,
																 0, -1,  0,  0,
																-1,  0,  0,  0 },
				MatrixNM<Complex, 4, 4>{              0,             0,             0, Complex(0, -1),
																              0,             0, Complex(0, 1),              0,
																              0, Complex(0, 1),             0,              0,
																 Complex(0, -1),             0,             0,              0 },
				MatrixNM<Complex, 4, 4>{ 0,  0,  1,  0,
																 0,  0,  0, -1,
																-1,  0,  0,  0,
																 0,  1,  0,  0 }
		};
		const static inline MatrixNM<Complex, 4, 4> DiracGamma5{ 0,  0,  1,  0,
																														 0,  0,  0,  1,
																														 1,  0,  0,  0,
																														 0,  1,  0,  0 };

		///////////////////             Creating Matrix from Vector              ///////////////////
		template<class Type> static Matrix<Type> RowMatrixFromVector(const Vector<Type>& b)
		{
			Matrix<Type> ret(1, (int)b.size());
			for (int i = 0; i < b.size(); i++)
				ret[0][i] = b[i];

			return ret;
		}
		template<class Type> static Matrix<Type> ColumnMatrixFromVector(const Vector<Type>& b)
		{
			Matrix<Type> ret((int)b.size(), 1);
			for (int i = 0; i < b.size(); i++)
				ret[i][0] = b[i];

			return ret;
		}
		template<class Type> static Matrix<Type> DiagonalMatrixFromVector(const Vector<Type>& b)
		{
			Matrix<Type> ret((int)b.size(), (int)b.size());
			for (int i = 0; i < b.size(); i++)
				ret[i][i] = b[i];

			return ret;
		}

		///////////////////                   Matrix helpers                     ///////////////////
		template<class Type> static Matrix<Type> Commutator(const Matrix<Type>& a, const Matrix<Type>& b)
		{
			return a * b - b * a;
		}
		template<class Type> static Matrix<Type> AntiCommutator(const Matrix<Type>& a, const Matrix<Type>& b)
		{
			return a * b + b * a;
		}

		template<class Type> static void MatrixDecomposeToSymAntisym(const Matrix<Type>& orig, Matrix<Type>& outSym, Matrix<Type>& outAntiSym)
		{
			if (orig.RowNum() != orig.ColNum())
				throw MatrixDimensionError("MatrixDecompose - matrix must be square", orig.RowNum(), orig.ColNum(), -1, -1);

			outSym = (orig + orig.GetTranspose()) * 0.5;
			outAntiSym = (orig - orig.GetTranspose()) * 0.5;
		}

		///////////////////                  Matrix functions                    ///////////////////
		template<class Type> static Matrix<Type> Exp(const Matrix<Type>& a, int n = 10)
		{
			Matrix<Type> ret(a.RowNum(), a.ColNum());
			Matrix<Type> a_pow_n(a);

			double fact = 1.0;
			ret.MakeUnitMatrix();

			for (int i = 1; i <= n; i++)
			{
				ret += a_pow_n / fact;

				a_pow_n = a_pow_n * a;
				fact *= (double)i;
			}

			return ret;
		}

		///////////////////                Real matrix helpers                   ///////////////////
		static bool IsOrthogonal(const Matrix<Real>& mat, double eps = Defaults::IsMatrixOrthogonalPrecision)
		{
			if (mat.RowNum() != mat.ColNum())
				throw MatrixDimensionError("IsOrthogonal - matrix must be square", mat.RowNum(), mat.ColNum(), -1, -1);

			Matrix<Real> matProd = mat * mat.GetTranspose();

			return matProd.IsUnit(eps);
		}

		///////////////////               Complex matrix helpers                 ///////////////////
		static Matrix<Real> GetRealPart(const Matrix<Complex>& a)
		{
			Matrix<Real> ret(a.RowNum(), a.ColNum());

			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++) {
					ret[i][j] = a[i][j].real();
				}

			return	ret;
		}
		static Matrix<Real> GetImagPart(const Matrix<Complex>& a)
		{
			Matrix<Real> ret(a.RowNum(), a.ColNum());

			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++) {
					ret[i][j] = a[i][j].imag();
				}

			return	ret;
		}

		static Matrix<Complex> GetConjugateTranspose(const Matrix<Complex>& mat)
		{
			Matrix<Complex> ret(mat.ColNum(), mat.RowNum());

			for (int i = 0; i < mat.RowNum(); i++)
				for (int j = 0; j < mat.ColNum(); j++)
					ret[j][i] = std::conj(mat[i][j]);

			return ret;
		}
		static Matrix<Complex> CmplxMatFromRealMat(const Matrix<Real>& mat)
		{
			Matrix<Complex> mat_cmplx(mat.RowNum(), mat.ColNum());

			for (int i = 0; i < mat.RowNum(); i++)
				for (int j = 0; j < mat.ColNum(); j++)
					mat_cmplx[i][j] = Complex(mat(i, j), 0.0);

			return mat_cmplx;
		}

		static bool IsComplexMatReal(const Matrix<Complex>& mat)
		{
			for (int i = 0; i < mat.RowNum(); i++)
				for (int j = 0; j < mat.ColNum(); j++)
					if (mat[i][j].imag() != 0.0)
						return false;

			return true;
		}
		static bool IsHermitian(const Matrix<Complex>& mat)
		{
			if (mat.RowNum() != mat.ColNum())
				throw MatrixDimensionError("IsHermitian - matrix must be square", mat.RowNum(), mat.ColNum(), -1, -1);

			for (int i = 0; i < mat.RowNum(); i++)
				for (int j = i + 1; j < mat.ColNum(); j++)
					if (mat[i][j] != std::conj(mat[j][i]))
						return false;
			return true;
		}
		static bool IsUnitary(const Matrix<Complex>& mat)
		{
			// IsUnitary - complex square matrix U is unitary if its conjugate transpose U* is also its inverse
			if (mat.RowNum() != mat.ColNum())
				throw MatrixDimensionError("IsUnitary - matrix must be square", mat.RowNum(), mat.ColNum(), -1, -1);

			Matrix<Complex> matProd = mat * GetConjugateTranspose(mat);

			return matProd.IsUnit();
		}

		///////////////////       Matrix<Complex> - Matrix<Real>  operations     ///////////////////
		static Matrix<Complex> AddMat(const Matrix<Complex>& a, const Matrix<Real>& b)
		{
			if (a.RowNum() != b.RowNum() || a.ColNum() != b.ColNum())
				throw MatrixDimensionError("AddMat(Complex, Real) - must be same dim", a.RowNum(), a.ColNum(), b.RowNum(), b.ColNum());

			Matrix<Complex> ret(a);
			for (int i = 0; i < b.RowNum(); i++)
				for (int j = 0; j < b.ColNum(); j++)
					ret[i][j] += b[i][j];
			return ret;
		}
		static Matrix<Complex> AddMat(const Matrix<Real>& a, const Matrix<Complex>& b)
		{
			if (a.RowNum() != b.RowNum() || a.ColNum() != b.ColNum())
				throw MatrixDimensionError("AddMat(Real, Complex) - must be same dim", a.RowNum(), a.ColNum(), b.RowNum(), b.ColNum());

			Matrix<Complex> ret(b);
			for (int i = 0; i < a.RowNum(); i++)
				for (int j = 0; j < a.ColNum(); j++)
					ret[i][j] += a[i][j];
			return ret;
		}

		static Matrix<Complex> SubMat(const Matrix<Complex>& a, const Matrix<Real>& b)
		{
			if (a.RowNum() != b.RowNum() || a.ColNum() != b.ColNum())
				throw MatrixDimensionError("AddMat(Complex, Real) - must be same dim", a.RowNum(), a.ColNum(), b.RowNum(), b.ColNum());

			Matrix<Complex> ret(a);
			for (int i = 0; i < b.RowNum(); i++)
				for (int j = 0; j < b.ColNum(); j++)
					ret[i][j] -= b[i][j];
			return ret;
		}
		static Matrix<Complex> SubMat(const Matrix<Real>& a, const Matrix<Complex>& b)
		{
			if (a.RowNum() != b.RowNum() || a.ColNum() != b.ColNum())
				throw MatrixDimensionError("AddMat(Real, Complex) - must be same dim", a.RowNum(), a.ColNum(), b.RowNum(), b.ColNum());

			Matrix<Complex> ret(b);
			for (int i = 0; i < a.RowNum(); i++)
				for (int j = 0; j < a.ColNum(); j++)
					ret[i][j] = a[i][j] - b[i][j];
			return ret;
		}

		static Matrix<Complex> MulMat(const Complex& a, const Matrix<Real>& b)
		{
			Matrix<Complex>	ret(b.RowNum(), b.ColNum());

			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++) {
					ret[i][j] = a * b[i][j];
				}

			return	ret;
		}
		static Matrix<Complex> MulMat(const Matrix<Real>& a, const Complex& b)
		{
			Matrix<Complex>	ret(a.RowNum(), a.ColNum());

			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++) {
					ret[i][j] = b * a[i][j];
				}

			return	ret;
		}
		
		static Matrix<Complex> MulMat(const Matrix<Complex>& a, const Matrix<Real>& b)
		{
			if (a.ColNum() != b.RowNum())
				throw MatrixDimensionError("Matrix::operator*(Complex, Real)) - a.colNum must be equal to b.rowNum", a.RowNum(), a.ColNum(), b.RowNum(), b.ColNum());

			Matrix<Complex>	ret(a.RowNum(), b.ColNum());
			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++) {
					ret[i][j] = 0.0;
					for (int k = 0; k < a.ColNum(); k++)
						ret[i][j] += a[i][k] * b[k][j];
				}

			return	ret;
		}
		static Matrix<Complex> MulMat(const Matrix<Real>& a, const Matrix<Complex>& b)
		{
			if (a.ColNum() != b.RowNum())
				throw MatrixDimensionError("Matrix::operator*(Real, Complex)) - a.colNum must be equal to b.rowNum", a.RowNum(), a.ColNum(), b.RowNum(), b.ColNum());

			Matrix<Complex>	ret(a.RowNum(), b.ColNum());
			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++) {
					ret[i][j] = 0.0;
					for (int k = 0; k < a.ColNum(); k++)
						ret[i][j] += a[i][k] * b[k][j];
				}

			return	ret;
		}

		static Vector<Complex> MulMatVec(const Matrix<Real>& a, const Vector<Complex>& b)
		{
			if (a.ColNum() != b.size())
				throw MatrixDimensionError("operator*(Mat a, Vec b) - a.colNum must be equal to vector size", a.RowNum(), a.ColNum(), (int)b.size(), -1);

			Vector<Complex>	ret(a.RowNum());
			for (int i = 0; i < a.RowNum(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < a.ColNum(); j++)
					ret[i] += a[i][j] * b[j];
			}
			return ret;
		}
		static Vector<Complex> MulMatVec(const Matrix<Complex>& a, const Vector<Real>& b)
		{
			if (a.ColNum() != b.size())
				throw MatrixDimensionError("operator*(Mat a, Vec b) - a.colNum must be equal to vector size", a.RowNum(), a.ColNum(), (int)b.size(), -1);

			Vector<Complex>	ret(a.RowNum());
			for (int i = 0; i < a.RowNum(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < a.ColNum(); j++)
					ret[i] += a[i][j] * b[j];
			}
			return ret;
		}

		static Vector<Complex> MulVecMat(const Vector<Complex>& a, const Matrix<Real>& b)
		{
			if (a.size() != b.RowNum())
				throw MatrixDimensionError("operator*(Vec a, Mat b) - vector size must be equal to b.rowNum", (int)a.size(), -1, b.RowNum(), b.ColNum());

			Vector<Complex>	ret(b.ColNum());
			for (int i = 0; i < b.ColNum(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < b.RowNum(); j++)
					ret[i] += a[j] * b[j][i];
			}
			return ret;
		}
		static Vector<Complex> MulVecMat(const Vector<Real>& a, const Matrix<Complex>& b)
		{
			if (a.size() != b.RowNum())
				throw MatrixDimensionError("operator*(Vec a, Mat b) - vector size must be equal to b.rowNum", (int)a.size(), -1, b.RowNum(), b.ColNum());

			Vector<Complex>	ret(b.ColNum());
			for (int i = 0; i < b.ColNum(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < b.RowNum(); j++)
					ret[i] += a[j] * b[j][i];
			}
			return ret;
		}
	}
}
///////////////////////////   ./include/core/LinAlgEqSolvers.h   ///////////////////////////


namespace MML
{
	///////////////////////   GAUSS-JORDAN SOLVER    /////////////////////////////
	template<class Type>
	class GaussJordanSolver
	{
	public:
		static bool Solve(Matrix<Type>& a, Matrix<Type>& b)
		{
			int i, icol, irow, j, k, l, ll;
			Real big;
			Type dum, pivinv;

			int n = a.RowNum();
			int m = b.ColNum();
			std::vector<int> indxc(n), indxr(n), ipiv(n);
			for (j = 0; j < n; j++) ipiv[j] = 0;
			for (i = 0; i < n; i++) {
				big = 0.0;
				for (j = 0; j < n; j++)
					if (ipiv[j] != 1)
						for (k = 0; k < n; k++) {
							if (ipiv[k] == 0) {
								if (Abs(a[j][k]) >= big) {
									big = Abs(a[j][k]);
									irow = j;
									icol = k;
								}
							}
						}
				++(ipiv[icol]);
				if (irow != icol) {
					for (l = 0; l < n; l++) std::swap(a[irow][l], a[icol][l]);
					for (l = 0; l < m; l++) std::swap(b[irow][l], b[icol][l]);
				}
				indxr[i] = irow;
				indxc[i] = icol;

				if (a[icol][icol] == Real{ 0.0 })
					return false;
				// throw SingularMatrixError("GaussJordanSolver::Solve - Singular Matrix");

				pivinv = Real{ 1.0 } / a[icol][icol];
				a[icol][icol] = 1.0;
				for (l = 0; l < n; l++) a[icol][l] *= pivinv;
				for (l = 0; l < m; l++) b[icol][l] *= pivinv;
				for (ll = 0; ll < n; ll++)
					if (ll != icol) {
						dum = a[ll][icol];
						a[ll][icol] = 0.0;
						for (l = 0; l < n; l++) a[ll][l] -= a[icol][l] * dum;
						for (l = 0; l < m; l++) b[ll][l] -= b[icol][l] * dum;
					}
			}
			for (l = n - 1; l >= 0; l--) {
				if (indxr[l] != indxc[l])
					for (k = 0; k < n; k++)
						std::swap(a[k][indxr[l]], a[k][indxc[l]]);
			}

			return true;
		}
		static Vector<Type> Solve(Matrix<Type>& a, const Vector<Type>& b)
		{
			Vector<Type> x(b.size());
			if (Solve(b, x) == true)
				return x;
			else
				throw SingularMatrixError("GaussJordanSolver::Solve - Singular Matrix");
		}
		static bool Solve(Matrix<Type>& a, Vector<Type>& b)
		{
			Matrix<Type> bmat(b.size(), 1);
			for (int i = 0; i < b.size(); i++)
				bmat[i][0] = b[i];

			bool ret = Solve(a, bmat);
			b = bmat.VectorFromColumn(0);
			return ret;
		}
	};
}
// end namespace
///////////////////////////   ./include/core/Derivation.h   ///////////////////////////



namespace MML
{
	class Derivation
	{
	public:
		static inline const Real NDer1_h = 2 * std::sqrt(Constants::Epsilon);
		static inline const Real NDer2_h = std::pow(3 * Constants::Epsilon, 1.0 / 3.0);
		static inline const Real NDer4_h = std::pow(11.25 * Constants::Epsilon, 1.0 / 5.0);     // 0.0012009323661373839 for double!
		static inline const Real NDer6_h = std::pow(Constants::Epsilon / 168.0, 1.0 / 7.0);
		static inline const Real NDer8_h = std::pow(551.25 * Constants::Epsilon, 1.0 / 9.0);

		/********************************************************************************************************************/
		/********                               Numerical derivatives of FIRST order                                 ********/
		/********************************************************************************************************************/
		//////////////////////////              RealFunction            //////////////////////////
		static Real NDer1(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = f(x + h);
			Real y0 = f(x);
			Real diff = yh - y0;
			if (error)
			{
				Real ym = f(x - h);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;

				// h*|f''(x)|*0.5 + (|f(x+h)+|f(x)|) * Constants::Epsilon/h
				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
			}
			return diff / h;
		}
		static Real NDer1(const IRealFunction& f, Real x, Real* error)
		{
			// Error bound ~eps^1/2
			// Note that this estimate of h differs from the best estimate by a factor of sqrt((|f(x)| + |f(x+h)|)/|f''(x)|).
			// Since this factor is invariant under the scaling f -> kf, then we are somewhat justified in approximating it by 1.
			// This approximation will get better as we move to higher orders of accuracy.
			return NDer1(f, x, NDer1_h, error);
		}
		static Real NDer1(const IRealFunction& f, Real x)
		{
			return NDer1(f, x, NDer1_h, nullptr);
		}
		
		static Real NDer1Left(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer1(f, x - 2 * NDer1_h, NDer1_h, error); }
		static Real NDer1Right(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer1(f, x + 2 * NDer1_h, NDer1_h, error); }
		static Real NDer1Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer1(f, x - 2 * h, h, error); }
		static Real NDer1Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer1(f, x + 2 * h, h, error); }

		static Real NSecDer1(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NSecDer1(f, x, NDer1_h, error);
		}
		static Real NSecDer1(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NDer2(f, x + h, h, error);
			Real y0 = NDer2(f, x, h, error);
			Real diff = yh - y0;
			if (error)
			{
				Real ym = NDer2(f, x - h, h, error);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;

				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
			}
			return diff / h;
		}

		static Real NThirdDer1(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NThirdDer1(f, x, NDer1_h, error);
		}
		static Real NThirdDer1(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NSecDer2(f, x + h, h, error);
			Real y0 = NSecDer2(f, x, h, error);
			Real diff = yh - y0;
			if (error)
			{
				Real ym = NSecDer2(f, x - h, h, error);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;

				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
			}
			return diff / h;
		}

		//////////////////////////             ScalarFunction           //////////////////////////
		template <int N>
		static Real NDer1Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer1Partial(f, deriv_index, point, NDer1_h, error);
		}

		template <int N>
		static Real NDer1Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x = point;
			Real y0 = f(x);

			x[deriv_index] = orig_x + h;
			Real yh = f(x);

			Real diff = yh - y0;
			if (error)
			{
				x[deriv_index] = orig_x - h;
				Real ym = f(x);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;
				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
			}
			return diff / h;
		}

		template <int N>
		static Real NSecDer1Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NSecDer1Partial(f, der_ind1, der_ind2, point, NDer1_h, error);
		}

		template <int N>
		static Real NSecDer1Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real x_orig_val = point[der_ind2];

			auto x_eval_pos = point;
			Real y0 = NDer2Partial(f, der_ind1, x_eval_pos, error);
			x_eval_pos[der_ind2] = x_orig_val + h;
			Real yh = NDer2Partial(f, der_ind1, x_eval_pos, error);

			Real diff = yh - y0;
			if (error)
			{
				x_eval_pos[der_ind2] = x_orig_val - h;

				Real ym = NDer2Partial(f, der_ind1, x_eval_pos, error);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;

				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
			}
			return diff / h;
		}

		template <int N>
		static VectorN<Real, N> NDer1PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer1PartialByAll(f, point, NDer1_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer1PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer1Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer1Partial(f, i, point, h);
			}

			return ret;
		}

		//////////////////////////             VectorFunction           //////////////////////////
		template <int N>
		static Real NDer1Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer1Partial(f, func_index, deriv_index, point, NDer1_h, error);
		}

		template <int N>
		static Real NDer1Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			auto x = point;

			Real x_orig = x[deriv_index];
			Real y0 = f(x)[func_index];

			x[deriv_index] = x_orig + h;
			Real yh = f(x)[func_index];

			Real diff = yh - y0;
			if (error)
			{
				x[deriv_index] = x_orig - h;
				Real ym = f(x)[func_index];
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;
				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
			}
			return diff / h;
		}

		template <int N>
		static VectorN<Real, N> NDer1PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer1PartialByAll(f, func_index, point, NDer1_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer1PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer1Partial(f, func_index, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer1Partial(f, func_index, i, point, h);
			}

			return ret;
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer1PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, MatrixNM<Real, N, N>* error = nullptr)
		{
			return NDer1PartialAllByAll(f, point, NDer1_h, error);
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer1PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, Real h, MatrixNM<Real, N, N>* error = nullptr)
		{
			MatrixNM<Real, N, N> ret;

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer1Partial(f, i, j, point, h, &((*error)(i, j)));
					else
						ret(i, j) = NDer1Partial(f, i, j, point, h);
				}

			return ret;
		}

		//////////////////////////             TensorField           //////////////////////////
		template <int N>
		static Real NDer1Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer1Partial(f, i, j, deriv_index, point, NDer1_h, error);
		}

		template <int N>
		static Real NDer1Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			auto x = point;

			Real x_orig = x[deriv_index];
			Real y0 = f.Component(i, j, x);

			x[deriv_index] = x_orig + h;
			Real yh = f.Component(i, j, x);

			Real diff = yh - y0;
			if (error)
			{
				x[deriv_index] = x_orig - h;
				Real ym = f.Component(i, j, x);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;
				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
			}
			return diff / h;
		}

		template <int N>
		static Real NDer1Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer1Partial(f, i, j, k, deriv_index, point, NDer1_h, error);
		}

		template <int N>
		static Real NDer1Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			auto x = point;

			Real x_orig = x[deriv_index];
			Real y0 = f.Component(i, j, k, x);

			x[deriv_index] = x_orig + h;
			Real yh = f.Component(i, j, k, x);

			Real diff = yh - y0;
			if (error)
			{
				x[deriv_index] = x_orig - h;
				Real ym = f.Component(i, j, k, x);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;
				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
			}
			return diff / h;
		}

		template <int N>
		static Real NDer1Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer1Partial(f, i, j, k, l, deriv_index, point, NDer1_h, error);
		}

		template <int N>
		static Real NDer1Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			auto x = point;

			Real x_orig = x[deriv_index];
			Real y0 = f.Component(i, j, k, l, x);

			x[deriv_index] = x_orig + h;
			Real yh = f.Component(i, j, k, l, x);

			Real diff = yh - y0;
			if (error)
			{
				x[deriv_index] = x_orig - h;
				Real ym = f.Component(i, j, k, l, x);
				Real ypph = std::abs(yh - 2 * y0 + ym) / h;
				*error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
			}
			return diff / h;
		}

		/////////////////////////             ParametricCurve           /////////////////////////
		template <int N>
		static VectorN<Real, N> NDer1(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer1(f, t, NDer1_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer1(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(t + h);
			VectorN<Real, N> y0 = f(t);
			VectorN<Real, N> diff = yh - y0;

			if (error)
			{
				VectorN<Real, N> ym = f(t - h);
				VectorN<Real, N> ypph_vec = yh - 2 * y0 + ym;

				Real ypph = ypph_vec.NormL2() / h;

				*error = ypph / 2 + (yh.NormL2() + y0.NormL2()) * Constants::Epsilon / h;
			}
			return diff / h;
		}

		template <int N>
		static VectorN<Real, N> NSecDer1(const IParametricCurve<N>& f, Real x, Real* error = nullptr)
		{
			return NSecDer1(f, x, NDer1_h, error);
		}

		template <int N>
		static VectorN<Real, N> NSecDer1(const IParametricCurve<N>& f, Real x, Real h, Real* error = nullptr)
		{
			VectorN<Real, N>  yh = NDer2(f, x + h, h, error);
			VectorN<Real, N>  y0 = NDer2(f, x, h, error);
			VectorN<Real, N>  diff = yh - y0;
			if (error)
			{
				VectorN<Real, N> ym = NDer2(f, x - h, h, error);
				VectorN<Real, N> ypph_vec = (yh - 2 * y0 + ym) / h;

				Real ypph = ypph_vec.NormL2();

				*error = ypph / 2 + (yh.NormL2() + y0.NormL2()) * Constants::Epsilon / h;
			}
			return diff / h;
		}

		template <int N>
		static VectorN<Real, N> NThirdDer1(const IParametricCurve<N>& f, Real x, Real* error = nullptr)
		{
			return NThirdDer1(f, x, NDer1_h, error);
		}

		template <int N>
		static VectorN<Real, N> NThirdDer1(const IParametricCurve<N>& f, Real x, Real h, Real* error = nullptr)
		{
			VectorN<Real, N>  yh = NSecDer2(f, x + h, h, error);
			VectorN<Real, N>  y0 = NSecDer2(f, x, h, error);
			VectorN<Real, N>  diff = yh - y0;
			if (error)
			{
				VectorN<Real, N> ym = NSecDer2(f, x - h, h, error);
				VectorN<Real, N> ypph_vec = (yh - 2 * y0 + ym) / h;

				Real ypph = ypph_vec.NormL2();

				*error = ypph / 2 + (yh.NormL2() + y0.NormL2()) * Constants::Epsilon / h;
			}
			return diff / h;
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SECOND order                                 ********/
		/********************************************************************************************************************/

		//////////////////////////              RealFunction            //////////////////////////
		static Real NDer2(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = f(x + h);
			Real ymh = f(x - h);
			Real diff = yh - ymh;
			if (error)
			{
				Real y2h = f(x + 2 * h);
				Real ym2h = f(x - 2 * h);
				*error = Constants::Epsilon * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}
		static Real NDer2(const IRealFunction& f, Real x, Real* error)
		{
			// Error bound ~eps^2/3
			// See the previous discussion to understand determination of h and the error bound.
			// Series[(f[x+h] - f[x-h])/(2*h), {h, 0, 4}]

			return NDer2(f, x, NDer2_h, error);
		}
		static Real NDer2(const IRealFunction& f, Real x)
		{
			return NDer2(f, x, NDer2_h, nullptr);
		}
		
		static Real NDer2Left(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer2(f, x - 2 * NDer2_h, NDer2_h, error); }
		static Real NDer2Right(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer2(f, x + 2 * NDer2_h, NDer2_h, error); }
		static Real NDer2Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer2(f, x - 3 * h, h, error); }
		static Real NDer2Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer2(f, x + 3 * h, h, error); }

		static Real NSecDer2(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NSecDer2(f, x, NDer2_h, error);
		}
		static Real NSecDer2(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NDer4(f, x + h, error);
			Real ymh = NDer4(f, x - h, error);
			Real diff = yh - ymh;
			if (error)
			{
				Real y2h = NDer4(f, x + 2 * h, error);
				Real ym2h = NDer4(f, x - 2 * h, error);

				*error = Constants::Epsilon * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}

		static Real NThirdDer2(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NThirdDer2(f, x, NDer2_h, error);
		}
		static Real NThirdDer2(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NSecDer4(f, x + h, error);
			Real ymh = NSecDer4(f, x - h, error);
			Real diff = yh - ymh;
			if (error)
			{
				Real y2h = NSecDer4(f, x + 2 * h, error);
				Real ym2h = NSecDer4(f, x - 2 * h, error);

				*error = Constants::Epsilon * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}

		//////////////////////////             ScalarFunction           //////////////////////////
		template <int N>
		static Real NDer2Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, deriv_index, point, NDer2_h, error);
		}

		template <int N>
		static Real NDer2Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			auto    x = point;
			x[deriv_index] = orig_x + h;
			Real yh = f(x);

			x[deriv_index] = orig_x - h;
			Real ymh = f(x);

			Real diff = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 2 * h;
				Real y2h = f(x);

				x[deriv_index] = orig_x - 2 * h;
				Real ym2h = f(x);

				*error = Constants::Epsilon * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}

		template <int N>
		static Real NSecDer2Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NSecDer2Partial(f, der_ind1, der_ind2, point, NDer2_h, error);
		}

		template <int N>
		static Real NSecDer2Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real orig_x = point[der_ind2];
			auto x_eval_pos = point;

			x_eval_pos[der_ind2] = orig_x + h;
			Real yh = NDer4Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - h;
			Real ymh = NDer4Partial(f, der_ind1, x_eval_pos, error);

			Real diff = yh - ymh;

			if (error)
			{
				x_eval_pos[der_ind2] = orig_x + 2 * h;
				Real y2h = NDer4Partial(f, der_ind1, x_eval_pos, error);

				x_eval_pos[der_ind2] = orig_x - 2 * h;
				Real ym2h = NDer4Partial(f, der_ind1, x_eval_pos, error);

				*error = Constants::Epsilon * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer2PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer2PartialByAll(f, point, NDer2_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer2PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer2Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer2Partial(f, i, point, h);
			}

			return ret;
		}

		//////////////////////////             VectorFunction           //////////////////////////
		template <int N>
		static Real NDer2Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, func_index, deriv_index, point, NDer2_h, error);
		}

		template <int N>
		static Real NDer2Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f(x)[func_index];

			x[deriv_index] = orig_x - h;
			Real ymh = f(x)[func_index];

			Real diff = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 2 * h;
				Real y2h = f(x)[func_index];

				x[deriv_index] = orig_x - 2 * h;
				Real ym2h = f(x)[func_index];

				*error = Constants::Epsilon * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer2PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer2PartialByAll(f, func_index, point, NDer2_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer2PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer2Partial(f, func_index, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer2Partial(f, func_index, i, point, h);
			}

			return ret;
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer2PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, MatrixNM<Real, N, N>* error = nullptr)
		{
			return NDer2PartialAllByAll(f, point, NDer2_h, error);
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer2PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, Real h, MatrixNM<Real, N, N>* error = nullptr)
		{
			MatrixNM<Real, N, N> ret;

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer2Partial(f, i, j, point, h, &((*error)(i, j)));
					else
						ret(i, j) = NDer2Partial(f, i, j, point, h);
				}

			return ret;
		}

		//////////////////////////             TensorField           //////////////////////////
		template <int N>
		static Real NDer2Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, i, j, deriv_index, point, NDer2_h, error);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f.Component(i, j, x);

			x[deriv_index] = orig_x - h;
			Real ymh = f.Component(i, j, x);

			Real diff = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 2 * h;
				Real y2h = f.Component(i, j, x);

				x[deriv_index] = orig_x - 2 * h;
				Real ym2h = f.Component(i, j, x);

				*error = Constants::Epsilon * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, i, j, k, deriv_index, point, NDer2_h, error);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f.Component(i, j, k, x);

			x[deriv_index] = orig_x - h;
			Real ymh = f.Component(i, j, k, x);

			Real diff = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 2 * h;
				Real y2h = f.Component(i, j, k, x);

				x[deriv_index] = orig_x - 2 * h;
				Real ym2h = f.Component(i, j, k, x);

				*error = Constants::Epsilon * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer2Partial(f, i, j, k, l, deriv_index, point, NDer2_h, error);
		}

		template <int N>
		static Real NDer2Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f.Component(i, j, k, l, x);

			x[deriv_index] = orig_x - h;
			Real ymh = f.Component(i, j, k, l, x);

			Real diff = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 2 * h;
				Real y2h = f.Component(i, j, k, l, x);

				x[deriv_index] = orig_x - 2 * h;
				Real ym2h = f.Component(i, j, k, l, x);

				*error = Constants::Epsilon * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
			}

			return diff / (2 * h);
		}
		/////////////////////////             ParametricCurve           /////////////////////////
		template <int N>
		static VectorN<Real, N> NDer2(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer2(f, t, NDer2_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer2(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(t + h);
			VectorN<Real, N> ymh = f(t - h);
			VectorN<Real, N> diff = yh - ymh;

			if (error)
			{
				VectorN<Real, N> yth = f(t + 2 * h);
				VectorN<Real, N> ymth = f(t - 2 * h);

				*error = Constants::Epsilon * ((yh + ymh) / (2 * h)).NormL2() + std::abs(((yth - ymth) / 2 - diff).NormL2()) / (6 * h);
			}
			return diff / (2 * h);
		}

		template <int N>
		static VectorN<Real, N> NSecDer2(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NSecDer2(f, t, NDer2_h, error);
		}

		template <int N>
		static VectorN<Real, N> NSecDer2(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = NDer4(f, t + h, error);
			VectorN<Real, N> ymh = NDer4(f, t - h, error);
			VectorN<Real, N> diff = yh - ymh;

			if (error)
			{
				VectorN<Real, N> yth = NDer4(f, t + 2 * h, error);
				VectorN<Real, N> ymth = NDer4(f, t - 2 * h, error);

				*error = Constants::Epsilon * ((yh + ymh) / (2 * h)).NormL2() + std::abs(((yth - ymth) / 2 - diff).NormL2()) / (6 * h);
			}
			return diff / (2 * h);
		}

		template <int N>
		static VectorN<Real, N> NThirdDer2(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NThirdDer2(f, t, NDer2_h, error);
		}

		template <int N>
		static VectorN<Real, N> NThirdDer2(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = NSecDer4(f, t + h, error);
			VectorN<Real, N> ymh = NSecDer4(f, t - h, error);
			VectorN<Real, N> diff = yh - ymh;

			if (error)
			{
				VectorN<Real, N> yth = NSecDer4(f, t + 2 * h, error);
				VectorN<Real, N> ymth = NSecDer4(f, t - 2 * h, error);

				*error = Constants::Epsilon * ((yh + ymh) / (2 * h)).NormL2() + std::abs(((yth - ymth) / 2 - diff).NormL2()) / (6 * h);
			}
			return diff / (2 * h);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of FOURTH order                                 ********/
		/********************************************************************************************************************/

		//////////////////////////              RealFunction            //////////////////////////
		static Real NDer4(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = f(x + h);
			Real ymh = f(x - h);
			Real y2h = f(x + 2 * h);
			Real ym2h = f(x - 2 * h);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				// Mathematica code to extract the remainder:
				// Series[(f[x-2*h]+ 8*f[x+h] - 8*f[x-h] - f[x+2*h])/(12*h), {h, 0, 7}]
				Real y3h = f(x + 3 * h);
				Real ym3h = f(x - 3 * h);

				// Error from fifth derivative:
				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				// Error from function evaluation:
				*error += Constants::Epsilon * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}
		static Real NDer4(const IRealFunction& f, Real x, Real* error)
		{
			// Error bound ~eps^4/5
			return NDer4(f, x, NDer4_h, error);
		}
		static Real NDer4(const IRealFunction& f, Real x)
		{
			return NDer4(f, x, NDer4_h, nullptr);
		}

		static Real NDer4Left(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer4(f, x - 4 * NDer4_h, NDer4_h, error); }
		static Real NDer4Right(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer4(f, x + 4 * NDer4_h, NDer4_h, error); }
		static Real NDer4Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer4(f, x - 4 * h, h, error); }
		static Real NDer4Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer4(f, x + 4 * h, h, error); }

		static Real NSecDer4(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NSecDer4(f, x, NDer4_h, error);
		}
		static Real NSecDer4(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NDer6(f, x + h, error);
			Real ymh = NDer6(f, x - h, error);
			Real y2h = NDer6(f, x + 2 * h, error);
			Real ym2h = NDer6(f, x - 2 * h, error);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				Real y3h = NDer6(f, x + 3 * h, error);
				Real ym3h = NDer6(f, x - 3 * h, error);

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Epsilon * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}

		static Real NThirdDer4(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NThirdDer4(f, x, NDer4_h, error);
		}
		static Real NThirdDer4(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NSecDer6(f, x + h, error);
			Real ymh = NSecDer6(f, x - h, error);
			Real y2h = NSecDer6(f, x + 2 * h, error);
			Real ym2h = NSecDer6(f, x - 2 * h, error);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				Real y3h = NSecDer6(f, x + 3 * h, error);
				Real ym3h = NSecDer6(f, x - 3 * h, error);

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Epsilon * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}

		//////////////////////////             ScalarFunction           //////////////////////////
		template <int N>
		static Real NDer4Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer4Partial(f, deriv_index, point, NDer4_h, error);
		}

		template <int N>
		static Real NDer4Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f(x);

			x[deriv_index] = orig_x - h;
			Real ymh = f(x);

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f(x);

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f(x);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 3 * h;
				Real y3h = f(x);

				x[deriv_index] = orig_x - 3 * h;
				Real ym3h = f(x);

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Epsilon * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}

		template <int N>
		static Real NSecDer4Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NSecDer4Partial(f, der_ind1, der_ind2, point, NDer4_h, error);
		}

		template <int N>
		static Real NSecDer4Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[der_ind2];
			auto x_eval_pos = point;

			x_eval_pos[der_ind2] = orig_x + h;
			Real yh = NDer6Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - h;
			Real ymh = NDer6Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x + 2 * h;
			Real y2h = NDer6Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - 2 * h;
			Real ym2h = NDer6Partial(f, der_ind1, x_eval_pos, error);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				x_eval_pos[der_ind2] = orig_x + 3 * h;
				Real y3h = NDer6Partial(f, der_ind1, x_eval_pos, error);

				x_eval_pos[der_ind2] = orig_x - 3 * h;
				Real ym3h = NDer6Partial(f, der_ind1, x_eval_pos, error);

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Epsilon * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer4PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer4PartialByAll(f, point, NDer4_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer4PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer4Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer4Partial(f, i, point, h);
			}

			return ret;
		}

		//////////////////////////             VectorFunction           //////////////////////////
		template <int N>
		static Real NDer4Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer4Partial(f, func_index, deriv_index, point, NDer4_h, error);
		}

		template <int N>
		static Real NDer4Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f(x)[func_index];

			x[deriv_index] = orig_x - h;
			Real ymh = f(x)[func_index];

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f(x)[func_index];

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f(x)[func_index];

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 3 * h;
				Real y3h = f(x)[func_index];

				x[deriv_index] = orig_x - 3 * h;
				Real ym3h = f(x)[func_index];

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Epsilon * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer4PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer4PartialByAll(f, func_index, point, NDer4_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer4PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer4Partial(f, func_index, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer4Partial(f, func_index, i, point, h);
			}

			return ret;
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer4PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, MatrixNM<Real, N, N>* error = nullptr)
		{
			return NDer4PartialAllByAll(f, point, NDer4_h, error);
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer4PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, Real h, MatrixNM<Real, N, N>* error = nullptr)
		{
			MatrixNM<Real, N, N> ret;

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer4Partial(f, i, j, point, h, &((*error)(i, j)));
					else
						ret(i, j) = NDer4Partial(f, i, j, point, h);
				}

			return ret;
		}

		//////////////////////////             TensorField           //////////////////////////
		template <int N>
		static Real NDer4Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer4Partial(f, i, j, deriv_index, point, NDer4_h, error);
		}

		template <int N>
		static Real NDer4Partial(const ITensorField2<N>& f, int i, int j, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f.Component(i, j, x);

			x[deriv_index] = orig_x - h;
			Real ymh = f.Component(i, j, x);

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f.Component(i, j, x);

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f.Component(i, j, x);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 3 * h;
				Real y3h = f.Component(i, j, x);

				x[deriv_index] = orig_x - 3 * h;
				Real ym3h = f.Component(i, j, x);

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Epsilon * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}

		template <int N>
		static Real NDer4Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer4Partial(f, i, j, k, deriv_index, point, NDer4_h, error);
		}

		template <int N>
		static Real NDer4Partial(const ITensorField3<N>& f, int i, int j, int k, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f.Component(i, j, k, x);

			x[deriv_index] = orig_x - h;
			Real ymh = f.Component(i, j, k, x);

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f.Component(i, j, k, x);

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f.Component(i, j, k, x);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 3 * h;
				Real y3h = f.Component(i, j, k, x);

				x[deriv_index] = orig_x - 3 * h;
				Real ym3h = f.Component(i, j, k, x);

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Epsilon * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}


		template <int N>
		static Real NDer4Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer4Partial(f, i, j, k, l, deriv_index, point, NDer4_h, error);
		}

		template <int N>
		static Real NDer4Partial(const ITensorField4<N>& f, int i, int j, int k, int l, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };
			x[deriv_index] = orig_x + h;
			Real yh = f.Component(i, j, k, l, x);

			x[deriv_index] = orig_x - h;
			Real ymh = f.Component(i, j, k, l, x);

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f.Component(i, j, k, l, x);

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f.Component(i, j, k, l, x);

			Real y2 = ym2h - y2h;
			Real y1 = yh - ymh;

			if (error)
			{
				x[deriv_index] = orig_x + 3 * h;
				Real y3h = f.Component(i, j, k, l, x);

				x[deriv_index] = orig_x - 3 * h;
				Real ym3h = f.Component(i, j, k, l, x);

				*error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
				*error += Constants::Epsilon * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}

		/////////////////////////             ParametricCurve           /////////////////////////
		template <int N>
		static VectorN<Real, N> NDer4(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer4(f, t, NDer4_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer4(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(t + h);
			VectorN<Real, N> ymh = f(t - h);
			VectorN<Real, N> y2h = f(t + 2 * h);
			VectorN<Real, N> ym2h = f(t - 2 * h);

			VectorN<Real, N> y2 = ym2h - y2h;
			VectorN<Real, N> y1 = yh - ymh;

			if (error)
			{
				VectorN<Real, N> y3h = f(t + 3 * h);
				VectorN<Real, N> ym3h = f(t - 3 * h);

				*error = std::abs((y3h - ym3h).NormL2() / 2 + 2 * (ym2h - y2h).NormL2() + 5 * (yh - ymh).NormL2() / 2) / (30 * h);
				*error += Constants::Epsilon * (y2h.NormL2() + ym2h.NormL2() + 8 * (ymh.NormL2() + yh.NormL2())) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}

		template <int N>
		static VectorN<Real, N> NSecDer4(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NSecDer4(f, t, NDer4_h, error);
		}

		template <int N>
		static VectorN<Real, N> NSecDer4(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = NDer6(f, t + h, error);
			VectorN<Real, N> ymh = NDer6(f, t - h, error);
			VectorN<Real, N> y2h = NDer6(f, t + 2 * h, error);
			VectorN<Real, N> ym2h = NDer6(f, t - 2 * h, error);

			VectorN<Real, N> y2 = ym2h - y2h;
			VectorN<Real, N> y1 = yh - ymh;

			if (error)
			{
				VectorN<Real, N> y3h = NDer6(f, t + 3 * h, error);
				VectorN<Real, N> ym3h = NDer6(f, t - 3 * h, error);

				*error = std::abs((y3h - ym3h).NormL2() / 2 + 2 * (ym2h - y2h).NormL2() + 5 * (yh - ymh).NormL2() / 2) / (30 * h);
				*error += Constants::Epsilon * (y2h.NormL2() + ym2h.NormL2() + 8 * (ymh.NormL2() + yh.NormL2())) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}


		template <int N>
		static VectorN<Real, N> NThirdDer4(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NThirdDer4(f, t, NDer4_h, error);
		}

		template <int N>
		static VectorN<Real, N> NThirdDer4(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = NSecDer6(f, t + h, error);
			VectorN<Real, N> ymh = NSecDer6(f, t - h, error);
			VectorN<Real, N> y2h = NSecDer6(f, t + 2 * h, error);
			VectorN<Real, N> ym2h = NSecDer6(f, t - 2 * h, error);

			VectorN<Real, N> y2 = ym2h - y2h;
			VectorN<Real, N> y1 = yh - ymh;

			if (error)
			{
				VectorN<Real, N> y3h = NSecDer6(f, t + 3 * h, error);
				VectorN<Real, N> ym3h = NSecDer6(f, t - 3 * h, error);

				*error = std::abs((y3h - ym3h).NormL2() / 2 + 2 * (ym2h - y2h).NormL2() + 5 * (yh - ymh).NormL2() / 2) / (30 * h);
				*error += Constants::Epsilon * (y2h.NormL2() + ym2h.NormL2() + 8 * (ymh.NormL2() + yh.NormL2())) / (12 * h);
			}
			return (y2 + 8 * y1) / (12 * h);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of SIXTH order                                 ********/
		/********************************************************************************************************************/

		//////////////////////////              RealFunction            //////////////////////////
		static Real NDer6(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			const Real eps = (std::numeric_limits<Real>::epsilon)();

			Real yh = f(x + h);
			Real ymh = f(x - h);
			Real y1 = yh - ymh;
			Real y2 = f(x - 2 * h) - f(x + 2 * h);
			Real y3 = f(x + 3 * h) - f(x - 3 * h);

			if (error)
			{
				// Mathematica code to generate fd scheme for 7th derivative:
				// Sum[(-1)^i*Binomial[7, i]*(f[x+(3-i)*h] + f[x+(4-i)*h])/2, {i, 0, 7}]
				// Mathematica to demonstrate that this is a finite difference formula for 7th derivative:
				// Series[(f[x+4*h]-f[x-4*h] + 6*(f[x-3*h] - f[x+3*h]) + 14*(f[x-h] - f[x+h] + f[x+2*h] - f[x-2*h]))/2, {h, 0, 15}]
				Real y7 = (f(x + 4 * h) - f(x - 4 * h) - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}
		static Real NDer6(const IRealFunction& f, Real x, Real* error)
		{
			// Error bound ~eps^6/7
			// Error: h^6f^(7)(x)/140 + 5|f(x)|eps/h
			return NDer6(f, x, NDer6_h, error);
		}
		static Real NDer6(const IRealFunction& f, Real x)
		{
			// Error bound ~eps^6/7
			// Error: h^6f^(7)(x)/140 + 5|f(x)|eps/h
			return NDer6(f, x, NDer6_h, nullptr);
		}

		static Real NDer6Left(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer6(f, x - 5 * NDer6_h, NDer6_h, error); }
		static Real NDer6Right(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer6(f, x + 5 * NDer6_h, NDer6_h, error); }
		static Real NDer6Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer6(f, x - 5 * h, h, error); }
		static Real NDer6Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer6(f, x + 5 * h, h, error); }

		static Real NSecDer6(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NSecDer6(f, x, NDer6_h, error);
		}
		static Real NSecDer6(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NDer8(f, x + h, error);
			Real ymh = NDer8(f, x - h, error);
			Real y1 = yh - ymh;
			Real y2 = NDer8(f, x - 2 * h, error) - NDer8(f, x + 2 * h, error);
			Real y3 = NDer8(f, x + 3 * h, error) - NDer8(f, x - 3 * h, error);

			if (error)
			{
				Real y7 = (NDer8(f, x + 4 * h, error) - NDer8(f, x - 4 * h, error) - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}

		static Real NThirdDer6(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NThirdDer6(f, x, NDer6_h, error);
		}
		static Real NThirdDer6(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NSecDer8(f, x + h, error);
			Real ymh = NSecDer8(f, x - h, error);
			Real y1 = yh - ymh;
			Real y2 = NSecDer8(f, x - 2 * h, error) - NSecDer8(f, x + 2 * h, error);
			Real y3 = NSecDer8(f, x + 3 * h, error) - NSecDer8(f, x - 3 * h, error);

			if (error)
			{
				Real y7 = (NSecDer8(f, x + 4 * h, error) - NSecDer8(f, x - 4 * h, error) - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}

		//////////////////////////             ScalarFunction           //////////////////////////
		template <int N>
		static Real NDer6Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer6Partial(f, deriv_index, point, NDer6_h, error);
		}

		template <int N>
		static Real NDer6Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };

			x[deriv_index] = orig_x + h;
			Real yh = f(x);

			x[deriv_index] = orig_x - h;
			Real ymh = f(x);

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f(x);

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f(x);

			x[deriv_index] = orig_x + 3 * h;
			Real y3h = f(x);

			x[deriv_index] = orig_x - 3 * h;
			Real ym3h = f(x);

			Real y1 = yh - ymh;
			Real y2 = ym2h - y2h;
			Real y3 = y3h - ym3h;

			if (error)
			{
				x[deriv_index] = orig_x + 4 * h;
				Real y4h = f(x);

				x[deriv_index] = orig_x - 4 * h;
				Real ym4h = f(x);

				Real y7 = (y4h - ym4h - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}

		template <int N>
		static Real NSecDer6Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NSecDer6Partial(f, der_ind1, der_ind2, point, NDer6_h, error);
		}

		template <int N>
		static Real NSecDer6Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[der_ind2];
			auto x_eval_pos = point;

			x_eval_pos[der_ind2] = orig_x + h;
			Real yh = NDer6Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - h;
			Real ymh = NDer6Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x + 2 * h;
			Real y2h = NDer6Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - 2 * h;
			Real ym2h = NDer6Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x + 3 * h;
			Real y3h = NDer6Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - 3 * h;
			Real ym3h = NDer6Partial(f, der_ind1, x_eval_pos, error);

			Real y1 = yh - ymh;
			Real y2 = ym2h - y2h;
			Real y3 = y3h - ym3h;

			if (error)
			{
				x_eval_pos[der_ind2] = orig_x + 4 * h;
				Real y4h = NDer6Partial(f, der_ind1, x_eval_pos, error);

				x_eval_pos[der_ind2] = orig_x - 4 * h;
				Real ym4h = NDer6Partial(f, der_ind1, x_eval_pos, error);

				Real y7 = (y4h - ym4h - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer6PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer6PartialByAll(f, point, NDer6_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer6PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer6Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer6Partial(f, i, point, h);
			}

			return ret;
		}

		//////////////////////////             VectorFunction           //////////////////////////
		template <int N>
		static Real NDer6Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer6Partial(f, func_index, deriv_index, point, NDer6_h, error);
		}

		template <int N>
		static Real NDer6Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };

			x[deriv_index] = orig_x + h;
			Real yh = f(x)[func_index];

			x[deriv_index] = orig_x - h;
			Real ymh = f(x)[func_index];

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f(x)[func_index];

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f(x)[func_index];

			x[deriv_index] = orig_x + 3 * h;
			Real y3h = f(x)[func_index];

			x[deriv_index] = orig_x - 3 * h;
			Real ym3h = f(x)[func_index];

			Real y1 = yh - ymh;
			Real y2 = ym2h - y2h;
			Real y3 = y3h - ym3h;

			if (error)
			{
				x[deriv_index] = orig_x + 4 * h;
				Real y4h = f(x)[func_index];

				x[deriv_index] = orig_x - 4 * h;
				Real ym4h = f(x)[func_index];

				Real y7 = (y4h - ym4h - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer6PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer6PartialByAll(f, func_index, point, NDer6_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer6PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer6Partial(f, func_index, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer6Partial(f, func_index, i, point, h);
			}

			return ret;
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer6PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, MatrixNM<Real, N, N>* error = nullptr)
		{
			return NDer6PartialAllByAll(f, point, NDer6_h, error);
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer6PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, Real h, MatrixNM<Real, N, N>* error = nullptr)
		{
			MatrixNM<Real, N, N> ret;

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer6Partial(f, i, j, point, h, &((*error)(i, j)));
					else
						ret(i, j) = NDer6Partial(f, i, j, point, h);
				}

			return ret;
		}

		/////////////////////////             ParametricCurve           /////////////////////////
		template <int N>
		static VectorN<Real, N> NDer6(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer6(f, t, NDer6_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer6(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(t + h);
			VectorN<Real, N> ymh = f(t - h);
			VectorN<Real, N> y1 = yh - ymh;
			VectorN<Real, N> y2 = f(t - 2 * h) - f(t + 2 * h);
			VectorN<Real, N> y3 = f(t + 3 * h) - f(t - 3 * h);

			if (error)
			{
				VectorN<Real, N> y7 = (f(t + 4 * h) - f(t - 4 * h) - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = y7.NormL2() / (140 * h) + 5 * (yh.NormL2() + ymh.NormL2()) * Constants::Epsilon / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}

		template <int N>
		static VectorN<Real, N> NSecDer6(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NSecDer6(f, t, NDer6_h, error);
		}

		template <int N>
		static VectorN<Real, N> NSecDer6(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = NDer8(f, t + h, error);
			VectorN<Real, N> ymh = NDer8(f, t - h, error);
			VectorN<Real, N> y1 = yh - ymh;
			VectorN<Real, N> y2 = NDer8(f, t - 2 * h, error) - NDer8(f, t + 2 * h, error);
			VectorN<Real, N> y3 = NDer8(f, t + 3 * h, error) - NDer8(f, t - 3 * h, error);

			if (error)
			{
				VectorN<Real, N> y7 = (NDer8(f, t + 4 * h, error) - NDer8(f, t - 4 * h, error) - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = y7.NormL2() / (140 * h) + 5 * (yh.NormL2() + ymh.NormL2()) * Constants::Epsilon / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}

		template <int N>
		static VectorN<Real, N> NThirdDer6(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NThirdDer6(f, t, NDer6_h, error);
		}

		template <int N>
		static VectorN<Real, N> NThirdDer6(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = NSecDer8(f, t + h, error);
			VectorN<Real, N> ymh = NSecDer8(f, t - h, error);
			VectorN<Real, N> y1 = yh - ymh;
			VectorN<Real, N> y2 = NSecDer8(f, t - 2 * h, error) - NSecDer8(f, t + 2 * h, error);
			VectorN<Real, N> y3 = NSecDer8(f, t + 3 * h, error) - NSecDer8(f, t - 3 * h, error);

			if (error)
			{
				VectorN<Real, N> y7 = (NSecDer8(f, t + 4 * h, error) - NSecDer8(f, t - 4 * h, error) - 6 * y3 - 14 * y1 - 14 * y2) / 2;

				*error = y7.NormL2() / (140 * h) + 5 * (yh.NormL2() + ymh.NormL2()) * Constants::Epsilon / h;
			}
			return (y3 + 9 * y2 + 45 * y1) / (60 * h);
		}

		/********************************************************************************************************************/
		/********                               Numerical derivatives of EIGHTH order                                ********/
		/********************************************************************************************************************/

		//////////////////////////              RealFunction            //////////////////////////
		static Real NDer8(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = f(x + h);
			Real ymh = f(x - h);
			Real y1 = yh - ymh;
			Real y2 = f(x - 2 * h) - f(x + 2 * h);
			Real y3 = f(x + 3 * h) - f(x - 3 * h);
			Real y4 = f(x - 4 * h) - f(x + 4 * h);

			Real tmp1 = 3 * y4 / 8 + 4 * y3;
			Real tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				// Mathematica code to generate fd scheme for 7th derivative:
				// Sum[(-1)^i*Binomial[9, i]*(f[x+(4-i)*h] + f[x+(5-i)*h])/2, {i, 0, 9}]
				// Mathematica to demonstrate that this is a finite difference formula for 7th derivative:
				// Series[(f[x+5*h]-f[x- 5*h])/2 + 4*(f[x-4*h] - f[x+4*h]) + 27*(f[x+3*h] - f[x-3*h])/2 + 24*(f[x-2*h]  - f[x+2*h]) + 21*(f[x+h] - f[x-h]), {h, 0, 15}]
				Real f9 = (f(x + 5 * h) - f(x - 5 * h)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}
			return (tmp1 + tmp2) / (105 * h);
		}
		static Real NDer8(const IRealFunction& f, Real x, Real* error)
		{
			// Error bound ~eps^8/9.
			// In Real precision, we only expect to lose two digits of precision while using this formula, at the cost of 8 function evaluations.
			// Error: h^8|f^(9)(x)|/630 + 7|f(x)|eps/h assuming 7 unstabilized additions.
			// Mathematica code to get the error:
			// Series[(f[x+h]-f[x-h])*(4/5) + (1/5)*(f[x-2*h] - f[x+2*h]) + (4/105)*(f[x+3*h] - f[x-3*h]) + (1/280)*(f[x-4*h] - f[x+4*h]), {h, 0, 9}]
			// If we used Kahan summation, we could get the max error down to h^8|f^(9)(x)|/630 + |f(x)|eps/h.

			return NDer8(f, x, NDer8_h, error);
		}
		static Real NDer8(const IRealFunction& f, Real x)
		{
			return NDer8(f, x, NDer8_h, nullptr);
		}

		static Real NDer8Left(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer8(f, x - 6 * NDer8_h, NDer8_h, error); }
		static Real NDer8Right(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer8(f, x + 6 * NDer8_h, NDer8_h, error); }
		static Real NDer8Left(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer8(f, x - 6 * h, h, error); }
		static Real NDer8Right(const IRealFunction& f, Real x, Real h, Real* error = nullptr) { return NDer8(f, x + 6 * h, h, error); }

		static Real NSecDer8(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NSecDer8(f, x, NDer8_h, error);
		}
		static Real NSecDer8(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NDer8(f, x + h, error);
			Real ymh = NDer8(f, x - h, error);
			Real y1 = yh - ymh;
			Real y2 = NDer8(f, x - 2 * h, error) - NDer8(f, x + 2 * h, error);
			Real y3 = NDer8(f, x + 3 * h, error) - NDer8(f, x - 3 * h, error);
			Real y4 = NDer8(f, x - 4 * h, error) - NDer8(f, x + 4 * h, error);

			Real tmp1 = 3 * y4 / 8 + 4 * y3;
			Real tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				Real f9 = (NDer8(f, x + 5 * h, error) - NDer8(f, x - 5 * h, error)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}
			return (tmp1 + tmp2) / (105 * h);
		}

		static Real NThirdDer8(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			return NThirdDer8(f, x, NDer8_h, error);
		}
		static Real NThirdDer8(const IRealFunction& f, Real x, Real h, Real* error = nullptr)
		{
			Real yh = NSecDer8(f, x + h, error);
			Real ymh = NSecDer8(f, x - h, error);
			Real y1 = yh - ymh;
			Real y2 = NSecDer8(f, x - 2 * h, error) - NSecDer8(f, x + 2 * h, error);
			Real y3 = NSecDer8(f, x + 3 * h, error) - NSecDer8(f, x - 3 * h, error);
			Real y4 = NSecDer8(f, x - 4 * h, error) - NSecDer8(f, x + 4 * h, error);

			Real tmp1 = 3 * y4 / 8 + 4 * y3;
			Real tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				Real f9 = (NSecDer8(f, x + 5 * h, error) - NSecDer8(f, x - 5 * h, error)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}
			return (tmp1 + tmp2) / (105 * h);
		}

		//////////////////////////             ScalarFunction           //////////////////////////
		template <int N>
		static Real NDer8Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer8Partial(f, deriv_index, point, NDer8_h, error);
		}

		template <int N>
		static Real NDer8Partial(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };

			x[deriv_index] = orig_x + h;
			Real yh = f(x);

			x[deriv_index] = orig_x - h;
			Real ymh = f(x);

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f(x);

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f(x);

			x[deriv_index] = orig_x + 3 * h;
			Real y3h = f(x);

			x[deriv_index] = orig_x - 3 * h;
			Real ym3h = f(x);

			x[deriv_index] = orig_x + 4 * h;
			Real y4h = f(x);

			x[deriv_index] = orig_x - 4 * h;
			Real ym4h = f(x);

			Real y1 = yh - ymh;
			Real y2 = ym2h - y2h;
			Real y3 = y3h - ym3h;
			Real y4 = ym4h - y4h;

			Real tmp1 = 3 * y4 / 8 + 4 * y3;
			Real tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				x[deriv_index] = orig_x + 5 * h;
				Real y5h = f(x);

				x[deriv_index] = orig_x - 5 * h;
				Real ym5h = f(x);

				Real f9 = (y5h - ym5h) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}

			return (tmp1 + tmp2) / (105 * h);
		}

		template <int N>
		static Real NSecDer8Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NSecDer8Partial(f, der_ind1, der_ind2, point, NDer8_h, error);
		}

		template <int N>
		static Real NSecDer8Partial(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[der_ind2];
			auto x_eval_pos = point;

			x_eval_pos[der_ind2] = orig_x + h;
			Real yh = NDer8Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - h;
			Real ymh = NDer8Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x + 2 * h;
			Real y2h = NDer8Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - 2 * h;
			Real ym2h = NDer8Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x + 3 * h;
			Real y3h = NDer8Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - 3 * h;
			Real ym3h = NDer8Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x + 4 * h;
			Real y4h = NDer8Partial(f, der_ind1, x_eval_pos, error);

			x_eval_pos[der_ind2] = orig_x - 4 * h;
			Real ym4h = NDer8Partial(f, der_ind1, x_eval_pos, error);

			Real y1 = yh - ymh;
			Real y2 = ym2h - y2h;
			Real y3 = y3h - ym3h;
			Real y4 = ym4h - y4h;

			Real tmp1 = 3 * y4 / 8 + 4 * y3;
			Real tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				x_eval_pos[der_ind2] = orig_x + 5 * h;
				Real y5h = NDer8Partial(f, der_ind1, x_eval_pos, error);

				x_eval_pos[der_ind2] = orig_x - 5 * h;
				Real ym5h = NDer8Partial(f, der_ind1, x_eval_pos, error);

				Real f9 = (y5h - ym5h) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}

			return (tmp1 + tmp2) / (105 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer8PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer8PartialByAll(f, point, NDer8_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer8PartialByAll(const IScalarFunction<N>& f, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer8Partial(f, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer8Partial(f, i, point, h);
			}

			return ret;
		}

		//////////////////////////             VectorFunction           //////////////////////////
		template <int N>
		static Real NDer8Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real* error = nullptr)
		{
			return NDer8Partial(f, func_index, deriv_index, point, NDer8_h, error);
		}

		template <int N>
		static Real NDer8Partial(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real h, Real* error = nullptr)
		{
			Real     orig_x = point[deriv_index];

			VectorN<Real, N> x{ point };

			x[deriv_index] = orig_x + h;
			Real yh = f(x)[func_index];

			x[deriv_index] = orig_x - h;
			Real ymh = f(x)[func_index];

			x[deriv_index] = orig_x + 2 * h;
			Real y2h = f(x)[func_index];

			x[deriv_index] = orig_x - 2 * h;
			Real ym2h = f(x)[func_index];

			x[deriv_index] = orig_x + 3 * h;
			Real y3h = f(x)[func_index];

			x[deriv_index] = orig_x - 3 * h;
			Real ym3h = f(x)[func_index];

			x[deriv_index] = orig_x + 4 * h;
			Real y4h = f(x)[func_index];

			x[deriv_index] = orig_x - 4 * h;
			Real ym4h = f(x)[func_index];

			Real y1 = yh - ymh;
			Real y2 = ym2h - y2h;
			Real y3 = y3h - ym3h;
			Real y4 = ym4h - y4h;

			Real tmp1 = 3 * y4 / 8 + 4 * y3;
			Real tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				x[deriv_index] = orig_x + 5 * h;
				Real y5h = f(x)[func_index];

				x[deriv_index] = orig_x - 5 * h;
				Real ym5h = f(x)[func_index];

				Real f9 = (y5h - ym5h) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh)) * Constants::Epsilon / h;
			}

			return (tmp1 + tmp2) / (105 * h);
		}

		template <int N>
		static VectorN<Real, N> NDer8PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, VectorN<Real, N>* error = nullptr)
		{
			return NDer8PartialByAll(f, func_index, point, NDer8_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer8PartialByAll(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, Real h, VectorN<Real, N>* error = nullptr)
		{
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++)
			{
				if (error)
					ret[i] = NDer8Partial(f, func_index, i, point, h, &(*error)[i]);
				else
					ret[i] = NDer8Partial(f, func_index, i, point, h);
			}

			return ret;
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer8PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, MatrixNM<Real, N, N>* error = nullptr)
		{
			return NDer8PartialAllByAll(f, point, NDer8_h, error);
		}

		template <int N>
		static MatrixNM<Real, N, N> NDer8PartialAllByAll(const IVectorFunction<N>& f, const VectorN<Real, N>& point, Real h, MatrixNM<Real, N, N>* error = nullptr)
		{
			MatrixNM<Real, N, N> ret;

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					if (error)
						ret(i, j) = NDer8Partial(f, i, j, point, h, &((*error)(i, j)));
					else
						ret(i, j) = NDer8Partial(f, i, j, point, h);
				}

			return ret;
		}

		/////////////////////////             ParametricCurve           /////////////////////////
		template <int N>
		static VectorN<Real, N> NDer8(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NDer8(f, t, NDer8_h, error);
		}

		template <int N>
		static VectorN<Real, N> NDer8(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = f(t + h);
			VectorN<Real, N> ymh = f(t - h);
			VectorN<Real, N> y1 = yh - ymh;
			VectorN<Real, N> y2 = f(t - 2 * h) - f(t + 2 * h);
			VectorN<Real, N> y3 = f(t + 3 * h) - f(t - 3 * h);
			VectorN<Real, N> y4 = f(t - 4 * h) - f(t + 4 * h);

			VectorN<Real, N> tmp1 = 3 * y4 / 8 + 4 * y3;
			VectorN<Real, N> tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				VectorN<Real, N> f9 = (f(t + 5 * h) - f(t - 5 * h)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = f9.NormL2() / (630 * h) + 7 * (yh.NormL2() + ymh.NormL2()) * Constants::Epsilon / h;
			}
			return (tmp1 + tmp2) / (105 * h);
		}

		template <int N>
		static VectorN<Real, N> NSecDer8(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NSecDer8(f, t, NDer8_h, error);
		}

		template <int N>
		static VectorN<Real, N> NSecDer8(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = NDer8(f, t + h, error);
			VectorN<Real, N> ymh = NDer8(f, t - h, error);
			VectorN<Real, N> y1 = yh - ymh;
			VectorN<Real, N> y2 = NDer8(f, t - 2 * h, error) - NDer8(f, t + 2 * h, error);
			VectorN<Real, N> y3 = NDer8(f, t + 3 * h, error) - NDer8(f, t - 3 * h, error);
			VectorN<Real, N> y4 = NDer8(f, t - 4 * h, error) - NDer8(f, t + 4 * h, error);

			VectorN<Real, N> tmp1 = 3 * y4 / 8 + 4 * y3;
			VectorN<Real, N> tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				VectorN<Real, N> f9 = (NDer8(f, t + 5 * h, error) - NDer8(f, t - 5 * h, error)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = f9.NormL2() / (630 * h) + 7 * (yh.NormL2() + ymh.NormL2()) * Constants::Epsilon / h;
			}
			return (tmp1 + tmp2) / (105 * h);
		}

		template <int N>
		static VectorN<Real, N> NThirdDer8(const IParametricCurve<N>& f, Real t, Real* error = nullptr)
		{
			return NThirdDer8(f, t, NDer8_h, error);
		}

		template <int N>
		static VectorN<Real, N> NThirdDer8(const IParametricCurve<N>& f, Real t, Real h, Real* error = nullptr)
		{
			VectorN<Real, N> yh = NSecDer8(f, t + h, error);
			VectorN<Real, N> ymh = NSecDer8(f, t - h, error);
			VectorN<Real, N> y1 = yh - ymh;
			VectorN<Real, N> y2 = NSecDer8(f, t - 2 * h, error) - NSecDer8(f, t + 2 * h, error);
			VectorN<Real, N> y3 = NSecDer8(f, t + 3 * h, error) - NSecDer8(f, t - 3 * h, error);
			VectorN<Real, N> y4 = NSecDer8(f, t - 4 * h, error) - NSecDer8(f, t + 4 * h, error);

			VectorN<Real, N> tmp1 = 3 * y4 / 8 + 4 * y3;
			VectorN<Real, N> tmp2 = 21 * y2 + 84 * y1;

			if (error)
			{
				VectorN<Real, N> f9 = (NSecDer8(f, t + 5 * h, error) - NSecDer8(f, t - 5 * h, error)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;

				*error = f9.NormL2() / (630 * h) + 7 * (yh.NormL2() + ymh.NormL2()) * Constants::Epsilon / h;
			}
			return (tmp1 + tmp2) / (105 * h);
		}


		/********************************************************************************************************************/
		/********                            Definitions of default derivation functions                             ********/
		/********************************************************************************************************************/
		static inline Real(*Derive)(const IRealFunction& f, Real x) = Derivation::NDer4;
		static inline Real(*DeriveErr)(const IRealFunction& f, Real x, Real* error) = Derivation::NDer4;

		static inline Real(*DeriveSec)(const IRealFunction& f, Real x, Real* error) = Derivation::NSecDer4;

		static inline Real(*DeriveThird)(const IRealFunction& f, Real x, Real* error) = Derivation::NThirdDer2;


		template<int N>
		static inline Real(*DerivePartial)(const IScalarFunction<N>& f, int deriv_index, const VectorN<Real, N>& point, Real* error) = Derivation::NDer4Partial;

		template<int N>
		static inline Real(*DeriveSecPartial)(const IScalarFunction<N>& f, int der_ind1, int der_ind2, const VectorN<Real, N>& point, Real* error) = Derivation::NSecDer4Partial;

		template<int N>
		static inline VectorN<Real, N>(*DerivePartialAll)(const IScalarFunction<N>& f, const VectorN<Real, N>& point, VectorN<Real, N>* error) = Derivation::NDer4PartialByAll;


		template<int N>
		static inline Real(*DeriveVecPartial)(const IVectorFunction<N>& f, int func_index, int deriv_index, const VectorN<Real, N>& point, Real* error) = Derivation::NDer4Partial;

		template<int N>
		static inline VectorN<Real, N>(*DeriveVecPartialAll)(const IVectorFunction<N>& f, int func_index, const VectorN<Real, N>& point, VectorN<Real, N>* error) = Derivation::NDer4PartialByAll;

		template<int N>
		static inline MatrixNM<Real, N, N>(*DeriveVecPartialAllByAll)(const IVectorFunction<N>& f, const VectorN<Real, N>& point, MatrixNM<Real, N, N>* error) = Derivation::NDer4PartialAllByAll;


		template<int N>
		static inline VectorN<Real, N>(*DeriveCurve)(const IParametricCurve<N>& f, Real x, Real* error) = Derivation::NDer4;

		template<int N>
		static inline VectorN<Real, N>(*DeriveCurveSec)(const IParametricCurve<N>& f, Real x, Real* error) = Derivation::NSecDer4;

		template<int N>
		static inline VectorN<Real, N>(*DeriveCurveThird)(const IParametricCurve<N>& f, Real x, Real* error) = Derivation::NThirdDer4;
	};


	template<int N, int M>
	class Jacobian
	{
	public:
		static MatrixNM<Real, N, N> calc(const IVectorFunction<N>& func, const VectorN<Real, N>& pos)
		{
			MatrixNM<Real, N, N> jac;

			for (int i = 0; i < N; ++i)
				for (int j = 0; j < N; ++j)
				{
					jac(i, j) = Derivation::NDer4Partial(func, i, j, pos);
				}

			return jac;
		}

		static MatrixNM<Real, M, N> calc(const IVectorFunctionNM<N, M>& func, const VectorN<Real, N>& pos)
		{
			MatrixNM<Real, M, N> jac;

			for (int i = 0; i < M; ++i)
				for (int j = 0; j < N; ++j)
				{
					jac(i, j) = Derivation::NDer4Partial(func, i, j, pos);
				}

			return jac;
		}
	};
}

///////////////////////////   ./include/core/Integration.h   ///////////////////////////



namespace MML
{
	enum IntegrationMethod { TRAP, SIMPSON, ROMBERG, GAUSS10 };


	static Real TrapRefine(const IRealFunction& func, const Real a, const Real b, const int n)
	{
		// This routine computes the nth stage of refinement of an extended trapezoidal rule. func is input
		// as a pointer to the function to be integrated between limits a and b, also input. When called with
		// n=1, the routine returns the crudest estimate of Rab f(x)dx. Subsequent calls with n=2,3,...
		// (in that sequential order) will improve the accuracy by adding 2n-2 additional interior points.
		Real x, tnm, sum, del;
		static Real s;
		int it, j;

		if (n == 1) {
			return (s = 0.5 * (b - a) * (func(a) + func(b)));
		}
		else
		{
			for (it = 1, j = 1; j < n - 1; j++)
				it <<= 1;

			tnm = it;
			del = (b - a) / tnm;
			x = a + 0.5 * del;

			for (sum = 0.0, j = 0; j < it; j++, x += del)
				sum += func(x);

			s = 0.5 * (s + (b - a) * sum / tnm);

			return s;
		}
	}

	static Real IntegrateTrap(const IRealFunction& func, const Real a, const Real b, Real req_eps, Real* achieved_precision)
	{
		// Returns the integral of the function func from a to b. The parameters EPS can be set to the
		// desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum allowed
		// number of steps. Integration is performed by the trapezoidal rule.

		// Unsophisticated as it is, routine qtrap is in fact a fairly robust way of doing
		// integrals of functions that are not very smooth. Increased sophistication will usually
		// translate into a higher-order method whose efficiency will be greater only for
		// sufficiently smooth integrands. qtrap is the method of choice, e.g., for an integrand
		// which is a function of a variable that is linearly interpolated between measured data
		// points. Be sure that you do not require too stringent an EPS, however: If qtrap takes
		// too many steps in trying to achieve your required accuracy, accumulated roundoff
		// errors may start increasing, and the routine may never converge. 
		// Value 1e-6 is just on the edge of trouble for most 32-bit machines; it is achievable when the
		// convergence is moderately rapid, but not otherwise.
		int j;
		Real s, olds = 0.0;
		Real diff = 0.0, threshold = 0.0;

		for (j = 0; j < Defaults::IntegrateTrapMaxSteps; j++)
		{
			s = TrapRefine(func, a, b, j + 1);

			if (j > 5)
			{
				diff = s - olds;
				threshold = req_eps * std::abs(olds);
				if (std::abs(diff) < threshold || (s == 0.0 && olds == 0.0))
				{
					// std::cout << "\ns : " << s << " olds : " << olds <<  " diff : " << diff << " threshold : " << threshold << std::endl;
					if (achieved_precision != nullptr)
						*achieved_precision = std::abs(diff);

					return s;
				}
			}
			olds = s;
		}
		if (achieved_precision != nullptr)
			*achieved_precision = std::abs(diff);

		return s;
	}
	static Real IntegrateTrap(const IRealFunction& func, const Real a, const Real b)
	{
		return IntegrateTrap(func, a, b, Defaults::IntegrateTrapEPS, nullptr);
	}


	static Real IntegrateGauss10(const IRealFunction& func, const Real a, const Real b)
	{
		// Returns the integral of the function func between a and b, by ten-point GaussLegendre integration: 
		// the function is evaluated exactly ten times at interior points in the range of integration.        
		static const Real x[] = { 0.1488743389816312,0.4333953941292472,
														0.6794095682990244,0.8650633666889845,0.9739065285171717 };
		static const Real w[] = { 0.2955242247147529,0.2692667193099963,
														0.2190863625159821,0.1494513491505806,0.0666713443086881 };
		Real xm = 0.5 * (b + a);
		Real xr = 0.5 * (b - a);
		Real s = 0;
		for (int j = 0; j < 5; j++) {
			Real dx = xr * x[j];
			s += w[j] * (func(xm + dx) + func(xm - dx));
		}
		return s *= xr;
	}

	struct SurfIntf2 : public IRealFunction
	{
		mutable Real xsav;
		IScalarFunction<2>& funcToInt;

		SurfIntf2(IScalarFunction<2>& func) : funcToInt(func) {}
		Real operator()(const Real y) const
		{
			VectorN<Real, 2> v{ xsav,y };
			return funcToInt(v);
		}
	};

	struct SurfIntf1 : public IRealFunction
	{
		mutable SurfIntf2 f2;
		IntegrationMethod method;

		IScalarFunction<2>& funcToInt;
		Real(*y1)(Real);
		Real(*y2)(Real);

		SurfIntf1(IScalarFunction<2>& func, IntegrationMethod inMethod, Real yy1(Real), Real yy2(Real)) : y1(yy1), y2(yy2), f2(func), funcToInt(func), method(inMethod)
		{}

		Real operator()(const Real x) const
		{
			f2.xsav = x;
			switch (method)
			{
				//case SIMPSON:
				//	return IntegrateSimpson(f2, y1(x), y2(x));
				//case ROMBERG:
				//	return IntegrateRomberg(f2, y1(x), y2(x));
			case GAUSS10:
				return IntegrateGauss10(f2, y1(x), y2(x));
			default:
				Real ret = IntegrateTrap(f2, y1(x), y2(x));
				//std::cout << "IntegrateTrap: " << ret << std::endl;
				return ret;
			}
		}
	};

	static Real IntegrateSurface(IScalarFunction<2>& func, IntegrationMethod method, const Real x1, const Real x2, Real y1(Real), Real y2(Real))
	{
		SurfIntf1 f1(func, method, y1, y2);

		switch (method)
		{
			//case SIMPSON:
			//	return IntegrateSimpson(f1, x1, x2);
			//case ROMBERG:
			//	return IntegrateRomberg(f1, x1, x2);
		case GAUSS10:
			return IntegrateGauss10(f1, x1, x2);
		default:
			return IntegrateTrap(f1, x1, x2);
		}
	}

	struct VolIntf3 : public IRealFunction
	{
		mutable Real xsav, ysav;
		IScalarFunction<3>& funcToInt;

		VolIntf3(IScalarFunction<3>& func) : funcToInt(func), xsav{ 0 }, ysav{ 0 } {}
		Real operator()(const Real z) const
		{
			VectorN<Real, 3> v{ xsav,ysav,z };
			return funcToInt(v);
		}
	};
	struct VolIntf2 : public IRealFunction
	{
		mutable VolIntf3 f3;

		IScalarFunction<3>& funcToInt;
		Real(*z1)(Real, Real);
		Real(*z2)(Real, Real);

		VolIntf2(IScalarFunction<3>& func, Real zz1(Real, Real), Real zz2(Real, Real)) : z1(zz1), z2(zz2), funcToInt(func), f3(func) {}

		Real operator()(const Real y) const
		{
			f3.ysav = y;
			return IntegrateGauss10(f3, z1(f3.xsav, y), z2(f3.xsav, y));
		}
	};
	struct VolIntf1 : public IRealFunction
	{
		mutable VolIntf2 f2;

		IScalarFunction<3>& funcToInt;
		Real(*y1)(Real);
		Real(*y2)(Real);

		VolIntf1(IScalarFunction<3>& func, Real yy1(Real), Real yy2(Real), Real z1(Real, Real),
			Real z2(Real, Real)) : y1(yy1), y2(yy2), f2(func, z1, z2), funcToInt(func)
		{}

		Real operator()(const Real x) const
		{
			f2.f3.xsav = x;
			return IntegrateGauss10(f2, y1(x), y2(x));
		}
	};

	static Real IntegrateVolume(IScalarFunction<3>& func, const Real x1, const Real x2, Real y1(Real), Real y2(Real),
		Real z1(Real, Real), Real z2(Real, Real))
	{
		VolIntf1 f1(func, y1, y2, z1, z2);

		return IntegrateGauss10(f1, x1, x2);
	}

	static inline Real(*Integrate)(const MML::IRealFunction& f, Real a, Real b, Real req_eps, Real* achieved_precision) = IntegrateTrap;

} // end namespace
///////////////////////////   ./include/core/Function.h   ///////////////////////////




namespace MML
{
	///////////////////////////     REAL FUNCTION      ////////////////////////////////////
	class RealFunction : public IRealFunction
	{
		Real(*_func)(const Real);
	public:
		RealFunction(Real(*inFunc)(const Real)) : _func(inFunc) {}

		Real operator()(const Real x) const { return _func(x); }
	};
	class RealFunctionFromStdFunc : public IRealFunction
	{
		std::function<Real(const Real)> _func;
	public:
		RealFunctionFromStdFunc(std::function<Real(const Real)> inFunc) : _func(inFunc) {}

		Real operator()(const Real x) const { return _func(x); }
	};

	///////////////////////////     SCALAR FUNCTION       //////////////////////////////////
	template<int N>
	class ScalarFunction : public IScalarFunction<N>
	{
		Real(*_func)(const VectorN<Real, N>&);
	public:
		ScalarFunction(Real(*inFunc)(const VectorN<Real, N>&)) : _func(inFunc) {}

		Real operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	template<int N>
	class ScalarFunctionFromStdFunc : public IScalarFunction<N>
	{
		std::function<Real(const VectorN<Real, N>&)> _func;
	public:
		ScalarFunctionFromStdFunc(std::function<Real(const VectorN<Real, N>&)> inFunc) : _func(inFunc) {}

		Real operator()(const VectorN<Real, N>& x) const { return _func(x); }
	};

	/////////////////////////    VECTOR FUNCTION N -> N      ///////////////////////////////////
	template<int N>
	class VectorFunction : public IVectorFunction<N>
	{
		VectorN<Real, N>(*_func)(const VectorN<Real, N>&);
	public:
		VectorFunction(VectorN<Real, N>(*inFunc)(const VectorN<Real, N>&)) : _func(inFunc) {}

		VectorN<Real, N>     operator()(const VectorN<Real, N>& x) const { return _func(x); }

		MatrixNM<Real, N, N> jacobian(const VectorN<Real, N>& x) const { return Jacobian<N, N>::calc(*this, x); }
	};

	template<int N>
	class VectorFunctionFromStdFunc : public IVectorFunction<N>
	{
		std::function<VectorN<Real, N>(const VectorN<Real, N>&)> _func;
	public:
		VectorFunctionFromStdFunc(std::function<VectorN<Real, N>(const VectorN<Real, N>&)>& inFunc) : _func(inFunc) {}

		VectorN<Real, N>     operator()(const VectorN<Real, N>& x) const { return _func(x); }

		MatrixNM<Real, N, N> jacobian(const VectorN<Real, N>& x) const { return Jacobian<N, N>::calc(*this, x); }
	};

	/////////////////////////    VECTOR FUNCTION M -> N      ///////////////////////////////////
	template<int N, int M>
	class VectorFunctionNM : public IVectorFunctionNM<N, M>
	{
		VectorN<Real, M>(*_func)(const VectorN<Real, N>&);
	public:
		VectorFunctionNM(VectorN<Real, M>(*inFunc)(const VectorN<Real, N>&)) : _func(inFunc) {}

		VectorN<Real, M>     operator()(const VectorN<Real, N>& x) const { return _func(x); }

		MatrixNM<Real, M, N> jacobian(const VectorN<Real, N>& x) const { return Jacobian<N, M>::calc(*this, x); }
	};

	template<int N, int M>
	class VectorFunctionNMFromStdFunc : public IVectorFunctionNM<N, M>
	{
		std::function<VectorN<Real, M>(const VectorN<Real, N>&)> _func;
	public:
		VectorFunctionNMFromStdFunc(std::function<VectorN<Real, M>(const VectorN<Real, N>&)>& inFunc) : _func(inFunc) {}

		VectorN<Real, M>     operator()(const VectorN<Real, N>& x) const { return _func(x); }

		MatrixNM<Real, M, N> jacobian(const VectorN<Real, N>& x) const { return Jacobian<N, M>::calc(*this, x); }
	};

	//////////////////////     PARAMETRIC CURVE             ///////////////////////////////////
	template<int N>
	class   ParametricCurve : public IParametricCurve<N>
	{
		Real _minT;
		Real _maxT;
		VectorN<Real, N>(*_func)(Real);
	public:
		ParametricCurve(VectorN<Real, N>(*inFunc)(Real)) : _func(inFunc), _minT(Constants::NegativeInf), _maxT(Constants::PositiveInf) {}
		ParametricCurve(Real minT, Real maxT, VectorN<Real, N>(*inFunc)(Real)) : _func(inFunc), _minT(minT), _maxT(maxT) {}

		Real getMinT() const { return _minT; }
		Real getMaxT() const { return _maxT; }

		virtual VectorN<Real, N> operator()(Real x) const { return _func(x); }
	};

	template<int N>
	class ParametricCurveFromStdFunc : public IParametricCurve<N>
	{
		Real _minT;
		Real _maxT;
		std::function<VectorN<Real, N>(Real)> _func;
	public:
		ParametricCurveFromStdFunc(std::function<VectorN<Real, N>(Real)>& inFunc) : _func(inFunc), _minT(Constants::NegativeInf), _maxT(Constants::PositiveInf) {}
		ParametricCurveFromStdFunc(Real minT, Real maxT, std::function<VectorN<Real, N>(Real)>& inFunc) : _func(inFunc), _minT(minT), _maxT(maxT) {}

		Real getMinT() const { return _minT; }
		Real getMaxT() const { return _maxT; }

		VectorN<Real, N> operator()(Real x) const { return _func(x); }
	};

	/////////////////////       PARAMETRIC SURFACE         //////////////////////////////////
	template<int N>
	class ParametricSurface : public IParametricSurface<N>
	{
		// TODO - ensure that N is at least 3!!!
		Real _minX;
		Real _maxX;
		Real _minY;
		Real _maxY;
		VectorN<Real, N>(*_func)(Real u, Real w);

	public:
		ParametricSurface(VectorN<Real, N>(*inFunc)(Real u, Real w)) : _func(inFunc), _minX(Constants::NegativeInf), _maxX(Constants::PositiveInf), _minY(Constants::NegativeInf), _maxY(Constants::PositiveInf) {}
		ParametricSurface(VectorN<Real, N>(*inFunc)(Real u, Real w), Real minX, Real maxX, Real minY, Real maxY) : _func(inFunc), _minX(minX), _maxX(maxX), _minY(minY), _maxY(maxY) {}

		VectorN<Real, N> operator()(Real u, Real w) const { return _func(u, w); }

		virtual Real getMinU() const { return _minX; }
		virtual Real getMaxU() const { return _maxX; }
		virtual Real getMinW() const { return _minY; }
		virtual Real getMaxW() const { return _maxY; }
	};

	// imati cemo i surface Discrete, kreiran od triangles?

	template<int N>
	class ParametricSurfaceFromStdFunc : public IParametricSurface<N>
	{
		Real _minX;
		Real _maxX;
		Real _minY;
		Real _maxY;
		std::function<VectorN<Real, N>(Real u, Real w)> _func;
	public:
		ParametricSurfaceFromStdFunc(std::function<VectorN<Real, N>(Real u, Real w)>& inFunc) : _func(inFunc), _minX(Constants::NegativeInf), _maxX(Constants::PositiveInf), _minY(Constants::NegativeInf), _maxY(Constants::PositiveInf) {}
		ParametricSurfaceFromStdFunc(std::function<VectorN<Real, N>(Real u, Real w)>& inFunc, Real minX, Real maxX, Real minY, Real maxY) : _func(inFunc), _minX(minX), _maxX(maxX), _minY(minY), _maxY(maxY) {}

		VectorN<Real, N> operator()(Real u, Real w) const { return _func(u, w); }

		virtual Real getMinU() const { return _minX; }
		virtual Real getMaxU() const { return _maxX; }
		virtual Real getMinW() const { return _minY; }
		virtual Real getMaxW() const { return _maxY; }
	};
} // end namespace

///////////////////////////   ./include/core/FunctionHelpers.h   ///////////////////////////




namespace MML
{

	static RealPolynom TaylorSeries2(IRealFunction& f, Real a)
	{
		RealPolynom ret(2);

		Real val = f(a);
		Real coef1 = Derivation::NDer6(f, a);
		Real coef2 = Derivation::NSecDer4(f, a) / 2.0;

		ret[0] = val - coef1 * a + coef2 * POW2(a);
		ret[1] = coef1 - 2.0 * coef2 * a;
		ret[2] = coef2;

		return ret;
	}
	static RealPolynom TaylorSeries3(IRealFunction& f, Real a)
	{
		RealPolynom ret(3);

		Real val = f(a);
		Real coef1 = Derivation::NDer6(f, a);
		Real coef2 = Derivation::NSecDer4(f, a) / 2.0;
		Real coef3 = Derivation::NThirdDer2(f, a) / 6.0;

		ret[0] = val - coef1 * a + coef2 * POW2(a) - coef3 * pow(a, 3);
		ret[1] = coef1 - 2.0 * coef2 * a + 3.0 * coef3 * POW2(a);
		ret[2] = coef2 - 3.0 * coef3 * a;
		ret[3] = coef3;

		return ret;
	}

	/////////////////////       FUNCTION HELPERS         //////////////////////////////////
	class RealFuncDerived1 : public IRealFunction
	{
		IRealFunction& _f1;
		Real _deriv_step;
	public:
		RealFuncDerived1(IRealFunction& f1) : _f1(f1), _deriv_step(0.0) {}
		RealFuncDerived1(IRealFunction& f1, Real deriv_step) : _f1(f1), _deriv_step(deriv_step) {}

		Real operator()(Real x) const {
			if (_deriv_step != 0.0)
				return Derivation::NDer1(_f1, x, _deriv_step);
			else
				return Derivation::NDer1(_f1, x);
		}
	};
	class RealFuncDerived2 : public IRealFunction
	{
		IRealFunction& _f1;
		Real _deriv_step;
	public:
		RealFuncDerived2(IRealFunction& f1) : _f1(f1), _deriv_step(0.0) {}
		RealFuncDerived2(IRealFunction& f1, Real deriv_step) : _f1(f1), _deriv_step(deriv_step) {}

		Real operator()(Real x) const {
			if (_deriv_step != 0.0)
				return Derivation::NDer2(_f1, x, _deriv_step);
			else
				return Derivation::NDer2(_f1, x);
		}
	};
	class RealFuncDerived4 : public IRealFunction
	{
		IRealFunction& _f1;
		Real _deriv_step;
	public:
		RealFuncDerived4(IRealFunction& f1) : _f1(f1), _deriv_step(0.0) {}
		RealFuncDerived4(IRealFunction& f1, Real deriv_step) : _f1(f1), _deriv_step(deriv_step) {}

		Real operator()(Real x) const {
			if (_deriv_step != 0.0)
				return Derivation::NDer4(_f1, x, _deriv_step);
			else
				return Derivation::NDer4(_f1, x);
		}
	};
	class RealFuncDerived6 : public IRealFunction
	{
		IRealFunction& _f1;
		Real _deriv_step;
	public:
		RealFuncDerived6(IRealFunction& f1) : _f1(f1), _deriv_step(0.0) {}
		RealFuncDerived6(IRealFunction& f1, Real deriv_step) : _f1(f1), _deriv_step(deriv_step) {}

		Real operator()(Real x) const {
			if (_deriv_step != 0.0)
				return Derivation::NDer6(_f1, x, _deriv_step);
			else
				return Derivation::NDer6(_f1, x);
		}
	};
	class RealFuncDerived8 : public IRealFunction
	{
		IRealFunction& _f1;
		Real _deriv_step;
	public:
		RealFuncDerived8(IRealFunction& f1) : _f1(f1), _deriv_step(0.0) {}
		RealFuncDerived8(IRealFunction& f1, Real deriv_step) : _f1(f1), _deriv_step(deriv_step) {}

		Real operator()(Real x) const {
			if (_deriv_step != 0.0)
				return Derivation::NDer8(_f1, x, _deriv_step);
			else
				return Derivation::NDer8(_f1, x);
		}
	};

	// TODO - integrate helper - ima i x1, a x2 je parametar funkcije

	class RealFuncDiffHelper : public IRealFunction
	{
		IRealFunction& _f1, & _f2;
	public:
		RealFuncDiffHelper(IRealFunction& f1, IRealFunction& f2) : _f1(f1), _f2(f2) {}
		Real operator()(Real x) const { return _f1(x) - _f2(x); }
	};

	class RealFuncDiffAbsHelper : public IRealFunction
	{
		IRealFunction& _f1, & _f2;
	public:
		RealFuncDiffAbsHelper(IRealFunction& f1, IRealFunction& f2) : _f1(f1), _f2(f2) {}
		Real operator()(Real x) const { return std::abs(_f1(x) - _f2(x)); }
	};

	class RealFuncDiffSqrHelper : public IRealFunction
	{
		IRealFunction& _f1, & _f2;
	public:
		RealFuncDiffSqrHelper(IRealFunction& f1, IRealFunction& f2) : _f1(f1), _f2(f2) {}
		Real operator()(Real x) const { return POW2(_f1(x) - _f2(x)); }
	};

} // end namespace

///////////////////////////   ./include/core/InterpolatedFunction.h   ///////////////////////////




namespace MML
{
	class RealFunctionInterpolated : public IRealFunction
	{
		double _xmin, _xmax;

	public:
		Real virtual rawinterp(int jlo, Real x) const = 0;
	};

	class RealFunctionInterpolatedBase : public RealFunctionInterpolated
	{
		// TODO - HIGH, LAKO, sve privatizirati!
	public:
		mutable int jsav, cor;

		int n, mm, dj;
		Vector<Real> xx, yy;

		RealFunctionInterpolatedBase(const Vector<Real>& x, const Vector<Real>& y, int m)
			: n((int)x.size()), mm(m), jsav(0), cor(0), xx(x), yy(y)
		{
			dj = std::min(1, (int)pow((Real)n, 0.25));
		}

		Real operator()(Real x) const
		{
			int jlo = cor ? hunt(x) : locate(x);
			return rawinterp(jlo, x);
		}

		// Given a value x, return a value j such that x is (insofar as possible) centered in the subrange
		// xx[j..j+mm-1], where xx is the stored pointer. The values in xx must be monotonic, either
		// increasing or decreasing. The returned value is not less than 0, nor greater than n-1.
		int locate(const Real x) const
		{
			int ju, jm, jl;
			if (n < 2 || mm < 2 || mm > n) throw("locate size error");
			bool ascnd = (xx[n - 1] >= xx[0]);
			jl = 0;
			ju = n - 1;
			while (ju - jl > 1) {
				jm = (ju + jl) >> 1;
				if (x >= xx[jm] == ascnd)
					jl = jm;
				else
					ju = jm;
			}
			cor = std::abs(jl - jsav) > dj ? 0 : 1;
			jsav = jl;
			return std::max(0, std::min(n - mm, jl - ((mm - 2) >> 1)));
		}

		// Given a value x, return a value j such that x is (insofar as possible) centered in the subrange
		// xx[j..j+mm-1], where xx is the stored pointer. The values in xx must be monotonic, either
		// increasing or decreasing. The returned value is not less than 0, nor greater than n-1.
		int hunt(const Real x) const
		{
			int jl = jsav, jm, ju, inc = 1;
			if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
			bool ascnd = (xx[n - 1] >= xx[0]);
			if (jl < 0 || jl > n - 1) {
				jl = 0;
				ju = n - 1;
			}
			else {
				if (x >= xx[jl] == ascnd) {
					for (;;) {
						ju = jl + inc;
						if (ju >= n - 1) { ju = n - 1; break; }
						else if (x < xx[ju] == ascnd) break;
						else {
							jl = ju;
							inc += inc;
						}
					}
				}
				else {
					ju = jl;
					for (;;) {
						jl = jl - inc;
						if (jl <= 0) { jl = 0; break; }
						else if (x >= xx[jl] == ascnd) break;
						else {
							ju = jl;
							inc += inc;
						}
					}
				}
			}
			while (ju - jl > 1) {
				jm = (ju + jl) >> 1;
				if (x >= xx[jm] == ascnd)
					jl = jm;
				else
					ju = jm;
			}
			cor = std::abs(jl - jsav) > dj ? 0 : 1;
			jsav = jl;
			return std::max(0, std::min(n - mm, jl - ((mm - 2) >> 1)));
		}
	};

	struct LinearInterpRealFunc : RealFunctionInterpolatedBase
	{
		LinearInterpRealFunc(Vector<Real>& xv, Vector<Real>& yv) : RealFunctionInterpolatedBase(xv, yv, 2) {}

		Real rawinterp(int j, Real x) const {
			if (xx[j] == xx[j + 1]) return yy[j];
			else return yy[j] + ((x - xx[j]) / (xx[j + 1] - xx[j])) * (yy[j + 1] - yy[j]);
		}
	};


	template<int N>
	class ParametricCurveInterpolated : public IParametricCurve<N>
	{
		// TODO 0.9 - HIGH, verify
		Real _minT;
		Real _maxT;
		std::shared_ptr<Vector<Real>> _xvals;
		std::shared_ptr<Matrix<Real>> _yvals;

	public:
		ParametricCurveInterpolated() {}
		ParametricCurveInterpolated(std::shared_ptr<Vector<Real>> xsave, std::shared_ptr<Matrix<Real>> ysave) : _xvals(xsave), _yvals(ysave) {}
		ParametricCurveInterpolated(int inPoints, Real* xsave, Real* ysave)
		{
			_minT = xsave[0];
			_maxT = xsave[inPoints - 1];

			_xvals = std::make_shared<Vector<Real>>(inPoints);
			_yvals = std::make_shared<Matrix<Real>>(N, inPoints);

			for (int i = 0; i < inPoints; i++)
			{
				(*_xvals)[i] = xsave[i];
				for (int j = 0; j < N; j++)
					(*_yvals)(j, i) = ysave[i * N + j];
			}
		}
		ParametricCurveInterpolated(const Vector<Real>& xsave, const Matrix<Real>& ysave)
		{
			_minT = xsave[0];
			_maxT = xsave[xsave.size() - 1];

			_xvals = std::make_shared<Vector<Real>>(xsave);
			_yvals = std::make_shared<Matrix<Real>>(ysave);
		}

		Real getMinT() const { return _minT; }
		Real getMaxT() const { return _maxT; }

		virtual VectorN<Real, N> operator()(Real x) const { return VectorN<Real, N>{0}; }
	};

	template<int N>
	class InterpolatedSurface : public IParametricSurface<N>
	{
	public:
		InterpolatedSurface() {}

		VectorN<Real, N> operator()(const VectorN<Real, 2>& x) const { return VectorN<Real, N>{}; }
	};
}

///////////////////////////   ./include/core/MetricTensor.h   ///////////////////////////



namespace MML
{
	template<int N>
	class MetricTensorField : public ITensorField2<N>
	{
	public:
		MetricTensorField() : ITensorField2<N>(2, 0) { }
		MetricTensorField(int numContra, int numCo) : ITensorField2<N>(numContra, numCo) { }

		Tensor2<N>   operator()(const VectorN<Real, N>& pos) const
		{
			Tensor2<N> ret(this->getNumContravar(), this->getNumCovar());

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					ret(i, j) = this->Component(i, j, pos);

			return ret;
		}

		Real GetChristoffelSymbolFirstKind(int i, int j, int k, const VectorN<Real, N>& pos) const
		{
			MetricTensorField<N>& g = *this;

			Real gamma_ijk = 0.0;
			for (int m = 0; m < N; m++)
			{
				gamma_ijk += g.Component(m, k, pos) * GetChristoffelSymbolSecondKind(m, i, j, pos);
			}
			return gamma_ijk;
		}

		Real GetChristoffelSymbolSecondKind(int i, int j, int k, const VectorN<Real, N>& pos) const
		{
			MetricTensorField<N>& g = *this;

			// depends on covar-contravar of metric tensor!
			Real gamma_ijk = 0.0;
			for (int l = 0; l < N; l++)
			{
				Real coef1 = Derivation::DerivePartial<N>(g, i, l, j, pos, nullptr);
				Real coef2 = Derivation::DerivePartial<N>(g, j, l, i, pos, nullptr);
				Real coef3 = Derivation::DerivePartial<N>(g, i, j, l, pos, nullptr);

				gamma_ijk += 0.5 * g.Component(i, l, pos) * (coef1 + coef2 - coef3);
			}
			return gamma_ijk;
		}

		VectorN<Real, N> CovariantDerivativeContravar(const IVectorFunction<N>& func, int j, const VectorN<Real, N>& pos) const
		{
			MetricTensorField<N>& g = *this;

			VectorN<Real, N> ret;
			VectorN<Real, N> vec_val = func(pos);

			for (int i = 0; i < N; i++) {
				Real comp_val = Derivation::DeriveVecPartial<N>(func, i, j, pos, nullptr);

				for (int k = 0; k < N; k++)
					comp_val += GetChristoffelSymbolSecondKind(i, k, j, pos) * vec_val[k];

				ret[i] = comp_val;
			}
			return ret;
		}

		Real CovariantDerivativeContravarComp(const IVectorFunction<N>& func, int i, int j, const VectorN<Real, N>& pos) const
		{
			MetricTensorField<N>& g = *this;

			Real ret = Derivation::DeriveVecPartial<N>(func, i, j, pos, nullptr);

			for (int k = 0; k < N; k++)
				ret += GetChristoffelSymbolSecondKind(i, k, j, pos) * func(pos)[k];

			return ret;
		}

		VectorN<Real, N> CovariantDerivativeCovar(const IVectorFunction<N>& func, int j, const VectorN<Real, N>& pos) const
		{
			MetricTensorField<N>& g = *this;

			VectorN<Real, N> ret;
			VectorN<Real, N> vec_val = func(pos);

			for (int i = 0; i < N; i++) {
				Real comp_val = Derivation::DeriveVecPartial<N>(func, i, j, pos, nullptr);

				for (int k = 0; k < N; k++)
					comp_val -= GetChristoffelSymbolSecondKind(k, i, j, pos) * vec_val[k];

				ret[i] = comp_val;
			}
			return ret;
		}

		Real CovariantDerivativeCovarComp(const IVectorFunction<N>& func, int i, int j, const VectorN<Real, N>& pos) const
		{
			MetricTensorField<N>& g = *this;

			Real comp_val = Derivation::DeriveVecPartial<N>(func, i, j, pos, nullptr);

			for (int k = 0; k < N; k++)
				comp_val -= GetChristoffelSymbolSecondKind(k, i, j, pos) * func(pos)[k];

			return comp_val;
		}
	};

	template<int N>
	class MetricTensorCartesian : public MetricTensorField<N>
	{
	public:
		MetricTensorCartesian() : MetricTensorField<N>(2, 0) { }

		Real Component(int i, int j, const VectorN<Real, N>& pos) const
		{
			if (i == j)
				return 1.0;
			else
				return 0.0;
		}
	};

	class MetricTensorSpherical : public MetricTensorField<3>
	{
	public:
		MetricTensorSpherical() : MetricTensorField<3>(0, 2) { }

		virtual  Real Component(int i, int j, const VectorN<Real, 3>& pos) const override
		{
			if (i == 0 && j == 0)
				return 1.0;
			else if (i == 1 && j == 1)
				return POW2(pos[0]);
			else if (i == 2 && j == 2)
			{
				auto d = pos[0] * pos[0] * sin(pos[1]) * sin(pos[1]);
				return d;
			}
			else
				return 0.0;
		}
	};
	class MetricTensorSphericalContravar : public MetricTensorField<3>
	{
	public:
		MetricTensorSphericalContravar() : MetricTensorField<3>(2, 0) { }

		virtual Real Component(int i, int j, const VectorN<Real, 3>& pos) const
		{
			if (i == 0 && j == 0)
				return 1.0;
			else if (i == 1 && j == 1)
				return 1 / (pos[0] * pos[0]);
			else if (i == 2 && j == 2)
				return 1 / (pos[0] * pos[0] * sin(pos[1]) * sin(pos[1]));
			else
				return 0.0;
		}
	};
	class MetricTensorCylindrical : public MetricTensorField<3>
	{
	public:
		MetricTensorCylindrical() : MetricTensorField<3>(2, 0) { }

		virtual Real Component(int i, int j, const VectorN<Real, 3>& pos) const
		{
			if (i == 0 && j == 0)
				return 1.0;
			else if (i == 1 && j == 1)
				return pos[0] * pos[0];
			else if (i == 2 && j == 2)
				return 1.0;
			else
				return 0.0;
		}
	};

	template<typename VectorFrom, typename VectorTo, int N>
	class MetricTensorFromCoordTransf : public MetricTensorField<N>
	{
		const ICoordTransfWithInverse<VectorFrom, VectorTo, N>& _coordTransf;

	public:
		MetricTensorFromCoordTransf(ICoordTransfWithInverse<VectorFrom, VectorTo, N>& inTransf) : _coordTransf(inTransf)
		{ }

		virtual Real Component(int i, int j, const VectorN<Real, N>& pos) const
		{
			Real g_ij = 0.0;
			for (int k = 0; k < N; k++)
			{
				auto der_k_by_i = Derivation::DerivePartial<N>(_coordTransf.coordTransfFunc(k), i, pos, nullptr);
				auto der_k_by_j = Derivation::DerivePartial<N>(_coordTransf.coordTransfFunc(k), j, pos, nullptr);

				g_ij += der_k_by_i * der_k_by_j;
			}
			return g_ij;
		}
	};
}
///////////////////////////   ./include/core/CoordTransf.h   ///////////////////////////




namespace MML
{
	template<typename VectorFrom, typename VectorTo, int N>
	class CoordTransf : public virtual ICoordTransf<VectorFrom, VectorTo, N>
	{
	public:
		virtual VectorTo   getBasisVec(int ind, const VectorFrom& pos)
		{
			VectorTo ret;

			for (int i = 0; i < N; i++)
				ret[i] = Derivation::NDer4Partial(this->coordTransfFunc(i), ind, pos);

			return ret;
		}
		virtual VectorFrom getInverseBasisVec(int ind, const VectorFrom& pos)
		{
			VectorFrom ret;

			for (int i = 0; i < N; i++)
				ret[i] = Derivation::NDer4Partial(this->coordTransfFunc(ind), i, pos);

			return ret;
		}

		MatrixNM<Real, N, N> jacobian(const VectorN<Real, N>& pos)
		{
			MatrixNM<Real, N, N> jac;

			for (int i = 0; i < N; ++i)
				for (int j = 0; j < N; ++j)
				{
					jac(i, j) = Derivation::NDer4Partial(this->coordTransfFunc(i), j, pos);
				}

			return jac;
		}

		VectorTo   transfVecContravariant(const VectorFrom& vec, const VectorFrom& pos)
		{
			VectorFrom ret;
			for (int i = 0; i < N; i++) {
				ret[i] = 0;
				for (int j = 0; j < N; j++)
					ret[i] += Derivation::NDer4Partial(this->coordTransfFunc(i), j, pos) * vec[j];
			}
			return ret;
		}

		VectorFrom transfInverseVecCovariant(const VectorTo& vec, const VectorFrom& pos)
		{
			VectorFrom ret;
			for (int i = 0; i < N; i++)
			{
				ret[i] = 0;
				for (int j = 0; j < N; j++)
					ret[i] += Derivation::NDer4Partial(this->coordTransfFunc(j), i, pos) * vec[j];
			}
			return ret;
		}
	};

	template<typename VectorFrom, typename VectorTo, int N>
	class CoordTransfWithInverse : public virtual CoordTransf<VectorFrom, VectorTo, N>,
		public virtual ICoordTransfWithInverse<VectorFrom, VectorTo, N>
	{
	public:
		virtual VectorFrom getContravarBasisVec(int ind, const VectorTo& pos)
		{
			VectorFrom ret;

			for (int i = 0; i < N; i++)
				ret[i] = Derivation::NDer4Partial(this->inverseCoordTransfFunc(ind), i, pos);

			return ret;
		}
		virtual VectorTo   getInverseContravarBasisVec(int ind, const VectorTo& pos)
		{
			VectorFrom ret;

			for (int i = 0; i < N; i++)
				ret[i] = Derivation::NDer4Partial(this->inverseCoordTransfFunc(i), ind, pos);

			return ret;
		}

		VectorTo   transfVecCovariant(const VectorFrom& vec, const VectorTo& pos)
		{
			VectorTo ret;
			for (int i = 0; i < N; i++)
			{
				ret[i] = 0;
				for (int j = 0; j < N; j++)
					ret[i] += Derivation::NDer4Partial(this->inverseCoordTransfFunc(j), i, pos) * vec[j];
			}

			return ret;
		}

		VectorFrom transfInverseVecContravariant(const VectorTo& vec, const VectorTo& pos)
		{
			VectorFrom ret;
			for (int i = 0; i < N; i++)
			{
				ret[i] = 0;
				for (int j = 0; j < N; j++)
					ret[i] += Derivation::NDer4Partial(this->inverseCoordTransfFunc(i), j, pos, 1e-8) * vec[j];
			}

			return ret;
		}

		Tensor2<N> transfTensor2(const Tensor2<N>& tensor, const VectorFrom& pos)
		{
			Tensor2<N> ret(tensor.NumContravar(), tensor.NumCovar());

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					ret(i, j) = 0;
					for (int k = 0; k < N; k++)
						for (int l = 0; l < N; l++)
						{
							double coef1, coef2;
							if (tensor._isContravar[0])
								coef1 = Derivation::NDer1Partial(this->coordTransfFunc(i), k, pos);
							else
								coef1 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(k), i, pos);

							if (tensor._isContravar[1])
								coef2 = Derivation::NDer1Partial(this->coordTransfFunc(j), l, pos);
							else
								coef2 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(l), j, pos);

							ret(i, j) += coef1 * coef2 * tensor(k, l);
						}
				}

			return ret;
		}

		Tensor3<N> transfTensor3(const Tensor3<N>& tensor, const VectorFrom& pos)
		{
			Tensor3<N> ret(tensor.NumContravar(), tensor.NumCovar());

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
					{
						ret.Component(i, j, k) = 0;
						for (int l = 0; l < N; l++)
							for (int m = 0; m < N; m++)
								for (int n = 0; n < N; n++)
								{
									double coef1, coef2, coef3;
									if (tensor._isContravar[0])
										coef1 = Derivation::NDer1Partial(this->coordTransfFunc(i), l, pos);
									else
										coef1 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(l), i, pos);

									if (tensor._isContravar[1])
										coef2 = Derivation::NDer1Partial(this->coordTransfFunc(j), m, pos);
									else
										coef2 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(m), j, pos);

									if (tensor._isContravar[2])
										coef3 = Derivation::NDer1Partial(this->coordTransfFunc(k), n, pos);
									else
										coef3 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(n), k, pos);

									ret.Component(i, j, k) += coef1 * coef2 * coef3 * tensor.Component(l, m, n);
								}
					}

			return ret;
		}

		Tensor4<N> transfTensor4(const Tensor4<N>& tensor, const VectorFrom& pos)
		{
			Tensor4<N> ret(tensor.NumContravar(), tensor.NumCovar());

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
						for (int l = 0; l < N; l++)
						{
							ret[i][j][k][l] = 0;
							for (int m = 0; m < N; m++)
								for (int n = 0; n < N; n++)
									for (int o = 0; o < N; o++)
										for (int p = 0; p < N; p++)
										{
											double coef1, coef2, coef3, coef4;
											if (tensor._isContravar[0])
												coef1 = Derivation::NDer1Partial(this->coordTransfFunc(i), m, pos);
											else
												coef1 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(m), i, pos);

											if (tensor._isContravar[1])
												coef2 = Derivation::NDer1Partial(this->coordTransfFunc(j), n, pos);
											else
												coef2 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(n), j, pos);

											if (tensor._isContravar[2])
												coef3 = Derivation::NDer1Partial(this->coordTransfFunc(k), o, pos);
											else
												coef3 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(o), k, pos);

											if (tensor._isContravar[3])
												coef4 = Derivation::NDer1Partial(this->coordTransfFunc(l), p, pos);
											else
												coef4 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(p), l, pos);

											ret[i][j][k][l] += coef1 * coef2 * coef3 * coef4 * tensor[m][n][o][p];
										}
						}

			return ret;
		}

		Tensor5<N> transfTensor5(const Tensor5<N>& tensor, const VectorFrom& pos)
		{
			Tensor5<N> ret(tensor.NumContravar(), tensor.NumCovar());

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
						for (int l = 0; l < N; l++)
							for (int m = 0; m < N; m++)
							{
								ret[i][j][k][l][m] = 0;
								for (int n = 0; n < N; n++)
									for (int o = 0; o < N; o++)
										for (int p = 0; p < N; p++)
											for (int q = 0; q < N; q++)
												for (int r = 0; r < N; r++)
												{
													double coef1, coef2, coef3, coef4, coef5;
													if (tensor._isContravar[0])
														coef1 = Derivation::NDer1Partial(this->coordTransfFunc(i), n, pos);
													else
														coef1 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(n), i, pos);

													if (tensor._isContravar[1])
														coef2 = Derivation::NDer1Partial(this->coordTransfFunc(j), o, pos);
													else
														coef2 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(o), j, pos);

													if (tensor._isContravar[2])
														coef3 = Derivation::NDer1Partial(this->coordTransfFunc(k), p, pos);
													else
														coef3 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(p), k, pos);

													if (tensor._isContravar[3])
														coef4 = Derivation::NDer1Partial(this->coordTransfFunc(l), q, pos);
													else
														coef4 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(q), l, pos);

													if (tensor._isContravar[4])
														coef4 = Derivation::NDer1Partial(this->coordTransfFunc(m), r, pos);
													else
														coef4 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(r), m, pos);

													ret[i][j][k][l][m] += coef1 * coef2 * coef3 * coef4 * coef5 * tensor[n][o][p][q][r];
												}
							}

			return ret;
		}
	};
}
///////////////////////////   ./include/core/CoordTransf/CoordTransf2D.h   ///////////////////////////



namespace MML
{
	class CoordTransfPolarToCartesian2D : public CoordTransfWithInverse<Vector2Polar, Vector2Cartesian, 2>
	{
		// q[0] = r     - radial distance
		// q[1] = phi   - polar angle
	public:
		static Real func1(const VectorN<Real, 2>& q) { return q[0] * cos(q[1]); }
		static Real func2(const VectorN<Real, 2>& q) { return q[0] * sin(q[1]); }

		// q[0] = x
		// q[1] = y
		static Real funcInverse1(const VectorN<Real, 2>& q) { return sqrt(q[0] * q[0] + q[1] * q[1]); }
		static Real funcInverse2(const VectorN<Real, 2>& q) { return atan2(q[1], q[0]); }

		inline static ScalarFunction<2> _func[2] = { ScalarFunction<2>{func1},
																								 ScalarFunction<2>{func2}
		};

		inline static ScalarFunction<2> _funcInverse[2] = { ScalarFunction<2>{funcInverse1},
																												ScalarFunction<2>{funcInverse2}
		};

		Vector2Cartesian     transf(const Vector2Polar& q) const { return Vector2Cartesian{ func1(q), func2(q) }; }
		Vector2Polar         transfInverse(const Vector2Cartesian& q) const { return Vector2Polar{ funcInverse1(q), funcInverse2(q) }; }

		const IScalarFunction<2>& coordTransfFunc(int i) const { return _func[i]; }
		const IScalarFunction<2>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
	};

	class CoordTransfCart2DRotation : public CoordTransfWithInverse<Vector2Cartesian, Vector2Cartesian, 2>
	{
	private:
		Real    _angle;
		MatrixNM<Real, 2, 2>  _transf;
		MatrixNM<Real, 2, 2>  _inverse;

		const ScalarFunctionFromStdFunc<2> _f1;
		const ScalarFunctionFromStdFunc<2> _f2;

		const ScalarFunctionFromStdFunc<2> _fInverse1;
		const ScalarFunctionFromStdFunc<2> _fInverse2;

	public:

		CoordTransfCart2DRotation(Real inAngle) : _angle(inAngle),
			_f1(std::function<Real(const VectorN<Real, 2>&)> { std::bind(&CoordTransfCart2DRotation::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 2>&)> { std::bind(&CoordTransfCart2DRotation::func2, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 2>&)> { std::bind(&CoordTransfCart2DRotation::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 2>&)> { std::bind(&CoordTransfCart2DRotation::funcInverse2, this, std::placeholders::_1) })
		{
			_transf[0][0] = cos(_angle);
			_transf[0][1] = -sin(_angle);
			_transf[1][0] = sin(_angle);
			_transf[1][1] = cos(_angle);

			_inverse[0][0] = cos(_angle);
			_inverse[0][1] = sin(_angle);
			_inverse[1][0] = -sin(_angle);
			_inverse[1][1] = cos(_angle);
		}

		Real func1(const VectorN<Real, 2>& q) const { return _transf[0][0] * q[0] + _transf[0][1] * q[1]; }
		Real func2(const VectorN<Real, 2>& q) const { return (_transf * q)[1]; }

		Real funcInverse1(const VectorN<Real, 2>& q) const { return (_inverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 2>& q) const { return (_inverse * q)[1]; }

		Vector2Cartesian    transf(const Vector2Cartesian& q) const { return Vector2Cartesian{ func1(q), func2(q) }; }
		Vector2Cartesian    transfInverse(const Vector2Cartesian& q) const { return Vector2Cartesian{ funcInverse1(q), funcInverse2(q) }; }

		const IScalarFunction<2>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else return _f2;
		}
		const IScalarFunction<2>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else return _fInverse2;
		}
	};
}

///////////////////////////   ./include/core/CoordTransf/CoordTransf3D.h   ///////////////////////////



namespace MML
{
	class CoordTransfCart3DRotationXAxis : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	private:
		Real    _angle;
		MatrixNM<Real, 3, 3>  _transf;
		MatrixNM<Real, 3, 3>  _inverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

	public:
		CoordTransfCart3DRotationXAxis(Real inAngle) : _angle(inAngle),
			_f1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationXAxis::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationXAxis::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationXAxis::func3, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationXAxis::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationXAxis::funcInverse2, this, std::placeholders::_1) }),
			_fInverse3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationXAxis::funcInverse3, this, std::placeholders::_1) })
		{
			_transf[0][0] = 1.0;
			_transf[1][1] = cos(_angle);
			_transf[1][2] = -sin(_angle);
			_transf[2][1] = sin(_angle);
			_transf[2][2] = cos(_angle);

			_inverse[0][0] = 1.0;
			_inverse[1][1] = cos(_angle);
			_inverse[1][2] = sin(_angle);
			_inverse[2][1] = -sin(_angle);
			_inverse[2][2] = cos(_angle);
		}

		Real func1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_inverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_inverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_inverse * q)[2]; }

		Vector3Cartesian    transf(const Vector3Cartesian& q) const { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
		Vector3Cartesian    transfInverse(const Vector3Cartesian& q) const { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i)const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}
	};
	class CoordTransfCart3DRotationYAxis : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	private:
		Real    _angle;
		MatrixNM<Real, 3, 3>  _transf;
		MatrixNM<Real, 3, 3>  _inverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

	public:
		CoordTransfCart3DRotationYAxis(Real inAngle) : _angle(inAngle),
			_f1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationYAxis::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationYAxis::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationYAxis::func3, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationYAxis::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationYAxis::funcInverse2, this, std::placeholders::_1) }),
			_fInverse3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationYAxis::funcInverse3, this, std::placeholders::_1) })
		{
			_transf[0][0] = cos(_angle);
			_transf[0][2] = sin(_angle);
			_transf[1][1] = 1.0;
			_transf[2][0] = -sin(_angle);
			_transf[2][2] = cos(_angle);

			_inverse[0][0] = cos(_angle);
			_inverse[0][2] = -sin(_angle);
			_inverse[1][1] = 1.0;
			_inverse[2][0] = sin(_angle);
			_inverse[2][2] = cos(_angle);
		}

		Real func1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_inverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_inverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_inverse * q)[2]; }

		Vector3Cartesian    transf(const Vector3Cartesian& q) const { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
		Vector3Cartesian    transfInverse(const Vector3Cartesian& q) const { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}
	};
	class CoordTransfCart3DRotationZAxis : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	private:
		Real    _angle;
		MatrixNM<Real, 3, 3>  _transf;
		MatrixNM<Real, 3, 3>  _inverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

	public:
		CoordTransfCart3DRotationZAxis(Real inAngle) : _angle(inAngle),
			_f1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationZAxis::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationZAxis::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationZAxis::func3, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationZAxis::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationZAxis::funcInverse2, this, std::placeholders::_1) }),
			_fInverse3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationZAxis::funcInverse3, this, std::placeholders::_1) })
		{
			_transf[0][0] = cos(_angle);
			_transf[0][1] = -sin(_angle);
			_transf[1][0] = sin(_angle);
			_transf[1][1] = cos(_angle);
			_transf[2][2] = 1.0;

			_inverse[0][0] = cos(_angle);
			_inverse[0][1] = sin(_angle);
			_inverse[1][0] = -sin(_angle);
			_inverse[1][1] = cos(_angle);
			_inverse[2][2] = 1.0;
		}

		Real func1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_inverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_inverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_inverse * q)[2]; }

		Vector3Cartesian    transf(const Vector3Cartesian& q) const { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
		Vector3Cartesian    transfInverse(const Vector3Cartesian& q) const { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}
	};

	// TODO 0.9 - LOW, TESKO, Euler angles, finish CoordTransfCart3DRotationGeneralAxis
	// class CoordTransfCart3DRotationGeneralAxis  : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	// {

	// };

	// Performs tranformation from original (Cartesian) system to orthogonal system defined by 
	// its base vectors expressed in original system.
	class CoordTransf3DCartOrthogonal : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	private:
		Vector3Cartesian _base[3];

		MatrixNM<Real, 3, 3> _transf;
		MatrixNM<Real, 3, 3> _transfInverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

		Real func1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[2]; }

	public:
		CoordTransf3DCartOrthogonal(const VectorN<Real, 3>& b1, const VectorN<Real, 3>& b2, const VectorN<Real, 3>& b3) :
			_f1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOrthogonal::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOrthogonal::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOrthogonal::func3, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOrthogonal::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOrthogonal::funcInverse2, this, std::placeholders::_1) }),
			_fInverse3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOrthogonal::funcInverse3, this, std::placeholders::_1) })
		{
			_base[0] = b1;
			_base[1] = b2;
			_base[2] = b3;

			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++)
				{
					_transf(i, j) = _base[i][j];
					_transfInverse(i, j) = _base[j][i];
				}
			}
		}

		MatrixNM<Real, 3, 3> getTransfMatrix() { return _transf; }
		MatrixNM<Real, 3, 3> getInvTransfMatrix() { return _transfInverse; }

		Vector3Cartesian    transf(const Vector3Cartesian& q) const { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
		Vector3Cartesian    transfInverse(const Vector3Cartesian& q) const { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}

		bool IsRightHanded()
		{
			Vector3Cartesian cross = VectorProd(_base[0], _base[1]);
			if (ScalarProd(cross, _base[2]) > 0.0)
				return true;
			else
				return false;
		}
	};

	// General 3D Cartesian transformation, given by matrix
	class CoordTransf3DCartGeneral : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	private:
		Vector3Cartesian _base[3];

		MatrixNM<Real, 3, 3> _transf;
		MatrixNM<Real, 3, 3> _transfInverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

		Real func1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[2]; }

	public:
		CoordTransf3DCartGeneral(const MatrixNM<Real, 3, 3>& transfMat) :
			_f1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartGeneral::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartGeneral::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartGeneral::func3, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartGeneral::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartGeneral::funcInverse2, this, std::placeholders::_1) }),
			_fInverse3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartGeneral::funcInverse3, this, std::placeholders::_1) })
		{
			_transf = transfMat;

			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++)
				{
					_base[i][j] = _transf(i, j);
					_transfInverse(i, j) = _transf[j][i];
				}
			}

			// check if it is orthogonal!
			// TODO - implement this

		}

		MatrixNM<Real, 3, 3> getTransfMatrix() { return _transf; }
		MatrixNM<Real, 3, 3> getInvTransfMatrix() { return _transfInverse; }

		Vector3Cartesian    transf(const Vector3Cartesian& q) const { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
		Vector3Cartesian    transfInverse(const Vector3Cartesian& q) const { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}

		bool IsRightHanded()
		{
			Vector3Cartesian cross = VectorProd(_base[0], _base[1]);
			if (ScalarProd(cross, _base[2]) > 0.0)
				return true;
			else
				return false;
		}
	};
	class CoordTransf3DCartOblique : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	private:
		Vector3Cartesian _base[3];
		Vector3Cartesian _dual[3];

		MatrixNM<Real, 3, 3> _alpha;
		MatrixNM<Real, 3, 3> _transf;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

		// TODO - ovo popraviti i napraviti kako spada
		Real func1(const VectorN<Real, 3>& q) const { return ScalarProd(q, MML::Vector3Cartesian(_dual[0])); }
		Real func2(const VectorN<Real, 3>& q) const { return ScalarProd(q, MML::Vector3Cartesian(_dual[1])); }
		Real func3(const VectorN<Real, 3>& q) const { return ScalarProd(q, MML::Vector3Cartesian(_dual[2])); }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

	public:
		CoordTransf3DCartOblique(VectorN<Real, 3> b1, VectorN<Real, 3> b2, VectorN<Real, 3> b3) :
			_f1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOblique::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOblique::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOblique::func3, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOblique::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOblique::funcInverse2, this, std::placeholders::_1) }),
			_fInverse3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOblique::funcInverse3, this, std::placeholders::_1) })
		{
			_base[0] = b1;
			_base[1] = b2;
			_base[2] = b3;

			Real V = ScalarProd(_base[0], VectorProd(_base[1], _base[2]));

			_dual[0] = (1 / V) * VectorProd(_base[1], _base[2]);
			_dual[1] = (1 / V) * VectorProd(_base[2], _base[0]);
			_dual[2] = (1 / V) * VectorProd(_base[0], _base[1]);

			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					_alpha(i, j) = _base[i][j];
					_transf(i, j) = _base[j][i];     // transponirano
				}
			}
		}

		MatrixNM<Real, 3, 3> getAlpha() { return _alpha; }
		MatrixNM<Real, 3, 3> getTransf() { return _alpha; }

		Vector3Cartesian    Base(int i) { return _base[i]; }
		Vector3Cartesian    Dual(int i) { return _dual[i]; }

		Vector3Cartesian    transf(const Vector3Cartesian& q) const { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
		Vector3Cartesian    transfInverse(const Vector3Cartesian& q) const { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}

		bool IsRightHanded()
		{
			Vector3Cartesian cross = VectorProd(_base[0], _base[1]);
			if (ScalarProd(cross, _base[2]) > 0.0)
				return true;
			else
				return false;
		}
	};
}

///////////////////////////   ./include/core/CoordTransf/CoordTransfSpherical.h   ///////////////////////////



namespace MML
{

	// TODO 0.9 - VIDJETI STO S MATH vs PHY konvencijama o redoslijedu koordinata
	class CoordTransfSphericalToCartesian : public CoordTransfWithInverse<Vector3Spherical, Vector3Cartesian, 3>
	{
	private:
		// q[0] = r     - radial distance
		// q[1] = theta - inclination
		// q[2] = phi   - azimuthal angle
		static Real x(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]) * cos(q[2]); }
		static Real y(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]) * sin(q[2]); }
		static Real z(const VectorN<Real, 3>& q) { return q[0] * cos(q[1]); }

		// q[0] = x
		// q[1] = y
		// q[2] = z
		static Real r(const VectorN<Real, 3>& q) { return sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]); }
		static Real theta(const VectorN<Real, 3>& q) { return acos(q[2] / sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2])); }
		static Real phi(const VectorN<Real, 3>& q) { return atan2(q[1], q[0]); }

		inline static ScalarFunction<3> _func[3] = { ScalarFunction<3>{x},
																								 ScalarFunction<3>{y},
																								 ScalarFunction<3>{z}
		};

		inline static ScalarFunction<3> _funcInverse[3] = { ScalarFunction<3>{r},
																												ScalarFunction<3>{theta},
																												ScalarFunction<3>{phi}
		};
	public:
		Vector3Cartesian     transf(const Vector3Spherical& q) const { return Vector3Cartesian{ x(q), y(q), z(q) }; }
		Vector3Spherical     transfInverse(const Vector3Cartesian& q) const { return Vector3Spherical{ r(q), theta(q), phi(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const { return _func[i]; }
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }

		Vector3Cartesian getBasisVectorExplicit(int ind, const Vector3Spherical& pos)
		{
			const Real r = pos[0];
			const Real theta = pos[1];
			const Real phi = pos[2];
			switch (ind)
			{
			case 0: return Vector3Cartesian{      sin(theta) * cos(phi),     sin(theta) * sin(phi),      cos(theta) };
			case 1: return Vector3Cartesian{  r * cos(theta) * cos(phi), r * cos(theta) * sin(phi), -r * sin(theta) };
			case 2: return Vector3Cartesian{ -r * sin(theta) * sin(phi), r * sin(theta) * cos(phi),						  0.0 };
			default:
				return Vector3Cartesian{ 0.0, 0.0, 0.0 };
			}
		}
		Vector3Cartesian getUnitBasisVectorExplicit(int ind, const Vector3Spherical& pos)
		{
			const Real r = pos[0];
			const Real theta = pos[1];
			const Real phi = pos[2];
			switch (ind)
			{
			case 0: return Vector3Cartesian{ sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta) };
			case 1: return Vector3Cartesian{ cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta) };
			case 2: return Vector3Cartesian{ -sin(phi), cos(phi), 0.0 };
			default:
				return Vector3Cartesian{ 0.0, 0.0, 0.0 };
			}
		}

		Vector3Spherical getInverseBasisVectorExplicit(int ind, const Vector3Spherical& pos)
		{
			const Real r = pos[0];
			const Real theta = pos[1];
			const Real phi = pos[2];
			switch(ind)
			{
			case 0: return Vector3Spherical{ sin(theta) * cos(phi), r * cos(theta) * cos(phi), -r * sin(theta) * sin(phi) };
			case 1: return Vector3Spherical{ sin(theta) * sin(phi), r * cos(theta) * sin(phi),  r * sin(theta) * cos(phi) };
			case 2: return Vector3Spherical{ cos(theta)           ,-r * sin(theta)           ,                        0.0 };
			default: 
				return Vector3Spherical{ 0.0, 0.0, 0.0 };
			}
		}
		Vector3Spherical getInverseUnitBasisVectorExplicit(int ind, const Vector3Spherical& pos)
		{
			const Real r = pos[0];
			const Real theta = pos[1];
			const Real phi = pos[2];
			switch(ind)
			{
			case 0: return Vector3Spherical{ sin(theta) * cos(phi),  cos(theta) * cos(phi), -sin(phi) };
			case 1: return Vector3Spherical{ sin(theta) * sin(phi),  cos(theta) * sin(phi),  cos(phi) };
			case 2: return Vector3Spherical{ cos(theta)           , -sin(theta)           ,  0.0 };
			default: 
				return Vector3Spherical{ 0.0, 0.0, 0.0 };
			}
		}
	};

	class CoordTransfCartesianToSpherical : public CoordTransfWithInverse<Vector3Cartesian, Vector3Spherical, 3>
	{
	private:
		// q[0] = x
		// q[1] = y
		// q[2] = z
		static Real r(const VectorN<Real, 3>& q) { return sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]); }
		static Real theta(const VectorN<Real, 3>& q) { return acos(q[2] / sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2])); }
		static Real phi(const VectorN<Real, 3>& q) { return atan2(q[1], q[0]); }

		// q[0] = r     - radial distance
		// q[1] = theta - inclination
		// q[2] = phi   - azimuthal angle
		static Real x(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]) * cos(q[2]); }
		static Real y(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]) * sin(q[2]); }
		static Real z(const VectorN<Real, 3>& q) { return q[0] * cos(q[1]); }

		inline static ScalarFunction<3> _func[3] = { ScalarFunction<3>{r},
																								 ScalarFunction<3>{theta},
																								 ScalarFunction<3>{phi}
		};

		inline static ScalarFunction<3> _funcInverse[3] = { ScalarFunction<3>{x},
																												ScalarFunction<3>{y},
																												ScalarFunction<3>{z}
		};
	public:
		Vector3Spherical     transf(const Vector3Cartesian& q) const { return Vector3Spherical{ r(q), theta(q), phi(q) }; }
		Vector3Cartesian     transfInverse(const Vector3Spherical& q) const { return Vector3Cartesian{ x(q), y(q), z(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const { return _func[i]; }
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
	};


	static CoordTransfSphericalToCartesian      CoordTransfSpherToCart;
	static CoordTransfCartesianToSpherical      CoordTransfCartToSpher;
}

///////////////////////////   ./include/core/CoordTransf/CoordTransfCylindrical.h   ///////////////////////////



namespace MML
{

	class CoordTransfCylindricalToCartesian : public CoordTransfWithInverse<Vector3Cylindrical, Vector3Cartesian, 3>
	{
	private:
		// q1 = r   - distance from symmetry axis
		// q2 = phi - angle to symmetry axis
		// q3 = z   - z
		static Real x(const VectorN<Real, 3>& q) { return q[0] * cos(q[1]); }
		static Real y(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]); }
		static Real z(const VectorN<Real, 3>& q) { return q[2]; }

		// q[0] = x
		// q[1] = y
		// q[2] = z
		static Real r(const VectorN<Real, 3>& q) { return sqrt(q[0] * q[0] + q[1] * q[1]); }
		static Real phi(const VectorN<Real, 3>& q) { return atan2(q[1], q[0]); }
		// z-coordinate is the same, so we'll use for inverse the same function as for forward transformation

		inline static ScalarFunction<3> _func[3] = { ScalarFunction<3>{x},
																								 ScalarFunction<3>{y},
																								 ScalarFunction<3>{z}
		};

		inline static ScalarFunction<3> _funcInverse[3] = { ScalarFunction<3>{r},
																												ScalarFunction<3>{phi},
																												ScalarFunction<3>{z}
		};
	public:
		Vector3Cartesian     transf(const Vector3Cylindrical& q)      const { return Vector3Cartesian{ x(q), y(q), z(q) }; }
		Vector3Cylindrical   transfInverse(const Vector3Cartesian& q) const { return Vector3Cylindrical{ r(q), phi(q), z(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i)        const { return _func[i]; }
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
	};

	class CoordTransfCartesianToCylindrical : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cylindrical, 3>
	{
	private:
		// q[0] = x
		// q[1] = y
		// q[2] = z
		static Real r(const VectorN<Real, 3>& q) { return sqrt(q[0] * q[0] + q[1] * q[1]); }
		static Real phi(const VectorN<Real, 3>& q) { return atan2(q[1], q[0]); }
		static Real z(const VectorN<Real, 3>& q) { return q[2]; }

		// q1 = r   - distance from symmetry axis
		// q2 = phi - angle to symmetry axis
		// q3 = z   - z
		static Real x(const VectorN<Real, 3>& q) { return q[0] * cos(q[1]); }
		static Real y(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]); }

		inline static ScalarFunction<3> _func[3] = { ScalarFunction<3>{r},
																								 ScalarFunction<3>{phi},
																								 ScalarFunction<3>{z}
		};

		inline static ScalarFunction<3> _funcInverse[3] = { ScalarFunction<3>{x},
																												ScalarFunction<3>{y},
																												ScalarFunction<3>{z}
		};
	public:
		Vector3Cylindrical transf(const Vector3Cartesian& q)          const { return Vector3Cylindrical{ r(q), phi(q), z(q) }; }
		Vector3Cartesian   transfInverse(const Vector3Cylindrical& q) const { return Vector3Cartesian{ x(q), y(q), z(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i)        const { return _func[i]; }
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
	};

	static CoordTransfCartesianToCylindrical    CoordTransfCartToCyl;
	static CoordTransfCylindricalToCartesian    CoordTransfCylToCart;

}

///////////////////////////   ./include/core/CoordTransf/CoordTransfLorentz.h   ///////////////////////////




namespace MML
{
	class CoordTransfLorentzXAxis : public CoordTransfWithInverse<Vector4Lorentz, Vector4Lorentz, 4>
	{
	private:
		Real    _velocity;			// expressed in units of c (speed of light)
		MatrixNM<Real, 4, 4>  _transf;
		MatrixNM<Real, 4, 4>  _inverse;

		const ScalarFunctionFromStdFunc<4> _f1, _f2, _f3, _f4;
		const ScalarFunctionFromStdFunc<4> _fInverse1, _fInverse2, _fInverse3, _fInverse4;

	public:
		CoordTransfLorentzXAxis(Real inVelocity) : _velocity(inVelocity),
			_f1(std::function<Real(const VectorN<Real, 4>&)> { std::bind(&CoordTransfLorentzXAxis::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 4>&)> { std::bind(&CoordTransfLorentzXAxis::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 4>&)> { std::bind(&CoordTransfLorentzXAxis::func3, this, std::placeholders::_1) }),
			_f4(std::function<Real(const VectorN<Real, 4>&)> { std::bind(&CoordTransfLorentzXAxis::func4, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 4>&)> { std::bind(&CoordTransfLorentzXAxis::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 4>&)> { std::bind(&CoordTransfLorentzXAxis::funcInverse2, this, std::placeholders::_1) }),
			_fInverse3(std::function<Real(const VectorN<Real, 4>&)> { std::bind(&CoordTransfLorentzXAxis::funcInverse3, this, std::placeholders::_1) }),
			_fInverse4(std::function<Real(const VectorN<Real, 4>&)> { std::bind(&CoordTransfLorentzXAxis::funcInverse4, this, std::placeholders::_1) })
		{
			double gamma = 1.0 / sqrt(1.0 - _velocity * _velocity);
			double beta = _velocity;

			_transf[0][0] = gamma;
			_transf[0][1] = -beta * gamma;
			_transf[1][0] = -beta * gamma;
			_transf[1][1] = gamma;
			_transf[2][2] = 1.0;
			_transf[3][3] = 1.0;

			_inverse[0][0] = gamma;
			_inverse[0][1] = -beta * gamma;
			_inverse[1][0] = -beta * gamma;
			_inverse[1][1] = gamma;
			_inverse[2][2] = 1.0;
			_inverse[3][3] = 1.0;
		}

		Real func1(const VectorN<Real, 4>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 4>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 4>& q) const { return (_transf * q)[2]; }
		Real func4(const VectorN<Real, 4>& q) const { return (_transf * q)[3]; }

		Real funcInverse1(const VectorN<Real, 4>& q) const { return (_inverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 4>& q) const { return (_inverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 4>& q) const { return (_inverse * q)[2]; }
		Real funcInverse4(const VectorN<Real, 4>& q) const { return (_inverse * q)[3]; }

		Vector4Lorentz    transf(const Vector4Lorentz& q) const { return Vector4Lorentz{ func1(q), func2(q), func3(q), func4(q) }; }
		Vector4Lorentz    transfInverse(const Vector4Lorentz& q) const { return Vector4Lorentz{ funcInverse1(q), funcInverse2(q), funcInverse3(q), funcInverse4(q) }; }

		const IScalarFunction<4>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else if (i == 2) return _f3;
			else return _f4;
		}
		const IScalarFunction<4>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else if (i == 2) return _fInverse3;
			else return _fInverse4;
		}
	};
}

///////////////////////////   ./include/core/FieldOperations.h   ///////////////////////////// grad
// - cart
// - spher
// - cyl






namespace MML
{
	namespace ScalarFieldOperations
	{
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////                 GENERAL COORD. FIELD OPERATIONS                    ///////////////////
		template<int N>
		static VectorN<Real, N> Gradient(IScalarFunction<N>& scalarField, const VectorN<Real, N>& pos, const MetricTensorField<N>& metricTensorField)
		{
			VectorN<Real, N> derivsAtPoint = Derivation::DerivePartialAll<N>(scalarField, pos, nullptr);

			// TODO - depends on covar-contravar of metric tensor!!!
			Tensor2<N> metricAtPoint(2, 0);
			metricTensorField.ValueAtPoint(pos, metricAtPoint);

			VectorN<Real, N> ret = metricAtPoint * derivsAtPoint;

			return ret;
		}

		template<int N>
		static Real Divergence(const IVectorFunction<N>& vectorField, const VectorN<Real, N>& pos, const MetricTensorField<N>& metricTensorField)
		{
			Real div = 0.0;
			VectorN<Real, N> vec_val = vectorField(pos);

			for (int i = 0; i < N; i++)
			{
				div += Derivation::DeriveVecPartial<N>(vectorField, i, i, pos, nullptr);

				// correction for general coordinates
				for (int k = 0; k < N; k++)
				{
					div += vec_val[k] * metricTensorField.GetChristoffelSymbolSecondKind(i, i, k, pos);
				}
			}
			return div;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////                   GRADIENT                     /////////////////////////////
		template<int N>
		static VectorN<Real, N> GradientCart(const IScalarFunction<N>& scalarField, const VectorN<Real, N>& pos)
		{
			return Derivation::DerivePartialAll<N>(scalarField, pos, nullptr);
		}
		template<int N>
		static VectorN<Real, N> GradientCart(const IScalarFunction<N>& scalarField, const VectorN<Real, N>& pos, int der_order)
		{
			switch (der_order)
			{
			case 1: return Derivation::NDer1PartialByAll<N>(scalarField, pos, nullptr);
			case 2: return Derivation::NDer2PartialByAll<N>(scalarField, pos, nullptr);
			case 4: return Derivation::NDer4PartialByAll<N>(scalarField, pos, nullptr);
			case 6: return Derivation::NDer6PartialByAll<N>(scalarField, pos, nullptr);
			case 8: return Derivation::NDer8PartialByAll<N>(scalarField, pos, nullptr);
			default:
				throw std::invalid_argument("GradientCart: der_order must be in 1, 2, 4, 6 or 8");
			}
		}

		static Vector3Spherical GradientSpher(const IScalarFunction<3>& scalarField, const Vector3Spherical& pos)
		{
			Vector3Spherical ret = Derivation::DerivePartialAll<3>(scalarField, pos, nullptr);

			ret[1] = ret[1] / pos[0];
			ret[2] = ret[2] / (pos[0] * sin(pos[1]));

			return ret;
		}
		static Vector3Spherical GradientSpher(const IScalarFunction<3>& scalarField, const Vector3Spherical& pos, int der_order)
		{
			Vector3Spherical ret;

			switch (der_order)
			{
			case 1: ret = Derivation::NDer1PartialByAll<3>(scalarField, pos, nullptr); break;
			case 2: ret = Derivation::NDer2PartialByAll<3>(scalarField, pos, nullptr); break;
			case 4: ret = Derivation::NDer4PartialByAll<3>(scalarField, pos, nullptr); break;
			case 6: ret = Derivation::NDer6PartialByAll<3>(scalarField, pos, nullptr); break;
			case 8: ret = Derivation::NDer8PartialByAll<3>(scalarField, pos, nullptr); break;
			default:
				throw std::invalid_argument("GradientSpher: der_order must be in 1, 2, 4, 6 or 8");
			}

			ret[1] = ret[1] / pos[0];
			ret[2] = ret[2] / (pos[0] * sin(pos[1]));

			return ret;
		}

		static Vector3Cylindrical GradientCyl(const IScalarFunction<3>& scalarField, const Vector3Cylindrical& pos)
		{
			Vector3Cylindrical ret = Derivation::DerivePartialAll<3>(scalarField, pos, nullptr);

			ret[1] = ret[1] / pos[0];

			return ret;
		}
		static Vector3Cylindrical GradientCyl(const IScalarFunction<3>& scalarField, const Vector3Cylindrical& pos, int der_order)
		{
			Vector3Cylindrical ret;

			switch (der_order)
			{
			case 1: ret = Derivation::NDer1PartialByAll<3>(scalarField, pos, nullptr); break;
			case 2: ret = Derivation::NDer2PartialByAll<3>(scalarField, pos, nullptr); break;
			case 4: ret = Derivation::NDer4PartialByAll<3>(scalarField, pos, nullptr); break;
			case 6: ret = Derivation::NDer6PartialByAll<3>(scalarField, pos, nullptr); break;
			case 8: ret = Derivation::NDer8PartialByAll<3>(scalarField, pos, nullptr); break;
			default:
				throw std::invalid_argument("GradientCyl: der_order must be in 1, 2, 4, 6 or 8");
			}
			ret[1] = ret[1] / pos[0];

			return ret;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////                  LAPLACIAN                     /////////////////////////////
		template<int N>
		static Real LaplacianCart(const IScalarFunction<N>& scalarField, const VectorN<Real, N>& pos)
		{
			Real lapl = 0.0;
			for (int i = 0; i < N; i++)
				lapl += Derivation::DeriveSecPartial<N>(scalarField, i, i, pos, nullptr);

			return lapl;
		}
		static Real LaplacianSpher(const IScalarFunction<3>& scalarField, const Vector3Spherical& pos)
		{
			const Real r = pos.R();
			const Real phi = pos.Phi();
			const Real theta = pos.Theta();

			Real first = Derivation::DeriveSecPartial<3>(scalarField, 0, 0, pos, nullptr);
			Real second = 2 / pos.R() * Derivation::DerivePartial<3>(scalarField, 0, pos, nullptr);
			Real third = 1 / (r * r * sin(theta)) * (cos(theta) * Derivation::DerivePartial<3>(scalarField, 1, pos, nullptr) + sin(theta) * Derivation::DeriveSecPartial<3>(scalarField, 1, 1, pos, nullptr));
			Real fourth = 1 / (r * r * sin(theta) * sin(theta));

			return first + second + third;
		}
		static Real LaplacianCyl(const IScalarFunction<3>& scalarField, const Vector3Cylindrical& pos)
		{
			const Real r = pos[0];

			Real first = 1 / r * (Derivation::DerivePartial<3>(scalarField, 0, pos, nullptr) + r * Derivation::DeriveSecPartial<3>(scalarField, 0, 0, pos, nullptr));
			Real second = 1 / (r * r) * Derivation::DeriveSecPartial<3>(scalarField, 1, 1, pos, nullptr);
			Real third = Derivation::DeriveSecPartial<3>(scalarField, 2, 2, pos, nullptr);

			return first + second + third;
		}
	};

	namespace VectorFieldOperations
	{
		////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////                  DIVERGENCE                    /////////////////////////////
		template<int N>
		static Real DivCart(const IVectorFunction<N>& vectorField, const VectorN<Real, N>& pos)
		{
			Real div = 0.0;
			for (int i = 0; i < N; i++)
				div += Derivation::DeriveVecPartial<N>(vectorField, i, i, pos, nullptr);

			return div;
		}

		static Real DivSpher(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& x)
		{
			VectorN<Real, 3> vals = vectorField(x);

			VectorN<Real, 3> derivs;
			for (int i = 0; i < 3; i++)
				derivs[i] = Derivation::DeriveVecPartial<3>(vectorField, i, i, x, nullptr);

			Real div = 0.0;
			div += 1 / (x[0] * x[0]) * (2 * x[0] * vals[0] + x[0] * x[0] * derivs[0]);
			div += 1 / (x[0] * sin(x[1])) * (cos(x[1]) * vals[1] + sin(x[1]) * derivs[1]);
			div += 1 / (x[0] * sin(x[1])) * derivs[2];

			return div;
		}

		static Real DivCyl(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& x)
		{
			VectorN<Real, 3> vals = vectorField(x);

			VectorN<Real, 3> derivs;
			for (int i = 0; i < 3; i++)
				derivs[i] = Derivation::DeriveVecPartial<3>(vectorField, i, i, x, nullptr);

			Real div = 0.0;
			div += 1 / x[0] * (vals[0] + x[0] * derivs[0]);
			div += 1 / x[0] * derivs[1];
			div += derivs[2];

			return div;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////                     CURL                       /////////////////////////////
		static Vector3Cartesian CurlCart(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& pos)
		{
			Real dzdy = Derivation::DeriveVecPartial<3>(vectorField, 2, 1, pos, nullptr);
			Real dydz = Derivation::DeriveVecPartial<3>(vectorField, 1, 2, pos, nullptr);

			Real dxdz = Derivation::DeriveVecPartial<3>(vectorField, 0, 2, pos, nullptr);
			Real dzdx = Derivation::DeriveVecPartial<3>(vectorField, 2, 0, pos, nullptr);

			Real dydx = Derivation::DeriveVecPartial<3>(vectorField, 1, 0, pos, nullptr);
			Real dxdy = Derivation::DeriveVecPartial<3>(vectorField, 0, 1, pos, nullptr);

			Vector3Cartesian curl{ dzdy - dydz, dxdz - dzdx, dydx - dxdy };

			return curl;
		}

		static Vector3Spherical CurlSpher(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& pos)
		{
			VectorN<Real, 3> vals = vectorField(pos);

			Real dphidtheta = Derivation::DeriveVecPartial<3>(vectorField, 2, 1, pos, nullptr);
			Real dthetadphi = Derivation::DeriveVecPartial<3>(vectorField, 1, 2, pos, nullptr);

			Real drdphi = Derivation::DeriveVecPartial<3>(vectorField, 0, 2, pos, nullptr);
			Real dphidr = Derivation::DeriveVecPartial<3>(vectorField, 2, 0, pos, nullptr);

			Real dthetadr = Derivation::DeriveVecPartial<3>(vectorField, 1, 0, pos, nullptr);
			Real drdtheta = Derivation::DeriveVecPartial<3>(vectorField, 0, 1, pos, nullptr);

			Vector3Spherical ret;
			const Real& r = pos[0];
			const Real& theta = pos[1];
			const Real& phi = pos[2];

			ret[0] = 1 / (r * sin(theta)) * (cos(theta) * vals[2] + sin(theta) * dphidtheta - dthetadphi);
			ret[1] = 1 / r * (1 / sin(theta) * drdphi - vals[2] - r * dphidr);
			ret[2] = 1 / r * (vals[1] + r * dthetadr - drdtheta);

			return ret;
		}

		static Vector3Cylindrical CurlCyl(const IVectorFunction<3>& vectorField, const VectorN<Real, 3>& pos)
		{
			VectorN<Real, 3> vals = vectorField(pos);

			Real dzdphi = Derivation::DeriveVecPartial<3>(vectorField, 2, 1, pos, nullptr);
			Real dphidz = Derivation::DeriveVecPartial<3>(vectorField, 1, 2, pos, nullptr);

			Real drdz = Derivation::DeriveVecPartial<3>(vectorField, 0, 2, pos, nullptr);
			Real dzdr = Derivation::DeriveVecPartial<3>(vectorField, 2, 0, pos, nullptr);

			Real dphidr = Derivation::DeriveVecPartial<3>(vectorField, 1, 0, pos, nullptr);
			Real drdphi = Derivation::DeriveVecPartial<3>(vectorField, 0, 1, pos, nullptr);

			Vector3Cylindrical ret{ (Real)(1.0 / pos[0] * dzdphi - dphidz), drdz - dzdr, 1 / pos[0] * (vals[1] + pos[0] * dphidr - drdphi) };

			return ret;
		}
	};
}
///////////////////////////   ./include/core/CurvesSurfaces.h   ///////////////////////////



namespace MML
{
	namespace Curves2D
	{
		/////////////////////////////              CARTESIAN PLANAR CURVES                  ///////////////////////////////
		class Circle2DCurve : public IParametricCurveParametrized<2>
		{
			Real _radius;
		public:
			Circle2DCurve() : _radius(1) {}
			Circle2DCurve(Real radius) : _radius(radius) {}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			int		getNumParam() { return 1; }
			Real	getParam(int i) const { return _radius; }
			void	setParam(int i, Real val) { _radius = val; }

			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{_radius* cos(t), _radius* sin(t)}; }
		};

		class LogSpiralCurve : public IParametricCurve<2>
		{
			Real _lambda, _c;
		public:
			LogSpiralCurve() : _lambda(-1), _c(1) {}
			LogSpiralCurve(Real lambda, Real c) : _lambda(lambda), _c(c) {
				if (lambda >= 0) throw std::invalid_argument("LogSpiralCurve: lambda must be negative.");
				if (c == 0) throw std::invalid_argument("LogSpiralCurve: c must not be zero.");
			}

			Real getMinT() const { return Constants::NegativeInf; }
			Real getMaxT() const { return Constants::PositiveInf; }

			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{exp(_lambda* t)* cos(t), exp(_lambda* t)* sin(t)}; }
		};

		class LemniscateCurve : public IParametricCurve<2>
		{
		public:
			Real getMinT() const { return Constants::NegativeInf; }
			Real getMaxT() const { return Constants::PositiveInf; }

			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{cos(t) / (1 + sin(t) * sin(t)), sin(t)* cos(t) / (1 + sin(t) * sin(t))}; }
		};

		class DeltoidCurve : public IParametricCurve<2>
		{
			int _n;
		public:
			DeltoidCurve() : _n(1) {}
			DeltoidCurve(int n) : _n(n) {}

			Real getMinT() const { return Constants::NegativeInf; }
			Real getMaxT() const { return Constants::PositiveInf; }

			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{2 * _n * cos(t) * (1 + cos(t)), 2 * _n * sin(t) * (1 - cos(t))}; }
		};

		class AstroidCurve : public IParametricCurve<2>
		{
			Real _c;
		public:
			AstroidCurve() : _c(1) {}
			AstroidCurve(Real c) : _c(c) {
				if (c <= 0) throw std::invalid_argument("AstroidCurve: c must be positive.");
			}

			Real getMinT() const { return Constants::NegativeInf; }
			Real getMaxT() const { return Constants::PositiveInf; }

			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{_c* cos(t)* cos(t)* cos(t), _c* sin(t)* sin(t)* sin(t)}; }
		};

		class EpitrochoidCurve : public IParametricCurve<2>
		{
			Real _radius, _c;
			int _n;
		public:
			EpitrochoidCurve() : _radius(1), _c(1), _n(1) {}
			EpitrochoidCurve(Real radius, Real c, int n) : _radius(radius), _c(c), _n(n) {}

			Real getMinT() const { return Constants::NegativeInf; }
			Real getMaxT() const { return Constants::PositiveInf; }

			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{cos(t) - _c * cos(_n * t), sin(t) - _c * sin(_n * t) }; }
		};

		class ArchimedeanSpiralCurve : public IParametricCurve<2>
		{
			Real _a;
		public:
			ArchimedeanSpiralCurve() : _a(1) {}
			ArchimedeanSpiralCurve(Real a) : _a(a) {}

			Real getMinT() const { return Constants::NegativeInf; }
			Real getMaxT() const { return Constants::PositiveInf; }

			VectorN<Real, 2> operator()(Real t) const { return MML::VectorN<Real, 2>{_a* t* cos(t), _a* t* sin(t)}; }
		};
	} // end namespace Curves2D

	namespace Curves3D
	{
		///////////////////////////////             CARTESIAN SPACE CURVES                  ///////////////////////////////
		class LineCurve : public IParametricCurve<3>
		{
			Line3D  _line;
			Real _minT;
			Real _maxT;
		public:
			LineCurve(Real minT, Real maxT, const Point3Cartesian& pnt, const Vector3Cartesian& dir) : _line(pnt, dir), _minT(minT), _maxT(maxT) {}
			LineCurve(Real t1, const Point3Cartesian& pnt1, Real t2, const Point3Cartesian& pnt2) 
			{
				// tocno samo ako je t1 = 0.0!!!
				_line.StartPoint() = pnt1;
				Vec3Cart dir = Vec3Cart(pnt1, pnt2);
				_line.Direction() = dir.NormL2() / (t2 - t1) * dir.GetAsUnitVector();
				_minT = t1;
				_maxT = t2;
			}

			Real getMinT() const { return _minT; }
			Real getMaxT() const { return _maxT; }

			VectorN<Real, 3> operator()(Real t) const
			{
				//if (t < _minT || t > _maxT)
				//	throw std::invalid_argument("LineCurve: t is out of range.");

				auto pnt = _line.PointOnLine(t);
				return VectorN<Real, 3>{pnt.X(), pnt.Y(), pnt.Z()};
			}
		};

		// TODO - add SquarePathXY, i SquarePathRoundCorner (with radius)

		class Circle3DXY : public IParametricCurve<3> {
			Real _radius;
		public:
			Circle3DXY() : _radius(1) {}
			Circle3DXY(Real radius) : _radius(radius) {}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real t) const { return MML::VectorN<Real, 3>{_radius* cos(t), _radius* sin(t), 0}; }
		};
		class Circle3DXZ : public IParametricCurve<3> {
			Real _radius;
		public:
			Circle3DXZ() : _radius(1) {}
			Circle3DXZ(Real radius) : _radius(radius) {}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real t) const { return MML::VectorN<Real, 3>{_radius* cos(t), 0, _radius* sin(t)}; }
		};
		class Circle3DYZ : public IParametricCurve<3> {
			Real _radius;
		public:
			Circle3DYZ() : _radius(1) {}
			Circle3DYZ(Real radius) : _radius(radius) {}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real t) const { return MML::VectorN<Real, 3>{0, _radius* cos(t), _radius* sin(t)}; }
		};

		// TODO - add general circle, lying in plane with given normal, with given center position

		class HelixCurve : public IParametricCurve<3>
		{
			Real _radius, _b;
		public:
			HelixCurve() : _radius(1.0), _b(1.0) {}
			HelixCurve(Real radius, Real b) : _radius(radius), _b(b) {}

			Real getMinT() const { return Constants::NegativeInf; }
			Real getMaxT() const { return Constants::PositiveInf; }

			VectorN<Real, 3> operator()(Real t) const { return MML::VectorN<Real, 3>{_radius* cos(t), _radius* sin(t), _b* t}; }

			Real getCurvature(Real t) const { return _radius / (POW2(_radius) + POW2(_b)); }
			Real getTorsion(Real t) const { return _b / (POW2(_radius) + POW2(_b)); }
		};

		class TwistedCubicCurve : public IParametricCurve<3>
		{
		public:
			Real getMinT() const { return Constants::NegativeInf; }
			Real getMaxT() const { return Constants::PositiveInf; }

			VectorN<Real, 3> operator()(Real t) const { return MML::VectorN<Real, 3>{t, t* t, t* t* t}; }
		};

		class ToroidalSpiralCurve : public IParametricCurve<3>
		{
			int _n;
			Real _scale = 1.0;
		public:
			ToroidalSpiralCurve() : _n(1) {}
			ToroidalSpiralCurve(int n) : _n(n) {}
			ToroidalSpiralCurve(Real scale) : _scale(scale) {}
			ToroidalSpiralCurve(int n, Real scale) : _n(n), _scale(scale) {}

			Real getMinT() const { return Constants::NegativeInf; }
			Real getMaxT() const { return Constants::PositiveInf; }

			VectorN<Real, 3> operator()(Real t) const { return MML::VectorN<Real, 3>{(_scale* (4 + sin(_n * t))* cos(t)), _scale* (4 + sin(_n * t))* sin(t), _scale* cos(_n* t)}; }
		};

		///////////////////////////////             SPHERICAL SPACE CURVES                  ///////////////////////////////
		class Circle3DXYSpherical : public IParametricCurve<3> {
			Real _radius;
		public:
			Circle3DXYSpherical() : _radius(1) {}
			Circle3DXYSpherical(Real radius) : _radius(radius) {}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real t) const { return MML::VectorN<Real, 3>{_radius, Constants::PI/2, t}; }
		};

	} // end namespace Curves3D

	namespace Surfaces
	{
		// TODO - SimplePlane3D, given by point and normal

		// MonkeySaddle
		class MonkeySaddle : public IParametricSurface<3>
		{
		public:
			Real getMinU() const { return -10; }
			Real getMaxU() const { return 10; }
			Real getMinW() const { return -10; }
			Real getMaxW() const { return 10; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{u, w, u* (u* u - 3 * w * w)}; }
		};
		// MobiusStrip
		class MobiusStrip : public IParametricSurface<3>
		{
		public:
			Real getMinU() const { return 0; }
			Real getMaxU() const { return 2 * Constants::PI; }
			Real getMinW() const { return -1; }
			Real getMaxW() const { return 1; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{(1 + w * cos(u / 2))* cos(u), (1 + w * cos(u / 2))* sin(u), w* sin(u / 2)}; }
		};
		// Torus
		class Torus : public IParametricSurface<3>
		{
			Real _R, _r;
		public:
			Torus() : _R(1), _r(0.5) {}
			Torus(Real R, Real r) : _R(R), _r(r) {}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return 2 * Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{(_R + _r * cos(w))* cos(u), (_R + _r * cos(w))* sin(u), _r* sin(w)}; }
		};
		// Sphere
		class Sphere : public IParametricSurface<3>
		{
			Real _R;
		public:
			Sphere() : _R(1) {}
			Sphere(Real R) : _R(R) {}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{_R* sin(u)* cos(w), _R* sin(u)* sin(w), _R* cos(u)}; }
		};
		// Ellipsoid
		class Ellipsoid : public IParametricSurface<3>
		{
			Real _a, _b, _c;
		public:
			Ellipsoid() : _a(1), _b(1), _c(1) {}
			Ellipsoid(Real a, Real b, Real c) : _a(a), _b(b), _c(c) {}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return 2 * Constants::PI; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{_a* sin(u)* cos(w), _b* sin(u)* sin(w), _c* cos(u)}; }
		};
		// Cylinder
		class Cylinder : public IParametricSurface<3>
		{
			Real _R, _H;
		public:
			Cylinder() : _R(1), _H(1) {}
			Cylinder(Real R, Real H) : _R(R), _H(H) {}

			Real getMinU() const { return 0; }
			Real getMaxU() const { return 2 * Constants::PI; }
			Real getMinW() const { return 0; }
			Real getMaxW() const { return _H; }

			VectorN<Real, 3> operator()(Real u, Real w) const { return MML::VectorN<Real, 3>{_R* cos(u), _R* sin(u), w}; }
		};
	} // end namespace Surfaces
}

///////////////////////////   ./include/core/DiracDeltaFunction.h   ///////////////////////////


namespace MML
{
	class DiracFunction : public IRealFunction
	{
	protected:
		int _N;
	public:
		DiracFunction(int N) : _N(N) {}
	};

	class DiracStep : public DiracFunction
	{
	public:
		DiracStep(int N) : DiracFunction(N) {}

		Real operator()(const Real x) const
		{
			if (x < -1.0 / (2 * _N) || x > 1.0 / (2 * _N))
				return 0.0;
			else
				return _N;
		}
	};
	class DiracExp : public DiracFunction
	{
	public:
		DiracExp(int N) : DiracFunction(N) {}

		Real operator()(const Real x) const { return _N / sqrt(2 * Constants::PI) * exp(-x * x * _N * _N); }
	};
	class DiracSqr : public DiracFunction
	{
	public:
		DiracSqr(int N) : DiracFunction(N) {}

		Real operator()(const Real x) const { return _N / Constants::PI / (1 + _N * _N * x * x); }
	};
	class DiracSin : public DiracFunction
	{
	public:
		DiracSin(int N) : DiracFunction(N) {}

		Real operator()(const Real x) const { return sin(_N * x) / (Constants::PI * x); }
	};
}

///////////////////////////   ./include/core/ODESystem.h   ///////////////////////////



namespace MML
{
	class ODESystem : public IODESystem
	{
	protected:
		int _dim;
		void (*_func)(Real, const Vector<Real>&, Vector<Real>&);

	public:
		ODESystem() : _dim(0), _func(nullptr) { }
		ODESystem(int n, void (*inFunc)(Real, const Vector<Real>&, Vector<Real>&)) : _dim(n), _func(inFunc) { }

		int getDim() const { return _dim; }
		void derivs(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const
		{
			_func(t, x, dxdt);
		}
	};

	class ODESystemWithJacobian : public ODESystem
	{
	private:
		void (*_funcJac)(const Real, const Vector<Real>&, Vector<Real>&, Matrix<Real>&);

	public:
		ODESystemWithJacobian() : _funcJac(nullptr) { }
		ODESystemWithJacobian(int n,
			void (*inFunc)(Real, const Vector<Real>&, Vector<Real>&),
			void (*inFuncJac)(const Real t, const Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx)
		) : ODESystem(n, inFunc), _funcJac(inFuncJac) { }

		void jacobian(const Real t, Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx)
		{
			if (_funcJac != nullptr)
				_funcJac(t, x, dxdt, dydx);
		}
	};

	// calculates needed Jacobian numerically
	class ODESystemWithNumJacobian : public ODESystem
	{
	public:
		ODESystemWithNumJacobian() { }
		ODESystemWithNumJacobian(int n, void (*inFunc)(Real, const Vector<Real>&, Vector<Real>&) ) 
			: ODESystem(n, inFunc) { }

		void jacobian(const Real t, Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx)
		{
			// formirati lokalnu vektorsku funkciju na bazi derivs()
//			_funcJac(t, x, dxdt, dydx);
		}
	};

	class ODESystemSolution
	{
		// first values are initial conditions
	public:
		int  _sys_dim;
		int  _count;
		Real _x1, _x2;
		Vector<Real> _xval;
		Matrix<Real> _yval;

		ODESystemSolution() {}
		ODESystemSolution(Real x1, Real x2, int dim, int maxSteps) : _sys_dim(dim), _count(maxSteps + 1), _x1(x1), _x2(x2)
		{
			_xval.resize(maxSteps + 1);
			_yval.Resize(dim, maxSteps + 1);
		}

		template<int N>
		ParametricCurveInterpolated<N> getSolutionAsParametricCurve() const
		{
			ParametricCurveInterpolated<N> curve(_xval, _yval);
			return curve;
		}

		LinearInterpRealFunc getSolutionAsLinearInterp(int component) const
		{
			Vector<Real> xsave = _xval;
			Vector<Real> ysave = _yval.VectorFromRow(component);

			return LinearInterpRealFunc(xsave, ysave);
		}

		bool Serialize(std::string fileName, std::string title) const
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "MULTI_REAL_FUNCTION_VARIABLE_SPACED" << std::endl;

			file << title << std::endl;
			file << _sys_dim << std::endl;
			file << _count << std::endl;
			file << _x1 << std::endl;
			file << _x2 << std::endl;

			for (int i = 0; i < _count; i++)
			{
				file << _xval[i] << " ";
				for (int j = 0; j < _sys_dim; j++)
				{
					file << _yval[j][i] << " ";
				}
				file << std::endl;
			}
			file.close();
			return true;
		}
		bool SerializeAsParametricCurve3D(std::string fileName, std::string title) const
		{
			if (_sys_dim != 3)
				return false;

			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "PARAMETRIC_CURVE_CARTESIAN_3D" << std::endl;
			file << "t1: " << _x1 << std::endl;
			file << "t2: " << _x2 << std::endl;
			file << "NumPoints: " << _count << std::endl;

			for (int i = 0; i < _count; i++)
			{
				file << _xval[i] << " ";
				for (int j = 0; j < _sys_dim; j++)
				{
					file << _yval[j][i] << " ";
				}
				file << std::endl;
			}

			file.close();
			return true;
		}
	};

	class ODESystemSolutionEqualSpacing
	{
		// first values are initial conditions
	public:
		int  _sys_dim;
		int  _count;
		Real x1, x2;
		Vector<Real> xval;
		Matrix<Real> yval;

		ODESystemSolutionEqualSpacing() {}
		ODESystemSolutionEqualSpacing(int dim, int numSteps) : _sys_dim(dim), _count(numSteps + 1)
		{
			xval.resize(numSteps + 1);
			yval.Resize(dim, numSteps + 1);
		}
	};
}
///////////////////////////   ./include/core/CoordSystem.h   ///////////////////////////



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
			if (_parentFrame != nullptr)
			{
				_parentFrame = parentFrame;
				parentFrame->AddChildFrame(this);
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
			
			// TODO
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
			double latitudeDeg = localPos[0];
			double longitudeDeg = localPos[1];
			double h = localPos[2];

			// formirati sferni vektor, i vidjeti koliko se zarotira za T
			Vector3Spherical spherePos({_radius + h, Utils::DegToRad(90 - latitudeDeg), Utils::DegToRad(longitudeDeg) });
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

///////////////////////////   ./include/algorithms/PathIntegration.h   ///////////////////////////




namespace MML
{
	class PathIntegration
	{
		template<int N>
		class HelperCurveLen : public IRealFunction
		{
			const IParametricCurve<N>& _curve;
		public:
			HelperCurveLen(const IParametricCurve<N>& curve) : _curve(curve) {}

			Real operator()(Real t) const
			{
				auto tangent_vec = Derivation::DeriveCurve<N>(_curve, t, nullptr);
				return tangent_vec.NormL2();
			}
		};

		template<int N>
		class HelperLineIntegralScalarFunc : public IRealFunction
		{
			const IScalarFunction<N>& _scalar_field;
			const IParametricCurve<N>& _curve;
		public:
			HelperLineIntegralScalarFunc(const IScalarFunction<N>& scalarField, const IParametricCurve<N>& curve) : _scalar_field(scalarField), _curve(curve) {}

			Real operator()(Real t) const
			{
				auto tangent_vec = Derivation::DeriveCurve<N>(_curve, t, nullptr);

				auto field_val = _scalar_field(_curve(t));
				auto ret = field_val * tangent_vec.NormL2();

				//auto gradient = ScalarFieldOperations::GradientCart(_potential, _curve(t));
				//auto ret = tangent_vec.ScalarProductCartesian(gradient);

				return ret;
			}
		};

		template<int N>
		class HelperLineIntegralVectorFunc : public IRealFunction
		{
			const IVectorFunction<N>& _vector_field;
			const IParametricCurve<N>& _curve;
		public:
			HelperLineIntegralVectorFunc(const IVectorFunction<N>& vectorField, const IParametricCurve<N>& curve) : _vector_field(vectorField), _curve(curve) {}

			Real operator()(Real t) const
			{
				auto tangent_vec = Derivation::DeriveCurve<N>(_curve, t, nullptr);
				auto field_vec = _vector_field(_curve(t));

				return tangent_vec.ScalarProductCartesian(field_vec);
			}
		};

	public:
		template<int N>
		static Real ParametricCurveLength(const IParametricCurve<N>& curve, const Real a, const Real b)
		{
			HelperCurveLen helper(curve);

			return IntegrateTrap(helper, a, b);
		}
		// TODO - dodati i Mass, prima funkciju gustoce kao param

		static Real LineIntegral(const IScalarFunction<3>& scalarField, const IParametricCurve<3>& curve, const Real a, const Real b, const Real eps = Defaults::WorkIntegralPrecision)
		{
			HelperLineIntegralScalarFunc helper(scalarField, curve);

			return IntegrateTrap(helper, a, b, eps, nullptr);
		}

		static Real LineIntegral(const IVectorFunction<3>& vectorField, const IParametricCurve<3>& curve, const Real a, const Real b, const Real eps = Defaults::LineIntegralPrecision)
		{
			HelperLineIntegralVectorFunc helper(vectorField, curve);

			return IntegrateTrap(helper, a, b, eps, nullptr);
		}
	};
} // end namespace
///////////////////////////   ./include/algorithms/SurfaceIntegration.h   ///////////////////////////﻿#if !defined MML_SURFACE_INTEGRATION_H





namespace MML
{
	class SurfaceIntegration
	{
	public:
		static Real SurfaceIntegral(const IVectorFunction<3>& vectorField, const IParametricSurface<3>& surface, const Real x1, const Real x2, const Real y1, const Real y2)
		{
			return 0.0;
		}
		static Real SurfaceIntegral(const IVectorFunction<3>& vectorField, const SolidSurfaces3D& solid)
		{
			Real total = 0.0;
			for (int i = 0; i < solid._surfaces.size(); i++)
				total += SurfaceIntegral(vectorField, solid._surfaces[i]);

			return total;
		}
		static Real CalcSurfaceContrib(const IVectorFunction<3>& vectorField, const RectSurface3D& surface)
		{
			Real result = 0.0;
			Vector3Cartesian    normal = surface.getNormal();
			Real area = surface.getArea();
			Point3Cartesian    center = surface.getCenter();

			Vector3Cartesian    value = vectorField(Vector3Cartesian(center.X(), center.Y(), center.Z()));

			Real dotProduct = ScalarProd(normal, value);

			result = dotProduct * area;

			return result;
		}

		static Real SurfaceIntegralImproveRecursively(const IVectorFunction<3>& vectorField, const RectSurface3D& surface, Real prev_value, int level)
		{
			// now we will calculate integral for surface divided to 4 equal parts
			Real result = 0.0;
			Point3Cartesian pnt_mid12 = (surface._pnt1 + surface._pnt2) / 2.0;
			Point3Cartesian pnt_mid23 = (surface._pnt2 + surface._pnt3) / 2.0;
			Point3Cartesian pnt_mid34 = (surface._pnt3 + surface._pnt4) / 2.0;
			Point3Cartesian pnt_mid41 = (surface._pnt4 + surface._pnt1) / 2.0;
			Point3Cartesian center = surface.getCenter();

			RectSurface3D s1(surface._pnt1, pnt_mid12, center, pnt_mid41);
			RectSurface3D s2(pnt_mid12, surface._pnt2, pnt_mid23, center);
			RectSurface3D s3(center, pnt_mid23, surface._pnt3, pnt_mid34);
			RectSurface3D s4(pnt_mid41, center, pnt_mid34, surface._pnt4);

			result += CalcSurfaceContrib(vectorField, s1);
			result += CalcSurfaceContrib(vectorField, s2);
			result += CalcSurfaceContrib(vectorField, s3);
			result += CalcSurfaceContrib(vectorField, s4);

			// compare to prev_value
			if (fabs(result - prev_value) < 0.0001 || level == 0)
			{
				return result;
			}
			else
			{
				Real new_result = 0.0;

				new_result += SurfaceIntegralImproveRecursively(vectorField, s1, result, level - 1);
				new_result += SurfaceIntegralImproveRecursively(vectorField, s2, result, level - 1);
				new_result += SurfaceIntegralImproveRecursively(vectorField, s3, result, level - 1);
				new_result += SurfaceIntegralImproveRecursively(vectorField, s4, result, level - 1);

				return new_result;
			}

			return result;
		}

		static Real SurfaceIntegral(const IVectorFunction<3>& vectorField, const RectSurface3D& surface)
		{
			Real result = CalcSurfaceContrib(vectorField, surface);

			return SurfaceIntegralImproveRecursively(vectorField, surface, result, 7);
		}
	};
} // end namespace
///////////////////////////   ./include/algorithms/VolumeIntegration.h   ///////////////////////////




namespace MML
{
	class VolumeIntegration
	{

		// moze se definirati od tijela, koe ima funkciju isWithin(), i onda komponirati

	public:
		static Real VolumeIntegral(const IScalarFunction<3>& scalarField, const SolidSurfaces3D& solid)
		{
			return 0.0;
		}
	};
} // end namespace
///////////////////////////   ./include/algorithms/ODESystemSolver.h   ///////////////////////////




namespace MML
{
	class RungeKuttaSolverDumb
	{
	public:
		void rk4(Vector<Real>& y, Vector<Real>& dydx, const Real x, const Real h, Vector<Real>& yout, IODESystem& sys)
		{
			int i;
			Real xh, hh, h6;

			int n = (int)y.size();
			Vector<Real> dym(n), dyt(n), yt(n);
			hh = h * 0.5;
			h6 = h / 6.0;
			xh = x + hh;

			for (i = 0; i < n; i++) 
				yt[i] = y[i] + hh * dydx[i];
			
			sys.derivs(xh, yt, dyt);
			
			for (i = 0; i < n; i++) 
				yt[i] = y[i] + hh * dyt[i];
			
			sys.derivs(xh, yt, dym);
			
			for (i = 0; i < n; i++) {
				yt[i] = y[i] + h * dym[i];
				dym[i] += dyt[i];
			}
			
			sys.derivs(x + h, yt, dyt);
			
			for (i = 0; i < n; i++)
				yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
		}

		ODESystemSolutionEqualSpacing integrate(IODESystem& sys, const Vector<Real>& vstart, const Real x1, const Real x2, int numSteps)
		{
			int i, k;
			Real x, h;
			int dim = sys.getDim();

			ODESystemSolutionEqualSpacing sol(dim, numSteps);

			Vector<Real> v(vstart), vout(dim), dv(dim);

			for (i = 0; i < dim; i++) {
				sol.yval[i][0] = v[i];
			}
			sol.xval[0] = x1;
			x = x1;
			h = (x2 - x1) / numSteps;
			
			for (k = 0; k < numSteps; k++) 
			{
				sys.derivs(x, v, dv);

				rk4(v, dv, x, h, vout, sys);

				x += h;

				sol.xval[k + 1] = x;

				for (i = 0; i < dim; i++) {
					v[i] = vout[i];
					sol.yval[i][k + 1] = v[i];
				}
			}

			return sol;
		}
	};

	class RungeKuttaSolverSimple
	{
		void rkck(Vector<Real>& y, Vector<Real>& dydx, const Real x,
			const Real h, Vector<Real>& yout, Vector<Real>& yerr,
			IODESystem& sys)
		{
			static const Real a2 = 0.2, a3 = 0.3, a4 = 0.6, a5 = 1.0, a6 = 0.875,
				b21 = 0.2, b31 = 3.0 / 40.0, b32 = 9.0 / 40.0, b41 = 0.3, b42 = -0.9,
				b43 = 1.2, b51 = -11.0 / 54.0, b52 = 2.5, b53 = -70.0 / 27.0,
				b54 = 35.0 / 27.0, b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0,
				b63 = 575.0 / 13824.0, b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0,
				c1 = 37.0 / 378.0, c3 = 250.0 / 621.0, c4 = 125.0 / 594.0, c6 = 512.0 / 1771.0,
				dc1 = c1 - 2825.0 / 27648.0, dc3 = c3 - 18575.0 / 48384.0,
				dc4 = c4 - 13525.0 / 55296.0, dc5 = -277.00 / 14336.0, dc6 = c6 - 0.25;
			int i;

			int n = y.size();
			Vector<Real> ak2(n), ak3(n), ak4(n), ak5(n), ak6(n), ytemp(n);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + b21 * h * dydx[i];
			sys.derivs(x + a2 * h, ytemp, ak2);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (b31 * dydx[i] + b32 * ak2[i]);
			sys.derivs(x + a3 * h, ytemp, ak3);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (b41 * dydx[i] + b42 * ak2[i] + b43 * ak3[i]);
			sys.derivs(x + a4 * h, ytemp, ak4);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (b51 * dydx[i] + b52 * ak2[i] + b53 * ak3[i] + b54 * ak4[i]);
			sys.derivs(x + a5 * h, ytemp, ak5);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (b61 * dydx[i] + b62 * ak2[i] + b63 * ak3[i] + b64 * ak4[i] + b65 * ak5[i]);
			sys.derivs(x + a6 * h, ytemp, ak6);
			for (i = 0; i < n; i++)
				yout[i] = y[i] + h * (c1 * dydx[i] + c3 * ak3[i] + c4 * ak4[i] + c6 * ak6[i]);
			for (i = 0; i < n; i++)
				yerr[i] = h * (dc1 * dydx[i] + dc3 * ak3[i] + dc4 * ak4[i] + dc5 * ak5[i] + dc6 * ak6[i]);
		}

		void rkqs(Vector<Real>& y, Vector<Real>& dydx, Real& x, const Real htry,
			const Real eps, Vector<Real>& yscal, Real& hdid, Real& hnext,
			IODESystem& sys)
		{
			const Real SAFETY = 0.9, PGROW = -0.2, PSHRNK = -0.25, ERRCON = 1.89e-4;
			int i;
			Real errmax, h, htemp, xnew;

			int n = (int)y.size();
			h = htry;
			Vector<Real> yerr(n), ytemp(n);

			for (;;) {
				rkck(y, dydx, x, h, ytemp, yerr, sys);

				errmax = 0.0;
				for (i = 0; i < n; i++)
					errmax = std::max(errmax, fabs(yerr[i] / yscal[i]));
				errmax /= eps;
				if (errmax <= 1.0)
					break;

				htemp = SAFETY * h * pow(errmax, PSHRNK);
				h = (h >= Real{ 0 } ? std::max<Real>(htemp, 0.1 * h) : std::min<Real>(htemp, 0.1 * h));
				xnew = x + h;

				if (xnew == x)
					throw("stepsize underflow in rkqs");
			}
			if (errmax > ERRCON)
				hnext = SAFETY * h * pow(errmax, PGROW);
			else
				hnext = 5.0 * h;

			x += (hdid = h);

			for (i = 0; i < n; i++)
				y[i] = ytemp[i];
		}

	public:
		ODESystemSolution integrate(IODESystem& sys, const Vector<Real>& ystart, const Real x1, const Real x2, int maxSteps, Real minSaveInterval,
																const Real eps, const Real h1, const Real hmin, int& nok, int& nbad)
		{
			const int MAXSTP = 10000;
			const Real TINY = 1.0e-30;

			int dim = sys.getDim();
			ODESystemSolution sol(x1, x2, dim, maxSteps);

			int kount = 0;
			int i, nstp;
			Real xsav, x, hnext, hdid, h;
			Vector<Real> yscal(dim), y(ystart), dydx(dim);

			x = x1;
			h = SIGN(h1, x2 - x1);
			nok = nbad = kount = 0;

			if (maxSteps > 0) xsav = x - minSaveInterval * 2.0;
			
			for (nstp = 0; nstp < MAXSTP; nstp++) 
			{
				sys.derivs(x, y, dydx);
			
				for (i = 0; i < dim; i++)
					yscal[i] = fabs(y[i]) + fabs(dydx[i] * h) + TINY;
				
				if (maxSteps > 0 && kount < maxSteps - 1 && fabs(x - xsav) > fabs(minSaveInterval)) {
					for (i = 0; i < dim; i++) sol._yval[i][kount] = y[i];
					sol._xval[kount++] = x;
					xsav = x;
				}
				
				if ((x + h - x2) * (x + h - x1) > 0.0) 
					h = x2 - x;
				
				rkqs(y, dydx, x, h, eps, yscal, hdid, hnext, sys);
				
				if (hdid == h) 
					++nok; 
				else 
					++nbad;
				
				if ((x - x2) * (x2 - x1) >= 0.0) {
					if (maxSteps != 0) {
						for (i = 0; i < dim; i++) sol._yval[i][kount] = y[i];
						sol._xval[kount++] = x;
					}
					return sol;
				}
				
				if (fabs(hnext) <= hmin) 
					throw("Step size too small in integrate");
				
				h = hnext;
			}
			throw("Too many steps in routine integrate");
		}
	};
}
///////////////////////////   ./include/algorithms/ParametricCurveAnalyzer.h   ///////////////////////////




namespace MML
{
	// TODO - LOW, Curve analyzer - cuspoid i slicne tocke

	class ParametricCurveAnalyzer
	{
	public:
		// helper class that returns unit tangent vector
		template<int N>
		class CurveTangentUnit : public IParametricCurve<N>
		{
			const IParametricCurve<N>& _curve;
		public:
			CurveTangentUnit(const IParametricCurve<N>& curve) : _curve(curve) {}

			VectorN<Real, N> operator()(Real t) const
			{
				auto tangent_vec = Derivation::DeriveCurve<N>(_curve, t, nullptr);
				return tangent_vec / tangent_vec.NormL2();
			}
		};

		// Functions for calculating curve parameters
		template<int N>
		static VectorN<Real, N> getTangent(const IParametricCurve<N>& curve, Real t)
		{
			return Derivation::DeriveCurve<N>(curve, t, nullptr);
		}
		template<int N>
		static VectorN<Real, N> getTangentUnit(const IParametricCurve<N>& curve, Real t)
		{
			auto tangent = getTangent(curve, t);
			return tangent / tangent.NormL2();
		}

		template<int N>
		static VectorN<Real, N> getNormal(const IParametricCurve<N>& curve, Real t)
		{
			return Derivation::DeriveCurveSec<N>(curve, t, nullptr);
		}

		template<int N>
		static VectorN<Real, N> getNormalUnit(const IParametricCurve<N>& curve, Real t)
		{
			auto normal = getNormal(curve, t);
			return normal / normal.NormL2();
		}
		template<int N>
		static VectorN<Real, N> getPrincipalNormal(const IParametricCurve<N>& curve, Real t)
		{
			auto y_der_1 = Vector3Cartesian(getTangent(curve, t));
			auto y_der_2 = Vector3Cartesian(Derivation::DeriveCurveSec<3>(curve, t, nullptr));

			Vector3Cartesian vec_prod1 = VectorProd(y_der_2, y_der_1);
			Vector3Cartesian res_vec = VectorProd(y_der_1, vec_prod1);

			return res_vec / (y_der_1.NormL2() * vec_prod1.NormL2());
		}

		template<int N>
		static VectorN<Real, N> getBinormal(const IParametricCurve<N>& curve, Real t)
		{
			auto y_der_1 = Vector3Cartesian(getTangent(curve, t));
			auto y_der_2 = Vector3Cartesian(Derivation::DeriveCurveSec<3>(curve, t, nullptr));

			Vector3Cartesian vec_prod1 = VectorProd(y_der_2, y_der_1);
			return vec_prod1 / vec_prod1.NormL2();
		}

		template<int N>
		static VectorN<Real, N> getCurvatureVector(const IParametricCurve<N>& curve, Real t)
		{
			auto y_der_1 = getTangent(curve, t);
			auto y_der_2 = Derivation::DeriveCurveSec<3>(curve, t, nullptr);

			Real res1 = pow(y_der_1.NormL2(), -2.0);
			auto   vec2 = y_der_2 - res1 * y_der_1.ScalarProductCartesian(y_der_2) * y_der_1;

			return vec2 / res1;
		}

		template<int N>
		static Real getCurvature(const IParametricCurve<N>& curve, Real t)
		{
			auto y_der_1 = getTangent(curve, t);
			auto y_der_2 = Derivation::DeriveCurveSec<3>(curve, t, nullptr);

			Real res1 = pow(y_der_1.NormL2(), -2.0);
			auto vec2 = y_der_2 - res1 * y_der_1.ScalarProductCartesian(y_der_2) * y_der_1;
			Real res2 = vec2.NormL2();

			return res1 * res2;
		}

		static Real getCurvature3(const IParametricCurve<3>& curve, Real t)
		{
			auto curve_first_der = Vector3Cartesian(getTangent(curve, t));
			auto curve_sec_der = Vector3Cartesian(Derivation::DeriveCurveSec<3>(curve, t, nullptr));

			auto prod = VectorProd(curve_first_der, curve_sec_der);

			return prod.NormL2() / pow(curve_first_der.NormL2(), 3);
		}

		static Real getTorsion3(const IParametricCurve<3>& curve, Real t)
		{
			auto curve_first_der = Vector3Cartesian(getTangent(curve, t));
			auto curve_sec_der = Vector3Cartesian(Derivation::DeriveCurveSec<3>(curve, t, nullptr));
			auto curve_third_der = Vector3Cartesian(Derivation::DeriveCurveThird<3>(curve, t, nullptr));

			auto prod = VectorProd(curve_first_der, curve_sec_der);

			Real temp = prod.ScalarProductCartesian(curve_third_der);

			return -temp / pow(prod.NormL2(), 2);
		}

		static Plane3D getOsculationPlane(const IParametricCurve<3>& curve, Real t)
		{
			return Plane3D(Vector3Cartesian(curve(t)).getAsPoint(), Vector3Cartesian(getNormal(curve, t)));
		}
		static Plane3D getNormalPlane(const IParametricCurve<3>& curve, Real t)
		{
			return Plane3D(Vector3Cartesian(curve(t)).getAsPoint(), Vector3Cartesian(getTangentUnit(curve, t)));
		}
		static Plane3D getRectifyingPlane(const IParametricCurve<3>& curve, Real t)
		{
			return Plane3D(Vector3Cartesian(curve(t)).getAsPoint(), Vector3Cartesian(getBinormal(curve, t)));
		}

		static void getMovingTrihedron(const IParametricCurve<3>& curve, Real t, Vector3Cartesian& tangent, Vector3Cartesian& normal, Vector3Cartesian& binormal)
		{
			tangent = Vector3Cartesian(getTangentUnit(curve, t));
			normal = Vector3Cartesian(getPrincipalNormal(curve, t));
			binormal = Vector3Cartesian(getBinormal(curve, t));
		}

		static bool isArcLengthParametrized(const IParametricCurve<3>& curve, Real t1, Real t2)
		{
			int numPnt = 100;
			Real delta = (t2 - t1) / numPnt;
			for (Real t = t1 + delta; t < t2; t += delta)
			{
				Real len = PathIntegration::ParametricCurveLength(curve, t1, t);
				if (fabs(len - (t - t1)) > 1e-03)
					return false;
			}

			return true;
		}
	};
}

///////////////////////////   ./include/algorithms/RootFinding.h   ///////////////////////////﻿#if !defined MML_ROOTFINDING_H





namespace MML
{
	namespace RootFinding
	{
		// Given a function or functor func and an initial guessed range x1 to x2, the routine expands
		// the range geometrically until a root is bracketed by the returned values x1 and x2(in which
		// case zbrac returns true) or until the range becomes unacceptably large(in which case zbrac
		// returns false).
		static bool zbrac(const IRealFunction& func, double& x1, double& x2)
		{
			const int NTRY = 50;
			const double FACTOR = 1.6;

			if (x1 == x2)
				throw("Bad initial range in zbrac");

			double f1 = func(x1);
			double f2 = func(x2);

			for (int j = 0; j < NTRY; j++)
			{
				if (f1 * f2 < 0.0)
					return true;

				if (std::abs(f1) < std::abs(f2))
					f1 = func(x1 += FACTOR * (x1 - x2));
				else
					f2 = func(x2 += FACTOR * (x2 - x1));
			}
			return false;
		}

		// Given a function or functor fx defined on the interval[x1, x2], subdivide the interval into
		// n equally spaced segments, and search for zero crossings of the function.nroot will be set
		// to the number of bracketing pairs found.If it is positive, the arrays xb1[0..nroot - 1] and
		// xb2[0..nroot - 1] will be filled sequentially with any bracketing pairs that are found.On input,
		// these vectors may have any size, including zero; they will be resized to   nroot.
		static void zbrak(const IRealFunction& fx, const Real x1, const Real x2, const int n, Vector<Real>& xb1,
			Vector<Real>& xb2, int& nroot)
		{
			int nb = 20;
			xb1.resize(nb);
			xb2.resize(nb);
			nroot = 0;
			Real dx = (x2 - x1) / n;
			Real x = x1;
			Real fp = fx(x1);

			for (int i = 0; i < n; i++)
			{
				Real fc = fx(x += dx);

				if (fc * fp <= 0.0)
				{
					xb1[nroot] = x - dx;
					xb2[nroot++] = x;
					if (nroot == nb)
					{
						Vector<Real> tempvec1(xb1), tempvec2(xb2);
						xb1.resize(2 * nb);
						xb2.resize(2 * nb);
						for (int j = 0; j < nb; j++)
						{
							xb1[j] = tempvec1[j];
							xb2[j] = tempvec2[j];
						}
						nb *= 2;
					}
				}
				fp = fc;
			}
		}

		// Using bisection, return the root of a function or functor func known to lie between x1 and x2.
		// The root will be refined until its accuracy is ˙xacc.
		static Real FindRootBisection(const IRealFunction& func, const Real x1, const Real x2, const Real xacc)
		{
			const int JMAX = 50;
			Real dx, xmid, rtb;

			Real f = func(x1);
			Real fmid = func(x2);

			if (f * fmid >= 0.0)
				throw("Root must be bracketed for bisection in FindRootBisection");

			rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
			for (int j = 0; j < JMAX; j++)
			{
				fmid = func(xmid = rtb + (dx *= 0.5));

				if (fmid <= 0.0)
					rtb = xmid;

				if (std::abs(dx) < xacc || fmid == 0.0)
					return rtb;
			}
			throw("Too many bisections in FindRootBisection");
		}

		// Using the false - position method, return the root of a function or functor func known to lie
		// between x1 and x2.The root is refined until its accuracy is ˙xacc
		static Real FindRootFalsePosition(const IRealFunction& func, const Real x1, const Real x2, const Real xacc)
		{
			const int MAXIT = 30;

			Real xl, xh, del;

			Real fl = func(x1);
			Real fh = func(x2);

			if (fl * fh > 0.0)
				throw("Root must be bracketed in FindRootFalsePosition");

			if (fl < 0.0) {
				xl = x1;
				xh = x2;
			}
			else {
				xl = x2;
				xh = x1;
				std::swap(fl, fh);
			}
			Real dx = xh - xl;
			for (int j = 0; j < MAXIT; j++)
			{
				Real rtf = xl + dx * fl / (fl - fh);
				Real f = func(rtf);

				if (f < 0.0) {
					del = xl - rtf;
					xl = rtf;
					fl = f;
				}
				else {
					del = xh - rtf;
					xh = rtf;
					fh = f;
				}
				dx = xh - xl;
				if (std::abs(del) < xacc || f == 0.0)
					return rtf;
			}
			throw("Maximum number of iterations exceeded in FindRootFalsePosition");
		}

		// Using the secant method, return the root of a function or functor func thought to lie between
		// x1 and x2.The root is refined until its accuracy is ˙xacc.
		static Real FindRootSecant(const IRealFunction& func, const Real x1, const Real x2, const Real xacc) {
			const int MAXIT = 30;
			Real xl, rts;
			Real fl = func(x1);
			Real f = func(x2);
			if (std::abs(fl) < std::abs(f)) {
				rts = x1;
				xl = x2;
				std::swap(fl, f);
			}
			else {
				xl = x1;
				rts = x2;
			}
			for (int j = 0; j < MAXIT; j++)
			{
				Real dx = (xl - rts) * f / (f - fl);
				xl = rts;
				fl = f;
				rts += dx;
				f = func(rts);

				if (std::abs(dx) < xacc || f == 0.0)
					return rts;
			}
			throw("Maximum number of iterations exceeded in FindRootSecant");
		}

		// Using Ridders’ method, return the root of a function or functor func known to lie between x1
		// and x2.The root will be refined to an approximate accuracy xacc.
		static Real FindRootRidders(const IRealFunction& func, const Real x1, const Real x2, const Real xacc) {
			const int MAXIT = 60;
			Real fl = func(x1);
			Real fh = func(x2);
			if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0))
			{
				Real xl = x1;
				Real xh = x2;
				Real ans = -9.99e99;

				for (int j = 0; j < MAXIT; j++)
				{
					Real xm = 0.5 * (xl + xh);
					Real fm = func(xm);
					Real s = sqrt(fm * fm - fl * fh);

					if (s == 0.0)
						return ans;

					Real xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s);

					if (std::abs(xnew - ans) <= xacc)
						return ans;

					ans = xnew;
					Real fnew = func(ans);

					if (fnew == 0.0)
						return ans;

					if (SIGN(fm, fnew) != fm) {
						xl = xm;
						fl = fm;
						xh = ans;
						fh = fnew;
					}
					else if (SIGN(fl, fnew) != fl) {
						xh = ans;
						fh = fnew;
					}
					else if (SIGN(fh, fnew) != fh) {
						xl = ans;
						fl = fnew;
					}
					else throw("never get here.");

					if (std::abs(xh - xl) <= xacc)
						return ans;
				}
				throw("FindRootRidders exceed maximum iterations");
			}
			else {
				if (fl == 0.0)
					return x1;
				if (fh == 0.0)
					return x2;

				throw("root must be bracketed in FindRootRidders.");
			}
		}

		// Using Brent’s method, return the root of a function or functor func known to lie between x1
		// and x2.The root will be refined until its accuracy is tol.
		static Real FindRootBrent(IRealFunction& func, const Real x1, const Real x2, const Real tol)
		{
			const int ITMAX = 100;
			const Real EPS = std::numeric_limits<Real>::epsilon();
			Real a = x1, b = x2, c = x2, d, e, fa = func(a), fb = func(b), fc, p, q, r, s, tol1, xm;

			if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
				throw("Root must be bracketed in FindRootBrent");

			fc = fb;
			for (int iter = 0; iter < ITMAX; iter++)
			{
				if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
					c = a;
					fc = fa;
					e = d = b - a;
				}
				if (std::abs(fc) < std::abs(fb)) {
					a = b;
					b = c;
					c = a;
					fa = fb;
					fb = fc;
					fc = fa;
				}

				tol1 = 2.0 * EPS * std::abs(b) + 0.5 * tol;
				xm = 0.5 * (c - b);

				if (std::abs(xm) <= tol1 || fb == 0.0)
					return b;

				if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
					s = fb / fa;
					if (a == c) {
						p = 2.0 * xm * s;
						q = 1.0 - s;
					}
					else {
						q = fa / fc;
						r = fb / fc;
						p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
						q = (q - 1.0) * (r - 1.0) * (s - 1.0);
					}

					if (p > 0.0)
						q = -q;
					p = std::abs(p);

					Real min1 = 3.0 * xm * q - std::abs(tol1 * q);
					Real min2 = std::abs(e * q);

					if (2.0 * p < (min1 < min2 ? min1 : min2)) {
						e = d;
						d = p / q;
					}
					else {
						d = xm;
						e = d;
					}
				}
				else {
					d = xm;
					e = d;
				}
				a = b;
				fa = fb;
				if (std::abs(d) > tol1)
					b += d;
				else
					b += SIGN(tol1, xm);
				fb = func(b);
			}
			throw("Maximum number of iterations exceeded in FindRootBrent");
		}

		// Using the Newton-Raphson method, return the root of a function known to lie in the interval
		// x1; x2.The root will be refined until its accuracy is known within ˙xacc.
		static Real FindRootNewton(const IRealFunction& funcd, const Real x1, const Real x2, const Real xacc) {
			const int JMAX = 20;
			Real rtn = 0.5 * (x1 + x2);
			for (int j = 0; j < JMAX; j++)
			{
				Real f = funcd(rtn);
				Real df = Derivation::NDer4(funcd, rtn);
				Real dx = f / df;
				rtn -= dx;
				if ((x1 - rtn) * (rtn - x2) < 0.0)
					throw("Jumped out of brackets in rtnewt");
				if (std::abs(dx) < xacc)
					return rtn;
			}
			throw("Maximum number of iterations exceeded in rtnewt");
		}

		//Using a combination of Newton - Raphson and bisection, return the root of a function bracketed
		//between x1 and x2.The root will be refined until its accuracy is known within ˙xacc.funcd
		//is a user - supplied struct that returns the function value as a functor and the first derivative of
		//the function at the point x as the function df(see text).
		static Real FindRootSafe(const IRealFunction& funcd, const Real x1, const Real x2, const Real xacc) {
			const int MAXIT = 100;
			Real xh, xl;
			Real fl = funcd(x1);
			Real fh = funcd(x2);

			if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
				throw("Root must be bracketed in rtsafe");

			if (fl == 0.0) return x1;
			if (fh == 0.0) return x2;

			if (fl < 0.0) {
				xl = x1;
				xh = x2;
			}
			else {
				xh = x1;
				xl = x2;
			}
			Real rts = 0.5 * (x1 + x2);
			Real dxold = std::abs(x2 - x1);
			Real dx = dxold;

			Real f = funcd(rts);
			Real df = Derivation::NDer4(funcd, rts);

			for (int j = 0; j < MAXIT; j++)
			{
				if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0)
					|| (std::abs(2.0 * f) > std::abs(dxold * df))) {
					dxold = dx;
					dx = 0.5 * (xh - xl);
					rts = xl + dx;
					if (xl == rts) return rts;
				}
				else {
					dxold = dx;
					dx = f / df;
					Real temp = rts;
					rts -= dx;
					if (temp == rts) return rts;
				}
				if (std::abs(dx) < xacc)
					return rts;

				f = funcd(rts);
				df = Derivation::NDer4(funcd, rts);

				if (f < 0.0)
					xl = rts;
				else
					xh = rts;
			}
			throw("Maximum number of iterations exceeded in rtsafe");
		}
	};
}


///////////////////////////   ./include/algorithms/Statistics.h   ///////////////////////////


namespace MML
{
	class Statistics
	{
	public:
		static Real Avg(const Vector<Real>& data)
		{
			Real outAvg = 0.0;
			int n = data.size();

			for (int j = 0; j < n; j++)
				outAvg += data[j];
			outAvg /= n;

			return outAvg;
		}

		static void AvgVar(const Vector<Real>& data, Real& outAvg, Real& outVar)
		{
			Real s, ep;
			int j, n = data.size();

			outAvg = 0.0;
			for (j = 0; j < n; j++)
				outAvg += data[j];
			outAvg /= n;

			outVar = ep = 0.0;
			for (j = 0; j < n; j++) {
				s = data[j] - outAvg;
				ep += s;
				outVar += s * s;
			}
			outVar = (outVar - ep * ep / n) / (n - 1);
		}

		static void Moments(const Vector<Real>& data, Real& ave, Real& adev, Real& sdev, Real& var, Real& skew, Real& curt)
		{
			int j, n = data.size();
			Real ep = 0.0, s, p;

			if (n <= 1)
				throw("n must be at least 2 in moment");

			s = 0.0;
			for (j = 0; j < n; j++)
				s += data[j];
			ave = s / n;

			adev = var = skew = curt = 0.0;
			for (j = 0; j < n; j++) {
				adev += std::abs(s = data[j] - ave);
				ep += s;
				var += (p = s * s);
				skew += (p *= s);
				curt += (p *= s);
			}
			adev /= n;
			var = (var - ep * ep / n) / (n - 1);
			sdev = sqrt(var);

			if (var != 0.0) {
				skew /= (n * var * sdev);
				curt = curt / (n * var * var) - 3.0;
			}
			else
				throw("No skew/kurtosis when variance = 0 (in moment)");
		}
	};
}
///////////////////////////   ./include/algorithms/FunctionsAnalyzer.h   ///////////////////////////






namespace MML
{
	// TODO - MED, function point analyzer at point and interval - is continuous, is derivation defined
	class RealFunctionAnalyzer
	{
		IRealFunction& _f;
		std::string _funcName;
	public:
		RealFunctionAnalyzer(IRealFunction& f) : _f(f) {}
		RealFunctionAnalyzer(IRealFunction& f, std::string inName) : _f(f), _funcName(inName) {}

		void PrintPointAnalysis(Real x, Real eps = 1e-6)
		{
			std::cout << "Function analysis at point: " << x << ":" << std::endl;
			std::cout << "  Defined at point:    " << (isDefinedAtPoint(x) ? "yes" : "no") << std::endl;
			std::cout << "  Continuous at point: " << (isContinuousAtPoint(x, eps) ? "yes" : "no") << std::endl;
			std::cout << "  Is inflection point: " << (isInflectionPoint(x, eps) ? "yes" : "no") << std::endl;
		}

		void PrintIntervalAnalysis(Real x1, Real x2, int numPoints, Real eps = 1e-6)
		{
			if (_funcName != "")
				std::cout << std::fixed << "f(x) = " << _funcName << " - Function analysis in interval [" << x1 << ", " << x2 << "] with " << numPoints << " points:" << std::endl;
			else
				std::cout << std::fixed << "Function analysis in interval [" << x1 << ", " << x2 << "] with " << numPoints << " points:" << std::endl;

			bool isDef = true;
			bool isCont = true;
			std::vector<Real> _notDefinedPoints;
			std::vector<Real> _notContinuousPoints;

			for (int i = 0; i < numPoints; i++)
			{
				Real x = x1 + i * (x2 - x1) / numPoints;
				if (!isDefinedAtPoint(x))
				{
					isDef = false;
					_notDefinedPoints.push_back(x);
				}
				else if (!isContinuousAtPoint(x, eps))
				{
					isCont = false;
					_notContinuousPoints.push_back(x);
				}
			}

			std::cout << "  Defined    : " << (isDef ? "yes" : "no");
			if (!_notDefinedPoints.empty())
			{
				std::cout << "  Not defined at points: ";
				for (int i = 0; i < _notDefinedPoints.size(); i++)
					std::cout << _notDefinedPoints[i] << " ";
				std::cout << std::endl;
			}
			else
				std::cout << std::endl;
			std::cout << "  Continuous : " << (isCont ? "yes" : "no");
			if (!_notContinuousPoints.empty())
			{
				std::cout << "  Not continuous at points: ";
				for (int i = 0; i < _notContinuousPoints.size(); i++)
					std::cout << _notContinuousPoints[i] << " ";
				std::cout << std::endl;
			}
			else
				std::cout << std::endl;
			std::cout << "  Monotonic  : " << (isMonotonic(x1, x2, numPoints) ? "yes" : "no") << std::endl;
			std::cout << "  Min        : " << MinInNPoints(x1, x2, numPoints) << std::endl;
			std::cout << "  Max        : " << MaxInNPoints(x1, x2, numPoints) << std::endl;
		}
		void PrintDetailedIntervalAnalysis(Real x1, Real x2, int numPoints, Real eps = 1e-6)
		{
			std::cout << std::fixed << "Function analysis in interval [" << x1 << ", " << x2 << "] with " << numPoints << " points:" << std::endl;
			std::cout << " Point " << "           Value     " << "         First.der.     " << "     Sec.der.     " << "Defined " << " Continuous " << " Inflection " << std::endl;
			for (int i = 0; i < numPoints; i++)
			{
				Real x = x1 + i * (x2 - x1) / numPoints;

				std::cout << std::setw(6) << x << " :  ";
				std::cout << std::setw(16) << _f(x) << " :  ";
				std::cout << std::setw(16) << Derivation::NDer1(_f, x) << " :  ";
				std::cout << std::setw(16) << Derivation::NSecDer2(_f, x) << " :  ";
				std::cout << (isDefinedAtPoint(x) ? "   yes   " : "    no   ");
				std::cout << (isContinuousAtPoint(x, eps) ? "   yes   " : "    no   ");
				std::cout << (isInflectionPoint(x, eps) ? "     yes   " : "      no   ") << std::endl;
			}
		}

		std::vector<Real> GetRoots(Real x1, Real x2, Real eps)
		{
			std::vector<Real> roots;
			Real step = (x2 - x1) / 1000;
			Real prev = _f(x1);
			for (int i = 1; i < 1000; i++)
			{
				Real curr = _f(x1 + i * step);
				if (prev * curr < 0)
				{
					Real root = RootFinding::FindRootBisection(_f, x1 + (i - 1) * step, x1 + i * step, eps);
					roots.push_back(root);
				}
				prev = curr;
			}
			return roots;
		}
		Vector<Real> GetLocalOptimums(Real x1, Real x2)
		{
			Vector<Real> optimums;

			return optimums;
		}
		Vector<Real> GetInflectionPoints(Real x1, Real x2)
		{
			Vector<Real> inflection_points;

			return inflection_points;
		}
		
		bool isDefinedAtPoint(Real x)
		{
			Real y = _f(x);
			return !std::isnan(y) && !std::isinf(y);
		}
		bool isContinuousAtPoint(Real x, Real eps)
		{
			if (!isDefinedAtPoint(x))
				return false;

			Real h = eps;
			Real val = _f(x);
			Real left = _f(x - h);
			Real right = _f(x + h);

			// handling case of constant function
			if (val == left && val == right)
				return true;

			Real abs_dif = std::abs(left - right);

			// smanji tu abs razliku na pola, i nadji h za koji to vrijedi
			Real req_new_abs = abs_dif / 2;
			h /= 2.0;
			while (h > eps / 1000.0)
			{
				left = _f(x - h);
				right = _f(x + h);

				if (std::abs(left - right) < req_new_abs)
					return true;

				h /= 2.0;
			}
			return false;
		}
		bool isLocalOptimum(Real x, Real eps)
		{
			Real left_sec_der = Derivation::NSecDer4(_f, x - 4 * eps);
			Real right_sec_der = Derivation::NSecDer4(_f, x + 4 * eps);

			return left_sec_der * right_sec_der > 0;
		}
		bool isInflectionPoint(Real x, Real eps)
		{
			// TODO - FIX - mora biti first der == 0!
			Real left_sec_der = Derivation::NSecDer4(_f, x - 4 * eps);
			Real right_sec_der = Derivation::NSecDer4(_f, x + 4 * eps);

			return std::abs(left_sec_der - right_sec_der) < eps;
		}
		bool isContinuous(Real x1, Real x2, int numPoints)
		{
			// kroz sve tocke
			// vidjeti gdje se derivacije razlikuju u znaku (potencijalni min/max, ili discontinuity)
			// istraziti dalje oko te tocke
			return false;
		}
		bool isMonotonic(Real x1, Real x2, int numPoints)
		{
			Real step = (x2 - x1) / numPoints;
			Real prev = _f(x1 + step);
			if (_f(x1) < _f(x1 + step))
			{
				for (int i = 2; i < numPoints; i++)
				{
					Real curr = _f(x1 + i * step);
					if (curr < prev)
						return false;
					prev = curr;
				}
			}
			else
			{
				for (int i = 2; i < numPoints; i++)
				{
					Real curr = _f(x1 + i * step);
					if (curr > prev)
						return false;
					prev = curr;
				}
			}
			return true;
		}
		Real MinInNPoints(Real x1, Real x2, int numPoints)
		{
			Real step = (x2 - x1) / numPoints;
			Real min = _f(x1);

			for (int i = 1; i < numPoints; i++)
			{
				Real curr = _f(x1 + i * step);
				if (curr < min)
					min = curr;
			}
			return min;
		}
		Real MaxInNPoints(Real x1, Real x2, int numPoints)
		{
			Real step = (x2 - x1) / numPoints;
			Real max = _f(x1);

			for (int i = 1; i < numPoints; i++)
			{
				Real curr = _f(x1 + i * step);
				if (curr > max)
					max = curr;
			}
			return max;
		}
	};

	class RealFunctionComparer
	{
		IRealFunction& _f1;
		IRealFunction& _f2;

	public:
		RealFunctionComparer(IRealFunction& f1, IRealFunction& f2) : _f1(f1), _f2(f2) {}

		Real getAbsDiffSum(Real a, Real b, int numPoints)
		{
			Real step = (b - a) / numPoints;
			Real sum = 0.0;

			for (int i = 0; i < numPoints; i++)
				sum += std::abs(_f1(a + i * step) - _f2(a + i * step));

			return sum;
		}
		Real getAbsDiffAvg(Real a, Real b, int numPoints)
		{
			return getAbsDiffSum(a, b, numPoints) / numPoints;
		}
		Real getAbsDiffMax(Real a, Real b, int numPoints)
		{
			Real step = (b - a) / numPoints;
			Real max = std::abs(_f1(a) - _f2(a));

			for (int i = 0; i < numPoints; i++)
			{
				Real diff = std::abs(_f1(a + i * step) - _f2(a + i * step));
				if (diff > max)
					max = diff;
			}
			return max;
		}
		Real getRelDiffSum(Real a, Real b, int numPoints)
		{
			Real step = (b - a) / numPoints;
			Real sum = 0.0;

			for (int i = 0; i < numPoints; i++)
				if (_f1(a + i * step) != 0.0)
					sum += std::abs(_f1(a + i * step) - _f2(a + i * step)) / std::abs(_f1(a + i * step));
				else
					--numPoints;

			return sum;
		}
		Real getRelDiffAvg(Real a, Real b, int numPoints)
		{
			return getRelDiffSum(a, b, numPoints) / numPoints;
		}
		Real getRelDiffMax(Real a, Real b, int numPoints)
		{
			Real step = (b - a) / numPoints;
			Real max = 0.0;

			for (int i = 0; i < numPoints; i++)
			{
				if (_f1(a + i * step) != 0.0)
				{
					Real diff = std::abs(_f1(a + i * step) - _f2(a + i * step)) / std::abs(_f1(a + i * step));
					if (diff > max)
						max = diff;
				}
			}
			return max;
		}

		///////////                  Integration measures                /////////
		Real getIntegratedDiff(Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP)
		{
			return getIntegratedDiff(_f1, _f2, a, b, method);
		}
		static Real getIntegratedDiff(IRealFunction& f1, IRealFunction& f2, Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP)
		{
			RealFuncDiffHelper helper(f1, f2);

			switch (method)
			{
			//case IntegrationMethod::SIMPSON:
			//	return IntegrateSimpson(helper, a, b);
			//case IntegrationMethod::ROMBERG:
			//	return IntegrateRomberg(helper, a, b);
			default:
				return IntegrateTrap(helper, a, b);
			}
		}

		Real getIntegratedAbsDiff(Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP)
		{
			return getIntegratedAbsDiff(_f1, _f2, a, b, method);
		}
		static Real getIntegratedAbsDiff(IRealFunction& f1, IRealFunction& f2, Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP)
		{
			RealFuncDiffAbsHelper helper(f1, f2);

			switch (method)
			{
			//case IntegrationMethod::SIMPSON:
			//	return IntegrateSimpson(helper, a, b);
			//case IntegrationMethod::ROMBERG:
			//	return IntegrateRomberg(helper, a, b);
			default:
				return IntegrateTrap(helper, a, b);
			}
		}

		Real getIntegratedSqrDiff(Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP)
		{
			return getIntegratedSqrDiff(_f1, _f2, a, b, method);
		}
		static Real getIntegratedSqrDiff(IRealFunction& f1, IRealFunction& f2, Real a, Real b, IntegrationMethod method = IntegrationMethod::TRAP)
		{
			RealFuncDiffSqrHelper helper(f1, f2);

			switch (method)
			{
			//case IntegrationMethod::SIMPSON:
			//	return IntegrateSimpson(helper, a, b);
			//case IntegrationMethod::ROMBERG:
			//	return IntegrateRomberg(helper, a, b);
			default:
				return IntegrateTrap(helper, a, b);
			}
		}
	};

	class ScalarFieldAnalyzer
	{
		IScalarFunction<3>& _f;
	public:
		ScalarFieldAnalyzer(IScalarFunction<3>& f) : _f(f) {}

		// kao out parametar vraca "mjeru" konzervativnosti
		bool IsConservative()
		{
			return false;
		}
	};

	class VectorFieldAnalyzer
	{
		IVectorFunction<3>& _f;
	public:
		VectorFieldAnalyzer(IVectorFunction<3>& f) : _f(f) {}

		bool IsSolenoidal()
		{
			return false;
		}
	};
}
///////////////////////////   ./include/tools/ConsolePrinter.h   ///////////////////////////

namespace MML
{
	// ConsolePrinterComplex
	// - Value je bazna klasa - Real, Vector, bool, string

	template<class _Tag, class _Value>
	class TablePrinter
	{
	public:
		std::string _tagName;
		std::vector<_Tag> _tags;

		int _tagWidth, _tagPrec;
		std::vector<std::tuple<int, int, char>> _listFormatSpec;

		std::vector<std::string> _valueNames;
		std::vector<std::vector<_Value>> _values;

		TablePrinter(std::string tagName, std::vector<std::string> valueNames) :
			_tagName(tagName), _valueNames(valueNames)
		{
			_tagWidth = 8;
			_tagPrec = 3;
			for (size_t i = 0; i < _valueNames.size(); i++)
			{
				_listFormatSpec.push_back(std::make_tuple(11, 5, ' '));
			}
		}

		TablePrinter(std::string tagName, int tagWidth, int tagPrec, std::vector<std::string> valueNames, std::vector<std::tuple<int, int, char>> listWidthPrec) :
			_tagName(tagName), _tagWidth(tagWidth), _tagPrec(tagPrec), _valueNames(valueNames), _listFormatSpec(listWidthPrec)
		{}

		void addRow(_Tag tag, std::vector<_Value> values)
		{
			if (values.size() != _valueNames.size())
				throw std::invalid_argument("Number of values does not match number of value names");

			_tags.push_back(tag);
			_values.push_back(values);
		}

		void Print()
		{
			std::cout << std::setw(_tagWidth) << _tagName << " ";
			for (size_t i = 0; i < _valueNames.size(); i++)
			{
				std::cout << std::setw(std::get<0>(_listFormatSpec[i])) << std::setprecision(std::get<1>(_listFormatSpec[i])) << _valueNames[i] << " ";
			}
			std::cout << std::endl;

			for (size_t i = 0; i < _tags.size(); i++)
			{

				if (std::get<2>(_listFormatSpec[0]) == 'S')
					std::cout << std::scientific << std::setw(_tagWidth) << std::setprecision(_tagPrec) << _tags[i] << " ";
				else
					std::cout << std::fixed << std::setw(_tagWidth) << std::setprecision(_tagPrec) << _tags[i] << " ";
				for (size_t j = 0; j < _values[i].size(); j++)
				{
					if (std::get<2>(_listFormatSpec[j]) == 'S')
						std::cout << std::scientific << std::setw(std::get<0>(_listFormatSpec[j])) << std::setprecision(std::get<1>(_listFormatSpec[j])) << _values[i][j] << " ";
					else
						std::cout << std::fixed << std::setw(std::get<0>(_listFormatSpec[j])) << std::setprecision(std::get<1>(_listFormatSpec[j])) << _values[i][j] << " ";
				}
				std::cout << std::endl;
			}
		}

		void Print(int tagWidth, int tagPrec, std::vector<std::pair<int, int>> listWidthPrec)
		{
			std::cout << std::fixed << std::setw(tagWidth) << _tagName << " ";
			for (size_t i = 0; i < _valueNames.size(); i++)
			{
				std::cout << std::setw(listWidthPrec[i].first) << _valueNames[i];
			}
			std::cout << std::endl;
			for (size_t i = 0; i < _tags.size(); i++)
			{
				std::cout << std::setw(tagWidth) << std::setprecision(tagPrec) << _tags[i] << " ";
				for (size_t j = 0; j < _values[i].size(); j++)
				{
					std::cout << std::setw(listWidthPrec[j].first) << std::setprecision(listWidthPrec[j].second) << _values[i][j];
				}
				std::cout << std::endl;
			}
		}
	};
}

///////////////////////////   ./include/tools/Serializer.h   ///////////////////////////


namespace MML
{
	class Serializer
	{
	public:
		// Real function serialization 
		static bool SaveRealMultiFunc(std::vector<IRealFunction*> funcs, std::string title, Real x1, Real x2, int count, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "MULTI_REAL_FUNCTION_VARIABLE_SPACED" << std::endl;

			file << title << std::endl;
			file << funcs.size() << std::endl;
			file << count << std::endl;
			file << x1 << std::endl;
			file << x2 << std::endl;

			for (int i = 0; i < count; i++)
			{
				double x = x1 + (x2 - x1) * i / (count - 1);
				file << x << " ";
				for (int j = 0; j < funcs.size(); j++)
				{
					file << (*funcs[j])(x) << " ";
				}
				file << std::endl;
			}

			file.close();
			return true;
		}

		static bool SaveRealFuncEquallySpaced(const IRealFunction& f, std::string title, Real x1, Real x2, int numPoints, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "REAL_FUNCTION_EQUALLY_SPACED" << std::endl;
			file << title << std::endl;
			file << "x1: " << x1 << std::endl;
			file << "x2: " << x2 << std::endl;
			file << "NumPoints: " << numPoints << std::endl;

			Real step = (x2 - x1) / (numPoints - 1);
			for (int i = 0; i < numPoints; i++)
			{
				Real x = x1 + i * step;
				file << f(x) << std::endl;
			}
			file.close();
			return true;
		}

		static bool SaveRealFuncEquallySpacedDetailed(const IRealFunction& f, std::string title, Real x1, Real x2, int numPoints, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "REAL_FUNCTION_EQUALLY_SPACED_DETAILED" << std::endl;
			file << title << std::endl;
			file << "x1: " << x1 << std::endl;
			file << "x2: " << x2 << std::endl;
			file << "NumPoints: " << numPoints << std::endl;

			Real step = (x2 - x1) / (numPoints - 1);
			for (int i = 0; i < numPoints; i++)
			{
				Real x = x1 + i * step;
				file << x << " " << f(x) << std::endl;
			}
			file.close();
			return true;
		}

		static bool SaveRealFuncVariableSpaced(const IRealFunction& f, std::string title, Vector<Real> points, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "REAL_FUNCTION_VARIABLE_SPACED" << std::endl;
			file << title << std::endl;

			for (int i = 0; i < points.size(); i++)
			{
				Real x = points[i];
				file << x << " " << f(x) << std::endl;
			}
			file.close();
			return true;
		}

		// Parametric curve serialization
		template<int N>
		static bool SaveParamCurve(const IRealToVectorFunction<N>& f, std::string inType, std::string title, Real t1, Real t2, int numPoints, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << inType << std::endl;
			file << title << std::endl;
			file << "t1: " << t1 << std::endl;
			file << "t2: " << t2 << std::endl;
			file << "NumPoints: " << numPoints << std::endl;

			Real delta = (t2 - t1) / (numPoints - 1);
			for (Real t = t1; t <= t2; t += delta)
			{
				file << t << " ";
				for (int i = 0; i < N; i++)
					file << f(t)[i] << " ";
				file << std::endl;
			}
			file.close();
			return true;
		}

		template<int N>
		static bool SaveParamCurve(const IRealToVectorFunction<N>& f, std::string inType, std::string title, Vector<Real> points, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << inType << std::endl;
			file << title << std::endl;
			file << "t1: " << points[0] << std::endl;
			file << "t2: " << points[points.size() - 1] << std::endl;
			file << "NumPoints: " << points.size() << std::endl;
			for (int i = 0; i < points.size(); i++)
			{
				Real t = points[i];
				file << t << " ";
				for (int i = 0; i < N; i++)
					file << f(t)[i] << " ";
				file << std::endl;
			}

			file.close();
			return true;
		}

		template<int N>
		static bool SaveAsParamCurve(std::vector<VectorN<Real, N>> vals, std::string inType, std::string title, Real t1, Real t2, int numPoints, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << inType << std::endl;
			file << title << std::endl;
			file << "t1: " << t1 << std::endl;
			file << "t2: " << t2 << std::endl;
			file << "NumPoints: " << numPoints << std::endl;

			Real delta = (t2 - t1) / (numPoints - 1);
			for (int i = 0; i < numPoints; i++)
			{
				Real t = t1 + i * delta;
				file << t << " ";
				for (int j = 0; j < N; j++)
					file << vals[i][j] << " ";
				file << std::endl;
			}
			file.close();
			return true;

		}

		//bool SerializeCartesian2D(Real t1, Real t2, int numPoints, std::string fileName) const
		//{
		//    return SaveParamCurve<2>("PARAMETRIC_CURVE_CARTESIAN_2D", t1, t2, numPoints, fileName);
		//}
		//bool SerializeCartesian2DAtPoints(Vector<Real> points, std::string fileName) const
		//{
		//    return SaveParamCurveAtPoints<2>("PARAMETRIC_CURVE_CARTESIAN_2D_AT_POINTS", points, fileName);
		//}
		static bool SaveParamCurveCartesian3D(const IRealToVectorFunction<3>& f, std::string title, Real t1, Real t2, int numPoints, std::string fileName)
		{
			return SaveParamCurve<3>(f, "PARAMETRIC_CURVE_CARTESIAN_3D", title, t1, t2, numPoints, fileName);
		}
		//bool SerializeCartesian3DAtPoints(Vector<Real> points, std::string fileName) const
		//{
		//    return SaveParamCurveAtPoints("PARAMETRIC_CURVE_CARTESIAN_3D", points, fileName);
		//}
		//bool SerializePolar(Real t1, Real t2, int numPoints, std::string fileName) const
		//{
		//    return SaveParamCurve("PARAMETRIC_CURVE_POLAR", t1, t2, numPoints, fileName);
		//}
		//bool SerializePolarAtPoints(Vector<Real> points, std::string fileName) const
		//{
		//    return SaveParamCurveAtPoints("PARAMETRIC_CURVE_POLAR", points, fileName);
		//}
		//bool SerializeSpherical(Real t1, Real t2, int numPoints, std::string fileName) const
		//{
		//    return SaveParamCurve("PARAMETRIC_CURVE_SPHERICAL", t1, t2, numPoints, fileName);
		//}
		//bool SerializeSphericalAtPoints(Vector<Real> points, std::string fileName) const
		//{
		//    return SaveParamCurveAtPoints("PARAMETRIC_CURVE_SPHERICAL", points, fileName);
		//}

		// Scalar function serialization
		static bool SaveScalarFunc2DCartesian(const IScalarFunction<2>& f, std::string title, Real x1, Real x2, int numPointsX, Real y1, Real y2, int numPointsY, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "SCALAR_FUNCTION_CARTESIAN_2D" << std::endl;
			file << title << std::endl;
			file << "x1: " << x1 << std::endl;
			file << "x2: " << x2 << std::endl;
			file << "NumPointsX: " << numPointsX << std::endl;
			file << "y1: " << y1 << std::endl;
			file << "y2: " << y2 << std::endl;
			file << "NumPointsY: " << numPointsY << std::endl;

			Real stepX = (x2 - x1) / (numPointsX - 1);
			Real stepY = (y2 - y1) / (numPointsY - 1);
			for (int i = 0; i < numPointsX; i++)
			{
				for (int j = 0; j < numPointsY; j++)
				{
					Real x = x1 + i * stepX;
					Real y = y1 + j * stepY;
					file << x << " " << y << " " << f(VectorN<Real, 2>{x, y}) << std::endl;
				}
			}
			return true;
		}

		static bool SaveScalarFunc3DCartesian(const IScalarFunction<3>& f, std::string title, Real x1, Real x2, int numPointsX, Real y1, Real y2, int numPointsY, Real z1, Real z2, int numPointsZ, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << "SCALAR_FUNCTION_CARTESIAN_3D" << std::endl;
			file << title << std::endl;
			file << "x1: " << x1 << std::endl;
			file << "x2: " << x2 << std::endl;
			file << "NumPointsX: " << numPointsX << std::endl;
			file << "y1: " << y1 << std::endl;
			file << "y2: " << y2 << std::endl;
			file << "NumPointsY: " << numPointsY << std::endl;
			file << "z1: " << z1 << std::endl;
			file << "z2: " << z2 << std::endl;
			file << "NumPointsZ: " << numPointsZ << std::endl;

			Real stepX = (x2 - x1) / (numPointsX - 1);
			Real stepY = (y2 - y1) / (numPointsY - 1);
			Real stepZ = (z2 - z1) / (numPointsZ - 1);
			for (int i = 0; i < numPointsX; i++)
			{
				for (int j = 0; j < numPointsY; j++)
				{
					for (int k = 0; k < numPointsZ; k++)
					{
						Real x = x1 + i * stepX;
						Real y = y1 + j * stepY;
						Real z = z1 + j * stepZ;
						file << x << " " << y << " " << z << " " << f(VectorN<Real, 3>{x, y, z}) << std::endl;
					}
				}
			}
			return true;
		}

		// vector function serialization
		static bool SaveVectorFunc3D(const IVectorFunction<3>& f, std::string inType, std::string title, Real x1_start, Real x1_end, int numPointsX1, Real x2_start, Real x2_end, int numPointsX2, Real x3_start, Real x3_end, int numPointsX3, std::string fileName)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << inType << std::endl;
			file << title << std::endl;

			Real stepX = (x1_end - x1_start) / (numPointsX1 - 1);
			Real stepY = (x2_end - x2_start) / (numPointsX2 - 1);
			Real stepZ = (x3_end - x3_start) / (numPointsX3 - 1);
			for (int i = 0; i < numPointsX1; i++)
				for (int j = 0; j < numPointsX2; j++)
					for (int k = 0; k < numPointsX3; k++)
					{
						Real x = x1_start + i * stepX;
						Real y = x2_start + j * stepY;
						Real z = x3_start + k * stepZ;
						auto val = f(VectorN<Real, 3>{x, y, z});
						file << x << " " << y << " " << z << " " << val[0] << " " << val[1] << " " << val[2] << std::endl;
					}

			file.close();
			return true;
		}

		static bool SaveVectorFunc3D(const IVectorFunction<3>& f, std::string inType, std::string title, Real x1_start, Real x1_end, int numPointsX1, Real x2_start, Real x2_end, int numPointsX2, Real x3_start, Real x3_end, int numPointsX3, std::string fileName, Real upper_threshold)
		{
			std::ofstream file(fileName);
			if (!file.is_open())
				return false;

			file << inType << std::endl;
			file << title << std::endl;

			Real stepX = (x1_end - x1_start) / (numPointsX1 - 1);
			Real stepY = (x2_end - x2_start) / (numPointsX2 - 1);
			Real stepZ = (x3_end - x3_start) / (numPointsX3 - 1);
			for (int i = 0; i < numPointsX1; i++)
				for (int j = 0; j < numPointsX2; j++)
					for (int k = 0; k < numPointsX3; k++)
					{
						Real x = x1_start + i * stepX;
						Real y = x2_start + j * stepY;
						Real z = x3_start + k * stepZ;
						auto val = f(VectorN<Real, 3>{x, y, z});

						if (val.NormL2() < upper_threshold)
							file << x << " " << y << " " << z << " " << val[0] << " " << val[1] << " " << val[2] << std::endl;
					}

			file.close();
			return true;
		}

		static bool SaveVectorFunc3DCartesian(const IVectorFunction<3>& f, std::string title, Real x1, Real x2, int numPointsX, Real y1, Real y2, int numPointsY, Real z1, Real z2, int numPointsZ, std::string fileName)
		{
			return SaveVectorFunc3D(f, "VECTOR_FIELD_3D_CARTESIAN", title, x1, x2, numPointsX, y1, y2, numPointsY, z1, z2, numPointsZ, fileName);
		}
		static bool SaveVectorFunc3DCartesian(const IVectorFunction<3>& f, std::string title, Real x1, Real x2, int numPointsX, Real y1, Real y2, int numPointsY, Real z1, Real z2, int numPointsZ, std::string fileName, Real upper_threshold)
		{
			return SaveVectorFunc3D(f, "VECTOR_FIELD_3D_CARTESIAN", title, x1, x2, numPointsX, y1, y2, numPointsY, z1, z2, numPointsZ, fileName, upper_threshold);
		}
		static bool SaveVectorFuncSpherical(const IVectorFunction<3>& f, std::string title, Real r1, Real r2, int numPointsR, Real theta1, Real theta2, int numPointsTheta, Real phi1, Real phi2, int numPointsPhi, std::string fileName)
		{
			return SaveVectorFunc3D(f, "VECTOR_FIELD_SPHERICAL", title, r1, r2, numPointsR, theta1, theta2, numPointsTheta, phi1, phi2, numPointsPhi, fileName);
		}
		static bool SaveVectorFuncSpherical(const IVectorFunction<3>& f, std::string title, Real r1, Real r2, int numPointsR, Real theta1, Real theta2, int numPointsTheta, Real phi1, Real phi2, int numPointsPhi, std::string fileName, Real upper_threshold)
		{
			return SaveVectorFunc3D(f, "VECTOR_FIELD_SPHERICAL", title, r1, r2, numPointsR, theta1, theta2, numPointsTheta, phi1, phi2, numPointsPhi, fileName, upper_threshold);
		}
	};
}
///////////////////////////   ./include/tools/Visualizer.h   ///////////////////////////



namespace MML
{
	class Visualizer
	{
		static inline std::string _pathResultFiles{ GLOB_PATH_ResultFiles };

		static inline std::string _pathRealFuncViz{ GLOB_PATH_RealFuncViz };
		static inline std::string _pathSurfaceViz{ GLOB_PATH_SurfaceViz };
		static inline std::string _pathParametricCurveViz{ GLOB_PATH_ParametricCurveViz };
		static inline std::string _pathVectorFieldViz{ GLOB_PATH_VectorFieldViz };

	public:
		static void VisualizeRealFunction(const IRealFunction& f, std::string title, Real x1, Real x2, int numPoints, std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			Serializer::SaveRealFuncEquallySpacedDetailed(f, title, x1, x2, numPoints, name);

			std::string command = _pathRealFuncViz + " " + name;
			system(command.c_str());
			std::cout << "VisualizeRealFunction: Not implemented for this OS" << std::endl;
		}

		static void VisualizeMultiRealFunction(std::vector<IRealFunction*> funcs, std::string title, Real x1, Real x2, int numPoints, std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			Serializer::SaveRealMultiFunc(funcs, title, x1, x2, numPoints, name);

			std::string command = _pathRealFuncViz + " " + name;
			system(command.c_str());
			std::cout << "VisualizeMultiRealFunction: Not implemented for this OS" << std::endl;
		}

		static void VisualizeScalarFunc2DCartesian(const IScalarFunction<2>& func, std::string title, Real x1, Real x2, int numPointsX, Real y1, Real y2, int numPointsY, std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			Serializer::SaveScalarFunc2DCartesian(func, title, x1, x2, numPointsX, y1, y2, numPointsY, name);

			std::string command = _pathSurfaceViz + " " + name;
			system(command.c_str());
			std::cout << "VisualizeScalarFunc2DCartesian: Not implemented for this OS" << std::endl;
		}

		static void VisualizeVectorField3DCartesian(const IVectorFunction<3>& func, std::string title, Real x1, Real x2, int numPointsX, Real y1, Real y2, int numPointsY, Real z1, Real z2, int numPointsZ, std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			Serializer::SaveVectorFunc3DCartesian(func, title, x1, x2, numPointsX, y1, y2, numPointsY, z1, z2, numPointsZ, name);

			std::string command = _pathVectorFieldViz + " " + name;
			system(command.c_str());
			std::cout << "VisualizeVectorField3DCartesian: Not implemented for this OS" << std::endl;
		}

		static void VisualizeParamCurve3D(const IRealToVectorFunction<3>& f, std::string title, Real t1, Real t2, int numPoints, std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			Serializer::SaveParamCurveCartesian3D(f, title, t1, t2, numPoints, name);

			std::string command = _pathParametricCurveViz + " " + name;
			system(command.c_str());
			std::cout << "VisualizeParamCurve3D: Not implemented for this OS" << std::endl;
		}

		static void VisualizeMultiParamCurve3D(std::vector<std::string> fileNames)
		{
			std::string params;
			for (auto& name : fileNames)
			{
				params = params + (_pathResultFiles + name + " ");
			}

			std::string command = _pathParametricCurveViz + " " + params;
			system(command.c_str());
			std::cout << "VisualizeMultiParamCurve3D: Not implemented for this OS" << std::endl;
		}

		static void VisualizeODESysSolAsMultiFunc(const ODESystemSolution& sol, std::string title, std::string fileName)
		{
			std::string name = _pathResultFiles + fileName;
			sol.Serialize(name, title);

			std::string command = _pathRealFuncViz + " " + name;
			system(command.c_str());
			std::cout << "VisualizeODESysSolAsMultiFunc: Not implemented for this OS" << std::endl;
		}

		static void VisualizeODESysSolAsParamCurve3(const ODESystemSolution& sol, std::string title, std::string fileName)
		{
			if (sol._sys_dim != 3)
				throw std::runtime_error("VisualizeODESysSolAsParamCurve3: system dimension must be 3");

			std::string name = _pathResultFiles + fileName;
			sol.SerializeAsParametricCurve3D(name, title);

			std::string command = _pathParametricCurveViz + " " + name;
			system(command.c_str());
			std::cout << "VisualizeODESysSolAsParamCurve3: Not implemented for this OS" << std::endl;
		}
	};
}

#endif
