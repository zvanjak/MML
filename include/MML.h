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
#  if (defined(TARGET_OS_OSX) && TARGET_OS_OSX == 1) || \
      (defined(TARGET_OS_MAC) && TARGET_OS_MAC == 1)
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



// Complex must have the same underlaying type as Real
typedef double               Real;      // default real type
typedef std::complex<double> Complex;   // default complex type

// Global paths for Visualizers
static const std::string GLOB_PATH_ResultFiles = "E:\\Projects\\MML\\results\\";
static const std::string GLOB_PATH_RealFuncViz = "E:\\Projects\\MML\\tools\\visualizers\\real_function_visualizer\\MML_RealFunctionVisualizer.exe";
static const std::string GLOB_PATH_SurfaceViz = "E:\\Projects\\MML\\tools\\visualizers\\scalar_function_2d_visualizer\\MML_ScalarFunction2Visualizer.exe";
static const std::string GLOB_PATH_ParametricCurveViz = "E:\\Projects\\MML\\tools\\visualizers\\parametric_curve_visualizer\\MML_ParametricCurveVisualizer.exe";
static const std::string GLOB_PATH_VectorFieldViz = "E:\\Projects\\MML\\tools\\visualizers\\vector_field_visualizer\\MML_VectorFieldVisualizer.exe";


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

	template<class T> inline T POW2(const T a) { return ((a) * (a)); }
	template<class T> inline T POW3(const T a) { return ((a) * (a) * (a)); }
	template<class T> inline T POW4(const T a) { return ((a) * (a) * (a) * (a)); }

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
		//////////               Default precisions             ///////////
		// TODO - make dependent on Real type (ie. different values for float, double and long double)
		static inline const double ComplexEqualityPrecision = 1e-15;
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

	//////////             Vector error exceptions            ///////////
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
///////////////////////////   ./include/interfaces/IAlgebra.h   ///////////////////////////

namespace MML
{
    // TODO - Monoid, Ring, CommutativeRing, CommutativeGroup
    template<class _ElemType>
    class IGroup
    {
        public:
        virtual _ElemType  identity() const = 0;
        virtual _ElemType  inverse(const _ElemType& a) const = 0;
        virtual _ElemType  op(const _ElemType& a, const _ElemType& b) const = 0;
    };

    template<class _ElemType>
    class Field
    {
        public:
        virtual _ElemType  identity() = 0;
        virtual _ElemType  zero() = 0;
        virtual _ElemType  inverse(const _ElemType& a) = 0;
        virtual _ElemType  add(const _ElemType& a, const _ElemType& b) = 0;
        virtual _ElemType  mul(const _ElemType& a, const _ElemType& b) = 0;
    };
}
///////////////////////////   ./include/interfaces/IVectorSpace.h   ///////////////////////////

namespace MML
{
	template<class _FieldType, class _VecType>
	class VectorSpace
	{
	public:
		virtual _FieldType  zero() const { return _FieldType(0); }
		virtual _FieldType  inverse(const _FieldType& b) const { return -b; }

		virtual _VecType    add(const _VecType& a, const _VecType& b) const { return a + b; }
		virtual _VecType    mul(const _FieldType& a, const _VecType& b) const { return a * b; }
	};

	template<class _FieldType, class _VecType>
	class NormedVectorSpace : public VectorSpace<_FieldType, _VecType>
	{
	public:
		virtual Real  norm(const _VecType& a) const = 0;
	};

	template<class _FieldType, class _VecType>
	class HilbertSpace : public NormedVectorSpace<_FieldType, _VecType>
	{
	public:
		virtual Real  norm(const _VecType& a) const
		{
			return sqrt(scal_prod(a, a));
		}
		virtual Real  scal_prod(const _VecType& a, const _VecType& b) const = 0;
	};

	template<class _VecTo, class _VecFrom>
	class ILinearOperator
	{
	public:
		virtual _VecTo  operator()(const _VecFrom& x) const = 0;
	};
}
///////////////////////////   ./include/interfaces/ITensor.h   ///////////////////////////

namespace MML
{
	template<int N>
	class ITensor2
	{
	public:
		virtual int   NumContravar() const = 0;
		virtual int   NumCovar() const = 0;

		virtual Real  Component(int i, int j) const = 0;
		virtual Real& Component(int i, int j) = 0;
	};

	template<int N>
	class ITensor3
	{
	public:
		virtual int   NumContravar() const = 0;
		virtual int   NumCovar() const = 0;

		virtual Real  Component(int i, int j, int k) const = 0;
		virtual Real& Component(int i, int j, int k) = 0;
	};

	template<int N>
	class ITensor4
	{
	public:
		virtual int   NumContravar() const = 0;
		virtual int   NumCovar() const = 0;

		virtual Real  Component(int i, int j, int k, int l) const = 0;
		virtual Real& Component(int i, int j, int k, int l) = 0;
	};

	template<int N>
	class ITensor5
	{
	public:
		virtual int   NumContravar() const = 0;
		virtual int   NumCovar() const = 0;

		virtual Real  Component(int i, int j, int k, int l, int m) const = 0;
		virtual Real& Component(int i, int j, int k, int l, int m) = 0;
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

		void GetEquidistantCovering(int numPoints) { }
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

///////////////////////////   ./include/base/Vector.h   ///////////////////////////

namespace MML
{
	template<class Type>
	class Vector
	{
	private:
		std::vector<Type> _elems;

	public:
		///////////////////////          Constructors and destructor       //////////////////////
		Vector() {}
		Vector(size_t n) : _elems(n) {}
		Vector(size_t n, Type val) : _elems(n, val) {}
		Vector(std::vector<Type> values) : _elems(values) {}
		Vector(std::initializer_list<Type> list) : _elems(list) {}
		Vector(int n, Type* vals) : _elems(n)
		{
			for (int i = 0; i < n; ++i)
				_elems[i] = vals[i];
		}

		void Resize(int newLen)
		{
			_elems.resize(newLen);
		}

		////////////////////////            Standard stuff             ////////////////////////
		int  size() const { return (int)_elems.size(); }

		static Vector GetUnitVector(int dimVec, int indUnit)
		{
			if (indUnit < 0 || indUnit >= dimVec)
				throw VectorDimensionError("Vector::GetUnitVector - wrong unit index", dimVec, indUnit);

			Vector ret(dimVec);
			ret[indUnit] = Type{ 1.0 };
			return ret;
		}

		bool IsEqual(const Vector& b, Real eps = Defaults::VectorEqualityPrecision) const
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
		static bool AreEqual(const Vector& a, const Vector& b, Type eps = Defaults::VectorEqualityPrecision)
		{
			return a.IsEqual(b, eps);
		}

		///////////////////////////            Operators             ///////////////////////////
		Type& operator[](int n) { return _elems[n]; }
		const Type& operator[](int n) const { return _elems[n]; }

		Vector operator-()          // unary minus
		{
			Vector ret(size());;
			for (int i = 0; i < size(); i++)
				ret._elems[i] = -(*this)[i];
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
		bool operator==(const Vector& b) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::operator==() - vectors must be equal size", size(), b.size());

			for (int i = 0; i < size(); i++)
			{
				if ((*this)[i] != b[i])
					return false;
			}
			return true;
		}

		friend Vector operator*(Type a, const Vector& b)
		{
			Vector ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret._elems[i] = a * b._elems[i];
			return ret;
		}
		friend Vector operator*(const Vector& a, Type b)
		{
			Vector ret(a.size());;
			for (int i = 0; i < a.size(); i++)
				ret._elems[i] = b * a._elems[i];
			return ret;
		}
		friend Vector operator/(const Vector& a, Type b)
		{
			Vector ret(a.size());
			for (int i = 0; i < a.size(); i++)
				ret._elems[i] = a._elems[i] / b;
			return ret;
		}

		//////////////////////                 Operations                 ///////////////////////
		Type ScalarProductCartesian(const Vector& b) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::ScalarProductCartesian - vectors must be equal size", size(), b.size());

			Type product = 0.0;
			for (int i = 0; i < size(); i++)
				product += (*this)[i] * b[i];
			return product;
		}
		Type NormL2() const
		{
			Type norm = 0.0;
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
		friend std::ostream& operator<<(std::ostream& stream, const Vector& a)
		{
			a.Print(stream, 15, 10);

			return stream;
		}
	};

	typedef Vector<int>     VectorInt;
	typedef Vector<float>   VectorFlt;
	typedef Vector<Real>    VectorDbl;
	typedef Vector<Complex> VectorComplex;

	typedef Vector<int>     VecI;
	typedef Vector<float>   VecF;
	typedef Vector<Real>    VecD;
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
		VectorN(const Type& init_val) {
			for (int i = 0; i < N; ++i)
				_val[i] = init_val;
		}
		VectorN(std::initializer_list<Type> list)
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
		VectorN(std::vector<Type> list)
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
		VectorN(Type* vals)
		{
			for (int i = 0; i < N; ++i)
				_val[i] = vals[i];
		}

		////////////////////////            Standard stuff             ////////////////////////
		int size() const { return N; }

		static VectorN GetUnitVector(int indUnit)
		{
			VectorN ret;
			ret[indUnit] = 1.0;
			return ret;
		}

		VectorN GetAsUnitVector() const
		{
			return VectorN{ (*this) / NormL2() };
		}
		VectorN GetAsUnitVectorAtPos(const VectorN& pos) const
		{
			return VectorN{ (*this) / NormL2() };
		}

		bool IsEqual(const VectorN& b, Type eps = Defaults::VectorEqualityPrecision) const
		{
			for (int i = 0; i < N; i++)
			{
				if (fabs((*this)[i] - b[i]) > eps)
					return false;
			}
			return true;
		}
		static bool AreEqual(const VectorN& a, const VectorN& b, Type eps = Defaults::VectorEqualityPrecision)
		{
			return a.IsEqual(b, eps);
		}

		///////////////////////////            Operators             ///////////////////////////
		Type& operator[](int n) { return _val[n]; }
		Type  operator[](int n) const { return _val[n]; }

		VectorN operator+(const VectorN& b) const
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = _val[i] + b._val[i];
			return ret;
		}
		VectorN operator-()         // unary minus
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = -_val[i];
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
			{
				if ((*this)[i] != b[i])
					return false;
			}
			return true;
		}

		friend VectorN operator*(const VectorN& a, Type b)
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = a[i] * b;
			return ret;
		}
		friend VectorN operator*(Type a, const VectorN& b)
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = a * b[i];
			return ret;
		}
		friend VectorN operator/(const VectorN& a, Type b)
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = a[i] / b;
			return ret;
		}

		//////////////////////                 Operations                 ///////////////////////
		Type ScalarProductCartesian(const VectorN& b) const
		{
			Type product = 0.0;
			for (int i = 0; i < N; i++)
				product += (*this)[i] * b[i];
			return product;
		}
		Type NormL2() const
		{
			Type norm = 0.0;
			for (int i = 0; i < size(); i++)
				norm += (*this)[i] * (*this)[i];
			return std::sqrt(norm);
		}
		Type AngleToVector(const VectorN& b) const
		{
			if (size() != b.size())
				throw VectorDimensionError("VectorN::AngleToVector - vectors must be equal size", size(), b.size());

			Type cosAngle = this->ScalarProductCartesian(b) / (NormL2() * b.NormL2());
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

        if( Abs(x) > zeroThreshold )
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

	typedef VectorN<float, 2> Vector2Flt;
	typedef VectorN<float, 3> Vector3Flt;
	typedef VectorN<float, 4> Vector4Flt;

	typedef VectorN<double, 2> Vector2Dbl;
	typedef VectorN<double, 3> Vector3Dbl;
	typedef VectorN<double, 4> Vector4Dbl;

	typedef VectorN<Complex, 2> Vector2Complex;
	typedef VectorN<Complex, 3> Vector3Complex;
	typedef VectorN<Complex, 4> Vector4Complex;

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
//#include <format>



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

				for (int i = 1; i < rows; i++)
					_data[i] = _data[i - 1] + cols;
			}
			else
				throw MatrixAllocationError("Matrix::Init - allocation error", rows, cols);
		}

	public:
		///////////////////////          Constructors and destructor       //////////////////////
		Matrix() : _rows(0), _cols(0), _data{ nullptr } {}
		Matrix(int rows, int cols) : _rows(rows), _cols(cols)
		{
			Init(rows, cols);
			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					_data[i][j] = 0;
		}
		Matrix(int rows, int cols, Type val) : _rows(rows), _cols(cols)
		{
			Init(rows, cols);
			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					_data[i][j] = val;
		}
		// useful if you have a pointer to continuous (row-wise) 2D array
		Matrix(int rows, int cols, Type* val) : _rows(rows), _cols(cols)
		{
			Init(rows, cols);
			for (int i = 0; i < _rows; ++i)
				for (int j = 0; j < _cols; ++j)
					_data[i][j] = *val++;
		}
		// in strict mode, you must supply ALL necessary values for complete matrix initialization
		Matrix(int rows, int cols, std::initializer_list<Type> values, bool strictMode = true) : _rows(rows), _cols(cols)
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
			if (rows <= 0 || cols < 0)
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

		typedef Type value_type;      // make T available externally

		////////////////////////            Standard stuff             ////////////////////////
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

		/////////////////////          Matrix to Vector conversions           ////////////////////
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

		/////////////////////          Matrix properties            ////////////////////
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

		bool IsDiagonal(double eps = Defaults::IsMatrixDiagonalPrecision) const
		{
			for (int i = 0; i < RowNum(); i++)
				for (int j = 0; j < ColNum(); j++)
					if (i != j && Abs((*this)[i][j]) > eps)
						return false;

			return true;
		}

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

		/////////////////////            Assignment operators             ////////////////////
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
		Matrix& operator=(Matrix&& m)
		{
			if (this == &m)
				return *this;

			std::swap(_data, m._data);
			std::swap(_rows, m._rows);
			std::swap(_cols, m._cols);

			return *this;
		}

		////////////////////            Access operators             ///////////////////////
		Type* operator[](int i) { return _data[i]; }
		const Type* operator[](const int i) const { return _data[i]; }

		Type  operator()(int i, int j) const { return _data[i][j]; }
		Type& operator()(int i, int j) { return _data[i][j]; }

		// version with checking bounds
		Type  ElemAt(int i, int j) const
		{
			if (i < 0 || i >= RowNum() || j < 0 || j >= ColNum())
				throw MatrixAccessBoundsError("Matrix::ElemAt", i, j, RowNum(), ColNum());

			return _data[i][j];
		}
		Type& ElemAt(int i, int j)
		{
			if (i < 0 || i >= RowNum() || j < 0 || j >= ColNum())
				throw MatrixAccessBoundsError("Matrix::ElemAt", i, j, RowNum(), ColNum());

			return _data[i][j];
		}

		////////////////////            Arithmetic operators             ////////////////////
		Matrix operator-()             // unary minus
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
			for (size_t i = 0; i < _rows; i++)
				for (size_t j = 0; j < _cols; j++)
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

		friend Matrix operator*(const Matrix<Type>& a, Type b)
		{
			int	i, j;
			Matrix	ret(a.RowNum(), a.ColNum());

			for (i = 0; i < a.RowNum(); i++)
				for (j = 0; j < a.ColNum(); j++)
					ret[i][j] = a._data[i][j] * b;

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
		friend Matrix operator/(const Matrix& a, Type b)
		{
			int	i, j;
			Matrix	ret(a.RowNum(), a.ColNum());

			for (i = 0; i < a.RowNum(); i++)
				for (j = 0; j < a.ColNum(); j++)
					ret[i][j] = a._data[i][j] / b;

			return ret;
		}
		friend Vector<Type> operator*(const Matrix<Type>& a, const Vector<Type>& b)
		{
			if (a.ColNum() != b.size())
				throw MatrixDimensionError("operator*(Mat a, Vec b) - a.colNum must be equal to vector size", a._rows, a._cols, (int)b.size(), -1);

			Vector<Type>	ret(a.RowNum());
			for (int i = 0; i < a.RowNum(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < a.ColNum(); j++)
					ret[i] += a._data[i][j] * b[j];
			}

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

		///////////////////////            Equality operations             //////////////////////
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

		bool IsEqual(const Matrix<Type>& b, Type eps = Defaults::MatrixEqualityPrecision) const
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
			return a.IsEqual(b, eps);
		}

		///////////////////            Trace, Inverse & Transpose             ///////////////////
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
				throw MatrixDimensionError("Matrix::Transpose - in-place Transpose possible only for square  matrix", _rows, _cols, -1, -1);

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

		///////////////////////////               I/O                 ///////////////////////////
		std::string to_string(int width, int precision) const
		{
			std::stringstream str;

			Print(str, width, precision);

			return str.str();
		}
		void   Print(std::ostream& stream, int width, int precision) const
		{
			stream << "Rows: " << RowNum() << " Cols: " << ColNum() << std::endl;

			for (size_t i = 0; i < RowNum(); i++)
			{
				stream << "[ ";
				for (size_t j = 0; j < ColNum(); j++)
				{
					stream << std::setw(width) << std::setprecision(precision) << _data[i][j] << ", ";
				}
				stream << " ]" << std::endl;
			}
		}
		void   Print(std::ostream& stream, int width, int precision, double epsZero) const
		{
			stream << "Rows: " << RowNum() << " Cols: " << ColNum() << std::endl;

			for (size_t i = 0; i < RowNum(); i++)
			{
				stream << "[ ";
				for (size_t j = 0; j < ColNum(); j++)
				{
					if (Abs(_data[i][j]) < epsZero)
						stream << std::setw(width) << std::setprecision(precision) << 0.0 << ", ";
					else
						stream << std::setw(width) << std::setprecision(precision) << _data[i][j] << ", ";
				}
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
				for (size_t i = 0; i < outMat.RowNum(); i++)
					for (size_t j = 0; j < outMat.ColNum(); j++)
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
				for (size_t i = 0; i < mat.RowNum(); i++)
				{
					for (size_t j = 0; j < mat.ColNum(); j++)
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

	////////////////////               Default Matrix typdefs                 ////////////////////
	typedef Matrix<int>     MatrixInt;
	typedef Matrix<float>   MatrixFlt;
	typedef Matrix<Real>    MatrixDbl;
	typedef Matrix<Complex> MatrixComplex;

	typedef Matrix<int>     MatI;
	typedef Matrix<float>   MatF;
	typedef Matrix<Real>    MatD;
	typedef Matrix<Complex> MatC;
}
///////////////////////////   ./include/base/MatrixBandDiag.h   ///////////////////////////


namespace MML
{
	template<class Type>
	class TridiagonalMatrix
	{
	private:
		int _dim;
		Vector<Type> _a;
		Vector<Type> _diag;
		Vector<Type> _c;

	public:
		TridiagonalMatrix(int dim, const Vector<Type>& a, const Vector<Type>& diag, const Vector<Type>& c) : _a(a), _diag(diag), _c(c), _dim(dim)
		{
			if (a.size() != dim || diag.size() != dim || c.size() != dim)
				throw("Error in TridiagonalMatrix constructor: wrong dimensions");
		}

		TridiagonalMatrix(int dim, std::initializer_list<Type> values) : _dim(dim), _a(dim), _diag(dim), _c(dim)
		{
			if (values.size() != dim * 3 - 2)
				throw("Error in TridiagonalMatrix constructor: wrong dimensions");

			auto val = values.begin();
			_a[0] = 0.0;
			_diag[0] = *val++;
			_c[0] = *val++;
			for (int i = 1; i < dim - 1; ++i)
			{
				_a[i] = *val++;
				_diag[i] = *val++;
				_c[i] = *val++;
			}
			_a[dim - 1] = *val++;
			_diag[dim - 1] = *val++;
			_c[dim - 1] = 0.0;
		}

		int RowNum() const { return _dim; }
		int ColNum() const { return _dim; }

		Type  operator()(int i, int j) const {
			if (i == j)
				return _diag[i];
			else if (i == j - 1)
				return _c[i];
			else if (i == j + 1 && j < _dim - 1)
				return _a[i];
			else
				return 0.0;
		}
		Type& operator()(int i, int j) {
			if (i == j)
				return _diag[i];
			else if (i == j - 1)
				return _c[i];
			else if (i == j + 1 && j < _dim - 1)
				return _a[i];
			else
				throw MatrixAccessBoundsError("TridiagonalMatrix::operator()", i, j, _dim, _dim);
		}

		// TODO 0.9 - dodati IsEqual
		// TODO 1.0 - imaju li smisla operacije s regularnim matricama? I BandDiag?
		TridiagonalMatrix operator+(const TridiagonalMatrix& b) const
		{
			if (_dim != b._dim)
				throw MatrixDimensionError("TridiagonalMatrix::operator+() - must be same dim", _dim, _dim, b._dim, b._dim);

			TridiagonalMatrix temp(*this);
			for (int i = 0; i < _a.size(); i++)
				temp._a[i] += b._a[i];
			for (int i = 0; i < _diag.size(); i++)
				temp._diag[i] += b._diag[i];
			for (int i = 0; i < _c.size(); i++)
				temp._c[i] += b._c[i];

			return temp;
		}
		TridiagonalMatrix operator-(const TridiagonalMatrix& b) const
		{
			if (_dim != b._dim)
				throw MatrixDimensionError("TridiagonalMatrix::operator+() - must be same dim", _dim, _dim, b._dim, b._dim);

			TridiagonalMatrix temp(*this);
			for (int i = 0; i < _a.size(); i++)
				temp._a[i] -= b._a[i];
			for (int i = 0; i < _diag.size(); i++)
				temp._diag[i] -= b._diag[i];
			for (int i = 0; i < _c.size(); i++)
				temp._c[i] -= b._c[i];

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
			int j, n = _a.size();
			Real bet;
			Vector<Type> gam(n);

			if (_diag[0] == 0.0)
				throw("Error 1 in tridag");

			sol[0] = rhs[0] / (bet = _diag[0]);

			for (j = 1; j < n; j++) {
				gam[j] = _c[j - 1] / bet;
				bet = _diag[j] - _a[j] * gam[j];

				if (bet == 0.0)
					throw("Error 2 in tridag");

				sol[j] = (rhs[j] - _a[j] * sol[j - 1]) / bet;
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

	class BandDiagonalMatrix
	{
		// The array a[0..n-1][0..m1+m2] stores A as follows: The diagonal elements are in a[0..n-1][m1].
		// Subdiagonal elements are in a[j..n-1][0..m1-1] with j > 0 appropriate to the number of
		// elements on each subdiagonal. Superdiagonal elements are in a[0..j][m1+1..m1+m2] with
		// j < n-1 appropriate to the number of elements on each superdiagonal.        
		int _dim;
		int _m1, _m2;
		Matrix<Real> _data;
	public:
		int RowNum() const { return _dim; }
		int ColNum() const { return _dim; }

		BandDiagonalMatrix(int dim, int m1, int m2, const Matrix<Real>& data) : _dim(dim), _m1(m1), _m2(m2), _data(data)
		{
			if (data.RowNum() != dim || data.ColNum() != m1 + m2 + 1)
				throw("Error in BandDiagonalMatrix constructor: wrong dimensions");
		}

		Real  operator()(int i, int j) const {
			if (i > j + _m1 || j > i + _m2)
				return 0.0;
			else
				return _data[i][j - i + _m1];
		}
		Real& operator()(int i, int j) {
			if (i > j + _m1 || j > i + _m2)
				throw MatrixAccessBoundsError("BandDiagonalMatrix::operator()", i, j, _dim, _dim);
			else
				return _data[i][j - i + _m1];
		}

		// Matrix multiply b = A * x, where A is band-diagonal with m1 rows below the diagonal and
		// m2 rows above. The input vector is x[0..n-1] and the output vector is b[0..n-1]. 
		void banmul(Matrix<Real>& a, const int m1, const int m2, Vector<Real>& x, Vector<Real>& b)
		{
			int i, j, k, tmploop, n = a.RowNum();
			for (i = 0; i < n; i++) {
				k = i - m1;
				tmploop = std::min(m1 + m2 + 1, int(n - k));
				b[i] = 0.0;
				for (j = std::max(0, -k); j < tmploop; j++) b[i] += a[i][j] * x[j + k];
			}
		}

		// TODO - HIGH, SREDNJE, implement basic operations between BandDiag and BandDiag


		Vector<Real> operator*(Vector<Real>& x)
		{
			Vector<Real> b(x.size());
			banmul(_data, _m1, _m2, x, b);
			return b;
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
		static MatrixNM<Type, 1, M>  RowMatrixFromVector(const VectorN<Type, M>& b)
		{
			MatrixNM<Type, 1, M>  ret;
			for (int j = 0; j < M; j++)
				ret._vals[0][j] = b[j];

			return ret;
		}
		static MatrixNM<Type, N, 1>  ColumnMatrixFromVector(const VectorN<Type, N>& b)
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

		friend MatrixNM operator*(const MatrixNM& a, Type b)
		{
			int	i, j;
			MatrixNM	ret(a);

			for (i = 0; i < a.RowNum(); i++)
				for (j = 0; j < a.ColNum(); j++)
					ret._vals[i][j] *= b;

			return ret;
		}
		friend MatrixNM operator/(const MatrixNM& a, Type b)
		{
			int	i, j;
			MatrixNM	ret(a);

			for (i = 0; i < a.RowNum(); i++)
				for (j = 0; j < a.ColNum(); j++)
					ret._vals[i][j] /= b;

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

		friend VectorN<Type, N> operator*(const MatrixNM<Type, N, M>& a, const VectorN<Type, M>& b)
		{
			int	i, j;
			VectorN<Type, N>	ret;

			for (i = 0; i < a.RowNum(); i++)
			{
				ret[i] = 0;
				for (j = 0; j < a.ColNum(); j++)
					ret[i] += a._vals[i][j] * b[j];
			}

			return ret;
		}
		friend VectorN<Type, M> operator*(const VectorN<Type, N>& a, const MatrixNM<Type, N, M>& b)
		{
			int	i, j;
			VectorN<Type, M>	ret;

			for (i = 0; i < b.ColNum(); i++)
			{
				ret[i] = 0;
				for (j = 0; j < b.RowNum(); j++)
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

			for (size_t i = 0; i < RowNum(); i++)
			{
				stream << "[ ";
				for (size_t j = 0; j < ColNum(); j++)
				{
					stream << std::setw(width) << std::setprecision(precision) << _vals[i][j] << ", ";
				}
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

///////////////////////////   ./include/base/MatrixSparse.h   ///////////////////////////

namespace MML
{
    // TODO - BIG!!! implement sparse matrix
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
				_ptrData[0] = numElem > 0 ? new Type[numElem] : nullptr;

			for (int i = 1; i < dim; i++)
				_ptrData[i] = _ptrData[i - 1] + i;
		}

	public:
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

		typedef Type value_type;      // make T available externally

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
///////////////////////////   ./include/base/Matrix3D.h   ///////////////////////////

namespace MML
{
	template <class Type>
	class Matrix3D {
	private:
		int _n;
		int _m;
		int _k;
		Type*** _v;

	public:
		Matrix3D() : _n(0), _m(0), _k(0), _v(nullptr) {}

		Matrix3D(int n, int m, int k) : _n(n), _m(m), _k(k), _v(new Type** [n])
		{
			int i, j;

			_v[0] = new Type * [n * m];
			_v[0][0] = new Type[n * m * k];

			for (j = 1; j < m; j++)
				_v[0][j] = _v[0][j - 1] + k;

			for (i = 1; i < n; i++) {
				_v[i] = _v[i - 1] + m;
				_v[i][0] = _v[i - 1][0] + m * k;

				for (j = 1; j < m; j++)
					_v[i][j] = _v[i][j - 1] + k;
			}
		}

		~Matrix3D()
		{
			if (_v != NULL) {
				delete[](_v[0][0]);
				delete[](_v[0]);
				delete[](_v);
			}
		}

		//subscripting: pointer to row i
		inline Type** operator[](const int i) { return _v[i]; }
		inline const Type* const* operator[](const int i) const { return _v[i]; }

		Type  operator()(int i, int j, int k) const { return _v[i][j][k]; }
		Type& operator()(int i, int j, int k) { return _v[i][j][k]; }

		inline int dim1() const { return _n; }
		inline int dim2() const { return _m; }
		inline int dim3() const { return _k; }
	};
}

///////////////////////////   ./include/base/Tensor.h   ///////////////////////////



namespace MML
{
	enum TensorIndexType { CONTRAVARIANT, COVARIANT };

	template <int N>
	class Tensor2 : public ITensor2<N>
	{
		// TODO - MED, LAKO, u MatrixNM!
		// TODO - MED, LAKO, osnovne operacije +, -, *
		// TODO - contraction
		Real _coeff[N][N] = { 0 };
	public:
		int _numContravar;
		int _numCovar;
		bool _isContravar[2];

		Tensor2(int nContra, int nCo) : _numContravar(nContra), _numCovar(nCo)
		{
			// mora biti covar + contra == 2
			if (_numContravar + _numCovar != 2)
				throw TensorCovarContravarNumError("Tensor2 ctor, wrong number of contravariant and covariant indices", nContra, nCo);

			for (int i = 0; i < nContra; i++)
				_isContravar[i] = true;

			for (int i = nContra; i < nContra + nCo; i++)
				_isContravar[i] = false;
		}
		Tensor2(TensorIndexType first, TensorIndexType second)
		{
			if (first == CONTRAVARIANT)
			{
				_numContravar++; _isContravar[0] = true;
			}
			else
			{
				_numCovar++;     _isContravar[0] = false;
			}

			if (second == CONTRAVARIANT)
			{
				_numContravar++; _isContravar[1] = true;
			}
			else
			{
				_numCovar++;     _isContravar[1] = false;
			}
		}

		int   NumContravar() const { return _numContravar; }
		int   NumCovar()     const { return _numCovar; }

		Real  Component(int i, int j) const { return _coeff[i][j]; }
		Real& Component(int i, int j) { return _coeff[i][j]; }

		void   Print(std::ostream& stream, int width, int precision) const
		{
			stream << "N = " << N << std::endl;

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
			stream << "N = " << N << std::endl;

			for (size_t i = 0; i < N; i++)
			{
				stream << "[ ";
				for (size_t j = 0; j < N; j++)
					stream << a._coeff[i][j] << ", ";
				stream << " ]" << std::endl;
			}
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

		Tensor3(int nContra, int nCo) : _numContravar(nContra), _numCovar(nCo) {
			// mora biti covar + contra == 3
			if (_numContravar + _numCovar != 3)
				throw TensorCovarContravarNumError("Tensor3 ctor, wrong number of contravariant and covariant indices", nContra, nCo);

			for (int i = 0; i < nContra; i++)
				_isContravar[i] = true;

			for (int i = nContra; i < nContra + nCo; i++)
				_isContravar[i] = false;
		}

		int   NumContravar() const { return _numContravar; }
		int   NumCovar()     const { return _numCovar; }

		Real  Component(int i, int j, int k) const { return _coeff[i][j][k]; }
		Real& Component(int i, int j, int k) { return _coeff[i][j][k]; }
	};

	template <int N>
	class Tensor4 : public ITensor4<N>
	{
		Real _coeff[N][N][N][N] = { 0 };
	public:
		int _numContravar;
		int _numCovar;
		bool _isContravar[4];

		Tensor4(int nContra, int nCo) : _numContravar(nContra), _numCovar(nCo) {
			if (_numContravar + _numCovar != 4)
				throw TensorCovarContravarNumError("Tensor4 ctor, wrong number of contravariant and covariant indices", nContra, nCo);

			for (int i = 0; i < nContra; i++)
				_isContravar[i] = true;

			for (int i = nContra; i < nContra + nCo; i++)
				_isContravar[i] = false;
		}

		int   NumContravar() const { return _numContravar; }
		int   NumCovar()     const { return _numCovar; }

		Real  Component(int i, int j, int k, int l) const { return _coeff[i][j][k][l]; }
		Real& Component(int i, int j, int k, int l) { return _coeff[i][j][k][l]; }
	};

	template <int N>
	class Tensor5 : public ITensor5<N>
	{
		Real _coeff[N][N][N][N][N];
	public:
		int _numContravar;
		int _numCovar;
		bool _isContravar[5];

		Tensor5(int nContra, int nCo) : _numContravar(nContra), _numCovar(nCo) {
			if (_numContravar + _numCovar != 5)
				throw TensorCovarContravarNumError("Tensor5 ctor, wrong number of contravariant and covariant indices", nContra, nCo);

			for (int i = 0; i < nContra; i++)
				_isContravar[i] = true;

			for (int i = nContra; i < nContra + nCo; i++)
				_isContravar[i] = false;

		}

		int   NumContravar() const { return _numContravar; }
		int   NumCovar()     const { return _numCovar; }

		Real  Component(int i, int j, int k, int l, int m) const { return _coeff[i][j][k][l][m]; }
		Real& Component(int i, int j, int k, int l, int m) { return _coeff[i][j][k][l][m]; }
	};
}
///////////////////////////   ./include/base/Polynom.h   ///////////////////////////


namespace MML
{
	// TODO - HIGH, TESKO?, dodati RealPolynom klasu, izvedenu iz IRealFunction implementirati derive (derivs - sve, first, sec, third zasebno)& integrate
	template <typename _Field, typename _CoefType = Real>
	class Polynom
	{
	private:
		std::vector<_CoefType> _vecCoef;
	public:
		Polynom() {}
		Polynom(int n) { _vecCoef.resize(n + 1); }
		Polynom(const std::vector<_CoefType>& vecCoef) : _vecCoef(vecCoef) {}
		Polynom(std::initializer_list<_CoefType> list) : _vecCoef(list) {}
		Polynom(const Polynom& Copy) : _vecCoef(Copy._vecCoef) {}
		~Polynom() {}

		Polynom& operator=(const Polynom& Copy) { _vecCoef = Copy._vecCoef; return *this; }

		int  GetDegree() const { return (int)_vecCoef.size() - 1; }
		void SetDegree(int newDeg) { _vecCoef.resize(newDeg + 1); }

		_Field  operator[] (int i) const { return _vecCoef[i]; }
		_Field& operator[] (int i) { return _vecCoef[i]; }

		_Field operator() (const _Field& x) {
			int j = GetDegree();
			_Field p = _vecCoef[j];

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

		bool operator==(const Polynom& b) const
		{
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

			while (nv >= 0 && v._vecCoef[nv] == 0.) nv--;

			if (nv < 0)
				throw("poldiv divide by zero polynomial");

			Polynom r = u;
			Polynom q(u.GetDegree());
			for (k = n - nv; k >= 0; k--)
			{
				q[k] = r[nv + k] / v[nv];
				for (j = nv + k - 1; j >= k; j--)
					r[j] -= q[k] * v[j - k];
			}
			for (j = nv; j <= n; j++) r[j] = 0.0;

			int nq = q.GetDegree();
			while (nq >= 0 && q[nq] == 0.)
				nq--;

			qout.SetDegree(nq - 1);
			for (j = 0; j < nq; j++)
				qout[j] = q[j];

			rout.SetDegree(nv - 1);
			for (j = 0; j < nv; j++)
				rout[j] = r[j];
		}

		friend Polynom operator*(const Polynom& a, _CoefType b)
		{
			Polynom ret;
			ret._vecCoef.resize(a.GetDegree() + 1);
			for (int i = 0; i <= a.GetDegree(); i++)
				ret._vecCoef[i] = a._vecCoef[i] * b;
			return ret;
		}

		friend Polynom operator*(_CoefType a, const Polynom& b)
		{
			Polynom ret;
			ret._vecCoef.resize(b.GetDegree() + 1);
			for (int i = 0; i <= b.GetDegree(); i++)
				ret._vecCoef[i] = a * b._vecCoef[i];
			return ret;
		}

		friend Polynom operator/(const Polynom& a, _CoefType b)
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
			for (int i = (int)_vecCoef.size() - 1; i >= 0; i--)
			{
				if (std::abs(_vecCoef[i]) == 0.0)
					continue;

				if (i != _vecCoef.size() - 1)
					stream << " + ";

				if (i == 0)
					stream << std::setw(width) << std::setprecision(precision) << _vecCoef[i];
				else if (i == 1)
					stream << std::setw(width) << std::setprecision(precision) << _vecCoef[i] << " * x" << i;
				else
					stream << std::setw(width) << std::setprecision(precision) << _vecCoef[i] << " * x^" << i;
			}

			return stream;
		}

		std::ostream& Print(std::ostream& stream) const
		{
			for (int i = (int)_vecCoef.size() - 1; i >= 0; i--)
			{
				if (std::abs(_vecCoef[i]) == 0.0)
					continue;

				if (i != _vecCoef.size() - 1)
					stream << " + ";

				if (i == 0)
					stream << _vecCoef[i];
				else if (i == 1)
					stream << _vecCoef[i] << " * x" << i;
				else
					stream << _vecCoef[i] << " * x^" << i;
			}

			return stream;
		}

		friend std::ostream& operator<<(std::ostream& stream, Polynom& a)
		{
			a.Print(stream);

			return stream;
		}
	};

	typedef Polynom<Real, Real>         RealPolynom;
	typedef Polynom<Complex, Complex>   ComplexPolynom;

	typedef Polynom<MatrixNM<Real, 2, 2>, Real>       Matrix2Polynom;
	typedef Polynom<MatrixNM<Real, 3, 3>, Real>       Matrix3Polynom;
	typedef Polynom<MatrixNM<Real, 4, 4>, Real>       Matrix4Polynom;
}

///////////////////////////   ./include/base/Algebra.h   ///////////////////////////



namespace MML
{
    // TODO - implement example groups - Z6, permutations
    // TODO - vector space + Gram Schmidt
    // TODO - linear transformations & OPERATORS
    
    // groups
    // finitve ima i numElem
    class GroupZ : public IGroup<int>
    {
        public:
        virtual int  identity() const { return 0;}
        virtual int  inverse(const int& a) const { return -a;}
        virtual int  op(const int& a, const int& b) const { return a + b;}
    };
    // FieldR
    // FieldC
}
///////////////////////////   ./include/base/Geometry.h   ///////////////////////////


namespace MML
{
	// TODO - MED, dodati Dist svugdje
	class Point2Cartesian
	{
	private:
		Real _x, _y;

	public:
		Real  X() const { return _x; }
		Real& X() { return _x; }
		Real  Y() const { return _y; }
		Real& Y() { return _y; }

		Point2Cartesian() {}
		Point2Cartesian(Real x, Real y) : _x(x), _y(y) {}

		Real Dist(const Point2Cartesian& b) const { return sqrt(POW2(b._x - _x) + POW2(b._y - _y)); }

		Point2Cartesian operator+(const Point2Cartesian& b) const { return Point2Cartesian(_x + b._x, _y + b._y); }
		Point2Cartesian operator-(const Point2Cartesian& b) const { return Point2Cartesian(_x - b._x, _y - b._y); }

		friend Point2Cartesian operator*(const Point2Cartesian& a, Real b) { return Point2Cartesian(a._x * b, a._y * b); }
		friend Point2Cartesian operator*(Real a, const Point2Cartesian& b) { return Point2Cartesian(a * b._x, a * b._y); }
		friend Point2Cartesian operator/(const Point2Cartesian& a, Real b) { return Point2Cartesian(a._x / b, a._y / b); }

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

		Point2Polar() {}
		Point2Polar(Real r, Real phi) : _r(r), _phi(phi) {}

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

		Point3Cartesian() {}
		Point3Cartesian(Real x, Real y, Real z) : _x(x), _y(y), _z(z) {}

		Real Dist(const Point3Cartesian& b) const { return sqrt(POW2(b._x - _x) + POW2(b._y - _y) + POW2(b._z - _z)); }

		Point3Cartesian operator+(const Point3Cartesian& b) const { return Point3Cartesian(_x + b._x, _y + b._y, _z + b._z); }
		Point3Cartesian operator-(const Point3Cartesian& b) const { return Point3Cartesian(_x - b._x, _y - b._y, _z - b._z); }

		friend Point3Cartesian operator*(const Point3Cartesian& a, Real b) { return Point3Cartesian(a._x * b, a._y * b, a._z * b); }
		friend Point3Cartesian operator*(Real a, const Point3Cartesian& b) { return Point3Cartesian(a * b._x, a * b._y, a * b._z); }
		friend Point3Cartesian operator/(const Point3Cartesian& a, Real b) { return Point3Cartesian(a._x / b, a._y / b, a._z / b); }
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

		Triangle() {}
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
		Vector3Cartesian(const VectorN<Real, 3>& b) : VectorN<Real, 3>{ b[0], b[1], b[2] } {}
		Vector3Cartesian(Real x, Real y, Real z) : VectorN<Real, 3>{ x, y, z } {}
		Vector3Cartesian(std::initializer_list<Real> list) : VectorN<Real, 3>(list) { }
		Vector3Cartesian(const Point3Cartesian& a, const Point3Cartesian& b)
		{
			_val[0] = b.X() - a.X();
			_val[1] = b.Y() - a.Y();
			_val[2] = b.Z() - a.Z();
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

		friend Real VectorsAngle(const Vector3Cartesian& a, const Vector3Cartesian& b)
		{
			Real cos_phi = ScalarProd(a, b) / (a.NormL2() * b.NormL2());

			return acos(cos_phi);
		}
	};

	class Vector3Spherical : public VectorN<Real, 3>
	{
	public:
		Real    R()     const { return _val[0]; }
		Real& R() { return _val[0]; }
		Real    Theta() const { return _val[1]; }
		Real& Theta() { return _val[1]; }
		Real    Phi()   const { return _val[2]; }
		Real& Phi() { return _val[2]; }

		Vector3Spherical() : VectorN<Real, 3>{ 0.0, 0.0, 0.0 } {}
		Vector3Spherical(const VectorN<Real, 3>& b) : VectorN<Real, 3>{ b[0], b[1], b[2] } {}
		Vector3Spherical(Real r, Real theta, Real phi) : VectorN<Real, 3>{ r, theta, phi } {}
		Vector3Spherical(std::initializer_list<Real> list) : VectorN<Real, 3>(list) { }

		// TODO - HIGH, HARD, osnovne operacije +, -
		Vector3Spherical GetAsUnitVector() const
		{
			return Vector3Spherical{ 1.0, Theta() , Phi() };
		}
		Vector3Spherical GetAsUnitVectorAtPos(const Vector3Spherical& pos) const
		{
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
		Vector3Cylindrical GetAsUnitVector() const
		{
			return Vector3Cylindrical{ (*this) / NormL2() };
		}
		Vector3Cylindrical GetAsUnitVectorAtPos(const Vector3Cylindrical& pos) const
		{
			return Vector3Cylindrical{ R(), Phi() / pos.R(), Z() };
		}
	};

	typedef Vector3Cartesian    Vec3Cart;
	typedef Vector3Spherical    Vec3Sph;
	typedef Vector3Cylindrical  Vec3Cyl;
}
///////////////////////////   ./include/base/StdFunctions.h   ///////////////////////////


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

///////////////////////////   ./include/interfaces/IFunction.h   ///////////////////////////


namespace MML
{
	//////////////////////////////////////////////////////////////////////
	template<typename _RetType, typename _ArgType>
	class IFunction
	{
	public:
		virtual _RetType operator()(_ArgType) const = 0;
	};

	//////////////////////////////////////////////////////////////////////
	class IRealFunction : public IFunction<Real, Real>
	{
	public:
		virtual Real operator()(Real) const = 0;

		void GetValues(Real x1, Real x2, int numPnt, Vector<Real>& outX, Vector<Real>& outY)
		{
			outX.Resize(numPnt);
			outY.Resize(numPnt);

			for (int i = 0; i < numPnt; i++) {
				outX[i] = x1 + i * (x2 - x1) / (numPnt - 1);
				outY[i] = (*this)(outX[i]);
			}
		}
		void Print(Real x1, Real x2, int numPnt)
		{
			for (int i = 0; i < numPnt; i++) {
				Real x = x1 + i * (x2 - x1) / (numPnt - 1);
				std::cout << x << " " << (*this)(x) << std::endl;
			}
		}
	};

	//////////////////////////////////////////////////////////////////////
	template<int N>
	class IScalarFunction : public IFunction<Real, const VectorN<Real, N>&>
	{
	public:
		virtual Real operator()(const VectorN<Real, N>& x) const = 0;
	};

	//////////////////////////////////////////////////////////////////////
	template<int N>
	class IRealToVectorFunction : public IFunction<VectorN<Real, N>, Real>
	{
	public:
		virtual VectorN<Real, N> operator()(Real x) const = 0;
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
	};

	//////////////////////////////////////////////////////////////////////
	template<int N>
	class IParametricCurve : public IRealToVectorFunction<N>
	{
	public:
		virtual VectorN<Real, N> operator()(Real x) const = 0;

		virtual Real getMinT() const = 0;
		virtual Real getMaxT() const = 0;

		// TODO - finish IsClosed() 
		bool isClosed(double tCycleStart, double tCycleEnd, double tVerificationIntervalEnd) const {
			// u t1 i t2 moraju biti isti
			// i onda do kraja intervala provjeriti da li se po tockama trace-a slaze
			return false;
		}

		std::vector<VectorN<Real, N>> GetTrace(double t1, double t2, int numPoints) const
		{
			std::vector<VectorN<Real, N>> ret;
			double deltaT = (t2 - t1) / (numPoints - 1);
			for (Real t = t1; t <= t2; t += deltaT)
				ret.push_back((*this)(t));
			return ret;
		}
	};

	//////////////////////////////////////////////////////////////////////
	// simple regular surface, defined on rectangular coordinate patch
	template<int N>
	class IParametricSurface : public IFunction<VectorN<Real, N>, const VectorN<Real, 2>&>
	{
	public:
		virtual VectorN<Real, N> operator()(Real u, Real w) const = 0;

		virtual Real getMinX() const = 0;
		virtual Real getMaxX() const = 0;
		virtual Real getMinY() const = 0;
		virtual Real getMaxY() const = 0;

		virtual VectorN<Real, N> operator()(const VectorN<Real, 2>& coord) const
		{
			return operator()(coord[0], coord[1]);
		}
		bool Serialize2DCartesian(Real x1, Real x2, int numPointsX, Real y1, Real y2, int numPointsY, std::string fileName) const
		{
			return false;
		}
	};

	// complex surface, with fixed u limits, but variable w limits (dependent on u)
	template<int N>
	class IParametricSurfaceComplex : public IFunction<VectorN<Real, N>, const VectorN<Real, 2>&>
	{
	public:
		virtual VectorN<Real, N> operator()(Real u, Real w) const = 0;

		virtual Real getMinX() const = 0;
		virtual Real getMaxX() const = 0;
		virtual Real getMinY(Real x) const = 0;
		virtual Real getMaxY(Real x) const = 0;

		virtual VectorN<Real, N> operator()(const VectorN<Real, 2>& coord) const
		{
			return operator()(coord[0], coord[1]);
		}
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
	};

	template<typename VectorFrom, typename VectorTo, int N>
	class ICoordTransfWithInverse : public virtual ICoordTransf<VectorFrom, VectorTo, N>
	{
	public:
		virtual       VectorFrom          transfInverse(const VectorTo& in) const = 0;
		virtual const IScalarFunction<N>& inverseCoordTransfFunc(int i) const = 0;
	};
}
///////////////////////////   ./include/interfaces/IODESystem.h   ///////////////////////////


namespace MML
{
	class IODESystem
	{
	public:
		virtual int     getDim() const = 0;
		virtual void    derivs(const Real, const Vector<Real>&, Vector<Real>&) const = 0;
		void operator()(const Real t, const Vector<Real>& x, Vector<Real>& dxdt) const { return derivs(t, x, dxdt); }
	};
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

		Point2Cartesian     StartPoint() const { return _point; }
		Point2Cartesian& StartPoint() { return _point; }

		Vector2Cartesian    Direction() const { return _direction; }
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

		Point2Cartesian     StartPoint() const { return _point1; }
		Point2Cartesian& StartPoint() { return _point1; }

		Point2Cartesian     EndPoint()  const { return _point2; }
		Point2Cartesian& EndPoint() { return _point2; }

		Point2Cartesian PointOnSegment(Real t)
		{
			if (t < 0.0 || t > 1.0)
				throw std::invalid_argument("t must be in [0,1]");

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

		// TODO - IsConvex()
		// TODO - Triangularization

		// TODO - IsInside
		bool IsInside(Point2Cartesian pnt) const
		{
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
		Line3D(const Point3Cartesian& pnt, const Vector3Cartesian dir)
		{
			_point = pnt;
			_direction = dir.GetAsUnitVector();
		}

		Line3D(const Point3Cartesian& a, const Point3Cartesian& b)
		{
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
			return ScalarProd(Direction(), b.Direction()) == 0.0f;
		}
		// TODO 0.9 - VERIFY distance Line - Point3
		Real Dist(const Point3Cartesian& pnt) const
		{
			Real dist = 0.0;

			const Real a = pnt.X();
			const Real b = pnt.Y();
			const Real c = pnt.Z();

			const Real x1 = StartPoint().X();
			const Real y1 = StartPoint().Y();
			const Real z1 = StartPoint().Z();

			const Real l = Direction().X();
			const Real m = Direction().Y();
			const Real n = Direction().Z();

			dist = POW2((a - x1) * m - (b - y1) * l) + POW2((b - y1) * n - (c - z1) * m) + POW2((c - z1) * l - (a - x1) * n);
			Real denom = l * l + m * m + n * n;

			return sqrt(dist / denom);
		}

		// Real DistChatGPT(const Point3Cartesian &pnt)
		// {
		//     Vector3Cartesian p0 = StartPoint();
		//     Vector3Cartesian p1 = p0 + Direction();
		//     Vector3Cartesian p = pnt;

		//     Vector3Cartesian numerator = (p - p0).CrossProduct(p - p1);
		//     Real denominator = (p1 - p0).Length();

		//     return numerator.Length() / denominator;
		// }

		Point3Cartesian NearestPoint(const Point3Cartesian& pnt) const
		{
			// https://math.stackexchange.com/questions/1521128/given-a-line-and-a-point-in-3d-how-to-find-the-closest-point-on-the-line         
			Point3Cartesian Q(StartPoint());
			Point3Cartesian R = PointOnLine(1.0);

			Vector3Cartesian RQ(Q, R);
			Vector3Cartesian QP(pnt, StartPoint());

			double t = QP.ScalarProductCartesian(RQ) / RQ.ScalarProductCartesian(RQ);

			return Q - t * RQ;
		}
		// TODO - pravac koji prolazi kroz tocku i sijece zadani pravac okomito
		// TODO - sjeciste dva pravca
	};

	class SegmentLine3D
	{
	private:
		Point3Cartesian _point1;
		Point3Cartesian _point2;

	public:
		SegmentLine3D(Point3Cartesian pnt1, Point3Cartesian pnt2) : _point1(pnt1), _point2(pnt2)
		{ }

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
			if (normal.NormL2() == 0.0)
				throw std::runtime_error("Plane3D ctor - normal is null vector");

			Vector3Cartesian unitNormal = normal.GetAsUnitVector();

			_A = unitNormal.X();
			_B = unitNormal.Y();
			_C = unitNormal.Z();
			_D = -(a.X() * unitNormal.X() + a.Y() * unitNormal.Y() + a.Z() * unitNormal.Z());
		}

		Plane3D(const Point3Cartesian& a, const Point3Cartesian& b, const Point3Cartesian& c)
			: Plane3D(a, VectorProd(Vector3Cartesian(a, b), Vector3Cartesian(a, c)))
		{ }

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
			Point3Cartesian x(seg_x, 0, 0);
			Point3Cartesian y(0, seg_y, 0);
			Point3Cartesian z(0, 0, seg_z);

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
			outseg_x = -_D / _A;
			outseg_y = -_D / _B;
			outseg_z = -_D / _C;
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
			// TODO - check if it is already ona the plane!
				// from given point and normal to plane form a Line3D
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

			// Calculate the distance between the line and the plane
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
			// TODO finishs
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
		{}

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
		RectSurface3D(Point3Cartesian pnt1, Point3Cartesian pnt2, Point3Cartesian pnt3, Point3Cartesian pnt4)
			: _pnt1(pnt1), _pnt2(pnt2), _pnt3(pnt3), _pnt4(pnt4)
		{
			// provjeriti da li su sve u ravnini!

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
		virtual Real getMinX() const { return _minX; }
		virtual Real getMaxX() const { return _maxX; }
		virtual Real getMinY() const { return _minY; }
		virtual Real getMaxY() const { return _maxY; }

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

	// TODO - IntegrableSolid? koji ima i potrebne funkcije kojima definira granice tijela?
	class IntegrableVolume3D
	{
		// osigurava da se znaju funkcije koje definiraju granice tijela
		// da se moze obaviti volume integracija
	};

	class SolidSurfaces3D
	{
	public:
		// solid body in 3D
		std::vector<RectSurface3D> _surfaces;

	public:
		bool IsInside(const Point3Cartesian& pnt) const
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

///////////////////////////   ./include/base/QuadraticForm.h   ///////////////////////////



namespace MML
{
	// TODO - MED, finish this
	template<int N>
	class QuadraticFormN : public IScalarFunction<N>
	{
		MatrixNM<Real, N, N> _mat;
	public:
		Real operator()(const VectorN<Real, N>& x) const
		{
			return x.ScalarProductCartesian(_mat * x);
		}

		// TODO - treat as polynom! Print

	};
}
///////////////////////////   ./include/base/VectorSpace.h   ///////////////////////////

//#include "interfaces/IFunction.h"


namespace MML
{
	// TODO - MED, vector space + Gram Schmidt
	// TODO - MED, linear transformations & OPERATORS
	template<int N>
	class RealVectorSpaceN : public HilbertSpace<Real, VectorN<Real, N>>
	{
	public:
		virtual Real  scal_prod(const VectorN<Real, N>& a, const VectorN<Real, N>& b) const
		{
			return a.ScalarProductCartesian(b);
		}
	};

	template<int N>
	class ComplexVectorSpaceN : public HilbertSpace<Complex, VectorN<Complex, N>>
	{
	public:
		virtual Real  scal_prod(const VectorN<Complex, N>& a, const VectorN<Complex, N>& b) const
		{
			Real product = 0.0;
			for (int i = 0; i < N; i++)
				product += (a[i] * std::conj(b[i])).real();
			return product;
		}
	};
}
///////////////////////////   ./include/base/LinearFunctional.h   ///////////////////////////



namespace MML
{
	template <int N, typename _Field = Real>
	class LinearFunctionalN
	{
	private:
		// TODO - LOW, this should be an array, not a vector
		VectorN<_Field, N> _vecCoef;
	public:
		LinearFunctionalN() {}
		LinearFunctionalN(const VectorN<_Field, N>& vecCoef) : _vecCoef(vecCoef) {}
		LinearFunctionalN(std::initializer_list<_Field> list) : _vecCoef(list) {}

		LinearFunctionalN(const LinearFunctionalN& Copy) : _vecCoef(Copy._vecCoef) {}
		~LinearFunctionalN() {}

		LinearFunctionalN& operator=(const LinearFunctionalN& Copy) { _vecCoef = Copy._vecCoef; return *this; }

		int Dim() const { return N; }

		_Field operator()(const VectorN<_Field, N>& vecX) const
		{
			_Field result = 0.0;
			for (int i = 0; i < N; i++)
				result += _vecCoef[i] * vecX[i];
			return result;
		}

		LinearFunctionalN operator+(const LinearFunctionalN& b) const
		{
			LinearFunctionalN result;
			for (int i = 0; i < N; i++)
				result._vecCoef[i] = _vecCoef[i] + b._vecCoef[i];
			return result;
		}

		LinearFunctionalN operator-(const LinearFunctionalN& b) const
		{
			LinearFunctionalN result;
			for (int i = 0; i < N; i++)
				result._vecCoef[i] = _vecCoef[i] - b._vecCoef[i];
			return result;
		}

		LinearFunctionalN operator*(_Field b) const
		{
			LinearFunctionalN result;
			for (int i = 0; i < _vecCoef.size(); i++)
				result._vecCoef[i] = _vecCoef[i] * b;
			return result;
		}
	};

	template <int N>
	class RealLinearFunctionalN : public IScalarFunction<N>
	{
	private:
		VectorN<Real, N> _vecCoef;
	public:
		RealLinearFunctionalN() {}
		RealLinearFunctionalN(const VectorN<Real, N>& vecCoef) : _vecCoef(vecCoef) {}
		RealLinearFunctionalN(std::initializer_list<Real> list) : _vecCoef(list) {}

		RealLinearFunctionalN(const RealLinearFunctionalN& Copy) = default;
		~RealLinearFunctionalN() {}

		RealLinearFunctionalN& operator=(const RealLinearFunctionalN& Copy) { _vecCoef = Copy._vecCoef; return *this; }

		int Dim() const { return N; }

		Real operator()(const VectorN<Real, N>& vecX) const
		{
			Real result = 0.0;
			for (int i = 0; i < N; i++)
				result += _vecCoef[i] * vecX[i];
			return result;
		}

		RealLinearFunctionalN operator+(const RealLinearFunctionalN& b) const
		{
			RealLinearFunctionalN result;
			for (int i = 0; i < N; i++)
				result._vecCoef[i] = _vecCoef[i] + b._vecCoef[i];
			return result;
		}

		RealLinearFunctionalN operator-(const RealLinearFunctionalN& b) const
		{
			RealLinearFunctionalN result;
			for (int i = 0; i < N; i++)
				result._vecCoef[i] = _vecCoef[i] - b._vecCoef[i];
			return result;
		}

		RealLinearFunctionalN operator*(Real b) const
		{
			RealLinearFunctionalN result;
			for (int i = 0; i < _vecCoef.size(); i++)
				result._vecCoef[i] = _vecCoef[i] * b;
			return result;
		}
	};

	class RealLinearFunctional
	{
	private:
		std::vector<Real> _vecCoef;
	public:
		RealLinearFunctional() {}
		RealLinearFunctional(int n) : _vecCoef(n, 0.0) {}
		RealLinearFunctional(const std::vector<Real>& vecCoef) : _vecCoef(vecCoef) {}
		RealLinearFunctional(std::initializer_list<Real> list) : _vecCoef(list) {}

		RealLinearFunctional(const RealLinearFunctional& Copy) : _vecCoef(Copy._vecCoef) {}
		~RealLinearFunctional() {}

		int Dim() const { return (int)_vecCoef.size(); }

		RealLinearFunctional& operator=(const RealLinearFunctional& Copy) { _vecCoef = Copy._vecCoef; return *this; }

		Real operator()(const std::vector<Real>& vecX) const
		{
			if (vecX.size() != Dim())
				throw std::runtime_error("RealLinearFunctional::operator() - incompatible vector size");

			Real result = 0.0;
			for (int i = 0; i < _vecCoef.size(); i++)
				result += _vecCoef[i] * vecX[i];
			return result;
		}

		RealLinearFunctional operator+(const RealLinearFunctional& b) const
		{
			if (b.Dim() != Dim())
				throw std::runtime_error("RealLinearFunctional::operator+() - incompatible vector size");

			RealLinearFunctional result;
			for (int i = 0; i < Dim(); i++)
				result._vecCoef[i] = _vecCoef[i] + b._vecCoef[i];
			return result;
		}

		RealLinearFunctional operator-(const RealLinearFunctional& b) const
		{
			if (b.Dim() != Dim())
				throw std::runtime_error("RealLinearFunctional::operator-() - incompatible vector size");

			RealLinearFunctional result;
			for (int i = 0; i < Dim(); i++)
				result._vecCoef[i] = _vecCoef[i] - b._vecCoef[i];
			return result;
		}

		RealLinearFunctional operator*(Real b) const
		{
			RealLinearFunctional result;
			for (int i = 0; i < Dim(); i++)
				result._vecCoef[i] = _vecCoef[i] + b;
			return result;
		}
	};
}

///////////////////////////   ./include/base/LinearOperator.h   ///////////////////////////



namespace MML
{
	// TODO - MED, TESKO, finish  this
	template<int N>
	class RealLinearOperatorN : public ILinearOperator<VectorN<Real, N>, VectorN<Real, N>>
	{
		MatrixNM<Real, N, N> m_matrix;
	public:
		VectorN<Real, N>  operator()(const VectorN<Real, N>& x) const
		{
			return m_matrix * x;
		}
	};

	template<int N, int M>
	class RealLinearOperatorNM : public ILinearOperator<VectorN<Real, N>, VectorN<Real, M>>
	{
		MatrixNM<Real, N, M> m_matrix;
	public:
		VectorN<Real, N>  operator()(const VectorN<Real, M>& x) const
		{
			return m_matrix * x;
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
		static bool AreEqual(const Vector<Complex>& a, const Vector<Complex>& b, double eps = Defaults::ComplexEqualityPrecision)
		{
			if (a.size() != b.size())
				return false;

			for (int i = 0; i < a.size(); i++)
				if (std::abs(a[i].real() - b[i].real()) > eps || std::abs(a[i].imag() - b[i].imag()) > eps)
					return false;

			return true;
		}

		///////////////////                     Vector helpers                     ///////////////////
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
				MatrixNM<Complex, 4, 4>{ 0,              0,              0,              Complex(0, -1),
																 0,              0,              Complex(0, 1),  0,
																 0,              Complex(0, 1),  0,              0,
																 Complex(0, -1), 0,              0,              0 },
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

		template<class Type> static void MatrixDecompose(const Matrix<Type>& orig, Matrix<Type>& outSym, Matrix<Type>& outAntiSym)
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

	///////////////////////   BAND-DIAGONAL SOLVER    ////////////////////////////
	class BandDiagLUSolver {
		int n, m1, m2;
		Matrix<Real> au, al;
		Vector<int> indx;
		Real d;

	public:
		// TODO 0.9 - HIGH, SREDNJE, parm mora biti BandDiagonalMatrix!
		// TODO - handling singularne matrice
		BandDiagLUSolver(Matrix<Real>& a, const int mm1, const int mm2)
			: n(a.RowNum()), au(a), m1(mm1), m2(mm2), al(n, m1), indx(n)
		{
			const Real TINY = 1.0e-40;
			int i, j, k, l, mm;
			Real dum;
			mm = m1 + m2 + 1;
			l = m1;
			for (i = 0; i < m1; i++) {
				for (j = m1 - i; j < mm; j++) au[i][j - l] = au[i][j];
				l--;
				for (j = mm - l - 1; j < mm; j++) au[i][j] = 0.0;
			}
			d = 1.0;
			l = m1;
			for (k = 0; k < n; k++) {
				dum = au[k][0];
				i = k;
				if (l < n) l++;
				for (j = k + 1; j < l; j++) {
					if (std::abs(au[j][0]) > std::abs(dum)) {
						dum = au[j][0];
						i = j;
					}
				}
				indx[k] = i + 1;
				if (dum == 0.0) au[k][0] = TINY;
				if (i != k) {
					d = -d;
					for (j = 0; j < mm; j++) std::swap(au[k][j], au[i][j]);
				}
				for (i = k + 1; i < l; i++) {
					dum = au[i][0] / au[k][0];
					al[k][i - k - 1] = dum;
					for (j = 1; j < mm; j++) au[i][j - 1] = au[i][j] - dum * au[k][j];
					au[i][mm - 1] = 0.0;
				}
			}
		}
		void Solve(const Vector<Real>& b, Vector<Real>& x)
		{
			int i, j, k, l, mm;
			Real dum;
			mm = m1 + m2 + 1;
			l = m1;
			for (k = 0; k < n; k++) x[k] = b[k];
			for (k = 0; k < n; k++) {
				j = indx[k] - 1;
				if (j != k) std::swap(x[k], x[j]);
				if (l < n) l++;
				for (j = k + 1; j < l; j++) x[j] -= al[k][j - k - 1] * x[k];
			}
			l = 1;
			for (i = n - 1; i >= 0; i--) {
				dum = x[i];
				for (k = 1; k < l; k++) dum -= au[i][k] * x[k + i];
				x[i] = dum / au[i][0];
				if (l < mm) l++;
			}
		}

		Vector<Real> Solve(const Vector<Real>& b)
		{
			Vector<Real> x(b.size());
			Solve(b, x);
			return x;
		}

		Real Det() const
		{
			Real dd = d;
			for (int i = 0; i < n; i++)
				dd *= au[i][0];
			return dd;
		}
	};

	///////////////////////      LU DECOMPOSITION       /////////////////////////////
	template<class Type>
	class LUDecompositionSolver
	{
	private:
		int n;
		const Matrix<Type>& refOrig;

		Matrix<Type> lu;
		std::vector<int> indx;
		Real d;

	public:
		// TODO - HIGH, HIGH, TESKO, napraviti da se može odraditi i inplace, a ne da se kao sad uvijek kreira kopija
		LUDecompositionSolver(const Matrix<Type>& inMatRef) : n(inMatRef.RowNum()), refOrig(inMatRef), lu(inMatRef), indx(n)
			// ne može        LUDecompositionSolver(const Matrix<Type>  &inMatRef) : n(inMatRef.RowNum()), refOrig(inMatRef), lu(Matrix<Type>(inMatRef.RowNum(), inMatRef.ColNum())), indx(n) 
			// lu mora postati pointer!
		{
			// Given a Matrix<Real> a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise
			// permutation of itself. a and n are input. a is output, arranged as in equation (NR 2.3.14);
			// indx[1..n] is an output Vector<Real> that records the row permutation effected by the partial
			// pivoting; d is output as ±1 depending on whether the number of row interchanges was even
			// or odd, respectively. This routine is used in combination with lubksb to solve linear equations
			// or invert a Matrix<Real>.
			const Real TINY = 1.0e-40;
			int i, imax, j, k;
			Real big, temp;
			Type temp2;
			Vector<Type> vv(n);
			d = 1.0;
			for (i = 0; i < n; i++) {
				big = 0.0;
				for (j = 0; j < n; j++)
					if ((temp = Abs(lu[i][j])) > big) big = temp;
				if (big == 0.0)
					throw SingularMatrixError("LUDecompositionSolver::ctor - Singular Matrix");

				vv[i] = 1.0 / big;
			}
			for (k = 0; k < n; k++) {
				big = 0.0;
				imax = k;
				for (i = k; i < n; i++) {
					temp = Abs(vv[i] * lu[i][k]);
					if (temp > big) {
						big = temp;
						imax = i;
					}
				}
				if (k != imax) {
					for (j = 0; j < n; j++) {
						temp2 = lu[imax][j];
						lu[imax][j] = lu[k][j];
						lu[k][j] = temp2;
					}
					d = -d;
					vv[imax] = vv[k];
				}
				indx[k] = imax;
				if (lu[k][k] == Real{ 0.0 }) lu[k][k] = TINY;
				for (i = k + 1; i < n; i++) {
					temp2 = lu[i][k] /= lu[k][k];
					for (j = k + 1; j < n; j++)
						lu[i][j] -= temp2 * lu[k][j];
				}
			}
		}

		bool Solve(const Vector<Type>& b, Vector<Type>& x)
		{
			// Solves the set of n linear equations A·X = B. Here a[1..n][1..n] is input, not as the Matrix<Real>
			// A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
			// as the permutation Vector<Real> returned by ludcmp. b[1..n] is input as the right-hand side Vector<Real>
			// B, and returns with the solution Vector<Real> X. a, n, and indx are not modified by this routine
			// and can be left in place for successive calls with different right-hand sides b. This routine takes
			// into account the possibility that b will begin with many zero elements, so it is efficient for use
			// in Matrix<Real> inversion
			int i, ii = 0, ip, j;
			Type sum;

			if (b.size() != n || x.size() != n)
				return false;

			for (i = 0; i < n; i++)
				x[i] = b[i];

			for (i = 0; i < n; i++) {
				ip = indx[i];
				sum = x[ip];
				x[ip] = x[i];
				if (ii != 0)
					for (j = ii - 1; j < i; j++) sum -= lu[i][j] * x[j];
				else if (sum != Real{ 0.0 })
					ii = i + 1;
				x[i] = sum;
			}

			for (i = n - 1; i >= 0; i--) {
				sum = x[i];
				for (j = i + 1; j < n; j++) sum -= lu[i][j] * x[j];
				x[i] = sum / lu[i][i];
			}
			return true;
		}

		Vector<Type> Solve(const Vector<Type>& b)
		{
			Vector<Type> x(b.size());
			if (Solve(b, x) == true)
				return x;
			else
				throw VectorDimensionError("LUDecompositionSolver::Solve - bad sizes", b.size(), n);
		}

		void Solve(Matrix<Type>& b, Matrix<Type>& x)
		{
			int i, j, m = b.ColNum();

			if (b.RowNum() != n || x.RowNum() != n || b.ColNum() != x.ColNum())
				throw("LUdcmp::solve bad sizes");

			Vector<Type> xx(n);

			for (j = 0; j < m; j++) {
				for (i = 0; i < n; i++)
					xx[i] = b[i][j];

				Solve(xx, xx);

				for (i = 0; i < n; i++)
					x[i][j] = xx[i];
			}
		}

		// Using the stored LU decomposition, return in ainv the matrix inverse 
		void inverse(Matrix<Type>& ainv)
		{
			int i, j;
			ainv.Resize(n, n);
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) ainv[i][j] = 0.;
				ainv[i][i] = 1.;
			}
			Solve(ainv, ainv);
		}

		Type det()
		{
			Real dd = d;
			for (int i = 0; i < n; i++)
				dd *= lu[i][i];
			return dd;
		}

		// Improves a solution Vector<Real> x[1..n] of the linear set of equations A · X = B. The Matrix<Real>
		// a[1..n][1..n], and the Vector<Real>s b[1..n] and x[1..n] are input, as is the dimension n.
		// Also input is alud[1..n][1..n], the LU decomposition of a as returned by ludcmp, and
		// the Vector<Real> indx[1..n] also returned by that routine. On output, only x[1..n] is modified,
		// to an improved set of values
		void mprove(const Vector<Real>& b, Vector<Real>& x)
		{
			int i, j;
			Vector<Real> r(n);

			for (i = 0; i < n; i++) {
				long double sdp = -b[i];
				for (j = 0; j < n; j++)
					sdp += (long double)refOrig[i][j] * (long double)x[j];
				r[i] = sdp;
			}

			Solve(r, r);

			for (i = 0; i < n; i++)
				x[i] -= r[i];
		}
	};

	///////////////////////   CHOLESKY DECOMPOSITION    /////////////////////////////
	class CholeskyDecompositionSolver
	{
	private:
		int n;
		Matrix<Real> el;

	public:
		CholeskyDecompositionSolver(Matrix<Real>& a) : n(a.RowNum()), el(a)
		{
			// Given a positive-definite symmetric Matrix<Real> a[1..n][1..n], this routine constructs its Cholesky
			// decomposition, A = L · LT . On input, only the upper triangle of a need be given; it is not
			// modified. The Cholesky factor L is returned in the lower triangle of a, except for its diagonal
			// elements which are returned in p[1..n]
			int i, j, k;
			Real sum;

			if (el.ColNum() != n)
				throw("need square Matrix<Real>");

			for (i = 0; i < n; i++) {
				for (j = i; j < n; j++) {
					for (sum = el[i][j], k = i - 1; k >= 0; k--) sum -= el[i][k] * el[j][k];
					if (i == j) {
						if (sum <= 0.0)
							throw("Cholesky failed");
						el[i][i] = sqrt(sum);
					}
					else el[j][i] = sum / el[i][i];
				}
			}
			for (i = 0; i < n; i++) for (j = 0; j < i; j++) el[j][i] = 0.;
		}
		void Solve(const Vector<Real>& b, Vector<Real>& x)
		{
			// Solves the set of n linear equations A · x = b, where a is a positive-definite symmetric Matrix<Real>.
			// a[1..n][1..n] and p[1..n] are input as the output of the routine choldc. Only the lower
			// triangle of a is accessed. b[1..n] is input as the right-hand side Vector<Real>. The solution Vector<Real> is
			// returned in x[1..n]. a, n, and p are not modified and can be left in place for successive calls
			// with different right-hand sides b. b is not modified unless you identify b and x in the calling
			// sequence, which is allowed.
			int i, k;
			Real sum;
			if (b.size() != n || x.size() != n) throw("bad lengths in Cholesky");
			for (i = 0; i < n; i++) {
				for (sum = b[i], k = i - 1; k >= 0; k--) sum -= el[i][k] * x[k];
				x[i] = sum / el[i][i];
			}
			for (i = n - 1; i >= 0; i--) {
				for (sum = x[i], k = i + 1; k < n; k++) sum -= el[k][i] * x[k];
				x[i] = sum / el[i][i];
			}
		}
		Vector<Real> Solve(const Vector<Real>& b)
		{
			Vector<Real> x(b.size());
			Solve(b, x);
			return x;
		}
		void inverse(Matrix<Real>& ainv) {
			int i, j, k;
			Real sum;
			ainv.Resize(n, n);
			for (i = 0; i < n; i++) for (j = 0; j <= i; j++) {
				sum = (i == j ? 1. : 0.);
				for (k = i - 1; k >= j; k--) sum -= el[i][k] * ainv[j][k];
				ainv[j][i] = sum / el[i][i];
			}
			for (i = n - 1; i >= 0; i--) for (j = 0; j <= i; j++) {
				sum = (i < j ? 0. : ainv[j][i]);
				for (k = i + 1; k < n; k++) sum -= el[k][i] * ainv[j][k];
				ainv[i][j] = ainv[j][i] = sum / el[i][i];
			}
		}
		Real logdet() {
			Real sum = 0.;
			for (int i = 0; i < n; i++) sum += log(el[i][i]);
			return 2. * sum;
		}

		void elmult(Vector<Real>& y, Vector<Real>& b) {
			int i, j;
			if (b.size() != n || y.size() != n) throw("bad lengths");
			for (i = 0; i < n; i++) {
				b[i] = 0.;
				for (j = 0; j <= i; j++) b[i] += el[i][j] * y[j];
			}
		}
		void elsolve(Vector<Real>& b, Vector<Real>& y) {
			int i, j;
			Real sum;
			if (b.size() != n || y.size() != n) throw("bad lengths");
			for (i = 0; i < n; i++) {
				for (sum = b[i], j = 0; j < i; j++) sum -= el[i][j] * y[j];
				y[i] = sum / el[i][i];
			}
		}
	};

	///////////////////////   QR DECOMPOSITION    /////////////////////////////
	class QRDecompositionSolver
	{
	public:         // TODO make private!
		int n;
		Matrix<Real> qt, r;
		bool sing;

	public:
		QRDecompositionSolver(const Matrix<Real>& a) : n(a.RowNum()), qt(n, n), r(a), sing(false)
		{
			// Constructs the QR decomposition of a[1..n][1..n]. The upper triangular Matrix<Real> R is returned in the upper triangle of a, 
			// except for the diagonal elements of R which are returned in d[1..n]. 
			int i, j, k;
			Vector<Real> c(n), d(n);
			Real scale, sigma, sum, tau;
			for (k = 0; k < n - 1; k++) {
				scale = 0.0;
				for (i = k; i < n; i++) scale = std::max(scale, std::abs(r[i][k]));
				if (scale == 0.0) {
					sing = true;
					c[k] = d[k] = 0.0;
				}
				else {
					for (i = k; i < n; i++) r[i][k] /= scale;
					for (sum = 0.0, i = k; i < n; i++) sum += POW2(r[i][k]);
					sigma = SIGN(sqrt(sum), r[k][k]);
					r[k][k] += sigma;
					c[k] = sigma * r[k][k];
					d[k] = -scale * sigma;
					for (j = k + 1; j < n; j++) {
						for (sum = 0.0, i = k; i < n; i++) sum += r[i][k] * r[i][j];
						tau = sum / c[k];
						for (i = k; i < n; i++) r[i][j] -= tau * r[i][k];
					}
				}
			}
			d[n - 1] = r[n - 1][n - 1];
			if (d[n - 1] == 0.0) sing = true;
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) qt[i][j] = 0.0;
				qt[i][i] = 1.0;
			}
			for (k = 0; k < n - 1; k++) {
				if (c[k] != 0.0) {
					for (j = 0; j < n; j++) {
						sum = 0.0;
						for (i = k; i < n; i++)
							sum += r[i][k] * qt[i][j];
						sum /= c[k];
						for (i = k; i < n; i++)
							qt[i][j] -= sum * r[i][k];
					}
				}
			}
			for (i = 0; i < n; i++) {
				r[i][i] = d[i];
				for (j = 0; j < i; j++) r[i][j] = 0.0;
			}
		}

		// Solves the set of n linear equations A · x = b. a[1..n][1..n], c[1..n], and d[1..n] are
		// input as the output of the routine qrdcmp and are not modified. b[1..n] is input as the
		// right-hand side Vector<Real>, and is overwritten with the solution Vector<Real> on output. 
		void Solve(const Vector<Real>& b, Vector<Real>& x)
		{
			qtmult(b, x);
			RSolve(x, x);
		}
		Vector<Real> Solve(const Vector<Real>& b)
		{
			Vector<Real> x(b.size());
			Solve(b, x);
			return x;
		}
		void RSolve(const Vector<Real>& b, Vector<Real>& x)
		{
			// Solves the set of n linear equations R · x = b, where R is an upper triangular Matrix<Real> stored in
			// a and d. a[1..n][1..n] and d[1..n] are input as the output of the routine qrdcmp and
			// are not modified. b[1..n] is input as the right-hand side Vector<Real>, and is overwritten with the
			// solution Vector<Real> on output            
			int i, j;
			Real sum;
			if (sing)
				throw SingularMatrixError("QRDecompositionSolver::rsolve - attempting solve in a singular QR");

			for (i = n - 1; i >= 0; i--) {
				sum = b[i];
				for (j = i + 1; j < n; j++)
					sum -= r[i][j] * x[j];
				x[i] = sum / r[i][i];
			}
		}

		void qtmult(const Vector<Real>& b, Vector<Real>& x) {
			int i, j;
			Real sum;
			for (i = 0; i < n; i++) {
				sum = 0.;
				for (j = 0; j < n; j++)
					sum += qt[i][j] * b[j];
				x[i] = sum;
			}
		}
		void update(Vector<Real>& u, Vector<Real>& v)
		{
			// Given the QR decomposition of some n × n Matrix<Real>, calculates the QR decomposition of the
			// Matrix<Real> Q·(R+ u x v). The quantities are dimensioned as r[1..n][1..n], qt[1..n][1..n],
			// u[1..n], and v[1..n]. Note that QT is input and returned in qt.            
			int i, k;
			Vector<Real> w(u);
			for (k = n - 1; k >= 0; k--)
				if (w[k] != 0.0) break;
			if (k < 0) k = 0;
			for (i = k - 1; i >= 0; i--) {
				rotate(i, w[i], -w[i + 1]);
				if (w[i] == 0.0)
					w[i] = std::abs(w[i + 1]);
				else if (std::abs(w[i]) > std::abs(w[i + 1]))
					w[i] = std::abs(w[i]) * sqrt(1.0 + POW2(w[i + 1] / w[i]));
				else w[i] = std::abs(w[i + 1]) * sqrt(1.0 + POW2(w[i] / w[i + 1]));
			}
			for (i = 0; i < n; i++) r[0][i] += w[0] * v[i];
			for (i = 0; i < k; i++)
				rotate(i, r[i][i], -r[i + 1][i]);
			for (i = 0; i < n; i++)
				if (r[i][i] == 0.0) sing = true;
		}
		void rotate(const int i, const Real a, const Real b)
		{
			// Given matrices r[1..n][1..n] and qt[1..n][1..n], carry out a Jacobi rotation on rows
			// i and i + 1 of each Matrix<Real>. a and b are the parameters of the rotation: cos phi = a=pa2 + b2,
			// sin phi = b=pa2 + b2.            
			int j;
			Real c, fact, s, w, y;
			if (a == 0.0) {
				c = 0.0;
				s = (b >= 0.0 ? 1.0 : -1.0);
			}
			else if (std::abs(a) > std::abs(b)) {
				fact = b / a;
				c = SIGN<Real>(1.0 / sqrt(1.0 + (fact * fact)), a);
				s = fact * c;
			}
			else {
				fact = a / b;
				s = SIGN<Real>(1.0 / sqrt(1.0 + (fact * fact)), b);
				c = fact * s;
			}
			for (j = i; j < n; j++) {
				y = r[i][j];
				w = r[i + 1][j];
				r[i][j] = c * y - s * w;
				r[i + 1][j] = s * y + c * w;
			}
			for (j = 0; j < n; j++) {
				y = qt[i][j];
				w = qt[i + 1][j];
				qt[i][j] = c * y - s * w;
				qt[i + 1][j] = s * y + c * w;
			}
		}
	};

	/////////////////////////////////   SVD DECOMPOSITION      /////////////////////////////
	class SVDecompositionSolver
	{
	private:
		int m, n;
		Matrix<Real> u, v;
		Vector<Real> w;
		Real eps, tsh;

	public:
		Vector<Real> getW() { return w; }
		Matrix<Real> getU() { return u; }
		Matrix<Real> getV() { return v; }

	public:
		SVDecompositionSolver(const Matrix<Real>& a) : m(a.RowNum()), n(a.ColNum()), u(a), v(n, n), w(n)
		{
			// Given a Matrix<Real> a[1..m][1..n], this routine computes its singular value decomposition, A = U·W ·V T . 
			// The Matrix<Real> U replaces a on output. 
			// The diagonal Matrix<Real> of singular values W is output as a Vector<Real> w[1..n]. 
			// The Matrix<Real> V (not the transpose V T ) is output as v[1..n][1..n].            
			eps = std::numeric_limits<Real>::epsilon();
			decompose();
			reorder();
			tsh = 0.5 * sqrt(m + n + 1.) * w[0] * eps;
		}

		Real inv_condition() {
			return (w[0] <= 0. || w[n - 1] <= 0.) ? 0. : w[n - 1] / w[0];
		}

		void Solve(const Vector<Real>& b, Vector<Real>& x, Real thresh = -1.)
		{
			// Solve A  x D b for a vector x using the pseudoinverse of A as obtained by SVD. If positive,
			// thresh is the threshold value below which singular values are considered as zero. If thresh is
			// negative, a default based on expected roundoff error is used.
			int i, j, jj;
			Real s;
			if (b.size() != m || x.size() != n) throw("solve bad sizes");
			Vector<Real> tmp(n);
			tsh = (thresh >= 0. ? thresh : 0.5 * sqrt(m + n + 1.) * w[0] * eps);
			for (j = 0; j < n; j++) {
				s = 0.0;
				if (w[j] > tsh) {
					for (i = 0; i < m; i++) s += u[i][j] * b[i];
					s /= w[j];
				}
				tmp[j] = s;
			}
			for (j = 0; j < n; j++) {
				s = 0.0;
				for (jj = 0; jj < n; jj++) s += v[j][jj] * tmp[jj];
				x[j] = s;
			}
		}

		Vector<Real> Solve(const Vector<Real>& b, Real thresh = -1.)
		{
			Vector<Real> x(b.size());
			Solve(b, x, thresh);
			return x;
		}

		// Solves m sets of n equations A  X D B using the pseudoinverse of A. The right-hand sides are
		// input as b[0..n-1][0..m-1], while x[0..n-1][0..m-1] returns the solutions. thresh as above.
		void Solve(const Matrix<Real>& b, Matrix<Real>& x, Real thresh = -1.)
		{
			int i, j, p = b.ColNum();
			if (b.RowNum() != m || x.RowNum() != n || x.ColNum() != p)
				throw("solve bad sizes");

			Vector<Real> xx(n), bcol(m);
			for (j = 0; j < p; j++) {
				for (i = 0; i < m; i++)
					bcol[i] = b[i][j];

				Solve(bcol, xx, thresh);

				for (i = 0; i < n; i++)
					x[i][j] = xx[i];
			}
		}

		// Return the rank of A, after zeroing any singular values smaller than thresh. If thresh is
		// negative, a default value based on estimated roundoff is used.        
		int Rank(Real thresh = -1.) {
			int j, nr = 0;
			tsh = (thresh >= 0. ? thresh : 0.5 * sqrt(m + n + 1.) * w[0] * eps);
			for (j = 0; j < n; j++) if (w[j] > tsh) nr++;
			return nr;
		}

		// Return the nullity of A, after zeroing any singular values smaller than thresh. Default value as above.
		int Nullity(Real thresh = -1.) {
			int j, nn = 0;
			tsh = (thresh >= 0. ? thresh : 0.5 * sqrt(m + n + 1.) * w[0] * eps);
			for (j = 0; j < n; j++) if (w[j] <= tsh) nn++;
			return nn;
		}

		// Gives an orthonormal basis for the range of A as the columns of a returned matrix. thresh as above.
		Matrix<Real> Range(Real thresh = -1.) {
			int i, j, nr = 0;
			Matrix<Real> rnge(m, Rank(thresh));
			for (j = 0; j < n; j++) {
				if (w[j] > tsh) {
					for (i = 0; i < m; i++) rnge[i][nr] = u[i][j];
					nr++;
				}
			}
			return rnge;
		}

		// Gives an orthonormal basis for the nullspace of A as the columns of a returned matrix. thresh as above
		Matrix<Real> Nullspace(Real thresh = -1.) {
			int j, jj, nn = 0;
			Matrix<Real> nullsp(n, Nullity(thresh));
			for (j = 0; j < n; j++) {
				if (w[j] <= tsh) {
					for (jj = 0; jj < n; jj++) nullsp[jj][nn] = v[jj][j];
					nn++;
				}
			}
			return nullsp;
		}

	private:
		void decompose() {
			bool flag;
			int i, its, j, jj, k, l, nm;
			Real anorm, c, f, g, h, s, scale, x, y, z;
			Vector<Real> rv1(n);
			g = scale = anorm = 0.0;
			for (i = 0; i < n; i++) {
				l = i + 2;
				rv1[i] = scale * g;
				g = s = scale = 0.0;
				if (i < m) {
					for (k = i; k < m; k++) scale += std::abs(u[k][i]);
					if (scale != 0.0) {
						for (k = i; k < m; k++) {
							u[k][i] /= scale;
							s += u[k][i] * u[k][i];
						}
						f = u[i][i];
						g = -SIGN(sqrt(s), f);
						h = f * g - s;
						u[i][i] = f - g;
						for (j = l - 1; j < n; j++) {
							for (s = 0.0, k = i; k < m; k++) s += u[k][i] * u[k][j];
							f = s / h;
							for (k = i; k < m; k++) u[k][j] += f * u[k][i];
						}
						for (k = i; k < m; k++) u[k][i] *= scale;
					}
				}
				w[i] = scale * g;
				g = s = scale = 0.0;
				if (i + 1 <= m && i + 1 != n) {
					for (k = l - 1; k < n; k++) scale += std::abs(u[i][k]);
					if (scale != 0.0) {
						for (k = l - 1; k < n; k++) {
							u[i][k] /= scale;
							s += u[i][k] * u[i][k];
						}
						f = u[i][l - 1];
						g = -SIGN(sqrt(s), f);
						h = f * g - s;
						u[i][l - 1] = f - g;
						for (k = l - 1; k < n; k++) rv1[k] = u[i][k] / h;
						for (j = l - 1; j < m; j++) {
							for (s = 0.0, k = l - 1; k < n; k++) s += u[j][k] * u[i][k];
							for (k = l - 1; k < n; k++) u[j][k] += s * rv1[k];
						}
						for (k = l - 1; k < n; k++) u[i][k] *= scale;
					}
				}
				anorm = std::max(anorm, (std::abs(w[i]) + std::abs(rv1[i])));
			}
			for (i = n - 1; i >= 0; i--) {
				if (i < n - 1) {
					if (g != 0.0) {
						for (j = l; j < n; j++)
							v[j][i] = (u[i][j] / u[i][l]) / g;
						for (j = l; j < n; j++) {
							for (s = 0.0, k = l; k < n; k++) s += u[i][k] * v[k][j];
							for (k = l; k < n; k++) v[k][j] += s * v[k][i];
						}
					}
					for (j = l; j < n; j++) v[i][j] = v[j][i] = 0.0;
				}
				v[i][i] = 1.0;
				g = rv1[i];
				l = i;
			}
			for (i = std::min(m, n) - 1; i >= 0; i--) {
				l = i + 1;
				g = w[i];
				for (j = l; j < n; j++) u[i][j] = 0.0;
				if (g != 0.0) {
					g = 1.0 / g;
					for (j = l; j < n; j++) {
						for (s = 0.0, k = l; k < m; k++) s += u[k][i] * u[k][j];
						f = (s / u[i][i]) * g;
						for (k = i; k < m; k++) u[k][j] += f * u[k][i];
					}
					for (j = i; j < m; j++) u[j][i] *= g;
				}
				else for (j = i; j < m; j++) u[j][i] = 0.0;
				++u[i][i];
			}
			for (k = n - 1; k >= 0; k--) {
				for (its = 0; its < 30; its++) {
					flag = true;
					for (l = k; l >= 0; l--) {
						nm = l - 1;
						if (l == 0 || std::abs(rv1[l]) <= eps * anorm) {
							flag = false;
							break;
						}
						if (std::abs(w[nm]) <= eps * anorm) break;
					}
					if (flag) {
						c = 0.0;
						s = 1.0;
						for (i = l; i < k + 1; i++) {
							f = s * rv1[i];
							rv1[i] = c * rv1[i];
							if (std::abs(f) <= eps * anorm) break;
							g = w[i];
							h = pythag(f, g);
							w[i] = h;
							h = 1.0 / h;
							c = g * h;
							s = -f * h;
							for (j = 0; j < m; j++) {
								y = u[j][nm];
								z = u[j][i];
								u[j][nm] = y * c + z * s;
								u[j][i] = z * c - y * s;
							}
						}
					}
					z = w[k];
					if (l == k) {
						if (z < 0.0) {
							w[k] = -z;
							for (j = 0; j < n; j++) v[j][k] = -v[j][k];
						}
						break;
					}
					if (its == 29) throw("no convergence in 30 svdcmp iterations");
					x = w[l];
					nm = k - 1;
					y = w[nm];
					g = rv1[nm];
					h = rv1[k];
					f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
					g = pythag(f, 1.0);
					f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
					c = s = 1.0;
					for (j = l; j <= nm; j++) {
						i = j + 1;
						g = rv1[i];
						y = w[i];
						h = s * g;
						g = c * g;
						z = pythag(f, h);
						rv1[j] = z;
						c = f / z;
						s = h / z;
						f = x * c + g * s;
						g = g * c - x * s;
						h = y * s;
						y *= c;
						for (jj = 0; jj < n; jj++) {
							x = v[jj][j];
							z = v[jj][i];
							v[jj][j] = x * c + z * s;
							v[jj][i] = z * c - x * s;
						}
						z = pythag(f, h);
						w[j] = z;
						if (z) {
							z = 1.0 / z;
							c = f * z;
							s = h * z;
						}
						f = c * g + s * y;
						x = c * y - s * g;
						for (jj = 0; jj < m; jj++) {
							y = u[jj][j];
							z = u[jj][i];
							u[jj][j] = y * c + z * s;
							u[jj][i] = z * c - y * s;
						}
					}
					rv1[l] = 0.0;
					rv1[k] = f;
					w[k] = x;
				}
			}
		}

		void reorder() {
			int i, j, k, s, inc = 1;
			Real sw;
			Vector<Real> su(m), sv(n);
			do { inc *= 3; inc++; } while (inc <= n);
			do {
				inc /= 3;
				for (i = inc; i < n; i++) {
					sw = w[i];
					for (k = 0; k < m; k++) su[k] = u[k][i];
					for (k = 0; k < n; k++) sv[k] = v[k][i];
					j = i;
					while (w[j - inc] < sw) {
						w[j] = w[j - inc];
						for (k = 0; k < m; k++) u[k][j] = u[k][j - inc];
						for (k = 0; k < n; k++) v[k][j] = v[k][j - inc];
						j -= inc;
						if (j < inc) break;
					}
					w[j] = sw;
					for (k = 0; k < m; k++) u[k][j] = su[k];
					for (k = 0; k < n; k++) v[k][j] = sv[k];

				}
			} while (inc > 1);
			for (k = 0; k < n; k++) {
				s = 0;
				for (i = 0; i < m; i++) if (u[i][k] < 0.) s++;
				for (j = 0; j < n; j++) if (v[j][k] < 0.) s++;
				if (s > (m + n) / 2) {
					for (i = 0; i < m; i++) u[i][k] = -u[i][k];
					for (j = 0; j < n; j++) v[j][k] = -v[j][k];
				}
			}
		}

		Real pythag(const Real a, const Real b) {
			Real absa = std::abs(a), absb = std::abs(b);
			return (absa > absb ? absa * sqrt(1.0 + POW2(absb / absa)) :
				(absb == 0.0 ? 0.0 : absb * sqrt(1.0 + POW2(absa / absb))));
		}
	};
} // end namespace
///////////////////////////   ./include/core/CoreUtils.h   ///////////////////////////




namespace MML
{
    namespace MatrixUtils
    {
        static Real Det(const Matrix<Real> &mat)
        {
            Matrix<Real> matCopy(mat);
            LUDecompositionSolver<Real> _solver(matCopy);

            return _solver.det();
        }

        static int Rank(const Matrix<Real> &mat)
        {
            Matrix<Real>     matcopy(mat);

            SVDecompositionSolver svdSolver(matcopy);

            return svdSolver.Rank();
        }
    }
}
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
		static Real NDer1(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			// Error bound ~eps^1/2
			// Note that this estimate of h differs from the best estimate by a factor of sqrt((|f(x)| + |f(x+h)|)/|f''(x)|).
			// Since this factor is invariant under the scaling f -> kf, then we are somewhat justified in approximating it by 1.
			// This approximation will get better as we move to higher orders of accuracy.
			return NDer1(f, x, NDer1_h, error);
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
		static Real NDer2(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			// Error bound ~eps^2/3
			// See the previous discussion to understand determination of h and the error bound.
			// Series[(f[x+h] - f[x-h])/(2*h), {h, 0, 4}]

			return NDer2(f, x, NDer2_h, error);
		}
		static Real NDer2Left(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer2(f, x - 2 * NDer2_h, NDer2_h, error); }
		static Real NDer2Right(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer2(f, x + 2 * NDer2_h, NDer2_h, error); }

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
		static Real NDer4(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			// Error bound ~eps^4/5
			return NDer4(f, x, NDer4_h, error);
		}
		static Real NDer4Left(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer4(f, x - 4 * NDer4_h, NDer4_h, error); }
		static Real NDer4Right(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer4(f, x + 4 * NDer4_h, NDer4_h, error); }

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
		static Real NDer6(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			// Error bound ~eps^6/7
			// Error: h^6f^(7)(x)/140 + 5|f(x)|eps/h
			return NDer6(f, x, NDer6_h, error);
		}
		static Real NDer6Left(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer6(f, x - 5 * NDer6_h, NDer6_h, error); }
		static Real NDer6Right(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer6(f, x + 5 * NDer6_h, NDer6_h, error); }

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
		static Real NDer8(const IRealFunction& f, Real x, Real* error = nullptr)
		{
			// Error bound ~eps^8/9.
			// In Real precision, we only expect to lose two digits of precision while using this formula, at the cost of 8 function evaluations.
			// Error: h^8|f^(9)(x)|/630 + 7|f(x)|eps/h assuming 7 unstabilized additions.
			// Mathematica code to get the error:
			// Series[(f[x+h]-f[x-h])*(4/5) + (1/5)*(f[x-2*h] - f[x+2*h]) + (4/105)*(f[x+3*h] - f[x-3*h]) + (1/280)*(f[x-4*h] - f[x+4*h]), {h, 0, 9}]
			// If we used Kahan summation, we could get the max error down to h^8|f^(9)(x)|/630 + |f(x)|eps/h.

			return NDer8(f, x, NDer8_h, error);
		}
		static Real NDer8Left(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer8(f, x - 6 * NDer8_h, NDer8_h, error); }
		static Real NDer8Right(const IRealFunction& f, Real x, Real* error = nullptr) { return NDer8(f, x + 6 * NDer8_h, NDer8_h, error); }

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
		static inline Real(*Derive)(const IRealFunction& f, Real x, Real* error) = Derivation::NDer4;

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
}

///////////////////////////   ./include/core/Jacobians.h   ///////////////////////////




namespace MML
{
	template<typename VectorFrom, typename VectorTo, int N>
	class Jacobian
	{
	public:
		static MatrixNM<Real, N, N> calc(const IVectorFunction<N>& func, const VectorN<Real, N>& x)
		{
			MatrixNM<Real, N, N> jac;

			for (int i = 0; i < N; ++i)
				for (int j = 0; j < N; ++j)
				{
					jac(i, j) = Derivation::NDer4Partial(func, i, j, x);
				}

			return jac;
		}

		template<int M>
		static MatrixNM<Real, M, N> calc(const IVectorFunctionNM<N, M>& func, const VectorN<Real, N>& x)
		{
			MatrixNM<Real, M, N> jac;

			for (int i = 0; i < M; ++i)
				for (int j = 0; j < N; ++j)
				{
					jac(i, j) = Derivation::NDer4Partial(func, i, j, x);
				}

			return jac;
		}

		static MatrixNM<Real, N, N> calc(const ICoordTransf<VectorFrom, VectorTo, N>& func, const VectorN<Real, N>& x)
		{
			MatrixNM<Real, N, N> jac;

			for (int i = 0; i < N; ++i)
				for (int j = 0; j < N; ++j)
				{
					jac(i, j) = Derivation::NDer4Partial(func.coordTransfFunc(i), j, x);
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

	static Real IntegrateTrap(const IRealFunction& func, const Real a, const Real b, Real req_eps)
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
					return s;
				}
			}
			olds = s;
		}
		return s;
		// TODO - add parameter to signal this situation
		throw IntegrationTooManySteps("qtrap");
	}
	static Real IntegrateTrap(const IRealFunction& func, const Real a, const Real b)
	{
		return IntegrateTrap(func, a, b, Defaults::IntegrateTrapEPS);
	}

	static Real IntegrateSimpson(const IRealFunction& func, const Real a, const Real b, Real req_eps)
	{
		// Returns the integral of the function func from a to b. The parameters EPS can be set to the
		// desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum allowed
		// number of steps. Integration is performed by Simpson’s rule.

		// The routine qsimp will in general be more efficient than qtrap (i.e., require
		// fewer function evaluations) when the function to be integrated has a finite 4th
		// derivative (i.e., a continuous 3rd derivative). The combination of qsimp and its
		// necessary workhorse trapzd is a good one for light-duty work.
		int j;
		Real s, st, ost = 0.0, os = 0.0;

		for (j = 0; j < Defaults::IntegrateSimpMaxSteps; j++)
		{
			st = TrapRefine(func, a, b, j + 1);
			s = (4.0 * st - ost) / 3.0;

			if (j > 5)
				if (std::abs(s - os) < req_eps * std::abs(os) ||
					(s == 0.0 && os == 0.0))
					return s;
			os = s;
			ost = st;
		}
		throw IntegrationTooManySteps("qsimp");
	}
	static Real IntegrateSimpson(const IRealFunction& func, const Real a, const Real b)
	{
		return IntegrateSimpson(func, a, b, Defaults::IntegrateSimpEPS);
	}

	static bool polint(Vector<Real>& xa, Vector<Real>& ya, const Real x, Real& y, Real& dy)
	{
		int i, m, ns = 0;
		Real den, dif, dift, ho, hp, w;

		int n = (int)xa.size();
		Vector<Real> c(n), d(n);
		dif = fabs(x - xa[0]);
		for (i = 0; i < n; i++) {
			if ((dift = fabs(x - xa[i])) < dif) {
				ns = i;
				dif = dift;
			}
			c[i] = ya[i];
			d[i] = ya[i];
		}
		y = ya[ns--];
		for (m = 1; m < n; m++) {
			for (i = 0; i < n - m; i++) {
				ho = xa[i] - x;
				hp = xa[i + m] - x;
				w = c[i + 1] - d[i];
				if ((den = ho - hp) == 0.0)
					// nrerror("Error in routine polint");
					return false;
				den = w / den;
				d[i] = hp * den;
				c[i] = ho * den;
			}
			y += (dy = (2 * (ns + 1) < (n - m) ? c[ns + 1] : d[ns--]));
		}
		return true;
	}

	static void polin2(Vector<Real>& x1a, Vector<Real>& x2a, Matrix<Real>& ya, const Real x1,
		const Real x2, Real& y, Real& dy)
		// Given arrays x1a[1..m] and x2a[1..n] of independent variables, and a submatrix of function
		// values ya[1..m][1..n], tabulated at the grid points defined by x1a and x2a; and given values
		// x1 and x2 of the independent variables; this routine returns an interpolated function value y,
		// and an accuracy indication dy (based only on the interpolation in the x1 direction, however).          
	{
		int j, k;

		int m = (int)x1a.size();
		int n = (int)x2a.size();
		Vector<Real> ymtmp(m), ya_t(n);
		for (j = 0; j < m; j++) {
			for (k = 0; k < n; k++)
				ya_t[k] = ya[j][k];

			polint(x2a, ya_t, x2, ymtmp[j], dy);
		}
		polint(x1a, ymtmp, x1, y, dy);
	}

	static Real IntegrateRomberg(const IRealFunction& func, const Real a, const Real b, Real req_eps)
	{
		// Returns the integral of the function func from a to b. Integration is performed by Romberg’s
		// method of order 2K, where, e.g., K=2 is Simpson’s rule.

		// The routine qromb, along with its required trapzd and polint, is quite
		// powerful for sufficiently smooth (e.g., analytic) integrands, integrated over intervals
		// which contain no singularities, and where the enRealoints are also nonsingular. qromb,
		// in such circumstances, takes many, many fewer function evaluations than either of
		// the routines in x4.2
		const int JMAXP = Defaults::IntegrateRombMaxSteps + 1, K = 5;
		Real ss, dss;
		Vector<Real> s(Defaults::IntegrateRombMaxSteps), h(JMAXP), s_t(K), h_t(K);

		h[0] = 1.0;
		for (int j = 1; j <= Defaults::IntegrateRombMaxSteps; j++)
		{
			s[j - 1] = TrapRefine(func, a, b, j);

			if (j >= K)
			{
				for (int i = 0; i < K; i++) {
					h_t[i] = h[j - K + i];
					s_t[i] = s[j - K + i];
				}

				polint(h_t, s_t, 0.0, ss, dss);

				if (std::abs(dss) <= req_eps * std::abs(ss))
					return ss;
			}

			h[j] = 0.25 * h[j - 1];
		}
		throw IntegrationTooManySteps("qromb");
	}
	static Real IntegrateRomberg(const IRealFunction& func, const Real a, const Real b)
	{
		return IntegrateRomberg(func, a, b, Defaults::IntegrateRombEPS);
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
			case SIMPSON:
				return IntegrateSimpson(f2, y1(x), y2(x));
			case ROMBERG:
				return IntegrateRomberg(f2, y1(x), y2(x));
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
		case SIMPSON:
			return IntegrateSimpson(f1, x1, x2);
		case ROMBERG:
			return IntegrateRomberg(f1, x1, x2);
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

	static inline Real(*Integrate)(const MML::IRealFunction& f, Real a, Real b, Real req_eps) = IntegrateSimpson;

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

		MatrixNM<Real, N, N> jacobian(const VectorN<Real, N>& x) const { return Jacobian<VectorN<Real, N>, VectorN<Real, N>, N>::calc(*this, x); }
	};

	template<int N>
	class VectorFunctionFromStdFunc : public IVectorFunction<N>
	{
		std::function<VectorN<Real, N>(const VectorN<Real, N>&)> _func;
	public:
		VectorFunctionFromStdFunc(std::function<VectorN<Real, N>(const VectorN<Real, N>&)>& inFunc) : _func(inFunc) {}

		VectorN<Real, N>     operator()(const VectorN<Real, N>& x) const { return _func(x); }

		MatrixNM<Real, N, N> jacobian(const VectorN<Real, N>& x) const { return Jacobian<VectorN<Real, N>, VectorN<Real, N>, N>::calc(*this, x); }
	};

	/////////////////////////    VECTOR FUNCTION N -> M      ///////////////////////////////////
	template<int N, int M>
	class VectorFunctionNM : public IVectorFunctionNM<N, M>
	{
		VectorN<Real, M>(*_func)(const VectorN<Real, N>&);
	public:
		VectorFunctionNM(VectorN<Real, M>(*inFunc)(const VectorN<Real, N>&)) : _func(inFunc) {}

		VectorN<Real, M>     operator()(const VectorN<Real, N>& x) const { return _func(x); }

		MatrixNM<Real, M, N> jacobian(const VectorN<Real, N>& x) const { return Jacobian<VectorN<Real, N>, VectorN<Real, M>, N>::calc(*this, x); }
	};

	template<int N, int M>
	class VectorFunctionNMFromStdFunc : public IVectorFunctionNM<N, M>
	{
		std::function<VectorN<Real, M>(const VectorN<Real, N>&)> _func;
	public:
		VectorFunctionNMFromStdFunc(std::function<VectorN<Real, M>(const VectorN<Real, N>&)>& inFunc) : _func(inFunc) {}

		VectorN<Real, M>     operator()(const VectorN<Real, N>& x) const { return _func(x); }

		MatrixNM<Real, M, N> jacobian(const VectorN<Real, N>& x) const { return Jacobian<VectorN<Real, N>, VectorN<Real, M>, N>::calc(*this, x); }
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

		virtual Real getMinX() const { return _minX; }
		virtual Real getMaxX() const { return _maxX; }
		virtual Real getMinY() const { return _minY; }
		virtual Real getMaxY() const { return _maxY; }

		// TODO - double getStartY(double x) const;     // ako surface patch  nije kvadratni
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

		virtual Real getMinX() const { return _minX; }
		virtual Real getMaxX() const { return _maxX; }
		virtual Real getMinY() const { return _minY; }
		virtual Real getMaxY() const { return _maxY; }
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

	// Polynomial interpolation object. Construct with x and y vectors, and the number M of points
	// to be used locally (polynomial order plus one), then call interp for interpolated values.
	struct PolynomInterpRealFunc : RealFunctionInterpolatedBase
	{
		mutable Real dy;

		// The user interface to Poly_interp is virtually the same as for Linear_interp
		// (end of ÷3.1), except that an additional argument in the constructor sets M, the number of points used (the order plus one). 
		PolynomInterpRealFunc(Vector<Real>& xv, Vector<Real>& yv, int m) : RealFunctionInterpolatedBase(xv, yv, m), dy(0.)
		{}

		// Given a value x, and using pointers to data xx and yy, this routine returns an interpolated
		// value y, and stores an error estimate dy. The returned value is obtained by mm-point polynomial
		// interpolation on the subrange xx[jl..jl+mm-1].
		Real rawinterp(int jl, Real x) const
		{
			int i, m, ns = 0;
			Real y, den, dif, dift, ho, hp, w;
			const Real* xa = &xx[jl];
			const Real* ya = &yy[jl];
			Vector<Real> c(mm), d(mm);
			dif = std::abs(x - xa[0]);
			for (i = 0; i < mm; i++) {
				if ((dift = std::abs(x - xa[i])) < dif) {
					ns = i;
					dif = dift;
				}
				c[i] = ya[i];
				d[i] = ya[i];
			}
			y = ya[ns--];
			for (m = 1; m < mm; m++) {
				for (i = 0; i < mm - m; i++) {
					ho = xa[i] - x;
					hp = xa[i + m] - x;
					w = c[i + 1] - d[i];
					if ((den = ho - hp) == 0.0) throw("Poly_interp error");
					den = w / den;
					d[i] = hp * den;
					c[i] = ho * den;
				}

				y += (dy = (2 * (ns + 1) < (mm - m) ? c[ns + 1] : d[ns--]));
				// After each column in the tableau is completed, we decide which correction, c or d, we
				// want to add to our accumulating value of y, i.e., which path to take through the tableau
				// — forking up or down. We do this in such a way as to take the most “straight line”
				// route through the tableau to its apex, updating ns accordingly to keep track of where
				// we are. This route keeps the partial approximations centered (insofar as possible) on
				// the target x. The last dy added is thus the error indication.                
			}
			return y;
		}
	};

	// Diagonal rational function interpolation object. Construct with x and y vectors, and the number
	// m of points to be used locally, then call interp for interpolated values.
	struct RationalInterpRealFunc : RealFunctionInterpolatedBase
	{
		mutable Real dy;
		RationalInterpRealFunc(Vector<Real>& xv, Vector<Real>& yv, int m) : RealFunctionInterpolatedBase(xv, yv, m), dy(0.)
		{}

		// Given a value x, and using pointers to data xx and yy, this routine returns an interpolated value
		// y, and stores an error estimate dy. The returned value is obtained by mm-point diagonal rational
		// function interpolation on the subrange xx[jl..jl+mm-1].        
		Real rawinterp(int jl, Real x) const
		{
			const Real TINY = 1.0e-99;
			int m, i, ns = 0;
			Real y, w, t, hh, h, dd;
			const Real* xa = &xx[jl], * ya = &yy[jl];
			Vector<Real> c(mm), d(mm);
			hh = std::abs(x - xa[0]);
			for (i = 0; i < mm; i++) {
				h = std::abs(x - xa[i]);
				if (h == 0.0) {
					dy = 0.0;
					return ya[i];
				}
				else if (h < hh) {
					ns = i;
					hh = h;
				}
				c[i] = ya[i];
				d[i] = ya[i] + TINY;
			}
			y = ya[ns--];
			for (m = 1; m < mm; m++) {
				for (i = 0; i < mm - m; i++) {
					w = c[i + 1] - d[i];
					h = xa[i + m] - x;
					t = (xa[i] - x) * d[i] / h;
					dd = t - c[i + 1];
					if (dd == 0.0)  // This error condition indicates that the interpolating function has a pole at the requested value of x.
						throw("Error in routine ratint");
					dd = w / dd;
					d[i] = c[i + 1] * dd;
					c[i] = t * dd;
				}
				y += (dy = (2 * (ns + 1) < (mm - m) ? c[ns + 1] : d[ns--]));
			}
			return y;
		}
	};

	// Cubic spline interpolation object. Construct with x and y vectors, and (optionally) values of
	// the first derivative at the endpoints, then call interp for interpolated values.
	struct SplineInterpRealFunc : RealFunctionInterpolatedBase
	{
		Vector<Real> y2;

		SplineInterpRealFunc(Vector<Real>& xv, Vector<Real>& yv, Real yp1 = 1.e99, Real ypn = 1.e99)
			: RealFunctionInterpolatedBase(xv, yv, 2), y2(xv.size())
		{
			sety2(&xv[0], &yv[0], yp1, ypn);
		}

		// This routine stores an array y2[0..n-1] with second derivatives of the interpolating function
		// at the tabulated points pointed to by xv, using function values pointed to by yv. If yp1 and/or
		// ypn are equal to 1  1099 or larger, the routine is signaled to set the corresponding boundary
		// condition for a natural spline, with zero second derivative on that boundary; otherwise, they are
		// the values of the first derivatives at the endpoints.
		void sety2(const Real* xv, const Real* yv, Real yp1, Real ypn)
		{
			int i, k;
			Real p, qn, sig, un;
			int n = (int)y2.size();
			Vector<Real> u(n - 1);
			if (yp1 > 0.99e99)
				y2[0] = u[0] = 0.0;
			else {
				y2[0] = -0.5;
				u[0] = (3.0 / (xv[1] - xv[0])) * ((yv[1] - yv[0]) / (xv[1] - xv[0]) - yp1);
			}
			for (i = 1; i < n - 1; i++) {
				sig = (xv[i] - xv[i - 1]) / (xv[i + 1] - xv[i - 1]);
				p = sig * y2[i - 1] + 2.0;
				y2[i] = (sig - 1.0) / p;
				u[i] = (yv[i + 1] - yv[i]) / (xv[i + 1] - xv[i]) - (yv[i] - yv[i - 1]) / (xv[i] - xv[i - 1]);
				u[i] = (6.0 * u[i] / (xv[i + 1] - xv[i - 1]) - sig * u[i - 1]) / p;
			}
			if (ypn > 0.99e99)
				qn = un = 0.0;
			else {
				qn = 0.5;
				un = (3.0 / (xv[n - 1] - xv[n - 2])) * (ypn - (yv[n - 1] - yv[n - 2]) / (xv[n - 1] - xv[n - 2]));
			}
			y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);
			for (k = n - 2; k >= 0; k--)
				y2[k] = y2[k] * y2[k + 1] + u[k];
		}

		// Given a value x, and using pointers to data xx and yy, and the stored vector of second derivatives
		// y2, this routine returns the cubic spline interpolated value y.        
		Real rawinterp(int jl, Real x) const
		{
			int klo = jl, khi = jl + 1;
			Real y, h, b, a;
			h = xx[khi] - xx[klo];
			if (h == 0.0) throw("Bad input to routine splint");
			a = (xx[khi] - x) / h;
			b = (x - xx[klo]) / h;
			y = a * yy[klo] + b * yy[khi] + ((a * a * a - a) * y2[klo]
				+ (b * b * b - b) * y2[khi]) * (h * h) / 6.0;
			return y;
		}
	};

	// Barycentric rational interpolation object. After constructing the object, call interp for interpolated values. Note that no error estimate dy is calculated.
	struct BaryRatInterpRealFunc : RealFunctionInterpolatedBase
	{
		Vector<Real> w;
		int d;

		// Constructor arguments are x and y vectors of length n, and order d of desired approximation.
		BaryRatInterpRealFunc(Vector<Real>& xv, Vector<Real>& yv, int dd)
			: RealFunctionInterpolatedBase(xv, yv, (int)xv.size()), w(n), d(dd)
		{
			if (n <= d) throw("d too large for number of points in BaryRat_interp");
			for (int k = 0; k < n; k++) {
				int imin = std::max(k - d, 0);
				int imax = k >= n - d ? n - d - 1 : k;
				Real temp = imin & 1 ? -1.0 : 1.0;
				Real sum = 0.0;
				for (int i = imin; i <= imax; i++) {
					int jmax = std::min(i + d, n - 1);
					Real term = 1.0;
					for (int j = i; j <= jmax; j++) {
						if (j == k) continue;
						term *= (xx[k] - xx[j]);
					}
					term = temp / term;
					temp = -temp;
					sum += term;
				}
				w[k] = sum;
			}
		}

		// Use equation (NR 3.4.9) to compute the barycentric rational interpolant. Note that jl is not used
		// since the approximation is global; it is included only for compatibility with Base_interp           
		Real rawinterp(int jl, Real x) const
		{
			Real num = 0, den = 0;
			for (int i = 0; i < n; i++) {
				Real h = x - xx[i];
				if (h == 0.0) {
					return yy[i];
				}
				else {
					Real temp = w[i] / h;
					num += temp * yy[i];
					den += temp;
				}
			}
			return num / den;
		}

		// No need to invoke hunt or locate since the interpolation is global, so override interp to simply
		// call rawinterp directly with a dummy value of jl.        
		Real interp(Real x) {
			return rawinterp(1, x);
		}
	};

	struct BilinInterpScalarFunction2D : public IScalarFunction<2>
	{
		int m, n;
		const Matrix<Real>& y;
		LinearInterpRealFunc x1terp, x2terp;

		BilinInterpScalarFunction2D(Vector<Real>& x1v, Vector<Real>& x2v, Matrix<Real>& ym)
			: m((int)x1v.size()), n((int)x2v.size()), y(ym),
			x1terp(x1v, x1v), x2terp(x2v, x2v) {}

		Real interp(Real x1p, Real x2p) const {
			int i, j;
			Real yy, t, u;
			i = x1terp.cor ? x1terp.hunt(x1p) : x1terp.locate(x1p);
			j = x2terp.cor ? x2terp.hunt(x2p) : x2terp.locate(x2p);
			t = (x1p - x1terp.xx[i]) / (x1terp.xx[i + 1] - x1terp.xx[i]);
			u = (x2p - x2terp.xx[j]) / (x2terp.xx[j + 1] - x2terp.xx[j]);
			yy = (1. - t) * (1. - u) * y[i][j] + t * (1. - u) * y[i + 1][j]
				+ (1. - t) * u * y[i][j + 1] + t * u * y[i + 1][j + 1];
			return yy;
		}

		Real operator()(const VectorN<Real, 2>& x) const
		{
			return interp(x[0], x[1]);
		}
	};

	struct PolynomInterpScalarFunction2D : public IScalarFunction<2>
	{
		int m, n, mm, nn;
		const Matrix<Real>& y;

		mutable Vector<Real> yv;
		mutable PolynomInterpRealFunc x1terp, x2terp;

		PolynomInterpScalarFunction2D(Vector<Real>& x1v, Vector<Real>& x2v, Matrix<Real>& ym, int mp, int np) 
			: m((int)x1v.size()), n((int)x2v.size()), mm(mp), nn(np), y(ym), yv(m),
				x1terp(x1v, yv, mm), x2terp(x2v, x2v, nn) {}

		Real interp(Real x1p, Real x2p) const {
			int i, j, k;
			i = x1terp.cor ? x1terp.hunt(x1p) : x1terp.locate(x1p);
			j = x2terp.cor ? x2terp.hunt(x2p) : x2terp.locate(x2p);
			for (k = i; k < i + mm; k++) {
				// TODO - FIX THIS!!!
				x2terp.yy = y.VectorFromRow(k); // &y[k][0];
				yv[k] = x2terp.rawinterp(j, x2p);
			}
			return x1terp.rawinterp(i, x1p);
		}

		Real operator()(const VectorN<Real, 2>& x) const
		{
			return interp(x[0], x[1]);
		}
	};

	struct SplineInterpScalarFunction2D : public IScalarFunction<2>
	{
		int m, n;
		const Matrix<Real>& y;

		Vector<Real>& x1;
		mutable Vector<Real> yv;

		// TODO 0.9 - HIGH, ovo popraviti
		Vector<SplineInterpRealFunc*> srp;
		std::vector<Vector<Real>> _yVals;

		SplineInterpScalarFunction2D(Vector<Real>& x1v, Vector<Real>& x2v, Matrix<Real>& ym)
			: m((int)x1v.size()), n((int)x2v.size()), y(ym), yv(m), x1(x1v), srp(m)
		{
			_yVals.resize(m);
			for (int i = 0; i < m; i++)
			{
				_yVals[i] = y.VectorFromRow(i);
				srp[i] = new SplineInterpRealFunc(x2v, _yVals[i]);
			}
		}

		~SplineInterpScalarFunction2D() {
			for (int i = 0; i < m; i++) delete srp[i];
		}

		Real interp(Real x1p, Real x2p) const
		{
			for (int i = 0; i < m; i++)
				yv[i] = (*srp[i])(x2p);

			SplineInterpRealFunc scol(x1, yv);

			return scol(x1p);
		}

		Real operator()(const VectorN<Real, 2>& x) const
		{
			return interp(x[0], x[1]);
		}
	};

	class InterpolatedScalarFunction3D : public IScalarFunction<3>
	{
	public:
		InterpolatedScalarFunction3D() {}

		Real operator()(const VectorN<Real, 3>& x) const { return 0.0; }
		virtual Real operator()(Real u, Real w, Real z)
		{
			VectorN<Real, 3> coord{ u,w,z };
			// TODO 0.9 - BIG
			return operator()(coord);
		}
	};

	// class ScalarFunction2Interpolated : public IScalarFunction<2>
	// {
	//     Real (*_func)(const VectorN<Real, 2> &);
	// public:
	//     ScalarFunction2Interpolated( Real (*inFunc)(const VectorN<Real, 2> &) ) : _func(inFunc)    {}

	//     Real operator()(const VectorN<Real, 2> &x) const  { return _func(x); }
	// };

	// class ScalarFunction3Interpolated : public IScalarFunction<3>
	// {
	//     Real (*_func)(const VectorN<Real, 3> &);
	// public:
	//     ScalarFunction3Interpolated( Real (*inFunc)(const VectorN<Real, 3> &) ) : _func(inFunc)    {}

	//     Real operator()(const VectorN<Real, 3> &x) const  { return _func(x); }
	// };

	// Object for interpolating a curve specified by n points in dim dimensions.
	template<int N>
	struct SplineInterpParametricCurve : public IParametricCurve<N>
	{
		int dim, n, in;
		bool cls;
		Matrix<Real> pts;
		Vector<Real> s;
		Vector<Real> ans;
		std::vector<SplineInterpRealFunc*> srp;

		// Constructor. The n   dim matrix ptsin inputs the data points. Input close as 0 for
		// an open curve, 1 for a closed curve. (For a closed curve, the last data point should not
		// duplicate the first — the algorithm will connect them.)
		SplineInterpParametricCurve(const Matrix<Real>& ptsin, bool close = 0)
			: n(ptsin.RowNum()), dim(ptsin.ColNum()), in(close ? 2 * n : n),
			cls(close), pts(dim, in), s(in), ans(dim), srp(dim)
		{
			int i, ii, im, j, ofs;
			Real ss, soff, db, de;
			ofs = close ? n / 2 : 0;
			s[0] = 0.;
			for (i = 0; i < in; i++) {
				ii = (i - ofs + n) % n;
				im = (ii - 1 + n) % n;
				for (j = 0; j < dim; j++) pts[j][i] = ptsin[ii][j];
				if (i > 0) {
					s[i] = s[i - 1] + rad(&ptsin[ii][0], &ptsin[im][0]);
					if (s[i] == s[i - 1]) throw("error in Curve_interp");
					// Consecutive points may not be identical. For a closed curve, the last data
					// point should not duplicate the first.                    
				}
			}
			ss = close ? s[ofs + n] - s[ofs] : s[n - 1] - s[0];
			soff = s[ofs];
			for (i = 0; i < in; i++) s[i] = (s[i] - soff) / ss;
			for (j = 0; j < dim; j++) {
				db = in < 4 ? 1.e99 : fprime(&s[0], &pts[j][0], 1);
				de = in < 4 ? 1.e99 : fprime(&s[in - 1], &pts[j][in - 1], -1);

				Vector<Real> vec = pts.VectorFromRow(j);
				srp[j] = new SplineInterpRealFunc(s, vec, db, de);
			}
		}
		~SplineInterpParametricCurve() {
			for (int j = 0; j < dim; j++) delete srp[j];
		}
		
		Real getMinT() const { return 0.0; }
		Real getMaxT() const { return 1.0; }

		// Interpolate a point on the stored curve. The point is parameterized by t, in the range [0,1].
		// For open curves, values of t outside this range will return extrapolations (dangerous!). For
		// closed curves, t is periodic with period 1
		VectorN<Real, N> interp(Real t) const
		{
			VectorN<Real, N> ans;

			if (cls)
				t = t - floor(t);
			for (int j = 0; j < dim; j++)
				ans[j] = (*srp[j])(t);

			return ans;
		}

		VectorN<Real, N> operator()(Real t) const
		{
			return interp(t);
		}

		// Utility for estimating the derivatives at the endpoints. x and y point to the abscissa and
		// ordinate of the endpoint. If pm is C1, points to the right will be used (left endpoint); if it
		// is -1, points to the left will be used (right endpoint). 
		Real fprime(Real* x, Real* y, int pm) {
			Real s1 = x[0] - x[pm * 1], s2 = x[0] - x[pm * 2], s3 = x[0] - x[pm * 3],
				s12 = s1 - s2, s13 = s1 - s3, s23 = s2 - s3;
			return -(s1 * s2 / (s13 * s23 * s3)) * y[pm * 3] + (s1 * s3 / (s12 * s2 * s23)) * y[pm * 2]
				- (s2 * s3 / (s1 * s12 * s13)) * y[pm * 1] + (1. / s1 + 1. / s2 + 1. / s3) * y[0];
		}

		Real rad(const Real* p1, const Real* p2) {
			Real sum = 0.;
			for (int i = 0; i < dim; i++)
				sum += POW2(p1[i] - p2[i]);
			return sqrt(sum);
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
					ret.Component(i, j) = this->Component(i, j, pos);

			return ret;
		}

		Real GetChristoffelSymbolFirstKind(int i, int j, int k, const VectorN<Real, N>& pos) const
		{
			MetricTensorField<N>& g = *this;

			Real gamma_ijk = 0.0;
			for (int m = 0; m < m; k++)
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

				gamma_ijk += 0.5 * g.Component(i, l) * (coef1 + coef2 - coef3);
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
		MetricTensorCartesian() : MetricTensorField<N>(N, 0) { }

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

		virtual Real Component(int i, int j, const VectorN<Real, 3>& pos) const
		{
			if (i == 0 && j == 0)
				return 1.0;
			else if (i == 1 && j == 1)
				return pos[0] * pos[0];
			else if (i == 2 && j == 2)
				return pos[0] * pos[0] * sin(pos[1]) * sin(pos[1]);
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
		ICoordTransfWithInverse<VectorFrom, VectorTo, N>& _coordTransf;

	public:
		MetricTensorFromCoordTransf(ICoordTransfWithInverse<VectorFrom, VectorTo, N>& inTransf) : _coordTransf(inTransf)
		{ }

		virtual Real Component(int i, int j, const VectorN<Real, N>& pos) const
		{
			Real g_ij = 0.0;
			for (int k = 0; k < N; k++)
			{
				g_ij += Derivation::DerivePartial<N>(_coordTransf.coordTransfFunc(k), i, pos, nullptr) * Derivation::DerivePartial<N>(_coordTransf.coordTransfFunc(k), j, pos, nullptr);
			}
			return g_ij;
		}
	};
}
///////////////////////////   ./include/core/CoordTransf.h   ///////////////////////////




namespace MML
{
	// ovdje dodati translational, rotational, galilean, lorentzian transf
	// SVE su to transformacije koordinata
	template<typename VectorFrom, typename VectorTo, int N>
	class CoordTransf : public virtual ICoordTransf<VectorFrom, VectorTo, N>
	{
	public:
		virtual VectorTo getCovariantBasisVec(int ind, const VectorFrom& pos)
		{
			VectorTo ret;

			for (int i = 0; i < N; i++)
				ret[i] = Derivation::NDer4Partial(this->coordTransfFunc(i), ind, pos);

			return ret;
		}

		virtual VectorTo getUnitVector(int ind, const VectorFrom& pos)
		{
			return getCovariantBasisVec(ind, pos).GetAsUnitVectorAtPos(pos);
		}

		MatrixNM<Real, N, N> jacobian(const VectorN<Real, N>& x)
		{
			return Jacobian<VectorFrom, VectorTo, N>::calc(*this, x);
		}

		VectorTo   transfVecContravariant(const VectorFrom& vec, const VectorFrom& pos)
		{
			VectorFrom ret;
			for (int i = 0; i < N; i++) {
				ret[i] = 0;
				for (int j = 0; j < N; j++)
					ret[i] += Derivation::NDer4Partial(this->coordTransfFunc(i), j, pos, 1e-8) * vec[j];
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
					ret[i] += Derivation::NDer4Partial(this->coordTransfFunc(j), i, pos, 1e-8) * vec[j];
			}
			return ret;
		}
	};

	template<typename VectorFrom, typename VectorTo, int N>
	class CoordTransfWithInverse : public virtual CoordTransf<VectorFrom, VectorTo, N>,
		public virtual ICoordTransfWithInverse<VectorFrom, VectorTo, N>
	{
	public:
		virtual VectorFrom getContravariantBasisVec(int ind, const VectorTo& pos)
		{
			VectorFrom ret;

			for (int i = 0; i < N; i++)
				ret[i] = Derivation::NDer4Partial(this->inverseCoordTransfFunc(ind), i, pos);

			return ret;
		}

		virtual VectorFrom getUnitVectorInverse(int ind, const VectorTo& pos)
		{
			return getContravariantBasisVec(ind, pos).GetAsUnitVectorAtPos(pos);
		}

		VectorTo transfVecCovariant(const VectorFrom& vec, const VectorTo& pos)
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
					ret.Component(i, j) = 0;
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

							ret.Component(i, j) += coef1 * coef2 * tensor.Component(k, l);
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

		IScalarFunction<2>& coordTransfFunc(int i) const { return _func[i]; }
		IScalarFunction<2>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
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

		Real func1(const VectorN<Real, 2>& q) const { return (_transf * q)[0]; }
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

	class CoordTransfCart3DRotationXAxis : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	public:
		Real    _angle;
		MatrixNM<Real, 3, 3>  _transf;
		MatrixNM<Real, 3, 3>  _inverse;

		const ScalarFunctionFromStdFunc<3>  _f1;
		const ScalarFunctionFromStdFunc<3>  _f2;
		const ScalarFunctionFromStdFunc<3>  _f3;

		const ScalarFunctionFromStdFunc<3> _fInverse1;
		const ScalarFunctionFromStdFunc<3> _fInverse2;
		const ScalarFunctionFromStdFunc<3> _fInverse3;

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

		Real func1(const VectorN<Real, 3>& q) const
		{
			return (_transf * q)[0];
		}
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

		const ScalarFunctionFromStdFunc<3> _f1;
		const ScalarFunctionFromStdFunc<3> _f2;
		const ScalarFunctionFromStdFunc<3> _f3;

		const ScalarFunctionFromStdFunc<3> _fInverse1;
		const ScalarFunctionFromStdFunc<3> _fInverse2;
		const ScalarFunctionFromStdFunc<3> _fInverse3;

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

		const ScalarFunctionFromStdFunc<3> _f1;
		const ScalarFunctionFromStdFunc<3> _f2;
		const ScalarFunctionFromStdFunc<3> _f3;

		const ScalarFunctionFromStdFunc<3> _fInverse1;
		const ScalarFunctionFromStdFunc<3> _fInverse2;
		const ScalarFunctionFromStdFunc<3> _fInverse3;

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

	class CoordTransfRectilinear : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	private:
		Vector3Cartesian _base[3];
		Vector3Cartesian _dual[3];

		MatrixNM<Real, 3, 3> _alpha;
		MatrixNM<Real, 3, 3> _transf;

		const ScalarFunctionFromStdFunc<3> _f1;
		const ScalarFunctionFromStdFunc<3> _f2;
		const ScalarFunctionFromStdFunc<3> _f3;

		const ScalarFunctionFromStdFunc<3> _fInverse1;
		const ScalarFunctionFromStdFunc<3> _fInverse2;
		const ScalarFunctionFromStdFunc<3> _fInverse3;

		// TODO - HIGH, TESKO, ovo popraviti i napraviti kako spada
		Real func1(const VectorN<Real, 3>& q) const { return ScalarProd(q, MML::Vector3Cartesian(_dual[0])); }
		Real func2(const VectorN<Real, 3>& q) const { return ScalarProd(q, MML::Vector3Cartesian(_dual[1])); }
		Real func3(const VectorN<Real, 3>& q) const { return ScalarProd(q, MML::Vector3Cartesian(_dual[2])); }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

	public:
		CoordTransfRectilinear(VectorN<Real, 3> b1,
			VectorN<Real, 3> b2,
			VectorN<Real, 3> b3) : _f1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfRectilinear::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfRectilinear::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfRectilinear::func3, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfRectilinear::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfRectilinear::funcInverse2, this, std::placeholders::_1) }),
			_fInverse3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfRectilinear::funcInverse3, this, std::placeholders::_1) })
		{
			_base[0] = b1;
			_base[1] = b2;
			_base[2] = b3;

			Vector3Cartesian cross1 = VectorProd(_base[1], _base[2]);
			_dual[0] = (1 / (ScalarProd(_base[0], cross1))) * cross1;

			Vector3Cartesian cross2 = VectorProd(_base[2], _base[0]);
			_dual[1] = (1 / (ScalarProd(_base[1], cross2))) * cross2;

			Vector3Cartesian cross3 = VectorProd(_base[0], _base[1]);
			_dual[2] = (1 / (ScalarProd(_base[2], cross3))) * cross3;

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

		IScalarFunction<3>& coordTransfFunc(int i) const { return _func[i]; }
		IScalarFunction<3>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }

		// TODO 0.9 - overload covar i contravar transf. funkcije s analiticki izracunatim jakobijanom derivacija
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

		IScalarFunction<3>& coordTransfFunc(int i) const { return _func[i]; }
		IScalarFunction<3>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
	};

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
		//static Real funcInverse3(const VectorN<Real, 3> &q) { return q[2]; }

		inline static ScalarFunction<3> _func[3] = { ScalarFunction<3>{x},
																								 ScalarFunction<3>{y},
																								 ScalarFunction<3>{z}
		};

		inline static ScalarFunction<3> _funcInverse[3] = { ScalarFunction<3>{r},
																												ScalarFunction<3>{phi},
																												ScalarFunction<3>{z}
		};
	public:
		Vector3Cartesian     transf(const Vector3Cylindrical& q) const { return Vector3Cartesian{ x(q), y(q), z(q) }; }
		Vector3Cylindrical   transfInverse(const Vector3Cartesian& q) const { return Vector3Cylindrical{ r(q), phi(q), z(q) }; }

		IScalarFunction<3>& coordTransfFunc(int i) const { return _func[i]; }
		IScalarFunction<3>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
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
		//static Real z(const VectorN<Real, 3> &q) { return q[2]; }

		inline static ScalarFunction<3> _func[3] = { ScalarFunction<3>{r},
																								 ScalarFunction<3>{phi},
																								 ScalarFunction<3>{z}
		};

		inline static ScalarFunction<3> _funcInverse[3] = { ScalarFunction<3>{x},
																												ScalarFunction<3>{y},
																												ScalarFunction<3>{z}
		};
	public:
		Vector3Cylindrical transf(const Vector3Cartesian& q) const { return Vector3Cylindrical{ r(q), phi(q), z(q) }; }
		Vector3Cartesian   transfInverse(const Vector3Cylindrical& q) const { return Vector3Cartesian{ x(q), y(q), z(q) }; }

		IScalarFunction<3>& coordTransfFunc(int i) const { return _func[i]; }
		IScalarFunction<3>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
	};

	static CoordTransfSphericalToCartesian      CoordTransfSpherToCart;
	static CoordTransfCylindricalToCartesian    CoordTransfCylToCart;
	static CoordTransfCartesianToSpherical      CoordTransfCartToSpher;
	static CoordTransfCartesianToCylindrical    CoordTransfCartToCyl;
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

			// TODO - depends on cocontra of metic tensor!!!
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
			case 1: ret = Derivation::NDer1PartialByAll<3>(scalarField, pos, nullptr);
			case 2: ret = Derivation::NDer2PartialByAll<3>(scalarField, pos, nullptr);
			case 4: ret = Derivation::NDer4PartialByAll<3>(scalarField, pos, nullptr);
			case 6: ret = Derivation::NDer6PartialByAll<3>(scalarField, pos, nullptr);
			case 8: ret = Derivation::NDer8PartialByAll<3>(scalarField, pos, nullptr);
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
			case 1: ret = Derivation::NDer1PartialByAll<3>(scalarField, pos, nullptr);
			case 2: ret = Derivation::NDer2PartialByAll<3>(scalarField, pos, nullptr);
			case 4: ret = Derivation::NDer4PartialByAll<3>(scalarField, pos, nullptr);
			case 6: ret = Derivation::NDer6PartialByAll<3>(scalarField, pos, nullptr);
			case 8: ret = Derivation::NDer8PartialByAll<3>(scalarField, pos, nullptr);
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
///////////////////////////   ./include/core/ChebyshevApproximation.h   ///////////////////////////



namespace MML
{
class ChebyshevApproximation {
	int n,m;
	Vector<Real> c;
	Real a,b;

    // Constructor from previously computed coefficients.
    ChebyshevApproximation(Vector<Real> &cc, Real aa, Real bb)
		: n(cc.size()), m(n), c(cc), a(aa), b(bb) {}

    // Set m, the number of coefficients after truncating to an error level thresh, and return the value set
	int setm(Real thresh) {while (m>1 && std::abs(c[m-1])<thresh) m--; return m;}
	
    // Chebyshev fit: Given a function func, lower and upper limits of the interval [a,b], compute and
    // save nn coefficients of the Chebyshev approximation such that func.(x) = sum( ... ), where y and x are related by (5.8.10). 
    // This routine is intended to be called with moderately large n (e.g., 30 or 50), the array of c’s subsequently to be truncated at the smaller value
    // m such that cm and subsequent elements are negligible.
    ChebyshevApproximation(const IRealFunction &func, Real aa, Real bb, int nn=50)
        : n(nn), m(nn), c(n), a(aa), b(bb)
    {
        const Real pi=3.141592653589793;
        int k,j;
        Real fac,bpa,bma,y,sum;
        Vector<Real> f(n);
        bma=0.5*(b-a);
        bpa=0.5*(b+a);
        for (k=0;k<n;k++) {
            y=cos(pi*(k+0.5)/n);
            f[k]=func(y*bma+bpa);
        }
        fac=2.0/n;
        for (j=0;j<n;j++) {
            sum=0.0;
            for (k=0;k<n;k++)
                sum += f[k]*cos(pi*j*(k+0.5)/n);
            c[j]=fac*sum;
        }
    }
    // The method eval has an argument for specifying how many leading coefficients
    // m should be used in the evaluation.
    Real eval(Real x, int m)
    {
        Real d=0.0,dd=0.0,sv,y,y2;
        int j;
        if ((x-a)*(x-b) > 0.0) throw("x not in range in Chebyshev::eval");
        y2=2.0*(y=(2.0*x-a-b)/(b-a));
        for (j=m-1;j>0;j--) {
            sv=d;
            d=y2*d-dd+c[j];
            dd=sv;
        }
        return y*d-dd+0.5*c[0];
    }

    // Return a new Chebyshev object that approximates the derivative of the existing function over
    // the same range [a,b].
    ChebyshevApproximation derivative()
    {
        int j;
        Real con;
        Vector<Real> cder(n);
        cder[n-1]=0.0;
        cder[n-2]=2*(n-1)*c[n-1];
        for (j=n-2;j>0;j--)
            cder[j-1]=cder[j+1]+2*j*c[j];
        con=2.0/(b-a);
        for (j=0;j<n;j++) cder[j] *= con;
        return ChebyshevApproximation(cder,a,b);
    }

    // Return a new Chebyshev object that approximates the indefinite integral of the existing function
    // over the same range [a,b]. The constant of integration is set so that the integral vanishes at a.
    ChebyshevApproximation integral()
    {
        int j;
        Real sum=0.0,fac=1.0,con;
        Vector<Real> cint(n);
        con=0.25*(b-a);
        for (j=1;j<n-1;j++) {
            cint[j]=con*(c[j-1]-c[j+1])/j;
            sum += fac*cint[j];
            fac = -fac;
        }
        cint[n-1]=con*c[n-2]/(n-1);
        sum += fac*cint[n-1];
        cint[0]=2.0*sum;
        return ChebyshevApproximation(cint,a,b);
    }

    // Inverse of routine polycofs in Chebyshev: Given an array of polynomial coefficients d[0..n-1],
    // construct an equivalent Chebyshev object.
    ChebyshevApproximation(Vector<Real> &d)
        : n(d.size()), m(n), c(n), a(-1.), b(1.)
    {
        c[n-1]=d[n-1];
        c[n-2]=2.0*d[n-2];
        for (int j=n-3;j>=0;j--) {
            c[j]=2.0*d[j]+c[j+2];
            for (int i=j+1;i<n-2;i++) {
                    c[i] = (c[i]+c[i+2])/2;
            }
            c[n-2] /= 2;
            c[n-1] /= 2;
        }
    }
    // Polynomial coefficients from a Chebyshev fit.
    Vector<Real> polycofs(int m)
    {
        int k,j;
        Real sv;
        Vector<Real> d(m),dd(m);
        for (j=0;j<m;j++) d[j]=dd[j]=0.0;
        d[0]=c[m-1];
        for (j=m-2;j>0;j--) {
            for (k=m-j;k>0;k--) {
                sv=d[k];
                d[k]=2.0*d[k-1]-dd[k];
                dd[k]=sv;
            }
            sv=d[0];
            d[0] = -dd[0]+c[j];
            dd[0]=sv;
        }
        for (j=m-1;j>0;j--) d[j]=d[j-1]-dd[j];
        d[0] = -dd[0]+0.5*c[0];
        return d;
    }
};

} // namespace MML

///////////////////////////   ./include/core/FunctionSpace.h   ///////////////////////////




namespace MML
{
	// razmisliti o complex verziji
	template<int N>
	class OrthogonalFunctionsSpaceN : public HilbertSpace<Real, VectorN<Real, N>>
	{
	protected:
		Real _low, _upp;
	public:
		OrthogonalFunctionsSpaceN(Real low, Real upp) : _low(low), _upp(upp) {}

		virtual Real getLeadCoef(int i) = 0;

		virtual const IRealFunction& getBasisFunc(int i) = 0;
		virtual const IRealFunction& getWeightFunc() = 0;

		virtual VectorN<Real, N> getRepr(const IRealFunction& f) {
			VectorN<Real, N> ret;

			for (int i = 0; i < N; i++) {
				RealFunctionFromStdFunc func_to_integrate([&](Real x) {
					return getWeightFunc()(x) * f(x) * getBasisFunc(i)(x);
					});

				ret[i] = getLeadCoef(i) * Integrate(func_to_integrate, _low, _upp, 1e-5);
			}

			return ret;
		}

		virtual Real  scal_prod(const VectorN<Real, N>& a, const  VectorN<Real, N>& b) const { return a.ScalarProductCartesian(b); }
	};

	class HermitianFunctionSpace5 : public OrthogonalFunctionsSpaceN<5>
	{
	private:
		const static inline RealFunction _weightFunc = RealFunction([](Real x) { return exp(-x * x); });

		const static inline RealFunction _basisFunctions[5] = {
				RealFunction([](Real x) { return (Real)std::hermite(0,x); }),
				RealFunction([](Real x) { return (Real)std::hermite(1,x); }),
				RealFunction([](Real x) { return (Real)std::hermite(2,x); }),
				RealFunction([](Real x) { return (Real)std::hermite(3,x); }),
				RealFunction([](Real x) { return (Real)std::hermite(4,x); })
		};
	public:
		HermitianFunctionSpace5() : OrthogonalFunctionsSpaceN(Constants::NegativeInf, Constants::PositiveInf) {}
		HermitianFunctionSpace5(Real low, Real upp) : OrthogonalFunctionsSpaceN(low, upp) {}

		Real getLeadCoef(int i) { return 1 / (std::pow(2, i) * Functions::Factorial(i) * sqrt(Constants::PI)); }

		const IRealFunction& getBasisFunc(int i) { return _basisFunctions[i]; }
		const IRealFunction& getWeightFunc() { return _weightFunc; }
	};

	class LegendreFunctionSpace5 : public OrthogonalFunctionsSpaceN<5>
	{
	private:
		const static inline RealFunction _weightFunc = RealFunction([](Real x) { return (Real)1.0; });

		const static inline RealFunction _basisFunctions[5] = {
				RealFunction([](Real x) { return (Real)std::legendre(0,x); }),
				RealFunction([](Real x) { return (Real)std::legendre(1,x); }),
				RealFunction([](Real x) { return (Real)std::legendre(2,x); }),
				RealFunction([](Real x) { return (Real)std::legendre(3,x); }),
				RealFunction([](Real x) { return (Real)std::legendre(4,x); })
		};
	public:
		LegendreFunctionSpace5() : OrthogonalFunctionsSpaceN(-1.0, 1.0) {}

		Real getLeadCoef(int i) { return i * 0.5; }

		const IRealFunction& getBasisFunc(int i) { return _basisFunctions[i]; }
		const IRealFunction& getWeightFunc() { return _weightFunc; }
	};

	class LaguerreFunctionSpace5 : public OrthogonalFunctionsSpaceN<5>
	{
	private:
		Real _alpha;

		Real weight_func(Real x) { return exp(-x) * pow(x, _alpha); }
		RealFunctionFromStdFunc _weightFunc = std::function<Real(Real)>{ std::bind(&LaguerreFunctionSpace5::weight_func, this, std::placeholders::_1) };

		const static inline RealFunction _basisFunctions[5] = {
				RealFunction([](Real x) { return (Real)std::laguerre(0,x); }),
				RealFunction([](Real x) { return (Real)std::laguerre(1,x); }),
				RealFunction([](Real x) { return (Real)std::laguerre(2,x); }),
				RealFunction([](Real x) { return (Real)std::laguerre(3,x); }),
				RealFunction([](Real x) { return (Real)std::laguerre(4,x); })
		};
	public:
		LaguerreFunctionSpace5(Real inAlpha) : OrthogonalFunctionsSpaceN(0.0, Constants::PositiveInf)
		{
			if (inAlpha <= -1.0)
				throw std::runtime_error("Alpha must bi bigger than -1");

			_alpha = inAlpha;
		}

		Real getLeadCoef(int i) { return Functions::Factorial(i) / std::tgamma(i + _alpha + 1); }

		const IRealFunction& getBasisFunc(int i) { return _basisFunctions[i]; }
		const IRealFunction& getWeightFunc() { return _weightFunc; }
	};

	// TODO - promijeniti u dinamicki vektor!!!, ali ostaviti N template param
	template<int N>
	class DiscretizedFunctionSpaceN : public VectorSpace<Real, VectorN<Real, N>>
	{
	private:
		Real _low, _upp;
	public:
		DiscretizedFunctionSpaceN(Real low, Real upp) : _low(low), _upp(upp) {}

		// each function is a vector of values, of that function at discretized points
		VectorN<Real, N> getRepr(const IRealFunction& f)
		{
			VectorN<Real, N> ret;
			for (int i = 0; i < N; i++)
			{
				Real x = _low + i * (_upp - _low) / (N - 1);
				ret[i] = f(x);
			}
			return ret;
		}

		// crucial, then you can define operators on functions, but actually vectors, and applying operator to function is multiplying matrix with vector
	};
}
///////////////////////////   ./include/core/Curves.h   ///////////////////////////



namespace MML
{
	namespace Curves2D
	{
		////////////////////////////////             PLANAR CURVES                  //////////////////////////////////
		class Circle2DCurve : public IParametricCurve<2>
		{
			Real _radius;
		public:
			Circle2DCurve() : _radius(1) {}
			Circle2DCurve(Real radius) : _radius(radius) {}

			Real getMinT() const { return 0.0; }
			Real getMaxT() const { return 2 * Constants::PI; }

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
		/////////////////////////////////             SPACE CURVES                  ///////////////////////////////////
		class LineCurve : public IParametricCurve<3>
		{
			Line3D  _line;
			Real _minT;
			Real _maxT;
		public:
			LineCurve(Real minT, Real maxT, const Point3Cartesian& pnt, const Vector3Cartesian dir) : _line(pnt, dir), _minT(minT), _maxT(maxT) {}

			Real getMinT() const { return _minT; }
			Real getMaxT() const { return _maxT; }

			VectorN<Real, 3> operator()(Real t) const
			{
				if (t < _minT || t > _maxT)
					throw std::invalid_argument("LineCurve: t is out of range.");

				auto pnt = _line.PointOnLine(t);
				return VectorN<Real, 3>{pnt.X(), pnt.Y(), pnt.Z()};
			}
		};

		// TODO 0.9 - add SquarePathXY

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
	}
}

///////////////////////////   ./include/core/Surfaces.h   ///////////////////////////


namespace MML
{
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


	}
}

///////////////////////////   ./include/core/DiracDeltaFunction.h   ///////////////////////////


namespace MML
{
	// TODO 0.9 - implement step function
	// TODO 0.9 - implement heaviside function  
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
		ODESystemWithJacobian() { }
		ODESystemWithJacobian(int n,
			void (*inFunc)(Real, const Vector<Real>&, Vector<Real>&),
			void (*inFuncJac)(const Real t, const Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx)
		) : ODESystem(n, inFunc), _funcJac(inFuncJac) { }

		void jacobian(const Real t, Vector<Real>& x, Vector<Real>& dxdt, Matrix<Real>& dydx)
		{
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
			_xval.Resize(maxSteps + 1);
			_yval.Resize(dim, maxSteps + 1);
		}

		template<int N>
		ParametricCurveInterpolated<N> getSolutionAsParametricCurve() const
		{
			ParametricCurveInterpolated<N> curve(_xval, _yval);
			return curve;
		}

		template<int N>
		SplineInterpParametricCurve<N> getSolutionAsSplineParametricCurve() const
		{
			return SplineInterpParametricCurve<N>(_yval);
		}

		LinearInterpRealFunc getSolutionAsLinearInterp(int component) const
		{
			Vector<Real> xsave = _xval;
			Vector<Real> ysave = _yval.VectorFromRow(component);

			return LinearInterpRealFunc(xsave, ysave);
		}
		PolynomInterpRealFunc getSolutionAsPolynomInterp(int component, int polynomDegree) const
		{
			Vector<Real> xsave = _xval;
			Vector<Real> ysave = _yval.VectorFromRow(component);

			return PolynomInterpRealFunc(xsave, ysave, polynomDegree);
		}
		SplineInterpRealFunc getSolutionAsSplineInterp(int component) const
		{
			Vector<Real> xsave = _xval;
			Vector<Real> ysave = _yval.VectorFromRow(component);

			return SplineInterpRealFunc(xsave, ysave);
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
			xval.Resize(numSteps + 1);
			yval.Resize(dim, numSteps + 1);
		}
	};
}
///////////////////////////   ./include/core/CoordSystem.h   ///////////////////////////



namespace MML
{
	// RectilinearCartesianOrthogonal
	class CoordSystemOrthogonalCartesian
	{
		Vector3Cartesian _base[3];

	public:
		CoordSystemOrthogonalCartesian(Vector3Cartesian b1, Vector3Cartesian b2, Vector3Cartesian b3)
		{
			_base[0] = b1;
			_base[1] = b2;
			_base[2] = b3;
		}

		bool isOrthogonal()
		{
			Real a01 = ScalarProd(_base[0], _base[1]);
			Real a02 = ScalarProd(_base[0], _base[2]);
			Real a12 = ScalarProd(_base[1], _base[2]);

			return sqrt(a01 * a01 + a02 * a02 + a12 * a12) < 1e-6;
		}
	};

	class CoordSystemObliqueCartesian
	{
	public:
		Vector3Cartesian _base[3];
		VectorN<Real, 3> _dual[3];

		MatrixNM<Real, 3, 3> _alpha;

		MatrixNM<Real, 3, 3> _transf;
		MatrixNM<Real, 3, 3> _inv;

	public:
		CoordSystemObliqueCartesian(Vector3Cartesian b1, Vector3Cartesian b2, Vector3Cartesian b3)
		{
			_base[0] = b1;
			_base[1] = b2;
			_base[2] = b3;

			Vector3Cartesian cross1 = VectorProd(_base[1], _base[2]);
			_dual[0] = (1 / (ScalarProd(_base[0], cross1))) * cross1;

			Vector3Cartesian cross2 = VectorProd(_base[2], _base[0]);
			_dual[1] = (1 / (ScalarProd(_base[1], cross2))) * cross2;

			Vector3Cartesian cross3 = VectorProd(_base[0], _base[1]);
			_dual[2] = (1 / (ScalarProd(_base[2], cross3))) * cross3;

			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					_alpha(i, j) = _base[i][j];
					_transf(i, j) = _base[j][i];     // transponirano
				}
			}

			_inv = _transf.GetInverse();

		}

		VectorN<Real, 3> Base(int i) { return _base[i]; }
		VectorN<Real, 3> Dual(int i) { return _dual[i]; }
	};

	class ReferentialFrame3D
		// get local coord in parent coord
	{ };

	class InertialFrame3D : public ReferentialFrame3D
	{
		// u odnosu na drugi referential frame ima samo konstantu brzinu
	};

	class MovingFrame3D : public ReferentialFrame3D
	{
		// getOriginPositionAtTime - vraca poziciju u odnosu na ReferentialFrame3D
		// getSpeedAtTime - vraca brzinu u odnosu na ReferentialFrame3D
	};

	// TODO - CircleOrbitingCartesianFrame3DToSpherical
	// ovo je referentni frame koji se rotira oko nekog centra mase, i treba ga zamisliti kao kocku koja rotira oko CM
	class CircleOrbitingCartesianFrame3DToCartesian : public MovingFrame3D
	{
		// koristimo Cartesian sustav - vraca pozicije u odnosu na CM oko kojeg orbitira U CARTESIAN KOORDINATAMA
		// što ukoliko parametre orbite ne zelim u Cartesian sustavu? - nova klasa
	public:
		Real _radius;
		Real _speed;
		Real _period;
		Real _angle_at_t0;
		// axis, za 3d slucaj, kad je ravnina orbite zakrenuta
		Vector3Cartesian _axis;

		ReferentialFrame3D _parentFrame;

		CircleOrbitingCartesianFrame3DToCartesian(ReferentialFrame3D parentFrame, Real radius, Real period)
		{
			_radius = radius;
			// _speed = speed; // izracunati
			_period = period;
			_parentFrame = parentFrame;
			_angle_at_t0 = 0;
		}
		CircleOrbitingCartesianFrame3DToCartesian(ReferentialFrame3D parentFrame, Real radius, Real period, Real angle_at_t0)
		{
			_radius = radius;
			_period = period;
			_angle_at_t0 = angle_at_t0;
			_parentFrame = parentFrame;
		}
		CircleOrbitingCartesianFrame3DToCartesian(ReferentialFrame3D parentFrame, Real radius, Real period, Real orbit_inclination, Real angle_at_t0)
		{
			_radius = radius;
			_period = period;
			_angle_at_t0 = angle_at_t0;
			// calc axis based on inclination
			_parentFrame = parentFrame;
		}

		// OVO JE KLJUCNA FUNKCIJA
		Vector3Cartesian GetOriginPositionAtTime(Real t)
		{
			// calculate rotational evolution of position of center of mass
			Real angle = _angle_at_t0 + 2 * Constants::PI * t / _period;
			// in z-plane!
			Vector3Cartesian CM_pos({ _radius * cos(angle), _radius * sin(angle), 0 });

			return CM_pos;
		}
		// ovo ce morati svaka!
		Vector3Cartesian GetLocalPositionAtTimeInCartesian(Vector3Cartesian pos, Real t)
		{
			// calculate evolution of position of center of mass
			Vector3Cartesian CM_pos = GetOriginPositionAtTime(t);

			// add local coordinates to CM position
			// BITNA PRETPOSTAVKA - kako naš sustav rotira oko CM, njegova apsolutna orijentacije se ne mijenja
			// ie, Zemljina (lokalna) os rotacije je jednom nagnuta OD Sunca, a za sest mjeseci nagnuta PREMA Suncu
			return CM_pos + pos;
		}
		// za danu tocku zadanu Cartesian koordinatama u lokalnom sustavu
		// vraca Cartesian poziciju (u LOKALNOM frameu u odnosu na CM), nakon vremena T
		Vector3Cartesian GetLocalPositionAtTime(Vector3Cartesian pos, Real t)
		{
			// calculate evolution of position of center of mass
			Vector3Cartesian CM_pos = GetOriginPositionAtTime(t);

			// add local coordinates to CM position
			// BITNA PRETPOSTAVKA - kako naš sustav rotira oko CM, njegova apsolutna orijentacije se ne mijenja
			// ie, Zemljina (lokalna) os rotacije je jednom nagnuta OD Sunca, a za sest mjeseci nagnuta PREMA Suncu
			return CM_pos + pos;
		}
		Vector3Cartesian GetLocalPositionAtTime(Vector3Spherical pos, Real t)
		{
			// calculate evolution of position of center of mass
			Vector3Cartesian CM_pos = GetOriginPositionAtTime(t);

			// transform given spherical coordinates to cartesian
			// add local coordinates to CM position        
			return Vector3Cartesian({ 0,0,0 });
		}
	};

	class RotatingFrame3D : public ReferentialFrame3D
	{
		// ovaj frame koristi cilindrični sustav (generalna rotacija oko osi)
	public:
		Real _period;
		Real _angle_at_t0;
		VectorN<Real, 3> _axis;     // pretpostavljamo z-axis za pocetak
		ReferentialFrame3D _parentFrame;

		RotatingFrame3D(ReferentialFrame3D parentFrame, Real period, VectorN<Real, 3> axis)
		{
			_period = period;
			_axis = axis;
			_parentFrame = parentFrame;
		}
		VectorN<Real, 3> GetPositionAtTime(Vector3Cylindrical pos, Real t)
		{
			Real angle = 2 * Constants::PI * t / _period;
			return VectorN<Real, 3>({ _axis[0] * cos(angle), _axis[1] * sin(angle), 0 });
			// TODO
		}
	};

	class SphericalRotatingFrame : public RotatingFrame3D
	{
		// ovaj frame radi sa spherical koordinatama
	public:
		SphericalRotatingFrame(ReferentialFrame3D parentFrame, Real period, VectorN<Real, 3> axis)
			: RotatingFrame3D(parentFrame, period, axis) {}

		VectorN<Real, 3> GetPositionAtTime(Vector3Spherical pos, Real t)
		{
			return VectorN<Real, 3>({ 0,0,0 });
		}
	};

	class HardSphereRotatingFrame : public RotatingFrame3D
	{
		// lokalne koordinate - lat, long, h 
		// vraca - spherical
		// ima svoj CENTAR MASE u sredini sfere, i u odnosu na njega vraća pozicije
		// koje su usglasene s axisom rotacije (lat, long)
	public:
		double _radius;

		HardSphereRotatingFrame(ReferentialFrame3D parentFrame, Real radius, Real period, VectorN<Real, 3> axis)
			: _radius(radius), RotatingFrame3D(parentFrame, period, axis) {}

		Vector3Cartesian GetOriginPositionAtTime(Real t)
		{
			return Vector3Cartesian({ 0,0,0 });       // ne mice se!!!!
		}

		// za danu tocku u lokalnom sustavu, vraca spherical poziciju (u LOKALNOM frameu u odnosu na CM), nakon vremena T
		Vector3Spherical GetLocalPositionAtTime(Real latitude, Real longitude, Real height, Real t)
		{
			// taking into consideration rotation of the sphere
			// calculate position after time t

			// pretvoriti u pravi spherical vektor
			// 
			return Vector3Spherical({ 0,0,0 });
		}

		Vector3Cartesian GetLocalPositionAtTimeInCartesianParentFrame(Real latitude, Real longitude, Real height, Real t)
		{
			// calculate evolution of position of center of mass
			Vector3Cartesian CM_pos = GetOriginPositionAtTime(t);

			Vector3Spherical localSphPos = GetLocalPositionAtTime(latitude, longitude, height, t);
			Vector3Cartesian localCartPos = CoordTransfSpherToCart.transf(localSphPos);
			
			return CM_pos + localCartPos;
		}
	};

	class HardSphereLocalCartesian : public InertialFrame3D
	{
		// u ctor dobije ref na HardSphereRotatingFrame, ITOCNO ODREDJENU TOCKU NA SFERI!!!
		// ima smisla - gleda nakon deltaT gdje je pozicija tocke u jednom i drugom
		// kosi hitac zadan u lokalnom kartezijevom, i izracunam
		// onda vidim gdje je taj lokalni kartezije u trenutku deltaT, i da li se 
		// slaze TRENUTNA tocka (x,y,z) di je sletio hitac, s onom kako sam izracunao
	};

	// da li mi treba Local3D koji za parenta ima HardSphereRotatingFrameToSpherical?
	// lokalni sustav, baziran na TOCNO ODREDJENOJ TOCKI SFERE, s x, y i z
	// za njega NE TREBA davati lat, long i h jer vec ima, a x, y i z transformira lokalno

 
}

///////////////////////////   ./include/core/Fields.h   ///////////////////////////




namespace MML::Fields
{
	////////////////////             INVERSE RADIAL FIELD                /////////////////
	// Potential fields
	static Real InverseRadialPotentialFieldCart(const VectorN<Real, 3>& x) { return 1.0 / x.NormL2(); }
	static Real InverseRadialPotentialFieldCart(Real constant, const VectorN<Real, 3>& x) { return constant / x.NormL2(); }
	static Real InverseRadialPotentialFieldSpher(const VectorN<Real, 3>& x) { return 1.0 / x[0]; }
	static Real InverseRadialPotentialFieldSpher(Real constant, const VectorN<Real, 3>& x) { return constant / x[0]; }
	static Real InverseRadialPotentialFieldCyl(const VectorN<Real, 3>& x) { return 1.0 / sqrt(x[0] * x[0] + x[2] * x[2]); }
	static Real InverseRadialPotentialFieldCyl(Real constant, const VectorN<Real, 3>& x) { return constant / sqrt(x[0] * x[0] + x[2] * x[2]); }

	// Force fields
	static VectorN<Real, 3> InverseRadialPotentialForceFieldCart(const VectorN<Real, 3>& x)     // TODO - ovo u Vector3Cartesian
	{
		return (-1) * x / std::pow((x).NormL2(), 3);
	}
	static VectorN<Real, 3> InverseRadialPotentialForceFieldCart(Real constant, const VectorN<Real, 3>& x)
	{
		return -constant * x / std::pow((x).NormL2(), 3);
	}
	// can't be Vector3Spherical, because in that case we can't form a ScalarFunction<3> out of it :()
	static VectorN<Real, 3> InverseRadialPotentialForceFieldSph(const VectorN<Real, 3>& x)
	{
		return VectorN<Real, 3>{-1 / (x[0] * x[0]), 0.0, 0.0};
	}
	static VectorN<Real, 3> InverseRadialPotentialForceFieldSph(Real constant, const VectorN<Real, 3>& x)
	{
		return VectorN<Real, 3>{-constant / (x[0] * x[0]), 0.0, 0.0};
	}

	class InverseRadialFieldCart : public IScalarFunction<3>
	{
	protected:
		Real _constant;
	public:
		InverseRadialFieldCart() : _constant(-1.0) {}
		InverseRadialFieldCart(Real constant) : _constant(constant) {}

		Real operator()(const VectorN<Real, 3>& x) const { return _constant * Fields::InverseRadialPotentialFieldCart(x); }
	};
	class InverseRadialFieldSpher : public IScalarFunction<3>
	{
	protected:
		Real _constant;
	public:
		InverseRadialFieldSpher() : _constant(-1.0) {}
		InverseRadialFieldSpher(Real constant) : _constant(constant) {}

		Real operator()(const VectorN<Real, 3>& x) const { return _constant * InverseRadialPotentialFieldSpher(x); }
	};

	class InverseRadialForceFieldCart : public IVectorFunction<3>
	{
	private:
		Real _constant;
	public:
		InverseRadialForceFieldCart() : _constant(-1.0) {}
		InverseRadialForceFieldCart(Real constant) : _constant(constant) {}

		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const { return _constant * Fields::InverseRadialPotentialForceFieldCart(x); }
	};
	class InverseRadialForceFieldSpher : public IVectorFunction<3>
	{
	private:
		Real _constant;
	public:
		InverseRadialForceFieldSpher() : _constant(-1.0) {}
		InverseRadialForceFieldSpher(Real constant) : _constant(constant) {}

		VectorN<Real, 3> operator()(const VectorN<Real, 3>& x) const { return _constant * InverseRadialPotentialForceFieldSph(x); }
	};
	// TODO - dodati multibody gravity
	// TODO - dodati InfiniteLineCurrentField
	// TODO - dodati InfinitePlaneCurrentField

}

///////////////////////////   ./include/core/PermutationGroup.h   ///////////////////////////



namespace MML
{
	// TODO - verify this!!!
	class PermutationGroup {
	public:
		PermutationGroup() {}

		PermutationGroup(const std::vector<std::vector<int>>& cycles) {
			for (const auto& cycle : cycles) {
				add_cycle(cycle);
			}
		}

		void add_cycle(const std::vector<int>& cycle) {
			m_cycles.push_back(cycle);
		}

		std::vector<int> apply(const std::vector<int>& perm) const {
			std::vector<int> result = perm;
			for (const auto& cycle : m_cycles) {
				for (size_t i = 0; i < cycle.size(); ++i) {
					result[cycle[i]] = perm[cycle[(i + 1) % cycle.size()]];
				}
			}
			return result;
		}

		PermutationGroup inverse() const {
			std::vector<std::vector<int>> inverse_cycles;
			for (const auto& cycle : m_cycles) {
				std::vector<int> inverse_cycle(cycle.size());
				for (size_t i = 0; i < cycle.size(); ++i) {
					inverse_cycle[cycle[i]] = cycle[(i + 1) % cycle.size()];
				}
				std::reverse(inverse_cycle.begin(), inverse_cycle.end());
				inverse_cycles.push_back(inverse_cycle);
			}
			return PermutationGroup(inverse_cycles);
		}

	private:
		std::vector<std::vector<int>> m_cycles;
	};
}
///////////////////////////   ./include/core/Serializer.h   ///////////////////////////


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
///////////////////////////   ./include/core/Visualizer.h   ///////////////////////////



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
///////////////////////////   ./include/core/ConsolePrinter.h   ///////////////////////////

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
		class HelperWorkIntegral : public IRealFunction
		{
			const IScalarFunction<N>& _potential;
			const IParametricCurve<N>& _curve;
		public:
			HelperWorkIntegral(const IScalarFunction<N>& potentialField, const IParametricCurve<N>& curve) : _potential(potentialField), _curve(curve) {}

			Real operator()(Real t) const
			{
				auto tangent_vec = Derivation::DeriveCurve<N>(_curve, t, nullptr);
				auto gradient = ScalarFieldOperations::GradientCart(_potential, _curve(t));

				return tangent_vec.ScalarProductCartesian(gradient);
			}
		};
		template<int N>
		class HelperLineIntegral : public IRealFunction
		{
			const IVectorFunction<N>& _vector_field;
			const IParametricCurve<N>& _curve;
		public:
			HelperLineIntegral(const IVectorFunction<N>& vectorField, const IParametricCurve<N>& curve) : _vector_field(vectorField), _curve(curve) {}

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

		static Real WorkIntegral(const IScalarFunction<3>& potentialField, const IParametricCurve<3>& curve, const Real a, const Real b, const Real eps = Defaults::WorkIntegralPrecision)
		{
			HelperWorkIntegral helper(potentialField, curve);

			return IntegrateTrap(helper, a, b, eps);
		}

		static Real LineIntegral(const IVectorFunction<3>& vectorField, const IParametricCurve<3>& curve, const Real a, const Real b, const Real eps = Defaults::LineIntegralPrecision)
		{
			HelperLineIntegral helper(vectorField, curve);

			return IntegrateTrap(helper, a, b, eps);
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

	public:
		static Real VolumeIntegral(const IScalarFunction<3>& scalarField, const SolidSurfaces3D& solid)
		{
			return 0.0;
		}
	};
} // end namespace
///////////////////////////   ./include/algorithms/EigenSystemSolvers.h   ///////////////////////////


// Given the eigenvalues d[0..n-1] and (optionally) the eigenvectors v[0..n-1][0..n-1] as determined by Jacobi (÷11.1) or tqli (÷11.4), 
// this routine sorts the eigenvalues into descending order and rearranges the columns of v correspondingly. The method is straight insertion.
static void eigsrt(MML::Vector<Real>& d, MML::Matrix<Real>* v = NULL)
{
	int k;
	int n = (int)d.size();
	for (int i = 0; i < n - 1; i++)
	{
		Real p = d[k = i];
		for (int j = i; j < n; j++)
			if (d[j] >= p)
				p = d[k = j];
		if (k != i)
		{
			d[k] = d[i];
			d[i] = p;
			if (v != NULL)
				for (int j = 0; j < n; j++)
				{
					p = (*v)[j][i];
					(*v)[j][i] = (*v)[j][k];
					(*v)[j][k] = p;
				}
		}
	}
}

namespace MML
{
	// Computes all eigenvalues and eigenvectors of a real symmetric matrix by Jacobi’s method.
	struct SymmMatEigenSolverJacobi
	{
		const int n;
		Matrix<Real> a, v;
		Vector<Real> d;
		int nrot;
		const Real EPS;

		// Computes all eigenvalues and eigenvectors of a real symmetric matrix a[0..n-1][0..n-1].
		// On output, d[0..n-1] contains the eigenvalues of a sorted into descending order, while
		// v[0..n-1][0..n-1] is a matrix whose columns contain the corresponding normalized eigenvectors. nrot contains the number of Jacobi rotations that were required. Only the upper
		// triangle of a is accessed.
		SymmMatEigenSolverJacobi(const Matrix<Real>& aa) : n(aa.RowNum()), a(aa), v(n, n), d(n), nrot(0),
			EPS(std::numeric_limits<Real>::epsilon())
		{
			int i, j, ip, iq;
			Real tresh, theta, tau, t, sm, s, h, g, c;
			Vector<Real> b(n), z(n);
			for (ip = 0; ip < n; ip++)
			{
				for (iq = 0; iq < n; iq++)
					v[ip][iq] = 0.0;
				v[ip][ip] = 1.0;
			}
			for (ip = 0; ip < n; ip++)
			{
				b[ip] = d[ip] = a[ip][ip];
				z[ip] = 0.0;
			}
			for (i = 1; i <= 50; i++)
			{
				sm = 0.0;
				for (ip = 0; ip < n - 1; ip++)
				{
					for (iq = ip + 1; iq < n; iq++)
						sm += std::abs(a[ip][iq]);
				}
				if (sm == 0.0)
				{
					eigsrt(d, &v);
					return;
				}
				if (i < 4)
					tresh = 0.2 * sm / (n * n);
				else
					tresh = 0.0;
				for (ip = 0; ip < n - 1; ip++)
				{
					for (iq = ip + 1; iq < n; iq++)
					{
						g = 100.0 * std::abs(a[ip][iq]);
						if (i > 4 && g <= EPS * std::abs(d[ip]) && g <= EPS * std::abs(d[iq]))
							a[ip][iq] = 0.0;
						else if (std::abs(a[ip][iq]) > tresh)
						{
							h = d[iq] - d[ip];
							if (g <= EPS * std::abs(h))
								t = (a[ip][iq]) / h;
							else
							{
								theta = 0.5 * h / (a[ip][iq]);
								t = 1.0 / (std::abs(theta) + sqrt(1.0 + theta * theta));
								if (theta < 0.0)
									t = -t;
							}
							c = 1.0 / sqrt(1 + t * t);
							s = t * c;
							tau = s / (1.0 + c);
							h = t * a[ip][iq];
							z[ip] -= h;
							z[iq] += h;
							d[ip] -= h;
							d[iq] += h;
							a[ip][iq] = 0.0;
							for (j = 0; j < ip; j++)
								rot(a, s, tau, j, ip, j, iq);
							for (j = ip + 1; j < iq; j++)
								rot(a, s, tau, ip, j, j, iq);
							for (j = iq + 1; j < n; j++)
								rot(a, s, tau, ip, j, iq, j);
							for (j = 0; j < n; j++)
								rot(v, s, tau, j, ip, j, iq);
							++nrot;
						}
					}
				}
				for (ip = 0; ip < n; ip++)
				{
					b[ip] += z[ip];
					d[ip] = b[ip];
					z[ip] = 0.0;
				}
			}
			throw("Too many iterations in routine jacobi");
		}
		inline void rot(Matrix<Real>& a, const Real s, const Real tau, const int i,
			const int j, const int k, const int l)
		{
			Real g = a[i][j];
			Real h = a[k][l];
			a[i][j] = g - s * (h + g * tau);
			a[k][l] = h + s * (g - h * tau);
		}
	};

	// Computes all eigenvalues and eigenvectors of a real symmetric matrix by reduction to tridiagonal
	// form followed by QL iteration.
	class SymmMatEigenSolver
	{
	public:
		int n;
		Matrix<Real> z;
		Vector<Real> d, e;
		bool yesvecs;

		Vector<Real> getEigenvalues() const { return d; }
		Vector<Real> getEigenvector(int i) const
		{
			Vector<Real> res = z.VectorFromColumn(i);
			return res;
		}

	public:
		// Computes all eigenvalues and eigenvectors of a real symmetric matrix a[0..n-1][0..n-1]
		// by reduction to tridiagonal form followed by QL iteration. On output, d[0..n-1] contains
		// the eigenvalues of a sorted into descending order, while z[0..n-1][0..n-1] is a matrix
		// whose columns contain the corresponding normalized eigenvectors. If yesvecs is input as
		// true (the default), then the eigenvectors are computed. If yesvecs is input as false, only
		// the eigenvalues are computed.
		SymmMatEigenSolver(const MatrixSym<Real>& a, bool yesvec = true) : n(a.RowNum()), d(n), e(n), yesvecs(yesvec)
		{
			z = a.GetAsMatrix();

			tred2();
			tqli();
			sort();
		}

		// Computes all eigenvalues and (optionally) eigenvectors of a real, symmetric, tridiagonal
		// matrix by QL iteration. On input, dd[0..n-1] contains the diagonal elements of the tridiagonal matrix. 
		// The vector ee[0..n-1] inputs the subdiagonal elements of the tridiagonal matrix, with ee[0] arbitrary. 
		// Output is the same as the constructor above.
		SymmMatEigenSolver(Vector<Real>& dd, Vector<Real>& ee, bool yesvec = true) : n((int)dd.size()), d(dd), e(ee), z(n, n), yesvecs(yesvec)
		{
			for (int i = 0; i < n; i++)
				z[i][i] = 1.0;
			tqli();
			sort();
		}
		void sort()
		{
			if (yesvecs)
				eigsrt(d, &z);
			else
				eigsrt(d);
		}

		// Householder reduction of a real symmetric matrix z[0..n-1][0..n-1]. (The input matrix A
		// to Symmeig is stored in z.) On output, z is replaced by the orthogonal matrix Q effecting
		// the transformation. d[0..n-1] contains the diagonal elements of the tridiagonal matrix and
		// e[0..n-1] the off-diagonal elements, with e[0]=0. If yesvecs is false, so that only eigenvalues
		// will subsequently be determined, several statements are omitted, in which case z contains no
		// useful information on output.
		void tred2()
		{
			int l, k, j, i;
			Real scale, hh, h, g, f;
			for (i = n - 1; i > 0; i--)
			{
				l = i - 1;
				h = scale = 0.0;
				if (l > 0)
				{
					for (k = 0; k < i; k++)
						scale += std::abs(z[i][k]);
					if (scale == 0.0)
						e[i] = z[i][l];
					else
					{
						for (k = 0; k < i; k++)
						{
							z[i][k] /= scale;
							h += z[i][k] * z[i][k];
						}
						f = z[i][l];
						g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
						e[i] = scale * g;
						h -= f * g;
						z[i][l] = f - g;
						f = 0.0;
						for (j = 0; j < i; j++)
						{
							if (yesvecs)
								z[j][i] = z[i][j] / h;
							g = 0.0;
							for (k = 0; k < j + 1; k++)
								g += z[j][k] * z[i][k];
							for (k = j + 1; k < i; k++)
								g += z[k][j] * z[i][k];
							e[j] = g / h;
							f += e[j] * z[i][j];
						}
						hh = f / (h + h);
						for (j = 0; j < i; j++)
						{
							f = z[i][j];
							e[j] = g = e[j] - hh * f;
							for (k = 0; k < j + 1; k++)
								z[j][k] -= (f * e[k] + g * z[i][k]);
						}
					}
				}
				else
					e[i] = z[i][l];
				d[i] = h;
			}
			if (yesvecs)
				d[0] = 0.0;
			e[0] = 0.0;
			for (i = 0; i < n; i++)
			{
				if (yesvecs)
				{
					if (d[i] != 0.0)
					{
						for (j = 0; j < i; j++)
						{
							g = 0.0;
							for (k = 0; k < i; k++)
								g += z[i][k] * z[k][j];
							for (k = 0; k < i; k++)
								z[k][j] -= g * z[k][i];
						}
					}
					d[i] = z[i][i];
					z[i][i] = 1.0;
					for (j = 0; j < i; j++)
						z[j][i] = z[i][j] = 0.0;
				}
				else
				{
					d[i] = z[i][i];
				}
			}
		}

		// QL algorithm with implicit shifts to determine the eigenvalues and (optionally) the eigenvectors
		// of a real, symmetric, tridiagonal matrix, or of a real symmetric matrix previously reduced by
		// tred2 (÷11.3). On input, d[0..n-1] contains the diagonal elements of the tridiagonal matrix.
		// On output, it returns the eigenvalues. The vector e[0..n-1] inputs the subdiagonal elements
		// of the tridiagonal matrix, with e[0] arbitrary. On output e is destroyed. If the eigenvectors of
		// a tridiagonal matrix are desired, the matrix z[0..n-1][0..n-1] is input as the identity matrix.
		// If the eigenvectors of a matrix that has been reduced by tred2 are required, then z is input as
		// the matrix output by tred2. In either case, column k of z returns the normalized eigenvector
		// corresponding to d[k]
		void tqli()
		{
			int m, l, iter, i, k;
			Real s, r, p, g, f, dd, c, b;
			const Real EPS = std::numeric_limits<Real>::epsilon();
			for (i = 1; i < n; i++)
				e[i - 1] = e[i];
			e[n - 1] = 0.0;
			for (l = 0; l < n; l++)
			{
				iter = 0;
				do
				{
					for (m = l; m < n - 1; m++)
					{
						dd = std::abs(d[m]) + std::abs(d[m + 1]);
						if (std::abs(e[m]) <= EPS * dd)
							break;
					}
					if (m != l)
					{
						if (iter++ == 30)
							throw("Too many iterations in tqli");
						g = (d[l + 1] - d[l]) / (2.0 * e[l]);
						r = pythag(g, 1.0);
						g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
						s = c = 1.0;
						p = 0.0;
						for (i = m - 1; i >= l; i--)
						{
							f = s * e[i];
							b = c * e[i];
							e[i + 1] = (r = pythag(f, g));
							if (r == 0.0)
							{
								d[i + 1] -= p;
								e[m] = 0.0;
								break;
							}
							s = f / r;
							c = g / r;
							g = d[i + 1] - p;
							r = (d[i] - g) * s + 2.0 * c * b;
							d[i + 1] = g + (p = s * r);
							g = c * r - b;
							if (yesvecs)
							{
								for (k = 0; k < n; k++)
								{
									f = z[k][i + 1];
									z[k][i + 1] = s * z[k][i] + c * f;
									z[k][i] = c * z[k][i] - s * f;
								}
							}
						}
						if (r == 0.0 && i >= l)
							continue;
						d[l] -= p;
						e[l] = g;
						e[m] = 0.0;
					}
				} while (m != l);
			}
		}
		Real pythag(const Real a, const Real b)
		{
			Real absa = std::abs(a), absb = std::abs(b);
			return (absa > absb ? absa * sqrt(1.0 + POW2(absb / absa)) : (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + POW2(absa / absb))));
		}
	};
	//////////////////////////////////////////////////////////////////////////////

	// Computes all eigenvalues and eigenvectors of a real nonsymmetric matrix by reduction to Hessenberg form followed by QR iteration.
	class EigenSolver
	{
	public:
		int n;
		Matrix<Real> a, zz;
		Vector<Complex> wri;
		Vector<Real> scale;
		Vector<int> perm;
		bool yesvecs, hessen;

		bool isRealEigenvalue(int ind) const { return wri[ind].imag() == 0.0; }
		Vector<Complex> getEigenvalues() const { return wri; }

		int getNumReal() const
		{
			int res = 0;
			for (int i = 0; i < n; i++)
				if (wri[i].imag() == 0.0)
					res++;
			return res;
		}
		int getNumComplex() const {
			return n - getNumReal();
		}

		Vector<Real> getRealEigenvalues() const
		{
			Vector<Real> res(getNumReal());
			int cnt = 0;
			for (int i = 0; i < n; i++)
				if (wri[i].imag() == 0.0)
					res[cnt++] = wri[i].real();
			return res;
		}
		Vector<Complex> getComplexEigenvalues() const
		{
			Vector<Complex> res(getNumComplex());
			int cnt = 0;
			for (int i = 0; i < n; i++)
				if (wri[i].imag() != 0.0)
					res[cnt++] = wri[i];
			return res;
		}

		Vector<Complex> getEigenvector(int ind) const
		{
			Vector<Complex> res(n);

			if (isRealEigenvalue(ind))
				for (int i = 0; i < n; i++)
					res[i] = zz(i, ind);
			else
			{
				// count how many real eigenvalues are there before this one
				int cnt = 0;
				for (int i = 0; i < ind; i++)
					if (isRealEigenvalue(i))
						cnt++;
				if (cnt % 2 == 0)
				{
					if (ind % 2 == 0)
					{
						for (int i = 0; i < n; i++)
							res[i] = Complex(zz(i, ind), zz(i, ind + 1));
					}
					else
					{
						for (int i = 0; i < n; i++)
							res[i] = Complex(zz(i, ind - 1), -zz(i, ind));
					}
				}
				else {
					if (ind % 2 == 0)
					{
						for (int i = 0; i < n; i++)
							res[i] = Complex(zz(i, ind - 1), -zz(i, ind));
					}
					else
					{
						for (int i = 0; i < n; i++)
							res[i] = Complex(zz(i, ind), zz(i, ind + 1));
					}
				}
			}

			return res;
		}
		Vector<Real> getRealPartEigenvector(int ind) const
		{
			Vector<Real> res(n);
			for (int i = 0; i < n; i++)
				res[i] = zz(i, ind);

			return res;
		}

		// Computes all eigenvalues and (optionally) eigenvectors of a real nonsymmetric matrix a[0..n-1][0..n-1] by reduction to Hessenberg form followed by QR iteration. 
		// If yesvecs is input as true (the default), then the eigenvectors are computed. Otherwise, only the eigenvalues are computed. 
		// If hessen is input as false (the default), the matrix is first reduced to Hessenberg form. Otherwise it is assumed that the matrix is already in Hessenberg from. 
		// On output, wri[0..n-1] contains the eigenvalues of a sorted into descending order, while zz[0..n-1][0..n-1] is a matrix whose columns contain the corresponding
		// eigenvectors. 
		// For a complex eigenvalue, only the eigenvector corresponding to the eigenvalue with a positive imaginary part is stored, with the real part in zz[0..n-1][i] and the
		// imaginary part in h.zz[0..n-1][i+1]. The eigenvectors are not normalized.

		EigenSolver(const Matrix<Real>& aa, bool yesvec = true, bool hessenb = false) : n(aa.RowNum()), a(aa), zz(n, n), wri(n), scale(n, 1.0), perm(n),
			yesvecs(yesvec), hessen(hessenb)
		{
			balance();
			if (!hessen)
				elmhes();
			if (yesvecs)
			{
				for (int i = 0; i < n; i++)
					zz[i][i] = 1.0;
				if (!hessen)
					eltran();
				hqr2();
				balbak();
				sortvecs();
			}
			else
			{
				hqr();
				sort();
			}
		}
		void balance()
		{
			const Real RADIX = std::numeric_limits<Real>::radix;
			bool done = false;
			Real sqrdx = RADIX * RADIX;
			while (!done)
			{
				done = true;
				for (int i = 0; i < n; i++)
				{
					Real r = 0.0, c = 0.0;
					for (int j = 0; j < n; j++)
						if (j != i)
						{
							c += std::abs(a[j][i]);
							r += std::abs(a[i][j]);
						}
					if (c != 0.0 && r != 0.0)
					{
						Real g = r / RADIX;
						Real f = 1.0;
						Real s = c + r;
						while (c < g)
						{
							f *= RADIX;
							c *= sqrdx;
						}
						g = r * RADIX;
						while (c > g)
						{
							f /= RADIX;
							c /= sqrdx;
						}
						if ((c + r) / f < 0.95 * s)
						{
							done = false;
							g = 1.0 / f;
							scale[i] *= f;
							for (int j = 0; j < n; j++)
								a[i][j] *= g;
							for (int j = 0; j < n; j++)
								a[j][i] *= f;
						}
					}
				}
			}
		}
		void elmhes()
		{
			for (int m = 1; m < n - 1; m++)
			{
				Real x = 0.0;
				int i = m;
				for (int j = m; j < n; j++)
				{
					if (std::abs(a[j][m - 1]) > std::abs(x))
					{
						x = a[j][m - 1];
						i = j;
					}
				}
				perm[m] = i;
				if (i != m)
				{
					for (int j = m - 1; j < n; j++)
						std::swap(a[i][j], a[m][j]);
					for (int j = 0; j < n; j++)
						std::swap(a[j][i], a[j][m]);
				}
				if (x != 0.0)
				{
					for (i = m + 1; i < n; i++)
					{
						Real y = a[i][m - 1];
						if (y != 0.0)
						{
							y /= x;
							a[i][m - 1] = y;
							for (int j = m; j < n; j++)
								a[i][j] -= y * a[m][j];
							for (int j = 0; j < n; j++)
								a[j][m] += y * a[j][i];
						}
					}
				}
			}
		}
		void eltran()
		{
			for (int mp = n - 2; mp > 0; mp--)
			{
				for (int k = mp + 1; k < n; k++)
					zz[k][mp] = a[k][mp - 1];
				int i = perm[mp];
				if (i != mp)
				{
					for (int j = mp; j < n; j++)
					{
						zz[mp][j] = zz[i][j];
						zz[i][j] = 0.0;
					}
					zz[i][mp] = 1.0;
				}
			}
		}

		void hqr()
		{
			int nn, m, l, k, j, its, i, mmin;
			Real z, y, x, w, v, u, t, s, r, q, p, anorm = 0.0;

			const Real EPS = std::numeric_limits<Real>::epsilon();
			for (i = 0; i < n; i++)
				for (j = std::max(i - 1, 0); j < n; j++)
					anorm += std::abs(a[i][j]);
			nn = n - 1;
			t = 0.0;
			while (nn >= 0)
			{
				its = 0;
				do
				{
					for (l = nn; l > 0; l--)
					{
						s = std::abs(a[l - 1][l - 1]) + std::abs(a[l][l]);
						if (s == 0.0)
							s = anorm;
						if (std::abs(a[l][l - 1]) <= EPS * s)
						{
							a[l][l - 1] = 0.0;
							break;
						}
					}
					x = a[nn][nn];
					if (l == nn)
					{
						wri[nn--] = x + t;
					}
					else
					{
						y = a[nn - 1][nn - 1];
						w = a[nn][nn - 1] * a[nn - 1][nn];
						if (l == nn - 1)
						{
							p = 0.5 * (y - x);
							q = p * p + w;
							z = sqrt(std::abs(q));
							x += t;
							if (q >= 0.0)
							{
								z = p + SIGN(z, p);
								wri[nn - 1] = wri[nn] = x + z;
								if (z != 0.0)
									wri[nn] = x - w / z;
							}
							else
							{
								wri[nn] = Complex(x + p, -z);
								wri[nn - 1] = conj(wri[nn]);
							}
							nn -= 2;
						}
						else
						{
							if (its == 30)
								throw("Too many iterations in hqr");
							if (its == 10 || its == 20)
							{
								t += x;
								for (i = 0; i < nn + 1; i++)
									a[i][i] -= x;
								s = std::abs(a[nn][nn - 1]) + std::abs(a[nn - 1][nn - 2]);
								y = x = 0.75 * s;
								w = -0.4375 * s * s;
							}
							++its;
							for (m = nn - 2; m >= l; m--)
							{
								z = a[m][m];
								r = x - z;
								s = y - z;
								p = (r * s - w) / a[m + 1][m] + a[m][m + 1];
								q = a[m + 1][m + 1] - z - r - s;
								r = a[m + 2][m + 1];
								s = std::abs(p) + std::abs(q) + std::abs(r);
								p /= s;
								q /= s;
								r /= s;
								if (m == l)
									break;
								u = std::abs(a[m][m - 1]) * (std::abs(q) + std::abs(r));
								v = std::abs(p) * (std::abs(a[m - 1][m - 1]) + std::abs(z) + std::abs(a[m + 1][m + 1]));
								if (u <= EPS * v)
									break;
							}
							for (i = m; i < nn - 1; i++)
							{
								a[i + 2][i] = 0.0;
								if (i != m)
									a[i + 2][i - 1] = 0.0;
							}
							for (k = m; k < nn; k++)
							{
								if (k != m)
								{
									p = a[k][k - 1];
									q = a[k + 1][k - 1];
									r = 0.0;
									if (k + 1 != nn)
										r = a[k + 2][k - 1];
									if ((x = std::abs(p) + std::abs(q) + std::abs(r)) != 0.0)
									{
										p /= x;
										q /= x;
										r /= x;
									}
								}
								if ((s = SIGN(sqrt(p * p + q * q + r * r), p)) != 0.0)
								{
									if (k == m)
									{
										if (l != m)
											a[k][k - 1] = -a[k][k - 1];
									}
									else
										a[k][k - 1] = -s * x;
									p += s;
									x = p / s;
									y = q / s;
									z = r / s;
									q /= p;
									r /= p;
									for (j = k; j < nn + 1; j++)
									{
										p = a[k][j] + q * a[k + 1][j];
										if (k + 1 != nn)
										{
											p += r * a[k + 2][j];
											a[k + 2][j] -= p * z;
										}
										a[k + 1][j] -= p * y;
										a[k][j] -= p * x;
									}
									mmin = nn < k + 3 ? nn : k + 3;
									for (i = l; i < mmin + 1; i++)
									{
										p = x * a[i][k] + y * a[i][k + 1];
										if (k + 1 != nn)
										{
											p += z * a[i][k + 2];
											a[i][k + 2] -= p * r;
										}
										a[i][k + 1] -= p * q;
										a[i][k] -= p;
									}
								}
							}
						}
					}
				} while (l + 1 < nn);
			}
		}
		void hqr2()
		{
			int nn, m, l, k, j, its, i, mmin, na;
			Real z, y, x, w, v, u, t, s, r, q, p, anorm = 0.0, ra, sa, vr, vi;

			const Real EPS = std::numeric_limits<Real>::epsilon();
			for (i = 0; i < n; i++)
				for (j = std::max(i - 1, 0); j < n; j++)
					anorm += std::abs(a[i][j]);
			nn = n - 1;
			t = 0.0;
			while (nn >= 0)
			{
				its = 0;
				do
				{
					for (l = nn; l > 0; l--)
					{
						s = std::abs(a[l - 1][l - 1]) + std::abs(a[l][l]);
						if (s == 0.0)
							s = anorm;
						if (std::abs(a[l][l - 1]) <= EPS * s)
						{
							a[l][l - 1] = 0.0;
							break;
						}
					}
					x = a[nn][nn];
					if (l == nn)
					{
						wri[nn] = a[nn][nn] = x + t;
						nn--;
					}
					else
					{
						y = a[nn - 1][nn - 1];
						w = a[nn][nn - 1] * a[nn - 1][nn];
						if (l == nn - 1)
						{
							p = 0.5 * (y - x);
							q = p * p + w;
							z = sqrt(std::abs(q));
							x += t;
							a[nn][nn] = x;
							a[nn - 1][nn - 1] = y + t;
							if (q >= 0.0)
							{
								z = p + SIGN(z, p);
								wri[nn - 1] = wri[nn] = x + z;
								if (z != 0.0)
									wri[nn] = x - w / z;
								x = a[nn][nn - 1];
								s = std::abs(x) + std::abs(z);
								p = x / s;
								q = z / s;
								r = sqrt(p * p + q * q);
								p /= r;
								q /= r;
								for (j = nn - 1; j < n; j++)
								{
									z = a[nn - 1][j];
									a[nn - 1][j] = q * z + p * a[nn][j];
									a[nn][j] = q * a[nn][j] - p * z;
								}
								for (i = 0; i <= nn; i++)
								{
									z = a[i][nn - 1];
									a[i][nn - 1] = q * z + p * a[i][nn];
									a[i][nn] = q * a[i][nn] - p * z;
								}
								for (i = 0; i < n; i++)
								{
									z = zz[i][nn - 1];
									zz[i][nn - 1] = q * z + p * zz[i][nn];
									zz[i][nn] = q * zz[i][nn] - p * z;
								}
							}
							else
							{
								wri[nn] = Complex(x + p, -z);
								wri[nn - 1] = conj(wri[nn]);
							}
							nn -= 2;
						}
						else
						{
							if (its == 30)
								throw("Too many iterations in hqr");
							if (its == 10 || its == 20)
							{
								t += x;
								for (i = 0; i < nn + 1; i++)
									a[i][i] -= x;
								s = std::abs(a[nn][nn - 1]) + std::abs(a[nn - 1][nn - 2]);
								y = x = 0.75 * s;
								w = -0.4375 * s * s;
							}
							++its;
							for (m = nn - 2; m >= l; m--)
							{
								z = a[m][m];
								r = x - z;
								s = y - z;
								p = (r * s - w) / a[m + 1][m] + a[m][m + 1];
								q = a[m + 1][m + 1] - z - r - s;
								r = a[m + 2][m + 1];
								s = std::abs(p) + std::abs(q) + std::abs(r);
								p /= s;
								q /= s;
								r /= s;
								if (m == l)
									break;
								u = std::abs(a[m][m - 1]) * (std::abs(q) + std::abs(r));
								v = std::abs(p) * (std::abs(a[m - 1][m - 1]) + std::abs(z) + std::abs(a[m + 1][m + 1]));
								if (u <= EPS * v)
									break;
							}
							for (i = m; i < nn - 1; i++)
							{
								a[i + 2][i] = 0.0;
								if (i != m)
									a[i + 2][i - 1] = 0.0;
							}
							for (k = m; k < nn; k++)
							{
								if (k != m)
								{
									p = a[k][k - 1];
									q = a[k + 1][k - 1];
									r = 0.0;
									if (k + 1 != nn)
										r = a[k + 2][k - 1];
									if ((x = std::abs(p) + std::abs(q) + std::abs(r)) != 0.0)
									{
										p /= x;
										q /= x;
										r /= x;
									}
								}
								if ((s = SIGN(sqrt(p * p + q * q + r * r), p)) != 0.0)
								{
									if (k == m)
									{
										if (l != m)
											a[k][k - 1] = -a[k][k - 1];
									}
									else
										a[k][k - 1] = -s * x;
									p += s;
									x = p / s;
									y = q / s;
									z = r / s;
									q /= p;
									r /= p;
									for (j = k; j < n; j++)
									{
										p = a[k][j] + q * a[k + 1][j];
										if (k + 1 != nn)
										{
											p += r * a[k + 2][j];
											a[k + 2][j] -= p * z;
										}
										a[k + 1][j] -= p * y;
										a[k][j] -= p * x;
									}
									mmin = nn < k + 3 ? nn : k + 3;
									for (i = 0; i < mmin + 1; i++)
									{
										p = x * a[i][k] + y * a[i][k + 1];
										if (k + 1 != nn)
										{
											p += z * a[i][k + 2];
											a[i][k + 2] -= p * r;
										}
										a[i][k + 1] -= p * q;
										a[i][k] -= p;
									}
									for (i = 0; i < n; i++)
									{
										p = x * zz[i][k] + y * zz[i][k + 1];
										if (k + 1 != nn)
										{
											p += z * zz[i][k + 2];
											zz[i][k + 2] -= p * r;
										}
										zz[i][k + 1] -= p * q;
										zz[i][k] -= p;
									}
								}
							}
						}
					}
				} while (l + 1 < nn);
			}
			if (anorm != 0.0)
			{
				for (nn = n - 1; nn >= 0; nn--)
				{
					p = real(wri[nn]);
					q = imag(wri[nn]);
					na = nn - 1;
					if (q == 0.0)
					{
						m = nn;
						a[nn][nn] = 1.0;
						for (i = nn - 1; i >= 0; i--)
						{
							w = a[i][i] - p;
							r = 0.0;
							for (j = m; j <= nn; j++)
								r += a[i][j] * a[j][nn];
							if (imag(wri[i]) < 0.0)
							{
								z = w;
								s = r;
							}
							else
							{
								m = i;

								if (imag(wri[i]) == 0.0)
								{
									t = w;
									if (t == 0.0)
										t = EPS * anorm;
									a[i][nn] = -r / t;
								}
								else
								{
									x = a[i][i + 1];
									y = a[i + 1][i];
									q = POW2(real(wri[i]) - p) + POW2(imag(wri[i]));
									t = (x * s - z * r) / q;
									a[i][nn] = t;
									if (std::abs(x) > std::abs(z))
										a[i + 1][nn] = (-r - w * t) / x;
									else
										a[i + 1][nn] = (-s - y * t) / z;
								}
								t = std::abs(a[i][nn]);
								if (EPS * t * t > 1)
									for (j = i; j <= nn; j++)
										a[j][nn] /= t;
							}
						}
					}
					else if (q < 0.0)
					{
						m = na;
						if (std::abs(a[nn][na]) > std::abs(a[na][nn]))
						{
							a[na][na] = q / a[nn][na];
							a[na][nn] = -(a[nn][nn] - p) / a[nn][na];
						}
						else
						{
							Complex temp = Complex(0.0, -a[na][nn]) / Complex(a[na][na] - p, q);
							a[na][na] = real(temp);
							a[na][nn] = imag(temp);
						}
						a[nn][na] = 0.0;
						a[nn][nn] = 1.0;
						for (i = nn - 2; i >= 0; i--)
						{
							w = a[i][i] - p;
							ra = sa = 0.0;
							for (j = m; j <= nn; j++)
							{
								ra += a[i][j] * a[j][na];
								sa += a[i][j] * a[j][nn];
							}
							if (imag(wri[i]) < 0.0)
							{
								z = w;
								r = ra;
								s = sa;
							}
							else
							{
								m = i;
								if (imag(wri[i]) == 0.0)
								{
									Complex temp = Complex(-ra, -sa) / Complex(w, q);
									a[i][na] = real(temp);
									a[i][nn] = imag(temp);
								}
								else
								{
									x = a[i][i + 1];
									y = a[i + 1][i];
									vr = POW2(real(wri[i]) - p) + POW2(imag(wri[i])) - q * q;
									vi = 2.0 * q * (real(wri[i]) - p);
									if (vr == 0.0 && vi == 0.0)
										vr = EPS * anorm * (std::abs(w) + std::abs(q) + std::abs(x) + std::abs(y) + std::abs(z));
									Complex temp = Complex(x * r - z * ra + q * sa, x * s - z * sa - q * ra) /
										Complex(vr, vi);
									a[i][na] = real(temp);
									a[i][nn] = imag(temp);
									if (std::abs(x) > std::abs(z) + std::abs(q))
									{
										a[i + 1][na] = (-ra - w * a[i][na] + q * a[i][nn]) / x;
										a[i + 1][nn] = (-sa - w * a[i][nn] - q * a[i][na]) / x;
									}
									else
									{
										Complex temp = Complex(-r - y * a[i][na], -s - y * a[i][nn]) /
											Complex(z, q);
										a[i + 1][na] = real(temp);
										a[i + 1][nn] = imag(temp);
									}
								}
							}
							t = std::max(std::abs(a[i][na]), std::abs(a[i][nn]));
							if (EPS * t * t > 1)
								for (j = i; j <= nn; j++)
								{
									a[j][na] /= t;
									a[j][nn] /= t;
								}
						}
					}
				}
				for (j = n - 1; j >= 0; j--)
					for (i = 0; i < n; i++)
					{
						z = 0.0;
						for (k = 0; k <= j; k++)
							z += zz[i][k] * a[k][j];
						zz[i][j] = z;
					}
			}
		}
		void balbak()
		{
			for (int i = 0; i < n; i++)
				for (int j = 0; j < n; j++)
					zz[i][j] *= scale[i];
		}

		void sort()
		{
			int i;
			for (int j = 1; j < n; j++)
			{
				Complex x = wri[j];
				for (i = j - 1; i >= 0; i--)
				{
					if (real(wri[i]) >= real(x))
						break;
					wri[i + 1] = wri[i];
				}
				wri[i + 1] = x;
			}
		}
		void sortvecs()
		{
			int i;
			Vector<Real> temp(n);
			for (int j = 1; j < n; j++)
			{
				Complex x = wri[j];
				for (int k = 0; k < n; k++)
					temp[k] = zz[k][j];
				for (i = j - 1; i >= 0; i--)
				{
					if (real(wri[i]) >= real(x))
						break;
					wri[i + 1] = wri[i];
					for (int k = 0; k < n; k++)
						zz[k][i + 1] = zz[k][i];
				}
				wri[i + 1] = x;
				for (int k = 0; k < n; k++)
					zz[k][i + 1] = temp[k];
			}
		}
	};

}
///////////////////////////   ./include/algorithms/ODESystemSteppers.h   ///////////////////////////



namespace MML
{
	struct StepperBase {
		IODESystem& _sys;
		int    n, neqn;
		bool   dense;

		Real& x;
		Real xold;
		Vector<Real>& y, & dydx;

		Real atol, rtol;
		Real EPS;

		Real hdid;
		Real hnext;
		Vector<Real> yout, yerr;

		StepperBase(IODESystem& sys, Vector<Real>& yy, Vector<Real>& dydxx, Real& xx, const Real atoll, const Real rtoll, bool dens)
			: _sys(sys), x(xx), y(yy), dydx(dydxx), atol(atoll),
			rtol(rtoll), dense(dens), n(sys.getDim()), neqn(n), yout(n), yerr(n) {}

		virtual Real dense_out(const int i, const Real x, const Real h) = 0;
	};

	// Dormand-Prince fifth-order Runge-Kutta step with monitoring of local truncation error to ensure
	// accuracy and adjust stepsize.
	struct StepperDopr5 : StepperBase {
		Vector<Real> k2, k3, k4, k5, k6;
		Vector<Real> rcont1, rcont2, rcont3, rcont4, rcont5;
		Vector<Real> dydxnew;

		// Input to the constructor are the dependent variable y[0..n-1] and its derivative dydx[0..n-1]
		// at the starting value of the independent variable x. Also input are the absolute and relative
		// tolerances, atol and rtol, and the boolean dense, which is true if dense output is required        
		StepperDopr5(IODESystem& sys, Vector<Real>& yy, Vector<Real>& dydx, Real& xx, const Real atoll, const Real rtoll, bool dens)
			: StepperBase(sys, yy, dydx, xx, atoll, rtoll, dens),
			k2(n), k3(n), k4(n), k5(n), k6(n),
			rcont1(n), rcont2(n), rcont3(n), rcont4(n), rcont5(n),
			dydxnew(n)
		{
			EPS = std::numeric_limits<Real>::epsilon();
		}

		// Attempts a step with stepsize htry. On output, y and x are replaced by their new values, hdid
		// is the stepsize that was actually accomplished, and hnext is the estimated next stepsize        
		void step(const Real htry) {
			Real h = htry;
			for (;;) {
				dy(h);
				Real err = error();
				if (con.success(err, h)) break;
				if (std::abs(h) <= std::abs(x) * EPS)
					throw("stepsize underflow in StepperDopr5");
			}
			if (dense)
				prepare_dense(h);
			dydx = dydxnew;
			y = yout;
			xold = x;
			x += (hdid = h);
			hnext = con.hnext;
		}

		// Given values for n variables y[0..n-1] and their derivatives dydx[0..n-1] known at x, use the
		// fifth-order Dormand-Prince Runge-Kutta method to advance the solution over an interval h and
		// store the incremented variables in yout[0..n-1]. Also store an estimate of the local truncation
		// error in yerr using the embedded fourth-order method        
		void dy(const Real h) {
			static const Real c2 = 0.2, c3 = 0.3, c4 = 0.8, c5 = 8.0 / 9.0, a21 = 0.2, a31 = 3.0 / 40.0,
				a32 = 9.0 / 40.0, a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0, a51 = 19372.0 / 6561.0,
				a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0, a61 = 9017.0 / 3168.0,
				a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0, a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0,
				a71 = 35.0 / 384.0, a73 = 500.0 / 1113.0, a74 = 125.0 / 192.0, a75 = -2187.0 / 6784.0,
				a76 = 11.0 / 84.0, e1 = 71.0 / 57600.0, e3 = -71.0 / 16695.0, e4 = 71.0 / 1920.0,
				e5 = -17253.0 / 339200.0, e6 = 22.0 / 525.0, e7 = -1.0 / 40.0;
			Vector<Real> ytemp(n);
			int i;
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * a21 * dydx[i];
			_sys.derivs(x + c2 * h, ytemp, k2);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (a31 * dydx[i] + a32 * k2[i]);
			_sys.derivs(x + c3 * h, ytemp, k3);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (a41 * dydx[i] + a42 * k2[i] + a43 * k3[i]);
			_sys.derivs(x + c4 * h, ytemp, k4);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (a51 * dydx[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
			_sys.derivs(x + c5 * h, ytemp, k5);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (a61 * dydx[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
			Real xph = x + h;
			_sys.derivs(xph, ytemp, k6);
			for (i = 0; i < n; i++)
				yout[i] = y[i] + h * (a71 * dydx[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
			_sys.derivs(xph, yout, dydxnew);
			for (i = 0; i < n; i++) {
				yerr[i] = h * (e1 * dydx[i] + e3 * k3[i] + e4 * k4[i] + e5 * k5[i] + e6 * k6[i] + e7 * dydxnew[i]);
			}
		}

		// Store coefficients of interpolating polynomial for dense output in rcont1...rcont5        
		void prepare_dense(const Real h) {
			Vector<Real> ytemp(n);
			static const Real d1 = -12715105075.0 / 11282082432.0,
				d3 = 87487479700.0 / 32700410799.0, d4 = -10690763975.0 / 1880347072.0,
				d5 = 701980252875.0 / 199316789632.0, d6 = -1453857185.0 / 822651844.0,
				d7 = 69997945.0 / 29380423.0;
			for (int i = 0; i < n; i++) {
				rcont1[i] = y[i];
				Real ydiff = yout[i] - y[i];
				rcont2[i] = ydiff;
				Real bspl = h * dydx[i] - ydiff;
				rcont3[i] = bspl;
				rcont4[i] = ydiff - h * dydxnew[i] - bspl;
				rcont5[i] = h * (d1 * dydx[i] + d3 * k3[i] + d4 * k4[i] + d5 * k5[i] + d6 * k6[i] +
					d7 * dydxnew[i]);
			}
		}

		// Evaluate interpolating polynomial for y[i] at location x, where xold <= x <= xold + h.
		Real dense_out(const int i, const Real x, const Real h) {
			Real s = (x - xold) / h;
			Real s1 = 1.0 - s;
			return rcont1[i] + s * (rcont2[i] + s1 * (rcont3[i] + s * (rcont4[i] + s1 * rcont5[i])));
		}

		// Use yerr to compute norm of scaled error estimate. A value less than one means the step was successful.        
		Real error() {
			Real err = 0.0, sk;
			for (int i = 0; i < n; i++) {
				sk = atol + rtol * std::max(std::abs(y[i]), std::abs(yout[i]));
				err += POW2(yerr[i] / sk);
			}
			return sqrt(err / n);
		}

		// Finally, the controller tests whether err <= 1 and adjusts the stepsize. The
		// default setting is beta D 0 (no PI control). Set beta to 0.04 or 0.08 to turn on PI control        
		struct Controller {
			Real hnext, errold;
			bool reject;

			Controller() : reject(false), errold(1.0e-4) {}

			// Returns true if err <= 1, false otherwise. If step was successful, sets hnext to the estimated
			// optimal stepsize for the next step. If the step failed, reduces h appropriately for another try.
			bool success(const Real err, Real& h) {
				static const Real beta = 0.0, alpha = 0.2 - beta * 0.75, safe = 0.9, minscale = 0.2, maxscale = 10.0;
				Real scale;
				if (err <= 1.0) {
					if (err == 0.0)
						scale = maxscale;
					else {
						scale = safe * pow(err, -alpha) * pow(errold, beta);
						if (scale < minscale) scale = minscale;
						if (scale > maxscale) scale = maxscale;
					}
					if (reject)
						hnext = h * std::min<Real>(scale, 1.0);
					else
						hnext = h * scale;
					errold = std::max<Real>(err, 1.0e-4);
					reject = false;
					return true;
				}
				else {
					scale = std::max(safe * pow(err, -alpha), minscale);
					h *= scale;
					reject = true;
					return false;
				}
			}
		};
		Controller con;
	};

	struct Dopr853_constants {
		const Real c2 = 0.526001519587677318785587544488e-01;
		const Real c3 = 0.789002279381515978178381316732e-01;
		const Real c4 = 0.118350341907227396726757197510e+00;
		const Real c5 = 0.281649658092772603273242802490e+00;
		const Real c6 = 0.333333333333333333333333333333e+00;
		const Real c7 = 0.25e+00;
		const Real c8 = 0.307692307692307692307692307692e+00;
		const Real c9 = 0.651282051282051282051282051282e+00;
		const Real c10 = 0.6e+00;
		const Real c11 = 0.857142857142857142857142857142e+00;
		const Real c14 = 0.1e+00;
		const Real c15 = 0.2e+00;
		const Real c16 = 0.777777777777777777777777777778e+00;

		const Real b1 = 5.42937341165687622380535766363e-2;
		const Real b6 = 4.45031289275240888144113950566e0;
		const Real b7 = 1.89151789931450038304281599044e0;
		const Real b8 = -5.8012039600105847814672114227e0;
		const Real b9 = 3.1116436695781989440891606237e-1;
		const Real b10 = -1.52160949662516078556178806805e-1;
		const Real b11 = 2.01365400804030348374776537501e-1;
		const Real b12 = 4.47106157277725905176885569043e-2;

		const Real bhh1 = 0.244094488188976377952755905512e+00;
		const Real bhh2 = 0.733846688281611857341361741547e+00;
		const Real bhh3 = 0.220588235294117647058823529412e-01;

		const Real er1 = 0.1312004499419488073250102996e-01;
		const Real er6 = -0.1225156446376204440720569753e+01;
		const Real er7 = -0.4957589496572501915214079952e+00;
		const Real er8 = 0.1664377182454986536961530415e+01;
		const Real er9 = -0.3503288487499736816886487290e+00;
		const Real er10 = 0.3341791187130174790297318841e+00;
		const Real er11 = 0.8192320648511571246570742613e-01;
		const Real er12 = -0.2235530786388629525884427845e-01;

		const Real a21 = 5.26001519587677318785587544488e-2;
		const Real a31 = 1.97250569845378994544595329183e-2;
		const Real a32 = 5.91751709536136983633785987549e-2;
		const Real a41 = 2.95875854768068491816892993775e-2;
		const Real a43 = 8.87627564304205475450678981324e-2;
		const Real a51 = 2.41365134159266685502369798665e-1;
		const Real a53 = -8.84549479328286085344864962717e-1;
		const Real a54 = 9.24834003261792003115737966543e-1;
		const Real a61 = 3.7037037037037037037037037037e-2;
		const Real a64 = 1.70828608729473871279604482173e-1;
		const Real a65 = 1.25467687566822425016691814123e-1;
		const Real a71 = 3.7109375e-2;
		const Real a74 = 1.70252211019544039314978060272e-1;
		const Real a75 = 6.02165389804559606850219397283e-2;
		const Real a76 = -1.7578125e-2;

		const Real a81 = 3.70920001185047927108779319836e-2;
		const Real a84 = 1.70383925712239993810214054705e-1;
		const Real a85 = 1.07262030446373284651809199168e-1;
		const Real a86 = -1.53194377486244017527936158236e-2;
		const Real a87 = 8.27378916381402288758473766002e-3;
		const Real a91 = 6.24110958716075717114429577812e-1;
		const Real a94 = -3.36089262944694129406857109825e0;
		const Real a95 = -8.68219346841726006818189891453e-1;
		const Real a96 = 2.75920996994467083049415600797e1;
		const Real a97 = 2.01540675504778934086186788979e1;
		const Real a98 = -4.34898841810699588477366255144e1;
		const Real a101 = 4.77662536438264365890433908527e-1;
		const Real a104 = -2.48811461997166764192642586468e0;
		const Real a105 = -5.90290826836842996371446475743e-1;
		const Real a106 = 2.12300514481811942347288949897e1;
		const Real a107 = 1.52792336328824235832596922938e1;
		const Real a108 = -3.32882109689848629194453265587e1;
		const Real a109 = -2.03312017085086261358222928593e-2;

		const Real a111 = -9.3714243008598732571704021658e-1;
		const Real a114 = 5.18637242884406370830023853209e0;
		const Real a115 = 1.09143734899672957818500254654e0;
		const Real a116 = -8.14978701074692612513997267357e0;
		const Real a117 = -1.85200656599969598641566180701e1;
		const Real a118 = 2.27394870993505042818970056734e1;
		const Real a119 = 2.49360555267965238987089396762e0;
		const Real a1110 = -3.0467644718982195003823669022e0;
		const Real a121 = 2.27331014751653820792359768449e0;
		const Real a124 = -1.05344954667372501984066689879e1;
		const Real a125 = -2.00087205822486249909675718444e0;
		const Real a126 = -1.79589318631187989172765950534e1;
		const Real a127 = 2.79488845294199600508499808837e1;
		const Real a128 = -2.85899827713502369474065508674e0;
		const Real a129 = -8.87285693353062954433549289258e0;
		const Real a1210 = 1.23605671757943030647266201528e1;
		const Real a1211 = 6.43392746015763530355970484046e-1;

		const Real a141 = 5.61675022830479523392909219681e-2;
		const Real a147 = 2.53500210216624811088794765333e-1;
		const Real a148 = -2.46239037470802489917441475441e-1;
		const Real a149 = -1.24191423263816360469010140626e-1;
		const Real a1410 = 1.5329179827876569731206322685e-1;
		const Real a1411 = 8.20105229563468988491666602057e-3;
		const Real a1412 = 7.56789766054569976138603589584e-3;
		const Real a1413 = -8.298e-3;

		const Real a151 = 3.18346481635021405060768473261e-2;
		const Real a156 = 2.83009096723667755288322961402e-2;
		const Real a157 = 5.35419883074385676223797384372e-2;
		const Real a158 = -5.49237485713909884646569340306e-2;
		const Real a1511 = -1.08347328697249322858509316994e-4;
		const Real a1512 = 3.82571090835658412954920192323e-4;
		const Real a1513 = -3.40465008687404560802977114492e-4;
		const Real a1514 = 1.41312443674632500278074618366e-1;
		const Real a161 = -4.28896301583791923408573538692e-1;
		const Real a166 = -4.69762141536116384314449447206e0;
		const Real a167 = 7.68342119606259904184240953878e0;
		const Real a168 = 4.06898981839711007970213554331e0;
		const Real a169 = 3.56727187455281109270669543021e-1;
		const Real a1613 = -1.39902416515901462129418009734e-3;
		const Real a1614 = 2.9475147891527723389556272149e0;
		const Real a1615 = -9.15095847217987001081870187138e0;

		const Real d41 = -0.84289382761090128651353491142e+01;
		const Real d46 = 0.56671495351937776962531783590e+00;
		const Real d47 = -0.30689499459498916912797304727e+01;
		const Real d48 = 0.23846676565120698287728149680e+01;
		const Real d49 = 0.21170345824450282767155149946e+01;
		const Real d410 = -0.87139158377797299206789907490e+00;
		const Real d411 = 0.22404374302607882758541771650e+01;
		const Real d412 = 0.63157877876946881815570249290e+00;
		const Real d413 = -0.88990336451333310820698117400e-01;
		const Real d414 = 0.18148505520854727256656404962e+02;
		const Real d415 = -0.91946323924783554000451984436e+01;
		const Real d416 = -0.44360363875948939664310572000e+01;

		const Real d51 = 0.10427508642579134603413151009e+02;
		const Real d56 = 0.24228349177525818288430175319e+03;
		const Real d57 = 0.16520045171727028198505394887e+03;
		const Real d58 = -0.37454675472269020279518312152e+03;
		const Real d59 = -0.22113666853125306036270938578e+02;
		const Real d510 = 0.77334326684722638389603898808e+01;
		const Real d511 = -0.30674084731089398182061213626e+02;
		const Real d512 = -0.93321305264302278729567221706e+01;
		const Real d513 = 0.15697238121770843886131091075e+02;
		const Real d514 = -0.31139403219565177677282850411e+02;
		const Real d515 = -0.93529243588444783865713862664e+01;
		const Real d516 = 0.35816841486394083752465898540e+02;

		const Real d61 = 0.19985053242002433820987653617e+02;
		const Real d66 = -0.38703730874935176555105901742e+03;
		const Real d67 = -0.18917813819516756882830838328e+03;
		const Real d68 = 0.52780815920542364900561016686e+03;
		const Real d69 = -0.11573902539959630126141871134e+02;
		const Real d610 = 0.68812326946963000169666922661e+01;
		const Real d611 = -0.10006050966910838403183860980e+01;
		const Real d612 = 0.77771377980534432092869265740e+00;
		const Real d613 = -0.27782057523535084065932004339e+01;
		const Real d614 = -0.60196695231264120758267380846e+02;
		const Real d615 = 0.84320405506677161018159903784e+02;
		const Real d616 = 0.11992291136182789328035130030e+02;

		const Real d71 = -0.25693933462703749003312586129e+02;
		const Real d76 = -0.15418974869023643374053993627e+03;
		const Real d77 = -0.23152937917604549567536039109e+03;
		const Real d78 = 0.35763911791061412378285349910e+03;
		const Real d79 = 0.93405324183624310003907691704e+02;
		const Real d710 = -0.37458323136451633156875139351e+02;
		const Real d711 = 0.10409964950896230045147246184e+03;
		const Real d712 = 0.29840293426660503123344363579e+02;
		const Real d713 = -0.43533456590011143754432175058e+02;
		const Real d714 = 0.96324553959188282948394950600e+02;
		const Real d715 = -0.39177261675615439165231486172e+02;
		const Real d716 = -0.14972683625798562581422125276e+03;
	};

	struct StepperDopr853 : StepperBase, Dopr853_constants {
		Vector<Real> yerr2;
		Vector<Real> k2, k3, k4, k5, k6, k7, k8, k9, k10;
		Vector<Real> rcont1, rcont2, rcont3, rcont4, rcont5, rcont6, rcont7, rcont8;

		StepperDopr853(IODESystem& sys, Vector<Real>& yy, Vector<Real>& dydxx, Real& xx,
			const Real atoll, const Real rtoll, bool dens)
			: StepperBase(sys, yy, dydxx, xx, atoll, rtoll, dens),
			yerr2(n), k2(n), k3(n), k4(n),
			k5(n), k6(n), k7(n), k8(n), k9(n), k10(n), rcont1(n), rcont2(n), rcont3(n),
			rcont4(n), rcont5(n), rcont6(n), rcont7(n), rcont8(n)
		{
			EPS = std::numeric_limits<Real>::epsilon();
		}

		void step(const Real htry) {
			Vector<Real> dydxnew(n);
			Real h = htry;
			for (;;) {
				dy(h);
				Real err = error(h);
				if (con.success(err, h)) break;
				if (std::abs(h) <= std::abs(x) * EPS)
					throw("stepsize underflow in StepperDopr853");
			}
			_sys.derivs(x + h, yout, dydxnew);
			if (dense)
				prepare_dense(h, dydxnew);
			dydx = dydxnew;
			y = yout;
			xold = x;
			x += (hdid = h);
			hnext = con.hnext;
		}

		void dy(const Real h) {
			Vector<Real> ytemp(n);
			int i;
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * a21 * dydx[i];
			_sys.derivs(x + c2 * h, ytemp, k2);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (a31 * dydx[i] + a32 * k2[i]);
			_sys.derivs(x + c3 * h, ytemp, k3);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (a41 * dydx[i] + a43 * k3[i]);
			_sys.derivs(x + c4 * h, ytemp, k4);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (a51 * dydx[i] + a53 * k3[i] + a54 * k4[i]);
			_sys.derivs(x + c5 * h, ytemp, k5);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (a61 * dydx[i] + a64 * k4[i] + a65 * k5[i]);
			_sys.derivs(x + c6 * h, ytemp, k6);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (a71 * dydx[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
			_sys.derivs(x + c7 * h, ytemp, k7);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (a81 * dydx[i] + a84 * k4[i] + a85 * k5[i] + a86 * k6[i] + a87 * k7[i]);
			_sys.derivs(x + c8 * h, ytemp, k8);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (a91 * dydx[i] + a94 * k4[i] + a95 * k5[i] + a96 * k6[i] + a97 * k7[i] +
					a98 * k8[i]);
			_sys.derivs(x + c9 * h, ytemp, k9);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (a101 * dydx[i] + a104 * k4[i] + a105 * k5[i] + a106 * k6[i] +
					a107 * k7[i] + a108 * k8[i] + a109 * k9[i]);
			_sys.derivs(x + c10 * h, ytemp, k10);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (a111 * dydx[i] + a114 * k4[i] + a115 * k5[i] + a116 * k6[i] +
					a117 * k7[i] + a118 * k8[i] + a119 * k9[i] + a1110 * k10[i]);
			_sys.derivs(x + c11 * h, ytemp, k2);
			Real xph = x + h;
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (a121 * dydx[i] + a124 * k4[i] + a125 * k5[i] + a126 * k6[i] +
					a127 * k7[i] + a128 * k8[i] + a129 * k9[i] + a1210 * k10[i] + a1211 * k2[i]);
			_sys.derivs(xph, ytemp, k3);
			for (i = 0; i < n; i++) {
				k4[i] = b1 * dydx[i] + b6 * k6[i] + b7 * k7[i] + b8 * k8[i] + b9 * k9[i] + b10 * k10[i] +
					b11 * k2[i] + b12 * k3[i];
				yout[i] = y[i] + h * k4[i];
			}
			for (i = 0; i < n; i++) {
				yerr[i] = k4[i] - bhh1 * dydx[i] - bhh2 * k9[i] - bhh3 * k3[i];
				yerr2[i] = er1 * dydx[i] + er6 * k6[i] + er7 * k7[i] + er8 * k8[i] + er9 * k9[i] +
					er10 * k10[i] + er11 * k2[i] + er12 * k3[i];
			}
		}
		void prepare_dense(const Real h, const Vector<Real>& dydxnew)
		{
			int i;
			Real ydiff, bspl;
			Vector<Real> ytemp(n);
			for (i = 0; i < n; i++) {
				rcont1[i] = y[i];
				ydiff = yout[i] - y[i];
				rcont2[i] = ydiff;
				bspl = h * dydx[i] - ydiff;
				rcont3[i] = bspl;
				rcont4[i] = ydiff - h * dydxnew[i] - bspl;
				rcont5[i] = d41 * dydx[i] + d46 * k6[i] + d47 * k7[i] + d48 * k8[i] +
					d49 * k9[i] + d410 * k10[i] + d411 * k2[i] + d412 * k3[i];
				rcont6[i] = d51 * dydx[i] + d56 * k6[i] + d57 * k7[i] + d58 * k8[i] +
					d59 * k9[i] + d510 * k10[i] + d511 * k2[i] + d512 * k3[i];
				rcont7[i] = d61 * dydx[i] + d66 * k6[i] + d67 * k7[i] + d68 * k8[i] +
					d69 * k9[i] + d610 * k10[i] + d611 * k2[i] + d612 * k3[i];
				rcont8[i] = d71 * dydx[i] + d76 * k6[i] + d77 * k7[i] + d78 * k8[i] +
					d79 * k9[i] + d710 * k10[i] + d711 * k2[i] + d712 * k3[i];
			}
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (a141 * dydx[i] + a147 * k7[i] + a148 * k8[i] + a149 * k9[i] +
					a1410 * k10[i] + a1411 * k2[i] + a1412 * k3[i] + a1413 * dydxnew[i]);
			_sys.derivs(x + c14 * h, ytemp, k10);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (a151 * dydx[i] + a156 * k6[i] + a157 * k7[i] + a158 * k8[i] +
					a1511 * k2[i] + a1512 * k3[i] + a1513 * dydxnew[i] + a1514 * k10[i]);
			_sys.derivs(x + c15 * h, ytemp, k2);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + h * (a161 * dydx[i] + a166 * k6[i] + a167 * k7[i] + a168 * k8[i] +
					a169 * k9[i] + a1613 * dydxnew[i] + a1614 * k10[i] + a1615 * k2[i]);
			_sys.derivs(x + c16 * h, ytemp, k3);
			for (i = 0; i < n; i++)
			{
				rcont5[i] = h * (rcont5[i] + d413 * dydxnew[i] + d414 * k10[i] + d415 * k2[i] + d416 * k3[i]);
				rcont6[i] = h * (rcont6[i] + d513 * dydxnew[i] + d514 * k10[i] + d515 * k2[i] + d516 * k3[i]);
				rcont7[i] = h * (rcont7[i] + d613 * dydxnew[i] + d614 * k10[i] + d615 * k2[i] + d616 * k3[i]);
				rcont8[i] = h * (rcont8[i] + d713 * dydxnew[i] + d714 * k10[i] + d715 * k2[i] + d716 * k3[i]);
			}
		}
		Real dense_out(const int i, const Real x, const Real h) {
			Real s = (x - xold) / h;
			Real s1 = 1.0 - s;
			return rcont1[i] + s * (rcont2[i] + s1 * (rcont3[i] + s * (rcont4[i] + s1 * (rcont5[i] +
				s * (rcont6[i] + s1 * (rcont7[i] + s * rcont8[i]))))));
		}
		Real error(const Real h) {
			Real err = 0.0, err2 = 0.0, sk, deno;
			for (int i = 0; i < n; i++) {
				sk = atol + rtol * std::max(std::abs(y[i]), std::abs(yout[i]));
				err2 += POW2(yerr[i] / sk);
				err += POW2(yerr2[i] / sk);
			}
			deno = err + 0.01 * err2;
			if (deno <= 0.0)
				deno = 1.0;
			return std::abs(h) * err * sqrt(1.0 / (n * deno));
		}
		struct Controller {
			Real hnext, errold;
			bool reject;
			Controller() : reject(false), errold(1.0e-4) {}
			bool success(const Real err, Real& h) {
				static const Real beta = 0.0, alpha = 1.0 / 8.0 - beta * 0.2, safe = 0.9, minscale = 0.333,
					maxscale = 6.0;
				Real scale;
				if (err <= 1.0) {
					if (err == 0.0)
						scale = maxscale;
					else {
						scale = safe * pow(err, -alpha) * pow(errold, beta);
						if (scale < minscale) scale = minscale;
						if (scale > maxscale) scale = maxscale;
					}
					if (reject)
						hnext = h * std::min<Real>(scale, 1.0);
					else
						hnext = h * scale;
					errold = std::max<Real>(err, 1.0e-4);
					reject = false;
					return true;
				}
				else {
					scale = std::max<Real>(safe * pow(err, -alpha), minscale);
					h *= scale;
					reject = true;
					return false;
				}
			}
		};
		Controller con;
	};

	// Bulirsch-Stoer step with monitoring of local truncation error to ensure accuracy and adjust stepsize
	struct StepperBS : StepperBase {
		static const int KMAXX = 8, IMAXX = KMAXX + 1;     // KMAXX is the maximum number of rows used in the extrapolation
		int k_targ;                 // Optimal row number for convergence.
		Vector<int> nseq;           // Stepsize sequence
		Vector<int> cost;
		Matrix<Real> table;         // Extrapolation tableau
		Vector<Real> dydxnew;
		int mu;
		Matrix<Real> coeff;         // Coefficients used in extrapolation tableau
		Vector<Real> errfac;        // Used to compute dense interpolation error.
		Matrix<Real> ysave;         // ysave and fsave store values and derivatives to be used for dense output
		Matrix<Real> fsave;
		Vector<int> ipoint;         // Keeps track of where values are stored in fsave.
		Vector<Real> dens;

		StepperBS(IODESystem& sys, Vector<Real>& yy, Vector<Real>& dydxx, Real& xx, const Real atoll, const Real rtoll, bool dens)
			: StepperBase(sys, yy, dydxx, xx, atoll, rtoll, dens),
			nseq(IMAXX), cost(IMAXX),
			table(KMAXX, n), dydxnew(n), coeff(IMAXX, IMAXX), errfac(2 * IMAXX + 2), ysave(IMAXX, n),
			fsave(IMAXX* (2 * IMAXX + 1), n), ipoint(IMAXX + 1), dens((2 * IMAXX + 5)* n)
		{
			EPS = std::numeric_limits<Real>::epsilon();
			if (dense)
				for (int i = 0; i < IMAXX; i++)
					nseq[i] = 4 * i + 2;
			else
				for (int i = 0; i < IMAXX; i++)
					nseq[i] = 2 * (i + 1);
			cost[0] = nseq[0] + 1;
			for (int k = 0; k < KMAXX; k++) cost[k + 1] = cost[k] + nseq[k + 1];
			hnext = -1.0e99;
			Real logfact = -log10(std::max<Real>(1.0e-12, rtol)) * 0.6 + 0.5;
			k_targ = std::max<int>(1, std::min<int>(KMAXX - 1, int(logfact)));
			for (int k = 0; k < IMAXX; k++) {
				for (int l = 0; l < k; l++) {
					Real ratio = Real(nseq[k]) / nseq[l];
					coeff[k][l] = 1.0 / (ratio * ratio - 1.0);
				}
			}
			for (int i = 0; i < 2 * IMAXX + 1; i++) {
				int ip5 = i + 5;
				errfac[i] = 1.0 / (ip5 * ip5);
				Real e = 0.5 * sqrt(Real(i + 1) / ip5);
				for (int j = 0; j <= i; j++) {
					errfac[i] *= e / (j + 1);
				}
			}
			ipoint[0] = 0;
			for (int i = 1; i <= IMAXX; i++) {
				int njadd = 4 * i - 2;
				if (nseq[i - 1] > njadd) njadd++;
				ipoint[i] = ipoint[i - 1] + njadd;
			}
		}

		void step(const Real htry) {
			const Real STEPFAC1 = 0.65, STEPFAC2 = 0.94, STEPFAC3 = 0.02, STEPFAC4 = 4.0,
				KFAC1 = 0.8, KFAC2 = 0.9;
			static bool first_step = true, last_step = false;
			static bool forward, reject = false, prev_reject = false;
			int i, k;
			Real fac, h, hnew, hopt_int, err;
			bool firstk;
			Vector<Real> hopt(IMAXX), work(IMAXX);
			Vector<Real> ysav(n), yseq(n);
			Vector<Real> ymid(n), scale(n);
			work[0] = 0;
			h = htry;
			forward = h > 0 ? true : false;
			for (i = 0; i < n; i++) ysav[i] = y[i];
			if (h != hnext && !first_step) {
				last_step = true;
			}
			if (reject) {
				prev_reject = true;
				last_step = false;
			}
			reject = false;
			firstk = true;
			hnew = std::abs(h);
		interp_error:
			while (firstk || reject) {
				h = forward ? hnew : -hnew;
				firstk = false;
				reject = false;
				if (std::abs(h) <= std::abs(x) * EPS)
					throw("step size underflow in StepperBS");
				int ipt = -1;
				for (k = 0; k <= k_targ + 1; k++) {
					dy(ysav, h, k, yseq, ipt);
					if (k == 0)
						y = yseq;
					else
						for (i = 0; i < n; i++)
							table[k - 1][i] = yseq[i];
					if (k != 0) {
						polyextr(k, table, y);
						err = 0.0;
						for (i = 0; i < n; i++) {
							scale[i] = atol + rtol * std::max(std::abs(ysav[i]), std::abs(y[i]));
							err += POW2((y[i] - table[0][i]) / scale[i]);
						}
						err = sqrt(err / n);
						Real expo = 1.0 / (2 * k + 1);
						Real facmin = pow(STEPFAC3, expo);
						if (err == 0.0)
							fac = 1.0 / facmin;
						else {
							fac = STEPFAC2 / pow(err / STEPFAC1, expo);
							fac = std::max<Real>(facmin / STEPFAC4, std::min<Real>(1.0 / facmin, fac));
						}
						hopt[k] = std::abs(h * fac);
						work[k] = cost[k] / hopt[k];
						if ((first_step || last_step) && err <= 1.0)
							break;
						if (k == k_targ - 1 && !prev_reject && !first_step && !last_step) {
							if (err <= 1.0)
								break;
							else if (err > POW2(nseq[k_targ] * nseq[k_targ + 1] / (nseq[0] * nseq[0]))) {
								reject = true;
								k_targ = k;
								if (k_targ > 1 && work[k - 1] < KFAC1 * work[k])
									k_targ--;
								hnew = hopt[k_targ];
								break;
							}
						}
						if (k == k_targ) {
							if (err <= 1.0)
								break;
							else if (err > POW2(nseq[k + 1] / nseq[0])) {
								reject = true;
								if (k_targ > 1 && work[k - 1] < KFAC1 * work[k])
									k_targ--;
								hnew = hopt[k_targ];
								break;
							}
						}
						if (k == k_targ + 1) {
							if (err > 1.0) {
								reject = true;
								if (k_targ > 1 && work[k_targ - 1] < KFAC1 * work[k_targ])
									k_targ--;
								hnew = hopt[k_targ];
							}
							break;
						}
					}
				}
				if (reject)
					prev_reject = true;
			}
			_sys.derivs(x + h, y, dydxnew);
			if (dense) {
				prepare_dense(h, dydxnew, ysav, scale, k, err);
				hopt_int = h / std::max<Real>(pow(err, 1.0 / (2 * k + 3)), 0.01);
				if (err > 10.0) {
					hnew = std::abs(hopt_int);
					reject = true;
					prev_reject = true;
					goto interp_error;
				}
			}
			dydx = dydxnew;
			xold = x;
			x += h;
			hdid = h;
			first_step = false;
			int kopt;
			if (k == 1)
				kopt = 2;
			else if (k <= k_targ) {
				kopt = k;
				if (work[k - 1] < KFAC1 * work[k])
					kopt = k - 1;
				else if (work[k] < KFAC2 * work[k - 1])
					kopt = std::min(k + 1, KMAXX - 1);
			}
			else {
				kopt = k - 1;
				if (k > 2 && work[k - 2] < KFAC1 * work[k - 1])
					kopt = k - 2;
				if (work[k] < KFAC2 * work[kopt])
					kopt = std::min(k, KMAXX - 1);
			}
			if (prev_reject) {
				k_targ = std::min(kopt, k);
				hnew = std::min(std::abs(h), hopt[k_targ]);
				prev_reject = false;
			}
			else {
				if (kopt <= k)
					hnew = hopt[kopt];
				else {
					if (k < k_targ && work[k] < KFAC2 * work[k - 1])
						hnew = hopt[k] * cost[kopt + 1] / cost[k];
					else
						hnew = hopt[k] * cost[kopt] / cost[k];
				}
				k_targ = kopt;
			}
			if (dense)
				hnew = std::min(hnew, std::abs(hopt_int));
			if (forward)
				hnext = hnew;
			else
				hnext = -hnew;
		}

		virtual void dy(const Vector<Real>& y, const Real htot, const int k, Vector<Real>& yend, int& ipt)
		{
			Vector<Real> ym(n), yn(n);
			int nstep = nseq[k];
			Real h = htot / nstep;
			for (int i = 0; i < n; i++) {
				ym[i] = y[i];
				yn[i] = y[i] + h * dydx[i];
			}
			Real xnew = x + h;
			_sys.derivs(xnew, yn, yend);
			Real h2 = 2.0 * h;
			for (int nn = 1; nn < nstep; nn++) {
				if (dense && nn == nstep / 2) {
					for (int i = 0; i < n; i++)
						ysave[k][i] = yn[i];
				}
				if (dense && std::abs(nn - nstep / 2) <= 2 * k + 1) {
					ipt++;
					for (int i = 0; i < n; i++)
						fsave[ipt][i] = yend[i];
				}
				for (int i = 0; i < n; i++) {
					Real swap = ym[i] + h2 * yend[i];
					ym[i] = yn[i];
					yn[i] = swap;
				}
				xnew += h;
				_sys.derivs(xnew, yn, yend);
			}
			if (dense && nstep / 2 <= 2 * k + 1) {
				ipt++;
				for (int i = 0; i < n; i++)
					fsave[ipt][i] = yend[i];
			}
			for (int i = 0; i < n; i++)
				yend[i] = 0.5 * (ym[i] + yn[i] + h * yend[i]);
		}

		void polyextr(const int k, Matrix<Real>& table, Vector<Real>& last)
		{
			int l = (int)last.size();
			for (int j = k - 1; j > 0; j--)
				for (int i = 0; i < l; i++)
					table[j - 1][i] = table[j][i] + coeff[k][j] * (table[j][i] - table[j - 1][i]);
			for (int i = 0; i < l; i++)
				last[i] = table[0][i] + coeff[k][0] * (table[0][i] - last[i]);
		}

		virtual void prepare_dense(const Real h, const Vector<Real>& dydxnew, const Vector<Real>& ysav, const Vector<Real>& scale, const int k, Real& error)
		{
			mu = 2 * k - 1;
			for (int i = 0; i < n; i++) {
				dens[i] = ysav[i];
				dens[n + i] = h * dydx[i];
				dens[2 * n + i] = y[i];
				dens[3 * n + i] = dydxnew[i] * h;
			}
			for (int j = 1; j <= k; j++) {
				Real dblenj = nseq[j];
				for (int l = j; l >= 1; l--) {
					Real factor = POW2(dblenj / nseq[l - 1]) - 1.0;
					for (int i = 0; i < n; i++)
						ysave[l - 1][i] = ysave[l][i] + (ysave[l][i] - ysave[l - 1][i]) / factor;
				}
			}
			for (int i = 0; i < n; i++)
				dens[4 * n + i] = ysave[0][i];
			for (int kmi = 1; kmi <= mu; kmi++) {
				int kbeg = (kmi - 1) / 2;
				for (int kk = kbeg; kk <= k; kk++) {
					Real facnj = pow(nseq[kk] / 2.0, kmi - 1);
					int ipt = ipoint[kk + 1] - 2 * kk + kmi - 3;
					for (int i = 0; i < n; i++)
						ysave[kk][i] = fsave[ipt][i] * facnj;
				}
				for (int j = kbeg + 1; j <= k; j++) {
					Real dblenj = nseq[j];
					for (int l = j; l >= kbeg + 1; l--) {
						Real factor = POW2(dblenj / nseq[l - 1]) - 1.0;
						for (int i = 0; i < n; i++)
							ysave[l - 1][i] = ysave[l][i] +
							(ysave[l][i] - ysave[l - 1][i]) / factor;
					}
				}
				for (int i = 0; i < n; i++)
					dens[(kmi + 4) * n + i] = ysave[kbeg][i] * h;
				if (kmi == mu) continue;
				for (int kk = kmi / 2; kk <= k; kk++) {
					int lbeg = ipoint[kk + 1] - 1;
					int lend = ipoint[kk] + kmi;
					if (kmi == 1) lend += 2;
					for (int l = lbeg; l >= lend; l -= 2)
						for (int i = 0; i < n; i++)
							fsave[l][i] = fsave[l][i] - fsave[l - 2][i];
					if (kmi == 1) {
						int l = lend - 2;
						for (int i = 0; i < n; i++)
							fsave[l][i] = fsave[l][i] - dydx[i];
					}
				}
				for (int kk = kmi / 2; kk <= k; kk++) {
					int lbeg = ipoint[kk + 1] - 2;
					int lend = ipoint[kk] + kmi + 1;
					for (int l = lbeg; l >= lend; l -= 2)
						for (int i = 0; i < n; i++)
							fsave[l][i] = fsave[l][i] - fsave[l - 2][i];
				}
			}
			dense_interp(n, dens, mu);
			error = 0.0;
			if (mu >= 1) {
				for (int i = 0; i < n; i++)
					error += POW2(dens[(mu + 4) * n + i] / scale[i]);
				error = sqrt(error / n) * errfac[mu - 1];
			}
		}
		virtual Real dense_out(const int i, const Real x, const Real h) {
			Real theta = (x - xold) / h;
			Real theta1 = 1.0 - theta;
			Real yinterp = dens[i] + theta * (dens[n + i] + theta1 * (dens[2 * n + i] * theta
				+ dens[3 * n + i] * theta1));
			if (mu < 0)
				return yinterp;
			Real theta05 = theta - 0.5;
			Real t4 = POW2(theta * theta1);
			Real c = dens[n * (mu + 4) + i];
			for (int j = mu; j > 0; j--)
				c = dens[n * (j + 3) + i] + c * theta05 / j;
			yinterp += t4 * c;
			return yinterp;
		}

		virtual void dense_interp(const int n, Vector<Real>& y, const int imit) {
			Real y0, y1, yp0, yp1, ydiff, aspl, bspl, ph0, ph1, ph2, ph3, fac1, fac2;
			Vector<Real> a(31);
			for (int i = 0; i < n; i++) {
				y0 = y[i];
				y1 = y[2 * n + i];
				yp0 = y[n + i];
				yp1 = y[3 * n + i];
				ydiff = y1 - y0;
				aspl = -yp1 + ydiff;
				bspl = yp0 - ydiff;
				y[n + i] = ydiff;
				y[2 * n + i] = aspl;
				y[3 * n + i] = bspl;
				if (imit < 0) continue;
				ph0 = (y0 + y1) * 0.5 + 0.125 * (aspl + bspl);
				ph1 = ydiff + (aspl - bspl) * 0.25;
				ph2 = -(yp0 - yp1);
				ph3 = 6.0 * (bspl - aspl);
				if (imit >= 1) {
					a[1] = 16.0 * (y[5 * n + i] - ph1);
					if (imit >= 3) {
						a[3] = 16.0 * (y[7 * n + i] - ph3 + 3 * a[1]);
						for (int im = 5; im <= imit; im += 2) {
							fac1 = im * (im - 1) / 2.0;
							fac2 = fac1 * (im - 2) * (im - 3) * 2.0;
							a[im] = 16.0 * (y[(im + 4) * n + i] + fac1 * a[im - 2] - fac2 * a[im - 4]);
						}
					}
				}
				a[0] = (y[4 * n + i] - ph0) * 16.0;
				if (imit >= 2) {
					a[2] = (y[n * 6 + i] - ph2 + a[0]) * 16.0;
					for (int im = 4; im <= imit; im += 2) {
						fac1 = im * (im - 1) / 2.0;
						fac2 = im * (im - 1) * (im - 2) * (im - 3);
						a[im] = (y[n * (im + 4) + i] + a[im - 2] * fac1 - a[im - 4] * fac2) * 16.0;
					}
				}
				for (int im = 0; im <= imit; im++)
					y[n * (im + 4) + i] = a[im];
			}
		}
	};
}
///////////////////////////   ./include/algorithms/ODESystemSteppers_Stiff.h   ///////////////////////////





namespace MML
{
	struct Ross_constants {
		const Real c2 = 0.386;
		const Real c3 = 0.21;
		const Real c4 = 0.63;
		const Real bet2p = 0.0317;
		const Real bet3p = 0.0635;
		const Real bet4p = 0.3438;
		const Real d1 = 0.2500000000000000e+00;
		const Real d2 = -0.1043000000000000e+00;
		const Real d3 = 0.1035000000000000e+00;
		const Real d4 = -0.3620000000000023e-01;
		const Real a21 = 0.1544000000000000e+01;
		const Real a31 = 0.9466785280815826e+00;
		const Real a32 = 0.2557011698983284e+00;
		const Real a41 = 0.3314825187068521e+01;
		const Real a42 = 0.2896124015972201e+01;
		const Real a43 = 0.9986419139977817e+00;
		const Real a51 = 0.1221224509226641e+01;
		const Real a52 = 0.6019134481288629e+01;
		const Real a53 = 0.1253708332932087e+02;
		const Real a54 = -0.6878860361058950e+00;
		const Real c21 = -0.5668800000000000e+01;
		const Real c31 = -0.2430093356833875e+01;
		const Real c32 = -0.2063599157091915e+00;
		const Real c41 = -0.1073529058151375e+00;
		const Real c42 = -0.9594562251023355e+01;
		const Real c43 = -0.2047028614809616e+02;
		const Real c51 = 0.7496443313967647e+01;
		const Real c52 = -0.1024680431464352e+02;
		const Real c53 = -0.3399990352819905e+02;
		const Real c54 = 0.1170890893206160e+02;
		const Real c61 = 0.8083246795921522e+01;
		const Real c62 = -0.7981132988064893e+01;
		const Real c63 = -0.3152159432874371e+02;
		const Real c64 = 0.1631930543123136e+02;
		const Real c65 = -0.6058818238834054e+01;
		const Real gam = 0.2500000000000000e+00;
		const Real d21 = 0.1012623508344586e+02;
		const Real d22 = -0.7487995877610167e+01;
		const Real d23 = -0.3480091861555747e+02;
		const Real d24 = -0.7992771707568823e+01;
		const Real d25 = 0.1025137723295662e+01;
		const Real d31 = -0.6762803392801253e+00;
		const Real d32 = 0.6087714651680015e+01;
		const Real d33 = 0.1643084320892478e+02;
		const Real d34 = 0.2476722511418386e+02;
		const Real d35 = -0.6594389125716872e+01;
	};

	// Fourth-order stiffly stable Rosenbrock step for integrating stiff ODEs, with monitoring of local
	// truncation error to adjust stepsize.
	struct StepperRoss : StepperBase, Ross_constants {
		Matrix<Real> dfdy;
		Vector<Real> dfdx;
		Vector<Real> k1, k2, k3, k4, k5, k6;
		Vector<Real> cont1, cont2, cont3, cont4;
		Matrix<Real> a;

		// Input to the constructor are the dependent variable y[0..n-1] and its derivative dydx[0..n-1]
		// at the starting value of the independent variable x. Also input are the absolute and relative
		// tolerances, atol and rtol, and the boolean dense, which is true if dense output is required
		StepperRoss(IODESystem& sys, Vector<Real>& yy, Vector<Real>& dydxx, Real& xx, const Real atoll, const Real rtoll, bool dens)
			: StepperBase(sys, yy, dydxx, xx, atoll, rtoll, dens),
			dfdy(n, n), dfdx(n), k1(n), k2(n),
			k3(n), k4(n), k5(n), k6(n), cont1(n), cont2(n), cont3(n), cont4(n), a(n, n)
		{
			EPS = std::numeric_limits<Real>::epsilon();
		}

		// Attempts a step with stepsize htry. On output, y and x are replaced by their new values, hdid
		// is the stepsize that was actually accomplished, and hnext is the estimated next stepsize.        
		void step(const Real htry) {
			Vector<Real> dydxnew(n);
			Real h = htry;
			dynamic_cast<ODESystemWithJacobian&>(_sys).jacobian(x, y, dfdx, dfdy);
			for (;;) {
				dy(h);
				Real err = error();
				if (con.success(err, h)) break;
				if (std::abs(h) <= std::abs(x) * EPS)
					throw("stepsize underflow in StepperRoss");
			}
			_sys.derivs(x + h, yout, dydxnew);
			if (dense)
				prepare_dense(h, dydxnew);
			dydx = dydxnew;
			y = yout;
			xold = x;
			x += (hdid = h);
			hnext = con.hnext;
		}

		// Given values for n variables y[0..n-1] and their derivatives dydx[0..n-1] known at x, use the
		// fourth-order stiffly stable Rosenbrock method to advance the solution over an interval h and
		// store the incremented variables in yout[0..n-1]. Also store an estimate of the local truncation
		// error in yerr using the embedded third-order method
		void dy(const Real h) {
			Vector<Real> ytemp(n), dydxnew(n);
			int i;
			for (i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) a[i][j] = -dfdy[i][j];
				a[i][i] += 1.0 / (gam * h);
			}
			LUDecompositionSolver<Real> alu(a);
			for (i = 0; i < n; i++)
				ytemp[i] = dydx[i] + h * d1 * dfdx[i];
			alu.Solve(ytemp, k1);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + a21 * k1[i];
			_sys.derivs(x + c2 * h, ytemp, dydxnew);
			for (i = 0; i < n; i++)
				ytemp[i] = dydxnew[i] + h * d2 * dfdx[i] + c21 * k1[i] / h;
			alu.Solve(ytemp, k2);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + a31 * k1[i] + a32 * k2[i];
			_sys.derivs(x + c3 * h, ytemp, dydxnew);
			for (i = 0; i < n; i++)
				ytemp[i] = dydxnew[i] + h * d3 * dfdx[i] + (c31 * k1[i] + c32 * k2[i]) / h;
			alu.Solve(ytemp, k3);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + a41 * k1[i] + a42 * k2[i] + a43 * k3[i];
			_sys.derivs(x + c4 * h, ytemp, dydxnew);
			for (i = 0; i < n; i++)
				ytemp[i] = dydxnew[i] + h * d4 * dfdx[i] + (c41 * k1[i] + c42 * k2[i] + c43 * k3[i]) / h;
			alu.Solve(ytemp, k4);
			for (i = 0; i < n; i++)
				ytemp[i] = y[i] + a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i];
			Real xph = x + h;
			_sys.derivs(xph, ytemp, dydxnew);
			for (i = 0; i < n; i++)
				k6[i] = dydxnew[i] + (c51 * k1[i] + c52 * k2[i] + c53 * k3[i] + c54 * k4[i]) / h;
			alu.Solve(k6, k5);
			for (i = 0; i < n; i++)
				ytemp[i] += k5[i];
			_sys.derivs(xph, ytemp, dydxnew);
			for (i = 0; i < n; i++)
				k6[i] = dydxnew[i] + (c61 * k1[i] + c62 * k2[i] + c63 * k3[i] + c64 * k4[i] + c65 * k5[i]) / h;
			alu.Solve(k6, yerr);
			for (i = 0; i < n; i++)
				yout[i] = ytemp[i] + yerr[i];
		}

		// Store coefficients of interpolating polynomial for dense output in cont1...cont4.
		void prepare_dense(const Real h, const Vector<Real>& dydxnew) {
			for (int i = 0; i < n; i++) {
				cont1[i] = y[i];
				cont2[i] = yout[i];
				cont3[i] = d21 * k1[i] + d22 * k2[i] + d23 * k3[i] + d24 * k4[i] + d25 * k5[i];
				cont4[i] = d31 * k1[i] + d32 * k2[i] + d33 * k3[i] + d34 * k4[i] + d35 * k5[i];
			}
		}

		// Evaluate interpolating polynomial for y[i] at location x, where xold <= x <= xold + h        
		Real dense_out(const int i, const Real x, const Real h) {
			Real s = (x - xold) / h;
			Real s1 = 1.0 - s;
			return cont1[i] * s1 + s * (cont2[i] + s1 * (cont3[i] + s * cont4[i]));
		}

		// Use yerr to compute norm of scaled error estimate. A value less than one means the step was successful        
		Real error() {
			Real err = 0.0, sk;
			for (int i = 0; i < n; i++) {
				sk = atol + rtol * std::max(std::abs(y[i]), std::abs(yout[i]));
				err += POW2(yerr[i] / sk);
			}
			return sqrt(err / n);
		}

		struct Controller {
			Real hnext;
			bool reject;
			bool first_step;
			Real errold;
			Real hold;
			Controller() : reject(false), first_step(true) {}
			bool success(Real err, Real& h) {
				static const Real safe = 0.9, fac1 = 5.0, fac2 = 1.0 / 6.0;
				Real fac = std::max<Real>(fac2, std::min<Real>(fac1, pow(err, 0.25) / safe));
				Real hnew = h / fac;
				if (err <= 1.0) {
					if (!first_step) {
						Real facpred = (hold / h) * pow(err * err / errold, 0.25) / safe;
						facpred = std::max(fac2, std::min(fac1, facpred));
						fac = std::max(fac, facpred);
						hnew = h / fac;
					}
					first_step = false;
					hold = h;
					errold = std::max<Real>(0.01, err);
					if (reject)
						hnew = (h >= 0.0 ? std::min(hnew, h) : std::max(hnew, h));
					hnext = hnew;
					reject = false;
					return true;
				}
				else {
					h = hnew;
					reject = true;
					return false;
				}
			}
		};
		Controller con;
	};

	struct StepperSemiImplExtr : StepperBase {
		static const int KMAXX = 12, IMAXX = KMAXX + 1;
		int k_targ;
		Vector<int> nseq;
		Vector<Real> cost;
		Matrix<Real> table;
		Matrix<Real> dfdy;
		Vector<Real> dfdx;
		Real jac_redo;
		bool calcjac;
		Real theta;
		Matrix<Real> a;
		int kright;
		Matrix<Real> coeff;
		Matrix<Real> fsave;
		Vector<Real> dens;
		Vector<Real> factrl;

		StepperSemiImplExtr(IODESystem& sys, Vector<Real>& yy, Vector<Real>& dydxx, Real& xx,
			const Real atoll, const Real rtoll, bool dens)
			: StepperBase(sys, yy, dydxx, xx, atoll, rtoll, dens), nseq(IMAXX), cost(IMAXX),
			table(KMAXX, n), dfdy(n, n), dfdx(n), calcjac(false),
			a(n, n), coeff(IMAXX, IMAXX),
			fsave((IMAXX - 1)* (IMAXX + 1) / 2 + 2, n), dens((IMAXX + 2)* n), factrl(IMAXX)
		{
			static const Real costfunc = 1.0, costjac = 5.0, costlu = 1.0, costsolve = 1.0;
			EPS = std::numeric_limits<Real>::epsilon();
			jac_redo = std::min<Real>(1.0e-4, rtol);
			theta = 2.0 * jac_redo;
			nseq[0] = 2;
			nseq[1] = 3;
			for (int i = 2; i < IMAXX; i++)
				nseq[i] = 2 * nseq[i - 2];
			cost[0] = costjac + costlu + nseq[0] * (costfunc + costsolve);
			for (int k = 0; k < KMAXX; k++)
				cost[k + 1] = cost[k] + (nseq[k + 1] - 1) * (costfunc + costsolve) + costlu;
			hnext = -1.0e99;
			Real logfact = -log10(rtol + atol) * 0.6 + 0.5;
			k_targ = std::max(1, std::min(KMAXX - 1, int(logfact)));
			for (int k = 0; k < IMAXX; k++) {
				for (int l = 0; l < k; l++) {
					Real ratio = Real(nseq[k]) / nseq[l];
					coeff[k][l] = 1.0 / (ratio - 1.0);
				}
			}
			factrl[0] = 1.0;
			for (int k = 0; k < IMAXX - 1; k++)
				factrl[k + 1] = (k + 1) * factrl[k];
		}
		void step(const Real htry) {
			const Real STEPFAC1 = 0.6, STEPFAC2 = 0.93, STEPFAC3 = 0.1, STEPFAC4 = 4.0,
				STEPFAC5 = 0.5, KFAC1 = 0.7, KFAC2 = 0.9;
			static bool first_step = true, last_step = false;
			static bool forward, reject = false, prev_reject = false;
			static Real errold;
			int i, k;
			Real fac, h, hnew, err;
			bool firstk;
			Vector<Real> hopt(IMAXX), work(IMAXX);
			Vector<Real> ysav(n), yseq(n);
			Vector<Real> ymid(n), scale(n);
			work[0] = 1.e30;
			h = htry;
			forward = h > 0 ? true : false;
			for (i = 0; i < n; i++) ysav[i] = y[i];
			if (h != hnext && !first_step) {
				last_step = true;
			}
			if (reject) {
				prev_reject = true;
				last_step = false;
				theta = 2.0 * jac_redo;
			}
			for (i = 0; i < n; i++)
				scale[i] = atol + rtol * std::abs(y[i]);
			reject = false;
			firstk = true;
			hnew = std::abs(h);
		compute_jac:
			if (theta > jac_redo && !calcjac) {
				dynamic_cast<ODESystemWithJacobian&>(_sys).jacobian(x, y, dfdx, dfdy);

				//_sys.jacobian(x,y,dfdx,dfdy);
				calcjac = true;
			}
			while (firstk || reject) {
				h = forward ? hnew : -hnew;
				firstk = false;
				reject = false;
				if (std::abs(h) <= std::abs(x) * EPS)
					throw("step size underflow in StepperSemiImplExtr");
				int ipt = -1;
				for (k = 0; k <= k_targ + 1; k++) {
					bool success = dy(ysav, h, k, yseq, ipt, scale);
					if (!success) {
						reject = true;
						hnew = std::abs(h) * STEPFAC5;
						break;
					}
					if (k == 0)
						y = yseq;
					else
						for (i = 0; i < n; i++)
							table[k - 1][i] = yseq[i];
					if (k != 0) {
						polyextr(k, table, y);
						err = 0.0;
						for (i = 0; i < n; i++) {
							scale[i] = atol + rtol * std::abs(ysav[i]);
							err += POW2((y[i] - table[0][i]) / scale[i]);
						}
						err = sqrt(err / n);
						if (err > 1.0 / EPS || (k > 1 && err >= errold)) {
							reject = true;
							hnew = std::abs(h) * STEPFAC5;
							break;
						}
						errold = std::max<Real>(4.0 * err, 1.0);
						Real expo = 1.0 / (k + 1);
						Real facmin = pow(STEPFAC3, expo);
						if (err == 0.0)
							fac = 1.0 / facmin;
						else {
							fac = STEPFAC2 / pow(err / STEPFAC1, expo);
							fac = std::max<Real>(facmin / STEPFAC4, std::min<Real>(1.0 / facmin, fac));
						}
						hopt[k] = std::abs(h * fac);
						work[k] = cost[k] / hopt[k];
						if ((first_step || last_step) && err <= 1.0)
							break;
						if (k == k_targ - 1 && !prev_reject && !first_step && !last_step) {
							if (err <= 1.0)
								break;
							else if (err > nseq[k_targ] * nseq[k_targ + 1] * 4.0) {
								reject = true;
								k_targ = k;
								if (k_targ > 1 && work[k - 1] < KFAC1 * work[k])
									k_targ--;
								hnew = hopt[k_targ];
								break;
							}
						}
						if (k == k_targ) {
							if (err <= 1.0)
								break;
							else if (err > nseq[k + 1] * 2.0) {
								reject = true;
								if (k_targ > 1 && work[k - 1] < KFAC1 * work[k])
									k_targ--;
								hnew = hopt[k_targ];
								break;
							}
						}
						if (k == k_targ + 1) {
							if (err > 1.0) {
								reject = true;
								if (k_targ > 1 && work[k_targ - 1] < KFAC1 * work[k_targ])
									k_targ--;
								hnew = hopt[k_targ];
							}
							break;
						}
					}
				}
				if (reject) {
					prev_reject = true;
					if (!calcjac) {
						theta = 2.0 * jac_redo;
						goto compute_jac;
					}
				}
			}
			calcjac = false;
			if (dense)
				prepare_dense(h, ysav, scale, k, err);
			xold = x;
			x += h;
			hdid = h;
			first_step = false;
			int kopt;
			if (k == 1)
				kopt = 2;
			else if (k <= k_targ) {
				kopt = k;
				if (work[k - 1] < KFAC1 * work[k])
					kopt = k - 1;
				else if (work[k] < KFAC2 * work[k - 1])
					kopt = std::min(k + 1, KMAXX - 1);
			}
			else {
				kopt = k - 1;
				if (k > 2 && work[k - 2] < KFAC1 * work[k - 1])
					kopt = k - 2;
				if (work[k] < KFAC2 * work[kopt])
					kopt = std::min(k, KMAXX - 1);
			}
			if (prev_reject) {
				k_targ = std::min(kopt, k);
				hnew = std::min(std::abs(h), hopt[k_targ]);
				prev_reject = false;
			}
			else {
				if (kopt <= k)
					hnew = hopt[kopt];
				else {
					if (k < k_targ && work[k] < KFAC2 * work[k - 1])
						hnew = hopt[k] * cost[kopt + 1] / cost[k];
					else
						hnew = hopt[k] * cost[kopt] / cost[k];
				}
				k_targ = kopt;
			}
			if (forward)
				hnext = hnew;
			else
				hnext = -hnew;
		}

		bool dy(const Vector<Real>& y, const Real htot, const int k, Vector<Real>& yend, int& ipt, const Vector<Real>& scale)
		{
			Vector<Real> del(n), ytemp(n), dytemp(n);
			int nstep = nseq[k];
			Real h = htot / nstep;
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) a[i][j] = -dfdy[i][j];
				a[i][i] += 1.0 / h;
			}
			LUDecompositionSolver<Real> alu(a);
			Real xnew = x + h;
			_sys.derivs(xnew, y, del);
			for (int i = 0; i < n; i++)
				ytemp[i] = y[i];
			alu.Solve(del, del);
			if (dense && nstep == k + 1) {
				ipt++;
				for (int i = 0; i < n; i++)
					fsave[ipt][i] = del[i];
			}
			for (int nn = 1; nn < nstep; nn++) {
				for (int i = 0; i < n; i++)
					ytemp[i] += del[i];
				xnew += h;
				_sys.derivs(xnew, ytemp, yend);
				if (nn == 1 && k <= 1) {
					Real del1 = 0.0;
					for (int i = 0; i < n; i++)
						del1 += POW2(del[i] / scale[i]);
					del1 = sqrt(del1);
					_sys.derivs(x + h, ytemp, dytemp);
					for (int i = 0; i < n; i++)
						del[i] = dytemp[i] - del[i] / h;
					alu.Solve(del, del);
					Real del2 = 0.0;
					for (int i = 0; i < n; i++)
						del2 += POW2(del[i] / scale[i]);
					del2 = sqrt(del2);
					theta = del2 / std::max<Real>(1.0, del1);
					if (theta > 1.0)
						return false;
				}
				alu.Solve(yend, del);
				if (dense && nn >= nstep - k - 1) {
					ipt++;
					for (int i = 0; i < n; i++)
						fsave[ipt][i] = del[i];
				}
			}
			for (int i = 0; i < n; i++)
				yend[i] = ytemp[i] + del[i];
			return true;
		}
		void polyextr(const int k, Matrix<Real>& table, Vector<Real>& last) {
			int l = (int)last.size();
			for (int j = k - 1; j > 0; j--)
				for (int i = 0; i < l; i++)
					table[j - 1][i] = table[j][i] + coeff[k][j] * (table[j][i] - table[j - 1][i]);
			for (int i = 0; i < l; i++)
				last[i] = table[0][i] + coeff[k][0] * (table[0][i] - last[i]);
		}

		void prepare_dense(const Real h, const Vector<Real>& ysav, const Vector<Real>& scale,
			const int k, Real& error) {
			kright = k;
			for (int i = 0; i < n; i++) {
				dens[i] = ysav[i];
				dens[n + i] = y[i];
			}
			for (int klr = 0; klr < kright; klr++) {
				if (klr >= 1) {
					for (int kk = klr; kk <= k; kk++) {
						int lbeg = ((kk + 3) * kk) / 2;
						int lend = lbeg - kk + 1;
						for (int l = lbeg; l >= lend; l--)
							for (int i = 0; i < n; i++)
								fsave[l][i] = fsave[l][i] - fsave[l - 1][i];
					}
				}
				for (int kk = klr; kk <= k; kk++) {
					Real facnj = nseq[kk];
					facnj = pow(facnj, klr + 1) / factrl[klr + 1];
					int ipt = ((kk + 3) * kk) / 2;
					int krn = (kk + 2) * n;
					for (int i = 0; i < n; i++) {
						dens[krn + i] = fsave[ipt][i] * facnj;
					}
				}
				for (int j = klr + 1; j <= k; j++) {
					Real dblenj = nseq[j];
					for (int l = j; l >= klr + 1; l--) {
						Real factor = dblenj / nseq[l - 1] - 1.0;
						for (int i = 0; i < n; i++) {
							int krn = (l + 2) * n + i;
							dens[krn - n] = dens[krn] + (dens[krn] - dens[krn - n]) / factor;
						}
					}
				}
			}
			for (int in = 0; in < n; in++) {
				for (int j = 1; j <= kright + 1; j++) {
					int ii = n * j + in;
					dens[ii] = dens[ii] - dens[ii - n];
				}
			}
		}
		Real dense_out(const int i, const Real x, const Real h) {
			Real theta = (x - xold) / h;
			int k = kright;
			Real yinterp = dens[(k + 1) * n + i];
			for (int j = 1; j <= k; j++)
				yinterp = dens[(k + 1 - j) * n + i] + yinterp * (theta - 1.0);
			return dens[i] + yinterp * theta;
		}
		void dense_interp(const int n, Vector<Real>& y, const int imit);
	};
}
///////////////////////////   ./include/algorithms/ODESystemSolver.h   ///////////////////////////




namespace MML
{
	struct  Output {
		// TODO - MED, LAKO, sve privatizirati
		int kmax;
		int _nvar;
		int nsave;

		int nok;
		int nbad;
		int nstp;

		bool dense;
		int count;

		Real x1, x2, xout, dxout;

		Vector<Real> xsave;
		Matrix<Real> ysave;

		Output() : kmax(-1), dense(false), count(0) {}

		// Constructor provides dense output at nsave equally spaced intervals. If nsave <= 0, output
		// is saved only at the actual integration steps.        
		Output(const int nsavee) : kmax(50), nsave(nsavee), count(0), xsave(kmax), nok(0), nbad(0) {
			dense = nsave > 0 ? true : false;
		}
		void init(const int neqn, const Real xlo, const Real xhi) {
			_nvar = neqn;
			if (kmax == -1) return;
			ysave.Resize(_nvar, kmax);
			if (dense) {
				x1 = xlo;
				x2 = xhi;
				xout = x1;
				dxout = (x2 - x1) / nsave;
			}
		}
		void resize() {
			int kold = kmax;
			kmax *= 2;

			Vector<Real> tempvec(xsave);
			xsave.Resize(kmax);
			for (int k = 0; k < kold; k++)
				xsave[k] = tempvec[k];

			Matrix<Real> tempmat(ysave);
			ysave.Resize(_nvar, kmax);
			for (int i = 0; i < _nvar; i++)
				for (int k = 0; k < kold; k++)
					ysave[i][k] = tempmat[i][k];
		}

		// Invokes dense_out function of stepper routine to produce output at xout. Normally called
		// by out rather than directly. Assumes that xout is between xold and xold+h, where the
		// stepper must keep track of xold, the location of the previous step, and x=xold+h, the
		// current step        
		void save_dense(StepperBase& s, const Real xout, const Real h) {
			if (count == kmax) resize();
			for (int i = 0; i < _nvar; i++)
				ysave[i][count] = s.dense_out(i, xout, h);
			xsave[count++] = xout;
		}

		// Saves values of current x and y.
		void save(const Real x, const Vector<Real>& y) {
			if (kmax <= 0) return;
			if (count == kmax) resize();
			for (int i = 0; i < _nvar; i++)
				ysave[i][count] = y[i];
			xsave[count++] = x;
		}

		// Typically called by Odeint to produce dense output. Input variables are nstp, the current
		// step number, the current values of x and y, the stepper s, and the stepsize h. A call with
		// nstp=-1 saves the initial values. The routine checks whether x is greater than the desired
		// output point xout. If so, it calls save_dense.        
		void out(const int nstp, const Real x, const Vector<Real>& y, StepperBase& s, const Real h) {
			if (!dense)
				throw("dense output not set in Output!");
			if (nstp == -1) {
				save(x, y);
				xout += dxout;
			}
			else {
				while ((x - xout) * (x2 - x1) > 0.0) {
					save_dense(s, xout, h);
					xout += dxout;
				}
			}
		}
	};

	template<class Stepper>
	class ODESystemSolver {
		static inline const Real EPS = std::numeric_limits<Real>::epsilon();
		static const int MAXSTP = 50000;

		Real         _curr_x;               // used as reference by stepper!
		Vector<Real> _curr_y;
		Vector<Real> _curr_dydx;
	public:
		IODESystem& _sys;
		Output&			_out;
		Stepper    _stepper;
		bool _dense;

		int getDim() { return _sys.getDim(); }

		ODESystemSolver(IODESystem& sys, const Real atol, const Real rtol, Output& out)
								: _sys(sys), _curr_y(sys.getDim()), _curr_dydx(sys.getDim()), _dense(out.dense), _out(out),
									_stepper(sys, _curr_y, _curr_dydx, _curr_x, atol, rtol, out.dense)
		{ }

		ODESystemSolution integrate(Vector<Real>& in_ystart, const Real x1, Real x2, Real h1, const Real hmin)
		{
			int  dim = _sys.getDim();
			Real h = SIGN(h1, x2 - x1);

			ODESystemSolution sol(x1, x2, dim, _out.nsave);

			_out.init(_stepper.neqn, x1, x2);

			_curr_x = x1;
			_curr_y = in_ystart;
			Vector<Real>& ystart = in_ystart;

			_sys.derivs(_curr_x, _curr_y, _curr_dydx);

			if (_dense)
				_out.out(-1, _curr_x, _curr_y, _stepper, h);
			else
				_out.save(_curr_x, _curr_y);

			for (_out.nstp = 0; _out.nstp < MAXSTP; _out.nstp++)
			{
				if ((_curr_x + h * 1.0001 - x2) * (x2 - x1) > 0.0)
					h = x2 - _curr_x;

				_stepper.step(h);

				if (_stepper.hdid == h) ++_out.nok;
				else ++_out.nbad;

				if (_dense)
					_out.out(_out.nstp, _curr_x, _curr_y, _stepper, _stepper.hdid);
				else
					_out.save(_curr_x, _curr_y);

				if ((_curr_x - x2) * (x2 - x1) >= 0.0) {
					// we are done
					for (int i = 0; i < getDim(); i++)
						ystart[i] = _curr_y[i];

					if (_out.kmax > 0 && std::abs(_out.xsave[_out.count - 1] - x2) > 100.0 * std::abs(x2) * EPS)
						_out.save(_curr_x, _curr_y);

					for (int i = 0; i <= _out.nsave; i++)
					{
						sol._xval[i] = _out.xsave[i];
						for (int j = 0; j < ystart.size(); j++)
						{
							sol._yval[j][i] = _out.ysave[j][i];
						}
					}

					return sol;
				}
				if (std::abs(_stepper.hnext) <= hmin)
					throw("Step size too small in ODESystemSolver::integrate");

				h = _stepper.hnext;
			}
			throw("Too many steps in routine ODESystemSolver::integrate");
		}
	};

	class RungeKuttaSolverDumb
	{
	public:
		void rk4(Vector<Real>& y, Vector<Real>& dydx, const Real x, const Real h,
			Vector<Real>& yout, IODESystem& sys)
		{
			int i;
			Real xh, hh, h6;

			int n = (int)y.size();
			Vector<Real> dym(n), dyt(n), yt(n);
			hh = h * 0.5;
			h6 = h / 6.0;
			xh = x + hh;
			for (i = 0; i < n; i++) yt[i] = y[i] + hh * dydx[i];
			sys.derivs(xh, yt, dyt);
			for (i = 0; i < n; i++) yt[i] = y[i] + hh * dyt[i];
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
			for (k = 0; k < numSteps; k++) {
				sys.derivs(x, v, dv);
				rk4(v, dv, x, h, vout, sys);
				if (x + h == x)
					throw("Step size too small in routine rkdumb");
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

			int n = (int)y.size();
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
			for (nstp = 0; nstp < MAXSTP; nstp++) {
				sys.derivs(x, y, dydx);
				for (i = 0; i < dim; i++)
					yscal[i] = fabs(y[i]) + fabs(dydx[i] * h) + TINY;
				if (maxSteps > 0 && kount < maxSteps - 1 && fabs(x - xsav) > fabs(minSaveInterval)) {
					for (i = 0; i < dim; i++) sol._yval[i][kount] = y[i];
					sol._xval[kount++] = x;
					xsav = x;
				}
				if ((x + h - x2) * (x + h - x1) > 0.0) h = x2 - x;
				rkqs(y, dydx, x, h, eps, yscal, hdid, hnext, sys);
				if (hdid == h) ++nok; else ++nbad;
				if ((x - x2) * (x2 - x1) >= 0.0) {
					if (maxSteps != 0) {
						for (i = 0; i < dim; i++) sol._yval[i][kount] = y[i];
						sol._xval[kount++] = x;
					}
					return sol;
				}
				if (fabs(hnext) <= hmin) throw("Step size too small in integrate");
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
		// template<int N>
		// static VectorN<Real, N> getNormalScaled(const IParametricCurve<N>& curve, Real t)
		// {
		// 	CurveTangentUnit  helper(curve);

		// 	return Derivation::DeriveCurve<N>(helper, t, nullptr);
		// }
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
			auto   vec2 = y_der_2 - res1 * y_der_1.ScalarProductCartesian(y_der_2) * y_der_1;
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

///////////////////////////   ./include/algorithms/Fourier.h   ///////////////////////////

namespace MML
{
	// TODO - viditi sto s Fourierom
	class Fourier
	{

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
			xb1.Resize(nb);
			xb2.Resize(nb);
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
						xb1.Resize(2 * nb);
						xb2.Resize(2 * nb);
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


///////////////////////////   ./include/algorithms/RootFindingMultidim.h   ///////////////////////////﻿#if !defined MML_ROOTFINDING_MULTIDIM_H



namespace MML::RootFinding
{
	// TODO 0.9 - HIGH, SREDNJE!!! finish
	//Given an n - dimensional point xold[0..n - 1], the value of the function and gradient there, fold
	//and g[0..n - 1], and a direction p[0..n - 1], finds a new point x[0..n - 1] along the direction
	//p from xold where the function or functor func has decreased “sufficiently.” The new function
	//value is returned in f.stpmax is an input quantity that limits the length of the steps so that you
	//do not try to evaluate the function in regions where it is undefined or subject to overflow.p is
	//usually the Newton direction.The output quantity check is false on a normal exit.It is true
	//when x is too close to xold.In a minimization algorithm, this usually signals convergence and
	//can be ignored.However, in a zero - finding algorithm the calling program should check whether
	//the convergence is spurious.
	template <class T>
	void lnsrch(Vector<Real>& xold, const Real fold, Vector<Real>& g, Vector<Real>& p,
		Vector<Real>& x, Real& f, const Real stpmax, bool& check, T& func)
	{
		const Real ALF = 1.0e-4, TOLX = std::numeric_limits<Real>::epsilon();
		Real a, alam, alam2 = 0.0, alamin, b, disc, f2 = 0.0;
		Real rhs1, rhs2, slope = 0.0, sum = 0.0, temp, test, tmplam;
		int i, n = xold.size();
		check = false;
		for (i = 0; i < n; i++) sum += p[i] * p[i];
		sum = sqrt(sum);
		if (sum > stpmax)
			for (i = 0; i < n; i++)
				p[i] *= stpmax / sum;
		for (i = 0; i < n; i++)
			slope += g[i] * p[i];
		if (slope >= 0.0) throw("Roundoff problem in lnsrch.");
		test = 0.0;
		for (i = 0; i < n; i++) {
			temp = std::abs(p[i]) / std::max<Real>(std::abs(xold[i]), 1.0);
			if (temp > test) test = temp;
		}
		alamin = TOLX / test;
		alam = 1.0;
		for (;;) {
			for (i = 0; i < n; i++) x[i] = xold[i] + alam * p[i];
			f = func(x);
			if (alam < alamin) {
				for (i = 0; i < n; i++) x[i] = xold[i];
				check = true;
				return;
			}
			else if (f <= fold + ALF * alam * slope) return;
			else {
				if (alam == 1.0)
					tmplam = -slope / (2.0 * (f - fold - slope));
				else {
					rhs1 = f - fold - alam * slope;
					rhs2 = f2 - fold - alam2 * slope;
					a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2)) / (alam - alam2);
					b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2 / (alam2 * alam2)) / (alam - alam2);
					if (a == 0.0) tmplam = -slope / (2.0 * b);
					else {
						disc = b * b - 3.0 * a * slope;
						if (disc < 0.0) tmplam = 0.5 * alam;
						else if (b <= 0.0) tmplam = (-b + sqrt(disc)) / (3.0 * a);
						else tmplam = -slope / (b + sqrt(disc));
					}
					if (tmplam > 0.5 * alam)
						tmplam = 0.5 * alam;
				}
			}
			alam2 = alam;
			f2 = f;
			alam = std::max(tmplam, 0.1 * alam);
		}
	}
	template <class T>
	struct NRfdjac {
		const Real EPS;
		T& func;
		NRfdjac(T& funcc) : EPS(1.0e-8), func(funcc) {}
		Matrix<Real> operator() (Vector<Real>& x, Vector<Real>& fvec) {
			int n = x.size();
			Matrix<Real> df(n, n);
			Vector<Real> xh = x;
			for (int j = 0; j < n; j++) {
				Real temp = xh[j];
				Real h = EPS * std::abs(temp);
				if (h == 0.0) h = EPS;
				xh[j] = temp + h;
				h = xh[j] - temp;
				Vector<Real> f = func(xh);
				xh[j] = temp;
				for (int i = 0; i < n; i++)
					df[i][j] = (f[i] - fvec[i]) / h;
			}
			return df;
		}
	};
	template <class T>
	struct NRfmin {
		Vector<Real> fvec;
		T& func;
		int n;
		NRfmin(T& funcc) : func(funcc) {}
		Real operator() (Vector<Real>& x) {
			n = x.size();
			Real sum = 0;
			fvec = func(x);
			for (int i = 0; i < n; i++) sum += POW2(fvec[i]);
			return 0.5 * sum;
		}
	};

	//Given an initial guess x[0..n - 1] for a root in n dimensions, find the root by a globally convergent
	//Newton’s method.The vector of functions to be zeroed, called fvec[0..n - 1] in the routine
	//below, is returned by the user - supplied function or functor vecfunc(see text).The output
	//quantity check is false on a normal return and true if the routine has converged to a local
	//minimum of the function fmin defined below.In this case try restarting from a different initial
	//guess
	template <class T>
	void newt(Vector<Real>& x, bool& check, T& vecfunc) {
		const int MAXITS = 200;
		const Real TOLF = 1.0e-8, TOLMIN = 1.0e-12, STPMX = 100.0;
		const Real TOLX = std::numeric_limits<Real>::epsilon();
		//Here MAXITS is the maximum number of iterations; TOLF sets the convergence criterion on
		//function values; TOLMIN sets the criterion for deciding whether spurious convergence to a
		//minimum of fmin has occurred; STPMX is the scaled maximum step length allowed in line
		//searches; and TOLX is the convergence criterion on delta x
		int i, j, its, n = x.size();
		Real den, f, fold, stpmax, sum, temp, test;
		Vector<Real> g(n), p(n), xold(n);
		Matrix<Real> fjac(n, n);
		NRfmin<T> fmin(vecfunc);
		NRfdjac<T> fdjac(vecfunc);
		Vector<Real>& fvec = fmin.fvec;
		f = fmin(x);
		test = 0.0;
		for (i = 0; i < n; i++)
			if (std::abs(fvec[i]) > test) test = std::abs(fvec[i]);
		if (test < 0.01 * TOLF) {
			check = false;
			return;
		}
		sum = 0.0;
		for (i = 0; i < n; i++) sum += POW2(x[i]);
		stpmax = STPMX * std::max(sqrt(sum), Real(n));
		for (its = 0; its < MAXITS; its++) {
			fjac = fdjac(x, fvec);
			for (i = 0; i < n; i++) {
				sum = 0.0;
				for (j = 0; j < n; j++) sum += fjac[j][i] * fvec[j];
				g[i] = sum;
			}
			for (i = 0; i < n; i++) xold[i] = x[i];
			fold = f;
			for (i = 0; i < n; i++) p[i] = -fvec[i];

			LUDecompositionSolver alu(fjac);

			alu.Solve(p, p);

			lnsrch(xold, fold, g, p, x, f, stpmax, check, fmin);

			test = 0.0;
			for (i = 0; i < n; i++)
				if (std::abs(fvec[i]) > test) test = std::abs(fvec[i]);
			if (test < TOLF) {
				check = false;
				return;
			}
			if (check) {
				test = 0.0;
				den = std::max(f, 0.5 * n);
				for (i = 0; i < n; i++) {
					temp = std::abs(g[i]) * std::max<Real>(std::abs(x[i]), 1.0) / den;
					if (temp > test) test = temp;
				}
				check = (test < TOLMIN);
				return;
			}
			test = 0.0;
			for (i = 0; i < n; i++) {
				temp = (std::abs(x[i] - xold[i])) / std::max<Real>(std::abs(x[i]), 1.0);
				if (temp > test) test = temp;
			}
			if (test < TOLX)
				return;
		}
		throw("MAXITS exceeded in newt");
	}

	//Given an initial guess x[0..n - 1] for a root in n dimensions, find the root by Broyden’s
	//method embedded in a globally convergent strategy.The vector of functions to be zeroed,
	//called fvec[0..n - 1] in the routine below, is returned by the user - supplied function or functor
	//vecfunc.The routines NRfdjac and NRfmin from newt are used.The output quantity check
	//is false on a normal return and true if the routine has converged to a local minimum of the
	//function fmin or if Broyden’s method can make no further progress.In this case try restarting
	//from a different initial guess
	template <class T>
	void broydn(Vector<Real>& x, bool& check, T& vecfunc)
	{
		const int MAXITS = 200;
		const Real EPS = std::numeric_limits<Real>::epsilon();
		const Real TOLF = 1.0e-8, TOLX = EPS, STPMX = 100.0, TOLMIN = 1.0e-12;
		//Here MAXITS is the maximum number of iterations; EPS is the machine precision; TOLF
		//is the convergence criterion on function values; TOLX is the convergence criterion on delta x;
		//STPMX is the scaled maximum step length allowed in line searches; and TOLMIN is used to
		//decide whether spurious convergence to a minimum of fmin has occurred.
		bool restrt, skip;
		int i, its, j, n = x.size();
		Real den, f, fold, stpmax, sum, temp, test;
		Vector<Real> fvcold(n), g(n), p(n), s(n), t(n), w(n), xold(n);
		QRDecompositionSolver* qr;
		NRfmin<T> fmin(vecfunc);
		NRfdjac<T> fdjac(vecfunc);
		Vector<Real>& fvec = fmin.fvec;
		f = fmin(x);
		test = 0.0;
		for (i = 0; i < n; i++)
			if (std::abs(fvec[i]) > test) test = std::abs(fvec[i]);
		if (test < 0.01 * TOLF) {
			check = false;
			return;
		}
		for (sum = 0.0, i = 0; i < n; i++) sum += POW2(x[i]);
		stpmax = STPMX * std::max(sqrt(sum), Real(n));
		restrt = true;
		for (its = 1; its <= MAXITS; its++) {
			if (restrt) {
				qr = new QRDecompositionSolver(fdjac(x, fvec));
				if (qr->sing) {
					Matrix<Real> one(n, n, 0.0);
					for (i = 0; i < n; i++) one[i][i] = 1.0;
					delete qr;
					qr = new QRDecompositionSolver(one);
				}
			}
			else {
				for (i = 0; i < n; i++) s[i] = x[i] - xold[i];
				for (i = 0; i < n; i++) {
					for (sum = 0.0, j = i; j < n; j++) sum += qr->r[i][j] * s[j];
					t[i] = sum;
				}
				skip = true;
				for (i = 0; i < n; i++) {
					for (sum = 0.0, j = 0; j < n; j++) sum += qr->qt[j][i] * t[j];
					w[i] = fvec[i] - fvcold[i] - sum;
					if (std::abs(w[i]) >= EPS * (std::abs(fvec[i]) + std::abs(fvcold[i]))) skip = false;
					else w[i] = 0.0;
				}
				if (!skip) {
					qr->qtmult(w, t);
					for (den = 0.0, i = 0; i < n; i++) den += POW2(s[i]);
					for (i = 0; i < n; i++) s[i] /= den;
					qr->update(t, s);
					if (qr->sing) throw("singular update in broydn");
				}
			}
			qr->qtmult(fvec, p);
			for (i = 0; i < n; i++)
				p[i] = -p[i];
			for (i = n - 1; i >= 0; i--) {
				for (sum = 0.0, j = 0; j <= i; j++) sum -= qr->r[j][i] * p[j];
				g[i] = sum;
			}
			for (i = 0; i < n; i++) {
				xold[i] = x[i];
				fvcold[i] = fvec[i];
			}
			fold = f;
			qr->RSolve(p, p);
			Real slope = 0.0;
			for (i = 0; i < n; i++) slope += g[i] * p[i];
			if (slope >= 0.0) {
				restrt = true;
				continue;
			}
			lnsrch(xold, fold, g, p, x, f, stpmax, check, fmin);
			test = 0.0;
			for (i = 0; i < n; i++)
				if (std::abs(fvec[i]) > test) test = std::abs(fvec[i]);
			if (test < TOLF) {
				check = false;
				delete qr;
				return;
			}
			if (check) {
				if (restrt) {
					delete qr;
					return;
				}
				else {
					test = 0.0;
					den = std::max(f, 0.5 * n);
					for (i = 0; i < n; i++) {
						temp = std::abs(g[i]) * std::max<Real>(std::abs(x[i]), 1.0) / den;
						if (temp > test) test = temp;
					}
					if (test < TOLMIN) {
						delete qr;
						return;
					}
					else restrt = true;
				}
			}
			else {
				restrt = false;
				test = 0.0;
				for (i = 0; i < n; i++) {
					temp = (std::abs(x[i] - xold[i])) / std::max<Real>(std::abs(x[i]), 1.0);
					if (temp > test) test = temp;
				}
				if (test < TOLX) {
					delete qr;
					return;
				}
			}
		}
		throw("MAXITS exceeded in broydn");
	}
}

///////////////////////////   ./include/algorithms/Statistics.h   ///////////////////////////


namespace MML
{
	class Statistics
	{
	public:
		static Real Avg(Vector<Real>& data)
		{
			Real outAvg = 0.0;
			int n = data.size();

			for (int j = 0; j < n; j++)
				outAvg += data[j];
			outAvg /= n;

			return outAvg;
		}

		static void AvgVar(Vector<Real>& data, Real& outAvg, Real& outVar)
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

		static void Moments(Vector<Real>& data, Real& ave, Real& adev, Real& sdev, Real& var, Real& skew, Real& curt)
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
///////////////////////////   ./include/algorithms/MatrixAlg.h   ///////////////////////////




namespace MML
{
	namespace MatrixUtils
	{
		static bool IsPositiveDefinite(const Matrix<Real>& mat)
		{
			EigenSolver eigenSolver(mat);

			if (eigenSolver.getNumReal() < mat.RowNum())
				return false;

			for (int i = 0; i < mat.RowNum(); i++)
				if (eigenSolver.getRealEigenvalues()[i] < 0.0)
					return false;

			return true;
		}

		static bool IsPositiveSemiDefinite(const Matrix<Real>& mat)
		{
			EigenSolver eigenSolver(mat);

			if (eigenSolver.getNumReal() < mat.RowNum())
				return false;

			for (int i = 0; i < mat.RowNum(); i++)
				if (eigenSolver.getRealEigenvalues()[i] <= 0.0)
					return false;

			return true;
		}
	}
}
///////////////////////////   ./include/algorithms/FunctionAnalyzers.h   ///////////////////////////






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
			case IntegrationMethod::SIMPSON:
				return IntegrateSimpson(helper, a, b);
			case IntegrationMethod::ROMBERG:
				return IntegrateRomberg(helper, a, b);
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
			case IntegrationMethod::SIMPSON:
				return IntegrateSimpson(helper, a, b);
			case IntegrationMethod::ROMBERG:
				return IntegrateRomberg(helper, a, b);
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
			case IntegrationMethod::SIMPSON:
				return IntegrateSimpson(helper, a, b);
			case IntegrationMethod::ROMBERG:
				return IntegrateRomberg(helper, a, b);
			default:
				return IntegrateTrap(helper, a, b);
			}
		}
	};
}

#endif
