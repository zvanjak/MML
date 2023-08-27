#ifndef MML_SINGLE_HEADER
#define MML_SINGLE_HEADER

#include <exception>
#include <stdexcept>

#include <cmath>
#include <limits>
#include <complex>

#include <string>
#include <vector>

#include <fstream>
#include <iostream>
#include <iomanip>

#include <initializer_list>
#include <algorithm>

#include <functional>///////////////////////////   ./include/MMLBase.h   ///////////////////////////







template<class T>
inline T SQR(const T a) {return ((a)*(a));}

template<class T>
inline T SIGN(const T &a, const T &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const float &a, const double &b)
	{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const double &a, const float &b)
	{return (float)(b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a));}
    
// template<class T>
// inline void SWAP(T &a, T &b)
// 	{T dum=a; a=b; b=dum;}

typedef double               Real;      // default real type
typedef std::complex<double> Complex;   // default complex type


namespace MML
{
    class Defaults
    {
        public:
        static inline const double MatrixEqualityPrecision = 1e-15;
        static inline const double VectorEqualityPrecision = 1e-15;

        static inline const double DerivationDefaultStep = 1e-6;
        
        static inline const int    IntegrateTrapMaxSteps = 30;
        static inline const double IntegrateTrapEPS      = 1.0e-10;

        static inline const int    IntegrateSimpMaxSteps = 30;
        static inline const double IntegrateSimpEPS      = 1.0e-8;

        static inline const int    IntegrateRombMaxSteps = 30;
        static inline const double IntegrateRombEPS      = 1.0e-10;

        static inline const double WorkIntegralPrecision = 1e-05;
        static inline const double LineIntegralPrecision = 1e-05;
    };

	class VectorDimensionError : public std::invalid_argument
	{
		public:
        int _size1, _size2;
		VectorDimensionError(std::string inMessage, int size1, int size2) : std::invalid_argument(inMessage), _size1(size1), _size2(size2)
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

    class IntegrationTooManySteps : public std::domain_error
	{
		public:
		std::string _operation;

		IntegrationTooManySteps(std::string inOperation) : std::domain_error(inOperation)
		{
			_operation = inOperation;
		}
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
        // bool contains(const IInterval &other) const = 0;
        // bool intersects(const IInterval &other) const = 0;

        virtual void GetEquidistantCovering(int numPoints) = 0;
    };

}
///////////////////////////   ./include/utilities/Intervals.h   ///////////////////////////



namespace MML
{
    // TODO - implement intervals properly (and use them in test beds)
    ///////////////////////////////////////////////   Interfaces    ///////////////////////////////////////////
    class BaseInterval : public IInterval
    {
    protected:
        Real _lower, _upper;

        BaseInterval(Real lower, Real upper) : _lower(lower), _upper(upper) { }
    public:
        Real getLowerBound() const { return _lower; }
        Real getUpperBound() const { return _upper; }
        Real getLength()     const { return _upper - _lower; }
        
        bool isContinuous()  const { return true; }

        // bool contains(const IInterval &other) const = 0;
        // bool intersects(const IInterval &other) const = 0;

        void GetEquidistantCovering(int numPoints) { }
    };

    class CompleteRInterval : public BaseInterval
    {
    public:
        CompleteRInterval() : BaseInterval(-std::numeric_limits<double>::max(), std::numeric_limits<double>::max()) { }

        bool contains(Real x) const {
            return true;
        }        
    };
    class OpenInterval : public BaseInterval
    {
        // for equidistant covering
        double _lowerRealDif = 0.0000001;
    public:
        OpenInterval(Real lower, Real upper) : BaseInterval(lower, upper) { }

        bool contains(Real x) const {
            return (x > _lower) && (x < _upper);
        }        
    };

    class OpenClosedInterval : public BaseInterval
    {
    public:
        OpenClosedInterval(Real lower, Real upper) : BaseInterval(lower, upper) { }

        bool contains(Real x) const {
            return (x > _lower) && (x <= _upper);
        }        
    };
    class ClosedInterval : public BaseInterval
    {
    public:
        ClosedInterval(Real lower, Real upper) : BaseInterval(lower, upper) { }
        
        bool contains(Real x) const {
            return (x >= _lower) && (x <= _upper);
        }        
    };   
    
    class ClosedOpenInterval : public BaseInterval
    {
    public:
        ClosedOpenInterval(Real lower, Real upper) : BaseInterval(lower, upper) { }
        
        bool contains(Real x) const {
            return (x >= _lower) && (x < _upper);
        }        
    };    
    class NegInfToOpenInterval : public BaseInterval
    {
    public:
        NegInfToOpenInterval(Real upper) : BaseInterval(-std::numeric_limits<double>::max(), upper) { }

        bool contains(Real x) const {
            return (x > _lower) && (x < _upper);
        }        
    };

    class NegInfToClosedInterval : public BaseInterval
    {
    public:
        NegInfToClosedInterval(Real upper) : BaseInterval(-std::numeric_limits<double>::max(), upper) { }

        bool contains(Real x) const {
            return (x > _lower) && (x <= _upper);
        }        
    };

    class OpenToInfInterval : public BaseInterval
    {
    public:
        OpenToInfInterval(Real lower) : BaseInterval(lower, std::numeric_limits<double>::max()) { }

        bool contains(Real x) const {
            return (x > _lower) && (x < _upper);
        }        
    };
    
    class ClosedToInfInterval : public BaseInterval
    {
    public:
        ClosedToInfInterval(Real lower) : BaseInterval(lower, std::numeric_limits<double>::max()) { }

        bool contains(Real x) const {
            return (x >= _lower) && (x < _upper);
        }        
    };

    class Interval : public IInterval
    {
        std::vector<BaseInterval *> _intervals;
    public:
        Interval(BaseInterval &interval)
        {
            //_intervals.push_back(interval);
        }
        Interval(std::initializer_list<BaseInterval> intervals)
        {
        }

        // TODO - ima initializer list s intervalima, i još provjeri i adjusta
    };

    class IntervalUnion : public IInterval
    {
        IInterval &_left, &_right;
    public:
        IntervalUnion(IInterval &left, IInterval &right) : _left(left), _right(right) 
        { 
            if( _left.getLowerBound() > _right.getLowerBound() )
            {
                IInterval &temp = _left;
                _left = _right;
                _right = temp;
                // std::swap(_left, _right);  // doesn't work???
            }
        }

        Real getLowerBound() const { return _left.getLowerBound(); }
        Real getUpperBound() const { return std::max(_left.getUpperBound(), _right.getUpperBound()); }
        Real getLength() const { 
            // FIXME - implement correctly
            // ako se ne sijeku
            if( _left.getUpperBound() < _right.getLowerBound() )
                return _left.getLength() + _right.getLength();
            else {
                // ili je jedan sadrzan cijeli u drugome

                // ili se standardno sijeku
            }
            return _left.getLength() + _right.getLength(); 

        }
        
        bool isContinuous() const { return false; }

        bool contains(Real x) const {
            return _left.contains(x) || _right.contains(x);
        }
        // bool contains(const IInterval &other) const = 0;
        // bool intersects(const IInterval &other) const = 0;

        void GetEquidistantCovering(int numPoints) { }
    };    
}

///////////////////////////   ./include/interfaces/IAlgebra.h   ///////////////////////////

namespace MML
{
    // group
    template<class _ElemType>
    class Group
    {
        public:
        virtual _ElemType  identity() = 0;
        virtual _ElemType  zero() = 0;
        virtual _ElemType  inverse(const _ElemType& a) = 0;
        virtual _ElemType  op(const _ElemType& a, const _ElemType& b) = 0;
    };

    // kako definirati zahtjeve na VecType? mora imati + i *
    class VectorSpaceElem
    {
        // TODO - define interfaces for vector space elements
        
    };

    template<int N, class _FieldType, class _VecType>
    class VectorSpace
    {
        public:
        virtual _FieldType  identity() = 0;
        virtual _VecType    zero() = 0;
        virtual _VecType    inverse(const _VecType &b) = 0;

        // zbrajanje dva vektora
        // množenje vektora skalarom
    };

    class RealNVectorSpace
    {
        public:
        virtual double  identity() = 0;
        virtual double  zero() = 0;
        virtual double  inverse(const double &b) = 0;

        // zbrajanje dva vektora
        // mnozenje vektora skalarom
    };

    // Hilbert space, kompleksni field!

    template<class _VecSpace>
    class VectorFromVecSpace
    {
        // vraca skalarni produkt definiran u template param
    };
}
///////////////////////////   ./include/core/Algebra.h   ///////////////////////////


namespace MML
{
    // TODO - implement Z6 group
    // TODO - permutation group
    // TODO - vector space + Gram Schmidt
    // TODO - linear transformations
    
    // groups
    class GroupZ : public Group<int>
    {
        public:
        virtual int  identity() { return 0;}
        virtual int  zero() { return 0;}
        virtual int  inverse(const int& a) { return -a;}
        virtual int  op(const int& a, const int& b) { return a + b;}
    };

}
///////////////////////////   ./include/core/Vector.h   ///////////////////////////

namespace MML
{
    template<class _Type>
    class Vector
    {
    private:
        std::vector<_Type> _elems;

    public:
        ///////////////////////          Constructors and destructor       //////////////////////
        Vector() {}
        Vector(size_t n) : _elems(n) {}
        Vector(size_t n, _Type val) : _elems(n, val) {}
        Vector(std::vector<_Type> values) : _elems(values) {}
        Vector(std::initializer_list<_Type> list) : _elems(list) {}
        Vector(int n, _Type *vals) : _elems(n) 
        {
            for(int i=0; i<n; ++i) 
                _elems[i] = vals[i];
        }

        void Resize(int newLen)
        {
            _elems.resize(newLen);
        }

        ////////////////////////            Standard stuff             ////////////////////////
        int  size() const { return (int) _elems.size(); }

        static Vector GetUnitVector(int dimVec, int indUnit)
        {
            // TODO - što ako je indUnit neispravan?
            Vector ret(dimVec);
            ret[indUnit] = _Type{1.0};
            return ret;
        }
        
        bool IsEqual(const Vector &b, _Type eps=Defaults::VectorEqualityPrecision) const
        {
            if (size() != b.size() )
                throw VectorDimensionError("Vector::IsEqual - vectors must be equal size", size(), b.size());

            for( int i=0; i<size(); i++ )
            {
                if( std::abs((*this)[i] - b[i]) > eps )
                    return false;
            }
            return true;
        }
        static bool AreEqual(const Vector &a,const Vector &b, _Type eps=Defaults::VectorEqualityPrecision)
        {
            return a.IsEqual(b, eps);
        }

        ///////////////////////////            Operators             ///////////////////////////
        _Type& operator[](int n)       { return _elems[n]; }
        _Type  operator[](int n) const { return _elems[n]; }

        Vector operator+(const Vector &b ) const
        {
            if (size() != b.size() )
                throw VectorDimensionError("Vector::operator+() - vectors must be equal size", size(), b.size());

            Vector ret(b.size());;
            for(int i=0; i<b.size(); i++)
                ret._elems[i] = (*this)[i] + b._elems[i];
            return ret;
        }
        Vector operator-(const Vector &b ) const
        {
            if (size() != b.size() )
                throw VectorDimensionError("Vector::operator-() - vectors must be equal size", size(), b.size());

            Vector ret(b.size());;
            for(int i=0; i<b.size(); i++)
                ret._elems[i] = (*this)[i] - b._elems[i];
            return ret;
        }
        friend Vector operator*(_Type a, const Vector &b )
        {
            Vector ret(b.size());;
            for(int i=0; i<b.size(); i++)
                ret._elems[i] = a * b._elems[i];
            return ret;
        }
        friend Vector operator*(const Vector &a, _Type b )
        {
            Vector ret(a.size());;
            for(int i=0; i<a.size(); i++)
                ret._elems[i] = b * a._elems[i];
            return ret;
        }
        friend Vector operator/(const Vector &a, _Type b)
        {
            Vector ret(a.size());
            for(int i=0; i<a.size(); i++)
                ret._elems[i] = a._elems[i] / b;
            return ret;
        }
        
        //////////////////////                 Operations                 ///////////////////////
        _Type ScalarProductCartesian(Vector &b)
        {
            if (size() != b.size() )
                throw VectorDimensionError("Vector::ScalarProductCartesian - vectors must be equal size", size(), b.size());

            _Type product = 0.0;
            for( int i=0; i<size(); i++ )
                product += (*this)[i] * b[i];
            return product;
        }
        _Type NormL2() const
        {
            _Type norm = 0.0;
            for( int i=0; i<size(); i++ )
                norm += (*this)[i] * (*this)[i];
            return std::sqrt(norm);
        }
        _Type AngleToVector(Vector &b)
        {
            if (size() != b.size() )
                throw VectorDimensionError("Vector::AngleToVector - vectors must be equal size", size(), b.size());

            _Type cosAngle = ScalarProductCartesian(b) / (NormL2() * b.NormL2());
            return std::acos(cosAngle);
        }

        ///////////////////////////               I/O                 ///////////////////////////
        std::string to_string(int width, int precision)
        {
            std::stringstream str;

            Print(str, width, precision);

            return str.str();
        }
        void Print(std::ostream& stream, int width, int precision) const
        {
            stream << "[";
            bool first = true;
            for(const _Type& x : _elems)
            {
                if( !first )
                    stream << ", ";
                else
                    first = false;

                stream << std::setw(width) << std::setprecision(precision) << x;
            }
            stream << "]";
        }   
        friend std::ostream& operator<<(std::ostream& stream, const Vector &a)
        {
            a.Print(stream, 15, 10);

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

///////////////////////////   ./include/core/VectorN.h   ///////////////////////////

namespace MML
{
    template<class _Type, int N> 
    class VectorN
    {
    protected:
        _Type  _val[N] = {0};

    public:
        ///////////////////////          Constructors and destructor       //////////////////////
        VectorN() {}
        VectorN(const _Type &init_val) {
            for(int i=0; i<N; ++i) 
                _val[i] = init_val;
        }
        VectorN(std::initializer_list<_Type> list) 
        { 
            int count{ 0 };
            for (auto element : list)
            {
                _val[count] = element;
                ++count;

                if( count >= N )
                    break;
            }
        }
        VectorN(std::vector<_Type> list) 
        { 
            int count{ 0 };
            for (auto element : list)
            {
                _val[count] = element;
                ++count;

                if( count >= N )
                    break;
            }
        }        
        VectorN(_Type *vals)
        {
            for(int i=0; i<N; ++i) 
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

        bool IsEqual(VectorN &b, _Type eps) const
        {
            for( int i=0; i<N; i++ )
            {
                if( fabs((*this)[i] - b[i]) > eps )
                    return false;
            }
            return true;
        }
        static bool AreEqual(const VectorN &a, const VectorN &b, _Type eps=Defaults::VectorEqualityPrecision)
        {
            return a.IsEqual(b, eps);
        }
        
        ///////////////////////////            Operators             ///////////////////////////
        _Type& operator[](int n)       { return _val[n]; }
        _Type  operator[](int n) const { return _val[n]; }

        VectorN operator+(const VectorN &b ) const
        {
            VectorN ret;
            for(int i=0; i<N; i++)
                ret._val[i] = _val[i] + b._val[i];
            return ret;
        }
        VectorN operator-(const VectorN &b ) const
        {
            VectorN ret;
            for(int i=0; i<N; i++)
                ret._val[i] = _val[i] - b._val[i];
            return ret;
        }
        friend VectorN operator*(const VectorN &a, _Type b )
        {
            VectorN ret;
            for(int i=0; i<N; i++)
                ret._val[i] = a[i] * b;
            return ret;
        }
        friend VectorN operator*(_Type a, const VectorN &b )
        {
            VectorN ret;
            for(int i=0; i<N; i++)
                ret._val[i] = a * b[i];
            return ret;
        }
        friend VectorN operator/(const VectorN &a, _Type b)
        {
            VectorN ret;
            for(int i=0; i<N; i++)
                ret._val[i] = a[i] / b;
            return ret;
        }

        //////////////////////                 Operations                 ///////////////////////
        _Type ScalarProductCartesian(const VectorN &b) const
        {
            _Type product = 0.0;
            for( int i=0; i<N; i++ )
                product += (*this)[i] * b[i];
            return product;
        }
        _Type NormL2() const
        {
            _Type norm = 0.0;
            for( int i=0; i<size(); i++ )
                norm += (*this)[i] * (*this)[i];
            return std::sqrt(norm);
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
            for(const _Type& x : _val)
            {
                if( !first )
                    stream << ", ";
                else
                    first = false;

                stream << std::setw(width) << std::setprecision(precision) << x;
            }
            stream << "]";

            return stream;
        }   
        friend std::ostream& operator<<(std::ostream& stream, const VectorN<_Type, N> &a)
        {
            a.Print(stream,15,10);

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

///////////////////////////   ./include/core/Matrix.h   ///////////////////////////
//#include <format>



namespace MML
{
    template<class _Type>
    class Matrix
    {
    private:
        int  _rows;
        int  _cols;
        _Type **_ptrData;

        void Init(int rows, int cols)
        {
            if( rows <= 0 || cols < 0 )
                throw MatrixDimensionError("Matrix::Invert - rowNum and colNum must be positive", rows, cols, -1, -1);

            _rows = rows;
            _cols = cols;
            // _ptrData = new _Type *[_rows];
            // for (int i = 0; i < _rows; ++i)
            //     _ptrData[i] = new _Type[_cols];

            _ptrData = new _Type*[rows];

            int numElem = rows * cols;
            if (_ptrData) 
                _ptrData[0] = numElem>0 ? new _Type[numElem] : nullptr;
            
            for (int i=1; i<rows; i++) 
                _ptrData[i] = _ptrData[i-1] + cols;
        }

    public:
        ///////////////////////          Constructors and destructor       //////////////////////
        Matrix() : _rows(0), _cols(0),  _ptrData{nullptr} {} 
        Matrix(int rows, int cols) : _rows(rows), _cols(cols)
        {
            Init(rows, cols);
            for (int i = 0; i < _rows; ++i)
                for (int j = 0; j < _cols; ++j)
                    _ptrData[i][j] = 0;
        }
        Matrix(int rows, int cols, _Type val) : _rows(rows), _cols(cols)
        {
            Init(rows, cols);
            for (int i = 0; i < _rows; ++i)
                for (int j = 0; j < _cols; ++j)
                    _ptrData[i][j] = val;
        }        
        Matrix(int rows, int cols, _Type *val) : _rows(rows), _cols(cols)
        {
            Init(rows, cols);
            for (int i = 0; i < _rows; ++i)
                for (int j = 0; j < _cols; ++j)
                    _ptrData[i][j] = *val++;
        }   
        Matrix(int rows, int cols, std::initializer_list<_Type> values) : _rows(rows), _cols(cols)
        {
            Init(rows, cols);
            
            auto val = values.begin();
            for (int i = 0; i < _rows; ++i)
                for (int j = 0; j < _cols; ++j)
                    if( val != values.end() )
                        _ptrData[i][j] = *val++;
                    else
                        _ptrData[i][j] = 0.0;
        }
        Matrix(const Matrix &m) : _rows(m._rows), _cols(m._cols)
        {
            Init(m._rows, m._cols);

            for (int i = 0; i < _rows; ++i)
                for (int j = 0; j < _cols; ++j)
                    _ptrData[i][j] = m._ptrData[i][j];
        }
        Matrix(Matrix &&m)
        {
            _ptrData = m._ptrData;

            _rows = m._rows;
            _cols = m._cols;

            m._rows = 0;
            m._cols = 0;
            m._ptrData = nullptr;
        }
        ~Matrix()
        {
            if (_ptrData != NULL) {
                delete[] (_ptrData[0]);
                delete[] (_ptrData);
            }
        }

        void Resize(int rows, int cols) 
        {
            for (int i = 0; i < _rows; ++i)
                if (_ptrData != nullptr && _ptrData[i] != nullptr)
                    delete[] _ptrData[i];
            if (_ptrData != nullptr)
                delete[] _ptrData;

            Init(rows, cols);

            for (int i = 0; i < _rows; ++i)
                for (int j = 0; j < _cols; ++j)
                    _ptrData[i][j] = 0;
        }

        typedef _Type value_type;      // make T available externally

        ////////////////////////            Standard stuff             ////////////////////////
        int RowNum() const { return (int) _rows; }
        int ColNum() const { return (int) _cols; }

        static Matrix GetUnitMatrix(int dim)
        {
            Matrix unitMat(dim, dim);
            
            for( int i=0; i<dim; i++ )
                unitMat._ptrData[i][i] = 1.0;

            return unitMat;            
        }
        void   MakeUnitMatrix(void)
        {
            if (_rows == _cols)
            {
                for (int i = 0; i < _rows; i++)
                    for (int j = 0; j < _cols; j++)
                        if (i == j)
                            _ptrData[i][j] = 1;
                        else
                            _ptrData[i][j] = 0;
            }
            else
                throw MatrixDimensionError("MakeUnitMatrix - must be square matrix", _rows, _cols, -1, -1);
        }

        bool IsEqual(const Matrix &b, _Type eps=Defaults::MatrixEqualityPrecision) const
        {
            if( RowNum() != b.RowNum() || ColNum() != b.ColNum() )
                return false;

            for( int i=0; i<RowNum(); i++ )
                for( int j=0; j<ColNum(); j++ )
                {
                    if( std::abs(_ptrData[i][j] - b._ptrData[i][j]) > eps )
                        return false;
                }
                
            return true;
        }
        static bool AreEqual(const Matrix &a, const Matrix &b, _Type eps=Defaults::MatrixEqualityPrecision) 
        {
            return a.IsEqual(b, eps);
        }        

        /////////////////////          Vector-Matrix conversion           /////////////////////
        static Matrix RowMatrixFromVector(const Vector<_Type> &b)
        {
            Matrix ret(1, (int) b.size());
            for( int i=0; i<b.size(); i++)
                ret[0][i] = b[i];

            return ret;
        }
        static Matrix ColumnMatrixFromVector(const Vector<_Type> &b)
        {
            Matrix ret((int) b.size(), 1);
            for( int i=0; i<b.size(); i++)
                ret[i][0] = b[i];

            return ret;
        }
        static Vector<_Type> VectorFromRow(const Matrix &a, int rowInd)
        {
           if( rowInd >= a.RowNum() )
                throw MatrixAccessBoundsError("VectorFromRow - row index must be less then a.RowNum()", rowInd, 0, a.RowNum(), a.ColNum());

            Vector<_Type> ret(a.ColNum());
            for( int i=0; i<a.ColNum(); i++)
                ret[i] = a(rowInd,i);

            return ret;
        }
        static Vector<_Type> VectorFromColumn(const Matrix &a, int colInd)
        {
           if( colInd >= a.ColNum() )
                throw MatrixAccessBoundsError("VectorFromColumn - column index must be less then a.ColNum()", 0, colInd, a.RowNum(), a.ColNum());

            Vector<_Type> ret(a.RowNum());
            for( int i=0; i<a.RowNum(); i++)
                ret[i] = a(i,colInd);

            return ret;
        }
        static Vector<_Type> VectorFromDiagonal(const Matrix &a)
        {
           if( a.RowNum() != a.ColNum() ) 
                throw MatrixDimensionError("VectorFromDiagonal - must be square matrix", a.RowNum(), a.ColNum(), -1, -1);

            Vector<_Type> ret(a.RowNum());
            for( int i=0; i<a.RowNum(); i++)
                ret[i] = a(i,i);

            return ret;
        }

        ///////////////////////////            Operators             ///////////////////////////
        Matrix &operator=(const Matrix &m)
        {
            if (this == &m)
                return *this;

            if (_rows != m._rows || _cols != m._cols)
            {
                for (size_t i = 0; i < _rows; ++i)
                    delete[] _ptrData[i];
                delete[] _ptrData;

                Init(m._rows, m._cols);
            }

            for (size_t i = 0; i < _rows; ++i)
                for (size_t j = 0; j < _cols; ++j)
                    _ptrData[i][j] = m._ptrData[i][j];

            return *this;
        }
        Matrix &operator=(Matrix &&m)
        {
            if (this == &m)
                return *this;

            std::swap(_ptrData, m._ptrData);
            std::swap(_rows, m._rows);
            std::swap(_cols, m._cols);

            return *this;
        }

        _Type* operator[](int i)                   { return _ptrData[i]; }
        const _Type* operator[](const int i) const { return _ptrData[i]; }
        
        _Type  operator()(int i, int j) const { return _ptrData[i][j]; }
        _Type& operator()(int i, int j)       { return _ptrData[i][j]; }        

        // version with checking bounds
        _Type  ElemAt(int i, int j) const 
        { 
            if( i<0 || i>=RowNum() || j<0 || j>=ColNum() )
                throw MatrixAccessBoundsError("Matrix::ElemAt", i, j, RowNum(), ColNum());

            return _ptrData[i][j]; 
        }
        _Type& ElemAt(int i, int j)       
        {
            if( i<0 || i>=RowNum() || j<0 || j>=ColNum() )
                throw MatrixAccessBoundsError("Matrix::ElemAt", i, j, RowNum(), ColNum());

            return _ptrData[i][j]; 
        }

        Matrix operator+(const Matrix &other) const
        {
            if (_rows != other._rows || _cols != other._cols)
                throw MatrixDimensionError("Matrix::operator+() - must be same dim", _rows, _cols, other._rows, other._cols);

            Matrix temp(_rows, _cols);
            for (size_t i = 0; i < _rows; i++)
                for (size_t j = 0; j < _cols; j++)
                    temp._ptrData[i][j] = other._ptrData[i][j] + _ptrData[i][j];

            return temp;
        }
        Matrix operator-(const Matrix &other) const
        {
            if (_rows != other._rows || _cols != other._cols)
                throw MatrixDimensionError("Matrix::operator-() - must be same dim", _rows, _cols, other._rows, other._cols);

            Matrix temp(_rows, _cols);
            for (int i = 0; i < _rows; i++)
                for (int j = 0; j < _cols; j++)
                    temp._ptrData[i][j] =  _ptrData[i][j] - other._ptrData[i][j];

            return temp;
        }
        Matrix operator*( const Matrix &b ) const
        {
            if( ColNum() != b.RowNum() )
                throw MatrixDimensionError("Matrix::operator*() - a.colNum must be equal to b.rowNum", _rows, _cols, b._rows, b._cols);

            Matrix	ret(RowNum(), b.ColNum());
            for( int i=0; i<ret.RowNum(); i++ )
                for( int j=0; j<ret.ColNum(); j++ )
                {
                    ret._ptrData[i][j] = 0;
                    for(int k=0; k<ColNum(); k++ )
                        ret._ptrData[i][j] += _ptrData[i][k] * b._ptrData[k][j];
                }

            return	ret;
        }

        friend Matrix operator*( const Matrix &a, _Type b )
        {
            int	i, j;
            Matrix	ret(a.RowNum(), a.ColNum());

            for( i=0; i<a.RowNum(); i++ )
                for( j=0; j<a.ColNum(); j++ )
                    ret[i][j] = a._ptrData[i][j] * b;

            return ret;
        }
        friend Matrix operator*( _Type a, const Matrix &b )
        {
            int	i, j;
            Matrix	ret(b.RowNum(), b.ColNum());

            for( i=0; i<b.RowNum(); i++ )
                for( j=0; j<b.ColNum(); j++ )
                    ret[i][j] = a * b._ptrData[i][j];

            return ret;
        }
        friend Matrix operator/(const Matrix &a, _Type b )
        {
            int	i, j;
            Matrix	ret(a.RowNum(), a.ColNum());

            for( i=0; i<a.RowNum(); i++ )
                for( j=0; j<a.ColNum(); j++ )
                    ret[i][j] = a._ptrData[i][j] / b;

            return ret;
        }

        friend Vector<_Type> operator*( const Matrix &a, const Vector<_Type> &b )
        {
            if( a.ColNum() != b.size() )
                throw MatrixDimensionError("operator*(Mat a, Vec b) - a.colNum must be equal to vector size", a._rows, a._cols, (int) b.size(), -1);

            Vector<_Type>	ret(a.RowNum());
            for( int i=0; i<a.RowNum(); i++ )
            {
                ret[i] = 0;
                for( int j=0; j<a.ColNum(); j++ )
                    ret[i] += a._ptrData[i][j] * b[j];
            }

            return ret;
        }
        friend Vector<_Type> operator*( const Vector<_Type> &a, const Matrix &b )
        {
            if( a.size() != b.RowNum() )
            {
                //std::string error = std::format("Hello {}!\n", "world");
                throw MatrixDimensionError("operator*(Vec a, Mat b) - vector size must be equal to b.rowNum", (int) a.size(), -1, b._rows, b._cols);
            }

            Vector<_Type>	ret(b.ColNum());
            for( int i=0; i<b.ColNum(); i++ )
            {
                ret[i] = 0;
                for( int j=0; j<b.RowNum(); j++ )
                    ret[i] += a[j] * b(j, i);
            }

            return ret;
        }

        //////////////////////            Inverse & Transpose             ///////////////////////
        void   Invert()
        {
            if( RowNum() != ColNum() ) 
                throw MatrixDimensionError("Matrix::Invert - must be square matrix", _rows, _cols, -1, -1);

            Matrix& a = *this;            
            Matrix  b(RowNum(), 1);      // dummy rhs

            b(0,0) = 1.0;

            int i,icol,irow,j,k,l,ll;
            _Type big,dum,pivinv;

            int n=a.RowNum();
            int m=b.ColNum();
            std::vector<int> indxc(n),indxr(n),ipiv(n);
            for (j=0;j<n;j++) ipiv[j]=0;
            for (i=0;i<n;i++) {
                big=0.0;
                for (j=0;j<n;j++)
                    if (ipiv[j] != 1)
                        for (k=0;k<n;k++) {
                            if (ipiv[k] == 0) {
                                if (std::abs(a[j][k]) >= big) {
                                    big=std::abs(a[j][k]);
                                    irow=j;
                                    icol=k;
                                }
                            }
                        }
                ++(ipiv[icol]);
                if (irow != icol) {
                    for (l=0;l<n;l++) std::swap(a[irow][l],a[icol][l]);
                    for (l=0;l<m;l++) std::swap(b[irow][l],b[icol][l]);
                }
                indxr[i]=irow;
                indxc[i]=icol;

                if (a[icol][icol] == 0.0) 
                    throw SingularMatrixError("Matrix::Invert, gaussj: Singular Matrix");

                pivinv=1.0/a[icol][icol];
                a[icol][icol]=1.0;
                for (l=0;l<n;l++) a[icol][l] *= pivinv;
                for (l=0;l<m;l++) b[icol][l] *= pivinv;
                for (ll=0;ll<n;ll++)
                    if (ll != icol) {
                        dum=a[ll][icol];
                        a[ll][icol]=0.0;
                        for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
                        for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
                    }
            }
            for (l=n-1;l>=0;l--) {
                if (indxr[l] != indxc[l])
                    for (k=0;k<n;k++)
                        std::swap(a[k][indxr[l]],a[k][indxc[l]]);
            }
        }   
        Matrix GetInverse() const
        {
            if( RowNum() != ColNum() ) 
                throw MatrixDimensionError("Matrix::GetInverse - must be square matrix", _rows, _cols, -1, -1);

            Matrix a(*this);              // making a copy, where inverse will be stored at the end
            
            a.Invert();
            
            return a;
        }     

        void   Transpose()
        {
            if( RowNum() != ColNum() ) 
                throw MatrixDimensionError("Matrix::Transpose - inplace Transpose possible only for square  matrix", _rows, _cols, -1, -1);


            for (size_t i = 0; i < RowNum(); i++)
                for (size_t j = i+1; j < ColNum(); j++)
                    std::swap(_ptrData[i][j], _ptrData[j][i]);
        }
        Matrix GetTranspose() const
        {
            Matrix ret(ColNum(), RowNum());

            for (size_t i = 0; i < RowNum(); i++)
                for (size_t j = 0; j < ColNum(); j++)
                    ret((int) i, (int) j) = _ptrData[j][i];

            return ret;
        }

        ///////////////////////////               I/O                 ///////////////////////////
        void   Print(std::ostream& stream, int width, int precision) const
        {
            stream << "Rows: " << RowNum() << " Cols: " << ColNum() << std::endl;

            for (size_t i = 0; i < RowNum(); i++)
            {
                stream << "[ ";
                for (size_t j = 0; j < ColNum(); j++)
                {
                    stream << std::setw(width) << std::setprecision(precision) << _ptrData[i][j] << ", ";
                }                
                stream << " ]" << std::endl;
            }
        }
        friend std::ostream& operator<<(std::ostream& stream, const Matrix &a)
        {
            a.Print(stream, 10, 3);

            return stream;
        }   
    
        // TODO - load matrix from file
        static Matrix LoadFromFile(std::string inFileName)
        {
            Matrix ret;
            return ret;
        }
    };

    typedef Matrix<int>     MatrixInt;
    typedef Matrix<float>   MatrixFlt;
    typedef Matrix<double>  MatrixDbl;
    typedef Matrix<Complex> MatrixComplex;

    typedef Matrix<int>     MatI;
    typedef Matrix<float>   MatF;
    typedef Matrix<double>  MatD;
    typedef Matrix<Complex> MatC;    
}

///////////////////////////   ./include/core/Matrix3D.h   ///////////////////////////

namespace MML
{
    template <class _Type>
    class Matrix3D {
    private:
        int _n;
        int _m;
        int _k;
        _Type ***v;

    public:
        Matrix3D(): _n(0), _m(0), _k(0), v(nullptr) {}

        Matrix3D(int n, int m, int k) : _n(n), _m(m), _k(k), v(new _Type**[n])
        {
            int i,j;
            
            v[0]    = new _Type*[n*m];
            v[0][0] = new _Type[n*m*k];
            
            for(j=1; j<m; j++) 
                v[0][j] = v[0][j-1] + k;

            for(i=1; i<n; i++) {
                v[i]    = v[i-1] + m;
                v[i][0] = v[i-1][0] + m*k;
                
                for(j=1; j<m; j++) 
                    v[i][j] = v[i][j-1] + k;
            }
        }        

        ~Matrix3D()
        {
            if (v != NULL) {
                delete[] (v[0][0]);
                delete[] (v[0]);
                delete[] (v);
            }
        }

        //subscripting: pointer to row i
        inline _Type**              operator[](const int i)       { return v[i]; }
        inline const _Type* const * operator[](const int i) const { return v[i]; }

        inline int dim1() const { return _n; }
        inline int dim2() const { return _m; }
        inline int dim3() const { return _k; }
    };
}

///////////////////////////   ./include/core/MatrixNM.h   ///////////////////////////


namespace MML
{
    template <class _Type, int N, int M>
    class MatrixNM
    {
    public:
        _Type _vals[N][M] = {{0}};

    public:
        //////////////////////////             Constructors           /////////////////////////
        MatrixNM() {}
        MatrixNM(std::initializer_list<_Type> values) 
        {
            auto val = values.begin();
            for (size_t i = 0; i < RowNum(); ++i)
                for (size_t j = 0; j < ColNum(); ++j)
                    if( val != values.end() )
                        _vals[i][j] = *val++;
                    else
                        _vals[i][j] = 0.0;
        }
        MatrixNM(const MatrixNM &m) 
        {
            for (size_t i = 0; i < RowNum(); ++i)
                for (size_t j = 0; j < ColNum(); ++j)
                    _vals[i][j] = m._vals[i][j];
        }
        MatrixNM(const _Type &m)        // initialize as diagonal matrix
        {
            for( int i=0; i<N; i++ )
                _vals[i][i] = _Type{m};
        }        

        typedef _Type value_type;      // make T available externally

        ////////////////////////            Standard stuff             ////////////////////////
        int RowNum() const { return N; }
        int ColNum() const { return M; }

        static MatrixNM GetUnitMatrix() 
        {
            MatrixNM unitMat;
            
            for( int i=0; i<N; i++ )
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
        }
        
        /////////////////////          Vector-Matrix conversion           /////////////////////
        static MatrixNM<_Type, 1, M>  RowMatrixFromVector(const VectorN<_Type, M> &b)
        {
            MatrixNM<_Type,1,M>  ret;
            for( int j=0; j<M; j++)
                ret._vals[0][j] = b[j];

            return ret;
        }      
        static MatrixNM<_Type, N, 1>  ColumnMatrixFromVector(const VectorN<_Type, N> &b)
        {
            MatrixNM<_Type,N, 1>  ret;
            for( int i=0; i<M; i++)
                ret._vals[i][0] = b[i];

            return ret;
        }
        static VectorN<_Type, M> VectorFromRow(const MatrixNM<_Type, N,M> &a, int rowInd)
        {
            VectorN<_Type, M> ret;
            for( int j=0; j<M; j++)
                ret[j] = a._vals[rowInd][j];

            return ret;
        }
        static VectorN<_Type, N> VectorFromColumn(const MatrixNM<_Type, N,M> &a, int colInd)
        {
            VectorN<_Type, N> ret;
            for( int i=0; i<N; i++)
                ret[i] = a._vals[i][colInd];

            return ret;
        }
        static VectorN<_Type, N> VectorFromDiagonal(const MatrixNM<_Type, N,N> &a)
        {
            VectorN<_Type, N> ret;
            for( int i=0; i<N; i++)
                ret[i] = a._vals[i][i];

            return ret;
        }        

        bool IsEqual(const MatrixNM &b, _Type eps=Defaults::MatrixEqualityPrecision) const
        {
            for( int i=0; i<RowNum(); i++ )
                for( int j=0; j<ColNum(); j++ )
                    if( std::abs(_vals[i][j] - b._vals[i][j]) > eps )
                        return false;
                
            return true;
        }
        bool AreEqual(const MatrixNM &a, const MatrixNM &b, _Type eps=Defaults::MatrixEqualityPrecision) const
        {
            return a.IsEqual(b, eps);
        }   
        
        ///////////////////////////            Operators             ///////////////////////////
        MatrixNM &operator=(const MatrixNM &m)
        {
            if (this == &m)
                return *this;

            for (size_t i = 0; i < RowNum(); ++i)
                for (size_t j = 0; j < ColNum(); ++j)
                    _vals[i][j] = m._vals[i][j];

            return *this;
        }
        MatrixNM &operator=(const _Type &m)
        {
            if (this == &m)
                return *this;

            for (size_t i = 0; i < RowNum(); ++i)
                for (size_t j = 0; j < ColNum(); ++j)
                    _vals[i][j] = m;

            return *this;
        }        

        _Type  operator()(int i, int j) const { return _vals[i][j]; }
        _Type& operator()(int i, int j)       { return _vals[i][j]; }        

        // version with checking bounds
        _Type  ElemAt(int i, int j) const 
        { 
            if( i<0 || i>=RowNum() || j<0 || j>=ColNum() )
                throw MatrixAccessBoundsError("MatrixNM::ElemAt", i, j, RowNum(), ColNum());

            return _vals[i][j]; 
        }
        _Type& ElemAt(int i, int j)       
        {
            if( i<0 || i>=RowNum() || j<0 || j>=ColNum() )
                throw MatrixAccessBoundsError("MatrixNM::ElemAt", i, j, RowNum(), ColNum());

            return _vals[i][j]; 
        }

        MatrixNM operator+(const MatrixNM &other) const
        {
            MatrixNM temp;
            for (size_t i = 0; i < RowNum(); i++)
                for (size_t j = 0; j < ColNum(); j++)
                    temp._vals[i][j] = other._vals[i][j] + _vals[i][j];
            return temp;
        }
        MatrixNM operator-(const MatrixNM &other) const
        {
            MatrixNM temp;
            for (size_t i = 0; i < RowNum(); i++)
                for (size_t j = 0; j < ColNum(); j++)
                    temp._vals[i][j] = _vals[i][j] - other._vals[i][j];
            return temp;
        }
        template<int K>
        MatrixNM<_Type, N, K>  operator*( const MatrixNM<_Type, M, K> &b ) const
        {
            MatrixNM<_Type, N, K>	ret;

            for( int i=0; i<ret.RowNum(); i++ )
                for( int j=0; j<ret.ColNum(); j++ )
                {
                    ret._vals[i][j] = 0;
                    for(int k=0; k<ColNum(); k++ )
                        ret._vals[i][j] += _vals[i][k] * b._vals[k][j];
                }

            return	ret;
        }        

        friend MatrixNM operator*( const MatrixNM &a, _Type b )
        {
            int	i, j;
            MatrixNM	ret(a);

            for( i=0; i<a.RowNum(); i++ )
                for( j=0; j<a.ColNum(); j++ )
                    ret._vals[i][j] *= b;

            return ret;
        }
        friend MatrixNM operator/( const MatrixNM &a, _Type b )
        {
            int	i, j;
            MatrixNM	ret(a);

            for( i=0; i<a.RowNum(); i++ )
                for( j=0; j<a.ColNum(); j++ )
                    ret._vals[i][j] /= b;

            return ret;
        }
        friend MatrixNM operator*( _Type a, const MatrixNM &b )
        {
            int	i, j;
            MatrixNM	ret;

            for( i=0; i<b.RowNum(); i++ )
                for( j=0; j<b.ColNum(); j++ )
                    ret._vals[i][j] = a * b._vals[i][j];

            return ret;
        }

        friend VectorN<_Type, N> operator*( const MatrixNM<_Type, N, M> &a, const VectorN<_Type, M> &b )
        {
            int	i, j;
            VectorN<_Type, N>	ret;

            for( i=0; i<a.RowNum(); i++ )
            {
                ret[i] = 0;
                for( j=0; j<a.ColNum(); j++ )
                    ret[i] += a._vals[i][j] * b[j];
            }

            return ret;
        }
        friend VectorN<_Type, M> operator*( const VectorN<_Type, N> &a, const MatrixNM<_Type, N, M> &b )
        {
            int	i, j;
            VectorN<_Type, M>	ret;

            for( i=0; i<b.ColNum(); i++ )
            {
                ret[i] = 0;
                for( j=0; j<b.RowNum(); j++ )
                    ret[i] += a[j] * b._vals[j][i];
            }

            return ret;
        }

        //////////////////////            Inverse & Transpose             ///////////////////////
        void Invert()
        {
            if( RowNum() != ColNum() ) 
                throw MatrixDimensionError("Matrix::Invert - must be square matrix", N, M, -1, -1);

            MatrixNM& a = *this;            
            MatrixNM<_Type,N,1>  b;      // dummy rhs

            b(0,0) = 1.0;

            int i,icol,irow,j,k,l,ll;
            _Type big,dum,pivinv;

            int n=RowNum();
            int m=b.ColNum();
            std::vector<int> indxc(n),indxr(n),ipiv(n);
            for (j=0;j<n;j++) ipiv[j]=0;
            for (i=0;i<n;i++) {
                big=0.0;
                for (j=0;j<n;j++)
                    if (ipiv[j] != 1)
                        for (k=0;k<n;k++) {
                            if (ipiv[k] == 0) {
                                if (std::abs(a._vals[j][k]) >= big) {
                                    big=std::abs(a._vals[j][k]);
                                    irow=j;
                                    icol=k;
                                }
                            }
                        }
                ++(ipiv[icol]);
                if (irow != icol) {
                    for (l=0;l<n;l++) std::swap(a._vals[irow][l],a._vals[icol][l]);
                    for (l=0;l<m;l++) std::swap(b._vals[irow][l],b._vals[icol][l]);
                }
                indxr[i]=irow;
                indxc[i]=icol;

                if (a._vals[icol][icol] == 0.0) 
                    throw SingularMatrixError("MatrixNM::Invert, gaussj: Singular Matrix");

                pivinv=1.0/a._vals[icol][icol];
                a._vals[icol][icol]=1.0;
                for (l=0;l<n;l++) a._vals[icol][l] *= pivinv;
                for (l=0;l<m;l++) b._vals[icol][l] *= pivinv;
                for (ll=0;ll<n;ll++)
                    if (ll != icol) {
                        dum=a._vals[ll][icol];
                        a._vals[ll][icol]=0.0;
                        for (l=0;l<n;l++) a._vals[ll][l] -= a._vals[icol][l]*dum;
                        for (l=0;l<m;l++) b._vals[ll][l] -= b._vals[icol][l]*dum;
                    }
            }
            for (l=n-1;l>=0;l--) {
                if (indxr[l] != indxc[l])
                    for (k=0;k<n;k++)
                        std::swap(a._vals[k][indxr[l]],a._vals[k][indxc[l]]);
            }
        }  
        MatrixNM GetInverse() const
        {
            if( RowNum() != ColNum() ) 
                throw MatrixDimensionError("Matrix::GetInverse - must be square matrix", N, M, -1, -1);

            MatrixNM a(*this);              // making a copy, where inverse will be stored at the end
            
            a.Invert();
            
            return a;
        }     

        void Transpose()
        {
            if( RowNum() != ColNum() ) 
                throw MatrixDimensionError("Matrix::Transpose - inplace Transpose possible only for square  matrix", N, M, -1, -1);

            for (size_t i = 0; i < RowNum(); i++)
                for (size_t j = i+1; j < ColNum(); j++)
                    std::swap(_vals[i][j], _vals[j][i]);
        }
        MatrixNM<_Type, M,N> GetTranspose() const
        {
            MatrixNM<_Type,M,N> ret;

            for (size_t i = 0; i < RowNum(); i++)
                for (size_t j = 0; j < ColNum(); j++)
                    ret._vals[i][j] = _vals[j][i];

            return ret;
        }             
 
        ///////////////////////////               I/O                 ///////////////////////////
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
        friend std::ostream& operator<<(std::ostream& stream, const MatrixNM &a)
        {
            a.Print(stream, 10, 3);

            return stream;
        }     
    };

    typedef MatrixNM<float, 2, 2> Matrix22Flt;
    typedef MatrixNM<float, 3, 3> Matrix33Flt;
    typedef MatrixNM<float, 4, 4> Matrix44Flt;

    typedef MatrixNM<double, 2, 2> Matrix22Dbl;
    typedef MatrixNM<double, 3, 3> Matrix33Dbl;
    typedef MatrixNM<double, 4, 4> Matrix44Dbl;

    typedef MatrixNM<Complex, 2, 2> Matrix22Complex;
    typedef MatrixNM<Complex, 3, 3> Matrix33Complex;
    typedef MatrixNM<Complex, 4, 4> Matrix44Complex;

    typedef MatrixNM<float, 2, 2> Mat22F;
    typedef MatrixNM<float, 3, 3> Mat33F;
    typedef MatrixNM<float, 4, 4> Mat44F;

    typedef MatrixNM<double, 2, 2> Mat22D;
    typedef MatrixNM<double, 3, 3> Mat33D;
    typedef MatrixNM<double, 4, 4> Mat44D;

    typedef MatrixNM<Complex, 2, 2> Mat22C;
    typedef MatrixNM<Complex, 3, 3> Mat33C;
    typedef MatrixNM<Complex, 4, 4> Mat44C;
}

///////////////////////////   ./include/core/MatrixSparse.h   ///////////////////////////

namespace MML
{

}

///////////////////////////   ./include/interfaces/ITensor.h   ///////////////////////////

namespace MML
{
    template<int N>
    class ITensor2
    {
    public:
        virtual double  Component(int i, int j) const = 0;
        virtual double& Component(int i, int j) = 0;
    };

    template<int N>
    class ITensorField2
    {
    public:
        virtual double Component(int i, int j, const VectorN<Real, N> &pos) const = 0;
        virtual void   ValueAtPoint(const VectorN<Real, N> &pos, ITensor2<N> &outVal) const = 0;
    };

    template<int N>
    class ITensor3
    {
    public:
        virtual double  Component(int i, int j, int k) const = 0;
        virtual double& Component(int i, int j, int k) = 0;
    };

    template<int N>
    class ITensorField3
    {
    public:
        virtual double Component(int i, int j, int k, const VectorN<Real, N> &pos) const = 0;
        virtual void   ValueAtPoint(const VectorN<Real, N> &pos, ITensor3<N> &outVal) const = 0;
    };

    template<int N>
    class ITensor4
    {
    public:
        virtual double  Component(int i, int j, int k, int l) const = 0;
        virtual double& Component(int i, int j, int k, int l) = 0;
    };

    template<int N>
    class ITensorField4
    {
    public:
        virtual double Component(int i, int j, int k, int l, const VectorN<Real, N> &pos) const = 0;
        virtual void   ValueAtPoint(const VectorN<Real, N> &pos, ITensor4<N> &outVal) const = 0;
    };        

}


///////////////////////////   ./include/core/Tensor.h   ///////////////////////////



namespace MML
{
    template <int N>
    class Tensor2 : public ITensor2<N>
    {
        double _coeff[N][N];
    public:
        int _numContravar;
        int _numCovar;

        Tensor2(int nContra, int nCo) : _numContravar(nContra), _numCovar(nCo) {}

        double  Component(int i, int j) const { return _coeff[i][j]; }
        double& Component(int i, int j)       { return _coeff[i][j];}    
    };

    template <int N>
    class Tensor3 : public ITensor3<N>
    {
        double _coeff[N][N][N];
    public:
        int _numContravar;
        int _numCovar;

        Tensor3(int nContra, int nCo) : _numContravar(nContra), _numCovar(nCo) {}
        
        double  Component(int i, int j, int k) const { return _coeff[i][j][k]; }
        double& Component(int i, int j, int k)       { return _coeff[i][j][k];}   
    };

    template <int N>
    class Tensor4 : public ITensor4<N>
    {
        double _coeff[N][N][N][N];
    public:
        int _numContravar;
        int _numCovar;

        Tensor4(int nContra, int nCo) : _numContravar(nContra), _numCovar(nCo) {}
        
        double  Component(int i, int j, int k, int l) const { return _coeff[i][j][k][l]; }
        double& Component(int i, int j, int k, int l)       { return _coeff[i][j][k][l];}   
    };
}
///////////////////////////   ./include/core/Constants.h   ///////////////////////////

namespace MML
{
    class Constants
    {
    public:
        static inline const Real PI = 3.141593;
        static inline const Real Epsilon = std::numeric_limits<Real>::epsilon();
        // static constexpr Real EpsilonSqrt = std::sqrt(Epsilon);
    };
}
///////////////////////////   ./include/core/CoreUtils.h   ///////////////////////////


namespace MML
{
    ///////////////////////////////////////////////   Interfaces    ///////////////////////////////////////////
    class Utils
    {
    public:
        static Matrix<Complex> CmplxMatFromRealMat(const Matrix<Real> &mat)
        {
            Matrix<Complex> mat_cmplx(mat.RowNum(), mat.ColNum());

            for(int i=0; i<mat.RowNum(); i++ )
                for(int j=0; j<mat.ColNum(); j++)
                    mat_cmplx[i][j] = Complex(mat(i,j), 0.0);
            
            return mat_cmplx;         
        }

        // IsComplexMatReal?

        // enabling Complex - Real vector operations
        static Vector<Complex> AddVec(const Vector<Complex> &a, const Vector<Real> &b)
        {
            Vector<Complex> ret(b.size());;
            for (int i = 0; i < b.size(); i++)
                ret[i] = a[i] + b[i];
            return ret;
        }
        static Vector<Complex> AddVec(const Vector<Real> &a, const Vector<Complex> &b)
        {
            Vector<Complex> ret(b.size());;
            for (int i = 0; i < b.size(); i++)
                ret[i] = a[i] + b[i];
            return ret;
        }

        static  Vector<Complex> SubVec(const Vector<Complex> &a, const Vector<Real> &b)
        {
            Vector<Complex> ret(b.size());;
            for (int i = 0; i < b.size(); i++)
                ret[i] = a[i] - b[i];
            return ret;
        }
        static Vector<Complex> SubVec(const Vector<Real> &a, const Vector<Complex> &b)
        {
            Vector<Complex> ret(b.size());;
            for (int i = 0; i < b.size(); i++)
                ret[i] = a[i] - b[i];
            return ret;
        }        

        // enabling Complex - Real Matrix operations
        static Matrix<Complex> AddMat(const Matrix<Complex> &a, const Matrix<Real> &b)
        {
            if (a.RowNum() != b.RowNum() || a.ColNum() != b.ColNum())
                throw MatrixDimensionError("AddMat(Complex, Real) - must be same dim", a.RowNum(), a.ColNum(), b.RowNum(), b.ColNum());

            Matrix<Complex> ret(a);
            for (int i = 0; i < b.RowNum(); i++)
                for (int j = 0; j < b.ColNum(); j++)
                    ret[i][j] += b[i][j];
            return ret;
        }
        static Matrix<Complex> AddMat(const Matrix<Real> &a, const Matrix<Complex> &b)
        {
            if (a.RowNum() != b.RowNum() || a.ColNum() != b.ColNum())
                throw MatrixDimensionError("AddMat(Real, Complex) - must be same dim", a.RowNum(), a.ColNum(), b.RowNum(), b.ColNum());

            Matrix<Complex> ret(b);
            for (int i = 0; i < a.RowNum(); i++)
                for (int j = 0; j < a.ColNum(); j++)
                    ret[i][j] += a[i][j];
            return ret;
        }

        static Matrix<Complex> SubMat(const Matrix<Complex> &a, const Matrix<Real> &b)
        {
            if (a.RowNum() != b.RowNum() || a.ColNum() != b.ColNum())
                throw MatrixDimensionError("AddMat(Complex, Real) - must be same dim", a.RowNum(), a.ColNum(), b.RowNum(), b.ColNum());

            Matrix<Complex> ret(a);
            for (int i = 0; i < b.RowNum(); i++)
                for (int j = 0; j < b.ColNum(); j++)
                    ret[i][j] -= b[i][j];
            return ret;
        }
        static Matrix<Complex> SubMat(const Matrix<Real> &a, const Matrix<Complex> &b)
        {
            if (a.RowNum() != b.RowNum() || a.ColNum() != b.ColNum())
                throw MatrixDimensionError("AddMat(Real, Complex) - must be same dim", a.RowNum(), a.ColNum(), b.RowNum(), b.ColNum());

            Matrix<Complex> ret(b);
            for (int i = 0; i < a.RowNum(); i++)
                for (int j = 0; j < a.ColNum(); j++)
                    ret[i][j] = a[i][j] - b[i][j];
            return ret;
        }

        static Matrix<Complex> MulMat(const Matrix<Complex> &a, const Matrix<Real> &b)
        {
            if( a.ColNum() != b.RowNum() )
                throw MatrixDimensionError("Matrix::operator*() - a.colNum must be equal to b.rowNum", a.RowNum(), a.ColNum(), b.RowNum(), b.ColNum());

            Matrix<Complex>	ret(a.RowNum(), b.ColNum());
            for( int i=0; i<ret.RowNum(); i++ )
                for( int j=0; j<ret.ColNum(); j++ ) {
                    ret[i][j] = 0.0;
                    for(int k=0; k<a.ColNum(); k++ )
                        ret[i][j] += a[i][k] * b[k][j];
                }

            return	ret;
        }
        static Matrix<Complex> MulMat(const Matrix<Real> &a, const Matrix<Complex> &b)
        {
            if( a.ColNum() != b.RowNum() )
                throw MatrixDimensionError("Matrix::operator*() - a.colNum must be equal to b.rowNum", a.RowNum(), a.ColNum(), b.RowNum(), b.ColNum());

            Matrix<Complex>	ret(a.RowNum(), b.ColNum());
            for( int i=0; i<ret.RowNum(); i++ )
                for( int j=0; j<ret.ColNum(); j++ ) {
                    ret[i][j] = 0.0;
                    for(int k=0; k<a.ColNum(); k++ )
                        ret[i][j] += a[i][k] * b[k][j];
                }

            return	ret;
        }  
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
    };

    //////////////////////////////////////////////////////////////////////
    class IRealFunction : public IFunction<double, double>
    {
        public:
        virtual double operator()(double) const = 0;
    };

    //////////////////////////////////////////////////////////////////////
    template<int N>
    class IScalarFunction : public IFunction<double, const VectorN<Real, N> &>
    {
        public:
        virtual double operator()(const VectorN<Real, N> &x) const = 0;
    };
    
    //////////////////////////////////////////////////////////////////////
    template<int N>
    class IRealToVectorFunction : public IFunction<VectorN<Real, N>, double>
    {
        public:
        virtual VectorN<Real, N> operator()(double x) const = 0;
    };

    //////////////////////////////////////////////////////////////////////
    template<int N>
    class IVectorFunction : public IFunction<VectorN<Real, N>, const VectorN<Real, N> &>
    {
        public:
        virtual VectorN<Real, N> operator()(const VectorN<Real, N> &x) const = 0;
        virtual double operator()(const VectorN<Real, N> &x, int component) const
        {
            VectorN<Real, N> val = (*this)(x);
            return val[component];
        }
    };

    //////////////////////////////////////////////////////////////////////
    template<int N, int M>
    class IVectorFunctionNM : public IFunction<VectorN<Real, M>, const VectorN<Real, N> &>
    {
        public:
        virtual VectorN<Real, M> operator()(const VectorN<Real, N> &x) const = 0;
        virtual double operator()(const VectorN<Real, N> &x, int component) const
        {
            VectorN<Real, M> val = (*this)(x);
            return val[component];
        }
    };

    //////////////////////////////////////////////////////////////////////
    template<int N>
    class IParametricCurve : public IRealToVectorFunction<N> // IFunction<VectorN<Real, N>, double>
    {
        public:
        virtual VectorN<Real, N> operator()(double x) const = 0;

        // GetMixX(), GetMaxX(), može vracati i infinity
    };

    //////////////////////////////////////////////////////////////////////
    template<int N>
    class IParametricSurface : public IFunction<VectorN<Real, N>, const VectorN<Real, 2> &>
    {
        public:
        virtual VectorN<Real, N> operator()(double u, double w) const = 0;
        virtual VectorN<Real, N> operator()(const VectorN<Real, 2> &coord) const 
        {
            return operator()(coord[0], coord[1]);
        }

        // GetMixX(), GetMaxX(), može vracati i infinity
        // GetMixY(), GetMaxY(), može vracati i infinity
        // da je povrsina omedjena

    };
}

///////////////////////////   ./include/interfaces/ICoordTransf.h   ///////////////////////////



namespace MML
{
    template<typename VectorFrom, typename VectorTo, int N>
    class ICoordTransf
    {
        public:
        virtual VectorTo    transf(const VectorFrom &in) = 0;       // transfCoord
        virtual VectorFrom  transfInverse(const VectorTo &in) = 0;

        virtual IScalarFunction<N>& coordTransfFunc(int i) = 0;
        virtual IScalarFunction<N>& inverseCoordTransfFunc(int i) = 0;
    };    
  
}

///////////////////////////   ./include/interfaces/IODESystem.h   ///////////////////////////


namespace MML
{
	class IODESystem
	{
	public:
        virtual int     getDim() = 0;
		virtual void    derivs(const double, const Vector<Real>&, Vector<Real> &) = 0;
		void operator()(const double t, const Vector<Real> &x, Vector<Real> &dxdt) { return derivs(t, x, dxdt); }
	};  
}

///////////////////////////   ./include/basic_types/Polynom.h   ///////////////////////////


namespace MML
{
    template <typename _Field, typename _CoefType = double>
    class Polynom
    {
    private:
        std::vector<_CoefType> _vecCoef;
    public:
        Polynom() {}
        Polynom(int n) { _vecCoef.resize(n+1); }
        Polynom(const std::vector<_CoefType> &vecCoef) : _vecCoef(vecCoef) {}
        Polynom(std::initializer_list<_CoefType> list) : _vecCoef(list) {}
        Polynom(const Polynom &Copy) : _vecCoef(Copy._vecCoef) {}
        ~Polynom() {}

        Polynom& operator=(const Polynom &Copy) { _vecCoef = Copy._vecCoef; return *this; }

        int  GetDegree() const     { return (int) _vecCoef.size() - 1; }
        void SetDegree(int newDeg) {  _vecCoef.resize(newDeg+1); }

        _Field  operator[] (int i) const { return _vecCoef[i]; }
        _Field& operator[] (int i)       { return _vecCoef[i]; }

        _Field operator() (const _Field &x) {
            int j = GetDegree();
            _Field p = _vecCoef[j];
            
            while (j>0) 
                p = p*x + _vecCoef[--j];
            return p;
        }        
        // Given the coefficients of a polynomial of degree nc as an array c[0..nc] of size nc+1 (with
        // c[0] being the constant term), and given a value x, this routine fills an output array pd of size
        // nd+1 with the value of the polynomial evaluated at x in pd[0], and the first nd derivatives at
        // x in pd[1..nd].
        void Derive(const Real x, Vector<Real> &pd)
        {
            int  nnd,j,i;
            int  nc=GetDegree();
            int  nd=pd.size()-1;
            Real cnst=1.0;

            pd[0]=(*this)[nc];
            for (j=1;j<nd+1;j++) 
                pd[j]=0.0;
            
            for (i=nc-1;i>=0;i--) 
            {
                nnd=(nd < (nc-i) ? nd : nc-i);
                for (j=nnd;j>0;j--) 
                    pd[j]=pd[j]*x+pd[j-1];

                pd[0]=pd[0]*x + (*this)[i];
            }
            for (i=2;i<nd+1;i++) {
                cnst *= i;
                pd[i] *= cnst;
            }
        }        

        bool operator==(const Polynom &b) const
        {
            if( _vecCoef.size() != b._vecCoef.size() )
                return false;
            for( int i=0; i<_vecCoef.size(); i++ )
                if( _vecCoef[i] != b._vecCoef[i] )
                    return false;
            return true;
        }

        Polynom operator+(const Polynom &b) const
        {
            Polynom result;
            int n = (int) std::max(_vecCoef.size(), b._vecCoef.size());
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

        Polynom operator-(const Polynom &b) const
        {
            Polynom result;
            int n = (int) std::max(_vecCoef.size(), b._vecCoef.size());
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

        Polynom operator*(const Polynom &b) const
        {
            Polynom result;

            int n = (int) (_vecCoef.size() + b._vecCoef.size() - 1);
            result._vecCoef.resize(n);
            for (int i = 0; i < _vecCoef.size(); i++)
                for (int j = 0; j < b._vecCoef.size(); j++)
                    result._vecCoef[i + j] += _vecCoef[i] * b._vecCoef[j];
            return result;
        }

        static void poldiv(const Polynom &u, const Polynom &v, Polynom &qout, Polynom &rout)
        {
            int k,j,n=u.GetDegree(),nv=v.GetDegree();

            while (nv >= 0 && v._vecCoef[nv] == 0.) nv--;

            if (nv < 0) 
                throw("poldiv divide by zero polynomial");
            
            Polynom r = u;
            Polynom q(u.GetDegree());
            for (k=n-nv;k>=0;k--) 
            {
                q[k]=r[nv+k]/v[nv];
                for (j=nv+k-1;j>=k;j--) 
                    r[j] -= q[k]*v[j-k];
            }
            for (j=nv;j<=n;j++) r[j]=0.0;

            int nq=q.GetDegree();            
            while (nq >= 0 && q[nv] == 0.) 
                nq--;

            qout.SetDegree(nq-1);
            for(j=0; j<nq; j++ )
                qout[j] = q[j];

            rout.SetDegree(nv-1);
            for(j=0; j<nv; j++ )
                rout[j] = r[j];

        }
        
        friend Polynom operator*(const Polynom &a, _CoefType b )
        {
            Polynom ret;
            ret._vecCoef.resize(a.GetDegree()+1);
            for(int i=0; i<=a.GetDegree(); i++)
                ret._vecCoef[i] = a._vecCoef[i] * b;
            return ret;
        }

        friend Polynom operator*(_CoefType a, const Polynom &b )
        {
            Polynom ret;
            ret._vecCoef.resize(b.GetDegree()+1);
            for(int i=0; i<=b.GetDegree(); i++)
                ret._vecCoef[i] = a * b._vecCoef[i];
            return ret;
        }

        friend Polynom operator/(const Polynom &a, _CoefType b)
        {
            Polynom ret;
            ret._vecCoef.resize(a.GetDegree()+1);
            for(int i=0; i<=a.GetDegree(); i++)
                ret._vecCoef[i] = a._vecCoef[i] / b;
            return ret;
        }        

        std::ostream& Print(std::ostream& stream, int width, int precision) const
        {
            for(int i=(int) _vecCoef.size()-1; i>=0; i--)
            {
                if( std::abs(_vecCoef[i]) == 0.0 )
                    continue;

                if( i != _vecCoef.size()-1 )
                    stream << " + ";

                if( i == 0 )
                    stream << std::setw(width) << std::setprecision(precision) << _vecCoef[i];
                else if ( i == 1 )
                    stream << std::setw(width) << std::setprecision(precision) << _vecCoef[i] << " * x" << i;
                else
                    stream << std::setw(width) << std::setprecision(precision) << _vecCoef[i] << " * x^" << i;
            }

            return stream;
        }

        std::ostream& Print(std::ostream& stream) const
        {
            for(int i=(int)_vecCoef.size()-1; i>=0; i--)
            {
                if( std::abs(_vecCoef[i]) == 0.0 )
                    continue;

                if( i != _vecCoef.size()-1 )
                    stream << " + ";
                
                if( i == 0 )
                    stream <<  _vecCoef[i];
                else if ( i == 1 )
                    stream <<  _vecCoef[i] << " * x" << i;
                else
                    stream <<  _vecCoef[i] << " * x^" << i;
            }

            return stream;
        }

        std::string to_string(int width, int precision) const
        {
            std::stringstream str;

            Print(str, width, precision);

            return str.str();
        }

        friend std::ostream& operator<<(std::ostream& stream, Polynom &a)
        {
            a.Print(stream);

            return stream;
        }        
    };
        
    typedef Polynom<Real, Real>         RealPolynom;
    typedef Polynom<Complex, Complex>   ComplexPolynom;

    typedef Polynom<MatrixNM<Real,2,2>, Real>       MatrixPolynomDim2;
    typedef Polynom<MatrixNM<Real,3,3>, Real>       MatrixPolynomDim3;
    typedef Polynom<MatrixNM<Real,4,4>, Real>       MatrixPolynomDim4;
}

///////////////////////////   ./include/basic_types/Geometry.h   ///////////////////////////


namespace MML
{
    class Triangle
    {
    private:
        double _a, _b, _c;

    public:
        double  A() const   { return _a; }
        double& A()         { return _a; }
        double  B() const   { return _b; }
        double& B()         { return _b; }
        double  C() const   { return _c; }
        double& C()         { return _c; }

        Triangle() {}
        Triangle(double a, double b, double c) : _a(a), _b(b), _c(c) {}

        double Area() const
        {
            double s = (_a + _b + _c) / 2.0;
            return sqrt(s * (s - _a) * (s - _b) * (s - _c));
        }
        bool IsRight() const
        {
            return (SQR(_a) + SQR(_b) == SQR(_c)) || (SQR(_a) + SQR(_c) == SQR(_b)) || (SQR(_b) + SQR(_c) == SQR(_a));
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

///////////////////////////   ./include/basic_types/Geometry2D.h   ///////////////////////////


namespace MML
{
    class Point2Cartesian
    {
    private:
        double _x, _y;

    public:
        double  X() const   { return _x; }
        double& X()         { return _x; }
        double  Y() const   { return _y; }
        double& Y()         { return _y; }

        Point2Cartesian() {}
        Point2Cartesian(double x, double y) : _x(x), _y(y) {}

        double Dist(const Point2Cartesian &b) const { return sqrt(SQR(b._x - _x) + SQR(b._y - _y) ); }

        Point2Cartesian operator+(const VectorN<Real, 2> &b) const { return Point2Cartesian(_x + b[0], _y + b[1]); }
        Point2Cartesian operator-(const VectorN<Real, 2> &b) const { return Point2Cartesian(_x - b[0], _y - b[1]); }
    };

    class Vector2Cartesian : public VectorN<Real, 2>
    {
    public:
        Vector2Cartesian() {}
        Vector2Cartesian(double x, double y) 
        {
            _val[0] = x;
            _val[1] = y;
        }
        Vector2Cartesian(const Point2Cartesian &a, const Point2Cartesian &b) 
        {
            _val[0] = b.X() - a.X();
            _val[1] = b.Y() - a.Y();
        }

        Real    X() const   { return _val[0]; }
        Real&   X()         { return _val[0]; }
        Real    Y() const   { return _val[1]; }
        Real&   Y()         { return _val[1]; }

        Vector2Cartesian GetUnitVector() const
        {
            VectorN<Real, 2> res = (*this) / NormL2();

            return Vector2Cartesian(res[0], res[1]);
        }
    };

    class Point2Polar
    {
    private:
        double _r, _phi;

    public:
        double  R() const   { return _r; }
        double& R()         { return _r; }
        double  Phi() const { return _phi; }
        double& Phi()       { return _phi; }
        
        Point2Polar() {}
        Point2Polar(double r, double phi) : _r(r), _phi(phi) {}
    };

    class Vector2Polar : public VectorN<Real, 2>
    {
    public:
        Real    R() const   { return _val[0]; }
        Real&   R()         { return _val[0]; }
        Real    Phi() const { return _val[1]; }
        Real&   Phi()       { return _val[1]; }

        Vector2Polar() {}
        Vector2Polar(double r, double phi) 
        {
            _val[0] = r;
            _val[1] = phi;
        }   
    };

    class Line2D
    {
    private:
        Point2Cartesian _point;
        Vector2Cartesian _direction; // unit vector in line direction

    public:
        Line2D(const Point2Cartesian &pnt, const Vector2Cartesian dir)
        {
            _point = pnt;
            _direction = dir.GetUnitVector();
        }

        Line2D(const Point2Cartesian &a, const Point2Cartesian &b)
        {
            Vector2Cartesian dir(a , b);
            _point = a;
            _direction = dir.GetUnitVector();
        }

        Point2Cartesian     StartPoint() const  { return _point; }
        Point2Cartesian&    StartPoint()        { return _point; }

        Vector2Cartesian    Direction() const   { return _direction; }
        Vector2Cartesian&   Direction()         { return _direction; }

        Point2Cartesian PointOnLine(double t)
        {
            VectorN<Real, 2> dist = t * _direction;
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

        SegmentLine2D(Point2Cartesian pnt1, Vector2Cartesian direction, double t) : _point1(pnt1)
        {
            _point2 = pnt1 + t * direction;
        }

        Point2Cartesian     StartPoint() const  { return _point1; }
        Point2Cartesian&    StartPoint()        { return _point1; }

        Point2Cartesian     EndPoint()  const  { return _point2; }
        Point2Cartesian&    EndPoint()         { return _point2; }

        Point2Cartesian PointOnSegment(double t)
        {
            // check t  u [0,1]
            VectorN<Real, 2> dist = t * Direction();
            Point2Cartesian ret = _point1 + dist;
            return ret;
        }

        double              Length()    const { return _point1.Dist(_point2); }
        Vector2Cartesian    Direction() const { return Vector2Cartesian(_point1, _point2); }
    };

    class Triangle2D
    {
    private:
        Point2Cartesian _pnt1, _pnt2, _pnt3;

    public:
        Triangle2D(Point2Cartesian pnt1, Point2Cartesian pnt2, Point2Cartesian pnt3) : _pnt1(pnt1), _pnt2(pnt2), _pnt3(pnt3)
        { }

        Point2Cartesian     Pnt1() const { return _pnt1; }
        Point2Cartesian&    Pnt1()       { return _pnt1; }
        Point2Cartesian     Pnt2() const { return _pnt2; }
        Point2Cartesian&    Pnt2()       { return _pnt2; }
        Point2Cartesian     Pnt3() const { return _pnt3; }
        Point2Cartesian&    Pnt3()       { return _pnt3; }

        double Area() const
        {
            double a = _pnt1.Dist(_pnt2);
            double b = _pnt2.Dist(_pnt3);
            double c = _pnt3.Dist(_pnt1);

            double s = (a + b + c) / 2.0;

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
        std::vector<Point2Cartesian>& Points()       { return _points; }

        double Area() const
        {
            double area = 0.0;
            int n = (int) _points.size();
            for (int i = 0; i < n; i++)
            {
                area += _points[i].X() * _points[(i + 1) % n].Y();
                area -= _points[i].Y() * _points[(i + 1) % n].X();
            }
            area /= 2.0;
            return area;
        }

        // TODO - IsConvex()
        

        // bool IsConvex() const
        // {
        //     int n = _points.size();
        //     if (n < 3) return false;

        //     bool got_negative = false;
        //     bool got_positive = false;
        //     int B, C;
        //     for (int A = 0; A < n; A++)
        //     {
        //         B = (A + 1) % n;
        //         C = (B + 1) % n;

        //         double cross_product = CrossProductLength(_points[A], _points[B], _points[C]);
        //         if (cross_product < 0)
        //         {
        //             got_negative = true;
        //         }
        //         else if (cross_product > 0)
        //         {
        //             got_positive = true;
        //         }
        //         if (got_negative && got_positive) return false;
        //     }

        //     return true;
        // }

        bool IsInside(Point2Cartesian pnt) const
        {
            return false;
        }                 
    };            
}
///////////////////////////   ./include/basic_types/Geometry3D.h   ///////////////////////////



namespace MML
{
    class Point3Cartesian
    {
    private:
        double _x, _y, _z;

    public:
        double  X() const   { return _x; }
        double& X()         { return _x; }
        double  Y() const   { return _y; }
        double& Y()         { return _y; }
        double  Z() const   { return _z; }
        double& Z()         { return _z; }

        Point3Cartesian() {}
        Point3Cartesian(double x, double y, double z) : _x(x), _y(y), _z(z) {}

        double Dist(const Point3Cartesian &b) const { return sqrt(SQR(b._x - _x) + SQR(b._y - _y) + SQR(b._z - _z)); }

        Point3Cartesian operator+(const VectorN<Real, 3> &b) const { return Point3Cartesian(_x + b[0], _y + b[1], _z + b[2]); }
        Point3Cartesian operator-(const VectorN<Real, 3> &b) const { return Point3Cartesian(_x - b[0], _y - b[1], _z - b[2]); }
    };

    class Vector3Cartesian : public VectorN<Real, 3>
    {
    public:
        Real    X() const   { return _val[0]; }
        Real&   X()         { return _val[0]; }
        Real    Y() const   { return _val[1]; }
        Real&   Y()         { return _val[1]; }
        Real    Z() const   { return _val[2]; }
        Real&   Z()         { return _val[2]; }

        Vector3Cartesian()                                  : VectorN<Real,3>{0.0, 0.0, 0.0} { }
        Vector3Cartesian(const VectorN<Real, 3> &b)         : VectorN<Real,3>{b[0], b[1], b[2]} { }
        Vector3Cartesian(double x, double y, double z)      : VectorN<Real,3>{x, y, z} { }
        Vector3Cartesian(std::initializer_list<Real> list)  : VectorN<Real,3>(list) { }
        Vector3Cartesian(const Point3Cartesian &a, const Point3Cartesian &b) 
        {
            _val[0] = b.X() - a.X();
            _val[1] = b.Y() - a.Y();
            _val[2] = b.Z() - a.Z();
        }

        Point3Cartesian getAsPoint()
        {
            return Point3Cartesian(_val[0], _val[1], _val[2]);
        }
        bool IsParallelTo(Vector3Cartesian &b, Real eps = 1e-15) const
        {
            Real norm1 = NormL2();
            Real norm2 = b.NormL2();

            return std::abs(X() / norm1 - b.X() / norm2) < eps &&
                   std::abs(Y() / norm1 - b.Y() / norm2) < eps &&
                   std::abs(Z() / norm1 - b.Z() / norm2) < eps;
        }

        bool IsPerpendicularTo(Vector3Cartesian &b, Real eps = 1e-15) const
        {
            if (std::abs(ScalarProd(*this, b)) < eps)
                return true;
            else
                return false;
        }

        Vector3Cartesian GetUnitVector() const
        {
            VectorN<Real, 3> res = (*this) / NormL2();

            return Vector3Cartesian{res};
        }

        friend Real ScalarProd(const Vector3Cartesian &a, const Vector3Cartesian &b)
        {
            return a.ScalarProductCartesian(b);
        }

        friend Vector3Cartesian VectorProd(const Vector3Cartesian &a, const Vector3Cartesian &b)
        {
            Vector3Cartesian ret;

            ret.X() = a.Y() * b.Z() - a.Z() * b.Y();
            ret.Y() = a.Z() * b.X() - a.X() * b.Z();
            ret.Z() = a.X() * b.Y() - a.Y() * b.X();

            return ret;
        }

        friend Real VectorsAngle(const Vector3Cartesian &a, const Vector3Cartesian &b)
        {
            Real cos_phi = ScalarProd(a, b) / (a.NormL2() * b.NormL2());

            return acos(cos_phi);
        }
    };

    class Vector3Spherical : public VectorN<Real, 3>
    {
    public:
        Real    R()     const   { return _val[0]; }
        Real&   R()             { return _val[0]; }
        Real    Theta() const   { return _val[1]; }
        Real&   Theta()         { return _val[1]; }
        Real    Phi()   const   { return _val[2]; }
        Real&   Phi()           { return _val[2]; }

        Vector3Spherical()                                   : VectorN<Real,3>{0.0, 0.0, 0.0} { }
        Vector3Spherical(const VectorN<Real, 3> &b)          : VectorN<Real,3>{b[0], b[1], b[2]} { }
        Vector3Spherical(double r, double theta, double phi) : VectorN<Real,3>{r, theta, phi} { }
        Vector3Spherical(std::initializer_list<Real> list)   : VectorN<Real,3>(list) { }
    };

    class Vector3Cylindrical : public VectorN<Real, 3>
    {
    public:
        Real    R()   const { return _val[0]; }
        Real&   R()         { return _val[0]; }
        Real    Phi() const { return _val[1]; }
        Real&   Phi()       { return _val[1]; }
        Real    Z()   const { return _val[2]; }
        Real&   Z()         { return _val[2]; }

        Vector3Cylindrical()                                 : VectorN<Real,3>{0.0, 0.0, 0.0} { }
        Vector3Cylindrical(const VectorN<Real, 3> &b)        : VectorN<Real,3>{b[0], b[1], b[2]} { }
        Vector3Cylindrical(double r, double phi, double z)   : VectorN<Real,3>{r, phi, z} { }
        Vector3Cylindrical(std::initializer_list<Real> list) : VectorN<Real,3>(list) { }
    };

    class Line3D
    {
    private:
        Point3Cartesian _point;
        Vector3Cartesian _direction; 

    public:
        Line3D() {}
        Line3D(const Point3Cartesian &pnt, const Vector3Cartesian dir)
        {
            _point = pnt;
            _direction = dir.GetUnitVector();
        }

        Line3D(const Point3Cartesian &a, const Point3Cartesian &b)
        {
            Vector3Cartesian dir(a, b);
            _point = a;
            _direction = dir.GetUnitVector();
        }

        Point3Cartesian StartPoint() const  { return _point; }
        Point3Cartesian &StartPoint()       { return _point; }

        Vector3Cartesian Direction() const  { return _direction; }
        Vector3Cartesian &Direction()       { return _direction; }

        Point3Cartesian PointOnLine(double t)
        {
            return _point + t * _direction;
        }

        bool IsPerpendicular(const Line3D &b ) const
        {
            return ScalarProd(Direction(), b.Direction()) == 0.0f;
        }
        // TODO distance Line - Point3
        // TODO nearest point on line
        // TODO pravac koji prolazi kroz tocku i sijece zadani pravac okomito
    };

    class SegmentLine3D
    {
    private:
        Point3Cartesian _point1;
        Point3Cartesian _point2;

    public:
        SegmentLine3D(Point3Cartesian pnt1, Point3Cartesian pnt2) : _point1(pnt1), _point2(pnt2)
        { }

        SegmentLine3D(Point3Cartesian pnt1, Vector3Cartesian direction, double t)
        {
            _point1 = pnt1;
            _point2 = pnt1 + t * direction;
        }

        Point3Cartesian StartPoint() const  { return _point1; }
        Point3Cartesian &StartPoint()       { return _point1; }

        Point3Cartesian EndPoint() const  { return _point2; }
        Point3Cartesian &EndPoint()       { return _point2; }

        Point3Cartesian PointOnSegment(double t)
        {
            // check t  u [0,1]
            VectorN<Real, 3> dist = t * Direction();
            Point3Cartesian ret = _point1 + dist;
            return ret;
        }

        double              Length()    const { return _point1.Dist(_point2); }
        Vector3Cartesian    Direction() const { return Vector3Cartesian(_point1, _point2); }
    };

    class Plane3D
    {
    private:
        double _A, _B, _C, _D;
    
    public:
        Plane3D(const Point3Cartesian &a, const Vector3Cartesian &normal)
        {
            // check for normal null vector
            Vector3Cartesian unitNormal = normal.GetUnitVector();

            _A = unitNormal.X();
            _B = unitNormal.Y();
            _C = unitNormal.Z();
            _D = -(a.X() * unitNormal.X() + a.Y() * unitNormal.Y() + a.Z() * unitNormal.Z());
        }

        Plane3D(const Point3Cartesian &a, const Point3Cartesian &b, const Point3Cartesian &c) : Plane3D(a, VectorProd(Vector3Cartesian(a, b), Vector3Cartesian(a, c)))
        { }

        Plane3D(double alpha, double beta, double gamma, double d)      // Hesseov (normalni) oblik
        {
            _A = cos(alpha);
            _B = cos(beta);
            _C = cos(gamma);
            _D = -d;
        }

        // tri segmenta na koord osima ctor
        Plane3D(double seg_x, double seg_y, double seg_z) 
        {
            Point3Cartesian x(seg_x, 0, 0);
            Point3Cartesian y(0, seg_y, 0);
            Point3Cartesian z(0, 0, seg_z);

            // TODO - if seg == 0
            _A = 1 / seg_x;
            _B = 1 / seg_y;
            _C = 1 / seg_z;
            _D = -1;
        }

        static Plane3D GetXYPlane() { return Plane3D(Point3Cartesian(0,0,0), Vector3Cartesian(0,0,1)); }
        static Plane3D GetXZPlane() { return Plane3D(Point3Cartesian(0,0,0), Vector3Cartesian(0,1,0)); }
        static Plane3D GetYZPlane() { return Plane3D(Point3Cartesian(0,0,0), Vector3Cartesian(1,0,0)); }

        double  A() const   { return _A; }
        double& A()         { return _A; }
        double  B() const   { return _B; }
        double& B()         { return _B; }
        double  C() const   { return _C; }
        double& C()         { return _C; }
        double  D() const   { return _D; }
        double& D()         { return _D; }

        Vector3Cartesian Normal() const { return Vector3Cartesian(_A, _B, _C); }

        void GetCoordAxisSegments(double &outseg_x, double& outseg_y, double& outseg_z)
        {
            outseg_x = - _D /  _A;
            outseg_y = - _D /  _B;
            outseg_z = - _D /  _C;
        }
        // TODO - GetHesseNormalFormParams


        bool IsPointOnPlane(const Point3Cartesian &pnt) const
        {
            return _A * pnt.X() + _B * pnt.Y() + _C * pnt.Z() + _D == 0.0;
        }
        double DistToPoint(const Point3Cartesian &pnt) const
        {
            double a = _A * pnt.X() + _B * pnt.Y() + _C * pnt.Z() + _D;
            double b = sqrt(_A * _A + _B * _B + _C * _C);

            return std::abs(a / b);
        }
        Point3Cartesian ProjectionToPlane(const Point3Cartesian &pnt) const
        {
            Point3Cartesian ret;
            return ret;
        }

        bool IsLineOnPlane(const Line3D &line) const
        {
            return false;
        }
        double AngleToLine(const Line3D &line) const
        {
            return 0.0;
        }        
        bool IntersectionWithLine(const Line3D &line, Point3Cartesian &out_inter_pnt) const
        {
            return false;
        }

        bool IsParallelToPlane(const Plane3D &plane) const
        {
            Vector3Cartesian norm1(_A, _B, _C);
            Vector3Cartesian norm2(plane._A, plane._B, plane._C);

            return norm1.IsParallelTo(norm2);
        }        
        bool IsPerpendicularToPlane(const Plane3D &plane) const
        {
            Vector3Cartesian norm1(_A, _B, _C);
            Vector3Cartesian norm2(plane._A, plane._B, plane._C);

            return norm1.IsPerpendicularTo(norm2);
        }
        double AngleToPlane(const Plane3D &plane) const
        {
            return 0.0;
        }     

        bool IntersectionWithPlane(const Plane3D &plane, Line3D &out_inter_line) const
        {
            return false;
        }
    };

    class Triangle3D
    {
    private:
        Point3Cartesian _pnt1, _pnt2, _pnt3;;
    };

    typedef Vector3Cartesian    Vec3Cart;
    typedef Vector3Spherical    Vec3Sph;
    typedef Vector3Cylindrical  Vec3Cyl;
}

///////////////////////////   ./include/basic_types/CoordSystem.h   ///////////////////////////



namespace MML
{
/*
SUSTINA
- imam koordinate u jednom, i (lokalne) koordinate u drugom koordinatnom sistemu
- kakva je veza?
- i kako izracunati jedne na osnovu drugih?

Kakvi koordinatni sistemi su interesantni
1. Ortogonal cartesian
2. Oblique cartesian
3. Cylindrical
4. Spherical
5. Generalized
6. Curvilinear
7. Generalized curvilinear
8. Orthogonal curvilinear
9. Oblique curvilinear
10. Rectilinear
11. Rectilinear orthogonal
12. Rectilinear oblique

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
            double a01 = ScalarProd(_base[0], _base[1]);
            double a02 = ScalarProd(_base[0], _base[2]);
            double a12 = ScalarProd(_base[1], _base[2]);

            return sqrt(a01*a01 + a02*a02 + a12*a12) < 1e-6;
        }
    };

    class CoordSystemObliqueCartesian
    {
        public:
        Vector3Cartesian _base[3];
        VectorN<double, 3> _dual[3];

        MatrixNM<Real,3,3> _alpha;

        MatrixNM<Real,3,3> _transf;
        MatrixNM<Real,3,3> _inv;

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

            for( int i=0; i<3; i++ )
            {
                for(int j=0; j<3; j++ )
                {
                    _alpha(i,j) =  _base[i][j];
                    _transf(i,j) = _base[j][i];     // transponirano
                }
            }

            _inv = _transf.GetInverse();

        }

        VectorN<double, 3> Base(int i) { return _base[i]; }
        VectorN<double, 3> Dual(int i) { return _dual[i]; }        
    };

    // RECTILINEAR - konstantni bazni vektori
    // CURVILINEAR - U SVAKOJ TOCKI IMA DRUGACIJE BAZNE VEKTORE
    template<int N>
    class CoordSystem
    {
        VectorN<double, N> _base[N];
        VectorN<double, N> _originCoord;

        // init s tri vektora baze
        // init s matricom transformacije?

        // isOrthogonal

        virtual void transf(VectorN<double, N> &x) = 0;
    };

    class CoordSystemTranslated
    {
        VectorN<double, 3> _originCoord;
    };

    class CoordSystemRotated
    {
        Vector3Cartesian _originCoord;
        // matrica transf

        // ima staticke membere za kreiranje bazicnih matrixa rotacija
    };

    template<int N>
    class MovingCoordSystem
    {
        VectorN<double, 3> _origin;
        virtual VectorN<Real, 3> transf(const VectorN<Real, 3> &x, double t) = 0;

    };    

    // rotacije - staticki
    // rotacije + translacije - staticki

    class InertialMovingFrame : public MovingCoordSystem<3>
    {
        public:
        // origin position at t = 0
        Vector3Cartesian _origin;

        Vector3Cartesian _speed;

        InertialMovingFrame() {}
        InertialMovingFrame(Vector3Cartesian origin, Vector3Cartesian speed)
        {
            _origin = origin;
            _speed = speed;
        }

        // transform to origin system
        VectorN<Real, 3> transf(const VectorN<Real, 3> &x, double t)
        {
            return x - _origin - _speed * t;
        }

    };

    class RotatingFrame : public MovingCoordSystem<3>
    {
        public:
        VectorN<double, 3> transf(const VectorN<double, 3> &x, double t)
        {
            return VectorN<double, 3>();
        }

        // os rotacije u odnosu na koordinatni sustav sredi8sta

        // transform to origin system
    };    

    class LorentzIntertialMovingFrame : public MovingCoordSystem<4>
    {
        public:
        VectorN<double, 3> transf(const VectorN<double, 3> &x, double t)
        {
            return VectorN<double, 3>();
        }

        // origin
        // brzina, odnosno boost

        // transform to origin system
    };

}

///////////////////////////   ./include/basic_types/Fields.h   ///////////////////////////


namespace MML
{
    static double InverseRadialFieldFuncCart(const VectorN<Real, 3> &x )   { return 1.0 / x.NormL2(); }
    static double InverseRadialFieldFuncSpher(const VectorN<Real, 3> &x )  { return 1.0 / x[0]; }
    static double InverseRadialFieldFuncCyl(const VectorN<Real, 3> &x )    { return 1.0 / sqrt(x[0]*x[0] + x[2]*x[2]); }

    class InverseRadialFieldCart : public IScalarFunction<3>
    {
    protected:
        double _constant;
    public:
        InverseRadialFieldCart() : _constant(-1.0) {}
        InverseRadialFieldCart(double constant) : _constant(constant) {}

        double operator()(const VectorN<double, 3> &x) const  { return _constant / x.NormL2(); }
    };

    class GravityPotentialFieldCart : public InverseRadialFieldCart
    {
    private:
        double _M;
        double _G;
    public:
        GravityPotentialFieldCart() : _M(1.0), _G(6.67408e-11), InverseRadialFieldCart(6.67408e-11) {}
        GravityPotentialFieldCart(double M) : _M(M), _G(6.67408e-11), InverseRadialFieldCart(M * 6.67408e-11) {}
        GravityPotentialFieldCart(double M, double G) : _M(M), _G(G), InverseRadialFieldCart(M * G) {}
    };
    
    class MultibodyGravityPotentialFieldCart : public IScalarFunction<3>
    {
    private:
        double _M;
        double _G;
        std::vector<double> _masses;
        std::vector<VectorN<double, 3>> _positions;
    public:
        MultibodyGravityPotentialFieldCart() : _M(1.0), _G(6.67408e-11) {}
        MultibodyGravityPotentialFieldCart(double M) : _M(M), _G(6.67408e-11) {}
        MultibodyGravityPotentialFieldCart(double M, double G) : _M(M), _G(G) {}

        double operator()(const VectorN<double, 3> &x) const  { return -_G * _M / x.NormL2(); }
    };      
    
    class CoulombPotentialFieldCart : public InverseRadialFieldCart
    {
    private:
        double _Q;
        double _C;
    public:
        CoulombPotentialFieldCart() : _Q(1.0), _C(6.67408e-11), InverseRadialFieldCart(8.987551787e9) {}
        CoulombPotentialFieldCart(double Q) : _Q(Q), _C(6.67408e-11), InverseRadialFieldCart(Q * 8.987551787e9) {}
        CoulombPotentialFieldCart(double Q, double C) : _Q(Q), _C(C), InverseRadialFieldCart(Q * C) {}
    };  
}

///////////////////////////   ./include/basic_types/Function.h   ///////////////////////////



namespace MML
{
    ///////////////////////////     REAL FUNCTION      ////////////////////////////////////
    class RealFunction : public IRealFunction
    {
        double (*_func)(const double) ;
    public:
        RealFunction(double (*inFunc)(const double) ) : _func(inFunc)    {}

        double operator()(const double x) const    { return _func(x); }

        // TODO - IsMonotone()? 
    };
    class RealFunctionFromStdFunc : public IRealFunction
    {
        std::function<double(const double)> _func;
    public:
        RealFunctionFromStdFunc(std::function<double(const double)> inFunc) : _func(inFunc)    {}

        double operator()(const double x) const    { return _func(x); }
    };

    ///////////////////////////     SCALAR FUNCTION       //////////////////////////////////
    template<int N>
    class ScalarFunction : public IScalarFunction<N>
    {
        double (*_func)(const VectorN<double, N> &);
    public:
        ScalarFunction( double (*inFunc)(const VectorN<double, N> &) ) : _func(inFunc)    {}

        double operator()(const VectorN<double, N> &x) const  { return _func(x); }
    };

    template<int N>
    class ScalarFunctionFromStdFunc : public IScalarFunction<N>
    {
        std::function<double(const VectorN<double, N> &)> _func;
    public:
        ScalarFunctionFromStdFunc(std::function<double(const VectorN<double, N> &)> inFunc) : _func(inFunc)     {}

        double operator()(const VectorN<double, N> &x) const  { return _func(x); }
    };
    
    /////////////////////////    VECTOR FUNCTION N -> N      ///////////////////////////////////
    template<int N>
    class VectorFunction : public IVectorFunction<N>
    {
        VectorN<double, N> (*_func)(const VectorN<double, N> &);
    public:
        VectorFunction( VectorN<double, N> (*inFunc)(const VectorN<double, N> &) ) : _func(inFunc)      {}

        VectorN<double, N> operator()(const VectorN<double, N> &x) const  { return _func(x); }
    };
    template<int N>
    class VectorFunctionFromStdFunc : public IVectorFunction<N>
    {
        std::function<VectorN<double, N>(const VectorN<double, N> &)> _func;
    public:
        VectorFunctionFromStdFunc(std::function<VectorN<double, N>(const VectorN<double, N> &)> &inFunc) : _func(inFunc)    {}

        VectorN<double, N> operator()(const VectorN<double, N> &x) const   { return _func(x); }
    };

   /////////////////////////    VECTOR FUNCTION N -> M      ///////////////////////////////////
    template<int N, int M>
    class VectorFunctionNM : public IVectorFunctionNM<N,M>
    {
        VectorN<double, M> (*_func)(const VectorN<double, N> &);
    public:
        VectorFunctionNM( VectorN<double, N> (*inFunc)(const VectorN<double, N> &) ) : _func(inFunc)      {}

        VectorN<double, M> operator()(const VectorN<double, N> &x) const  { return _func(x); }
    };
    template<int N, int M>
    class VectorFunctionNMFromStdFunc : public IVectorFunctionNM<N,M>
    {
        std::function<VectorN<double, M>(const VectorN<double, N> &)> _func;
    public:
        VectorFunctionNMFromStdFunc(std::function<VectorN<double, M>(const VectorN<double, N> &)> &inFunc) : _func(inFunc)    {}

        VectorN<double, M> operator()(const VectorN<double, N> &x) const   { return _func(x); }
    };

    //////////////////////     PARAMETRIC CURVE             ///////////////////////////////////
    template<int N>
    class ParametricCurve : public IParametricCurve<N>
    {
        VectorN<double, N> (*_func)(double);
    public:
        ParametricCurve( VectorN<double, N> (*inFunc)(double) ) : _func(inFunc)    {}

        virtual VectorN<double, N> operator()(double x) const  { return _func(x); }
    };

    template<int N>
    class ParametricCurveFromStdFunc : public IParametricCurve<N>
    {
        std::function<VectorN<double, N>(double)> _func;
    public:
        ParametricCurveFromStdFunc(std::function<VectorN<double, N>(double)> &inFunc) : _func(inFunc)    {}

        VectorN<double, N> operator()(double x) const   { return _func(x); }
    };

    /////////////////////       PARAMETRIC SURFACE         //////////////////////////////////
    template<int N>
    class ParametricSurface : public IParametricSurface<N>
    {
        VectorN<double, N> (*_func)(double u, double w);
    public:
        ParametricSurface( VectorN<double, N> (*inFunc)(double u, double w) ) : _func(inFunc)    {}

        VectorN<double, N> operator()(double u, double w) const  { return _func(u,w); }
    };

    template<int N>
    class ParametricSurfaceFromStdFunc : public IParametricSurface<N>
    {
        std::function<VectorN<double, N>(double u, double w)> _func;
    public:
        ParametricSurfaceFromStdFunc(std::function<VectorN<double, N>(double u, double w)> &inFunc) : _func(inFunc)    {}

        VectorN<double, N> operator()(double u, double w) const   { return _func(u,w); }
    };    

    /////////////////////       FUNCTION HELPERS         //////////////////////////////////
    class RealFuncDiffHelper : public IRealFunction
    {
        IRealFunction &_f1, &_f2;
    public:
        RealFuncDiffHelper(IRealFunction &f1, IRealFunction &f2) : _f1(f1), _f2(f2) {}
        double operator()(double x) const { return _f1(x) - _f2(x); }
    };
    
    class RealFuncDiffAbsHelper : public IRealFunction
    {
        IRealFunction &_f1, &_f2;
    public:
        RealFuncDiffAbsHelper(IRealFunction &f1, IRealFunction &f2) : _f1(f1), _f2(f2) {}
        double operator()(double x) const { return std::abs(_f1(x) - _f2(x)); }
    };

    class RealFuncDiffSqrHelper : public IRealFunction
    {
        IRealFunction &_f1, &_f2;
    public:
        RealFuncDiffSqrHelper(IRealFunction &f1, IRealFunction &f2) : _f1(f1), _f2(f2) {}
        double operator()(double x) const { return SQR(_f1(x) - _f2(x)); }
    };    

} // end namespace

///////////////////////////   ./include/basic_types/InterpolatedFunction.h   ///////////////////////////



namespace MML
{
    class Base_interp : public IRealFunction
    {
    public:
        mutable int jsav, cor;

        int n, mm, dj;
        const Real *xx, *yy;

        Base_interp(Vector<Real> &x, const Real *y, int m)
            : n((int)x.size()), mm(m), jsav(0), cor(0), xx(&x[0]), yy(y) 
        {
            dj = std::min(1,(int)pow((Real)n,0.25));
        }

        Real interp(Real x) const {
            int jlo = cor ? hunt(x) : locate(x);
            return rawinterp(jlo,x);
        }

        double operator()(Real x) const
        {
            return interp(x);
        }

        // Given a value x, return a value j such that x is (insofar as possible) centered in the subrange
        // xx[j..j+mm-1], where xx is the stored pointer. The values in xx must be monotonic, either
        // increasing or decreasing. The returned value is not less than 0, nor greater than n-1.
        int locate(const Real x) const
        {
            int ju,jm,jl;
            if (n < 2 || mm < 2 || mm > n) throw("locate size error");
            bool ascnd=(xx[n-1] >= xx[0]);
            jl=0;
            ju=n-1;
            while (ju-jl > 1) {
                jm = (ju+jl) >> 1;
                if (x >= xx[jm] == ascnd)
                    jl=jm;
                else
                    ju=jm;
            }
            cor = abs(jl-jsav) > dj ? 0 : 1;
            jsav = jl;
            return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));
        }

        // Given a value x, return a value j such that x is (insofar as possible) centered in the subrange
        // xx[j..j+mm-1], where xx is the stored pointer. The values in xx must be monotonic, either
        // increasing or decreasing. The returned value is not less than 0, nor greater than n-1.
        int hunt(const Real x) const
        {
            int jl=jsav, jm, ju, inc=1;
            if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
            bool ascnd=(xx[n-1] >= xx[0]);
            if (jl < 0 || jl > n-1) {
                jl=0;
                ju=n-1;
            } else {
                if (x >= xx[jl] == ascnd) {
                    for (;;) {
                        ju = jl + inc;
                        if (ju >= n-1) { ju = n-1; break;}
                        else if (x < xx[ju] == ascnd) break;
                        else {
                            jl = ju;
                            inc += inc;
                        }
                    }
                } else {
                    ju = jl;
                    for (;;) {
                        jl = jl - inc;
                        if (jl <= 0) { jl = 0; break;}
                        else if (x >= xx[jl] == ascnd) break;
                        else {
                            ju = jl;
                            inc += inc;
                        }
                    }
                }
            }
            while (ju-jl > 1) {
                jm = (ju+jl) >> 1;
                if (x >= xx[jm] == ascnd)
                    jl=jm;
                else
                    ju=jm;
            }
            cor = abs(jl-jsav) > dj ? 0 : 1;
            jsav = jl;
            return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));
        }
        
        Real virtual rawinterp(int jlo, Real x) const = 0;
    };
    
    struct LinearInterpRealFunc : Base_interp
    {
        LinearInterpRealFunc(Vector<Real> &xv, Vector<Real> &yv) : Base_interp(xv,&yv[0],2)  {}

        Real rawinterp(int j, Real x) const {
            if (xx[j]==xx[j+1]) return yy[j];
            else return yy[j] + ((x-xx[j])/(xx[j+1]-xx[j]))*(yy[j+1]-yy[j]);
        }
    };

    // Polynomial interpolation object. Construct with x and y vectors, and the number M of points
    // to be used locally (polynomial order plus one), then call interp for interpolated values.
    struct PolynomInterpRealFunc : Base_interp
    {
        mutable Real dy;

        // The user interface to Poly_interp is virtually the same as for Linear_interp
        // (end of ÷3.1), except that an additional argument in the constructor sets M, the number of points used (the order plus one). 
        PolynomInterpRealFunc(Vector<Real> &xv, Vector<Real> &yv, int m) : Base_interp(xv,&yv[0],m), dy(0.) 
        {}
        
        // Given a value x, and using pointers to data xx and yy, this routine returns an interpolated
        // value y, and stores an error estimate dy. The returned value is obtained by mm-point polynomial
        // interpolation on the subrange xx[jl..jl+mm-1].
        Real rawinterp(int jl, Real x) const
        {
            int i,m,ns=0;
            Real y,den,dif,dift,ho,hp,w;
            const Real *xa = &xx[jl], *ya = &yy[jl];
            Vector<Real> c(mm),d(mm);
            dif=abs(x-xa[0]);
            for (i=0;i<mm;i++) {
                if ((dift=abs(x-xa[i])) < dif) {
                    ns=i;
                    dif=dift;
                }
                c[i]=ya[i];
                d[i]=ya[i];
            }
            y=ya[ns--];
            for (m=1;m<mm;m++) {
                for (i=0;i<mm-m;i++) {
                    ho=xa[i]-x;
                    hp=xa[i+m]-x;
                    w=c[i+1]-d[i];
                    if ((den=ho-hp) == 0.0) throw("Poly_interp error");
                    den=w/den;
                    d[i]=hp*den;
                    c[i]=ho*den;
                }

                y += (dy=(2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--]));
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
    struct RationalInterpRealFunc : Base_interp
    {
        mutable Real dy;
        RationalInterpRealFunc(Vector<Real> &xv, Vector<Real> &yv, int m) : Base_interp(xv,&yv[0],m), dy(0.) 
        {}
        
        // Given a value x, and using pointers to data xx and yy, this routine returns an interpolated value
        // y, and stores an error estimate dy. The returned value is obtained by mm-point diagonal rational
        // function interpolation on the subrange xx[jl..jl+mm-1].        
        Real rawinterp(int jl, Real x) const
        {
            const Real TINY=1.0e-99;
            int m,i,ns=0;
            Real y,w,t,hh,h,dd;
            const Real *xa = &xx[jl], *ya = &yy[jl];
            Vector<Real> c(mm),d(mm);
            hh=abs(x-xa[0]);
            for (i=0;i<mm;i++) {
                h=abs(x-xa[i]);
                if (h == 0.0) {
                    dy=0.0;
                    return ya[i];
                } else if (h < hh) {
                    ns=i;
                    hh=h;
                }
                c[i]=ya[i];
                d[i]=ya[i]+TINY;
            }
            y=ya[ns--];
            for (m=1;m<mm;m++) {
                for (i=0;i<mm-m;i++) {
                    w=c[i+1]-d[i];
                    h=xa[i+m]-x;
                    t=(xa[i]-x)*d[i]/h;
                    dd=t-c[i+1];
                    if (dd == 0.0)  // This error condition indicates that the interpolating function has a pole at the requested value of x.
                        throw("Error in routine ratint");
                    dd=w/dd;
                    d[i]=c[i+1]*dd;
                    c[i]=t*dd;
                }
                y += (dy=(2*(ns+1) < (mm-m) ? c[ns+1] : d[ns--]));
            }
            return y;
        }
    };

    // Cubic spline interpolation object. Construct with x and y vectors, and (optionally) values of
    // the first derivative at the endpoints, then call interp for interpolated values.
    struct SplineInterpRealFunc : Base_interp
    {
        Vector<Real> y2;
        
        SplineInterpRealFunc(Vector<Real> &xv, Vector<Real> &yv, Real yp1=1.e99, Real ypn=1.e99)
            : Base_interp(xv,&yv[0],2), y2(xv.size())
        {
            sety2(&xv[0],&yv[0],yp1,ypn);
        }

        SplineInterpRealFunc(Vector<Real> &xv, const Real *yv, Real yp1=1.e99, Real ypn=1.e99)
            : Base_interp(xv,yv,2), y2(xv.size())
        {
            sety2(&xv[0],yv,yp1,ypn);
        }

        // This routine stores an array y2[0..n-1] with second derivatives of the interpolating function
        // at the tabulated points pointed to by xv, using function values pointed to by yv. If yp1 and/or
        // ypn are equal to 1  1099 or larger, the routine is signaled to set the corresponding boundary
        // condition for a natural spline, with zero second derivative on that boundary; otherwise, they are
        // the values of the first derivatives at the endpoints.
        void sety2(const Real *xv, const Real *yv, Real yp1, Real ypn)
        {
            int i,k;
            Real p,qn,sig,un;
            int n=(int) y2.size();
            Vector<Real> u(n-1);
            if (yp1 > 0.99e99)
                y2[0]=u[0]=0.0;
            else {
                y2[0] = -0.5;
                u[0]=(3.0/(xv[1]-xv[0]))*((yv[1]-yv[0])/(xv[1]-xv[0])-yp1);
            }
            for (i=1;i<n-1;i++) {
                sig=(xv[i]-xv[i-1])/(xv[i+1]-xv[i-1]);
                p=sig*y2[i-1]+2.0;
                y2[i]=(sig-1.0)/p;
                u[i]=(yv[i+1]-yv[i])/(xv[i+1]-xv[i]) - (yv[i]-yv[i-1])/(xv[i]-xv[i-1]);
                u[i]=(6.0*u[i]/(xv[i+1]-xv[i-1])-sig*u[i-1])/p;
            }
            if (ypn > 0.99e99)
                qn=un=0.0;
            else {
                qn=0.5;
                un=(3.0/(xv[n-1]-xv[n-2]))*(ypn-(yv[n-1]-yv[n-2])/(xv[n-1]-xv[n-2]));
            }
            y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
            for (k=n-2;k>=0;k--)
                y2[k]=y2[k]*y2[k+1]+u[k];
        }

        // Given a value x, and using pointers to data xx and yy, and the stored vector of second derivatives
        // y2, this routine returns the cubic spline interpolated value y.        
        Real rawinterp(int jl, Real x) const
        {
            int klo=jl,khi=jl+1;
            Real y,h,b,a;
            h=xx[khi]-xx[klo];
            if (h == 0.0) throw("Bad input to routine splint");
            a=(xx[khi]-x)/h;
            b=(x-xx[klo])/h;
            y=a*yy[klo]+b*yy[khi]+((a*a*a-a)*y2[klo]
                +(b*b*b-b)*y2[khi])*(h*h)/6.0;
            return y;
        }
    };

    // Barycentric rational interpolation object. After constructing the object, call interp for interpolated values. Note that no error estimate dy is calculated.
    struct BaryRatInterpRealFunc : Base_interp
    {
        Vector<Real> w;
        int d;

        // Constructor arguments are x and y vectors of length n, and order d of desired approximation.
        BaryRatInterpRealFunc(Vector<Real> &xv, Vector<Real> &yv, int dd) 
            : Base_interp(xv,&yv[0], (int) xv.size()), w(n), d(dd)
        {
            if (n<=d) throw("d too large for number of points in BaryRat_interp");
            for (int k=0;k<n;k++) {
                int imin=std::max(k-d,0);
                int imax = k >= n-d ? n-d-1 : k;
                Real temp = imin & 1 ? -1.0 : 1.0;
                Real sum=0.0;
                for (int i=imin;i<=imax;i++) {
                    int jmax=std::min(i+d,n-1);
                    Real term=1.0;
                    for (int j=i;j<=jmax;j++) {
                        if (j==k) continue;
                        term *= (xx[k]-xx[j]);
                    }
                    term=temp/term;
                    temp=-temp;
                    sum += term;
                }
                w[k]=sum;
            }
        }

        // Use equation (NR 3.4.9) to compute the barycentric rational interpolant. Note that jl is not used
        // since the approximation is global; it is included only for compatibility with Base_interp           
        Real rawinterp(int jl, Real x) const
        {
            Real num=0,den=0;
            for (int i=0;i<n;i++) {
                Real h=x-xx[i];
                if (h == 0.0) {
                    return yy[i];
                } else {
                    Real temp=w[i]/h;
                    num += temp*yy[i];
                    den += temp;
                }
            }
            return num/den;
        }

        // No need to invoke hunt or locate since the interpolation is global, so override interp to simply
        // call rawinterp directly with a dummy value of jl.        
        Real interp(Real x) {
            return rawinterp(1,x);
        }
    };   


    struct BilinInterpScalarFunction2D : public IScalarFunction<2>
    {
        int m,n;
        const Matrix<Real> &y;
        LinearInterpRealFunc x1terp, x2terp;

        BilinInterpScalarFunction2D(Vector<Real> &x1v, Vector<Real> &x2v, Matrix<Real> &ym)
            : m((int) x1v.size()), n( (int) x2v.size()), y(ym),
            x1terp(x1v,x1v), x2terp(x2v,x2v) {}

        double interp(double x1p, double x2p) const {
            int i,j;
            double yy, t, u;
            i = x1terp.cor ? x1terp.hunt(x1p) : x1terp.locate(x1p);
            j = x2terp.cor ? x2terp.hunt(x2p) : x2terp.locate(x2p);
            t = (x1p-x1terp.xx[i])/(x1terp.xx[i+1]-x1terp.xx[i]);
            u = (x2p-x2terp.xx[j])/(x2terp.xx[j+1]-x2terp.xx[j]);
            yy = (1.-t)*(1.-u)*y[i][j] + t*(1.-u)*y[i+1][j]
                + (1.-t)*u*y[i][j+1] + t*u*y[i+1][j+1];
            return yy;
        }

        double operator()(const VectorN<Real, 2> &x) const    
        { 
            return interp(x[0], x[1]); 
        }
    };
    
    struct PolynomInterpScalarFunction2D : public IScalarFunction<2>
    {
        int m,n,mm,nn;
        const Matrix<Real> &y;
        
        mutable Vector<Real> yv;
        mutable PolynomInterpRealFunc x1terp, x2terp;

        PolynomInterpScalarFunction2D(Vector<Real> &x1v, Vector<Real> &x2v, Matrix<Real> &ym,
            int mp, int np) : m((int) x1v.size()), n( (int) x2v.size()),
            mm(mp), nn(np), y(ym), yv(m),
            x1terp(x1v,yv,mm), x2terp(x2v,x2v,nn) {}

        double interp(double x1p, double x2p) const {
            int i,j,k;
            i = x1terp.cor ? x1terp.hunt(x1p) : x1terp.locate(x1p);
            j = x2terp.cor ? x2terp.hunt(x2p) : x2terp.locate(x2p);
            for (k=i;k<i+mm;k++) {
                x2terp.yy = &y[k][0];
                yv[k] = x2terp.rawinterp(j,x2p);
            }
            return x1terp.rawinterp(i,x1p);
        }

        double operator()(const VectorN<Real, 2> &x) const    
        { 
            return interp(x[0], x[1]); 
        }        
    };
    
    struct SplineInterpScalarFunction2D : public IScalarFunction<2>
    {
        int m,n;
        const Matrix<Real> &y;
        
        Vector<Real> &x1;
        mutable Vector<Real> yv;
        Vector<SplineInterpRealFunc*> srp;

        SplineInterpScalarFunction2D(Vector<Real> &x1v, Vector<Real> &x2v, Matrix<Real> &ym)
            : m((int) x1v.size()), n((int) x2v.size()), y(ym), yv(m), x1(x1v), srp(m) 
            {
            for (int i=0;i<m;i++) 
                srp[i] = new SplineInterpRealFunc(x2v,&y[i][0]);
        }

        ~SplineInterpScalarFunction2D(){
            for (int i=0;i<m;i++) delete srp[i];
        }

        double interp(double x1p, double x2p) const 
        {
            for (int i=0;i<m;i++) 
                yv[i] = (*srp[i]).interp(x2p);
            
            SplineInterpRealFunc scol(x1,yv);
            
            return scol.interp(x1p);
        }

        double operator()(const VectorN<Real, 2> &x) const    
        { 
            return interp(x[0], x[1]); 
        }        
    };
    
    class InterpolatedScalarFunction3D : public IScalarFunction<3>
    {
        public:
        InterpolatedScalarFunction3D() {}

        double operator()(const VectorN<Real, 3> &x) const    { return 0.0; }
        virtual double operator()(double u, double w, double z)
        {
            VectorN<Real, 3> coord{u,w,z};

            return operator()(coord);
        }            
    };

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

        SplineInterpParametricCurve(Matrix<Real> &ptsin, bool close=0)
        : n(ptsin.RowNum()), dim(ptsin.ColNum()), in(close ? 2*n : n),
        cls(close), pts(dim,in), s(in), ans(dim), srp(dim) 
        {
            int i,ii,im,j,ofs;
            Real ss,soff,db,de;
            ofs = close ? n/2 : 0;
            s[0] = 0.;
            for (i=0;i<in;i++) {
                ii = (i-ofs+n) % n;
                im = (ii-1+n) % n;
                for (j=0;j<dim;j++) pts[j][i] = ptsin[ii][j];
                if (i>0) {
                    s[i] = s[i-1] + rad(&ptsin[ii][0],&ptsin[im][0]);
                    if (s[i] == s[i-1]) throw("error in Curve_interp");
                }
            }
            ss = close ? s[ofs+n]-s[ofs] : s[n-1]-s[0];
            soff = s[ofs];
            for (i=0;i<in;i++) s[i] = (s[i]-soff)/ss;
            for (j=0;j<dim;j++) {
                db = in < 4 ? 1.e99 : fprime(&s[0],&pts[j][0],1);
                de = in < 4 ? 1.e99 : fprime(&s[in-1],&pts[j][in-1],-1);
                srp[j] = new SplineInterpRealFunc(s,&pts[j][0],db,de);
            }
        }
        ~SplineInterpParametricCurve() {
            for (int j=0;j<dim;j++) delete srp[j];
        }
        
        VectorN<Real, N> &interp(Real t) const 
        {
            VectorN<Real, N> ans;

            if (cls) 
                t = t - floor(t);
            for (int j=0;j<dim;j++) 
                ans[j] = (*srp[j]).interp(t);
            
            return ans;
        }

        VectorN<Real, N> operator()(Real t) const    
        { 
            return interp(t); 
        }

        Real fprime(Real *x, Real *y, int pm) {
            Real s1 = x[0]-x[pm*1], s2 = x[0]-x[pm*2], s3 = x[0]-x[pm*3],
                s12 = s1-s2, s13 = s1-s3, s23 = s2-s3;
            return -(s1*s2/(s13*s23*s3))*y[pm*3]+(s1*s3/(s12*s2*s23))*y[pm*2]
                -(s2*s3/(s1*s12*s13))*y[pm*1]+(1./s1+1./s2+1./s3)*y[0];
        }

        Real rad(const Real *p1, const Real *p2) {
            Real sum = 0.;
            for (int i=0;i<dim;i++) 
                sum += SQR(p1[i]-p2[i]);
            return sqrt(sum);
        }

    };
    
    template<int N>
    class InterpolatedSurface : public IParametricSurface<N>
    {
        public:
        InterpolatedSurface() {}

        VectorN<Real, N> operator()(const VectorN<Real, 2> &x) const    { return VectorN<Real, N>{}; }
    };       
}

///////////////////////////   ./include/basic_types/Functionals.h   ///////////////////////////


namespace MML
{
    template <int _Dim, typename _Field = Real>
    class LinearFunctionalN
    {
    private:
        VectorN<_Field, _Dim> _vecCoef;
    public:
        LinearFunctionalN() {}
        LinearFunctionalN(const VectorN<_Field, _Dim> &vecCoef) : _vecCoef(vecCoef) {}
        LinearFunctionalN(std::initializer_list<_Field> list) : _vecCoef(list) {}

        LinearFunctionalN(const LinearFunctionalN &Copy) : _vecCoef(Copy._vecCoef) {}
        ~LinearFunctionalN() {}

        LinearFunctionalN& operator=(const LinearFunctionalN &Copy) { _vecCoef = Copy._vecCoef; return *this; }

        _Field operator()(const VectorN<_Field, _Dim> &vecX) const
        {
            _Field result = 0.0;
            for (int i = 0; i < _Dim; i++)
                result += _vecCoef[i] * vecX[i];
            return result;
        }

        LinearFunctionalN operator+(const LinearFunctionalN &b) const
        {
            LinearFunctionalN result;
            for (int i = 0; i < _Dim; i++)
            {
                if (i < _vecCoef.size())
                    result._vecCoef[i] += _vecCoef[i];
                if (i < b._vecCoef.size())
                    result._vecCoef[i] += b._vecCoef[i];
            }
            return result;
        }

        LinearFunctionalN operator-(const LinearFunctionalN &b) const
        {
            LinearFunctionalN result;
            for (int i = 0; i < _Dim; i++)
            {
                if (i < _vecCoef.size())
                    result._vecCoef[i] += _vecCoef[i];
                if (i < b._vecCoef.size())
                    result._vecCoef[i] -= b._vecCoef[i];
            }
            return result;
        }

        LinearFunctionalN operator*(double b) const
        {
            LinearFunctionalN result;
            for (int i = 0; i < _vecCoef.size();    i++)
                result._vecCoef[i] = _vecCoef[i] * b;
            return result;
        }

    };

    class LinearFunctional
    {
    private:
        std::vector<double> _vecCoef;
    public:
        LinearFunctional() {}
        LinearFunctional(const std::vector<double> &vecCoef) : _vecCoef(vecCoef) {}
        LinearFunctional(std::initializer_list<double> list) : _vecCoef(list) {}

        LinearFunctional(const LinearFunctional &Copy) : _vecCoef(Copy._vecCoef) {}
        ~LinearFunctional() {}

        LinearFunctional& operator=(const LinearFunctional &Copy) { _vecCoef = Copy._vecCoef; return *this; }

        double operator()(const std::vector<double> &vecX) const
        {
            double result = 0.0;
            for (int i = 0; i < _vecCoef.size(); i++)
                result += _vecCoef[i] * vecX[i];
            return result;
        }

        LinearFunctional operator+(const LinearFunctional &b) const
        {
            LinearFunctional result;
            int n = (int) std::max(_vecCoef.size(), b._vecCoef.size());
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

        LinearFunctional operator-(const LinearFunctional &b) const
        {
            LinearFunctional result;
            int n = (int) std::max(_vecCoef.size(), b._vecCoef.size());
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

        LinearFunctional operator*(double b) const
        {
            LinearFunctional result;
            result._vecCoef.resize(_vecCoef.size());
            for (int i = 0; i < _vecCoef.size(); i++)
                result._vecCoef[i] = _vecCoef[i] * b;
            return result;
        }
    };
}

///////////////////////////   ./include/basic_types/Functions.h   ///////////////////////////



namespace MML
{
    namespace StdFunctions
    {
        static inline Real Sin(Real x) { return sin(x); }
        static inline Real Cos(Real x) { return cos(x); }
        static inline Real Tan(Real x) { return tan(x); }
        static inline Real Exp(Real x) { return exp(x); }
        static inline Real Log(Real x) { return log(x); }
        static inline Real Sqrt(Real x) { return sqrt(x); }
        static inline Real Pow(Real x, Real y) { return pow(x, y); }
        static inline Real Sinh(Real x) { return sinh(x); }
        static inline Real Cosh(Real x) { return cosh(x); }
        static inline Real Tanh(Real x) { return tanh(x); }
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
        
        static inline Complex Sin(Complex x) { return sin(x); }
        static inline Complex Cos(Complex x) { return cos(x); }
        static inline Complex Tan(Complex x) { return tan(x); }
        static inline Complex Exp(Complex x) { return exp(x); }
        static inline Complex Log(Complex x) { return log(x); }
        static inline Complex Sqrt(Complex x) { return sqrt(x); }
        static inline Complex Pow(Complex x, Complex y) { return pow(x, y); }
        static inline Complex Sinh(Complex x) { return sinh(x); }
        static inline Complex Cosh(Complex x) { return cosh(x); }
        static inline Complex Tanh(Complex x) { return tanh(x); }
        static inline Complex Asin(Complex x) { return asin(x); }
        static inline Complex Acos(Complex x) { return acos(x); }
        static inline Complex Atan(Complex x) { return atan(x); }
        static inline Complex Asinh(Complex x) { return asinh(x); }
        static inline Complex Acosh(Complex x) { return acosh(x); }
        static inline Complex Atanh(Complex x) { return atanh(x); }        

    }
}

///////////////////////////   ./include/basic_types/Curves.h   ///////////////////////////



namespace MML
{
    // TODO - dodati još par planarnih i prostornih krivulja
    // TODO - dodati polarne krivulje r = r(phi)
    namespace Curves
    {
        ////////////////////////////////             PLANAR CURVES                  //////////////////////////////////
        class Circle2DCurve : public IParametricCurve<2> 
        {
            Real _radius;
        public:
            Circle2DCurve() : _radius(1) {}
            Circle2DCurve(Real radius) : _radius(radius) {}

            VectorN<Real, 2> operator()(double t) const  { return MML::VectorN<Real, 2>{_radius * cos(t), _radius * sin(t)}; }
        };

        class LogSpiralCurve : public IParametricCurve<2> 
        {
            Real _lambda, _c;       // lambda < 0;   c != 0
        public:
            LogSpiralCurve() : _lambda(-1), _c(1) {}
            LogSpiralCurve(Real lambda, Real c) : _lambda(lambda), _c(c) {}

            VectorN<Real, 2> operator()(double t) const  { return MML::VectorN<Real, 2>{exp(_lambda * t) * cos(t), exp(_lambda * t) * sin(t)}; }
        };

        class LemniscateCurve : public IParametricCurve<2> 
        {
        public:
            VectorN<Real, 2> operator()(double t) const  { return MML::VectorN<Real, 2>{cos(t) / (1 + sin(t)*sin(t)), sin(t) * cos(t) / (1 + sin(t)*sin(t))}; }
        };    

        class DeltoidCurve : public IParametricCurve<2> 
        {
            int _n;
        public:
            DeltoidCurve() : _n(1) {}
            DeltoidCurve(int n) : _n(n) {}

            VectorN<Real, 2> operator()(double t) const  { return MML::VectorN<Real, 2>{2 * _n * cos(t) * (1 + cos(t)), 2 * _n * sin(t) * (1 - cos(t))}; }
        };

        class AstroidCurve : public IParametricCurve<2> 
        {
            Real _c;            // c > 0
        public:
            AstroidCurve() : _c(1) {}
            AstroidCurve(Real c) : _c(c) {}

            VectorN<Real, 2> operator()(double t) const  { return MML::VectorN<Real, 2>{_c * cos(t)* cos(t)* cos(t), _c * sin(t)* sin(t)* sin(t)}; }
        };

        class EpitrochoidCurve : public IParametricCurve<2> 
        {
            Real _radius, _c;
            int _n;
        public:
            EpitrochoidCurve() : _radius(1), _c(1), _n(1) {}
            EpitrochoidCurve(Real radius, Real c, int n) : _radius(radius), _c(c), _n(n) {}

            VectorN<Real, 2> operator()(double t) const  { return MML::VectorN<Real, 2>{cos(t) - _c * cos(_n*t), sin(t) - _c * sin(_n*t) }; }
        };

        /////////////////////////////             PLANAR POLAR CURVES                  ////////////////////////////////

        /////////////////////////////////             SPACE CURVES                  ///////////////////////////////////
        class Circle3DXY : public IParametricCurve<3> {
            Real _radius;
        public:
            Circle3DXY() : _radius(1) {}
            Circle3DXY(Real radius) : _radius(radius) {}

            VectorN<Real, 3> operator()(double t) const  { return MML::VectorN<Real, 3>{_radius * cos(t), _radius * sin(t), 0}; }
        };
        class Circle3DXZ : public IParametricCurve<3> {
            Real _radius;
        public:
            Circle3DXZ() : _radius(1) {}
            Circle3DXZ(Real radius) : _radius(radius) {}

            VectorN<Real, 3> operator()(double t) const  { return MML::VectorN<Real, 3>{_radius * cos(t), 0, _radius * sin(t)}; }
        };
        class Circle3DYZ : public IParametricCurve<3> {
            Real _radius;
        public:
            Circle3DYZ() : _radius(1) {}
            Circle3DYZ(Real radius) : _radius(radius) {}
            
            VectorN<Real, 3> operator()(double t) const  { return MML::VectorN<Real, 3>{0, _radius * cos(t), _radius * sin(t)}; }
        };
        
        class HelixCurve : public IParametricCurve<3> 
        {
            Real _radius, _b;
        public:
            HelixCurve() : _radius(1.0), _b(1.0) {}
            HelixCurve(Real radius, Real b) : _radius(radius), _b(b) {}
    
            VectorN<Real, 3> operator()(double t) const  { return MML::VectorN<Real, 3>{_radius * cos(t), _radius * sin(t), _b * t}; }

            double getCurvature(double t) const  { return _radius / (SQR(_radius) + SQR(_b)); }
            double getTorsion(double t) const    { return _b / (SQR(_radius) + SQR(_b)); }
        };

        class TwistedCubicCurve : public IParametricCurve<3> 
        {
        public:
            VectorN<Real, 3> operator()(double t) const  { return MML::VectorN<Real, 3>{t, t * t, t * t * t}; }
        };

        class ToroidalSpiralCurve : public IParametricCurve<3> 
        {
            int _n;
        public:
            ToroidalSpiralCurve() : _n(1) {}
            ToroidalSpiralCurve(int n) : _n(n) {}

            VectorN<Real, 3> operator()(double t) const  { return MML::VectorN<Real, 3>{(4 + sin(_n*t)) * cos(t), (4 + sin(_n*t)) * sin(t), cos(_n*t)}; }
        };
    }
}

///////////////////////////   ./include/basic_types/Surfaces.h   ///////////////////////////


namespace MML
{
    namespace Surfaces
    {
        static ParametricSurface<3> test1([](double u, double w) { return MML::VectorN<Real, 3>{u, w, u+w}; });
    }
}

///////////////////////////   ./include/basic_types/DeltaFunction.h   ///////////////////////////

namespace MML
{

    class DeltaFunction
    {
        double _value;
        double _eps;
        // zajedno odredjuju sirinu i visinu pravokutnika
    };
}

///////////////////////////   ./include/basic_types/ODESystem.h   ///////////////////////////



namespace MML
{
	class ODESystem : public IODESystem
	{
    protected:
        int _dim;
        void (*_func)(double, const Vector<Real>&, Vector<Real> &);

	public:
        ODESystem() : _dim(0), _func(nullptr) { }
        ODESystem(int n, void (*inFunc)(double, const Vector<Real>&, Vector<Real> &)) : _dim(n), _func(inFunc) { }
        
        int getDim() { return _dim; }
		void derivs(const double t, const Vector<Real> &x, Vector<Real> &dxdt)
        {
            _func(t, x, dxdt);
        }
	};

	class ODESystemWithJacobian : public ODESystem
	{
    private:
        int _dim;
        void (*_funcJac)(const double, const Vector<Real>&, Vector<Real> &, Matrix<Real> &);

	public:
        ODESystemWithJacobian() { }
        ODESystemWithJacobian(int n, 
                              void (*inFunc)(double, const Vector<Real>&, Vector<Real> &),
                              void (*inFuncJac)(const double t, const Vector<Real> &x, Vector<Real> &dxdt, Matrix<Real> &dydx)
                              ) : ODESystem(n, inFunc), _funcJac(inFuncJac) { }
        
        void jacobian(const double t, Vector<Real> &x, Vector<Real> &dxdt, Matrix<Real> &dydx)
        {
            _funcJac(t, x, dxdt, dydx);
        }
	};

    class ODESystemSolution
    {
        // first values are initial conditions
    public:
        int  _sys_dim;
        int  _count;
        Real x1,x2;
        Vector<Real> xval;
        Matrix<Real> yval;
        
        ODESystemSolution() {}
        ODESystemSolution(int dim, int maxSteps) : _sys_dim(dim)
        {
            xval.Resize(maxSteps+1);
            yval.Resize(dim, maxSteps+1);
        }
    };

    class ODESystemSolutionEqualSpacing
    {
        // first values are initial conditions
    public:
        int  _sys_dim;
        int  _count;
        Real x1,x2;
        Vector<Real> xval;
        Matrix<Real> yval;
        
        ODESystemSolutionEqualSpacing() {}
        ODESystemSolutionEqualSpacing(int dim, int numSteps) : _sys_dim(dim), _count(numSteps+1)
        {
            xval.Resize(numSteps+1);
            yval.Resize(dim, numSteps+1);
        }
        // InterpRealFunctionLinear getSolution(int component)
        // {
        //     // na osnovu xsave i ysave napravi interpoliranu funkciju

        //     return InterpRealFunctionLinear();
        // }

    };
}
///////////////////////////   ./include/algorithms/Derivation.h   ///////////////////////////




namespace MML
{
    // TODO - finish 2nd deriv for Scalar function, 2, 4, 6, and 8th order
    class Derivation
    {
    public:
        /********************************************************************************************************************/
        /********                               Numerical derivatives of FIRST order                                 ********/
        /********************************************************************************************************************/
        //////////////////////////              RealFunction            //////////////////////////
        static Real NDer1(const MML::IRealFunction &f, Real x, Real* error = nullptr)
        {
            // Error bound ~eps^1/2
            // Note that this estimate of h differs from the best estimate by a factor of sqrt((|f(x)| + |f(x+h)|)/|f''(x)|).
            // Since this factor is invariant under the scaling f -> kf, then we are somewhat justified in approximating it by 1.
            // This approximation will get better as we move to higher orders of accuracy.
            const Real h = 2 * std::sqrt(Constants::Epsilon);

            return NDer1(f, x, h, error);
        }

        static Real NDer1(const MML::IRealFunction &f, Real x, Real h, Real* error = nullptr)
        {
            Real yh   = f(x + h);
            Real y0   = f(x);
            Real diff = yh - y0;
            if (error)
            {
                Real ym   = f(x - h);
                Real ypph = std::abs(yh - 2 * y0 + ym) / h;

                // h*|f''(x)|*0.5 + (|f(x+h)+|f(x)|)*eps/h
                *error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
            }
            return diff / h;
        }

        static Real NSecDer1(const MML::IRealFunction &f, Real x, Real* error = nullptr)
        {
            Real h = 2 * std::sqrt(Constants::Epsilon);

            return NSecDer1(f, x, h, error);
        }

        static Real NSecDer1(const MML::IRealFunction &f, Real x, Real h, Real* error = nullptr)
        {
            Real yh   = NDer2(f, x + h, h, error);
            Real y0   = NDer2(f, x, h, error);
            Real diff = yh - y0;
            if (error)
            {
                Real ym   = NDer2(f, x - h, h, error);
                Real ypph = std::abs(yh - 2 * y0 + ym) / h;

                *error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
            }
            return diff / h;
        }
        
        static Real NThirdDer1(const MML::IRealFunction &f, Real x, Real* error = nullptr)
        {
            Real h = 2 * std::sqrt(Constants::Epsilon);

            return NThirdDer1(f, x, h, error);
        }

        static Real NThirdDer1(const MML::IRealFunction &f, Real x, Real h, Real* error = nullptr)
        {
            Real yh   = NSecDer2(f, x + h, h, error);
            Real y0   = NSecDer2(f, x, h, error);
            Real diff = yh - y0;
            if (error)
            {
                Real ym   = NSecDer2(f, x - h, h, error);
                Real ypph = std::abs(yh - 2 * y0 + ym) / h;

                *error = ypph / 2 + (std::abs(yh) + std::abs(y0)) * Constants::Epsilon / h;
            }
            return diff / h;
        }
        
        //////////////////////////             ScalarFunction           //////////////////////////
        template <int N>
        static Real NDer1Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, Real *error = nullptr)
        {
            Real h = 2 * sqrt(std::numeric_limits<Real>::epsilon());

            return NDer1Partial(f, deriv_index, point, h, error);
        }

        template <int N>
        static Real NDer1Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, Real h, Real *error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();

            Real     orig_x     = point[deriv_index];

            VectorN<Real, N> x  = point;
            Real y0 = f(x);

            x[deriv_index] = orig_x + h;
            Real yh = f(x);

            Real diff = yh - y0;
            if (error)
            {
                x[deriv_index] = orig_x - h;
                Real ym   = f(x);
                Real ypph = std::abs(yh - 2 * y0 + ym) / h;
                *error = ypph / 2 + (std::abs(yh) + std::abs(y0))*eps / h;
            }
            return diff / h;
        } 

        template <int N>
        static Real NSecDer1Partial(const IScalarFunction<N> &f, int der_ind1, int der_ind2, const VectorN<Real, N> &point, Real *error = nullptr)
        {
            Real h = 2 * sqrt(std::numeric_limits<Real>::epsilon());

            return NSecDer1Partial(f, der_ind1, der_ind2, point, h, error);
        }

        template <int N>
        static Real NSecDer1Partial(const IScalarFunction<N> &f, int der_ind1, int der_ind2, const VectorN<Real, N> &point, Real h, Real *error = nullptr)
        {
            Real x_orig_val      = point[der_ind2];

            auto x_eval_pos      = point;
            Real y0              = NDer2Partial(f, der_ind1, x_eval_pos, error);
            x_eval_pos[der_ind2] = x_orig_val + h;
            Real yh              = NDer2Partial(f, der_ind1, x_eval_pos, error);

            Real diff = yh - y0;
            if (error)
            {
                x_eval_pos[der_ind2] = x_orig_val - h;
                
                Real ym   = NDer2Partial(f, der_ind1, x_eval_pos, error);
                Real ypph = std::abs(yh - 2 * y0 + ym) / h;
                
                *error    = ypph / 2 + (std::abs(yh) + std::abs(y0)) * std::numeric_limits<Real>::epsilon() / h;
            }
            return diff / h;
        } 
        
        template <int N>
        static VectorN<Real, N> NDer1PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            Real h = 2 * sqrt(std::numeric_limits<Real>::epsilon());

            return NDer1PartialByAll(f, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer1PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, Real h, VectorN<Real, N> *error = nullptr)
        {
            VectorN<Real, N> ret;

            for( int i=0; i<N; i++)
            {
                if( error )
                    ret[i] = NDer1Partial(f, i, point, h, &(*error)[i]);
                else
                    ret[i] = NDer1Partial(f, i, point, h);
            }

            return ret;
        }

        //////////////////////////             VectorFunction           //////////////////////////
        template <int N>
        static Real NDer1Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, Real *error = nullptr)
        {
            const Real h = 2 * sqrt(std::numeric_limits<Real>::epsilon());

            return NDer1Partial(f, func_index, deriv_index, point, h, error);
        }

        template <int N>
        static Real NDer1Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, Real h, Real *error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();

            auto   x      = point;

            Real x_orig = x[deriv_index];
            Real y0     = f(x)[func_index];

            x[deriv_index] = x_orig + h;
            Real yh      = f(x)[func_index];

            Real diff = yh - y0;
            if (error)
            {
                x[deriv_index] = x_orig - h;
                Real ym      = f(x)[func_index];
                Real ypph    = std::abs(yh - 2 * y0 + ym) / h;
                *error = ypph / 2 + (std::abs(yh) + std::abs(y0))*eps / h;
            }
            return diff / h;
        }

        template <int N>
        static VectorN<Real, N> NDer1PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            Real h = 2 * sqrt(std::numeric_limits<Real>::epsilon());

            return NDer1PartialByAll(f, func_index, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer1PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, Real h, VectorN<Real, N> *error = nullptr)
        {
            VectorN<Real, N> ret;

            for( int i=0; i<N; i++)
            {
                if( error )
                    ret[i] = NDer1Partial(f, func_index, i, point, h, &(*error)[i]);
                else
                    ret[i] = NDer1Partial(f, func_index, i, point, h);
            }

            return ret;         
        }

        template <int N>
        static MatrixNM<Real, N,N> NDer1PartialAllByAll(const IVectorFunction<N> &f, const VectorN<Real, N> &point, MatrixNM<Real, N,N> *error = nullptr)
        {
            Real h = 2 * sqrt(std::numeric_limits<Real>::epsilon());

            return NDer1PartialAllByAll(f, point, h, error);
        }

        template <int N>
        static MatrixNM<Real, N,N> NDer1PartialAllByAll(const IVectorFunction<N> &f, const VectorN<Real, N> &point, Real h, MatrixNM<Real, N,N> *error = nullptr)
        {
            MatrixNM<Real, N,N> ret;

            for( int i=0; i<N; i++)
                for( int j=0; j<N; j++)
                {
                    if( error )
                        ret(i,j) = NDer1Partial(f, i, j, point, h, &((*error)(i,j)));
                    else
                        ret(i,j) = NDer1Partial(f, i, j, point, h);
                }

            return ret;         
        }

        /////////////////////////             ParametricCurve           /////////////////////////
        template <int N>
        static VectorN<Real, N> NDer1(const MML::IParametricCurve<N> &f, Real t, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();
            Real h = 2 * std::sqrt(eps);

            return NDer1(f, t, h, error);
        }        

        template <int N>
        static VectorN<Real, N> NDer1(const MML::IParametricCurve<N> &f, Real t, Real h, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();

            VectorN<Real, N> yh   = f(t + h);
            VectorN<Real, N> y0   = f(t);
            VectorN<Real, N> diff = yh - y0;
            
            if (error)
            {
                VectorN<Real, N> ym = f(t - h);
                VectorN<Real, N> ypph_vec = yh - 2 * y0 + ym;

                Real ypph = ypph_vec.NormL2() / h;

                *error = ypph / 2 + (yh.NormL2() + y0.NormL2()) * eps / h;
            }
            return diff / h;
        }

        template <int N>        
        static VectorN<Real, N> NSecDer1(const MML::IParametricCurve<N> &f, Real x, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();
            Real h = 2 * std::sqrt(eps);

            return NSecDer1(f, x, h, error);
        }

        template <int N>
        static VectorN<Real, N> NSecDer1(const MML::IParametricCurve<N> &f, Real x, Real h, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();

            VectorN<Real, N>  yh   = NDer2(f, x + h, h, error);
            VectorN<Real, N>  y0   = NDer2(f, x, h, error);
            VectorN<Real, N>  diff = yh - y0;
            if (error)
            {
                VectorN<Real, N> ym       = NDer2(f, x - h, h, error);
                VectorN<Real, N> ypph_vec = (yh - 2 * y0 + ym) / h;
                
                Real ypph = ypph_vec.NormL2();

                *error = ypph / 2 + (yh.NormL2() + y0.NormL2()) * eps / h;
            }
            return diff / h;
        }

        template <int N>        
        static VectorN<Real, N> NThirdDer1(const MML::IParametricCurve<N> &f, Real x, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();
            Real h = 2 * std::sqrt(eps);

            return NThirdDer1(f, x, h, error);
        }

        template <int N>
        static VectorN<Real, N> NThirdDer1(const MML::IParametricCurve<N> &f, Real x, Real h, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();

            VectorN<Real, N>  yh   = NSecDer2(f, x + h, h, error);
            VectorN<Real, N>  y0   = NSecDer2(f, x, h, error);
            VectorN<Real, N>  diff = yh - y0;
            if (error)
            {
                VectorN<Real, N> ym       = NSecDer2(f, x - h, h, error);
                VectorN<Real, N> ypph_vec = (yh - 2 * y0 + ym) / h;
                
                Real ypph = ypph_vec.NormL2();

                *error = ypph / 2 + (yh.NormL2() + y0.NormL2()) * eps / h;
            }
            return diff / h;
        }        
        
        /********************************************************************************************************************/
        /********                               Numerical derivatives of SECOND order                                 ********/
        /********************************************************************************************************************/
        
        //////////////////////////              RealFunction            //////////////////////////
        static Real NDer2(const MML::IRealFunction &f, Real x, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            // Error bound ~eps^2/3
            // See the previous discussion to understand determination of h and the error bound.
            // Series[(f[x+h] - f[x-h])/(2*h), {h, 0, 4}]
            Real h = std::pow(3 * eps, static_cast<Real>(1) / static_cast<Real>(3));
            //h = detail::make_xph_representable(x, h);

            return NDer2(f, x, h, error);
        } 

        static Real NDer2(const MML::IRealFunction &f, Real x, Real h, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            Real yh   = f(x + h);
            Real ymh  = f(x - h);
            Real diff = yh - ymh;
            if (error)
            {
                Real y2h  = f(x + 2 * h);
                Real ym2h = f(x - 2 * h);
                *error = eps * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
            }

            return diff / (2 * h);
        }

        static Real NSecDer2(const MML::IRealFunction &f, Real x, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(3 * eps, static_cast<Real>(1) / static_cast<Real>(3));

            return NSecDer2(f, x, h, error);
        }

        static Real NSecDer2(const MML::IRealFunction &f, Real x, Real h, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            Real yh   = NDer4(f, x + h, error);
            Real ymh  = NDer4(f, x - h, error);
            Real diff = yh - ymh;
            if (error)
            {
                Real y2h   = NDer4(f, x + 2 * h, error);
                Real ym2h  = NDer4(f, x - 2 * h, error);

                *error = eps * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
            }

            return diff / (2 * h);
        }          

        static Real NThirdDer2(const MML::IRealFunction &f, Real x, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(3 * eps, static_cast<Real>(1) / static_cast<Real>(3));

            return NThirdDer2(f, x, h, error);
        }

        static Real NThirdDer2(const MML::IRealFunction &f, Real x, Real h, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            Real yh   = NSecDer4(f, x + h, error);
            Real ymh  = NSecDer4(f, x - h, error);
            Real diff = yh - ymh;
            if (error)
            {
                Real y2h   = NSecDer4(f, x + 2 * h, error);
                Real ym2h  = NSecDer4(f, x - 2 * h, error);

                *error = eps * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
            }

            return diff / (2 * h);
        }        
        
        //////////////////////////             ScalarFunction           //////////////////////////
        template <int N>
        static Real NDer2Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(3 * eps, static_cast<Real>(1) / static_cast<Real>(3));

            return NDer2Partial(f, deriv_index, point, h, error);
        }

        template <int N>
        static Real NDer2Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, Real h, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

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
                *error = eps * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
            }

            return diff / (2 * h);
        }

        template <int N>
        static Real NSecDer2Partial(const IScalarFunction<N> &f, int der_ind1, int der_ind2, const VectorN<Real, N> &point, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(3 * eps, static_cast<Real>(1) / static_cast<Real>(3));

            return NSecDer2Partial(f, der_ind1, der_ind2, point, h, error);
        }

        template <int N>
        static Real NSecDer2Partial(const IScalarFunction<N> &f, int der_ind1, int der_ind2, const VectorN<Real, N> &point, Real h, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            Real orig_x          = point[der_ind2];
            auto x_eval_pos      = point;
            
            x_eval_pos[der_ind2] = orig_x + h;
            Real yh              = NDer4Partial(f, der_ind1, x_eval_pos, error);

            x_eval_pos[der_ind2] = orig_x - h;
            Real ymh             = NDer4Partial(f, der_ind1, x_eval_pos, error);

            Real diff            = yh - ymh;

            if (error)
            {
                x_eval_pos[der_ind2] = orig_x + 2 * h;
                Real y2h             = NDer4Partial(f, der_ind1, x_eval_pos, error);

                x_eval_pos[der_ind2] = orig_x - 2 * h;
                Real ym2h            = NDer4Partial(f, der_ind1, x_eval_pos, error);
                
                *error = eps * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
            }

            return diff / (2 * h);
        }

        template <int N>
        static VectorN<Real, N> NDer2PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(3 * eps, static_cast<Real>(1) / static_cast<Real>(3));

            return NDer2PartialByAll(f, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer2PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, Real h, VectorN<Real, N> *error = nullptr)
        {
            VectorN<Real, N> ret;

            for( int i=0; i<N; i++)
            {
                if( error )
                    ret[i] = NDer2Partial(f, i, point, h, &(*error)[i]);
                else
                    ret[i] = NDer2Partial(f, i, point, h);
            }

            return ret;
        }

        //////////////////////////             VectorFunction           //////////////////////////
        template <int N>
        static Real NDer2Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(3 * eps, static_cast<Real>(1) / static_cast<Real>(3));

            return NDer2Partial(f, func_index, deriv_index, point, h, error);
        }        

        template <int N>
        static Real NDer2Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, Real h, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            Real     orig_x = point[deriv_index];

            VectorN<Real, N> x{point};
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

                *error = eps * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
            }

            return diff / (2 * h);
        }

        template <int N>
        static VectorN<Real, N> NDer2PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(3 * eps, static_cast<Real>(1) / static_cast<Real>(3));

            return NDer2PartialByAll(f, func_index, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer2PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, Real h, VectorN<Real, N> *error = nullptr)
        {
            VectorN<Real, N> ret;

            for( int i=0; i<N; i++)
            {
                if( error )
                    ret[i] = NDer2Partial(f, func_index, i, point, h, &(*error)[i]);
                else
                    ret[i] = NDer2Partial(f, func_index, i, point, h);
            }

            return ret;         
        }

        template <int N>
        static MatrixNM<Real, N,N> NDer2PartialAllByAll(const IVectorFunction<N> &f, const VectorN<Real, N> &point, MatrixNM<Real, N,N> *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(3 * eps, static_cast<Real>(1) / static_cast<Real>(3));

            return NDer2PartialAllByAll(f, point, h, error);
        }

        template <int N>
        static MatrixNM<Real, N,N> NDer2PartialAllByAll(const IVectorFunction<N> &f, const VectorN<Real, N> &point, Real h, MatrixNM<Real, N,N> *error = nullptr)
        {
            MatrixNM<Real, N,N> ret;

            for( int i=0; i<N; i++)
                for( int j=0; j<N; j++)
                {
                    if( error )
                        ret(i,j) = NDer2Partial(f, i, j, point, h, &((*error)(i,j)));
                    else
                        ret(i,j) = NDer2Partial(f, i, j, point, h);
                }

            return ret;         
        }

        /////////////////////////             ParametricCurve           /////////////////////////
        template <int N>
        static VectorN<Real, N> NDer2(const MML::IParametricCurve<N> &f, Real t, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();
            Real h = std::pow(3 * eps, static_cast<Real>(1) / static_cast<Real>(3));

            return NDer2(f, t, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer2(const MML::IParametricCurve<N> &f, Real t, Real h, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();

            VectorN<Real, N> yh   = f(t + h);
            VectorN<Real, N> ymh  = f(t - h);
            VectorN<Real, N> diff = yh - ymh;
            
            if (error)
            {
                VectorN<Real, N> yth  = f(t + 2 * h);
                VectorN<Real, N> ymth = f(t - 2 * h);
                *error = eps * ((yh + ymh) / (2 * h)).NormL2() + std::abs(((yth - ymth) / 2 - diff).NormL2()) / (6 * h);
            }
            return diff / (2 * h);
        } 

        template <int N>
        static VectorN<Real, N> NSecDer2(const MML::IParametricCurve<N> &f, Real t, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();
            Real h = std::pow(3 * eps, static_cast<Real>(1) / static_cast<Real>(3));

            return NSecDer2(f, t, h, error);
        }

        template <int N>
        static VectorN<Real, N> NSecDer2(const MML::IParametricCurve<N> &f, Real t, Real h, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();

            VectorN<Real, N> yh   = NDer4(f, t + h, error);
            VectorN<Real, N> ymh  = NDer4(f, t - h, error);
            VectorN<Real, N> diff = yh - ymh;
            
            if (error)
            {
                VectorN<Real, N> yth  = NDer4(f, t + 2 * h, error);
                VectorN<Real, N> ymth = NDer4(f, t - 2 * h, error);
                *error = eps * ((yh + ymh) / (2 * h)).NormL2() + std::abs(((yth - ymth) / 2 - diff).NormL2()) / (6 * h);
            }
            return diff / (2 * h);
        }

        template <int N>
        static VectorN<Real, N> NThirdDer2(const MML::IParametricCurve<N> &f, Real t, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();
            Real h = std::pow(3 * eps, static_cast<Real>(1) / static_cast<Real>(3));

            return NThirdDer2(f, t, h, error);
        }

        template <int N>
        static VectorN<Real, N> NThirdDer2(const MML::IParametricCurve<N> &f, Real t, Real h, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();

            VectorN<Real, N> yh   = NSecDer4(f, t + h, error);
            VectorN<Real, N> ymh  = NSecDer4(f, t - h, error);
            VectorN<Real, N> diff = yh - ymh;
            
            if (error)
            {
                VectorN<Real, N> yth  = NSecDer4(f, t + 2 * h, error);
                VectorN<Real, N> ymth = NSecDer4(f, t - 2 * h, error);
                *error = eps * ((yh + ymh) / (2 * h)).NormL2() + std::abs(((yth - ymth) / 2 - diff).NormL2()) / (6 * h);
            }
            return diff / (2 * h);
        }

        /********************************************************************************************************************/
        /********                               Numerical derivatives of FOURTH order                                 ********/
        /********************************************************************************************************************/
        
        //////////////////////////              RealFunction            //////////////////////////
        static Real NDer4(const MML::IRealFunction &f, Real x, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            // Error bound ~eps^4/5
            Real h = std::pow(11.25*eps, (Real)1 / (Real)5);
            //h = detail::make_xph_representable(x, h);

            return NDer4(f, x, h, error);
        }

        static Real NDer4(const MML::IRealFunction &f, Real x, Real h, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();
            
            Real yh   = f(x + h);
            Real ymh  = f(x - h);
            Real y2h  = f(x + 2 * h);
            Real ym2h = f(x - 2 * h);
            
            Real y2 = ym2h - y2h;
            Real y1 = yh - ymh;
            
            if (error)
            {
                // Mathematica code to extract the remainder:
                // Series[(f[x-2*h]+ 8*f[x+h] - 8*f[x-h] - f[x+2*h])/(12*h), {h, 0, 7}]
                Real y3h  = f(x + 3 * h);
                Real ym3h = f(x - 3 * h);

                // Error from fifth derivative:
                *error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
                // Error from function evaluation:
                *error += eps * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
            }
            return (y2 + 8 * y1) / (12 * h);
        }

        static Real NSecDer4(const MML::IRealFunction &f, Real x, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(11.25*eps, (Real)1 / (Real)5);

            return NSecDer4(f, x, h, error);
        }

        static Real NSecDer4(const MML::IRealFunction &f, Real x, Real h, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();
            
            Real yh   = NDer6(f, x + h, error);
            Real ymh  = NDer6(f, x - h, error);
            Real y2h  = NDer6(f, x + 2 * h, error);
            Real ym2h = NDer6(f, x - 2 * h, error);
            
            Real y2 = ym2h - y2h;
            Real y1 = yh - ymh;
            
            if (error)
            {
                Real y3h  = NDer6(f, x + 3 * h, error);
                Real ym3h = NDer6(f, x - 3 * h, error);

                *error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
                *error += eps * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
            }
            return (y2 + 8 * y1) / (12 * h);
        }

        static Real NThirdDer4(const MML::IRealFunction &f, Real x, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(11.25*eps, (Real)1 / (Real)5);

            return NThirdDer4(f, x, h, error);
        }

        static Real NThirdDer4(const MML::IRealFunction &f, Real x, Real h, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();
            
            Real yh   = NSecDer6(f, x + h, error);
            Real ymh  = NSecDer6(f, x - h, error);
            Real y2h  = NSecDer6(f, x + 2 * h, error);
            Real ym2h = NSecDer6(f, x - 2 * h, error);
            
            Real y2 = ym2h - y2h;
            Real y1 = yh - ymh;
            
            if (error)
            {
                Real y3h  = NSecDer6(f, x + 3 * h, error);
                Real ym3h = NSecDer6(f, x - 3 * h, error);

                *error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
                *error += eps * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
            }
            return (y2 + 8 * y1) / (12 * h);
        }

        //////////////////////////             ScalarFunction           //////////////////////////
        template <int N>
        static Real NDer4Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(11.25*eps, (Real)1 / (Real)5);

            return NDer4Partial(f, deriv_index, point, h, error);
        }

        template <int N>
        static Real NDer4Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, Real h, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            Real     orig_x = point[deriv_index];

            VectorN<Real, N> x{point};
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

                *error += eps * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
            }
            return (y2 + 8 * y1) / (12 * h);
        }

        template <int N>
        static Real NSecDer4Partial(const IScalarFunction<N> &f, int der_ind1, int der_ind2, const VectorN<Real, N> &point, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(11.25*eps, (Real)1 / (Real)5);

            return NSecDer4Partial(f, der_ind1, der_ind2, point, h, error);
        }

        template <int N>
        static Real NSecDer4Partial(const IScalarFunction<N> &f, int der_ind1, int der_ind2, const VectorN<Real, N> &point, Real h, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            Real     orig_x = point[der_ind2];
            auto x_eval_pos = point;
            
            x_eval_pos[der_ind2] = orig_x + h;
            Real yh              = NDer6Partial(f, der_ind1, x_eval_pos, error);

            x_eval_pos[der_ind2] = orig_x - h;
            Real ymh             = NDer6Partial(f, der_ind1, x_eval_pos, error);

            x_eval_pos[der_ind2] = orig_x + 2 * h;
            Real y2h             = NDer6Partial(f, der_ind1, x_eval_pos, error);

            x_eval_pos[der_ind2] = orig_x - 2 * h;
            Real ym2h            = NDer6Partial(f, der_ind1, x_eval_pos, error);

            Real y2 = ym2h - y2h;
            Real y1 = yh - ymh;
            
            if (error)
            {
                x_eval_pos[der_ind2] = orig_x + 3 * h;
                Real y3h             = NDer6Partial(f, der_ind1, x_eval_pos, error);

                x_eval_pos[der_ind2] = orig_x - 3 * h;
                Real ym3h            = NDer6Partial(f, der_ind1, x_eval_pos, error);

                *error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);

                *error += eps * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
            }
            return (y2 + 8 * y1) / (12 * h);
        }

        template <int N>
        static VectorN<Real, N> NDer4PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(11.25*eps, (Real)1 / (Real)5);

            return NDer4PartialByAll(f, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer4PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, Real h, VectorN<Real, N> *error = nullptr)
        {
            VectorN<Real, N> ret;

            for( int i=0; i<N; i++)
            {
                if( error )
                    ret[i] = NDer4Partial(f, i, point, h, &(*error)[i]);
                else
                    ret[i] = NDer4Partial(f, i, point, h);
            }

            return ret;
        }

        //////////////////////////             VectorFunction           //////////////////////////
        template <int N>
        static Real NDer4Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(11.25*eps, (Real)1 / (Real)5);

            return NDer4Partial(f, func_index, deriv_index, point, h, error);
        }

        template <int N>
        static Real NDer4Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, Real h, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            Real     orig_x = point[deriv_index];

            VectorN<Real, N> x{point};
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

                *error += eps * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
            }
            return (y2 + 8 * y1) / (12 * h);
        }

        template <int N>
        static VectorN<Real, N> NDer4PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(11.25*eps, (Real)1 / (Real)5);

            return NDer4PartialByAll(f, func_index, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer4PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, Real h, VectorN<Real, N> *error = nullptr)
        {
            VectorN<Real, N> ret;

            for( int i=0; i<N; i++)
            {
                if( error )
                    ret[i] = NDer4Partial(f, func_index, i, point, h, &(*error)[i]);
                else
                    ret[i] = NDer4Partial(f, func_index, i, point, h);
            }

            return ret;         
        }

        template <int N>
        static MatrixNM<Real, N,N> NDer4PartialAllByAll(const IVectorFunction<N> &f, const VectorN<Real, N> &point, MatrixNM<Real, N,N> *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(3 * eps, static_cast<Real>(1) / static_cast<Real>(3));

            return NDer4PartialAllByAll(f, point, h, error);
        }

        template <int N>
        static MatrixNM<Real, N,N> NDer4PartialAllByAll(const IVectorFunction<N> &f, const VectorN<Real, N> &point, Real h, MatrixNM<Real, N,N> *error = nullptr)
        {
            MatrixNM<Real, N,N> ret;

            for( int i=0; i<N; i++)
                for( int j=0; j<N; j++)
                {
                    if( error )
                        ret(i,j) = NDer4Partial(f, i, j, point, h, &((*error)(i,j)));
                    else
                        ret(i,j) = NDer4Partial(f, i, j, point, h);
                }

            return ret;         
        }

        /////////////////////////             ParametricCurve           /////////////////////////
        template <int N>
        static VectorN<Real, N> NDer4(const MML::IParametricCurve<N> &f, Real t, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();
            Real h = std::pow(11.25*eps, 1.0 / 5.0);

            return NDer4(f, t, h, error);
        }        

        template <int N>
        static VectorN<Real, N> NDer4(const MML::IParametricCurve<N> &f, Real t, Real h, Real* error = nullptr)
        {
            VectorN<Real, N> yh   = f(t + h);
            VectorN<Real, N> ymh  = f(t - h);
            VectorN<Real, N> y2h  = f(t + 2 * h);
            VectorN<Real, N> ym2h = f(t - 2 * h);
            
            VectorN<Real, N> y2 = ym2h - y2h;
            VectorN<Real, N> y1 = yh - ymh;

            if (error)
            {
                VectorN<Real, N> y3h = f(t + 3 * h);
                VectorN<Real, N> ym3h = f(t - 3 * h);
                
                const Real eps = std::numeric_limits<Real>::epsilon();

                *error = std::abs((y3h - ym3h).NormL2() / 2 + 2 * (ym2h - y2h).NormL2() + 5 * (yh - ymh).NormL2() / 2) / (30 * h);
                *error += eps * (y2h.NormL2() + ym2h.NormL2() + 8 * (ymh.NormL2() + yh.NormL2())) / (12 * h);
            }
            return (y2 + 8 * y1) / (12 * h);
        }

        template <int N>
        static VectorN<Real, N> NSecDer4(const MML::IParametricCurve<N> &f, Real t, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();
            Real h = std::pow(11.25*eps, 1.0 / 5.0);

            return NSecDer4(f, t, h, error);
        }        

        template <int N>
        static VectorN<Real, N> NSecDer4(const MML::IParametricCurve<N> &f, Real t, Real h, Real* error = nullptr)
        {
            VectorN<Real, N> yh   = NDer6(f, t + h, error);
            VectorN<Real, N> ymh  = NDer6(f, t - h, error);
            VectorN<Real, N> y2h  = NDer6(f, t + 2 * h, error);
            VectorN<Real, N> ym2h = NDer6(f, t - 2 * h, error);
            
            VectorN<Real, N> y2 = ym2h - y2h;
            VectorN<Real, N> y1 = yh - ymh;

            if (error)
            {
                VectorN<Real, N> y3h  = NDer6(f, t + 3 * h, error);
                VectorN<Real, N> ym3h = NDer6(f, t - 3 * h, error);
                
                const Real eps = std::numeric_limits<Real>::epsilon();

                *error = std::abs((y3h - ym3h).NormL2() / 2 + 2 * (ym2h - y2h).NormL2() + 5 * (yh - ymh).NormL2() / 2) / (30 * h);
                *error += eps * (y2h.NormL2() + ym2h.NormL2() + 8 * (ymh.NormL2() + yh.NormL2())) / (12 * h);
            }
            return (y2 + 8 * y1) / (12 * h);
        }


        template <int N>
        static VectorN<Real, N> NThirdDer4(const MML::IParametricCurve<N> &f, Real t, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();
            Real h = std::pow(11.25*eps, 1.0 / 5.0);

            return NThirdDer4(f, t, h, error);
        }        

        template <int N>
        static VectorN<Real, N> NThirdDer4(const MML::IParametricCurve<N> &f, Real t, Real h, Real* error = nullptr)
        {
            VectorN<Real, N> yh   = NSecDer6(f, t + h, error);
            VectorN<Real, N> ymh  = NSecDer6(f, t - h, error);
            VectorN<Real, N> y2h  = NSecDer6(f, t + 2 * h, error);
            VectorN<Real, N> ym2h = NSecDer6(f, t - 2 * h, error);
            
            VectorN<Real, N> y2 = ym2h - y2h;
            VectorN<Real, N> y1 = yh - ymh;

            if (error)
            {
                VectorN<Real, N> y3h  = NSecDer6(f, t + 3 * h, error);
                VectorN<Real, N> ym3h = NSecDer6(f, t - 3 * h, error);
                
                const Real eps = std::numeric_limits<Real>::epsilon();

                *error = std::abs((y3h - ym3h).NormL2() / 2 + 2 * (ym2h - y2h).NormL2() + 5 * (yh - ymh).NormL2() / 2) / (30 * h);
                *error += eps * (y2h.NormL2() + ym2h.NormL2() + 8 * (ymh.NormL2() + yh.NormL2())) / (12 * h);
            }
            return (y2 + 8 * y1) / (12 * h);
        }        
        
        /********************************************************************************************************************/
        /********                               Numerical derivatives of SIXTH order                                 ********/
        /********************************************************************************************************************/
        
        //////////////////////////              RealFunction            //////////////////////////
        static Real NDer6(const MML::IRealFunction &f, Real x, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            // Error bound ~eps^6/7
            // Error: h^6f^(7)(x)/140 + 5|f(x)|eps/h
            Real h = std::pow(eps / 168, (Real)1 / (Real)7);
            //h = detail::make_xph_representable(x, h);

            return NDer6(f, x, h, error);
        }

        static Real NDer6(const MML::IRealFunction &f, Real x, Real h, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            Real yh  = f(x + h);
            Real ymh = f(x - h);
            Real y1  = yh - ymh;
            Real y2  = f(x - 2 * h) - f(x + 2 * h);
            Real y3  = f(x + 3 * h) - f(x - 3 * h);

            if (error)
            {
                // Mathematica code to generate fd scheme for 7th derivative:
                // Sum[(-1)^i*Binomial[7, i]*(f[x+(3-i)*h] + f[x+(4-i)*h])/2, {i, 0, 7}]
                // Mathematica to demonstrate that this is a finite difference formula for 7th derivative:
                // Series[(f[x+4*h]-f[x-4*h] + 6*(f[x-3*h] - f[x+3*h]) + 14*(f[x-h] - f[x+h] + f[x+2*h] - f[x-2*h]))/2, {h, 0, 15}]
                Real y7 = (f(x + 4 * h) - f(x - 4 * h) - 6 * y3 - 14 * y1 - 14 * y2) / 2;
                *error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh))*eps / h;
            }
            return (y3 + 9 * y2 + 45 * y1) / (60 * h);
        }               

        static Real NSecDer6(const MML::IRealFunction &f, Real x, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(eps / 168, (Real)1 / (Real)7);

            return NSecDer6(f, x, h, error);
        }

        static Real NSecDer6(const MML::IRealFunction &f, Real x, Real h, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            Real yh  = NDer8(f, x + h, error);
            Real ymh = NDer8(f, x - h, error);
            Real y1  = yh - ymh;
            Real y2  = NDer8(f, x - 2 * h, error) - NDer8(f, x + 2 * h, error);
            Real y3  = NDer8(f, x + 3 * h, error) - NDer8(f, x - 3 * h, error);

            if (error)
            {
                Real y7 = (NDer8(f, x + 4 * h, error) - NDer8(f, x - 4 * h, error) - 6 * y3 - 14 * y1 - 14 * y2) / 2;
                *error    = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh))*eps / h;
            }
            return (y3 + 9 * y2 + 45 * y1) / (60 * h);
        }

        static Real NThirdDer6(const MML::IRealFunction &f, Real x, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(eps / 168, (Real)1 / (Real)7);

            return NThirdDer6(f, x, h, error);
        }

        static Real NThirdDer6(const MML::IRealFunction &f, Real x, Real h, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            Real yh  = NSecDer8(f, x + h, error);
            Real ymh = NSecDer8(f, x - h, error);
            Real y1  = yh - ymh;
            Real y2  = NSecDer8(f, x - 2 * h, error) - NSecDer8(f, x + 2 * h, error);
            Real y3  = NSecDer8(f, x + 3 * h, error) - NSecDer8(f, x - 3 * h, error);

            if (error)
            {
                Real y7 = (NSecDer8(f, x + 4 * h, error) - NSecDer8(f, x - 4 * h, error) - 6 * y3 - 14 * y1 - 14 * y2) / 2;
                *error    = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh))*eps / h;
            }
            return (y3 + 9 * y2 + 45 * y1) / (60 * h);
        }

        //////////////////////////             ScalarFunction           //////////////////////////
        template <int N>
        static Real NDer6Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(eps / 168, (Real)1 / (Real)7);

            return NDer6Partial(f, deriv_index, point, h, error);
        }

        template <int N>
        static Real NDer6Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, Real h, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            Real     orig_x = point[deriv_index];

            VectorN<Real, N> x{point};

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
                *error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh))*eps / h;
            }
            return (y3 + 9 * y2 + 45 * y1) / (60 * h);
        }

        template <int N>
        static Real NSecDer6Partial(const IScalarFunction<N> &f, int der_ind1, int der_ind2,const VectorN<Real, N> &point, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(eps / 168, (Real)1 / (Real)7);

            return NSecDer6Partial(f, der_ind1, der_ind2, point, h, error);
        }

        template <int N>
        static Real NSecDer6Partial(const IScalarFunction<N> &f, int der_ind1, int der_ind2, const VectorN<Real, N> &point, Real h, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            Real     orig_x = point[der_ind2];
            auto x_eval_pos = point;

            x_eval_pos[der_ind2] = orig_x + h;
            Real yh              = NDer6Partial(f, der_ind1, x_eval_pos, error);

            x_eval_pos[der_ind2] = orig_x - h;
            Real ymh             = NDer6Partial(f, der_ind1, x_eval_pos, error);

            x_eval_pos[der_ind2] = orig_x + 2 * h;
            Real y2h             = NDer6Partial(f, der_ind1, x_eval_pos, error);

            x_eval_pos[der_ind2] = orig_x - 2 * h;
            Real ym2h            = NDer6Partial(f, der_ind1, x_eval_pos, error);

            x_eval_pos[der_ind2] = orig_x + 3 * h;
            Real y3h             = NDer6Partial(f, der_ind1, x_eval_pos, error);

            x_eval_pos[der_ind2] = orig_x - 3 * h;
            Real ym3h            = NDer6Partial(f, der_ind1, x_eval_pos, error);

            Real y1 = yh - ymh;
            Real y2 = ym2h - y2h;
            Real y3 = y3h - ym3h;

            if (error)
            {
                x_eval_pos[der_ind2] = orig_x + 4 * h;
                Real y4h             = NDer6Partial(f, der_ind1, x_eval_pos, error);

                x_eval_pos[der_ind2] = orig_x - 4 * h;
                Real ym4h            = NDer6Partial(f, der_ind1, x_eval_pos, error);

                Real y7 = (y4h - ym4h - 6 * y3 - 14 * y1 - 14 * y2) / 2;
                *error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh))*eps / h;
            }
            return (y3 + 9 * y2 + 45 * y1) / (60 * h);
        }

        template <int N>
        static VectorN<Real, N> NDer6PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(eps / 168, (Real)1 / (Real)7);

            return NDer6PartialByAll(f, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer6PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, Real h, VectorN<Real, N> *error = nullptr)
        {
            VectorN<Real, N> ret;

            for( int i=0; i<N; i++)
            {
                if( error )
                    ret[i] = NDer6Partial(f, i, point, h, &(*error)[i]);
                else
                    ret[i] = NDer6Partial(f, i, point, h);
            }

            return ret;
        }

        //////////////////////////             VectorFunction           //////////////////////////
        template <int N>
        static Real NDer6Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(eps / 168, (Real)1 / (Real)7);

            return NDer6Partial(f, func_index, deriv_index, point, h, error);
        }

        template <int N>
        static Real NDer6Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, Real h, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            Real     orig_x = point[deriv_index];

            VectorN<Real, N> x{point};

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
                *error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh))*eps / h;
            }
            return (y3 + 9 * y2 + 45 * y1) / (60 * h);
        }

        template <int N>
        static VectorN<Real, N> NDer6PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(eps / 168, (Real)1 / (Real)7);

            return NDer6PartialByAll(f, func_index, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer6PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, Real h, VectorN<Real, N> *error = nullptr)
        {
            VectorN<Real, N> ret;

            for( int i=0; i<N; i++)
            {
                if( error )
                    ret[i] = NDer6Partial(f, func_index, i, point, h, &(*error)[i]);
                else
                    ret[i] = NDer6Partial(f, func_index, i, point, h);
            }

            return ret;         
        }

        template <int N>
        static MatrixNM<Real, N,N> NDer6PartialAllByAll(const IVectorFunction<N> &f, const VectorN<Real, N> &point, MatrixNM<Real, N,N> *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(3 * eps, static_cast<Real>(1) / static_cast<Real>(3));

            return NDer6PartialAllByAll(f, point, h, error);
        }

        template <int N>
        static MatrixNM<Real, N,N> NDer6PartialAllByAll(const IVectorFunction<N> &f, const VectorN<Real, N> &point, Real h, MatrixNM<Real, N,N> *error = nullptr)
        {
            MatrixNM<Real, N,N> ret;

            for( int i=0; i<N; i++)
                for( int j=0; j<N; j++)
                {
                    if( error )
                        ret(i,j) = NDer6Partial(f, i, j, point, h, &((*error)(i,j)));
                    else
                        ret(i,j) = NDer6Partial(f, i, j, point, h);
                }

            return ret;         
        }

        /////////////////////////             ParametricCurve           /////////////////////////
        template <int N>
        static VectorN<Real, N> NDer6(const MML::IParametricCurve<N> &f, Real t, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(3 * eps, static_cast<Real>(1) / static_cast<Real>(3));;

            return NDer6(f, t, h, error);
        }        

        template <int N>
        static VectorN<Real, N> NDer6(const MML::IParametricCurve<N> &f, Real t, Real h, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            VectorN<Real, N> yh  = f(t + h);
            VectorN<Real, N> ymh = f(t - h);
            VectorN<Real, N> y1  = yh - ymh;
            VectorN<Real, N> y2  = f(t - 2 * h) - f(t + 2 * h);
            VectorN<Real, N> y3  = f(t + 3 * h) - f(t - 3 * h);

            if (error)
            {
                VectorN<Real, N> y7 = (f(t + 4 * h) - f(t - 4 * h) - 6 * y3 - 14 * y1 - 14 * y2) / 2;
                *error = y7.NormL2() / (140 * h) + 5 * (yh.NormL2() + ymh.NormL2())*eps / h;
            }
            return (y3 + 9 * y2 + 45 * y1) / (60 * h);
        }

        template <int N>
        static VectorN<Real, N> NSecDer6(const MML::IParametricCurve<N> &f, Real t, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(3 * eps, static_cast<Real>(1) / static_cast<Real>(3));;

            return NSecDer6(f, t, h, error);
        }             

        template <int N>
        static VectorN<Real, N> NSecDer6(const MML::IParametricCurve<N> &f, Real t, Real h, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            VectorN<Real, N> yh  = NDer8(f, t + h, error);
            VectorN<Real, N> ymh = NDer8(f, t - h, error);
            VectorN<Real, N> y1  = yh - ymh;
            VectorN<Real, N> y2  = NDer8(f, t - 2 * h, error) - NDer8(f, t + 2 * h, error);
            VectorN<Real, N> y3  = NDer8(f, t + 3 * h, error) - NDer8(f, t - 3 * h, error);

            if (error)
            {
                VectorN<Real, N> y7 = (NDer8(f, t + 4 * h, error) - NDer8(f, t - 4 * h, error) - 6 * y3 - 14 * y1 - 14 * y2) / 2;
                *error = y7.NormL2() / (140 * h) + 5 * (yh.NormL2() + ymh.NormL2())*eps / h;
            }
            return (y3 + 9 * y2 + 45 * y1) / (60 * h);
        }

        template <int N>
        static VectorN<Real, N> NThirdDer6(const MML::IParametricCurve<N> &f, Real t, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real       h   = std::pow(3 * eps, static_cast<Real>(1) / static_cast<Real>(3));;

            return NThirdDer6(f, t, h, error);
        }             

        template <int N>
        static VectorN<Real, N> NThirdDer6(const MML::IParametricCurve<N> &f, Real t, Real h, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            VectorN<Real, N> yh  = NSecDer8(f, t + h, error);
            VectorN<Real, N> ymh = NSecDer8(f, t - h, error);
            VectorN<Real, N> y1  = yh - ymh;
            VectorN<Real, N> y2  = NSecDer8(f, t - 2 * h, error) - NSecDer8(f, t + 2 * h, error);
            VectorN<Real, N> y3  = NSecDer8(f, t + 3 * h, error) - NSecDer8(f, t - 3 * h, error);

            if (error)
            {
                VectorN<Real, N> y7 = (NSecDer8(f, t + 4 * h, error) - NSecDer8(f, t - 4 * h, error) - 6 * y3 - 14 * y1 - 14 * y2) / 2;
                *error = y7.NormL2() / (140 * h) + 5 * (yh.NormL2() + ymh.NormL2())*eps / h;
            }
            return (y3 + 9 * y2 + 45 * y1) / (60 * h);
        }
        
        /********************************************************************************************************************/
        /********                               Numerical derivatives of EIGHTH order                                ********/
        /********************************************************************************************************************/
        
        //////////////////////////              RealFunction            //////////////////////////
        static Real NDer8(const MML::IRealFunction &f, Real x, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();
            // Error bound ~eps^8/9.
            // In Real precision, we only expect to lose two digits of precision while using this formula, at the cost of 8 function evaluations.
            // Error: h^8|f^(9)(x)|/630 + 7|f(x)|eps/h assuming 7 unstabilized additions.
            // Mathematica code to get the error:
            // Series[(f[x+h]-f[x-h])*(4/5) + (1/5)*(f[x-2*h] - f[x+2*h]) + (4/105)*(f[x+3*h] - f[x-3*h]) + (1/280)*(f[x-4*h] - f[x+4*h]), {h, 0, 9}]
            // If we used Kahan summation, we could get the max error down to h^8|f^(9)(x)|/630 + |f(x)|eps/h.
            Real h = std::pow(551.25*eps, (Real)1 / (Real)9);
            //h = detail::make_xph_representable(x, h);

            return NDer8(f, x, h, error);
        }

        static Real NDer8(const MML::IRealFunction &f, Real x, Real h, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            Real yh  = f(x + h);
            Real ymh = f(x - h);
            Real y1  = yh - ymh;
            Real y2  = f(x - 2 * h) - f(x + 2 * h);
            Real y3  = f(x + 3 * h) - f(x - 3 * h);
            Real y4  = f(x - 4 * h) - f(x + 4 * h);

            Real tmp1 = 3 * y4 / 8 + 4 * y3;
            Real tmp2 = 21 * y2 + 84 * y1;

            if (error)
            {
                // Mathematica code to generate fd scheme for 7th derivative:
                // Sum[(-1)^i*Binomial[9, i]*(f[x+(4-i)*h] + f[x+(5-i)*h])/2, {i, 0, 9}]
                // Mathematica to demonstrate that this is a finite difference formula for 7th derivative:
                // Series[(f[x+5*h]-f[x- 5*h])/2 + 4*(f[x-4*h] - f[x+4*h]) + 27*(f[x+3*h] - f[x-3*h])/2 + 24*(f[x-2*h]  - f[x+2*h]) + 21*(f[x+h] - f[x-h]), {h, 0, 15}]
                Real f9 = (f(x + 5 * h) - f(x - 5 * h)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;
                *error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh))*eps / h;
            }
            return (tmp1 + tmp2) / (105 * h);
        }

        static Real NSecDer8(const MML::IRealFunction &f, Real x, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();
            Real h = std::pow(551.25*eps, (Real)1 / (Real)9);

            return NSecDer8(f, x, h, error);
        }

        static Real NSecDer8(const MML::IRealFunction &f, Real x, Real h, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();

            Real yh  = NDer8(f, x + h, error);
            Real ymh = NDer8(f, x - h, error);
            Real y1  = yh - ymh;
            Real y2  = NDer8(f, x - 2 * h, error) - NDer8(f, x + 2 * h, error);
            Real y3  = NDer8(f, x + 3 * h, error) - NDer8(f, x - 3 * h, error);
            Real y4  = NDer8(f, x - 4 * h, error) - NDer8(f, x + 4 * h, error);

            Real tmp1 = 3 * y4 / 8 + 4 * y3;
            Real tmp2 = 21 * y2 + 84 * y1;

            if (error)
            {
                Real f9 = (NDer8(f, x + 5 * h, error) - NDer8(f, x - 5 * h, error)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;
                *error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh))*eps / h;
            }
            return (tmp1 + tmp2) / (105 * h);
        }

        static Real NThirdDer8(const MML::IRealFunction &f, Real x, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();
            Real h = std::pow(551.25*eps, (Real)1 / (Real)9);

            return NThirdDer8(f, x, h, error);
        }

        static Real NThirdDer8(const MML::IRealFunction &f, Real x, Real h, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();

            Real yh  = NSecDer8(f, x + h, error);
            Real ymh = NSecDer8(f, x - h, error);
            Real y1  = yh - ymh;
            Real y2  = NSecDer8(f, x - 2 * h, error) - NSecDer8(f, x + 2 * h, error);
            Real y3  = NSecDer8(f, x + 3 * h, error) - NSecDer8(f, x - 3 * h, error);
            Real y4  = NSecDer8(f, x - 4 * h, error) - NSecDer8(f, x + 4 * h, error);

            Real tmp1 = 3 * y4 / 8 + 4 * y3;
            Real tmp2 = 21 * y2 + 84 * y1;

            if (error)
            {
                Real f9 = (NSecDer8(f, x + 5 * h, error) - NSecDer8(f, x - 5 * h, error)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;
                *error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh))*eps / h;
            }
            return (tmp1 + tmp2) / (105 * h);
        }        

        //////////////////////////             ScalarFunction           //////////////////////////
        template <int N>
        static Real NDer8Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(551.25*eps, 1.0 / 9.0);

            return NDer8Partial(f, deriv_index, point, h, error);
        }

        template <int N>
        static Real NDer8Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, Real h, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            Real     orig_x = point[deriv_index];

            VectorN<Real, N> x{point};

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
                Real y5h     = f(x);

                x[deriv_index] = orig_x - 5 * h;
                Real ym5h    = f(x);

                Real f9 = (y5h - ym5h) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;
                *error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh))*eps / h;

            }

            return (tmp1 + tmp2) / (105 * h);            
        }

        template <int N>
        static Real NSecDer8Partial(const IScalarFunction<N> &f, int der_ind1, int der_ind2,  const VectorN<Real, N> &point, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(551.25*eps, 1.0 / 9.0);

            return NSecDer8Partial(f, der_ind1, der_ind2, point, h, error);
        }

        template <int N>
        static Real NSecDer8Partial(const IScalarFunction<N> &f, int der_ind1, int der_ind2,  const VectorN<Real, N> &point, Real h, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            Real     orig_x = point[der_ind2];
            auto x_eval_pos = point;

            x_eval_pos[der_ind2] = orig_x + h;
            Real yh              = NDer8Partial(f, der_ind1, x_eval_pos, error);

            x_eval_pos[der_ind2] = orig_x - h;
            Real ymh             = NDer8Partial(f, der_ind1, x_eval_pos, error);

            x_eval_pos[der_ind2] = orig_x + 2 * h;
            Real y2h             = NDer8Partial(f, der_ind1, x_eval_pos, error);

            x_eval_pos[der_ind2] = orig_x - 2 * h;
            Real ym2h            = NDer8Partial(f, der_ind1, x_eval_pos, error);

            x_eval_pos[der_ind2] = orig_x + 3 * h;
            Real y3h             = NDer8Partial(f, der_ind1, x_eval_pos, error);

            x_eval_pos[der_ind2] = orig_x - 3 * h;
            Real ym3h            = NDer8Partial(f, der_ind1, x_eval_pos, error);

            x_eval_pos[der_ind2] = orig_x + 4 * h;
            Real y4h             = NDer8Partial(f, der_ind1, x_eval_pos, error);

            x_eval_pos[der_ind2] = orig_x - 4 * h;
            Real ym4h            = NDer8Partial(f, der_ind1, x_eval_pos, error);

            Real y1 = yh - ymh;
            Real y2 = ym2h - y2h;
            Real y3 = y3h - ym3h;
            Real y4 = ym4h - y4h;

            Real tmp1 = 3 * y4 / 8 + 4 * y3;
            Real tmp2 = 21 * y2 + 84 * y1;

            if (error)
            {
                x_eval_pos[der_ind2] = orig_x + 5 * h;
                Real y5h             = NDer8Partial(f, der_ind1, x_eval_pos, error);

                x_eval_pos[der_ind2] = orig_x - 5 * h;
                Real ym5h            = NDer8Partial(f, der_ind1, x_eval_pos, error);

                Real f9 = (y5h - ym5h) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;
                *error  = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh))*eps / h;
            }

            return (tmp1 + tmp2) / (105 * h);            
        }

        template <int N>
        static VectorN<Real, N> NDer8PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(551.25*eps, (Real)1 / (Real)9);

            return NDer8PartialByAll(f, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer8PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, Real h, VectorN<Real, N> *error = nullptr)
        {
            VectorN<Real, N> ret;

            for( int i=0; i<N; i++)
            {
                if( error )
                    ret[i] = NDer8Partial(f, i, point, h, &(*error)[i]);
                else
                    ret[i] = NDer8Partial(f, i, point, h);
            }

            return ret;
        }

        //////////////////////////             VectorFunction           //////////////////////////
        template <int N>
        static Real NDer8Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(551.25*eps, (Real)1 / (Real)9);

            return NDer8Partial(f, func_index, deriv_index, point, h, error);
        }

        template <int N>
        static Real NDer8Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, Real h, Real *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            Real     orig_x = point[deriv_index];

            VectorN<Real, N> x{point};

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
                *error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh))*eps / h;
            }

            return (tmp1 + tmp2) / (105 * h);            
        }        

        template <int N>
        static VectorN<Real, N> NDer8PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(551.25*eps, (Real)1 / (Real)9);

            return NDer8PartialByAll(f, func_index, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer8PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, Real h, VectorN<Real, N> *error = nullptr)
        {
            VectorN<Real, N> ret;

            for( int i=0; i<N; i++)
            {
                if( error )
                    ret[i] = NDer8Partial(f, func_index, i, point, h, &(*error)[i]);
                else
                    ret[i] = NDer8Partial(f, func_index, i, point, h);
            }

            return ret;         
        }

        template <int N>
        static MatrixNM<Real, N,N> NDer8PartialAllByAll(const IVectorFunction<N> &f, const VectorN<Real, N> &point, MatrixNM<Real, N,N> *error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h = std::pow(551.25*eps, (Real)1 / (Real)9);

            return NDer8PartialAllByAll(f, point, h, error);
        }

        template <int N>
        static MatrixNM<Real, N,N> NDer8PartialAllByAll(const IVectorFunction<N> &f, const VectorN<Real, N> &point, Real h, MatrixNM<Real, N,N> *error = nullptr)
        {
            MatrixNM<Real, N,N> ret;

            for( int i=0; i<N; i++)
                for( int j=0; j<N; j++)
                {
                    if( error )
                        ret(i,j) = NDer8Partial(f, i, j, point, h, &((*error)(i,j)));
                    else
                        ret(i,j) = NDer8Partial(f, i, j, point, h);
                }

            return ret;         
        }

        /////////////////////////             ParametricCurve           /////////////////////////
        template <int N>
        static VectorN<Real, N> NDer8(const MML::IParametricCurve<N> &f, Real t, Real* error = nullptr)
        {
            const Real eps = std::numeric_limits<Real>::epsilon();
            Real h         = std::pow(551.25*eps, (Real)1 / (Real)9);

            return NDer8(f, t, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer8(const MML::IParametricCurve<N> &f, Real t, Real h, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            VectorN<Real, N> yh  = f(t + h);
            VectorN<Real, N> ymh = f(t - h);
            VectorN<Real, N> y1  = yh - ymh;
            VectorN<Real, N> y2  = f(t - 2 * h) - f(t + 2 * h);
            VectorN<Real, N> y3  = f(t + 3 * h) - f(t - 3 * h);
            VectorN<Real, N> y4  = f(t - 4 * h) - f(t + 4 * h);

            VectorN<Real, N> tmp1 = 3 * y4 / 8 + 4 * y3;
            VectorN<Real, N> tmp2 = 21 * y2 + 84 * y1;

            if (error)
            {
                VectorN<Real, N> f9 = (f(t + 5 * h) - f(t - 5 * h)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;
                *error = f9.NormL2() / (630 * h) + 7 * (yh.NormL2() + ymh.NormL2())*eps / h;
            }
            return (tmp1 + tmp2) / (105 * h);
        }

        template <int N>
        static VectorN<Real, N> NSecDer8(const MML::IParametricCurve<N> &f, Real t, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h         = std::pow(551.25*eps, (Real)1 / (Real)9);

            return NSecDer8(f, t, h, error);
        }

        template <int N>
        static VectorN<Real, N> NSecDer8(const MML::IParametricCurve<N> &f, Real t, Real h, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            VectorN<Real, N> yh  = NDer8(f, t + h, error);
            VectorN<Real, N> ymh = NDer8(f, t - h, error);
            VectorN<Real, N> y1  = yh - ymh;
            VectorN<Real, N> y2  = NDer8(f, t - 2 * h, error) - NDer8(f, t + 2 * h, error);
            VectorN<Real, N> y3  = NDer8(f, t + 3 * h, error) - NDer8(f, t - 3 * h, error);
            VectorN<Real, N> y4  = NDer8(f, t - 4 * h, error) - NDer8(f, t + 4 * h, error);

            VectorN<Real, N> tmp1 = 3 * y4 / 8 + 4 * y3;
            VectorN<Real, N> tmp2 = 21 * y2 + 84 * y1;

            if (error)
            {
                VectorN<Real, N> f9 = (NDer8(f, t + 5 * h, error) - NDer8(f, t - 5 * h, error)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;
                *error = f9.NormL2() / (630 * h) + 7 * (yh.NormL2() + ymh.NormL2())*eps / h;
            }
            return (tmp1 + tmp2) / (105 * h);
        }        

        template <int N>
        static VectorN<Real, N> NThirdDer8(const MML::IParametricCurve<N> &f, Real t, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();
            Real h         = std::pow(551.25*eps, (Real)1 / (Real)9);

            return NThirdDer8(f, t, h, error);
        }

        template <int N>
        static VectorN<Real, N> NThirdDer8(const MML::IParametricCurve<N> &f, Real t, Real h, Real* error = nullptr)
        {
            const Real eps = (std::numeric_limits<Real>::epsilon)();

            VectorN<Real, N> yh  = NSecDer8(f, t + h, error);
            VectorN<Real, N> ymh = NSecDer8(f, t - h, error);
            VectorN<Real, N> y1  = yh - ymh;
            VectorN<Real, N> y2  = NSecDer8(f, t - 2 * h, error) - NSecDer8(f, t + 2 * h, error);
            VectorN<Real, N> y3  = NSecDer8(f, t + 3 * h, error) - NSecDer8(f, t - 3 * h, error);
            VectorN<Real, N> y4  = NSecDer8(f, t - 4 * h, error) - NSecDer8(f, t + 4 * h, error);

            VectorN<Real, N> tmp1 = 3 * y4 / 8 + 4 * y3;
            VectorN<Real, N> tmp2 = 21 * y2 + 84 * y1;

            if (error)
            {
                VectorN<Real, N> f9 = (NSecDer8(f, t + 5 * h, error) - NSecDer8(f, t - 5 * h, error)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;
                *error = f9.NormL2() / (630 * h) + 7 * (yh.NormL2() + ymh.NormL2())*eps / h;
            }
            return (tmp1 + tmp2) / (105 * h);
        }     
        
        /********************************************************************************************************************/
        /********                            Definitions of default derivation functions                             ********/
        /********************************************************************************************************************/
        static inline Real(*Derive)(const MML::IRealFunction &f, Real x, Real* error) = Derivation::NDer4;

        static inline Real(*DeriveSec)(const MML::IRealFunction &f, Real x, Real* error) = Derivation::NSecDer4;

        static inline Real(*DeriveThird)(const MML::IRealFunction &f, Real x, Real* error) = Derivation::NThirdDer2;
        
        template<int N>
        static inline Real(*DerivePartial)(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, Real *error) = Derivation::NDer4Partial;

        template<int N>
        static inline Real(*DeriveSecPartial)(const IScalarFunction<N> &f, int der_ind1, int der_ind2, const VectorN<Real, N> &point, Real *error) = Derivation::NSecDer4Partial;

        template<int N>
        static inline VectorN<Real, N>(*DerivePartialAll)(const IScalarFunction<N> &f, const VectorN<Real, N> &point, VectorN<Real, N> *error) = Derivation::NDer4PartialByAll;

        template<int N>
        static inline Real(*DeriveVecPartial)(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, Real *error) = Derivation::NDer4Partial;

        template<int N>
        static inline VectorN<Real, N>(*DeriveVecPartialAll)(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, VectorN<Real, N> *error) = Derivation::NDer4PartialByAll;

        template<int N>
        static inline MatrixNM<Real, N, N>(*DeriveVecPartialAllByAll)(const IVectorFunction<N> &f, const VectorN<Real, N> &point, MatrixNM<Real, N,N> *error) = Derivation::NDer4PartialAllByAll;

        template<int N>
        static inline VectorN<Real, N>(*DeriveCurve)(const IParametricCurve<N> &f, Real x, Real* error) = Derivation::NDer4;

        template<int N>
        static inline VectorN<Real, N>(*DeriveCurveSec)(const IParametricCurve<N> &f, Real x, Real* error) = Derivation::NSecDer4;

        template<int N>
        static inline VectorN<Real, N>(*DeriveCurveThird)(const IParametricCurve<N> &f, Real x, Real* error) = Derivation::NThirdDer4;
    };
}

///////////////////////////   ./include/algorithms/Integration.h   ///////////////////////////



namespace MML
{
    static bool polint(Vector<Real> &xa, Vector<Real> &ya, const double x, double &y, double &dy)
    {
        int i,m,ns=0;
        double den,dif,dift,ho,hp,w;

        int n=(int) xa.size();
        Vector<Real> c(n),d(n);
        dif=fabs(x-xa[0]);
        for (i=0;i<n;i++) {
            if ((dift=fabs(x-xa[i])) < dif) {
                ns=i;
                dif=dift;
            }
            c[i]=ya[i];
            d[i]=ya[i];
        }
        y=ya[ns--];
        for (m=1;m<n;m++) {
            for (i=0;i<n-m;i++) {
                ho=xa[i]-x;
                hp=xa[i+m]-x;
                w=c[i+1]-d[i];
                if ((den=ho-hp) == 0.0) 
                    // nrerror("Error in routine polint");
                    return false;
                den=w/den;
                d[i]=hp*den;
                c[i]=ho*den;
            }
            y += (dy=(2*(ns+1) < (n-m) ? c[ns+1] : d[ns--]));
        }
        return true;
    }

    static void polin2(Vector<Real> &x1a, Vector<Real> &x2a, Matrix<Real> &ya, const double x1,
                const double x2, double &y, double &dy)
    // Given arrays x1a[1..m] and x2a[1..n] of independent variables, and a submatrix of function
    // values ya[1..m][1..n], tabulated at the grid points defined by x1a and x2a; and given values
    // x1 and x2 of the independent variables; this routine returns an interpolated function value y,
    // and an accuracy indication dy (based only on the interpolation in the x1 direction, however).          
    {
        int j,k;

        int m = (int) x1a.size();
        int n = (int) x2a.size();
        Vector<Real> ymtmp(m),ya_t(n);
        for (j=0;j<m;j++) {
            for (k=0;k<n;k++) 
                ya_t[k]=ya[j][k];

            polint(x2a,ya_t,x2,ymtmp[j],dy);
        }
        polint(x1a,ymtmp,x1,y,dy);
    }
	class Integration
	{
		public:

        static double TrapRefine(const IRealFunction &func, const double a, const double b, const int n)
        {
            // This routine computes the nth stage of refinement of an extended trapezoidal rule. func is input
            // as a pointer to the function to be integrated between limits a and b, also input. When called with
            // n=1, the routine returns the crudest estimate of Rab f(x)dx. Subsequent calls with n=2,3,...
            // (in that sequential order) will improve the accuracy by adding 2n-2 additional interior points.
            double x,tnm,sum,del;
            static double s;
            int it,j;

            if (n == 1) {
                return (s=0.5*(b-a)*(func(a)+func(b)));
            } else 
            {
                for (it=1,j=1;j<n-1;j++) 
                    it <<= 1;
            
                tnm=it;
                del=(b-a)/tnm;
                x=a+0.5*del;
            
                for (sum=0.0,j=0;j<it;j++,x+=del) 
                    sum += func(x);
            
                s=0.5*(s+(b-a)*sum/tnm);
            
                return s;
            }
        }

        static double IntegrateTrap(const IRealFunction &func, const double a, const double b, double req_eps)
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
            double s,olds=0.0;
            double diff = 0.0, threshold = 0.0;

            for (j=0;j<Defaults::IntegrateTrapMaxSteps;j++) 
            {
                s=TrapRefine(func,a,b,j+1);

                if (j > 5)
                {
                    diff = s-olds;
                    threshold = req_eps * std::abs(olds);
                    //std::cout << "\ns : " << s << " olds : " << olds <<  " diff : " << diff << " threshold : " << threshold << std::endl;
                    if (std::abs(diff) < threshold || (s == 0.0 && olds == 0.0)) 
                        return s;
                }
                olds=s;
            }
            throw IntegrationTooManySteps("qtrap");
        }
        static double IntegrateTrap(const IRealFunction &func, const double a, const double b)
        {
            return IntegrateTrap(func, a, b, Defaults::IntegrateTrapEPS);
        }

        static double IntegrateSimpson(const IRealFunction &func, const double a, const double b, double req_eps)
        {
            // Returns the integral of the function func from a to b. The parameters EPS can be set to the
            // desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum allowed
            // number of steps. Integration is performed by Simpson’s rule.

            // The routine qsimp will in general be more efficient than qtrap (i.e., require
            // fewer function evaluations) when the function to be integrated has a finite 4th
            // derivative (i.e., a continuous 3rd derivative). The combination of qsimp and its
            // necessary workhorse trapzd is a good one for light-duty work.
            int j;
            double s,st,ost=0.0,os=0.0;

            for (j=0;j<Defaults::IntegrateSimpMaxSteps;j++) 
            {
                st=TrapRefine(func,a,b,j+1);
                s=(4.0*st-ost)/3.0;
                
                if (j > 5)
                    if (std::abs(s-os) < req_eps * std::abs(os) ||
                        (s == 0.0 && os == 0.0)) 
                        return s;
                os=s;
                ost=st;
            }
            throw IntegrationTooManySteps("qsimp");
        }

        static double IntegrateSimpson(const IRealFunction &func, const double a, const double b)
        {
            return IntegrateSimpson(func, a, b, Defaults::IntegrateSimpEPS);
        }

        static double IntegrateRomberg(const IRealFunction &func, double a, double b, double req_eps)
        {
            // Returns the integral of the function func from a to b. Integration is performed by Romberg’s
            // method of order 2K, where, e.g., K=2 is Simpson’s rule.

            // The routine qromb, along with its required trapzd and polint, is quite
            // powerful for sufficiently smooth (e.g., analytic) integrands, integrated over intervals
            // which contain no singularities, and where the endoubleoints are also nonsingular. qromb,
            // in such circumstances, takes many, many fewer function evaluations than either of
            // the routines in x4.2
            const int JMAXP=Defaults::IntegrateRombMaxSteps+1, K=5;
            double ss,dss;
            Vector<double> s(Defaults::IntegrateRombMaxSteps),h(JMAXP),s_t(K),h_t(K);

            h[0]=1.0;
            for (int j=1;j<=Defaults::IntegrateRombMaxSteps;j++) 
            {
                s[j-1]=TrapRefine(func,a,b,j);
            
                if (j >= K) 
                {
                    for (int i=0;i<K;i++) {
                        h_t[i]=h[j-K+i];
                        s_t[i]=s[j-K+i];
                    }
                
                    polint(h_t,s_t,0.0,ss,dss);
                
                    if (std::abs(dss) <= req_eps * std::abs(ss)) 
                        return ss;
                }

                h[j]=0.25*h[j-1];
            }
            throw IntegrationTooManySteps("qromb");
        }

        static double IntegrateRomberg(const IRealFunction &func, const double a, const double b)
        {
            return IntegrateRomberg(func, a, b, Defaults::IntegrateRombEPS);
        }

        static inline double(*Integrate)(const MML::IRealFunction &f, double a, double b, double req_eps) = Integration::IntegrateSimpson;

	};
} // end namespace
///////////////////////////   ./include/basic_types/CoordTransf.h   ///////////////////////////





namespace MML
{
    // ovdje dodati translational, rotational, galilean, lorentzian transf
    // SVE su to transformacije koordinata
    // TODO - dodati Cart2DToPolar kao primjer 2D transformacije
    // TODO - dodati CartesianRotation kao primjer 3D transformacije
    template<typename VectorFrom, typename VectorTo, int N>
    class CoordTransf : ICoordTransf<VectorFrom, VectorTo, N>
    {
        public:
        virtual VectorTo getUnitVector(int ind, const VectorFrom &pos)
        {
            VectorTo ret;

            return ret;
        }

        // transf contravariant vector
        VectorTo contravariantTransf(const VectorFrom &vec, const VectorFrom &pos) 
        {
            VectorFrom ret;

            for( int j=0; j<N; j++ )
            {
                ret[j] = 0;
                for( int k=0; k<N; k++)
                {
                    ret[j] += Derivation::NDer1Partial(this->coordTransfFunc(j), k, pos, 1e-8) * vec[k];
                }
            }

            return ret;
        }
        VectorFrom contravariantTransfInverse(const VectorTo &vec, const VectorTo &pos) 
        {
            VectorFrom ret;

            for( int k=0; k<N; k++ )
            {
                ret[k] = 0;
                for( int j=0; j<N; j++)
                {
                    ret[k] += Derivation::NDer1Partial(this->inverseCoordTransfFunc(k), j, pos, 1e-8) * vec[j];
                }
            }

            return ret;
        }

        // transform covariant vector
        VectorTo covariantTransf(const VectorFrom &vec, const VectorTo &pos)
        {
            VectorTo ret;

            for( int j=0; j<N; j++ )
            {
                ret[j] = 0;
                for( int k=0; k<N; k++)
                {
                    ret[j] += Derivation::NDer1Partial(this->inverseCoordTransfFunc(k), j, pos, 1e-8) * vec[k];
                }
            }

            return ret;
        }

        VectorFrom covariantTransfInverse(const VectorTo &vec, const VectorFrom &pos)
        {
            VectorFrom ret;

            for( int k=0; k<N; k++ )
            {
                ret[k] = 0;
                for( int j=0; j<N; j++)
                {
                    ret[k] += Derivation::NDer1Partial(this->coordTransfFunc(j), k, pos, 1e-8) * vec[j];
                }
            }

            return ret;
        }        

        // transform tensor
    }; 

    class CoordTransfRectilinear : public CoordTransf<VectorN<Real, 3>, VectorN<Real, 3>, 3>
    {
    private:
        Vector3Cartesian _base[3];
        Vector3Cartesian _dual[3];

        ScalarFunctionFromStdFunc<3> _f1;
        ScalarFunctionFromStdFunc<3> _f2;
        ScalarFunctionFromStdFunc<3> _f3;

        ScalarFunctionFromStdFunc<3> _fInverse1;
        ScalarFunctionFromStdFunc<3> _fInverse2;
        ScalarFunctionFromStdFunc<3> _fInverse3;

    public:
        MatrixNM<Real,3,3> _alpha;
        MatrixNM<Real,3,3> _transf;

        CoordTransfRectilinear( VectorN<Real, 3> b1, 
                                VectorN<Real, 3> b2, 
                                VectorN<Real, 3> b3) : _f1( std::function<double(const VectorN<Real, 3>&)> { std::bind( &CoordTransfRectilinear::func1, this, std::placeholders::_1 ) } ),
                                                       _f2( std::function<double(const VectorN<Real, 3>&)> { std::bind( &CoordTransfRectilinear::func2, this, std::placeholders::_1 ) } ),
                                                       _f3( std::function<double(const VectorN<Real, 3>&)> { std::bind( &CoordTransfRectilinear::func3, this, std::placeholders::_1 ) } ),
                                                       _fInverse1( std::function<double(const VectorN<Real, 3>&)> { std::bind( &CoordTransfRectilinear::funcInverse1, this, std::placeholders::_1 ) } ),
                                                       _fInverse2( std::function<double(const VectorN<Real, 3>&)> { std::bind( &CoordTransfRectilinear::funcInverse2, this, std::placeholders::_1 ) } ),
                                                       _fInverse3( std::function<double(const VectorN<Real, 3>&)> { std::bind( &CoordTransfRectilinear::funcInverse3, this, std::placeholders::_1 ) } )
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

            for( int i=0; i<3; i++ )
            {
                for(int j=0; j<3; j++ )
                {
                    _alpha(i,j) =  _base[i][j];
                    _transf(i,j) = _base[j][i];     // transponirano
                }
            }
        }

        double func1(const VectorN<Real, 3> &q) { return ScalarProd(q, MML::Vector3Cartesian(_dual[0])); }
        double func2(const VectorN<Real, 3> &q) { return ScalarProd(q, MML::Vector3Cartesian(_dual[1])); }
        double func3(const VectorN<Real, 3> &q) { return ScalarProd(q, MML::Vector3Cartesian(_dual[2])); }

        double funcInverse1(const VectorN<Real, 3> &q) { return (_transf * q)[0]; }
        double funcInverse2(const VectorN<Real, 3> &q) { return (_transf * q)[1]; }
        double funcInverse3(const VectorN<Real, 3> &q) { return (_transf * q)[2]; }

        Vector3Cartesian    Base(int i) { return _dual[i]; }
        Vector3Cartesian    Dual(int i) { return _dual[i]; }

        VectorN<Real, 3>    transf(const VectorN<Real, 3> &q)           { return VectorN<Real, 3>{ func1(q), func2(q), func3(q) }; }
        VectorN<Real, 3>    transfInverse(const VectorN<Real, 3> &q)    { return VectorN<Real, 3>{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }
        
        IScalarFunction<3>& coordTransfFunc(int i)         
        { 
            if( i == 0 ) return _f1;
            else if( i == 1 ) return _f2;
            else return _f3; 
        }
        IScalarFunction<3>&  inverseCoordTransfFunc(int i)
        {
            if( i == 0 ) return _fInverse1;
            else if( i == 1 ) return _fInverse2;
            else return _fInverse3; 
        }

        bool IsRightHanded()
        {
            Vector3Cartesian cross = VectorProd(_base[0], _base[1]);
            if( ScalarProd(cross, _base[2]) > 0.0 )
                return true;
            else
                return false;
        }
    };

    class CoordTransfSphericalToCartesian : public CoordTransf<Vector3Spherical, Vector3Cartesian, 3>
    {
        // q[0] = r     - radial distance
        // q[1] = theta - inclination
        // q[2] = phi   - azimuthal angle
        public:
        static double func1(const VectorN<Real, 3> &q) { return q[0] * sin(q[1]) * cos(q[2]); }
        static double func2(const VectorN<Real, 3> &q) { return q[0] * sin(q[1]) * sin(q[2]); }
        static double func3(const VectorN<Real, 3> &q) { return q[0] * cos(q[1]); }

        // q[0] = x
        // q[1] = y
        // q[2] = z
        static double funcInverse1(const VectorN<Real, 3> &q) { return sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]); }
        static double funcInverse2(const VectorN<Real, 3> &q) { return atan2(sqrt(q[0]*q[0] + q[1]*q[1]), q[2]); }
        static double funcInverse3(const VectorN<Real, 3> &q) { return atan2(q[1], q[0]); }

        inline static ScalarFunction<3> _func[3] = { ScalarFunction<3>{func1},
                                                                ScalarFunction<3>{func2},
                                                                ScalarFunction<3>{func3}
                                                              };

        inline static ScalarFunction<3> _funcInverse[3] = { ScalarFunction<3>{funcInverse1},
                                                                       ScalarFunction<3>{funcInverse2},
                                                                       ScalarFunction<3>{funcInverse3}
                                                                     };

        Vector3Cartesian     transf(const Vector3Spherical &q)           { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
        Vector3Spherical     transfInverse(const Vector3Cartesian &q)    { return Vector3Spherical{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }
        
        IScalarFunction<3>&  coordTransfFunc(int i)         { return _func[i]; }
        IScalarFunction<3>&  inverseCoordTransfFunc(int i)  { return _funcInverse[i]; }
    };

    class CoordTransfCartesianToSpherical : public CoordTransf<Vector3Cartesian, Vector3Spherical, 3>
    {
        // q[0] = x
        // q[1] = y
        // q[2] = z
        public:
        static double func1(const VectorN<Real, 3> &q) { return sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]); }
        static double func2(const VectorN<Real, 3> &q) { return atan2(sqrt(q[0]*q[0] + q[1]*q[1]), q[2]); }
        static double func3(const VectorN<Real, 3> &q) { return atan2(q[1], q[0]); }

        // q[0] = r     - radial distance
        // q[1] = theta - inclination
        // q[2] = phi   - azimuthal angle
        public:
        static double funcInverse1(const VectorN<Real, 3> &q) { return q[0] * sin(q[1]) * cos(q[2]); }
        static double funcInverse2(const VectorN<Real, 3> &q) { return q[0] * sin(q[1]) * sin(q[2]); }
        static double funcInverse3(const VectorN<Real, 3> &q) { return q[0] * cos(q[1]); }

        inline static ScalarFunction<3> _func[3] = { 
                                                    ScalarFunction<3>{func1},
                                                    ScalarFunction<3>{func2},
                                                    ScalarFunction<3>{func3}
                                                };

        inline static ScalarFunctionFromStdFunc<3> _funcInverse[3] = { 
                                                                        ScalarFunctionFromStdFunc<3>{std::function<double(const VectorN<Real, 3>&)>{funcInverse1}},
                                                                        ScalarFunctionFromStdFunc<3>{std::function<double(const VectorN<Real, 3>&)>{funcInverse2}},
                                                                        ScalarFunctionFromStdFunc<3>{std::function<double(const VectorN<Real, 3>&)>{funcInverse3}}
                                                                    };

        Vector3Spherical     transf(const Vector3Cartesian &q)           { return Vector3Spherical{ func1(q), func2(q), func3(q) }; }
        Vector3Cartesian     transfInverse(const Vector3Spherical &q)    { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

        IScalarFunction<3>&  coordTransfFunc(int i)         { return _func[i]; }
        IScalarFunction<3>&  inverseCoordTransfFunc(int i)  { return _funcInverse[i]; }
    };

    class CoordTransfCylindricalToCartesian : public CoordTransf<Vector3Cylindrical, Vector3Cartesian, 3>
    {
        // q1 = r   - distance from symmetry axis
        // q2 = phi - angle to symmetry axis
        // q3 = z   - z
        public:
        static double func1(const VectorN<Real, 3> &q) { return q[0] * cos(q[1]); }
        static double func2(const VectorN<Real, 3> &q) { return q[0] * sin(q[1]); }
        static double func3(const VectorN<Real, 3> &q) { return q[2]; }

        // q[0] = x
        // q[1] = y
        // q[2] = z
        static double funcInverse1(const VectorN<Real, 3> &q) { return sqrt(q[0]*q[0] + q[1]*q[1]); }
        static double funcInverse2(const VectorN<Real, 3> &q) { return atan2(q[1], q[0]); }
        static double funcInverse3(const VectorN<Real, 3> &q) { return q[2]; }

        inline static ScalarFunction<3> _func[3] = { ScalarFunction<3>{func1},
                                                    ScalarFunction<3>{func2},
                                                    ScalarFunction<3>{func3}
                                                    };

        inline static ScalarFunction<3> _funcInverse[3] = { ScalarFunction<3>{funcInverse1},
                                                            ScalarFunction<3>{funcInverse2},
                                                            ScalarFunction<3>{funcInverse3}
                                                            };

        Vector3Cartesian     transf(const Vector3Cylindrical &q)         { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
        Vector3Cylindrical   transfInverse(const Vector3Cartesian &q)    { return Vector3Cylindrical{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }
        
        IScalarFunction<3>&  coordTransfFunc(int i)         { return _func[i]; }
        IScalarFunction<3>&  inverseCoordTransfFunc(int i)  { return _funcInverse[i]; }
    };    

    class CoordTransfCartesianToCylindrical : public CoordTransf<Vector3Cartesian, Vector3Cylindrical, 3>
    {
        // q[0] = x
        // q[1] = y
        // q[2] = z
        public:
        static double func1(const VectorN<Real, 3> &q) { return sqrt(q[0]*q[0] + q[1]*q[1]); }
        static double func2(const VectorN<Real, 3> &q) { return atan2(q[1], q[0]); }
        static double func3(const VectorN<Real, 3> &q) { return q[2]; }

        // q1 = r   - distance from symmetry axis
        // q2 = phi - angle to symmetry axis
        // q3 = z   - z
        static double funcInverse1(const VectorN<Real, 3> &q) { return q[0] * cos(q[1]); }
        static double funcInverse2(const VectorN<Real, 3> &q) { return q[0] * sin(q[1]); }
        static double funcInverse3(const VectorN<Real, 3> &q) { return q[2]; }

        inline static ScalarFunction<3> _func[3] = { ScalarFunction<3>{func1},
                                                                ScalarFunction<3>{func2},
                                                                ScalarFunction<3>{func3}
                                                              };

        inline static ScalarFunction<3> _funcInverse[3] = { ScalarFunction<3>{funcInverse1},
                                                                       ScalarFunction<3>{funcInverse2},
                                                                       ScalarFunction<3>{funcInverse3}
                                                                     };

        Vector3Cylindrical transf(const Vector3Cartesian &q)            { return Vector3Cylindrical{ func1(q), func2(q), func3(q) }; }
        Vector3Cartesian   transfInverse(const Vector3Cylindrical &q)   { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }
        
        IScalarFunction<3>&  coordTransfFunc(int i)         { return _func[i]; }
        IScalarFunction<3>&  inverseCoordTransfFunc(int i)  { return _funcInverse[i]; }
    };

    static CoordTransfSphericalToCartesian      CoordTransfSpherToCart;
    static CoordTransfCylindricalToCartesian    CoordTransfCylToCart;
    static CoordTransfCartesianToSpherical      CoordTransfCartToSpher;
    static CoordTransfCartesianToCylindrical    CoordTransfCartToCyl;
}

///////////////////////////   ./include/basic_types/MetricTensor.h   ///////////////////////////



namespace MML
{
    template<int N>
    class MetricTensor : public ITensorField2<N>
    {
    public:
        void   ValueAtPoint(const VectorN<Real, N> &pos, ITensor2<N> &outVal) const 
        {
            for( int i=0; i<N; i++ )
                for( int j=0; j<N; j++ )
                    outVal.Component(i,j) = this->Component(i,j, pos);
        }
    };

    template<int N>
    class MetricTensorCartesian: public MetricTensor<N>
    {
        public:
        double Component(int i, int j, const VectorN<Real, N> &pos) const
        {
            if( i == j )
                return 1.0;
            else
                return 0.0;
        }
    };

    class MetricTensorSpherical: public MetricTensor<3>
    {
        public:
        virtual double Component(int i, int j, const VectorN<Real, 3> &pos) const
        {
            if( i == 0 && j == 0 )
                return 1.0;
            else if( i == 1 && j == 1 )
                return pos[0] * pos[0];
            else if( i == 2 && j == 2 )
                return pos[0] * pos[0] * sin(pos[1]) * sin(pos[1]);
            else
                return 0.0;
        }
    };

    class MetricTensorCylindrical: public MetricTensor<3>
    {
        public:
        virtual double Component(int i, int j, const VectorN<Real, 3> &pos) const
        {
            if( i == 0 && j == 0 )
                return 1.0;
            else if( i == 1 && j == 1 )
                return pos[0] * pos[0];
            else if( i == 2 && j == 2 )
                return 1.0;
            else
                return 0.0;
        }
    };

    template<typename VectorFrom, typename VectorTo, int N>
    class MetricTensorFromCoordTransf: public MetricTensor<N>
    {
        ICoordTransf<VectorFrom, VectorTo, N> &_coordTransf;

        public:
        MetricTensorFromCoordTransf(ICoordTransf<VectorFrom, VectorTo, N> &inTransf) : _coordTransf(inTransf)
        { }

        virtual double Component(int i, int j, const VectorN<Real, N> &pos) const
        {
            double g_ij = 0.0;
            for(int l=0; l<N; l++)
            {
                g_ij += Derivation::DerivePartial<N>(_coordTransf.coordTransfFunc(l), i, pos, nullptr) * Derivation::DerivePartial<N>(_coordTransf.coordTransfFunc(l), j, pos, nullptr);
            }
            return g_ij;
        }
    };
}
///////////////////////////   ./include/algorithms/FieldOperations.h   ///////////////////////////// grad
// - cart
// - spher
// - cyl







namespace MML
{
    // TODO - Vec. field op. - Laplacian
    class ScalarFieldOperations
    {
        public:
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////                   GRADIENT                     /////////////////////////////
        template<int N>
        static VectorN<Real, N> Gradient(IScalarFunction<N> &scalarField, const MetricTensor<N>& metricTensor, const VectorN<Real, N> &pos)
        {
            VectorN<Real, N> derivsAtPoint = Derivation::DerivePartialAll<N>(scalarField, pos, nullptr);
            
            Tensor2<N> metricAtPoint(2,0);
            metricTensor.ValueAtPoint(pos, metricAtPoint);

            VectorN<Real, N> ret = metricAtPoint * derivsAtPoint;

            return ret;
        }

        template<int N>
        static VectorN<Real, N> GradientCart(const IScalarFunction<N> &scalarField, const VectorN<Real, N> &pos)
        {
            return Derivation::DerivePartialAll<N>(scalarField, pos, nullptr);
        }

        static Vector3Spherical GradientSpher(const IScalarFunction<3> &scalarField, const Vector3Spherical &pos)
        {
            Vector3Spherical ret = Derivation::DerivePartialAll<3>(scalarField, pos, nullptr);

            ret[1] = ret[1] / pos[0];
            ret[2] = ret[2] / (pos[0] * sin(pos[1]));

            return ret;
        }

        static Vector3Cylindrical GradientCyl(const IScalarFunction<3> &scalarField, const Vector3Cylindrical &pos)
        {
            Vector3Cylindrical ret = Derivation::DerivePartialAll<3>(scalarField, pos, nullptr);

            ret[1] = ret[1] / pos[0];

            return ret;
        }            
    };

    class VectorFieldOperations
    {
    public:
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////                  DIVERGENCE                    /////////////////////////////
        template<int N>
        static double DivCart(const IVectorFunction<N> &vectorField, const VectorN<Real, N> &pos)
        {
            double div = 0.0;
            for( int i=0; i<N; i++ )
                div += Derivation::DeriveVecPartial<N>(vectorField, i, i, pos, nullptr);

            return div;
        }    

        static double DivSpher(const IVectorFunction<3> &vectorField, const VectorN<Real, 3> &x)
        {
            VectorN<Real, 3> vals = vectorField(x);

            VectorN<Real, 3> derivs;
            for( int i=0; i<3; i++ )
                derivs[i] = Derivation::DeriveVecPartial<3>(vectorField, i, i, x, nullptr);
            
            double div = 0.0;
            div += 1 / (x[0]*x[0]) * (2 * x[0] * vals[0] + x[0]*x[0] * derivs[0]);
            div += 1 / (x[0] * sin(x[1])) * (cos(x[1]) * vals[1] + sin(x[1]) * derivs[1]);
            div += 1 / (x[0] * sin(x[1])) * derivs[2];

            return div;
        }           

        static double DivCyl(const IVectorFunction<3> &vectorField, const VectorN<Real, 3> &x)
        {
            VectorN<Real, 3> vals = vectorField(x);

            VectorN<Real, 3> derivs;
            for( int i=0; i<3; i++ )
                derivs[i] = Derivation::DeriveVecPartial<3>(vectorField, i, i, x, nullptr);
            
            double div = 0.0;
            div += 1 / x[0] * (vals[0] + x[0] * derivs[0]);
            div += 1 / x[0] * derivs[1];
            div += derivs[2];

            return div;
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////                     CURL                       /////////////////////////////
        static Vector3Cartesian CurlCart(const IVectorFunction<3> &vectorField, const VectorN<Real, 3> &pos)
        {
            double dzdy = Derivation::DeriveVecPartial<3>(vectorField, 2, 1, pos, nullptr);
            double dydz = Derivation::DeriveVecPartial<3>(vectorField, 1, 2, pos, nullptr);

            double dxdz = Derivation::DeriveVecPartial<3>(vectorField, 0, 2, pos, nullptr);
            double dzdx = Derivation::DeriveVecPartial<3>(vectorField, 2, 0, pos, nullptr);

            double dydx = Derivation::DeriveVecPartial<3>(vectorField, 1, 0, pos, nullptr);
            double dxdy = Derivation::DeriveVecPartial<3>(vectorField, 0, 1, pos, nullptr);

            Vector3Cartesian curl{dzdy - dydz, dxdz - dzdx, dydx - dxdy};

            return curl;
        }    

        static Vector3Spherical CurlSpher(const IVectorFunction<3> &vectorField, const VectorN<Real, 3> &pos)
        {
            VectorN<Real, 3> vals = vectorField(pos);

            double dphidtheta = Derivation::DeriveVecPartial<3>(vectorField, 2, 1, pos, nullptr);
            double dthetadphi = Derivation::DeriveVecPartial<3>(vectorField, 1, 2, pos, nullptr);

            double drdphi = Derivation::DeriveVecPartial<3>(vectorField, 0, 2, pos, nullptr);
            double dphidr = Derivation::DeriveVecPartial<3>(vectorField, 2, 0, pos, nullptr);

            double dthetadr = Derivation::DeriveVecPartial<3>(vectorField, 1, 0, pos, nullptr);
            double drdtheta = Derivation::DeriveVecPartial<3>(vectorField, 0, 1, pos, nullptr);

            Vector3Spherical ret;
            const double &r     = pos[0];
            const double &theta = pos[1];
            const double &phi   = pos[2];

            ret[0] = 1 / (r * sin(theta)) * (cos(theta) * vals[2] + sin(theta) * dphidtheta - dthetadphi);
            ret[1] = 1 / r * (1 / sin(theta)  * drdphi - vals[2] - r * dphidr);
            ret[2] = 1 / r * (vals[1] + r * dthetadr - drdtheta);

            return ret;
        }           

        static Vector3Cylindrical CurlCyl(const IVectorFunction<3> &vectorField, const VectorN<Real, 3> &pos)
        {
            VectorN<Real, 3> vals = vectorField(pos);

            double dzdphi = Derivation::DeriveVecPartial<3>(vectorField, 2, 1, pos, nullptr);
            double dphidz = Derivation::DeriveVecPartial<3>(vectorField, 1, 2, pos, nullptr);

            double drdz = Derivation::DeriveVecPartial<3>(vectorField, 0, 2, pos, nullptr);
            double dzdr = Derivation::DeriveVecPartial<3>(vectorField, 2, 0, pos, nullptr);

            double dphidr = Derivation::DeriveVecPartial<3>(vectorField, 1, 0, pos, nullptr);
            double drdphi = Derivation::DeriveVecPartial<3>(vectorField, 0, 1, pos, nullptr);

            Vector3Cylindrical ret{1.0 / pos[0] * dzdphi - dphidz, drdz - dzdr, 1 / pos[0] * (vals[1] + pos[0] * dphidr - drdphi)};

            return ret;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////                  LAPLACIAN                     /////////////////////////////
        template<int N>
        static VectorN<Real, N> LaplacianCart(const IScalarFunction<N> &scalarField, const VectorN<Real, N> &pos)
        {
            double lapl = 0.0;
            for( int i=0; i<N; i++ )
                lapl += Derivation::NSecDer1Partial<N>(scalarField, i, i, pos, nullptr);

            return div;
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
            const IParametricCurve<N> &_curve;
        public:
            HelperCurveLen(const IParametricCurve<N> &curve) : _curve(curve) {}

            double operator()(double t) const 
            {
                auto tangent_vec = Derivation::DeriveCurve<N>(_curve, t, nullptr);
                return tangent_vec.NormL2();
            }
        };
        template<int N>
        class HelperWorkIntegral : public IRealFunction
        {
            const IScalarFunction<N>  &_potential;
            const IParametricCurve<N> &_curve;
        public:
            HelperWorkIntegral(const IScalarFunction<N> &potentialField, const IParametricCurve<N> &curve) : _potential(potentialField), _curve(curve) {}

            double operator()(double t) const 
            {
                auto tangent_vec = Derivation::DeriveCurve<N>(_curve, t, nullptr);
                auto gradient    = ScalarFieldOperations::GradientCart(_potential, _curve(t));

                return tangent_vec.ScalarProductCartesian(gradient);
            }
        };
        template<int N>
        class HelperLineIntegral : public IRealFunction
        {
            const IVectorFunction<N>  &_vector_field;
            const IParametricCurve<N> &_curve;
        public:
            HelperLineIntegral(const IVectorFunction<N> &vectorField, const IParametricCurve<N> &curve) : _vector_field(vectorField), _curve(curve) {}

            double operator()(double t) const 
            {
                auto tangent_vec = Derivation::DeriveCurve<N>(_curve, t, nullptr);
                auto field_vec   = _vector_field(_curve(t));

                return tangent_vec.ScalarProductCartesian(field_vec);
            }
        };        
    public:            
        template<int N>
        static double ParametricCurveLength(const IParametricCurve<N> &curve, const double a, const double b)
        {
            HelperCurveLen helper(curve);
            
            return Integration::IntegrateTrap(helper, a, b);
        }

        static double WorkIntegral(const IScalarFunction<3> &potentialField, const IParametricCurve<3> &curve, const double a, const double b, const double eps=Defaults::WorkIntegralPrecision)
        {
            HelperWorkIntegral helper(potentialField, curve);
            
            return Integration::IntegrateTrap(helper, a, b, eps);
        }

        static double LineIntegral(const IVectorFunction<3> &vectorField, const IParametricCurve<3> &curve, const double a, const double b, const double eps=Defaults::LineIntegralPrecision)
        {
            HelperLineIntegral helper(vectorField, curve);
            
            return Integration::IntegrateTrap(helper, a, b, eps);            
        }

        static double SurfaceIntegral(const IVectorFunction<3> &vectorField, const IParametricSurface<3> &surface, const double x1, const double x2, const double y1, const double y2)
        {
            return 0.0;
        }
	};
} // end namespace
///////////////////////////   ./include/algorithms/LinAlgEqSolvers.h   ///////////////////////////


namespace MML
{
    class GaussJordanSolver
    {
        public:

        static bool Solve(Matrix<Real> &a, Matrix<Real> &b)
        {
            int i,icol,irow,j,k,l,ll;
            double big,dum,pivinv;

            int n=a.RowNum();
            int m=b.ColNum();
            std::vector<int> indxc(n),indxr(n),ipiv(n);
            for (j=0;j<n;j++) ipiv[j]=0;
            for (i=0;i<n;i++) {
                big=0.0;
                for (j=0;j<n;j++)
                    if (ipiv[j] != 1)
                        for (k=0;k<n;k++) {
                            if (ipiv[k] == 0) {
                                if (std::abs(a[j][k]) >= big) {
                                    big=std::abs(a[j][k]);
                                    irow=j;
                                    icol=k;
                                }
                            }
                        }
                ++(ipiv[icol]);
                if (irow != icol) {
                    for (l=0;l<n;l++) std::swap(a[irow][l],a[icol][l]);
                    for (l=0;l<m;l++) std::swap(b[irow][l],b[icol][l]);
                }
                indxr[i]=irow;
                indxc[i]=icol;

                if (a[icol][icol] == 0.0) 
                    throw SingularMatrixError("GaussJordanSolver::Solve - Singular Matrix");

                pivinv=1.0/a[icol][icol];
                a[icol][icol]=1.0;
                for (l=0;l<n;l++) a[icol][l] *= pivinv;
                for (l=0;l<m;l++) b[icol][l] *= pivinv;
                for (ll=0;ll<n;ll++)
                    if (ll != icol) {
                        dum=a[ll][icol];
                        a[ll][icol]=0.0;
                        for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
                        for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
                    }
            }
            for (l=n-1;l>=0;l--) {
                if (indxr[l] != indxc[l])
                    for (k=0;k<n;k++)
                        std::swap(a[k][indxr[l]],a[k][indxc[l]]);
            }

            return true;
        }

        static bool Solve(Matrix<Real> &a, Vector<Real> &b)
        {
            Matrix<Real> bmat = Matrix<Real>::ColumnMatrixFromVector(b);
            return Solve(a, bmat);
        }
    };

    class LUDecompositionSolver
    {
    private:
        int n;
        Matrix<Real> &refOrig;

        Matrix<Real> lu;
        std::vector<int> indx;
        double d;
    
    public:
        LUDecompositionSolver(Matrix<Real>  &a) : n(a.RowNum()), refOrig(a), lu(a), indx(n) 
        {
            // Given a Matrix<Real> a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise
            // permutation of itself. a and n are input. a is output, arranged as in equation (NR 2.3.14);
            // indx[1..n] is an output Vector<Real> that records the row permutation effected by the partial
            // pivoting; d is output as ±1 depending on whether the number of row interchanges was even
            // or odd, respectively. This routine is used in combination with lubksb to solve linear equations
            // or invert a Matrix<Real>.
            const double TINY=1.0e-40;
            int i,imax,j,k;
            double big,temp;
            Vector<Real> vv(n);
            d=1.0;
            for (i=0;i<n;i++) {
                big=0.0;
                for (j=0;j<n;j++)
                    if ((temp=std::abs(lu[i][j])) > big) big=temp;
                if (big == 0.0) 
                    throw SingularMatrixError("LUDecompositionSolver::ctor - Singular Matrix<Real>");

                vv[i]=1.0/big;
            }
            for (k=0;k<n;k++) {
                big=0.0;
                imax=k;
                for (i=k;i<n;i++) {
                    temp=vv[i]*std::abs(lu[i][k]);
                    if (temp > big) {
                        big=temp;
                        imax=i;
                    }
                }
                if (k != imax) {
                    for (j=0;j<n;j++) {
                        temp=lu[imax][j];
                        lu[imax][j]=lu[k][j];
                        lu[k][j]=temp;
                    }
                    d = -d;
                    vv[imax]=vv[k];
                }
                indx[k]=imax;
                if (lu[k][k] == 0.0) lu[k][k]=TINY;
                for (i=k+1;i<n;i++) {
                    temp=lu[i][k] /= lu[k][k];
                    for (j=k+1;j<n;j++)
                        lu[i][j] -= temp*lu[k][j];
                }
            }
        }
        void Solve(Vector<Real> &b, Vector<Real> &x)
        {
            // Solves the set of n linear equations A·X = B. Here a[1..n][1..n] is input, not as the Matrix<Real>
            // A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
            // as the permutation Vector<Real> returned by ludcmp. b[1..n] is input as the right-hand side Vector<Real>
            // B, and returns with the solution Vector<Real> X. a, n, and indx are not modified by this routine
            // and can be left in place for successive calls with different right-hand sides b. This routine takes
            // into account the possibility that b will begin with many zero elements, so it is efficient for use
            // in Matrix<Real> inversion
            int i,ii=0,ip,j;
            double sum;
            if (b.size() != n || x.size() != n)
                throw("LUdcmp::solve bad sizes");
            for (i=0;i<n;i++) 
                x[i] = b[i];
            for (i=0;i<n;i++) {
                ip=indx[i];
                sum=x[ip];
                x[ip]=x[i];
                if (ii != 0)
                    for (j=ii-1;j<i;j++) sum -= lu[i][j]*x[j];
                else if (sum != 0.0)
                    ii=i+1;
                x[i]=sum;
            }
            for (i=n-1;i>=0;i--) {
                sum=x[i];
                for (j=i+1;j<n;j++) sum -= lu[i][j]*x[j];
                x[i]=sum/lu[i][i];
            }
        }

        void Solve(Matrix<Real> &b, Matrix<Real> &x)
        {
            int i,j,m=b.ColNum();
            
            if (b.RowNum() != n || x.RowNum() != n || b.ColNum() != x.ColNum())
                throw("LUdcmp::solve bad sizes");
            
            Vector<Real> xx(n);
            
            for (j=0;j<m;j++) {
                for (i=0;i<n;i++) 
                    xx[i] = b[i][j];
                
                Solve(xx,xx);
                
                for (i=0;i<n;i++) 
                    x[i][j] = xx[i];
            }
        }        

        // Using the stored LU decomposition, return in ainv the matrix inverse 
        void inverse(Matrix<Real> &ainv)
        {
            int i,j;
            ainv.Resize(n,n);
            for (i=0;i<n;i++) {
                for (j=0;j<n;j++) ainv[i][j] = 0.;
                ainv[i][i] = 1.;
            }
            Solve(ainv,ainv);
        }
        double det()
        {
            double dd = d;
            for (int i=0;i<n;i++) 
                dd *= lu[i][i];
            return dd;
        }
        
        // Improves a solution Vector<Real> x[1..n] of the linear set of equations A · X = B. The Matrix<Real>
        // a[1..n][1..n], and the Vector<Real>s b[1..n] and x[1..n] are input, as is the dimension n.
        // Also input is alud[1..n][1..n], the LU decomposition of a as returned by ludcmp, and
        // the Vector<Real> indx[1..n] also returned by that routine. On output, only x[1..n] is modified,
        // to an improved set of values
        void mprove(Vector<Real> &b, Vector<Real> &x)
        {
            int i,j;
            Vector<Real> r(n);
            
            for (i=0;i<n;i++) {
                long double  sdp = -b[i];
                for (j=0;j<n;j++)
                    sdp += (long double)refOrig[i][j] * (long double)x[j];
                r[i]=sdp;
            }
            
            Solve(r,r);
            
            for (i=0;i<n;i++) 
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
        CholeskyDecompositionSolver(Matrix<Real> &a) : n(a.RowNum()), el(a) 
        {
            // Given a positive-definite symmetric Matrix<Real> a[1..n][1..n], this routine constructs its Cholesky
            // decomposition, A = L · LT . On input, only the upper triangle of a need be given; it is not
            // modified. The Cholesky factor L is returned in the lower triangle of a, except for its diagonal
            // elements which are returned in p[1..n]
            int i,j,k;
            double sum;
            
            if (el.ColNum() != n) 
                throw("need square Matrix<Real>");
            
            for (i=0;i<n;i++) {
                for (j=i;j<n;j++) {
                    for (sum=el[i][j],k=i-1;k>=0;k--) sum -= el[i][k]*el[j][k];
                    if (i == j) {
                        if (sum <= 0.0)
                            throw("Cholesky failed");
                        el[i][i]=sqrt(sum);
                    } else el[j][i]=sum/el[i][i];
                }
            }
            for (i=0;i<n;i++) for (j=0;j<i;j++) el[j][i] = 0.;
        }
        void Solve(Vector<Real> &b, Vector<Real> &x) 
        {
            // Solves the set of n linear equations A · x = b, where a is a positive-definite symmetric Matrix<Real>.
            // a[1..n][1..n] and p[1..n] are input as the output of the routine choldc. Only the lower
            // triangle of a is accessed. b[1..n] is input as the right-hand side Vector<Real>. The solution Vector<Real> is
            // returned in x[1..n]. a, n, and p are not modified and can be left in place for successive calls
            // with different right-hand sides b. b is not modified unless you identify b and x in the calling
            // sequence, which is allowed.
            int i,k;
            double sum;
            if (b.size() != n || x.size() != n) throw("bad lengths in Cholesky");
            for (i=0;i<n;i++) {
                for (sum=b[i],k=i-1;k>=0;k--) sum -= el[i][k]*x[k];
                x[i]=sum/el[i][i];
            }
            for (i=n-1;i>=0;i--) {
                for (sum=x[i],k=i+1;k<n;k++) sum -= el[k][i]*x[k];
                x[i]=sum/el[i][i];
            }		
        }
        void elmult(Vector<Real> &y, Vector<Real> &b) {
            int i,j;
            if (b.size() != n || y.size() != n) throw("bad lengths");
            for (i=0;i<n;i++) {
                b[i] = 0.;
                for (j=0;j<=i;j++) b[i] += el[i][j]*y[j];
            }
        }
        void elsolve(Vector<Real> &b, Vector<Real> &y) {
            int i,j;
            double sum;
            if (b.size() != n || y.size() != n) throw("bad lengths");
            for (i=0;i<n;i++) {
                for (sum=b[i],j=0; j<i; j++) sum -= el[i][j]*y[j];
                y[i] = sum/el[i][i];
            }
        }
        void inverse(Matrix<Real> &ainv) {
            int i,j,k;
            double sum;
            ainv.Resize(n,n);
            for (i=0;i<n;i++) for (j=0;j<=i;j++){
                sum = (i==j? 1. : 0.);
                for (k=i-1;k>=j;k--) sum -= el[i][k]*ainv[j][k];
                ainv[j][i]= sum/el[i][i];
            }
            for (i=n-1;i>=0;i--) for (j=0;j<=i;j++){
                sum = (i<j? 0. : ainv[j][i]);
                for (k=i+1;k<n;k++) sum -= el[k][i]*ainv[j][k];
                ainv[i][j] = ainv[j][i] = sum/el[i][i];
            }				
        }
        double logdet() {
            double sum = 0.;
            for (int  i=0; i<n; i++) sum += log(el[i][i]);
            return 2.*sum;
        }
    };

    ///////////////////////   QR DECOMPOSITION    /////////////////////////////
    class QRDecompositionSolver
    {
    private:
        int n;
        Matrix<Real> qt, r;
        bool sing;    

    public:

        QRDecompositionSolver(Matrix<Real> &a) : n(a.RowNum()), qt(n,n), r(a), sing(false) 
        {
            // Constructs the QR decomposition of a[1..n][1..n]. The upper triangular Matrix<Real> R is returned in the upper triangle of a, 
            // except for the diagonal elements of R which are returned in d[1..n]. 
            int i,j,k;
            Vector<Real> c(n), d(n);
            double scale,sigma,sum,tau;
            for (k=0;k<n-1;k++) {
                scale=0.0;
                for (i=k;i<n;i++) scale=std::max(scale,std::abs(r[i][k]));
                if (scale == 0.0) {
                    sing=true;
                    c[k]=d[k]=0.0;
                } else {
                    for (i=k;i<n;i++) r[i][k] /= scale;
                    for (sum=0.0,i=k;i<n;i++) sum += SQR(r[i][k]);
                    sigma=SIGN(sqrt(sum),r[k][k]);
                    r[k][k] += sigma;
                    c[k]=sigma*r[k][k];
                    d[k] = -scale*sigma;
                    for (j=k+1;j<n;j++) {
                        for (sum=0.0,i=k;i<n;i++) sum += r[i][k]*r[i][j];
                        tau=sum/c[k];
                        for (i=k;i<n;i++) r[i][j] -= tau*r[i][k];
                    }
                }
            }
            d[n-1]=r[n-1][n-1];
            if (d[n-1] == 0.0) sing=true;
            for (i=0;i<n;i++) {
                for (j=0;j<n;j++) qt[i][j]=0.0;
                qt[i][i]=1.0;
            }
            for (k=0;k<n-1;k++) {
                if (c[k] != 0.0) {
                    for (j=0;j<n;j++) {
                        sum=0.0;
                        for (i=k;i<n;i++)
                            sum += r[i][k]*qt[i][j];
                        sum /= c[k];
                        for (i=k;i<n;i++)
                            qt[i][j] -= sum*r[i][k];
                    }
                }
            }
            for (i=0;i<n;i++) {
                r[i][i]=d[i];
                for (j=0;j<i;j++) r[i][j]=0.0;
            }
        }

        // Solves the set of n linear equations A · x = b. a[1..n][1..n], c[1..n], and d[1..n] are
        // input as the output of the routine qrdcmp and are not modified. b[1..n] is input as the
        // right-hand side Vector<Real>, and is overwritten with the solution Vector<Real> on output. 
        void Solve(Vector<Real> &b, Vector<Real> &x) 
        {           
            qtmult(b,x);
            rsolve(x,x);
        }

        void qtmult(Vector<Real> &b, Vector<Real> &x) {
            int i,j;
            double sum;
            for (i=0;i<n;i++) {
                sum = 0.;
                for (j=0;j<n;j++) sum += qt[i][j]*b[j];
                x[i] = sum;
            }
        }

        void rsolve(Vector<Real> &b, Vector<Real> &x) 
        {
            // Solves the set of n linear equations R · x = b, where R is an upper triangular Matrix<Real> stored in
            // a and d. a[1..n][1..n] and d[1..n] are input as the output of the routine qrdcmp and
            // are not modified. b[1..n] is input as the right-hand side Vector<Real>, and is overwritten with the
            // solution Vector<Real> on output            
            int i,j;
            double sum;
            if (sing) 
                throw SingularMatrixError("QRDecompositionSolver::rsolve - attempting solve in a singular QR");

            for (i=n-1;i>=0;i--) {
                sum=b[i];
                for (j=i+1;j<n;j++) sum -= r[i][j]*x[j];
                x[i]=sum/r[i][i];
            }
        }
        void update(Vector<Real> &u, Vector<Real> &v) 
        {
            // Given the QR decomposition of some n × n Matrix<Real>, calculates the QR decomposition of the
            // Matrix<Real> Q·(R+ u x v). The quantities are dimensioned as r[1..n][1..n], qt[1..n][1..n],
            // u[1..n], and v[1..n]. Note that QT is input and returned in qt.            
            int i,k;
            Vector<Real> w(u);
            for (k=n-1;k>=0;k--)
                if (w[k] != 0.0) break;
            if (k < 0) k=0;
            for (i=k-1;i>=0;i--) {
                rotate(i,w[i],-w[i+1]);
                if (w[i] == 0.0)
                    w[i]=std::abs(w[i+1]);
                else if (std::abs(w[i]) > std::abs(w[i+1]))
                    w[i]=std::abs(w[i])*sqrt(1.0+SQR(w[i+1]/w[i]));
                else w[i]=std::abs(w[i+1])*sqrt(1.0+SQR(w[i]/w[i+1]));
            }
            for (i=0;i<n;i++) r[0][i] += w[0]*v[i];
            for (i=0;i<k;i++)
                rotate(i,r[i][i],-r[i+1][i]);
            for (i=0;i<n;i++)
                if (r[i][i] == 0.0) sing=true;
        }

        void rotate(const int i, const double a, const double b)
        {
            // Given matrices r[1..n][1..n] and qt[1..n][1..n], carry out a Jacobi rotation on rows
            // i and i + 1 of each Matrix<Real>. a and b are the parameters of the rotation: cos phi = a=pa2 + b2,
            // sin phi = b=pa2 + b2.            
            int j;
            double c,fact,s,w,y;
            if (a == 0.0) {
                c=0.0;
                s=(b >= 0.0 ? 1.0 : -1.0);
            } else if (std::abs(a) > std::abs(b)) {
                fact=b/a;
                c=SIGN(1.0/sqrt(1.0+(fact*fact)),a);
                s=fact*c;
            } else {
                fact=a/b;
                s=SIGN(1.0/sqrt(1.0+(fact*fact)),b);
                c=fact*s;
            }
            for (j=i;j<n;j++) {
                y=r[i][j];
                w=r[i+1][j];
                r[i][j]=c*y-s*w;
                r[i+1][j]=s*y+c*w;
            }
            for (j=0;j<n;j++) {
                y=qt[i][j];
                w=qt[i+1][j];
                qt[i][j]=c*y-s*w;
                qt[i+1][j]=s*y+c*w;
            }
        }    

    };

    // void tridag(Vector<Real> &a, Vector<Real> &b, Vector<Real> &c, Vector<Real> &r, Vector<Real> &u)
    // {
    //     // Solves for a Vector<Real> u[1..n] the tridiagonal linear set given by equation (2.4.1). a[1..n],
    //     // b[1..n], c[1..n], and r[1..n] are input Vector<Real>s and are not modified.

    //     int j;
    //     double bet;

    //     int n=(int)a.size();
    //     Vector<Real> gam(n);

    //     if (b[0] == 0.0) 
    //         //nrerror("Error 1 in tridag");
    //         return;

    //     u[0]=r[0]/(bet=b[0]);
    //     for (j=1;j<n;j++) {
    //         gam[j]=c[j-1]/bet;
    //         bet=b[j]-a[j]*gam[j];
    //         if (bet == 0.0) 
    //             // nrerror("Error 2 in tridag");
    //             return;

    //         u[j]=(r[j]-a[j]*u[j-1])/bet;
    //     }
    //     for (j=(n-2);j>=0;j--)
    //         u[j] -= gam[j+1]*u[j+1];
    // }

    /////////////////////////////////   SVD DECOMPOSITION      /////////////////////////////
    class SVDecompositionSolver 
    {
    private:
        int m,n;
        Matrix<Real> u,v;
        Vector<Real> w;
        double eps, tsh;
    
    public:
        SVDecompositionSolver(Matrix<Real> &a) : m(a.RowNum()), n(a.ColNum()), u(a), v(n,n), w(n) 
        {
            // Given a Matrix<Real> a[1..m][1..n], this routine computes its singular value decomposition, A = U·W ·V T . 
            // The Matrix<Real> U replaces a on output. 
            // The diagonal Matrix<Real> of singular values W is output as a Vector<Real> w[1..n]. 
            // The Matrix<Real> V (not the transpose V T ) is output as v[1..n][1..n].            
            eps = std::numeric_limits<double>::epsilon();
            decompose();
            reorder();
            tsh = 0.5*sqrt(m+n+1.)*w[0]*eps;
        }

        double inv_condition() {
            return (w[0] <= 0. || w[n-1] <= 0.) ? 0. : w[n-1]/w[0];
        }
        
        void solve(Vector<Real> &b, Vector<Real> &x, double thresh = -1.) 
        {
            // Solve A  x D b for a vector x using the pseudoinverse of A as obtained by SVD. If positive,
            // thresh is the threshold value below which singular values are considered as zero. If thresh is
            // negative, a default based on expected roundoff error is used.
            int i,j,jj;
            double s;
            if (b.size() != m || x.size() != n) throw("solve bad sizes");
            Vector<Real> tmp(n);
            tsh = (thresh >= 0. ? thresh : 0.5*sqrt(m+n+1.)*w[0]*eps);
            for (j=0;j<n;j++) {
                s=0.0;
                if (w[j] > tsh) {
                    for (i=0;i<m;i++) s += u[i][j]*b[i];
                    s /= w[j];
                }
                tmp[j]=s;
            }
            for (j=0;j<n;j++) {
                s=0.0;
                for (jj=0;jj<n;jj++) s += v[j][jj]*tmp[jj];
                x[j]=s;
            }
        }

        // Solves m sets of n equations A  X D B using the pseudoinverse of A. The right-hand sides are
        // input as b[0..n-1][0..m-1], while x[0..n-1][0..m-1] returns the solutions. thresh as above.
        void solve(Matrix<Real> &b, Matrix<Real> &x, double thresh = -1.)
        {
            int i,j,p=b.ColNum();
            if (b.RowNum() != m || x.RowNum() != n || x.ColNum() != p)
                throw("solve bad sizes");
            Vector<Real> xx(n),bcol(m);
            for (j=0;j<p;j++) {
                for (i=0;i<m;i++) bcol[i] = b[i][j];
                solve(bcol,xx,thresh);
                for (i=0;i<n;i++) x[i][j] = xx[i];
            }
        }

        // Return the rank of A, after zeroing any singular values smaller than thresh. If thresh is
        // negative, a default value based on estimated roundoff is used.        
        int rank(double thresh = -1.) {
            int j,nr=0;
            tsh = (thresh >= 0. ? thresh : 0.5*sqrt(m+n+1.)*w[0]*eps);
            for (j=0;j<n;j++) if (w[j] > tsh) nr++;
            return nr;
        }

        // Return the nullity of A, after zeroing any singular values smaller than thresh. Default value as above.
        int nullity(double thresh = -1.) {
            int j,nn=0;
            tsh = (thresh >= 0. ? thresh : 0.5*sqrt(m+n+1.)*w[0]*eps);
            for (j=0;j<n;j++) if (w[j] <= tsh) nn++;
            return nn;
        }

        // Give an orthonormal basis for the range of A as the columns of a returned matrix. thresh as above.
        Matrix<Real> range(double thresh = -1.){
            int i,j,nr=0;
            Matrix<Real> rnge(m,rank(thresh));
            for (j=0;j<n;j++) {
                if (w[j] > tsh) {
                    for (i=0;i<m;i++) rnge[i][nr] = u[i][j];
                    nr++;
                }
            }
            return rnge;
        }

        // Give an orthonormal basis for the nullspace of A as the columns of a returned matrix. thresh as above
        Matrix<Real> nullspace(double thresh = -1.){
            int j,jj,nn=0;
            Matrix<Real> nullsp(n,nullity(thresh));
            for (j=0;j<n;j++) {
                if (w[j] <= tsh) {
                    for (jj=0;jj<n;jj++) nullsp[jj][nn] = v[jj][j];
                    nn++;
                }
            }
            return nullsp;
        }
        void decompose() {
            bool flag;
            int i,its,j,jj,k,l,nm;
            double anorm,c,f,g,h,s,scale,x,y,z;
            Vector<Real> rv1(n);
            g = scale = anorm = 0.0;
            for (i=0;i<n;i++) {
                l=i+2;
                rv1[i]=scale*g;
                g=s=scale=0.0;
                if (i < m) {
                    for (k=i;k<m;k++) scale += std::abs(u[k][i]);
                    if (scale != 0.0) {
                        for (k=i;k<m;k++) {
                            u[k][i] /= scale;
                            s += u[k][i]*u[k][i];
                        }
                        f=u[i][i];
                        g = -SIGN(sqrt(s),f);
                        h=f*g-s;
                        u[i][i]=f-g;
                        for (j=l-1;j<n;j++) {
                            for (s=0.0,k=i;k<m;k++) s += u[k][i]*u[k][j];
                            f=s/h;
                            for (k=i;k<m;k++) u[k][j] += f*u[k][i];
                        }
                        for (k=i;k<m;k++) u[k][i] *= scale;
                    }
                }
                w[i]=scale *g;
                g=s=scale=0.0;
                if (i+1 <= m && i+1 != n) {
                    for (k=l-1;k<n;k++) scale += std::abs(u[i][k]);
                    if (scale != 0.0) {
                        for (k=l-1;k<n;k++) {
                            u[i][k] /= scale;
                            s += u[i][k]*u[i][k];
                        }
                        f=u[i][l-1];
                        g = -SIGN(sqrt(s),f);
                        h=f*g-s;
                        u[i][l-1]=f-g;
                        for (k=l-1;k<n;k++) rv1[k]=u[i][k]/h;
                        for (j=l-1;j<m;j++) {
                            for (s=0.0,k=l-1;k<n;k++) s += u[j][k]*u[i][k];
                            for (k=l-1;k<n;k++) u[j][k] += s*rv1[k];
                        }
                        for (k=l-1;k<n;k++) u[i][k] *= scale;
                    }
                }
                anorm=std::max(anorm,(std::abs(w[i])+std::abs(rv1[i])));
            }
            for (i=n-1;i>=0;i--) {
                if (i < n-1) {
                    if (g != 0.0) {
                        for (j=l;j<n;j++)
                            v[j][i]=(u[i][j]/u[i][l])/g;
                        for (j=l;j<n;j++) {
                            for (s=0.0,k=l;k<n;k++) s += u[i][k]*v[k][j];
                            for (k=l;k<n;k++) v[k][j] += s*v[k][i];
                        }
                    }
                        for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
                    }
                    v[i][i]=1.0;
                    g=rv1[i];
                    l=i;
                }
                for (i=std::min(m,n)-1;i>=0;i--) {
                    l=i+1;
                    g=w[i];
                    for (j=l;j<n;j++) u[i][j]=0.0;
                    if (g != 0.0) {
                        g=1.0/g;
                        for (j=l;j<n;j++) {
                            for (s=0.0,k=l;k<m;k++) s += u[k][i]*u[k][j];
                            f=(s/u[i][i])*g;
                            for (k=i;k<m;k++) u[k][j] += f*u[k][i];
                        }
                        for (j=i;j<m;j++) u[j][i] *= g;
                    } else for (j=i;j<m;j++) u[j][i]=0.0;
                    ++u[i][i];
                }
                for (k=n-1;k>=0;k--) {
                    for (its=0;its<30;its++) {
                        flag=true;
                        for (l=k;l>=0;l--) {
                            nm=l-1;
                            if (l == 0 || std::abs(rv1[l]) <= eps*anorm) {
                                flag=false;
                                break;
                            }
                            if (std::abs(w[nm]) <= eps*anorm) break;
                        }
                        if (flag) {
                            c=0.0;
                            s=1.0;
                            for (i=l;i<k+1;i++) {
                                f=s*rv1[i];
                                rv1[i]=c*rv1[i];
                                if (std::abs(f) <= eps*anorm) break;
                                g=w[i];
                                h=pythag(f,g);
                                w[i]=h;
                                h=1.0/h;
                                c=g*h;
                                s = -f*h;
                                for (j=0;j<m;j++) {
                                    y=u[j][nm];
                                    z=u[j][i];
                                    u[j][nm]=y*c+z*s;
                                    u[j][i]=z*c-y*s;
                                }
                            }
                        }
                        z=w[k];
                        if (l == k) {
                            if (z < 0.0) {
                                w[k] = -z;
                                for (j=0;j<n;j++) v[j][k] = -v[j][k];
                            }
                            break;
                        }
                        if (its == 29) throw("no convergence in 30 svdcmp iterations");
                        x=w[l];
                        nm=k-1;
                        y=w[nm];
                        g=rv1[nm];
                        h=rv1[k];
                        f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
                        g=pythag(f,1.0);
                        f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
                        c=s=1.0;
                        for (j=l;j<=nm;j++) {
                            i=j+1;
                            g=rv1[i];
                            y=w[i];
                            h=s*g;
                            g=c*g;
                            z=pythag(f,h);
                            rv1[j]=z;
                            c=f/z;
                            s=h/z;
                            f=x*c+g*s;
                            g=g*c-x*s;
                            h=y*s;
                            y *= c;
                            for (jj=0;jj<n;jj++) {
                                x=v[jj][j];
                                z=v[jj][i];
                                v[jj][j]=x*c+z*s;
                                v[jj][i]=z*c-x*s;
                            }
                            z=pythag(f,h);
                            w[j]=z;
                            if (z) {
                                z=1.0/z;
                                c=f*z;
                                s=h*z;
                            }
                            f=c*g+s*y;
                            x=c*y-s*g;
                            for (jj=0;jj<m;jj++) {
                                y=u[jj][j];
                                z=u[jj][i];
                                u[jj][j]=y*c+z*s;
                                u[jj][i]=z*c-y*s;
                            }
                        }
                        rv1[l]=0.0;
                        rv1[k]=f;
                        w[k]=x;
                    }
                }
            }

        void reorder() {
            int i,j,k,s,inc=1;
            double sw;
            Vector<Real> su(m), sv(n);
            do { inc *= 3; inc++; } while (inc <= n);
            do {
                inc /= 3;
                for (i=inc;i<n;i++) {
                    sw = w[i];
                    for (k=0;k<m;k++) su[k] = u[k][i];
                    for (k=0;k<n;k++) sv[k] = v[k][i];
                    j = i;
                    while (w[j-inc] < sw) {
                        w[j] = w[j-inc];
                        for (k=0;k<m;k++) u[k][j] = u[k][j-inc];
                        for (k=0;k<n;k++) v[k][j] = v[k][j-inc];
                        j -= inc;
                        if (j < inc) break;
                    }
                    w[j] = sw;
                    for (k=0;k<m;k++) u[k][j] = su[k];
                    for (k=0;k<n;k++) v[k][j] = sv[k];

                }
            } while (inc > 1);
            for (k=0;k<n;k++) {
                s=0;
                for (i=0;i<m;i++) if (u[i][k] < 0.) s++;
                for (j=0;j<n;j++) if (v[j][k] < 0.) s++;
                if (s > (m+n)/2) {
                    for (i=0;i<m;i++) u[i][k] = -u[i][k];
                    for (j=0;j<n;j++) v[j][k] = -v[j][k];
                }
            }
        }

        double pythag(const double a, const double b) {
            double absa=std::abs(a), absb=std::abs(b);
            return (absa > absb ? absa*sqrt(1.0+SQR(absb/absa)) :
                (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb))));
        }
    };
} // end namespace
///////////////////////////   ./include/algorithms/EigenSystemSolvers.h   ///////////////////////////


// Given the eigenvalues d[0..n-1] and (optionally) the eigenvectors v[0..n-1][0..n-1] as determined by Jacobi (÷11.1) or tqli (÷11.4), this routine sorts the eigenvalues into descending
// order and rearranges the columns of v correspondingly. The method is straight insertion.
static void eigsrt(MML::Vector<Real> &d, MML::Matrix<Real> *v = NULL)
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
    struct Jacobi
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
        Jacobi(Matrix<Real> &aa) : n(aa.RowNum()), a(aa), v(n, n), d(n), nrot(0),
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
                        sm += abs(a[ip][iq]);
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
                        g = 100.0 * abs(a[ip][iq]);
                        if (i > 4 && g <= EPS * abs(d[ip]) && g <= EPS * abs(d[iq]))
                            a[ip][iq] = 0.0;
                        else if (abs(a[ip][iq]) > tresh)
                        {
                            h = d[iq] - d[ip];
                            if (g <= EPS * abs(h))
                                t = (a[ip][iq]) / h;
                            else
                            {
                                theta = 0.5 * h / (a[ip][iq]);
                                t = 1.0 / (abs(theta) + sqrt(1.0 + theta * theta));
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
        inline void rot(Matrix<Real> &a, const Real s, const Real tau, const int i,
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

    public:
        // Computes all eigenvalues and eigenvectors of a real symmetric matrix a[0..n-1][0..n-1]
        // by reduction to tridiagonal form followed by QL iteration. On output, d[0..n-1] contains
        // the eigenvalues of a sorted into descending order, while z[0..n-1][0..n-1] is a matrix
        // whose columns contain the corresponding normalized eigenvectors. If yesvecs is input as
        // true (the default), then the eigenvectors are computed. If yesvecs is input as false, only
        // the eigenvalues are computed.
        SymmMatEigenSolver(Matrix<Real> &a, bool yesvec = true) : n(a.RowNum()), z(a), d(n), e(n), yesvecs(yesvec)
        {
            tred2();
            tqli();
            sort();
        }

        // Computes all eigenvalues and (optionally) eigenvectors of a real, symmetric, tridiagonal
        // matrix by QL iteration. On input, dd[0..n-1] contains the diagonal elements of the tridiagonal matrix. The vector ee[0..n-1] inputs the subdiagonal elements of the tridiagonal
        // matrix, with ee[0] arbitrary. Output is the same as the constructor above.
        SymmMatEigenSolver(Vector<Real> &dd, Vector<Real> &ee, bool yesvec = true) : n((int)dd.size()), d(dd), e(ee), z(n, n), yesvecs(yesvec)
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
                        scale += abs(z[i][k]);
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
                        dd = abs(d[m]) + abs(d[m + 1]);
                        if (abs(e[m]) <= EPS * dd)
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
            return (absa > absb ? absa * sqrt(1.0 + SQR(absb / absa)) : (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + SQR(absa / absb))));
        }
    };
    //////////////////////////////////////////////////////////////////////////////

    // Computes all eigenvalues and eigenvectors of a real nonsymmetric matrix by reduction to Hessenberg form followed by QR iteration.
    class Unsymmeig
    {
    public:
        int n;
        Matrix<Real> a, zz;
        Vector<Complex> wri;
        Vector<Real> scale;
        Vector<int> perm;
        bool yesvecs, hessen;

        bool isRealEigenvalue(int ind) const   { return wri[ind].imag() == 0.0; }
        Vector<Complex> getEigenvalues() const { return wri; }

        int getNumReal() const 
        {
            int res = 0;
            for (int i = 0; i < n; i++)
                if(wri[i].imag() == 0.0 )
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
                if(wri[i].imag() == 0.0 )
                    res[cnt++] = wri[i].real();
            return res;
        }
        Vector<Complex> getComplexEigenvalues() const
        {
            Vector<Complex> res(getNumComplex());
            int cnt = 0;
            for (int i = 0; i < n; i++)
                if(wri[i].imag() != 0.0 )
                    res[cnt++] = wri[i];
            return res;
        }

        Vector<Complex> getEigenvector(int ind) const
        {
            Vector<Complex> res(n);

            if( isRealEigenvalue(ind) )
                for (int i = 0; i < n; i++)
                    res[i] = zz(i,ind);
            else
            {
                // count how many real eigenvalues are ther ebefore this one
                int cnt = 0;
                for (int i = 0; i < ind; i++)
                    if( isRealEigenvalue(i) )
                        cnt++;
                if( cnt % 2 == 0)
                {
                    if( ind % 2 == 0)
                    {
                        for (int i = 0; i < n; i++)
                            res[i] = Complex(zz(i,ind) , zz(i,ind+1));
                    }
                    else
                    {
                        for (int i = 0; i < n; i++)
                            res[i] = Complex(zz(i,ind-1) , -zz(i,ind));
                    }
                }
                else {
                    if( ind % 2 == 0 )
                    {
                        for (int i = 0; i < n; i++)
                            res[i] = Complex(zz(i,ind-1) , -zz(i,ind));
                    }
                    else
                    {
                        for (int i = 0; i < n; i++)
                            res[i] = Complex(zz(i,ind) , zz(i,ind+1));
                    }
                }
            }

            return res;
        }
        Vector<Real> getRealPartEigenvector(int ind) const
        {
            Vector<Real> res(n);
            for (int i = 0; i < n; i++)
                res[i] = zz(i,ind);

            return res;
        }        

        // Computes all eigenvalues and (optionally) eigenvectors of a real nonsymmetric matrix a[0..n-1][0..n-1] by reduction to Hessenberg form followed by QR iteration. 
        // If yesvecs is input as true (the default), then the eigenvectors are computed. Otherwise, only the eigenvalues are computed. 
        // If hessen is input as false (the default), the matrix is first reduced to Hessenberg form. Otherwise it is assumed that the matrix is already in Hessenberg from. 
        // On output, wri[0..n-1] contains the eigenvalues of a sorted into descending order, while zz[0..n-1][0..n-1] is a matrix whose columns contain the corresponding
        // eigenvectors. 
        // For a complex eigenvalue, only the eigenvector corresponding to the eigenvalue with a positive imaginary part is stored, with the real part in zz[0..n-1][i] and the
        // imaginary part in h.zz[0..n-1][i+1]. The eigenvectors are not normalized.

        Unsymmeig(Matrix<Real> &aa, bool yesvec = true, bool hessenb = false) : n(aa.RowNum()), a(aa), zz(n, n), wri(n), scale(n, 1.0), perm(n),
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
                            c += abs(a[j][i]);
                            r += abs(a[i][j]);
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
                    if (abs(a[j][m - 1]) > abs(x))
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
                    anorm += abs(a[i][j]);
            nn = n - 1;
            t = 0.0;
            while (nn >= 0)
            {
                its = 0;
                do
                {
                    for (l = nn; l > 0; l--)
                    {
                        s = abs(a[l - 1][l - 1]) + abs(a[l][l]);
                        if (s == 0.0)
                            s = anorm;
                        if (abs(a[l][l - 1]) <= EPS * s)
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
                            z = sqrt(abs(q));
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
                                s = abs(a[nn][nn - 1]) + abs(a[nn - 1][nn - 2]);
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
                                s = abs(p) + abs(q) + abs(r);
                                p /= s;
                                q /= s;
                                r /= s;
                                if (m == l)
                                    break;
                                u = abs(a[m][m - 1]) * (abs(q) + abs(r));
                                v = abs(p) * (abs(a[m - 1][m - 1]) + abs(z) + abs(a[m + 1][m + 1]));
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
                                    if ((x = abs(p) + abs(q) + abs(r)) != 0.0)
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
                    anorm += abs(a[i][j]);
            nn = n - 1;
            t = 0.0;
            while (nn >= 0)
            {
                its = 0;
                do
                {
                    for (l = nn; l > 0; l--)
                    {
                        s = abs(a[l - 1][l - 1]) + abs(a[l][l]);
                        if (s == 0.0)
                            s = anorm;
                        if (abs(a[l][l - 1]) <= EPS * s)
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
                            z = sqrt(abs(q));
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
                                s = abs(x) + abs(z);
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
                                s = abs(a[nn][nn - 1]) + abs(a[nn - 1][nn - 2]);
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
                                s = abs(p) + abs(q) + abs(r);
                                p /= s;
                                q /= s;
                                r /= s;
                                if (m == l)
                                    break;
                                u = abs(a[m][m - 1]) * (abs(q) + abs(r));
                                v = abs(p) * (abs(a[m - 1][m - 1]) + abs(z) + abs(a[m + 1][m + 1]));
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
                                    if ((x = abs(p) + abs(q) + abs(r)) != 0.0)
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
                                    q = SQR(real(wri[i]) - p) + SQR(imag(wri[i]));
                                    t = (x * s - z * r) / q;
                                    a[i][nn] = t;
                                    if (abs(x) > abs(z))
                                        a[i + 1][nn] = (-r - w * t) / x;
                                    else
                                        a[i + 1][nn] = (-s - y * t) / z;
                                }
                                t = abs(a[i][nn]);
                                if (EPS * t * t > 1)
                                    for (j = i; j <= nn; j++)
                                        a[j][nn] /= t;
                            }
                        }
                    }
                    else if (q < 0.0)
                    {
                        m = na;
                        if (abs(a[nn][na]) > abs(a[na][nn]))
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
                                    vr = SQR(real(wri[i]) - p) + SQR(imag(wri[i])) - q * q;
                                    vi = 2.0 * q * (real(wri[i]) - p);
                                    if (vr == 0.0 && vi == 0.0)
                                        vr = EPS * anorm * (abs(w) + abs(q) + abs(x) + abs(y) + abs(z));
                                    Complex temp = Complex(x * r - z * ra + q * sa, x * s - z * sa - q * ra) /
                                                   Complex(vr, vi);
                                    a[i][na] = real(temp);
                                    a[i][nn] = imag(temp);
                                    if (abs(x) > abs(z) + abs(q))
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
        ODESystem   &_sys;
        int    n, neqn;
        bool   dense;

        double &x;
        double xold;
        Vector<Real> &y, &dydx;

        double atol,rtol;
        double EPS;

        double hdid;
        double hnext;
        Vector<Real> yout,yerr;

        StepperBase(ODESystem &sys, Vector<Real> &yy, Vector<Real> &dydxx, double &xx, const double atoll, const double rtoll, bool dens) 
            : _sys(sys), x(xx),y(yy),dydx(dydxx),atol(atoll),
              rtol(rtoll),dense(dens),n(sys.getDim()),neqn(n),yout(n),yerr(n) {}

        virtual Real dense_out(const int i, const Real x, const Real h) = 0;
    };

    // Dormand-Prince fifth-order Runge-Kutta step with monitoring of local truncation error to ensure
    // accuracy and adjust stepsize.
    struct StepperDopr5 : StepperBase {
        Vector<Real> k2,k3,k4,k5,k6;
        Vector<Real> rcont1,rcont2,rcont3,rcont4,rcont5;
        Vector<Real> dydxnew;
        
        // Input to the constructor are the dependent variable y[0..n-1] and its derivative dydx[0..n-1]
        // at the starting value of the independent variable x. Also input are the absolute and relative
        // tolerances, atol and rtol, and the boolean dense, which is true if dense output is required        
        StepperDopr5(ODESystem &sys, Vector<Real> &yy, Vector<Real> &dydx, Real &xx, const Real atoll, const Real rtoll, bool dens) 
                : StepperBase(sys, yy, dydx, xx, atoll, rtoll, dens), 
                  k2(n),k3(n),k4(n),k5(n),k6(n),
                  rcont1(n),rcont2(n),rcont3(n),rcont4(n),rcont5(n),
                  dydxnew(n) 
        {
	        EPS=std::numeric_limits<Real>::epsilon();
        }
        
        // Attempts a step with stepsize htry. On output, y and x are replaced by their new values, hdid
        // is the stepsize that was actually accomplished, and hnext is the estimated next stepsize        
        void step(const Real htry) {
            Real h=htry;
            for (;;) {
                dy(h);
                Real err=error();
                if (con.success(err,h)) break;
                if (abs(h) <= abs(x)*EPS)
                    throw("stepsize underflow in StepperDopr5");
            }
            if (dense)
                prepare_dense(h);
            dydx=dydxnew;
            y=yout;
            xold=x;
            x += (hdid=h);
            hnext=con.hnext;
        }  

        // Given values for n variables y[0..n-1] and their derivatives dydx[0..n-1] known at x, use the
        // fifth-order Dormand-Prince Runge-Kutta method to advance the solution over an interval h and
        // store the incremented variables in yout[0..n-1]. Also store an estimate of the local truncation
        // error in yerr using the embedded fourth-order method        
        void dy(const Real h) {
            static const Real c2=0.2,c3=0.3,c4=0.8,c5=8.0/9.0,a21=0.2,a31=3.0/40.0,
            a32=9.0/40.0,a41=44.0/45.0,a42=-56.0/15.0,a43=32.0/9.0,a51=19372.0/6561.0,
            a52=-25360.0/2187.0,a53=64448.0/6561.0,a54=-212.0/729.0,a61=9017.0/3168.0,
            a62=-355.0/33.0,a63=46732.0/5247.0,a64=49.0/176.0,a65=-5103.0/18656.0,
            a71=35.0/384.0,a73=500.0/1113.0,a74=125.0/192.0,a75=-2187.0/6784.0,
            a76=11.0/84.0,e1=71.0/57600.0,e3=-71.0/16695.0,e4=71.0/1920.0,
            e5=-17253.0/339200.0,e6=22.0/525.0,e7=-1.0/40.0;
            Vector<Real> ytemp(n);
            int i;
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*a21*dydx[i];
            _sys.derivs(x+c2*h,ytemp,k2);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(a31*dydx[i]+a32*k2[i]);
            _sys.derivs(x+c3*h,ytemp,k3);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(a41*dydx[i]+a42*k2[i]+a43*k3[i]);
            _sys.derivs(x+c4*h,ytemp,k4);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(a51*dydx[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]);
            _sys.derivs(x+c5*h,ytemp,k5);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(a61*dydx[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]);
            Real xph=x+h;
            _sys.derivs(xph,ytemp,k6);
            for (i=0;i<n;i++)
                yout[i]=y[i]+h*(a71*dydx[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);
            _sys.derivs(xph,yout,dydxnew);
            for (i=0;i<n;i++) {
                yerr[i]=h*(e1*dydx[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*dydxnew[i]);
            }
        }

        // Store coefficients of interpolating polynomial for dense output in rcont1...rcont5        
        void prepare_dense(const Real h) {
            Vector<Real> ytemp(n);
            static const Real d1=-12715105075.0/11282082432.0,
            d3=87487479700.0/32700410799.0, d4=-10690763975.0/1880347072.0,
            d5=701980252875.0/199316789632.0, d6=-1453857185.0/822651844.0,
            d7=69997945.0/29380423.0;
            for (int i=0;i<n;i++) {
                rcont1[i]=y[i];
                Real ydiff=yout[i]-y[i];
                rcont2[i]=ydiff;
                Real bspl=h*dydx[i]-ydiff;
                rcont3[i]=bspl;
                rcont4[i]=ydiff-h*dydxnew[i]-bspl;
                rcont5[i]=h*(d1*dydx[i]+d3*k3[i]+d4*k4[i]+d5*k5[i]+d6*k6[i]+
                    d7*dydxnew[i]);
            }
        }

        // Evaluate interpolating polynomial for y[i] at location x, where xold <= x <= xold + h.
        Real dense_out(const int i, const Real x, const Real h) {
            Real s=(x-xold)/h;
            Real s1=1.0-s;
            return rcont1[i]+s*(rcont2[i]+s1*(rcont3[i]+s*(rcont4[i]+s1*rcont5[i])));
        }

        // Use yerr to compute norm of scaled error estimate. A value less than one means the step was successful.        
        Real error() {
            Real err=0.0,sk;
            for (int i=0;i<n;i++) {
                sk=atol+rtol*std::max(abs(y[i]),abs(yout[i]));
                err += SQR(yerr[i]/sk);
            }
            return sqrt(err/n);
        }
        
        // Finally, the controller tests whether err <= 1 and adjusts the stepsize. The
        // default setting is beta D 0 (no PI control). Set beta to 0.04 or 0.08 to turn on PI control        
        struct Controller {
            Real hnext,errold;
            bool reject;

            Controller() : reject(false), errold(1.0e-4) {}

            // Returns true if err <= 1, false otherwise. If step was successful, sets hnext to the estimated
            // optimal stepsize for the next step. If the step failed, reduces h appropriately for another try.
            bool success(const Real err,Real &h) {
                static const Real beta=0.0,alpha=0.2-beta*0.75,safe=0.9,minscale=0.2, maxscale=10.0;
                Real scale;
                if (err <= 1.0) {
                    if (err == 0.0)
                        scale=maxscale;
                    else {
                        scale=safe*pow(err,-alpha)*pow(errold,beta);
                        if (scale<minscale) scale=minscale;
                        if (scale>maxscale) scale=maxscale;
                    }
                    if (reject)
                        hnext=h*std::min(scale,1.0);
                    else
                        hnext=h*scale;
                    errold=std::max(err,1.0e-4);
                    reject=false;
                    return true;
                } else {
                    scale=std::max(safe*pow(err,-alpha),minscale);
                    h *= scale;
                    reject=true;
                    return false;
                }
            }
        };
        Controller con;
    }; 

    struct Dopr853_constants {
        const Real c2  = 0.526001519587677318785587544488e-01;
        const Real c3  = 0.789002279381515978178381316732e-01;
        const Real c4  = 0.118350341907227396726757197510e+00;
        const Real c5  = 0.281649658092772603273242802490e+00;
        const Real c6  = 0.333333333333333333333333333333e+00;
        const Real c7  = 0.25e+00;
        const Real c8  = 0.307692307692307692307692307692e+00;
        const Real c9  = 0.651282051282051282051282051282e+00;
        const Real c10 = 0.6e+00;
        const Real c11 = 0.857142857142857142857142857142e+00;
        const Real c14 = 0.1e+00;
        const Real c15 = 0.2e+00;
        const Real c16 = 0.777777777777777777777777777778e+00;

        const Real b1 =   5.42937341165687622380535766363e-2;
        const Real b6 =   4.45031289275240888144113950566e0;
        const Real b7 =   1.89151789931450038304281599044e0;
        const Real b8 =  -5.8012039600105847814672114227e0;
        const Real b9 =   3.1116436695781989440891606237e-1;
        const Real b10 = -1.52160949662516078556178806805e-1;
        const Real b11 =  2.01365400804030348374776537501e-1;
        const Real b12 =  4.47106157277725905176885569043e-2;

        const Real bhh1 = 0.244094488188976377952755905512e+00;
        const Real bhh2 = 0.733846688281611857341361741547e+00;
        const Real bhh3 = 0.220588235294117647058823529412e-01;

        const Real er1  =  0.1312004499419488073250102996e-01;
        const Real er6  = -0.1225156446376204440720569753e+01;
        const Real er7  = -0.4957589496572501915214079952e+00;
        const Real er8  =  0.1664377182454986536961530415e+01;
        const Real er9  = -0.3503288487499736816886487290e+00;
        const Real er10 =  0.3341791187130174790297318841e+00;
        const Real er11 =  0.8192320648511571246570742613e-01;
        const Real er12 = -0.2235530786388629525884427845e-01;

        const Real a21 =    5.26001519587677318785587544488e-2;
        const Real a31 =    1.97250569845378994544595329183e-2;
        const Real a32 =    5.91751709536136983633785987549e-2;
        const Real a41 =    2.95875854768068491816892993775e-2;
        const Real a43 =    8.87627564304205475450678981324e-2;
        const Real a51 =    2.41365134159266685502369798665e-1;
        const Real a53 =   -8.84549479328286085344864962717e-1;
        const Real a54 =    9.24834003261792003115737966543e-1;
        const Real a61 =    3.7037037037037037037037037037e-2;
        const Real a64 =    1.70828608729473871279604482173e-1;
        const Real a65 =    1.25467687566822425016691814123e-1;
        const Real a71 =    3.7109375e-2;
        const Real a74 =    1.70252211019544039314978060272e-1;
        const Real a75 =    6.02165389804559606850219397283e-2;
        const Real a76 =   -1.7578125e-2;

        const Real a81 =    3.70920001185047927108779319836e-2;
        const Real a84 =    1.70383925712239993810214054705e-1;
        const Real a85 =    1.07262030446373284651809199168e-1;
        const Real a86 =   -1.53194377486244017527936158236e-2;
        const Real a87 =    8.27378916381402288758473766002e-3;
        const Real a91 =    6.24110958716075717114429577812e-1;
        const Real a94 =   -3.36089262944694129406857109825e0;
        const Real a95 =   -8.68219346841726006818189891453e-1;
        const Real a96 =    2.75920996994467083049415600797e1;
        const Real a97 =    2.01540675504778934086186788979e1;
        const Real a98 =   -4.34898841810699588477366255144e1;
        const Real a101 =   4.77662536438264365890433908527e-1;
        const Real a104 =  -2.48811461997166764192642586468e0;
        const Real a105 =  -5.90290826836842996371446475743e-1;
        const Real a106 =   2.12300514481811942347288949897e1;
        const Real a107 =   1.52792336328824235832596922938e1;
        const Real a108 =  -3.32882109689848629194453265587e1;
        const Real a109 =  -2.03312017085086261358222928593e-2;

        const Real a111 =  -9.3714243008598732571704021658e-1;
        const Real a114 =   5.18637242884406370830023853209e0;
        const Real a115 =   1.09143734899672957818500254654e0;
        const Real a116 =  -8.14978701074692612513997267357e0;
        const Real a117 =  -1.85200656599969598641566180701e1;
        const Real a118 =   2.27394870993505042818970056734e1;
        const Real a119 =   2.49360555267965238987089396762e0;
        const Real a1110 = -3.0467644718982195003823669022e0;
        const Real a121 =   2.27331014751653820792359768449e0;
        const Real a124 =  -1.05344954667372501984066689879e1;
        const Real a125 =  -2.00087205822486249909675718444e0;
        const Real a126 =  -1.79589318631187989172765950534e1;
        const Real a127 =   2.79488845294199600508499808837e1;
        const Real a128 =  -2.85899827713502369474065508674e0;
        const Real a129 =  -8.87285693353062954433549289258e0;
        const Real a1210 =  1.23605671757943030647266201528e1;
        const Real a1211 =  6.43392746015763530355970484046e-1;

        const Real a141 =  5.61675022830479523392909219681e-2;
        const Real a147 =  2.53500210216624811088794765333e-1;
        const Real a148 = -2.46239037470802489917441475441e-1;
        const Real a149 = -1.24191423263816360469010140626e-1;
        const Real a1410 =  1.5329179827876569731206322685e-1;
        const Real a1411 =  8.20105229563468988491666602057e-3;
        const Real a1412 =  7.56789766054569976138603589584e-3;
        const Real a1413 = -8.298e-3;

        const Real a151 =  3.18346481635021405060768473261e-2;
        const Real a156 =  2.83009096723667755288322961402e-2;
        const Real a157 =  5.35419883074385676223797384372e-2;
        const Real a158 = -5.49237485713909884646569340306e-2;
        const Real a1511 = -1.08347328697249322858509316994e-4;
        const Real a1512 =  3.82571090835658412954920192323e-4;
        const Real a1513 = -3.40465008687404560802977114492e-4;
        const Real a1514 =  1.41312443674632500278074618366e-1;
        const Real a161 = -4.28896301583791923408573538692e-1;
        const Real a166 = -4.69762141536116384314449447206e0;
        const Real a167 =  7.68342119606259904184240953878e0;
        const Real a168 =  4.06898981839711007970213554331e0;
        const Real a169 =  3.56727187455281109270669543021e-1;
        const Real a1613 = -1.39902416515901462129418009734e-3;
        const Real a1614 =  2.9475147891527723389556272149e0;
        const Real a1615 = -9.15095847217987001081870187138e0;

        const Real d41  = -0.84289382761090128651353491142e+01;
        const Real d46  =  0.56671495351937776962531783590e+00;
        const Real d47  = -0.30689499459498916912797304727e+01;
        const Real d48  =  0.23846676565120698287728149680e+01;
        const Real d49  =  0.21170345824450282767155149946e+01;
        const Real d410 = -0.87139158377797299206789907490e+00;
        const Real d411 =  0.22404374302607882758541771650e+01;
        const Real d412 =  0.63157877876946881815570249290e+00;
        const Real d413 = -0.88990336451333310820698117400e-01;
        const Real d414 =  0.18148505520854727256656404962e+02;
        const Real d415 = -0.91946323924783554000451984436e+01;
        const Real d416 = -0.44360363875948939664310572000e+01;

        const Real d51  =  0.10427508642579134603413151009e+02;
        const Real d56  =  0.24228349177525818288430175319e+03;
        const Real d57  =  0.16520045171727028198505394887e+03;
        const Real d58  = -0.37454675472269020279518312152e+03;
        const Real d59  = -0.22113666853125306036270938578e+02;
        const Real d510 =  0.77334326684722638389603898808e+01;
        const Real d511 = -0.30674084731089398182061213626e+02;
        const Real d512 = -0.93321305264302278729567221706e+01;
        const Real d513 =  0.15697238121770843886131091075e+02;
        const Real d514 = -0.31139403219565177677282850411e+02;
        const Real d515 = -0.93529243588444783865713862664e+01;
        const Real d516 =  0.35816841486394083752465898540e+02;

        const Real d61 =  0.19985053242002433820987653617e+02;
        const Real d66 = -0.38703730874935176555105901742e+03;
        const Real d67 = -0.18917813819516756882830838328e+03;
        const Real d68 =  0.52780815920542364900561016686e+03;
        const Real d69 = -0.11573902539959630126141871134e+02;
        const Real d610 =  0.68812326946963000169666922661e+01;
        const Real d611 = -0.10006050966910838403183860980e+01;
        const Real d612 =  0.77771377980534432092869265740e+00;
        const Real d613 = -0.27782057523535084065932004339e+01;
        const Real d614 = -0.60196695231264120758267380846e+02;
        const Real d615 =  0.84320405506677161018159903784e+02;
        const Real d616 =  0.11992291136182789328035130030e+02;

        const Real d71  = -0.25693933462703749003312586129e+02;
        const Real d76  = -0.15418974869023643374053993627e+03;
        const Real d77  = -0.23152937917604549567536039109e+03;
        const Real d78  =  0.35763911791061412378285349910e+03;
        const Real d79  =  0.93405324183624310003907691704e+02;
        const Real d710 = -0.37458323136451633156875139351e+02;
        const Real d711 =  0.10409964950896230045147246184e+03;
        const Real d712 =  0.29840293426660503123344363579e+02;
        const Real d713 = -0.43533456590011143754432175058e+02;
        const Real d714 =  0.96324553959188282948394950600e+02;
        const Real d715 = -0.39177261675615439165231486172e+02;
        const Real d716 = -0.14972683625798562581422125276e+03;
    };
    
    struct StepperDopr853 : StepperBase, Dopr853_constants {
        Vector<Real> yerr2;
        Vector<Real> k2,k3,k4,k5,k6,k7,k8,k9,k10;
        Vector<Real> rcont1,rcont2,rcont3,rcont4,rcont5,rcont6,rcont7,rcont8;
        
        StepperDopr853(ODESystem &sys, Vector<Real> &yy,Vector<Real> &dydxx,Real &xx,
                       const Real atoll,const Real rtoll,bool dens) 
                        : StepperBase(sys, yy,dydxx,xx,atoll,rtoll,dens),
                            yerr2(n),k2(n),k3(n),k4(n),
                            k5(n),k6(n),k7(n),k8(n),k9(n),k10(n),rcont1(n),rcont2(n),rcont3(n),
                            rcont4(n),rcont5(n),rcont6(n),rcont7(n),rcont8(n) 
        {
            EPS=std::numeric_limits<Real>::epsilon();
        }

        void step(const Real htry) {
            Vector<Real> dydxnew(n);
            Real h=htry;
            for (;;) {
                dy(h);
                Real err=error(h);
                if (con.success(err,h)) break;
                if (abs(h) <= abs(x)*EPS)
                    throw("stepsize underflow in StepperDopr853");
            }
            _sys.derivs(x+h,yout,dydxnew);
            if (dense)
                prepare_dense(h,dydxnew);
            dydx=dydxnew;
            y=yout;
            xold=x;
            x += (hdid=h);
            hnext=con.hnext;
        }

        void dy(const Real h) {
            Vector<Real> ytemp(n);
            int i;
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*a21*dydx[i];
            _sys.derivs(x+c2*h,ytemp,k2);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(a31*dydx[i]+a32*k2[i]);
            _sys.derivs(x+c3*h,ytemp,k3);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(a41*dydx[i]+a43*k3[i]);
            _sys.derivs(x+c4*h,ytemp,k4);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(a51*dydx[i]+a53*k3[i]+a54*k4[i]);
            _sys.derivs(x+c5*h,ytemp,k5);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(a61*dydx[i]+a64*k4[i]+a65*k5[i]);
            _sys.derivs(x+c6*h,ytemp,k6);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(a71*dydx[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);
            _sys.derivs(x+c7*h,ytemp,k7);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(a81*dydx[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i]);
            _sys.derivs(x+c8*h,ytemp,k8);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(a91*dydx[i]+a94*k4[i]+a95*k5[i]+a96*k6[i]+a97*k7[i]+
                    a98*k8[i]);
            _sys.derivs(x+c9*h,ytemp,k9);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(a101*dydx[i]+a104*k4[i]+a105*k5[i]+a106*k6[i]+
                    a107*k7[i]+a108*k8[i]+a109*k9[i]);
            _sys.derivs(x+c10*h,ytemp,k10);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(a111*dydx[i]+a114*k4[i]+a115*k5[i]+a116*k6[i]+
                    a117*k7[i]+a118*k8[i]+a119*k9[i]+a1110*k10[i]);
            _sys.derivs(x+c11*h,ytemp,k2);
            Real xph=x+h;
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(a121*dydx[i]+a124*k4[i]+a125*k5[i]+a126*k6[i]+
                    a127*k7[i]+a128*k8[i]+a129*k9[i]+a1210*k10[i]+a1211*k2[i]);
            _sys.derivs(xph,ytemp,k3);
            for (i=0;i<n;i++) {
                k4[i]=b1*dydx[i]+b6*k6[i]+b7*k7[i]+b8*k8[i]+b9*k9[i]+b10*k10[i]+
                    b11*k2[i]+b12*k3[i];
                yout[i]=y[i]+h*k4[i];
            }
            for (i=0;i<n;i++) {
                yerr[i]=k4[i]-bhh1*dydx[i]-bhh2*k9[i]-bhh3*k3[i];
                yerr2[i]=er1*dydx[i]+er6*k6[i]+er7*k7[i]+er8*k8[i]+er9*k9[i]+
                    er10*k10[i]+er11*k2[i]+er12*k3[i];
            }
        }
        void prepare_dense(const Real h,const Vector<Real> &dydxnew)
        {
            int i;
            Real ydiff,bspl;
            Vector<Real> ytemp(n);
            for (i=0;i<n;i++) {
                rcont1[i]=y[i];
                ydiff=yout[i]-y[i];
                rcont2[i]=ydiff;
                bspl=h*dydx[i]-ydiff;
                rcont3[i]=bspl;
                rcont4[i]=ydiff-h*dydxnew[i]-bspl;
                rcont5[i]=d41*dydx[i]+d46*k6[i]+d47*k7[i]+d48*k8[i]+
                    d49*k9[i]+d410*k10[i]+d411*k2[i]+d412*k3[i];
                rcont6[i]=d51*dydx[i]+d56*k6[i]+d57*k7[i]+d58*k8[i]+
                    d59*k9[i]+d510*k10[i]+d511*k2[i]+d512*k3[i];
                rcont7[i]=d61*dydx[i]+d66*k6[i]+d67*k7[i]+d68*k8[i]+
                    d69*k9[i]+d610*k10[i]+d611*k2[i]+d612*k3[i];
                rcont8[i]=d71*dydx[i]+d76*k6[i]+d77*k7[i]+d78*k8[i]+
                    d79*k9[i]+d710*k10[i]+d711*k2[i]+d712*k3[i];
            }
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(a141*dydx[i]+a147*k7[i]+a148*k8[i]+a149*k9[i]+
                    a1410*k10[i]+a1411*k2[i]+a1412*k3[i]+a1413*dydxnew[i]);
            _sys.derivs(x+c14*h,ytemp,k10);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(a151*dydx[i]+a156*k6[i]+a157*k7[i]+a158*k8[i]+
                    a1511*k2[i]+a1512*k3[i]+a1513*dydxnew[i]+a1514*k10[i]);
            _sys.derivs(x+c15*h,ytemp,k2);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(a161*dydx[i]+a166*k6[i]+a167*k7[i]+a168*k8[i]+
                    a169*k9[i]+a1613*dydxnew[i]+a1614*k10[i]+a1615*k2[i]);
            _sys.derivs(x+c16*h,ytemp,k3);
            for (i=0;i<n;i++)
            {
                rcont5[i]=h*(rcont5[i]+d413*dydxnew[i]+d414*k10[i]+d415*k2[i]+d416*k3[i]);
                rcont6[i]=h*(rcont6[i]+d513*dydxnew[i]+d514*k10[i]+d515*k2[i]+d516*k3[i]);
                rcont7[i]=h*(rcont7[i]+d613*dydxnew[i]+d614*k10[i]+d615*k2[i]+d616*k3[i]);
                rcont8[i]=h*(rcont8[i]+d713*dydxnew[i]+d714*k10[i]+d715*k2[i]+d716*k3[i]);
            }
        }
        Real dense_out(const int i,const Real x,const Real h) {
            Real s=(x-xold)/h;
            Real s1=1.0-s;
            return rcont1[i]+s*(rcont2[i]+s1*(rcont3[i]+s*(rcont4[i]+s1*(rcont5[i]+
                s*(rcont6[i]+s1*(rcont7[i]+s*rcont8[i]))))));
        }
        Real error(const Real h) {
            Real err=0.0,err2=0.0,sk,deno;
            for (int i=0;i<n;i++) {
                sk=atol+rtol*std::max(abs(y[i]),abs(yout[i]));
                err2 += SQR(yerr[i]/sk);
                err += SQR(yerr2[i]/sk);
            }
            deno=err+0.01*err2;
            if (deno <= 0.0)
                deno=1.0;
            return abs(h)*err*sqrt(1.0/(n*deno));
        }
        struct Controller {
            Real hnext,errold;
            bool reject;
            Controller() : reject(false), errold(1.0e-4) {}
            bool success(const Real err, Real &h) {
                static const Real beta=0.0,alpha=1.0/8.0-beta*0.2,safe=0.9,minscale=0.333,
                    maxscale=6.0;
                Real scale;
                if (err <= 1.0) {
                    if (err == 0.0)
                        scale=maxscale;
                    else {
                        scale=safe*pow(err,-alpha)*pow(errold,beta);
                        if (scale<minscale) scale=minscale;
                        if (scale>maxscale) scale=maxscale;
                    }
                    if (reject)
                        hnext=h*std::min(scale,1.0);
                    else
                        hnext=h*scale;
                    errold=std::max(err,1.0e-4);
                    reject=false;
                    return true;
                } else {
                    scale=std::max(safe*pow(err,-alpha),minscale);
                    h *= scale;
                    reject=true;
                    return false;
                }
            }
        };
        Controller con;
    };

    // Bulirsch-Stoer step with monitoring of local truncation error to ensure accuracy and adjust stepsize
    struct StepperBS : StepperBase {
        static const int KMAXX=8,IMAXX=KMAXX+1;     // KMAXX is the maximum number of rows used in the extrapolation
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

        StepperBS(ODESystem &sys, Vector<Real> &yy,Vector<Real> &dydxx,Real &xx, const Real atoll,const Real rtoll, bool dens) 
            : StepperBase(sys, yy,dydxx,xx,atoll,rtoll,dens),
                nseq(IMAXX),cost(IMAXX),
                table(KMAXX,n),dydxnew(n),coeff(IMAXX,IMAXX),errfac(2*IMAXX+2),ysave(IMAXX,n),
                fsave(IMAXX*(2*IMAXX+1),n),ipoint(IMAXX+1),dens((2*IMAXX+5)*n) 
        {
            EPS=std::numeric_limits<Real>::epsilon();
            if (dense)
                for (int i=0;i<IMAXX;i++)
                    nseq[i]=4*i+2;
            else
                for (int i=0;i<IMAXX;i++)
                    nseq[i]=2*(i+1);
            cost[0]=nseq[0]+1;
            for (int k=0;k<KMAXX;k++) cost[k+1]=cost[k]+nseq[k+1];
            hnext=-1.0e99;
            Real logfact=-log10(std::max(1.0e-12,rtol))*0.6+0.5;
            k_targ=std::max(1,std::min(KMAXX-1,int(logfact)));
            for (int k = 0; k<IMAXX; k++) {
                for (int l=0; l<k; l++) {
                    Real ratio=Real(nseq[k])/nseq[l];
                    coeff[k][l]=1.0/(ratio*ratio-1.0);
                }
            }
            for (int i=0; i<2*IMAXX+1; i++) {
                int ip5=i+5;
                errfac[i]=1.0/(ip5*ip5);
                Real e = 0.5*sqrt(Real(i+1)/ip5);
                for (int j=0; j<=i; j++) {
                    errfac[i] *= e/(j+1);
                }
            }
            ipoint[0]=0;
            for (int i=1; i<=IMAXX; i++) {
                int njadd=4*i-2;
                if (nseq[i-1] > njadd) njadd++;
                ipoint[i]=ipoint[i-1]+njadd;
            }
        }

        void step(const Real htry) {
            const Real STEPFAC1=0.65,STEPFAC2=0.94,STEPFAC3=0.02,STEPFAC4=4.0,
                KFAC1=0.8,KFAC2=0.9;
            static bool first_step=true,last_step=false;
            static bool forward,reject=false,prev_reject=false;
            int i,k;
            Real fac,h,hnew,hopt_int,err;
            bool firstk;
            Vector<Real> hopt(IMAXX),work(IMAXX);
            Vector<Real> ysav(n),yseq(n);
            Vector<Real> ymid(n),scale(n);
            work[0]=0;
            h=htry;
            forward = h>0 ? true : false;
            for (i=0;i<n;i++) ysav[i]=y[i];
            if (h != hnext && !first_step) {
                last_step=true;
            }
            if (reject) {
                prev_reject=true;
                last_step=false;
            }
            reject=false;
            firstk=true;
            hnew=abs(h);
            interp_error:
            while (firstk || reject) {
                h = forward ? hnew : -hnew;
                firstk=false;
                reject=false;
                if (abs(h) <= abs(x)*EPS)
                    throw("step size underflow in StepperBS");
                int ipt=-1;
                for (k=0; k<=k_targ+1;k++) {
                    dy(ysav,h,k,yseq,ipt);
                    if (k == 0)
                        y=yseq;
                    else
                        for (i=0;i<n;i++)
                            table[k-1][i]=yseq[i];
                    if (k != 0) {
                        polyextr(k,table,y);
                        err=0.0;
                        for (i=0;i<n;i++) {
                            scale[i]=atol+rtol*std::max(abs(ysav[i]),abs(y[i]));
                            err+=SQR((y[i]-table[0][i])/scale[i]);
                        }
                        err=sqrt(err/n);
                        Real expo=1.0/(2*k+1);
                        Real facmin=pow(STEPFAC3,expo);
                        if (err == 0.0)
                            fac=1.0/facmin;
                        else {
                            fac=STEPFAC2/pow(err/STEPFAC1,expo);
                            fac=std::max(facmin/STEPFAC4,std::min(1.0/facmin,fac));
                        }
                        hopt[k]=abs(h*fac);
                        work[k]=cost[k]/hopt[k];
                        if ((first_step || last_step) && err <= 1.0)
                            break;
                        if (k == k_targ-1 && !prev_reject && !first_step && !last_step) {
                            if (err <= 1.0)
                                break;
                            else if (err>SQR(nseq[k_targ]*nseq[k_targ+1]/(nseq[0]*nseq[0]))) {
                                reject=true;
                                k_targ=k;
                                if (k_targ>1 && work[k-1]<KFAC1*work[k])
                                    k_targ--;
                                hnew=hopt[k_targ];
                                break;
                            }
                        }
                        if (k == k_targ) {
                            if (err <= 1.0)
                                break;
                            else if (err>SQR(nseq[k+1]/nseq[0])) {
                                reject=true;
                                if (k_targ>1 && work[k-1]<KFAC1*work[k])
                                    k_targ--;
                                hnew=hopt[k_targ];
                                break;
                            }
                        }
                        if (k == k_targ+1) {
                            if (err > 1.0) {
                                reject=true;
                                if (k_targ>1 && work[k_targ-1]<KFAC1*work[k_targ])
                                    k_targ--;
                                hnew=hopt[k_targ];
                            }
                            break;
                        }
                    }
                }
                if (reject)
                    prev_reject=true;
            }
            _sys.derivs(x+h,y,dydxnew);
            if (dense) {
                prepare_dense(h,dydxnew,ysav,scale,k,err);
                hopt_int=h/std::max(pow(err,1.0/(2*k+3)),0.01);
                if (err > 10.0) {
                    hnew=abs(hopt_int);
                    reject=true;
                    prev_reject=true;
                    goto interp_error;
                }
            }
            dydx=dydxnew;
            xold=x;
            x+=h;
            hdid=h;
            first_step=false;
            int kopt;
            if (k == 1)
                kopt=2;
            else if (k <= k_targ) {
                kopt=k;
                if (work[k-1] < KFAC1*work[k])
                    kopt=k-1;
                else if (work[k] < KFAC2*work[k-1])
                    kopt=std::min(k+1,KMAXX-1);
            } else {
                kopt=k-1;
                if (k > 2 && work[k-2] < KFAC1*work[k-1])
                    kopt=k-2;
                if (work[k] < KFAC2*work[kopt])
                    kopt=std::min(k,KMAXX-1);
            }
            if (prev_reject) {
                k_targ=std::min(kopt,k);
                hnew=std::min(std::abs(h),hopt[k_targ]);
                prev_reject=false;
            }
            else {
                if (kopt <= k)
                    hnew=hopt[kopt];
                else {
                    if (k<k_targ && work[k]<KFAC2*work[k-1])
                        hnew=hopt[k]*cost[kopt+1]/cost[k];
                    else
                        hnew=hopt[k]*cost[kopt]/cost[k];
                }
                k_targ=kopt;
            }
            if (dense)
                hnew=std::min(hnew,std::abs(hopt_int));
            if (forward)
                hnext=hnew;
            else
                hnext=-hnew;	
        }

        virtual void dy(const Vector<Real> &y,const Real htot,const int k,Vector<Real> &yend, int &ipt) 
        {
            Vector<Real> ym(n),yn(n);
            int nstep=nseq[k];
            Real h=htot/nstep;
            for (int i=0;i<n;i++) {
                ym[i]=y[i];
                yn[i]=y[i]+h*dydx[i];
            }
            Real xnew=x+h;
            _sys.derivs(xnew,yn,yend);
            Real h2=2.0*h;
            for (int nn=1;nn<nstep;nn++) {
                if (dense && nn == nstep/2) {
                        for (int i=0;i<n;i++)
                            ysave[k][i]=yn[i];
                }
                if (dense && abs(nn-nstep/2) <= 2*k+1) {
                    ipt++;
                    for (int i=0;i<n;i++)
                        fsave[ipt][i]=yend[i];
                }
                for (int i=0;i<n;i++) {
                    Real swap=ym[i]+h2*yend[i];
                    ym[i]=yn[i];
                    yn[i]=swap;
                }
                xnew += h;
                _sys.derivs(xnew,yn,yend);
            }
            if (dense && nstep/2 <= 2*k+1) {
                ipt++;
                for (int i=0;i<n;i++)
                    fsave[ipt][i]=yend[i];
            }
            for (int i=0;i<n;i++)
                yend[i]=0.5*(ym[i]+yn[i]+h*yend[i]);
        }

        void polyextr(const int k, Matrix<Real> &table, Vector<Real> &last) 
        {
            int l=(int)last.size();
            for (int j=k-1; j>0; j--)
                for (int i=0; i<l; i++)
                    table[j-1][i]=table[j][i]+coeff[k][j]*(table[j][i]-table[j-1][i]);
            for (int i=0; i<l; i++)
                last[i]=table[0][i]+coeff[k][0]*(table[0][i]-last[i]);
        }

        virtual void prepare_dense(const Real h,const Vector<Real> &dydxnew, const Vector<Real> &ysav,const Vector<Real> &scale,const int k,Real &error) 
        {
            mu=2*k-1;
            for (int i=0; i<n; i++) {
                dens[i]=ysav[i];
                dens[n+i]=h*dydx[i];
                dens[2*n+i]=y[i];
                dens[3*n+i]=dydxnew[i]*h;
            }
            for (int j=1; j<=k; j++) {
                Real dblenj=nseq[j];
                for (int l=j; l>=1; l--) {
                    Real factor=SQR(dblenj/nseq[l-1])-1.0;
                    for (int i=0; i<n; i++)
                        ysave[l-1][i]=ysave[l][i]+(ysave[l][i]-ysave[l-1][i])/factor;
                }
            }
            for (int i=0; i<n; i++)
                dens[4*n+i]=ysave[0][i];
            for (int kmi=1; kmi<=mu; kmi++) {
                int kbeg=(kmi-1)/2;
                for (int kk=kbeg; kk<=k; kk++) {
                    Real facnj=pow(nseq[kk]/2.0,kmi-1);
                    int ipt=ipoint[kk+1]-2*kk+kmi-3;
                    for (int i=0; i<n; i++)
                        ysave[kk][i]=fsave[ipt][i]*facnj;
                }
                for (int j=kbeg+1; j<=k; j++) {
                    Real dblenj=nseq[j];
                    for (int l=j; l>=kbeg+1; l--) {
                        Real factor=SQR(dblenj/nseq[l-1])-1.0;
                        for (int i=0; i<n; i++)
                            ysave[l-1][i]=ysave[l][i]+
                                (ysave[l][i]-ysave[l-1][i])/factor;
                    }
                }
                for (int i=0; i<n; i++)
                    dens[(kmi+4)*n+i]=ysave[kbeg][i]*h;
                if (kmi == mu) continue;
                for (int kk=kmi/2; kk<=k; kk++) {
                    int lbeg=ipoint[kk+1]-1;
                    int lend=ipoint[kk]+kmi;
                    if (kmi == 1) lend += 2;
                    for (int l=lbeg; l>=lend; l-=2)
                        for (int i=0; i<n; i++)
                            fsave[l][i]=fsave[l][i]-fsave[l-2][i];
                    if (kmi == 1) {
                        int l=lend-2;
                        for (int i=0; i<n; i++)
                            fsave[l][i]=fsave[l][i]-dydx[i];
                    }
                }
                for (int kk=kmi/2; kk<=k; kk++) {
                    int lbeg=ipoint[kk+1]-2;
                    int lend=ipoint[kk]+kmi+1;
                    for (int l=lbeg; l>=lend; l-=2)
                        for (int i=0; i<n; i++)
                            fsave[l][i]=fsave[l][i]-fsave[l-2][i];
                }
            }
            dense_interp(n,dens,mu);
            error=0.0;
            if (mu >= 1) {
                for (int i=0; i<n; i++)
                    error += SQR(dens[(mu+4)*n+i]/scale[i]);
                error=sqrt(error/n)*errfac[mu-1];
            }
        }
        virtual Real dense_out(const int i,const Real x,const Real h) {
            Real theta=(x-xold)/h;
            Real theta1=1.0-theta;
            Real yinterp=dens[i]+theta*(dens[n+i]+theta1*(dens[2*n+i]*theta
                +dens[3*n+i]*theta1));
            if (mu<0)
                return yinterp;
            Real theta05=theta-0.5;
            Real t4=SQR(theta*theta1);
            Real c=dens[n*(mu+4)+i];
            for (int j=mu;j>0; j--)
                c=dens[n*(j+3)+i]+c*theta05/j;
            yinterp += t4*c;
            return yinterp;
        }

        virtual void dense_interp(const int n, Vector<Real> &y, const int imit) {
            Real y0,y1,yp0,yp1,ydiff,aspl,bspl,ph0,ph1,ph2,ph3,fac1,fac2;
            Vector<Real> a(31);
            for (int i=0; i<n; i++) {
                y0=y[i];
                y1=y[2*n+i];
                yp0=y[n+i];
                yp1=y[3*n+i];
                ydiff=y1-y0;
                aspl=-yp1+ydiff;
                bspl=yp0-ydiff;
                y[n+i]=ydiff;
                y[2*n+i]=aspl;
                y[3*n+i]=bspl;
                if (imit < 0) continue;
                ph0=(y0+y1)*0.5+0.125*(aspl+bspl);
                ph1=ydiff+(aspl-bspl)*0.25;
                ph2=-(yp0-yp1);
                ph3=6.0*(bspl-aspl);
                if (imit >= 1) {
                    a[1]=16.0*(y[5*n+i]-ph1);
                    if (imit >= 3) {
                        a[3]=16.0*(y[7*n+i]-ph3+3*a[1]);
                        for (int im=5; im <=imit; im+=2) {
                            fac1=im*(im-1)/2.0;
                            fac2=fac1*(im-2)*(im-3)*2.0;
                            a[im]=16.0*(y[(im+4)*n+i]+fac1*a[im-2]-fac2*a[im-4]);
                        }
                    }
                }
                a[0]=(y[4*n+i]-ph0)*16.0;
                if (imit >= 2) {
                    a[2]=(y[n*6+i]-ph2+a[0])*16.0;
                    for (int im=4; im <=imit; im+=2) {
                        fac1=im*(im-1)/2.0;
                        fac2=im*(im-1)*(im-2)*(im-3);
                        a[im]=(y[n*(im+4)+i]+a[im-2]*fac1-a[im-4]*fac2)*16.0;
                    }
                }
                for (int im=0; im<=imit; im++)
                    y[n*(im+4)+i]=a[im];
            }
        }
    };
}
///////////////////////////   ./include/algorithms/ODESystemSteppers_Stiff.h   ///////////////////////////





namespace MML
{
    struct Ross_constants {
        const Real c2=0.386;
        const Real c3=0.21;
        const Real c4=0.63;
        const Real bet2p=0.0317;
        const Real bet3p=0.0635;
        const Real bet4p=0.3438;
        const Real d1= 0.2500000000000000e+00;
        const Real d2=-0.1043000000000000e+00;
        const Real d3= 0.1035000000000000e+00;
        const Real d4=-0.3620000000000023e-01;
        const Real a21= 0.1544000000000000e+01;
        const Real a31= 0.9466785280815826e+00;
        const Real a32= 0.2557011698983284e+00;
        const Real a41= 0.3314825187068521e+01;
        const Real a42= 0.2896124015972201e+01;
        const Real a43= 0.9986419139977817e+00;
        const Real a51= 0.1221224509226641e+01;
        const Real a52= 0.6019134481288629e+01;
        const Real a53= 0.1253708332932087e+02;
        const Real a54=-0.6878860361058950e+00;
        const Real c21=-0.5668800000000000e+01;
        const Real c31=-0.2430093356833875e+01;
        const Real c32=-0.2063599157091915e+00;
        const Real c41=-0.1073529058151375e+00;
        const Real c42=-0.9594562251023355e+01;
        const Real c43=-0.2047028614809616e+02;
        const Real c51= 0.7496443313967647e+01;
        const Real c52=-0.1024680431464352e+02;
        const Real c53=-0.3399990352819905e+02;
        const Real c54= 0.1170890893206160e+02;
        const Real c61= 0.8083246795921522e+01;
        const Real c62=-0.7981132988064893e+01;
        const Real c63=-0.3152159432874371e+02;
        const Real c64= 0.1631930543123136e+02;
        const Real c65=-0.6058818238834054e+01;
        const Real gam= 0.2500000000000000e+00;
        const Real d21= 0.1012623508344586e+02;
        const Real d22=-0.7487995877610167e+01;
        const Real d23=-0.3480091861555747e+02;
        const Real d24=-0.7992771707568823e+01;
        const Real d25= 0.1025137723295662e+01;
        const Real d31=-0.6762803392801253e+00;
        const Real d32= 0.6087714651680015e+01;
        const Real d33= 0.1643084320892478e+02;
        const Real d34= 0.2476722511418386e+02;
        const Real d35=-0.6594389125716872e+01;
    };

    // Fourth-order stiffly stable Rosenbrock step for integrating stiff ODEs, with monitoring of local
    // truncation error to adjust stepsize.
    struct StepperRoss : StepperBase, Ross_constants {
        Matrix<Real> dfdy;
        Vector<Real> dfdx;
        Vector<Real> k1,k2,k3,k4,k5,k6;
        Vector<Real> cont1,cont2,cont3,cont4;
        Matrix<Real> a;

        // Input to the constructor are the dependent variable y[0..n-1] and its derivative dydx[0..n-1]
        // at the starting value of the independent variable x. Also input are the absolute and relative
        // tolerances, atol and rtol, and the boolean dense, which is true if dense output is required
        StepperRoss(ODESystem &sys, Vector<Real> &yy, Vector<Real> &dydxx, Real &xx, const Real atoll,const Real rtoll, bool dens) 
            : StepperBase(sys, yy,dydxx,xx,atoll,rtoll,dens),
              dfdy(n,n),dfdx(n),k1(n),k2(n),
              k3(n),k4(n),k5(n),k6(n),cont1(n),cont2(n),cont3(n),cont4(n),a(n,n) 
        {
            EPS=std::numeric_limits<Real>::epsilon();
        }
        
        // Attempts a step with stepsize htry. On output, y and x are replaced by their new values, hdid
        // is the stepsize that was actually accomplished, and hnext is the estimated next stepsize.        
        void step(const Real htry) {
            Vector<Real> dydxnew(n);
            Real h=htry;
            dynamic_cast<ODESystemWithJacobian&>(_sys).jacobian(x,y,dfdx,dfdy);
            for (;;) {
                dy(h);
                Real err=error();
                if (con.success(err,h)) break;
                if (abs(h) <= abs(x)*EPS)
                    throw("stepsize underflow in StepperRoss");
            }
            _sys.derivs(x+h,yout,dydxnew);
            if (dense)
                prepare_dense(h,dydxnew);
            dydx=dydxnew;
            y=yout;
            xold=x;
            x += (hdid=h);
            hnext=con.hnext;
        }

        // Given values for n variables y[0..n-1] and their derivatives dydx[0..n-1] known at x, use the
        // fourth-order stiffly stable Rosenbrock method to advance the solution over an interval h and
        // store the incremented variables in yout[0..n-1]. Also store an estimate of the local truncation
        // error in yerr using the embedded third-order method
        void dy(const Real h) {
            Vector<Real> ytemp(n),dydxnew(n);
            int i;
            for (i=0;i<n;i++) {
                for (int j=0;j<n;j++) a[i][j] = -dfdy[i][j];
                a[i][i] += 1.0/(gam*h);
            }
            LUDecompositionSolver alu(a);
            for (i=0;i<n;i++)
                ytemp[i]=dydx[i]+h*d1*dfdx[i];
            alu.Solve(ytemp,k1);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+a21*k1[i];
            _sys.derivs(x+c2*h,ytemp,dydxnew);
            for (i=0;i<n;i++)
                ytemp[i]=dydxnew[i]+h*d2*dfdx[i]+c21*k1[i]/h;
            alu.Solve(ytemp,k2);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+a31*k1[i]+a32*k2[i];
            _sys.derivs(x+c3*h,ytemp,dydxnew);
            for (i=0;i<n;i++)
                ytemp[i]=dydxnew[i]+h*d3*dfdx[i]+(c31*k1[i]+c32*k2[i])/h;
            alu.Solve(ytemp,k3);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+a41*k1[i]+a42*k2[i]+a43*k3[i];
            _sys.derivs(x+c4*h,ytemp,dydxnew);
            for (i=0;i<n;i++)
                ytemp[i]=dydxnew[i]+h*d4*dfdx[i]+(c41*k1[i]+c42*k2[i]+c43*k3[i])/h;
            alu.Solve(ytemp,k4);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+a51*k1[i]+a52*k2[i]+a53*k3[i]+a54*k4[i];
            Real xph=x+h;
            _sys.derivs(xph,ytemp,dydxnew);
            for (i=0;i<n;i++)
                k6[i]=dydxnew[i]+(c51*k1[i]+c52*k2[i]+c53*k3[i]+c54*k4[i])/h;
            alu.Solve(k6,k5);
            for (i=0;i<n;i++)
                ytemp[i] += k5[i];
            _sys.derivs(xph,ytemp,dydxnew);
            for (i=0;i<n;i++)
                k6[i]=dydxnew[i]+(c61*k1[i]+c62*k2[i]+c63*k3[i]+c64*k4[i]+c65*k5[i])/h;
            alu.Solve(k6,yerr);
            for (i=0;i<n;i++)
                yout[i]=ytemp[i]+yerr[i];
        }

        // Store coefficients of interpolating polynomial for dense output in cont1...cont4.
        void prepare_dense(const Real h,const Vector<Real> &dydxnew) {
            for (int i=0;i<n;i++) {
                cont1[i]=y[i];
                cont2[i]=yout[i];
                cont3[i]=d21*k1[i]+d22*k2[i]+d23*k3[i]+d24*k4[i]+d25*k5[i];
                cont4[i]=d31*k1[i]+d32*k2[i]+d33*k3[i]+d34*k4[i]+d35*k5[i];
            }
        }

        // Evaluate interpolating polynomial for y[i] at location x, where xold <= x <= xold + h        
        Real dense_out(const int i,const Real x,const Real h) {
            Real s=(x-xold)/h;
            Real s1=1.0-s;
            return cont1[i]*s1+s*(cont2[i]+s1*(cont3[i]+s*cont4[i]));
        }

        // Use yerr to compute norm of scaled error estimate. A value less than one means the step was successful        
        Real error() {
            Real err=0.0,sk;
            for (int i=0;i<n;i++) {
                sk=atol+rtol*std::max(abs(y[i]),abs(yout[i]));
                err += SQR(yerr[i]/sk);
            }
            return sqrt(err/n);
        }

        struct Controller {
            Real hnext;
            bool reject;
            bool first_step;
            Real errold;
            Real hold;
            Controller() : reject(false), first_step(true) {}
            bool success(Real err, Real &h) {
                static const Real safe=0.9,fac1=5.0,fac2=1.0/6.0;
                Real fac=std::max(fac2,std::min(fac1,pow(err,0.25)/safe));
                Real hnew=h/fac;
                if (err <= 1.0) {
                    if (!first_step) {
                        Real facpred=(hold/h)*pow(err*err/errold,0.25)/safe;
                        facpred=std::max(fac2,std::min(fac1,facpred));
                        fac=std::max(fac,facpred);
                        hnew=h/fac;
                    }
                    first_step=false;
                    hold=h;
                    errold=std::max(0.01,err);
                    if (reject)
                        hnew=(h >= 0.0 ? std::min(hnew,h) : std::max(hnew,h));
                    hnext=hnew;
                    reject=false;
                    return true;
                } else {
                    h=hnew;
                    reject=true;
                    return false;
                }
            }
        };
        Controller con;
    };

    struct StepperSemiImplExtr : StepperBase {
        static const int KMAXX=12,IMAXX=KMAXX+1;
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

        StepperSemiImplExtr(ODESystem &sys, Vector<Real> &yy, Vector<Real> &dydxx, Real &xx,
        const Real atoll,const Real rtoll, bool dens)
            : StepperBase(sys, yy,dydxx,xx,atoll,rtoll,dens),nseq(IMAXX),cost(IMAXX),
            table(KMAXX,n),dfdy(n,n),dfdx(n),calcjac(false),
            a(n,n),coeff(IMAXX,IMAXX),
            fsave((IMAXX-1)*(IMAXX+1)/2+2,n),dens((IMAXX+2)*n),factrl(IMAXX) 
        {
            static const Real costfunc=1.0,costjac=5.0,costlu=1.0,costsolve=1.0;
            EPS=std::numeric_limits<Real>::epsilon();
            jac_redo=std::min(1.0e-4,rtol);
            theta=2.0*jac_redo;
            nseq[0]=2;
            nseq[1]=3;
            for (int i=2;i<IMAXX;i++)
                nseq[i]=2*nseq[i-2];
            cost[0]=costjac+costlu+nseq[0]*(costfunc+costsolve);
            for (int k=0;k<KMAXX;k++)
                cost[k+1]=cost[k]+(nseq[k+1]-1)*(costfunc+costsolve)+costlu;
            hnext=-1.0e99;
            Real logfact=-log10(rtol+atol)*0.6+0.5;
            k_targ=std::max(1,std::min(KMAXX-1,int(logfact)));
            for (int k=0; k<IMAXX; k++) {
                for (int l=0; l<k; l++) {
                    Real ratio=Real(nseq[k])/nseq[l];
                    coeff[k][l]=1.0/(ratio-1.0);
                }
            }
            factrl[0]=1.0;
            for (int k=0; k<IMAXX-1; k++)
                factrl[k+1]=(k+1)*factrl[k];
        }
        void step(const Real htry) {
            const Real STEPFAC1=0.6,STEPFAC2=0.93,STEPFAC3=0.1,STEPFAC4=4.0,
                STEPFAC5=0.5,KFAC1=0.7,KFAC2=0.9;
            static bool first_step=true,last_step=false;
            static bool forward,reject=false,prev_reject=false;
            static Real errold;
            int i,k;
            Real fac,h,hnew,err;
            bool firstk;
            Vector<Real> hopt(IMAXX),work(IMAXX);
            Vector<Real> ysav(n),yseq(n);
            Vector<Real> ymid(n),scale(n);
            work[0]=1.e30;
            h=htry;
            forward = h>0 ? true : false;
            for (i=0;i<n;i++) ysav[i]=y[i];
            if (h != hnext && !first_step) {
                last_step=true;
            }
            if (reject) {
                prev_reject=true;
                last_step=false;
                theta=2.0*jac_redo;
            }
            for (i=0;i<n;i++)
                scale[i]=atol+rtol*abs(y[i]);
            reject=false;
            firstk=true;
            hnew=abs(h);
            compute_jac:
            if (theta > jac_redo && !calcjac) {
                dynamic_cast<ODESystemWithJacobian&>(_sys).jacobian(x,y,dfdx,dfdy);

                //_sys.jacobian(x,y,dfdx,dfdy);
                calcjac=true;
            }
            while (firstk || reject) {
                h = forward ? hnew : -hnew;
                firstk=false;
                reject=false;
                if (abs(h) <= abs(x)*EPS)
                    throw("step size underflow in StepperSemiImplExtr");
                int ipt=-1;
                for (k=0; k<=k_targ+1;k++) {
                    bool success=dy(ysav,h,k,yseq,ipt,scale);
                    if (!success) {
                        reject=true;
                        hnew=abs(h)*STEPFAC5;
                        break;
                    }
                    if (k == 0)
                        y=yseq;
                    else
                        for (i=0;i<n;i++)
                            table[k-1][i]=yseq[i];
                    if (k != 0) {
                        polyextr(k,table,y);
                        err=0.0;
                        for (i=0;i<n;i++) {
                            scale[i]=atol+rtol*abs(ysav[i]);
                            err+=SQR((y[i]-table[0][i])/scale[i]);
                        }
                        err=sqrt(err/n);
                        if (err > 1.0/EPS || (k > 1 && err >= errold)) {
                            reject=true;
                            hnew=abs(h)*STEPFAC5;
                            break;
                        }
                        errold=std::max(4.0*err,1.0);
                        Real expo=1.0/(k+1);
                        Real facmin=pow(STEPFAC3,expo);
                        if (err == 0.0)
                            fac=1.0/facmin;
                        else {
                            fac=STEPFAC2/pow(err/STEPFAC1,expo);
                            fac=std::max(facmin/STEPFAC4,std::min(1.0/facmin,fac));
                        }
                        hopt[k]=abs(h*fac);
                        work[k]=cost[k]/hopt[k];
                        if ((first_step || last_step) && err <= 1.0)
                            break;
                        if (k == k_targ-1 && !prev_reject && !first_step && !last_step) {
                            if (err <= 1.0)
                                break;
                            else if (err>nseq[k_targ]*nseq[k_targ+1]*4.0) {
                                reject=true;
                                k_targ=k;
                                if (k_targ>1 && work[k-1]<KFAC1*work[k])
                                    k_targ--;
                                hnew=hopt[k_targ];
                                break;
                            }
                        }
                        if (k == k_targ) {
                            if (err <= 1.0)
                                break;
                            else if (err>nseq[k+1]*2.0) {
                                reject=true;
                                if (k_targ>1 && work[k-1]<KFAC1*work[k])
                                    k_targ--;
                                hnew=hopt[k_targ];
                                break;
                            }
                        }
                        if (k == k_targ+1) {
                            if (err > 1.0) {
                                reject=true;
                                if (k_targ>1 && work[k_targ-1]<KFAC1*work[k_targ])
                                    k_targ--;
                                hnew=hopt[k_targ];
                            }
                            break;
                        }
                    }
                }
                if (reject) {
                    prev_reject=true;
                    if (!calcjac) {
                        theta=2.0*jac_redo;
                        goto compute_jac;
                    }
                }
            }
            calcjac=false;
            if (dense)
                prepare_dense(h,ysav,scale,k,err);
            xold=x;
            x+=h;
            hdid=h;
            first_step=false;
            int kopt;
            if (k == 1)
                kopt=2;
            else if (k <= k_targ) {
                kopt=k;
                if (work[k-1] < KFAC1*work[k])
                    kopt=k-1;
                else if (work[k] < KFAC2*work[k-1])
                    kopt=std::min(k+1,KMAXX-1);
            } else {
                kopt=k-1;
                if (k > 2 && work[k-2] < KFAC1*work[k-1])
                    kopt=k-2;
                if (work[k] < KFAC2*work[kopt])
                    kopt=std::min(k,KMAXX-1);
            }
            if (prev_reject) {
                k_targ=std::min(kopt,k);
                hnew=std::min(std::abs(h),hopt[k_targ]);
                prev_reject=false;
            }
            else {
                if (kopt <= k)
                    hnew=hopt[kopt];
                else {
                    if (k<k_targ && work[k]<KFAC2*work[k-1])
                        hnew=hopt[k]*cost[kopt+1]/cost[k];
                    else
                        hnew=hopt[k]*cost[kopt]/cost[k];
                }
                k_targ=kopt;
            }
            if (forward)
                hnext=hnew;
            else
                hnext=-hnew;	
        }

        bool dy(const Vector<Real> &y,const Real htot,const int k,Vector<Real> &yend, int &ipt,const Vector<Real> &scale) 
        {
            Vector<Real> del(n),ytemp(n),dytemp(n);
            int nstep=nseq[k];
            Real h=htot/nstep;
            for (int i=0;i<n;i++) {
                for (int j=0;j<n;j++) a[i][j] = -dfdy[i][j];
                a[i][i] += 1.0/h;
            }
            LUDecompositionSolver alu(a);
            Real xnew=x+h;
            _sys.derivs(xnew,y,del);
            for (int i=0;i<n;i++)
                ytemp[i]=y[i];
            alu.Solve(del,del);
            if (dense && nstep==k+1) {
                ipt++;
                for (int i=0;i<n;i++)
                    fsave[ipt][i]=del[i];
            }
            for (int nn=1;nn<nstep;nn++) {
                for (int i=0;i<n;i++)
                    ytemp[i] += del[i];
                xnew += h;
                _sys.derivs(xnew,ytemp,yend);
                if (nn ==1 && k<=1) {
                    Real del1=0.0;
                    for (int i=0;i<n;i++)
                        del1 += SQR(del[i]/scale[i]);
                    del1=sqrt(del1);
                    _sys.derivs(x+h,ytemp,dytemp);
                    for (int i=0;i<n;i++)
                        del[i]=dytemp[i]-del[i]/h;
                    alu.Solve(del,del);
                    Real del2=0.0;
                    for (int i=0;i<n;i++)
                        del2 += SQR(del[i]/scale[i]);
                    del2=sqrt(del2);
                    theta=del2/std::max(1.0,del1);
                    if (theta > 1.0)
                        return false;
                }
                alu.Solve(yend,del);
                if (dense && nn >= nstep-k-1) {
                    ipt++;
                    for (int i=0;i<n;i++)
                        fsave[ipt][i]=del[i];
                }
            }
            for (int i=0;i<n;i++)
                yend[i]=ytemp[i]+del[i];
            return true;
        }
        void polyextr(const int k,Matrix<Real> &table,Vector<Real> &last) {
            int l=(int) last.size();
            for (int j=k-1; j>0; j--)
                for (int i=0; i<l; i++)
                    table[j-1][i]=table[j][i]+coeff[k][j]*(table[j][i]-table[j-1][i]);
            for (int i=0; i<l; i++)
                last[i]=table[0][i]+coeff[k][0]*(table[0][i]-last[i]);
        }

        void prepare_dense(const Real h,const Vector<Real> &ysav,const Vector<Real> &scale,
                            const int k,Real &error) {
            kright=k;
            for (int i=0; i<n; i++) {
                dens[i]=ysav[i];
                dens[n+i]=y[i];
            }
            for (int klr=0; klr < kright; klr++) {
                if (klr >= 1) {
                    for (int kk=klr; kk<=k; kk++) {
                        int lbeg=((kk+3)*kk)/2;
                        int lend=lbeg-kk+1;
                        for (int l=lbeg; l>=lend; l--)
                            for (int i=0; i<n; i++)
                                fsave[l][i]=fsave[l][i]-fsave[l-1][i];
                    }
                }
                for (int kk=klr; kk<=k; kk++) {
                    Real facnj=nseq[kk];
                    facnj=pow(facnj,klr+1)/factrl[klr+1];
                    int ipt=((kk+3)*kk)/2;
                        int krn=(kk+2)*n;
                    for (int i=0; i<n; i++) {
                        dens[krn+i]=fsave[ipt][i]*facnj;
                    }
                }
                for (int j=klr+1; j<=k; j++) {
                    Real dblenj=nseq[j];
                    for (int l=j; l>=klr+1; l--) {
                        Real factor=dblenj/nseq[l-1]-1.0;
                        for (int i=0; i<n; i++) {
                            int krn=(l+2)*n+i;
                            dens[krn-n]=dens[krn]+(dens[krn]-dens[krn-n])/factor;
                        }
                    }
                }
            }
            for (int in=0; in<n; in++) {
                for (int j=1; j<=kright+1; j++) {
                    int ii=n*j+in;
                    dens[ii]=dens[ii]-dens[ii-n];
                }
            }
        }
        Real dense_out(const int i,const Real x,const Real h) {
            Real theta=(x-xold)/h;
            int k=kright;
            Real yinterp=dens[(k+1)*n+i];
            for (int j=1; j<=k; j++)
                yinterp=dens[(k+1-j)*n+i]+yinterp*(theta-1.0);
            return dens[i]+yinterp*theta;
        }
        void dense_interp(const int n, Vector<Real> &y, const int imit);
    };    
}
///////////////////////////   ./include/algorithms/ODESystemSolvers.h   ///////////////////////////



namespace MML
{

    struct Output {
        int kmax;
        int _nvar;
        int nsave;

        int nok;
        int nbad;
        int nstp;        

        bool dense;
        int count;

        Real x1,x2,xout,dxout;

        Vector<Real> xsave;
        Matrix<Real> ysave;
        
        Output() : kmax(-1),dense(false),count(0) {}

        // Constructor provides dense output at nsave equally spaced intervals. If nsave <= 0, output
        // is saved only at the actual integration steps.        
        Output(const int nsavee) : kmax(50),nsave(nsavee),count(0),xsave(kmax), nok(0), nbad(0) {
            dense = nsave > 0 ? true : false;
        }
        void init(const int neqn, const Real xlo, const Real xhi) {
            _nvar=neqn;
            if (kmax == -1) return;
            ysave.Resize(_nvar,kmax);
            if (dense) {
                x1=xlo;
                x2=xhi;
                xout=x1;
                dxout=(x2-x1)/nsave;
            }
        }
        void resize() {
            int kold=kmax;
            kmax *= 2;
            
            Vector<Real> tempvec(xsave);
            xsave.Resize(kmax);
            for (int k=0; k<kold; k++)
                xsave[k]=tempvec[k];
            
            Matrix<Real> tempmat(ysave);
            ysave.Resize(_nvar,kmax);
            for (int i=0; i<_nvar; i++)
                for (int k=0; k<kold; k++)
                    ysave[i][k]=tempmat[i][k];
        }

        // Invokes dense_out function of stepper routine to produce output at xout. Normally called
        // by out rather than directly. Assumes that xout is between xold and xold+h, where the
        // stepper must keep track of xold, the location of the previous step, and x=xold+h, the
        // current step        
        void save_dense(StepperBase &s, const Real xout, const Real h) {
            if (count == kmax) resize();
            for (int i=0;i<_nvar;i++)
                ysave[i][count]=s.dense_out(i,xout,h);
            xsave[count++]=xout;
        }

        // Saves values of current x and y.
        void save(const Real x, const Vector<Real> &y) {
            if (kmax <= 0) return;
            if (count == kmax) resize();
            for (int i=0;i<_nvar;i++)
                ysave[i][count]=y[i];
            xsave[count++]=x;
        }

        // Typically called by Odeint to produce dense output. Input variables are nstp, the current
        // step number, the current values of x and y, the stepper s, and the stepsize h. A call with
        // nstp=-1 saves the initial values. The routine checks whether x is greater than the desired
        // output point xout. If so, it calls save_dense.        
        void out(const int nstp,const Real x,const Vector<Real> &y,StepperBase &s,const Real h) {
            if (!dense)
                throw("dense output not set in Output!");
            if (nstp == -1) {
                save(x,y);
                xout += dxout;
            } else {
                while ((x-xout)*(x2-x1) > 0.0) {
                    save_dense(s,xout,h);
                    xout += dxout;
                }
            }
        }
    };

    template<class Stepper>
    class ODESystemSolver {
        static inline const double EPS=std::numeric_limits<Real>::epsilon();
        static const int MAXSTP=50000;
        
        Real         _curr_x;               // used as reference by stepper!
        Vector<Real> _curr_y;
        Vector<Real> _curr_dydx;
    public:
        ODESystem &_sys;
        Output    &_out;
        Stepper   _stepper;
        bool _dense;

        int getDim() { return _sys.getDim(); }

        ODESystemSolver(ODESystem &sys, const Real atol, const Real rtol, Output &out) 
                            : _sys(sys), _curr_y(sys.getDim()), _curr_dydx(sys.getDim()), _dense(out.dense), _out(out),
                              _stepper(sys, _curr_y, _curr_dydx, _curr_x, atol, rtol, _dense) 
        { }

        ODESystemSolution integrate(Vector<Real> &in_ystart, const Real x1, Real x2, Real h1, const Real hmin)
        {
            _out.init(_stepper.neqn,x1,x2);

            Real h   = SIGN(h1,x2-x1);    
            int  dim = _sys.getDim();

            _curr_x = x1;
            _curr_y = in_ystart;
            Vector<Real> &ystart = in_ystart;
            ODESystemSolution sol(dim, _out.nsave);

            _sys.derivs(_curr_x, _curr_y, _curr_dydx);

            if (_dense)
                _out.out(-1, _curr_x, _curr_y, _stepper, h);
            else
                _out.save(_curr_x, _curr_y);

            for (_out.nstp=0; _out.nstp<MAXSTP; _out.nstp++) 
            {
                if ((_curr_x+h*1.0001-x2)*(x2-x1) > 0.0)
                    h=x2-_curr_x;

                _stepper.step(h);

                if (_stepper.hdid == h) ++_out.nok; 
                else ++_out.nbad;

                if (_dense)
                    _out.out(_out.nstp, _curr_x, _curr_y, _stepper, _stepper.hdid);
                else
                    _out.save(_curr_x, _curr_y);

                if ((_curr_x-x2)*(x2-x1) >= 0.0) {
                    for (int i=0;i<getDim();i++) ystart[i]=_curr_y[i];

                    if (_out.kmax > 0 && abs(_out.xsave[_out.count-1]-x2) > 100.0*abs(x2)*EPS)
                        _out.save(_curr_x, _curr_y);
                    
                    for(int i=0; i<=_out.nsave; i++)
                    {
                        sol.xval[i] = _out.xsave[i];
                        for(int j=0; j<ystart.size(); j++)
                        {
                            sol.yval[j][i] = _out.ysave[j][i];
                        }
                    }

                    return sol;
                }
                if (abs(_stepper.hnext) <= hmin) 
                    throw("Step size too small in Odeint");
                
                h=_stepper.hnext;
            }
            throw("Too many steps in routine Odeint");
        }        
    };

}
///////////////////////////   ./include/algorithms/ODESystemSolversLegacy.h   ///////////////////////////




namespace MML
{
    class RungeKuttaSolverDumb
    {
    public:        
        void rk4(Vector<Real> &y, Vector<Real> &dydx, const double x, const double h,
            Vector<Real> &yout, ODESystem &sys)
        {
            int i;
            double xh,hh,h6;

            int n= (int) y.size();
            Vector<Real> dym(n),dyt(n),yt(n);
            hh=h*0.5;
            h6=h/6.0;
            xh=x+hh;
            for (i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
            sys.derivs(xh,yt,dyt);
            for (i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
            sys.derivs(xh,yt,dym);
            for (i=0;i<n;i++) {
                yt[i]=y[i]+h*dym[i];
                dym[i] += dyt[i];
            }
            sys.derivs(x+h,yt,dyt);
            for (i=0;i<n;i++)
                yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
        }

        ODESystemSolutionEqualSpacing integrate(ODESystem &sys, const Vector<Real> &vstart, const double x1, const double x2, int numSteps)
        {
            int i,k;
            double x,h;
            int dim=sys.getDim();
            
            ODESystemSolutionEqualSpacing sol(dim, numSteps);

            Vector<Real> v(vstart),vout(dim),dv(dim);
            for (i=0;i<dim;i++) {
                sol.yval[i][0]=v[i];
            }
            sol.xval[0]=x1;
            x=x1;
            h=(x2-x1)/numSteps;
            for (k=0;k<numSteps;k++) {
                sys.derivs(x,v,dv);
                rk4(v,dv,x,h,vout,sys);
                if (x+h == x)
                    throw("Step size too small in routine rkdumb");
                x += h;
                sol.xval[k+1]=x;
                for (i=0;i<dim;i++) {
                    v[i]=vout[i];
                    sol.yval[i][k+1]=v[i];
                }
            }

            return sol;
        }    
    };    

    class RungeKuttaNR2
    {
        void rkck(Vector<Real> &y, Vector<Real> &dydx, const double x,
            const double h, Vector<Real> &yout, Vector<Real> &yerr,
            ODESystem &sys)
        {
            static const double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875,
                b21=0.2, b31=3.0/40.0, b32=9.0/40.0, b41=0.3, b42 = -0.9,
                b43=1.2, b51 = -11.0/54.0, b52=2.5, b53 = -70.0/27.0,
                b54=35.0/27.0, b61=1631.0/55296.0, b62=175.0/512.0,
                b63=575.0/13824.0, b64=44275.0/110592.0, b65=253.0/4096.0,
                c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0,
                dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0,
                dc4=c4-13525.0/55296.0, dc5 = -277.00/14336.0, dc6=c6-0.25;
            int i;

            int n= (int) y.size();
            Vector<Real> ak2(n),ak3(n),ak4(n),ak5(n),ak6(n),ytemp(n);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+b21*h*dydx[i];
            sys.derivs(x+a2*h,ytemp,ak2);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
            sys.derivs(x+a3*h,ytemp,ak3);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
            sys.derivs(x+a4*h,ytemp,ak4);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
            sys.derivs(x+a5*h,ytemp,ak5);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
            sys.derivs(x+a6*h,ytemp,ak6);
            for (i=0;i<n;i++)
                yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
            for (i=0;i<n;i++)
                yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
        }

        void rkqs(Vector<Real> &y, Vector<Real> &dydx, double &x, const double htry,
            const double eps, Vector<Real> &yscal, double &hdid, double &hnext,
            ODESystem &sys)
        {
            const double SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4;
            int i;
            double errmax,h,htemp,xnew;

            int n= (int) y.size();
            h=htry;
            Vector<Real> yerr(n),ytemp(n);
            for (;;) {
                rkck(y,dydx,x,h,ytemp,yerr,sys);
                errmax=0.0;
                for (i=0;i<n;i++) errmax=std::max(errmax,fabs(yerr[i]/yscal[i]));
                errmax /= eps;
                if (errmax <= 1.0) break;
                htemp=SAFETY*h*pow(errmax,PSHRNK);
                h=(h >= 0.0 ? std::max(htemp,0.1*h) : std::min(htemp,0.1*h));
                xnew=x+h;
                if (xnew == x) throw("stepsize underflow in rkqs");
            }
            if (errmax > ERRCON) hnext=SAFETY*h*pow(errmax,PGROW);
            else hnext=5.0*h;
            x += (hdid=h);
            for (i=0;i<n;i++) y[i]=ytemp[i];
        }    

    public:

        ODESystemSolution integrate(ODESystem &sys, const Vector<Real> &ystart, const double x1, const double x2, int maxSteps, double minSaveInterval,
                                    const double eps, const double h1, const double hmin, int &nok, int &nbad)
        {
            const int MAXSTP=10000;
            const double TINY=1.0e-30;

            int dim = sys.getDim();
            ODESystemSolution sol(dim,maxSteps);

            int kount = 0;
            int i,nstp;
            double xsav,x,hnext,hdid,h;
            Vector<Real> yscal(dim),y(ystart),dydx(dim);

            x=x1;
            h=SIGN(h1,x2-x1);
            nok = nbad = kount = 0;

            if (maxSteps > 0) xsav=x-minSaveInterval*2.0;
            for (nstp=0;nstp<MAXSTP;nstp++) {
                sys.derivs(x,y,dydx);
                for (i=0;i<dim;i++)
                    yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
                if (maxSteps > 0 && kount < maxSteps-1 && fabs(x-xsav) > fabs(minSaveInterval)) {
                    for (i=0;i<dim;i++) sol.yval[i][kount]=y[i];
                    sol.xval[kount++]=x;
                    xsav=x;
                }
                if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
                rkqs(y,dydx,x,h,eps,yscal,hdid,hnext,sys);
                if (hdid == h) ++nok; else ++nbad;
                if ((x-x2)*(x2-x1) >= 0.0) {
                    if (maxSteps != 0) {
                        for (i=0;i<dim;i++) sol.yval[i][kount]=y[i];
                        sol.xval[kount++]=x;
                    }
                    return sol;
                }
                if (fabs(hnext) <= hmin) throw("Step size too small in odeint");
                h=hnext;
            }
            throw("Too many steps in routine odeint");
        }

    };
}

///////////////////////////   ./include/algorithms/DiffGeometryAlgorithms.h   ///////////////////////////




namespace MML
{  
    class DiffGeometry
    {
    public:
        template<int N>
        class CurveTangentUnit : public IParametricCurve<N>
        {
            const IParametricCurve<N> &_curve;
        public:
            CurveTangentUnit(const IParametricCurve<N> &curve) : _curve(curve) {}

            VectorN<Real, N> operator()(double t) const 
            {
                auto tangent_vec = Derivation::DeriveCurve<N>(_curve, t, nullptr);
                return tangent_vec / tangent_vec.NormL2();
            }
        };

        template<int N>
        static VectorN<Real, N> getTangent(const IParametricCurve<N> &curve, double t)
        {
            return Derivation::DeriveCurve<N>(curve, t, nullptr);
        }
        template<int N>
        static VectorN<Real, N> getTangentUnit(const IParametricCurve<N> &curve, double t)
        {
            auto tangent = getTangent(curve, t);
            return tangent / tangent.NormL2();
        }

        template<int N>
        static VectorN<Real, N> getNormal(const IParametricCurve<N> &curve, double t)
        {
            return Derivation::DeriveCurveSec<N>(curve, t, nullptr);
        }
        template<int N>
        static VectorN<Real, N> getNormalScaled(const IParametricCurve<N> &curve, double t)
        {
            CurveTangentUnit  helper(curve);

            return Derivation::DeriveCurve<N>(helper, t, nullptr);
        }        
        template<int N>
        static VectorN<Real, N> getNormalUnit(const IParametricCurve<N> &curve, double t)
        {
            auto normal = getNormal(curve, t);
            return normal / normal.NormL2();
        }
        template<int N>
        static VectorN<Real, N> getPrincipalNormal(const IParametricCurve<N> &curve, double t)
        {
            auto y_der_1 = Vector3Cartesian(getTangent(curve, t));
            auto y_der_2 = Vector3Cartesian(Derivation::DeriveCurveSec<3>(curve, t, nullptr));

            Vector3Cartesian vec_prod1 = VectorProd(y_der_2, y_der_1);
            Vector3Cartesian res_vec   = VectorProd(y_der_1, vec_prod1);

            return res_vec / (y_der_1.NormL2() * vec_prod1.NormL2());
        }        

        template<int N>
        static VectorN<Real, N> getBinormal(const IParametricCurve<N> &curve, double t)
        {
            auto y_der_1 = Vector3Cartesian(getTangent(curve, t));
            auto y_der_2 = Vector3Cartesian(Derivation::DeriveCurveSec<3>(curve, t, nullptr));

            Vector3Cartesian vec_prod1 = VectorProd(y_der_2, y_der_1);
            return vec_prod1 / vec_prod1.NormL2();
        }  

        template<int N>
        static VectorN<Real, N> getCurvatureVector(const IParametricCurve<N> &curve, double t)
        {
            auto y_der_1 = getTangent(curve, t);
            auto y_der_2 = Derivation::DeriveCurveSec<3>(curve, t, nullptr);

            double res1 = pow(y_der_1.NormL2(), -2.0);
            auto   vec2 = y_der_2 - res1 * y_der_1.ScalarProductCartesian(y_der_2) * y_der_1;

            return vec2 / res1;
        }  

        template<int N>
        static Real getCurvature(const IParametricCurve<N> &curve, double t)
        {
            auto y_der_1 = getTangent(curve, t);
            auto y_der_2 = Derivation::DeriveCurveSec<3>(curve, t, nullptr);

            double res1 = pow(y_der_1.NormL2(), -2.0);
            auto   vec2 = y_der_2 - res1 * y_der_1.ScalarProductCartesian(y_der_2) * y_der_1;
            double res2 = vec2.NormL2();

            return res1 * res2;
        }  
    
        static Real getCurvature3(const IParametricCurve<3> &curve, double t)
        {
            auto curve_first_der = Vector3Cartesian( getTangent(curve, t) );
            auto curve_sec_der   = Vector3Cartesian( Derivation::DeriveCurveSec<3>(curve, t, nullptr) );

            auto prod = VectorProd(curve_first_der, curve_sec_der);

            return prod.NormL2() / pow(curve_first_der.NormL2(), 3);
        }  

        static Real getTorsion3(const IParametricCurve<3> &curve, double t)
        {
            auto curve_first_der = Vector3Cartesian( getTangent(curve, t) );
            auto curve_sec_der   = Vector3Cartesian( Derivation::DeriveCurveSec<3>(curve, t, nullptr) );
            auto curve_third_der = Vector3Cartesian( Derivation::DeriveCurveThird<3>(curve, t, nullptr) );

            auto prod = VectorProd(curve_first_der, curve_sec_der);

            Real temp = prod.ScalarProductCartesian(curve_third_der);

            return -temp / pow(prod.NormL2(), 2);
        }  

        static Plane3D getOsculationPlane(const IParametricCurve<3> &curve, double t)
        {
            Vector3Cartesian vec_pnt( curve(t) );
            
            Plane3D ret( vec_pnt.getAsPoint(), Vector3Cartesian(getNormal(curve, t)) );

            return ret;
        }

        static Plane3D getNormalPlane(const IParametricCurve<3> &curve, double t)
        {
            Vector3Cartesian vec_pnt( curve(t) );
            
            Plane3D ret( vec_pnt.getAsPoint(), Vector3Cartesian(getTangentUnit(curve, t)) );

            return ret;
        }

        static Plane3D getRectifyingPlane(const IParametricCurve<3> &curve, double t)
        {
            Vector3Cartesian vec_pnt( curve(t) );
            
            Plane3D ret( vec_pnt.getAsPoint(), Vector3Cartesian(getBinormal(curve, t)) );

            return ret;
        }

        static void getMovingTrihedron(const IParametricCurve<3> &curve, double t, Vector3Cartesian &tangent, Vector3Cartesian &normal, Vector3Cartesian &binormal)
        {
            tangent   = Vector3Cartesian(getTangentUnit(curve, t));
            normal    = Vector3Cartesian(getPrincipalNormal(curve, t));
            binormal  = Vector3Cartesian(getBinormal(curve, t));
        }

        static bool isArcLengthParametrized(const IParametricCurve<3> &curve, double t1, double t2)
        {
            int numPnt = 100;
            double delta = (t2 - t1) / numPnt;
            for(double t=t1+delta; t < t2; t += delta)
            {
                double len = PathIntegration::ParametricCurveLength(curve, t1, t);
                if( fabs(len - (t - t1)) > 1e-03 )
                    return false;
            }

            return true;
        }        
    };
}

///////////////////////////   ./include/algorithms/FunctionAnalyzer.h   ///////////////////////////



namespace MML
{
    // TODO - function point analyzer at point
    class FunctionAnalyzer
    {
    public:
        enum IntegrationMethod { TRAP, SIMPSON, ROMBERG } ;

        static double FuncDiff(IRealFunction &f1, IRealFunction &f2, double a, double b, IntegrationMethod method = TRAP)
        {
            RealFuncDiffHelper helper(f1, f2);

            switch (method)
            {
                case SIMPSON:
                    return Integration::IntegrateSimpson(helper, a, b);
                case ROMBERG:
                    return Integration::IntegrateRomberg(helper, a, b);
                default:
                    return Integration::IntegrateTrap(helper, a, b);
            }
        }

        static double FuncDiffAbs(IRealFunction &f1, IRealFunction &f2, double a, double b, IntegrationMethod method = TRAP)
        {
            RealFuncDiffAbsHelper helper(f1, f2);

            switch (method)
            {
                case SIMPSON:
                    return Integration::IntegrateSimpson(helper, a, b);
                case ROMBERG:
                    return Integration::IntegrateRomberg(helper, a, b);
                default:
                    return Integration::IntegrateTrap(helper, a, b);
            }
        }
        
        static double FuncDiffSqr(IRealFunction &f1, IRealFunction &f2, double a, double b, IntegrationMethod method = TRAP)
        {
            RealFuncDiffSqrHelper helper(f1, f2);

            switch (method)
            {
                case SIMPSON:
                    return Integration::IntegrateSimpson(helper, a, b);
                case ROMBERG:
                    return Integration::IntegrateRomberg(helper, a, b);
                default:
                    return Integration::IntegrateTrap(helper, a, b);
            }
        }
    };
}
///////////////////////////   ./include/algorithms/Fourier.h   ///////////////////////////

namespace MML
{
    class Fourier
    {

    };
}


///////////////////////////   ./include/algorithms/RootFinding.h   ///////////////////////////

namespace MML
{
    class RootFinding
    {

    };
}


///////////////////////////   ./include/algorithms/Statistics.h   ///////////////////////////

namespace MML
{
    class Statistics
    {
    public:
        static void moment(Vector<Real> &data, Real &ave, Real &adev, Real &sdev, Real &var,
            Real &skew, Real &curt) {
            int j,n=data.size();
            Real ep=0.0,s,p;
            if (n <= 1) throw("n must be at least 2 in moment");
            s=0.0;
            for (j=0;j<n;j++) s += data[j];
            ave=s/n;
            adev=var=skew=curt=0.0;
            for (j=0;j<n;j++) {
                adev += abs(s=data[j]-ave);
                ep += s;
                var += (p=s*s);
                skew += (p *= s);
                curt += (p *= s);
            }
            adev /= n;
            var=(var-ep*ep/n)/(n-1);
            sdev=sqrt(var);
            if (var != 0.0) {
                skew /= (n*var*sdev);
                curt=curt/(n*var*var)-3.0;
            } else throw("No skew/kurtosis when variance = 0 (in moment)");
        }
        static void avevar(Vector<Real> &data, Real &ave, Real &var) {
            Real s,ep;
            int j,n=data.size();
            ave=0.0;
            for (j=0;j<n;j++) ave += data[j];
            ave /= n;
            var=ep=0.0;
            for (j=0;j<n;j++) {
                s=data[j]-ave;
                ep += s;
                var += s*s;
            }
            var=(var-ep*ep/n)/(n-1);
        }
    };
}


///////////////////////////   ./include/systems/LinAlgEqSystem.h   ///////////////////////////

namespace MML
{
    class LinAlgEqSystem
    {
        // dva ctora - Matrix, Matrix i Matrix, Vector

        // solve by GJ

        // solve by LU - dvije verzije - in place, i verzija koja vraca dekompoziciju
        // perform LU
        //
        // QR, Cholesky
        
        // SVD decomposition
        
        // find eigenvalues

        // ima i Verify - das mu solution i onda vidis kolika je norma razlike u odnosu na zadani rhs

    };
}


///////////////////////////   ./include/systems/DiffEqSystem.h   ///////////////////////////

namespace MML
{
    class DiffEqEqSystem
    {
        // zdaje mu se skup funkcija i pocetni uvjeti
        // dvije vrste - initial & boundary conditions

    };
}



#endif
