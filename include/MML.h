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

#include <functional>

///////////////////////////   ./include/MMLBase.h   ///////////////////////////







template<class T>
inline T SQR(const T a) {return a*a;}

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
    class DefaultParams
    {
        public:
        static inline const double DerivationDefaultStep = 1e-6;
        
        static inline const int    IntegrateTrapJMAX = 20;
        static inline const double IntegrateTrapEPS = 1.0e-10;

        static inline const int    IntegrateSimpJMAX = 20;
        static inline const double IntegrateSimpEPS = 1.0e-10;

        static inline const int    IntegrateRombJMAX = 20;
        static inline const double IntegrateRombEPS = 1.0e-10;
    };

	class VectorDimensionError : public std::invalid_argument
	{
		public:
		std::string _operation;

		VectorDimensionError(std::string inOperation, std::string inMessage) : std::invalid_argument(inMessage)
		{
			_operation = inOperation;
		}
	};	
	class MatrixAccessBoundsError : public std::out_of_range 
	{
		public:
		std::string _operation;

		MatrixAccessBoundsError(std::string inOperation, std::string inMessage) : std::out_of_range(inMessage)
		{
			_operation = inOperation; 
		}
	};

	class MatrixDimensionError : public std::invalid_argument
	{
		public:
		std::string _operation;

		MatrixDimensionError(std::string inOperation, std::string inMessage) : std::invalid_argument(inMessage)
		{
			_operation = inOperation;
		}
	};

	class SingularMatrixError : public std::domain_error
	{
		public:
		std::string _operation;

		SingularMatrixError(std::string inOperation) : std::domain_error(inOperation)
		{
			_operation = inOperation;
		}
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

///////////////////////////   ./include/utilities/TabulatedValues.h   ///////////////////////////

namespace MML
{
    ///////////////////////////////////////////////   Interfaces    ///////////////////////////////////////////
    class ITabulatedValues1D
    {
    public:
        virtual int    NumValues() = 0;

        virtual double Value(int ind) = 0;
        virtual double Value(int ind, double &outX) = 0;

        virtual ~ITabulatedValues1D() {}
    };

    class ITabulatedValues2D
    {
    public:
        virtual int    NumValues1() = 0;
        virtual int    NumValues2() = 0;

        virtual double Value(int ind1, int in2) = 0;
        virtual double Value(int ind1, int in2, double &outX, double &outY) = 0;

        virtual ~ITabulatedValues2D() {}
    };

    class ITabulatedValues3D
    {
    public:
        virtual int    NumValues1() = 0;
        virtual int    NumValues2() = 0;
        virtual int    NumValues3() = 0;

        virtual double Value(int ind1, int ind2, int ind3) = 0;
        virtual double Value(int ind1, int ind2, int ind3, double &outX, double &outY, double &outZ) = 0;

        virtual ~ITabulatedValues3D() {}
    };

    ///////////////////////////////////////////////   Implementations    ///////////////////////////////////////////
    class TabulatedValues1DEqualSpacing : public ITabulatedValues1D
    {
    private:
        double _x1, _x2;
        std::vector<double> _funcValues;

    public:
        TabulatedValues1DEqualSpacing() : _x1{0.0}, _x2{0.0}  { }
        TabulatedValues1DEqualSpacing(double x1, double x2, std::vector<double> values) :  _x1{x1}, _x2{x2}, _funcValues{values}  { }

        int    NumValues() { return (int) _funcValues.size(); }

        double Value(int ind)
        {
            if( !(0 <= ind && ind < NumValues()) )
                throw std::out_of_range("wrong index");
            else
                return _funcValues[ind];
        }

        double Value(int ind, double &outX) 
        {
            if( !(0 <= ind && ind < NumValues()) )
                throw std::out_of_range("wrong index");
            else
            {
                outX = _x1 + (_x2 - _x1) / NumValues() * ind;
                return _funcValues[ind];
            }
        }
    };

    class TabulatedValues1D : public ITabulatedValues1D
    {
    private:
        double _x1, _x2;
        std::vector<double> _xValues;
        std::vector<double> _funcValues;

    public:
        TabulatedValues1D() : _x1{0.0}, _x2{0.0}  { }
        TabulatedValues1D(std::vector<double> xValues, std::vector<double> funcValues) : _xValues{xValues}, _funcValues{funcValues}  
        { 
            if( _xValues.size() != _funcValues.size() )
                throw std::out_of_range("must be equal");

            _x1 = _xValues[0];
            _x2 = _xValues[NumValues()-1];
        }

        int    NumValues() { return (int) _funcValues.size(); }

        double Value(int ind)
        {
            if( !(0 <= ind && ind < NumValues()) )
                throw std::out_of_range("wrong index");
            else
                return _funcValues[ind];
        }

        double Value(int ind, double &outX) 
        {
            if( !(0 <= ind && ind < NumValues()) )
                throw std::out_of_range("wrong index");
            else
            {
                outX = _xValues[ind];
                return _funcValues[ind];
            }
        }        
    };
}

///////////////////////////   ./include/basic_types/Algebra.h   ///////////////////////////
// group

// Z6 group

// permutation group

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
    // množenje vektora skalarom
};

// Hilbert space, kompleksni field!

template<class _VecSpace>
class VectorFromVecSpace
{
    // vraća skalarni produkt definiran u template param
};

void f()
{

}///////////////////////////   ./include/basic_types/Vector.h   ///////////////////////////


namespace MML
{
    template<class _Type>
    class Vector
    {
        private:
            std::vector<_Type> _elems;

        public:
            Vector(size_t n) : _elems(n) {}
            Vector(size_t n, _Type val) : _elems(n, val) {}
            Vector(std::vector<_Type> values) : _elems(values) {}
            Vector(std::initializer_list<_Type> list) : _elems(list) {}
            Vector(int n, _Type *vals) : _elems(n) 
            {
                for(int i=0; i<n; ++i) 
                    _elems[i] = vals[i];
            }

            static Vector GetUnitVector(int dimVec, int indUnit)
            {
                // TODO - što ako je indUnit neispravan?
                Vector ret(dimVec);
                ret[indUnit] = _Type{1.0};
                return ret;
            }

            _Type& operator[](int n)       { return _elems[n]; }
            _Type  operator[](int n) const { return _elems[n]; }
            
            size_t size() const { return _elems.size(); }

            _Type ScalarProductCartesian(Vector &b)
            {
                if (size() != b.size() )
                    throw VectorDimensionError("Vector::ScalarProductCartesian", "wrong dim");

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

            bool IsEqual(const Vector &b, _Type eps=1e-15) const
            {
                if (size() != b.size() )
                    throw VectorDimensionError("Vector::IsEqual", "wrong dim");

                for( int i=0; i<size(); i++ )
                {
                    if( std::abs((*this)[i] - b[i]) > eps )
                        return false;
                }
                return true;
            }

            static bool AreEqual(const Vector &a,const Vector &b, _Type eps=1e-15)
            {
                if (a.size() != b.size() )
                    throw VectorDimensionError("Vector::IsEqual", "wrong dim");

                for( int i=0; i<a.size(); i++ )
                {
                    if( std::abs(a[i] - b[i]) > eps )
                        return false;
                }
                return true;
            }
            Vector operator+(const Vector &b ) const
            {
                if (size() != b.size() )
                    throw VectorDimensionError("Vector::op+", "wrong dim");

                Vector ret(b.size());;
                for(int i=0; i<b.size(); i++)
                    ret._elems[i] = (*this)[i] + b._elems[i];
                return ret;
            }

            Vector operator-(const Vector &b ) const
            {
                if (size() != b.size() )
                    throw VectorDimensionError("Vector::op-", "wrong dim");

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

            // enabling Complex - Double/Float operations
            // friend Vector<Complex> operator+(const Vector<Complex>& a, const Vector<double>& b)
            // {
            //   Vector<Complex> ret(b.size());;
            //   for (int i = 0; i < b.size(); i++)
            //     ret[i] = a[i] + b[i];
            //   return ret;
            // }

            // friend Vector<Complex> operator-(const Vector<Complex>& a, const Vector<double>& b)
            // {
            //   Vector<Complex> ret(b.size());;
            //   for (int i = 0; i < b.size(); i++)
            //     ret[i] = a[i] - b[i];
            //   return ret;
            // }

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

///////////////////////////   ./include/basic_types/VectorN.h   ///////////////////////////

namespace MML
{
    template<class _Type, int N> 
    class VectorN
    {
        protected:
            _Type  _val[N] = {0};

        public:
        VectorN() {}
        VectorN(_Type init_val) : _val{init_val}     {}

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

        static VectorN GetUnitVector(int indUnit)
        {
            VectorN ret;
            ret[indUnit] = 1.0;
            return ret;
        }

        _Type& operator[](int n)       { return _val[n]; }
        _Type  operator[](int n) const { return _val[n]; }

        int size() const { return N; }

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

        bool IsEqual(VectorN &b, _Type eps) const
        {
            for( int i=0; i<N; i++ )
            {
                if( fabs((*this)[i] - b[i]) > eps )
                    return false;
            }
            return true;
        }

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

///////////////////////////   ./include/basic_types/Matrix.h   ///////////////////////////


namespace MML
{
    template<class _Type>
    class Matrix
    {
    private:
        size_t _rows;
        size_t _cols;
        _Type **_ptrData;

    public:
        Matrix() : _rows(0), _cols(0),  _ptrData{nullptr} {} 
        Matrix(size_t rows, size_t cols) : _rows(rows), _cols(cols)
        {
            _ptrData = new _Type *[_rows];
            for (size_t i = 0; i < _rows; ++i)
                _ptrData[i] = new _Type[_cols];
            for (size_t i = 0; i < _rows; ++i)
                for (size_t j = 0; j < _cols; ++j)
                    _ptrData[i][j] = 0;
        }
        Matrix(size_t rows, size_t cols, std::initializer_list<_Type> values) : _rows(rows), _cols(cols)
        {
            _ptrData = new _Type *[_rows];
            for (size_t i = 0; i < _rows; ++i)
                _ptrData[i] = new _Type[_cols];
            
            auto val = values.begin();
            for (size_t i = 0; i < _rows; ++i)
                for (size_t j = 0; j < _cols; ++j)
                    if( val != values.end() )
                    {
                        _ptrData[i][j] = *val;
                        ++val;
                    }
                    else
                        _ptrData[i][j] = 0.0;
        }
        Matrix(const Matrix &m) : _rows(m._rows), _cols(m._cols)
        {
            _ptrData = new _Type *[_rows];
            for (size_t i = 0; i < _rows; ++i)
                _ptrData[i] = new _Type[_cols];
            for (size_t i = 0; i < _rows; ++i)
                for (size_t j = 0; j < _cols; ++j)
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
            for (size_t i = 0; i < _rows; ++i)
                if (_ptrData != nullptr && _ptrData[i] != nullptr)
                    delete[] _ptrData[i];
            if (_ptrData != nullptr)
                delete[] _ptrData;
        }

        static Matrix GetUnitMatrix(int dim)
        {
            Matrix unitMat(dim, dim);
            
            for( int i=0; i<dim; i++ )
                unitMat._ptrData[i][i] = 1.0;

            return unitMat;            
        }

        int RowNum() const { return (int) _rows; }
        int ColNum() const { return (int) _cols; }

        void MakeUnitMatrix(void)
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
                throw MatrixDimensionError("MakeUnitMatrix", "must be square");
        }

        void Resize(size_t rows, size_t cols) 
        {
            for (size_t i = 0; i < _rows; ++i)
                if (_ptrData != nullptr && _ptrData[i] != nullptr)
                    delete[] _ptrData[i];
            if (_ptrData != nullptr)
                delete[] _ptrData;

            _rows = rows;
            _cols = cols;
            _ptrData = new _Type *[_rows];
            for (size_t i = 0; i < _rows; ++i)
                _ptrData[i] = new _Type[_cols];
            for (size_t i = 0; i < _rows; ++i)
                for (size_t j = 0; j < _cols; ++j)
                    _ptrData[i][j] = 0;
        }

        static Matrix RowMatrixFromVector(const Vector<_Type> &b)
        {
            Matrix ret(1, b.size());
            for( int i=0; i<b.size(); i++)
                ret[0][i] = b[i];

            return ret;
        }

        static Matrix ColumnMatrixFromVector(const Vector<_Type> &b)
        {
            Matrix ret(b.size(), 1);
            for( int i=0; i<b.size(); i++)
                ret[i][0] = b[i];

            return ret;
        }

        static Vector<_Type> VectorFromRow(const Matrix &a, int rowInd)
        {
            Vector<_Type> ret(a.ColNum());
            for( int i=0; i<a.ColNum(); i++)
                ret[i] = a(rowInd,i);

            return ret;
        }

        static Vector<_Type> VectorFromColumn(const Matrix &a, int colInd)
        {
            Vector<_Type> ret(a.RowNum());
            for( int i=0; i<a.RowNum(); i++)
                ret[i] = a(i,colInd);

            return ret;
        }

        static Vector<_Type> VectorFromDiagonal(const Matrix &a)
        {
            // TODO verify square matrix
            Vector<_Type> ret(a.RowNum());
            for( int i=0; i<a.RowNum(); i++)
                ret[i] = a(i,i);

            return ret;
        }

        bool IsEqual(const Matrix &b, _Type eps=1e-15) const
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

        bool AreEqual(const Matrix &a, const Matrix &b, _Type eps=1e-15) const
        {
            if( a.RowNum() != b.RowNum() || a.ColNum() != b.ColNum() )
                return false;

            for( int i=0; i<a.RowNum(); i++ )
                for( int j=0; j<a.ColNum(); j++ )
                {
                    if( std::abs(a._ptrData[i][j] - b._ptrData[i][j]) > eps )
                        return false;
                }
                
            return true;
        }        

        Matrix &operator=(const Matrix &m)
        {
            if (this == &m)
                return *this;

            if (_rows != m._rows || _cols != m._cols)
            {
                for (size_t i = 0; i < _rows; ++i)
                    delete[] _ptrData[i];
                delete[] _ptrData;

                _rows = m._rows;
                _cols = m._cols;
                _ptrData = new _Type *[_rows];
                for (size_t i = 0; i < _rows; ++i)
                    _ptrData[i] = new _Type[_cols];
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

        _Type* operator[](int i)        { return _ptrData[i]; }
        
        _Type  operator()(int i, int j) const { return _ptrData[i][j]; }
        _Type& operator()(int i, int j)       { return _ptrData[i][j]; }        

        // version with checking bounds
        _Type  ElemAt(int i, int j) const 
        { 
            if( i<0 || i>=RowNum() || j<0 || j>=ColNum() )
                throw MatrixAccessBoundsError("Matrix::ElemAt", "i=, j=, rows=, cols=");

            return _ptrData[i][j]; 
        }
        _Type& ElemAt(int i, int j)       
        {
            if( i<0 || i>=RowNum() || j<0 || j>=ColNum() )
                throw MatrixAccessBoundsError("MatrixElemAt", "i=, j=, rows=, cols=");

            return _ptrData[i][j]; 
        }

        Matrix operator+(const Matrix &other) const
        {
            if (_rows != other._rows || _cols != other._cols)
                throw MatrixDimensionError("op+", "wrong dim");

            Matrix temp(_rows, _cols);
            for (size_t i = 0; i < _rows; i++)
                for (size_t j = 0; j < _cols; j++)
                    temp._ptrData[i][j] = other._ptrData[i][j] + _ptrData[i][j];

            return temp;
        }

        Matrix operator-(const Matrix &other) const
        {
            if (_rows != other._rows || _cols != other._cols)
                throw MatrixDimensionError("op-", "wrong dim");

            Matrix temp(_rows, _cols);
            for (int i = 0; i < _rows; i++)
                for (int j = 0; j < _cols; j++)
                    temp._ptrData[i][j] =  _ptrData[i][j] - other._ptrData[i][j];

            return temp;
        }

        Matrix  operator*( const Matrix &b ) const
        {
            if( ColNum() != b.RowNum() )
                throw MatrixDimensionError("op*", "wrong dim");

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
                throw MatrixDimensionError("Matrix * Vector<_Type>", "wrong dim");

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
                //std::string error = std::format("Hello {}!\n", "world");
                throw MatrixDimensionError("Vector<_Type> * Matrix", "Vector<_Type> dim = N, Matrix row num = M");

            Vector<_Type>	ret(b.ColNum());
            for( int i=0; i<b.ColNum(); i++ )
            {
                ret[i] = 0;
                for( int j=0; j<b.RowNum(); j++ )
                    ret[i] += a[i] * b(i,j);
            }

            return ret;
        }

        void Invert()
        {
            if( RowNum() != ColNum() ) 
                throw MatrixDimensionError("Matrix::Invert", "");

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
                throw MatrixDimensionError("Matrix::GetInvert", "");

            Matrix a(*this);              // making a copy, where inverse will be stored at the end
            
            a.Invert();
            
            return a;
        }     

        void Transpose()
        {
            // check dimensions - in place Transpose only for square matrices
            if( RowNum() != ColNum() ) 
                throw MatrixDimensionError("Matrix::Transpose", "");

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

        void Print(std::ostream& stream, int width, int precision) const
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
    };
}

///////////////////////////   ./include/basic_types/Mat3D.h   ///////////////////////////

namespace MML
{
    template <class T>
    class Mat3D {
    private:
        int nn;
        int mm;
        int kk;
        T ***v;
    public:
        Mat3D(): nn(0), mm(0), kk(0), v(NULL) {}

        Mat3D(int n, int m, int k) : nn(n), mm(m), kk(k), v(new T**[n])
        {
            int i,j;
            v[0] = new T*[n*m];
            v[0][0] = new T[n*m*k];
            for(j=1; j<m; j++) v[0][j] = v[0][j-1] + k;
            for(i=1; i<n; i++) {
                v[i] = v[i-1] + m;
                v[i][0] = v[i-1][0] + m*k;
                for(j=1; j<m; j++) v[i][j] = v[i][j-1] + k;
            }
        }        

        ~Mat3D()
        {
            if (v != NULL) {
                delete[] (v[0][0]);
                delete[] (v[0]);
                delete[] (v);
            }
        }

        //subscripting: pointer to row i
        inline T**              operator[](const int i)       { return v[i]; }
        inline const T* const * operator[](const int i) const { return v[i]; }

        inline int dim1() const { return nn; }
        inline int dim2() const { return mm; }
        inline int dim3() const { return kk; }
    };
}

///////////////////////////   ./include/basic_types/MatrixNM.h   ///////////////////////////


namespace MML
{
    template <class _Type, int N, int M>
    class MatrixNM
    {
    public:
        _Type _vals[N][M] = {{0}};

    public:
        MatrixNM() {}
        MatrixNM(std::initializer_list<_Type> values) 
        {
            auto val = values.begin();
            for (size_t i = 0; i < RowNum(); ++i)
                for (size_t j = 0; j < ColNum(); ++j)
                    if( val != values.end() )
                    {
                        _vals[i][j] = *val;
                        ++val;
                    }
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

        static MatrixNM GetUnitMatrix() 
        {
            MatrixNM unitMat;
            
            for( int i=0; i<N; i++ )
                unitMat._vals[i][i] = 1.0;

            return unitMat;
        }

        int RowNum() const { return (int) N; }
        int ColNum() const { return (int) M; }

        void MakeUnitMatrix(void)
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
        
        static MatrixNM<_Type, 1,M> RowMatrixFromVector(const VectorN<_Type, M> &b)
        {
            MatrixNM<_Type,1,M>  ret;
            for( int j=0; j<M; j++)
                ret._vals[0][j] = b[j];

            return ret;
        }
        
        static MatrixNM<_Type, N, 1> ColumnMatrixFromVector(const VectorN<_Type, N> &b)
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

        bool IsEqual(const MatrixNM &b, _Type eps=1e-15) const
        {
            for( int i=0; i<RowNum(); i++ )
                for( int j=0; j<ColNum(); j++ )
                    if( fabs(_vals[i][j] - b._vals[i][j]) > eps )
                        return false;
                
            return true;
        }

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
                throw MatrixAccessBoundsError("MatrixElemAt", "i=, j=, rows=, cols=");

            return _vals[i][j]; 
        }
        _Type& ElemAt(int i, int j)       
        {
            if( i<0 || i>=RowNum() || j<0 || j>=ColNum() )
                throw MatrixAccessBoundsError("MatrixElemAt", "i=, j=, rows=, cols=");

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

        friend VectorN<_Type, N> operator*( const MatrixNM<_Type, N,M> &a, const VectorN<_Type, M> &b )
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

        friend VectorN<_Type, M> operator*( const VectorN<_Type, N> &a, const MatrixNM<_Type, N,M> &b )
        {
            int	i, j;
            VectorN<_Type, M>	ret;

            for( i=0; i<b.ColNum(); i++ )
            {
                ret[i] = 0;
                for( j=0; j<b.RowNum(); j++ )
                    ret[i] += a[i] * b._vals[j][i];
            }

            return ret;
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

        friend std::ostream& operator<<(std::ostream& stream, const MatrixNM &a)
        {
            a.Print(stream, 10, 3);

            return stream;
        }     

        void Invert()
        {
            if( RowNum() != ColNum() ) 
                throw MatrixDimensionError("MatrixNM::Invert", "");

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
                throw MatrixDimensionError("MatrixNM::GetInvert", "");

            MatrixNM a(*this);              // making a copy, where inverse will be stored at the end
            
            a.Invert();
            
            return a;
        }     

        void Transpose()
        {
            // check dimensions - in place Transpose only for square matrices
            if( RowNum() != ColNum() ) 
                throw MatrixDimensionError("MatrixNM::Transpose", "");

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
    };
}

///////////////////////////   ./include/basic_types/Polynom.h   ///////////////////////////


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

    // ovo ce biti generalni polinom
    template <typename _Field, typename _CoefType = double>
    class Polynom
    {
    private:
        std::vector<_CoefType> _vecCoef;
    public:
        Polynom() {}
        Polynom(const std::vector<_CoefType> &vecCoef) : _vecCoef(vecCoef) {}
        Polynom(std::initializer_list<_CoefType> list) : _vecCoef(list) {}

        Polynom(const Polynom &Copy) : _vecCoef(Copy._vecCoef) {}
        ~Polynom() {}

        Polynom& operator=(const Polynom &Copy) { _vecCoef = Copy._vecCoef; return *this; }

        int GetDegree() const { return (int) _vecCoef.size() - 1; }

        _Field operator()(const _Field &x) const
        {
            _Field result = 1.0; 
            _Field power = x;
            for (int i = 1; i < _vecCoef.size(); i++)
            {
                result = result + _vecCoef[i] * power;
                power = power * x;
            }
            return result;
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

        friend Polynom operator*(const Polynom &a, _CoefType b )
        {
            Polynom ret;
            for(int i=0; i<=a.GetDegree(); i++)
                ret._vecCoef[i] = a._vecCoef[i] * b;
            return ret;
        }

        friend Polynom operator*(_CoefType a, const Polynom &b )
        {
            Polynom ret;
            for(int i=0; i<=b.GetDegree(); i++)
                ret._vecCoef[i] = a * b._vecCoef[i];
            return ret;
        }

        friend Polynom operator/(const Polynom &a, _CoefType b)
        {
            Polynom ret;
            for(int i=0; i<=a.GetDegree(); i++)
                ret._vecCoef[i] = a._vecCoef[i] / b;
            return ret;
        }        

        std::ostream& Print(std::ostream& stream, int width, int precision) const
        {
            stream << "[";
            bool first = true;
            for(const _CoefType& x : _vecCoef)
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

        std::string to_string(int width, int precision) const
        {
            std::stringstream str;

            Print(str, width, precision);

            return str.str();
        }

        friend std::ostream& operator<<(std::ostream& stream, Polynom &a)
        {
            a.Print(stream,15,10);

            return stream;
        }        
    };
        
    typedef Polynom<Real, Real>         RealPolynom;
    typedef Polynom<Complex, Complex>   ComplexPolynom;

    typedef Polynom<MatrixNM<Real,2,2>, Real>       MatrixPolynomDim2;
    typedef Polynom<MatrixNM<Real,3,3>, Real>       MatrixPolynomDim3;
    typedef Polynom<MatrixNM<Real,4,4>, Real>       MatrixPolynomDim4;
}

///////////////////////////   ./include/basic_types/MatrixSparse.h   ///////////////////////////

namespace MML
{

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

        Point2Cartesian operator+(const VectorN<Real, 2> &b) const
        {
            return Point2Cartesian(_x + b[0], _y + b[1]);
        }
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

        // TODO - tocka i double, koji je kut prema pozitivnoj x osi

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
        SegmentLine2D(Point2Cartesian pnt1, Point2Cartesian pnt2)
        {
            _point1 = pnt1;
            _point2 = pnt2;
        }

        SegmentLine2D(Point2Cartesian pnt1, Vector2Cartesian direction, double t)
        {
            _point1 = pnt1;
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
        Triangle2D(Point2Cartesian pnt1, Point2Cartesian pnt2, Point2Cartesian pnt3)
        {
            _pnt1 = pnt1;
            _pnt2 = pnt2;
            _pnt3 = pnt3;
        }

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

        std::vector<Point2Cartesian> Points() const { return _points; }
        std::vector<Point2Cartesian>& Points() { return _points; }

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

        Point3Cartesian operator+(const VectorN<Real, 3> &b) const
        {
            return Point3Cartesian(_x + b[0], _y + b[1], _z + b[2]);
        }
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
        // distance Line - Point3
        // nearest point on line
        // pravac koji prolazi kroz tocku i sijece zadani pravac okomito
    };

    class SegmentLine3D
    {
    private:
        Point3Cartesian _point1;
        Point3Cartesian _point2;

    public:
        SegmentLine3D(Point3Cartesian pnt1, Point3Cartesian pnt2)
        {
            _point1 = pnt1;
            _point2 = pnt2;
        }

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
        Plane3D(double seg_x, double seg_y, double seg_z) // Hesseov (normalni) oblik
        {
            Point3Cartesian x(seg_x, 0, 0);
            Point3Cartesian y(0, seg_y, 0);
            Point3Cartesian z(0, 0, seg_z);

            // TODO
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
        bool IntersectionWithTwoPlanes(const Plane3D &plane1, const Plane3D &plane2, Line3D &out_inter_line) const
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
SUŠTINA
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

MovingDynamicalSytem3D earth_around_sun(funkcija ovisnosti pozicije u odnosu na GLOBALNI KARETEZIJEV sustav);
RotatingSystem3D earth_rotation(earth_around_sun);
- za dane dvije koord, lat i long, daje poziciju u odnosu na dani koord sustav
LocalCartesian3D earth_surface(earth_rotation, lat, long);

LorentzInertialMovingFrame observer_moving_frame(vektor smjera, ovisnost pozicije o t); /// moze i (0,0,0,0) - stoi na mjestu
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

///////////////////////////   ./include/interfaces/ILinAlgEqSystemSolver.h   ///////////////////////////

namespace MML
{
    class ILinAlgEqSystemSolver
    {
        // dvije Solve, za Vector i Matrix
    };
}


///////////////////////////   ./include/algorithms/LinAlgEqSolvers.h   ///////////////////////////


namespace MML
{
    class LinAlgEqSolvers
    {
        // staticke funkcije za low level rjesavanje
    };

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
            // permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
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
        void mprove(Vector<Real> &b, Vector<Real> &x)
        {
            // Improves a solution Vector<Real> x[1..n] of the linear set of equations A · X = B. The Matrix<Real>
            // a[1..n][1..n], and the Vector<Real>s b[1..n] and x[1..n] are input, as is the dimension n.
            // Also input is alud[1..n][1..n], the LU decomposition of a as returned by ludcmp, and
            // the Vector<Real> indx[1..n] also returned by that routine. On output, only x[1..n] is modified,
            // to an improved set of values
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
            // Constructs the QR decomposition of a[1..n][1..n]. The upper triangular Matrix<Real> R is returned in the upper triangle of a, except for the diagonal elements of R which are returned in
            // d[1..n]. The orthogonal Matrix<Real> Q is represented as a product of n − 1 Householder matrices
            // Q1 : : : Qn−1, where Qj = 1 − uj ⊗ uj=cj. The ith component of uj is zero for i = 1; : : : ; j − 1
            // while the nonzero components are returned in a[i][j] for i = j; : : : ; n. sing returns as
            // true (1) if singularity is encountered during the decomposition, but the decomposition is still
            // completed in this case; otherwise it returns false (0).
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
        void Solve(Vector<Real> &b, Vector<Real> &x) 
        {
            // Solves the set of n linear equations A · x = b. a[1..n][1..n], c[1..n], and d[1..n] are
            // input as the output of the routine qrdcmp and are not modified. b[1..n] is input as the
            // right-hand side Vector<Real>, and is overwritten with the solution Vector<Real> on output.            
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
            // Matrix<Real> Q·(R+u⊗v). The quantities are dimensioned as r[1..n][1..n], qt[1..n][1..n],
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
            // i and i + 1 of each Matrix<Real>. a and b are the parameters of the rotation: cos θ = a=pa2 + b2,
            // sin θ = b=pa2 + b2.            
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
            // Given a Matrix<Real> a[1..m][1..n], this routine computes its singular value decomposition, A =
            // U·W ·V T . The Matrix<Real> U replaces a on output. 
            // The diagonal Matrix<Real> of singular values W is output as a Vector<Real> w[1..n]. The Matrix<Real> V (not the transpose V T ) is output as v[1..n][1..n].            
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
            // Solves A·X = B for a Vector<Real> X, where A is specified by the arrays u[1..m][1..n], w[1..n],
            // v[1..n][1..n] as returned by svdcmp. m and n are the dimensions of a, and will be equal for
            // square matrices. b[1..m] is the input right-hand side. x[1..n] is the output solution Vector<Real>.
            // No input quantities are destroyed, so the routine may be called sequentially with different b’s.        
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
        int rank(double thresh = -1.) {
            int j,nr=0;
            tsh = (thresh >= 0. ? thresh : 0.5*sqrt(m+n+1.)*w[0]*eps);
            for (j=0;j<n;j++) if (w[j] > tsh) nr++;
            return nr;
        }

        int nullity(double thresh = -1.) {
            int j,nn=0;
            tsh = (thresh >= 0. ? thresh : 0.5*sqrt(m+n+1.)*w[0]*eps);
            for (j=0;j<n;j++) if (w[j] <= tsh) nn++;
            return nn;
        }

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

void eigsrt(MML::Vector<Real> &d, MML::Matrix<Real> *v = NULL)
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
    class IVectorFunction : public IFunction<VectorN<Real, N>, const VectorN<Real, N> &>
    {
        public:
        virtual VectorN<Real, N> operator()(const VectorN<Real, N> &x) const = 0;
        virtual double operator()(VectorN<Real, N> &x, int component) const
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
        virtual double operator()(VectorN<Real, N> &x, int component) const
        {
            VectorN<Real, M> val = (*this)(x);
            return val[component];
        }
    };

    //////////////////////////////////////////////////////////////////////
    template<int N>
    class IParametricCurve : public IFunction<VectorN<Real, N>, double>
    {
        public:
        virtual VectorN<Real, N> operator()(double x) const = 0;

        // GetMixX(), GetMaxX(), može vraćati i infinity
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

        // GetMixX(), GetMaxX(), može vraćati i infinity
        // GetMixY(), GetMaxY(), može vraćati i infinity
        // da je površina omeđena

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
    class ScalarFunctionFromFuncPtr : public IScalarFunction<N>
    {
        public:
        double (*_func)(const VectorN<double, N> &);

        ScalarFunctionFromFuncPtr( double (*inFunc)(const VectorN<double, N> &) ) : _func(inFunc)    {}

        double operator()(const VectorN<double, N> &x) const  { return _func(x); }
    };

    template<int N>
    class ScalarFunctionFromStdFunc : public IScalarFunction<N>
    {
        public:
        std::function<double(const VectorN<double, N> &)> _func;

        ScalarFunctionFromStdFunc(std::function<double(const VectorN<double, N> &)> inFunc) : _func(inFunc)     {}

        double operator()(const VectorN<double, N> &x) const  { return _func(x); }
    };
    
    /////////////////////////    VECTOR FUNCTION N -> N      ///////////////////////////////////
    template<int N>
    class VectorFunctionFromFuncPtr : public IVectorFunction<N>
    {
        public:
        VectorN<double, N> (*_func)(const VectorN<double, N> &);

        VectorFunctionFromFuncPtr( VectorN<double, N> (*inFunc)(const VectorN<double, N> &) ) : _func(inFunc)      {}

        VectorN<double, N> operator()(const VectorN<double, N> &x) const  { return _func(x); }
    };
    template<int N>
    class VectorFunctionFromStdFunc : public IVectorFunction<N>
    {
        public:
        std::function<VectorN<double, N>(const VectorN<double, N> &)> _func;

        VectorFunctionFromStdFunc(std::function<VectorN<double, N>(const VectorN<double, N> &)> &inFunc) : _func(inFunc)    {}

        VectorN<double, N> operator()(const VectorN<double, N> &x) const   { return _func(x); }
    };

   /////////////////////////    VECTOR FUNCTION N -> M      ///////////////////////////////////
    template<int N, int M>
    class VectorFunctionNMFromFuncPtr : public IVectorFunctionNM<N,M>
    {
        public:
        VectorN<double, M> (*_func)(const VectorN<double, N> &);

        VectorFunctionNMFromFuncPtr( VectorN<double, N> (*inFunc)(const VectorN<double, N> &) ) : _func(inFunc)      {}

        VectorN<double, M> operator()(const VectorN<double, N> &x) const  { return _func(x); }
    };
    template<int N, int M>
    class VectorFunctionNMFromStdFunc : public IVectorFunctionNM<N,M>
    {
        public:
        std::function<VectorN<double, M>(const VectorN<double, N> &)> _func;

        VectorFunctionNMFromStdFunc(std::function<VectorN<double, M>(const VectorN<double, N> &)> &inFunc) : _func(inFunc)    {}

        VectorN<double, M> operator()(const VectorN<double, N> &x) const   { return _func(x); }
    };

    //////////////////////////////////////////////////////////////////////
    template<int N>
    class ParametricCurveFromFuncPtr : public IParametricCurve<N>
    {
        public:
        VectorN<double, N> (*_func)(double);

        ParametricCurveFromFuncPtr( VectorN<double, N> (*inFunc)(double) ) : _func(inFunc)    {}

        VectorN<double, N> operator()(double x) const  { return _func(x); }
    };

    template<int N>
    class ParametricCurveFromStdFunc : public IParametricCurve<N>
    {
        public:
        std::function<VectorN<double, N>(double)> _func;

        ParametricCurveFromStdFunc(std::function<VectorN<double, N>(double)> &inFunc) : _func(inFunc)    {}

        VectorN<double, N> operator()(double x) const   { return _func(x); }
    };

    //////////////////////////////////////////////////////////////////////
    template<int N>
    class ParametricSurfaceFromFuncPtr : public IParametricSurface<N>
    {
        public:
        VectorN<double, N> (*_func)(double u, double w);

        ParametricSurfaceFromFuncPtr( VectorN<double, N> (*inFunc)(double u, double w) ) : _func(inFunc)    {}

        VectorN<double, N> operator()(double u, double w) const  { return _func(u,w); }
    };

    template<int N>
    class ParametricSurfaceFromStdFunc : public IParametricSurface<N>
    {
        public:
        std::function<VectorN<double, N>(double u, double w)> _func;

        ParametricSurfaceFromStdFunc(std::function<VectorN<double, N>(double u, double w)> &inFunc) : _func(inFunc)    {}

        VectorN<double, N> operator()(double u, double w) const   { return _func(u,w); }
    };    
} // end namespace

///////////////////////////   ./include/basic_types/Functions.h   ///////////////////////////


namespace MML
{
    namespace Curves
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
        // static inline Abs(Real x) { return abs(x); }
        // static inline Floor(Real x) { return floor(x); }
        // static inline Ceil(Real x) { return ceil(x); }
        // static inline Round(Real x) { return round(x); }
        // static inline Sign(Real x) { return sign(x); }

    }
}

///////////////////////////   ./include/basic_types/Curves.h   ///////////////////////////


namespace MML
{
    namespace Curves
    {
        static ParametricCurveFromFuncPtr<3> helix_curve([](double t) { return MML::VectorN<Real, 3>{cos(t), sin(t), t}; });
    }
}

///////////////////////////   ./include/basic_types/InterpolatedFunction.h   ///////////////////////////


namespace MML
{
    class InterpolatedRealFunction : public IRealFunction
    {
        std::unique_ptr<ITabulatedValues1D> _values;

        public:
        InterpolatedRealFunction() 
        {

        }

        double operator()(double x) const    { return 0.0; }
    };

    class InterpolatedScalarFunction2D : public IScalarFunction<2>
    {
        public:
        InterpolatedScalarFunction2D() {}

        double operator()(const VectorN<Real, 2> &x) const    { return 0.0; }
        virtual double operator()(double u, double w)
        {
            VectorN<Real, 2> coord{u,w};

            return operator()(coord);
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

    template<int N>
    class InterpolatedCurve : public IParametricCurve<N>
    {
        public:
        InterpolatedCurve() {}

        VectorN<Real, N> operator()(double t) const    { return VectorN<Real, N>{}; }
    };

    template<int N>
    class InterpolatedSurface : public IParametricSurface<N>
    {
        public:
        InterpolatedSurface() {}

        VectorN<Real, N> operator()(const VectorN<Real, 2> &x) const    { return VectorN<Real, N>{}; }
    };       
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

///////////////////////////   ./include/algorithms/Derivation.h   ///////////////////////////



namespace MML
{
    class Derivation
    {
    public:
        static double make_xph_representable(double x, double  h)
        {
            using std::numeric_limits;
            // Redefine h so that x + h is representable. Not using this trick leads to large error.
            // The compiler flag -ffast-math evaporates these operations . . .
            double temp = x + h;
            h = temp - x;
//            double g = boost::math::nextafter(x, (numeric_limits<double>::max)()) - x;
            // Handle the case x + h == x:
            if (h == 0)
            {
  //              h = boost::math::nextafter(x, (numeric_limits<double>::max)()) - x;
            }
            return h;
        }

        /********************************************************************************************************************/
        /********                               Numerical derivatives of FIRST order                                 ********/
        /********************************************************************************************************************/
        static double NDer1(const MML::IRealFunction &f, double x, double* error = nullptr)
        {
            const double eps = std::numeric_limits<double>::epsilon();
            // Error bound ~eps^1/2
            // Note that this estimate of h differs from the best estimate by a factor of sqrt((|f(x)| + |f(x+h)|)/|f''(x)|).
            // Since this factor is invariant under the scaling f -> kf, then we are somewhat justified in approximating it by 1.
            // This approximation will get better as we move to higher orders of accuracy.
            double h = 2 * std::sqrt(eps);
            //h = make_xph_representable(x, h);

            double yh = f(x + h);
            double y0 = f(x);
            double diff = yh - y0;
            if (error)
            {
                double ym = f(x - h);
                double ypph = std::abs(yh - 2 * y0 + ym) / h;
                // h*|f''(x)|*0.5 + (|f(x+h)+|f(x)|)*eps/h
                *error = ypph / 2 + (std::abs(yh) + std::abs(y0))*eps / h;
            }
            return diff / h;
        }
     
        static double NDer1(const MML::IRealFunction &f, double x, double h, double* error = nullptr)
        {
            const double eps = std::numeric_limits<double>::epsilon();

            double yh = f(x + h);
            double y0 = f(x);
            double diff = yh - y0;
            if (error)
            {
                double ym = f(x - h);
                double ypph = std::abs(yh - 2 * y0 + ym) / h;
                *error = ypph / 2 + (std::abs(yh) + std::abs(y0))*eps / h;
            }
            return diff / h;
        }

        template <int N>
        static double NDer1Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, double *error = nullptr)
        {
            double h = 2 * sqrt(std::numeric_limits<double>::epsilon());

            return NDer1Partial(f, deriv_index, point, h, error);
        }

        template <int N>
        static double NDer1Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, double h, double *error = nullptr)
        {
            const double eps = std::numeric_limits<double>::epsilon();

            double     orig_x = point[deriv_index];

            VectorN<Real, N> x = point;
            double y0 = f(x);

            x[deriv_index] = orig_x + h;
            double yh = f(x);

            double diff = yh - y0;
            if (error)
            {
                x[deriv_index] = orig_x - h;
                double ym = f(x);
                double ypph = std::abs(yh - 2 * y0 + ym) / h;
                *error = ypph / 2 + (std::abs(yh) + std::abs(y0))*eps / h;
            }
            return diff / h;
        } 

        template <int N>
        static VectorN<Real, N> NDer1PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            double h = 2 * sqrt(std::numeric_limits<double>::epsilon());

            return NDer1PartialByAll(f, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer1PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, double h, VectorN<Real, N> *error = nullptr)
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

        template <int N>
        static double NDer1Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, double *error = nullptr)
        {
            double h = 2 * sqrt(std::numeric_limits<double>::epsilon());

            return NDer1Partial(f, func_index, deriv_index, point, h, error);
        }

        template <int N>
        static double NDer1Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, double h, double *error = nullptr)
        {
            const double eps = std::numeric_limits<double>::epsilon();

            VectorN<Real, N> x{point};
            double x_orig = x[deriv_index];
            double y0 = f(x)[func_index];

            x[deriv_index] = x_orig + h;
            double yh = f(x)[func_index];

            double diff = yh - y0;
            if (error)
            {
                x[deriv_index] = x_orig - h;
                double ym = f(x)[func_index];
                double ypph = std::abs(yh - 2 * y0 + ym) / h;
                *error = ypph / 2 + (std::abs(yh) + std::abs(y0))*eps / h;
            }
            return diff / h;
        }

        template <int N>
        static VectorN<Real, N> NDer1PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            double h = 2 * sqrt(std::numeric_limits<double>::epsilon());

            return NDer1PartialByAll(f, func_index, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer1PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, double h, VectorN<Real, N> *error = nullptr)
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
            double h = 2 * sqrt(std::numeric_limits<double>::epsilon());

            return NDer1PartialAllByAll(f, point, h, error);
        }

        template <int N>
        static MatrixNM<Real, N,N> NDer1PartialAllByAll(const IVectorFunction<N> &f, const VectorN<Real, N> &point, double h, MatrixNM<Real, N,N> *error = nullptr)
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

        /********************************************************************************************************************/
        /********                               Numerical derivatives of SECOND order                                 ********/
        /********************************************************************************************************************/
        static double NDer2(const MML::IRealFunction &f, double x, double* error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            // Error bound ~eps^2/3
            // See the previous discussion to understand determination of h and the error bound.
            // Series[(f[x+h] - f[x-h])/(2*h), {h, 0, 4}]
            double h = std::pow(3 * eps, static_cast<double>(1) / static_cast<double>(3));
            //h = detail::make_xph_representable(x, h);

            double yh = f(x + h);
            double ymh = f(x - h);
            double diff = yh - ymh;
            if (error)
            {
                double yth = f(x + 2 * h);
                double ymth = f(x - 2 * h);
                *error = eps * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((yth - ymth) / 2 - diff) / (6 * h);
            }

            return diff / (2 * h);
        } 

        static double NDer2(const MML::IRealFunction &f, double x, double h, double* error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();

            double yh = f(x + h);
            double ymh = f(x - h);
            double diff = yh - ymh;
            if (error)
            {
                double y2h = f(x + 2 * h);
                double ym2h = f(x - 2 * h);
                *error = eps * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
            }

            return diff / (2 * h);
        }  

        template <int N>
        static double NDer2Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, double *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(3 * eps, static_cast<double>(1) / static_cast<double>(3));

            return NDer2Partial(f, deriv_index, point, h, error);
        }

        template <int N>
        static double NDer2Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, double h, double *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();

            double     orig_x = point[deriv_index];

            VectorN<Real, N> x{point};
            x[deriv_index] = orig_x + h;
            double yh = f(x);

            x[deriv_index] = orig_x - h;
            double ymh = f(x);

            double diff = yh - ymh;

            if (error)
            {
                x[deriv_index] = orig_x + 2 * h;
                double y2h = f(x);

                x[deriv_index] = orig_x - 2 * h;
                double ym2h = f(x);
                *error = eps * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
            }

            return diff / (2 * h);
        }

        template <int N>
        static VectorN<Real, N> NDer2PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(3 * eps, static_cast<double>(1) / static_cast<double>(3));

            return NDer2PartialByAll(f, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer2PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, double h, VectorN<Real, N> *error = nullptr)
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

        template <int N>
        static double NDer2Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, double *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(3 * eps, static_cast<double>(1) / static_cast<double>(3));

            return NDer2Partial(f, func_index, deriv_index, point, h, error);
        }        

        template <int N>
        static double NDer2Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, double h, double *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();

            double     orig_x = point[deriv_index];

            VectorN<Real, N> x{point};
            x[deriv_index] = orig_x + h;
            double yh = f(x)[func_index];

            x[deriv_index] = orig_x - h;
            double ymh = f(x)[func_index];

            double diff = yh - ymh;

            if (error)
            {
                x[deriv_index] = orig_x + 2 * h;
                double y2h = f(x)[func_index];

                x[deriv_index] = orig_x - 2 * h;
                double ym2h = f(x)[func_index];                

                *error = eps * (std::abs(yh) + std::abs(ymh)) / (2 * h) + std::abs((y2h - ym2h) / 2 - diff) / (6 * h);
            }

            return diff / (2 * h);
        }

        template <int N>
        static VectorN<Real, N> NDer2PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(3 * eps, static_cast<double>(1) / static_cast<double>(3));

            return NDer2PartialByAll(f, func_index, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer2PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, double h, VectorN<Real, N> *error = nullptr)
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
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(3 * eps, static_cast<double>(1) / static_cast<double>(3));

            return NDer2PartialAllByAll(f, point, h, error);
        }

        template <int N>
        static MatrixNM<Real, N,N> NDer2PartialAllByAll(const IVectorFunction<N> &f, const VectorN<Real, N> &point, double h, MatrixNM<Real, N,N> *error = nullptr)
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

        /********************************************************************************************************************/
        /********                               Numerical derivatives of FOURTH order                                 ********/
        /********************************************************************************************************************/
        static double NDer4(const MML::IRealFunction &f, double x, double* error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            // Error bound ~eps^4/5
            double h = std::pow(11.25*eps, (double)1 / (double)5);
            //h = detail::make_xph_representable(x, h);
            
            double yh = f(x + h);
            double ymh = f(x - h);
            double y2h = f(x + 2 * h);
            double ym2h = f(x - 2 * h);
            
            double y2 = ym2h - y2h;
            double y1 = yh - ymh;
            
            if (error)
            {
                // Mathematica code to extract the remainder:
                // Series[(f[x-2*h]+ 8*f[x+h] - 8*f[x-h] - f[x+2*h])/(12*h), {h, 0, 7}]
                double y3h = f(x + 3 * h);
                double ym3h = f(x - 3 * h);
                
                // Error from fifth derivative:
                *error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
                // Error from function evaluation:
                *error += eps * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
            }
            return (y2 + 8 * y1) / (12 * h);
        }

        static double NDer4(const MML::IRealFunction &f, double x, double h, double* error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            
            double yh = f(x + h);
            double ymh = f(x - h);
            double y2h = f(x + 2 * h);
            double ym2h = f(x - 2 * h);
            
            double y2 = ym2h - y2h;
            double y1 = yh - ymh;
            
            if (error)
            {
                double y3h = f(x + 3 * h);
                double ym3h = f(x - 3 * h);

                *error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);
                *error += eps * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
            }
            return (y2 + 8 * y1) / (12 * h);
        }

        template <int N>
        static double NDer4Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, double *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(11.25*eps, (double)1 / (double)5);

            return NDer4Partial(f, deriv_index, point, h, error);
        }

        template <int N>
        static double NDer4Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, double h, double *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();

            double     orig_x = point[deriv_index];

            VectorN<Real, N> x{point};
            x[deriv_index] = orig_x + h;
            double yh = f(x);

            x[deriv_index] = orig_x - h;
            double ymh = f(x);

            x[deriv_index] = orig_x + 2 * h;
            double y2h = f(x);

            x[deriv_index] = orig_x - 2 * h;
            double ym2h = f(x);

            double y2 = ym2h - y2h;
            double y1 = yh - ymh;
            
            if (error)
            {
                x[deriv_index] = orig_x + 3 * h;
                double y3h = f(x);

                x[deriv_index] = orig_x - 3 * h;
                double ym3h = f(x);

                *error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);

                *error += eps * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
            }
            return (y2 + 8 * y1) / (12 * h);
        }

        template <int N>
        static VectorN<Real, N> NDer4PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(11.25*eps, (double)1 / (double)5);

            return NDer4PartialByAll(f, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer4PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, double h, VectorN<Real, N> *error = nullptr)
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

        template <int N>
        static double NDer4Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, double *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(11.25*eps, (double)1 / (double)5);

            return NDer4Partial(f, func_index, deriv_index, point, h, error);
        }

        template <int N>
        static double NDer4Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, double h, double *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();

            double     orig_x = point[deriv_index];

            VectorN<Real, N> x{point};
            x[deriv_index] = orig_x + h;
            double yh = f(x)[func_index];

            x[deriv_index] = orig_x - h;
            double ymh = f(x)[func_index];

            x[deriv_index] = orig_x + 2 * h;
            double y2h = f(x)[func_index];

            x[deriv_index] = orig_x - 2 * h;
            double ym2h = f(x)[func_index];

            double y2 = ym2h - y2h;
            double y1 = yh - ymh;
            
            if (error)
            {
                x[deriv_index] = orig_x + 3 * h;
                double y3h = f(x)[func_index];

                x[deriv_index] = orig_x - 3 * h;
                double ym3h = f(x)[func_index];

                *error = std::abs((y3h - ym3h) / 2 + 2 * (ym2h - y2h) + 5 * (yh - ymh) / 2) / (30 * h);

                *error += eps * (std::abs(y2h) + std::abs(ym2h) + 8 * (std::abs(ymh) + std::abs(yh))) / (12 * h);
            }
            return (y2 + 8 * y1) / (12 * h);
        }

        template <int N>
        static VectorN<Real, N> NDer4PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(11.25*eps, (double)1 / (double)5);

            return NDer4PartialByAll(f, func_index, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer4PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, double h, VectorN<Real, N> *error = nullptr)
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
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(3 * eps, static_cast<double>(1) / static_cast<double>(3));

            return NDer4PartialAllByAll(f, point, h, error);
        }

        template <int N>
        static MatrixNM<Real, N,N> NDer4PartialAllByAll(const IVectorFunction<N> &f, const VectorN<Real, N> &point, double h, MatrixNM<Real, N,N> *error = nullptr)
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

        /********************************************************************************************************************/
        /********                               Numerical derivatives of SIXTH order                                 ********/
        /********************************************************************************************************************/
        static double NDer6(const MML::IRealFunction &f, double x, double* error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            // Error bound ~eps^6/7
            // Error: h^6f^(7)(x)/140 + 5|f(x)|eps/h
            double h = std::pow(eps / 168, (double)1 / (double)7);
            //h = detail::make_xph_representable(x, h);

            double yh = f(x + h);
            double ymh = f(x - h);
            double y1 = yh - ymh;
            double y2 = f(x - 2 * h) - f(x + 2 * h);
            double y3 = f(x + 3 * h) - f(x - 3 * h);

            if (error)
            {
                // Mathematica code to generate fd scheme for 7th derivative:
                // Sum[(-1)^i*Binomial[7, i]*(f[x+(3-i)*h] + f[x+(4-i)*h])/2, {i, 0, 7}]
                // Mathematica to demonstrate that this is a finite difference formula for 7th derivative:
                // Series[(f[x+4*h]-f[x-4*h] + 6*(f[x-3*h] - f[x+3*h]) + 14*(f[x-h] - f[x+h] + f[x+2*h] - f[x-2*h]))/2, {h, 0, 15}]
                double y7 = (f(x + 4 * h) - f(x - 4 * h) - 6 * y3 - 14 * y1 - 14 * y2) / 2;
                *error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh))*eps / h;
            }
            return (y3 + 9 * y2 + 45 * y1) / (60 * h);
        }

        static double NDer6(const MML::IRealFunction &f, double x, double h, double* error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();

            double yh = f(x + h);
            double ymh = f(x - h);
            double y1 = yh - ymh;
            double y2 = f(x - 2 * h) - f(x + 2 * h);
            double y3 = f(x + 3 * h) - f(x - 3 * h);

            if (error)
            {
                double y7 = (f(x + 4 * h) - f(x - 4 * h) - 6 * y3 - 14 * y1 - 14 * y2) / 2;
                *error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh))*eps / h;
            }
            return (y3 + 9 * y2 + 45 * y1) / (60 * h);
        }        

        template <int N>
        static double NDer6Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, double *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(eps / 168, (double)1 / (double)7);

            return NDer6Partial(f, deriv_index, point, h, error);
        }

        template <int N>
        static double NDer6Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, double h, double *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();

            double     orig_x = point[deriv_index];

            VectorN<Real, N> x{point};

            x[deriv_index] = orig_x + h;
            double yh = f(x);

            x[deriv_index] = orig_x - h;
            double ymh = f(x);

            x[deriv_index] = orig_x + 2 * h;
            double y2h = f(x);

            x[deriv_index] = orig_x - 2 * h;
            double ym2h = f(x);

            x[deriv_index] = orig_x + 3 * h;
            double y3h = f(x);

            x[deriv_index] = orig_x - 3 * h;
            double ym3h = f(x);

            double y1 = yh - ymh;
            double y2 = ym2h - y2h;
            double y3 = y3h - ym3h;

            if (error)
            {
                x[deriv_index] = orig_x + 4 * h;
                double y4h = f(x);

                x[deriv_index] = orig_x - 4 * h;
                double ym4h = f(x);

                double y7 = (y4h - ym4h - 6 * y3 - 14 * y1 - 14 * y2) / 2;
                *error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh))*eps / h;
            }
            return (y3 + 9 * y2 + 45 * y1) / (60 * h);
        }

        template <int N>
        static VectorN<Real, N> NDer6PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(eps / 168, (double)1 / (double)7);

            return NDer6PartialByAll(f, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer6PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, double h, VectorN<Real, N> *error = nullptr)
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

        template <int N>
        static double NDer6Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, double *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(eps / 168, (double)1 / (double)7);

            return NDer6Partial(f, func_index, deriv_index, point, h, error);
        }

        template <int N>
        static double NDer6Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, double h, double *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();

            double     orig_x = point[deriv_index];

            VectorN<Real, N> x{point};

            x[deriv_index] = orig_x + h;
            double yh = f(x)[func_index];

            x[deriv_index] = orig_x - h;
            double ymh = f(x)[func_index];

            x[deriv_index] = orig_x + 2 * h;
            double y2h = f(x)[func_index];

            x[deriv_index] = orig_x - 2 * h;
            double ym2h = f(x)[func_index];

            x[deriv_index] = orig_x + 3 * h;
            double y3h = f(x)[func_index];

            x[deriv_index] = orig_x - 3 * h;
            double ym3h = f(x)[func_index];

            double y1 = yh - ymh;
            double y2 = ym2h - y2h;
            double y3 = y3h - ym3h;

            if (error)
            {
                x[deriv_index] = orig_x + 4 * h;
                double y4h = f(x)[func_index];

                x[deriv_index] = orig_x - 4 * h;
                double ym4h = f(x)[func_index];

                double y7 = (y4h - ym4h - 6 * y3 - 14 * y1 - 14 * y2) / 2;
                *error = std::abs(y7) / (140 * h) + 5 * (std::abs(yh) + std::abs(ymh))*eps / h;
            }
            return (y3 + 9 * y2 + 45 * y1) / (60 * h);
        }

        template <int N>
        static VectorN<Real, N> NDer6PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(eps / 168, (double)1 / (double)7);

            return NDer6PartialByAll(f, func_index, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer6PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, double h, VectorN<Real, N> *error = nullptr)
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
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(3 * eps, static_cast<double>(1) / static_cast<double>(3));

            return NDer6PartialAllByAll(f, point, h, error);
        }

        template <int N>
        static MatrixNM<Real, N,N> NDer6PartialAllByAll(const IVectorFunction<N> &f, const VectorN<Real, N> &point, double h, MatrixNM<Real, N,N> *error = nullptr)
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

        /********************************************************************************************************************/
        /********                               Numerical derivatives of EIGHTH order                                ********/
        /********************************************************************************************************************/
        static double NDer8(const MML::IRealFunction &f, double x, double* error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            // Error bound ~eps^8/9.
            // In double precision, we only expect to lose two digits of precision while using this formula, at the cost of 8 function evaluations.
            // Error: h^8|f^(9)(x)|/630 + 7|f(x)|eps/h assuming 7 unstabilized additions.
            // Mathematica code to get the error:
            // Series[(f[x+h]-f[x-h])*(4/5) + (1/5)*(f[x-2*h] - f[x+2*h]) + (4/105)*(f[x+3*h] - f[x-3*h]) + (1/280)*(f[x-4*h] - f[x+4*h]), {h, 0, 9}]
            // If we used Kahan summation, we could get the max error down to h^8|f^(9)(x)|/630 + |f(x)|eps/h.
            double h = std::pow(551.25*eps, (double)1 / (double)9);
            //h = detail::make_xph_representable(x, h);

            double yh = f(x + h);
            double ymh = f(x - h);
            double y1 = yh - ymh;
            double y2 = f(x - 2 * h) - f(x + 2 * h);
            double y3 = f(x + 3 * h) - f(x - 3 * h);
            double y4 = f(x - 4 * h) - f(x + 4 * h);

            double tmp1 = 3 * y4 / 8 + 4 * y3;
            double tmp2 = 21 * y2 + 84 * y1;

            if (error)
            {
                // Mathematica code to generate fd scheme for 7th derivative:
                // Sum[(-1)^i*Binomial[9, i]*(f[x+(4-i)*h] + f[x+(5-i)*h])/2, {i, 0, 9}]
                // Mathematica to demonstrate that this is a finite difference formula for 7th derivative:
                // Series[(f[x+5*h]-f[x- 5*h])/2 + 4*(f[x-4*h] - f[x+4*h]) + 27*(f[x+3*h] - f[x-3*h])/2 + 24*(f[x-2*h]  - f[x+2*h]) + 21*(f[x+h] - f[x-h]), {h, 0, 15}]
                double f9 = (f(x + 5 * h) - f(x - 5 * h)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;
                *error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh))*eps / h;
            }
            return (tmp1 + tmp2) / (105 * h);
        }

        static double NDer8(const MML::IRealFunction &f, double x, double h, double* error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();

            double yh = f(x + h);
            double ymh = f(x - h);
            double y1 = yh - ymh;
            double y2 = f(x - 2 * h) - f(x + 2 * h);
            double y3 = f(x + 3 * h) - f(x - 3 * h);
            double y4 = f(x - 4 * h) - f(x + 4 * h);

            double tmp1 = 3 * y4 / 8 + 4 * y3;
            double tmp2 = 21 * y2 + 84 * y1;

            if (error)
            {
                double f9 = (f(x + 5 * h) - f(x - 5 * h)) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;
                *error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh))*eps / h;
            }
            return (tmp1 + tmp2) / (105 * h);
        }       

        template <int N>
        static double NDer8Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, double *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(551.25*eps, (double)1 / (double)9);

            return NDer8Partial(f, deriv_index, point, h, error);
        }

        template <int N>
        static double NDer8Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, double h, double *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();

            double     orig_x = point[deriv_index];

            VectorN<Real, N> x{point};

            x[deriv_index] = orig_x + h;
            double yh = f(x);

            x[deriv_index] = orig_x - h;
            double ymh = f(x);

            x[deriv_index] = orig_x + 2 * h;
            double y2h = f(x);

            x[deriv_index] = orig_x - 2 * h;
            double ym2h = f(x);

            x[deriv_index] = orig_x + 3 * h;
            double y3h = f(x);

            x[deriv_index] = orig_x - 3 * h;
            double ym3h = f(x);

            x[deriv_index] = orig_x + 4 * h;
            double y4h = f(x);

            x[deriv_index] = orig_x - 4 * h;
            double ym4h = f(x);

            double y1 = yh - ymh;
            double y2 = ym2h - y2h;
            double y3 = y3h - ym3h;
            double y4 = ym4h - y4h;

            double tmp1 = 3 * y4 / 8 + 4 * y3;
            double tmp2 = 21 * y2 + 84 * y1;

            if (error)
            {
                x[deriv_index] = orig_x + 5 * h;
                double y5h = f(x);

                x[deriv_index] = orig_x - 5 * h;
                double ym5h = f(x);

                double f9 = (y5h - ym5h) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;
                *error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh))*eps / h;

            }

            return (tmp1 + tmp2) / (105 * h);            
        }

        template <int N>
        static VectorN<Real, N> NDer8PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(551.25*eps, (double)1 / (double)9);

            return NDer8PartialByAll(f, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer8PartialByAll(const IScalarFunction<N> &f, const VectorN<Real, N> &point, double h, VectorN<Real, N> *error = nullptr)
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

        template <int N>
        static double NDer8Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, double *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(551.25*eps, (double)1 / (double)9);

            return NDer8Partial(f, func_index, deriv_index, point, h, error);
        }

        template <int N>
        static double NDer8Partial(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, double h, double *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();

            double     orig_x = point[deriv_index];

            VectorN<Real, N> x{point};

            x[deriv_index] = orig_x + h;
            double yh = f(x)[func_index];

            x[deriv_index] = orig_x - h;
            double ymh = f(x)[func_index];

            x[deriv_index] = orig_x + 2 * h;
            double y2h = f(x)[func_index];

            x[deriv_index] = orig_x - 2 * h;
            double ym2h = f(x)[func_index];

            x[deriv_index] = orig_x + 3 * h;
            double y3h = f(x)[func_index];

            x[deriv_index] = orig_x - 3 * h;
            double ym3h = f(x)[func_index];

            x[deriv_index] = orig_x + 4 * h;
            double y4h = f(x)[func_index];

            x[deriv_index] = orig_x - 4 * h;
            double ym4h = f(x)[func_index];

            double y1 = yh - ymh;
            double y2 = ym2h - y2h;
            double y3 = y3h - ym3h;
            double y4 = ym4h - y4h;

            double tmp1 = 3 * y4 / 8 + 4 * y3;
            double tmp2 = 21 * y2 + 84 * y1;

            if (error)
            {
                x[deriv_index] = orig_x + 5 * h;
                double y5h = f(x)[func_index];

                x[deriv_index] = orig_x - 5 * h;
                double ym5h = f(x)[func_index];

                double f9 = (y5h - ym5h) / 2 + 4 * y4 + 27 * y3 / 2 + 24 * y2 + 21 * y1;
                *error = std::abs(f9) / (630 * h) + 7 * (std::abs(yh) + std::abs(ymh))*eps / h;
            }

            return (tmp1 + tmp2) / (105 * h);            
        }        

        template <int N>
        static VectorN<Real, N> NDer8PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, VectorN<Real, N> *error = nullptr)
        {
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(551.25*eps, (double)1 / (double)9);

            return NDer8PartialByAll(f, func_index, point, h, error);
        }

        template <int N>
        static VectorN<Real, N> NDer8PartialByAll(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, double h, VectorN<Real, N> *error = nullptr)
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
            const double eps = (std::numeric_limits<double>::epsilon)();
            double h = std::pow(3 * eps, static_cast<double>(1) / static_cast<double>(3));

            return NDer8PartialAllByAll(f, point, h, error);
        }

        template <int N>
        static MatrixNM<Real, N,N> NDer8PartialAllByAll(const IVectorFunction<N> &f, const VectorN<Real, N> &point, double h, MatrixNM<Real, N,N> *error = nullptr)
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

        /********************************************************************************************************************/
        /********                            Definitions of default derivation functions                             ********/
        /********************************************************************************************************************/
        static inline double(*Derive)(const MML::IRealFunction &f, double x, double* error) = Derivation::NDer4;
        
        template<int N>
        static inline double(*DerivePartial)(const IScalarFunction<N> &f, int deriv_index, const VectorN<Real, N> &point, double *error) = Derivation::NDer4Partial;

        template<int N>
        static inline VectorN<Real, N>(*DerivePartialAll)(const IScalarFunction<N> &f, const VectorN<Real, N> &point, VectorN<Real, N> *error) = Derivation::NDer4PartialByAll;

        template<int N>
        static inline double(*DeriveVecPartial)(const IVectorFunction<N> &f, int func_index, int deriv_index, const VectorN<Real, N> &point, double *error) = Derivation::NDer4Partial;

        template<int N>
        static inline VectorN<Real, N>(*DeriveVecPartialAll)(const IVectorFunction<N> &f, int func_index, const VectorN<Real, N> &point, VectorN<Real, N> *error) = Derivation::NDer4PartialByAll;

        template<int N>
        static inline MatrixNM<Real, N,N>(*DeriveVecPartialAllByAll)(const IVectorFunction<N> &f, const VectorN<Real, N> &point, MatrixNM<Real, N,N> *error) = Derivation::NDer4PartialAllByAll;
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
                    ret[j] += Derivation::NDer1Partial(coordTransfFunc(j), k, pos, 1e-8) * vec[k];
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
                    ret[k] += Derivation::NDer1Partial(inverseCoordTransfFunc(k), j, pos, 1e-8) * vec[j];
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
                    ret[j] += Derivation::NDer1Partial(inverseCoordTransfFunc(k), j, pos, 1e-8) * vec[k];
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
                    ret[k] += Derivation::NDer1Partial(coordTransfFunc(j), k, pos, 1e-8) * vec[j];
                }
            }

            return ret;
        }        

        // transform tensor
    };    
  
}

///////////////////////////   ./include/basic_types/CoordTransf.h   ///////////////////////////



namespace MML
{
    // ovdje dodati translational, rotational, galilean, lorentzian transf
    // SVE su to transformacije koordinata
    class CoordTransfRectilinear : public ICoordTransf<VectorN<Real, 3>, VectorN<Real, 3>, 3>
    {
        public:
        Vector3Cartesian _base[3];
        VectorN<Real, 3> _dual[3];

        MatrixNM<Real,3,3> _alpha;
        MatrixNM<Real,3,3> _transf;

        ScalarFunctionFromStdFunc<3> _f1;
        ScalarFunctionFromStdFunc<3> _f2;
        ScalarFunctionFromStdFunc<3> _f3;

        ScalarFunctionFromStdFunc<3> _fInverse1;
        ScalarFunctionFromStdFunc<3> _fInverse2;
        ScalarFunctionFromStdFunc<3> _fInverse3;

        public:
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

        VectorN<Real, 3> Dual(int i)
        {
            return _dual[i];
        }
        VectorN<Real, 3>           transf(const VectorN<Real, 3> &q)           { return VectorN<Real, 3>{ func1(q), func2(q), func3(q) }; }
        VectorN<Real, 3>           transfInverse(const VectorN<Real, 3> &q)    { return VectorN<Real, 3>{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }
        
        IScalarFunction<3>&  coordTransfFunc(int i)         
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

    class CoordTransfSphericalToCartesian : public ICoordTransf<Vector3Spherical, Vector3Cartesian, 3>
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

        inline static ScalarFunctionFromFuncPtr<3> _func[3] = { ScalarFunctionFromFuncPtr<3>{func1},
                                                                ScalarFunctionFromFuncPtr<3>{func2},
                                                                ScalarFunctionFromFuncPtr<3>{func3}
                                                              };

        inline static ScalarFunctionFromFuncPtr<3> _funcInverse[3] = { ScalarFunctionFromFuncPtr<3>{funcInverse1},
                                                                       ScalarFunctionFromFuncPtr<3>{funcInverse2},
                                                                       ScalarFunctionFromFuncPtr<3>{funcInverse3}
                                                                     };

        Vector3Cartesian     transf(const Vector3Spherical &q)           { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
        Vector3Spherical     transfInverse(const Vector3Cartesian &q)    { return Vector3Spherical{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }
        
        IScalarFunction<3>&  coordTransfFunc(int i)         { return _func[i]; }
        IScalarFunction<3>&  inverseCoordTransfFunc(int i)  { return _funcInverse[i]; }
    };

    class CoordTransfCartesianToSpherical : public ICoordTransf<Vector3Cartesian, Vector3Spherical, 3>
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

        inline static ScalarFunctionFromFuncPtr<3> _func[3] = { 
                                                    ScalarFunctionFromFuncPtr<3>{func1},
                                                    ScalarFunctionFromFuncPtr<3>{func2},
                                                    ScalarFunctionFromFuncPtr<3>{func3}
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

    class CoordTransfCylindricalToCartesian : public ICoordTransf<Vector3Cylindrical, Vector3Cartesian, 3>
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

        inline static ScalarFunctionFromFuncPtr<3> _func[3] = { ScalarFunctionFromFuncPtr<3>{func1},
                                                                ScalarFunctionFromFuncPtr<3>{func2},
                                                                ScalarFunctionFromFuncPtr<3>{func3}
                                                              };

        inline static ScalarFunctionFromFuncPtr<3> _funcInverse[3] = { ScalarFunctionFromFuncPtr<3>{funcInverse1},
                                                                       ScalarFunctionFromFuncPtr<3>{funcInverse2},
                                                                       ScalarFunctionFromFuncPtr<3>{funcInverse3}
                                                                     };

        Vector3Cartesian     transf(const Vector3Cylindrical &q)         { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
        Vector3Cylindrical   transfInverse(const Vector3Cartesian &q)    { return Vector3Cylindrical{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }
        
        IScalarFunction<3>&  coordTransfFunc(int i)         { return _func[i]; }
        IScalarFunction<3>&  inverseCoordTransfFunc(int i)  { return _funcInverse[i]; }
    };    

    class CoordTransfCartesianToCylindrical : public ICoordTransf<Vector3Cartesian, Vector3Cylindrical, 3>
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

        inline static ScalarFunctionFromFuncPtr<3> _func[3] = { ScalarFunctionFromFuncPtr<3>{func1},
                                                                ScalarFunctionFromFuncPtr<3>{func2},
                                                                ScalarFunctionFromFuncPtr<3>{func3}
                                                              };

        inline static ScalarFunctionFromFuncPtr<3> _funcInverse[3] = { ScalarFunctionFromFuncPtr<3>{funcInverse1},
                                                                       ScalarFunctionFromFuncPtr<3>{funcInverse2},
                                                                       ScalarFunctionFromFuncPtr<3>{funcInverse3}
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

///////////////////////////   ./include/basic_types/Tensor.h   ///////////////////////////


namespace MML
{
    // imamo TensorData - to su koeficijenti
    // imamo ponašanje
    template <int N>
    class TensorRank2
    {
        private:
            double _coeff[N][N];
            int _numContravar;
            int _numCovar;

            double ScalarProduct(VectorN<Real, N> a, VectorN<Real, N> b)
            {
                return 0.0;
            }

            // TensorRank2 Transform(ICoordTransfTest<N> &a)
            // {
            //     return TensorRank2{};
            // }
    };

    template <int N>
    class TensorRank3
    {
        private:
            double _coeff[N][N][N];
    };

    template <int N>
    class TensorRank4
    {
        private:
            double _coeff[N][N][N][N];
    };
}

///////////////////////////   ./include/basic_types/MetricTensor.h   ///////////////////////////

namespace MML
{
    
template<int N>
class MetricTensor
{
    public:
    virtual double Component(int i, int j, const VectorN<Real, N> &pos) const = 0;
    virtual MatrixNM<Real, N, N> MetricAtPoint(const VectorN<Real, N> &pos) const = 0;
};

template<int N>
class MetricTensorCartesian: public MetricTensor<N>
{
    public:
    virtual double Component(int i, int j, const VectorN<Real, N> &pos) const
    {
        if( i == j )
            return 1.0;
        else
            return 0.0;
    }

    virtual MatrixNM<Real, N, N> MetricAtPoint(const VectorN<Real, N> &pos) const
    {
        MatrixNM<Real, N, N> ret;

        for( int i=0; i<3; i++ )
            for( int j=0; j<N; j++ )
                ret(i,j) = Component(i,j, pos);

        return ret;
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

    virtual MatrixNM<Real, 3, 3> MetricAtPoint(const VectorN<Real, 3> &pos) const
    {
        MatrixNM<Real, 3, 3> ret;

        for( int i=0; i<3; i++ )
            for( int j=0; j<3; j++ )
                ret(i,j) = Component(i,j, pos);

        return ret;
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

    virtual MatrixNM<Real, 3, 3> MetricAtPoint(const VectorN<Real, 3> &pos) const
    {
        MatrixNM<Real, 3, 3> ret;

        for( int i=0; i<3; i++ )
            for( int j=0; j<3; j++ )
                ret(i,j) = Component(i,j, pos);

        return ret;
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

    virtual MatrixNM<Real, N, N> MetricAtPoint(const VectorN<Real, N> &pos) const 
    {
        MatrixNM<Real, N, N> ret;

        for( int i=0; i<3; i++ )
            for( int j=0; j<N; j++ )
                ret(i,j) = Component(i,j, pos);

        return ret;
    }
};

}

///////////////////////////   ./include/basic_types/Fields.h   ///////////////////////////


namespace MML
{
    static double InverseRadialFieldFuncCart(const VectorN<Real, 3> &x )   { return 1.0 / x.NormL2(); }

    static double InverseRadialFieldFuncSpher(const VectorN<Real, 3> &x )
    {
        return 1.0 / x[0];
    }

    static double InverseRadialFieldFuncCyl(const VectorN<Real, 3> &x )
    {
        double r = x.NormL2();

        return 10.0 / sqrt(x[0]*x[0] + x[2]*x[2]);
    }

    class InverseRadialFieldCart : public IScalarFunction<3>
    {
    protected:
        double _constant;
    public:
        InverseRadialFieldCart() : _constant(1.0) {}
        InverseRadialFieldCart(double constant) : _constant(constant) {}

        double operator()(const VectorN<double, 3> &x) const  { return -_constant / x.NormL2(); }
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

///////////////////////////   ./include/algorithms/DiffEqSolvers.h   ///////////////////////////



namespace MML
{
    class RKSolver
    {
    public:

        void rk4(Vector<Real> &y, Vector<Real> &dydx, const double x, const double h,
            Vector<Real> &yout, void derivs(const double, Vector<Real> &, Vector<Real> &))
        {
            int i;
            double xh,hh,h6;

            int n= (int) y.size();
            Vector<Real> dym(n),dyt(n),yt(n);
            hh=h*0.5;
            h6=h/6.0;
            xh=x+hh;
            for (i=0;i<n;i++) yt[i]=y[i]+hh*dydx[i];
            derivs(xh,yt,dyt);
            for (i=0;i<n;i++) yt[i]=y[i]+hh*dyt[i];
            derivs(xh,yt,dym);
            for (i=0;i<n;i++) {
                yt[i]=y[i]+h*dym[i];
                dym[i] += dyt[i];
            }
            derivs(x+h,yt,dyt);
            for (i=0;i<n;i++)
                yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
        }


        void rkdumb(Vector<Real> &vstart, const double x1, const double x2, Vector<Real> &xx, Matrix<Real> &y,
            void derivs(const double, Vector<Real> &, Vector<Real> &))
        {
            int i,k;
            double x,h;

            int nvar=y.RowNum();
            int nstep=y.ColNum()-1;
            Vector<Real> v(nvar),vout(nvar),dv(nvar);
            for (i=0;i<nvar;i++) {
                v[i]=vstart[i];
                y[i][0]=v[i];
            }
            xx[0]=x1;
            x=x1;
            h=(x2-x1)/nstep;
            for (k=0;k<nstep;k++) {
                derivs(x,v,dv);
                rk4(v,dv,x,h,vout,derivs);
                if (x+h == x)
                    throw("Step size too small in routine rkdumb");
                x += h;
                xx[k+1]=x;
                for (i=0;i<nvar;i++) {
                    v[i]=vout[i];
                    y[i][k+1]=v[i];
                }
            }
        }        

        void rkck(Vector<Real> &y, Vector<Real> &dydx, const double x,
            const double h, Vector<Real> &yout, Vector<Real> &yerr,
            void derivs(const double, Vector<Real> &, Vector<Real> &))
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
            derivs(x+a2*h,ytemp,ak2);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
            derivs(x+a3*h,ytemp,ak3);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
            derivs(x+a4*h,ytemp,ak4);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
            derivs(x+a5*h,ytemp,ak5);
            for (i=0;i<n;i++)
                ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
            derivs(x+a6*h,ytemp,ak6);
            for (i=0;i<n;i++)
                yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
            for (i=0;i<n;i++)
                yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
        }


        void rkqs(Vector<Real> &y, Vector<Real> &dydx, double &x, const double htry,
            const double eps, Vector<Real> &yscal, double &hdid, double &hnext,
            void derivs(const double, Vector<Real> &, Vector<Real> &))
        {
            const double SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4;
            int i;
            double errmax,h,htemp,xnew;

            int n= (int) y.size();
            h=htry;
            Vector<Real> yerr(n),ytemp(n);
            for (;;) {
                rkck(y,dydx,x,h,ytemp,yerr,derivs);
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


        double dxsav;
        int kmax,kount;
        Vector<Real> *xp_p;
        Matrix<Real> *yp_p;

        void odeint(Vector<Real> &ystart, const double x1, const double x2, const double eps,
            const double h1, const double hmin, int &nok, int &nbad,
            void derivs(const double, Vector<Real> &, Vector<Real> &),
            void rkqs(Vector<Real> &, Vector<Real> &, double &, const double, const double,
            Vector<Real> &, double &, double &, void (*)(const double, Vector<Real> &, Vector<Real> &)))
        {
            const int MAXSTP=10000;
            const double TINY=1.0e-30;
            int i,nstp;
            double xsav,x,hnext,hdid,h;

            int nvar= (int) ystart.size();
            Vector<Real> yscal(nvar),y(nvar),dydx(nvar);
            Vector<Real> &xp=*xp_p;
            Matrix<Real> &yp=*yp_p;
            x=x1;
            h=SIGN(h1,x2-x1);
            nok = nbad = kount = 0;
            for (i=0;i<nvar;i++) y[i]=ystart[i];
            if (kmax > 0) xsav=x-dxsav*2.0;
            for (nstp=0;nstp<MAXSTP;nstp++) {
                derivs(x,y,dydx);
                for (i=0;i<nvar;i++)
                    yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
                if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
                    for (i=0;i<nvar;i++) yp[i][kount]=y[i];
                    xp[kount++]=x;
                    xsav=x;
                }
                if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
                rkqs(y,dydx,x,h,eps,yscal,hdid,hnext,derivs);
                if (hdid == h) ++nok; else ++nbad;
                if ((x-x2)*(x2-x1) >= 0.0) {
                    for (i=0;i<nvar;i++) ystart[i]=y[i];
                    if (kmax != 0) {
                        for (i=0;i<nvar;i++) yp[i][kount]=y[i];
                        xp[kount++]=x;
                    }
                    return;
                }
                if (fabs(hnext) <= hmin) throw("Step size too small in odeint");
                h=hnext;
            }
            throw("Too many steps in routine odeint");
        }

    };
}

///////////////////////////   ./include/algorithms/Interpolators.h   ///////////////////////////


namespace MML
{
    static bool polint(Vector<Real> &xa, Vector<Real> &ya, const double x, double &y, double &dy)
    // Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns a value y, and
    // an error estimate dy. If P (x) is the polynomial of degree N − 1 such that P (xai) = yai; i =
    // 1; : : : ; n, then the returned value y = P (x).
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

    // TODO - unused function
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

    // TODO - unused function
    static bool ratint(Vector<Real> &xa, Vector<Real> &ya, const double x, double &y, double &dy)
    {
        // Given arrays xa[1..n] and ya[1..n], and given a value of x, this routine returns a value of
        // y and an accuracy estimate dy. The value returned is that of the diagonal rational function,
        // evaluated at x, which passes through the n points (xai; yai), i = 1:::n.

        const double TINY=1.0e-25;
        int m,i,ns=0;
        double w,t,hh,h,dd;

        int n=(int)xa.size();
        Vector<Real> c(n),d(n);
        hh=fabs(x-xa[0]);
        for (i=0;i<n;i++) {
            h=fabs(x-xa[i]);
            if (h == 0.0) {
                y=ya[i];
                dy=0.0;
                return true;
            } else if (h < hh) {
                ns=i;
                hh=h;
            }
            c[i]=ya[i];
            d[i]=ya[i]+TINY;
        }
        y=ya[ns--];
        for (m=1;m<n;m++) {
            for (i=0;i<n-m;i++) {
                w=c[i+1]-d[i];
                h=xa[i+m]-x;
                t=(xa[i]-x)*d[i]/h;
                dd=t-c[i+1];
                if (dd == 0.0) 
                    //nrerror("Error in routine ratint");
                    return false;
                dd=w/dd;
                d[i]=c[i+1]*dd;
                c[i]=t*dd;
            }
            y += (dy=(2*(ns+1) < (n-m) ? c[ns+1] : d[ns--]));
        }
        return true;
    }
}  // end namespace

///////////////////////////   ./include/algorithms/Integration.h   ///////////////////////////


namespace MML
{
	class Integration
	{
		public:

        static double TrapRefine(IRealFunction &func, const double a, const double b, const int n)
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

        static double IntegrateTrap(IRealFunction &func, const double a, const double b)
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
            // errors may start increasing, and the routine may never converge. A value 10−6
            // is just on the edge of trouble for most 32-bit machines; it is achievable when the
            // convergence is moderately rapid, but not otherwise.
            int j;
            double s,olds=0.0;

            for (j=0;j<DefaultParams::IntegrateTrapJMAX;j++) 
            {
                s=TrapRefine(func,a,b,j+1);

                if (j > 5)
                    if (std::abs(s-olds) < DefaultParams::IntegrateTrapEPS * std::abs(olds) ||
                        (s == 0.0 && olds == 0.0)) 
                        return s;
                olds=s;
            }
            throw IntegrationTooManySteps("qtrap");
        }

        static double IntegrateSimpson(IRealFunction &func, const double a, const double b)
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

            for (j=0;j<DefaultParams::IntegrateSimpJMAX;j++) 
            {
                st=TrapRefine(func,a,b,j+1);
                s=(4.0*st-ost)/3.0;
                
                if (j > 5)
                    if (std::abs(s-os) < DefaultParams::IntegrateSimpEPS * std::abs(os) ||
                        (s == 0.0 && os == 0.0)) 
                        return s;
                os=s;
                ost=st;
            }
            throw IntegrationTooManySteps("qsimp");
        }

        static double IntegrateRomberg(IRealFunction &func, double a, double b)
        {
            // Returns the integral of the function func from a to b. Integration is performed by Romberg’s
            // method of order 2K, where, e.g., K=2 is Simpson’s rule.

            // The routine qromb, along with its required trapzd and polint, is quite
            // powerful for sufficiently smooth (e.g., analytic) integrands, integrated over intervals
            // which contain no singularities, and where the endoubleoints are also nonsingular. qromb,
            // in such circumstances, takes many, many fewer function evaluations than either of
            // the routines in x4.2
            const int JMAXP=DefaultParams::IntegrateRombJMAX+1, K=5;
            double ss,dss;
            Vector<double> s(DefaultParams::IntegrateRombJMAX),h(JMAXP),s_t(K),h_t(K);

            h[0]=1.0;
            for (int j=1;j<=DefaultParams::IntegrateRombJMAX;j++) 
            {
                s[j-1]=TrapRefine(func,a,b,j);
            
                if (j >= K) 
                {
                    for (int i=0;i<K;i++) {
                        h_t[i]=h[j-K+i];
                        s_t[i]=s[j-K+i];
                    }
                
                    polint(h_t,s_t,0.0,ss,dss);
                
                    if (std::abs(dss) <= DefaultParams::IntegrateRombEPS * std::abs(ss)) 
                        return ss;
                }

                h[j]=0.25*h[j-1];
            }
            throw IntegrationTooManySteps("qromb");
        }
	};
} // end namespace
///////////////////////////   ./include/algorithms/PathIntegration.h   ///////////////////////////


namespace MML
{
    // LineIntegralScalar
    // LineIntegralVector
    // SurfaceIntegral - flux
	class PathIntegration
	{
		public:
        template<int N>
        static double PathIntegrateTrap(IScalarFunction<N> &func, IParametricCurve<N> &path, const double a, const double b)
        {
            return 0.0;
        }


	};
} // end namespace
///////////////////////////   ./include/algorithms/FieldOperations.h   ///////////////////////////// grad
// - cart
// - spher
// - cyl




namespace MML
{
    class ScalarFieldOperations
    {
        public:
                // Generalni gradijent
        template<int N>
        static VectorN<Real, N> Gradient(IScalarFunction<N> &scalarField, const MetricTensor<N>& metric, const VectorN<Real, N> &pos)
        {
            VectorN<Real, N> derivsAtPoint = Derivation::DerivePartialAll<N>(scalarField, pos, nullptr);

            VectorN<Real, N> ret = metric.MetricAtPoint(pos) * derivsAtPoint;

            return ret;
        }

        template<int N>
        static VectorN<Real, N> GradientCart(IScalarFunction<N> &scalarField, const VectorN<Real, N> &pos)
        {
            return Derivation::DerivePartialAll<N>(scalarField, pos, nullptr);

            // const MetricTensor<N>& m = Metric();

            // MatrixNM<N,N> metricAtPoint = m.MetricAtPoint(pos);
            // ret = m.MetricAtPoint(pos) * derivsAtPoint;
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
        template<int N>
        static double DivCart(const IVectorFunction<N> &vectorField, const VectorN<Real, N> &pos)
        {
            // VectorN<Real, N> v1 = Derivation::DeriveVecPartialAll<N>(vectorField, 0, pos, nullptr);
            // MatrixNM<N,N> v2 = Derivation::DeriveVecPartialAllByAll<N>(vectorField, pos, nullptr);

            // Derivation::DeriveVecPartial<N> = Derivation::NDer6Partial;

            double div = 0.0;
            for( int i=0; i<N; i++ )
            {
                div += Derivation::DeriveVecPartial<N>(vectorField, i, i, pos, nullptr);
            }

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
    };
}

///////////////////////////   ./include/algorithms/RootFinding.h   ///////////////////////////// bit će toga


namespace MML
{
    class RootFinding
    {

    };
}


///////////////////////////   ./include/algorithms/Statistics.h   ///////////////////////////// bit će toga


namespace MML
{
    class Statistics
    {

    };
}


///////////////////////////   ./include/systems/LinAlgEqSystem.h   ///////////////////////////

namespace MML
{
    class LinAlgEqSystem
    {
        // dva ctora - Matrix, Matrix i Matrix, Vector

        // solve by GJ

        // solve by LU - dvije verzije - in place, i verzija koja vraća dekompoziciju
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
