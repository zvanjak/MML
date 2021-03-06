#if !defined  __MML_BASIC_TYPES_SINGLE_HEADER
#define __MML_BASIC_TYPES_SINGLE_HEADER

#include <exception>

#include <string>
#include <vector>

#include <cmath>

#include <iostream>
#include <iomanip>

#include <functional>
#include <limits>

#include <initializer_list>

#define SIGN(a,b) ( (b) >= 0.0 ? fabs(a): -fabs(a) )

#define SQR(a) ( (a) * (a) )///////////////////////////   ./include/MMLStandardHeaders.h   ///////////////////////////






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
    };

    class ITabulatedValues2D
    {
    public:
        virtual int    NumValues1() = 0;
        virtual int    NumValues2() = 0;

        virtual double Value(int ind1, int in2) = 0;
        virtual double Value(int ind1, int in2, double &outX, double &outY) = 0;
    };

    class ITabulatedValues3D
    {
    public:
        virtual int    NumValues1() = 0;
        virtual int    NumValues2() = 0;
        virtual int    NumValues3() = 0;

        virtual double Value(int ind1, int in2, int ind3) = 0;
        virtual double Value(int ind1, int in2, int ind3, double &outX, double &outY, double &outZ) = 0;
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

    class TabulatedValues1D
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

///////////////////////////   ./include/basic_types/Vector.h   ///////////////////////////

namespace MML
{
    class Vector
    {
        private:
            std::vector<double> _elems;

        public:
            Vector(size_t n) : _elems(n) {}
            Vector(std::initializer_list<double> list) : _elems(list) { }

            double& operator[](int n) { return _elems[n]; }
            double  operator[](int n) const { return _elems[n]; }
            
            size_t size() const { return _elems.size(); }

            double ScalarProductCartesian(Vector &b)
            {
                double product = 0.0;
                for( int i=0; i<size(); i++ )
                    product += (*this)[i] * b[i];
                return product;
            }

            double NormCartesian()
            {
                double norm = 0.0;
                for( int i=0; i<size(); i++ )
                    norm += (*this)[i] * (*this)[i];
                return norm;
            }

            bool IsEqual(const Vector &b, double eps)
            {
                for( int i=0; i<size(); i++ )
                {
                    if( fabs((*this)[i] - b[i]) > eps )
                        return false;
                }
                return true;
            }

            std::string to_string(std::string format)
            {
                std::string a{"0.0"};

                return a;
            }

            Vector operator+(Vector &b )
            {
                Vector ret(b.size());;
                for(int i=0; i<b.size(); i++)
                    ret._elems[i] = (*this)[i] + b._elems[i];
                return ret;
            }

            Vector operator-(Vector &b )
            {
                Vector ret(b.size());;
                for(int i=0; i<b.size(); i++)
                    ret._elems[i] = (*this)[i] - b._elems[i];
                return ret;
            }

            friend Vector operator*(double a, Vector &b )
            {
                Vector ret(b.size());;
                for(int i=0; i<b.size(); i++)
                    ret._elems[i] = a * b._elems[i];
                return ret;
            }

            friend Vector operator*(Vector &a, double b )
            {
                Vector ret(a.size());;
                for(int i=0; i<a.size(); i++)
                    ret._elems[i] = b * a._elems[i];
                return ret;
            }

            friend Vector operator/(Vector &a, double b)
            {
                Vector ret(a.size());
                for(int i=0; i<a.size(); i++)
                    ret._elems[i] = a._elems[i] / b;
                return ret;
            }

            friend std::ostream& operator<<(std::ostream& stream, const Vector &a)
            {
                stream << "[";
                bool first = true;
                for(const double& x : a._elems)
                {
                    if( !first )
                        stream << ", ";
                    else
                        first = false;

                    stream << std::setprecision(15) << x;
                }
                stream << "]";

                return stream;
            }
    };
}

///////////////////////////   ./include/basic_types/VectorN.h   ///////////////////////////

namespace MML
{
    template<int N> 
    class VectorN
    {
        protected:
            double  _val[N];

        public:
        VectorN() {}

        VectorN(std::initializer_list<double> list) 
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

        double& operator[](int n)       { return _val[n]; }
        double  operator[](int n) const { return _val[n]; }

        int size() const { return N; }

        double NormL2()
        {
            double norm = 0.0;
            for( int i=0; i<size(); i++ )
                norm += (*this)[i] * (*this)[i];
            return norm;
        }

        bool IsEqual(VectorN &b, double eps)
        {
            return false;
        }

        VectorN operator+(VectorN &b )
        {
            VectorN ret;
            for(int i=0; i<N; i++)
                ret._val[i] = _val[i] + b._val[i];
            return ret;
        }

        VectorN operator-(VectorN &b )
        {
            VectorN ret;
            for(int i=0; i<N; i++)
                ret._val[i] = _val[i] - b._val[i];
            return ret;
        }

        friend VectorN operator*(VectorN &a, double b )
        {
            VectorN ret;
            for(int i=0; i<N; i++)
                ret._val[i] = a[i] * b;
            return ret;
        }

        friend VectorN operator*(double a, VectorN &b )
        {
            VectorN ret;
            for(int i=0; i<N; i++)
                ret._val[i] = a * b[i];
            return ret;
        }

        friend VectorN operator/(VectorN &a, double b)
        {
            VectorN ret;
            for(int i=0; i<N; i++)
                ret._val[i] = a[i] / b;
            return ret;
        }

        friend std::ostream& operator<<(std::ostream& stream, const VectorN<N> &a)
        {
            stream << "[";
            bool first = true;
            for(const double& x : a._val)
            {
                if( !first )
                    stream << ", ";
                else
                    first = false;

                stream << std::setw(8) << x;
            }
            stream << "]";

            return stream;
        }
    };

    class Vector2Cartesian : public VectorN<2>
    {
        public:
            double  X() const   { return _val[0]; }
            double& X()         { return _val[0]; }
            double  Y() const   { return _val[1]; }
            double& Y()         { return _val[1]; }
    };

    class Vector2Polar : public VectorN<2>
    {
        public:
            double  R() const   { return _val[0]; }
            double& R()         { return _val[0]; }
            double  Phi() const { return _val[1]; }
            double& Phi()       { return _val[1]; }
    };

    class Vector3Cartesian : public VectorN<3>
    {
        public:
            double  X() const   { return _val[0]; }
            double& X()         { return _val[0]; }
            double  Y() const   { return _val[1]; }
            double& Y()         { return _val[1]; }
            double  Z() const   { return _val[2]; }
            double& Z()         { return _val[2]; }

            bool IsParallelTo(Vector3Cartesian &b, double eps=1e-15)
            {
                double norm1 = NormL2();
                double norm2 = b.NormL2();

                return  fabs(X()/norm1 - b.X()/norm2) < eps &&
                        fabs(Y()/norm1 - b.Y()/norm2) < eps &&
                        fabs(Z()/norm1 - b.Z()/norm2) < eps;
            }
            bool IsPerpendicularTo(Vector3Cartesian &b, double eps=1e-15)
            {
                if( fabs(ScalarProd(*this,b)) < eps)
                    return true;
                else
                    return false;
            }

            Vector3Cartesian GetUnitVector()
            {
                VectorN<3> res = (*this) / NormL2();

                return Vector3Cartesian{res};
            }

            friend double ScalarProd(Vector3Cartesian &a, Vector3Cartesian &b)
            {
                return a.X() * b.X() + a.Y() * b.Y() + a.Z() * b.Z();
            }

            friend Vector3Cartesian VectorProd(Vector3Cartesian &a, Vector3Cartesian &b)
            {
                Vector3Cartesian ret;

                ret.X() = a.Y() * b.Z() - a.Z() * b.Y();
                ret.Y() = a.Z() * b.X() - a.X() * b.Z();
                ret.Z() = a.X() * b.Y() - a.Y() * b.X();

                return ret;
            }

            friend double VectorsAngle(Vector3Cartesian &a, Vector3Cartesian &b)
            {
                double cos_phi = ScalarProd(a,b) / (a.NormL2() * b.NormL2());

                return acos(cos_phi);
            }            
    };

    class Vector3Spherical : public VectorN<3>
    {
        public:
            double  R() const       { return _val[0]; }
            double& R()             { return _val[0]; }
            double  Theta() const   { return _val[1]; }
            double& Theta()         { return _val[1]; }
            double  Phi() const     { return _val[2]; }
            double& Phi()           { return _val[2]; }
    };    
}



typedef MML::VectorN<2> Vector2;
typedef MML::VectorN<3> Vector3;
typedef MML::VectorN<4> Vector4;

// SPECIJALIZACIJE ZA 2, 3 I 4
// template<>
// class VectorN<3>
// {

// };

///////////////////////////   ./include/basic_types/Matrix.h   ///////////////////////////


namespace MML
{
    class Matrix
    {
    public:
        Matrix() : _rows(0), _cols(0)
        {
            _ptrData = nullptr;
            //std::cout << "Matrix::default constructor \n";
        }
        Matrix(size_t rows, size_t cols) : _rows(rows), _cols(cols)
        {
            _ptrData = new double *[_rows];
            for (size_t i = 0; i < _rows; ++i)
                _ptrData[i] = new double[_cols];
            for (size_t i = 0; i < _rows; ++i)
                for (size_t j = 0; j < _cols; ++j)
                    _ptrData[i][j] = 0;
            
            //std::cout << "Matrix::constructor \n";
        }
        Matrix(size_t rows, size_t cols, std::initializer_list<double> values) : _rows(rows), _cols(cols)
        {
            _ptrData = new double *[_rows];
            for (size_t i = 0; i < _rows; ++i)
                _ptrData[i] = new double[_cols];
            
            size_t cnt = 0;
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
            
            //std::cout << "Matrix::constructor \n";
        }
        Matrix(const Matrix &m) : _rows(m._rows), _cols(m._cols)
        {
            _ptrData = new double *[_rows];
            for (size_t i = 0; i < _rows; ++i)
                _ptrData[i] = new double[_cols];
            for (size_t i = 0; i < _rows; ++i)
                for (size_t j = 0; j < _cols; ++j)
                    _ptrData[i][j] = m._ptrData[i][j];

            //std::cout << "Matrix::copy constructor \n";
        }
        Matrix(Matrix &&m)
        {
            _ptrData = m._ptrData;

            _rows = m._rows;
            _cols = m._cols;

            m._rows = 0;
            m._cols = 0;
            m._ptrData = nullptr;

            //std::cout << "Matrix::move constructor \n";
        }
        ~Matrix()
        {
            for (size_t i = 0; i < _rows; ++i)
                if (_ptrData != nullptr && _ptrData[i] != nullptr)
                    delete[] _ptrData[i];
            if (_ptrData != nullptr)
                delete[] _ptrData;
            //std::cout << "Matrix::destructor \n";
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
        }

        static Matrix RowMatrixFromVector(Vector &b)
        {
            Matrix ret(1, b.size());
            for( int i=0; i<b.size(); i++)
                ret[0][i] = b[i];

            return ret;
        }

        static Matrix ColumnMatrixFromVector(Vector &b)
        {
            Matrix ret(b.size(), 1);
            for( int i=0; i<b.size(); i++)
                ret[i][0] = b[i];

            return ret;
        }

        static Vector VectorFromRow(Matrix &a, int rowInd)
        {
            Vector ret(a.ColNum());
            for( int i=0; i<a.ColNum(); i++)
                ret[i] = a[rowInd][i];

            return ret;
        }

        static Vector VectorFromColumn(Matrix &a, int colInd)
        {
            Vector ret(a.RowNum());
            for( int i=0; i<a.RowNum(); i++)
                ret[i] = a[i][colInd];

            return ret;
        }
        
        bool IsEqual(Matrix &b, double eps)
        {
            return false;
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
                _ptrData = new double *[_rows];
                for (size_t i = 0; i < _rows; ++i)
                    _ptrData[i] = new double[_cols];
            }

            for (size_t i = 0; i < _rows; ++i)
                for (size_t j = 0; j < _cols; ++j)
                    _ptrData[i][j] = m._ptrData[i][j];

            //std::cout << "Matrix::operator = \n";
            return *this;
        }
        Matrix &operator=(Matrix &&m)
        {
            if (this == &m)
                return *this;

            std::swap(_ptrData, m._ptrData);
            std::swap(_rows, m._rows);
            std::swap(_cols, m._cols);

            //std::cout << "Matrix::move operator = \n";
            return *this;
        }

        double *operator[](int i)
        {
            return _ptrData[i];
        }

        Matrix operator+(const Matrix &other)
        {
            //std::cout << "Matrix::operator + START \n";

            Matrix temp(_rows, _cols);
            if (_rows != other._rows || _cols != other._cols)
            {
                for (size_t i = 0; i < _rows; i++)
                    for (size_t j = 0; j < _cols; j++)
                        temp._ptrData[i][j] = _ptrData[i][j];
                return temp;
            }
            else
            {
                for (size_t i = 0; i < _rows; i++)
                    for (size_t j = 0; j < _cols; j++)
                        temp._ptrData[i][j] = other._ptrData[i][j] + _ptrData[i][j];
            }
            //std::cout << "Matrix::operator + END \n";
            return temp;
        }

        Matrix operator-(const Matrix &other)
        {
            //std::cout << "Matrix::operator + START \n";

            Matrix temp(_rows, _cols);
            if (_rows != other._rows || _cols != other._cols)
            {
                for (int i = 0; i < _rows; i++)
                    for (int j = 0; j < _cols; j++)
                        temp._ptrData[i][j] = _ptrData[i][j];
                return temp;
            }
            else
            {
                for (int i = 0; i < _rows; i++)
                    for (int j = 0; j < _cols; j++)
                        temp._ptrData[i][j] =  _ptrData[i][j] - other._ptrData[i][j];
            }
            //std::cout << "Matrix::operator + END \n";
            return temp;
        }

        Matrix  operator*( const Matrix &b ) const
        {
            Matrix	ret(RowNum(), b.ColNum());

            if( ColNum() == b.RowNum() )
            {
                for( int i=0; i<ret.RowNum(); i++ )
                    for( int j=0; j<ret.ColNum(); j++ )
                    {
                        ret._ptrData[i][j] = 0;
                        for(int k=0; k<ColNum(); k++ )
                            ret._ptrData[i][j] += _ptrData[i][k] * b._ptrData[k][j];
                    }
            }
            else
            {				// krive dimenzije matrice
            }

            return	ret;
        }

       friend Matrix operator*( Matrix &a, double b )
        {
            int	i, j;
            Matrix	ret(a.RowNum(), a.ColNum());

            for( i=0; i<a.RowNum(); i++ )
                for( j=0; j<a.ColNum(); j++ )
                    ret[i][j] = a._ptrData[i][j] * b;

            return ret;
        }

        friend Matrix operator*( double a, Matrix &b )
        {
            int	i, j;
            Matrix	ret(b.RowNum(), b.ColNum());

            for( i=0; i<b.RowNum(); i++ )
                for( j=0; j<b.ColNum(); j++ )
                    ret[i][j] = a * b._ptrData[i][j];

            return ret;
        }

        friend Vector operator*( Matrix &a, Vector &b )
        {
            int	i, j;
            Vector	ret(a.RowNum());

            if( a.ColNum() == b.size() )
            {
                for( i=0; i<a.RowNum(); i++ )
                {
                    ret[i] = 0;
                    for( j=0; j<a.ColNum(); j++ )
                        ret[i] += a._ptrData[i][j] * b[j];
                }
            }
            else
            {
                std::cout << "Krive dimenzije kod mnozenja matrice i vektora !!!";
            }

            return ret;
        }

        friend Vector operator*( Vector &a, Matrix &b )
        {
            int	i, j;
            Vector	ret(b.ColNum());

            if( a.size() == b.RowNum() )
            {
                for( i=0; i<b.ColNum(); i++ )
                {
                    ret[i] = 0;
                    for( j=0; j<b.RowNum(); j++ )
                        ret[i] += a[i] * b[i][j];
                }
            }
            else
            {
                std::cout << "Krive dimenzije kod mnozenja matrice i vektora !!!";
            }

            return ret;
        }

        void Print()
        {
            std::cout << "Rows: " << _rows << "  Cols: " << _cols << std::endl;

            for (size_t i = 0; i < _rows; i++)
            {
                for (size_t j = 0; j < _cols; j++)
                {
                    std::cout << std::setw(10) << std::setprecision(3) << _ptrData[i][j] << ", ";
                }
                std::cout << std::endl;
            }
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
            _ptrData = new double *[_rows];
            for (size_t i = 0; i < _rows; ++i)
                _ptrData[i] = new double[_cols];
            for (size_t i = 0; i < _rows; ++i)
                for (size_t j = 0; j < _cols; ++j)
                    _ptrData[i][j] = 0;
            
            //std::cout << "Matrix::constructor \n";
        }

        friend std::ostream& operator<<(std::ostream& stream, const Matrix &a)
        {
            stream << "Rows: " << a.RowNum() << "  Cols: " << a.ColNum() << std::endl;

            for (size_t i = 0; i < a.RowNum(); i++)
            {
                stream << "[ ";
                for (size_t j = 0; j < a.ColNum(); j++)
                {
                    stream << std::setw(10) << std::setprecision(3) << a._ptrData[i][j] << ", ";
                }                
                stream << " ]" << std::endl;
            }

            return stream;
        }        

    private:
        size_t _rows;
        size_t _cols;
        double **_ptrData;
    };
}

///////////////////////////   ./include/basic_types/MatrixNM.h   ///////////////////////////


namespace MML
{

    template <int N, int M>
    class MatrixNM
    {
    public:
        double _vals[N][M];

    public:
        MatrixNM() 
        {
            for (size_t i = 0; i < N; ++i)
                for (size_t j = 0; j < M; ++j)
                    _vals[i][j] = 0.0;
        }

        MatrixNM(std::initializer_list<double> values) 
        {
            size_t cnt = 0;
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

        MatrixNM &operator=(const MatrixNM &m)
        {
            if (this == &m)
                return *this;

            for (size_t i = 0; i < RowNum(); ++i)
                for (size_t j = 0; j < ColNum(); ++j)
                    _vals[i][j] = m._vals[i][j];

            return *this;
        }

        double &operator()(int i, int j) { return _vals[i][j]; }

        MatrixNM operator+(const MatrixNM &other)
        {
            MatrixNM temp;
            for (size_t i = 0; i < RowNum(); i++)
                for (size_t j = 0; j < ColNum(); j++)
                    temp._vals[i][j] = other._vals[i][j] + _vals[i][j];
            return temp;
        }

        MatrixNM operator-(const MatrixNM &other)
        {
            MatrixNM temp;
            for (size_t i = 0; i < RowNum(); i++)
                for (size_t j = 0; j < ColNum(); j++)
                    temp._vals[i][j] = _vals[i][j] - other._vals[i][j];
            return temp;
        }

        template<int K>
        MatrixNM<N, K>  operator*( const MatrixNM<M, K> &b ) const
        {
            MatrixNM<N, K>	ret;

            for( int i=0; i<ret.RowNum(); i++ )
                for( int j=0; j<ret.ColNum(); j++ )
                {
                    ret._vals[i][j] = 0;
                    for(int k=0; k<ColNum(); k++ )
                        ret._vals[i][j] += _vals[i][k] * b._vals[k][j];
                }

            return	ret;
        }        


       friend MatrixNM operator*( MatrixNM &a, double b )
        {
            int	i, j;
            MatrixNM	ret(a.RowNum(), a.ColNum());

            for( i=0; i<a.RowNum(); i++ )
                for( j=0; j<a.ColNum(); j++ )
                    ret[i][j] = a._ptrData[i][j] * b;

            return ret;
        }

        friend MatrixNM operator*( double a, MatrixNM &b )
        {
            int	i, j;
            MatrixNM	ret(b.RowNum(), b.ColNum());

            for( i=0; i<b.RowNum(); i++ )
                for( j=0; j<b.ColNum(); j++ )
                    ret[i][j] = a * b._ptrData[i][j];

            return ret;
        }

        friend VectorN<N> operator*( MatrixNM &a, VectorN<M> &b )
        {
            int	i, j;
            VectorN<N>	ret;

            for( i=0; i<a.RowNum(); i++ )
            {
                ret[i] = 0;
                for( j=0; j<a.ColNum(); j++ )
                    ret[i] += a._vals[i][j] * b[j];
            }

            return ret;
        }

        friend VectorN<M> operator*( VectorN<N> &a, MatrixNM &b )
        {
            int	i, j;
            VectorN<M>	ret;

            for( i=0; i<b.ColNum(); i++ )
            {
                ret[i] = 0;
                for( j=0; j<b.RowNum(); j++ )
                    ret[i] += a[i] * b._vals[j][i];
            }

            return ret;
        }

        friend std::ostream& operator<<(std::ostream& stream, const MatrixNM &a)
        {
            stream << "Rows: " << a.RowNum() << "  Cols: " << a.ColNum() << std::endl;

            for (size_t i = 0; i < a.RowNum(); i++)
            {
                stream << "[ ";
                for (size_t j = 0; j < a.ColNum(); j++)
                {
                    stream << std::setw(10) << std::setprecision(3) << a._vals[i][j] << ", ";
                }                
                stream << " ]" << std::endl;
            }

            return stream;
        }   
        
        void Print()
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    std::cout << std::setw(10) << std::setprecision(3) << _vals[i][j] << ", ";
                }
                std::cout << std::endl;
            }
        }               
    };
}

///////////////////////////   ./include/basic_types/Polynom.h   ///////////////////////////
namespace MML
{
	class	Polynom
	{
	private:
		int			m_nDegree;
		double	*m_pdCoef;

	public:
		Polynom();
		Polynom( int n );
		Polynom( const Polynom &Copy );
		~Polynom();

		friend	int			GetDegree( const Polynom &a );
		friend	int			GetRealDegree( const Polynom & a );
		friend	double	GetCoef( const Polynom &a, int CoefNum );
		friend	double	GetLeadingCoef( const Polynom &a );
		friend	void		ReducePolynom( Polynom *a );

		double		Val( double x );
		Polynom		operator+( const Polynom &b ) const;
		Polynom		operator-( const Polynom &a ) const;
		Polynom		operator*( const Polynom &b ) const;
		Polynom		operator/( const Polynom &b ) const;
		Polynom		operator%( const Polynom &b ) const;	// dijeljenje po modulu (ostatak)
		Polynom&	operator=( const Polynom &b );

		Polynom		operator+=( const Polynom &b );
		Polynom		operator-=( const Polynom &b );
		Polynom		operator*=( const Polynom &b );
		Polynom		operator*=( const double b );
		Polynom		operator/=( const double &b );
	//	Polynom		operator/=( const Polynom &b );

		double&		operator[]( int i );

		friend	Polynom	operator*( double a, const Polynom &b );
		friend	Polynom	operator*( const Polynom &a, double b );

		friend	void	Input( Polynom *a );
		friend	void	Print( const Polynom &a );
	};
}
///////////////////////////   ./include/basic_types/MatrixSparse.h   ///////////////////////////
namespace MML
{

}

///////////////////////////   ./include/basic_types/Geometry.h   ///////////////////////////

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
    };
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

            Point3Cartesian operator+(VectorN<3> &b)
            {
                Point3Cartesian ret;
                ret.X() = X() + b[0];
                ret.Y() = Y() + b[1];
                ret.Z() = Z() + b[2];
                return ret;
            }
            Vector3Cartesian operator-(Point3Cartesian &b)
            {
                Vector3Cartesian ret;
                ret.X() = b.X() - X();
                ret.Y() = b.Y() - Y();
                ret.Z() = b.Z() - Z();
                return ret;
            }
    };    

    class Line2D
    {
    };

    class Line3D
    {
        public:
            Point3Cartesian _point;
            Vector3Cartesian _direction;        // da bude jedinicni vektor? normiranje u ctoru?

            Line3D(Point3Cartesian &pnt, Vector3Cartesian dir)
            {
                _point = pnt;
                _direction = dir.GetUnitVector();
            }

            Line3D(Point3Cartesian &a, Point3Cartesian &b)
            {
                Vector3Cartesian dir = b - a;
                _point = a;
                _direction = dir.GetUnitVector();
            }

            Point3Cartesian PointOnLine(double t)
            {
                VectorN<3> dist = t * _direction;
                Point3Cartesian ret = _point + dist;
                return ret;
            }

            // distance Line - Point3
            // nearest point on line
    };

    class SegmentLine2D
    {
        // reprezentacija - dvije Point2D

        // ctors
        //  dviej tocke
        //  tocka i koeficijne smjera
    };

    class SegmentLine3D
    {
        // reprezentacija - dvije Point3D

        // ctors
        //  dviej tocke
        //  tocka i koeficijne smjera
    };

    class Plane
    {
        public:
            double _A, _B, _C, _D;

            Plane(Point3Cartesian &a, Vector3Cartesian &normal)
            {

            }

            Plane(Point3Cartesian &a, Point3Cartesian &b, Point3Cartesian &C)
            {

            }

            Plane(double alpha, double beta, double gamma, double d)        // Hesseov (normalni) oblik
            {
                _A = cos(alpha);
                _B = cos(beta);
                _C = cos(gamma);
                _D = -d;
            }

            // tri segmenta na koord osima ctor

            Vector3Cartesian Normal()
            {
                Vector3Cartesian ret;

                return ret;
            }

            // GetCoordAxisSegments
            // GetHesseNormalFormParams

            // IsPointOnPlane

            // IsLineOnPlane

            // Angle with line

            // Intersection with line

            // 2 plane intersection
    };

    class Triangle2D
    {

    };

    class Polygon2D
    {

    };
    
    class Triangle3D
    {

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
    class IScalarFunction : public IFunction<double, VectorN<N> &>
    {
        public:
        virtual double operator()(VectorN<N> &x) const = 0;
    };

    //////////////////////////////////////////////////////////////////////
    template<int N>
    class IVectorFunction : public IFunction<VectorN<N>, VectorN<N> &>
    {
        public:
        virtual VectorN<N> operator()(VectorN<N> &x) const = 0;
    };
}

///////////////////////////   ./include/basic_types/Function.h   ///////////////////////////


namespace MML
{
    //////////////////////////////////////////////////////////////////////
    class RealFunction : public IRealFunction
    {
        double (*_func)(double) ;

        public:
        RealFunction(double (*inFunc)(double) ) : _func(inFunc)    {}

        double operator()(double x) const    { return _func(x); }
    };
    class RealFunctionFromStdFunc : public IRealFunction
    {
        std::function<double(double)> _func;

        public:
        RealFunctionFromStdFunc(std::function<double(double)> &inFunc) : _func(inFunc)    {}

        double operator()(double x) const    { return _func(x); }
    };

    //////////////////////////////////////////////////////////////////////
    template<int N>
    class ScalarFunctionFromFuncPtr : public IScalarFunction<N>
    {
        public:
        double (*_func)(VectorN<N> &);

        ScalarFunctionFromFuncPtr( double (*inFunc)(VectorN<N> &) ) : _func(inFunc)
        {}

        double operator()(VectorN<N> &x) const  { return _func(x); }
    };

    template<int N>
    class ScalarFunctionFromStdFunc : public IScalarFunction<N>
    {
        public:
        std::function<double(VectorN<N> &)> _func;

        ScalarFunctionFromStdFunc(std::function<double(VectorN<N> &)> inFunc) : _func(inFunc)
        {}

        double operator()(VectorN<N> &x) const  { return _func(x); }
    };
    
    //////////////////////////////////////////////////////////////////////
    template<int N>
    class VectorFunctionFromStdFunc : public IVectorFunction<N>
    {
        public:
        std::function<VectorN<N>(VectorN<N> &)> _func;

        VectorFunctionFromStdFunc(std::function<VectorN<N>(VectorN<N> &)> &inFunc) : _func(inFunc)    {}

        VectorN<N> operator()(VectorN<N> &x) const   { return _func(x); }
    };
} // end namespace

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

    template<int N>    
    class InterpolatedScalarFunction : public IScalarFunction<N>
    {
        public:
        InterpolatedScalarFunction() 
        {

        }

        double operator()(VectorN<N> x) const    { return 0.0; }
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

///////////////////////////   ./include/interfaces/ICoordSystemTransf.h   ///////////////////////////


namespace MML
{
    template<int N>
    class ICoordSystemTransf
    {
        public:
        virtual VectorN<N> transf(VectorN<N> in) = 0;
        virtual IScalarFunction<N>& coordTransfFunc(int i) = 0;
        virtual IScalarFunction<N>& inverseCoordTransfFunc(int i) = 0;
    };
}

///////////////////////////   ./include/basic_types/Tensors.h   ///////////////////////////

namespace MML
{
    template <int N>
    class TensorRank2
    {
        private:
            double _coeff[N][N];
            int _numContravar;
            int _numCovar;

            double ScalarProduct(VectorN<N> a, VectorN<N> b)
            {
                return 0.0;
            }

            TensorRank2 Transform(ICoordSystemTransf<N> &a)
            {
                return TensorRank2{};
            }
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

///////////////////////////   ./include/algorithms/Derivation.h   ///////////////////////////


namespace MML
{
    class Derivation
    {
    public:
        static double NDer1(const MML::IRealFunction &f, double x, double* error = nullptr)
        {
            using std::sqrt;
            using std::pow;
            using std::abs;
            using std::numeric_limits;

            const double eps = std::numeric_limits<double>::epsilon();
            // Error bound ~eps^1/2
            // Note that this estimate of h differs from the best estimate by a factor of sqrt((|f(x)| + |f(x+h)|)/|f''(x)|).
            // Since this factor is invariant under the scaling f -> kf, then we are somewhat justified in approximating it by 1.
            // This approximation will get better as we move to higher orders of accuracy.
            double h = 2 * sqrt(eps);
            //h = detail::make_xph_representable(x, h);

            double yh = f(x + h);
            double y0 = f(x);
            double diff = yh - y0;
            if (error)
            {
                double ym = f(x - h);
                double ypph = abs(yh - 2 * y0 + ym) / h;
                // h*|f''(x)|*0.5 + (|f(x+h)+|f(x)|)*eps/h
                *error = ypph / 2 + (abs(yh) + abs(y0))*eps / h;
            }
            return diff / h;
        }
            
        static double NDer1(const MML::IRealFunction &f, double x, double h, double* error = nullptr)
        {
            using std::abs;

            const double eps = std::numeric_limits<double>::epsilon();

            double yh = f(x + h);
            double y0 = f(x);
            double diff = yh - y0;
            if (error)
            {
                double ym = f(x - h);
                double ypph = abs(yh - 2 * y0 + ym) / h;
                // h*|f''(x)|*0.5 + (|f(x+h)+|f(x)|)*eps/h
                *error = ypph / 2 + (abs(yh) + abs(y0))*eps / h;
            }
            return diff / h;
        }

        template <int N>
        class HelperFunc
        {
            const IScalarFunction<N> &_func;
            const VectorN<N> &_x;
            int _der_ind;

        public:
            HelperFunc(const IScalarFunction<N> &inFunc, const VectorN<N> &x) : _func(inFunc), _x(x)
            { }
            HelperFunc(const IScalarFunction<N> &inFunc, const VectorN<N> &x, int der_ind) : _func(inFunc), _x(x), _der_ind(der_ind)
            { }

            void setDerInd(int ind) { _der_ind = ind; }

            double operator()(double x) const
            {
                VectorN<N> copy{_x};

                copy[_der_ind] = x;

                return _func(copy);
            }
        };

        template <int N>
        static double NDer1Partial(const IScalarFunction<N> &f, int deriv_index, const VectorN<N> &point, double h, double *error = nullptr)
        {
            double x = point[deriv_index];
            HelperFunc<N> helper{f, point, deriv_index};

            const double eps = std::numeric_limits<double>::epsilon();

            double yh = helper(x + h);
            double y0 = helper(x);
            double diff = yh - y0;
            if (error)
            {
                double ym = helper(x - h);
                double ypph = abs(yh - 2 * y0 + ym) / h;
                // h*|f''(x)|*0.5 + (|f(x+h)+|f(x)|)*eps/h
                *error = ypph / 2 + (abs(yh) + abs(y0))*eps / h;
            }
            return diff / h;
        }

        // BIG TODO - derivacija vektorskog polja PO ODREDJENOM SMJERU
        // template <int N>
        // static double DeriveDirectional(const VectorNFunction<N> &f, VectorN<N> &deriv_direction, VectorN<N> &point)
        // {
        //     double x = point[deriv_index];
        //     HelperFunc<N> helper{f, point, deriv_index};

        //     double dfdx = boost::math::differentiation::finite_difference_derivative(helper, x);

        //     return dfdx;
        // }

        template <int N>
        static VectorN<N> NDer1PartialByAll(const IScalarFunction<N> &f, const VectorN<N> &point, double h, VectorN<N> *error = nullptr)
        {
            VectorN<N> ret;

            for( int i=0; i<N; i++)
            {
                int deriv_index = i;

                double x = point[deriv_index];
                HelperFunc<N> helper{f, point, deriv_index};

                const double eps = std::numeric_limits<double>::epsilon();

                double yh = helper(x + h);
                double y0 = helper(x);
                double diff = yh - y0;
                if (error)
                {
                    double ym = helper(x - h);
                    double ypph = abs(yh - 2 * y0 + ym) / h;
                    // h*|f''(x)|*0.5 + (|f(x+h)+|f(x)|)*eps/h
                    (*error)[i] = ypph / 2 + (abs(yh) + abs(y0))*eps / h;
                }

                ret[i] = diff / h;
            }

            return ret;
        }
    };

}

///////////////////////////   ./include/algorithms/DiffEqSolvers.h   ///////////////////////////





namespace MML
{
    class RKSolver
    {
    public:

        void rk4(Vector &y, Vector &dydx, const double x, const double h,
            Vector &yout, void derivs(const double, Vector &, Vector &))
        {
            int i;
            double xh,hh,h6;

            int n= (int) y.size();
            Vector dym(n),dyt(n),yt(n);
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


        void rkdumb(Vector &vstart, const double x1, const double x2, Vector &xx, Matrix &y,
            void derivs(const double, Vector &, Vector &))
        {
            int i,k;
            double x,h;

            int nvar=y.RowNum();
            int nstep=y.ColNum()-1;
            Vector v(nvar),vout(nvar),dv(nvar);
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

        void rkck(Vector &y, Vector &dydx, const double x,
            const double h, Vector &yout, Vector &yerr,
            void derivs(const double, Vector &, Vector &))
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
            Vector ak2(n),ak3(n),ak4(n),ak5(n),ak6(n),ytemp(n);
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


        void rkqs(Vector &y, Vector &dydx, double &x, const double htry,
            const double eps, Vector &yscal, double &hdid, double &hnext,
            void derivs(const double, Vector &, Vector &))
        {
            const double SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4;
            int i;
            double errmax,h,htemp,xnew;

            int n= (int) y.size();
            h=htry;
            Vector yerr(n),ytemp(n);
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
        Vector *xp_p;
        Matrix *yp_p;

        void odeint(Vector &ystart, const double x1, const double x2, const double eps,
            const double h1, const double hmin, int &nok, int &nbad,
            void derivs(const double, Vector &, Vector &),
            void rkqs(Vector &, Vector &, double &, const double, const double,
            Vector &, double &, double &, void (*)(const double, Vector &, Vector &)))
        {
            const int MAXSTP=10000;
            const double TINY=1.0e-30;
            int i,nstp;
            double xsav,x,hnext,hdid,h;

            int nvar= (int) ystart.size();
            Vector yscal(nvar),y(nvar),dydx(nvar);
            Vector &xp=*xp_p;
            Matrix &yp=*yp_p;
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

///////////////////////////   ./include/algorithms/EigenSystemSolvers.h   ///////////////////////////
namespace MML
{
    class EigenSystemSolvers
    {

    };
}


///////////////////////////   ./include/algorithms/Interpolators.h   ///////////////////////////


namespace MML
{
    static bool polint(Vector &xa, Vector &ya, const double x, double &y, double &dy)
    // Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns a value y, and
    // an error estimate dy. If P (x) is the polynomial of degree N ??? 1 such that P (xai) = yai; i =
    // 1; : : : ; n, then the returned value y = P (x).
    {
        int i,m,ns=0;
        double den,dif,dift,ho,hp,w;

        int n=(int) xa.size();
        Vector c(n),d(n);
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

    static void polin2(Vector &x1a, Vector &x2a, Matrix &ya, const double x1,
                const double x2, double &y, double &dy)
    // Given arrays x1a[1..m] and x2a[1..n] of independent variables, and a submatrix of function
    // values ya[1..m][1..n], tabulated at the grid points defined by x1a and x2a; and given values
    // x1 and x2 of the independent variables; this routine returns an interpolated function value y,
    // and an accuracy indication dy (based only on the interpolation in the x1 direction, however).          
    {
        int j,k;

        int m = (int) x1a.size();
        int n = (int) x2a.size();
        Vector ymtmp(m),ya_t(n);
        for (j=0;j<m;j++) {
            for (k=0;k<n;k++) 
                ya_t[k]=ya[j][k];

            polint(x2a,ya_t,x2,ymtmp[j],dy);
        }
        polint(x1a,ymtmp,x1,y,dy);
    }

    static bool ratint(Vector &xa, Vector &ya, const double x, double &y, double &dy)
    {
        // Given arrays xa[1..n] and ya[1..n], and given a value of x, this routine returns a value of
        // y and an accuracy estimate dy. The value returned is that of the diagonal rational function,
        // evaluated at x, which passes through the n points (xai; yai), i = 1:::n.

        const double TINY=1.0e-25;
        int m,i,ns=0;
        double w,t,hh,h,dd;

        int n=(int)xa.size();
        Vector c(n),d(n);
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
	static double trapzd(double func(const double), const double a, const double b, const int n)
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
		} else {
			for (it=1,j=1;j<n-1;j++) it <<= 1;
			tnm=it;
			del=(b-a)/tnm;
			x=a+0.5*del;
			for (sum=0.0,j=0;j<it;j++,x+=del) sum += func(x);
			s=0.5*(s+(b-a)*sum/tnm);
			return s;
		}
	}

	static double qtrap(double func(const double), const double a, const double b)
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
		// errors may start increasing, and the routine may never converge. A value 10???6
		// is just on the edge of trouble for most 32-bit machines; it is achievable when the
		// convergence is moderately rapid, but not otherwise.

		const int JMAX=20;
		const double EPS=1.0e-10;
		int j;
		double s,olds=0.0;

		for (j=0;j<JMAX;j++) {
			s=trapzd(func,a,b,j+1);
			if (j > 5)
				if (fabs(s-olds) < EPS*fabs(olds) ||
					(s == 0.0 && olds == 0.0)) 
					return s;
			olds=s;
		}
		// nrerror("Too many steps in routine qtrap");
		return 0.0;
	}

	static double qsimp(double func(const double), const double a, const double b)
	{
		// Returns the integral of the function func from a to b. The parameters EPS can be set to the
		// desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum allowed
		// number of steps. Integration is performed by Simpson???s rule.

		// The routine qsimp will in general be more efficient than qtrap (i.e., require
		// fewer function evaluations) when the function to be integrated has a finite 4th
		// derivative (i.e., a continuous 3rd derivative). The combination of qsimp and its
		// necessary workhorse trapzd is a good one for light-duty work.

		const int JMAX=20;
		const double EPS=1.0e-10;
		int j;
		double s,st,ost=0.0,os=0.0;

		for (j=0;j<JMAX;j++) {
			st=trapzd(func,a,b,j+1);
			s=(4.0*st-ost)/3.0;
			if (j > 5)
				if (fabs(s-os) < EPS*fabs(os) ||
					(s == 0.0 && os == 0.0)) 
					return s;
			os=s;
			ost=st;
		}
		//nrerror("Too many steps in routine qsimp");
		return 0.0;
	}

	static double qromb(double func(const double), double a, double b)
	{
		// Returns the integral of the function func from a to b. Integration is performed by Romberg???s
		// method of order 2K, where, e.g., K=2 is Simpson???s rule.

		// The routine qromb, along with its required trapzd and polint, is quite
		// powerful for sufficiently smooth (e.g., analytic) integrands, integrated over intervals
		// which contain no singularities, and where the endoubleoints are also nonsingular. qromb,
		// in such circumstances, takes many, many fewer function evaluations than either of
		// the routines in x4.2

		const int JMAX=20, JMAXP=JMAX+1, K=5;
		const double EPS=1.0e-10;
		double ss,dss;
		Vector s(JMAX),h(JMAXP),s_t(K),h_t(K);
		int i,j;

		h[0]=1.0;
		for (j=1;j<=JMAX;j++) {
			s[j-1]=trapzd(func,a,b,j);
			if (j >= K) {
				for (i=0;i<K;i++) {
					h_t[i]=h[j-K+i];
					s_t[i]=s[j-K+i];
				}
				polint(h_t,s_t,0.0,ss,dss);
				if (fabs(dss) <= EPS*fabs(ss)) return ss;
			}
			h[j]=0.25*h[j-1];
		}
		//nrerror("Too many steps in routine qromb");
		return 0.0;
	}

} // end namespace
///////////////////////////   ./include/interfaces/ILinAlgEqSystemSolver.h   ///////////////////////////
namespace MML
{
    class ILinAlgEqSystemSolver
    {
        // zadaje mu se matrica sistema
        // dvije Solve, za Vector i Matrix
    };
}


///////////////////////////   ./include/algorithms/LinAlgEqSolvers.h   ///////////////////////////




namespace MML
{
    class GaussJordanSolver
    {
        public:

        static bool Solve(Matrix &a, Matrix &b)
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
                                if (fabs(a[j][k]) >= big) {
                                    big=fabs(a[j][k]);
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
                    //nrerror("gaussj: Singular Matrix");
                    return false;

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
    };

    class LUDecompositionSolver
    {
    private:
        int n;
        Matrix &refOrig;

        Matrix lu;
        std::vector<int> indx;
        double d;
    
    public:
        LUDecompositionSolver(Matrix  &a) : n(a.RowNum()), lu(a), refOrig(a), indx(n) 
        {
            // Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise
            // permutation of itself. a and n are input. a is output, arranged as in equation (2.3.14) above;
            // indx[1..n] is an output vector that records the row permutation effected by the partial
            // pivoting; d is output as ??1 depending on whether the number of row interchanges was even
            // or odd, respectively. This routine is used in combination with lubksb to solve linear equations
            // or invert a matrix.
            const double TINY=1.0e-40;
            int i,imax,j,k;
            double big,temp;
            Vector vv(n);
            d=1.0;
            for (i=0;i<n;i++) {
                big=0.0;
                for (j=0;j<n;j++)
                    if ((temp=abs(lu[i][j])) > big) big=temp;
                if (big == 0.0) 
                    throw("Singular matrix in LUdcmp");
                vv[i]=1.0/big;
            }
            for (k=0;k<n;k++) {
                big=0.0;
                imax=k;
                for (i=k;i<n;i++) {
                    temp=vv[i]*abs(lu[i][k]);
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
        void Solve(Vector &b, Vector &x)
        {
            // Solves the set of n linear equations A??X = B. Here a[1..n][1..n] is input, not as the matrix
            // A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input
            // as the permutation vector returned by ludcmp. b[1..n] is input as the right-hand side vector
            // B, and returns with the solution vector X. a, n, and indx are not modified by this routine
            // and can be left in place for successive calls with different right-hand sides b. This routine takes
            // into account the possibility that b will begin with many zero elements, so it is efficient for use
            // in matrix inversion
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

        void Solve(Matrix &b, Matrix &x)
        {
            int i,j,m=b.ColNum();
            
            if (b.RowNum() != n || x.RowNum() != n || b.ColNum() != x.ColNum())
                throw("LUdcmp::solve bad sizes");
            
            Vector xx(n);
            
            for (j=0;j<m;j++) {
                for (i=0;i<n;i++) 
                    xx[i] = b[i][j];
                
                Solve(xx,xx);
                
                for (i=0;i<n;i++) 
                    x[i][j] = xx[i];
            }
        }        

        void inverse(Matrix &ainv)
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
        void mprove(Vector &b, Vector &x)
        {
            // Improves a solution vector x[1..n] of the linear set of equations A ?? X = B. The matrix
            // a[1..n][1..n], and the vectors b[1..n] and x[1..n] are input, as is the dimension n.
            // Also input is alud[1..n][1..n], the LU decomposition of a as returned by ludcmp, and
            // the vector indx[1..n] also returned by that routine. On output, only x[1..n] is modified,
            // to an improved set of values
            int i,j;
            Vector r(n);
            
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
        Matrix el;
    
    public:    
        CholeskyDecompositionSolver(Matrix &a) : n(a.RowNum()), el(a) 
        {
            // Given a positive-definite symmetric matrix a[1..n][1..n], this routine constructs its Cholesky
            // decomposition, A = L ?? LT . On input, only the upper triangle of a need be given; it is not
            // modified. The Cholesky factor L is returned in the lower triangle of a, except for its diagonal
            // elements which are returned in p[1..n]
            int i,j,k;
            double sum;
            
            if (el.ColNum() != n) 
                throw("need square matrix");
            
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
        void solve(Vector &b, Vector &x) 
        {
            // Solves the set of n linear equations A ?? x = b, where a is a positive-definite symmetric matrix.
            // a[1..n][1..n] and p[1..n] are input as the output of the routine choldc. Only the lower
            // triangle of a is accessed. b[1..n] is input as the right-hand side vector. The solution vector is
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
        void elmult(Vector &y, Vector &b) {
            int i,j;
            if (b.size() != n || y.size() != n) throw("bad lengths");
            for (i=0;i<n;i++) {
                b[i] = 0.;
                for (j=0;j<=i;j++) b[i] += el[i][j]*y[j];
            }
        }
        void elsolve(Vector &b, Vector &y) {
            int i,j;
            double sum;
            if (b.size() != n || y.size() != n) throw("bad lengths");
            for (i=0;i<n;i++) {
                for (sum=b[i],j=0; j<i; j++) sum -= el[i][j]*y[j];
                y[i] = sum/el[i][i];
            }
        }
        void inverse(Matrix &ainv) {
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
        Matrix qt, r;
        bool sing;    

    public:

        QRDecompositionSolver(Matrix &a) : n(a.RowNum()), qt(n,n), r(a), sing(false) 
        {
            // Constructs the QR decomposition of a[1..n][1..n]. The upper triangular matrix R is returned in the upper triangle of a, except for the diagonal elements of R which are returned in
            // d[1..n]. The orthogonal matrix Q is represented as a product of n ??? 1 Householder matrices
            // Q1 : : : Qn???1, where Qj = 1 ??? uj ??? uj=cj. The ith component of uj is zero for i = 1; : : : ; j ??? 1
            // while the nonzero components are returned in a[i][j] for i = j; : : : ; n. sing returns as
            // true (1) if singularity is encountered during the decomposition, but the decomposition is still
            // completed in this case; otherwise it returns false (0).
            int i,j,k;
            Vector c(n), d(n);
            double scale,sigma,sum,tau;
            for (k=0;k<n-1;k++) {
                scale=0.0;
                for (i=k;i<n;i++) scale=std::max(scale,abs(r[i][k]));
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
        void Solve(Vector &b, Vector &x) 
        {
            // Solves the set of n linear equations A ?? x = b. a[1..n][1..n], c[1..n], and d[1..n] are
            // input as the output of the routine qrdcmp and are not modified. b[1..n] is input as the
            // right-hand side vector, and is overwritten with the solution vector on output.            
            qtmult(b,x);
            rsolve(x,x);
        }

        void qtmult(Vector &b, Vector &x) {
            int i,j;
            double sum;
            for (i=0;i<n;i++) {
                sum = 0.;
                for (j=0;j<n;j++) sum += qt[i][j]*b[j];
                x[i] = sum;
            }
        }

        void rsolve(Vector &b, Vector &x) 
        {
            // Solves the set of n linear equations R ?? x = b, where R is an upper triangular matrix stored in
            // a and d. a[1..n][1..n] and d[1..n] are input as the output of the routine qrdcmp and
            // are not modified. b[1..n] is input as the right-hand side vector, and is overwritten with the
            // solution vector on output            
            int i,j;
            double sum;
            if (sing) throw("attempting solve in a singular QR");
            for (i=n-1;i>=0;i--) {
                sum=b[i];
                for (j=i+1;j<n;j++) sum -= r[i][j]*x[j];
                x[i]=sum/r[i][i];
            }
        }
        void update(Vector &u, Vector &v) 
        {
            // Given the QR decomposition of some n ?? n matrix, calculates the QR decomposition of the
            // matrix Q??(R+u???v). The quantities are dimensioned as r[1..n][1..n], qt[1..n][1..n],
            // u[1..n], and v[1..n]. Note that QT is input and returned in qt.            
            int i,k;
            Vector w(u);
            for (k=n-1;k>=0;k--)
                if (w[k] != 0.0) break;
            if (k < 0) k=0;
            for (i=k-1;i>=0;i--) {
                rotate(i,w[i],-w[i+1]);
                if (w[i] == 0.0)
                    w[i]=abs(w[i+1]);
                else if (abs(w[i]) > abs(w[i+1]))
                    w[i]=abs(w[i])*sqrt(1.0+SQR(w[i+1]/w[i]));
                else w[i]=abs(w[i+1])*sqrt(1.0+SQR(w[i]/w[i+1]));
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
            // i and i + 1 of each matrix. a and b are the parameters of the rotation: cos ?? = a=pa2 + b2,
            // sin ?? = b=pa2 + b2.            
            int j;
            double c,fact,s,w,y;
            if (a == 0.0) {
                c=0.0;
                s=(b >= 0.0 ? 1.0 : -1.0);
            } else if (abs(a) > abs(b)) {
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

    // void tridag(Vector &a, Vector &b, Vector &c, Vector &r, Vector &u)
    // {
    //     // Solves for a vector u[1..n] the tridiagonal linear set given by equation (2.4.1). a[1..n],
    //     // b[1..n], c[1..n], and r[1..n] are input vectors and are not modified.

    //     int j;
    //     double bet;

    //     int n=(int)a.size();
    //     Vector gam(n);

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
        Matrix u,v;
        Vector w;
        double eps, tsh;
    
    public:
        SVDecompositionSolver(Matrix &a) : m(a.RowNum()), n(a.ColNum()), u(a), v(n,n), w(n) 
        {
            // Given a matrix a[1..m][1..n], this routine computes its singular value decomposition, A =
            // U??W ??V T . The matrix U replaces a on output. 
            // The diagonal matrix of singular values W is output as a vector w[1..n]. The matrix V (not the transpose V T ) is output as v[1..n][1..n].            
            eps = std::numeric_limits<double>::epsilon();
            decompose();
            reorder();
            tsh = 0.5*sqrt(m+n+1.)*w[0]*eps;
        }


        double inv_condition() {
            return (w[0] <= 0. || w[n-1] <= 0.) ? 0. : w[n-1]/w[0];
        }

    
    void solve(Vector &b, Vector &x, double thresh = -1.) 
    {
        // Solves A??X = B for a vector X, where A is specified by the arrays u[1..m][1..n], w[1..n],
        // v[1..n][1..n] as returned by svdcmp. m and n are the dimensions of a, and will be equal for
        // square matrices. b[1..m] is the input right-hand side. x[1..n] is the output solution vector.
        // No input quantities are destroyed, so the routine may be called sequentially with different b???s.        
        int i,j,jj;
        double s;
        if (b.size() != m || x.size() != n) throw("solve bad sizes");
        Vector tmp(n);
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

    void solve(Matrix &b, Matrix &x, double thresh = -1.)
    {
        int i,j,p=b.ColNum();
        if (b.RowNum() != m || x.RowNum() != n || x.ColNum() != p)
            throw("solve bad sizes");
        Vector xx(n),bcol(m);
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

    Matrix range(double thresh = -1.){
        int i,j,nr=0;
        Matrix rnge(m,rank(thresh));
        for (j=0;j<n;j++) {
            if (w[j] > tsh) {
                for (i=0;i<m;i++) rnge[i][nr] = u[i][j];
                nr++;
            }
        }
        return rnge;
    }

    Matrix nullspace(double thresh = -1.){
        int j,jj,nn=0;
        Matrix nullsp(n,nullity(thresh));
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
        Vector rv1(n);
        g = scale = anorm = 0.0;
        for (i=0;i<n;i++) {
            l=i+2;
            rv1[i]=scale*g;
            g=s=scale=0.0;
            if (i < m) {
                for (k=i;k<m;k++) scale += abs(u[k][i]);
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
                for (k=l-1;k<n;k++) scale += abs(u[i][k]);
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
            anorm=std::max(anorm,(abs(w[i])+abs(rv1[i])));
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
                        if (l == 0 || abs(rv1[l]) <= eps*anorm) {
                            flag=false;
                            break;
                        }
                        if (abs(w[nm]) <= eps*anorm) break;
                    }
                    if (flag) {
                        c=0.0;
                        s=1.0;
                        for (i=l;i<k+1;i++) {
                            f=s*rv1[i];
                            rv1[i]=c*rv1[i];
                            if (abs(f) <= eps*anorm) break;
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
            Vector su(m), sv(n);
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
            double absa=abs(a), absb=abs(b);
            return (absa > absb ? absa*sqrt(1.0+SQR(absb/absa)) :
                (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb))));
        }
    };
} // end namespace
///////////////////////////   ./include/algorithms/RootFinding.h   ///////////////////////////// bit ??e toga

namespace MML
{
    class RootFinding
    {

    };
}


///////////////////////////   ./include/systems/CoordSystem.h   ///////////////////////////



namespace MML
{
    class CoordTransfSphericalToCartesian : public ICoordSystemTransf<3>
    {
        // q[0] = r     - radial distance
        // q[1] = theta - inclination
        // q[2] = phi   - azimuthal angle
        public:
        static double func1(VectorN<3> &q) { return q[0] * sin(q[1]) * cos(q[2]); }
        static double func2(VectorN<3> &q) { return q[0] * sin(q[1]) * sin(q[2]); }
        static double func3(VectorN<3> &q) { return q[0] * cos(q[1]); }

        // q[0] = x
        // q[1] = y
        // q[2] = z
        static double funcInverse1(VectorN<3> &q) { return sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]); }
        static double funcInverse2(VectorN<3> &q) { return atan2(sqrt(q[0]*q[0] + q[1]*q[1]), q[2]); }
        static double funcInverse3(VectorN<3> &q) { return atan2(q[1], q[0]); }

        inline static ScalarFunctionFromStdFunc<3> _func[3] = { 
                                                    ScalarFunctionFromStdFunc<3>{std::function<double(VectorN<3>&)>{func1}},
                                                    ScalarFunctionFromStdFunc<3>{std::function<double(VectorN<3>&)>{func2}},
                                                    ScalarFunctionFromStdFunc<3>{std::function<double(VectorN<3>&)>{func3}}
                                                };

        inline static ScalarFunctionFromStdFunc<3> _funcInverse[3] = { 
                                                    ScalarFunctionFromStdFunc<3>{std::function<double(VectorN<3>&)>{funcInverse1}},
                                                    ScalarFunctionFromStdFunc<3>{std::function<double(VectorN<3>&)>{funcInverse2}},
                                                    ScalarFunctionFromStdFunc<3>{std::function<double(VectorN<3>&)>{funcInverse3}}
                                                };

        VectorN<3>           transf(VectorN<3> q)           { return VectorN<3>{ func1(q), func2(q), func3(q) }; }
        IScalarFunction<3>&  coordTransfFunc(int i)         { return _func[i]; }
        IScalarFunction<3>&  inverseCoordTransfFunc(int i)  { return _funcInverse[i]; }
    };

    class CoordTransfCartesianToSpherical : public ICoordSystemTransf<3>
    {
        // q[0] = x
        // q[1] = y
        // q[2] = z
        public:
        static double func1(VectorN<3> &q) { return sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]); }
        static double func2(VectorN<3> &q) { return atan2(sqrt(q[0]*q[0] + q[1]*q[1]), q[2]); }
        static double func3(VectorN<3> &q) { return atan2(q[1], q[0]); }

        // q[0] = r     - radial distance
        // q[1] = theta - inclination
        // q[2] = phi   - azimuthal angle
        public:
        static double funcInverse1(VectorN<3> &q) { return q[0] * sin(q[1]) * cos(q[2]); }
        static double funcInverse2(VectorN<3> &q) { return q[0] * sin(q[1]) * sin(q[2]); }
        static double funcInverse3(VectorN<3> &q) { return q[0] * cos(q[1]); }

        inline static ScalarFunctionFromStdFunc<3> _func[3] = { 
                                                    ScalarFunctionFromStdFunc<3>{std::function<double(VectorN<3>&)>{func1}},
                                                    ScalarFunctionFromStdFunc<3>{std::function<double(VectorN<3>&)>{func2}},
                                                    ScalarFunctionFromStdFunc<3>{std::function<double(VectorN<3>&)>{func3}}
                                                };

        inline static ScalarFunctionFromStdFunc<3> _funcInverse[3] = { 
                                                    ScalarFunctionFromStdFunc<3>{std::function<double(VectorN<3>&)>{funcInverse1}},
                                                    ScalarFunctionFromStdFunc<3>{std::function<double(VectorN<3>&)>{funcInverse2}},
                                                    ScalarFunctionFromStdFunc<3>{std::function<double(VectorN<3>&)>{funcInverse3}}
                                                };

        VectorN<3>           transf(VectorN<3> q)           { return VectorN<3>{ func1(q), func2(q), func3(q) }; }
        IScalarFunction<3>&  coordTransfFunc(int i)         { return _func[i]; }
        IScalarFunction<3>&  inverseCoordTransfFunc(int i)  { return _funcInverse[i]; }
    };
}

///////////////////////////   ./include/systems/LinAlgEqSystem.h   ///////////////////////////
namespace MML
{
    class LinAlgEqSystem
    {
        // solve by GJ
        // perform LU, QR, Cholesky, SVD decomposition
        // find eigenvalues

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
