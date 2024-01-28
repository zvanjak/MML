# Polynom class

Class representing general polynom

~~~ c++
template <typename _Field, typename _CoefType = Real>
class Polynom
{
private:
    std::vector<_CoefType> _vecCoef;

public:
    int  GetDegree() const     { return (int) _vecCoef.size() - 1; }
    void SetDegree(int newDeg) {  _vecCoef.resize(newDeg+1); }

    _Field  operator[] (int i) const { return _vecCoef[i]; }
    _Field& operator[] (int i)       { return _vecCoef[i]; }

    _Field operator() (const _Field &x);

    // Given the coefficients of a polynomial of degree nc as an array c[0..nc] of size nc+1 (with
    // c[0] being the constant term), and given a value x, this routine fills an output array pd of size
    // nd+1 with the value of the polynomial evaluated at x in pd[0], and the first nd derivatives at
    // x in pd[1..nd].
    void Derive(const Real x, Vector<Real> &pd);

    bool operator==(const Polynom &b) const;

    Polynom operator+(const Polynom &b) const;
    Polynom operator-(const Polynom &b) const;
    Polynom operator*(const Polynom &b) const;

    static void poldiv(const Polynom &u, const Polynom &v, Polynom &qout, Polynom &rout);

    friend Polynom operator*(const Polynom &a, _CoefType b );
    friend Polynom operator*(_CoefType a, const Polynom &b );
    friend Polynom operator/(const Polynom &a, _CoefType b);

    std::string to_string(int width, int precision) const;
    std::ostream& Print(std::ostream& stream, int width, int precision) const;
    std::ostream& Print(std::ostream& stream) const;
    friend std::ostream& operator<<(std::ostream& stream, Polynom &a);
};

// predefined typedefs
typedef Polynom<Real, Real>         RealPolynom;
typedef Polynom<Complex, Complex>   ComplexPolynom;

typedef Polynom<MatrixNM<Real,2,2>, Real>       MatrixPolynomDim2;
typedef Polynom<MatrixNM<Real,3,3>, Real>       MatrixPolynomDim3;
typedef Polynom<MatrixNM<Real,4,4>, Real>       MatrixPolynomDim4;
~~~
