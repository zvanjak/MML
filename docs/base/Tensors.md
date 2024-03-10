# Tensor classes

Four classes representing tensors of different order, with different number of dimensions.
Capable of handling various combination of covariant and contravariant indices.

## Tensor types

**Tensor2\<int Dim>>**

~~~ c++
enum TensorIndexType { CONTRAVARIANT, COVARIANT };

template <int N>
class Tensor2 : public ITensor2<N>
{
    Real _coeff[N][N];
public:
    int _numContravar;
    int _numCovar;
    bool _isContravar[2];
    Tensor2(int nContra, int nCo) : _numContravar(nContra), _numCovar(nCo) 
    Tensor2(TensorIndexType first, TensorIndexType second) 
    
    int   NumContravar() const { return _numContravar; }
    int   NumCovar()     const { return _numCovar;}

    Real  Component(int i, int j) const { return _coeff[i][j]; }
    Real& Component(int i, int j)       { return _coeff[i][j]; }

    void   Print(std::ostream& stream, int width, int precision) const;
    friend std::ostream& operator<<(std::ostream& stream, const Tensor2 &a);
};
~~~

Similar definitions for:

**Tensor3\<int Dim>**

**Tensor4\<int Dim>**

**Tensor5\<int Dim>**

