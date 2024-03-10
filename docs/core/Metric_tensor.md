# Metric tensor

Class inheriteed from ITensorField2, representing metric tensor of a given space. 

Base (abstract) class

~~~C++
template<int N>
class MetricTensorField : public ITensorField2<N>
{
public:
	MetricTensorField() : ITensorField2<N>(2, 0) { }
	MetricTensorField(int numContra, int numCo) : ITensorField2<N>(numContra, numCo) { }
  
	// this function, defined in ITensorField2, is not implemented, and should be implemented in derived class
	// virtual Real    Component(int i, int j, const VectorN<Real, N>& pos) const = 0;

	Tensor2<N>   operator()(const VectorN<Real, N>& pos) const;

	Real GetChristoffelSymbolFirstKind(int i, int j, int k, const VectorN<Real, N>& pos) const
	Real GetChristoffelSymbolSecondKind(int i, int j, int k, const VectorN<Real, N>& pos) const

	VectorN<Real, N> CovariantDerivativeContravar(const IVectorFunction<N>& func, int j, const VectorN<Real, N>& pos) const
	Real             CovariantDerivativeContravarComp(const IVectorFunction<N>& func, int i, int j, const VectorN<Real, N>& pos) const

	VectorN<Real, N> CovariantDerivativeCovar(const IVectorFunction<N>& func, int j, const VectorN<Real, N>& pos) const
	Real             CovariantDerivativeCovarComp(const IVectorFunction<N>& func, int i, int j, const VectorN<Real, N>& pos) const
};
~~~

## Defined metric tensors

~~~C++
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
~~~

## General metric tensor, defined with coord. transformation

~~~C++
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
~~~

## Example usage



