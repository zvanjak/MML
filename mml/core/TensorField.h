///////////////////////////////////////////////////////////////////////////////////////////
///                         MinimalMathLibrary (MML)                                  ///
///                                                                                   ///
///  File:        TensorField.h                                                       ///
///  Description: Tensor field coordinate transformation adapters                     ///
///               Express tensor fields in different coordinate systems               ///
///                                                                                   ///
///  Copyright:   (c) 2024-2026 Zvonimir Vanjak                                       ///
///  License:     MIT License (see LICENSE.md)                                         ///
///                                                                                   ///
///////////////////////////////////////////////////////////////////////////////////////////
#if !defined MML_TENSOR_FIELD_H
#define MML_TENSOR_FIELD_H

#ifdef MML_USE_SINGLE_HEADER
#include "MML.h"
#else
#include "interfaces/ITensorField.h"
#include "core/CoordTransf.h"
#endif

namespace MML
{
	/// @brief Adapter that expresses a rank-2 tensor field in a different coordinate system
	///
	/// Given a rank-2 tensor field T defined in "From" coordinates and an invertible
	/// coordinate transformation From→To, produces a new ITensorField2 that can be
	/// evaluated directly in "To" coordinates.
	///
	/// At each evaluation point p_to in the target coordinate system:
	///   1. Maps p_to back to source coordinates: p_from = transfInverse(p_to)
	///   2. Evaluates the original field:         T = originalField(p_from)
	///   3. Applies the tensor transformation law: T' = J · T · J^T  (Jacobian-based)
	///
	/// Example: Express an EM field tensor from Cartesian to Spherical coordinates
	/// @code
	///   SomeEMField  F_cart(/* ... */);           // ITensorField2<3> in Cartesian
	///   CoordTransfCartesianToSpherical transf;   // with inverse
	///
	///   TransformedTensorField2<Vector3Cartesian, Vector3Spherical, 3>
	///       F_sph(F_cart, transf);
	///
	///   Tensor2<3> T = F_sph(VectorN<Real, 3>{2.0, PI/4, PI/3});  // evaluate in spherical
	/// @endcode
	///
	/// @tparam VectorFrom Source coordinate system vector type
	/// @tparam VectorTo   Target coordinate system vector type
	/// @tparam N          Dimension of coordinate space
	template<typename VectorFrom, typename VectorTo, int N>
	class TransformedTensorField2 : public ITensorField2<N>
	{
		const ITensorField2<N>& _originalField;
		const CoordTransfWithInverse<VectorFrom, VectorTo, N>& _coordTransf;

	public:
		TransformedTensorField2(
			const ITensorField2<N>& originalField,
			const CoordTransfWithInverse<VectorFrom, VectorTo, N>& coordTransf)
			: ITensorField2<N>(originalField.getNumContravar(), originalField.getNumCovar())
			, _originalField(originalField)
			, _coordTransf(coordTransf)
		{}

		Tensor2<N> operator()(const VectorN<Real, N>& pos) const override
		{
			VectorTo   posTo(pos);
			VectorFrom posFrom = _coordTransf.transfInverse(posTo);
			Tensor2<N> tensor  = _originalField(posFrom);
			return _coordTransf.transfTensor2(tensor, posFrom);
		}

		Real Component(int i, int j, const VectorN<Real, N>& pos) const override
		{
			return (*this)(pos)(i, j);
		}
	};

	/// @brief Adapter that expresses a rank-3 tensor field in a different coordinate system
	///
	/// Same principle as TransformedTensorField2 but for rank-3 tensor fields
	/// (e.g., Christoffel symbols).
	///
	/// @tparam VectorFrom Source coordinate system vector type
	/// @tparam VectorTo   Target coordinate system vector type
	/// @tparam N          Dimension of coordinate space
	template<typename VectorFrom, typename VectorTo, int N>
	class TransformedTensorField3 : public ITensorField3<N>
	{
		const ITensorField3<N>& _originalField;
		const CoordTransfWithInverse<VectorFrom, VectorTo, N>& _coordTransf;

	public:
		TransformedTensorField3(
			const ITensorField3<N>& originalField,
			const CoordTransfWithInverse<VectorFrom, VectorTo, N>& coordTransf)
			: ITensorField3<N>(originalField.getNumContravar(), originalField.getNumCovar())
			, _originalField(originalField)
			, _coordTransf(coordTransf)
		{}

		Tensor3<N> operator()(const VectorN<Real, N>& pos) const override
		{
			VectorTo   posTo(pos);
			VectorFrom posFrom = _coordTransf.transfInverse(posTo);
			Tensor3<N> tensor  = _originalField(posFrom);
			return _coordTransf.transfTensor3(tensor, posFrom);
		}

		Real Component(int i, int j, int k, const VectorN<Real, N>& pos) const override
		{
			return (*this)(pos).Component(i, j, k);
		}
	};

	/// @brief Adapter that expresses a rank-4 tensor field in a different coordinate system
	///
	/// Same principle as TransformedTensorField2 but for rank-4 tensor fields
	/// (e.g., Riemann curvature tensor).
	///
	/// @tparam VectorFrom Source coordinate system vector type
	/// @tparam VectorTo   Target coordinate system vector type
	/// @tparam N          Dimension of coordinate space
	template<typename VectorFrom, typename VectorTo, int N>
	class TransformedTensorField4 : public ITensorField4<N>
	{
		const ITensorField4<N>& _originalField;
		const CoordTransfWithInverse<VectorFrom, VectorTo, N>& _coordTransf;

	public:
		TransformedTensorField4(
			const ITensorField4<N>& originalField,
			const CoordTransfWithInverse<VectorFrom, VectorTo, N>& coordTransf)
			: ITensorField4<N>(originalField.getNumContravar(), originalField.getNumCovar())
			, _originalField(originalField)
			, _coordTransf(coordTransf)
		{}

		Tensor4<N> operator()(const VectorN<Real, N>& pos) const override
		{
			VectorTo   posTo(pos);
			VectorFrom posFrom = _coordTransf.transfInverse(posTo);
			Tensor4<N> tensor  = _originalField(posFrom);
			return _coordTransf.transfTensor4(tensor, posFrom);
		}

		Real Component(int i, int j, int k, int l, const VectorN<Real, N>& pos) const override
		{
			return (*this)(pos)[i][j][k][l];
		}
	};

	/// @brief Adapter that expresses a rank-5 tensor field in a different coordinate system
	///
	/// Same principle as TransformedTensorField2 but for rank-5 tensor fields.
	///
	/// @tparam VectorFrom Source coordinate system vector type
	/// @tparam VectorTo   Target coordinate system vector type
	/// @tparam N          Dimension of coordinate space
	template<typename VectorFrom, typename VectorTo, int N>
	class TransformedTensorField5 : public ITensorField5<N>
	{
		const ITensorField5<N>& _originalField;
		const CoordTransfWithInverse<VectorFrom, VectorTo, N>& _coordTransf;

	public:
		TransformedTensorField5(
			const ITensorField5<N>& originalField,
			const CoordTransfWithInverse<VectorFrom, VectorTo, N>& coordTransf)
			: ITensorField5<N>(originalField.getNumContravar(), originalField.getNumCovar())
			, _originalField(originalField)
			, _coordTransf(coordTransf)
		{}

		Tensor5<N> operator()(const VectorN<Real, N>& pos) const override
		{
			VectorTo   posTo(pos);
			VectorFrom posFrom = _coordTransf.transfInverse(posTo);
			Tensor5<N> tensor  = _originalField(posFrom);
			return _coordTransf.transfTensor5(tensor, posFrom);
		}

		Real Component(int i, int j, int k, int l, int m, const VectorN<Real, N>& pos) const override
		{
			return (*this)(pos)[i][j][k][l][m];
		}
	};
}

#endif // MML_TENSOR_FIELD_H
