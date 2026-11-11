#include <mml/base/Matrix/Matrix.h>
#include <mml/base/Matrix/MatrixNM.h>
#include <mml/base/Matrix/MatrixSym.h>
#include <mml/base/Matrix/MatrixTriDiag.h>
#include <mml/base/DifferentialGeometry_base.h>

#if defined(MML_LINEAR_ALG_DIRECT_H) || defined(MML_MATRIX_INVERSION_H) || defined(MML_CONSOLE_PRINTER_H)
#error Base containers must not depend on inversion algorithms or tools.
#endif

template<class Type>
concept ContainerInverse = requires(const Type& matrix) { matrix.Inverse(); };
template<class Type>
concept ContainerInvert = requires(Type& matrix) { matrix.Invert(); };
template<class Type>
concept LegacyContainerInverse = requires(const Type& matrix) { matrix.GetInverse(); };

static_assert(!ContainerInverse<MML::Matrix<MML::Real>>);
static_assert(!ContainerInverse<MML::MatrixNM<MML::Real, 2, 2>>);
static_assert(!ContainerInverse<MML::MatrixSym<MML::Real>>);
static_assert(!ContainerInverse<MML::TridiagonalMatrix<MML::Real>>);
static_assert(!ContainerInvert<MML::Matrix<MML::Real>>);
static_assert(!ContainerInvert<MML::MatrixNM<MML::Real, 2, 2>>);
static_assert(!LegacyContainerInverse<MML::MatrixNM<MML::Real, 2, 2>>);
static_assert(!LegacyContainerInverse<MML::MatrixSym<MML::Real>>);
static_assert(!LegacyContainerInverse<MML::TridiagonalMatrix<MML::Real>>);
