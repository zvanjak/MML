# Algebra Boundaries

MML core retains the SO2/SO3 rotation groups and SE2/SE3 rigid-motion value types,
including typed point/tangent-vector bridges. Include `mml/base/Algebra_base.h` or
the individual headers below `mml/base/Algebra/LieGroups/`.

Finite/discrete algebra is supplied by `MML_Ext` in the private MML-Packages-Private
repository, not by core or its single-header distribution. Its new aggregate is
`mml_ext/base/Algebra_base.h`; all 17 relocated component headers keep MML namespaces
and include guards. The two aggregate headers have distinct guards and can coexist.

Generic exact field-matrix determinant and inverse remain core algorithms at
`mml/algorithms/MatrixAlg/FieldMatrixAlgorithms.h`. Their namespace and arithmetic
are unchanged. Both core Basis and optional finite-field matrices reuse them.

Numerical `Polynom` stays core; exact `Algebra::Polynomial` belongs to the extension.
No compatibility headers or upward dependencies remain.

See [Core Algebra](../base/Algebra.md) for the API boundary and retained Lie-group
conventions. Run `MML_DocsApp --algebra` for the core example, or the receiving
repository's `MML_DocsApp --discrete-algebra` for the transferred group/action and
finite-field examples plus all six cookbook recipes.
