# Function Interfaces, Bases, and the FunctionSpaces Extension

MML core retains the reusable mathematical primitives:

- Function interfaces in `mml/interfaces/IFunction.h`.
- Function wrappers in `mml/base/Function.h`.
- Orthogonal bases in `mml/core/OrthogonalBasis/`.

The finite FunctionSpaces, operator-discretization and collocation BVP framework
has moved to **MML-Packages-Private**, its canonical implementation home:

- Aggregate: `include/mml_ext/core/FunctionSpaces.h`.
- Components: `include/mml_ext/core/FunctionSpaces/`.
- CMake target: `MML_Ext` (also `MML::Ext`).
- Namespaces and include guards are unchanged.

Switch the include prefix from `mml/core/FunctionSpaces` to
`mml_ext/core/FunctionSpaces` and link the extension target. The old core headers
and single-header definitions are removed, without forwarding headers or a core
dependency on packages.

Implementation, tests and runnable examples are maintained together in the
[private extension documentation](https://github.com/zvanjak/MML-Packages-Private/blob/master/docs/mml_ext/Function_spaces.md).
This migration does not publish private package source into the public companion
repositories.

For core primitives, see the [function interfaces](../../mml/interfaces/IFunction.h) and
[OrthogonalBases](OrthogonalBases.md). Historical size reports describe their
dated snapshots, not the post-relocation core.
