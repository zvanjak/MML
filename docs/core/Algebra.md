# Computational Algebra And Symmetry

MML's algebra layer makes symmetry computational without attempting to become a full
computer algebra system. Algebraic value types, finite algorithms, and symmetry
operations live in `mml/base`.

## Entry Points

```cpp
#include <mml/base/Algebra_base.h>
```

The implemented capability includes:

- permutations, finite groups, $C_n$, and $D_n$;
- orbits, stabilizers, Burnside counting, Cayley tables, and conjugacy classes;
- exact modular rings, prime fields, polynomials, and $GF(p^n)$;
- representations, characters, invariant projections, and group averaging;
- $SO(2)$, $SO(3)$, $SE(2)$, and $SE(3)$;
- lattices, reciprocal lattices, affine symmetries, and motif expansion.

All compositions use `compose(left, right)`, with `right` applied first for
transformations. Exact structures use exact equality; numerical algorithms expose
tolerance or equality policies.

## Worked Demo

`Docs_Demo_Algebra()` in `src/docs_demos/core/docs_demo_algebra.cpp` demonstrates:

1. cyclic and dihedral groups;
2. polygon orbits and stabilizers;
3. binary-necklace counting with Burnside's lemma;
4. an $SO(3)$ axis-angle rotation;
5. lattice motif expansion modulo unit-cell translations.

Build it through the normal examples gate:

```powershell
cmake --build build --config Release --target MML_DocsApp --parallel
```

Detailed API conventions and examples are in [../base/Algebra.md](../base/Algebra.md).
The original scope and implementation ledger are in
[../improvements_in_beads/Algebra_improvement.md](../improvements_in_beads/Algebra_improvement.md).

## Boundaries

- Named space-group databases, lattice reduction, and large computational-group
  algorithms remain out of scope.
- Numerical `Polynom` remains the interpolation/calculus type; exact
  `Algebra::Polynomial` serves finite algebra.
- Lie and rigid-motion types integrate with typed points and tangent vectors while
  remaining independent of Core coordinate-map machinery.
- Cayley graph and interpolation-frame data are renderer-neutral; visualization tools
  may serialize them without introducing visualization dependencies into Base/Core.