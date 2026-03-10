# MML Release 1.2 Plan

## Overview

**Release 1.2 - Making MML complete and cross-platform numerical computing toolkit**

Release 1.2 is a major release focusing on:
- Technical debt cleanup
- Cross-platform improvements  
- Feature enhancements across all modules
- Documentation completeness
- Visualization improvements

**Status:** In Progress  
**Priority:** P0 (Critical)  
**Epic ID:** `MinimalMathLibrary-1isb`

---

## Release Summary

| Metric | Count |
|--------|-------|
| **Total Sub-Epics** | 8 |
| **Total Tasks** | ~170 |
| **Completed Tasks** | ~55 |
| **Remaining Tasks** | ~115 |
| **Estimated Effort** | 300+ hours |

---

## Sub-Epics Overview

### 1. Claude Sonnet 4.5 Improvements (`MinimalMathLibrary-4e0b`)

**Status:** In Progress | **Priority:** P0

Epic for implementing high-priority improvements identified in comprehensive Claude Sonnet 4.5 analysis.

**Goal:** Elevate MML from A (89/100) to A+ (92+/100)

| Status | Count |
|--------|-------|
| ✅ Completed | 16 |
| 🔄 Open | 7 |

#### Completed Tasks
- ✅ Standardize size accessors (size() everywhere)
- ✅ Standardize error handling (result types everywhere)
- ✅ Split BaseUtils into focused utility modules
- ✅ Implement concrete Tensor2 and Tensor3 classes
- ✅ Add config objects to all algorithms for consistency
- ✅ Standardize parameter passing (const ref for input, value for sink)
- ✅ Add Big-O complexity to all algorithm documentation
- ✅ Standardize result structs with full diagnostics
- ✅ Add LoadResult type for DataLoader consistency
- ✅ Fix unnecessary copies (add nodiscard, use move)
- ✅ Add sparse matrix support (CSR/CSC formats)
- ✅ Add DAE support (IODESystemDAE interface)
- ✅ Add event detection (IODESystemWithEvents)
- ✅ Split ComputationalGeometry.h into smaller headers

#### Remaining Tasks
- 🔄 Reduce header blast radius (forward declarations) - 8h
- 🔄 Add thread-safety documentation to all public APIs - 6h
- 🔄 Fix CSV escape sequences in Serializer - 6h
- 🔄 Add vcpkg and conan package manager support - 10h
- 🔄 Expand coordinate systems (elliptical, parabolic, toroidal) - 16h
- 🔄 Add streaming I/O for large files - 14h
- 🔄 Expand statistical distributions library - 12h

---

### 2. GPT-5.2 Analysis Improvements (`MinimalMathLibrary-om8d`)

**Status:** In Progress | **Priority:** P0

Implementation work derived from comprehensive GPT-5.2 deep-dive analysis.

| Status | Count |
|--------|-------|
| ✅ Completed | 30 |
| 🔄 Open | 13 |

#### Completed P0 Tasks (Correctness/Safety)
- ✅ Base header self-sufficiency (include what you use)
- ✅ Fix Complex + mixed-type numeric routines
- ✅ Remove header-scope using-namespace directives
- ✅ Fix incorrect base math formulas
- ✅ Fix GaussJordanSolver type narrowing
- ✅ Replace string-literal throws with MML exceptions
- ✅ Define singularity-handling policy
- ✅ Fix ODE step calculator output contract
- ✅ Harden BVPShootingMethod failure handling
- ✅ Clarify computational geometry boolean ops semantics
- ✅ Fix Timer Start/Reset semantics
- ✅ Unify Serializer return style

#### Completed P1 Tasks (Consistency/Maintainability)
- ✅ Establish library-wide equality/tolerance policy
- ✅ Normalize exception usage for base primitives
- ✅ Document Real configurability and guarantees
- ✅ Make ReferenceFrame3D parent/child relationships safe
- ✅ Audit CoordTransf return types
- ✅ Fix layering drift in MatrixUtils.h
- ✅ Deprecate implicit scalar conversions on Result structs
- ✅ Standardize algorithm result structs
- ✅ Extract shared numeric validation helpers
- ✅ Remove Hessenberg reduction duplication
- ✅ Clarify IInterval getLength semantics
- ✅ Document coordTransfFunc lifetime/ownership
- ✅ Make IParametricCurveParametrized evaluation const
- ✅ Contain Visualizer output paths
- ✅ Make ConsolePrinter borders robust
- ✅ Separate compute from printing in algorithms

#### Remaining Tasks
- 🔄 Reduce copying in hot base APIs (P2)
- 🔄 Add non-owning span/view APIs for vectors (P2)
- 🔄 Add constexpr/noexcept where safe (P2)
- 🔄 Add high-ROI core contract tests (P2)
- 🔄 Cache/batch expensive metric computations (P2)
- 🔄 Make MinimumEnclosingCircle deterministic (P2)
- 🔄 Create shared IParametrized base (P2)
- 🔄 Guard interface sampling helpers against numPoints < 2 (P2)
- 🔄 Add optional component-wise evaluation for IVectorFunctionNM (P1)
- 🔄 Reduce header blast radius of Serializer/DataLoader (P1)
- 🔄 Reduce include weight/coupling in interfaces (P2)
- 🔄 Standardize error handling boundary across Tools (P2)
- 🔄 Add robust escaping for ConsolePrinter exports (P2)

---

### 3. Sigma Engine Implementation (`MinimalMathLibrary-t886`)

**Status:** Open | **Priority:** P0

Transform ExpressionEvaluator into a full interactive mathematical environment.

| Status | Count |
|--------|-------|
| ✅ Completed | 3 |
| 🔄 Open | 4 |

#### Completed Phases
- ✅ Phase 1: Session & Variables Foundation
- ✅ Phase 2: User-Defined Functions

#### Remaining Phases
- 🔄 Phase 3: AST Layer for Symbolic Operations (2-3 weeks)
- 🔄 Phase 4: Symbolic Differentiation and Simplification (2-3 weeks)
- 🔄 Phase 5: Visualization Integration (2 weeks)
- 🔄 Phase 6: Polish, Documentation and Distribution (1-2 weeks)

**Key Features:**
- Variable persistence across expressions
- User-defined functions with parameters
- Symbolic differentiation and simplification
- Integrated visualization/plotting
- Session save/load (.sigma format)

---

### 4. float128 Support on GCC (`MinimalMathLibrary-uwwl`)

**Status:** Open | **Priority:** P0

Add comprehensive support for GCC's `__float128` type as the Real precision type.

| Status | Count |
|--------|-------|
| 🔄 In Progress | 1 |
| 🔄 Open | 6 |

#### Tasks
- 🔄 Add quadmath header and float128 typedef infrastructure (in progress)
- 🔄 Handle Complex type limitations with float128
- 🔄 Create MMLQuadMath.h adapter for quadmath functions
- 🔄 Fix numeric_limits and constants for float128
- 🔄 Update CMakeLists.txt for quadmath library linking
- 🔄 Implement I/O operations for float128 types
- 🔄 Add float128 testing and validation suite

**Key Challenges:**
- `std::complex<__float128>` doesn't exist
- `std::numeric_limits<__float128>` doesn't exist  
- Standard math functions don't work with `__float128`
- Custom I/O operators needed

---

### 5. Adjust docs_demos with Documentation (`MinimalMathLibrary-1isb.1`)

**Status:** Open | **Priority:** P1

Align all docs_demos with corresponding header files.

| Sub-Epic | Tasks |
|----------|-------|
| mml/base - Core data structures | 29 |
| mml/core - Mathematical operations | 17 |
| mml/algorithms - Numerical algorithms | 24 |
| mml/systems - Dynamical systems | 2 |
| mml/tools - Utilities | 6 |
| mml/interfaces - Abstract interfaces | 10 |
| **Total** | **88** |

Each task involves:
- Verifying existing demo covers all header functionality
- Creating new demos where missing
- Ensuring consistent documentation style

---

### 6. Preparing MML Visualizers for 1.2 (`MinimalMathLibrary-1isb.2`)

**Status:** Open | **Priority:** P1

Update and improve all MML visualizers across platforms.

#### Tasks
- 🔄 Update all MML visualizers with latest platform versions
  - WPF, FLTK, Qt on Windows, Linux, Mac
  - Unified look-and-feel
- 🔄 Define nice visualizations for each visualizer
  - 2-5 sample visualizations per visualizer
  - Select best two for visualization data folder
- 🔄 Make FLTK default backend for 2D on Linux and Mac
- 🔄 Add demo functions to visualization_examples app
  - Cross-platform demo functions
  - Handle FLTK 3D limitations gracefully

**Note:** Some tasks involve work in separate MML_Visualizers project.

---

### 7. Documentation Style Alignment (`MinimalMathLibrary-er1b`)

**Status:** ✅ Closed

Aligned all header documentation to `///` single-line Doxygen style.

---

### 8. PDE Solvers Documentation (`MinimalMathLibrary-uf4h`)

**Status:** ✅ Closed

Comprehensive documentation for PDE Solver Module covering:
- Sparse matrices
- Grid infrastructure
- Iterative solvers
- Elliptic, parabolic, hyperbolic PDEs

---

## Priority Matrix

### P0 - Critical (Must Complete for 1.2)

| Epic | Remaining Tasks | Est. Hours |
|------|-----------------|------------|
| Claude Sonnet 4.5 Improvements | 7 | 72 |
| GPT-5.2 Improvements (P0/P1) | 2 | 14 |
| Sigma Engine (Phases 3-6) | 4 | 200+ |
| float128 Support | 6 | 60 |

### P1 - High Priority

| Epic | Remaining Tasks | Est. Hours |
|------|-----------------|------------|
| docs_demos Alignment | 88 | 176 |
| Visualizers Preparation | 4 | 40 |
| GPT-5.2 Improvements (P1) | 2 | 16 |

### P2 - Medium Priority (Nice to Have)

| Epic | Remaining Tasks | Est. Hours |
|------|-----------------|------------|
| GPT-5.2 Improvements (P2) | 9 | 50 |
| Claude Sonnet 4.5 (P2 tasks) | 3 | 42 |

---

## Recommended Execution Order

### Phase 1: Foundation (Complete Current Work)
1. Complete float128 infrastructure task (in progress)
2. Finish remaining P0 tasks from GPT-5.2
3. Complete Claude Sonnet 4.5 remaining P1 tasks

### Phase 2: Core Features
1. Sigma Engine Phase 3 (AST Layer)
2. Sigma Engine Phase 4 (Symbolic Differentiation)
3. float128 remaining tasks

### Phase 3: Visualization & Documentation
1. Visualizers preparation (all 4 tasks)
2. Sigma Engine Phase 5 (Visualization Integration)
3. Begin docs_demos alignment (prioritize most-used modules)

### Phase 4: Polish & Release
1. Complete docs_demos alignment
2. Sigma Engine Phase 6 (Polish)
3. P2 improvements as time permits
4. Package manager support (vcpkg/conan)

---

## Key Milestones

| Milestone | Target |
|-----------|--------|
| P0 critical tasks complete | Q1 2026 |
| Sigma Engine Phase 4 complete | Q2 2026 |
| Visualization improvements complete | Q2 2026 |
| docs_demos 50% complete | Q2 2026 |
| float128 support complete | Q2 2026 |
| Release 1.2 candidate | Q3 2026 |

---

## Risk Areas

1. **Sigma Engine Complexity** - AST and symbolic differentiation are complex features
2. **float128 Cross-Platform** - GCC-only feature, may complicate builds
3. **docs_demos Volume** - 88 tasks is substantial documentation work
4. **Visualizer Dependencies** - Work spans multiple projects (MML_Visualizers)

---

## Success Criteria

- [ ] All P0 tasks completed
- [ ] Sigma Engine fully functional with sessions, user functions, differentiation, visualization
- [ ] float128 support working on GCC with proper fallback
- [ ] All visualizers updated and working on all platforms
- [ ] At least 70% of docs_demos aligned
- [ ] All 3293+ existing tests continue to pass
- [ ] Cross-platform builds working (Windows/Linux/macOS)
- [ ] Package manager support (vcpkg or conan)

---

*Document generated: January 27, 2026*  
*Based on Beads issue tracker analysis*
