import os
from pathlib import Path

# Get the script directory and project root
script_dir = Path(__file__).parent
project_root = script_dir.parent.parent

print(f"Script dir: {script_dir}")
print(f"Project root: {project_root}")

# =============================================================================
# MML Single Header File List
#
# Files are ordered by dependencies - foundational types first, then interfaces,
# base classes, core algorithms, and finally high-level systems.
#
# AGGREGATE HEADERS are commented out with [AGGREGATE] marker - they only include
# other headers that are already listed individually. They exist for user convenience
# but are not needed in the single-header build.
#
# Last updated: 2026-07-03 (SparseMatrix and SparseSolvers added)
# =============================================================================

listFiles = [
    # =========================================================================
    # FOUNDATION - Core definitions and precision
    # =========================================================================
    # IMPORTANT: MMLTypeDefs.h must come first because it defines Real.
    # MMLPrecision.h must come before MMLBase.h because MMLBase.h uses
    # PrecisionValues<Real> for default tolerance values in AlgorithmContext
    # and Defaults struct once modular includes are stripped.
    "mml/MMLTypeDefs.h",
    "mml/MMLExceptions.h",
    "mml/MMLPrecision.h",
    "mml/MMLConcepts.h",
    "mml/MMLBase.h",
    "mml/MMLNumericValidation.h",
    "mml/MMLSingularityHandling.h",
    "mml/MMLVisualizators.h",

    # =========================================================================
    # INTERFACES - Abstract base classes (part 1 - no dependencies)
    # =========================================================================
    "mml/interfaces/IInterval.h",

    # =========================================================================
    # BASE - Fundamental data types and structures
    # =========================================================================
    # Intervals
    "mml/base/Intervals.h",

    # Algebra metadata and foundational contracts
    "mml/base/Algebra/AlgebraTraits.h",
    "mml/base/Algebra/Permutation.h",
    "mml/base/Algebra/FiniteGroup.h",
    "mml/base/Algebra/CyclicGroup.h",
    "mml/base/Algebra/DihedralGroup.h",
    "mml/base/Algebra/GroupAction.h",
    "mml/base/Algebra/ModInt.h",
    "mml/base/Algebra/ModularStructures.h",
    "mml/base/Algebra/Polynomial.h",
    "mml/base/Algebra/FiniteField.h",

    # Matrix print formatting (standalone, no MML deps)
    "mml/base/Matrix/MatrixPrintFormat.h",

    # Standard math functions
    "mml/base/StandardFunctions.h",
    "mml/base/SpecialFunctions.h",

    # Number theory, exact rationals, combinatorics (NumberTheory before Rational before Combinatorics)
    "mml/base/NumberTheory.h",
    "mml/base/Rational.h",
    "mml/base/Combinatorics.h",

    # Linear algebra - Vectors FIRST (needed by ITensor, IParametrized, Geometry, etc.)
    "mml/base/Vector/Vector.h",
    "mml/base/Vector/VectorN.h",

    # Typed geometry frame tags (needed by rank-1 tensor aliases)
    "mml/base/DifferentialGeometry/FrameTags.h",

    # IParametrized needs Vector (but not VectorN)
    "mml/interfaces/IParametrized.h",

    # Geometry fundamentals - Core geometry components (in GeometryBase/ subfolder)
    # [AGGREGATE] "mml/base/Geometry/Geometry.h",  # Only includes GeometryBase/* headers below
    "mml/base/Geometry/GeometryBase/GeometryPoints.h",
    "mml/base/Geometry/GeometryBase/Geometry2DShapes.h",
    "mml/base/Geometry/GeometryBase/Geometry3DShapes.h",

    # Linear algebra - Matrices (now in Matrix/ subfolder)
    "mml/base/Matrix/Matrix.h",
    "mml/base/Matrix/Matrix3D.h",
    "mml/base/Matrix/MatrixNM.h",

    # Algebra representations depend on fixed-size matrices and vectors
    "mml/base/Algebra/Representation.h",
    "mml/base/Matrix/MatrixSym.h",
    "mml/base/Matrix/MatrixTriDiag.h",
    "mml/base/Matrix/MatrixBandDiag.h",

    # Sparse matrices - Individual components (in SparseMatrix/ subfolder)
    "mml/base/SparseMatrix/SparseMatrixCOO.h",
    "mml/base/SparseMatrix/SparseMatrixCSR.h",
    "mml/base/SparseMatrix/SparseMatrixCSC.h",
    "mml/base/SparseMatrix/SparseMatrix.h",

    # Typed rank-1 tensor and typed geometry primitives
    "mml/base/Tensor/Rank1Tensor.h",
    "mml/base/DifferentialGeometry/Point.h",
    "mml/base/DifferentialGeometry/Metric.h",
    "mml/base/DifferentialGeometry/DifferentialForm.h",
    "mml/base/DifferentialGeometry/Hodge.h",

    # ITensor interface (needs VectorN)
    "mml/interfaces/ITensor.h",

    # Symbol utilities are used by tensor Levi-Civita helpers
    "mml/base/BaseUtils/SymbolUtils.h",

    # Tensors
    "mml/base/Tensor/MatrixNDim.h",
    "mml/base/Tensor/TensorCommon.h",
    "mml/base/Tensor/Tensor2.h",
    "mml/base/Tensor/Tensor3.h",
    "mml/base/Tensor/Tensor4.h",
    "mml/base/Tensor/Tensor5.h",

    # Shared algorithm config/result contracts used by base-level detailed APIs
    "mml/base/AlgorithmTypes.h",

    # Interfaces dependent on base objects
    "mml/interfaces/IFunction.h",
    "mml/base/DifferentialGeometry/TypedFields.h",
    "mml/interfaces/IComplexFunction.h",
    "mml/interfaces/IODESystem.h",
    "mml/interfaces/IODESystemWithEvents.h",
    "mml/interfaces/IDynamicalSystem.h",
    "mml/interfaces/IDiscreteMap.h",
    "mml/interfaces/IODESystemDAE.h",
    "mml/interfaces/IODESystemDAEWithEvents.h",
    "mml/interfaces/ICoordTransf.h",
    "mml/interfaces/ITensorField.h",
    "mml/interfaces/IODESystemStepCalculator.h",
    "mml/interfaces/IODESystemStepper.h",

    # Basic utilities (MUST come before VectorTypes which uses Utils::SphericalToCartesian)
    # [AGGREGATE] "mml/base/BaseUtils.h",  # Only includes BaseUtils/* headers below
    "mml/base/BaseUtils/AngleCoordUtils.h",
    "mml/base/BaseUtils/ComparisonUtils.h",
    "mml/base/BaseUtils/VectorOps.h",
    "mml/base/BaseUtils/MatrixOps.h",
    "mml/base/BaseUtils/MixedTypeOps.h",
    "mml/base/Random.h",
    "mml/base/QuasiRandom.h",

    # Dimensional vector types (3D uses coordinate conversion utilities)
    "mml/base/Vector/VectorTypes2D.h",
    "mml/base/Vector/VectorTypes3D.h",
    "mml/base/Vector/VectorTypes4D.h",

    "mml/base/Geometry/GeometrySpherical.h",

    # Quaternions
    "mml/base/Quaternions.h",

    # Geometry Lie groups depend on fixed-size vectors/matrices and quaternions
    "mml/base/Algebra/LieGroups/SO2.h",
    "mml/base/Algebra/LieGroups/SO3.h",
    "mml/base/Algebra/LieGroups/SE2.h",
    "mml/base/Algebra/LieGroups/SE3.h",

    # Polynomials
    "mml/base/Polynom.h",
    "mml/base/ChebyshevPolynom.h",

    # 2D Geometry - Individual components (in Geometry2D/ subfolder)
    # [AGGREGATE] "mml/base/Geometry/Geometry2D.h",  # Only includes Geometry2D/* headers below
    "mml/base/Geometry/Geometry2D/Geometry2DLines.h",
    "mml/base/Geometry/Geometry2D/Geometry2DTriangle.h",
    "mml/base/Geometry/Geometry2D/Geometry2DPolygon.h",
    "mml/base/Geometry/Geometry2D/Geometry2DCircleBox.h",

    # 3D Geometry - Individual components (in Geometry3D/ subfolder)
    # [AGGREGATE] "mml/base/Geometry/Geometry3D.h",  # Only includes Geometry3D/* headers below
    "mml/base/Geometry/Geometry3D/Geometry3DLines.h",
    "mml/base/Geometry/Geometry3D/Geometry3DPlane.h",
    "mml/base/Geometry/Geometry3D/Geometry3DSurfaces.h",

    # 3D Bodies - Individual components (in Geometry3DBodies/ subfolder)
    # [AGGREGATE] "mml/base/Geometry/Geometry3DBodies.h",  # Only includes Geometry3DBodies/* headers below
    "mml/base/Geometry/Geometry3DBodies/Geometry3DBounding.h",
    "mml/base/Geometry/Geometry3DBodies/Geometry3DBodyBase.h",
    "mml/base/Geometry/Geometry3DBodies/Geometry3DCubes.h",
    "mml/base/Geometry/Geometry3DBodies/Geometry3DTorusCylinder.h",
    "mml/base/Geometry/Geometry3DBodies/Geometry3DSpherePyramid.h",

    # Functions
    "mml/base/Function.h",
    "mml/base/ComplexFunction.h",

    # Interpolation - Individual components (in InterpolatedFunctions/ subfolder)
    # [AGGREGATE] "mml/base/InterpolatedFunction.h",  # Only includes InterpolatedFunctions/* headers below
    "mml/base/InterpolatedFunctions/InterpolationTypes.h",
    "mml/base/InterpolatedFunctions/InterpolatedRealFunctionLinear.h",
    "mml/base/InterpolatedFunctions/InterpolatedRealFunctionPolynomial.h",
    "mml/base/InterpolatedFunctions/InterpolatedRealFunctionRational.h",
    "mml/base/InterpolatedFunctions/InterpolatedRealFunctionBarycentric.h",
    "mml/base/InterpolatedFunctions/InterpolatedFunctionSpline.h",
    "mml/base/InterpolatedFunctions/InterpolatedRealFunctionAdvanced.h",
    "mml/base/InterpolatedFunctions/Interpolation2DFunction.h",
    "mml/base/InterpolatedFunctions/InterpolationParametricCurve.h",

    # Dirac delta and special functions
    "mml/base/DiracDeltaFunction.h",

    # ODE System base
    "mml/base/ODESystem.h",
    "mml/base/ODESystemSolution.h",
    "mml/base/DAESystem.h",

    # Graph
    "mml/base/Graph.h",

    # Chebyshev approximation
    "mml/base/ChebyshevApproximation.h",

    # =========================================================================
    # CORE - Core algorithms and operations
    # =========================================================================
    # Algebra law checks and operations over base algebraic structures
    "mml/base/Algebra/Algorithms/AlgebraLawChecks.h",
    "mml/base/Algebra/Algorithms/FiniteGroupAlgorithms.h",
    "mml/base/Algebra/Algorithms/GroupActionAlgorithms.h",
    "mml/base/Algebra/Algorithms/ModularArithmeticAlgorithms.h",
    "mml/base/Algebra/Algorithms/FieldMatrixAlgorithms.h",
    "mml/base/Algebra/Algorithms/PolynomialAlgorithms.h",
    "mml/base/Algebra/Algorithms/RepresentationAlgorithms.h",

    # Linear algebra solvers - Individual components (in LinAlgEqSolvers/ subfolder)
    # [AGGREGATE] "mml/core/LinAlgEqSolvers.h",  # Only includes LinAlgEqSolvers/* headers below
    "mml/core/LinAlgEqSolvers/LinAlgDirect.h",
    "mml/core/LinAlgEqSolvers/LinAlgQR.h",
    "mml/core/LinAlgEqSolvers/LinAlgSVD.h",
    "mml/core/LinAlgEqSolvers/LinAlgComplexSVD.h",
    "mml/core/LinAlgEqSolvers/LinAlgEqSolvers_iterative.h",

    # Sparse linear solvers - Individual components (in SparseSolvers/ subfolder)
    "mml/core/SparseSolvers/IterativeSolverBase.h",
    "mml/core/SparseSolvers/ConjugateGradient.h",
    "mml/core/SparseSolvers/BiCGSTAB.h",
    "mml/core/SparseSolvers/GMRES.h",
    "mml/core/SparseSolvers/Preconditioners.h",
    "mml/core/SparseSolvers/IterativeSolvers.h",

    # Vector spaces - Individual components (in VectorSpaces/ subfolder)
    "mml/core/VectorSpaces/Basis.h",
    "mml/core/VectorSpaces/Subspace.h",
    "mml/core/VectorSpaces/LinearMap.h",
    "mml/core/VectorSpaces/DualSpace.h",
    "mml/core/VectorSpaces/InnerProductSpace.h",
    "mml/core/VectorSpaces/AffineSpace.h",
    # [AGGREGATE] "mml/core/VectorSpaces.h",  # Includes VectorSpaces/* headers above

    # Richardson extrapolation and coordinate singularities
    "mml/base/RichardsonExtrapolation.h",
    "mml/core/DifferentialGeometry/CoordinateSingularities.h",

    # Derivation - Individual components (in Derivation/ subfolder)
    # [AGGREGATE] "mml/core/Derivation.h",  # Only includes Derivation/* headers below
    "mml/core/Derivation/DerivationBase.h",
    "mml/core/Derivation/FirstDerivativeStencil.h",
    "mml/core/Derivation/DerivationRealFunction.h",
    "mml/core/Derivation/DerivationScalarFunction.h",
    "mml/core/Derivation/DerivationVectorFunction.h",
    "mml/core/Derivation/DerivationParametricCurve.h",
    "mml/core/Derivation/DerivationParametricSurface.h",
    "mml/core/Derivation/DerivationTensorField.h",
    "mml/core/Derivation/Jacobians.h",
    "mml/core/Derivation/Hessians.h",
    "mml/core/Derivation/DerivationComplexStep.h",
    "mml/core/Derivation/DerivationComplex.h",
    "mml/core/Derivation/ForwardAD.h",
    "mml/core/Derivation/ReverseAD.h",
    "mml/core/Derivation/ADJacobians.h",

    # Integration - Individual components (in Integration/ subfolder)
    # [AGGREGATE] "mml/core/Integration.h",  # Only includes Integration/* headers below
    # IntegrationBase MUST come first - defines IntegrationResult
    # GaussKronrod.h MUST come before Integration1D.h because Integration1D.h
    # calls Integration::IntegrateGK21 in its IntegrateGK21 wrapper function
    "mml/core/Integration/IntegrationBase.h",
    "mml/core/Integration/GaussKronrod.h",
    "mml/core/Integration/Integration1D.h",
    "mml/core/Integration/Integration2D.h",
    "mml/core/Integration/Integration3D.h",
    "mml/core/Integration/GaussianQuadrature.h",
    "mml/core/Integration/IntegrationImproper.h",
    "mml/core/Integration/Integration2DAdaptive.h",
    "mml/core/Integration/Integration3DAdaptive.h",
    "mml/core/Integration/PathIntegration.h",
    "mml/core/Integration/SurfaceIntegration.h",
    "mml/core/Integration/MonteCarloIntegration.h",

    # Function helpers and fields
    "mml/core/FunctionHelpers.h",
    "mml/core/CoordTransf/CoordTransfBase.h",
    "mml/core/MetricTensor.h",
    "mml/core/DifferentialGeometry/MetricAdapters.h",
    "mml/core/Fields/FieldOperationsCommon.h",
    "mml/core/Fields/ScalarFieldOperations.h",
    "mml/core/Fields/VectorFieldOperations.h",
    "mml/core/DifferentialGeometry/FieldOperations.h",
    "mml/core/Fields/Fields.h",

    # Coordinate transformations - CoordTransfBase.h is emitted before MetricTensor.h
    "mml/core/CoordTransf/CoordTransf2D.h",
    "mml/core/CoordTransf/CoordTransf3D.h",
    "mml/core/CoordTransf/CoordTransfSpherical.h",
    "mml/core/CoordTransf/CoordTransfCylindrical.h",
    "mml/core/DifferentialGeometry/CoordinateMap.h",

    # Field operations, curves, surfaces
    "mml/core/Curves.h",
    "mml/core/Surfaces.h",
    "mml/core/Fields/TensorFieldAdapters.h",

    # Function spaces and bases
    "mml/core/OrthogonalBasis.h",
    "mml/core/OrthogonalBasis/LegendreBasis.h",
    "mml/core/OrthogonalBasis/ChebyshevBasis.h",
    "mml/core/OrthogonalBasis/HermiteBasis.h",
    "mml/core/OrthogonalBasis/LaguerreBasis.h",

    # Function spaces - Individual components (in FunctionSpaces/ subfolder)
    "mml/core/FunctionSpaces/FunctionSpacesBase.h",
    "mml/core/FunctionSpaces/FunctionSpace1D.h",
    "mml/core/FunctionSpaces/TrialSpace1D.h",
    "mml/core/FunctionSpaces/FunctionExpansion1D.h",
    "mml/core/FunctionSpaces/FunctionSpaceResult.h",
    "mml/core/FunctionSpaces/OrthogonalBasisTrialSpace1D.h",
    "mml/core/FunctionSpaces/Projection.h",
    "mml/core/FunctionSpaces/Interpolation.h",
    "mml/core/FunctionSpaces/ChebyshevCollocationSpace1D.h",
    "mml/core/FunctionSpaces/LinearDifferentialOperator1D.h",
    "mml/core/FunctionSpaces/BoundaryCondition1D.h",
    "mml/core/FunctionSpaces/OperatorAssembly.h",
    "mml/core/FunctionSpaces/DenseBVPSolver1D.h",
    "mml/core/FunctionSpaces/BVPDiagnostics1D.h",
    "mml/core/FunctionSpaces/LinearOperator.h",
    "mml/core/FunctionSpaces/MatrixFreeOperators1D.h",
    # [AGGREGATE] "mml/core/FunctionSpaces.h",  # Includes FunctionSpaces/* headers above

    # Complex analysis (contour integration, residues, argument principle)
    "mml/core/ComplexAnalysis.h",

    # =========================================================================
    # ALGORITHMS - Numerical algorithms
    # =========================================================================
    # Matrix algorithms
    # [AGGREGATE] "mml/algorithms/EigenSystemSolvers.h",  # Includes Eigen/* headers below
    "mml/algorithms/Eigen/EigenSolverConfig.h",
    "mml/algorithms/Eigen/HessenbergReduction.h",
    "mml/algorithms/Eigen/detail/RealSchurAnalysis.h",
    "mml/algorithms/Eigen/detail/HessenbergQRIteration.h",
    "mml/algorithms/Eigen/detail/SchurEigenvectors.h",
    "mml/algorithms/Eigen/SymmMatEigenSolverJacobi.h",
    "mml/algorithms/Eigen/HermitianMatEigenSolverJacobi.h",
    "mml/algorithms/Eigen/SymmMatEigenSolverQR.h",
    "mml/algorithms/Eigen/ComplexEigenSolver.h",
    "mml/algorithms/Eigen/EigenSolver.h",
    "mml/algorithms/MatrixAnalysisTypes.h",
    # [AGGREGATE] "mml/algorithms/MatrixAlg.h",  # Includes MatrixAlg/* headers below
    "mml/algorithms/MatrixAlg/Properties.h",
    "mml/algorithms/MatrixAlg/Measures.h",
    "mml/algorithms/MatrixAlg/Algebra.h",
    "mml/algorithms/MatrixAlg/Decompositions.h",
    "mml/algorithms/MatrixAlg/Subspaces.h",
    "mml/algorithms/MatrixAlg/Eigensystems.h",
    "mml/algorithms/MatrixAlg/Definiteness.h",

    # Curve fitting and interpolation
    "mml/algorithms/CurveFitting.h",

    # Root finding - Individual components (in RootFinding/ subfolder)
    # [AGGREGATE] "mml/algorithms/RootFinding.h",  # Only includes RootFinding/* headers below
    "mml/algorithms/RootFinding/RootFindingBase.h",
    "mml/algorithms/RootFinding/RootIsolation.h",
    "mml/algorithms/RootFinding/RootFindingBracketing.h",
    "mml/algorithms/RootFinding/RootFindingMethods.h",
    "mml/algorithms/RootFinding/RootFindingAllRoots.h",
    "mml/algorithms/RootFinding/NonlinearSystemSolvers.h",
    "mml/algorithms/RootFinding/RootFindingPolynoms.h",
    "mml/algorithms/RootFinding/RootFindingComplex.h",

    # ODE Solvers - Individual components (in ODESolvers/ subfolder)
    # [AGGREGATE] "mml/algorithms/ODESolvers.h",  # Only includes ODESolvers/* headers below
    "mml/algorithms/ODESolvers/ODERKCoefficients.h",  # Must come before shared stepper infrastructure
    "mml/algorithms/ODESolvers/ODEStepperInfrastructure.h",
    "mml/algorithms/ODESolvers/ODEStepCalculators.h",
    "mml/algorithms/ODESolvers/ODESteppers.h",
    "mml/algorithms/ODESolvers/ODESolverFixedStep.h",
    "mml/algorithms/ODESolvers/ODESolverAdaptive.h",
    "mml/algorithms/ODESolvers/ODESolverStiff.h",
    "mml/algorithms/ODESolvers/ODESolverEventDetection.h",
    "mml/algorithms/ODESolvers/BVPShootingMethod.h",

    # DAE solvers (in DAESolvers/ subfolder)
    # [AGGREGATE] "mml/algorithms/DAESolvers.h",  # Only includes DAESolvers/* headers below
    "mml/algorithms/DAESolvers/DAESolverBase.h",
    "mml/algorithms/DAESolvers/DAENumericalJacobian.h",
    "mml/algorithms/DAESolvers/DAEBackwardEuler.h",
    "mml/algorithms/DAESolvers/DAEAdaptive.h",
    "mml/algorithms/DAESolvers/DAEBDF2.h",
    "mml/algorithms/DAESolvers/DAEBDF4.h",
    "mml/algorithms/DAESolvers/DAERODAS.h",
    "mml/algorithms/DAESolvers/DAERadauIIA.h",
    "mml/algorithms/DAESolvers/DAEEventDetection.h",

    # Geodesic helpers depend on metric tensors and ODE fixed-step solvers
    "mml/core/DifferentialGeometry/InducedMetric.h",
    "mml/algorithms/Geodesic.h",

    # Statistics core
    "mml/algorithms/Statistics.h",
    # Statistics - Individual components (in Statistics/ subfolder)
    "mml/algorithms/Statistics/StatisticsBase.h",
    "mml/algorithms/Statistics/Histogram.h",
    "mml/algorithms/Statistics/Distributions.h",
    "mml/algorithms/Statistics/DiscreteDistributions.h",

    # Optimization core
    "mml/algorithms/Optimization/Optimization.h",
    "mml/algorithms/Optimization/Multidim/MultidimTypes.h",
    "mml/algorithms/Optimization/Multidim/NelderMead.h",
    "mml/algorithms/Optimization/Multidim/LineSearch.h",
    "mml/algorithms/Optimization/Multidim/Powell.h",
    "mml/algorithms/Optimization/Multidim/QuasiNewton.h",
    "mml/algorithms/Optimization/Multidim/MultidimSolvers.h",
    "mml/algorithms/Optimization/Constraints/BoundConstraints.h",
    "mml/algorithms/Optimization/Constraints/BoxConstrainedNelderMead.h",
    "mml/algorithms/Optimization/Constraints/BoxConstrainedPowell.h",
    "mml/algorithms/Optimization/Constraints/ProjectedGradient.h",
    "mml/algorithms/Optimization/Constraints/AugmentedLagrangian.h",
    "mml/algorithms/Optimization/LP/LPTypes.h",
    "mml/algorithms/Optimization/LP/LinearProgram.h",
    "mml/algorithms/Optimization/LP/SimplexTableau.h",
    "mml/algorithms/Optimization/LP/SimplexSolver.h",
    "mml/algorithms/Optimization/LP/LPSolvers.h",

    # Function analysis
    "mml/algorithms/Analyzers/FunctionsAnalyzer.h",
    "mml/algorithms/Analyzers/FieldAnalyzers.h",
    "mml/algorithms/Analyzers/MatrixAnalyzer.h",

    # Fourier analysis
    "mml/algorithms/Fourier/Fourier.h",
    "mml/algorithms/Fourier/FourierRealFFT.h",
    "mml/algorithms/Fourier/FourierWindowing.h",

    # Graph algorithms
    # [AGGREGATE] "mml/algorithms/GraphAlg.h"
    "mml/algorithms/Graphs/detail/GraphAlgorithmUtils.h",
    "mml/algorithms/Graphs/GraphTraversal.h",
    "mml/algorithms/Graphs/GraphShortestPaths.h",
    "mml/algorithms/Graphs/GraphStructure.h",
    "mml/algorithms/Graphs/GraphFlow.h",
    "mml/algorithms/Graphs/GraphSpanningTree.h",

    # Computational geometry - Individual components (in CompGeometry/ subfolder)
    # [AGGREGATE] "mml/algorithms/ComputationalGeometry.h",  # Only includes CompGeometry/* headers below
    "mml/algorithms/CompGeometry/CompGeometryBase.h",
    "mml/algorithms/CompGeometry/RobustPredicates.h",
    "mml/algorithms/CompGeometry/ConvexHull.h",
    "mml/algorithms/CompGeometry/ConvexHull3D.h",
    "mml/algorithms/CompGeometry/Triangulation.h",
    "mml/algorithms/CompGeometry/VoronoiDiagram.h",
    "mml/algorithms/CompGeometry/Circles.h",
    "mml/algorithms/CompGeometry/Intersections.h",
    "mml/algorithms/CompGeometry/KDTree.h",
    "mml/algorithms/CompGeometry/PolygonOps.h",

    # =========================================================================
    # TOOLS - Utility classes
    # =========================================================================
    "mml/tools/CsvUtils.h",
    "mml/tools/ConsolePrinter.h",
	"mml/tools/IOResult.h",

	# Persistence - Durable object save/load components
	# [AGGREGATE] "mml/tools/Persistence.h",  # Only includes persistence/* headers below
	"mml/tools/persistence/PersistenceBase.h",
	"mml/tools/persistence/JSON.h",
    "mml/tools/persistence/BinaryBase.h",
    "mml/tools/persistence/VectorBinary.h",
    "mml/tools/persistence/MatrixBinary.h",
	"mml/tools/persistence/VectorJSON.h",
	"mml/tools/persistence/MatrixJSON.h",
	"mml/tools/persistence/FunctionJSON.h",
	"mml/tools/persistence/MatrixIO.h",

    # Serializer - Visualizer presentation-data exporters
    # [AGGREGATE] "mml/tools/Serializer.h",  # Only includes serializer/* headers below
    "mml/tools/serializer/SerializerBase.h",
    "mml/tools/serializer/SerializerVectorFields.h",
    "mml/tools/serializer/SerializerFunctions.h",
    "mml/tools/serializer/SerializerCurves.h",
    "mml/tools/serializer/SerializerSurfaces.h",
    "mml/tools/serializer/SerializerODE.h",
    "mml/tools/serializer/SerializerSimulation.h",

    "mml/tools/Visualizer.h",
    "mml/tools/ThreadPool.h",
    "mml/tools/Timer.h",

    # DataLoader - Individual components (in data_loader/ subfolder)
    # [AGGREGATE] "mml/tools/DataLoader.h",  # Only includes data_loader/* headers below
    "mml/tools/data_loader/DataLoaderTypes.h",
    "mml/tools/data_loader/DataLoaderParsing.h",
    "mml/tools/data_loader/DataLoaderCSV.h",
    "mml/tools/data_loader/DataLoaderJSON.h",

    # =========================================================================
    # SYSTEMS - High-level mathematical systems
    # =========================================================================
    "mml/systems/LinearSystem.h",
    "mml/systems/DynamicalSystem/DynamicalSystemTypes.h",
    "mml/systems/DynamicalSystem/DynamicalSystemBase.h",
    "mml/systems/ContinuousSystems.h",
    "mml/systems/DiscreteMaps.h",
    # [AGGREGATE] "mml/systems/DynamicalSystem.h",
    "mml/systems/DynamicalSystem/DynamicalAnalysisCommon.h",
    "mml/systems/DynamicalSystem/FixedPointAnalysis.h",
    "mml/systems/DynamicalSystem/LyapunovAnalysis.h",
    "mml/systems/DynamicalSystem/BifurcationAnalysis.h",
    "mml/systems/DynamicalSystem/PhaseSpaceAnalysis.h",
    "mml/systems/DynamicalSystemAnalyzer.h",
]

std_headers_file = script_dir / "MMLStandardHeaders.h"
with open(std_headers_file, "r", encoding='utf-8') as fStdHeaders:
    stdHeaderText = fStdHeaders.read()

header_variants = [
    {
        "file_name": "MML.h",
        "guard_name": "MML_SINGLE_HEADER",
        "precision_macro": "MML_USE_DOUBLE",
    },
    {
        "file_name": "MML_float.h",
        "guard_name": "MML_SINGLE_HEADER_FLOAT",
        "precision_macro": "MML_USE_FLOAT",
    },
    {
        "file_name": "MML_long_double.h",
        "guard_name": "MML_SINGLE_HEADER_LONG_DOUBLE",
        "precision_macro": "MML_USE_LONG_DOUBLE",
    },
]

def build_standard_header_text(guard_name, precision_macro):
    text = stdHeaderText.replace("MML_SINGLE_HEADER", guard_name)

    guard_define = f"#define {guard_name}\n"
    variant_guard = (
        guard_define +
        "\n#ifndef MML_SINGLE_HEADER_VARIANT_INCLUDED\n" +
        "#define MML_SINGLE_HEADER_VARIANT_INCLUDED\n" +
        "#else\n" +
        "#error Only one MML single-header variant can be included in one translation unit.\n" +
        "#endif\n"
    )
    if precision_macro:
        variant_guard += (
            "\n#undef MML_USE_FLOAT\n" +
            "#undef MML_USE_DOUBLE\n" +
            "#undef MML_USE_LONG_DOUBLE\n" +
            f"#define {precision_macro}\n"
        )

    return text.replace(guard_define, variant_guard, 1).rstrip() + '\n'

def strip_header_comment(content):
    """
    Remove the MML header comment block from the start of a file.
    The header block consists of lines starting with '///' at the beginning of the file.
    Only removes the initial header block - other comments in the file are preserved.
    """
    lines = content.splitlines(keepends=True)

    # Skip empty lines at the start
    start_idx = 0
    while start_idx < len(lines) and lines[start_idx].strip() == '':
        start_idx += 1

    # Check if we have a header block (lines starting with ///)
    if start_idx < len(lines) and lines[start_idx].strip().startswith('///'):
        # Skip all consecutive /// comment lines
        while start_idx < len(lines) and lines[start_idx].strip().startswith('///'):
            start_idx += 1

    # Return the remaining content
    return ''.join(lines[start_idx:])

def should_skip_preprocessor_line(line):
    """
    Determine if a preprocessor line should be skipped in the single-header output.

    Skip:
    - #include directives (consolidated at top of MML.h)
    - #pragma once (single header has its own guard)
    - Include guards in any form:
      - #ifndef MML_..._H / #define MML_..._H
      - #if !defined MML_..._H / #if !defined(MML_..._H)

    Keep (important for conditional compilation):
    - #ifdef / #if / #if defined / #elif (non-guard)
    - #else
    - #endif (except for include guard endings)
    - #define (except include guard defines)
    - #undef
    - #error / #warning
    - #pragma (except 'once')
    """
    stripped = line.strip()

    # Always skip #include
    if stripped.startswith("#include"):
        return True

    # Always skip #pragma once
    if stripped.startswith("#pragma") and "once" in stripped:
        return True

    # Skip include guards: #ifndef MML_..._H or #ifndef _MML_..._H
    if stripped.startswith("#ifndef"):
        guard_name = stripped[7:].strip()
        if guard_name.startswith("MML_") or guard_name.startswith("_MML_"):
            if guard_name.endswith("_H") or guard_name.endswith("_H_"):
                return True

    # Skip include guards: #if !defined MML_..._H or #if !defined(MML_..._H)
    if stripped.startswith("#if"):
        # Check for "#if !defined" or "#if ! defined" patterns
        rest = stripped[3:].strip()
        if rest.startswith("!"):
            rest = rest[1:].strip()
            if rest.startswith("defined"):
                rest = rest[7:].strip()
                # Handle both "defined MML_..." and "defined(MML_...)"
                if rest.startswith("("):
                    rest = rest[1:].strip()
                    if rest.endswith(")"):
                        rest = rest[:-1].strip()
                guard_name = rest
                if guard_name.startswith("MML_") or guard_name.startswith("_MML_"):
                    if guard_name.endswith("_H") or guard_name.endswith("_H_"):
                        return True

    # Skip include guard defines: #define MML_..._H
    if stripped.startswith("#define"):
        parts = stripped.split()
        if len(parts) >= 2:
            define_name = parts[1]
            if (define_name.startswith("MML_") or define_name.startswith("_MML_")):
                if define_name.endswith("_H") or define_name.endswith("_H_"):
                    return True

    # Keep all other preprocessor directives
    return False


def write_single_header(output_file, guard_name, precision_macro):
    with open(output_file, "w", encoding='utf-8') as fSingleHeaderFile:
        fSingleHeaderFile.write(build_standard_header_text(guard_name, precision_macro))

        for fileName in listFiles :
            file_path = project_root / fileName
            # Try UTF-8 first, fallback to latin-1 if that fails
            try:
                with open(file_path,'r', encoding='utf-8-sig') as f:  # utf-8-sig strips BOM
                    content = f.read()
            except UnicodeDecodeError:
                with open(file_path,'r', encoding='latin-1') as f:
                    content = f.read()

            # Remove the header comment block from the start of the file
            content = strip_header_comment(content)

            # dodati liniju komentara na početku
            fSingleHeaderFile.write("\n///////////////////////////   " + fileName + "   ///////////////////////////\n")

            # Track if we're at the very end of file (for skipping trailing #endif of include guard)
            lines = content.splitlines(keepends=True)

            # Find the last non-empty line index
            last_content_idx = len(lines) - 1
            while last_content_idx >= 0 and lines[last_content_idx].strip() == '':
                last_content_idx -= 1

            output_lines = []
            skipped_preprocessor_line = False
            for idx, line in enumerate(lines):
                stripped = line.strip()
                if stripped.startswith("#"):
                    # Check if this is the final #endif (include guard closing)
                    if idx == last_content_idx and stripped == "#endif":
                        # Skip trailing include guard #endif
                        skipped_preprocessor_line = True
                        continue
                    # Also check for #endif with comment like "#endif // MML_..."
                    if idx == last_content_idx and stripped.startswith("#endif"):
                        guard_comment = stripped[6:].strip()
                        if guard_comment.startswith("//") and ("MML_" in guard_comment or "_H" in guard_comment):
                            skipped_preprocessor_line = True
                            continue

                    if should_skip_preprocessor_line(line):
                        skipped_preprocessor_line = True
                        continue

                if skipped_preprocessor_line and stripped == '':
                    continue

                line_without_eol = line.rstrip('\r\n')
                has_eol = len(line_without_eol) != len(line)
                output_lines.append(line_without_eol.rstrip(' \t') + ('\n' if has_eol else ''))
                skipped_preprocessor_line = False

            while output_lines and output_lines[0].strip() == '':
                output_lines.pop(0)
            while output_lines and output_lines[-1].strip() == '':
                output_lines.pop()

            for line in output_lines:
                line_without_eol = line.rstrip('\r\n')
                has_eol = len(line_without_eol) != len(line)
                fSingleHeaderFile.write(line_without_eol)
                if has_eol:
                    fSingleHeaderFile.write('\n')

        fSingleHeaderFile.write(f"\n#endif   //{guard_name}\n")

for variant in header_variants:
    output_file = project_root / "mml" / "single_header" / variant["file_name"]
    write_single_header(output_file, variant["guard_name"], variant["precision_macro"])
    print(f"Generated {output_file}")
