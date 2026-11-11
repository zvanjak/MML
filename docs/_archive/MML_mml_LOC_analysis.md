# MML LOC Analysis: Files Over 300 Lines

Generated: 2026-07-18

Scope: all regular files under `mml/`.

Counting method: physical lines from `wc -l`, including blank lines and comments. The threshold is strictly greater than 300 LOC.

## Summary

| Metric | Count |
|--------|------:|
| Total files scanned | 250 |
| Total LOC under `mml/` | 185,979 |
| Files longer than 300 LOC | 128 |
| LOC in files longer than 300 LOC | 167,155 |
| Share of files longer than 300 LOC | 51.2% |
| Share of LOC in files longer than 300 LOC | 89.9% |
| Generated single-header LOC | 89,887 |
| Long-file LOC excluding generated single header | 77,268 |

## Notes

- [mml/single_header/MML.h](../mml/single_header/MML.h) dominates the report at 89,887 LOC. It is included because it is a regular file under `mml/`, but it should be treated as generated packaging output rather than a refactoring target.
- Excluding the single header, the largest concentration of long files is in [mml/algorithms](../mml/algorithms), [mml/core](../mml/core), and [mml/base](../mml/base).
- Most long files are header-only implementation units, which is expected for this library, but the top non-generated files are still good candidates for future modularization or documentation-focused review.

## Long-File LOC By Area

| Area | Files > 300 LOC | LOC in long files | Share of long-file LOC |
|------|----------------:|------------------:|-----------------------:|
| [mml/single_header](../mml/single_header) | 1 | 89,887 | 53.8% |
| [mml/algorithms](../mml/algorithms) | 35 | 25,819 | 15.4% |
| [mml/core](../mml/core) | 36 | 20,473 | 12.2% |
| [mml/base](../mml/base) | 34 | 19,507 | 11.7% |
| [mml/tools](../mml/tools) | 9 | 5,340 | 3.2% |
| [mml/systems](../mml/systems) | 6 | 3,595 | 2.2% |
| [mml/interfaces](../mml/interfaces) | 3 | 1,164 | 0.7% |
| mml root headers | 4 | 1,370 | 0.8% |

## Largest Non-Generated Files

| Rank | LOC | File |
|-----:|----:|------|
| 1 | 2,325 | [mml/algorithms/Optimization/LinearProgramming.h](../mml/algorithms/Optimization/LinearProgramming.h) |
| 2 | 2,068 | [mml/algorithms/Optimization/OptimizationMultidim.h](../mml/algorithms/Optimization/OptimizationMultidim.h) |
| 3 | 1,502 | [mml/algorithms/Statistics/Distributions.h](../mml/algorithms/Statistics/Distributions.h) |
| 4 | 1,387 | [mml/base/Tensor.h](../mml/base/Tensor.h) |
| 5 | 1,294 | [mml/base/Matrix/Matrix.h](../mml/base/Matrix/Matrix.h) |
| 6 | 1,277 | [mml/core/FieldOperations.h](../mml/core/FieldOperations.h) |
| 7 | 1,218 | [mml/algorithms/ODESolvers/ODESteppers.h](../mml/algorithms/ODESolvers/ODESteppers.h) |
| 8 | 1,213 | [mml/core/Derivation/DerivationRealFunction.h](../mml/core/Derivation/DerivationRealFunction.h) |
| 9 | 1,161 | [mml/tools/Visualizer.h](../mml/tools/Visualizer.h) |
| 10 | 1,139 | [mml/algorithms/FunctionsAnalyzer.h](../mml/algorithms/FunctionsAnalyzer.h) |

## All Files Longer Than 300 LOC

| LOC | File | Area |
|----:|------|------|
| 89,887 | [mml/single_header/MML.h](../mml/single_header/MML.h) | mml/single_header |
| 2,325 | [mml/algorithms/Optimization/LinearProgramming.h](../mml/algorithms/Optimization/LinearProgramming.h) | mml/algorithms |
| 2,068 | [mml/algorithms/Optimization/OptimizationMultidim.h](../mml/algorithms/Optimization/OptimizationMultidim.h) | mml/algorithms |
| 1,502 | [mml/algorithms/Statistics/Distributions.h](../mml/algorithms/Statistics/Distributions.h) | mml/algorithms |
| 1,387 | [mml/base/Tensor.h](../mml/base/Tensor.h) | mml/base |
| 1,294 | [mml/base/Matrix/Matrix.h](../mml/base/Matrix/Matrix.h) | mml/base |
| 1,277 | [mml/core/FieldOperations.h](../mml/core/FieldOperations.h) | mml/core |
| 1,218 | [mml/algorithms/ODESolvers/ODESteppers.h](../mml/algorithms/ODESolvers/ODESteppers.h) | mml/algorithms |
| 1,213 | [mml/core/Derivation/DerivationRealFunction.h](../mml/core/Derivation/DerivationRealFunction.h) | mml/core |
| 1,161 | [mml/tools/Visualizer.h](../mml/tools/Visualizer.h) | mml/tools |
| 1,139 | [mml/algorithms/FunctionsAnalyzer.h](../mml/algorithms/FunctionsAnalyzer.h) | mml/algorithms |
| 1,124 | [mml/core/LinAlgEqSolvers/LinAlgDirect.h](../mml/core/LinAlgEqSolvers/LinAlgDirect.h) | mml/core |
| 1,120 | [mml/algorithms/Optimization/Optimization.h](../mml/algorithms/Optimization/Optimization.h) | mml/algorithms |
| 1,118 | [mml/systems/LinearSystem.h](../mml/systems/LinearSystem.h) | mml/systems |
| 1,099 | [mml/algorithms/EigenSystemSolvers.h](../mml/algorithms/EigenSystemSolvers.h) | mml/algorithms |
| 1,077 | [mml/algorithms/EigenSolverHelpers.h](../mml/algorithms/EigenSolverHelpers.h) | mml/algorithms |
| 1,008 | [mml/base/InterpolatedFunctions/InterpolatedRealFunction.h](../mml/base/InterpolatedFunctions/InterpolatedRealFunction.h) | mml/base |
| 634 | [mml/base/Vector/VectorTypes3D.h](../../mml/base/Vector/VectorTypes3D.h) | mml/base |
| 301 | [mml/base/Vector/VectorTypes2D.h](../../mml/base/Vector/VectorTypes2D.h) | mml/base |
| 968 | [mml/core/MatrixUtils.h](../mml/core/MatrixUtils.h) | mml/core |
| 959 | [mml/tools/ConsolePrinter.h](../mml/tools/ConsolePrinter.h) | mml/tools |
| 906 | [mml/algorithms/Fourier/Fourier.h](../mml/algorithms/Fourier/Fourier.h) | mml/algorithms |
| 867 | [mml/base/Matrix/MatrixNM.h](../mml/base/Matrix/MatrixNM.h) | mml/base |
| 866 | [mml/base/Matrix/MatrixSym.h](../mml/base/Matrix/MatrixSym.h) | mml/base |
| 855 | [mml/base/Intervals.h](../mml/base/Intervals.h) | mml/base |
| 849 | [mml/tools/README_Serialization.md](../mml/tools/README_Serialization.md) | mml/tools |
| 832 | [mml/algorithms/RootFinding/RootFindingMethods.h](../mml/algorithms/RootFinding/RootFindingMethods.h) | mml/algorithms |
| 816 | [mml/algorithms/GraphAlgorithms.h](../mml/algorithms/GraphAlgorithms.h) | mml/algorithms |
| 814 | [mml/systems/DynamicalSystemAnalyzers.h](../mml/systems/DynamicalSystemAnalyzers.h) | mml/systems |
| 806 | [mml/base/Quaternions.h](../mml/base/Quaternions.h) | mml/base |
| 802 | [mml/core/CoordTransf/CoordTransf3D.h](../mml/core/CoordTransf/CoordTransf3D.h) | mml/core |
| 766 | [mml/algorithms/Statistics.h](../mml/algorithms/Statistics.h) | mml/algorithms |
| 751 | [mml/core/Integration/GaussianQuadrature.h](../mml/core/Integration/GaussianQuadrature.h) | mml/core |
| 739 | [mml/algorithms/ODESolvers/ODESolverStiff.h](../mml/algorithms/ODESolvers/ODESolverStiff.h) | mml/algorithms |
| 733 | [mml/core/Curves.h](../mml/core/Curves.h) | mml/core |
| 727 | [mml/algorithms/CompGeometry/KDTree.h](../mml/algorithms/CompGeometry/KDTree.h) | mml/algorithms |
| 714 | [mml/core/Integration/Integration1D.h](../mml/core/Integration/Integration1D.h) | mml/core |
| 709 | [mml/core/Integration/IntegrationImproper.h](../mml/core/Integration/IntegrationImproper.h) | mml/core |
| 708 | [mml/algorithms/Statistics/StatisticsHypothesis.h](../mml/algorithms/Statistics/StatisticsHypothesis.h) | mml/algorithms |
| 639 | [mml/core/Surfaces.h](../mml/core/Surfaces.h) | mml/core |
| 635 | [mml/algorithms/CurveFitting.h](../mml/algorithms/CurveFitting.h) | mml/algorithms |
| 634 | [mml/core/LinAlgEqSolvers/LinAlgEqSolvers_iterative.h](../mml/core/LinAlgEqSolvers/LinAlgEqSolvers_iterative.h) | mml/core |
| 627 | [mml/core/Derivation/DerivationScalarFunction.h](../mml/core/Derivation/DerivationScalarFunction.h) | mml/core |
| 626 | [mml/base/Vector/Vector.h](../mml/base/Vector/Vector.h) | mml/base |
| 624 | [mml/base/Polynom.h](../mml/base/Polynom.h) | mml/base |
| 617 | [mml/algorithms/ODESolvers/BVPShootingMethod.h](../mml/algorithms/ODESolvers/BVPShootingMethod.h) | mml/algorithms |
| 600 | [mml/base/Graph.h](../mml/base/Graph.h) | mml/base |
| 600 | [mml/core/MetricTensor.h](../mml/core/MetricTensor.h) | mml/core |
| 588 | [mml/base/Geometry/GeometryCore/Geometry2DShapes.h](../mml/base/Geometry/GeometryCore/Geometry2DShapes.h) | mml/base |
| 585 | [mml/core/Derivation/DerivationVectorFunction.h](../mml/core/Derivation/DerivationVectorFunction.h) | mml/core |
| 582 | [mml/algorithms/GraphSpectral.h](../mml/algorithms/GraphSpectral.h) | mml/algorithms |
| 575 | [mml/base/SparseMatrix/SparseMatrixCSR.h](../mml/base/SparseMatrix/SparseMatrixCSR.h) | mml/base |
| 558 | [mml/tools/data_loader/DataLoaderTypes.h](../mml/tools/data_loader/DataLoaderTypes.h) | mml/tools |
| 536 | [mml/core/CoordSystem.h](../mml/core/CoordSystem.h) | mml/core |
| 533 | [mml/algorithms/CompGeometry/ConvexHull3D.h](../mml/algorithms/CompGeometry/ConvexHull3D.h) | mml/algorithms |
| 530 | [mml/base/DAESystem.h](../mml/base/DAESystem.h) | mml/base |
| 521 | [mml/algorithms/RootFinding/RootFindingBracketing.h](../mml/algorithms/RootFinding/RootFindingBracketing.h) | mml/algorithms |
| 519 | [mml/core/Integration/GaussKronrod.h](../mml/core/Integration/GaussKronrod.h) | mml/core |
| 508 | [mml/algorithms/ODESolvers/ODEStepCalculators.h](../mml/algorithms/ODESolvers/ODEStepCalculators.h) | mml/algorithms |
| 508 | [mml/systems/DynamicalSystemAnalyzer.h](../mml/systems/DynamicalSystemAnalyzer.h) | mml/systems |
| 505 | [mml/interfaces/IFunction.h](../mml/interfaces/IFunction.h) | mml/interfaces |
| 495 | [mml/algorithms/ODESolvers/ODESolverAdaptive.h](../mml/algorithms/ODESolvers/ODESolverAdaptive.h) | mml/algorithms |
| 494 | [mml/core/Integration/MonteCarloIntegration.h](../mml/core/Integration/MonteCarloIntegration.h) | mml/core |
| 493 | [mml/algorithms/RootFinding/RootFindingPolynoms.h](../mml/algorithms/RootFinding/RootFindingPolynoms.h) | mml/algorithms |
| 492 | [mml/algorithms/CompGeometry/Triangulation.h](../mml/algorithms/CompGeometry/Triangulation.h) | mml/algorithms |
| 492 | [mml/core/LinAlgEqSolvers/LinAlgSVD.h](../mml/core/LinAlgEqSolvers/LinAlgSVD.h) | mml/core |
| 486 | [mml/base/Matrix/Matrix3D.h](../mml/base/Matrix/Matrix3D.h) | mml/base |
| 480 | [mml/base/Matrix/MatrixBandDiag.h](../mml/base/Matrix/MatrixBandDiag.h) | mml/base |
| 478 | [mml/core/Fields.h](../mml/core/Fields.h) | mml/core |
| 475 | [mml/core/Integration/Integration3DAdaptive.h](../mml/core/Integration/Integration3DAdaptive.h) | mml/core |
| 473 | [mml/algorithms/CompGeometry/VoronoiDiagram.h](../mml/algorithms/CompGeometry/VoronoiDiagram.h) | mml/algorithms |
| 471 | [mml/base/Geometry/Geometry2DCore/Geometry2DCircleBox.h](../mml/base/Geometry/Geometry2DCore/Geometry2DCircleBox.h) | mml/base |
| 465 | [mml/base/ChebyshevApproximation.h](../mml/base/ChebyshevApproximation.h) | mml/base |
| 460 | [mml/base/Geometry/Geometry2DCore/Geometry2DLines.h](../mml/base/Geometry/Geometry2DCore/Geometry2DLines.h) | mml/base |
| 458 | [mml/core/ComplexAnalysis.h](../mml/core/ComplexAnalysis.h) | mml/core |
| 453 | [mml/base/BaseUtils/MatrixOps.h](../mml/base/BaseUtils/MatrixOps.h) | mml/base |
| 449 | [mml/base/Geometry/Geometry3DBodiesCore/Geometry3DBounding.h](../mml/base/Geometry/Geometry3DBodiesCore/Geometry3DBounding.h) | mml/base |
| 442 | [mml/base/Vector/VectorN.h](../mml/base/Vector/VectorN.h) | mml/base |
| 441 | [mml/base/Geometry/Geometry3DCore/Geometry3DLines.h](../mml/base/Geometry/Geometry3DCore/Geometry3DLines.h) | mml/base |
| 440 | [mml/core/Integration/Integration2DAdaptive.h](../mml/core/Integration/Integration2DAdaptive.h) | mml/core |
| 432 | [mml/systems/ContinuousSystems.h](../mml/systems/ContinuousSystems.h) | mml/systems |
| 431 | [mml/tools/data_loader/DataLoaderCSV.h](../mml/tools/data_loader/DataLoaderCSV.h) | mml/tools |
| 418 | [mml/core/FunctionHelpers.h](../mml/core/FunctionHelpers.h) | mml/core |
| 415 | [mml/algorithms/FieldLineTracer.h](../mml/algorithms/FieldLineTracer.h) | mml/algorithms |
| 412 | [mml/core/Integration/SurfaceIntegration.h](../mml/core/Integration/SurfaceIntegration.h) | mml/core |
| 410 | [mml/base/Geometry/GeometryCore/GeometryPoints.h](../mml/base/Geometry/GeometryCore/GeometryPoints.h) | mml/base |
| 405 | [mml/base/Function.h](../mml/base/Function.h) | mml/base |
| 401 | [mml/base/AlgorithmTypes.h](../mml/base/AlgorithmTypes.h) | mml/base |
| 401 | [mml/core/CoordTransf.h](../mml/core/CoordTransf.h) | mml/core |
| 397 | [mml/core/SparseSolvers/Preconditioners.h](../mml/core/SparseSolvers/Preconditioners.h) | mml/core |
| 393 | [mml/core/Derivation/ForwardAD.h](../mml/core/Derivation/ForwardAD.h) | mml/core |
| 378 | [mml/core/SparseSolvers/GMRES.h](../mml/core/SparseSolvers/GMRES.h) | mml/core |
| 374 | [mml/base/Geometry/GeometryCore/Geometry3DShapes.h](../mml/base/Geometry/GeometryCore/Geometry3DShapes.h) | mml/base |
| 374 | [mml/core/Derivation/ReverseAD.h](../mml/core/Derivation/ReverseAD.h) | mml/core |
| 369 | [mml/core/Derivation/Jacobians.h](../mml/core/Derivation/Jacobians.h) | mml/core |
| 366 | [mml/base/DifferentialGeometry/DifferentialForm.h](../mml/base/DifferentialGeometry/DifferentialForm.h) | mml/base |
| 365 | [mml/systems/README.md](../mml/systems/README.md) | mml/systems |
| 364 | [mml/base/Geometry/Geometry2DCore/Geometry2DPolygon.h](../mml/base/Geometry/Geometry2DCore/Geometry2DPolygon.h) | mml/base |
| 364 | [mml/core/LinAlgEqSolvers/LinAlgQR.h](../mml/core/LinAlgEqSolvers/LinAlgQR.h) | mml/core |
| 360 | [mml/tools/serializer/SerializerSurfaces.h](../mml/tools/serializer/SerializerSurfaces.h) | mml/tools |
| 359 | [mml/MMLExceptions.h](../mml/MMLExceptions.h) | mml root |
| 358 | [mml/systems/DiscreteMaps.h](../mml/systems/DiscreteMaps.h) | mml/systems |
| 357 | [mml/tools/serializer/SerializerFunctions.h](../mml/tools/serializer/SerializerFunctions.h) | mml/tools |
| 356 | [mml/MMLVisualizators.h](../mml/MMLVisualizators.h) | mml root |
| 354 | [mml/tools/data_loader/DataLoaderJSON.h](../mml/tools/data_loader/DataLoaderJSON.h) | mml/tools |
| 350 | [mml/base/BaseUtils/MixedTypeOps.h](../mml/base/BaseUtils/MixedTypeOps.h) | mml/base |
| 349 | [mml/core/Derivation/DerivationTensorField.h](../mml/core/Derivation/DerivationTensorField.h) | mml/core |
| 342 | [mml/algorithms/CompGeometry/PolygonOps.h](../mml/algorithms/CompGeometry/PolygonOps.h) | mml/algorithms |
| 339 | [mml/base/Geometry/Geometry3DBodiesCore/Geometry3DTorusCylinder.h](../mml/base/Geometry/Geometry3DBodiesCore/Geometry3DTorusCylinder.h) | mml/base |
| 339 | [mml/interfaces/IODESystemDAE.h](../mml/interfaces/IODESystemDAE.h) | mml/interfaces |
| 337 | [mml/base/Geometry/Geometry3DBodiesCore/Geometry3DSpherePyramid.h](../mml/base/Geometry/Geometry3DBodiesCore/Geometry3DSpherePyramid.h) | mml/base |
| 336 | [mml/algorithms/Statistics/StatisticsConfidence.h](../mml/algorithms/Statistics/StatisticsConfidence.h) | mml/algorithms |
| 330 | [mml/MMLPrecision.h](../mml/MMLPrecision.h) | mml root |
| 330 | [mml/core/Derivation/DerivationParametricCurve.h](../mml/core/Derivation/DerivationParametricCurve.h) | mml/core |
| 327 | [mml/base/Geometry/Geometry3DCore/Geometry3DSurfaces.h](../mml/base/Geometry/Geometry3DCore/Geometry3DSurfaces.h) | mml/base |
| 325 | [mml/MMLBase.h](../mml/MMLBase.h) | mml root |
| 320 | [mml/base/Geometry/Geometry3DBodiesCore/Geometry3DBodyBase.h](../mml/base/Geometry/Geometry3DBodiesCore/Geometry3DBodyBase.h) | mml/base |
| 320 | [mml/interfaces/ITensor.h](../mml/interfaces/ITensor.h) | mml/interfaces |
| 319 | [mml/algorithms/RootFinding/RootFindingComplex.h](../mml/algorithms/RootFinding/RootFindingComplex.h) | mml/algorithms |
| 317 | [mml/algorithms/Statistics/StatisticsRank.h](../mml/algorithms/Statistics/StatisticsRank.h) | mml/algorithms |
| 315 | [mml/base/Matrix/MatrixTriDiag.h](../mml/base/Matrix/MatrixTriDiag.h) | mml/base |
| 313 | [mml/algorithms/DAESolvers/DAEBDF4.h](../mml/algorithms/DAESolvers/DAEBDF4.h) | mml/algorithms |
| 313 | [mml/core/Derivation/DerivationComplex.h](../mml/core/Derivation/DerivationComplex.h) | mml/core |
| 312 | [mml/algorithms/DAESolvers/DAENumericalJacobian.h](../mml/algorithms/DAESolvers/DAENumericalJacobian.h) | mml/algorithms |
| 311 | [mml/tools/ThreadPool.h](../mml/tools/ThreadPool.h) | mml/tools |
| 308 | [mml/base/Geometry/Geometry3DBodiesCore/Geometry3DCubes.h](../mml/base/Geometry/Geometry3DBodiesCore/Geometry3DCubes.h) | mml/base |
| 306 | [mml/core/VectorSpaces/AffineSpace.h](../mml/core/VectorSpaces/AffineSpace.h) | mml/core |
| 305 | [mml/algorithms/ODESolvers/ODERKCoefficients.h](../mml/algorithms/ODESolvers/ODERKCoefficients.h) | mml/algorithms |
| 304 | [mml/algorithms/DAESolvers/DAERadauIIA.h](../mml/algorithms/DAESolvers/DAERadauIIA.h) | mml/algorithms |