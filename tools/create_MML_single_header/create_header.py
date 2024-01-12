import os


print(os.getcwd())

listFiles = ["./include/MMLBase.h",
             
            "./include/utilities/Constants.h",
            "./include/utilities/StdFunctions.h",
            "./include/utilities/DataContainers.h",
            "./include/utilities/Matrix3D.h",

            "./include/interfaces/IInterval.h",
            "./include/interfaces/IAlgebra.h",
            "./include/interfaces/IVectorSpace.h",
            "./include/interfaces/ITensor.h",
            "./include/utilities/Intervals.h",
            
            "./include/base/Algebra.h",
            "./include/base/Geometry.h",
            "./include/base/Vector.h",
            "./include/base/VectorN.h",
            "./include/base/VectorTypes.h",
            "./include/base/Matrix.h",
            "./include/base/MatrixBandDiag.h",
            "./include/base/MatrixNM.h",
            "./include/base/MatrixSparse.h",
            "./include/base/MatrixSym.h",
            "./include/base/Tensor.h",
            "./include/base/BaseUtils.h",

            "./include/interfaces/IFunction.h",
            "./include/interfaces/ICoordTransf.h",
            "./include/interfaces/IODESystem.h",

            "./include/core/LinAlgEqSolvers.h",
            "./include/core/VectorSpace.h",
            "./include/core/Polynom.h",
            "./include/core/MatrixUtils.h",

            "./include/core/Derivation.h",
            "./include/core/Integration.h",
            "./include/core/Function.h",
            "./include/core/FunctionHelpers.h",
            "./include/core/Jacobians.h",
            "./include/core/InterpolatedFunction.h",
            "./include/core/LinearFunctional.h",
            "./include/core/LinearOperator.h",
            "./include/core/QuadraticForm.h",
            "./include/core/MetricTensor.h",
            "./include/core/CoordTransf.h",
            "./include/core/FieldOperations.h",
            "./include/core/ChebyshevApproximation.h",

            "./include/basic_types/FunctionSpace.h",
            "./include/basic_types/Curves.h",
            "./include/basic_types/Surfaces.h",
            "./include/basic_types/DiracDeltaFunction.h",
            "./include/basic_types/ODESystem.h",
            "./include/basic_types/Geometry2D.h",
            "./include/basic_types/Geometry3D.h",
            "./include/basic_types/CoordSystem.h",
            "./include/basic_types/Fields.h",
            "./include/basic_types/PermutationGroup.h",

            "./include/algorithms/PathIntegration.h",
            "./include/algorithms/SurfaceIntegration.h",
            "./include/algorithms/EigenSystemSolvers.h",
            "./include/algorithms/ODESystemSteppers.h",
            "./include/algorithms/ODESystemSteppers_Stiff.h",
            "./include/algorithms/ODESystemSolver.h",
            "./include/algorithms/DiffGeometryAlgorithms.h",
            "./include/algorithms/FunctionAnalyzer.h",
            "./include/algorithms/Fourier.h",
            "./include/algorithms/RootFinding.h",
            "./include/algorithms/RootFindingMultidim.h",
            "./include/algorithms/Statistics.h",

            "./include/systems/LinAlgEqSystem.h",
            "./include/systems/DiffEqSystem.h"]

fSingleHeaderFile = open("./include/single_header/MML.h", "w")

# definirati #ifdef -> prebaceno u standard headers
# fSingleHeaderFile.write("#ifndef MML_SINGLE_HEADER\n")
# fSingleHeaderFile.write("#define MML_SINGLE_HEADER\n\n")

# iskopirati standard headers 
fStdHeaders = open("./tools/create_MML_single_header/MMLStandardHeaders.h", "r")
stdHeaderLines = fStdHeaders.readlines()
for stdLine in stdHeaderLines:
    fSingleHeaderFile.write(stdLine)

for fileName in listFiles :
    f=open(fileName,'r', encoding='utf-8') 

    # dodati liniju komentara na početku
    fSingleHeaderFile.write("///////////////////////////   " + fileName + "   ///////////////////////////")

    # lines = f.readlines()
    # for line in lines:
    while (line := f.readline()) :
        if line.startswith("#") == False :
            #print (line)
            fSingleHeaderFile.write(line)

fSingleHeaderFile.write("\n#endif\n")
