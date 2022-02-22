# add_library(fatrop SHARED FatropMemory.cpp LinearAlgebraBlasfeo.cpp InterfaceMUMPS.cpp)
set(FATROP_INTERNAL
blasfeo_wrapper/LinearAlgebraBlasfeo.hpp    blasfeo_wrapper/LinearAlgebraBlasfeo.cpp
function_evaluation/CasadiCodegen.hpp       function_evaluation/CasadiCodegen.cpp
aux/Aux.hpp                                 aux/Aux.cpp
aux/Common.hpp                              aux/Common.cpp
aux/DynamicLib.hpp                          aux/DynamicLib.cpp
aux/FatropVector.hpp                        aux/FatropVector.cpp
aux/LinearAlgebra.hpp                       aux/LinearAlgebra.cpp
aux/SmartPtr.hpp                            aux/SmartPtr.cpp
debug/FatropDebugTools.hpp                  debug/FatropDebugTools.cpp
debug/LinearAlgebraEigen.hpp                debug/LinearAlgebraEigen.cpp
debug/RandomOCP.hpp                         debug/RandomOCP.cpp
function_evaluation/FunctionEvaluation.hpp  function_evaluation/FunctionEvaluation.cpp
ocp/BFOCP.hpp                               ocp/BFOCP.cpp
ocp/BFOCPAdapter.hpp                        ocp/BFOCPAdapter.cpp
ocp/BFOCPBasic.hpp                          ocp/BFOCPBasic.cpp
ocp/DuInfEvaluator.hpp                      ocp/DuInfEvaluator.cpp
ocp/FatropOCP.hpp                           ocp/FatropOCP.cpp
ocp/FatropOCPBuilder.hpp                    ocp/FatropOCPBuilder.cpp
ocp/OCP.hpp                                 ocp/OCP.cpp
ocp/OCPDims.hpp                             ocp/OCPDims.cpp
ocp/OCPInitializer.hpp                      ocp/OCPInitializer.cpp
ocp/OCPKKT.hpp                              ocp/OCPKKT.cpp
ocp/OCPLinearSolver.hpp                     ocp/OCPLinearSolver.cpp
ocp/OCPLSRiccati.hpp                        ocp/OCPLSRiccati.cpp
ocp/OCPNoScaling.hpp                        ocp/OCPNoScaling.cpp
ocp/OCPScalingMethod.hpp                    ocp/OCPScalingMethod.cpp
solver/AlgorithmQuantities.hpp              solver/AlgorithmQuantities.cpp
solver/AlgStrategy.hpp                      solver/AlgStrategy.cpp
solver/FatropAlg.hpp                        solver/FatropAlg.cpp
solver/FatropData.hpp                       solver/FatropData.cpp
solver/FatropParams.hpp                     solver/FatropParams.cpp
solver/FatropPrinter.hpp                    solver/FatropPrinter.cpp
solver/Filter.hpp                           solver/Filter.cpp
solver/IterationData.hpp                    solver/IterationData.cpp
solver/LineSearch.hpp                       solver/LineSearch.cpp
solver/StepAcceptor.hpp                     solver/StepAcceptor.cpp
solver/Timer.hpp                            solver/Timer.cpp
sparse/FatropSparse.hpp                     sparse/FatropSparse.cpp
sparse/Interface.hpp                        sparse/Interface.cpp
sparse/InterfaceMA57.hpp                    sparse/InterfaceMA57.cpp
sparse/InterfaceMUMPS.hpp                   sparse/InterfaceMUMPS.cpp
sparse/InterfacePardiso.hpp                 sparse/InterfacePardiso.cpp
sparse/SparseOCP.hpp                        sparse/SparseOCP.cpp
templates/NLPAlg.hpp                        templates/NLPAlg.cpp
json/json.h json/json.cpp
)
set(FATROPHEADERS
blasfeo_wrapper/LinearAlgebraBlasfeo.hpp
function_evaluation/CasadiCodegen.hpp
aux/Aux.hpp
aux/Common.hpp
aux/DynamicLib.hpp
aux/FatropVector.hpp
aux/LinearAlgebra.hpp
aux/SmartPtr.hpp
debug/FatropDebugTools.hpp
debug/LinearAlgebraEigen.hpp
debug/RandomOCP.hpp
function_evaluation/FunctionEvaluation.hpp
ocp/BFOCP.hpp
ocp/BFOCPAdapter.hpp
ocp/BFOCPBasic.hpp
ocp/DuInfEvaluator.hpp
ocp/FatropOCP.hpp
ocp/FatropOCPBuilder.hpp
ocp/OCP.hpp
ocp/OCPDims.hpp
ocp/OCPInitializer.hpp
ocp/OCPKKT.hpp
ocp/OCPLinearSolver.hpp
ocp/OCPLSRiccati.hpp
ocp/OCPNoScaling.hpp
ocp/OCPScalingMethod.hpp
solver/AlgorithmQuantities.hpp
solver/AlgStrategy.hpp
solver/FatropAlg.hpp
solver/FatropData.hpp
solver/FatropParams.hpp
solver/FatropPrinter.hpp
solver/Filter.hpp
solver/IterationData.hpp
solver/LineSearch.hpp
solver/StepAcceptor.hpp
solver/Timer.hpp
sparse/FatropSparse.hpp
sparse/Interface.hpp
sparse/InterfaceMA57.hpp
sparse/InterfaceMUMPS.hpp
sparse/InterfacePardiso.hpp
sparse/SparseOCP.hpp
templates/NLPAlg.hpp)
add_library(fatrop SHARED ${FATROP_INTERNAL})
target_link_libraries(fatrop m dl blasfeo) 
find_package(MUMPS)
target_link_libraries(fatrop MUMPS::MPISEQ)
target_link_libraries(fatrop MUMPS::DMUMPS)
target_include_directories(fatrop
    PRIVATE
        ${PROJECT_SOURCE_DIR}/fatrop
)