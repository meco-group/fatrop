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
# debug/FatropDebugTools.hpp                  debug/FatropDebugTools.cpp
# debug/LinearAlgebraEigen.hpp                debug/LinearAlgebraEigen.cpp
# debug/RandomOCP.hpp                         debug/RandomOCP.cpp
function_evaluation/FunctionEvaluation.hpp  function_evaluation/FunctionEvaluation.cpp
ocp/OCPAbstract.hpp                         ocp/OCPAbstract.cpp
ocp/BasicOCPApplication.hpp
# ocp/BFOCPAL.hpp                             ocp/BFOCPAL.cpp
ocp/BFOCPAdapter.hpp                        ocp/BFOCPAdapter.cpp
# ocp/BFOCPAdapterAL.hpp                      ocp/BFOCPAdapterAL.cpp
ocp/BFOCPBasic.hpp                          ocp/BFOCPBasic.cpp
# ocp/OCPBuilder.hpp                          ocp/OCPBuilder.cpp
# ocp/OCPBuilderAL.hpp                        ocp/OCPBuilderAL.cpp
ocp/DuInfEvaluator.hpp                      ocp/DuInfEvaluator.cpp
ocp/FatropOCP.hpp                           ocp/FatropOCP.cpp
ocp/BasicOCPSamplers.hpp                    ocp/BasicOCPSamplers.cpp
# ocp/FatropOCPAL.hpp                         ocp/FatropOCPAL.cpp
ocp/FatropOCPBuilder.hpp                    ocp/FatropOCPBuilder.cpp
ocp/OCP.hpp                                 ocp/OCP.cpp
# ocp/OCPAL.hpp                               ocp/OCPAL.cpp
ocp/OCPDims.hpp                             ocp/OCPDims.cpp
ocp/OCPInitializer.hpp                      ocp/OCPInitializer.cpp
ocp/OCPKKT.hpp                              ocp/OCPKKT.cpp
ocp/OCPLinearSolver.hpp                     ocp/OCPLinearSolver.cpp
ocp/OCPLSRiccati.hpp                        ocp/OCPLSRiccati.cpp
ocp/OCPNoScaling.hpp                        ocp/OCPNoScaling.cpp
ocp/OCPScalingMethod.hpp                    ocp/OCPScalingMethod.cpp
# ocp/LineSearchDDP.hpp                       ocp/LineSearchDDP.cpp
solver/AlgBuilder.hpp              #solver/AlgorithmQuantities.cpp
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
solver/FatropStats.hpp                            
# alm_solver/FatropALMAlg.hpp                 alm_solver/FatropALMAlg.cpp
# alm_solver/FatropALMData.hpp                alm_solver/FatropALMData.cpp
# sparse/FatropSparse.hpp                     sparse/FatropSparse.cpp
# sparse/Interface.hpp                        sparse/Interface.cpp
# sparse/InterfaceMA57.hpp                    sparse/InterfaceMA57.cpp
# sparse/InterfaceMUMPS.hpp                   sparse/InterfaceMUMPS.cpp
# sparse/InterfacePardiso.hpp                 sparse/InterfacePardiso.cpp
# sparse/SparseOCP.hpp                        sparse/SparseOCP.cpp
templates/NLPAlg.hpp                        templates/NLPAlg.cpp
# templates/FatropNLPAL.hpp                   templates/FatropNLPAL.cpp
templates/FatropApplication.hpp             templates/FatropApplication.cpp
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
# debug/FatropDebugTools.hpp
# debug/LinearAlgebraEigen.hpp
# debug/RandomOCP.hpp
function_evaluation/FunctionEvaluation.hpp
ocp/BasicOCPApplication.hpp
ocp/OCPAbstract.hpp
ocp/BFOCPAdapter.hpp
ocp/BFOCPBasic.hpp
ocp/DuInfEvaluator.hpp
ocp/FatropOCP.hpp
# ocp/FatropOCPAL.hpp
ocp/FatropOCPBuilder.hpp
ocp/BasicOCPSamplers.hpp
ocp/OCP.hpp
# ocp/OCPAL.hpp
ocp/OCPDims.hpp
ocp/OCPInitializer.hpp
ocp/OCPKKT.hpp
ocp/OCPLinearSolver.hpp
ocp/OCPLSRiccati.hpp
ocp/OCPNoScaling.hpp
ocp/OCPScalingMethod.hpp
# ocp/LineSearchDDP.hpp
solver/AlgorithmQuantities.hpp
solver/AlgBuilder.hpp
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
solver/FatropStats.hpp                            
# sparse/FatropSparse.hpp
# sparse/Interface.hpp
# sparse/InterfaceMA57.hpp
# sparse/InterfaceMUMPS.hpp
# sparse/InterfacePardiso.hpp
# sparse/SparseOCP.hpp
templates/FatropApplication.hpp
templates/NLPAlg.hpp)
add_library(fatrop SHARED ${FATROP_INTERNAL})
target_link_libraries(fatrop m dl blasfeo) 

# find_package(MUMPS)
# target_link_libraries(fatrop MUMPS::MPISEQ)
# target_link_libraries(fatrop MUMPS::DMUMPS)
target_include_directories(fatrop
    PRIVATE
        ${PROJECT_SOURCE_DIR}/fatrop
)

add_executable(RunFatrop executables/RunFatrop.cpp)
target_link_libraries(RunFatrop fatrop) 
target_include_directories(RunFatrop
    PRIVATE
        ${PROJECT_SOURCE_DIR}/fatrop
)

# set_target_properties(fatrop PROPERTIES PUBLIC_HEADER ${FATROPHEADERS})
# set_target_properties(fatrop PROPERTIES PUBLIC_HEADER "${FATROPHEADERS}")
list(APPEND INCLUDEDIRSFATROP
.
# aux
# blasfeo_wrapper
# debug
# function_evaluation
# json
# ocp
# solver
# sparse
# templates)
)
foreach(incl IN LISTS INCLUDEDIRSFATROP) 
install(DIRECTORY "${CMAKE_SOURCE_DIR}/fatrop/${incl}" # source directory
        DESTINATION "include/fatrop/${incl}" # target directory
        FILES_MATCHING # install only matched files
        PATTERN "*.hpp" # select C++ header files
        PATTERN "*.h" # select C header files
)
endforeach()
#install(TARGETS fatrop
#PUBLIC_HEADER DESTINATION fatrop)

set(INCLUDE_INSTALL_DIR "include/fatrop")
set(LIBRARY_INSTALL_DIR "lib")
set(RUNTIME_INSTALL_DIR "bin")
set(SYSCONFIG_INSTALL_DIR "etc/fatrop")

# Create the distribution package(s)
set(CPACK_PACKAGE_VERSION ${VERSION})
set(CPACK_INCLUDE_TOPLEVEL_DIRECTORY 0)
set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})

set(CPACK_PACKAGE_NAME "fatrop")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")

INSTALL(TARGETS fatrop blasfeo RunFatrop
  EXPORT fatropTargets
  LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR}
  ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR}
  INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR}
  RUNTIME DESTINATION ${RUNTIME_INSTALL_DIR}
  )

include(CMakePackageConfigHelpers)
set(ConfigFileInstallDir cmake/fatrop)
# to do: make the install dirs not-hard coded
# set(INCLUDE_INSTALL_DIR ${INCLUDE_INSTALL_DIR} CACHE PATH "install path for include files")
# set(LIBRARY_INSTALL_DIR ${LIBRARY_INSTALL_DIR} CACHE PATH "install path for libraries")
# set(RUNTIME_INSTALL_DIR ${RUNTIME_INSTALL_DIR} CACHE PATH "install path for runtimes")
configure_package_config_file(fatropConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/fatropConfig.cmake"
  INSTALL_DESTINATION "${ConfigFileInstallDir}"
  PATH_VARS INCLUDE_INSTALL_DIR LIBRARY_INSTALL_DIR
  )
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/fatropConfigVersion.cmake"
  VERSION "${VERSION}"
  COMPATIBILITY SameMajorVersion)

EXPORT(EXPORT fatropTargets
  FILE fatropTargets.cmake)

INSTALL(FILES
  "${CMAKE_CURRENT_BINARY_DIR}/fatropConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/fatropConfigVersion.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/fatropTargets.cmake"
  DESTINATION "${ConfigFileInstallDir}")

include(CPack)
