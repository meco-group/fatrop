#[[
Fatrop - A fast trajectory optimization solver
Copyright (C) 2022, 2023 Lander Vanroye, KU Leuven. All rights reserved.

This file is part of Fatrop.

Fatrop is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Fatrop is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Fatrop.  If not, see <http://www.gnu.org/licenses/>.]]
set(FATROP_INTERNAL
blasfeo_wrapper/LinearAlgebraBlasfeo.hpp    blasfeo_wrapper/LinearAlgebraBlasfeo.cpp
function_evaluation/CasadiCodegen.hpp       function_evaluation/CasadiCodegen.cpp
auxiliary/VectorUtils.hpp                   auxiliary/VectorUtils.cpp
auxiliary/Common.hpp                        auxiliary/Common.cpp
auxiliary/DynamicLib.hpp                    auxiliary/DynamicLib.cpp
auxiliary/FatropVector.hpp                  auxiliary/FatropVector.cpp
auxiliary/LinearAlgebra.hpp                 auxiliary/LinearAlgebra.cpp
function_evaluation/FunctionEvaluation.hpp  function_evaluation/FunctionEvaluation.cpp
ocp/OCPAbstract.hpp                         ocp/OCPAbstract.cpp
ocp/StageOCPApplication.hpp                 ocp/StageOCPApplication.cpp
ocp/OCPAdapter.hpp                          ocp/OCPAdapter.cpp
ocp/StageOCP.hpp                            ocp/StageOCP.cpp
ocp/DuInfEvaluator.hpp                      ocp/DuInfEvaluator.cpp
ocp/FatropOCP.hpp                           ocp/FatropOCP.cpp
ocp/StageOCPExpressions.hpp                 ocp/StageOCPExpressions.cpp
ocp/FatropOCPBuilder.hpp                    ocp/FatropOCPBuilder.cpp
ocp/OCP.hpp                                 ocp/OCP.cpp
ocp/OCPDims.hpp                             ocp/OCPDims.cpp
ocp/OCPInitializer.hpp                      ocp/OCPInitializer.cpp
ocp/OCPKKT.hpp                              ocp/OCPKKT.cpp
ocp/OCPLinearSolver.hpp                     ocp/OCPLinearSolver.cpp
ocp/OCPLSRiccati.hpp                        ocp/OCPLSRiccati.cpp
ocp/OCPNoScaling.hpp                        ocp/OCPNoScaling.cpp
ocp/OCPScalingMethod.hpp                    ocp/OCPScalingMethod.cpp
solver/AlgBuilder.hpp              
solver/AlgStrategy.hpp                      solver/AlgStrategy.cpp
solver/FatropAlg.hpp                        solver/FatropAlg.cpp
solver/FatropData.hpp                       solver/FatropData.cpp
solver/FatropOptions.hpp                    solver/FatropOptions.cpp
solver/FatropPrinter.hpp                    solver/FatropPrinter.cpp
solver/Filter.hpp                           solver/Filter.cpp
solver/IterationData.hpp                    solver/IterationData.cpp
solver/LineSearch.hpp                       solver/LineSearch.cpp
solver/FatropStats.hpp                            
templates/NLPAlg.hpp                        templates/NLPAlg.cpp
json/json.h json/json.cpp
)
set(FATROPHEADERS
blasfeo_wrapper/LinearAlgebraBlasfeo.hpp
function_evaluation/CasadiCodegen.hpp
auxiliary/VectorUtils.hpp
auxiliary/Common.hpp
auxiliary/DynamicLib.hpp
auxiliary/FatropVector.hpp
auxiliary/LinearAlgebra.hpp
function_evaluation/FunctionEvaluation.hpp
ocp/StageOCPApplication.hpp
ocp/OCPAbstract.hpp
ocp/OCPAdapter.hpp
ocp/Stage.hpp
ocp/DuInfEvaluator.hpp
ocp/FatropOCP.hpp
ocp/FatropOCPBuilder.hpp
ocp/StageOCPExpressions.hpp
ocp/OCP.hpp
ocp/OCPDims.hpp
ocp/OCPInitializer.hpp
ocp/OCPKKT.hpp
ocp/OCPLinearSolver.hpp
ocp/OCPLSRiccati.hpp
ocp/OCPNoScaling.hpp
ocp/OCPScalingMethod.hpp
solver/AlgorithmQuantities.hpp
solver/AlgBuilder.hpp
solver/AlgStrategy.hpp
solver/FatropAlg.hpp
solver/FatropData.hpp
solver/FatropOptions.hpp
solver/FatropPrinter.hpp
solver/Filter.hpp
solver/IterationData.hpp
solver/LineSearch.hpp
solver/StepAcceptor.hpp
solver/Timer.hpp
solver/FatropStats.hpp                            
templates/FatropApplication.hpp
templates/NLPAlg.hpp)
# add_library(fatrop SHARED ${FATROP_INTERNAL})
add_library(fatrop SHARED ${FATROP_INTERNAL})

if (WIN32)
  target_link_libraries(fatrop dl blasfeo)
else()
  target_link_libraries(fatrop m dl blasfeo) 
endif (WIN32)


add_executable(RunFatrop executables/RunFatrop.cpp)
target_link_libraries(RunFatrop fatrop) 

target_include_directories(RunFatrop
    PRIVATE
        ${PROJECT_SOURCE_DIR}/fatrop
)

target_include_directories(fatrop PUBLIC  
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/fatrop>  
    $<INSTALL_INTERFACE:include/fatrop> # <prefix>/include/mylib
)

list(APPEND INCLUDEDIRSFATROP
.
)
foreach(incl IN LISTS INCLUDEDIRSFATROP) 
install(DIRECTORY "${CMAKE_SOURCE_DIR}/fatrop/${incl}" # source directory
        DESTINATION "include/fatrop/${incl}" # target directory
        FILES_MATCHING # install only matched files
        PATTERN "*.hpp" # select C++ header files
        PATTERN "*.h" # select C header files
)
endforeach()

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

INSTALL(TARGETS fatrop RunFatrop
  EXPORT fatropConfig
  LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR}
  ARCHIVE DESTINATION ${LIBRARY_INSTALL_DIR}
  INCLUDES DESTINATION ${INCLUDE_INSTALL_DIR}
  RUNTIME DESTINATION ${RUNTIME_INSTALL_DIR}
  )

INSTALL(EXPORT fatropConfig DESTINATION cmake)