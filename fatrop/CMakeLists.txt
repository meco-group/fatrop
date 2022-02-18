# add_library(fatrop SHARED FatropMemory.cpp LinearAlgebraBlasfeo.cpp InterfaceMUMPS.cpp)
add_library(fatrop SHARED blasfeo_wrapper/LinearAlgebraBlasfeo.cpp function_evaluation/CasadiCodegen.cpp)
target_link_libraries(fatrop m dl blasfeo) 
target_include_directories(fatrop
    PRIVATE
        ${PROJECT_SOURCE_DIR}/fatrop
)

# target_link_libraries(fatrop esmumps) 
# target_include_directories(fatrop
#     PRIVATE
#         $ENV{HOME}/installations/MUMPS_5.4.1/include
# )