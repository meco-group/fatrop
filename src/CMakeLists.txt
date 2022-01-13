# add_library(fatrop SHARED FatropMemory.cpp LinearAlgebraBlasfeo.cpp InterfaceMUMPS.cpp)
add_library(fatrop SHARED FatropMemory.cpp LinearAlgebraBlasfeo.cpp CasadiCodegen.cpp)
target_link_libraries(fatrop m dl blasfeo) 
target_include_directories(fatrop
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)


# target_link_libraries(fatrop esmumps) 
# target_include_directories(fatrop
#     PRIVATE
#         /home/lander/installations/MUMPS_5.4.1/include
# )