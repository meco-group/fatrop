# add_library(fatrop SHARED FatropMemory.cpp LinearAlgebraBlasfeo.cpp InterfaceMUMPS.cpp)
set(FATROP_INTERNAL
blasfeo_wrapper/LinearAlgebraBlasfeo.hpp  blasfeo_wrapper/LinearAlgebraBlasfeo.cpp
function_evaluation/CasadiCodegen.hpp     function_evaluation/CasadiCodegen.cpp
)
add_library(fatrop SHARED ${FATROP_INTERNAL})
target_link_libraries(fatrop m dl blasfeo) 
target_include_directories(fatrop
    PRIVATE
        ${PROJECT_SOURCE_DIR}/fatrop
)

