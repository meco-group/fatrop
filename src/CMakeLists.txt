add_library(fatrop SHARED FatropMemory.cpp FatropLinearAlgebraBlasfeo.cpp)
target_link_libraries(fatrop m dl blasfeo) 
target_include_directories(fatrop
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)