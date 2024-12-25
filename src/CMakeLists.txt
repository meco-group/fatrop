# Define source files for the fatrop library
set(FATROP_SOURCES
    linear_algebra/lu_factorization.cpp
    linear_algebra/vector.cpp
    ocp/dims.cpp
    ocp/jacobian.cpp
    ocp/hessian.cpp
    ocp/problem_info.cpp
)

# Create the fatrop shared library
add_library(fatrop SHARED ${FATROP_SOURCES})

# Link fatrop with blasfeo
target_link_libraries(fatrop PUBLIC blasfeo)
