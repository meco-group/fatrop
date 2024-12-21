# Define source files for the fatrop library
set(FATROP_SOURCES
    linear_algebra/lu_factorization.cpp
    # Add other source files here
)

# Create the fatrop shared library
add_library(fatrop SHARED ${FATROP_SOURCES})

# Link fatrop with blasfeo
target_link_libraries(fatrop PUBLIC blasfeo)