# Define legacy source files
set(LEGACY_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/OCPCInterface.cpp
    # Add other legacy source files here as needed
)

# Add legacy sources to the main fatrop target
target_sources(fatrop PRIVATE ${LEGACY_SOURCES})

# Add legacy include directory to the main fatrop target
target_include_directories(fatrop
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# Install header files (only .hpp)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
        DESTINATION ${INCLUDE_INSTALL_DIR}
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")