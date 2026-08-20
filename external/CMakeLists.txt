include(FetchContent)

option(WITH_BUILD_BLASFEO "Build BLASFEO from source" OFF)
if(WITH_BUILD_BLASFEO AND NOT FATROP_LA_BACKEND STREQUAL "eigen")
    FetchContent_Declare(
        blasfeo
        GIT_REPOSITORY https://github.com/giaf/blasfeo.git
        GIT_TAG master
    )
    FetchContent_MakeAvailable(blasfeo)
endif()
