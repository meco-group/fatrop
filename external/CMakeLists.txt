include(FetchContent)

option(WITH_BUILD_BLASFEO "Build BLASFEO from source" OFF)
if(WITH_BUILD_BLASFEO)
    FetchContent_Declare(
        blasfeo
        GIT_REPOSITORY https://github.com/giaf/blasfeo.git
        GIT_TAG master
    )
    FetchContent_MakeAvailable(blasfeo)
endif()
