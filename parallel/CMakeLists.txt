add_library(WHAMparallel SHARED)

target_sources(WHAMparallel PRIVATE
    "OpenMP.cpp"
)

set(WHAMPARALLEL_HEADER
    "MPI_environment.h"
    "OpenMP_buffer.h"
    "OpenMP.h")

target_include_directories(WHAMparallel PUBLIC 
                         "${CMAKE_SOURCE_DIR}")

# target_compile_options(parallel PRIVATE -fPIC)

message("Required package of OpenMP is ${OpenMP_CXX_FOUND}")

# resource: https://gitlab.kitware.com/cmake/cmake/-/issues/21818
if(OpenMP_CXX_FOUND)
    target_link_libraries(WHAMparallel PUBLIC OpenMP::OpenMP_CXX)
endif()

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(WHAMparallel 
    INCLUDE_SUBDIR WHAMparallel 
    INCLUDES ${WHAMPARALLEL_HEADER})
    