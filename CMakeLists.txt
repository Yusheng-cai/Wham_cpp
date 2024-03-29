cmake_minimum_required(VERSION 3.18.1)

# Include the cmake scripts directory
LIST( APPEND CMAKE_MODULE_PATH
  ${CMAKE_SOURCE_DIR}/cmake/
)

# create the new project
project(Wham)

# Require C++11
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# reference : https://stackoverflow.com/questions/8591762/ifdef-debug-with-cmake-independent-from-platform
message("build type = ${CMAKE_BUILD_TYPE}")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -Wall -g -DMY_DEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -ffast-math -march=native")

# define my own debug defintion 
# set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DMY_DEBUG")

# Requires OpenMP
find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Requires FFTW3
find_package(FFTW3)
if (FFTW3_FOUND)
  message(STATUS "FFTW3 library: ${FFTW3_LIBRARIES}")
endif()

# install the libraries as well as the executables in lib & bin folder in build/
# This has to be after project is defined
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")

enable_testing()

add_subdirectory(src)
add_subdirectory(tools)
add_subdirectory(parallel)
add_subdirectory(test)

add_executable(Wham main.cpp)

target_include_directories(Wham PUBLIC ${CMAKE_SOURCE_DIR})

target_link_libraries(Wham PUBLIC WHAMsrc WHAMtools WHAMparallel)

install(TARGETS Wham DESTINATION bin) 
