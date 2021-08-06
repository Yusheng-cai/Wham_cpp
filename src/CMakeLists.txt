add_library(src SHARED)

target_sources(src PRIVATE
    "DataFileParser.cpp"
    "TimeSeries.cpp"
    "Wham.cpp"
)

set(src_HEADER
    "DataFileParser.h"
    "TimeSeries.h"
    "Wham.h"
)

target_include_directories(src
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

target_link_libraries(src PUBLIC tools parallel)

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(src 
    INCLUDE_SUBDIR src
    INCLUDES ${src_HEADER})