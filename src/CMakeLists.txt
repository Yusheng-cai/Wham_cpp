add_library(WHAMsrc SHARED)

target_sources(WHAMsrc PRIVATE
    "DataFileParser.cpp"
    "TimeSeries.cpp"
    "Wham.cpp"
    "Uwham.cpp"
    "Bias.cpp"
    "SimpleBias.cpp"
    "UwhamCalculationStrategy.cpp"
    "UwhamAdaptiveMethods.cpp"
    "Driver.cpp"
    "UwhamLBFGS.cpp"
    "Bin.cpp"
)

set(WHAMsrc_HEADER
    "DataFileParser.h"
    "TimeSeries.h"
    "Wham.h"
    "Uwham.h"
    "Bias.h"
    "SimpleBias.h"
    "UwhamCalculationStrategy.h"
    "UwhamAdaptiveMethods.h"
    "Driver.h"
    "UwhamLBFGS.h"
    "Bin.h"
)

target_include_directories(WHAMsrc
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

target_link_libraries(WHAMsrc PUBLIC WHAMtools WHAMparallel)

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(WHAMsrc 
    INCLUDE_SUBDIR WHAMsrc
    INCLUDES ${WHAMsrc_HEADER})