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
    "FFT.cpp"
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
    "FFT.h"
)

target_include_directories(WHAMsrc
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)

target_include_directories(WHAMsrc PUBLIC "${FFTW3_INCLUDE_DIRS}")
message("In src, ${FFTW3_LIBRARIES}")
target_link_libraries(WHAMsrc PUBLIC WHAMtools WHAMparallel ${FFTW3_LIBRARIES})

include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(WHAMsrc 
    INCLUDE_SUBDIR WHAMsrc
    INCLUDES ${WHAMsrc_HEADER})