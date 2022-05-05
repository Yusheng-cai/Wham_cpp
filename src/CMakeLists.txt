add_library(WHAMsrc SHARED)

target_sources(WHAMsrc PRIVATE
    "DataFileParser.cpp"
    "TimeSeries.cpp"
    "Wham.cpp"
    "Uwham.cpp"
    "BWham.cpp"
    "Bias.cpp"
    "SimpleBias.cpp"
    "BwhamCalculationStrategy.cpp"
    "BwhamLBFGS.cpp"
    "UWhamCalculationStrategy.cpp"
    "UwhamAdaptiveMethods.cpp"
    "Driver.cpp"
    "UwhamLBFGS.cpp"
    "Bin.cpp"
    "FFT.cpp"
    "UwhamReweight.cpp"
    "SparseSampling.cpp"
    "TSoperation.cpp"
    "HistogramOP.cpp"
    "CombineData.cpp"
    "Reweight.cpp"
    "UwhamConditionalReweight.cpp"
)

set(WHAMsrc_HEADER
    "DataFileParser.h"
    "TimeSeries.h"
    "Wham.h"
    "Uwham.h"
    "BWham.h"
    "Bias.h"
    "SimpleBias.h"
    "BwhamCalculationStrategy.h"
    "BwhamLBFGS.h"
    "UwhamCalculationStrategy.h"
    "UwhamAdaptiveMethods.h"
    "Driver.h"
    "UwhamLBFGS.h"
    "Bin.h"
    "FFT.h"
    "UwhamReweight.h"
    "SparseSampling.h"
    "HistogramOP.h"
    "TSoperation.h"
    "CombineData.h"
    "Reweight.h"
    "UwhamConditionalReweight.h"
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
