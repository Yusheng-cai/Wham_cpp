add_library(WHAMtools SHARED)

target_sources(WHAMtools PRIVATE
    "Algorithm.cpp"
    "RandomNumber.cpp"
    "InputParser.cpp"
    "CommandLineArguments.cpp"
    "FileSystem.cpp"
    "OutputFunction.cpp"
)

set(WHAMTOOLS_HEADER
    "Assert.h"
    "Algorithm.h"
    "CommonTypes.h"
    "GenericFactory.h"
    "CommonOperations.h"
    "InputParser.h"
    "RandomNumber.h"
    "CommandLineArguments.h"
    "Registry.h"
    "FileSystem.h"
    "OutputFunction.h"
)

target_include_directories(WHAMtools
                        PUBLIC 
                        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
                        $<INSTALL_INTERFACE:${PROJECT_BINARY_DIR}/include>
)
# target_compile_options(tools PRIVATE -fPIC)


include(InstallLibraryWithStandardSetup)
InstallLibraryWithStandardSetup(WHAMtools 
    INCLUDE_SUBDIR WHAMtools
    INCLUDES ${WHAMTOOLS_HEADER})

