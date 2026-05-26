include(TestSetup)

# Set the libraries needed for test
set(_core_test_lib src tools parallel)

# compile the list of names of the tests
set(TESTS_NAMES "")

# find the bash program, if not then crash the process 
find_program(BASH_PROGRAM bash)
if (NOT BASH_PROGRAM)
    message(FATAL_ERROR "Bash program is not found on this computer.")
endif()

# define the test script
set(TEST_RUN_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh")

# Add unit tests for lattice testing

set(PROGRAM_PATH "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Wham")
set(TESTDATA "${CMAKE_CURRENT_LIST_DIR}/testdata/ModelPotential1d")

# Now add the test which tests the main program (only test if they run) 
set(OMP_THREADS 1 4 8)

# Current official CTest suite. Other fixture directories under test/ are
# intentionally not registered until their expected behavior is reviewed.
set(MAIN_TESTS "")
list(APPEND MAIN_TESTS "testAdaptive")
list(APPEND MAIN_TESTS "testLBFGS")

# Add the test and run with the bash script run_test.sh
foreach(THREADS IN LISTS OMP_THREADS)
    foreach(TEST_NAME IN LISTS MAIN_TESTS)
        set(TEST_SCRIPT "${TEST_NAME}/input.dat")
        set(TEST_DIR "${CMAKE_CURRENT_LIST_DIR}/${TEST_NAME}")
        set(TEST_CMD ${BASH_PROGRAM} ${TEST_RUN_SCRIPT} -p ${PROGRAM_PATH} -d ${TEST_DIR} 
        -o ${CMAKE_CURRENT_BINARY_DIR}/pji.out -o ${CMAKE_CURRENT_BINARY_DIR}/norm.out -o ${CMAKE_CURRENT_BINARY_DIR}/kl.out -r pjiref.out -r normref.out -r klref.out 
             -i ${CMAKE_CURRENT_LIST_DIR}/${TEST_SCRIPT} -n ${THREADS} -a ${TESTDATA})

        # set working directory to be the directory processed currently
        add_test(NAME "${TEST_NAME}::OMP_${THREADS}"
                COMMAND ${TEST_CMD}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
        ) 
    endforeach()
endforeach()

add_executable(testWhamTools testWhamTools/test_wham_tools.cpp)
target_include_directories(testWhamTools PUBLIC ${CMAKE_SOURCE_DIR})
target_link_libraries(testWhamTools PUBLIC WHAMsrc WHAMtools WHAMparallel)
add_test(NAME testWhamTools COMMAND testWhamTools)

add_custom_target(build_test COMMAND ${CMAKE_CTEST_COMMAND} 
                             DEPENDS ${TESTS_NAMES})
