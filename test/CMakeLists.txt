file(GLOB TEST_FILES *.in)
foreach(loop_var IN LISTS TEST_FILES)
    get_filename_component(FNAME ${loop_var} NAME_WE)
    add_test(NAME ${FNAME}
        COMMAND ${CMAKE_COMMAND}
        -DTEST_PROG=$<TARGET_FILE:ZZ_polynomial>
        -DTEST_NAME=${FNAME}
        -DSOURCEDIR=${CMAKE_HOME_DIRECTORY}/test
        -P ${CMAKE_HOME_DIRECTORY}/cmake/script/runtest.cmake)
endforeach()
foreach(loop_var IN LISTS TEST_FILES)
    get_filename_component(FNAME ${loop_var} NAME_WE)
    add_test(NAME ${FNAME}_imstr
        COMMAND ${CMAKE_COMMAND}
        -DTEST_PROG=$<TARGET_FILE:ZZ_polynomial>
        -DTEST_NAME=${FNAME}
        -DTEST_IMSTR=ON
        -DSOURCEDIR=${CMAKE_HOME_DIRECTORY}/test
        -P ${CMAKE_HOME_DIRECTORY}/cmake/script/runtest.cmake)
endforeach()
file(GLOB TEST_FILES *.bin)
foreach(loop_var IN LISTS TEST_FILES)
    get_filename_component(FNAME ${loop_var} NAME_WE)
    add_test(NAME ${FNAME}
        COMMAND ${CMAKE_COMMAND}
        -DTEST_PROG=$<TARGET_FILE:ZZ_polynomial>
        -DTEST_NAME=${FNAME}
        -DTEST_BL=ON
        -DSOURCEDIR=${CMAKE_HOME_DIRECTORY}/test
        -P ${CMAKE_HOME_DIRECTORY}/cmake/script/runtest.cmake)
endforeach()

