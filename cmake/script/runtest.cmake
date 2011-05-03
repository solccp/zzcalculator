

if (TEST_IMSTR)
    execute_process(COMMAND ${TEST_PROG} -P -f ${SOURCEDIR}/${TEST_NAME}.in TIMEOUT 120 RESULT_VARIABLE HAD_ERROR OUTPUT_FILE output.txt ERROR_FILE error.txt) 
else ()
    execute_process(COMMAND ${TEST_PROG} ${SOURCEDIR}/${TEST_NAME}.in TIMEOUT 120 RESULT_VARIABLE HAD_ERROR OUTPUT_FILE output.txt ERROR_FILE error.txt) 
endif()

if(HAD_ERROR)
    message(FATAL_ERROR "Test failed")
endif()


execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files output.txt ${SOURCEDIR}/${TEST_NAME}.out RESULT_VARIABLE DIFFERENT1)
if(DIFFERENT1)
    message(FATAL_ERROR "Test failed - ZZ_polynomial files differ")
endif()

if (TEST_IMSTR)
    execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files intermediate_strs.yaml ${SOURCEDIR}/${TEST_NAME}.yaml RESULT_VARIABLE DIFFERENT2)
    if(DIFFERENT2)
        message(FATAL_ERROR "Test failed - intermediate structures files differ")
    endif()                 
endif()


