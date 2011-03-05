
execute_process(COMMAND ${TEST_PROG} ${SOURCEDIR}/${TEST_NAME}.in TIMEOUT 120 RESULT_VARIABLE HAD_ERROR OUTPUT_FILE output.txt ERROR_FILE error.txt) 

if(HAD_ERROR)
    message(FATAL_ERROR "Test failed")
endif()

execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files output.txt ${SOURCEDIR}/${TEST_NAME}.out RESULT_VARIABLE DIFFERENT)
if(DIFFERENT)
    message(FATAL_ERROR "Test failed - files differ")
endif()                 
