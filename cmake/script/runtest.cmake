

if (TEST_IMSTR)
    execute_process(COMMAND ${TEST_PROG} -P -f ${SOURCEDIR}/${TEST_NAME}.in TIMEOUT 120 RESULT_VARIABLE HAD_ERROR OUTPUT_FILE output.txt ERROR_FILE error.txt) 
elseif(TEST_BL)
    execute_process(COMMAND ${TEST_PROG} -b ${SOURCEDIR}/${TEST_NAME}.bl ${SOURCEDIR}/${TEST_NAME}.bin TIMEOUT 120 RESULT_VARIABLE HAD_ERROR OUTPUT_FILE output.txt ERROR_FILE error.txt) 
elseif(TEST_XML)
    execute_process(COMMAND ${TEST_PROG} -Q ${SOURCEDIR}/${TEST_NAME}.in TIMEOUT 120 RESULT_VARIABLE HAD_ERROR OUTPUT_FILE output.txt ERROR_FILE error.txt) 
elseif(TEST_CML)
    execute_process(COMMAND ${TEST_PROG} -p ${LEVEL} ${SOURCEDIR}/${TEST_NAME}.in TIMEOUT 120 RESULT_VARIABLE HAD_ERROR OUTPUT_FILE output.txt ERROR_FILE error.txt) 
else ()
    execute_process(COMMAND ${TEST_PROG} ${SOURCEDIR}/${TEST_NAME}.in TIMEOUT 120 RESULT_VARIABLE HAD_ERROR OUTPUT_FILE output.txt ERROR_FILE error.txt) 
endif()

if(HAD_ERROR)
    message(FATAL_ERROR "Test failed")
endif()

if (TEST_CML)
    execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files substructures.cml ${SOURCEDIR}/${TEST_NAME}.cml RESULT_VARIABLE DIFFERENT1)
elseif (TEST_XML)
    execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files output.txt ${SOURCEDIR}/${TEST_NAME}.xml RESULT_VARIABLE DIFFERENT1)
else ()
    execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files output.txt ${SOURCEDIR}/${TEST_NAME}.out RESULT_VARIABLE DIFFERENT1)
endif()
if(DIFFERENT1)
        message(FATAL_ERROR "Test failed - ZZ_polynomial files differ")
endif()

if (TEST_IMSTR)
    execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files intermediate_strs.yaml ${SOURCEDIR}/${TEST_NAME}.yaml RESULT_VARIABLE DIFFERENT2)
    if(DIFFERENT2)
        message(FATAL_ERROR "Test failed - intermediate structures files differ")
    endif()                 
endif()




