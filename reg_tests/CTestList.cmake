#=============================================================================
# Functions for adding tests / Categories of tests
#=============================================================================

# Standard regression test
function(add_test_r testname np)
    add_test(${testname} sh -c "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/${nalu_ex_name} ${MPIEXEC_POSTFLAGS} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${NALU_GOLD_NORMS_DIR}/test_files/${testname}/${testname}.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 10800 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "regression")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_r)

# Standard performance test
function(add_test_p testname np)
    add_test(${testname} sh -c "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/${nalu_ex_name} ${MPIEXEC_POSTFLAGS} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${NALU_GOLD_NORMS_DIR}/test_files/${testname}/${testname}.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 10800 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "performance")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_p)

# Regression test with single restart
function(add_test_r_rst testname np)
    add_test(${testname} sh -c "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/${nalu_ex_name} ${MPIEXEC_POSTFLAGS} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${NALU_GOLD_NORMS_DIR}/test_files/${testname}/${testname}.norm.gold ${TOLERANCE}; ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/${nalu_ex_name} ${MPIEXEC_POSTFLAGS} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_rst.i -o ${testname}_rst.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname}_rst ${NALU_GOLD_NORMS_DIR}/test_files/${testname}/${testname}_rst.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 10800 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "regression")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_r_rst)

# Regression test with input
function(add_test_r_inp testname np)
    add_test(${testname} sh -c "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/${nalu_ex_name} ${MPIEXEC_POSTFLAGS} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname} ${NALU_GOLD_NORMS_DIR}/test_files/${testname}/${testname}.norm.gold ${TOLERANCE}; ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/${nalu_ex_name} ${MPIEXEC_POSTFLAGS} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_Input.i -o ${testname}_Input.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname}_Input ${NALU_GOLD_NORMS_DIR}/test_files/${testname}/${testname}_Input.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 10800 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "regression")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_r_inp)

# Verification test with three resolutions
function(add_test_v3 testname np)
    add_test(${testname} sh -c "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/${nalu_ex_name} ${MPIEXEC_POSTFLAGS} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R0.i -o ${testname}_R0.log && ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/${nalu_ex_name} ${MPIEXEC_POSTFLAGS} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R1.i -o ${testname}_R1.log && ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/${nalu_ex_name} ${MPIEXEC_POSTFLAGS} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R2.i -o ${testname}_R2.log && python ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/norms.py")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 10800 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "verification")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_v3)

# Verification test with two resolutions
function(add_test_v2 testname np)
    add_test(${testname} sh -c "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/${nalu_ex_name} ${MPIEXEC_POSTFLAGS} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R0.i -o ${testname}_R0.log && ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/${nalu_ex_name} ${MPIEXEC_POSTFLAGS} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}_R1.i -o ${testname}_R1.log && python ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/norms.py")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 10800 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "verification")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_v2)

# Regression test that runs with different numbers of processes
function(add_test_r_np testname np)
    add_test(${testname}Np${np} sh -c "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/${nalu_ex_name} ${MPIEXEC_POSTFLAGS} -i ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i -o ${testname}Np${np}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail.sh ${testname}Np${np} ${NALU_GOLD_NORMS_DIR}/test_files/${testname}/${testname}Np${np}.norm.gold ${TOLERANCE}")
    set_tests_properties(${testname}Np${np} PROPERTIES TIMEOUT 10800 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "regression")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_r_np)

# Standard unit test
function(add_test_u testname np)
    if(${np} EQUAL 1)
      set(GTEST_SHUFFLE "--gtest_shuffle")
    else()
      unset(GTEST_SHUFFLE)
    endif()
    add_test(${testname} sh -c "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/${utest_ex_name} ${GTEST_SHUFFLE}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 2000 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "unit")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_u)

# GPU unit test
function(add_test_u_gpu testname np)
    set(FILTER "--gtest_filter=BasicKokkos.discover_execution_space:*.NGP*")
    add_test(${testname} sh -c "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/${utest_ex_name} ${FILTER}")
    set_tests_properties(${testname} PROPERTIES TIMEOUT 2000 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "unit")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
endfunction(add_test_u_gpu)

# Regression test with catalyst capability
function(add_test_r_cat testname np ncat)
    if(ENABLE_PARAVIEW_CATALYST)
      if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i.in)
        add_test(${testname} sh -c "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${np} ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/${nalu_ex_catalyst_name} -i ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}/${testname}_catalyst.i -o ${testname}.log && ${CMAKE_CURRENT_SOURCE_DIR}/pass_fail_catalyst.sh ${testname} ${ncat}")
        set_tests_properties(${testname} PROPERTIES TIMEOUT 10800 PROCESSORS ${np} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}" LABELS "regression")
        file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
        set(CATALYST_FILE_INPUT_DECK_COMMAND "catalyst_file_name: catalyst.txt")
        configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/${testname}.i.in
                       ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname}/${testname}_catalyst.i @ONLY)
        file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_files/${testname}/catalyst.txt
             DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test_files/${testname})
      endif()
    else()
      add_test_r(${testname} ${np})
    endif()
endfunction(add_test_r_cat)

if(NOT ENABLE_CUDA)


  #=============================================================================
  # Unit tests
  #=============================================================================
  add_test_u(unitTest1 1)
  add_test_u(unitTest2 2)


else(NOT ENABLE_CUDA)

  #=============================================================================
  # GPU unit tests
  #=============================================================================
  add_test_u_gpu(unitTestGPU 1)

endif(NOT ENABLE_CUDA)
