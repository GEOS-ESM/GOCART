function(copy_directory source destination)
  if(EXISTS ${source})
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${source} ${destination})
  else()
    message(FATAL_ERROR "Source directory not found: ${source}")
  endif()
endfunction()

function(link_directory source destination)
  if(EXISTS ${source})
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${source} ${destination})
  else()
    message(FATAL_ERROR "Source directory not found: ${source}")
  endif()
endfunction()

function(copy_restarts root_dir expdir)
  file(GLOB checkpoint_subdirs LIST_DIRECTORIES true ${root_dir}/checkpoints/*)
  foreach(ckpt_dir IN LISTS checkpoint_subdirs)
    get_filename_component(ckpt_name ${ckpt_dir} NAME)
    if(IS_DIRECTORY ${ckpt_dir} AND NOT ckpt_name STREQUAL "last")
      copy_directory(${ckpt_dir} ${expdir}/checkpoints/${ckpt_name})
    endif()
  endforeach()
endfunction()

function(run_geos num_procs case_name expdir)
  execute_process(
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${num_procs} -prepend-rank ${MPIEXEC_PREFLAGS} ${MY_BINARY_DIR}/GEOS.x cap.yaml
    RESULT_VARIABLE CMD_RESULT
    WORKING_DIRECTORY ${expdir}
    COMMAND_ECHO STDOUT
  )
  if(EXISTS ${expdir}/PET0.ESMF_LogFile)
    execute_process(COMMAND ${CMAKE_COMMAND} -E cat ${expdir}/PET0.ESMF_LogFile)
  endif()
  if(CMD_RESULT)
    message(FATAL_ERROR "Error running ${case_name}")
  endif()
endfunction()

function(compare_results baseline_dir current_dir)
  file(GLOB baseline_files ${baseline_dir}/*.nc)
  foreach(baseline_file IN LISTS baseline_files)
    get_filename_component(fname ${baseline_file} NAME)
    message(STATUS "Comparing ${fname}")
    execute_process(
      COMMAND cmp ${baseline_file} ${current_dir}/${fname}
      RESULT_VARIABLE CMP_RESULT
      OUTPUT_VARIABLE CMP_OUTPUT
      ERROR_VARIABLE CMP_OUTPUT
    )
    if(CMP_RESULT)
      message(FATAL_ERROR "Files differ: ${fname}\n${CMP_OUTPUT}")
    endif()
  endforeach()
endfunction()

function(run_case case_name regression_data_dir)
  string(RANDOM LENGTH 24 expdir)
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E make_directory ${expdir}
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_LIST_DIR}/${case_name} ${expdir}
  )

  set(root_dir ${regression_data_dir}/${case_name})
  set(num_procs "6")
  set(checkpoints_dir ${root_dir}/checkpoints/last)

  if(NOT EXISTS ${root_dir})
    message(STATUS "Regression data not found for ${case_name}: ${root_dir} -- skipping")
    return()
  endif()

  copy_restarts(${root_dir} ${expdir})
  link_directory(${regression_data_dir}/ExtData ${expdir}/ExtData)
  run_geos(${num_procs} ${case_name} ${expdir})

  if(NOT EXISTS ${checkpoints_dir})
    message(WARNING "Baseline directory [${checkpoints_dir}] not found. Skipping comparison.")
    return()
  endif()
  compare_results(${checkpoints_dir} ${expdir}/checkpoints/last)

  # execute_process(
  #   COMMAND ${CMAKE_COMMAND} -E rm -rf ${expdir}
  # )
endfunction()

run_case(${TEST_CASE} ${REGRESSION_DATA_DIR})
