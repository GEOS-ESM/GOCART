include("${ESMA_REGRESSION_HELPERS}")

function(copy_restarts root_dir expdir)
  file(GLOB checkpoint_subdirs LIST_DIRECTORIES true ${root_dir}/checkpoints/*)
  foreach(ckpt_dir IN LISTS checkpoint_subdirs)
    get_filename_component(ckpt_name ${ckpt_dir} NAME)
    if(IS_DIRECTORY ${ckpt_dir} AND NOT ckpt_name STREQUAL "last")
      copy_directory(${ckpt_dir} ${expdir}/checkpoints/${ckpt_name})
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
