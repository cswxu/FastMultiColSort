list(APPEND experiments
    main_tworound_exhaustive
    main_run_an_instance
    main_tpch_exhaustive_actual
    main_tpch_complete_scenario
    main_tpch_costmodel_accuracy
	main_run_an_instance_hashing
    main_tworound_exhaustive_hashing
    collect_rank
    collect_mcs_time
    collect_whole_query_time
    collect_transpose
    )

foreach(ee ${experiments})
    add_executable(${ee} "${ee}.cpp")
    target_link_libraries(${ee} pcm_wrapper multiAttrSort-core )
endforeach()

