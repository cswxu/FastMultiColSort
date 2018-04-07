include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

list(APPEND test_list
		sortnetwork_test
        mergekernel_test
        composer_test
		misc_test
		project_profiling_test
		enum_plan_test
		rrs_test
		mergesort_profiling_test
#		subgroup_test
    )

# find_program(MEMCHECK_CMD valgrind )

foreach(tt ${test_list})
    add_executable(${tt}  "${tt}.cpp")
    target_link_libraries(${tt} multiAttrSort-core gtest gtest_main)
#    add_test(NAME ${tt} COMMAND ${MEMCHECK_CMD} --tool=memcheck --leak-check=yes ./${tt} --gtest_color=yes)
    add_test(NAME ${tt} COMMAND ${tt} --gtest_color=yes)
    add_dependencies(check-build ${tt})
endforeach()



