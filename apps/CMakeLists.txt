add_executable(main main.c)
set_target_properties(main PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
target_link_libraries(main PRIVATE argtable3 m ljs)
target_include_directories(main PUBLIC "${INCLUDE_DIR}")
target_compile_features(main PRIVATE c_std_99 )



if (USE_ASAN)
  target_compile_options(main PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_libraries(main PRIVATE asan ubsan)
  target_link_options(main PRIVATE -fsanitize=address,undefined)
endif()
