# set(IMAGE_INCLUDE_DIR ${INCLUDE_DIR}/image)

add_library( ljs STATIC
        # store.c ${CURRENT_INCLUDE_DIR}/store.h
        system.c ${INCLUDE_DIR}/system.h
        simulation.c ${INCLUDE_DIR}/simulation.h
        # velocity_verlet.c ${INCLUDE_DIR}/velocity_verlet.h
        io.c ${INCLUDE_DIR}/io.h

        # ${INCLUDE_DIR}/type.h
        # ${INCLUDE_DIR}/global.h
        )

set_target_properties( ljs  PROPERTIES 
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

target_include_directories( ljs PUBLIC "${INCLUDE_DIR}" )

target_link_libraries( ljs PRIVATE m )


if (USE_ASAN)
  target_compile_options(ljs PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
  target_link_libraries(ljs PRIVATE asan ubsan)
  target_link_options(ljs PRIVATE -fsanitize=address,undefined)
endif()

# add_subdirectory(convolution_layer)
# add_subdirectory(neural_network)
