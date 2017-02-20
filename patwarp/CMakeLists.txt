# All our sources
# CMake actually handles modules correctly!!!
list(APPEND pw_main
      "patwarp_module.f"
      "patwarp.f"
      )

add_executable(${PATWARP_NAME} ${pw_main})
target_link_libraries(${PATWARP_NAME} ${METIS_LIBRARIES})
