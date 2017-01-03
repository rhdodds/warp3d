# Try to find HYPRE

# Defines the following:
#     HYPRE_FOUND
#     HYPRE_INCLUDE_DIRS
#     HYPRE_LIBRARIES

FIND_PATH(HYPRE_INCLUDE_DIR HYPRE.h /usr/include /usr/local/include /opt/packages/hypre/include)
FIND_LIBRARY(HYPRE_LIBRARY NAMES HYPRE PATHS /usr/lib /usr/local/lib /opt/packages/hypre/lib)

if (HYPRE_INCLUDE_DIR AND HYPRE_LIBRARY)
      SET(HYPRE_FOUND TRUE)
endif()

if (HYPRE_FOUND)
      set(HYPRE_LIBRARIES ${HYPRE_LIBRARY})
      set(HYPRE_INCLUDE_DIRS ${HYPRE_INCLUDE_DIR})
else()
      if (HYPRE_FIND_REQUIRED)
            message(FATAL_ERROR "Could not find HYPRE")
      endif()
endif()
