# Try to find the current version of MKL

# Defines the following:
#     MKL_FOUND
#     MKL_INCLUDE_DIRS
#     MKL_LIBRARIES

if(NOT DEFINED ENV{MKLROOT})
      set(MKL_FOUND FALSE)
      if(MKL_FIND_REQUIRED)
            message(FATAL_ERROR "MKLROOT not set, cannot find MKL")
      endif()
      return()
endif()

set(MKL_INCLUDE_DIRS "$ENV{MKLROOT}/include")
set(MKL_LIBRARIES "-L$ENV{MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl")

set(MKL_FOUND TRUE)
