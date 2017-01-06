cmake_minimum_required(VERSION 3.2)

# Custom finders
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

# Project name and languages
project(warp3d Fortran)

# Version
add_definitions(-DBuild=3250)

# Some general structure
if (WIN32)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/run_windows_64)
      add_definitions(-DWindows)
elseif(APPLE)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/run_linux_mac_os_x)
      add_definitions(-DOSX)
else()
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/run_linux_em64t)
      add_definitions(-DLinux)
endif()

# Configuration options: MPI and Hypre
option(USE_MPI "Compile the hybrid MPI/OpenMP version" OFF)
option(USE_HYPRE "Compile with the Hypre AMG solver" OFF)
option(USE_NEML "Include the NEML material library" OFF)

if(USE_MPI)
      find_package(MPI REQUIRED)
      include_directories(${MPI_Fortran_INCLUDE_PATH})
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_COMPILE_FLAGS}")
      link_directories(${MPI_Fortran_LIBRARIES})
      
      add_definitions(-DUSE_MPI)
      set(EXE_NAME "warp3d.mpi_omp")
      
      # This will require some thought for non-Intel compilers
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mt_mpi")
else()
      set(EXE_NAME "warp3d.omp")
endif()

if(USE_HYPRE)
      find_package(HYPRE REQUIRED)
      include_directories(${HYPRE_INCLUDE_DIRS})
      add_definitions(-DUSE_HYPRE)
endif()

if(USE_NEML)
      set(NEML_LIBRARY_PATH "/home/messner/projects/neml/lib" CACHE STRING "")
      set(NEML_LIBRARIES "-L${NEML_LIBRARY_PATH} -lneml")
      add_definitions(-DUSE_NEML)
endif()

# Force the preprocessor for all fortran files
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")

# Everyone needs openmp
FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif()

# Ditto MKL
FIND_PACKAGE(MKL REQUIRED)
if (MKL_FOUND)
      include_directories(${MKL_INCLUDE_DIRS})
endif()

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
      # General flags
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpconstant -fixed -72 -traceback")

      # Optimized build
      set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -no-inline-factor -ftz -ip -xHOST")
      # Debug build
      set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -check all,noarg_temp_created -fpe0")
else()
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -m64")
      # Please note the manual override to the Intel threading library
      # otherwise you don't get Pardiso threads
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed -liomp5")
      message(WARNING "WARP3D not tested on compilers other than Intel!")
      # Optimized build
      set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
      # Debug build
      set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -ffpe-trap=invalid")
endif()


add_subdirectory(src)
