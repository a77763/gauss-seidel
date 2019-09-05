# CMakeList.txt : Top-level CMake project file, do global configuration
# and include sub-projects here.
#
cmake_minimum_required (VERSION 3.10)

set(CMAKE_C_COMPILER_WORKS 1)
set(CMAKE_CXX_COMPILER_WORKS 1)

project ("GS" LANGUAGES CXX)

find_package(MPI REQUIRED)
#message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_STANDARD 11)

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    SET(CMAKE_INSTALL_PREFIX
        "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "Installation prefix, default '${CMAKE_CURRENT_SOURCE_DIR}'" FORCE
    )
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

set(THIS_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin")

# Include sub-projects.
add_subdirectory ("source")

install(TARGETS main
    RUNTIME DESTINATION bin
)
