# CMakeList.txt : CMake project for TSP, include source and define
# project specific logic here.


include_directories(${PROJECT_SOURCE_DIR})

add_executable (main main.cpp GS.cpp GS.hpp)
target_link_libraries(main PUBLIC MPI::MPI_CXX)
