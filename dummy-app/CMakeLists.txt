cmake_minimum_required(VERSION 3.2)
project(dummy_app)

set(DASH_BASE_DIR "$ENV{HOME}/opt/dash-0.3.0/")

find_package(DASH-MPI REQUIRED HINTS ${DASH_BASE_DIR})
find_package(DART-BASE REQUIRED HINTS ${DASH_BASE_DIR})
find_package(DART-MPI REQUIRED HINTS ${DASH_BASE_DIR})

if(DASH-MPI_FOUND)
  message("DASH MPI found")
  include_directories("${DASH_BASE_DIR}/include")
else (DASH-MPI_FOUND)
  message("DASH MPI not found")
endif (DASH-MPI_FOUND)

if(DART-MPI_FOUND)
  message("DART MPI found")
endif(DART-MPI_FOUND)

if(DART-BASE_FOUND)
  message("DART BASE found")
endif(DART-BASE_FOUND)

add_executable(dummy_app main.cpp)
target_link_libraries(dummy_app dash-mpi dart-base dart-mpi)
