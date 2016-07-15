cmake_minimum_required(VERSION 3.2)
project(lulesh)

## Adjust to your needs

## DART Implementation
set(DASH_DART_IMPL "MPI")
## Project Sources
file(GLOB SOURCES *.cc)

set(CXX_FLAGS "-O3 -Ofast -DUSE_DASH")

## Do not edit below this line unless you know what you are doing

if("$ENV{DASH_BASE}" STREQUAL "")
  set(DASH_BASE_DIR "$ENV{HOME}/opt/dash-0.3.0/")
else()
  set(DASH_BASE_DIR $ENV{DASH_BASE})
endif()

find_package("DASH-${DASH_DART_IMPL}" REQUIRED HINTS ${DASH_BASE_DIR})
find_package(DART-BASE REQUIRED HINTS ${DASH_BASE_DIR})
find_package("DART-${DASH_DART_IMPL}" REQUIRED HINTS ${DASH_BASE_DIR})

message("Looking for DASH installation in: ${DASH_BASE_DIR}")
if(DASH-MPI_FOUND)
  message("DASH ${DASH_DART_IMPL} found")
  include_directories("${DASH_BASE_DIR}/include")
endif()

if(DART-MPI_FOUND)
  message("DART ${DASH_DART_IMPL} found")
endif()

if(DART-BASE_FOUND)
  message("DART BASE found")
endif()

set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${CXX_FLAGS}" )

## set dependencies and build options
add_executable(lulesh ${SOURCES})
target_link_libraries(lulesh dash-mpi dart-base dart-mpi)

## Print summary
message("Sources:     ${SOURCES}")
message("CXX Flags:   ${CMAKE_CXX_FLAGS}")
