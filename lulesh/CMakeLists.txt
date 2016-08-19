cmake_minimum_required(VERSION 3.0)
project(lulesh)

## Adjust to your needs

## DART Implementation
set(DASH_DART_IMPL "MPI")
## Project Sources
file(GLOB SOURCES *.cc)

set(CXX_FLAGS "-O3 -Ofast -DUSE_DASH")
set(CMAKE_CXX_STANDARD 11)

## Do not edit below this line unless you know what you are doing

if("$ENV{DASH_BASE}" STREQUAL "")
  set(DASH_BASE_DIR "$ENV{HOME}/opt/dash-0.3.0")
  message(STATUS "DASH_BASE not set, using ${DASH_BASE_DIR}")
else()
  set(DASH_BASE_DIR $ENV{DASH_BASE})
  message(STATUS "DASH_BASE set to ${DASH_BASE_DIR}")
endif()

set(DASH_BASE_CMAKE_DIR "${DASH_BASE_DIR}/cmake")

find_package("DASH-${DASH_DART_IMPL}" REQUIRED HINTS ${DASH_BASE_CMAKE_DIR})
find_package(MPI)

message(STATUS "Looking for DASH installation in: ${DASH_BASE_DIR}")
if(DASH-MPI_FOUND)
  message(STATUS "DASH ${DASH_DART_IMPL} found")
endif()

if(DART-MPI_FOUND)
  message(STATUS "DART ${DASH_DART_IMPL} found")
endif()

message(STATUS "MPI include path: ${MPI_INCLUDE_PATH}")
message(STATUS "MPI libraries:")
foreach (MPI_C_LIB ${MPI_C_LIBRARIES})
  message(STATUS "    " ${MPI_C_LIB})
endforeach()
message(STATUS "MPI link flags:")
foreach (MPI_LINK_FLAG ${MPI_LINK_FLAGS})
  message(STATUS "   " ${MPI_LINK_FLAG})
endforeach()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${MPI_INCLUDE_PATH}" )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DASH_CXX_FLAGS} ${CXX_FLAGS}" )
string(TOLOWER ${DASH_DART_IMPL} DASH_DART_IMPL_LC)

## set dependencies and build options
add_executable(
  ${PROJECT_NAME}
  ${SOURCES})

target_link_libraries(
  ${PROJECT_NAME}
  ${DASH_LIBRARIES}
  ${MPI_C_LIBRARIES})

if(MPI_LINK_FLAGS)
  set_target_properties(
    ${PROJECT_NAME} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

## Print summary
message(STATUS "Sources:          ${SOURCES}")
message(STATUS "DASH Libraries:   ${DASH_LIBRARIES}")
message(STATUS "CXX Flags:        ${CMAKE_CXX_FLAGS}")
