cmake_minimum_required(VERSION 3.0)
project(lulesh)

## Adjust to your needs

## DART Implementation
set(DASH_DART_IMPL "MPI")
## Project Sources
file(GLOB SOURCES *.cc)

#set(CXX_FLAGS "-O3 -Ofast -DUSE_DASH")
#set(CXX_FLAGS "-O0 -g -DUSE_DASH")
set(CXX_FLAGS "-O3 -g -DUSE_DASH")
#set(CXX_FLAGS "-O3 -g -DUSE_DASH -p")
#set(CXX_FLAGS "${CXX_FLAGS} -DDASH_USE_GET")
set(CMAKE_CXX_STANDARD 14)

## Do not edit below this line unless you know what you are doing

if("$ENV{DASH_BASE}" STREQUAL "")
  set(DASH_BASE_DIR "$ENV{HOME}/opt/dash-0.3.0")
  message(STATUS "DASH_BASE not set, using ${DASH_BASE_DIR}")
else()
  set(DASH_BASE_DIR $ENV{DASH_BASE})
  message(STATUS "DASH_BASE set to ${DASH_BASE_DIR}")
endif()

set(DASH_BASE_CMAKE_DIR "${DASH_BASE_DIR}/share/cmake")

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


foreach (MPI_INC_ ${MPI_INCLUDE_PATH})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I${MPI_INC_}" )
endforeach()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DASH_CXX_FLAGS} ${CXX_FLAGS}" )
set(CXX_FLAGS ${CMAKE_CXX_FLAGS})

string(TOLOWER ${DASH_DART_IMPL} DASH_DART_IMPL_LC)

## set dependencies and build options
add_executable(
  ${PROJECT_NAME}
  ${SOURCES})

#set(DASH_LIBRARIES "${DASH_LIBRARIES} -ldart-tasking -L$ENV{HOME}/opt/extrae-3.5.4-gnu7.3.0-shared/lib -lptmpitrace")
#set(DASH_LIBRARIES "${DASH_LIBRARIES} -ldart-tasking")
#set(DASH_LIBRARIES "${DASH_LIBRARIES}")

#set(EXTRAE_LIBRARIES "-L$ENV{HOME}/opt-cray/extrae-3.5.4-gnu7.3.0-shared/lib -lptmpitrace")

if("$ENV{TCMALLOC_ROOT}" STREQUAL "")
  set(TCMALLOC_LIBRARIES "")
else()
  set(TCMALLOC_LIBRARIES "-L$ENV{TCMALLOC_ROOT}/lib/ -ltcmalloc_minimal")
  message(STATUS "Using TCMalloc libraries: ${TCMALLOC_LIBRARIES}")
endif()

#set (DMAPP_LIBRARIES "lulesh-tasks/CMakeLists.txt")

target_link_libraries(
  ${PROJECT_NAME}
  ${EXTRAE_LIBRARIES}
  ${TCMALLOC_LIBRARIES}
  ${DASH_LIBRARIES}
  "${DASH_BASE_DIR}/lib/libdart-tasking.a"
#  "$ENV{HOME}/opt-cray/extrae-3.5.4-gnu7.3.0/lib/libptmpitrace.a"
  ${DMAPP_LIBRARIES}
  ${MPI_C_LIBRARIES})

if(MPI_LINK_FLAGS)
  set_target_properties(
    ${PROJECT_NAME} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

set_target_properties(
  ${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "${CMAKE_CXX_FLAGS}")

#set(CXX_FLAGS "${CXX_FLAGS} -DDASH_USE_GET")
#set(CMAKE_CXX_FLAGS "${CXX_FLAGS}")
add_executable(
  ${PROJECT_NAME}_get
  ${SOURCES})
target_link_libraries(
  ${PROJECT_NAME}_get
  ${TCMALLOC_LIBRARIES}
  ${EXTRAE_LIBRARIES}
  ${DASH_LIBRARIES}
  "${DASH_BASE_DIR}/lib/libdart-tasking.a"
#  "$ENV{HOME}/opt-cray/extrae-3.5.4-gnu7.3.0/lib/libptmpitrace.so"
  ${DMAPP_LIBRARIES}
  ${MPI_C_LIBRARIES})
set_target_properties(
  ${PROJECT_NAME}_get PROPERTIES COMPILE_FLAGS "-DDASH_USE_GET")
if(MPI_LINK_FLAGS)
  set_target_properties(
    ${PROJECT_NAME}_get PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

add_executable(
  ${PROJECT_NAME}_copyin
  ${SOURCES})
target_link_libraries(
  ${PROJECT_NAME}_copyin
  ${TCMALLOC_LIBRARIES}
  ${EXTRAE_LIBRARIES}
  ${DASH_LIBRARIES}
  "${DASH_BASE_DIR}/lib/libdart-tasking.a"
#  "$ENV{HOME}/opt-cray/extrae-3.5.4-gnu7.3.0/lib/libptmpitrace.so"
  ${DMAPP_LIBRARIES}
  ${MPI_C_LIBRARIES})
set_target_properties(
  ${PROJECT_NAME}_copyin PROPERTIES COMPILE_FLAGS "-DDASH_USE_COPYIN")
if(MPI_LINK_FLAGS)
  set_target_properties(
    ${PROJECT_NAME}_copyin PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()




## Print summary
message(STATUS "Sources:          ${SOURCES}")
message(STATUS "DASH Libraries:   ${DASH_LIBRARIES}")
message(STATUS "CXX Flags:        ${CMAKE_CXX_FLAGS}")
