cmake_minimum_required(VERSION 3.2)
project(pattern-visualizer)

## Adjust to your needs

set(INSTALL_PREFIX "$ENV{HOME}/public_html/${PROJECT_NAME}"
    CACHE STRING "Installation path")

## Install path
set(CMAKE_INSTALL_PREFIX ${INSTALL_PREFIX})

## DART Implementation
set(DASH_DART_IMPL "MPI")
## Project Sources
set(SOURCES "src/main")

set(CXX_FLAGS "-O3 -Ofast")

## Do not edit below this line unless you know what you are doing

if("$ENV{DASH_BASE}" STREQUAL "")
  set(DASH_BASE_DIR "$ENV{HOME}/opt/dash-0.3.0")
else()
  set(DASH_BASE_DIR $ENV{DASH_BASE})
endif()

set(DASH_BASE_CMAKE_DIR "${DASH_BASE_DIR}/cmake")

find_package("DASH-${DASH_DART_IMPL}" REQUIRED HINTS ${DASH_BASE_CMAKE_DIR})

message("Looking for DASH installation in: ${DASH_BASE_DIR}")
if(DASH-MPI_FOUND)
  message("DASH ${DASH_DART_IMPL} found")
endif()

if(DART-MPI_FOUND)
  message("DART ${DASH_DART_IMPL} found")
endif()

set( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${DASH_CXX_FLAGS} ${CXX_FLAGS}" )
string(TOLOWER ${DASH_DART_IMPL} DASH_DART_IMPL_LC)

## set dependencies and build options
add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} ${DASH_LIBRARIES})

## Install
file(COPY htdocs DESTINATION .)
add_custom_command(
  TARGET ${PROJECT_NAME}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${PROJECT_NAME}> ${CMAKE_BINARY_DIR}/htdocs/
)

install(DIRECTORY htdocs/
        DESTINATION .)
install(TARGETS ${PROJECT_NAME}
        DESTINATION .)

## Print summary
message("Sources:         ${SOURCES}")
message("DASH Libraries:  ${DASH_LIBRARIES}")
message("CXX Flags:       ${CMAKE_CXX_FLAGS}")
