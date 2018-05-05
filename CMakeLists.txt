cmake_minimum_required(VERSION 3.9)

project(smallpt)

set(ISPC "ispc" CACHE FILEPATH "ispc executable")

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
	set(ISPC_FLAGS -g -O0)
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo")
  set(ISPC_FLAGS -g)
else()
  set(ISPC_FLAGS )
endif()

add_custom_command(
  OUTPUT smallpt_ispc.obj
  COMMAND ${ISPC} --arch=x86-64 -DISPC_USE_CONCRT ${CMAKE_SOURCE_DIR}/smallpt.ispc ${ISPC_FLAGS} --opt=fast-math -o smallpt_ispc.obj -h smallpt_ispc.h
  DEPENDS smallpt.ispc)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(smallpt smallpt_ispc.obj smallpt.cpp smallpt_c.cpp tasksys.cpp)

target_link_libraries(smallpt)
