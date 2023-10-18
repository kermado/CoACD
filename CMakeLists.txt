cmake_minimum_required(VERSION 3.24)
project(CoACD LANGUAGES C CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W0 -fpermissive -D_USE_MATH_DEFINES")
    set(CMAKE_SHARED_LIBRARY_PREFIX "lib")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fpermissive -D_USE_MATH_DEFINES")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g3")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -msimd128 -mavx")
set(CMAKE_CXX_FLAGS "-pthread -fno-exceptions")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wno-int-in-bool-context -Wsign-compare")
    set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fsanitize=address")
endif()

add_compile_definitions(DISABLE_SPDLOG)
include(eigen)

file(GLOB_RECURSE COACD_SRC "src/*.cc" "src/*.cpp" "public/*.cpp")
add_library(coacd STATIC ${COACD_SRC})
target_include_directories(coacd PUBLIC public)
target_include_directories(coacd PRIVATE 3rd/triangle)

if (NOT DEFINED WASM)
    add_executable(main main.cpp)
    target_link_libraries(main coacd)
endif()