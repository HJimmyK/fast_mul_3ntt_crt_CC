cmake_minimum_required(VERSION 3.10)
project(test LANGUAGES C)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
add_compile_options(-O3)

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-g)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -ffast-math -fomit-frame-pointer")
elseif(CMAKE_C_COMPILER_ID MATCHES "MSVC")
    add_compile_options(/Zi)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /Ot /GL")
endif()

add_executable(test main.c)