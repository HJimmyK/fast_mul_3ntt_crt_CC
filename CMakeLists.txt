cmake_minimum_required(VERSION 3.10)
project(test LANGUAGES C)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
add_compile_options(-O3)

# 核心新增：添加调试信息（让VTune识别函数名），分编译器配置
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    # GCC/Clang：-g 生成调试符号，与-O3共存不影响性能
    add_compile_options(-g)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -ffast-math -fomit-frame-pointer")
elseif(CMAKE_C_COMPILER_ID MATCHES "MSVC")
    # MSVC：/Zi 生成调试符号（.pdb文件），/Od是关闭优化（此处保留你的/Ot优化，仅加/Zi）
    add_compile_options(/Zi)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /Ot /GL")
endif()

add_executable(test main.c)