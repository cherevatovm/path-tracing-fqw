cmake_minimum_required(VERSION 3.15)
project(PathTracer LANGUAGES CXX)

# General settings
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(PATH_TRACER_FAST_MATH "Enable unsafe fast floating point math" OFF)
option(PATH_TRACER_NATIVE_ARCH "Optimize for native CPU (SIMD, tuning)" ON)
option(PATH_TRACER_LTO "Enable link-time optimization" ON)


# Dependency search
find_package(OpenMP REQUIRED)
find_package(OpenGL REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
find_package(GLEW REQUIRED)
find_package(Threads REQUIRED)
if(WIN32)
    find_package(imgui CONFIG REQUIRED)
else()
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(IMGUI REQUIRED imgui)
endif()


# Project's source files
file(GLOB_RECURSE PROJECT_SOURCES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
)

set(THIRDPARTY_SOURCES
    thirdparty/ImFileDialog/ImFileDialog.cpp
)

add_executable(PathTracer
    ${PROJECT_SOURCES}
    ${THIRDPARTY_SOURCES}
)


# Include directories
target_include_directories(PathTracer PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/ImFileDialog
)

if(UNIX AND NOT APPLE)
    target_include_directories(PathTracer SYSTEM PRIVATE
        ${IMGUI_INCLUDE_DIRS}
    )
    list(GET IMGUI_INCLUDE_DIRS 0 IMGUI_ROOT_DIR)
    target_include_directories(PathTracer SYSTEM PRIVATE
        ${IMGUI_ROOT_DIR}/backends
)
endif()


# Compiler flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(PathTracer PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Wno-unused-parameter
        -Wno-sign-compare

        $<$<CONFIG:Debug>:-O0 -g3 -fno-omit-frame-pointer>
        $<$<CONFIG:Release>:-O3>
    )
endif()

if(MSVC)
    target_compile_options(PathTracer PRIVATE
        /W4
        $<$<CONFIG:Release>:/O2>
        $<$<CONFIG:Debug>:/Od /Zi>
    )

    target_compile_definitions(PathTracer PRIVATE
        NOMINMAX
        WIN32_LEAN_AND_MEAN
        _CRT_SECURE_NO_WARNINGS
    )
endif()


# Additional optimization options
target_compile_options(PathTracer PRIVATE
    # GCC/Clang
    $<$<AND:$<BOOL:${PATH_TRACER_FAST_MATH}>,$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang>>:
        -ffast-math
        -fno-math-errno
        -ffp-contract=fast
    >

    # MSVC
    $<$<AND:$<BOOL:${PATH_TRACER_FAST_MATH}>,$<CONFIG:Release>,$<CXX_COMPILER_ID:MSVC>>:
        /fp:fast
        /fp:contract
    >
)

target_compile_options(PathTracer PRIVATE
    # GCC/Clang
    $<$<AND:$<BOOL:${PATH_TRACER_NATIVE_ARCH}>,$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang>>:
        -march=native
    >

    # MSVC
    $<$<AND:$<BOOL:${PATH_TRACER_NATIVE_ARCH}>,$<CONFIG:Release>,$<CXX_COMPILER_ID:MSVC>>:
        /arch:AVX2
    >
)

include(CheckIPOSupported)
check_ipo_supported(RESULT IPO_SUPPORTED OUTPUT IPO_ERROR)

if(PATH_TRACER_LTO AND IPO_SUPPORTED)
    set_property(TARGET PathTracer PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()


# Linking
target_link_libraries(PathTracer PRIVATE
    OpenMP::OpenMP_CXX
    OpenGL::GL
    glfw
    GLEW::GLEW
    Threads::Threads
)

if(WIN32)
    target_link_libraries(PathTracer PRIVATE imgui::imgui)
else()
    target_link_libraries(PathTracer PRIVATE ${IMGUI_LIBRARIES})
endif()

if(UNIX AND NOT APPLE)
    target_link_options(PathTracer PRIVATE
        $<$<CONFIG:Debug>:-rdynamic>
    )
endif()


message(STATUS "Project: PathTracer")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Fast-math: ${PATH_TRACER_FAST_MATH}")
message(STATUS "Native arch: ${PATH_TRACER_NATIVE_ARCH}")
message(STATUS "Link-time optimization: ${PATH_TRACER_LTO}")
message(STATUS "OpenMP found: ${OpenMP_CXX_FOUND}")
