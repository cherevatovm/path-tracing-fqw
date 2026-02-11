cmake_minimum_required(VERSION 3.15)
project(PathTracer LANGUAGES CXX)
include(FetchContent)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(PATH_TRACER_FAST_MATH "Enable fast math optimizations" OFF)
option(PATH_TRACER_NATIVE_ARCH "Optimize for native CPU" ON)
option(PATH_TRACER_SANITIZERS "Enable address/UB sanitizers (Debug only)" OFF)

# Search for system libraries
find_package(OpenMP REQUIRED)
find_package(OpenGL REQUIRED)
find_package(glfw3 3.3 REQUIRED)
find_package(GLEW REQUIRED)
find_package(Threads REQUIRED)

# Search for ImGui via pkg-config
find_package(PkgConfig REQUIRED)
pkg_check_modules(IMGUI QUIET imgui)

if(IMGUI_FOUND)
    message(STATUS "ImGui found: ${IMGUI_INCLUDE_DIRS}")
else()
    message(WARNING "ImGui not found via pkg-config")
endif()

# Project's source files
file(GLOB_RECURSE PROJECT_SOURCES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
)

# Thirdparty sources
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
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty
    ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/ImFileDialog
)

if(IMGUI_FOUND)
    target_include_directories(PathTracer SYSTEM PRIVATE
        ${IMGUI_INCLUDE_DIRS}
    )
endif()

# Set warning parameters and build configuration
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

if(PATH_TRACER_FAST_MATH)
    target_compile_options(PathTracer PRIVATE
        $<$<CONFIG:Release>:-ffast-math>
    )
endif()

if(PATH_TRACER_NATIVE_ARCH)
    target_compile_options(PathTracer PRIVATE
        $<$<CONFIG:Release>:-march=native>
    )
endif()

# Sanitizers
if(PATH_TRACER_SANITIZERS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(PathTracer PRIVATE
        -fsanitize=address,undefined
        -fno-omit-frame-pointer
        -O1
    )
    target_link_options(PathTracer PRIVATE
        -fsanitize=address,undefined
    )
endif()

# Linking
target_link_libraries(PathTracer PRIVATE
    OpenMP::OpenMP_CXX
    OpenGL::GL
    glfw
    GLEW::GLEW
    Threads::Threads
    ${IMGUI_LIBRARIES}
)

target_link_options(PathTracer PRIVATE
    $<$<CONFIG:Debug>:-rdynamic>
)

message(STATUS "Project: PathTracer")
message(STATUS "Compiler: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "C++ standard: ${CMAKE_CXX_STANDARD}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}") 
message(STATUS "Fast-math: ${PATH_TRACER_FAST_MATH}")
message(STATUS "Native arch: ${PATH_TRACER_NATIVE_ARCH}")
message(STATUS "Sanitizers: ${PATH_TRACER_SANITIZERS}")
message(STATUS "OpenMP found: ${OpenMP_CXX_FOUND}")
