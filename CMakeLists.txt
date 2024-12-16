cmake_minimum_required(VERSION 3.14)

# In-source guard
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not supported")
endif()

# Init project
project(
    gelocus
    VERSION 0.0.0
    DESCRIPTION "Earth-centered reference frame transformations library"
    HOMEPAGE_URL "https://github.com/gunvirranu/gelocus"
    LANGUAGES C
)

# For CMake 3.21+, variable is set by default by project()
if(CMAKE_VERSION VERSION_LESS 3.21.0)
    string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}" PROJECT_IS_TOP_LEVEL)
endif()

# Create Developer Mode option, only if top-level project
if(PROJECT_IS_TOP_LEVEL)
    option(gelocus_DEVELOPER_MODE "Enable developer mode" OFF)
endif()

# ---- Warning Guard ----
# target_include_directories with SYSTEM modifier will request compiler
# to omit warnings from provided paths, if the compiler supports that.
# This is to provide a user experience similar to find_package when
# add_subdirectory or FetchContent is used to consume this project.
set(warning_guard "")
if(NOT PROJECT_IS_TOP_LEVEL)
    option(
        gelocus_INCLUDES_WITH_SYSTEM
        "Use SYSTEM modifier for gelocus's includes, disabling warnings"
        ON
    )
    mark_as_advanced(gelocus_INCLUDES_WITH_SYSTEM)
    if(gelocus_INCLUDES_WITH_SYSTEM)
        set(warning_guard SYSTEM)
    endif()
endif()

# Declare library and specify source files
add_library(
    gelocus
    src/gelocus.c
)

# Include directory
target_include_directories(
    gelocus ${warning_guard} PUBLIC
    "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
)

# Specify C99 dependency
target_compile_features(gelocus PUBLIC c_std_99)

# Developer mode
if(gelocus_DEVELOPER_MODE)
    if(NOT PROJECT_IS_TOP_LEVEL)
        message(AUTHOR_WARNING "Developer mode enabled")
    endif()

    # Hop on over to tests/CMakeLists.txt
    add_subdirectory(tests)
endif()
