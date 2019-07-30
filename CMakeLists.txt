# Copyright (c) 2019 Reed A. Cartwright <reed@cartwright.ht>
# 
# This file is part of the Ultimate Source Code Project.
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

cmake_minimum_required(VERSION 3.8.0 FATAL_ERROR)
cmake_policy(SET CMP0025 NEW)
cmake_policy(SET CMP0048 NEW)

# guard against in-source builds
set(the_source_dir "${CMAKE_CURRENT_SOURCE_DIR}")
set(the_binary_dir "${CMAKE_CURRENT_BINARY_DIR}")

if("${the_source_dir}" STREQUAL "${the_binary_dir}")
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

################################################################################
#  Mutation Toolkit Project Configuration
# 

## Important Variables
set(CMAKE_MODULE_PATH "${the_source_dir}/modules")
set(CMAKE_BUILD_TYPE_INIT "RelWithDebInfo")
set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)

## Compiling options
option(USE_STATIC_LIBS "Link against static libraries.")
option(USE_STATIC_FLAG "Use -static flag to create completely static binaries.")
option(DEVEL_MODE "Enable developer mode.")

################################################################################
# Init Project

project(MUTK VERSION 0.1
  DESCRIPTION "Mutation Toolkit"
  HOMEPAGE_URL "https://github.com/ultimatesource/mutk/"
  LANGUAGES CXX C
)

## Enable Developer CMake Flags
include("conf/devopt.cmake")

# Gitmeta data and semantic versioning
include("conf/gitmeta.cmake")
include("conf/semver.cmake")

add_semver_file(${CMAKE_PROJECT_NAME} METADATA_GIT)

## Guard against bad build-type strings
if(CMAKE_BUILD_TYPE)
  string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
  if(    NOT cmake_build_type_tolower STREQUAL "debug"
     AND NOT cmake_build_type_tolower STREQUAL "release"
     AND NOT cmake_build_type_tolower STREQUAL "relwithdebinfo"
     AND NOT cmake_build_type_tolower STREQUAL "none")
    message(FATAL_ERROR "Unsupported build type \"${CMAKE_BUILD_TYPE}\". Allowed values are Debug, Release, RelWithDebInfo, and None (case-insensitive).")
  endif()
  message(STATUS "Compiling with build type '${CMAKE_BUILD_TYPE}'.")
  if(cmake_build_type_tolower STREQUAL "debug")
    if(NOT DEFINED CMAKE_EXPORT_COMPILE_COMMANDS)
      set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
    endif()
  endif()
endif()

## Turn on unit and executable testing
enable_testing()
include("conf/doctest.cmake")

set(pretest_all ALL)
if(DEVOPT_EXCLUDE_PRETEST_FROM_ALL)
  set(pretest_all)
endif()

add_custom_target(pretest ${pretest_all})
add_dependencies(pretest configure)

configure_file(mutk.sublime-project mutk.sublime-project COPYONLY)

################################################################################
# Installation Information

STRING(TOLOWER ${CMAKE_PROJECT_NAME} basename_for_install)

SET(CMAKE_DATA_DIR "share/${basename_for_install}" CACHE STRING
  "Install location for data (relative to prefix).")
SET(CMAKE_DOC_DIR "doc/${basename_for_install}" CACHE STRING
  "Install location for documentation (relative to prefix).")
SET(CMAKE_MAN_DIR "man" CACHE STRING
  "Install location for man pages (relative to prefix).")
MARK_AS_ADVANCED(CMAKE_DATA_DIR CMAKE_DOC_DIR CMAKE_MAN_DIR)

INCLUDE(GNUInstallDirs)
SET(CMAKE_INSTALL_LIBEXECDIR  "${CMAKE_INSTALL_LIBEXECDIR}/${basename_for_install}")

INSTALL(FILES README.md LICENSE DESTINATION ${CMAKE_DOC_DIR})

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(tests)
ADD_SUBDIRECTORY(doc)

########### Add uninstall target ###############
CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/conf/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake") 
