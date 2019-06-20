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

CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0 FATAL_ERROR)
CMAKE_POLICY(SET CMP0025 NEW)
CMAKE_POLICY(SET CMP0048 NEW)

################################################################################
#  Mutation Toolkit Project Configuration
# 

## Important Variables
SET(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/modules")
SET(CMAKE_BUILD_TYPE_INIT "RelWithDebInfo")
SET_PROPERTY(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS ON)

## Compiling options
SET(USE_STATIC_LIBS OFF CACHE BOOL "Link against static libraries.")
SET(USE_STATIC_FLAG OFF CACHE BOOL "Use -static flag to create completely static binaries.")
SET(DEVEL_MODE OFF CACHE BOOL "Enable developer mode.")

################################################################################
# Init Project

PROJECT(MUTK VERSION 0.1.0
  DESCRIPTION "Mutation Toolkit"
  HOMEPAGE_URL "https://github.com/ultimatesource/mutk/"
  LANGUAGES CXX C)
SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_EXTENSIONS ON)

# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

## Guard against bad build-type strings
IF(CMAKE_BUILD_TYPE)
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
ENDIF()


################################################################################
# Always make a full statement RPATH
# See http://www.itk.org/Wiki/CMake_RPATH_handling

SET(CMAKE_SKIP_BUILD_RPATH FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}" is_system_dir)
IF("${is_system_dir}" STREQUAL "-1")
  SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}")
ENDIF("${is_system_dir}" STREQUAL "-1")
################################################################################

IF(DEVEL_MODE)
  ADD_DEFINITIONS(-DMUTK_DEVELOPER_MODE)
ENDIF()

## Turn on unit and executable testing
ENABLE_TESTING()

## Enable Developer CMake Flags
INCLUDE(devopt)

################################################################################
# Version Information

## Detect Project Version Information
FIND_PACKAGE(Git)
SET(COMMIT_ID "unknown_commit_id")
IF(GIT_FOUND)
  EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} describe --always --dirty
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
    RESULT_VARIABLE res_var OUTPUT_VARIABLE git_commit_id
    ERROR_QUIET)
  IF(${res_var} EQUAL 0 )
    STRING( REPLACE "\n" "" COMMIT_ID "${git_commit_id}" )
  ENDIF()
ENDIF()

MESSAGE(STATUS "${CMAKE_PROJECT_NAME} version: ${CMAKE_PROJECT_VERSION}-${COMMIT_ID}.")

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
  "${CMAKE_CURRENT_SOURCE_DIR}/modules/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)

ADD_CUSTOM_TARGET(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake") 
