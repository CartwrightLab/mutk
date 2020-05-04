# Copyright (c) 2020 Reed A. Cartwright <reed@cartwright.ht>
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

include(CheckIncludeFile)
include(CheckFunctionExists)
include(CheckLibraryExists)

add_subdirectory(lib)

# Setup a target that will build version.h.in -> version.h
set(configure_script_code "
include(\"${CMAKE_CURRENT_BINARY_DIR}/../${CMAKE_PROJECT_NAME}.version\")\n
")
add_configured_file(version.h
    INPUT "${CMAKE_CURRENT_SOURCE_DIR}/version.h.in"
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/version.h"
    EXTRA_DEPS "${CMAKE_CURRENT_BINARY_DIR}/../${CMAKE_PROJECT_NAME}.version"
    VARIABLES ${CMAKE_PROJECT_NAME}_SEMVER_VERSION
              ${CMAKE_PROJECT_NAME}_SEMVER_INTEGER
              ${CMAKE_PROJECT_NAME}_SEMVER_MAJOR
              ${CMAKE_PROJECT_NAME}_SEMVER_MINOR
              ${CMAKE_PROJECT_NAME}_SEMVER_PATCH
)
add_dependencies(configure-version.h configure-${CMAKE_PROJECT_NAME}.version)
unset(configure_script_code)

#####################################################################
# mutk binary

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/mutk.cpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/mutk.cpp"
  @ONLY)
add_executable(mutk "${CMAKE_CURRENT_BINARY_DIR}/mutk.cpp")
target_link_libraries(mutk libmutk)
install(TARGETS mutk RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

#####################################################################
# mutk subcommands

set(mutk_subcommands version genseed modelfit
  #graph
  )

foreach(cmd ${mutk_subcommands})
  add_executable(mutk-${cmd} mutk-${cmd}.cpp)
  target_link_libraries(mutk-${cmd} libmutk CLI11::CLI11)

  install(TARGETS mutk-${cmd} RUNTIME DESTINATION ${CMAKE_INSTALL_LIBEXECDIR})

  clang_format_target(mutk-${cmd})
  clang_tidy_target(mutk-${cmd})  
endforeach()

# target-specific options
target_link_libraries(mutk-genseed minionrng::minionrng)

###########################################################
# Unit Test support
add_executable(mutk-doctest EXCLUDE_FROM_ALL mutk-doctest.cpp)

target_link_libraries(mutk-doctest libmutk-doctest doctest::doctest)

doctest_discover_tests(mutk-doctest PROPERTIES "FIXTURES_REQUIRED;doctest_fixture")

add_test(NAME [=["[pretest] build_mutk-doctest"]=] COMMAND "${CMAKE_COMMAND}"
  --build "${CMAKE_BINARY_DIR}"
  --config "$<CONFIG>"
  --target mutk-doctest
)
set_tests_properties([=["[pretest] build_mutk-doctest"]=] PROPERTIES FIXTURES_SETUP doctest_fixture)

add_dependencies(pretest mutk-doctest)
