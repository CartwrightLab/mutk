# This CMake File defines several useful developer options

SET(DEVOPT_ENABLE_GPERFTOOLS OFF CACHE BOOL "Enable profiling with gperftools.")
SET(DEVOPT_ENABLE_COVERAGE_REPORT OFF CACHE BOOL "Enable code coverage reporting.")
SET(DEVOPT_EXCLUDE_PRETEST_FROM_ALL OFF CACHE BOOL "Do not build the pretest target by default.")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

SET(devopt_LIBRARIES)

if(DEVOPT_ENABLE_GPERFTOOLS)
  find_package(Gperftools COMPONENTS profiler)
  if(GPERFTOOLS_FOUND)
    message(STATUS "DEVOPT: Profiling with gperftools enabled. Use CPUPROFILE environmental variable to turn on profiling and specify output file.")
    set(devopt_LIBRARIES ${devopt_LIBRARIES} GPERFTOOLS::GPERFTOOLS)
  else()
    message(FATAL_ERROR "Gperftools was not found. Please disable the flag DEVOPT_ENABLE_GPERFTOOLS and try again.")
  endif()
endif()

#####################################################################
# COVERAGE REPORT

add_library(devopt_coverage INTERFACE)

if(DEVOPT_ENABLE_COVERAGE_REPORT AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options(devopt_coverage INTERFACE
    -O0        # no optimization
    -g         # generate debug info
    --coverage # sets all required flags
  )
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options(devopt_coverage INTERFACE --coverage)
  else()
    target_link_libraries(devopt_coverage INTERFACE --coverage)
  endif()

  # setup targets for generating reports
  add_subdirectory("${CMAKE_SOURCE_DIR}/coverage")
endif()

#####################################################################
# CLANG FORMAT

# Setup code formatting if possible
find_package(ClangFormat)

function(clang_format_files name)
  if(NOT ClangFormat_FOUND)
    return() # a NOOP 
  endif()
  
  foreach(source ${ARGN})
    get_filename_component(source "${source}" ABSOLUTE)
    list(APPEND sources "${source}")
  endforeach()

  add_custom_target(format_${name}
    COMMAND "${CLANG_FORMAT_EXECUTABLE}"
      -style=file -i ${sources}
    COMMENT
      "Formating sources of ${name} ..."
  )
  if(TARGET format)
    add_dependencies(format format_${name})
  else()
    add_custom_target(format DEPENDS format_${name})
  endif()

  add_custom_target(check_format_${name}
    COMMAND "${CLANG_FORMAT_EXECUTABLE}"
      -style=file -output-replacements-xml
      ${sources} > ${CMAKE_BINARY_DIR}/check_format_${name}.xml
    COMMAND ! grep -q "'\\breplacement\\b'" ${CMAKE_BINARY_DIR}/check_format_${name}.xml
    BYPRODUCTS 
      ${CMAKE_BINARY_DIR}/check_format_${name}.xml
    COMMENT
      "Checking format of sources of ${name} ..."
  )
  if(TARGET check_format)
    add_dependencies(check_format check_format_${name})
  else()
    add_custom_target(check_format DEPENDS check_format_${name})
  endif()

endfunction()

function(clang_format_target target) 
  get_target_property(target_sources "${target}" SOURCES)
  clang_format_files("${target}" ${target_sources})
endfunction()

#####################################################################
# CPPCHECK

find_package(Cppcheck)

if(Cppcheck_FOUND AND CMAKE_EXPORT_COMPILE_COMMANDS)

add_custom_target(cppcheck
  COMMAND "${CPPCHECK_EXECUTABLE}"
    "--project=${CMAKE_BINARY_DIR}/compile_commands.json"
    "--suppressions-list=${CMAKE_CURRENT_LIST_DIR}/cppcheck_suppressions.txt"
    --enable=all
    --inline-suppr
  COMMENT "Looking for programming errors with Cppcheck ..."
)

add_custom_target(check_cppcheck
  COMMAND "${CPPCHECK_EXECUTABLE}"
    "--project=${CMAKE_BINARY_DIR}/compile_commands.json"
    "--suppressions-list=${CMAKE_CURRENT_LIST_DIR}/cppcheck_suppressions.txt"
    --enable=all
    --inline-suppr
    --error-exitcode=2
    --xml
    "--output-file=${CMAKE_BINARY_DIR}/cppcheck_results.xml"
  BYPRODUCTS
    "${CMAKE_BINARY_DIR}/cppcheck_results.xml"
  COMMENT "Checking if source code passes Cppcheck ..."
)

endif()


#--suppressions-list=../../conf/cppcheck_suppressions.txt --project=compile_commands.json