set(version_name "version.h")
set(version_src_file "${the_source_dir}/src/version.h.in")
set(version_dest_file "${the_binary_dir}/src/version.h")

set(the_version_major "${CMAKE_PROJECT_VERSION_MAJOR}")
set(the_version_minor "${CMAKE_PROJECT_VERSION_MINOR}")
set(the_version_patch "${CMAKE_PROJECT_VERSION_PATCH}")
set(the_version_tweak "${CMAKE_PROJECT_VERSION_TWEAK}") # not used

# Are we are in master/ a prerelease?
if("${the_version_patch}" STREQUAL "")
  set(the_version_str "${the_version_major}.${the_version_minor}-prerelease")
  set(the_version_release 0)
else()
  set(the_version_str "${the_version_major}.${the_version_minor}.${the_version_patch}")
  set(the_version_release 1)
endif()

math(EXPR the_version_num
  "(${the_version_major})*1000000 + (${the_version_minor})*1000 + (${the_version_patch}${the_version_release})"
  OUTPUT_FORMAT DECIMAL)

message(STATUS "Version: ${the_version_str}")

# Process Git Hash information
set(THE_SOURCE_BUILT_FROM_GIT)
if("$Format:$" STREQUAL "")
  set(the_version_git_hash       "$Format:%H$")
  set(the_version_git_hash_short "$Format:%h$")

  option(THE_SOURCE_IS_PATCHED "Set to ON if patches have been applied." OFF)
  if(THE_SOURCE_IS_PATCHED)
    set(the_version_git_dirty "dirty")
  endif()
elseif(GIT_FOUND)
  set(configure_code "
if (IS_DIRECTORY \"${the_source_dir}/.git\")
  set(THE_SOURCE_BUILT_FROM_GIT TRUE)
  execute_process(
    COMMAND           \"${GIT_EXECUTABLE}\"
                      rev-parse
                      HEAD
    WORKING_DIRECTORY \"${the_source_dir}\"
    RESULT_VARIABLE   git_return
    OUTPUT_VARIABLE   the_version_git_hash)
  execute_process(
    COMMAND           \"${GIT_EXECUTABLE}\"
                      rev-parse
                      --short
                      HEAD
    WORKING_DIRECTORY \"${the_source_dir}\"
    RESULT_VARIABLE   git_return
    OUTPUT_VARIABLE   the_version_git_hash_short)
  execute_process(
    COMMAND           \"${GIT_EXECUTABLE}\"
                      diff
                      --no-ext-diff
                      --quiet
                      --exit-code
    WORKING_DIRECTORY \"${the_source_dir}\"
    RESULT_VARIABLE   git_return)
  string(STRIP \"\${the_version_git_hash}\" the_version_git_hash)
  string(STRIP \"\${the_version_git_hash_short}\" the_version_git_hash_short)
  if(git_return)
    set(the_version_git_dirty \"dirty\")
  endif()

  message(STATUS \"version: \${the_version_str}\")
  message(STATUS \"git hash: \${the_version_git_hash}\")
  message(STATUS \"git short hash: \${the_version_git_hash_short}\")
  message(STATUS \"git dirty: \${the_version_git_dirty}\")
endif()
")
else()
  set(the_version_git_hash       "<unknown>")
  set(the_version_git_hash_short "<unknown>")
  set(the_version_git_dirty      "<unknown>")
endif()

set(configure_script "${version_dest_file}.cmake")

file(WRITE "${configure_script}"
  "# configure script for ${version_h_dest_file}\n"
)

set(configure_vars 
  the_version_major
  the_version_minor
  the_version_patch
  the_version_tweak
  the_version_str
  the_version_num
  the_version_git_hash
  the_version_git_hash_short
  the_version_git_dirty
)

foreach(arg IN LISTS configure_vars)
  file(APPEND "${configure_script}"
    "set(${arg} \"${${arg}}\")\n")
endforeach()

file(APPEND "${configure_script}" "${configure_code}")

file(APPEND "${configure_script}" "
set(the_version_str_long \"\${the_version_str}-g\${the_version_git_hash_short}\")
if(the_version_git_dirty)
  string(APPEND the_version_str_long \"-\${the_version_git_dirty}\")
endif()
")

file(APPEND "${configure_script}" "
configure_file(
  \"${version_src_file}\"
  \"${version_dest_file}.output\"
  @ONLY)\n")

file(APPEND "${configure_script}" "
configure_file(
  \"${version_dest_file}.output\"
  \"${version_dest_file}\"
  COPYONLY)\n")

set(clean_files
  "${version_dest_file}"
  "${version_dest_file}.output"
  )

set_directory_properties(PROPERTIES
      ADDITIONAL_MAKE_CLEAN_FILES "${clean_files}")

add_custom_command(
  OUTPUT  "${version_dest_file}.noexist"
  BYPRODUCTS "${version_dest_file}"
  COMMAND "${CMAKE_COMMAND}" -P "${configure_script}"
  MAIN_DEPENDENCY "${version_src_file}"
  DEPENDS "${version_src_file}"
          "${configure_script}"
  WORKING_DIRECTORY
          "${CMAKE_CURRENT_BINARY_DIR}"
  COMMENT "Configuring ${version_name} file \"${version_src_file}\" -> \"${version_dest_file}\"")

add_custom_target(configure-${version_name}
  DEPENDS "${version_dest_file}"
  SOURCES "${version_src_file}")
