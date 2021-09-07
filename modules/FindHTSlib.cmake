include(LibFindMacros)

find_package(ZLIB)

libfind_pkg_detect(HTSLIB htslib FIND_PATH htslib/hts.h FIND_LIBRARY hts)

libfind_version_header(HTSLIB htslib/hts.h HTS_VERSION QUIET)

if(NOT HTSLIB_VERSION AND NOT CMAKE_CROSSCOMPILING)
  find_program(HTSLIB_EXECUTABLE NAMES htsfile)
  execute_process(
    COMMAND ${HTSLIB_EXECUTABLE} "--version"
    OUTPUT_VARIABLE hts_version_raw
  )
  if(hts_version_raw MATCHES "htsfile \\(htslib\\) ([0-9.]+)")
    set(HTSLIB_VERSION ${CMAKE_MATCH_1})
  endif()
endif()

libfind_process(HTSLIB)

if(HTSLIB_FOUND AND NOT TARGET HTSLIB::HTSLIB)
  add_library(HTSLIB::HTSLIB UNKNOWN IMPORTED)
  set_target_properties(HTSLIB::HTSLIB PROPERTIES
    IMPORTED_LOCATION ${HTSLIB_LIBRARIES}
    INTERFACE_INCLUDE_DIRECTORIES ${HTSLIB_INCLUDE_DIRS}
    INTERFACE_LINK_LIBRARIES ZLIB::ZLIB
  )
endif()
