include(CheckIncludeFile)
include(CheckFunctionExists)
include(CheckLibraryExists)

add_subdirectory(lib)

include_directories(BEFORE "${CMAKE_CURRENT_BINARY_DIR}")

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/mutk.sh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/mutk.sh"
  @ONLY)
install(PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/mutk.sh" DESTINATION ${CMAKE_INSTALL_BINDIR} RENAME mutk)
