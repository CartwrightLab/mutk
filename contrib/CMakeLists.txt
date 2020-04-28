add_library(doctest INTERFACE)
add_library(doctest::doctest ALIAS doctest)
target_include_directories(doctest INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/doctest")

add_library(minionrng STATIC minionrng/minion.cc)
add_library(minionrng::minionrng ALIAS minionrng)
target_include_directories(minionrng PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/minionrng")

add_library(CLI11 INTERFACE)
add_library(CLI11::CLI11 ALIAS CLI11)
target_include_directories(CLI11 INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/CLI")
