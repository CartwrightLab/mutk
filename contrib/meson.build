doctest_dep = declare_dependency(include_directories: include_directories('doctest'))
cli_dep = declare_dependency(include_directories : include_directories('CLI'))

minionrng = static_library(
    'minionrng', ['minionrng/minion.cc'],
    include_directories : include_directories('minionrng'),
)

minionrng_dep = declare_dependency(
    include_directories : include_directories('minionrng'),
    link_with : minionrng
)
