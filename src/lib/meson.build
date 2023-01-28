libmutk_sources = files([
  'version.cpp',
  'call.cpp',
  'pedigree.cpp',
  'utility.cpp',
  'newick.cpp',
  'mutation.cpp',
  'graph_builder.cpp',
  'graph_peeler.cpp',
  'junction_tree.cpp',
  'potential.cpp',
  'potential-cloning.cpp',
  'potential-selfing.cpp',
  'mutation_builder.cpp'
])

libmutk_deps = [boost_dep, doctest_dep, eigen_dep, xtensor_dep, xblas_dep]

libmutk = static_library('mutk', [libmutk_sources, version_file],
  include_directories : inc,
  dependencies : libmutk_deps,
  cpp_args : ['-DDOCTEST_CONFIG_DISABLE']
)
