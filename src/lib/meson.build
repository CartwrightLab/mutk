
libmutk_sources = [
  'version.cpp',
  'call.cpp',
  'pedigree.cpp',
  'utility.cpp',
  'newick.cpp',
  'relationship_graph.cpp',
  'mutation.cpp'
]

libmutk = static_library('mutk', libmutk_sources,
	include_directories : inc,
	dependencies : [boost_dep, doctest_dep, eigen_dep]
)
