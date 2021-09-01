
inc = include_directories('include')

add_global_arguments('-D', language : 'cpp')

boost_dep = dependency('boost')
eigen_dep = dependency('eigen3')
doctest_dep = declare_dependency(include_directories: include_directories('../contrib/doctest'))

subdir('include')
subdir('lib')

