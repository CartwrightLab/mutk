# Process version information
semver = meson.project_version().split('.')
configure_file(input : 'version.h.in.in',
  output : 'version.h.in',
  configuration : {
    'version_major' : semver[0],
    'version_minor' : semver[1],
    'version_patch' : semver[2],
    'version_int'   : (semver[0].to_int()*1000+semver[1].to_int())*10000+semver[2].to_int(),
    'version_vcs'   : '@VCS_TAG@'
  })
version_file = vcs_tag(
    input : meson.current_build_dir() / 'version.h.in',
    output : 'version.h'
    )

inc = include_directories(['.','include'])

boost_dep = dependency('boost')
eigen_dep = dependency('eigen3')
htslib_dep = dependency('htslib')

subdir('include')
subdir('lib')

progs=['version', 'genseed', 'modelfit', 'graph']

foreach p : progs
  exe = executable('mutk-@0@'.format(p), ['mutk-@0@.cpp'.format(p), version_file],
    link_with : [libmutk],
    include_directories : inc,
    dependencies : [eigen_dep, cli_dep, htslib_dep, minionrng_dep],
    cpp_args : ['-DDOCTEST_CONFIG_DISABLE']
  )
endforeach

doctest_exe = executable('mutk-doctest', ['mutk-doctest.cpp', version_file],
  link_whole : [libmutk_doctest],
  include_directories : inc,
  dependencies : [doctest_dep, eigen_dep, cli_dep, htslib_dep, minionrng_dep]
)
