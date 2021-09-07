# Process version information
semver = meson.project_version().split('.')
version_src = configure_file(input : 'version.h.in.in',
  output : 'version.h.in',
  configuration : {
    'version_major' : semver[0],
    'version_minor' : semver[1],
    'version_patch' : semver[2],
    'version_int'   : (semver[0].to_int()*1000+semver[1].to_int())*10000+semver[2].to_int(),
    'version_vcs'   : '@VCS_TAG@'
  })

version_file = vcs_tag( input : version_src, output : 'version.h' )

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
    cpp_args : ['-DDOCTEST_CONFIG_DISABLE'],
    install : true,
    install_dir : get_option('libexecdir')
  )
endforeach

mutk_src = configure_file(input : 'mutk.cpp.in',
  output : 'mutk.cpp',
  configuration : {
    'source_dir' : meson.current_source_dir(),
    'prefix' : get_option('prefix'),
    'libexecdir' : get_option('libexecdir')
  })

executable('mutk', [mutk_src, version_file],
  include_directories : inc,
  link_with : [libmutk],
  install : true
  )