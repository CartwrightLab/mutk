project('mutk', 'cpp',
  version : '0.1.9000',
  license : 'MIT',
  default_options : ['cpp_std=c++17']
)

subdir('contrib')
subdir('src')
subdir('tests')

configure_file(input : 'mutk.sublime-project.in',
  output : 'mutk.sublime-project',
  configuration : {
    'source_dir' : meson.project_source_root(),
    'project_name' : meson.project_name()
  }
)
