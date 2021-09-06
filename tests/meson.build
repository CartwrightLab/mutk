doctest_exe = executable('libmutk-doctest', ['libmutk-doctest.cpp', version_file, libmutk_sources],
  include_directories : inc,
  dependencies : [doctest_dep, eigen_dep, cli_dep, htslib_dep, minionrng_dep],
  build_by_default : false
)

fs = import('fs')

tests = fs.read('libmutk-doctest-tests.txt').strip().split('\n')

foreach test_case : tests
  test(f'[libmutk] @test_case@', doctest_exe, args :[f'--test-case=@test_case@'])
endforeach
