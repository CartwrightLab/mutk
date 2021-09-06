doctest_exe = executable('libmutk-doctest', ['libmutk-doctest.cpp', version_file, libmutk_sources],
  include_directories : inc,
  dependencies : [doctest_dep, eigen_dep, cli_dep, htslib_dep, minionrng_dep],
  build_by_default : false
)

test_list = files('libmutk-doctest-tests.txt')

doctest_list_wrapper = find_program('doctest-list-wrapper.sh')

run_target('update-tests',
  command : [doctest_list_wrapper, 'update', doctest_exe, test_list],
)

test('[pretest] Checking test list completeness', doctest_list_wrapper,
  args : ['test', doctest_exe, test_list],
  is_parallel : false,
  priority : 1000)

fs = import('fs')
tests = fs.read(test_list).strip().split('\n')
foreach test_case : tests
  test(f'[libmutk] @test_case@', doctest_exe, args : [f'--test-case=@test_case@'])
endforeach
