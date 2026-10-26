# Test executables for the RS module. Each returns 0 on success.
#   ctest                    -> fast tests
#   ctest -L slow            -> only the n = 4096 tests
#   ctest -LE slow           -> everything except slow tests

function(rs_add_test name)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} PRIVATE rs_decoders rs_harness)
  target_include_directories(${name} PRIVATE "../../../HEwithRScodes/tests")
  rs_set_warnings(${name})
  add_test(NAME ${name} COMMAND ${name})
  set_tests_properties(${name} PROPERTIES LABELS "rs;fast" TIMEOUT 120)
endfunction()

rs_add_test(test_field)
rs_add_test(test_rs_code)
rs_add_test(test_welch_berlekamp)
rs_add_test(test_poly_fast)
rs_add_test(test_fast_gcd)

# Same executable, target length n = 4096
add_test(NAME test_welch_berlekamp_n4096 COMMAND test_welch_berlekamp --large)
set_tests_properties(test_welch_berlekamp_n4096 PROPERTIES LABELS "rs;slow" TIMEOUT 600)
add_test(NAME test_fast_gcd_n4096 COMMAND test_fast_gcd --large)
set_tests_properties(test_fast_gcd_n4096 PROPERTIES LABELS "rs;slow" TIMEOUT 600)

# Smoke test: the benchmark app runs end to end on a small code
add_test(NAME rs_bench_smoke COMMAND rs_bench 256 128 3 0 --errors 10 --errors 65 --sigma 0.2)
set_tests_properties(rs_bench_smoke PROPERTIES LABELS "rs;fast" TIMEOUT 120)
