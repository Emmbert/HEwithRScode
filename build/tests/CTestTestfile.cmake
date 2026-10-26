# CMake generated Testfile for 
# Source directory: /Users/emma/Desktop/HEandECC/HEwithRScodes/tests
# Build directory: /Users/emma/Desktop/HEandECC/HEwithRScodes/build/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_field "/Users/emma/Desktop/HEandECC/HEwithRScodes/build/tests/test_field")
set_tests_properties(test_field PROPERTIES  LABELS "rs;fast" TIMEOUT "120" _BACKTRACE_TRIPLES "/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;11;add_test;/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;15;rs_add_test;/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;0;")
add_test(test_rs_code "/Users/emma/Desktop/HEandECC/HEwithRScodes/build/tests/test_rs_code")
set_tests_properties(test_rs_code PROPERTIES  LABELS "rs;fast" TIMEOUT "120" _BACKTRACE_TRIPLES "/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;11;add_test;/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;16;rs_add_test;/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;0;")
add_test(test_welch_berlekamp "/Users/emma/Desktop/HEandECC/HEwithRScodes/build/tests/test_welch_berlekamp")
set_tests_properties(test_welch_berlekamp PROPERTIES  LABELS "rs;fast" TIMEOUT "120" _BACKTRACE_TRIPLES "/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;11;add_test;/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;17;rs_add_test;/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;0;")
add_test(test_poly_fast "/Users/emma/Desktop/HEandECC/HEwithRScodes/build/tests/test_poly_fast")
set_tests_properties(test_poly_fast PROPERTIES  LABELS "rs;fast" TIMEOUT "120" _BACKTRACE_TRIPLES "/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;11;add_test;/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;18;rs_add_test;/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;0;")
add_test(test_fast_gcd "/Users/emma/Desktop/HEandECC/HEwithRScodes/build/tests/test_fast_gcd")
set_tests_properties(test_fast_gcd PROPERTIES  LABELS "rs;fast" TIMEOUT "120" _BACKTRACE_TRIPLES "/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;11;add_test;/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;19;rs_add_test;/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;0;")
add_test(test_welch_berlekamp_n4096 "/Users/emma/Desktop/HEandECC/HEwithRScodes/build/tests/test_welch_berlekamp" "--large")
set_tests_properties(test_welch_berlekamp_n4096 PROPERTIES  LABELS "rs;slow" TIMEOUT "600" _BACKTRACE_TRIPLES "/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;22;add_test;/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;0;")
add_test(test_fast_gcd_n4096 "/Users/emma/Desktop/HEandECC/HEwithRScodes/build/tests/test_fast_gcd" "--large")
set_tests_properties(test_fast_gcd_n4096 PROPERTIES  LABELS "rs;slow" TIMEOUT "600" _BACKTRACE_TRIPLES "/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;24;add_test;/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;0;")
add_test(rs_bench_smoke "/Users/emma/Desktop/HEandECC/HEwithRScodes/build/rs_bench" "256" "128" "3" "0" "--errors" "10" "--errors" "65" "--sigma" "0.2")
set_tests_properties(rs_bench_smoke PROPERTIES  LABELS "rs;fast" TIMEOUT "120" _BACKTRACE_TRIPLES "/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;28;add_test;/Users/emma/Desktop/HEandECC/HEwithRScodes/tests/CMakeLists.txt;0;")
