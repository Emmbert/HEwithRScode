cmake_minimum_required(VERSION 3.16)
project(HEwithRScodes LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)  # for clangd / CLion

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(HE_BUILD_EXPERIMENTS "Build the OpenFHE noise experiments and main (requires OpenFHE)" ON)
option(RS_BUILD_TESTS "Build the RS unit tests and register them with CTest" ON)

function(rs_set_warnings target)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${target} PRIVATE -Wall -Wextra)
  endif()
endfunction()

# =============================================================================
# Reed-Solomon module (independent of OpenFHE)
#   rs_core      field, polynomials, negacyclic RS code, channels, Decoder interface
#   rs_decoders  one source file per decoder + registry
#   rs_harness   generic decoder test + benchmark
# =============================================================================
add_library(rs_core
  "../../HEwithRScodes/src/rs/core/field.cpp"
        "../../HEwithRScodes/src/rs/core/rs_code.cpp"
        "../../HEwithRScodes/src/rs/core/poly_fast.cpp"
        "../../HEwithRScodes/src/rs/channel.cpp"
)
target_include_directories(rs_core PUBLIC include)
rs_set_warnings(rs_core)

add_library(rs_decoders
  "../../HEwithRScodes/src/rs/decoders/welch_berlekamp.cpp"
        "../../HEwithRScodes/src/rs/decoders/fast_gcd.cpp"
        "../../HEwithRScodes/src/rs/decoders/registry.cpp"
        # src/rs/decoders/berlekamp_massey.cpp
)
target_link_libraries(rs_decoders PUBLIC rs_core)
rs_set_warnings(rs_decoders)

add_library(rs_harness
  "../../HEwithRScodes/src/rs/harness/decoder_test.cpp"
        "../../HEwithRScodes/src/rs/harness/benchmark.cpp"
)
target_link_libraries(rs_harness PUBLIC rs_core)
rs_set_warnings(rs_harness)

add_executable(rs_bench "../../HEwithRScodes/apps/rs_bench.cpp")
target_link_libraries(rs_bench PRIVATE rs_decoders rs_harness)
rs_set_warnings(rs_bench)

# =============================================================================
# Homomorphic encryption noise experiments (OpenFHE)
# =============================================================================
if(HE_BUILD_EXPERIMENTS)
  find_package(OpenFHE REQUIRED)

  # All OpenFHE settings bundled in one interface target, so they only apply to
  # code that uses OpenFHE (and not, e.g., -Werror to the RS module).
  add_library(openfhe_deps INTERFACE)
  separate_arguments(OPENFHE_CXX_FLAGS_LIST UNIX_COMMAND "${OpenFHE_CXX_FLAGS}")
  target_compile_options(openfhe_deps INTERFACE ${OPENFHE_CXX_FLAGS_LIST})
  target_include_directories(openfhe_deps SYSTEM INTERFACE
    ${OPENMP_INCLUDES}
    ${OpenFHE_INCLUDE}
    ${OpenFHE_INCLUDE}/third-party/include
    ${OpenFHE_INCLUDE}/core
    ${OpenFHE_INCLUDE}/pke
    ${OpenFHE_INCLUDE}/binfhe
  )
  target_link_directories(openfhe_deps INTERFACE ${OpenFHE_LIBDIR} ${OPENMP_LIBRARIES})
  target_link_libraries(openfhe_deps INTERFACE ${OpenFHE_SHARED_LIBRARIES})
  find_package(OpenMP)  # OpenFHE is usually built with OpenMP; -fopenmp is also needed at link time
  if(OpenMP_CXX_FOUND)
    target_link_libraries(openfhe_deps INTERFACE OpenMP::OpenMP_CXX)
  endif()

  add_library(noise_sim
          src/noise_simulations/crypto_context.cpp
          src/noise_simulations/encryption.cpp
          src/noise_simulations/noise_experiment.cpp
          src/noise_simulations/multiplication_experiment.cpp
          src/noise_simulations/noise_growth_experiment.cpp   # new
          src/noise_simulations/param_sweep.cpp
  )
  # Both styles work: #include "noise_simulations/encryption.h" and #include "encryption.h"
  target_include_directories(noise_sim PUBLIC
          include
          include/noise_simulations
  )
  target_link_libraries(noise_sim PUBLIC openfhe_deps)

  add_executable(HEwithRScodes main.cpp)
  target_link_libraries(HEwithRScodes PRIVATE noise_sim rs_decoders rs_harness)
endif()

# =============================================================================
# Tests (CTest)
# =============================================================================
if(RS_BUILD_TESTS)
  enable_testing()
  add_subdirectory("../HEwithRScodes/tests")
endif()
