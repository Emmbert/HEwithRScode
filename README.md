# Running Experiments on HE combined with Reed-Solomon codes

This instructions are working on Apple silicon M3.

## 1. Prerequisites
If no c++ compiler installe:
```bash
xcode-select --install      # Xcode Command Line Tools (clang, make, git)
brew install cmake          # if not already installed
```

## 2. Build and install OpenFHE
#### 2.1. Install OpenMP support

Apple clang doesn't include OpenMP by default, so install it via Homebrew:

```bash
brew install libomp
```

Confirm it installed correctly:
```bash
brew --prefix libomp
ls $(brew --prefix libomp)/lib
```
You should see `libomp.dylib` listed.

### 2.2 Install OpenFHE

```bash
git clone https://github.com/openfheorg/openfhe-development.git
cd openfhe-development
mkdir build && cd build

cmake .. \
  -DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include" \
  -DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include" \
  -DOpenMP_C_LIB_NAMES="omp" \
  -DOpenMP_CXX_LIB_NAMES="omp" \
  -DOpenMP_omp_LIBRARY="$(brew --prefix libomp)/lib/libomp.dylib"

make -j$(sysctl -n hw.ncpu)
sudo make install
```

**After this finishes**, confirm where things landed:
```bash
find / -name "openfhe.h" 2>/dev/null
find / -name "libOPENFHEpke*" 2>/dev/null
```
Typically:
- Headers → `/usr/local/include/openfhe/...`
- Libraries → `/usr/local/lib/libOPENFHEpke.dylib` (and `core`, `binfhe`)


## 3. Compile your program

```bash
clang++ -std=c++17 main.cpp -o tests \
  -I/usr/local/include/openfhe \
  -I/usr/local/include/openfhe/pke \
  -I/usr/local/include/openfhe/core \
  -I/usr/local/include/openfhe/binfhe \
  -L/usr/local/lib \
  -lOPENFHEpke -lOPENFHEcore -lOPENFHEbinfhe \
  -Xpreprocessor -fopenmp \
  -I$(brew --prefix libomp)/include \
  -L$(brew --prefix libomp)/lib -lomp \
  -Wl,-rpath,/usr/local/lib
```

## 4. Run it

```bash
./test
```

---

## X. (Optional) Using CMake instead of manual compile commands

Once you're compiling more than one file regularly, this avoids retyping long commands:

**CMakeLists.txt:**
```cmake
cmake_minimum_required(VERSION 3.5.1)
project(MyOpenFHEApp)

find_package(OpenFHE REQUIRED)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${OpenFHE_INCLUDE_DIR})
link_directories(${OpenFHE_LIBDIR})

add_executable(myapp main.cpp)
target_link_libraries(myapp PUBLIC ${OpenFHE_SHARED_LIBRARIES})
```

Then:
```bash
mkdir build && cd build
cmake ..
make
./myapp
```

# Negacyclic RS decoder test bench

Code: R_p = Z_p[x]/(x^n + 1), n = 2^d, p prime with 2n | p - 1.
Locators x_i = psi^(2i+1) (the roots of x^n + 1, i.e. the primitive 2n-th roots of unity),
codeword = negacyclic NTT of the message polynomial m(x), deg m < k.

```
include/rs/      field.hpp, poly.hpp     GF(p) + polynomial arithmetic
                 rs_code.hpp             RSCode (negacyclic NTT encode/interpolate, parameters)
                 channel.hpp             ReceivedWord, symbol-error and additive Gaussian channels
                 decoder.hpp             abstract Decoder interface
                 decoder_test.hpp        generic correctness test for any decoder
                 benchmark.hpp           runAllDecoders(), printComparison(), writeCsv()
src/             implementations of the above
decoders/        one .hpp/.cpp per decoder (+ selfTest<Name>()), registry.cpp lists all
apps/main.cpp    rs_bench: self tests + comparison
```

Build:  `cmake -S . -B build && cmake --build build`
(or `g++ -std=c++17 -O2 -Iinclude -Idecoders src/*.cpp decoders/*.cpp apps/main.cpp -o rs_bench`)

Run:    `./build/rs_bench 4096 3072 50 40961 --errors 512 --errors 513 --sigma 0.2 --csv decoder_benchmarks_out.csv`
(arguments: n k trials p; p = 0 picks the smallest suitable prime; `--psi` sets your own root)

Add a decoder: create `decoders/<name>.hpp/.cpp` deriving from `rs::Decoder`,
write `selfTest<Name>()` using `rs::testDecoder()`, add one line each in `registry.cpp`
(and the .cpp to `CMakeLists.txt`).


### Run the scripts to find parameter sets that are close to the multiplication depth broder

```python
source ../HEmeetsECC/venv/bin/activate
cd simulation_results
python3 find_candidates.py calibrate --p 65537 --K 10 --n 16384 32768 --out calib_params.csv          
```

```bash
cmake --build build -j$(sysctl -n hw.ncpu)
./build/HEwithRScodes simulation_results/calib_params.csv simulation_results/calib 10
```

```python
python3 find_candidates.py select --p 65537 --K 10 --calib calib --targets 0.1 0.25 0.5 --out candidates.csv
```

```bash
./build/HEwithRScodes simulation_results/candidates.csv simulation_results/final 11
```

```python
python3 find_candidates.py select --p 65537 --K 10 --calib calib --targets 0.1 0.25 0.5 --out candidates.csv
```