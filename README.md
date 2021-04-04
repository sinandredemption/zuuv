# zuuv
Multiprecision Floating-point to Continued Fraction Cruncher

Zuuv is a program that can compute <link>simple continued fraction of a number</link> from it's floating-point representation. It is reasonably fast and can compute billions of terms on a regular PC. It's features include:
* Advanced algorithms for log-linear runtime
* Ability to utilize hard-disk to perform extremely large calculations
* Ability to automatically save and resume partially computations
* Basic multi-threading

## Installation
Windows build:
Linux build:
The program requires a minimum of four threads to run efficiently. I haven't personally tested it on a dual-core setup, but I would expect performance issues.
After downloading, to check if everything is alright, do a quick continuant benchmark (option 4) to see if everything is running fine.

## Benchmarking
There are two default benchmarks provided:
### Benchmark Continued Fraction Cruncher
Zuuv attempts to calculate the terms of regular continued fraction of a random fraction. This benchmark timing isn't very scalable beyond 4-cores (the procedure is essentially sequential in nature, and hence very difficult to parallalize beyond 4-cores. However, I am open to ideas.).
### Benchmark Continuant Cruncher
Zuuv attempts to calculate the continued of a list of random numbers (generated according to the Gauss-Kuzmin distribution to simulate continued fraction terms). The benchmark timing should benefit with increase in processor cores, as the procedure is well-parallelized.

## Doing an actual computation
A file pi_100k_hex.txt is provided to experiment.
## Compiling
Zuuv is written in C++17 and depends on the mpir-3.0.0 library.
### Visual Studio
You will need to download and extract the mpir-3.0.0 library, and enable C++17 standard.
#### mpir
1. Download and extract the mpir-3.0.0 library.
2. In Project -> Properties -> C/C++ -> General -> Additional Include Directories, add the path to mpir forlder inside mpir-3.0.0
3. In Project -> Properties -> Linker -> Input -> Additional Dependencies, add path to mpir.lib
#### C++17
Goto Project -> Properties -> C/C++ -> Language -> C++ Language Standard, and choose "ISO C++17 Standard"
### g++
1. Install the mpir library with C++ enabled (`./configure --enable-cxx`)
2. Use the following command line: `g++ *.cpp -std=c++17 -pthread -lmpir` (of course, you can further optimize by `-O2` etc.)

In case of problems, feel free to write to me.

