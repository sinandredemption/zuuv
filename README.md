# Zuuv
Multiprecision Floating-point to Continued Fraction Cruncher

Zuuv is a program that can compute [regular continued fraction](https://en.wikipedia.org/wiki/Continued_fraction) of a number from it's floating-point representation (a file containing hex or dec digits). It's features include:
* Advanced algorithms for log-linear runtime
* Ability to utilize hard-disk to perform extremely large calculations
* Ability to automatically save and resume partially computations
* Basic multi-threading
With enough dedication, Zuuv can be utilized to break world records for computation of regular continued fraction terms of important mathematical constants such as Pi, Euler-Mascheroni etc, even on regular PCs.

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
### Benchmark Disk-based Multiplication
Zuuv attempts to benchmark a critical multiplication routine required in disk-based compuatations. This procedure again is well parallelized and should benefit well from scaling.

## Doing an actual computation
Zuuv takes in a file containing the stream of digits of a number and calculates it's regular continued fraction terms in multiple iterations, each stored in the subdirectory `iterations`. A file pi_1m_hex.txt, containing 1 million hexadecimal digits of pi, is provided to experiment with Zuuv before launching a serious computation.

There are two compuatation modes.

| Computation Mode | RAM requirements | Disk requirements | Resumable |  Speed  | Comments |
| ---------------- | ---------------- | ----------------- | --------- | ------- | -------- |
|     RAM based    |       High       |    Output only    |     No    | Fastest | Use for small compuations that can fit in RAM |
|    Disk based    | Adjustible (low) |       High        |    Yes    |   Fast  | Use for large compuations that can't entirely be done in RAM |

### Disk based computations
Disk based computations are the way to go for large computations. Disk-based computations can be used to control RAM requirements through adjusting `bytes per file`. Also, disk based computations are resumable and fault tolerant: if the program is terminated for any reason, the computation can begin again from the same iteration given that the values in folders `disk_mpz` and `iterations` are preserved.

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

## Contributing
Here is an incomplete list of contribution ideas:
* Break the world record for calculation of regular fractions terms of some important mathematical constant. This is actually easy to do with Zuuv even on a regular PC, given that you're dedicated enough.
* Redesign the UI of Zuuv. Set defaults everywhere so that the user doesn't have to think too much.
* Better parallelization. Zuuv can't parallelize well enough beyond 4-cores. More scalable changes would probably require a algorithmic changes. I'm open to ideas.
* Make RAM usage estimates more precise.
If you're interested in contributing to Zuuv or are considering breaking a world record, do let me know and I will try to help you.
