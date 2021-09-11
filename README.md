# Zuuv
Multiprecision Floating-point to Continued Fraction Cruncher

## Introduction
Zuuv is a program that can compute the [regular continued fraction](https://en.wikipedia.org/wiki/Continued_fraction) of a number from it's floating-point representation (a file containing hex or dec digits). Its features include:
* Advanced algorithms for log-linear runtime
* Ability to utilize hard-disk to perform extremely large calculations
* Ability to automatically save and resume partially computations
* Multi-threading

## Notable Computations
Zuuv has been utilized to break world records for computation of regular continued fraction terms of important mathematical constants such as Pi, Euler-Mascheroni etc, even on regular PCs. Here are the current record computations done using Zuuv:

| Constant | Terms | Date | Time | System Details | Link | Comments |
| -------- | ----- | ---- | ---- | -------------- | ---- | ---- |
| Pi       | 30,113,021,586 | Apr 2021 | ~6 days | Intel Core i5 9300H, 8GB DDR4 @ 2666MHz | [All of the terms](https://drive.google.com/drive/folders/1--Qh9Xxq1i6oeHnTXzKrQ9FoguHreBKy?usp=sharing) | The computation added a new term to the OEIS sequence [A033089](https://oeis.org/A033089) |
| Euler's constant | 16,695,279,010 | Apr 2021 | ~2 days | Intel Core i5 9300H, 8GB DDR4 @ 2666MHz | [All of the terms](https://drive.google.com/drive/folders/15LjUnGcZiDAJLptqwl_Wf69zFLV0HkIK?usp=sharing) | See comments on [A033091](https://oeis.org/A033091) |

## Installation
Windows builds, along with the required DLLs can be found in the `builds` directory.
- `modern` subdirectory contains builds for newer processors
- `legacy` subdirectory contains builds for older processors

### Linux
_IMPORTANT:_ Starting Release 2.0, **LINUX IS NO LONGER SUPPORTED**
For older versions of Zuuv, Linux build might require [mpir-3.0.0](https://mpir.org/mpir-3.0.0.zip) library to be installed.
1. Download and install **mpir** library: `wget https://mpir.org/mpir-3.0.0.zip`
2. Navigate to `mpir-3.0.0` folder and run `./configure --enable-cxx`. Install all the requirements with `sudo apt install <program>` in the case that `./configure` fails due to some missing program.
3. Run `make install`

### Note on performance
The program requires a minimum of **four threads** to run efficiently. I haven't personally tested it on a dual-core setup, but I would expect performance issues.
After downloading, to check if everything is alright, do a quick continuant benchmark (option 4) to see if everything is running fine.

## Doing an actual computation
Zuuv takes in a file containing the stream of digits of a number and calculates its regular continued fraction terms in multiple iterations, each stored in the subdirectory `iterations`. A file `pi_1m_hex.txt`, containing one million hexadecimal digits of pi, is provided in the `source/constants` subdirectory to experiment with Zuuv before launching a serious computation.
It is recommended to estimate the resources and the time a computation would need before starting a computation. The option "Estimate RAM usage and computation times" in Zuuv would provide you with a very rough but practical estimate.

#### How many terms can I compute from X digits of pi?
As a rule of thumb, you can compute `1.167 x <number of hex digits>` terms or `0.97 x <number of dec digits>` terms of simple continued fraction. For example, a file containing five billion hex digits of pi can be used to compute `5 billion x 1.167 = 5.835 billion` terms of simple continued fraction of pi.

## Computation modes
There are two computation modes.

| Computation Mode | RAM requirements | Disk requirements | Resumable |  Speed  | Comments |
| ---------------- | ---------------- | ----------------- | --------- | ------- | -------- |
|     RAM based    |       High       |    Output only    |     No    | Fastest | Use for small computations that can fit in RAM |
|    Disk based    | Adjustable (low) |       High        |    Yes    |   Fast  | Use for large computations that can't entirely be done in RAM |

### Disk based computations
Disk based computations are the way to go for large computations. Disk-based computations can be used to control RAM requirements through adjusting `bytes per file`. Also, disk based computations are resumable and fault tolerant: if the program is terminated for any reason, the computation can begin again from the same iteration given that the values in folders `disk_mpz` and `iterations` are preserved.

## Benchmarking
There are four default benchmarks provided:
### Benchmark Continued Fraction Cruncher
Zuuv attempts to calculate the terms of the regular continued fraction of a random fraction. This benchmark isn't _very_ scalable beyond 4-cores (the procedure is essentially sequential in nature, and hence difficult to parallelize beyond 4-cores. However, I am open to ideas.).
### Benchmark Continuant Cruncher
Zuuv attempts to calculate the continued of a list of random numbers (generated according to the Gauss-Kuzmin distribution to simulate continued fraction terms). The benchmark timing should benefit with increase in processor cores, as the procedure is well-parallelized.
### Benchmark Disk-based Multiplication
Zuuv attempts to benchmark a critical multiplication routine required in disk-based computations. This procedure again is well parallelized and should benefit well from scaling.
### Benchmark RAM-based Multiplication
Zuuv compares the speed of MPIR-based multiplication with YMP-based multiplication.

## Compiling the newer versions (Release 2.0 or greater)
This repository contains all the dependencies needed for a compilation under Microsoft Visual Studio.
1. Clone this repository
2. Open the project file "zuuv.vcxproj"
3. Run the compilation

## Compiling the older versions (Release 1.0)
Zuuv is written in **C++17** and depends on the **mpir-3.0.0 library.**
### Visual Studio
You will need to download and extract the mpir-3.0.0 library, and enable C++17 standard.
1. Clone this repository.
2. Goto _File -> New -> Project from Existing Code_. Choose the `source` folder
#### mpir
3. Download and extract the [mpir-3.0.0](https://mpir.org/mpir-3.0.0.zip) library.
4. In _Project -> Properties -> C/C++ -> General -> Additional Include Directories_, add the path to `mpir` sub-directory inside `mpir-3.0.0` folder
5. In _Project -> Properties -> Linker -> Input -> Additional Dependencies_, add path to `mpir.lib`
#### C++17
6. Goto _Project -> Properties -> C/C++ -> Language -> C++ Language Standard_, and choose "_ISO C++17 Standard_"
### g++
1. Install the **mpir-3.0.0 library** following the instructions given under _"[Installation](#installation)" -> "Linux"_
2. Navigate to `source` folder and use the following command-line

`g++ *.cpp -I /usr/lib/local/include -L /usr/lib/local/lib -lmpir -pthread -std=c++17 -Ofast -DNDEBUG -o zuuv`

In case of problems, feel free to write to me.

## Contributing
Here is an incomplete list of contribution ideas:
* **Break some world records with Zuuv.** Calculate regular continued fraction expansions of some important mathematical constant. This is actually easy to do with Zuuv even on a regular PC, given that you're dedicated enough. Given the attention computation of decimal digits of constants gets, it is a shame that continued fraction exapansions haven't catched up, especially considering that continued fraction expansions are more mathematically interesting.
* **Redesign the UI of Zuuv.** Set defaults everywhere so that the user doesn't have to think too much.
* **Better parallelization.** Crunching continued fraction expansions is inherently sequential in nature. As such, parallelization might require algorithmic changes. Zuuv still doesn't use often use more than 50% of the available CPU multithreading power on average. However, maybe certain routines can be parallized via a GPU (e.g. basecase calculations of continuants of a list of numbers).

If you're interested in contributing to Zuuv or are considering breaking a world record, do let me know, and I will try to help you.
