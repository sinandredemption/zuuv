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

## Compiling
Zuuv is written in C++17 and depends on the mpir-3.0.0 library. You will need to download and extract it.
### Visual Studio
