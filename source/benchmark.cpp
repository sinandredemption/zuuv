#include "benchmark.h"
#include <mpirxx.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <thread>
#include <random>
#include "continuant.h"
#include "cruncher.h"
#include "utils.h"

void benchmark(bool verify) {
  double ram_target, f = 0;
  int iters = 0;
  std::cout << "Enter target RAM usage (in MB): ";
  std::cin >> ram_target;
  std::cout << "Enter # iterations: ";
  std::cin >> iters;
  size_t frac_size = ram_target * 12825;  // fixme
  std::cout << "\nBenchmarking against mpz_gcd() using target RAM = " << ram_target
    << "MB (frac size = " << frac_size << ")..." << std::endl;

  for (int i = 0; i < iters; ++i) {
    mpz_class a(mpz_rand(frac_size)), b(mpz_rand(frac_size));
    std::cout << "Iteration \t" << i << "/" << iters << std::endl;

    std::cout << "reg_cf_terms()...";
    double start = wall_clock();
    auto c = reg_cf_terms(a, b);
    double end = wall_clock() - start;
    std::cout << "\t" << end << "ms for " << c.list.size() << " convergents (largest: "
      << *std::max_element(c.list.begin(), c.list.end()) << ")" << std::endl;

    std::cout << "mpz_gcd()...";
    start = wall_clock();
    mpz_class g(0);
    mpz_gcd(g.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    double end2 = wall_clock() - start;
    std::cout << "    \t" << end2 << "ms for gcd(" << g.get_d() << ") (f: " << end / end2 << ")" << std::endl;

    f += end / end2;

    if (verify) {
      std::cout << "Verifying...";
      if (c.updated == false) c.update();
      if (a != g * c.p_k || b != g * c.q_k) {
        std::cerr << " ERROR. Incorrect calculation.";
        std::cerr << "Press 'y' to dump the calculations: ";

        char yes;
        std::cin >> yes;
        if (yes == 'y') {
          CFTerms actual_cf;
          basecase_reg_cf_terms(a, b, actual_cf, false);

          std::ofstream error_file("errors.txt");

          error_file << "a: " << a.get_str() << "\nb: " << b.get_str() << "\ng: "
            << g.get_str() << "\np_k" << c.p_k.get_str() << "\nq_k: " << c.q_k.get_str() << std::endl;
          error_file << "Correct CF: ";
          for (auto corr : actual_cf)
            error_file << corr << ' ';
          error_file << "\nCF: ";
          for (auto incorr : c.list)
            error_file << incorr << ' ';
          error_file.close();
        }
      }

      std::cout << "CORRECT";
    }
    std::cout << std::endl;
  }
  f /= iters;
  std::cout << "Average f: " << f << std::endl;
}

void benchmark_continuant() {
  const int Iters = 10;
  const int StartTerms = 50000;
  double f = 0;

  std::cout << "\nIterations: " << Iters << " | Cores: "
    << std::thread::hardware_concurrency() / 2 << std::endl << std::endl;

  for (int i = 0; i < Iters; ++i) {
    int terms = (StartTerms << i);
    std::default_random_engine rng;
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    CFTerms cfterms;

    std::wcout << "Iteration \t" << i + 1 << "/" << Iters << " \t(" << terms << " terms)...";
    cfterms.reserve(terms);
    for (int j = 0; j < terms; ++j) {
      // Use inverse transform sampling to generate pseudo-rands over the Gauss-Kurmin distribution

      double d = unif(rng);
      d = exp2(d);
      d = (-2. + d) / (1 - d);   // Inverse transform sampling

      mpir_ui num = (mpir_ui)std::ceil(d);
      cfterms.push_back(num);
    }
    std::cout << std::endl;

    std::cout << "single-core...";
    double start = wall_clock();
    ContinuantCache::clear();
    auto c1 = ContinuantCache::cached_continuant(0, cfterms.size() - 1, cfterms.size() / 2, cfterms);
    ContinuantCache::clear();
    double end = wall_clock() - start;
    std::cout << "\t" << end << "ms " << std::endl;

    std::cout << "multi-core...";
    start = wall_clock();
    ContinuantCache::clear();
    auto c2 = parallel_continuant(0, cfterms.size() - 1, cfterms);
    ContinuantCache::clear();
    double end2 = wall_clock() - start;
    std::cout << "\t" << end2 << "ms" << std::endl;

    if (end / end2 > 1)
      f += end / end2;

    if (c1 == c2) {
      std::cout << "EQUAL | " << "result size: \t" << mpz_sizeinbase(c1.get_mpz_t(), 2) << " bits | multi-core util: \t"
        << end / end2 << "x" << std::endl;
    }
    else {
      std::cout << "UNEQUAL" << std::endl;
    }
    std::cout << std::endl;
  }

  f /= Iters;
  std::cout << "Multi-core utilization: " << f << "x" << std::endl;
}

void benchmark_basic_continuant() {
  CFTerms cl;
  for (int j = 0; j < Params::ThresholdUseBasicContProc * 2; ++j)
    cl.push_back(10. * rand() / double(RAND_MAX));
  for (int i = 0; i < 13; ++i) {
    int n = (10000 << i);
    std::cout << "Round #" << i << " n: " << n << "...";
    std::cout << std::endl;

    double start = wall_clock();
    size_t checksum = 0;
    for (int j = 0; j < n; ++j) {
      //auto c1 = basic_continuant(0, cfterms.size() - 1, cfterms);
      ContinuantCache::clear();
      auto c1 = ContinuantCache::cached_continuant(0, cl.size() - 1, cl.size() / 2, cl);
      ContinuantCache::clear();
      checksum += mpz_sizeinbase(c1.get_mpz_t(), 2);
    }
    double end = wall_clock() - start;
    std::cout << "\n" << end << "ms" << std::endl;

    //start = wall_clock();
    //checksum = 0;
    //for (int j = 0; j < n; ++j) {
    //	//auto c1 = parallel_basic_continuant(0, cfterms.size() - 1, cfterms);
    //	ContinuantCache::clear();
    //	auto c1 = ContinuantCache::parallel_continuant(0, 0, cfterms.size() - 1, cfterms.size() / 2, cfterms);
    //	ContinuantCache::clear();
    //	checksum += mpz_sizeinbase(c1.get_mpz_t(), 2);
    //}
    //double end2 = wall_clock() - start;
    //std::cout << end2 << "ms (f: " << end / end2 << ")" << std::endl;

    std::cout << std::endl;
  }
}