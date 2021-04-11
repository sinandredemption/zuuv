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
#include "disk_mpz.h"

void benchmark_cf_cruncher(bool verify) {
  double frac_size, f = 0;
  int iters = 0;
  std::cout << "Enter fraction size (in MB): ";
  std::cin >> frac_size;
  std::cout << "Enter # iterations: ";
  std::cin >> iters;
  frac_size *= 1024 * 1024;
  std::cout << "\nBenchmarking against mpz_gcd() using fraction size = " << frac_size / (1024.*1024.)
    << "MB..." << std::endl;

  for (int i = 0; i < iters; ++i) {
    mpz_class a(mpz_rand(frac_size)), b(mpz_rand(frac_size));
    std::cout << "Iteration \t" << i << "/" << iters << std::endl;

    std::cout << "reg_cf_terms()...";
    double start = wall_clock();
    auto c = reg_cf_terms(a, b);
    double end = wall_clock() - start;

    if (c.list.empty())
      c.onload();

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
        std::cerr << " ERROR. Incorrect calculation. ";
        std::cerr << "Do you want to dump the calculations (yes/no)? ";

        std::string choice;
        std::cin >> choice;
        if (choice == "yes") {
          std::cout << "Dumping..." << std::flush;
          CFTerms actual_cf;
          basecase_reg_cf_terms(a, b, actual_cf, false);

          std::ofstream error_file("errors.txt");

          error_file << "a: " << a.get_str() << "\nb: " << b.get_str() << "\ng: "
            << g.get_str() << "\np_k" << c.p_k.get_str() << "\nq_k: " << c.q_k.get_str() << std::endl;
          error_file << "Correct CF: ";
          int i = 0;
          for (auto corr : actual_cf)
            error_file << corr << ((++i % 100 == 0) ? '\n' : ' ');
          i = 0;
          error_file << "\nCF: ";
          for (auto incorr : c.list)
            error_file << incorr << ((++i % 100 == 0) ? '\n' : ' ');
          error_file.close();
          std::cout << " completed" << std::endl;
        }
      } else std::cout << "CORRECT";
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

    std::cout << "type 1...";
    double start = wall_clock();
    ContinuantCache::clear();
    auto c1 = ContinuantCache::cached_continuant(0, cfterms.size() - 1, cfterms.size() / 2, cfterms);
    ContinuantCache::clear();
    double end = wall_clock() - start;
    std::cout << "\t" << end << "ms " << std::endl;

    std::cout << "type 2...";
    start = wall_clock();
    ContinuantCache::clear();
    ContinuantCache::dummy_run(0, cfterms.size() - 1, cfterms.size() / 2);
    auto c2 = ContinuantCache::continuant_new(0, cfterms.size() - 1, cfterms.size() / 2, cfterms);
    ContinuantCache::clear();
    double end2 = wall_clock() - start;
    std::cout << "\t" << end2 << "ms" << std::endl;

    if (end / end2 > 1)
      f += end / end2;

    if (c1 == c2) {
      std::cout << "EQUAL | " << "result size: \t" << mpz_sizeinbase(c1.get_mpz_t(), 2)
        << " bits | multi-core util: \t" << end / end2 << "x" << std::endl;
    }
    else {
      std::cout << "UNEQUAL" << std::endl;
    }
    std::cout << std::endl;
  }

  f /= Iters;
  std::cout << "Multi-core utilization: " << f << "x" << std::endl;
}

void benchmark_mult() {
  size_t bytes_per_file = 0;
  size_t nfiles = 0;
  size_t nthreads = 0;
  size_t iters = 0;

  while (nfiles == 0) {
    std::cout << "# of files: ";
    std::cin >> nfiles;
  }
  while (bytes_per_file == 0) {
    std::cout << "bytes per file: ";
    std::cin >> bytes_per_file;
  }
  while (nthreads == 0) {
    std::cout << "# of threads to use: ";
    std::cin >> nthreads;
  }
  while (iters == 0) {
    std::cout << "# of iterations: ";
    std::cin >> iters;
  }

  std::cout << "Benchmarking disk_mpz multiplication... (frac size = "
    << (2 * nfiles * bytes_per_file) / (1024.*1024.) << "MB, " << nthreads << " threads)...";

  disk_mpq frac("bench", "", bytes_per_file);

  for (size_t i = 0; i < nfiles; ++i) {
    frac.get_den().pushback_mpz(mpz_rand(bytes_per_file));
    frac.get_num().pushback_mpz(mpz_rand(bytes_per_file));
  }

  std::cout << std::endl << std::endl;
  
  for (size_t i = 1; i <= iters; ++i) {
    std::cout << "Iteration " << i << "/" << iters << "...";

    mpz_class a = mpz_rand(bytes_per_file), b = mpz_rand(bytes_per_file);

    double start = wall_clock();
    disk_mpz res(disk_mpz::mt_cross_mult_sub("bench_result", bytes_per_file, nthreads,
      frac.get_num(), a, frac.get_den(), b));
    double end = wall_clock() - start;

    std::cout << "\n" << end << "ms | " <<
      (nfiles * 2. * bytes_per_file * 1000.) / (end * 1024. * 1024.) << "MB/s" << std::endl;

    while (res.files() > 0)
      res.pop_top();

    std::cout << std::endl;
  }
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