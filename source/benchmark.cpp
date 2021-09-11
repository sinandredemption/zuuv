#include "benchmark.h"
#include <mpirxx.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <thread>
#include <random>
#include <ymp/ymp.h>
#include "continuant.h"
#include "cruncher.h"
#include "utils.h"
#include "disk_mpz.h"
#include "mul.h"

void benchmark_cf_cruncher(bool verify)
{
  multiplication::init();

  double frac_size = 0;
  int iterations = 0;

  while (frac_size == 0)
  {
    std::cout << "Enter fraction size in MB (e.g. 0.5): ";
    std::cin >> frac_size;
  }

  while (iterations == 0)
  {
    std::cout << "Enter # iterations: ";
    std::cin >> iterations;
  }

  frac_size *= 1024 * 1024;
  //size_t mem = multiplication::allocate_mem(frac_size * 2);
  //std::cout << "\nusing " << mem / (1024. * 1024.) << "MB of RAM for multiplication" << std::endl;

  std::cout << "benchmarking against mpz_gcd() using fraction size = " << frac_size / (1024.*1024.) << "MB...\n" << std::endl;

  double total_reg_cf_time = 0, total_mpz_gcd_time = 0;

  for (int i = 0; i < iterations; ++i) {
    mpz_class a(mpz_rand(frac_size)), b(mpz_rand(frac_size));
    double start;


    std::cout << "Iteration \t" << i+1 << "/" << iterations << std::endl;

    std::cout << "reg_cf_expansion()...";
    start = wall_clock();

    auto c1 = reg_cf_expansion(a, b);

    double end1 = wall_clock() - start;
    total_reg_cf_time += end1;

    if (c1.list.empty())
      c1.onload();

    std::cout << "\t" << end1 << "ms for " << c1.list.size() << " convergents (largest: "
      << *std::max_element(c1.list.begin(), c1.list.end()) << ")" << std::endl;

    std::cout << "mpz_gcd()...";

    start = wall_clock();

    mpz_class g(0);
    mpz_gcd(g.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());

    double end2 = wall_clock() - start;
    total_mpz_gcd_time += end2;

    std::cout << "    \t" << end2 << "ms for gcd = " << g.get_d() << std::endl;

    if (verify) {
      std::cout << "verifying...";

      if (c1.updated == false) c1.update_continuants();

      if (a != g * c1.continuants.p_k || b != g * c1.continuants.q_k)
      {
        std::cerr << " ERROR. Incorrect calculation. ";
        std::cerr << "Do you want to dump the calculations (yes/no)? ";

        std::string choice;
        std::cin >> choice;

        if (choice == "yes") {
          std::cout << "Dumping..." << std::flush;
          CFTerms actual_cf = basecase_calc_reg_cf_expansion(a, b, false);

          std::ofstream error_file("errors.txt");

          error_file << "a: " << a.get_str() << "\nb: " << b.get_str() << "\ng: "
            << g.get_str() << "\np_k" << c1.continuants.p_k.get_str() << "\nq_k: " << c1.continuants.q_k.get_str() << std::endl;

          error_file << "Correct CF: ";
          
          int i = 0;
          for (auto corr : actual_cf)
            error_file << corr << ((++i % 100 == 0) ? '\n' : ' ');
          
          i = 0;
          error_file << "\nCF: ";
          for (auto incorr : c1.list)
            error_file << incorr << ((++i % 100 == 0) ? '\n' : ' ');
          
          error_file.close();
          std::cout << " completed" << std::endl;
        }
      } else std::cout << " ok";
    }
    std::cout << std::endl;
  }

  std::cout << "Total time taken: " << total_reg_cf_time << "ms [reg_cf_expansion()] | " << total_mpz_gcd_time << "ms [mpz_gcd()]" << std::endl;
  std::cout << "Factor: " << total_reg_cf_time / total_mpz_gcd_time << "x" << std::endl;
  std::cout << "Scratch RAM used for multiplication: " << multiplication::scratch_mem_size / (1024. * 1024.) << "MB" << std::endl;
}

void benchmark_continuant() {
  multiplication::init();

  int iterations = 0, total_terms = 0;

  while (iterations == 0) {
    std::cout << "Enter # iterations: ";
    std::cin >> iterations;
  }

  while (total_terms == 0) {
    std::cout << "Enter total number of random terms (e.g. 100000): ";

    std::string s;
    std::cin >> s;

    total_terms = parse_shorthand_num(s);
  }

  double singlethread_time, multithread_time;
  singlethread_time = multithread_time = 0;

  std::cout << "\nIterations: " << iterations << " | Cores: "
    << std::thread::hardware_concurrency() / 2 << std::endl << std::endl;

  for (int i = 0; i < iterations; ++i) {
    double start;
    CFTerms cfterms;
    cfterms.reserve(total_terms);

    std::cout << "Iteration \t" << i + 1 << "/" << iterations << " \t(" << total_terms << " terms)...";

    // Use inverse transform sampling to generate pseudo-rands over the Gauss-Kurmin distribution
    std::default_random_engine rng;
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    for (int j = 0; j < total_terms; ++j) {
      double d = unif(rng);
      d = exp2(d);
      d = (-2. + d) / (1 - d);   // Inverse transform sampling

      mpir_ui num = (mpir_ui)std::ceil(d);
      cfterms.push_back(num);
    }

    std::cout << std::endl;

    std::cout << "single-threaded...";
    start = wall_clock();

    ContinuantCache::clear();
    ContinuantCache::build_cache_indices(0, cfterms.size() - 1, cfterms.size() / 2);
    auto c1 = ContinuantCache::single_threaded_continuant(0, cfterms.size() - 1, cfterms.size() / 2, cfterms);
    ContinuantCache::clear();

    double end1 = wall_clock() - start;
    std::cout << "\t" << end1 << "ms (" << float(mpz_sizeinbase(c1.get_mpz_t(), 2)) / (end1 * 0.125 * 1.024 * 1024) << "Mb/s)" << std::endl;
    singlethread_time += end1;

    std::cout << "multi-threaded...";
    start = wall_clock();
    
    ContinuantCache::clear();
    ContinuantCache::build_cache_indices(0, cfterms.size() - 1, cfterms.size() / 2);
    auto c2 = ContinuantCache::multi_threaded_continuant(0, cfterms.size() - 1, cfterms.size() / 2, cfterms);
    ContinuantCache::clear();
    
    double end2 = wall_clock() - start;
    std::cout << "\t" << end2 << "ms (" << float(mpz_sizeinbase(c2.get_mpz_t(), 2)) / (end2 * 0.125 * 1.024 * 1024) << "Mb/s)" << std::endl;
    multithread_time += end2;

    if (c1 == c2)
    {
      std::cout << "EQUAL | " << "result size: \t" << mpz_sizeinbase(c1.get_mpz_t(), 2)
        << " bits | multi-core util: " << end1 / end2 << "x" << std::endl;
    }
    else std::cout << "UNEQUAL" << std::endl;
    std::cout << std::endl;
  }

  std::cout << "Total time: " << singlethread_time << "ms (single-thread) | " << multithread_time << "ms (multi-thread)" << std::endl;
  std::cout << "multi-core utilization factor: " << (singlethread_time / multithread_time) << "x" << std::endl;

}

void benchmark_disk_mul() {
  multiplication::init();

  size_t bytes_per_file = 0;
  size_t nfiles = 0;
  size_t threads = 0;
  size_t iters = 0;
  size_t ymp_or_mpir = 0;

  while (nfiles == 0) {
    std::cout << "# of files: ";
    std::cin >> nfiles;
  }
  while (bytes_per_file < 1024 * 1024) {
    std::cout << "file size (MBs): ";
    std::cin >> bytes_per_file;
    bytes_per_file *= 1024 * 1024;
  }
  while (threads == 0) {
    std::cout << "# of threads to use: ";
    std::cin >> threads;
  }
  while (iters == 0) {
    std::cout << "# of iterations: ";
    std::cin >> iters;
  }

  std::cout << "Benchmarking disk_mpz multiplication... ("<< threads << " threads)...";
  std::cout << std::endl << std::endl;

  double total_time = 0;

  for (size_t i = 1; i <= iters; ++i) {
    std::cout << "Iteration " << i << "/" << iters << "...";

    disk_mpq frac("bench", "", bytes_per_file);

    for (size_t i = 0; i < nfiles; ++i) {
      frac.get_den().pushback_mpz(mpz_rand(bytes_per_file));
      frac.get_num().pushback_mpz(mpz_rand(bytes_per_file));
    }

    mpz_class a = mpz_rand(bytes_per_file), b = mpz_rand(bytes_per_file);

    double start = wall_clock();
    disk_mpz res(disk_mpz::cross_mul_sub("bench_result", frac.get_num(), a, frac.get_den(), b, bytes_per_file, threads));
    double end = wall_clock() - start;

    total_time += end;
    std::cout << "\n" << end << "ms | " <<
      (nfiles * 2. * bytes_per_file * 1000.) / (end * 1024. * 1024.) << "MB/s" << std::endl;

    while (res.files() > 0)
      res.pop_top();

    std::cout << std::endl;

    frac.destroy();
  }

  std::cout << "average = " << (iters * nfiles * 2. * bytes_per_file * 1000.) / (total_time * 1024. * 1024.) << "MB/s" << std::endl;
}

int benchmark_ram_mul() {
  multiplication::init();

  size_t mulsize = 0;
  size_t iters   = 0;

  while (mulsize == 0) {
    std::cout << "size of operands: ";

    std::string s;
    std::cin >> s;

    mulsize = parse_shorthand_num(s);
  }
  
  std::cout << "Memory required: " << ymp::LowLevel::mul_iPsize<mp_limb_t>(mulsize, std::thread::hardware_concurrency()) / (1024. * 1024.);
  std::cout << "MB" << std::endl;

  while (iters == 0) {
    std::cout << "# of repititions: ";
    std::cin >> iters;
  }

  double total_mpir_time, total_ymp1_time, total_ymp2_time;
  total_mpir_time = total_ymp1_time = total_ymp2_time = 0;

  bool is_ok = true;
  mpz_class a, b, c1, c2, c3;
  for (int i = 1; i <= iters; ++i)
  {
    if (i % 2) mulsize *= 1. + float(rand()) / (10 * float(RAND_MAX));
    else       mulsize /= 1. + float(rand()) / (10 * float(RAND_MAX));

    a = mpz_rand(mulsize);
    b = mpz_rand(mulsize);

    std::cout << "Iteration " << i << "/" << iters << " (size = " << mulsize * 8. / (1024. * 1024.) << "MB)...\n";

    // MPIR
    auto mpir_time = wall_clock();
    mpz_mul(c1.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
    mpir_time = wall_clock() - mpir_time;

    std::cout << "mpir                  : " << mpir_time << "ms" << std::endl;

    total_mpir_time += mpir_time;

    // YMP - single-threaded
    auto ymp1_time = wall_clock();
    multiplication::mul(c2.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t(), 1);
    ymp1_time = wall_clock() - ymp1_time;

    std::cout << "ymp (single-threaded) : " << ymp1_time << "ms";
    std::cout << " (mem = " << (ymp::LowLevel::mul_iPsize<mp_limb_t>(mulsize * 2, 1)) / (8 * 1024 * 1024) << "Mb)" << std::endl;

    total_ymp1_time += ymp1_time;

    // YMP - multi-threaded
    auto ymp2_time = wall_clock();
    multiplication::mul(c3.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t(), std::thread::hardware_concurrency());
    ymp2_time = wall_clock() - ymp2_time;

    std::cout << "ymp (multi-threaded)  : " << ymp2_time << "ms";
    std::cout << " (mem = " << (ymp::LowLevel::mul_iPsize<mp_limb_t>(mulsize * 2, std::thread::hardware_concurrency())) / (8*1024 * 1024) << "Mb)" << std::endl;

    total_ymp2_time += ymp2_time;

    is_ok = is_ok && (c1 == c2) && (c2 == c3);

    std::cout << (is_ok ? "ok" : "ERROR") << std::endl;
    
    std::cout << std::endl;

    if (!is_ok) break;
  }

  std::cout << "mpir                 : " << total_mpir_time << "ms" << std::endl;
  std::cout << "ymp (single-threaded): " << total_ymp1_time << "ms" << std::endl;
  std::cout << "ymp (multi-threaded) : " << total_ymp2_time << "ms" << std::endl;

  std::cout << "done" << std::endl;
  return is_ok ? 0 : -1;
}

void benchmark_estimate_time() {
  double fraction_size = 0;

  while (fraction_size == 0) {
    std::cout << "Enter fraction size in MBs (e.g. 50.0): ";
    std::cin >> fraction_size;
  }

  auto prev_flags = std::cout.flags();
  std::cout.precision(2);
  std::cout << std::fixed;

  std::cout << "Estimated Peak RAM usage: " << peak_ram_usage(fraction_size) << "MB" << std::endl;
  std::cout << std::endl;

  char c = 0;
  while ((c != 'y') && (c != 'n')) {
    std::cout << "Do you also want to estimate computation time? (y/n): ";
    std::cin >> c;
  }

  if (c == 'y')
  {
    multiplication::init();

    size_t digits = 0, terms = 0;

    while (digits < 1) {
      std::cout << "Number of hex digits in the input file: ";

      std::string s;
      std::cin >> s;

      digits = parse_shorthand_num(s);
    }

    while (terms < 1) {
      std::cout << "Number of output terms: ";
      
      std::string s;
      std::cin >> s;

      terms = parse_shorthand_num(s);
    }

    std::cout << "Estimating time required for reg_cf_expansion()...";

    const int Repititions = 3;
    const int DefaultFractionSize1 = 1 * 1024 * 1024; // 1MB
    const int DefaultFractionSize2 = 10 * 1024 * 1024; // 1MB
    double average_time = 0, estimated_time = 0, estimated_time1, estimated_time2;

    for (int i = 0; i < Repititions; ++i)
    {
      mpz_class num = mpz_rand(DefaultFractionSize1), den = mpz_rand(DefaultFractionSize1);

      double t = wall_clock();
      reg_cf_expansion(num, den, true, true);
      t = wall_clock() - t;

      average_time += t;
    }

    average_time /= Repititions;

    estimated_time = (average_time / ((DefaultFractionSize1/(1024.*1024.)) * 10 * std::sqrt(std::log(10 * (DefaultFractionSize1 / (1024. * 1024.))))));
    estimated_time *= fraction_size * 10.0 * std::sqrt(std::log(fraction_size));
    estimated_time1 = estimated_time;

    std::cout << " " << estimated_time1 / 1000. << "s per iteration";
    std::cout << std::endl;

    average_time = 0;

    std::cout << "Estimating time required for cross-multiplication...";
    for (int i = 0; i < Repititions; ++i) {
      disk_mpq frac("bench", "", DefaultFractionSize2);

      for (size_t j = 0; j < 5; ++j) {
        frac.get_den().pushback_mpz(mpz_rand(DefaultFractionSize2));
        frac.get_num().pushback_mpz(mpz_rand(DefaultFractionSize2));
      }

      mpz_class a = mpz_rand(DefaultFractionSize2), b = mpz_rand(DefaultFractionSize2);

      double t = wall_clock();
      disk_mpz res(disk_mpz::cross_mul_sub("bench_result", frac.get_num(), a, frac.get_den(), b, DefaultFractionSize2));
      average_time += wall_clock() - t;

      while (res.files() > 0)
        res.pop_top();

      frac.destroy();
    }

    average_time /= Repititions;

    size_t nfiles = (digits / (2. * 1024. * 1024.)) / fraction_size;
    estimated_time = (nfiles / 5) * (fraction_size / (DefaultFractionSize2 / (1024. * 1024.))) * average_time * 2;
    estimated_time2 = estimated_time;

    std::cout << " " << estimated_time2 / 1000 << "s per iteration";
    std::cout << std::endl;

    size_t total_iterations = terms / (fraction_size * 4.9e6);

    std::cout << "Total time = " << format_time(estimated_time1 + estimated_time2) << " x " << total_iterations;
    std::cout << " = " << format_time(total_iterations * (estimated_time1 + estimated_time2));
    std::cout << std::endl;
  }
}